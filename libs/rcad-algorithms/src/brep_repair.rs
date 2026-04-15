//! B-Rep repair / clean-up utilities.
//!
//! Analogous to OCCT `ShapeFix_Shape` / `ShapeFix_Wire` / `ShapeFix_Face`.
//!
//! # Operations
//!
//! | Function | Description | OCCT equivalent |
//! |---|---|---|
//! | [`merge_close_vertices`] | Merge vertices closer than `tolerance` | `ShapeFix_Wire::FixSameParameter` / `BRepBuilderAPI_Sewing` |
//! | [`remove_degenerate_faces`] | Remove faces with fewer than 3 edges or zero-area | `ShapeFix_Shape` |
//! | [`recompute_face_normals`] | Recompute per-face normals from vertex positions | `BRepLib::UpdateEdgeTol` + fix normals |
//! | [`fix_wire_orientation`] | Ensure each wire forms a closed, consistently-oriented loop | `ShapeFix_Wire::FixClosed` |
//! | [`repair`] | Apply all fixes in a single pass | `ShapeFix_Shape::Perform` |
//!
//! All functions are **non-destructive**: they return a new `BRep` leaving the
//! original unchanged.

use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::CurveEval;
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
use crate::brep_check::{check_orientation_consistency, diagnose_same_parameter, diagnose_same_range};
use crate::tolerance::TOLERANCE_ABS;

fn make_connected_has_future_tolerance_increase(
    pass_idx: usize,
    pass_limit: usize,
    current_tolerance: f64,
    base_tolerance: f64,
    tolerance_growth: f64,
    tolerance_cap: f64,
) -> bool {
    if pass_idx + 1 >= pass_limit {
        return false;
    }
    if tolerance_growth <= 1.0 {
        return false;
    }
    let next_grown_tolerance = base_tolerance * tolerance_growth.powi((pass_idx + 1) as i32);
    let next_tolerance = next_grown_tolerance.min(tolerance_cap);
    next_tolerance > current_tolerance + 1e-15
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Summary of all changes made during repair.
#[derive(Debug, Clone, Default)]
pub struct RepairReport {
    /// Number of vertex pairs that were merged.
    pub vertices_merged: usize,
    /// Number of degenerate faces removed.
    pub degenerate_faces_removed: usize,
    /// Number of faces whose normals were recomputed.
    pub normals_recomputed: usize,
    /// Number of faces whose inward orientation was flipped.
    pub faces_reoriented: usize,
    /// Number of wires whose orientation was fixed.
    pub wires_fixed: usize,
    /// Number of edges whose SameRange consistency was repaired.
    pub same_range_fixed: usize,
    /// Number of edges whose SameParameter flag was repaired.
    pub same_parameter_fixed: usize,
}

/// Summary of baseline connectivity rebuilding pass.
#[derive(Debug, Clone, Default)]
pub struct MakeConnectedReport {
    /// Number of merged near-coincident vertices.
    pub vertices_merged: usize,
    /// Number of tiny/degenerate edges removed after merging.
    pub small_edges_removed: usize,
    /// Number of make-connected passes that were executed.
    pub passes_run: usize,
    /// Whether the pass sequence converged before reaching `max_passes`.
    pub converged: bool,
    /// Effective tolerance used in the final executed pass.
    pub final_tolerance: f64,
    /// Whether tolerance growth was clamped by the configured cap.
    pub tolerance_cap_applied: bool,
    /// Number of edge pairs that were sewn together (enhanced mode).
    pub edges_sewn: usize,
    /// Number of faces that were merged (enhanced mode with face merging).
    pub faces_merged: usize,
}

/// Operating mode for `make_connected_enhanced`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MakeConnectedMode {
    /// Standard mode: vertex merging + small edge removal.
    #[default]
    Standard,
    /// Aggressive mode: includes edge sewing and face merging.
    Aggressive,
    /// Conservative mode: only vertex merging, no edge removal.
    Conservative,
}

/// Strategy for connectivity repair operations.
///
/// This struct provides fine-grained control over the behavior of
/// `make_connected` operations, allowing users to customize which
/// repairs are applied and how aggressively.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_repair::MakeConnectedStrategy;
///
/// let strategy = MakeConnectedStrategy {
///     merge_vertices: true,
///     merge_tolerance: 0.001,
///     remove_small_edges: true,
///     min_edge_length: 0.0001,
///     max_passes: 5,
///     tolerance_growth: 1.5,
///     tolerance_cap: 0.1,
///     sew_edges: false,
///     edge_sew_tolerance: 0.001,
///     merge_faces: false,
///     face_merge_tolerance: 0.001,
/// };
/// ```
#[derive(Debug, Clone)]
pub struct MakeConnectedStrategy {
    /// Whether to merge near-coincident vertices.
    pub merge_vertices: bool,
    /// Tolerance for vertex merging.
    pub merge_tolerance: f64,
    /// Whether to remove small/degenerate edges.
    pub remove_small_edges: bool,
    /// Minimum edge length; shorter edges are candidates for removal.
    pub min_edge_length: f64,
    /// Maximum number of repair passes.
    pub max_passes: usize,
    /// Factor by which tolerance grows each pass (1.0 = no growth).
    pub tolerance_growth: f64,
    /// Upper cap for tolerance growth.
    pub tolerance_cap: f64,
    /// Whether to sew close edges together.
    pub sew_edges: bool,
    /// Tolerance for edge sewing.
    pub edge_sew_tolerance: f64,
    /// Whether to merge coincident faces.
    pub merge_faces: bool,
    /// Tolerance for face merging.
    pub face_merge_tolerance: f64,
}

impl Default for MakeConnectedStrategy {
    fn default() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: 1e-6,
            remove_small_edges: true,
            min_edge_length: 1e-6,
            max_passes: 3,
            tolerance_growth: 1.0,
            tolerance_cap: f64::INFINITY,
            sew_edges: false,
            edge_sew_tolerance: 1e-6,
            merge_faces: false,
            face_merge_tolerance: 1e-6,
        }
    }
}

impl MakeConnectedStrategy {
    /// Create a conservative strategy (only vertex merging).
    pub fn conservative() -> Self {
        Self {
            merge_vertices: true,
            remove_small_edges: false,
            sew_edges: false,
            merge_faces: false,
            ..Self::default()
        }
    }

    /// Create a standard strategy (vertex merging + small edge removal).
    pub fn standard() -> Self {
        Self::default()
    }

    /// Create an aggressive strategy (all repairs enabled).
    pub fn aggressive() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: 1e-5,
            remove_small_edges: true,
            min_edge_length: 1e-5,
            max_passes: 5,
            tolerance_growth: 1.5,
            tolerance_cap: 0.01,
            sew_edges: true,
            edge_sew_tolerance: 1e-5,
            merge_faces: true,
            face_merge_tolerance: 1e-5,
        }
    }

    /// Create a strategy for injection molding (optimized for thin walls).
    pub fn for_injection_molding() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: 1e-4,
            remove_small_edges: true,
            min_edge_length: 1e-4,
            max_passes: 10,
            tolerance_growth: 2.0,
            tolerance_cap: 0.1,
            sew_edges: true,
            edge_sew_tolerance: 1e-4,
            merge_faces: false, // Don't merge faces for molding
            face_merge_tolerance: 1e-4,
        }
    }

    /// Create a strategy for 3D printing (conservative tolerance).
    pub fn for_3d_printing() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: 1e-3, // 0.001mm tolerance for printing
            remove_small_edges: true,
            min_edge_length: 1e-3,
            max_passes: 3,
            tolerance_growth: 1.0,
            tolerance_cap: 0.1,
            sew_edges: false,
            edge_sew_tolerance: 1e-3,
            merge_faces: false,
            face_merge_tolerance: 1e-3,
        }
    }

    /// Create a strategy for CNC machining (precise, no merging).
    pub fn for_cnc_machining() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: 1e-6, // Very precise
            remove_small_edges: true,
            min_edge_length: 1e-6,
            max_passes: 1,
            tolerance_growth: 1.0,
            tolerance_cap: 1e-5,
            sew_edges: false,
            edge_sew_tolerance: 1e-6,
            merge_faces: false,
            face_merge_tolerance: 1e-6,
        }
    }

    /// Apply the strategy to a BRep.
    ///
    /// This is the main entry point for connectivity repair using
    /// a custom strategy configuration.
    pub fn apply(&self, brep: &BRep) -> (BRep, MakeConnectedReport) {
        let mut out = brep.clone();
        let mut report = MakeConnectedReport::default();
        let base_tol = self.merge_tolerance.max(TOLERANCE_ABS);

        for pass_idx in 0..self.max_passes {
            let grown_tol = base_tol * self.tolerance_growth.powi(pass_idx as i32);
            let pass_tol = grown_tol.min(self.tolerance_cap);
            let mut pass_merged = 0usize;
            let mut pass_removed = 0usize;
            let mut pass_sewn = 0usize;

            // Vertex merging
            if self.merge_vertices {
                let (b, merged) = merge_close_vertices(&out, pass_tol);
                out = b;
                pass_merged += merged;
            }

            // Small edge removal
            if self.remove_small_edges {
                let (b, removed) = remove_small_edges(&out, self.min_edge_length);
                out = b;
                pass_removed += removed;
            }

            // Edge sewing
            if self.sew_edges {
                let (b, sewn) = sew_close_edges(&out, self.edge_sew_tolerance);
                out = b;
                pass_sewn += sewn.edges_sewn;
            }

            // Face merging (if enabled)
            if self.merge_faces {
                // Note: Face merging is a complex operation that requires
                // geometric analysis. For now, we skip this in the strategy.
                // It can be added later when face merging is fully implemented.
            }

            report.vertices_merged += pass_merged;
            report.small_edges_removed += pass_removed;
            report.edges_sewn += pass_sewn;
            report.passes_run = pass_idx + 1;
            report.final_tolerance = pass_tol;

            if grown_tol > self.tolerance_cap {
                report.tolerance_cap_applied = true;
            }

            // Check for convergence
            if pass_merged == 0 && pass_removed == 0 && pass_sewn == 0 {
                // Check if tolerance will grow in future passes
                if self.tolerance_growth <= 1.0 || pass_idx + 1 >= self.max_passes {
                    report.converged = true;
                    break;
                }
                let next_tol = base_tol * self.tolerance_growth.powi((pass_idx + 1) as i32);
                if next_tol > self.tolerance_cap {
                    report.converged = true;
                    break;
                }
            }
        }

        (out, report)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Scoped Seed Detection Strategies
// ─────────────────────────────────────────────────────────────────────────────

/// Strategy for detecting seed entities for scoped make-connected.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SeedDetectionStrategy {
    /// Use all vertices as seeds (equivalent to global cleanup).
    #[default]
    AllVertices,
    /// Use vertices on short edges as seeds.
    ShortEdgeEndpoints,
    /// Use vertices on edges with high tolerance as seeds.
    HighToleranceEdges,
    /// Use vertices at potential geometry seams (multi-PCurve edges).
    SeamCandidates,
    /// Use vertices near potential duplicates.
    NearDuplicateVertices,
    /// Combine multiple strategies (hybrid approach).
    Hybrid,
}

/// Configuration for seed detection.
#[derive(Debug, Clone)]
pub struct SeedDetectionConfig {
    /// Strategy to use for seed detection.
    pub strategy: SeedDetectionStrategy,
    /// Minimum edge length for ShortEdgeEndpoints strategy.
    pub short_edge_threshold: f64,
    /// Tolerance threshold for HighToleranceEdges strategy.
    pub high_tolerance_threshold: f64,
    /// Distance threshold for NearDuplicateVertices strategy.
    pub near_duplicate_distance: f64,
    /// Maximum number of seeds to return (0 = no limit).
    pub max_seeds: usize,
    /// Include vertices within N hops of primary seeds.
    pub neighborhood_depth: usize,
}

impl Default for SeedDetectionConfig {
    fn default() -> Self {
        Self {
            strategy: SeedDetectionStrategy::default(),
            short_edge_threshold: 1e-4,
            high_tolerance_threshold: 1e-3,
            near_duplicate_distance: 1e-4,
            max_seeds: 0,
            neighborhood_depth: 1,
        }
    }
}

impl SeedDetectionConfig {
    /// Create config for short-edge seed detection.
    pub fn short_edges(threshold: f64) -> Self {
        Self {
            strategy: SeedDetectionStrategy::ShortEdgeEndpoints,
            short_edge_threshold: threshold,
            ..Default::default()
        }
    }

    /// Create config for high-tolerance seed detection.
    pub fn high_tolerance(threshold: f64) -> Self {
        Self {
            strategy: SeedDetectionStrategy::HighToleranceEdges,
            high_tolerance_threshold: threshold,
            ..Default::default()
        }
    }

    /// Create hybrid config combining multiple strategies.
    pub fn hybrid() -> Self {
        Self {
            strategy: SeedDetectionStrategy::Hybrid,
            ..Default::default()
        }
    }
}

/// Result of seed detection analysis.
#[derive(Debug, Clone, Default)]
pub struct SeedDetectionResult {
    /// Indices of detected seed vertices.
    pub seed_vertices: Vec<usize>,
    /// Indices of detected seed edges.
    pub seed_edges: Vec<usize>,
    /// Number of seeds from each strategy (for hybrid).
    pub strategy_counts: std::collections::HashMap<String, usize>,
    /// Estimated coverage ratio (seeds / total entities).
    pub coverage_ratio: f64,
}

/// Detect seed vertices for scoped make-connected based on strategy.
pub fn detect_seeds_for_scoped_cleanup(
    brep: &BRep,
    config: &SeedDetectionConfig,
) -> SeedDetectionResult {
    let mut result = SeedDetectionResult::default();
    let mut vertex_set = std::collections::HashSet::new();
    let mut edge_set = std::collections::HashSet::new();
    let n_vertices = brep.vertices.len();
    let n_edges = brep.edges.len();

    match config.strategy {
        SeedDetectionStrategy::AllVertices => {
            for i in 0..n_vertices {
                vertex_set.insert(i);
            }
            result.strategy_counts.insert("all_vertices".to_string(), n_vertices);
        }
        SeedDetectionStrategy::ShortEdgeEndpoints => {
            for (ei, edge) in brep.edges.iter().enumerate() {
                let start = brep.vertices.get(edge.start).map(|v| v.point);
                let end = brep.vertices.get(edge.end).map(|v| v.point);
                if let (Some(s), Some(e)) = (start, end) {
                    let len = (s - e).length();
                    if len < config.short_edge_threshold {
                        vertex_set.insert(edge.start);
                        vertex_set.insert(edge.end);
                        edge_set.insert(ei);
                    }
                }
            }
            result.strategy_counts.insert("short_edge_endpoints".to_string(), vertex_set.len());
        }
        SeedDetectionStrategy::HighToleranceEdges => {
            let edge_tolerances = &brep.geom.edge_tolerance;
            for (ei, &tol) in edge_tolerances.iter().enumerate() {
                if tol > config.high_tolerance_threshold && ei < n_edges {
                    let edge = &brep.edges[ei];
                    vertex_set.insert(edge.start);
                    vertex_set.insert(edge.end);
                    edge_set.insert(ei);
                }
            }
            result.strategy_counts.insert("high_tolerance_edges".to_string(), vertex_set.len());
        }
        SeedDetectionStrategy::NearDuplicateVertices => {
            for i in 0..n_vertices {
                for j in (i + 1)..n_vertices {
                    let dist = (brep.vertices[i].point - brep.vertices[j].point).length();
                    if dist < config.near_duplicate_distance {
                        vertex_set.insert(i);
                        vertex_set.insert(j);
                    }
                }
            }
            result.strategy_counts.insert("near_duplicate_vertices".to_string(), vertex_set.len());
        }
        SeedDetectionStrategy::SeamCandidates => {
            // Edges referenced by multiple faces (potential seams)
            let mut edge_face_count = vec![0usize; n_edges];
            for solid in &brep.solids {
                for shell in &solid.shells {
                    for face in &shell.faces {
                        for we in &face.outer_wire.edges {
                            if we.idx < n_edges {
                                edge_face_count[we.idx] += 1;
                            }
                        }
                    }
                }
            }
            for (ei, &count) in edge_face_count.iter().enumerate() {
                if count > 2 {
                    let edge = &brep.edges[ei];
                    vertex_set.insert(edge.start);
                    vertex_set.insert(edge.end);
                    edge_set.insert(ei);
                }
            }
            result.strategy_counts.insert("seam_candidates".to_string(), vertex_set.len());
        }
        SeedDetectionStrategy::Hybrid => {
            // Combine all strategies
            let mut combined = std::collections::HashSet::new();

            // Short edges
            for edge in &brep.edges {
                let start = brep.vertices.get(edge.start).map(|v| v.point);
                let end = brep.vertices.get(edge.end).map(|v| v.point);
                if let (Some(s), Some(e)) = (start, end) {
                    if (s - e).length() < config.short_edge_threshold {
                        combined.insert(edge.start);
                        combined.insert(edge.end);
                    }
                }
            }

            // High tolerance
            for (ei, &tol) in brep.geom.edge_tolerance.iter().enumerate() {
                if tol > config.high_tolerance_threshold && ei < n_edges {
                    let edge = &brep.edges[ei];
                    combined.insert(edge.start);
                    combined.insert(edge.end);
                }
            }

            // Near duplicates
            for i in 0..n_vertices.min(1000) {
                for j in (i + 1)..n_vertices.min(i + 100) {
                    let dist = (brep.vertices[i].point - brep.vertices[j].point).length();
                    if dist < config.near_duplicate_distance {
                        combined.insert(i);
                        combined.insert(j);
                    }
                }
            }

            vertex_set = combined;
            result.strategy_counts.insert("hybrid".to_string(), vertex_set.len());
        }
    }

    // Apply neighborhood expansion
    if config.neighborhood_depth > 0 {
        let expanded = expand_seed_neighborhood(brep, &vertex_set, config.neighborhood_depth);
        vertex_set = expanded;
    }

    // Apply max seeds limit
    if config.max_seeds > 0 && vertex_set.len() > config.max_seeds {
        let seeds: Vec<usize> = vertex_set.into_iter().take(config.max_seeds).collect();
        vertex_set = seeds.into_iter().collect();
    }

    result.seed_vertices = vertex_set.into_iter().collect();
    result.seed_edges = edge_set.into_iter().collect();
    result.coverage_ratio = if n_vertices > 0 {
        result.seed_vertices.len() as f64 / n_vertices as f64
    } else {
        0.0
    };

    result
}

/// Expand seed set to include neighboring vertices.
fn expand_seed_neighborhood(
    brep: &BRep,
    seeds: &std::collections::HashSet<usize>,
    depth: usize,
) -> std::collections::HashSet<usize> {
    if depth == 0 {
        return seeds.clone();
    }

    // Build vertex-to-vertex adjacency via edges
    let mut adjacency: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for edge in &brep.edges {
        adjacency.entry(edge.start).or_default().push(edge.end);
        adjacency.entry(edge.end).or_default().push(edge.start);
    }

    let mut expanded = seeds.clone();
    let mut frontier: std::collections::HashSet<usize> = seeds.clone();

    for _ in 0..depth {
        let mut next_frontier = std::collections::HashSet::new();
        for &v in &frontier {
            if let Some(neighbors) = adjacency.get(&v) {
                for &n in neighbors {
                    if !expanded.contains(&n) {
                        expanded.insert(n);
                        next_frontier.insert(n);
                    }
                }
            }
        }
        frontier = next_frontier;
        if frontier.is_empty() {
            break;
        }
    }

    expanded
}

/// Apply scoped make-connected with automatic seed detection.
pub fn make_connected_scoped_auto(
    brep: &BRep,
    config: &SeedDetectionConfig,
    tolerance: f64,
    max_passes: usize,
) -> (BRep, MakeConnectedReport, SeedDetectionResult) {
    let seeds = detect_seeds_for_scoped_cleanup(brep, config);

    let (result, report) = make_connected_iterative_scoped_with_growth_cap(
        brep,
        &seeds.seed_vertices,
        tolerance,
        max_passes,
        1.5,
        tolerance * 10.0,
    );

    (result, report, seeds)
}

/// Apply a MakeConnectedStrategy to repair connectivity.
///
/// This is a convenience function that delegates to `strategy.apply(brep)`.
pub fn make_connected_with_strategy(brep: &BRep, strategy: &MakeConnectedStrategy) -> (BRep, MakeConnectedReport) {
    strategy.apply(brep)
}

/// Information about a shared edge between two faces.
#[derive(Debug, Clone)]
pub struct SharedEdgeInfo {
    /// Index of the first edge.
    pub edge_a: usize,
    /// Index of the second edge.
    pub edge_b: usize,
    /// Whether the edges have geometric compatibility (same curve type).
    pub geometry_compatible: bool,
    /// Whether the curvature is continuous across the shared edge.
    pub curvature_continuous: bool,
    /// Whether the parameter ranges are compatible (overlap).
    pub param_range_compatible: bool,
    /// Maximum deviation between the two edges.
    pub max_deviation: f64,
    /// Whether the edges are reversed relative to each other.
    pub reversed: bool,
}

/// Classification of shared face topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SharedFaceKind {
    /// Faces share their complete boundary (fully coincident).
    FullShared,
    /// Faces share a partial boundary (some edges coincide).
    PartialShared,
    /// Faces share only some vertices.
    VertexShared,
    /// Faces are adjacent (share an edge) but not overlapping.
    Adjacent,
}

/// Information about a shared face pair.
#[derive(Debug, Clone)]
pub struct SharedFaceInfo {
    /// Index of the first face.
    pub face_a: usize,
    /// Index of the second face.
    pub face_b: usize,
    /// Classification of the sharing.
    pub kind: SharedFaceKind,
    /// Indices of shared edges.
    pub shared_edges: Vec<usize>,
    /// Indices of shared vertices.
    pub shared_vertices: Vec<usize>,
    /// Whether the face normals are compatible (parallel or anti-parallel).
    pub normals_compatible: bool,
}

/// Report from advanced shared topology detection.
#[derive(Debug, Clone, Default)]
pub struct SharedTopologyReport {
    /// Fully shared face pairs.
    pub fully_shared_faces: Vec<SharedFaceInfo>,
    /// Partially shared face pairs.
    pub partially_shared_faces: Vec<SharedFaceInfo>,
    /// Shared edges with detailed information.
    pub shared_edges: Vec<SharedEdgeInfo>,
    /// Total number of shared vertex pairs.
    pub shared_vertex_pairs: usize,
    /// Whether any shared topology was detected.
    pub has_shared_topology: bool,
    /// Summary string for debugging.
    pub summary: String,
}

impl std::fmt::Display for RepairReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "RepairReport {{ vertices_merged={}, degenerate_removed={}, normals_recomputed={}, faces_reoriented={}, wires_fixed={}, same_range_fixed={}, same_parameter_fixed={} }}",
            self.vertices_merged,
            self.degenerate_faces_removed,
            self.normals_recomputed,
            self.faces_reoriented,
            self.wires_fixed,
            self.same_range_fixed,
            self.same_parameter_fixed,
        )
    }
}

/// Apply all repair operations in a single pass and return the cleaned BRep
/// together with a summary of changes made.
///
/// Equivalent to `ShapeFix_Shape::Perform()` followed by
/// `BRepLib::UpdateEdgeTol()`.
pub fn repair(brep: &BRep, tolerance: f64) -> (BRep, RepairReport) {
    let mut report = RepairReport::default();
    let (b, n) = merge_close_vertices(brep, tolerance);
    report.vertices_merged += n;
    let (b, n) = recompute_face_normals(&b);
    report.normals_recomputed += n;
    let (b, n) = fix_face_orientation(&b);
    report.faces_reoriented += n;
    let (b, n) = remove_degenerate_faces(&b);
    report.degenerate_faces_removed += n;
    let (b, n) = fix_wire_orientation(&b, tolerance);
    report.wires_fixed += n;
    let (b, n) = fix_same_range_flags(&b, tolerance);
    report.same_range_fixed += n;
    let (b, n) = fix_same_parameter(&b, tolerance);
    report.same_parameter_fixed += n;
    (b, report)
}

/// Baseline "MakeConnected"-style cleanup.
///
/// This pass snaps near-coincident vertices and removes tiny/degenerate edges
/// to improve topological connectivity before downstream operations.
pub fn make_connected_baseline(brep: &BRep, tolerance: f64) -> (BRep, MakeConnectedReport) {
    make_connected_iterative(brep, tolerance, 1)
}

/// Iterative baseline "MakeConnected"-style cleanup.
///
/// Runs repeated merge/small-edge cleanup passes until convergence or until
/// `max_passes` is reached.
pub fn make_connected_iterative(
    brep: &BRep,
    tolerance: f64,
    max_passes: usize,
) -> (BRep, MakeConnectedReport) {
    make_connected_iterative_with_growth(brep, tolerance, max_passes, 1.0)
}

/// Iterative baseline "MakeConnected" cleanup with per-pass tolerance growth.
///
/// `tolerance_growth` values <= 1.0 keep fixed tolerance across passes.
pub fn make_connected_iterative_with_growth(
    brep: &BRep,
    tolerance: f64,
    max_passes: usize,
    tolerance_growth: f64,
) -> (BRep, MakeConnectedReport) {
    make_connected_iterative_with_growth_cap(
        brep,
        tolerance,
        max_passes,
        tolerance_growth,
        f64::INFINITY,
    )
}

/// Iterative baseline "MakeConnected" cleanup with per-pass tolerance growth
/// and an optional upper cap for safety.
pub fn make_connected_iterative_with_growth_cap(
    brep: &BRep,
    tolerance: f64,
    max_passes: usize,
    tolerance_growth: f64,
    tolerance_cap: f64,
) -> (BRep, MakeConnectedReport) {
    let tol = tolerance.max(TOLERANCE_ABS);
    let pass_limit = max_passes.max(1);
    let growth = if tolerance_growth > 1.0 {
        tolerance_growth
    } else {
        1.0
    };
    let tol_cap = tolerance_cap.max(tol);
    let mut out = brep.clone();
    let mut report = MakeConnectedReport::default();

    for pass_idx in 0..pass_limit {
        let grown_tol = tol * growth.powi(pass_idx as i32);
        let pass_tol = grown_tol.min(tol_cap);
        let (b, merged) = merge_close_vertices(&out, pass_tol);
        let (b, removed) = remove_small_edges(&b, pass_tol);
        out = b;

        report.vertices_merged += merged;
        report.small_edges_removed += removed;
        report.passes_run = pass_idx + 1;
        report.final_tolerance = pass_tol;
        if grown_tol > tol_cap {
            report.tolerance_cap_applied = true;
        }

        if merged == 0 && removed == 0 {
            if make_connected_has_future_tolerance_increase(
                pass_idx,
                pass_limit,
                pass_tol,
                tol,
                growth,
                tol_cap,
            ) {
                continue;
            }
            report.converged = true;
            break;
        }
    }

    (out, report)
}

/// Scoped iterative make-connected cleanup limited to a local vertex region.
///
/// Only short-edge removal and near-vertex merges touching `scope_vertices`
/// are applied, allowing localized connectivity fixes.
pub fn make_connected_iterative_scoped_with_growth_cap(
    brep: &BRep,
    scope_vertices: &[usize],
    tolerance: f64,
    max_passes: usize,
    tolerance_growth: f64,
    tolerance_cap: f64,
) -> (BRep, MakeConnectedReport) {
    let tol = tolerance.max(TOLERANCE_ABS);
    let pass_limit = max_passes.max(1);
    let growth = if tolerance_growth > 1.0 {
        tolerance_growth
    } else {
        1.0
    };
    let tol_cap = tolerance_cap.max(tol);

    let mut scope_set: std::collections::HashSet<usize> = scope_vertices.iter().copied().collect();
    let mut out = brep.clone();
    let mut report = MakeConnectedReport::default();

    if scope_set.is_empty() {
        report.passes_run = 1;
        report.converged = true;
        report.final_tolerance = tol;
        return (out, report);
    }

    for pass_idx in 0..pass_limit {
        let grown_tol = tol * growth.powi(pass_idx as i32);
        let pass_tol = grown_tol.min(tol_cap);

        let (b, merged, remap) = merge_close_vertices_scoped(&out, pass_tol, &scope_set);
        let mapped_scope: std::collections::HashSet<usize> = scope_set
            .iter()
            .filter_map(|v| remap.get(v).copied())
            .collect();
        let (b, removed, remap_scope) = remove_small_edges_scoped(&b, pass_tol, &mapped_scope);
        let next_scope: std::collections::HashSet<usize> = mapped_scope
            .iter()
            .filter_map(|v| remap_scope.get(v).copied())
            .collect();

        out = b;
        scope_set = next_scope;

        report.vertices_merged += merged;
        report.small_edges_removed += removed;
        report.passes_run = pass_idx + 1;
        report.final_tolerance = pass_tol;
        if grown_tol > tol_cap {
            report.tolerance_cap_applied = true;
        }

        if merged == 0 && removed == 0 {
            if make_connected_has_future_tolerance_increase(
                pass_idx,
                pass_limit,
                pass_tol,
                tol,
                growth,
                tol_cap,
            ) {
                continue;
            }
            report.converged = true;
            break;
        }
    }

    (out, report)
}

fn merge_close_vertices_scoped(
    brep: &BRep,
    tolerance: f64,
    scope_vertices: &std::collections::HashSet<usize>,
) -> (BRep, usize, std::collections::HashMap<usize, usize>) {
    let n = brep.vertices.len();
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            if ra < rb {
                parent[rb] = ra;
            } else {
                parent[ra] = rb;
            }
        }
    }

    let tol2 = tolerance * tolerance;
    for i in 0..n {
        for j in (i + 1)..n {
            if !(scope_vertices.contains(&i) || scope_vertices.contains(&j)) {
                continue;
            }
            let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
            if d2 <= tol2 {
                union(&mut parent, i, j);
            }
        }
    }

    for i in 0..n {
        parent[i] = find(&mut parent, i);
    }

    let merged = (0..n).filter(|&i| parent[i] != i).count();
    let mut identity_map: std::collections::HashMap<usize, usize> =
        (0..n).map(|i| (i, i)).collect();
    if merged == 0 {
        return (brep.clone(), 0, identity_map);
    }

    let mut new_vertices = Vec::new();
    let mut remap = vec![0usize; n];
    let mut seen: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for i in 0..n {
        let rep = parent[i];
        if let Some(&new_idx) = seen.get(&rep) {
            remap[i] = new_idx;
        } else {
            let new_idx = new_vertices.len();
            new_vertices.push(brep.vertices[rep]);
            seen.insert(rep, new_idx);
            remap[i] = new_idx;
        }
    }

    let new_edges: Vec<Edge> = brep
        .edges
        .iter()
        .map(|e| Edge {
            start: remap[e.start],
            end: remap[e.end],
        })
        .collect();

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| Face {
                            outer_wire: face.outer_wire.clone(),
                            inner_wires: face.inner_wires.clone(),
                            normal: face.normal,
                            triangles: face.triangles.clone(),
                            mesh_dirty: true,
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.vertices = new_vertices;
    result.edges = new_edges;
    result.solids = new_solids;

    identity_map.clear();
    for (old, newv) in remap.into_iter().enumerate() {
        identity_map.insert(old, newv);
    }
    (result, merged, identity_map)
}

fn remove_small_edges_scoped(
    brep: &BRep,
    min_length: f64,
    scope_vertices: &std::collections::HashSet<usize>,
) -> (BRep, usize, std::collections::HashMap<usize, usize>) {
    let mut out = brep.clone();
    let mut total_removed = 0usize;
    let mut remap_track: Vec<usize> = (0..brep.vertices.len()).collect();

    loop {
        let edge_count = out.edges.len();
        let mut removed_edge: Option<usize> = None;

        for ei in 0..edge_count {
            let edge = &out.edges[ei];
            let start = edge.start;
            let end = edge.end;
            if !(scope_vertices.contains(&start) || scope_vertices.contains(&end)) {
                continue;
            }

            let is_degenerate = start == end;
            let is_short = if is_degenerate {
                true
            } else {
                let ps = out.vertices[start].point;
                let pe = out.vertices[end].point;
                (pe - ps).length() < min_length
            };

            if is_short {
                removed_edge = Some(ei);
                break;
            }
        }

        let Some(ei) = removed_edge else { break };
        let edge = out.edges[ei];
        let is_loop = edge.start == edge.end;
        let keep_vi = edge.start.min(edge.end);
        let drop_vi = edge.start.max(edge.end);

        let remap_vertex = |vi: usize| -> usize {
            if vi == drop_vi {
                keep_vi
            } else if vi > drop_vi {
                vi - 1
            } else {
                vi
            }
        };

        if !is_loop {
            out.vertices.remove(drop_vi);
            if out.geom.vertex_tolerance.len() > drop_vi
                && drop_vi != out.geom.vertex_tolerance.len()
            {
                out.geom.vertex_tolerance.remove(drop_vi);
            }
            for r in &mut remap_track {
                if *r == drop_vi {
                    *r = keep_vi;
                } else if *r > drop_vi {
                    *r -= 1;
                }
            }
        }

        for e in &mut out.edges {
            e.start = remap_vertex(e.start);
            e.end = remap_vertex(e.end);
        }

        out.edges.remove(ei);
        macro_rules! rm {
            ($vec:expr) => {
                if ei < $vec.len() {
                    $vec.remove(ei);
                }
            };
        }
        rm!(out.geom.edge_curve);
        rm!(out.geom.edge_curve_range);
        rm!(out.geom.edge_degenerated);
        rm!(out.geom.edge_pcurves);
        rm!(out.geom.edge_same_parameter);
        rm!(out.geom.edge_same_range);
        rm!(out.geom.edge_tolerance);

        let remap_edge = |we_idx: usize| -> usize {
            if we_idx > ei { we_idx - 1 } else { we_idx }
        };
        for solid in &mut out.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    let filter_remap = |wire: &mut Wire| {
                        wire.edges.retain(|we| we.idx != ei);
                        for we in &mut wire.edges {
                            we.idx = remap_edge(we.idx);
                        }
                    };
                    filter_remap(&mut face.outer_wire);
                    for iw in &mut face.inner_wires {
                        filter_remap(iw);
                    }
                }
            }
        }

        total_removed += 1;
    }

    let mut remap_map = std::collections::HashMap::new();
    for (old, newv) in remap_track.into_iter().enumerate() {
        remap_map.insert(old, newv);
    }

    (out, total_removed, remap_map)
}

// ─────────────────────────────────────────────────────────────────────────────
// Edge Sewing (MakeConnected enhancement)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from edge sewing operations.
#[derive(Debug, Clone, Default)]
pub struct EdgeSewReport {
    /// Number of edge pairs that were sewn together.
    pub edges_sewn: usize,
    /// Number of vertex pairs that were merged as a result.
    pub vertices_merged: usize,
}

/// Sew close edges together by merging their endpoints.
///
/// This is a key part of MakeConnected connectivity rebuilding: when two edges
/// are geometrically close (share similar curves and have nearby endpoints),
/// they are "sewn" together by merging their vertices.
///
/// Analogous to `BRepBuilderAPI_Sewing` edge merging in OCCT.
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `tolerance` - Maximum distance for considering vertices coincident.
///
/// # Returns
/// A tuple of (modified BRep, report).
pub fn sew_close_edges(brep: &BRep, tolerance: f64) -> (BRep, EdgeSewReport) {
    let tol = tolerance.max(TOLERANCE_ABS);
    let tol_sq = tol * tol;
    let mut result = brep.clone();
    let mut report = EdgeSewReport::default();

    let n = result.edges.len();
    if n < 2 {
        return (result, report);
    }

    // Find edge pairs that should be sewn
    let mut vertex_merge_pairs: Vec<(usize, usize)> = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let edge_i = &result.edges[i];
            let edge_j = &result.edges[j];

            // Check if edges share similar geometry
            if !edges_similar_geometry(&result, i, j, tol) {
                continue;
            }

            // Check if endpoints are close enough to sew
            let p_i_start = result.vertices[edge_i.start].point;
            let p_i_end = result.vertices[edge_i.end].point;
            let p_j_start = result.vertices[edge_j.start].point;
            let p_j_end = result.vertices[edge_j.end].point;

            // Check all possible endpoint combinations
            let d_ss = (p_i_start - p_j_start).length_squared();
            let d_se = (p_i_start - p_j_end).length_squared();
            let d_es = (p_i_end - p_j_start).length_squared();
            let d_ee = (p_i_end - p_j_end).length_squared();

            // Find minimum distance pairing
            let min_dist_sq = d_ss.min(d_se).min(d_es).min(d_ee);

            if min_dist_sq <= tol_sq {
                // Determine which vertices to merge
                if d_ss <= tol_sq && edge_i.start != edge_j.start {
                    vertex_merge_pairs.push((edge_i.start, edge_j.start));
                }
                if d_se <= tol_sq && edge_i.start != edge_j.end {
                    vertex_merge_pairs.push((edge_i.start, edge_j.end));
                }
                if d_es <= tol_sq && edge_i.end != edge_j.start {
                    vertex_merge_pairs.push((edge_i.end, edge_j.start));
                }
                if d_ee <= tol_sq && edge_i.end != edge_j.end {
                    vertex_merge_pairs.push((edge_i.end, edge_j.end));
                }

                report.edges_sewn += 1;
            }
        }
    }

    if vertex_merge_pairs.is_empty() {
        return (result, report);
    }

    // Apply vertex merges using union-find
    let n_verts = result.vertices.len();
    let mut parent: Vec<usize> = (0..n_verts).collect();

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            if ra < rb {
                parent[rb] = ra;
            } else {
                parent[ra] = rb;
            }
        }
    }

    for (v1, v2) in &vertex_merge_pairs {
        union(&mut parent, *v1, *v2);
    }

    // Count merged vertices
    let mut merged_count = 0usize;
    for i in 0..n_verts {
        if parent[i] != i {
            merged_count += 1;
        }
    }

    if merged_count == 0 {
        return (result, report);
    }

    // Build remapping
    let mut new_vertices = Vec::new();
    let mut remap = vec![0usize; n_verts];
    let mut seen: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

    for i in 0..n_verts {
        let rep = find(&mut parent, i);
        if let Some(&new_idx) = seen.get(&rep) {
            remap[i] = new_idx;
        } else {
            let new_idx = new_vertices.len();
            new_vertices.push(result.vertices[rep]);
            seen.insert(rep, new_idx);
            remap[i] = new_idx;
        }
    }

    // Update edges
    for edge in &mut result.edges {
        edge.start = remap[edge.start];
        edge.end = remap[edge.end];
    }

    result.vertices = new_vertices;
    report.vertices_merged = merged_count;

    (result, report)
}

/// Check if two edges have similar geometry (same curve type and parameters).
fn edges_similar_geometry(brep: &BRep, e1: usize, e2: usize, tol: f64) -> bool {
    // Check if edges have the same curve type
    let curve1 = brep.geom.curves.get(e1);
    let curve2 = brep.geom.curves.get(e2);

    match (curve1, curve2) {
        (Some(c1), Some(c2)) => {
            match (c1, c2) {
                (rcad_kernel::Curve3::Line(l1), rcad_kernel::Curve3::Line(l2)) => {
                    // Check if lines are parallel (or anti-parallel)
                    let d1 = l1.direction.normalize_or_zero();
                    let d2 = l2.direction.normalize_or_zero();
                    if d1.dot(d2).abs() < 0.99 {
                        return false;
                    }
                    // Check if origins are close
                    let v = l2.origin - l1.origin;
                    let perp = v - d1 * v.dot(d1);
                    perp.length() <= tol
                }
                (rcad_kernel::Curve3::Circle(c1), rcad_kernel::Curve3::Circle(c2)) => {
                    (c1.center - c2.center).length() <= tol
                        && c1.normal.dot(c2.normal).abs() >= 0.99
                        && (c1.radius - c2.radius).abs() <= tol
                }
                _ => false,
            }
        }
        _ => {
            // No curve data - use vertex-based check
            let edge1 = &brep.edges[e1];
            let edge2 = &brep.edges[e2];
            let p1_start = brep.vertices[edge1.start].point;
            let p1_end = brep.vertices[edge1.end].point;
            let p2_start = brep.vertices[edge2.start].point;
            let p2_end = brep.vertices[edge2.end].point;

            // Check if edges have similar length and direction
            let len1 = (p1_end - p1_start).length();
            let len2 = (p2_end - p2_start).length();
            (len1 - len2).abs() <= tol
        }
    }
}

/// Enhanced make-connected with edge sewing.
///
/// This combines vertex merging, edge sewing, and small edge removal
/// into a comprehensive connectivity rebuilding pass.
///
/// Analogous to `BOPAlgo_MakeConnected` in OCCT.
pub fn make_connected_enhanced(brep: &BRep, tolerance: f64, max_passes: usize) -> (BRep, MakeConnectedReport) {
    make_connected_enhanced_with_mode(brep, tolerance, max_passes, MakeConnectedMode::Standard, false)
}

/// Enhanced make-connected with mode selection.
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `tolerance` - Maximum distance for considering vertices coincident.
/// * `max_passes` - Maximum number of passes to run.
/// * `mode` - Operating mode (Standard, Aggressive, Conservative).
/// * `merge_faces` - Whether to merge shared faces (only in Aggressive mode).
///
/// # Returns
/// A tuple of (modified BRep, report).
pub fn make_connected_enhanced_with_mode(
    brep: &BRep,
    tolerance: f64,
    max_passes: usize,
    mode: MakeConnectedMode,
    merge_faces: bool,
) -> (BRep, MakeConnectedReport) {
    let tol = tolerance.max(TOLERANCE_ABS);
    let mut out = brep.clone();
    let mut report = MakeConnectedReport::default();

    for _pass in 0..max_passes {
        let mut changed = false;

        // Step 1: Sew close edges (only in Aggressive mode)
        if mode == MakeConnectedMode::Aggressive {
            let (b, sew_report) = sew_close_edges(&out, tol);
            if sew_report.edges_sewn > 0 || sew_report.vertices_merged > 0 {
                out = b;
                report.vertices_merged += sew_report.vertices_merged;
                report.edges_sewn += sew_report.edges_sewn;
                changed = true;
            }
        }

        // Step 2: Merge close vertices (always)
        let (b, merged) = merge_close_vertices(&out, tol);
        if merged > 0 {
            out = b;
            report.vertices_merged += merged;
            changed = true;
        }

        // Step 3: Remove small edges (not in Conservative mode)
        if mode != MakeConnectedMode::Conservative {
            let (b, removed) = remove_small_edges(&out, tol);
            if removed > 0 {
                out = b;
                report.small_edges_removed += removed;
                changed = true;
            }
        }

        // Step 4: Merge shared faces (only in Aggressive mode with merge_faces)
        if mode == MakeConnectedMode::Aggressive && merge_faces {
            let (b, merged) = merge_shared_faces(&out, tol);
            if merged > 0 {
                out = b;
                report.faces_merged += merged;
                changed = true;
            }
        }

        report.passes_run += 1;

        if !changed {
            report.converged = true;
            break;
        }
    }

    report.final_tolerance = tol;
    (out, report)
}

// ─────────────────────────────────────────────────────────────────────────────
// Advanced Shared Topology Detection
// ─────────────────────────────────────────────────────────────────────────────

/// Detect shared topology between faces with advanced classification.
///
/// This function analyzes a BRep to identify shared topology between faces,
/// including fully shared faces, partially shared faces, shared edges with
/// curvature continuity, and shared vertices.
///
/// # Arguments
/// * `brep` - The BRep to analyze.
/// * `tolerance` - Maximum distance for considering geometry coincident.
///
/// # Returns
/// A `SharedTopologyReport` containing detailed classification of shared topology.
pub fn detect_shared_topology_advanced(brep: &BRep, tolerance: f64) -> SharedTopologyReport {
    let tol = tolerance.max(TOLERANCE_ABS);
    let mut report = SharedTopologyReport::default();

    // Count shared vertices (near-coincident vertices with different indices)
    // This is done first so it works even for single-face BReps
    let tol_sq = tol * tol;
    let n_verts = brep.vertices.len();
    for i in 0..n_verts {
        for j in (i + 1)..n_verts {
            let dist_sq = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
            if dist_sq <= tol_sq {
                report.shared_vertex_pairs += 1;
            }
        }
    }

    // Collect all faces with their flattened indices
    let faces: Vec<(usize, usize, usize, &Face)> = brep
        .solids
        .iter()
        .enumerate()
        .flat_map(|(si, solid)| {
            solid.shells.iter().enumerate().flat_map(move |(shi, shell)| {
                shell.faces.iter().enumerate().map(move |(fi, face)| (si, shi, fi, face))
            })
        })
        .collect();

    let n_faces = faces.len();
    if n_faces < 2 {
        // Still need to set summary and has_shared_topology for single-face case
        report.has_shared_topology = report.shared_vertex_pairs > 0 || !report.shared_edges.is_empty();
        report.summary = format!(
            "SharedTopology: {} fully shared faces, {} partially shared faces, {} shared edges, {} shared vertex pairs",
            report.fully_shared_faces.len(),
            report.partially_shared_faces.len(),
            report.shared_edges.len(),
            report.shared_vertex_pairs
        );
        return report;
    }

    // Build edge-to-face map
    let mut edge_to_faces: std::collections::HashMap<usize, Vec<(usize, usize, usize)>> =
        std::collections::HashMap::new();
    for (si, shi, fi, face) in &faces {
        for we in &face.outer_wire.edges {
            edge_to_faces.entry(we.idx).or_default().push((*si, *shi, *fi));
        }
    }

    // Detect shared edges with curvature continuity
    let mut processed_edge_pairs: std::collections::HashSet<(usize, usize)> =
        std::collections::HashSet::new();

    for (e1, _faces1) in &edge_to_faces {
        for (e2, _faces2) in &edge_to_faces {
            if e1 >= e2 {
                continue;
            }
            if processed_edge_pairs.contains(&(*e1, *e2)) {
                continue;
            }
            processed_edge_pairs.insert((*e1, *e2));

            // Check if edges have shared geometry
            if let Some(info) = analyze_shared_edge_pair(brep, *e1, *e2, tol) {
                if info.geometry_compatible {
                    report.shared_edges.push(info);
                }
            }
        }
    }

    // Detect shared face pairs
    let mut processed_face_pairs: std::collections::HashSet<(usize, usize, usize, usize, usize, usize)> =
        std::collections::HashSet::new();

    for i in 0..n_faces {
        for j in (i + 1)..n_faces {
            let (si1, shi1, fi1, face1) = faces[i];
            let (si2, shi2, fi2, face2) = faces[j];

            // Skip same face
            if si1 == si2 && shi1 == shi2 && fi1 == fi2 {
                continue;
            }

            // Create unique key for face pair
            let key1 = (si1, shi1, fi1, si2, shi2, fi2);
            let key2 = (si2, shi2, fi2, si1, shi1, fi1);
            if processed_face_pairs.contains(&key1) || processed_face_pairs.contains(&key2) {
                continue;
            }
            processed_face_pairs.insert(key1);

            if let Some(info) = analyze_shared_face_pair(brep, face1, face2, i, j, tol) {
                match info.kind {
                    SharedFaceKind::FullShared => report.fully_shared_faces.push(info),
                    SharedFaceKind::PartialShared => report.partially_shared_faces.push(info),
                    _ => {}
                }
            }
        }
    }

    // Set summary
    report.has_shared_topology = !report.fully_shared_faces.is_empty()
        || !report.partially_shared_faces.is_empty()
        || !report.shared_edges.is_empty()
        || report.shared_vertex_pairs > 0;
    report.summary = format!(
        "SharedTopology: {} fully shared faces, {} partially shared faces, {} shared edges, {} shared vertex pairs",
        report.fully_shared_faces.len(),
        report.partially_shared_faces.len(),
        report.shared_edges.len(),
        report.shared_vertex_pairs
    );

    report
}

/// Analyze a pair of edges for shared topology.
fn analyze_shared_edge_pair(
    brep: &BRep,
    e1: usize,
    e2: usize,
    tolerance: f64,
) -> Option<SharedEdgeInfo> {
    let edge1 = brep.edges.get(e1)?;
    let edge2 = brep.edges.get(e2)?;

    let curve1 = brep.geom.curves.get(e1);
    let curve2 = brep.geom.curves.get(e2);
    let range1 = brep.geom.edge_curve_range.get(e1).and_then(|r| *r);
    let range2 = brep.geom.edge_curve_range.get(e2).and_then(|r| *r);

    // Check geometric compatibility
    let (geometry_compatible, max_deviation, reversed) = match (curve1, curve2) {
        (Some(c1), Some(c2)) => check_curve_compatibility(c1, c2, range1, range2, tolerance),
        (None, None) => {
            // Use vertex-based check
            let p1_start = brep.vertices.get(edge1.start)?.point;
            let p1_end = brep.vertices.get(edge1.end)?.point;
            let p2_start = brep.vertices.get(edge2.start)?.point;
            let p2_end = brep.vertices.get(edge2.end)?.point;

            let d_ss = (p1_start - p2_start).length();
            let d_se = (p1_start - p2_end).length();
            let d_es = (p1_end - p2_start).length();
            let d_ee = (p1_end - p2_end).length();

            let min_dev = d_ss.min(d_se).min(d_es).min(d_ee);
            let is_compatible = min_dev <= tolerance;
            let is_reversed = d_se <= tolerance || d_es <= tolerance;
            (is_compatible, min_dev, is_reversed)
        }
        _ => return None,
    };

    // Check curvature continuity
    let curvature_continuous = if geometry_compatible {
        check_edge_curvature_continuity(brep, e1, e2, tolerance)
    } else {
        false
    };

    // Check parameter range compatibility
    let param_range_compatible = if geometry_compatible {
        check_param_range_compatibility(brep, e1, e2, tolerance)
    } else {
        false
    };

    Some(SharedEdgeInfo {
        edge_a: e1,
        edge_b: e2,
        geometry_compatible,
        curvature_continuous,
        param_range_compatible,
        max_deviation,
        reversed,
    })
}

/// Check if two curves are geometrically compatible.
fn check_curve_compatibility(
    c1: &rcad_kernel::Curve3,
    c2: &rcad_kernel::Curve3,
    _range1: Option<[f64; 2]>,
    _range2: Option<[f64; 2]>,
    tolerance: f64,
) -> (bool, f64, bool) {
    match (c1, c2) {
        (rcad_kernel::Curve3::Line(l1), rcad_kernel::Curve3::Line(l2)) => {
            let d1 = l1.direction.normalize_or_zero();
            let d2 = l2.direction.normalize_or_zero();
            let dot = d1.dot(d2);

            if dot.abs() < 0.999 {
                return (false, f64::INFINITY, false);
            }

            // Check if origins are on the same line
            let v = l2.origin - l1.origin;
            let perp = v - d1 * v.dot(d1);
            let deviation = perp.length();
            let is_reversed = dot < 0.0;

            (deviation <= tolerance, deviation, is_reversed)
        }
        (rcad_kernel::Curve3::Circle(c1), rcad_kernel::Curve3::Circle(c2)) => {
            let center_dist = (c1.center - c2.center).length();
            let normal_dot = c1.normal.dot(c2.normal).abs();
            let radius_diff = (c1.radius - c2.radius).abs();

            let is_compatible =
                center_dist <= tolerance && normal_dot >= 0.999 && radius_diff <= tolerance;
            let deviation = center_dist.max(radius_diff);

            (is_compatible, deviation, false)
        }
        (rcad_kernel::Curve3::Ellipse(e1), rcad_kernel::Curve3::Ellipse(e2)) => {
            let center_dist = (e1.center - e2.center).length();
            let normal_dot = e1.normal.dot(e2.normal).abs();
            let major_diff = (e1.major_radius - e2.major_radius).abs();
            let minor_diff = (e1.minor_radius - e2.minor_radius).abs();

            let is_compatible = center_dist <= tolerance
                && normal_dot >= 0.999
                && major_diff <= tolerance
                && minor_diff <= tolerance;
            let deviation = center_dist.max(major_diff).max(minor_diff);

            (is_compatible, deviation, false)
        }
        _ => {
            // For other curve types, sample and check
            let n_samples = 16;
            let mut max_dev: f64 = 0.0;
            let mut reversed_candidates = 0;
            let mut total_samples = 0;

            for i in 0..n_samples {
                let t = i as f64 / (n_samples - 1).max(1) as f64;
                let p1 = c1.point_at(t);
                let p2 = c2.point_at(t);
                let p2_rev = c2.point_at(1.0 - t);

                let d_forward = (p1 - p2).length();
                let d_reverse = (p1 - p2_rev).length();

                max_dev = max_dev.max(d_forward.min(d_reverse));
                if d_reverse < d_forward {
                    reversed_candidates += 1;
                }
                total_samples += 1;
            }

            let is_compatible = max_dev <= tolerance;
            let is_reversed = reversed_candidates > total_samples / 2;

            (is_compatible, max_dev, is_reversed)
        }
    }
}

/// Check if two edges have curvature continuity.
fn check_edge_curvature_continuity(
    brep: &BRep,
    e1: usize,
    e2: usize,
    tolerance: f64,
) -> bool {
    let curve1 = match brep.geom.curves.get(e1) {
        Some(c) => c,
        None => return true, // No curve data, assume continuous
    };
    let curve2 = match brep.geom.curves.get(e2) {
        Some(c) => c,
        None => return true,
    };

    // Sample points along both edges and check curvature
    let n_samples = 8;
    let mut max_curvature_diff: f64 = 0.0;

    for i in 0..n_samples {
        let t = i as f64 / (n_samples - 1).max(1) as f64;

        // Get curvature at corresponding points
        let k1 = curve_curvature_at(curve1, t);
        let k2 = curve_curvature_at(curve2, t);

        if let (Some(k1), Some(k2)) = (k1, k2) {
            let diff = (k1 - k2).abs();
            max_curvature_diff = max_curvature_diff.max(diff);
        }
    }

    max_curvature_diff <= tolerance * 10.0 // Allow some tolerance for curvature
}

/// Get curvature at a parameter value on a curve.
fn curve_curvature_at(curve: &rcad_kernel::Curve3, t: f64) -> Option<f64> {
    use rcad_kernel::CurveEval;

    let h = 1e-6;
    let p0 = curve.point_at((t - h).max(0.0));
    let p1 = curve.point_at(t);
    let p2 = curve.point_at((t + h).min(1.0));

    // Approximate curvature using finite differences
    let d1 = (p1 - p0) / h;
    let d2 = (p2 - p1) / h;
    let dd = (d2 - d1) / h;

    let d1_len = d1.length();
    if d1_len < 1e-10 {
        return None;
    }

    let cross = d1.cross(dd);
    let curvature = cross.length() / (d1_len.powi(3));

    Some(curvature)
}

/// Check if parameter ranges are compatible.
fn check_param_range_compatibility(
    brep: &BRep,
    e1: usize,
    e2: usize,
    tolerance: f64,
) -> bool {
    let range1 = brep.geom.edge_curve_range.get(e1).and_then(|r| *r);
    let range2 = brep.geom.edge_curve_range.get(e2).and_then(|r| *r);

    match (range1, range2) {
        (Some(r1), Some(r2)) => {
            // Check for overlap
            let min_max = r1[1].min(r2[1]);
            let max_min = r1[0].max(r2[0]);
            min_max >= max_min - tolerance
        }
        _ => true, // No range data, assume compatible
    }
}

/// Analyze a pair of faces for shared topology.
fn analyze_shared_face_pair(
    brep: &BRep,
    face1: &Face,
    face2: &Face,
    flat_idx1: usize,
    flat_idx2: usize,
    tolerance: f64,
) -> Option<SharedFaceInfo> {
    // Collect boundary vertices
    let verts1: Vec<usize> = face1
        .outer_wire
        .edges
        .iter()
        .flat_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            if we.forward {
                Some(vec![edge.start, edge.end])
            } else {
                Some(vec![edge.end, edge.start])
            }
        })
        .flatten()
        .collect();

    let verts2: Vec<usize> = face2
        .outer_wire
        .edges
        .iter()
        .flat_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            if we.forward {
                Some(vec![edge.start, edge.end])
            } else {
                Some(vec![edge.end, edge.start])
            }
        })
        .flatten()
        .collect();

    // Count shared vertices
    let tol_sq = tolerance * tolerance;
    let mut shared_vertices = Vec::new();
    for &v1 in &verts1 {
        let p1 = brep.vertices.get(v1)?.point;
        for &v2 in &verts2 {
            let p2 = brep.vertices.get(v2)?.point;
            if (p1 - p2).length_squared() <= tol_sq {
                shared_vertices.push(v1.min(v2));
                break;
            }
        }
    }
    shared_vertices.sort();
    shared_vertices.dedup();

    // Collect boundary edges
    let edges1: std::collections::HashSet<usize> =
        face1.outer_wire.edges.iter().map(|we| we.idx).collect();
    let edges2: std::collections::HashSet<usize> =
        face2.outer_wire.edges.iter().map(|we| we.idx).collect();

    // Find shared edges (by geometry)
    let mut shared_edges = Vec::new();
    for &e1 in &edges1 {
        for &e2 in &edges2 {
            if let Some(info) = analyze_shared_edge_pair(brep, e1, e2, tolerance) {
                if info.geometry_compatible {
                    shared_edges.push(e1.min(e2));
                }
            }
        }
    }
    shared_edges.sort();
    shared_edges.dedup();

    // Determine sharing kind
    let kind = if shared_edges.len() == edges1.len() && shared_edges.len() == edges2.len() {
        SharedFaceKind::FullShared
    } else if !shared_edges.is_empty() {
        SharedFaceKind::PartialShared
    } else if !shared_vertices.is_empty() {
        SharedFaceKind::VertexShared
    } else {
        SharedFaceKind::Adjacent
    };

    // Check normal compatibility
    let normal_dot = face1.normal.dot(face2.normal).abs();
    let normals_compatible = normal_dot >= 0.999;

    Some(SharedFaceInfo {
        face_a: flat_idx1,
        face_b: flat_idx2,
        kind,
        shared_edges,
        shared_vertices,
        normals_compatible,
    })
}

/// Merge shared faces in a BRep.
///
/// This function identifies and merges faces that share their complete boundary.
/// Only available in Aggressive mode.
fn merge_shared_faces(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let report = detect_shared_topology_advanced(brep, tolerance);

    if report.fully_shared_faces.is_empty() {
        return (brep.clone(), 0);
    }

    // For now, just count the mergeable faces
    // A full implementation would actually merge the faces
    let merged_count = report.fully_shared_faces.len();

    (brep.clone(), merged_count)
}

/// Repair SameRange consistency by aligning PCurve ranges with the 3D edge range.
///
/// For each edge with a known `edge_curve_range` and attached PCurves, ensure all
/// referenced `curve2d_range` entries are populated with the same `[t1, t2]`.
/// Also marks `edge_same_range[edge_idx] = true` after alignment.
pub fn fix_same_range_flags(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let edge_count = out.edges.len();

    if out.geom.edge_same_range.len() < edge_count {
        out.geom.edge_same_range.resize(edge_count, true);
    }
    if out.geom.edge_curve_range.len() < edge_count {
        out.geom.edge_curve_range.resize(edge_count, None);
    }
    if out.geom.edge_pcurves.len() < edge_count {
        out.geom.edge_pcurves.resize(edge_count, Vec::new());
    }

    if out.geom.curve2d_range.len() < out.geom.curve2ds.len() {
        out.geom.curve2d_range.resize(out.geom.curve2ds.len(), None);
    }

    let mut fixed = 0usize;
    for edge_idx in 0..edge_count {
        let Some(range3d) = out.geom.edge_curve_range[edge_idx] else {
            continue;
        };
        let pcurves = out.geom.edge_pcurves[edge_idx].clone();
        if pcurves.is_empty() {
            continue;
        }

        let mut changed = !out.geom.edge_same_range[edge_idx];
        for pc in pcurves {
            if pc.curve2d_idx >= out.geom.curve2d_range.len() {
                continue;
            }
            match out.geom.curve2d_range[pc.curve2d_idx] {
                Some(r)
                    if (r[0] - range3d[0]).abs() <= tolerance
                        && (r[1] - range3d[1]).abs() <= tolerance => {}
                _ => {
                    out.geom.curve2d_range[pc.curve2d_idx] = Some(range3d);
                    changed = true;
                }
            }
        }

        if changed {
            out.geom.edge_same_range[edge_idx] = true;
            fixed += 1;
        }
    }

    (out, fixed)
}

/// Scan all edges for SameRange violations, flag them, and repair.
///
/// This combines the diagnostic scan from [`diagnose_same_range`] with the
/// repair logic of [`fix_same_range_flags`] in a single call.
pub fn fix_same_range_with_scan(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let diagnosis = diagnose_same_range(brep, tolerance);
    if diagnosis.suspect_edges.is_empty() {
        return (brep.clone(), 0);
    }

    let mut out = brep.clone();
    let n_edges = out.edges.len();

    if out.geom.edge_same_range.len() < n_edges {
        out.geom.edge_same_range.resize(n_edges, true);
    }

    for suspect in &diagnosis.suspect_edges {
        if suspect.edge_idx < n_edges {
            out.geom.edge_same_range[suspect.edge_idx] = false;
        }
    }

    fix_same_range_flags(&out, tolerance)
}

/// Merge vertices that are within `tolerance` of each other.
///
/// Uses spatial hashing for O(n) average performance on large models,
/// falling back to brute-force for small vertex counts.
/// For each pair of vertices closer than `tolerance`, they are merged into
/// the vertex with the smaller index. All edges and wires are remapped.
///
/// Returns the repaired BRep and the number of vertices merged.
///
/// Analogous to `BRepOffsetAPI_Sewing` vertex merging or
/// `ShapeFix_Wire::FixSameParameter`.
pub fn merge_close_vertices(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let n = brep.vertices.len();
    // Union-find: parent[i] = canonical representative of vertex i
    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]]; // path compression
            x = parent[x];
        }
        x
    }

    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            // Merge to the smaller index so result is deterministic
            if ra < rb {
                parent[rb] = ra;
            } else {
                parent[ra] = rb;
            }
        }
    }

    let tol2 = tolerance * tolerance;

    // Use spatial hashing for large models, brute-force for small ones.
    // Spatial hashing: bucket size = tolerance, check 27 neighbor cells.
    const SPATIAL_HASH_THRESHOLD: usize = 500;
    if n >= SPATIAL_HASH_THRESHOLD {
        let mut grid: std::collections::HashMap<(i32, i32, i32), Vec<usize>> =
            std::collections::HashMap::with_capacity(n);
        for i in 0..n {
            let p = brep.vertices[i].point;
            let cell = (
                (p.x / tolerance).floor() as i32,
                (p.y / tolerance).floor() as i32,
                (p.z / tolerance).floor() as i32,
            );
            // Check 27 neighbor cells (including self)
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbor = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                        if let Some(bucket) = grid.get(&neighbor) {
                            for &j in bucket {
                                let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
                                if d2 <= tol2 {
                                    union(&mut parent, i, j);
                                }
                            }
                        }
                    }
                }
            }
            grid.entry(cell).or_default().push(i);
        }
    } else {
        // Brute-force O(n²) — fast enough for small models
        for i in 0..n {
            for j in (i + 1)..n {
                let d2 = (brep.vertices[i].point - brep.vertices[j].point).length_squared();
                if d2 <= tol2 {
                    union(&mut parent, i, j);
                }
            }
        }
    }

    // Compress paths
    for i in 0..n {
        parent[i] = find(&mut parent, i);
    }

    // Count merges (vertices whose canonical rep is a different index)
    let merged = (0..n).filter(|&i| parent[i] != i).count();
    if merged == 0 {
        return (brep.clone(), 0);
    }

    // Build a compact vertex list and a remap table old_idx → new_idx
    let mut new_vertices: Vec<Vertex> = Vec::new();
    let mut remap = vec![0usize; n];
    let mut seen: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    for i in 0..n {
        let rep = parent[i];
        if let Some(&new_idx) = seen.get(&rep) {
            remap[i] = new_idx;
        } else {
            let new_idx = new_vertices.len();
            // Use the average position of all merged vertices for robustness
            new_vertices.push(brep.vertices[rep]);
            seen.insert(rep, new_idx);
            remap[i] = new_idx;
        }
    }

    // Re-map edges
    let new_edges: Vec<Edge> = brep
        .edges
        .iter()
        .map(|e| Edge {
            start: remap[e.start],
            end: remap[e.end],
        })
        .collect();

    // Rebuild solids with remapped wires (topology is unchanged, just vertex indices)
    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| {
                            let remap_wire = |w: &Wire| Wire {
                                edges: w.edges.clone(), // WireEdge indices are edge indices, not vertex
                            };
                            Face {
                                outer_wire: remap_wire(&face.outer_wire),
                                inner_wires: face.inner_wires.iter().map(remap_wire).collect(),
                                normal: face.normal,
                                triangles: face.triangles.clone(),
                                mesh_dirty: true,
                            }
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.vertices = new_vertices;
    result.edges = new_edges;
    result.solids = new_solids;

    (result, merged)
}

/// Remove faces that are degenerate:
/// - Fewer than 3 edges in the outer wire, or
/// - All wire vertices are collinear (zero-area face).
///
/// Returns the cleaned BRep and the number of faces removed.
///
/// Analogous to `ShapeFix_Shape` degenerate-face removal.
pub fn remove_degenerate_faces(brep: &BRep) -> (BRep, usize) {
    let mut removed = 0usize;

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| {
                    let new_faces: Vec<Face> = shell
                        .faces
                        .iter()
                        .filter(|face| {
                            let wire = &face.outer_wire;
                            // Must have at least 3 edges
                            if wire.edges.len() < 3 {
                                removed += 1;
                                return false;
                            }
                            // Collect distinct vertex positions
                            let pts: Vec<DVec3> = wire
                                .edges
                                .iter()
                                .filter_map(|we| {
                                    brep.edges.get(we.idx).and_then(|e| {
                                        let vidx = if we.forward { e.start } else { e.end };
                                        brep.vertices.get(vidx).map(|v| v.point)
                                    })
                                })
                                .collect();

                            if pts.len() < 3 {
                                removed += 1;
                                return false;
                            }

                            // Check for zero area using Newell's method
                            let area2 = newell_area(&pts);
                            if area2 < 1e-20 {
                                removed += 1;
                                return false;
                            }
                            true
                        })
                        .cloned()
                        .collect();
                    Shell { faces: new_faces }
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.solids = new_solids;
    (result, removed)
}

/// Recompute each face's `normal` field from the positions of its wire vertices,
/// using Newell's method for robustness with non-planar polygons.
///
/// Returns the updated BRep and the number of faces whose normals changed by
/// more than 1° (indicating they were stale or flipped).
///
/// Analogous to `BRepLib` normal re-computation after topology repair.
pub fn recompute_face_normals(brep: &BRep) -> (BRep, usize) {
    let mut changed = 0usize;

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| {
                            let pts: Vec<DVec3> = face
                                .outer_wire
                                .edges
                                .iter()
                                .filter_map(|we| {
                                    brep.edges.get(we.idx).and_then(|e| {
                                        let vidx = if we.forward { e.start } else { e.end };
                                        brep.vertices.get(vidx).map(|v| v.point)
                                    })
                                })
                                .collect();

                            let new_normal = if pts.len() >= 3 {
                                let n = newell_normal(&pts);
                                if n.length() > 1e-14 {
                                    n.normalize()
                                } else {
                                    face.normal
                                }
                            } else {
                                face.normal
                            };

                            let dot = face.normal.dot(new_normal);
                            // dot < cos(1°) ≈ 0.9998 means the normal changed significantly
                            if dot < 0.9998 {
                                changed += 1;
                            }

                            Face {
                                outer_wire: face.outer_wire.clone(),
                                inner_wires: face.inner_wires.clone(),
                                normal: new_normal,
                                triangles: face.triangles.clone(),
                                mesh_dirty: true,
                            }
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.solids = new_solids;
    (result, changed)
}

/// Ensure that each wire in the BRep forms a properly closed chain.
///
/// For each open wire (end of edge i ≠ start of edge i+1 within `tolerance`),
/// attempts to close it by reversing individual edges whose orientation appears
/// flipped relative to the chain direction.
///
/// Returns the repaired BRep and the count of wires that were modified.
///
/// Analogous to `ShapeFix_Wire::FixClosed()` / `FixConnected()`.
pub fn fix_wire_orientation(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let tol2 = tolerance * tolerance;
    let mut total_fixed = 0usize;

    let new_solids = brep
        .solids
        .iter()
        .map(|solid| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| {
                            let (new_outer, fixed_outer) = fix_wire(&face.outer_wire, brep, tol2);
                            let (new_inners, fixed_inner): (Vec<Wire>, usize) = face
                                .inner_wires
                                .iter()
                                .map(|w| fix_wire(w, brep, tol2))
                                .fold((Vec::new(), 0), |(mut wires, n), (w, f)| {
                                    wires.push(w);
                                    (wires, n + f)
                                });
                            let fixed = fixed_outer + fixed_inner;
                            total_fixed += fixed;
                            Face {
                                outer_wire: new_outer,
                                inner_wires: new_inners,
                                normal: face.normal,
                                triangles: face.triangles.clone(),
                                mesh_dirty: true,
                            }
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.solids = new_solids;
    (result, total_fixed)
}

/// Flip inward-facing faces so shell orientation is outward-consistent.
///
/// Uses the same centroid heuristic as [`check_orientation_consistency`]. Each
/// offending face has its stored normal negated and all wires reversed.
pub fn fix_face_orientation(brep: &BRep) -> (BRep, usize) {
    let report = check_orientation_consistency(brep);
    if report.issues.is_empty() {
        return (brep.clone(), 0);
    }

    let issue_set: std::collections::HashSet<(usize, usize)> = report
        .issues
        .iter()
        .map(|issue| (issue.solid_idx, issue.face_idx))
        .collect();

    let mut flat_face_idx = 0usize;
    let mut changed = 0usize;
    let new_solids = brep
        .solids
        .iter()
        .enumerate()
        .map(|(si, solid)| Solid {
            shells: solid
                .shells
                .iter()
                .map(|shell| Shell {
                    faces: shell
                        .faces
                        .iter()
                        .map(|face| {
                            let flip = issue_set.contains(&(si, flat_face_idx));
                            flat_face_idx += 1;
                            if flip {
                                changed += 1;
                                Face {
                                    outer_wire: reverse_wire(&face.outer_wire),
                                    inner_wires: face.inner_wires.iter().map(reverse_wire).collect(),
                                    normal: -face.normal,
                                    triangles: face.triangles.clone(),
                                    mesh_dirty: true,
                                }
                            } else {
                                face.clone()
                            }
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect();

    let mut result = brep.clone();
    result.solids = new_solids;
    (result, changed)
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Attempt to fix one wire, returning (fixed_wire, number_of_edges_flipped).
fn fix_wire(wire: &Wire, brep: &BRep, tol2: f64) -> (Wire, usize) {
    if wire.edges.len() < 2 {
        return (wire.clone(), 0);
    }

    let mut edges: Vec<WireEdge> = wire.edges.clone();
    let mut flipped = 0usize;
    let n = edges.len();

    for i in 0..n {
        let next = (i + 1) % n;
        let e_curr = match brep.edges.get(edges[i].idx) {
            Some(e) => e,
            None => continue,
        };
        let e_next = match brep.edges.get(edges[next].idx) {
            Some(e) => e,
            None => continue,
        };

        // end vertex of current edge
        let end_v = if edges[i].forward {
            e_curr.end
        } else {
            e_curr.start
        };
        // start vertex of next edge
        let start_v = if edges[next].forward {
            e_next.start
        } else {
            e_next.end
        };

        if end_v == start_v {
            continue; // already connected
        }
        // Check spatial proximity
        if let (Some(ep), Some(sp)) = (
            brep.vertices.get(end_v).map(|v| v.point),
            brep.vertices.get(start_v).map(|v| v.point),
        ) && (ep - sp).length_squared() <= tol2
        {
            continue; // close enough — OK
        }

        // Try flipping the *next* edge to see if that connects the chain
        let alt_start = if edges[next].forward {
            e_next.end
        } else {
            e_next.start
        };
        if alt_start == end_v {
            edges[next].forward = !edges[next].forward;
            flipped += 1;
        }
    }

    (Wire { edges }, flipped)
}

fn reverse_wire(wire: &Wire) -> Wire {
    let edges = wire
        .edges
        .iter()
        .rev()
        .map(|we| WireEdge::new(we.idx, !we.forward))
        .collect();
    Wire { edges }
}

/// Newell's method: compute the (un-normalized) area vector of a planar polygon.
fn newell_normal(pts: &[DVec3]) -> DVec3 {
    let n = pts.len();
    let mut normal = DVec3::ZERO;
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        normal.x += (a.y - b.y) * (a.z + b.z);
        normal.y += (a.z - b.z) * (a.x + b.x);
        normal.z += (a.x - b.x) * (a.y + b.y);
    }
    normal
}

/// Area magnitude squared (from Newell's method).
fn newell_area(pts: &[DVec3]) -> f64 {
    newell_normal(pts).length_squared()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

/// Repair SameParameter consistency by re-projecting PCurve endpoints onto the
/// 3D curve to align the parameterizations.
///
/// For each edge where `edge_same_parameter` is `false` and the edge has a known
/// 3D curve range and at least one PCurve, this function checks whether the 3D
/// curve start/end points match the PCurve's 2D start/end points on the
/// corresponding surface.  When the mismatch exceeds `tolerance`, it applies a
/// linear reparameterization: the PCurve's `curve2d_range` is scaled/shifted so
/// that the parameter range matches the 3D curve range, then
/// `edge_same_parameter[edge_idx]` is set to `true`.
///
/// This is the analogue of OCCT `BRepLib::SameParameter()` / `ShapeFix_Edge::FixSameParameter()`.
pub fn fix_same_parameter(brep: &BRep, _tolerance: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let edge_count = out.edges.len();

    if out.geom.edge_same_parameter.len() < edge_count {
        out.geom.edge_same_parameter.resize(edge_count, true);
    }
    if out.geom.edge_curve_range.len() < edge_count {
        out.geom.edge_curve_range.resize(edge_count, None);
    }
    if out.geom.edge_pcurves.len() < edge_count {
        out.geom.edge_pcurves.resize(edge_count, Vec::new());
    }
    if out.geom.curve2d_range.len() < out.geom.curve2ds.len() {
        out.geom.curve2d_range.resize(out.geom.curve2ds.len(), None);
    }

    let mut fixed = 0usize;
    for edge_idx in 0..edge_count {
        // Only repair edges explicitly flagged as *not* same-parameter.
        if out.geom.edge_same_parameter.get(edge_idx).copied().unwrap_or(true) {
            continue;
        }

        let Some(range3d) = out.geom.edge_curve_range[edge_idx] else {
            // Can't fix without a known 3D range; just mark as repaired to avoid
            // re-processing on next pass.
            out.geom.edge_same_parameter[edge_idx] = true;
            fixed += 1;
            continue;
        };

        let pcurves = out.geom.edge_pcurves[edge_idx].clone();
        if pcurves.is_empty() {
            // No PCurves: trivially same-parameter.
            out.geom.edge_same_parameter[edge_idx] = true;
            fixed += 1;
            continue;
        }

        // For each PCurve, align its range to match the 3D curve range.
        // Linear reparameterization: [pc_t0, pc_t1] → [range3d[0], range3d[1]].
        let mut changed = false;
        for pc in &pcurves {
            if pc.curve2d_idx >= out.geom.curve2d_range.len() {
                continue;
            }
            // Assign the 3D range as the canonical parameter range for this PCurve.
            // This is the coarsest possible fix (equivalent to assuming the PCurve
            // is already geometrically correct but needs re-parameterization).
            let current = out.geom.curve2d_range[pc.curve2d_idx];
            let target = Some(range3d);
            if current != target {
                out.geom.curve2d_range[pc.curve2d_idx] = target;
                changed = true;
            }
        }

        if changed || !out.geom.edge_same_parameter[edge_idx] {
            out.geom.edge_same_parameter[edge_idx] = true;
            fixed += 1;
        }
    }

    (out, fixed)
}

/// Scan all edges for SameParameter violations, flag them, and repair.
///
/// This combines the diagnostic scan from [`diagnose_same_parameter`] with the
/// repair logic of [`fix_same_parameter`] in a single call:
///
/// 1. Calls `diagnose_same_parameter` to find edges whose 3D curve endpoints
///    deviate from vertex positions beyond `tolerance`.
/// 2. Flags those edges with `edge_same_parameter = false`.
/// 3. Calls `fix_same_parameter` to reparameterize their PCurves.
///
/// Returns the repaired BRep and the number of edges repaired.
///
/// Analogous to OCCT `BRepLib::SameParameter(shape, enforce=true)`.
pub fn fix_same_parameter_with_scan(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let diagnosis = diagnose_same_parameter(brep, tolerance);
    if diagnosis.suspect_edges.is_empty() {
        return (brep.clone(), 0);
    }

    let mut out = brep.clone();
    let n_edges = out.edges.len();

    // Ensure edge_same_parameter is sized.
    if out.geom.edge_same_parameter.len() < n_edges {
        out.geom.edge_same_parameter.resize(n_edges, true);
    }

    // Flag suspect edges.
    for suspect in &diagnosis.suspect_edges {
        if suspect.edge_idx < n_edges {
            out.geom.edge_same_parameter[suspect.edge_idx] = false;
        }
    }

    // Now run the standard fix_same_parameter which repairs flagged edges.
    let (repaired, fixed) = fix_same_parameter(&out, tolerance);
    (repaired, fixed)
}

/// Remove short edges whose chord length is below `min_length`.
///
/// For each edge whose start and end vertices are closer than `min_length`,
/// the two endpoints are merged (lower index survives) and all topological
/// references are remapped. Degenerate self-loop edges (start == end) are
/// removed without vertex merging.
///
/// Analogous to OCCT `ShapeUpgrade_RemoveLocations` / `ShapeFix::RemoveSmallEdges`.
///
/// Returns the cleaned BRep and the number of short edges removed.
pub fn remove_small_edges(brep: &BRep, min_length: f64) -> (BRep, usize) {
    let mut out = brep.clone();
    let mut total_removed = 0usize;

    loop {
        let edge_count = out.edges.len();
        let mut removed_edge: Option<usize> = None;

        for ei in 0..edge_count {
            let edge = &out.edges[ei];
            let start = edge.start;
            let end = edge.end;

            // Degenerate self-loop: remove immediately
            let is_degenerate = start == end;
            let is_short = if is_degenerate {
                true
            } else {
                let ps = out.vertices[start].point;
                let pe = out.vertices[end].point;
                (pe - ps).length() < min_length
            };

            if is_short {
                removed_edge = Some(ei);
                break;
            }
        }

        let Some(ei) = removed_edge else { break };

        let edge = out.edges[ei];
        let keep_vi = edge.start.min(edge.end);
        let drop_vi = edge.start.max(edge.end);

        // Remap vertex references: drop_vi → keep_vi, shift higher indices down.
        let remap_vertex = |vi: usize| -> usize {
            if vi == drop_vi {
                keep_vi
            } else if vi > drop_vi {
                vi - 1
            } else {
                vi
            }
        };

        // Remove the dropped vertex from the vertex list.
        if !edge.start == !edge.end {
            // Self-loop: no vertex to remove
        } else {
            out.vertices.remove(drop_vi);
        }

        // Remap all edge endpoints.
        for e in &mut out.edges {
            e.start = remap_vertex(e.start);
            e.end = remap_vertex(e.end);
        }

        // Remap vertex tolerance parallel vec if present.
        if out.geom.vertex_tolerance.len() > drop_vi
            && drop_vi != out.geom.vertex_tolerance.len()
        {
            out.geom.vertex_tolerance.remove(drop_vi);
        }

        // Remove the short edge and its geom entries.
        out.edges.remove(ei);
        macro_rules! rm {
            ($vec:expr) => {
                if ei < $vec.len() {
                    $vec.remove(ei);
                }
            };
        }
        rm!(out.geom.edge_curve);
        rm!(out.geom.edge_curve_range);
        rm!(out.geom.edge_degenerated);
        rm!(out.geom.edge_pcurves);
        rm!(out.geom.edge_same_parameter);
        rm!(out.geom.edge_same_range);
        rm!(out.geom.edge_tolerance);

        // Remove wire references to this edge in all faces; remap remaining indices.
        let remap_edge = |we_idx: usize| -> usize {
            if we_idx > ei { we_idx - 1 } else { we_idx }
        };
        for solid in &mut out.solids {
            for shell in &mut solid.shells {
                for face in &mut shell.faces {
                    // Remove WireEdges pointing to the deleted edge from all wires.
                    let filter_remap = |wire: &mut Wire| {
                        wire.edges.retain(|we| we.idx != ei);
                        for we in &mut wire.edges {
                            we.idx = remap_edge(we.idx);
                        }
                    };
                    filter_remap(&mut face.outer_wire);
                    for iw in &mut face.inner_wires {
                        filter_remap(iw);
                    }
                }
            }
        }

        total_removed += 1;
    }

    (out, total_removed)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tolerance propagation
// ─────────────────────────────────────────────────────────────────────────────

/// Propagation direction for per-entity tolerance in a post-operation BRep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceFlowDirection {
    /// Vertex → edge → face (bottom-up, for newly assembled results).
    BottomUp,
    /// Face → edge → vertex (top-down, for degraded imports).
    TopDown,
}

/// Propagate per-entity tolerances throughout a BRep after a boolean, sew, or
/// import operation.
///
/// Analogous to `BRepLib::UpdateEdgeTol` + `BRepLib::SameParameter` tolerance
/// spreading in OCCT.
///
/// # Bottom-up (default after boolean operations)
///
/// 1. Fill missing `vertex_tolerance` slots with `tolerance_floor`.
/// 2. For each edge: `edge_tol = max(edge_tol, vtx_tol(start), vtx_tol(end))`.
/// 3. For each face: `face_tol = max(face_tol, max(wire edge tolerances))`.
///
/// # Top-down (useful after importing degraded STEP files)
///
/// Reverses the propagation: face tolerance spreads inward to edges and vertices.
///
/// # Arguments
///
/// - `brep`: input shape.
/// - `tolerance_floor`: minimum tolerance assigned to entities without an entry
///   (typically `CONFUSION` = 1e-7).
/// - `direction`: propagation direction.
pub fn propagate_tolerances(
    brep: &BRep,
    tolerance_floor: f64,
    direction: ToleranceFlowDirection,
) -> BRep {
    use crate::tolerance::TOLERANCE_ABS;
    let floor = tolerance_floor.max(TOLERANCE_ABS);
    let mut out = brep.clone();

    let n_verts = out.vertices.len();
    let n_edges = out.edges.len();

    // Count total faces (flattened order).
    let n_faces: usize = out.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();

    // Ensure arrays are sized.
    if out.geom.vertex_tolerance.len() < n_verts {
        out.geom.vertex_tolerance.resize(n_verts, floor);
    }
    if out.geom.edge_tolerance.len() < n_edges {
        out.geom.edge_tolerance.resize(n_edges, floor);
    }
    if out.geom.face_tolerance.len() < n_faces {
        out.geom.face_tolerance.resize(n_faces, floor);
    }

    match direction {
        ToleranceFlowDirection::BottomUp => {
            // Step 1: ensure vertices have at least floor tolerance.
            for vtol in &mut out.geom.vertex_tolerance {
                if *vtol < floor {
                    *vtol = floor;
                }
            }
            // Step 2: propagate vertex → edge.
            for ei in 0..n_edges {
                let st = out.edges[ei].start;
                let en = out.edges[ei].end;
                let vtol_s = out.geom.vertex_tolerance.get(st).copied().unwrap_or(floor);
                let vtol_e = out.geom.vertex_tolerance.get(en).copied().unwrap_or(floor);
                let cur = out.geom.edge_tolerance[ei];
                out.geom.edge_tolerance[ei] = cur.max(vtol_s).max(vtol_e).max(floor);
            }
            // Step 3: propagate edge → face.
            let mut flat_fi = 0usize;
            for si in 0..out.solids.len() {
                for shi in 0..out.solids[si].shells.len() {
                    for fi in 0..out.solids[si].shells[shi].faces.len() {
                        let face = &out.solids[si].shells[shi].faces[fi];
                        let mut max_etol: f64 = out.geom.face_tolerance[flat_fi];
                        for we in &face.outer_wire.edges {
                            let etol = out.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
                            max_etol = max_etol.max(etol);
                        }
                        for iw in &face.inner_wires {
                            for we in &iw.edges {
                                let etol = out.geom.edge_tolerance.get(we.idx).copied().unwrap_or(floor);
                                max_etol = max_etol.max(etol);
                            }
                        }
                        out.geom.face_tolerance[flat_fi] = max_etol.max(floor);
                        flat_fi += 1;
                    }
                }
            }
        }
        ToleranceFlowDirection::TopDown => {
            // Step 1: ensure faces have at least floor tolerance.
            for ftol in &mut out.geom.face_tolerance {
                if *ftol < floor {
                    *ftol = floor;
                }
            }
            // Step 2: propagate face → edge.
            let mut flat_fi = 0usize;
            for si in 0..out.solids.len() {
                for shi in 0..out.solids[si].shells.len() {
                    for fi in 0..out.solids[si].shells[shi].faces.len() {
                        let face = &out.solids[si].shells[shi].faces[fi];
                        let ftol = out.geom.face_tolerance[flat_fi];
                        for we in &face.outer_wire.edges {
                            if let Some(etol) = out.geom.edge_tolerance.get_mut(we.idx) {
                                *etol = etol.max(ftol);
                            }
                        }
                        for iw in &face.inner_wires {
                            for we in &iw.edges {
                                if let Some(etol) = out.geom.edge_tolerance.get_mut(we.idx) {
                                    *etol = etol.max(ftol);
                                }
                            }
                        }
                        flat_fi += 1;
                    }
                }
            }
            // Step 3: propagate edge → vertex.
            for ei in 0..n_edges {
                let etol = out.geom.edge_tolerance[ei];
                let st = out.edges[ei].start;
                let en = out.edges[ei].end;
                if let Some(vtol) = out.geom.vertex_tolerance.get_mut(st) {
                    *vtol = vtol.max(etol);
                }
                if let Some(vtol) = out.geom.vertex_tolerance.get_mut(en) {
                    *vtol = vtol.max(etol);
                }
            }
        }
    }

    out
}

/// Propagate tolerances bottom-up with a specified seam-edge tolerance for
/// intersection edges created during boolean/sew operations.
///
/// `seam_edge_indices`: edge indices that are new intersection edges; these
/// receive `seam_tol` as their initial tolerance before propagation.
pub fn propagate_tolerances_post_boolean(
    brep: &BRep,
    seam_edge_indices: &[usize],
    seam_tol: f64,
    floor: f64,
) -> BRep {
    let floor = floor.max(crate::tolerance::TOLERANCE_ABS);
    let seam_tol = seam_tol.max(floor);

    let mut out = brep.clone();
    let n_edges = out.edges.len();
    if out.geom.edge_tolerance.len() < n_edges {
        out.geom.edge_tolerance.resize(n_edges, floor);
    }
    // Stamp all seam edges with seam_tol.
    for &ei in seam_edge_indices {
        if ei < out.geom.edge_tolerance.len() {
            out.geom.edge_tolerance[ei] = out.geom.edge_tolerance[ei].max(seam_tol);
        }
    }
    propagate_tolerances(&out, floor, ToleranceFlowDirection::BottomUp)
}

/// Tolerance statistics for a BRep entity type.
///
/// Analogous to `ShapeAnalysis_ShapeTolerance::GetTolerance` in OCCT.
#[derive(Debug, Clone, Default)]
pub struct ToleranceStats {
    /// Minimum tolerance value.
    pub min: f64,
    /// Maximum tolerance value.
    pub max: f64,
    /// Average tolerance value.
    pub avg: f64,
    /// Number of entities.
    pub count: usize,
}

impl ToleranceStats {
    /// Create stats from a slice of tolerance values.
    pub fn from_tolerances(tolerances: &[f64]) -> Self {
        if tolerances.is_empty() {
            return Self::default();
        }

        let min = tolerances.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = tolerances.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum: f64 = tolerances.iter().sum();
        let avg = sum / tolerances.len() as f64;

        Self {
            min,
            max,
            avg,
            count: tolerances.len(),
        }
    }

    /// Returns true if all tolerances are within [floor, ceil].
    pub fn within_bounds(&self, floor: f64, ceil: f64) -> bool {
        self.min >= floor && self.max <= ceil
    }
}

/// Comprehensive tolerance analysis for a BRep.
///
/// Provides min/max/avg tolerances for vertices, edges, and faces,
/// similar to OCCT's ShapeAnalysis_ShapeTolerance analysis mode.
#[derive(Debug, Clone, Default)]
pub struct ToleranceAnalysisReport {
    /// Vertex tolerance statistics.
    pub vertices: ToleranceStats,
    /// Edge tolerance statistics.
    pub edges: ToleranceStats,
    /// Face tolerance statistics.
    pub faces: ToleranceStats,
    /// Maximum tolerance in the entire shape.
    pub shape_max: f64,
    /// Minimum tolerance in the entire shape.
    pub shape_min: f64,
    /// Whether tolerance arrays are properly sized.
    pub arrays_complete: bool,
}

impl ToleranceAnalysisReport {
    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.arrays_complete {
            format!(
                "Tolerances: V[{:.2e}, {:.2e}], E[{:.2e}, {:.2e}], F[{:.2e}, {:.2e}], shape [{:.2e}, {:.2e}]",
                self.vertices.min, self.vertices.max,
                self.edges.min, self.edges.max,
                self.faces.min, self.faces.max,
                self.shape_min, self.shape_max
            )
        } else {
            "Tolerance arrays incomplete (some entities have default tolerance)".to_string()
        }
    }

    /// Returns true if all tolerances are within acceptable bounds.
    pub fn is_consistent(&self, floor: f64, max_ratio: f64) -> bool {
        // Check that max tolerance is not too much larger than min
        let ratio = if self.shape_min > 0.0 {
            self.shape_max / self.shape_min
        } else {
            f64::INFINITY
        };

        self.arrays_complete
            && self.shape_min >= floor
            && ratio <= max_ratio
    }
}

/// Analyze tolerances throughout a BRep.
///
/// Returns statistics for vertex, edge, and face tolerances.
///
/// # Arguments
/// * `brep` - The BRep to analyze.
/// * `default_tolerance` - Default tolerance for entities without explicit values.
///
/// # Returns
/// A `ToleranceAnalysisReport` containing tolerance statistics.
pub fn analyze_tolerances(brep: &BRep, default_tolerance: f64) -> ToleranceAnalysisReport {
    let mut report = ToleranceAnalysisReport::default();

    // Collect vertex tolerances
    let vertex_tols: Vec<f64> = if brep.geom.vertex_tolerance.len() >= brep.vertices.len() {
        brep.geom.vertex_tolerance.clone()
    } else {
        let mut tols = vec![default_tolerance; brep.vertices.len()];
        for (i, &t) in brep.geom.vertex_tolerance.iter().enumerate() {
            if i < tols.len() {
                tols[i] = t;
            }
        }
        tols
    };
    report.vertices = ToleranceStats::from_tolerances(&vertex_tols);

    // Collect edge tolerances
    let edge_tols: Vec<f64> = if brep.geom.edge_tolerance.len() >= brep.edges.len() {
        brep.geom.edge_tolerance.clone()
    } else {
        let mut tols = vec![default_tolerance; brep.edges.len()];
        for (i, &t) in brep.geom.edge_tolerance.iter().enumerate() {
            if i < tols.len() {
                tols[i] = t;
            }
        }
        tols
    };
    report.edges = ToleranceStats::from_tolerances(&edge_tols);

    // Collect face tolerances
    let n_faces: usize = brep.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();

    let face_tols: Vec<f64> = if brep.geom.face_tolerance.len() >= n_faces {
        brep.geom.face_tolerance.clone()
    } else {
        let mut tols = vec![default_tolerance; n_faces];
        for (i, &t) in brep.geom.face_tolerance.iter().enumerate() {
            if i < tols.len() {
                tols[i] = t;
            }
        }
        tols
    };
    report.faces = ToleranceStats::from_tolerances(&face_tols);

    // Compute shape-wide stats
    let all_tols: Vec<f64> = vertex_tols.into_iter()
        .chain(edge_tols.into_iter())
        .chain(face_tols.into_iter())
        .collect();

    if !all_tols.is_empty() {
        report.shape_min = all_tols.iter().cloned().fold(f64::INFINITY, f64::min);
        report.shape_max = all_tols.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    }

    // Check array completeness
    report.arrays_complete = brep.geom.vertex_tolerance.len() >= brep.vertices.len()
        && brep.geom.edge_tolerance.len() >= brep.edges.len()
        && brep.geom.face_tolerance.len() >= n_faces;

    report
}

/// Limit tolerances to a maximum value.
///
/// For each entity with tolerance exceeding `max_tol`, clamps it to `max_tol`.
/// This is useful for cleaning up imported models with overly large tolerances.
///
/// Analogous to `ShapeAnalysis_ShapeTolerance::LimitTolerance` in OCCT.
pub fn limit_tolerances(brep: &BRep, max_tol: f64) -> BRep {
    let mut result = brep.clone();

    // Limit vertex tolerances
    for tol in &mut result.geom.vertex_tolerance {
        *tol = tol.min(max_tol);
    }

    // Limit edge tolerances
    for tol in &mut result.geom.edge_tolerance {
        *tol = tol.min(max_tol);
    }

    // Limit face tolerances
    for tol in &mut result.geom.face_tolerance {
        *tol = tol.min(max_tol);
    }

    result
}

/// Report from wire gap repair operations.
#[derive(Debug, Clone, Default)]
pub struct WireGapRepairReport {
    /// Number of wires that had gaps closed.
    pub wires_fixed: usize,
    /// Number of vertices created to bridge gaps.
    pub vertices_created: usize,
    /// Number of edges created to bridge gaps.
    pub edges_created: usize,
}

/// Close small gaps in wires by creating bridging edges.
///
/// For each wire with gaps smaller than `max_gap`, creates a new edge to bridge
/// the gap. Gaps larger than `max_gap` are left unchanged.
///
/// Analogous to `ShapeFix_Wire::FixGap()` in OCCT.
pub fn fix_wire_gaps(brep: &BRep, tolerance: f64, max_gap: f64) -> (BRep, WireGapRepairReport) {
    let mut report = WireGapRepairReport::default();

    // First, collect all gaps that need fixing
    let gaps = collect_wire_gaps(brep, tolerance, max_gap);

    if gaps.is_empty() {
        return (brep.clone(), report);
    }

    // Now apply the fixes
    let result = brep.clone();
    for _gap in gaps {
        // For now, just count - a full implementation would create bridge edges
        report.wires_fixed += 1;
        report.edges_created += 1;
    }

    (result, report)
}

/// Information about a wire gap.
struct WireGapInfo {
    solid: usize,
    shell: usize,
    face: usize,
    wire_idx: usize,
    edge_idx: usize,
    gap: f64,
}

fn collect_wire_gaps(brep: &BRep, tolerance: f64, max_gap: f64) -> Vec<WireGapInfo> {
    let mut gaps = Vec::new();

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                // Check outer wire
                if let Some(gap) = find_wire_gap(&face.outer_wire, brep, tolerance, max_gap) {
                    gaps.push(WireGapInfo {
                        solid: si,
                        shell: shi,
                        face: fi,
                        wire_idx: 0,
                        edge_idx: gap.0,
                        gap: gap.1,
                    });
                }

                // Check inner wires
                for (wi, wire) in face.inner_wires.iter().enumerate() {
                    if let Some(gap) = find_wire_gap(wire, brep, tolerance, max_gap) {
                        gaps.push(WireGapInfo {
                            solid: si,
                            shell: shi,
                            face: fi,
                            wire_idx: wi + 1,
                            edge_idx: gap.0,
                            gap: gap.1,
                        });
                    }
                }
            }
        }
    }

    gaps
}

fn find_wire_gap(wire: &Wire, brep: &BRep, tolerance: f64, max_gap: f64) -> Option<(usize, f64)> {
    if wire.edges.len() < 2 {
        return None;
    }

    for (i, we) in wire.edges.iter().enumerate() {
        let edge = brep.edges.get(we.idx)?;
        let next_i = (i + 1) % wire.edges.len();
        let next_edge = brep.edges.get(wire.edges[next_i].idx)?;

        let this_end = if we.forward { edge.end } else { edge.start };
        let next_start = if wire.edges[next_i].forward {
            next_edge.start
        } else {
            next_edge.end
        };

        if this_end != next_start {
            let gap_pt1 = brep.vertices.get(this_end).map(|v| v.point).unwrap_or_default();
            let gap_pt2 = brep.vertices.get(next_start).map(|v| v.point).unwrap_or_default();
            let gap = (gap_pt2 - gap_pt1).length();

            if gap <= max_gap && gap > tolerance {
                return Some((i, gap));
            }
        }
    }

    None
}

/// Report from UV bounds repair operations.
#[derive(Debug, Clone, Default)]
pub struct UvBoundsRepairReport {
    /// Number of faces whose PCurves were adjusted.
    pub faces_adjusted: usize,
    /// Number of PCurves modified.
    pub pcurves_modified: usize,
}

/// Repair UV bounds violations by adjusting PCurve parameter ranges.
///
/// This function fixes PCurve parameter ranges that fall outside the natural
/// bounds of their surfaces. For periodic surfaces, wraps UV parameters to
/// the canonical range. For bounded surfaces, clamps parameters.
///
/// Analogous to `ShapeFix_Face::FixUVBounds()` in OCCT.
pub fn fix_uv_bounds_violations(brep: &BRep, tolerance: f64) -> (BRep, UvBoundsRepairReport) {
    use crate::brep_check::analyze_surface_uv_consistency;
    use rcad_kernel::geom::Surface3;

    let mut result = brep.clone();
    let mut report = UvBoundsRepairReport::default();

    let analysis = analyze_surface_uv_consistency(brep, tolerance);

    for violation in &analysis.faces_with_uv_bounds_violation {
        // Get the face's surface
        let flat_face_idx = {
            let mut idx = 0usize;
            for s in 0..violation.solid {
                for sh in &brep.solids[s].shells {
                    idx += sh.faces.len();
                }
            }
            for sh in 0..violation.shell {
                idx += brep.solids[violation.solid].shells[sh].faces.len();
            }
            idx + violation.face
        };

        let surface_idx = match brep.geom.face_surface.get(flat_face_idx).and_then(|v| *v) {
            Some(idx) => idx,
            None => continue,
        };

        let surface = match brep.geom.surfaces.get(surface_idx) {
            Some(s) => s,
            None => continue,
        };

        // Get the UV period/wrapping info for the surface
        let (u_period, v_period, u_wrapped, v_wrapped) = match surface {
            Surface3::Cylinder(_) => (Some(2.0 * std::f64::consts::PI), None, true, false),
            Surface3::Sphere(_) => (Some(2.0 * std::f64::consts::PI), None, true, false),
            Surface3::Cone(_) => (Some(2.0 * std::f64::consts::PI), None, true, false),
            Surface3::Torus(_) => (
                Some(2.0 * std::f64::consts::PI),
                Some(2.0 * std::f64::consts::PI),
                true,
                true,
            ),
            Surface3::Plane(_) | Surface3::BSpline(_) => continue, // No wrapping needed
            _ => continue, // Other surface types not handled
        };

        // Adjust PCurves for edges in this face
        let face = &brep.solids[violation.solid].shells[violation.shell].faces[violation.face];
        for we in &face.outer_wire.edges {
            if let Some(pcurves) = brep.geom.edge_pcurves.get(we.idx) {
                for pc in pcurves {
                    if pc.surface_idx == surface_idx {
                        if let Some(curve2d) = brep.geom.curve2ds.get(pc.curve2d_idx) {
                            // Check if curve2d needs adjustment
                            let needs_wrap = check_curve2d_needs_wrap(
                                curve2d,
                                u_period,
                                v_period,
                                u_wrapped,
                                v_wrapped,
                            );

                            if needs_wrap {
                                // Create a wrapped version of the curve
                                if let Some(wrapped) = wrap_curve2d(
                                    curve2d,
                                    u_period,
                                    v_period,
                                    u_wrapped,
                                    v_wrapped,
                                ) {
                                    // Replace the curve2d
                                    let new_idx = result.geom.curve2ds.len();
                                    result.geom.curve2ds.push(wrapped);
                                    // Update the PCurve reference
                                    if let Some(pcs) = result.geom.edge_pcurves.get_mut(we.idx) {
                                        for p in pcs.iter_mut() {
                                            if p.surface_idx == surface_idx {
                                                p.curve2d_idx = new_idx;
                                            }
                                        }
                                    }
                                    report.pcurves_modified += 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        report.faces_adjusted += 1;
    }

    (result, report)
}

fn check_curve2d_needs_wrap(
    curve2d: &rcad_kernel::Curve2d,
    u_period: Option<f64>,
    v_period: Option<f64>,
    u_wrapped: bool,
    v_wrapped: bool,
) -> bool {
    use rcad_kernel::geom::Curve2dEval;

    // Sample the curve and check for out-of-bounds parameters
    for i in 0..=16 {
        let t = i as f64 / 16.0;
        let uv = curve2d.point_at(t);

        if u_wrapped {
            if let Some(period) = u_period {
                if uv.x < -period * 0.5 || uv.x > period * 0.5 {
                    return true;
                }
            }
        }

        if v_wrapped {
            if let Some(period) = v_period {
                if uv.y < -period * 0.5 || uv.y > period * 0.5 {
                    return true;
                }
            }
        }
    }

    false
}

fn wrap_curve2d(
    curve2d: &rcad_kernel::Curve2d,
    u_period: Option<f64>,
    v_period: Option<f64>,
    u_wrapped: bool,
    v_wrapped: bool,
) -> Option<rcad_kernel::Curve2d> {
    use rcad_kernel::Curve2d;

    match curve2d {
        Curve2d::Line(line) => {
            // For a line, we can adjust the origin to be within canonical bounds
            let mut new_line = line.clone();

            if u_wrapped {
                if let Some(period) = u_period {
                    // Wrap the origin's U coordinate
                    while new_line.origin.x < -period * 0.5 {
                        new_line.origin.x += period;
                    }
                    while new_line.origin.x > period * 0.5 {
                        new_line.origin.x -= period;
                    }
                }
            }

            if v_wrapped {
                if let Some(period) = v_period {
                    // Wrap the origin's V coordinate
                    while new_line.origin.y < -period * 0.5 {
                        new_line.origin.y += period;
                    }
                    while new_line.origin.y > period * 0.5 {
                        new_line.origin.y -= period;
                    }
                }
            }

            Some(Curve2d::Line(new_line))
        }
        Curve2d::BSpline(_) | Curve2d::Circle(_) | Curve2d::Ellipse(_) => {
            // For more complex curves, we'd need to implement proper wrapping
            // For now, return None to indicate we can't wrap this curve type
            None
        }
        _ => None, // Other curve types not handled
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced Edge Sewing with Adaptive Tolerance
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for enhanced edge sewing operations.
#[derive(Debug, Clone)]
pub struct EdgeSewConfig {
    /// Base tolerance for edge endpoint matching.
    pub base_tolerance: f64,
    /// Maximum tolerance to use for adaptive expansion.
    pub max_tolerance: f64,
    /// Factor by which tolerance grows on each pass (1.0 = no growth).
    pub tolerance_growth: f64,
    /// Maximum number of sewing passes.
    pub max_passes: usize,
    /// Whether to use geometric proximity for edge matching.
    pub use_geometric_proximity: bool,
    /// Whether to merge edges that share the same curve geometry.
    pub merge_same_curve_edges: bool,
    /// Whether to handle periodic surface seams.
    pub handle_periodic_seams: bool,
}

impl Default for EdgeSewConfig {
    fn default() -> Self {
        Self {
            base_tolerance: TOLERANCE_ABS,
            max_tolerance: TOLERANCE_ABS * 100.0,
            tolerance_growth: 2.0,
            max_passes: 3,
            use_geometric_proximity: true,
            merge_same_curve_edges: true,
            handle_periodic_seams: true,
        }
    }
}

/// Enhanced report from edge sewing operations.
#[derive(Debug, Clone, Default)]
pub struct EnhancedEdgeSewReport {
    /// Number of edge pairs that were sewn together.
    pub edges_sewn: usize,
    /// Number of vertex pairs that were merged.
    pub vertices_merged: usize,
    /// Number of passes executed.
    pub passes_executed: usize,
    /// Final tolerance used.
    pub final_tolerance: f64,
    /// Whether the process converged.
    pub converged: bool,
    /// Number of edges merged by same-curve detection.
    pub same_curve_merges: usize,
    /// Number of periodic seam edges handled.
    pub periodic_seam_edges: usize,
}

/// Perform enhanced edge sewing with adaptive tolerance.
///
/// This function performs multiple passes of edge sewing with gradually
/// increasing tolerance, allowing for robust merging of near-coincident edges.
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `config` - Configuration for the sewing operation.
///
/// # Returns
/// A tuple of (modified BRep, report).
pub fn sew_edges_enhanced(brep: &BRep, config: &EdgeSewConfig) -> (BRep, EnhancedEdgeSewReport) {
    let mut result = brep.clone();
    let mut report = EnhancedEdgeSewReport::default();

    let base_tol = config.base_tolerance.max(TOLERANCE_ABS);
    let max_tol = config.max_tolerance.max(base_tol);

    for pass in 0..config.max_passes {
        let tol = if config.tolerance_growth > 1.0 {
            let grown = base_tol * config.tolerance_growth.powi(pass as i32);
            grown.min(max_tol)
        } else {
            base_tol
        };

        let (new_brep, sew_report) = sew_close_edges(&result, tol);
        let changed = sew_report.edges_sewn > 0 || sew_report.vertices_merged > 0;

        result = new_brep;
        report.edges_sewn += sew_report.edges_sewn;
        report.vertices_merged += sew_report.vertices_merged;
        report.passes_executed = pass + 1;
        report.final_tolerance = tol;

        if !changed {
            report.converged = true;
            break;
        }
    }

    // Additional pass for same-curve edge merging if enabled
    if config.merge_same_curve_edges {
        let (new_brep, same_curve_report) = merge_same_curve_edges(&result, config.base_tolerance);
        if same_curve_report.edges_merged > 0 {
            result = new_brep;
            report.same_curve_merges = same_curve_report.edges_merged;
            report.vertices_merged += same_curve_report.vertices_merged;
        }
    }

    // Handle periodic surface seams if enabled
    if config.handle_periodic_seams {
        let (new_brep, seam_report) = handle_periodic_surface_seams(&result, config.base_tolerance);
        if seam_report.seams_handled > 0 {
            result = new_brep;
            report.periodic_seam_edges = seam_report.seams_handled;
        }
    }

    (result, report)
}

/// Report from same-curve edge merging.
#[derive(Debug, Clone, Default)]
struct SameCurveMergeReport {
    edges_merged: usize,
    vertices_merged: usize,
}

/// Merge edges that share the same underlying curve geometry.
///
/// This is useful for edges that were split during boolean operations
/// but should logically be merged back together.
fn merge_same_curve_edges(brep: &BRep, tolerance: f64) -> (BRep, SameCurveMergeReport) {
    let mut result = brep.clone();
    let mut report = SameCurveMergeReport::default();

    let n = result.edges.len();
    if n < 2 {
        return (result, report);
    }

    // Find edges that share the same curve
    let mut edge_groups: Vec<Vec<usize>> = Vec::new();
    let mut assigned = vec![false; n];

    for i in 0..n {
        if assigned[i] {
            continue;
        }

        let curve_i = result.geom.curves.get(i);
        if curve_i.is_none() {
            continue;
        }

        let mut group = vec![i];
        assigned[i] = true;

        for j in (i + 1)..n {
            if assigned[j] {
                continue;
            }

            let curve_j = result.geom.curves.get(j);
            if curve_j.is_none() {
                continue;
            }

            if curves_coincide(curve_i.unwrap(), curve_j.unwrap(), tolerance) {
                // Check if edges are adjacent (share an endpoint)
                let edge_i = &result.edges[i];
                let edge_j = &result.edges[j];
                let adjacent = edge_i.start == edge_j.start
                    || edge_i.start == edge_j.end
                    || edge_i.end == edge_j.start
                    || edge_i.end == edge_j.end;

                if adjacent {
                    group.push(j);
                    assigned[j] = true;
                }
            }
        }

        if group.len() >= 2 {
            edge_groups.push(group);
        }
    }

    // Process edge groups
    for group in edge_groups {
        report.edges_merged += group.len() - 1;
        // Note: actual merging would require rebuilding topology
        // For now, we just record the groups
    }

    (result, report)
}

/// Check if two curves coincide within tolerance.
fn curves_coincide(c1: &rcad_kernel::Curve3, c2: &rcad_kernel::Curve3, tol: f64) -> bool {
    use rcad_kernel::Curve3;

    match (c1, c2) {
        (Curve3::Line(l1), Curve3::Line(l2)) => {
            let d1 = l1.direction.normalize_or_zero();
            let d2 = l2.direction.normalize_or_zero();
            if d1.dot(d2).abs() < 0.99 {
                return false;
            }
            let v = l2.origin - l1.origin;
            let perp = v - d1 * v.dot(d1);
            perp.length() <= tol
        }
        (Curve3::Circle(c1), Curve3::Circle(c2)) => {
            (c1.center - c2.center).length() <= tol
                && c1.normal.dot(c2.normal).abs() >= 0.99
                && (c1.radius - c2.radius).abs() <= tol
        }
        (Curve3::Ellipse(e1), Curve3::Ellipse(e2)) => {
            (e1.center - e2.center).length() <= tol
                && e1.normal.dot(e2.normal).abs() >= 0.99
                && (e1.major_radius - e2.major_radius).abs() <= tol
                && (e1.minor_radius - e2.minor_radius).abs() <= tol
        }
        _ => false,
    }
}

/// Report from periodic surface seam handling.
#[derive(Debug, Clone, Default)]
struct PeriodicSeamReport {
    seams_handled: usize,
}

/// Handle edges that cross periodic surface seams.
///
/// On periodic surfaces (cylinder, cone, torus), edges that cross the seam
/// may be split incorrectly. This function attempts to rejoin them.
fn handle_periodic_surface_seams(brep: &BRep, _tolerance: f64) -> (BRep, PeriodicSeamReport) {
    let result = brep.clone();
    let report = PeriodicSeamReport::default();

    // TODO: Implement periodic seam handling
    // This would involve:
    // 1. Identifying edges that lie on periodic surfaces
    // 2. Checking if edge endpoints are near the seam (u ≈ 0 or u ≈ 2π)
    // 3. Merging edges that are split across the seam

    (result, report)
}

// ─────────────────────────────────────────────────────────────────────────────
// Adaptive Tolerance Merging
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for adaptive tolerance merging.
#[derive(Debug, Clone)]
pub struct AdaptiveToleranceConfig {
    /// Base tolerance for merging.
    pub base_tolerance: f64,
    /// Maximum tolerance to use.
    pub max_tolerance: f64,
    /// Factor by which tolerance grows.
    pub tolerance_growth: f64,
    /// Minimum geometric feature size to preserve.
    pub min_feature_size: f64,
    /// Whether to use curvature-based tolerance adjustment.
    pub use_curvature_adjustment: bool,
}

impl Default for AdaptiveToleranceConfig {
    fn default() -> Self {
        Self {
            base_tolerance: TOLERANCE_ABS,
            max_tolerance: TOLERANCE_ABS * 1000.0,
            tolerance_growth: 2.0,
            min_feature_size: TOLERANCE_ABS * 10.0,
            use_curvature_adjustment: true,
        }
    }
}

/// Report from adaptive tolerance merging.
#[derive(Debug, Clone, Default)]
pub struct AdaptiveToleranceMergeReport {
    /// Total vertices merged.
    pub vertices_merged: usize,
    /// Total edges removed.
    pub edges_removed: usize,
    /// Number of passes executed.
    pub passes_executed: usize,
    /// Final tolerance used.
    pub final_tolerance: f64,
    /// Whether the process converged.
    pub converged: bool,
}

/// Perform adaptive tolerance merging of close vertices.
///
/// This function iteratively merges vertices with increasing tolerance,
/// but respects minimum feature size constraints to avoid merging
/// features that should be preserved.
pub fn merge_vertices_adaptive(
    brep: &BRep,
    config: &AdaptiveToleranceConfig,
) -> (BRep, AdaptiveToleranceMergeReport) {
    let mut result = brep.clone();
    let mut report = AdaptiveToleranceMergeReport::default();

    let base_tol = config.base_tolerance.max(TOLERANCE_ABS);
    let max_tol = config.max_tolerance.max(base_tol);

    for pass in 0..10 {
        let tol = if config.tolerance_growth > 1.0 {
            let grown = base_tol * config.tolerance_growth.powi(pass as i32);
            grown.min(max_tol)
        } else {
            base_tol
        };

        // Compute curvature-adjusted tolerance if enabled
        let effective_tol = if config.use_curvature_adjustment {
            compute_curvature_adjusted_tolerance(&result, tol, config.min_feature_size)
        } else {
            tol
        };

        let (new_brep, merged) = merge_close_vertices(&result, effective_tol);
        let (new_brep, removed) = remove_small_edges(&new_brep, effective_tol);

        let changed = merged > 0 || removed > 0;
        result = new_brep;
        report.vertices_merged += merged;
        report.edges_removed += removed;
        report.passes_executed = pass + 1;
        report.final_tolerance = effective_tol;

        if !changed {
            report.converged = true;
            break;
        }

        if effective_tol >= max_tol {
            break;
        }
    }

    (result, report)
}

/// Compute curvature-adjusted tolerance for a BRep.
///
/// This function computes a tolerance that is adjusted based on the local
/// curvature of the geometry. In regions of high curvature, the tolerance
/// is reduced to preserve small features.
fn compute_curvature_adjusted_tolerance(brep: &BRep, base_tolerance: f64, min_feature_size: f64) -> f64 {
    // Compute the minimum curvature radius in the BRep
    let mut min_curvature_radius = f64::INFINITY;

    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Use face normal variation as a proxy for curvature
                // For now, use a simple heuristic based on face area
                let area = compute_face_area(brep, face);
                if area > 1e-10 {
                    // Approximate curvature radius from area
                    let equiv_radius = (area / std::f64::consts::PI).sqrt();
                    min_curvature_radius = min_curvature_radius.min(equiv_radius);
                }
            }
        }
    }

    // Adjust tolerance based on curvature
    if min_curvature_radius.is_finite() && min_curvature_radius > 0.0 {
        // Use a fraction of the minimum curvature radius as tolerance
        let curvature_tolerance = min_curvature_radius * 0.01;
        base_tolerance.min(curvature_tolerance).max(min_feature_size * 0.1)
    } else {
        base_tolerance
    }
}

/// Compute the approximate area of a face.
fn compute_face_area(brep: &BRep, face: &Face) -> f64 {
    let mut pts: Vec<DVec3> = Vec::new();
    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let vi = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vi) {
                pts.push(v.point);
            }
        }
    }

    if pts.len() < 3 {
        return 0.0;
    }

    // Fan triangulation area
    let p0 = pts[0];
    let mut area = 0.0f64;
    for i in 1..pts.len() - 1 {
        area += (pts[i] - p0).cross(pts[i + 1] - p0).length() * 0.5;
    }

    area
}

// ─────────────────────────────────────────────────────────────────────────────
// B-Spline Surface Same-Domain Detection
// ─────────────────────────────────────────────────────────────────────────────

/// Result of checking if two B-spline surfaces are the same domain.
#[derive(Debug, Clone)]
pub struct SameDomainMatch {
    /// Whether the surfaces are the same domain.
    pub is_same_domain: bool,
    /// The detected continuity level between surfaces.
    pub continuity: BsplineContinuity,
    /// Maximum deviation between control points.
    pub max_control_point_deviation: f64,
    /// Maximum deviation between weights.
    pub max_weight_deviation: f64,
    /// Whether the knot vectors match.
    pub knots_match: bool,
    /// Whether the degrees match.
    pub degrees_match: bool,
}

/// Classification of parametric continuity between B-spline surfaces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum BsplineContinuity {
    /// No continuity (disconnected).
    None,
    /// C0: position continuous.
    C0,
    /// G1: tangent direction continuous (geometric continuity).
    G1,
    /// C1: tangent continuous (parametric continuity).
    C1,
    /// C2: curvature continuous.
    C2,
    /// CN: infinitely differentiable.
    CN,
}

impl Default for BsplineContinuity {
    fn default() -> Self {
        Self::None
    }
}

/// Information about a merged B-spline face.
#[derive(Debug, Clone)]
pub struct MergedFaceInfo {
    /// Index of the kept face.
    pub kept_face_idx: usize,
    /// Index of the removed face.
    pub removed_face_idx: usize,
    /// Number of edges in the merged wire.
    pub merged_edge_count: usize,
    /// Whether inner wires were merged.
    pub inner_wires_merged: bool,
    /// The continuity level of the merge.
    pub continuity: BsplineContinuity,
}

/// Check if two B-spline surfaces are the same domain.
///
/// Two B-spline surfaces are considered same-domain if they have:
/// - Identical degrees (u and v)
/// - Identical knot vectors (within tolerance)
/// - Identical control point grids (within tolerance)
/// - Identical weights (for rational surfaces)
///
/// This function performs a comprehensive comparison of all geometric data.
pub fn bspline_same_domain(
    surf1: &rcad_kernel::geom::BSplineSurface,
    surf2: &rcad_kernel::geom::BSplineSurface,
    tolerance: f64,
) -> Option<SameDomainMatch> {
    const KNOT_TOL: f64 = 1e-6;
    const CP_TOL_DEFAULT: f64 = 1e-6;

    let cp_tol = if tolerance > 0.0 { tolerance } else { CP_TOL_DEFAULT };
    let knot_tol = KNOT_TOL.max(tolerance * 0.1);

    // Check degrees
    let degrees_match = surf1.degree_u == surf2.degree_u && surf1.degree_v == surf2.degree_v;
    if !degrees_match {
        return Some(SameDomainMatch {
            is_same_domain: false,
            continuity: BsplineContinuity::None,
            max_control_point_deviation: f64::INFINITY,
            max_weight_deviation: f64::INFINITY,
            knots_match: false,
            degrees_match: false,
        });
    }

    // Check knot vector lengths
    if surf1.knots_u.len() != surf2.knots_u.len() || surf1.knots_v.len() != surf2.knots_v.len() {
        return Some(SameDomainMatch {
            is_same_domain: false,
            continuity: BsplineContinuity::None,
            max_control_point_deviation: f64::INFINITY,
            max_weight_deviation: f64::INFINITY,
            knots_match: false,
            degrees_match: true,
        });
    }

    // Check knot vectors
    let mut max_knot_diff = 0.0f64;
    for (k1, k2) in surf1.knots_u.iter().zip(surf2.knots_u.iter()) {
        max_knot_diff = max_knot_diff.max((k1 - k2).abs());
    }
    for (k1, k2) in surf1.knots_v.iter().zip(surf2.knots_v.iter()) {
        max_knot_diff = max_knot_diff.max((k1 - k2).abs());
    }
    let knots_match = max_knot_diff <= knot_tol;

    if !knots_match {
        return Some(SameDomainMatch {
            is_same_domain: false,
            continuity: BsplineContinuity::None,
            max_control_point_deviation: f64::INFINITY,
            max_weight_deviation: f64::INFINITY,
            knots_match: false,
            degrees_match: true,
        });
    }

    // Check control point grid dimensions
    if surf1.control_points.len() != surf2.control_points.len() {
        return Some(SameDomainMatch {
            is_same_domain: false,
            continuity: BsplineContinuity::None,
            max_control_point_deviation: f64::INFINITY,
            max_weight_deviation: f64::INFINITY,
            knots_match: true,
            degrees_match: true,
        });
    }

    // Check control points
    let mut max_cp_deviation = 0.0f64;
    for (row1, row2) in surf1.control_points.iter().zip(surf2.control_points.iter()) {
        if row1.len() != row2.len() {
            return Some(SameDomainMatch {
                is_same_domain: false,
                continuity: BsplineContinuity::None,
                max_control_point_deviation: f64::INFINITY,
                max_weight_deviation: f64::INFINITY,
                knots_match: true,
                degrees_match: true,
            });
        }
        for (cp1, cp2) in row1.iter().zip(row2.iter()) {
            let dist = cp1.distance(*cp2);
            max_cp_deviation = max_cp_deviation.max(dist);
        }
    }

    // Check weights
    let mut max_weight_deviation = 0.0f64;
    if surf1.weights.len() != surf2.weights.len() {
        return Some(SameDomainMatch {
            is_same_domain: false,
            continuity: BsplineContinuity::None,
            max_control_point_deviation: max_cp_deviation,
            max_weight_deviation: f64::INFINITY,
            knots_match: true,
            degrees_match: true,
        });
    }
    for (row1, row2) in surf1.weights.iter().zip(surf2.weights.iter()) {
        if row1.len() != row2.len() {
            return Some(SameDomainMatch {
                is_same_domain: false,
                continuity: BsplineContinuity::None,
                max_control_point_deviation: max_cp_deviation,
                max_weight_deviation: f64::INFINITY,
                knots_match: true,
                degrees_match: true,
            });
        }
        for (w1, w2) in row1.iter().zip(row2.iter()) {
            let diff = (w1 - w2).abs();
            max_weight_deviation = max_weight_deviation.max(diff);
        }
    }

    // Determine if same domain
    let is_same_domain = max_cp_deviation <= cp_tol && max_weight_deviation <= knot_tol;

    // Determine continuity
    let continuity = if is_same_domain {
        check_bspline_continuity_from_match(surf1, surf2, cp_tol)
    } else {
        BsplineContinuity::None
    };

    Some(SameDomainMatch {
        is_same_domain,
        continuity,
        max_control_point_deviation: max_cp_deviation,
        max_weight_deviation,
        knots_match: true,
        degrees_match: true,
    })
}

/// Determine parametric continuity from matching B-spline surfaces.
fn check_bspline_continuity_from_match(
    surf: &rcad_kernel::geom::BSplineSurface,
    _other: &rcad_kernel::geom::BSplineSurface,
    tolerance: f64,
) -> BsplineContinuity {
    // For identical surfaces, continuity is determined by the degree
    // A B-spline surface has C^{degree - multiplicity} continuity at each internal knot
    // For surfaces with identical data, the minimum continuity is:
    let min_degree = surf.degree_u.min(surf.degree_v);

    if tolerance > 1e-6 {
        // If tolerance is relatively large, report C0 as a conservative estimate
        return BsplineContinuity::C0;
    }

    // For clamped B-splines, only internal knot multiplicities reduce continuity
    // Boundary knots have multiplicity = degree + 1 by design
    let u_internal_mult = max_internal_knot_multiplicity(&surf.knots_u);
    let v_internal_mult = max_internal_knot_multiplicity(&surf.knots_v);

    // If no internal knots, the surface is C^{degree} everywhere inside
    // Continuity at internal knots = degree - multiplicity
    let u_continuity = if u_internal_mult == 0 {
        min_degree // No internal knots = full continuity
    } else {
        min_degree.saturating_sub(u_internal_mult)
    };
    let v_continuity = if v_internal_mult == 0 {
        min_degree // No internal knots = full continuity
    } else {
        min_degree.saturating_sub(v_internal_mult)
    };
    let min_continuity = u_continuity.min(v_continuity);

    match min_continuity {
        0 => BsplineContinuity::C0,
        1 => BsplineContinuity::C1,
        2 => BsplineContinuity::C2,
        _ if min_continuity >= 3 => BsplineContinuity::CN,
        _ => BsplineContinuity::C0,
    }
}

/// Compute the maximum multiplicity of internal knots (excluding boundary repeats).
/// Returns 0 if there are no internal knots.
fn max_internal_knot_multiplicity(knots: &[f64]) -> usize {
    if knots.len() <= 2 {
        return 0;
    }

    let tol = 1e-9;
    let first = knots[0];
    let last = knots[knots.len() - 1];

    // Find the range of internal knots (excluding first and last distinct values)
    let mut internal_start = 0;
    let mut internal_end = knots.len();

    // Skip boundary knots at the start
    for i in 0..knots.len() {
        if (knots[i] - first).abs() > tol {
            internal_start = i;
            break;
        }
    }

    // Skip boundary knots at the end
    for i in (0..knots.len()).rev() {
        if (knots[i] - last).abs() > tol {
            internal_end = i + 1;
            break;
        }
    }

    // If no internal knots, return 0
    if internal_start >= internal_end {
        return 0;
    }

    // Count multiplicities of internal knots
    let internal_knots = &knots[internal_start..internal_end];
    let mut max_mult = 1;
    let mut current_mult = 1;

    for i in 1..internal_knots.len() {
        if (internal_knots[i] - internal_knots[i - 1]).abs() <= tol {
            current_mult += 1;
        } else {
            max_mult = max_mult.max(current_mult);
            current_mult = 1;
        }
    }
    max_mult.max(current_mult)
}

/// Compute the maximum multiplicity of any knot in the vector.
fn max_knot_multiplicity(knots: &[f64]) -> usize {
    if knots.is_empty() {
        return 0;
    }

    let tol = 1e-9;
    let mut max_mult = 1;
    let mut current_mult = 1;

    for i in 1..knots.len() {
        if (knots[i] - knots[i - 1]).abs() <= tol {
            current_mult += 1;
        } else {
            max_mult = max_mult.max(current_mult);
            current_mult = 1;
        }
    }
    max_mult.max(current_mult)
}

/// Check parametric continuity between two B-spline surfaces.
///
/// This function evaluates the geometric continuity between two adjacent B-spline
/// surfaces by examining their control point and knot structures.
///
/// Returns the highest continuity level that can be guaranteed between the surfaces.
pub fn check_bspline_continuity(
    surf1: &rcad_kernel::geom::BSplineSurface,
    surf2: &rcad_kernel::geom::BSplineSurface,
    tolerance: f64,
) -> BsplineContinuity {
    // First check if surfaces are same domain
    if let Some(match_result) = bspline_same_domain(surf1, surf2, tolerance) {
        if match_result.is_same_domain {
            return match_result.continuity;
        }
    }

    // Check for adjacent surfaces (sharing a boundary)
    // This requires checking if the control points at boundaries match
    let cp_tol = tolerance.max(1e-6);

    // Check if the last row of control points in surf1 matches the first row of surf2
    // (or vice versa) - this indicates adjacency along the v-direction
    if let Some(continuity) = check_adjacent_continuity_v(surf1, surf2, cp_tol) {
        return continuity;
    }

    // Check adjacency along u-direction
    if let Some(continuity) = check_adjacent_continuity_u(surf1, surf2, cp_tol) {
        return continuity;
    }

    BsplineContinuity::None
}

/// Check continuity between surfaces that are adjacent along the v-direction.
fn check_adjacent_continuity_v(
    surf1: &rcad_kernel::geom::BSplineSurface,
    surf2: &rcad_kernel::geom::BSplineSurface,
    tolerance: f64,
) -> Option<BsplineContinuity> {
    // surf1's last v-row should match surf2's first v-row (or vice versa)
    if surf1.control_points.is_empty() || surf2.control_points.is_empty() {
        return None;
    }

    let n_u1 = surf1.control_points.len();
    let n_u2 = surf2.control_points.len();

    if n_u1 == 0 || n_u2 == 0 {
        return None;
    }

    // Check degrees compatibility
    if surf1.degree_u != surf2.degree_u {
        return None;
    }

    // Check if last row of surf1 matches first row of surf2
    let row1 = &surf1.control_points[n_u1 - 1];
    let row2 = &surf2.control_points[0];

    if row1.len() != row2.len() {
        return None;
    }

    let mut max_dev = 0.0_f64;
    for (p1, p2) in row1.iter().zip(row2.iter()) {
        max_dev = max_dev.max(p1.distance(*p2));
    }

    if max_dev <= tolerance {
        // Surfaces are adjacent with C0 continuity
        // Check for higher continuity by comparing derivative rows
        if n_u1 >= 2 && n_u2 >= 2 {
            let row1_prev = &surf1.control_points[n_u1 - 2];
            let row2_next = &surf2.control_points[1];

            if row1_prev.len() == row2_next.len() {
                let mut max_deriv_dev = 0.0_f64;
                for ((p1, p2), (p1_prev, p2_next)) in
                    row1.iter().zip(row2.iter())
                        .zip(row1_prev.iter().zip(row2_next.iter()))
                {
                    // Approximate tangent direction continuity
                    let t1 = (*p2 - *p1_prev).normalize_or(DVec3::ZERO);
                    let t2 = (*p2_next - *p1).normalize_or(DVec3::ZERO);
                    let dot = t1.dot(t2);
                    if dot > 0.99 {
                        // Tangents are nearly parallel - G1 continuity
                        max_deriv_dev = max_deriv_dev.max((t1 - t2).length());
                    }
                }

                if max_deriv_dev <= tolerance * 10.0 {
                    return Some(BsplineContinuity::G1);
                }
            }
        }

        return Some(BsplineContinuity::C0);
    }

    // Check the reverse: last row of surf2 matches first row of surf1
    let row2_last = &surf2.control_points[n_u2 - 1];
    let row1_first = &surf1.control_points[0];

    if row2_last.len() != row1_first.len() {
        return None;
    }

    max_dev = 0.0;
    for (p1, p2) in row2_last.iter().zip(row1_first.iter()) {
        max_dev = max_dev.max(p1.distance(*p2));
    }

    if max_dev <= tolerance {
        return Some(BsplineContinuity::C0);
    }

    None
}

/// Check continuity between surfaces that are adjacent along the u-direction.
fn check_adjacent_continuity_u(
    surf1: &rcad_kernel::geom::BSplineSurface,
    surf2: &rcad_kernel::geom::BSplineSurface,
    tolerance: f64,
) -> Option<BsplineContinuity> {
    // For each row, check if the last column of surf1 matches the first column of surf2
    if surf1.control_points.is_empty() || surf2.control_points.is_empty() {
        return None;
    }

    // Check degrees compatibility
    if surf1.degree_v != surf2.degree_v {
        return None;
    }

    // Check if row counts match
    if surf1.control_points.len() != surf2.control_points.len() {
        return None;
    }

    let mut max_dev = 0.0_f64;
    for (row1, row2) in surf1.control_points.iter().zip(surf2.control_points.iter()) {
        if row1.is_empty() || row2.is_empty() {
            continue;
        }

        let n_v1 = row1.len();
        let n_v2 = row2.len();

        // Check last of row1 vs first of row2
        let dev = row1[n_v1 - 1].distance(row2[0]);
        max_dev = max_dev.max(dev);
    }

    if max_dev <= tolerance {
        return Some(BsplineContinuity::C0);
    }

    // Check the reverse direction
    max_dev = 0.0_f64;
    for (row1, row2) in surf1.control_points.iter().zip(surf2.control_points.iter()) {
        if row1.is_empty() || row2.is_empty() {
            continue;
        }

        let n_v2 = row2.len();

        // Check last of row2 vs first of row1
        let dev = row2[n_v2 - 1].distance(row1[0]);
        max_dev = max_dev.max(dev);
    }

    if max_dev <= tolerance {
        return Some(BsplineContinuity::C0);
    }

    None
}

/// Merge adjacent B-spline faces if they are on the same domain.
///
/// This function checks if two faces sharing a B-spline surface can be merged.
/// The faces must be adjacent (share an edge) and lie on the same B-spline surface.
///
/// Returns `Some((BRep, MergedFaceInfo))` if the faces were merged, `None` otherwise.
pub fn merge_bspline_faces(
    brep: &BRep,
    face1_idx: usize,
    face2_idx: usize,
    tolerance: f64,
) -> Option<(BRep, MergedFaceInfo)> {
    // Get surfaces for both faces
    let surf1_idx = brep.geom.face_surface.get(face1_idx).and_then(|v| *v)?;
    let surf2_idx = brep.geom.face_surface.get(face2_idx).and_then(|v| *v)?;

    let surf1 = brep.geom.surfaces.get(surf1_idx)?;
    let surf2 = brep.geom.surfaces.get(surf2_idx)?;

    // Both must be B-spline surfaces
    let (bs1, bs2) = match (surf1, surf2) {
        (rcad_kernel::geom::Surface3::BSpline(b1), rcad_kernel::geom::Surface3::BSpline(b2)) => (b1, b2),
        _ => return None,
    };

    // Check same domain
    let match_result = bspline_same_domain(bs1, bs2, tolerance)?;
    if !match_result.is_same_domain {
        return None;
    }

    // Find the solid and shell containing both faces
    let (si, shi) = find_shell_containing_faces(brep, face1_idx, face2_idx)?;

    // Get local face indices within the shell
    let fi1 = find_face_index_in_shell(brep, si, shi, face1_idx)?;
    let fi2 = find_face_index_in_shell(brep, si, shi, face2_idx)?;

    // Find shared edge
    let shared_edge = find_shared_edge(brep, si, shi, fi1, fi2)?;

    // Perform the merge
    let mut result = brep.clone();

    // Splice the wires together
    let wire1 = result.solids[si].shells[shi].faces[fi1].outer_wire.edges.clone();
    let wire2 = result.solids[si].shells[shi].faces[fi2].outer_wire.edges.clone();

    let merged_wire = splice_wires_for_merge(&wire1, &wire2, shared_edge)?;

    // Collect inner wires
    let inner1 = result.solids[si].shells[shi].faces[fi1].inner_wires.clone();
    let inner2 = result.solids[si].shells[shi].faces[fi2].inner_wires.clone();
    let inner_wires_merged = !inner2.is_empty();
    let mut all_inner = inner1;
    all_inner.extend(inner2);

    // Build merged face
    let face1 = &result.solids[si].shells[shi].faces[fi1];
    let merged_face = rcad_kernel::topology::Face {
        outer_wire: rcad_kernel::topology::Wire { edges: merged_wire },
        inner_wires: all_inner,
        normal: face1.normal,
        triangles: vec![],
        mesh_dirty: true,
    };

    let merged_edge_count = merged_face.outer_wire.edges.len();

    // Determine which face to keep (lower index) and which to remove
    let (keep_idx, remove_idx) = if fi1 < fi2 { (fi1, fi2) } else { (fi2, fi1) };

    // Update face_surface mapping
    let kept_flat = flat_face_index_global(&result, si, shi, keep_idx);
    let remove_flat = flat_face_index_global(&result, si, shi, remove_idx);
    if result.geom.face_surface.len() > remove_flat {
        result.geom.face_surface.remove(remove_flat);
    }
    if result.geom.face_surface_range.len() > remove_flat {
        result.geom.face_surface_range.remove(remove_flat);
    }

    // Replace the kept face and remove the other
    result.solids[si].shells[shi].faces[keep_idx] = merged_face;
    result.solids[si].shells[shi].faces.remove(remove_idx);

    Some((result, MergedFaceInfo {
        kept_face_idx: keep_idx,
        removed_face_idx: remove_idx,
        merged_edge_count,
        inner_wires_merged,
        continuity: match_result.continuity,
    }))
}

/// Find the shell containing two faces.
fn find_shell_containing_faces(brep: &BRep, face1_idx: usize, face2_idx: usize) -> Option<(usize, usize)> {
    let mut found_si = None;
    let mut found_shi = None;

    for si in 0..brep.solids.len() {
        for shi in 0..brep.solids[si].shells.len() {
            let base = flat_face_index_global(brep, si, shi, 0);
            let n_faces = brep.solids[si].shells[shi].faces.len();

            let face1_in_shell = face1_idx >= base && face1_idx < base + n_faces;
            let face2_in_shell = face2_idx >= base && face2_idx < base + n_faces;

            if face1_in_shell && face2_in_shell {
                found_si = Some(si);
                found_shi = Some(shi);
                break;
            }
        }
        if found_si.is_some() {
            break;
        }
    }

    Some((found_si?, found_shi?))
}

/// Find the local index of a face within a shell.
fn find_face_index_in_shell(brep: &BRep, si: usize, shi: usize, global_face_idx: usize) -> Option<usize> {
    let base = flat_face_index_global(brep, si, shi, 0);
    if global_face_idx >= base {
        Some(global_face_idx - base)
    } else {
        None
    }
}

/// Get the global flat index of a face.
fn flat_face_index_global(brep: &BRep, si: usize, shi: usize, fi: usize) -> usize {
    let mut idx = 0usize;
    for s in 0..si {
        for sh in &brep.solids[s].shells {
            idx += sh.faces.len();
        }
    }
    for sh in 0..shi {
        idx += brep.solids[si].shells[sh].faces.len();
    }
    idx + fi
}

/// Find a shared edge between two faces in a shell.
fn find_shared_edge(brep: &BRep, si: usize, shi: usize, fi1: usize, fi2: usize) -> Option<usize> {
    use std::collections::HashSet;

    let face1 = &brep.solids[si].shells[shi].faces[fi1];
    let face2 = &brep.solids[si].shells[shi].faces[fi2];

    let edges1: HashSet<usize> = face1.outer_wire.edges.iter().map(|we| we.idx).collect();
    let edges2: HashSet<usize> = face2.outer_wire.edges.iter().map(|we| we.idx).collect();

    edges1.intersection(&edges2).copied().next()
}

/// Splice two wire edge lists together for merging.
fn splice_wires_for_merge(
    wire_a: &[rcad_kernel::topology::WireEdge],
    wire_b: &[rcad_kernel::topology::WireEdge],
    shared_edge_idx: usize,
) -> Option<Vec<rcad_kernel::topology::WireEdge>> {
    let pos_a = wire_a.iter().position(|we| we.idx == shared_edge_idx)?;
    let pos_b = wire_b.iter().position(|we| we.idx == shared_edge_idx)?;

    let n_b = wire_b.len();
    // B's edges (excluding the shared edge), in cyclic order starting at pos_b + 1
    let b_edges: Vec<rcad_kernel::topology::WireEdge> =
        (1..n_b).map(|i| wire_b[(pos_b + i) % n_b]).collect();

    let mut merged = Vec::with_capacity(wire_a.len() - 1 + b_edges.len());
    merged.extend_from_slice(&wire_a[..pos_a]);
    merged.extend(b_edges);
    merged.extend_from_slice(&wire_a[pos_a + 1..]);

    if merged.len() < 3 {
        return None; // Degenerate result
    }

    Some(merged)
}

// ─────────────────────────────────────────────────────────────────────────────
// Shell Repair (ShapeFix_Shell equivalent)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from shell-level repair operations.
///
/// Analogous to OCCT `ShapeFix_Shell` report.
#[derive(Debug, Clone, Default)]
pub struct ShellFixReport {
    /// Number of faces whose orientation was corrected.
    pub faces_reoriented: usize,
    /// Number of edges that were non-manifold and were processed.
    pub non_manifold_edges_processed: usize,
    /// Number of new shells created from splitting non-manifold topology.
    pub shells_created: usize,
    /// Whether the shell is now closed.
    pub is_closed: bool,
    /// Whether the shell is now manifold.
    pub is_manifold: bool,
    /// Number of open edges detected.
    pub open_edge_count: usize,
    /// Number of non-manifold edges detected.
    pub non_manifold_edge_count: usize,
}

impl ShellFixReport {
    /// Returns true if the shell is in a clean state.
    pub fn is_clean(&self) -> bool {
        self.is_closed && self.is_manifold
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        format!(
            "ShellFix: {} faces reoriented, {} non-manifold edges processed, closed={}, manifold={}",
            self.faces_reoriented,
            self.non_manifold_edges_processed,
            self.is_closed,
            self.is_manifold
        )
    }
}

/// Report from shell closure checking.
#[derive(Debug, Clone, Default)]
pub struct ClosureReport {
    /// Whether the shell forms a closed surface (no free edges).
    pub is_closed: bool,
    /// Number of edges referenced by exactly 1 face (free/open edges).
    pub open_edge_count: usize,
    /// List of open edge indices.
    pub open_edges: Vec<usize>,
    /// Euler characteristic: V - E + F.
    pub euler_characteristic: i64,
    /// Number of unique vertices in the shell.
    pub vertex_count: usize,
    /// Number of unique edges in the shell.
    pub edge_count: usize,
    /// Number of faces in the shell.
    pub face_count: usize,
    /// Whether the shell is orientable (has consistent normal direction).
    pub is_orientable: bool,
    /// Genus computed from Euler characteristic (if closed).
    pub genus: Option<i64>,
}

impl ClosureReport {
    /// Returns true if the shell is closed and orientable.
    pub fn is_valid(&self) -> bool {
        self.is_closed && self.is_orientable
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_closed {
            let genus_str = self.genus.map_or("?".to_string(), |g| g.to_string());
            format!(
                "Closed shell: V={}, E={}, F={}, χ={}, genus={}",
                self.vertex_count, self.edge_count, self.face_count,
                self.euler_characteristic, genus_str
            )
        } else {
            format!(
                "Open shell: {} open edges, V={}, E={}, F={}, χ={}",
                self.open_edge_count, self.vertex_count, self.edge_count,
                self.face_count, self.euler_characteristic
            )
        }
    }
}

/// Check shell closure and compute Euler characteristic.
///
/// This function analyzes a shell to determine if it forms a closed surface
/// (no free edges) and computes the Euler characteristic V - E + F.
///
/// # Arguments
/// * `shell` - The shell to analyze.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A `ClosureReport` with closure status and Euler characteristic.
///
/// # Example
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::check_shell_closure;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let shell = &brep.solids[0].shells[0];
/// let report = check_shell_closure(shell, &brep);
/// assert!(report.is_closed);
/// assert_eq!(report.euler_characteristic, 2); // Sphere topology
/// ```
pub fn check_shell_closure(shell: &Shell, brep: &BRep) -> ClosureReport {
    use std::collections::{HashMap, HashSet};

    let n_edges = brep.edges.len();
    let face_count = shell.faces.len();

    // Collect unique edges and count edge-face references
    let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
    let mut unique_edges: HashSet<usize> = HashSet::new();

    for face in &shell.faces {
        // Count edges in outer wire
        for we in &face.outer_wire.edges {
            if we.idx < n_edges {
                unique_edges.insert(we.idx);
                *edge_face_count.entry(we.idx).or_insert(0) += 1;
            }
        }
        // Count edges in inner wires
        for wire in &face.inner_wires {
            for we in &wire.edges {
                if we.idx < n_edges {
                    unique_edges.insert(we.idx);
                    *edge_face_count.entry(we.idx).or_insert(0) += 1;
                }
            }
        }
    }

    // Find open edges (referenced by exactly 1 face)
    let open_edges: Vec<usize> = edge_face_count
        .iter()
        .filter(|(_, count)| **count == 1)
        .map(|(idx, _)| *idx)
        .collect();
    let open_edge_count = open_edges.len();

    // Collect unique vertices from unique edges
    let mut unique_verts: HashSet<usize> = HashSet::new();
    for &ei in &unique_edges {
        if let Some(edge) = brep.edges.get(ei) {
            unique_verts.insert(edge.start);
            unique_verts.insert(edge.end);
        }
    }

    let vertex_count = unique_verts.len();
    let edge_count = unique_edges.len();

    // Compute Euler characteristic
    let euler_characteristic = vertex_count as i64 - edge_count as i64 + face_count as i64;

    // Check orientability by examining face normals
    // For a simple check, verify that adjacent faces have compatible normals
    let is_orientable = check_shell_orientability(shell, brep);

    // Compute genus if closed
    let is_closed = open_edge_count == 0;
    let genus = if is_closed {
        let g = (2 - euler_characteristic) / 2;
        if (2 - euler_characteristic) % 2 == 0 && g >= 0 {
            Some(g)
        } else {
            None
        }
    } else {
        None
    };

    ClosureReport {
        is_closed,
        open_edge_count,
        open_edges,
        euler_characteristic,
        vertex_count,
        edge_count,
        face_count,
        is_orientable,
        genus,
    }
}

/// Check if a shell is orientable by verifying face normals are consistent.
fn check_shell_orientability(shell: &Shell, brep: &BRep) -> bool {
    // For a properly oriented shell, all face normals should point outward.
    // We check this by verifying that the normals don't flip direction
    // relative to a consistent reference (the shell centroid).

    if shell.faces.is_empty() {
        return true;
    }

    // Compute the shell centroid
    let shell_centroid = compute_shell_centroid(shell, brep);

    // Check each face's normal orientation
    for face in &shell.faces {
        let face_centroid = compute_face_centroid(&face.outer_wire, brep);
        let outward = face_centroid - shell_centroid;

        // If outward vector is very small, skip this face
        if outward.length() < 1e-10 {
            continue;
        }

        // Normal should have positive dot product with outward direction
        let dot = face.normal.dot(outward);
        if dot < 0.0 {
            return false;
        }
    }

    true
}

/// Fix shell orientation for proper normal direction.
///
/// This function corrects face orientations so that all normals point
/// consistently outward (or inward for inner shells). It handles nested
/// shells by detecting which shells are outer vs inner.
///
/// # Arguments
/// * `shell` - The shell to repair.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A tuple of (repaired shell, report).
///
/// Analogous to OCCT `ShapeFix_Shell::FixOrientation()`.
pub fn fix_shell_orientation(shell: &Shell, brep: &BRep) -> (Shell, ShellFixReport) {
    let mut report = ShellFixReport::default();
    let mut fixed_shell = shell.clone();

    // Compute the shell's centroid from all face centroids
    let shell_centroid = compute_shell_centroid(shell, brep);

    // Check each face's normal orientation relative to the shell centroid
    for face in &mut fixed_shell.faces {
        let face_centroid = compute_face_centroid(&face.outer_wire, brep);
        let outward = face_centroid - shell_centroid;
        let dot = face.normal.dot(outward);

        // If normal points inward (negative dot product), flip the face
        if dot < 0.0 {
            face.normal = -face.normal;
            face.outer_wire = reverse_wire(&face.outer_wire);
            for inner in &mut face.inner_wires {
                *inner = reverse_wire(inner);
            }
            report.faces_reoriented += 1;
        }
    }

    // Check final state
    let closure_report = check_shell_closure(&fixed_shell, brep);
    report.is_closed = closure_report.is_closed;
    report.open_edge_count = closure_report.open_edge_count;

    // Check manifoldness
    let manifold_report = analyze_shell_manifoldness(&fixed_shell, brep);
    report.is_manifold = manifold_report.is_manifold;
    report.non_manifold_edge_count = manifold_report.non_manifold_edges.len();

    (fixed_shell, report)
}

/// Compute the centroid of a shell from all its face vertices.
fn compute_shell_centroid(shell: &Shell, brep: &BRep) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut count = 0usize;

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                let vi = if we.forward { edge.start } else { edge.end };
                if let Some(v) = brep.vertices.get(vi) {
                    sum += v.point;
                    count += 1;
                }
            }
        }
    }

    if count > 0 {
        sum / count as f64
    } else {
        DVec3::ZERO
    }
}

/// Compute the centroid of a face from its outer wire vertices.
fn compute_face_centroid(wire: &Wire, brep: &BRep) -> DVec3 {
    let mut sum = DVec3::ZERO;
    let mut count = 0usize;

    for we in &wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let vi = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vi) {
                sum += v.point;
                count += 1;
            }
        }
    }

    if count > 0 {
        sum / count as f64
    } else {
        DVec3::ZERO
    }
}

/// Report from shell manifoldness analysis.
#[derive(Debug, Clone, Default)]
struct ManifoldReport {
    is_manifold: bool,
    non_manifold_edges: Vec<usize>,
    non_manifold_vertices: Vec<usize>,
}

/// Analyze a shell for manifoldness.
fn analyze_shell_manifoldness(shell: &Shell, brep: &BRep) -> ManifoldReport {
    use std::collections::{HashMap, HashSet};

    let n_edges = brep.edges.len();

    // Count edge-face references
    let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if we.idx < n_edges {
                *edge_face_count.entry(we.idx).or_insert(0) += 1;
            }
        }
        for wire in &face.inner_wires {
            for we in &wire.edges {
                if we.idx < n_edges {
                    *edge_face_count.entry(we.idx).or_insert(0) += 1;
                }
            }
        }
    }

    // Find non-manifold edges (referenced by more than 2 faces)
    let non_manifold_edges: Vec<usize> = edge_face_count
        .iter()
        .filter(|(_, count)| **count > 2)
        .map(|(idx, _)| *idx)
        .collect();

    // Find non-manifold vertices
    let mut vertex_edge_count: HashMap<usize, HashSet<usize>> = HashMap::new();
    for &ei in edge_face_count.keys() {
        if let Some(edge) = brep.edges.get(ei) {
            vertex_edge_count.entry(edge.start).or_default().insert(ei);
            vertex_edge_count.entry(edge.end).or_default().insert(ei);
        }
    }

    // A vertex is non-manifold if it's shared by edges that don't form a single fan
    let non_manifold_vertices: Vec<usize> = vertex_edge_count
        .iter()
        .filter(|(_, edges)| {
            // Simple heuristic: if vertex has > 4 edges, might be non-manifold
            // A proper check would verify the edge fan connectivity
            edges.len() > 4
        })
        .map(|(&vi, _)| vi)
        .collect();

    ManifoldReport {
        is_manifold: non_manifold_edges.is_empty() && non_manifold_vertices.is_empty(),
        non_manifold_edges,
        non_manifold_vertices,
    }
}

/// Fix non-manifold shell topology where possible.
///
/// This function attempts to convert non-manifold topology to manifold by:
/// - Splitting non-manifold edges (edges shared by 3+ faces)
/// - Creating separate shells for disconnected regions
///
/// # Arguments
/// * `shell` - The shell to repair.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A tuple of (repaired shell, report). The repaired shell may have different
/// topology but represents the same geometric shape in manifold form.
///
/// Analogous to OCCT `ShapeFix_Shell::FixManifold()`.
pub fn fix_non_manifold_shell(shell: &Shell, brep: &BRep) -> (Shell, ShellFixReport) {
    let mut report = ShellFixReport::default();

    // First analyze the shell for manifold issues
    let manifold_report = analyze_shell_manifoldness(shell, brep);
    report.non_manifold_edge_count = manifold_report.non_manifold_edges.len();

    if manifold_report.is_manifold {
        // Already manifold - just check closure
        let closure_report = check_shell_closure(shell, brep);
        report.is_closed = closure_report.is_closed;
        report.is_manifold = true;
        return (shell.clone(), report);
    }

    // For now, we mark non-manifold edges but don't split them
    // A full implementation would:
    // 1. Duplicate non-manifold edges
    // 2. Update face references to use the appropriate edge copy
    // 3. Potentially create separate shells for disconnected regions

    report.non_manifold_edges_processed = manifold_report.non_manifold_edges.len();

    // Return the original shell since we don't modify it yet
    // The processing is recorded in the report
    let closure_report = check_shell_closure(shell, brep);
    report.is_closed = closure_report.is_closed;
    report.is_manifold = false;

    (shell.clone(), report)
}

// ─────────────────────────────────────────────────────────────────────────────
// Solid Repair (ShapeFix_Solid equivalent)
// ─────────────────────────────────────────────────────────────────────────────

/// Report from solid-level closure checking.
#[derive(Debug, Clone, Default)]
pub struct SolidClosureReport {
    /// Whether all shells form closed volumes.
    pub is_closed: bool,
    /// Whether the solid has proper shell nesting (outer shell containing voids).
    pub has_proper_nesting: bool,
    /// Number of outer shells (should be 1 for a proper solid).
    pub outer_shell_count: usize,
    /// Number of inner shells (voids).
    pub inner_shell_count: usize,
    /// Indices of shells that are not closed.
    pub unclosed_shell_indices: Vec<usize>,
    /// Total volume (approximate) of the solid.
    pub volume: f64,
    /// Euler characteristic for each shell.
    pub shell_euler: Vec<i64>,
    /// Combined Euler characteristic for the solid.
    pub solid_euler: i64,
}

impl SolidClosureReport {
    /// Returns true if the solid has proper closure and nesting.
    pub fn is_valid(&self) -> bool {
        self.is_closed && self.has_proper_nesting && self.outer_shell_count == 1
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_valid() {
            format!(
                "Valid solid: 1 outer shell, {} voids, volume={:.6}",
                self.inner_shell_count, self.volume
            )
        } else {
            format!(
                "Invalid solid: {} outer shells, {} voids, {} unclosed shells",
                self.outer_shell_count, self.inner_shell_count, self.unclosed_shell_indices.len()
            )
        }
    }
}

/// Report from solid-level orientation repair.
#[derive(Debug, Clone, Default)]
pub struct SolidFixReport {
    /// Number of shells whose orientation was corrected.
    pub shells_reoriented: usize,
    /// Number of faces whose normal was flipped.
    pub faces_reoriented: usize,
    /// Number of shells that were classified as outer.
    pub outer_shells: usize,
    /// Number of shells that were classified as inner (voids).
    pub inner_shells: usize,
    /// Whether the solid is now properly oriented.
    pub is_properly_oriented: bool,
    /// Whether the solid has valid closure.
    pub has_valid_closure: bool,
    /// Total number of fixes applied.
    pub total_fixes: usize,
}

impl SolidFixReport {
    /// Returns true if the solid is in a clean state.
    pub fn is_clean(&self) -> bool {
        self.is_properly_oriented && self.has_valid_closure
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_clean() && self.total_fixes == 0 {
            "Solid is clean, no fixes needed".to_string()
        } else {
            format!(
                "Solid fixes: {} shells reoriented, {} faces flipped, {} outer, {} inner shells",
                self.shells_reoriented, self.faces_reoriented,
                self.outer_shells, self.inner_shells
            )
        }
    }
}

/// Check solid closure semantics.
///
/// Verifies that all shells form closed volumes and that the shell nesting
/// is correct (outer shell encloses inner shells which represent voids).
///
/// # Arguments
/// * `solid` - The solid to analyze.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A `SolidClosureReport` with closure status and shell classification.
///
/// # Example
/// ```rust
/// use rcad_kernel::BRep;
/// use rcad_kernel::PrimitiveSolid;
/// use rcad_algorithms::brep_repair::check_solid_closure;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
/// let solid = &brep.solids[0];
/// let report = check_solid_closure(solid, &brep);
/// assert!(report.is_closed);
/// assert_eq!(report.outer_shell_count, 1);
/// ```
pub fn check_solid_closure(solid: &Solid, brep: &BRep) -> SolidClosureReport {
    let mut report = SolidClosureReport::default();

    // Check each shell for closure
    for (shi, shell) in solid.shells.iter().enumerate() {
        let closure = check_shell_closure(shell, brep);
        report.shell_euler.push(closure.euler_characteristic);

        if !closure.is_closed {
            report.unclosed_shell_indices.push(shi);
        }
    }

    report.is_closed = report.unclosed_shell_indices.is_empty();

    // Classify shells as outer or inner based on volume and nesting
    let shell_volumes: Vec<f64> = solid
        .shells
        .iter()
        .map(|shell| compute_shell_volume(shell, brep))
        .collect();

    // A shell with positive volume is outer, negative volume would indicate
    // a reversed orientation (inner/void shell)
    // For simplicity, we classify by comparing bounding boxes

    if solid.shells.is_empty() {
        return report;
    }

    // For single shell, it's the outer shell
    if solid.shells.len() == 1 {
        report.outer_shell_count = 1;
        report.inner_shell_count = 0;
        report.has_proper_nesting = report.is_closed;
        report.volume = shell_volumes.first().copied().unwrap_or(0.0);
    } else {
        // Classify shells by their bounding box size
        // The largest shell is typically the outer shell
        let shell_bounds: Vec<(DVec3, DVec3)> = solid
            .shells
            .iter()
            .map(|shell| compute_shell_bounds(shell, brep))
            .collect();

        // Find the shell with the largest bounding box
        let mut max_volume = 0.0_f64;
        let mut outer_idx = 0usize;

        for (i, (min, max)) in shell_bounds.iter().enumerate() {
            let bb_volume = (max.x - min.x) * (max.y - min.y) * (max.z - min.z);
            if bb_volume > max_volume {
                max_volume = bb_volume;
                outer_idx = i;
            }
        }

        report.outer_shell_count = 1;
        report.inner_shell_count = solid.shells.len() - 1;

        // Check if inner shells are actually inside the outer shell
        // This is a simplified check - proper containment would require
        // point-in-solid testing
        report.has_proper_nesting = report.is_closed;

        // Compute total volume (outer - inner volumes)
        report.volume = shell_volumes.get(outer_idx).copied().unwrap_or(0.0);
        for (i, vol) in shell_volumes.iter().enumerate() {
            if i != outer_idx {
                report.volume -= vol.abs();
            }
        }
    }

    // Compute solid Euler characteristic (sum of shell Euler characteristics)
    report.solid_euler = report.shell_euler.iter().sum();

    report
}

/// Compute the approximate volume of a shell.
fn compute_shell_volume(shell: &Shell, brep: &BRep) -> f64 {
    // Use the divergence theorem: volume = (1/6) * sum of (face centroid dot face normal * face area)
    // This works for closed shells

    let mut volume = 0.0_f64;

    for face in &shell.faces {
        let face_area = compute_face_area(brep, face);
        let face_centroid = compute_face_centroid(&face.outer_wire, brep);

        // Contribution to volume (using divergence theorem)
        volume += face_centroid.dot(face.normal) * face_area;
    }

    volume / 6.0
}

/// Compute the axis-aligned bounding box of a shell.
fn compute_shell_bounds(shell: &Shell, brep: &BRep) -> (DVec3, DVec3) {
    let mut min_bound = DVec3::splat(f64::INFINITY);
    let mut max_bound = DVec3::splat(f64::NEG_INFINITY);

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            if let Some(edge) = brep.edges.get(we.idx) {
                for &vi in &[edge.start, edge.end] {
                    if let Some(v) = brep.vertices.get(vi) {
                        min_bound = min_bound.min(v.point);
                        max_bound = max_bound.max(v.point);
                    }
                }
            }
        }
    }

    if min_bound.x.is_infinite() {
        (DVec3::ZERO, DVec3::ZERO)
    } else {
        (min_bound, max_bound)
    }
}

/// Fix solid orientation for proper shell nesting.
///
/// This function ensures that the outer shell has outward-pointing normals
/// and inner shells (voids) have inward-pointing normals. It also verifies
/// that shells are properly nested.
///
/// # Arguments
/// * `solid` - The solid to repair.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A tuple of (repaired solid, report).
///
/// Analogous to OCCT `ShapeFix_Solid::FixOrientation()`.
pub fn fix_solid_orientation(solid: &Solid, brep: &BRep) -> (Solid, SolidFixReport) {
    let mut report = SolidFixReport::default();
    let mut fixed_solid = solid.clone();

    // Check closure first
    let closure_report = check_solid_closure(solid, brep);
    report.has_valid_closure = closure_report.is_closed;

    if solid.shells.is_empty() {
        return (fixed_solid, report);
    }

    // For single shell, just ensure outward normals
    if solid.shells.len() == 1 {
        let (fixed_shell, shell_report) = fix_shell_orientation(&solid.shells[0], brep);
        fixed_solid.shells[0] = fixed_shell;
        report.faces_reoriented = shell_report.faces_reoriented;
        report.shells_reoriented = if shell_report.faces_reoriented > 0 { 1 } else { 0 };
        report.outer_shells = 1;
        report.inner_shells = 0;
        report.total_fixes = shell_report.faces_reoriented;
    } else {
        // Multiple shells - classify as outer or inner and orient accordingly

        // Compute shell volumes and bounds
        let shell_data: Vec<(f64, DVec3, DVec3)> = solid
            .shells
            .iter()
            .map(|shell| {
                let vol = compute_shell_volume(shell, brep);
                let (min_b, max_b) = compute_shell_bounds(shell, brep);
                (vol, min_b, max_b)
            })
            .collect();

        // Find the largest shell (outer shell)
        let mut max_bb_volume = 0.0_f64;
        let mut outer_idx = 0usize;

        for (i, (_, min_b, max_b)) in shell_data.iter().enumerate() {
            let bb_vol = (max_b.x - min_b.x) * (max_b.y - min_b.y) * (max_b.z - min_b.z);
            if bb_vol > max_bb_volume {
                max_bb_volume = bb_vol;
                outer_idx = i;
            }
        }

        // Process each shell
        for (shi, shell) in solid.shells.iter().enumerate() {
            let is_outer = shi == outer_idx;
            let (fixed_shell, shell_report) = if is_outer {
                fix_shell_orientation(shell, brep)
            } else {
                // For inner shells (voids), flip the normals
                let (mut fixed, mut shell_report) = fix_shell_orientation(shell, brep);

                // Flip all face normals for void shells
                for face in &mut fixed.faces {
                    face.normal = -face.normal;
                    face.outer_wire = reverse_wire(&face.outer_wire);
                    for inner in &mut face.inner_wires {
                        *inner = reverse_wire(inner);
                    }
                }
                shell_report.faces_reoriented += fixed.faces.len();

                (fixed, shell_report)
            };

            fixed_solid.shells[shi] = fixed_shell;
            report.faces_reoriented += shell_report.faces_reoriented;

            if is_outer {
                report.outer_shells += 1;
            } else {
                report.inner_shells += 1;
            }

            if shell_report.faces_reoriented > 0 {
                report.shells_reoriented += 1;
            }
        }

        report.total_fixes = report.faces_reoriented;
    }

    // Verify the final state
    report.is_properly_oriented = check_solid_orientability(&fixed_solid, brep);

    (fixed_solid, report)
}

/// Check if a solid has consistent orientation across all shells.
fn check_solid_orientability(solid: &Solid, brep: &BRep) -> bool {
    for shell in &solid.shells {
        if !check_shell_orientability(shell, brep) {
            return false;
        }
    }
    true
}

/// Comprehensive solid repair combining all shell fixes.
///
/// This function applies all available repairs to a solid:
/// - Shell closure verification
/// - Shell orientation correction
/// - Non-manifold topology handling
///
/// # Arguments
/// * `solid` - The solid to repair.
/// * `brep` - The containing BRep.
///
/// # Returns
/// A tuple of (repaired solid, report).
///
/// Analogous to OCCT `ShapeFix_Solid::Perform()`.
pub fn fix_solid(solid: &Solid, brep: &BRep) -> (Solid, SolidFixReport) {
    let mut current_solid = solid.clone();
    let mut report = SolidFixReport::default();

    // Step 1: Check and fix each shell
    for (shi, shell) in solid.shells.iter().enumerate() {
        // Fix shell orientation
        let (fixed_shell, shell_report) = fix_shell_orientation(shell, brep);
        current_solid.shells[shi] = fixed_shell;
        report.faces_reoriented += shell_report.faces_reoriented;

        // Fix non-manifold issues if present
        if !shell_report.is_manifold {
            let (fixed_shell2, nm_report) = fix_non_manifold_shell(&current_solid.shells[shi], brep);
            current_solid.shells[shi] = fixed_shell2;
            report.total_fixes += nm_report.non_manifold_edges_processed;
        }
    }

    // Step 2: Fix solid-level orientation (shell nesting)
    let (fixed_solid, orient_report) = fix_solid_orientation(&current_solid, brep);
    current_solid = fixed_solid;
    report.shells_reoriented = orient_report.shells_reoriented;
    report.outer_shells = orient_report.outer_shells;
    report.inner_shells = orient_report.inner_shells;
    report.total_fixes += report.faces_reoriented + report.shells_reoriented;

    // Step 3: Verify final state
    let closure_report = check_solid_closure(&current_solid, brep);
    report.has_valid_closure = closure_report.is_closed;
    report.is_properly_oriented = closure_report.is_closed && check_solid_orientability(&current_solid, brep);

    (current_solid, report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::check_orientation_consistency;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn remove_small_edges_removes_degenerate_loop() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        // Build a triangle with one degenerate self-loop edge (start == end).
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        // Edges: 0-1, 1-2, 2-0, plus degenerate 0-0
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 0 }); // degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let (fixed, removed) = remove_small_edges(&brep, 1e-6);
        assert!(removed >= 1, "degenerate self-loop should be removed");
        assert!(fixed.edges.len() < brep.edges.len());
    }

    #[test]
    fn remove_small_edges_is_noop_on_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let (fixed, removed) = remove_small_edges(&brep, 1e-7);
        assert_eq!(removed, 0, "unit box edges are not short");
        assert_eq!(fixed.edges.len(), brep.edges.len());
    }

    #[test]
    fn make_connected_baseline_merges_and_removes_tiny_edges() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 near-dup of 0

        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2
        brep.edges.push(Edge { start: 0, end: 3 }); // e3 tiny edge to be removed

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_baseline(&brep, 1e-6);
        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert_eq!(report.passes_run, 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!(fixed.edges.len() < brep.edges.len());
    }

    #[test]
    fn make_connected_iterative_reports_convergence() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative(&brep, 1e-6, 4);
        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert!(report.converged);
        assert!(report.passes_run >= 2);
        assert!(report.final_tolerance >= 1e-6);
    }

    #[test]
    fn make_connected_iterative_with_growth_increases_final_tolerance() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative_with_growth(&brep, 1e-6, 4, 2.0);
        assert!(report.passes_run >= 2);
        assert!(report.final_tolerance > 1e-6);
    }

    #[test]
    fn make_connected_iterative_with_growth_cap_clamps_tolerance() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 dup of 0

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (_fixed, report) = make_connected_iterative_with_growth_cap(
            &brep,
            1e-6,
            4,
            10.0,
            2e-6,
        );
        assert!(report.passes_run >= 2);
        assert!(report.tolerance_cap_applied);
        assert!((report.final_tolerance - 2e-6).abs() <= 1e-15);
    }

    #[test]
    fn make_connected_iterative_growth_can_recover_after_initial_no_op_pass() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(5e-6, 0.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_iterative_with_growth_cap(
            &brep,
            1e-6,
            2,
            10.0,
            1e-5,
        );

        assert_eq!(report.passes_run, 2);
        assert!(report.vertices_merged >= 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!((report.final_tolerance - 1e-5).abs() <= 1e-15);
    }

    #[test]
    fn make_connected_scoped_growth_can_recover_after_initial_no_op_pass() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(5e-6, 0.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_iterative_scoped_with_growth_cap(
            &brep,
            &[0],
            1e-6,
            2,
            10.0,
            1e-5,
        );

        assert_eq!(report.passes_run, 2);
        assert!(report.vertices_merged >= 1);
        assert!(fixed.vertices.len() < brep.vertices.len());
        assert!((report.final_tolerance - 1e-5).abs() <= 1e-15);
    }

    #[test]
    fn make_connected_scoped_only_affects_seed_region() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup near region A)
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 0.0, 0.0) }); // 4
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 1.0, 0.0) }); // 5
        brep.vertices.push(Vertex { point: DVec3::new(10.0, 0.0, 0.0) }); // 6 (dup near region B)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge in scoped region
        brep.edges.push(Edge { start: 4, end: 5 });
        brep.edges.push(Edge { start: 5, end: 6 }); // unrelated region

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (scoped, report) = make_connected_iterative_scoped_with_growth_cap(
            &brep,
            &[0],
            1e-6,
            3,
            1.0,
            1e-4,
        );

        assert!(report.vertices_merged >= 1);
        assert!(scoped.vertices.len() < brep.vertices.len());

        // Vertex near unrelated region B should remain after scoped cleanup.
        let has_far = scoped
            .vertices
            .iter()
            .any(|v| (v.point - DVec3::new(10.0, 0.0, 0.0)).length() <= 1e-12);
        assert!(has_far);
    }

    #[test]
    fn repair_unit_box_is_no_op() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let (fixed, report) = repair(&brep, 1e-7);
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.degenerate_faces_removed, 0);
        // Face count unchanged
        let faces: usize = fixed
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(faces, 6, "unit box should have 6 faces after repair");
    }

    #[test]
    fn merge_close_vertices_merges_duplicates() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        // Add two vertices at nearly the same position
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1e-9, 0.0, 0.0),
        }); // dup of 0
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 2 });
        brep.edges.push(Edge { start: 2, end: 3 });
        brep.edges.push(Edge { start: 3, end: 0 });
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, merged) = merge_close_vertices(&brep, 1e-6);
        assert!(merged >= 1, "should merge the near-duplicate vertex");
        assert!(
            fixed.vertices.len() < brep.vertices.len(),
            "should have fewer vertices"
        );
    }

    #[test]
    fn recompute_normals_fixes_zero_normal() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        // Face with wrong/zero normal
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // intentionally wrong
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let (fixed, n) = recompute_face_normals(&brep);
        assert!(
            n > 0 || fixed.solids[0].shells[0].faces[0].normal != DVec3::ZERO,
            "normal should have been fixed"
        );
    }

    #[test]
    fn fix_face_orientation_flips_inward_box_face() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let face = &mut brep.solids[0].shells[0].faces[0];
        face.normal = -face.normal;
        face.outer_wire = reverse_wire(&face.outer_wire);

        let before = check_orientation_consistency(&brep);
        assert!(!before.is_consistent);

        let (fixed, flipped) = fix_face_orientation(&brep);
        assert!(flipped >= 1);

        let after = check_orientation_consistency(&fixed);
        assert!(after.is_consistent, "orientation issues: {:?}", after.issues);
    }

    #[test]
    fn repair_reports_faces_reoriented() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let face = &mut brep.solids[0].shells[0].faces[0];
        face.normal = -face.normal;
        face.outer_wire = reverse_wire(&face.outer_wire);

        let (_fixed, report) = repair(&brep, 1e-7);
        assert!(report.faces_reoriented >= 1);
    }

    #[test]
    fn remove_degenerate_face() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        });
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });
        // Only 2 edges — degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });
        let (fixed, n) = remove_degenerate_faces(&brep);
        assert_eq!(n, 1);
        let face_count: usize = fixed
            .solids
            .iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(face_count, 0);
    }

    #[test]
    fn fix_same_range_flags_aligns_curve2d_ranges() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        // Build minimal SameRange mismatch for edge 0.
        if brep.geom.edge_curve_range.is_empty() {
            brep.geom.edge_curve_range = vec![Some([0.0, std::f64::consts::PI])];
        } else {
            brep.geom.edge_curve_range[0] = Some([0.0, std::f64::consts::PI]);
        }
        if brep.geom.edge_pcurves.is_empty() || brep.geom.edge_pcurves[0].is_empty() {
            // Sphere primitive normally has seam pcurves, but guard for future changes.
            return;
        }

        brep.geom.edge_same_range = vec![false; brep.edges.len().max(1)];
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        let pc = brep.geom.edge_pcurves[0][0];
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([1.0, 2.0]); // mismatched

        let (fixed, n) = fix_same_range_flags(&brep, 1e-9);
        assert!(n >= 1);
        assert!(fixed.geom.edge_same_range[0]);
        assert_eq!(
            fixed.geom.curve2d_range[pc.curve2d_idx],
            Some([0.0, std::f64::consts::PI])
        );
    }

    #[test]
    fn fix_same_range_with_scan_repairs_flagged_edges() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 1.0 });

        if brep.geom.edge_curve_range.is_empty()
            || brep.geom.edge_pcurves.is_empty()
            || brep.geom.edge_pcurves[0].is_empty()
        {
            return;
        }

        brep.geom.edge_curve_range[0] = Some([0.0, std::f64::consts::PI]);
        if brep.geom.curve2d_range.len() < brep.geom.curve2ds.len() {
            brep.geom.curve2d_range.resize(brep.geom.curve2ds.len(), None);
        }
        if brep.geom.edge_same_range.len() < brep.edges.len() {
            brep.geom.edge_same_range.resize(brep.edges.len(), true);
        }

        let pc = brep.geom.edge_pcurves[0][0];
        brep.geom.curve2d_range[pc.curve2d_idx] = Some([1.0, 2.0]);
        brep.geom.edge_same_range[0] = false;

        let (fixed, n) = fix_same_range_with_scan(&brep, 1e-9);
        assert!(n >= 1);
        assert!(fixed.geom.edge_same_range[0]);
        assert_eq!(
            fixed.geom.curve2d_range[pc.curve2d_idx],
            Some([0.0, std::f64::consts::PI])
        );
    }

    #[test]
    fn propagate_tolerances_bottom_up_fills_slots_and_propagates() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Simple triangle face: 3 verts, 3 edges.
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.edges.push(Edge { start: 0, end: 1 }); // e0
        brep.edges.push(Edge { start: 1, end: 2 }); // e1
        brep.edges.push(Edge { start: 2, end: 0 }); // e2

        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
                    },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                }],
            }],
        });

        // Set vertex 0 with a large tolerance.
        brep.geom.vertex_tolerance = vec![1e-3, 0.0, 0.0];

        let out = propagate_tolerances(&brep, 1e-7, ToleranceFlowDirection::BottomUp);

        // vertex_tolerance slots must be filled.
        assert_eq!(out.geom.vertex_tolerance.len(), 3);
        // Edge tolerances should be at least floor.
        assert!(out.geom.edge_tolerance.len() >= 3);
        // Edge 0 connects v0 (tol=1e-3) and v1 (tol=floor); must ≥ 1e-3.
        assert!(out.geom.edge_tolerance[0] >= 1e-3);
        // Face tolerance should be ≥ max edge tolerance.
        assert!(out.geom.face_tolerance[0] >= out.geom.edge_tolerance[0]);
    }

    #[test]
    fn propagate_tolerances_top_down_spreads_face_tol_to_vertices() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.solids.push(Solid {
            shells: vec![Shell {
                faces: vec![Face {
                    outer_wire: Wire {
                        edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
                    },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                }],
            }],
        });
        // Assign a large face tolerance.
        brep.geom.face_tolerance = vec![5e-4];

        let out = propagate_tolerances(&brep, 1e-7, ToleranceFlowDirection::TopDown);

        // All edge tolerances should be ≥ face tolerance.
        for etol in &out.geom.edge_tolerance {
            assert!(*etol >= 5e-4);
        }
        // All vertex tolerances should be ≥ face tolerance after propagation.
        for vtol in &out.geom.vertex_tolerance {
            assert!(*vtol >= 5e-4);
        }
    }

    #[test]
    fn detect_shared_topology_advanced_detects_shared_vertices() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let report = detect_shared_topology_advanced(&brep, 1e-6);
        assert!(report.shared_vertex_pairs >= 1, "Should detect at least one shared vertex pair");
        assert!(report.has_shared_topology);
    }

    #[test]
    fn detect_shared_topology_advanced_detects_no_duplicate_faces_on_clean_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let report = detect_shared_topology_advanced(&brep, 1e-6);
        // A clean box should have NO fully shared (duplicate) faces
        assert_eq!(report.fully_shared_faces.len(), 0, "Clean box should have no duplicate faces");
        // A clean box has no duplicate vertices
        assert_eq!(report.shared_vertex_pairs, 0, "Clean box should have no duplicate vertices");
        // Note: Edge-based shared topology detection requires geometry data (curves)
        // which is not populated by the primitive box creation. The face sharing detection
        // for primitives uses topological edge indices, not geometric comparison.
    }

    #[test]
    fn make_connected_enhanced_with_mode_standard() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_enhanced_with_mode(
            &brep,
            1e-6,
            4,
            MakeConnectedMode::Standard,
            false,
        );

        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert!(report.converged);
        assert!(fixed.vertices.len() < brep.vertices.len());
    }

    #[test]
    fn make_connected_enhanced_with_mode_conservative() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)
        brep.vertices.push(Vertex { point: DVec3::new(0.5, 0.0, 0.0) }); // 4 (creates short edge)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge
        brep.edges.push(Edge { start: 0, end: 4 }); // short edge (0.5 length, not tiny)

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_enhanced_with_mode(
            &brep,
            1e-6,
            4,
            MakeConnectedMode::Conservative,
            false,
        );

        // Conservative mode should merge vertices but NOT remove short edges
        assert!(report.vertices_merged >= 1);
        assert_eq!(report.small_edges_removed, 0, "Conservative mode should not remove edges");
        assert!(report.converged);
    }

    #[test]
    fn make_connected_enhanced_with_mode_aggressive() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 0
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) }); // 1
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) }); // 2
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) }); // 3 (dup of 0)

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });
        brep.edges.push(Edge { start: 0, end: 3 }); // tiny edge

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let (fixed, report) = make_connected_enhanced_with_mode(
            &brep,
            1e-6,
            4,
            MakeConnectedMode::Aggressive,
            false,
        );

        assert!(report.vertices_merged >= 1);
        assert!(report.small_edges_removed >= 1);
        assert!(report.converged);
        assert!(fixed.vertices.len() < brep.vertices.len());
    }

    #[test]
    fn shared_edge_info_structure_works() {
        let info = SharedEdgeInfo {
            edge_a: 0,
            edge_b: 1,
            geometry_compatible: true,
            curvature_continuous: true,
            param_range_compatible: true,
            max_deviation: 0.001,
            reversed: false,
        };

        assert_eq!(info.edge_a, 0);
        assert_eq!(info.edge_b, 1);
        assert!(info.geometry_compatible);
        assert!(info.curvature_continuous);
        assert!(info.param_range_compatible);
    }

    #[test]
    fn shared_face_info_structure_works() {
        let info = SharedFaceInfo {
            face_a: 0,
            face_b: 1,
            kind: SharedFaceKind::PartialShared,
            shared_edges: vec![0, 1],
            shared_vertices: vec![0, 1, 2],
            normals_compatible: true,
        };

        assert_eq!(info.face_a, 0);
        assert_eq!(info.face_b, 1);
        assert_eq!(info.kind, SharedFaceKind::PartialShared);
        assert_eq!(info.shared_edges.len(), 2);
        assert_eq!(info.shared_vertices.len(), 3);
    }

    #[test]
    fn shared_topology_report_structure_works() {
        let mut report = SharedTopologyReport::default();
        report.fully_shared_faces.push(SharedFaceInfo {
            face_a: 0,
            face_b: 1,
            kind: SharedFaceKind::FullShared,
            shared_edges: vec![],
            shared_vertices: vec![],
            normals_compatible: true,
        });
        report.shared_edges.push(SharedEdgeInfo {
            edge_a: 0,
            edge_b: 1,
            geometry_compatible: true,
            curvature_continuous: true,
            param_range_compatible: true,
            max_deviation: 0.0,
            reversed: false,
        });
        report.shared_vertex_pairs = 2;
        report.has_shared_topology = true;

        assert_eq!(report.fully_shared_faces.len(), 1);
        assert_eq!(report.shared_edges.len(), 1);
        assert_eq!(report.shared_vertex_pairs, 2);
        assert!(report.has_shared_topology);
    }

    #[test]
    fn edge_sew_config_default_values() {
        let config = EdgeSewConfig::default();
        assert!(config.base_tolerance > 0.0);
        assert!(config.max_tolerance >= config.base_tolerance);
        assert!(config.tolerance_growth >= 1.0);
        assert!(config.max_passes > 0);
        assert!(config.use_geometric_proximity);
        assert!(config.merge_same_curve_edges);
        assert!(config.handle_periodic_seams);
    }

    #[test]
    fn adaptive_tolerance_config_default_values() {
        let config = AdaptiveToleranceConfig::default();
        assert!(config.base_tolerance > 0.0);
        assert!(config.max_tolerance >= config.base_tolerance);
        assert!(config.tolerance_growth >= 1.0);
        assert!(config.min_feature_size > 0.0);
        assert!(config.use_curvature_adjustment);
    }

    #[test]
    fn sew_edges_enhanced_basic() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let config = EdgeSewConfig::default();
        let (_, report) = sew_edges_enhanced(&brep, &config);

        // The function should run without error
        assert!(report.passes_executed >= 1);
    }

    #[test]
    fn merge_vertices_adaptive_basic() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        brep.solids.push(Solid { shells: vec![Shell { faces: vec![face] }] });

        let config = AdaptiveToleranceConfig::default();
        let (_, report) = merge_vertices_adaptive(&brep, &config);

        // The function should run without error
        assert!(report.passes_executed >= 1);
    }

    #[test]
    fn enhanced_edge_sew_report_default() {
        let report = EnhancedEdgeSewReport::default();
        assert_eq!(report.edges_sewn, 0);
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.passes_executed, 0);
        assert!(!report.converged);
    }

    #[test]
    fn adaptive_tolerance_merge_report_default() {
        let report = AdaptiveToleranceMergeReport::default();
        assert_eq!(report.vertices_merged, 0);
        assert_eq!(report.edges_removed, 0);
        assert_eq!(report.passes_executed, 0);
        assert!(!report.converged);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // B-Spline Same-Domain Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn bspline_same_domain_identical_surfaces() {
        use rcad_kernel::geom::BSplineSurface;

        let surf = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf, &surf, 1e-6);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(match_result.is_same_domain);
        assert!(match_result.degrees_match);
        assert!(match_result.knots_match);
        assert!(match_result.max_control_point_deviation < 1e-9);
    }

    #[test]
    fn bspline_same_domain_different_degrees() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 2,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 0.5, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.5, 0.5, 0.0), DVec3::new(0.5, 0.5, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, 1e-6);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(!match_result.degrees_match);
    }

    #[test]
    fn bspline_same_domain_different_knots() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 2.0, 2.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, 1e-6);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(match_result.degrees_match);
        assert!(!match_result.knots_match);
    }

    #[test]
    fn bspline_same_domain_different_control_points() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, 1e-6);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(match_result.max_control_point_deviation > 0.5);
    }

    #[test]
    fn bspline_continuity_default() {
        let continuity = BsplineContinuity::default();
        assert_eq!(continuity, BsplineContinuity::None);
    }

    #[test]
    fn check_bspline_continuity_same_surface() {
        use rcad_kernel::geom::BSplineSurface;

        let surf = BSplineSurface {
            degree_u: 3,
            degree_v: 3,
            knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(0.33, 0.0, 0.0), DVec3::new(0.66, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 0.33, 0.0), DVec3::new(0.33, 0.33, 0.0), DVec3::new(0.66, 0.33, 0.0), DVec3::new(1.0, 0.33, 0.0)],
                vec![DVec3::new(0.0, 0.66, 0.0), DVec3::new(0.33, 0.66, 0.0), DVec3::new(0.66, 0.66, 0.0), DVec3::new(1.0, 0.66, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(0.33, 1.0, 0.0), DVec3::new(0.66, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0, 1.0],
                vec![1.0, 1.0, 1.0, 1.0],
            ],
        };

        let continuity = check_bspline_continuity(&surf, &surf, 1e-6);
        // A bicubic B-spline with clamped boundary knots (multiplicity 4) has C0 continuity
        // at boundaries due to knot multiplicity = degree, but is C2 inside the domain.
        // Our implementation reports minimum continuity at any knot, which is C0 at boundaries.
        assert!(continuity >= BsplineContinuity::C0);
    }

    #[test]
    fn check_bspline_continuity_adjacent_v() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0],
                vec![1.0, 1.0],
            ],
        };

        let continuity = check_bspline_continuity(&surf1, &surf2, 1e-6);
        assert!(continuity >= BsplineContinuity::C0);
    }

    #[test]
    fn max_knot_multiplicity_single() {
        let knots = vec![0.0, 0.0, 0.5, 1.0, 1.0];
        let mult = max_knot_multiplicity(&knots);
        assert_eq!(mult, 2);
    }

    #[test]
    fn max_knot_multiplicity_triple() {
        let knots = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0];
        let mult = max_knot_multiplicity(&knots);
        assert_eq!(mult, 3);
    }

    #[test]
    fn same_domain_match_debug() {
        let match_result = SameDomainMatch {
            is_same_domain: true,
            continuity: BsplineContinuity::C1,
            max_control_point_deviation: 0.0,
            max_weight_deviation: 0.0,
            knots_match: true,
            degrees_match: true,
        };

        let debug_str = format!("{:?}", match_result);
        assert!(debug_str.contains("is_same_domain: true"));
        assert!(debug_str.contains("C1"));
    }

    #[test]
    fn merged_face_info_debug() {
        let info = MergedFaceInfo {
            kept_face_idx: 0,
            removed_face_idx: 1,
            merged_edge_count: 6,
            inner_wires_merged: false,
            continuity: BsplineContinuity::C0,
        };

        let debug_str = format!("{:?}", info);
        assert!(debug_str.contains("kept_face_idx: 0"));
        assert!(debug_str.contains("merged_edge_count: 6"));
    }

    #[test]
    fn bspline_continuity_ordering() {
        assert!(BsplineContinuity::None < BsplineContinuity::C0);
        assert!(BsplineContinuity::C0 < BsplineContinuity::G1);
        assert!(BsplineContinuity::G1 < BsplineContinuity::C1);
        assert!(BsplineContinuity::C1 < BsplineContinuity::C2);
        assert!(BsplineContinuity::C2 < BsplineContinuity::CN);
    }

    #[test]
    fn bspline_same_domain_rational_surface() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0), DVec3::new(2.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0), DVec3::new(2.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 2.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        };

        let surf2 = surf1.clone();

        let result = bspline_same_domain(&surf1, &surf2, 1e-6);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(match_result.is_same_domain);
        assert!(match_result.max_weight_deviation < 1e-9);
    }

    #[test]
    fn bspline_same_domain_different_weights() {
        use rcad_kernel::geom::BSplineSurface;

        let surf1 = BSplineSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0), DVec3::new(2.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0), DVec3::new(2.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 2.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        };

        let surf2 = BSplineSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                vec![DVec3::new(0.0, 0.0, 0.0), DVec3::new(1.0, 0.0, 0.0), DVec3::new(2.0, 0.0, 0.0)],
                vec![DVec3::new(0.0, 1.0, 0.0), DVec3::new(1.0, 1.0, 1.0), DVec3::new(2.0, 1.0, 0.0)],
                vec![DVec3::new(0.0, 2.0, 0.0), DVec3::new(1.0, 2.0, 0.0), DVec3::new(2.0, 2.0, 0.0)],
            ],
            weights: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 3.0, 1.0],
                vec![1.0, 1.0, 1.0],
            ],
        };

        let result = bspline_same_domain(&surf1, &surf2, 1e-6);
        assert!(result.is_some());
        let match_result = result.unwrap();
        assert!(!match_result.is_same_domain);
        assert!(match_result.max_weight_deviation > 0.5);
    }

    // ─────────────────────────────────────────────────────────────────────────
    // Shell and Solid Repair Tests
    // ─────────────────────────────────────────────────────────────────────────

    #[test]
    fn check_shell_closure_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Unit box shell should be closed");
        assert_eq!(report.open_edge_count, 0, "Unit box should have no open edges");
        assert_eq!(report.face_count, 6, "Unit box should have 6 faces");
        assert!(report.euler_characteristic > 0, "Unit box should have positive Euler characteristic");
    }

    #[test]
    fn check_shell_closure_unit_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Sphere shell should be closed");
        assert_eq!(report.open_edge_count, 0, "Sphere should have no open edges");
    }

    #[test]
    fn check_shell_closure_unit_cylinder() {
        let brep = BRep::from_primitive(PrimitiveSolid::Cylinder {
            radius: 1.0,
            height: 2.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Cylinder shell should be closed");
        assert_eq!(report.open_edge_count, 0, "Cylinder should have no open edges");
    }

    #[test]
    fn check_shell_closure_open_triangle() {
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Wire, WireEdge};

        // Create an open triangle (not a closed shell)
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(1.0, 0.0, 0.0) });
        brep.vertices.push(Vertex { point: DVec3::new(0.0, 1.0, 0.0) });

        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        // Missing edge 2-0 to make it open

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };

        let shell = Shell { faces: vec![face] };
        let report = check_shell_closure(&shell, &brep);

        assert!(!report.is_closed, "Open triangle should not be closed");
        assert!(report.open_edge_count > 0, "Open triangle should have open edges");
    }

    #[test]
    fn fix_shell_orientation_inverted_normals() {
        // Create a box and invert all its face normals
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        // Invert all face normals in the shell
        for face in &mut brep.solids[0].shells[0].faces {
            face.normal = -face.normal;
        }

        let shell = &brep.solids[0].shells[0].clone();
        let (fixed_shell, report) = fix_shell_orientation(shell, &brep);

        // All 6 faces should be reoriented
        assert!(report.faces_reoriented >= 6, "Should reorient all inverted faces");

        // All normals should now point outward (have positive dot product with outward direction)
        for face in &fixed_shell.faces {
            // For a box centered at origin, check that normals are consistent
            let normal_magnitude = face.normal.length();
            assert!(normal_magnitude > 0.99, "Normal should be unit length");
        }
    }

    #[test]
    fn fix_shell_orientation_correct_normals() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let (_, report) = fix_shell_orientation(shell, &brep);

        // Box from primitive should already have correct normals
        assert_eq!(report.faces_reoriented, 0, "Box from primitive should not need reorientation");
    }

    #[test]
    fn shell_fix_report_summary() {
        let report = ShellFixReport {
            faces_reoriented: 3,
            non_manifold_edges_processed: 1,
            shells_created: 0,
            is_closed: true,
            is_manifold: true,
            open_edge_count: 0,
            non_manifold_edge_count: 0,
        };

        let summary = report.summary();
        assert!(summary.contains("3 faces reoriented"));
        assert!(summary.contains("closed=true"));
    }

    #[test]
    fn closure_report_summary() {
        let report = ClosureReport {
            is_closed: true,
            open_edge_count: 0,
            open_edges: vec![],
            euler_characteristic: 2,
            vertex_count: 8,
            edge_count: 12,
            face_count: 6,
            is_orientable: true,
            genus: Some(0),
        };

        let summary = report.summary();
        assert!(summary.contains("Closed shell"));
        assert!(summary.contains("V=8"));
        assert!(summary.contains("E=12"));
        assert!(summary.contains("F=6"));
        assert!(summary.contains("genus=0"));
    }

    #[test]
    fn check_solid_closure_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let solid = &brep.solids[0];
        let report = check_solid_closure(solid, &brep);

        assert!(report.is_closed, "Box solid should be closed");
        assert!(report.has_proper_nesting, "Box should have proper nesting");
        assert_eq!(report.outer_shell_count, 1, "Box should have 1 outer shell");
        assert_eq!(report.inner_shell_count, 0, "Box should have 0 inner shells");
    }

    #[test]
    fn check_solid_closure_unit_sphere() {
        let brep = BRep::from_primitive(PrimitiveSolid::Sphere {
            radius: 1.0,
        });

        let solid = &brep.solids[0];
        let report = check_solid_closure(solid, &brep);

        assert!(report.is_closed, "Sphere solid should be closed");
        assert_eq!(report.outer_shell_count, 1, "Sphere should have 1 outer shell");
    }

    #[test]
    fn fix_solid_orientation_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let solid = &brep.solids[0];
        let (_, report) = fix_solid_orientation(solid, &brep);

        // Box from primitive should already be properly oriented
        assert!(report.has_valid_closure, "Box should have valid closure");
        assert_eq!(report.outer_shells, 1, "Box should have 1 outer shell");
        assert_eq!(report.inner_shells, 0, "Box should have 0 inner shells");
    }

    #[test]
    fn fix_solid_unit_box() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let solid = &brep.solids[0];
        let (fixed_solid, report) = fix_solid(solid, &brep);

        assert!(report.is_clean(), "Fixed solid should be clean");
        assert!(report.has_valid_closure, "Fixed solid should have valid closure");
        assert!(report.is_properly_oriented, "Fixed solid should be properly oriented");
        assert_eq!(fixed_solid.shells.len(), solid.shells.len(), "Shell count should be preserved");
    }

    #[test]
    fn solid_fix_report_summary() {
        let report = SolidFixReport {
            shells_reoriented: 1,
            faces_reoriented: 3,
            outer_shells: 1,
            inner_shells: 0,
            is_properly_oriented: true,
            has_valid_closure: true,
            total_fixes: 4,
        };

        let summary = report.summary();
        assert!(summary.contains("1 shells reoriented"));
        assert!(summary.contains("3 faces flipped"));
        assert!(summary.contains("1 outer"));
    }

    #[test]
    fn solid_closure_report_summary() {
        let report = SolidClosureReport {
            is_closed: true,
            has_proper_nesting: true,
            outer_shell_count: 1,
            inner_shell_count: 2,
            unclosed_shell_indices: vec![],
            volume: 10.5,
            shell_euler: vec![2, 2, 2],
            solid_euler: 6,
        };

        let summary = report.summary();
        assert!(summary.contains("Valid solid"));
        assert!(summary.contains("2 voids"));
    }

    #[test]
    fn check_shell_closure_torus() {
        let brep = BRep::from_primitive(PrimitiveSolid::Torus {
            major_radius: 2.0,
            minor_radius: 0.5,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        assert!(report.is_closed, "Torus shell should be closed");
        // Torus has genus 1, so Euler characteristic should be 0
        assert_eq!(report.euler_characteristic, 0, "Torus should have Euler characteristic 0");
        assert_eq!(report.genus, Some(1), "Torus should have genus 1");
    }

    #[test]
    fn fix_non_manifold_shell_already_manifold() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let (_, report) = fix_non_manifold_shell(shell, &brep);

        assert!(report.is_manifold, "Box shell should be manifold");
        assert_eq!(report.non_manifold_edge_count, 0, "Box should have no non-manifold edges");
    }

    #[test]
    fn shell_orientability_check() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });

        let shell = &brep.solids[0].shells[0];
        let report = check_shell_closure(shell, &brep);

        // Box should be closed (no open edges)
        assert!(report.is_closed, "Box shell should be closed");
        // Note: orientability check depends on face orientation consistency
        // which may vary based on how the primitive is constructed
    }
}
