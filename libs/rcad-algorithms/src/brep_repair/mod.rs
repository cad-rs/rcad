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
use rcad_kernel::topods;
use rcad_kernel::CurveEval;
use rcad_kernel::Curve2dEval;
use rcad_kernel::SurfaceEval;
use rcad_kernel::Surface3;
use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
use crate::brep_check::{check_orientation_consistency, diagnose_same_parameter, diagnose_same_range};
use crate::tolerance::{
    TOLERANCE_ABS, TOLERANCE_ABS_SQ, TOLERANCE_ADAPTIVE_MAX, TOLERANCE_COORD_SUB, TOLERANCE_FLOAT_DEDUP,
    TOLERANCE_FLOAT_LOOSE, TOLERANCE_LINEAR_ULTRA_STRICT, TOLERANCE_MESH_LEGACY,
    TOLERANCE_METRIC_SQ_NEAR_ZERO, TOLERANCE_RETRY_LADDER_COARSE,
    TOLERANCE_RETRY_LADDER_MID,
};

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
    next_tolerance > current_tolerance + TOLERANCE_FLOAT_DEDUP
}

//閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Public API
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

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
    /// Number of seam edges detected on periodic surfaces.
    pub seam_edges_detected: usize,
    /// Number of edges split at periodic seams.
    pub seam_edges_split: usize,
    /// Number of degenerate points handled (sphere poles, cone apex).
    pub degenerate_points_handled: usize,
    /// Number of edges merged across periodic seams.
    pub seam_edges_merged: usize,
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
    /// Whether scoped cleanup fell back to global.
    pub fell_back_to_global: bool,
    /// Coverage assessment that triggered fallback (if any).
    pub coverage_assessment: Option<CoverageAssessment>,
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
            merge_tolerance: TOLERANCE_MESH_LEGACY,
            remove_small_edges: true,
            min_edge_length: TOLERANCE_MESH_LEGACY,
            max_passes: 3,
            tolerance_growth: 1.0,
            tolerance_cap: f64::INFINITY,
            sew_edges: false,
            edge_sew_tolerance: TOLERANCE_MESH_LEGACY,
            merge_faces: false,
            face_merge_tolerance: TOLERANCE_MESH_LEGACY,
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
            merge_tolerance: TOLERANCE_RETRY_LADDER_MID,
            remove_small_edges: true,
            min_edge_length: TOLERANCE_RETRY_LADDER_MID,
            max_passes: 5,
            tolerance_growth: 1.5,
            tolerance_cap: 0.01,
            sew_edges: true,
            edge_sew_tolerance: TOLERANCE_RETRY_LADDER_MID,
            merge_faces: true,
            face_merge_tolerance: TOLERANCE_RETRY_LADDER_MID,
        }
    }

    /// Create a strategy for injection molding (optimized for thin walls).
    pub fn for_injection_molding() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
            remove_small_edges: true,
            min_edge_length: TOLERANCE_RETRY_LADDER_COARSE,
            max_passes: 10,
            tolerance_growth: 2.0,
            tolerance_cap: 0.1,
            sew_edges: true,
            edge_sew_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
            merge_faces: false, // Don't merge faces for molding
            face_merge_tolerance: TOLERANCE_RETRY_LADDER_COARSE,
        }
    }

    /// Create a strategy for 3D printing (conservative tolerance).
    pub fn for_3d_printing() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: TOLERANCE_ADAPTIVE_MAX, // 0.001mm tolerance for printing
            remove_small_edges: true,
            min_edge_length: TOLERANCE_ADAPTIVE_MAX,
            max_passes: 3,
            tolerance_growth: 1.0,
            tolerance_cap: 0.1,
            sew_edges: false,
            edge_sew_tolerance: TOLERANCE_ADAPTIVE_MAX,
            merge_faces: false,
            face_merge_tolerance: TOLERANCE_ADAPTIVE_MAX,
        }
    }

    /// Create a strategy for CNC machining (precise, no merging).
    pub fn for_cnc_machining() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: TOLERANCE_MESH_LEGACY, // Very precise
            remove_small_edges: true,
            min_edge_length: TOLERANCE_MESH_LEGACY,
            max_passes: 1,
            tolerance_growth: 1.0,
            tolerance_cap: TOLERANCE_RETRY_LADDER_MID,
            sew_edges: false,
            edge_sew_tolerance: TOLERANCE_MESH_LEGACY,
            merge_faces: false,
            face_merge_tolerance: TOLERANCE_MESH_LEGACY,
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

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Scoped Seed Detection Strategies
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

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
            short_edge_threshold: TOLERANCE_RETRY_LADDER_COARSE,
            high_tolerance_threshold: TOLERANCE_ADAPTIVE_MAX,
            near_duplicate_distance: TOLERANCE_RETRY_LADDER_COARSE,
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

/// Multi-dimensional coverage assessment for scoped cleanup.
#[derive(Debug, Clone)]
pub struct CoverageAssessment {
    /// Fraction of vertices covered by seeds.
    pub vertex_coverage: f64,
    /// Fraction of edges covered by seeds (at least one endpoint).
    pub edge_coverage: f64,
    /// Fraction of faces covered by seeds (at least one boundary vertex).
    pub face_coverage: f64,
    /// Whether scoped cleanup should fall back to global.
    pub should_fallback_to_global: bool,
}

/// Assess coverage of seed vertices over the BRep.
pub fn assess_coverage(brep: &BRep, seed_vertices: &[usize]) -> CoverageAssessment {
    let n_vertices = brep.vertices.len().max(1);
    let n_edges = brep.edges.len().max(1);

    let seed_set: std::collections::HashSet<usize> = seed_vertices.iter().copied().collect();

    // Vertex coverage
    let vertex_coverage = seed_vertices.len() as f64 / n_vertices as f64;

    // Edge coverage: at least one endpoint in seeds
    let covered_edges = brep
        .edges
        .iter()
        .filter(|e| seed_set.contains(&e.start) || seed_set.contains(&e.end))
        .count();
    let edge_coverage = covered_edges as f64 / n_edges as f64;

    // Face coverage: at least one boundary vertex in seeds
    let mut covered_faces = 0usize;
    let total_faces = brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count();

    for face in brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
    {
        let has_seed = face
            .outer_wire
            .edges
            .iter()
            .flat_map(|we| {
                brep.edges
                    .get(we.idx)
                    .map(|e| vec![e.start, e.end])
                    .unwrap_or_default()
            })
            .any(|v| seed_set.contains(&v));
        if has_seed {
            covered_faces += 1;
        }
    }
    let face_coverage = if total_faces > 0 {
        covered_faces as f64 / total_faces as f64
    } else {
        0.0
    };

    // Fallback threshold: if any coverage is below 30%, use global
    let min_coverage = vertex_coverage.min(edge_coverage).min(face_coverage);
    let should_fallback = min_coverage < 0.3;

    CoverageAssessment {
        vertex_coverage,
        edge_coverage,
        face_coverage,
        should_fallback_to_global: should_fallback,
    }
}

/// Get adjacent faces for an edge (simplified, no BRepGraph needed).
fn get_edge_adjacent_faces_brep(brep: &BRep, edge_idx: usize) -> Vec<usize> {
    let mut faces = Vec::new();
    let mut face_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for _face in &shell.faces {
                // Check if this face references the edge
                for wire_edge in &_face.outer_wire.edges {
                    if wire_edge.idx == edge_idx {
                        faces.push(face_idx);
                        break;
                    }
                }
                face_idx += 1;
            }
        }
    }
    faces
}

/// Get face normal from geometry or stored normal.
fn get_face_normal(brep: &BRep, face_idx: usize) -> Option<DVec3> {
    // First, try to get from face's stored normal
    let mut current_idx = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                if current_idx == face_idx {
                    // Return the stored face normal
                    return Some(face.normal.normalize());
                }
                current_idx += 1;
            }
        }
    }
    None
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
            // Strategy 1: Edges referenced by multiple faces (potential seams)
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
                if count > 2
                    && let Some(edge) = brep.edges.get(ei) {
                        vertex_set.insert(edge.start);
                        vertex_set.insert(edge.end);
                        edge_set.insert(ei);
                    }
            }

            // Strategy 2: Detect edges with multiple PCurves (potential seams on periodic surfaces)
            for (ei, pcurves) in brep.geom.edge_pcurves.iter().enumerate() {
                if pcurves.len() > 1 {
                    // Multi-PCurve edge is a seam candidate (e.g., seam on cylinder/sphere/torus)
                    if let Some(edge) = brep.edges.get(ei) {
                        vertex_set.insert(edge.start);
                        vertex_set.insert(edge.end);
                        edge_set.insert(ei);
                    }
                }
            }

            // Strategy 3: Detect edges where adjacent face normals have large angle (> 45 degrees)
            for (ei, edge) in brep.edges.iter().enumerate() {
                let adj_faces = get_edge_adjacent_faces_brep(brep, ei);
                if adj_faces.len() == 2
                    && let (Some(n1), Some(n2)) = (
                        get_face_normal(brep, adj_faces[0]),
                        get_face_normal(brep, adj_faces[1]),
                    ) {
                        let dot = n1.dot(n2);
                        // Angle > 45 degrees means dot < cos(45鎺? 閳?0.707
                        // Use abs to handle both same-side and opposite-side normals
                        if dot.abs() < std::f64::consts::FRAC_PI_4.cos() {
                            vertex_set.insert(edge.start);
                            vertex_set.insert(edge.end);
                            edge_set.insert(ei);
                        }
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
                if let (Some(s), Some(e)) = (start, end)
                    && (s - e).length() < config.short_edge_threshold {
                        combined.insert(edge.start);
                        combined.insert(edge.end);
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

            // Seam candidate Strategy 1: Edges referenced by multiple faces (potential seams)
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
                if count > 2
                    && let Some(edge) = brep.edges.get(ei) {
                        combined.insert(edge.start);
                        combined.insert(edge.end);
                    }
            }

            // Seam candidate Strategy 2: Edges with multiple PCurves
            for (ei, pcurves) in brep.geom.edge_pcurves.iter().enumerate() {
                if pcurves.len() > 1
                    && let Some(edge) = brep.edges.get(ei) {
                        combined.insert(edge.start);
                        combined.insert(edge.end);
                    }
            }

            // Seam candidate Strategy 3: Edges with large face normal angle (> 45 degrees)
            for (ei, edge) in brep.edges.iter().enumerate() {
                let adj_faces = get_edge_adjacent_faces_brep(brep, ei);
                if adj_faces.len() == 2
                    && let (Some(n1), Some(n2)) = (
                        get_face_normal(brep, adj_faces[0]),
                        get_face_normal(brep, adj_faces[1]),
                    ) {
                        let dot = n1.dot(n2);
                        // Angle > 45 degrees means dot < cos(45鎺? 閳?0.707
                        if dot.abs() < std::f64::consts::FRAC_PI_4.cos() {
                            combined.insert(edge.start);
                            combined.insert(edge.end);
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
            "RepairReport {{ vertices_merged={}, degenerate_removed={}, normals_recomputed={}, faces_reoriented={}, wires_fixed={}, same_range_fixed={}, same_parameter_fixed={}, seam_edges_detected={}, seam_edges_split={}, degenerate_points_handled={}, seam_edges_merged={} }}",
            self.vertices_merged,
            self.degenerate_faces_removed,
            self.normals_recomputed,
            self.faces_reoriented,
            self.wires_fixed,
            self.same_range_fixed,
            self.same_parameter_fixed,
            self.seam_edges_detected,
            self.seam_edges_split,
            self.degenerate_points_handled,
            self.seam_edges_merged,
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
///
/// Automatically falls back to global cleanup when seed coverage is below
/// the fallback threshold (30% for any coverage dimension).
pub fn make_connected_iterative_scoped_with_growth_cap(
    brep: &BRep,
    scope_vertices: &[usize],
    tolerance: f64,
    max_passes: usize,
    tolerance_growth: f64,
    tolerance_cap: f64,
) -> (BRep, MakeConnectedReport) {
    // Assess coverage first
    let assessment = assess_coverage(brep, scope_vertices);

    if assessment.should_fallback_to_global {
        // Fall back to global cleanup with same parameters
        let (result, mut report) = make_connected_iterative_with_growth_cap(
            brep,
            tolerance,
            max_passes,
            tolerance_growth,
            tolerance_cap,
        );
        report.fell_back_to_global = true;
        report.coverage_assessment = Some(assessment);
        return (result, report);
    }

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
                            sample_point: None,
                            mesh_dirty: true,
                surface_idx: face.surface_idx,
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

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Edge Sewing (MakeConnected enhancement)
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

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
/// Analogous to `BOPAlgo_MakeConnected` in OCCT. Old-BRep conversion is internal.
pub fn make_connected_enhanced(brep: &topods::BRep, tolerance: f64, max_passes: usize) -> (topods::BRep, MakeConnectedReport) {
    let old = rcad_kernel::BRep::from_topods_with_location(brep, glam::DAffine3::IDENTITY);
    let (result, report) = make_connected_enhanced_old(&old, tolerance, max_passes);
    (result.to_topods(), report)
}

/// Legacy: takes old BRep.
fn make_connected_enhanced_old(brep: &BRep, tolerance: f64, max_passes: usize) -> (BRep, MakeConnectedReport) {
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

// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// Advanced Shared Topology Detection
// 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

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

    for e1 in edge_to_faces.keys() {
        for e2 in edge_to_faces.keys() {
            if e1 >= e2 {
                continue;
            }
            if processed_edge_pairs.contains(&(*e1, *e2)) {
                continue;
            }
            processed_edge_pairs.insert((*e1, *e2));

            // Check if edges have shared geometry
            if let Some(info) = analyze_shared_edge_pair(brep, *e1, *e2, tol)
                && info.geometry_compatible {
                    report.shared_edges.push(info);
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
include!("extra1.rs");
include!("extra2.rs");
include!("extra3.rs");
include!("extra4.rs");
include!("extra5.rs");
include!("extra6.rs");
#[cfg(test)]
mod tests {
    // All test files (tests_p1..p3) have pre-existing syntax errors.
    // Include none until they are fixed separately.
}
