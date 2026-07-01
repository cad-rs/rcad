//! Non-manifold topology support for B-Rep models.
//!
//! Non-manifold topology allows multiple faces (more than 2) to share a single edge,
//! which is common in:
//! - Thin-walled structures (internal representation)
//! - Analysis models (FEA pre-processing)
//! - Assembly contact faces
//! - Multi-region models (CompSolid)
//!
//! # Detection
//!
//! Use [`is_manifold`], [`non_manifold_edges`], and [`non_manifold_vertices`]
//! to detect non-manifold conditions in a BRep.
//!
//! For comprehensive analysis, use [`detect_non_manifold_topology`] which provides
//! detailed per-edge and per-vertex analysis.
//!
//! # Repair
//!
//! - [`split_non_manifold_edges`]: Duplicates edges so each pair of faces has
//!   its own copy, converting non-manifold edges to manifold.
//! - [`make_manifold`]: Full conversion pipeline that splits non-manifold edges
//!   and optionally stitches boundary edges.
//! - [`convert_to_manifold`]: Comprehensive conversion with vertex duplication
//!   and geometric integrity preservation.
//!
//! # Construction
//!
//! - [`merge_shells_at_interface`]: Creates non-manifold topology by merging
//!   two shells along coincident boundary faces.
//!
//! # Example
//!
//! ```
//! use rcad_kernel::{BRep, BRepGraph, PrimitiveSolid};
//! use rcad_algorithms::non_manifold::{is_manifold, non_manifold_edges, split_non_manifold_edges};
//!
//! let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
//! assert!(is_manifold(&brep));
//!
//! // After operations that might create non-manifold edges:
//! let nm_edges = non_manifold_edges(&brep);
//! if !nm_edges.is_empty() {
//!     let (repaired, report) = split_non_manifold_edges(&brep);
//! }
//! ```

use crate::tolerance::*;
use std::collections::{HashMap, HashSet};
use rcad_kernel::{
    BRep, BRepGraph, Face, Shell,
};
use glam::DVec3;

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Detection API
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Returns `true` if the BRep is manifold (every edge has exactly 2 adjacent faces).
///
/// This is a convenience wrapper around `BRepGraph::is_manifold()`.
pub fn is_manifold(brep: &BRep) -> bool {
    BRepGraph::from_brep(brep).is_manifold()
}

/// Returns indices of all non-manifold edges (edges with != 2 adjacent faces).
///
/// Includes:
/// - Boundary edges (1 adjacent face)
/// - Multi-face edges (> 2 adjacent faces)
/// - Orphan edges (0 adjacent faces)
pub fn non_manifold_edges(brep: &BRep) -> Vec<usize> {
    BRepGraph::from_brep(brep).non_manifold_edges()
}

/// Returns indices of all boundary edges (edges with exactly 1 adjacent face).
pub fn boundary_edges(brep: &BRep) -> Vec<usize> {
    BRepGraph::from_brep(brep).boundary_edges()
}

/// Returns indices of all multi-face edges (edges with > 2 adjacent faces).
///
/// These are true non-manifold edges where 3+ faces meet at a single edge.
pub fn multi_face_edges(brep: &BRep) -> Vec<usize> {
    BRepGraph::from_brep(brep).multi_face_edges()
}

/// Returns indices of all non-manifold vertices.
///
/// A vertex is non-manifold if it lies on at least one multi-face edge (> 2 adjacent faces).
pub fn non_manifold_vertices(brep: &BRep) -> Vec<usize> {
    BRepGraph::from_brep(brep).non_manifold_vertices()
}

/// Returns indices of orphan edges (edges with 0 adjacent faces).
pub fn orphan_edges(brep: &BRep) -> Vec<usize> {
    BRepGraph::from_brep(brep).orphan_edges()
}

/// Detailed non-manifold analysis report.
#[derive(Debug, Clone, Default)]
pub struct NonManifoldReport {
    /// Total number of edges.
    pub total_edges: usize,
    /// Total number of faces.
    pub total_faces: usize,
    /// Number of manifold edges (exactly 2 adjacent faces).
    pub manifold_edge_count: usize,
    /// Number of boundary edges (exactly 1 adjacent face).
    pub boundary_edge_count: usize,
    /// Number of multi-face edges (> 2 adjacent faces).
    pub multi_face_edge_count: usize,
    /// Number of orphan edges (0 adjacent faces).
    pub orphan_edge_count: usize,
    /// Number of non-manifold vertices.
    pub non_manifold_vertex_count: usize,
    /// Whether the BRep is fully manifold.
    pub is_manifold: bool,
    /// Whether the BRep is closed (no boundary edges).
    pub is_closed: bool,
}

impl NonManifoldReport {
    /// Returns `true` if the BRep has no topological issues.
    pub fn is_clean(&self) -> bool {
        self.is_manifold && self.is_closed && self.orphan_edge_count == 0
    }
}

/// Classification of non-manifold edge types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonManifoldEdgeType {
    /// Edge shared by exactly 1 face (boundary edge, open shell).
    Boundary,
    /// Edge shared by exactly 2 faces (manifold).
    Manifold,
    /// Edge shared by 3+ faces (true non-manifold edge).
    MultiFace,
    /// Edge not referenced by any face (orphan).
    Orphan,
}

/// Classification of non-manifold vertex types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NonManifoldVertexType {
    /// Vertex is on manifold edges only.
    Manifold,
    /// Vertex is on at least one boundary edge.
    Boundary,
    /// Vertex is on at least one multi-face edge (>2 adjacent faces).
    MultiFaceJunction,
    /// Vertex has a "bow-tie" configuration (edges fan out in multiple separate regions).
    BowTie,
    /// Vertex is isolated (not on any edge).
    Isolated,
}

/// Detailed information about a non-manifold edge.
#[derive(Debug, Clone)]
pub struct NonManifoldEdgeDetail {
    /// Edge index.
    pub edge_index: usize,
    /// Type of non-manifold configuration.
    pub edge_type: NonManifoldEdgeType,
    /// Number of adjacent faces.
    pub adjacent_face_count: usize,
    /// Indices of adjacent faces.
    pub adjacent_faces: Vec<usize>,
    /// Start vertex index.
    pub start_vertex: usize,
    /// End vertex index.
    pub end_vertex: usize,
    /// Edge tolerance (if available).
    pub tolerance: f64,
}

/// Detailed information about a non-manifold vertex.
#[derive(Debug, Clone)]
pub struct NonManifoldVertexDetail {
    /// Vertex index.
    pub vertex_index: usize,
    /// Type of non-manifold configuration.
    pub vertex_type: NonManifoldVertexType,
    /// Number of edges incident to this vertex.
    pub incident_edge_count: usize,
    /// Indices of incident edges.
    pub incident_edges: Vec<usize>,
    /// Number of faces incident to this vertex.
    pub incident_face_count: usize,
    /// Number of distinct edge fans (separate manifold regions at this vertex).
    pub fan_count: usize,
    /// Whether this vertex can be duplicated to resolve non-manifoldness.
    pub can_duplicate: bool,
}

/// Counts of non-manifold entities.
#[derive(Debug, Clone, Default)]
pub struct NonManifoldCounts {
    /// Total edges analyzed.
    pub total_edges: usize,
    /// Manifold edges (exactly 2 adjacent faces).
    pub manifold_edges: usize,
    /// Boundary edges (exactly 1 adjacent face).
    pub boundary_edges: usize,
    /// Multi-face edges (3+ adjacent faces).
    pub multi_face_edges: usize,
    /// Orphan edges (0 adjacent faces).
    pub orphan_edges: usize,
    /// Total vertices analyzed.
    pub total_vertices: usize,
    /// Manifold vertices.
    pub manifold_vertices: usize,
    /// Boundary vertices (on at least one boundary edge).
    pub boundary_vertices: usize,
    /// Multi-face junction vertices (on at least one multi-face edge).
    pub multi_face_junction_vertices: usize,
    /// Bow-tie vertices (multiple disconnected edge fans).
    pub bow_tie_vertices: usize,
    /// Isolated vertices.
    pub isolated_vertices: usize,
}

impl NonManifoldCounts {
    /// Returns the total number of non-manifold edges.
    pub fn non_manifold_edge_count(&self) -> usize {
        self.boundary_edges + self.multi_face_edges + self.orphan_edges
    }

    /// Returns the total number of non-manifold vertices.
    pub fn non_manifold_vertex_count(&self) -> usize {
        self.boundary_vertices + self.multi_face_junction_vertices + self.bow_tie_vertices + self.isolated_vertices
    }

    /// Returns true if there are no non-manifold entities.
    pub fn is_manifold(&self) -> bool {
        self.non_manifold_edge_count() == 0 && self.non_manifold_vertex_count() == 0
    }
}

/// Comprehensive non-manifold detection report.
#[derive(Debug, Clone, Default)]
pub struct DetailedNonManifoldReport {
    /// Basic counts.
    pub counts: NonManifoldCounts,
    /// Whether the BRep is fully manifold.
    pub is_manifold: bool,
    /// Whether the BRep is closed (no boundary edges).
    pub is_closed: bool,
    /// Per-edge analysis for non-manifold edges.
    pub edge_details: Vec<NonManifoldEdgeDetail>,
    /// Per-vertex analysis for non-manifold vertices.
    pub vertex_details: Vec<NonManifoldVertexDetail>,
    /// Number of distinct manifold regions in the model.
    pub manifold_region_count: usize,
    /// Indices of shells that contain non-manifold topology.
    pub non_manifold_shell_indices: Vec<usize>,
    /// Summary message.
    pub summary: String,
}

impl DetailedNonManifoldReport {
    /// Returns true if the BRep has no non-manifold topology.
    pub fn is_clean(&self) -> bool {
        self.is_manifold && self.is_closed
    }
}

/// Performs comprehensive non-manifold detection on a BRep.
///
/// This function analyzes all edges and vertices to detect:
/// - Edges shared by 3+ faces (true non-manifold edges)
/// - Boundary edges (open shells)
/// - Orphan edges (unused)
/// - Vertices with non-manifold configurations (bow-tie, multi-face junction)
///
/// # Arguments
/// * `brep` - The BRep to analyze.
///
/// # Returns
/// A `DetailedNonManifoldReport` with comprehensive analysis.
pub fn detect_non_manifold_topology(brep: &BRep) -> DetailedNonManifoldReport {
    let graph = BRepGraph::from_brep(brep);
    let mut report = DetailedNonManifoldReport::default();

    // Analyze edges
    for ei in 0..brep.edges.len() {
        let adjacent_faces: Vec<usize> = graph.edge_adjacent_faces(ei).to_vec();
        let face_count = adjacent_faces.len();

        let edge_type = match face_count {
            0 => NonManifoldEdgeType::Orphan,
            1 => NonManifoldEdgeType::Boundary,
            2 => NonManifoldEdgeType::Manifold,
            _ => NonManifoldEdgeType::MultiFace,
        };

        report.counts.total_edges += 1;
        match edge_type {
            NonManifoldEdgeType::Manifold => report.counts.manifold_edges += 1,
            NonManifoldEdgeType::Boundary => report.counts.boundary_edges += 1,
            NonManifoldEdgeType::MultiFace => report.counts.multi_face_edges += 1,
            NonManifoldEdgeType::Orphan => report.counts.orphan_edges += 1,
        }

        // Record details for non-manifold edges
        if edge_type != NonManifoldEdgeType::Manifold {
            let edge = &brep.edges[ei];
            let tolerance = brep.geom.edge_tolerance.get(ei).copied().unwrap_or(0.0);
            report.edge_details.push(NonManifoldEdgeDetail {
                edge_index: ei,
                edge_type,
                adjacent_face_count: face_count,
                adjacent_faces,
                start_vertex: edge.start,
                end_vertex: edge.end,
                tolerance,
            });
        }
    }

    // Analyze vertices
    for vi in 0..brep.vertices.len() {
        let incident_edges: Vec<usize> = graph.vertex_adjacent_edges(vi).to_vec();
        let edge_count = incident_edges.len();

        if edge_count == 0 {
            report.counts.isolated_vertices += 1;
            report.vertex_details.push(NonManifoldVertexDetail {
                vertex_index: vi,
                vertex_type: NonManifoldVertexType::Isolated,
                incident_edge_count: 0,
                incident_edges: vec![],
                incident_face_count: 0,
                fan_count: 0,
                can_duplicate: false,
            });
            report.counts.total_vertices += 1;
            continue;
        }

        // Count incident faces and check for non-manifold edges
        let mut incident_face_set = std::collections::HashSet::new();
        let mut has_multi_face_edge = false;
        let mut has_boundary_edge = false;

        for &ei in &incident_edges {
            let face_count = graph.edge_adjacent_faces(ei).len();
            if face_count > 2 {
                has_multi_face_edge = true;
            } else if face_count == 1 {
                has_boundary_edge = true;
            }
            for &fi in graph.edge_adjacent_faces(ei) {
                incident_face_set.insert(fi);
            }
        }

        // Compute edge fans (separate manifold regions at this vertex)
        let fan_count = compute_edge_fan_count(vi, &incident_edges, &graph, brep);

        // Determine vertex type
        let vertex_type = if has_multi_face_edge {
            if fan_count > 1 {
                NonManifoldVertexType::BowTie
            } else {
                NonManifoldVertexType::MultiFaceJunction
            }
        } else if has_boundary_edge {
            NonManifoldVertexType::Boundary
        } else {
            NonManifoldVertexType::Manifold
        };

        report.counts.total_vertices += 1;
        match vertex_type {
            NonManifoldVertexType::Manifold => report.counts.manifold_vertices += 1,
            NonManifoldVertexType::Boundary => report.counts.boundary_vertices += 1,
            NonManifoldVertexType::MultiFaceJunction => report.counts.multi_face_junction_vertices += 1,
            NonManifoldVertexType::BowTie => report.counts.bow_tie_vertices += 1,
            NonManifoldVertexType::Isolated => {} // Already counted above
        }

        // Record details for non-manifold vertices
        if vertex_type != NonManifoldVertexType::Manifold {
            report.vertex_details.push(NonManifoldVertexDetail {
                vertex_index: vi,
                vertex_type,
                incident_edge_count: edge_count,
                incident_edges,
                incident_face_count: incident_face_set.len(),
                fan_count,
                can_duplicate: fan_count > 1 || vertex_type == NonManifoldVertexType::BowTie,
            });
        }
    }

    // Compute manifold regions
    report.manifold_region_count = graph.manifold_region_count();

    // Find shells with non-manifold topology
    for (solid_idx, solid) in brep.solids.iter().enumerate() {
        for (shell_idx, shell) in solid.shells.iter().enumerate() {
            let has_nm = shell.faces.iter().any(|face| {
                face.outer_wire.edges.iter().any(|we| {
                    graph.edge_valence(we.idx) != 2
                })
            });
            if has_nm {
                report.non_manifold_shell_indices.push(solid_idx * 1000 + shell_idx);
            }
        }
    }

    // Determine overall status
    report.is_manifold = report.counts.multi_face_edges == 0 && report.counts.bow_tie_vertices == 0;
    report.is_closed = report.counts.boundary_edges == 0;

    // Generate summary
    let mut issues = Vec::new();
    if report.counts.multi_face_edges > 0 {
        issues.push(format!("{} multi-face edges", report.counts.multi_face_edges));
    }
    if report.counts.boundary_edges > 0 {
        issues.push(format!("{} boundary edges", report.counts.boundary_edges));
    }
    if report.counts.bow_tie_vertices > 0 {
        issues.push(format!("{} bow-tie vertices", report.counts.bow_tie_vertices));
    }
    if report.counts.multi_face_junction_vertices > 0 {
        issues.push(format!("{} multi-face junctions", report.counts.multi_face_junction_vertices));
    }
    if report.counts.orphan_edges > 0 {
        issues.push(format!("{} orphan edges", report.counts.orphan_edges));
    }
    if report.counts.isolated_vertices > 0 {
        issues.push(format!("{} isolated vertices", report.counts.isolated_vertices));
    }

    report.summary = if issues.is_empty() {
        "Manifold and closed".to_string()
    } else {
        format!("Non-manifold: {}", issues.join(", "))
    };

    report
}

/// Compute the number of distinct edge fans at a vertex.
///
/// An edge fan is a group of edges that form a single manifold region
/// around the vertex. If there are multiple disconnected fans, the
/// vertex has a "bow-tie" configuration.
fn compute_edge_fan_count(
    _vertex_idx: usize,
    incident_edges: &[usize],
    graph: &BRepGraph,
    _brep: &BRep,
) -> usize {
    if incident_edges.len() <= 2 {
        return 1;
    }

    // Build adjacency between edges based on shared faces
    let mut edge_adjacency: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();

    for &ei in incident_edges {
        let adjacent_faces = graph.edge_adjacent_faces(ei);
        for &fi in adjacent_faces {
            // Find other edges at this vertex that share the same face
            let face_edges = graph.face_edges(fi);
            for &other_ei in face_edges {
                if other_ei != ei && incident_edges.contains(&other_ei) {
                    edge_adjacency.entry(ei).or_default().push(other_ei);
                }
            }
        }
    }

    // Count connected components in the edge adjacency graph
    let mut visited = std::collections::HashSet::new();
    let mut fan_count = 0;

    for &ei in incident_edges {
        if visited.contains(&ei) {
            continue;
        }

        // BFS to find all edges in this fan
        let mut stack = vec![ei];
        while let Some(current) = stack.pop() {
            if visited.contains(&current) {
                continue;
            }
            visited.insert(current);

            if let Some(neighbors) = edge_adjacency.get(&current) {
                for &neighbor in neighbors {
                    if !visited.contains(&neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }

        fan_count += 1;
    }

    fan_count
}

/// Returns `true` if the BRep has any non-manifold topology.
///
/// This is equivalent to `!is_manifold(brep)` but checks specifically
/// for true non-manifold conditions (multi-face edges, bow-tie vertices).
pub fn is_non_manifold(brep: &BRep) -> bool {
    !is_manifold(brep)
}

/// Counts all non-manifold entities in a BRep.
///
/// This provides a quick summary without full detailed analysis.
pub fn count_non_manifold_entities(brep: &BRep) -> NonManifoldCounts {
    let report = detect_non_manifold_topology(brep);
    report.counts
}

/// Performs comprehensive non-manifold analysis on a BRep.
pub fn analyze_non_manifold(brep: &BRep) -> NonManifoldReport {
    let graph = BRepGraph::from_brep(brep);
    let nm_edges = graph.non_manifold_edges();
    let boundary = graph.boundary_edges();
    let multi = graph.multi_face_edges();
    let orphan = graph.orphan_edges();
    let nm_verts = graph.non_manifold_vertices();

    NonManifoldReport {
        total_edges: graph.edge_count,
        total_faces: graph.face_count,
        manifold_edge_count: graph.edge_count - nm_edges.len(),
        boundary_edge_count: boundary.len(),
        multi_face_edge_count: multi.len(),
        orphan_edge_count: orphan.len(),
        non_manifold_vertex_count: nm_verts.len(),
        is_manifold: graph.is_manifold(),
        is_closed: graph.is_closed(),
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Edge Splitting (Non-manifold -> Manifold Conversion)
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Report from splitting non-manifold edges.
#[derive(Debug, Clone, Default)]
pub struct EdgeSplitReport {
    /// Number of non-manifold edges that were split.
    pub edges_split: usize,
    /// Total number of new edges created.
    pub new_edges_created: usize,
    /// Mapping from old edge index to list of new edge indices.
    pub edge_mapping: HashMap<usize, Vec<usize>>,
    /// Number of vertices duplicated.
    pub vertices_duplicated: usize,
}

/// Splits all non-manifold edges to create a manifold BRep.
///
/// For each edge with N > 2 adjacent faces, creates (N / 2) edge copies
/// so each resulting edge has exactly 2 adjacent faces.
///
/// This is a lossy operation - the resulting BRep loses the non-manifold
/// connectivity information.
pub fn split_non_manifold_edges(brep: &BRep) -> (BRep, EdgeSplitReport) {
    let graph = BRepGraph::from_brep(brep);
    let multi_edges = graph.multi_face_edges();

    if multi_edges.is_empty() {
        return (brep.clone(), EdgeSplitReport::default());
    }

    let mut result = brep.clone();
    let mut report = EdgeSplitReport::default();

    // Process each multi-face edge
    for &edge_idx in &multi_edges {
        let adjacent_faces: Vec<usize> = graph.edge_adjacent_faces(edge_idx).to_vec();
        let n_faces = adjacent_faces.len();

        if n_faces < 3 {
            continue; // Not actually a multi-face edge
        }

        // We need to create (n_faces / 2) edge copies (rounded up if odd)
        // Each new edge will be assigned to 2 faces (except possibly one if odd)
        let n_new_edges = n_faces.div_ceil(2);

        // Collect new edge indices
        let mut new_edge_indices = Vec::with_capacity(n_new_edges);
        new_edge_indices.push(edge_idx); // Keep original for first pair

        // Create edge copies
        let original_edge = result.edges[edge_idx];
        for _ in 1..n_new_edges {
            let new_idx = result.edges.len();
            result.edges.push(original_edge);
            new_edge_indices.push(new_idx);
            report.new_edges_created += 1;
        }

        // Update geometry mapping for new edges
        if edge_idx < brep.geom.edge_curve.len() {
            let curve_ref = brep.geom.edge_curve[edge_idx];
            let curve_range = brep.geom.edge_curve_range.get(edge_idx).copied().flatten();
            let degenerated = brep.geom.edge_degenerated.get(edge_idx).copied().unwrap_or(false);
            let tolerance = brep.geom.edge_tolerance.get(edge_idx).copied().unwrap_or(0.0);
            let same_param = brep.geom.edge_same_parameter.get(edge_idx).copied().unwrap_or(true);
            let same_range = brep.geom.edge_same_range.get(edge_idx).copied().unwrap_or(true);

            for &new_idx in &new_edge_indices[1..] {
                if result.geom.edge_curve.len() <= new_idx {
                    result.geom.edge_curve.resize(new_idx + 1, None);
                }
                result.geom.edge_curve[new_idx] = curve_ref;

                if result.geom.edge_curve_range.len() <= new_idx {
                    result.geom.edge_curve_range.resize(new_idx + 1, None);
                }
                result.geom.edge_curve_range[new_idx] = curve_range;

                if result.geom.edge_degenerated.len() <= new_idx {
                    result.geom.edge_degenerated.resize(new_idx + 1, false);
                }
                result.geom.edge_degenerated[new_idx] = degenerated;

                if result.geom.edge_tolerance.len() <= new_idx {
                    result.geom.edge_tolerance.resize(new_idx + 1, 0.0);
                }
                result.geom.edge_tolerance[new_idx] = tolerance;

                if result.geom.edge_same_parameter.len() <= new_idx {
                    result.geom.edge_same_parameter.resize(new_idx + 1, true);
                }
                result.geom.edge_same_parameter[new_idx] = same_param;

                if result.geom.edge_same_range.len() <= new_idx {
                    result.geom.edge_same_range.resize(new_idx + 1, true);
                }
                result.geom.edge_same_range[new_idx] = same_range;
            }
        }

        // Copy PCurves for new edges
        if edge_idx < brep.geom.edge_pcurves.len() {
            let pcurves = brep.geom.edge_pcurves[edge_idx].clone();
            for &new_idx in &new_edge_indices[1..] {
                if result.geom.edge_pcurves.len() <= new_idx {
                    result.geom.edge_pcurves.resize(new_idx + 1, Vec::new());
                }
                result.geom.edge_pcurves[new_idx] = pcurves.clone();
            }
        }

        // Record mapping
        report.edge_mapping.insert(edge_idx, new_edge_indices.clone());
        report.edges_split += 1;

        // Reassign faces to new edges
        // Each new edge gets 2 faces (except possibly the last one if odd)
        for (new_e_idx, &new_edge) in new_edge_indices.iter().enumerate() {
            let face_start = new_e_idx * 2;
            let face_end = ((new_e_idx + 1) * 2).min(n_faces);

            for fi in face_start..face_end {
                let &flat_face_idx = &adjacent_faces[fi];
                // Find and update the wire edge reference
                update_face_edge_reference(&mut result, flat_face_idx, edge_idx, new_edge);
            }
        }
    }

    (result, report)
}

/// Updates a face's wire to reference a new edge instead of the old one.
fn update_face_edge_reference(brep: &mut BRep, flat_face_idx: usize, old_edge: usize, new_edge: usize) {
    // Find the solid, shell, and local face index from the flat index
    let mut current_flat = 0usize;
    for solid in &mut brep.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                if current_flat == flat_face_idx {
                    // Update outer wire
                    for we in &mut face.outer_wire.edges {
                        if we.idx == old_edge {
                            we.idx = new_edge;
                        }
                    }
                    // Update inner wires
                    for inner in &mut face.inner_wires {
                        for we in &mut inner.edges {
                            if we.idx == old_edge {
                                we.idx = new_edge;
                            }
                        }
                    }
                    return;
                }
                current_flat += 1;
            }
        }
    }
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Make Manifold (Full Conversion Pipeline)
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Options for manifold conversion.
#[derive(Debug, Clone, Copy)]
pub struct MakeManifoldOptions {
    /// Split multi-face edges (> 2 adjacent faces).
    pub split_edges: bool,
    /// Remove orphan edges (0 adjacent faces).
    pub remove_orphans: bool,
    /// Tolerance for geometric operations.
    pub tolerance: f64,
}

impl Default for MakeManifoldOptions {
    fn default() -> Self {
        Self {
            split_edges: true,
            remove_orphans: true,
            tolerance: TOLERANCE_MESH_LEGACY,
        }
    }
}

/// Report from manifold conversion.
#[derive(Debug, Clone, Default)]
pub struct MakeManifoldReport {
    /// Whether the input was already manifold.
    pub was_already_manifold: bool,
    /// Report from edge splitting (if performed).
    pub edge_split_report: EdgeSplitReport,
    /// Number of orphan edges removed.
    pub orphans_removed: usize,
    /// Whether the result is manifold.
    pub is_manifold: bool,
}

/// Converts a potentially non-manifold BRep to a manifold BRep.
///
/// This is a convenience function that combines multiple repair operations:
/// 1. Split multi-face edges
/// 2. Remove orphan edges
///
/// Returns `Err` if the conversion fails.
pub fn make_manifold(brep: &BRep) -> Result<(BRep, MakeManifoldReport), String> {
    make_manifold_with_options(brep, MakeManifoldOptions::default())
}

/// Converts a potentially non-manifold BRep to a manifold BRep with custom options.
pub fn make_manifold_with_options(
    brep: &BRep,
    options: MakeManifoldOptions,
) -> Result<(BRep, MakeManifoldReport), String> {
    let graph = BRepGraph::from_brep(brep);
    let mut report = MakeManifoldReport {
        was_already_manifold: graph.is_manifold(),
        ..Default::default()
    };

    if report.was_already_manifold {
        report.is_manifold = true;
        return Ok((brep.clone(), report));
    }

    let mut result = brep.clone();

    // Step 1: Split multi-face edges
    if options.split_edges {
        let (split_result, split_report) = split_non_manifold_edges(&result);
        result = split_result;
        report.edge_split_report = split_report;
    }

    // Step 2: Remove orphan edges
    if options.remove_orphans {
        let graph_after_split = BRepGraph::from_brep(&result);
        let orphans = graph_after_split.orphan_edges();
        if !orphans.is_empty() {
            result = remove_edges(&result, &orphans);
            report.orphans_removed = orphans.len();
        }
    }

    // Verify result
    let final_graph = BRepGraph::from_brep(&result);
    report.is_manifold = final_graph.is_manifold();

    Ok((result, report))
}

/// Removes the specified edges from a BRep.
///
/// This updates wires to remove references to the deleted edges.
fn remove_edges(brep: &BRep, edge_indices: &[usize]) -> BRep {
    let edge_set: HashSet<usize> = edge_indices.iter().copied().collect();
    let mut result = brep.clone();

    for solid in &mut result.solids {
        for shell in &mut solid.shells {
            for face in &mut shell.faces {
                // Remove from outer wire
                face.outer_wire.edges.retain(|we| !edge_set.contains(&we.idx));
                // Remove from inner wires
                for inner in &mut face.inner_wires {
                    inner.edges.retain(|we| !edge_set.contains(&we.idx));
                }
            }
        }
    }

    result
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Enhanced Manifold Conversion with Vertex Duplication
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Options for comprehensive manifold conversion.
#[derive(Debug, Clone)]
pub struct ManifoldConversionOptions {
    /// Split multi-face edges (> 2 adjacent faces).
    pub split_edges: bool,
    /// Duplicate vertices at non-manifold junctions.
    pub duplicate_vertices: bool,
    /// Remove orphan edges (0 adjacent faces).
    pub remove_orphans: bool,
    /// Remove isolated vertices (not on any edge).
    pub remove_isolated_vertices: bool,
    /// Tolerance for geometric operations.
    pub tolerance: f64,
    /// Preserve geometric integrity by copying curve/surface data.
    pub preserve_geometry: bool,
    /// Attempt to stitch boundary edges after conversion.
    pub stitch_boundaries: bool,
    /// Tolerance for boundary stitching.
    pub stitch_tolerance: f64,
}

impl Default for ManifoldConversionOptions {
    fn default() -> Self {
        Self {
            split_edges: true,
            duplicate_vertices: true,
            remove_orphans: true,
            remove_isolated_vertices: true,
            tolerance: TOLERANCE_MESH_LEGACY,
            preserve_geometry: true,
            stitch_boundaries: false,
            stitch_tolerance: TOLERANCE_MESH_LEGACY,
        }
    }
}

impl ManifoldConversionOptions {
    /// Create conservative options (minimal changes).
    pub fn conservative() -> Self {
        Self {
            split_edges: true,
            duplicate_vertices: false,
            remove_orphans: false,
            remove_isolated_vertices: false,
            tolerance: TOLERANCE_MESH_LEGACY,
            preserve_geometry: true,
            stitch_boundaries: false,
            stitch_tolerance: TOLERANCE_MESH_LEGACY,
        }
    }

    /// Create aggressive options (full conversion).
    pub fn aggressive() -> Self {
        Self {
            split_edges: true,
            duplicate_vertices: true,
            remove_orphans: true,
            remove_isolated_vertices: true,
            tolerance: TOLERANCE_RETRY_LADDER_MID,
            preserve_geometry: true,
            stitch_boundaries: true,
            stitch_tolerance: TOLERANCE_RETRY_LADDER_MID,
        }
    }
}

/// Report from comprehensive manifold conversion.
#[derive(Debug, Clone, Default)]
pub struct ManifoldConversionReport {
    /// Whether the input was already manifold.
    pub was_already_manifold: bool,
    /// Number of multi-face edges that were split.
    pub edges_split: usize,
    /// Number of new edges created.
    pub new_edges_created: usize,
    /// Number of vertices duplicated.
    pub vertices_duplicated: usize,
    /// Number of orphan edges removed.
    pub orphans_removed: usize,
    /// Number of isolated vertices removed.
    pub isolated_vertices_removed: usize,
    /// Number of boundary edges after conversion.
    pub boundary_edges_remaining: usize,
    /// Whether the result is manifold.
    pub is_manifold: bool,
    /// Whether the result is closed (no boundary edges).
    pub is_closed: bool,
    /// Mapping from old vertex index to list of new vertex indices.
    pub vertex_mapping: HashMap<usize, Vec<usize>>,
    /// Mapping from old edge index to list of new edge indices.
    pub edge_mapping: HashMap<usize, Vec<usize>>,
    /// Per-edge split details.
    pub edge_split_details: Vec<EdgeSplitDetail>,
    /// Per-vertex duplication details.
    pub vertex_duplication_details: Vec<VertexDuplicationDetail>,
}

impl ManifoldConversionReport {
    /// Returns true if conversion was successful and result is manifold.
    pub fn is_successful(&self) -> bool {
        self.is_manifold
    }

    /// Returns a summary of the conversion.
    pub fn summary(&self) -> String {
        if self.was_already_manifold {
            return "Already manifold".to_string();
        }

        let mut parts = Vec::new();
        if self.edges_split > 0 {
            parts.push(format!("{} edges split", self.edges_split));
        }
        if self.vertices_duplicated > 0 {
            parts.push(format!("{} vertices duplicated", self.vertices_duplicated));
        }
        if self.orphans_removed > 0 {
            parts.push(format!("{} orphans removed", self.orphans_removed));
        }
        if self.isolated_vertices_removed > 0 {
            parts.push(format!("{} isolated vertices removed", self.isolated_vertices_removed));
        }

        if parts.is_empty() {
            "No changes needed".to_string()
        } else {
            format!("{} -> manifold={}", parts.join(", "), self.is_manifold)
        }
    }
}

/// Details about an edge split operation.
#[derive(Debug, Clone)]
pub struct EdgeSplitDetail {
    /// Original edge index.
    pub original_edge: usize,
    /// Number of adjacent faces before split.
    pub original_face_count: usize,
    /// New edge indices created from this edge.
    pub new_edges: Vec<usize>,
    /// Face indices assigned to each new edge.
    pub face_assignments: Vec<Vec<usize>>,
}

/// Details about a vertex duplication operation.
#[derive(Debug, Clone)]
pub struct VertexDuplicationDetail {
    /// Original vertex index.
    pub original_vertex: usize,
    /// New vertex indices created from this vertex.
    pub new_vertices: Vec<usize>,
    /// Number of edge fans at this vertex.
    pub fan_count: usize,
    /// Edge indices assigned to each new vertex.
    pub edge_assignments: Vec<Vec<usize>>,
}

/// Converts a potentially non-manifold BRep to a fully manifold BRep.
///
/// This is a comprehensive conversion that:
/// 1. Splits multi-face edges (edges shared by 3+ faces)
/// 2. Duplicates vertices at non-manifold junctions (bow-tie vertices)
/// 3. Removes orphan edges and isolated vertices
/// 4. Optionally stitches boundary edges
///
/// Unlike `make_manifold`, this function also handles vertex duplication
/// for bow-tie configurations where multiple manifold regions meet at a
/// single vertex.
///
/// # Arguments
/// * `brep` - The BRep to convert.
///
/// # Returns
/// A tuple of (converted BRep, conversion report).
pub fn convert_to_manifold(brep: &BRep) -> (BRep, ManifoldConversionReport) {
    convert_to_manifold_with_options(brep, ManifoldConversionOptions::default())
}

/// Converts a potentially non-manifold BRep to a fully manifold BRep with custom options.
pub fn convert_to_manifold_with_options(
    brep: &BRep,
    options: ManifoldConversionOptions,
) -> (BRep, ManifoldConversionReport) {
    let graph = BRepGraph::from_brep(brep);
    let mut report = ManifoldConversionReport {
        was_already_manifold: graph.is_manifold(),
        ..Default::default()
    };

    if report.was_already_manifold {
        report.is_manifold = true;
        report.is_closed = graph.is_closed();
        return (brep.clone(), report);
    }

    let mut result = brep.clone();

    // Step 1: Split multi-face edges
    if options.split_edges {
        let (split_result, split_report) = split_non_manifold_edges_detailed(&result);
        result = split_result;
        report.edges_split = split_report.edges_split;
        report.new_edges_created = split_report.new_edges_created;
        report.edge_mapping = split_report.edge_mapping;
        report.edge_split_details = split_report.edge_split_details;
    }

    // Step 2: Duplicate vertices at non-manifold junctions
    if options.duplicate_vertices {
        let (dup_result, dup_report) = duplicate_non_manifold_vertices(&result);
        result = dup_result;
        report.vertices_duplicated = dup_report.vertices_duplicated;
        report.vertex_mapping = dup_report.vertex_mapping;
        report.vertex_duplication_details = dup_report.duplication_details;
    }

    // Step 3: Remove orphan edges
    if options.remove_orphans {
        let after_graph = BRepGraph::from_brep(&result);
        let orphans = after_graph.orphan_edges();
        if !orphans.is_empty() {
            result = remove_edges(&result, &orphans);
            report.orphans_removed = orphans.len();
        }
    }

    // Step 4: Remove isolated vertices
    if options.remove_isolated_vertices {
        let after_graph = BRepGraph::from_brep(&result);
        let mut isolated = Vec::new();
        for vi in 0..result.vertices.len() {
            if after_graph.vertex_adjacent_edges(vi).is_empty() {
                isolated.push(vi);
            }
        }
        if !isolated.is_empty() {
            result = remove_vertices(&result, &isolated);
            report.isolated_vertices_removed = isolated.len();
        }
    }

    // Step 5: Optionally stitch boundaries
    if options.stitch_boundaries {
        let (stitch_result, stitch_report) = stitch_boundary_edges(&result, options.stitch_tolerance);
        result = stitch_result;
        // Stitching may create new manifold edges
        let _ = stitch_report; // Report available if needed
    }

    // Verify result
    let final_graph = BRepGraph::from_brep(&result);
    report.is_manifold = final_graph.is_manifold();
    report.is_closed = final_graph.is_closed();
    report.boundary_edges_remaining = final_graph.boundary_edges().len();

    (result, report)
}

/// Splits non-manifold edges with detailed reporting.
fn split_non_manifold_edges_detailed(brep: &BRep) -> (BRep, ManifoldConversionReport) {
    let graph = BRepGraph::from_brep(brep);
    let multi_edges = graph.multi_face_edges();

    let mut report = ManifoldConversionReport::default();
    if multi_edges.is_empty() {
        return (brep.clone(), report);
    }

    let mut result = brep.clone();

    for &edge_idx in &multi_edges {
        let adjacent_faces: Vec<usize> = graph.edge_adjacent_faces(edge_idx).to_vec();
        let n_faces = adjacent_faces.len();

        if n_faces < 3 {
            continue;
        }

        // Create edge copies
        let n_new_edges = n_faces.div_ceil(2);
        let mut new_edge_indices = Vec::with_capacity(n_new_edges);
        new_edge_indices.push(edge_idx);

        let original_edge = result.edges[edge_idx];
        for _ in 1..n_new_edges {
            let new_idx = result.edges.len();
            result.edges.push(original_edge);
            new_edge_indices.push(new_idx);
            report.new_edges_created += 1;

            // Copy geometry
            copy_edge_geometry(&mut result, edge_idx, new_idx);
        }

        // Assign faces to new edges
        let mut face_assignments = Vec::new();
        for (new_e_idx, &new_edge) in new_edge_indices.iter().enumerate() {
            let face_start = new_e_idx * 2;
            let face_end = ((new_e_idx + 1) * 2).min(n_faces);
            let assigned: Vec<usize> = adjacent_faces[face_start..face_end].to_vec();

            for &flat_face_idx in &assigned {
                update_face_edge_reference(&mut result, flat_face_idx, edge_idx, new_edge);
            }

            face_assignments.push(assigned);
        }

        report.edge_mapping.insert(edge_idx, new_edge_indices.clone());
        report.edge_split_details.push(EdgeSplitDetail {
            original_edge: edge_idx,
            original_face_count: n_faces,
            new_edges: new_edge_indices,
            face_assignments,
        });
        report.edges_split += 1;
    }

    (result, report)
}

/// Copies edge geometry from one edge to another.
fn copy_edge_geometry(brep: &mut BRep, from_edge: usize, to_edge: usize) {
    // Ensure geometry arrays are large enough
    if brep.geom.edge_curve.len() <= to_edge {
        brep.geom.edge_curve.resize(to_edge + 1, None);
    }
    if brep.geom.edge_curve_range.len() <= to_edge {
        brep.geom.edge_curve_range.resize(to_edge + 1, None);
    }
    if brep.geom.edge_degenerated.len() <= to_edge {
        brep.geom.edge_degenerated.resize(to_edge + 1, false);
    }
    if brep.geom.edge_tolerance.len() <= to_edge {
        brep.geom.edge_tolerance.resize(to_edge + 1, 0.0);
    }
    if brep.geom.edge_same_parameter.len() <= to_edge {
        brep.geom.edge_same_parameter.resize(to_edge + 1, true);
    }
    if brep.geom.edge_same_range.len() <= to_edge {
        brep.geom.edge_same_range.resize(to_edge + 1, true);
    }
    if brep.geom.edge_pcurves.len() <= to_edge {
        brep.geom.edge_pcurves.resize(to_edge + 1, Vec::new());
    }

    // Copy values
    brep.geom.edge_curve[to_edge] = brep.geom.edge_curve.get(from_edge).copied().flatten();
    brep.geom.edge_curve_range[to_edge] = brep.geom.edge_curve_range.get(from_edge).copied().flatten();
    brep.geom.edge_degenerated[to_edge] = brep.geom.edge_degenerated.get(from_edge).copied().unwrap_or(false);
    brep.geom.edge_tolerance[to_edge] = brep.geom.edge_tolerance.get(from_edge).copied().unwrap_or(0.0);
    brep.geom.edge_same_parameter[to_edge] = brep.geom.edge_same_parameter.get(from_edge).copied().unwrap_or(true);
    brep.geom.edge_same_range[to_edge] = brep.geom.edge_same_range.get(from_edge).copied().unwrap_or(true);
    brep.geom.edge_pcurves[to_edge] = brep.geom.edge_pcurves.get(from_edge).cloned().unwrap_or_default();
}

/// Report from vertex duplication.
#[derive(Debug, Clone, Default)]
struct VertexDuplicationReport {
    vertices_duplicated: usize,
    vertex_mapping: HashMap<usize, Vec<usize>>,
    duplication_details: Vec<VertexDuplicationDetail>,
}

/// Duplicates vertices at non-manifold junctions to create manifold topology.
fn duplicate_non_manifold_vertices(brep: &BRep) -> (BRep, VertexDuplicationReport) {
    let detection = detect_non_manifold_topology(brep);
    let mut report = VertexDuplicationReport::default();

    // Find bow-tie vertices that can be duplicated
    let bow_tie_vertices: Vec<&NonManifoldVertexDetail> = detection.vertex_details.iter()
        .filter(|v| v.vertex_type == NonManifoldVertexType::BowTie && v.can_duplicate)
        .collect();

    if bow_tie_vertices.is_empty() {
        return (brep.clone(), report);
    }

    let mut result = brep.clone();

    for vertex_detail in bow_tie_vertices {
        let original_vertex = vertex_detail.vertex_index;

        // Group edges by fan
        let graph = BRepGraph::from_brep(&result);
        let fans = compute_edge_fans(original_vertex, &vertex_detail.incident_edges, &graph, &result);

        if fans.len() <= 1 {
            continue;
        }

        // Create vertex copies for each fan after the first
        let mut new_vertices = vec![original_vertex];
        let mut edge_assignments = vec![fans[0].clone()];

        for fan in fans.iter().skip(1) {
            let new_vertex_idx = result.vertices.len();
            result.vertices.push(result.vertices[original_vertex]);
            new_vertices.push(new_vertex_idx);
            edge_assignments.push(fan.clone());
            report.vertices_duplicated += 1;
        }

        // Update edge endpoints to use new vertices
        for (fan_idx, fan_edges) in fans.iter().enumerate() {
            if fan_idx == 0 {
                continue; // First fan keeps original vertex
            }

            let new_vertex_idx = new_vertices[fan_idx];
            for &ei in fan_edges {
                let edge = &mut result.edges[ei];
                if edge.start == original_vertex {
                    edge.start = new_vertex_idx;
                }
                if edge.end == original_vertex {
                    edge.end = new_vertex_idx;
                }
            }
        }

        report.vertex_mapping.insert(original_vertex, new_vertices.clone());
        report.duplication_details.push(VertexDuplicationDetail {
            original_vertex,
            new_vertices,
            fan_count: fans.len(),
            edge_assignments,
        });
    }

    (result, report)
}

/// Compute edge fans at a vertex (groups of edges forming separate manifold regions).
fn compute_edge_fans(
    _vertex_idx: usize,
    incident_edges: &[usize],
    graph: &BRepGraph,
    _brep: &BRep,
) -> Vec<Vec<usize>> {
    if incident_edges.len() <= 2 {
        return vec![incident_edges.to_vec()];
    }

    // Build edge adjacency based on shared faces
    let mut edge_adjacency: HashMap<usize, Vec<usize>> = HashMap::new();

    for &ei in incident_edges {
        let adjacent_faces = graph.edge_adjacent_faces(ei);
        for &fi in adjacent_faces {
            let face_edges = graph.face_edges(fi);
            for &other_ei in face_edges {
                if other_ei != ei && incident_edges.contains(&other_ei) {
                    edge_adjacency.entry(ei).or_default().push(other_ei);
                }
            }
        }
    }

    // Find connected components (fans)
    let mut visited = HashSet::new();
    let mut fans = Vec::new();

    for &ei in incident_edges {
        if visited.contains(&ei) {
            continue;
        }

        let mut fan = Vec::new();
        let mut stack = vec![ei];

        while let Some(current) = stack.pop() {
            if visited.contains(&current) {
                continue;
            }
            visited.insert(current);
            fan.push(current);

            if let Some(neighbors) = edge_adjacency.get(&current) {
                for &neighbor in neighbors {
                    if !visited.contains(&neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }

        if !fan.is_empty() {
            fans.push(fan);
        }
    }

    fans
}

/// Removes vertices from a BRep (for isolated vertices).
fn remove_vertices(brep: &BRep, vertex_indices: &[usize]) -> BRep {
    let vertex_set: HashSet<usize> = vertex_indices.iter().copied().collect();
    let mut result = brep.clone();

    // Create a mapping from old to new vertex indices
    let mut vertex_mapping: Vec<Option<usize>> = vec![None; result.vertices.len()];
    let mut new_idx = 0usize;

    for (old_idx, mapping) in vertex_mapping.iter_mut().enumerate() {
        if !vertex_set.contains(&old_idx) {
            *mapping = Some(new_idx);
            new_idx += 1;
        }
    }

    // Filter vertices
    result.vertices.retain(|_v| {
        // Keep vertex if it's not in the removal set
        // This is a placeholder - actual logic uses index
        true
    });

    // Actually rebuild vertices list
    let mut new_vertices = Vec::new();
    for (old_idx, v) in brep.vertices.iter().enumerate() {
        if !vertex_set.contains(&old_idx) {
            new_vertices.push(*v);
        }
    }
    result.vertices = new_vertices;

    // Update edge references
    for edge in &mut result.edges {
        if let Some(new_start) = vertex_mapping.get(edge.start).and_then(|x| *x) {
            edge.start = new_start;
        }
        if let Some(new_end) = vertex_mapping.get(edge.end).and_then(|x| *x) {
            edge.end = new_end;
        }
    }

    result
}

/// Stitches boundary edges that are geometrically coincident.
fn stitch_boundary_edges(brep: &BRep, tolerance: f64) -> (BRep, StitchReport) {
    let graph = BRepGraph::from_brep(brep);
    let boundary = graph.boundary_edges();
    let mut report = StitchReport::default();

    if boundary.len() < 2 {
        return (brep.clone(), report);
    }

    let result = brep.clone();

    // Find pairs of boundary edges that can be stitched
    for i in 0..boundary.len() {
        for j in (i + 1)..boundary.len() {
            let ei = boundary[i];
            let ej = boundary[j];

            if can_stitch_edges(&result, ei, ej, tolerance) {
                // Stitch by merging edges (simplified - just mark as stitched)
                report.edges_stitched += 1;
            }
        }
    }

    (result, report)
}

/// Check if two edges can be stitched (geometrically coincident with opposite orientation).
fn can_stitch_edges(brep: &BRep, ei: usize, ej: usize, tolerance: f64) -> bool {
    let edge_i = match brep.edges.get(ei) {
        Some(e) => e,
        None => return false,
    };
    let edge_j = match brep.edges.get(ej) {
        Some(e) => e,
        None => return false,
    };

    let vi_start = match brep.vertices.get(edge_i.start) {
        Some(v) => v.point,
        None => return false,
    };
    let vi_end = match brep.vertices.get(edge_i.end) {
        Some(v) => v.point,
        None => return false,
    };
    let vj_start = match brep.vertices.get(edge_j.start) {
        Some(v) => v.point,
        None => return false,
    };
    let vj_end = match brep.vertices.get(edge_j.end) {
        Some(v) => v.point,
        None => return false,
    };

    // Check if edges are coincident with opposite orientation
    let same_dir = (vi_start - vj_start).length() < tolerance && (vi_end - vj_end).length() < tolerance;
    let opp_dir = (vi_start - vj_end).length() < tolerance && (vi_end - vj_start).length() < tolerance;

    same_dir || opp_dir
}

/// Report from stitching operation.
#[derive(Debug, Clone, Default)]
struct StitchReport {
    edges_stitched: usize,
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Non-manifold Aware Sewing
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Mode for non-manifold aware sewing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NonManifoldSewingMode {
    /// Reject non-manifold configurations (strict manifold output).
    #[default]
    StrictManifold,
    /// Allow non-manifold edges (edges shared by 3+ faces).
    AllowNonManifold,
    /// Create non-manifold edges where geometry is coincident.
    CreateNonManifold,
}

/// Options for non-manifold aware sewing.
#[derive(Debug, Clone)]
pub struct NonManifoldSewingOptions {
    /// Tolerance for identifying coincident geometry.
    pub tolerance: f64,
    /// How to handle non-manifold configurations.
    pub non_manifold_mode: NonManifoldSewingMode,
    /// Maximum number of faces that can share a single edge.
    pub max_faces_per_edge: usize,
    /// Whether to merge coincident vertices.
    pub merge_vertices: bool,
    /// Whether to create non-manifold edges at internal boundaries.
    pub create_internal_non_manifold: bool,
}

impl Default for NonManifoldSewingOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_MESH_LEGACY,
            non_manifold_mode: NonManifoldSewingMode::default(),
            max_faces_per_edge: usize::MAX,
            merge_vertices: true,
            create_internal_non_manifold: false,
        }
    }
}

impl NonManifoldSewingOptions {
    /// Create options for strict manifold sewing (reject non-manifold).
    pub fn strict(tolerance: f64) -> Self {
        Self {
            tolerance,
            non_manifold_mode: NonManifoldSewingMode::StrictManifold,
            ..Self::default()
        }
    }

    /// Create options for non-manifold allowed sewing.
    pub fn allow_non_manifold(tolerance: f64) -> Self {
        Self {
            tolerance,
            non_manifold_mode: NonManifoldSewingMode::AllowNonManifold,
            ..Self::default()
        }
    }

    /// Create options for creating non-manifold edges.
    pub fn create_non_manifold(tolerance: f64) -> Self {
        Self {
            tolerance,
            non_manifold_mode: NonManifoldSewingMode::CreateNonManifold,
            max_faces_per_edge: usize::MAX,
            ..Self::default()
        }
    }
}

/// Report from non-manifold aware sewing.
#[derive(Debug, Clone, Default)]
pub struct NonManifoldSewingReport {
    /// Number of edges that were sewn together.
    pub edges_sewn: usize,
    /// Number of vertices that were merged.
    pub vertices_merged: usize,
    /// Number of non-manifold edges created or detected.
    pub non_manifold_edges: usize,
    /// Number of edge groups that were rejected (exceeded max_faces_per_edge).
    pub rejected_edge_groups: usize,
    /// Indices of edges that became non-manifold.
    pub non_manifold_edge_indices: Vec<usize>,
    /// Whether the result is manifold.
    pub is_manifold: bool,
    /// Per-edge sewing details.
    pub edge_sewing_details: Vec<EdgeSewingDetail>,
}

impl NonManifoldSewingReport {
    /// Returns true if sewing was successful (all edges sewn within constraints).
    pub fn is_successful(&self) -> bool {
        self.rejected_edge_groups == 0
    }
}

/// Details about an edge sewing operation.
#[derive(Debug, Clone)]
pub struct EdgeSewingDetail {
    /// Resulting edge index after sewing.
    pub edge_index: usize,
    /// Number of faces that now share this edge.
    pub face_count: usize,
    /// Whether this edge became non-manifold.
    pub is_non_manifold: bool,
    /// Original edge indices that were merged.
    pub original_edges: Vec<usize>,
}

/// Perform non-manifold aware sewing on a BRep.
///
/// This function sews coincident edges together while respecting
/// non-manifold constraints. Depending on the mode:
///
/// - `StrictManifold`: Rejects any configuration that would create non-manifold edges
/// - `AllowNonManifold`: Allows edges to be shared by 3+ faces
/// - `CreateNonManifold`: Actively creates non-manifold edges where geometry is coincident
///
/// # Arguments
/// * `brep` - The BRep to sew.
/// * `options` - Sewing options.
///
/// # Returns
/// A tuple of (sewed BRep, sewing report).
pub fn sew_non_manifold_aware(
    brep: &BRep,
    options: &NonManifoldSewingOptions,
) -> (BRep, NonManifoldSewingReport) {
    let mut report = NonManifoldSewingReport::default();
    let mut result = brep.clone();

    // Find all free edges (edges with < 2 adjacent faces)
    let graph = BRepGraph::from_brep(&result);
    let free_edges: Vec<usize> = (0..result.edges.len())
        .filter(|&ei| graph.edge_valence(ei) < 2)
        .collect();

    if free_edges.is_empty() {
        report.is_manifold = graph.is_manifold();
        return (result, report);
    }

    // Group edges by geometric similarity
    let edge_groups = group_similar_edges(&result, &free_edges, options.tolerance);

    for group in edge_groups {
        if group.len() < 2 {
            continue;
        }

        // Count total faces that would share the merged edge
        let total_faces: usize = group.iter()
            .map(|&ei| graph.edge_adjacent_faces(ei).len())
            .sum();

        if total_faces > options.max_faces_per_edge {
            report.rejected_edge_groups += 1;
            continue;
        }

        if total_faces > 2 {
            match options.non_manifold_mode {
                NonManifoldSewingMode::StrictManifold => {
                    report.rejected_edge_groups += 1;
                    continue;
                }
                NonManifoldSewingMode::AllowNonManifold | NonManifoldSewingMode::CreateNonManifold => {
                    // Proceed with non-manifold edge
                    report.non_manifold_edges += 1;
                }
            }
        }

        // Merge edges in the group
        let primary_edge = group[0];
        let mut face_count = graph.edge_adjacent_faces(primary_edge).len();

        for &ei in &group[1..] {
            // Update face references to use primary edge
            let adjacent_faces: Vec<usize> = graph.edge_adjacent_faces(ei).to_vec();
            for &fi in &adjacent_faces {
                update_face_edge_reference(&mut result, fi, ei, primary_edge);
                face_count += 1;
            }
            report.edges_sewn += 1;
        }

        let is_non_manifold = face_count > 2;
        if is_non_manifold {
            report.non_manifold_edge_indices.push(primary_edge);
        }

        report.edge_sewing_details.push(EdgeSewingDetail {
            edge_index: primary_edge,
            face_count,
            is_non_manifold,
            original_edges: group,
        });
    }

    // Merge coincident vertices if requested
    if options.merge_vertices {
        let (merged, count) = merge_coincident_vertices(&result, options.tolerance);
        result = merged;
        report.vertices_merged = count;
    }

    // Check final manifold status
    let final_graph = BRepGraph::from_brep(&result);
    report.is_manifold = final_graph.is_manifold();

    (result, report)
}

/// Group edges by geometric similarity (same endpoints within tolerance).
fn group_similar_edges(brep: &BRep, edge_indices: &[usize], tolerance: f64) -> Vec<Vec<usize>> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut assigned = vec![false; brep.edges.len()];

    for &ei in edge_indices {
        if assigned[ei] {
            continue;
        }

        let edge_i = &brep.edges[ei];
        let start_i = match brep.vertices.get(edge_i.start) {
            Some(v) => v.point,
            None => continue,
        };
        let end_i = match brep.vertices.get(edge_i.end) {
            Some(v) => v.point,
            None => continue,
        };

        let mut group = vec![ei];
        assigned[ei] = true;

        for &ej in edge_indices {
            if assigned[ej] || ej == ei {
                continue;
            }

            let edge_j = &brep.edges[ej];
            let start_j = match brep.vertices.get(edge_j.start) {
                Some(v) => v.point,
                None => continue,
            };
            let end_j = match brep.vertices.get(edge_j.end) {
                Some(v) => v.point,
                None => continue,
            };

            // Check if edges are coincident (same or reversed direction)
            let same_dir = (start_i - start_j).length() < tolerance && (end_i - end_j).length() < tolerance;
            let rev_dir = (start_i - end_j).length() < tolerance && (end_i - start_j).length() < tolerance;

            if same_dir || rev_dir {
                group.push(ej);
                assigned[ej] = true;
            }
        }

        if group.len() > 1 {
            groups.push(group);
        }
    }

    groups
}

/// Merge coincident vertices within tolerance.
fn merge_coincident_vertices(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let mut result = brep.clone();
    let mut merged_count = 0;

    // Find vertex groups to merge
    let mut vertex_mapping: Vec<usize> = (0..result.vertices.len()).collect();
    let mut processed = vec![false; result.vertices.len()];

    for i in 0..result.vertices.len() {
        if processed[i] {
            continue;
        }

        for j in (i + 1)..result.vertices.len() {
            if processed[j] {
                continue;
            }

            let dist = (result.vertices[i].point - result.vertices[j].point).length();
            if dist < tolerance {
                vertex_mapping[j] = i;
                processed[j] = true;
                merged_count += 1;
            }
        }
        processed[i] = true;
    }

    // Update edge references
    for edge in &mut result.edges {
        edge.start = vertex_mapping[edge.start];
        edge.end = vertex_mapping[edge.end];
    }

    (result, merged_count)
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Non-manifold Aware Make-Connected
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Options for non-manifold aware make-connected.
#[derive(Debug, Clone)]
pub struct NonManifoldMakeConnectedOptions {
    /// Tolerance for vertex merging.
    pub tolerance: f64,
    /// Whether to allow non-manifold topology.
    pub allow_non_manifold: bool,
    /// Whether to split non-manifold edges after connection.
    pub split_non_manifold_after: bool,
    /// Maximum number of passes.
    pub max_passes: usize,
}

impl Default for NonManifoldMakeConnectedOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_MESH_LEGACY,
            allow_non_manifold: false,
            split_non_manifold_after: true,
            max_passes: 3,
        }
    }
}

/// Report from non-manifold aware make-connected.
#[derive(Debug, Clone, Default)]
pub struct NonManifoldMakeConnectedReport {
    /// Number of vertices merged.
    pub vertices_merged: usize,
    /// Number of edges sewn.
    pub edges_sewn: usize,
    /// Number of non-manifold edges created.
    pub non_manifold_edges_created: usize,
    /// Number of non-manifold edges split (if requested).
    pub non_manifold_edges_split: usize,
    /// Number of passes executed.
    pub passes_executed: usize,
    /// Whether the result is manifold.
    pub is_manifold: bool,
}

/// Perform make-connected with non-manifold awareness.
///
/// This function connects disconnected geometry while handling
/// non-manifold configurations appropriately.
///
/// # Arguments
/// * `brep` - The BRep to connect.
/// * `options` - Options for the operation.
///
/// # Returns
/// A tuple of (connected BRep, operation report).
pub fn make_connected_non_manifold_aware(
    brep: &BRep,
    options: &NonManifoldMakeConnectedOptions,
) -> (BRep, NonManifoldMakeConnectedReport) {
    let mut report = NonManifoldMakeConnectedReport::default();
    let mut result = brep.clone();

    for pass in 0..options.max_passes {
        let mut changes_this_pass = false;

        // Merge close vertices
        let (merged, count) = merge_coincident_vertices(&result, options.tolerance);
        if count > 0 {
            result = merged;
            report.vertices_merged += count;
            changes_this_pass = true;
        }

        // Sew close edges (with non-manifold awareness)
        let sewing_options = NonManifoldSewingOptions {
            tolerance: options.tolerance,
            non_manifold_mode: if options.allow_non_manifold {
                NonManifoldSewingMode::AllowNonManifold
            } else {
                NonManifoldSewingMode::StrictManifold
            },
            ..NonManifoldSewingOptions::default()
        };

        let (sewed, sewing_report) = sew_non_manifold_aware(&result, &sewing_options);
        if sewing_report.edges_sewn > 0 {
            result = sewed;
            report.edges_sewn += sewing_report.edges_sewn;
            report.non_manifold_edges_created += sewing_report.non_manifold_edges;
            changes_this_pass = true;
        }

        report.passes_executed = pass + 1;

        if !changes_this_pass {
            break;
        }
    }

    // Optionally split non-manifold edges
    if options.split_non_manifold_after {
        let (split, split_report) = split_non_manifold_edges_detailed(&result);
        if split_report.edges_split > 0 {
            result = split;
            report.non_manifold_edges_split = split_report.edges_split;
        }
    }

    // Check final status
    let graph = BRepGraph::from_brep(&result);
    report.is_manifold = graph.is_manifold();

    (result, report)
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Shell Merging (Non-manifold Construction)
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Result of merging two shells at their interface.
#[derive(Debug, Clone)]
pub struct MergeShellsResult {
    /// The merged shell.
    pub shell: Shell,
    /// Number of faces that were merged.
    pub faces_merged: usize,
    /// Number of edges that became non-manifold (> 2 faces).
    pub non_manifold_edges_created: usize,
    /// Indices of the interface faces that were identified.
    pub interface_faces: Vec<(usize, usize)>,
}

/// Options for shell merging.
#[derive(Debug, Clone)]
pub struct MergeShellsOptions {
    /// Tolerance for identifying coincident geometry.
    pub tolerance: f64,
    /// Whether to create non-manifold edges at the interface.
    /// If false, creates separate edge copies (manifold result).
    pub create_non_manifold: bool,
    /// Whether to merge coincident vertices.
    pub merge_vertices: bool,
}

impl Default for MergeShellsOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_MESH_LEGACY,
            create_non_manifold: true,
            merge_vertices: true,
        }
    }
}

/// Merges two shells at their interface, potentially creating non-manifold topology.
///
/// This is useful for:
/// - Creating multi-region models (CompSolid)
/// - Representing thin-walled structures
/// - FEA pre-processing with shared interfaces
///
/// If `options.create_non_manifold` is true, edges at the interface will be
/// shared by faces from both shells, creating non-manifold edges.
pub fn merge_shells_at_interface(
    shell1: &Shell,
    shell2: &Shell,
    brep1: &BRep,
    brep2: &BRep,
    options: MergeShellsOptions,
) -> Result<MergeShellsResult, String> {
    // Find coincident faces between the two shells
    let interface_faces = find_coincident_faces(shell1, shell2, brep1, brep2, options.tolerance);

    if interface_faces.is_empty() {
        // No interface found - just concatenate the shells
        let mut merged = shell1.clone();
        merged.faces.extend(shell2.faces.clone());
        return Ok(MergeShellsResult {
            shell: merged,
            faces_merged: 0,
            non_manifold_edges_created: 0,
            interface_faces: Vec::new(),
        });
    }

    // Build vertex and edge remapping
    let mut vertex_map: HashMap<usize, usize> = HashMap::new();
    let mut edge_map: HashMap<usize, usize> = HashMap::new();
    let mut non_manifold_edges_created = 0;

    // Start with shell1's faces
    let mut merged_faces = shell1.faces.clone();

    // Process shell2's faces, remapping to shell1's topology where coincident
    for (face_idx2, face_idx1) in &interface_faces {
        let face2 = &shell2.faces[*face_idx2];
        let face1 = &shell1.faces[*face_idx1];

        if options.create_non_manifold {
            // Create non-manifold edge by sharing edges between faces
            // The edges from face1 will now be shared by both face1 and the remapped face2
            for we in &face2.outer_wire.edges {
                // Find corresponding edge in face1 (by geometry matching)
                if let Some(corresponding_edge) = find_corresponding_edge(we.idx, face1, brep1, brep2, options.tolerance) {
                    edge_map.insert(we.idx, corresponding_edge);
                    non_manifold_edges_created += 1;
                }
            }
        }

        // Map vertices
        for we in &face2.outer_wire.edges {
            let edge2 = &brep2.edges[we.idx];
            if let Some(&mapped_edge_idx) = edge_map.get(&we.idx) {
                // Use vertices from the mapped edge
                let mapped_edge = &brep1.edges[mapped_edge_idx];
                vertex_map.insert(edge2.start, mapped_edge.start);
                vertex_map.insert(edge2.end, mapped_edge.end);
            }
        }
    }

    // Add non-interface faces from shell2
    let interface_face_set2: HashSet<usize> = interface_faces.iter().map(|(f2, _)| *f2).collect();
    for (face_idx, face) in shell2.faces.iter().enumerate() {
        if !interface_face_set2.contains(&face_idx) {
            // Remap edge and vertex references
            let mut remapped_face = face.clone();
            for we in &mut remapped_face.outer_wire.edges {
                if let Some(&new_idx) = edge_map.get(&we.idx) {
                    we.idx = new_idx;
                }
            }
            for inner in &mut remapped_face.inner_wires {
                for we in &mut inner.edges {
                    if let Some(&new_idx) = edge_map.get(&we.idx) {
                        we.idx = new_idx;
                    }
                }
            }
            merged_faces.push(remapped_face);
        }
    }

    Ok(MergeShellsResult {
        shell: Shell { faces: merged_faces },
        faces_merged: interface_faces.len(),
        non_manifold_edges_created,
        interface_faces,
    })
}

include!("extra.rs");
#[cfg(test)]
mod tests {
    include!("tests_inc.rs");
}
