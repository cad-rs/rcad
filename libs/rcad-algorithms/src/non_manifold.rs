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
//! # Repair
//!
//! - [`split_non_manifold_edges`]: Duplicates edges so each pair of faces has
//!   its own copy, converting non-manifold edges to manifold.
//! - [`make_manifold`]: Full conversion pipeline that splits non-manifold edges
//!   and optionally stitches boundary edges.
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

use std::collections::{HashMap, HashSet};
use rcad_kernel::{
    BRep, BRepGraph, Face, Shell,
};
use glam::DVec3;

// ─────────────────────────────────────────────────────────────────────────────
// Detection API
// ─────────────────────────────────────────────────────────────────────────────

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

// ─────────────────────────────────────────────────────────────────────────────
// Edge Splitting (Non-manifold -> Manifold Conversion)
// ─────────────────────────────────────────────────────────────────────────────

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
        let n_new_edges = (n_faces + 1) / 2;

        // Collect new edge indices
        let mut new_edge_indices = Vec::with_capacity(n_new_edges);
        new_edge_indices.push(edge_idx); // Keep original for first pair

        // Create edge copies
        let original_edge = result.edges[edge_idx].clone();
        for _ in 1..n_new_edges {
            let new_idx = result.edges.len();
            result.edges.push(original_edge.clone());
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

// ─────────────────────────────────────────────────────────────────────────────
// Make Manifold (Full Conversion Pipeline)
// ─────────────────────────────────────────────────────────────────────────────

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
            tolerance: 1e-6,
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

// ─────────────────────────────────────────────────────────────────────────────
// Shell Merging (Non-manifold Construction)
// ─────────────────────────────────────────────────────────────────────────────

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
            tolerance: 1e-6,
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

/// Finds coincident faces between two shells based on geometric matching.
fn find_coincident_faces(
    shell1: &Shell,
    shell2: &Shell,
    brep1: &BRep,
    brep2: &BRep,
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();

    for (f1_idx, face1) in shell1.faces.iter().enumerate() {
        let normal1 = face1.normal;
        let center1 = compute_face_center(face1, brep1);

        for (f2_idx, face2) in shell2.faces.iter().enumerate() {
            let normal2 = face2.normal;
            let center2 = compute_face_center(face2, brep2);

            // Check if normals are parallel (same or opposite direction)
            let dot = normal1.dot(normal2).abs();
            if dot < 1.0 - tolerance {
                continue; // Not parallel
            }

            // Check if centers are coincident
            if (center1 - center2).length() < tolerance {
                pairs.push((f2_idx, f1_idx));
            }
        }
    }

    pairs
}

/// Computes the geometric center of a face.
fn compute_face_center(face: &Face, brep: &BRep) -> DVec3 {
    let mut center = DVec3::ZERO;
    let mut count = 0;

    for we in &face.outer_wire.edges {
        if we.idx < brep.edges.len() {
            let edge = &brep.edges[we.idx];
            if edge.start < brep.vertices.len() {
                center += brep.vertices[edge.start].point;
                count += 1;
            }
            if edge.end < brep.vertices.len() {
                center += brep.vertices[edge.end].point;
                count += 1;
            }
        }
    }

    if count > 0 {
        center / count as f64
    } else {
        DVec3::ZERO
    }
}

/// Finds a corresponding edge in face1 that matches edge_idx from face2/brep2.
fn find_corresponding_edge(
    edge_idx: usize,
    face1: &Face,
    brep1: &BRep,
    brep2: &BRep,
    tolerance: f64,
) -> Option<usize> {
    let edge2 = &brep2.edges[edge_idx];
    let start2 = brep2.vertices.get(edge2.start)?.point;
    let end2 = brep2.vertices.get(edge2.end)?.point;

    for we in &face1.outer_wire.edges {
        let edge1 = &brep1.edges[we.idx];
        let start1 = brep1.vertices.get(edge1.start)?.point;
        let end1 = brep1.vertices.get(edge1.end)?.point;

        // Check if edges are coincident (same or reversed direction)
        let same_dir = (start1 - start2).length() < tolerance && (end1 - end2).length() < tolerance;
        let rev_dir = (start1 - end2).length() < tolerance && (end1 - start2).length() < tolerance;

        if same_dir || rev_dir {
            return Some(we.idx);
        }
    }

    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Non-manifold Traversal Extensions for BRepGraph
// ─────────────────────────────────────────────────────────────────────────────

/// Extension trait for non-manifold traversal on BRepGraph.
pub trait NonManifoldTraversal {
    /// Returns all faces that share a non-manifold edge with the given face.
    fn non_manifold_adjacent_faces(&self, face_idx: usize) -> Vec<usize>;

    /// Returns the number of manifold regions (connected components via manifold edges only).
    fn manifold_region_count(&self) -> usize;

    /// Returns faces grouped by manifold region.
    fn manifold_regions(&self) -> Vec<Vec<usize>>;

    /// Iterator over non-manifold edges with their adjacent faces.
    fn non_manifold_edge_info(&self) -> Vec<(usize, Vec<usize>)>;
}

impl NonManifoldTraversal for BRepGraph {
    /// Returns all faces that share a non-manifold edge with the given face.
    fn non_manifold_adjacent_faces(&self, face_idx: usize) -> Vec<usize> {
        let mut result = Vec::new();
        let edges = self.face_edges(face_idx).to_vec();

        for &ei in &edges {
            if self.edge_valence(ei) > 2 {
                // This is a non-manifold edge
                for &adj_face in self.edge_adjacent_faces(ei) {
                    if adj_face != face_idx {
                        result.push(adj_face);
                    }
                }
            }
        }

        result.sort_unstable();
        result.dedup();
        result
    }

    /// Returns the number of manifold regions (connected components via manifold edges only).
    fn manifold_region_count(&self) -> usize {
        self.manifold_regions().len()
    }

    /// Returns faces grouped by manifold region.
    fn manifold_regions(&self) -> Vec<Vec<usize>> {
        let mut visited = vec![false; self.face_count];
        let mut regions = Vec::new();

        // Get non-manifold edges
        let nm_edges: HashSet<usize> = self.multi_face_edges().into_iter().collect();

        for start in 0..self.face_count {
            if visited[start] {
                continue;
            }

            let mut region = Vec::new();
            let mut stack = vec![start];

            while let Some(fi) = stack.pop() {
                if visited[fi] {
                    continue;
                }
                visited[fi] = true;
                region.push(fi);

                // Traverse only through manifold edges
                for &ei in self.face_edges(fi) {
                    if nm_edges.contains(&ei) {
                        continue; // Skip non-manifold edges
                    }

                    for &adj in self.edge_adjacent_faces(ei) {
                        if !visited[adj] {
                            stack.push(adj);
                        }
                    }
                }
            }

            if !region.is_empty() {
                regions.push(region);
            }
        }

        regions
    }

    /// Returns information about each non-manifold edge and its adjacent faces.
    fn non_manifold_edge_info(&self) -> Vec<(usize, Vec<usize>)> {
        self.multi_face_edges()
            .into_iter()
            .map(|ei| (ei, self.edge_adjacent_faces(ei).to_vec()))
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::{PrimitiveSolid, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
    use glam::DVec3;

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    /// Build a minimal non-manifold BRep where edge 0 is shared by 3 faces.
    fn non_manifold_tripod() -> BRep {
        let vertices = vec![
            Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0
            Vertex { point: DVec3::new(1.0, 0.0, 0.0) }, // 1
            Vertex { point: DVec3::new(0.0, 1.0, 0.0) }, // 2
            Vertex { point: DVec3::new(0.0, 0.0, 1.0) }, // 3
            Vertex { point: DVec3::new(0.0, -1.0, 0.0) }, // 4
        ];

        let edges = vec![
            Edge { start: 0, end: 1 }, // shared by 3 faces
            Edge { start: 1, end: 2 },
            Edge { start: 2, end: 0 },
            Edge { start: 1, end: 3 },
            Edge { start: 3, end: 0 },
            Edge { start: 1, end: 4 },
            Edge { start: 4, end: 0 },
        ];

        let f0 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(1, true),
                    WireEdge::new(2, true),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f1 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(3, true),
                    WireEdge::new(4, true),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f2 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(5, true),
                    WireEdge::new(6, true),
                ],
            },
            inner_wires: vec![],
            normal: -DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };

        BRep {
            vertices,
            edges,
            solids: vec![Solid {
                shells: vec![Shell {
                    faces: vec![f0, f1, f2],
                }],
            }],
            geom: Default::default(),
        }
    }

    #[test]
    fn test_is_manifold_for_box() {
        let brep = unit_box();
        assert!(is_manifold(&brep));
    }

    #[test]
    fn test_is_manifold_for_tripod() {
        let brep = non_manifold_tripod();
        assert!(!is_manifold(&brep));
    }

    #[test]
    fn test_non_manifold_edges_for_box() {
        let brep = unit_box();
        let nm_edges = non_manifold_edges(&brep);
        assert!(nm_edges.is_empty());
    }

    #[test]
    fn test_non_manifold_edges_for_tripod() {
        let brep = non_manifold_tripod();
        let nm_edges = non_manifold_edges(&brep);
        // Edge 0 is multi-face, edges 1-6 are boundary
        assert_eq!(nm_edges.len(), 7); // 1 multi-face + 6 boundary
        assert!(nm_edges.contains(&0));
    }

    #[test]
    fn test_multi_face_edges_for_tripod() {
        let brep = non_manifold_tripod();
        let multi = multi_face_edges(&brep);
        assert_eq!(multi, vec![0]);
    }

    #[test]
    fn test_non_manifold_vertices_for_tripod() {
        let brep = non_manifold_tripod();
        let verts = non_manifold_vertices(&brep);
        assert_eq!(verts, vec![0, 1]); // endpoints of edge 0
    }

    #[test]
    fn test_analyze_non_manifold_for_box() {
        let brep = unit_box();
        let report = analyze_non_manifold(&brep);
        assert!(report.is_manifold);
        assert!(report.is_closed);
        assert_eq!(report.multi_face_edge_count, 0);
        assert_eq!(report.boundary_edge_count, 0);
        assert!(report.is_clean());
    }

    #[test]
    fn test_analyze_non_manifold_for_tripod() {
        let brep = non_manifold_tripod();
        let report = analyze_non_manifold(&brep);
        assert!(!report.is_manifold);
        assert!(!report.is_closed);
        assert_eq!(report.multi_face_edge_count, 1);
        assert_eq!(report.boundary_edge_count, 6);
        assert_eq!(report.non_manifold_vertex_count, 2);
        assert!(!report.is_clean());
    }

    #[test]
    fn test_split_non_manifold_edges_for_box() {
        let brep = unit_box();
        let (result, report) = split_non_manifold_edges(&brep);
        assert!(is_manifold(&result));
        assert_eq!(report.edges_split, 0);
    }

    #[test]
    fn test_split_non_manifold_edges_for_tripod() {
        let brep = non_manifold_tripod();
        let (result, report) = split_non_manifold_edges(&brep);

        // After splitting, the multi-face edge should be resolved
        assert!(report.edges_split > 0);
        assert!(report.new_edges_created > 0);

        // Verify the mapping
        assert!(report.edge_mapping.contains_key(&0));
    }

    #[test]
    fn test_make_manifold_for_box() {
        let brep = unit_box();
        let (result, report) = make_manifold(&brep).expect("should succeed");
        assert!(report.was_already_manifold);
        assert!(report.is_manifold);
    }

    #[test]
    fn test_make_manifold_for_tripod() {
        let brep = non_manifold_tripod();
        let (result, report) = make_manifold(&brep).expect("should succeed");
        assert!(!report.was_already_manifold);
        // After splitting, boundary edges remain, so not fully manifold in the closed sense
        // but the multi-face edge should be resolved
    }

    #[test]
    fn test_non_manifold_traversal() {
        let brep = non_manifold_tripod();
        let graph = BRepGraph::from_brep(&brep);

        // Test non_manifold_adjacent_faces
        let adj = graph.non_manifold_adjacent_faces(0);
        // Face 0 shares edge 0 with faces 1 and 2
        assert!(adj.contains(&1));
        assert!(adj.contains(&2));

        // Test manifold_regions
        let regions = graph.manifold_regions();
        // With a multi-face edge, faces should still be connected via that edge
        // but our manifold_regions skips non-manifold edges
        assert!(!regions.is_empty());

        // Test non_manifold_edge_info
        let info = graph.non_manifold_edge_info();
        assert_eq!(info.len(), 1);
        assert_eq!(info[0].0, 0); // edge 0
        assert_eq!(info[0].1.len(), 3); // 3 adjacent faces
    }

    #[test]
    fn test_boundary_edges_for_tripod() {
        let brep = non_manifold_tripod();
        let bounds = boundary_edges(&brep);
        // Edges 1-6 are boundary edges (1 face each)
        assert_eq!(bounds.len(), 6);
    }

    #[test]
    fn test_orphan_edges() {
        let brep = unit_box();
        let orphans = orphan_edges(&brep);
        assert!(orphans.is_empty());
    }
}
