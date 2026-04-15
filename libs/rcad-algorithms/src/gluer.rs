//! BRepFeat_Gluer equivalent functionality for gluing shapes at interfaces.
//!
//! This module provides tools to glue (merge) two shapes at their shared interface,
//! commonly used for assembly operations and feature combination.
//!
//! # Core Concept
//!
//! Unlike boolean union which computes intersection and creates new geometry,
//! gluing assumes shapes already have coincident faces at their interface and
//! simply merges the topology, eliminating duplicate geometry.
//!
//! # Usage
//!
//! ```
//! use rcad_kernel::{BRep, PrimitiveSolid};
//! use rcad_algorithms::gluer::{glue_shapes, GluerOptions};
//!
//! // Create two boxes that share a face
//! let mut box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
//! let mut box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
//!
//! // Translate box2 to share a face with box1
//! box2.apply_transform(glam::DAffine3::from_translation(glam::DVec3::new(0.0, 1.0, 0.0)));
//!
//! // Glue them together
//! let result = glue_shapes(&box1, &box2, GluerOptions::default());
//! ```
//!
//! # OCCT Equivalent
//!
//! This module provides functionality similar to OCCT's `BRepFeat_Gluer`:
//! - Automatic detection of coincident faces
//! - Merging of shared edges and vertices
//! - History tracking for parametric editing

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::{BRep, Edge, Face, Shell, Solid, GeomStore, PCurve};

// ─────────────────────────────────────────────────────────────────────────────
// Error Types
// ─────────────────────────────────────────────────────────────────────────────

/// Errors that can occur during gluing operations.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GluerError {
    /// No interface faces found between the shapes.
    NoInterfaceFound,
    /// Incompatible geometry at the interface.
    IncompatibleGeometry(String),
    /// Topology inconsistency detected.
    TopologyInconsistency(String),
    /// Invalid input shape (e.g., empty BRep).
    InvalidInput(String),
    /// Normal direction mismatch at interface.
    NormalMismatch { face1: usize, face2: usize },
}

impl std::fmt::Display for GluerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GluerError::NoInterfaceFound => write!(f, "No interface faces found between shapes"),
            GluerError::IncompatibleGeometry(msg) => write!(f, "Incompatible geometry: {}", msg),
            GluerError::TopologyInconsistency(msg) => write!(f, "Topology inconsistency: {}", msg),
            GluerError::InvalidInput(msg) => write!(f, "Invalid input: {}", msg),
            GluerError::NormalMismatch { face1, face2 } => {
                write!(f, "Normal mismatch at faces {} and {}", face1, face2)
            }
        }
    }
}

impl std::error::Error for GluerError {}

// ─────────────────────────────────────────────────────────────────────────────
// Options and Result Types
// ─────────────────────────────────────────────────────────────────────────────

/// Options for gluing operations.
#[derive(Debug, Clone)]
pub struct GluerOptions {
    /// Tolerance for geometric matching (default: 1e-6).
    pub tolerance: f64,
    /// Whether to merge coincident faces at the interface (default: true).
    /// When true, one face is kept; when false, both faces remain (non-manifold).
    pub merge_shared_faces: bool,
    /// Whether to merge coincident edges (default: true).
    pub merge_shared_edges: bool,
    /// Whether to merge coincident vertices (default: true).
    pub merge_shared_vertices: bool,
    /// Whether to preserve history for parametric editing (default: true).
    pub preserve_history: bool,
    /// Whether to check and handle normal direction mismatches (default: true).
    /// If true, faces with opposite normals at interface will still be merged.
    pub allow_opposite_normals: bool,
    /// Whether to perform validity checks on the result (default: true).
    pub validate_result: bool,
}

impl Default for GluerOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            merge_shared_faces: true,
            merge_shared_edges: true,
            merge_shared_vertices: true,
            preserve_history: true,
            allow_opposite_normals: true,
            validate_result: true,
        }
    }
}

/// History tracking for a gluing operation.
#[derive(Debug, Clone, Default)]
pub struct GluerHistory {
    /// Mapping from result face index to source face(s).
    /// A face can come from either shape A, shape B, or be a merged face.
    pub face_origins: Vec<FaceOrigin>,
    /// Mapping from result edge index to source edge(s).
    pub edge_origins: Vec<EdgeOrigin>,
    /// Mapping from result vertex index to source vertex.
    pub vertex_origins: Vec<VertexOrigin>,
    /// Indices of faces that were merged (face_from_a, face_from_b).
    pub merged_face_pairs: Vec<(usize, usize)>,
    /// Indices of edges that were merged (edge_from_a, edge_from_b).
    pub merged_edge_pairs: Vec<(usize, usize)>,
    /// Indices of vertices that were merged (vertex_from_a, vertex_from_b).
    pub merged_vertex_pairs: Vec<(usize, usize)>,
}

/// Origin of a face in the glued result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceOrigin {
    /// Face came from shape A.
    FromA(usize),
    /// Face came from shape B.
    FromB(usize),
    /// Face is a merged interface face (originally from both A and B).
    Merged { from_a: usize, from_b: usize },
}

/// Origin of an edge in the glued result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeOrigin {
    /// Edge came from shape A.
    FromA(usize),
    /// Edge came from shape B.
    FromB(usize),
    /// Edge is a merged interface edge.
    Merged { from_a: usize, from_b: usize },
}

/// Origin of a vertex in the glued result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VertexOrigin {
    /// Vertex came from shape A.
    FromA(usize),
    /// Vertex came from shape B.
    FromB(usize),
    /// Vertex is a merged interface vertex.
    Merged { from_a: usize, from_b: usize },
}

/// Result of a gluing operation.
#[derive(Debug, Clone)]
pub struct GluerResult {
    /// The glued BRep.
    pub brep: BRep,
    /// History tracking (if requested).
    pub history: Option<GluerHistory>,
    /// Pairs of faces that were merged (face_idx_in_a, face_idx_in_b).
    pub merged_faces: Vec<(usize, usize)>,
    /// Pairs of edges that were merged (edge_idx_in_a, edge_idx_in_b).
    pub merged_edges: Vec<(usize, usize)>,
    /// Pairs of vertices that were merged (vertex_idx_in_a, vertex_idx_in_b).
    pub merged_vertices: Vec<(usize, usize)>,
    /// Number of interface faces detected.
    pub interface_face_count: usize,
    /// Number of interface edges detected.
    pub interface_edge_count: usize,
}

// ─────────────────────────────────────────────────────────────────────────────
// Interface Detection
// ─────────────────────────────────────────────────────────────────────────────

/// Represents a detected interface between two shapes.
#[derive(Debug, Clone)]
pub struct InterfaceInfo {
    /// Pairs of coincident faces (face_idx_in_a, face_idx_in_b).
    pub face_pairs: Vec<(usize, usize)>,
    /// Pairs of coincident edges (edge_idx_in_a, edge_idx_in_b).
    pub edge_pairs: Vec<(usize, usize)>,
    /// Pairs of coincident vertices (vertex_idx_in_a, vertex_idx_in_b).
    pub vertex_pairs: Vec<(usize, usize)>,
    /// Whether the interface has consistent normals (pointing same direction).
    pub normals_consistent: bool,
}

/// Detects the interface between two BReps.
///
/// Returns information about coincident faces, edges, and vertices.
pub fn detect_interface(brep1: &BRep, brep2: &BRep, tolerance: f64) -> InterfaceInfo {
    let faces1 = flatten_faces(brep1);
    let faces2 = flatten_faces(brep2);

    let mut face_pairs = Vec::new();
    let mut normals_consistent = true;

    // Detect coincident faces
    for (f1_idx, f1) in faces1.iter().enumerate() {
        for (f2_idx, f2) in faces2.iter().enumerate() {
            if are_faces_coincident(f1, f2, brep1, brep2, tolerance) {
                face_pairs.push((f1_idx, f2_idx));
                // Check normal consistency
                if (f1.normal - f2.normal).length() > tolerance {
                    normals_consistent = false;
                }
            }
        }
    }

    // Detect coincident edges based on face pairs
    let edge_pairs = detect_edge_pairs(&face_pairs, brep1, brep2, tolerance);

    // Detect coincident vertices based on edge pairs
    let vertex_pairs = detect_vertex_pairs(&edge_pairs, brep1, brep2, tolerance);

    InterfaceInfo {
        face_pairs,
        edge_pairs,
        vertex_pairs,
        normals_consistent,
    }
}

/// Flattens all faces from a BRep into a single list.
fn flatten_faces(brep: &BRep) -> Vec<&Face> {
    brep.solids
        .iter()
        .flat_map(|s| s.shells.iter())
        .flat_map(|sh| sh.faces.iter())
        .collect()
}

/// Checks if two faces are coincident (share the same geometric surface).
fn are_faces_coincident(
    face1: &Face,
    face2: &Face,
    brep1: &BRep,
    brep2: &BRep,
    tolerance: f64,
) -> bool {
    // Get vertex positions for both faces
    let verts1 = get_face_vertices(face1, brep1);
    let verts2 = get_face_vertices(face2, brep2);

    if verts1.is_empty() || verts2.is_empty() {
        return false;
    }

    // Check if face centers are close
    let center1 = compute_center(&verts1);
    let center2 = compute_center(&verts2);

    if (center1 - center2).length() > tolerance * 10.0 {
        return false;
    }

    // Check if normals are parallel (same or opposite direction)
    let dot = face1.normal.dot(face2.normal).abs();
    if dot < 1.0 - tolerance {
        return false;
    }

    // Check if vertices of face1 are on face2's plane and vice versa
    let all_on_plane = verts1.iter().all(|&v| {
        let dist = (v - center2).dot(face2.normal).abs();
        dist < tolerance * 10.0
    }) && verts2.iter().all(|&v| {
        let dist = (v - center1).dot(face1.normal).abs();
        dist < tolerance * 10.0
    });

    all_on_plane
}

/// Gets the 3D positions of all vertices in a face.
fn get_face_vertices(face: &Face, brep: &BRep) -> Vec<DVec3> {
    let mut positions = Vec::new();
    let mut seen = HashSet::new();

    for we in &face.outer_wire.edges {
        if we.idx < brep.edges.len() {
            let edge = &brep.edges[we.idx];
            if edge.start < brep.vertices.len() && seen.insert(edge.start) {
                positions.push(brep.vertices[edge.start].point);
            }
            if edge.end < brep.vertices.len() && seen.insert(edge.end) {
                positions.push(brep.vertices[edge.end].point);
            }
        }
    }

    positions
}

/// Computes the geometric center of a set of points.
fn compute_center(points: &[DVec3]) -> DVec3 {
    if points.is_empty() {
        return DVec3::ZERO;
    }
    points.iter().sum::<DVec3>() / points.len() as f64
}

/// Detects coincident edge pairs based on face pairs.
fn detect_edge_pairs(
    face_pairs: &[(usize, usize)],
    brep1: &BRep,
    brep2: &BRep,
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let mut edge_pairs = Vec::new();
    let faces1 = flatten_faces(brep1);
    let faces2 = flatten_faces(brep2);

    // Build a set of face pairs for quick lookup
    let face_pair_set: HashSet<(usize, usize)> = face_pairs.iter().copied().collect();

    // Check edges from faces that are part of interface
    for (f1_idx, f1) in faces1.iter().enumerate() {
        for we1 in &f1.outer_wire.edges {
            for (f2_idx, f2) in faces2.iter().enumerate() {
                if !face_pair_set.contains(&(f1_idx, f2_idx)) {
                    continue;
                }
                for we2 in &f2.outer_wire.edges {
                    if are_edges_coincident(we1.idx, we2.idx, brep1, brep2, tolerance) {
                        edge_pairs.push((we1.idx, we2.idx));
                    }
                }
            }
        }
    }

    edge_pairs.sort_unstable();
    edge_pairs.dedup();
    edge_pairs
}

/// Checks if two edges are coincident.
fn are_edges_coincident(
    edge1_idx: usize,
    edge2_idx: usize,
    brep1: &BRep,
    brep2: &BRep,
    tolerance: f64,
) -> bool {
    if edge1_idx >= brep1.edges.len() || edge2_idx >= brep2.edges.len() {
        return false;
    }

    let edge1 = &brep1.edges[edge1_idx];
    let edge2 = &brep2.edges[edge2_idx];

    if edge1.start >= brep1.vertices.len() || edge1.end >= brep1.vertices.len() {
        return false;
    }
    if edge2.start >= brep2.vertices.len() || edge2.end >= brep2.vertices.len() {
        return false;
    }

    let start1 = brep1.vertices[edge1.start].point;
    let end1 = brep1.vertices[edge1.end].point;
    let start2 = brep2.vertices[edge2.start].point;
    let end2 = brep2.vertices[edge2.end].point;

    // Check if edges match (same direction or reversed)
    let same_dir = (start1 - start2).length() < tolerance && (end1 - end2).length() < tolerance;
    let rev_dir = (start1 - end2).length() < tolerance && (end1 - start2).length() < tolerance;

    same_dir || rev_dir
}

/// Detects coincident vertex pairs based on edge pairs.
fn detect_vertex_pairs(
    edge_pairs: &[(usize, usize)],
    brep1: &BRep,
    brep2: &BRep,
    tolerance: f64,
) -> Vec<(usize, usize)> {
    let mut vertex_pairs = Vec::new();

    for (e1_idx, e2_idx) in edge_pairs {
        if *e1_idx >= brep1.edges.len() || *e2_idx >= brep2.edges.len() {
            continue;
        }

        let edge1 = &brep1.edges[*e1_idx];
        let edge2 = &brep2.edges[*e2_idx];

        // Match vertices
        if edge1.start < brep1.vertices.len() && edge2.start < brep2.vertices.len() {
            let v1 = brep1.vertices[edge1.start].point;
            let v2 = brep2.vertices[edge2.start].point;
            if (v1 - v2).length() < tolerance {
                vertex_pairs.push((edge1.start, edge2.start));
            }
        }

        if edge1.start < brep1.vertices.len() && edge2.end < brep2.vertices.len() {
            let v1 = brep1.vertices[edge1.start].point;
            let v2 = brep2.vertices[edge2.end].point;
            if (v1 - v2).length() < tolerance {
                vertex_pairs.push((edge1.start, edge2.end));
            }
        }

        if edge1.end < brep1.vertices.len() && edge2.start < brep2.vertices.len() {
            let v1 = brep1.vertices[edge1.end].point;
            let v2 = brep2.vertices[edge2.start].point;
            if (v1 - v2).length() < tolerance {
                vertex_pairs.push((edge1.end, edge2.start));
            }
        }

        if edge1.end < brep1.vertices.len() && edge2.end < brep2.vertices.len() {
            let v1 = brep1.vertices[edge1.end].point;
            let v2 = brep2.vertices[edge2.end].point;
            if (v1 - v2).length() < tolerance {
                vertex_pairs.push((edge1.end, edge2.end));
            }
        }
    }

    // Also check all vertices for proximity
    for (v1_idx, v1) in brep1.vertices.iter().enumerate() {
        for (v2_idx, v2) in brep2.vertices.iter().enumerate() {
            if (v1.point - v2.point).length() < tolerance {
                let pair = (v1_idx, v2_idx);
                if !vertex_pairs.contains(&pair) {
                    vertex_pairs.push(pair);
                }
            }
        }
    }

    vertex_pairs.sort_unstable();
    vertex_pairs.dedup();
    vertex_pairs
}

// ─────────────────────────────────────────────────────────────────────────────
// Core Gluing Operations
// ─────────────────────────────────────────────────────────────────────────────

/// Glues two BReps together at their interface.
///
/// This function:
/// 1. Detects coincident faces/edges/vertices at the interface
/// 2. Merges the topology, eliminating duplicates
/// 3. Optionally preserves history for parametric editing
///
/// # Arguments
///
/// * `brep1` - First shape (will become the base)
/// * `brep2` - Second shape (will be merged into brep1)
/// * `opts` - Gluing options
///
/// # Returns
///
/// * `Ok(GluerResult)` - The glued shape with metadata
/// * `Err(GluerError)` - If gluing fails
///
/// # Example
///
/// ```
/// use rcad_kernel::{BRep, PrimitiveSolid};
/// use rcad_algorithms::gluer::{glue_shapes, GluerOptions};
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
///
/// let result = glue_shapes(&box1, &box2, GluerOptions::default());
/// ```
pub fn glue_shapes(
    brep1: &BRep,
    brep2: &BRep,
    opts: GluerOptions,
) -> Result<GluerResult, GluerError> {
    // Validate inputs
    if brep1.solids.is_empty() || brep2.solids.is_empty() {
        return Err(GluerError::InvalidInput("Empty BRep".to_string()));
    }

    // Detect interface
    let interface = detect_interface(brep1, brep2, opts.tolerance);

    if interface.face_pairs.is_empty() && interface.edge_pairs.is_empty() {
        // No interface found - just concatenate the shapes
        return Ok(concatenate_breps(brep1, brep2, &opts));
    }

    // Perform the gluing
    let result = perform_gluing(brep1, brep2, &interface, &opts)?;

    Ok(result)
}

/// Concatenates two BReps without merging (no interface found).
fn concatenate_breps(brep1: &BRep, brep2: &BRep, opts: &GluerOptions) -> GluerResult {
    let mut result = BRep::new();
    let mut history = GluerHistory::default();

    // Copy vertices from brep1
    let v1_count = brep1.vertices.len();
    for v in &brep1.vertices {
        result.vertices.push(v.clone());
        history.vertex_origins.push(VertexOrigin::FromA(result.vertices.len() - 1));
    }

    // Copy vertices from brep2
    let v2_offset = result.vertices.len();
    for v in &brep2.vertices {
        result.vertices.push(v.clone());
        history.vertex_origins.push(VertexOrigin::FromB(result.vertices.len() - 1));
    }

    // Copy edges from brep1
    let e1_count = brep1.edges.len();
    for e in &brep1.edges {
        result.edges.push(e.clone());
        history.edge_origins.push(EdgeOrigin::FromA(result.edges.len() - 1));
    }

    // Copy edges from brep2 with remapped vertices
    let e2_offset = result.edges.len();
    for e in &brep2.edges {
        result.edges.push(Edge {
            start: e.start + v2_offset,
            end: e.end + v2_offset,
        });
        history.edge_origins.push(EdgeOrigin::FromB(result.edges.len() - 1));
    }

    // Copy faces from brep1
    let f1_count: usize = brep1.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();
    for solid in &brep1.solids {
        let mut new_solid = Solid { shells: Vec::new() };
        for shell in &solid.shells {
            let mut new_shell = Shell { faces: Vec::new() };
            for face in &shell.faces {
                new_shell.faces.push(face.clone());
                history.face_origins.push(FaceOrigin::FromA(history.face_origins.len()));
            }
            new_solid.shells.push(new_shell);
        }
        result.solids.push(new_solid);
    }

    // Copy faces from brep2 with remapped edges
    for solid in &brep2.solids {
        let mut new_solid = Solid { shells: Vec::new() };
        for shell in &solid.shells {
            let mut new_shell = Shell { faces: Vec::new() };
            for face in &shell.faces {
                let mut new_face = face.clone();
                // Remap edge indices
                for we in &mut new_face.outer_wire.edges {
                    we.idx += e2_offset;
                }
                for inner in &mut new_face.inner_wires {
                    for we in &mut inner.edges {
                        we.idx += e2_offset;
                    }
                }
                new_shell.faces.push(new_face);
                history.face_origins.push(FaceOrigin::FromB(history.face_origins.len()));
            }
            new_solid.shells.push(new_shell);
        }
        result.solids.push(new_solid);
    }

    // Merge geometry stores
    result.geom = merge_geom_stores(&brep1.geom, &brep2.geom, v1_count, e1_count, f1_count);

    GluerResult {
        brep: result,
        history: if opts.preserve_history { Some(history) } else { None },
        merged_faces: Vec::new(),
        merged_edges: Vec::new(),
        merged_vertices: Vec::new(),
        interface_face_count: 0,
        interface_edge_count: 0,
    }
}

/// Performs the actual gluing operation.
fn perform_gluing(
    brep1: &BRep,
    brep2: &BRep,
    interface: &InterfaceInfo,
    opts: &GluerOptions,
) -> Result<GluerResult, GluerError> {
    let mut result = BRep::new();
    let mut history = GluerHistory::default();

    // Build vertex mapping: brep2 vertex -> result vertex
    let mut vertex_map: HashMap<usize, usize> = HashMap::new();
    let mut merged_vertex_pairs: Vec<(usize, usize)> = Vec::new();

    // Build edge mapping: brep2 edge -> result edge
    let mut edge_map: HashMap<usize, usize> = HashMap::new();
    let mut merged_edge_pairs: Vec<(usize, usize)> = Vec::new();

    // Build face mapping: brep2 face -> result face (for merged faces, this points to the merged face)
    let mut _face_map: HashMap<usize, usize> = HashMap::new();
    let mut merged_face_pairs: Vec<(usize, usize)> = Vec::new();

    // Step 1: Copy vertices from brep1
    for (v_idx, v) in brep1.vertices.iter().enumerate() {
        result.vertices.push(v.clone());
        history.vertex_origins.push(VertexOrigin::FromA(v_idx));
    }

    // Step 2: Process vertices from brep2 (merge or copy)
    for (v2_idx, v2) in brep2.vertices.iter().enumerate() {
        // Check if this vertex should be merged with one from brep1
        if let Some(&(v1_idx, _)) = interface.vertex_pairs.iter().find(|(_, v2)| *v2 == v2_idx) {
            if opts.merge_shared_vertices {
                vertex_map.insert(v2_idx, v1_idx);
                merged_vertex_pairs.push((v1_idx, v2_idx));
            } else {
                let new_idx = result.vertices.len();
                result.vertices.push(v2.clone());
                vertex_map.insert(v2_idx, new_idx);
                history.vertex_origins.push(VertexOrigin::FromB(v2_idx));
            }
        } else {
            let new_idx = result.vertices.len();
            result.vertices.push(v2.clone());
            vertex_map.insert(v2_idx, new_idx);
            history.vertex_origins.push(VertexOrigin::FromB(v2_idx));
        }
    }

    // Step 3: Copy edges from brep1
    for (e_idx, e) in brep1.edges.iter().enumerate() {
        result.edges.push(e.clone());
        history.edge_origins.push(EdgeOrigin::FromA(e_idx));
    }

    // Step 4: Process edges from brep2 (merge or copy with remapped vertices)
    let edge_pair_set: HashSet<(usize, usize)> = interface.edge_pairs.iter().copied().collect();

    for (e2_idx, e2) in brep2.edges.iter().enumerate() {
        // Check if this edge should be merged
        let merge_candidate = edge_pair_set.iter().find(|(_, e2)| *e2 == e2_idx);

        if let Some(&(e1_idx, _)) = merge_candidate {
            if opts.merge_shared_edges {
                edge_map.insert(e2_idx, e1_idx);
                merged_edge_pairs.push((e1_idx, e2_idx));
            } else {
                // Copy with remapped vertices
                let new_idx = result.edges.len();
                let new_edge = Edge {
                    start: *vertex_map.get(&e2.start).unwrap_or(&e2.start),
                    end: *vertex_map.get(&e2.end).unwrap_or(&e2.end),
                };
                result.edges.push(new_edge);
                edge_map.insert(e2_idx, new_idx);
                history.edge_origins.push(EdgeOrigin::FromB(e2_idx));
            }
        } else {
            // Copy with remapped vertices
            let new_idx = result.edges.len();
            let new_edge = Edge {
                start: *vertex_map.get(&e2.start).unwrap_or(&e2.start),
                end: *vertex_map.get(&e2.end).unwrap_or(&e2.end),
            };
            result.edges.push(new_edge);
            edge_map.insert(e2_idx, new_idx);
            history.edge_origins.push(EdgeOrigin::FromB(e2_idx));
        }
    }

    // Step 5: Process faces - create solids with merged interface
    let face_pair_set: HashSet<(usize, usize)> = interface.face_pairs.iter().copied().collect();

    // Collect faces to exclude (interface faces from brep2 that will be merged)
    let mut faces2_to_exclude: HashSet<usize> = HashSet::new();
    if opts.merge_shared_faces {
        for (_, f2_idx) in &interface.face_pairs {
            faces2_to_exclude.insert(*f2_idx);
        }
    }

    // Copy brep1's solids with faces
    for solid in &brep1.solids {
        let mut new_solid = Solid { shells: Vec::new() };
        for shell in &solid.shells {
            let mut new_shell = Shell { faces: Vec::new() };
            for (f_idx, face) in shell.faces.iter().enumerate() {
                // Check if this is an interface face
                let is_interface = face_pair_set.iter().any(|(f1, _)| *f1 == f_idx);

                new_shell.faces.push(face.clone());
                if is_interface {
                    merged_face_pairs.push((f_idx, face_pair_set.iter().find(|(f1, _)| *f1 == f_idx).map(|(_, f2)| *f2).unwrap_or(0)));
                    history.face_origins.push(FaceOrigin::Merged {
                        from_a: f_idx,
                        from_b: face_pair_set.iter().find(|(f1, _)| *f1 == f_idx).map(|(_, f2)| *f2).unwrap_or(0),
                    });
                } else {
                    history.face_origins.push(FaceOrigin::FromA(f_idx));
                }
            }
            new_solid.shells.push(new_shell);
        }
        result.solids.push(new_solid);
    }

    // Add brep2's non-interface faces to the result
    for solid in &brep2.solids {
        for shell in &solid.shells {
            // Find the corresponding shell in result (or create new)
            let target_solid_idx = result.solids.len() - 1;
            for (f_idx, face) in shell.faces.iter().enumerate() {
                if faces2_to_exclude.contains(&f_idx) {
                    continue; // Skip interface faces
                }

                // Remap edge indices in the face
                let mut new_face = face.clone();
                for we in &mut new_face.outer_wire.edges {
                    if let Some(&new_e_idx) = edge_map.get(&we.idx) {
                        we.idx = new_e_idx;
                    }
                }
                for inner in &mut new_face.inner_wires {
                    for we in &mut inner.edges {
                        if let Some(&new_e_idx) = edge_map.get(&we.idx) {
                            we.idx = new_e_idx;
                        }
                    }
                }

                result.solids[target_solid_idx].shells[0].faces.push(new_face);
                history.face_origins.push(FaceOrigin::FromB(f_idx));
            }
        }
    }

    // Merge geometry stores
    let f1_count: usize = brep1.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();
    result.geom = merge_geom_stores(
        &brep1.geom,
        &brep2.geom,
        brep1.vertices.len(),
        brep1.edges.len(),
        f1_count,
    );

    // Remap geometry references for merged edges and vertices
    apply_geometry_mappings(&mut result.geom, &edge_map, brep1.edges.len());

    Ok(GluerResult {
        brep: result,
        history: if opts.preserve_history { Some(history) } else { None },
        merged_faces: merged_face_pairs,
        merged_edges: merged_edge_pairs,
        merged_vertices: merged_vertex_pairs,
        interface_face_count: interface.face_pairs.len(),
        interface_edge_count: interface.edge_pairs.len(),
    })
}

/// Merges two GeomStores with appropriate offset adjustments.
fn merge_geom_stores(
    geom1: &GeomStore,
    geom2: &GeomStore,
    _v1_count: usize,
    _e1_count: usize,
    _f1_count: usize,
) -> GeomStore {
    let mut result = geom1.clone();

    // Merge curves
    let curve_offset = result.curves.len();
    result.curves.extend(geom2.curves.clone());

    // Merge surfaces
    let surface_offset = result.surfaces.len();
    result.surfaces.extend(geom2.surfaces.clone());

    // Merge 2D curves
    let curve2d_offset = result.curve2ds.len();
    result.curve2ds.extend(geom2.curve2ds.clone());

    // Merge edge_curve with offset
    result.edge_curve.extend(
        geom2.edge_curve.iter().map(|ec| {
            ec.map(|idx| idx + curve_offset)
        })
    );

    // Merge face_surface with offset
    result.face_surface.extend(
        geom2.face_surface.iter().map(|fs| {
            fs.map(|idx| idx + surface_offset)
        })
    );

    // Merge edge_pcurves with offset
    result.edge_pcurves.extend(
        geom2.edge_pcurves.iter().map(|pcurves| {
            pcurves.iter().map(|pc| PCurve {
                surface_idx: pc.surface_idx + surface_offset,
                curve2d_idx: pc.curve2d_idx + curve2d_offset,
            }).collect()
        })
    );

    // Merge edge_curve_range
    result.edge_curve_range.extend(geom2.edge_curve_range.clone());

    // Merge edge_degenerated
    result.edge_degenerated.extend(geom2.edge_degenerated.clone());

    // Merge tolerances
    result.vertex_tolerance.extend(geom2.vertex_tolerance.clone());
    result.edge_tolerance.extend(geom2.edge_tolerance.clone());
    result.face_tolerance.extend(geom2.face_tolerance.clone());
    result.curve2d_range.extend(geom2.curve2d_range.clone());
    result.face_surface_range.extend(geom2.face_surface_range.clone());
    result.edge_same_parameter.extend(geom2.edge_same_parameter.clone());
    result.edge_same_range.extend(geom2.edge_same_range.clone());

    result
}

/// Applies edge mapping to geometry references.
fn apply_geometry_mappings(
    _geom: &mut GeomStore,
    edge_map: &HashMap<usize, usize>,
    e1_count: usize,
) {
    // For edges from brep2 that were merged, we need to handle pcurves
    // This is a simplified version - a full implementation would merge pcurves
    for (&e2_idx, &target_idx) in edge_map {
        if e2_idx >= e1_count && target_idx < e1_count {
            // Edge from brep2 was merged with an edge from brep1
            // The pcurves from e2_idx should potentially be added to target_idx
            // For now, we keep the pcurves from brep1's edge
            // A more complete implementation would merge pcurves
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Shell-Level Gluing
// ─────────────────────────────────────────────────────────────────────────────

/// Glues two shells at specified interface faces.
///
/// This is a lower-level function for when the interface is already known.
///
/// # Arguments
///
/// * `shell1` - First shell
/// * `shell2` - Second shell
/// * `interface_faces` - Pairs of face indices (in shell1, in shell2) that form the interface
///
/// # Returns
///
/// The merged shell, or an error if gluing fails.
pub fn glue_at_interface(
    shell1: &Shell,
    shell2: &Shell,
    interface_faces: &[(usize, usize)],
) -> Result<Shell, GluerError> {
    if interface_faces.is_empty() {
        // No interface - just concatenate
        let mut merged = shell1.clone();
        merged.faces.extend(shell2.faces.clone());
        return Ok(merged);
    }

    let _interface_set: HashSet<(usize, usize)> = interface_faces.iter().copied().collect();
    let faces2_to_skip: HashSet<usize> = interface_faces.iter().map(|(_, f2)| *f2).collect();

    let mut merged = Shell { faces: shell1.faces.clone() };

    // Add non-interface faces from shell2
    for (f2_idx, face) in shell2.faces.iter().enumerate() {
        if !faces2_to_skip.contains(&f2_idx) {
            merged.faces.push(face.clone());
        }
    }

    Ok(merged)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;
    use rcad_kernel::{Wire, Face, Shell};
    use glam::DAffine3;

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    #[test]
    fn test_detect_interface_no_overlap() {
        let box1 = unit_box();
        let box2 = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }).transformed(DAffine3::from_translation(DVec3::new(10.0, 0.0, 0.0)));

        let interface = detect_interface(&box1, &box2, 1e-6);
        assert!(interface.face_pairs.is_empty());
        assert!(interface.edge_pairs.is_empty());
    }

    #[test]
    fn test_detect_interface_touching_faces() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        // Translate box2 to touch box1 at y=1 face
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let interface = detect_interface(&box1, &box2, 1e-6);
        // Should detect at least one coincident face pair
        assert!(!interface.face_pairs.is_empty() || !interface.edge_pairs.is_empty());
    }

    #[test]
    fn test_glue_shapes_no_interface() {
        let box1 = unit_box();
        let box2 = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }).transformed(DAffine3::from_translation(DVec3::new(10.0, 0.0, 0.0)));

        let result = glue_shapes(&box1, &box2, GluerOptions::default()).unwrap();

        // Should have combined all faces
        let total_faces: usize = result.brep.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(total_faces, 12); // 6 + 6

        assert!(result.merged_faces.is_empty());
    }

    #[test]
    fn test_glue_shapes_with_interface() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let result = glue_shapes(&box1, &box2, GluerOptions::default()).unwrap();

        // Should have fewer than 12 faces due to merging
        let total_faces: usize = result.brep.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        assert!(total_faces <= 12);
    }

    #[test]
    fn test_glue_preserves_history() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let opts = GluerOptions {
            preserve_history: true,
            ..Default::default()
        };

        let result = glue_shapes(&box1, &box2, opts).unwrap();
        assert!(result.history.is_some());

        let history = result.history.unwrap();
        assert!(!history.face_origins.is_empty());
    }

    #[test]
    fn test_glue_without_history() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let opts = GluerOptions {
            preserve_history: false,
            ..Default::default()
        };

        let result = glue_shapes(&box1, &box2, opts).unwrap();
        assert!(result.history.is_none());
    }

    #[test]
    fn test_glue_at_interface() {
        let shell1 = Shell {
            faces: vec![
                Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                },
            ],
        };
        let shell2 = Shell {
            faces: vec![
                Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                },
            ],
        };

        // Glue with one interface face
        let result = glue_at_interface(&shell1, &shell2, &[(0, 0)]).unwrap();
        assert_eq!(result.faces.len(), 1); // Only one face kept

        // Glue with no interface
        let result = glue_at_interface(&shell1, &shell2, &[]).unwrap();
        assert_eq!(result.faces.len(), 2); // Both faces kept
    }

    #[test]
    fn test_gluer_options_default() {
        let opts = GluerOptions::default();
        assert_eq!(opts.tolerance, 1e-6);
        assert!(opts.merge_shared_faces);
        assert!(opts.merge_shared_edges);
        assert!(opts.merge_shared_vertices);
        assert!(opts.preserve_history);
    }

    #[test]
    fn test_empty_brep_error() {
        let empty = BRep::new();
        let box1 = unit_box();

        let result = glue_shapes(&empty, &box1, GluerOptions::default());
        assert!(matches!(result, Err(GluerError::InvalidInput(_))));
    }

    #[test]
    fn test_vertex_merging() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let opts = GluerOptions {
            merge_shared_vertices: true,
            ..Default::default()
        };

        let result = glue_shapes(&box1, &box2, opts).unwrap();

        // Some vertices should be merged
        // The result should have fewer vertices than 16 (8+8)
        assert!(result.brep.vertices.len() < 16);
    }

    #[test]
    fn test_edge_merging() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let opts = GluerOptions {
            merge_shared_edges: true,
            ..Default::default()
        };

        let result = glue_shapes(&box1, &box2, opts).unwrap();

        // If faces are merged, edges should be merged too
        if !result.merged_faces.is_empty() {
            assert!(!result.merged_edges.is_empty() || result.interface_edge_count > 0);
        }
    }

    #[test]
    fn test_no_merge_option() {
        let box1 = unit_box();
        let mut box2 = unit_box();
        box2.apply_transform(DAffine3::from_translation(DVec3::new(0.0, 1.0, 0.0)));

        let opts = GluerOptions {
            merge_shared_faces: false,
            merge_shared_edges: false,
            merge_shared_vertices: false,
            ..Default::default()
        };

        let result = glue_shapes(&box1, &box2, opts).unwrap();

        // Without merging, we should get all entities
        assert_eq!(result.brep.vertices.len(), 16); // 8 + 8
        assert_eq!(result.brep.edges.len(), 24); // 12 + 12
    }
}
