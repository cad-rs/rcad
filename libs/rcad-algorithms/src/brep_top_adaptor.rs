//! BRepTopAdaptor-style topology adapters for high-level topology traversal.
//!
//! This module provides high-level adapters for exploring BRep topology,
//! analogous to OCCT's `TopExp_Explorer` and `BRepAdaptor` classes.
//!
//! # Overview
//!
//! - **Explorers**: `FaceExplorer`, `EdgeExplorer`, `VertexExplorer`, `WireExplorer`
//!   provide forward-only iteration over topology elements.
//! - **ShapeIterator**: Generic iterator implementing `std::iter::Iterator` for all shape types.
//! - **Topology queries**: Helper functions for adjacency queries.
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::brep_top_adaptor::*;
//! use rcad_kernel::BRep;
//!
//! let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
//!     width: 1.0, height: 1.0, depth: 1.0
//! });
//!
//! // Count faces using FaceExplorer
//! let mut explorer = FaceExplorer::new(&brep);
//! let mut face_count = 0;
//! while explorer.next().is_some() {
//!     face_count += 1;
//! }
//! assert_eq!(face_count, 6);
//! ```

pub use crate::brep_tools::ShapeType;
use crate::tolerance::*;
use rcad_kernel::topods::{self, TShape};
use rcad_kernel::topology::Face;

// =============================================================================
// Face Adaptor
// =============================================================================

/// Adapter providing convenient access to face data.
///
/// Analogous to OCCT `BRepAdaptor_Face`.
#[derive(Debug, Clone)]
pub struct FaceAdaptor<'a> {
    brep: &'a rcad_kernel::BRep,
    face_idx: usize, // tshape index of TShape::Face
}

impl<'a> FaceAdaptor<'a> {
    /// Create a new face adaptor.
    pub fn new(brep: &'a rcad_kernel::BRep, face_idx: usize) -> Self {
        Self { brep, face_idx }
    }

    /// Returns the face index (tshape index of the TShape::Face).
    pub fn index(&self) -> usize {
        self.face_idx
    }

    /// Returns the shell index containing this face (always 0 in flat TShape model).
    pub fn shell_index(&self) -> usize {
        0
    }

    /// Returns the solid index containing this face (always 0 in flat TShape model).
    pub fn solid_index(&self) -> usize {
        0
    }

    fn get_face_data(&self) -> Option<&topods::TFaceData> {
        let ts = self.brep.tshapes.get(self.face_idx)?;
        match ts.as_ref() {
            TShape::Face(fd) => Some(fd),
            _ => None,
        }
    }

    /// Returns a reference to the face topology (always None in the new TShape model).
    pub fn face(&self) -> Option<&'a Face> {
        None
    }

    /// Returns the surface index for this face, if available (now returns tshape index itself).
    pub fn surface_index(&self) -> Option<usize> {
        self.get_face_data()
            .and_then(|fd| fd.surface.as_ref().map(|_| self.face_idx))
    }

    /// Returns the number of edges in the outer wire.
    pub fn edge_count(&self) -> usize {
        self.get_face_data()
            .and_then(|fd| {
                let wts = self.brep.tshapes.get(fd.outer_wire.index)?;
                match wts.as_ref() {
                    TShape::Wire(wd) => Some(wd.edges.len()),
                    _ => None,
                }
            })
            .unwrap_or(0)
    }

    /// Returns the number of inner wires (holes).
    pub fn inner_wire_count(&self) -> usize {
        self.get_face_data()
            .map(|fd| fd.inner_wires.len())
            .unwrap_or(0)
    }

    /// Returns the face tolerance.
    pub fn tolerance(&self) -> f64 {
        TOLERANCE_MESH_LEGACY
    }
}

// =============================================================================
// Edge Adaptor
// =============================================================================

/// Adapter providing convenient access to edge data.
///
/// Analogous to OCCT `BRepAdaptor_Curve` / `BRepAdaptor_Edge`.
#[derive(Debug, Clone)]
pub struct EdgeAdaptor<'a> {
    brep: &'a rcad_kernel::BRep,
    edge_idx: usize, // tshape index of TShape::Edge
}

impl<'a> EdgeAdaptor<'a> {
    fn get_edge_data(&self) -> Option<&topods::TEdgeData> {
        let ts = self.brep.tshapes.get(self.edge_idx)?;
        match ts.as_ref() {
            TShape::Edge(ed) => Some(ed),
            _ => None,
        }
    }

    /// Create a new edge adaptor.
    pub fn new(brep: &'a rcad_kernel::BRep, edge_idx: usize) -> Self {
        Self { brep, edge_idx }
    }

    /// Returns the edge index (tshape index).
    pub fn index(&self) -> usize {
        self.edge_idx
    }

    /// Returns the edge topology (start and end vertex indices) — not available in new model.
    pub fn edge(&self) -> Option<rcad_kernel::topology::Edge> {
        None
    }

    /// Returns the start vertex index.
    pub fn start_vertex(&self) -> Option<usize> {
        self.get_edge_data().map(|ed| ed.first.index)
    }

    /// Returns the end vertex index.
    pub fn end_vertex(&self) -> Option<usize> {
        self.get_edge_data().map(|ed| ed.last.index)
    }

    /// Returns the 3D curve reference, if available.
    pub fn curve_index(&self) -> Option<usize> {
        // In the new model, curves are stored directly on TEdgeData.
        // Return the edge index if a curve exists, None otherwise.
        self.get_edge_data()
            .and_then(|ed| ed.curve.as_ref().map(|_| self.edge_idx))
    }

    /// Returns the parameter range for this edge.
    pub fn parameter_range(&self) -> Option<[f64; 2]> {
        self.get_edge_data().map(|ed| ed.range)
    }

    /// Returns true if this edge is degenerate.
    pub fn is_degenerate(&self) -> bool {
        self.get_edge_data()
            .map(|ed| ed.degenerated)
            .unwrap_or(false)
    }

    /// Returns true if this edge is closed (start == end vertex).
    pub fn is_closed(&self) -> bool {
        self.get_edge_data()
            .map(|ed| ed.first.index == ed.last.index)
            .unwrap_or(false)
    }

    /// Returns the edge tolerance.
    pub fn tolerance(&self) -> f64 {
        self.get_edge_data()
            .map(|ed| ed.tolerance)
            .filter(|&t| t > 0.0)
            .unwrap_or(TOLERANCE_MESH_LEGACY)
    }
}

// =============================================================================
// Vertex Adaptor
// =============================================================================

/// Adapter providing convenient access to vertex data.
///
/// Analogous to OCCT `BRepAdaptor_Point` / `BRep_Tool` for vertices.
#[derive(Debug, Clone)]
pub struct VertexAdaptor<'a> {
    brep: &'a rcad_kernel::BRep,
    vertex_idx: usize, // tshape index of TShape::Vertex
}

impl<'a> VertexAdaptor<'a> {
    /// Create a new vertex adaptor.
    pub fn new(brep: &'a rcad_kernel::BRep, vertex_idx: usize) -> Self {
        Self { brep, vertex_idx }
    }

    /// Returns the vertex index (tshape index).
    pub fn index(&self) -> usize {
        self.vertex_idx
    }

    /// Returns the 3D point location of the vertex.
    pub fn point(&self) -> Option<glam::DVec3> {
        self.brep.vertex_point(self.vertex_idx)
    }

    /// Returns the vertex tolerance.
    pub fn tolerance(&self) -> f64 {
        self.brep
            .tshapes
            .get(self.vertex_idx)
            .and_then(|ts| {
                if let TShape::Vertex(vd) = ts.as_ref() {
                    Some(vd.tolerance)
                } else {
                    None
                }
            })
            .filter(|&t| t > 0.0)
            .unwrap_or(TOLERANCE_MESH_LEGACY)
    }
}

// =============================================================================
// Face Explorer
// =============================================================================

/// Forward-only explorer for faces in a BRep.
///
/// Analogous to OCCT `TopExp_Explorer(shape, TopAbs_FACE)`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::FaceExplorer;
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// let mut explorer = FaceExplorer::new(&brep);
/// let mut faces = Vec::new();
/// while let Some(idx) = explorer.next() {
///     faces.push(idx);
/// }
/// assert_eq!(faces.len(), 6);
/// ```
#[derive(Debug, Clone)]
pub struct FaceExplorer<'a> {
    brep: &'a rcad_kernel::BRep,
    tshape_idx: usize,
    current: Option<FaceAdaptor<'a>>,
}

impl<'a> FaceExplorer<'a> {
    /// Create a new face explorer.
    pub fn new(brep: &'a rcad_kernel::BRep) -> Self {
        Self {
            brep,
            tshape_idx: 0,
            current: None,
        }
    }

    /// Advance to the next face and return its tshape index.
    ///
    /// Returns `None` when all faces have been visited.
    pub fn next(&mut self) -> Option<usize> {
        while let Some(ts) = self.brep.tshapes.get(self.tshape_idx) {
            let idx = self.tshape_idx;
            self.tshape_idx += 1;
            if matches!(ts.as_ref(), TShape::Face(_)) {
                self.current = Some(FaceAdaptor::new(self.brep, idx));
                return Some(idx);
            }
        }
        self.current = None;
        None
    }

    /// Returns the adaptor for the current face.
    ///
    /// Only valid after a successful call to `next()`.
    pub fn current_adaptor(&self) -> Option<&FaceAdaptor<'a>> {
        self.current.as_ref()
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.tshape_idx = 0;
        self.current = None;
    }
}

// =============================================================================
// Edge Explorer
// =============================================================================

/// Forward-only explorer for edges in a BRep.
///
/// Analogous to OCCT `TopExp_Explorer(shape, TopAbs_EDGE)`.
#[derive(Debug, Clone)]
pub struct EdgeExplorer<'a> {
    brep: &'a rcad_kernel::BRep,
    tshape_idx: usize,
    current: Option<EdgeAdaptor<'a>>,
}

impl<'a> EdgeExplorer<'a> {
    /// Create a new edge explorer.
    pub fn new(brep: &'a rcad_kernel::BRep) -> Self {
        Self {
            brep,
            tshape_idx: 0,
            current: None,
        }
    }

    /// Advance to the next edge and return its tshape index.
    ///
    /// Returns `None` when all edges have been visited.
    pub fn next(&mut self) -> Option<usize> {
        while let Some(ts) = self.brep.tshapes.get(self.tshape_idx) {
            let idx = self.tshape_idx;
            self.tshape_idx += 1;
            if matches!(ts.as_ref(), TShape::Edge(_)) {
                self.current = Some(EdgeAdaptor::new(self.brep, idx));
                return Some(idx);
            }
        }
        self.current = None;
        None
    }

    /// Returns the adaptor for the current edge.
    ///
    /// Only valid after a successful call to `next()`.
    pub fn current_adaptor(&self) -> Option<&EdgeAdaptor<'a>> {
        self.current.as_ref()
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.tshape_idx = 0;
        self.current = None;
    }
}

// =============================================================================
// Vertex Explorer
// =============================================================================

/// Forward-only explorer for vertices in a BRep.
///
/// Analogous to OCCT `TopExp_Explorer(shape, TopAbs_VERTEX)`.
#[derive(Debug, Clone)]
pub struct VertexExplorer<'a> {
    brep: &'a rcad_kernel::BRep,
    tshape_idx: usize,
}

impl<'a> VertexExplorer<'a> {
    /// Create a new vertex explorer.
    pub fn new(brep: &'a rcad_kernel::BRep) -> Self {
        Self {
            brep,
            tshape_idx: 0,
        }
    }

    /// Advance to the next vertex and return its tshape index.
    ///
    /// Returns `None` when all vertices have been visited.
    pub fn next(&mut self) -> Option<usize> {
        while let Some(ts) = self.brep.tshapes.get(self.tshape_idx) {
            let idx = self.tshape_idx;
            self.tshape_idx += 1;
            if matches!(ts.as_ref(), TShape::Vertex(_)) {
                return Some(idx);
            }
        }
        None
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.tshape_idx = 0;
    }
}

// =============================================================================
// Wire Explorer
// =============================================================================

/// An edge reference with orientation from wire traversal.
#[derive(Debug, Clone, Copy)]
pub struct OrientedEdge {
    /// Edge index in `BRep.edges`.
    pub edge_idx: usize,
    /// True if traversed in forward direction (start to end).
    pub forward: bool,
}

impl OrientedEdge {
    /// Create a new oriented edge.
    pub fn new(edge_idx: usize, forward: bool) -> Self {
        Self { edge_idx, forward }
    }
}

/// Forward-only explorer for edges in a wire (face boundary).
///
/// Analogous to OCCT `TopExp_Explorer(face, TopAbs_EDGE)` or `BRepTools_WireExplorer`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::WireExplorer;
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Explore the wire of the first face
/// let mut explorer = WireExplorer::new(&brep, 0);
/// let mut edge_count = 0;
/// while explorer.next().is_some() {
///     edge_count += 1;
/// }
/// assert_eq!(edge_count, 4); // Each box face has 4 edges
/// ```
#[derive(Debug, Clone)]
pub struct WireExplorer<'a> {
    brep: &'a rcad_kernel::BRep,
    face_idx: usize,
    wire_idx: usize, // 0 = outer wire, 1+ = inner wires
    edge_idx: usize,
    current: Option<OrientedEdge>,
}

impl<'a> WireExplorer<'a> {
    /// Create a new wire explorer for a specific face.
    ///
    /// `face_idx` is the tshape index of the TShape::Face.
    pub fn new(brep: &'a rcad_kernel::BRep, face_idx: usize) -> Self {
        Self {
            brep,
            face_idx,
            wire_idx: 0,
            edge_idx: 0,
            current: None,
        }
    }

    /// Get the edge refs for a given wire ShapeRef.
    fn get_wire_edges(&self, wire_ref: &topods::ShapeRef) -> Vec<(usize, bool)> {
        let Some(wts) = self.brep.tshapes.get(wire_ref.index) else {
            return Vec::new();
        };
        let TShape::Wire(wd) = wts.as_ref() else {
            return Vec::new();
        };
        wd.edges
            .iter()
            .map(|er| {
                let forward = er.orientation == topods::Orientation::Forward;
                (er.index, forward)
            })
            .collect()
    }

    /// Get all wire edge groups for the face: first outer wire, then inner wires.
    fn get_face_wires(&self) -> Vec<Vec<(usize, bool)>> {
        let Some(ts) = self.brep.tshapes.get(self.face_idx) else {
            return Vec::new();
        };
        let TShape::Face(fd) = ts.as_ref() else {
            return Vec::new();
        };
        let mut wires = Vec::new();
        wires.push(self.get_wire_edges(&fd.outer_wire));
        for inner in &fd.inner_wires {
            wires.push(self.get_wire_edges(inner));
        }
        wires
    }

    /// Advance to the next edge and return its oriented reference.
    ///
    /// Returns `None` when all edges (outer and inner wires) have been visited.
    pub fn next(&mut self) -> Option<OrientedEdge> {
        let wires = self.get_face_wires();
        if wires.is_empty() {
            return None;
        }

        loop {
            let wire = wires.get(self.wire_idx)?;

            if self.edge_idx < wire.len() {
                let (idx, forward) = wire[self.edge_idx];
                self.current = Some(OrientedEdge::new(idx, forward));
                self.edge_idx += 1;
                return self.current;
            } else {
                // Move to next wire
                self.wire_idx += 1;
                self.edge_idx = 0;
            }
        }
    }

    /// Returns the current oriented edge.
    pub fn current(&self) -> Option<OrientedEdge> {
        self.current
    }

    /// Returns true if currently iterating over the outer wire.
    pub fn is_outer_wire(&self) -> bool {
        self.wire_idx == 0
    }

    /// Returns the current wire index (0 = outer, 1+ = inner).
    pub fn wire_index(&self) -> usize {
        self.wire_idx
    }

    /// Reset the explorer to start from the beginning.
    pub fn reset(&mut self) {
        self.wire_idx = 0;
        self.edge_idx = 0;
        self.current = None;
    }
}

// =============================================================================
// Shape Iterator
// =============================================================================

/// Internal state for shape iteration.
#[derive(Debug, Clone)]
struct ShapeIterState {
    solid_idx: usize,
    shell_idx: usize,
    face_idx: usize,
    wire_idx: usize,
    edge_idx: usize,
    vertex_idx: usize,
}

impl ShapeIterState {
    fn new() -> Self {
        Self {
            solid_idx: 0,
            shell_idx: 0,
            face_idx: 0,
            wire_idx: 0,
            edge_idx: 0,
            vertex_idx: 0,
        }
    }
}

/// Generic iterator over shapes of a specific type in a BRep.
///
/// Implements `std::iter::Iterator` for idiomatic iteration.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::{ShapeIterator, ShapeType};
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Iterate over all faces
/// let faces: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Face).collect();
/// assert_eq!(faces.len(), 6);
///
/// // Iterate over all edges
/// let edges: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Edge).collect();
/// assert_eq!(edges.len(), 12);
///
/// // Iterate over all vertices
/// let vertices: Vec<usize> = ShapeIterator::new(&brep, ShapeType::Vertex).collect();
/// assert_eq!(vertices.len(), 8);
/// ```
#[derive(Debug, Clone)]
pub struct ShapeIterator<'a> {
    brep: &'a rcad_kernel::BRep,
    shape_type: ShapeType,
    state: ShapeIterState,
    done: bool,
}

impl<'a> ShapeIterator<'a> {
    /// Create a new shape iterator for the given shape type.
    pub fn new(brep: &'a rcad_kernel::BRep, shape_type: ShapeType) -> Self {
        Self {
            brep,
            shape_type,
            state: ShapeIterState::new(),
            done: false,
        }
    }
}

impl<'a> Iterator for ShapeIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        // Fast path: iterate tshapes by shape type, returning tshape indices.
        let candidates: Vec<usize> = self
            .brep
            .tshapes
            .iter()
            .enumerate()
            .filter(|(_, ts)| match self.shape_type {
                ShapeType::Vertex => matches!(ts.as_ref(), TShape::Vertex(_)),
                ShapeType::Edge => matches!(ts.as_ref(), TShape::Edge(_)),
                ShapeType::Wire => matches!(ts.as_ref(), TShape::Wire(_)),
                ShapeType::Face => matches!(ts.as_ref(), TShape::Face(_)),
                ShapeType::Shell => matches!(ts.as_ref(), TShape::Shell(_)),
                ShapeType::Solid => matches!(ts.as_ref(), TShape::Solid(_)),
                ShapeType::Compound | ShapeType::CompSolid | ShapeType::Empty => false,
            })
            .map(|(idx, _)| idx)
            .collect();

        // Use the appropriate counter based on shape type
        let next_idx = match self.shape_type {
            ShapeType::Vertex => {
                if self.state.vertex_idx < candidates.len() {
                    let idx = candidates[self.state.vertex_idx];
                    self.state.vertex_idx += 1;
                    Some(idx)
                } else {
                    None
                }
            }
            ShapeType::Edge => {
                if self.state.edge_idx < candidates.len() {
                    let idx = candidates[self.state.edge_idx];
                    self.state.edge_idx += 1;
                    Some(idx)
                } else {
                    None
                }
            }
            ShapeType::Face => {
                if self.state.face_idx < candidates.len() {
                    let idx = candidates[self.state.face_idx];
                    self.state.face_idx += 1;
                    Some(idx)
                } else {
                    None
                }
            }
            ShapeType::Shell => {
                if self.state.shell_idx < candidates.len() {
                    let idx = candidates[self.state.shell_idx];
                    self.state.shell_idx += 1;
                    Some(idx)
                } else {
                    None
                }
            }
            ShapeType::Solid => {
                if self.state.solid_idx < candidates.len() {
                    let idx = candidates[self.state.solid_idx];
                    self.state.solid_idx += 1;
                    Some(idx)
                } else {
                    None
                }
            }
            ShapeType::Wire => {
                // Wires are TShape::Wire entries in tshapes.
                if self.state.wire_idx < candidates.len() {
                    let idx = candidates[self.state.wire_idx];
                    self.state.wire_idx += 1;
                    Some(idx)
                } else {
                    None
                }
            }
            ShapeType::Compound | ShapeType::CompSolid | ShapeType::Empty => None,
        };

        if next_idx.is_none() {
            self.done = true;
        }
        next_idx
    }
}

// =============================================================================
// Topology Queries (tshape-based)
// =============================================================================

/// Returns the edge indices from all wires of a face (by tshape index).
fn collect_face_edge_indices(brep: &rcad_kernel::BRep, face_idx: usize) -> Vec<usize> {
    let mut result = Vec::new();
    let ts = match brep.tshapes.get(face_idx) {
        Some(ts) => ts,
        None => return result,
    };
    let TShape::Face(fd) = ts.as_ref() else {
        return result;
    };

    // Collect from outer wire
    if let Some(wts) = brep.tshapes.get(fd.outer_wire.index) {
        if let TShape::Wire(wd) = wts.as_ref() {
            for er in &wd.edges {
                result.push(er.index);
            }
        }
    }
    // Collect from inner wires
    for iw in &fd.inner_wires {
        if let Some(wts) = brep.tshapes.get(iw.index) {
            if let TShape::Wire(wd) = wts.as_ref() {
                for er in &wd.edges {
                    result.push(er.index);
                }
            }
        }
    }
    result
}

/// Returns all edge indices referenced by a face (including inner wires).
///
/// Duplicate edge indices are preserved as they appear in the wire.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::edges_of_face;
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// let edges = edges_of_face(&brep, 0);
/// assert_eq!(edges.len(), 4);
/// ```
pub fn edges_of_face(brep: &rcad_kernel::BRep, face_idx: usize) -> Vec<usize> {
    collect_face_edge_indices(brep, face_idx)
}

/// Returns all face indices that reference the given edge.
///
/// For a manifold solid, each edge is typically shared by 2 faces.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::faces_of_edge;
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Each edge of a box is shared by exactly 2 faces
/// let faces = faces_of_edge(&brep, 0);
/// assert_eq!(faces.len(), 2);
/// ```
pub fn faces_of_edge(brep: &rcad_kernel::BRep, edge_idx: usize) -> Vec<usize> {
    let mut result = Vec::new();
    for (fi, ts) in brep.tshapes.iter().enumerate() {
        let TShape::Face(fd) = ts.as_ref() else {
            continue;
        };
        // Check outer wire edges
        if let Some(wts) = brep.tshapes.get(fd.outer_wire.index) {
            if let TShape::Wire(wd) = wts.as_ref() {
                if wd.edges.iter().any(|er| er.index == edge_idx) {
                    result.push(fi);
                    continue;
                }
            }
        }
        // Check inner wire edges
        for iw in &fd.inner_wires {
            if let Some(wts) = brep.tshapes.get(iw.index) {
                if let TShape::Wire(wd) = wts.as_ref() {
                    if wd.edges.iter().any(|er| er.index == edge_idx) {
                        result.push(fi);
                        break;
                    }
                }
            }
        }
    }
    result
}

/// Returns the (start, end) vertex indices (tshape indices) of an edge.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::vertices_of_edge;
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// let (start, end) = vertices_of_edge(&brep, 0);
/// assert!(start < 8); // Box has 8 vertices (0-7)
/// assert!(end < 8);
/// ```
pub fn vertices_of_edge(brep: &rcad_kernel::BRep, edge_idx: usize) -> (usize, usize) {
    match brep.tshapes.get(edge_idx) {
        Some(ts) => match ts.as_ref() {
            TShape::Edge(ed) => (ed.first.index, ed.last.index),
            _ => (usize::MAX, usize::MAX),
        },
        None => (usize::MAX, usize::MAX),
    }
}

/// Returns all edge indices that reference the given vertex.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::edges_of_vertex;
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Each vertex of a box has 3 incident edges
/// let edges = edges_of_vertex(&brep, 0);
/// assert_eq!(edges.len(), 3);
/// ```
pub fn edges_of_vertex(brep: &rcad_kernel::BRep, vertex_idx: usize) -> Vec<usize> {
    brep.tshapes
        .iter()
        .enumerate()
        .filter_map(|(ei, ts)| {
            if let TShape::Edge(ed) = ts.as_ref() {
                if ed.first.index == vertex_idx || ed.last.index == vertex_idx {
                    return Some(ei);
                }
            }
            None
        })
        .collect()
}

/// Returns all faces that share the given vertex (through their edges).
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_top_adaptor::faces_of_vertex;
/// use rcad_kernel::BRep;
///
/// let brep = rcad_kernel::BRep::from_primitive(rcad_kernel::PrimitiveSolid::Box {
///     width: 1.0, height: 1.0, depth: 1.0
/// });
///
/// // Each vertex of a box is shared by 3 faces
/// let faces = faces_of_vertex(&brep, 0);
/// assert_eq!(faces.len(), 3);
/// ```
pub fn faces_of_vertex(brep: &rcad_kernel::BRep, vertex_idx: usize) -> Vec<usize> {
    let mut result = Vec::new();

    // Get all edges that reference this vertex
    let vertex_edges = edges_of_vertex(brep, vertex_idx);

    // Get all faces that reference these edges
    for edge_idx in vertex_edges {
        for face_idx in faces_of_edge(brep, edge_idx) {
            if !result.contains(&face_idx) {
                result.push(face_idx);
            }
        }
    }

    result
}

/// Returns the number of faces in a BRep.
pub fn face_count(brep: &rcad_kernel::BRep) -> usize {
    brep.tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Face(_)))
        .count()
}

/// Returns the number of shells in a BRep.
pub fn shell_count(brep: &rcad_kernel::BRep) -> usize {
    brep.tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Shell(_)))
        .count()
}

/// Returns the number of wires in a BRep (including inner wires).
pub fn wire_count(brep: &rcad_kernel::BRep) -> usize {
    brep.tshapes
        .iter()
        .filter(|ts| matches!(ts.as_ref(), TShape::Wire(_)))
        .count()
}

// =============================================================================
// Tests
// =============================================================================
