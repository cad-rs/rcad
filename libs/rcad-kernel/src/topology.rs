use glam::DVec3;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Vertex {
    pub point: DVec3,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Edge {
    pub start: usize,
    pub end: usize,
}

/// An edge reference with explicit traversal direction inside a Wire.
///
/// `forward = true`  → traverse edge from `edge.start` to `edge.end`.
/// `forward = false` → traverse edge from `edge.end`   to `edge.start`.
///
/// Analogous to OCCT `TopoDS_Edge` with `FORWARD` / `REVERSED` orientation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WireEdge {
    /// Index into `BRep.edges`.
    pub idx: usize,
    /// Traversal direction: `true` = forward (start→end), `false` = reversed.
    pub forward: bool,
}

impl WireEdge {
    pub const fn new(idx: usize, forward: bool) -> Self {
        Self { idx, forward }
    }
    /// Shorthand: forward reference.
    pub const fn fwd(idx: usize) -> Self {
        Self { idx, forward: true }
    }
    /// Shorthand: reversed reference.
    pub const fn rev(idx: usize) -> Self {
        Self {
            idx,
            forward: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Wire {
    pub edges: Vec<WireEdge>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face {
    pub outer_wire: Wire,
    pub inner_wires: Vec<Wire>,
    pub normal: DVec3,
    /// Pre-triangulated vertex index triples (into BRep.vertices)
    pub triangles: Vec<[usize; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    pub faces: Vec<Face>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Solid {
    pub shells: Vec<Shell>,
}
