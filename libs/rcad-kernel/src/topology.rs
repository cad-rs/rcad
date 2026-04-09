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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wire_edge_fwd_rev() {
        let fwd = WireEdge::fwd(3);
        assert_eq!(fwd.idx, 3);
        assert!(fwd.forward);

        let rev = WireEdge::rev(5);
        assert_eq!(rev.idx, 5);
        assert!(!rev.forward);
    }

    #[test]
    fn wire_contains_edges() {
        let w = Wire {
            edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::rev(2)],
        };
        assert_eq!(w.edges.len(), 3);
        assert!(!w.edges[2].forward);
    }

    #[test]
    fn face_has_outer_wire_and_no_inner_wires_by_default() {
        let f = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
        };
        assert!(f.inner_wires.is_empty());
        assert_eq!(f.normal, DVec3::Z);
    }

    #[test]
    fn face_with_inner_wire() {
        let f = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![Wire {
                edges: vec![WireEdge::fwd(3), WireEdge::fwd(4)],
            }],
            normal: DVec3::Y,
            triangles: vec![],
        };
        assert_eq!(f.inner_wires.len(), 1);
        assert_eq!(f.inner_wires[0].edges.len(), 2);
    }

    #[test]
    fn shell_contains_faces() {
        let shell = Shell {
            faces: vec![
                Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::X,
                    triangles: vec![],
                },
                Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::NEG_X,
                    triangles: vec![],
                },
            ],
        };
        assert_eq!(shell.faces.len(), 2);
    }

    #[test]
    fn solid_contains_shells() {
        let solid = Solid {
            shells: vec![Shell { faces: vec![] }],
        };
        assert_eq!(solid.shells.len(), 1);
    }
}
