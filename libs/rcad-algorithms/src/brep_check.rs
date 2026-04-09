//! BRep validity checker.
//!
//! Analogous to OCCT `BRepCheck_Analyzer`. Checks structural and geometric
//! consistency of a BRep without modifying it.
//!
//! # Checks performed
//!
//! - **C1 Wire closure**: every wire must form a closed chain — the end vertex of
//!   each edge must equal the start vertex of the next edge.
//! - **C2 Face normal consistency**: each face's stored normal must not be a zero
//!   vector.
//! - **C3 Degenerate face**: faces with fewer than 3 wire edges are degenerate.
//! - **C4 Edge index validity**: WireEdge indices must be within bounds of
//!   `brep.edges`.
//! - **C5 Vertex index validity**: each edge's start/end indices must be within
//!   bounds of `brep.vertices`.

use glam::DVec3;
use rcad_kernel::BRep;

/// A single validity issue found during checking.
#[derive(Debug, Clone, PartialEq)]
pub enum CheckIssue {
    /// Wire is not closed: end vertex of edge `edge_idx` does not match start
    /// vertex of the next edge in the wire (solid `solid`, shell `shell`,
    /// face `face`, position `wire_pos`).
    OpenWire {
        solid: usize,
        shell: usize,
        face: usize,
        /// Index of the edge within the wire where the gap occurs.
        wire_pos: usize,
    },
    /// Face normal is a zero vector.
    ZeroNormal {
        solid: usize,
        shell: usize,
        face: usize,
    },
    /// Face outer wire has fewer than 3 edges.
    DegenerateFace {
        solid: usize,
        shell: usize,
        face: usize,
    },
    /// A WireEdge references an edge index that is out of bounds.
    InvalidEdgeIndex {
        solid: usize,
        shell: usize,
        face: usize,
        edge_idx: usize,
    },
    /// An edge references a vertex index that is out of bounds.
    InvalidVertexIndex { edge: usize, vertex_idx: usize },
}

impl std::fmt::Display for CheckIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CheckIssue::OpenWire {
                solid,
                shell,
                face,
                wire_pos,
            } => write!(
                f,
                "OpenWire: solid={solid} shell={shell} face={face} at wire pos {wire_pos}"
            ),
            CheckIssue::ZeroNormal { solid, shell, face } => {
                write!(f, "ZeroNormal: solid={solid} shell={shell} face={face}")
            }
            CheckIssue::DegenerateFace { solid, shell, face } => {
                write!(f, "DegenerateFace: solid={solid} shell={shell} face={face}")
            }
            CheckIssue::InvalidEdgeIndex {
                solid,
                shell,
                face,
                edge_idx,
            } => write!(
                f,
                "InvalidEdgeIndex: solid={solid} shell={shell} face={face} edge={edge_idx}"
            ),
            CheckIssue::InvalidVertexIndex { edge, vertex_idx } => {
                write!(f, "InvalidVertexIndex: edge={edge} vertex={vertex_idx}")
            }
        }
    }
}

/// Result of a BRep validity check.
#[derive(Debug, Clone)]
pub struct CheckResult {
    pub issues: Vec<CheckIssue>,
}

impl CheckResult {
    /// Returns `true` if no issues were found.
    pub fn is_valid(&self) -> bool {
        self.issues.is_empty()
    }
}

/// Check the validity of a BRep and return a `CheckResult` with any issues found.
///
/// Analogous to OCCT `BRepCheck_Analyzer::Perform()`.
pub fn check(brep: &BRep) -> CheckResult {
    let mut issues = Vec::new();
    let n_edges = brep.edges.len();
    let n_verts = brep.vertices.len();

    // C5: edge vertex bounds
    for (eidx, edge) in brep.edges.iter().enumerate() {
        if edge.start >= n_verts {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: edge.start,
            });
        }
        if edge.end >= n_verts {
            issues.push(CheckIssue::InvalidVertexIndex {
                edge: eidx,
                vertex_idx: edge.end,
            });
        }
    }

    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                let wire = &face.outer_wire;

                // C2: zero normal
                if face.normal == DVec3::ZERO {
                    issues.push(CheckIssue::ZeroNormal {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                }

                // C3: degenerate face
                if wire.edges.len() < 3 {
                    issues.push(CheckIssue::DegenerateFace {
                        solid: si,
                        shell: shi,
                        face: fi,
                    });
                    // Can't check wire closure for degenerate face
                    continue;
                }

                // C4: edge index bounds + collect start/end vertices for wire closure check
                let mut valid = true;
                let mut wire_verts: Vec<(usize, usize)> = Vec::new(); // (start_vidx, end_vidx)
                for we in &wire.edges {
                    if we.idx >= n_edges {
                        issues.push(CheckIssue::InvalidEdgeIndex {
                            solid: si,
                            shell: shi,
                            face: fi,
                            edge_idx: we.idx,
                        });
                        valid = false;
                    } else {
                        let edge = &brep.edges[we.idx];
                        let (sv, ev) = if we.forward {
                            (edge.start, edge.end)
                        } else {
                            (edge.end, edge.start)
                        };
                        wire_verts.push((sv, ev));
                    }
                }

                if !valid {
                    continue;
                }

                // C1: wire closure — end of edge[i] must match start of edge[i+1]
                let n = wire_verts.len();
                for i in 0..n {
                    let next = (i + 1) % n;
                    let end_v = wire_verts[i].1;
                    let start_v = wire_verts[next].0;
                    if end_v != start_v {
                        // Tolerance check: allow same position even if different vertex objects
                        let end_pt = brep.vertices[end_v].point;
                        let start_pt = brep.vertices[start_v].point;
                        if (end_pt - start_pt).length() > 1e-6 {
                            issues.push(CheckIssue::OpenWire {
                                solid: si,
                                shell: shi,
                                face: fi,
                                wire_pos: i,
                            });
                        }
                    }
                }
            }
        }
    }

    CheckResult { issues }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn unit_box_is_valid() {
        let brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        let result = check(&brep);
        assert!(
            result.is_valid(),
            "unit box should pass all checks; issues: {:?}",
            result.issues
        );
    }

    #[test]
    fn open_wire_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        // Build a BRep with a deliberately open wire (gap between edge 1 end and edge 0 start)
        let mut brep = BRep::new();
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 0.0, 0.0),
        }); // 0
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 0.0, 0.0),
        }); // 1
        brep.vertices.push(Vertex {
            point: DVec3::new(1.0, 1.0, 0.0),
        }); // 2
        brep.vertices.push(Vertex {
            point: DVec3::new(0.0, 1.0, 0.0),
        }); // 3 (gap: wire goes 0→1→2 then 2→0 skips 3)

        // Edge 0: v0 → v1; Edge 1: v1 → v2; Edge 2: v2 → v0 (skips v3 — would close)
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 3, end: 0 }); // intentional gap: starts at v3 not v2

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(1),
                    WireEdge::fwd(2), // e2 starts at v3, but e1 ends at v2 → open
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(!result.is_valid(), "open wire should be detected");
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::OpenWire { .. }))
        );
    }

    #[test]
    fn degenerate_face_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 0 });

        // Face with only 2 edges — degenerate
        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1)],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::DegenerateFace { .. }))
        );
    }

    #[test]
    fn zero_normal_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        for p in [DVec3::ZERO, DVec3::X, DVec3::Y, DVec3::Z] {
            brep.vertices.push(Vertex { point: p });
        }
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.edges.push(Edge { start: 1, end: 2 });
        brep.edges.push(Edge { start: 2, end: 0 });

        let face = Face {
            outer_wire: Wire {
                edges: vec![WireEdge::fwd(0), WireEdge::fwd(1), WireEdge::fwd(2)],
            },
            inner_wires: vec![],
            normal: DVec3::ZERO, // zero normal — invalid
            triangles: vec![],
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::ZeroNormal { .. })),
            "expected ZeroNormal issue"
        );
    }

    #[test]
    fn invalid_edge_index_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.vertices.push(Vertex { point: DVec3::X });
        brep.vertices.push(Vertex { point: DVec3::Y });
        brep.edges.push(Edge { start: 0, end: 1 }); // only edge 0 exists

        let face = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::fwd(0),
                    WireEdge::fwd(99), // out-of-bounds
                    WireEdge::fwd(0),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
        };
        brep.solids.push(Solid {
            shells: vec![Shell { faces: vec![face] }],
        });

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::InvalidEdgeIndex { .. })),
            "expected InvalidEdgeIndex issue"
        );
    }

    #[test]
    fn invalid_vertex_index_is_detected() {
        use glam::DVec3;
        use rcad_kernel::topology::{Edge, Vertex};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: DVec3::ZERO });
        brep.edges.push(Edge { start: 0, end: 99 }); // vertex 99 doesn't exist

        let result = check(&brep);
        assert!(
            result
                .issues
                .iter()
                .any(|i| matches!(i, CheckIssue::InvalidVertexIndex { .. })),
            "expected InvalidVertexIndex issue"
        );
    }
}
