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

/// Returns `true` as the serde default for the `mesh_dirty` field.
fn face_mesh_dirty_default() -> bool {
    true
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Face {
    pub outer_wire: Wire,
    pub inner_wires: Vec<Wire>,
    pub normal: DVec3,
    /// Pre-triangulated vertex index triples (into BRep.vertices).
    pub triangles: Vec<[usize; 3]>,
    /// When `true` the cached `triangles` are stale and should be recomputed
    /// before rendering.  Set to `false` by [`mesh_brep`] after tessellation,
    /// and restored to `true` by [`Face::invalidate_mesh`].
    ///
    /// This field is not serialised (transient rendering state).
    #[serde(skip, default = "face_mesh_dirty_default")]
    pub mesh_dirty: bool,
}

impl Face {
    /// Mark this face's cached mesh as stale so it will be re-tessellated on
    /// the next [`mesh_brep`] call.
    pub fn invalidate_mesh(&mut self) {
        self.mesh_dirty = true;
    }

    /// Returns `true` if the cached triangulation is up-to-date.
    pub fn mesh_is_clean(&self) -> bool {
        !self.mesh_dirty && !self.triangles.is_empty()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Shell {
    pub faces: Vec<Face>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Solid {
    pub shells: Vec<Shell>,
}

/// A connected solid made of multiple adjacent solids that share boundary faces.
///
/// Analogous to OCCT `TopoDS_CompSolid`. All contained solids must form a
/// topologically connected manifold body. CompSolid allows expressing structures
/// like multi-region models (e.g. a solid that is split into sub-regions by
/// internal surfaces) without performing a full boolean union.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CompSolid {
    /// The constituent connected solids.
    pub solids: Vec<Solid>,
    /// Optional label for this CompSolid (e.g. from an assembly tree).
    #[serde(default)]
    pub label: Option<String>,
}

impl CompSolid {
    /// Create an empty CompSolid.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a CompSolid from a list of solids.
    pub fn from_solids(solids: Vec<Solid>) -> Self {
        Self { solids, label: None }
    }

    /// Total number of faces across all constituent solids.
    pub fn face_count(&self) -> usize {
        self.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    }
}

/// A heterogeneous collection of shapes (solids, shells, wires, etc.).
///
/// Analogous to OCCT `TopoDS_Compound`. A Compound can hold any mix of:
/// - Complete solids (`BRep`)
/// - Connected solid groups (`CompSolid`)
/// - Free shells
/// - Free wires / edges
///
/// Compounds are the top-level shape type for assemblies and imported STEP
/// files that contain multiple disconnected bodies.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Compound {
    /// Named sub-shapes (solids).
    ///
    /// Each entry is `(label, shape)`. The label is optional and may be
    /// empty — it is used for assembly-tree bookkeeping and STEP name mapping.
    pub solids: Vec<(Option<String>, Solid)>,
    /// Named CompSolids (multi-region connected solid groups).
    pub comp_solids: Vec<(Option<String>, CompSolid)>,
    /// Loose shells not attached to any solid.
    pub shells: Vec<(Option<String>, Shell)>,
    /// Nested sub-compounds (for deeply hierarchical assemblies).
    pub compounds: Vec<(Option<String>, Compound)>,
}

impl Compound {
    /// Create an empty Compound.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a solid with an optional label.
    pub fn add_solid(&mut self, label: Option<String>, solid: Solid) {
        self.solids.push((label, solid));
    }

    /// Add a CompSolid with an optional label.
    pub fn add_comp_solid(&mut self, label: Option<String>, comp_solid: CompSolid) {
        self.comp_solids.push((label, comp_solid));
    }

    /// Flatten all constituent solids into a single list (discards compound hierarchy).
    pub fn flatten_solids(&self) -> Vec<&Solid> {
        let mut out = Vec::new();
        for (_, s) in &self.solids {
            out.push(s);
        }
        for (_, cs) in &self.comp_solids {
            for s in &cs.solids {
                out.push(s);
            }
        }
        for (_, sub) in &self.compounds {
            out.extend(sub.flatten_solids());
        }
        out
    }

    /// Total face count across all constituent shapes.
    pub fn face_count(&self) -> usize {
        self.flatten_solids()
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    }
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
            mesh_dirty: true,
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
            mesh_dirty: true,
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
                    mesh_dirty: true,
                },
                Face {
                    outer_wire: Wire { edges: vec![] },
                    inner_wires: vec![],
                    normal: DVec3::NEG_X,
                    triangles: vec![],
                    mesh_dirty: true,
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
