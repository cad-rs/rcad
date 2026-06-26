use std::sync::Arc;
use std::collections::HashMap;
use glam::DVec3;
use serde::{Deserialize, Serialize};
use crate::geom::{Curve3, Surface3};

/// OCCT TopAbs_Orientation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Orientation {
    Forward,
    Reversed,
    Internal,
    External,
}

impl Orientation {
    pub const fn is_forward(self) -> bool {
        matches!(self, Orientation::Forward)
    }
}

/// OCCT TopAbs_ShapeEnum
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ShapeType {
    Vertex,
    Edge,
    Wire,
    Face,
    Shell,
    Solid,
    CompSolid,
    Compound,
}

/// TopoDS_Shape equivalent: index into the TShape pool + Orientation.
/// This is a value type (Copy, small) — the geometric data lives in the TShape Arc.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ShapeRef {
    /// Index into BRep.tshapes[]
    pub index: usize,
    pub orientation: Orientation,
}

impl ShapeRef {
    pub const fn new(index: usize) -> Self {
        Self { index, orientation: Orientation::Forward }
    }
    pub const fn with_orientation(index: usize, orientation: Orientation) -> Self {
        Self { index, orientation }
    }
}

/// TShape — shared geometric/topological data (analogous to TopoDS_TShape + subclasses).
/// Stored in Arc<TShape> within BRep so multiple ShapeRefs share the same data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TShape {
    Vertex(TVertexData),
    Edge(TEdgeData),
    Wire(TWireData),
    Face(TFaceData),
    Shell(TShellData),
    Solid(TSolidData),
    CompSolid(Vec<ShapeRef>),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TVertexData {
    pub point: DVec3,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TEdgeData {
    pub curve: Option<usize>,
    pub first: ShapeRef,
    pub last: ShapeRef,
    pub range: [f64; 2],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TWireData {
    pub edges: Vec<ShapeRef>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TFaceData {
    pub surface: Option<usize>,
    pub outer_wire: ShapeRef,
    pub inner_wires: Vec<ShapeRef>,
    pub sample_point: Option<DVec3>,
    /// UV domain [umin, umax, vmin, vmax] — used by surface area calculation.
    pub uv_domain: Option<[f64; 4]>,
    /// INTERNAL vertices (OCCT: TopAbs_INTERNAL sub-shapes, BRep_Builder.Add(aF, aV)).
    pub internal_vertices: Vec<ShapeRef>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TShellData {
    pub faces: Vec<ShapeRef>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TSolidData {
    pub shells: Vec<ShapeRef>,
}

/// BRep top-level shape container — all TShapes in a single pool with shared Arc ownership.
/// Analogous to OCCT's Doc/assembly structure where all TShapes live in a shared scope.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BRep {
    pub tshapes: Vec<Arc<TShape>>,
    pub curves: Vec<Curve3>,
    pub surfaces: Vec<Surface3>,
    pub curve2ds: Vec<crate::geom::Curve2d>,
}

impl Default for BRep {
    fn default() -> Self { Self::new() }
}

impl BRep {
    pub fn new() -> Self {
        Self { tshapes: Vec::new(), curves: Vec::new(), surfaces: Vec::new(), curve2ds: Vec::new() }
    }

    pub fn add_tvertex(&mut self, point: DVec3) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Vertex(TVertexData { point })));
        ShapeRef::new(index)
    }

    pub fn add_tedge(&mut self, curve: Option<usize>, first: ShapeRef, last: ShapeRef, range: [f64; 2]) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Edge(TEdgeData { curve, first, last, range })));
        ShapeRef::new(index)
    }

    pub fn add_twire(&mut self, edges: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Wire(TWireData { edges })));
        ShapeRef::new(index)
    }

    pub fn add_tface(&mut self, surface: Option<usize>, outer_wire: ShapeRef, inner_wires: Vec<ShapeRef>, sample_point: Option<DVec3>, uv_domain: Option<[f64; 4]>, internal_vertices: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Face(TFaceData { surface, outer_wire, inner_wires, sample_point, uv_domain, internal_vertices })));
        ShapeRef::new(index)
    }

    pub fn add_tshell(&mut self, faces: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Shell(TShellData { faces })));
        ShapeRef::new(index)
    }

    pub fn add_tsolid(&mut self, shells: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::Solid(TSolidData { shells })));
        ShapeRef::new(index)
    }

    pub fn add_tcompsolid(&mut self, solids: Vec<ShapeRef>) -> ShapeRef {
        let index = self.tshapes.len();
        self.tshapes.push(Arc::new(TShape::CompSolid(solids)));
        ShapeRef::new(index)
    }

    pub fn vertex(&self, r: ShapeRef) -> &TVertexData {
        match &*self.tshapes[r.index] {
            TShape::Vertex(v) => v,
            _ => panic!("ShapeRef {} is not a Vertex", r.index),
        }
    }

    pub fn edge(&self, r: ShapeRef) -> &TEdgeData {
        match &*self.tshapes[r.index] {
            TShape::Edge(e) => e,
            _ => panic!("ShapeRef {} is not an Edge", r.index),
        }
    }

    pub fn wire(&self, r: ShapeRef) -> &TWireData {
        match &*self.tshapes[r.index] {
            TShape::Wire(w) => w,
            _ => panic!("ShapeRef {} is not a Wire", r.index),
        }
    }

    pub fn face(&self, r: ShapeRef) -> &TFaceData {
        match &*self.tshapes[r.index] {
            TShape::Face(f) => f,
            _ => panic!("ShapeRef {} is not a Face", r.index),
        }
    }

    pub fn shell(&self, r: ShapeRef) -> &TShellData {
        match &*self.tshapes[r.index] {
            TShape::Shell(s) => s,
            _ => panic!("ShapeRef {} is not a Shell", r.index),
        }
    }

    pub fn solid(&self, r: ShapeRef) -> &TSolidData {
        match &*self.tshapes[r.index] {
            TShape::Solid(s) => s,
            _ => panic!("ShapeRef {} is not a Solid", r.index),
        }
    }

    /// Build a minimal Cube (for testing).
    pub fn build_unit_cube() -> (Self, ShapeRef) {
        let mut brep = BRep::new();
        let mut bld = BRepBuilder::new();
        bld.build_unit_cube(&mut brep);
        let root = bld.root.expect("cube should have root solid");
        (brep, root)
    }
}

// ---------------------------------------------------------------------------
// BRepBuilder — OCCT BRep_Builder equivalent for incrementally constructing BRep
// ---------------------------------------------------------------------------

pub struct BRepBuilder {
    pub root: Option<ShapeRef>,
    vertex_cache: Vec<[f64; 3]>,
}

impl BRepBuilder {
    pub fn new() -> Self {
        Self { root: None, vertex_cache: Vec::new() }
    }

    fn find_or_add_vertex(&mut self, brep: &mut BRep, pt: DVec3) -> ShapeRef {
        for (i, &cached) in self.vertex_cache.iter().enumerate() {
            let dp = DVec3::new(cached[0], cached[1], cached[2]) - pt;
            if dp.length_squared() < 1e-30 {
                return ShapeRef::new(i);
            }
        }
        let r = brep.add_tvertex(pt);
        self.vertex_cache.push([pt.x, pt.y, pt.z]);
        r
    }

    pub fn build_unit_cube(&mut self, brep: &mut BRep) {
        // 8 vertices
        let v: Vec<ShapeRef> = [
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
        ].iter().map(|&p| self.find_or_add_vertex(brep, DVec3::new(p[0], p[1], p[2]))).collect();

        // 12 edges (unit box: 4 around each of 3 axis-aligned rings)
        let faces_verts: [[usize; 4]; 6] = [
            [0, 1, 2, 3], // Z- (bottom)
            [4, 5, 6, 7], // Z+ (top)
            [0, 1, 5, 4], // Y- (front)
            [2, 3, 7, 6], // Y+ (back)
            [0, 3, 7, 4], // X- (left)
            [1, 2, 6, 5], // X+ (right)
        ];

        // Build 12 edges, reuse shared vertex references between faces
        let mut edge_map = std::collections::HashMap::new();
        let mut edge_for_face = Vec::new();

        for &fv in &faces_verts {
            let mut face_edges = Vec::new();
            for j in 0..4 {
                let a = fv[j];
                let b = fv[(j + 1) % 4];
                let key = if a < b { (a, b, true) } else { (b, a, false) };
                let e_idx = *edge_map.entry((a.min(b), a.max(b))).or_insert_with(|| {
                    let va = v[a];
                    let vb = v[b];
                    brep.add_tedge(None, va, ShapeRef::with_orientation(vb.index, Orientation::Reversed), [0.0, 1.0]).index
                });
                let orient = if a < b { Orientation::Forward } else { Orientation::Reversed };
                face_edges.push(ShapeRef::with_orientation(e_idx, orient));
            }
            edge_for_face.push(face_edges);
        }

        // Build 6 faces, collecting their refs for shell building
        let mut face_refs = Vec::new();
        for (i, face_edges) in edge_for_face.into_iter().enumerate() {
            let wire = brep.add_twire(face_edges);
            let face = brep.add_tface(None, wire, vec![], None, None, vec![]);
            let face_ref = ShapeRef::new(face.index);

            // Outer normal: orient face based on the face index
            // For a unit cube at origin, faces 0,2,4 have inward normals, orient REVERSED
            // faces 1,3,5 have outward normals, orient FORWARD
            let orient = match i {
                0 => Orientation::Reversed,  // Z-
                1 => Orientation::Forward,   // Z+
                2 => Orientation::Reversed,  // Y-
                3 => Orientation::Forward,   // Y+
                4 => Orientation::Reversed,  // X-
                5 => Orientation::Forward,   // X+
                _ => unreachable!(),
            };
            face_refs.push(ShapeRef::with_orientation(face.index, orient));
        }

        // Build shell from the 6 oriented faces
        let shell = brep.add_tshell(face_refs);
        // Build solid from the shell
        let solid = brep.add_tsolid(vec![shell]);
        self.root = Some(solid);
    }
}

// ---------------------------------------------------------------------------
// ShapeType helpers
// ---------------------------------------------------------------------------

impl TShape {
    pub fn shape_type(&self) -> ShapeType {
        match self {
            TShape::Vertex(_) => ShapeType::Vertex,
            TShape::Edge(_) => ShapeType::Edge,
            TShape::Wire(_) => ShapeType::Wire,
            TShape::Face(_) => ShapeType::Face,
            TShape::Shell(_) => ShapeType::Shell,
            TShape::Solid(_) => ShapeType::Solid,
            TShape::CompSolid(_) => ShapeType::CompSolid,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orientation_values() {
        assert!(Orientation::Forward.is_forward());
        assert!(!Orientation::Reversed.is_forward());
        assert!(!Orientation::Internal.is_forward());
        assert!(!Orientation::External.is_forward());
    }

    #[test]
    fn test_shape_ref_construction() {
        let r = ShapeRef::new(5);
        assert_eq!(r.index, 5);
        assert_eq!(r.orientation, Orientation::Forward);

        let r2 = ShapeRef::with_orientation(3, Orientation::Reversed);
        assert_eq!(r2.index, 3);
        assert_eq!(r2.orientation, Orientation::Reversed);
    }

    #[test]
    fn test_empty_brep() {
        let brep = BRep::new();
        assert!(brep.tshapes.is_empty());
    }

    #[test]
    fn test_vertex_creation() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(brep.tshapes.len(), 1);
        assert_eq!(v.index, 0);
        assert_eq!(brep.vertex(v).point, DVec3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn test_vertex_sharing() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::ZERO);
        // Two separate calls create two distinct TShapes (OCCT: each call adds a unique TShape)
        assert_eq!(brep.tshapes.len(), 2);
        assert_eq!(brep.vertex(v0).point, brep.vertex(v1).point);
        // Same Arc — points are equal but TShape identity differs (OCCT behavior)
        assert!(!Arc::ptr_eq(&brep.tshapes[v0.index], &brep.tshapes[v1.index]));
    }

    #[test]
    fn test_edge_creation() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let e = brep.add_tedge(None, v0, v1, [0.0, 1.0]);
        assert_eq!(brep.tshapes.len(), 3);
        let ed = brep.edge(e);
        assert_eq!(ed.first.index, v0.index);
        assert_eq!(ed.last.index, v1.index);
    }

    #[test]
    fn test_edge_shares_vertex_tshape() {
        // Two edges sharing the same vertex TShape
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        let v2 = brep.add_tvertex(DVec3::new(1.0, 1.0, 0.0));

        let e0 = brep.add_tedge(None, v0, v1, [0.0, 1.0]);
        let e1 = brep.add_tedge(None, v1, v2, [0.0, 1.0]);

        // Both edges reference v1 at index 1
        assert_eq!(brep.edge(e0).last.index, v1.index);
        assert_eq!(brep.edge(e1).first.index, v1.index);
        // Same TShape identity (v1)
        assert!(Arc::ptr_eq(&brep.tshapes[brep.edge(e0).last.index], &brep.tshapes[brep.edge(e1).first.index]));
    }

    #[test]
    fn test_wire_and_face() {
        let mut brep = BRep::new();
        let v = (0..4).map(|i| brep.add_tvertex(DVec3::new(i as f64, 0.0, 0.0))).collect::<Vec<_>>();
        let e0 = brep.add_tedge(None, v[0], v[1], [0.0, 1.0]);
        let e1 = brep.add_tedge(None, v[1], v[2], [0.0, 1.0]);
        let e2 = brep.add_tedge(None, v[2], v[3], [0.0, 1.0]);
        let wire = brep.add_twire(vec![e0, e1, e2]);
        let face = brep.add_tface(None, wire, vec![], None, None, vec![]);
        let fd = brep.face(face);
        assert_eq!(brep.tshapes.len(), 9); // 4V + 3E + 1W + 1F
        assert!(fd.inner_wires.is_empty());
        assert!(fd.sample_point.is_none());
    }

    #[test]
    fn test_shape_type_discrimination() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::ZERO);
        assert_eq!(brep.tshapes[v.index].shape_type(), ShapeType::Vertex);

        let e = brep.add_tedge(None, v, v, [0.0, 1.0]);
        assert_eq!(brep.tshapes[e.index].shape_type(), ShapeType::Edge);

        let w = brep.add_twire(vec![e]);
        assert_eq!(brep.tshapes[w.index].shape_type(), ShapeType::Wire);

        let f = brep.add_tface(None, w, vec![], None, None, vec![]);
        assert_eq!(brep.tshapes[f.index].shape_type(), ShapeType::Face);

        let sh = brep.add_tshell(vec![f]);
        assert_eq!(brep.tshapes[sh.index].shape_type(), ShapeType::Shell);

        let so = brep.add_tsolid(vec![sh]);
        assert_eq!(brep.tshapes[so.index].shape_type(), ShapeType::Solid);
    }

    #[test]
    fn test_orientation_on_shape_reference() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let v1 = brep.add_tvertex(DVec3::X);
        // Edge references v1 with REVERSED orientation (last vertex traversed backward)
        let e = brep.add_tedge(None, v0, ShapeRef::with_orientation(v1.index, Orientation::Reversed), [0.0, 1.0]);
        let ed = brep.edge(e);
        assert_eq!(ed.last.orientation, Orientation::Reversed);
    }

    #[test]
    fn test_clone_preserves_shape_identity() {
        let mut brep = BRep::new();
        let v0 = brep.add_tvertex(DVec3::ZERO);
        let _v1 = brep.add_tvertex(DVec3::X);
        let e = brep.add_tedge(None, v0, v0, [0.0, 1.0]);

        let cloned = brep.clone();
        // Same TShape identity in clone (Arc::ptr_eq across clone)
        assert_eq!(cloned.tshapes.len(), brep.tshapes.len());
        assert!(Arc::ptr_eq(&cloned.tshapes[e.index], &brep.tshapes[e.index]));
    }

    #[test]
    fn test_serialize_roundtrip() {
        let mut brep = BRep::new();
        let v = brep.add_tvertex(DVec3::new(1.0, 2.0, 3.0));
        let e = brep.add_tedge(None, v, v, [0.0, 1.0]);

        let json = serde_json::to_string(&brep).unwrap();
        let restored: BRep = serde_json::from_str(&json).unwrap();
        assert_eq!(restored.tshapes.len(), 2);
        assert_eq!(restored.vertex(ShapeRef::new(0)).point, DVec3::new(1.0, 2.0, 3.0));
        assert_eq!(restored.edge(ShapeRef::new(1)).range, [0.0, 1.0]);
    }

    #[test]
    fn test_brep_builder_unit_cube() {
        let (_brep, root) = BRep::build_unit_cube();
        assert_eq!(root.orientation, Orientation::Forward);
    }

    #[test]
    fn test_wire_edge_orientation() {
        let mut brep = BRep::new();
        let v = (0..4).map(|i| brep.add_tvertex(DVec3::new(i as f64, 0.0, 0.0))).collect::<Vec<_>>();
        let e0 = brep.add_tedge(None, v[0], v[1], [0.0, 1.0]);

        // Wire with forward edge
        let wire_fwd = brep.add_twire(vec![e0]);
        assert_eq!(brep.tshapes[wire_fwd.index].shape_type(), ShapeType::Wire);

        // Wire with reversed edge
        let e0_rev = ShapeRef::with_orientation(e0.index, Orientation::Reversed);
        let wire_rev = brep.add_twire(vec![e0_rev]);
        if let TShape::Wire(ref wd) = *brep.tshapes[wire_rev.index] {
            assert_eq!(wd.edges[0].orientation, Orientation::Reversed);
        } else {
            panic!("expected Wire");
        }

        // Same TShape for the edge, different orientation on the reference
        assert!(Arc::ptr_eq(&brep.tshapes[e0.index], &brep.tshapes[e0_rev.index]));
    }

    #[test]
    fn test_to_topods_roundtrip_sa() {
        // Build a simple box BRep, convert to topods and back, check SA
        let mut orig = crate::BRep::new();
        // 8 vertices of a unit cube
        let v: Vec<_> = [
            [0.0,0.0,0.0],[1.0,0.0,0.0],[1.0,1.0,0.0],[0.0,1.0,0.0],
            [0.0,0.0,1.0],[1.0,0.0,1.0],[1.0,1.0,1.0],[0.0,1.0,1.0],
        ].iter().map(|&p| orig.vertices.push(crate::topology::Vertex { point: DVec3::new(p[0],p[1],p[2]) })).collect::<Vec<_>>();
        // Simple box with 6 faces, 12 edges
        let face_vert_idxs: [[usize;4]; 6] = [
            [0,1,2,3],[4,5,6,7],[0,1,5,4],[2,3,7,6],[0,3,7,4],[1,2,6,5],
        ];
        for &fv in &face_vert_idxs {
            for j in 0..4 {
                let a = fv[j]; let b = fv[(j+1)%4];
                // Check if edge already exists
                let exists = orig.edges.iter().position(|e| (e.start==a&&e.end==b)||(e.start==b&&e.end==a));
                if exists.is_none() {
                    orig.edges.push(crate::topology::Edge { start: a, end: b });
                }
            }
        }
        // Convert to topods and back
        let t = orig.to_topods();
        let back = crate::BRep::from_topods(&t);
        assert_eq!(back.vertices.len(), 8, "vertex count");
        assert_eq!(back.edges.len(), 12, "edge count");
        // Check surface area (unit cube SA = 6.0)
        let sa = crate::surface_area(&back);
        assert!((sa - 6.0).abs() < 0.01, "unit cube SA: expected 6.0, got {}", sa);
    }
}
