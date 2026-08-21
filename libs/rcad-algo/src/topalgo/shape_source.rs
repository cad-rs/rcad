// Shape source abstraction for the face/solid classifiers.
//
// OCCT: BRepClass / BRepClass3d (TKTopAlgo) classify points against
// TopoDS_Shape using the BRepAdaptor_* adaptors (TKBRep). They never touch
// BOPDS (TKBO). rcad's boolean pipeline stores its shapes in the DS (BOPDS)
// at flat indices, so the classifier code needs a way to read shape data
// without depending on the DS type.
//
// This trait exposes exactly the BOPDS-style shape access the classifiers use
// (NbShapes, ShapeInfo::Shape/ShapeType/SubShapes, ShapeIndex, vertex-edge
// map, BRep_Tool-style Surface/Tolerance). The bop::ds::DS implements it, so
// the dependency direction matches OCCT: topalgo (TKTopAlgo) depends only on
// rcad-kernel (TKBRep), and bop (TKBO) depends on topalgo.

use rcad_kernel::geom::{Curve2d, Surface3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{CurveRepresentation, Orientation, ShapeType, TShape};
use std::collections::HashMap;

/// OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L354-360) — the edge's pcurve
/// on the face, with the seam (CurveOnClosedSurface) pcurve selected by the
/// edge orientation in the face: FORWARD → pcurve1 (the u=2*PI image of the
/// seam), REVERSED → pcurve2 (the u=0 image). The returned range is the edge's
/// forward parameter range; the caller reverses it when the edge is REVERSED.
/// rcad keys the edge pcurve map by the DS face index, but input-shape pcurves
/// (make_cylinder etc.) are keyed by the face's preserved BRep `index` — both
/// keys are consulted.
pub fn edge_pcurve_on_face(
    ds: &dyn ShapeSource,
    edge_idx: usize,
    face_idx: usize,
    ori: Orientation,
) -> Option<(Curve2d, f64, f64)> {
    let f = ds.shape_at(face_idx);
    // OCCT BRep_Tool::CurveOnSurface (BRep_Tool.cxx L345): the pcurve key
    // location is L.Predivided(E.Location()) — the face location divided by
    // the edge's location. A located edge (prism top edge) shares its TShape
    // with the base edge but has its own pcurve key.
    let e = ds.shape_at(edge_idx);
    let key_loc = crate::bop::algo::compose_face_edge_pcurve_location(
        f.location, e.location, ds.locations());
    let face_key = (f.ptr_id(), key_loc);
    match &*e.data {
        TShape::Edge(ed) => {
            if let Some((pc1, pc2, range)) = ed.representations.iter().find_map(|r| match r {
                CurveRepresentation::CurveOnClosedSurface {
                    face,
                    pcurve1,
                    pcurve2,
                    range,
                } if *face == face_key => Some((pcurve1.clone(), pcurve2.clone(), *range)),
                _ => None,
            }) {
                let pc = if ori == Orientation::Reversed { pc2 } else { pc1 };
                return Some((pc, range[0], range[1]));
            }
            ed.pcurves.get(&face_key).cloned()
        }
        _ => None,
    }
}

/// BOPDS-style shape data source for the classifiers.
pub trait ShapeSource: Send + Sync {
    /// BOPDS_DS::NbShapes — total number of shapes in the DS.
    fn nb_shapes(&self) -> usize;
    /// BOPDS_ShapeInfo::Shape at a flat DS index.
    fn shape_at(&self, i: usize) -> Shape;
    /// BOPDS_ShapeInfo::ShapeType at a flat DS index.
    fn shape_type(&self, i: usize) -> ShapeType;
    /// BOPDS_ShapeInfo::SubShapes at a flat DS index (child shape indices).
    fn sub_shapes(&self, i: usize) -> &[usize];
    /// BOPDS_DS::ShapeIndex — map a (ptr_id, location) to its flat index.
    fn map_shape_index(&self, ptr_id: u64, location: u32) -> Option<usize>;
    /// BOPDS_DS vertex → incident edges map.
    fn map_ve(&self, vertex: usize) -> Option<&Vec<usize>>;
    /// BRep_Tool::Surface(face) — face surface at a flat DS index.
    fn face_surface(&self, i: usize) -> Option<Surface3>;
    /// BRep_Tool::Tolerance(vertex) — vertex tolerance at a flat DS index.
    fn vertex_tolerance(&self, i: usize) -> f64;
    /// BRep_Tool::Degenerated(edge) — true when the edge at a flat DS index is
    /// degenerated (no 3D curve, coincident vertices).
    fn is_edge_degenerated(&self, i: usize) -> bool;
    /// TopLoc_Location by index (0 = identity). Used by the edge-vertex
    /// traversal to compose the edge Location into its vertices
    /// (TopoDS_Iterator cumLoc semantics).
    fn get_location(&self, idx: u32) -> glam::DAffine3;
    /// The full TopLoc_Location table (index 0 = identity), for composing
    /// edge+vertex Locations (TopoDS_Iterator cumLoc semantics).
    fn locations(&self) -> &[glam::DAffine3];
}

/// A ShapeSource adapter exposing a single face and its wire edges under the
/// flat DS indexing (index 0 = the face, 1..N = the wire edges in order).
/// The synthetic draft-solid faces have no DS registration, but
/// IntTools_FClass2d::Init (OCCT IntTools_FClass2d.cxx L77-621) builds its
/// classifier from the face and the edge pcurves alone — OCCT receives a
/// TopoDS_Face, not a BOPDS DS. This adapter makes the rcad FClass2d
/// translation usable for such faces.
pub struct FaceShapeSource<'a> {
    face: &'a Shape,
    surf: Option<Surface3>,
    edges: Vec<Shape>,
    edge_index: HashMap<(u64, u32), usize>,
    locations: &'a [glam::DAffine3],
}

impl<'a> FaceShapeSource<'a> {
    /// Build the adapter: index 0 is the face, the wire edges (outer + inner,
    /// in traversal order) follow. `surf` is the face surface (already
    /// location-transformed by the caller, matching the DS convention).
    pub fn new(face: &'a Shape, surf: Surface3, locations: &'a [glam::DAffine3]) -> Self {
        let mut edges: Vec<Shape> = Vec::new();
        if let TShape::Face(fd) = &*face.data {
            for w in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                if let TShape::Wire(wd) = &*w.data {
                    for e in &wd.edges {
                        edges.push(e.clone());
                    }
                }
            }
        }
        let mut edge_index: HashMap<(u64, u32), usize> = HashMap::new();
        for (i, e) in edges.iter().enumerate() {
            edge_index.insert((e.ptr_id(), e.location), i + 1);
        }
        FaceShapeSource {
            face,
            surf: Some(surf),
            edges,
            edge_index,
            locations,
        }
    }
}

impl ShapeSource for FaceShapeSource<'_> {
    fn nb_shapes(&self) -> usize {
        1 + self.edges.len()
    }
    fn shape_at(&self, i: usize) -> Shape {
        if i == 0 {
            self.face.clone()
        } else {
            self.edges.get(i - 1).cloned().unwrap_or_else(Shape::null)
        }
    }
    fn shape_type(&self, i: usize) -> ShapeType {
        self.shape_at(i).shape_type()
    }
    fn sub_shapes(&self, _i: usize) -> &[usize] {
        &[]
    }
    fn map_shape_index(&self, ptr_id: u64, location: u32) -> Option<usize> {
        if self.face.ptr_id() == ptr_id && self.face.location == location {
            Some(0)
        } else {
            self.edge_index.get(&(ptr_id, location)).copied()
        }
    }
    fn map_ve(&self, _vertex: usize) -> Option<&Vec<usize>> {
        None
    }
    fn face_surface(&self, i: usize) -> Option<Surface3> {
        if i == 0 {
            self.surf.clone()
        } else {
            None
        }
    }
    fn vertex_tolerance(&self, _i: usize) -> f64 {
        0.0
    }
    fn is_edge_degenerated(&self, i: usize) -> bool {
        match &*self.shape_at(i).data {
            TShape::Edge(ed) => ed.degenerated,
            _ => true,
        }
    }
    fn get_location(&self, idx: u32) -> glam::DAffine3 {
        self.locations
            .get(idx as usize)
            .copied()
            .unwrap_or(glam::DAffine3::IDENTITY)
    }
    fn locations(&self) -> &[glam::DAffine3] {
        self.locations
    }
}
