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
    let face_key = (f.ptr_id(), f.location);
    match &*ds.shape_at(edge_idx).data {
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
}
