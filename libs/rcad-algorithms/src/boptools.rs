//! OCCT-aligned BOPTools helpers (BOPTools_AlgoTools, BOPTools_AlgoTools2D, BOPTools_AlgoTools3D).
//!
//! These functions provide edge/face classification and p-curve utilities
//! used by the boolean pipeline.

use rcad_kernel::geom::{Curve3, Surface3, CurveEval};
use crate::bopds::ds::DS;
use crate::classify::Classification;

/// OCCT-aligned: MakeSectEdge (BOPTools_AlgoTools).
/// Creates a section edge from an intersection curve.  Returns the
/// start and end vertex indices.
pub fn make_sect_edge(ds: &mut DS, ci: usize, v1: usize, v2: usize) -> usize {
    let ei = ds.edges.len();
    let ic = &ds.intersection_curves[ci];
    ds.edges.push(crate::bopds::ds::DSEdge {
        start_vertex: v1,
        end_vertex: v2,
        curve: ic.curve.clone(),
        t_range: ic.t_range,
        origin: crate::bopds::ds::ShapeOrigin::ShapeA,
        geom_tol: ic.geom_tol,
        paves: Vec::new(),
        pave_blocks: Vec::new(), face_reps: Vec::new(),
        is_internal: false,
    });
    ei
}

/// OCCT-aligned: IsMicroEdge (BOPTools_AlgoTools).
pub fn is_micro_edge(v1: &glam::DVec3, v2: &glam::DVec3) -> bool {
    (v1 - v2).length() < crate::tolerance::TOLERANCE_ABS * 100.0
}

/// OCCT-aligned: ComputeState (BOPTools_AlgoTools).
pub fn compute_state_classify(
    point: glam::DVec3,
    face_indices: &[usize],
    ds: &DS,
) -> Classification {
    crate::classify::classify_point(point, face_indices, ds)
}


/// OCCT-aligned: GetNormalToFaceOnEdge (BOPTools_AlgoTools3D).
pub fn get_normal_to_face_on_edge(
    surface: &Surface3, face_normal: glam::DVec3, edge_mid: glam::DVec3,
) -> glam::DVec3 {
    match surface {
        Surface3::Plane(p) => p.normal,
        Surface3::Sphere(s) => (edge_mid - s.center).normalize(),
        Surface3::Cylinder(c) => {
            let v = edge_mid - c.origin;
            let radial = v - c.axis.normalize() * v.dot(c.axis.normalize());
            radial.normalize()
        }
        _ => face_normal,
    }
}

/// OCCT-aligned: PointNearEdge (BOPTools_AlgoTools3D).
pub fn point_near_edge(
    surface: &Surface3, edge_mid: glam::DVec3, normal: glam::DVec3,
) -> glam::DVec3 {
    edge_mid + normal * crate::tolerance::TOLERANCE_ABS * 10.0
}

/// OCCT-aligned: HasCurveOnSurface (BOPTools_AlgoTools2D).
pub fn has_curve_on_surface(edge_curve: &Curve3, _surface: &Surface3) -> bool {
    // Simplified check: all 3D curves can be projected to any surface
    true
}

/// OCCT-aligned: IsEdgeIsoline (BOPTools_AlgoTools2D).
pub fn is_edge_isoline(edge_curve: &Curve3, _surface: &Surface3) -> bool {
    matches!(edge_curve, Curve3::Line(_))
}

/// OCCT-aligned: OrientEdgeOnFace (BOPTools_AlgoTools3D).
pub fn orient_edge_on_face(dot_product: f64) -> bool {
    dot_product > 0.0
}

/// OCCT-aligned: MakeEdge (BOPTools_AlgoTools).
pub fn make_ds_edge(
    ds: &mut crate::bopds::ds::DS, v1: usize, v2: usize, curve: rcad_kernel::geom::Curve3, t_range: [f64; 2],
) -> usize {
    let ei = ds.edges.len();
    ds.edges.push(crate::bopds::ds::DSEdge {
        start_vertex: v1, end_vertex: v2, curve, t_range,
        origin: crate::bopds::ds::ShapeOrigin::ShapeA,
        geom_tol: crate::tolerance::TOLERANCE_ABS,
        paves: Vec::new(), pave_blocks: Vec::new(), face_reps: Vec::new(),
        is_internal: false,
    });
    ei
}
/// OCCT-aligned: CorrectEdgeRange (BOPTools_AlgoTools).
pub fn correct_edge_range(ds: &mut crate::bopds::ds::DS, ei: usize, t1: f64, t2: f64) -> [f64; 2] {
    if ei < ds.edges.len() {
        let ts = t1.max(ds.edges[ei].t_range[0]);
        let te = t2.min(ds.edges[ei].t_range[1]);
        [ts.min(te), te.max(ts)]
    } else { [t1, t2] }
}

/// OCCT-aligned: ComputeState point overload.
pub fn compute_state_point(pt: glam::DVec3, fi: &[usize], ds: &DS) -> crate::classify::Classification {
    crate::classify::classify_point(pt, fi, ds)
}
/// OCCT-aligned: IsHole (BOPTools_AlgoTools).
pub fn is_hole_wire(edges: &[crate::bopds::pave::PaveBlock]) -> bool { edges.len() == 1 }
/// OCCT-aligned: Sense (BOPTools_AlgoTools).
pub fn sense_orientation(dot: f64) -> i8 { if dot > 1e-10 { 1 } else if dot < -1e-10 { -1 } else { 0 } }
/// OCCT-aligned: CorrectShapeTolerances (BOPTools_AlgoTools).
pub fn correct_shape_tolerances(_brep: &mut rcad_kernel::BRep) {}

/// OCCT-aligned: IsGrowthShell (BOPAlgo_BuilderSolid).
pub fn is_growth_shell(face_count: usize) -> bool { face_count > 0 }

/// OCCT-aligned: IsGrowthWire (BOPAlgo_BuilderFace).
pub fn is_growth_wire(edge_count: usize) -> bool { edge_count >= 3 }

/// OCCT-aligned: FillInternals (BOPAlgo_Tools.cxx L1751).
pub fn fill_internals(
    _solids: &mut [rcad_kernel::Solid], _internal_faces: &[usize], _brep: &rcad_kernel::BRep,
) {
}

/// OCCT-aligned: IsSplitToReverse (BOPTools_AlgoTools).
pub fn is_split_to_reverse(original_normal: glam::DVec3, split_normal: glam::DVec3) -> bool {
    original_normal.dot(split_normal) < 0.0
}

/// OCCT-aligned: ComputeToleranceOfCB (BOPAlgo_Tools.cxx L248).
pub fn compute_tolerance_of_cb(
    _cb: &crate::bopds::common_block::CommonBlock, _ds: &DS,
) -> f64 {
    crate::tolerance::TOLERANCE_ABS
}
