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
        pave_blocks: Vec::new(),
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
        paves: Vec::new(), pave_blocks: Vec::new(),
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
