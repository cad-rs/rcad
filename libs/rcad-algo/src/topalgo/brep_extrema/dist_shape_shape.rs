// OCCT BRepExtrema_DistShapeShape (BRepExtrema_DistShapeShape.cxx) — 1:1
// translation of the edge-edge MIN case, as used by IntTools_EdgeEdge::Perform()
// (BOPAlgo_PaveFiller fast rejection: "if (d > 1.1 * myTol) return;").
//
// OCCT computes the minimum distance between two edges as the minimum over:
//   - the curve-curve interior extrema clipped to the ranges
//     (Extrema_ExtCC -> Extrema_ExtElC for line-conic; OCCT L177-317,
//      L832-901 PrepareResults),
//   - the edge vertex sub-shapes (vertex-edge and vertex-vertex distances,
//     OCCT DistanceSS / DistShapeShape Perform).
//
// rcad: the boolean DS represents an edge as Curve3 + parameter range, so the
// edge-edge distance is a function of two curve segments.

use glam::DVec3;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_cc::ExtremaExtCC;
use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{Curve3, CurveEval};

/// OCCT BRepExtrema_DistShapeShape(edge1, edge2, Extrema_ExtFlag_MIN).Value()
/// for the edge-edge case: the exact minimum 3D distance between the curve
/// segment [t1, t2] of `c1` and [u1, u2] of `c2`.
pub fn min_distance_edge_segments(
    c1: &Curve3,
    t1: f64,
    t2: f64,
    c2: &Curve3,
    u1: f64,
    u2: f64,
) -> f64 {
    // OCCT: BRepAdaptor_Curve C1(edge1)/C2(edge2) with the edge ranges.
    let a_c1 = GeomCurveAdaptor::with_range(c1.clone(), t1, t2);
    let a_c2 = GeomCurveAdaptor::with_range(c2.clone(), u1, u2);
    let a_tool1 = CurveToolHandle::for_curve3(c1, &a_c1, &a_c1);
    let a_tool2 = CurveToolHandle::for_curve3(c2, &a_c2, &a_c2);

    // OCCT: Extrema_ExtCC ExtCC(C1, C2, U1, U2, V1, V2, TolC1, TolC2)
    // (ExtCC.cxx L177-317). The analytical ExtElC branches run for a line
    // against an elementary curve (cxx L247-294), the general branch
    // otherwise; PrepareResults clips to the ranges (cxx L832-901).
    let an_ext_cc = ExtremaExtCC::new_curves_ranged(
        &a_tool1, &a_tool2, t1, t2, u1, u2, 1.0e-10, 1.0e-10,
    );

    let mut best = f64::INFINITY;

    // OCCT Extrema_ExtCC::NbExt() / SquareDistance(i) (cxx L351-358 / L340-347).
    if an_ext_cc.is_done() {
        let a_nb_ext = an_ext_cc.nb_ext();
        for an_idx in 1..=a_nb_ext {
            best = best.min(an_ext_cc.square_distance(an_idx));
        }
    }
    // OCCT Extrema_ExtCC::TrimmedSquareDistances (mydist11/12/21/22,
    // cxx L375-393).
    let (a_d11, a_d12, a_d21, a_d22, _, _, _, _) = an_ext_cc.trimmed_square_distances();
    best = best.min(a_d11);
    best = best.min(a_d12);
    best = best.min(a_d21);
    best = best.min(a_d22);

    // OCCT BRepExtrema_DistShapeShape: the edge vertex sub-shapes — endpoint of
    // one edge to the other edge's curve (ExtPC / ExtPElC).
    for &te in &[t1, t2] {
        let p = c1.point_at(te);
        let a_d = closest_point_on_curve_range(c2, p, u1, u2, 64).distance;
        best = best.min(a_d * a_d);
    }
    for &ue in &[u1, u2] {
        let p = c2.point_at(ue);
        let a_d = closest_point_on_curve_range(c1, p, t1, t2, 64).distance;
        best = best.min(a_d * a_d);
    }

    best.sqrt()
}

/// OCCT BRepExtrema_DistShapeShape(edge, vertex).Value() — the minimum 3D
/// distance from the curve segment [t1, t2] to the point (ExtP / ExtPC
/// clipped to the range, BRepExtrema_DistShapeShape.cxx Perform(Edge, Vertex)).
pub fn min_distance_edge_vertex(c: &Curve3, t1: f64, t2: f64, p: DVec3) -> f64 {
    closest_point_on_curve_range(c, p, t1, t2, 64).distance
}
