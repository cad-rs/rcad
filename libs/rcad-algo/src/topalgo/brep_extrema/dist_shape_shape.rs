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

use rcad_kernel::base::geom_api::project::closest_point_on_curve_range;
use rcad_kernel::base::extrema::{ext_cc_line_conic, line_line_extrema};
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
    let mut best = f64::INFINITY;

    // OCCT Extrema_ExtCC dispatch (L247-294): line-conic handled analytically
    // by Extrema_ExtElC. Normalize so a Line comes first.
    let (line, conic, lt1, lt2, cu1, cu2) = match (c1, c2) {
        (Curve3::Line(l), c2) => (l, c2, t1, t2, u1, u2),
        (c1, Curve3::Line(l)) => (l, c1, u1, u2, t1, t2),
        (Curve3::Line(a), Curve3::Line(b)) => {
            // Line-Line (OCCT ExtElC L268-357): interior closest pair + vertices.
            for (d, tp, up) in line_line_extrema(a, b) {
                if tp >= t1 - f64::EPSILON && tp <= t2 + f64::EPSILON
                    && up >= u1 - f64::EPSILON && up <= u2 + f64::EPSILON
                {
                    best = best.min(d);
                }
            }
            for &te in &[t1, t2] {
                best = best.min(closest_point_on_curve_range(c2, a.point_at(te), u1, u2, 64).distance);
            }
            for &ue in &[u1, u2] {
                best = best.min(closest_point_on_curve_range(c1, b.point_at(ue), t1, t2, 64).distance);
            }
            return best;
        }
        // Both non-line: not reachable from the IntTools_EdgeEdge fast-reject
        // (its condition guarantees one curve is a Line). Fall back to the
        // sampling curve-curve extrema.
        _ => return rcad_kernel::base::extrema::extrema_curve_curve(c1, c2, 64).min_distance(),
    };

    // OCCT Extrema_ExtCC for line-conic: interior extrema in range + corners.
    let cc = ext_cc_line_conic(line, lt1, lt2, conic, cu1, cu2);
    for (d, _, _) in &cc.interior {
        best = best.min(*d);
    }
    // OCCT TrimmedSquareDistances (mydist11/12/21/22, L375-393).
    best = best.min(cc.corners.dist11.sqrt());
    best = best.min(cc.corners.dist12.sqrt());
    best = best.min(cc.corners.dist21.sqrt());
    best = best.min(cc.corners.dist22.sqrt());

    // OCCT vertex sub-shapes of the edges: endpoint of one edge to the other
    // edge's curve (ExtPC / ExtPElC).
    for &te in &[lt1, lt2] {
        let p = line.point_at(te);
        best = best.min(closest_point_on_curve_range(conic, p, cu1, cu2, 64).distance);
    }
    for &ue in &[cu1, cu2] {
        let p = conic.point_at(ue);
        best = best.min(closest_point_on_curve_range(&Curve3::Line(*line), p, lt1, lt2, 64).distance);
    }

    best
}
