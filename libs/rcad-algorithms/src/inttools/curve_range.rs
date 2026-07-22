use rcad_kernel::geom::{Curve3, CurveEval};
use crate::tolerance::TOLERANCE_CLAMP_MIN;

/// Curve parameter step: parameter increment needed to move `tol` distance along curve.
///
/// OCCT: Adaptor3d_Curve::Resolution(tol) = tol / |dP/dt|
/// (BRepLib_1.cxx L61, IntTools_ShrunkRange.cxx L162)
///
/// For a curve parameterized by `t`, the step `resolution` satisfies:
/// `|P(t + resolution) - P(t)|  ?tol` (first-order approximation using tangent speed).
///
/// When the tangent speed is nearly zero (singularity), the resolution is clamped to `tol`.
pub fn curve_resolution(curve: &Curve3, t: f64, tol: f64) -> f64 {
    let speed = curve.tangent_at(t).length();
    if speed < TOLERANCE_CLAMP_MIN {
        tol
    } else {
        tol / speed
    }
}

/// Compute the shrunk (valid) range for a curve segment, excluding the tolerance
/// spheres around the endpoint vertices.
///
/// IntTools_ShrunkRange::Perform() (IntTools_ShrunkRange.cxx L107-191)
///              + BRepLib::FindValidRange (BRepLib_1.cxx L173-258)
///
/// The shrunk range [t_start, t_end] is the portion of [t1, t2] where the curve point
/// is outside both vertex tolerance spheres. This is critical for:
/// - Edge-Face intersection: the shrunk range tells us where the edge is "truly away"
///   from its endpoint vertices, avoiding false intersections near vertices
///   (OCCT IntTools_ShrunkRange purpose).
/// - Micro-edge detection: if the shrunk range is empty, the PaveBlock is too short
///   and should be removed (OCCT `!IsSplittable()` check).
///
/// OCCT algorithm (IntTools_ShrunkRange::Perform):
/// 1. Guard: return None if (t2 - t1) < Precision::PConfusion() (L117-120)
/// 2. Get vertex tolerances; for each vertex:
///    aTolV = max(aTolV, aTolE) + Precision::Confusion() (L129-142)
/// 3. Call BRepLib::FindValidRange which:
///    a. Computes eps = max(curve.Resolution(aTolE) * 0.1, Epsilon(aMaxPar),
///                          Precision::PConfusion())
///    b. For each endpoint, calls findNearestValidPoint that steps along the curve
///       until outside the tolerance sphere, then binary-search refines.
///    c. Returns (theFirst, theLast) with theFirst < theLast.
/// 4. Guard: return None if (myTS2 - myTS1) < Precision::PConfusion() (L152-155)
/// 5. Compute edge length on shrunk range (L159-170)
/// 6. Guard: return None if length < Precision::Confusion() (L171-174)
/// 7. Set is_splittable if length > 2*aTolE + 2*Precision::Confusion() (L184-187)
///
/// This implementation provides the core shrunk range computation (steps 1-4).
/// It uses curve_resolution at the endpoint parameters as an O(1) approximation
/// of the full stepping + binary search approach, which is accurate when the curve
/// speed is nearly constant near the endpoints.
///
/// # Arguments
/// * `curve` - The 3D curve
/// * `t_range` - The full parameter range [t1, t2] (t1 < t2)
/// * `v1_tol` - Geometric tolerance at the start vertex
/// * `v2_tol` - Geometric tolerance at the end vertex
/// * `edge_tol` - Geometric tolerance of the edge
///
/// # Returns
/// * `Some([t_start, t_end])` - Valid shrunk range where t_start < t_end
/// * `None` - Micro-edge; the entire range is covered by tolerance spheres
use glam::DVec3;

/// Step along curve, return first parameter outside tolerance sphere.
fn find_nearest_valid_point(curve: &Curve3, t_start: f64, t_end: f64, center: DVec3, tol: f64, step: f64) -> Option<f64> {
    let tol_sq = tol * tol;
    let mut t = t_start;
    let dir = if t_end > t_start { 1.0 } else { -1.0 };
    // OCCT BRepLib_1.cxx L70: bounds iterations to ~100.
    let max_iter = 4096;
    let mut iter = 0usize;
    while iter < max_iter && (t - t_start) * (t_end - t_start) >= 0.0 {
        iter += 1;
        if (curve.point_at(t) - center).length_squared() >= tol_sq {
            if (t - t_start).abs() < step * 0.5 { return None; }
            let mut lo = if dir > 0.0 { t - step } else { t + step };
            let mut hi = t;
            for _ in 0..12 { let mid = (lo+hi)*0.5;
                if (curve.point_at(mid)-center).length_squared() >= tol_sq { hi = mid; } else { lo = mid; } }
            return Some(hi);
        }
        let next_t = t + step * dir;
        if next_t == t { break; } // step too small for fp resolution
        t = next_t;
    }
    None
}

pub fn shrunk_range(
    curve: &Curve3, t_range: [f64; 2], v1_tol: f64, v2_tol: f64, edge_tol: f64,
) -> Option<[f64; 2]> {
    let [t1, t2] = t_range;
    // OCCT IntTools_ShrunkRange.cxx L117-120: if range < PConfusion → micro-edge
    if (t2 - t1).abs() < rcad_kernel::tolerance::CONFUSION { return None; }
    // OCCT L124-142: tolerance adjustments
    let confusion = rcad_kernel::tolerance::CONFUSION;
    let a_tol_v1 = v1_tol.max(edge_tol) + confusion;
    let a_tol_v2 = v2_tol.max(edge_tol) + confusion;
    // Compute parametric steps for vertex spheres (OCCT: BRepLib_1.cxx L61: Resolution)
    // OCCT L162-169: parametric tolerance for edge
    let step1 = curve_resolution(curve, t1, a_tol_v1);
    let step2 = curve_resolution(curve, t2, a_tol_v2);
    let step1 = step1 * 5.0;
    let step2 = step2 * 5.0;
    // OCCT L146-151: FindValidRange — find parametric bounds outside vertex spheres
    let p1 = curve.point_at(t1);
    let p2 = curve.point_at(t2);
    let ts1 = find_nearest_valid_point(curve, t1, t2, p1, a_tol_v1, step1);
    let ts2 = find_nearest_valid_point(curve, t2, t1, p2, a_tol_v2, step2);
    let (the_first, the_last) = match (ts1, ts2) {
        (Some(f), Some(l)) => (f, l), (Some(f), None) => (f, t2),
        (None, Some(l)) => (t1, l), (None, None) => return None,
    };
    // OCCT L152-156: if shrunk range < PConfusion → micro-edge
    if the_first >= the_last || (the_last - the_first) < rcad_kernel::tolerance::CONFUSION { return None; }
    // OCCT L162-169: parametric tolerance for length computation
    // OCCT L170-175: GCPnts_AbscissaPoint::Length — compute arc length of shrunk range
    // rcad: skip length check — edges with valid FindValidRange are accepted.
    // OCCT L184-187: check splittable (length > 2*TolE + 2*Confusion)
    // rcad: splittable check done by caller (analyze_shrunk_data).
    Some([the_first, the_last])
}

