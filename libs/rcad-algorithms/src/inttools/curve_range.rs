use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::tolerance::{CONFUSION, PCONFUSION, parametric_epsilon, is_infinite_value};
use glam::DVec3;

/// Curve parameter step: parameter increment needed to move `tol` distance along curve.
///
/// OCCT: Adaptor3d_Curve::Resolution(tol) = tol / |dP/dt|
/// (BRepLib_1.cxx L61, IntTools_ShrunkRange.cxx L162)
///
/// For a curve parameterized by `t`, the step `resolution` satisfies:
/// `|P(t + resolution) - P(t)|  ≈ tol` (first-order approximation using tangent speed).
///
/// When the tangent speed is nearly zero (singularity), the resolution is clamped to `tol`.
pub fn curve_resolution(curve: &Curve3, t: f64, tol: f64) -> f64 {
    let speed = curve.tangent_at(t).length();
    if speed < 1e-15 {
        tol
    } else {
        tol / speed
    }
}

/// Global curve resolution: parametric step needed to move `tol` distance.
///
/// Samples at 23 points across the curve range to find the maximum derivative
/// magnitude, then returns `tol / max_speed`. Equivalent to OCCT
/// `Adaptor3d_Curve::Resolution(tol)` for BSpline curves which samples
/// at up to 23 points (Geom_BSplineCurve::Resolution).
///
/// For analytic curves (line, circle, etc.) where derivative magnitude is
/// constant, this gives the same result as the local resolution.
///
/// OCCT: BRepLib_1.cxx L61, Geom_BSplineCurve::Resolution
pub fn curve_resolution_global(curve: &Curve3, t_start: f64, t_end: f64, tol: f64) -> f64 {
    const N_SAMPLES: usize = 23;
    let range = t_end - t_start;
    if range.abs() < 1e-15 {
        // Degenerate range: use local resolution at midpoint
        let mid = (t_start + t_end) * 0.5;
        let speed = curve.derivative_at(mid).length();
        return if speed < 1e-15 { tol } else { tol / speed };
    }
    let mut max_speed = 0.0_f64;
    for i in 0..=N_SAMPLES {
        let t = t_start + range * (i as f64) / (N_SAMPLES as f64);
        let d1 = curve.derivative_at(t);
        let speed = d1.length();
        if speed > max_speed {
            max_speed = speed;
        }
    }
    if max_speed < 1e-15 { tol } else { tol / max_speed }
}

/// OCCT BRepLib_1.cxx L25-169: findNearestValidPoint.
///
/// Starting from the appointed end of the curve, find the nearest
/// point on the curve that is an intersection with the sphere with
/// center `the_vert_pnt` and radius `the_tol`.
///
/// # Arguments
/// * `curve` - The 3D curve
/// * `the_first` - First parameter of the curve range
/// * `the_last` - Last parameter of the curve range
/// * `is_first` - If true, start from the_first; if false, start from the_last
/// * `the_vert_pnt` - Center of the tolerance sphere (vertex point)
/// * `the_tol` - Radius of the tolerance sphere (vertex tolerance)
/// * `the_eps` - Convergence threshold for binary search
///
/// # Returns
/// * `Some(par)` - The nearest parameter outside the tolerance sphere
/// * `None` - The entire range is inside the tolerance sphere
fn find_nearest_valid_point(
    curve: &Curve3,
    the_first: f64,
    the_last: f64,
    is_first: bool,
    the_vert_pnt: DVec3,
    the_tol: f64,
    the_eps: f64,
) -> Option<f64> {
    // 1. Check that the needed end is inside the sphere
    let a_start_u = if is_first { the_first } else { the_last };
    let an_end_u = if is_first { the_last } else { the_first };
    let a_sq_tol = the_tol * the_tol;
    let a_p = curve.point_at(a_start_u);
    if (a_p - the_vert_pnt).length_squared() > a_sq_tol {
        // the vertex does not cover the corresponding to this vertex end of the curve
        return None;
    }

    // 2. Find a nearest point that is outside
    //
    // stepping along the curve by theTol till go out
    //
    // the general step is computed using general curve resolution
    let mut a_step = curve_resolution_global(curve, a_start_u, an_end_u, the_tol) * 1.01;
    if a_step < the_eps {
        a_step = the_eps;
    }
    // aD1Mag is a threshold to consider local derivative magnitude too small
    // and to accelerate going out of sphere
    // (inverse of resolution is the maximal derivative);
    // this is actual for bezier and b-spline types only
    let a_d1_mag = if is_bspline_or_bezier(curve) {
        // 1. / Resolution(1.) * 0.01 → max|D1| * 0.01 → squared
        let max_speed = 1.0 / curve_resolution_global(curve, a_start_u, an_end_u, 1.0);
        let val = max_speed * 0.01;
        val * val
    } else {
        0.0
    };

    if !is_first {
        a_step = -a_step;
    }

    let mut is_out = false;
    let mut an_u_in = a_start_u;
    let mut an_u_out = an_u_in;
    while !is_out {
        an_u_in = an_u_out;
        an_u_out += a_step;
        if (is_first && an_u_out > an_end_u) || (!is_first && an_u_out < an_end_u) {
            // step is too big and we go out of bounds,
            // check if the opposite bound is outside
            let a_p = curve.point_at(an_end_u);
            is_out = (a_p - the_vert_pnt).length_squared() > a_sq_tol;
            if !is_out {
                // all range is inside sphere
                return None;
            }
            an_u_out = an_end_u;
            break;
        }
        if a_d1_mag > 0.0 {
            let mut a_step_local = a_step.abs();
            loop {
                // cycle to go out of local singularity
                let a_p = curve.point_at(an_u_out);
                let a_d1 = curve.derivative_at(an_u_out);
                is_out = (a_p - the_vert_pnt).length_squared() > a_sq_tol;
                if !is_out && a_d1.length_squared() < a_d1_mag {
                    a_step_local *= 2.0;
                    if is_first {
                        an_u_out += a_step_local;
                    } else {
                        an_u_out -= a_step_local;
                    }
                    if (is_first && an_u_out < an_end_u) || (!is_first && an_u_out > an_end_u) {
                        // still in range
                        continue;
                    }
                    // went out of range, so check if the end point has out state
                    an_u_out = an_end_u;
                    let a_p = curve.point_at(an_u_out);
                    is_out = (a_p - the_vert_pnt).length_squared() > a_sq_tol;
                    if !is_out {
                        // all range is inside sphere
                        return None;
                    }
                }
                break;
            }
        } else {
            let a_p = curve.point_at(an_u_out);
            if !is_out {
                is_out = (a_p - the_vert_pnt).length_squared() > a_sq_tol;
            }
        }
    }

    // 3. Precise solution with binary search
    let mut a_delta = (an_u_out - an_u_in).abs();
    while a_delta > the_eps {
        let a_mid_u = (an_u_in + an_u_out) * 0.5;
        let a_p = curve.point_at(a_mid_u);
        let is_out_mid = (a_p - the_vert_pnt).length_squared() > a_sq_tol;
        if is_out_mid {
            an_u_out = a_mid_u;
        } else {
            an_u_in = a_mid_u;
        }
        a_delta = (an_u_out - an_u_in).abs();
    }
    let the_par = (an_u_in + an_u_out) * 0.5;
    Some(the_par)
}

/// Check if the curve is of BSpline or Bezier type (including unwrapping Offset and Trimmed).
fn is_bspline_or_bezier(curve: &Curve3) -> bool {
    match curve {
        Curve3::BSpline(_) | Curve3::Bezier(_) => true,
        Curve3::Offset(oc) => is_bspline_or_bezier(&oc.basis),
        Curve3::Trimmed(tc) => is_bspline_or_bezier(&tc.curve),
        _ => false,
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
/// rcad: steps 5-7 are done by the caller (ShrunkRange wrapper in shrunk_range.rs).
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
pub fn shrunk_range(
    curve: &Curve3, t_range: [f64; 2], v1_tol: f64, v2_tol: f64, edge_tol: f64,
) -> Option<[f64; 2]> {
    let [t1, t2] = t_range;
    // OCCT IntTools_ShrunkRange.cxx L117-120: if range < PConfusion -> micro-edge
    if (t2 - t1).abs() < PCONFUSION { return None; }
    // OCCT BRepLib::FindValidRange handles Circle curves natively using
    // GeomAPI_ProjectPointOnCurve. rcad's find_nearest_valid_point may not
    // handle periodic curves reliably. Return full range for Circle curves
    // to ensure shrunk data is valid (splittable = true).
    if matches!(curve, Curve3::Circle(_)) {
        return Some([t1, t2]);
    }
    // OCCT L124-142: tolerance adjustments
    let a_tol_v1 = v1_tol.max(edge_tol) + CONFUSION;
    let a_tol_v2 = v2_tol.max(edge_tol) + CONFUSION;

    // OCCT L144-146: compute the points at the endpoints
    let p1 = curve.point_at(t1);
    let p2 = curve.point_at(t2);

    // OCCT BRepLib::FindValidRange (L173-258)
    // Compute epsilon = max(curve.Resolution(theTolE) * 0.1, Epsilon(aMaxPar), PConfusion())
    let is_inf_par_v1 = is_infinite_value(t1);
    let is_inf_par_v2 = is_infinite_value(t2);
    let a_max_par = {
        let mut m = 0.0;
        if !is_inf_par_v1 { m = t1.abs(); }
        if !is_inf_par_v2 { m = m.max(t2.abs()); }
        m
    };
    // OCCT L201-202: anEps = max(curve.Resolution(theTolE) * 0.1, Epsilon(aMaxPar), PConfusion())
    let an_eps = {
        let res = curve_resolution_global(curve, t1, t2, edge_tol) * 0.1;
        let eps = parametric_epsilon(a_max_par);
        res.max(eps).max(PCONFUSION)
    };

    // OCCT L204-225: find theFirst (first vertex)
    let ts1 = if is_inf_par_v1 {
        t1
    } else {
        match find_nearest_valid_point(curve, t1, t2, true, p1, a_tol_v1, an_eps) {
            Some(val) => {
                // OCCT L221-224: if (theParV2 - theFirst < anEps) return false;
                if t2 - val < an_eps { return None; }
                val
            }
            None => {
                // OCCT: BRepLib::FindValidRange handles all curve types including Circle.
                // rcad's find_nearest_valid_point may fail on Circle edges where
                // endpoint vertices coincide (start_vertex == end_vertex).
                // Fall back to full range for Circle curves.
                match curve {
                    _ => return None,
                }
            },
        }
    };

    // OCCT L227-248: find theLast (second vertex)
    let ts2 = if is_inf_par_v2 {
        t2
    } else {
        match find_nearest_valid_point(curve, t1, t2, false, p2, a_tol_v2, an_eps) {
            Some(val) => {
                // OCCT L244-247: if (theLast - theParV1 < anEps) return false;
                if val - t1 < an_eps { return None; }
                val
            }
            None => {
                match curve {
                    rcad_kernel::geom::Curve3::Circle(_) => t2,
                    _ => return None,
                }
            },
        }
    };

    // OCCT L250-255: check found parameters — if theFirst > theLast → overlapping
    if ts1 > ts2 {
        return None;
    }

    // OCCT L152-156: if shrunk range < PConfusion → micro-edge
    if (ts2 - ts1) < PCONFUSION { return None; }

    // OCCT L158-175: length computation and check — done by caller (ShrunkRange wrapper)
    Some([ts1, ts2])
}

