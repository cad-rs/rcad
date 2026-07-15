/// Greedy nearest-neighbor ordering of a point cloud into one or more chains.
///
/// Returns a `Vec` of chains (each chain is `Vec<DVec3>`).
/// Points that can't be extended within `gap_tol` start a new chain.
/// After chain formation, chains are stitched together when their endpoints
/// are close enough, producing fewer, longer chains (typically one closed loop
/// for a single intersection curve).
fn greedy_order_points(pts: Vec<DVec3>, gap_floor: f64) -> Vec<Vec<DVec3>> {
 if pts.is_empty() {
 return vec![];
 }

 let gap_floor = if gap_floor.is_finite() && gap_floor > 0.0 {
 gap_floor
 } else {
 0.0
 };

 // Estimate gap tolerance from average nearest-neighbor distance
 // (rough: use 3x the median distance between sorted x-coordinates).
 // Also compute the bounding-box diagonal to prevent gap_tol from
 // shrinking excessively when dense analytic-distance crossing points
 // reduce the median nn-distance (e.g. cone-cylinder pairs).
 let gap_tol = {
 let mut dists: Vec<f64> = Vec::with_capacity(pts.len());
 let mut bbox_min = DVec3::splat(f64::INFINITY);
 let mut bbox_max = DVec3::splat(f64::NEG_INFINITY);
 for i in 0..pts.len() {
 let pi = pts[i];
 bbox_min = bbox_min.min(pi);
 bbox_max = bbox_max.max(pi);
 let mut best = f64::INFINITY;
 for j in 0..pts.len() {
 if i != j {
 let d = (pi - pts[j]).length();
 if d < best {
 best = d;
 }
 }
 }
 dists.push(best);
 }
 dists.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
 let median = dists[dists.len() / 2];
 let bbox_diag = (bbox_max - bbox_min).length();
 // Floor by 2% of the bounding-box diagonal so gap_tol doesn't collapse
 // below  ?% of the curve's spatial extent (handles analytic-distance
 // oversampling while keeping Steinmetz and other small intersections intact).
 let bbox_floor = bbox_diag * 0.02;
 (median * 5.0).max(bbox_floor).max(TOLERANCE_COORD_SUB).max(gap_floor)
 };

 // Stitch gap tolerance: more generous than within-chain growth.
 // Allow up to 3  the within-chain gap to merge separate chains.
 let stitch_tol = gap_tol * 3.0;

 let mut used = vec![false; pts.len()];
 let mut chains: Vec<Vec<DVec3>> = Vec::new();

 loop {
 // Find first unused point
 let start = match used.iter().position(|&u| !u) {
 Some(i) => i,
 None => break,
 };
 used[start] = true;
 let mut chain = vec![pts[start]];

 loop {
 let last = *chain.last().expect("chain is non-empty (starts with 1 element)");
 // Find nearest unused point within gap_tol
 let mut best_dist = gap_tol;
 let mut best_idx = None;
 for (i, &used_i) in used.iter().enumerate() {
 if !used_i {
 let d = (pts[i] - last).length();
 if d < best_dist {
 best_dist = d;
 best_idx = Some(i);
 }
 }
 }
 match best_idx {
 Some(idx) => {
 used[idx] = true;
 chain.push(pts[idx]);
 }
 None => break,
 }
 }

 chains.push(chain);
 }

 //  € € Chain stitching  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
 // Repeatedly merge pairs of chains whose endpoints are within stitch_tol.
 // This turns fragmented arc segments into a single closed loop.
 let mut changed = true;
 while changed && chains.len() > 1 {
 changed = false;
 'outer: for i in 0..chains.len() {
 for j in (i + 1)..chains.len() {
 let end_i = *chains[i].last().expect("chains[i] is non-empty");
 let start_j = chains[j][0];
 let end_j = *chains[j].last().expect("chains[j] is non-empty");
 let start_i = chains[i][0];

 // Determine merge direction
 let (merge_rev_j, close_enough) =
 if (end_i - start_j).length() <= stitch_tol {
 (false, true) // i + j
 } else if (end_i - end_j).length() <= stitch_tol {
 (true, true)  // i + reversed j
 } else if (end_j - start_i).length() <= stitch_tol {
 // j + i: handled next iteration via swapped roles
 (false, false)
 } else if (start_j - start_i).length() <= stitch_tol {
 (false, false)
 } else {
 (false, false)
 };

 if close_enough {
 let chain_j = chains.remove(j);
 let appended: Vec<DVec3> = if merge_rev_j {
 chain_j.into_iter().rev().collect()
 } else {
 chain_j
 };
 chains[i].extend(appended);
 changed = true;
 break 'outer;
 }

 // Also handle: j ends near i start  ?prepend j to i
 if (end_j - start_i).length() <= stitch_tol {
 let chain_j = chains.remove(j);
 let mut merged = chain_j;
 merged.append(&mut chains[i]);
 chains[i] = merged;
 changed = true;
 break 'outer;
 }
 // j start near i start  ?prepend reversed j
 if (start_j - start_i).length() <= stitch_tol {
 let chain_j = chains.remove(j);
 let mut merged: Vec<DVec3> = chain_j.into_iter().rev().collect();
 merged.append(&mut chains[i]);
 chains[i] = merged;
 changed = true;
 break 'outer;
 }
 }
 }
 }

 chains
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// UV-space Newton-Raphson refinement (IntPatch_TheSearchInside)
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// UV-space Newton-Raphson refinement of an intersection point.
///
/// Given an initial guess `(u, v)` on `s1`, find `(u', v')` where
/// F(u',v')  ?0, where F(u,v) = signed_distance(s2, s1.point_at(u,v)).
///
/// OCCT equivalent: `IntPatch_TheSurfFunction::Value` + `math_FunctionSetRoot`
/// in `IntStart_SearchInside::Perform` (IntStart_SearchInside.gxx L211).
/// OCCT uses analytic gradients; rcad uses finite differences for generality
/// across all surface types (Surface3 does not expose partial derivatives).
fn refine_uv_intersection(
 s1: &Surface3,
 s2: &Surface3,
 u: f64,
 v: f64,
 tol: f64,
) -> Option<(DVec3, f64, f64)> {
 let tol = tol.abs().max(TOLERANCE_ABS);
 let eps = 1e-6; // finite-difference step in UV
 let max_iter = 20;
 let mut u = u;
 let mut v = v;

 for _ in 0..max_iter {
 let p = s1.point_at(u, v);
 if !p.is_finite() {
 break;
 }
 let f = surface_implicit(s2, p);

 if f.abs() < tol {
 return Some((p, u, v));
 }

 // Finite-difference Jacobian: dF/du, dF/dv
 let pu = s1.point_at(u + eps, v);
 let pv = s1.point_at(u, v + eps);
 let fu = if pu.is_finite() { surface_implicit(s2, pu) } else { f };
 let fv = if pv.is_finite() { surface_implicit(s2, pv) } else { f };

 let df_du = (fu - f) / eps;
 let df_dv = (fv - f) / eps;
 let grad2 = df_du * df_du + df_dv * df_dv;

 if grad2 < TOLERANCE_LEN_SQ_DIV_SAFE {
 break;
 }

 // Gauss-Newton step: (u,v) += -F * / || || 
 let du = -f * df_du / grad2;
 let dv = -f * df_dv / grad2;

 if du.abs() < 1e-10 && dv.abs() < 1e-10 {
 // Converged to machine precision
 let pf = s1.point_at(u, v);
 let ff = if pf.is_finite() { surface_implicit(s2, pf) } else { f };
 if ff.abs() < tol {
 return Some((pf, u, v));
 }
 break;
 }

 u += du;
 v += dv;
 }

 None
}

/// IntStart_SearchInside  ?constrained Newton-Raphson on F(u,v)=0
/// within UV bounds [u0, u1] [v0, v1] (matching OCCT's Binf/Bsup box).
/// After convergence, checks |F(u,v)| <= func_tol.
fn refine_uv_intersection_bounded(
 s1: &Surface3,
 s2: &Surface3,
 u: f64,
 v: f64,
 u0: f64,
 u1: f64,
 v0: f64,
 v1: f64,
 func_tol: f64,
) -> Option<(DVec3, f64, f64)> {
 let func_tol = func_tol.abs().max(TOLERANCE_ABS);
 let eps = 1e-6;
 let max_iter = 20;
 let mut u = u.clamp(u0, u1);
 let mut v = v.clamp(v0, v1);

 for _ in 0..max_iter {
 let p = s1.point_at(u, v);
 if !p.is_finite() {
 break;
 }
 let f = surface_implicit(s2, p);

 if f.abs() < func_tol {
 return Some((p, u, v));
 }

 let pu = s1.point_at((u + eps).min(u1), v);
 let pv = s1.point_at(u, (v + eps).min(v1));
 let fu = if pu.is_finite() { surface_implicit(s2, pu) } else { f };
 let fv = if pv.is_finite() { surface_implicit(s2, pv) } else { f };

 let df_du = (fu - f) / eps;
 let df_dv = (fv - f) / eps;
 let grad2 = df_du * df_du + df_dv * df_dv;

 if grad2 < TOLERANCE_LEN_SQ_DIV_SAFE {
 break;
 }

 let du = -f * df_du / grad2;
 let dv = -f * df_dv / grad2;

 if du.abs() < 1e-10 && dv.abs() < 1e-10 {
 let pf = s1.point_at(u, v);
 let ff = if pf.is_finite() { surface_implicit(s2, pf) } else { f };
 if ff.abs() < func_tol {
 return Some((pf, u, v));
 }
 break;
 }

 u += du;
 v += dv;
 // clamp to search box (Binf/Bsup)
 u = u.clamp(u0, u1);
 v = v.clamp(v0, v1);
 }

 None
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// polyline  ?BSpline curve fitting
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// polyline-to-BSpline curve fitting for surface-surface intersection
/// results.
///
/// OCCT reference: `GeomInt_IntSS::MakeBSpline` (GeomInt_IntSS.cxx, lines ~180-280).
///
/// OCCT algorithm (simplified):
/// 1. Approximate the polyline with a BSpline using least-squares fitting
/// 2. Check the approximation error at each polyline point
/// 3. If max error > tolerance, increase control points and re-fit
/// 4. The result is a C2-continuous (degree  ?3) BSpline curve
///
/// This implementation:
/// - Uses chord-length parameterization (matching OCCT's approach)
/// - Fits a cubic (degree-3, or  ?2 for small inputs) BSpline via least-squares
/// - Progressively increases the number of control points when deviation
/// exceeds `max_tol`, up to exact interpolation (zero deviation)
/// - Returns `Some(Curve3::BSpline(...))` on success, or `None` when there
/// are too few points or fitting fails
///
///  ?Same least-squares + adaptive refinement strategy as
/// `GeomInt_IntSS::MakeBSpline`. The chord-length parameterization and
/// clamped cubic knot vector match OCCT's internal approach.
pub fn polyline_to_bspline(points: &[DVec3], max_tol: f64) -> Option<Curve3> {
 let n = points.len();
 if n < 4 {
 // Too few points for a meaningful BSpline fit.
 return None;
 }

 // Compute chord-length parameters for deviation checking.
 // These map each input point to its parameter value in [0, 1].
 let params = chord_length_params_3d(points)?;

 // Strategy: progressively refine by increasing the number of control points.
 // Start with ~half the points as control points, up to n.
 let mut n_ctrl = (n / 2).clamp(4, n);

 for _ in 0..8 {
 // When n_ctrl >= n, fall back to exact interpolation (zero deviation).
 let bspline = if n_ctrl >= n {
 rcad_kernel::fit::interpolate_points(points).ok()?
 } else {
 rcad_kernel::fit::approximate_points(points, n_ctrl).ok()?
 };

 // Compute max deviation at polyline points.
 let max_dev = max_bspline_deviation(&bspline, points, &params);

 if max_dev <= max_tol {
 //  ?Deviation within tolerance  ?accept the BSpline.
 return Some(Curve3::BSpline(bspline));
 }

 // Exact interpolation with n_ctrl == n should achieve near-zero
 // deviation (limited only by floating-point precision). If we are
 // already at n_ctrl >= n, return the exact fit anyway.
 if n_ctrl >= n {
 return Some(Curve3::BSpline(bspline));
 }

 // Increase control points: move halfway from current toward n.
 let remaining = n - n_ctrl;
 let increment = remaining / 2;
 n_ctrl = n.min(n_ctrl + increment.max(1));
 }

 None
}

/// Compute chord-length parameterization for 3D points, normalized to [0, 1].
///
/// Returns `None` when the points are degenerate (all coincident).
fn chord_length_params_3d(pts: &[DVec3]) -> Option<Vec<f64>> {
 let n = pts.len();
 let mut params = Vec::with_capacity(n);
 params.push(0.0);
 let mut total = 0.0;
 for i in 1..n {
 total += (pts[i] - pts[i - 1]).length();
 params.push(total);
 }
 if total < 1e-14 {
 return None;
 }
 for p in &mut params {
 *p /= total;
 }
 Some(params)
}

/// Compute the maximum deviation between a BSpline curve and the input data
/// points at their chord-length parameter values.
fn max_bspline_deviation(bspline: &BSplineCurve3, data_pts: &[DVec3], params: &[f64]) -> f64 {
 let mut max_dev = 0.0;
 for (i, pt) in data_pts.iter().enumerate() {
 let eval_pt = bspline.point_at(params[i]);
 let dev = (*pt - eval_pt).length();
 if dev > max_dev {
 max_dev = dev;
 }
 }
 max_dev
}

/// Convert all `SurfaceCurve::Polyline` entries in a
/// `SurfaceSurfaceIntersection` to BSpline approximations when beneficial.
///
///  ?corresponds to the post-processing step in
/// `GeomInt_IntSS::Perform` that replaces raw polylines with BSpline curves.
///
/// Call this after [intersect_surfaces_with_density] or any surface-surface
/// intersection that may produce polyline segments.
pub fn convert_polylines_to_bsplines(
 result: &mut SurfaceSurfaceIntersection,
 max_tol: f64,
) {
 for entry in &mut result.curves {
 if let SurfaceCurve::Polyline(pts) = &entry.curve_3d {
 if pts.len() >= 4 {
 if let Some(bspline) = polyline_to_bspline(pts, max_tol) {
 if let Curve3::BSpline(b) = bspline {
 entry.curve_3d = SurfaceCurve::BSplineCurve(Box::new(b));
 }
 }
 }
 }
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Tests
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

