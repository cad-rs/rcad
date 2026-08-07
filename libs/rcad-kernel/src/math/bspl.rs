//! OCCT BSplCLib + BSplSLib: BSpline curve and surface evaluation.
//!
//! Core algorithms:
//! - Cox-De Boor evaluation (rational NURBS, 3D and 2D)
//! - BSpline derivative via homogeneous quotient rule
//! - Knot span binary search (same knot vector convention as OCCT)
//!
//! OCCT source: src/FoundationClasses/TKMath/BSplCLib/BSplCLib.cxx
//!             src/FoundationClasses/TKMath/BSplSLib/BSplSLib.cxx

use crate::geom::{BezierCurve3, BezierSurface, BSplineCurve3, BSplineSurface};
use glam::{DVec2, DVec3};

/// Find the knot span index `k` such that `knots[k] <= t < knots[k+1]`.
/// OCCT: BSplCLib::BinSearch — used by all BSpline evaluators.
pub fn find_knot_span(degree: usize, knots: &[f64], t: f64) -> usize {
    let t_min = knots[degree];
    let t_max = knots[knots.len() - degree - 1];
    let t_clamped = t.clamp(t_min, t_max);
    let mut span = degree;
    for (i, &knot) in knots.iter().enumerate().take(knots.len() - degree - 1).skip(degree) {
        if knot <= t_clamped { span = i; } else { break; }
    }
    span
}

// ══════════════════════════════════════════════════════════════════════════
// BSplCLib — curve evaluation
// ══════════════════════════════════════════════════════════════════════════

/// OCCT BSplCLib: weight accessor — a NULL weight array (non-rational curve)
/// is treated as all weights equal to 1.0.
#[inline]
fn wgt(weights: &[f64], i: usize) -> f64 {
    if weights.is_empty() {
        1.0
    } else {
        weights[i]
    }
}

/// Cox-De Boor evaluation of a rational BSpline curve (3D).
/// OCCT: BSplCLib::Eval.
pub fn de_boor(degree: usize, knots: &[f64], points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 { return DVec3::ZERO; }
    let k = find_knot_span(degree, knots, t);
    let mut r = vec![DVec3::ZERO; degree + 1];
    let mut w = vec![0.0f64; degree + 1];
    for j in 0..=degree {
        let idx = k - degree + j;
        let wgt_j = wgt(weights, idx);
        r[j] = points[idx] * wgt_j;
        w[j] = wgt_j;
    }
    for level in 1..=degree {
        for j in 0..=(degree - level) {
            let idx_a = k - degree + j + level;
            let a = (t - knots[idx_a]) / (knots[idx_a + degree - level + 1] - knots[idx_a]);
            let a = a.clamp(0.0, 1.0);
            r[j] = r[j] * (1.0 - a) + r[j + 1] * a;
            w[j] = w[j] * (1.0 - a) + w[j + 1] * a;
        }
    }
    if w[0].abs() > 1e-15 { r[0] / w[0] } else { r[0] }
}

/// Cox-De Boor evaluation of a rational BSpline curve (2D).
/// OCCT: BSplCLib::Eval (2D overload).
pub fn de_boor_2d(degree: usize, knots: &[f64], points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 { return DVec2::ZERO; }
    let k = find_knot_span(degree, knots, t);
    let mut r = vec![DVec2::ZERO; degree + 1];
    let mut w = vec![0.0f64; degree + 1];
    for j in 0..=degree {
        let idx = k - degree + j;
        let wgt_j = wgt(weights, idx);
        r[j] = points[idx] * wgt_j;
        w[j] = wgt_j;
    }
    for level in 1..=degree {
        for j in 0..=(degree - level) {
            let idx_a = k - degree + j + level;
            let a = (t - knots[idx_a]) / (knots[idx_a + degree - level + 1] - knots[idx_a]);
            let a = a.clamp(0.0, 1.0);
            r[j] = r[j] * (1.0 - a) + r[j + 1] * a;
            w[j] = w[j] * (1.0 - a) + w[j + 1] * a;
        }
    }
    if w[0].abs() > 1e-15 { r[0] / w[0] } else { r[0] }
}

/// Cox-De Boor evaluation returning a homogeneous 4-vector `(wx, wy, wz, w)`.
/// Used by the surface evaluator (tensor product) to postpone division.
/// OCCT: BSplCLib::Eval with homogeneous output.
pub fn de_boor_homo(degree: usize, knots: &[f64], points: &[DVec3], weights: &[f64], t: f64) -> [f64; 4] {
    let n = points.len();
    if n == 0 { return [0.0; 4]; }
    let k = find_knot_span(degree, knots, t);
    let mut r = vec![[0.0f64; 4]; degree + 1];
    for j in 0..=degree {
        let idx = k - degree + j;
        let p = points[idx];
        let w = wgt(weights, idx);
        r[j] = [p.x * w, p.y * w, p.z * w, w];
    }
    for level in 1..=degree {
        for j in 0..=(degree - level) {
            let idx_a = k - degree + j + level;
            let denom = knots[idx_a + degree - level + 1] - knots[idx_a];
            let a = if denom.abs() > 1e-15 { ((t - knots[idx_a]) / denom).clamp(0.0, 1.0) } else { 0.0 };
            for c in 0..4 {
                r[j][c] = r[j][c] * (1.0 - a) + r[j + 1][c] * a;
            }
        }
    }
    r[0]
}

/// Cox-De Boor evaluation (2D, homogeneous 3-vector `(ux, uy, u)`).
pub fn de_boor_homo_2d(degree: usize, knots: &[f64], points: &[DVec2], weights: &[f64], t: f64) -> [f64; 3] {
    let n = points.len();
    if n == 0 { return [0.0; 3]; }
    let k = find_knot_span(degree, knots, t);
    let mut r = vec![[0.0f64; 3]; degree + 1];
    for j in 0..=degree {
        let idx = k - degree + j;
        let p = points[idx];
        let w = wgt(weights, idx);
        r[j] = [p.x * w, p.y * w, w];
    }
    for level in 1..=degree {
        for j in 0..=(degree - level) {
            let idx_a = k - degree + j + level;
            let denom = knots[idx_a + degree - level + 1] - knots[idx_a];
            let a = if denom.abs() > 1e-15 { ((t - knots[idx_a]) / denom).clamp(0.0, 1.0) } else { 0.0 };
            for c in 0..3 {
                r[j][c] = r[j][c] * (1.0 - a) + r[j + 1][c] * a;
            }
        }
    }
    r[0]
}

/// BSpline derivative via homogeneous quotient rule.
/// OCCT: BSplCLib::EvalDerivative.
pub fn bspline_tangent(degree: usize, knots: &[f64], points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n < 2 || degree == 0 { return DVec3::ZERO; }
    let p = degree as f64;
    let m = n - 1;
    let mut a_prime = Vec::with_capacity(m);
    let mut w_prime = vec![0.0f64; m];
    for i in 0..m {
        let denom = knots[i + degree + 1] - knots[i + 1];
        if denom.abs() < 1e-15 {
            a_prime.push(DVec3::ZERO);
        } else {
            let s = p / denom;
            let w_im1 = wgt(weights, i + 1);
            let w_i = wgt(weights, i);
            a_prime.push(s * (w_im1 * points[i + 1] - w_i * points[i]));
            w_prime[i] = s * (w_im1 - w_i);
        }
    }
    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0f64; m];
    let cp_prime = de_boor(degree - 1, deriv_knots, &a_prime, &unit, t);
    let w_val = de_boor(degree, knots, points, weights, t);
    let w_deriv = if w_prime.iter().any(|&w| w.abs() > 1e-15) {
        de_boor(degree - 1, deriv_knots, &from_vec_scalar(&w_prime), &unit, t).x
    } else { 0.0 };

    // Quotient rule: (C' W - C W') / W²
    let w0 = de_boor_homo(degree, knots, points, weights, t);
    let ww = w0[3];
    if ww.abs() > 1e-15 {
        let c_val = DVec3::new(w0[0], w0[1], w0[2]) / ww;
        (cp_prime - c_val * w_deriv) / ww
    } else {
        cp_prime
    }
}

/// BSpline derivative for 2D curves.
pub fn bspline_tangent_2d(degree: usize, knots: &[f64], points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n < 2 || degree == 0 { return DVec2::ZERO; }
    let p = degree as f64;
    let m = n - 1;
    let mut a_prime = Vec::with_capacity(m);
    let mut w_prime = vec![0.0f64; m];
    for i in 0..m {
        let denom = knots[i + degree + 1] - knots[i + 1];
        if denom.abs() < 1e-15 {
            a_prime.push(DVec2::ZERO);
        } else {
            let s = p / denom;
            let w_im1 = wgt(weights, i + 1);
            let w_i = wgt(weights, i);
            a_prime.push(s * (w_im1 * points[i + 1] - w_i * points[i]));
            w_prime[i] = s * (w_im1 - w_i);
        }
    }
    let deriv_knots = &knots[1..knots.len() - 1];
    let unit = vec![1.0f64; m];
    let cp_prime = de_boor_2d(degree - 1, deriv_knots, &a_prime, &unit, t);
    let w_deriv = if w_prime.iter().any(|&w| w.abs() > 1e-15) {
        de_boor_2d(degree - 1, deriv_knots, &from_vec_scalar_2d(&w_prime), &unit, t).x
    } else { 0.0 };
    let homo = de_boor_homo_2d(degree, knots, points, weights, t);
    let ww = homo[2];
    if ww.abs() > 1e-15 {
        let c_val = DVec2::new(homo[0], homo[1]) / ww;
        (cp_prime - c_val * w_deriv) / ww
    } else {
        cp_prime
    }
}

/// Helper: convert Vec<f64> to Vec<DVec3> with scalar.x = value.
fn from_vec_scalar(v: &[f64]) -> Vec<DVec3> {
    v.iter().map(|&x| DVec3::new(x, 0.0, 0.0)).collect()
}

/// Helper: convert Vec<f64> to Vec<DVec2> with scalar.x = value.
fn from_vec_scalar_2d(v: &[f64]) -> Vec<DVec2> {
    v.iter().map(|&x| DVec2::new(x, 0.0)).collect()
}

// ══════════════════════════════════════════════════════════════════════════
// BSplCLib — knot insertion and segmentation (Boehm)
// ══════════════════════════════════════════════════════════════════════════

/// Count the multiplicity of knot value `u` in the expanded knot vector.
/// OCCT: multiplicity of a knot in the Knots array.
fn knot_multiplicity(knots: &[f64], u: f64) -> usize {
    knots.iter().filter(|&&k| (k - u).abs() < 1e-12).count()
}

/// Distinct (compressed) knot values with their multiplicities.
/// OCCT: Knots() + Multiplicities().
fn compress_knots(knots: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut vals: Vec<f64> = Vec::new();
    let mut mults: Vec<usize> = Vec::new();
    for &k in knots {
        match vals.last() {
            Some(&last) if (k - last).abs() < 1e-15 => {
                *mults.last_mut().unwrap() += 1;
            }
            _ => {
                vals.push(k);
                mults.push(1);
            }
        }
    }
    (vals, mults)
}

/// OCCT BSplCLib::FirstUKnotIndex — index (0-based) into the compressed knot
/// array of the first knot whose cumulative multiplicity from the start
/// exceeds the degree.
///
/// For a clamped (non-periodic) curve the first boundary has multiplicity
/// Degree+1, so this returns 1 (the first unique knot). For a periodic curve
/// (boundary multiplicity == Degree) this returns 2 (skipping the seam knot).
pub fn first_uknot_index(knots: &[f64], degree: usize) -> usize {
    let (_, mults) = compress_knots(knots);
    let mut idx = 0;
    let mut sigma = mults[0];
    while sigma <= degree && idx + 1 < mults.len() {
        idx += 1;
        sigma += mults[idx];
    }
    idx
}

/// OCCT BSplCLib::LastUKnotIndex — index (0-based) into the compressed knot
/// array of the last knot whose cumulative multiplicity from the end exceeds
/// the degree.
pub fn last_uknot_index(knots: &[f64], degree: usize) -> usize {
    let (_, mults) = compress_knots(knots);
    let mut idx = mults.len() - 1;
    let mut sigma = mults[idx];
    while sigma <= degree && idx > 0 {
        idx -= 1;
        sigma += mults[idx];
    }
    idx
}

/// OCCT Geom_BSplineCurve::IsPeriodic() derived from the knot structure.
///
/// A periodic (unclamped) B-spline has first and last knot multiplicities
/// equal to the degree (not Degree+1), and the effective parameter range is
/// [Knots(FirstUKnotIndex), Knots(LastUKnotIndex)].
pub fn bspline_is_periodic(knots: &[f64], degree: usize) -> bool {
    let (_, mults) = compress_knots(knots);
    if mults.len() < 3 {
        return false;
    }
    mults[0] == degree && *mults.last().unwrap() == degree
}

/// Boehm knot insertion — insert knot `u` once into a rational B-spline.
///
/// Operates in homogeneous coordinates `(P*w, w)` (NURBS Book Algorithm A5.1).
/// OCCT: BSplCLib::InsertKnot (single insertion).
fn insert_knot_once(
    degree: usize,
    knots: &[f64],
    wpts: &[DVec3],
    wts: &[f64],
    u: f64,
) -> (Vec<f64>, Vec<DVec3>, Vec<f64>) {
    let k = find_knot_span(degree, knots, u);
    let p = degree;

    // New knots: insert u after index k.
    let mut new_knots = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.push(u);
    new_knots.extend_from_slice(&knots[k + 1..]);

    let n = wpts.len();
    let mut new_wpts = Vec::with_capacity(n + 1);
    let mut new_wts = Vec::with_capacity(n + 1);
    for i in 0..=n {
        if i <= k - p {
            new_wpts.push(wpts[i]);
            new_wts.push(wts[i]);
        } else if i >= k + 1 {
            new_wpts.push(wpts[i - 1]);
            new_wts.push(wts[i - 1]);
        } else {
            // k - p + 1 <= i <= k
            let denom = knots[i + p] - knots[i];
            let alpha = if denom.abs() > 1e-15 {
                ((u - knots[i]) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            new_wpts.push(alpha * wpts[i] + (1.0 - alpha) * wpts[i - 1]);
            new_wts.push(alpha * wts[i] + (1.0 - alpha) * wts[i - 1]);
        }
    }
    (new_knots, new_wpts, new_wts)
}

/// OCCT Geom_BSplineCurve::Segment(U1, U2, Tol) — extract the sub-curve on
/// [u1, u2] via knot insertion.
///
/// Inserts u1 and u2 to multiplicity Degree+1 (clamping), then extracts the
/// sub-curve between them. `tol` is used to snap u1/u2 to an existing knot
/// within tolerance (mirroring OCCT's tolerance handling in Segment).
pub fn segment_bspline_curve(
    b: &crate::geom::BSplineCurve3,
    u1: f64,
    u2: f64,
    tol: f64,
) -> Option<crate::geom::BSplineCurve3> {
    use crate::geom::BSplineCurve3;
    let degree = b.degree;
    let n = b.knots.len();
    if degree == 0 || n < 2 * degree + 2 {
        return None;
    }
    // Effective parameter domain via the UKnot indices — correct for both
    // clamped (boundary mult == Degree+1) and periodic (boundary mult == Degree)
    // knot vectors (OCCT BSplCLib::FirstUKnotIndex/LastUKnotIndex).
    let first = {
        let (vals, _) = compress_knots(&b.knots);
        let idx = first_uknot_index(&b.knots, degree);
        vals[idx]
    };
    let last = {
        let (vals, _) = compress_knots(&b.knots);
        let idx = last_uknot_index(&b.knots, degree);
        vals[idx]
    };

    // OCCT: snap to an existing knot within tolerance, then clamp to bounds.
    let snap = |v: f64| -> f64 {
        let mut r = v;
        for &k in &b.knots {
            if (k - v).abs() <= tol {
                r = k;
                break;
            }
        }
        r
    };
    let mut a = snap(u1);
    let mut c = snap(u2);
    if a < first {
        a = first;
    }
    if c > last {
        c = last;
    }
    if a >= c {
        return None;
    }

    // Work in homogeneous coordinates.
    let mut knots = b.knots.clone();
    let mut wpts: Vec<DVec3> = b
        .control_points
        .iter()
        .zip(&b.weights)
        .map(|(p, w)| *p * *w)
        .collect();
    let mut wts: Vec<f64> = b.weights.clone();

    // Clamp a to multiplicity Degree+1.
    let target = degree + 1;
    while knot_multiplicity(&knots, a) < target {
        (knots, wpts, wts) = insert_knot_once(degree, &knots, &wpts, &wts, a);
    }
    // Clamp c to multiplicity Degree+1.
    while knot_multiplicity(&knots, c) < target {
        (knots, wpts, wts) = insert_knot_once(degree, &knots, &wpts, &wts, c);
    }

    // Extract [a, c]: a occupies [i1, i1+degree], c occupies [i2-degree, i2].
    let i1 = knots.iter().position(|&k| (k - a).abs() < 1e-12)?;
    let i2 = knots.iter().rposition(|&k| (k - c).abs() < 1e-12)?;
    if i2 <= i1 {
        return None;
    }

    let sub_knots = knots[i1..=i2].to_vec();
    // sub-poles: wpts[i1 .. i2-degree] (exclusive end)
    let sub_wpts = &wpts[i1..i2 - degree];
    let sub_wts = &wts[i1..i2 - degree];
    let sub_pts: Vec<DVec3> = sub_wpts
        .iter()
        .zip(sub_wts)
        .map(|(wp, w)| if w.abs() > 1e-15 { *wp / *w } else { *wp })
        .collect();

    Some(BSplineCurve3 {
        degree,
        knots: sub_knots,
        control_points: sub_pts,
        weights: sub_wts.to_vec(),
        is_periodic: false,
    })
}

// ══════════════════════════════════════════════════════════════════════════
// BSplSLib::Resolution — max derivative bound and parametric resolution
// ══════════════════════════════════════════════════════════════════════════

/// OCCT `Standard_Real::Epsilon` (Standard_Real.hxx L242-246): the absolute
/// difference between `x` and the next representable double of the same sign
/// (one ULP). Used as the weight-variation threshold in `Rational()`.
#[inline]
fn occt_epsilon(x: f64) -> f64 {
    if x >= 0.0 {
        x.next_up() - x
    } else {
        x - x.next_down()
    }
}

/// OCCT `Rational()` weight-variation detection (Geom_BSplineSurface.cxx
/// L110-138 / Geom_BezierSurface.cxx L443-469). Returns `(URational,
/// VRational)` with OCCT's naming: `URational` is set when the weights vary
/// along the columns (the V direction), `VRational` when they vary along the
/// rows (the U direction). Note OCCT names are opposite to intuition — these
/// flags are passed straight into `BSplSLib::Resolution` where `URational`
/// guards the U-direction loop and `VRational` the V-direction loop.
pub(crate) fn surface_rational_flags(weights: &[Vec<f64>]) -> (bool, bool) {
    // OCCT L110-124: VRational = weights vary along rows (I, I+1, fixed J).
    let mut v_rational = false;
    'v_outer: for v in 0..weights[0].len() {
        for u in 0..weights.len() - 1 {
            let w = weights[u][v];
            if (w - weights[u + 1][v]).abs() > occt_epsilon(w.abs()) {
                v_rational = true;
                break 'v_outer;
            }
        }
    }
    // OCCT L126-137: URational = weights vary along columns (J, J+1, fixed I).
    let mut u_rational = false;
    'u_outer: for u in 0..weights.len() {
        for v in 0..weights[u].len() - 1 {
            let w = weights[u][v];
            if (w - weights[u][v + 1]).abs() > occt_epsilon(w.abs()) {
                u_rational = true;
                break 'u_outer;
            }
        }
    }
    (u_rational, v_rational)
}

/// OCCT `BSplSLib::Resolution` (BSplSLib.cxx L3519-3858): bounds the maximum
/// surface derivative over the knot spans and returns the parametric
/// resolution from a 3D tolerance:
///
///   U/VTolerance = Tolerance3D / (MaxDerivative * sqrt(2))
///
/// `poles` is the control point grid `[u][v]` and `weights` the weight grid
/// (only read when `u_rational`/`v_rational` is set; the rational branch also
/// divides by the minimum weight). `flat_knots_u`/`flat_knots_v` are the fully
/// expanded knot sequences (clamped, non-periodic — rcad surfaces carry no
/// periodic flag). `u_rational`/`v_rational` are the stored `myURational`/
/// `myVRational` flags from `surface_rational_flags`.
///
/// Indexing follows OCCT 1-based `NCollection_Array1/Array2` access: a flat
/// knot index `i` maps to `flat_knots[i - 1]`, a pole `Poles.Value(ii, jj)`
/// to `poles[ii - 1][jj - 1]`. All modulos over `poles_length` are replicated
/// verbatim (they are identities for a clamped grid, but OCCT applies them).
/// Returns `(utol, vtol)`, both 0.0 when a max derivative is 0.
pub fn bspl_surface_resolution(
    poles: &[Vec<DVec3>],
    weights: &[Vec<f64>],
    flat_knots_u: &[f64],
    flat_knots_v: &[f64],
    u_degree: usize,
    v_degree: usize,
    u_rational: bool,
    v_rational: bool,
    tol3d: f64,
) -> (f64, f64) {
    let mut max_derivative = [0.0f64; 2];

    let p_row_length = poles[0].len();
    let p_col_length = poles.len();

    // OCCT L3556-3571: min weight over the whole grid (only when rational).
    let mut min_weights = 0.0;
    if u_rational || v_rational {
        min_weights = weights[0][0];
        for u in 0..p_col_length {
            for v in 0..p_row_length {
                let w = weights[u][v];
                if w < min_weights {
                    min_weights = w;
                }
            }
        }
    }

    let ud1 = u_degree + 1;
    let vd1 = v_degree + 1;
    let num_poles_u = flat_knots_u.len() - ud1;
    let num_poles_v = flat_knots_v.len() - vd1;
    let poles_length_u = p_col_length;
    let poles_length_v = p_row_length;

    // OCCT L3578-3710: U direction (URational guards the U loop).
    if u_rational {
        let ud2 = u_degree * 2;
        let vd2 = v_degree * 2;
        for ii in 2..=num_poles_u {
            let ii_index = (ii - 1) % poles_length_u + 1;
            let ii_minus = (ii - 2) % poles_length_u + 1;
            let mut inverse = flat_knots_u[ii + u_degree - 1] - flat_knots_u[ii - 1];
            inverse = 1.0 / inverse;
            let lower0 = (ii as isize - ud1 as isize).max(1) as usize;
            let upper0 = (ii + ud2 + 1).min(num_poles_u);
            for jj in 1..=num_poles_v {
                let jj_index = (jj - 1) % poles_length_v + 1;
                let lower1 = (jj as isize - vd1 as isize).max(1) as usize;
                let upper1 = (jj + vd2 + 1).min(num_poles_v);
                let pij = poles[ii_index - 1][jj_index - 1];
                let wij = weights[ii_index - 1][jj_index - 1];
                let pmj = poles[ii_minus - 1][jj_index - 1];
                let wmj = weights[ii_minus - 1][jj_index - 1];
                let (xij, yij, zij) = (pij.x, pij.y, pij.z);
                let (xmj, ymj, zmj) = (pmj.x, pmj.y, pmj.z);
                for pp in lower0..=upper0 {
                    let pp_index = (pp - 1) % poles_length_u + 1;
                    for qq in lower1..=upper1 {
                        let qq_index = (qq - 1) % poles_length_v + 1;
                        let ppq = poles[pp_index - 1][qq_index - 1];
                        let (xpq, ypq, zpq) = (ppq.x, ppq.y, ppq.z);
                        let mut value = 0.0;
                        let mut factor = (xpq - xij) * wij;
                        factor -= (xpq - xmj) * wmj;
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                        factor = (ypq - yij) * wij;
                        factor -= (ypq - ymj) * wmj;
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                        factor = (zpq - zij) * wij;
                        factor -= (zpq - zmj) * wmj;
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                        value *= inverse;
                        if max_derivative[0] < value {
                            max_derivative[0] = value;
                        }
                    }
                }
            }
        }
        max_derivative[0] /= min_weights;
    } else {
        // OCCT L3671-3708: non-rational U direction (pure adjacent difference).
        for ii in 2..=num_poles_u {
            let ii_index = (ii - 1) % poles_length_u + 1;
            let ii_minus = (ii - 2) % poles_length_u + 1;
            let mut inverse = flat_knots_u[ii + u_degree - 1] - flat_knots_u[ii - 1];
            inverse = 1.0 / inverse;
            for jj in 1..=num_poles_v {
                let jj_index = (jj - 1) % poles_length_v + 1;
                let pij = poles[ii_index - 1][jj_index - 1];
                let pmj = poles[ii_minus - 1][jj_index - 1];
                let mut value = 0.0;
                let mut factor = pij.x - pmj.x;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pij.y - pmj.y;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pij.z - pmj.z;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                value *= inverse;
                if max_derivative[0] < value {
                    max_derivative[0] = value;
                }
            }
        }
    }
    max_derivative[0] *= u_degree as f64;

    // OCCT L3711-3843: V direction (VRational guards the V loop; ii is now
    // the V index, jj the U index).
    if v_rational {
        let ud2 = u_degree * 2;
        let vd2 = v_degree * 2;
        for ii in 2..=num_poles_v {
            let ii_index = (ii - 1) % poles_length_v + 1;
            let ii_minus = (ii - 2) % poles_length_v + 1;
            let mut inverse = flat_knots_v[ii + v_degree - 1] - flat_knots_v[ii - 1];
            inverse = 1.0 / inverse;
            let lower0 = (ii as isize - vd1 as isize).max(1) as usize;
            let upper0 = (ii + vd2 + 1).min(num_poles_v);
            for jj in 1..=num_poles_u {
                let jj_index = (jj - 1) % poles_length_u + 1;
                let lower1 = (jj as isize - ud1 as isize).max(1) as usize;
                let upper1 = (jj + ud2 + 1).min(num_poles_u);
                let pji = poles[jj_index - 1][ii_index - 1];
                let wji = weights[jj_index - 1][ii_index - 1];
                let pjm = poles[jj_index - 1][ii_minus - 1];
                let wjm = weights[jj_index - 1][ii_minus - 1];
                let (xji, yji, zji) = (pji.x, pji.y, pji.z);
                let (xjm, yjm, zjm) = (pjm.x, pjm.y, pjm.z);
                // OCCT L3757-3796: pp iterates the U window slice but the
                // pole is read as Value(qq_index, pp_index) with the modulos
                // crossed (L3759/L3764) — replicated verbatim.
                for pp in lower1..=upper1 {
                    let pp_index = (pp - 1) % poles_length_v + 1;
                    for qq in lower0..=upper0 {
                        let qq_index = (qq - 1) % poles_length_u + 1;
                        let pqp = poles[qq_index - 1][pp_index - 1];
                        let (xqp, yqp, zqp) = (pqp.x, pqp.y, pqp.z);
                        let mut value = 0.0;
                        let mut factor = (xqp - xji) * wji;
                        factor -= (xqp - xjm) * wjm;
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                        factor = (yqp - yji) * wji;
                        factor -= (yqp - yjm) * wjm;
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                        factor = (zqp - zji) * wji;
                        factor -= (zqp - zjm) * wjm;
                        if factor < 0.0 {
                            factor = -factor;
                        }
                        value += factor;
                        value *= inverse;
                        if max_derivative[1] < value {
                            max_derivative[1] = value;
                        }
                    }
                }
            }
        }
        max_derivative[1] /= min_weights;
    } else {
        // OCCT L3804-3841: non-rational V direction.
        for ii in 2..=num_poles_v {
            let ii_index = (ii - 1) % poles_length_v + 1;
            let ii_minus = (ii - 2) % poles_length_v + 1;
            let mut inverse = flat_knots_v[ii + v_degree - 1] - flat_knots_v[ii - 1];
            inverse = 1.0 / inverse;
            for jj in 1..=num_poles_u {
                let jj_index = (jj - 1) % poles_length_u + 1;
                let pji = poles[jj_index - 1][ii_index - 1];
                let pjm = poles[jj_index - 1][ii_minus - 1];
                let mut value = 0.0;
                let mut factor = pji.x - pjm.x;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pji.y - pjm.y;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = pji.z - pjm.z;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                value *= inverse;
                if max_derivative[1] < value {
                    max_derivative[1] = value;
                }
            }
        }
    }
    max_derivative[1] *= v_degree as f64;

    // OCCT L3844-3857.
    max_derivative[0] *= std::f64::consts::SQRT_2;
    max_derivative[1] *= std::f64::consts::SQRT_2;
    if max_derivative[0] != 0.0 && max_derivative[1] != 0.0 {
        (tol3d / max_derivative[0], tol3d / max_derivative[1])
    } else {
        (0.0, 0.0)
    }
}

/// OCCT `Geom_BSplineSurface::Resolution` (Geom_BSplineSurface_1.cxx
/// L2197-2222): `BSplSLib::Resolution` over the stored (already expanded)
/// knot vectors with the construction-time rational flags.
pub fn bspline_surface_resolution(s: &BSplineSurface, tol3d: f64) -> (f64, f64) {
    let u_deg = s.degree_u;
    let v_deg = s.degree_v;
    let (u_rational, v_rational) = surface_rational_flags(&s.weights);
    bspl_surface_resolution(
        &s.control_points,
        &s.weights,
        &s.knots_u,
        &s.knots_v,
        u_deg,
        v_deg,
        u_rational,
        v_rational,
        tol3d,
    )
}

/// OCCT `Geom_BezierSurface::Resolution` (Geom_BezierSurface.cxx L1991-2037):
/// `BSplSLib::Resolution` over the flat knot sequences `[0;deg+1] ++ [1;deg+1]`
/// in each direction (UKnots()/UMultiplicities() both pass the unit interval).
pub fn bezier_surface_resolution(s: &BezierSurface, tol3d: f64) -> (f64, f64) {
    let u_deg = s.control_points.len() - 1;
    let v_deg = s.control_points[0].len() - 1;
    let flat_u: Vec<f64> = std::iter::repeat_n(0.0, u_deg + 1)
        .chain(std::iter::repeat_n(1.0, u_deg + 1))
        .collect();
    let flat_v: Vec<f64> = std::iter::repeat_n(0.0, v_deg + 1)
        .chain(std::iter::repeat_n(1.0, v_deg + 1))
        .collect();
    let (u_rational, v_rational) = surface_rational_flags(&s.weights);
    bspl_surface_resolution(
        &s.control_points,
        &s.weights,
        &flat_u,
        &flat_v,
        u_deg,
        v_deg,
        u_rational,
        v_rational,
        tol3d,
    )
}

// ══════════════════════════════════════════════════════════════════════════
// 1D curve Resolution — GeomAdaptor_Curve::Resolution subchain
// (GeomAdaptor_Curve.cxx L1116-1148 -> Geom_BezierCurve/Geom_BSplineCurve
//  ::Resolution -> BSplCLib::Resolution 1D flat-array overload)
// ══════════════════════════════════════════════════════════════════════════

/// Rational() weight-variation detection for curves.
/// OCCT `Geom_BSplineCurve::Rational` (Geom_BSplineCurve.cxx L98-108) and
/// `Geom_BezierCurve::Rational` (Geom_BezierCurve.cxx L60-73): a curve is
/// rational iff some adjacent weight difference exceeds `gp::Resolution()`
/// which equals `RealSmall()` = DBL_MIN (gp.hxx L60) = `f64::MIN_POSITIVE`.
/// NOTE: this threshold differs from the surface `Rational()` detection
/// (1-ULP Epsilon in `surface_rational_flags`).
fn curve_rational(weights: &[f64]) -> bool {
    if weights.is_empty() {
        return false;
    }
    let mut current = weights[0];
    for &w in &weights[1..] {
        let delta = w - current;
        if delta.abs() > f64::MIN_POSITIVE {
            return true;
        }
        current = w;
    }
    false
}

/// OCCT `BSplCLib::PrepareUnperiodize` (BSplCLib.cxx L2967-3020): for a
/// periodic curve, count the knots/poles of the unperiodized curve.
/// `mults` is 0-based (rcad `compress_knots` output); OCCT's 1-based
/// `Mults(Lower()+k)` maps to `mults[k]`. Returns (nb_knots, nb_poles).
fn prepare_unperiodize(degree: usize, mults: &[usize]) -> (usize, usize) {
    // OCCT L2972-2980: NbKnots = Mults.Length(), NbPoles = -Degree - 1 + sum(Mults).
    let mut nb_knots = mults.len();
    let mut nb_poles = mults.iter().sum::<usize>() - (degree + 1);

    // OCCT L2983-3000: add knots at the beginning to raise the multiplicities
    // to Degree + 1.  k starts at Mults.Upper() - 1 (1-based) = 0-based len - 2.
    let mut sigma = mults[0];
    let mut k = mults.len() - 2;
    while sigma < degree + 1 {
        sigma += mults[k];
        nb_poles += mults[k];
        k -= 1;
        nb_knots += 1;
    }
    if sigma > degree + 1 {
        nb_poles -= sigma - degree - 1;
    }

    // OCCT L3002-3018: add knots at the end to raise the multiplicities to
    // Degree + 1.  k starts at Mults.Lower() + 1 (1-based = 2) = 0-based 1.
    let mut sigma = mults[mults.len() - 1];
    let mut k = 1;
    while sigma < degree + 1 {
        sigma += mults[k];
        nb_poles += mults[k];
        k += 1;
        nb_knots += 1;
    }
    if sigma > degree + 1 {
        nb_poles -= sigma - degree - 1;
    }

    (nb_knots, nb_poles)
}

/// OCCT `BSplCLib::Resolution` 1D flat-array overload (BSplCLib.cxx L4316-4820),
/// ArrayDimension = 3 case (BSplCLib_CurveComputation.pxx L1947-1965).
/// `flat_knots` is the fully expanded knot vector (rcad `BSplineCurve3.knots`,
/// or `[0;deg+1] ++ [1;deg+1]` for Bezier); `weights` mirrors the OCCT
/// `Weights()` accessor — NULL (non-rational) maps to `None`.
/// The flat `Poles` double array maps to `[DVec3]` (semantic, bottom-level).
/// OCCT tail (L4811-4820): NO `sqrt(2)` factor, guard `RealSmall()` = DBL_MIN.
pub fn bspl_curve_resolution(
    poles: &[DVec3],
    weights: Option<&[f64]>,
    flat_knots: &[f64],
    degree: usize,
    tol3d: f64,
) -> f64 {
    let deg1 = degree + 1;
    let deg2 = (degree << 1) + 1;
    let num_poles = flat_knots.len() - deg1; // local, from FlatKnots.Length()
    let num_poles_count = poles.len(); // the NumPoles argument (== num_poles for clamped)
    let mut max_derivative = 0.0f64;

    if let Some(wg) = weights {
        // OCCT L4446-4459: rational branch — minimum weight scan.
        let mut min_weights = wg[0];
        for &w in &wg[1..num_poles_count] {
            if w < min_weights {
                min_weights = w;
            }
        }
        // OCCT L4461-4527: windowed derivative estimate.
        for ii in 1..num_poles {
            let ii_index = ii % num_poles_count;
            let ii_minus = (ii - 1) % num_poles_count;
            let p_ii = poles[ii_index];
            let p_mi = poles[ii_minus];
            let wg_ii_index = wg[ii_index];
            let wg_ii_minus = wg[ii_minus];
            let mut inverse = flat_knots[ii + degree] - flat_knots[ii];
            inverse = 1.0 / inverse;
            let lower = if ii >= deg1 { ii - deg1 } else { 0 };
            let upper = if deg2 + ii > num_poles { num_poles } else { deg2 + ii };
            for jj in lower..upper {
                let p_jj = poles[jj % num_poles_count];
                let mut value = 0.0;
                let mut factor =
                    (p_jj.x - p_ii.x) * wg_ii_index - (p_jj.x - p_mi.x) * wg_ii_minus;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = (p_jj.y - p_ii.y) * wg_ii_index - (p_jj.y - p_mi.y) * wg_ii_minus;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                factor = (p_jj.z - p_ii.z) * wg_ii_index - (p_jj.z - p_mi.z) * wg_ii_minus;
                if factor < 0.0 {
                    factor = -factor;
                }
                value += factor;
                value *= inverse;
                if max_derivative < value {
                    max_derivative = value;
                }
            }
        }
        // OCCT L4527: divide by the minimum weight.
        max_derivative /= min_weights;
    } else {
        // OCCT L4529-4569: non-rational branch — adjacent pole differences.
        for ii in 1..num_poles {
            let ii_index = ii % num_poles_count;
            let ii_minus = (ii - 1) % num_poles_count;
            let p_ii = poles[ii_index];
            let p_mi = poles[ii_minus];
            let mut inverse = flat_knots[ii + degree] - flat_knots[ii];
            inverse = 1.0 / inverse;
            let mut value = 0.0;
            let mut factor = p_ii.x - p_mi.x;
            if factor < 0.0 {
                factor = -factor;
            }
            value += factor;
            factor = p_ii.y - p_mi.y;
            if factor < 0.0 {
                factor = -factor;
            }
            value += factor;
            factor = p_ii.z - p_mi.z;
            if factor < 0.0 {
                factor = -factor;
            }
            value += factor;
            value *= inverse;
            if max_derivative < value {
                max_derivative = value;
            }
        }
    }

    // OCCT L4811-4820: tail.
    max_derivative *= degree as f64;
    if max_derivative > f64::MIN_POSITIVE {
        tol3d / max_derivative
    } else {
        tol3d / f64::MIN_POSITIVE
    }
}

/// OCCT `Geom_BezierCurve::Resolution` (Geom_BezierCurve.cxx L743-758):
/// `BSplCLib::Resolution` over `KnotSequence()` = FlatBezierKnots(deg) =
/// `[0;deg+1] ++ [1;deg+1]` (L857-870), with `Weights()` (nullptr if
/// non-rational) and `Tolerance3D = 1.`; result scaled by the caller.
pub fn bezier_curve_resolution(c: &BezierCurve3, tol3d: f64) -> f64 {
    let degree = c.control_points.len() - 1;
    let flat_knots: Vec<f64> = std::iter::repeat_n(0.0, degree + 1)
        .chain(std::iter::repeat_n(1.0, degree + 1))
        .collect();
    let rational = curve_rational(&c.weights);
    let weights = if rational { Some(&c.weights[..]) } else { None };
    bspl_curve_resolution(&c.control_points, weights, &flat_knots, degree, tol3d)
}

/// OCCT `Geom_BSplineCurve::Resolution` (Geom_BSplineCurve_1.cxx L756-798):
/// periodic branch unperiodizes the poles (PrepareUnperiodize + modulo wrap)
/// and keeps `myFlatKnots`; non-periodic passes the poles directly.  Both pass
/// `Weights()` (nullptr if non-rational) and `Tolerance3D = 1.`.  OCCT builds
/// `new_weights` in the periodic branch but the `Resolution` call passes
/// `Weights()` (the original array), so only `new_poles` is needed in rcad.
pub fn bspline_curve_resolution(c: &BSplineCurve3, tol3d: f64) -> f64 {
    let rational = curve_rational(&c.weights);
    let weights = if rational { Some(&c.weights[..]) } else { None };
    if c.is_periodic {
        let (_, mults) = compress_knots(&c.knots);
        let (_nb_knots, nb_poles) = prepare_unperiodize(c.degree, &mults);
        let mut new_poles = Vec::with_capacity(nb_poles);
        for ii in 1..=nb_poles {
            new_poles.push(c.control_points[(ii - 1) % c.control_points.len()]);
        }
        bspl_curve_resolution(&new_poles, weights, &c.knots, c.degree, tol3d)
    } else {
        bspl_curve_resolution(&c.control_points, weights, &c.knots, c.degree, tol3d)
    }
}

#[cfg(test)]
mod tests {
    use glam::DVec3;
    use crate::geom::{BSplineCurve3, CurveEval};

    fn sample_cubic() -> BSplineCurve3 {
        // Clamped cubic BSpline, 6 control points -> 10 knots (degree + np + 1).
        BSplineCurve3 {
            degree: 3,
            knots: vec![0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0, 1.0],
            control_points: vec![
                DVec3::ZERO,
                DVec3::new(1.0, 0.0, 0.0),
                DVec3::new(2.0, 1.0, 0.0),
                DVec3::new(3.0, 1.0, 0.0),
                DVec3::new(4.0, 0.0, 0.0),
                DVec3::new(5.0, 0.0, 0.0),
            ],
            weights: vec![1.0; 6],
            is_periodic: false,
        }
    }

    #[test]
    fn segment_matches_original_on_subrange() {
        let c = sample_cubic();
        let seg = super::segment_bspline_curve(&c, 0.2, 0.8, 1e-9).expect("segment");
        assert_eq!(seg.degree, 3);
        // Segmented curve must coincide with the original on [0.2, 0.8].
        for i in 0..=20 {
            let t = 0.2 + (0.8 - 0.2) * (i as f64) / 20.0;
            let p_orig = c.point_at(t);
            let p_seg = seg.point_at(t);
            assert!(
                (p_orig - p_seg).length() < 1e-6,
                "t={t}: orig={p_orig:?} seg={p_seg:?}"
            );
        }
        // Endpoints must match exactly.
        assert!((seg.point_at(0.2) - c.point_at(0.2)).length() < 1e-9);
        assert!((seg.point_at(0.8) - c.point_at(0.8)).length() < 1e-9);
    }

    #[test]
    fn segment_full_range_is_identity() {
        let c = sample_cubic();
        let seg = super::segment_bspline_curve(&c, 0.0, 1.0, 1e-9).expect("segment");
        assert_eq!(seg.control_points.len(), c.control_points.len());
        assert!((seg.point_at(0.5) - c.point_at(0.5)).length() < 1e-9);
    }

    #[test]
    fn uknot_index_clamped_vs_periodic() {
        // Clamped cubic: boundary mult = 4 = degree+1.
        let clamped = vec![0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0];
        assert_eq!(super::first_uknot_index(&clamped, 3), 0);
        assert_eq!(super::last_uknot_index(&clamped, 3), 4);
        assert!(!super::bspline_is_periodic(&clamped, 3));

        // Periodic cubic: boundary mult = 3 = degree.
        let periodic = vec![0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0];
        assert_eq!(super::first_uknot_index(&periodic, 3), 1);
        assert_eq!(super::last_uknot_index(&periodic, 3), 3);
        assert!(super::bspline_is_periodic(&periodic, 3));
    }
}
