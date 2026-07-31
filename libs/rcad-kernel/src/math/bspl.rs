//! OCCT BSplCLib + BSplSLib: BSpline curve and surface evaluation.
//!
//! Core algorithms:
//! - Cox-De Boor evaluation (rational NURBS, 3D and 2D)
//! - BSpline derivative via homogeneous quotient rule
//! - Knot span binary search (same knot vector convention as OCCT)
//!
//! OCCT source: src/FoundationClasses/TKMath/BSplCLib/BSplCLib.cxx
//!             src/FoundationClasses/TKMath/BSplSLib/BSplSLib.cxx

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
        r[j] = points[idx] * weights[idx];
        w[j] = weights[idx];
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
        r[j] = points[idx] * weights[idx];
        w[j] = weights[idx];
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
        let w = weights[idx];
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
        let w = weights[idx];
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
            a_prime.push(s * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
            w_prime[i] = s * (weights[i + 1] - weights[i]);
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
            a_prime.push(s * (weights[i + 1] * points[i + 1] - weights[i] * points[i]));
            w_prime[i] = s * (weights[i + 1] - weights[i]);
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
