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
