//! OCCT PLib: Polynomial Library — polynomial evaluation and basis conversion.
//!
//! Evaluation in Bezier form (de Casteljau), power basis (Horner), and
//! conversion between polynomial bases used by BSpline algorithms.
//!
//! OCCT source: src/FoundationClasses/TKMath/PLib/PLib.cxx

use glam::{DVec2, DVec3};

// ══════════════════════════════════════════════════════════════════════════
// De Casteljau — Bezier curve evaluation
// ══════════════════════════════════════════════════════════════════════════

/// De Casteljau evaluation of a rational Bezier curve (3D).
/// Repeated linear interpolation of weighted control points.
/// OCCT: PLib::Eval — Bezier basis.
pub fn de_casteljau_3d(points: &[DVec3], weights: &[f64], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 { return DVec3::ZERO; }
    if n == 1 { return points[0]; }
    let mut w: Vec<f64> = weights.to_vec();
    let mut r: Vec<DVec3> = points.iter().zip(weights.iter())
        .map(|(p, &wi)| *p * wi).collect();
    let mut m = n;
    while m > 1 {
        for i in 0..m - 1 {
            r[i] = r[i] * (1.0 - t) + r[i + 1] * t;
            w[i] = w[i] * (1.0 - t) + w[i + 1] * t;
        }
        m -= 1;
    }
    if w[0].abs() > 1e-15 { r[0] / w[0] } else { r[0] }
}

/// De Casteljau evaluation of a rational Bezier curve (2D).
pub fn de_casteljau_2d(points: &[DVec2], weights: &[f64], t: f64) -> DVec2 {
    let n = points.len();
    if n == 0 { return DVec2::ZERO; }
    if n == 1 { return points[0]; }
    let mut w: Vec<f64> = weights.to_vec();
    let mut r: Vec<DVec2> = points.iter().zip(weights.iter())
        .map(|(p, &wi)| *p * wi).collect();
    let mut m = n;
    while m > 1 {
        for i in 0..m - 1 {
            r[i] = r[i] * (1.0 - t) + r[i + 1] * t;
            w[i] = w[i] * (1.0 - t) + w[i + 1] * t;
        }
        m -= 1;
    }
    if w[0].abs() > 1e-15 { r[0] / w[0] } else { r[0] }
}

/// Non-rational de Casteljau (weights all 1). Evaluates control polygon directly.
pub fn de_casteljau_linear(points: &[DVec3], t: f64) -> DVec3 {
    let n = points.len();
    if n == 0 { return DVec3::ZERO; }
    if n == 1 { return points[0]; }
    let mut p = points.to_vec();
    for _k in 0..n - 1 {
        for i in 0..n - 1 {
            if i + 1 < p.len() {
                p[i] = p[i] * (1.0 - t) + p[i + 1] * t;
            }
        }
    }
    p[0]
}

// ══════════════════════════════════════════════════════════════════════════
// Power basis — Horner evaluation
// ══════════════════════════════════════════════════════════════════════════

/// Horner evaluation: a₀ + a₁·t + a₂·t² + ... + aₙ·tⁿ
/// OCCT: PLib::EvalPolynomial.
pub fn eval_polynomial(coeffs: &[f64], t: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * t + c)
}

/// Evaluate polynomial and its derivative simultaneously.
/// Returns (value, derivative).
pub fn eval_polynomial_d1(coeffs: &[f64], t: f64) -> (f64, f64) {
    let n = coeffs.len();
    if n == 0 { return (0.0, 0.0); }
    let mut val = coeffs[n - 1];
    let mut deriv = 0.0;
    for i in (1..n).rev() {
        deriv = deriv * t + val;
        val = val * t + coeffs[i - 1];
    }
    (val, deriv)
}

// ══════════════════════════════════════════════════════════════════════════
// Basis conversion utilities
// ══════════════════════════════════════════════════════════════════════════

/// Compute binomial coefficient C(n, k).
pub fn binomial(n: usize, k: usize) -> f64 {
    if k > n { return 0.0; }
    let k = k.min(n - k);
    let mut r = 1.0f64;
    for i in 1..=k {
        r = r * (n - k + i) as f64 / i as f64;
    }
    r
}

/// Convert power basis coefficients to Bezier control points (1D).
/// OCCT: `PLib::CoefficientsPoles` (dim = 1).
/// Input: power coefficients [a₀, a₁, ..., aₙ] for Σ aᵢ·tⁱ on t∈[0,1].
/// Output: Bezier control points [c₀, ..., cₙ] for Σ cᵢ·Bⁿᵢ(t).
///
/// Uses the identity `t^j = Σ_{i=j..n} (C(i,j)/C(n,j)) · Bⁿᵢ(t)`, giving
/// `c_i = Σ_{j=0..i} a_j · C(i,j)/C(n,j)`.
pub fn power_to_bezier(coeffs: &[f64]) -> Vec<f64> {
    let n = coeffs.len() - 1;
    let mut bez = vec![0.0; n + 1];
    for i in 0..=n {
        let mut sum = 0.0;
        for j in 0..=i {
            sum += coeffs[j] * binomial(i, j) / binomial(n, j);
        }
        bez[i] = sum;
    }
    bez
}

// ══════════════════════════════════════════════════════════════════════════
// Tests
// ══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eval_horner_constant() {
        assert!((eval_polynomial(&[5.0], 10.0) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn eval_horner_linear() {
        assert!((eval_polynomial(&[3.0, 2.0], 4.0) - 11.0).abs() < 1e-12);
    }

    #[test]
    fn eval_horner_quadratic() {
        assert!((eval_polynomial(&[1.0, 2.0, 1.0], 3.0) - 16.0).abs() < 1e-12);
    }

    #[test]
    fn eval_d1_linear() {
        let (v, d) = eval_polynomial_d1(&[1.0, 2.0], 3.0); // 1 + 2x
        assert!((v - 7.0).abs() < 1e-12);
        assert!((d - 2.0).abs() < 1e-12);
    }

    #[test]
    fn binomial_basic() {
        assert!((binomial(5, 2) - 10.0).abs() < 1e-12);
        assert!((binomial(4, 0) - 1.0).abs() < 1e-12);
        assert!((binomial(4, 4) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn de_casteljau_line_midpoint() {
        let pts = vec![DVec3::ZERO, DVec3::new(2.0, 0.0, 0.0)];
        let w = vec![1.0, 1.0];
        let p = de_casteljau_3d(&pts, &w, 0.5);
        assert!((p - DVec3::new(1.0, 0.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn de_casteljau_linear_midpoint() {
        let pts = vec![DVec3::ZERO, DVec3::new(4.0, 0.0, 0.0)];
        let p = de_casteljau_linear(&pts, 0.5);
        assert!((p - DVec3::new(2.0, 0.0, 0.0)).length() < 1e-12);
    }

    #[test]
    fn power_to_bezier_linear() {
        // f(t) = 1 + 2t, on [0,1]
        // Bezier: c0 = 1, c1 = 3
        let bez = power_to_bezier(&[1.0, 2.0]);
        assert_eq!(bez.len(), 2);
        assert!((bez[0] - 1.0).abs() < 1e-12);
        assert!((bez[1] - 3.0).abs() < 1e-12);
    }
}
