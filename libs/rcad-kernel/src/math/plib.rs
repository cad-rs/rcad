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
// OCCT PLib.cxx multi-dimensional polynomial evaluation and Hermite
// interpolation (used by AdvApprox / Convert_CompPolynomialToPoles / the
// helix pipeline).
// ══════════════════════════════════════════════════════════════════════════

/// OCCT `PLib::EvalPolynomial(Par, DerivativeRequest, Degree, Dimension,
/// PolynomialCoeff, Results)` (PLib.cxx L945-1029).
///
/// Coefficients are stored degree-major: `[c0(dim), c1(dim), ...,
/// cDegree(dim)]`.  On return `results` holds `(1 + derivative_request) *
/// dimension` values laid out `[value(dim), d1(dim), ..., dN(dim)]`.
/// The fused Horner recursion matches the OCCT optimized `eval_poly1/2`
/// templates operation-for-operation.
pub fn eval_polynomial_flat(
    par: f64,
    derivative_request: i32,
    degree: i32,
    dimension: i32,
    polynomial_coeff: &[f64],
    results: &mut [f64],
) {
    let dim = dimension as usize;
    match derivative_request {
        1 => {
            let mut local0 = vec![0.0f64; dim];
            let mut local1 = vec![0.0f64; dim];
            local0.copy_from_slice(&polynomial_coeff[..dim]);
            let mut coeffs = 0usize;
            for _a_deg in 0..degree {
                coeffs += dim;
                for i in 0..dim {
                    let val = local0[i];
                    local1[i] = local1[i] * par + val;
                    local0[i] = val * par + polynomial_coeff[coeffs + i];
                }
            }
            results[..dim].copy_from_slice(&local0);
            results[dim..2 * dim].copy_from_slice(&local1);
        }
        2 => {
            let mut local0 = vec![0.0f64; dim];
            let mut local1 = vec![0.0f64; dim];
            let mut local2 = vec![0.0f64; dim];
            local0.copy_from_slice(&polynomial_coeff[..dim]);
            let mut coeffs = 0usize;
            for _a_deg in 0..degree {
                coeffs += dim;
                for i in 0..dim {
                    let d1 = local1[i];
                    let val = local0[i];
                    local2[i] = local2[i] * par + d1 * 2.0;
                    local1[i] = d1 * par + val;
                    local0[i] = val * par + polynomial_coeff[coeffs + i];
                }
            }
            results[..dim].copy_from_slice(&local0);
            results[dim..2 * dim].copy_from_slice(&local1);
            results[2 * dim..3 * dim].copy_from_slice(&local2);
        }
        _ => {
            // General case for DerivativeRequest > 2 (and 0).
            let res_size = (1 + derivative_request) as usize * dim;
            for v in results[..res_size].iter_mut() {
                *v = 0.0;
            }
            let mut a_coeffs = degree as usize * dim;
            for _a_deg in 0..=degree {
                // aPtr walks from the highest derivative slot down to slot 0.
                let mut a_ptr = res_size - dim;
                for a_deriv in (1..=derivative_request).rev() {
                    let an_original = a_ptr - dim;
                    for ind in 0..dim {
                        results[a_ptr + ind] =
                            results[a_ptr + ind] * par + results[an_original + ind] * a_deriv as f64;
                    }
                    a_ptr = an_original;
                }
                for ind in 0..dim {
                    results[a_ptr + ind] = results[a_ptr + ind] * par + polynomial_coeff[a_coeffs + ind];
                }
                if a_coeffs >= dim {
                    a_coeffs -= dim;
                }
            }
        }
    }
}

/// OCCT `PLib::NoDerivativeEvalPolynomial` (value only).  `degree_dimension`
/// is the offset of the HIGHEST-degree coefficient (`Deg * Dimension`); the
/// Horner walk descends with stride `dimension` (OCCT eval_poly0).
pub fn no_derivative_eval_polynomial_flat(
    par: f64,
    degree: i32,
    dimension: i32,
    degree_dimension: i32,
    polynomial_coeff: &[f64],
    results: &mut [f64],
) {
    let dim = dimension as usize;
    let mut coeffs = degree_dimension as usize;
    let mut local = vec![0.0f64; dim];
    local.copy_from_slice(&polynomial_coeff[coeffs..coeffs + dim]);
    for _a_deg in 0..degree {
        coeffs -= dim;
        for i in 0..dim {
            local[i] = local[i] * par + polynomial_coeff[coeffs + i];
        }
    }
    results[..dim].copy_from_slice(&local);
}

/// OCCT `PLib::HermiteInterpolate(Dimension, FirstParameter, LastParameter,
/// FirstOrder, LastOrder, FirstConstr, LastConstr, Coefficients)`
/// (PLib.cxx L1931-2027) — constrained Hermite interpolation on [-1, 1]
/// (here on the caller-supplied parameters), solved via `math_Gauss`.
/// `first_constr` / `last_constr` are (dimension x (order+1)) row-major
/// matrices addressed as (idim, order) 1-based.  Returns false when the
/// system is singular.
#[allow(clippy::too_many_arguments)]
pub fn hermite_interpolate(
    dimension: usize,
    first_parameter: f64,
    last_parameter: f64,
    first_order: usize,
    last_order: usize,
    first_constr: &crate::math::MatD,
    last_constr: &crate::math::MatD,
    coefficients: &mut [f64],
) -> bool {
    let pattern: [[f64; 6]; 3] = [
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        [0.0, 0.0, 2.0, 6.0, 12.0, 20.0],
    ];

    let n = first_order + last_order + 1; // matrix order (0-based extent)
    let mut a = crate::math::MatD::new(n + 1, n + 1);

    for irow in 0..=first_order {
        let mut first_val = 1.0f64;
        for icol in 0..=n {
            a.set(irow + 1, icol + 1, pattern[irow][icol] * first_val);
            if irow <= icol {
                first_val *= first_parameter;
            }
        }
    }
    for irow in 0..=last_order {
        let mut last_val = 1.0f64;
        for icol in 0..=n {
            a.set(irow + first_order + 2, icol + 1, pattern[irow][icol] * last_val);
            if irow <= icol {
                last_val *= last_parameter;
            }
        }
    }

    let equations = crate::math::math_gauss::MathGauss::new(&a);

    for idim in 1..=dimension {
        let mut b = crate::math::VecD::new(n + 1);
        for icol in 0..=first_order {
            b.set(icol + 1, first_constr.get(idim, icol + 1));
        }
        for icol in 0..=last_order {
            b.set(first_order + 2 + icol, last_constr.get(idim, icol + 1));
        }

        if !equations.is_done() {
            return false;
        }
        equations.solve(&mut b);

        for icol in 0..=n {
            coefficients[dimension * icol + idim - 1] = b.get(icol + 1);
        }
    }
    true
}

/// OCCT `PLib::JacobiParameters(ConstraintOrder, MaxDegree, Code,
/// NbGaussPoints, WorkDegree)` (PLib.cxx L2049-2192).
pub fn jacobi_parameters(
    constraint_order: crate::math::GeomAbsShape,
    max_degree: usize,
    code: i32,
    nb_gauss_points: &mut usize,
    work_degree: &mut usize,
) {
    const NDEG8: usize = 8;
    const NDEG10: usize = 10;
    const NDEG15: usize = 15;
    const NDEG20: usize = 20;
    const NDEG25: usize = 25;
    const NDEG30: usize = 30;
    const NDEG40: usize = 40;
    const NDEG50: usize = 50;
    const NDEG61: usize = 61;

    let niv_constr = crate::math::p_lib_jacobi::niv_constr(constraint_order);
    assert!(
        max_degree >= 2 * niv_constr + 1,
        "Invalid MaxDegree"
    );

    if code >= 1 {
        *work_degree = max_degree + 9;
    } else {
        *work_degree = max_degree + 6;
    }

    // Nbre mini de points nécessaires.
    let mut ipmin = 0usize;
    let wd = *work_degree;
    if wd < NDEG8 {
        ipmin = NDEG8;
    } else if wd < NDEG10 {
        ipmin = NDEG10;
    } else if wd < NDEG15 {
        ipmin = NDEG15;
    } else if wd < NDEG20 {
        ipmin = NDEG20;
    } else if wd < NDEG25 {
        ipmin = NDEG25;
    } else if wd < NDEG30 {
        ipmin = NDEG30;
    } else if wd < NDEG40 {
        ipmin = NDEG40;
    } else if wd < NDEG50 {
        ipmin = NDEG50;
    } else if wd < NDEG61 {
        ipmin = NDEG61;
    } else {
        panic!("Invalid MaxDegree");
    }

    // Nbre de points voulus.
    let iwant = match code {
        -5 => NDEG8,
        -4 => NDEG10,
        -3 => NDEG15,
        -2 => NDEG20,
        -1 => NDEG25,
        1 => NDEG30,
        2 => NDEG40,
        3 => NDEG50,
        4 => NDEG61,
        _ => panic!("Invalid Code"),
    };

    *nb_gauss_points = ipmin.max(iwant);
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

/// OCCT `PLib::EvalLagrange(Parameter, DerivativeRequest, Degree, Dimension,
/// Values, Parameters, Results)` (PLib.cxx L1122-1249) — evaluates the
/// Lagrange interpolation polynomial through `Degree + 1` points (each of
/// `Dimension` coordinates, stored interleaved in `values`, with the point
/// parameters in `parameters`) and its derivatives up to `DerivativeRequest`
/// at `parameter`.  `results` receives the value row followed by the
/// derivative rows (`(DerivativeRequest + 1) * Dimension` entries).  Returns
/// 1 when two parameters coincide (division by zero, results undefined).
pub fn eval_lagrange(
    parameter: f64,
    derivative_request: usize,
    degree: usize,
    dimension: usize,
    values: &[f64],
    parameters: &[f64],
    results: &mut [f64],
) -> i32 {
    let mut local_request = derivative_request;
    if local_request >= degree {
        local_request = degree;
    }

    // Build the divided differences array (copied from the point values).
    let mut divided_differences: Vec<f64> = vec![0.0; (degree + 1) * dimension];
    for (i, v) in values.iter().enumerate().take((degree + 1) * dimension) {
        divided_differences[i] = *v;
    }

    let mut return_code = 0;
    let mut ii = degree as i64;
    while ii >= 0 {
        let mut jj = degree as i64;
        while jj > degree as i64 - ii {
            let index = jj * dimension as i64;
            let index1 = index - dimension as i64;
            for kk in 0..dimension {
                divided_differences[(index + kk as i64) as usize] -=
                    divided_differences[(index1 + kk as i64) as usize];
            }
            let mut difference =
                parameters[jj as usize] - parameters[(jj - degree as i64 - 1 + ii) as usize];
            if difference.abs() < f64::MIN_POSITIVE {
                return_code = 1;
                ii = -1;
                break;
            }
            difference = 1.0e0 / difference;
            for kk in 0..dimension {
                divided_differences[(index + kk as i64) as usize] *= difference;
            }
            jj -= 1;
        }
        ii -= 1;
    }
    if return_code != 0 {
        return return_code;
    }

    // Evaluate the Newton form: P(t) = [t1]P + (t-t1)[t1,t2]P + ...
    let index = degree * dimension;
    for kk in 0..dimension {
        results[kk] = divided_differences[index + kk];
    }
    for (i, r) in results.iter_mut().enumerate().take((local_request + 1) * dimension) {
        if i >= dimension {
            *r = 0.0e0;
        }
    }

    let mut ii = degree as i64;
    while ii >= 1 {
        let difference = parameter - parameters[(ii - 1) as usize];
        let mut jj = local_request as i64;
        while jj > 0 {
            let index = (jj * dimension as i64) as usize;
            let index1 = index - dimension;
            for kk in 0..dimension {
                results[index + kk] *= difference;
                results[index + kk] += results[index1 + kk] * jj as f64;
            }
            jj -= 1;
        }
        let index = ((ii - 1) * dimension as i64) as usize;
        for kk in 0..dimension {
            results[kk] *= difference;
            results[kk] += divided_differences[index + kk];
        }
        ii -= 1;
    }
    return_code
}

/// OCCT PLib::CoefficientsPoles for gp_Pnt arrays (PLib.cxx L1482-1493),
/// dispatching to the dim = 3 version (L1522-1608).  Converts the power
/// coefficients of a degree `n = len - 1` polynomial into the Bezier poles.
/// Only the non-rational path (WCoefs == nullptr) is used by GeomFill_Coons.
pub fn coefficients_poles(coefs: &[DVec3]) -> Vec<DVec3> {
    let reflen = coefs.len();
    let mut poles = vec![DVec3::ZERO; reflen];
    // Les Extremites.
    poles[0] = coefs[0];
    poles[reflen - 1] = coefs[reflen - 1];

    for i in 2..reflen {
        let cnp = binomial(reflen - 1, i - 1);
        poles[i - 1] = coefs[i - 1] / cnp;
    }

    for i in 1..=reflen - 1 {
        for j in (i..reflen).rev() {
            let prev = poles[j - 1];
            poles[j] += prev;
        }
    }
    poles
}
