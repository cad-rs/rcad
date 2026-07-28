//! OCCT MathPoly: polynomial solvers.
//!
//! Corresponds to OCCT `math_DirectPolynomialRoots` and `math_Laguerre`.
//!
//! - `math_DirectPolynomialRoots` — solve_linear, solve_quadratic, solve_cubic, solve_quartic
//! - `math_Laguerre` — laguerre_roots (general polynomial root finding)
//! - `MathUtils_Poly` — poly_eval (Horner's method)

const TOL_FLOAT_DEDUP: f64 = 1e-15;
const TOL_LEN_MIN: f64 = 1e-12;
const TOL_LINEAR_ULTRA_STRICT: f64 = 1e-10;

// =============================================================================
// math_DirectPolynomialRoots
// =============================================================================

/// Solve linear equation ax + b = 0
pub fn solve_linear(a: f64, b: f64) -> Option<f64> {
    if a.abs() < TOL_FLOAT_DEDUP {
        if b.abs() < TOL_FLOAT_DEDUP { Some(0.0) } else { None }
    } else {
        Some(-b / a)
    }
}

/// Solve quadratic equation ax^2 + bx + c = 0
///
/// Returns real roots in ascending order.
pub fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < TOL_FLOAT_DEDUP {
        return solve_linear(b, c).into_iter().collect();
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return Vec::new();
    }
    if disc.abs() < TOL_FLOAT_DEDUP {
        return vec![-b / (2.0 * a)];
    }
    let sqrt_disc = disc.sqrt();
    let q = if b >= 0.0 {
        -0.5 * (b + sqrt_disc)
    } else {
        -0.5 * (b - sqrt_disc)
    };
    let mut roots = vec![q / a, c / q];
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

/// Real cube root
fn cube_root(x: f64) -> f64 {
    if x >= 0.0 { x.cbrt() } else { -(-x).cbrt() }
}

/// Solve cubic equation ax^3 + bx^2 + cx + d = 0
///
/// Uses Cardano's formula. Returns real roots in ascending order.
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < TOL_FLOAT_DEDUP {
        return solve_quadratic(b, c, d);
    }
    // Normalize to x^3 + px^2 + qx + r = 0
    let p = b / a;
    let q = c / a;
    let r = d / a;
    // Substitute x = t - p/3 to get t^3 + at + b = 0
    let a_coef = q - p * p / 3.0;
    let b_coef = 2.0 * p * p * p / 27.0 - p * q / 3.0 + r;
    let disc = b_coef * b_coef / 4.0 + a_coef * a_coef * a_coef / 27.0;
    let offset = p / 3.0;

    if disc.abs() < TOL_FLOAT_DEDUP {
        if a_coef.abs() < TOL_FLOAT_DEDUP {
            return vec![-offset];
        }
        let t = 3.0 * b_coef / a_coef;
        let mut roots = vec![-offset + t, -offset - t / 2.0];
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots
    } else if disc > 0.0 {
        let sqrt_disc = disc.sqrt();
        let u = cube_root(-b_coef / 2.0 + sqrt_disc);
        let v = cube_root(-b_coef / 2.0 - sqrt_disc);
        vec![u + v - offset]
    } else {
        let m = 2.0 * (-a_coef / 3.0).sqrt();
        let theta = (-b_coef / 2.0) / ((-a_coef * a_coef * a_coef / 27.0).sqrt());
        let theta = theta.clamp(-1.0, 1.0);
        let theta = theta.acos() / 3.0;
        let mut roots = vec![
            m * theta.cos() - offset,
            m * (theta + std::f64::consts::TAU / 3.0).cos() - offset,
            m * (theta + 2.0 * std::f64::consts::TAU / 3.0).cos() - offset,
        ];
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots
    }
}

/// Solve quartic equation ax^4 + bx^3 + cx^2 + dx + e = 0
///
/// Uses Ferrari's method. Returns real roots in ascending order.
pub fn solve_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Vec<f64> {
    if a.abs() < TOL_FLOAT_DEDUP {
        return solve_cubic(b, c, d, e);
    }

    let p = b / a;
    let q = c / a;
    let r = d / a;
    let s = e / a;

    let a1 = q - 3.0 * p * p / 8.0;
    let b1 = r + p * p * p / 8.0 - p * q / 2.0;
    let c1 = s - 3.0 * p * p * p * p / 256.0 + p * p * q / 16.0 - p * r / 4.0;

    if b1.abs() < TOL_LEN_MIN {
        let disc = a1 * a1 - 4.0 * c1;
        if disc < -TOL_LINEAR_ULTRA_STRICT {
            return Vec::new();
        }
        if disc.abs() < TOL_LINEAR_ULTRA_STRICT {
            let y = (-a1 / 2.0).sqrt();
            return vec![y - p / 4.0, -y - p / 4.0];
        }
        let sqrt_disc = disc.sqrt();
        let y1_sq = (-a1 + sqrt_disc) / 2.0;
        let y2_sq = (-a1 - sqrt_disc) / 2.0;

        let mut roots = Vec::new();
        if y1_sq >= -TOL_LINEAR_ULTRA_STRICT {
            let y1 = y1_sq.max(0.0).sqrt();
            roots.push(y1 - p / 4.0);
            roots.push(-y1 - p / 4.0);
        }
        if y2_sq >= -TOL_LINEAR_ULTRA_STRICT {
            let y2 = y2_sq.max(0.0).sqrt();
            roots.push(y2 - p / 4.0);
            roots.push(-y2 - p / 4.0);
        }
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return roots;
    }

    let resolvent_roots = solve_cubic(1.0, 2.0 * a1, a1 * a1 - 4.0 * c1, -b1 * b1);
    if resolvent_roots.is_empty() {
        return Vec::new();
    }

    let t = resolvent_roots
        .iter()
        .find(|&&t| t > TOL_LINEAR_ULTRA_STRICT)
        .copied()
        .unwrap_or(resolvent_roots[0]);

    let sqrt_t = t.max(0.0).sqrt();
    let mut roots = Vec::new();

    if sqrt_t > TOL_LINEAR_ULTRA_STRICT {
        let inner1 = -(a1 + t + b1 / sqrt_t);
        let inner2 = -(a1 + t - b1 / sqrt_t);

        if inner1 >= -TOL_LINEAR_ULTRA_STRICT {
            let s1 = inner1.max(0.0).sqrt();
            roots.push((sqrt_t + s1) / 2.0 - p / 4.0);
            roots.push((sqrt_t - s1) / 2.0 - p / 4.0);
        }
        if inner2 >= -TOL_LINEAR_ULTRA_STRICT {
            let s2 = inner2.max(0.0).sqrt();
            roots.push((-sqrt_t + s2) / 2.0 - p / 4.0);
            roots.push((-sqrt_t - s2) / 2.0 - p / 4.0);
        }
    } else {
        let inner = -(a1 + t);
        if inner >= -TOL_LINEAR_ULTRA_STRICT {
            let s = inner.max(0.0).sqrt();
            roots.push(s / 2.0 - p / 4.0);
            roots.push(-s / 2.0 - p / 4.0);
        }
    }

    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots.dedup_by(|a, b| (*a - *b).abs() < TOL_LINEAR_ULTRA_STRICT);
    roots
}

// =============================================================================
// math_Laguerre — general polynomial root finding via Laguerre's method
// =============================================================================

/// Find all real roots of polynomial a₀ + a₁x + ... + a_nx^n using Laguerre's method.
pub fn laguerre_roots(coeffs: &[f64]) -> Vec<f64> {
    let mut a: Vec<f64> = coeffs.iter().copied().collect();
    while a.len() > 1 && a.last().map_or(false, |&c| c.abs() < TOL_FLOAT_DEDUP) {
        a.pop();
    }
    let n = a.len() - 1;
    if n == 0 {
        return vec![];
    }
    if n <= 4 {
        return match n {
            1 => vec![-a[0] / a[1]],
            2 => solve_quadratic(a[2], a[1], a[0]),
            3 => solve_cubic(a[3], a[2], a[1], a[0]),
            _ => solve_quartic(a[4], a[3], a[2], a[1], a[0]),
        };
    }

    let mut roots = Vec::new();
    let mut deg = n;
    while deg >= 2 {
        let mut x = 0.0;
        for _ in 0..200 {
            let mut p = vec![0.0; deg + 1];
            p[deg] = a[deg];
            for i in (0..deg).rev() {
                p[i] = a[i] + x * p[i + 1];
            }
            let fv = p[0];
            if fv.abs() < 1e-14 {
                break;
            }
            let mut p1 = vec![0.0; deg];
            for i in 0..deg {
                p1[i] = a[i + 1] * (i as f64 + 1.0);
            }
            let mut fp = 0.0;
            for i in (0..deg - 1).rev() {
                fp = fp * x + p1[i];
            }
            let mut fp2 = 0.0;
            for i in (0..deg - 2).rev() {
                fp2 = fp2 * x + p1[i + 1] * (i as f64 + 1.0);
            }
            let g = fp / fv;
            let h = g * g - fp2 / fv;
            let d_ = ((deg as f64 - 1.0) * (deg as f64 * h - g * g)).sqrt();
            let s = if (g + d_).abs() > (g - d_).abs() {
                deg as f64 / (g + d_)
            } else {
                deg as f64 / (g - d_)
            };
            if s.abs() < 1e-14 {
                break;
            }
            x -= s;
            if s.abs() < 1e-12 {
                break;
            }
        }
        roots.push(x);
        let mut b = vec![0.0; deg];
        b[deg - 1] = a[deg];
        for i in (0..deg - 1).rev() {
            b[i] = a[i + 1] + x * b[i + 1];
        }
        a = b;
        deg = a.len() - 1;
    }
    if deg == 1 {
        roots.push(-a[0] / a[1]);
    }
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

// =============================================================================
// MathUtils_Poly — polynomial evaluation (Horner's method)
// =============================================================================

/// Evaluate a polynomial using Horner's method.
pub fn poly_eval(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}
