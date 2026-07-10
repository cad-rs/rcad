//! OCCT math-style mathematical utilities.
//!
//! Provides mathematical algorithms and solvers including:
//! - Root finding (Newton-Raphson, bisection, secant)
//! - Multi-dimensional Newton methods
//! - Polynomial solvers (linear through quartic)
//! - Eigenvalue/matrix utilities
//! - Numerical integration
//! - Optimization (golden section search)

use glam::{DMat2, DMat3, DVec2, DVec3};
use std::f64::consts::FRAC_1_SQRT_2;
use crate::geom::{BSplineSurface, Surface3};

// Tolerance constants (local copies — moved from rcad-algorithms tolerance)
const TOLERANCE_FLOAT_DEDUP: f64 = 1e-15;
const TOLERANCE_LEN_MIN: f64 = 1e-12;
const TOLERANCE_LINEAR_ULTRA_STRICT: f64 = 1e-10;

// =============================================================================
// Root Finding
// =============================================================================

/// Newton-Raphson method for finding roots of a function.
///
/// # Arguments
/// * `f` - The function to find root of
/// * `df` - The derivative of the function
/// * `x0` - Initial guess
/// * `tol` - Tolerance for convergence
/// * `max_iter` - Maximum number of iterations
///
/// # Returns
/// The root if found within tolerance and iteration limit
pub fn newton_raphson(
    f: fn(f64) -> f64,
    df: fn(f64) -> f64,
    x0: f64,
    tol: f64,
    max_iter: usize,
) -> Option<f64> {
    let mut x = x0;

    for _ in 0..max_iter {
        let fx = f(x);
        if fx.abs() < tol {
            return Some(x);
        }

        let dfx = df(x);
        if dfx.abs() < TOLERANCE_FLOAT_DEDUP {
            return None; // Derivative too small
        }

        let x_new = x - fx / dfx;

        // Check for convergence
        if (x_new - x).abs() < tol {
            return Some(x_new);
        }

        x = x_new;
    }

    // Check final value
    if f(x).abs() < tol {
        Some(x)
    } else {
        None
    }
}

/// Bisection method for finding roots of a function.
///
/// # Arguments
/// * `f` - The function to find root of
/// * `a` - Lower bound of interval
/// * `b` - Upper bound of interval
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root if found within interval and tolerance
pub fn bisection(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> Option<f64> {
    let mut lo = a;
    let mut hi = b;

    let f_lo = f(lo);
    let f_hi = f(hi);

    // Check if interval brackets a root
    if f_lo * f_hi > 0.0 {
        return None;
    }

    // Check if bounds are already roots
    if f_lo.abs() < tol {
        return Some(lo);
    }
    if f_hi.abs() < tol {
        return Some(hi);
    }

    let max_iter = ((hi - lo) / tol).ceil() as usize;

    for _ in 0..max_iter {
        let mid = (lo + hi) / 2.0;
        let f_mid = f(mid);

        if f_mid.abs() < tol || (hi - lo) / 2.0 < tol {
            return Some(mid);
        }

        if f_lo * f_mid < 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }

    Some((lo + hi) / 2.0)
}

/// Secant method for finding roots of a function.
///
/// # Arguments
/// * `f` - The function to find root of
/// * `x0` - First initial guess
/// * `x1` - Second initial guess
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root if found within tolerance
pub fn secant(f: fn(f64) -> f64, x0: f64, x1: f64, tol: f64) -> Option<f64> {
    let mut x_prev = x0;
    let mut x_curr = x1;

    let max_iter = 100;

    for _ in 0..max_iter {
        let f_prev = f(x_prev);
        let f_curr = f(x_curr);

        if f_curr.abs() < tol {
            return Some(x_curr);
        }

        let denom = f_curr - f_prev;
        if denom.abs() < TOLERANCE_FLOAT_DEDUP {
            return None;
        }

        let x_new = x_curr - f_curr * (x_curr - x_prev) / denom;

        if (x_new - x_curr).abs() < tol {
            return Some(x_new);
        }

        x_prev = x_curr;
        x_curr = x_new;
    }

    if f(x_curr).abs() < tol {
        Some(x_curr)
    } else {
        None
    }
}

// =============================================================================
// Multi-dimensional Newton Methods
// =============================================================================

/// Newton-Raphson method for 2D systems of equations.
///
/// Solves the system F(x) = 0 where F: R^2 -> R^2
///
/// # Arguments
/// * `f` - The function vector F(x)
/// * `jacobian` - The Jacobian matrix of F
/// * `x0` - Initial guess
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root vector if found
pub fn newton_2d(
    f: fn(DVec2) -> DVec2,
    jacobian: fn(DVec2) -> DMat2,
    x0: DVec2,
    tol: f64,
) -> Option<DVec2> {
    let mut x = x0;
    let max_iter = 50;

    for _ in 0..max_iter {
        let fx = f(x);

        if fx.length() < tol {
            return Some(x);
        }

        let j = jacobian(x);
        let det = j.determinant();

        if det.abs() < TOLERANCE_FLOAT_DEDUP {
            return None; // Singular Jacobian
        }

        let j_inv = j.inverse();
        let delta = j_inv * fx;

        let x_new = x - delta;

        if delta.length() < tol {
            return Some(x_new);
        }

        x = x_new;
    }

    if f(x).length() < tol {
        Some(x)
    } else {
        None
    }
}

/// Newton-Raphson method for 3D systems of equations.
///
/// Solves the system F(x) = 0 where F: R^3 -> R^3
///
/// # Arguments
/// * `f` - The function vector F(x)
/// * `jacobian` - The Jacobian matrix of F
/// * `x0` - Initial guess
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The root vector if found
pub fn newton_3d(
    f: fn(DVec3) -> DVec3,
    jacobian: fn(DVec3) -> DMat3,
    x0: DVec3,
    tol: f64,
) -> Option<DVec3> {
    let mut x = x0;
    let max_iter = 50;

    for _ in 0..max_iter {
        let fx = f(x);

        if fx.length() < tol {
            return Some(x);
        }

        let j = jacobian(x);
        let det = j.determinant();

        if det.abs() < TOLERANCE_FLOAT_DEDUP {
            return None; // Singular Jacobian
        }

        if let Some(j_inv) = inverse_3x3(j) {
            let delta = j_inv * fx;
            let x_new = x - delta;

            if delta.length() < tol {
                return Some(x_new);
            }

            x = x_new;
        } else {
            return None;
        }
    }

    if f(x).length() < tol {
        Some(x)
    } else {
        None
    }
}

// =============================================================================
// Polynomial Solvers
// =============================================================================

/// Solve linear equation ax + b = 0
pub fn solve_linear(a: f64, b: f64) -> Option<f64> {
    if a.abs() < TOLERANCE_FLOAT_DEDUP {
        if b.abs() < TOLERANCE_FLOAT_DEDUP {
            Some(0.0) // Infinite solutions, return 0
        } else {
            None // No solution
        }
    } else {
        Some(-b / a)
    }
}

/// Solve quadratic equation ax^2 + bx + c = 0
///
/// Returns real roots in ascending order
pub fn solve_quadratic(a: f64, b: f64, c: f64) -> Vec<f64> {
    if a.abs() < TOLERANCE_FLOAT_DEDUP {
        // Linear case
        return solve_linear(b, c).into_iter().collect();
    }

    let disc = b * b - 4.0 * a * c;

    if disc < 0.0 {
        return Vec::new(); // No real roots
    }

    if disc.abs() < TOLERANCE_FLOAT_DEDUP {
        // Single root (double root)
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

/// Solve cubic equation ax^3 + bx^2 + cx + d = 0
///
/// Uses Cardano's formula and returns real roots in ascending order
pub fn solve_cubic(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    if a.abs() < TOLERANCE_FLOAT_DEDUP {
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

    if disc.abs() < TOLERANCE_FLOAT_DEDUP {
        // One or two roots (discriminant near zero)
        if a_coef.abs() < TOLERANCE_FLOAT_DEDUP {
            // Triple root
            return vec![-offset];
        }
        let t = 3.0 * b_coef / a_coef;
        let mut roots = vec![-offset + t, -offset - t / 2.0];
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        roots
    } else if disc > 0.0 {
        // One real root
        let sqrt_disc = disc.sqrt();
        let u = cube_root(-b_coef / 2.0 + sqrt_disc);
        let v = cube_root(-b_coef / 2.0 - sqrt_disc);
        vec![u + v - offset]
    } else {
        // Three real roots (trigonometric solution)
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

/// Real cube root
fn cube_root(x: f64) -> f64 {
    if x >= 0.0 {
        x.cbrt()
    } else {
        -(-x).cbrt()
    }
}

/// Solve quartic equation ax^4 + bx^3 + cx^2 + dx + e = 0
///
/// Uses Ferrari's method and returns real roots in ascending order
pub fn solve_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Vec<f64> {
    if a.abs() < TOLERANCE_FLOAT_DEDUP {
        return solve_cubic(b, c, d, e);
    }

    // Normalize to x^4 + px^3 + qx^2 + rx + s = 0
    let p = b / a;
    let q = c / a;
    let r = d / a;
    let s = e / a;

    // Substitute x = y - p/4 to get depressed quartic y^4 + py^2 + qy + r = 0
    let a1 = q - 3.0 * p * p / 8.0;
    let b1 = r + p * p * p / 8.0 - p * q / 2.0;
    let c1 = s - 3.0 * p * p * p * p / 256.0 + p * p * q / 16.0 - p * r / 4.0;

    // Handle b1 = 0 case (quartic with only even powers)
    if b1.abs() < TOLERANCE_LEN_MIN {
        let disc = a1 * a1 - 4.0 * c1;
        if disc < -TOLERANCE_LINEAR_ULTRA_STRICT {
            return Vec::new();
        }
        if disc.abs() < TOLERANCE_LINEAR_ULTRA_STRICT {
            let y = (-a1 / 2.0).sqrt();
            return vec![y - p / 4.0, -y - p / 4.0];
        }
        let sqrt_disc = disc.sqrt();
        let y1_sq = (-a1 + sqrt_disc) / 2.0;
        let y2_sq = (-a1 - sqrt_disc) / 2.0;

        let mut roots = Vec::new();
        if y1_sq >= -TOLERANCE_LINEAR_ULTRA_STRICT {
            let y1 = y1_sq.max(0.0).sqrt();
            roots.push(y1 - p / 4.0);
            roots.push(-y1 - p / 4.0);
        }
        if y2_sq >= -TOLERANCE_LINEAR_ULTRA_STRICT {
            let y2 = y2_sq.max(0.0).sqrt();
            roots.push(y2 - p / 4.0);
            roots.push(-y2 - p / 4.0);
        }
        roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
        return roots;
    }

    // Solve resolvent cubic: t^3 + 2*a1*t^2 + (a1^2 - 4*c1)*t - b1^2 = 0
    let resolvent_roots = solve_cubic(1.0, 2.0 * a1, a1 * a1 - 4.0 * c1, -b1 * b1);

    if resolvent_roots.is_empty() {
        return Vec::new();
    }

    // Find a positive root of the resolvent
    let t = resolvent_roots
        .iter()
        .find(|&&t| t > TOLERANCE_LINEAR_ULTRA_STRICT)
        .copied()
        .unwrap_or(resolvent_roots[0]);

    let sqrt_t = t.max(0.0).sqrt();

    let mut roots = Vec::new();

    if sqrt_t > TOLERANCE_LINEAR_ULTRA_STRICT {
        let inner1 = -(a1 + t + b1 / sqrt_t);
        let inner2 = -(a1 + t - b1 / sqrt_t);

        if inner1 >= -TOLERANCE_LINEAR_ULTRA_STRICT {
            let s1 = inner1.max(0.0).sqrt();
            roots.push((sqrt_t + s1) / 2.0 - p / 4.0);
            roots.push((sqrt_t - s1) / 2.0 - p / 4.0);
        }
        if inner2 >= -TOLERANCE_LINEAR_ULTRA_STRICT {
            let s2 = inner2.max(0.0).sqrt();
            roots.push((-sqrt_t + s2) / 2.0 - p / 4.0);
            roots.push((-sqrt_t - s2) / 2.0 - p / 4.0);
        }
    } else {
        // t is nearly zero, use alternative formula
        let inner = -(a1 + t);
        if inner >= -TOLERANCE_LINEAR_ULTRA_STRICT {
            let s = inner.max(0.0).sqrt();
            roots.push(s / 2.0 - p / 4.0);
            roots.push(-s / 2.0 - p / 4.0);
        }
    }

    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots.dedup_by(|a, b| (*a - *b).abs() < TOLERANCE_LINEAR_ULTRA_STRICT);
    roots
}

// =============================================================================
// Eigenvalue/Matrix Utilities
// =============================================================================

/// Compute eigenvalues of a 2x2 matrix.
///
/// Returns eigenvalues in descending order (largest first)
pub fn eigenvalues_2x2(m: DMat2) -> (f64, f64) {
    // For matrix [[a, b], [c, d]]:
    // eigenvalues satisfy: lambda^2 - (a+d)*lambda + (ad-bc) = 0
    let trace = m.x_axis.x + m.y_axis.y; // trace = sum of diagonal
    let det = m.determinant();

    let disc = trace * trace - 4.0 * det;

    if disc < 0.0 {
        // Complex eigenvalues - return real parts
        let real_part = trace / 2.0;
        (real_part, real_part)
    } else {
        let sqrt_disc = disc.sqrt();
        let e1 = (trace + sqrt_disc) / 2.0;
        let e2 = (trace - sqrt_disc) / 2.0;
        if e1 >= e2 {
            (e1, e2)
        } else {
            (e2, e1)
        }
    }
}

/// Compute eigenvalues of a 3x3 matrix using characteristic polynomial.
///
/// Returns eigenvalues sorted in descending order (largest first)
pub fn eigenvalues_3x3(m: DMat3) -> (f64, f64, f64) {
    // Characteristic polynomial: det(A - lambda*I) = 0
    // For 3x3: -lambda^3 + tr(A)*lambda^2 - S*lambda + det(A) = 0
    // where S = sum of principal minors

    let trace = m.x_axis.x + m.y_axis.y + m.z_axis.z; // trace = sum of diagonal
    let det = m.determinant();

    // Sum of principal 2x2 minors:
    // M11 = (a22*a33 - a23*a32), M22 = (a11*a33 - a13*a31), M33 = (a11*a22 - a12*a21)
    let s = (m.y_axis.y * m.z_axis.z - m.y_axis.z * m.z_axis.y)  // M11
          + (m.x_axis.x * m.z_axis.z - m.x_axis.z * m.z_axis.x)  // M22
          + (m.x_axis.x * m.y_axis.y - m.x_axis.y * m.y_axis.x); // M33

    // Coefficients of characteristic polynomial: lambda^3 - trace*lambda^2 + s*lambda - det = 0
    let roots = solve_cubic(1.0, -trace, s, -det);

    match roots.len() {
        0 => (0.0, 0.0, 0.0),
        1 => (roots[0], roots[0], roots[0]),
        2 => {
            if roots[0] >= roots[1] {
                (roots[0], roots[1], roots[1])
            } else {
                (roots[1], roots[0], roots[0])
            }
        }
        3 => {
            let mut r = roots;
            r.sort_by(|a, b| b.partial_cmp(a).unwrap()); // Descending
            (r[0], r[1], r[2])
        }
        _ => (0.0, 0.0, 0.0),
    }
}

/// Compute the inverse of a 3x3 matrix.
///
/// Returns None if the matrix is singular
pub fn inverse_3x3(m: DMat3) -> Option<DMat3> {
    let det = determinant_3x3(m);
    if det.abs() < TOLERANCE_FLOAT_DEDUP {
        return None;
    }

    // Cofactor matrix (transpose of adjugate)
    let cofactor00 = m.y_axis.y * m.z_axis.z - m.y_axis.z * m.z_axis.y;
    let cofactor01 = -(m.y_axis.x * m.z_axis.z - m.y_axis.z * m.z_axis.x);
    let cofactor02 = m.y_axis.x * m.z_axis.y - m.y_axis.y * m.z_axis.x;

    let cofactor10 = -(m.x_axis.y * m.z_axis.z - m.x_axis.z * m.z_axis.y);
    let cofactor11 = m.x_axis.x * m.z_axis.z - m.x_axis.z * m.z_axis.x;
    let cofactor12 = -(m.x_axis.x * m.z_axis.y - m.x_axis.y * m.z_axis.x);

    let cofactor20 = m.x_axis.y * m.y_axis.z - m.x_axis.z * m.y_axis.y;
    let cofactor21 = -(m.x_axis.x * m.y_axis.z - m.x_axis.z * m.y_axis.x);
    let cofactor22 = m.x_axis.x * m.y_axis.y - m.x_axis.y * m.y_axis.x;

    let inv_det = 1.0 / det;

    Some(DMat3::from_cols(
        DVec3::new(cofactor00 * inv_det, cofactor10 * inv_det, cofactor20 * inv_det),
        DVec3::new(cofactor01 * inv_det, cofactor11 * inv_det, cofactor21 * inv_det),
        DVec3::new(cofactor02 * inv_det, cofactor12 * inv_det, cofactor22 * inv_det),
    ))
}

/// Compute the determinant of a 3x3 matrix.
pub fn determinant_3x3(m: DMat3) -> f64 {
    m.determinant()
}

// =============================================================================
// Numerical Integration
// =============================================================================

/// Simpson's rule for numerical integration.
///
/// Integrates f from a to b using n subintervals (n must be even)
pub fn simpson_integrate(f: fn(f64) -> f64, a: f64, b: f64, n: usize) -> f64 {
    let n = if !n.is_multiple_of(2) { n + 1 } else { n }; // Ensure even
    let h = (b - a) / n as f64;

    let mut sum = f(a) + f(b);

    for i in 1..n {
        let x = a + i as f64 * h;
        let coef = if i % 2 == 0 { 2.0 } else { 4.0 };
        sum += coef * f(x);
    }

    sum * h / 3.0
}

/// Gaussian quadrature nodes and weights for various orders
fn gaussian_nodes_weights(n: usize) -> Vec<(f64, f64)> {
    match n {
        1 => vec![(0.0, 2.0)],
        2 => vec![
            (-FRAC_1_SQRT_2, 1.0),
            (FRAC_1_SQRT_2, 1.0),
        ],
        3 => vec![
            (0.0, 8.0 / 9.0),
            (-0.7745966692414834, 5.0 / 9.0),
            (0.7745966692414834, 5.0 / 9.0),
        ],
        4 => vec![
            (-0.8611363115940526, 0.3478548451374538),
            (-0.3399810435848563, 0.6521451548625461),
            (0.3399810435848563, 0.6521451548625461),
            (0.8611363115940526, 0.3478548451374538),
        ],
        5 => vec![
            (0.0, 0.5688888888888889),
            (-0.5384693101056831, 0.4786286704993665),
            (0.5384693101056831, 0.4786286704993665),
            (-0.906_179_845_938_664, 0.2369268850561891),
            (0.906_179_845_938_664, 0.2369268850561891),
        ],
        6 => vec![
            (-0.932_469_514_203_152, 0.1713244923791704),
            (-0.6612093864662645, 0.3607615730481386),
            (-0.2386191860831969, 0.467_913_934_572_691),
            (0.2386191860831969, 0.467_913_934_572_691),
            (0.6612093864662645, 0.3607615730481386),
            (0.932_469_514_203_152, 0.1713244923791704),
        ],
        _ => gaussian_nodes_weights(6), // Default to 6-point rule
    }
}

/// Gaussian quadrature for numerical integration.
///
/// Integrates f from a to b using n-point Gaussian quadrature
pub fn gaussian_quadrature(f: fn(f64) -> f64, a: f64, b: f64, n_points: usize) -> f64 {
    let nodes_weights = gaussian_nodes_weights(n_points);

    // Transform from [-1, 1] to [a, b]
    let scale = (b - a) / 2.0;
    let shift = (a + b) / 2.0;

    let mut sum = 0.0;
    for (node, weight) in nodes_weights {
        let x = shift + scale * node;
        sum += weight * f(x);
    }

    sum * scale
}

// =============================================================================
// Optimization
// =============================================================================

/// Golden ratio constant
const PHI: f64 = 1.618033988749895;
const RESPHI: f64 = 0.3819660112501051; // 1/phi^2

/// Golden section search for finding minimum.
///
/// # Arguments
/// * `f` - Function to minimize
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The x value that minimizes f in [a, b]
pub fn golden_section_min<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> f64 {
    let mut lo = a;
    let mut hi = b;

    let mut c = lo + RESPHI * (hi - lo);
    let mut d = hi - RESPHI * (hi - lo);

    let mut fc = f(c);
    let mut fd = f(d);

    while (hi - lo).abs() > tol {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = lo + RESPHI * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = hi - RESPHI * (hi - lo);
            fd = f(d);
        }
    }

    (lo + hi) / 2.0
}

/// Golden section search for finding maximum.
///
/// # Arguments
/// * `f` - Function to maximize
/// * `a` - Lower bound
/// * `b` - Upper bound
/// * `tol` - Tolerance for convergence
///
/// # Returns
/// The x value that maximizes f in [a, b]
pub fn golden_section_max<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> f64 {
    golden_section_min(|x| -f(x), a, b, tol)
}

// =============================================================================
// SVD — Singular Value Decomposition (3x3)
// =============================================================================

/// Solve A*x = b for 3x3 matrix A.
///
/// Uses Gaussian elimination with partial pivoting (via inverse). If A is
/// singular, falls back to a damped least-squares approach.
pub fn svd_solve_3x3(a: DMat3, b: DVec3) -> Option<DVec3> {
    // For non-singular matrices, use direct inverse (most efficient).
    if let Some(inv) = inverse_3x3(a) {
        return Some(inv * b);
    }
    // For singular matrices, use SVD-based pseudo-inverse.
    let (u, s, vt) = svd_jacobi_3x3(a)?;
    let s_max = s.x.max(s.y.max(s.z));
    let tol = s_max * 1e-14;
    let utb = u.transpose() * b;
    let y = DVec3::new(
        if s.x > tol { utb.x / s.x } else { 0.0 },
        if s.y > tol { utb.y / s.y } else { 0.0 },
        if s.z > tol { utb.z / s.z } else { 0.0 },
    );
    Some(vt.transpose() * y)
}

/// Jacobi SVD for 3x3 matrices using one-sided Jacobi.
///
/// Returns (U, singular_values, V^T) where A = U * diag(s) * V^T.
fn svd_jacobi_3x3(mut a: DMat3) -> Option<(DMat3, DVec3, DMat3)> {
    let tol = 1e-14;
    let mut v = DMat3::IDENTITY;

    for _ in 0..100 {
        let mut max_gamma = 0.0;
        let mut best_p = 0;
        let mut best_q = 1;

        // Find pair with largest off-diagonal correlation
        for p in 0..3 {
            for q in (p + 1)..3 {
                let gamma = a.col(p).dot(a.col(q));
                let alpha = a.col(p).length_squared();
                let beta = a.col(q).length_squared();
                let denom = (alpha * beta).sqrt();
                if denom > tol {
                    let ratio = gamma.abs() / denom;
                    if ratio > max_gamma {
                        max_gamma = ratio;
                        best_p = p;
                        best_q = q;
                    }
                }
            }
        }

        if max_gamma < tol {
            break;
        }

        let p = best_p;
        let q = best_q;
        let alpha = a.col(p).length_squared();
        let beta = a.col(q).length_squared();
        let gamma = a.col(p).dot(a.col(q));

        // Compute Jacobi rotation
        let tau = (beta - alpha) / (2.0 * gamma);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        // Update columns p, q of A: [a_p' a_q'] = [a_p a_q] * J
        let ap = a.col(p);
        let aq = a.col(q);
        let new_ap = ap * c + aq * s;
        let new_aq = -ap * s + aq * c;
        a = DMat3::from_cols(
            if p == 0 { new_ap } else if q == 0 { new_aq } else { a.col(0) },
            if p == 1 { new_ap } else if q == 1 { new_aq } else { a.col(1) },
            if p == 2 { new_ap } else if q == 2 { new_aq } else { a.col(2) },
        );

        // Update V similarly
        let vp = v.col(p);
        let vq = v.col(q);
        let new_vp = vp * c + vq * s;
        let new_vq = -vp * s + vq * c;
        v = DMat3::from_cols(
            if p == 0 { new_vp } else if q == 0 { new_vq } else { v.col(0) },
            if p == 1 { new_vp } else if q == 1 { new_vq } else { v.col(1) },
            if p == 2 { new_vp } else if q == 2 { new_vq } else { v.col(2) },
        );
    }

    // Singular values = column norms of A
    let mut s = DVec3::new(a.col(0).length(), a.col(1).length(), a.col(2).length());
    let tol_s = 1e-14;
    if s.x < tol_s { s.x = 0.0; }
    if s.y < tol_s { s.y = 0.0; }
    if s.z < tol_s { s.z = 0.0; }

    // Normalize columns of A to get U
    let u0 = if s.x > tol_s { a.col(0) / s.x } else { DVec3::X };
    let u1 = if s.y > tol_s { a.col(1) / s.y } else { DVec3::Y };
    let u2 = if s.z > tol_s { a.col(2) / s.z } else { DVec3::Z };
    let u = DMat3::from_cols(u0, u1, u2);

    Some((u, s, v.transpose()))
}

// =============================================================================
// Least Squares
// =============================================================================

/// Fit a line y = a + b*x to data points using linear least squares.
///
/// Returns (intercept, slope) or None if insufficient data.
pub fn least_squares_linear(x: &[f64], y: &[f64]) -> Option<(f64, f64)> {
    let n = x.len().min(y.len());
    if n < 2 { return None; }

    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;

    for i in 0..n {
        sx += x[i];
        sy += y[i];
        sxx += x[i] * x[i];
        sxy += x[i] * y[i];
    }

    let denom = (n as f64) * sxx - sx * sx;
    if denom.abs() < TOLERANCE_FLOAT_DEDUP { return None; }

    let b = ((n as f64) * sxy - sx * sy) / denom;
    let a = (sy - b * sx) / (n as f64);
    Some((a, b))
}

// =============================================================================
// BFGS — Broyden-Fletcher-Goldfarb-Shanno optimization
// =============================================================================

/// Minimize a function with gradient using the BFGS quasi-Newton method.
///
/// * `x0` - Initial guess
/// * `f_grad` - Function returning (value) and filling gradient
/// * `tol` - Gradient norm convergence tolerance
/// * `max_iter` - Maximum iterations
///
/// Returns the minimizer location or None on failure.
pub fn bfgs_minimize(
    x0: &[f64],
    f_grad: impl Fn(&[f64], &mut [f64]) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut grad = vec![0.0; n];
    let mut f = f_grad(&x, &mut grad);

    // Initial inverse Hessian approximation: identity
    let mut h_inv = vec![0.0; n * n];
    for i in 0..n {
        h_inv[i * n + i] = 1.0;
    }

    for _ in 0..max_iter {
        // Check convergence: ||grad|| < tol
        let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gn < tol {
            return Some(x);
        }

        // Search direction: p = -H * grad
        let mut p = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                p[i] -= h_inv[i * n + j] * grad[j];
            }
        }

        // Line search: find step size alpha
        let alpha = line_search_backtracking(&x, &p, &grad, f, &f_grad);

        // Update
        let mut s = vec![0.0; n];
        for i in 0..n {
            s[i] = alpha * p[i];
            x[i] += s[i];
        }

        let mut new_grad = vec![0.0; n];
        let new_f = f_grad(&x, &mut new_grad);

        // Gradient difference: y = g_{k+1} - g_k
        let mut y = vec![0.0; n];
        for i in 0..n {
            y[i] = new_grad[i] - grad[i];
        }

        grad = new_grad;
        f = new_f;

        // BFGS update of inverse Hessian approximation
        // H_{k+1} = (I - ρ*s*y^T) * H_k * (I - ρ*y*s^T) + ρ*s*s^T
        // where ρ = 1/(y^T * s)
        let sy = s.iter().zip(y.iter()).map(|(s, y)| s * y).sum::<f64>();
        if sy.abs() < TOLERANCE_FLOAT_DEDUP { continue; }

        let rho = 1.0 / sy;

        // Compute H*y (store as hy)
        let mut hy = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                hy[i] += h_inv[i * n + j] * y[j];
            }
        }

        // Compute y^T * H * y
        let ythy = y.iter().zip(hy.iter()).map(|(y, hy)| y * hy).sum::<f64>();

        // Update: H = H + ρ * ( (1 + ρ*y^T*H*y)*s*s^T - s*y^T*H - H*y*s^T )
        let factor = 1.0 + rho * ythy;
        let mut h_new = h_inv.clone();
        for i in 0..n {
            for j in 0..n {
                h_new[i * n + j] += rho * (factor * s[i] * s[j] - s[i] * hy[j] - hy[i] * s[j]);
            }
        }
        h_inv = h_new;
    }

    // Check final gradient
    let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
    if gn < tol { Some(x) } else { None }
}

/// Simple backtracking line search.
fn line_search_backtracking(
    x: &[f64],
    p: &[f64],
    _grad: &[f64],
    f_current: f64,
    f_grad: impl Fn(&[f64], &mut [f64]) -> f64,
) -> f64 {
    let c = 1e-4;
    let rho = 0.5;
    let mut alpha = 1.0;
    let mut trial = vec![0.0; x.len()];
    let mut g = vec![0.0; x.len()];

    for _ in 0..20 {
        for i in 0..x.len() {
            trial[i] = x[i] + alpha * p[i];
        }
        let f_trial = f_grad(&trial, &mut g);
        if f_trial <= f_current + c * alpha * p.iter().zip(_grad.iter()).map(|(p, g)| p * g).sum::<f64>() {
            return alpha;
        }
        alpha *= rho;
    }
    alpha
}

// =============================================================================
// Newton Minimization (with Hessian)
// =============================================================================

/// Minimize a function using Newton's method with gradient and Hessian.
///
/// * `x0` - Initial guess
/// * `f_grad_hess` - Function returning value, filling gradient and Hessian
///   (Hessian stored as flat array: row-major, n×n)
/// * `tol` - Gradient norm convergence tolerance
/// * `max_iter` - Maximum iterations
///
/// Returns the minimizer location or None on failure.
pub fn newton_minimize(
    x0: &[f64],
    f_grad_hess: impl Fn(&[f64], &mut [f64], &mut [f64]) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut grad = vec![0.0; n];
    let mut hess = vec![0.0; n * n];
    let mut f = f_grad_hess(&x, &mut grad, &mut hess);

    for _ in 0..max_iter {
        let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gn < tol {
            return Some(x);
        }

        // Solve H * p = -g for p using simple Gaussian elimination
        if let Some(p) = solve_linear_system(&hess, &grad.iter().map(|g| -g).collect::<Vec<_>>(), n) {
            // Line search: use simple bisection with gradient evaluation
            let mut alpha = 1.0;
            let mut trial = x.clone();
            let mut g_trial = vec![0.0; n];
            for _ in 0..20 {
                for i in 0..n {
                    trial[i] = x[i] + alpha * p[i];
                }
                let f_trial = f_grad_hess(&trial, &mut g_trial, &mut hess);
                let directional = p.iter().zip(grad.iter()).map(|(p, g)| p * g).sum::<f64>();
                if f_trial <= f + 1e-4 * alpha * directional {
                    break;
                }
                alpha *= 0.5;
            }

            for i in 0..n {
                x[i] += alpha * p[i];
            }
            f = f_grad_hess(&x, &mut grad, &mut hess);
        } else {
            return None;
        }
    }

    let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
    if gn < tol { Some(x) } else { None }
}

/// Solve A*x = b for small n×n system using Gaussian elimination.
fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    // Augmented matrix
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }

    // Forward elimination
    for col in 0..n {
        // Partial pivot
        let mut max_row = col;
        for row in (col + 1)..n {
            if aug[row * (n + 1) + col].abs() > aug[max_row * (n + 1) + col].abs() {
                max_row = row;
            }
        }
        if aug[max_row * (n + 1) + col].abs() < TOLERANCE_FLOAT_DEDUP {
            return None;
        }
        if max_row != col {
            for j in col..=n {
                aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
        }

        // Eliminate below
        for row in (col + 1)..n {
            let factor = aug[row * (n + 1) + col] / aug[col * (n + 1) + col];
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i * (n + 1) + n];
        for j in (i + 1)..n {
            sum -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (n + 1) + i];
    }
    Some(x)
}

// =============================================================================
// FRPR — Fletcher-Reeves Polak-Ribiere conjugate gradient optimization
// =============================================================================

/// Minimize function with gradient using FRPR conjugate gradient method.
pub fn frpr_minimize(
    x0: &[f64],
    f_grad: impl Fn(&[f64], &mut [f64]) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut grad = vec![0.0; n];
    let _f = f_grad(&x, &mut grad);
    let mut d = grad.iter().map(|g| -g).collect::<Vec<_>>();
    let mut prev_norm2 = grad.iter().map(|g| g * g).sum::<f64>();

    for _ in 0..max_iter {
        if prev_norm2.sqrt() < tol { return Some(x); }
        let alpha = line_search_backtracking(&x, &d, &grad, _f, &f_grad);
        for i in 0..n { x[i] += alpha * d[i]; }
        let mut new_grad = vec![0.0; n];
        let _f = f_grad(&x, &mut new_grad);
        let new_norm2 = new_grad.iter().map(|g| g * g).sum::<f64>();
        let gg_diff: f64 = new_grad.iter().zip(grad.iter()).map(|(ng, g)| ng * (ng - g)).sum();
        let beta = if prev_norm2 > TOLERANCE_FLOAT_DEDUP { (gg_diff / prev_norm2).max(0.0) } else { 0.0 };
        for i in 0..n { d[i] = -new_grad[i] + beta * d[i]; }
        grad = new_grad;
        prev_norm2 = new_norm2;
    }
    if prev_norm2.sqrt() < tol { Some(x) } else { None }
}

// =============================================================================
// Powell — Derivative-free direction-set optimization
// =============================================================================

/// Minimize using Powell's derivative-free conjugate direction method.
pub fn powell_minimize(x0: &[f64], f: impl Fn(&[f64]) -> f64, tol: f64, max_iter: usize) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut dirs: Vec<Vec<f64>> = (0..n).map(|i| { let mut d = vec![0.0; n]; d[i] = 1.0; d }).collect();
    let mut prev_f = f(&x);
    for _ in 0..max_iter {
        let x_start = x.clone();
        let mut delta = 0.0;
        let mut best_dir = 0;
        for i in 0..n {
            // Minimize along direction i via golden section
            let (xn, fn_) = minimize_1d(&x, &dirs[i], &f);
            let dec = prev_f - fn_;
            if dec > delta { delta = dec; best_dir = i; }
            x = xn; prev_f = fn_;
        }
        if x.iter().zip(x_start.iter()).map(|(a, b)| (a - b).abs()).sum::<f64>() < tol { return Some(x); }
        // New conjugate direction
        let mut nd = Vec::with_capacity(n);
        for i in 0..n { nd.push(x[i] - x_start[i]); }
        let nn = nd.iter().map(|d| d * d).sum::<f64>().sqrt();
        if nn > TOLERANCE_FLOAT_DEDUP {
            for i in 0..n { nd[i] /= nn; }
            for i in best_dir..n - 1 { dirs.swap(i, i + 1); }
            *dirs.last_mut().unwrap() = nd;
        }
    }
    None
}

fn minimize_1d(x: &[f64], dir: &[f64], f: impl Fn(&[f64]) -> f64) -> (Vec<f64>, f64) {
    let n = x.len();
    let at = |alpha: f64| -> f64 { let mut t = x.to_vec(); for i in 0..n { t[i] += alpha * dir[i]; } f(&t) };
    let alpha = golden_section_min(at, 0.0, 10.0, 1e-10);
    let mut xn = x.to_vec();
    for i in 0..n { xn[i] += alpha * dir[i]; }
    let fxn = f(&xn);
    (xn, fxn)
}

// =============================================================================
// Householder — QR via Householder reflections
// =============================================================================

/// Solve A*x = b using Householder QR decomposition.
pub fn householder_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    if n == 0 { return None; }
    let mut r = a.to_vec();
    let mut rhs = b.to_vec();
    for k in 0..n - 1 { // Skip last column (no rows below to zero out)
        let mut nrm2 = 0.0;
        for i in k..n { nrm2 += r[i * n + k] * r[i * n + k]; }
        if nrm2 < TOLERANCE_FLOAT_DEDUP { continue; }
        let nrm = nrm2.sqrt();
        let sign = if r[k * n + k] >= 0.0 { -1.0 } else { 1.0 };
        let beta = -sign * nrm;
        let vk = r[k * n + k] - sign * nrm;
        for i in (k + 1)..n { r[i * n + k] /= vk; }
        r[k * n + k] = beta;
        for j in (k + 1)..n {
            let mut sp = r[k * n + k] * r[k * n + j];
            for i in (k + 1)..n { sp += r[i * n + k] * r[i * n + j]; }
            let tau = sp / beta;
            r[k * n + j] -= tau * r[k * n + k];
            for i in (k + 1)..n { r[i * n + j] -= tau * r[i * n + k]; }
        }
        let mut sp = r[k * n + k] * rhs[k];
        for i in (k + 1)..n { sp += r[i * n + k] * rhs[i]; }
        let tau = sp / beta;
        rhs[k] -= tau * r[k * n + k];
        for i in (k + 1)..n { rhs[i] -= tau * r[i * n + k]; }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..n { s -= r[i * n + j] * x[j]; }
        if r[i * n + i].abs() < TOLERANCE_FLOAT_DEDUP { return None; }
        x[i] = s / r[i * n + i];
    }
    Some(x)
}

// =============================================================================
// Crout — LU decomposition with partial pivoting
// =============================================================================

/// Solve A*x = b using Crout LU decomposition.
pub fn crout_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    if n == 0 { return None; }
    let mut lu = a.to_vec();
    let mut p: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut max_val = lu[p[k] * n + k].abs();
        let mut max_r = k;
        for i in (k + 1)..n { let v = lu[p[i] * n + k].abs(); if v > max_val { max_val = v; max_r = i; } }
        p.swap(k, max_r);
        if lu[p[k] * n + k].abs() < TOLERANCE_FLOAT_DEDUP { return None; }
        for i in (k + 1)..n {
            let f = lu[p[i] * n + k] / lu[p[k] * n + k];
            lu[p[i] * n + k] = f;
            for j in (k + 1)..n { lu[p[i] * n + j] -= f * lu[p[k] * n + j]; }
        }
    }
    let mut y = vec![0.0; n];
    for i in 0..n { let mut s = b[p[i]]; for j in 0..i { s -= lu[p[i] * n + j] * y[j]; } y[i] = s; }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for j in (i + 1)..n { s -= lu[p[i] * n + j] * x[j]; }
        if lu[p[i] * n + i].abs() < TOLERANCE_FLOAT_DEDUP { return None; }
        x[i] = s / lu[p[i] * n + i];
    }
    Some(x)
}

// =============================================================================
// BissecNewton — Hybrid bisection + Newton root finding
// =============================================================================

/// Find root using hybrid bisection-Newton method (safe + fast).
pub fn biss_newton(f: impl Fn(f64) -> f64, df: impl Fn(f64) -> f64,
    a: f64, b: f64, tol: f64) -> Option<f64> {
    let mut lo = a.min(b); let mut hi = a.max(b);
    let mut x = (lo + hi) / 2.0;
    for _ in 0..100 {
        let fx = f(x);
        if fx.abs() < tol { return Some(x); }
        let dfx = df(x);
        let xn = if dfx.abs() > TOLERANCE_FLOAT_DEDUP {
            let nv = x - fx / dfx;
            if nv > lo && nv < hi { nv } else { (lo + hi) / 2.0 }
        } else { (lo + hi) / 2.0 };
        if (xn - x).abs() < tol { return Some(xn); }
        if f(lo) * f(xn) <= 0.0 { hi = xn; } else { lo = xn; }
        x = xn;
    }
    if f(x).abs() < tol { Some(x) } else { None }
}

// =============================================================================
// TrigonometricFunctionRoots
// Solve: a*cos²(x) + 2b*cos(x)sin(x) + c*cos(x) + d*sin(x) + e = 0
// =============================================================================

fn trig_f(a: f64, b: f64, c: f64, d: f64, e: f64) -> impl Fn(f64) -> f64 {
    move |x: f64| { let (s, c_) = x.sin_cos(); a * c_ * c_ + 2.0 * b * c_ * s + c * c_ + d * s + e }
}

/// Find roots of trigonometric equation in [x_min, x_max].
pub fn trig_roots(a: f64, b: f64, c: f64, d: f64, e: f64, x_min: f64, x_max: f64) -> Vec<f64> {
    let f = trig_f(a, b, c, d, e);
    find_roots_in(f, x_min, x_max, 200)
}

/// Solve d*sin(x) + e = 0 in [x_min, x_max].
pub fn trig_roots_sin_only(d: f64, e: f64, x_min: f64, x_max: f64) -> Vec<f64> {
    trig_roots(0.0, 0.0, 0.0, d, e, x_min, x_max)
}

/// Solve c*cos(x) + d*sin(x) + e = 0 in [x_min, x_max].
pub fn trig_roots_cos_sin(c: f64, d: f64, e: f64, x_min: f64, x_max: f64) -> Vec<f64> {
    trig_roots(0.0, 0.0, c, d, e, x_min, x_max)
}

// =============================================================================
// Laguerre — Polynomial root finding via Laguerre's method
// =============================================================================

/// Find all roots of polynomial a₀ + a₁x + ... + a_nx^n using Laguerre's method.
pub fn laguerre_roots(coeffs: &[f64]) -> Vec<f64> {
    let mut a: Vec<f64> = coeffs.iter().copied().collect();
    while a.len() > 1 && a.last().map_or(false, |&c| c.abs() < TOLERANCE_FLOAT_DEDUP) { a.pop(); }
    let n = a.len() - 1;
    if n == 0 { return vec![]; }
    if n <= 4 { return match n { 1 => vec![-a[0] / a[1]], 2 => solve_quadratic(a[2], a[1], a[0]), 3 => solve_cubic(a[3], a[2], a[1], a[0]), _ => solve_quartic(a[4], a[3], a[2], a[1], a[0]) }; }

    let mut roots = Vec::new();
    let mut deg = n;
    while deg >= 2 {
        let mut x = 0.0;
        for _ in 0..200 {
            let mut p = vec![0.0; deg + 1]; p[deg] = a[deg];
            for i in (0..deg).rev() { p[i] = a[i] + x * p[i + 1]; }
            let fv = p[0];
            if fv.abs() < 1e-14 { break; }
            let mut p1 = vec![0.0; deg]; for i in 0..deg { p1[i] = a[i + 1] * (i as f64 + 1.0); }
            let mut fp = 0.0; for i in (0..deg - 1).rev() { fp = fp * x + p1[i]; }
            let mut fp2 = 0.0; for i in (0..deg - 2).rev() { fp2 = fp2 * x + p1[i + 1] * (i as f64 + 1.0); }
            let g = fp / fv; let h = g * g - fp2 / fv;
            let d_ = ((deg as f64 - 1.0) * (deg as f64 * h - g * g)).sqrt();
            let s = if (g + d_).abs() > (g - d_).abs() { deg as f64 / (g + d_) } else { deg as f64 / (g - d_) };
            if s.abs() < 1e-14 { break; }
            x -= s;
            if s.abs() < 1e-12 { break; }
        }
        roots.push(x);
        let mut b = vec![0.0; deg]; b[deg - 1] = a[deg];
        for i in (0..deg - 1).rev() { b[i] = a[i + 1] + x * b[i + 1]; }
        a = b; deg = a.len() - 1;
    }
    if deg == 1 { roots.push(-a[0] / a[1]); }
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

// =============================================================================
// BrentMinimum — Brent's method for 1D minimization
// =============================================================================

/// Minimize a 1D function using Brent's method.
pub fn brent_minimize(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> f64 {
    const PHI: f64 = 0.3819660112501051;
    let (mut lo, mut hi) = if a < b { (a, b) } else { (b, a) };
    let (mut x, mut w, mut v) = ((a + b) / 2.0, (a + b) / 2.0, (a + b) / 2.0);
    let (mut fx, mut fw, mut fv) = (f(x), f(x), f(x));
    let (mut d, mut e): (f64, f64) = (0.0, 0.0);
    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        let tol1 = tol * x.abs() + 1e-12;
        let tol2 = 2.0 * tol1;
        if (x - mid).abs() <= tol2 - (hi - lo) / 2.0 { return x; }
        let mut use_para = false;
        let mut u = 0.0;
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv); let qq = (x - v) * (fx - fw);
            let p = (x - v) * qq - (x - w) * r; let q = 2.0 * (qq - r);
            if q.abs() > tol1 { u = x - p / q; if u > lo + tol1 && u < hi - tol1 && (u - x).abs() < e { use_para = true; } }
        }
        if !use_para { u = if x >= mid { x - PHI * (x - lo) } else { x + PHI * (hi - x) }; e = d; d = u - x; }
        let fu = f(u);
        if fu <= fx {
            if u >= x { lo = x; } else { hi = x; }
            v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
        } else {
            if u >= x { hi = u; } else { lo = u; }
            if fu <= fw || w == x { v = w; fv = fw; w = u; fw = fu; }
            else if fu <= fv || v == x || v == w { v = u; fv = fu; }
        }
    }
    x
}

// =============================================================================
// MultipleRoots — Find all roots of f(x) = 0 in [a, b]
// =============================================================================

/// Find all roots in [a, b] by scanning for sign changes.
pub fn find_roots_in(f: impl Fn(f64) -> f64, a: f64, b: f64, n_intervals: usize) -> Vec<f64> {
    let step = (b - a) / n_intervals.max(1) as f64;
    let mut roots = Vec::new();
    for i in 0..n_intervals {
        let x1 = a + i as f64 * step; let x2 = x1 + step;
        let f1 = f(x1); let f2 = f(x2);
        if f1 * f2 < 0.0 { if let Some(r) = bisection(&f, x1, x2, 1e-10) { roots.push(r); } }
        else if f1.abs() < 1e-10 && !roots.iter().any(|r| (r - x1).abs() < 1e-8) { roots.push(x1); }
    }
    if f(b).abs() < 1e-10 && !roots.iter().any(|r| (r - b).abs() < 1e-8) { roots.push(b); }
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

/// Find a bracket [a,b] where f(a) and f(b) have opposite signs.
pub fn bracket_root(f: impl Fn(f64) -> f64, x0: f64, step: f64, max_steps: usize) -> Option<(f64, f64)> {
    let mut a = x0; let mut fa = f(a);
    for _ in 0..max_steps {
        let b = a + step; let fb = f(b);
        if fa * fb <= 0.0 { return Some((a, b)); }
        a = b; fa = fb;
    }
    None
}

// =============================================================================
// MathFunctor — General function evaluation utility
// =============================================================================

/// Evaluate a polynomial using Horner's method.
pub fn poly_eval(coeffs: &[f64], x: f64) -> f64 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

// =============================================================================
// GlobOptMin — Global optimization via grid + local refinement
// =============================================================================

/// Global optimizer: coarse grid evaluation + local coordinate descent.
pub fn glob_opt_min(
    f: impl Fn(&[f64]) -> f64,
    lower: &[f64],
    upper: &[f64],
    n_cells: usize,
    n_local: usize,
) -> Vec<f64> {
    let n = lower.len();
    if n == 0 { return vec![]; }
    let mut candidates: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut current = vec![0.0; n];
    fn grid_eval(f: &impl Fn(&[f64]) -> f64, lower: &[f64], upper: &[f64], nc: usize,
        cand: &mut Vec<(f64, Vec<f64>)>, cur: &mut Vec<f64>, dim: usize) {
        if dim == cur.len() { cand.push((f(cur), cur.clone())); return; }
        let step = (upper[dim] - lower[dim]) / nc as f64;
        for i in 0..=nc { cur[dim] = lower[dim] + i as f64 * step; grid_eval(f, lower, upper, nc, cand, cur, dim + 1); }
    }
    grid_eval(&f, lower, upper, n_cells, &mut candidates, &mut current, 0);
    candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    candidates.truncate(n_local.max(1));
    let mut best_x = candidates[0].1.clone();
    let mut best_f = candidates[0].0;
    for (_, mut x) in candidates {
        for _ in 0..5 {
            for i in 0..n {
                let step = (upper[i] - lower[i]) * 0.01;
                let fx = f(&x);
                let mut xt = x.clone(); xt[i] += step;
                if xt[i] <= upper[i] && f(&xt) < fx { x[i] = xt[i]; continue; }
                xt[i] = x[i] - step;
                if xt[i] >= lower[i] && f(&xt) < fx { x[i] = xt[i]; }
            }
        }
        let fv = f(&x);
        if fv < best_f { best_f = fv; best_x = x; }
    }
    best_x
}

// =============================================================================
// PSO — Particle Swarm Optimization
// =============================================================================

/// Minimize a function using Particle Swarm Optimization.
pub fn pso_minimize(
    f: impl Fn(&[f64]) -> f64,
    lower: &[f64],
    upper: &[f64],
    n_particles: usize,
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = lower.len();
    let mut rng = fastrand::Rng::new();
    let mut positions: Vec<Vec<f64>> = Vec::with_capacity(n_particles);
    let mut velocities: Vec<Vec<f64>> = Vec::with_capacity(n_particles);
    let mut pbest: Vec<Vec<f64>> = Vec::with_capacity(n_particles);
    let mut pbest_val: Vec<f64> = Vec::with_capacity(n_particles);
    for _ in 0..n_particles {
        let mut pos = Vec::with_capacity(n);
        for i in 0..n { pos.push(lower[i] + rng.f64() * (upper[i] - lower[i])); }
        let fv = f(&pos);
        positions.push(pos.clone()); velocities.push(vec![0.0; n]);
        pbest.push(pos); pbest_val.push(fv);
    }
    let mut gbest = pbest[0].clone();
    let mut gbest_val = pbest_val[0];
    for i in 1..n_particles { if pbest_val[i] < gbest_val { gbest = pbest[i].clone(); gbest_val = pbest_val[i]; } }
    for _ in 0..max_iter {
        let prev = gbest_val;
        for i in 0..n_particles {
            for j in 0..n {
                velocities[i][j] = 0.72 * velocities[i][j] + 1.49 * rng.f64() * (pbest[i][j] - positions[i][j])
                    + 1.49 * rng.f64() * (gbest[j] - positions[i][j]);
                let vmax = (upper[j] - lower[j]) * 0.2;
                velocities[i][j] = velocities[i][j].clamp(-vmax, vmax);
                positions[i][j] = (positions[i][j] + velocities[i][j]).clamp(lower[j], upper[j]);
            }
            let fv = f(&positions[i]);
            if fv < pbest_val[i] { pbest_val[i] = fv; pbest[i] = positions[i].clone(); }
            if fv < gbest_val { gbest_val = fv; gbest = positions[i].clone(); }
        }
        if (prev - gbest_val).abs() < tol && gbest_val.abs() > 1e-15 { break; }
    }
    gbest
}

// =============================================================================
// LM — Levenberg-Marquardt nonlinear least squares
// =============================================================================

/// Levenberg-Marquardt solver for nonlinear least squares.
///
/// Minimizes 0.5 * Σ f_i(x)² where f_i are residuals.
/// `func` fills the residual vector `f` (n_eq) and Jacobian `J` (n_eq × n, column-major),
/// and returns the sum-of-squares value 0.5*Σf_i².
pub fn lm_solve(
    x0: &[f64],
    mut func: impl FnMut(&[f64], &mut [f64], &mut [f64]) -> f64,
    n_eq: usize,
    max_iter: usize,
    tol: f64,
) -> Option<Vec<f64>> {
    let n = x0.len();
    if n == 0 || n_eq == 0 { return None; }
    let mut x = x0.to_vec();
    let mut f = vec![0.0; n_eq];
    let mut jac = vec![0.0; n_eq * n];
    let mut lambda = 1.0;
    let mut cost = func(&x, &mut f, &mut jac);
    for _ in 0..max_iter {
        let mut jtj = vec![0.0; n * n];
        let mut jtf = vec![0.0; n];
        for i in 0..n {
            for k in 0..n {
                let mut s = 0.0;
                for r in 0..n_eq { s += jac[r * n + i] * jac[r * n + k]; }
                jtj[i * n + k] = s;
            }
        }
        for i in 0..n {
            let mut s = 0.0;
            for r in 0..n_eq { s += jac[r * n + i] * f[r]; }
            jtf[i] = s;
        }
        if jtf.iter().map(|v| v * v).sum::<f64>().sqrt() < tol { return Some(x); }
        let mut h = jtj.clone();
        for i in 0..n { h[i * n + i] += lambda; }
        let rhs: Vec<f64> = jtf.iter().map(|v| -v).collect();
        let delta = solve_linear_system(&h, &rhs, n);
        let (cost_new, gain_ratio) = match delta {
            Some(ref d) => {
                let mut xn = vec![0.0; n];
                for i in 0..n { xn[i] = x[i] + d[i]; }
                let mut fn_ = vec![0.0; n_eq];
                let mut jn = vec![0.0; n_eq * n];
                let cn = func(&xn, &mut fn_, &mut jn);
                let pred: f64 = 0.5 * d.iter().zip(rhs.iter()).map(|(d, r)| d * r).sum::<f64>();
                let gr = if pred.abs() > 1e-15 { (cost - cn) / pred } else { 0.0 };
                (cn, gr)
            }
            None => (cost, -1.0)
        };
        if gain_ratio > 0.0 {
            let d = delta.as_ref().unwrap();
            for i in 0..n { x[i] += d[i]; }
            let _ = func(&x, &mut f, &mut jac);
            cost = cost_new;
            lambda *= if gain_ratio > 0.75 { 0.5 } else if gain_ratio < 0.25 { 2.0 } else { 1.0 };
            lambda = lambda.max(1e-12).min(1e12);
            let dn: f64 = d.iter().map(|d| d * d).sum::<f64>().sqrt();
            if dn < tol { return Some(x); }
        } else {
            lambda *= 2.0;
            if lambda > 1e12 { break; }
        }
    }
    if cost < tol { Some(x) } else { None }
}

// =============================================================================
// GeomPlate — Thin-plate spline surface from constraint points
// =============================================================================

/// Thin-plate spline radial basis function: φ(r) = r² * log(r)
fn tps_rbf(r: f64) -> f64 {
    if r < 1e-15 { 0.0 } else { r * r * r.ln() }
}

/// Solve a thin-plate spline interpolation problem.
///
/// Given N constraint points (xᵢ, yᵢ, zᵢ), finds the TPS function:
///   f(x,y) = a₀ + a₁*x + a₂*y + Σ wᵢ * φ(‖(x,y) - (xᵢ,yᵢ)‖)
///
/// Returns (w, a) where w is N coefficient weights and a is [a₀, a₁, a₂].
pub fn thin_plate_spline(points: &[DVec3]) -> Option<(Vec<f64>, Vec<f64>)> {
    let n = points.len();
    if n < 3 { return None; }

    // Build the (N+3)×(N+3) system:
    // [ K   T ] [w]   [z]
    // [ T^T 0 ] [a] = [0]
    let m = n + 3;
    let mut mat = vec![0.0; m * m];
    let mut rhs = vec![0.0; m];

    for i in 0..n {
        for j in 0..n {
            let r = (points[i].truncate() - points[j].truncate()).length();
            mat[i * m + j] = tps_rbf(r);
        }
        mat[i * m + n + 0] = 1.0;      // T: [1, x, y]
        mat[i * m + n + 1] = points[i].x;
        mat[i * m + n + 2] = points[i].y;
        mat[n * m + i] = 1.0;           // T^T
        mat[(n + 1) * m + i] = points[i].x;
        mat[(n + 2) * m + i] = points[i].y;
        rhs[i] = points[i].z;
    }

    solve_linear_system(&mat, &rhs, m).map(|sol| {
        let w = sol[0..n].to_vec();
        let a = vec![sol[n], sol[n + 1], sol[n + 2]];
        (w, a)
    })
}

/// Evaluate a thin-plate spline at position (x, y).
///
/// * `w` — N weights from `thin_plate_spline`
/// * `a` — [a₀, a₁, a₂] affine coefficients
/// * `points` — the N constraint points used in the TPS solve
pub fn evaluate_tps(x: f64, y: f64, w: &[f64], a: &[f64], points: &[DVec3]) -> f64 {
    let mut f = a[0] + a[1] * x + a[2] * y;
    for (i, &wi) in w.iter().enumerate() {
        let r = DVec2::new(x - points[i].x, y - points[i].y).length();
        f += wi * tps_rbf(r);
    }
    f
}

/// Build a plate surface (BSplineSurface) from constraint points.
///
/// Uses thin-plate spline interpolation, evaluates on a regular grid,
/// and constructs a BSpline surface. Returns `Surface3::BSpline`.
///
/// * `constraints` — constraint points (x, y, z)
/// * `n_u` — number of control points in u-direction
/// * `n_v` — number of control points in v-direction
/// * `tol` — TPS solver tolerance
pub fn build_plate_surface(constraints: &[DVec3], n_u: usize, n_v: usize) -> Option<Surface3> {
    if constraints.len() < 3 || n_u < 2 || n_v < 2 { return None; }

    // Find bounding box of constraints
    let mut x_min = f64::INFINITY; let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY; let mut y_max = f64::NEG_INFINITY;
    for p in constraints {
        x_min = x_min.min(p.x); x_max = x_max.max(p.x);
        y_min = y_min.min(p.y); y_max = y_max.max(p.y);
    }
    if (x_max - x_min).abs() < 1e-10 { x_max = x_min + 1.0; }
    if (y_max - y_min).abs() < 1e-10 { y_max = y_min + 1.0; }
    let pad_x = (x_max - x_min) * 0.1;
    let pad_y = (y_max - y_min) * 0.1;
    x_min -= pad_x; x_max += pad_x;
    y_min -= pad_y; y_max += pad_y;

    // Solve TPS
    let (w, a) = thin_plate_spline(constraints)?;

    // Evaluate on grid
    let mut cp = vec![vec![DVec3::ZERO; n_u]; n_v];
    for vi in 0..n_v {
        let y = y_min + (y_max - y_min) * vi as f64 / (n_v - 1) as f64;
        for ui in 0..n_u {
            let x = x_min + (x_max - x_min) * ui as f64 / (n_u - 1) as f64;
            let z = evaluate_tps(x, y, &w, &a, constraints);
            cp[vi][ui] = DVec3::new(x, y, z);
        }
    }

    // Build clamped cubic BSpline surface
    let degree = 3.min(n_u - 1).min(n_v - 1);
    let (knots_u, knots_v) = build_bspline_knots(n_u, n_v, degree);

    Some(Surface3::BSpline(BSplineSurface {
        degree_u: degree,
        degree_v: degree,
        knots_u,
        knots_v,
        control_points: cp,
        weights: vec![vec![1.0; n_u]; n_v],
    }))
}

/// Build clamped knot vectors for BSpline surface of given degree.
fn build_bspline_knots(n_u: usize, n_v: usize, degree: usize) -> (Vec<f64>, Vec<f64>) {
    let build = |n: usize| -> Vec<f64> {
        let nk = n + degree + 1;
        if n <= degree {
            return vec![0.0; nk];
        }
        let mut k = Vec::with_capacity(nk);
        for _ in 0..=degree { k.push(0.0); }
        let n_int = n - degree - 1;
        if n_int > 0 {
            let step = 1.0 / (n_int + 1) as f64;
            for i in 1..=n_int { k.push(i as f64 * step); }
        }
        for _ in 0..=degree { k.push(1.0); }
        k
    };
    (build(n_u), build(n_v))
}
