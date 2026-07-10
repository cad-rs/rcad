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
pub fn bisection(f: fn(f64) -> f64, a: f64, b: f64, tol: f64) -> Option<f64> {
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
// Tests
// =============================================================================


