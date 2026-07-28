// =============================================================================
// Internal helper functions
// =============================================================================

use super::*;

/// Get the domain for a 2D curve, handling special cases.
pub(crate) fn curve2d_domain(curve: &Curve2d) -> [f64; 2] {
    match curve {
        Curve2d::Line(_) => [-1e6, 1e6], // Clamp infinite lines
        Curve2d::Circle(_) => [0.0, 2.0 * PI],
        Curve2d::Ellipse(_) => [0.0, 2.0 * PI],
        Curve2d::CircleInvolute(_) => [-10.0, 10.0], // Practical range
        Curve2d::ArchimedeanSpiral(_) => [0.0, 6.0 * PI], // ~3 turns
        Curve2d::LogarithmicSpiral(_) => [0.0, 4.0 * PI], // ~2 turns
        Curve2d::SineWave(_) => [-10.0, 10.0],
        Curve2d::Parabola(_) => [-1e4, 1e4],
        Curve2d::Hyperbola(_) => [-1e4, 1e4],
        Curve2d::BSpline(bspline) => {
            let n = bspline.knots.len();
            if n < 2 {
                return [0.0, 1.0];
            }
            [
                bspline.knots[bspline.degree],
                bspline.knots[n - bspline.degree - 1],
            ]
        }
        Curve2d::Bezier(_) => [0.0, 1.0],
        Curve2d::Trimmed(tc) => curve2d_domain(tc.curve.as_ref()),
        Curve2d::Offset(c) => curve2d_domain(c.basis.as_ref()),
        Curve2d::AHTBezier(_) => [0.0, 1.0],
        Curve2d::TBezier(c) => [0.0, std::f64::consts::PI / c.alpha],
    }
}

/// Compute the first derivative of a 2D curve using finite differences.
pub(crate) fn curve2d_derivative(curve: &Curve2d, t: f64) -> DVec2 {
    curve.derivative_at(t)
}

/// Compute the second derivative of a 2D curve via the Curve2dEval trait.
pub(crate) fn curve2d_second_derivative(curve: &Curve2d, t: f64) -> DVec2 {
    curve.derivative2_at(t)
}

/// Compute the unit tangent vector of a 2D curve.
pub(crate) fn curve2d_tangent(curve: &Curve2d, t: f64) -> DVec2 {
    let d = curve2d_derivative(curve, t);
    let len = d.length();
    if len < TOLERANCE_FLOAT_DEDUP {
        DVec2::X
    } else {
        d / len
    }
}

/// Newton refinement for curve-curve intersection.
pub(crate) fn refine_curve2d_intersection(
    curve1: &Curve2d,
    curve2: &Curve2d,
    domain1: [f64; 2],
    domain2: [f64; 2],
    t1: f64,
    t2: f64,
) -> (f64, f64) {
    let mut t1 = t1;
    let mut t2 = t2;

    const MAX_ITER: usize = 30;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let p1 = curve1.point_at(t1);
        let p2 = curve2.point_at(t2);

        let d1 = curve2d_derivative(curve1, t1);
        let d2 = curve2d_derivative(curve2, t2);

        let diff = p1 - p2;

        // Gradient of distance squared
        let f1 = diff.dot(d1);
        let f2 = -diff.dot(d2);

        // Hessian (second derivatives)
        let d1_2 = curve2d_second_derivative(curve1, t1);
        let d2_2 = curve2d_second_derivative(curve2, t2);

        let h11 = d1.dot(d1) + diff.dot(d1_2);
        let h22 = d2.dot(d2) - diff.dot(d2_2);
        let h12 = -d1.dot(d2);

        let det = h11 * h22 - h12 * h12;
        if det.abs() < TOL {
            break;
        }

        let dt1 = (-f1 * h22 + f2 * h12) / det;
        let dt2 = (-f2 * h11 + f1 * h12) / det;

        t1 += dt1;
        t2 += dt2;

        t1 = t1.clamp(domain1[0], domain1[1]);
        t2 = t2.clamp(domain2[0], domain2[1]);

        if dt1.abs() < TOL && dt2.abs() < TOL {
            break;
        }
    }

    (t1, t2)
}

/// Newton refinement for point-to-curve distance.
pub(crate) fn refine_point_curve2d_distance(
    curve: &Curve2d,
    domain: [f64; 2],
    point: DVec2,
    initial_t: f64,
) -> f64 {
    let mut t = initial_t;

    const MAX_ITER: usize = 20;
    const TOL: f64 = TOLERANCE_LINEAR_ULTRA_STRICT;

    for _ in 0..MAX_ITER {
        let p = curve.point_at(t);
        let d = curve2d_derivative(curve, t);

        let diff = p - point;
        let f = diff.dot(d);

        let d2 = curve2d_second_derivative(curve, t);
        let df = d.dot(d) + diff.dot(d2);

        if df.abs() < TOL {
            break;
        }

        let delta = -f / df;
        t += delta;

        t = t.clamp(domain[0], domain[1]);

        if delta.abs() < TOL {
            break;
        }
    }

    t
}

/// Newton refinement for curve-to-curve distance.
pub(crate) fn refine_curve2d_distance(
    curve1: &Curve2d,
    curve2: &Curve2d,
    domain1: [f64; 2],
    domain2: [f64; 2],
    t1: f64,
    t2: f64,
) -> (f64, f64) {
    // Reuse intersection refinement (same mathematics)
    refine_curve2d_intersection(curve1, curve2, domain1, domain2, t1, t2)
}

// =============================================================================
// Interpolation helpers (from kernel fit.rs)
// =============================================================================

/// Chord-length parameterization for 2D points, normalized to [0, 1].
pub(crate) fn chord_length_params_2d(pts: &[DVec2]) -> Vec<f64> {
    let n = pts.len();
    let mut params = Vec::with_capacity(n);
    params.push(0.0_f64);
    let mut total = 0.0_f64;
    for i in 1..n {
        total += (pts[i] - pts[i - 1]).length();
        params.push(total);
    }
    if total < TOLERANCE_FLOAT_LOOSE {
        return vec![0.0; n];
    }
    for p in &mut params {
        *p /= total;
    }
    params
}

/// Clamped knot vector derived from parameters.
pub(crate) fn clamped_knots_from_params(params: &[f64], degree: usize) -> Vec<f64> {
    let n = params.len();
    let m = n + degree + 1;
    let mut knots = vec![0.0_f64; m];

    // First degree+1 knots = 0
    for knot in knots.iter_mut().take(degree + 1) {
        *knot = 0.0;
    }
    // Last degree+1 knots = 1
    for knot in knots.iter_mut().skip(m - degree - 1) {
        *knot = 1.0;
    }
    // Interior knots: average of degree consecutive params
    if degree < n {
        for j in 1..(n - degree) {
            let mut avg = 0.0;
            for param in params.iter().skip(j).take(degree) {
                avg += param;
            }
            knots[j + degree] = avg / degree as f64;
        }
    }
    knots
}

/// Solve the interpolation system for 2D points.
pub(crate) fn solve_interpolation_2d(
    params: &[f64],
    knots: &[f64],
    degree: usize,
    pts: &[DVec2],
) -> Vec<DVec2> {
    let n = pts.len();
    let a = collocation_matrix_2d(params, knots, degree, n, n);

    let rhs_x: Vec<f64> = pts.iter().map(|p| p.x).collect();
    let rhs_y: Vec<f64> = pts.iter().map(|p| p.y).collect();

    let cx = gauss_solve_2d(&a, &rhs_x);
    let cy = gauss_solve_2d(&a, &rhs_y);

    (0..n).map(|i| DVec2::new(cx[i], cy[i])).collect()
}

/// Build collocation matrix for B-spline interpolation.
fn collocation_matrix_2d(
    params: &[f64],
    knots: &[f64],
    degree: usize,
    n_data: usize,
    n_ctrl: usize,
) -> Vec<Vec<f64>> {
    params[..n_data]
        .iter()
        .map(|&t| all_basis_fns_2d(t, knots, degree, n_ctrl))
        .collect()
}

/// Find the knot span index.
fn find_span_2d(n_ctrl: usize, degree: usize, t: f64, knots: &[f64]) -> usize {
    let n = n_ctrl - 1;
    if t >= knots[n + 1] {
        return n;
    }
    if t <= knots[degree] {
        return degree;
    }
    let mut lo = degree;
    let mut hi = n + 1;
    let mut mid = (lo + hi) / 2;
    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

/// Evaluate all basis functions at parameter t.
fn basis_fns_2d(span: usize, t: f64, degree: usize, knots: &[f64]) -> Vec<f64> {
    let mut n = vec![0.0_f64; degree + 1];
    let mut left = vec![0.0_f64; degree + 1];
    let mut right = vec![0.0_f64; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0_f64;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Evaluate all n_ctrl basis functions at t (dense).
fn all_basis_fns_2d(t: f64, knots: &[f64], degree: usize, n_ctrl: usize) -> Vec<f64> {
    let span = find_span_2d(n_ctrl, degree, t, knots);
    let local = basis_fns_2d(span, t, degree, knots);
    let mut result = vec![0.0_f64; n_ctrl];
    for (k, &val) in local.iter().enumerate().take(degree + 1) {
        let idx = span - degree + k;
        if idx < n_ctrl {
            result[idx] = val;
        }
    }
    result
}

/// Gaussian elimination with partial pivoting.
fn gauss_solve_2d(a: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut mat: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.push(rhs[i]);
            r
        })
        .collect();

    for col in 0..n {
        let mut max_row = col;
        let mut max_val = mat[col][col].abs();
        for (row, row_data) in mat.iter().enumerate().skip(col + 1) {
            if row_data[col].abs() > max_val {
                max_val = row_data[col].abs();
                max_row = row;
            }
        }
        mat.swap(col, max_row);

        let pivot = mat[col][col];
        if pivot.abs() < TOLERANCE_FLOAT_LOOSE {
            continue;
        }

        for row in (col + 1)..n {
            let factor = mat[row][col] / pivot;
            let (lower, upper) = mat.split_at_mut(row);
            let pivot_row = &lower[col];
            let elim_row = &mut upper[0];
            for (elim_val, &pivot_val) in
                elim_row[col..=n].iter_mut().zip(pivot_row[col..=n].iter())
            {
                *elim_val -= pivot_val * factor;
            }
        }
    }

    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut sum = mat[i][n];
        for j in (i + 1)..n {
            sum -= mat[i][j] * x[j];
        }
        let diag = mat[i][i];
        x[i] = if diag.abs() > TOLERANCE_FLOAT_LOOSE {
            sum / diag
        } else {
            0.0
        };
    }
    x
}
