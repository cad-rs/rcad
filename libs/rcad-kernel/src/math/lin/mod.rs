//! OCCT MathLin: linear algebra operations.
//!
//! Corresponds to OCCT `math_SVD`, `math_Gauss`, `math_Crout`,
//! `math_Householder`, `math_Matrix`, `math_LeastSquare`.

use glam::{DMat2, DMat3, DVec3};

const TOL_FLOAT_DEDUP: f64 = 1e-15;

// =============================================================================
// math_Matrix — matrix operations
// =============================================================================

/// Compute eigenvalues of a 2x2 matrix. Returns in descending order.
pub fn eigenvalues_2x2(m: DMat2) -> (f64, f64) {
    let trace = m.x_axis.x + m.y_axis.y;
    let det = m.determinant();
    let disc = trace * trace - 4.0 * det;
    if disc < 0.0 {
        let real_part = trace / 2.0;
        (real_part, real_part)
    } else {
        let sqrt_disc = disc.sqrt();
        let e1 = (trace + sqrt_disc) / 2.0;
        let e2 = (trace - sqrt_disc) / 2.0;
        if e1 >= e2 { (e1, e2) } else { (e2, e1) }
    }
}

/// Compute eigenvalues of a 3x3 matrix. Returns in descending order.
pub fn eigenvalues_3x3(m: DMat3) -> (f64, f64, f64) {
    let trace = m.x_axis.x + m.y_axis.y + m.z_axis.z;
    let det = m.determinant();
    let s = (m.y_axis.y * m.z_axis.z - m.y_axis.z * m.z_axis.y)
          + (m.x_axis.x * m.z_axis.z - m.x_axis.z * m.z_axis.x)
          + (m.x_axis.x * m.y_axis.y - m.x_axis.y * m.y_axis.x);
    let roots = crate::math::math_poly::solve_cubic(1.0, -trace, s, -det);
    match roots.len() {
        0 => (0.0, 0.0, 0.0),
        1 => (roots[0], roots[0], roots[0]),
        2 => {
            if roots[0] >= roots[1] { (roots[0], roots[1], roots[1]) }
            else { (roots[1], roots[0], roots[0]) }
        }
        3 => {
            let mut r = roots;
            r.sort_by(|a, b| b.partial_cmp(a).unwrap());
            (r[0], r[1], r[2])
        }
        _ => (0.0, 0.0, 0.0),
    }
}

/// Compute the inverse of a 3x3 matrix. Returns None if singular.
pub fn inverse_3x3(m: DMat3) -> Option<DMat3> {
    let det = m.determinant();
    if det.abs() < TOL_FLOAT_DEDUP { return None; }
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
// math_SVD — Singular Value Decomposition
// =============================================================================

/// Solve A*x = b for 3x3 matrix A using SVD-based pseudo-inverse when singular.
pub fn svd_solve_3x3(a: DMat3, b: DVec3) -> Option<DVec3> {
    if let Some(inv) = inverse_3x3(a) {
        return Some(inv * b);
    }
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
fn svd_jacobi_3x3(mut a: DMat3) -> Option<(DMat3, DVec3, DMat3)> {
    let tol = 1e-14;
    let mut v = DMat3::IDENTITY;
    for _ in 0..100 {
        let mut max_gamma = 0.0;
        let mut best_p = 0;
        let mut best_q = 1;
        for p in 0..3 {
            for q in (p + 1)..3 {
                let gamma = a.col(p).dot(a.col(q));
                let alpha = a.col(p).length_squared();
                let beta = a.col(q).length_squared();
                let denom = (alpha * beta).sqrt();
                if denom > tol {
                    let ratio = gamma.abs() / denom;
                    if ratio > max_gamma {
                        max_gamma = ratio; best_p = p; best_q = q;
                    }
                }
            }
        }
        if max_gamma < tol { break; }
        let (p, q) = (best_p, best_q);
        let alpha = a.col(p).length_squared();
        let beta = a.col(q).length_squared();
        let gamma = a.col(p).dot(a.col(q));
        let tau = (beta - alpha) / (2.0 * gamma);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let ap = a.col(p);
        let aq = a.col(q);
        let new_ap = ap * c + aq * s;
        let new_aq = -ap * s + aq * c;
        a = DMat3::from_cols(
            if p == 0 { new_ap } else if q == 0 { new_aq } else { a.col(0) },
            if p == 1 { new_ap } else if q == 1 { new_aq } else { a.col(1) },
            if p == 2 { new_ap } else if q == 2 { new_aq } else { a.col(2) },
        );
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
    let mut s = DVec3::new(a.col(0).length(), a.col(1).length(), a.col(2).length());
    let tol_s = 1e-14;
    if s.x < tol_s { s.x = 0.0; }
    if s.y < tol_s { s.y = 0.0; }
    if s.z < tol_s { s.z = 0.0; }
    let u0 = if s.x > tol_s { a.col(0) / s.x } else { DVec3::X };
    let u1 = if s.y > tol_s { a.col(1) / s.y } else { DVec3::Y };
    let u2 = if s.z > tol_s { a.col(2) / s.z } else { DVec3::Z };
    let u = DMat3::from_cols(u0, u1, u2);
    Some((u, s, v.transpose()))
}

// =============================================================================
// math_Gauss — Gaussian elimination
// =============================================================================

/// Solve A*x = b for small n×n system using Gaussian elimination.
pub fn solve_linear_system(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }
    for col in 0..n {
        let mut max_row = col;
        for row in (col + 1)..n {
            if aug[row * (n + 1) + col].abs() > aug[max_row * (n + 1) + col].abs() {
                max_row = row;
            }
        }
        if aug[max_row * (n + 1) + col].abs() < TOL_FLOAT_DEDUP { return None; }
        if max_row != col {
            for j in col..=n { aug.swap(col * (n + 1) + j, max_row * (n + 1) + j); }
        }
        for row in (col + 1)..n {
            let factor = aug[row * (n + 1) + col] / aug[col * (n + 1) + col];
            for j in col..=n { aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j]; }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i * (n + 1) + n];
        for j in (i + 1)..n { sum -= aug[i * (n + 1) + j] * x[j]; }
        x[i] = sum / aug[i * (n + 1) + i];
    }
    Some(x)
}

// =============================================================================
// math_Crout — LU decomposition with partial pivoting
// =============================================================================

/// Solve A*x = b using Crout LU decomposition.
pub fn crout_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    if n == 0 { return None; }
    let mut lu = a.to_vec();
    let mut p: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut max_val = lu[p[k] * n + k].abs();
        let mut max_r = k;
        for i in (k + 1)..n {
            let v = lu[p[i] * n + k].abs();
            if v > max_val { max_val = v; max_r = i; }
        }
        p.swap(k, max_r);
        if lu[p[k] * n + k].abs() < TOL_FLOAT_DEDUP { return None; }
        for i in (k + 1)..n {
            let f = lu[p[i] * n + k] / lu[p[k] * n + k];
            lu[p[i] * n + k] = f;
            for j in (k + 1)..n { lu[p[i] * n + j] -= f * lu[p[k] * n + j]; }
        }
    }
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[p[i]];
        for j in 0..i { s -= lu[p[i] * n + j] * y[j]; }
        y[i] = s;
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for j in (i + 1)..n { s -= lu[p[i] * n + j] * x[j]; }
        if lu[p[i] * n + i].abs() < TOL_FLOAT_DEDUP { return None; }
        x[i] = s / lu[p[i] * n + i];
    }
    Some(x)
}

// =============================================================================
// math_Householder — QR via Householder reflections
// =============================================================================

/// Solve A*x = b using Householder QR decomposition.
pub fn householder_solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    if n == 0 { return None; }
    let mut r = a.to_vec();
    let mut rhs = b.to_vec();
    for k in 0..n - 1 {
        let mut nrm2 = 0.0;
        for i in k..n { nrm2 += r[i * n + k] * r[i * n + k]; }
        if nrm2 < TOL_FLOAT_DEDUP { continue; }
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
        if r[i * n + i].abs() < TOL_FLOAT_DEDUP { return None; }
        x[i] = s / r[i * n + i];
    }
    Some(x)
}

// =============================================================================
// math_LeastSquare — least squares fitting
// =============================================================================

/// Fit a line y = a + b*x to data points using linear least squares.
pub fn least_squares_linear(x: &[f64], y: &[f64]) -> Option<(f64, f64)> {
    let n = x.len().min(y.len());
    if n < 2 { return None; }
    let mut sx = 0.0; let mut sy = 0.0;
    let mut sxx = 0.0; let mut sxy = 0.0;
    for i in 0..n {
        sx += x[i]; sy += y[i];
        sxx += x[i] * x[i]; sxy += x[i] * y[i];
    }
    let denom = (n as f64) * sxx - sx * sx;
    if denom.abs() < TOL_FLOAT_DEDUP { return None; }
    let b = ((n as f64) * sxy - sx * sy) / denom;
    let a = (sy - b * sx) / (n as f64);
    Some((a, b))
}
