//! OCCT math_ThinPlateSpline — thin-plate spline surface from constraint points.
//!
//! (Moved here from math_utils.rs as part of OCCT MathPoly/GProp reorganization.)

use crate::geom::{BSplineSurface, Surface3};
use glam::DVec3;

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

    let m = n + 3;
    let mut mat = vec![0.0; m * m];
    let mut rhs = vec![0.0; m];

    for i in 0..n {
        for j in 0..n {
            let r = (points[i].truncate() - points[j].truncate()).length();
            mat[i * m + j] = tps_rbf(r);
        }
        mat[i * m + n + 0] = 1.0;
        mat[i * m + n + 1] = points[i].x;
        mat[i * m + n + 2] = points[i].y;
        mat[n * m + i] = 1.0;
        mat[(n + 1) * m + i] = points[i].x;
        mat[(n + 2) * m + i] = points[i].y;
        rhs[i] = points[i].z;
    }

    crate::math::lin::solve_linear_system(&mat, &rhs, m).map(|sol| {
        let w = sol[0..n].to_vec();
        let a = vec![sol[n], sol[n + 1], sol[n + 2]];
        (w, a)
    })
}

/// Evaluate a thin-plate spline at position (x, y).
pub fn evaluate_tps(x: f64, y: f64, w: &[f64], a: &[f64], points: &[DVec3]) -> f64 {
    let mut f = a[0] + a[1] * x + a[2] * y;
    for (i, &wi) in w.iter().enumerate() {
        let r = DVec2::new(x - points[i].x, y - points[i].y).length();
        f += wi * tps_rbf(r);
    }
    f
}

use glam::DVec2;

/// Build a plate surface (BSplineSurface) from constraint points.
pub fn build_plate_surface(constraints: &[DVec3], n_u: usize, n_v: usize) -> Option<Surface3> {
    if constraints.len() < 3 || n_u < 2 || n_v < 2 { return None; }

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

    let (w, a) = thin_plate_spline(constraints)?;

    let mut cp = vec![vec![DVec3::ZERO; n_u]; n_v];
    for vi in 0..n_v {
        let y = y_min + (y_max - y_min) * vi as f64 / (n_v - 1) as f64;
        for ui in 0..n_u {
            let x = x_min + (x_max - x_min) * ui as f64 / (n_u - 1) as f64;
            let z = evaluate_tps(x, y, &w, &a, constraints);
            cp[vi][ui] = DVec3::new(x, y, z);
        }
    }

    let degree = 3.min(n_u - 1).min(n_v - 1);
    let (knots_u, knots_v) = build_bspline_knots(n_u, n_v, degree);

    Some(Surface3::BSpline(BSplineSurface {
        degree_u: degree, degree_v: degree,
        knots_u, knots_v,
        control_points: cp,
        weights: vec![vec![1.0; n_u]; n_v],
    }))
}

/// Build clamped knot vectors for BSpline surface of given degree.
fn build_bspline_knots(n_u: usize, n_v: usize, degree: usize) -> (Vec<f64>, Vec<f64>) {
    let build = |n: usize| -> Vec<f64> {
        let nk = n + degree + 1;
        if n <= degree { return vec![0.0; nk]; }
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
