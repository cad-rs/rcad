//! OCCT MathSys: nonlinear system of equations solvers.
//!
//! Corresponds to OCCT `math_FunctionSetRoot` / `MathSys`.
//! Functions: newton_2d, newton_3d.

use glam::{DMat2, DMat3, DVec2, DVec3};

const TOL_FLOAT_DEDUP: f64 = 1e-15;

// =============================================================================
// MathSys_Newton2D — Newton-Raphson for 2D systems
// =============================================================================

/// Newton-Raphson for 2D systems. Solves F(x) = 0 where F: R^2 -> R^2.
pub fn newton_2d(
    f: fn(DVec2) -> DVec2,
    jacobian: fn(DVec2) -> DMat2,
    x0: DVec2,
    tol: f64,
) -> Option<DVec2> {
    let mut x = x0;
    for _ in 0..50 {
        let fx = f(x);
        if fx.length() < tol { return Some(x); }
        let j = jacobian(x);
        let det = j.determinant();
        if det.abs() < TOL_FLOAT_DEDUP { return None; }
        let delta = j.inverse() * fx;
        let x_new = x - delta;
        if delta.length() < tol { return Some(x_new); }
        x = x_new;
    }
    if f(x).length() < tol { Some(x) } else { None }
}

// =============================================================================
// MathSys_Newton3D — Newton-Raphson for 3D systems
// =============================================================================

fn inverse_3x3(m: DMat3) -> Option<DMat3> {
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

/// Newton-Raphson for 3D systems. Solves F(x) = 0 where F: R^3 -> R^3.
pub fn newton_3d(
    f: fn(DVec3) -> DVec3,
    jacobian: fn(DVec3) -> DMat3,
    x0: DVec3,
    tol: f64,
) -> Option<DVec3> {
    let mut x = x0;
    for _ in 0..50 {
        let fx = f(x);
        if fx.length() < tol { return Some(x); }
        let j = jacobian(x);
        let det = j.determinant();
        if det.abs() < TOL_FLOAT_DEDUP { return None; }
        if let Some(j_inv) = inverse_3x3(j) {
            let delta = j_inv * fx;
            let x_new = x - delta;
            if delta.length() < tol { return Some(x_new); }
            x = x_new;
        } else { return None; }
    }
    if f(x).length() < tol { Some(x) } else { None }
}
