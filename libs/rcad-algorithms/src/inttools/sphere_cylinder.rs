//! Analytic intersection of a sphere and a cylinder.
//!
//! # Case classification
//!
//! ## Axis-aligned case (sphere centre on cylinder axis)
//!
//! When the sphere centre **C** lies on the cylinder axis (`d_perp ≈ 0`), the
//! intersection degenerates to one or two circles perpendicular to the axis:
//!
//! ```text
//! h = h_c ± √(R² − r²)
//! ```
//!
//! where
//! - `h_c = (C − O) · â`  (height of sphere centre above cylinder origin),
//! - `R` = sphere radius,
//! - `r` = cylinder radius,
//! - `â` = unit cylinder axis,
//! - `h` = height of the intersection circle on the cylinder axis.
//!
//! Each such `h` yields a circle of radius `r` centred at `O + h · â` with
//! normal `â`.
//!
//! ## Parallel-axis offset case
//!
//! When the sphere centre is off-axis but the sphere and cylinder axes are
//! parallel (or the cylinder has no preferred axis direction relative to the
//! sphere), we can still decide:
//!
//! - Let `d` = perpendicular distance from sphere centre to cylinder axis.
//! - The sphere surface is at radial distances `[d − R, d + R]` from the axis.
//! - The cylinder surface is at radial distance `r` from the axis.
//!
//! Therefore:
//! - If `d − R > r` or `d + R < r` (and `d > r` for the latter): **no intersection**.
//! - If `|d − r| ≤ R`: the sphere surface intersects the cylinder surface;
//!   the exact intersection is a quartic (Viviani-type) curve — return `General`.
//!
//! ## General / skew case
//!
//! For arbitrary off-axis configurations the intersection is a quartic space
//! curve (Viviani-type).  We solve it analytically by substituting the
//! cylinder parametrisation into the sphere equation, yielding a quadratic
//! in `v` for each cylinder azimuth `u` — see [`intersect_skew_sphere_cylinder`].

use glam::DVec3;
use rcad_kernel::geom::{any_perpendicular, Circle3, CylindricalSurface, SphericalSurface};
use rcad_kernel::SurfaceEval;
use std::f64::consts::TAU;

use crate::tolerance::*;
use super::pcurve_derive::refine_polyline;

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Analytic result of sphere × cylinder intersection.
#[derive(Debug, Clone)]
pub enum SphereCylinderResult {
    /// Sphere and cylinder do not intersect (disjoint or one fully inside the
    /// other when the axis-aligned condition also holds).
    NoIntersection,
    /// Exactly one tangent circle (R = r and sphere centre on axis, or the two
    /// roots coincide).
    TangentCircle(Circle3),
    /// Two distinct intersection circles (axis-aligned, `R > r`).
    TwoCircles(Circle3, Circle3),
    /// The intersection is a quartic space curve.  The caller should fall back
    /// to numeric marching.
    General,
    /// Skew (off-axis) configuration solved analytically via cylinder-azimuth
    /// sampling.  Each inner Vec is a polyline branch of the intersection curve
    /// in 3D (at most two branches, from the ± sqrt of the quadratic).
    SkewQuartic(Vec<Vec<DVec3>>),
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `sphere` and `cyl`.
///
/// Returns one of [`SphereCylinderResult`]'s variants:
///
/// - [`NoIntersection`](SphereCylinderResult::NoIntersection) — disjoint.
/// - [`TangentCircle`](SphereCylinderResult::TangentCircle) — one tangent circle
///   (axis-aligned case, discriminant = 0).
/// - [`TwoCircles`](SphereCylinderResult::TwoCircles) — two circles (axis-aligned).
/// - [`General`](SphereCylinderResult::General) — quartic; fall back to marching.
/// - [`SkewQuartic`](SphereCylinderResult::SkewQuartic) — analytic quartic branches
///   (off-axis, solved via cylinder parametrisation).
///
/// The axis-aligned tolerance is ten times the absolute position tolerance.
pub fn intersect_sphere_cylinder(
    sphere: &SphericalSurface,
    cyl: &CylindricalSurface,
) -> SphereCylinderResult {
    intersect_sphere_cylinder_with_tolerance(sphere, cyl, 0.0)
}

/// Compute sphere-cylinder intersection with additional fuzzy tolerance.
///
/// This relaxes axis-aligned and distance early-out checks by `fuzzy_tol` so
/// near-coincident cases can still classify into analytic branches.
pub fn intersect_sphere_cylinder_with_tolerance(
    sphere: &SphericalSurface,
    cyl: &CylindricalSurface,
    fuzzy_tol: f64,
) -> SphereCylinderResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);
    let axis = cyl.axis.normalize();
    let d = sphere.center - cyl.origin;
    let d_along = d.dot(axis);
    let d_perp_vec = d - axis * d_along;
    let d_perp = d_perp_vec.length(); // perpendicular distance: sphere centre → cyl axis

    let r = cyl.radius;
    let big_r = sphere.radius;

    // ── Axis-aligned case ─────────────────────────────────────────────────────
    if d_perp < tol * 10.0 {
        // Sphere centre is on (or extremely close to) the cylinder axis.
        let disc = big_r * big_r - r * r;

        if disc < -tol {
            return SphereCylinderResult::NoIntersection;
        }

        let h_c = d_along;

        if disc.abs() < tol {
            let center = cyl.origin + axis * h_c;
            return SphereCylinderResult::TangentCircle(Circle3::new(center, axis, r));
        }

        let delta_h = disc.sqrt();
        let c1 = Circle3::new(cyl.origin + axis * (h_c - delta_h), axis, r,
        );
        let c2 = Circle3::new(cyl.origin + axis * (h_c + delta_h), axis, r,
        );
        return SphereCylinderResult::TwoCircles(c1, c2);
    }

    // ── Off-axis: early-out distance test ─────────────────────────────────────
    //
    // The cylinder lateral surface is everywhere at distance `r` from the axis.
    // The sphere surface spans radial distances [d_perp − R, d_perp + R] from
    // the axis (considering all points on the sphere surface).
    //
    // No intersection when:
    //   (a) d_perp - R > r  →  sphere is entirely outside the cylinder (far side)
    //   (b) d_perp + R < r  →  sphere is entirely inside the cylinder (near side)
    //       but only when the sphere is smaller than the cylinder radius + offset
    //
    // Case (a): closest radial approach of sphere to axis exceeds cylinder radius.
    if d_perp - big_r > r + tol {
        return SphereCylinderResult::NoIntersection;
    }
    // Case (b): sphere is fully enclosed inside the cylinder laterally.
    if d_perp + big_r < r - tol {
        return SphereCylinderResult::NoIntersection;
    }

    // ── Quartic (Viviani-type) intersection ───────────────────────────────────
    // Try analytic quartic solver first; fall back to General (marching) if
    // it returns no branches (e.g. near-tangent configurations).
    let skew_result = intersect_skew_sphere_cylinder(sphere, cyl);
    if !skew_result.is_empty() {
        SphereCylinderResult::SkewQuartic(skew_result)
    } else {
        SphereCylinderResult::General
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Skew-axis analytic solver
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the intersection of a sphere and cylinder with skew (off-axis)
/// configuration using an analytic quartic solver.
///
/// # Theory
///
/// Cylinder parametrisation (u = azimuth [0, 2π), v = height along axis):
///
/// ```text
/// P(u,v) = O_cyl + v·a_cyl + r_cyl·(cos(u)·x_cyl + sin(u)·y_cyl)
/// ```
///
/// The sphere surface satisfies:
///
/// ```text
/// |P - O_sph|² = R²
/// ```
///
/// Substituting P(u,v) gives a quadratic in v for each fixed u:
///
/// ```text
/// a_v·v² + b_v(u)·v + c_v(u) = 0
///
/// a_v   = 1                                              (always, |a_cyl| = 1)
/// b_v(u) = 2·(D0·a_cyl)
/// c_v(u) = |D0|² - R²
/// D0(u)  = O_cyl - O_sph + r_cyl·(cos(u)·x_cyl + sin(u)·y_cyl)
/// ```
///
/// Since a_v = 1, the quadratic never degenerates.  For each u we compute
/// `v = (-b_v ± sqrt(b_v² - 4·c_v)) / 2`, giving up to two branches.
fn intersect_skew_sphere_cylinder(
    sphere: &SphericalSurface,
    cyl: &CylindricalSurface,
) -> Vec<Vec<DVec3>> {
    let a_cyl = cyl.axis.normalize();
    let o_cyl = cyl.origin;
    let o_sph = sphere.center;
    let r_cyl = cyl.radius;
    let r_sph = sphere.radius;

    // Perpendicular basis for cylinder (must match CylindricalSurface::point_at).
    let x_cyl = any_perpendicular(a_cyl);
    let y_cyl = a_cyl.cross(x_cyl).normalize();

    let delta_o = o_cyl - o_sph; // O_cyl - O_sph

    const N_SAMPLES: usize = 128;
    const CHORD_TOL: f64 = crate::inttools::CHORD_TOLERANCE;
    const REFINE_DEPTH: usize = crate::inttools::CHORD_REFINE_DEPTH;
    let mut branch_plus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);
    let mut branch_minus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cos_u, sin_u) = (u.cos(), u.sin());

        // D0(u) = (O_cyl - O_sph) + r_cyl·(cos(u)·x_cyl + sin(u)·y_cyl)
        let d0 = delta_o + r_cyl * (cos_u * x_cyl + sin_u * y_cyl);

        // b_v(u) = 2·(D0·a_cyl)
        let b_v = 2.0 * d0.dot(a_cyl);

        // c_v(u) = |D0|² - R²
        let c_v = d0.length_squared() - r_sph * r_sph;

        // v² + b_v·v + c_v = 0  (a_v = 1, always non-degenerate)
        let disc = b_v * b_v - 4.0 * c_v;

        if disc < 0.0 {
            // No intersection at this u azimuth.
            continue;
        }

        let sqrt_disc = disc.sqrt();

        let v_plus = (-b_v + sqrt_disc) * 0.5;
        let v_minus = (-b_v - sqrt_disc) * 0.5;

        if v_plus.is_finite() {
            let p = cyl.point_at(u, v_plus);
            if p.is_finite() {
                branch_plus.push((u, p));
            }
        }
        if v_minus.is_finite() {
            let p = cyl.point_at(u, v_minus);
            if p.is_finite() {
                branch_minus.push((u, p));
            }
        }
    }

    // Adaptive refinement: subdivide segments where chord error exceeds tolerance
    let eval_plus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d0 = delta_o + r_cyl * (cos_u * x_cyl + sin_u * y_cyl);
        let b_v = 2.0 * d0.dot(a_cyl);
        let c_v = d0.length_squared() - r_sph * r_sph;
        let disc = b_v * b_v - 4.0 * c_v;
        if disc < 0.0 { return None; }
        let v = (-b_v + disc.sqrt()) * 0.5;
        if v.is_finite() { let p = cyl.point_at(u, v); if p.is_finite() { return Some(p); } }
        None
    };
    let eval_minus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d0 = delta_o + r_cyl * (cos_u * x_cyl + sin_u * y_cyl);
        let b_v = 2.0 * d0.dot(a_cyl);
        let c_v = d0.length_squared() - r_sph * r_sph;
        let disc = b_v * b_v - 4.0 * c_v;
        if disc < 0.0 { return None; }
        let v = (-b_v - disc.sqrt()) * 0.5;
        if v.is_finite() { let p = cyl.point_at(u, v); if p.is_finite() { return Some(p); } }
        None
    };

    let (mut branch_plus, mut branch_minus): (Vec<DVec3>, Vec<DVec3>) = (
        refine_polyline(&branch_plus, eval_plus, CHORD_TOL, REFINE_DEPTH)
            .into_iter().map(|(_, p)| p).collect(),
        refine_polyline(&branch_minus, eval_minus, CHORD_TOL, REFINE_DEPTH)
            .into_iter().map(|(_, p)| p).collect(),
    );

    // Dedup: remove trailing points that nearly duplicate the first point
    // (closed curve degeneracy).
    let dedup = |pts: &mut Vec<DVec3>| {
        while pts.len() >= 3 {
            let n = pts.len();
            if (pts[n - 1] - pts[0]).length_squared() < TOLERANCE_VEC_SQ_MIN {
                pts.pop();
            } else {
                break;
            }
        }
    };

    let mut branches = Vec::new();
    if branch_plus.len() >= 2 {
        dedup(&mut branch_plus);
        if branch_plus.len() >= 2 {
            branches.push(branch_plus);
        }
    }
    if branch_minus.len() >= 2 {
        // Check the minus branch is distinct from the plus branch.  If they're
        // nearly the same (tangent intersection), keep only one.
        let is_distinct = branches.is_empty()
            || (branch_minus[0] - branches[0][0]).length_squared() > TOLERANCE_VEC_SQ_MIN;
        if is_distinct {
            dedup(&mut branch_minus);
            if branch_minus.len() >= 2 {
                branches.push(branch_minus);
            }
        }
    }

    branches
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────


