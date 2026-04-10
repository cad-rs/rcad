//! Analytic intersection of a sphere and a cylinder.
//!
//! # General theory
//!
//! The intersection of a sphere and a cylinder is generically a quartic space
//! curve (Viviani's curve and its generalisations).  Full analytic parametric
//! representation of the quartic is implemented as a Bernstein/rational form in
//! specialised software; here we handle the **axis-aligned special case** that
//! covers the vast majority of practical CAD situations.
//!
//! ## Axis-aligned case
//!
//! When the sphere centre **C** lies on the cylinder axis, the intersection
//! degenerates to one or two circles perpendicular to the axis:
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
//! ## Non-axis-aligned case
//!
//! When `d_perp > 0` (sphere centre is off the axis) the intersection is a
//! genuine quartic and we return [`SphereCylinderResult::General`], signalling
//! the caller to fall back to numeric marching.

use rcad_kernel::geom::{Circle3, CylindricalSurface, SphericalSurface};

use crate::tolerance::TOLERANCE_ABS;

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
    /// Sphere centre is off the cylinder axis — the intersection is a quartic
    /// space curve.  The caller should fall back to numeric marching.
    General,
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the analytic intersection of `sphere` and `cyl`.
///
/// Returns one of [`SphereCylinderResult`]'s variants:
///
/// - [`NoIntersection`](SphereCylinderResult::NoIntersection) — no real
///   intersection (axis-aligned case only; the `General` path always returns
///   [`General`](SphereCylinderResult::General)).
/// - [`TangentCircle`](SphereCylinderResult::TangentCircle) — one circle.
/// - [`TwoCircles`](SphereCylinderResult::TwoCircles) — two circles.
/// - [`General`](SphereCylinderResult::General) — quartic; fall back to
///   numeric marching.
///
/// The axis-aligned tolerance is ten times the absolute position tolerance.
pub fn intersect_sphere_cylinder(
    sphere: &SphericalSurface,
    cyl: &CylindricalSurface,
) -> SphereCylinderResult {
    let axis = cyl.axis.normalize();
    let d = sphere.center - cyl.origin;
    let d_along = d.dot(axis);
    let d_perp_vec = d - axis * d_along;
    let d_perp = d_perp_vec.length();

    // If the sphere centre is not on the cylinder axis, fall back to marching.
    if d_perp > TOLERANCE_ABS * 10.0 {
        return SphereCylinderResult::General;
    }

    // Axis-aligned case.
    let r = cyl.radius;
    let big_r = sphere.radius;

    // Discriminant: big_R² − r²
    let disc = big_r * big_r - r * r;

    if disc < -TOLERANCE_ABS {
        // Sphere too small to reach the cylinder surface along the axis.
        return SphereCylinderResult::NoIntersection;
    }

    // Height of sphere centre on the cylinder axis
    let h_c = d_along;

    if disc.abs() < TOLERANCE_ABS {
        // Tangent: one circle at h_c
        let center = cyl.origin + axis * h_c;
        return SphereCylinderResult::TangentCircle(Circle3 {
            center,
            normal: axis,
            radius: r,
        });
    }

    // Two distinct circles
    let delta_h = disc.sqrt();
    let h1 = h_c - delta_h;
    let h2 = h_c + delta_h;

    let c1 = Circle3 {
        center: cyl.origin + axis * h1,
        normal: axis,
        radius: r,
    };
    let c2 = Circle3 {
        center: cyl.origin + axis * h2,
        normal: axis,
        radius: r,
    };

    SphereCylinderResult::TwoCircles(c1, c2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;

    fn sphere(center: DVec3, radius: f64) -> SphericalSurface {
        SphericalSurface { center, axis: DVec3::Z, radius }
    }

    fn cyl(origin: DVec3, radius: f64) -> CylindricalSurface {
        CylindricalSurface { origin, axis: DVec3::Z, radius }
    }

    /// R > r, centre on axis → two circles
    #[test]
    fn two_circles_axis_aligned() {
        let sph = sphere(DVec3::new(0.0, 0.0, 3.0), 5.0);
        let c = cyl(DVec3::ZERO, 4.0);
        match intersect_sphere_cylinder(&sph, &c) {
            SphereCylinderResult::TwoCircles(c1, c2) => {
                // delta_h = sqrt(25 - 16) = 3
                // h1 = 3 - 3 = 0, h2 = 3 + 3 = 6
                assert!((c1.center.z - 0.0).abs() < 1e-9, "c1.z={}", c1.center.z);
                assert!((c2.center.z - 6.0).abs() < 1e-9, "c2.z={}", c2.center.z);
                assert!((c1.radius - 4.0).abs() < 1e-9);
                assert!((c2.radius - 4.0).abs() < 1e-9);
            }
            other => panic!("expected TwoCircles, got {other:?}"),
        }
    }

    /// R = r, centre on axis → tangent circle
    #[test]
    fn tangent_circle_equal_radii() {
        let sph = sphere(DVec3::new(0.0, 0.0, 5.0), 3.0);
        let c = cyl(DVec3::ZERO, 3.0);
        match intersect_sphere_cylinder(&sph, &c) {
            SphereCylinderResult::TangentCircle(tc) => {
                // h_c = 5, disc = 0, so circle at z = 5
                assert!((tc.center.z - 5.0).abs() < 1e-9, "tc.z={}", tc.center.z);
                assert!((tc.radius - 3.0).abs() < 1e-9);
            }
            other => panic!("expected TangentCircle, got {other:?}"),
        }
    }

    /// R < r, centre on axis → no intersection
    #[test]
    fn no_intersection_sphere_smaller() {
        let sph = sphere(DVec3::new(0.0, 0.0, 0.0), 1.0);
        let c = cyl(DVec3::ZERO, 2.0);
        assert!(matches!(
            intersect_sphere_cylinder(&sph, &c),
            SphereCylinderResult::NoIntersection
        ));
    }

    /// Centre NOT on axis → General
    #[test]
    fn general_off_axis() {
        // Sphere centre at (1, 0, 0) — far from Z axis (d_perp = 1)
        let sph = sphere(DVec3::new(1.0, 0.0, 0.0), 5.0);
        let c = cyl(DVec3::ZERO, 2.0);
        assert!(matches!(intersect_sphere_cylinder(&sph, &c), SphereCylinderResult::General));
    }
}
