//! Analytic intersection of two cylinders.
//!
//! # Case classification
//!
//! ## Parallel axes
//!
//! When the two cylinder axes are parallel (cross-product length 鈮?0):
//!
//! - **Coaxial**: axes coincide 鈫?intersection is the full cylinder surface of
//!   the smaller radius (degenerate; we return `Coaxial`).
//! - **Offset parallel**: axes are parallel but distinct.  The gap between the
//!   two axis lines is `d`.
//!   - `d 鈮?r1 + r2`: no intersection.
//!   - `d = r1 + r2` (within tolerance): external tangent, one generator line.
//!   - `|r1 鈭?r2| < d < r1 + r2`: two generator lines (cross-section chords).
//!   - `d = |r1 鈭?r2|`: internal tangent, one generator line.
//!   - `d < |r1 鈭?r2|`: one cylinder inside the other, no surface intersection.
//!
//! ## Perpendicular axes (Steinmetz intersection)
//!
//! When the axes are perpendicular and the cross-section distance equals zero
//! (axes actually intersect), the intersection curves are two ellipses 鈥?//! specifically the classic Steinmetz configuration.  We return
//! `Perpendicular(TwoEllipses(...))` for this sub-case.
//!
//! ## General skew axes
//!
//! For all other orientations we return `General`, signalling the caller to
//! fall back to numeric marching.

use glam::DVec3;
use rcad_kernel::SurfaceEval;
use rcad_kernel::geom::{Circle3, CylindricalSurface, Ellipse3, Line3, any_perpendicular};
use std::f64::consts::TAU;

use super::pcurve_derive::refine_polyline;


// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Result type
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Analytic result of cylinder 脳 cylinder intersection.
#[derive(Debug, Clone)]
pub enum CylinderCylinderResult {
    /// The cylinders do not intersect (disjoint or fully nested).
    NoIntersection,
    /// Cylinders are coaxial; the intersection is the full lateral surface of
    /// the smaller cylinder.
    Coaxial,
    /// Parallel axes: the intersection is exactly one generator line (external
    /// or internal tangent).
    OneGeneratorLine(Line3),
    /// Parallel axes: the intersection consists of two generator lines.
    TwoGeneratorLines(Line3, Line3),
    /// Perpendicular intersecting axes (Steinmetz): two ellipses.
    TwoEllipses(Ellipse3, Ellipse3),
    /// Perpendicular intersecting axes, equal radii: two circles.
    TwoCircles(Circle3, Circle3),
    /// Perpendicular axes with non-zero offset (axes do not intersect),
    /// but cylinders overlap (dist 鈮?r1 + r2). Produces one or two
    /// space curves parametrized on cyl1's surface:
    /// `P(胃) = O1 + v*A1 + R1*(cos(胃)*U1 + sin(胃)*V1)`
    /// where `v = dz 卤 sqrt(R2虏 - (dy + R1*sin(胃-伪))虏)`.
    PerpendicularOffsetCurves {
        cyl1: CylindricalSurface,
        cyl2: CylindricalSurface,
        /// Distance between axes
        dist: f64,
    },
    /// Skew axes (non-parallel, non-perpendicular): analytic quartic solution.
    ///
    /// The two cylinders intersect in a quartic space curve.  For each cylinder
    /// azimuth u 鈭?[0, 2蟺) the second cylinder's equation reduces to a quadratic
    /// in the height v, solved analytically.  Two branches (卤 sqrt) are returned
    /// as polylines.
    SkewQuartic(Vec<Vec<DVec3>>),
    /// General case (skew axes or oblique angle not handled analytically).
    /// The caller should fall back to numeric marching.
    General,
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Main function
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the analytic intersection of `cyl1` and `cyl2`.
pub fn intersect_cylinder_cylinder(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
) -> CylinderCylinderResult {
    intersect_cylinder_cylinder_with_eps(cyl1, cyl2, rcad_kernel::rcad_kernel::CONFUSION, rcad_kernel::ANGULAR)
}

/// Compute cylinder-cylinder intersection with additional fuzzy tolerance.
pub fn intersect_cylinder_cylinder_with_tolerance(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    fuzzy_tol: f64,
) -> CylinderCylinderResult {
    let linear_tol = rcad_kernel::rcad_kernel::CONFUSION + fuzzy_tol.max(0.0);
    let angular_tol = rcad_kernel::ANGULAR + fuzzy_tol.max(0.0);
    intersect_cylinder_cylinder_with_eps(cyl1, cyl2, linear_tol, angular_tol)
}

fn intersect_cylinder_cylinder_with_eps(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    linear_tol: f64,
    angular_tol: f64,
) -> CylinderCylinderResult {
    let a1 = cyl1.axis.normalize();
    let a2 = cyl2.axis.normalize();

    let cross = a1.cross(a2);
    let sin_angle = cross.length(); // |sin 胃|

    // 鈹€鈹€ Parallel axes 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    if sin_angle < angular_tol {
        return intersect_parallel_cylinders(cyl1, cyl2, a1, linear_tol);
    }

    // 鈹€鈹€ Perpendicular axes 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    let cos_angle = a1.dot(a2).abs();
    if cos_angle < angular_tol {
        return intersect_perpendicular_cylinders(cyl1, cyl2, a1, a2, linear_tol);
    }

    // 鈹€鈹€ Skew axes (analytic quartic solver) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    // Parameterize cyl1, substitute into cyl2 equation 鈫?quadratic in v per u.
    let skew_result = intersect_skew_cylinder_cylinder(cyl1, cyl2);
    if !skew_result.is_empty() {
        return CylinderCylinderResult::SkewQuartic(skew_result);
    }

    // 鈹€鈹€ General / oblique (fallback to marching) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    CylinderCylinderResult::General
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Parallel axes
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn intersect_parallel_cylinders(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    axis: DVec3,
    linear_tol: f64,
) -> CylinderCylinderResult {
    let r1 = cyl1.radius;
    let r2 = cyl2.radius;

    // Perpendicular distance between the two parallel axes.
    let delta = cyl2.origin - cyl1.origin;
    let delta_perp = delta - axis * delta.dot(axis);
    let d = delta_perp.length();

    // Coaxial check
    if d < linear_tol {
        if (r1 - r2).abs() < linear_tol {
            return CylinderCylinderResult::Coaxial;
        }
        // One inside the other along the same axis
        return CylinderCylinderResult::NoIntersection;
    }

    let sum = r1 + r2;
    let diff = (r1 - r2).abs();

    if d > sum + linear_tol {
        return CylinderCylinderResult::NoIntersection;
    }
    if d < diff - linear_tol {
        // One cylinder fully inside the other
        return CylinderCylinderResult::NoIntersection;
    }

    // Direction from cyl1 axis to cyl2 axis (perpendicular)
    let dir_perp = delta_perp.normalize();

    // External tangent
    if (d - sum).abs() < linear_tol {
        let point = cyl1.origin + dir_perp * r1;
        return CylinderCylinderResult::OneGeneratorLine(Line3 {
            origin: point,
            direction: axis,
        });
    }
    // Internal tangent
    if (d - diff).abs() < linear_tol {
        // The tangent line is on the side of the smaller cylinder that is
        // closest to the larger cylinder's axis.
        let point = if r1 >= r2 {
            cyl1.origin + dir_perp * r1
        } else {
            cyl1.origin - dir_perp * r1
        };
        return CylinderCylinderResult::OneGeneratorLine(Line3 {
            origin: point,
            direction: axis,
        });
    }

    // Two generator lines: find the two intersection points in the
    // perpendicular cross-section.
    //
    // In 2D: circle1 centred at origin radius r1, circle2 centred at (d, 0)
    // radius r2.  The intersection x-coordinate:
    //   x = (d虏 + r1虏 - r2虏) / (2d)
    //   y = 卤sqrt(r1虏 - x虏)
    let x = (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let y_sq = r1 * r1 - x * x;
    if y_sq < 0.0 {
        return CylinderCylinderResult::NoIntersection;
    }
    let y = y_sq.sqrt();

    // Orthogonal unit vector in the cross-section plane
    let v_perp = axis.cross(dir_perp).normalize();

    let p1 = cyl1.origin + dir_perp * x + v_perp * y;
    let p2 = cyl1.origin + dir_perp * x - v_perp * y;

    CylinderCylinderResult::TwoGeneratorLines(
        Line3 {
            origin: p1,
            direction: axis,
        },
        Line3 {
            origin: p2,
            direction: axis,
        },
    )
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Perpendicular intersecting axes (Steinmetz)
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

fn intersect_perpendicular_cylinders(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    a1: DVec3,
    a2: DVec3,
    linear_tol: f64,
) -> CylinderCylinderResult {
    let r1 = cyl1.radius;
    let r2 = cyl2.radius;

    // Find the closest point between the two axes (skew lines).
    // Parametric form: P = O1 + t*a1  and  Q = O2 + s*a2
    // The connecting vector at closest approach is perpendicular to both axes.
    let w = cyl1.origin - cyl2.origin;
    let b = a1.dot(a2);
    let denom = 1.0 - b * b;

    if denom.abs() < 1e-12 {
        // Degenerate (should not reach here since we checked perpendicularity)
        return CylinderCylinderResult::General;
    }

    let d1 = a1.dot(w);
    let d2 = a2.dot(w);
    let t = (b * d2 - d1) / denom;
    let s = (d2 - b * d1) / denom;

    let closest1 = cyl1.origin + a1 * t;
    let closest2 = cyl2.origin + a2 * s;

    // Perpendicular distance between axes
    let dist = (closest1 - closest2).length();

    if dist > r1 + r2 + linear_tol {
        return CylinderCylinderResult::NoIntersection;
    }

    // For the Steinmetz case the axes must actually cross (dist 鈮?0).
    // For larger dist we can still handle overlapping perpendicular cylinders
    // analytically using a 胃-parametrization on cyl1's surface.
    if dist > linear_tol * 10.0 {
        if dist <= r1 + r2 + linear_tol {
            return CylinderCylinderResult::PerpendicularOffsetCurves {
                cyl1: *cyl1,
                cyl2: *cyl2,
                dist,
            };
        }
        return CylinderCylinderResult::NoIntersection;
    }

    // Intersection point of the two axes
    let origin = (closest1 + closest2) * 0.5;

    // Third axis = a1 脳 a2  (normal to both, the "viewing" direction)
    let _a3 = a1.cross(a2).normalize();

    if (r1 - r2).abs() < linear_tol {
        // Equal radii: the Steinmetz intersection of two perpendicular cylinders
        // produces two ellipses, each lying in a diagonal plane.
        //
        // Derivation: from Cyl1 eq |P 脳 a1| = r1 and Cyl2 eq |P 脳 a2| = r2,
        // with a1路a2 = 0 and r1 = r2 = r, we get (P路(a1+a2))(P路(a1-a2)) = 0.
        // So P lies in plane n1路P = 0 or n2路P = 0 where:
        //   n1 = (a1 + a2).normalize()
        //   n2 = (a1 - a2).normalize()
        //
        // In each plane, the intersection with Cyl1 is an ellipse with:
        //   major_radius = r路鈭? (along direction (a1鈭揳2) in the plane)
        //   minor_radius = r   (along the direction normal脳major_dir)
        let sqrt2 = std::f64::consts::SQRT_2;
        let n1 = (a1 + a2).normalize();
        let n2 = (a1 - a2).normalize();
        let u1 = (a1 - a2).normalize();
        let u2 = (a1 + a2).normalize();

        let ellipse1 = Ellipse3 {
            center: origin,
            normal: n1,
            major_dir: u1,
            major_radius: r1 * sqrt2,
            minor_radius: r1,
        };
        let ellipse2 = Ellipse3 {
            center: origin,
            normal: n2,
            major_dir: u2,
            major_radius: r1 * sqrt2,
            minor_radius: r1,
        };
        return CylinderCylinderResult::TwoEllipses(ellipse1, ellipse2);
    }

    // Unequal radii: the Steinmetz two-ellipse derivation only works
    // when r1 == r2 (the difference of squares x虏 鈭?z虏 = r1虏 鈭?r2虏
    // factorises into planes n1路P = 0 and n2路P = 0).  For unequal
    // radii the intersection curves are general space curves 鈥?not
    // planar ellipses 鈥?but we can still use the analytic
    // 胃-parametrization of PerpendicularOffsetCurves (same formula
    // works for both offset and intersecting axes).
    if dist <= r1 + r2 + linear_tol {
        return CylinderCylinderResult::PerpendicularOffsetCurves {
            cyl1: *cyl1,
            cyl2: *cyl2,
            dist,
        };
    }
    CylinderCylinderResult::NoIntersection
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Skew-axis analytic solver
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Compute the intersection of two cylinders with skew (non-parallel,
/// non-perpendicular) axes using an analytic quartic solver.
///
/// # Theory
///
/// Cylinder 1 parametrisation (u = azimuth [0, 2蟺), v = height along axis):
///
/// ```text
/// P(u,v) = O1 + v路a1 + r1路(cos(u)路x1 + sin(u)路y1)
/// ```
///
/// Cylinder 2 implicit equation (for any point P, distance from axis a2 is r2):
///
/// ```text
/// |P - O2|虏 - ((P - O2)路a2)虏 = r2虏
/// ```
///
/// Substituting P(u,v) gives F(v) = 0 where
///
/// ```text
/// a_v路v虏 + b_v(u)路v + c_v(u) = 0
///
/// a_v   = 1 - (a1路a2)虏                                           (constant)
/// b_v(u) = 2路(D0路a1) - 2路(D0路a2)路(a1路a2)
/// c_v(u) = |D0|虏 - (D0路a2)虏 - r2虏
/// D0(u)  = O1 - O2 + r1路(cos(u)路x1 + sin(u)路y1)
/// ```
///
/// For each u we solve the quadratic for v, giving two branches (卤 sqrt).
fn intersect_skew_cylinder_cylinder(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
) -> Vec<Vec<DVec3>> {
    let a1 = cyl1.axis.normalize();
    let a2 = cyl2.axis.normalize();
    let o1 = cyl1.origin;
    let o2 = cyl2.origin;
    let r1 = cyl1.radius;
    let r2_sq = cyl2.radius * cyl2.radius;

    // Perpendicular basis for cyl1 (must match CylindricalSurface::point_at).
    let x1 = any_perpendicular(a1);
    let y1 = a1.cross(x1).normalize();

    // Constant coefficient a_v = 1 - (a1路a2)虏.
    let a1_dot_a2 = a1.dot(a2);
    let a_v = 1.0 - a1_dot_a2 * a1_dot_a2;

    let delta = o1 - o2; // O1 - O2

    const N_SAMPLES: usize = 128;
    const CHORD_TOL: f64 = crate::bop::int_tools::CHORD_TOLERANCE;
    const REFINE_DEPTH: usize = crate::bop::int_tools::CHORD_REFINE_DEPTH;
    let mut branch_plus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);
    let mut branch_minus: Vec<(f64, DVec3)> = Vec::with_capacity(N_SAMPLES + 1);

    for i in 0..=N_SAMPLES {
        let u = (i as f64 / N_SAMPLES as f64) * TAU;
        let (cos_u, sin_u) = (u.cos(), u.sin());

        // D0(u) = (O1 - O2) + r1路(cos(u)路x1 + sin(u)路y1)
        let d0 = delta + r1 * (cos_u * x1 + sin_u * y1);
        let d0_a1 = d0.dot(a1);
        let d0_a2 = d0.dot(a2);
        let d0_sq = d0.length_squared();

        // b_v(u) = 2路(D0路a1) - 2路(D0路a2)路(a1路a2)
        let b_v = 2.0 * d0_a1 - 2.0 * d0_a2 * a1_dot_a2;

        // c_v(u) = |D0|虏 - (D0路a2)虏 - r2虏
        let c_v = d0_sq - d0_a2 * d0_a2 - r2_sq;

        if a_v.abs() > 1e-12 {
            // Quadratic: a_v路v虏 + b_v路v + c_v = 0
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                continue;
            }
            let sqrt_disc = disc.sqrt();
            let two_a_v = 2.0 * a_v;

            let v_plus = (-b_v + sqrt_disc) / two_a_v;
            let v_minus = (-b_v - sqrt_disc) / two_a_v;

            if v_plus.is_finite() {
                let p = cyl1.point_at(u, v_plus);
                if p.is_finite() {
                    branch_plus.push((u, p));
                }
            }
            if v_minus.is_finite() {
                let p = cyl1.point_at(u, v_minus);
                if p.is_finite() {
                    branch_minus.push((u, p));
                }
            }
        } else if b_v.abs() > 1e-12 {
            // a_v 鈮?0 (near-parallel axes): solve linear b_v路v + c_v = 0
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cyl1.point_at(u, v);
                if p.is_finite() {
                    branch_plus.push((u, p));
                    branch_minus.push((u, p));
                }
            }
        }
    }

    // Adaptive refinement: subdivide segments where chord error exceeds tolerance
    let eval_plus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d0 = delta + r1 * (cos_u * x1 + sin_u * y1);
        let d0_a1 = d0.dot(a1);
        let d0_a2 = d0.dot(a2);
        let d0_sq = d0.length_squared();
        let b_v = 2.0 * d0_a1 - 2.0 * d0_a2 * a1_dot_a2;
        let c_v = d0_sq - d0_a2 * d0_a2 - r2_sq;
        if a_v.abs() > 1e-12 {
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                return None;
            }
            let v = (-b_v + disc.sqrt()) / (2.0 * a_v);
            if v.is_finite() {
                let p = cyl1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else if b_v.abs() > 1e-12 {
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cyl1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else {
            None
        }
    };
    let eval_minus = |u: f64| -> Option<DVec3> {
        let (cos_u, sin_u) = (u.cos(), u.sin());
        let d0 = delta + r1 * (cos_u * x1 + sin_u * y1);
        let d0_a1 = d0.dot(a1);
        let d0_a2 = d0.dot(a2);
        let d0_sq = d0.length_squared();
        let b_v = 2.0 * d0_a1 - 2.0 * d0_a2 * a1_dot_a2;
        let c_v = d0_sq - d0_a2 * d0_a2 - r2_sq;
        if a_v.abs() > 1e-12 {
            let disc = b_v * b_v - 4.0 * a_v * c_v;
            if disc < 0.0 {
                return None;
            }
            let v = (-b_v - disc.sqrt()) / (2.0 * a_v);
            if v.is_finite() {
                let p = cyl1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else if b_v.abs() > 1e-12 {
            let v = -c_v / b_v;
            if v.is_finite() {
                let p = cyl1.point_at(u, v);
                if p.is_finite() {
                    return Some(p);
                }
            }
            None
        } else {
            None
        }
    };

    let (mut branch_plus, mut branch_minus): (Vec<DVec3>, Vec<DVec3>) = (
        refine_polyline(&branch_plus, eval_plus, CHORD_TOL, REFINE_DEPTH)
            .into_iter()
            .map(|(_, p)| p)
            .collect(),
        refine_polyline(&branch_minus, eval_minus, CHORD_TOL, REFINE_DEPTH)
            .into_iter()
            .map(|(_, p)| p)
            .collect(),
    );

    // Dedup: remove trailing points that nearly duplicate the first point
    // (closed curve degeneracy).
    let dedup = |pts: &mut Vec<DVec3>| {
        while pts.len() >= 3 {
            let n = pts.len();
            if (pts[n - 1] - pts[0]).length_squared() < 1e-24 {
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
        // Check the minus branch is distinct from the plus branch.
        let is_distinct = branches.is_empty()
            || (branch_minus[0] - branches[0][0]).length_squared() > 1e-24;
        if is_distinct {
            dedup(&mut branch_minus);
            if branch_minus.len() >= 2 {
                branches.push(branch_minus);
            }
        }
    }

    branches
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Perpendicular offset curves sampling
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€

/// Sample the perpendicular offset curves parameterization and return polylines.
///
/// For two perpendicular cylinders with offset (non-intersecting) axes,
/// the intersection curve(s) can be parameterized on cyl1's surface:
///
/// ```text
/// P(胃) = O1 + v(胃)路a1 + R1路(cos(胃)路u1 + sin(胃)路v1)
/// v(胃) = dz 卤 鈭?R2虏 - (R1路cos(胃) - dx)虏)
/// ```
///
/// Returns two polylines (one per closed loop), each combining the + and -
/// branches. The loops meet at tangent points where the discriminant is zero.
pub fn sample_perpendicular_offset_curves(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    _dist: f64,
    n_samples: usize,
) -> Vec<Vec<DVec3>> {
    let a1 = cyl1.axis.normalize();
    let a2 = cyl2.axis.normalize();
    let r1 = cyl1.radius;
    let r2 = cyl2.radius;
    let r2_sq = r2 * r2;
    let w = cyl1.origin - cyl2.origin;
    let denom = 1.0 - a1.dot(a2) * a1.dot(a2);
    if denom.abs() < 1e-12 {
        return vec![];
    }
    let d1 = a1.dot(w);
    let d2 = a2.dot(w);
    let t = (a1.dot(a2) * d2 - d1) / denom;
    let s = (d2 - a1.dot(a2) * d1) / denom;
    let conn = (cyl1.origin + a1 * t) - (cyl2.origin + a2 * s);
    let conn_len = conn.length();
    let u1 = if conn_len < 1e-12 {
        a1.cross(a2).normalize()
    } else {
        conn / conn_len
    };
    let v1 = a1.cross(u1).normalize();
    let delta = cyl2.origin - cyl1.origin;
    let dx = delta.dot(u1);
    let dz = delta.dot(a1);
    let cos_min = ((dx - r2) / r1).clamp(-1.0, 1.0);
    let cos_max = ((dx + r2) / r1).clamp(-1.0, 1.0);
    if cos_min > cos_max {
        return vec![];
    }
    let t_low = cos_max.acos();
    let t_high = cos_min.acos();

    let mut branches = Vec::new();

    for (t_start, t_end) in [(t_low, t_high), (TAU - t_high, TAU - t_low)] {
        // Forward: branch = +1, 胃 = t_start 鈫?t_end
        // Backward: branch = -1, 胃 = t_end 鈫?t_start (reversed in the loop)
        // Combined they form a single closed curve.
        let n_pts = n_samples * 2 + 1;
        let mut pts: Vec<DVec3> = Vec::with_capacity(n_pts);

        // Forward: positive sqrt branch
        for i in 0..=n_samples {
            let theta = t_start + (t_end - t_start) * i as f64 / n_samples as f64;
            let (ct, st) = (theta.cos(), theta.sin());
            let diff = r1 * ct - dx;
            let disc = (r2_sq - diff * diff).max(0.0).sqrt();
            let v_z = dz + disc;
            pts.push(cyl1.origin + v_z * a1 + r1 * (ct * u1 + st * v1));
        }
        // Backward: negative sqrt branch (reversed, skip the already-sampled t_end)
        for i in 1..=n_samples {
            let theta = t_end - (t_end - t_start) * i as f64 / n_samples as f64;
            let (ct, st) = (theta.cos(), theta.sin());
            let diff = r1 * ct - dx;
            let disc = (r2_sq - diff * diff).max(0.0).sqrt();
            let v_z = dz - disc;
            pts.push(cyl1.origin + v_z * a1 + r1 * (ct * u1 + st * v1));
        }

        if pts.len() >= 2 {
            branches.push(pts);
        }
    }

    branches
}

// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
// Tests
// 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
