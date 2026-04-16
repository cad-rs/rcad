//! Analytic intersection of a plane with a torus.
//!
//! # Cases
//!
//! - **Perpendicular to axis**: Two circles (inner and outer equator)
//! - **Parallel to axis**: Two circles or figure-8 depending on offset
//! - **Oblique**: Complex curve, fall back to numerical marching

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Plane, ToroidalSurface};

use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_ANG};

/// Result of plane x torus intersection.
#[derive(Debug, Clone)]
pub enum PlaneTorusResult {
    /// The plane does not intersect the torus.
    NoIntersection,
    /// Single tangent circle.
    TangentCircle(Circle3),
    /// Two circles (perpendicular case).
    TwoCircles(Circle3, Circle3),
    /// Complex intersection, fall back to numerical marching.
    General,
}

/// Compute the analytic intersection of `plane` and `torus`.
pub fn intersect_plane_torus(plane: &Plane, torus: &ToroidalSurface) -> PlaneTorusResult {
    intersect_plane_torus_with_tolerance(plane, torus, 0.0)
}

/// Plane x torus intersection with fuzzy tolerance.
pub fn intersect_plane_torus_with_tolerance(
    plane: &Plane,
    torus: &ToroidalSurface,
    fuzzy_tol: f64,
) -> PlaneTorusResult {
    let tol = TOLERANCE_ABS + fuzzy_tol.max(0.0);

    // Normalize plane normal and torus axis
    let n = plane.normal.normalize();
    let a = torus.axis.normalize();

    // Check if plane is perpendicular to torus axis
    let dot_na = n.dot(a).abs();

    if dot_na > 1.0 - TOLERANCE_ANG {
        // Plane perpendicular to axis: circular cross-section
        return intersect_plane_torus_perpendicular(plane, torus, tol);
    }

    // Check if plane is parallel to torus axis
    if dot_na < TOLERANCE_ANG {
        // Plane parallel to axis: may produce two circles
        return intersect_plane_torus_parallel(plane, torus, tol);
    }

    // General oblique case: fall back to numerical
    PlaneTorusResult::General
}

fn intersect_plane_torus_perpendicular(
    plane: &Plane,
    torus: &ToroidalSurface,
    tol: f64,
) -> PlaneTorusResult {
    // Distance from torus center to plane along axis
    let signed_dist = (torus.center - plane.origin).dot(torus.axis);
    let abs_dist = signed_dist.abs();

    // Maximum distance for intersection is the minor radius
    if abs_dist > torus.minor_radius + tol {
        return PlaneTorusResult::NoIntersection;
    }

    // Tangent case: one circle
    if (abs_dist - torus.minor_radius).abs() < tol {
        let center = torus.center - torus.axis * signed_dist;
        return PlaneTorusResult::TangentCircle(Circle3 {
            center,
            normal: torus.axis,
            radius: torus.major_radius,
        });
    }

    // Two circles at height signed_dist from torus center
    // Circle radius on the tube: sqrt(r^2 - d^2) where r = minor_radius, d = distance
    let tube_circle_r = (torus.minor_radius * torus.minor_radius - signed_dist * signed_dist).sqrt();

    // Two circles at major_radius +/- tube_circle_r from axis
    let r1 = torus.major_radius + tube_circle_r;
    let r2 = (torus.major_radius - tube_circle_r).max(0.0);

    let center = torus.center - torus.axis * signed_dist;

    if r2 < tol {
        // Inner circle degenerates to point
        PlaneTorusResult::TangentCircle(Circle3 {
            center,
            normal: torus.axis,
            radius: r1,
        })
    } else {
        PlaneTorusResult::TwoCircles(
            Circle3 { center, normal: torus.axis, radius: r1 },
            Circle3 { center, normal: torus.axis, radius: r2 },
        )
    }
}

fn intersect_plane_torus_parallel(
    plane: &Plane,
    torus: &ToroidalSurface,
    tol: f64,
) -> PlaneTorusResult {
    // Distance from torus axis to plane
    let to_plane = plane.origin - torus.center;
    let along_axis = to_plane.dot(torus.axis);
    let in_plane = to_plane - torus.axis * along_axis;
    let dist_to_axis = in_plane.length();

    // Check if plane passes through the torus tube
    // The torus tube sweeps from (R-r) to (R+r) from axis
    let r_min = torus.major_radius - torus.minor_radius;
    let r_max = torus.major_radius + torus.minor_radius;

    if dist_to_axis > r_max + tol {
        return PlaneTorusResult::NoIntersection;
    }

    // For now, return General for parallel case (complex geometry)
    // TODO: Implement analytic circles for specific parallel configurations
    PlaneTorusResult::General
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plane_perpendicular_to_torus_axis_produces_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane perpendicular to Y axis, slicing through center
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // Outer circle at major_radius + minor_radius
                assert!((c1.radius - 6.0).abs() < 1e-6, "Outer circle radius expected 6.0, got {}", c1.radius);
                // Inner circle at major_radius - minor_radius
                assert!((c2.radius - 4.0).abs() < 1e-6, "Inner circle radius expected 4.0, got {}", c2.radius);
                // Both circles should have the same center
                assert!((c1.center - DVec3::ZERO).length() < 1e-6);
                assert!((c2.center - DVec3::ZERO).length() < 1e-6);
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_parallel_to_torus_axis_returns_general() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::General));
    }

    #[test]
    fn plane_perpendicular_tangent_to_torus() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane tangent to torus (at top of tube)
        let plane = Plane {
            origin: DVec3::new(0.0, 1.0, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TangentCircle(c) => {
                // Tangent circle at the major radius
                assert!((c.radius - 5.0).abs() < 1e-6);
            }
            other => panic!("Expected TangentCircle, got {:?}", other),
        }
    }

    #[test]
    fn plane_perpendicular_no_intersection() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane outside torus
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::NoIntersection));
    }

    #[test]
    fn plane_perpendicular_offset_produces_two_circles() {
        // Torus centered at origin with axis along Y
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane offset by 0.5 from center
        let plane = Plane {
            origin: DVec3::new(0.0, 0.5, 0.0),
            normal: DVec3::Y,
        };

        let result = intersect_plane_torus(&plane, &torus);

        match result {
            PlaneTorusResult::TwoCircles(c1, c2) => {
                // tube_circle_r = sqrt(1 - 0.25) = sqrt(0.75) = 0.866025...
                let expected_tube_r = (1.0_f64 * 1.0 - 0.5_f64 * 0.5).sqrt();
                let expected_r1 = 5.0 + expected_tube_r;
                let expected_r2 = 5.0 - expected_tube_r;

                assert!((c1.radius - expected_r1).abs() < 1e-6, "Outer circle radius mismatch");
                assert!((c2.radius - expected_r2).abs() < 1e-6, "Inner circle radius mismatch");
            }
            other => panic!("Expected TwoCircles, got {:?}", other),
        }
    }

    #[test]
    fn plane_oblique_returns_general() {
        let torus = ToroidalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Y,
            major_radius: 5.0,
            minor_radius: 1.0,
        };
        // Plane at 45 degrees to torus axis
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::new(1.0, 1.0, 0.0).normalize(),
        };

        let result = intersect_plane_torus(&plane, &torus);
        assert!(matches!(result, PlaneTorusResult::General));
    }
}
