use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

#[derive(Debug, Clone)]
pub enum PlaneConicalResult {
    NoIntersection,
    Point(DVec3),
    SingleLine(Line3),
    TwoLines(Line3, Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
    // Parabola/Hyperbola would need new Curve3 variants — deferred to Phase 6.
}

pub fn intersect_plane_cone(plane: &Plane, cone: &ConicalSurface) -> PlaneConicalResult {
    let cos_angle = plane.normal.dot(cone.axis).abs();
    let sin_angle = (1.0 - cos_angle * cos_angle).sqrt().max(0.0);

    // Distance from apex to plane
    let apex_to_plane = (plane.origin - cone.apex).dot(plane.normal);

    if cos_angle < TOLERANCE_ANG {
        // Plane parallel to axis
        if apex_to_plane.abs() < TOLERANCE_ABS {
            // Plane passes through apex
            return PlaneConicalResult::Point(cone.apex);
        }
        // May intersect cone in hyperbola — approximate as no intersection for now
        return PlaneConicalResult::NoIntersection;
    }

    if (cos_angle - 1.0).abs() < TOLERANCE_ANG {
        // Plane perpendicular to axis → circle
        if apex_to_plane.abs() < TOLERANCE_ABS {
            return PlaneConicalResult::Point(cone.apex);
        }
        let t = apex_to_plane / cone.axis.dot(plane.normal);
        let center = cone.apex + cone.axis * t;
        let radius = (t * cone.half_angle_rad.tan()).abs();
        if radius < TOLERANCE_ABS {
            return PlaneConicalResult::Point(center);
        }
        return PlaneConicalResult::Circle(Circle3 {
            center,
            normal: cone.axis,
            radius,
        });
    }

    // Check if plane passes through apex
    if apex_to_plane.abs() < TOLERANCE_ABS {
        // Plane through apex — result is point, one line, or two lines
        let angle_between = sin_angle.atan2(cos_angle);
        let half_angle = cone.half_angle_rad;

        if (angle_between - half_angle).abs() < TOLERANCE_ANG {
            // Tangent — single line
            let dir = plane.normal.cross(cone.axis).normalize();
            let gen_dir = (cone.axis * half_angle.cos() + dir * half_angle.sin()).normalize();
            return PlaneConicalResult::SingleLine(Line3 {
                origin: cone.apex,
                direction: gen_dir,
            });
        }

        if angle_between < half_angle {
            // Two generatrices
            let cross = plane.normal.cross(cone.axis);
            if is_zero_vec(cross) {
                return PlaneConicalResult::Point(cone.apex);
            }
            let perp_in_plane = cross.normalize();
            let projected_axis =
                (cone.axis - plane.normal * cone.axis.dot(plane.normal)).normalize();

            let d1 =
                (projected_axis * half_angle.cos() + perp_in_plane * half_angle.sin()).normalize();
            let d2 =
                (projected_axis * half_angle.cos() - perp_in_plane * half_angle.sin()).normalize();

            return PlaneConicalResult::TwoLines(
                Line3 {
                    origin: cone.apex,
                    direction: d1,
                },
                Line3 {
                    origin: cone.apex,
                    direction: d2,
                },
            );
        }

        return PlaneConicalResult::Point(cone.apex);
    }

    // General oblique case — ellipse (when angle between plane and axis < cone half-angle)
    // This is an approximation; full conic section analysis would produce parabolas/hyperbolas too.
    let t = apex_to_plane / cone.axis.dot(plane.normal);
    let center = cone.apex + cone.axis * t;
    let base_radius = (t * cone.half_angle_rad.tan()).abs();

    if base_radius < TOLERANCE_ABS {
        return PlaneConicalResult::Point(center);
    }

    let major_radius = base_radius / cos_angle;
    let minor_radius = base_radius;
    let major_dir = (cone.axis - plane.normal * cone.axis.dot(plane.normal)).normalize();

    PlaneConicalResult::Ellipse(Ellipse3 {
        center,
        normal: plane.normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_cone() -> ConicalSurface {
        ConicalSurface {
            apex: DVec3::ZERO,
            axis: DVec3::Y,
            radius: 1.0,
            half_angle_rad: std::f64::consts::FRAC_PI_4, // 45°
        }
    }

    #[test]
    fn perpendicular_plane_circle() {
        let plane = Plane {
            origin: DVec3::new(0.0, 2.0, 0.0),
            normal: DVec3::Y,
        };
        match intersect_plane_cone(&plane, &test_cone()) {
            PlaneConicalResult::Circle(c) => {
                assert!((c.center.y - 2.0).abs() < TOLERANCE_ABS);
                assert!((c.radius - 2.0).abs() < 0.01); // tan(45°)*2 = 2
            }
            other => panic!("Expected Circle, got {other:?}"),
        }
    }

    #[test]
    fn plane_through_apex() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::X,
        };
        let result = intersect_plane_cone(&plane, &test_cone());
        // Plane through apex, perpendicular to X → two lines or point
        assert!(matches!(
            result,
            PlaneConicalResult::TwoLines(_, _) | PlaneConicalResult::Point(_)
        ));
    }
}
