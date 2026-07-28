#[allow(unused_imports)]
use glam::DVec3;
use rcad_kernel::geom::*;



#[derive(Debug, Clone)]
pub enum PlaneCylinderResult {
    NoIntersection,
    TangentLine(Line3),
    TwoLines(Line3, Line3),
    Circle(Circle3),
    Ellipse(Ellipse3),
}

pub fn intersect_plane_cylinder(plane: &Plane, cyl: &CylindricalSurface) -> PlaneCylinderResult {
    let cos_angle = plane.normal.dot(cyl.axis).abs();

    if cos_angle < rcad_kernel::ANGULAR {
        // Plane parallel to cylinder axis
        let axis_to_plane = (plane.origin - cyl.origin).dot(plane.normal);
        let dist = axis_to_plane.abs();

        if dist > cyl.radius + rcad_kernel::rcad_kernel::CONFUSION {
            return PlaneCylinderResult::NoIntersection;
        }
        if (dist - cyl.radius).abs() < rcad_kernel::rcad_kernel::CONFUSION {
            let tang_point = cyl.origin + plane.normal * axis_to_plane;
            return PlaneCylinderResult::TangentLine(Line3 {
                origin: tang_point,
                direction: cyl.axis,
            });
        }
        let offset_dir = plane.normal.cross(cyl.axis).normalize();
        let half_chord = (cyl.radius * cyl.radius - dist * dist).sqrt();
        let center_on_plane = cyl.origin + plane.normal * axis_to_plane;

        let l1_origin = center_on_plane + offset_dir * half_chord;
        let l2_origin = center_on_plane - offset_dir * half_chord;

        return PlaneCylinderResult::TwoLines(
            Line3 {
                origin: l1_origin,
                direction: cyl.axis,
            },
            Line3 {
                origin: l2_origin,
                direction: cyl.axis,
            },
        );
    }

    if (cos_angle - 1.0).abs() < rcad_kernel::ANGULAR {
        // Plane perpendicular to cylinder axis 鈫?circle
        let t = (plane.origin - cyl.origin).dot(cyl.axis);
        let center = cyl.origin + cyl.axis * t;
        return PlaneCylinderResult::Circle(Circle3::new(center, cyl.axis, cyl.radius));
    }

    // General oblique case 鈫?ellipse
    let major_radius = cyl.radius / cos_angle;
    let minor_radius = cyl.radius;

    let t = (plane.origin - cyl.origin).dot(plane.normal) / cyl.axis.dot(plane.normal);
    let center = cyl.origin + cyl.axis * t;

    let major_dir = (cyl.axis - plane.normal * cyl.axis.dot(plane.normal)).normalize();

    PlaneCylinderResult::Ellipse(Ellipse3 {
        center,
        normal: plane.normal,
        major_dir,
        major_radius,
        minor_radius,
    })
}
