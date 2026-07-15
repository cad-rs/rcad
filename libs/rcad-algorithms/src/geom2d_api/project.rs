//! 3D curve → 2D curve projection onto a plane surface.
//! GeomAPI::To2d (geomapi.cxx L311-340).

use rcad_kernel::geom::*;
use glam::{DVec2, DVec3};
use crate::tolerance::{TOLERANCE_ABS, TOLERANCE_LEN_SQ_DIV_SAFE};

/// Build a 2D coordinate frame from a plane surface.
/// Returns (origin, u_axis, v_axis) using the plane's stored axes.
pub(crate) fn plane_frame(plane: &Plane) -> (DVec3, DVec3, DVec3) {
    (plane.origin, plane.u_dir, plane.v_dir)
}

/// Project a 3D point onto a plane's 2D coordinate system.
pub(crate) fn project_point_to_plane_uv(point: DVec3, plane: &Plane) -> DVec2 {
    let (origin, u_axis, v_axis) = plane_frame(plane);
    let d = point - origin;
    DVec2::new(d.dot(u_axis), d.dot(v_axis))
}

/// GeomAPI::To2d — project a 3D curve onto a plane surface.
///
/// Supports Line3, Circle3, Ellipse3. Returns None for unsupported curve types.
/// The plane's coordinate frame defines the 2D space:
///   u = (P - origin)·u_axis, v = (P - origin)·v_axis
pub fn project_curve_to_plane(curve: &Curve3, surface: &Surface3) -> Option<Curve2d> {
    let plane = match surface {
        Surface3::Plane(p) => p,
        _ => return None, // Only planar projection is supported
    };
    let (origin, u_axis, v_axis) = plane_frame(plane);

    match curve {
        Curve3::Line(l) => {
            let p_start = project_point_to_plane_uv(l.origin, plane);
            let p_end = project_point_to_plane_uv(l.origin + l.direction, plane);
            let dir = (p_end - p_start).normalize_or_zero();
            if dir.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                return None;
            }
            Some(Curve2d::Line(Line2d { origin: p_start, direction: dir }))
        }
        Curve3::Circle(c) => {
            let center_2d = project_point_to_plane_uv(c.center, plane);
            // Project x_dir axis endpoint to find 2D axes
            let p_u = project_point_to_plane_uv(c.center + c.x_dir * c.radius, plane);
            let dir_u = (p_u - center_2d).normalize_or_zero();
            if dir_u.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                return None;
            }
            let perp = DVec2::new(-dir_u.y, dir_u.x);
            let uv_radius = (p_u - center_2d).length();
            Some(Curve2d::Circle(Circle2d {
                center: center_2d,
                x_dir: dir_u,
                y_dir: perp,
                radius: uv_radius,
            }))
        }
        Curve3::Ellipse(e) => {
            let center_2d = project_point_to_plane_uv(e.center, plane);
            // Project major_dir axis endpoint to get 2D major direction and radius
            let p_major = project_point_to_plane_uv(e.center + e.major_dir * e.major_radius, plane);
            let dir_major = (p_major - center_2d).normalize_or_zero();
            if dir_major.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                return None;
            }
            let uv_major_radius = (p_major - center_2d).length();
            // Project minor axis endpoint
            let minor_dir_3d = if e.major_dir.cross(e.normal).length_squared() > TOLERANCE_LEN_SQ_DIV_SAFE {
                e.major_dir.cross(e.normal).normalize()
            } else {
                any_perpendicular(e.major_dir).normalize()
            };
            let p_minor = project_point_to_plane_uv(e.center + minor_dir_3d * e.minor_radius, plane);
            let dir_minor = (p_minor - center_2d).normalize_or_zero();
            let uv_minor_radius = (p_minor - center_2d).length();

            // Ensure axes are perpendicular in 2D
            let y_dir = DVec2::new(-dir_major.y, dir_major.x);
            let minor_comp = dir_minor.dot(y_dir);
            let corrected_minor = if minor_comp.abs() > 0.5 {
                uv_minor_radius
            } else {
                uv_minor_radius.max(uv_major_radius * 0.01) // degenerate ellipse guard
            };

            Some(Curve2d::Ellipse(Ellipse2d {
                center: center_2d,
                major_dir: dir_major,
                major_radius: uv_major_radius.max(TOLERANCE_ABS),
                minor_radius: corrected_minor.max(TOLERANCE_ABS),
            }))
        }
        _ => None,
    }
}
