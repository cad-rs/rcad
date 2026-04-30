use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

/// V-V: test if two vertices coincide.
pub fn vertex_vertex_coincide(p1: DVec3, p2: DVec3) -> bool {
    points_coincide(p1, p2)
}

/// V-E: test if a vertex lies on a line edge within its parametric range.
/// Returns the parameter if it does.
pub fn vertex_on_line(point: DVec3, line: &Line3, t_range: [f64; 2]) -> Option<f64> {
    vertex_on_line_with_tol(point, line, t_range, TOLERANCE_ABS)
}

/// Same as [`vertex_on_line`], with a caller-supplied coincidence / param-margin tolerance
/// (clamped to at least [`TOLERANCE_ABS`]). Used by [`crate::pave_filler::PaveFiller`] with
/// [`bopds::ds::DS::fuzzy_tol`].
pub fn vertex_on_line_with_tol(
    point: DVec3,
    line: &Line3,
    t_range: [f64; 2],
    coincident_tol: f64,
) -> Option<f64> {
    let tol = coincident_tol.max(TOLERANCE_ABS);
    let tol_sq = tol * tol;
    let v = point - line.origin;
    let t = v.dot(line.direction);
    let closest = line.origin + line.direction * t;
    if (closest - point).length_squared() > tol_sq {
        return None;
    }
    if t < t_range[0] - tol || t > t_range[1] + tol {
        return None;
    }
    Some(t)
}

/// V-F: test if a vertex lies on a plane (within tolerance).
pub fn vertex_on_plane(point: DVec3, plane: &Plane) -> bool {
    vertex_on_plane_with_tol(point, plane, TOLERANCE_ABS)
}

/// Same as [`vertex_on_plane`] with explicit distance tolerance (minimum [`TOLERANCE_ABS`]).
pub fn vertex_on_plane_with_tol(point: DVec3, plane: &Plane, distance_tol: f64) -> bool {
    let tol = distance_tol.max(TOLERANCE_ABS);
    let d = (point - plane.origin).dot(plane.normal);
    d.abs() < tol
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_on_line_segment() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let t = vertex_on_line(DVec3::new(0.5, 0.0, 0.0), &line, [0.0, 1.0]);
        assert!(t.is_some());
        assert!((t.unwrap() - 0.5).abs() < TOLERANCE_ABS);
    }

    #[test]
    fn vertex_off_line() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        assert!(vertex_on_line(DVec3::new(0.5, 1.0, 0.0), &line, [0.0, 1.0]).is_none());
    }

    #[test]
    fn vertex_on_plane_test() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        assert!(vertex_on_plane(DVec3::new(5.0, 3.0, 0.0), &plane));
        assert!(!vertex_on_plane(DVec3::new(5.0, 3.0, 1.0), &plane));
    }

    #[test]
    fn vertex_on_line_with_tol_allows_near_miss_inside_margin() {
        let line = Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        };
        let pt = DVec3::new(0.5, 2e-7, 0.0);
        assert!(vertex_on_line(pt, &line, [0.0, 1.0]).is_none());
        assert!(vertex_on_line_with_tol(pt, &line, [0.0, 1.0], TOLERANCE_MESH_LEGACY).is_some());
    }
}
