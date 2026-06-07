use glam::DVec3;
use rcad_kernel::geom::*;

use crate::tolerance::*;

pub struct EdgeFaceHit {
    pub point: DVec3,
    pub edge_param: f64,
}

/// Intersect a line segment (bounded by t_range) with a plane.
/// Does NOT check face boundary containment — caller must do that.
pub fn intersect_line_plane(line: &Line3, t_range: [f64; 2], plane: &Plane) -> Option<EdgeFaceHit> {
    intersect_line_plane_with_tol(line, t_range, plane, TOLERANCE_ABS)
}

/// Same as [`intersect_line_plane`] with explicit edge-parameter margin (minimum [`TOLERANCE_ABS`]).
/// Parallel/near-parallel denom threshold stays strict at [`TOLERANCE_ABS`].
pub fn intersect_line_plane_with_tol(
    line: &Line3,
    t_range: [f64; 2],
    plane: &Plane,
    param_tol: f64,
) -> Option<EdgeFaceHit> {
    let ptol = param_tol.max(TOLERANCE_ABS);
    let denom = line.direction.dot(plane.normal);
    if denom.abs() < TOLERANCE_ABS {
        return None;
    }
    let t = (plane.origin - line.origin).dot(plane.normal) / denom;
    if t < t_range[0] - ptol || t > t_range[1] + ptol {
        return None;
    }
    let point = line.origin + line.direction * t;
    Some(EdgeFaceHit {
        point,
        edge_param: t,
    })
}

/// Build a local 2D basis on the plane (u, v axes).
pub fn plane_local_basis(plane: &Plane) -> (DVec3, DVec3) {
    let n = plane.normal;
    let ref_dir = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let u = n.cross(ref_dir).normalize();
    let v = n.cross(u).normalize();
    (u, v)
}

/// Check if `point` lies inside a planar face whose boundary vertices are given
/// in order. Uses 2D projection + ray-casting.
pub fn point_in_planar_face(point: DVec3, plane: &Plane, face_verts: &[DVec3]) -> bool {
    point_in_planar_face_with_tol(point, plane, face_verts, TOLERANCE_ABS)
}

/// Same as [`point_in_planar_face`], with a 2D ray-cast margin (minimum [`TOLERANCE_ABS`]).
/// Use the same magnitude as pave [`bopds::ds::DS::fuzzy_tol`] for consistent V–F containment.
pub fn point_in_planar_face_with_tol(
    point: DVec3,
    plane: &Plane,
    face_verts: &[DVec3],
    geom_tol: f64,
) -> bool {
    if face_verts.len() < 3 {
        return false;
    }
    let eps = geom_tol.max(TOLERANCE_ABS);
    let (u_axis, v_axis) = plane_local_basis(plane);

    let project = |p: DVec3| -> (f64, f64) {
        let d = p - plane.origin;
        (d.dot(u_axis), d.dot(v_axis))
    };

    let (px, py) = project(point);
    let poly: Vec<(f64, f64)> = face_verts.iter().map(|v| project(*v)).collect();

    ray_cast_contains_with_tol(px, py, &poly, eps)
}

fn ray_cast_contains_with_tol(px: f64, py: f64, poly: &[(f64, f64)], eps: f64) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = poly[i];
        let (xj, yj) = poly[j];
        if (yi >= py - eps) != (yj >= py - eps) {
            let dy = yj - yi;
            if dy.abs() < eps {
                j = i;
                continue;
            }
            let xint = (xj - xi) * (py - yi) / dy + xi;
            if px < xint + eps {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Clip an infinite line to a convex polygon on a plane.
/// Returns the parametric interval `(t_min, t_max)` of the line inside the polygon,
/// or None if the line doesn't cross the polygon.
///
/// This is a wrapper around [`clip_line_to_polygon_with_tol`] that returns only the
/// first inside interval for backward compatibility with convex polygons.
pub fn clip_line_to_convex_polygon(
    line: &Line3,
    plane: &Plane,
    face_verts: &[DVec3],
) -> Option<(f64, f64)> {
    let intervals = clip_line_to_polygon_with_tol(line, plane, face_verts, TOLERANCE_ABS);
    intervals.first().copied()
}

/// Clip an infinite line to a (possibly non-convex) polygon on a plane.
/// Returns all parametric intervals along the line that lie inside the polygon,
/// or an empty vec if the line doesn't cross the polygon.
///
/// Algorithm: find all edge-intersection t-values along the line, sort and dedup,
/// then test each interval midpoint via ray-cast point-in-polygon.
pub fn clip_line_to_polygon_with_tol(
    line: &Line3,
    plane: &Plane,
    face_verts: &[DVec3],
    geom_tol: f64,
) -> Vec<(f64, f64)> {
    let eps = geom_tol.max(TOLERANCE_ABS);
    if face_verts.len() < 3 {
        return vec![];
    }
    let (u_axis, v_axis) = plane_local_basis(plane);

    // Project line direction and origin onto 2D
    let line_u = line.direction.dot(u_axis);
    let line_v = line.direction.dot(v_axis);
    let d = line.origin - plane.origin;
    let origin_u = d.dot(u_axis);
    let origin_v = d.dot(v_axis);

    // Project polygon vertices to 2D
    let pts_2d: Vec<(f64, f64)> = face_verts
        .iter()
        .map(|v| {
            let d = *v - plane.origin;
            (d.dot(u_axis), d.dot(v_axis))
        })
        .collect();

    let n = pts_2d.len();
    let mut t_vals = Vec::new();

    // Find all edge intersection t-values: for each polygon edge, compute
    // where the line crosses the edge's supporting line and check if the
    // intersection lies within the edge segment.
    for i in 0..n {
        let j = (i + 1) % n;
        let (ax, ay) = pts_2d[i];
        let (bx, by) = pts_2d[j];
        let ex = bx - ax;
        let ey = by - ay;

        // line_dir × edge_dir — zero means parallel
        let denom = line_u * ey - line_v * ex;
        if denom.abs() < eps {
            // OCCT-aligned: when the line coincides with a polygon edge (parallel
            // AND zero distance), use the edge endpoints as intersection t-values.
            // Without this, boundary-coincident intersection lines lose the portion
            // of the line that runs along the polygon edge (e.g. bfuse_simple B3
            // face[6] intersection at z=0, which is on the face boundary).
            let dist = (origin_u - ax) * ey - (origin_v - ay) * ex;
            if dist.abs() < eps {
                // Line coincides with this edge — add t at both endpoints
                let dir_len2 = line_u * line_u + line_v * line_v;
                if dir_len2 > 1e-30 {
                    let t_a = ((ax - origin_u) * line_u + (ay - origin_v) * line_v) / dir_len2;
                    let t_b = ((bx - origin_u) * line_u + (by - origin_v) * line_v) / dir_len2;
                    t_vals.push(t_a);
                    t_vals.push(t_b);
                }
            }
            continue;
        }

        // t where line crosses the edge's supporting line:
        //   line(t) = edge(i) + s * (edge(i+1) - edge(i))
        //   (p[i] - origin) × edge_dir / (line_dir × edge_dir)
        let t = ((ax - origin_u) * ey - (ay - origin_v) * ex) / denom;

        // Edge parameter s — check the intersection is within the segment
        let s = if ex.abs() > eps {
            (origin_u + t * line_u - ax) / ex
        } else {
            (origin_v + t * line_v - ay) / ey
        };

        if s >= -eps && s <= 1.0 + eps {
            t_vals.push(t);
        }
    }

    if t_vals.is_empty() {
        // Line doesn't cross any edge — check if origin is inside
        if point_in_planar_face_with_tol(line.origin, plane, face_verts, geom_tol) {
            return vec![(f64::NEG_INFINITY, f64::INFINITY)];
        }
        return vec![];
    }

    // Sort by t-value
    t_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Remove duplicates within tolerance
    let mut deduped: Vec<f64> = Vec::new();
    for &t in &t_vals {
        if deduped.is_empty() {
            deduped.push(t);
        } else {
            let last = deduped[deduped.len() - 1];
            if (t - last).abs() > eps {
                deduped.push(t);
            }
        }
    }

    // Compute polygon's t-range along the line direction to choose a safe
    // "far" offset for testing unbounded intervals.
    let line_dir_len_sq = line_u * line_u + line_v * line_v;
    let (min_vert_t, max_vert_t) = if line_dir_len_sq > eps {
        let mut min_t = f64::INFINITY;
        let mut max_t = f64::NEG_INFINITY;
        for &(u, v) in &pts_2d {
            let t = ((u - origin_u) * line_u + (v - origin_v) * line_v) / line_dir_len_sq;
            min_t = min_t.min(t);
            max_t = max_t.max(t);
        }
        (min_t, max_t)
    } else {
        (0.0, 0.0)
    };
    let far_offset = (max_vert_t - min_vert_t + 1.0).max(1.0) * 2.0;

    let mut result = Vec::new();

    // Unbounded interval before first t-value
    {
        let test_t = deduped[0] - far_offset;
        let test_pt = line.origin + line.direction * test_t;
        if point_in_planar_face_with_tol(test_pt, plane, face_verts, geom_tol) {
            result.push((f64::NEG_INFINITY, deduped[0]));
        }
    }

    // Intervals between consecutive t-values
    for k in 0..deduped.len() - 1 {
        let t_mid = (deduped[k] + deduped[k + 1]) / 2.0;
        let mid_pt = line.origin + line.direction * t_mid;
        if point_in_planar_face_with_tol(mid_pt, plane, face_verts, geom_tol) {
            result.push((deduped[k], deduped[k + 1]));
        } else {
            // The midpoint might land exactly on a boundary edge (coincident
            // line case).  Nudge perpendicular to the line direction (try
            // BOTH directions) and retry.
            let perp = line.direction.cross(plane.normal).normalize_or_zero() * (eps * 10.0);
            if perp.length_squared() > 0.0 {
                if point_in_planar_face_with_tol(mid_pt + perp, plane, face_verts, geom_tol)
                    || point_in_planar_face_with_tol(mid_pt - perp, plane, face_verts, geom_tol)
                {
                    result.push((deduped[k], deduped[k + 1]));
                }
            }
        }
    }

    // Unbounded interval after last t-value
    {
        let test_t = deduped[deduped.len() - 1] + far_offset;
        let test_pt = line.origin + line.direction * test_t;
        if point_in_planar_face_with_tol(test_pt, plane, face_verts, geom_tol) {
            result.push((deduped[deduped.len() - 1], f64::INFINITY));
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_through_plane() {
        let line = Line3 {
            origin: DVec3::new(0.5, 0.5, -1.0),
            direction: DVec3::Z,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let hit = intersect_line_plane(&line, [-10.0, 10.0], &plane).unwrap();
        assert!((hit.edge_param - 1.0).abs() < TOLERANCE_ABS);
        assert!(points_coincide(hit.point, DVec3::new(0.5, 0.5, 0.0)));
    }

    #[test]
    fn line_parallel_to_plane() {
        let line = Line3 {
            origin: DVec3::new(0.0, 0.0, 1.0),
            direction: DVec3::X,
        };
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        assert!(intersect_line_plane(&line, [-10.0, 10.0], &plane).is_none());
    }

    #[test]
    fn point_inside_square() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        assert!(point_in_planar_face(
            DVec3::new(0.5, 0.5, 0.0),
            &plane,
            &verts
        ));
        assert!(!point_in_planar_face(
            DVec3::new(1.5, 0.5, 0.0),
            &plane,
            &verts
        ));
    }

    #[test]
    fn clip_line_to_square() {
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let line = Line3 {
            origin: DVec3::new(0.5, -1.0, 0.0),
            direction: DVec3::Y,
        };
        let (t_min, t_max) = clip_line_to_convex_polygon(&line, &plane, &verts).unwrap();
        assert!((t_min - 1.0).abs() < TOLERANCE_MESH_LEGACY, "t_min={t_min}");
        assert!((t_max - 2.0).abs() < TOLERANCE_MESH_LEGACY, "t_max={t_max}");

        // New API also returns the same single interval for convex polygons
        let intervals = clip_line_to_polygon_with_tol(&line, &plane, &verts, TOLERANCE_ABS);
        assert_eq!(intervals.len(), 1);
        assert!((intervals[0].0 - 1.0).abs() < TOLERANCE_MESH_LEGACY);
        assert!((intervals[0].1 - 2.0).abs() < TOLERANCE_MESH_LEGACY);
    }

    #[test]
    fn clip_line_to_l_shape() {
        // L-shaped polygon matching the H1/H2 top cap profile:
        //   (0,0)-(2,0)-(2,1)-(1,1)-(1,3)-(0,3)
        // The reflex vertex at (1,1) makes this non-convex.
        let plane = Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
        };
        let verts = vec![
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(2.0, 0.0, 0.0),
            DVec3::new(2.0, 1.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(1.0, 3.0, 0.0),
            DVec3::new(0.0, 3.0, 0.0),
        ];

        // Line through the wide region (y = 0.5): should clip from x=0 to x=2
        let line1 = Line3 {
            origin: DVec3::new(-1.0, 0.5, 0.0),
            direction: DVec3::X,
        };
        let intervals1 = clip_line_to_polygon_with_tol(&line1, &plane, &verts, TOLERANCE_ABS);
        assert_eq!(intervals1.len(), 1, "line through wide region should have 1 interval");
        assert!((intervals1[0].0 - 1.0).abs() < TOLERANCE_MESH_LEGACY, "t_min={}", intervals1[0].0);
        assert!((intervals1[0].1 - 3.0).abs() < TOLERANCE_MESH_LEGACY, "t_max={}", intervals1[0].1);

        // Line through the narrow arm (y = 1.5): should clip from x=0 to x=1
        let line2 = Line3 {
            origin: DVec3::new(-1.0, 1.5, 0.0),
            direction: DVec3::X,
        };
        let intervals2 = clip_line_to_polygon_with_tol(&line2, &plane, &verts, TOLERANCE_ABS);
        assert_eq!(intervals2.len(), 1, "line through narrow arm should have 1 interval");
        assert!((intervals2[0].0 - 1.0).abs() < TOLERANCE_MESH_LEGACY, "t_min={}", intervals2[0].0);
        assert!((intervals2[0].1 - 2.0).abs() < TOLERANCE_MESH_LEGACY, "t_max={}", intervals2[0].1);

        // Line through the notch (x = 1.5): only one interval through the lower arm
        // At x=1.5 the upper arm (x ∈ [0,1]) doesn't cover this x, so the line
        // enters at y=0, exits at y=1, and stays outside from y=1 onward.
        let line3 = Line3 {
            origin: DVec3::new(1.5, -1.0, 0.0),
            direction: DVec3::Y,
        };
        let intervals3 = clip_line_to_polygon_with_tol(&line3, &plane, &verts, TOLERANCE_ABS);
        assert_eq!(intervals3.len(), 1, "line through notch should have 1 interval");
        assert!((intervals3[0].0 - 1.0).abs() < TOLERANCE_MESH_LEGACY, "t_min={}", intervals3[0].0);
        assert!((intervals3[0].1 - 2.0).abs() < TOLERANCE_MESH_LEGACY, "t_max={}", intervals3[0].1);

        // Line that misses the polygon entirely
        let line4 = Line3 {
            origin: DVec3::new(5.0, 0.0, 0.0),
            direction: DVec3::Y,
        };
        let intervals4 = clip_line_to_polygon_with_tol(&line4, &plane, &verts, TOLERANCE_ABS);
        assert!(intervals4.is_empty(), "line outside should have 0 intervals");
    }
}
