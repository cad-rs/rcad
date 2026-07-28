use crate::bop::ds::ds::DS;
use crate::tolerance::*;
use glam::DVec3;
use rcad_kernel::geom::*;

pub struct EdgeFaceHit {
    pub point: DVec3,
    pub edge_param: f64,
}

/// Intersect a line segment (bounded by t_range) with a plane.
/// Does NOT check face boundary containment  ?caller must do that.
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
/// Use the same magnitude as pave [`crate::bop::ds::ds::DS::fuzzy_tol`] for consistent V–F containment.
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

///  ?compute_edge_face_criteria (IntTools_EdgeFace.cxx L528-548).
/// Computes the tolerance sum for edge-face intersection.
/// For BSpline/Bezier curves with large tolerance ratio, uses max.
pub fn compute_edge_face_criteria(edge_tol: f64, face_tol: f64, curve_type: &Curve3) -> f64 {
    let fuzz = 0.0;
    let a_tol_f = face_tol + fuzz;
    let a_tol_e = edge_tol + fuzz;
    match curve_type {
        Curve3::BSpline(_) | Curve3::Bezier(_) => {
            let diff1 = a_tol_e / a_tol_f.max(TOLERANCE_LEN_SQ_DIV_SAFE);
            let diff2 = a_tol_f / a_tol_e.max(TOLERANCE_LEN_SQ_DIV_SAFE);
            if diff1 > 100.0 || diff2 > 100.0 {
                a_tol_e.max(a_tol_f)
            } else {
                1.5 * a_tol_e + a_tol_f
            }
        }
        _ => a_tol_e + a_tol_f,
    }
}

///  ?IsEqDistance (IntTools_EdgeFace.cxx L240-299).
/// Checks if point is near the axis of a cylindrical/conical/toroidal surface,
/// returning the surface's radius at that point.
pub fn is_eq_distance(p: DVec3, surface: &Surface3, tol: f64) -> Option<f64> {
    match surface {
        Surface3::Cylinder(c) => {
            let v = p - c.origin;
            let proj = v.dot(c.axis);
            let radial = v - c.axis * proj;
            let dist = radial.length();
            if dist < tol { Some(c.radius) } else { None }
        }
        Surface3::Cone(cn) => {
            let axis = cn.axis_dir();
            let v = p - cn.apex;
            let proj = v.dot(axis);
            let radial = v - axis * proj;
            let dist = radial.length();
            if dist < tol {
                let r_at_z = cn.radius + proj * cn.half_angle_rad.tan();
                Some(r_at_z.abs())
            } else {
                None
            }
        }
        Surface3::Torus(t) => {
            let d_center = (p - t.center).length();
            let dc = (d_center - t.major_radius).abs();
            if dc < tol { Some(t.minor_radius) } else { None }
        }
        _ => None,
    }
}

///  ?IsCoincident (IntTools_EdgeFace.cxx L62-163).
/// Checks if an edge is coincident with a face by sampling points along the
/// edge, projecting onto the face, and classifying them.
/// Returns true if >50% of sample points project within criteria AND are IN.
pub fn is_coincident_edge_face(
    curve: &Curve3,
    t_range: [f64; 2],
    surface: &Surface3,
    face_tol: f64,
    edge_tol: f64,
    context: &mut crate::bop::int_tools::context::Context,
    ds: &DS,
    face_idx: usize,
) -> bool {
    let criteria = compute_edge_face_criteria(edge_tol, face_tol, curve);
    let a_tresh = 0.5;
    let a_nb_seg = if matches!(curve, Curve3::Line(_)) && matches!(surface, Surface3::Plane(_)) {
        2
    } else {
        23
    };
    let a_tresh_idx_f = ((a_nb_seg + 1) as f64 * 0.25) as i32;
    let a_tresh_idx_l = ((a_nb_seg + 1) as f64 * 0.75) as i32;

    let mut t1 = t_range[0];
    let mut t2 = t_range[1];
    let bnd_shift = 0.01 * (t2 - t1);
    t1 += bnd_shift;
    t2 -= bnd_shift;
    let dt = (t2 - t1) / a_nb_seg as f64;

    let mut i_cnt = 0i32;
    let mut is_classified = false;

    for i in 0..=a_nb_seg {
        let t = t1 + i as f64 * dt;
        let p = curve.point_at(t);
        let proj = context.proj_ps(ds, face_idx, p);
        let Some((uv, _pt3d, dist)) = proj else {
            continue;
        };
        if dist > criteria {
            if dist > 100.0 * criteria {
                return false;
            }
            continue;
        }
        i_cnt += 1;
        if ((0 < i) && (i < a_tresh_idx_f)) || ((a_tresh_idx_l < i) && (i < a_nb_seg)) {
            continue;
        }
        if is_classified && (i != a_nb_seg) {
            continue;
        }
        let state = context.fclass2d(ds, face_idx).perform(uv, true);
        use crate::bop::int_tools::fclass2d::State;
        if state == State::Out {
            return false;
        }
        if i != 0 {
            is_classified = true;
        }
    }

    let coeff = i_cnt as f64 / (a_nb_seg + 1) as f64;
    coeff > a_tresh
}

///  ?IsCoplanar (IntTools_EdgeFace.cxx L788-813).
/// Checks if a curve lies in the plane of a planar surface.
pub fn is_coplanar(curve: &Curve3, surface: &Surface3) -> bool {
    let Surface3::Plane(pl) = surface else {
        return false;
    };
    match curve {
        Curve3::Line(l) => l.direction.dot(pl.normal).abs() < 1e-12,
        Curve3::Circle(c) => c.normal.dot(pl.normal).abs() > 1.0 - 1e-12,
        Curve3::Ellipse(e) => e.normal.dot(pl.normal).abs() > 1.0 - 1e-12,
        _ => false,
    }
}

///  ?IsRadius (IntTools_EdgeFace.cxx L815-843).
/// Checks if a curve's radius matches the surface's curvature radius.
pub fn is_radius(curve: &Curve3, surface: &Surface3, criteria: f64) -> bool {
    match (curve, surface) {
        (Curve3::Circle(c), Surface3::Sphere(s)) => {
            let dist = (c.center - s.center).length();
            (dist - s.radius).abs() < criteria
        }
        (Curve3::Circle(c), Surface3::Cylinder(cyl)) => {
            let v = c.center - cyl.origin;
            let proj = v.dot(cyl.axis);
            let radial = v - cyl.axis * proj;
            let axis_dist = radial.length();
            (axis_dist - c.radius).abs() < criteria
        }
        _ => false,
    }
}

///  ?MakeType (IntTools_EdgeFace.cxx L304-359).
/// Determines whether a common part is EDGE or VERTEX type.
pub fn make_edge_face_type(
    edge_t_range: [f64; 2],
    common_t_range: [f64; 2],
    curve: &Curve3,
    criteria: f64,
    is_whole_range: bool,
) -> (i8, f64) {
    let [af1, al1] = common_t_range;
    let [ef1, el1] = edge_t_range;
    let pf = curve.point_at(af1);
    let pl = curve.point_at(al1);
    let df1 = (pf - pl).length();

    if (df1 > criteria * 2.0) && is_whole_range {
        return (1, 0.0); // TopAbs_EDGE
    }

    if is_whole_range {
        let tm = (af1 + al1) * 0.5;
        let dist = (pf - curve.point_at(tm)).length();
        if dist > criteria * 2.0 {
            return (1, 0.0); // TopAbs_EDGE
        }
    }

    let tm = (af1 + al1) * 0.5;
    (0, tm) // TopAbs_VERTEX
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

        // line_dir × edge_dir  ?zero means parallel
        let denom = line_u * ey - line_v * ex;
        if denom.abs() < eps {
            // when the line coincides with a polygon edge (parallel
            // AND zero distance), use the edge endpoints as intersection t-values.
            // Without this, boundary-coincident intersection lines lose the portion
            // of the line that runs along the polygon edge (e.g. bfuse_simple B3
            // face[6] intersection at z=0, which is on the face boundary).
            let dist = (origin_u - ax) * ey - (origin_v - ay) * ex;
            if dist.abs() < eps {
                // Line coincides with this edge  ?add t at both endpoints
                let dir_len2 = line_u * line_u + line_v * line_v;
                if dir_len2 > TOLERANCE_LEN_SQ_DIV_SAFE {
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

        // Edge parameter s  ?check the intersection is within the segment
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
        // Line doesn't cross any edge  ?check if origin is inside
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
