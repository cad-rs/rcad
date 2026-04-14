use glam::DVec3;
use rcad_kernel::geom::*;

use crate::bopds::ds::*;
use crate::inttools;
use crate::tolerance::{AdaptiveTolerance, ToleranceLevel};

/// Classification of a point relative to a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Classification {
    In,
    Out,
    On,
}

/// Classify a point relative to a solid (defined by its face indices in DS).
/// Uses ray-casting with multiple ray directions for robustness.
pub fn classify_point(point: DVec3, solid_face_indices: &[usize], ds: &DS) -> Classification {
    if solid_face_indices.is_empty() {
        return Classification::Out;
    }

    // Compute adaptive tolerance based on model scale
    let tol = AdaptiveTolerance::from_scale(ds.model_scale());
    let on_surface_tol = tol.tolerance(ToleranceLevel::Relaxed);

    if let Some(class) = classify_analytic_cylinder_solid(point, solid_face_indices, ds, on_surface_tol) {
        return class;
    }

    if let Some(class) = classify_analytic_cone_solid(point, solid_face_indices, ds, on_surface_tol) {
        return class;
    }

    // Fast path for an analytic torus solid: if all faces belong to the same torus
    // surface, classify by the torus implicit equation instead of ray casting.
    let mut analytic_torus: Option<rcad_kernel::geom::ToroidalSurface> = None;
    let mut all_same_torus = true;
    for &fi in solid_face_indices {
        match &ds.faces[fi].surface {
            Surface3::Torus(t) => {
                if let Some(base) = analytic_torus {
                    let same = (base.center - t.center).length() <= on_surface_tol * 10.0
                        && (base.axis.normalize_or_zero().dot(t.axis.normalize_or_zero())).abs() >= 0.9999
                        && (base.major_radius - t.major_radius).abs() <= on_surface_tol * 10.0
                        && (base.minor_radius - t.minor_radius).abs() <= on_surface_tol * 10.0;
                    if !same {
                        all_same_torus = false;
                        break;
                    }
                } else {
                    analytic_torus = Some(*t);
                }
            }
            _ => {
                all_same_torus = false;
                break;
            }
        }
    }
    if all_same_torus {
        if let Some(torus) = analytic_torus {
            let axis = torus.axis.normalize_or_zero();
            let delta = point - torus.center;
            let z = delta.dot(axis);
            let radial = delta - z * axis;
            let rho = radial.length();
            let tube_dist = ((rho - torus.major_radius).powi(2) + z * z).sqrt();
            let signed = tube_dist - torus.minor_radius;
            if signed.abs() <= on_surface_tol {
                return Classification::On;
            }
            return if signed < 0.0 {
                Classification::In
            } else {
                Classification::Out
            };
        }
    }

    // First check: is the point ON any face surface within face bounds?
    for &fi in solid_face_indices {
        let face = &ds.faces[fi];
        match &face.surface {
            Surface3::Plane(plane) => {
                let d = (point - plane.origin).dot(plane.normal);
                if d.abs() < on_surface_tol {
                    let face_verts = ds.face_boundary_points(fi);
                    if inttools::edge_face::point_in_planar_face(point, plane, &face_verts) {
                        return Classification::On;
                    }
                }
            }
            Surface3::Cylinder(c) => {
                let v = point - c.origin;
                let along = v.dot(c.axis);
                let perp = (v - c.axis * along).length();
                if (perp - c.radius).abs() < on_surface_tol {
                    return Classification::On;
                }
            }
            Surface3::Sphere(s) => {
                if ((point - s.center).length() - s.radius).abs() < on_surface_tol {
                    return Classification::On;
                }
            }
            Surface3::Cone(c) => {
                let axis = c.axis_dir();
                let v = point - c.apex_point();
                let along = v.dot(axis);
                let perp = (v - axis * along).length();
                if along >= -on_surface_tol && (perp - along.max(0.0) * c.half_angle_rad.tan()).abs() < on_surface_tol {
                    return Classification::On;
                }
            }
            _ => {}
        }
    }

    // Try multiple ray directions to avoid edge/vertex hits.
    // Use non-axis-aligned directions to minimize hitting polygon edges.
    let rays = [
        DVec3::new(0.8017, 0.2673, 0.5345).normalize(),
        DVec3::new(-0.3333, 0.6667, 0.6667).normalize(),
        DVec3::new(0.5774, -0.5774, 0.5774).normalize(),
    ];

    for ray_dir in &rays {
        match ray_cast_classify(point, *ray_dir, solid_face_indices, ds) {
            Some(class) => return class,
            None => continue, // ambiguous hit, try next ray
        }
    }

    // Fallback
    Classification::Out
}

fn classify_analytic_cylinder_solid(
    point: DVec3,
    solid_face_indices: &[usize],
    ds: &DS,
    on_surface_tol: f64,
) -> Option<Classification> {
    let mut cylinder: Option<rcad_kernel::geom::CylindricalSurface> = None;
    for &fi in solid_face_indices {
        match &ds.faces[fi].surface {
            Surface3::Cylinder(c) => {
                if let Some(base) = cylinder {
                    let same = (base.origin - c.origin).length() <= on_surface_tol * 10.0
                        && (base.axis.normalize_or_zero().dot(c.axis.normalize_or_zero())).abs() >= 0.9999
                        && (base.radius - c.radius).abs() <= on_surface_tol * 10.0;
                    if !same {
                        return None;
                    }
                } else {
                    cylinder = Some(*c);
                }
            }
            Surface3::Plane(_) => {}
            _ => return None,
        }
    }

    let cylinder = cylinder?;
    let axis = cylinder.axis.normalize_or_zero();
    let mut h_min = f64::INFINITY;
    let mut h_max = f64::NEG_INFINITY;
    for &fi in solid_face_indices {
        for v in ds.face_boundary_points(fi) {
            let h = (v - cylinder.origin).dot(axis);
            h_min = h_min.min(h);
            h_max = h_max.max(h);
        }
    }
    if !h_min.is_finite() || !h_max.is_finite() {
        return None;
    }

    let local = point - cylinder.origin;
    let along = local.dot(axis);
    let radial = (local - axis * along).length();

    if along < h_min - on_surface_tol || along > h_max + on_surface_tol {
        return Some(Classification::Out);
    }
    if radial > cylinder.radius + on_surface_tol {
        return Some(Classification::Out);
    }
    if (radial - cylinder.radius).abs() <= on_surface_tol
        || (along - h_min).abs() <= on_surface_tol
        || (along - h_max).abs() <= on_surface_tol
    {
        return Some(Classification::On);
    }
    Some(Classification::In)
}

fn classify_analytic_cone_solid(
    point: DVec3,
    solid_face_indices: &[usize],
    ds: &DS,
    on_surface_tol: f64,
) -> Option<Classification> {
    let mut cone: Option<rcad_kernel::geom::ConicalSurface> = None;
    for &fi in solid_face_indices {
        match &ds.faces[fi].surface {
            Surface3::Cone(c) => {
                if let Some(base) = cone {
                    let same = (base.apex_point() - c.apex_point()).length() <= on_surface_tol * 10.0
                        && (base.axis_dir().dot(c.axis_dir())).abs() >= 0.9999
                        && (base.half_angle_rad - c.half_angle_rad).abs() <= 1e-9;
                    if !same {
                        return None;
                    }
                } else {
                    cone = Some(*c);
                }
            }
            Surface3::Plane(_) => {}
            _ => return None,
        }
    }

    let cone = cone?;
    let axis = cone.axis_dir();
    let apex = cone.apex_point();
    let tan_half = cone.half_angle_rad.tan();
    if tan_half.abs() < 1e-12 {
        return None;
    }

    let mut along_min = f64::INFINITY;
    let mut along_max = f64::NEG_INFINITY;
    for &fi in solid_face_indices {
        for v in ds.face_boundary_points(fi) {
            let along = (v - apex).dot(axis);
            along_min = along_min.min(along);
            along_max = along_max.max(along);
        }
    }
    if !along_min.is_finite() || !along_max.is_finite() {
        return None;
    }

    let local = point - apex;
    let along = local.dot(axis);
    let radial = (local - axis * along).length();
    let allowed = along.max(0.0) * tan_half;

    if along < along_min - on_surface_tol || along > along_max + on_surface_tol {
        return Some(Classification::Out);
    }
    if radial > allowed + on_surface_tol {
        return Some(Classification::Out);
    }
    if (radial - allowed).abs() <= on_surface_tol
        || (along - along_min).abs() <= on_surface_tol
        || (along - along_max).abs() <= on_surface_tol
    {
        return Some(Classification::On);
    }
    Some(Classification::In)
}

/// Cast a single ray and count face crossings. Returns None if the ray hits
/// a face edge/vertex (ambiguous).
fn ray_cast_classify(
    point: DVec3,
    ray_dir: DVec3,
    solid_face_indices: &[usize],
    ds: &DS,
) -> Option<Classification> {
    let mut crossings = 0u32;
    let tol = AdaptiveTolerance::from_scale(ds.model_scale());
    let ray_tol = tol.tolerance(ToleranceLevel::Strict);
    let boundary_tol = tol.tolerance(ToleranceLevel::Relaxed);
    let parallel_tol_sq = tol.tolerance_sq(ToleranceLevel::Strict);

    for &fi in solid_face_indices {
        let face = &ds.faces[fi];
        match &face.surface {
            Surface3::Plane(plane) => {
                let denom = ray_dir.dot(plane.normal);
                if denom.abs() < ray_tol {
                    continue; // parallel to ray
                }
                let t = (plane.origin - point).dot(plane.normal) / denom;
                if t < ray_tol {
                    continue; // behind ray origin
                }

                let hit = point + ray_dir * t;

                // Check if hit is near a face boundary edge/vertex — if so, this ray is ambiguous
                let face_verts = ds.face_boundary_points(fi);
                if is_near_polygon_boundary(&hit, &face_verts, plane, boundary_tol) {
                    return None; // ambiguous, try different ray
                }

                if inttools::edge_face::point_in_planar_face(hit, plane, &face_verts) {
                    crossings += 1;
                }
            }
            Surface3::Cylinder(c) => {
                // Ray-cylinder intersection: |perp(ray_origin + t*ray_dir - origin)|² = r²
                let oc = point - c.origin;
                let axis = c.axis.normalize();
                let d = ray_dir - axis * ray_dir.dot(axis);
                let f = oc - axis * oc.dot(axis);
                let a = d.length_squared();
                if a < parallel_tol_sq {
                    continue; // ray parallel to cylinder axis
                }
                let b = 2.0 * d.dot(f);
                let cc = f.length_squared() - c.radius * c.radius;
                let disc = b * b - 4.0 * a * cc;
                if disc < 0.0 {
                    continue;
                }
                // Compute height range along cylinder axis from boundary verts.
                let face_verts = ds.face_boundary_points(fi);
                let (h_min, h_max) = if face_verts.len() >= 2 {
                    let mut mn = f64::INFINITY;
                    let mut mx = f64::NEG_INFINITY;
                    for &v in &face_verts {
                        let h = (v - c.origin).dot(axis);
                        mn = mn.min(h);
                        mx = mx.max(h);
                    }
                    (mn, mx)
                } else {
                    (-1e9, 1e9) // unbounded fallback
                };
                // Compute angular range of this cylinder face (finite arc vs full cylinder).
                // Build a reference perpendicular to the axis for angle measurement.
                let (angle_min, angle_max) = cylinder_face_angle_range(c, &face_verts, axis);
                let slack = boundary_tol;
                let sq = disc.sqrt();
                for &t in &[(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
                    if t > ray_tol {
                        let hit = point + ray_dir * t;
                        let h = (hit - c.origin).dot(axis);
                        if h >= h_min - slack && h <= h_max + slack {
                            // Check angular containment if the face is a partial arc
                            if angle_max - angle_min < std::f64::consts::TAU - 0.01 {
                                let radial = hit - c.origin - axis * h;
                                let angle = cylinder_angle(c, radial);
                                if !angle_in_range(angle, angle_min, angle_max, slack / c.radius) {
                                    continue;
                                }
                            }
                            crossings += 1;
                        }
                    }
                }
            }
            Surface3::Sphere(s) => {
                // Ray-sphere intersection
                let oc = point - s.center;
                let a = ray_dir.length_squared();
                let b = 2.0 * oc.dot(ray_dir);
                let cc = oc.length_squared() - s.radius * s.radius;
                let disc = b * b - 4.0 * a * cc;
                if disc < 0.0 {
                    continue;
                }
                let sq = disc.sqrt();
                for &t in &[(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
                    if t > ray_tol {
                        let hit = point + ray_dir * t;
                        let face_verts = ds.face_boundary_points(fi);
                        // A sphere primitive has at most 2 wire vertices (poles),
                        // which form a degenerate AABB. Treat the whole sphere
                        // surface as the face whenever we don't have a proper polygon.
                        let in_face = if face_verts.len() < 3 {
                            true // whole sphere is one face
                        } else {
                            point_in_face_aabb(hit, &face_verts, boundary_tol)
                        };
                        if in_face {
                            crossings += 1;
                        }
                    }
                }
            }
            Surface3::Cone(c) => {
                // Ray-cone intersection (finite cone approximated by AABB test)
                let axis = c.axis_dir();
                let tan_a = c.half_angle_rad.tan();
                let apex = c.apex_point();
                let co = point - apex;
                let d_along = ray_dir.dot(axis);
                let co_along = co.dot(axis);
                let d_perp = ray_dir - axis * d_along;
                let co_perp = co - axis * co_along;
                let a = d_perp.length_squared() - tan_a * tan_a * d_along * d_along;
                let b = 2.0 * (d_perp.dot(co_perp) - tan_a * tan_a * d_along * co_along);
                let cc = co_perp.length_squared() - tan_a * tan_a * co_along * co_along;
                let disc = b * b - 4.0 * a * cc;
                if a.abs() < parallel_tol_sq || disc < 0.0 {
                    continue;
                }
                let sq = disc.sqrt();
                for &t in &[(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
                    if t > ray_tol {
                        let hit = point + ray_dir * t;
                        let along = (hit - apex).dot(axis);
                        if along >= 0.0 {
                            let face_verts = ds.face_boundary_points(fi);
                            if point_in_face_aabb(hit, &face_verts, boundary_tol) {
                                crossings += 1;
                            }
                        }
                    }
                }
            }
            Surface3::Torus(t) => {
                // Algebraic ray-torus intersection via quartic polynomial.
                // For torus at center C, axis A, major radius R, minor radius r:
                //   (|X|² + R² − r²)² = 4R²·|X_radial|²   where X = P_local + tt·D
                let p_local = point - t.center;
                let axis = t.axis.normalize();
                let za = ray_dir.dot(axis);
                let zp = p_local.dot(axis);
                let p2 = p_local.length_squared();
                let hf = p_local.dot(ray_dir);
                let r_maj = t.major_radius;
                let r_min = t.minor_radius;

                let r2 = r_maj * r_maj;
                let e1 = p2 + r2 - r_min * r_min;
                let a1c = 1.0 - za * za;
                let b1 = 2.0 * (hf - za * zp);
                let c1 = p2 - zp * zp;

                let coeff3 = 4.0 * hf;
                let coeff2 = 4.0 * hf * hf + 2.0 * e1 - 4.0 * r2 * a1c;
                let coeff1_q = 4.0 * hf * e1 - 4.0 * r2 * b1;
                let coeff0 = e1 * e1 - 4.0 * r2 * c1;

                let q = |tt: f64| -> f64 {
                    (((tt + coeff3) * tt + coeff2) * tt + coeff1_q) * tt + coeff0
                };

                let face_verts = ds.face_boundary_points(fi);
                let t_max = if face_verts.is_empty() {
                    (r_maj + r_min) * 6.0
                } else {
                    face_verts
                        .iter()
                        .map(|&v| (v - point).length())
                        .fold(0.0_f64, f64::max)
                        + r_min * 2.0
                };

                const N_SCAN: usize = 64;
                let step = (t_max - ray_tol) / N_SCAN as f64;
                let mut t_prev = ray_tol;
                let mut q_prev = q(t_prev);

                for i in 1..=N_SCAN {
                    let tt = ray_tol + step * i as f64;
                    let q_curr = q(tt);
                    if q_prev * q_curr <= 0.0 {
                        let mut lo = t_prev;
                        let mut hi = tt;
                        for _ in 0..32 {
                            let mid = 0.5 * (lo + hi);
                            if q(lo) * q(mid) <= 0.0 { hi = mid; } else { lo = mid; }
                        }
                        let root = 0.5 * (lo + hi);
                        let hit = point + ray_dir * root;
                        let h_loc = hit - t.center;
                        let z_hit = h_loc.dot(axis);
                        let radial_sq = h_loc.length_squared() - z_hit * z_hit;
                        if radial_sq >= -1e-9
                            && point_in_face_aabb(hit, &face_verts, boundary_tol)
                        {
                            crossings += 1;
                        }
                    }
                    t_prev = tt;
                    q_prev = q_curr;
                }
            }
            _ => {}
        }
    }

    Some(if crossings % 2 == 1 {
        Classification::In
    } else {
        Classification::Out
    })
}

/// Compute the angular range [min, max] of a cylinder face in radians.
/// Returns (0, TAU) for a full cylinder face.
fn cylinder_face_angle_range(
    c: &rcad_kernel::geom::CylindricalSurface,
    face_verts: &[DVec3],
    axis: DVec3,
) -> (f64, f64) {
    if face_verts.len() < 2 {
        return (0.0, std::f64::consts::TAU);
    }
    let angles: Vec<f64> = face_verts
        .iter()
        .map(|&v| {
            let radial = v - c.origin - axis * (v - c.origin).dot(axis);
            cylinder_angle(c, radial)
        })
        .collect();
    let mut min_a = angles[0];
    let mut max_a = angles[0];
    for &a in &angles[1..] {
        if a < min_a { min_a = a; }
        if a > max_a { max_a = a; }
    }
    // If the angle span is essentially zero (all boundary vertices at same angle,
    // e.g. on the seam), we cannot determine the arc range — treat as full cylinder.
    if max_a - min_a < 1e-6 {
        return (0.0, std::f64::consts::TAU);
    }
    // If span > π, the face might wrap around. Detect wraparound:
    // if max-min > π, try wrapping angles and recompute.
    if max_a - min_a > std::f64::consts::PI {
        // Try normalizing angles relative to the first angle
        let ref_a = angles[0];
        let wrapped: Vec<f64> = angles
            .iter()
            .map(|&a| {
                let mut d = a - ref_a;
                while d < 0.0 { d += std::f64::consts::TAU; }
                while d > std::f64::consts::TAU { d -= std::f64::consts::TAU; }
                d
            })
            .collect();
        let span = wrapped.iter().cloned().fold(0.0_f64, f64::max);
        if span < std::f64::consts::TAU - 0.01 {
            // Not a full cylinder
            let wmin = 0.0_f64;
            let wmax = span;
            return (ref_a + wmin, ref_a + wmax);
        } else {
            return (0.0, std::f64::consts::TAU);
        }
    }
    // If the angular span is near zero, the boundary vertices all lie on the same
    // generator line (e.g. the seam), which means we cannot determine the actual
    // arc extent from vertex angles alone — treat as a full cylinder.
    if max_a - min_a < 1e-6 {
        return (0.0, std::f64::consts::TAU);
    }
    (min_a, max_a)
}

/// Compute the angle of a radial vector relative to the cylinder's reference direction.
fn cylinder_angle(c: &rcad_kernel::geom::CylindricalSurface, radial: DVec3) -> f64 {
    // Use any_perpendicular to get a reference direction
    let axis = c.axis.normalize();
    let ref_dir = rcad_kernel::geom::any_perpendicular(axis).normalize();
    let perp_dir = axis.cross(ref_dir).normalize();
    let x = radial.dot(ref_dir);
    let y = radial.dot(perp_dir);
    x.atan2(y) // in [-π, π]
}

/// Check if angle is within [min, max] range (with angular slack).
fn angle_in_range(angle: f64, min_a: f64, max_a: f64, slack: f64) -> bool {
    angle >= min_a - slack && angle <= max_a + slack
}

/// Conservative face containment check using AABB of the face boundary vertices.
fn point_in_face_aabb(point: DVec3, face_verts: &[DVec3], slack: f64) -> bool {
    if face_verts.is_empty() {
        return false;
    }
    let mut mn = face_verts[0];
    let mut mx = face_verts[0];
    for &v in face_verts.iter().skip(1) {
        mn = mn.min(v);
        mx = mx.max(v);
    }
    point.cmpge(mn - DVec3::splat(slack)).all() && point.cmple(mx + DVec3::splat(slack)).all()
}

/// Check if a point is close to any edge of a polygon (within tolerance).
fn is_near_polygon_boundary(point: &DVec3, verts: &[DVec3], plane: &Plane, boundary_tol: f64) -> bool {
    let (u_axis, v_axis) = inttools::edge_face::plane_local_basis(plane);
    let project = |p: DVec3| -> (f64, f64) {
        let d = p - plane.origin;
        (d.dot(u_axis), d.dot(v_axis))
    };

    let (px, py) = project(*point);
    let n = verts.len();
    let tol_sq = boundary_tol * boundary_tol;

    for i in 0..n {
        let j = (i + 1) % n;
        let (ax, ay) = project(verts[i]);
        let (bx, by) = project(verts[j]);

        // Distance from point to line segment (ax,ay)-(bx,by)
        let dx = bx - ax;
        let dy = by - ay;
        let len_sq = dx * dx + dy * dy;
        if len_sq < tol_sq {
            continue;
        }

        let t = ((px - ax) * dx + (py - ay) * dy) / len_sq;
        let t = t.clamp(0.0, 1.0);
        let cx = ax + t * dx;
        let cy = ay + t * dy;
        let dist_sq = (px - cx) * (px - cx) + (py - cy) * (py - cy);

        if dist_sq < tol_sq {
            return true;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom_populate::populate_box_geom;
    use rcad_kernel::{BRep, PrimitiveSolid};

    #[test]
    fn point_inside_box() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        populate_box_geom(&mut brep);
        let ds = DS::new(&brep, &BRep::new());
        let face_indices: Vec<usize> = (0..ds.faces.len())
            .filter(|&i| ds.faces[i].origin == ShapeOrigin::ShapeA)
            .collect();

        assert_eq!(
            classify_point(DVec3::new(0.5, 0.5, 0.5), &face_indices, &ds),
            Classification::In
        );
        assert_eq!(
            classify_point(DVec3::new(2.0, 0.5, 0.5), &face_indices, &ds),
            Classification::Out
        );
    }
}
