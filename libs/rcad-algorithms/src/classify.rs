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
                let v = point - c.apex;
                let along = v.dot(c.axis);
                let perp = (v - c.axis * along).length();
                if (perp - along * c.half_angle_rad.tan()).abs() < on_surface_tol {
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
                let d = ray_dir - c.axis * ray_dir.dot(c.axis);
                let f = oc - c.axis * oc.dot(c.axis);
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
                // Cylinder primitives only have 2 seam vertices but their axis
                // projections give the correct height extent for the finite face.
                let face_verts = ds.face_boundary_points(fi);
                let axis = c.axis.normalize();
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
                let slack = boundary_tol;
                let sq = disc.sqrt();
                for &t in &[(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)] {
                    if t > ray_tol {
                        let hit = point + ray_dir * t;
                        let h = (hit - c.origin).dot(axis);
                        if h >= h_min - slack && h <= h_max + slack {
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
                let tan_a = c.half_angle_rad.tan();
                let co = point - c.apex;
                let d_along = ray_dir.dot(c.axis);
                let co_along = co.dot(c.axis);
                let d_perp = ray_dir - c.axis * d_along;
                let co_perp = co - c.axis * co_along;
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
                        let along = (hit - c.apex).dot(c.axis);
                        if along >= 0.0 {
                            let face_verts = ds.face_boundary_points(fi);
                            if point_in_face_aabb(hit, &face_verts, boundary_tol) {
                                crossings += 1;
                            }
                        }
                    }
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
