use glam::DVec3;
use rcad_kernel::geom::*;

use crate::bopds::ds::*;
use crate::inttools;
use crate::tolerance::*;

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

    // First check: is the point ON any face surface within face bounds?
    for &fi in solid_face_indices {
        let face = &ds.faces[fi];
        match &face.surface {
            Surface3::Plane(plane) => {
                let d = (point - plane.origin).dot(plane.normal);
                if d.abs() < TOLERANCE_ABS * 100.0 {
                    let face_verts = ds.face_boundary_points(fi);
                    if inttools::edge_face::point_in_planar_face(point, plane, &face_verts) {
                        return Classification::On;
                    }
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

    for &fi in solid_face_indices {
        let face = &ds.faces[fi];
        match &face.surface {
            Surface3::Plane(plane) => {
                let denom = ray_dir.dot(plane.normal);
                if denom.abs() < TOLERANCE_ABS {
                    continue; // parallel to ray
                }
                let t = (plane.origin - point).dot(plane.normal) / denom;
                if t < TOLERANCE_ABS {
                    continue; // behind ray origin
                }

                let hit = point + ray_dir * t;

                // Check if hit is near a face boundary edge/vertex — if so, this ray is ambiguous
                let face_verts = ds.face_boundary_points(fi);
                if is_near_polygon_boundary(&hit, &face_verts, plane) {
                    return None; // ambiguous, try different ray
                }

                if inttools::edge_face::point_in_planar_face(hit, plane, &face_verts) {
                    crossings += 1;
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

/// Check if a point is close to any edge of a polygon (within tolerance).
fn is_near_polygon_boundary(point: &DVec3, verts: &[DVec3], plane: &Plane) -> bool {
    let (u_axis, v_axis) = inttools::edge_face::plane_local_basis(plane);
    let project = |p: DVec3| -> (f64, f64) {
        let d = p - plane.origin;
        (d.dot(u_axis), d.dot(v_axis))
    };

    let (px, py) = project(*point);
    let n = verts.len();
    let boundary_tol = TOLERANCE_ABS * 1000.0; // generous boundary tolerance

    for i in 0..n {
        let j = (i + 1) % n;
        let (ax, ay) = project(verts[i]);
        let (bx, by) = project(verts[j]);

        // Distance from point to line segment (ax,ay)-(bx,by)
        let dx = bx - ax;
        let dy = by - ay;
        let len_sq = dx * dx + dy * dy;
        if len_sq < TOLERANCE_ABS * TOLERANCE_ABS {
            continue;
        }

        let t = ((px - ax) * dx + (py - ay) * dy) / len_sq;
        let t = t.clamp(0.0, 1.0);
        let cx = ax + t * dx;
        let cy = ay + t * dy;
        let dist_sq = (px - cx) * (px - cx) + (py - cy) * (py - cy);

        if dist_sq < boundary_tol * boundary_tol {
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
