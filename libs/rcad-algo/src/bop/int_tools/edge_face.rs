// OCCT IntTools_EdgeFace — edge-face intersection
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Surface3};
#[derive(Debug, Clone)]
pub struct EdgeFaceHit { pub point: DVec3, pub edge_param: f64 }

pub fn intersect_line_plane(line: &rcad_kernel::geom::Line3, t_range: [f64; 2], plane: &rcad_kernel::geom::Plane) -> Option<EdgeFaceHit> {
    let denom = line.direction.dot(plane.normal);
    if denom.abs() < rcad_kernel::CONFUSION { return None; }
    let t = (plane.origin - line.origin).dot(plane.normal) / denom;
    if t < t_range[0] - rcad_kernel::CONFUSION || t > t_range[1] + rcad_kernel::CONFUSION { return None; }
    Some(EdgeFaceHit { point: line.origin + line.direction * t, edge_param: t })
}
pub fn plane_local_basis(plane: &rcad_kernel::geom::Plane) -> (DVec3, DVec3) {
    let n = plane.normal;
    let ref_dir = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    (n.cross(ref_dir).normalize(), n.cross(n.cross(ref_dir)).normalize())
}
pub fn point_in_planar_face(point: DVec3, plane: &rcad_kernel::geom::Plane, face_verts: &[DVec3]) -> bool {
    if face_verts.len() < 3 { return false; }
    let (u_axis, v_axis) = plane_local_basis(plane);
    let project = |p: DVec3| { let d = p - plane.origin; (d.dot(u_axis), d.dot(v_axis)) };
    let (px, py) = project(point);
    let poly: Vec<(f64, f64)> = face_verts.iter().map(|v| project(*v)).collect();
    let n = poly.len(); let mut inside = false; let mut j = n - 1;
    for i in 0..n {
        let (xi, yi) = poly[i]; let (xj, yj) = poly[j];
        if (yi > 0.0) != (yj > 0.0) {
            let xint = (xj - xi) * (0.0 - yi) / (yj - yi) + xi;
            if px < xint { inside = !inside; }
        }
        j = i;
    }
    inside
}
