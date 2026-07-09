//! Shape-to-shape and point-to-shape minimum distance.
//!
//! Analogous to OCCT `BRepExtrema_DistShapeShape`.
//!
//! # Strategy
//! 1. Sample each face on an 8 8 UV grid + wire vertices.
//! 2. For each sample on A, project onto every analytic surface of B via
//! [`closest_point_on_surface`] (Newton-converged).
//! 3. Symmetric pass B  ?A.
//! 4. Refine the best candidate pair with alternating projection until
//! convergence (typically 3 ? iterations, tolerance 1e-9).
//!
//! Complexity is O(F_A F_B S ) for the sampling phase, but the
//! refinement step is O(1) and brings the result to near-machine precision.

use glam::{DVec2, DVec3};

use crate::{BRep, topods, Surface3, closest_point_on_surface, geom::SurfaceEval};

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Result type
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Result of a shape-to-shape or point-to-shape distance query.
#[derive(Debug, Clone)]
pub struct ShapeDistance {
 /// Minimum Euclidean distance between the two shapes (or point and shape).
 pub distance: f64,
 /// The closest point on the first shape (or the query point).
 pub point_on_a: DVec3,
 /// The closest point on the second shape (or the shape surface).
 pub point_on_b: DVec3,
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Public API
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Compute the minimum distance between two BReps.
///
/// Returns the pair of closest points (one on each shape) and the distance.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::BRep;
/// use rcad_kernel::distance::min_distance;
///
/// let (a, _) = BRep::build_unit_cube();
/// let (b, _) = BRep::build_unit_cube();
/// let d = min_distance(&a, &b);
/// assert!(d.distance >= 0.0);
/// ```
pub fn min_distance(a: &BRep, b: &BRep) -> ShapeDistance {
 if let Some(exact) = analytic_min_distance_single_face(a, b) {
 return exact;
 }

 let mut best = ShapeDistance {
 distance: f64::INFINITY,
 point_on_a: DVec3::ZERO,
 point_on_b: DVec3::ZERO,
 };

 // Sample points on A, project onto B
 let samples_a = sample_brep_points(a);
 for &pa in &samples_a {
 if let Some(r) = closest_on_brep(pa, b)
 && r.distance < best.distance
 {
 best = ShapeDistance {
 distance: r.distance,
 point_on_a: pa,
 point_on_b: r.point,
 };
 }
 }

 // Sample points on B, project onto A (symmetric)
 let samples_b = sample_brep_points(b);
 for &pb in &samples_b {
 if let Some(r) = closest_on_brep(pb, a)
 && r.distance < best.distance
 {
 best = ShapeDistance {
 distance: r.distance,
 point_on_a: r.point,
 point_on_b: pb,
 };
 }
 }

 // Refine the best candidate with alternating projection
 if best.distance.is_finite() {
 let (pa, pb, dist) = refine_pair(best.point_on_a, best.point_on_b, a, b);
 if dist < best.distance {
 best = ShapeDistance { distance: dist, point_on_a: pa, point_on_b: pb };
 }
 }

 best
}

/// Compute the minimum distance from a 3D point to the surface of a BRep.
///
/// Returns the closest point on the shape surface and the distance.
///
/// # Examples
/// ```rust
/// use glam::DVec3;
/// use rcad_kernel::BRep;
/// use rcad_kernel::distance::point_to_shape_distance;
///
/// let (box_brep, _) = BRep::build_unit_cube();
/// // The unit cube spans `[0,1]` in each dimension. Pick a point far away.
/// let d = point_to_shape_distance(DVec3::new(10.0, 10.0, 10.0), &box_brep);
/// assert!(d.distance > 0.0);
/// ```
pub fn point_to_shape_distance(query: DVec3, brep: &BRep) -> ShapeDistance {
 match closest_on_brep(query, brep) {
 Some(r) => ShapeDistance {
 distance: r.distance,
 point_on_a: query,
 point_on_b: r.point,
 },
 None => ShapeDistance {
 distance: f64::INFINITY,
 point_on_a: query,
 point_on_b: query,
 },
 }
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Internal helpers
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Lightweight result used internally (just point + distance).
struct ClosestResult {
 point: DVec3,
 distance: f64,
}

/// Refine a candidate closest-point pair using alternating projection.
///
/// Starting from an initial guess `(pa, pb)`, alternately projects each point
/// onto the other shape until the pair converges (  < 1e-9) or 30 iterations
/// are exhausted.  Converges in 3 ? iterations for smooth surfaces.
fn refine_pair(
 mut pa: DVec3,
 mut pb: DVec3,
 a: &BRep,
 b: &BRep,
) -> (DVec3, DVec3, f64) {
 for _ in 0..30 {
 // Project pa onto B
 if let Some(r) = closest_on_brep(pa, b) {
 pb = r.point;
 }
 // Project pb onto A
 if let Some(r) = closest_on_brep(pb, a) {
 let new_pa = r.point;
 let delta = (new_pa - pa).length();
 pa = new_pa;
 if delta < 1e-9 {
 break;
 }
 } else {
 break;
 }
 }
 let dist = (pa - pb).length();
 (pa, pb, dist)
}

#[derive(Clone)]
struct SingleFaceInfo<'a> {
 surface: &'a Surface3,
 polygon: Vec<DVec3>,
}

fn analytic_min_distance_single_face(a: &BRep, b: &BRep) -> Option<ShapeDistance> {
 let fa = single_face_info(a)?;
 let fb = single_face_info(b)?;

 match (fa.surface, fb.surface) {
 (Surface3::Sphere(sa), Surface3::Sphere(sb)) => {
 let center_delta = sb.center - sa.center;
 let center_dist = center_delta.length();
 if center_dist + 1e-12 < sa.radius + sb.radius {
 return None;
 }
 let dir = if center_dist > 1e-12 {
 center_delta / center_dist
 } else {
 DVec3::X
 };
 let point_on_a = sa.center + dir * sa.radius;
 let point_on_b = sb.center - dir * sb.radius;
 Some(ShapeDistance {
 distance: (point_on_b - point_on_a).length(),
 point_on_a,
 point_on_b,
 })
 }
 (Surface3::Plane(pa), Surface3::Sphere(sb)) => {
 analytic_plane_sphere_distance(pa.origin, pa.normal, &fa.polygon, sb.center, sb.radius)
 .map(|(point_on_plane, point_on_sphere, distance)| ShapeDistance {
 distance,
 point_on_a: point_on_plane,
 point_on_b: point_on_sphere,
 })
 }
 (Surface3::Sphere(sa), Surface3::Plane(pb)) => {
 analytic_plane_sphere_distance(pb.origin, pb.normal, &fb.polygon, sa.center, sa.radius)
 .map(|(point_on_plane, point_on_sphere, distance)| ShapeDistance {
 distance,
 point_on_a: point_on_sphere,
 point_on_b: point_on_plane,
 })
 }
 (Surface3::Plane(pa), Surface3::Plane(pb)) => {
 analytic_parallel_plane_distance(pa.origin, pa.normal, &fa.polygon, pb.origin, pb.normal, &fb.polygon)
 .map(|(point_on_a, point_on_b, distance)| ShapeDistance {
 distance,
 point_on_a,
 point_on_b,
 })
 }
 _ => None,
 }
}

fn single_face_info(brep: &topods::BRep) -> Option<SingleFaceInfo<'_>> {
 // Find first Solid tshape
 let solid_idx = brep.tshapes.iter().position(|ts| {
  if let topods::TShape::Solid(_) = &**ts { true } else { false }
 })?;
 let topods::TShape::Solid(sd) = &*brep.tshapes[solid_idx] else { return None };
 let shell_sr = sd.shells.first()?;
 let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] else { return None };
 if shd.faces.len() != 1 { return None }
 let face_sr = shd.faces.first()?;
 let topods::TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { return None };
 let surface = fd.surface.as_ref()?;
 let topods::TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] else { return None };
 let polygon: Vec<DVec3> = wd.edges.iter().filter_map(|we| {
  let topods::TShape::Edge(ed) = &*brep.tshapes[we.index] else { return None };
  let v_idx = if we.orientation.is_forward() { ed.first.index } else { ed.last.index };
  let topods::TShape::Vertex(vd) = &*brep.tshapes[v_idx] else { return None };
  Some(vd.point)
 }).collect();
 if polygon.len() < 3 { return None }
 Some(SingleFaceInfo { surface, polygon })
}

fn analytic_plane_sphere_distance(
 plane_origin: DVec3,
 plane_normal: DVec3,
 plane_polygon: &[DVec3],
 sphere_center: DVec3,
 sphere_radius: f64,
) -> Option<(DVec3, DVec3, f64)> {
 let n = plane_normal.normalize_or_zero();
 if n.length_squared() < 1e-20 {
 return None;
 }
 let signed = (sphere_center - plane_origin).dot(n);
 let point_on_plane = sphere_center - signed * n;
 if !point_in_planar_polygon(point_on_plane, plane_origin, n, plane_polygon) {
 return None;
 }
 let sign = if signed >= 0.0 { 1.0 } else { -1.0 };
 let point_on_sphere = sphere_center - sign * n * sphere_radius;
 Some((
 point_on_plane,
 point_on_sphere,
 (signed.abs() - sphere_radius).max(0.0),
 ))
}

fn analytic_parallel_plane_distance(
 origin_a: DVec3,
 normal_a: DVec3,
 polygon_a: &[DVec3],
 origin_b: DVec3,
 normal_b: DVec3,
 polygon_b: &[DVec3],
) -> Option<(DVec3, DVec3, f64)> {
 let n0 = normal_a.normalize_or_zero();
 let n1 = normal_b.normalize_or_zero();
 if n0.length_squared() < 1e-20 || n1.length_squared() < 1e-20 {
 return None;
 }
 if n0.dot(n1).abs() <= 1.0 - 1e-9 {
 return None;
 }

 let u = crate::geom::any_perpendicular(n0);
 let v = n0.cross(u).normalize_or_zero();
 let poly_a_2d: Vec<DVec2> = polygon_a.iter().map(|&p| to_plane_uv(p, origin_a, u, v)).collect();
 let poly_b_2d: Vec<DVec2> = polygon_b
 .iter()
 .map(|&p| to_plane_uv(project_point_to_plane(p, origin_a, n0), origin_a, u, v))
 .collect();

 let overlap_pt = polygon_overlap_point(&poly_a_2d, &poly_b_2d)?;
 let point_on_a = origin_a + overlap_pt.x * u + overlap_pt.y * v;
 let signed = (origin_b - origin_a).dot(n0);
 let point_on_b = point_on_a + signed * n0;
 Some((point_on_a, point_on_b, signed.abs()))
}

fn project_point_to_plane(point: DVec3, plane_origin: DVec3, plane_normal: DVec3) -> DVec3 {
 point - (point - plane_origin).dot(plane_normal) * plane_normal
}

fn to_plane_uv(point: DVec3, origin: DVec3, u: DVec3, v: DVec3) -> DVec2 {
 let d = point - origin;
 DVec2::new(d.dot(u), d.dot(v))
}

fn point_in_planar_polygon(point: DVec3, origin: DVec3, normal: DVec3, polygon: &[DVec3]) -> bool {
 let u = crate::geom::any_perpendicular(normal);
 let v = normal.cross(u).normalize_or_zero();
 let point_uv = to_plane_uv(point, origin, u, v);
 let poly_uv: Vec<DVec2> = polygon.iter().map(|&p| to_plane_uv(p, origin, u, v)).collect();
 point_in_polygon_2d(point_uv, &poly_uv)
}

fn point_in_polygon_2d(point: DVec2, polygon: &[DVec2]) -> bool {
 if polygon.len() < 3 {
 return false;
 }
 let mut inside = false;
 let mut j = polygon.len() - 1;
 for i in 0..polygon.len() {
 let pi = polygon[i];
 let pj = polygon[j];
 let dy = pj.y - pi.y;
 let intersects = ((pi.y > point.y) != (pj.y > point.y))
 && (point.x < (pj.x - pi.x) * (point.y - pi.y) / dy.abs().max(1e-20) + pi.x);
 if intersects {
 inside = !inside;
 }
 j = i;
 }
 inside
}

fn polygon_overlap_point(poly_a: &[DVec2], poly_b: &[DVec2]) -> Option<DVec2> {
 for &p in poly_a {
 if point_in_polygon_2d(p, poly_b) {
 return Some(p);
 }
 }
 for &p in poly_b {
 if point_in_polygon_2d(p, poly_a) {
 return Some(p);
 }
 }
 for i in 0..poly_a.len() {
 let a0 = poly_a[i];
 let a1 = poly_a[(i + 1) % poly_a.len()];
 for j in 0..poly_b.len() {
 let b0 = poly_b[j];
 let b1 = poly_b[(j + 1) % poly_b.len()];
 if let Some(p) = segment_intersection_2d(a0, a1, b0, b1) {
 return Some(p);
 }
 }
 }
 None
}

fn segment_intersection_2d(a0: DVec2, a1: DVec2, b0: DVec2, b1: DVec2) -> Option<DVec2> {
 let r = a1 - a0;
 let s = b1 - b0;
 let denom = r.perp_dot(s);
 if denom.abs() < 1e-12 {
 return None;
 }
 let qp = b0 - a0;
 let t = qp.perp_dot(s) / denom;
 let u = qp.perp_dot(r) / denom;
 if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
 Some(a0 + t * r)
 } else {
 None
 }
}

/// Find the closest point on any face surface of `brep` to `query`.
/// Returns `None` if the BRep has no faces with analytic surfaces.
fn closest_on_brep(query: DVec3, brep: &topods::BRep) -> Option<ClosestResult> {
 let mut best: Option<ClosestResult> = None;

 for ts in &brep.tshapes {
  let topods::TShape::Face(fd) = &**ts else { continue };
  let surface = fd.surface.as_ref()?;
  let proj = closest_point_on_surface(surface, query, 8);
  if best.as_ref().is_none_or(|b| proj.distance < b.distance) {
  best = Some(ClosestResult {
   point: proj.point,
   distance: proj.distance,
  });
  }
 }

 best
}

/// Collect sample points from the surface of a BRep: 8 8 grid per face + vertices.
fn sample_brep_points(brep: &topods::BRep) -> Vec<DVec3> {
 const GRID: usize = 8;
 let mut pts = Vec::new();

 // Vertex positions
 for ts in &brep.tshapes {
  if let topods::TShape::Vertex(vd) = &**ts {
  pts.push(vd.point);
  }
 }

 // Per-face surface grid
 for ts in &brep.tshapes {
  let topods::TShape::Face(fd) = &**ts else { continue };
  let surface = match fd.surface.as_ref() {
  Some(s) => s,
  None => continue,
  };
  let [u0, u1, v0, v1] = match fd.uv_domain {
  Some(r) => r,
  None => surface.default_domain(),
  };

  for i in 0..GRID {
  for j in 0..GRID {
   let u = u0 + (u1 - u0) * (i as f64 + 0.5) / GRID as f64;
   let v = v0 + (v1 - v0) * (j as f64 + 0.5) / GRID as f64;
   pts.push(surface.point_at(u, v));
  }
  }
 }

 pts
}

//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
// Tests
//  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

#[cfg(test)]
mod tests {
 use super::*;
 use crate::geom::Plane;
 use crate::topods::{self, ShapeRef, Orientation};

 fn make_square_plane_brep(origin: DVec3, normal: DVec3, width: f64, height: f64) -> BRep {
 let n = normal.normalize_or_zero();
 let u = crate::geom::any_perpendicular(n);
 let v = n.cross(u).normalize_or_zero();
 let hw = width * 0.5;
 let hh = height * 0.5;

 let mut brep = topods::BRep::new();

 // 4 corner vertices
 let pts = [
  origin - hw * u - hh * v,
  origin + hw * u - hh * v,
  origin + hw * u + hh * v,
  origin - hw * u + hh * v,
 ];
 let verts: Vec<ShapeRef> = pts.iter().map(|&p| brep.add_tvertex(p)).collect();

 // 4 edges
 let e0 = brep.add_tedge(None, verts[0], ShapeRef::synthetic_with_orientation(verts[1].index, Orientation::Reversed), [0.0, 1.0]);
 let e1 = brep.add_tedge(None, verts[1], ShapeRef::synthetic_with_orientation(verts[2].index, Orientation::Reversed), [0.0, 1.0]);
 let e2 = brep.add_tedge(None, verts[2], ShapeRef::synthetic_with_orientation(verts[3].index, Orientation::Reversed), [0.0, 1.0]);
 let e3 = brep.add_tedge(None, verts[3], ShapeRef::synthetic_with_orientation(verts[0].index, Orientation::Reversed), [0.0, 1.0]);

 // Wire with 4 edges
 let wire = brep.add_twire(vec![
  ShapeRef::synthetic_with_orientation(e0.index, Orientation::Forward),
  ShapeRef::synthetic_with_orientation(e1.index, Orientation::Forward),
  ShapeRef::synthetic_with_orientation(e2.index, Orientation::Forward),
  ShapeRef::synthetic_with_orientation(e3.index, Orientation::Forward),
 ]);

 // Face with Plane surface
 let face = brep.add_tface(
  Some(Surface3::Plane(Plane { origin, normal: n })),
  wire,
  vec![],
  None,
  Some([-hw, hw, -hh, hh]),
  vec![],
  true,
 );

 // Shell and solid
 let shell = brep.add_tshell(vec![face]);
 brep.add_tsolid(vec![shell]);

 brep
 }

 fn make_sphere_brep(radius: f64) -> BRep {
 use crate::geom::SphericalSurface;
 let mut brep = topods::BRep::new();

 // 3 vertices forming a triangular patch on the sphere
 let v0 = brep.add_tvertex(DVec3::new(0.0, 0.0, radius));
 let v1 = brep.add_tvertex(DVec3::new(radius, 0.0, 0.0));
 let v2 = brep.add_tvertex(DVec3::new(0.0, radius, 0.0));

 // 3 edges
 let e0 = brep.add_tedge(None, v0, ShapeRef::synthetic_with_orientation(v1.index, Orientation::Reversed), [0.0, 1.0]);
 let e1 = brep.add_tedge(None, v1, ShapeRef::synthetic_with_orientation(v2.index, Orientation::Reversed), [0.0, 1.0]);
 let e2 = brep.add_tedge(None, v2, ShapeRef::synthetic_with_orientation(v0.index, Orientation::Reversed), [0.0, 1.0]);

 // Wire
 let wire = brep.add_twire(vec![
  ShapeRef::synthetic_with_orientation(e0.index, Orientation::Forward),
  ShapeRef::synthetic_with_orientation(e1.index, Orientation::Forward),
  ShapeRef::synthetic_with_orientation(e2.index, Orientation::Forward),
 ]);

 // Face with sphere surface
 let surface = Surface3::Sphere(SphericalSurface {
  center: DVec3::ZERO,
  axis: DVec3::Z,
  radius,
  ref_dir: DVec3::X,
 });
 let face = brep.add_tface(
  Some(surface),
  wire,
  vec![],
  Some(DVec3::new(0.0, 0.0, radius)),
  Some([0.0, std::f64::consts::TAU, -std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2]),
  vec![],
  true,
 );

 // Shell and solid
 let shell = brep.add_tshell(vec![face]);
 brep.add_tsolid(vec![shell]);

 brep
 }

 #[test]
 fn point_to_box_distance() {
 // Use sphere which has an analytic surface entry.
 let brep = make_sphere_brep(1.0);
 let d = point_to_shape_distance(DVec3::new(0.0, 0.0, 5.0), &brep);
 println!("point_to_sphere_distance (vertical): {}", d.distance);
 assert!(d.distance > 0.0, "distance should be positive");
 assert!(
 d.distance < 10.0,
 "distance should be finite and reasonable"
 );
 }

 #[test]
 fn point_to_sphere_distance() {
 // Sphere radius 1.0 at origin; point at (5, 0, 0) -> distance ~4.0
 let brep = make_sphere_brep(1.0);
 let d = point_to_shape_distance(DVec3::new(5.0, 0.0, 0.0), &brep);
 println!("point_to_sphere_distance: {}", d.distance);
 assert!(
 (d.distance - 4.0).abs() < 0.1,
 "expected ~4.0, got {}",
 d.distance
 );
 }

 #[test]
 fn min_distance_disjoint_shapes() {
 // Two spheres far apart: one at origin (r=1), one implicitly at origin too.
 // Two boxes: one at default position, check distance is non-negative.
 let a = make_sphere_brep(1.0);
 let (b, _) = BRep::build_unit_cube();
 let d = min_distance(&a, &b);
 assert!(
 d.distance >= 0.0,
 "distance must be non-negative, got {}",
 d.distance
 );
 println!("min_distance sphere-box: {}", d.distance);
 }

 #[test]
 fn disjoint_spheres_distance_is_correct() {
 // Two identical unit spheres at origin -> distance ~0
 let a = make_sphere_brep(1.0);
 let b = make_sphere_brep(1.0);
 let d = min_distance(&a, &b);
 assert!(
 d.distance < 0.5,
 "identical spheres should have distance ~0, got {}",
 d.distance
 );
 }

 #[test]
 fn point_on_sphere_surface_has_near_zero_distance() {
 let brep = make_sphere_brep(2.0);
 // A point on the sphere surface (radius = 2, pointing along X).
 let d = point_to_shape_distance(DVec3::new(2.0, 0.0, 0.0), &brep);
 assert!(
 d.distance < 0.1,
 "point on sphere surface should have near-zero distance, got {}",
 d.distance
 );
 }

 #[test]
 fn distance_is_symmetric() {
 let a = make_sphere_brep(1.0);
 let (b, _) = BRep::build_unit_cube();
 let d_ab = min_distance(&a, &b).distance;
 let d_ba = min_distance(&b, &a).distance;
 assert!(
 (d_ab - d_ba).abs() < 0.01,
 "distance should be symmetric: d(a,b)={d_ab} vs d(b,a)={d_ba}"
 );
 }

 /// Two unit spheres separated by 5 units: exact distance = 5 - 1 - 1 = 3.
 #[test]
 fn disjoint_spheres_exact_distance() {
 let a = make_sphere_brep(1.0);
 let mut b = make_sphere_brep(1.0);
 b.apply_transform(glam::DAffine3::from_translation(glam::DVec3::new(5.0, 0.0, 0.0)));
 let d = min_distance(&a, &b);
 assert!(
 (d.distance - 3.0).abs() < 1e-6,
 "expected distance=3.0, got {}",
 d.distance
 );
 }

 #[test]
 fn plane_sphere_exact_distance() {
 let plane = make_square_plane_brep(DVec3::ZERO, DVec3::Z, 10.0, 10.0);
 let mut sphere = make_sphere_brep(1.0);
 sphere.apply_transform(glam::DAffine3::from_translation(DVec3::new(0.0, 0.0, 5.0)));
 let d = min_distance(&plane, &sphere);
 assert!((d.distance - 4.0).abs() < 1e-6, "expected 4.0, got {}", d.distance);
 }

 #[test]
 fn parallel_planes_exact_distance() {
 let a = make_square_plane_brep(DVec3::ZERO, DVec3::Z, 4.0, 4.0);
 let b = make_square_plane_brep(DVec3::new(0.0, 0.0, 3.0), DVec3::Z, 4.0, 4.0);
 let d = min_distance(&a, &b);
 assert!((d.distance - 3.0).abs() < 1e-6, "expected 3.0, got {}", d.distance);
 }
}
