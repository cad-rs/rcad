use super::*;

impl<'a> super::PaveFiller<'a> {
 pub(crate) fn intersect_torus_cylinder_faces(
 &mut self,
 f1: usize,
 f2: usize,
 torus: &ToroidalSurface,
 cylinder: &CylindricalSurface,
 ) {
 let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
 let s1 = Surface3::Torus(*torus);
 let s2 = Surface3::Cylinder(*cylinder);
 self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
 }

 pub(crate) fn intersect_torus_cone_faces(
 &mut self,
 f1: usize,
 f2: usize,
 torus: &ToroidalSurface,
 cone: &ConicalSurface,
 ) {
 let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
 let s1 = Surface3::Torus(*torus);
 let s2 = Surface3::Cone(*cone);
 self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
 }

 pub(crate) fn intersect_torus_torus_faces(
 &mut self,
 f1: usize,
 f2: usize,
 torus1: &ToroidalSurface,
 torus2: &ToroidalSurface,
 ) {
 let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
 let s1 = Surface3::Torus(*torus1);
 let s2 = Surface3::Torus(*torus2);
 self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
 }

 pub(crate) fn intersect_torus_sphere_faces(
 &mut self,
 f1: usize,
 f2: usize,
 torus: &ToroidalSurface,
 sphere: &SphericalSurface,
 ) {
 let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
 let s1 = Surface3::Torus(*torus);
 let s2 = Surface3::Sphere(*sphere);
 self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
 }

 pub(crate) fn trim_curve_to_faces(
 ds: &DS,
 curve: &Curve3,
 search_range: [f64; 2],
 f1: usize,
 f2: usize,
 ) -> Option<[f64; 2]> {
 use crate::medial_axis::point_in_polygon_2d;
 use rcad_kernel::projection::closest_point_on_surface;
 use std::f64::consts::TAU;

 const N: usize = 256;

 let face1 = &ds.faces[f1];
 let face2 = &ds.faces[f2];
 let uv_bnd1 = face1.uv_boundary.as_ref()?;
 let uv_bnd2 = face2.uv_boundary.as_ref()?;
 let s1 = &face1.surface;
 let s2 = &face2.surface;

 // UV from 3D point on a surface, normalising u �?[0, 2 ].
 let uv_on_surface = |surface: &Surface3, p: DVec3| -> DVec2 {
 match surface {
 Surface3::Cone(cone) => {
 let uv = cone.world_to_uv(p);
 DVec2::new(if uv.x < 0.0 { uv.x + TAU } else { uv.x }, uv.y)
 }
 Surface3::Sphere(sph) => sph.world_to_uv(p),
 Surface3::Cylinder(cyl) => {
 let x_ax = cyl.ref_dir.normalize();
 let y_ax = cyl.axis.cross(x_ax).normalize();
 let local = p - cyl.origin;
 let u = local.dot(y_ax).atan2(local.dot(x_ax));
 DVec2::new(if u < 0.0 { u + TAU } else { u }, local.dot(cyl.axis))
 }
 _ => {
 let proj = closest_point_on_surface(surface, p, 16);
 DVec2::new(proj.params.0, proj.params.1)
 }
 }
 };

 // True when the curve point at t is inside *both* faces' UV boundaries.
 // For planar faces the 3D point must actually lie on the plane (not just
 // project there), else an off-surface point would be a false positive.
 let point_in_both = |t: f64| -> bool {
 let pt = curve.point_at(t);
 for (sf, bnd) in &[(s1, uv_bnd1), (s2, uv_bnd2)] {
 if let Surface3::Plane(pl) = sf {
 if (pt - pl.origin).dot(pl.normal).abs() > TOLERANCE_COORD_SUB {
 return false;
 }
 }
 let uv = uv_on_surface(sf, pt);
 if !point_in_polygon_2d(uv, bnd) {
 return false;
 }
 }
 true
 };

 let [t0, t1] = search_range;
 let step = (t1 - t0) / N as f64;
 let mut seg_start: Option<(usize, f64)> = None;
 let mut segments: Vec<(usize, usize, f64, f64)> = Vec::new();

 for i in 0..=N {
 let t = t0 + step * i as f64;
 let inside = point_in_both(t);

 if inside {
 if seg_start.is_none() {
 seg_start = Some((i, t));
 }
 } else if let Some((si, st)) = seg_start.take() {
 if t - st > TOLERANCE_LINEAR_ULTRA_STRICT {
 segments.push((si, i, st, t));
 }
 }
 }
 if let Some((si, st)) = seg_start.take() {
 if t1 - st > TOLERANCE_LINEAR_ULTRA_STRICT {
 segments.push((si, N, st, t1));
 }
 }

 // Longest segment
 let (si, ei, rough_start, rough_end) = segments.into_iter().max_by(|a, b| {
 (a.3 - a.2)
 .partial_cmp(&(b.3 - b.2))
 .unwrap_or(std::cmp::Ordering::Equal)
 })?;

 //  € € binary-search refinement of both endpoints  € €
 // Start: between sample (si-1, outside) and sample (si, inside).
 let refined_start = if si > 0 {
 let t_out = t0 + step * (si - 1) as f64;
 let mut lo = t_out; // outside
 let mut hi = rough_start; // inside
 for _ in 0..48 {
 let mid = 0.5 * (lo + hi);
 if point_in_both(mid) {
 hi = mid;
 } else {
 lo = mid;
 }
 }
 hi
 } else {
 rough_start
 };

 // End: between sample (ei-1, inside) and sample (ei, outside).
 let refined_end = if ei < N {
 let t_in = t0 + step * (ei - 1) as f64;
 let t_out = t0 + step * ei as f64;
 let mut lo = t_in;  // inside
 let mut hi = t_out; // outside
 for _ in 0..48 {
 let mid = 0.5 * (lo + hi);
 if point_in_both(mid) {
 lo = mid;
 } else {
 hi = mid;
 }
 }
 lo
 } else {
 rough_end
 };

 if refined_end - refined_start > TOLERANCE_LINEAR_ULTRA_STRICT {
 Some([refined_start, refined_end])
 } else {
 None
 }
 }

 //  € € Plane �?Cone analytic face-face intersection  € € € € € € € € € € € € € € € € € € € € € € € € € €

 pub(crate) fn register_torus_intersection(
 &mut self,
 f1: usize,
 f2: usize,
 s1: &Surface3,
 s2: &Surface3,
 torus_is_f1: bool,
 ) {
 use inttools::intss::{intersect_surfaces_with_density_tol, SurfaceCurve};
 use inttools::pcurve_derive::polyline_pcurve_by_projection;

 let result = intersect_surfaces_with_density_tol(s1, s2, 48, self.ff_tol(f1, f2));
 if result.is_empty() {
 return;
 }

 for sir in &result.curves {
 match &sir.curve_3d {
 SurfaceCurve::Circle(circle) => {
 // Only split into half-circles for torus cylinder intersections where the
 // full circle spans 100% of cylinder U (triggers has_full_wrap fallback).
 // For other surface types the full circle is simpler and more robust.
 // Note: s1 is always Torus by calling convention, s2 is the other surface.
 if matches!(s2, Surface3::Cylinder(_)) {
 let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
 if torus_is_f1 { (Some(a.clone()), Some(b.clone())) }
 else { (Some(b.clone()), Some(a.clone())) }
 } else { (None, None) };

 let mut curve_indices = Vec::new();
 for (t0, t1) in [(0.0, std::f64::consts::PI), (std::f64::consts::PI, std::f64::consts::TAU)] {
 let pts = sample_circle_arc(circle, t0, t1, 16);
 if pts.len() < 2 { continue; }
 let v_start = self.ds.add_vertex(pts[0]);
 let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

 let curve_idx = self.ds.intersection_curves.len();
 eprintln!("[IC] CREATE ci={} f1={} f2={} sv={} ev={}", curve_idx, f1, f2, v_start, v_end);
 self.ds.intersection_curves.push(IntersectionCurve {
 curve: Curve3::Circle(*circle),
 polyline: vec![],
 start_vertex: v_start,
 end_vertex: v_end,
 t_range: [t0, t1],
 pcurve_on_a: pca.clone(),
 pcurve_on_b: pcb.clone(),
 geom_tol: crate::tolerance::TOLERANCE_ABS,
 pave_blocks: Vec::new(),
 curve_extra: crate::bopds::ds::CurveExtra::default(),
 });

 self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f1].face_info.vertices_in.insert(v_start);
 self.ds.faces[f1].face_info.vertices_in.insert(v_end);
 self.ds.faces[f2].face_info.vertices_in.insert(v_start);
 self.ds.faces[f2].face_info.vertices_in.insert(v_end);

 curve_indices.push(curve_idx);
 }

 if !curve_indices.is_empty() {
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF{
 f1, f2, curves: curve_indices,  points: vec![],
  tangent_faces: false,
  });
 }
 } else {
 let pts = sample_circle_arc(circle, 0.0, std::f64::consts::TAU, 32);
 if pts.len() < 2 { continue; }
 let v_start = self.ds.add_vertex(pts[0]);
 let v_end = self.ds.add_vertex(pts[pts.len() - 1]);
 let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
 if torus_is_f1 { (Some(a.clone()), Some(b.clone())) }
 else { (Some(b.clone()), Some(a.clone())) }
 } else { (None, None) };

 let curve_idx = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(IntersectionCurve {
 curve: Curve3::Circle(*circle),
 polyline: vec![],
 start_vertex: v_start,
 end_vertex: v_end,
 t_range: [0.0, std::f64::consts::TAU],
 pcurve_on_a: pca,
 pcurve_on_b: pcb,
 geom_tol: crate::tolerance::TOLERANCE_ABS,
 pave_blocks: Vec::new(),
 curve_extra: crate::bopds::ds::CurveExtra::default(),
 });

 self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f1].face_info.vertices_in.insert(v_start);
 self.ds.faces[f1].face_info.vertices_in.insert(v_end);
 self.ds.faces[f2].face_info.vertices_in.insert(v_start);
 self.ds.faces[f2].face_info.vertices_in.insert(v_end);

 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF{
 f1, f2, curves: vec![curve_idx],  points: vec![],
  tangent_faces: false,
  });
 }
 }
 SurfaceCurve::Polyline(pts) => {
 if pts.len() < 2 {
 continue;
 }
 let v_start = self.ds.add_vertex(pts[0]);
 let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

 let arc_len: f64 = pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
 let dir = (pts[pts.len() - 1] - pts[0]).normalize_or_zero();

 let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
 if torus_is_f1 {
 (Some(a.clone()), Some(b.clone()))
 } else {
 (Some(b.clone()), Some(a.clone()))
 }
 } else {
 (
 polyline_pcurve_by_projection(pts, s1),
 polyline_pcurve_by_projection(pts, s2),
 )
 };

 let curve_idx = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(IntersectionCurve {
 curve: Curve3::Line(Line3 {
 origin: pts[0],
 direction: if dir.length_squared() > 0.5 { dir } else { DVec3::X },
 }),
 polyline: pts.clone(),
 start_vertex: v_start,
 end_vertex: v_end,
 t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
 pcurve_on_a: pca,
 pcurve_on_b: pcb,
 geom_tol: crate::tolerance::TOLERANCE_ABS,
 pave_blocks: Vec::new(),
 curve_extra: crate::bopds::ds::CurveExtra::default(),
 });

 self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f1].face_info.vertices_in.insert(v_start);
 self.ds.faces[f1].face_info.vertices_in.insert(v_end);
 self.ds.faces[f2].face_info.vertices_in.insert(v_start);
 self.ds.faces[f2].face_info.vertices_in.insert(v_end);

 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF{
 f1,
 f2,
 curves: vec![curve_idx],
  points: vec![],
  tangent_faces: false,
  });
 }
 SurfaceCurve::Ellipse(ellipse) => {
 let pts = sample_circle_arc(
 &Circle3::new(ellipse.center, ellipse.normal, ellipse.major_radius,
 ),
 0.0,
 std::f64::consts::TAU,
 32,
 );
 if pts.len() < 2 {
 continue;
 }
 let v_start = self.ds.add_vertex(pts[0]);
 let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

 let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
 if torus_is_f1 {
 (Some(a.clone()), Some(b.clone()))
 } else {
 (Some(b.clone()), Some(a.clone()))
 }
 } else {
 (None, None)
 };

 let curve_idx = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(IntersectionCurve {
 curve: Curve3::Ellipse(*ellipse),
 polyline: vec![],
 start_vertex: v_start,
 end_vertex: v_end,
 t_range: [0.0, std::f64::consts::TAU],
 pcurve_on_a: pca,
 pcurve_on_b: pcb,
 geom_tol: crate::tolerance::TOLERANCE_ABS,
 pave_blocks: Vec::new(),
 curve_extra: crate::bopds::ds::CurveExtra::default(),
 });

 self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f1].face_info.vertices_in.insert(v_start);
 self.ds.faces[f1].face_info.vertices_in.insert(v_end);
 self.ds.faces[f2].face_info.vertices_in.insert(v_start);
 self.ds.faces[f2].face_info.vertices_in.insert(v_end);

 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF{
 f1,
 f2,
 curves: vec![curve_idx],
  points: vec![],
  tangent_faces: false,
  });
 }
 SurfaceCurve::Line(line) => {
 let pts = self.ds.face_boundary_points(f1);
 let pts2 = self.ds.face_boundary_points(f2);
 let bbox1_min = pts.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
 let bbox1_max = pts.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));
 let bbox2_min = pts2.iter().fold(DVec3::INFINITY, |a, &b| a.min(b));
 let bbox2_max = pts2.iter().fold(DVec3::NEG_INFINITY, |a, &b| a.max(b));

 let lo = bbox1_min.min(bbox2_min);
 let hi = bbox1_max.max(bbox2_max);
 let extent = (hi - lo).length() * 0.5 + 1.0;

 let p_start = line.origin + line.direction * (-extent);
 let p_end = line.origin + line.direction * extent;

 let v_start = self.ds.add_vertex(p_start);
 let v_end = self.ds.add_vertex(p_end);

 let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
 if torus_is_f1 {
 (Some(a.clone()), Some(b.clone()))
 } else {
 (Some(b.clone()), Some(a.clone()))
 }
 } else {
 (None, None)
 };

 let curve_idx = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(IntersectionCurve {
 curve: Curve3::Line(*line),
 polyline: vec![p_start, p_end],
 start_vertex: v_start,
 end_vertex: v_end,
 t_range: [-extent, extent],
 pcurve_on_a: pca,
 pcurve_on_b: pcb,
 geom_tol: crate::tolerance::TOLERANCE_ABS,
 pave_blocks: Vec::new(),
 curve_extra: crate::bopds::ds::CurveExtra::default(),
 });

 self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f1].face_info.vertices_in.insert(v_start);
 self.ds.faces[f1].face_info.vertices_in.insert(v_end);
 self.ds.faces[f2].face_info.vertices_in.insert(v_start);
 self.ds.faces[f2].face_info.vertices_in.insert(v_end);

 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF{
 f1,
 f2,
 curves: vec![curve_idx],
  points: vec![],
  tangent_faces: false,
  });
 }
 SurfaceCurve::BSplineCurve(b) => {
 // Sample the BSpline to produce a polyline for face splitting.
 use rcad_kernel::geom::CurveEval;
 let n_samples = 33_usize;
 let mut pts: Vec<DVec3> = Vec::with_capacity(n_samples);
 for i in 0..n_samples {
 let t = i as f64 / (n_samples - 1) as f64;
 pts.push(b.point_at(t));
 }
 if pts.len() < 2 {
 continue;
 }
 let v_start = self.ds.add_vertex(pts[0]);
 let v_end = self.ds.add_vertex(pts[pts.len() - 1]);

 let arc_len: f64 = pts.windows(2).map(|w| (w[1] - w[0]).length()).sum();
 let dir = (pts[pts.len() - 1] - pts[0]).normalize_or_zero();

 let (pca, pcb) = if let (Some(a), Some(b)) = (&sir.pcurve_on_a, &sir.pcurve_on_b) {
 if torus_is_f1 {
 (Some(a.clone()), Some(b.clone()))
 } else {
 (Some(b.clone()), Some(a.clone()))
 }
 } else {
 (
 polyline_pcurve_by_projection(&pts, s1),
 polyline_pcurve_by_projection(&pts, s2),
 )
 };

 let curve_idx = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(IntersectionCurve {
 curve: Curve3::BSpline((**b).clone()),
 polyline: pts.clone(),
 start_vertex: v_start,
 end_vertex: v_end,
 t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
 pcurve_on_a: pca,
 pcurve_on_b: pcb,
 geom_tol: crate::tolerance::TOLERANCE_ABS,
 pave_blocks: Vec::new(),
 curve_extra: crate::bopds::ds::CurveExtra::default(),
 });

 self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f1].face_info.vertices_in.insert(v_start);
 self.ds.faces[f1].face_info.vertices_in.insert(v_end);
 self.ds.faces[f2].face_info.vertices_in.insert(v_start);
 self.ds.faces[f2].face_info.vertices_in.insert(v_end);

 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF{
 f1,
 f2,
 curves: vec![curve_idx],
  points: vec![],
  tangent_faces: false,
  });
 }
 SurfaceCurve::Point(_) | SurfaceCurve::Parabola(_) | SurfaceCurve::Hyperbola(_) => {
 // Skip degenerate / unsupported curve types for now
 }
 }
 }
 }

 pub(crate) fn intersect_torus_plane_faces(
 &mut self,
 f1: usize,
 f2: usize,
 torus: &ToroidalSurface,
 plane: &Plane,
 ) {
 let torus_is_f1 = matches!(self.ds.faces[f1].surface, Surface3::Torus(_));
 let s1 = Surface3::Torus(*torus);
 let s2 = Surface3::Plane(*plane);
 self.register_torus_intersection(f1, f2, &s1, &s2, torus_is_f1);
 }

}
