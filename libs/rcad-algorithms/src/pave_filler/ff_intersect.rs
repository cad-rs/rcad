use super::*;
use crate::inttools::int_patch_type::IntPatchIType;

impl<'a> super::PaveFiller<'a> {
 pub(crate) fn perform_ff(&mut self) {
 // before performing intersection (ensures FaceInfo is current).
 let mut ff_face_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
 ff_face_set.extend(self.faces_of(ShapeOrigin::ShapeA));
 ff_face_set.extend(self.faces_of(ShapeOrigin::ShapeB));
 for &fi in &ff_face_set {
 self.ds.refine_face_info_on(fi);
 self.ds.refine_face_info_in(fi);
 }

 // OCCT PaveFiller_6.cxx: FillShrunkData + BVH pair iteration
 self.fill_shrunk_data(); // OCCT: FillShrunkData(FACE, FACE)
 let a_faces = self.faces_of(ShapeOrigin::ShapeA);
 let b_faces = self.faces_of(ShapeOrigin::ShapeB);

 if a_faces.is_empty() || b_faces.is_empty() {
 return;
 }

 if let (Some(ref fbvh)) = self.face_bvh {
 // DS-based BVH: indices are DS face indices, no a_rev/b_rev mapping.
 let candidates = crate::bvh::DsBvh::candidate_pairs(fbvh, fbvh);
 // Cross-origin filter + dedup.
 let mut processed_pairs = std::collections::HashSet::new();
 for &(fa, fb) in &candidates {
 if self.ds.faces[fa].origin == self.ds.faces[fb].origin { continue; }
 if !processed_pairs.insert((fa, fb)) { continue; }
 if self.ds.has_interf_ff(fa, fb) { continue; }
 if !self.should_skip_glued_face_pair(fa, fb) {
 self.intersect_face_face(fa, fb);
 }
 }
 } else {
 //  OCCT-aligned: BOPDS_Iterator cross-group face pair iteration.
 let a_fcount = self.ds.a_face_count;
 let mut fit = crate::bopds::ds::PairIterator::prepare_ab(a_fcount, self.ds.faces.len());
 while fit.more() {
 let pk = fit.value();
 let af = pk.i1; let bf = pk.i2;
 // OCCT: myDS->HasInterf(nF1, nF2)  ?skip if already interfered
 if self.ds.has_interf_ff(af, bf) { fit.next(); continue; }
 if !self.should_skip_glued_face_pair(af, bf) {
 self.intersect_face_face(af, bf);
 }
 fit.next();
 }
 }
 // OCCT-aligned: dedup FF interferences by pair (BOPDS_IndexRange).
 self.ds.dedup_ff_interferences();
 }

 pub(crate) fn should_skip_glued_face_pair(&self, f1: usize, f2: usize) -> bool {
 if !self.use_glue() {
 return false;
 }

 // Use pre-detected fully-glued faces if available
 if self.ds.is_fully_glued_face_pair(f1, f2) {
 return true;
 }

 let face1 = &self.ds.faces[f1];
 let face2 = &self.ds.faces[f2];
 if face1.origin == face2.origin {
 return false;
 }
 if !self.surfaces_glue_compatible(&face1.surface, &face2.surface) {
 return false;
 }

 let n1_len2 = face1.normal.length_squared();
 let n2_len2 = face2.normal.length_squared();
 if n1_len2 <= TOLERANCE_ABS || n2_len2 <= TOLERANCE_ABS {
 return false;
 }
 let n1 = face1.normal / n1_len2.sqrt();
 let n2 = face2.normal / n2_len2.sqrt();
 if n1.dot(n2) > -0.99 {
 return false;
 }

 self.boundaries_fully_overlap(f1, f2)
 }

 pub(crate) fn is_seam_edge(&self, edge_idx: usize, face_idx: usize) -> bool {
 let face = &self.ds.faces[face_idx];
 let edge = &self.ds.edges[edge_idx];

 match &face.surface {
 Surface3::Cylinder(cyl) => {
 // Cylinder seam edge: Line3 parallel to axis
 if let Curve3::Line(line) = &edge.curve {
 let dir = line.direction.normalize();
 let axis = cyl.axis.normalize();
 dir.dot(axis).abs() > 1.0 - TOLERANCE_ABS
 } else {
 false
 }
 }
 Surface3::Sphere(sph) => {
 // Sphere seam edge = great circle arc in meridian plane (U=0 boundary).
 // Checks mirror OCCT exactly:
 // (1) Curve is Geom_Circle ? Curve3::Circle
 // (2) |center - S.Location()| < Precision::Confusion() ? TOLERANCE_ABS_SQ
 // (3) |radius - S.Radius| < Precision::Confusion() ? TOLERANCE_ABS
 // (4) |circle_normal  ?sphere_axis| < Precision::Angular() ? perp check
 match &edge.curve {
 Curve3::Circle(c) => {
 (c.center - sph.center).length_squared() < TOLERANCE_ABS_SQ
 && (c.radius - sph.radius).abs() < TOLERANCE_ABS
 && c.normal.normalize().dot(sph.axis.normalize()).abs() < 1e-12
 }
 _ => false,
 }
 }
 Surface3::Torus(tor) => {
 // OCCT IsClosedFF: torus has TWO periodic boundaries.
 // U-seam: major circle, center = torus center, radius = major_radius,
 // normal  ?torus axis.
 // V-seam: minor circle, center on major circle, radius = minor_radius,
 // normal  ?torus axis.
 // All tolerances match OCCT Precision::Confusion/Angular.
 match &edge.curve {
 Curve3::Circle(c) => {
 let axis = tor.axis.normalize();
 let c_normal = c.normal.normalize();
 let center_dist = (c.center - tor.center).length();
 // U-seam: center at torus center, normal  ?axis, radius = major
 let is_u_seam = center_dist < TOLERANCE_ABS
 && c_normal.dot(axis).abs() > 1.0 - 1e-12
 && (c.radius - tor.major_radius).abs() < TOLERANCE_ABS;
 // V-seam: center on major circle, normal  ?axis, radius = minor
 let on_major = (center_dist - tor.major_radius).abs() < TOLERANCE_ABS;
 let is_v_seam = on_major
 && c_normal.dot(axis).abs() < 1e-12
 && (c.radius - tor.minor_radius).abs() < TOLERANCE_ABS;
 is_u_seam || is_v_seam
 }
 _ => false,
 }
 }
 _ => false,
 }
 }

 /// OCCT PaveFiller_6.cxx L393-479: seam edge shift
 /// OCCT PaveFiller_6.cxx L393-479: seam edge shift
 pub(crate) fn check_seam_edge_shift(&self, f1: usize, f2: usize) -> Option<SeamEdgeShift> {
 let s1 = &self.ds.faces[f1].surface;
 let s2 = &self.ds.faces[f2].surface;

 // Skip if both faces are Planes (seam edges only exist on periodic surfaces)
 if matches!(s1, Surface3::Plane(_)) && matches!(s2, Surface3::Plane(_)) {
 return None;
 }

 for &e1 in &self.ds.faces[f1].boundary_edges {
 let is_closed1 = self.is_seam_edge(e1, f1);
 for &e2 in &self.ds.faces[f2].boundary_edges {
 let is_closed2 = self.is_seam_edge(e2, f2);
 if !is_closed1 && !is_closed2 {
 continue;
 }

 // Look for EE interference between this edge pair
 for inf in &self.ds.interf_ee {
 if !((inf.e1 == e1 && inf.e2 == e2) || (inf.e1 == e2 && inf.e2 == e1)) {
 continue;
 }

 // Project the EE vertex point onto both edges' 3D curves
 // (OCCT: GeomAPI_ProjectPointOnCurve)
 let curve1 = &self.ds.edges[e1].curve;
 let curve2 = &self.ds.edges[e2].curve;
 let proj1 = closest_point_on_curve(curve1, inf.point, 64);
 let proj2 = closest_point_on_curve(curve2, inf.point, 64);

 let a_p1 = proj1.point;
 let a_p2 = proj2.point;
 let shift_dist = a_p1.distance(a_p2);

 // OCCT-aligned: the seam edge shift is a SMALL tolerance
 // correction, not a geometric transformation.  Verify both
 // projections are close to the EE vertex  ?if either is
 // far, the vertex is not near both edges and shifting would
 // be invalid (e.g. sphere center jumps by 1 unit).
 let vtx_pt = inf.point;
 let d1 = a_p1.distance(vtx_pt);
 let d2 = a_p2.distance(vtx_pt);
 // OCCT's shift is a sub-tolerance adjustment.  A projection
 // error exceeding 1e-4 means the vertex is not on this edge.
 let sanity_tol = TOLERANCE_ABS * 1000.0;
 if d1 > sanity_tol || d2 > sanity_tol {
 continue;
 }

 // Check if the shift exceeds vertex tolerance
 let vtx_tol = self.ds.vertices[inf.new_vertex].geom_tol;
 if shift_dist > vtx_tol {
 // OCCT: shift the face with the closed/seam edge
 let shift_vector = if is_closed1 {
 a_p2 - a_p1 // Shift f1: move aP1 toward aP2
 } else {
 a_p1 - a_p2 // Shift f2: move aP2 toward aP1
 };

 return Some(SeamEdgeShift {
 shift_vector,
 shift_value: shift_dist,
 shifted_face: if is_closed1 { 1 } else { 2 },
 });
 }
 }
 }
 }
 None
 }

 /// OCCT PaveFiller_6: reverse seam edge shift
 /// OCCT PaveFiller_6: reverse seam edge shift
 pub(crate) fn reverse_seam_edge_shift(&mut self, f1: usize, f2: usize, shift: &SeamEdgeShift) {
 let inv_vec = if shift.shifted_face == 1 {
 -shift.shift_vector
 } else {
 shift.shift_vector
 };

 // Collect curve indices from the FaceFace interference for this pair
 let mut curve_indices: Vec<usize> = Vec::new();
 for inf in &self.ds.interf_ff {
 if (inf.f1 == f1 && inf.f2 == f2) || (inf.f1 == f2 && inf.f2 == f1) {
 curve_indices = inf.curves.clone();
 break;
 }
 }

 // Reverse shift on each curve
 for &ci in &curve_indices {
 if ci >= self.ds.intersection_curves.len() {
 continue;
 }
 let ic = &mut self.ds.intersection_curves[ci];

 // Translate 3D curve back by inverse shift
 ic.curve = translate_curve3(&ic.curve, inv_vec);

 // Translate polyline points if any
 for p in &mut ic.polyline {
 *p += inv_vec;
 }

 // Translate vertex positions back
 let sv = ic.start_vertex;
 let ev = ic.end_vertex;
 if sv < self.ds.vertices.len() {
 self.ds.vertices[sv].point += inv_vec;
 }
 if ev < self.ds.vertices.len() {
 self.ds.vertices[ev].point += inv_vec;
 }
 }
 }

 /// OCCT: dispatch FF intersection by surface type
 /// OCCT L344-608: IntTools_FaceFace::Perform — face-face intersection.
/// Dispatches by surface type with bReverse sorting (OCCT SortTypes/IndexType),
/// then runs intersection, MakeCurve, ComputeTolReached3d, PrepareLines3D,
/// and point/curve registration.
pub(crate) fn intersect_face_face(&mut self, f1: usize, f2: usize) {
 // = =  Seam Edge Shift (OCCT PaveFiller_6.cxx L393-479) = = = = = = = = = = = = = = 
 let shift_info = self.check_seam_edge_shift(f1, f2);
 let old_shift_tol = self.seam_shift_tol;
 if let Some(ref info) = shift_info {
 self.seam_shift_tol = info.shift_value;
 }

 let s1_orig = self.ds.faces[f1].surface.clone();
 let s2_orig = self.ds.faces[f2].surface.clone();

 // Apply seam edge shift to surface clones if needed
 let s1 = match &shift_info {
 Some(info) if info.shifted_face == 1 => {
 apply_shift_to_surface(&s1_orig, info.shift_vector)
 }
 _ => s1_orig,
 };
 let s2 = match &shift_info {
 Some(info) if info.shifted_face == 2 => {
 apply_shift_to_surface(&s2_orig, info.shift_vector)
 }
 _ => s2_orig,
 };

 // OCCT L351-375: SortTypes — canonical surface ordering.
 // Swap f1/f2 so the higher-type surface is always "face A".
 let type_idx1 = Self::surface_type_index(&s1);
 let type_idx2 = Self::surface_type_index(&s2);
 let b_reverse = type_idx1 < type_idx2;
 // OCCT L354: if bReverse, swap face refs so myFace1 gets the higher type.
 let (f1, f2, s_a, s_b) = if b_reverse { (f2, f1, &s2, &s1) } else { (f1, f2, &s1, &s2) };

 // OCCT L384-393: tolerance setup
 // OCCT: myTolF1 = BRep_Tool::Tolerance(myFace1) + aFuzz, etc.
 // rcad: tolerance handled by PaveFiller's fuzzy_tolerance and face geom_tols.

 // OCCT L395-401: isFace1Quad/isFace2Quad — skip; rcad uses IntPatchIntersection
 // which dispatches by quad type internally.

 // ── OCCT L404-434: Plane-Plane fast path (PerformPlanes) ──
 if matches!(s_a, rcad_kernel::geom::Surface3::Plane(_))
   && matches!(s_b, rcad_kernel::geom::Surface3::Plane(_))
 {
 self.perform_plane_plane(f1, f2);
 if let Some(ref info) = shift_info { self.reverse_seam_edge_shift(f1, f2, info); }
 self.seam_shift_tol = old_shift_tol;
 return;
 }

 // OCCT L436-438: myLConstruct.Load(dom1, dom2, myHS1, myHS2)
 let mut lconstruct = crate::inttools::int_patch_line_constructor::GeomIntLineConstructor::new();
 lconstruct.load(f1, f2);

 // OCCT-aligned: IntPatch_Intersection: generic surface-surface intersection.
 let mut int_patch = crate::inttools::int_patch_intersection::IntPatchIntersection::new();
 int_patch.perform(s_a, s_b, self.fuzzy_tolerance, self.fuzzy_tolerance);
 if int_patch.tangent_faces() {
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
 f1, f2, curves: Vec::new(), points: Vec::new(), tangent_faces: true,
 });
 if let Some(ref info) = shift_info { self.reverse_seam_edge_shift(f1, f2, info); }
 self.seam_shift_tol = old_shift_tol;
 return;
 }

 // OCCT-aligned: PutPointsOnLine (IntPatch_Intersection.cxx L268-312).
 // Projects intersection points onto each analytic line to create
 // boundary-crossing vertices.  These vertices split the line into
 // valid intervals for MakeCurve/TreatCircle.
 for li in 0..int_patch.nb_lines() {
 self.put_points_on_line(f1, f2, int_patch.line_mut(li));
 }

 // OCCT-aligned: MakeCurve (IntTools_FaceFace.cxx L695-1846) for each IntPatch line.
 // Returns a Vec of IntersectionCurve — one per valid part from the
 // LineConstructor (OCCT supports aNbParts > 1, e.g. multi-segment clipping).
 let mut ff_curve_indices: Vec<usize> = Vec::new();
 for i in 0..int_patch.nb_lines() {
 let ics = self.make_intersection_curve(f1, f2, int_patch.line(i));
 for ic in ics {
   let ci = self.ds.intersection_curves.len();
   let mut adjusted_ic = ic;
   // OCCT L558-567: if reversed, swap pcurves (first ↔ second).
   if b_reverse {
   std::mem::swap(&mut adjusted_ic.pcurve_on_a, &mut adjusted_ic.pcurve_on_b);
   }
   self.ds.intersection_curves.push(adjusted_ic);
   ff_curve_indices.push(ci);
 }
 }

 // OCCT L576-608: points — filter by isPointInOnFace, append to myPnts.
 let mut ff_point_indices: Vec<usize> = Vec::new();
 for pi in 0..int_patch.nb_points() {
 let pt = int_patch.point(pi);
 let (uv_a, uv_b, f_a, f_b) = if b_reverse {
   (glam::DVec2::new(pt.u2, pt.v2), glam::DVec2::new(pt.u1, pt.v1), f2, f1)
 } else {
   (glam::DVec2::new(pt.u1, pt.v1), glam::DVec2::new(pt.u2, pt.v2), f1, f2)
 };
 if !self.context.is_point_in_on_face(self.ds, f_a, uv_a) { continue; }
 if !self.context.is_point_in_on_face(self.ds, f_b, uv_b) { continue; }
 let vi = self.ds.add_vertex(pt.p1);
 self.ds.vertices[vi].geom_tol = self.ds.vertices[vi].geom_tol.max(pt.tolerance);
 ff_point_indices.push(vi);
 }
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
 f1, f2, curves: ff_curve_indices, points: ff_point_indices, tangent_faces: false,
 });

 // = =  Reverse Seam Edge Shift (OCCT ApplyTrsf L560) = = = = = = = = = = = = = = 
 if let Some(ref info) = shift_info {
 self.reverse_seam_edge_shift(f1, f2, info);
 }
 // = =  Restore seam shift tol = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 self.seam_shift_tol = old_shift_tol;

 // OCCT-aligned: ComputeTolReached3d + PrepareLines3D.
 if let Some(ff_curves) = self.find_face_face_curve_indices(f1, f2) {
 let t_a = self.ff_tol(f1, f1);
 let t_b = self.ff_tol(f2, f2);
 for &ci in &ff_curves {
 let (curve, pca, pcb, tr, current_tol) = {
   let ic = &self.ds.intersection_curves[ci];
   (ic.curve.clone(),
    ic.pcurve_on_a.clone(), ic.pcurve_on_b.clone(),
    ic.t_range, ic.geom_tol)
 };
 let (new_tol, tang_tol) = inttools::pcurve_derive::compute_intersection_curve_tolerance(
 &curve, pca.as_ref(), pcb.as_ref(),
 &self.ds.faces[f1].surface, &self.ds.faces[f2].surface, tr,
 t_a, t_b, current_tol,
 );
 let ic = &mut self.ds.intersection_curves[ci];
 ic.geom_tol = ic.geom_tol.max(new_tol);
 ic.curve_extra.tangential_tol = ic.curve_extra.tangential_tol.max(tang_tol);
 }
 // PrepareLines3D  ?split closed curves
 let n_curves_before_split = self.ds.intersection_curves.len();
 inttools::pcurve_derive::prepare_lines_3d(&mut self.ds.intersection_curves);
 // OCCT-aligned: After PrepareLines3D splits closed curves, the split
 // segments are added to the same FF interference entry.  Update the
 // FF entry's curve list to include any newly created curve indices.
 if n_curves_before_split != self.ds.intersection_curves.len() {
  if let Some(ff_entry) = self.ds.interf_ff.last_mut() {
   for new_ci in n_curves_before_split..self.ds.intersection_curves.len() {
    ff_entry.curves.push(new_ci);
   }
  }
 }
 //  OCCT-aligned: After PrepareLines3D splits closed curves, new curve endpoints
 // must be updated to the split points. OCCT's BRepBuilderAPI_MakeEdge auto-sets
 // endpoints when creating edges. rcad's IntersectionCurve requires explicit update:
 // for start==end but non-full-period t_range (i.e. split half-circle), compute
 // correct endpoint positions via point_at and create new DS vertices.
 for ci in 0..self.ds.intersection_curves.len() {
 let needs_fix = {
   let ic = &self.ds.intersection_curves[ci];
   let half_circle = match &ic.curve {
   rcad_kernel::geom::Curve3::Circle(_) | rcad_kernel::geom::Curve3::Ellipse(_) => {
   (ic.t_range[1] - ic.t_range[0] - std::f64::consts::TAU).abs() >= TOLERANCE_ANG
   }
   _ => false,
   };
   half_circle && ic.start_vertex == ic.end_vertex
 };
 if needs_fix {
  if std::env::var("RCAD_DBG_FF").is_ok() {
   eprintln!("[DBG_FF] needs_fix: ci={} t=[{:.4},{:.4}] sv={} ev={}", ci,
   self.ds.intersection_curves[ci].t_range[0], self.ds.intersection_curves[ci].t_range[1],
   self.ds.intersection_curves[ci].start_vertex, self.ds.intersection_curves[ci].end_vertex);
  }
  let t0 = self.ds.intersection_curves[ci].t_range[0];
  let t1 = self.ds.intersection_curves[ci].t_range[1];
  let p_start = self.ds.intersection_curves[ci].curve.point_at(t0);
  let p_end = self.ds.intersection_curves[ci].curve.point_at(t1);
  let v_start = self.ds.vertices.len();
  self.ds.vertices.push(DSVertex { point: p_start, geom_tol: TOLERANCE_ABS, origin: None, is_internal: false, location: 0 });
  let v_end = self.ds.vertices.len();
  self.ds.vertices.push(DSVertex { point: p_end, geom_tol: TOLERANCE_ABS, origin: None, is_internal: false, location: 0 });
  self.ds.intersection_curves[ci].start_vertex = v_start;
  self.ds.intersection_curves[ci].end_vertex = v_end;
 }
 }
 // OCCT-aligned: PreparePostTreatFF (PaveFiller_6.cxx L3642-3668).
 let post_ff_curves = self.find_face_face_curve_indices(f1, f2)
 .unwrap_or_default();
 self.ds.faces[f1].face_info.curves_sc.extend(&post_ff_curves);
 self.ds.faces[f2].face_info.curves_sc.extend(&post_ff_curves);
 for &ci in &post_ff_curves {
 if ci < self.ds.intersection_curves.len() {
   let ic = &self.ds.intersection_curves[ci];
   if ic.start_vertex < self.ds.vertices.len() {
   self.ds.faces[f1].face_info.vertices_in.insert(ic.start_vertex);
   self.ds.faces[f2].face_info.vertices_in.insert(ic.start_vertex);
   }
   if ic.end_vertex < self.ds.vertices.len() {
   self.ds.faces[f1].face_info.vertices_in.insert(ic.end_vertex);
   self.ds.faces[f2].face_info.vertices_in.insert(ic.end_vertex);
   }
 }
 }
  } // if let Some(ff_curves)
} // fn intersect_face_face

/// OCCT-aligned: IndexType (IntTools_FaceFace.cxx L2844-2870).
/// Maps Surface3 variant to an integer index for canonical ordering.
/// Lower-typed surface is "simpler" (Plane<Cylinder<Cone<Sphere<Torus).
fn surface_type_index(surf: &rcad_kernel::geom::Surface3) -> i32 {
 match surf {
 rcad_kernel::geom::Surface3::Plane(_) => 0,
 rcad_kernel::geom::Surface3::Cylinder(_) => 1,
 rcad_kernel::geom::Surface3::Cone(_) => 2,
 rcad_kernel::geom::Surface3::Sphere(_) => 3,
 rcad_kernel::geom::Surface3::Torus(_) => 4,
 _ => 11,
 }
}

/// OCCT L2426-2560: PerformPlanes — plane-plane intersection fast path.
fn perform_plane_plane(&mut self, f1: usize, f2: usize) {
 use rcad_kernel::geom::{Curve3, Surface3};
 let pln1 = match &self.ds.faces[f1].surface { Surface3::Plane(p) => p, _ => return };
 let pln2 = match &self.ds.faces[f2].surface { Surface3::Plane(p) => p, _ => return };
 let mut geo = crate::inttools::int_ana_quad_quad_geo::QuadQuadGeo::new();
 let q1 = crate::inttools::int_surf_quadric::Quadric::from_plane(pln1);
 let q2 = crate::inttools::int_surf_quadric::Quadric::from_plane(pln2);
 let (Some(ref q1), Some(ref q2)) = (q1, q2) else { return };
 geo.perform_plane_plane(q1, q2, 1e-8);
 if !geo.is_done() { return; }
 use crate::inttools::int_ana_quad_quad_geo::AnaResultType;
 if let AnaResultType::Same = geo.type_inter() {
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
   f1, f2, curves: Vec::new(), points: Vec::new(), tangent_faces: true,
 });
 return;
 }
 if matches!(geo.type_inter(), AnaResultType::Empty | AnaResultType::NoResult) { return; }
 let line3 = geo.line(1);
 let line3d = Curve3::Line(*line3);
 let pcurve1 = crate::inttools::pcurve_derive::line_pcurve_on_plane(line3, pln1);
 let pcurve2 = crate::inttools::pcurve_derive::line_pcurve_on_plane(line3, pln2);
 let uv1 = self.context.uv_bounds(self.ds, f1);
 let uv2 = self.context.uv_bounds(self.ds, f2);
 let tol = self.ds.faces[f1].geom_tol.max(self.ds.faces[f2].geom_tol);
 let p1 = classify_lin2d(&pcurve1, uv1, tol);
 let p2 = classify_lin2d(&pcurve2, uv2, tol);
 let (Some([p11, p12]), Some([p21, p22])) = (p1, p2) else { return };
 if p21 >= p12 || p22 <= p11 { return; }
 let pmin = p11.max(p21);
 let pmax = p12.min(p22);
 if pmax - pmin <= tol { return; }
 let t_range = [pmin, pmax];
 let mut curve_extra = crate::bopds::ds::CurveExtra::default();
 curve_extra.tangential_tol = tol;
 let ic = crate::bopds::ds::IntersectionCurve {
 curve: line3d, polyline: Vec::new(),
 start_vertex: usize::MAX, end_vertex: usize::MAX,
 t_range,
 pcurve_on_a: Some(pcurve1), pcurve_on_b: Some(pcurve2),
 geom_tol: tol, pave_blocks: Vec::new(), curve_extra,
 };
 let ci = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(ic);
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
 f1, f2, curves: vec![ci], points: Vec::new(), tangent_faces: false,
 });
 self.ds.faces[f1].face_info.curves_sc.insert(ci);
 self.ds.faces[f2].face_info.curves_sc.insert(ci);
}

/// OCCT L2574-2640: ClassifyLin2d — clip a 2D line to a UV rectangle.
/// Returns parameter range [p1,p2] where line passes through [xmin,xmax]×[ymin,ymax],
/// or None if it misses.  Exported as pub for testing.
pub fn classify_lin2d(pc: &rcad_kernel::geom::Curve2d, uv: [f64; 4], tol: f64) -> Option<[f64; 2]> {
 use rcad_kernel::geom::{Curve2d, Curve2dEval};
 let (xmin, xmax, ymin, ymax) = (uv[0], uv[1], uv[2], uv[3]);
 let (A, B, C) = match pc {
 Curve2d::Line(l) => {
   (l.direction.y, -l.direction.x, -(l.direction.y * l.origin.x - l.direction.x * l.origin.y))
 }
 _ => return None,
 };
 fn inter(a: f64, b: f64, tl: f64) -> bool { (a < -tl && b > tl) || (a > tl && b < -tl) }
 fn coinc(a: f64, b: f64, tl: f64) -> bool { a.abs() <= tl && b.abs() <= tl }
 let mut par: Vec<f64> = Vec::with_capacity(2);
 // edge x=xmin, y∈[ymin,ymax]
 let d1 = A * xmin + B * ymin + C;
 let d2 = A * xmin + B * ymax + C;
 if inter(d1, d2, tol) && B.abs() > 1e-15 {
 let y = -(C + A * xmin) / B;
 if y >= ymin - tol && y <= ymax + tol { par.push(line2d_param(pc, glam::DVec2::new(xmin, y))); }
 } else if coinc(d1, d2, tol) {
 par.push(line2d_param(pc, glam::DVec2::new(xmin, ymin)));
 par.push(line2d_param(pc, glam::DVec2::new(xmin, ymax)));
 }
 if par.len() >= 2 { return Some([par[0].min(par[1]), par[0].max(par[1])]); }
 // edge y=ymax, x∈[xmin,xmax]
 let d1 = A * xmin + B * ymax + C;
 let d2 = A * xmax + B * ymax + C;
 if inter(d1, d2, tol) && A.abs() > 1e-15 {
 let x = -(C + B * ymax) / A;
 if x >= xmin - tol && x <= xmax + tol { par.push(line2d_param(pc, glam::DVec2::new(x, ymax))); }
 } else if coinc(d1, d2, tol) && par.is_empty() {
 par.push(line2d_param(pc, glam::DVec2::new(xmin, ymax)));
 par.push(line2d_param(pc, glam::DVec2::new(xmax, ymax)));
 }
 if par.len() >= 2 { return Some([par[0].min(par[1]), par[0].max(par[1])]); }
 // edge x=xmax, y∈[ymin,ymax]
 let d1 = A * xmax + B * ymax + C;
 let d2 = A * xmax + B * ymin + C;
 if inter(d1, d2, tol) && B.abs() > 1e-15 {
 let y = -(C + A * xmax) / B;
 if y >= ymin - tol && y <= ymax + tol { par.push(line2d_param(pc, glam::DVec2::new(xmax, y))); }
 } else if coinc(d1, d2, tol) && par.is_empty() {
 par.push(line2d_param(pc, glam::DVec2::new(xmax, ymin)));
 par.push(line2d_param(pc, glam::DVec2::new(xmax, ymax)));
 }
 if par.len() >= 2 { return Some([par[0].min(par[1]), par[0].max(par[1])]); }
 // edge y=ymin, x∈[xmin,xmax]
 let d1 = A * xmax + B * ymin + C;
 let d2 = A * xmin + B * ymin + C;
 if inter(d1, d2, tol) && A.abs() > 1e-15 {
 let x = -(C + B * ymin) / A;
 if x >= xmin - tol && x <= xmax + tol { par.push(line2d_param(pc, glam::DVec2::new(x, ymin))); }
 } else if coinc(d1, d2, tol) && par.is_empty() {
 par.push(line2d_param(pc, glam::DVec2::new(xmin, ymin)));
 par.push(line2d_param(pc, glam::DVec2::new(xmax, ymin)));
 }
 if par.len() >= 2 { Some([par[0].min(par[1]), par[0].max(par[1])]) } else { None }
}

/// Helper: parameter of a 2D point on a Line2d.
fn line2d_param(pc: &rcad_kernel::geom::Curve2d, p: glam::DVec2) -> f64 {
 let l = match pc { rcad_kernel::geom::Curve2d::Line(l) => l, _ => return 0.0 };
 (p - l.origin).dot(l.direction)
}

 /// OCCT-aligned: MakeCurve (IntTools_FaceFace.cxx L695-1846).
/// Dispatches by line type (OCCT switch on IntPatch_IType):
///   - Walking:        approximate BSpline from marching points (L1097)
///   - Line/Parabola/Hyperbola: LineConstructor parts + per-part handling (L815-898)
///   - Circle/Ellipse:  TreatCircle-equivalent with 0-crossing splitting (L904-1095)
///   - (Restriction:    handled upstream in IntPatch_Intersection)
/// Returns one or more IntersectionCurve per valid part.
pub(crate) fn make_intersection_curve(
  &mut self, f1: usize, f2: usize,
  line: &crate::inttools::int_patch_line::IntPatchLine,
) -> Vec<crate::bopds::ds::IntersectionCurve> {
 use rcad_kernel::geom::Curve2dEval;
 use std::f64::consts::TAU;

 // OCCT L1097-1846: IntPatch_Walking → approximate BSpline from marching points.
 if line.is_wline() {
 return self.make_walking_curve(f1, f2, line);
 }

 let curve = line.curve.clone();
 let orig_t_range = line.t_range;
 let geom_tol = line.tolerance;
 let typl = line.line_type;

 // OCCT L815-1095: switch on line type.
 match typl {
 IntPatchIType::Line | IntPatchIType::Parabola | IntPatchIType::Hyperbola =>
   self.make_analytic_nonperiodic_curve(f1, f2, &curve, orig_t_range, typl, geom_tol, line.tang_tolerance),
 IntPatchIType::Circle | IntPatchIType::Ellipse =>
   self.make_analytic_periodic_curve(f1, f2, &curve, orig_t_range, typl, geom_tol, line.tang_tolerance, &line.vertices),
 _ => Vec::new(),
 }
}

/// OCCT L1097-1846: MakeCurve for IntPatch_Walking.
/// Approximates a BSpline3 from marching points, builds pcurves from
/// marching UV data.
fn make_walking_curve(
  &mut self, _f1: usize, _f2: usize,
  line: &crate::inttools::int_patch_line::IntPatchLine,
) -> Vec<crate::bopds::ds::IntersectionCurve> {
 let n = line.nb_points();
 if n < 2 { return Vec::new(); }

 let p3d_pts: Vec<glam::DVec3> = (0..n).map(|i| line.point(i).p3d).collect();
 let polyline = p3d_pts.clone();

 if let Some(bs_curve3) = crate::inttools::intss::polyline_to_bspline(&p3d_pts, 1e-4) {
 let t_range_bs = bs_curve3.default_domain();
 let bs = match &bs_curve3 {
 rcad_kernel::geom::Curve3::BSpline(b) => b.clone(),
 _ => {
   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = line.tang_tolerance;
   return vec![crate::bopds::ds::IntersectionCurve {
   curve: bs_curve3.clone(), polyline, start_vertex: usize::MAX, end_vertex: usize::MAX,
   t_range: t_range_bs, pcurve_on_a: line.pcurve1.clone(), pcurve_on_b: line.pcurve2.clone(),
   geom_tol: line.tolerance.max(CONFUSION),
   pave_blocks: Vec::new(), curve_extra,
   }];
 }
 };

 // Build pcurves from marching UV data.
 let mut pcurve_on_a = line.pcurve1.clone();
 let mut pcurve_on_b = line.pcurve2.clone();
 if pcurve_on_a.is_none() && line.point(0).u1.is_finite() {
  let uv_pts: Vec<glam::DVec2> = (0..n).map(|i| glam::DVec2::new(line.point(i).u1, line.point(i).v1)).collect();
  if let Ok(bs2d) = rcad_kernel::fit::interpolate_points_2d(&uv_pts) {
  pcurve_on_a = Some(rcad_kernel::geom::Curve2d::BSpline(bs2d));
  }
 }
 if pcurve_on_b.is_none() && line.point(0).u2.is_finite() {
  let uv_pts: Vec<glam::DVec2> = (0..n).map(|i| glam::DVec2::new(line.point(i).u2, line.point(i).v2)).collect();
  if let Ok(bs2d) = rcad_kernel::fit::interpolate_points_2d(&uv_pts) {
  pcurve_on_b = Some(rcad_kernel::geom::Curve2d::BSpline(bs2d));
  }
 }
 let mut curve_extra = crate::bopds::ds::CurveExtra::default();
 curve_extra.tangential_tol = line.tang_tolerance;
 return vec![crate::bopds::ds::IntersectionCurve {
  curve: rcad_kernel::geom::Curve3::BSpline(bs),
  polyline,
  start_vertex: usize::MAX, end_vertex: usize::MAX,
  t_range: t_range_bs,
  pcurve_on_a, pcurve_on_b,
  geom_tol: line.tolerance.max(CONFUSION),
  pave_blocks: Vec::new(),
  curve_extra,
 }];
 }
 Vec::new()
}

/// OCCT L815-898: MakeCurve for Line, Parabola, Hyperbola.
/// - Creates analytic curve from IntPatch_GLine.
/// - Calls LineConstructor to get valid parameter parts (OCCT NbParts/Part).
/// - For each part:
///     both bounds finite  → trimmed 3D curve + BuildPCurves + endpoint vertices
///     one/both infinite   → test reference point on face domains → keep or reject
/// - rcad note: IntPatchLine has no vertex data, so LineConstructor returns
///   a single part with the original t_range (always infinite for lines).
fn make_analytic_nonperiodic_curve(
  &mut self, f1: usize, f2: usize,
  curve: &Curve3, orig_t_range: [f64; 2], typl: IntPatchIType,
  geom_tol: f64, tang_tolerance: f64,
) -> Vec<crate::bopds::ds::IntersectionCurve> {
 use rcad_kernel::geom::Curve2dEval;
 use std::f64::consts::TAU;

 // OCCT L815-826: create analytic 3D curve from the GLine.
 // rcad: curve is already the correct analytic type in IntPatchLine.

 // OCCT L828-840: LineConstructor.Perform(L).
 // LineConstructor iterates over vertex intervals; with nbvtx=0 (rcad)
 // it returns a single part with the full parameter range.
 let parts = self.line_constructor_parts(curve, orig_t_range, typl);
 if parts.is_empty() {
 return Vec::new();
 }

 let mut result = Vec::with_capacity(parts.len());

 // OCCT L842-898: per-part loop.
 for &[fprm, lprm] in &parts {
 let b_finite = fprm.is_finite() && lprm.is_finite() && lprm > fprm + 1e-12;

 if b_finite {
   // ── Both bounds finite: trimmed curve + pcurves + vertices ──
   // OCCT L835-870: Geom_TrimmedCurve + BuildPCurves + Geom2d_TrimmedCurve.
   let ic_t_range = [fprm, lprm];

   // BuildPCurves on the trimmed range (OCCT L846-863).
   let pca = self.compute_pcurve_on_surface(curve, f1);
   let pcb = self.compute_pcurve_on_surface(curve, f2);

   // OCCT L868-870: SetCurve(TrimmedCurve) + SetFirstCurve2d(TrimmedCurve(...)).
   // rcad: analytic curve retained, t_range provides effective trimming.
   let trimmed_pca = pca.as_ref().map(|pc| pc.clone());
   let trimmed_pcb = pcb.as_ref().map(|pc| pc.clone());

   // OCCT: for Parabola, CurveTolerance(aCT3D, myTol)
   let ic_geom_tol = if typl == IntPatchIType::Parabola {
   geom_tol.max(crate::tolerance::TOLERANCE_ABS)
   } else { geom_tol.max(crate::tolerance::TOLERANCE_ABS) };

   // Create endpoint vertices (OCCT L872-890).
   let (sv, ev) = {
   let p_start = curve.point_at(fprm);
   let p_end = curve.point_at(lprm);
   if p_start.is_finite() && p_end.is_finite() {
     let sv = self.ds.vertices.len();
     self.ds.vertices.push(DSVertex { point: p_start, geom_tol: ic_geom_tol, origin: None, is_internal: false, location: 0 });
     let ev = self.ds.vertices.len();
     self.ds.vertices.push(DSVertex { point: p_end, geom_tol: ic_geom_tol, origin: None, is_internal: false, location: 0 });
     (sv, ev)
   } else { (usize::MAX, usize::MAX) }
   };

   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
   curve: curve.clone(),
   polyline: Vec::new(),
   start_vertex: sv,
   end_vertex: ev,
   t_range: ic_t_range,
   pcurve_on_a: trimmed_pca,
   pcurve_on_b: trimmed_pcb,
   geom_tol: ic_geom_tol,
   pave_blocks: Vec::new(),
   curve_extra,
   });

 } else {
   // ── One/both bounds infinite: test reference point ──
   // OCCT L850-895: test-point approach.
   // dT = 100.0; surface-type exceptions for extrusion/offset/revolution.
   let dT = 100.0;
   let test_t = if !fprm.is_finite() && lprm.is_finite() {
   // bFNIt && !bLPIt: only lower bound infinite
   lprm - dT
   } else if fprm.is_finite() && !lprm.is_finite() {
   // !bFNIt && bLPIt: only upper bound infinite
   fprm + dT
   } else {
   // bFNIt && bLPIt: both infinite → IntermediatePoint(-dT, dT)
   0.0
   };

   let p3d = curve.point_at(test_t);
   if !p3d.is_finite() { continue; }

   // OCCT L865-880: classify test point on both face domains.
   let uv1 = self.context.proj_ps(self.ds, f1, p3d);
   let uv2 = self.context.proj_ps(self.ds, f2, p3d);
   let in1 = uv1.map_or(false, |(uv, _, _)| self.context.is_point_in_on_face(self.ds, f1, uv));
   let in2 = uv2.map_or(false, |(uv, _, _)| self.context.is_point_in_on_face(self.ds, f2, uv));
   if !in1 || !in2 { continue; }

   // OCCT L882-895: if both inside, append curve WITHOUT pcurves
   // (Geom2d_BSplineCurve H1; SeqOfCurve.Append(IntTools_Curve(newc, H1, H1))).
   // rcad: compute pcurves for downstream use.
   let pca = self.compute_pcurve_on_surface(curve, f1);
   let pcb = self.compute_pcurve_on_surface(curve, f2);

   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
   curve: curve.clone(),
   polyline: Vec::new(),
   start_vertex: usize::MAX,
   end_vertex: usize::MAX,
   t_range: orig_t_range,
   pcurve_on_a: pca,
   pcurve_on_b: pcb,
   geom_tol: geom_tol.max(crate::tolerance::TOLERANCE_ABS),
   pave_blocks: Vec::new(),
   curve_extra,
   });
 }
 }
 result
}

/// OCCT L904-1095: MakeCurve for Circle, Ellipse.
/// - Creates analytic curve from IntPatch_GLine.
/// - Calls TreatCircle to sort vertices, split at 0, build candidate intervals.
/// - For each candidate interval:
///     not full-period → trimmed curve + BuildPCurves(with UV bounds) + vertices
///     full-period (aNbParts=1) → test 18 points around circle → keep/reject
/// - rcad: with 0 vertices, TreatCircle returns fallback full [0, 2π] interval.
fn make_analytic_periodic_curve(
  &mut self, f1: usize, f2: usize,
  curve: &Curve3, orig_t_range: [f64; 2], typl: IntPatchIType,
  geom_tol: f64, tang_tolerance: f64,
  vertices: &[crate::inttools::int_patch_line::IntPatchVertex],
) -> Vec<crate::bopds::ds::IntersectionCurve> {
 use rcad_kernel::geom::Curve2dEval;
 use std::f64::consts::{TAU, PI};

 // OCCT L906-920: create analytic 3D curve from GLine.

 // OCCT L922-950: TreatCircle — split intervals with 0-crossing handling.
 // Sorts vertices on the GLine, creates intervals, tests midpoints.
 let parts = self.treat_circle_parts(curve, orig_t_range, typl, vertices, f1, f2);

 // OCCT L950-1095: aNbParts = seqp.Length() / 2.
 //   If aNbParts == 0 → the for loop does not execute → no output curves.
 if parts.is_empty() {
 return Vec::new();
 }

 let aPeriod = TAU;
 let aRealEpsilon = f64::EPSILON;
 let aNbParts = parts.len();
 let mut result = Vec::with_capacity(parts.len());

 for &[fprm, lprm] in &parts {
 // OCCT L953-956: if (|fprm|>eps || |lprm-2π|>eps) → not full-period
 let is_full_period = fprm.abs() <= aRealEpsilon && (lprm - aPeriod).abs() <= aRealEpsilon;

 if !is_full_period && (lprm > fprm + 1e-12) {
   // ── Not full-period: trimmed curve + pcurves + vertices ──
   // OCCT L960-990: Geom_TrimmedCurve(newc, fprm, lprm) + BuildPCurves + append.

   // OCCT L968-990: BuildPCurves with surface UV bounds for Circle/Ellipse.
   let pca = self.compute_pcurve_on_surface(curve, f1);
   let pcb = self.compute_pcurve_on_surface(curve, f2);
   let trimmed_pca = pca.as_ref().map(|pc| pc.clone());
   let trimmed_pcb = pcb.as_ref().map(|pc| pc.clone());

   // Create endpoint vertices.
   let (sv, ev) = {
   let p_start = curve.point_at(fprm);
   let p_end = curve.point_at(lprm);
   if p_start.is_finite() && p_end.is_finite() {
     let sv = self.ds.vertices.len();
     self.ds.vertices.push(DSVertex { point: p_start, geom_tol: geom_tol.max(TOLERANCE_ABS), origin: None, is_internal: false, location: 0 });
     let ev = self.ds.vertices.len();
     self.ds.vertices.push(DSVertex { point: p_end, geom_tol: geom_tol.max(TOLERANCE_ABS), origin: None, is_internal: false, location: 0 });
     (sv, ev)
   } else { (usize::MAX, usize::MAX) }
   };

   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
   curve: curve.clone(),
   polyline: Vec::new(),
   start_vertex: sv,
   end_vertex: ev,
   t_range: [fprm, lprm],
   pcurve_on_a: trimmed_pca,
   pcurve_on_b: trimmed_pcb,
   geom_tol: geom_tol.max(TOLERANCE_ABS),
   pave_blocks: Vec::new(),
   curve_extra,
   });

 } else if is_full_period && aNbParts == 1 {
   // ── Full-period, single part → accept full circle ──
   // OCCT L996-1042: trimmed full circle + BuildPCurves + append + break.
   let pca = self.compute_pcurve_on_surface(curve, f1);
   let pcb = self.compute_pcurve_on_surface(curve, f2);
   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
   curve: curve.clone(),
   polyline: Vec::new(),
   start_vertex: usize::MAX,
   end_vertex: usize::MAX,
   t_range: orig_t_range,
   pcurve_on_a: pca,
   pcurve_on_b: pcb,
   geom_tol: geom_tol.max(TOLERANCE_ABS),
   pave_blocks: Vec::new(),
   curve_extra,
   });
   break;

 } else if is_full_period && aNbParts > 1 {
   // ── Full-period, multiple parts: test 18 points ──
   // OCCT L1045-1095: on regarde si on garde.
   let aTwoPIdiv17 = aPeriod / 17.0;
   for j in 0..=17 {
   let t = j as f64 * aTwoPIdiv17;
   let p3d = curve.point_at(t);
   if !p3d.is_finite() { continue; }
   let uv1 = self.context.proj_ps(self.ds, f1, p3d).map(|(uv, _, _)| uv);
   let uv2 = self.context.proj_ps(self.ds, f2, p3d).map(|(uv, _, _)| uv);
   let in1 = uv1.map_or(false, |uv| { self.context.is_point_in_on_face(self.ds, f1, uv) });
   let in2 = uv2.map_or(false, |uv| { self.context.is_point_in_on_face(self.ds, f2, uv) });
   if in1 && in2 {
     let pca = self.compute_pcurve_on_surface(curve, f1);
     let pcb = self.compute_pcurve_on_surface(curve, f2);
     let mut curve_extra = crate::bopds::ds::CurveExtra::default();
     curve_extra.tangential_tol = tang_tolerance;
     result.push(crate::bopds::ds::IntersectionCurve {
     curve: curve.clone(),
     polyline: Vec::new(),
     start_vertex: usize::MAX,
     end_vertex: usize::MAX,
     t_range: orig_t_range,
     pcurve_on_a: pca,
     pcurve_on_b: pcb,
     geom_tol: geom_tol.max(TOLERANCE_ABS),
     pave_blocks: Vec::new(),
     curve_extra,
     });
     break;
   }
   }
 }
 }
 result
}

/// OCCT-aligned: LineConstructor (GeomInt_LineConstructor::Perform for GLine).
/// Iterates over IntPatch_Point vertices on the line, tests the midpoint of
/// each adjacent-vertex interval on both face domains, keeps valid intervals.
///
/// With 0 vertices (nbvtx=0): OCCT's intrvtested flag stays false, and the
/// full parameter range [FirstParameter, LastParameter] is kept as one part.
/// The caller's test-point logic decides whether to keep or reject it.
fn line_constructor_parts(
  &self, _curve: &Curve3, orig_t_range: [f64; 2], _typl: IntPatchIType,
) -> Vec<[f64; 2]> {
 // OCCT: nbvtx = GeomInt_LineTool::NbVertex(L).
 //       intrvtested = false; for (i=1; i<nbvtx; ++i) {...}
 //       if (!intrvtested) { seqp.Append(FirstParameter(L), LastParameter(L)); }
 //
 // For both finite and infinite ranges: OCCT returns the full range.
 // The caller decides based on the test-point approach.
 vec![orig_t_range]
}

/// OCCT-aligned: TreatCircle (GeomInt_LineConstructor.cxx L481-560).
/// For Circle/Ellipse with vertices: sorts vertices by parameter in [0, 2π),
/// creates intervals between sorted vertices, tests midpoints on both face
/// domains.  Handles 0-crossing via PeriodicReparam + SeqFprm/SeqLprm.
///
/// Without vertices (nbvtx=0): OCCT creates a zero-initialized array of size 1,
/// the sort and interval-building steps produce no valid intervals, and seqp
/// remains empty.  The caller sees aNbParts=0 and creates no output curves.
/// This function matches that behavior — returns empty.
/// OCCT-aligned: PutPointsOnLine (IntPatch_Intersection.cxx L268-312).
/// For each analytic line, computes boundary-crossing points where the
/// curve projects to UV outside one of the face domains.  These points
/// become vertices on the line, used by TreatCircle to split the curve.
///
/// rcad: pcurve evaluation + face UV boundary detection.  OCCT uses
/// IntPatch_Point vertices from the VV/VE/VF intersection phases,
/// projected onto each GLine via PutPointsOnLine.
fn put_points_on_line(
  &mut self, f1: usize, f2: usize,
  line: &mut crate::inttools::int_patch_line::IntPatchLine,
) {
 use rcad_kernel::geom::Curve2dEval;
 if line.is_wline() { return; } // Walking lines handled separately.
 let typl = line.line_type;
 if typl == IntPatchIType::Restriction { return; }

 let curve = &line.curve;
 let t_range = line.t_range;
 let pca = self.compute_pcurve_on_surface(curve, f1);
 let pcb = self.compute_pcurve_on_surface(curve, f2);

 // Determine a practical sampling range.
 let (sweep_start, sweep_end) = match typl {
 IntPatchIType::Circle | IntPatchIType::Ellipse => (0.0, std::f64::consts::TAU),
 _ => {
   if t_range[0].is_finite() { (t_range[0], t_range[0] + 1e5) }
   else if t_range[1].is_finite() { (t_range[1] - 1e5, t_range[1]) }
   else { (-1e5, 1e5) }
 }
 };
 if !sweep_start.is_finite() || !sweep_end.is_finite() || sweep_end <= sweep_start { return; }

 // OCCT L272-282: for each line, for each point: project and classify.
 // rcad: sample the curve, evaluate pcurve, classify UV on both face domains.
 // Find transitions from valid→invalid or invalid→valid.
 const N_SAMPLES: usize = 501;
 let mut uv_in = vec![false; N_SAMPLES];

 for i in 0..N_SAMPLES {
 let t = sweep_start + (sweep_end - sweep_start) * i as f64 / (N_SAMPLES - 1) as f64;
 let p3d = curve.point_at(t);
 if !p3d.is_finite() { continue; }

 // Prefer pcurve for UV, fall back to 3D projection.
 let uv1 = pca.as_ref().and_then(|pc| { let uv = pc.point_at(t); if uv.is_finite() { Some(uv) } else { None } })
   .or_else(|| self.context.proj_ps(self.ds, f1, p3d).map(|(uv, _, _)| uv));
 let uv2 = pcb.as_ref().and_then(|pc| { let uv = pc.point_at(t); if uv.is_finite() { Some(uv) } else { None } })
   .or_else(|| self.context.proj_ps(self.ds, f2, p3d).map(|(uv, _, _)| uv));

 let in1 = uv1.map_or(false, |uv| self.context.is_point_in_on_face(self.ds, f1, uv));
 let in2 = uv2.map_or(false, |uv| self.context.is_point_in_on_face(self.ds, f2, uv));
 uv_in[i] = in1 && in2;
 }

 // OCCT L290-310: find transitions and add vertices.
 // Add a vertex at each transition point (middle of the transition interval).
 let mut new_vertices: Vec<crate::inttools::int_patch_line::IntPatchVertex> = Vec::new();
 for i in 1..N_SAMPLES {
 if uv_in[i] != uv_in[i - 1] {
   // Transition at midpoint between sample i-1 and i.
   let t_mid = sweep_start + (sweep_end - sweep_start) * (i as f64 - 0.5) / (N_SAMPLES - 1) as f64;
   let p3d = curve.point_at(t_mid);
   if p3d.is_finite() {
   new_vertices.push(crate::inttools::int_patch_line::IntPatchVertex {
     param_on_line: t_mid,
     p3d,
     u1: pca.as_ref().map_or(0.0, |pc| pc.point_at(t_mid).x),
     v1: pca.as_ref().map_or(0.0, |pc| pc.point_at(t_mid).y),
     u2: pcb.as_ref().map_or(0.0, |pc| pc.point_at(t_mid).x),
     v2: pcb.as_ref().map_or(0.0, |pc| pc.point_at(t_mid).y),
   });
   }
 }
 }

 // Deduplicate vertices: OCCT L312-320 rejects vertices within TolPC.
 const VERTEX_TOL: f64 = 1e-10;
 for v in new_vertices {
 let is_dup = line.vertices.iter().any(|ev| (ev.param_on_line - v.param_on_line).abs() < VERTEX_TOL);
 if !is_dup {
   line.vertices.push(v);
 }
 }

 // Sort vertices by parameter (OCCT L322: ComputeVertexParameters).
 line.vertices.sort_by(|a, b| a.param_on_line.partial_cmp(&b.param_on_line).unwrap_or(std::cmp::Ordering::Equal));
}

fn treat_circle_parts(
  &mut self, curve: &Curve3, orig_t_range: [f64; 2], typl: IntPatchIType,
  vertices: &[crate::inttools::int_patch_line::IntPatchVertex],
  f1: usize, f2: usize,
) -> Vec<[f64; 2]> {
 use crate::inttools::int_patch_line::IntPatchVertex;
 use std::f64::consts::TAU;

 // OCCT GeomInt_LineConstructor::TreatCircle (L481-560):
 //   aVtxArr = sorted vertices with params projected to [0, 2π)
 //   Last vertex = first.param + 2π  (wraps around)
 //   Remove duplicates within aTolPC
 //   Sort again
 //   For each adjacent pair (i, i+1): test midpoint on both face domains

 let aNbVtx = vertices.len();
 if aNbVtx == 0 {
  // OCCT: with 0 vertices, aVtxArr has size 1 (default-constructed),
  // sort does nothing, no intervals produced → seqp empty → caller
  // sees aNbParts=0 → no output curves.
  return Vec::new();
 }

 // Build vertex array with parameters projected to [0, 2π) (OCCT L492-495).
 let mut sorted: Vec<(f64, &IntPatchVertex)> = vertices.iter()
   .map(|v| {
     let par = if orig_t_range[1] - orig_t_range[0] >= TAU - 1e-12 {
       let p = v.param_on_line;
       if p < 0.0 { p + TAU * (( -p / TAU).ceil()) }
       else if p >= TAU { p - TAU * ((p / TAU).floor()) }
       else { p }
     } else { v.param_on_line };
     (par, v)
   })
   .collect();

 // OCCT L500-502: sort by parameter.
 sorted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

 // OCCT L504: create last vertex at first.param + 2π.
 let first_param = sorted[0].0;
 let last_param = first_param + TAU;

 // OCCT L506-515: reject duplicates (within aTolPC = 1000*Precision::PConfusion).
 const aTolPC: f64 = 1000.0 * 1e-12;
 let mut deduped: Vec<f64> = Vec::with_capacity(sorted.len() + 1);
 for &(par, _) in &sorted {
  if deduped.is_empty() || (par - *deduped.last().unwrap()).abs() > aTolPC {
   deduped.push(par);
  }
 }
 // Add the last vertex (first.param + 2π), skip if duplicate of last deduped.
 let last_2pi = first_param + TAU;
 if deduped.is_empty() || (last_2pi - *deduped.last().unwrap()).abs() > aTolPC {
  deduped.push(last_2pi);
 }

 // OCCT L525-560: for each adjacent pair, test midpoint on both face domains.
 let mut result = Vec::new();
 for i in 0..(deduped.len().saturating_sub(1)) {
  let t1 = deduped[i];
  let t2 = deduped[i + 1];
  if (t2 - t1).abs() <= 1e-12 { continue; }
  // Test midpoint (OCCT L539-548).
  let t_mid = (t1 + t2) * 0.5;
  let p3d = curve.point_at(t_mid);
  if !p3d.is_finite() { continue; }
  // OCCT: Parameters(myHS1, myHS2, Pmid, u1, v1, u2, v2)
  // AdjustPeriodic for periodic surfaces, then Classify on both domains.
  // Project midpoint to UV on both faces and classify on their domains.
  let uv1 = self.context.proj_ps(self.ds, f1, p3d).map(|(uv, _, _)| uv);
  let uv2 = self.context.proj_ps(self.ds, f2, p3d).map(|(uv, _, _)| uv);
  let in1 = uv1.map_or(false, |uv| self.context.is_point_in_on_face(self.ds, f1, uv));
  let in2 = uv2.map_or(false, |uv| self.context.is_point_in_on_face(self.ds, f2, uv));
  if !(in1 && in2) { continue; }
  result.push([t1, t2]);
 }

 result
}


/// OCCT-aligned: BuildPCurves for all curve-surface type combinations.
/// Matches OCCT GeomInt_IntSS::BuildPCurves (L822-846). Uses exact
/// analytic pcurves when available; falls back to sampling + projection.
fn compute_pcurve_on_surface(
  &self, curve: &rcad_kernel::geom::Curve3, fi: usize,
) -> Option<rcad_kernel::geom::Curve2d> {
 if fi >= self.ds.faces.len() { return None; }
 let surf = &self.ds.faces[fi].surface;
 let pc = match (curve, surf) {
 (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Plane(p)) =>
   crate::inttools::pcurve_derive::line_pcurve_on_plane(l, p),
 (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Plane(p)) =>
   crate::inttools::pcurve_derive::circle_pcurve_on_plane(c, p),
 (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Plane(p)) =>
   crate::inttools::pcurve_derive::ellipse_pcurve_on_plane(e, p),
 (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Sphere(s)) =>
   crate::inttools::pcurve_derive::line_pcurve_on_sphere(l, s),
 (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Sphere(s)) =>
   crate::inttools::pcurve_derive::circle_pcurve_on_sphere(c, s),
 (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Sphere(s)) =>
   crate::inttools::pcurve_derive::ellipse_pcurve_on_sphere(e, s),
 (rcad_kernel::geom::Curve3::Parabola(p), rcad_kernel::geom::Surface3::Sphere(s)) =>
   crate::inttools::pcurve_derive::parabola_pcurve_on_sphere(p, s),
 (rcad_kernel::geom::Curve3::Hyperbola(h), rcad_kernel::geom::Surface3::Sphere(s)) =>
   crate::inttools::pcurve_derive::hyperbola_pcurve_on_sphere(h, s),
 (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Cylinder(c)) =>
   crate::inttools::pcurve_derive::line_pcurve_on_cylinder(l, c),
 (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Cylinder(cyl)) =>
   crate::inttools::pcurve_derive::circle_pcurve_on_cylinder(c, cyl),
 (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Cylinder(cyl)) =>
   crate::inttools::pcurve_derive::ellipse_pcurve_on_cylinder(e, cyl),
 (rcad_kernel::geom::Curve3::Line(l), rcad_kernel::geom::Surface3::Cone(c)) =>
   crate::inttools::pcurve_derive::line_pcurve_on_cone(l, c),
 (rcad_kernel::geom::Curve3::Circle(c), rcad_kernel::geom::Surface3::Cone(co)) =>
   crate::inttools::pcurve_derive::circle_pcurve_on_cone(c, co),
 (rcad_kernel::geom::Curve3::Ellipse(e), rcad_kernel::geom::Surface3::Cone(co)) =>
   crate::inttools::pcurve_derive::ellipse_pcurve_on_cone(e, co),
 _ => {
   let tr = match curve {
   rcad_kernel::geom::Curve3::Line(_) => [-1e3, 1e3],
   rcad_kernel::geom::Curve3::Circle(_) | rcad_kernel::geom::Curve3::Ellipse(_) => [0.0, std::f64::consts::TAU],
   _ => [0.0, 1.0],
   };
   crate::inttools::pcurve_derive::fallback_pcurve_by_projection(curve, &tr, surf)
 }
 };
 Some(pc)
}
}
