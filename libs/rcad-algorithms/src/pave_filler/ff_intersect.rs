use super::*;

impl<'a> super::PaveFiller<'a> {
 pub(crate) fn perform_ff(&mut self) {
 // OCCT PaveFiller_6.cxx L288-314: UpdateFaceInfoOn/In for all FF participant faces
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

 if let (Some(bvh_a), Some(bvh_b)) = (self.bvh_a, self.bvh_b) {
 // Build reverse maps: BRep face index  ?position in a_faces/b_faces
 let a_max_idx = a_faces.iter().map(|&dsi| self.ds.faces[dsi].source_face_idx).max().unwrap_or(0);
 let b_max_idx = b_faces.iter().map(|&dsi| self.ds.faces[dsi].source_face_idx).max().unwrap_or(0);
 let mut a_rev = vec![usize::MAX; a_max_idx + 1];
 for (pos, &dsi) in a_faces.iter().enumerate() {
 a_rev[self.ds.faces[dsi].source_face_idx] = pos;
 }
 let mut b_rev = vec![usize::MAX; b_max_idx + 1];
 for (pos, &dsi) in b_faces.iter().enumerate() {
 b_rev[self.ds.faces[dsi].source_face_idx] = pos;
 }

 let candidates = Bvh::candidate_pairs(bvh_a, bvh_b);
 let mut processed_pairs = std::collections::HashSet::new();
 for (fa_brep, fb_brep) in candidates {
 if let (Some(&ai), Some(&bi)) = (a_rev.get(fa_brep), b_rev.get(fb_brep))
 && ai != usize::MAX && bi != usize::MAX {
 //  ?OCCT-aligned: BVH may produce duplicate candidate pairs when a face appears
 // in multiple intersecting leaf nodes, causing duplicate intersection curves.
 // OCCT PaveFiller processes each face pair once (FF matrix uses BOPDS_IndexRange
 // to mark pairs as already processed).
 if !processed_pairs.insert((ai, bi)) { continue; }
 let af = a_faces[ai];
 let bf = b_faces[bi];
 if self.should_skip_glued_face_pair(af, bf) {
 continue;
 }
 self.intersect_face_face(af, bf);
 }
 }
 } else {
 //  ?OCCT-aligned: BOPDS_Iterator cross-group face pair iteration.
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
 // OCCT BOPAlgo_PaveFiller_6.cxx L106-134 IsClosedFF:
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
 /// OCCT: dispatch FF intersection by surface type
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

  // OCCT-aligned: IntPatch_Intersection: generic surface-surface intersection.
  let mut int_patch = crate::inttools::int_patch_intersection::IntPatchIntersection::new();
  int_patch.perform(&s1, &s2, self.fuzzy_tolerance, self.fuzzy_tolerance);
  // OCCT L536-540: tangent faces -> early return, no curves.
  if int_patch.tangent_faces() {
  self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
  f1, f2, curves: Vec::new(), points: Vec::new(), tangent_faces: true,
  });
  if let Some(ref info) = shift_info { self.reverse_seam_edge_shift(f1, f2, info); }
  self.seam_shift_tol = old_shift_tol;
  return;
  }
  // OCCT-aligned: MakeCurve (IntTools_FaceFace.cxx L695-1846) for each IntPatch line.
  let mut ff_curve_indices: Vec<usize> = Vec::new();
  for i in 0..int_patch.nb_lines() {
  let ci = self.ds.intersection_curves.len();
  let ic = self.make_intersection_curve(f1, f2, int_patch.line(i));
  self.ds.intersection_curves.push(ic);
  ff_curve_indices.push(ci);
  }
  // OCCT L540-550: one InterferenceFF per face pair, all curves in one entry.
  // OCCT L565-607: also process intersection points (PntOn2Faces).
  let mut ff_point_indices: Vec<usize> = Vec::new();
  for pi in 0..int_patch.nb_points() {
  let pt = int_patch.point(pi);
  // OCCT L578-586: validate point is within both face domains
  if !self.context.is_point_in_on_face(self.ds, f1, glam::DVec2::new(pt.u1, pt.v1)) { continue; }
  if !self.context.is_point_in_on_face(self.ds, f2, glam::DVec2::new(pt.u2, pt.v2)) { continue; }
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
 // OCCT L629-638: get curve, pcurves, starting tolerance, parameter range
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
 // OCCT L679-689: SetTolerance + SetTangentialTolerance
 let ic = &mut self.ds.intersection_curves[ci];
 ic.geom_tol = ic.geom_tol.max(new_tol);
 ic.curve_extra.tangential_tol = ic.curve_extra.tangential_tol.max(tang_tol);
 }
 // PrepareLines3D  ?split closed curves
 inttools::pcurve_derive::prepare_lines_3d(&mut self.ds.intersection_curves);
 //  ?OCCT-aligned: After PrepareLines3D splits closed curves, new curve endpoints
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
 // Push new vertices directly  ?do NOT use add_vertex which merges
 // near-coincident points.  Multiple sphere-plane intersection circles
 // pass through the same spatial positions (e.g. (1,0,0), (0,1,0),
 // (0,0, 1)), and each half-circle needs its OWN start/end vertices
 // to form distinct arcs.  Merging them collapses different
 // face-boundary arcs onto the same vertex pair, creating the
 // same section edges for different curves and breaking the loop.
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
 // Only add valid vertex indices to vertices_in.  Unclipped
 // analytical curves may still have usize::MAX endpoints.
 if ic.start_vertex < self.ds.vertices.len() {
 self.ds.faces[f1].face_info.vertices_in.insert(ic.start_vertex);
 self.ds.faces[f2].face_info.vertices_in.insert(ic.start_vertex);
 }
 if ic.end_vertex < self.ds.vertices.len() {
 self.ds.faces[f1].face_info.vertices_in.insert(ic.end_vertex);
 self.ds.faces[f2].face_info.vertices_in.insert(ic.end_vertex);
 }
 // Register curve endpoints as vertices_on if they match face boundary vertices
 let bv1 = self.ds.faces[f1].boundary_verts.clone();
 let bv2 = self.ds.faces[f2].boundary_verts.clone();
 for &bvi in &bv1 {
 if bvi == ic.start_vertex || bvi == ic.end_vertex {
 self.ds.faces[f1].face_info.vertices_on.insert(bvi);
 }
 }
 for &bvi in &bv2 {
 if bvi == ic.start_vertex || bvi == ic.end_vertex {
 self.ds.faces[f2].face_info.vertices_on.insert(bvi);
 }
 }
 }
 }
 //  ?OCCT-aligned: InitPaveBlock1 for all curves (PaveFiller_6.cxx L800).
 // Creates an initial PaveBlock on each curve for ext_pave tracking.
 for ci in 0..self.ds.intersection_curves.len() {
 
 }
 }
 }

 /// OCCT: plane-plane intersection
 /// OCCT: plane-plane intersection

 ///  ?OCCT-aligned: CheckSelfInterference (BOPAlgo_PaveFiller_11.cxx L28-221).
 /// Builds vertex aces and edge aces connection maps per source range,
 /// detects acquired self-intersections (same vertex/edge used by >1 face from same operand).
 pub(crate) fn check_self_interference(&self) -> Vec<String> {
 // OCCT L30-34: single argument == self-interference mode, skip.
 if self.ds.a_vertex_count == 0 {
 return Vec::new();
 }

 // OCCT L38-41: iterate ranges (A and B operands).
 let a_end = self.ds.a_vertex_count;
 let ranges: [(usize, usize, &str); 2] = [
 (0, a_end, "A"),
 (a_end, self.ds.vertices.len(), "B"),
 ];

 let mut warnings: Vec<String> = Vec::new();

 for &(range_start, range_end, _name) in &ranges {
 // OCCT L43-48: aMCSI  ?map of connections: vertex/edge  ?list of faces.
 let mut v_to_faces: std::collections::HashMap<usize, Vec<usize>> =
 std::collections::HashMap::new();
 let mut e_to_faces: std::collections::HashMap<usize, Vec<usize>> =
 std::collections::HashMap::new();
 // OCCT L48: aMCBFence  ?skip already-processed CommonBlocks.
 let mut cb_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();

 // OCCT L51-197: iterate shapes in this range.
 // rcad: process faces whose source solids are in this range.
 let origin = if range_start == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
 for fi in 0..self.ds.faces.len() {
 let face = &self.ds.faces[fi];
 if face.origin != origin { continue; }

 // OCCT L151-173: FACE  ?analyze IN and SC vertices.
 for &vi in &face.face_info.vertices_in {
 v_to_faces.entry(vi).or_default().push(fi);
 }
 // OCCT L156-172: VerticesSc
 for &vi in &face.face_info.vertices_in {
 if !face.face_info.vertices_sc.is_empty() {
 // vertices_sc is the set of vertices of section curves
 }
 }

 // OCCT L175-195: PBsIn / PBsSc  ?edge ace connections.
 for pb_idx in face.face_info.pave_blocks_in.iter()
 .chain(face.face_info.pave_blocks_sc.iter())
 {
 if *pb_idx >= self.ds.pave_blocks.len() { continue; }
 let pb = &self.ds.pave_blocks[*pb_idx];
  let ei = pb.0.read().unwrap().original_edge;
  if ei >= self.ds.edges.len() { continue; }
  e_to_faces.entry(ei).or_default().push(fi);

  // OCCT L112-148: CommonBlock analysis  ?check if same CB
  // contains edges from the same argument.
  if let Some(cb_idx) = pb.0.read().unwrap().common_block_idx {
 if cb_fence.insert(cb_idx) {
 if cb_idx < self.ds.common_blocks.len() {
 let cb = &self.ds.common_blocks[cb_idx];
 let same_arg_edges: Vec<usize> = cb.pave_blocks().iter()
 .filter_map(|&(pbi, _)| {
 let pb2 = &self.ds.pave_blocks[pbi];
 let e = pb2.0.read().unwrap().original_edge;
 if e < self.ds.edges.len()
 && self.ds.edges[e].origin == face.origin
 {
 Some(e)
 } else { None }
 })
 .collect();
 if same_arg_edges.len() > 1 {
 warnings.push(format!(
 "Acquired self-intersection: CommonBlock {:?} contains {} edges from same argument",
 cb_idx, same_arg_edges.len()
 ));
 }
 }
 }
 }
 }
 }

 // OCCT L198-219: Analyze connections  ?if any vertex/edge connects
 // >1 face from the same argument  ?self-interference.
 for (_vi, faces) in &v_to_faces {
 if faces.len() > 1 {
 warnings.push(format!(
 "Self-interference: vertex {:?} belongs to {} faces from same argument",
 _vi, faces.len()
 ));
 }
 }
 for (_ei, faces) in &e_to_faces {
 if faces.len() > 1 {
 warnings.push(format!(
 "Self-interference: edge {:?} belongs to {} faces from same argument",
 _ei, faces.len()
 ));
 }
 }
 }

 warnings
 }

 pub(crate) fn make_split_edges(&mut self) {
 // OCCT L392: UpdateCommonBlocksWithSDVertices  ?before creating split edges,
 // ensure CommonBlocks reference correct (SD-deduplicated) vertex indices.
 self.ds.update_common_blocks_with_sd_vertices();

 // Phase 1: collect PaveBlock data without creating new edges (avoids
 // mutable borrow conflict with self.ds.edges iteration).
 struct BlockData {
 ei: usize,
 sv: usize, ev: usize,
 t_start: f64, t_end: f64,
 curve: Curve3,
 origin: ShapeOrigin,
 geom_tol: f64,
 face_reps: Vec<DSCurveRepOnFace>,
 }
 let mut all_blocks: Vec<BlockData> = Vec::new();
 let n_orig_edges = self.ds.edges.len();

 //  ?OCCT-aligned: MakeSplitEdges (PaveFiller_7.cxx) only creates split
 // edges and sets PaveBlock->Edge() (pb.0.read().unwrap().new_edge).  rcad also initializes
 // pave_blocks on source edges here so downstream FillImagesEdges can
 // read pb.0.read().unwrap().new_edge.  my_images / my_origins are NOT populated here  ?
 // that is FillImagesEdges' responsibility (build_edge_images in ds.rs).

 for ei in 0..n_orig_edges {
 let edge = &self.ds.edges[ei];
 if edge.paves.is_empty() {
 continue;
 }

 // OCCT L408-414: skip degenerated edges (HasFlag).
 if self.ds.is_edge_degenerated(ei) {
 continue;
 }

 let mut all_paves = vec![
 Pave { vertex_idx: edge.start_vertex, param: edge.t_range[0] },
 Pave { vertex_idx: edge.end_vertex, param: edge.t_range[1] },
 ];
 all_paves.extend_from_slice(&edge.paves);
 all_paves.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));
 all_paves.dedup_by(|a, b| params_equal(a.param, b.param));

 // OCCT L406-421: CommonBlock dedup  ?skip PBs whose CommonBlock was already processed.
 let mut processed_common_blocks: std::collections::HashSet<usize> = std::collections::HashSet::new();

 for w in all_paves.windows(2) {
  let pb = PaveBlock::new(ei, w[0], w[1]);

  // OCCT L416-421: skip PaveBlock whose CommonBlock already processed.
  if let Some(cb_idx) = pb.common_block_idx {
  if !processed_common_blocks.insert(cb_idx) {
  continue;
  }
  }

  // OCCT L425-430: skip if no new vertices (both vertices are from source shapes).
  if !self.ds.is_new_vertex(pb.pave1.vertex_idx)
  && !self.ds.is_new_vertex(pb.pave2.vertex_idx)
  {
  continue;
  }

  let t1 = pb.pave1.param;
  let t2 = pb.pave2.param;
  let (t_start, t_end) = if t1 < t2 { (t1, t2) } else { (t2, t1) };
  let split_curve = pb.curve.clone().unwrap_or_else(|| edge.curve.clone());
  all_blocks.push(BlockData {
  ei,
  sv: pb.pave1.vertex_idx,
  ev: pb.pave2.vertex_idx,
 t_start, t_end,
 curve: split_curve,
 origin: edge.origin,
 geom_tol: edge.geom_tol,
 face_reps: edge.face_reps.clone(),
 });
 }
 }

 // Phase 2: create new DSEdges for each collected block + set pave_blocks
 // on source edges (MakeSplitEdges).  my_images / my_origins are NOT
 // populated here  ?that is FillImagesEdges' job (build_edge_images in ds.rs).
 let mut edge_pbs: std::collections::HashMap<usize, Vec<(usize, usize, f64, f64, usize)>> =
 std::collections::HashMap::new();

 for data in &all_blocks {
 let new_ei = self.ds.edges.len();
 self.ds.edges.push(DSEdge {
 start_vertex: data.sv,
 end_vertex: data.ev,
 curve: data.curve.clone(),
 t_range: [data.t_start, data.t_end],
 origin: data.origin,
 geom_tol: data.geom_tol,
 paves: vec![],
 pave_blocks: vec![],
 face_reps: data.face_reps.clone(),
 is_internal: false,
 vertex_params: {
 let mut vp = std::collections::HashMap::new();
 vp.insert(data.sv, data.t_start);
 vp.insert(data.ev, data.t_end);
 vp
 },
 face_tolerances: Vec::new(),
  is_geometric: true,
  location: 0,
  });

  // Track for pave_blocks assignment on source edge
 edge_pbs.entry(data.ei).or_default().push((
 data.sv, data.ev, data.t_start, data.t_end, new_ei,
 ));
 }

 //  ?OCCT-aligned: Set pave_blocks on source edges that were split,
 // so Builder::fill_images_edges can read pb.0.read().unwrap().new_edge.
 for (ei, blocks) in &edge_pbs {
 let pbs: Vec<PaveBlock> = blocks.iter().map(|&(sv, ev, t_start, t_end, new_ei)| {
 let mut pb = PaveBlock::new(*ei,
 Pave { vertex_idx: sv, param: t_start },
 Pave { vertex_idx: ev, param: t_end },
 );
  pb.new_edge = Some(new_ei);
 pb
 }).collect();
 self.ds.edges[*ei].pave_blocks = pbs.into_iter().map(|pb| crate::bopds::pave::SharedPB::new(pb)).collect();
 }
 }

 // = = =  Helpers = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

 pub(crate) fn verts_of(&self, origin: ShapeOrigin) -> Vec<usize> {
 self.ds
 .vertices
 .iter()
 .enumerate()
 .filter(|(_, v)| v.origin == Some(origin))
 .map(|(i, _)| i)
 .collect()
 }

 pub(crate) fn edges_of(&self, origin: ShapeOrigin) -> Vec<usize> {
 self.ds
 .edges
 .iter()
 .enumerate()
 .filter(|(_, e)| e.origin == origin)
 .map(|(i, _)| i)
 .collect()
 }

 pub(crate) fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
 self.ds
 .faces
 .iter()
 .enumerate()
 .filter(|(_, f)| f.origin == origin)
 .map(|(i, _)| i)
 .collect()
 }

 /// OCCT-aligned: MakeCurve for IntPatch analytic lines
/// (IntTools_FaceFace.cxx L695-1846). Builds an IntersectionCurve with
/// proper t_range, pcurves, and endpoint vertices from an IntPatchLine.
/// For infinite-range Line3 curves: clips to face UV bounds using the
/// boundary-vertex bounding box.
pub(crate) fn make_intersection_curve(
  &mut self, f1: usize, f2: usize,
  line: &crate::inttools::int_patch_line::IntPatchLine,
) -> crate::bopds::ds::IntersectionCurve {
 let curve = line.curve.clone();
 let mut t_range = line.t_range;
 let mut pcurve_on_a: Option<rcad_kernel::geom::Curve2d> = line.pcurve1.clone();
 let mut pcurve_on_b: Option<rcad_kernel::geom::Curve2d> = line.pcurve2.clone();
 let mut start_vertex = usize::MAX;
 let mut end_vertex = usize::MAX;
 let mut geom_tol = line.tolerance;
 let mut polyline: Vec<glam::DVec3> = Vec::new();

 // Step 0: For walking lines (IntPatch_Walking), create BSpline from points.
 // OCCT MakeCurve L724-744: handles IntPatch_Walking by approximating a
 // BSpline from the walking points.  The pcurves are already known from
 // the marching process (wline_pnts stores UV on both surfaces).
 if line.is_wline() {
 let n = line.nb_points();
 if n >= 2 {
 let p3d_pts: Vec<glam::DVec3> = (0..n).map(|i| line.point(i).p3d).collect();
 polyline = p3d_pts.clone();
 if let Some(bs_curve3) = crate::inttools::intss::polyline_to_bspline(&p3d_pts, 1e-4) {
 let t_range_bs = bs_curve3.default_domain();
 let bs = match &bs_curve3 {
 rcad_kernel::geom::Curve3::BSpline(b) => b.clone(),
 _ => { return crate::bopds::ds::IntersectionCurve {
   curve: bs_curve3.clone(), polyline, start_vertex: usize::MAX, end_vertex: usize::MAX,
   t_range: t_range_bs, pcurve_on_a, pcurve_on_b, geom_tol: line.tolerance.max(1e-7),
   pave_blocks: Vec::new(), curve_extra: crate::bopds::ds::CurveExtra { tangential_tol: line.tang_tolerance, ..Default::default() },
 }; }
 };
 // For walking lines with pcurve info, build pcurves from point UV data
 if n >= 2 && pcurve_on_a.is_none() && line.point(0).u1.is_finite() {
   let uv_samples: Vec<(f64, f64, f64)> = (0..n).map(|i| {
   let p = line.point(i);
   let t = t_range_bs[0] + (t_range_bs[1] - t_range_bs[0]) * i as f64 / (n - 1) as f64;
   (t, p.u1, p.v1)
   }).collect();
   let knots: Vec<f64> = uv_samples.iter().map(|(t, _, _)| *t).collect();
   let uv_pts: Vec<glam::DVec2> = uv_samples.iter().map(|(_, u, v)| glam::DVec2::new(*u, *v)).collect();
   if let Ok(bs2d) = rcad_kernel::fit::interpolate_points_2d(&uv_pts) {
   let mut bs2d_clone = bs2d.clone();
   for (k, &tk) in knots.iter().enumerate() {
   if k < bs2d_clone.knots.len() { bs2d_clone.knots[k] = tk; }
   }
   pcurve_on_a = Some(rcad_kernel::geom::Curve2d::BSpline(bs2d_clone));
   }
 }
 if n >= 2 && pcurve_on_b.is_none() && line.point(0).u2.is_finite() {
   let uv_samples: Vec<(f64, f64, f64)> = (0..n).map(|i| {
   let p = line.point(i);
   let t = t_range_bs[0] + (t_range_bs[1] - t_range_bs[0]) * i as f64 / (n - 1) as f64;
   (t, p.u2, p.v2)
   }).collect();
   let uv_pts: Vec<glam::DVec2> = uv_samples.iter().map(|(_, u, v)| glam::DVec2::new(*u, *v)).collect();
   if let Ok(bs2d) = rcad_kernel::fit::interpolate_points_2d(&uv_pts) {
   pcurve_on_b = Some(rcad_kernel::geom::Curve2d::BSpline(bs2d));
   }
 }
 return crate::bopds::ds::IntersectionCurve {
   curve: rcad_kernel::geom::Curve3::BSpline(bs),
   polyline,
   start_vertex: usize::MAX, end_vertex: usize::MAX,
   t_range: t_range_bs,
   pcurve_on_a, pcurve_on_b,
   geom_tol: line.tolerance.max(1e-7),
   pave_blocks: Vec::new(),
   curve_extra: crate::bopds::ds::CurveExtra { tangential_tol: line.tang_tolerance, ..Default::default() },
 };
 }
 }
 }

 // Step 1: For infinite-range Line3 curves, clip to face UV bounds.
 let needs_clipping = matches!(&curve, rcad_kernel::geom::Curve3::Line(_))
   && t_range[0] <= -1e9 && t_range[1] >= 1e9;
 if needs_clipping {
 let dir = match &curve {
 rcad_kernel::geom::Curve3::Line(l) => l.direction,
 _ => unreachable!(),
 };
 if !dir.is_finite() || dir.length_squared() < 1e-30 {
 return self.make_raw_intersection_curve(curve, t_range, pcurve_on_a, pcurve_on_b, geom_tol, line.tang_tolerance);
 }
 let line_origin = match &curve {
 rcad_kernel::geom::Curve3::Line(l) => l.origin,
 _ => unreachable!(),
 };
 let mut t_min = f64::NEG_INFINITY;
 let mut t_max = f64::INFINITY;
 for &fi in &[f1, f2] {
 if fi >= self.ds.faces.len() { continue; }
 let face = &self.ds.faces[fi];
 let (base, x_axis, y_axis) = match &face.surface {
 rcad_kernel::geom::Surface3::Plane(p) => {
 let abs = p.normal.abs();
 let candidate = if abs.x <= abs.y && abs.x <= abs.z { glam::DVec3::X }
 else if abs.y <= abs.z { glam::DVec3::Y }
 else { glam::DVec3::Z };
 let x = p.normal.cross(candidate).normalize();
 let y = p.normal.cross(x);
 (p.origin, x, y)
 }
 _ => continue,
 };
 let mut u_min = f64::MAX; let mut u_max = f64::NEG_INFINITY;
 let mut v_min = f64::MAX; let mut v_max = f64::NEG_INFINITY;
 let mut has_uv = false;
 for &vi in &face.boundary_verts {
 if vi < self.ds.vertices.len() {
 let pt = self.ds.vertices[vi].point;
 let u = (pt - base).dot(x_axis);
 let v = (pt - base).dot(y_axis);
 u_min = u_min.min(u); u_max = u_max.max(u);
 v_min = v_min.min(v); v_max = v_max.max(v);
 has_uv = true;
 }
 }
 if !has_uv || !(u_max > u_min + 1e-12) || !(v_max > v_min + 1e-12) { continue; }
 let d_base = line_origin - base;
 let u0 = d_base.dot(x_axis);
 let du = dir.dot(x_axis);
 let v0 = d_base.dot(y_axis);
 let dv = dir.dot(y_axis);
 if du.abs() > 1e-30 {
 let t_lo = (u_min - u0) / du;
 let t_hi = (u_max - u0) / du;
 t_min = t_min.max(t_lo.min(t_hi)); t_max = t_max.min(t_lo.max(t_hi));
 } else if u0 < u_min - 1e-12 || u0 > u_max + 1e-12 {
 return self.make_raw_intersection_curve(curve, t_range, pcurve_on_a, pcurve_on_b, geom_tol, line.tang_tolerance);
 }
 if dv.abs() > 1e-30 {
 let t_lo = (v_min - v0) / dv;
 let t_hi = (v_max - v0) / dv;
 t_min = t_min.max(t_lo.min(t_hi)); t_max = t_max.min(t_lo.max(t_hi));
 } else if v0 < v_min - 1e-12 || v0 > v_max + 1e-12 {
 return self.make_raw_intersection_curve(curve, t_range, pcurve_on_a, pcurve_on_b, geom_tol, line.tang_tolerance);
 }
 }
 if t_min.is_finite() && t_max.is_finite() && t_max > t_min + 1e-12 {
 let p_start = curve.point_at(t_min);
 let p_end = curve.point_at(t_max);
 if p_start.is_finite() && p_end.is_finite() {
 t_range = [t_min, t_max];
 start_vertex = self.ds.vertices.len();
 self.ds.vertices.push(crate::bopds::ds::DSVertex { point: p_start, geom_tol: crate::tolerance::TOLERANCE_ABS, origin: None, is_internal: false, location: 0 });
 end_vertex = self.ds.vertices.len();
 self.ds.vertices.push(crate::bopds::ds::DSVertex { point: p_end, geom_tol: crate::tolerance::TOLERANCE_ABS, origin: None, is_internal: false, location: 0 });
 }
 }
 }

 // Step 2: Compute pcurves for analytic curves on Plane surfaces.
 // OCCT BuildPCurves (L822-846): projects 3D curve to each surface UV.
 if pcurve_on_a.is_none() {
 pcurve_on_a = self.compute_pcurve_on_surface(&curve, f1);
 }
 if pcurve_on_b.is_none() {
 pcurve_on_b = self.compute_pcurve_on_surface(&curve, f2);
 }

 // Step 3: Create endpoint vertices at finite t_range boundaries.
 // Skip Circle/Ellipse -- prepare_lines_3d splits closed curves
 // and needs_fix creates the correct endpoint vertices.
 let is_closed_arc = matches!(&curve, rcad_kernel::geom::Curve3::Circle(_) | rcad_kernel::geom::Curve3::Ellipse(_));
 if !is_closed_arc && start_vertex == usize::MAX && end_vertex == usize::MAX
   && t_range[0].is_finite() && t_range[1].is_finite()
   && t_range[1] > t_range[0] + 1e-12
 {
 let p_start = curve.point_at(t_range[0]);
 let p_end = curve.point_at(t_range[1]);
 if p_start.is_finite() && p_end.is_finite() {
 start_vertex = self.ds.vertices.len();
 self.ds.vertices.push(crate::bopds::ds::DSVertex { point: p_start, geom_tol: crate::tolerance::TOLERANCE_ABS, origin: None, is_internal: false, location: 0 });
 end_vertex = self.ds.vertices.len();
 self.ds.vertices.push(crate::bopds::ds::DSVertex { point: p_end, geom_tol: crate::tolerance::TOLERANCE_ABS, origin: None, is_internal: false, location: 0 });
 }
 }

 let mut curve_extra = crate::bopds::ds::CurveExtra::default();
 curve_extra.tangential_tol = line.tang_tolerance;
 crate::bopds::ds::IntersectionCurve {
 curve,
 polyline,
 start_vertex,
 end_vertex,
 t_range,
 pcurve_on_a,
 pcurve_on_b,
 geom_tol,
 pave_blocks: Vec::new(),
 curve_extra,
 }
}

/// Fallback: create raw IntersectionCurve from the given parameters
/// without any clipping/pcurve computation (for unclippable curves).
fn make_raw_intersection_curve(
  &self,
  curve: rcad_kernel::geom::Curve3,
  t_range: [f64; 2],
  pcurve_on_a: Option<rcad_kernel::geom::Curve2d>,
  pcurve_on_b: Option<rcad_kernel::geom::Curve2d>,
  geom_tol: f64,
  tang_tolerance: f64,
) -> crate::bopds::ds::IntersectionCurve {
 let mut curve_extra = crate::bopds::ds::CurveExtra::default();
 curve_extra.tangential_tol = tang_tolerance;
 crate::bopds::ds::IntersectionCurve {
 curve, polyline: Vec::new(), start_vertex: usize::MAX, end_vertex: usize::MAX,
 t_range, pcurve_on_a, pcurve_on_b, geom_tol, pave_blocks: Vec::new(), curve_extra,
 }
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
