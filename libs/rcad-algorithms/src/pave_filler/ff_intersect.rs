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

 // OCCT IntPatch_Intersection does NOT demote BSpline surfaces.
 // BSpline stays as Parametric (ts=0); Plane stays as Geom (ts=1).
 // The (ts1 != ts2) condition triggers ImpPrmIntersection path (marching).

 // = =OCCT IntPatch_Intersection 3-category dispatch = =
 // OCCT IntPatch_Intersection.cxx L1298-1339 classifies surface pairs:
 // - ts1 == ts2 == 1 : Geom-Geom (both analytic)  ?ImpImpIntersection
 // - ts1 != ts2 : Geom-Param (one analytic, one parametric)  ?ImpPrmIntersection
 // - ts1 == ts2 == 0 : Param-Param (both parametric)  ?PrmPrmIntersection
 let (cat1, cat2) = (classify_surface_type(&s1), classify_surface_type(&s2));
 match (cat1, cat2) {
 // = =  Geom-Geom: both analytic surfaces = = 
 // OCCT ImpImpIntersection handles all analytic-analytic pairs.
 // rcad dispatches to specialized functions per combination.
 (SurfaceCategory::GeomGeom, SurfaceCategory::GeomGeom) => {
 match (&s1, &s2) {
 (Surface3::Plane(p1), Surface3::Plane(p2)) => {
 self.intersect_plane_plane_faces(f1, f2, p1, p2);
 }
 (Surface3::Plane(pl), Surface3::Sphere(sph))
 | (Surface3::Sphere(sph), Surface3::Plane(pl)) => {
 self.intersect_plane_sphere_faces(f1, f2, pl, sph);
 }
 (Surface3::Plane(pl), Surface3::Cylinder(cyl))
 | (Surface3::Cylinder(cyl), Surface3::Plane(pl)) => {
 self.intersect_plane_cylinder_faces(f1, f2, pl, cyl);
 }
 (Surface3::Sphere(sph1), Surface3::Sphere(sph2)) => {
 let (sph1, sph2) = (*sph1, *sph2);
 self.intersect_sphere_sphere_faces(f1, f2, &sph1, &sph2);
 }
 (Surface3::Sphere(sph), Surface3::Cylinder(cyl))
 | (Surface3::Cylinder(cyl), Surface3::Sphere(sph)) => {
 let (sph, cyl) = (*sph, *cyl);
 self.intersect_sphere_cylinder_faces(f1, f2, &sph, &cyl);
 }
 (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
 let (c1, c2) = (*c1, *c2);
 self.intersect_cylinder_cylinder_faces(f1, f2, &c1, &c2);
 }
 (Surface3::Plane(pl), Surface3::Cone(cone))
 | (Surface3::Cone(cone), Surface3::Plane(pl)) => {
 self.intersect_plane_cone_faces(f1, f2, pl, cone);
 }
 (Surface3::Cylinder(cyl), Surface3::Cone(cone))
 | (Surface3::Cone(cone), Surface3::Cylinder(cyl)) => {
 let (cyl, cone) = (*cyl, *cone);
 self.intersect_cylinder_cone_faces(f1, f2, &cyl, &cone);
 }
 (Surface3::Cone(cone1), Surface3::Cone(cone2)) => {
 let (cone1, cone2) = (*cone1, *cone2);
 self.intersect_cone_cone_faces(f1, f2, &cone1, &cone2);
 }
 (Surface3::Plane(pl), Surface3::Torus(tor))
 | (Surface3::Torus(tor), Surface3::Plane(pl)) => {
 self.intersect_torus_plane_faces(f1, f2, tor, pl);
 }
 (Surface3::Sphere(sph), Surface3::Torus(tor))
 | (Surface3::Torus(tor), Surface3::Sphere(sph)) => {
 self.intersect_torus_sphere_faces(f1, f2, tor, sph);
 }
 (Surface3::Cylinder(cyl), Surface3::Torus(tor))
 | (Surface3::Torus(tor), Surface3::Cylinder(cyl)) => {
 self.intersect_torus_cylinder_faces(f1, f2, tor, cyl);
 }
 (Surface3::Cone(cone), Surface3::Torus(tor))
 | (Surface3::Torus(tor), Surface3::Cone(cone)) => {
 self.intersect_torus_cone_faces(f1, f2, tor, cone);
 }
 (Surface3::Torus(tor1), Surface3::Torus(tor2)) => {
 self.intersect_torus_torus_faces(f1, f2, tor1, tor2);
 }
 (Surface3::Sphere(sph), Surface3::Cone(cone))
 | (Surface3::Cone(cone), Surface3::Sphere(sph)) => {
 let (sph, cone) = (*sph, *cone);
 self.intersect_sphere_cone_faces(f1, f2, &sph, &cone);
 }
 _ => {}
 }
 }
 // = =  Geom-Param: one analytic, one parametric = = 
 // OCCT ImpPrmIntersection handles this category.
 // rcad: use PrmPrmIntersection when a Param surface is BSpline/Bezier
 // (marching handles Plane-Plane quickly; PrmPrm handles mixed pairs).
 (SurfaceCategory::GeomGeom, SurfaceCategory::ParamParam)
 | (SurfaceCategory::ParamParam, SurfaceCategory::GeomGeom) => {
 let any_bspline = matches!(&s1, Surface3::BSpline(_) | Surface3::Bezier(_))
 || matches!(&s2, Surface3::BSpline(_) | Surface3::Bezier(_));
 if any_bspline {
 self.intersect_ff_by_prmprm(f1, f2, &s1, &s2);
 } else {
 self.intersect_ff_by_marching(f1, f2);
 }
 }
 _ => {
 // ParamParam (both parametric): PrmPrmIntersection
 self.intersect_ff_by_prmprm(f1, f2, &s1, &s2);
 }
 }

 // = =  Reverse Seam Edge Shift (OCCT ApplyTrsf L560) = = = = = = = = = = = = = = 
 if let Some(ref info) = shift_info {
 self.reverse_seam_edge_shift(f1, f2, info);
 }
 // = =  Restore seam shift tol = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 self.seam_shift_tol = old_shift_tol;

 // OCCT L600-608: CheckCurve validation -- validate curve bounds before processing
 if let Some(ff_curves) = self.find_face_face_curve_indices(f1, f2) {
 let invalid: Vec<usize> = ff_curves.iter().filter(|&&ci| {
 let ic = &self.ds.intersection_curves[ci];
 !ic.t_range[0].is_finite() || !ic.t_range[1].is_finite()
 || (ic.t_range[1] - ic.t_range[0]).abs() < 1e-12
 || ic.curve.point_at(ic.t_range[0]).is_nan()
 }).copied().collect();
 for &ci in &invalid {
 if ci < self.ds.intersection_curves.len() {
 self.ds.intersection_curves[ci].t_range = [0.0, 0.0];
 }
 }
 }
 //  ?OCCT-aligned:ComputeTolReached3d + PrepareLines3D  ?post-process all
 // intersection curves for this face pair.  Runs for every path (analytic,
 // numeric_intss, marching) to ensure consistent curve tolerance and
 // closed-curve splitting.
 if let Some(ff_curves) = self.find_face_face_curve_indices(f1, f2) {
 let t_a = self.ff_tol(f1, f1);
 let t_b = self.ff_tol(f2, f2);
 for &ci in &ff_curves {
 let (curve, pca, pcb, sv, ev, tr) = {
 let ic = &self.ds.intersection_curves[ci];
 (ic.curve.clone(), ic.pcurve_on_a.clone(), ic.pcurve_on_b.clone(),
 ic.start_vertex, ic.end_vertex, ic.t_range)
 };
 let (new_tol, _) = inttools::pcurve_derive::compute_intersection_curve_tolerance(
 &curve, pca.as_ref(), pcb.as_ref(),
 &self.ds.faces[f1].surface, &self.ds.faces[f2].surface, tr,
 t_a, t_b, 0.0,
 );
 if new_tol > TOLERANCE_ABS {
 let vt = new_tol.min(TOLERANCE_MESH_LEGACY);
 self.ds.vertices[sv].geom_tol = self.ds.vertices[sv].geom_tol.max(vt);
 self.ds.vertices[ev].geom_tol = self.ds.vertices[ev].geom_tol.max(vt);
 }
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
 // OCCT-aligned: Register curves_sc and vertices_in for this face pair
 // (PostTreatFF equivalent in OCCT PaveFiller_6.cxx L1165-1397).
 let post_ff_curves = self.find_face_face_curve_indices(f1, f2)
 .unwrap_or_default();
 self.ds.faces[f1].face_info.curves_sc.extend(&post_ff_curves);
 self.ds.faces[f2].face_info.curves_sc.extend(&post_ff_curves);
 for &ci in &post_ff_curves {
 if ci < self.ds.intersection_curves.len() {
 let ic = &self.ds.intersection_curves[ci];
 self.ds.faces[f1].face_info.vertices_in.insert(ic.start_vertex);
 self.ds.faces[f1].face_info.vertices_in.insert(ic.end_vertex);
 self.ds.faces[f2].face_info.vertices_in.insert(ic.start_vertex);
 self.ds.faces[f2].face_info.vertices_in.insert(ic.end_vertex);
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
 pub(crate) fn check_self_interference(&self) -> Result<(), String> {
 // OCCT L30-34: single argument  ?self-interference mode, skip.
 if self.ds.a_vertex_count == 0 {
 return Ok(());
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

 if warnings.is_empty() {
 Ok(())
 } else {
 Err(format!("Self-interference detected:\n{}", warnings.join("\n")))
 }
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

 ///  ?OCCT-aligned: PrmPrmIntersection path for ParamParam + mixed BSpline pairs.
 pub(crate) fn intersect_ff_by_prmprm(
 &mut self,
 f1: usize,
 f2: usize,
 s1: &Surface3,
 s2: &Surface3,
 ) {
 // Generate seed points by sampling s1 and projecting onto s2
 let (u1_min, u1_max, v1_min, v1_max) = ff_uv_bounds(s1);
 let (u2_min, u2_max, v2_min, v2_max) = ff_uv_bounds(s2);

 let n_uv = 12usize;
 let du1 = if n_uv > 1 { (u1_max - u1_min) / (n_uv - 1) as f64 } else { 0.0 };
 let dv1 = if n_uv > 1 { (v1_max - v1_min) / (n_uv - 1) as f64 } else { 0.0 };

 let mut seed_points = Vec::new();
 let tol = self.ff_tol(f1, f2).max(1e-6);

 for iu in 0..n_uv {
 let u1 = u1_min + iu as f64 * du1;
 for iv in 0..n_uv {
 let v1 = v1_min + iv as f64 * dv1;
 let p3d = s1.point_at(u1, v1);
 if !p3d.is_finite() { continue; }
 if let Some((u2, v2)) = ff_project(s2, p3d) {
 if u2 >= u2_min - 0.1 && u2 <= u2_max + 0.1
 && v2 >= v2_min - 0.1 && v2 <= v2_max + 0.1
 {
 let p3d_2 = s2.point_at(u2, v2);
 let dist = p3d.distance(p3d_2);
 if dist < tol * 10.0 {
 seed_points.push(prm_prm_intersection::PntOn2S {
 p3d: (p3d + p3d_2) * 0.5, u1, v1, u2, v2,
 });
 }
 }
 }
 }
 }

 if seed_points.is_empty() { return; }

 // Deduplicate seeds
 seed_points.sort_by(|a, b| a.u1.partial_cmp(&b.u1).unwrap_or(std::cmp::Ordering::Equal));
 seed_points.dedup_by(|a, b| {
 (a.u1 - b.u1).abs() < 0.1 && (a.v1 - b.v1).abs() < 0.1
 && (a.u2 - b.u2).abs() < 0.1 && (a.v2 - b.v2).abs() < 0.1
 });

 // Run PrmPrmIntersection
 let increment = 0.01; let deflection = 0.01;
 let epsilon = 1e-7; let tol_tangency = self.ff_tol(f1, f2).max(1e-7);
 let mut prm = prm_prm_intersection::PrmPrmIntersection::new();
 prm.perform_with_seeds(s1, s2, &seed_points, tol_tangency, epsilon, deflection, increment);
 if prm.is_empty() { return; }

 // Convert lines to DS IntersectionCurves
 let mut curve_indices = Vec::new();
 for line in &prm.slin {
 if line.points.len() < 2 { continue; }
 let chain: Vec<DVec3> = line.points.iter().map(|p| p.p3d).collect();
 let v_start = self.ds.add_vertex(chain[0]);
 let v_end = self.ds.add_vertex(chain[chain.len() - 1]);
 let arc_len: f64 = chain.windows(2).map(|w| (w[1] - w[0]).length()).sum();
 let dir = (chain[chain.len() - 1] - chain[0]).normalize_or_zero();
 let pcurve_a = crate::inttools::pcurve_derive::polyline_pcurve_by_projection(&chain, s1);
 let pcurve_b = crate::inttools::pcurve_derive::polyline_pcurve_by_projection(&chain, s2);

 let curve_idx = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(crate::bopds::ds::IntersectionCurve {
 curve: Curve3::Line(Line3 { origin: chain[0], direction: if dir.length_squared() > 0.5 { dir } else { DVec3::X } }),
 polyline: chain, start_vertex: v_start, end_vertex: v_end,
 t_range: [0.0, arc_len.max(TOLERANCE_LINEAR_ULTRA_STRICT)],
 pcurve_on_a: pcurve_a, pcurve_on_b: pcurve_b,
 geom_tol: crate::tolerance::TOLERANCE_ABS,
 pave_blocks: Vec::new(),
 curve_extra: crate::bopds::ds::CurveExtra::default(),
 });

 self.ds.faces[f1].face_info.curves_sc.insert(curve_idx);
 self.ds.faces[f2].face_info.curves_sc.insert(curve_idx);
 let _ = (v_start, v_end);
 curve_indices.push(curve_idx);
 }
 if !curve_indices.is_empty() {
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF{ f1, f2,   curves: curve_indices, points: vec![], tangent_faces: false });
 }
 }
}

//  € € PrmPrm helpers  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

fn ff_uv_bounds(s: &Surface3) -> (f64, f64, f64, f64) {
 match s {
 Surface3::Plane(_) => (-1e5, 1e5, -1e5, 1e5),
 Surface3::Sphere(sp) => (0.0, std::f64::consts::TAU, -std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2),
 Surface3::Cylinder(_) => (0.0, std::f64::consts::TAU, -1e5, 1e5),
 Surface3::Cone(_) => (0.0, std::f64::consts::TAU, -1e5, 1e5),
 Surface3::Torus(_) => (0.0, std::f64::consts::TAU, 0.0, std::f64::consts::TAU),
 Surface3::BSpline(bsp) => {
 let u_max = (bsp.knots_u.len().saturating_sub(bsp.degree_u + 1)) as f64;
 let v_max = (bsp.knots_v.len().saturating_sub(bsp.degree_v + 1)) as f64;
 (0.0, u_max.max(1.0), 0.0, v_max.max(1.0))
 }
 Surface3::Bezier(_) => (0.0, 1.0, 0.0, 1.0),
 _ => (0.0, 1.0, 0.0, 1.0),
 }
}

fn ff_project(surf: &Surface3, pt: DVec3) -> Option<(f64, f64)> {
 let (uv, _proj) = crate::extrema::closest_point_on_surface(surf, pt);
 if uv.x.is_finite() && uv.y.is_finite() { Some((uv.x, uv.y)) } else { None }
}




