use super::*;
use crate::inttools::int_patch_type::IntPatchIType;

/// Work item for Phase 2 of PerformFF (OCCT BOPAlgo_PaveFiller_6.cxx L365-518).
/// Each FF pair's setup data, computed in Phase 1, consumed in Phase 2.
struct FFWork {
 f1: usize,
 f2: usize,
 shift_info: Option<super::SeamEdgeShift>,
}

impl<'a> super::PaveFiller<'a> {
 pub(crate) fn perform_ff(&mut self) {
  // OCCT BOPAlgo_PaveFiller_6.cxx L286-623: PerformFF

  // ===== Phase 0: pairs + UpdateFaceInfo (OCCT L291-320) =====

  // OCCT L291: myIterator->Initialize(TopAbs_FACE, TopAbs_FACE)
  let pairs = self.ff_candidate_pairs();
  let i_size = pairs.len();

  // OCCT L295-302: collect touched faces from intersection pairs
  let mut a_mi_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
  for &(n_f1, n_f2) in &pairs {
    a_mi_fence.insert(n_f1);
    a_mi_fence.insert(n_f2);
  }

  // OCCT L304-311: collect rest of touched faces (those with HasReference)
  for fi in 0..self.ds.faces.len() {
    if fi < self.ds.shape_info.len() && self.ds.shape_info[fi].has_reference() {
      a_mi_fence.insert(fi);
    }
  }

  // OCCT L313-314: UpdateFaceInfoOn/In
  for &fi in &a_mi_fence {
    self.ds.refine_face_info_on(fi);
    self.ds.refine_face_info_in(fi);
  }

  // OCCT L316-320: early return if no intersection pairs
  if i_size == 0 {
    return;
  }

  // ===== Phase 1 (OCCT L365-518): prepare work items =====
  // OCCT collects BOPAlgo_FaceFace objects; rcad collects FFWork structs.

  let mut a_work: Vec<FFWork> = Vec::new();

  for &(n_f1, n_f2) in &pairs {
    if self.ds.has_interf_ff(n_f1, n_f2) {
      continue;
    }

    if !self.use_glue() {
      // ---- Non-glue path (OCCT L374-508) ----

      let is_plane1 = matches!(self.ds.faces[n_f1].surface, Surface3::Plane(_));
      let is_plane2 = matches!(self.ds.faces[n_f2].surface, Surface3::Plane(_));

      // OCCT L381-392: CheckPlanes - skip parallel planes that don't truly intersect
      if is_plane1 && is_plane2 {
        if !self.check_planes(n_f1, n_f2) {
          // OCCT L387-391: register empty FF interference for non-intersecting planes
          self.ds.interf_ff.push(InterferenceFF {
            f1: n_f1, f2: n_f2,
            curves: Vec::new(),
            points: Vec::new(),
            tangent_faces: false,
          });
          continue;
        }
      }

      // OCCT L394-487: seam edge shift computation (not applied until Phase 2)
      let shift_info = self.check_seam_edge_shift(n_f1, n_f2);

      // OCCT L499-504: GetEFPnts — collect EF intersection points as seeds
      // for IntPatch_Intersection. (rcad: stored on work item for Phase 2.)

      a_work.push(FFWork { f1: n_f1, f2: n_f2, shift_info });
    } else {
      // ---- Glue mode (OCCT L511-517) ----
      self.ds.interf_ff.push(InterferenceFF {
        f1: n_f1, f2: n_f2,
        curves: Vec::new(),
        points: Vec::new(),
        tangent_faces: false,
      });
    }
  }

  // ===== Phase 2 (OCCT L538-622): execute + process results =====
  // OCCT: parallel Perform on all BOPAlgo_FaceFace objects, then process.
  // rcad: sequential for now.


  for work in a_work {
    if self.check_stop("PerformFF") { return; }


    // OCCT L546-554: BOPAlgo_FaceFace::Perform (IntPatch + MakeCurve).
    // intersect_face_face handles: seam shift application, IntPatch, MakeCurve,
    // PutPointsOnLine, PrepareLines3D, ApplyTrsf, ComputeTolReached3d, PostTreatFF.
    self.intersect_face_face(work.f1, work.f2);

    // OCCT L600-601: IntTools_Tools::CheckCurve — discard curves too thin
    // to form valid edges (size < 3*Precision::Confusion()).
    if let Some(ff_entry) = self.ds.interf_ff.last_mut() {
      if ff_entry.f1 == work.f1 && ff_entry.f2 == work.f2 {
        ff_entry.curves.retain(|&ci| {
          if ci >= self.ds.intersection_curves.len() { return false; }
          let ic = &self.ds.intersection_curves[ci];
          // OCCT L558-572: IntTools_Tools::CheckCurve — build 3D bounding box
          // via BndLib_Add3dCurve::Add, then check !IsThin(3*Precision::Confusion()).
          let tol_cmp = 3.0 * crate::TOLERANCE_ABS;
          let mut bb_min = DVec3::splat(f64::INFINITY);
          let mut bb_max = DVec3::splat(f64::NEG_INFINITY);
          let n_samples = 100usize.max(2);
          for si in 0..n_samples {
            let t = ic.t_range[0] + (ic.t_range[1] - ic.t_range[0]) * (si as f64) / ((n_samples - 1) as f64);
            let p = ic.curve.point_at(t);
            bb_min = bb_min.min(p);
            bb_max = bb_max.max(p);
          }
          // Expand box by tolerance (OCCT BndLib_Add3dCurve::Add uses max(tol, tang_tol))
          let expand = ic.geom_tol.max(ic.curve_extra.tangential_tol).max(crate::TOLERANCE_ABS);
          bb_min -= DVec3::splat(expand);
          bb_max += DVec3::splat(expand);
          // IsThin check: reject if box is thin in ALL three directions
          let dx = bb_max.x - bb_min.x;
          let dy = bb_max.y - bb_min.y;
          let dz = bb_max.z - bb_min.z;
          let valid = dx > tol_cmp || dy > tol_cmp || dz > tol_cmp;
          if !valid && std::env::var("RCAD_DBG_FF").is_ok() {
            eprintln!("[FF] CheckCurve discard IC[{}] box=({:.2e},{:.2e},{:.2e}) tol={:.2e}", ci, dx, dy, dz, tol_cmp);
          }
          valid
        });
        // OCCT L568-570: myDS->AddInterf(nF1, nF2) — mark FF pair as processed
        // if it has valid curves or points.
        if !ff_entry.curves.is_empty() || !ff_entry.points.is_empty() {
          self.ds.try_add_interf(work.f1, work.f2);
        }
      }
    }
  }

  // Dedup FF interferences by pair (BOPDS_IndexRange equivalent)
  self.ds.dedup_ff_interferences();
 }

 pub(crate) fn ff_candidate_pairs(&self) -> Vec<(usize, usize)> {
  // Build face BVH locally (OCCT builds BOPTools_BoxTree inside PerformFF)
  // Equivalent to BOPDS_Iterator::Initialize(TopAbs_FACE, TopAbs_FACE)
  let face_bvh = {
      let mut indices = Vec::new();
      let mut aabbs = Vec::new();
      for (fi, _f) in self.ds.faces.iter().enumerate() {
          indices.push(fi);
          aabbs.push(crate::bopds::ds::face_aabb::face_aabb(self.ds, fi));
      }
      if indices.len() >= 20 { Some(crate::bvh::DsBvh::build(indices, aabbs)) } else { None }
  };

  // Get FF pair candidates from BVH or pair iterator.
  // OCCT equivalent: myIterator->Initialize(TopAbs_FACE, TopAbs_FACE) loop.
  if let Some(ref fbvh) = face_bvh {
    let candidates = crate::bvh::DsBvh::candidate_pairs(fbvh, fbvh);
    candidates
      .into_iter()
      .filter(|&(fa, fb)| self.ds.face_origin(fa) != self.ds.face_origin(fb))
      .collect()
  } else {
    let a_fcount = self.ds.a_face_count;
    let mut result = Vec::new();
    let mut fit = crate::bopds::ds::PairIterator::prepare_ab(a_fcount, self.ds.faces.len());
    while fit.more() {
      let pk = fit.value();
      result.push((pk.i1, pk.i2));
      fit.next();
    }
    result
  }
 }

 /// OCCT PaveFiller_6.cxx L419-486: seam edge shift using aEEMap.
 /// Checks EE intersections between seam (closed) edges of a face pair
 /// and returns the shift needed to align them, or None.


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

 // the seam edge shift is a SMALL tolerance
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
 self.ds.vertex_data_mut(sv).point += inv_vec;
 }
 if ev < self.ds.vertices.len() {
 self.ds.vertex_data_mut(ev).point += inv_vec;
 }
 }
 }

 /// OCCT: dispatch FF intersection by surface type
 /// OCCT L344-608: IntTools_FaceFace::Perform  -- face-face intersection.
/// Dispatches by surface type with bReverse sorting (OCCT SortTypes/IndexType),
/// then runs intersection, MakeCurve, ComputeTolReached3d, PrepareLines3D,
/// and point/curve registration.
pub(crate) fn intersect_face_face(&mut self, f1: usize, f2: usize) {
 let dbg_ff = std::env::var("RCAD_DBG_FF").is_ok();
 if dbg_ff { eprintln!("[FF] intersect_face_face: f1={} f2={}", f1, f2); }
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

 // OCCT L351-375: SortTypes  -- canonical surface ordering.
 // Swap f1/f2 so the higher-type surface is always "face A".
 let type_idx1 = Self::surface_type_index(&s1);
 let type_idx2 = Self::surface_type_index(&s2);
 let b_reverse = type_idx1 < type_idx2;
 // OCCT L354: if bReverse, swap face refs so myFace1 gets the higher type.
 let (f1, f2, s_a, s_b) = if b_reverse { (f2, f1, &s2, &s1) } else { (f1, f2, &s1, &s2) };

 // OCCT L384-393: tolerance setup
 // OCCT: myTolF1 = BRep_Tool::Tolerance(myFace1) + aFuzz, etc.
 // rcad: tolerance handled by PaveFiller's fuzzy_tolerance and face geom_tols.
// Compute TolFF = max(face tolerances) per OCCT ToleranceFF (BOPAlgo_PaveFiller_6.cxx L3918-3942).
let tol1 = self.ds.faces.get(f1).map_or(1e-7, |f| f.geom_tol);
let tol2 = self.ds.faces.get(f2).map_or(1e-7, |f| f.geom_tol);
let mut a_tol_ff = tol1.max(tol2);
fn is_analytic_ff(surf: &Surface3) -> bool {
    matches!(surf, Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_))
}
if !is_analytic_ff(&self.ds.faces[f1].surface) || !is_analytic_ff(&self.ds.faces[f2].surface) {
    a_tol_ff = a_tol_ff.max(5e-6);
}
// Ensure minimum tolerance for IntPatch to work
a_tol_ff = a_tol_ff.max(1e-7);
if dbg_ff { eprintln!("[FF] ToleranceFF: f1={} tol={:.2e} f2={} tol={:.2e} -> a_tol_ff={:.2e}", f1, tol1, f2, tol2, a_tol_ff); }

 // OCCT L395-401: isFace1Quad/isFace2Quad  -- skip; rcad uses IntPatchIntersection
 // which dispatches by quad type internally.

 //  OCCT L404-434: Plane-Plane fast path (PerformPlanes) 
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

 // IntPatch_Intersection: generic surface-surface intersection.
 let mut int_patch = crate::inttools::int_patch_intersection::IntPatchIntersection::new();
 int_patch.perform(s_a, s_b, a_tol_ff, a_tol_ff);
 if int_patch.tangent_faces() {
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
 f1, f2, curves: Vec::new(), points: Vec::new(), tangent_faces: true,
 });
 if let Some(ref info) = shift_info { self.reverse_seam_edge_shift(f1, f2, info); }
 self.seam_shift_tol = old_shift_tol;
 return;
 }

 // PutPointsOnLine (IntPatch_Intersection.cxx L268-312).
 // Projects intersection points onto each analytic line to create
 // boundary-crossing vertices.  These vertices split the line into
 // valid intervals for MakeCurve/TreatCircle.
 for li in 0..int_patch.nb_lines() {
 self.put_points_on_line(f1, f2, int_patch.line_mut(li));
 }

 // OCCT L498-504: GetEFPnts → SetList passes EF points to IntPatch's PutPointsOnLine.
 // rcad: IntPatch skips PutPointsOnLine; EF=0 for sphere-sphere (PerformEF gap).
 // EF projection here would require EF>0. Currently EF=0, so no points to project.

 // MakeCurve (IntTools_FaceFace.cxx L695-1846) for each IntPatch line.
 // Returns a Vec of IntersectionCurve  -- one per valid part from the
 // LineConstructor (OCCT supports aNbParts > 1, e.g. multi-segment clipping).
 let mut ff_curve_indices: Vec<usize> = Vec::new();
 for i in 0..int_patch.nb_lines() {
 let ics = self.make_intersection_curve(f1, f2, int_patch.line(i));
 for ic in ics {
   let ci = self.ds.intersection_curves.len();
   let mut adjusted_ic = ic;
   // OCCT L558-567: if reversed, swap pcurves (first  -- second).
   if b_reverse {
   std::mem::swap(&mut adjusted_ic.pcurve_on_a, &mut adjusted_ic.pcurve_on_b);
   }
   self.ds.intersection_curves.push(adjusted_ic);
   // OCCT: vertices created after MakeCurve (in Process/PerformFF).
   // BRepBuilderAPI_MakeVertex(P3D) + myDS->Index for each endpoint.
   let (p_start, p_end, t0, t1) = {
     let ic_ref = &self.ds.intersection_curves[ci];
     let t0 = ic_ref.t_range[0];
     let t1 = ic_ref.t_range[1];
     (ic_ref.curve.point_at(t0), ic_ref.curve.point_at(t1), t0, t1)
   };
   let sv = if p_start.is_finite() { self.ds.add_vertex(p_start) } else { usize::MAX };
   let ev = if p_end.is_finite() { self.ds.add_vertex(p_end) } else { usize::MAX };
   self.ds.intersection_curves[ci].start_vertex = sv;
   self.ds.intersection_curves[ci].end_vertex = ev;
   // OCCT: init IC pave_blocks so MakeSplitEdges can create section edges
   {
    use crate::bopds::pave::{Pave, PaveBlock, SharedPB};
    let sv = self.ds.intersection_curves[ci].start_vertex;
    let ev = self.ds.intersection_curves[ci].end_vertex;
    let t0 = self.ds.intersection_curves[ci].t_range[0];
    let t1 = self.ds.intersection_curves[ci].t_range[1];
    let pb = PaveBlock::new(0, Pave{vertex_idx:sv, param:t0}, Pave{vertex_idx:ev, param:t1});
    let spb = SharedPB::new(pb);
    self.ds.intersection_curves[ci].pave_blocks.push(spb.clone());
    self.ds.pave_blocks.push(spb);
   }
   if std::env::var("RCAD_DBG_MB").is_ok() {
    let ic2 = &self.ds.intersection_curves[ci];
    eprintln!("[DBG_IC3] PUSHED IC[{}]: geom_tol={:.6e} t_range={:.6} {:.6}", ci, ic2.geom_tol, ic2.t_range[0], ic2.t_range[1]);
   }
   ff_curve_indices.push(ci);
 }
 }

 // OCCT L576-608: points  -- filter by isPointInOnFace, append to myPnts.
 let mut ff_point_indices: Vec<crate::bopds::ds::types::FFPoint> = Vec::new();
 for pi in 0..int_patch.nb_points() {
 let pt = int_patch.point(pi);
 let (uv_a, uv_b, f_a, f_b) = if b_reverse {
   (glam::DVec2::new(pt.u2, pt.v2), glam::DVec2::new(pt.u1, pt.v1), f2, f1)
 } else {
   (glam::DVec2::new(pt.u1, pt.v1), glam::DVec2::new(pt.u2, pt.v2), f1, f2)
 };
 if !self.context.is_point_in_on_face(self.ds, f_a, uv_a) { continue; }
 if !self.context.is_point_in_on_face(self.ds, f_b, uv_b) { continue; }
 // FFPoint stores point data inline (OCCT BOPDS_Point). No DS vertex created yet.
 ff_point_indices.push(crate::bopds::ds::types::FFPoint::new(pt.p1, uv_a, uv_b));
 }
 if std::env::var("RCAD_DBG_FF").is_ok() { eprintln!("[FF]   -> curves={} nLines={}", ff_curve_indices.len(), int_patch.nb_lines()); }
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
 f1, f2, curves: ff_curve_indices, points: ff_point_indices, tangent_faces: false,
 });

 // = =  Reverse Seam Edge Shift (OCCT ApplyTrsf L560) = = = = = = = = = = = = = = 
 if let Some(ref info) = shift_info {
 self.reverse_seam_edge_shift(f1, f2, info);
 }
 // = =  Restore seam shift tol = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 self.seam_shift_tol = old_shift_tol;

 // ComputeTolReached3d + PrepareLines3D.
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
 inttools::pcurve_derive::prepare_lines_3d(&mut self.ds.intersection_curves, false);
 // After PrepareLines3D splits closed curves, the split
 // segments are added to the same FF interference entry.  Update the
 // FF entry's curve list to include any newly created curve indices.
 if n_curves_before_split != self.ds.intersection_curves.len() {
  if let Some(ff_entry) = self.ds.interf_ff.last_mut() {
   for new_ci in n_curves_before_split..self.ds.intersection_curves.len() {
    ff_entry.curves.push(new_ci);
   }
  }
 }
 //  After PrepareLines3D splits closed curves, new curve endpoints
 // must be updated to the split points. For start==end but
 // non-full-period t_range (i.e. split half-circle), compute
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
   half_circle && ic.start_vertex != usize::MAX && ic.start_vertex == ic.end_vertex
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
  self.ds.push_vertex(DSVertex { point: p_start, geom_tol: TOLERANCE_ABS, origin: None, is_internal: true, location: 0 }, None);
  let v_end = self.ds.vertices.len();
  self.ds.push_vertex(DSVertex { point: p_end, geom_tol: TOLERANCE_ABS, origin: None, is_internal: true, location: 0 }, None);
  self.ds.intersection_curves[ci].start_vertex = v_start;
  self.ds.intersection_curves[ci].end_vertex = v_end;
 }
 }
 // PreparePostTreatFF (PaveFiller_6.cxx L3642-3668).
 let post_ff_curves = self.find_face_face_curve_indices(f1, f2)
 .unwrap_or_default();
 self.ds.face_info_mut(f1).curves_sc.extend(&post_ff_curves);
 self.ds.face_info_mut(f2).curves_sc.extend(&post_ff_curves);
 for &ci in &post_ff_curves {
 if ci < self.ds.intersection_curves.len() {
   let ic = &self.ds.intersection_curves[ci];
   let sv = ic.start_vertex;
   let ev = ic.end_vertex;
   if sv < self.ds.vertices.len() {
   self.ds.face_info_mut(f1).vertices_in.insert(sv);
   self.ds.face_info_mut(f2).vertices_in.insert(sv);
   }
   if ev < self.ds.vertices.len() {
   self.ds.face_info_mut(f1).vertices_in.insert(ev);
   self.ds.face_info_mut(f2).vertices_in.insert(ev);
   }
 }
 }
  } // if let Some(ff_curves)
} // fn intersect_face_face

/// IndexType (IntTools_FaceFace.cxx L2844-2870).
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

/// OCCT L2426-2560: PerformPlanes  -- plane-plane intersection fast path.
fn perform_plane_plane(&mut self, f1: usize, f2: usize) {
 use rcad_kernel::geom::{Curve3, Surface3};
 let pln1 = match &self.ds.faces[f1].surface { Surface3::Plane(p) => p, _ => return };
 let pln2 = match &self.ds.faces[f2].surface { Surface3::Plane(p) => p, _ => return };
 let mut geo = crate::inttools::int_ana_quad_quad_geo::QuadQuadGeo::new();
 let q1 = crate::inttools::int_surf_quadric::Quadric::from_plane(pln1);
 let q2 = crate::inttools::int_surf_quadric::Quadric::from_plane(pln2);
 geo.perform_plane_plane(&q1, &q2, 1e-8, self.fuzzy_tolerance);
 if !geo.is_done() { return; }
 use crate::inttools::int_ana_quad_quad_geo::AnaResultType;
 if let AnaResultType::Same = geo.type_inter() {
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
   f1, f2, curves: Vec::new(), points: Vec::new(), tangent_faces: true,
 });
 return;
 }
 if matches!(geo.type_inter(), AnaResultType::Empty) { return; }
 let line3 = geo.line(1);
 let line3d = Curve3::Line(line3.clone());
 let pcurve1 = crate::inttools::pcurve_derive::line_pcurve_on_plane(&line3, pln1);
 let pcurve2 = crate::inttools::pcurve_derive::line_pcurve_on_plane(&line3, pln2);
 // OCCT L2514: new Geom_TrimmedCurve(aGLin, pmin, pmax)
 // OCCT L2521: new Geom2d_TrimmedCurve(C2d, pmin, pmax)
 let uv1 = self.context.uv_bounds(self.ds, f1);
 let uv2 = self.context.uv_bounds(self.ds, f2);
 let tol = self.ds.face_tolerance(f1).max(self.ds.face_tolerance(f2));
 let p1 = crate::inttools::classify_lin2d::classify_lin2d(&pcurve1, uv1, tol);
 let p2 = crate::inttools::classify_lin2d::classify_lin2d(&pcurve2, uv2, tol);
 let (Some([p11, p12]), Some([p21, p22])) = (p1, p2) else { return };
 if p21 >= p12 || p22 <= p11 { return; }
 let pmin = p11.max(p21);
 let pmax = p12.min(p22);
 if pmax - pmin <= tol { return; }
 let t_range = [pmin, pmax];
 let mut curve_extra = crate::bopds::ds::CurveExtra::default();
 curve_extra.tangential_tol = tol;
 // OCCT L2514: new Geom_TrimmedCurve(aGLin, pmin, pmax)
 // OCCT L2521: new Geom2d_TrimmedCurve(C2d, pmin, pmax)
 let trimmed_curve = Curve3::Trimmed(Box::new(TrimmedCurve3::new(line3d.clone(), pmin, pmax)));
 let trimmed_pca = Some(Curve2d::Trimmed(TrimmedCurve2 { curve: Box::new(pcurve1), t_min: pmin, t_max: pmax }));
 let trimmed_pcb = Some(Curve2d::Trimmed(TrimmedCurve2 { curve: Box::new(pcurve2), t_min: pmin, t_max: pmax }));
 let ic = crate::bopds::ds::IntersectionCurve {
 curve: trimmed_curve, polyline: Vec::new(),
 start_vertex: usize::MAX, end_vertex: usize::MAX,
 t_range,
 pcurve_on_a: trimmed_pca, pcurve_on_b: trimmed_pcb,
 geom_tol: tol, pave_blocks: Vec::new(), curve_extra,
 };
 // OCCT: vertices created by BRepBuilderAPI_MakeVertex + myDS->Index
 // after the curve is stored (matching the MakeCurve caller pattern).
 let sv = self.ds.add_vertex(ic.curve.point_at(t_range[0]));
 let ev = self.ds.add_vertex(ic.curve.point_at(t_range[1]));
 let mut ic = ic;
 ic.start_vertex = sv;
 ic.end_vertex = ev;
 let ci = self.ds.intersection_curves.len();
 self.ds.intersection_curves.push(ic);
 self.ds.interf_ff.push(crate::bopds::ds::InterferenceFF {
 f1, f2, curves: vec![ci], points: Vec::new(), tangent_faces: false,
 });
 self.ds.face_info_mut(f1).curves_sc.insert(ci);
 self.ds.face_info_mut(f2).curves_sc.insert(ci);
}



 /// MakeCurve (IntTools_FaceFace.cxx L695-1846).
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

 // ===== OCCT IntTools_FaceFace.cxx L695-751 =====
 // OCCT L700-714: local vars
 // OCCT L717: reapprox label (not needed in sequential rcad)
 // OCCT L719: Tolpc = myTolApprox
 // OCCT L720: bAvoidLineConstructor = false
 let mut b_avoid_line_constructor = false;

 // OCCT L721-722: L = myIntersector.Line(Index); typl = L->ArcType();
 let typl = line.line_type;

 // OCCT L724-744: IntPatch_Walking special handling
 if line.is_wline() {
   let nbp = line.nb_points();
   if nbp >= 2 {
     let p1 = line.point(0).p3d;
     let p2 = line.point(nbp - 1).p3d;
     // OCCT L740-743: if endpoints are nearly coincident, use LineConstructor
     if p1.distance_squared(p2) < 1e-14 {
       b_avoid_line_constructor = false;
     }
   }
 }

 // OCCT L748-751: IntPatch_Restriction — skip LineConstructor
 if typl == crate::inttools::int_patch_type::IntPatchIType::Restriction {
   b_avoid_line_constructor = true;
 }

 // OCCT L755-773: LineConstructor.Perform(L)
 // If !IsDone → return empty. If NbParts <= 0 → return empty.
 let parts: Vec<[f64; 2]> = if !b_avoid_line_constructor {
   let p = self.line_constructor_parts(
     &line.curve, line.t_range, typl, &line.vertices, f1, f2);
   if p.is_empty() { return Vec::new(); }
   p
 } else {
   // OCCT L748-750: for Restriction, skip LineConstructor, use full range
   // rcad: use the full t_range as a single part
   vec![line.t_range]
 };

 // OCCT L776-1846: switch(typl)
 match typl {
 crate::inttools::int_patch_type::IntPatchIType::Line
 | crate::inttools::int_patch_type::IntPatchIType::Parabola
 | crate::inttools::int_patch_type::IntPatchIType::Hyperbola =>
   self.make_analytic_nonperiodic_curve(f1, f2, &line.curve, &parts, typl, line.tolerance, line.tang_tolerance),
 crate::inttools::int_patch_type::IntPatchIType::Circle
 | crate::inttools::int_patch_type::IntPatchIType::Ellipse =>
   self.make_analytic_periodic_curve(f1, f2, &line.curve, line.t_range, typl, line.tolerance, line.tang_tolerance, &line.vertices),
 crate::inttools::int_patch_type::IntPatchIType::Walking
 | crate::inttools::int_patch_type::IntPatchIType::Restriction =>
   self.make_walking_curve(f1, f2, line),
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
///     both bounds finite   -- trimmed 3D curve + BuildPCurves + endpoint vertices
///     one/both infinite    -- test reference point on face domains  -- keep or reject
/// - rcad note: IntPatchLine has no vertex data, so LineConstructor returns
///   a single part with the original t_range (always infinite for lines).
fn make_analytic_nonperiodic_curve(
  &mut self, f1: usize, f2: usize,
  curve: &Curve3, parts: &Vec<[f64; 2]>, typl: IntPatchIType,
  geom_tol: f64, tang_tolerance: f64,
) -> Vec<crate::bopds::ds::IntersectionCurve> {
 use rcad_kernel::geom::Curve2dEval;
 use std::f64::consts::TAU;

 // OCCT L815-826: create analytic 3D curve from the GLine.
 // rcad: curve is already the correct analytic type in IntPatchLine.

 // OCCT L828-840: LineConstructor.Perform(L) already done upstream.
 // parts already computed by make_intersection_curve.

 if parts.is_empty() {
 return Vec::new();
 }

 let mut result = Vec::with_capacity(parts.len());

 // OCCT L842-898: per-part loop.
 for part in parts.iter() {
 let &[fprm, lprm] = part;
 let b_finite = fprm.is_finite() && lprm.is_finite() && lprm > fprm + 1e-12;

 if b_finite {
   //  Both bounds finite: trimmed curve + pcurves + vertices 
   // OCCT L835-870: Geom_TrimmedCurve + BuildPCurves + Geom2d_TrimmedCurve.
   let ic_t_range = [fprm, lprm];

   // OCCT L816-820: for Parabola, CurveTolerance(aCT3D, myTol)
   let ic_geom_tol = if typl == IntPatchIType::Parabola {
     // OCCT: IntTools_Tools::CurveTolerance(aCT3D, myTol)
     crate::boptools::curve_tolerance(curve, geom_tol)
   } else {
     geom_tol.max(crate::tolerance::TOLERANCE_ABS)
   };

   // OCCT L822-846: BuildPCurves on the trimmed range.
   // OCCT: GeomInt_IntSS::BuildPCurves(fprm, lprm, Tolpc, surface, newc, C2d)
   // OCCT L822-832: if (myApprox1) { ... }
   let pca;
   let pcb;
   // myApprox1 (always true in rcad)
   {
     let raw = self.compute_pcurve_on_surface(curve, f1);
     if raw.is_none() { continue; }
     // OCCT L832: aCurve.SetFirstCurve2d(new Geom2d_TrimmedCurve(C2d, fprm, lprm))
     pca = Some(Curve2d::Trimmed(TrimmedCurve2 {
       curve: Box::new(raw.unwrap()),
       t_min: fprm, t_max: lprm,
     }));
   }
   // myApprox2 (always true in rcad)
   {
     let raw = self.compute_pcurve_on_surface(curve, f2);
     if raw.is_none() { continue; }
     pcb = Some(Curve2d::Trimmed(TrimmedCurve2 {
       curve: Box::new(raw.unwrap()),
       t_min: fprm, t_max: lprm,
     }));
   }

   // OCCT L814-815: new Geom_TrimmedCurve(newc, fprm, lprm)
   let trimmed_curve = Curve3::Trimmed(Box::new(TrimmedCurve3::new(curve.clone(), fprm, lprm)));

   // OCCT: no vertex creation in MakeCurve (vertices created later in caller).
   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
   curve: trimmed_curve,
   polyline: Vec::new(),
   start_vertex: usize::MAX,
   end_vertex: usize::MAX,
   t_range: [fprm, lprm],
   pcurve_on_a: pca,
   pcurve_on_b: pcb,
   geom_tol: ic_geom_tol,
   pave_blocks: Vec::new(),
   curve_extra,
   });

 } else {
   //  One/both bounds infinite: test reference point 
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
   // bFNIt && bLPIt: both infinite  -- OCCT IntTools_Tools::IntermediatePoint(-dT, dT)
   crate::boptools::intermediate_point_occt(-dT, dT)
   };

   let p3d = curve.point_at(test_t);
   if !p3d.is_finite() { continue; }

   // OCCT L865-867: get surface types for the test-point branch
   let surf1 = &self.ds.faces[f1].surface;
   let surf2 = &self.ds.faces[f2].surface;
   let is_extrusion_rev_offset = |s: &Surface3| -> bool {
     matches!(s, Surface3::LinearExtrusion(_) | Surface3::Revolution(_) | Surface3::Offset(_))
   };

   // OCCT L875-882: if either surface is extrusion/offset/revolution,
   // append curve with empty pcurves (= H1, H1 in OCCT) and skip Classify.
   if is_extrusion_rev_offset(surf1) || is_extrusion_rev_offset(surf2) {
     let mut curve_extra = crate::bopds::ds::CurveExtra::default();
     curve_extra.tangential_tol = tang_tolerance;
     result.push(crate::bopds::ds::IntersectionCurve {
       curve: curve.clone(),
       polyline: Vec::new(),
       start_vertex: usize::MAX,
       end_vertex: usize::MAX,
       t_range: [fprm, lprm],
       pcurve_on_a: None,
       pcurve_on_b: None,
       geom_tol: geom_tol.max(crate::tolerance::TOLERANCE_ABS),
       pave_blocks: Vec::new(),
       curve_extra,
     });
     continue;
   }

   // OCCT L886-892: Parameters + Classify on both face domains.
   // OCCT: Tol = Precision::Confusion()
   // OCCT: Parameters(myHS1, myHS2, ptref, u1, v1, u2, v2)
   let uv1 = self.context.proj_ps(self.ds, f1, p3d);
   let uv2 = self.context.proj_ps(self.ds, f2, p3d);
   // OCCT L888: ok = (dom1->Classify(gp_Pnt2d(u1, v1), Tol) != TopAbs_OUT)
   // rcad uses fclass2d with CONFUION tolerance for 1:1 match
   let in1 = uv1.map_or(false, |(uv, _, _)| {
     let fc = self.context.fclass2d(self.ds, f1);
     fc.perform(uv, true) != crate::inttools::fclass2d::State::Out
   });
   if !in1 { continue; }
   // OCCT L890-891: if (ok) { ok = dom2->Classify(...) }
   let in2 = uv2.map_or(false, |(uv, _, _)| {
     let fc = self.context.fclass2d(self.ds, f2);
     fc.perform(uv, true) != crate::inttools::fclass2d::State::Out
   });
   if !in2 { continue; }

   // OCCT L893-896: append curve with empty pcurves (H1, H1 in OCCT).
   // rcad note: OCCT does not compute pcurves here; keeping empty to match.
   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
     curve: curve.clone(),
     polyline: Vec::new(),
     start_vertex: usize::MAX,
     end_vertex: usize::MAX,
     t_range: [fprm, lprm],
     pcurve_on_a: None,
     pcurve_on_b: None,
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
///     not full-period  -- trimmed curve + BuildPCurves(with UV bounds) + vertices
///     full-period (aNbParts=1)  -- test 18 points around circle  -- keep/reject
/// - rcad: with 0 vertices, TreatCircle returns fallback full [0, 2閿滅 interval.
fn make_analytic_periodic_curve(
  &mut self, f1: usize, f2: usize,
  curve: &Curve3, orig_t_range: [f64; 2], typl: IntPatchIType,
  geom_tol: f64, tang_tolerance: f64,
  vertices: &[crate::inttools::int_patch_line::IntPatchVertex],
) -> Vec<crate::bopds::ds::IntersectionCurve> {
 if std::env::var("RCAD_DBG_MB").is_ok() {
  eprintln!("[DBG_IC] make_analytic_periodic_curve: geom_tol={:.6e} t_range={:.6} {:.6}", geom_tol, orig_t_range[0], orig_t_range[1]);
 }
 use rcad_kernel::geom::Curve2dEval;
 use std::f64::consts::{TAU, PI};

 // OCCT L906-920: create analytic 3D curve from GLine.

 // OCCT L922-950: TreatCircle  -- split intervals with 0-crossing handling.
 // Sorts vertices on the GLine, creates intervals, tests midpoints.
 let parts = self.treat_circle_parts(curve, orig_t_range, typl, vertices, f1, f2);

 // OCCT L950-1095: aNbParts = seqp.Length() / 2.
 //   If aNbParts == 0  -- the for loop does not execute  -- no output curves.
 if parts.is_empty() {
 return Vec::new();
 }

 let aPeriod = TAU;
 let aRealEpsilon = f64::EPSILON;
 let aNbParts = parts.len();
 let mut result = Vec::with_capacity(parts.len());

 for &[fprm, lprm] in &parts {
 // OCCT L953-956: if (|fprm|>eps || |lprm-2閿滅皸>eps)  -- not full-period
 let is_full_period = fprm.abs() <= aRealEpsilon && (lprm - aPeriod).abs() <= aRealEpsilon;

 if !is_full_period && (lprm > fprm + 1e-12) {
   //  Not full-period: trimmed curve + pcurves + vertices 
   // OCCT L960-990: Geom_TrimmedCurve(newc, fprm, lprm) + BuildPCurves + append.

   // OCCT L968-990: BuildPCurves with surface UV bounds for Circle/Ellipse.
   let pca = self.compute_pcurve_on_surface(curve, f1);
   let pcb = self.compute_pcurve_on_surface(curve, f2);
   let trimmed_pca = pca.clone();
   let trimmed_pcb = pcb.clone();
   let trimmed_curve = Curve3::Trimmed(Box::new(TrimmedCurve3::new(curve.clone(), fprm, lprm)));

   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
   curve: trimmed_curve,
   polyline: Vec::new(),
   start_vertex: usize::MAX,
   end_vertex: usize::MAX,
   t_range: [fprm, lprm],
   pcurve_on_a: trimmed_pca,
   pcurve_on_b: trimmed_pcb,
   geom_tol: geom_tol.max(TOLERANCE_ABS),
   pave_blocks: Vec::new(),
   curve_extra,
   });

 } else if is_full_period && aNbParts == 1 {
   // Full-period, single part — accept full circle
   // OCCT L996-1042: trimmed full circle + BuildPCurves + append + break.
   let pca = self.compute_pcurve_on_surface(curve, f1);
   let pcb = self.compute_pcurve_on_surface(curve, f2);
   // OCCT: no vertex creation in MakeCurve (vertices created later in caller).
   let trimmed_curve = Curve3::Trimmed(Box::new(TrimmedCurve3::new(curve.clone(), fprm, lprm)));
   let mut curve_extra = crate::bopds::ds::CurveExtra::default();
   curve_extra.tangential_tol = tang_tolerance;
   result.push(crate::bopds::ds::IntersectionCurve {
   curve: trimmed_curve,
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
   //  Full-period, multiple parts: test 18 points 
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

/// LineConstructor (GeomInt_LineConstructor::Perform for GLine).
/// Iterates over IntPatch_Point vertices on the line, tests the midpoint of
/// each adjacent-vertex interval on both face domains, keeps valid intervals.
///
/// With 0 vertices (nbvtx=0): OCCT's intrvtested flag stays false, and the
/// full parameter range [FirstParameter, LastParameter] is kept as one part.
/// The caller's test-point logic decides whether to keep or reject it.
fn line_constructor_parts(
  &mut self, curve: &Curve3, orig_t_range: [f64; 2], _typl: IntPatchIType,
  vertices: &[crate::inttools::int_patch_line::IntPatchVertex],
  f1: usize, f2: usize,
) -> Vec<[f64; 2]> {
 // OCCT GeomInt_LineConstructor: split line at vertex positions,
 // test each interval's midpoint on both face domains.
 if vertices.is_empty() {
   return vec![orig_t_range];
 }
 let a_tol_pc: f64 = 1000.0 * 1e-12;
 let mut sorted: Vec<f64> = vertices.iter()
   .map(|v| v.param_on_line)
   .filter(|&t| t >= orig_t_range[0] - a_tol_pc && t <= orig_t_range[1] + a_tol_pc)
   .collect();
 sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
 // Dedup
 let mut deduped: Vec<f64> = Vec::new();
 for &p in &sorted {
   if deduped.is_empty() || (p - *deduped.last().unwrap()).abs() > a_tol_pc {
     deduped.push(p);
   }
 }
 // Add range boundaries if needed
 let t0 = if deduped.is_empty() || (deduped[0] - orig_t_range[0]).abs() > a_tol_pc {
   orig_t_range[0]
 } else {
   deduped[0]
 };
 let t1 = if deduped.is_empty() || (orig_t_range[1] - *deduped.last().unwrap()).abs() > a_tol_pc {
   orig_t_range[1]
 } else {
   *deduped.last().unwrap()
 };
 let mut all_breaks: Vec<f64> = Vec::new();
 all_breaks.push(t0);
 all_breaks.extend(deduped.iter().filter(|&&p| p > t0 + a_tol_pc && p < t1 - a_tol_pc));
 if (t1 - *all_breaks.last().unwrap()).abs() > a_tol_pc {
   all_breaks.push(t1);
 }
 // For each adjacent pair, test midpoint on both face domains
 let mut result = Vec::new();
 for i in 0..(all_breaks.len().saturating_sub(1)) {
   let b1 = all_breaks[i];
   let b2 = all_breaks[i + 1];
   if (b2 - b1).abs() <= 1e-12 { continue; }
   let t_mid = (b1 + b2) * 0.5;
   let p3d = curve.point_at(t_mid);
   if !p3d.is_finite() { continue; }
   let uv1 = self.context.proj_ps(self.ds, f1, p3d).map(|(uv, _, _)| uv);
   let uv2 = self.context.proj_ps(self.ds, f2, p3d).map(|(uv, _, _)| uv);
   let in1 = uv1.map_or(false, |uv| self.context.is_point_in_on_face(self.ds, f1, uv));
   let in2 = uv2.map_or(false, |uv| self.context.is_point_in_on_face(self.ds, f2, uv));
   if in1 && in2 {
     result.push([b1, b2]);
   }
 }
 if result.is_empty() {
   // Fallback: no valid intervals found, try the full range
   result.push(orig_t_range);
 }
 result
}

/// TreatCircle (GeomInt_LineConstructor.cxx L481-560).
/// For Circle/Ellipse with vertices: sorts vertices by parameter in [0, 2閿?,
/// creates intervals between sorted vertices, tests midpoints on both face
/// domains.  Handles 0-crossing via PeriodicReparam + SeqFprm/SeqLprm.
///
/// Without vertices (nbvtx=0): OCCT creates a zero-initialized array of size 1,
/// the sort and interval-building steps produce no valid intervals, and seqp
/// remains empty.  The caller sees aNbParts=0 and creates no output curves.
/// This function matches that behavior  -- returns empty.
/// PutPointsOnLine (IntPatch_Intersection.cxx L268-312).
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
 use crate::pave_filler::helpers::project_vertex_to_curve;
 if line.is_wline() { return; }
 let typl = line.line_type;
 if typl == IntPatchIType::Restriction { return; }
 let curve = &line.curve;
 let t_range = line.t_range;
 let a_tol_c = line.tolerance;
 for &fi in &[f1, f2] {
  if fi >= self.ds.faces.len() { continue; }
  let verts: Vec<usize> = self.ds.faces[fi].boundary_verts.clone();
  for &vi in &verts {
   if vi >= self.ds.vertices.len() { continue; }
   let v_pt = self.ds.vertex_point(vi);
   let t_opt = project_vertex_to_curve(v_pt, curve, a_tol_c);
   let t = match t_opt {
    Some(t) if t >= t_range[0] - a_tol_c && t <= t_range[1] + a_tol_c => t,
    _ => continue,
   };
   let is_dup = line.vertices.iter().any(|ev| (ev.param_on_line - t).abs() < 1e-10);
   if is_dup { continue; }
   let pca = self.compute_pcurve_on_surface(curve, f1);
   let pcb = self.compute_pcurve_on_surface(curve, f2);
   let uv1 = pca.as_ref().map(|pc| pc.point_at(t)).filter(|uv| uv.is_finite());
   let uv2 = pcb.as_ref().map(|pc| pc.point_at(t)).filter(|uv| uv.is_finite());
   let p3d = curve.point_at(t);
   line.vertices.push(crate::inttools::int_patch_line::IntPatchVertex {
    param_on_line: t, p3d,
    u1: uv1.map_or(0.0, |uv| uv.x),
    v1: uv1.map_or(0.0, |uv| uv.y),
    u2: uv2.map_or(0.0, |uv| uv.x),
    v2: uv2.map_or(0.0, |uv| uv.y),
   });
  }
 }
 line.vertices.sort_by(|a, b| a.param_on_line.partial_cmp(&b.param_on_line).unwrap_or(std::cmp::Ordering::Equal));
 // OCCT PutPointsOnLine projects parametric boundary points (seam/pole) onto the line.
 // rcad's IntPatch skips this; detect via proj_ps UV seam wrap (periodic u-jump).
 if matches!(typl, IntPatchIType::Circle | IntPatchIType::Ellipse)
   || matches!(typl, IntPatchIType::Line | IntPatchIType::Parabola | IntPatchIType::Hyperbola)
 {
   self.project_parametric_boundary(f1, f2, line, curve.clone(), t_range);
 }
}

/// OCCT PutPointsOnLine equivalent: projects each face's parametric boundary (seam)
/// onto the intersection line. Detects the u-seam wrap via proj_ps UV sampling.
/// For closed curves, adds a complementary splitting point when only one seam
/// crossing is found (mimicking EF GetEFPnts which rcad lacks).
fn project_parametric_boundary(
  &mut self, f1: usize, f2: usize,
  line: &mut crate::inttools::int_patch_line::IntPatchLine,
  curve: Curve3, t_range: [f64; 2],
) {
 let n_samples = 200usize;
 let mut candidates: Vec<f64> = Vec::new();
 for &fi in &[f1, f2] {
  if fi >= self.ds.faces.len() { continue; }
  let mut prev_uv: Option<DVec2> = None;
  for si in 0..=n_samples {
   let t = t_range[0] + (t_range[1] - t_range[0]) * (si as f64) / (n_samples as f64);
   let p3d = curve.point_at(t);
   if !p3d.is_finite() { prev_uv = None; continue; }
   let uv = self.context.proj_ps(self.ds, fi, p3d).map(|(u, _, _)| u);
   let Some(uv) = uv else { prev_uv = None; continue; };
   if let Some(puv) = prev_uv {
    // Seam crossing: UV jumps across the periodic boundary (π→-π or -π→π).
    if (uv.x - puv.x).abs() > std::f64::consts::PI {
     candidates.push(t.clamp(t_range[0], t_range[1]));
    }
   }
   prev_uv = Some(uv);
  }
 }
 candidates.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
 candidates.dedup_by(|a, b| (*a - *b).abs() < 1e-4);
 // OCCT also projects EF points onto the line (GetEFPnts). Without EF>0, add
 // diametrically opposite point so the closed curve is split into two intervals.
 if candidates.len() == 1 && (t_range[1] - t_range[0]).abs() >= std::f64::consts::TAU - 1e-8 {
  let t_opp = candidates[0] + std::f64::consts::PI;
  let t_opp = if t_opp > t_range[1] { t_opp - std::f64::consts::TAU } else { t_opp };
  candidates.push(t_opp);
 }
 candidates.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
 candidates.dedup_by(|a, b| (*a - *b).abs() < 1e-4);
 for &t in &candidates {
  if t < t_range[0] || t > t_range[1] { continue; }
  let is_dup = line.vertices.iter().any(|ev| (ev.param_on_line - t).abs() < 1e-10);
  if is_dup { continue; }
  let p3d = curve.point_at(t);
  let pca = self.compute_pcurve_on_surface(&curve, f1);
  let pcb = self.compute_pcurve_on_surface(&curve, f2);
  let uv1 = pca.as_ref().map(|pc| pc.point_at(t)).unwrap_or(DVec2::ZERO);
  let uv2 = pcb.as_ref().map(|pc| pc.point_at(t)).unwrap_or(DVec2::ZERO);
  line.vertices.push(crate::inttools::int_patch_line::IntPatchVertex {
    param_on_line: t, p3d, u1: uv1.x, v1: uv1.y, u2: uv2.x, v2: uv2.y,
  });
 }
 line.vertices.sort_by(|a, b| a.param_on_line.partial_cmp(&b.param_on_line).unwrap_or(std::cmp::Ordering::Equal));
 line.vertices.dedup_by(|a, b| (a.param_on_line - b.param_on_line).abs() < 1e-10);
}

fn treat_circle_parts(
  &mut self, curve: &Curve3, orig_t_range: [f64; 2], typl: IntPatchIType,
  vertices: &[crate::inttools::int_patch_line::IntPatchVertex],
  f1: usize, f2: usize,
) -> Vec<[f64; 2]> {
 use crate::inttools::int_patch_line::IntPatchVertex;
 use std::f64::consts::TAU;
 if std::env::var("RCAD_DBG_FF").is_ok() {
  eprintln!("[FF] treat_circle_parts: f1={} f2={} nVtx={}", f1, f2, vertices.len());
 }

  // OCCT GeomInt_LineConstructor::TreatCircle (L674-733):
 //   RejectMicroCircle, sort, RejectDuplicates, midpoint test

 // OCCT L679-681: RejectMicroCircle -- skip circles/ellipses smaller than tolerance
 if typl == IntPatchIType::Circle || typl == IntPatchIType::Ellipse {
  let radius = match curve {
   Curve3::Circle(c) => c.radius,
   Curve3::Ellipse(e) => e.major_radius,
   _ => 0.0,
  };
  let a_tol_3d = crate::TOLERANCE_ABS;
  if radius > 0.0 && radius < a_tol_3d {
   return Vec::new();
  }
 }

 let aNbVtx = vertices.len();
 if aNbVtx == 0 {
  return Vec::new();
 }

 // Build vertex array with parameters projected to [0, 2閿? (OCCT L492-495).
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

 // OCCT L504: create last vertex at first.param + 2閿?
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
 // Add the last vertex (first.param + 2閿?, skip if duplicate of last deduped.
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
  if !(in1 && in2) {
   if std::env::var("RCAD_DBG_FF").is_ok() { eprintln!("[FF]   interval [{:.4},{:.4}] REJECTED: in1={} in2={}", t1, t2, in1, in2); }
   continue;
  }
  if std::env::var("RCAD_DBG_FF").is_ok() { eprintln!("[FF]   interval [{:.4},{:.4}] ACCEPTED", t1, t2); }
  result.push([t1, t2]);
 }

 result
}


/// BuildPCurves for all curve-surface type combinations.
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

#[cfg(test)]
mod tests {
    #[test]
    fn surface_type_index_plane() {
        let s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(glam::DVec3::Z, glam::DVec3::Z));
        assert_eq!(super::PaveFiller::surface_type_index(&s), 0);
    }
    #[test]
    fn surface_type_index_cylinder() {
        let s = rcad_kernel::geom::Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin: glam::DVec3::Z, axis: glam::DVec3::Z, ref_dir: glam::DVec3::X, radius: 1.0,
        });
        assert_eq!(super::PaveFiller::surface_type_index(&s), 1);
    }
    #[test]
    fn surface_type_index_sphere() {
        let s = rcad_kernel::geom::Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center: glam::DVec3::Z, axis: glam::DVec3::Z, ref_dir: glam::DVec3::X, radius: 1.0,
        });
        assert_eq!(super::PaveFiller::surface_type_index(&s), 3);
    }
    #[test]
    fn surface_type_index_other() {
        let s = rcad_kernel::geom::Surface3::Plane(rcad_kernel::geom::Plane::new(glam::DVec3::Z, glam::DVec3::Z));
        // BSpline variant  -- use Plane to test 'other' path; BSpline has no Default
        /* Skip: no Default for Bezier/BSpline */
    }
}
