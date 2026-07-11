use std::collections::{HashMap, HashSet};

use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::PCurve;

use crate::bvh::{Aabb, DsBvh};
use crate::bopds::ds::{
 DS, DSEdge, DSCurveRepOnFace, Interference, IntersectionCurve, ShapeOrigin,
};
use crate::bopds::pave::*;
use crate::inttools;
use crate::tolerance::*;
use super::helpers::*;

impl<'a> super::PaveFiller<'a> {
 pub(super) fn make_blocks(&mut self) {
 if std::env::var("RCAD_DEBUG_MB").is_ok() {
 eprintln!("[MB] ENTER make_blocks");
 }
 // OCCT L652-655: GlueOff guard
 if self.use_glue() {
 return;
 }

 // OCCT L657-659: Collect FF interferences (InterfFF)
 let ff_interfs: Vec<(usize, usize, Vec<usize>, Vec<usize>)> = self.ds.interf_ff.iter()
 .map(|ff| (ff.f1, ff.f2, ff.curves.clone(), ff.points.clone()))
 .collect();

 // OCCT L660-663: Early return when no FF interferences
 if ff_interfs.is_empty() {
 return;
 }

 // OCCT L666-669: Local variables
 let a_nb_ff = ff_interfs.len();
 let mut n_f1: usize;
 let mut n_f2: usize;
 let mut n_v1: usize;
 let mut n_v2: usize;
 let mut a_t1: f64;
 let mut a_t2: f64;
 let mut b_exist: bool; // OCCT L668: bExist, reused across points/curves

 // OCCT L681-683: Edge shape (skip  ?rcad uses DSEdge, not TopoDS_Edge)

 // OCCT L685-718: Per-iteration collections (simplified  ?no IncAllocator in rcad)
 // OCCT L687: aLSE  ?shared edges between the two faces
 let mut a_lse: Vec<usize> = Vec::new();
 // OCCT L689-694: Vertex maps for ON/IN/Common/Stick/EF/Bounds
 let mut a_mv_on_in: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mv_common: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mv_stick: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mv_ef: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
 // OCCT L695-696: PaveBlock maps
 let mut a_mpb_on_in: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mpb_common: std::collections::HashSet<usize> = std::collections::HashSet::new();
 // OCCT L699: aMVTol  ?vertex tolerance map (needs UnBind, use Vec<(usize, f64)>)
 let mut a_mv_tol: Vec<(usize, f64)> = Vec::new();
 // OCCT L704: aLPB  ?temporary list of PaveBlocks from update()
 let mut a_lpb: Vec<PaveBlock> = Vec::new();
 // OCCT L706: aMSCPB  — map from section edge index to (interf_index, curve_index)
 let mut a_mscpb: std::collections::HashMap<usize, (usize, usize)> = std::collections::HashMap::new();
 // OCCT L707: aMVI  — map from DS vertex index -> DS vertex index (identity)
 let mut a_mvi: std::collections::HashSet<usize> = std::collections::HashSet::new();
 // OCCT L708-709: aDMExEdges  ?map PB -> list of existing edges
 // OCCT L710: aDMNewSD  ?map old vertex -> new SD vertex
 let mut a_dm_new_sd: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 // OCCT L712: aDMVLV  ?vertex-vertex coincidence map
 let mut a_dm_vlv: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
 // OCCT L695: aMVBounds  ?bound vertices from PutBoundPaveOnCurve
 let mut a_mv_bounds: std::collections::HashSet<usize> = std::collections::HashSet::new();
 // OCCT L714: aMicroPB  ?micro PBs (too short for valid range)
 let mut a_micro_pb: Vec<PaveBlock> = Vec::new();
 // OCCT L715-716: aVertsOnRejectedPB
 let mut a_verts_on_rejected_pb: Vec<usize> = Vec::new();
 // OCCT L717: aPBFacesMap  ?map PB -> list of faces to add it to
 let mut a_pb_faces_map: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
 // OCCT L704: aMPBAdd  ?set of existing PBs already processed (cross-iteration)
 let mut a_mpb_add: std::collections::HashSet<usize> = std::collections::HashSet::new();
 // OCCT L720: aFFToRecheck  ?indices of FF pairs needing recheck
 let mut a_ff_to_recheck: Vec<usize> = Vec::new();
 let a_nb_ff_prev = a_nb_ff;

 // Cross-loop state: section edge tracking
 let mut existing_edge_map: std::collections::HashMap<(usize, usize, usize, usize), usize> = std::collections::HashMap::new();
 let mut reg_sec_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
 self.ds.section_edge_refs = vec![Vec::new(); self.ds.intersection_curves.len()];

 // OCCT L725-1107: Loop over FF pairs
 for i in 0..a_nb_ff {
 // OCCT L731-733: Recheck logic
 let cur_ind = if i < a_nb_ff_prev { i } else { a_ff_to_recheck[i - a_nb_ff_prev] };

 // OCCT L735-736: Get FF pair indices
 let (n_f1_val, n_f2_val) = (ff_interfs[cur_ind].0, ff_interfs[cur_ind].1);
 n_f1 = n_f1_val;
 n_f2 = n_f2_val;
 let curves_of_ff = &ff_interfs[cur_ind].2;
 let points_of_ff = &ff_interfs[cur_ind].3;

 // OCCT L738-745: Get points and curves of this FF pair
 let a_nb_c = curves_of_ff.len();
 let a_nb_p = points_of_ff.len();
 // OCCT L738: aVP = aFF.ChangePoints
 // OCCT L740: aVC = aFF.ChangeCurves
 if a_nb_p == 0 && a_nb_c == 0 {
 continue; // OCCT L742-744: skip if no points AND no curves
 }

 // OCCT L747-748: Get face references (skip  ?rcad uses DS indices)
 // OCCT L750: aTolFF
 let a_tol_ff = self.ff_tol(n_f1, n_f2);

 // OCCT L752-753: FaceInfo references
 // OCCT L755-767: Clear per-iteration collections
 a_mv_on_in.clear();
 a_mv_common.clear();
 a_mpb_on_in.clear();
 a_mpb_common.clear();
 a_mv_stick.clear();
 a_mv_ef.clear();
 a_mi.clear();
 a_mv_tol.clear();
 a_lse.clear();
 a_mv_bounds.clear();

 // OCCT L770: SubShapesOnIn(nF1, nF2, aMVOnIn, aMVCommon, aMPBOnIn, aMPBCommon)
 {
 let f1_face = &self.ds.faces[n_f1];
 let f2_face = &self.ds.faces[n_f2];
 // Step 1: PBs ON/IN — add PB + endpoints to aMVOnIn + aMPBOnIn
 let pb_sets_fn = |pb_set: &indexmap::IndexSet<usize>,
                   mpb_set: &mut std::collections::HashSet<usize>,
                   mv_set: &mut std::collections::HashSet<usize>| {
  for &pb_idx in pb_set {
   mpb_set.insert(pb_idx);
   if pb_idx < self.ds.pave_blocks.len() {
    let pb = &self.ds.pave_blocks[pb_idx];
    let (v1, v2) = pb.0.read().unwrap().indices();
    mv_set.insert(v1);
    mv_set.insert(v2);
   }
  }
 };
 pb_sets_fn(&f1_face.face_info.pave_blocks_on, &mut a_mpb_on_in, &mut a_mv_on_in);
 pb_sets_fn(&f1_face.face_info.pave_blocks_in, &mut a_mpb_on_in, &mut a_mv_on_in);
 pb_sets_fn(&f2_face.face_info.pave_blocks_on, &mut a_mpb_on_in, &mut a_mv_on_in);
 pb_sets_fn(&f2_face.face_info.pave_blocks_in, &mut a_mpb_on_in, &mut a_mv_on_in);
 // Step 2: Common PBs — PBsOn1 also in PBsOn2/PBsIn2
 for &pb_idx in &f1_face.face_info.pave_blocks_on {
  if f2_face.face_info.pave_blocks_on.contains(&pb_idx)
  || f2_face.face_info.pave_blocks_in.contains(&pb_idx) {
   a_mpb_common.insert(pb_idx);
   if pb_idx < self.ds.pave_blocks.len() {
    let pb = &self.ds.pave_blocks[pb_idx];
    let (v1, v2) = pb.0.read().unwrap().indices();
    a_mv_common.insert(v1);
    a_mv_common.insert(v2);
   }
  }
 }
 // Step 3: f1 VerticesOn/In also in f2 → aMVOnIn + aMVCommon
 for &vi in &f1_face.face_info.vertices_on {
  if f2_face.face_info.vertices_on.contains(&vi) || f2_face.face_info.vertices_in.contains(&vi) {
   a_mv_on_in.insert(vi);
   a_mv_common.insert(vi);
  } else {
   a_mv_on_in.insert(vi);
  }
 }
 for &vi in &f1_face.face_info.vertices_in {
  if f2_face.face_info.vertices_on.contains(&vi) || f2_face.face_info.vertices_in.contains(&vi) {
   a_mv_on_in.insert(vi);
   a_mv_common.insert(vi);
  } else {
   a_mv_on_in.insert(vi);
  }
 }
 // Step 3b: f2 VerticesOn/In — process opposite face too (OCCT missing in rcad)
 for &vi in &f2_face.face_info.vertices_on {
  if f1_face.face_info.vertices_on.contains(&vi) || f1_face.face_info.vertices_in.contains(&vi) {
   a_mv_on_in.insert(vi);
   a_mv_common.insert(vi);
  } else {
   a_mv_on_in.insert(vi);
  }
 }
 for &vi in &f2_face.face_info.vertices_in {
  if f1_face.face_info.vertices_on.contains(&vi) || f1_face.face_info.vertices_in.contains(&vi) {
   a_mv_on_in.insert(vi);
   a_mv_common.insert(vi);
  } else {
   a_mv_on_in.insert(vi);
  }
 }
 }

 // OCCT L771: SharedEdges(nF1, nF2, aLSE)
 {
 let f1_edge_set: std::collections::HashSet<usize> = self.ds.faces[n_f1].boundary_edges.iter().copied().collect();
 for &ei in &self.ds.faces[n_f2].boundary_edges {
 if f1_edge_set.contains(&ei) {
 a_lse.push(ei);
 }
 }
 }

 // OCCT L775-793: Treat Points
 for j in 0..a_nb_p {
  // OCCT L781-782: BOPDS_Point& aNP = aVP.ChangeValue(j);
  // OCCT: const gp_Pnt& aP = aNP.Pnt();
  if points_of_ff[j] >= self.ds.ff_points.len() { continue; }
  let a_p = self.ds.ff_points[points_of_ff[j]];

  // OCCT L784: IsExistingVertex(aP, aTolFF, aMVOnIn)
  b_exist = self.is_existing_vertex_at_point(a_p, a_tol_ff, &a_mv_on_in);

  if !b_exist {
   // OCCT L787: BOPTools_AlgoTools::MakeNewVertex(aP, aTolFF, aV)
   let n_v = self.ds.add_vertex(a_p);
   self.ds.vertices[n_v].geom_tol = a_tol_ff;

   // SKIP: rcad has no aMSCPB (CoupleOfPaveBlocks map) equivalent yet.
   // OCCT L789-791: aCPB.SetIndexInterf(i); aCPB.SetIndex(j); aMSCPB.Add(aV, aCPB)
  }
 }

 // OCCT L793: GetStickVertices(nF1, nF2, aMVStick, aMVEF, aMI)
 // OCCT BOPAlgo_PaveFiller_6.cxx L2879-2937
 {
 // Build aMI = all sub-shapes of nF1  ?nF2 (GetFullShapeMap)
 // OCCT L2896-2897: GetFullShapeMap(nF1, aMI); GetFullShapeMap(nF2, aMI)
 // rcad: use build_face_shape_map which collects all edges/vertices of a face
 let mut a_mi = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
 let a_mi_b = crate::pave_filler::build_face_shape_map(self.ds, n_f2);
 // OCCT L2900-2920: VV, VE, EE, VF with HasIndexNew
 // VE: stick vertices from VE
 for inf in &self.ds.interf_ve {
 let belongs = a_mi.contains(&inf.vertex) || a_mi_b.contains(&inf.vertex);
 if !belongs { continue; }
 a_mv_stick.insert(inf.vertex);
 a_mi.insert(inf.vertex);
 }
 // VF: stick vertices from VF
 for inf in &self.ds.interf_vf {
 let belongs = a_mi.contains(&inf.vertex) || a_mi_b.contains(&inf.vertex);
 if !belongs { continue; }
 a_mv_stick.insert(inf.vertex);
 a_mi.insert(inf.vertex);
 }
 // EE: collect new vertices
 for inf in &self.ds.interf_ee {
 if inf.new_vertex == usize::MAX { continue; }
 let s1_in_pair = a_mi.contains(&inf.e1) || a_mi_b.contains(&inf.e1);
 let s2_in_pair = a_mi.contains(&inf.e2) || a_mi_b.contains(&inf.e2);
 if !s1_in_pair || !s2_in_pair { continue; }
 if let Some(n_v) = self.ds.has_shape_sd(inf.new_vertex) {
 a_mv_stick.insert(n_v);
 a_mi.insert(n_v);
 } else {
 a_mv_stick.insert(inf.new_vertex);
 a_mi.insert(inf.new_vertex);
 }
 }
 // VV: collect merged vertices
 for inf in &self.ds.interf_vv {
 if inf.merged_vertex == usize::MAX { continue; }
 let s1_in_pair = a_mi.contains(&inf.v1) || a_mi_b.contains(&inf.v1);
 let s2_in_pair = a_mi.contains(&inf.v2) || a_mi_b.contains(&inf.v2);
 if !s1_in_pair || !s2_in_pair { continue; }
 if let Some(n_v) = self.ds.has_shape_sd(inf.merged_vertex) {
 a_mv_stick.insert(n_v);
 a_mi.insert(n_v);
 } else {
 a_mv_stick.insert(inf.merged_vertex);
 a_mi.insert(inf.merged_vertex);
 }
 }
 // OCCT L2921-2936: EF interferences
 for inf in &self.ds.interf_ef {
 if inf.new_vertex == usize::MAX { continue; }
 // EF vertices were already validated by the EF pass; add all directly.
 // rcad: build_face_shape_map doesn't include face index, so skip pair check.
 if let Some(n_v) = self.ds.has_shape_sd(inf.new_vertex) {
 a_mv_stick.insert(n_v);
 a_mv_ef.insert(n_v);
 a_mi.insert(n_v);
 } else {
 a_mv_stick.insert(inf.new_vertex);
 a_mv_ef.insert(inf.new_vertex);
 a_mi.insert(inf.new_vertex);
 }
 }
 }

 // OCCT L796-809: Loop over curves  ?PutPavesOnCurve
 let aMI = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
 let aMI_ref = &aMI;
 for &ci in curves_of_ff {
 if ci >= self.ds.intersection_curves.len() { continue; }
 // OCCT L799-800: aNC.InitPaveBlock1()
 // ensures at least one PaveBlock exists for PutPavesOnCurve to add ext_paves to.
 if self.ds.intersection_curves[ci].pave_blocks.is_empty() {
  let pb = crate::bopds::pave::PaveBlock::new_curve_block();
  self.ds.intersection_curves[ci].pave_blocks.push(crate::bopds::pave::SharedPB::new(pb));
 }

 // OCCT L802-808: PutPavesOnCurve(aMVOnIn, aMVCommon, aNC, aMI, aMVEF, aMVTol, aDMVLV)
 self.put_paves_on_curve(&a_mv_on_in, &a_mv_common, ci, &aMI, &a_mv_ef);
 }

 // OCCT L814: FilterPavesOnCurves  ?remove bad paves across all curves
 self.filter_paves_on_curves(curves_of_ff);

 // Second loop over curves
 for (j, &ci) in curves_of_ff.iter().enumerate() {
 if ci >= self.ds.intersection_curves.len() { continue; }
 // OCCT L821: PutStickPavesOnCurve(aF1, aF2, aMI, aVC, j, aMVStick, aMVTol, aDMVLV)
 self.put_stick_paves_on_curve(ci, &aMI, &a_mv_stick);

 // OCCT L823-826: PutEFPavesOnCurve (single curve case)
 if a_nb_c == 1 {
  self.put_ef_paves_on_curve(ci, &aMI, &a_mv_ef);
 }

 // OCCT L828-843: PutBoundPaveOnCurve(aF1, aF2, aNC, aLBV) + aDMBV
 // OCCT BOPAlgo_PaveFiller_6.cxx L2340-2400
 // For each un-vertexed endpoint of the IC, check if it lies on both face
 // surfaces (3D distance check, not UV domain), and if so create a vertex
 // and append it as an ext_pave on the IC's PB.
 {
 let mut a_lbv: Vec<usize> = Vec::new();
 // Clone IC data to avoid borrow conflict
 let ic_data = {
 let ic = &self.ds.intersection_curves[ci];
 (ic.curve.clone(), ic.t_range, ic.geom_tol, ic.start_vertex, ic.end_vertex)
 };
 let (ic_curve, ic_t_range, ic_geom_tol, ic_sv, ic_ev) = ic_data;
 let a_tol_r3d = ic_geom_tol.max(crate::tolerance::TOLERANCE_ABS);
 let a_t = [ic_t_range[0], ic_t_range[1]];
 let a_p = [ic_curve.point_at(a_t[0]), ic_curve.point_at(a_t[1])];
 // Ensure pave_block1 exists (InitPaveBlock1 equivalent).
 // OCCT BOPAlgo_PaveFiller_6.cxx L2437-2456: InitPaveBlock1 allocates a
 // PaveBlock with the curve's t_range and sv/ev endpoints, enabling
 // subsequent put_pave_on_curve to append ext_paves.  Without this,
 // analytical curves never acquire pave blocks and the section-edge
 // pipeline processes nothing.
 if self.ds.intersection_curves[ci].pave_blocks.is_empty() {
 let local_sv = self.ds.intersection_curves[ci].start_vertex;
 let local_ev = self.ds.intersection_curves[ci].end_vertex;
 if local_sv < self.ds.vertices.len() && local_ev < self.ds.vertices.len() {
 let mut pb = crate::bopds::pave::PaveBlock::new(
 crate::bopds::pave::NO_EDGE,
 crate::bopds::pave::Pave { vertex_idx: local_sv, param: ic_t_range[0] },
 crate::bopds::pave::Pave { vertex_idx: local_ev, param: ic_t_range[1] },
 );
 pb.curve = Some(ic_curve.clone());
 self.ds.intersection_curves[ci].pave_blocks.push(crate::bopds::pave::SharedPB::new(pb));
 } else if std::env::var("RCAD_DBG_MB").is_ok() {
 eprintln!("[MB] InitPaveBlock1 skip: sv={} ev={} nV={} ci={}",
 local_sv, local_ev, self.ds.vertices.len(), ci);
 }
 }
 // getBoundPaves: check which endpoints already have vertices
 // by comparing ext_pave positions and IC sv/ev against bound points.
 let mut a_bnd_nv = [usize::MAX; 2];
 {
 let ic = &self.ds.intersection_curves[ci];
 if let Some(pb) = ic.pave_blocks.first() {
 for ep in &pb.0.read().unwrap().ext_paves {
 let pt = ic_curve.point_at(ep.param);
 for j in 0..2 {
 if (pt - a_p[j]).length_squared() < a_tol_r3d * a_tol_r3d {
 a_bnd_nv[j] = ep.vertex_idx;
 }
 }
 }
 }
 if ic_sv < self.ds.vertices.len() {
 let sv_pt = self.ds.vertices[ic_sv].point;
 for j in 0..2 {
 if (sv_pt - a_p[j]).length_squared() < a_tol_r3d * a_tol_r3d {
 a_bnd_nv[j] = ic_sv;
 }
 }
 }
 if ic_ev < self.ds.vertices.len() {
 let ev_pt = self.ds.vertices[ic_ev].point;
 for j in 0..2 {
 if (ev_pt - a_p[j]).length_squared() < a_tol_r3d * a_tol_r3d {
 a_bnd_nv[j] = ic_ev;
 }
 }
 }
 }
 let a_tol_v_new = crate::tolerance::TOLERANCE_ABS;
 let is_closed = a_p[1].distance(a_p[0]) < a_tol_v_new;
 if is_closed && (a_bnd_nv[0] != usize::MAX || a_bnd_nv[1] != usize::MAX) {
 // OCCT L2357-2360: closed curve with endpoints  ?nothing to do
 }
 for j in 0..2 {
 if a_bnd_nv[j] != usize::MAX { continue; }
 if j == 1 && is_closed { continue; }
 // OCCT L2372: IsValidPointForFaces  ?3D distance to each face surface
 let mut bvf = true;
 for &fi in &[n_f1, n_f2] {
 if fi == usize::MAX { continue; }
 let surf = &self.ds.faces[fi].surface;
 let dist = match surf {
 Surface3::Plane(p) => (a_p[j] - p.origin).dot(p.normal).abs(),
 Surface3::Sphere(s) => ((a_p[j] - s.center).length() - s.radius).abs(),
 Surface3::Cylinder(c) => {
 let v = a_p[j] - c.origin;
 let axis = c.axis.normalize();
 let radial = v - axis * v.dot(axis);
 (radial.length() - c.radius).abs()
 }
 _ => f64::MAX,
 };
 if dist > a_tol_r3d { bvf = false; break; }
 }
 if !bvf { continue; }
 // OCCT L2377-2396: create vertex + add DS + append ext_pave
 let n_vn = self.ds.add_vertex(a_p[j]);
 self.ds.vertices[n_vn].geom_tol = a_tol_r3d;
 if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first_mut() {
 pb.0.write().unwrap().append_ext_pave(crate::bopds::pave::Pave {
 vertex_idx: n_vn, param: a_t[j],
 });
 }
 a_lbv.push(n_vn);
 a_mv_bounds.insert(n_vn);
 }
 }
 } // OCCT L844: end second curve loop

 // OCCT L847-851: PutClosingPaveOnCurve for each curve
 for &ci in curves_of_ff {
 if ci >= self.ds.intersection_curves.len() { continue; }
 self.put_closing_pave_on_curve(ci);
 }

 // OCCT L874-894: Build aPBTree (BOPTools_BoxTree) from ON/IN PBs
 let mut a_pb_indices: Vec<usize> = Vec::new();
 let mut a_pb_aabbs: Vec<Aabb> = Vec::new();
 for &pb_idx in &a_mpb_on_in {
  if pb_idx >= self.ds.pave_blocks.len() { continue; }
  let pb = &self.ds.pave_blocks[pb_idx];
  if pb.0.read().unwrap().new_edge.is_none() && pb.0.read().unwrap().original_edge == NO_EDGE { continue; }
  let ei = pb.0.read().unwrap().new_edge.unwrap_or(pb.0.read().unwrap().original_edge);
  if ei >= self.ds.edges.len() { continue; }
  if self.ds.is_edge_degenerated(ei) { continue; }
  // Compute edge AABB from start/end vertices + tolerance
  let sv = self.ds.edges[ei].start_vertex;
  let ev = self.ds.edges[ei].end_vertex;
  let aabb = if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
   let tol = self.ds.edges[ei].geom_tol;
   let mn = self.ds.vertices[sv].point.min(self.ds.vertices[ev].point) - DVec3::splat(tol);
   let mx = self.ds.vertices[sv].point.max(self.ds.vertices[ev].point) + DVec3::splat(tol);
   Aabb { min: mn, max: mx }
  } else { continue; };
  a_pb_indices.push(pb_idx);
  a_pb_aabbs.push(aabb);
 }
 let a_pb_tree = if !a_pb_indices.is_empty() {
  Some(DsBvh::build(a_pb_indices, a_pb_aabbs))
 } else { None };

 // OCCT L877-879: Check if this FF pair needs rechecking
 let mut is_to_recheck = a_nb_c > 0 && i < a_nb_ff_prev;

 // OCCT L882-1066: Make section edges (third loop over curves)
 for (j, &ci) in curves_of_ff.iter().enumerate() {
 if ci >= self.ds.intersection_curves.len() { continue; }

 // OCCT L884-886: Get curve data
 // aTolR3D = max(aNC.Tolerance(), aNC.TangentialTolerance())
 let a_tol_r3d = {
 let ic = &self.ds.intersection_curves[ci];
 ic.geom_tol.max(ic.curve_extra.tangential_tol)
 };

 // OCCT L888-892: aLPBC = aNC.ChangePaveBlocks(); aPB1 = aNC.ChangePaveBlock1()
 // aLPB.Clear(); aPB1->Update(aLPB, false);
 a_lpb.clear();
 {
 let ic = &mut self.ds.intersection_curves[ci];
  if let Some(mut pb1) = ic.pave_blocks.first().map(|spb| spb.0.write().unwrap()) {
 let sub_pbs = pb1.update(false);
 a_lpb = sub_pbs;
 }
 }

 // OCCT L894-897: if (aLPB.Extent() != 0) { isToRecheck = false; }
 if !a_lpb.is_empty() {
 is_to_recheck = false;
 }

 // OCCT L899-1063: Process each sub-PB
 for a_pb in &a_lpb {
 // OCCT L903-904: aPB->Indices(nV1, nV2); aPB->Range(aT1, aT2);
  (n_v1, n_v2) = a_pb.indices();
  (a_t1, a_t2) = a_pb.range();

 // OCCT L906-909: fabs(aT1-aT2) < Precision::PConfusion()  ?continue
 // Precision::PConfusion() = Confusion() * 0.01 = 1e-9
 if (a_t2 - a_t1).abs() < 1e-9 {
 continue;
 }

 // OCCT L912-918: IsValidBlockForFaces
 // OCCT L733: aMidPar = IntTools_Tools::IntermediatePoint(theT1, theT2)
 // PAR_T = 10*e^(-PI) = 0.43213918, mid = (1-PAR_T)*T1 + PAR_T*T2
 let ic = &self.ds.intersection_curves[ci];
 let curve = &ic.curve;
 let mid_t = 0.56786082 * a_t1 + 0.43213918 * a_t2;
 let mid_pt = curve.point_at(mid_t);
 let mut ok = true;
 for (k, &fi) in [n_f1, n_f2].iter().enumerate() {
 if fi == usize::MAX { continue; }
 // OCCT L746-756: if pcurve exists, D0(aMidPar, uv), IsPointInOnFace
 let pcurve = if k == 0 { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
 if let Some(pc) = pcurve {
 let uv = pc.point_at(mid_t);
 // OCCT L752: IsPointInOnFace  ?true if State IN or ON
 if !self.context.is_point_in_on_face(self.ds, fi, uv) { ok = false; break; }
 } else {
 // OCCT L759: IsValidPointForFace(aP, aF, theTol)  ?project 3D point onto surface
 let surf = if k == 0 { &self.ds.faces[n_f1].surface } else { &self.ds.faces[n_f2].surface };
 if !self.context.is_valid_point_for_face(mid_pt, fi, a_tol_r3d) { ok = false; break; }
 }
 }
 if !ok { continue; } // OCCT L755: bFlag false  ?skip this PB

 // OCCT L920-930: IsExistingPaveBlock via aLSE (shared edges)
 // OCCT BOPAlgo_PaveFiller_6.cxx L2020-2075
 // Uses geometry-based detection: intermediate point  ?bounding box  ?ComputePE
 let mut n_e_out: usize = usize::MAX;
 let mut a_tol_new: f64 = -1.0;
 let b_exist_lse = {
 if a_lse.is_empty() {
 false
 } else {
 let a_tm = 0.56786082 * a_t1 + 0.43213918 * a_t2;
 let a_pm = {
 let ic = &self.ds.intersection_curves[ci];
 ic.curve.point_at(a_tm)
 };
 let a_tol = {
 let v1_tol = if n_v1 < self.ds.vertices.len() { self.ds.vertices[n_v1].geom_tol } else { a_tol_r3d };
 let v2_tol = if n_v2 < self.ds.vertices.len() { self.ds.vertices[n_v2].geom_tol } else { a_tol_r3d };
 v1_tol.max(v2_tol)
 };
 let mut found = false;
 let mut best_dist = f64::MAX;
 for &sei in &a_lse {
 if sei >= self.ds.edges.len() { continue; }
 let se = &self.ds.edges[sei];
 let a_tol_e = se.geom_tol;
 let a_tol_check = a_tol_e.max(a_tol);
 // ComputePE: project a_pm onto edge curve, check distance
 let (_t, a_proj) = crate::extrema::closest_point_on_curve(&se.curve, a_pm);
 let dist = (a_proj - a_pm).length();
 if dist <= a_tol_check && dist < best_dist {
 found = true;
 n_e_out = sei;
 a_tol_new = dist;
 best_dist = dist;
 }
 }
 found
 }
 };
 if b_exist_lse {
 // OCCT L926-930: UpdateEdgeTolerance + UpdateSavedTolerance
 if a_tol_new > 0.0 {
 // Update edge tolerance
 if n_e_out < self.ds.edges.len() {
 self.ds.edges[n_e_out].geom_tol = self.ds.edges[n_e_out].geom_tol.max(a_tol_new);
 }
 // Save vertex tolerances
 for &vi in &[n_v1, n_v2] {
 if vi < self.ds.vertices.len() {
 a_mv_tol.push((vi, self.ds.vertices[vi].geom_tol));
 }
 }
 }
 continue;
 }

 // OCCT L936-960: FindValidRange check
 let has_valid_range = {
 if n_v1 < self.ds.vertices.len() && n_v2 < self.ds.vertices.len() {
 let v1_pt = self.ds.vertices[n_v1].point;
 let v2_pt = self.ds.vertices[n_v2].point;
 let v1_tol = a_tol_r3d.max(self.ds.vertices[n_v1].geom_tol);
 let v2_tol = a_tol_r3d.max(self.ds.vertices[n_v2].geom_tol);
 let ic = &self.ds.intersection_curves[ci];
 find_valid_range(&ic.curve, a_t1, a_t2, a_tol_r3d, v1_pt, v1_tol, v2_pt, v2_tol).is_some()
 } else { false }
 };
 if !has_valid_range {
 // OCCT L951-959: aMicroPB.Add(aPB); aMVI.Bind
 // But only if neither vertex is a bound vertex (aMVBounds guard).
 // Bound vertices are IC endpoints  ?their PBs are handled
 // separately in post-treatment, not as micro edges.
 if !a_mv_bounds.contains(&n_v1) && !a_mv_bounds.contains(&n_v2) {
 a_micro_pb.push(a_pb.clone());
 }
 continue;
 }

 // OCCT L962-1021: IsExistingPaveBlock via aMPBOnIn + aPBTree
 let a_tm = 0.56786082 * a_t1 + 0.43213918 * a_t2;
 let a_pm = curve.point_at(a_tm);
 let a_p1 = curve.point_at(a_t1);
 let a_p2 = curve.point_at(a_t2);
 let a_tol_v11 = if n_v1 < self.ds.vertices.len() { self.ds.vertices[n_v1].geom_tol } else { a_tol_r3d };
 let a_tol_v12 = if n_v2 < self.ds.vertices.len() { self.ds.vertices[n_v2].geom_tol } else { a_tol_r3d };
 let a_tol_v1 = a_tol_v11.max(a_tol_v12);
 let a_tol_check = a_tol_r3d;
 let a_max_tol_add = 0.001_f64.min(10.0 * a_tol_check);
 // Query BVH tree for candidate PBs near sub-PB mid-point
 let candidates: Vec<usize> = if let Some(ref pb_tree) = a_pb_tree {
  let query_box = Aabb {
   min: a_pm - DVec3::splat(a_tol_v1 + a_tol_check),
   max: a_pm + DVec3::splat(a_tol_v1 + a_tol_check),
  };
  pb_tree.query_aabb(&query_box)
 } else { Vec::new() };
 let b_exist_on_in = {
  if candidates.is_empty() {
  false
  } else {
 let mut found_pb_idx = usize::MAX;
 let mut best_dist = f64::MAX;
 let mut best_a_tol_new = -1.0;

 for &pb_idx in &candidates {
 if pb_idx >= self.ds.pave_blocks.len() { continue; }
 let existing_pb = &self.ds.pave_blocks[pb_idx];
 // OCCT L2154-2155: nV21, nV22
 let (n_v21, n_v22) = existing_pb.0.read().unwrap().indices();
 // OCCT L2157-2159: aTolV21, aTolV22, aTolV2
 let a_tol_v21 = if n_v21 < self.ds.vertices.len() { self.ds.vertices[n_v21].geom_tol } else { a_tol_r3d };
 let a_tol_v22 = if n_v22 < self.ds.vertices.len() { self.ds.vertices[n_v22].geom_tol } else { a_tol_r3d };
 let a_tol_v2 = a_tol_v21.max(a_tol_v22);
 // OCCT L2165: iFlag1  ?vertex match for start
 let i_flag1 = n_v1 == n_v21 || n_v1 == n_v22;
 // OCCT L2166: iFlag2  ?vertex match for end
 // OR edge AABB overlaps end-point AABB (!aBoxSp.IsOut(aBoxP2))
 let edge_ei = existing_pb.0.read().unwrap().new_edge.unwrap_or(existing_pb.0.read().unwrap().original_edge);
 let i_flag2 = if n_v2 == n_v21 || n_v2 == n_v22 {
 true
 } else if edge_ei < self.ds.edges.len() {
 // OCCT: aBoxSp (edge AABB, from ShapeInfo) vs aBoxP2
 let sv = self.ds.edges[edge_ei].start_vertex;
 let ev = self.ds.edges[edge_ei].end_vertex;
 let e_min = if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
 self.ds.vertices[sv].point.min(self.ds.vertices[ev].point)
 } else { a_p2 };
 let e_max = if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
 self.ds.vertices[sv].point.max(self.ds.vertices[ev].point)
 } else { a_p2 };
 let e_tol = a_tol_v21.max(a_tol_v22);
 let sp_min = e_min - DVec3::splat(e_tol);
 let sp_max = e_max + DVec3::splat(e_tol);
 // aBoxP2 enlarged by aTolV12
 let p2_min = a_p2 - DVec3::splat(a_tol_v12);
 let p2_max = a_p2 + DVec3::splat(a_tol_v12);
 // AABB overlap (OCCT: !aBoxSp.IsOut(aBoxP2))
 !(sp_max.x < p2_min.x || sp_min.x > p2_max.x
 || sp_max.y < p2_min.y || sp_min.y > p2_max.y
 || sp_max.z < p2_min.z || sp_min.z > p2_max.z)
 } else { false };
 if !i_flag2 { continue; }

 let edge_idx = existing_pb.0.read().unwrap().new_edge.unwrap_or(existing_pb.0.read().unwrap().original_edge);
 if edge_idx >= self.ds.edges.len() { continue; }
 let existing_edge = &self.ds.edges[edge_idx];

 // OCCT L2173-2176: init aCoeff, aDistm1m2, aPEStatus
 let mut a_coeff = 1.0;
 let mut a_dist_m1m2 = 0.0;
 let mut a_pe_status = 1;
 // OCCT L2178: aRealTol
 let mut a_real_tol = a_tol_check;
 // OCCT L2179-2187: IsCommonBlock
 if a_mpb_common.contains(&pb_idx) {
 a_real_tol = a_real_tol.max(a_tol_v1.max(a_tol_v2));
 // theMPBCommon check  ?rcad uses a_mpb_common (already true)
 a_real_tol *= 2.0;
 }
 // OCCT L2189-2230: iFlag1==2 && iFlag2==2 (both vertices match)
 // Tangent-based tolerance increase for non-linear edges
 // rcad: simplified  ?skip the tangent check for now (linear approximation)
 // OCCT L2232-2252: Mid-point projection
 let (_t, proj) = crate::extrema::closest_point_on_curve(
 &existing_edge.curve, a_pm);
 let dist_to_sp = (proj - a_pm).length();
 if dist_to_sp > a_real_tol { continue; }

 // OCCT L2254-2261: if iFlag1==1, project P1 onto edge
 let mut dist_p1 = f64::MAX;
 if !i_flag1 {
 let (_t1, p1_proj) = crate::extrema::closest_point_on_curve(
 &existing_edge.curve, a_p1);
 dist_p1 = (p1_proj - a_p1).length();
 }
 // OCCT L2263-2270: if iFlag2==1 (bbox-only, not vertex), project P2
 let mut dist_to_use = dist_to_sp;
 if n_v2 != n_v21 && n_v2 != n_v22 {
 // iFlag2 was bbox-only (not vertex match)  ?project P2
 let (_t2, p2_proj) = crate::extrema::closest_point_on_curve(
 &existing_edge.curve, a_p2);
 let dist_p2 = (p2_proj - a_p2).length();
 if dist_to_use < dist_p2 {
 dist_to_use = dist_p2;
 }
 }

 // OCCT L2272-2280: select best candidate
 let i_flag1_ok = i_flag1 || dist_p1 <= a_real_tol;
 if i_flag1_ok && dist_to_use < best_dist {
 found_pb_idx = pb_idx;
 best_a_tol_new = a_coeff * dist_to_use;
 best_dist = dist_to_use;
 }
 }
 if found_pb_idx != usize::MAX {
 n_e_out = self.ds.pave_blocks[found_pb_idx].0.read().unwrap().new_edge.unwrap_or(
 self.ds.pave_blocks[found_pb_idx].0.read().unwrap().original_edge);
 a_tol_new = best_a_tol_new;
 true
 } else {
 false
 }
 }
 };
 if b_exist_on_in {
 // OCCT L964-1021: Existing PB found, may need to add to other face
 let existing_pb = &self.ds.pave_blocks[candidates.iter().find_map(
 |&p| if self.ds.pave_blocks[p].0.read().unwrap().new_edge.unwrap_or(
 self.ds.pave_blocks[p].0.read().unwrap().original_edge) == n_e_out
 { Some(p) } else { None }
 ).unwrap_or(usize::MAX)];
 if existing_pb.0.read().unwrap().new_edge.is_some() || existing_pb.0.read().unwrap().original_edge < self.ds.edges.len() {
  // Find the matching PB index in candidates
  let pb_idx_f1 = candidates.iter().find(|&&p| {
  let e = self.ds.pave_blocks[p].0.read().unwrap().new_edge
  .unwrap_or(self.ds.pave_blocks[p].0.read().unwrap().original_edge);
  e == n_e_out
  }).copied().unwrap_or(usize::MAX);
  let b_in_f1 = {
  self.ds.faces[n_f1].face_info.pave_blocks_on.contains(&pb_idx_f1)
  || self.ds.faces[n_f1].face_info.pave_blocks_in.contains(&pb_idx_f1)
  };
  let pb_idx_f2 = candidates.iter().find(|&&p| {
  let e = self.ds.pave_blocks[p].0.read().unwrap().new_edge
  .unwrap_or(self.ds.pave_blocks[p].0.read().unwrap().original_edge);
  e == n_e_out
  }).copied().unwrap_or(usize::MAX);
  let b_in_f2 = {
  self.ds.faces[n_f2].face_info.pave_blocks_on.contains(&pb_idx_f2)
  || self.ds.faces[n_f2].face_info.pave_blocks_in.contains(&pb_idx_f2)
  };
 if !b_in_f1 || !b_in_f2 {
 // Update edge tolerance: OCCT L968-985
 if n_e_out < self.ds.edges.len() {
 self.ds.edges[n_e_out].geom_tol = self.ds.edges[n_e_out].geom_tol.max(a_tol_new);
 }
 // aPBFacesMap: OCCT L988-993
 let n_f = if b_in_f1 { n_f2 } else { n_f1 };
 a_pb_faces_map.entry(candidates.iter().find(|&&p| {
 let e = self.ds.pave_blocks[p].0.read().unwrap().new_edge.unwrap_or(self.ds.pave_blocks[p].0.read().unwrap().original_edge);
 e == n_e_out
 }).copied().unwrap_or(usize::MAX))
 .or_default()
 .push(n_f);
 // PreparePostTreatFF: OCCT L1015-1021
 // Append PB to aLPBC, register in aMSCPB/aMVI
 // rcad: register PB in both faces' pave_blocks_sc
 // OCCT L1046: aMPBAdd.Add(aPBOut)  ?only process once
 if let Some(&pb_idx) = candidates.iter().find(|&&p| {
 let e = self.ds.pave_blocks[p].0.read().unwrap().new_edge.unwrap_or(self.ds.pave_blocks[p].0.read().unwrap().original_edge);
 e == n_e_out
 }) {
 if a_mpb_add.insert(pb_idx) {
 let ic_curves = &mut self.ds.intersection_curves[ci];
 ic_curves.pave_blocks.push(self.ds.pave_blocks[pb_idx].clone());
 for &fi in &[n_f1, n_f2] {
 if fi != usize::MAX {
 self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
 }
 }
 }
 }
 }
 }
 continue;
 }

 // OCCT L1023-1044: MakeEdge + MakePCurve
 let ic = &self.ds.intersection_curves[ci];
 let pca = ic.pcurve_on_a.clone();
 let pcb = ic.pcurve_on_b.clone();
 let new_ei = crate::boptools::make_edge(self.ds, ci, n_v1, n_v2, a_t1, a_t2, a_tol_r3d);
 crate::boptools::make_pcurve(
  self.ds, new_ei, n_f1, n_f2, ci,
  self.section_attribute.pcurve_on_s1,
  self.section_attribute.pcurve_on_s2,
  pca.as_ref(), pcb.as_ref(),
  Some([a_t1, a_t2]), Some([a_t1, a_t2]),
 );
 // set PB edge and register in section_edge_refs
 if new_ei < self.ds.edges.len() {
  if let Some(epb) = self.ds.edges[new_ei].pave_blocks.first_mut() {
   epb.0.write().unwrap().new_edge = Some(new_ei);
  }
 }
 self.ds.section_edge_refs[ci].push(new_ei);
 // OCCT L1066-1067: aLPBC.Append(aPB)
 let mut sub_pb = a_pb.clone();
 sub_pb.new_edge = Some(new_ei);
 // OCCT L1069-1075: aMSCPB.Add(aES, aCPB) + aMVI.Bind(aV1, nV1)
 a_mscpb.insert(new_ei, (cur_ind, j));
 a_mvi.insert(n_v1);
 a_mvi.insert(n_v2);
 // rcad: allocate a global PB and register on both faces' pave_blocks_sc.
 let g_pb_idx = self.ds.allocate_pave_block(sub_pb.clone());
 for &fi in &[n_f1, n_f2] {
 if fi != usize::MAX {
 self.ds.faces[fi].face_info.pave_blocks_sc.insert(g_pb_idx);
 }
 }
 // OCCT L1079-1080: aMVTol.UnBind(nV1/nV2)
 a_mv_tol.retain(|&(v, _)| v != n_v1 && v != n_v2);
 // OCCT L1082-1094: ProcessExistingPaveBlocks
 for &pb_idx in &candidates {
 if pb_idx >= self.ds.pave_blocks.len() { continue; }
 // OCCT L3139: theMPB.Contains(aPBF)  ?skip already-processed
 if a_mpb_add.contains(&pb_idx) { continue; }
 let a_pbf = &self.ds.pave_blocks[pb_idx];
 let (pbsv, pbev) = a_pbf.0.read().unwrap().indices();
 // Check if PB shares vertices with the new edge
 if pbsv == n_v1 || pbsv == n_v2 || pbev == n_v1 || pbev == n_v2 {
 a_mpb_add.insert(pb_idx);
 // Check if this PB is already in both faces' ON/IN
 let b_in_f1 = self.ds.faces[n_f1].face_info.pave_blocks_on.contains(&pb_idx)
 || self.ds.faces[n_f1].face_info.pave_blocks_in.contains(&pb_idx);
 let b_in_f2 = self.ds.faces[n_f2].face_info.pave_blocks_on.contains(&pb_idx)
 || self.ds.faces[n_f2].face_info.pave_blocks_in.contains(&pb_idx);
 if b_in_f1 && b_in_f2 {
 // Register in curve PB list + pave_blocks_sc
 self.ds.intersection_curves[ci].pave_blocks.push(a_pbf.clone());
 for &fi in &[n_f1, n_f2] {
 if fi != usize::MAX {
 self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
 }
 }
 } else {
 // Add to PBFacesMap for the missing face
 let n_f = if b_in_f1 { n_f2 } else { n_f1 };
 a_pb_faces_map.entry(pb_idx).or_default().push(n_f);
 // Register in curve PB list + pave_blocks_sc for both faces
 self.ds.intersection_curves[ci].pave_blocks.push(a_pbf.clone());
 for &fi in &[n_f1, n_f2] {
 if fi != usize::MAX {
 self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
 }
 }
 }
 }
 }
 } // OCCT L1063: end sub-PB loop

 // OCCT L1065: aLPBC.RemoveFirst()  ?remove the parent PB
 if ci < self.ds.intersection_curves.len() {
 let ic = &mut self.ds.intersection_curves[ci];
 if !ic.pave_blocks.is_empty() {
 ic.pave_blocks.remove(0);
 }
 }
 } // OCCT L1066: end Make section edges loop

 // OCCT L1067-1071: Recheck logic
 if is_to_recheck {
 a_ff_to_recheck.push(cur_ind);
 }

 // OCCT L1073-1095: Restore vertex tolerances + reset bounding boxes
 a_mv_tol.sort_by(|a, b| a.0.cmp(&b.0));
 a_mv_tol.dedup_by_key(|a| a.0);
 for &(n_v, saved_tol) in &a_mv_tol {
 if n_v < self.ds.vertices.len() {
 // OCCT L1112-1116: Restore ORIGINAL tolerance (not max)
 self.ds.vertices[n_v].geom_tol = saved_tol;
 }
 }
 // OCCT L1091-1094: UnBind from aDMVLV (separate loop)
 for &(n_v, _) in &a_mv_tol {
 a_dm_vlv.remove(&n_v);
 }

 // OCCT L1097-1106: ProcessExistingPaveBlocks (post-section-edge)
 // Registers existing PBs from ON/IN sets that overlap with new section edges
 // into the section curve for each face.
 {
 // Collect PB indices to add to each face's pave_blocks_sc
 let mut pbs_to_add: Vec<(usize, usize)> = Vec::new(); // (fi, pb_idx)
 for &ci in curves_of_ff {
 if ci >= self.ds.intersection_curves.len() { continue; }
 for &sei in &self.ds.section_edge_refs[ci] {
 let se = &self.ds.edges[sei];
 let (sv, ev) = (se.start_vertex, se.end_vertex);
 for &fi in &[n_f1, n_f2] {
 let face = &self.ds.faces[fi];
 for &pb_idx in face.face_info.pave_blocks_on.iter()
 .chain(face.face_info.pave_blocks_in.iter())
 {
 if pb_idx < self.ds.pave_blocks.len() {
 let pb = &self.ds.pave_blocks[pb_idx];
 let (pbsv, pbev) = pb.0.read().unwrap().indices();
 if pbsv == sv || pbsv == ev || pbev == sv || pbev == ev {
 pbs_to_add.push((fi, pb_idx));
 }
 }
 }
 }
 }
 }
 // Apply collected adds
 for &(fi, pb_idx) in &pbs_to_add {
 self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
 }
 }
 } // OCCT L1107: end FF pair loop

 // ===== Post-loop phases (OCCT L1109-1136) =====

 // OCCT L1109-1110: RemoveMicroSectionEdges
 // Micro section edges are PBs whose FindValidRange failed (too short).
 // The a_micro_pb list has been populated during section edge creation.
 // OCCT BOPAlgo_PaveFiller_6.cxx L4341-4419
 // rcad: match section edges to micro PBs by (vertex_idx) pair, since
 // the sub-PBs from update() don't have new_edge set (unlike OCCT where
 // aMSCPB maps edge shape -> PB directly).
 {
 for ci in 0..self.ds.intersection_curves.len() {
 let mut keep: Vec<usize> = Vec::new();
 for &sei in &self.ds.section_edge_refs[ci] {
 let is_micro = if sei < self.ds.edges.len() {
 let e = &self.ds.edges[sei];
 a_micro_pb.iter().any(|pb: &PaveBlock| {
  let (pv1, pv2) = pb.indices();
 (pv1 == e.start_vertex && pv2 == e.end_vertex)
 || (pv1 == e.end_vertex && pv2 == e.start_vertex)
 })
 } else { false };
 if !is_micro {
 keep.push(sei);
 } else if sei < self.ds.edges.len() {
 self.ds.edges[sei].pave_blocks.clear();
 }
 }
 self.ds.section_edge_refs[ci] = keep;
 }
 }

 // OCCT L1112: MakeSDVerticesFF(aDMVLV, aDMNewSD)
 // Create SD vertices for coinciding VV/VE/VF vertex groups
 // OCCT BOPAlgo_PaveFiller_6.cxx L1173-1193
 {
 // a_dm_vlv: map from vertex index to list of coincident vertices
 let key_list: Vec<(usize, Vec<usize>)> = a_dm_vlv.iter()
 .map(|(k, v)| (*k, v.clone()))
 .collect();
 for (_, verts) in &key_list {
 if verts.len() < 2 { continue; }
 // Create one SD vertex for the group.  Use the centroid of all
 // vertices in the group as the SD point.
 let mut sum = glam::DVec3::ZERO;
 for &v in verts {
 if v < self.ds.vertices.len() {
 sum += self.ds.vertices[v].point;
 }
 }
 let sd_pt = sum / verts.len() as f64;
 // Find the "best" existing vertex in the group (the one closest to centroid)
 let mut best_v = verts[0];
 let mut best_dist = f64::MAX;
 for &v in verts {
 if v < self.ds.vertices.len() {
 let d = (self.ds.vertices[v].point - sd_pt).length_squared();
 if d < best_dist { best_v = v; best_dist = d; }
 }
 }
 // Map all other vertices to the best one
 for &v in verts {
 if v != best_v && v < self.ds.vertices.len() {
 a_dm_new_sd.insert(v, best_v);
 }
 }
 }
 }

 // OCCT L1114-1120: PostTreatFF(aMSCPB, aDMExEdges, aDMNewSD, aMicroPB,
 // aVertsOnRejectedPB, aAllocator, theRange)
 // Post-process section edges: create missing PBs, register in face info.
 self.post_treat_ff();

 // OCCT L1125-1126: CorrectToleranceOfSE()
 // Reduce tolerance of section edges where appropriate.
 for ci in 0..self.ds.intersection_curves.len() {
 for &sei in &self.ds.section_edge_refs[ci] {
 if sei < self.ds.edges.len() {
 let edge_tol = self.ds.edges[sei].geom_tol;
 let curve_tol = if ci < self.ds.intersection_curves.len() {
 self.ds.intersection_curves[ci].geom_tol
 } else { edge_tol };
 // Use the smaller of edge and curve tolerance (OCCT CorrectToleranceOfSE)
 self.ds.edges[sei].geom_tol = edge_tol.min(curve_tol).max(TOLERANCE_ABS);
 }
 }
 }

 // OCCT L1127-1128: UpdateFaceInfo(aDMExEdges, aDMNewSD, aPBFacesMap)
 // Register section edge PBs in their missing faces via aPBFacesMap.
 // Then recompute vertices_in for each face from curve endpoints.
 for (&pb_idx, faces) in &a_pb_faces_map {
 if pb_idx < self.ds.pave_blocks.len() {
 for &fi in faces {
 if fi < self.ds.faces.len() {
 self.ds.faces[fi].face_info.pave_blocks_sc.insert(pb_idx);
 }
 }
 }
 }
 for fi in 0..self.ds.faces.len() {
 for &ci in self.ds.faces[fi].face_info.curves_sc_only().iter() {
 if ci < self.ds.intersection_curves.len() {
 let ic = &self.ds.intersection_curves[ci];
 self.ds.faces[fi].face_info.vertices_in.insert(ic.start_vertex);
 self.ds.faces[fi].face_info.vertices_in.insert(ic.end_vertex);
 }
 }
 }

 // OCCT L1129-1130: UpdatePaveBlocks(aDMNewSD)
 // Update PB vertex indices for SD vertices
 for (old_v, new_v) in &a_dm_new_sd {
 for ei in 0..self.ds.edges.len() {
  for spb in &mut self.ds.edges[ei].pave_blocks {
  let mut pb = spb.0.write().unwrap();
  if pb.pave1.vertex_idx == *old_v { pb.pave1.vertex_idx = *new_v; }
  if pb.pave2.vertex_idx == *old_v { pb.pave2.vertex_idx = *new_v; }
 }
 }
 for fi in 0..self.ds.faces.len() {
 if self.ds.faces[fi].face_info.vertices_in.contains(old_v) {
 self.ds.faces[fi].face_info.vertices_in.remove(old_v);
 self.ds.faces[fi].face_info.vertices_in.insert(*new_v);
 }
 }
 }

 // OCCT L1133-1136: PutSEInOtherFaces
 self.put_se_in_other_faces();

 // OCCT-aligned: Build edge images
 self.ds.build_edge_images();

 if std::env::var("RCAD_DEBUG_SPLIT").is_ok() {
 let n_circle = self.ds.intersection_curves.iter().filter(|ic| matches!(ic.curve, Curve3::Circle(_))).count();
 let n_total = self.ds.intersection_curves.len();
 eprintln!("[SPLIT] END_MAKE_BLOCKS total_curves={} circle_curves={}", n_total, n_circle);
 for fi in 0..self.ds.faces.len() {
 let face = &self.ds.faces[fi];
 eprintln!("[SPLIT] face[{}] curves_sc={} vertices_in={}", fi, face.face_info.curves_sc.len(), face.face_info.vertices_in.len());
 }
 }

 // OCCT-aligned: InitPaveBlock1 is called per-curve in the first loop above.
 // The empty loop here is removed -- its purpose was served.
 }
}

/// OCCT-aligned: PostTreatFF (BOPAlgo_PaveFiller_6.cxx L1197-?).
/// Post-processes section edges created by MakeBlocks.
/// Creates missing PBs, registers in face info, and handles
/// technological vertices.  Currently a stub -- the full ~2000-line
/// OCCT function body will be translated in a subsequent alignment pass.
impl<'a> super::PaveFiller<'a> {
 pub(super) fn post_treat_ff(&mut self) {
  // OCCT L1208-1212: aNbS = theMSCPB.Extent(); if (!aNbS) return;
  // rcad: section edges are registered in section_edge_refs.
  let total_section_edges: usize = self.ds.section_edge_refs.iter().map(|v| v.len()).sum();
  if total_section_edges == 0 {
   return;
  }
  // OCCT PostTreatFF body (PaveFiller_6.cxx L1197-?).
  // Handles: registering section edges in face info, creating
  // SD vertices for micro edges, cleaning technological vertices,
  // and fixing up pave blocks for same-domain face pairs.
  //
  // rcad: the core registration is already done in the third loop
  // (Make section edges) above via section_edge_refs and
  // pave_blocks_sc.  Additional passes (technological vertices,
  // multi-face section edge propagation) are deferred.
 }
}


