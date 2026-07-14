use std::collections::{HashMap, HashSet};

use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::PCurve;

use crate::bopalgo;
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
 if self.use_glue() {
 return;
 }
 let ff_interfs: Vec<(usize, usize, Vec<usize>, Vec<usize>)> = self.ds.interf_ff.iter()
 .map(|ff| (ff.f1, ff.f2, ff.curves.clone(), ff.points.clone()))
 .collect();
 if ff_interfs.is_empty() {
 return;
 }
 let a_nb_ff = ff_interfs.len();
 let mut n_f1: usize;
 let mut n_f2: usize;
 let mut n_v1: usize;
 let mut n_v2: usize;
 let mut a_t1: f64;
 let mut a_t2: f64;
 let mut b_exist: bool;
 let mut a_lse: Vec<usize> = Vec::new();
 let mut a_mv_on_in: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mv_common: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mv_stick: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mv_ef: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mi: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mpb_on_in: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mpb_common: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_mv_tol: Vec<(usize, f64)> = Vec::new();
 let mut a_lpb: Vec<PaveBlock> = Vec::new();
 let mut a_mscpb: std::collections::HashMap<usize, (usize, usize)> = std::collections::HashMap::new();
 let mut a_mvi: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_dm_new_sd: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
 let mut a_dm_vlv: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
 let mut a_mv_bounds: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_micro_pb: Vec<PaveBlock> = Vec::new();
 let mut a_verts_on_rejected_pb: Vec<usize> = Vec::new();
 let mut a_pb_faces_map: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
 let mut a_mpb_add: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_ff_to_recheck: Vec<usize> = Vec::new();
 let a_nb_ff_prev = a_nb_ff;

 // Cross-loop state: section edge tracking
 let mut existing_edge_map: std::collections::HashMap<(usize, usize, usize, usize), usize> = std::collections::HashMap::new();
 let mut reg_sec_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
 self.ds.section_edge_refs = vec![Vec::new(); self.ds.intersection_curves.len()];
 for i in 0..a_nb_ff {
 let cur_ind = if i < a_nb_ff_prev { i } else { a_ff_to_recheck[i - a_nb_ff_prev] };
 let (n_f1_val, n_f2_val) = (ff_interfs[cur_ind].0, ff_interfs[cur_ind].1);
 n_f1 = n_f1_val;
 n_f2 = n_f2_val;
 let curves_of_ff = &ff_interfs[cur_ind].2;
 let points_of_ff = &ff_interfs[cur_ind].3;
 let a_nb_c = curves_of_ff.len();
 let a_nb_p = points_of_ff.len();
 if a_nb_p == 0 && a_nb_c == 0 {
 continue;
 }
 let a_tol_ff = self.ff_tol(n_f1, n_f2);
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
 {
 let f1_face = &self.ds.faces[n_f1];
 let f2_face = &self.ds.faces[n_f2];
 // Step 1: PBs ON/IN  ?add PB + endpoints to aMVOnIn + aMPBOnIn
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
 // Step 2: Common PBs  ?PBsOn1 also in PBsOn2/PBsIn2
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
 // OCCT-aligned: SubShapesOnIn  ?Face1 vertices in Face2  ?aMVOnIn + aMVCommon
 // OCCT BOPDS_DS.cxx L1124-1142: only Face1's VOn/VIn checked against Face2's VOn/VIn
 for &vi in &f1_face.face_info.vertices_on {
  if f2_face.face_info.vertices_on.contains(&vi) || f2_face.face_info.vertices_in.contains(&vi) {
   a_mv_on_in.insert(vi);
   a_mv_common.insert(vi);
  }
 }
 for &vi in &f1_face.face_info.vertices_in {
  if f2_face.face_info.vertices_on.contains(&vi) || f2_face.face_info.vertices_in.contains(&vi) {
   a_mv_on_in.insert(vi);
   a_mv_common.insert(vi);
  }
 }
 }
 // OCCT-aligned: SharedEdges (BOPDS_DS.cxx L1147-1208)
 {
 let mut a_first_face_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for &ei in &self.ds.faces[n_f1].boundary_edges {
  let pbs = self.ds.pave_blocks(ei);
  if pbs.is_empty() {
   a_first_face_edges.insert(ei);
  } else {
   for pb in pbs {
    let re = self.ds.real_pave_block_edge(ei, pb)
     .or(pb.0.read().unwrap().new_edge)
     .unwrap_or(ei);
    a_first_face_edges.insert(re);
   }
  }
 }
 for &ei in &self.ds.faces[n_f2].boundary_edges {
  let pbs = self.ds.pave_blocks(ei);
  if pbs.is_empty() {
   if a_first_face_edges.contains(&ei) {
    a_lse.push(ei);
   }
  } else {
   for pb in pbs {
    let re = self.ds.real_pave_block_edge(ei, pb)
     .or(pb.0.read().unwrap().new_edge)
     .unwrap_or(ei);
    if a_first_face_edges.contains(&re) {
     a_lse.push(re);
    }
   }
  }
 }
 }
 // OCCT L775-793: 1. Treat Points
 for j in 0..a_nb_p {
  // OCCT: const gp_Pnt& aP = aNP.Pnt();
  if points_of_ff[j] >= self.ds.ff_points.len() { continue; }
  let a_p = self.ds.ff_points[points_of_ff[j]];
  b_exist = self.is_existing_vertex_at_point(a_p, a_tol_ff, &a_mv_on_in);
  if !b_exist {
   let n_v = self.ds.add_vertex(a_p);
   self.ds.vertices[n_v].geom_tol = a_tol_ff;
   // OCCT: aCPB.SetIndexInterf(i) + aCPB.SetIndex(j) + aMSCPB.Add(aV, aCPB)
   a_mscpb.insert(n_v, (cur_ind, j));
  }
 }
 // OCCT BOPAlgo_PaveFiller_6.cxx L2879-2937
 {
 // Build aMI = all sub-shapes of nF1  ?nF2 (GetFullShapeMap)
 // rcad: use build_face_shape_map which collects all edges/vertices of a face
 let mut a_mi = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
 let a_mi_b = crate::pave_filler::build_face_shape_map(self.ds, n_f2);
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
 for inf in &self.ds.interf_ef {
 if inf.new_vertex == usize::MAX { continue; }
 let e_in_pair = a_mi.contains(&inf.edge) || a_mi_b.contains(&inf.edge);
 let f_in_pair = a_mi.contains(&inf.face) || a_mi_b.contains(&inf.face);
 if !e_in_pair || !f_in_pair { continue; }
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
 let aMI = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
 let aMI_ref = &aMI;
 for &ci in curves_of_ff {
 if ci >= self.ds.intersection_curves.len() { continue; }
 // ensures at least one PaveBlock exists for PutPavesOnCurve to add ext_paves to.
 if self.ds.intersection_curves[ci].pave_blocks.is_empty() {
  let pb = crate::bopds::pave::PaveBlock::new_curve_block();
  self.ds.intersection_curves[ci].pave_blocks.push(crate::bopds::pave::SharedPB::new(pb));
 }
 self.put_paves_on_curve(&a_mv_on_in, &a_mv_common, ci, &aMI, &a_mv_ef);
 }
 self.filter_paves_on_curves(curves_of_ff);

 // Second loop over curves
 for (j, &ci) in curves_of_ff.iter().enumerate() {
 if ci >= self.ds.intersection_curves.len() { continue; }
 self.put_stick_paves_on_curve(ci, &aMI, &a_mv_stick);
 if a_nb_c == 1 {
  self.put_ef_paves_on_curve(ci, &aMI, &a_mv_ef);
 }
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
 }
 for j in 0..2 {
 if a_bnd_nv[j] != usize::MAX { continue; }
 if j == 1 && is_closed { continue; }
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
 // OCCT-aligned: UV boundary transition detection.  Sample the curve on each
 // face's FClass2d to find where the pcurve crosses the trimmed face boundary.
 // OCCT achieves this via PutBoundPaveOnCurve + aIC.Bounds() (clipped curve).
 {
 let a_faces = [n_f1, n_f2];
 for (k, &fi) in a_faces.iter().enumerate() {
  if fi == usize::MAX { continue; }
  let pc_opt = if k == 0 { self.ds.intersection_curves[ci].pcurve_on_a.clone() }
               else { self.ds.intersection_curves[ci].pcurve_on_b.clone() };
  let Some(pc) = pc_opt else { continue };
  let tt0 = ic_t_range[0];
  let tt1 = ic_t_range[1];
  let span = tt1 - tt0;
  if span <= TOLERANCE_CLAMP_MIN { continue; }
  let n_samp = 129;
  let first_in = self.context.is_point_in_on_face(self.ds, fi, pc.point_at(tt0));
  if std::env::var("RCAD_DEBUG_MB").is_ok() {
   eprintln!("[MB]  UV_bd ci={} k={} fi={} first_in={} n_samp={}", ci, k, fi, first_in, n_samp);
  }
  let mut prev_t = tt0;
  let mut prev_in = first_in;
  for i in 1..=n_samp {
   let t = tt0 + span * i as f64 / n_samp as f64;
   let in_on = self.context.is_point_in_on_face(self.ds, fi, pc.point_at(t));
   if in_on != prev_in {
    let mut lo = prev_t;
    let mut hi = t;
    for _ in 0..20 {
     let mid = (lo + hi) * 0.5;
     let mid_in = self.context.is_point_in_on_face(self.ds, fi, pc.point_at(mid));
     if mid_in == prev_in { lo = mid; } else { hi = mid; }
    }
    let ct = (lo + hi) * 0.5;
    let cp = ic_curve.point_at(ct);
    let nv = self.ds.add_vertex(cp);
    self.ds.vertices[nv].geom_tol = a_tol_r3d;
    if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first_mut() {
     pb.0.write().unwrap().append_ext_pave(Pave { vertex_idx: nv, param: ct });
    }
    a_lbv.push(nv);
    prev_in = in_on;
   }
   prev_t = t;
  }
 }
 }
 }
 }
 for &ci in curves_of_ff {
 if ci >= self.ds.intersection_curves.len() { continue; }
 self.put_closing_pave_on_curve(ci);
 }
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
 let mut is_to_recheck = a_nb_c > 0 && i < a_nb_ff_prev;
 for (j, &ci) in curves_of_ff.iter().enumerate() {
 if ci >= self.ds.intersection_curves.len() { continue; }
 // aTolR3D = max(aNC.Tolerance(), aNC.TangentialTolerance())
 let a_tol_r3d = {
 let ic = &self.ds.intersection_curves[ci];
 ic.geom_tol.max(ic.curve_extra.tangential_tol)
 };
 // aLPB.Clear(); aPB1->Update(aLPB, false);
 a_lpb.clear();
 {
 let ic = &mut self.ds.intersection_curves[ci];
  if let Some(mut pb1) = ic.pave_blocks.first().map(|spb| spb.0.write().unwrap()) {
  // OCCT-aligned: ensure pave1/pave2 use the curve's clipped endpoint vertices.
  // clip_curve_to_face_uv sets t_range and start/end_vertex; propagate to PB.
  if pb1.original_edge == NO_EDGE && ic.start_vertex < self.ds.vertices.len()
   && ic.end_vertex < self.ds.vertices.len() {
   pb1.pave1 = crate::bopds::pave::Pave { vertex_idx: ic.start_vertex, param: ic.t_range[0] };
   pb1.pave2 = crate::bopds::pave::Pave { vertex_idx: ic.end_vertex, param: ic.t_range[1] };
  }
  // OCCT-aligned: Update with theFlag=false (BOPAlgo_PaveFiller_6.cxx L912).
  // The flag=false means ext_paves alone define the sub-PB boundaries.
  // The InitPaveBlock1 + PutPavesOnCurve handling of IC endpoints ensures
  // the sv/ev are present as ext_paves, so Update(false) still produces
  // valid sub-PBs covering the full curve range.
 let sub_pbs = pb1.update(false);
 a_lpb = sub_pbs;
 }
 }
 if !a_lpb.is_empty() {
 is_to_recheck = false;
 }
 for a_pb in &a_lpb {
  (n_v1, n_v2) = a_pb.indices();
  (a_t1, a_t2) = a_pb.range();
 // Precision::PConfusion() = Confusion() * 0.01 = 1e-9
 if (a_t2 - a_t1).abs() < 1e-9 {
  if std::env::var("RCAD_DEBUG_MB").is_ok() { eprintln!("[MB]  skip_zero_range t1={} t2={}", a_t1, a_t2); }
 continue;
 }
 // PAR_T = 10*e^(-PI) = 0.43213918, mid = (1-PAR_T)*T1 + PAR_T*T2
 let ic = &self.ds.intersection_curves[ci];
 let curve = &ic.curve;
 let mid_t = 0.56786082 * a_t1 + 0.43213918 * a_t2;
 let mid_pt = curve.point_at(mid_t);
 let mut ok = true;
 for (k, &fi) in [n_f1, n_f2].iter().enumerate() {
 if fi == usize::MAX { continue; }
 let pcurve = if k == 0 { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
 if let Some(pc) = pcurve {
  let uv = pc.point_at(mid_t);
  let in_on = self.context.is_point_in_on_face(self.ds, fi, uv);
  if !in_on {
   // OCCT-aligned: fallback to 3D projection when pcurve-based check fails.
   // This handles sub-PBs whose midpoint lies outside the trimmed face domain
   // but whose 3D point is valid on the face surface. OCCT achieves this with
   // finer sub-PB subdivision (more ext_paves); rcad falls back to projection.
   let valid = self.context.is_valid_point_for_face(mid_pt, fi, a_tol_r3d);
   if !valid {
    if std::env::var("RCAD_DEBUG_MB").is_ok() { eprintln!("[MB]   fail_face k={} fi={} pcurve_uv=({:.4},{:.4}) proj_fail", k, fi, uv.x, uv.y); }
    ok = false; break;
   }
  }
 } else {
  let surf = if k == 0 { &self.ds.faces[n_f1].surface } else { &self.ds.faces[n_f2].surface };
  let valid = self.context.is_valid_point_for_face(mid_pt, fi, a_tol_r3d);
  if !valid {
   if std::env::var("RCAD_DEBUG_MB").is_ok() { eprintln!("[MB]   fail_face k={} fi={} no_pcurve proj_fail", k, fi); }
   ok = false; break;
  }
 }
 }
 if !ok {
  if std::env::var("RCAD_DEBUG_MB").is_ok() { eprintln!("[MB]  skip_valid_block nv1={} nv2={} t1={} t2={}", n_v1, n_v2, a_t1, a_t2); }
  continue; }
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
 // But only if neither vertex is a bound vertex (aMVBounds guard).
 // Bound vertices are IC endpoints  ?their PBs are handled
 // separately in post-treatment, not as micro edges.
 if !a_mv_bounds.contains(&n_v1) && !a_mv_bounds.contains(&n_v2) {
 a_micro_pb.push(a_pb.clone());
 }
 continue;
 }
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
 let (n_v21, n_v22) = existing_pb.0.read().unwrap().indices();
 let a_tol_v21 = if n_v21 < self.ds.vertices.len() { self.ds.vertices[n_v21].geom_tol } else { a_tol_r3d };
 let a_tol_v22 = if n_v22 < self.ds.vertices.len() { self.ds.vertices[n_v22].geom_tol } else { a_tol_r3d };
 let a_tol_v2 = a_tol_v21.max(a_tol_v22);
 let i_flag1 = n_v1 == n_v21 || n_v1 == n_v22;
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
 let mut a_coeff = 1.0;
 let mut a_dist_m1m2 = 0.0;
 let mut a_pe_status = 1;
 let mut a_real_tol = a_tol_check;
 if a_mpb_common.contains(&pb_idx) {
 a_real_tol = a_real_tol.max(a_tol_v1.max(a_tol_v2));
 // theMPBCommon check  ?rcad uses a_mpb_common (already true)
 a_real_tol *= 2.0;
 }
 // Tangent-based tolerance increase for non-linear edges
 // rcad: simplified  ?skip the tangent check for now (linear approximation)
 let (_t, proj) = crate::extrema::closest_point_on_curve(
 &existing_edge.curve, a_pm);
 let dist_to_sp = (proj - a_pm).length();
 if dist_to_sp > a_real_tol { continue; }
 let mut dist_p1 = f64::MAX;
 if !i_flag1 {
 let (_t1, p1_proj) = crate::extrema::closest_point_on_curve(
 &existing_edge.curve, a_p1);
 dist_p1 = (p1_proj - a_p1).length();
 }
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
 let mut sub_pb = a_pb.clone();
 sub_pb.new_edge = Some(new_ei);
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
 a_mv_tol.retain(|&(v, _)| v != n_v1 && v != n_v2);
 for &pb_idx in &candidates {
 if pb_idx >= self.ds.pave_blocks.len() { continue; }
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
 }
 // OCCT-aligned: aLPBC.RemoveFirst() (BOPAlgo_PaveFiller_6.cxx L1097).
 // Remove the InitPaveBlock1 PB from the curve's PaveBlocks list now that
 // all sub-PBs have been processed.
 if !self.ds.intersection_curves[ci].pave_blocks.is_empty() {
  self.ds.intersection_curves[ci].pave_blocks.remove(0);
 }
 }
 if is_to_recheck {
 a_ff_to_recheck.push(cur_ind);
 }
 a_mv_tol.sort_by(|a, b| a.0.cmp(&b.0));
 a_mv_tol.dedup_by_key(|a| a.0);
 for &(n_v, saved_tol) in &a_mv_tol {
 if n_v < self.ds.vertices.len() {
 self.ds.vertices[n_v].geom_tol = saved_tol;
 }
 }
 for &(n_v, _) in &a_mv_tol {
 a_dm_vlv.remove(&n_v);
 }
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
 }

 // ===== Post-loop phase =====
 // OCCT-aligned: RemoveMicroSectionEdges (BOPAlgo_PaveFiller_6.cxx L1142, L4341-4417).
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
 // OCCT-aligned: MakeSDVerticesFF (BOPAlgo_PaveFiller_6.cxx L1145, L1173-1193).
 // Create SD vertices for coinciding VV/VE/VF vertex groups.
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
 // OCCT-aligned: PostTreatFF (BOPAlgo_PaveFiller_6.cxx L1146-1152, L1197-1701).
 // rcad uses crate::bopalgo::intersect_vertices + self.make_sd_vertices_vv
 // to fuse section-edge vertices, replacing OCCT's nested PaveFiller.
 {
 // Phase 1 (OCCT L1235-1263): Find unused vertices from FF pairs.
 // For each FF, get stick vertices, remove those used by curves' paves.
 let mut a_verts_unused: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut a_ind_map: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for ff in &self.ds.interf_ff {
  let n_f1 = ff.f1;
  let n_f2 = ff.f2;
  if n_f1 >= self.ds.faces.len() || n_f2 >= self.ds.faces.len() { continue; }
  // Build stick vertices (OCCT: GetStickVertices)
  let mut a_mv: std::collections::HashSet<usize> = std::collections::HashSet::new();
  let mut a_mi = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
  let a_mi_b = crate::pave_filler::build_face_shape_map(self.ds, n_f2);
  for v in &a_mi_b { a_mi.insert(*v); }
  // Collect vertex indices from all interference types belonging to this FF's faces
  for inf in &self.ds.interf_ve {
   if !a_mi.contains(&inf.vertex) { continue; }
   a_mv.insert(inf.vertex);
  }
  for inf in &self.ds.interf_vf {
   if !a_mi.contains(&inf.vertex) { continue; }
   a_mv.insert(inf.vertex);
  }
  for inf in &self.ds.interf_ee {
   if inf.new_vertex == usize::MAX { continue; }
   if !a_mi.contains(&inf.e1) || !a_mi.contains(&inf.e2) { continue; }
   let n_v = self.ds.has_shape_sd(inf.new_vertex).unwrap_or(inf.new_vertex);
   a_mv.insert(n_v);
  }
  for inf in &self.ds.interf_vv {
   if inf.merged_vertex == usize::MAX { continue; }
   if !a_mi.contains(&inf.v1) || !a_mi.contains(&inf.v2) { continue; }
   let n_v = self.ds.has_shape_sd(inf.merged_vertex).unwrap_or(inf.merged_vertex);
   a_mv.insert(n_v);
  }
  for inf in &self.ds.interf_ef {
   if inf.new_vertex == usize::MAX { continue; }
   let e_in = a_mi.contains(&inf.edge);
   let f_in = a_mi.contains(&inf.face);
   if !e_in || !f_in { continue; }
   let n_v = self.ds.has_shape_sd(inf.new_vertex).unwrap_or(inf.new_vertex);
   a_mv.insert(n_v);
  }
  // Remove used vertices (OCCT: RemoveUsedVertices — iterates curves' PB paves)
  for &ci in &ff.curves {
   if ci >= self.ds.intersection_curves.len() { continue; }
   let ic = &self.ds.intersection_curves[ci];
   for spb in &ic.pave_blocks {
    let pb = spb.0.read().unwrap();
    a_mv.remove(&pb.pave1.vertex_idx);
    a_mv.remove(&pb.pave2.vertex_idx);
    for ep in &pb.ext_paves {
     a_mv.remove(&ep.vertex_idx);
    }
   }
  }
  // OCCT: IndMap fence — vertices appearing once go to VertsUnused;
  // appearing twice get removed.
  for &vi in &a_mv {
   if a_ind_map.contains(&vi) {
    a_verts_unused.remove(&vi);
   } else {
    a_ind_map.insert(vi);
    a_verts_unused.insert(vi);
   }
  }
 }

 // Phase 2 (OCCT L1269-1308): Early return for single-entry case (skip, no TopoDS shapes).
 let a_nb_s = a_mscpb.len();
 let a_nb_me = a_micro_pb.len();
 if a_nb_s > 0 && a_nb_me == 0 {
  // Collect all unique vertex indices from section edges
  let mut a_post_verts: Vec<usize> = a_mvi.iter().copied().collect();
  // Add vertices from micro PBs
  for pb in &a_micro_pb {
   let (pv1, pv2) = pb.indices();
   a_post_verts.push(pv1);
   a_post_verts.push(pv2);
  }
  // Add unused vertices (Phase 1)
  for &vi in &a_verts_unused {
   a_post_verts.push(vi);
  }
  // Phase 3 (OCCT L1310-1348)+Phase 6 (OCCT L1421-1430): Fuse vertices.
  // Use OCCT-aligned intersect_vertices to group close vertices,
  // then make_sd_vertices_vv to create SD entries for each group.
  let fuzzy = self.fuzzy_tolerance;
  let blocks = crate::bopalgo::intersect_vertices(&a_post_verts, self.ds, fuzzy);
  for block in &blocks {
   if block.len() >= 2 {
    self.make_sd_vertices_vv(block);
   }
  }
  // Phase 8 (OCCT L1690-1701): Extract SD entries from shape_sd into aDMNewSD.
  for &vi in &a_post_verts {
   if let Some(sd) = self.ds.has_shape_sd(vi) {
    if vi != sd {
     a_dm_new_sd.entry(vi).or_insert(sd);
    }
   }
  }
  // Follow SD chains: resolve old→new through any intermediate mappings.
  let sd_keys: Vec<usize> = a_dm_new_sd.keys().copied().collect();
  for &k in &sd_keys {
   if let Some(&v) = a_dm_new_sd.get(&k) {
    let mut chain = v;
    let mut depth = 0;
    while let Some(&next) = a_dm_new_sd.get(&chain) {
     if next == chain || depth > 100 { break; }
     chain = next;
     depth += 1;
    }
    if chain != v {
     a_dm_new_sd.insert(k, chain);
    }
   }
  }
 }
 }

 // ⏳ OCCT-aligned: CorrectToleranceOfSE (BOPAlgo_PaveFiller_6.cxx L1158, L4105-4306).
 // rcad: only min(curve_tol, edge_tol); OCCT also reduces vertex tolerances.
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
 // ⏳ OCCT-aligned: UpdateFaceInfo (BOPAlgo_PaveFiller_6.cxx L1161, L1705-1978).
 // rcad: simple pave_blocks_sc + vertices_in; OCCT also handles CommonBlocks,
 // existing-edge PB replacement via aDMExEdges, and SD vertex remapping in
 // PaveBlocksOn/PaveBlocksIn/PaveBlocksSc.
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
 // ⏳ OCCT-aligned: UpdatePaveBlocks (BOPAlgo_PaveFiller_6.cxx L1163, L3712-3750).
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
 // OCCT-aligned: PutSEInOtherFaces (BOPAlgo_PaveFiller_6.cxx L1167-1168).
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

 }
}


