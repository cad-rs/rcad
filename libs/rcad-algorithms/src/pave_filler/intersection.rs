use std::collections::{HashSet, BTreeMap};

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve3, Surface3, any_perpendicular};
use rcad_kernel::CurveEval;

use crate::bopalgo::{fill_map, make_blocks};
use crate::bopds::ds::{DS, ShapeOrigin, InterferenceVV, InterferenceVE, InterferenceVF, InterferenceEE, InterferenceEF};
use crate::bopds::pave::Pave;
use crate::inttools;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::pave_filler::helpers::*;
use crate::tolerance::*;

impl<'a> super::PaveFiller<'a> {
 /// OCCT PaveFiller L145: BOPDS_Iterator  ?BVH pair enumeration
 pub(crate) fn build_ds_bvh(&self, is_a: bool, is_edge: bool) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let (ds_start, end) = if is_edge {
 if is_a { (0, self.ds.a_edge_count) }
 else { (self.ds.a_edge_count, self.ds.edges.len()) }
 } else {
 if is_a { (0, self.ds.a_vertex_count) }
 else { (self.ds.a_vertex_count, self.ds.vertices.len()) }
 };
 let n = end - ds_start;
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);

 for local_i in 0..n {
 let ds_i = ds_start + local_i;
 indices.push(ds_i);
 let aabb = if is_edge {
 let e = &self.ds.edges[ds_i];
 let pts = [self.ds.vertices[e.start_vertex].point,
 self.ds.vertices[e.end_vertex].point];
 let mut a = Aabb::empty();
 for &p in &pts { a.expand_point(p); }
 // Expand for edge tolerance
 let tol = e.geom_tol.max(CONFUSION);
 a.min -= DVec3::splat(tol);
 a.max += DVec3::splat(tol);
 a
 } else {
 let pt = self.ds.vertices[ds_i].point;
 let tol = self.ds.vertices[ds_i].geom_tol.max(CONFUSION);
 Aabb { min: pt - DVec3::splat(tol), max: pt + DVec3::splat(tol) }
 };
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 ///  OCCT-aligned: BOPDS_Iterator  ?build a single BVH for all elements
 /// of the given shape type (both operands A and B combined), used for
 /// single-pass cross-operand pair traversal.
 pub(crate) fn build_ds_bvh_combined(&self, is_edge: bool) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let n = if is_edge { self.ds.edges.len() } else { self.ds.vertices.len() };
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);
 for ds_i in 0..n {
 indices.push(ds_i);
 let aabb = if is_edge {
 let e = &self.ds.edges[ds_i];
 let pts = [self.ds.vertices[e.start_vertex].point,
 self.ds.vertices[e.end_vertex].point];
 let mut a = Aabb::empty();
 for &p in &pts { a.expand_point(p); }
 let tol = e.geom_tol.max(CONFUSION);
 a.min -= DVec3::splat(tol);
 a.max += DVec3::splat(tol);
 a
 } else {
 let pt = self.ds.vertices[ds_i].point;
 let tol = self.ds.vertices[ds_i].geom_tol.max(CONFUSION);
 Aabb { min: pt - DVec3::splat(tol), max: pt + DVec3::splat(tol) }
 };
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 /// OCCT PaveFiller_2.cxx L141-206: PerformVE
 /// ✅ OCCT-aligned: BOPDS_Iterator::Initialize(VERTEX, EDGE) — single pass.
 /// Cross-operand filtering is done by BOPDS_Iterator; this function
 /// applies remaining type-specific filters and calls intersect_ve.
 pub(crate) fn perform_ve_bvh(&mut self, pairs: &[(usize, usize)]) {
 use rayon::prelude::*;
 self.fill_shrunk_data();
 // pairs come pre-computed from BOPDS_Iterator
 let ds = &self.ds;
 let a_vc = ds.a_vertex_count;
 let a_ec = ds.a_edge_count;
 if pairs.is_empty() { return; }
 let filtered: Vec<(usize, usize)> = pairs.par_iter()
 .filter(|&(vi, ei)| {
 if (*vi < a_vc) == (*ei < a_ec) { return false; }
 !ds.edge_has_vertex(*vi, *ei) && !ds.edge_has_flag(*ei)
 && !ds.has_interf_ve(*vi, *ei) && !ds.has_interf_ve_via_faces(*vi, *ei)
 && !ds.is_edge_degenerated(*ei)
 && !ds.edges[*ei].pave_blocks.is_empty()
 && ds.edges[*ei].pave_blocks[0].0.read().unwrap().is_splittable
 })
 .copied()
 .collect();
 self.intersect_ve(&filtered);
 }

 ///  OCCT-aligned: IntersectVE (PaveFiller_2.cxx L212-394).
 /// Processes vertex-edge pairs with SD vertex resolution, PB endpoint
 /// dedup (aMVPB), and aDMVSD fence map, matching OCCT's structure.
 fn intersect_ve(&mut self, pairs: &[(usize, usize)]) {
 if pairs.is_empty() { return; }
 // rcad: interferences Vec (no pre-allocation needed)
 // Group vertices by edge, then resolve SD and dedup per (nVSD, nE)
 let mut edge_verts: std::collections::HashMap<usize, Vec<usize>> =
 std::collections::HashMap::new();
 for &(vi, ei) in pairs {
 edge_verts.entry(ei).or_default().push(vi);
 }

 // aDMVSD map: (nVSD, nE)  ?list of original vertices deduped to same SD root
 let mut a_dmv_sd: std::collections::HashMap<(usize, usize), Vec<usize>> =
 std::collections::HashMap::new();
 let mut a_m_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for (&ei, verts) in &edge_verts {
 let a_mv_pb: std::collections::HashSet<usize> = self.ds.edges[ei].paves.iter()
 .map(|p| p.vertex_idx)
 .collect();

 for &vi in verts {
 let n_vsd = self.ds.has_shape_sd(vi).unwrap_or(vi);
 if a_mv_pb.contains(&n_vsd) { continue; }
 let key = (n_vsd, ei);
 a_dmv_sd.entry(key).or_default().push(vi);
 }
 }
 for (&(n_vsd, ei), original_verts) in &a_dmv_sd {
 let point = self.ds.vertices[n_vsd].point;
 let edge = &self.ds.edges[ei];
 let te = self.ve_tol(n_vsd, ei);

 let t_opt = crate::pave_filler::helpers::project_vertex_to_curve(
 point, &edge.curve, te);
 let t = match t_opt {
 Some(t) if t >= edge.t_range[0] && t <= edge.t_range[1] => t,
 _ => continue,
 };
 let dist_3d = edge.curve.point_at(t).distance(point);
 if dist_3d > self.ds.vertices[n_vsd].geom_tol {
 self.ds.vertices[n_vsd].geom_tol = dist_3d;
 self.ds.increased_ss.insert(n_vsd);
 }
 // OCCT adds pave via aPave.SetIndex(nVx) using the UpdateVertex result.
 // rcad: push Pave directly to edge's pave list.
 let edge_had_paves = !self.ds.edges[ei].paves.is_empty();
 for &vi in original_verts {
 let has_vertex_at_t = self.ds.edges[ei].paves.iter()
 .any(|p| (p.param - t).abs() < TOLERANCE_ABS && p.vertex_idx == vi);
 if !has_vertex_at_t {
 self.ds.edges[ei].paves.push(Pave { vertex_idx: vi, param: t });
 }
 }
 if !edge_had_paves || self.ds.edges[ei].paves.len() > 1 {
 a_m_edges.insert(ei);
 }
 for &vi in original_verts {
 let already_interfered = self.ds.interf_ve.iter().any(|inf| {
 inf.vertex == vi && inf.edge == ei
 });
 if !already_interfered {
 self.ds.interf_ve.push(InterferenceVE{
 vertex: vi,
 edge: ei,
 param: t,
 });
 }
 }
 }
 if !a_m_edges.is_empty() {
 self.split_pave_blocks(&a_m_edges, true);
 }
 }
 /// OCCT PaveFiller_3.cxx L145-244: PerformEE
 /// ✅ OCCT-aligned: BOPDS_Iterator::Initialize(EDGE, EDGE) — single pass.
 /// Cross-operand filtering via a_edge_count.
 /// OCCT-aligned: PerformEE (PaveFiller_3.cxx L145-590).
 /// Returns the set of edges that were modified (got new paves).
 /// The caller should call split_pave_blocks for remaining edges
 /// after treat_new_vertices has processed the new vertex edges.
 pub(crate) fn perform_ee_bvh(&mut self, pairs: &[(usize, usize)])
 -> std::collections::HashSet<usize>
 {
 use rayon::prelude::*;
 self.fill_shrunk_data();
 // pairs come pre-computed from BOPDS_Iterator
 let ds = &self.ds;
 let a_ec = ds.a_edge_count;
 let blocks: Vec<(usize, usize, [f64; 2], [f64; 2])> = pairs.par_iter()
 .filter(|&(ae, be)| {
 if (*ae < a_ec) == (*be < a_ec) { return false; }
 !ds.edge_has_flag(*ae) && !ds.edge_has_flag(*be)
 && !ds.has_interf_ee(*ae, *be)
 && !ds.is_edge_degenerated(*ae) && !ds.is_edge_degenerated(*be)
 })
 .flat_map(|&(ae, be)| {
 let ra = Self::get_pb_boxes(ds, ae, ds.edges[ae].t_range);
 let rb = Self::get_pb_boxes(ds, be, ds.edges[be].t_range);
 let mut v = Vec::new();
 for &r1 in &ra { for &r2 in &rb { v.push((ae, be, r1, r2)); } }
 v
 })
 .collect();
 let mut modified: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for &(ae, be, r1, r2) in &blocks {
 self.intersect_ee(ae, be, r1, r2, &mut modified);
 }
 // SplitPaveBlocks) is handled in the caller (perform() in mod.rs):
 // - treat_new_vertices()  ?PerformNewVertices
 // - split_pave_blocks() for remaining modified edges
 modified
 }
 /// OCCT: PaveBlock range extraction (GetPBBox equivalent)
 pub(crate) fn get_pb_boxes(ds: &DS, edge_idx: usize, edge_t_range: [f64; 2]) -> Vec<[f64; 2]> {
 let paves = &ds.edges[edge_idx].paves;
 if paves.is_empty() { return vec![edge_t_range]; }
 let mut params: Vec<f64> = paves.iter().map(|p| p.param).filter(|p| p.is_finite()).collect();
 params.sort_by(|a, b| a.partial_cmp(b).unwrap());
 params.dedup();
 let tol = ds.edges[edge_idx].geom_tol.max(crate::tolerance::TOLERANCE_ABS);
 let mut ranges = Vec::new();
 let mut prev = edge_t_range[0];
 for &p in &params {
 if (p - prev).abs() > tol { ranges.push([prev, p]); }
 prev = p;
 }
 if (edge_t_range[1] - prev).abs() > tol { ranges.push([prev, edge_t_range[1]]); }
 ranges
 }
 /// ✅ OCCT-aligned: PerformVF (PaveFiller_4.cxx L165-298).
 /// BOPDS_Iterator::Initialize(VERTEX, FACE) — single BVH pass.
 /// SD vertex resolution + aMVFPairs dedup matching OCCT form.
 pub(crate) fn perform_vf_bvh(&mut self, pairs: &[(usize, usize)]) {
 use rayon::prelude::*;
 self.fill_shrunk_data();
 // pairs come pre-computed from BOPDS_Iterator
 let ds = &self.ds;
 let a_vc = ds.a_vertex_count;
 let a_fc = ds.a_face_count;
 if pairs.is_empty() { return; }
 // Skip already-interfered pairs; resolve SD vertices.
 let filtered: Vec<(usize, usize)> = pairs.par_iter()
 .filter(|&(vi, fi)| {
 if (*vi < a_vc) == (*fi < a_fc) { return false; }
 if ds.has_interf_vf(*vi, *fi) { return false; }
 if ds.has_interf_ve_via_faces(*vi, *fi) { return false; }
 true
 })
 .copied()
 .collect();
 let mut a_mvf_pairs: std::collections::HashSet<(usize, usize)> =
 std::collections::HashSet::new();
 for &(vi, fi) in &filtered {
 let n_vsd = ds.has_shape_sd(vi).unwrap_or(vi);
 a_mvf_pairs.insert((n_vsd, fi));
 }
 for &(vi, fi) in &a_mvf_pairs {
 self.check_vertex_face(vi, fi);
 }
 }
 /// OCCT PaveFiller_5.cxx L165-300: PerformEF
 pub(crate) fn perform_ef_bvh(&mut self, pairs: &[(usize, usize)]) {
 use rayon::prelude::*;
 self.fill_shrunk_data();
 // pairs come pre-computed from BOPDS_Iterator
 let ds = &self.ds;
 let a_edge_count = ds.a_edge_count;
 let a_face_count = ds.a_face_count;
 let blocks: Vec<(usize, usize, [f64; 2])> = pairs.par_iter()
 .filter(|&(ei, fi)| {
 let same_range = (*ei < a_edge_count && *fi < a_face_count)
 || (*ei >= a_edge_count && *fi >= a_face_count);
 if same_range { return false; }
 !ds.edge_has_flag(*ei) && !ds.is_edge_degenerated(*ei)
 && !ds.has_interf_ef(*ei, *fi)
 })
 .flat_map(|&(ei, fi)| {
 //  OCCT-aligned: iterate edge's PaveBlocks (L246-248: ChangePaveBlocks(nE)).
 // Skip PBs already in face's PaveBlocksOn (L257-260: aMPBF.Contains(aPBR)).
 let face_pbon: Vec<usize> = ds.faces[fi].face_info.pave_blocks_on.iter().copied().collect();
 let mut results = Vec::new();
 for pb_idx in 0..ds.edges[ei].pave_blocks.len() {
 let pb = &ds.edges[ei].pave_blocks[pb_idx];
 let real_original = pb.0.read().unwrap().original_edge;
 if face_pbon.contains(&real_original) { continue; }
 let aT1 = pb.0.read().unwrap().pave1.param;
 let aT2 = pb.0.read().unwrap().pave2.param;
 let range = [aT1.min(aT2), aT1.max(aT2)];
 results.push((ei, fi, range));
 }
 if results.is_empty() {
 // Fallback: no PBs  ?use full edge range (OCCT L240: aLPB outer iteration)
 let r = Self::get_pb_boxes(ds, ei, ds.edges[ei].t_range);
 r.into_iter().map(move |range| (ei, fi, range)).collect::<Vec<_>>()
 } else { results }
 })
 .collect();
 for &(ei, fi, r) in &blocks {
 self.intersect_ef(ei, fi, &r);
 }
 }
 /// OCCT BOPDS_Iterator: face BVH construction
 pub(crate) fn build_ds_bvh_face(&self, is_a: bool) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let (start, end) = if is_a {
 (0, self.ds.a_face_count)
 } else {
 (self.ds.a_face_count, self.ds.faces.len())
 };
 let n = end - start;
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);
 for local_i in 0..n {
 let fi = start + local_i;
 indices.push(fi);
 let f = &self.ds.faces[fi];
 let mut aabb = Aabb::empty();
 // Boundary vertices
 for &vi in &f.boundary_verts {
 if vi < self.ds.vertices.len() {
 aabb.expand_point(self.ds.vertices[vi].point);
 }
 }
 // OCCT BndLib_AddSurface: expand AABB for curved surfaces.
 // Sphere: full sphere AABB = center  ?radius (face boundary
 // vertices only cover a patch, not the whole sphere volume).
 // Cylinder/Cone: boundary vertices already span the full
 // parametric extent  ?no extra expansion needed.
 if let Surface3::Sphere(s) = &f.surface {
 let r = s.radius.abs();
 aabb.expand_point(s.center + DVec3::splat(r));
 aabb.expand_point(s.center - DVec3::splat(r));
 }
 let tol = f.geom_tol.max(CONFUSION);
 aabb.min -= DVec3::splat(tol);
 aabb.max += DVec3::splat(tol);
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 ///  OCCT-aligned: BOPDS_Iterator  ?combined face BVH (both operands).
 pub(crate) fn build_ds_bvh_face_all(&self) -> crate::bvh::DsBvh {
 use crate::bvh::{Aabb, DsBvh};
 let n = self.ds.faces.len();
 let mut indices = Vec::with_capacity(n);
 let mut aabbs = Vec::with_capacity(n);
 for fi in 0..n {
 indices.push(fi);
 let f = &self.ds.faces[fi];
 let mut aabb = Aabb::empty();
 for &vi in &f.boundary_verts {
 if vi < self.ds.vertices.len() {
 aabb.expand_point(self.ds.vertices[vi].point);
 }
 }
 if let Surface3::Sphere(s) = &f.surface {
 let r = s.radius.abs();
 aabb.expand_point(s.center + DVec3::splat(r));
 aabb.expand_point(s.center - DVec3::splat(r));
 }
 let tol = f.geom_tol.max(CONFUSION);
 aabb.min -= DVec3::splat(tol);
 aabb.max += DVec3::splat(tol);
 aabbs.push(aabb);
 }
 DsBvh::build(indices, aabbs)
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ve_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all vertices are shared, skip V-E pass
 let shared_verts = &self.ds.shared_topology.shared_vertices;
 if shared_verts.is_empty() {
 return false;
 }

 // Check if all vertices from shape A have matches in shape B
 let a_verts: std::collections::HashSet<usize> = self.ds.vertices
 .iter()
 .enumerate()
 .filter(|(_, v)| v.origin == Some(ShapeOrigin::ShapeA))
 .map(|(i, _)| i)
 .collect();

 let matched_a: std::collections::HashSet<usize> = shared_verts
 .iter()
 .map(|(a, _)| *a)
 .collect();

 a_verts == matched_a && !a_verts.is_empty()
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ee_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 let shared_edges = &self.ds.shared_topology.shared_edges;
 if shared_edges.is_empty() {
 return false;
 }

 // Check if all edges from shape A have matches in shape B
 let a_edges: std::collections::HashSet<usize> = self.ds.edges
 .iter()
 .enumerate()
 .filter(|(_, e)| e.origin == ShapeOrigin::ShapeA)
 .map(|(i, _)| i)
 .collect();

 let matched_a: std::collections::HashSet<usize> = shared_edges
 .iter()
 .map(|(a, _)| *a)
 .collect();

 a_edges == matched_a && !a_edges.is_empty()
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_vf_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all faces are fully glued, skip V-F pass
 !self.ds.shared_topology.fully_glued_faces.is_empty()
 && self.ds.shared_topology.fully_glued_faces.len()
 == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ef_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all faces are fully glued, skip E-F pass
 !self.ds.shared_topology.fully_glued_faces.is_empty()
 && self.ds.shared_topology.fully_glued_faces.len()
 == self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count)
 }
 /// rcad glue-mode acceleration (no OCCT equivalent)
 pub(crate) fn should_skip_ff_pass(&self) -> bool {
 if !self.use_glue() {
 return false;
 }

 // If all faces are fully glued, skip F-F pass
 let total_face_pairs = self.ds.a_face_count * (self.ds.faces.len() - self.ds.a_face_count);
 self.ds.shared_topology.fully_glued_faces.len() == total_face_pairs && total_face_pairs > 0
 }
 /// ✅ OCCT-aligned: PerformVV (PaveFiller_1.cxx L45-132).
 /// Builds vertex-vertex connection map (FillMap), groups connected
 /// vertices (MakeBlocks), then creates SD vertices for each group.
 /// Pairs come pre-computed from BOPDS_Iterator (cross-operand, AABB-filtered).
  pub(crate) fn perform_vv(&mut self, pairs: &[(usize, usize)]) {
  // OCCT L47-56: n1, n2, iFlag, aSize; iterator init + early return
  let a_vc = self.ds.a_vertex_count;
  let a_size = a_vc * (self.ds.vertices.len() - a_vc);
  if a_size == 0 { return; }

  // ✅ OCCT-aligned: BOPDS_Iterator(VERTEX, VERTEX) — BVH-based pair enumeration.
  // OCCT L68-76: myIterator->Initialize(VERTEX, VERTEX) returns overlapping AABB pairs.
  let mut a_mili: BTreeMap<usize, Vec<usize>> = BTreeMap::new();

  // OCCT L68-98: 1. Map V/LV — build connection map of close vertex pairs.
  for &(n1, n2) in pairs {
  // Skip same-operand pairs (OCCT: cross-operand only)
  if (n1 < a_vc) == (n2 < a_vc) { continue; }

  // OCCT L77-81: if HasInterf — FillMap + continue
  if self.ds.has_interf_vv(n1, n2) {
  fill_map(&mut a_mili, n1, n2);
  continue;
  }

  // OCCT L84-91: Resolve SD vertices (HasShapeSD) + ComputeVV
  let n1sd = self.ds.has_shape_sd(n1).unwrap_or(n1);
  let n2sd = self.ds.has_shape_sd(n2).unwrap_or(n2);

  // OCCT L93: ComputeVV(aV1, aV2, myFuzzyValue) — tolerance-based distance check
  let tol = self.vv_pair_tol(n1, n2);
  let dist = (self.ds.vertices[n1sd].point - self.ds.vertices[n2sd].point).length();
  let i_flag = if dist <= tol { 0 } else { 1 };

  // OCCT L94-97: if !iFlag (vertices interfere) — FillMap
  if i_flag == 0 {
  fill_map(&mut a_mili, n1, n2);
  }
  }
 let a_m_blocks = make_blocks(&a_mili);
 for block in &a_m_blocks {
 if block.len() < 2 { continue; }
 self.make_sd_vertices_vv(block);
 }
  // OCCT ShapesSD is a DataMap<source, target> (one direction).
  // rcad stores (source, target) bidirectionally; dedup via HashSet.
  let a_dmii: std::collections::HashSet<usize> =
    self.ds.shape_sd.sd_vertices_iter().map(|&(k, _)| k).collect();
  for &n1 in &a_dmii {
    self.ds.init_pave_blocks_for_vertex(n1);
  }
  }

 ///  OCCT-aligned: MakeSDVertices (PaveFiller_1.cxx L136-233).
 /// Merges a connected group of vertices into a single SD vertex.
 /// The first vertex in the block becomes the merge target; all other
 /// vertices in the block are remapped to it via AddShapeSD and
 /// Interference::VertexVertex entries.
 pub(super) fn make_sd_vertices_vv(&mut self, block: &[usize]) {
 if block.len() < 2 { return; }
 // nSD tracks the existing SD vertex index; others go into aLV.
 let mut n_sd: Option<usize> = None;
 let mut a_lv: Vec<usize> = Vec::with_capacity(block.len());
 for &n_x in block {
 if let Some(n_sd1) = self.ds.has_shape_sd(n_x) {
 if n_sd.is_none() {
 n_sd = Some(n_sd1);
 }
 }
 a_lv.push(n_x);
 }
 // rcad: use the minimum-index vertex as the merged target.
 let n_v = n_sd.unwrap_or_else(|| *a_lv.iter().min().unwrap());
 // rcad: geom_tol of the merged vertex is max of all members' tolerances.
 if let Some(&target) = block.iter().max_by(|&&a, &&b| {
 self.ds.vertices[a].geom_tol.partial_cmp(&self.ds.vertices[b].geom_tol).unwrap()
 }) {
 if n_v < self.ds.vertices.len() {
 self.ds.vertices[n_v].geom_tol = self.ds.vertices[n_v].geom_tol
 .max(self.ds.vertices[target].geom_tol);
 }
 }
 for i in 0..block.len() {
 let n1 = block[i];
 self.ds.add_shape_sd(n1, n_v);
 // rcad: ShapeOrigin check.
 // (OCCT L208-218: self-interfering shape warning  ?skipped for brevity)
 for j in (i + 1)..block.len() {
 let n2 = block[j];
 // rcad: push VertexVertex interference
 self.ds.interf_vv.push(InterferenceVV{
 v1: n1,
 v2: n2,
 merged_vertex: n_v,
 });
 }
 }
 }
 /// OCCT PaveFiller_2.cxx L141-206: PerformVE
 pub(crate) fn perform_ve(&mut self) {
 // with HasSubShape / HasFlag / HasInterf / HasInterfShapeSubShapes skips.
 self.fill_shrunk_data();
 let a_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeA);
 let b_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeB);
 // shrink data is computed on-the-fly in check_vertex_edge via ve_tol().
 //
 // rcad: manual O(n ? loop (see PairIterator in perform_ee for BVH pattern).

 for &vi in &a_verts {
 for &ei in &b_edges {
 if self.ds.edge_has_vertex(vi, ei) { continue; }
 if self.ds.edge_has_flag(ei) { continue; }
 if self.ds.has_interf_ve(vi, ei) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, ei) { continue; }
 if self.ds.is_edge_degenerated(ei) { continue; }
 if self.ds.edges[ei].pave_blocks.is_empty() { continue; }
 if !self.ds.edges[ei].pave_blocks[0].0.read().unwrap().is_splittable { continue; }
 self.check_vertex_edge(vi, ei);
 }
 }

 let b_verts: Vec<usize> = self.verts_of(ShapeOrigin::ShapeB);
 let a_edges: Vec<usize> = self.edges_of(ShapeOrigin::ShapeA);

 for &vi in &b_verts {
 for &ei in &a_edges {
 if self.ds.edge_has_vertex(vi, ei) { continue; }
 if self.ds.edge_has_flag(ei) { continue; }
 if self.ds.has_interf_ve(vi, ei) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, ei) { continue; }
 if self.ds.is_edge_degenerated(ei) { continue; }
 if self.ds.edges[ei].pave_blocks.is_empty() { continue; }
 if !self.ds.edges[ei].pave_blocks[0].0.read().unwrap().is_splittable { continue; }
 self.check_vertex_edge(vi, ei);
 }
 }
 }
 /// OCCT PaveFiller_2.cxx L104-121: ComputeVE
 pub(crate) fn check_vertex_edge(&mut self, vi: usize, ei: usize) {
 let point = self.ds.vertices[vi].point;
 let edge_curve = self.ds.edges[ei].curve.clone();
 let t_range = self.ds.edges[ei].t_range;
 let te = self.ve_tol(vi, ei);
 match &edge_curve {
 Curve3::Line(line) => {
 if let Some(t) = inttools::vertex_ops::vertex_on_line_with_tol(
 point,
 line,
 t_range,
 te,
 ) {
 self.ds.interf_ve.push(InterferenceVE{
 vertex: vi,
 edge: ei,
 param: t,
 });
 self.ds.edges[ei].paves.push(Pave {
 vertex_idx: vi,
 param: t,
 });
 }
 }
 Curve3::Circle(circle) => {
 // Check if point lies on the circle arc
 let v = point - circle.center;
 let dist = v.length();
 if (dist - circle.radius).abs() < te {
 let on_plane = v.dot(circle.normal).abs() < te;
 if on_plane {
 // Compute angular parameter
 let u = if circle.normal.x.abs() < 0.9 {
 circle.normal.cross(DVec3::X).normalize()
 } else {
 circle.normal.cross(DVec3::Y).normalize()
 };
 let w = circle.normal.cross(u);
 let theta = w.dot(v).atan2(u.dot(v));
 if theta >= t_range[0] - te && theta <= t_range[1] + te {
 //  OCCT-aligned: only create VE interference if the vertex is
 // within tolerance of the edge's 3D curve at the computed param.
 let on_edge_3d = edge_curve.point_at(theta).distance(point) <= te;
 if on_edge_3d {
 self.ds.interf_ve.push(InterferenceVE{
 vertex: vi,
 edge: ei,
 param: theta,
 });
 self.ds.edges[ei].paves.push(Pave {
 vertex_idx: vi,
 param: theta,
 });
 }
 }
 }
 }
 }
 _ => {
 //  OCCT-aligned: general curve projection (IntTools_Context:
 // GeomAPI_ProjectPointOnCurve for arbitrary curve types).
 // rcad: coarse 21-sample grid to find closest approach.
 let mut best_t = t_range[0];
 let mut best_d = f64::MAX;
 for si in 0..21 {
 let t = t_range[0] + (t_range[1] - t_range[0]) * (si as f64 / 20.0);
 let d = edge_curve.point_at(t).distance(point);
 if d < best_d { best_d = d; best_t = t; }
 }
 if best_d <= te {
 self.ds.interf_ve.push(InterferenceVE{
 vertex: vi,
 edge: ei,
 param: best_t,
 });
 self.ds.edges[ei].paves.push(Pave {
 vertex_idx: vi,
 param: best_t,
 });
 }
 }
 }
 }
 /// OCCT PaveFiller_3.cxx L145-244: PerformEE
 pub(crate) fn perform_ee(&mut self) {
 // with HasFlag / PaveBlock emptiness / GetPBBox skip conditions.
 self.fill_shrunk_data();
 let a_count = self.ds.a_edge_count;

 // Build a set of shared edge pairs for fast lookup when glue is enabled
 let shared_edge_set: std::collections::HashSet<(usize, usize)> = if self.use_glue() {
 self.ds
 .shared_topology
 .shared_edges
 .iter()
 .map(|(e1, e2)| (*e1, *e2))
 .collect()
 } else {
 std::collections::HashSet::new()
 };
 // rcad: PairIterator for cross-group pairs (A-edges  ?B-edges).
 // For PaveBlock-level precision (OCCT L200-232), iterate sub-ranges
 // of each edge defined by existing paves (from VE or prior intersections).
 // Each sub-range = one logical PaveBlock.
 let mut it = crate::bopds::ds::PairIterator::prepare_ab(a_count, self.ds.edges.len());
 while it.more() {
 let pk = it.value();
 let ae = pk.i1; let be = pk.i2;
 if self.ds.edge_has_flag(ae) || self.ds.edge_has_flag(be) {
 it.next(); continue;
 }
 if self.ds.has_interf_ee(ae, be) {
 it.next(); continue;
 }

 if self.ds.is_edge_degenerated(ae) || self.ds.is_edge_degenerated(be) {
 it.next(); continue;
 }
 // PaveBlocks of each edge).  rcad: build sub-ranges from existing
 // paves to limit intersection to relevant sub-segments.
 let ranges_a = self.collect_paveblock_ranges(ae, self.ds.edges[ae].t_range);
 let ranges_b = self.collect_paveblock_ranges(be, self.ds.edges[be].t_range);

 if ranges_a.is_empty() || ranges_b.is_empty() {
 it.next(); continue;
 }

 if self.use_glue() && shared_edge_set.contains(&(ae, be)) {
 // Glue: use first pave point as shared vertex
 let pv = self.ds.edges[ae].start_vertex;
 if !self.ds.has_interf_ee(ae, be) {
 self.ds.interf_ee.push(InterferenceEE{
 e1: ae, e2: be,
 point: self.ds.vertices[pv].point,
 param1: self.ds.edges[ae].t_range[0],
 param2: self.ds.edges[be].t_range[0],
 new_vertex: pv,
 });
 }
 } else {
 let mut _ee_modified: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for ra in &ranges_a {
 for rb in &ranges_b {
 self.intersect_ee(ae, be, *ra, *rb, &mut _ee_modified);
 }
 }
 }
 it.next();
 }
 }
 /// OCCT PaveFiller_3.cxx L580-640: CheckEdgeEdge
 pub(crate) fn intersect_ee(&mut self, e1: usize, e2: usize,
 range1: [f64; 2], range2: [f64; 2],
 modified: &mut std::collections::HashSet<usize>) {
 let edge1 = &self.ds.edges[e1];
 let edge2 = &self.ds.edges[e2];
 let tol = self.ee_tol(e1, e2);
 // Capture all edge data before mutable borrow
 let e1_curve = edge1.curve.clone();
 let e2_curve = edge2.curve.clone();
 drop(edge1);
 drop(edge2);
 // ✅ OCCT-aligned: FillShrunkData computes shrunk ranges for each pave block.
 // If shrunk_range fails (edge too short), skip this pair entirely
 // (=OCCT BOPAlgo_PaveFiller_3: !aPB->IsSplittable() → continue).
 let sr1 = match crate::inttools::curve_range::shrunk_range(
  &e1_curve, range1, tol, tol, tol) {
  Some(sr) => sr,
  None => return,
 };
 let sr2 = match crate::inttools::curve_range::shrunk_range(
  &e2_curve, range2, tol, tol, tol) {
  Some(sr) => sr,
  None => return,
 };

 // Compute intersections restricted to shrunk sub-ranges.
 let hits: Vec<(f64, f64, DVec3)> = match (&e1_curve, &e2_curve) {
 (Curve3::Line(l1), Curve3::Line(l2)) => {
 intersect_line_line(l1, sr1, l2, sr2, tol)
 .into_iter().map(|(t1, t2, p)| (t1, t2, p)).collect()
 }
 (Curve3::Line(l), Curve3::Circle(c)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, sr1, tol) && in_range(*t_circle, sr2, tol)
 })
 .map(|(t_line, t_circle, p)| (t_line, t_circle, p))
 .collect(),
 (Curve3::Circle(c), Curve3::Line(l)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, sr2, tol) && in_range(*t_circle, sr1, tol)
 })
 .map(|(t_line, t_circle, p)| (t_circle, t_line, p))
 .collect(),
  (Curve3::Circle(c1), Curve3::Circle(c2)) => intersect_circle_circle(c1, c2, tol)
  .into_iter()
  .filter_map(|p| {
  let t1 = circle_param(p, c1);
  let t2 = circle_param(p, c2);
  if in_range(t1, sr1, tol) && in_range(t2, sr2, tol) {
  Some((t1, t2, p))
  } else { None }
  })
  .collect(),
  _ => {
  // Fallback: use EdgeEdgeIntersector for non-analytic curve pairs
  let mut ee = crate::inttools::edge_edge::EdgeEdgeIntersector::new();
  ee.set_edges(e1, sr1, e2, sr2, self.ds);
  ee.set_fuzzy_value(tol);
  ee.perform();
  ee.common_parts().iter().map(|cp| {
  let t1 = cp.vertex_param1;
  let t2 = cp.vertex_param2;
  Some((t1, t2, cp.bounding_point1))
  }).flatten().collect()
  }
 };

 //  OCCT-aligned: Process each intersection result (PaveFiller_3.cxx L682-750).
 // For each valid intersection, create a new vertex and record EE interference.
 // OCCT's UpdateVertex handles proximity via tolerance merging; rcad creates
 // vertices directly (architecture diff: rcad DSVertex has no UpdateVertex).
 for (t1, t2, point) in hits {
 // ✅ OCCT-aligned: restrict to shrunk range.  IntTools_EdgeEdge computes
 // within the shrunk range; results at/outside the boundary are
 // endpoint-coincident (handled by VV/VE/VF) or coincide with an existing
 // pave vertex — neither should create a new EE interference.
 if t1 < sr1[0] || t1 > sr1[1] || t2 < sr2[0] || t2 > sr2[1] { continue; }
 // ✅ OCCT-aligned: skip tangent/colinear edge pairs.  OCCT
 // (PaveFiller_3.cxx) checks aEECP.TangentEdges() after EE computation
 // and defers the entire range to ForceInterfEE (CommonBlocks).
 let is_parallel = match (&e1_curve, &e2_curve) {
  (Curve3::Line(l1), Curve3::Line(l2)) =>
   l1.direction.cross(l2.direction).length() <= tol,
  _ => false,
 };
 if is_parallel { continue; }
 let new_v = self.ds.add_vertex(point);
 let new_v = self.ds.add_vertex(point);
 self.ds.interf_ee.push(InterferenceEE{
 e1, e2, point, param1: t1, param2: t2, new_vertex: new_v,
 });
 self.ds.edges[e1].paves.push(Pave { vertex_idx: new_v, param: t1 });
 self.ds.edges[e2].paves.push(Pave { vertex_idx: new_v, param: t2 });
 modified.insert(e1);
 modified.insert(e2);
 }
 }
 /// OCCT PaveFiller_3.cxx L580-640: CheckEdgeEdge
 pub(crate) fn check_edge_edge(&mut self, e1: usize, e2: usize) {
 let edge1 = &self.ds.edges[e1];
 let edge2 = &self.ds.edges[e2];
 let tol = self.ee_tol(e1, e2);

 let hits: Vec<(f64, f64, DVec3)> = match (&edge1.curve, &edge2.curve) {
 (Curve3::Line(l1), Curve3::Line(l2)) => {
 intersect_line_line(l1, edge1.t_range, l2, edge2.t_range, tol)
 .into_iter()
 .map(|(t1, t2, p)| (t1, t2, p))
 .collect()
 }
 (Curve3::Line(l), Curve3::Circle(c)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, edge1.t_range, tol)
 && in_range(*t_circle, edge2.t_range, tol)
 })
 .map(|(t_line, t_circle, p)| (t_line, t_circle, p))
 .collect(),
 (Curve3::Circle(c), Curve3::Line(l)) => intersect_line_circle(l, c, tol)
 .into_iter()
 .filter(|(t_line, t_circle, _)| {
 in_range(*t_line, edge2.t_range, tol)
 && in_range(*t_circle, edge1.t_range, tol)
 })
 .map(|(t_line, t_circle, p)| (t_circle, t_line, p))
 .collect(),
 (Curve3::Circle(c1), Curve3::Circle(c2)) => intersect_circle_circle(c1, c2, tol)
 .into_iter()
 .filter_map(|p| {
 let t1 = circle_param(p, c1);
 let t2 = circle_param(p, c2);
 if in_range(t1, edge1.t_range, tol) && in_range(t2, edge2.t_range, tol) {
 Some((t1, t2, p))
 } else {
 None
 }
 })
 .collect(),
 _ => vec![],
 };

 for (t1, t2, point) in hits {
 let new_v = self.ds.add_vertex(point);
 self.ds.interf_ee.push(InterferenceEE{
 e1,
 e2,
 point,
 param1: t1,
 param2: t2,
 new_vertex: new_v,
 });
 self.ds.edges[e1].paves.push(Pave {
 vertex_idx: new_v,
 param: t1,
 });
 self.ds.edges[e2].paves.push(Pave {
 vertex_idx: new_v,
 param: t2,
 });
 }
 }
 /// OCCT PaveFiller L575-590: TreatNewVertices (merge EE/EF new vertices)
 pub(crate) fn treat_new_vertices(&mut self) -> Vec<usize> {
 // = =  Phase 1: Collect new vertices (OCCT L696-702) = = = = = = = = = = = = = = = 
 #[derive(Clone, Copy)]
 struct NewVertInfo { idx: usize, pos: DVec3, tol: f64 }
 let mut new_verts: Vec<NewVertInfo> = Vec::new();
 let mut seen = std::collections::BTreeSet::new();
 for inf in &self.ds.interf_ee {
 if seen.insert(inf.new_vertex) {
 let v_tol = self.ds.vertices[inf.new_vertex].geom_tol.max(self.ds.fuzzy_tol);
 new_verts.push(NewVertInfo { idx: inf.new_vertex, pos: inf.point, tol: v_tol });
 }
 }
 for inf in &self.ds.interf_ef {
 if seen.insert(inf.new_vertex) {
 let v_tol = self.ds.vertices[inf.new_vertex].geom_tol.max(self.ds.fuzzy_tol);
 new_verts.push(NewVertInfo { idx: inf.new_vertex, pos: inf.point, tol: v_tol });
 }
 }
 if new_verts.len() < 2 { return vec![]; }

 // = =  Phase 2: IntersectVertices (BOPAlgo_Tools.hxx L1119-1205) = = =
 let gap = self.ds.fuzzy_tol / 2.0;

 use crate::bvh::{Aabb, DsBvh};
 let nv = new_verts.len();
 let mut bvh_indices: Vec<usize> = Vec::with_capacity(nv);
 let mut bvh_aabbs: Vec<Aabb> = Vec::with_capacity(nv);
 for i in 0..nv {
 let v = &new_verts[i];
 bvh_indices.push(i);
 let half = v.tol + gap;
 bvh_aabbs.push(Aabb {
 min: v.pos - DVec3::splat(half),
 max: v.pos + DVec3::splat(half),
 });
 }
 let bvh = DsBvh::build(bvh_indices, bvh_aabbs);

 let pairs = DsBvh::candidate_pairs(&bvh, &bvh);
 let mut a_mili: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
 for &(ia, ib) in &pairs {
   let geom_vi = self.ds.vertices[new_verts[ia].idx].geom_tol;
   let geom_vj = self.ds.vertices[new_verts[ib].idx].geom_tol;
   let merge_tol = geom_vi + geom_vj + self.ds.fuzzy_tol;
   if (new_verts[ia].pos - new_verts[ib].pos).length() <= merge_tol {
     fill_map(&mut a_mili, ia, ib);
   }
 }
 let a_blocks = make_blocks(&a_mili);
 let mut groups: Vec<Vec<usize>> = a_blocks.iter()
   .map(|block| block.iter().map(|&i| new_verts[i].idx).collect())
   .collect();
 {
   let mut taken: HashSet<usize> = HashSet::new();
   for group in &groups {
     for &vi in group {
       taken.insert(vi);
     }
   }
   for i in 0..nv {
     if !taken.contains(&new_verts[i].idx) {
       groups.push(vec![new_verts[i].idx]);
     }
   }
 }

 // = =  Phase 3: MakeVertex for each chain (BOPTools_AlgoTools::MakeVertex) = = =
 // Single element -> reuse the vertex.
 // Multiple elements -> BRepLib::BoundingVertex computes center + tolerance.
 let mut survivors: Vec<usize> = Vec::new();
 for members in &groups {
 if members.len() < 2 {
 survivors.push(members[0]);
 continue;
 }

 // OCCT: BRepLib::BoundingVertex computes center + tolerance.
 // rcad: compute centroid of all member vertex positions.
 let centroid = members.iter()
 .map(|&vi| self.ds.vertices[vi].point)
 .sum::<DVec3>() / members.len() as f64;
 let max_tol = members.iter()
 .map(|&vi| self.ds.vertices[vi].geom_tol)
 .max_by(|a, b| a.partial_cmp(b).unwrap())
 .unwrap_or(self.ds.fuzzy_tol);

 // OCCT BRep_Builder::MakeVertex: create new vertex at centroid.
 let new_vi = self.ds.add_vertex(centroid);
 self.ds.vertices[new_vi].geom_tol = max_tol;
 // OCCT: myIncreasedSS.Add(nV)  ?mark tolerance as increased.
 self.ds.increased_ss.insert(new_vi);

 // Update interferences and paves to point to the new vertex
 for &old_vi in members {
 if old_vi == new_vi { continue; }
 for edge in &mut self.ds.edges {
 for pave in &mut edge.paves {
 if pave.vertex_idx == old_vi { pave.vertex_idx = new_vi; }
 }
 }
 for inf in &mut self.ds.interf_ee {
 if inf.new_vertex == old_vi { inf.new_vertex = new_vi; }
 }
 for inf in &mut self.ds.interf_ef {
 if inf.new_vertex == old_vi { inf.new_vertex = new_vi; }
 }
 for face in &mut self.ds.faces {
 if face.face_info.vertices_on.remove(&old_vi) {
 face.face_info.vertices_on.insert(new_vi);
 }
 if face.face_info.vertices_in.remove(&old_vi) {
 face.face_info.vertices_in.insert(new_vi);
 }
 }
 }
 survivors.push(new_vi);
 }
 survivors
 }
 /// OCCT PaveFiller_5.cxx L359-420: RepeatIntersection
 pub(crate) fn repeat_intersection(&mut self) {
 if self.ds.increased_ss.is_empty() { return; }
 let candidates: Vec<usize> = self.ds.increased_ss.iter().copied().collect();

 // Build set of existing interferences for dedup
 //  OCCT L398-413: PerformVV  ?PerformVE  ?PerformVF
 use std::collections::BTreeSet;
 let mut ve_done: BTreeSet<(usize, usize)> = BTreeSet::new();
 let mut vf_done: BTreeSet<(usize, usize)> = BTreeSet::new();
 for inf in &self.ds.interf_ve {
 ve_done.insert((inf.vertex, inf.edge));
 }
 for inf in &self.ds.interf_vf {
 vf_done.insert((inf.vertex, inf.face));
 }

 // = =  VV: check survivors against vertices on the other side = = = = = = 
 //  OCCT L398: PerformVV(aPS.Next())
 // VV safe: if pair already in interferences, add_vertex will dedup
 for &vi in &candidates {
 let vi_origin = self.ds.vertices[vi].origin;
 let other_verts: Vec<usize> = self.ds.vertices.iter().enumerate()
 .filter(|(j, v)| {
 if *j == vi { return false; }
 match (vi_origin, v.origin) {
 (Some(ShapeOrigin::ShapeA), Some(ShapeOrigin::ShapeB)) => true,
 (Some(ShapeOrigin::ShapeB), Some(ShapeOrigin::ShapeA)) => true,
 _ => false,
 }
 })
 .map(|(j, _)| j)
 .collect();
 for &vj in &other_verts {
 let tol = self.vv_pair_tol(vi, vj);
 let dist = (self.ds.vertices[vi].point - self.ds.vertices[vj].point).length();
 if dist <= tol {
 self.ds.interf_vv.push(InterferenceVV{
 v1: vi, v2: vj, merged_vertex: vi,
 });
 }
 }
 }

 // = =  VE: check survivors against edges on the other side = = = = = = 
 //  OCCT L403: PerformVE(aPS.Next())
 for &vi in &candidates {
 let vi_origin = self.ds.vertices[vi].origin;
 let other_edges: Vec<usize> = match vi_origin {
 Some(ShapeOrigin::ShapeA) => self.edges_of(ShapeOrigin::ShapeB),
 Some(ShapeOrigin::ShapeB) => self.edges_of(ShapeOrigin::ShapeA),
 _ => continue,
 };
 for &ei in &other_edges {
 if ve_done.contains(&(vi, ei)) { continue; }
 self.check_vertex_edge(vi, ei);
 }
 }

 // = =  VF: check survivors against faces on the other side = = = = = = 
 //  OCCT L408: PerformVF(aPS.Next())
 for &vi in &candidates {
 let vi_origin = self.ds.vertices[vi].origin;
 let other_faces: Vec<usize> = match vi_origin {
 Some(ShapeOrigin::ShapeA) => self.faces_of(ShapeOrigin::ShapeB),
 Some(ShapeOrigin::ShapeB) => self.faces_of(ShapeOrigin::ShapeA),
 _ => continue,
 };
 for &fi in &other_faces {
 if vf_done.contains(&(vi, fi)) { continue; }
 self.check_vertex_face(vi, fi);
 }
 }
 }
 /// OCCT PaveFiller_4.cxx: PerformVF
 pub(crate) fn perform_vf(&mut self) {
 // OCCT PaveFiller_4.cxx: FillShrunkData + BVH pair iteration
 // with HasInterf skip condition.
 self.fill_shrunk_data(); // OCCT: FillShrunkData(VERTEX, FACE)
 let a_verts = self.verts_of(ShapeOrigin::ShapeA);
 let b_faces = self.faces_of(ShapeOrigin::ShapeB);
 for &vi in &a_verts {
 for &fi in &b_faces {
 if self.ds.has_interf_vf(vi, fi) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, fi) { continue; }
 self.check_vertex_face(vi, fi);
 }
 }
 let b_verts = self.verts_of(ShapeOrigin::ShapeB);
 let a_faces = self.faces_of(ShapeOrigin::ShapeA);
 for &vi in &b_verts {
 for &fi in &a_faces {
 if self.ds.has_interf_vf(vi, fi) { continue; }
 if self.ds.has_interf_ve_via_faces(vi, fi) { continue; }
 self.check_vertex_face(vi, fi);
 }
 }
 }

 ///  OCCT-aligned: CheckVertexFace (PaveFiller_4.cxx L249-298).
 /// Vertex/Face proximity check with SD vertex resolution.
 /// OCCT: BOPAlgo_VertexFace parallel solver + result processing;
 /// rcad: sequential equivalent with same projection logic.
 pub(crate) fn check_vertex_face(&mut self, vi: usize, fi: usize) {
 let n_vsd = self.ds.has_shape_sd(vi).unwrap_or(vi);
 let point = self.ds.vertices[n_vsd].point;
 let face = &self.ds.faces[fi];
 let tf = self.vf_tol(n_vsd, fi);
    let (is_on, proj_dist, proj_u, proj_v): (bool, f64, f64, f64) = match &face.surface {
        Surface3::Plane(plane) => {
            if inttools::vertex_ops::vertex_on_plane_with_tol(point, plane, tf) {
                let face_verts = self.ds.face_boundary_points(fi);
                let on_face = inttools::edge_face::point_in_planar_face_with_tol(
                    point, plane, &face_verts, tf);
                // UV: project point onto plane UV coordinate system
                let u_axis = plane.u_dir;
                let v_axis = plane.v_dir;
                let diff = point - plane.origin;
                let u = diff.dot(u_axis);
                let v = diff.dot(v_axis);
                (on_face, 0.0, u, v)
            } else {
                (false, f64::MAX, 0.0, 0.0)
            }
        }
        surface => {
            let proj = rcad_kernel::projection::closest_point_on_surface(
                surface, point, 16);
            let a_tol_v = self.ds.vertices[n_vsd].geom_tol;
            let a_tol_f = face.geom_tol;
            let a_tol_sum = a_tol_v + a_tol_f + self.ds.fuzzy_tol.max(tf);
            if proj.distance <= a_tol_sum {
                let uv = DVec2::new(proj.params.0, proj.params.1);
                let fclass = FClass2d::new(self.ds, fi, tf);
                let inside = fclass.perform(uv, false) == State::In;
                (inside, proj.distance, proj.params.0, proj.params.1)
            } else {
                (false, f64::MAX, 0.0, 0.0)
            }
        }
    };

    if is_on {
        self.ds.interf_vf.push(InterferenceVF{
            vertex: n_vsd,
            face: fi,
            u: proj_u,
            v: proj_v,
        });
 if proj_dist > 0.0 && proj_dist < f64::MAX
 && proj_dist > self.ds.vertices[n_vsd].geom_tol
 {
 self.ds.vertices[n_vsd].geom_tol = proj_dist;
 self.ds.increased_ss.insert(n_vsd);
 }

 //  OCCT-aligned: ALL VF vertices go to VerticesIn (OCCT L297: aMVIn.Add)
 self.ds.faces[fi].face_info.vertices_in.insert(n_vsd);
 }
 }
 ///  OCCT-aligned: PerformEF (PaveFiller_5.cxx L165-300).
 /// Iterates (edge, face) pairs with HasFlag / HasInterf skip conditions.
 /// Uses full edge range (not sub-ranges)  ?matching OCCT's original PB iteration.
 pub(crate) fn perform_ef(&mut self) {
 self.fill_shrunk_data();
 let a_edges = self.edges_of(ShapeOrigin::ShapeA);
 let b_faces = self.faces_of(ShapeOrigin::ShapeB);

 for &ei in &a_edges {
 if self.ds.edge_has_flag(ei) { continue; }
 if self.ds.is_edge_degenerated(ei) { continue; }
 let etr = self.ds.edges[ei].t_range;
 for &fi in &b_faces {
 if self.ds.has_interf_ef(ei, fi) { continue; }
 self.intersect_ef(ei, fi, &etr);
 }
 }

 let b_edges = self.edges_of(ShapeOrigin::ShapeB);
 let a_faces = self.faces_of(ShapeOrigin::ShapeA);

 for &ei in &b_edges {
 if self.ds.edge_has_flag(ei) { continue; }
 if self.ds.is_edge_degenerated(ei) { continue; }
 let etr = self.ds.edges[ei].t_range;
 for &fi in &a_faces {
 if self.ds.has_interf_ef(ei, fi) { continue; }
 self.intersect_ef(ei, fi, &etr);
 }
 }
 }
 /// OCCT PaveFiller_3.cxx L222-228: GetPBBox (PaveBlock range)
 pub(crate) fn collect_paveblock_ranges(&self, edge_idx: usize, edge_t_range: [f64; 2]) -> Vec<[f64; 2]> {
 let paves = &self.ds.edges[edge_idx].paves;
 if paves.is_empty() {
 return vec![edge_t_range];
 }
 let mut params: Vec<f64> = paves.iter().map(|p| p.param).filter(|p| p.is_finite()).collect();
 params.sort_by(|a, b| a.partial_cmp(b).unwrap());
 let edge_tol = self.ds.edges[edge_idx].geom_tol.max(self.tol());
 // Deduplicate
 params.dedup_by(|a, b| (*a - *b).abs() < edge_tol);
 // Include endpoints
 let mut bounds = vec![edge_t_range[0]];
 bounds.extend(params);
 bounds.push(edge_t_range[1]);
 // Build ranges
 let mut ranges = Vec::new();
 for w in bounds.windows(2) {
 if w[1] - w[0] > edge_tol {
 ranges.push([w[0], w[1]]);
 }
 }
 ranges
 }
 /// OCCT: shrunk range correction for face tolerance
 pub(crate) fn correct_range_for_face(edge_curve: &Curve3, etf: f64, range: [f64; 2]) -> [f64; 2] {
 const DT: f64 = 1e-12;
 match edge_curve {
 Curve3::Line(_) => range,
 Curve3::Circle(c) => {
 let a_res = etf / c.radius.max(TOLERANCE_ABS);
 let new_first = range[0] + a_res;
 let new_last = range[1] - a_res;
 if new_last - new_first < DT { range } else { [new_first, new_last] }
 }
 Curve3::Ellipse(e) => {
 let a_res = etf / e.major_radius.max(TOLERANCE_ABS);
 let new_first = range[0] + a_res;
 let new_last = range[1] - a_res;
 if new_last - new_first < DT { range } else { [new_first, new_last] }
 }
 _ => {
 let new_first = range[0] + etf;
 let new_last = range[1] - etf;
 if new_last - new_first < DT { range } else { [new_first, new_last] }
 }
 }
 }
 /// OCCT IntTools_FClass2d: point-in-face check
 pub(crate) fn is_point_in_face(&self, point: DVec3, face_idx: usize, tol: f64) -> bool {
 let face = &self.ds.faces[face_idx];
 match &face.surface {
 Surface3::Plane(plane) => {
 let verts = self.ds.face_boundary_points(face_idx);
 inttools::edge_face::point_in_planar_face_with_tol(point, plane, &verts, tol)
 }
 Surface3::Sphere(sphere) => {
 let uv = sphere.world_to_uv(point);
 uv.x >= -tol && uv.x <= std::f64::consts::TAU + tol
 && uv.y >= -std::f64::consts::FRAC_PI_2 - tol
 && uv.y <= std::f64::consts::FRAC_PI_2 + tol
 }
 Surface3::Cylinder(cyl) => {
 let axis = cyl.axis.normalize_or_zero();
 if axis.length_squared() < 0.5 { return true; }
 let local = point - cyl.origin;
 let v = local.dot(axis);
 let radial = local - axis * v;
 let u = radial.y.atan2(radial.x);
 let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
 u >= -tol && u <= std::f64::consts::TAU + tol
 }
 Surface3::Cone(cone) => {
 let axis = cone.axis_dir();
 let ap = point - cone.apex;
 let v = ap.dot(axis);
 let radial = ap - axis * v;
 let u = radial.y.atan2(radial.x);
 let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
 u >= -tol && u <= std::f64::consts::TAU + tol
 }
 Surface3::Torus(torus) => {
 let axis = torus.axis.normalize_or_zero();
 if axis.length_squared() < 0.5 { return true; }
 let local = point - torus.center;
 let v = local.dot(axis);
 let radial = local - axis * v;
 let r = radial.length();
 let u = radial.y.atan2(radial.x);
 let u = if u < 0.0 { u + std::f64::consts::TAU } else { u };
 let v_angle = ((r - torus.major_radius) / torus.minor_radius.max(tol)).acos();
 let v = if (r - torus.major_radius).abs() <= torus.minor_radius + tol { v_angle } else { 0.0 };
 u >= -tol && u <= std::f64::consts::TAU + tol
 && v >= -tol && v <= std::f64::consts::TAU + tol
 }
 _ => true,
 }
 }
 /// OCCT PaveFiller_5.cxx L340-480: IntersectEdgeFace
 pub(crate) fn intersect_ef(&mut self, edge_idx: usize, face_idx: usize, pb_range: &[f64; 2]) {
 let edge_curve = self.ds.edges[edge_idx].curve.clone();
 let edge_t_range = self.ds.edges[edge_idx].t_range;

 // Use PaveBlock range to constrain intersection interval (OCCT L262: SetRange(aPBRange))
 let ef_range = [
 pb_range[0].max(edge_t_range[0]),
 pb_range[1].min(edge_t_range[1]),
 ];
 let etf = self.ef_tol(edge_idx, face_idx);
 if ef_range[1] - ef_range[0] <= etf {
 return;
 }
 let ef_range = Self::correct_range_for_face(&edge_curve, etf, ef_range);
 if ef_range[1] - ef_range[0] <= etf {
 return;
 }
 let face_surface = self.ds.faces[face_idx].surface.clone();

 // Dispatch based on curve type  ?surface type
 let hits: Vec<(DVec3, f64)> = match (&edge_curve, &face_surface) {
 (Curve3::Line(line), Surface3::Plane(plane)) => {
 inttools::edge_face::intersect_line_plane_with_tol(
 line,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.edge_param))
 .collect()
 }
 (Curve3::Line(line), Surface3::Cylinder(cyl)) => {
 inttools::curve_surface::intersect_line_cylinder_with_tol(
 line,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Line(line), Surface3::Sphere(sph)) => {
 inttools::curve_surface::intersect_line_sphere_with_tol(
 line,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Line(line), Surface3::Cone(cone)) => {
 inttools::curve_surface::intersect_line_cone_with_tol(
 line,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Circle(circle), Surface3::Plane(plane)) => {
 // Use edge start vertex as reference direction for  ?0
 let sv = self.ds.edges[edge_idx].start_vertex;
 let ref_dir = (self.ds.vertices[sv].point - circle.center).normalize();
 inttools::curve_surface::intersect_circle_plane_with_ref(
 circle, ef_range, plane, etf, Some(ref_dir),
 )
 .into_iter().map(|h| (h.point, h.curve_param)).collect()
 }
 (Curve3::Circle(circle), Surface3::Cylinder(cyl)) => {
 inttools::curve_surface::intersect_circle_cylinder_with_tol(
 circle,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Circle(circle), Surface3::Sphere(sph)) => {
 inttools::curve_surface::intersect_circle_sphere_with_tol(
 circle,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Circle(circle), Surface3::Cone(cone)) => {
 inttools::curve_surface::intersect_circle_cone_with_tol(
 circle,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.curve_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Plane(plane)) => {
 //  OCCT-aligned: IntAna_IntConicQuad Ellipse  ?Plane
 inttools::ellipse_intersection::intersect_ellipse_plane_with_tol(
 ellipse,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Cylinder(cyl)) => {
 //  ?Partially aligned: numeric fallback, same as OCCT for rare cases
 inttools::ellipse_intersection::intersect_ellipse_cylinder_with_tol(
 ellipse,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Sphere(sph)) => {
 inttools::ellipse_intersection::intersect_ellipse_sphere_with_tol(
 ellipse,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Ellipse(ellipse), Surface3::Cone(cone)) => {
 inttools::ellipse_intersection::intersect_ellipse_cone_with_tol(
 ellipse,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.ellipse_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Plane(plane)) => {
 //  OCCT-aligned: IntAna_IntConicQuad Parabola  ?Plane
 inttools::parabola_intersection::intersect_parabola_plane_with_tol(
 parabola,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Cylinder(cyl)) => {
 //  ?Partially aligned: numeric fallback
 inttools::parabola_intersection::intersect_parabola_cylinder_with_tol(
 parabola,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Sphere(sph)) => {
 inttools::parabola_intersection::intersect_parabola_sphere_with_tol(
 parabola,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Parabola(parabola), Surface3::Cone(cone)) => {
 inttools::parabola_intersection::intersect_parabola_cone_with_tol(
 parabola,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.parabola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Plane(plane)) => {
 //  OCCT-aligned: IntAna_IntConicQuad Hyperbola  ?Plane
 inttools::hyperbola_intersection::intersect_hyperbola_plane_with_tol(
 hyperbola,
 ef_range,
 plane,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Cylinder(cyl)) => {
 //  ?Partially aligned: numeric fallback
 inttools::hyperbola_intersection::intersect_hyperbola_cylinder_with_tol(
 hyperbola,
 ef_range,
 cyl,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Sphere(sph)) => {
 inttools::hyperbola_intersection::intersect_hyperbola_sphere_with_tol(
 hyperbola,
 ef_range,
 sph,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 (Curve3::Hyperbola(hyperbola), Surface3::Cone(cone)) => {
 inttools::hyperbola_intersection::intersect_hyperbola_cone_with_tol(
 hyperbola,
 ef_range,
 cone,
 etf,
 )
 .into_iter()
 .map(|h| (h.point, h.hyperbola_param))
 .collect()
 }
 _ => {
 // Numeric fallback: sample the curve, find sign changes of the
 // surface implicit function. Works for any Curve3  ?Surface3 pair.
 intersect_edge_face_numeric(&edge_curve, &face_surface, ef_range, etf)
 }
 };

 for (point, edge_param) in hits {
 //  OCCT-aligned: IsPointInFace check for ALL surface types (PaveFiller_5.cxx L523)
 let in_face = self.is_point_in_face(point, face_idx, etf);
 if !in_face {
 let near_face_vert = match &face_surface {
 Surface3::Plane(_) => {
 self.ds.face_boundary_points(face_idx).iter().any(|&vp| {
 (vp - point).length() <= etf
 })
 }
 _ => false,
 };
 if !near_face_vert { continue; }
 }

 //  OCCT-aligned: Always create EF interference for intersection hits.
 // OCCT IntTools_EdgeFace creates a new vertex for each hit, even when
 // the hit coincides with an existing edge endpoint.  SD vertex merging
 // handles near-coincident vertices later (MakeSDVerticesFF in PostTreat).
 // rcad: do NOT skip endpoint-coincident hits  ?they are needed for
 // PutPaveOnCurve to split intersection curve pave blocks.
 let new_v = self.ds.add_vertex(point);
 // Register vertices_on for the new vertex if it's near the edge boundary
 let sv = self.ds.edges[edge_idx].start_vertex;
 let ev = self.ds.edges[edge_idx].end_vertex;
 let tol = etf
 .max(self.ds.vertices[sv].geom_tol)
 .max(self.ds.vertices[ev].geom_tol);
 if (point - self.ds.vertices[sv].point).length() <= tol
 || (point - self.ds.vertices[ev].point).length() <= tol
 {
 self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
 }
 //  OCCT-aligned: Create EF interference for EVERY hit, even at edge endpoints.
 // OCCT IntTools_EdgeFace creates a new vertex for each hit (no dedup).
 // rcad: remove the vertices_on skip check  ?always push interference.
 self.ds.interf_ef.push(InterferenceEF{
 edge: edge_idx,
 face: face_idx,
 point,
 edge_param,
 new_vertex: new_v,
 });
 // Only mark vertices_on if actually inserted (avoid duplicate insert msg)
 if !self.ds.faces[face_idx].face_info.vertices_on.contains(&new_v) {
 self.ds.faces[face_idx].face_info.vertices_on.insert(new_v);
 }
 self.ds.edges[edge_idx].paves.push(Pave {
 vertex_idx: new_v,
 param: edge_param,
 });
 }
 }
 /// OCCT HasInterf: skip already-processed pairs
 pub fn skip_redundant_interferences(&self) -> std::collections::HashSet<(usize, usize, u8)> {
 let mut skip_set = std::collections::HashSet::new();

 if !self.use_glue() {
 return skip_set;
 }

 // Skip V-V for shared vertices
 for &(va, vb) in &self.ds.shared_topology.shared_vertices {
 skip_set.insert((va, vb, 0)); // 0 = V-V
 }

 // Skip E-E for shared edges
 for &(ea, eb) in &self.ds.shared_topology.shared_edges {
 skip_set.insert((ea, eb, 2)); // 2 = E-E
 }

 // Skip F-F for fully glued faces
 for &(fa, fb) in &self.ds.shared_topology.fully_glued_faces {
 skip_set.insert((fa, fb, 5)); // 5 = F-F
 }

 skip_set
 }

 /// OCCT myDS->HasInterf: check existing EE interference
 pub(crate) fn has_ee_interf(&self, e1: usize, e2: usize) -> bool {
 self.ds.interf_ee.iter().any(|inf| {
 (inf.e1 == e1 && inf.e2 == e2) || (inf.e1 == e2 && inf.e2 == e1)
 })
 }
}

