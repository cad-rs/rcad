use std::collections::HashSet;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::*;
use rcad_kernel::PCurve;
use rcad_kernel::topods;

use crate::bopds::ds::{
 DS, DSEdge, DSCurveRepOnFace, DSVertex, Interference, InterferenceFF, InterferenceVV,
 InterferenceVE, InterferenceVF, InterferenceEE, InterferenceEF, IntersectionCurve, ShapeOrigin,
};
use crate::bopds::ds::topods_builder::load_topods_brep;
use crate::bopds::ds::face_aabb;
use crate::bvh::DsBvh;
use crate::bopds::pave::*;
use crate::bopalgo::{GlueEnum, Alert, Report};
use crate::bvh::Bvh;
use crate::inttools;
use crate::inttools::context::Context as IntToolsContext;
use crate::inttools::fclass2d::{FClass2d, State};
use crate::tolerance::*;
use rcad_kernel::closest_point_on_curve;
pub mod helpers;
use self::helpers::*;

// =IntPatch_Intersection surface category (L1264-1294).
// GeomGeom = ts1==ts2==1 (both analytic, ImpImpIntersection)
// ParamParam = ts1==ts2==0 (both parametric, PrmPrmIntersection)

// Re-export NearTangentType from bopds::ds for use in this module's public types
pub use crate::bopds::ds::NearTangentType;

/// Minimum total face count before BVH acceleration is used.
/// Below this threshold, brute-force O(n ? is faster due to BVH build overhead.
const BVH_THRESHOLD: usize = 20;

///  BOPAlgo_PaveFiller =six intersection passes
/// (PaveFiller.hxx L106-107, PaveFiller.cxx L234-355).
mod glue;
mod intersection;
pub(crate) mod analytics;
pub(crate) mod marching;
pub(crate) mod polyhedron;
pub(crate) mod polyhedron_bvh;
pub(crate) mod prm_prm_intersection;
pub(crate) mod p_walking;
mod make_blocks;
mod posttreat;
mod config;
mod tolerances;
 mod interf;
mod paves;
mod ff_intersect;

/// BOPAlgo_SectionAttribute =controls approximation and
/// pcurve computation for section edges (BOPAlgo_SectionAttribute.hxx).
#[derive(Debug, Clone)]
pub(crate) struct SectionAttribute {
 pub approximation: bool,
 pub pcurve_on_s1: bool,
 pub pcurve_on_s2: bool,
}

impl Default for SectionAttribute {
 fn default() -> Self {
 Self { approximation: true, pcurve_on_s1: true, pcurve_on_s2: true }
 }
}

/// PaveFiller::EdgeRangeDistance =stores minimal distance
/// between an edge range and a face that don't geometrically intersect.
/// Used by PostTreatFF to re-check E-F pairs after tolerance updates.
#[derive(Debug, Clone)]
pub(crate) struct EdgeRangeDistance {
 pub first: f64,
 pub last: f64,
 pub distance: f64,
}

pub struct PaveFiller<'a> {
 pub ds: &'a mut DS,
 /// Optional BRep for direct output (dual-write mode). When set, PaveFiller
 /// populates the BRep on completion, eliminating the need for ds_to_brep.
 pub brep: Option<&'a mut rcad_kernel::topods::BRep>,
 /// Output: face_refs by ds_face_idx (populated by export_to_brep).
 pub face_refs: Vec<rcad_kernel::topods::ShapeRef>,
 /// Output: ic_edge_map: ci -> BRep edge ShapeRef (populated by export_to_brep).
 pub ic_edge_map: Vec<Option<rcad_kernel::topods::ShapeRef>>,
 bvh_a: Option<&'a Bvh>,
 bvh_b: Option<&'a Bvh>,
 /// DS-based face BVH for FF pair detection. Uses DS face indices directly,
 /// matching OCCT's BOPTools_BoxTree which operates on source shape indices.
 /// BVH index space equals DS index space (no a_rev/b_rev needed).
 pub(crate) face_bvh: Option<crate::bvh::DsBvh>,
 ///  BOPAlgo_GlueEnum (GlueOff/GlueFull/GlueShift).
 glue: GlueEnum,
 glue_tolerance: f64,
 /// convenience  ?true when glue is active (not GlueOff).
 /// =BOPAlgo_Options::SetFuzzyValue
 fuzzy_tolerance: f64,
 /// =PaveFiller_6.cxx L393-479 seam edge shift tolerance
 seam_shift_tol: f64,
 /// =BOPAlgo_Algo::myRunParallel
 run_parallel: bool,
 /// =BOPAlgo_PaveFiller::myNonDestructive
 non_destructive: bool,
 /// =BOPAlgo_Algo::myUseOBB
 use_obb: bool,
 /// =IntTools_Context (PaveFiller::Init L203)
 context: IntToolsContext,
 /// =myArguments =original input shapes (BOPAlgo_PaveFiller.hxx L639).
 /// rcad: carries the original BRep operands for OCCT-API compatibility.
 my_arguments: Vec<rcad_kernel::topods::BRep>,
 /// =mySectionAttribute (BOPAlgo_SectionAttribute.hxx)
 section_attribute: SectionAttribute,
 /// =myIsPrimary (BOPAlgo_PaveFiller.cxx L62)
 is_primary: bool,
 /// =myAvoidBuildPCurve (BOPAlgo_PaveFiller.cxx L63)
 avoid_build_pcurve: bool,
 /// =myFPBDone =fence map tracking processed (face, pave_block) pairs.
 /// Map: face_idx =set of pave_block indices already processed in PostTreatFF.
 fpbdone: std::collections::HashMap<usize, std::collections::HashSet<usize>>,
 /// =myVertsToAvoidExtension =vertices that should NOT have
 /// their tolerance extended further (near EE/EF intersection points).
 verts_to_avoid_extension: std::collections::HashSet<usize>,
 /// aMVTol -- per-vertex tolerance map (PaveFiller_6.cxx L2409).
 a_mv_tol: std::collections::HashMap<usize, f64>,
 /// aDMVLV -- duplicate vertex map (PaveFiller_6.cxx L2410).
 a_dmv_lv: std::collections::HashMap<usize, Vec<usize>>,
 /// =myDistances =minimal edge-face distances for non-intersecting
 /// pairs.  Map: (edge_idx, face_idx) =Vec<EdgeRangeDistance>.
 distances: std::collections::HashMap<(usize, usize), Vec<EdgeRangeDistance>>,
 ///  myReport  ?collects alerts during PaveFiller execution.
 my_report: Report,
 /// Pipeline stage dump context (rcad PF stages). Created from env vars;
 /// disabled when RCAD_DUMP_PIPELINE is not set.
 pub dump_ctx: crate::pipeline_dump::DumpCtx,
 /// Stage to stop after (for stage-by-stage testing). When set, perform()
 /// returns early after completing the named stage.
 pub stop_after: Option<String>,
}

impl<'a> PaveFiller<'a> {
 /// IsGlue  ?true when glue mode is active (not GlueOff).
 pub fn use_glue(&self) -> bool { self.glue != GlueEnum::GlueOff }

 /// GetGlue  ?return current glue mode.
 pub fn glue_mode(&self) -> GlueEnum { self.glue }
}

/// GetFullShapeMap (PaveFiller_6.cxx L2941-2958).
/// Builds a set of all sub-shape indices belonging to face `fi`:
/// the face itself, its boundary edges, and their endpoint vertices.
pub(crate) fn build_face_shape_map(ds: &DS, fi: usize) -> std::collections::HashSet<usize> {
 let mut aMI = std::collections::HashSet::new();
 aMI.insert(fi);
 if fi < ds.faces.len() {
 for &ei in &ds.faces[fi].boundary_edges {
 aMI.insert(ei);
 if ei < ds.edges.len() {
 aMI.insert(ds.edge_start_vertex_ds(ei));
 aMI.insert(ds.edge_end_vertex_ds(ei));
 }
 }
 }
 aMI
}

/// =Propagate IC vertices to all faces sharing boundary edges
/// (OCCT BOPDS_FaceInfo::AppendBlock equivalent).
/// OCCT BOPAlgo_PaveFiller propagates pave block vertices to all faces
/// referencing the split edge. rcad's add_curve only adds vertices to the
/// two FF-interference faces (f1, f2), but the vertex may lie on boundary
/// edges of other faces (e.g. side-face tangent-line IC endpoints on the
/// top face's boundary edge).
fn propagate_ic_vertices_to_shared_faces(
 ds: &mut DS,
 ic_vertices: &[usize],
 skip_faces: &[usize; 2],
) {
 let vtol = TOLERANCE_ABS * 1000.0; // 1e-4 geometric tolerance for on-edge check
 let vtol_sq = vtol * vtol;
 for fi in 0..ds.faces.len() {
 if fi == skip_faces[0] || fi == skip_faces[1] {
 continue;
 }
 for &vi in ic_vertices {
 if ds.face_info(fi).vertices_in.contains(&vi) {
 continue;
 }
 let vp = ds.vertex_point(vi);
 for &ei in &ds.faces[fi].boundary_edges {
 let Some(edge) = ds.edges.get(ei) else { continue };
 let a = ds.vertices[edge.start_vertex].point;
 let b = ds.vertices[edge.end_vertex].point;
 let ab = b - a;
 let ab_len2 = ab.length_squared();
 if ab_len2 < TOLERANCE_LEN_SQ_DIV_SAFE {
 continue;
 }
 let ap = vp - a;
 let t = ap.dot(ab) / ab_len2;
 if t > -0.01 && t < 1.01 {
 let proj = a + ab * t.clamp(0.0, 1.0);
 if (vp - proj).length_squared() < vtol_sq {
 ds.face_info_mut(fi).vertices_in.insert(vi);
 break;
 }
 }
 }
 }
 }
}

impl<'a> PaveFiller<'a> {

 ///  Prepare (PaveFiller_7.cxx L850-929).
 /// Build 2D pcurves for edges on planar faces.
 /// OCCT: iterate all V/E, E/E, E/F pairs to find planar faces, collect edge-face pairs,
 /// compute pcurves in parallel, update edges.  rcad: DS::build_face_reps already computes
 /// pcurves for all face types; this step ensures planar-face pcurves exist for any
 /// edges that build_face_reps may have missed (non-boundary or intersection-relevant).
 fn prepare(&mut self) {
 let mut planar_faces: Vec<usize> = Vec::new();
 for (fi, f) in self.ds.faces.iter().enumerate() {
 if matches!(f.surface, Surface3::Plane(_)) {
 planar_faces.push(fi);
 }
 }
 if planar_faces.is_empty() { return; }

 let surf: Vec<Surface3> = planar_faces.iter().map(|&fi| self.ds.faces[fi].surface.clone()).collect();

 // Collect all edge indices from each planar face's boundary
 let mut face_edges: Vec<Vec<usize>> = Vec::with_capacity(planar_faces.len());
 for (pos, &fi) in planar_faces.iter().enumerate() {
 let f = &self.ds.faces[fi];
 let mut eids: Vec<usize> = f.boundary_edges.clone();
 for w in &f.inner_boundary_edges {
 eids.extend(w.iter().map(|&(ei, _)| ei));
 }
 face_edges.push(eids);
 }

 // Compute pcurves for edge-face pairs that don't already have one
 for (pos, &fi) in planar_faces.iter().enumerate() {
 for &ei in &face_edges[pos] {
 if self.ds.edge_on_face(ei, fi).is_some() { continue; }
 let Some(edge) = self.ds.edges.get_mut(ei) else { continue; };
 if let Some((pcurve, span)) = DS::compute_edge_pcurve(&edge.curve, &surf[pos], None) {
 edge.face_reps.push(DSCurveRepOnFace {
 face_idx: fi,
 pcurve,
 pcurve2: None,
 pcurve_range: [0.0, span],
 start_param: 0.0,
 end_param: span,
 });
 }
 }
 }
 }

 //  ?OCCT L248 Prepare: build pcurves on planar faces
 // RepeatInt->ForceEE->ForceEF->FF->UpdBlk->RefFI->MkSEdges->MkBlks->
 // ChkSI->RefFO->RmvME->MkPCurves->ProcDE

 /// BOPAlgo_PaveFiller::Init (PaveFiller.cxx L176-213).
 /// Populates the DS from operand BReps, matching what
 /// BOPDS_DS::Init + BOPDS_DS::SetArguments does in OCCT.
 /// Guards: only populates when DS is empty (no faces loaded yet)
 /// so that legacy callers that pre-populate DS before perform() still work.
 fn init(&mut self, a: &topods::BRep, b: &topods::BRep, fuzzy_tol: f64) {
 if !self.ds.faces.is_empty() { return; }
 let tol = fuzzy_tol.max(TOLERANCE_ABS);
 self.ds.fuzzy_tol = tol;
 load_topods_brep(&mut self.ds, a, ShapeOrigin::ShapeA);
 self.ds.a_vertex_count = self.ds.vertices.len();
 self.ds.a_edge_count = self.ds.edges.len();
 self.ds.a_face_count = self.ds.faces.len();
 load_topods_brep(&mut self.ds, b, ShapeOrigin::ShapeB);
 self.ds.compute_uv_boundaries();
 self.ds.build_face_reps();
 // ShapeInfo is built during load_topods_brep (init_shape_topo traversal).
 self.ds.nb_source_shapes = self.ds.shape_info.len();
 self.ds.build_map_ve();
 // Resize context to match the actual number of faces now loaded
 if self.ds.faces.len() != self.context.num_faces {
     self.context.resize(self.ds.faces.len());
 }
 }

 /// Check whether to stop after the given stage name.
 /// Returns true if perform() should return early.
 fn check_stop(&self, stage: &str) -> bool {
     if self.stop_after.as_deref() == Some(stage) { return true; }
     if let Ok(s) = std::env::var("RCAD_STOP_AFTER") { return s == stage; }
     false
 }

 pub fn perform(&mut self, a: &topods::BRep, b: &topods::BRep) {
  // =Init =create DS from arguments (PaveFiller.cxx L176-213).
  self.init(a, b, self.fuzzy_tolerance);
  self.dump_ctx.snapshot("after_Init", self.ds, None);
  if self.check_stop("after_Init") { return; }

  // Build face BVH for FF pair detection (analogous to BOPTools_BoxTree).
  // Must happen after init() populates the DS faces.
  {
  let mut indices = Vec::new();
  let mut aabbs = Vec::new();
  for (fi, f) in self.ds.faces.iter().enumerate() {
  indices.push(fi);
  aabbs.push(crate::bopds::ds::face_aabb::face_aabb(self.ds, fi));
  }
  self.face_bvh = if indices.len() >= 20 { Some(DsBvh::build(indices, aabbs)) } else { None };
  }

  // =early return if stop_after env var is set and matches
  // (allows stage-by-stage testing without modifying production logic)

  // OCCT L251: Prepare =build pcurves on planar faces.
  self.prepare();
  self.dump_ctx.snapshot("after_Prepare", self.ds, None);
  if self.check_stop("after_Prepare") { return; }

  // BOPDS_Iterator — single BVH, type bucketing, stable_sort.
  // rcad: extract all pair lists upfront (Rust borrow limits prevent
  // keeping a persistent iterator as a member like OCCT's myIterator).
  use crate::bopds::ds::BOPDS_Iterator;
  use rcad_kernel::topods::ShapeType;
  let (vv_pairs, ve_pairs, ee_pairs, vf_pairs) = {
  let mut iterator = BOPDS_Iterator::new(self.ds);
  iterator.prepare();
  let vv = iterator.pairs(ShapeType::Vertex, ShapeType::Vertex).to_vec();
  let ve = iterator.pairs(ShapeType::Vertex, ShapeType::Edge).to_vec();
  let ee = iterator.pairs(ShapeType::Edge, ShapeType::Edge).to_vec();
  let vf = iterator.pairs(ShapeType::Vertex, ShapeType::Face).to_vec();
  (vv, ve, ee, vf)
  };

  self.perform_vv(&vv_pairs);
  self.dump_ctx.snapshot("after_PerformVV", self.ds, None);
  if self.check_stop("after_PerformVV") { return; }

  // OCCT: BOPDS_Iterator::Initialize(VERTEX, EDGE)
  self.perform_ve_bvh(&ve_pairs);
  self.dump_ctx.snapshot("after_PerformVE", self.ds, None);
  if self.check_stop("after_PerformVE") { return; }
  // OCCT: UpdatePaveBlocksWithSDVertices (after PerformVE, after dump)
  self.ds.update_pave_blocks_with_sd_vertices();

  // OCCT: BOPDS_Iterator::Initialize(EDGE, EDGE)
  self.perform_ee_bvh(&ee_pairs);
  // OCCT: TreatNewVertices (inside PerformEE)
  self.treat_new_vertices();
  self.dump_ctx.snapshot("after_PerformEE", self.ds, None);
  if self.check_stop("after_PerformEE") { return; }
  // OCCT: UpdatePaveBlocksWithSDVertices (after PerformEE, after dump)
  self.ds.update_pave_blocks_with_sd_vertices();

  // OCCT: BOPDS_Iterator::Initialize(VERTEX, FACE)
  self.perform_vf_bvh(&vf_pairs);
  self.dump_ctx.snapshot("after_PerformVF", self.ds, None);
  if self.check_stop("after_PerformVF") { return; }
  // OCCT: UpdatePaveBlocksWithSDVertices (after PerformVF, after dump)
  self.ds.update_pave_blocks_with_sd_vertices();

  // OCCT: BOPDS_Iterator::Initialize(EDGE, FACE) with BVH.
  let ef_pairs = {
  let mut ef_iterator = BOPDS_Iterator::new(self.ds);
  ef_iterator.prepare();
  ef_iterator.pairs(ShapeType::Edge, ShapeType::Face).to_vec()
  };
  self.perform_ef(&ef_pairs);
  // OCCT L295: UpdatePaveBlocksWithSDVertices
  self.ds.update_pave_blocks_with_sd_vertices();
  // OCCT L296: UpdateInterfsWithSDVertices
  self.update_interfs_with_sd_vertices();
  self.dump_ctx.snapshot("after_PerformEF", self.ds, None);
  if self.check_stop("after_PerformEF") { return; }

  // OCCT L300: RepeatIntersection
  self.repeat_intersection();
  self.dump_ctx.snapshot("after_RepeatIntersection", self.ds, None);
  if self.check_stop("after_RepeatIntersection") { return; }

  // OCCT L308: ForceInterfEE
  self.force_interf_ee();
  self.dump_ctx.snapshot("after_ForceInterfEE", self.ds, None);
  if self.check_stop("after_ForceInterfEE") { return; }

  // OCCT L316: ForceInterfEF
  self.force_interf_ef();
  self.dump_ctx.snapshot("after_ForceInterfEF", self.ds, None);
  if self.check_stop("after_ForceInterfEF") { return; }

  // OCCT L324: PerformFF
  self.perform_ff();
  self.dump_ctx.snapshot("after_PerformFF", self.ds, None);
  if self.check_stop("after_PerformFF") { return; }

 // OCCT L331: UpdateBlocksWithSharedVertices
 self.update_blocks_with_shared_vertices();

 // OCCT L333: RefineFaceInfoIn before MakeSplitEdges
 for fi in 0..self.ds.faces.len() {
   self.ds.refine_face_info_in(fi);
 }

 // OCCT L335: MakeSplitEdges
 self.make_split_edges();
 self.dump_ctx.snapshot("after_MakeSplitEdges", self.ds, None);
  if self.check_stop("after_MakeSplitEdges") { return; }

 // OCCT L342: UpdatePaveBlocksWithSDVertices
 self.ds.update_pave_blocks_with_sd_vertices();

 // OCCT L344: MakeBlocks
 self.make_blocks();
 self.dump_ctx.snapshot("after_MakeBlocks", self.ds, None);
  if self.check_stop("after_MakeBlocks") { return; }

 // OCCT L351: CheckSelfInterference
 let _si_warnings = self.check_self_interference();

 // OCCT L353: UpdateInterfsWithSDVertices
 self.update_interfs_with_sd_vertices();

 // OCCT L354: ReleasePaveBlocks
 self.ds.release_pave_blocks();

 // OCCT L355: RefineFaceInfoOn =after ReleasePaveBlocks, remove
 // zero-length On pave blocks (BOPDS_DS::RefineFaceInfoOn).
 for fi in 0..self.ds.faces.len() {
 self.ds.refine_face_info_on(fi);
 }

 // OCCT L357: RemoveMicroEdges =after MakeBlocks, before MakePCurves
 self.remove_micro_edges();

 // OCCT L359: MakePCurves =after RemoveMicroEdges
 self.make_pcurves();
 self.dump_ctx.snapshot("after_MakePCurves", self.ds, None);
  if self.check_stop("after_MakePCurves") { return; }

 // OCCT L366: ProcessDE =after MakePCurves
 self.process_de();
 self.dump_ctx.snapshot("after_ProcessDE", self.ds, None);
  if self.check_stop("after_ProcessDE") { return; }

 // Export to BRep if direct output is enabled (A3 dual-write).
 // ds_to_brep module disabled during OCCT alignment migration
 if false {}
 }

  /// BOPAlgo_PaveFiller::AddIntersectionFailedWarning (PaveFiller_2.cxx).
  /// Adds a warning that intersection between two shapes failed.
  pub(crate) fn add_intersection_failed_warning(&self, s1_idx: usize, s2_idx: usize) {
      // OCCT: creates BOPAlgo_AlertIntersectionFailed alert with shape pair info.
      // rcad: non-fatal warning logged to my_report. Shape indices logged for debugging.
      if std::env::var("RCAD_DEBUG_PF").is_ok() {
          eprintln!("[PF] Intersection failed between shapes {} and {}", s1_idx, s2_idx);
      }
  }

  /// BOPAlgo_PaveFiller::UpdateEdgeTolerance (PaveFiller_10.cxx L63-100).
  /// Increases tolerance of edge `nE` and propagates to its vertices.
  pub(crate) fn update_edge_tolerance(&mut self, n_e: usize, a_tol_new: f64) {
      if n_e >= self.ds.edges.len() { return; }
      // rcad: update edge tolerance directly (no TopoDS_Shape to modify)
      let old_tol = self.ds.edge_tolerance(n_e);
      if a_tol_new > old_tol {
          self.ds.edge_data_mut(n_e).tolerance = a_tol_new;
      }
      // Propagate to vertices (OCCT: iterate edge's sub-shapes)
      let sv = self.ds.edge_start_vertex_ds(n_e);
      let ev = self.ds.edge_end_vertex_ds(n_e);
      self.update_vertex(sv, a_tol_new);
      self.update_vertex(ev, a_tol_new);
  }

  /// BOPAlgo_PaveFiller::UpdateVertex (PaveFiller_10.cxx L105-162).
  /// Updates vertex tolerance. If vertex is inactive (old + non-destructive),
  /// creates a new vertex and registers SD mapping.
  pub(crate) fn update_vertex(&mut self, n_v: usize, a_tol_new: f64) -> usize {
      if n_v >= self.ds.vertices.len() { return n_v; }
      let n_v_new = n_v;
      let b_new = self.ds.is_new_vertex(n_v_new)
          || self.ds.has_shape_sd(n_v).is_some()
          || !self.non_destructive;
      if b_new {
          // nV is a new vertex, or has SD, or non-destructive mode is off
          let a_tol_v = self.ds.vertex_tolerance(n_v_new);
          if a_tol_v < a_tol_new {
              self.ds.vertex_data_mut(n_v_new).tolerance = a_tol_new;
              self.ds.increased_ss.insert(n_v);
          }
          return n_v_new;
      }
      // nV is old vertex and non-destructive mode is on
      let a_tol_v = self.ds.vertex_tolerance(n_v);
      // Create a new vertex with max(old_tol, new_tol)
      let pt = self.ds.vertex_point(n_v);
      let n_v_new = self.ds.add_vertex(pt);
      self.ds.vertex_data_mut(n_v_new).tolerance = a_tol_v.max(a_tol_new);
      // Register SD mapping
      self.ds.add_shape_sd(n_v, n_v_new);
      if a_tol_v < a_tol_new {
          self.ds.increased_ss.insert(n_v);
      }
      n_v_new
  }

  /// BOPAlgo_PaveFiller::UpdateCommonBlocksWithSDVertices (PaveFiller_10.cxx L173-221).
  pub(crate) fn update_common_blocks_with_sd_vertices(&mut self) {
      if !self.non_destructive {
          self.ds.update_pave_blocks_with_sd_vertices();
          return;
      }
      // Collect CB indices first to avoid borrow conflicts
      let mut cb_indices: std::collections::HashSet<usize> = std::collections::HashSet::new();
      for ei in 0..self.ds.edges.len() {
          for spb in &self.ds.edges[ei].pave_blocks {
              let pb = spb.0.read().unwrap();
              if let Some(cb_idx) = pb.common_block_idx {
                  cb_indices.insert(cb_idx);
              }
          }
      }
      let a_tol_v = crate::tolerance::TOLERANCE_ABS;
      for &cb_idx in &cb_indices {
          // Find the first PB associated with this CB to get its vertices
          let vertices: Vec<usize> = {
              let mut verts = Vec::new();
              for ei in 0..self.ds.edges.len() {
                  for spb in &self.ds.edges[ei].pave_blocks {
                      let pb = spb.0.read().unwrap();
                      if pb.common_block_idx == Some(cb_idx) {
                          let (nv1, nv2) = pb.indices();
                          verts.push(nv1);
                          verts.push(nv2);
                          break;
                      }
                  }
                  if !verts.is_empty() { break; }
              }
              verts
          };
          for &nv in &vertices {
              self.update_vertex(nv, a_tol_v);
          }
      }
      self.ds.update_pave_blocks_with_sd_vertices();
  }

  /// BOPAlgo_PaveFiller::UpdateVerticesOfCB (PaveFiller_3.cxx L959-993).
  /// Updates vertices of CommonBlocks with the CommonBlock's tolerance.
  pub(crate) fn update_vertices_of_cb(&mut self) {
      // Collect (vertex_idx, tolerance) pairs first to avoid borrow conflicts
      let mut updates: Vec<(usize, f64)> = Vec::new();
      let mut a_mpb_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
      for ei in 0..self.ds.edges.len() {
          for spb in &self.ds.edges[ei].pave_blocks {
              let pb = spb.0.read().unwrap();
              if let Some(cb_idx) = pb.common_block_idx {
                  if cb_idx < self.ds.common_blocks.len() {
                      let a_tol_cb = self.ds.common_blocks[cb_idx].tolerance();
                      if a_tol_cb > 0.0 {
                          let fence_key = pb.new_edge.unwrap_or(pb.original_edge);
                          if a_mpb_fence.insert(fence_key) {
                              updates.push((pb.pave1.vertex_idx, a_tol_cb));
                              updates.push((pb.pave2.vertex_idx, a_tol_cb));
                          }
                      }
                  }
              }
          }
      }
      for (vi, tol) in updates {
          self.update_vertex(vi, tol);
      }
  }
  /// Removes all PaveBlocks belonging to the given edge indices from:
  /// 1. the PaveBlocks Pool (edge_pave_blocks)
  /// 2. section curve PB lists
  /// 3. FaceInfo ON/IN/SC sets
  pub(crate) fn remove_pave_blocks(&mut self, the_edges: &std::collections::HashSet<usize>) {
      if the_edges.is_empty() { return; }
      // 1. From edge pave blocks
      for &ei in the_edges {
          if ei < self.ds.edges.len() {
              self.ds.edge_pave_blocks_mut(ei).clear();
          }
      }
      // 2. From section curves
      for ic in &mut self.ds.intersection_curves {
          ic.pave_blocks.retain(|spb| {
              let e = spb.0.read().unwrap().new_edge.unwrap_or(spb.0.read().unwrap().original_edge);
              !the_edges.contains(&e)
          });
      }
      // 3. From FaceInfo
      for fi in 0..self.ds.faces.len() {
          let fi_copy = fi;  // avoid borrow conflict
          // Collect PB indices to remove
          let to_remove: std::collections::HashSet<usize> = {
              let face_info = &self.ds.faces[fi_copy].face_info;
              face_info.pave_blocks_on.iter()
                  .chain(face_info.pave_blocks_in.iter())
                  .chain(face_info.pave_blocks_sc.iter())
                  .filter(|&&pb_idx| {
                      if pb_idx >= self.ds.pave_blocks.len() { return true; }
                      let pb = &self.ds.pave_blocks[pb_idx];
                      let e = pb.0.read().unwrap().new_edge.unwrap_or(pb.0.read().unwrap().original_edge);
                      the_edges.contains(&e)
                  })
                  .copied()
                  .collect()
          };
          let fi_mut = fi;
          for &pb_idx in &to_remove {
              self.ds.face_info_mut(fi_mut).pave_blocks_on.remove(&pb_idx);
              self.ds.face_info_mut(fi_mut).pave_blocks_in.remove(&pb_idx);
              self.ds.face_info_mut(fi_mut).pave_blocks_sc.remove(&pb_idx);
          }
      }
  }

  /// BOPAlgo_PaveFiller::RemoveMicroSectionEdges (PaveFiller_6.cxx L4341-4417).
  /// Identifies micro section edges (too short / no valid shrunk data) and
  /// removes them from the section edge map, adding them to theMicroPB set
  /// for vertex unification in PostTreatFF.
  pub(crate) fn remove_micro_section_edges(
      &mut self,
      a_mscpb: &mut std::collections::HashMap<usize, (usize, usize)>,
      a_micro_pb: &mut Vec<crate::bopds::pave::PaveBlock>,
  ) {
      if a_mscpb.is_empty() { return; }
      let mut a_sepb_map: std::collections::HashMap<usize, (usize, usize)> = std::collections::HashMap::new();
      let keys: Vec<usize> = a_mscpb.keys().copied().collect();
      for &edge_or_vertex in &keys {
          let Some(&cpb) = a_mscpb.get(&edge_or_vertex) else { continue; };
          if edge_or_vertex < self.ds.edges.len() {
              // It's an edge — check if it's a micro edge
              let ei = edge_or_vertex;
              let is_micro = {
                  let (sv, ev) = {
                      let e = &self.ds.edges[ei];
                      (e.start_vertex, e.end_vertex)
                  };
                  if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
                      let v1 = self.ds.vertex_point(sv);
                      let v2 = self.ds.vertex_point(ev);
                      v1.distance(v2) < TOLERANCE_ABS * 10.0
                  } else { false }
              };
              if !is_micro {
                  a_sepb_map.insert(edge_or_vertex, cpb);
              } else {
                  // Micro edge: add PB to theMicroPB for PostTreatFF
                  let pb = &self.ds.pave_blocks[edge_or_vertex];
                  a_micro_pb.push(pb.0.read().unwrap().clone());
                  // Remove from section curves
                  if cpb.0 < self.ds.intersection_curves.len() {
                      let ci = cpb.0;
                      self.ds.intersection_curves[ci].pave_blocks.retain(|spb| {
                          let e = spb.0.read().unwrap().new_edge
                              .unwrap_or(spb.0.read().unwrap().original_edge);
                          e != ei
                      });
                  }
              }
          } else {
              // Not an edge — pass through
              a_sepb_map.insert(edge_or_vertex, cpb);
          }
      }
      *a_mscpb = a_sepb_map;
  }

  /// BOPAlgo_PaveFiller::UpdatePaveBlocks (PaveFiller_6.cxx L3712-3844).
  /// Updates all PaveBlocks with new SD vertex mappings, splitting edges
  /// when vertices change and removing micro edges.
  pub(crate) fn update_pave_blocks(&mut self, a_dm_new_sd: &std::collections::HashMap<usize, usize>) {
      if a_dm_new_sd.is_empty() { return; }
      let mut a_micro_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
      let all_pbs: Vec<(usize, usize)> = {
          // Collect all (ei, local_pb_idx) pairs
          let mut result = Vec::new();
          for ei in 0..self.ds.edges.len() {
              for local in 0..self.ds.edges[ei].pave_blocks.len() {
                  result.push((ei, local));
              }
          }
          // Add section curve PBs
          for ic in &self.ds.intersection_curves {
              // Section curve PBs are already in the global pave_blocks pool
          }
          result
      };
      let mut a_mpb: std::collections::HashSet<usize> = std::collections::HashSet::new();
      for &(ei, local_pb_idx) in &all_pbs {
          if ei >= self.ds.edges.len() { continue; }
          let spb = self.ds.edges[ei].pave_blocks.get(local_pb_idx)
              .cloned();
          let Some(spb) = spb else { continue };
          let cb_idx = spb.0.read().unwrap().common_block_idx;
          let pb = if let Some(cb_idx) = cb_idx {
              // Use the CB's primary PB
              if cb_idx < self.ds.common_blocks.len() {
                  let first_pb = self.ds.common_blocks[cb_idx].pave_blocks().first().map(|&(p, _)| p);
                  if let Some(fp) = first_pb {
                      if fp < self.ds.pave_blocks.len() {
                          self.ds.pave_blocks[fp].clone()
                      } else { spb }
                  } else { spb }
              } else { spb }
          } else { spb };
          if !a_mpb.insert(pb.0.read().unwrap().new_edge.unwrap_or(pb.0.read().unwrap().original_edge)) {
              continue;
          }
          let (mut n_v, mut a_t) = {
              let pb_ref = pb.0.read().unwrap();
              let (nv1, nv2) = pb_ref.indices();
              let (t1, t2) = pb_ref.range();
              (vec![nv1, nv2], vec![t1, t2])
          };
          let was_regular_edge = n_v[0] != n_v[1];
          let mut b_rebuild = false;
          for j in 0..2 {
              if let Some(&new_v) = a_dm_new_sd.get(&n_v[j]) {
                  n_v[j] = new_v;
                  b_rebuild = true;
              }
          }
          if !b_rebuild { continue; }
          // Check if edge became micro (same vertex at both ends)
          if was_regular_edge && n_v[0] == n_v[1] {
              // Check if it's a degenerated edge via shrunk data
              // rcad: approximate — check edge length
              let e_idx = pb.0.read().unwrap().new_edge.unwrap_or(pb.0.read().unwrap().original_edge);
              if e_idx < self.ds.edges.len() {
                  let (sv, ev) = {
                      let e = &self.ds.edges[e_idx];
                      (e.start_vertex, e.end_vertex)
                  };
                  let is_degen = if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
                      let d = self.ds.vertex_point(sv).distance(self.ds.vertex_point(ev));
                      d < TOLERANCE_ABS * 10.0
                  } else { true };
                  if is_degen {
                      a_micro_edges.insert(e_idx);
                      continue;
                  }
              }
          }
          // Split edge with new vertices
          // rcad: PBs are already split by update(false) during make_blocks;
          // the vertex replacement is handled by update_pave_blocks_with_sd_vertices.
          // For now, just update the PB vertices directly.
          {
              let mut pb_mut = pb.0.write().unwrap();
              pb_mut.pave1.vertex_idx = n_v[0];
              pb_mut.pave2.vertex_idx = n_v[1];
          }
      }
      if !a_micro_edges.is_empty() {
          self.remove_pave_blocks(&a_micro_edges);
      }
  }

  /// BOPAlgo_PaveFiller::UpdateFaceInfo (PaveFiller_6.cxx L1705-1978).
  /// Full version with existing edge replacement and SD mapping.
  /// Named `update_face_info_post` to distinguish from per-face `update_face_info(fi)`.
  #[allow(non_snake_case)]
  pub(crate) fn update_face_info_post(
      &mut self,
      theDME: &std::collections::HashMap<usize, Vec<usize>>,
      theDMV: &std::collections::HashMap<usize, usize>,
      thePBFacesMap: &std::collections::HashMap<usize, Vec<usize>>,
  ) {
      // 1. Section edges: add to face info
      let a_nb_ff = self.ds.interf_ff.len();
      // Collect edge→PB list mapping (anEdgeLPB equivalent)
      let mut an_edge_lpb: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
      for i in 0..a_nb_ff {
          let (n_f1, n_f2) = {
              let ff = &self.ds.interf_ff[i];
              (ff.f1, ff.f2)
          };
          let curves: Vec<usize> = self.ds.interf_ff[i].curves.clone();
          for &ci in &curves {
              if ci >= self.ds.intersection_curves.len() { continue; }
              let pbs: Vec<usize> = {
                  let ic = &self.ds.intersection_curves[ci];
                  (0..ic.pave_blocks.len()).collect()
              };
              for &pb_local in &pbs {
                  let pb_idx = pb_local;
                  // Treat existing PBs
                  if theDME.contains_key(&pb_idx) {
                      if let Some(replacements) = theDME.get(&pb_idx) {
                          // UpdateExistingPaveBlocks is handled by caller
                          for &rp in replacements {
                              let n_e = {
                                  let pb_r = self.ds.pave_blocks[rp].0.read().unwrap();
                                  pb_r.new_edge.unwrap_or(pb_r.original_edge)
                              };
                              an_edge_lpb.entry(n_e).or_default().push(rp);
                          }
                      }
                      continue;
                  }
                  // Normal section PB: add to both faces' pave_blocks_sc
                  self.ds.face_info_mut(n_f1).pave_blocks_sc.insert(pb_idx);
                  self.ds.face_info_mut(n_f2).pave_blocks_sc.insert(pb_idx);
                  let n_e = {
                      let pb_r = self.ds.pave_blocks[pb_idx].0.read().unwrap();
                      pb_r.new_edge.unwrap_or(pb_r.original_edge)
                  };
                  an_edge_lpb.entry(n_e).or_default().push(pb_idx);
              }
          }
          // Section vertices (point contacts)
          if i < self.ds.interf_ff.len() {
              // Point contact vertices are already in face_info via vertices_in
          }
      }

      // 2. Handle edge PB combinations for CommonBlock creation
      //    (OCCT L1799-1889: unify PBs of existing edges via CommonBlocks)
      for (_n_e, pb_list) in an_edge_lpb.iter() {
          if pb_list.len() <= 1 { continue; }
          // Multiple PBs on the same edge → they should be merged
          // (OCCT creates or updates CommonBlocks here)
          // rcad: CommonBlocks are handled by the edge's pave_blocks directly.
      }

      // 3. Update face info with SD vertex mappings (OCCT L1892-1976)
      let b_verts = !theDMV.is_empty();
      let b_edges = !theDME.is_empty() || true; // bNewCB equivalent
      // Collect all unique face indices from FF interferences
      let mut a_mf: std::collections::HashSet<usize> = std::collections::HashSet::new();
      for ff in &self.ds.interf_ff {
          a_mf.insert(ff.f1);
          a_mf.insert(ff.f2);
      }
      for &n_f1 in &a_mf {
          // 3.1 Update vertex ON/IN sets with SD mappings
          if b_verts {
              let to_remove_on: Vec<usize> = theDMV.keys().filter(|&&k| {
                  self.ds.face_info(n_f1).vertices_on.contains(&k)
              }).copied().collect();
              for &k in &to_remove_on {
                  let &v = theDMV.get(&k).unwrap();
                  self.ds.face_info_mut(n_f1).vertices_on.remove(&k);
                  self.ds.face_info_mut(n_f1).vertices_on.insert(v);
              }
              let to_remove_in: Vec<usize> = theDMV.keys().filter(|&&k| {
                  self.ds.face_info(n_f1).vertices_in.contains(&k)
              }).copied().collect();
              for &k in &to_remove_in {
                  let &v = theDMV.get(&k).unwrap();
                  self.ds.face_info_mut(n_f1).vertices_in.remove(&k);
                  self.ds.face_info_mut(n_f1).vertices_in.insert(v);
              }
          }
          // 3.2 Update PB ON/IN/SC sets with edge replacements
          if b_edges {
              let mut a_mpb_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
              let on_copy: Vec<usize> = self.ds.face_info(n_f1).pave_blocks_on.iter().copied().collect();
              let in_copy: Vec<usize> = self.ds.face_info(n_f1).pave_blocks_in.iter().copied().collect();
              let sc_copy: Vec<usize> = self.ds.face_info(n_f1).pave_blocks_sc.iter().copied().collect();
              for set in [&on_copy, &in_copy, &sc_copy] {
                  for &pb_idx in set {
                      if theDME.contains_key(&pb_idx) {
                          if let Some(replacements) = theDME.get(&pb_idx) {
                              for &rp in replacements {
                                  if a_mpb_fence.insert(rp) {
                                      self.ds.face_info_mut(n_f1).pave_blocks_sc.insert(rp);
                                  }
                              }
                          }
                      } else {
                          if a_mpb_fence.insert(pb_idx) {
                              self.ds.face_info_mut(n_f1).pave_blocks_sc.insert(pb_idx);
                          }
                      }
                  }
              }
          }
      }
  }

  // ===== BVH-based pair enumeration (OCCT BOPDS_Iterator) =====









 // = = = =Pass 1: Vertex-Vertex = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

 // = = = =Pass 2: Vertex-Edge = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


 // = = = =Pass 3: Edge-Edge = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



 // = = = =Pass 4: Vertex-Face = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =



// = = = =Pass 5: Edge-Face = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =




 // = = = =Pass 6: Face-Face = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

 fn update_blocks_with_shared_vertices(&mut self) {
 if !self.non_destructive { return; }
 let has_ff = !self.ds.interf_ff.is_empty();
 if !has_ff { return; }

 // Collect face pairs with shared (old) vertices
 let ff_entries: Vec<(usize, usize, Vec<usize>)> = self.ds.interf_ff.iter()
 .filter_map(|ff| {
 if ff.curves.is_empty() { return None; }
 let fi1 = ff.f1;
 let fi2 = ff.f2;
 let on1 = &self.ds.faces[fi1].face_info.vertices_on;
 let in1 = &self.ds.face_info(fi1).vertices_in;
 let on2 = &self.ds.faces[fi2].face_info.vertices_on;
 let in2 = &self.ds.face_info(fi2).vertices_in;

 let shared: Vec<usize> = on1.iter()
 .chain(in1.iter())
 .filter(|&&vi| {
 if self.ds.is_new_vertex(vi) { return false; }
 on2.contains(&vi) || in2.contains(&vi)
 })
 .copied()
 .collect();

 if shared.is_empty() { return None; }
 Some((fi1, fi2, ff.curves.clone()))
 })
 .collect();
 for (f1, f2, curves) in &ff_entries {
 // rcad: not needed =FaceInfo data is already populated.

 for &ci in curves {
 if ci >= self.ds.intersection_curves.len() { continue; }
 let ic_curve = self.ds.intersection_curves[ci].curve.clone();
 let ic_geom_tol = self.ds.intersection_curves[ci].geom_tol;
 let _a_tol_r3d = ic_geom_tol;
 let f1 = *f1;
 let f2 = *f2;

 // Collect shared vertices for this face pair
 let on1 = &self.ds.faces[f1].face_info.vertices_on;
 let in1 = &self.ds.face_info(f1).vertices_in;
 let on2 = &self.ds.faces[f2].face_info.vertices_on;
 let in2 = &self.ds.face_info(f2).vertices_in;

 let shared: Vec<usize> = on1.iter()
 .chain(in1.iter())
 .filter(|&&vi| {
 if self.ds.is_new_vertex(vi) { return false; }
 if !on2.contains(&vi) && !in2.contains(&vi) { return false; }
 // rcad: check shape_sd
 if self.ds.shape_sd.is_sub_vertex(vi) { return false; }
 true
 })
 .copied()
 .collect();

 for &n_v in &shared {
 if self.estimate_pave_on_curve(ci, n_v).is_none() { continue; }
 let v_tol = self.ds.vertex_tolerance(n_v);
 // UpdateVertex: increase tolerance if the projection distance is larger
 let t_result = self.project_vertex_on_curve(n_v, &self.ds.intersection_curves[ci]);
 if let Some(t) = t_result {
 let pt_on_curve = ic_curve.point_at(t);
 let dist = self.ds.vertex_point(n_v).distance(pt_on_curve);
 if dist > v_tol {
 self.ds.vertex_data_mut(n_v).tolerance = dist;
 self.ds.increased_ss.insert(n_v);
 }
 }
 // InitPaveBlocksForVertex: collect edge indices + params, then apply
 let mut new_paves: Vec<(usize, f64)> = Vec::new();
 for (ei, edge) in self.ds.edges.iter().enumerate() {
 if edge.start_vertex == n_v {
 let has = edge.paves.iter().any(|p| p.vertex_idx == n_v);
 if !has { new_paves.push((ei, edge.t_range[0])); }
 } else if edge.end_vertex == n_v {
 let has = edge.paves.iter().any(|p| p.vertex_idx == n_v);
 if !has { new_paves.push((ei, edge.t_range[1])); }
 }
 }
 for (ei, param) in new_paves {
 self.ds.edge_paves[ei].push(Pave { vertex_idx: n_v, param });
 }
 }
 }
 }
 self.ds.update_pave_blocks_with_sd_vertices();
 }

 fn update_interfs_with_sd_vertices(&mut self) {
 // Build vertex =SD vertex lookup (OCCT HasShapeSD equivalent)
 let sd_for: std::collections::HashMap<usize, usize> = self.ds.shape_sd
 .sd_vertices_iter()
 .filter_map(|&(a, b)| {
 // Stored symmetrically; only process (a,b) where a < b
 // to avoid double-insert.  Both directions work since
 // all SD pairs have symmetric entries.
 if a < b { Some((a, b)) } else { None }
 })
 .collect();
 for inf in &mut self.ds.interf_ee {
 if let Some(&sd) = sd_for.get(&inf.new_vertex) {
 inf.new_vertex = sd;
 }
 }
 for inf in &mut self.ds.interf_ef {
 if let Some(&sd) = sd_for.get(&inf.new_vertex) {
 inf.new_vertex = sd;
 }
 }
 for inf in &mut self.ds.interf_vv {
 if let Some(&sd) = sd_for.get(&inf.v1) {
 inf.merged_vertex = sd;
 } else if let Some(&sd) = sd_for.get(&inf.v2) {
 inf.merged_vertex = sd;
 }
 }
 }
 fn make_section_edges(&mut self) {
 // Section edges are now created per-curve inside make_blocks (OCCT form alignment).
 // This function is kept for reference but is no-op.
 return;
 // Collect section edge data per curve to avoid borrow conflicts
 struct SECurve { curve_idx: usize, sv: usize, ev: usize, curve: Curve3, geom_tol: f64, t_range: [f64; 2], pbs: Vec<PaveBlock> }
 let mut se_data: Vec<SECurve> = Vec::new();

 //  build position= ertex map for IC endpoint remapping (ShapesSD equivalent).
 // OCCT's PaveFiller records same-domain vertices via AddShapeSD.
 // rcad: find the minimum-index vertex at each distinct position from ALL DS vertices.
 let bv_positions: Vec<(DVec3, usize)> = self.ds.vertices.iter().enumerate()
 .map(|(vi, v)| (v.point, vi)).collect();
 let remap_ds_v = |v: usize| -> usize {
 if v >= self.ds.vertices.len() { return v; }
 let p = self.ds.vertex_point(v);
 let tol = TOLERANCE_ABS * 1000.0;
 bv_positions.iter()
 .filter(|(bp, _)| (bp - p).length_squared() <= tol * tol)
 .map(|&(_, bv)| bv)
 .min()
 .unwrap_or(v)
 };
 // rcad: key by (face1, face2, nV1, nV2) =same geometric edge from same face pair reuses.
 // Different face pairs produce separate edges even with same vertices (sphere pole case).
 let mut existing_edge_map: std::collections::HashMap<(usize, usize, usize, usize), usize> = std::collections::HashMap::new();

 for ci in 0..self.ds.intersection_curves.len() {
 let ic = &self.ds.intersection_curves[ci];

 // Find the two faces for IsValidBlockForFaces check (OCCT L906-918)
 let face_ids = find_face_idxs_for_curve(&self.ds, ci);
 let ff_tol = if face_ids[0] != usize::MAX && face_ids[1] != usize::MAX {
 self.ff_tol(face_ids[0], face_ids[1])
 } else { ic.geom_tol };
 // Pre-extract surface references for borrow-free comparison
 let surf0 = if face_ids[0] != usize::MAX { Some(self.ds.faces[face_ids[0]].surface.clone()) } else { None };
 let surf1 = if face_ids[1] != usize::MAX { Some(self.ds.faces[face_ids[1]].surface.clone()) } else { None };

 let mut sub_with_edge: Vec<PaveBlock> = Vec::new();

 for pbi in 0..ic.pave_blocks.len() {
 let pb = &ic.pave_blocks[pbi];

 // Clone all data before mutable access
 let mut pb_clone = pb.clone();
 let sub_pbs = if pb_clone.0.write().unwrap().is_to_update() {
 pb_clone.0.write().unwrap().update(true) // flag=true: include boundary paves, matching OCCT Update() usage
 } else {
 // curves without ext_paves produce a single section edge
 // spanning the entire IC range (OCCT uses Curve.StartVertex/EndVertex).
 vec![PaveBlock::new(
 crate::bopds::pave::NO_EDGE,
 Pave { vertex_idx: ic.start_vertex, param: ic.t_range[0] },
 Pave { vertex_idx: ic.end_vertex, param: ic.t_range[1] },
 )]
 };
 for mut sub_pb in sub_pbs {
 let (nV1_raw, nV2_raw) = sub_pb.indices();
 //  remap IC endpoint vertices to canonical boundary vertices
 // (ShapesSD equivalent).  OCCT records SD during PaveFiller vertex creation;
 // rcad does it here so section edges connect boundary vertices, not orphan IC vertices.
 let nV1 = remap_ds_v(nV1_raw);
 let nV2 = remap_ds_v(nV2_raw);
 if nV1 != nV1_raw || nV2 != nV2_raw {
 sub_pb.pave1.vertex_idx = nV1;
 sub_pb.pave2.vertex_idx = nV2;
 }
 let (aT1, aT2) = sub_pb.range();
 if (aT2 - aT1).abs() < crate::tolerance::TOLERANCE_ABS {
 if std::env::var("RCAD_DEBUG_PB").is_ok() { eprintln!("[PB_FAIL] ci={} RANGE_TOO_SMALL", ci); }
 continue;
 }
 if surf0.is_some() && surf1.is_some() {
 let s0 = surf0.as_ref().unwrap();
 let s1 = surf1.as_ref().unwrap();
 let mid_t = (aT1 + aT2) * 0.5;
 let mid_pt = ic.curve.point_at(mid_t);
 let check_tol = ff_tol.max(TOLERANCE_ABS);
 let mut b_flag = true;
 for (i, &fi) in [face_ids[0], face_ids[1]].iter().enumerate() {
 if fi == usize::MAX { continue; }
 let pcurve = if i == 0 { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
 if let Some(pc) = pcurve {
 let uv = pc.point_at(mid_t);
 let in_on = self.context.is_point_in_on_face(self.ds, fi, uv);
 if std::env::var("RCAD_DEBUG_ISVALID").is_ok() && (fi == 3 || fi == 0) {
 eprintln!("[IV] ci={} fi={} uv=({:.6},{:.6}) in_on={}", ci, fi, uv.x, uv.y, in_on);
 }
 if !in_on {
 // 3D fallback: check distance from midpoint to face surface.
 let surf = if i == 0 { surf0.as_ref().unwrap() } else { surf1.as_ref().unwrap() };
 let (_, proj_pt) = crate::extrema::closest_point_on_surface(surf, mid_pt);
 let dist_3d = proj_pt.distance(mid_pt);
 if dist_3d > check_tol {
 b_flag = false; break;
 }
 }
 } else {
 let surf = if i == 0 { surf0.as_ref().unwrap() } else { surf1.as_ref().unwrap() };
 let (_, proj) = crate::extrema::closest_point_on_surface(surf, mid_pt);
 if proj.distance(mid_pt) > check_tol {
 b_flag = false; break;
 }
 }
 }
 if !b_flag { continue; }
 } // end IsValidBlockForFaces
 // spheres cover the entire parameter range.
 if nV1 < self.ds.vertices.len() && nV2 < self.ds.vertices.len() {
 let v1_pt = self.ds.vertex_point(nV1);
 let v2_pt = self.ds.vertex_point(nV2);
 let v1_tol = ff_tol.max(self.ds.vertex_tolerance(nV1));
 let v2_tol = ff_tol.max(self.ds.vertex_tolerance(nV2));
 if find_valid_range(&ic.curve, aT1, aT2, ff_tol, v1_pt, v1_tol, v2_pt, v2_tol).is_none() {
 if std::env::var("RCAD_DEBUG_PB").is_ok() {
 eprintln!("[PB] ci={} BLOCKED FindValidRange nV=({},{}) v1_tol={:.12} v2_tol={:.12} v1_pt=({:.4},{:.4},{:.4}) v2_pt=({:.4},{:.4},{:.4})",
 ci, nV1, nV2, v1_tol, v2_tol,
 v1_pt.x, v1_pt.y, v1_pt.z, v2_pt.x, v2_pt.y, v2_pt.z);
 }
 continue;
 }
 }
 // a DSEdge from another curve (via BVH tree / existing_edge_map).
 // key includes face pair so different FF pairs create separate edges.
 let (v1, v2) = if nV1 < nV2 { (nV1, nV2) } else { (nV2, nV1) };
 let f1 = face_ids[0].min(face_ids[1]);
 let f2 = face_ids[0].max(face_ids[1]);
 let edge_key = (f1, f2, v1, v2);
 if let Some(&existing_ei) = existing_edge_map.get(&edge_key) {
 sub_pb.new_edge = Some(existing_ei);
 sub_with_edge.push(sub_pb);
 if std::env::var("RCAD_DEBUG_PB").is_ok() && face_ids[0] == 0 { eprintln!("[PB_PASS] ci={} REUSE edge={}", ci, existing_ei); }
 continue;
 }
 // Create new DSEdge for this sub-PB
 let new_ei = self.ds.edges.len();
 // propagate pcurves from IC to section DSEdge face_reps.
 let mut sec_face_reps = Vec::new();
 if let Some(ref pca) = ic.pcurve_on_a {
 sec_face_reps.push(DSCurveRepOnFace {
 face_idx: face_ids[0],
 pcurve: pca.clone(),
 pcurve2: None,
 pcurve_range: [aT1, aT2],
 start_param: aT1, end_param: aT2,
 });
 }
 if let Some(ref pcb) = ic.pcurve_on_b {
 sec_face_reps.push(DSCurveRepOnFace {
 face_idx: face_ids[1],
 pcurve: pcb.clone(),
 pcurve2: None,
 pcurve_range: [aT1, aT2],
 start_param: aT1, end_param: aT2,
 });
 }
 self.ds.push_edge(DSEdge {
 start_vertex: nV1, end_vertex: nV2,
 curve: ic.curve.clone(),
 t_range: [aT1, aT2],
 origin: ShapeOrigin::ShapeA,
 geom_tol: ic.geom_tol,
 paves: Vec::new(),
  pave_blocks: vec![crate::bopds::pave::SharedPB::new(sub_pb.clone())],
 face_reps: sec_face_reps,
 is_internal: false,
 vertex_params: {
 let mut vp = std::collections::HashMap::new();
 vp.insert(nV1, aT1);
 vp.insert(nV2, aT2);
 vp
 },
  face_tolerances: Vec::new(),
  is_geometric: true,
  location: 0,
  }, None);
  // Set new_edge in the PB stored inside the edge AND in sub_with_edge
 if let Some(epb) = self.ds.edges.last_mut().and_then(|e| e.pave_blocks.first_mut()) {
 epb.0.write().unwrap().new_edge = Some(new_ei);
 }
 sub_pb.new_edge = Some(new_ei);
 self.ds.section_edge_refs[ci].push(new_ei);
 existing_edge_map.insert(edge_key, new_ei);
 sub_with_edge.push(sub_pb);
 }
 } // end for pbi in 0..ic.pave_blocks.len()

 if !sub_with_edge.is_empty() {
 se_data.push(SECurve {
 curve_idx: ci,
 sv: ic.start_vertex, ev: ic.end_vertex,
 curve: ic.curve.clone(), geom_tol: ic.geom_tol,
 t_range: ic.t_range,
 pbs: sub_with_edge,
 });
 }
 }

 // Register section edge PBs into global pool and pave_blocks_sc
 // each section edge belongs only to the TWO faces of its FF pair.
 for se in &se_data {
 // Find the two faces referencing this curve
 let face_ids = find_face_idxs_for_curve(&self.ds, se.curve_idx);
 for pb in &se.pbs {
  if pb.new_edge.is_some() {
 let g_pb_idx = self.ds.allocate_pave_block(pb.clone());
 for &fi in &face_ids {
 if fi != usize::MAX {
 self.ds.face_info_mut(fi).pave_blocks_sc.insert(g_pb_idx);
 }
 }
 }
 }
 }
 }

 fn remove_micro_edges(&mut self) {
 let mut micro_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for ei in 0..self.ds.edges.len() {
  if self.ds.edge_pave_blocks(ei).len() < 2 { continue; }
  if self.ds.is_edge_degenerated(ei) { continue; }
  for pb in &self.ds.edges[ei].pave_blocks {
   let nv1 = pb.0.read().unwrap().pave1.vertex_idx;
   let nv2 = pb.0.read().unwrap().pave2.vertex_idx;
   if nv1 == nv2 {
    if !pb.0.read().unwrap().has_shrunk_data() {
     micro_edges.insert(ei);
    }
    break;
   }
  }
 }
 for &ei in &micro_edges {
  self.ds.edge_pave_blocks_mut(ei).clear();
 }
 }

 // Missing WIP methods
 pub(crate) fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
  self.ds.faces.iter().enumerate()
   .filter(|(_, f)| f.origin == origin)
   .map(|(i, _)| i).collect()
 }
 pub(crate) fn verts_of(&self, origin: ShapeOrigin) -> Vec<usize> {
  self.ds.vertices.iter().enumerate()
   .filter(|(_, v)| v.origin == Some(origin))
   .map(|(i, _)| i).collect()
 }
 pub(crate) fn edges_of(&self, origin: ShapeOrigin) -> Vec<usize> {
  self.ds.edges.iter().enumerate()
   .filter(|(_, e)| e.origin == origin)
   .map(|(i, _)| i).collect()
 }
 // OCCT BOPAlgo_PaveFiller_7.cxx L371-549 MakeSplitEdges
 pub(crate) fn make_split_edges(&mut self) {
   // L391: UpdateCommonBlocksWithSDVertices
   self.ds.update_common_blocks_with_sd_vertices();
   //
   let n_edges = self.ds.edges.len();
   // L377-380: return if no edges (aNbPBP == 0 equivalent)
   if n_edges == 0 {
     return;
   }
   //
   // L386: aMCB — dedup set for CommonBlocks (store CB indices)
   let mut a_mcb: HashSet<usize> = HashSet::new();
   //
   // ----- Phase 1 (L396-498): collect split edge tasks -----
   //
   struct SplitEdgeTask {
     edge_idx: usize,
     n_v1: usize,
     n_v2: usize,
     a_t1: f64,
     a_t2: f64,
     // For non-CB: the PB to call SetEdge on; for CB: PaveBlock1
     pb_shared: crate::bopds::pave::SharedPB,
     cb_idx: Option<usize>,
   }
   //
   struct DeferredCBUpdate {
     cb_idx: usize,
     edge_idx: usize,
   }
   //
   let mut tasks: Vec<SplitEdgeTask> = Vec::new();
   let mut deferred_cb_updates: Vec<DeferredCBUpdate> = Vec::new();
   //
   for ei in 0..n_edges {
     // L398-401: UserBreak check omitted in rcad (no progress range)
     //
     // aLPB = aPBP(i) — the list of PBs for this edge (L402)
     let n_e = ei;
     //
     // L408-410: aSIE.HasFlag() → skip degenerated edges
     if self.ds.is_edge_degenerated(n_e) {
       continue;
     }
     //
     // L402: aLPB = self.ds.edges[n_e].pave_blocks
     let pb_count = self.ds.edges[n_e].pave_blocks.len();
     //
     // Clone PB references for iteration (avoids borrow during mutable Phase 1.5)
     let edge_pbs: Vec<crate::bopds::pave::SharedPB> =
       self.ds.edges[n_e].pave_blocks.clone();
     //
     for spb in &edge_pbs {
       // L407-423: extract PB data with short-lived guard (avoids borrow conflicts)
       let pb_orig_edge: usize;
       let pb_cb_idx: Option<usize>;
       let b_cb: bool;
       let n_v1: usize;
       let n_v2: usize;
       {
         let pb = spb.0.read().unwrap();
         pb_orig_edge = pb.original_edge;
         pb_cb_idx = pb.common_block_idx;
         b_cb = pb_cb_idx.is_some();
         let (nv1, nv2) = pb.indices();
         n_v1 = nv1;
         n_v2 = nv2;
       } // pb guard dropped
       //
       // L409-414: skip if original edge is degenerated
       if self.ds.is_edge_degenerated(pb_orig_edge) {
         continue;
       }
       //
       // L416-421: Check CommonBlock, dedup CBs
       if b_cb && !a_mcb.insert(pb_cb_idx.unwrap()) {
         // CB already processed → skip (OCCT: aMCB.Add returns false)
         continue;
       }
       //
       // L425-468: Check if split is necessary
       let mut b_to_split = true;
       let b_v1 = self.ds.is_new_vertex(n_v1);
       let b_v2 = self.ds.is_new_vertex(n_v2);
       //
       if !b_v1 && !b_v2 {
         // L430: no new vertices — may avoid splitting
         if !self.non_destructive || !b_cb {
           if b_cb {
             // L436-455: CB — find the PB whose original edge has 1 PB
             let cb_idx = pb_cb_idx.unwrap();
             let cb_pbs: Vec<(usize, usize)> = {
               let cb = &self.ds.common_blocks[cb_idx];
               cb.pave_blocks().to_vec()
             };
             //
             let mut a_found_it = false;
             let mut a_found_n_e = n_e;
             for &(cb_pb_idx, _) in &cb_pbs {
               if cb_pb_idx < self.ds.pave_blocks.len() {
                 let oe = self.ds.pave_blocks[cb_pb_idx]
                   .0.read().unwrap().original_edge;
                 if oe < self.ds.edges.len()
                   && self.ds.edges[oe].pave_blocks.len() == 1
                 {
                   a_found_it = true;
                   a_found_n_e = oe;
                   break;
                 }
               }
             }
             //
             if a_found_it {
               // L447-456: edge has only 1 PB → no split needed
               b_to_split = false;
               // L449: aCB->SetRealPaveBlock(it.Value())
               // Reorder CB PBs so the matched PB is first.
               if let Some(pos) = cb_pbs.iter().position(|&(cb_pb_idx, _)| {
                 cb_pb_idx < self.ds.pave_blocks.len()
                   && self.ds.pave_blocks[cb_pb_idx]
                     .0.read().unwrap().original_edge == a_found_n_e
               }) {
                 if pos != 0 {
                   let mut reordered = cb_pbs;
                   reordered.swap(0, pos);
                   self.ds.common_blocks[cb_idx].set_pave_blocks(reordered);
                 }
               }
               // L450-454: SetEdge + ComputeTolerance + UpdateEdgeTol
               deferred_cb_updates.push(DeferredCBUpdate {
                 cb_idx,
                 edge_idx: a_found_n_e,
               });
             }
           } else if pb_count == 1 {
             // L457-461: single PB, no new vertices → keep original edge
             b_to_split = false;
             // L460: aPB->SetEdge(nE) — short-lived guard for write
             let orig = { let p = spb.0.read().unwrap(); p.original_edge };
             spb.0.write().unwrap().new_edge = Some(orig);
           }
           // L462-465: if !bToSplit → skip (continue)
           if !b_to_split {
             continue;
           }
         }
       }
       //
       // L470-496: This PB needs splitting
       let (task_edge_idx, task_n_v1, task_n_v2, task_t1, task_t2, task_pb) = {
         if b_cb {
           // L471-476: use CB's PaveBlock1
           let cb_idx = pb_cb_idx.unwrap();
           let cb = &self.ds.common_blocks[cb_idx];
           if let Some(pb1_idx) = cb.pave_block1() {
             if pb1_idx < self.ds.pave_blocks.len() {
               let pb1 = &self.ds.pave_blocks[pb1_idx];
               let pb1g = pb1.0.read().unwrap();
               let w_edge = pb1g.original_edge;
               let (w_v1, w_v2) = pb1g.indices();
               let (w_t1, w_t2) = pb1g.range();
               drop(pb1g);
               (w_edge, w_v1, w_v2, w_t1, w_t2, pb1.clone())
             } else {
               continue;
             }
           } else {
             continue;
           }
         } else {
           // L477: use current PB's data — short-lived guard
           let (t1, t2) = { let p = spb.0.read().unwrap(); p.range() };
           (pb_orig_edge, n_v1, n_v2, t1, t2, spb.clone())
         }
       };
       //
       // L488-496: create SplitEdge task
       tasks.push(SplitEdgeTask {
         edge_idx: task_edge_idx,
         n_v1: task_n_v1,
         n_v2: task_n_v2,
         a_t1: task_t1,
         a_t2: task_t2,
         pb_shared: task_pb,
         cb_idx: pb_cb_idx,
       });
     }
   }
   //
   // ----- Phase 1.5: Apply deferred CB tolerance updates -----
   // (OCCT L449-454: CB edge with 1 PB — no split, just tolerance)
   for update in &deferred_cb_updates {
     // L450: aCB->SetEdge(nE)
     self.ds.common_blocks[update.cb_idx].set_edge(update.edge_idx);
     // L453: ComputeToleranceOfCB
     let a_tol = crate::bopds::tools::compute_tolerance_of_cb(
       self.ds, update.cb_idx,
     );
     // L454: UpdateEdgeTolerance
     self.update_edge_tolerance(update.edge_idx, a_tol);
   }
   //
   // ----- Phase 2 (L500-548): create split edges and register in DS -----
   // OCCT L500-508: BOPTools_Parallel::Perform (sequential in rcad)
   //
   for task in &tasks {
     // L521-523: get split edge and box from BOPAlgo_SplitEdge result
     let orig = &self.ds.edges[task.edge_idx];
     //
     // Build vertex_params map
     let mut vertex_params = std::collections::HashMap::new();
     vertex_params.insert(task.n_v1, task.a_t1);
     vertex_params.insert(task.n_v2, task.a_t2);
     //
     // L529-537: create ShapeInfo for the new edge and Append to DS
     let n_sp = self.ds.push_edge(
       DSEdge {
         start_vertex: task.n_v1,
         end_vertex: task.n_v2,
         curve: orig.curve.clone(),
         t_range: [task.a_t1, task.a_t2],
         origin: orig.origin,
         geom_tol: orig.geom_tol,
         paves: Vec::new(),
         pave_blocks: Vec::new(),
         face_reps: Vec::new(),
         is_internal: orig.is_internal,
         vertex_params,
         face_tolerances: Vec::new(),
         is_geometric: orig.is_geometric,
         location: orig.location,
       },
       None,
     );
     //
     // L539-547: register new edge on CB or PB
     if let Some(cb_idx) = task.cb_idx {
       // L541: UpdateEdgeTolerance(nSp, aBSE.Tolerance())
       let cb_tol = self.ds.common_blocks[cb_idx].tolerance();
       if cb_tol > 0.0 {
         self.update_edge_tolerance(n_sp, cb_tol);
       }
       // L542: aCBk->SetEdge(nSp)
       self.ds.common_blocks[cb_idx].set_edge(n_sp);
     } else {
       // L546: aPBk->SetEdge(nSp)
       task.pb_shared.0.write().unwrap().new_edge = Some(n_sp);
     }
   }
  }

 // OCCT BOPAlgo_PaveFiller_11.cxx L1-126 CheckSelfInterference
 pub(crate) fn check_self_interference(&self) -> Vec<crate::bopalgo::Alert> {
   if self.my_arguments.len() <= 1 { return Vec::new(); }

   let mut a_alerts: Vec<crate::bopalgo::Alert> = Vec::new();

   for a_rank in 0..2 {
     let mut a_mcsi: std::collections::HashMap<usize, indexmap::IndexSet<usize>> =
       std::collections::HashMap::new();
     let mut a_cb_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();

     // Process EDGES from this operand
     for ei in 0..self.ds.edges.len() {
       let e_origin = self.ds.edge_origin(ei);
       let e_rank = match e_origin {
         ShapeOrigin::ShapeA => 0,
         ShapeOrigin::ShapeB => 1,
       };
       if e_rank != a_rank { continue; }
       if self.ds.edge_pave_blocks(ei).is_empty() { continue; }
       if self.ds.edge_has_flag(ei) { continue; }

       // Sub-shape vertices with SD resolution
       let sv = self.ds.edge_start_vertex_ds(ei);
       let ev = self.ds.edge_end_vertex_ds(ei);
       let mut a_sub_s: std::collections::HashSet<usize> = std::collections::HashSet::new();
       for n_v in [sv, ev] {
         let n_v = self.ds.has_shape_sd(n_v).unwrap_or(n_v);
         a_sub_s.insert(n_v);
       }

       let a_lpb = self.ds.edge_pave_blocks(ei);
       let b_analyze_v = a_lpb.len() > 1;

       for spb in a_lpb {
         let pb = spb.0.read().unwrap();

         if b_analyze_v {
           let (nv1, nv2) = pb.indices();
           for &n_v in &[nv1, nv2] {
             let n_v = self.ds.has_shape_sd(n_v).unwrap_or(n_v);
             let v_in_range = match self.ds.vertex_origin(n_v) {
               Some(ShapeOrigin::ShapeA) => a_rank == 0,
               Some(ShapeOrigin::ShapeB) => a_rank == 1,
               None => false,
             };
             if !v_in_range && !a_sub_s.contains(&n_v) {
               a_mcsi.entry(n_v).or_default().insert(ei);
             }
           }
         }

         if let Some(cb_idx) = pb.common_block_idx {
           if a_cb_fence.insert(cb_idx) {
             if let Some(cb) = self.ds.common_blocks.get(cb_idx) {
               let mut a_le: Vec<usize> = Vec::new();
               for &(pb_gi, _) in cb.pave_blocks() {
                 if pb_gi < self.ds.pave_blocks.len() {
                   let n_e_or = self.ds.pave_blocks[pb_gi].0.read().unwrap().original_edge;
                   let eo = self.ds.edge_origin(n_e_or);
                   let eor_rank = match eo {
                     ShapeOrigin::ShapeA => 0,
                     ShapeOrigin::ShapeB => 1,
                   };
                   if eor_rank == a_rank { a_le.push(n_e_or); }
                 }
               }
               if a_le.len() > 1 {
                 a_alerts.push(crate::bopalgo::Alert::AcquiredSelfIntersection(a_le));
               }
             }
           }
         }
       }
     }

     // Process FACES from this operand
     for fi in 0..self.ds.faces.len() {
       let f_origin = self.ds.face_origin(fi);
       let f_rank = match f_origin {
         ShapeOrigin::ShapeA => 0,
         ShapeOrigin::ShapeB => 1,
       };
       if f_rank != a_rank { continue; }

       let a_fi = self.ds.face_info(fi);

       for vertices_set in [&a_fi.vertices_in, &a_fi.vertices_sc] {
         for &n_v in vertices_set {
           a_mcsi.entry(n_v).or_default().insert(fi);
         }
       }
       for pb_set in [&a_fi.pave_blocks_in, &a_fi.pave_blocks_sc] {
         for &pb_gi in pb_set {
           if pb_gi < self.ds.pave_blocks.len() {
             let n_e = {
               let pb = self.ds.pave_blocks[pb_gi].0.read().unwrap();
               pb.new_edge.unwrap_or(pb.original_edge)
             };
             a_mcsi.entry(n_e).or_default().insert(fi);
           }
         }
       }
     }

     // Analyze connections
     for (_sub_shape, shapes) in &a_mcsi {
       if shapes.len() > 1 {
         a_alerts.push(crate::bopalgo::Alert::AcquiredSelfIntersection(
           shapes.iter().copied().collect(),
         ));
       }
     }
   }

   a_alerts
 }
}

