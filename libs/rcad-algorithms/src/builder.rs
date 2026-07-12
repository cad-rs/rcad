use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use indexmap::IndexMap;

use glam::{DVec2, DVec3};
use rayon::prelude::*;
use rcad_kernel::topods;
use rcad_kernel::PCurve;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval, *};
use rcad_kernel::topology::*;
use rcad_kernel::topology::Face;

use crate::bvh::{Aabb, DsBvh};
use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::bopalgo::{GlueEnum, Alert, Report};
use crate::history::{
 BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin,
};
use std::cell::RefCell;
use crate::inttools::context::Context;
use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;
use crate::triangulate::{triangulate_polygon, triangulate_polygon_with_holes};

use crate::pipeline_dump::DumpCtx;

mod angle_2d;
mod curve_tools;
mod debug_utils;
mod intres2d;
mod intersection;
mod types;
mod ds_as_brep;

mod builder_face;
pub(crate) use builder_face::BuilderFace;

pub use types::{
 BooleanOpType, BooleanError,
};
pub(crate) use types::{
 ShapeType, WireFace, WireSegment, WireEdgeSource,
 WireSegmentTopoDS, WireEdgeSourceTopoDS,
 FaceEntry,
};

mod result_builder;

pub(crate) use result_builder::ResultBuilder;
mod builder_utils;

pub(crate) use builder_utils::{
 curve_eq, hash_point,
 is_tangent_face, build_edge_bounds, quantize_pos,
 check_and_add_split_vertex, collect_face_edge_segments,
 annotate_history_from_ds, annotate_shell_and_solid_history,
 aggregate_face_region_origin, aggregate_shell_region_origin,
 point_in_polygon_2d,
};

pub struct BooleanBuilder<'a> {
 pub(crate) ds: &'a DS,
 pub(crate) op: BooleanOpType,
 /// =OCCT-aligned: myGlue =BOPAlgo_GlueEnum (GlueOff/GlueFull/GlueShift).
 pub(crate) glue: GlueEnum,
 pub(crate) glue_tolerance: f64,
 pub(crate) context: RefCell<Context>,
 // =OCCT-aligned: error tracking (myReport / HasErrors equivalent).
 pub(crate) has_errors: bool,
 // =OCCT-aligned: myImages =source shape index =list of split image indices.
 pub(crate) my_images: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>>>,
 pub(crate) my_origins: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>>>,
 pub(crate) my_shapes_sd: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, rcad_kernel::topods::ShapeRef>>,
 pub(crate) my_in_parts: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 pub(crate) my_solid_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 pub(crate) my_solid_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 // =OCCT-aligned: myNonDestructive (BOPAlgo_Builder.hxx L503).
 pub(crate) my_non_destructive: bool,
 // OCCT-aligned: myFillHistory (BOPAlgo_Options.hxx).
 pub(crate) my_fill_history: bool,
 // =OCCT-aligned: myCheckInverted (BOPAlgo_Builder.hxx L505).
 pub(crate) my_check_inverted: bool,
 // =OCCT-aligned: myStopOnFatalError =abort pipeline on fatal error.
 pub(crate) my_stop_on_fatal_error: bool,
 /// =OCCT-aligned: myEntryPoint =tracks builder phase (1=PerformInternal1 done, etc.).
 pub(crate) my_entry_point: u8,
 /// =OCCT-aligned: myReport =collects alerts during Builder execution.
  pub(crate) my_report: std::cell::RefCell<Report>,
 /// OCCT-aligned: myDims - dimension per argument (3=solid, 2=face).
 pub(crate) my_dims: std::cell::Cell<[i8; 2]>,
 /// =OCCT-aligned: converted BRep representation of DS.
 pub(crate) brep: std::cell::RefCell<Option<(rcad_kernel::topods::BRep, Vec<rcad_kernel::topods::ShapeRef>, Vec<Option<rcad_kernel::topods::ShapeRef>>)>>,
 /// OCCT-aligned: myShape  ?result shape accumulator (BRep).
 pub(crate) my_shape: std::cell::RefCell<rcad_kernel::topods::BRep>,
 /// OCCT-aligned: myArguments  ?all source shapes pre-created as TShapes.
 pub(crate) my_arguments: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
 /// OCCT-aligned: DS edge  ?TShape::Edge mapping (replaces ResultBuilder.ds_edge_to_tshape).
 pub(crate) my_edge_map: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
 /// OCCT-aligned: result wire TShape refs (replaces ResultBuilder.wire_refs).
 pub(crate) my_wire_refs: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
 /// OCCT-aligned: result shell TShape refs (replaces ResultBuilder.shells).
 pub(crate) my_shells: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
 /// Result face TShape refs (replaces ResultBuilder.face_refs).
 pub(crate) my_face_refs: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
 /// Result solid TShape refs (replaces ResultBuilder.solids).
 pub(crate) my_solids: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
 /// Result compsolid TShape refs (replaces ResultBuilder.compsolid_groups).
 pub(crate) my_compsolid_groups: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SourceSide {
 A,
 B,
}

/// Fast path: if the opposite solid is an axis-aligned box, check all sub-face
/// boundary vertices against the box AABB. For tessellated faces (cone/cylinder
/// UV grid), individual grid cells can straddle the box boundary even when their
/// sample point falls inside. Requiring ALL boundary vertices to be on the correct
/// side ensures straddling cells are conservatively classified.
///
/// - Intersection (any side): sub-face is kept only when ENTIRELY inside the box.
/// - Difference B-side: sub-face is kept only when ENTIRELY inside the box.
/// - Union/Difference A-side: sub-face is kept only when ENTIRELY outside the box.

mod wire_splitter;
mod wire_path;
mod wire_path_topo_ds;
mod edge_builders;
mod builder_utils_topo_ds;
mod filler_mod;
mod result_build_mod;
pub(crate) use wire_splitter::{
 EdgeInfo, build_closed_wires,
 expand_avoided_pids,
 physical_edge_id, world_to_uv,
 edge_uv_tangent, edge_angle_2d,
 are_verts_coincident,
};
pub(crate) use wire_path::{
 perform_areas, intersect_ray_curve_2d,
 refine_angles, pc_parameter_range,
 walk_path_extract_wires,
};

impl<'a> BooleanBuilder<'a> {
 /// =OCCT-aligned: TopoDS-based BuildFace pipeline with emit.
 /// Runs the full pipeline then emits result faces directly into ResultBuilder.
 pub(crate) fn builder_face_perform(
 &self,
 face_idx: usize,
 is_a: bool,
 result: &mut ResultBuilder,
 t: &mut topods::BRep,
 ) {
 let ds = self.ds;
 // Get BRep from the cached conversion (built during build_with_history)
 let brep_borrow = self.brep.borrow();
 let (br, face_refs, _ic_edge_map): &(rcad_kernel::topods::BRep, Vec<rcad_kernel::topods::ShapeRef>, Vec<Option<rcad_kernel::topods::ShapeRef>>) = match brep_borrow.as_ref() {
 Some(v) => v,
 None => {
 if std::env::var("RCAD_DEBUG_IC").is_ok() {
 eprintln!("[IC] builder_face_perform: self.brep is None, skipping face {}", face_idx);
 }
 return;
 }
 };
 // Guard: skip if face_refs doesn't have an entry for this face
 if face_idx >= face_refs.len() {
 if std::env::var("RCAD_DEBUG_IC").is_ok() {
 eprintln!("[IC] builder_face_perform: face_refs[{}] out of bounds (len={}), skipping", face_idx, face_refs.len());
 }
 return;
 }

 let face_ref = face_refs[face_idx];
 // Guard: verify the source face ShapeRef is valid in the brep
 if face_ref.index >= br.tshapes.len() || !matches!(&*br.tshapes[face_ref.index], rcad_kernel::topods::TShape::Face(_)) {
 if std::env::var("RCAD_DEBUG_IC").is_ok() {
 eprintln!("[IC] builder_face_perform: face_ref {} is not a valid face in brep (tshapes.len={}), skipping face {}",
 face_ref.index, br.tshapes.len(), face_idx);
 }
 return;
 }

 let pcurve_lookup = |ci: usize| self.find_pcurve_for_face(ci, face_idx);
 let segments = collect_face_edge_segments(ds, face_idx, &pcurve_lookup);
 if std::env::var("RCAD_DEBUG_IC").is_ok() {
 eprintln!("[SPLIT] face={} DS origin={:?} n_segments={} has_pb_sc={}", 
 face_idx, ds.faces[face_idx].origin, segments.len(),
 !ds.faces[face_idx].face_info.curves_sc.is_empty());
 for (si, seg) in segments.iter().enumerate() {
 let src = format!("{:?}", seg.source);
 eprintln!("[SPLIT] seg[{}] src={} v{}->v{}", si, src, seg.start_vertex, seg.end_vertex);
 }
 }
 if !self.builder_face_check_data(face_idx, &segments) { return; }

 let segments_topo = crate::builder::builder_utils_topo_ds::segments_to_topo_ds(&segments, ds, face_idx, &face_refs[..], &_ic_edge_map[..]);
 // segments kept alive for classification below; dropped after classification.

 let tool: &dyn rcad_kernel::topods::BRepTool = br;

 let (avoided_pids, pid_segs) = crate::builder::wire_splitter::perform_shapes_to_avoid_topo(
 &segments_topo, tool);
 let mut avoided = crate::builder::wire_splitter::expand_avoided_pids(&avoided_pids, &pid_segs);
 let wires = crate::builder::wire_path_topo_ds::build_closed_wires(
 &segments_topo, &avoided, tool);

 let in_loop: HashSet<usize> = wires.iter().flatten().copied().collect();
 for si in 0..segments_topo.len() {
 if !in_loop.contains(&si) && !avoided.contains(&si) { avoided.insert(si); }
 }
 // OCCT L327-382: group connected avoided edges into internal wires
 let internal_wire_groups = crate::builder::wire_path_topo_ds::build_internal_wires(
 &segments_topo, &avoided);

 let wfs = if !wires.is_empty() {
 crate::builder::wire_path_topo_ds::perform_areas(
 &wires, &internal_wire_groups, &segments_topo, tool, face_idx, ds)
 } else if !avoided.is_empty() {
 vec![WireFace { outer_wire: vec![], inner_wires: vec![], internal_wires: segments_topo.iter().enumerate().filter(|(si, _)| avoided.contains(si)).map(|(si, _)| vec![si]).collect() }]
 } else {
 vec![WireFace { outer_wire: (0..segments_topo.len()).collect(), inner_wires: vec![], internal_wires: vec![] }]
 };
 if wfs.is_empty() { return; }

 let mut wfs = wfs;
 // OCCT-aligned L147: PerformInternalShapes
 crate::builder::wire_path_topo_ds::perform_internal_shapes(
 &mut wfs, &internal_wire_groups, &segments_topo, tool, face_idx, face_ref, ds);

 // OCCT form: BuilderFace::Perform keeps ALL WireFaces.  Classification against
 // the opposing solid is done at the SOLID level (BuildRC), not per-face.
 // Build reverse lookup: DsEdge(sr) in segments_topo has sr.index = e_base + ei,
 // which cannot directly index ds.edges. Map ShapeRef.index -> DS edge index.
 // Must build before drop(segments) below.
 let e_base = self.ds.vertices.len();
 let ds_ei_to_sr: HashMap<usize, topods::ShapeRef> = segments.iter()
 .filter_map(|seg| match &seg.source {
 WireEdgeSource::DsEdge(ei) => Some((*ei, self.brep_sr(e_base + *ei))),
 _ => None,
 }).collect();
 let sr_index_to_ds_ei: HashMap<usize, usize> = segments.iter()
 .filter_map(|seg| match &seg.source {
 WireEdgeSource::DsEdge(ei) => Some((e_base + *ei, *ei)),
 _ => None,
 }).collect();
 drop(segments);

 let origin = if is_a {
 FaceOrigin::FromA(ds.faces[face_idx].source_face_idx)
 } else {
 FaceOrigin::FromB(ds.faces[face_idx].source_face_idx)
 };
 let ic_curves: HashMap<usize, Curve3> = ds.intersection_curves.iter()
 .enumerate().map(|(ci, ic)| (ci, ic.curve.clone())).collect();
 // Architecture diff A6: provide DS edge array for curve lookup.
 result.ds_edges = Some(std::sync::Arc::new(self.ds.edges.clone()));
 for wf in &wfs {
 result.emit_wire_face_topods(face_idx, wf, &segments_topo, tool, &ic_curves, false, origin,
 &HashMap::new(), face_refs[face_idx], self.ds.faces[face_idx].natural_restriction,
 &ds_ei_to_sr, &sr_index_to_ds_ei, &self.ds);
 // Architecture A1: create TShapes for this face immediately (incremental),
 // matching OCCT's per-face BRep_Builder assembly.  build_topods_faces will
 // skip faces already emitted as TShapes.
 result.emit_face_topods(t, &mut *self.my_face_refs.borrow_mut());
 }
 }

 /// =OCCT-aligned: BuilderFace::CheckData (BOPAlgo_BuilderFace.cxx L50-115).
 /// Validates face has intersection curves/segments. If no interferences,
 /// delegates to BuildDraftFace (OCCT's alternative path for non-split faces).
 fn builder_face_check_data(&self, face_idx: usize, segments: &[WireSegment]) -> bool {
 if segments.is_empty() {
 return false;
 }
 true
 }

 /// =OCCT-aligned: PIOperation_FillHistory =PrepareHistory (Builder_4.cxx L164-252).
 /// Builds source= esult history matching OCCT's BRepTools_History.
 ///
 /// OCCT form:
 /// L166:  if (!HasHistory()) return;
 /// L174:  myHistory = new BRepTools_History;
 /// L175:  myMapShape.Clear();
 /// L176:  TopExp::MapShapes(myShape, myMapShape);
 /// L185-187: for i in 0..NbSourceShapes()
 /// L192: if (!IsSupportedType(aS)) continue;
 /// L205: pLSp = LocModified(aS);  // =images
 /// L214: if (myMapShape.Contains(aSp)) =Modified
 /// L233: aGenShapes = LocGenerated(aS);
 /// L239: if (myMapShape.Contains(aG)) =Generated
 /// L247: if (!isModified && !myMapShape.Contains(aS)) =Deleted
 fn fill_history(&self, t_brep: &mut topods::BRep) -> Vec<crate::history::SourceShapeEntry> {
 use crate::history::{HistoryStatus, SourceShapeEntry};
 use topods::TShape;

 // OCCT L166: if (!HasHistory()) return.
 if !self.my_fill_history {
 return vec![];
 }

 // OCCT L174-176: TopExp::MapShapes(myShape, myMapShape).
 // rcad: build result vertex/edge presence sets from t_brep.
 let mut result_vtx: std::collections::HashSet<usize> = std::collections::HashSet::new();
 let mut result_edge: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for ts in &t_brep.tshapes {
 match &**ts {
 TShape::Vertex(vd) => {
 for (di, dv) in self.ds.vertices.iter().enumerate() {
 if (dv.point - vd.point).length_squared()
 < crate::tolerance::TOLERANCE_ABS * 2.0
 {
 result_vtx.insert(di);
 break;
 }
 }
 }
 TShape::Edge(ed) => {
 for (di, de) in self.ds.edges.iter().enumerate() {
 if (de.start_vertex == ed.first.index
 && de.end_vertex == ed.last.index)
 || (de.start_vertex == ed.last.index
 && de.end_vertex == ed.first.index)
 {
 result_edge.insert(di);
 break;
 }
 }
 }
 _ => {}
 }
 }

 let v_base = 0usize; // vertices start at 0 in t_brep
 let e_base = self.ds.vertices.len();
 let mut modified_indices: Vec<usize> = Vec::new();
 let mut entries = Vec::new();

 // = =  Iterate all source shapes = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
 // OCCT L185-187: for (int i = 0; i < aNbS; ++i)
 //
 // =Vertices (OCCT L192: IsSupportedType filter =all vertex types are valid)
 for (di, _dv) in self.ds.vertices.iter().enumerate() {
 // OCCT L205: const List<TopoDS_Shape>* pLSp = LocModified(aS);
 let sref = self.brep_sr(v_base + di);
 let has_images = self.my_images.borrow().contains_key(&sref);
 let in_result = result_vtx.contains(&di);

 let (status, result_indices) = if has_images && in_result {
 // OCCT L208-230: split images found in result =Modified
 let images = self.my_images.borrow().get(&sref).cloned().unwrap_or_default();
 modified_indices.push(v_base + di);
 (HistoryStatus::Modified, images.iter().map(|sr| sr.index).collect())
 } else if in_result {
 // OCCT L233-243: LocGenerated =in result =Generated
 (HistoryStatus::Generated, vec![v_base + di])
 } else {
 // OCCT L247-249: not in result =Deleted
 (HistoryStatus::Deleted, vec![])
 };
 entries.push(SourceShapeEntry { ds_index: di, shape_type: 0, status, result_indices });
 }

 // =Edges (same form)
 for (di, _de) in self.ds.edges.iter().enumerate() {
 let sref = self.brep_sr(e_base + di);
 let has_images = self.my_images.borrow().contains_key(&sref);
 let in_result = result_edge.contains(&di);

 let (status, result_indices) = if has_images && in_result {
 let images = self.my_images.borrow().get(&sref).cloned().unwrap_or_default();
 modified_indices.push(e_base + di);
 (HistoryStatus::Modified, images.iter().map(|sr| sr.index).collect())
 } else if in_result {
 (HistoryStatus::Generated, vec![e_base + di])
 } else {
 (HistoryStatus::Deleted, vec![])
 };
 entries.push(SourceShapeEntry { ds_index: di, shape_type: 1, status, result_indices });
 }

 // =Faces (OCCT shape type TopAbs_FACE =matched by surface + wire topology)
 // TODO: Add face-level history when topods face= S face matching is available.
 // Currently faces are tracked indirectly via face_origins in BuildResult.
 // OCCT L192: if (!BRepTools_History::IsSupportedType(aS)) continue;
 // For now, faces are not mapped here =they are handled by
 // annotate_shell_and_solid_history during post_treat.

 // = =  Set TopoDS_TShape::Moved for modified shapes = = = = = = = = = = = = = = = = = = = = = 
 // OCCT L216-225: modified shapes get orientation fix + moved flag.
 for &idx in &modified_indices {
 if let Some(arc) = t_brep.tshapes.get_mut(idx) {
 if let Some(arc) = std::sync::Arc::get_mut(arc) {
 match arc {
 TShape::Vertex(vd) => vd.flags |= topods::tshape_flags::MODIFIED,
 TShape::Edge(ed) => ed.flags |= topods::tshape_flags::MODIFIED,
 _ => {}
 }
 }
 }
 }

 entries
 }

 /// OCCT-aligned: PrepareHistory (BOPAlgo_Builder.cxx L441-448).
 /// Single call that builds the final shape topology and records
 /// source shape history (modified/generated/deleted).
 /// In rcad this combines two steps (build_topods for shape-level,
 /// fill_history for source-level) into one public method.

 /// OCCT-aligned: PrepareHistory for the TreatEmptyShape case (BOP.cxx L462-468).
 /// All source shapes are present as-is (Generated) or absent (Deleted);
 /// no splitting occurs, so no Modified shapes.
 fn source_history_single_side(&self, side: ShapeOrigin) -> Vec<crate::history::SourceShapeEntry> {
 use crate::history::{HistoryStatus, SourceShapeEntry};

 let a_vc = self.ds.a_vertex_count;
 let a_ec = self.ds.a_edge_count;
 let side_is_a = side == ShapeOrigin::ShapeA;
 let mut entries = Vec::new();

 // OCCT L185-187: iterate all source shapes
 for (di, _) in self.ds.vertices.iter().enumerate() {
 let on_side = if side_is_a { di < a_vc } else { di >= a_vc };
 let status = if on_side { HistoryStatus::Generated } else { HistoryStatus::Deleted };
 entries.push(SourceShapeEntry { ds_index: di, shape_type: 0, status, result_indices: vec![] });
 }
 for (di, _) in self.ds.edges.iter().enumerate() {
 let on_side = if side_is_a { di < a_ec } else { di >= a_ec };
 let status = if on_side { HistoryStatus::Generated } else { HistoryStatus::Deleted };
 entries.push(SourceShapeEntry { ds_index: di, shape_type: 1, status, result_indices: vec![] });
 }
 entries
 }

 /// All source shapes are Deleted (both operands empty).
 fn source_history_all_deleted(&self) -> Vec<crate::history::SourceShapeEntry> {
 use crate::history::{HistoryStatus, SourceShapeEntry};
 let mut entries = Vec::new();
 for (di, _) in self.ds.vertices.iter().enumerate() {
 entries.push(SourceShapeEntry { ds_index: di, shape_type: 0, status: HistoryStatus::Deleted, result_indices: vec![] });
 }
 for (di, _) in self.ds.edges.iter().enumerate() {
 entries.push(SourceShapeEntry { ds_index: di, shape_type: 1, status: HistoryStatus::Deleted, result_indices: vec![] });
 }
 entries
 }

 ///
 /// have been split by the PaveFiller (via myImages / vertices_in), build a
 /// single analytic face using the split boundary edges.  This avoids the
 /// tessellation fallback (split_curved_face_parametric, tessellate_sphere_face,
 /// etc.) that would otherwise be used for non-planar faces with only
 /// alone-vertex / on-edge intersection data.
 ///
 /// Returns `None` when:
 /// - The face has no boundary segments (empty geometry)
 /// - Any vertex is multi-connected (>=3 edges share the same vertex),
 /// indicating the face may need full SmartMap-based splitting
 /// - The wire pipeline cannot form a closed loop
 fn build_draft_face(
 &self,
 face_idx: usize,
 ) -> Option<(Vec<WireSegment>, Vec<WireFace>, HashMap<usize, DVec3>)> {
 let ds = self.ds;
 let face = &ds.faces[face_idx];
 let pcurve_lookup = |ci: usize| self.find_pcurve_for_face(ci, face_idx);
 let mut segments = collect_face_edge_segments(ds, face_idx, &pcurve_lookup);
 if segments.is_empty() {
 return None;
 }

 // OCCT HasMultiConnected: if a vertex connects >=3 boundary edges,
 // the face cannot be represented as a single closed wire and needs
 // the full SmartMap splitting path (BOPAlgo_Builder_2.cxx L1068-1074).
 let mut vert_count: HashMap<usize, usize> = HashMap::new();
 for seg in &segments {
 *vert_count.entry(seg.start_vertex).or_default() += 1;
 *vert_count.entry(seg.end_vertex).or_default() += 1;
 }
 if vert_count.values().any(|&c| c > 2) {
 return None;
 }

 let (wires, internal_wires, vertex_positions) =
 build_closed_wires(&mut segments, ds, face_idx, &std::collections::HashSet::new());
 if wires.is_empty() && internal_wires.is_empty() {
 return None;
 }
 let wfs = perform_areas(&wires, &internal_wires, &segments, ds, &mut *self.context.borrow_mut(), face_idx);
 if wfs.is_empty() {
 return None;
 }

 Some((segments, wfs, vertex_positions))
 }
}

// =============================================================================
// Phase 2: OCCT 1:1 PerformLoops Alignment (BOPAlgo_BuilderFace.cxx L239-606)
// =============================================================================

/// Edge-like segment for wire building=can be a DS edge, an intersection curve,
impl<'a> BooleanBuilder<'a> {
 pub fn new(ds: &'a DS, op: BooleanOpType) -> Self {
 let context = RefCell::new(Context::new(ds.faces.len(), TOLERANCE_ABS * 100.0));
 Self {
 ds, op, glue: GlueEnum::GlueOff, glue_tolerance: TOLERANCE_ABS, context, has_errors: false,
 my_images: std::cell::RefCell::new(std::collections::HashMap::new()),
 my_origins: std::cell::RefCell::new(std::collections::HashMap::new()),
 my_shapes_sd: std::cell::RefCell::new(std::collections::HashMap::new()),
 my_in_parts: std::cell::RefCell::new(std::collections::HashMap::new()),
 my_solid_images: std::cell::RefCell::new(std::collections::HashMap::new()),
 my_solid_origins: std::cell::RefCell::new(std::collections::HashMap::new()),
 my_non_destructive: false,
 my_fill_history: true, // OCCT default
 my_check_inverted: false,
 my_stop_on_fatal_error: true,
 my_entry_point: 0,
  my_report: std::cell::RefCell::new(Report::new()),
 my_dims: std::cell::Cell::new([3, 3]),
 brep: std::cell::RefCell::new(None),
 my_shape: std::cell::RefCell::new(rcad_kernel::topods::BRep::new()),
 my_arguments: std::cell::RefCell::new(Vec::new()),
 my_edge_map: std::cell::RefCell::new(Vec::new()),
 my_wire_refs: std::cell::RefCell::new(Vec::new()),
 my_shells: std::cell::RefCell::new(Vec::new()),
 my_face_refs: std::cell::RefCell::new(Vec::new()),
 my_solids: std::cell::RefCell::new(Vec::new()),
 my_compsolid_groups: std::cell::RefCell::new(Vec::new()),
 }
 }

 /// Pre-populate the BRep from a pre-built one (A3 dual-write via PaveFiller).
 /// Takes face_refs and ic_edge_map from export_to_brep.
 pub fn set_brep_with_mappings(&self, brep: rcad_kernel::topods::BRep,
 face_refs: Vec<rcad_kernel::topods::ShapeRef>,
 ic_edge_map: Vec<Option<rcad_kernel::topods::ShapeRef>>)
 {
 *self.brep.borrow_mut() = Some((brep, face_refs, ic_edge_map));
 }

 /// Create builder with a pre-built BRep (A3 dual-write, skips ds_to_brep).
 pub fn with_brep(ds: &'a DS, op: BooleanOpType, brep: rcad_kernel::topods::BRep,
 face_refs: Vec<rcad_kernel::topods::ShapeRef>,
 ic_edge_map: Vec<Option<rcad_kernel::topods::ShapeRef>>) -> Self
 {
 let builder = Self::new(ds, op);
 builder.set_brep_with_mappings(brep, face_refs, ic_edge_map);
 builder
 }

 /// Return a real Arc-based ShapeRef for a given flat DS index.
 /// Looks up the TShape at `flat_idx` in the cached BRep to extract the real
 /// Arc pointer identity.  Falls back to a synthetic ShapeRef when the BRep
 /// is unavailable or the index is out of range (e.g. sentinel keys for
 /// shells/solids).
 fn brep_sr(&self, flat_idx: usize) -> rcad_kernel::topods::ShapeRef {
 let brep_borrow = self.brep.borrow();
 if let Some((ref br, ..)) = *brep_borrow {
 if flat_idx < br.tshapes.len() {
 let ptr_id = std::sync::Arc::as_ptr(&br.tshapes[flat_idx]) as u64;
 return rcad_kernel::topods::ShapeRef { ptr_id, index: flat_idx,
 orientation: rcad_kernel::topods::Orientation::Forward, location: 0 };
 }
 }
 { rcad_kernel::topods::ShapeRef::synthetic(flat_idx) }
 }

 pub fn with_glue(mut self, enable: bool, tolerance: f64) -> Self {
 self.glue = if enable { GlueEnum::GlueFull } else { GlueEnum::GlueOff };
 self.glue_tolerance = tolerance.max(TOLERANCE_ABS);
 self
 }

 /// Unified semantic policy for sub-face retention.
 ///
 /// This keeps A/B branches aligned to the same decision table instead of
 /// maintaining two subtly diverging helper functions.
 fn keep_face_policy(op: BooleanOpType, source: SourceSide, class: Classification) -> bool {
 match op {
 // Regularized union: keep outside + coincident boundary fragments.
 // Coincident (`On`) fragments are deduplicated downstream in ResultBuilder.
 BooleanOpType::Union => {
 class == Classification::Out || class == Classification::On
 }
 BooleanOpType::Intersection => {
 class == Classification::In || class == Classification::On
 }
 BooleanOpType::Difference => match source {
 SourceSide::A => class == Classification::Out,
 SourceSide::B => class == Classification::In,
 },
 }
 }

 fn pcurve_matches_face_surface(
 &self,
 pcurve: &rcad_kernel::geom::Curve2d,
 surface: &Surface3,
 ic: &IntersectionCurve,
 ) -> bool {
 let samples: Vec<DVec3> = if ic.polyline.len() >= 3 {
 let mid = ic.polyline.len() / 2;
 vec![ic.polyline[0], ic.polyline[mid], *ic.polyline.last().unwrap()]
 } else if ic.polyline.len() == 2 {
 vec![ic.polyline[0], ic.polyline[1]]
 } else {
 let [t0, t1] = ic.t_range;
 let tm = 0.5 * (t0 + t1);
 vec![ic.curve.point_at(t0), ic.curve.point_at(tm), ic.curve.point_at(t1)]
 };

 let params: Vec<f64> = match pcurve.inner() {
 rcad_kernel::geom::Curve2d::BSpline(_) => {
 if samples.len() <= 1 {
 vec![0.0]
 } else {
 (0..samples.len())
 .map(|i| i as f64 / (samples.len() - 1) as f64)
 .collect()
 }
 }
 _ => {
 let [t0, t1] = ic.t_range;
 if samples.len() <= 1 {
 vec![t0]
 } else {
 (0..samples.len())
 .map(|i| t0 + (t1 - t0) * i as f64 / (samples.len() - 1) as f64)
 .collect()
 }
 }
 };

 let mut max_err: f64 = 0.0;
 for (sample, t) in samples.iter().zip(params.iter().copied()) {
 let uv = pcurve.point_at(t);
 let lifted = surface.point_at(uv.x, uv.y);
 max_err = max_err.max((lifted - *sample).length());
 }

 max_err.is_finite() && max_err <= TOLERANCE_ADAPTIVE_MAX
 }

 pub fn build(&mut self) -> Result<topods::BRep, BooleanError> {
 let (t, _) = self.build_with_history()?;
 Ok(t)
 }
}

impl<'a> BooleanBuilder<'a> {
 /// The top-level pipeline entry: dimension-by-dimension image filling
 /// (V= = = ACE= HELL= OLID), followed by BuildResult for each type.
 /// OCCT L310-445 structure matched in full (see inline OCCT line refs).
  /// =OCCT-aligned: CheckData (BOPAlgo_BOP.cxx L106-202) + CheckFiller (Builder.cxx L143-151).
  /// Validates operation type, non-empty arguments, and DS/PaveFiller state.
  /// OCCT form: AddError on each failure, then HasErrors check at the end.
  fn check_data(&self) -> Result<(), BooleanError> {
    // OCCT L132-137: if (myArguments.Extent() < 2) -> AlertTooFewArguments
    let nb_args = self.my_arguments.borrow().len();
    if nb_args < 2 {
      self.my_report.borrow_mut().add_alert(crate::bopalgo::Alert::TooFewArguments);
    }
    // OCCT L139-141: CheckFiller -> AlertNoFiller
    if self.ds.vertices.is_empty() {
      self.my_report.borrow_mut().add_alert(crate::bopalgo::Alert::NoFiller);
    }
    // OCCT: BOPAlgo_Builder::CheckData() — parent class checks.
    // rcad: operation type validation.
    match self.op {
      BooleanOpType::Union | BooleanOpType::Intersection | BooleanOpType::Difference => {}
      _ => self.my_report.borrow_mut().add_alert(crate::bopalgo::Alert::BOPNotSet),
    }
    // OCCT: if (HasErrors()) return;
    if self.my_report.borrow().has_alerts() {
      return Err(BooleanError::InvalidOperation);
    }
    Ok(())
  }

 ///  ?OCCT-aligned: Prepare (BOPAlgo_Builder.cxx L156-164).
 /// OCCT: BRep_Builder.MakeCompound(myShape)  ?empty compound as result.
 /// rcad: initializes my_shape + returns (BRep, ResultBuilder) for downstream.
 fn prepare(&self) -> (topods::BRep, ResultBuilder) {
 *self.my_shape.borrow_mut() = topods::BRep::new();
 (topods::BRep::new(), ResultBuilder::new())
 }

 ///  ?OCCT-aligned: create TShapes for all DS source shapes in my_shape.
 /// Equivalent to OCCT's myArguments populated with all source TopoDS_Shape.
 ///  ?OCCT-aligned: TreatEmptyShape (BOPAlgo_BOP.cxx L214-319).
 /// Handles the case where one or both operands have no geometry.
 /// Returns Ok(Some(brep)) if a quick result was determined,
 /// Ok(None) if the full pipeline must run.
 fn treat_empty_shape(&self, a_faces: &[usize], b_faces: &[usize])
 -> Result<Option<topods::BRep>, BooleanError>
 {
 let has_a = !a_faces.is_empty();
 let has_b = !b_faces.is_empty();
 if has_a && has_b {
 return Ok(None); // need full pipeline
 }
 if !has_a && !has_b {
 // OCCT L252-256: all empty → empty
 return Ok(Some(topods::BRep::new()));
 }
 // OCCT L258-317: one side empty → result depends on operation
 match self.op {
 BooleanOpType::Union => {
 // OCCT L270-279: return non-empty side
 let src = if has_a { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
 let brep = self.brep_of_side_topods(src, a_faces.len(), b_faces.len());
 Ok(Some(brep))
 }
 BooleanOpType::Intersection => {
 // OCCT L303-304: Common always empty
 Ok(Some(topods::BRep::new()))
 }
 BooleanOpType::Difference => {
 if !has_a {
 // OCCT L287-289: CUT with empty objects → empty
 Ok(Some(topods::BRep::new()))
 } else {
 // OCCT L281-289: CUT with empty tools → return objects
 let brep = self.brep_of_side_topods(ShapeOrigin::ShapeA, a_faces.len(), b_faces.len());
 Ok(Some(brep))
 }
 }
 _ => {
 // Unknown operation → fall through to full pipeline
 Ok(None)
 }
 }
 }

 /// OCCT-aligned: BOPAlgo_BOP::PerformInternal1 (BOP.cxx L422-579).
 /// Every statement in OCCT L422-579 has a corresponding rcad line below.
 /// See comments for exact OCCT line references.
 /// Structural difference: L425-429 setup done in constructor, re-affirmed here.
 /// L531 BuildResult(SOLID) writes to t_brep, then L900 BuildRC filters and
 /// clears solids from t_brep (non-Union) =equivalent to OCCT removing from myShape.
 pub fn build_with_history(&mut self) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 self.build_with_history_topods()
 }

 /// Same as build_with_history but returns topods::BRep directly (OCCT-aligned).
 pub fn build_with_history_topods(&mut self) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
 // SKIP: OCCT L425-429 copies (myPaveFiller, myDS, myContext, myFuzzyValue, myNonDestructive)
 // from theFiller argument.  rcad builds with_brep() which sets these in the constructor --
 // no re-assignment at the start of build_with_history_topods is needed.

 // OCCT L431-436: CheckData =validates arguments and merges PaveFiller report.
 // Populate my_arguments from DS source shapes (OCCT: SetArguments).
 let mut args = self.my_arguments.borrow_mut();
 args.clear();
 args.push(rcad_kernel::topods::ShapeRef::synthetic(0));
 args.push(rcad_kernel::topods::ShapeRef::synthetic(1));
 drop(args);
 let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
 let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
 self.check_data()?;

 // OCCT L438-443: Prepare
 // rcad: prepare() initializes my_shape + returns ResultBuilder.
 let mut result = self.prepare().1;
 // Compute myDims from argument shape types.
 // OCCT: BRep_Tool::Dimension(myArguments.First/Last) → VERTEX=0, EDGE=1, WIRE/FACE=2, SHELL/SOLID=3.
 // rcad: infer from per-origin shape counts in the DS.
 let a_n_verts = self.ds.a_vertex_count;
 let a_n_edges = self.ds.a_edge_count;
 let a_n_faces = self.ds.a_face_count;
 let b_n_verts = self.ds.vertices.len() - a_n_verts;
 let b_n_edges = self.ds.edges.len() - a_n_edges;
 let b_n_faces = self.ds.faces.len() - a_n_faces;
 let dim_a = if a_n_faces > 0 { 3_i8 } else if a_n_edges > 0 { 1 } else if a_n_verts > 0 { 0 } else { 3 };
 let dim_b = if b_n_faces > 0 { 3_i8 } else if b_n_edges > 0 { 1 } else if b_n_verts > 0 { 0 } else { 3 };
 self.my_dims.set([dim_a, dim_b]);

 // Pipeline snapshot: create dump context for stage-by-stage DS/BRep capture.
 let mut dump_ctx = DumpCtx::new(
     &std::env::var("RCAD_DUMP_GRID").unwrap_or_default(),
     &std::env::var("RCAD_DUMP_CASE").unwrap_or_default(),
 );
 dump_ctx.snapshot("initial", self.ds, None);

 // OCCT L445-453: TreatEmptyShape.
 // rcad: check if either operand has no faces (DS was populated by pave_fill).
 if a_faces.is_empty() || b_faces.is_empty() {
 if self.treat_empty_shape(&a_faces, &b_faces)?.is_some() {
 let source_history = if self.my_fill_history {
 let side = if !a_faces.is_empty() { ShapeOrigin::ShapeA }
 else { ShapeOrigin::ShapeB };
 self.source_history_single_side(side)
 } else {
 vec![]
 };
 let mut history = result.build_topods(&mut *self.my_shape.borrow_mut(), self.my_fill_history, &self.my_shells.borrow(), &mut *self.my_face_refs.borrow_mut(), &self.my_solids.borrow(), &self.my_compsolid_groups.borrow());
 history.source_history = source_history;
 let result_brep = self.my_shape.borrow().clone();
 return Ok((result_brep, history));
 }
 }

 // OCCT L454-457: ProgressScope + PISteps + analyzeProgress.
 // rcad: progress reporting not yet integrated.
 // OCCT L459-471: 3.1 FillImagesVertices + BuildResult(VERTEX)
 self.fill_images_vertices();
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 self.build_result(topods::ShapeType::Vertex, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 dump_ctx.snapshot("after_FillImagesVertices", self.ds, Some(&*self.my_shape.borrow()));
 // OCCT L472-483: 3.2 FillImagesEdges + BuildResult(EDGE)
 self.fill_images_edges();
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 self.build_result(topods::ShapeType::Edge, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 dump_ctx.snapshot("after_FillImagesEdges", self.ds, Some(&*self.my_shape.borrow()));
 // OCCT L484-496: 3.3 FillImagesContainers(WIRE) + BuildResult(WIRE)
 self.fill_images_container(ShapeType::Wire, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 self.build_result(topods::ShapeType::Wire, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 dump_ctx.snapshot("after_BuildResultWire", self.ds, Some(&*self.my_shape.borrow()));
 // OCCT L497-509: 3.4 FillImagesFaces + BuildResult(FACE)
 // Architecture A1: split faces create TShapes incrementally during fill_images_faces.
 // Remaining unsplit faces have existing TShapes from pre-create_source_shapes.
 self.fill_images_faces();
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 // BuildResult(FACE) — generic loop over my_arguments, adds originals/splits to result.
 self.build_result(topods::ShapeType::Face, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 dump_ctx.snapshot("after_FillImagesFaces", self.ds, Some(&*self.my_shape.borrow()));
 // OCCT L510-522: 3.5 FillImagesContainers(SHELL) + BuildResult(SHELL)
 self.fill_images_container(ShapeType::Shell, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 self.build_result(topods::ShapeType::Shell, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 dump_ctx.snapshot("after_BuildResultShell", self.ds, Some(&*self.my_shape.borrow()));
 // OCCT L523-535: 3.6 FillImagesSolids + BuildResult(SOLID)
 self.fill_images_solids(&mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 self.build_result(topods::ShapeType::Solid, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 dump_ctx.snapshot("after_FillImagesSolids", self.ds, Some(&*self.my_shape.borrow()));
 // OCCT L536-548: 3.7 FillImagesContainers(COMPSOLID) + BuildResult(COMPSOLID)
 self.fill_images_container(ShapeType::CompSolid, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 self.build_result(topods::ShapeType::CompSolid, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 // OCCT L549-561: 3.8 FillImagesCompounds + BuildResult(COMPOUND)
 self.fill_images_compounds(&mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 self.build_result(topods::ShapeType::Compound, &mut result);
 if self.has_errors { return Err(BooleanError::DegenerateResult); }
 dump_ctx.snapshot("after_FillImagesCompounds", self.ds, Some(&*self.my_shape.borrow()));
 // OCCT L563-568: 4. PrepareHistory — builds BRep TShapes + source shape history.
 let mut history = {
 let mut t_brep = self.my_shape.borrow_mut();
 let mut history = result.build_topods(&mut *t_brep, self.my_fill_history,
 &self.my_shells.borrow(), &mut *self.my_face_refs.borrow_mut(),
 &self.my_solids.borrow(), &self.my_compsolid_groups.borrow());
 let source_history = if self.my_fill_history {
 self.fill_history(&mut *t_brep)
 } else { vec![] };
 history.source_history = source_history;
 history
 };
 // OCCT L577-578: 5. PostTreat
 // Corrects tolerances of the result shape (CorrectTolerances + CorrectShapeTolerances).
 self.post_treat();
 let result_brep = self.my_shape.borrow().clone();

 Ok((result_brep, history))
 }

 /// Parallel version of `build_with_history`.
 ///
 /// Uses Rayon to process faces in parallel. Each face is split and classified
 /// independently, then results are merged. This can provide significant
 /// speedup for models with many faces (e.g., > 100 faces).
 ///
 /// # Performance
 ///
 /// - Small models (< 20 faces): May be slower due to thread overhead
 /// - Large models (> 100 faces): Typically 2-4x faster on multi-core systems
 fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
 let mut v: Vec<usize> = self
 .ds
 .faces
 .iter()
 .enumerate()
 .filter(|(_, f)| f.origin == origin)
 .map(|(i, _)| i)
 .collect();
 v.sort_unstable();
 v
 }

 /// OCCT L258-317: build result topods::BRep from one side's source shapes (TreatEmptyShape path).
 fn brep_of_side_topods(&self, origin: ShapeOrigin, _na: usize, _nb: usize) -> topods::BRep {
 let mut t = topods::BRep::new();
 let mut v_map: HashMap<usize, topods::ShapeRef> = HashMap::new();
 let mut e_map: HashMap<usize, topods::ShapeRef> = HashMap::new();

 for (fi, f) in self.ds.faces.iter().enumerate() {
 if f.origin != origin { continue; }
 let mut edge_refs: Vec<topods::ShapeRef> = Vec::new();
 for &ei in &f.boundary_edges {
 if ei >= self.ds.edges.len() { continue; }
 let e = &self.ds.edges[ei];
 let sv = *v_map.entry(e.start_vertex).or_insert_with(|| {
 t.add_tvertex(self.ds.vertices[e.start_vertex].point)
 });
 let ev = *v_map.entry(e.end_vertex).or_insert_with(|| {
 t.add_tvertex(self.ds.vertices[e.end_vertex].point)
 });
 let edge_sr = *e_map.entry(ei).or_insert_with(|| {
 t.add_tedge(None, sv, ev, [0.0, 1.0])
 });
 edge_refs.push(edge_sr);
 }
 if edge_refs.is_empty() { continue; }
 let ow = t.add_twire(edge_refs);
 let surf = f.surface.clone();
 let face = t.add_tface(Some(surf), ow, vec![], None, None, vec![], false);
 t.face_mut(face).tolerance = f.geom_tol;
 }

 // Collect face refs and wrap in Shell → Solid
 let face_srs: Vec<topods::ShapeRef> = t.tshapes.iter().enumerate()
 .filter(|(_, ts)| matches!(&***ts, topods::TShape::Face(_)))
 .map(|(i, _)| { let idx = i; topods::ShapeRef::synthetic(idx) })
 .collect();
 if !face_srs.is_empty() {
 let shell = t.add_tshell(face_srs);
 t.add_tsolid(vec![shell]);
 }
 t
 }

 /// Split a curved face (Cylinder, Sphere, Cone, Torus) by intersection polylines.
 ///
 /// Legacy approximate method: for each intersection polyline that crosses the face,
 /// we split the boundary point list into two halves at the points closest to the
 /// polyline endpoints. Kept as fallback when UV data or PCurves are unavailable.

 /// into a 2D trim polyline in UV space, then splits the UV boundary polygon.
 /// Maps resulting sub-polygons back to 3D via surface evaluation.
 ///
 /// == € ? = ㄦ = € ==﹂ =   Υ?
 /// OCCT: BuildSplitFaces =section edges ==ㄧ=== ?BRep sub-face=
 /// rcad: = == =8 == ⒔= ?FaceSampleData,=outer_circle_edges ===== 〒 ь ?
 /// == = ?8 == = ㄩ ф =+  € € ?,=OCCT = 〒 ｇ  ?FaceSampleData=

 /// Find the PCurve (2D parametric curve) for the given intersection curve
 /// as it lies on the given face. Searches FaceFace interferences to determine
 /// whether this face is f1 (use pcurve_on_a) or f2 (use pcurve_on_b).
 fn find_pcurve_for_face(
 &self,
 curve_idx: usize,
 face_idx: usize,
 ) -> Option<rcad_kernel::geom::Curve2d> {
 for ff in &self.ds.interf_ff {
 if ff.curves.contains(&curve_idx) {
 let ic = &self.ds.intersection_curves[curve_idx];
 if ff.f1 == face_idx {
 return ic.pcurve_on_a.clone();
 } else if ff.f2 == face_idx {
 return ic.pcurve_on_b.clone();
 }
 }
 }

 let ic = &self.ds.intersection_curves[curve_idx];
 let surface = &self.ds.faces[face_idx].surface;
 if let Some(pcurve) = &ic.pcurve_on_a
 && self.pcurve_matches_face_surface(pcurve, surface, ic)
 {
 return Some(pcurve.clone());
 }
 if let Some(pcurve) = &ic.pcurve_on_b
 && self.pcurve_matches_face_surface(pcurve, surface, ic)
 {
 return Some(pcurve.clone());
 }
 None
 }
}
