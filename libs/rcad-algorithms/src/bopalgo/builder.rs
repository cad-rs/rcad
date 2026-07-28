use indexmap::IndexMap;
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use glam::{DVec2, DVec3};
use rayon::prelude::*;
use rcad_kernel::PCurve;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval, *};
use rcad_kernel::topods;
use rcad_kernel::topology::Face;
use rcad_kernel::topology::*;

use crate::bopalgo::{Alert, GlueEnum, Report};
use crate::bopds::ds::*;
use crate::boptools::bvh::{Aabb, BoxTree};
use crate::classify::{Classification, classify_point};
use crate::history::{
    BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin,
};
use crate::inttools::context::Context;
use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;
use crate::triangulate::{triangulate_polygon, triangulate_polygon_with_holes};
use std::cell::RefCell;

use crate::pipeline_dump::DumpCtx;

mod angle_2d;
mod curve_tools;
mod debug_utils;
mod ds_as_brep;
pub(crate) mod intersection;
pub(crate) mod intres2d;
mod types;

mod builder_face;
pub(crate) use builder_face::BuilderFace;

pub use types::{BooleanError, BooleanOpType};
pub(crate) use types::{
    FaceEntry, ShapeType, WireEdgeSource, WireEdgeSourceTopoDS, WireFace, WireOrientation,
    WireSegment, WireSegmentTopoDS,
};

mod result_builder;

pub(crate) use result_builder::ResultBuilder;
mod builder_utils;

pub(crate) use builder_utils::{
    aggregate_face_region_origin, aggregate_shell_region_origin, annotate_history_from_ds,
    annotate_shell_and_solid_history, build_edge_bounds, check_and_add_split_vertex, curve_eq,
    hash_point, is_tangent_face, point_in_polygon_2d, quantize_pos,
};

// OCCT BOPAlgo_Builder.hxx L75-507 + parent class fields (BOPAlgo_BuilderShape, BOPAlgo_Options, BOPAlgo_BOP).
// Flattened into one Rust struct because Rust has no C++ inheritance.
pub struct BooleanBuilder<'a> {
    // BOPAlgo_Builder.hxx L495 =myPaveFiller — NOT stored; rcad runs PaveFiller before Builder, passes DS only.
    // BOPAlgo_Builder.hxx L496 =myDS
    pub(crate) ds: &'a DS,
    // BOPAlgo_Builder.hxx L497 =myContext
    pub(crate) context: RefCell<Context>,
    // BOPAlgo_Builder.hxx L492 =myArguments
    pub(crate) my_arguments: std::cell::RefCell<Vec<rcad_kernel::topods::Shape>>,
    // BOPAlgo_Builder.hxx L494 =myMapFence - fence map for argument uniqueness
    pub(crate) my_map_fence:
        std::cell::RefCell<std::collections::HashSet<rcad_kernel::topods::Shape>>,
    // BOPAlgo_Builder.hxx L498 =myEntryPoint - controls deletion of PaveFiller
    pub(crate) my_entry_point: u8,
    // BOPAlgo_Builder.hxx L499 =myImages - map of Images of the sub-shapes of arguments
    pub(crate) my_images: std::cell::RefCell<
        std::collections::HashMap<
            rcad_kernel::topods::Shape,
            Vec<rcad_kernel::topods::Shape>,
        >,
    >,
    // BOPAlgo_Builder.hxx L500 =myShapesSD - map of SD Shapes
    pub(crate) my_shapes_sd: std::cell::RefCell<
        std::collections::HashMap<rcad_kernel::topods::Shape, rcad_kernel::topods::Shape>,
    >,
    // BOPAlgo_Builder.hxx L501 =myOrigins - back map of Images
    pub(crate) my_origins: std::cell::RefCell<
        std::collections::HashMap<
            rcad_kernel::topods::Shape,
            Vec<rcad_kernel::topods::Shape>,
        >,
    >,
    // BOPAlgo_Builder.hxx L502 =myInParts - map of own and acquired IN faces of the arguments solids
    // rcad: keyed by source solid index (usize) ← OCCT keys by TopoDS_Shape; Rust adaptation for DS-indexed access.
    pub(crate) my_in_parts: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    /// rcad: ptr_ids of split face images classified as IN in fill_in_3d_parts.
    pub(crate) my_split_images_in:
        std::cell::RefCell<std::collections::HashMap<usize, std::collections::HashSet<u64>>>,
    // BOPAlgo_Builder.hxx L503 =myNonDestructive
    pub(crate) my_non_destructive: bool,
    // BOPAlgo_Builder.hxx L504 =myGlue
    pub(crate) my_glue: GlueEnum,
    // BOPAlgo_Builder.hxx L505 =myCheckInverted
    pub(crate) my_check_inverted: bool,

    // BOPAlgo_Options.hxx L143-147 (parent class)
    pub(crate) my_report: std::cell::RefCell<Report>,
    pub(crate) my_run_parallel: bool,
    pub(crate) my_fuzzy_value: f64,
    pub(crate) my_use_obb: bool,

    // BOPAlgo_BuilderShape.hxx L143-150 (parent class)
    pub(crate) my_shape: std::cell::RefCell<rcad_kernel::topods::BRep>,
    pub(crate) my_fill_history: bool,

    // BOPAlgo_BOP.hxx L124-126
    pub(crate) my_operation: BooleanOpType,
    pub(crate) my_dims: std::cell::Cell<[i8; 2]>,

    // --- rcad-specific (BRep flat-storage tracking, no OCCT equivalent) ---
    pub(crate) brep: std::cell::RefCell<
        Option<(
            rcad_kernel::topods::BRep,
            Vec<rcad_kernel::topods::Shape>,
            Vec<Option<rcad_kernel::topods::Shape>>,
        )>,
    >,
    pub(crate) my_edge_map: std::cell::RefCell<Vec<rcad_kernel::topods::Shape>>,
    pub(crate) my_wire_refs: std::cell::RefCell<Vec<rcad_kernel::topods::Shape>>,
    pub(crate) my_shells: std::cell::RefCell<Vec<rcad_kernel::topods::Shape>>,
    pub(crate) my_face_refs: std::cell::RefCell<Vec<rcad_kernel::topods::Shape>>,
    pub(crate) my_solids: std::cell::RefCell<Vec<rcad_kernel::topods::Shape>>,
    pub(crate) my_compsolid_groups: std::cell::RefCell<Vec<rcad_kernel::topods::Shape>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SourceSide {
    A,
    B,
}

mod builder_utils_topo_ds;
mod edge_builders;
mod filler_mod;
mod result_build_mod;
mod wire_path;
mod wire_path_topo_ds;
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
pub(crate) use wire_path::{
    intersect_ray_curve_2d, pc_parameter_range, refine_angles, walk_path_extract_wires,
};
pub(crate) use wire_splitter::{
    EdgeInfo, are_verts_coincident, build_closed_wires, edge_angle_2d, edge_uv_tangent,
    expand_avoided_pids, physical_edge_id, world_to_uv,
};

impl<'a> BooleanBuilder<'a> {
    /// OCCT BOPAlgo_Options::HasErrors - checks for fail alerts in myReport.
    fn has_errors(&self) -> bool {
        self.my_report.borrow().has_errors()
    }

    /// =PIOperation_FillHistory =PrepareHistory (Builder_4.cxx L164-252).
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
    /// OCCT BOPAlgo_Builder::PrepareHistory (BOPAlgo_Builder_4.cxx L164-252).
    fn prepare_history(&self, t_brep: &mut topods::BRep) -> Vec<crate::history::SourceShapeEntry> {
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
                    for di in 0..self.ds.vertex_count() {
                        if (self.ds.vertex_point(di) - vd.point).length_squared()
                            < crate::tolerance::TOLERANCE_ABS * 2.0
                        {
                            result_vtx.insert(di);
                            break;
                        }
                    }
                }
                TShape::Edge(ed) => {
                    for di in 0..self.ds.edge_count() {
                        let sv = self.ds.edge_start_vertex_ds(di);
                        let ev = self.ds.edge_end_vertex_ds(di);
                        if (sv == ed.first.index && ev == ed.last.index)
                            || (sv == ed.last.index && ev == ed.first.index)
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
        let e_base = self.ds.vertex_count();
        let mut modified_indices: Vec<usize> = Vec::new();
        let mut entries = Vec::new();

        // = =  Iterate all source shapes = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        // OCCT L185-187: for (int i = 0; i < aNbS; ++i)
        //
        // =Vertices (OCCT L192: IsSupportedType filter =all vertex types are valid)
        for di in 0..self.ds.vertex_count() {
            // OCCT L205: const List<TopoDS_Shape>* pLSp = LocModified(aS);
            let sref = self.brep_sr(v_base + di);
            let has_images = self.my_images.borrow().contains_key(&sref);
            let in_result = result_vtx.contains(&di);

            let (status, result_indices) = if has_images && in_result {
                // OCCT L208-230: split images found in result =Modified
                let images = self.my_images
                    .borrow()
                    .get(&sref)
                    .cloned()
                    .unwrap_or_default();
                modified_indices.push(v_base + di);
                (
                    HistoryStatus::Modified,
                    images.iter().map(|sr| sr.index).collect(),
                )
            } else if in_result {
                // OCCT L233-243: LocGenerated =in result =Generated
                (HistoryStatus::Generated, vec![v_base + di])
            } else {
                // OCCT L247-249: not in result =Deleted
                (HistoryStatus::Deleted, vec![])
            };
            entries.push(SourceShapeEntry {
                ds_index: di,
                shape_type: 0,
                status,
                result_indices,
            });
        }

        // =Edges (same form)
        for di in 0..self.ds.edge_count() {
            let sref = self.brep_sr(e_base + di);
            let has_images = self.my_images.borrow().contains_key(&sref);
            let in_result = result_edge.contains(&di);

            let (status, result_indices) = if has_images && in_result {
                let images = self.my_images
                    .borrow()
                    .get(&sref)
                    .cloned()
                    .unwrap_or_default();
                modified_indices.push(e_base + di);
                (
                    HistoryStatus::Modified,
                    images.iter().map(|sr| sr.index).collect(),
                )
            } else if in_result {
                (HistoryStatus::Generated, vec![e_base + di])
            } else {
                (HistoryStatus::Deleted, vec![])
            };
            entries.push(SourceShapeEntry {
                ds_index: di,
                shape_type: 1,
                status,
                result_indices,
            });
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

    /// PrepareHistory (BOPAlgo_Builder.cxx L441-448).
    /// Single call that builds the final shape topology and records
    /// source shape history (modified/generated/deleted).
    /// In rcad this combines two steps (build_topods for shape-level,
    /// fill_history for source-level) into one public method.

    /// PrepareHistory for the TreatEmptyShape case (BOP.cxx L462-468).
    /// All source shapes are present as-is (Generated) or absent (Deleted);
    /// no splitting occurs, so no Modified shapes.
    fn source_history_single_side(
        &self,
        side: ShapeOrigin,
    ) -> Vec<crate::history::SourceShapeEntry> {
        use crate::history::{HistoryStatus, SourceShapeEntry};

        let a_vc = self.ds.a_vertex_count();
        let a_ec = self.ds.a_edge_count();
        let side_is_a = side == ShapeOrigin::ShapeA;
        let mut entries = Vec::new();

        // OCCT L185-187: iterate all source shapes
        for di in 0..self.ds.vertex_count() {
            let on_side = if side_is_a { di < a_vc } else { di >= a_vc };
            let status = if on_side {
                HistoryStatus::Generated
            } else {
                HistoryStatus::Deleted
            };
            entries.push(SourceShapeEntry {
                ds_index: di,
                shape_type: 0,
                status,
                result_indices: vec![],
            });
        }
        for di in 0..self.ds.edge_count() {
            let on_side = if side_is_a { di < a_ec } else { di >= a_ec };
            let status = if on_side {
                HistoryStatus::Generated
            } else {
                HistoryStatus::Deleted
            };
            entries.push(SourceShapeEntry {
                ds_index: di,
                shape_type: 1,
                status,
                result_indices: vec![],
            });
        }
        entries
    }

    /// All source shapes are Deleted (both operands empty).
    fn source_history_all_deleted(&self) -> Vec<crate::history::SourceShapeEntry> {
        use crate::history::{HistoryStatus, SourceShapeEntry};
        let mut entries = Vec::new();
        for di in 0..self.ds.vertex_count() {
            entries.push(SourceShapeEntry {
                ds_index: di,
                shape_type: 0,
                status: HistoryStatus::Deleted,
                result_indices: vec![],
            });
        }
        for di in 0..self.ds.edge_count() {
            entries.push(SourceShapeEntry {
                ds_index: di,
                shape_type: 1,
                status: HistoryStatus::Deleted,
                result_indices: vec![],
            });
        }
        entries
    }

    ///
    /// have been split by the PaveFiller (via myImages / vertices_in), build a
    /// single analytic face using the split boundary edges.  This avoids the
    /// tessellation fallback (split_curved_face_parametric, tessellate_sphere_face,
    /// etc.) that would otherwise be used for non-planar faces with only
    /// ✅ OCCT-aligned: BuildDraftFace (Builder_2.cxx L1052-1189).
    /// Builds a draft face by substituting split images into the original wire.
    /// Returns None when OCCT would return a null face (INTERNAL edges,
    /// multi-connected vertices, edge unification failures).
    fn build_draft_face(
        &self,
        face_idx: usize,
    ) -> Option<(Vec<WireSegment>, Vec<WireFace>, HashMap<usize, DVec3>)> {
        let ds = self.ds;
        let e_base = ds.vertex_count();

        // OCCT L1073-1078: vertex counter / edge unification fence
        let mut edge_fence: HashSet<usize> = HashSet::new();
        let mut vert_count: HashMap<usize, usize> = HashMap::new();
        let mut segments: Vec<WireSegment> = Vec::new();
        let mut vertex_positions: HashMap<usize, DVec3> = HashMap::new();

        // Helper: add segment from an edge (original or split image)
        let mut add_segment = |sp_ei: usize,
                               sp_fwd: bool,
                               pcurve_from: Option<usize>,
                               segs: &mut Vec<WireSegment>,
                               vp: &mut HashMap<usize, DVec3>| {
            if sp_ei >= ds.edge_count() {
                return;
            }
            let sv = ds.edge_start_vertex_ds(sp_ei);
            let ev = ds.edge_end_vertex_ds(sp_ei);
            if !vp.contains_key(&sv) {
                if let Some(v) = ds.vertices.get(sv) {
                    vp.insert(sv, v.point);
                }
            }
            if !vp.contains_key(&ev) {
                if let Some(v) = ds.vertices.get(ev) {
                    vp.insert(ev, v.point);
                }
            }
            let (seg_sv, seg_ev) = if sp_fwd { (sv, ev) } else { (ev, sv) };
            let rep = ds.edge_on_face(sp_ei, face_idx)
                .or_else(|| pcurve_from.and_then(|ei| ds.edge_on_face(ei, face_idx)));
            segs.push(WireSegment {
                start_vertex: seg_sv,
                end_vertex: seg_ev,
                source: WireEdgeSource::DsEdge(sp_ei),
                orientation: WireOrientation::Forward,
                is_closed_on_face: false,
                second_pcurve: None,
                first_pcurve: rep.map(|r| r.pcurve.clone()),
                t_range: rep.map(|r| r.pcurve_range).unwrap_or(ds.edge_range(sp_ei)),
            });
        };

        // OCCT L1080-1181: iterate original face wires (outer + inner)
        let face_boundary_edges = ds.face_boundary_edges(face_idx);
        let face_boundary_edge_forwards = ds.face_boundary_edge_forwards(face_idx);
        for i in 0..face_boundary_edges.len() {
            let ei = face_boundary_edges[i];
            if ei >= ds.edge_count() {
                continue;
            }
            let forward = face_boundary_edge_forwards.get(i).copied().unwrap_or(true);

            // OCCT L1104-1110: INTERNAL edge -> return null
            if ds.edge_is_internal(ei) {
                return None;
            }
            let b_is_degenerated = ds.is_edge_degenerated(ei);

            // OCCT L1118: theImages.Seek(aE)
            let e_sr = self.brep_sr(e_base + ei);
            let has_images = self.my_images.borrow().contains_key(&e_sr);

            if !has_images {
                // OCCT L1120-1135: edge without split images
                if !b_is_degenerated {
                    *vert_count.entry(ds.edge_start_vertex_ds(ei)).or_default() += 1;
                    *vert_count.entry(ds.edge_end_vertex_ds(ei)).or_default() += 1;
                }
                // OCCT L1128: edge unification (aMEdges.Add)
                if !edge_fence.insert(ei) {
                    return None;
                }
                // OCCT L1133: aBB.Add(aNewWire, aE)
                add_segment(ei, forward, None, &mut segments, &mut vertex_positions);
            } else {
                // OCCT L1137-1175: edge has split images
                let imgs = self.my_images
                    .borrow()
                    .get(&e_sr)
                    .cloned()
                    .unwrap_or_default();
                for sp_sr in &imgs {
                    let sp_ei = sp_sr.index.saturating_sub(e_base);
                    if sp_ei >= ds.edge_count() {
                        continue;
                    }
                    // OCCT L1143: HasMultiConnected
                    if !b_is_degenerated {
                        *vert_count.entry(ds.edge_start_vertex_ds(sp_ei))
                            .or_default() += 1;
                        *vert_count.entry(ds.edge_end_vertex_ds(sp_ei)).or_default() += 1;
                    }
                    // OCCT L1149: edge unification
                    if !edge_fence.insert(sp_ei) {
                        return None;
                    }
                    // OCCT L1154: aSp.Orientation(anOriE)
                    // OCCT L1155-1159: degenerated -> add as-is
                    if b_is_degenerated {
                        add_segment(sp_ei, forward, None, &mut segments, &mut vertex_positions);
                        continue;
                    }
                    // OCCT L1161-1166: closed on face -> DoSplitSEAMOnFace (simplified: skip for draft)
                    let needs_rev =
                        crate::bopalgo::builder::edge_builders::is_split_to_reverse(ds, sp_ei, ei);
                    add_segment(
                        sp_ei,
                        forward != needs_rev,
                        Some(ei),
                        &mut segments,
                        &mut vertex_positions,
                    );
                }
            }
        }

        // OCCT: inner wires processed by same TopoDS_Iterator
        let face_inner_boundary = ds.face_inner_boundary(face_idx);
        for inner_wire in face_inner_boundary {
            for &(ei, forward) in inner_wire {
                if ei >= ds.edge_count() {
                    continue;
                }
                if ds.edge_is_internal(ei) {
                    return None;
                }
                let b_is_degenerated = ds.is_edge_degenerated(ei);
                let e_sr = self.brep_sr(e_base + ei);
                let has_images = self.my_images.borrow().contains_key(&e_sr);

                if !has_images {
                    if !b_is_degenerated {
                        *vert_count.entry(ds.edge_start_vertex_ds(ei)).or_default() += 1;
                        *vert_count.entry(ds.edge_end_vertex_ds(ei)).or_default() += 1;
                    }
                    if !edge_fence.insert(ei) {
                        return None;
                    }
                    add_segment(ei, forward, None, &mut segments, &mut vertex_positions);
                } else {
                    let imgs = self.my_images
                        .borrow()
                        .get(&e_sr)
                        .cloned()
                        .unwrap_or_default();
                    for sp_sr in &imgs {
                        let sp_ei = sp_sr.index.saturating_sub(e_base);
                        if sp_ei >= ds.edge_count() {
                            continue;
                        }
                        if !b_is_degenerated {
                            *vert_count.entry(ds.edge_start_vertex_ds(sp_ei))
                                .or_default() += 1;
                            *vert_count.entry(ds.edge_end_vertex_ds(sp_ei)).or_default() += 1;
                        }
                        if !edge_fence.insert(sp_ei) {
                            return None;
                        }
                        if b_is_degenerated {
                            add_segment(sp_ei, forward, None, &mut segments, &mut vertex_positions);
                            continue;
                        }
                        // OCCT L1161-1166: closed on face -> DoSplitSEAMOnFace (simplified: skip for draft)
                        let needs_rev = crate::bopalgo::builder::edge_builders::is_split_to_reverse(
                            ds, sp_ei, ei,
                        );
                        add_segment(
                            sp_ei,
                            forward != needs_rev,
                            Some(ei),
                            &mut segments,
                            &mut vertex_positions,
                        );
                    }
                }
            }
        }

        // OCCT L1082: check multi-connected vertices
        if vert_count.values().any(|&c| c > 2) {
            return None;
        }
        if segments.is_empty() {
            return None;
        }

        // Single WireFace with one outer wire (OCCT: one draft face per face)
        let wf = WireFace {
            outer_wire: (0..segments.len()).collect(),
            inner_wires: vec![],
            internal_wires: vec![],
        };
        Some((segments, vec![wf], vertex_positions))
    }
}

// =============================================================================
// Phase 2: OCCT 1:1 PerformLoops Alignment (BOPAlgo_BuilderFace.cxx L239-606)
// =============================================================================

/// Edge-like segment for wire building=can be a DS edge, an intersection curve,
impl<'a> BooleanBuilder<'a> {
    // OCCT BOPAlgo_Builder::BOPAlgo_Builder() (empty constructor)
    pub fn new(ds: &'a DS, op: BooleanOpType) -> Self {
        let context = RefCell::new(Context::new(ds.face_count(), TOLERANCE_ABS * 100.0));
        Self {
            ds,
            context,
            my_arguments: std::cell::RefCell::new(Vec::new()),
            my_map_fence: std::cell::RefCell::new(std::collections::HashSet::new()),
            my_entry_point: 0,
            my_images: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_shapes_sd: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_origins: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_in_parts: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_split_images_in: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_non_destructive: false,
            my_glue: GlueEnum::GlueOff,
            my_check_inverted: false,
            my_report: std::cell::RefCell::new(Report::new()),
            my_run_parallel: false,
            my_fuzzy_value: TOLERANCE_ABS,
            my_use_obb: false,
            my_shape: std::cell::RefCell::new(rcad_kernel::topods::BRep::new()),
            my_fill_history: true, // OCCT default (BOPAlgo_BuilderShape.hxx L122)
            my_operation: op,
            my_dims: std::cell::Cell::new([3, 3]),
            // --- rcad-specific ---
            brep: std::cell::RefCell::new(None),
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
    pub fn set_brep_with_mappings(
        &self,
        brep: rcad_kernel::topods::BRep,
        face_refs: Vec<rcad_kernel::topods::Shape>,
        ic_edge_map: Vec<Option<rcad_kernel::topods::Shape>>,
    ) {
        *self.brep.borrow_mut() = Some((brep, face_refs, ic_edge_map));
    }

    /// Create builder with a pre-built BRep (A3 dual-write, skips ds_to_brep).
    pub fn with_brep(
        ds: &'a DS,
        op: BooleanOpType,
        brep: rcad_kernel::topods::BRep,
        face_refs: Vec<rcad_kernel::topods::Shape>,
        ic_edge_map: Vec<Option<rcad_kernel::topods::Shape>>,
    ) -> Self {
        let builder = Self::new(ds, op);
        builder.set_brep_with_mappings(brep, face_refs, ic_edge_map);
        builder
    }

    /// Shape backed by shared Arc in ds.shapes (myDS->Shape(n) identity).
    /// Falls back to synthetic for out-of-range indices (shell/solid sentinel keys).
    fn brep_sr(&self, flat_idx: usize) -> rcad_kernel::topods::Shape {
        if flat_idx < self.ds.shapes.len() {
            return rcad_kernel::topods::Shape {
                data: self.ds.shapes[flat_idx].clone(),
                index: flat_idx,
                orientation: rcad_kernel::topods::Orientation::Forward,
                location: 0,
            };
        }
        { rcad_kernel::topods::Shape::synthetic(flat_idx, rcad_kernel::topods::Orientation::Forward) }
    }

    pub fn with_glue(mut self, enable: bool, tolerance: f64) -> Self {
        self.my_glue = if enable {
            GlueEnum::GlueFull
        } else {
            GlueEnum::GlueOff
        };
        self.my_fuzzy_value = tolerance.max(TOLERANCE_ABS);
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
            BooleanOpType::Union => class == Classification::Out || class == Classification::On,
            BooleanOpType::Intersection => {
                class == Classification::In || class == Classification::On
            }
            BooleanOpType::Difference => match source {
                SourceSide::A => class == Classification::Out,
                SourceSide::B => class == Classification::In,
            },
        }
    }
}

impl<'a> BooleanBuilder<'a> {
    /// The top-level pipeline entry: dimension-by-dimension image filling
    /// (V= = = ACE= HELL= OLID), followed by BuildResult for each type.
    /// OCCT BOPAlgo_Builder::CheckData (Builder.cxx L132-142).
    /// Validates operation type and non-empty arguments.
    /// OCCT form: AddError on each failure, then HasErrors check at the end.
    fn check_data(&self) -> Result<(), BooleanError> {
        // OCCT L132-137: if (myArguments.Extent() < 2) -> AlertTooFewArguments
        let nb_args = self.my_arguments.borrow().len();
        if nb_args < 2 {
            self.my_report.borrow_mut()
                .add_alert(crate::bopalgo::Alert::TooFewArguments);
        }
        // OCCT: BOPAlgo_BuilderShape::CheckData() — base class checks.
        // rcad: operation type validation (BOPAlgo_BOP::CheckData adds this).
        match self.my_operation {
            BooleanOpType::Union | BooleanOpType::Intersection | BooleanOpType::Difference => {}
            _ => self.my_report
                .borrow_mut()
                .add_alert(crate::bopalgo::Alert::BOPNotSet),
        }
        // OCCT: if (HasErrors()) return;
        if self.has_errors() {
            return Err(BooleanError::InvalidOperation);
        }
        Ok(())
    }

    /// OCCT BOPAlgo_Builder::CheckFiller (Builder.cxx L144-152).
    /// Checks if the PaveFiller has been set and merges its report.
    fn check_filler(&self) {
        // OCCT L146-149: if (!myPaveFiller) -> AlertNoFiller
        // rcad: DS must be populated by PaveFiller (vertices exist = filler ran).
        // OCCT L151: GetReport()->Merge(myPaveFiller->GetReport());
        if self.ds.vertex_count() == 0 {
            self.my_report.borrow_mut()
                .add_alert(crate::bopalgo::Alert::NoFiller);
        }
    }

    ///  ?Prepare (BOPAlgo_Builder.cxx L156-164).
    /// OCCT: BRep_Builder.MakeCompound(myShape)  ?empty compound as result.
    /// rcad: initializes my_shape + returns (BRep, ResultBuilder) for downstream.
    fn prepare(&self) -> (topods::BRep, ResultBuilder) {
        *self.my_shape.borrow_mut() = topods::BRep::new();
        (topods::BRep::new(), ResultBuilder::new())
    }

    ///  ?create TShapes for all DS source shapes in my_shape.
    /// Equivalent to OCCT's myArguments populated with all source TopoDS_Shape.
    ///  ?TreatEmptyShape (BOPAlgo_BOP.cxx L214-319).
    /// Handles the case where one or both operands have no geometry.
    /// Returns Ok(Some(brep)) if a quick result was determined,
    /// Ok(None) if the full pipeline must run.
    fn treat_empty_shape(
        &self,
        a_faces: &[usize],
        b_faces: &[usize],
    ) -> Result<Option<topods::BRep>, BooleanError> {
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
        match self.my_operation {
            BooleanOpType::Union => {
                // OCCT L270-279: return non-empty side
                let src = if has_a {
                    ShapeOrigin::ShapeA
                } else {
                    ShapeOrigin::ShapeB
                };
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
                    let brep =
                        self.brep_of_side_topods(ShapeOrigin::ShapeA, a_faces.len(), b_faces.len());
                    Ok(Some(brep))
                }
            }
            _ => {
                // Unknown operation → fall through to full pipeline
                Ok(None)
            }
        }
    }

    /// BOPAlgo_BOP::PerformInternal1 (BOP.cxx L422-579).
    /// Every statement in OCCT L422-579 has a corresponding rcad line below.
    /// See comments for exact OCCT line references.
    /// Structural difference: L425-429 setup done in constructor, re-affirmed here.
    /// L531 BuildResult(SOLID) writes to t_brep, then L900 BuildRC filters and
    /// clears solids from t_brep (non-Union) =equivalent to OCCT removing from myShape.
    pub fn build_with_history(&mut self) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
        self.build_with_history_topods()
    }

    /// Same as build_with_history but returns topods::BRep directly ().
    pub fn build_with_history_topods(
        &mut self,
    ) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
        // SKIP: OCCT L425-429 copies (myPaveFiller, myDS, myContext, myFuzzyValue, myNonDestructive)
        // from theFiller argument.  rcad builds with_brep() which sets these in the constructor --
        // no re-assignment at the start of build_with_history_topods is needed.

        // OCCT L431-436: CheckData =validates arguments and merges PaveFiller report.
        // Populate my_arguments from DS source shapes (OCCT: SetArguments).
        let mut args = self.my_arguments.borrow_mut();
        args.clear();
        args.push(rcad_kernel::topods::Shape::synthetic(0, rcad_kernel::topods::Orientation::Forward));
        args.push(rcad_kernel::topods::Shape::synthetic(1, rcad_kernel::topods::Orientation::Forward));
        drop(args);
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
        self.check_data()?;

        // OCCT BOPAlgo_BOP::CheckData L138: CheckFiller
        self.check_filler();
        if self.has_errors() {
            return Err(BooleanError::InvalidOperation);
        }

        // OCCT L438-443: Prepare
        // rcad: prepare() initializes my_shape + returns ResultBuilder.
        let mut result = self.prepare().1;
        // Compute myDims from argument shape types.
        // OCCT: BRep_Tool::Dimension(myArguments.First/Last) → VERTEX=0, EDGE=1, WIRE/FACE=2, SHELL/SOLID=3.
        // rcad: infer from per-origin shape counts in the DS.
        let a_n_verts = self.ds.a_vertex_count();
        let a_n_edges = self.ds.a_edge_count();
        let a_n_faces = self.ds.a_face_count();
        let b_n_verts = self.ds.vertex_count() - a_n_verts;
        let b_n_edges = self.ds.edge_count() - a_n_edges;
        let b_n_faces = self.ds.face_count() - a_n_faces;
        let dim_a = if a_n_faces > 0 {
            3_i8
        } else if a_n_edges > 0 {
            1
        } else if a_n_verts > 0 {
            0
        } else {
            3
        };
        let dim_b = if b_n_faces > 0 {
            3_i8
        } else if b_n_edges > 0 {
            1
        } else if b_n_verts > 0 {
            0
        } else {
            3
        };
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
                    let side = if !a_faces.is_empty() {
                        ShapeOrigin::ShapeA
                    } else {
                        ShapeOrigin::ShapeB
                    };
                    self.source_history_single_side(side)
                } else {
                    vec![]
                };
                let mut history = result.build_topods(
                    &mut *self.my_shape.borrow_mut(),
                    self.my_fill_history,
                    &self.my_shells.borrow(),
                    &mut *self.my_face_refs.borrow_mut(),
                    &self.my_solids.borrow(),
                    &self.my_compsolid_groups.borrow(),
                );
                history.source_history = source_history;
                let result_brep = self.my_shape.borrow().clone();
                return Ok((result_brep, history));
            }
        }

        // OCCT L454-457: ProgressScope + PISteps + analyzeProgress.
        // rcad: progress reporting not yet integrated.
        // OCCT L459-471: 3.1 FillImagesVertices + BuildResult(VERTEX)
        self.fill_images_vertices();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Vertex, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_FillImagesVertices",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L472-483: 3.2 FillImagesEdges + BuildResult(EDGE)
        self.fill_images_edges();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Edge, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_FillImagesEdges",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L484-496: 3.3 FillImagesContainers(WIRE) + BuildResult(WIRE)
        self.fill_images_container(ShapeType::Wire, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Wire, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_BuildResultWire",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L497-509: 3.4 FillImagesFaces + BuildResult(FACE)
        // Architecture A1: split faces create TShapes incrementally during fill_images_faces.
        // Remaining unsplit faces have existing TShapes from pre-create_source_shapes.
        self.fill_images_faces();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        // BuildResult(FACE) — generic loop over my_arguments, adds originals/splits to result.
        self.build_result(topods::ShapeType::Face, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_FillImagesFaces",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L510-522: 3.5 FillImagesContainers(SHELL) + BuildResult(SHELL)
        self.fill_images_container(ShapeType::Shell, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Shell, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_BuildResultShell",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L523-535: 3.6 FillImagesSolids + BuildResult(SOLID)
        self.fill_images_solids(&mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Solid, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_FillImagesSolids",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L536-548: 3.7 FillImagesContainers(COMPSOLID) + BuildResult(COMPSOLID)
        self.fill_images_container(ShapeType::CompSolid, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::CompSolid, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_BuildResultCompSolid",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L549-561: 3.8 FillImagesCompounds + BuildResult(COMPOUND)
        self.fill_images_compounds(&mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Compound, &mut result);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        dump_ctx.snapshot(
            "after_FillImagesCompounds",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        // OCCT L563-568: 4. BuildShape (BOPAlgo_BOP.cxx L885+)
        // OCCT L563-568: 4. BuildShape (BOPAlgo_BOP.cxx L885+)
        // BuildRC: filter by building element membership (BOPAlgo_BOP.cxx L597-881).
        if self.my_operation != BooleanOpType::Union {
            let mut t = self.my_shape.borrow_mut();
            self.build_rc(&mut result, &mut *t);
        }
        // OCCT L569-573: 5. PrepareHistory — builds BRep TShapes + source shape history.
        let mut history = {
            let mut t_brep = self.my_shape.borrow_mut();
            let mut history = result.build_topods(
                &mut *t_brep,
                self.my_fill_history,
                &self.my_shells.borrow(),
                &mut *self.my_face_refs.borrow_mut(),
                &self.my_solids.borrow(),
                &self.my_compsolid_groups.borrow(),
            );
            let source_history = if self.my_fill_history {
                self.prepare_history(&mut *t_brep)
            } else {
                vec![]
            };
            history.source_history = source_history;
            history
        };
        // OCCT L577-578: 5. PostTreat
        // Corrects tolerances of the result shape (CorrectTolerances + CorrectShapeTolerances).
        dump_ctx.snapshot(
            "after_PrepareHistory",
            self.ds,
            Some(&*self.my_shape.borrow()),
        );
        self.post_treat();
        dump_ctx.snapshot("after_PostTreat", self.ds, Some(&*self.my_shape.borrow()));
        let result_brep = self.my_shape.borrow().clone();

        Ok((result_brep, history))
    }
}

/// Stage snapshot: counts of DS and BRep entities at a pipeline boundary.
#[derive(Debug, Clone)]
pub(crate) struct StageSnapshot {
    pub stage: u32,
    pub stage_name: &'static str,
    pub n_ds_vertices: usize,
    pub n_ds_edges: usize,
    pub n_ds_faces: usize,
    pub n_ds_pave_blocks: usize,
    pub n_ds_intersection_curves: usize,
    pub n_ds_interf_ff: usize,
    pub n_brep_vertices: usize,
    pub n_brep_edges: usize,
    pub n_brep_faces: usize,
    pub n_brep_shells: usize,
    pub n_brep_solids: usize,
}

/// Run pipeline stage by stage, collecting a snapshot after each stage.
/// Returns the final result + Vec of per-stage snapshots for test analysis.
/// Stops early and returns Ok((partial_result, snapshots)) if any stage
/// triggers `has_errors`, so the caller can inspect the failure point.
impl<'a> BooleanBuilder<'a> {
    pub(crate) fn build_with_history_stage_by_stage(
        &mut self,
    ) -> Result<(topods::BRep, BooleanHistory, Vec<StageSnapshot>), BooleanError> {
        let mut snapshots: Vec<StageSnapshot> = Vec::with_capacity(12);

        // Helper macro: capture current DS + BRep state into a snapshot.
        // Uses a macro instead of closure to avoid borrow conflicts.
        macro_rules! snap {
            ($stage:expr, $name:expr) => {{
                let b = self.my_shape.borrow();
                let (nv, ne, nf, nsh, nso) = count_brep_entities(&b);
                snapshots.push(StageSnapshot {
                    stage: $stage,
                    stage_name: $name,
                    n_ds_vertices: self.ds.vertex_count(),
                    n_ds_edges: self.ds.edge_count(),
                    n_ds_faces: self.ds.face_count(),
                    n_ds_pave_blocks: self.ds.pave_blocks.len(),
                    n_ds_intersection_curves: self.ds.intersection_curves.len(),
                    n_ds_interf_ff: self.ds.interf_ff.len(),
                    n_brep_vertices: nv,
                    n_brep_edges: ne,
                    n_brep_faces: nf,
                    n_brep_shells: nsh,
                    n_brep_solids: nso,
                });
            }};
        }

        // OCCT L431-436: CheckData
        let mut args = self.my_arguments.borrow_mut();
        args.clear();
        args.push(rcad_kernel::topods::Shape::synthetic(0, rcad_kernel::topods::Orientation::Forward));
        args.push(rcad_kernel::topods::Shape::synthetic(1, rcad_kernel::topods::Orientation::Forward));
        drop(args);
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
        self.check_data()?;
        self.check_filler();
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }

        // OCCT L438-443: Prepare
        let mut result = self.prepare().1;
        let a_n_verts = self.ds.a_vertex_count();
        let a_n_edges = self.ds.a_edge_count();
        let a_n_faces = self.ds.a_face_count();
        let b_n_verts = self.ds.vertex_count() - a_n_verts;
        let b_n_edges = self.ds.edge_count() - a_n_edges;
        let b_n_faces = self.ds.face_count() - a_n_faces;
        let dim_a = if a_n_faces > 0 {
            3_i8
        } else if a_n_edges > 0 {
            1
        } else if a_n_verts > 0 {
            0
        } else {
            3
        };
        let dim_b = if b_n_faces > 0 {
            3_i8
        } else if b_n_edges > 0 {
            1
        } else if b_n_verts > 0 {
            0
        } else {
            3
        };
        self.my_dims.set([dim_a, dim_b]);

        snap!(1, "Prepare");

        // TreatEmptyShape
        if a_faces.is_empty() || b_faces.is_empty() {
            if self.treat_empty_shape(&a_faces, &b_faces)?.is_some() {
                let source_history = if self.my_fill_history {
                    let side = if !a_faces.is_empty() {
                        ShapeOrigin::ShapeA
                    } else {
                        ShapeOrigin::ShapeB
                    };
                    self.source_history_single_side(side)
                } else {
                    vec![]
                };
                let mut history = result.build_topods(
                    &mut *self.my_shape.borrow_mut(),
                    self.my_fill_history,
                    &self.my_shells.borrow(),
                    &mut *self.my_face_refs.borrow_mut(),
                    &self.my_solids.borrow(),
                    &self.my_compsolid_groups.borrow(),
                );
                history.source_history = source_history;
                let result_brep = self.my_shape.borrow().clone();
                snapshots.push(StageSnapshot {
                    stage: 1,
                    stage_name: "TreatEmptyShape_early_return",
                    n_ds_vertices: self.ds.vertex_count(),
                    n_ds_edges: self.ds.edge_count(),
                    n_ds_faces: self.ds.face_count(),
                    n_ds_pave_blocks: self.ds.pave_blocks.len(),
                    n_ds_intersection_curves: self.ds.intersection_curves.len(),
                    n_ds_interf_ff: self.ds.interf_ff.len(),
                    n_brep_vertices: result_brep.tshapes
                        .iter()
                        .filter(|ts| matches!(&***ts, topods::TShape::Vertex(_)))
                        .count(),
                    n_brep_edges: result_brep.tshapes
                        .iter()
                        .filter(|ts| matches!(&***ts, topods::TShape::Edge(_)))
                        .count(),
                    n_brep_faces: result_brep.tshapes
                        .iter()
                        .filter(|ts| matches!(&***ts, topods::TShape::Face(_)))
                        .count(),
                    n_brep_shells: result_brep.tshapes
                        .iter()
                        .filter(|ts| matches!(&***ts, topods::TShape::Shell(_)))
                        .count(),
                    n_brep_solids: result_brep.tshapes
                        .iter()
                        .filter(|ts| matches!(&***ts, topods::TShape::Solid(_)))
                        .count(),
                });
                return Ok((result_brep, history, snapshots));
            }
        }

        // Stage 1: FillImagesVertices + BuildResult(Vertex)
        snap!(2, "before_FillImagesVertices");
        self.fill_images_vertices();
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::Vertex, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(3, "after_FillImagesVertices");

        // Stage 2: FillImagesEdges + BuildResult(Edge)
        snap!(4, "before_FillImagesEdges");
        self.fill_images_edges();
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::Edge, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(5, "after_FillImagesEdges");

        // Stage 3: FillImagesContainers(WIRE) + BuildResult(WIRE)
        snap!(6, "before_FillImagesContainers_Wire");
        self.fill_images_container(ShapeType::Wire, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::Wire, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(7, "after_BuildResultWire");

        // Stage 4: FillImagesFaces + BuildResult(FACE)
        snap!(8, "before_FillImagesFaces");
        self.fill_images_faces();
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::Face, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(9, "after_FillImagesFaces");

        // Stage 5: FillImagesContainers(SHELL) + BuildResult(SHELL)
        snap!(10, "before_FillImagesContainers_Shell");
        self.fill_images_container(ShapeType::Shell, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::Shell, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(11, "after_BuildResultShell");

        // Stage 6: FillImagesSolids + BuildResult(SOLID)
        snap!(12, "before_FillImagesSolids");
        self.fill_images_solids(&mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::Solid, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(13, "after_FillImagesSolids");

        // Stage 7: FillImagesContainers(COMPSOLID) + BuildResult(COMPSOLID)
        snap!(14, "before_FillImagesContainers_CompSolid");
        self.fill_images_container(ShapeType::CompSolid, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::CompSolid, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(15, "after_BuildResultCompSolid");

        // Stage 8: FillImagesCompounds + BuildResult(COMPOUND)
        snap!(16, "before_FillImagesCompounds");
        self.fill_images_compounds(&mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        self.build_result(topods::ShapeType::Compound, &mut result);
        if self.has_errors() {
            return Ok((
                self.my_shape.borrow().clone(),
                BooleanHistory::default(),
                snapshots,
            ));
        }
        snap!(17, "after_FillImagesCompounds");

        // PrepareHistory
        let mut history = {
            let mut t_brep = self.my_shape.borrow_mut();
            let mut history = result.build_topods(
                &mut *t_brep,
                self.my_fill_history,
                &self.my_shells.borrow(),
                &mut *self.my_face_refs.borrow_mut(),
                &self.my_solids.borrow(),
                &self.my_compsolid_groups.borrow(),
            );
            let source_history = if self.my_fill_history {
                self.prepare_history(&mut *t_brep)
            } else {
                vec![]
            };
            history.source_history = source_history;
            history
        };
        snap!(18, "after_PrepareHistory");

        // PostTreat
        self.post_treat();
        let result_brep = self.my_shape.borrow().clone();
        snap!(19, "after_PostTreat");

        Ok((result_brep, history, snapshots))
    }
}

/// Count V/E/F/Shell/Solid TShapes in a BRep (owned, no borrow).
fn count_brep_entities(b: &topods::BRep) -> (usize, usize, usize, usize, usize) {
    let mut nv = 0;
    let mut ne = 0;
    let mut nf = 0;
    let mut nsh = 0;
    let mut nso = 0;
    for ts in &b.tshapes {
        match &**ts {
            topods::TShape::Vertex(_) => nv += 1,
            topods::TShape::Edge(_) => ne += 1,
            topods::TShape::Face(_) => nf += 1,
            topods::TShape::Shell(_) => nsh += 1,
            topods::TShape::Solid(_) => nso += 1,
            _ => {}
        }
    }
    (nv, ne, nf, nsh, nso)
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
impl<'a> BooleanBuilder<'a> {
    fn faces_of(&self, origin: ShapeOrigin) -> Vec<usize> {
        let mut v: Vec<usize> = self.ds
            .faces.iter()
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
        let mut v_map: HashMap<usize, topods::Shape> = HashMap::new();
        let mut e_map: HashMap<usize, topods::Shape> = HashMap::new();

        for fi in 0..self.ds.face_count() {
            if self.ds.face_origin(fi) != origin {
                continue;
            }
            let mut edge_refs: Vec<topods::Shape> = Vec::new();
            for &ei in self.ds.face_boundary_edges(fi) {
                if ei >= self.ds.edge_count() {
                    continue;
                }
                let sv_idx = self.ds.edge_start_vertex_ds(ei);
                let ev_idx = self.ds.edge_end_vertex_ds(ei);
                let sv = v_map.entry(sv_idx)
                    .or_insert_with(|| t.add_tvertex(self.ds.vertex_point(sv_idx)))
                    .clone();
                let ev = v_map.entry(ev_idx)
                    .or_insert_with(|| t.add_tvertex(self.ds.vertex_point(ev_idx)))
                    .clone();
                let edge_sr = e_map.entry(ei)
                    .or_insert_with(|| t.add_tedge(None, sv, ev, [0.0, 1.0]))
                    .clone();
                edge_refs.push(edge_sr);
            }
            if edge_refs.is_empty() {
                continue;
            }
            let ow = t.add_twire(edge_refs);
            let surf = self.ds.face_surface(fi).cloned().unwrap_or_else(|| {
                Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z))
            });
            let face = t.add_tface(Some(surf), ow, vec![], None, None, vec![], false);
            t.face_mut(face).tolerance = self.ds.face_tolerance(fi);
        }

        // Collect face refs and wrap in Shell → Solid
        let face_srs: Vec<topods::Shape> = t.tshapes
            .iter()
            .enumerate()
            .filter(|(_, ts)| matches!(&***ts, topods::TShape::Face(_)))
            .map(|(i, _)| {
                let idx = i;
                topods::Shape::synthetic(idx, topods::Orientation::Forward)
            })
            .collect();
        if !face_srs.is_empty() {
            let shell = t.add_tshell(face_srs);
            t.add_tsolid(vec![shell]);
        }
        t
    }
}
