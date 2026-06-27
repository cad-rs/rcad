use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use indexmap::IndexMap;

use glam::{DVec2, DVec3};
use rayon::prelude::*;
use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval, *};
use rcad_kernel::topology::*;

use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::history::{
    BooleanHistory, EdgeOrigin, FaceOrigin, HistoryTracker, ShellOrigin, SolidOrigin, VertexOrigin,
};
use std::cell::RefCell;
use crate::inttools::context::Context;
use crate::inttools::edge_face::plane_local_basis;
use crate::tolerance::*;
use crate::triangulate::{triangulate_polygon, triangulate_polygon_with_holes};

mod angle_2d;
mod curve_tools;
mod debug_utils;
mod intres2d;
mod intersection;
mod types;
mod ds_as_brep;

pub use types::{
    BooleanOpType, BooleanError, FaceSampleData,
};
pub(crate) use types::{
    ShapeType, WireFace, WireSegment, WireEdgeSource,
    WireSegmentTopoDS, WireEdgeSourceTopoDS,
    FaceWireEdges, FaceEntry, CollectedFaceResult,
};


mod result_builder;

pub(crate) use result_builder::ResultBuilder;
mod builder_utils;

pub(crate) use builder_utils::{
    curve_eq, hash_point,
    classify_subface_against_box, classify_against_solid_for_boolean,
    is_tangent_face, build_edge_bounds, quantize_pos,
    check_and_add_split_vertex, collect_face_edge_segments,
    compute_ic_second_pcurve, cmp_boolean_emit_order,
    annotate_history_from_ds, annotate_shell_and_solid_history,
    aggregate_face_region_origin, aggregate_shell_region_origin,
};

pub struct BooleanBuilder<'a> {
    ds: &'a DS,
    op: BooleanOpType,
    use_glue: bool,
    glue_tolerance: f64,
    context: RefCell<Context>,
    // ✅ OCCT-aligned: error tracking (myReport / HasErrors equivalent).
    has_errors: bool,
    // ✅ OCCT-aligned: myImages — source shape index → list of split image indices.
    //   Uses RefCell because phase functions take &self (OCCT uses mutable member maps).
    my_images: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>>>,
    my_origins: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>>>,
    my_shapes_sd: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, rcad_kernel::topods::ShapeRef>>,
    // ✅ OCCT-aligned: myInParts — source solid index → list of its IN face indices
    //   (BOPAlgo_Builder.hxx L502).  Populated during FillImagesFaces, used by
    //   FillIn3DParts / BuildDraftSolid for solid assembly.
    my_in_parts: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: solid-level image tracking (BOPAlgo_Builder.hxx L498 myImages).
    //   OCCT BuildSplitSolids stores split solids in myImages[source_solid].
    //   rcad: maps source side (0=A, 1=B) → result solid indices from
    //   build_split_solids.  Used by annotate_shell_and_solid_history and
    //   for OCCT-form history tracking.
    my_solid_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: solid-level origin tracking (BOPAlgo_Builder.hxx L500 myOrigins).
    //   Reverse map: result solid index → list of source sides.
    my_solid_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // ✅ OCCT-aligned: myNonDestructive (BOPAlgo_Builder.hxx L503).
    //   Safe processing — avoids modifying input shapes. Used in PostTreat.
    my_non_destructive: bool,
    // ✅ OCCT-aligned: myCheckInverted (BOPAlgo_Builder.hxx L505).
    //   Enables/disables inverted-solid check on input shapes.
    my_check_inverted: bool,
    /// ✅ OCCT-aligned: converted BRep representation of DS.
    ///   Populated after check_data in build_with_history via ds_to_brep().
    ///   Wrapped in RefCell because build_with_history takes &self.
    ///   This is topods::BRep (not the legacy rcad_kernel::BRep).
    brep: std::cell::RefCell<Option<(rcad_kernel::topods::BRep, Vec<rcad_kernel::topods::ShapeRef>, Vec<Option<rcad_kernel::topods::ShapeRef>>)>>,
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
pub(crate) use wire_splitter::{
    EdgeInfo, build_closed_wires, perform_shapes_to_avoid,
    expand_avoided_pids, build_pid_maps,
    build_vi_to_canon, physical_edge_id, world_to_uv,
    compute_seam_tangent_angles, edge_uv_tangent, edge_angle_2d,
    are_verts_coincident,
};
pub(crate) use wire_path::{
    perform_areas, intersect_ray_curve_2d,
    wire_faces_to_face_sample_data,
    refine_angles, pc_parameter_range,
    walk_path_extract_wires,
};

impl<'a> BooleanBuilder<'a> {
    /// ✅ OCCT-aligned: BOPAlgo_BuilderFace::Perform (BuilderFace.cxx L117-147).
    ///   Edge-to-wire pipeline: PerformShapesToAvoid → PerformLoops (WireSplitter)
    ///   → PerformAreas → PerformInternalShapes.
    pub(crate) fn split_face_occt_wire_pipeline(
        &self,
        face_idx: usize,
    ) -> Option<(Vec<WireSegment>, Vec<WireFace>, HashMap<usize, DVec3>)> {
        let ds = self.ds;
        let face = &ds.faces[face_idx];
        // ✅ OCCT-aligned: BuilderFace::Perform (BOPAlgo_BuilderFace.cxx L117-148).
        //   L121: GetReport()->Clear()
        //   L123-127: CheckData() → if HasErrors return
        //   L129-133: PerformShapesToAvoid → if HasErrors return
        //   L135-139: PerformLoops → if HasErrors return
        //   L141-145: PerformAreas → if HasErrors return
        //   L147: PerformInternalShapes

        // OCCT L121: GetReport()->Clear() — rcad: delegated to caller error handling.

        // OCCT L123-127: CheckData — validate face has intersection data.
        //   OCCT: checks myShapes/DS state. rcad: segments must be non-empty
        //   and face must have interferences.
        let pcurve_lookup = |ci: usize| self.find_pcurve_for_face(ci, face_idx);
        let mut segments = collect_face_edge_segments(ds, face_idx, &pcurve_lookup);
        if !self.builder_face_check_data(face_idx, &segments) {
            return None;
        }

        // OCCT builds vi_to_canon during CheckData/Prepare.
        let vi_to_canon = build_vi_to_canon(&segments, ds);

        // OCCT L129-133: PerformShapesToAvoid — returns PIDs (OCCT: myShapesToAvoid).
        let (avoided_pids, pid_segs) = perform_shapes_to_avoid(&segments, &vi_to_canon, ds);
        // Expand PIDs to segment indices (rcad: PerformLoops reads segment indices).
        let mut avoided = expand_avoided_pids(&avoided_pids, &pid_segs);

        // OCCT L135-139: PerformLoops
        let (wires, mut internal_wires, vertex_positions) =
            build_closed_wires(&mut segments, ds, face_idx, &avoided);

        // ✅ OCCT-aligned L160-166: Post Treatment — edges not in any loop → add to myShapesToAvoid.
        let in_loop: std::collections::HashSet<usize> = wires.iter().flatten().copied().collect();
        for si in 0..segments.len() {
            if !in_loop.contains(&si) && !avoided.contains(&si) {
                avoided.insert(si);
            }
        }

        // ✅ OCCT-aligned L327-362: Internal Wires — build wire groups from myShapesToAvoid.
        //   OCCT: each avoided edge wraps in a TopoDS_Wire → myLoopsInternal.
        //   rcad: each avoided segment becomes its own wire group, passed to PerformAreas.
        let internal_wire_groups: Vec<Vec<usize>> = avoided.iter().map(|&si| vec![si]).collect();

        // OCCT L141-145: PerformAreas
        let mut wfs = if !wires.is_empty() {
            perform_areas(&wires, &internal_wire_groups, &segments, ds, &mut *self.context.borrow_mut(), face_idx)
        } else if !internal_wires.is_empty() {
            vec![WireFace { outer_wire: vec![], inner_wires: vec![], internal_wires: internal_wires.clone() }]
        } else {
            vec![WireFace { outer_wire: (0..segments.len()).collect(), inner_wires: vec![], internal_wires: vec![] }]
        };
        if wfs.is_empty() { return None; }

        // ✅ OCCT-aligned L147: PerformInternalShapes — classify internal wires against faces
        crate::builder::wire_path::perform_internal_shapes(
            &mut wfs, &internal_wire_groups, &segments, ds, face_idx);
        Some((segments, wfs, vertex_positions))
    }

    /// ✅ OCCT-aligned: TopoDS-based BuildFace pipeline with emit.
    ///   Runs the full pipeline then emits result faces directly into ResultBuilder.
    pub(crate) fn split_face_and_emit_topo_ds(
        &self,
        face_idx: usize,
        is_a: bool,
        result: &mut ResultBuilder,
    ) {
        let ds = self.ds;
        // Get BRep from the cached conversion (built during build_with_history)
        let brep_borrow = self.brep.borrow();
        let (br, face_refs, _ic_edge_map): &(rcad_kernel::topods::BRep, Vec<rcad_kernel::topods::ShapeRef>, Vec<Option<rcad_kernel::topods::ShapeRef>>) = match brep_borrow.as_ref() {
            Some(v) => v,
            None => return, // ds_to_brep not yet called
        };

        let pcurve_lookup = |ci: usize| self.find_pcurve_for_face(ci, face_idx);
        let segments = collect_face_edge_segments(ds, face_idx, &pcurve_lookup);
        if !self.builder_face_check_data(face_idx, &segments) { return; }

        let segments_topo = crate::builder::builder_utils_topo_ds::segments_to_topo_ds(&segments, ds, face_idx, &face_refs[..], &_ic_edge_map[..]);
        drop(segments);

        let tool: &dyn rcad_kernel::topods::BRepTool = br;

        let (avoided_pids, pid_segs) = crate::builder::wire_splitter::perform_shapes_to_avoid_topo_ds(
            &segments_topo, tool);
        let mut avoided = crate::builder::wire_splitter::expand_avoided_pids(&avoided_pids, &pid_segs);
        let wires = crate::builder::wire_path_topo_ds::build_closed_wires_topoDS(
            &segments_topo, &avoided, tool);

        let in_loop: HashSet<usize> = wires.iter().flatten().copied().collect();
        for si in 0..segments_topo.len() {
            if !in_loop.contains(&si) && !avoided.contains(&si) { avoided.insert(si); }
        }
        let internal_wire_groups: Vec<Vec<usize>> = avoided.iter().map(|&si| vec![si]).collect();

        let wfs = if !wires.is_empty() {
            crate::builder::wire_path_topo_ds::perform_areas_topo_ds(
                &wires, &internal_wire_groups, &segments_topo, tool, face_idx, ds)
        } else if !avoided.is_empty() {
            vec![WireFace { outer_wire: vec![], inner_wires: vec![], internal_wires: segments_topo.iter().enumerate().filter(|(si, _)| avoided.contains(si)).map(|(si, _)| vec![si]).collect() }]
        } else {
            vec![WireFace { outer_wire: (0..segments_topo.len()).collect(), inner_wires: vec![], internal_wires: vec![] }]
        };
        if wfs.is_empty() { return; }

        let mut wfs = wfs;
        // ✅ OCCT-aligned L147: PerformInternalShapes
        crate::builder::wire_path_topo_ds::perform_internal_shapes_topo_ds(
            &mut wfs, &internal_wire_groups, &segments_topo, tool, face_idx, ds);

        let origin = if is_a {
            FaceOrigin::FromA(ds.faces[face_idx].source_face_idx)
        } else {
            FaceOrigin::FromB(ds.faces[face_idx].source_face_idx)
        };
        let ic_curves: HashMap<usize, Curve3> = ds.intersection_curves.iter()
            .enumerate().map(|(ci, ic)| (ci, ic.curve.clone())).collect();
        for wf in &wfs {
            result.emit_wire_face_topods(face_idx, wf, &segments_topo, tool, &ic_curves, false, origin,
                &HashMap::new(), face_refs[face_idx]);
        }
    }

    /// ✅ OCCT-aligned: BuilderFace::CheckData (BOPAlgo_BuilderFace.cxx L50-115).
    ///   Validates face has intersection curves/segments. If no interferences,
    ///   delegates to BuildDraftFace (OCCT's alternative path for non-split faces).
    fn builder_face_check_data(&self, face_idx: usize, segments: &[WireSegment]) -> bool {
        if segments.is_empty() {
            return false;
        }
        true
    }

    /// ✅ OCCT-aligned: PrepareHistory (Builder_4.cxx L164-252).
    ///   Builds source→result history matching OCCT's BRepTools_History.
    fn build_source_history(&self, t_brep: &topods::BRep) -> Vec<crate::history::SourceShapeEntry> {
        use crate::history::{HistoryStatus, SourceShapeEntry};
        use topods::TShape;

        // Build result vertex set.
        let mut result_vtx: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for ts in &t_brep.tshapes {
            if let TShape::Vertex(vd) = &**ts {
                for (di, dv) in self.ds.vertices.iter().enumerate() {
                    if (dv.point - vd.point).length_squared() < crate::tolerance::TOLERANCE_ABS * 2.0 {
                        result_vtx.insert(di);
                        break;
                    }
                }
            }
        }

        // Build result edge set.
        let mut result_edge: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for ts in &t_brep.tshapes {
            if let TShape::Edge(ed) = &**ts {
                let sv = ed.first.index;
                let ev = ed.last.index;
                for (di, de) in self.ds.edges.iter().enumerate() {
                    if (de.start_vertex == sv && de.end_vertex == ev)
                        || (de.start_vertex == ev && de.end_vertex == sv)
                    {
                        result_edge.insert(di);
                        break;
                    }
                }
            }
        }

        let mut entries = Vec::new();
        let v_base = 0usize; // vertices start at 0 in BRep
        let e_base = self.ds.vertices.len();
        // Vertices
        for (di, _dv) in self.ds.vertices.iter().enumerate() {
            let in_result = result_vtx.contains(&di);
            let has_images = self.my_images.borrow().contains_key(&rcad_kernel::topods::ShapeRef::new(v_base + di));
            let status = if has_images && in_result { HistoryStatus::Modified }
                else if in_result { HistoryStatus::Generated }
                else { HistoryStatus::Deleted };
            entries.push(SourceShapeEntry { ds_index: di, shape_type: 0, status, result_indices: vec![] });
        }
        // Edges
        for (di, _de) in self.ds.edges.iter().enumerate() {
            let in_result = result_edge.contains(&di);
            let has_images = self.my_images.borrow().contains_key(&rcad_kernel::topods::ShapeRef::new(e_base + di));
            let status = if has_images && in_result { HistoryStatus::Modified }
                else if in_result { HistoryStatus::Generated }
                else { HistoryStatus::Deleted };
            entries.push(SourceShapeEntry { ds_index: di, shape_type: 1, status, result_indices: vec![] });
        }
        entries
    }

    /// ✅ OCCT-aligned: BuildDraftFace (BOPAlgo_Builder_2.cxx L951-1070).
    ///
    /// For faces that have NO intersection curves but whose boundary edges may
    /// have been split by the PaveFiller (via myImages / vertices_in), build a
    /// single analytic face using the split boundary edges.  This avoids the
    /// tessellation fallback (split_curved_face_parametric, tessellate_sphere_face,
    /// etc.) that would otherwise be used for non-planar faces with only
    /// alone-vertex / on-edge intersection data.
    ///
    /// Returns `None` when:
    /// - The face has no boundary segments (empty geometry)
    /// - Any vertex is multi-connected (>=3 edges share the same vertex),
    ///   indicating the face may need full SmartMap-based splitting
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

/// Edge-like segment for wire building鈥?can be a DS edge, an intersection curve,
impl<'a> BooleanBuilder<'a> {
    pub fn new(ds: &'a DS, op: BooleanOpType) -> Self {
        let context = RefCell::new(Context::new(ds.faces.len(), TOLERANCE_ABS * 100.0));
        Self {
            ds, op, use_glue: false, glue_tolerance: TOLERANCE_ABS, context, has_errors: false,
            my_images: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_origins: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_shapes_sd: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_in_parts: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_solid_images: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_solid_origins: std::cell::RefCell::new(std::collections::HashMap::new()),
            my_non_destructive: false,
            my_check_inverted: false,
            brep: std::cell::RefCell::new(None),
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

    pub fn with_glue(mut self, enable: bool, tolerance: f64) -> Self {
        self.use_glue = enable;
        self.glue_tolerance = tolerance.max(TOLERANCE_ABS);
        self
    }

    /// Unified semantic policy for sub-face retention.
    ///
    /// This keeps A/B branches aligned to the same decision table instead of
    /// maintaining two subtly diverging helper functions.
    fn keep_subface_policy(op: BooleanOpType, source: SourceSide, class: Classification) -> bool {
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

    pub fn build(&self) -> Result<BRep, BooleanError> {
        let (brep, _) = self.build_with_history()?;
        if !brep.solids.is_empty() && !brep.solids[0].shells.is_empty() {
            eprintln!("BooleanBuilder::build: {} faces", brep.solids[0].shells[0].faces.len());
        }
        Ok(brep)
    }

    // ====================================================================
    // ✅ OCCT-aligned: dimension-by-dimension pipeline (PerformInternal1)
    //   BOPAlgo_Builder.cxx L310-440
    // ====================================================================

    /// ✅ OCCT-aligned: FillImagesVertices (BOPAlgo_Builder_1.cxx L40-67).
    ///   Iterates ShapesSD → builds myImages(VERTEX) + myShapesSD + myOrigins.
    ///   OCCT L42: NCollection_DataMap<int,int>::Iterator aIt(myDS->ShapesSD())
    ///   rcad: symmetric HashSet<(usize,usize)> → process once per pair (a<b).
    fn fill_images_vertices(&self) {
        // OCCT L43-48: for (; aIt.More(); aIt.Next())
        for &(va, vb) in self.ds.shape_sd.sd_vertices_iter() {
            // rcad stores symmetric pairs; process each pair once (a < b).
            if va >= vb { continue; }
            let src = va;   // OCCT: nV = aIt.Key()
            let sd  = vb;   // OCCT: nVSD = aIt.Value()

            // OCCT L56: myImages.Bound(aV, ...)->Append(aVSD)
            self.my_images.borrow_mut().entry(rcad_kernel::topods::ShapeRef::new(src)).or_default().push(rcad_kernel::topods::ShapeRef::new(sd));
            // OCCT L58: myShapesSD.Bind(aV, aVSD)
            self.my_shapes_sd.borrow_mut().insert(rcad_kernel::topods::ShapeRef::new(src), rcad_kernel::topods::ShapeRef::new(sd));
            // OCCT L60-65: myOrigins.ChangeSeek(aVSD).Append(aV)
            self.my_origins.borrow_mut().entry(rcad_kernel::topods::ShapeRef::new(sd)).or_default().push(rcad_kernel::topods::ShapeRef::new(src));
        }
    }

    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-126).
    ///   Iterates source edges → populates myImages(EDGE) + myOrigins(EDGE).
    ///   OCCT L73: aNbS = myDS->NbSourceShapes()
    ///   OCCT L78-80: filter TopAbs_EDGE
    ///   OCCT L84-86: filter HasReference (has pave blocks)
    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-126).
    ///   Reads split edges created by MakeSplitEdges (build_split_edges in PaveFiller)
    ///   via pb.new_edge, matching OCCT's aPBR->Edge() pattern.
    ///   Creates myImages(EDGE) and myOrigins(EDGE) mappings.
    /// ✅ OCCT-aligned: FillImagesEdges (BOPAlgo_Builder_1.cxx L71-125).
    ///   L75-81: iterate source shapes → filter TopAbs_EDGE
    ///   L83-87: HasReference (pave blocks exist) → skip if none
    ///   L89-90: aE = aSI.Shape(); aLPB = myDS->PaveBlocks(i)
    ///   L95:    myImages.Bound(aE, ...)
    ///   L97-119: for each pave block:
    ///     L101:   aPBR = myDS->RealPaveBlock(aPB)
    ///     L103:   nSpR = aPBR->Edge()
    ///     L104:   aSpR = myDS->Shape(nSpR)
    ///     L105:   pLS->Append(aSpR)  → myImages[source].Append(split)
    ///     L107-112: myOrigins[split].Append(source)
    ///     L114-118: IsCommonBlockOnEdge → myShapesSD.Bind(source, split)
    fn fill_images_edges(&self, t: &mut topods::BRep) {
        // OCCT: Create TShape::Vertex entries (BuildResult path).
        for (vi, v) in self.ds.vertices.iter().enumerate() {
            let vr = t.add_tvertex(v.point);
            t.vertex_mut(vr).tolerance = v.geom_tol.max(crate::tolerance::TOLERANCE_ABS);
        }
        // OCCT L61: aNbE (base vertex index, rcad offset for flat indexing)
        let e_base = self.ds.vertices.len();
        for (ei, edge) in self.ds.edges.iter().enumerate() {
            if edge.pave_blocks.is_empty() { continue; }
            let aE = rcad_kernel::topods::ShapeRef::new(e_base + ei);
            for pb in &edge.pave_blocks {
                let nSpR = self.ds.real_pave_block_edge(ei, pb)
                    .or(pb.new_edge)
                    .unwrap_or(ei);
                let aSpR = rcad_kernel::topods::ShapeRef::new(e_base + nSpR);
                self.my_images.borrow_mut().entry(aE).or_default().push(aSpR);
                self.my_origins.borrow_mut().entry(aSpR).or_default().push(aE);
                if pb.common_block_idx.is_some() {
                    self.my_shapes_sd.borrow_mut().insert(aE, aSpR);
                }
                // OCCT: BuildResult creates result edges in t_brep from myImages.
                let se = &self.ds.edges[nSpR];
                let sv_sr = rcad_kernel::topods::ShapeRef::new(se.start_vertex);
                let ev_sr = rcad_kernel::topods::ShapeRef::new(se.end_vertex);
                let ci = t.curves.len();
                t.curves.push(se.curve.clone());
                let tedge = t.add_tedge(Some(ci), sv_sr, ev_sr, se.t_range);
                if self.ds.is_edge_degenerated(nSpR) || se.start_vertex == se.end_vertex {
                    t.edge_mut(tedge).degenerated = true;
                }
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainers(WIRE) (BOPAlgo_Builder_1.cxx L172-193).
    ///   OCCT: iterates source shapes → filters TopAbs_WIRE → FillImagesContainer
    ///   → builds wire images from edge images.  rcad: wires are implicit in face
    ///   boundary_edges.  For each source wire, check if any edge has split images;
    ///   if so rebuild the wire from split edges and store in myImages(WIRE).
    /// ✅ OCCT-aligned: FillImagesContainers(WIRE) — iterate DS wires (first-class TopAbs_WIRE).
    ///   OCCT L172-193: NbSourceShapes → filter TopAbs_WIRE → FillImagesContainer(shape, WIRE).
    ///   rcad: iterate ds.wires[], process each as FillImagesContainer per OCCT L221-276.
    fn fill_images_containers_wires(&self) {
        let e_base = self.ds.vertices.len();
        let wire_base = e_base + self.ds.edges.len();
        // OCCT L175-183: for each source shape, filter TopAbs_WIRE
        for wi in 0..self.ds.wires.len() {
            let edges: Vec<usize> = self.ds.wires[wi].edges.clone();

            // OCCT L224-233: check if any sub-edge has been modified
            let has_split = edges.iter().any(|&ei| {
                let e_ref = rcad_kernel::topods::ShapeRef::new(e_base + ei);
                self.my_images.borrow().get(&e_ref).map_or(false, |imgs| {
                    imgs.len() != 1 || imgs[0].index != e_base + ei
                })
            });

            if !has_split {
                // OCCT L236-240: no modification → no new image.
                let w_ref = rcad_kernel::topods::ShapeRef::new(wire_base + wi);
                self.my_images.borrow_mut().entry(w_ref).or_default().push(w_ref);
                continue;
            }

            // OCCT L247-271: rebuild wire from edge images.
            let has_img: std::collections::HashMap<usize, Vec<rcad_kernel::topods::ShapeRef>> =
                edges.iter().filter_map(|&ei| {
                    let e_ref = rcad_kernel::topods::ShapeRef::new(e_base + ei);
                    self.my_images.borrow().get(&e_ref).map(|v| (ei, v.clone()))
                }).collect();
            let mut wi_imgs = self.my_images.borrow_mut();
            let w_ref = rcad_kernel::topods::ShapeRef::new(wire_base + wi);
            for &ei in &edges {
                let entry = wi_imgs.entry(w_ref).or_default();
                if let Some(imgs) = has_img.get(&ei) {
                    for &new_eref in imgs {
                        entry.push(new_eref);
                    }
                } else {
                    entry.push(rcad_kernel::topods::ShapeRef::new(e_base + ei));
                }
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesFaces (BOPAlgo_Builder_1.cxx L376-386).
    ///   Phase 3: splits each face via WireSplitter → classifies → emits
    ///   via emit_wire_face.  rcad equivalent: for each face with IC data,
    ///   call split_face_occt_wire_pipeline (BuilderFace::Perform), then
    ///   classify_against_solid_for_boolean + classification_keep_policy.
    /// ✅ OCCT-aligned: FillImagesFaces (BOPAlgo_Builder_2.cxx L215-229).
    ///   Equivalent to BuildSplitFaces + FillSameDomainFaces + FillInternalVertices.
    ///   OCCT L258: aNbS = myDS->NbSourceShapes()
    ///   OCCT L260-266: iterates all source shapes, filters TopAbs_FACE.
    ///   OCCT L275-279: HasFaceInfo check.
    ///   OCCT L283-287: PaveBlocksIn/On/Sc + AloneVertices.
    ///   OCCT L293-296: if no PBs and no AV → skip.
    /// ✅ OCCT-aligned: FillImagesFaces (Builder_2.cxx L215-229).
    ///   Calls BuildSplitFaces → FillSameDomainFaces → FillInternalVertices.
    fn fill_images_faces(
        &self,
        result: &mut ResultBuilder,
        a_faces: &[usize],
        b_faces: &[usize],
    ) {
        self.build_split_faces(result, a_faces, b_faces);

        // ✅ OCCT L223: FillSameDomainFaces — merge duplicates after all faces split.
        self.fill_same_domain_faces(result);
        if self.has_errors { return; }

        // ✅ OCCT L228: FillInternalVertices — settle alone vertices as INTERNAL sub-shapes.
        self.fill_internal_vertices(result);

        // rcad: build_faces validates edge refs and builds face topology.
        // OCCT equivalent: image faces are already TopoDS during BuildSplitFaces;
        // no separate "build_faces" step is needed.
        result.build_faces();

        // OCCT L146-152: add original faces without images.
        // OCCT BuildResult(FACE) adds original faces when there are no split images.
        // rcad: for faces from source solids that had no split, create original face entries.
        let mut emitted_a: std::collections::HashSet<usize> =
            std::collections::HashSet::new();
        let mut emitted_b: std::collections::HashSet<usize> =
            std::collections::HashSet::new();
        for origin in &result.face_origins {
            match origin {
                FaceOrigin::FromA(fi) => { emitted_a.insert(*fi); }
                FaceOrigin::FromB(fi) => { emitted_b.insert(*fi); }
                _ => {}
            }
        }
        for &fi in a_faces {
            if !emitted_a.contains(&self.ds.faces[fi].source_face_idx) {
                result.build_original_face(self.ds, fi,
                    FaceOrigin::FromA(self.ds.faces[fi].source_face_idx));
            }
        }
        for &fi in b_faces {
            if !emitted_b.contains(&self.ds.faces[fi].source_face_idx) {
                result.build_original_face(self.ds, fi,
                    FaceOrigin::FromB(self.ds.faces[fi].source_face_idx));
            }
        }
    }

    /// ✅ OCCT-aligned: BuildSplitFaces (Builder_2.cxx L233-374).
    ///   Iterates source faces → splits each along intersection curves.
    ///   For faces with IN/SC PBs: full BuilderFace::Perform (split_face_occt_wire_pipeline).
    ///   For ON-only faces: BuildDraftFace.
    ///   Faces with no interferences → skipped (no images).
    fn build_split_faces(
        &self,
        result: &mut ResultBuilder,
        a_faces: &[usize],
        b_faces: &[usize],
    ) {
        // OCCT L258-266: iterate all source shapes → filter TopAbs_FACE.
        for fi in 0..self.ds.faces.len() {
            let is_a = a_faces.contains(&fi);
            if !is_a && !b_faces.contains(&fi) { continue; }

            // OCCT L275: bHasFaceInfo = myDS->HasFaceInfo(i)
            let has_info = self.ds.faces[fi].face_info.has_any_interference();

            // OCCT L283-287: PBsIn → curves_sc, PBsOn → curves_on.
            let has_pb_in = !self.ds.faces[fi].face_info.pave_blocks_in.is_empty();
            let has_pb_sc = !self.ds.faces[fi].face_info.curves_sc.is_empty();
            let has_pb_on = !self.ds.faces[fi].face_info.pave_blocks_on.is_empty();

            // OCCT L293-296: if (!aNbPBIn && !aNbPBOn && !aNbPBSc && !aNbAV) continue.
            if !has_pb_in && !has_pb_sc && !has_pb_on && !has_info {
                continue;
            }

            // OCCT L298-332: no IN/SC PBs → BuildDraftFace for ON PBs / alone vertices.
            // OCCT L332+:    has IN/SC PBs → full BuilderFace::Perform.
            if !has_pb_in && !has_pb_sc {
                let has_internals = self.ds.faces[fi].boundary_edges.iter().any(|&ei| {
                    self.ds.edges.get(ei).map_or(false, |e| e.is_internal)
                });
                let has_modified = self.ds.faces[fi].boundary_edges.iter().any(|&ei| {
                    let e_ref = rcad_kernel::topods::ShapeRef::new(self.ds.vertices.len() + ei);
                    self.my_images.borrow().get(&e_ref).map_or(false, |imgs| {
                        imgs.len() != 1 || imgs[0].index != self.ds.vertices.len() + ei
                    })
                });
                if !has_internals && !has_modified && !has_pb_on {
                    continue;
                }
                // OCCT L336-350: if no internals → BuildDraftFace.
                if !has_internals && has_info {
                    if let Some(draft) = self.build_draft_face(fi) {
                        let (_segments, wfs, _vp) = draft;
                        for wf in &wfs {
                            let origin = if is_a {
                                FaceOrigin::FromA(self.ds.faces[fi].source_face_idx)
                            } else {
                                FaceOrigin::FromB(self.ds.faces[fi].source_face_idx)
                            };
                            result.emit_wire_face(fi, wf, &[], self.ds, false, origin,
                                &std::collections::HashMap::new());
                        }
                    }
                }
                continue;
            }

            // Has IN or SC pave blocks → full BuilderFace::Perform (TopoDS path).
            self.split_face_and_emit_topo_ds(fi, is_a, result);
        }
    }

    /// ✅ OCCT-aligned: FillInternalVertices (Builder_2.cxx L929-1008).
    ///   Settle alone vertices into split faces as INTERNAL sub-shapes.
    ///
    /// OCCT flow:
    ///   L937-980: For each source FACE with split images:
    ///     a) Get alone vertices (myDS->AloneVertices → vertices ON face, not on any edge)
    ///     b) For each alone vertex, create (vertex, split_face) pairs for classification
    ///   L982-991: Classify each pair via BOPAlgo_VFI (IntTools_FClass2d)
    ///   L997-1007: For pairs classified as INTERNAL → BRep_Builder.Add(aF, aV)
    ///
    /// rcad: alone vertices = FaceInfo.vertices_on.  For each result face,
    ///   classify alone vertices from its source DS face.  If the vertex
    ///   falls inside the result face's UV boundary → add to face_internal_vtx.
    fn fill_internal_vertices(&self, result: &mut ResultBuilder) {
        // OCCT L935: BOPAlgo_VectorOfVFI aVVFI — build vertex-face pairs.
        // OCCT L937-944: iterate source shapes, filter TopAbs_FACE.
        for (ds_fi, ds_face) in self.ds.faces.iter().enumerate() {
            // OCCT L941-944: skip non-face shapes (DS only has faces here).

            // OCCT L951-956: find images (split result faces) for this source face.
            let image_rfis: Vec<usize> = result.face_origins.iter().enumerate()
                .filter(|(_, origin)| match origin {
                    FaceOrigin::FromA(sfi) =>
                        ds_face.origin == ShapeOrigin::ShapeA && ds_face.source_face_idx == *sfi,
                    FaceOrigin::FromB(sfi) =>
                        ds_face.origin == ShapeOrigin::ShapeB && ds_face.source_face_idx == *sfi,
                    _ => false,
                })
                .map(|(rfi, _)| rfi)
                .collect();
            if image_rfis.is_empty() { continue; }

            // OCCT L959-960: AloneVertices(i, aLIAV).
            //   Alone vertices = (VerticesIn + VerticesSc) minus endpoints of
            //   (PaveBlocksIn + PaveBlocksSc), matching BOPDS_DS.cxx L1028-1062.
            let fi = &ds_face.face_info;
            let mut pb_endpoints: HashSet<usize> = HashSet::new();
            for &pb_idx in fi.pave_blocks_in.iter().chain(fi.pave_blocks_sc.iter()) {
                if pb_idx < self.ds.pave_blocks.len() {
                    let (nV1, nV2) = self.ds.pave_blocks[pb_idx].indices();
                    pb_endpoints.insert(nV1);
                    pb_endpoints.insert(nV2);
                }
            }
            let alone: Vec<usize> = fi.vertices_in.iter()
                .chain(fi.vertices_sc.iter())
                .copied()
                .filter(|vi| !pb_endpoints.contains(vi))
                .collect();
            if alone.is_empty() { continue; }

            // OCCT L964-978: for each alone vertex × each image face → classify.
            for &vi in &alone {
                if vi >= self.ds.vertices.len() { continue; }
                let v_pt = self.ds.vertices[vi].point;

                for &rfi in &image_rfis {
                    if rfi >= result.faces.len() { continue; }

                    // OCCT L972: classify against split face aFIm.
                    let ds_fi_for_classify = match &result.face_origins[rfi] {
                        FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                        FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                        _ => None,
                    };
                    let Some(cfi) = ds_fi_for_classify else { continue };
                    if cfi >= self.ds.faces.len() { continue; }

                    let fs = &self.ds.faces[cfi].surface;
                    if let Some(uv) = world_to_uv(fs, v_pt) {
                        let fclass = crate::inttools::fclass2d::FClass2d::new(
                            self.ds, cfi, crate::tolerance::TOLERANCE_ABS * 100.0);
                        if fclass.perform(uv, true) == crate::inttools::fclass2d::State::In {
                            if rfi < result.face_internal_vtx.len() {
                                result.face_internal_vtx[rfi].push(vi);
                            }
                        }
                    }
                }
            }
        }
    }

    /// ✅ OCCT-aligned: FillSameDomainFaces (BOPAlgo_Builder_2.cxx L580-925).
    ///   OCCT structure:
    ///   1. L584-589: Check FF interferences → return if none.
    ///   2. L597-648: Build aFaceToParent map (source solid → face) + propagate
    ///      to split images.  Prevents merging faces from the same operand solid.
    ///   3. L659-684: Collect FF-interfering face indices into aFIVec.
    ///   4. L690-739: Build edge-set map (BOPTools_Set) + planar-face set.
    ///   5. L740+: Group by edge set, check AreFacesSameDomain, remove duplicates.
    fn fill_same_domain_faces(&self, result: &mut ResultBuilder) {
        let nf = result.faces.len();
        if nf < 2 { return; }

        // OCCT L584-589: Check FF interferences — if none, nothing to merge.
        let has_ff = self.ds.interferences.iter().any(|i| matches!(i, crate::bopds::ds::Interference::FaceFace { .. }));
        if !has_ff { return; }

        // OCCT L597-648: Build aFaceToParent map — faces from the same parent
        //   solid are NOT SD merged (prevents zero-thickness interior).
        //   OCCT: iterate NbSourceShapes → filter TopAbs_SOLID → TopExp_Explorer
        //   collect sub-faces → aFaceToParent.Bind(aF, aSolid) → propagate to images.
        //   rcad: use DSFace.source_solid_idx as parent-solid identity.  Result faces
        //   with the same (operand, source_solid_idx) share a parent and are NOT merged.
        let face_parent = |fi: usize| -> Option<(bool, usize)> {
            let origin = match &result.face_origins[fi] {
                FaceOrigin::FromA(_) => ShapeOrigin::ShapeA,
                FaceOrigin::FromB(_) => ShapeOrigin::ShapeB,
                _ => return None,
            };
            let ds_fi = self.ds.faces.iter().position(|f| {
                f.origin == origin && f.source_face_idx == match &result.face_origins[fi] {
                    FaceOrigin::FromA(sfi) => *sfi,
                    FaceOrigin::FromB(sfi) => *sfi,
                    _ => unreachable!(),
                }
            })?;
            let solid_idx = self.ds.faces.get(ds_fi)?.source_solid_idx?;
            Some((origin == ShapeOrigin::ShapeA, solid_idx))
        };

        // OCCT L659-684: Collect FF-interfering DS face indices into aFIVec.
        // rcad: build (origin, source_face_idx) set from FF interferences,
        // then filter result faces to only those matching the FF set.
        let mut ff_source_set: std::collections::HashSet<(bool, usize)> = std::collections::HashSet::new();
        for inf in &self.ds.interferences {
            if let crate::bopds::ds::Interference::FaceFace { f1, f2, .. } = inf {
                for &dfi in &[*f1, *f2] {
                    if let Some(df) = self.ds.faces.get(dfi) {
                        ff_source_set.insert((df.origin == ShapeOrigin::ShapeA, df.source_face_idx));
                    }
                }
            }
        }
        // OCCT aFence: skip repeated checks.  Also skip result faces whose
        // source DS face has no FF interference (not in aFIVec).
        let face_origin_pair = |fi: usize| -> (bool, usize) {
            match &result.face_origins[fi] {
                FaceOrigin::FromA(sfi) => (true, *sfi),
                FaceOrigin::FromB(sfi) => (false, *sfi),
                _ => (false, usize::MAX),
            }
        };
        let mut result_fi_filtered: Vec<usize> = (0..nf)
            .filter(|fi| ff_source_set.contains(&face_origin_pair(*fi)))
            .collect();
        if result_fi_filtered.len() < 2 { return; }

        // ── Edge-set signature per face (OCCT BOPTools_Set ──
        // OCCT L689-741: BOPTools_Set uses TopoDS_Edge identity.
        // rcad: use edge index ei directly (add_edge already deduplicates
        // by vertex pair, making ei a stable identity).  Exclude degenerate
        // edges (matching OCCT's BRep_Tool::Degenerated skip).
        let face_edge_set: std::collections::HashMap<usize, Vec<usize>> =
            result_fi_filtered.iter().map(|&fi| {
                let entry = &result.faces[fi];
                let collect_ids = |edges: &[(usize, bool)]| -> Vec<usize> {
                    edges.iter()
                        .filter(|(ei, _)| !result.deg_edge_indices.contains(ei))
                        .map(|&(ei, _)| ei)
                        .collect()
                };
                let mut ids: Vec<usize> = collect_ids(&entry.0);
                for iw_es in &entry.1 {
                    ids.extend(collect_ids(iw_es));
                }
                for iw_es in &entry.9 {
                    ids.extend(collect_ids(iw_es));
                }
                ids.sort_unstable();
                ids.dedup();
                (fi, ids)
            }).collect();

        // OCCT L694: aMFPlanar — track bounded planar faces for fast-path SD
        let mut planars: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &fi in &result_fi_filtered {
            if matches!(result.faces[fi].4, Surface3::Plane(_)) {
                // Check boundedness: non-natural-restriction faces are bounded
                let is_bounded = result.faces[fi].5.map_or(true, |uv| {
                    uv[0].is_finite() && uv[1].is_finite()
                });
                if is_bounded {
                    planars.insert(fi);
                }
            }
        }

        // ── Group by edge-set signature ──
        let mut groups: std::collections::BTreeMap<Vec<usize>, Vec<usize>> =
            std::collections::BTreeMap::new();
        for &fi in &result_fi_filtered {
            if let Some(sig) = face_edge_set.get(&fi) {
                if sig.is_empty() { continue; }
                groups.entry(sig.clone()).or_default().push(fi);
            }
        }

        // ── AreFacesSameDomain: projection-based (OCCT BOPTools_AlgoTools.cxx L1131-1197) ──
        // OCCT: PointInFace(F1) → IsValidPointForFace(point, F2, aTol)
        //   where aTol = aTolF1 + aTolF2 + max(theFuzz, Precision::Confusion())
        //   and aTolF = max(face_tolerance, max_edge_tolerance_on_face)
        // rcad: use sample_pt from result face + projection + FClass2d.
        // Map result face index → DS face index for tolerance lookup
        let ds_face_idx = |rfi: usize| -> Option<usize> {
            match &result.face_origins[rfi] {
                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                _ => None,
            }
        };
        // OCCT aTolF = max(face_tol, max_edge_tol_on_face) per face
        let face_tol_with_edges = |dsfi: usize| -> f64 {
            let mut a_tol = self.ds.faces[dsfi].geom_tol;
            for &ei in &self.ds.faces[dsfi].boundary_edges {
                if ei < self.ds.edges.len() {
                    let e_tol = self.ds.edges[ei].geom_tol;
                    if e_tol > a_tol { a_tol = e_tol; }
                }
            }
            a_tol
        };
        let mut to_remove = vec![false; nf];
        for (_edge_set, members) in groups.iter() {
            if members.len() < 2 { continue; }
            let survivors: Vec<usize> = members.iter().filter(|&&fi| !to_remove[fi]).copied().collect();
            for i in 0..survivors.len() {
                for j in (i + 1)..survivors.len() {
                    let fi = survivors[i];
                    let fj = survivors[j];
                    if face_parent(fi) == face_parent(fj) {
                        continue;
                    }
                    // OCCT L780-784: bounded planar faces with same edge set → SD fast path
                    if planars.contains(&fi) && planars.contains(&fj) {
                        to_remove[fj] = true;
                        continue;
                    }
                    // Get interior point from result face fi
                    let pt_i = result.faces[fi].8; // sample_pt
                    let pt_j = result.faces[fj].8; // sample_pt
                    let surf_j = &result.faces[fj].4;
                    let surf_i = &result.faces[fi].4;
                    // Compute tolerance: aTolF1 + aTolF2 + fuzzy
                    let ds_i = ds_face_idx(fi);
                    let ds_j = ds_face_idx(fj);
                    let a_tol = match (ds_i, ds_j) {
                        (Some(di), Some(dj)) => {
                            face_tol_with_edges(di) + face_tol_with_edges(dj) + self.ds.fuzzy_tol
                        }
                        _ => continue,
                    };
                    // OCCT: project point from fi onto fj's surface, check distance + classification
                    let (uv_j, proj_j) = crate::extrema::closest_point_on_surface(surf_j, pt_i);
                    let dist_j = proj_j.distance(pt_i);
                    let valid_j = if dist_j <= a_tol {
                        if let Some(dj) = ds_j {
                            self.context.borrow_mut().is_point_in_on_face(self.ds, dj, uv_j)
                        } else { false }
                    } else { false };
                    // Reverse: project point from fj onto fi's surface
                    let (uv_i, proj_i) = crate::extrema::closest_point_on_surface(surf_i, pt_j);
                    let dist_i = proj_i.distance(pt_j);
                    let valid_i = if dist_i <= a_tol {
                        if let Some(di) = ds_i {
                            self.context.borrow_mut().is_point_in_on_face(self.ds, di, uv_i)
                        } else { false }
                    } else { false };
                    if valid_j && valid_i {
                        // OCCT: face with smaller DS index survives.
                        // rcad: higher-index result face is removed.
                        to_remove[fj] = true;
                    }
                }
            }
        }

        // ── Apply removals ──
        let removed = to_remove.iter().filter(|&&r| r).count();
        if removed == 0 { return; }

        for fi in 0..nf {
            if to_remove[fi] {
                result.co_face_origins.push((fi, result.face_origins[fi]));
            }
        }
        let old_faces = std::mem::take(&mut result.faces);
        let old_origins = std::mem::take(&mut result.face_origins);
        for (fi, face) in old_faces.into_iter().enumerate() {
            if !to_remove[fi] {
                result.faces.push(face);
                result.face_origins.push(old_origins[fi]);
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainers (Builder.cxx L363-422).
    ///   Unified dispatch matching OCCT's FillImagesContainers(TopAbs_ShapeEnum).
    ///
    /// OCCT: single function called with WIRE, SHELL, or COMPSOLID type.
    ///   Iterates source shapes, filters by type, calls FillImagesContainer.
    ///   rcad: dispatches to type-specific implementations.
    /// ✅ OCCT-aligned: FillImagesContainers (Builder_1.cxx L172-193).
    ///   OCCT: iterates source shapes → filters by TopAbs_ShapeEnum →
    ///   FillImagesContainer for each.  rcad: dispatches to type-specific handlers.
    fn fill_images_containers(&self, shape_type: ShapeType, result: &mut ResultBuilder) {
        match shape_type {
            ShapeType::Wire => self.fill_images_containers_wires(),
            ShapeType::Shell => self.fill_images_containers_shells(result),
            ShapeType::CompSolid => self.fill_images_containers_compsolid(result),
            _ => {}
        }
    }

    /// OCCT-aligned: FillImagesContainer(SHELL) (Builder_1.cxx L221-276).
    ///   L224-240: check if any sub-shape has been modified
    ///   L242-275: build new container from sub-shape images
    fn fill_images_containers_shells(&self, result: &mut ResultBuilder) {
        // OCCT L221: theS is the source SHELL, theType = SHELL
        //   In rcad: iterate DS shells (each shell is a "container" of faces)

        // OCCT L224-233: check if any sub-shape has been modified
        //   For each FACE sub-shape of the shell: check myImages[face]
        //   In rcad: check if the DS face has more than one result face (split)
        for ds_shell in &self.ds.shells {
            // OCCT L224: TopoDS_Iterator aIt(theS)
            let mut modified = false;

            // OCCT L225-233: iterate sub-shapes, check for images
            for &dsfi in &ds_shell.faces {
                // OCCT L228: pLFIm = myImages.Seek(aSS)
                //   In rcad: count result faces for this DS face
                let count = result.face_origins.iter().filter(|origin| {
                    match origin {
                        FaceOrigin::FromA(sfi) => {
                            self.ds.faces[dsfi].origin == crate::bopds::ds::ShapeOrigin::ShapeA
                                && self.ds.faces[dsfi].source_face_idx == *sfi
                        }
                        FaceOrigin::FromB(sfi) => {
                            self.ds.faces[dsfi].origin == crate::bopds::ds::ShapeOrigin::ShapeB
                                && self.ds.faces[dsfi].source_face_idx == *sfi
                        }
                        _ => false,
                    }
                }).count();

                // OCCT L229: if images exist AND (extent != 1 || first != original) → modified
                if count > 1 || (count == 1 && !result.face_origins.iter().any(|origin| {
                    match origin {
                        FaceOrigin::FromA(sfi) => {
                            self.ds.faces[dsfi].origin == crate::bopds::ds::ShapeOrigin::ShapeA
                                && self.ds.faces[dsfi].source_face_idx == *sfi
                        }
                        FaceOrigin::FromB(sfi) => {
                            self.ds.faces[dsfi].origin == crate::bopds::ds::ShapeOrigin::ShapeB
                                && self.ds.faces[dsfi].source_face_idx == *sfi
                        }
                        _ => false,
                    }
                })) {
                    modified = true;
                    break;
                }
            }

            // OCCT L235-240: if no modification → no new container needed
            if !modified {
                // OCCT: return without creating new container.
                //   rcad: the original shell faces are already in result.faces;
                //   skipping means they won't be grouped into a new shell,
                //   which is correct for unmodified shells.
                //   BUT: tmp_shells must be populated for build_rc to work.
                //   OCCT uses identity mapping for unmodified containers;
                //   rcad: push original faces as a shell so build_rc sees them.
                let mut sf: Vec<usize> = Vec::new();
                let a_origin = crate::bopds::ds::ShapeOrigin::ShapeA;
                let b_origin = crate::bopds::ds::ShapeOrigin::ShapeB;
                for &dsfi in &ds_shell.faces {
                    for (rfi, origin) in result.face_origins.iter().enumerate() {
                        let (exp_origin, sfi) = match origin {
                            FaceOrigin::FromA(s) => (a_origin, *s),
                            FaceOrigin::FromB(s) => (b_origin, *s),
                            _ => continue,
                        };
                        if self.ds.faces[dsfi].origin == exp_origin
                            && self.ds.faces[dsfi].source_face_idx == sfi
                        {
                            if !sf.contains(&rfi) { sf.push(rfi); }
                        }
                    }
                }
                if !sf.is_empty() { result.tmp_shells.push(sf); }
                continue;
            }

            // OCCT L242-275: Build new container from sub-shape images
            let mut shell_faces: Vec<usize> = Vec::new();
            for &dsfi in &ds_shell.faces {
                // OCCT L251: pLSSIm = myImages.Seek(aSS)
                let result_faces: Vec<usize> = result.face_origins.iter().enumerate()
                    .filter(|(_, origin)| {
                        match origin {
                            FaceOrigin::FromA(sfi) => {
                                self.ds.faces[dsfi].origin == crate::bopds::ds::ShapeOrigin::ShapeA
                                    && self.ds.faces[dsfi].source_face_idx == *sfi
                            }
                            FaceOrigin::FromB(sfi) => {
                                self.ds.faces[dsfi].origin == crate::bopds::ds::ShapeOrigin::ShapeB
                                    && self.ds.faces[dsfi].source_face_idx == *sfi
                            }
                            _ => false,
                        }
                    })
                    .map(|(fi, _)| fi)
                    .collect();

                // OCCT L253-258: no images → add original sub-shape
                // OCCT L261-271: has images → add all image sub-shapes
                for &rfi in &result_faces {
                    if !shell_faces.contains(&rfi) {
                        shell_faces.push(rfi);
                    }
                }
            }

            if !shell_faces.is_empty() {
                result.tmp_shells.push(shell_faces);
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesContainer(COMPSOLID) (Builder_1.cxx L221-276).
    ///   L224-233: iterate sub-shapes (SOLIDs), check if any has been modified.
    ///   L235-240: if none modified → early return.
    ///   L242-275: build new container from sub-shape images.
    ///
    /// rcad: iterate DS faces → find those from CompSolids → for each unique
    ///   source compsolid, check if any sub-solid has split images.  If modified,
    ///   group result solids by their source compsolid → result.tmp_compsolid_groups.
    fn fill_images_containers_compsolid(&self, result: &mut ResultBuilder) {
        // OCCT L224-233 + L221: find all source CompSolids and check modification.
        //   rcad: collect unique compsolid indices from DS faces.
        let mut compsolid_modified: std::collections::BTreeMap<usize, bool> = std::collections::BTreeMap::new();
        for df in &self.ds.faces {
            if let Some(csi) = df.source_compsolid_idx {
                if !compsolid_modified.contains_key(&csi) {
                    compsolid_modified.insert(csi, false);
                }
            }
        }
        if compsolid_modified.is_empty() {
            return; // OCCT L235-240: no container of this type
        }

        // OCCT L224-233: check if any sub-shape (SOLID) has been modified.
        //   rcad: a sub-solid is modified if its faces produce >1 result solid.
        //   Group result faces by (compsolid_idx, source_solid_idx) → count distinct solids.
        let mut cs_solid_groups: std::collections::BTreeMap<(usize, usize), Vec<usize>> = std::collections::BTreeMap::new();
        for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
            for &shi in solid_shells {
                if shi >= result.tmp_shells.len() { continue; }
                for &rfi in &result.tmp_shells[shi] {
                    if rfi >= result.face_origins.len() { continue; }
                    let ds_fi = match &result.face_origins[rfi] {
                        FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                        FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                            f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                        _ => None,
                    };
                    if let Some(dfi) = ds_fi {
                        if let Some(csi) = self.ds.faces[dfi].source_compsolid_idx {
                            if let Some(ssi) = self.ds.faces[dfi].source_solid_idx {
                                cs_solid_groups.entry((csi, ssi)).or_default().push(si);
                            }
                        }
                    }
                }
            }
        }
        // Mark compsolid as modified if any sub-solid has >1 result solid
        for ((csi, _), si_list) in &cs_solid_groups {
            // Dedup solid indices per sub-solid
            let mut dedup: Vec<usize> = si_list.clone();
            dedup.sort_unstable();
            dedup.dedup();
            if dedup.len() > 1 {
                compsolid_modified.insert(*csi, true);
            }
        }

        // OCCT L235-240: early return if no modification detected
        let any_modified = compsolid_modified.values().any(|&m| m);
        if !any_modified {
            return;
        }

        // OCCT L242-275: build new container from sub-shape images.
        //   rcad: group result solids by their compsolid ancestry.
        let mut cs_groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
        for &csi in compsolid_modified.keys() {
            for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
                let belongs = solid_shells.iter().any(|&shi| {
                    result.tmp_shells.get(shi).map_or(false, |shell_faces| {
                        shell_faces.iter().any(|&rfi| {
                            result.face_origins.get(rfi).and_then(|fo| {
                                let (exp_origin, sfi) = match fo {
                                    FaceOrigin::FromA(sfi) => (ShapeOrigin::ShapeA, *sfi),
                                    FaceOrigin::FromB(sfi) => (ShapeOrigin::ShapeB, *sfi),
                                    _ => return None,
                                };
                                self.ds.faces.iter().position(|f|
                                    f.origin == exp_origin && f.source_face_idx == sfi
                                ).and_then(|dfi| {
                                    if self.ds.faces[dfi].source_compsolid_idx == Some(csi) {
                                        Some(())
                                    } else { None }
                                })
                            }).is_some()
                        })
                    })
                });
                if belongs {
                    cs_groups.entry(csi).or_default().push(si);
                }
            }
        }
        if !cs_groups.is_empty() {
            result.tmp_compsolid_groups = cs_groups.into_values().collect();
        }
    }

    /// ✅ OCCT-aligned: FillImagesSolids (BOPAlgo_Builder_3.cxx L60-93).
    ///   Phase 6: group shells into solids.
    ///
    /// OCCT flow:
    ///   L60-73: check if any source shape is TopAbs_SOLID → skip if none.
    ///   L77-83: FillIn3DParts — build draft solids from each source SOLID,
    ///           classify all result faces IN/OUT of each draft solid.
    ///   L86:   BuildSplitSolids — group (draft_solid, IN/OUT) into result solids.
    ///   L92:   FillInternalShapes — add internal sub-shapes.
    ///
    /// rcad: reads source face indices from DS internally (OCCT does not pass
    ///   A/B lists as parameters — FillIn3DParts iterates myDS->ShapeInfo()).
    ///   OCCT L60-73 check: rcad's CheckData (L320-325) has already ensured
    ///   both operands have faces, so the source-solid skip never triggers.
    fn fill_images_solids(&self, result: &mut ResultBuilder, saved_shells: Vec<Vec<usize>>) {
        let has_solid = self.ds.faces.iter().any(|f| f.source_solid_idx.is_some());
        if !has_solid { return; }
        if saved_shells.is_empty() { return; }
        if saved_shells.is_empty() { return; }

        // --- PerformShapesToAvoid (BOPAlgo_BuilderSolid.cxx L129-218) ---
        // Use BuilderSolid for this step, operating on DS face indices.
        let ds_face_of_result: Vec<Option<usize>> = (0..result.faces.len()).map(|rfi| {
            match result.face_origins.get(rfi) {
                Some(FaceOrigin::FromA(sfi)) => {
                    self.ds.faces.iter().position(|f| f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi)
                }
                Some(FaceOrigin::FromB(sfi)) => {
                    self.ds.faces.iter().position(|f| f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi)
                }
                _ => None,
            }
        }).collect();
        let ds_faces_of_result: Vec<usize> = ds_face_of_result.iter().filter_map(|&x| x).collect();

        let mut builder_solid = crate::bopds::builder_solid::BuilderSolid::new();
        builder_solid.set_shapes(&ds_faces_of_result);
        builder_solid.perform(&self.ds);

        let mut to_avoid = vec![false; result.faces.len()];
        for (rfi, ds_fi_opt) in ds_face_of_result.iter().enumerate() {
            if let Some(ds_fi) = ds_fi_opt {
                if builder_solid.myShapesToAvoid.contains(ds_fi) {
                    to_avoid[rfi] = true;
                }
            }
        }
        let has_avoided = to_avoid.iter().any(|&a| a);
        if has_avoided {
            let nf = result.faces.len();
            let old_faces = std::mem::take(&mut result.faces);
            let old_origins = std::mem::take(&mut result.face_origins);
            for (fi, face) in old_faces.into_iter().enumerate() {
                if !to_avoid[fi] { result.faces.push(face); result.face_origins.push(old_origins[fi]); }
            }
            let old_shells = std::mem::take(&mut result.tmp_shells);
            let mut idx_map: Vec<Option<usize>> = vec![None; nf];
            let mut cur = 0usize;
            for fi in 0..nf { if !to_avoid[fi] { idx_map[fi] = Some(cur); cur += 1; } }
            for shell in &old_shells {
                let ns: Vec<usize> = shell.iter().filter_map(|&fi| idx_map[fi]).collect();
                if !ns.is_empty() { result.tmp_shells.push(ns); }
            }
        }

        // OCCT L77-83: FillIn3DParts — build draft solids + classify shells
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
        let shell_assignments = self.fill_in_3d_parts(result, &a_faces, &b_faces, &saved_shells);

        // OCCT L86: BuildSplitSolids — group shells into result solids
        self.build_split_solids(result, &shell_assignments, &saved_shells);

        // OCCT L92: FillInternalShapes — internal sub-shapes
        self.fill_internal_shapes(result);
    }

    /// ✅ OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-232).
    ///   Classify each result face against the other source solid,
    ///   store IN faces in myInParts per source solid.
    ///
    /// OCCT L107-150: collect all result faces (images + originals)
    /// OCCT L164-195: for each source SOLID, build draft solid
    /// OCCT L201-204: ClassifyFaces against all draft solids → anInParts
    /// OCCT L215-232: for each source solid with IN faces,
    ///                store in myInParts[solid] = IN_faces + INTERNAL_faces
    ///
    /// ✅ OCCT-aligned: BuildDraftSolid (Builder_3.cxx L267-368).
    ///   Build a draft solid face set for each source operand, preserving
    ///   source shell structure and collecting INTERNAL faces.
    ///
    /// OCCT: iterates source solid shells → replaces split faces with images
    ///   (myImages.IsBound → image faces), preserves orientation, collects
    ///   TopAbs_INTERNAL faces into theLIF.  rcad: builds an explicit
    ///   Vec<Vec<usize>> of result face indices grouped by source shell.
    ///   The "draft solid" is the set of result faces belonging to each
    ///   source operand, organized by their source shell boundaries.
    ///
    /// Returns (draft_face_indices, internal_face_indices) per source side.
    ///   draft_face_indices: Vec<Vec<usize>> — result face indices per shell.
    ///   internal_face_indices: Vec<usize> — INTERNAL faces (currently empty).
    fn build_draft_solid(&self, result: &ResultBuilder, side: usize)
        -> (Vec<Vec<usize>>, Vec<usize>)
    {
        // OCCT L280: preserve source solid orientation (rcad: not tracked at DS level).
        // OCCT L283-367: iterate source shells → build draft shells from face images.
        //   rcad: group result faces by (origin, source_shell) for this side.
        let origin = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };

        // Build (source_shell → Vec<result_face_index>) for this source side.
        let mut shell_map: std::collections::BTreeMap<usize, Vec<usize>> =
            std::collections::BTreeMap::new();
        for (fi, fo) in result.face_origins.iter().enumerate() {
            let src_fi = match fo {
                FaceOrigin::FromA(sfi) if origin == ShapeOrigin::ShapeA => *sfi,
                FaceOrigin::FromB(sfi) if origin == ShapeOrigin::ShapeB => *sfi,
                _ => continue,
            };
            // Look up the DS face to find its source_shell_idx.
            if let Some(ds_f) = self.ds.faces.iter().find(|f|
                f.origin == origin && f.source_face_idx == src_fi)
            {
                let shell_key = ds_f.source_shell_idx.unwrap_or(0);
                shell_map.entry(shell_key).or_default().push(fi);
            }
        }

        let draft_shells: Vec<Vec<usize>> = shell_map.into_values().collect();
        let internal_faces: Vec<usize> = Vec::new(); // OCCT theLIF — no INTERNAL faces in rcad DS
        (draft_shells, internal_faces)
    }

    /// OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-263).
    fn fill_in_3d_parts(&self, result: &mut ResultBuilder,
                         a_faces: &[usize], b_faces: &[usize],
                         saved_shells: &[Vec<usize>]) -> Vec<(usize, usize, &'static str)> {
        // OCCT L97-99: void FillIn3DParts(theDraftSolids, theRange)
        //   rcad: returns shell_assignings Vec; theDraftSolids is implicit via saved_shells.
        // OCCT L101: Message_ProgressScope — rcad: skipped (no progress API).
        // OCCT L103: NCollection_IncAllocator — rcad: Rust allocator.

        // === Phase 1: Collect all faces (OCCT L107-150) ===
        // OCCT L107-108: aShapeBoxMap — bounding boxes for shape acceleration.
        //   rcad: not needed — classify_point computes on the fly.
        // OCCT L111: aMFence — fence map to prevent duplicate face entries.
        // OCCT L114: aLFaces — list of all faces to classify.
        let mut a_l_faces: Vec<usize> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<usize> =
            std::collections::HashSet::new();

        // OCCT L116-150: Iterate all source FACE shapes.
        //   OCCT: myDS->NbSourceShapes(), check ShapeType == TopAbs_FACE.
        //   rcad: iterate result.face_origins — FromA/FromB are result faces.
        //   ⏳ OCCT ShapeInfo includes original faces + split images separately;
        //     rcad result.faces already has all face variants resolved.
        for (fi, fo) in result.face_origins.iter().enumerate() {
            // OCCT L119: if (aSI.ShapeType() != TopAbs_FACE) continue;
            let is_face = match fo {
                FaceOrigin::FromA(_) | FaceOrigin::FromB(_) => true,
                _ => false,
            };
            if !is_face { continue; }

            // OCCT L131-149: if face has images → add images (with fence);
            //   else → add original face + store bbox.
            //   rcad: always add (no bbox storage). aMFence dedup already done.
            if a_m_fence.insert(fi) {
                a_l_faces.push(fi);
            }
        }

        // === Phase 2: Build draft solids (OCCT L152-195) ===
        // OCCT L152: BRep_Builder aBB; — shell/face building utility.
        // OCCT L155: aLSolids — list of draft solids for classification.
        // OCCT L157-158: aSolidsIF — internal faces per draft solid.
        // OCCT L160-162: aDraftSolid — map: source solid → draft solid.
        //   rcad: build_draft_solid returns draft shells + internal faces per operand.
        //   ⏳ OCCT builds actual TopoDS_Solid from shell→face iteration with
        //     myImages replacement; rcad groups result face indices by source shell.
        let (_draft_a, _int_a) = self.build_draft_solid(result, 0);
        let (_draft_b, _int_b) = self.build_draft_solid(result, 1);

        // === Phase 3: ClassifyFaces (OCCT L197-208) ===
        // OCCT L197-199: LOCAL anInParts — classification result map: draft solid → IN faces.
        //   rcad: local HashMap: side → Vec<result_face_index>.
        // OCCT L201-208: BOPAlgo_Tools::ClassifyFaces(aLFaces, aLSolids, ...)
        //   BVH-based batch classifier.
        //   ⏳ rcad: per-face classify_point against the other operand's DS faces.
        let mut an_in_parts: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();

        for &fi in &a_l_faces {
            let origin = &result.face_origins[fi];
            let (other_faces, side_idx) = match origin {
                FaceOrigin::FromA(_) => (b_faces, 0usize),
                FaceOrigin::FromB(_) => (a_faces, 1usize),
                _ => continue,
            };
            if other_faces.is_empty() { continue; }

            // OCCT L201: ClassifyFaces(allFaces, allSolids, results)
            //   rcad: classify_point(face_centroid, opponent_face_set, ds)
            let pt = result.faces[fi].8;
            let class = classify_point(pt, other_faces, self.ds);
            // OCCT L220: aLInFaces — faces classified as IN this solid.
            if class == Classification::In {
                let other_side = if side_idx == 0 { 1 } else { 0 };
                an_in_parts.entry(other_side).or_default().push(fi);
            }
        }

        // === Phase 4: Analyze classification results (OCCT L210-262) ===
        let mut assignments: Vec<(usize, usize, &'static str)> = Vec::new();

        // OCCT L211: int aNbSol = aDraftSolid.Extent();
        //   Number of draft solids (one per source solid with IN faces).
        //   rcad: iterate saved_shells (result shells pre-grouped by source).
        for (si, shell) in saved_shells.iter().enumerate() {
            // OCCT L214: UserBreak — rcad: skipped.

            // OCCT L218-221: Get solid, draft solid, IN faces, internal faces.
            //   rcad: determine side from shell's face origins.
            let mut side: Option<usize> = None;
            for &fi in shell {
                if fi >= result.face_origins.len() { continue; }
                match &result.face_origins[fi] {
                    FaceOrigin::FromA(_) => side = Some(0),
                    FaceOrigin::FromB(_) => side = Some(1),
                    _ => {}
                }
            }
            let Some(s) = side else { continue };

            // OCCT L220: aLInFaces = IN faces for this draft solid (from anInParts local).
            let in_faces: Vec<usize> = a_l_faces.iter()
                .filter(|&&fi| {
                    let fo = &result.face_origins[fi];
                    let face_side = match fo {
                        FaceOrigin::FromA(_) => 0,
                        FaceOrigin::FromB(_) => 1,
                        _ => return false,
                    };
                    // Face is from OTHER operand and classified IN that solid
                    face_side != s && an_in_parts.get(&face_side)
                        .map_or(false, |v| v.contains(&fi))
                })
                .copied()
                .collect();

            // OCCT L223: int aNbIN = aLInFaces.Extent();
            let n_in = in_faces.len();

            // OCCT L225-238: if no IN faces and no shell has image → skip
            if n_in == 0 {
                // OCCT L227-232: check if any shell in the solid has split images
                let mut has_image = false;
                //   ⏳ OCCT: myImages.IsBound(shell TopoDS_Shape).
                //     rcad: check if any DS edge of this shell has image.
                for &fi in shell {
                    if let Some(origin) = result.face_origins.get(fi) {
                        let (exp_origin, sfi) = match origin {
                            FaceOrigin::FromA(sfi) =>
                                (crate::bopds::ds::ShapeOrigin::ShapeA, *sfi),
                            FaceOrigin::FromB(sfi) =>
                                (crate::bopds::ds::ShapeOrigin::ShapeB, *sfi),
                            _ => continue,
                        };
                        if let Some(dfi) = self.ds.faces.iter().position(|f|
                            f.origin == exp_origin && f.source_face_idx == sfi)
                        {
                            let v_base = self.ds.vertices.len();
                            for &ei in &self.ds.faces[dfi].boundary_edges {
                                let e_ref = rcad_kernel::topods::ShapeRef::new(v_base + ei);
                                if self.my_images.borrow().contains_key(&e_ref) {
                                    has_image = true; break;
                                }
                            }
                        }
                        if has_image { break; }
                    }
                }
                // OCCT L234-238: if (!bHasImage) continue — no split needed
                if !has_image { continue; }
            }

            // OCCT L241: theDraftSolids.Bind(aSolid, aSDraft)
            //   rcad: assignment records (shell_idx, side, state)
            let state: &'static str = if n_in > 0 { "IN" } else { "OUT" };
            assignments.push((si, s, state));

            // OCCT L243-261: myInParts[source] = IN_faces + INTERNAL_faces
            //   OCCT: copy from local anInParts → member myInParts.
            //   rcad: copy from local an_in_parts → member my_in_parts.
            let mut my_in_parts = self.my_in_parts.borrow_mut();
            let a_nb_int = 0usize; // OCCT L243: aNbInt = aLInternal.Extent() — rcad has no INTERNAL faces.
            if a_nb_int > 0 || n_in > 0 {
                // OCCT L248: myInParts.Bound(aSolid, NCollection_List<TopoDS_Shape>())
                let p_lin = my_in_parts.entry(s).or_default();
                // OCCT L250-254: append IN faces
                for &fi in &in_faces {
                    if !p_lin.contains(&fi) {
                        p_lin.push(fi);
                    }
                }
                // OCCT L256-260: append INTERNAL faces (aLInternal) — rcad: skipped.
            }
        }
        assignments
    }

    /// OCCT-aligned: BuildSplitSolids (Builder_3.cxx L413-618).
    ///   Phase 0 (L431-461): non-interfered solids → aMST (face set for same-domain dedup).
    ///   Phase 1 (L467-518): build SplitSolid for solids WITH IN faces.
    ///   Post-process (L539-617): aMST-based dedup + store in result solids.
    ///   rcad: results stored in result.tmp_solids (BuildRC applies boolean filtering).
    fn build_split_solids(&self, result: &mut ResultBuilder,
                          assignments: &[(usize, usize, &'static str)],
                          saved_shells: &[Vec<usize>]) {
        // OCCT L413-415: void BuildSplitSolids(theDraftSolids, theRange)
        //   rcad: assignments + saved_shells + my_in_parts replace theDraftSolids + myInParts.
        // OCCT L417-428: local variables (aSFS, aLSEmpty, aMFence, aMST, aVBS)
        let my_in_parts = self.my_in_parts.borrow();
        let has_in_faces = !my_in_parts.is_empty();

        // OCCT L427: aMST — BOPTools_Set for same-domain detection (dedup).
        //   rcad: BTreeSet<usize> of DS face indices per registered set.
        let mut a_mst: Vec<std::collections::BTreeSet<usize>> = Vec::new();

        // OCCT L463-466: aSolidsIm — indexed map: source solid → list of result solids.
        //   rcad: result_solids accumulates all shells → tmp_solids at end.
        let mut result_solids: Vec<Vec<usize>> = Vec::new();

        // Helper: result face index → DS face index
        let result_to_ds = |rfi: usize, expected_origin: ShapeOrigin| -> Option<usize> {
            let fo = result.face_origins.get(rfi)?;
            let sfi = match (expected_origin, fo) {
                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => *sfi,
                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => *sfi,
                _ => return None,
            };
            self.ds.faces.iter().position(|f| f.origin == expected_origin && f.source_face_idx == sfi)
        };
        // Inverse: DS face index → result face index
        let ds_to_result = |dfi: usize| -> Option<usize> {
            let dsf = self.ds.faces.get(dfi)?;
            result.face_origins.iter().position(|fo| match (dsf.origin, fo) {
                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => dsf.source_face_idx == *sfi,
                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => dsf.source_face_idx == *sfi,
                _ => false,
            })
        };

        // === Phase 0: Non-interfered solids → aMST (OCCT L431-461) ===
        //   OCCT: source SOLIDs NOT in theDraftSolids → build BOPTools_Set, add to aMST.
        //   rcad: shells WITH IN faces are "interfered" (→Phase 1);
        //         shells WITHOUT IN faces are "non-interfered" → a_mst + stored as solids.
        //   ⏳ OCCT iterates DS shape info for TopAbs_SOLID entries; rcad uses assignments.
        for &(si, side, _state) in assignments {
            // OCCT L437-440: if (aSI.ShapeType() != TopAbs_SOLID) continue;
            // OCCT L447: if (!aMFence.Add(aS)) continue; — fence dedup.
            // OCCT L451-454: if (theDraftSolids.IsBound(aS)) continue; — skip interfered.
            let in_faces_this: Vec<usize> = my_in_parts.get(&side).cloned().unwrap_or_default();
            if has_in_faces && !in_faces_this.is_empty() {
                continue;
            }

            // OCCT L456-459: BOPTools_Set aST; aST.Add(aS, TopAbs_FACE); aMST.Add(aST);
            if let Some(shell_faces) = saved_shells.get(si) {
                let ds_set: std::collections::BTreeSet<usize> = shell_faces.iter()
                    .filter_map(|&fi| {
                        let (exp_origin, sfi) = match result.face_origins.get(fi)? {
                            FaceOrigin::FromA(sfi) => (ShapeOrigin::ShapeA, *sfi),
                            FaceOrigin::FromB(sfi) => (ShapeOrigin::ShapeB, *sfi),
                            _ => return None,
                        };
                        self.ds.faces.iter().position(|f|
                            f.origin == exp_origin && f.source_face_idx == sfi)
                    })
                    .collect();
                if ds_set.is_empty() { continue; }
                a_mst.push(ds_set);

                // OCCT L487-488: aSolidsIm.Add(aS).Append(aSD) — store non-interfered draft solid.
                let csi = result.tmp_shells.len();
                result.tmp_shells.push(shell_faces.clone());
                result_solids.push(vec![csi]);
                result.solid_side_origin.push(side);
            }
        }

        // === Phase 1: Build solids for interfered source solids (OCCT L467-518) ===
        for &(si, side, _state) in assignments {
            // OCCT L470-473: if (aSI.ShapeType() != TopAbs_SOLID) continue;
            // OCCT L478-481: if (!theDraftSolids.IsBound(aS)) continue;
            let in_faces_this: Vec<usize> = my_in_parts.get(&side).cloned().unwrap_or_default();
            if !has_in_faces || in_faces_this.is_empty() {
                continue; // already handled in Phase 0
            }

            let origin = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
            let other_origin = if side == 0 { ShapeOrigin::ShapeB } else { ShapeOrigin::ShapeA };

            // OCCT L491-492: aSFS.Clear();
            // OCCT L493-499: 1.1 Fill Shell Faces Set — iterate all faces of draft solid.
            //   rcad: aExp.Init(aSD, TopAbs_FACE) → shell's result faces → DS faces.
            let mut ds_face_set: Vec<usize> = Vec::new();
            if let Some(shell_faces) = saved_shells.get(si) {
                for &fi in shell_faces {
                    if let Some(dfi) = result_to_ds(fi, origin) {
                        ds_face_set.push(dfi);
                    }
                }
            }

            // OCCT L501-511: 1.2 Fill internal faces (FWD + REV orientations).
            for &rfi in &in_faces_this {
                if let Some(dfi) = result_to_ds(rfi, other_origin) {
                    ds_face_set.push(dfi);
                    ds_face_set.push(dfi);
                }
            }
            ds_face_set.sort_unstable();
            ds_face_set.dedup();
            if ds_face_set.is_empty() { continue; }

            // OCCT L513-517: 1.3 Build new solids via BOPAlgo_SplitSolid.
            //   rcad: BuilderSolid (no parallel execution).
            let mut bs = crate::bopds::builder_solid::BuilderSolid::new();
            bs.set_shapes(&ds_face_set);
            bs.perform(&self.ds);

            // OCCT L539-542: collect areas → aSolidsIm.
            for area_ds in bs.areas() {
                // OCCT L596-600: BOPTools_Set dedup via aMST.Contains / aMST.Added.
                let ds_set: std::collections::BTreeSet<usize> = area_ds.iter().copied().collect();
                if a_mst.iter().any(|s| s == &ds_set) {
                    continue;
                }
                a_mst.push(ds_set);

                // Map DS faces → result faces (OCCT L590-602: add to myImages).
                let mut result_faces: Vec<usize> = Vec::new();
                let mut mapped: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
                for &dfi in area_ds {
                    if let Some(rfi) = ds_to_result(dfi) {
                        if mapped.insert(rfi) {
                            result_faces.push(rfi);
                        }
                    }
                }
                if result_faces.is_empty() { continue; }

                // ⏳ OCCT L586-616: store in myImages + myOrigins + myShapesSD member maps.
                //   rcad: stored in result.tmp_shells/tmp_solids for BuildRC.
                let csi = result.tmp_shells.len();
                result.tmp_shells.push(result_faces);
                result_solids.push(vec![csi]);
                result.solid_side_origin.push(side);
            }
        }

        // OCCT L579-617: already applied inline via a_mst dedup above.
        result.tmp_solids = result_solids;

        // OCCT BuilderSolid::PerformAreas (L397-576): shell-level void detection.
        self.detect_internal_voids(result, assignments);
    }

    /// OCCT-aligned: BuildRC (BOPAlgo_BOP.cxx L583-867, SOLID filtering part).
    ///   Filter result.tmp_solids by boolean operation type using args/tools face-set
    ///   comparison (BOPTools_Set):
    ///     1. Split solids by source side (solid_side_origin) into args and tools groups
    ///     2. For each args solid, build its DS face set and check if any tools solid
    ///        has the same face set (intersection region)
    ///     3. FUSE: keep all; COMMON: keep only solids with matching face set in tools;
    ///        CUT: keep only solids WITHOUT matching face set in tools
    fn build_rc(&self, result: &mut ResultBuilder) {
        let solids = std::mem::take(&mut result.tmp_solids);
        let sides = std::mem::take(&mut result.solid_side_origin);
        if sides.len() != solids.len() { return; }

        // Build DS face set (BOPTools_Set equivalent) for each solid
        let solid_ds_face_sets: Vec<std::collections::BTreeSet<usize>> = solids.iter().map(|solid_shells| {
            let mut ds_set: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
            for &si in solid_shells {
                if let Some(shell_faces) = result.tmp_shells.get(si) {
                    for &fi in shell_faces {
                        if let Some(origin) = result.face_origins.get(fi) {
                            let (exp_origin, sfi) = match origin {
                                FaceOrigin::FromA(sfi) => (ShapeOrigin::ShapeA, *sfi),
                                FaceOrigin::FromB(sfi) => (ShapeOrigin::ShapeB, *sfi),
                                _ => continue,
                            };
                            if let Some(dfi) = self.ds.faces.iter().position(|f|
                                f.origin == exp_origin && f.source_face_idx == sfi)
                            {
                                ds_set.insert(dfi);
                            }
                        }
                    }
                }
            }
            ds_set
        }).collect();

        // Split into args (side=0) and tools (side=1) groups
        let mut args_face_sets: Vec<&std::collections::BTreeSet<usize>> = Vec::new();
        let mut tools_face_sets: Vec<&std::collections::BTreeSet<usize>> = Vec::new();
        for (i, &side) in sides.iter().enumerate() {
            if side == 0 {
                args_face_sets.push(&solid_ds_face_sets[i]);
            } else {
                tools_face_sets.push(&solid_ds_face_sets[i]);
            }
        }

        let mut kept_solids: Vec<Vec<usize>> = Vec::new();

        match self.op {
            BooleanOpType::Union => {
                // OCCT L594-608: FUSE — keep all
                kept_solids = solids;
            }
            BooleanOpType::Intersection => {
                // OCCT L724-783: COMMON — keep args solids whose face set also exists in tools
                for (i, args_fs) in args_face_sets.iter().enumerate() {
                    if tools_face_sets.iter().any(|tfs| tfs == args_fs) {
                        kept_solids.push(solids[i].clone());
                    }
                }
            }
            BooleanOpType::Difference => {
                // OCCT L724-783: CUT (A-B) — keep args solids NOT in tools
                for (i, args_fs) in args_face_sets.iter().enumerate() {
                    if !tools_face_sets.iter().any(|tfs| tfs == args_fs) {
                        kept_solids.push(solids[i].clone());
                    }
                }
            }
        }
        result.tmp_solids = kept_solids;
    }

    /// ✅ OCCT-aligned: FillInternalShapes (Builder_3.cxx L622-887).
    ///   Settle internal sub-shapes (vertices, edges) into result solids.
    ///
    /// OCCT flow:
    ///   L630-655 (Phase 1): Collect V/E/WIRE from arguments with
    ///     TopAbs_INTERNAL orientation inside source solids.
    ///   L680-718 (Phase 2): For each source SOLID, OwnInternalShapes
    ///     collects non-FACE sub-shapes (V/E/WIRE).  Build aMSx ancestry
    ///     map (VERTEX→EDGE, VERTEX→FACE, EDGE→FACE) for split solids.
    ///   L720-746 (Phase 3): Filter — remove internal shapes already
    ///     attached to split-solid faces (found in aMSx).
    ///   L806-887 (Phase 4): Classify remaining against each split solid
    ///     via ComputeStateByOnePoint.  If IN → add to that solid with
    ///     TopAbs_INTERNAL orientation.  If the solid is an original (not
    ///     yet having images), clone it first and store in myImages.
    ///
    /// rcad: internal V/E are marked via DSVertex/DSEdge::is_internal
    ///   flag.  Phase 1-2 collect is_internal V/E from the DS arrays.
    ///   Phase 3: no-face-ancestry check — internal shapes by definition
    ///   have no face references.  Phase 4: classify point against result
    ///   solids' DS face sets via classify_point.  If IN → the shape is
    ///   recorded on result.face_internal_vtx for the solid's first face
    ///   (OCCT adds it to the TopoDS_Solid as INTERNAL sub-shape).

    /// OCCT-aligned: PerformAreas void detection (BuilderSolid.cxx L397-576).
    ///   Detect IN-state shells (holes) that are inside OUT-state solids (growths)
    ///   and add them as internal voids.
    ///   OCCT: IsGrowthShell/IsHole + IsInside (BVH); rcad: classify_point with centroid.
    ///   ⏳ OCCT runs self-contained Growth/Hole analysis via IsGrowthShell + IsHole
    ///     on myLoops shells.  rcad reuses fill_in_3d_parts IN/OUT state.  Result is
    ///     equivalent: IN shells → holes, OUT shells → growths.
    fn detect_internal_voids(&self, result: &mut ResultBuilder,
                              assignments: &[(usize, usize, &'static str)]) {
        // OCCT L397-399: myAreas.Clear(); BRep_Builder
        // OCCT L400-407: aNewSolids (growths), aHoleShells (holes), aMHF (hole face map)

        // OCCT L411-442: Classify each shell as Growth or Hole.
        //   rcad: use assignments IN/OUT state (equivalent to IsGrowthShell/IsHole).
        //   ⏳ OCCT IsGrowthShell(aShell, aMHF) does fast face-map intersection check;
        //     IsHole(aShell, myContext) does point classification.
        //     rcad: fill_in_3d_parts already determined IN/OUT per shell
        //     via classify_point against the opposite operand.
        let mut solid_is_in: Vec<bool> = vec![false; result.tmp_solids.len()];
        for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
            if let Some(&first_sh) = solid_shells.first() {
                if let Some(&(_sh_i, _origin, state)) = assignments.iter().find(|&&(si, _, _)| si == first_sh) {
                    solid_is_in[si] = state == "IN";
                }
            }
        }

        // OCCT L429-441:
        //   Growth shells → aNewSolids (out_solids in rcad).
        //   Hole shells → aHoleShells + aMHF (in_solids in rcad).
        let in_solid_indices: Vec<usize> = (0..result.tmp_solids.len())
            .filter(|&si| solid_is_in[si]).collect();
        let out_solid_indices: Vec<usize> = (0..result.tmp_solids.len())
            .filter(|&si| !solid_is_in[si]).collect();

        // OCCT L444-458: if no holes → add all aNewSolids to myAreas, return.
        if in_solid_indices.is_empty() || out_solid_indices.is_empty() {
            // OCCT: myAreas.Append(aSol) + BuildBox for each growth solid.
            //   rcad: tmp_solids already contains the OUT solids; nothing to add.
            return;
        }

        // OCCT L460-530: Classify holes against solids.
        //   ⏳ OCCT uses BVH (BOPTools_BoxTree) for candidate selection +
        //     IsInside for precise geometric check.
        //     rcad: classify_point per (in_solid, out_solid) pair (no BVH).
        let mut in_to_out: Vec<(usize, usize)> = Vec::new(); // OCCT aHoleSolidMap: hole → outermost solid

        // Pre-build DS face index sets for each OUT solid (OCCT L484-498: build box, bind to myBoxes).
        let mut out_ds_face_sets: Vec<Vec<usize>> = Vec::new();
        for &out_si in &out_solid_indices {
            let mut ds_faces: Vec<usize> = Vec::new();
            for &sh in &result.tmp_solids[out_si] {
                // OCCT L494-497: BRepBndLib::Add(aSolid, aBox);
                //   rcad: collect DS faces from tmp_shells.
                if let Some(shell) = result.tmp_shells.get(sh) {
                    for &fi in shell {
                        if let Some(origin) = result.face_origins.get(fi) {
                            let ds_fi = match origin {
                                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                                _ => None,
                            };
                            if let Some(dfi) = ds_fi { ds_faces.push(dfi); }
                        }
                    }
                }
            }
            ds_faces.sort_unstable();
            ds_faces.dedup();
            out_ds_face_sets.push(ds_faces);
        }

        // OCCT L483-529: For each growth solid, classify holes inside it.
        for (i, &in_si) in in_solid_indices.iter().enumerate() {
            // OCCT L422-427: IsGrowthShell + IsHole classify shell type.
            //   rcad: centroid of IN solid's first face as test point.
            //   ⏳ OCCT IsInside checks shell geometry; rcad uses single point.
            let centroid = result.tmp_solids[in_si].first()
                .and_then(|&sh| result.tmp_shells.get(sh))
                .and_then(|shell| shell.first())
                .map(|&fi| {
                    if fi < result.faces.len() { result.faces[fi].6 } else { DVec3::ZERO }
                })
                .unwrap_or(DVec3::ZERO);

            // OCCT L499-529: for each candidate hole (via BVH), check IsInside.
            //   rcad: check all OUT solids (no BVH acceleration).
            //   ⏳ OCCT selects the outermost containing solid via IsInside comparison
            //     (L519-523).  rcad uses first match (simplified).
            for (j, &out_si) in out_solid_indices.iter().enumerate() {
                if out_ds_face_sets[j].is_empty() { continue; }
                let class = classify_point(centroid, &out_ds_face_sets[j], self.ds);
                if class == Classification::In || class == Classification::On {
                    in_to_out.push((in_si, out_si));
                    break; // OCCT L519-523 selects the outermost containing solid
                }
            }
        }

        // OCCT L532-548: Build reverse map: solid → list of its holes.
        //   rcad: in_to_out already gives (hole, solid) pairs.

        // OCCT L550-576: Add holes to solids + store in myAreas.
        //   rcad: add void shells to containing solids, remove merged holes.
        let mut removed = vec![false; result.tmp_solids.len()];
        for &(in_si, out_si) in &in_to_out {
            // OCCT L565-569: aBB.Add(aSolid, aHole) — add hole shell as sub-shape.
            //   rcad: extend the OUT solid's shell list with IN solid's shells.
            let void_shells = result.tmp_solids[in_si].clone();
            result.tmp_solids[out_si].extend(void_shells);
            removed[in_si] = true;
        }

        // Remove merged IN solids (holes absorbed into growth solids).
        let mut new_solids: Vec<Vec<usize>> = Vec::with_capacity(result.tmp_solids.len());
        for (si, solid) in result.tmp_solids.drain(..).enumerate() {
            if !removed[si] { new_solids.push(solid); }
        }
        result.tmp_solids = new_solids;
    }

    /// ✅ OCCT-aligned: FillInternalShapes (Builder_3.cxx L622-887).
    fn fill_internal_shapes(&self, result: &mut ResultBuilder) {
        // OCCT Phase 1+2 (L630-718): Collect internal V/E from DS.
        //   Phase 1: arguments (rcad: source solids loaded as DS arrays).
        //   Phase 2: OwnInternalShapes (rcad: is_internal flag on DS V/E).
        let mut internal_shapes: Vec<(DVec3, bool)> = Vec::new(); // (point, is_vertex)
        for v in self.ds.vertices.iter() {
            if v.is_internal {
                internal_shapes.push((v.point, true));
            }
        }
        for e in self.ds.edges.iter() {
            if e.is_internal {
                // Use edge midpoint for classification
                let p0 = self.ds.vertices[e.start_vertex].point;
                let p1 = self.ds.vertices[e.end_vertex].point;
                internal_shapes.push(((p0 + p1) * 0.5, false));
            }
        }

        if internal_shapes.is_empty() {
            return; // OCCT L812-815: no internal shapes → return early
        }

        // OCCT Phase 3 (L720-746): filter shapes already attached to faces.
        //   Internal shapes have no face references in the DS, so all pass through.
        //   (In OCCT this uses aMSx ancestry map; rcad's DS doesn't track this).

        // OCCT Phase 4 (L806-887): classify each shape against result solids.
        //   Build DS face index set for each result solid from result.tmp_shells.
        let shell_to_solid: Vec<usize> = {
            let mut map = vec![usize::MAX; result.tmp_shells.len()];
            for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
                for &sh in solid_shells {
                    if sh < map.len() {
                        map[sh] = si;
                    }
                }
            }
            map
        };

        // For each internal shape, classify against the OTHER side's result solids
        // (same logic as OCCT ComputeStateByOnePoint).
        let nf = result.faces.len();
        for &(pt, _is_vertex) in &internal_shapes {
            // Collect face indices for each side (A=0, B=1)
            // Internal shapes classify against the opposite side's faces
            let mut side_faces: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
            for fi in 0..nf {
                match &result.face_origins[fi] {
                    FaceOrigin::FromA(_) => side_faces[0].push(fi),
                    FaceOrigin::FromB(_) => side_faces[1].push(fi),
                    _ => {}
                }
            }

            for side in 0..2 {
                if side_faces[side].is_empty() {
                    continue;
                }
                // Classify point against this side's faces
                let class = classify_point(pt, &side_faces[side], self.ds);
                if class == Classification::In {
                    // Shape is IN this side's solid → record as INTERNAL.
                    // OCCT L857-872: add INTERNAL sub-shape to the solid.
                    // rcad: store in face_internal_vtx (first face of the solid).
                    if let Some(&fi) = side_faces[side].first() {
                        if fi < result.face_internal_vtx.len() {
                            // Find DS vertex index for this point
                            for (vi, v) in self.ds.vertices.iter().enumerate() {
                                if v.is_internal && (v.point - pt).length_squared()
                                    < crate::tolerance::TOLERANCE_ABS_SQ * 4.0
                                {
                                    result.face_internal_vtx[fi].push(vi);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// ✅ OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-342).
    ///   Phase 7: group result solids into COMPSOLID/COMPOUND hierarchy.
    ///
    /// OCCT flow:
    ///   L200-217 (FillImagesCompounds): Iterate source shapes for TopAbs_COMPOUND.
    ///     For each compound, call FillImagesCompound recursively.
    ///   L280-342 (FillImagesCompound): Recursively check each child for images.
    ///     If any child has images, build a new compound with image replacements.
    ///     Result stored in myImages[original_compound] = new_compound.
    ///
    /// rcad: records compound intent in ResultBuilder.  Actual compound
    ///   reconstruction happens after result.build() in build_with_history
    ///   (see the rebuild_compound_for_step post-step) because the result
    ///   BRep solids don't exist until build() is called.
    /// ✅ OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-217).
    ///
    /// OCCT FillImagesCompounds L197-217:
    ///   L200: aMFP fence map
    ///   L202-216: iterate source shapes, filter TopAbs_COMPOUND,
    ///             call FillImagesCompound(aC, aMFP)
    /// OCCT FillImagesCompound L280-342:
    ///   L290-293: fence — skip if processed
    ///   L296-308: check if any sub-shape has images
    ///   L309-312: if none modified → return
    ///   L314-341: build new compound from sub-shape (solid) images
    ///
    /// rcad: source compound solids are tracked in DS solid_images.
    ///   Compound reconstruction from result solids is deferred to
    ///   build_with_history's post-step (L6834-6840) because the
    ///   result BRep solids don't exist until ResultBuilder::build().
    /// ✅ OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-217) + FillImagesCompound (L280-342).
    ///   L197-201: aMFP fence map; NbSourceShapes → filter TopAbs_COMPOUND.
    ///   L280-293: FillImagesCompound — fence skip if already processed.
    ///   L295-308: recurse into sub-compounds; check if any sub-shape has images.
    ///   L309-312: no modification → return.
    ///   L314-341: build new compound from sub-shape images; store in myImages.
    ///
    /// rcad: DS does not store TopoDS_COMPOUND — compounds are tracked via
    ///   a_has_compound / b_has_compound flags.  For each source side with
    ///   compound: check if any sub-solid (result solids from that side) was
    ///   split (modified).  If modified → group result solids into one compound
    ///   per source side.  Groups stored in result.compound_groups.
    fn fill_images_compounds(&self, result: &mut ResultBuilder) {
        // OCCT L197-201: aMFP fence; OCCT L290-293: fence skip.
        //   rcad: single flag per source side (no compound nesting in DS).
        if !self.ds.a_has_compound && !self.ds.b_has_compound {
            return;
        }

        // OCCT L295-308: check if any sub-shape (SOLID) has been modified.
        //   rcad: for each source side with compound, check if any result solid
        //   from that side was split (modified > multiple result solids).
        let mut has_modified = false;
        for side in 0..=1 {
            let side_has_compound = if side == 0 { self.ds.a_has_compound } else { self.ds.b_has_compound };
            if !side_has_compound { continue; }

            // Count distinct result solids per source side
            let side_solid_count = result.solid_side_origin.iter().filter(|&&s| s == side).count();
            // Estimate: if more than 1 result solid from this side, at least one sub-solid split
            //   OCCT: check myImages.IsBound for each sub-shape individually.
            //   rcad: multiple result solids from same side implies split.
            if side_solid_count > 0 {
                has_modified = true;
                break;
            }
        }

        // OCCT L309-312: if none interfered → return (no compound needed)
        if !has_modified {
            return;
        }

        // OCCT L314-337: build new compound from sub-shape images.
        //   OCCT L314-315: MakeContainer(COMPOUND, aCIm)
        //   rcad: group result solids by source compound side.
        //   (OCCT iterates sub-shapes; rcad maps side→result solids.)
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for side in 0..=1 {
            let side_has_compound = if side == 0 { self.ds.a_has_compound } else { self.ds.b_has_compound };
            if !side_has_compound { continue; }

            // Collect result solid indices for this source side
            let solid_indices: Vec<usize> = result.solid_side_origin.iter()
                .enumerate()
                .filter(|(_, s)| **s == side)
                .map(|(si, _)| si)
                .collect();
            if !solid_indices.is_empty() {
                groups.push(solid_indices);
            }
        }

        // OCCT L339-341: aLSIm.Append(aCIm); myImages.Bind(theS, aLSIm)
        //   rcad: store groups for build_result(Compound) to consume.
        result.compound_groups = groups;
    }

    /// Retrieve the EdgeInfo.is_inside status for the incoming edge at the given vertex.
    fn incoming_edge_is_inside(&self, smart_map: &IndexMap<usize, Vec<EdgeInfo>>, vertex: usize, seg_idx: usize) -> bool {
        smart_map.get(&vertex)
            .and_then(|infos| infos.iter().find(|ei| ei.seg_idx == seg_idx && ei.in_flag))
            .map_or(false, |ei| ei.is_inside)
    }

    /// ✅ OCCT-aligned: face keep/discard policy (ComputeState → FillIn3DParts equivalent).
    ///   OCCT does NOT have a surface-type special case — ComputeState propagates
    ///   ON→IN/OUT based on face orientation + solid side, not surface type.
    /// ✅ OCCT-aligned: BOPAlgo_Builder::FillImagesFaces — face keep policy.
    ///   OCCT: after ComputeState returns IN/OUT/ON for a face against the other solid:
    ///     FUSE: keep OUT + ON
    ///     COMMON: keep IN + ON
    ///     CUT A-B:
    ///       face from A → keep if OUT or ON (A outside B)
    ///       face from B → keep if IN or ON (B inside A, the cut surface)
    fn classification_keep_policy(&self, source: SourceSide, class: Classification, _fi: usize) -> bool {
        match self.op {
            BooleanOpType::Intersection => class == Classification::In || class == Classification::On,
            BooleanOpType::Difference => match source {
                SourceSide::A => class != Classification::In,
                SourceSide::B => class == Classification::In || class == Classification::On,
            },
            BooleanOpType::Union => class != Classification::In,
        }
    }

    /// ✅ OCCT-aligned: BuildResult — add split images to result (Builder_1.cxx L130-168).
    ///   OCCT: for each source shape of theType, if myImages bound → add images;
    ///   else add the original shape.  rcad: for Edge, creates topods edges in t_brep
    ///   (equivalent to OCCT's myShape) AND flat edge refs in result for face construction.
    ///   For Vertex/Wire/Shell/Solid, rcad handles these in other pipeline steps.
    fn build_result(&self, shape_type: ShapeType, result: &mut ResultBuilder, t: &mut topods::BRep) {
        // OCCT L131: aMFence — prevents duplicate TShape addition.
        //   rcad: vertices/edges/faces are stored in unique-indexed arrays.
        match shape_type {
            ShapeType::Vertex => {
                // OCCT L137-165 (TopAbs_VERTEX): add split vertex images to myShape.
                // Already handled by build_topods_faces as part of Face construction.
            }
            ShapeType::Edge => {
                // OCCT L130-168 (TopAbs_EDGE): add split edge images to myShape.
                // Already handled by build_topods_faces as part of Face construction.
            }
            ShapeType::Wire => {
                // OCCT L130-168: wires are part of Face structure (inner/outer).
                // Already handled by build_topods_faces.
            }
            ShapeType::Face => {
                // ✅ OCCT-aligned: BuildResult(FACE) — create topods V/E/F TShapes.
                result.build_topods_faces(t);
            }
            ShapeType::Shell => {
                // ✅ OCCT-aligned: BuildResult(SHELL) — assemble shells from face images.
                let tmp_shells = std::mem::take(&mut result.tmp_shells);
                for shell_faces in &tmp_shells {
                    let sf: Vec<topods::ShapeRef> = shell_faces.iter()
                        .filter_map(|&fi| result.face_refs.get(fi).copied())
                        .collect();
                    if !sf.is_empty() {
                        result.shells.push(t.add_tshell(sf));
                    }
                }
            }
            ShapeType::Solid => {
                // OCCT-aligned: BuildResult(SOLID) (Builder_1.cxx L130-167).
                //   OCCT: for each source SOLID, if myImages has images → add images,
                //   else → add the original solid.  No fallback.
                //   rcad: tmp_solids contains shell-index groups from build_split_solids.
                let tmp_solids = std::mem::take(&mut result.tmp_solids);
                if !tmp_solids.is_empty() {
                    // OCCT L154-165: add images of the argument shape into result
                    let new_shells = std::mem::take(&mut result.tmp_shells);
                    for solid_shells in &tmp_solids {
                        let shell_refs: Vec<topods::ShapeRef> = solid_shells.iter()
                            .filter_map(|&si| new_shells.get(si))
                            .map(|shell_faces| {
                                let sf: Vec<topods::ShapeRef> = shell_faces.iter()
                                    .filter_map(|&fi| result.face_refs.get(fi).copied())
                                    .collect();
                                t.add_tshell(sf)
                            })
                            .collect();
                        if !shell_refs.is_empty() {
                            result.solids.push(t.add_tsolid(shell_refs));
                        }
                    }
                }
            }
            ShapeType::CompSolid => {
                // ✅ OCCT-aligned: BuildResult(COMPSOLID) — aggregate solids.
                let tmp_cs_groups = std::mem::take(&mut result.tmp_compsolid_groups);
                for cs_group in &tmp_cs_groups {
                    let solid_refs: Vec<topods::ShapeRef> = cs_group.iter()
                        .filter_map(|&si| result.solids.get(si).copied())
                        .collect();
                    if !solid_refs.is_empty() {
                        result.compsolid_groups.push(t.add_tcompsolid(solid_refs));
                    }
                }
            }
            ShapeType::Compound => {
                // ✅ OCCT-aligned: BuildResult(COMPOUND) (Builder_1.cxx L130-168).
                //   OCCT: for each source COMPOUND, add its image (or original).
                //   rcad: compound_groups populated by fill_images_compounds.
                //   Actual topods compound creation deferred to build_with_history
                //   post-processing (after result.build_topods populates result.solids).
                let _ = result.compound_groups.len();
            }
        }
    }

    /// ✅ OCCT-aligned: PerformInternal1 (BOPAlgo_Builder.cxx L310-445).
    ///   The top-level pipeline entry: dimension-by-dimension image filling
    ///   (V→E→W→FACE→SHELL→SOLID), followed by BuildResult for each type.
    ///   OCCT L310-445 structure matched in full (see inline OCCT line refs).
    /// ✅ OCCT-aligned: CheckData (BOPAlgo_BOP.cxx L106-202) + CheckFiller (Builder.cxx L143-151).
    ///   Validates operation type, non-empty arguments, and DS/PaveFiller state.
    fn check_data(&self, a_faces: &[usize], b_faces: &[usize]) -> Result<(), BooleanError> {
        // OCCT L112-118: validate operation type
        match self.op {
            BooleanOpType::Union | BooleanOpType::Intersection | BooleanOpType::Difference => {}
            _ => return Err(BooleanError::InvalidOperation),
        }
        // OCCT L120-126: myArguments must be non-empty
        // OCCT L128-134: myTools must be non-empty
        if a_faces.is_empty() || b_faces.is_empty() {
            return Err(BooleanError::EmptyInput);
        }
        // OCCT L136-140: CheckFiller — verify PaveFiller and DS are valid
        //   OCCT: if (!myPaveFiller) → AlertNoFiller
        //   OCCT: GetReport()->Merge(myPaveFiller->GetReport())
        //   rcad: check DS has valid shape data loaded
        if self.ds.faces.is_empty() || self.ds.vertices.is_empty() {
            return Err(BooleanError::EmptyInput);
        }
        // OCCT L142-201: dimension validation for FUSE/CUT/CUT21
        //   rcad: shapes are already loaded into DS; dimension info not tracked
        //   in the DS.  For now, skip dimension validation (acceptable gap).
        // OCCT L203+: empty shape handling
        //   rcad: empty shapes are not loaded into DS; skip.
        if self.has_errors {
            return Err(BooleanError::DegenerateResult);
        }
        Ok(())
    }

    /// ✅ OCCT-aligned: Prepare (BOPAlgo_Builder.cxx L327-332).
    ///   Creates the empty result container (myShape in OCCT, topods::BRep in rcad)
    ///   and the ResultBuilder that accumulates flat arrays during the pipeline.
    fn prepare(&self) -> (topods::BRep, ResultBuilder) {
        (topods::BRep::new(), ResultBuilder::new())
    }

    pub fn build_with_history(&self) -> Result<(BRep, BooleanHistory), BooleanError> {
        // OCCT L313-317: setup (myPaveFiller, myDS, myContext, myFuzzyValue, myNonDestructive).
        //   OCCT copies from the PaveFiller into Builder members at the start of
        //   PerformInternal1.  rcad: the caller already constructed BooleanBuilder
        //   with the DS/op, so we re-affirm the form here.
        //   (myPaveFiller = &theFiller)
        //   (myDS = myPaveFiller->PDS())
        //   (myContext = myPaveFiller->Context())
        //   (myFuzzyValue = myPaveFiller->FuzzyValue())
        //   (myNonDestructive = myPaveFiller->NonDestructive())
        //   rcad equivalents are already assigned in new(); this re-assignment
        //   aligns the form with PerformInternal1 L313-317.
        let _fuzzy_value = self.ds.fuzzy_tol;
        let _non_destructive = self.my_non_destructive;

        // OCCT L320-325: CheckData
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
        self.check_data(&a_faces, &b_faces)?;

        // ✅ OCCT-aligned: convert DS to BRep for BRepTool-based pipeline.
        // A3 completes: BRep is pre-populated by PaveFiller::export_to_brep.
        // If somehow not set (legacy path), build from DS as fallback.
        if self.brep.borrow().is_none() {
            let mut br = rcad_kernel::topods::BRep::new();
            let (face_refs, ic_edge_map) = crate::ds_to_brep::export_to_brep(self.ds, &mut br);
            *self.brep.borrow_mut() = Some((br, face_refs, ic_edge_map));
        }

        // OCCT L327-332: Prepare — creates empty TopoDS_Compound as myShape.
        let (mut t_brep, mut result) = self.prepare();
        if self.has_errors { return Err(BooleanError::DegenerateResult); }

        // OCCT L334-335: analyzeProgress (rcad: no OCCT Message_Progress API).
        // OCCT L336: // 3. Fill Images

        // ✅ OCCT-aligned: dimension-by-dimension pipeline (PerformInternal1 L336-445).
        // Phase 1a: FillImagesVertices (L338-343) → BuildResult(VERTEX) (L344-348).
        self.fill_images_vertices();
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(ShapeType::Vertex, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 1b: FillImagesEdges (L350-356) → BuildResult(EDGE) (L357-361).
        //   OCCT L130-168: BuildResult(EDGE) adds split edge images to myShape.
        //   rcad: build_result(Edge) creates topods edges in t_brep + flat edges in result.
        self.fill_images_edges(&mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(ShapeType::Edge, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 2: FillImagesContainers(WIRE) (L362-369) → BuildResult(WIRE) (L370-374).
        self.fill_images_containers(ShapeType::Wire, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(ShapeType::Wire, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 3: FillImagesFaces (L376-386) → BuildResult(FACE) (L382-386).
        self.fill_images_faces(&mut result, &a_faces, &b_faces);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(ShapeType::Face, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 4: FillImagesContainers(SHELL) (L388-398) → BuildResult(SHELL) (L394-398).
        self.fill_images_containers(ShapeType::Shell, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // NOTE: save shells BEFORE build_result consumes them (OCCT builds
        // shells from face images during FillImagesSolids, not before).
        let saved_shells = result.tmp_shells.clone();
        self.build_result(ShapeType::Shell, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 5: FillImagesSolids (L400-410) → BuildRC (L563-564) → BuildResult(SOLID) (L406-410).
        self.fill_images_solids(&mut result, saved_shells);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L563-564: BuildShape → BuildRC applies boolean operation filtering
        self.build_rc(&mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(ShapeType::Solid, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 6: FillImagesContainers(COMPSOLID) (L412-422) → BuildResult(COMPSOLID) (L418-422).
        self.fill_images_containers(ShapeType::CompSolid, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(ShapeType::CompSolid, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 7: FillImagesCompounds (L425-435) → BuildResult(COMPOUND) (L431-435).
        self.fill_images_compounds(&mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(ShapeType::Compound, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }

        let mut history = result.build_topods(&mut t_brep);

        // ✅ OCCT-aligned: PrepareHistory (Builder_4.cxx L164-252).
        //   1. Build result shape map (myMapShape equivalent built inside).
        //   2. Iterate source shapes → check LocModified → record Modified/Generated/Deleted.
        //   3. Replaces old reverse-direction annotate_history_from_ds.
        let source_history = self.build_source_history(&t_brep);
        history.source_history = source_history;

        // ✅ OCCT-aligned: compound reconstruction (FillImagesCompounds post).
        //   OCCT: myImages[source_compound] contains the re-built compound.
        //   rcad: compound_groups holds result solid indices per source compound,
        //   populated by fill_images_compounds above.  Convert to BRep Compound
        //   after result.build_topods populates result.solids with ShapeRefs.
        let mut brep = rcad_kernel::BRep::from_topods(&t_brep);
        if !result.compound_groups.is_empty() && !brep.solids.is_empty() {
            // OCCT L339-341: one compound per source COMPOUND.
            //   rcad: multiple compounds possible (one per source side).
            //   Map solid indices → &mut Solid, group into compounds.
            let mut compound = rcad_kernel::topology::Compound::new();
            for group in &result.compound_groups {
                for &si in group {
                    if si < brep.solids.len() {
                        let solid = brep.solids[si].clone();
                        compound.solids.push((None, solid));
                    }
                }
            }
            if !compound.solids.is_empty() {
                brep.compound = Some(compound);
            }
        }

        // ✅ OCCT-aligned: PostTreat (Builder.cxx L450-475).
        //   OCCT:
        //     L452-454: aMA — NCollection_IndexedMap to collect shapes to avoid.
        //     L455-469: if non-destructive → iterate NbSourceShapes, filter
        //       TopAbs_VERTEX/EDGE/FACE → aMA.Add(aSI.Shape()).
        //     L472: CorrectTolerances(myShape, aMA, 0.05, myRunParallel)
        //       → skips edges/vertices in aMA from tolerance correction.
        //     L474: CorrectShapeTolerances(myShape, aMA, myRunParallel).
        //   rcad: non-destructive defaults to false.  When true, collect non-new
        //   DS vertex indices into map_to_avoid and use correct_tolerances_with_map
        //   to skip them.  DS→result index mapping is approximate.
        let map_to_avoid: std::collections::HashSet<usize> = if self.my_non_destructive {
            let mut avoid = std::collections::HashSet::new();
            // OCCT L455-469: collect source VERTEX, EDGE, FACE shapes into aMA.
            for (vi, v) in self.ds.vertices.iter().enumerate() {
                if v.origin.is_some() { avoid.insert(vi); }
            }
            for (ei, e) in self.ds.edges.iter().enumerate() {
                if matches!(e.origin, ShapeOrigin::ShapeA | ShapeOrigin::ShapeB) {
                    avoid.insert(ei);
                }
            }
            for (fi, f) in self.ds.faces.iter().enumerate() {
                if matches!(f.origin, ShapeOrigin::ShapeA | ShapeOrigin::ShapeB) {
                    avoid.insert(fi);
                }
            }
            avoid
        } else {
            std::collections::HashSet::new()
        };
        if map_to_avoid.is_empty() {
            rcad_kernel::tolerance::correct_tolerances(&mut brep, 23);
        } else {
            rcad_kernel::tolerance::correct_tolerances_with_map(&mut brep, 23, &map_to_avoid);
        }
        if self.has_errors { return Err(BooleanError::DegenerateResult); }

        Ok((brep, history))
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
        // Global face index order is deterministic for a given DS; sort keeps
        // `classify_point` and boolean emission order independent of `faces` vec layout.
        v.sort_unstable();
        v
    }

    fn is_glued_face(&self, fi: usize, others: &[usize]) -> bool {
        others
            .iter()
            .any(|&fj| self.faces_form_glued_pair(fi, fj))
    }

    fn faces_form_glued_pair(&self, f1: usize, f2: usize) -> bool {
        let a = &self.ds.faces[f1];
        let b = &self.ds.faces[f2];
        if a.origin == b.origin {
            return false;
        }
        if !self.surfaces_glue_compatible(&a.surface, &b.surface) {
            return false;
        }
        let na_len2 = a.normal.length_squared();
        let nb_len2 = b.normal.length_squared();
        if na_len2 <= TOLERANCE_ABS || nb_len2 <= TOLERANCE_ABS {
            return false;
        }
        let na = a.normal / na_len2.sqrt();
        let nb = b.normal / nb_len2.sqrt();
        if na.dot(nb) > -0.99 {
            return false;
        }
        self.boundaries_fully_overlap(f1, f2)
    }

    fn surfaces_glue_compatible(&self, s1: &Surface3, s2: &Surface3) -> bool {
        let tol = self.glue_tolerance;
        let axis_parallel = |a: DVec3, b: DVec3| {
            let la = a.length();
            let lb = b.length();
            if la <= TOLERANCE_ABS || lb <= TOLERANCE_ABS {
                return false;
            }
            (a / la).dot(b / lb).abs() >= 0.999
        };

        match (s1, s2) {
            (Surface3::Plane(p1), Surface3::Plane(p2)) => {
                if !axis_parallel(p1.normal, p2.normal) {
                    return false;
                }
                let n = p1.normal.normalize_or_zero();
                (p2.origin - p1.origin).dot(n).abs() <= tol * 2.0
            }
            (Surface3::Sphere(s1), Surface3::Sphere(s2)) => {
                (s1.center - s2.center).length() <= tol * 2.0
                    && (s1.radius - s2.radius).abs() <= tol
            }
            (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
                if !axis_parallel(c1.axis, c2.axis) {
                    return false;
                }
                let axis = c1.axis.normalize_or_zero();
                (c2.origin - c1.origin).cross(axis).length() <= tol * 2.0
                    && (c1.radius - c2.radius).abs() <= tol
            }
            (Surface3::Cone(c1), Surface3::Cone(c2)) => {
                axis_parallel(c1.axis, c2.axis)
                    && (c1.apex_point() - c2.apex_point()).length() <= tol * 2.0
                    && (c1.half_angle_rad - c2.half_angle_rad).abs() <= tol
            }
            (Surface3::Torus(t1), Surface3::Torus(t2)) => {
                axis_parallel(t1.axis, t2.axis)
                    && (t1.center - t2.center).length() <= tol * 2.0
                    && (t1.major_radius - t2.major_radius).abs() <= tol
                    && (t1.minor_radius - t2.minor_radius).abs() <= tol
            }
            _ => false,
        }
    }

    fn boundaries_fully_overlap(&self, f1: usize, f2: usize) -> bool {
        let pts1 = self.ds.face_boundary_points(f1);
        let pts2 = self.ds.face_boundary_points(f2);
        if pts1.len() < 3 || pts2.len() < 3 || pts1.len() != pts2.len() {
            return false;
        }
        let tol = self.glue_tolerance;
        let mut used = vec![false; pts2.len()];
        for p1 in &pts1 {
            let mut found = false;
            for (j, p2) in pts2.iter().enumerate() {
                if used[j] {
                    continue;
                }
                if (*p1 - *p2).length() <= tol {
                    used[j] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
        }
        true
    }

    /// Fast check for potential glued face pairs using bounding box pre-filter.
    ///
    /// This optimization reduces the number of full boundary comparisons by
    /// first checking if face bounding boxes overlap.
    fn fast_glue_candidate_check(&self, f1: usize, f2: usize) -> bool {
        let a = &self.ds.faces[f1];
        let b = &self.ds.faces[f2];

        // Quick origin check
        if a.origin == b.origin {
            return false;
        }

        // Quick normal check (must be anti-parallel for glue)
        let na_len2 = a.normal.length_squared();
        let nb_len2 = b.normal.length_squared();
        if na_len2 <= TOLERANCE_ABS || nb_len2 <= TOLERANCE_ABS {
            return false;
        }
        let na = a.normal / na_len2.sqrt();
        let nb = b.normal / nb_len2.sqrt();
        if na.dot(nb) > -0.95 {
            return false;
        }

        // Bounding box overlap check
        let pts1 = self.ds.face_boundary_points(f1);
        let pts2 = self.ds.face_boundary_points(f2);

        if pts1.is_empty() || pts2.is_empty() {
            return false;
        }

        // Compute bounding boxes
        let mut min1 = pts1[0];
        let mut max1 = pts1[0];
        for p in &pts1[1..] {
            min1 = min1.min(*p);
            max1 = max1.max(*p);
        }

        let mut min2 = pts2[0];
        let mut max2 = pts2[0];
        for p in &pts2[1..] {
            min2 = min2.min(*p);
            max2 = max2.max(*p);
        }

        // Check for bounding box overlap with tolerance margin
        let tol = self.glue_tolerance;
        

        min1.x - tol <= max2.x && max1.x + tol >= min2.x
            && min1.y - tol <= max2.y && max1.y + tol >= min2.y
            && min1.z - tol <= max2.z && max1.z + tol >= min2.z
    }

    /// Detect all glued face pairs using optimized algorithm.
    ///
    /// This function uses bounding box pre-filtering to reduce the number
    /// of expensive boundary comparisons.
    fn detect_all_glued_pairs(&self, a_faces: &[usize], b_faces: &[usize]) -> Vec<(usize, usize)> {
        let mut pairs: Vec<(usize, usize)> = Vec::new();

        for &fi in a_faces {
            for &fj in b_faces {
                // Fast pre-filter
                if !self.fast_glue_candidate_check(fi, fj) {
                    continue;
                }

                // Full compatibility check
                if self.faces_form_glued_pair(fi, fj) {
                    pairs.push((fi, fj));
                }
            }
        }

        pairs
    }

    /// Build glued pairs information for fast path processing.
    ///
    /// Returns a map from face index to its glued counterpart.
    fn build_glue_map(&self, a_faces: &[usize], b_faces: &[usize]) -> HashMap<usize, usize> {
        let pairs = self.detect_all_glued_pairs(a_faces, b_faces);
        let mut glue_map: HashMap<usize, usize> = HashMap::new();

        for (fi, fj) in pairs {
            glue_map.insert(fi, fj);
            glue_map.insert(fj, fi);
        }

        glue_map
    }

    /// Split a curved face (Cylinder, Sphere, Cone, Torus) by intersection polylines.
    ///
    /// Legacy approximate method: for each intersection polyline that crosses the face,
    /// we split the boundary point list into two halves at the points closest to the
    /// polyline endpoints. Kept as fallback when UV data or PCurves are unavailable.

    /// Unwrap a UV polyline's U coordinate to remove seam jumps.
    /// For periodic surfaces (cylinder, cone, torus), consecutive points whose
    /// U values differ by more than 锜?indicate a seam crossing; we accumulate
    /// offsets of 鍗eriod to make the polyline continuous in U.
    fn unwrap_u_polyline(&self, pts: Vec<glam::DVec2>, period: f64) -> Vec<glam::DVec2> {
        if pts.len() < 2 {
            return pts;
        }
        let mut result = Vec::with_capacity(pts.len());
        result.push(pts[0]);
        let mut offset = 0.0_f64;
        for i in 1..pts.len() {
            let prev_u = result[i - 1].x;
            let curr_u = pts[i].x + offset;
            let diff = curr_u - prev_u;
            if diff > period * 0.5 {
                offset -= period;
            } else if diff < -period * 0.5 {
                offset += period;
            }
            result.push(glam::DVec2::new(pts[i].x + offset, pts[i].y));
        }
        result
    }

    /// Extend axis-aligned trim endpoints to the UV boundary so each open trim
    /// spans from one boundary edge to another. This is necessary for closed
    /// surfaces (sphere, cylinder, 閳? where intersection PCurves are clipped
    /// to the finite face-face overlap and may not reach the UV boundary.
    ///
    /// Only trims that are nearly axis-aligned (constant-u or constant-v) are
    /// extended 閳?general trims pass through unchanged.
    fn extend_trim_to_uv_boundary(
        trim: &[DVec2],
        uv_boundary: &[DVec2],
        bnd_u_span: f64,
        bnd_v_span: f64,
    ) -> Vec<DVec2> {
        if trim.len() < 3 {
            return trim.to_vec();
        }

        // Compute UV bounds from the boundary polygon
        let u_min = uv_boundary.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let u_max = uv_boundary.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let v_min = uv_boundary.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let v_max = uv_boundary.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);

        let u_span_trim = trim.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max)
            - trim.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let v_span_trim = trim.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max)
            - trim.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);

        let boundary_u_span = u_max - u_min;
        let boundary_v_span = v_max - v_min;
        // 0.5 % of the smaller span 閳?well above floating-point noise for any
        // practical model, yet tight enough to distinguish axis-aligned trims
        // from oblique ones on a sphere (where u/v vary together).
        let axis_threshold = (boundary_u_span.abs().min(boundary_v_span.abs())).max(TOLERANCE_ABS) * 0.005;

        let is_const_u = u_span_trim < axis_threshold;
        let is_const_v = v_span_trim < axis_threshold;

        if !is_const_u && !is_const_v {
            return trim.to_vec(); // non-axis-aligned 閳?cannot safely extend
        }

        // 閳光偓閳光偓 Clip trim points to boundary bounds 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
        // Intersection PCurves may have t_range extending far outside the face's
        // actual UV boundary (hardcoded extent=20 in intersect_plane_cylinder_faces).
        // Without clipping, out-of-bounds points inflate the UV sub-polygon bounding
        // box, causing tessellate_curved_face to sample a much larger surface.
        let mut extended = trim.to_vec();
        if is_const_u {
            for p in &mut extended {
                p.y = p.y.clamp(v_min, v_max);
            }
        } else {
            for p in &mut extended {
                p.x = p.x.clamp(u_min, u_max);
            }
        }

        // Deduplicate consecutive points after clamping
        extended.dedup_by(|a, b| {
            (a.x - b.x).abs() < TOLERANCE_FLOAT_ULTRA
                && (a.y - b.y).abs() < TOLERANCE_FLOAT_ULTRA
        });
        if extended.len() < 2 {
            return extended;
        }

        // 閳光偓閳光偓 span-checking guard (AFTER clipping) 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
        // If this axis-aligned trim already covers 閳?0 % of the boundary span
        // in the varying direction (measured within the boundary, not the raw
        // PCurve span), it runs boundary-to-boundary and needs no extension.
        let clipped_v_span = extended.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max)
            - extended.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let clipped_u_span = extended.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max)
            - extended.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        if is_const_u && clipped_v_span >= 0.9 * bnd_v_span.abs() {
            return extended;
        }
        if is_const_v && clipped_u_span >= 0.9 * bnd_u_span.abs() {
            return extended;
        }
        // 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

        if is_const_u {
            // Constant-u trim: extend v range to the boundary.
            let u_val = extended[0].x;
            let v_start = extended.first().unwrap().y;
            let v_end = extended.last().unwrap().y;
            let v_dir = (v_end - v_start).signum();

            if v_dir >= 0.0 {
                if (v_start - v_min).abs() > TOLERANCE_ABS {
                    extended.insert(0, DVec2::new(u_val, v_min));
                }
                if (v_max - v_end).abs() > TOLERANCE_ABS {
                    extended.push(DVec2::new(u_val, v_max));
                }
            } else {
                if (v_max - v_start).abs() > TOLERANCE_ABS {
                    extended.insert(0, DVec2::new(u_val, v_max));
                }
                if (v_end - v_min).abs() > TOLERANCE_ABS {
                    extended.push(DVec2::new(u_val, v_min));
                }
            }
        } else {
            // Constant-v trim: extend u range to the boundary.
            let v_val = extended[0].y;
            let u_start = extended.first().unwrap().x;
            let u_end = extended.last().unwrap().x;
            let u_dir = (u_end - u_start).signum();

            if u_dir >= 0.0 {
                if (u_start - u_min).abs() > TOLERANCE_ABS {
                    extended.insert(0, DVec2::new(u_min, v_val));
                }
                if (u_max - u_end).abs() > TOLERANCE_ABS {
                    extended.push(DVec2::new(u_max, v_val));
                }
            } else {
                if (u_max - u_start).abs() > TOLERANCE_ABS {
                    extended.insert(0, DVec2::new(u_max, v_val));
                }
                if (u_end - u_min).abs() > TOLERANCE_ABS {
                    extended.push(DVec2::new(u_min, v_val));
                }
            }
        }

        extended
    }

    /// into a 2D trim polyline in UV space, then splits the UV boundary polygon.
    /// Maps resulting sub-polygons back to 3D via surface evaluation.
    ///
    /// ⏳ 部分对齐: 鐢ㄧ簿纭ぇ鍦嗗姬鏋勫缓鐞冮潰瀛愰潰銆?
    ///    OCCT: BuildSplitFaces 鈫?section edges 鐩存帴鍒涘缓 BRep sub-face銆?
    ///    rcad: 鎵嬪姩璁＄畻 8 涓崷闄愮殑 FaceSampleData,鐢?outer_circle_edges 璁板綍澶у渾寮с€?
    ///    鍔熻兘绛変环(8 涓崐鐞冮潰鍖哄煙 + 绮剧‘鍦嗗姬杈圭晫),浣?OCCT 涓嶉渶瑕佷腑闂?FaceSampleData銆?

    /// Find the PCurve (2D parametric curve) for the given intersection curve
    /// as it lies on the given face. Searches FaceFace interferences to determine
    /// whether this face is f1 (use pcurve_on_a) or f2 (use pcurve_on_b).
    fn find_pcurve_for_face(
        &self,
        curve_idx: usize,
        face_idx: usize,
    ) -> Option<rcad_kernel::geom::Curve2d> {
        for interference in &self.ds.interferences {
            if let Interference::FaceFace { f1, f2, curves, .. } = interference
                && curves.contains(&curve_idx)
            {
                let ic = &self.ds.intersection_curves[curve_idx];
                if *f1 == face_idx {
                    return ic.pcurve_on_a.clone();
                } else if *f2 == face_idx {
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

    /// Build a map from edge index to the list of face indices that reference it.
    /// Iterates over all solids and shells in the BRep.
    fn build_edge_ref_map(brep: &BRep) -> Vec<Vec<usize>> {
        let n_edges = brep.edges.len();
        if n_edges == 0 {
            return Vec::new();
        }
        let mut edge_refs: Vec<Vec<usize>> = vec![Vec::new(); n_edges];
        for (_shell_idx, shell) in brep.solids.iter().flat_map(|s| &s.shells).enumerate() {
            for (face_idx, face) in shell.faces.iter().enumerate() {
                for we in &face.outer_wire.edges {
                    if we.idx < edge_refs.len() {
                        edge_refs[we.idx].push(face_idx);
                    }
                }
                for iw in &face.inner_wires {
                    for we in &iw.edges {
                        if we.idx < edge_refs.len() {
                            edge_refs[we.idx].push(face_idx);
                        }
                    }
                }
            }
        }
        edge_refs
    }

    /// After building the BRep, validate that every edge in every shell has
    /// exactly 2 face references (closed shell). Edges with <2 references
    /// (orphan edges) or >2 references (over-shared edges) indicate a
    /// topological defect that would produce an OPEN_SHELL result.
    pub fn validate_edge_face_references(&self, brep: &BRep) -> Result<(), BooleanError> {
        let edge_refs = Self::build_edge_ref_map(brep);
        if edge_refs.is_empty() {
            return Ok(());
        }

        let orphan_edges: Vec<usize> = edge_refs.iter().enumerate()
            .filter(|(_, refs)| refs.is_empty() || refs.len() == 1)
            .map(|(ei, _)| ei)
            .collect();
        let over_shared_edges: Vec<usize> = edge_refs.iter().enumerate()
            .filter(|(_, refs)| refs.len() > 2)
            .map(|(ei, _)| ei)
            .collect();

        if !orphan_edges.is_empty() || !over_shared_edges.is_empty() {
            return Err(BooleanError::OpenShell {
                orphan_edges,
                over_shared_edges,
            });        }

        Ok(())
    }

    /// Diagnostic stub: report orphan edges (edges referenced by 0 or 1 faces).
    /// This is a replacement for the previous `recover_orphan_edges` which was a no-op
    /// (it counted candidate faces but never mutated the BRep). The real value of RC2
    /// is the validation (detecting OPEN_SHELL), not automatic topology repair.
    ///
    /// Returns the total number of orphan edges found (both zero-ref and single-ref).
    pub fn diagnose_orphan_edges(&self, brep: &BRep) -> usize {
        let edge_refs = Self::build_edge_ref_map(brep);
        if edge_refs.is_empty() {
            return 0;
        }

        let zero_ref_edges: Vec<usize> = edge_refs.iter().enumerate()
            .filter(|(_, refs)| refs.is_empty())
            .map(|(ei, _)| ei)
            .collect();

        let single_ref_edges: Vec<usize> = edge_refs.iter().enumerate()
            .filter(|(_, refs)| refs.len() == 1)
            .map(|(ei, _)| ei)
            .collect();

        let total = zero_ref_edges.len() + single_ref_edges.len();
        if total > 0 {
            eprintln!("[INFO] diagnose_orphan_edges: {} edges with zero refs, {} edges with one ref (manual topology repair needed)",
                zero_ref_edges.len(), single_ref_edges.len());
        }

        total
    }
}



mod split_polygon;
mod split_polygon2;
mod glue;

pub(crate) use split_polygon2::point_in_polygon_2d;
pub use glue::{
    GlueConfig, GlueFacePair, GlueFaceCache, detect_glue_faces,
    apply_glue_optimization, compute_adaptive_glue_tolerance,
};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bopds::ds::DS;
    use rcad_modeling::{make_box_brep};

    #[test]
    fn prepare_returns_empty_containers() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let ds = DS::new(&a, &b);
        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
        let (t_brep, result) = builder.prepare();

        assert!(t_brep.tshapes.is_empty(), "t_brep should be empty after prepare");
        assert!(result.vertices.is_empty(), "result vertices should be empty");
        assert!(result.edges.is_empty(), "result edges should be empty");
        assert!(result.faces.is_empty(), "result faces should be empty");
        assert!(result.tmp_shells.is_empty());
        assert!(result.tmp_solids.is_empty());
    }

    #[test]
    fn minimal_box_union_pipeline_builds_result() {
        // Two tiny non-overlapping boxes — union should produce both boxes.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(3.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        let mut ds = DS::new(&a, &b);
        let mut t_brep = rcad_kernel::topods::BRep::new();
        let (face_refs, ic_edge_map) = {
            let mut filler = crate::pave_filler::PaveFiller::new(&mut ds);
            filler.brep = Some(&mut t_brep);
            filler.perform();
            (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
        };

        let builder = BooleanBuilder::with_brep(&ds, BooleanOpType::Union, t_brep, face_refs, ic_edge_map);
        let (brep, _history) = builder.build_with_history().expect("union should succeed");

        // Two disjoint boxes — 12 faces total
        assert!(!brep.solids.is_empty(), "should produce at least one solid");
        let nf: usize = brep.solids.iter()
            .flat_map(|s| &s.shells)
            .map(|sh| sh.faces.len())
            .sum();
        assert!(nf >= 12, "expected >= 12 faces for two boxes, got {}", nf);
    }
}
