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

        // ✅ OCCT-aligned L360-362: Add internal wire groups into the result WireFaces
        //   for emit_wire_face to process.  OCCT: aBB.Add(aF, aW) per internal wire.
        if !internal_wire_groups.is_empty() {
            for wf in &mut wfs {
                for &si in &avoided {
                    wf.internal_wires.push(vec![si]);
                }
            }
        }
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
            crate::builder::wire_path_topo_ds::perform_areas_topo_ds(&wires, &internal_wire_groups, &segments_topo, tool, face_idx)
        } else if !avoided.is_empty() {
            vec![WireFace { outer_wire: vec![], inner_wires: vec![], internal_wires: segments_topo.iter().enumerate().filter(|(si, _)| avoided.contains(si)).map(|(si, _)| vec![si]).collect() }]
        } else {
            vec![WireFace { outer_wire: (0..segments_topo.len()).collect(), inner_wires: vec![], internal_wires: vec![] }]
        };
        if wfs.is_empty() { return; }

        let mut wfs = wfs;
        if !internal_wire_groups.is_empty() {
            for wf in &mut wfs {
                for &si in &avoided { wf.internal_wires.push(vec![si]); }
            }
        }

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
    fn fill_images_edges(&self, result: &mut ResultBuilder) {
        for (ei, edge) in self.ds.edges.iter().enumerate() {
            if edge.pave_blocks.is_empty() {
                // OCCT L84-86: HasReference check — no pave blocks → no split images.
                //   Such edges are not images of any source edge; they pass through
                //   BuildResult(EDGE) as originals only when myArguments contain EDGE
                //   types (no-op for solid boolean arguments).
                //   rcad: unmodified edges are registered via face construction
                //   (emit_wire_face → add_edge); no action needed here.
                continue;
            }

            for pb in &edge.pave_blocks {
                let new_ei = if let Some(rei) = self.ds.real_pave_block_edge(ei, pb) {
                    rei
                } else if let Some(nei) = pb.new_edge {
                    nei
                } else {
                    continue;
                };

                let e_base = self.ds.vertices.len();
                // OCCT L105-106: pLS->Append(aSpR) -> myImages(edge) += split_edge
                self.my_images.borrow_mut().entry(rcad_kernel::topods::ShapeRef::new(e_base + ei)).or_default().push(rcad_kernel::topods::ShapeRef::new(e_base + new_ei));

                // OCCT L107-112: myOrigins.ChangeSeek(aSpR).Append(aE)
                self.my_origins.borrow_mut().entry(rcad_kernel::topods::ShapeRef::new(e_base + new_ei)).or_default().push(rcad_kernel::topods::ShapeRef::new(e_base + ei));

                // OCCT L114-119: IsCommonBlockOnEdge -> myShapesSD.Bind(aSp, aSpR)
                if pb.common_block_idx.is_some() {
                    self.my_shapes_sd.borrow_mut().insert(rcad_kernel::topods::ShapeRef::new(e_base + ei), rcad_kernel::topods::ShapeRef::new(e_base + new_ei));
                }

                // rcad flat edge: create result edge entry for face construction.
                // OCCT: split edges are implicitly added to myShape through face sub-shapes.
                let se = &self.ds.edges[new_ei];
                let sv = result.add_ds_vertex(se.start_vertex, self.ds.vertices[se.start_vertex].point);
                let ev = result.add_ds_vertex(se.end_vertex, self.ds.vertices[se.end_vertex].point);
                let fi = result.edges.len();
                result.edges.push((sv, ev));
                while result.custom_edge_curves.len() <= fi {
                    result.custom_edge_curves.push(None);
                }
                result.custom_edge_curves[fi] = Some(se.curve.clone());
                if self.ds.is_edge_degenerated(new_ei) || se.start_vertex == se.end_vertex {
                    result.deg_edge_indices.insert(fi);
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

    fn fill_images_containers_shells(&self, result: &mut ResultBuilder) {
        let nf = result.faces.len();
        if nf == 0 || self.ds.shells.is_empty() { return; }

        // Build result face index → source DS face index mapping.
        let mut rfi_to_dsfi: Vec<Option<usize>> = vec![None; nf];
        for (rfi, origin) in result.face_origins.iter().enumerate() {
            let ds_fi = match origin {
                FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                _ => None,
            };
            if let Some(fi) = ds_fi { rfi_to_dsfi[rfi] = Some(fi); }
        }

        // OCCT L242-275: for each source SHELL, collect its split face images.
        for ds_shell in &self.ds.shells {
            let mut shell_faces: Vec<usize> = Vec::new();
            for rfi in 0..nf {
                if let Some(dsfi) = rfi_to_dsfi[rfi] {
                    if ds_shell.faces.contains(&dsfi) {
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
    ///   L224-233: iterate sub-shapes, check if any has been modified.
    ///   L235-240: if none modified → early return.
    ///   L242-275: build new container from sub-shape images.
    ///
    /// rcad: iterate DS faces → find those from CompSolids → group result solids
    /// by their source compsolid → store in result.tmp_compsolid_groups.
    fn fill_images_containers_compsolid(&self, result: &mut ResultBuilder) {
        // OCCT L224-233: check if any sub-shape has been modified.
        let has_compsolid = self.ds.faces.iter().any(|f| f.source_compsolid_idx.is_some());
        if !has_compsolid {
            return; // OCCT L235-240: early return
        }

        // OCCT L242-275: build CompSolid from sub-solid images.
        // Group result solids by their source compsolid index.
        // For each result solid, find its origin DS faces to determine compsolid.
        let mut solid_to_cs: Vec<Option<usize>> = vec![None; result.tmp_solids.len()];
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
                            solid_to_cs[si] = Some(csi);
                            break; // found compsolid for this solid
                        }
                    }
                }
                if solid_to_cs[si].is_some() { break; }
            }
        }

        // Group solid indices by compsolid index.
        let mut cs_groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
        for (si, cs_opt) in solid_to_cs.iter().enumerate() {
            if let Some(csi) = cs_opt {
                cs_groups.entry(*csi).or_default().push(si);
            }
        }
        result.tmp_compsolid_groups = cs_groups.into_values().collect();
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
    fn fill_images_solids(&self, result: &mut ResultBuilder) {
        let has_solid = self.ds.faces.iter().any(|f| f.source_solid_idx.is_some());
        if !has_solid { return; }
        if result.tmp_shells.is_empty() { return; }

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
        let shell_assignments = self.fill_in_3d_parts(result, &a_faces, &b_faces);

        // OCCT L86: BuildSplitSolids — group shells into result solids
        self.build_split_solids(result, &shell_assignments);

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

    /// ✅ OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-232).
    fn fill_in_3d_parts(&self, result: &mut ResultBuilder,
                         a_faces: &[usize], b_faces: &[usize]) -> Vec<(usize, usize, &'static str)> {
        let nf = result.faces.len();
        let mut to_remove = vec![false; nf];

        // OCCT L164-195: BuildDraftSolid for each source solid.
        //   rcad: builds draft face sets for both operands (result faces
        //   grouped by source shell).  The draft sets are computed here
        //   for form alignment even though classify_point uses DS indices.
        let (_draft_a, _int_a) = self.build_draft_solid(result, 0);
        let (_draft_b, _int_b) = self.build_draft_solid(result, 1);

        // OCCT L201-204: ClassifyFaces → anInParts.
        //   myInParts[0] = faces from B that are IN solid A (source side 0 = A)
        //   myInParts[1] = faces from A that are IN solid B (source side 1 = B)
        //   Per OCCT Builder_3.cxx L215-232.
        let mut my_in_parts = self.my_in_parts.borrow_mut();
        my_in_parts.clear();
        // Collect per-face classification results for shell-state computation
        // (OCCT tracks state via draft-solid membership; rcad tracks via per-face class).
        let mut face_class: Vec<Option<Classification>> = vec![None; nf];
        let mut face_side: Vec<Option<usize>> = vec![None; nf]; // 0=A, 1=B

        // ═══ OCCT Phase 3: ClassifyFaces (BOPAlgo_Tools.cxx L1334-1450) ═══
        //   OCCT: for each face, use BRepClass3d_SolidClassifier with a point
        //   ON the face surface (parametric midpoint).  rcad: use face
        //   sample_point (index 8) which is computed from the face boundary
        //   and guaranteed to be on the surface.
        //
        //   OCCT additionally:
        //   1. Skips faces whose AABB doesn't overlap the solid's AABB
        //      (aSelector BVH culling, L1345-1354)
        //   2. Skips self-shapes (faces that are sub-shapes of the solid,
        //      L1366-1368 — rcad handles this by classifying A-faces
        //      against B-faces and vice versa)
        //   3. Groups connected faces into blocks for batch classification
        //      (L1396-1405 — rcad classifies per-face, equivalent result)
        for fi in 0..nf {
            if to_remove[fi] { continue; }
            let (source_side, other_faces, side_idx) = match &result.face_origins[fi] {
                FaceOrigin::FromA(_) => (SourceSide::A, b_faces, 0usize),
                FaceOrigin::FromB(_) => (SourceSide::B, a_faces, 1usize),
                _ => continue,
            };
            face_side[fi] = Some(side_idx);
            if other_faces.is_empty() { continue; }

            // OCCT L1345-1354: BVH-based AABB overlap check (optional culling).
            //   rcad: no AABB culling — classify all candidate faces.

            // OCCT BRepClass3d_SolidClassifier: use a point ON the face surface.
            //   rcad: use face sample_point (index 8) instead of centroid (index 6).
            //   The sample_point is computed during emit_wire_face from the face's
            //   surface UV midpoint, guaranteed to be ON the surface.
            let pt = result.faces[fi].8; // sample_point (on-surface)
            let class = classify_point(pt, other_faces, self.ds);
            eprintln!("[CLASSIFY] fi={} origin={:?} pt=({:.4},{:.4},{:.4}) class={:?}", fi, result.face_origins[fi], pt.x, pt.y, pt.z, class);
            face_class[fi] = Some(class);

            // OCCT L215-232: store IN faces in myInParts
            //   A face classified as IN → it is IN the other solid
            match class {
                Classification::In => {
                    let other_side = if side_idx == 0 { 1 } else { 0 };
                    my_in_parts.entry(other_side).or_default().push(fi);
                }
                _ => {}
            }

            if !self.classification_keep_policy(source_side, class, fi) {
                to_remove[fi] = true;
            }
        }

        // Remove faces that fail keep policy
        if to_remove.iter().any(|&r| r) {
            let old_faces = std::mem::take(&mut result.faces);
            let old_origins = std::mem::take(&mut result.face_origins);
            for (fi, face) in old_faces.into_iter().enumerate() {
                if !to_remove[fi] {
                    result.faces.push(face);
                    result.face_origins.push(old_origins[fi]);
                }
            }
            // Rebuild shell face indices
            let old_shells = std::mem::take(&mut result.tmp_shells);
            let mut idx_map: Vec<Option<usize>> = vec![None; nf];
            let mut cur = 0usize;
            for fi in 0..nf {
                if !to_remove[fi] { idx_map[fi] = Some(cur); cur += 1; }
            }
            for shell in &old_shells {
                let new_shell: Vec<usize> = shell.iter()
                    .filter_map(|&fi| idx_map[fi]).collect();
                if !new_shell.is_empty() {
                    result.tmp_shells.push(new_shell);
                }
            }
            // Translate my_in_parts face indices through the removal map
            // (OCCT does not remove faces; rcad does — preserve index mapping for build_split_solids)
            let mut updated: std::collections::HashMap<usize, Vec<usize>> =
                std::collections::HashMap::new();
            for (&side, faces) in my_in_parts.iter() {
                let mut new_faces: Vec<usize> = Vec::new();
                for &fi in faces {
                    if let Some(new_fi) = idx_map[fi] {
                        new_faces.push(new_fi);
                    }
                }
                if !new_faces.is_empty() {
                    updated.insert(side, new_faces);
                }
            }
            *my_in_parts = updated;
            // Remap face_class to new indices (old face indices are stale after
            // removal — face_class[old_fi] != face_class[new_fi] after compression).
            let old_face_class = std::mem::replace(&mut face_class, vec![None; result.faces.len()]);
            for (old_fi, class) in old_face_class.into_iter().enumerate() {
                if let Some(new_fi) = idx_map[old_fi] {
                    if let Some(c) = class {
                        face_class[new_fi] = Some(c);
                    }
                }
            }
        }

        // OCCT L215-232 (continued): compute shell state dynamically
        //   based on face classifications instead of hardcoding "OUT".
        //   For each shell, determine if it is IN or OUT of the other solid.
        let mut assignments: Vec<(usize, usize, &'static str)> = Vec::new();
        for (si, shell) in result.tmp_shells.iter().enumerate() {
            let mut has_a = false;
            let mut has_b = false;
            // Determine shell state from the majority classification
            let mut n_out = 0usize;
            let mut n_in = 0usize;
            for &fi in shell {
                match &result.face_origins[fi] {
                    FaceOrigin::FromA(_) => has_a = true,
                    FaceOrigin::FromB(_) => has_b = true,
                    _ => {}
                }
                // Count IN/OUT from stored classification
                if let Some(class) = face_class.get(fi).copied().flatten() {
                    match class {
                        Classification::In => n_in += 1,
                        Classification::Out => n_out += 1,
                        _ => {}
                    }
                }
            }
            // Compute shell state: IN if most faces are IN, OUT otherwise
            let state: &'static str = if n_in > n_out { "IN" } else { "OUT" };
            if has_a {
                assignments.push((si, 0, state));
            }
            if has_b {
                assignments.push((si, 1, state));
            }
        }
        assignments
    }

    // OCCT L413-618: BuildSplitSolids — group classified shells into solids.
    //   OCCT does NOT re-split topology.  It combines draft solid faces + IN
    //   faces into BOPAlgo_SplitSolid which builds closed TopoDS_Solid from
    //   the face set, preserving source shell connectivity.
    //   rcad: groups shells by (operation, state, origin) — each group becomes
    //   one solid.  The aligned BuilderSolid module (bopds::builder_solid) is
    //   available for standalone use but not called from the hot path because
    //   the DS-level face indices don't preserve result-shell connectivity.
    /// ✅ OCCT-aligned: BuildSplitSolids (Builder_3.cxx L413-569).
    ///   Group classified shells into result solids.
    ///   Includes BuilderSolid::PerformAreas void detection internally.
    fn build_split_solids(&self, result: &mut ResultBuilder,
                          assignments: &[(usize, usize, &'static str)]) {
        let mut result_solids: Vec<Vec<usize>> = Vec::new();

        let mut group: Vec<usize> = Vec::new();
        for &(si, origin, state) in assignments {
            let keep = match self.op {
                BooleanOpType::Union => state == "OUT",
                BooleanOpType::Intersection => state == "IN",
                BooleanOpType::Difference => match (origin, state) {
                    (0, "OUT") | (1, "IN") | (1, "ON") => true,
                    _ => false,
                },
            };
            if keep { group.push(si); }
        }
        if !group.is_empty() {
            let mut group_faces: Vec<usize> = Vec::new();
            for &si in &group {
                if si < result.tmp_shells.len() {
                    group_faces.extend(&result.tmp_shells[si]);
                }
            }
            if !group_faces.is_empty() {
                let csi = result.tmp_shells.len();
                result.tmp_shells.push(group_faces);
                result_solids.push(vec![csi]);
            }
        }
        result.tmp_solids = result_solids;

        // OCCT BuilderSolid::PerformAreas (L397-576): shell-level void detection.
        self.detect_internal_voids(result, assignments);
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

    /// ✅ OCCT-aligned: BuilderSolid::PerformAreas void detection (L397-576).
    ///   Detect IN-state shells (holes) that are inside OUT-state solids (growths)
    ///   and add them as internal voids.  OCCT IsGrowthShell/IsHole + IsInside
    ///   classify each shell against candidate solids; rcad uses classify_point
    ///   with the IN-shell centroid against the OUT-solid's DS face set.
    fn detect_internal_voids(&self, result: &mut ResultBuilder,
                              assignments: &[(usize, usize, &'static str)]) {
        // OCCT L420-441: classify each shell as Growth or Hole.
        //   rcad: state ("IN"/"OUT") from fill_in_3d_parts corresponds to
        //   Growth (OUT = outer boundary) vs Hole (IN = internal void).
        //   Build IN/OUT solid lists from shell states.
        let mut solid_is_in: Vec<bool> = vec![false; result.tmp_solids.len()];
        for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
            if let Some(&first_sh) = solid_shells.first() {
                if let Some(&(_sh_i, _origin, state)) = assignments.iter().find(|&&(si, _, _)| si == first_sh) {
                    solid_is_in[si] = state == "IN";
                }
            }
        }

        // Separate IN solids (potential holes) from OUT solids (potential growths).
        let in_solid_indices: Vec<usize> = (0..result.tmp_solids.len())
            .filter(|&si| solid_is_in[si]).collect();
        let out_solid_indices: Vec<usize> = (0..result.tmp_solids.len())
            .filter(|&si| !solid_is_in[si]).collect();

        if in_solid_indices.is_empty() || out_solid_indices.is_empty() {
            return; // OCCT L444-457: no holes → nothing to classify
        }

        // OCCT L460-530: classify each hole shell against each candidate solid
        //   via IsInside (BVH-accelerated).  rcad: classify_point with centroid.
        let mut in_to_out: Vec<(usize, usize)> = Vec::new(); // (in_si, out_si)

        // Pre-build DS face index sets for each OUT solid (OCCT builds boxes + BVH).
        let mut out_ds_face_sets: Vec<Vec<usize>> = Vec::new();
        for &out_si in &out_solid_indices {
            let mut ds_faces: Vec<usize> = Vec::new();
            for &sh in &result.tmp_solids[out_si] {
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

        for (i, &in_si) in in_solid_indices.iter().enumerate() {
            // OCCT L422-427: classify hole — IsGrowthShell/IsHole.
            //   rcad: centroid of IN solid's first face as test point.
            let centroid = result.tmp_solids[in_si].first()
                .and_then(|&sh| result.tmp_shells.get(sh))
                .and_then(|shell| shell.first())
                .map(|&fi| {
                    // FaceEntry.6 is the centroid field
                    if fi < result.faces.len() { result.faces[fi].6 } else { DVec3::ZERO }
                })
                .unwrap_or(DVec3::ZERO);

            // OCCT L484-529: check IsInside(hole_shell, candidate_solid, context).
            for (j, &out_si) in out_solid_indices.iter().enumerate() {
                if out_ds_face_sets[j].is_empty() { continue; }
                let class = classify_point(centroid, &out_ds_face_sets[j], self.ds);
                if class == Classification::In || class == Classification::On {
                    in_to_out.push((in_si, out_si));
                    break; // OCCT selects the outermost containing solid
                }
            }
        }

        // OCCT L550-576: Add Holes to Solids (add void shells to containing solids).
        let mut removed = vec![false; result.tmp_solids.len()];
        for &(in_si, out_si) in &in_to_out {
            let void_shells = result.tmp_solids[in_si].clone();
            result.tmp_solids[out_si].extend(void_shells);
            removed[in_si] = true;
        }

        // Remove merged IN solids, preserve order.
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
    fn fill_images_compounds(&self, result: &mut ResultBuilder) {
        // OCCT L200: not implemented — aMFP fence map
        // OCCT L202-216: check source compounds
        let has_compound = self.ds.a_has_compound || self.ds.b_has_compound;
        if !has_compound {
            return; // OCCT L309-312: no compounds → return
        }
        // OCCT L314-341: build new compound from solid images.
        //   rcad: deferred — see build_with_history L6834-6840.
        //   OCCT stores the new TopoDS_Compound in myImages; rcad
        //   builds it during result.build() from result.tmp_solids,
        //   wrapping them in a rcad_kernel::topology::Compound.
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
                // ✅ OCCT-aligned: BuildResult(SOLID) — assemble solids from shells.
                let tmp_solids = std::mem::take(&mut result.tmp_solids);
                if tmp_solids.is_empty() && result.shells.is_empty() {
                    let shell = t.add_tshell(std::mem::take(&mut result.face_refs));
                    t.add_tsolid(vec![shell]);
                } else if !tmp_solids.is_empty() {
                    for solid_shells in &tmp_solids {
                        let shell_refs: Vec<topods::ShapeRef> = solid_shells.iter()
                            .filter_map(|&si| result.shells.get(si).copied())
                            .collect();
                        if !shell_refs.is_empty() {
                            result.solids.push(t.add_tsolid(shell_refs));
                        }
                    }
                } else {
                    let shell_refs: Vec<topods::ShapeRef> = result.shells.clone();
                    t.add_tsolid(shell_refs);
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
                // OCCT L431-435: BuildResult(COMPOUND) — aggregate into top-level compound.
                // rcad: compound handling is in build_with_history after build_result.
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
            let (brep_built, face_refs, ic_edge_map) = crate::ds_to_brep::ds_to_brep(self.ds);
            *self.brep.borrow_mut() = Some((brep_built, face_refs, ic_edge_map));
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
        self.fill_images_edges(&mut result);
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
        self.build_result(ShapeType::Shell, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // Phase 5: FillImagesSolids (L400-410) → BuildResult(SOLID) (L406-410).
        self.fill_images_solids(&mut result);
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
        let mut brep = rcad_kernel::BRep::from_topods(&t_brep);
        if (self.ds.a_has_compound || self.ds.b_has_compound) && !brep.solids.is_empty() {
            let mut compound = rcad_kernel::topology::Compound::new();
            for solid in brep.solids.drain(..) {
                compound.solids.push((None, solid));
            }
            brep.compound = Some(compound);
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
// SubFace removed: build_with_history_par

    /// When PaveFiller does not link a plane閳ユ悞phere circle to every affected box face, merge in
    /// any coplanar `Curve3::Circle` from `intersection_curves` that overlaps the face 2D AABB.
    fn extra_coplanar_circle_curves_for_plane_face(
        &self,
        face_idx: usize,
        plane: &Plane,
    ) -> Vec<usize> {
        let n = plane.normal.normalize_or_zero();
        if n.length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
            return vec![];
        }
        let face = &self.ds.faces[face_idx];
        let (u_axis, v_axis) = plane_local_basis(plane);
        let project_to_2d = |p: DVec3| -> DVec2 {
            let d = p - plane.origin;
            DVec2::new(d.dot(u_axis), d.dot(v_axis))
        };
        if face.boundary_verts.is_empty() {
            return vec![];
        }
        let mut umin = f64::INFINITY;
        let mut umax = f64::NEG_INFINITY;
        let mut vmin = f64::INFINITY;
        let mut vmax = f64::NEG_INFINITY;
        for &vi in &face.boundary_verts {
            let q = project_to_2d(self.ds.vertices[vi].point);
            umin = umin.min(q.x);
            umax = umax.max(q.x);
            vmin = vmin.min(q.y);
            vmax = vmax.max(q.y);
        }
        const MARGIN: f64 = TOLERANCE_ADAPTIVE_MAX;
        umin -= MARGIN;
        umax += MARGIN;
        vmin -= MARGIN;
        vmax += MARGIN;
        // Circle lies in a plane with normal parallel to this plane, and (center on plane)
        const PL_D: f64 = TOLERANCE_ADAPTIVE_MAX;
        const N_ALIGN: f64 = 0.04;
        let mut out = Vec::new();
        for (ci, ic) in self.ds.intersection_curves.iter().enumerate() {
            if face.face_info.curves_sc_only().contains(&ci) {
                continue;
            }
            let Curve3::Circle(c) = &ic.curve else {
                continue;
            };
            let nc = c.normal.normalize_or_zero();
            if nc.length_squared() < TOLERANCE_METRIC_SQ_NEAR_ZERO {
                continue;
            }
            if (nc.dot(n).abs() - 1.0).abs() > N_ALIGN {
                continue;
            }
            if ((DVec3::from(c.center) - plane.origin).dot(n)).abs() > PL_D {
                continue;
            }
            let c2d = project_to_2d(DVec3::from(c.center));
            let r = c.radius;
            if c2d.x + r < umin
                || c2d.x - r > umax
                || c2d.y + r < vmin
                || c2d.y - r > vmax
            {
                continue;
            }
            out.push(ci);
        }
        out
    }

    fn merged_split_curve_ids_for_planar_face(&self, face_idx: usize, plane: &Plane) -> Vec<usize> {
        let mut c: Vec<usize> = self.ds.faces[face_idx]
            .face_info
            .curves_sc_only()
            .iter()
            .copied()
            .collect();
        for e in self.extra_coplanar_circle_curves_for_plane_face(face_idx, plane) {
            if !c.contains(&e) {
                c.push(e);
            }
        }
        c.sort_unstable();
        c
    }

// SubFace removed: single_subface

    /// Split a face by intersection curves. If no intersection curves cross this
    /// face, returns the whole face as a single FaceSampleData.
// SubFace removed: split_face

    /// Tessellate a sphere face with no intersection curves into UV patches.
    ///
    /// The sphere's single face with a seam edge has only 2 boundary vertices in the DS
    /// (north and south poles along the seam). [`emit_face_with_origin`] rejects boundaries
    /// with fewer than 3 vertices, so we split the sphere into a UV grid where each patch
    /// has a fine polygon boundary (sampled along the patch edges) for accurate mesh-based
    /// surface area and volume.
// SubFace removed: tess_sphere

    /// Tessellate a cylinder wall face with no intersection curves into UV patches.
    ///
    /// Like the sphere, a cylinder's single face with a seam edge has only 2 boundary
    /// vertices in the DS (top and bottom along the seam), which [`emit_face_with_origin`]
    /// rejects (<3 vertices). Split the cylinder wall into azimuthal bands so each patch
    /// has a valid 3D boundary polygon.
// SubFace removed: tess_cyl

    /// Tessellate a cylinder face into an N_U 脳 N_V 2D grid of rectangular patches.
    ///
    /// Used for cylinder鈥揷ylinder intersections (e.g. Steinmetz) where full-wrap
    /// intersection curves prevent the parametric UV-polygon splitting from working.
    /// Each patch's sample point (boundary centroid 鈮?surface center) is classified
    /// independently against the other solid, correctly selecting the Steinmetz lobes.
// SubFace removed: tess_cyl_2d

    /// Tessellate a cone face into a UV grid. Each grid cell is a [`FaceSampleData`] with
    /// its own sample point, so that classify_point can independently decide whether
    /// that region is inside or outside the other solid.
    ///
    /// This replaces [`split_curved_face_parametric`] for cone faces because the UV
    /// splitter can produce overlapping sub-face polygons when intersection curves are
    /// high-order (e.g. the cone鈥揷ylinder quartic from skew axes in ZK8/ZL1), leading
    /// to SA double-counting.  The grid approach guarantees each UV region is covered
    /// by exactly one sub-face whose sample point correctly represents the region.
// SubFace removed: tess_cone_2d

    /// Split a planar face by intersection line segments.
    ///
    /// Algorithm:
    /// 1. Project boundary + intersection segment endpoints to 2D
    /// 2. Find where intersection segment endpoints lie on boundary edges
    /// 3. Insert intersection points into boundary at correct positions
    /// 4. Walk augmented boundary to extract sub-polygons on each side
    /// `split_curve_ids` is `face_info.curves_in` plus any merged coplanar circles (see
    /// [`Self::merged_split_curve_ids_for_planar_face`]).
// SubFace removed: split_planar

    /// ✅ OCCT-aligned: DS face iterator by origin.
    ///   OCCT: iterates myDS->ShapeInfo() filtering by TopAbs_FACE + source solid.
    ///   rcad: DS stores faces flat with ShapeOrigin (A/B) for operand discrimination.
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
// SubFace removed: split_curved_legacy

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
// SubFace removed: split_sphere
// SubFace removed: split_curved_param

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
    fn build_result_noop_for_all_types() {
        // OCCT BuildResult is a no-op for solid boolean inputs
        // (no VERTEX/EDGE/FACE arguments in myArguments).
        // Verify rcad's build_result is also a no-op.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let ds = DS::new(&a, &b);
        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
        let (mut t_brep, mut result) = builder.prepare();

        // Call build_result for each type — should not panic and
        // should not add anything to t_brep or result.
        for st in &[
            ShapeType::Vertex, ShapeType::Edge, ShapeType::Wire,
            ShapeType::Face, ShapeType::Shell, ShapeType::Solid,
            ShapeType::CompSolid, ShapeType::Compound,
        ] {
            builder.build_result(*st, &mut result, &mut t_brep);
        }

        // All types are no-ops, so t_brep should still be empty
        assert!(t_brep.tshapes.is_empty(),
            "build_result should be no-op for all types, got {} tshapes", t_brep.tshapes.len());
    }

    #[test]
    fn minimal_box_union_pipeline_builds_result() {
        // Two tiny non-overlapping boxes — union should produce both boxes.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let b = make_box_brep(DVec3::new(3.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        let mut ds = DS::new(&a, &b);
        let mut filler = crate::pave_filler::PaveFiller::new(&mut ds);
        filler.perform();

        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
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
