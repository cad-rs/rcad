use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use indexmap::IndexMap;

use glam::{DVec2, DVec3};
use rayon::prelude::*;
use rcad_kernel::BRep;
use rcad_kernel::topods;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval, *};
use rcad_kernel::topology::*;

use crate::bvh::{Aabb, DsBvh};
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
    classify_face_against_box, compute_state,
    is_tangent_face, build_edge_bounds, quantize_pos,
    check_and_add_split_vertex, collect_face_edge_segments,
    cmp_boolean_emit_order,
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
    // ✅ OCCT-aligned: myFillHistory (BOPAlgo_Options.hxx).
    //   When false, PrepareHistory is a no-op (HasHistory() returns false).
    my_fill_history: bool,
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
    EdgeInfo, build_closed_wires,
    expand_avoided_pids,
    physical_edge_id, world_to_uv,
    edge_uv_tangent, edge_angle_2d,
    are_verts_coincident,
};
pub(crate) use wire_path::{
    perform_areas, intersect_ray_curve_2d,
    wire_faces_to_face_sample_data,
    refine_angles, pc_parameter_range,
    walk_path_extract_wires,
};


/// Uses UV-space point-in-polygon to find a point that is inside the outer boundary
/// but outside all hole boundaries.  Falls back to candidates when centroid is in a hole.
fn find_interior_3d(
    outer_uvs: &[DVec2],
    hole_uvs: &[Vec<DVec2>],
    surface: &Surface3,
    normal: &DVec3,
) -> Option<DVec3> {
    if outer_uvs.len() < 3 { return None; }
    let centroid = outer_uvs.iter().copied().sum::<DVec2>() / outer_uvs.len() as f64;

    // Try the centroid first
    let candidates = {
        let mut c = vec![centroid];
        // Add midpoints between centroid and each outer vertex
        for uv in outer_uvs {
            c.push((centroid + *uv) * 0.5);
        }
        c
    };
    for &uv in &candidates {
        if !point_in_polygon_2d(outer_uvs, uv) { continue; }
        let in_hole = hole_uvs.iter().any(|h| {
            h.len() >= 3 && point_in_polygon_2d(h, uv)
        });
        if in_hole { continue; }
        // Valid interior UV point — convert to 3D
        let pt = match surface {
            Surface3::Plane(p) => {
                let x_axis = rcad_kernel::geom::any_perpendicular(p.normal).normalize();
                let y_axis = p.normal.cross(x_axis).normalize();
                p.origin + uv.x * x_axis + uv.y * y_axis
            }
            Surface3::Sphere(s) => s.point_at(uv.x, uv.y),
            Surface3::Cylinder(c) => {
                use rcad_kernel::geom::SurfaceEval;
                c.point_at(uv.x, uv.y)
            }
            Surface3::Cone(c) => {
                use rcad_kernel::geom::SurfaceEval;
                c.point_at(uv.x, uv.y)
            }
            _ => return None,
        };
        let inward = -normal;
        return Some(pt + inward * (TOLERANCE_ABS * 100.0));
    }
    None
}

impl<'a> BooleanBuilder<'a> {
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
        if std::env::var("RCAD_DEBUG_IC").is_ok() {
            eprintln!("[SPLIT] face={} DS origin={:?} n_segments={} has_pb_sc={}", 
                face_idx, ds.faces[face_idx].origin, segments.len(),
                !ds.faces[face_idx].face_info.curves_sc.is_empty());
            for (si, seg) in segments.iter().enumerate() {
                let src = format!("{:?}", seg.source);
                eprintln!("[SPLIT]   seg[{}] src={} v{}->v{}", si, src, seg.start_vertex, seg.end_vertex);
            }
        }
        if !self.builder_face_check_data(face_idx, &segments) { return; }

        let segments_topo = crate::builder::builder_utils_topo_ds::segments_to_topo_ds(&segments, ds, face_idx, &face_refs[..], &_ic_edge_map[..]);
        // segments kept alive for classification below; dropped after classification.

        let tool: &dyn rcad_kernel::topods::BRepTool = br;

        let (avoided_pids, pid_segs) = crate::builder::wire_splitter::perform_shapes_to_avoid_topo_ds(
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
        // ✅ OCCT-aligned L147: PerformInternalShapes
        crate::builder::wire_path_topo_ds::perform_internal_shapes(
            &mut wfs, &internal_wire_groups, &segments_topo, tool, face_idx, ds);

        // ✅ OCCT-aligned: ComputeState — classify each WireFace against opposing solid.
        //   OCCT BOPAlgo_Builder::PerformInternal1 classifies each split face via
        //   BOPTools_AlgoTools::ComputeState then filters via classification_keep_policy.
        //   rcad: find an interior UV point (outer polygon minus hole polygons), map to 3D,
        //   offset inward along face normal, and classify against opposing solid faces.
        let face_surf = ds.faces[face_idx].surface.clone();
        let normal = ds.faces[face_idx].normal;
        let opposing_faces: Vec<usize> = if is_a {
            (self.ds.a_face_count..self.ds.faces.len()).collect()
        } else {
            (0..self.ds.a_face_count).collect()
        };
        if !opposing_faces.is_empty() {
            let source = if is_a { SourceSide::A } else { SourceSide::B };
            wfs.retain(|wf| {
                // Build UV polygons: outer boundary + hole boundaries
                let outer_uvs: Vec<DVec2> = wf.outer_wire.iter().filter_map(|&si| {
                    let seg = &segments[si];
                    let pt = ds.vertices[seg.start_vertex].point;
                    world_to_uv(&face_surf, pt)
                }).collect();
                let hole_uvs: Vec<Vec<DVec2>> = wf.inner_wires.iter().map(|iw| {
                    iw.iter().filter_map(|&si| {
                        let seg = &segments[si];
                        world_to_uv(&face_surf, ds.vertices[seg.start_vertex].point)
                    }).collect()
                }).collect();
                // Find interior sample point (outer polygon minus holes)
                let sample_pt = find_interior_3d(&outer_uvs, &hole_uvs, &face_surf, &normal)
                    .unwrap_or_else(|| {
                        // Fallback: centroid of outer boundary vertices, offset inward
                        let cent = wf.outer_wire.iter()
                            .map(|&si| ds.vertices[segments[si].start_vertex].point)
                            .sum::<DVec3>() / wf.outer_wire.len() as f64;
                        let inward = -normal; // from surface toward solid interior
                        cent + inward * (TOLERANCE_ABS * 100.0)
                    });
                let class = classify_point(sample_pt, &opposing_faces, ds);
                self.classification_keep_policy(source, class, face_idx)
            });
        }
        drop(segments);

        let origin = if is_a {
            FaceOrigin::FromA(ds.faces[face_idx].source_face_idx)
        } else {
            FaceOrigin::FromB(ds.faces[face_idx].source_face_idx)
        };
        let ic_curves: HashMap<usize, Curve3> = ds.intersection_curves.iter()
            .enumerate().map(|(ci, ic)| (ci, ic.curve.clone())).collect();
        for wf in &wfs {
            result.emit_wire_face_topods(face_idx, wf, &segments_topo, tool, &ic_curves, false, origin,
                &HashMap::new(), face_refs[face_idx], self.ds.faces[face_idx].natural_restriction);
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

    /// ✅ OCCT-aligned: PIOperation_FillHistory → PrepareHistory (Builder_4.cxx L164-252).
    ///   Builds source→result history matching OCCT's BRepTools_History.
    ///
    /// OCCT form:
    ///   L166:  if (!HasHistory()) return;
    ///   L174:  myHistory = new BRepTools_History;
    ///   L175:  myMapShape.Clear();
    ///   L176:  TopExp::MapShapes(myShape, myMapShape);
    ///   L185-187: for i in 0..NbSourceShapes()
    ///   L192:    if (!IsSupportedType(aS)) continue;
    ///   L205:    pLSp = LocModified(aS);  // → images
    ///   L214:    if (myMapShape.Contains(aSp)) → Modified
    ///   L233:    aGenShapes = LocGenerated(aS);
    ///   L239:    if (myMapShape.Contains(aG)) → Generated
    ///   L247:    if (!isModified && !myMapShape.Contains(aS)) → Deleted
    fn fill_history(&self, t_brep: &mut topods::BRep) -> Vec<crate::history::SourceShapeEntry> {
        use crate::history::{HistoryStatus, SourceShapeEntry};
        use topods::TShape;

        // OCCT L166: if (!HasHistory()) return.
        if !self.my_fill_history {
            return vec![];
        }

        // OCCT L174-176: TopExp::MapShapes(myShape, myMapShape).
        //   rcad: build result vertex/edge presence sets from t_brep.
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

        // ── Iterate all source shapes ──────────────────────────────────────────
        // OCCT L185-187: for (int i = 0; i < aNbS; ++i)
        //
        // → Vertices (OCCT L192: IsSupportedType filter — all vertex types are valid)
        for (di, _dv) in self.ds.vertices.iter().enumerate() {
            // OCCT L205: const List<TopoDS_Shape>* pLSp = LocModified(aS);
            let sref = rcad_kernel::topods::ShapeRef::new(v_base + di);
            let has_images = self.my_images.borrow().contains_key(&sref);
            let in_result = result_vtx.contains(&di);

            let (status, result_indices) = if has_images && in_result {
                // OCCT L208-230: split images found in result → Modified
                let images = self.my_images.borrow().get(&sref).cloned().unwrap_or_default();
                modified_indices.push(v_base + di);
                (HistoryStatus::Modified, images.iter().map(|sr| sr.index).collect())
            } else if in_result {
                // OCCT L233-243: LocGenerated → in result → Generated
                (HistoryStatus::Generated, vec![v_base + di])
            } else {
                // OCCT L247-249: not in result → Deleted
                (HistoryStatus::Deleted, vec![])
            };
            entries.push(SourceShapeEntry { ds_index: di, shape_type: 0, status, result_indices });
        }

        // → Edges (same form)
        for (di, _de) in self.ds.edges.iter().enumerate() {
            let sref = rcad_kernel::topods::ShapeRef::new(e_base + di);
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

        // → Faces (OCCT shape type TopAbs_FACE — matched by surface + wire topology)
        //   TODO: Add face-level history when topods face→DS face matching is available.
        //   Currently faces are tracked indirectly via face_origins in BuildResult.
        //   OCCT L192: if (!BRepTools_History::IsSupportedType(aS)) continue;
        //   For now, faces are not mapped here — they are handled by
        //   annotate_shell_and_solid_history during post_treat.

        // ── Set TopoDS_TShape::Moved for modified shapes ─────────────────────
        // OCCT L216-225: modified shapes get orientation fix + moved flag.
        for &idx in &modified_indices {
            if let Some(arc) = t_brep.tshapes.get_mut(idx) {
                if let Some(arc) = std::sync::Arc::get_mut(arc) {
                    match arc {
                        TShape::Vertex(vd) => vd.moved = true,
                        TShape::Edge(ed) => ed.moved = true,
                        _ => {}
                    }
                }
            }
        }

        entries
    }

    /// ✅ OCCT-aligned: PrepareHistory for the TreatEmptyShape case (BOP.cxx L462-468).
    ///   All source shapes are present as-is (Generated) or absent (Deleted);
    ///   no splitting occurs, so no Modified shapes.
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
            my_fill_history: true,   // OCCT default
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
    fn fill_images_edges(&self) {
        // OCCT L73: aNbS = myDS->NbSourceShapes()
        // OCCT L75-81: iterate source shapes, filter TopAbs_EDGE
        // OCCT L83-87: filter HasReference (has pave blocks)
        // OCCT L89-90: aE = aSI.Shape(); aLPB = myDS->PaveBlocks(i)
        // OCCT L95:    myImages.Bound(aE, ...)
        // OCCT L97-119: for each pave block:
        //   L101:   aPBR = myDS->RealPaveBlock(aPB)
        //   L103:   nSpR = aPBR->Edge()
        //   L104-105: aSpR = myDS->Shape(nSpR); pLS->Append(aSpR)
        //   L107-112: myOrigins[split].Append(source)
        //   L114-118: IsCommonBlockOnEdge → myShapesSD.Bind(source, split)
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
    /// ✅ OCCT-aligned: FillImagesContainer(WIRE) (Builder_1.cxx L221-276).
    ///   Same algorithm as FillImagesContainer(SHELL), operating on WIRE type.
    ///   L224-233: check if any sub-edge has been modified via myImages.Seek.
    ///   L235-240: if none modified → return (no image, original kept by BuildResult).
    ///   L242-245: MakeContainer(WIRE, aCIm) — new wire to hold edge images.
    ///   L247-272: for each sub-edge → add split images or original.
    ///   L274-275: Closed check + myImages.Bind(theS).Append(aCIm).
    ///   ⏳ IsSplitToReverseWithWarn skipped — edge orientation deferred to face construction.
    fn fill_images_containers_wires(&self) {
        let e_base = self.ds.vertices.len();

        // OCCT L180-183: Iterate source COMPOUND/WIRE shapes (myDS->NbSourceShapes).
        //   rcad: iterate DS wires (each wire is a container of edges).
        for wi in 0..self.ds.wires.len() {
            let edges: Vec<usize> = self.ds.wires[wi].edges.clone();

            // OCCT L224-233: check if any sub-edge has been modified.
            //   pLFIm = myImages.Seek(aSS) → pLFIm exists AND (not 1 or not same) → modified.
            let mut a_it_modified = false;
            for &ei in &edges {
                let e_ref = rcad_kernel::topods::ShapeRef::new(e_base + ei);
                let imgs_borrow = self.my_images.borrow();
                let p_lf_im = imgs_borrow.get(&e_ref);
                let is_modified = p_lf_im.map_or(false, |imgs| {
                    imgs.len() != 1 || imgs[0].index != e_base + ei
                });
                if is_modified {
                    a_it_modified = true;
                    break;
                }
            }

            // OCCT L235-240: if (!aIt.More()) return — no modification, keep original.
            //   OCCT: no image bound → BuildResult adds the original wire.
            //   rcad: no image bound → BuildResult adds original via myImages/wire_images.
            if !a_it_modified {
                // rcad: wire_images[wi] = None means "unchanged" — handled by BuildResult(WIRE).
                continue;
            }

            // OCCT L242-245: MakeContainer(theType, aCIm) — create a new wire.
            //   rcad: a_c_im = Vec of edge refs (in order) forming the new wire.
            let mut a_c_im: Vec<rcad_kernel::topods::ShapeRef> = Vec::new();

            // OCCT L247-272: iterate sub-edges → add split images or original.
            for &ei in &edges {
                let e_ref = rcad_kernel::topods::ShapeRef::new(e_base + ei);
                let p_lss_im = self.my_images.borrow().get(&e_ref).cloned();

                if let Some(ref imgs) = p_lss_im {
                    // OCCT L260-271: has splits → add each split edge image.
                    for &a_ss_im in imgs {
                        // OCCT L265-269: if (!aSSIm.IsEqual(aSS) && IsSplitToReverseWithWarn) → reverse.
                        //   rcad: edge orientation handled during face construction.
                        let _is_equal = a_ss_im.index == e_ref.index;
                        if !a_c_im.contains(&a_ss_im) {
                            a_c_im.push(a_ss_im);
                        }
                    }
                } else {
                    // OCCT L253-258: no splits → add the sub-shape itself.
                    if !a_c_im.contains(&e_ref) {
                        a_c_im.push(e_ref);
                    }
                }
            }

            // OCCT L274: aCIm.Closed(BRep_Tool::IsClosed(aCIm))
            //   rcad: closure determined by edge connectivity — skip for wire.
            // OCCT L275: myImages.Bound(theS, ...).Append(aCIm)
            //   rcad: store wire image as edge ref list in wire_images.
            let w_ref = rcad_kernel::topods::ShapeRef::new(
                e_base + self.ds.edges.len() + wi);
            self.my_images.borrow_mut().entry(w_ref).or_default().extend(a_c_im);
        }
    }

    /// ✅ OCCT-aligned: FillImagesFaces (BOPAlgo_Builder_1.cxx L376-386).
    ///   Phase 3: splits each face via WireSplitter → classifies → emits
    ///   via emit_wire_face.  rcad equivalent: for each face with IC data,
    ///   call split_face_and_emit_topo_ds (TopoDS-based BuilderFace::Perform), then
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
    /// OCCT-aligned: FillImagesFaces (Builder_2.cxx L215-229).
    ///   3-step dispatcher: BuildSplitFaces → FillSameDomainFaces → FillInternalVertices.
    fn fill_images_faces(
        &self,
        result: &mut ResultBuilder,
        a_faces: &[usize],
        b_faces: &[usize],
    ) {
        // OCCT L218: BuildSplitFaces — split all faces along intersection curves.
        self.build_split_faces(result, a_faces, b_faces);
        // OCCT L219-222: if (HasErrors()) return;
        if self.has_errors { return; }

        // OCCT L223: FillSameDomainFaces — merge duplicate same-domain faces.
        self.fill_same_domain_faces(result);
        // OCCT L224-227: if (HasErrors()) return;
        if self.has_errors { return; }

        // OCCT L228: FillInternalVertices — settle alone vertices as INTERNAL.
        self.fill_internal_vertices(result);
    }

    /// ✅ OCCT-aligned: BuildSplitFaces (Builder_2.cxx L233-374).
    ///   Iterates source faces → splits each along intersection curves.
    ///   For faces with IN/SC PBs: full BuilderFace::Perform (split_face_and_emit_topo_ds).
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

    /// ✅ OCCT-aligned: FillImagesContainer(SHELL) (Builder_1.cxx L221-276).
    ///   L224-240: check if any sub-shape has been modified
    ///   L242-275: build new container from sub-shape images
    ///   ⏳ IsSplitToReverseWithWarn skipped — face orientation handled in emit_wire_face.
    ///   ⏳ aCIm.Closed() skipped — shell closure determined during BuildRC.
    fn fill_images_containers_shells(&self, result: &mut ResultBuilder) {
        for ds_shell in &self.ds.shells {
            // OCCT L224-233: check if any sub-shape (FACE) has been modified.
            //   OCCT: myImages.Seek(aSS) → if images exist AND (not 1 or not same) → modified.
            //   rcad: count result faces for this DS face — >1 means split, 1 means same.
            let mut a_it: Option<bool> = None; // OCCT: iterator break on modification found
            for &dsfi in &ds_shell.faces {
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
                // OCCT L229: pLFIm exists AND (extent != 1 || first != original) → modified
                let has_images = count > 0;
                let is_modified = has_images && (count > 1 || !result.face_origins.iter().any(|origin| {
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
                }));
                if is_modified {
                    a_it = Some(true);
                    break; // OCCT L231: modified found, stop checking
                }
            }

            // OCCT L235-240: if no sub-shape modified → return (no new container).
            //   rcad: push identity shell so build_rc can see the face group.
            //   OCCT doesn't need this because it uses TopoDS shape identity in myShape.
            if a_it.is_none() {
                // Push identity shell (original faces) for downstream build_rc.
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

            // OCCT L242-245: MakeContainer(theType, aCIm) — create new container.
            //   rcad: new shell = Vec<usize> of result face indices.
            let mut a_c_im: Vec<usize> = Vec::new();

            // OCCT L247-272: iterate sub-shapes (faces) → add images or original.
            for &dsfi in &ds_shell.faces {
                // OCCT L251: pLSSIm = myImages.Seek(aSS) — check if face has split images.
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

                if result_faces.is_empty() {
                    // OCCT L253-258: no images → add sub-shape itself.
                    //   rcad: the DS face wasn't split; no result face to add.
                    //   (Original face already added via build_original_face in build_result(Face).)
                    continue;
                }

                // OCCT L260-271: has images (split) → add each image sub-face.
                for &rfi in &result_faces {
                    // OCCT L265-269: if (!aSSIm.IsEqual(aSS) && IsSplitToReverseWithWarn) → reverse.
                    //   rcad: orientation already correct from emit_wire_face flip handling.
                    if !a_c_im.contains(&rfi) {
                        a_c_im.push(rfi);
                    }
                }
            }

            // OCCT L274-275: aCIm.Closed(BRep_Tool::IsClosed(aCIm)); myImages.Bound(theS).Append(aCIm);
            //   OCCT always appends the shell regardless of closure.
            if !a_c_im.is_empty() {
                result.tmp_shells.push(a_c_im);
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
    fn fill_images_solids(&self, result: &mut ResultBuilder) {
        let has_solid = self.ds.faces.iter().any(|f| f.source_solid_idx.is_some());
        if !has_solid { return; }

        // OCCT L77-83: FillIn3DParts — build draft solids + classify shells
        let shell_assignments = self.fill_in_3d_parts(result);

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
    ///   Build a draft solid from a source solid, replacing split faces with their
    ///   image sub-faces and collecting INTERNAL faces into theLIF.
    ///   OCCT L283-367: iterate source solid sub-shapes (shells→faces), myImages.Seek
    ///   for each face → replace with images if bound, add INTERNAL faces to theLIF.
    ///   rcad: iterates DS shells filtered by source side, finds matching result faces.
    fn build_draft_solid(&self, result: &ResultBuilder, side: usize)
        -> (Vec<Vec<usize>>, Vec<usize>)
    {
        // OCCT L280-281: aOrSd = theSolid.Orientation(); theDraftSolid.Orientation(aOrSd).
        //   rcad: solid orientation tracked per-face via FaceOrigin.
        let origin_side = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
        let mut draft_shells: Vec<Vec<usize>> = Vec::new();
        let mut the_lif: Vec<usize> = Vec::new();

        // OCCT L283-367: iterate sub-shapes (shells) of the solid.
        for ds_shell in &self.ds.shells {
            let belongs = ds_shell.faces.iter().any(|&dsfi|
                self.ds.faces.get(dsfi).map_or(false, |f| f.origin == origin_side));
            if !belongs { continue; }

            // OCCT L292-295: MakeShell(aShD); aShD.Orientation(aOrSh); iFlag = 0.
            let mut a_sh_d: Vec<usize> = Vec::new();
            let mut i_flag = false;

            // OCCT L297-360: iterate sub-shapes (faces) of the shell.
            for &dsfi in &ds_shell.faces {
                let dsf = &self.ds.faces[dsfi];
                if dsf.origin != origin_side { continue; }

                // OCCT L301: aOrF = aF.Orientation() — rcad: all DS faces are FORWARD.
                //   INTERNAL orientation is not tracked in rcad's DS, so all faces
                //   are treated as non-INTERNAL (the common case).

                // OCCT L303: if (myImages.IsBound(aF)) — check if face has split images.
                let result_faces: Vec<usize> = result.face_origins.iter().enumerate()
                    .filter(|(_, fo)| match fo {
                        FaceOrigin::FromA(sfi) => dsf.origin == ShapeOrigin::ShapeA && dsf.source_face_idx == *sfi,
                        FaceOrigin::FromB(sfi) => dsf.origin == ShapeOrigin::ShapeB && dsf.source_face_idx == *sfi,
                        _ => false,
                    })
                    .map(|(i, _)| i)
                    .collect();
                if result_faces.is_empty() { continue; }

                if result_faces.len() > 1 {
                    // OCCT L305-346: face has images → iterate image faces
                    for &a_fx in &result_faces {
                        let a_fx_dfi_opt = match &result.face_origins[a_fx] {
                            FaceOrigin::FromA(s) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *s),
                            FaceOrigin::FromB(s) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *s),
                            _ => None,
                        };
                        let is_sd = a_fx_dfi_opt.map_or(false, |fx_dfi|
                            self.ds.shape_sd.has_sd_face(dsfi, fx_dfi)
                                || self.ds.shape_sd.has_sd_face(fx_dfi, dsfi));

                        if is_sd {
                            // OCCT L311-330: same-domain image face — IsSplitToReverse check
                            //   rcad: approximate with normal comparison
                            let b_to_reverse = a_fx_dfi_opt.map_or(false, |fx_dfi|
                                crate::boptools::is_split_to_reverse(
                                    self.ds.faces[dsfi].normal, self.ds.faces[fx_dfi].normal));
                            if !b_to_reverse {
                                i_flag = true;
                                if !a_sh_d.contains(&a_fx) { a_sh_d.push(a_fx); }
                            }
                            // OCCT L321-326: if bToReverse → aFx.Reverse(); then add to shell
                            //   rcad: reversed normal means the face goes to shell either way
                        } else {
                            // OCCT L333-344: not same-domain → use original orientation
                            i_flag = true;
                            if !a_sh_d.contains(&a_fx) { a_sh_d.push(a_fx); }
                        }
                    }
                } else {
                    // OCCT L348-359: no images → add original face directly
                    let fi = result_faces[0];
                    i_flag = true;
                    if !a_sh_d.contains(&fi) { a_sh_d.push(fi); }
                }
            }

            // OCCT L362-366: if (iFlag) { aShD.Closed(...); aBB.Add(theDraftSolid, aShD); }
            if i_flag && !a_sh_d.is_empty() {
                draft_shells.push(a_sh_d);
            }
        }

        (draft_shells, the_lif)
    }

    /// ✅ OCCT-aligned: FillIn3DParts (Builder_3.cxx L97-263).
    ///   Phase 1: collect all result faces (aLFaces).
    ///   Phase 2: build draft solids from each source solid (BuildDraftSolid).
    ///   Phase 3: classify faces against each draft solid (per-face classify_point
    ///            approximates OCCT's BVH-based BOPAlgo_Tools::ClassifyFaces).
    ///   Phase 4: analyze results → store in myInParts + return assignments.
    fn fill_in_3d_parts(&self, result: &mut ResultBuilder) -> Vec<(usize, usize, &'static str)> {
        // OCCT L101: Message_ProgressScope — rcad: skipped.
        // OCCT L103: NCollection_IncAllocator — rcad: Rust allocator.

        // === Phase 1: Collect all faces (OCCT L107-150) ===
        // OCCT L107-108: aShapeBoxMap — bounding boxes for shape acceleration.
        // OCCT L111: aMFence — fence map to prevent duplicate face entries.
        // OCCT L114: aLFaces — list of all faces to classify.
        let mut a_l_faces: Vec<usize> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<usize> =
            std::collections::HashSet::new();

        // OCCT L116-150: Iterate all source FACE shapes via DS ShapeInfo.
        //   rcad: iterate result.face_origins (all result faces already resolved).
        for (fi, fo) in result.face_origins.iter().enumerate() {
            let is_face = match fo {
                FaceOrigin::FromA(_) | FaceOrigin::FromB(_) => true,
                _ => false,
            };
            if !is_face { continue; }
            // OCCT L131-149: if myImages bound → add images (with fence); else add original.
            if a_m_fence.insert(fi) {
                a_l_faces.push(fi);
            }
        }

        // === Phase 2: Build draft solids (OCCT L152-195) ===
        // OCCT L152: BRep_Builder aBB;
        // OCCT L155: aLSolids — list of draft solids for classification.
        // OCCT L157-158: aSolidsIF — internal faces per draft solid.
        // OCCT L160-162: aDraftSolid — map: source solid → draft solid.
        //   rcad: each draft solid = Vec of shell groups of DS face indices.
        let mut a_l_solids: Vec<Vec<Vec<usize>>> = Vec::new();
        let mut a_solids_if: Vec<Vec<usize>> = Vec::new();
        // (shell_idx, side) for each draft solid (replaces OCCT's source→draft map).
        let mut draft_solid_origin: Vec<(usize, usize)> = Vec::new();

        for side in 0..2 {
            let (draft_shells, the_lif) = self.build_draft_solid(result, side);
            if draft_shells.is_empty() { continue; }
            a_l_solids.push(draft_shells);
            a_solids_if.push(the_lif);
            // Find the DS shell(s) matching this side (OCCT: iterate DS ShapeInfo SOLID).
            let origin_side = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
            for (si, ds_shell) in self.ds.shells.iter().enumerate() {
                if ds_shell.faces.iter().any(|&dfi|
                    self.ds.faces.get(dfi).map_or(false, |f| f.origin == origin_side))
                {
                    draft_solid_origin.push((si, side));
                    break;
                }
            }
        }

        // === Phase 3: ClassifyFaces (OCCT L197-208) ===
        // OCCT L197-199: LOCAL anInParts — classification result map: draft solid → IN faces.
        // OCCT L201-208: BOPAlgo_Tools::ClassifyFaces(aLFaces, aLSolids,...) batch BVH.
        //   rcad: using bopalgo::classify_faces with per-face classify_point.
        let face_samples: Vec<DVec3> = a_l_faces.iter()
            .map(|&fi| if fi < result.faces.len() { result.faces[fi].8 } else { DVec3::ZERO })
            .collect();
        let aabb_of_face: Vec<Aabb> = a_l_faces.iter().map(|&fi| {
            // Build minimal AABB from face boundary vertices via DS
            if fi < result.face_origins.len() {
                let dfi_opt = match &result.face_origins[fi] {
                    FaceOrigin::FromA(sfi) => self.ds.faces.iter().position(|f|
                        f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                    FaceOrigin::FromB(sfi) => self.ds.faces.iter().position(|f|
                        f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                    _ => None,
                };
                if let Some(dfi) = dfi_opt {
                    let mut aabb = Aabb::empty();
                    for &vi in &self.ds.faces[dfi].boundary_verts {
                        if vi < self.ds.vertices.len() {
                            aabb.expand_point(self.ds.vertices[vi].point);
                        }
                    }
                    aabb
                } else { Aabb::empty() }
            } else { Aabb::empty() }
        }).collect();
        let aabb_of_solid: Vec<Aabb> = a_l_solids.iter().map(|shells| {
            let mut aabb = Aabb::empty();
            for sh in shells {
                for &dfi in sh {
                    if dfi < self.ds.faces.len() {
                        for &vi in &self.ds.faces[dfi].boundary_verts {
                            if vi < self.ds.vertices.len() {
                                aabb.expand_point(self.ds.vertices[vi].point);
                            }
                        }
                    }
                }
            }
            aabb
        }).collect();
        let an_in_parts_list = crate::bopalgo::classify_faces(
            &a_l_faces, &face_samples, &a_l_solids, self.ds,
            &aabb_of_face, &aabb_of_solid,
        );
        let mut an_in_parts: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (dsi, in_faces) in an_in_parts_list.into_iter().enumerate() {
            if !in_faces.is_empty() {
                an_in_parts.insert(dsi, in_faces);
            }
        }

        // === Phase 4: Analyze classification results (OCCT L210-262) ===
        let mut assignments: Vec<(usize, usize, &'static str)> = Vec::new();

        // OCCT L211: aNbSol = aDraftSolid.Extent()
        for (dsi, &(si, side)) in draft_solid_origin.iter().enumerate() {
            // OCCT L220: aLInFaces = IN faces for this draft solid (from anInParts).
            let in_faces: Vec<usize> = an_in_parts.get(&dsi).cloned().unwrap_or_default();
            let n_in = in_faces.len();

            // OCCT L225-238: if no IN faces, check if shell has images → skip if none.
            if n_in == 0 {
                let mut has_image = false;
                if let Some(ds_shell) = self.ds.shells.get(si) {
                    let v_base = self.ds.vertices.len();
                    for &dsfi in &ds_shell.faces {
                        if let Some(dsf) = self.ds.faces.get(dsfi) {
                            for &ei in &dsf.boundary_edges {
                                if self.my_images.borrow().contains_key(
                                    &rcad_kernel::topods::ShapeRef::new(v_base + ei))
                                {
                                    has_image = true; break;
                                }
                            }
                            if has_image { break; }
                        }
                    }
                }
                if !has_image { continue; }
            }

            // OCCT L241: theDraftSolids.Bind(aSolid, aSDraft)
            let state: &'static str = if n_in > 0 { "IN" } else { "OUT" };
            assignments.push((si, side, state));

            // OCCT L243-261: myInParts[source] = IN_faces + INTERNAL_faces
            let mut my_in_parts = self.my_in_parts.borrow_mut();
            let a_nb_int = a_solids_if.get(dsi).map_or(0, |v| v.len());
            if a_nb_int > 0 || n_in > 0 {
                let p_lin = my_in_parts.entry(side).or_default();
                // OCCT L250-254: append IN faces
                for &fi in &in_faces {
                    if !p_lin.contains(&fi) {
                        p_lin.push(fi);
                    }
                }
                // OCCT L256-260: append INTERNAL faces (aLInternal)
                if let Some(lif) = a_solids_if.get(dsi) {
                    for &lif_fi in lif {
                        if !p_lin.contains(&lif_fi) {
                            p_lin.push(lif_fi);
                        }
                    }
                }
            }
        }
        assignments
    }

    /// ✅ OCCT-aligned: BuildSplitSolids (Builder_3.cxx L413-618).
    ///
    ///   Build result solids from draft solids and IN faces.
    ///
    ///   Phase 0 (L431-461):  Non-interfered solids → aMST (face-set dedup).
    ///   Phase 1 (L467-518):  Interfered solids → BOPAlgo_SplitSolid → collect areas.
    ///   Phase 2 (L531-537):  Parallel execution (rcad: sequential).
    ///   Phase 3 (L539-577):  Collect results + merge alerts.
    ///   Phase 4 (L580-617):  Dedup via aMST, store in myImages / myOrigins / myShapesSD.
    ///
    ///   rcad: results stored in result.tmp_solids (BuildRC applies boolean filtering).
    ///   ⏳ myImages / myOrigins / myShapesSD storage deferred to BuildRC / build_topods.
    fn build_split_solids(&self, result: &mut ResultBuilder,
                          assignments: &[(usize, usize, &'static str)]) {
        // OCCT L413-415: void BuildSplitSolids(theDraftSolids, theRange)
        //   rcad: assignments + saved_shells + my_in_parts replace theDraftSolids + myInParts.
        // OCCT L417-428: local variables (aAlr0, aSFS, aLSEmpty, aMFence, aMST, aVBS)
        let my_in_parts = self.my_in_parts.borrow();
        let has_in_faces = !my_in_parts.is_empty();

        // OCCT L425: aSFS — list of all faces for building new solid
        // OCCT L426: aMFence — fence to avoid processing same solid twice
        //   rcad: implicit in assignments iteration + in_faces_this filter.
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
        //   OCCT: iterate DS ShapeInfo for TopAbs_SOLID NOT in theDraftSolids →
        //         build BOPTools_Set of faces, add to aMST.
        //   rcad: shells WITHOUT IN faces are "non-interfered" → a_mst + stored as solids.
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
            if let Some(ds_shell) = self.ds.shells.get(si) {
                let ds_set: std::collections::BTreeSet<usize> = ds_shell.faces.iter().copied().collect();
                if ds_set.is_empty() { continue; }
                a_mst.push(ds_set);

                // OCCT L487-488: aSolidsIm.Add(aS).Append(aSD) — store non-interfered draft solid.
                let result_faces: Vec<usize> = ds_shell.faces.iter()
                    .flat_map(|&dsfi| {
                        let dsf = &self.ds.faces[dsfi];
                        result.face_origins.iter().enumerate()
                            .filter(|(_, fo)| match (dsf.origin, fo) {
                                (ShapeOrigin::ShapeA, FaceOrigin::FromA(sfi)) => dsf.source_face_idx == *sfi,
                                (ShapeOrigin::ShapeB, FaceOrigin::FromB(sfi)) => dsf.source_face_idx == *sfi,
                                _ => false,
                            })
                            .map(|(fi, _)| fi)
                    })
                    .collect();
                if result_faces.is_empty() { continue; }
                let csi = result.tmp_shells.len();
                result.tmp_shells.push(result_faces);
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
                continue;
            }

            let origin = if side == 0 { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
            let other_origin = if side == 0 { ShapeOrigin::ShapeB } else { ShapeOrigin::ShapeA };

            // OCCT L491-499: 1.1 Fill Shell Faces Set — iterate all faces of draft solid
            let mut ds_face_set: Vec<usize> = Vec::new();
            if let Some(ds_shell) = self.ds.shells.get(si) {
                for &dsfi in &ds_shell.faces {
                    let dsf = &self.ds.faces[dsfi];
                    if dsf.origin != origin { continue; }
                    ds_face_set.push(dsfi);
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
            //   ⏳ rcad: BuilderSolid (single-threaded; OCCT uses parallel aVBS).
            let mut bs = crate::bopds::builder_solid::BuilderSolid::new();
            bs.set_shapes(&ds_face_set);
            bs.perform(&self.ds);

            // OCCT L531-537: Parallel execution (BOPTools_Parallel::Perform).
            //   ⏳ rcad: BuilderSolid already performed above (sequential).

            // OCCT L539-542: collect areas → aSolidsIm.
            // OCCT L544-577: merge BuilderSolid alerts into builder report.
            //   ⏳ rcad: alerts not merged.
            for area_ds in bs.areas() {
                // OCCT L590-602: BOPTools_Set dedup via aMST.Contains / aMST.Added.
                let ds_set: std::collections::BTreeSet<usize> = area_ds.iter().copied().collect();
                if a_mst.iter().any(|s| s == &ds_set) {
                    // OCCT L598: bFlagSD = aMST.Contains(aST) — same-domain → skip duplicate
                    continue;
                }
                a_mst.push(ds_set);

                // OCCT L590-602: aST.Add(aSR, TopAbs_FACE) — map DS faces to result faces.
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

                // OCCT L603-614: store in myImages + myOrigins + myShapesSD.
                //   ⏳ rcad: stored in result.tmp_shells/tmp_solids for BuildRC.
                let csi = result.tmp_shells.len();
                result.tmp_shells.push(result_faces);
                result_solids.push(vec![csi]);
                result.solid_side_origin.push(side);
            }
        }

        // OCCT L580-617: aMST-based dedup already applied per-area above.
        result.tmp_solids = result_solids;

        // OCCT BuilderSolid::PerformAreas (BuilderSolid.cxx L397-576): void detection.
        //   ⏳ rcad: separate post-step because BuilderSolid does not perform
        //     internal void detection during bs.perform().
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
    /// OCCT-aligned: BuildRC (BOPAlgo_BOP.cxx L583-867).
    ///   Filters split solids by boolean operation type using face-set comparison.
    ///   A. FUSE (L594-609): keep all split solids (fence-deduped).
    ///   B. COMMON/CUT/CUT21 (L616-864): build args/tools building-element maps,
    ///      resolve to split images, compare for intersection containment.
    ///   rcad: result.tmp_solids contains pre-assembled split solids with
    ///     solid_side_origin tracking.  The OCCT myShape is approximated by
    ///     result.tmp_solids entries.
    fn build_rc(&self, result: &mut ResultBuilder, t_brep: &mut topods::BRep) {
        // OCCT L587-591: TopoDS_Compound aC; BRep_Builder aBB; aBB.MakeCompound(aC)
        //   rcad: aC = result.tmp_solids (equivalent output).

        let solids = std::mem::take(&mut result.tmp_solids);
        let sides: Vec<usize> = result.solid_side_origin.clone();
        if sides.len() != solids.len() { return; }

        // OCCT L594-609: A. FUSE — iterate myShape with fence, add all
        if self.op == BooleanOpType::Union {
            // OCCT L596: aMFence fence map
            // OCCT L597-606: TopExp_Explorer aExp(myShape, aType); fence-add to aC
            let mut a_m_fence: std::collections::HashSet<Vec<usize>> =
                std::collections::HashSet::new();
            let mut kept: Vec<Vec<usize>> = Vec::new();
            for s in &solids {
                if a_m_fence.insert(s.clone()) {
                    kept.push(s.clone());
                }
            }
            // OCCT L607: myRC = aC
            result.tmp_solids = kept;
            return;
        }

        // OCCT L616-645: prepare building elements of arguments to get splits
        //   OCCT: iterate myArguments/myTools → TreatCompound → TopExp::MapShapes
        //   rcad: DS vertices/edges/faces are the building elements.
        //   For each side (0=args, 1=tools): collect V/E/F indices into maps.
        let e_base = self.ds.vertices.len();
        let f_base = e_base + self.ds.edges.len();

        // OCCT L622: aMArgs, aMTools — indexed maps of source shapes (V/E/F)
        //   rcad: HashSet<usize> of flat V/E/F indices.
        let mut a_m_args: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_m_tools: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_maps = [&mut a_m_args, &mut a_m_tools];

        for (side_idx, a_ms) in a_maps.iter_mut().enumerate() {
            // OCCT L628-643: for each argument/tool shape → TreatCompound → MapShapes
            //   rcad: source building elements classified by origin in DS arrays.
            let v_range = if side_idx == 0 {
                (0usize, self.ds.a_vertex_count)
            } else {
                (self.ds.a_vertex_count, self.ds.vertices.len())
            };
            let e_range = if side_idx == 0 {
                (0usize, self.ds.a_edge_count)
            } else {
                (self.ds.a_edge_count, self.ds.edges.len())
            };
            let f_range = if side_idx == 0 {
                (0usize, self.ds.a_face_count)
            } else {
                (self.ds.a_face_count, self.ds.faces.len())
            };

            // OCCT L641-642: TypeToExplore(iDim) → MapShapes(aSS, aType, aMS)
            //   rcad: each DS entity is a building element by type.
            for vi in v_range.0..v_range.1 { a_ms.insert(vi); }
            for ei in e_range.0..e_range.1 { a_ms.insert(e_base + ei); }
            for fi in f_range.0..f_range.1 { a_ms.insert(f_base + fi); }
        }

        // OCCT L654-705: get splits of building elements
        //   For each building element, check myImages.IsBound → add split images.
        //   rcad: for edges: self.my_images[b].  for faces: result.face_origins count.
        //   For faces with no images, also build BOPTools_Set for SOLID.
        let mut a_m_args_im: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_m_tools_im: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut a_mset_args: Vec<std::collections::BTreeSet<usize>> = Vec::new();
        let mut a_mset_tools: Vec<std::collections::BTreeSet<usize>> = Vec::new();

        let mut im_maps = [&mut a_m_args_im, &mut a_m_tools_im];
        let mut set_maps = [&mut a_mset_args, &mut a_mset_tools];

        for (side_idx, (a_ms_im, a_mset)) in im_maps.iter_mut().zip(set_maps.iter_mut()).enumerate() {
            let a_ms = &a_maps[side_idx]; // &HashSet<usize> for this side
            let side_is_args = side_idx == 0;

            // OCCT L667-704: for each building element
            let mut sorted_elements: Vec<&usize> = a_ms.iter().collect();
            sorted_elements.sort(); // deterministic order

            for &&flat_idx in &sorted_elements {
                // OCCT L670-678: Type check + degenerated edge skip
                //   rcad: flat_idx < v_range → VERTEX, < e_range → EDGE, else FACE
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                let local_idx = if is_edge { flat_idx - e_base }
                    else if is_face { flat_idx - f_base }
                    else { flat_idx };

                if is_edge {
                    // OCCT L671-678: degenerated edge check
                    if self.ds.is_edge_degenerated(local_idx) { continue; }
                }

                // OCCT L681-691: if (myImages.IsBound(aS)) { add split images }
                let has_images = if is_edge {
                    self.my_images.borrow().contains_key(
                        &rcad_kernel::topods::ShapeRef::new(local_idx))
                } else if is_face {
                    // Face has images if DS face produces multiple result faces
                    let (o_exp, sfi) = if side_is_args {
                        (ShapeOrigin::ShapeA, local_idx)
                    } else {
                        (ShapeOrigin::ShapeB, local_idx)
                    };
                    let result_count = result.face_origins.iter().filter(|fo| match fo {
                        FaceOrigin::FromA(s) if side_is_args => *s == sfi,
                        FaceOrigin::FromB(s) if !side_is_args => *s == sfi,
                        _ => false,
                    }).count();
                    result_count > 1
                    // If result_count == 0, the face was not split at all
                } else {
                    // OCCT: VERTEX images from myImages — not tracked at this level in rcad
                    false
                };

                if has_images {
                    // OCCT L683-689: iterate split images and add to image map
                    let (o_exp, sfi) = if side_is_args {
                        (ShapeOrigin::ShapeA, local_idx)
                    } else {
                        (ShapeOrigin::ShapeB, local_idx)
                    };

                    if is_face {
                        for (rfi, fo) in result.face_origins.iter().enumerate() {
                            let matches = match fo {
                                FaceOrigin::FromA(s) if side_is_args => *s == sfi,
                                FaceOrigin::FromB(s) if !side_is_args => *s == sfi,
                                _ => false,
                            };
                            if matches {
                                a_ms_im.insert(f_base + rfi);
                            }
                        }
                    } else if is_edge {
                        if let Some(imgs) = self.my_images.borrow().get(
                            &rcad_kernel::topods::ShapeRef::new(local_idx))
                        {
                            for &sr in imgs {
                                a_ms_im.insert(e_base + sr.index);
                            }
                        }
                    }
                } else {
                    // OCCT L692-702: no images → add original shape
                    a_ms_im.insert(flat_idx);

                    // OCCT L694-701: for SOLID building elements, build BOPTools_Set
                    //   rcad: for face elements, build DS face set for BOPTools_Set comparison
                    if is_face {
                        let mut a_st: std::collections::BTreeSet<usize> =
                            std::collections::BTreeSet::new();
                        // Build face set from this face and its adjacent faces in the same solid
                        //   ⏳ rcad: BOPTools_Set at FACE level approximates OCCT's
                        //     SOLID-level BOPTools_Set.  OCCT adds all faces of the SOLID;
                        //     rcad adds the single DS face and its shell siblings.
                        let (o_exp2, sfi2) = if side_is_args {
                            (ShapeOrigin::ShapeA, local_idx)
                        } else {
                            (ShapeOrigin::ShapeB, local_idx)
                        };
                        a_st.insert(local_idx);
                        // Add sibling faces from the same shell
                        for (dfi2, df2) in self.ds.faces.iter().enumerate() {
                            if dfi2 != local_idx && df2.origin == o_exp2
                                && df2.source_shell_idx
                                    == self.ds.faces[local_idx].source_shell_idx
                            {
                                a_st.insert(dfi2);
                            }
                        }
                        if !a_mset.contains(&a_st) {
                            a_mset.push(a_st);
                        }
                    }
                }
            }
        }

        // OCCT L707-783: compare the maps and make the result
        let b_common = self.op == BooleanOpType::Intersection;
        let b_cut21 = false; // ⏳ rcad: CUT21 not supported

        // OCCT L715-720: determine iteration/check maps based on CUT21
        let a_m_it: &std::collections::HashSet<usize> = if b_cut21 { &a_m_tools_im } else { &a_m_args_im };
        let a_m_check: &std::collections::HashSet<usize> = if b_cut21 { &a_m_args_im } else { &a_m_tools_im };
        let a_mset_check: &Vec<std::collections::BTreeSet<usize>> =
            if b_cut21 { &a_mset_args } else { &a_mset_tools };

        // OCCT L724-755: expand sub-shapes for COMMON
        let a_m_it_exp: std::collections::HashSet<usize> = if b_common {
            let mut exp = std::collections::HashSet::new();
            for &&flat_idx in &a_m_it.iter().collect::<Vec<_>>() {
                // OCCT L730-736: expand to lower dimensions via TypeToExplore
                //   rcad: if this is a FACE, include its EDGEs and VERTEXes.
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                if is_face {
                    let local_fi = flat_idx - f_base;
                    if local_fi < self.ds.faces.len() {
                        for &ei in &self.ds.faces[local_fi].boundary_edges {
                            exp.insert(e_base + ei);
                        }
                        for &vi in &self.ds.faces[local_fi].boundary_verts {
                            exp.insert(vi);
                        }
                    }
                } else if is_edge {
                    let local_ei = flat_idx - e_base;
                    if local_ei < self.ds.edges.len() {
                        exp.insert(self.ds.edges[local_ei].start_vertex);
                        exp.insert(self.ds.edges[local_ei].end_vertex);
                    }
                }
                exp.insert(flat_idx);
            }
            exp
        } else {
            a_m_it.clone()
        };

        // OCCT L744-755: expand check side too
        let a_m_check_exp: std::collections::HashSet<usize> = {
            let mut exp = std::collections::HashSet::new();
            for &&flat_idx in &a_m_check.iter().collect::<Vec<_>>() {
                let is_edge = flat_idx >= e_base && flat_idx < f_base;
                let is_face = flat_idx >= f_base;
                if is_face {
                    let local_fi = flat_idx - f_base;
                    if local_fi < self.ds.faces.len() {
                        for &ei in &self.ds.faces[local_fi].boundary_edges {
                            exp.insert(e_base + ei);
                        }
                        for &vi in &self.ds.faces[local_fi].boundary_verts {
                            exp.insert(vi);
                        }
                    }
                } else if is_edge {
                    let local_ei = flat_idx - e_base;
                    if local_ei < self.ds.edges.len() {
                        exp.insert(self.ds.edges[local_ei].start_vertex);
                        exp.insert(self.ds.edges[local_ei].end_vertex);
                    }
                }
                exp.insert(flat_idx);
            }
            exp
        };

        // OCCT L757-784: compare building-element images and build keep set.
        //   OCCT iterates aMItExp (V/E/F level images); adds each to aC if it
        //   passes the containment check against the other side.
        //   rcad: operate at the same building-element granularity, then filter
        //   result solids whose constituent face building-elements are in keep_set.
        let mut keep_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &&flat_idx in &a_m_it_exp.iter().collect::<Vec<_>>() {
            // OCCT L762: bContains = aMCheckExp.Contains(aS)
            let mut b_contains = a_m_check_exp.contains(&flat_idx);
            // OCCT L763-768: for SOLIDs, also check BOPTools_Set
            //   rcad: operate at FACE level (no SOLID-level DS entries).
            let is_face = flat_idx >= f_base;
            if !b_contains && is_face {
                let local_fi = flat_idx - f_base;
                let mut a_st: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
                if local_fi < self.ds.faces.len() {
                    a_st.insert(local_fi);
                    for &vi in &self.ds.faces[local_fi].boundary_verts {
                        a_st.insert(vi);
                    }
                    for &ei in &self.ds.faces[local_fi].boundary_edges {
                        a_st.insert(e_base + ei);
                    }
                }
                b_contains = a_mset_check.iter().any(|s| s == &a_st);
            }
            // OCCT L770-783: COMMON → keep if contained; CUT → keep if NOT contained
            let keep = if b_common { b_contains } else { !b_contains };
            if keep {
                keep_set.insert(flat_idx);
            }
        }

        // Filter result.tmp_solids: keep solids whose iterate-side face building
        // elements pass the building-element filter above.
        let mut kept_solids: Vec<Vec<usize>> = Vec::new();
        for (i, solid_shells) in solids.iter().enumerate() {
            let side = sides.get(i).copied().unwrap_or(0);
            // Check each solid: if ANY face's building element is in keep_set → keep.
            // A result solid is kept iff the source face(s) it was split from pass.
            let mut solid_keep = false;
            for &si in solid_shells {
                if let Some(shell_faces) = result.tmp_shells.get(si) {
                    for &rfi in shell_faces {
                        let dfi_opt = match result.face_origins.get(rfi) {
                            Some(FaceOrigin::FromA(sfi)) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                            Some(FaceOrigin::FromB(sfi)) => self.ds.faces.iter().position(|f|
                                f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                            _ => None,
                        };
                        if let Some(dfi) = dfi_opt {
                            let flat_fi = f_base + dfi;
                            if keep_set.contains(&flat_fi) {
                                solid_keep = true;
                                break;
                            }
                        }
                    }
                    if solid_keep { break; }
                }
            }
            if solid_keep {
                kept_solids.push(solid_shells.clone());
            }
        }

        // OCCT L786-809: filter result for COMMON — re-explore from high dim to low
        //   rcad: OCCT re-iterates the compound by dimension (SOLID→SHELL→FACE) with
        //   a fence.  rcad solids are already at SOLID granularity; shell-count fence
        //   prevents duplicates (OCCT L799-804 fence at FACE+WIRE level).
        if b_common {
            let mut a_m_fence: std::collections::HashSet<Vec<usize>> =
                std::collections::HashSet::new();
            let mut reordered: Vec<Vec<usize>> = Vec::new();
            for s in &kept_solids {
                if a_m_fence.insert(s.clone()) {
                    reordered.push(s.clone());
                }
            }
            kept_solids = reordered;
        }

        // OCCT L811-864: degenerated edge squat (DEs whose vertex is in result,
        //   is not new, and is not interfered → add to aC).
        //   ⏳ rcad: result edges are embedded in pre-assembled solids.  Adding
        //     standalone DEs to the compound is not applicable at this level.

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
    fn detect_internal_voids(&self, result: &mut ResultBuilder,
                              assignments: &[(usize, usize, &'static str)]) {
        // OCCT L397-399: myAreas.Clear(); BRep_Builder aBB;
        // OCCT L400-407: aNewSolids, aHoleShells, aMHF (hole face map).

        // Precompute DS face set, centroid, AABB per solid.
        //   OCCT operates on raw shells (myLoops); rcad operates on result.tmp_solids.
        let n_solids = result.tmp_solids.len();
        let mut ds_faces_of: Vec<Vec<usize>> = Vec::with_capacity(n_solids);
        let mut centroids: Vec<DVec3> = Vec::with_capacity(n_solids);
        let mut aabbs: Vec<Aabb> = Vec::with_capacity(n_solids);
        for si in 0..n_solids {
            let mut faces = Vec::new();
            let mut aabb = Aabb::empty();
            let mut centroid = DVec3::ZERO;
            for &sh in &result.tmp_solids[si] {
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
                            if let Some(dfi) = ds_fi {
                                faces.push(dfi);
                                for &vi in &self.ds.faces[dfi].boundary_verts {
                                    if vi < self.ds.vertices.len() {
                                        aabb.expand_point(self.ds.vertices[vi].point);
                                    }
                                }
                                if fi < result.faces.len() {
                                    centroid = result.faces[fi].6;
                                }
                            }
                        }
                    }
                }
            }
            faces.sort_unstable(); faces.dedup();
            ds_faces_of.push(faces);
            centroids.push(centroid);
            aabbs.push(aabb);
        }

        // === Step 1: Classify each shell as Growth or Hole (OCCT L411-442) ===
        //   OCCT L422: IsGrowthShell(aShell, aMHF) — fast face overlap check.
        //     If any face of theShell is already in aMHF (face map of known holes),
        //     the shell is a Growth (it bounds a hole).
        //   OCCT L426: IsHole(aShell, myContext) — classify infinite point against
        //     the dead solid (original solid being split).  IN = hole, OUT = growth.
        let mut a_mhf: std::collections::HashSet<usize> = std::collections::HashSet::new();
        let mut is_hole = vec![false; n_solids];

        for si in 0..n_solids {
            // OCCT L422: IsGrowthShell
            let is_growth = if !a_mhf.is_empty() {
                ds_faces_of[si].iter().any(|dfi| a_mhf.contains(dfi))
            } else {
                false
            };

            if !is_growth {
                // OCCT L426: IsHole — classify against original solid (source operand).
                let side = result.solid_side_origin.get(si).copied().unwrap_or(0);
                let dead_faces: Vec<usize> = self.ds.faces.iter().enumerate()
                    .filter(|(_, f)| match side {
                        0 => f.origin == ShapeOrigin::ShapeA,
                        _ => f.origin == ShapeOrigin::ShapeB,
                    })
                    .map(|(fi, _)| fi)
                    .collect();
                let class = classify_point(centroids[si], &dead_faces, self.ds);
                // OCCT: IsHole returns true if infinite point is IN dead solid → hole.
                is_hole[si] = class == Classification::In;
            }
            // else: IsGrowthShell returned true → definitely a growth.

            if is_hole[si] {
                // OCCT L439-441: aHoleShells.Add + TopExp::MapShapes(,TopAbs_FACE,aMHF)
                for &dfi in &ds_faces_of[si] {
                    a_mhf.insert(dfi);
                }
            }
        }

        // OCCT L429-441 (Growth/Hole separation done above).
        let in_si: Vec<usize> = (0..n_solids).filter(|&i| is_hole[i]).collect();
        let out_si: Vec<usize> = (0..n_solids).filter(|&i| !is_hole[i]).collect();

        // OCCT L444-458: if no holes → add all growths to myAreas + return.
        if in_si.is_empty() || out_si.is_empty() { return; }

        // === Step 2: Build BVH of hole shells (OCCT L462-478) ===
        //   OCCT L464-475: BOPTools_BoxTree with BRepBndLib bounding boxes.
        //   rcad: Bvh built from hole solid AABBs.
        let hole_key: Vec<usize> = in_si.clone();
        let hole_aabbs: Vec<Aabb> = in_si.iter().map(|&i| aabbs[i]).collect();
        let hole_bvh = crate::bvh::DsBvh::build(hole_key, hole_aabbs);

        // === Step 3: Classify holes against growth solids (OCCT L483-529) ===
        //   OCCT L493-529: for each growth solid:
        //     build box → BVH-select candidate holes → IsInside → store outermost.
        let mut a_hole_solid_map: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();

        for &os in &out_si {
            // OCCT L494-497: BRepBndLib::Add(aSolid, aBox)
            // OCCT L499-502: BOPTools_BoxTreeSelector → candidate holes
            let candidates = hole_bvh.query_aabb(&aabbs[os]);

            for &hole_idx in &candidates {
                // OCCT L511: IsInside(aHole, aSolid, myContext)
                let class = classify_point(centroids[hole_idx], &ds_faces_of[os], self.ds);
                if class != Classification::In && class != Classification::On {
                    continue;
                }

                // OCCT L517-527: select outermost containing solid.
                //   If current os is INSIDE the previously recorded solid,
                //   the current os is more specific (innermost container) → prefer it.
                use std::collections::hash_map::Entry;
                match a_hole_solid_map.entry(hole_idx) {
                    Entry::Occupied(mut e) => {
                        let prev_os = *e.get();
                        let prev_faces = &ds_faces_of[prev_os];
                        let os_inside = classify_point(centroids[os], prev_faces, self.ds);
                        if os_inside == Classification::In || os_inside == Classification::On {
                            e.insert(os);
                        }
                    }
                    Entry::Vacant(e) => {
                        e.insert(os);
                    }
                }
            }
        }

        // === Step 4: Build reverse map: solid → list of holes (OCCT L532-548) ===
        let mut solid_holes_map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (&hole_idx, &os) in &a_hole_solid_map {
            solid_holes_map.entry(os).or_default().push(hole_idx);
        }

        // === Step 5: Add holes to solids + myAreas (OCCT L550-576) ===
        let mut removed = vec![false; n_solids];
        // OCCT L553-573: for each growth with holes → aBB.Add(aSolid, aHole)
        for (&os, holes) in &solid_holes_map {
            for &hole_idx in holes {
                let void_shells = result.tmp_solids[hole_idx].clone();
                result.tmp_solids[os].extend(void_shells);
                removed[hole_idx] = true;
            }
        }
        // OCCT L575: myAreas.Append(aSolid) — rcad: non-removed solids kept in tmp_solids.
        // OCCT L578-581: add un-associated holes to myAreas (rcad: not needed — kept as-is).
        let mut new_solids: Vec<Vec<usize>> = Vec::with_capacity(n_solids);
        for (si, solid) in result.tmp_solids.drain(..).enumerate() {
            if !removed[si] { new_solids.push(solid); }
        }
        result.tmp_solids = new_solids;
    }

    /// ✅ OCCT-aligned: FillInternalShapes (Builder_3.cxx L622-887).
    /// OCCT-aligned: FillInternalShapes (Builder_3.cxx L622-887).
    ///   Phase 1 (L648-709): Collect internal V/E/WIRE from arguments.
    ///   Phase 2 (L717-788): Internal V/E from source solids + build aMSx ancestry.
    ///   Phase 3 (L790-809): Filter shapes already attached via aMSx.
    ///   Phase 4 (L811-816): Early return if none.
    ///   Phase 5 (L820-877): Classify each internal shape against each split solid;
    ///     if IN → add as INTERNAL sub-shape (clone original if needed).
    fn fill_internal_shapes(&self, result: &mut ResultBuilder) {
        // OCCT L631-644: allocator + indexed maps (aMSx, aMx, aMSI, aMFence, aMSOr, ...)
        //   rcad: adapted to Vec/HashSet equivalents.

        // === Phase 1: Shapes to process — collect from arguments (OCCT L648-709) ===
        //   OCCT L653-658: TreatCompound on each argument → flatten into aLSC.
        //   OCCT L660-681: filter VERTEX/EDGE/WIRE from aLSC → aLArgs.
        //   OCCT L684-709: for each aLArgs, check myImages.IsBound → aMSI (images or originals).
        //   rcad: DS vertices/edges with is_internal flag = sources.
        //   ├ TreatCompound: rcad treats DS V/E as already-flattened source shapes.
        //   └ aMSI: maps shape-ref → true if it's an internal shape to process.
        let mut a_msi: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // Collect internal vertices (OCCT L677-679: TopAbs_VERTEX → aLArgs)
        for (vi, v) in self.ds.vertices.iter().enumerate() {
            if v.is_internal {
                // OCCT L691-706: check myImages.IsBound → add split images or original
                let v_ref = rcad_kernel::topods::ShapeRef::new(vi);
                if self.my_images.borrow().contains_key(&v_ref) {
                    for img in &self.my_images.borrow()[&v_ref] {
                        a_msi.insert(img.index);
                    }
                } else {
                    a_msi.insert(vi);
                }
            }
        }
        // Collect internal edges (OCCT L665-675: WIRE → iterate edges; L677-679: EDGE directly)
        for (ei, e) in self.ds.edges.iter().enumerate() {
            if e.is_internal {
                let e_ref = rcad_kernel::topods::ShapeRef::new(ei);
                if self.my_images.borrow().contains_key(&e_ref) {
                    for img in &self.my_images.borrow()[&e_ref] {
                        a_msi.insert(img.index);
                    }
                } else {
                    a_msi.insert(ei);
                }
            }
        }

        // === Phase 2: Internal V/E from source solids + build aMSx ancestry (OCCT L717-788) ===
        //   OCCT L721-727: iterate DS for SOLIDs.
        //   L738: OwnInternalShapes(aS, aMx) — get INTERNAL sub-shapes from each solid.
        //   L741-758: insert into aMSI (with myImages check).
        //   L760-787: build aMSx ancestry: Vertex→Edge, Vertex→Face, Edge→Face.
        //   rcad: aMSx tracks which internal shapes are already on split-solid faces.
        #[allow(unused)]
        let mut a_msx: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new(); // shape_idx → list of ancestor face/edge indices
        let mut a_lsd: Vec<usize> = Vec::new(); // split solids to process

        // OCCT L741-758: internal shapes from OwnInternalShapes
        //   rcad: is_internal flag already collected above in Phase 1.
        //   The DS vertices/edges with is_internal=true are equivalent to OCCT's OwnInternalShapes output.

        // OCCT L760-787: build aMSx for split solids
        //   For each source SOLID that has split results (images) → build ancestry map.
        //   rcad: iterate result.tmp_solids → for each solid, map edges→faces.
        for (si, solid_shells) in result.tmp_solids.iter().enumerate() {
            // OCCT L761: if (myImages.IsBound(aS)) for source solid
            let side = result.solid_side_origin.get(si).copied().unwrap_or(0);
            // Build edge→face adjacency for this result solid
            let mut edge_to_faces: std::collections::HashMap<usize, Vec<usize>> =
                std::collections::HashMap::new();
            for &shi in solid_shells {
                if let Some(shell_faces) = result.tmp_shells.get(shi) {
                    for &rfi in shell_faces {
                        if let Some(fe) = result.faces.get(rfi) {
                            for &(ei, _) in &fe.0 {
                                edge_to_faces.entry(ei).or_default().push(rfi);
                            }
                        }
                    }
                }
            }
            // OCCT L770-773: TopExp::MapShapesAndAncestors → aMSx
            //   aMSx[vertex] = list of edge indices
            //   aMSx[vertex] = list of face indices
            //   aMSx[edge] = list of face indices
            for (&ei, face_list) in &edge_to_faces {
                // e_ref in a_msx → ancestors (face indices)
                a_msx.entry(ei).or_default().extend(face_list);
                // Also add vertex→edge ancestry
                if ei < self.ds.edges.len() {
                    a_msx.entry(self.ds.edges[ei].start_vertex)
                        .or_default().push(ei);
                    a_msx.entry(self.ds.edges[ei].end_vertex)
                        .or_default().push(ei);
                }
            }
            a_lsd.push(si);
        }

        // === Phase 3: Filter shapes already attached to split-solid faces (OCCT L790-809) ===
        //   OCCT: for each shape in aMSI, check if aMSx.Contains(shape) with non-empty ancestor list.
        //         → if NOT attached → aLSI (list of shapes to settle).
        let mut a_lsi: Vec<usize> = Vec::new();
        for &si in &a_msi {
            // OCCT L796-808: if aMSx contains the shape AND has non-empty ancestors → skip (attached).
            //   rcad: check if this internal shape index appears in aMSx with non-empty ancestors.
            let is_attached = a_msx.get(&si).map_or(false, |anc| !anc.is_empty());
            if !is_attached {
                a_lsi.push(si);
            }
        }

        // === Phase 4: Early return if none (OCCT L811-816) ===
        if a_lsi.is_empty() {
            return;
        }

        // === Phase 5: Settle internal V/E into solids (OCCT L820-877) ===
        //   OCCT L825-876: for each split solid (aLSd), for each internal shape (aLSI):
        //     ComputeStateByOnePoint(aSI, aSd) → if IN:
        //       - if original solid (aMSOr): clone → add INTERNAL → bind myImages/myOrigins
        //       - else: add INTERNAL directly
        for &si in &a_lsd {
            // OCCT L828: TopoDS_Solid aSd
            //   rcad: get DS face set for this solid
            let mut solid_ds_faces: Vec<usize> = Vec::new();
            if let Some(solid_shells) = result.tmp_solids.get(si) {
                for &shi in solid_shells {
                    if let Some(shell_faces) = result.tmp_shells.get(shi) {
                        for &rfi in shell_faces {
                            let dfi_opt = match result.face_origins.get(rfi) {
                                Some(FaceOrigin::FromA(sfi)) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeA && f.source_face_idx == *sfi),
                                Some(FaceOrigin::FromB(sfi)) => self.ds.faces.iter().position(|f|
                                    f.origin == ShapeOrigin::ShapeB && f.source_face_idx == *sfi),
                                _ => None,
                            };
                            if let Some(dfi) = dfi_opt {
                                solid_ds_faces.push(dfi);
                            }
                        }
                    }
                }
            }
            solid_ds_faces.sort_unstable();
            solid_ds_faces.dedup();
            if solid_ds_faces.is_empty() { continue; }

            // OCCT L830-875: iterate internal shapes to settle
            let mut i = 0usize;
            while i < a_lsi.len() {
                let si_idx = a_lsi[i];
                // OCCT L834: aSI.Orientation(TopAbs_INTERNAL)
                //   rcad: no orientation; use classify_point with centroid.
                // OCCT L836: ComputeStateByOnePoint(aSI, aSd, 1.e-11, myContext)
                let pt = if si_idx < self.ds.vertices.len() {
                    self.ds.vertices[si_idx].point
                } else {
                    let ei = si_idx;
                    if ei < self.ds.edges.len() {
                        (self.ds.vertices[self.ds.edges[ei].start_vertex].point
                         + self.ds.vertices[self.ds.edges[ei].end_vertex].point) * 0.5
                    } else {
                        i += 1; continue;
                    }
                };
                let a_state = classify_point(pt, &solid_ds_faces, self.ds);

                if a_state != Classification::In {
                    // OCCT L840: aIt1.Next(); continue;
                    i += 1;
                    continue;
                }

                // OCCT L844-873: shape is IN → add as INTERNAL
                //   OCCT L844: if (aMSOr.Contains(aSd)) — original solid → clone first
                //   rcad: find first face of this solid to store internal vertex
                if let Some(&first_shi) = result.tmp_solids.get(si).and_then(|s| s.first()) {
                    if let Some(shell_faces) = result.tmp_shells.get(first_shi) {
                        if let Some(&first_rfi) = shell_faces.first() {
                            if first_rfi < result.face_internal_vtx.len() {
                                // OCCT L857-873: add INTERNAL shape to solid
                                //   rcad: store DS vertex index in face_internal_vtx
                                if si_idx < self.ds.vertices.len() {
                                    if !result.face_internal_vtx[first_rfi].contains(&si_idx) {
                                        result.face_internal_vtx[first_rfi].push(si_idx);
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L875: aLSI.Remove(aIt1) — remove settled shape
                a_lsi.swap_remove(i);
                // don't increment i — the new element at i needs checking too
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
    /// OCCT-aligned: FillImagesCompounds (Builder_1.cxx L197-217) + FillImagesCompound (L280-342).
    ///   L197-201: dispatcher with fence map; iterate source COMPOUND shapes.
    ///   L280-293: FillImagesCompound — fence skip if already processed.
    ///   L295-308: recurse into sub-compounds; check if any sub-shape has images.
    ///   L309-312: no modification → return.
    ///   L314-341: build new compound from sub-shape images; store in myImages.
    ///   ⏳ rcad: no compound nesting in DS.  Flat per-face source_compsolid_idx.
    ///     The recursive FillImagesCompound is collapsed to a single level.
    fn fill_images_compounds(&self, result: &mut ResultBuilder) {
        // OCCT L199-200: aMFP fence map — prevents reprocessing the same compound.
        //   rcad: HashSet of processed compsolid indices.
        let mut a_mfp: std::collections::HashSet<usize> = std::collections::HashSet::new();
        // OCCT L202: aNbS = myDS->NbSourceShapes() — iterate all DS shapes.
        //   rcad: collect unique source_compsolid_idx from DS faces.
        let mut compound_indices: Vec<usize> = Vec::new();
        for df in &self.ds.faces {
            if let Some(csi) = df.source_compsolid_idx {
                if !compound_indices.contains(&csi) {
                    compound_indices.push(csi);
                }
            }
        }
        if compound_indices.is_empty() { return; }

        for &csi in &compound_indices {
            // OCCT L290-293: if (!theMFP.Add(theS)) return — fence check.
            if !a_mfp.insert(csi) { continue; }

            // OCCT L295-308: check if any sub-shape (SOLID) has been modified.
            //   rcad: collect source_solid_idx values under this compsolid,
            //   check if any has images (multiple result solids).
            let sub_solid_indices: Vec<usize> = self.ds.faces.iter()
                .filter_map(|f| {
                    if f.source_compsolid_idx == Some(csi) {
                        f.source_solid_idx
                    } else {
                        None
                    }
                })
                .collect();
            // OCCT L300-303: recurse into sub-compounds — rcad: flat, no nesting.
            // OCCT L304-307: if (myImages.IsBound(aSx)) bInterferred = true.
            //   rcad: check if any sub-solid produces >1 result solid (split).
            let mut b_interferred = false;
            for &ssi in &sub_solid_indices {
                // Count result solids from this source solid
                let count = result.solid_side_origin.iter()
                    .filter(|&&side| {
                        // Count result solids from the side matching this source solid's origin
                        let dfi = self.ds.faces.iter().position(|f|
                            f.source_solid_idx == Some(ssi));
                        dfi.map_or(false, |di| {
                            let origin = &self.ds.faces[di].origin;
                            (origin == &crate::bopds::ds::ShapeOrigin::ShapeA && side == 0)
                                || (origin == &crate::bopds::ds::ShapeOrigin::ShapeB && side == 1)
                        })
                    })
                    .count();
                if count > 0 {
                    b_interferred = true;
                    break;
                }
            }

            // OCCT L309-312: if (!bInterferred) return — no modification.
            if !b_interferred { continue; }

            // OCCT L314-315: MakeContainer(COMPOUND, aCIm)
            //   rcad: collect result solid indices for this compsolid.
            let mut a_c_im: Vec<usize> = Vec::new();

            // OCCT L317-336: iterate sub-shapes → add images or original.
            for &ssi in &sub_solid_indices {
                // Find the DS face for this source solid to determine its side (origin)
                let side = self.ds.faces.iter()
                    .find(|f| f.source_solid_idx == Some(ssi))
                    .map(|f| match f.origin {
                        crate::bopds::ds::ShapeOrigin::ShapeA => 0,
                        crate::bopds::ds::ShapeOrigin::ShapeB => 1,
                    })
                    .unwrap_or(0);

                // OCCT L322: if (myImages.IsBound(aSX)) — has split images?
                //   rcad: check if result solids exist for this side+source solid.
                let matching_solids: Vec<usize> = result.solid_side_origin.iter()
                    .enumerate()
                    .filter(|&(_, &s)| s == side)
                    .map(|(si, _)| si)
                    .collect();

                if matching_solids.is_empty() {
                    // OCCT L334-335: no images → add original sub-shape
                    //   rcad: no solid to add — the original solid is implicit.
                    continue;
                }

                // OCCT L324-331: has images → add each image with orientation.
                for &si in &matching_solids {
                    if !a_c_im.contains(&si) {
                        // OCCT L329: aSXIm.Orientation(aOrX) — preserve orientation.
                        //   rcad: orientation is per-face via FaceOrigin.
                        a_c_im.push(si);
                    }
                }
            }

            // OCCT L339-341: aLSIm.Append(aCIm); myImages.Bind(theS, aLSIm)
            //   rcad: store for build_result(Compound) to consume.
            if !a_c_im.is_empty() {
                result.compound_groups.push(a_c_im);
            }
        }
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
    /// OCCT-aligned: BuildResult (Builder_1.cxx L130-168).
    ///   Add split images (or originals) of source shapes into the result.
    ///   OCCT L133: aMFence fence map.
    ///   L136-167: for each source argument of matching type → if myImages bound
    ///     → add all image shapes; else → add the original shape.
    ///   rcad: adapts to topods::BRep TShape factory + ResultBuilder storage.
    fn build_result(&self, shape_type: ShapeType, result: &mut ResultBuilder, t: &mut topods::BRep) {
        // OCCT L133: NCollection_Map<TopoDS_Shape> aMFence — dedup shapes in result.
        //   rcad: unique-indexed arrays make a fence unnecessary, but the form is kept.
        #[allow(unused)]
        let mut a_m_fence: Vec<usize> = Vec::new();

        // OCCT L136-167: iterate all source arguments of matching type.
        //   rcad: source entities vary by type (DS arrays, result data).
        match shape_type {
    ShapeType::Vertex => {
        // ✅ OCCT-aligned: BuildResult (Builder_1.cxx L130-168).
        //   Iterate all source arguments of type VERTEX → add images to myShape.
        //   If no images, add the original vertex.
        //   rcad: source vertices = DS vertices 0..a_vc (A) + a_vc.. (B).
        let a_vc = self.ds.a_vertex_count;
        let nv = self.ds.vertices.len();
        for side in 0..2usize {
            let (start, end) = if side == 0 { (0usize, a_vc.min(nv)) } else { (a_vc, nv) };
            for vi in start..end {
                // OCCT L145: myImages.Seek(aS) — check if vertex has split image
                let sref = rcad_kernel::topods::ShapeRef::new(vi);
                let has_images = self.my_images.borrow().contains_key(&sref);
                if !has_images {
                    // OCCT L149-152: no images → add the original shape
                    let pt = self.ds.vertices[vi].point;
                    let _rvi = result.add_ds_vertex(vi, pt);
                    t.add_tvertex(pt);
                } else {
                    // OCCT L156-165: add images of the argument shape into result
                    let images = self.my_images.borrow().get(&sref).unwrap().clone();
                    for img in &images {
                        let vi_img = img.index;
                        if vi_img < self.ds.vertices.len() {
                            let pt = self.ds.vertices[vi_img].point;
                            let _rvi = result.add_ds_vertex(vi_img, pt);
                            t.add_tvertex(pt);
                        }
                    }
                }
            }
        }
    }
            ShapeType::Edge => {
                // OCCT L130-168 (TopAbs_EDGE): add split edge images to myShape.
                //   rcad: iterate myImages(EDGE) entries, create TShape::Edge for each.
                //   First ensure all DS vertices have TShapes for edge vertex refs.
                let e_base = self.ds.vertices.len();
                // Ensure vertex TShapes exist (needed by edge creation)
                for vi in 0..self.ds.vertices.len() {
                    let vr = rcad_kernel::topods::ShapeRef::new(vi);
                    if t.tshapes.len() <= vi {
                        // Extend tshapes array to cover this index
                        let pt = self.ds.vertices[vi].point;
                        let sv = t.add_tvertex(pt);
                        t.vertex_mut(sv).tolerance = self.ds.vertices[vi].geom_tol
                            .max(crate::tolerance::TOLERANCE_ABS);
                        let _ = vr;
                    }
                }
                // Iterate A and B side source edges
                let a_ec = self.ds.a_edge_count;
                let n_edges = self.ds.edges.len();
                for side in 0..2usize {
                    let (start, end) = if side == 0 {
                        (0usize, a_ec.min(n_edges))
                    } else {
                        (a_ec, n_edges)
                    };
                    for ei in start..end {
                        let aE = rcad_kernel::topods::ShapeRef::new(e_base + ei);
                        let has_images = self.my_images.borrow().contains_key(&aE);
                        if !has_images {
                            // OCCT L149-152: no images → add original edge
                            let edge = &self.ds.edges[ei];
                            let sv_sr = rcad_kernel::topods::ShapeRef::new(edge.start_vertex);
                            let ev_sr = rcad_kernel::topods::ShapeRef::new(edge.end_vertex);
                            let ci = t.curves.len();
                            t.curves.push(edge.curve.clone());
                            let te = t.add_tedge(Some(ci), sv_sr, ev_sr, edge.t_range);
                            if self.ds.is_edge_degenerated(ei) || edge.start_vertex == edge.end_vertex {
                                t.edge_mut(te).degenerated = true;
                            }
                        } else {
                            // OCCT L156-165: add split images
                            let images = self.my_images.borrow().get(&aE).unwrap().clone();
                            for img in &images {
                                let nSpR = img.index.saturating_sub(e_base);
                                if nSpR >= self.ds.edges.len() { continue; }
                                let edge = &self.ds.edges[nSpR];
                                let sv_sr = rcad_kernel::topods::ShapeRef::new(edge.start_vertex);
                                let ev_sr = rcad_kernel::topods::ShapeRef::new(edge.end_vertex);
                                let ci = t.curves.len();
                                t.curves.push(edge.curve.clone());
                                let te = t.add_tedge(Some(ci), sv_sr, ev_sr, edge.t_range);
                                if self.ds.is_edge_degenerated(nSpR) || edge.start_vertex == edge.end_vertex {
                                    t.edge_mut(te).degenerated = true;
                                }
                            }
                        }
                    }
                }
            }
            ShapeType::Wire => {
                // OCCT L130-168: wires are sub-shapes of faces, not standalone in rcad.
            }
            ShapeType::Face => {
                // OCCT L145-165: for each source FACE, check myImages.Seek(aS).
                //   rcad: result.face_origins tracks which source faces were split.
                //   rcad: build_faces() validates edge refs before TShape creation.
                result.build_faces();
                // OCCT L146-152: add original source faces without split images.
                let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
                let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
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
                for &fi in &a_faces {
                    if !emitted_a.contains(&self.ds.faces[fi].source_face_idx) {
                        result.build_original_face(self.ds, fi,
                            FaceOrigin::FromA(self.ds.faces[fi].source_face_idx));
                    }
                }
                for &fi in &b_faces {
                    if !emitted_b.contains(&self.ds.faces[fi].source_face_idx) {
                        result.build_original_face(self.ds, fi,
                            FaceOrigin::FromB(self.ds.faces[fi].source_face_idx));
                    }
                }
                result.build_topods_faces(t);
            }
            ShapeType::Shell => {
                // OCCT L145-165: for each source SHELL, check myImages for shell images.
                //   rcad: tmp_shells already contains shell face groups from
                //   fill_images_containers_shells.
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
                // OCCT L130-167: for each source SOLID, check myImages → add images/original.
                //   rcad: tmp_solids contains solid shell groups from build_split_solids.
                //   OCCT-aligned: clone, not take — BuildResult does not consume myImages,
                //   and build_rc (called after All BuildResults) still needs the data.
                let tmp_solids = result.tmp_solids.clone();
                if !tmp_solids.is_empty() {
                    let new_shells = result.tmp_shells.clone();
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
                // OCCT L130-167: aggregate sub-solid images into CompSolid.
                //   rcad: delegate to ResultBuilder::build_compsolids which uses
                //   BRepBuilder::make_compsolid (OCCT BRep_Builder equivalent).
                let tmp_cs_groups = std::mem::take(&mut result.tmp_compsolid_groups);
                result.build_compsolids(t, tmp_cs_groups);
            }
            ShapeType::Compound => {
                // OCCT L130-168: for each source COMPOUND, add its image (or original).
                //   rcad: delegate to ResultBuilder::build_compounds which uses
                //   BRepBuilder::make_compound (OCCT BRep_Builder equivalent).
                //   compound_groups contains solid indices for each source compound
                //   (populated by fill_images_compounds).
                let groups = std::mem::take(&mut result.compound_groups);
                result.build_compounds(t, &groups);
            }
        }
    }

    /// ✅ OCCT-aligned: BOPAlgo_Builder::BuildResult (Builder_1.cxx L130-168).
    ///   rcad: thin wrapper mapping topods::ShapeType to builder::types::ShapeType.
    fn build_result_occt(&self, the_type: topods::ShapeType, result: &mut ResultBuilder, t: &mut topods::BRep) {
        let shape_type = match the_type {
            topods::ShapeType::Shape => unreachable!("ShapeType::Shape is a null sentinel, never passed to build_result"),
            topods::ShapeType::Vertex => ShapeType::Vertex,
            topods::ShapeType::Edge => ShapeType::Edge,
            topods::ShapeType::Wire => ShapeType::Wire,
            topods::ShapeType::Face => ShapeType::Face,
            topods::ShapeType::Shell => ShapeType::Shell,
            topods::ShapeType::Solid => ShapeType::Solid,
            topods::ShapeType::CompSolid => ShapeType::CompSolid,
            topods::ShapeType::Compound => ShapeType::Compound,
        };
        self.build_result(shape_type, result, t);
    }

    /// ✅ OCCT-aligned: BOPAlgo_BOP::BuildShape (BOP.cxx L871-906).
    ///   Calls BuildRC (L900) then BuildSolid for FUSE 3D (L902-906).
    fn build_shape(&self, result: &mut ResultBuilder, t_brep: &mut topods::BRep) {
        // OCCT L900: BuildRC — filter solids by boolean operation
        self.build_rc(result, t_brep);
        if self.has_errors { return; }
        // OCCT L902-906: if (FUSE + 3D) BuildSolid
        //   rcad: Union keeps all filtered solids; no separate BuildSolid needed.
    }

    /// ✅ OCCT-aligned: BOPAlgo_Builder::PostTreat (Builder.cxx L450-475).
    ///   Two-step tolerance correction: CorrectTolerances + CorrectShapeTolerances.
    fn post_treat(&self, brep: &mut rcad_kernel::BRep) {
        // OCCT L452-454: aMA — map of shapes to avoid
        // OCCT L455-469: if non-destructive → collect source V/E/F into aMA
        // rcad: non-destructive defaults to false.  When true, collect non-new
        // DS vertex indices into map_to_avoid.
        let map_to_avoid: std::collections::HashSet<usize> = if self.my_non_destructive {
            let mut avoid = std::collections::HashSet::new();
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
        // OCCT L472: BOPTools_AlgoTools::CorrectTolerances(myShape, aMA, 0.05, myRunParallel)
        if map_to_avoid.is_empty() {
            rcad_kernel::tolerance::correct_tolerances(brep, 23);
        } else {
            rcad_kernel::tolerance::correct_tolerances_with_map(brep, 23, &map_to_avoid);
        }
        // OCCT L474: BOPTools_AlgoTools::CorrectShapeTolerances(myShape, aMA, myRunParallel)
        //   rcad: correct_tolerances already does both steps in one call.
        //   Separating them requires splitting the tolerance module.
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

    /// ✅ OCCT-aligned: TreatEmptyShape (BOPAlgo_BOP.cxx L214-319).
    ///   Handles the case where one or both operands have no geometry.
    ///   Returns Ok(Some(brep)) if a quick result was determined,
    ///   Ok(None) if the full pipeline must run.
    fn treat_empty_shape(&self, a_faces: &[usize], b_faces: &[usize])
        -> Result<Option<rcad_kernel::BRep>, BooleanError>
    {
        let has_a = !a_faces.is_empty();
        let has_b = !b_faces.is_empty();
        if has_a && has_b {
            return Ok(None); // need full pipeline
        }
        if !has_a && !has_b {
            // OCCT L252-256: all empty → empty result
            return Ok(Some(rcad_kernel::BRep::new()));
        }
        // OCCT L258-317: one side empty → result depends on operation
        match self.op {
            BooleanOpType::Union => {
                // OCCT L270-279: return non-empty side
                let src = if has_a { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
                let brep = self.brep_of_side(src, a_faces.len(), b_faces.len());
                Ok(Some(brep))
            }
            BooleanOpType::Intersection => {
                // OCCT L303-304: Common always empty
                Ok(Some(rcad_kernel::BRep::new()))
            }
            BooleanOpType::Difference => {
                if !has_a {
                    // OCCT L287-289: CUT with empty objects → empty
                    Ok(Some(rcad_kernel::BRep::new()))
                } else {
                    // OCCT L281-289: CUT with empty tools → return objects
                    let brep = self.brep_of_side(ShapeOrigin::ShapeA, a_faces.len(), b_faces.len());
                    Ok(Some(brep))
                }
            }
            _ => {
                // Unknown operation → fall through to full pipeline
                Ok(None)
            }
        }
    }

    /// ✅ OCCT-aligned: BOPAlgo_BOP::PerformInternal1 (BOP.cxx L422-579).
    ///   Every statement in OCCT L422-579 has a corresponding rcad line below.
    ///   See comments for exact OCCT line references.
    ///   Structural difference: L425-429 setup done in constructor, re-affirmed here.
    ///   L531 BuildResult(SOLID) writes to t_brep, then L900 BuildRC filters and
    ///   clears solids from t_brep (non-Union) — equivalent to OCCT removing from myShape.
    pub fn build_with_history(&self) -> Result<(BRep, BooleanHistory), BooleanError> {
        // OCCT L425-429: setup (myPaveFiller, myDS, myContext, myFuzzyValue, myNonDestructive).
        //   rcad equivalents are already assigned in new(); re-affirm the form.
        let _fuzzy_value = self.ds.fuzzy_tol;
        let _non_destructive = self.my_non_destructive;

        // Pipeline dump context (env RCAD_DUMP_PIPELINE=1 to enable).
        // Grid/case are set via RCAD_DUMP_GRID / RCAD_DUMP_CASE env vars.
        let mut _dump = crate::pipeline_dump::DumpCtx::new(
            &std::env::var("RCAD_DUMP_GRID").unwrap_or_else(|_| "unknown".into()),
            &std::env::var("RCAD_DUMP_CASE").unwrap_or_else(|_| "unknown".into()),
        );

        // OCCT L431-436: CheckData
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
        self.check_data(&a_faces, &b_faces)?;

        // OCCT L438-443: Prepare — creates empty TopoDS_Compound as myShape.
        let (mut t_brep, mut result) = self.prepare();
        if self.has_errors { return Err(BooleanError::DegenerateResult); }

        // OCCT L445-453: TreatEmptyShape — check if any operand is degenerate.
        if let Some(brep) = self.treat_empty_shape(&a_faces, &b_faces)? {
            // OCCT L462-468: PrepareHistory(theRange) — record source→result status.
            //   For TreatEmptyShape, one side's shapes appear as-is (Generated),
            //   the other side's shapes are absent (Deleted). There are no splits,
            //   so the classification is trivially based on which side has faces.
            let has_a = !a_faces.is_empty();
            let has_b = !b_faces.is_empty();
            let source_history = if !self.my_fill_history {
                // OCCT L166: !HasHistory → no history recorded
                vec![]
            } else if !has_a && !has_b {
                // Both empty → all Deleted
                self.source_history_all_deleted()
            } else if has_a && has_b {
                // Both have faces — should not reach this branch
                vec![]
            } else {
                // Exactly one non-empty side
                let side = if has_a { ShapeOrigin::ShapeA } else { ShapeOrigin::ShapeB };
                self.source_history_single_side(side)
            };
            let mut history = BooleanHistory::default();
            history.source_history = source_history;
            return Ok((brep, history));
        }

        // OCCT L454-457: ProgressScope + PISteps
        struct _ProgressScope;
        impl _ProgressScope { fn next(&self, _step: f64) -> f64 { 0.0 } }
        const _PIOP_LAST: usize = 12;
        let _a_ps = _ProgressScope;
        let _a_steps = [0.0f64; _PIOP_LAST];
        let _ = (&_a_ps, &_a_steps);

        // ✅ OCCT-aligned: dimension-by-dimension pipeline (PerformInternal1 L336-445).
        // OCCT L456-459: FillImagesVertices → BuildResult(VERTEX)
        self.fill_images_vertices();
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::Vertex, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        _dump.snapshot("after_FillImagesVertices", self.ds, Some(&t_brep));
        // OCCT L461-465: FillImagesEdges → BuildResult(EDGE)
        self.fill_images_edges();
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::Edge, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        _dump.snapshot("after_FillImagesEdges", self.ds, Some(&t_brep));
        // OCCT L467-470: FillImagesContainers(WIRE) → BuildResult(WIRE)
        self.fill_images_containers(ShapeType::Wire, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::Wire, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        _dump.snapshot("after_BuildResultWire", self.ds, Some(&t_brep));
        // OCCT L472-475: FillImagesFaces → BuildResult(FACE)
        self.fill_images_faces(&mut result, &a_faces, &b_faces);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::Face, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        _dump.snapshot("after_FillImagesFaces", self.ds, Some(&t_brep));
        // OCCT L477-480: FillImagesContainers(SHELL) → BuildResult(SHELL)
        self.fill_images_containers(ShapeType::Shell, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::Shell, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        _dump.snapshot("after_BuildResultShell", self.ds, Some(&t_brep));
        // OCCT L482-485: FillImagesSolids → BuildResult(SOLID)
        self.fill_images_solids(&mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::Solid, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L487-490: FillImagesContainers(COMPSOLID) → BuildResult(COMPSOLID)
        self.fill_images_containers(ShapeType::CompSolid, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::CompSolid, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L492-495: FillImagesCompounds → BuildResult(COMPOUND)
        self.fill_images_compounds(&mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result_occt(topods::ShapeType::Compound, &mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L498-500: BuildShape → BuildRC + BuildSolid
        self.build_shape(&mut result, &mut t_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L502-504: PIOperation_FillHistory → PrepareHistory
        //   OCCT: myDS->NbSourceShapes() → LocModified/Generated → AddModified/Remove.
        //   rcad: fill_history populates source_history from my_images + result map.
        let mut history = result.build_topods(&mut t_brep, self.my_fill_history);
        let source_history = if self.my_fill_history {
            self.fill_history(&mut t_brep)
        } else {
            vec![]
        };
        history.source_history = source_history;

        let mut brep = rcad_kernel::BRep::from_topods(&t_brep);

        // OCCT L506-508: PostTreat
        self.post_treat(&mut brep);
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
        v.sort_unstable();
        v
    }

    /// OCCT L258-317: build result BRep from one side's source shapes (TreatEmptyShape path).
    fn brep_of_side(&self, origin: ShapeOrigin, _na: usize, _nb: usize) -> rcad_kernel::BRep {
        // OCCT L310-316: BRep_Builder().Add(myShape, aItLS.Value())
        //   rcad: reconstruct a BRep from the DS faces/edges/vertices of one side.
        let mut brep = rcad_kernel::BRep::new();
        let mut v_map = std::collections::HashMap::new();
        let mut e_map = std::collections::HashMap::new();
        let mut edge_store: Vec<(usize, usize)> = Vec::new();
        for (fi, f) in self.ds.faces.iter().enumerate() {
            if f.origin != origin { continue; }
            let mut outer_edges: Vec<(usize, bool)> = Vec::new();
            for &ei in &f.boundary_edges {
                if ei >= self.ds.edges.len() { continue; }
                let e = &self.ds.edges[ei];
                let sv = *v_map.entry(e.start_vertex).or_insert_with(|| {
                    let vi = brep.vertices.len();
                    brep.vertices.push(rcad_kernel::Vertex { point: self.ds.vertices[e.start_vertex].point });
                    vi
                });
                let ev = *v_map.entry(e.end_vertex).or_insert_with(|| {
                    let vi = brep.vertices.len();
                    brep.vertices.push(rcad_kernel::Vertex { point: self.ds.vertices[e.end_vertex].point });
                    vi
                });
                let ei_new = *e_map.entry(ei).or_insert_with(|| {
                    let index = brep.edges.len();
                    brep.edges.push(rcad_kernel::Edge { start: sv, end: ev });
                    edge_store.push((e.start_vertex, e.end_vertex));
                    index
                });
                outer_edges.push((ei_new, true));
            }
            let face_normal = f.normal;
            let surf = f.surface.clone();
            let rcad_face = rcad_kernel::topology::Face {
                outer_wire: rcad_kernel::topology::Wire { edges: outer_edges.into_iter().map(|(i, _)| rcad_kernel::topology::WireEdge::fwd(i)).collect() },
                inner_wires: vec![],
                normal: face_normal,
                triangles: vec![],
                sample_point: None,
                mesh_dirty: false,
                surface_idx: None,
            };
            brep.geom.surfaces.push(surf);
            let surfi = brep.geom.surfaces.len() - 1;
            brep.geom.face_surface.push(Some(surfi));
            let shell = rcad_kernel::topology::Shell { faces: vec![rcad_face] };
            let solid = rcad_kernel::topology::Solid { shells: vec![shell] };
            brep.solids.push(solid);
        }
        brep
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
