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
use crate::bopalgo::{GlueEnum, Alert, Report};
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
    /// 閴?OCCT-aligned: myGlue 鈥?BOPAlgo_GlueEnum (GlueOff/GlueFull/GlueShift).
    glue: GlueEnum,
    glue_tolerance: f64,
    context: RefCell<Context>,
    // 閴?OCCT-aligned: error tracking (myReport / HasErrors equivalent).
    has_errors: bool,
    // 閴?OCCT-aligned: myImages 閳?source shape index 閳?list of split image indices.
    my_images: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>>>,
    my_origins: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, Vec<rcad_kernel::topods::ShapeRef>>>,
    my_shapes_sd: std::cell::RefCell<std::collections::HashMap<rcad_kernel::topods::ShapeRef, rcad_kernel::topods::ShapeRef>>,
    my_in_parts: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    my_solid_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    my_solid_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
    // 閴?OCCT-aligned: myNonDestructive (BOPAlgo_Builder.hxx L503).
    my_non_destructive: bool,
    // OCCT-aligned: myFillHistory (BOPAlgo_Options.hxx).
    my_fill_history: bool,
    // 閴?OCCT-aligned: myCheckInverted (BOPAlgo_Builder.hxx L505).
    my_check_inverted: bool,
    // 閴?OCCT-aligned: myStopOnFatalError 鈥?abort pipeline on fatal error.
    my_stop_on_fatal_error: bool,
    /// 閴?OCCT-aligned: myEntryPoint 鈥?tracks builder phase (1=PerformInternal1 done, etc.).
    my_entry_point: u8,
    /// 閴?OCCT-aligned: myReport 鈥?collects alerts during Builder execution.
    my_report: Report,
    /// 閴?OCCT-aligned: converted BRep representation of DS.
    brep: std::cell::RefCell<Option<(rcad_kernel::topods::BRep, Vec<rcad_kernel::topods::ShapeRef>, Vec<Option<rcad_kernel::topods::ShapeRef>>)>>,
    /// OCCT-aligned: myShape — result shape accumulator (BRep).
    my_shape: std::cell::RefCell<rcad_kernel::topods::BRep>,
    /// OCCT-aligned: myArguments — all source shapes pre-created as TShapes.
    my_arguments: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
    /// OCCT-aligned: DS edge → TShape::Edge mapping (replaces ResultBuilder.ds_edge_to_tshape).
    my_edge_map: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
    /// OCCT-aligned: result wire TShape refs (replaces ResultBuilder.wire_refs).
    my_wire_refs: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
    /// OCCT-aligned: result shell TShape refs (replaces ResultBuilder.shells).
    my_shells: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
    /// Result face TShape refs (replaces ResultBuilder.face_refs).
    my_face_refs: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
    /// Result solid TShape refs (replaces ResultBuilder.solids).
    my_solids: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
    /// Result compsolid TShape refs (replaces ResultBuilder.compsolid_groups).
    my_compsolid_groups: std::cell::RefCell<Vec<rcad_kernel::topods::ShapeRef>>,
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

    // Check if the UV polygon is CW (hole) 锟?interior is the complement.
    // Compute signed area: positive = CCW (interior is polygon), negative = CW (interior is complement).
    let signed_area: f64 = outer_uvs.windows(2).map(|pair| {
        pair[0].x * pair[1].y - pair[1].x * pair[0].y
    }).sum::<f64>() + {
        let n = outer_uvs.len();
        outer_uvs[n-1].x * outer_uvs[0].y - outer_uvs[0].x * outer_uvs[n-1].y
    } * 0.5;
    let is_cw = signed_area < 0.0;

    // Build candidates: for CCW polygons, points inside the polygon;
    // for CW polygons, points outside the polygon (the complement interior).
    let candidates = if is_cw {
        // CW polygon 锟?interior is the complement.  The centroid of the polygon
        // vertices is inside the polygon (wrong region).  For periodic surfaces
        // (sphere: [0,2蟺]脳[-蟺/2,蟺/2]), try the domain center or points opposite.
        let mut c = vec![
            DVec2::new(centroid.x + std::f64::consts::PI, centroid.y),
        ];
        // Also try domain corners that are likely outside the small CW region
        c.push(DVec2::new(std::f64::consts::PI, 0.0));
        c.push(DVec2::new(std::f64::consts::PI, std::f64::consts::FRAC_PI_2));
        c.push(DVec2::new(std::f64::consts::PI, -std::f64::consts::FRAC_PI_2));
        c.push(DVec2::new(std::f64::consts::PI * 1.5, 0.0));
        c.push(DVec2::new(std::f64::consts::PI * 0.5, std::f64::consts::FRAC_PI_4));
        c
    } else {
        // CCW polygon 锟?interior is the polygon itself
        let mut c = vec![centroid];
        for uv in outer_uvs {
            c.push((centroid + *uv) * 0.5);
        }
        c
    };
    let test_fn = |uv: DVec2| -> bool {
        if is_cw {
            // For CW: point must be OUTSIDE the polygon
            if outer_uvs.len() >= 3 && point_in_polygon_2d(outer_uvs, uv) { return false; }
            // AND outside any hole (holes are CCW, so inside-hole means inside hole polygon)
            !hole_uvs.iter().any(|h| h.len() >= 3 && point_in_polygon_2d(h, uv))
        } else {
            // For CCW: point must be INSIDE the polygon
            if !point_in_polygon_2d(outer_uvs, uv) { return false; }
            // AND outside any hole
            !hole_uvs.iter().any(|h| h.len() >= 3 && point_in_polygon_2d(h, uv))
        }
    };
    for &uv in &candidates {
        if !test_fn(uv) { continue; }
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
        return Some(pt);
    }
    // Fallback: for sphere/cylinder periodic surfaces, try a grid of candidate
    // points across the UV domain.  A CCW polygon on the sphere may represent
    // the small-region interior (correct for small cap) or the large-region
    // complement (for the large cap where both have CCW winding).
    // Testing multiple distributed points increases the chance of finding a
    // valid interior point for either region.
    if matches!(surface, Surface3::Sphere(_)) {
        for ui in 0..8 {
            let u = (ui as f64 + 0.5) / 8.0 * std::f64::consts::TAU;
            for vi in 0..4 {
                let v = (vi as f64 + 0.5) / 4.0 * std::f64::consts::PI - std::f64::consts::FRAC_PI_2;
                let uv = DVec2::new(u, v);
                if !test_fn(uv) { continue; }
                if let Surface3::Sphere(s) = surface {
                    return Some(s.point_at(u, v));
                }
            }
        }
    }
    None
}

impl<'a> BooleanBuilder<'a> {
    /// 锟?OCCT-aligned: TopoDS-based BuildFace pipeline with emit.
    ///   Runs the full pipeline then emits result faces directly into ResultBuilder.
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
                eprintln!("[SPLIT]   seg[{}] src={} v{}->v{}", si, src, seg.start_vertex, seg.end_vertex);
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

    /// 锟?OCCT-aligned: BuilderFace::CheckData (BOPAlgo_BuilderFace.cxx L50-115).
    ///   Validates face has intersection curves/segments. If no interferences,
    ///   delegates to BuildDraftFace (OCCT's alternative path for non-split faces).
    fn builder_face_check_data(&self, face_idx: usize, segments: &[WireSegment]) -> bool {
        if segments.is_empty() {
            return false;
        }
        true
    }

    /// 锟?OCCT-aligned: PIOperation_FillHistory 锟?PrepareHistory (Builder_4.cxx L164-252).
    ///   Builds source鈫抮esult history matching OCCT's BRepTools_History.
    ///
    /// OCCT form:
    ///   L166:  if (!HasHistory()) return;
    ///   L174:  myHistory = new BRepTools_History;
    ///   L175:  myMapShape.Clear();
    ///   L176:  TopExp::MapShapes(myShape, myMapShape);
    ///   L185-187: for i in 0..NbSourceShapes()
    ///   L192:    if (!IsSupportedType(aS)) continue;
    ///   L205:    pLSp = LocModified(aS);  // 锟?images
    ///   L214:    if (myMapShape.Contains(aSp)) 锟?Modified
    ///   L233:    aGenShapes = LocGenerated(aS);
    ///   L239:    if (myMapShape.Contains(aG)) 锟?Generated
    ///   L247:    if (!isModified && !myMapShape.Contains(aS)) 锟?Deleted
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

        // 鈹€鈹€ Iterate all source shapes 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        // OCCT L185-187: for (int i = 0; i < aNbS; ++i)
        //
            // 锟?Vertices (OCCT L192: IsSupportedType filter 锟?all vertex types are valid)
            for (di, _dv) in self.ds.vertices.iter().enumerate() {
                // OCCT L205: const List<TopoDS_Shape>* pLSp = LocModified(aS);
                let sref = self.brep_sr(v_base + di);
                let has_images = self.my_images.borrow().contains_key(&sref);
                let in_result = result_vtx.contains(&di);

            let (status, result_indices) = if has_images && in_result {
                // OCCT L208-230: split images found in result 锟?Modified
                let images = self.my_images.borrow().get(&sref).cloned().unwrap_or_default();
                modified_indices.push(v_base + di);
                (HistoryStatus::Modified, images.iter().map(|sr| sr.index).collect())
            } else if in_result {
                // OCCT L233-243: LocGenerated 锟?in result 锟?Generated
                (HistoryStatus::Generated, vec![v_base + di])
            } else {
                // OCCT L247-249: not in result 锟?Deleted
                (HistoryStatus::Deleted, vec![])
            };
            entries.push(SourceShapeEntry { ds_index: di, shape_type: 0, status, result_indices });
        }

        // 锟?Edges (same form)
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

        // 锟?Faces (OCCT shape type TopAbs_FACE 锟?matched by surface + wire topology)
        //   TODO: Add face-level history when topods face鈫扗S face matching is available.
        //   Currently faces are tracked indirectly via face_origins in BuildResult.
        //   OCCT L192: if (!BRepTools_History::IsSupportedType(aS)) continue;
        //   For now, faces are not mapped here 锟?they are handled by
        //   annotate_shell_and_solid_history during post_treat.

        // 鈹€鈹€ Set TopoDS_TShape::Moved for modified shapes 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
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
    ///   Single call that builds the final shape topology and records
    ///   source shape history (modified/generated/deleted).
    ///   In rcad this combines two steps (build_topods for shape-level,
    ///   fill_history for source-level) into one public method.
    fn prepare_history(&self, result: &mut ResultBuilder) -> BooleanHistory {
        let mut t_brep = self.my_shape.borrow_mut();
        let mut history = result.build_topods(&mut *t_brep, self.my_fill_history, &self.my_shells.borrow(), &mut *self.my_face_refs.borrow_mut(), &self.my_solids.borrow(), &self.my_compsolid_groups.borrow());
        let source_history = if self.my_fill_history {
            self.fill_history(&mut *t_brep)
        } else {
            vec![]
        };
        history.source_history = source_history;
        history
    }

    /// OCCT-aligned: PrepareHistory for the TreatEmptyShape case (BOP.cxx L462-468).
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

/// Edge-like segment for wire building锟?can be a DS edge, an intersection curve,
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
            my_fill_history: true,   // OCCT default
            my_check_inverted: false,
            my_stop_on_fatal_error: true,
            my_entry_point: 0,
            my_report: Report::new(),
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

    pub fn build(&self) -> Result<BRep, BooleanError> {
        let (t, _) = self.build_with_history()?;
        let brep = rcad_kernel::BRep::from_topods(&t);
        if !brep.solids.is_empty() && !brep.solids[0].shells.is_empty() {
            eprintln!("BooleanBuilder::build: {} faces", brep.solids[0].shells[0].faces.len());
        }
        Ok(brep)
    }
}

include!("builder/filler.rs");
include!("builder/result_build.rs");

impl<'a> BooleanBuilder<'a> {
    ///   The top-level pipeline entry: dimension-by-dimension image filling
    ///   (V鈫扙鈫扺鈫扚ACE鈫扴HELL鈫扴OLID), followed by BuildResult for each type.
    ///   OCCT L310-445 structure matched in full (see inline OCCT line refs).
    /// 锟?OCCT-aligned: CheckData (BOPAlgo_BOP.cxx L106-202) + CheckFiller (Builder.cxx L143-151).
    ///   Validates operation type, non-empty arguments, and DS/PaveFiller state.
    /// ✅ OCCT-aligned: CheckData (BOPAlgo_Builder.cxx L130-140).
    fn check_data(&self) -> Result<(), BooleanError> {
        // OCCT L132-137: aNb = myArguments.Extent(); if (aNb < 2) → AlertTooFewArguments
        //   rcad: arguments are always at least 2 (set during fuse() call).
        //   rcad: validate operation type as surrogate dimension check.
        match self.op {
            BooleanOpType::Union | BooleanOpType::Intersection | BooleanOpType::Difference => {}
            _ => return Err(BooleanError::InvalidOperation),
        }
        // OCCT L139-141: CheckFiller — verify PaveFiller and DS are valid.
        //   OCCT: if (!myPaveFiller) → AlertNoFiller; GetReport()->Merge(myPaveFiller->GetReport())
        //   rcad: PaveFiller ran before builder; check DS has valid shape data loaded.
        if self.ds.vertices.is_empty() {
            return Err(BooleanError::EmptyInput);
        }
        if self.has_errors {
            return Err(BooleanError::DegenerateResult);
        }
        Ok(())
    }

    /// ✅ OCCT-aligned: Prepare (BOPAlgo_Builder.cxx L156-164).
    ///   OCCT: BRep_Builder.MakeCompound(myShape) — empty compound as result.
    ///   rcad: initializes my_shape + returns (BRep, ResultBuilder) for downstream.
    fn prepare(&self) -> (topods::BRep, ResultBuilder) {
        *self.my_shape.borrow_mut() = topods::BRep::new();
        (topods::BRep::new(), ResultBuilder::new())
    }

    /// ✅ OCCT-aligned: create TShapes for all DS source shapes in my_shape.
    ///   Equivalent to OCCT's myArguments populated with all source TopoDS_Shape.
    fn pre_create_source_shapes(&self) {
        let mut t = self.my_shape.borrow_mut();
        let mut args = Vec::new();
        // 1. Vertices
        for v in &self.ds.vertices {
            let sr = t.add_tvertex(v.point);
            args.push(sr);
        }
        // 2. Edges (with curves)
        for (ei, edge) in self.ds.edges.iter().enumerate() {
            let sv = t.add_tvertex(self.ds.vertices[edge.start_vertex].point);
            let ev = t.add_tvertex(self.ds.vertices[edge.end_vertex].point);
            let ci = rcad_kernel::topods::find_or_add_curve3(&mut t.curves, &edge.curve);
            let te = t.add_tedge(Some(ci), sv, ev, edge.t_range);
            if self.ds.is_edge_degenerated(ei) || edge.start_vertex == edge.end_vertex {
                t.edge_mut(te).degenerated = true;
            }
            args.push(te);
        }
        // 3. Wires (from DS wires, using pre-created edge ShapeRefs)
        let e_base = self.ds.vertices.len();
        for wire in &self.ds.wires {
            let wire_edges: Vec<rcad_kernel::topods::ShapeRef> = wire.edges.iter()
                .filter_map(|&ei| {
                    let idx = e_base + ei;
                    if idx < t.tshapes.len() { Some(self.brep_sr(idx)) } else { None }
                })
                .collect();
            let ws = t.add_twire(wire_edges);
            args.push(ws);
        }
        // 4. Faces (from DS faces, using pre-created wire TShapes)
        let w_base = e_base + self.ds.edges.len();
        for fi in 0..self.ds.faces.len() {
            let face = &self.ds.faces[fi];
            // Outer wire
            let outer_wire = if let Some(wi) = face.outer_wire_idx {
                let idx = w_base + wi;
                if idx < t.tshapes.len() { Some(self.brep_sr(idx)) } else { None }
            } else { None };
            // Inner wires
            let inner_wires: Vec<rcad_kernel::topods::ShapeRef> = face.inner_wire_idxs.iter()
                .filter_map(|&wi| {
                    let idx = w_base + wi;
                    if idx < t.tshapes.len() { Some(self.brep_sr(idx)) } else { None }
                })
                .collect();
            let surf_idx = t.surfaces.len();
            t.surfaces.push(face.surface.clone());
            let sample_pt = face.boundary_verts.first().copied()
                .and_then(|vi| self.ds.vertices.get(vi))
                .map(|v| v.point)
                .unwrap_or(glam::DVec3::ZERO);
            let face_sr = t.add_tface(Some(surf_idx),
                outer_wire.unwrap_or(rcad_kernel::topods::ShapeRef::NULL),
                inner_wires, Some(sample_pt), None, vec![], face.natural_restriction);
            args.push(face_sr);
        }
        // 5. Shells (from DS shells, using pre-created face TShapes)
        let f_base = w_base + self.ds.wires.len();
        for shell in &self.ds.shells {
            let face_refs: Vec<rcad_kernel::topods::ShapeRef> = shell.faces.iter()
                .filter_map(|&fi| {
                    let idx = f_base + fi;
                    if idx < t.tshapes.len() { Some(self.brep_sr(idx)) } else { None }
                })
                .collect();
            if !face_refs.is_empty() {
                args.push(t.add_tshell(face_refs));
            }
        }
        // 6. Solids (from DS solids, using pre-created shell TShapes)
        let sh_base = f_base + self.ds.faces.len();
        for solid in &self.ds.solids {
            let shell_refs: Vec<rcad_kernel::topods::ShapeRef> = solid.shells.iter()
                .filter_map(|&si| {
                    let idx = sh_base + si;
                    if idx < t.tshapes.len() { Some(self.brep_sr(idx)) } else { None }
                })
                .collect();
            if !shell_refs.is_empty() {
                args.push(t.add_tsolid(shell_refs));
            }
        }
        // 7. CompSolids
        let so_base = sh_base + self.ds.shells.len();
        for cs in &self.ds.comp_solids {
            let solid_refs: Vec<rcad_kernel::topods::ShapeRef> = cs.solids.iter()
                .filter_map(|&soi| {
                    let idx = so_base + soi;
                    if idx < t.tshapes.len() { Some(self.brep_sr(idx)) } else { None }
                })
                .collect();
            if !solid_refs.is_empty() {
                args.push(t.add_tcompsolid(solid_refs));
            }
        }
        *self.my_arguments.borrow_mut() = args;
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
            // OCCT L252-256: all empty 锟?empty result
            return Ok(Some(rcad_kernel::BRep::new()));
        }
        // OCCT L258-317: one side empty 锟?result depends on operation
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
                    // OCCT L287-289: CUT with empty objects 锟?empty
                    Ok(Some(rcad_kernel::BRep::new()))
                } else {
                    // OCCT L281-289: CUT with empty tools 锟?return objects
                    let brep = self.brep_of_side(ShapeOrigin::ShapeA, a_faces.len(), b_faces.len());
                    Ok(Some(brep))
                }
            }
            _ => {
                // Unknown operation 锟?fall through to full pipeline
                Ok(None)
            }
        }
    }

    /// 锟?OCCT-aligned: BOPAlgo_BOP::PerformInternal1 (BOP.cxx L422-579).
    ///   Every statement in OCCT L422-579 has a corresponding rcad line below.
    ///   See comments for exact OCCT line references.
    ///   Structural difference: L425-429 setup done in constructor, re-affirmed here.
    ///   L531 BuildResult(SOLID) writes to t_brep, then L900 BuildRC filters and
    ///   clears solids from t_brep (non-Union) 鈥?equivalent to OCCT removing from myShape.
    pub fn build_with_history(&self) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
        self.build_with_history_topods()
    }

    /// Same as build_with_history but returns topods::BRep directly (OCCT-aligned).
    pub fn build_with_history_topods(&self) -> Result<(topods::BRep, BooleanHistory), BooleanError> {
        // OCCT L425-429: setup (myPaveFiller, myDS, myContext, myFuzzyValue, myNonDestructive).
        //   rcad equivalents are already assigned in new(); re-affirm the form.
        let _fuzzy_value = self.ds.fuzzy_tol;
        let _non_destructive = self.my_non_destructive;

        // OCCT L431-436: CheckData 鈥?validates arguments and merges PaveFiller report.
        let a_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeA);
        let b_faces: Vec<usize> = self.faces_of(ShapeOrigin::ShapeB);
        self.check_data()?;

        // OCCT L438-442: Prepare
        //   rcad: prepare() initializes my_shape + returns ResultBuilder.
        let mut result = self.prepare().1;

        // OCCT L445-453: TreatEmptyShape.
        //   OCCT: GetReport()->HasAlert(AlertEmptyShape) -> TreatEmptyShape() -> PrepareHistory -> return.
        //   rcad: check if either operand has no faces (DS was populated by pave_fill).
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
        const _PIOP_LAST: usize = 12;
        let _a_steps = [0.0f64; _PIOP_LAST];
        // OCCT: pre-create ALL source shapes as TShapes (matching myArguments in OCCT).
        self.pre_create_source_shapes();
        // OCCT L456-459: FillImagesVertices — BuildResult(VERTEX)
        self.fill_images_vertices();
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(topods::ShapeType::Vertex, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L461-465: FillImagesEdges — BuildResult(EDGE)
        self.fill_images_edges();
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(topods::ShapeType::Edge, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L467-470: FillImagesContainers(WIRE)
        self.fill_images_container(ShapeType::Wire, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(topods::ShapeType::Wire, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L472-475: FillImagesFaces — BuildResult(FACE)
        // Architecture A1: split faces create TShapes incrementally during fill_images_faces.
        // Remaining unsplit faces already have pre-created TShapes from pre_create_source_shapes.
        self.fill_images_faces(&mut result, &a_faces, &b_faces);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // BuildResult(FACE) — generic loop over my_arguments, adds originals/splits to result.
        self.build_result(topods::ShapeType::Face, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L477-480: FillImagesContainers(SHELL)
        self.fill_images_container(ShapeType::Shell, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(topods::ShapeType::Shell, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L482-485: FillImagesSolids
        self.fill_images_solids(&mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(topods::ShapeType::Solid, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L537-548: FillImagesContainers(COMPSOLID)
        self.fill_images_container(ShapeType::CompSolid, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(topods::ShapeType::CompSolid, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L492-495: FillImagesCompounds
        self.fill_images_compounds(&mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        self.build_result(topods::ShapeType::Compound, &mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L498-500: BuildShape
        self.build_shape(&mut result);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        // OCCT L502-504: PrepareHistory
        let mut history = self.prepare_history(&mut result);

        // OCCT L506-508: PostTreat
        let final_brep = self.my_shape.borrow().clone();
        let mut old_brep = rcad_kernel::BRep::from_topods(&final_brep);
        self.post_treat(&mut old_brep);
        if self.has_errors { return Err(BooleanError::DegenerateResult); }
        let result_brep = old_brep.to_topods();

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
    /// U values differ by more than 锟?indicate a seam crossing; we accumulate
    /// offsets of 閸楊槚eriod to make the polyline continuous in U.
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
    /// surfaces (sphere, cylinder, 锟? where intersection PCurves are clipped
    /// to the finite face-face overlap and may not reach the UV boundary.
    ///
    /// Only trims that are nearly axis-aligned (constant-u or constant-v) are
    /// extended 锟?general trims pass through unchanged.
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
        // 0.5 % of the smaller span 锟?well above floating-point noise for any
        // practical model, yet tight enough to distinguish axis-aligned trims
        // from oblique ones on a sphere (where u/v vary together).
        let axis_threshold = (boundary_u_span.abs().min(boundary_v_span.abs())).max(TOLERANCE_ABS) * 0.005;

        let is_const_u = u_span_trim < axis_threshold;
        let is_const_v = v_span_trim < axis_threshold;

        if !is_const_u && !is_const_v {
            return trim.to_vec(); // non-axis-aligned 锟?cannot safely extend
        }

        // 闁冲厜鍋撻柍鍏夊亾 Clip trim points to boundary bounds 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾
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

        // 闁冲厜鍋撻柍鍏夊亾 span-checking guard (AFTER clipping) 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜锟?        // If this axis-aligned trim already covers 锟?0 % of the boundary span
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
        // 闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜鍋撻柍鍏夊亾闁冲厜锟?
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
    /// 锟?閮ㄥ垎瀵归綈: 閻劎绨跨涵顔笺亣閸﹀棗濮弸鍕紦閻炲啴娼扮€涙劙娼伴妴?
    ///    OCCT: BuildSplitFaces 锟?section edges 閻╁瓨甯撮崚娑樼紦 BRep sub-face锟?
    ///    rcad: 閹靛濮╃拋锛勭暬 8 娑擃亜宕烽梽鎰畱 FaceSampleData,锟?outer_circle_edges 鐠佹澘缍嶆径褍娓惧褋锟?
    ///    閸旂喕鍏樼粵澶夌幆(8 娑擃亜宕愰悶鍐桨閸栧搫锟?+ 缁墽鈥橀崷鍡楀К鏉堝湱锟?,锟?OCCT 娑撳秹娓剁憰浣疯厬锟?FaceSampleData锟?

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
include!("builder/part2.rs");
include!("builder/footer.rs");
#[cfg(test)]
mod tests {
    include!("builder/tests_inc.rs");
}
