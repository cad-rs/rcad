use std::collections::{HashMap, HashSet, BTreeSet};
use glam::DVec2; use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::topods;
use std::cell::RefCell;
use crate::bopds::ds::*;
use crate::classify::{Classification, classify_point};
use crate::history::{BooleanHistory, FaceOrigin, ShellOrigin, SolidOrigin, VertexOrigin, EdgeOrigin, HistoryTracker};
use crate::inttools::context::Context;
use crate::tolerance::*;
use crate::builder::types::{BooleanOpType, WireSegment, WireEdgeSource, WireOrientation};
use crate::builder::SourceSide;

use crate::builder::angle_2d::angle_2d;
use crate::builder::wire_splitter::{world_to_uv, edge_uv_tangent, edge_angle_2d, are_verts_coincident, is_edge_isoline};
use crate::builder::edge_builders::{build_sphere_seam_segments, build_cylinder_seam_segments, is_split_to_reverse};

///  ?OCCT-aligned: compare two Curve3 for identity (same TShape).
pub(crate) fn curve_eq(a: &Curve3, b: &Curve3) -> bool {
 match (a, b) {
 (Curve3::Circle(ca), Curve3::Circle(cb)) => {
 (ca.center - cb.center).length_squared() < TOLERANCE_ABS_SQ
 && (ca.normal - cb.normal).length_squared() < TOLERANCE_ABS_SQ
 && (ca.radius - cb.radius).abs() < TOLERANCE_ABS
 }
 (Curve3::Line(la), Curve3::Line(lb)) => {
 (la.origin - lb.origin).length_squared() < TOLERANCE_ABS_SQ
 && (la.direction - lb.direction).length_squared() < TOLERANCE_ABS_SQ
 }
 _ => false,
 }
}

pub(crate) fn hash_point(p: DVec3) -> u64 {
 // Quantize to tolerance grid for spatial hashing
 let scale = 1.0 / TOLERANCE_ABS;
 let ix = (p.x * scale).round() as i64;
 let iy = (p.y * scale).round() as i64;
 let iz = (p.z * scale).round() as i64;
 // FNV-1a style hash
 let mut h: u64 = 14695981039346656037;
 for v in [ix, iy, iz] {
 h ^= v as u64;
 h = h.wrapping_mul(1099511628211);
 }
 h
}

/// Annotate a `BooleanHistory` with per-edge and per-vertex origins by
/// matching result BRep positions against the DS vertex/edge pool.
///
/// Both `edge_origins` and `vertex_origins` are filled in-place.
/// OCCT PostTreat equivalent: builds shape-to-origin maps for history tracking.
///
/// OCCT ref: BOPAlgo_Builder_3.cxx  ?`BOPAlgo_Builder::PostTreat`
/// (L1-250: builds `myLocModified` and `myLocGenerated` maps from DS images).
///
/// OCCT PostTreat algorithm (line-by-line mapping):
/// L20-40:  For each original shape, iterate sub-shapes (vertices, edges, faces).
/// L42-80:  Check `myImages[ei]` on each edge  ?if non-empty, record as Modified.
/// L82-110: For edges without images but present in result  ?record as Preserved.
/// L112-130: Generated edges (intersection edges)  ?record in myGenerated.
/// L132-170: For faces, check if wire edges were split  ?Modified; if not in
/// result  ?IsDeleted.
/// L172-200: Generated faces  ?myGenerated.
/// L202-230: Vertex tracking (fromA/fromB/intersection).
/// L232-250: Compute IsDeleted for entities absent from the result shape.
///
/// Differences from OCCT PostTreat:
/// - OCCT's PostTreat builds two maps: *myLocModified* (original -> last-modified
/// shape, for tracking splits and merges) and *myLocGenerated* (original -> list of
/// generated sub-shapes).  rcad's `annotate_history_from_ds` builds a simpler
/// `BooleanHistory` with flat `VertexOrigin`/`EdgeOrigin` arrays indexed by result
/// BRep position.
/// - OCCT PostTreat processes vertices, edges, and faces by iterating the DS images
/// (`myImages`, `myOrigins`, `myShapesSD`) and copying images from the source DS.
/// rcad uses spatial proximity (vertex point comparison) to match result vertices
/// to DS vertices, then traces edge origin from matched endpoints.
/// - OCCT PostTreat sets `myModified` for faces that were split (maps old -> new faces
/// via `myImages`).  rcad builds `FaceOrigin` separately (in `aggregate_face_origin`).
/// - OCCT PostTreat is called once at the end of `BOPAlgo_Builder::Build`.  rcad calls
/// `annotate_history_from_ds` inside `boolean_op_with_retry` after result assembly.
///
/// See also `BooleanHistory::update_with_post_treat()` for a more OCCT-aligned
/// implementation that uses `ds.my_images` instead of spatial proximity.
///
///  ?OCCT-aligned: core concept (history tracking from DS) matches OCCT's
/// image-map-based approach, adapted for rcad's flat-array data model.
///  ?OCCT-aligned: TopExp::MapShapes(myShape, myMapShape)  ?build result S index map.
/// OCCT maps TopoDS_Shape  ?identity for myMapShape lookup.
/// rcad: maps result vertex index  ?DS vertex index, result edge index  ?(DS vertices).
/// Used by PrepareHistory to determine Modified/Generated/Deleted provenance.
#[allow(dead_code)]
pub(crate) fn map_result_shapes(brep: &topods::BRep, ds: &DS) -> (Vec<usize>, Vec<(usize, usize)>) {
 // Collect flat vertex list from topods in ShapeRef.index order
 let topo_vertices: Vec<DVec3> = brep.tshapes.iter()
 .filter_map(|ts| match &**ts { topods::TShape::Vertex(v) => Some(v.point), _ => None })
 .collect();
 let mut result_to_ds: Vec<usize> = vec![usize::MAX; topo_vertices.len()];
 for (ri, pt) in topo_vertices.iter().enumerate() {
 for (di, dv) in ds.vertices.iter().enumerate() {
 if (dv.point - *pt).length_squared() < crate::tolerance::TOLERANCE_ABS * crate::tolerance::TOLERANCE_ABS * 4.0 {
 result_to_ds[ri] = di;
 break;
 }
 }
 }
 // Edge pairs from topods edges: map ShapeRef.index -> flat position
 let topo_edges: Vec<(usize, usize)> = brep.tshapes.iter()
 .filter_map(|ts| match &**ts { topods::TShape::Edge(e) => Some((e.first.index, e.last.index)), _ => None })
 .collect();
 let edge_pairs: Vec<(usize, usize)> = topo_edges.iter()
 .map(|&(s, e)| {
 let ds_s = result_to_ds.get(s).copied().unwrap_or(usize::MAX);
 let ds_e = result_to_ds.get(e).copied().unwrap_or(usize::MAX);
 (ds_s, ds_e)
 })
 .collect();
 (result_to_ds, edge_pairs)
}

///  ?OCCT-aligned: PrepareHistory (Builder_4.cxx L164-252).
/// OCCT iterates source shapes  ?LocModified  ?AddModified / AddGenerated / Remove.
/// rcad: uses pre-built result_to_ds map to annotate vertex/edge provenance.
#[allow(dead_code)]
pub(crate) fn annotate_history_from_ds(brep: &topods::BRep, history: &mut BooleanHistory, ds: &DS) {
 let (result_to_ds, _) = map_result_shapes(brep, ds);

 // OCCT L176: MapShapes done.  Annotate vertex origins (FromA/FromB/Intersection).
 let a_vc = ds.a_vertex_count;
 let n_result_verts = brep.tshapes.iter().filter(|ts| std::matches!(ts.as_ref(), topods::TShape::Vertex(_))).count();
 let mut vertex_origins: Vec<VertexOrigin> = Vec::with_capacity(n_result_verts);
 for ri in 0..n_result_verts {
 let di = result_to_ds[ri];
 let origin = if di == usize::MAX {
 VertexOrigin::Intersection
 } else if di < a_vc {
 VertexOrigin::FromA(di)
 } else {
 VertexOrigin::FromB(di - a_vc)
 };
 vertex_origins.push(origin);
 }
 history.vertex_origins = vertex_origins;

 // --- edge origins ---
 let a_vc = ds.a_vertex_count;
 let n_result_edges = brep.tshapes.iter().filter(|ts| std::matches!(ts.as_ref(), topods::TShape::Edge(_))).count();
 let mut edge_origins: Vec<EdgeOrigin> = Vec::with_capacity(n_result_edges);
 let a_ec = ds.a_edge_count;
 let total_ds_edges = ds.edges.len();

 for ts in &brep.tshapes {
 if let topods::TShape::Edge(ed) = &**ts {
 let ds_s = result_to_ds.get(ed.first.index).copied().unwrap_or(usize::MAX);
 let ds_e = result_to_ds.get(ed.last.index).copied().unwrap_or(usize::MAX);

 let origin = if ds_s == usize::MAX || ds_e == usize::MAX {
 EdgeOrigin::Generated
 } else if ds_s < a_vc && ds_e < a_vc {
 // Both endpoints are A vertices =look for a DS edge in A range.
 let found = (0..a_ec.min(total_ds_edges)).find(|&dei| {
 let de = &ds.edges[dei];
 (de.start_vertex == ds_s && de.end_vertex == ds_e)
 || (de.start_vertex == ds_e && de.end_vertex == ds_s)
 }); match found {
 Some(dei) => EdgeOrigin::FromA(dei),
 None => EdgeOrigin::SplitFromA(ds_s.min(a_vc - 1)),
 }
 } else if ds_s >= a_vc && ds_e >= a_vc {
 // Both endpoints are B vertices =look for a DS edge in B range.
 let found = (a_ec..total_ds_edges).find(|&dei| {
 let de = &ds.edges[dei];
 (de.start_vertex == ds_s && de.end_vertex == ds_e)
 || (de.start_vertex == ds_e && de.end_vertex == ds_s)
 }); match found {
 Some(dei) => EdgeOrigin::FromB(dei - a_ec),
 None => EdgeOrigin::SplitFromB(ds_s.min(ds.vertices.len().saturating_sub(1)) - a_vc),
 }
 } else {
 EdgeOrigin::Generated
 };
 edge_origins.push(origin);
 }
 }
 history.edge_origins = edge_origins;
}

pub(crate) fn aggregate_face_region_origin(face_origins: &[FaceOrigin]) -> ShellOrigin {
 let mut has_a = false;
 let mut has_b = false;
 let mut has_generated = false;
 for origin in face_origins {
 match origin {
 FaceOrigin::FromA(_) => has_a = true,
 FaceOrigin::FromB(_) => has_b = true,
 FaceOrigin::Generated => has_generated = true,
 }
 }

 match (has_a, has_b, has_generated) {
 (true, false, false) => ShellOrigin::FromA,
 (false, true, false) => ShellOrigin::FromB,
 (false, false, true) => ShellOrigin::Generated,
 _ => ShellOrigin::Mixed,
 }
}

pub(crate) fn aggregate_shell_region_origin(shell_origins: &[ShellOrigin]) -> SolidOrigin {
 let mut has_a = false;
 let mut has_b = false;
 let mut has_generated = false;
 let mut has_mixed = false;
 for origin in shell_origins {
 match origin {
 ShellOrigin::FromA => has_a = true,
 ShellOrigin::FromB => has_b = true,
 ShellOrigin::Generated => has_generated = true,
 ShellOrigin::Mixed => has_mixed = true,
 }
 }

 if has_mixed {
 return SolidOrigin::Mixed;
 }

 match (has_a, has_b, has_generated) {
 (true, false, false) => SolidOrigin::FromA,
 (false, true, false) => SolidOrigin::FromB,
 (false, false, true) => SolidOrigin::Generated,
 _ => SolidOrigin::Mixed,
 }
}

///  ?OCCT-aligned: PrepareHistory shell/solid provenance (Builder_4.cxx L164-252).
/// OCCT iterates source shapes  ?LocModified  ?AddModified/AddGenerated/Remove.
/// rcad: aggregates per-face origins to shell/solid level via face_region  ?shell  ?solid.
pub(crate) fn annotate_shell_and_solid_history(brep: &topods::BRep, history: &mut BooleanHistory) {
 let mut face_cursor = 0;
 let mut shell_origins = Vec::new();
 let mut solid_origins = Vec::new();

 for ts in &brep.tshapes {
 if let topods::TShape::Solid(sd) = &**ts {
 let solid_shell_start = shell_origins.len();
 for shell_sr in &sd.shells {
 if let topods::TShape::Shell(shd) = &*brep.tshapes[shell_sr.index] {
 let shell_face_count = shd.faces.len();
 let shell_face_origins = history
 .face_origins
 .get(face_cursor..face_cursor + shell_face_count)
 .unwrap_or(&[]);
 shell_origins.push(aggregate_face_region_origin(shell_face_origins));
 face_cursor += shell_face_count;
 }
 }
 solid_origins.push(aggregate_shell_region_origin(&shell_origins[solid_shell_start..]));
 }
 }

 if face_cursor != history.face_origins.len() {
 // Face count mismatch: BRep has more/fewer faces than history tracks.
 // This happens when compound reconstruction adds/removes faces or when
 // the face order in BRep differs from the emission order.  OCCT's
 // history tracking works with TopoDS shape identity  ?rcad's index-based
 // tracking is inherently more fragile.  Pad shell_origins to match.
 eprintln!("[HISTORY] face_cursor={} != history={}",
 face_cursor, history.face_origins.len());
 }
 history.shell_origins = shell_origins;
 history.solid_origins = solid_origins;
}

/// Boolean result builder (OCCT: BOPAlgo_BOP).
/// Tracks face splice origins and participates in `BooleanHistory`.
pub struct BooleanBuilder<'a> {
 ds: &'a DS,
 op: BooleanOpType,
 use_glue: bool,
 glue_tolerance: f64,
 context: RefCell<Context>,
 //  ?OCCT-aligned: error tracking (myReport / HasErrors equivalent).
 has_errors: bool,
 //  ?OCCT-aligned: myImages  ?source shape index  ?list of split image indices.
 // Uses RefCell because phase functions take &self (OCCT uses mutable member maps).
 my_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 //  ?OCCT-aligned: myOrigins  ?split shape index  ?list of source origin indices.
 my_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 //  ?OCCT-aligned: myShapesSD  ?source shape index  ?same-domain shape index.
 my_shapes_sd: std::cell::RefCell<std::collections::HashMap<usize, usize>>,
 //  ?OCCT-aligned: split edges created by FillImagesEdges (PaveBlock  ?new DSEdge).
 // Stored here because DS is immutable (rcad uses &'a DS); their indices start
 // at ds.edges.len() and are referenced by my_images(EDGE) / my_origins(EDGE).
 split_edges: std::cell::RefCell<Vec<crate::bopds::ds::DSEdge>>,
 //  ?OCCT-aligned: myInParts  ?source solid index  ?list of its IN face indices
 // (BOPAlgo_Builder.hxx L502).  Populated during FillImagesFaces, used by
 // FillIn3DParts / BuildDraftSolid for solid assembly.
 my_in_parts: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 //  ?OCCT-aligned: solid-level image tracking (BOPAlgo_Builder.hxx L498 myImages).
 // OCCT BuildSplitSolids stores split solids in myImages[source_solid].
 // rcad: maps source side (0=A, 1=B)  ?result solid indices from
 // build_split_solids.  Used by annotate_shell_and_solid_history and
 // for OCCT-form history tracking.
 my_solid_images: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 //  ?OCCT-aligned: solid-level origin tracking (BOPAlgo_Builder.hxx L500 myOrigins).
 // Reverse map: result solid index  ?list of source sides.
 my_solid_origins: std::cell::RefCell<std::collections::HashMap<usize, Vec<usize>>>,
 //  ?OCCT-aligned: myNonDestructive (BOPAlgo_Builder.hxx L503).
 // Safe processing  ?avoids modifying input shapes. Used in PostTreat.
 my_non_destructive: bool,
 //  ?OCCT-aligned: myCheckInverted (BOPAlgo_Builder.hxx L505).
 // Enables/disables inverted-solid check on input shapes.
 my_check_inverted: bool,
}

// =============================================================================
// OCCT 1:1  ? IsInternalFace (BOPTools_AlgoTools.cxx L791-872)
// =============================================================================

///  ?OCCT-aligned:  ?MEF (Map Edge= aces) = 椤?椤?閳??
/// OCCT BOPAlgo_FillIn3DParts::MapEdgesAndFaces (BOPAlgo_Tools.cxx L1479-1503)
/// OCCT-aligned: IsTangentFace (BOPTools_AlgoTools).
/// Checks if two faces are tangent (parallel normals + close distance).
pub fn is_tangent_face(fi_a: usize, fi_b: usize, ds: &crate::bopds::ds::DS, angle_tol: f64, dist_tol: f64) -> bool {
 let face_a = &ds.faces[fi_a];
 let face_b = &ds.faces[fi_b];
 let n_dot = face_a.normal.dot(face_b.normal).abs();
 if n_dot < angle_tol.cos() { return false; }
 let sample_a = if !face_a.boundary_verts.is_empty() {
 ds.vertices[face_a.boundary_verts[0]].point
 } else { return false; };
 let dist = match &face_b.surface {
 rcad_kernel::geom::Surface3::Plane(p) => (sample_a - p.origin).dot(p.normal).abs(),
 rcad_kernel::geom::Surface3::Sphere(s) => ((sample_a - s.center).length() - s.radius).abs(),
 _ => return false,
 };
 dist < dist_tol
}

pub(crate) fn build_edge_bounds(face_indices: &[usize], ds: &DS) -> std::collections::BTreeSet<usize> {
 let mut bounds: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
 for &fi in face_indices {
 let face = &ds.faces[fi];
 for &ei in &face.boundary_edges {
 bounds.insert(ei);
 }
 }
 bounds
}

///  ?OCCT-aligned: PointInFace  椤??= ?FaceSampleData =UV domain  = ?
/// OCCT BOPTools_AlgoTools3D.cxx L885-917
///
/// rcad  ? FaceSampleData  ?uv_domain =uv_centroid,= ?UV centroid
///  閿?=(OCCT =Hatcher =2D point-in-face, ?rcad =FaceSampleData
///  椤?椤? ?UV centroid = ? ?
// (point_in_face, classify_by_off_solid_edge removed  ?dead after ComputeState alignment)

/// = =3D  閿??u64 key,= 椤掑倵鍋??
pub(crate) fn quantize_pos(p: DVec3, tolerance: f64) -> u64 {
 let scale = 1.0 / tolerance;
 let x = (p.x * scale).round() as i64;
 let y = (p.y * scale).round() as i64;
 let z = (p.z * scale).round() as i64;
 //  = ?u64
 let xb = (x as u64) & 0x3FFFFF;
 let yb = (y as u64) & 0x3FFFFF;
 let zb = (z as u64) & 0x3FFFFF;
 (xb << 42) | (yb << 21) | zb
}

///  ?OCCT-aligned: IsInternalFace  椤??(BOPTools_AlgoTools.cxx L791-872)
///
///  椤?
/// Level 1:  椤?= 閳?=solid  閿?椤? 1  = 閵??
/// 椤? 閳╁唭? ?solid = ?
/// Level 2: ComputeState == =  solid  閿?= 椤愩儺鍞??
/// = ?PointInFace =classify_point ?
///
///  閳?? Some(true) =  ?solid = ?(IN)
/// Some(false) =  ?solid = ?(OUT)
/// None =  
/// Check if a DS vertex lies on the boundary edge between sv/ev, and if so add it
/// to split_verts with its parametric position t.
///  ?OCCT-aligned: FillImagesEdges checks pave blocks per edge (global scope).
pub(crate) fn check_and_add_split_vertex(
 ds: &DS,
 sv: usize,
 ev: usize,
 vi: usize,
 p_a: DVec3,
 ab: DVec3,
 ab_len2: f64,
 split_verts: &mut Vec<(usize, f64)>,
) {
 if vi == sv || vi == ev {
 return;
 }
 let p = ds.vertices[vi].point;
 let ap = p - p_a;
 let t = ap.dot(ab) / ab_len2;
 if t > 1e-8 && t < 1.0 - 1e-8 {
 let proj = p_a + ab * t;
 if (p - proj).length_squared() < 1e-10 {
 split_verts.push((vi, t));
 }
 }
}

///  ?OCCT-aligned: BuildSplitFaces edge assembly (L357-489) + DoSplitSEAMOnFace (L58-227).
pub(crate) fn collect_face_edge_segments(ds: &DS, face_idx: usize, pcurve_lookup: &impl Fn(usize) -> Option<Curve2d>) -> Vec<WireSegment> {
 let face = &ds.faces[face_idx];
 let mut segments: Vec<WireSegment> = Vec::new();
 let mut processed_seam_ds_edges: std::collections::HashSet<usize> = std::collections::HashSet::new();

 //  ?OCCT-aligned: boundary vertex position map (ShapesSD equivalent).
 // OCCT's DS shares vertices via ShapesSD during PaveFiller.
 // rcad: vertex remapping is done in make_section_edges_from_curve_pbs,
 // so IC endpoints already reference canonical vertices by this point.

 // Check if surface is closed (U/V)  for seam edge detection
 // OCCT L383-388: GeomLib::IsClosed  U/V
 let (is_u_closed, is_v_closed) = match &face.surface {
 Surface3::Sphere(_) => (true, true),
 Surface3::Cylinder(_) => (true, false),
 Surface3::Cone(_) => (true, false),
 _ => (false, false),
 };

 // ================================================================
 // 1. Original boundary edges (OCCT L357-460)
 // ================================================================
 // OCCT-aligned: orient boundary edges consistently for closed loop.
 // OCCT's TopExp_Explorer returns edges with the orientation they have
 // in the face's wire  ?each edge's end vertex matches the next edge's
 // start vertex.  rcad DS stores edges with arbitrary orientation.
 // Without this fix, a box face may have boundary edges like [2 ?, 3 ?,
 // 6 ?, 2 ?] where BOTH 3 ? and 6 ? end at vertex 7 (no outgoing edge
 // from 7), making the SmartMap connectivity wrong and preventing the
 // wire splitter from forming closed loops (fi=3 was failing).
 let mut prev_end: Option<usize> = None;
 for &ei in &face.boundary_edges {
 let edge = &ds.edges[ei];
 let (sv, ev) = match prev_end {
 Some(pe) if edge.start_vertex == pe => (edge.start_vertex, edge.end_vertex),
 Some(pe) if edge.end_vertex == pe => (edge.end_vertex, edge.start_vertex),
 _ => (edge.start_vertex, edge.end_vertex),
 };
 prev_end = Some(ev);

 //  ?OCCT L369: check if edge was split by intersection (myImages.IsBound).
 let edge_is_split = ei < ds.my_images.len() && ds.my_images[ei].len() > 1;

 if !edge_is_split {
 //  ?OCCT L395-404: seam detection for unsplit edges on periodic surfaces.
 // OCCT iterates all wire edges uniformly (no split/unsplit distinction);
 // rcad processes unsplit edges here  ?must detect seam before adding.
 let b_is_degenerated = ds.is_edge_degenerated(ei);
 let b_is_seam = !b_is_degenerated && (is_u_closed || is_v_closed)
 && ds.edge_on_face(ei, face_idx).map_or(false, |rep| {
 let (is_uiso, is_v_iso) = is_edge_isoline(&rep.pcurve, rep.pcurve_range);
 (is_u_closed && is_uiso) || (is_v_closed && is_v_iso)
 });
 if b_is_seam {
 if matches!(face.surface, Surface3::Sphere(_)) {
 if !processed_seam_ds_edges.insert(ei) { continue; }
 segments.extend(build_sphere_seam_segments(ds, ei, sv, ev, face, face_idx));
 } else {
 segments.extend(build_cylinder_seam_segments(ds, ei, sv, ev, face));
 }
 continue;
 }
 //  ?OCCT L371-382: unsplit edge  ?add directly.
 // OCCT L371-377: INTERNAL orientation  ?FWD+REV.
 // OCCT L379-381: FORWARD/REVERSED  ?add with orientation.
 let is_internal = ds.edges[ei].is_internal;
 let rep = ds.edge_on_face(ei, face_idx);
 let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
 Some(&edge.curve), Some(edge.t_range));
 let src = WireEdgeSource::DsEdge(ei);
 if is_internal {
 // OCCT L373-377: INTERNAL unsplit  ?FWD + REV
 segments.push(WireSegment {
 start_vertex: sv, end_vertex: ev, source: src.clone(),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
 });
 segments.push(WireSegment {
 start_vertex: ev, end_vertex: sv, source: src,
 orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: None, t_range: [0.0, 1.0],
 });
 } else {
 // OCCT L379-381: non-INTERNAL  ?add with orientation
 segments.push(WireSegment {
 start_vertex: sv, end_vertex: ev, source: src,
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
 });
 }
 continue;
 }

 //  ?OCCT L395-404: bIsClosed via IsClosed + IsEdgeIsoline.
 // On U-closed periodic surfaces (Sphere, Cylinder, Cone), seam edges
 // are U-isolines. Vertex coincidence NOT required (sphere pole-to-pole).
 let b_is_degenerated = ds.is_edge_degenerated(ei);
 let b_is_seam = if !b_is_degenerated && (is_u_closed || is_v_closed) {
 if let Some(rep) = ds.edge_on_face(ei, face_idx) {
 let (is_uiso, is_v_iso) = is_edge_isoline(&rep.pcurve, rep.pcurve_range);
 (is_u_closed && is_uiso) || (is_v_closed && is_v_iso)
 } else {
 false
 }
 } else {
 false
 };

 //  ?OCCT L408-464: iterate split sub-edges (aLIE from myImages.Find).
 if b_is_degenerated {
 // OCCT L413-417: iterate sub-edges, set orientation, append
 for &sub_ei in &ds.my_images[ei] {
 let sub_edge = &ds.edges[sub_ei];
 let sv_seg = sub_edge.start_vertex;
 let ev_seg = sub_edge.end_vertex;
 if sv_seg == ev_seg { continue; }
 segments.push(WireSegment {
 start_vertex: sv_seg, end_vertex: ev_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Forward,
 is_closed_on_face: true, second_pcurve: None, first_pcurve: None,
 t_range: [0.0, 1.0],
 });
 }
 } else if b_is_seam && matches!(face.surface, Surface3::Sphere(_)) {
 if !processed_seam_ds_edges.insert(ei) { continue; }
 segments.extend(build_sphere_seam_segments(ds, ei, sv, ev, face, face_idx));
 } else if b_is_seam {
 segments.extend(build_cylinder_seam_segments(ds, ei, sv, ev, face));
 } else {
 //  ?OCCT-aligned L408-464: three-branch split edge processing.
 // For each sub-edge from my_images, after degenerated (handled above):
 // 1. INTERNAL original (L420-426) -> FWD+REV
 // 2. Seam bIsClosed (L429-455) -> FWD+REV with fence (handled above)
 // 3. Normal (L457-462) -> orientation + IsSplitToReverseWithWarn
 if !ds.my_images.is_empty() && ei < ds.my_images.len() && !ds.my_images[ei].is_empty() {
 let is_original_internal = ds.edges[ei].is_internal;
 for &sub_ei in &ds.my_images[ei] {
 let sub_edge = &ds.edges[sub_ei];
 let sv_seg = sub_edge.start_vertex;
 let ev_seg = sub_edge.end_vertex;
 if sv_seg == ev_seg { continue; }
 // OCCT L420-426: INTERNAL original -> each sub-edge FWD+REV
 if is_original_internal {
 let (t_start, t_end) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
 Some(&sub_edge.curve), Some(sub_edge.t_range));
 let rep = ds.edge_on_face(sub_ei, face_idx)
 .or_else(|| ds.edge_on_face(ei, face_idx));
 segments.push(WireSegment {
 start_vertex: sv_seg, end_vertex: ev_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 segments.push(WireSegment {
 start_vertex: ev_seg, end_vertex: sv_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: None, t_range: [0.0, 1.0],
 });
 continue;
 }
 // OCCT L457-462: normal split -> orientation + IsSplitToReverseWithWarn
 let needs_reverse = is_split_to_reverse(ds, sub_ei, ei);
 let (t_fwd, t_rev) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
 Some(&sub_edge.curve), Some(sub_edge.t_range));
 let rep = ds.edge_on_face(sub_ei, face_idx)
 .or_else(|| ds.edge_on_face(ei, face_idx));
 if needs_reverse {
 segments.push(WireSegment {
 start_vertex: ev_seg, end_vertex: sv_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 } else {
 segments.push(WireSegment {
 start_vertex: sv_seg, end_vertex: ev_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 } }
 } else {
 // OCCT: edge not split 閳?add as single segment (Builder_2.cxx L374-378).
 let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
 Some(&ds.edges[ei].curve), Some(ds.edges[ei].t_range));
 let rep = ds.edge_on_face(ei, face_idx);
 segments.push(WireSegment {
 start_vertex: sv, end_vertex: ev,
 source: WireEdgeSource::DsEdge(ei),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 }
 }
 }

 // ================================================================
 //  ?OCCT-aligned: inner wire edges  ?same processing as outer boundary.
 // OCCT TopExp_Explorer iterates all wires' edges in one loop.
 // rcad stores them separately, so we apply identical logic here.
 // ================================================================
 for inner_wire in &face.inner_boundary_edges {
 for &(ei, forward_in_wire) in inner_wire {
 let edge = &ds.edges[ei];
 let (sv, ev) = if forward_in_wire {
 (edge.start_vertex, edge.end_vertex)
 } else {
 (edge.end_vertex, edge.start_vertex)
 };
 if sv == ev { continue; }

 let edge_is_split = ei < ds.my_images.len() && ds.my_images[ei].len() > 1;

 if !edge_is_split {
 let is_internal = ds.edges[ei].is_internal;
 let rep = ds.edge_on_face(ei, face_idx);
 let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
 Some(&edge.curve), Some(edge.t_range));
 let src = WireEdgeSource::DsEdge(ei);
 if is_internal {
 segments.push(WireSegment {
 start_vertex: sv, end_vertex: ev, source: src.clone(),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
 });
 segments.push(WireSegment {
 start_vertex: ev, end_vertex: sv, source: src,
 orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: None, t_range: [0.0, 1.0],
 });
 } else {
 segments.push(WireSegment {
 start_vertex: sv, end_vertex: ev, source: src,
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or(edge.t_range),
 });
 }
 continue;
 }

 let b_is_degenerated = ds.is_edge_degenerated(ei);
 let b_is_seam = if !b_is_degenerated && (is_u_closed || is_v_closed) {
 if let Some(rep) = ds.edge_on_face(ei, face_idx) {
 let (is_uiso, is_v_iso) = is_edge_isoline(&rep.pcurve, rep.pcurve_range);
 (is_u_closed && is_uiso) || (is_v_closed && is_v_iso)
 } else {
 false
 }
 } else {
 false
 };

 if b_is_degenerated {
 // OCCT L413-417: iterate sub-edges, set orientation, append
 for &sub_ei in &ds.my_images[ei] {
 let sub_edge = &ds.edges[sub_ei];
 let sv_seg = sub_edge.start_vertex;
 let ev_seg = sub_edge.end_vertex;
 if sv_seg == ev_seg { continue; }
 segments.push(WireSegment {
 start_vertex: sv_seg, end_vertex: ev_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Forward,
 is_closed_on_face: true, second_pcurve: None, first_pcurve: None,
 t_range: [0.0, 1.0],
 });
 }
 } else if b_is_seam && matches!(face.surface, Surface3::Sphere(_)) {
 if !processed_seam_ds_edges.insert(ei) { continue; }
 segments.extend(build_sphere_seam_segments(ds, ei, sv, ev, face, face_idx));
 } else if b_is_seam {
 segments.extend(build_cylinder_seam_segments(ds, ei, sv, ev, face));
 } else {
 if !ds.my_images.is_empty() && ei < ds.my_images.len() && !ds.my_images[ei].is_empty() {
 let is_original_internal = ds.edges[ei].is_internal;
 for &sub_ei in &ds.my_images[ei] {
 let sub_edge = &ds.edges[sub_ei];
 let sv_seg = sub_edge.start_vertex;
 let ev_seg = sub_edge.end_vertex;
 if sv_seg == ev_seg { continue; }
 if is_original_internal {
 let (t_start, t_end) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
 Some(&sub_edge.curve), Some(sub_edge.t_range));
 let rep = ds.edge_on_face(sub_ei, face_idx)
 .or_else(|| ds.edge_on_face(ei, face_idx));
 segments.push(WireSegment {
 start_vertex: sv_seg, end_vertex: ev_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 segments.push(WireSegment {
 start_vertex: ev_seg, end_vertex: sv_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: None, t_range: [0.0, 1.0],
 });
 continue;
 }
 let needs_reverse = is_split_to_reverse(ds, sub_ei, ei);
 let (t_fwd, t_rev) = edge_uv_tangent(ds, sv_seg, ev_seg, &face.surface,
 Some(&sub_edge.curve), Some(sub_edge.t_range));
 let rep = ds.edge_on_face(sub_ei, face_idx)
 .or_else(|| ds.edge_on_face(ei, face_idx));
 if needs_reverse {
 segments.push(WireSegment {
 start_vertex: ev_seg, end_vertex: sv_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Reversed, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 } else {
 segments.push(WireSegment {
 start_vertex: sv_seg, end_vertex: ev_seg,
 source: WireEdgeSource::DsEdge(sub_ei),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 } }
 } else {
 let (t_start, t_end) = edge_uv_tangent(ds, sv, ev, &face.surface,
 Some(&ds.edges[ei].curve), Some(ds.edges[ei].t_range));
 let rep = ds.edge_on_face(ei, face_idx);
 segments.push(WireSegment {
 start_vertex: sv, end_vertex: ev,
 source: WireEdgeSource::DsEdge(ei),
 orientation: WireOrientation::Forward, is_closed_on_face: false, second_pcurve: None,
 first_pcurve: rep.map(|r| r.pcurve.clone()),
 t_range: rep.map(|r| r.pcurve_range).unwrap_or([0.0, 1.0]),
 });
 }
 }
 }  // end inner_wire edge loop
 }  // end inner_wire loop

 // ================================================================
 // IN edge PBs (OCCT BOPAlgo_Builder_2.cxx L467-480).
 // Each IN PaveBlock contributes its split edge as FWD+REV.
 let boundary_set: std::collections::HashSet<usize> =
 face.boundary_edges.iter().copied().collect();
 let mut pb_dedup: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for &pb_idx in &face.face_info.pave_blocks_in {
 if pb_idx >= ds.pave_blocks.len() { continue; }
 let pb = &ds.pave_blocks[pb_idx];
 if boundary_set.contains(&pb.0.read().unwrap().original_edge) { continue; }
 if !pb_dedup.insert(pb.0.read().unwrap().original_edge) { continue; }
 let ei = pb.0.read().unwrap().new_edge.unwrap_or(pb.0.read().unwrap().original_edge);
 if ei >= ds.edges.len() { continue; }
 let edge = &ds.edges[ei];
 let face_surf = &ds.faces[face_idx].surface;
 let t_start = edge_angle_2d(&edge.curve, edge.t_range[0], edge.t_range, face_surf, false, ds.vertices[edge.start_vertex].geom_tol);
 let t_end = edge_angle_2d(&edge.curve, edge.t_range[1], edge.t_range, face_surf, true, ds.vertices[edge.end_vertex].geom_tol);
 // OCCT: aLE.Append(aSp) with FORWARD orientation.
 segments.push(WireSegment {
 start_vertex: edge.start_vertex, end_vertex: edge.end_vertex,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
 is_closed_on_face: false, second_pcurve: None, first_pcurve: None,
 t_range: edge.t_range,
 });
 // OCCT: aLE.Append(aSp) with REVERSED orientation.
 segments.push(WireSegment {
 start_vertex: edge.end_vertex, end_vertex: edge.start_vertex,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Reversed,
 is_closed_on_face: false, second_pcurve: None, first_pcurve: None,
 t_range: edge.t_range,
 });
 }

 // Section edges from pave_blocks_sc (OCCT: aMSCPB entries from MakeBlocks).
 // These edges are NOT in boundary_edges 閳?they're registered only in
 // pave_blocks_sc.  OCCT adds each edge once to myShapes; the WireSplitter
 // determines the required orientation during loop walking.
 let mut sc_dedup: std::collections::HashSet<usize> = std::collections::HashSet::new();
 for &pb_idx in &face.face_info.pave_blocks_sc {
 if pb_idx >= ds.pave_blocks.len() { continue; }
 let pb = &ds.pave_blocks[pb_idx];
 let ei = pb.0.read().unwrap().new_edge.unwrap_or(pb.0.read().unwrap().original_edge);
 if ei >= ds.edges.len() { continue; }
 if boundary_set.contains(&ei) { continue; }
 if !sc_dedup.insert(ei) { continue; }
 if ds.is_edge_degenerated(ei) { continue; }
 let edge = &ds.edges[ei];
 let mut sv = edge.start_vertex;
 let mut ev = edge.end_vertex;
 let is_deg = sv == ev || (ds.vertices[sv].point - ds.vertices[ev].point).length_squared() < TOLERANCE_ABS_SQ;
 if is_deg { continue; }
 // OCCT: orient section edge so its end_vertex doesn't get starved
 // of outgoing edges.  The vertex connectivity from previously added
 // segments determines the direction that balances in/out degree.
 {
 let ev_in = segments.iter().filter(|s| s.end_vertex == ev).count();
 let ev_out = segments.iter().filter(|s| s.start_vertex == ev).count();
 if ev_in >= ev_out + 1 {
 std::mem::swap(&mut sv, &mut ev);
 }
 }
 let sec_pcurve = edge.face_reps.iter().find(|r| r.face_idx == face_idx).map(|r| r.pcurve.clone());
 segments.push(WireSegment {
 start_vertex: sv, end_vertex: ev,
 source: WireEdgeSource::DsEdge(ei), orientation: WireOrientation::Forward,
 is_closed_on_face: false, second_pcurve: None, first_pcurve: sec_pcurve,
 t_range: edge.t_range,
 });
 }
  segments
}

/// Check if a 2D point is inside a 2D polygon using ray casting (OCCT-aligned utility).
pub(crate) fn point_in_polygon_2d(poly: &[DVec2], pt: DVec2) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let vi = poly[i];
        let vj = poly[j];
        if ((vi.y > pt.y) != (vj.y > pt.y))
            && (pt.x < (vj.x - vi.x) * (pt.y - vi.y) / (vj.y - vi.y) + vi.x)
        {
            inside = !inside;
        }
        j = i;
    }
    inside
}
