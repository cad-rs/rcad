use crate::bopds::face_info::FaceInfo;
use crate::bopds::ds::types::*;
use crate::tolerance::*;
use rcad_kernel::topods;
use rcad_kernel::geom::{Curve2dEval, Curve3, Line3, Plane, Surface3};
use glam::DVec3;
use std::collections::HashMap;

/// Build DS from two topods::BRep shapes.
///
/// Analogous to [`DS::new_with_fuzzy`] but reads directly from the TopoDS-aligned
/// representation, bypassing the deprecated `rcad_kernel::BRep`.
pub fn new_from_topods(a: &topods::BRep, b: &topods::BRep, fuzzy_tol: f64) -> DS {
 let tol = fuzzy_tol.max(TOLERANCE_ABS);
 let mut ds = DS {
 vertices: Vec::new(),
 edges: Vec::new(),
 wires: Vec::new(),
 shells: Vec::new(),
 solids: Vec::new(),
 comp_solids: Vec::new(),
 faces: Vec::new(),
 interf_vv: Vec::new(),
 interf_ve: Vec::new(),
 interf_vf: Vec::new(),
 interf_ee: Vec::new(),
 interf_ef: Vec::new(),
 interf_ff: Vec::new(),
 interf_vz: Vec::new(),
 interf_ez: Vec::new(),
 interf_fz: Vec::new(),
 interf_zz: Vec::new(),
 intersection_curves: Vec::new(),
 ff_points: Vec::new(),
 section_edge_refs: Vec::new(),
 fuzzy_tol: tol,
 a_vertex_count: 0,
 a_edge_count: 0,
 a_face_count: 0,
 shared_topology: SharedTopologyInfo::default(),
 shape_sd: ShapeSD::new(0, &SharedTopologyInfo::default()),
 same_domain_overlaps: Vec::new(),
 common_blocks: Vec::new(),
 my_images: Vec::new(),
 my_origins: Vec::new(),
 wire_images: Vec::new(),
 shell_images: Vec::new(),
 solid_images: Vec::new(),
 pave_blocks: Vec::new(),
 locations: Vec::new(),
 increased_ss: std::collections::HashSet::new(),
 interf_tb: std::collections::HashSet::new(),
 map_ve: std::collections::HashMap::new(),
 shape_info: Vec::new(),
 nb_source_shapes: 0,
 };

 load_topods_brep(&mut ds, a, ShapeOrigin::ShapeA);
 ds.a_vertex_count = ds.vertices.len();
 ds.a_edge_count = ds.edges.len();
 ds.a_face_count = ds.faces.len();
 load_topods_brep(&mut ds, b, ShapeOrigin::ShapeB);
 ds.compute_uv_boundaries();
 ds.build_face_reps();
 ds.init_shape_info();
 ds
}

/// Load a topods::BRep into the DS.
///
/// Walks the tshape pool to extract vertices, edges, faces, shells, solids,
/// and compsolids, building the flat DS arrays with proper index mapping.
pub fn load_topods_brep(ds: &mut DS, brep: &topods::BRep, origin: ShapeOrigin) {
 // ----- Step 1: Collect Vertex TShapes -> DS vertices -----
 let mut v_map: HashMap<usize, usize> = HashMap::new();
 for (ti, ts) in brep.tshapes.iter().enumerate() {
 if let topods::TShape::Vertex(vd) = ts.as_ref() {
 let tol = vd.tolerance.max(TOLERANCE_ABS);
 if let Some(existing) = ds.find_vertex_near(vd.point, tol) {
 v_map.insert(ti, existing);
 } else {
 let vi = ds.vertices.len();
 ds.vertices.push(DSVertex {
 point: vd.point,
 origin: Some(origin),
 geom_tol: vd.tolerance,
 is_internal: false,
 location: 0,
 });
 v_map.insert(ti, vi);
 }
 }
 }

 // ----- Step 2: Collect Edge TShapes -> DS edges -----
 let mut e_map: HashMap<usize, usize> = HashMap::new();
 for (ti, ts) in brep.tshapes.iter().enumerate() {
 if let topods::TShape::Edge(ed) = ts.as_ref() {
 let start = v_map.get(&ed.first.index).copied().unwrap_or(0);
 let end = v_map.get(&ed.last.index).copied().unwrap_or(0);
 let curve = ed.curve.clone().unwrap_or_else(|| {
 let p0 = ds.vertices[start].point;
 let p1 = ds.vertices[end].point;
 let dir = (p1 - p0).normalize_or_zero();
 Curve3::Line(Line3 { origin: p0, direction: dir })
 });
 let t_range = ed.range;
 let mut vertex_params = std::collections::HashMap::new();
 for (&vi, &param) in &ed.vertex_params {
 if let Some(&ds_vi) = v_map.get(&vi) {
 vertex_params.insert(ds_vi, param);
 }
 }
 if !vertex_params.contains_key(&start) { vertex_params.insert(start, t_range[0]); }
 if !vertex_params.contains_key(&end) { vertex_params.insert(end, t_range[1]); }
 let is_geometric = ed.curve.is_some() && (matches!(&curve, Curve3::Line(_) | Curve3::Circle(_)
 | Curve3::Ellipse(_) | Curve3::BSpline(_) | Curve3::Bezier(_)));
 let ds_ei = ds.edges.len();
 ds.edges.push(DSEdge {
 start_vertex: start,
 end_vertex: end,
 curve,
 t_range,
 origin,
 geom_tol: ed.tolerance,
 paves: Vec::new(),
 pave_blocks: Vec::new(),
 face_reps: Vec::new(),
 is_internal: false,
 vertex_params,
 face_tolerances: Vec::new(),
 is_geometric,
 location: 0,
 });
 e_map.insert(ti, ds_ei);
 ds.init_pave_blocks_for_edge(ds_ei);
 }
 }

 // ----- Step 3: Build compsolid-solid map -----
 let mut cs_solid_map: HashMap<usize, usize> = HashMap::new();
 for ts in brep.tshapes.iter() {
 if let topods::TShape::CompSolid(cs_solids) = ts.as_ref() {
 for (pos, sr) in cs_solids.iter().enumerate() {
 cs_solid_map.insert(sr.index, pos);
 }
 }
 }

 // ----- Step 4: Hierarchy �?Solid -> Shell -> Face -> Wire -> Edge -----
 let mut face_flat_idx = 0usize;
 let mut shell_counter = 0usize;
 let mut solid_counter = 0usize;
 let mut s_map: HashMap<usize, usize> = HashMap::new();
 let mut f_map: HashMap<usize, usize> = HashMap::new();

 for (ti, ts) in brep.tshapes.iter().enumerate() {
 if let topods::TShape::Solid(sd) = ts.as_ref() {
 s_map.insert(ti, solid_counter);
 let compsolid_idx = cs_solid_map.get(&ti).copied();
 let shell_start = ds.shells.len();

 for shell_sr in &sd.shells {
 let shell_data = brep.shell(*shell_sr);
 let prev_face_count = ds.faces.len();

 for face_sr in &shell_data.faces {
 let face_data = brep.face(*face_sr);
 let surface = face_data.surface.clone().unwrap_or_else(|| {
 let origin = DVec3::ZERO;
 Surface3::Plane(Plane::new(origin, DVec3::Z))
 });
 let outer_wire_data = brep.wire(face_data.outer_wire);
 let boundary_edges_ordered = reorder_wire_topods(&outer_wire_data.edges, brep, &e_map);
 let boundary_edges: Vec<usize> = boundary_edges_ordered.iter().map(|&(ei, _)| ei).collect();
 let boundary_edge_forwards: Vec<bool> = boundary_edges_ordered.iter().map(|&(_, fwd)| fwd).collect();
 let boundary_verts = {
 let edges = &outer_wire_data.edges;
 if edges.is_empty() {
 Vec::new()
 } else if edges.len() == 1 {
 let ed = brep.edge(edges[0]);
 vec![
 v_map.get(&ed.first.index).copied().unwrap_or(0),
 v_map.get(&ed.last.index).copied().unwrap_or(0),
 ]
 } else {
 let mut verts = Vec::with_capacity(edges.len());
 for i in 0..edges.len() {
 let next_i = (i + 1) % edges.len();
 let e = brep.edge(edges[i]);
 let en = brep.edge(edges[next_i]);
 let shared = if e.first.index == en.first.index || e.first.index == en.last.index {
 e.first.index
 } else {
 e.last.index
 };
 let non_shared = if shared == e.first.index { e.last.index } else { e.first.index };
 verts.push(v_map.get(&non_shared).copied().unwrap_or(0));
 }
 verts
 }
 };
 let outer_wire_idx = Some(ds.wires.len());
 ds.wires.push(DSWire { edges: boundary_edges.clone() });
 let inner_boundary_edges: Vec<Vec<(usize, bool)>> = face_data.inner_wires.iter()
 .map(|iw_sr| {
 let iw_data = brep.wire(*iw_sr);
 iw_data.edges.iter()
 .map(|we_sr| (e_map.get(&we_sr.index).copied().unwrap_or(0), we_sr.orientation.is_forward()))
 .collect()
 })
 .collect();
 let inner_wire_idxs: Vec<usize> = (0..face_data.inner_wires.len())
 .map(|_| { let wi = ds.wires.len(); ds.wires.push(DSWire { edges: Vec::new() }); wi })
 .collect();
 for (ii, iw_sr) in face_data.inner_wires.iter().enumerate() {
 let iw_data = brep.wire(*iw_sr);
 ds.wires[inner_wire_idxs[ii]].edges = iw_data.edges.iter()
 .map(|we_sr| e_map.get(&we_sr.index).copied().unwrap_or(0))
 .collect();
 }
 let normal = match &surface {
 Surface3::Plane(p) => p.normal,
 _ => {
 if boundary_verts.len() >= 3 {
 let p0 = ds.vertices.get(boundary_verts[0]).map(|v| v.point).unwrap_or(DVec3::ZERO);
 let p1 = ds.vertices.get(boundary_verts[1]).map(|v| v.point).unwrap_or(DVec3::ZERO);
 let p2 = ds.vertices.get(boundary_verts[2]).map(|v| v.point).unwrap_or(DVec3::ZERO);
 (p1 - p0).cross(p2 - p0).normalize_or_zero()
 } else { DVec3::Z }
 }
 };
 let ds_fi = ds.faces.len();
 f_map.insert(face_sr.index, ds_fi);
 ds.faces.push(DSFace {
 surface,
 boundary_verts,
 boundary_edges,
 boundary_edge_forwards,
 inner_boundary_edges,
 outer_wire_idx,
 inner_wire_idxs,
 normal,
 origin,
 face_info: FaceInfo::default(),
 source_face_idx: face_flat_idx,
 geom_tol: face_data.tolerance,
 location: 0,
 uv_boundary: None,
 natural_restriction: face_data.natural_restriction,
 source_shell_idx: Some(shell_counter),
 source_solid_idx: Some(solid_counter),
 source_compsolid_idx: compsolid_idx,
 });
 face_flat_idx += 1;
 }
 let shell_face_idxs: Vec<usize> = (prev_face_count..ds.faces.len()).collect();
 if !shell_face_idxs.is_empty() {
 ds.shells.push(DSShell { faces: shell_face_idxs });
 }
 shell_counter += 1;
 }
 let solid_shells: Vec<usize> = (shell_start..ds.shells.len()).collect();
 if !solid_shells.is_empty() {
 ds.solids.push(DSSolid { shells: solid_shells });
 }
 solid_counter += 1;
 }
 }

 // ----- Step 5: CompSolids -----
 for (_, ts) in brep.tshapes.iter().enumerate() {
 if let topods::TShape::CompSolid(cs_solids) = ts.as_ref() {
 let ds_solid_indices: Vec<usize> = cs_solids.iter()
 .filter_map(|sr| s_map.get(&sr.index).copied())
 .collect();
 if !ds_solid_indices.is_empty() {
 ds.comp_solids.push(DSCompSolid { solids: ds_solid_indices });
 }
 }
 }

 // ----- Step 6: Transfer pcurves from topods edge data -----
 for (ti, ts) in brep.tshapes.iter().enumerate() {
 if let topods::TShape::Edge(ed) = ts.as_ref() {
 let Some(&ds_ei) = e_map.get(&ti) else { continue; };
 for (&face_ti, (pcurve, param_start, param_end)) in &ed.pcurves {
 let Some(&ds_fi) = f_map.get(&face_ti) else { continue; };
 if ds.edges[ds_ei].face_reps.iter().any(|r| r.face_idx == ds_fi) { continue; }
 let uv_start = pcurve.point_at(*param_start);
 let uv_end = pcurve.point_at(*param_end);
 let span = (uv_end - uv_start).length();
 if span < TOLERANCE_CLAMP_MIN || !span.is_finite() { continue; }
 ds.edges[ds_ei].face_reps.push(DSCurveRepOnFace {
 face_idx: ds_fi,
 pcurve: pcurve.clone(),
 pcurve2: None,
 pcurve_range: [*param_start, *param_end],
 start_param: *param_start,
 end_param: *param_end,
 });
 }
 for rep in &ed.representations {
 match rep {
 topods::CurveRepresentation::CurveOnSurface { face, pcurve, range } => {
 let Some(&ds_fi) = f_map.get(face) else { continue; };
 if ds.edges[ds_ei].face_reps.iter().any(|r| r.face_idx == ds_fi) { continue; }
 let span = range[1] - range[0];
 if span < TOLERANCE_CLAMP_MIN { continue; }
 ds.edges[ds_ei].face_reps.push(DSCurveRepOnFace {
 face_idx: ds_fi,
 pcurve: pcurve.clone(),
 pcurve2: None,
 pcurve_range: *range,
 start_param: range[0],
 end_param: range[1],
 });
 }
 topods::CurveRepresentation::CurveOnClosedSurface { face, pcurve1, pcurve2, range } => {
 let Some(&ds_fi) = f_map.get(face) else { continue; };
 if ds.edges[ds_ei].face_reps.iter().any(|r| r.face_idx == ds_fi) { continue; }
 let span = range[1] - range[0];
 if span < TOLERANCE_CLAMP_MIN { continue; }
 ds.edges[ds_ei].face_reps.push(DSCurveRepOnFace {
 face_idx: ds_fi,
 pcurve: pcurve1.clone(),
 pcurve2: Some(pcurve2.clone()),
 pcurve_range: *range,
 start_param: range[0],
 end_param: range[1],
 });
 }
 _ => {}
 }
 }
 }
 }
}

/// OCCT-aligned: reorder wire edges by traversal order (TopExp_Explorer).
/// TopoDS version: reads edge vertex adjacency from topods::BRep.
/// Returns (DS edge index, forward_in_wire) pairs in traversal order.
pub fn reorder_wire_topods(
 wire_edges: &[topods::ShapeRef],
 brep: &topods::BRep,
 e_map: &HashMap<usize, usize>,
) -> Vec<(usize, bool)> {
 if wire_edges.len() <= 1 {
 return wire_edges.iter().map(|we| (e_map.get(&we.index).copied().unwrap_or(0), we.orientation.is_forward())).collect();
 }
 let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
 for (i, we_sr) in wire_edges.iter().enumerate() {
 let ed = brep.edge(*we_sr);
 adj.entry(ed.first.index).or_default().push(i);
 adj.entry(ed.last.index).or_default().push(i);
 }
 let first_ed = brep.edge(wire_edges[0]);
 let mut cur = first_ed.first.index;
 let mut used = vec![false; wire_edges.len()];
 let mut ordered = Vec::with_capacity(wire_edges.len());
 for _ in 0..wire_edges.len() {
 let next_i = adj.entry(cur).or_default().iter().copied()
 .find(|&i| !used[i])
 .expect("wire is not closed -- broken topology");
 used[next_i] = true;
 let we_sr = &wire_edges[next_i];
 ordered.push((e_map.get(&we_sr.index).copied().unwrap_or(0), we_sr.orientation.is_forward()));
 let ed = brep.edge(*we_sr);
 cur = if ed.first.index == cur { ed.last.index } else { ed.first.index };
 }
 ordered
}
