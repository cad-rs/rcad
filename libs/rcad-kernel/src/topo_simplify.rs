use std::collections::HashMap;
use crate::{BRep, Edge, WireEdge, topology::Wire};

fn oriented_edge_vertices(brep: &BRep, we: WireEdge) -> Option<(usize, usize)> {
 let e = brep.edges.get(we.idx)?;
 if we.forward {
 Some((e.start, e.end))
 } else {
 Some((e.end, e.start))
 }
}

fn points_are_collinear_forward(
 a: glam::DVec3,
 b: glam::DVec3,
 c: glam::DVec3,
 linear_tol: f64,
) -> bool {
 let ab = b - a;
 let bc = c - b;
 let ab_len = ab.length();
 let bc_len = bc.length();
 if ab_len <= 1e-12 || bc_len <= 1e-12 {
 return false;
 }
 let cross = ab.cross(bc).length();
 let dot = ab.dot(bc);
 cross <= linear_tol * (ab_len + bc_len) && dot > 0.0
}

fn find_existing_edge_between_vertices(brep: &BRep, from: usize, to: usize) -> Option<WireEdge> {
 for (idx, e) in brep.edges.iter().enumerate() {
 if e.start == from && e.end == to {
 return Some(WireEdge::fwd(idx));
 }
 if e.start == to && e.end == from {
 return Some(WireEdge::rev(idx));
 }
 }
 None
}

fn push_topology_edge(brep: &mut BRep, start: usize, end: usize) -> WireEdge {
 let idx = brep.edges.len();
 brep.edges.push(Edge { start, end });

 if !brep.geom.edge_curve.is_empty() {
 brep.geom.edge_curve.push(None);
 }
 if !brep.geom.edge_pcurves.is_empty() {
 brep.geom.edge_pcurves.push(Vec::new());
 }
 if !brep.geom.edge_curve_range.is_empty() {
 brep.geom.edge_curve_range.push(None);
 }
 if !brep.geom.edge_degenerated.is_empty() {
 brep.geom.edge_degenerated.push(false);
 }
 if !brep.geom.edge_tolerance.is_empty() {
 brep.geom.edge_tolerance.push(0.0);
 }
 if !brep.geom.edge_same_parameter.is_empty() {
 brep.geom.edge_same_parameter.push(false);
 }
 if !brep.geom.edge_same_range.is_empty() {
 brep.geom.edge_same_range.push(false);
 }

 WireEdge::fwd(idx)
}

fn merge_collinear_consecutive_edges_in_wire(
 brep: &mut BRep,
 wire: &mut Vec<WireEdge>,
 linear_tol: f64,
) -> usize {
 if wire.len() < 4 {
 return 0;
 }

 let mut merged = 0usize;
 loop {
 let n = wire.len();
 if n < 4 {
 break;
 }

 let mut changed = false;
 for i in 0..n {
 let j = (i + 1) % n;
 let Some((u, v1)) = oriented_edge_vertices(brep, wire[i]) else {
 continue;
 };
 let Some((v2, w)) = oriented_edge_vertices(brep, wire[j]) else {
 continue;
 };
 if v1 != v2 || u == w {
 continue;
 }

 let Some(p_u) = brep.vertices.get(u).map(|v| v.point) else {
 continue;
 };
 let Some(p_v) = brep.vertices.get(v1).map(|v| v.point) else {
 continue;
 };
 let Some(p_w) = brep.vertices.get(w).map(|v| v.point) else {
 continue;
 };
 if !points_are_collinear_forward(p_u, p_v, p_w, linear_tol) {
 continue;
 }

 let bridge = find_existing_edge_between_vertices(brep, u, w)
 .unwrap_or_else(|| push_topology_edge(brep, u, w));

 if i + 1 < n {
 wire.splice(i..=i + 1, [bridge]);
 } else {
 wire.pop();
 wire.remove(0);
 wire.insert(0, bridge);
 }
 merged += 1;
 changed = true;
 break;
 }

 if !changed {
 break;
 }
 }
 merged
}

/// Merge fragmented collinear edges inside all face wires (outer + inner).
///
/// This pass is intentionally local/topological: it only collapses *consecutive*
/// collinear edge pairs in a wire into a direct bridge edge, preserving wire
/// order and loop closure.
///
/// Returns the number of edge-pair merges performed.
pub fn merge_collinear_edges_in_wires(brep: &mut BRep, linear_tol: f64) -> usize {
 let tol = linear_tol.max(1e-12);
 let mut merged_total = 0usize;

 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 for fi in 0..brep.solids[si].shells[shi].faces.len() {
 let mut outer = brep.solids[si].shells[shi].faces[fi].outer_wire.edges.clone();
 merged_total += merge_collinear_consecutive_edges_in_wire(brep, &mut outer, tol);
 brep.solids[si].shells[shi].faces[fi].outer_wire.edges = outer;

 let inner_count = brep.solids[si].shells[shi].faces[fi].inner_wires.len();
 for wi in 0..inner_count {
 let mut inner =
 brep.solids[si].shells[shi].faces[fi].inner_wires[wi].edges.clone();
 merged_total += merge_collinear_consecutive_edges_in_wire(brep, &mut inner, tol);
 brep.solids[si].shells[shi].faces[fi].inner_wires[wi].edges = inner;
 }
 }
 }
 }

 merged_total
}

/// ✅ OCCT : BRep — OCCT BuildSolid + MakeBlocks  ,
/// rcad build_split_edges  , 。
///
/// (OCCT BOPAlgo_BuilderSolid / BOPAlgo_Tools::MakeBlocks):
/// 1. vertex → edges  
/// 2. degree=2 vertex ( 2 ):
/// a. (same edge_curve index)
/// b. Line  :  
/// c.  :  ,  wire  , 
/// 3.  
///
///  。
pub fn merge_collinear_brep_edges(brep: &mut BRep, linear_tol: f64) -> usize {
 let tol = linear_tol.max(1e-12);
 let mut total_merged = 0usize;

 //  
 loop {
 // 1. vertex → edges  
 let mut vert_edges: HashMap<usize, Vec<usize>> = HashMap::new();
 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 for fi in 0..brep.solids[si].shells[shi].faces.len() {
 let face = &brep.solids[si].shells[shi].faces[fi];
 let mut collect = |wire: &Wire| {
 for we in &wire.edges {
 if let Some(e) = brep.edges.get(we.idx) {
 vert_edges.entry(e.start).or_default().push(we.idx);
 vert_edges.entry(e.end).or_default().push(we.idx);
 }
 }
 };
 collect(&face.outer_wire);
 for w in &face.inner_wires {
 collect(w);
 }
 }
 }
 }

 let mut merged_any = false;

 // 2. degree=2 vertex  
 let mut to_merge: Vec<(usize, usize, usize)> = Vec::new(); // (keep_ei, remove_ei, shared_v)

 for (vi, edges) in &vert_edges {
 if edges.len() != 2 { continue; }
 let e1 = edges[0];
 let e2 = edges[1];
 if e1 == e2 { continue; }

 //  
 let c1 = brep.geom.edge_curve.get(e1).copied().flatten();
 let c2 = brep.geom.edge_curve.get(e2).copied().flatten();
 let same_curve = match (c1, c2) {
 (Some(ci1), Some(ci2)) => ci1 == ci2,
 _ => false,
 };
 if !same_curve { continue; }

 //  
 let edge1 = &brep.edges[e1];
 let edge2 = &brep.edges[e2];
 let shared_v = *vi;
 let (e1_other, e2_other) = if edge1.start == shared_v {
 (edge1.end, if edge2.start == shared_v { edge2.end } else { edge2.start })
 } else if edge1.end == shared_v {
 (edge1.start, if edge2.start == shared_v { edge2.end } else { edge2.start })
 } else {
 continue;
 };
 if e1_other == e2_other { continue; } //  

 //  
 let p1 = brep.vertices.get(e1_other).map(|v| v.point);
 let ps = brep.vertices.get(shared_v).map(|v| v.point);
 let p2 = brep.vertices.get(e2_other).map(|v| v.point);
 let (Some(p_a), Some(p_b), Some(p_c)) = (p1, ps, p2) else { continue; };

 //  : AB × BC  
 let ab = p_b - p_a;
 let bc = p_c - p_b;
 let ab_len = ab.length();
 let bc_len = bc.length();
 if ab_len <= 1e-12 || bc_len <= 1e-12 { continue; }
 let cross = ab.cross(bc).length();
 if cross > tol * (ab_len + bc_len) { continue; }
 //  : AB · BC > 0 ( )
 if ab.dot(bc) <= 0.0 { continue; }

 // (face_internal_vertices  )
 let used_as_corner = brep.geom.face_internal_vertices.iter().any(|fiv| fiv.contains(&shared_v));
 if used_as_corner { continue; }

 to_merge.push((e1, e2, shared_v));
 }

 // 3. ( )
 to_merge.sort_unstable();
 to_merge.dedup_by(|a, b| a.1 == b.1 || a.0 == b.1 || a.1 == b.0);

 for &(keep_ei, remove_ei, _shared_v) in &to_merge {
 if keep_ei == remove_ei { continue; }
 // A→B→C
 let ek = &brep.edges[keep_ei];
 let er = &brep.edges[remove_ei];
 let (new_start, new_end) = if ek.start == er.start || ek.start == er.end {
 // keep edge starts at shared vertex, need to extend backward
 let other_v = if ek.start == er.start { er.end } else { er.start };
 (other_v, ek.end)
 } else if ek.end == er.start || ek.end == er.end {
 let other_v = if ek.end == er.start { er.end } else { er.start };
 (ek.start, other_v)
 } else {
 continue;
 };

 // keep  
 if let Some(e) = brep.edges.get_mut(keep_ei) {
 e.start = new_start;
 e.end = new_end;
 }

 // wire  : remove_ei → keep_ei
 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 for fi in 0..brep.solids[si].shells[shi].faces.len() {
 let face = &mut brep.solids[si].shells[shi].faces[fi];
 fn remap_in_wire(wire: &mut Wire, from: usize, to: usize) {
 for we in &mut wire.edges {
 if we.idx == from { we.idx = to; }
 }
 }
 remap_in_wire(&mut face.outer_wire, remove_ei, keep_ei);
 for w in &mut face.inner_wires {
 remap_in_wire(w, remove_ei, keep_ei);
 }
 }
 }
 }

 merged_any = true;
 total_merged += 1;
 }

 if !merged_any { break; }

 //  
 let mut keep_set: Vec<bool> = (0..brep.edges.len()).map(|_| true).collect();
 // wire  ,  to_merge  
 for &(_, remove_ei, _) in &to_merge {
 if remove_ei < keep_set.len() {
 keep_set[remove_ei] = false;
 }
 }
 // edge  
 let mut remap: Vec<Option<usize>> = vec![None; brep.edges.len()];
 let mut new_edges: Vec<Edge> = Vec::new();
 for (i, e) in brep.edges.iter().enumerate() {
 if i < keep_set.len() && keep_set[i] {
 remap[i] = Some(new_edges.len());
 new_edges.push(*e);
 }
 }
 // wire  
 for si in 0..brep.solids.len() {
 for shi in 0..brep.solids[si].shells.len() {
 for fi in 0..brep.solids[si].shells[shi].faces.len() {
 let face = &mut brep.solids[si].shells[shi].faces[fi];
 fn remap_wire(wire: &mut Wire, remap: &[Option<usize>]) {
 for we in &mut wire.edges {
 if let Some(new) = remap.get(we.idx).copied().flatten() {
 we.idx = new;
 }
 }
 }
 remap_wire(&mut face.outer_wire, &remap);
 for w in &mut face.inner_wires {
 remap_wire(w, &remap);
 }
 }
 }
 }
 // geom  
 let mut new_ec: Vec<Option<usize>> = Vec::new();
 for (i, ec) in brep.geom.edge_curve.iter().enumerate() {
 if i < keep_set.len() && keep_set[i] {
 new_ec.push(*ec);
 }
 }
 brep.edges = new_edges;
 brep.geom.edge_curve = new_ec;
 // edge_pcurves compute_face_pcurves  
 // edge_tolerance, edge_degenerated, edge_same_parameter, edge_curve_range  
 if brep.geom.edge_tolerance.len() == remap.len() {
 let mut new_et: Vec<f64> = Vec::new();
 for (i, et) in brep.geom.edge_tolerance.iter().enumerate() {
 if i < keep_set.len() && keep_set[i] { new_et.push(*et); }
 }
 brep.geom.edge_tolerance = new_et;
 }
 if brep.geom.edge_degenerated.len() == remap.len() {
 let mut new_ed: Vec<bool> = Vec::new();
 for (i, ed) in brep.geom.edge_degenerated.iter().enumerate() {
 if i < keep_set.len() && keep_set[i] { new_ed.push(*ed); }
 }
 brep.geom.edge_degenerated = new_ed;
 }
 }

 total_merged
}

#[cfg(test)]
mod tests {
 use super::*;
 use crate::topology::{Face, Shell, Solid, Vertex, Wire};

 #[test]
 fn merge_collinear_edges_in_wires_collapses_fragmented_chain() {
 let mut brep = BRep::new();
 brep.vertices = vec![
 Vertex { point: glam::DVec3::new(0.0, 0.0, 0.0) }, // 0
 Vertex { point: glam::DVec3::new(1.0, 0.0, 0.0) }, // 1 split point
 Vertex { point: glam::DVec3::new(2.0, 0.0, 0.0) }, // 2
 Vertex { point: glam::DVec3::new(2.0, 1.0, 0.0) }, // 3
 Vertex { point: glam::DVec3::new(0.0, 1.0, 0.0) }, // 4
 ];
 brep.edges = vec![
 Edge { start: 0, end: 1 }, // e0
 Edge { start: 1, end: 2 }, // e1
 Edge { start: 2, end: 3 }, // e2
 Edge { start: 3, end: 4 }, // e3
 Edge { start: 4, end: 0 }, // e4
 ];
 brep.solids = vec![Solid {
 shells: vec![Shell {
 faces: vec![Face {
 outer_wire: Wire {
 edges: vec![
 WireEdge::fwd(0),
 WireEdge::fwd(1),
 WireEdge::fwd(2),
 WireEdge::fwd(3),
 WireEdge::fwd(4),
 ],
 },
 inner_wires: vec![],
 normal: glam::DVec3::Z,
 triangles: vec![],
 sample_point: None,
 mesh_dirty: true,
 surface_idx: None,
 }],
 }],
 }];

 let merged = merge_collinear_edges_in_wires(&mut brep, 1e-6);
 assert!(merged >= 1);
 let edge_count = brep.solids[0].shells[0].faces[0].outer_wire.edges.len();
 assert_eq!(edge_count, 4);
 }
}
