use std::collections::HashMap;
use crate::{BRep, WireEdge, topology::Wire};

fn oriented_edge_vertices(brep: &BRep, we: WireEdge) -> Option<(usize, usize)> {
 let edges = brep.flat_edges();
 let e = edges.get(we.idx).copied()?;
 if we.forward {
 Some((e.0, e.1))
 } else {
 Some((e.1, e.0))
 }
}

fn points_are_collinear_forward(
 a: glam::DVec3,
 b: glam::DVec3,
 c: glam::DVec3,
 linear_tol: f64,
) -> bool {
 let tol = linear_tol.max(1e-12);
 let ab = b - a;
 let bc = c - b;
 let ab_len = ab.length();
 let bc_len = bc.length();
 if ab_len <= 1e-12 || bc_len <= 1e-12 {
 return false;
 }
 let cross = ab.cross(bc).length();
 let dot = ab.dot(bc);
 cross <= tol * (ab_len + bc_len) && dot > 0.0
}

fn find_existing_edge_between_vertices(brep: &BRep, from: usize, to: usize) -> Option<WireEdge> {
 let edges = brep.flat_edges();
 for (idx, (s, e)) in edges.iter().enumerate() {
  if *s == from && *e == to {
  return Some(WireEdge::fwd(idx));
  }
  if *s == to && *e == from {
  return Some(WireEdge::rev(idx));
  }
 }
 None
}

fn push_topology_edge(brep: &mut BRep, start: usize, end: usize) -> WireEdge {
 let idx = brep.add_edge_flat(start, end, None, [0.0, 1.0]);
 WireEdge::fwd(idx)
}

fn merge_collinear_consecutive_edges_in_wire(
 _brep: &mut BRep,
 wire: &mut Vec<WireEdge>,
 linear_tol: f64,
) -> usize {
 if wire.len() < 4 {
 return 0;
 }
 let tol = linear_tol.max(1e-12);
 let mut merged = 0usize;
 loop {
  let mut i = 0;
  let len = wire.len();
  let mut changed = false;
  while i + 2 < len {
  let u = wire[i].idx;
  let v1 = wire[(i + 1) % len].idx;
  let w = wire[(i + 2) % len].idx;
  if v1 != u && u == w {
   // Collinear check uses vertex positions from BRep
   let p_u = _brep.vertex_point(u).unwrap_or_default();
   let p_v = _brep.vertex_point(v1).unwrap_or_default();
   let p_w = _brep.vertex_point(w).unwrap_or_default();
   if points_are_collinear_forward(p_u, p_v, p_w, tol) {
   // Remove middle edge
   wire.remove((i + 1) % wire.len());
   merged += 1;
   changed = true;
   continue;
   }
  }
  i += 1;
  }
  if !changed {
   break;
  }
 }
 merged
}

/// Merge collinear consecutive edge pairs in all wires of a BRep.
///
/// Scans every face's outer and inner wires for pairs of consecutive edges
/// that share a common vertex and whose 3D geometry is near-collinear.
/// When found, the middle vertex and its two incident edges are replaced by
/// a single direct edge bridging the outer vertices.
///
/// Returns the total number of edge-pair merges performed.
pub fn merge_collinear_edges_in_wires(_brep: &mut BRep, _linear_tol: f64) -> usize {
 // TODO: rewrite for topods::BRep (wire mutation via tshape replacement)
 0
}

/// Merge collinear edges in a BRep by consolidating degree-2 vertices.
///
/// Unlike `merge_collinear_edges_in_wires`, this operates on the full BRep
/// edge graph, not per-wire. It finds vertices with exactly two incident
/// edges that share the same 3D curve and are collinear, then replaces the
/// two edges with a single edge spanning the pair.
///
/// Returns the total number of edge merges performed.
pub fn merge_collinear_brep_edges(_brep: &mut BRep, _linear_tol: f64) -> usize {
 // TODO: rewrite for topods::BRep (edge mutation via tshape replacement)
 0
}
