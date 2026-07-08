//! Topology query helpers — analogous to OCCT `TopExp_Explorer` and
//! `TopExp::MapShapesAndAncestors`.
//!
//! All functions operate on `topods::BRep` and use ShapeRef-based access.

use std::collections::HashSet;

use crate::topology::WireEdge;
use crate::topods;

/// Returns wire edges in order, keeping only the **first** occurrence of each
/// `edge_idx` (semantic edge — one row per edge index in tree UIs).
///
/// Closed periodic faces often list the same seam edge twice (e.g. forward then
/// reverse) so the boundary is closed; the duplicate [`WireEdge`] entries are
/// not separate semantic edges.
pub fn wire_edges_unique_by_index(edges: &[WireEdge]) -> Vec<&WireEdge> {
    let mut seen = HashSet::new();
    edges
        .iter()
        .filter(|we| seen.insert(we.idx))
        .collect()
}

/// Helper: find all vertex tshape indices in a topods::BRep.
fn collect_vertex_indices(brep: &topods::BRep) -> Vec<usize> {
    brep.tshapes
        .iter()
        .enumerate()
        .filter(|(_, ts)| matches!(&**ts, topods::TShape::Vertex(_)))
        .map(|(i, _)| i)
        .collect()
}

/// Helper: find all edge tshape indices in a topods::BRep.
fn collect_edge_indices(brep: &topods::BRep) -> Vec<usize> {
    brep.tshapes
        .iter()
        .enumerate()
        .filter(|(_, ts)| matches!(&**ts, topods::TShape::Edge(_)))
        .map(|(i, _)| i)
        .collect()
}

/// Helper: find all face tshape indices in a topods::BRep.
fn collect_face_indices(brep: &topods::BRep) -> Vec<usize> {
    brep.tshapes
        .iter()
        .enumerate()
        .filter(|(_, ts)| matches!(&**ts, topods::TShape::Face(_)))
        .map(|(i, _)| i)
        .collect()
}

/// Helper: find all solid tshape indices in a topods::BRep.
fn collect_solid_indices(brep: &topods::BRep) -> Vec<usize> {
    brep.tshapes
        .iter()
        .enumerate()
        .filter(|(_, ts)| matches!(&**ts, topods::TShape::Solid(_)))
        .map(|(i, _)| i)
        .collect()
}

/// Returns the flat face tshape indices that reference `edge_tshape_idx`
/// in their outer wire.
///
/// For a well-formed manifold solid each edge is shared by exactly two faces.
/// Returns fewer if the edge is a boundary edge, more if non-manifold.
pub fn edge_adjacent_faces(brep: &topods::BRep, edge_tshape_idx: usize) -> Vec<usize> {
    collect_face_indices(brep)
        .into_iter()
        .filter(|&fi| {
            let topods::TShape::Face(fd) = &*brep.tshapes[fi] else { return false };
            let outer_wire = &brep.tshapes[fd.outer_wire.index];
            let topods::TShape::Wire(wd) = &**outer_wire else { return false };
            wd.edges.iter().any(|sr| sr.index == edge_tshape_idx)
        })
        .collect()
}

/// Returns all edge tshape indices referenced in the outer wire of face
/// `face_tshape_idx`.
///
/// Duplicate edge indices (e.g. seam edges) are preserved as they appear
/// in the wire.
pub fn face_edges(brep: &topods::BRep, face_tshape_idx: usize) -> Vec<usize> {
    let topods::TShape::Face(fd) = &*brep.tshapes[face_tshape_idx] else { return vec![] };
    let topods::TShape::Wire(wd) = &*brep.tshapes[fd.outer_wire.index] else { return vec![] };
    wd.edges.iter().map(|sr| sr.index).collect()
}

/// Returns all edge tshape indices `ei` where the edge has `vertex_tshape_idx`
/// as first or last vertex.
pub fn vertex_adjacent_edges(brep: &topods::BRep, vertex_tshape_idx: usize) -> Vec<usize> {
    collect_edge_indices(brep)
        .into_iter()
        .filter(|&ei| {
            let topods::TShape::Edge(ed) = &*brep.tshapes[ei] else { return false };
            ed.first.index == vertex_tshape_idx || ed.last.index == vertex_tshape_idx
        })
        .collect()
}

/// Number of face TShapes in the BRep.
pub fn face_count(brep: &topods::BRep) -> usize {
    collect_face_indices(brep).len()
}

/// Number of edge TShapes in the BRep.
pub fn edge_count(brep: &topods::BRep) -> usize {
    collect_edge_indices(brep).len()
}

/// Number of vertex TShapes in the BRep.
pub fn vertex_storage_len(brep: &topods::BRep) -> usize {
    vertex_indices(brep).len()
}

/// Number of topological vertices referenced by edges.
pub fn topological_vertex_count(brep: &topods::BRep) -> usize {
    let vi = vertex_indices(brep);
    vi.len()
}

/// Returns `true` if `edge_tshape_idx` is a degenerate edge — i.e. its start
/// and end vertices are the same point, or it is explicitly flagged degenerate.
///
/// Analogous to `BRep_Tool::Degenerated(edge)` in OCCT.
pub fn is_degenerate_edge(brep: &topods::BRep, edge_tshape_idx: usize) -> bool {
    let topods::TShape::Edge(ed) = &*brep.tshapes[edge_tshape_idx] else { return false };
    if ed.degenerated {
        return true;
    }
    let Some(v_start) = brep.tshapes.get(ed.first.index) else { return false };
    let Some(v_end) = brep.tshapes.get(ed.last.index) else { return false };
    let topods::TShape::Vertex(vd_start) = &**v_start else { return false };
    let topods::TShape::Vertex(vd_end) = &**v_end else { return false };
    (vd_end.point - vd_start.point).length_squared() < 1e-20
}

/// Returns all edge tshape indices that are candidates for seam edges — edges
/// whose start and end vertices are the same (or nearly the same) 3D point,
/// but that are *not* flagged as degenerate.
///
/// Seam edges appear on surfaces of revolution (cylinder, sphere, cone, torus)
/// where the UV seam is a real geometric curve even though it starts and ends
/// at the same 3D location.
pub fn periodic_seam_edge_indices(brep: &topods::BRep) -> Vec<usize> {
    collect_edge_indices(brep)
        .into_iter()
        .filter(|&ei| {
            let topods::TShape::Edge(ed) = &*brep.tshapes[ei] else { return false };
            if ed.degenerated {
                return false;
            }
            let Some(v_start) = brep.tshapes.get(ed.first.index) else { return false };
            let Some(v_end) = brep.tshapes.get(ed.last.index) else { return false };
            let topods::TShape::Vertex(vd_start) = &**v_start else { return false };
            let topods::TShape::Vertex(vd_end) = &**v_end else { return false };
            (vd_end.point - vd_start.point).length_squared() < 1e-10
        })
        .collect()
}

/// Returns semantic vertex tshape indices (topological vertices).
///
/// Includes vertices referenced by any topological edge, de-duplicated.
pub fn vertex_indices(brep: &topods::BRep) -> Vec<usize> {
    let mut out = Vec::new();
    for (ei, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Edge(ed) = &**ts {
            out.push(ed.first.index);
            out.push(ed.last.index);
        }
    }
    out.sort_unstable();
    out.dedup();
    let v_count = brep.tshapes.len();
    out.retain(|&vi| vi < v_count);
    out
}

/// Topological vertex indices suitable for creator display/picking overlays.
pub fn salient_vertex_indices(brep: &topods::BRep) -> Vec<usize> {
    const COLLINEAR_CHAIN_DOT_THRESHOLD: f64 = -0.965_925_826;

    vertex_indices(brep)
        .into_iter()
        .filter(|&vi| {
            let adj = vertex_adjacent_edges(brep, vi);
            if adj.len() != 2 {
                return true;
            }
            // Center vertex point
            let topods::TShape::Vertex(vd_center) = &*brep.tshapes[vi] else { return true };
            let center = vd_center.point;

            // Get the two adjacent edges and find the opposite-end vertices
            let topods::TShape::Edge(ed_a) = &*brep.tshapes[adj[0]] else { return true };
            let topods::TShape::Edge(ed_b) = &*brep.tshapes[adj[1]] else { return true };

            let other_a = if ed_a.first.index == vi { ed_a.last.index } else { ed_a.first.index };
            let other_b = if ed_b.first.index == vi { ed_b.last.index } else { ed_b.first.index };

            let topods::TShape::Vertex(vd_a) = &*brep.tshapes[other_a] else { return true };
            let topods::TShape::Vertex(vd_b) = &*brep.tshapes[other_b] else { return true };

            let dir_a = (vd_a.point - center).normalize();
            let dir_b = (vd_b.point - center).normalize();
            dir_a.dot(dir_b) > COLLINEAR_CHAIN_DOT_THRESHOLD
        })
        .collect()
}
