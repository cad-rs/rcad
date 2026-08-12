//! Open shell sewing: merge multiple BReps into a single watertight solid
//! by identifying and stitching near-coincident edges.
//!
//! Analogous to OCCT `BRepOffsetAPI_Sewing`.
//!
//! # Algorithm
//! 1. Concatenate all vertices, edges, and faces from every input BRep into
//!    a single pool, reindexing everything.
//! 2. **Vertex merging** (union-find): vertices within `tolerance` of each
//!    other are merged **only when they originate from different input BReps**.
//!    Same-shell vertices are never merged even if duplicated numerically — this
//!    avoids collapsing quad corners when stitching independent shells (`RCAD ZP3`).
//! 3. **Edge matching**: after vertex merging, edges that share both endpoint
//!    vertices (in either orientation) and originated from different input
//!    BReps are considered "stitched" — they represent the same boundary edge.
//!    Duplicate edges are removed and wire references updated; when the merged
//!    duplicate uses the opposite `(start,end)` order from the canonical edge,
//!    [`WireEdge.forward`] is toggled so loops remain consistently oriented.
//! 4. **Shell assembly**: all faces are collected into a single shell.
//!    Free edges (with only one incident face) are reported.
//!
//! # Limitations
//! - Surfaces are not preserved from input faces (flat-index intermediate model).
//! - Only the outer wire of each face is considered during edge matching
//!   (inner wires are preserved as-is, with reindexed edge refs).

use glam::DVec3;
use rcad_kernel::{
    topods,
    topods::{Orientation, TShape},
    topology::{Edge, Face, Vertex, Wire, WireEdge},
};

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// Result returned by [`sew_shells`].
#[derive(Debug)]
pub struct SewingResult {
    /// The merged BRep containing all input faces in a single shell.
    pub brep: topods::BRep,
    /// Number of edge pairs that were stitched (shared boundary resolved).
    pub stitched_pairs: usize,
    /// Indices of edges in the result BRep that have only one incident face
    /// (open boundary edges).
    pub free_edges: Vec<usize>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Merge multiple BReps into a single BRep, stitching near-coincident edges.
///
/// # Arguments
/// * `breps` — slice of BReps to merge (must have at least one solid).
/// * `tolerance` — vertex proximity threshold for merging.
///
/// # Examples
/// ```rust
/// use rcad_kernel::topods;
/// use rcad_modeling::sew_shells;
///
/// let (a, _) = topods::BRep::build_unit_cube();
/// let (b, _) = topods::BRep::build_unit_cube();
/// let result = sew_shells(&[a, b], 1e-6);
/// assert!(result.brep.has_solids());
/// ```
pub fn sew_shells(breps: &[topods::BRep], tolerance: f64) -> SewingResult {
    if breps.is_empty() {
        return SewingResult {
            brep: topods::BRep::new(),
            stitched_pairs: 0,
            free_edges: Vec::new(),
        };
    }

    // ── Step 1: concatenate everything ────────────────────────────────────

    let mut all_vertices: Vec<Vertex> = Vec::new();
    let mut all_edges: Vec<Edge> = Vec::new();
    let mut all_faces: Vec<Face> = Vec::new();

    // Track the vertex/edge offset for each input BRep
    let mut vertex_offsets: Vec<usize> = Vec::with_capacity(breps.len());
    let mut edge_offsets: Vec<usize> = Vec::with_capacity(breps.len());
    // Which input shell (`breps` index) each concatenated vertex came from.
    let mut vertex_src_brep: Vec<usize> = Vec::new();

    // Reusable mapping buffers for each input brep
    let mut ts_to_flat_v: Vec<usize> = Vec::new();
    let mut ts_to_flat_e: Vec<usize> = Vec::new();

    for (bi, brep) in breps.iter().enumerate() {
        let v_off = all_vertices.len();
        let e_off = all_edges.len();
        vertex_offsets.push(v_off);
        edge_offsets.push(e_off);

        // Build tshape-index → flat-vertex-index map for this input BRep.
        // Flat vertex indices include the running offset so edges can use
        // the value directly.
        ts_to_flat_v.clear();
        ts_to_flat_v.resize(brep.tshapes.len(), usize::MAX);
        {
            let mut vi = 0usize;
            for (tsi, ts) in brep.tshapes.iter().enumerate() {
                if let TShape::Vertex(_) = &**ts {
                    ts_to_flat_v[tsi] = v_off + vi;
                    vi += 1;
                }
            }
        }

        // Build tshape-index → flat-edge-index map (also offset).
        ts_to_flat_e.clear();
        ts_to_flat_e.resize(brep.tshapes.len(), usize::MAX);
        {
            let mut ei = 0usize;
            for (tsi, ts) in brep.tshapes.iter().enumerate() {
                if let TShape::Edge(_) = &**ts {
                    ts_to_flat_e[tsi] = e_off + ei;
                    ei += 1;
                }
            }
        }

        // Vertices
        for ts in &brep.tshapes {
            if let TShape::Vertex(vd) = &**ts {
                all_vertices.push(Vertex { point: vd.point });
                vertex_src_brep.push(bi);
            }
        }

        // Edges (reindexed)
        for ts in &brep.tshapes {
            if let TShape::Edge(ed) = &**ts {
                all_edges.push(Edge {
                    start: ts_to_flat_v[ed.first.index],
                    end: ts_to_flat_v[ed.last.index],
                });
            }
        }

        // Faces (reindex edge refs in wires)
        let solid_ref = brep.tshapes.iter().enumerate().find_map(|(i, ts)| {
            if matches!(&**ts, TShape::Solid(_)) {
                Some(topods::Shape::synthetic(i, topods::Orientation::Forward))
            } else {
                None
            }
        });
        if let Some(sr) = solid_ref {
            let sd = brep.solid(sr);
            for shell_sr in &sd.shells {
                let shd = brep.shell(shell_sr.clone());
                for face_sr in &shd.faces {
                    let fd = brep.face(face_sr.clone());

                    let reindex_wire = |w_sr: topods::Shape| -> Wire {
                        let wd = brep.wire(w_sr);
                        Wire {
                            edges: wd
                                .edges
                                .iter()
                                .map(|we| WireEdge {
                                    idx: ts_to_flat_e[we.index],
                                    forward: we.orientation == Orientation::Forward,
                                    location: we.location,
                                })
                                .collect(),
                        }
                    };

                    all_faces.push(Face {
                        outer_wire: reindex_wire(fd.outer_wire.clone()),
                        inner_wires: fd.inner_wires.iter().map(|w| reindex_wire(w.clone())).collect(),
                        normal: DVec3::ZERO,
                        triangles: Vec::new(),
                        sample_point: fd.sample_point,
                        mesh_dirty: true,
                        surface_idx: None,
                    });
                }
            }
        }
    }

    let n_verts = all_vertices.len();
    let n_edges = all_edges.len();

    // ── Step 2: union-find vertex merge ───────────────────────────────────

    let mut parent: Vec<usize> = (0..n_verts).collect();

    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[rb] = ra;
        }
    }

    let tol2 = tolerance * tolerance;
    for i in 0..n_verts {
        for j in (i + 1)..n_verts {
            if vertex_src_brep[i] == vertex_src_brep[j] {
                continue;
            }
            let d2 = (all_vertices[i].point - all_vertices[j].point).length_squared();
            if d2 <= tol2 {
                union(&mut parent, i, j);
            }
        }
    }

    // Build canonical index map: old index → canonical representative
    let canon: Vec<usize> = (0..n_verts).map(|i| find(&mut parent, i)).collect();

    // Remap edge endpoints
    for e in &mut all_edges {
        e.start = canon[e.start];
        e.end = canon[e.end];
    }

    // ── Step 3: edge deduplication ─────────────────────────────────────────
    // Build a map from (min_v, max_v) → first edge index.
    // For each duplicate, record the pair (kept, duplicate) as stitched.

    use std::collections::HashMap;
    let mut edge_key_to_idx: HashMap<(usize, usize), usize> = HashMap::new();
    // Maps old edge index → canonical edge index (for dedup)
    let mut edge_canon: Vec<usize> = (0..n_edges).collect();
    let mut stitched_pairs = 0usize;

    for i in 0..n_edges {
        let e = &all_edges[i];
        let key = (e.start.min(e.end), e.start.max(e.end));
        if let Some(&existing) = edge_key_to_idx.get(&key) {
            // This edge is a duplicate of `existing` — stitch
            edge_canon[i] = existing;
            stitched_pairs += 1;
        } else {
            edge_key_to_idx.insert(key, i);
        }
    }

    // When a duplicate edge is merged onto the canonical edge, orientations may differ:
    // `all_edges[i]` and `all_edges[edge_canon[i]]` share endpoints but might be reversed.
    // Wires must flip `WireEdge.forward` when remapping idx so face boundaries stay closed.
    let mut dup_edge_flip: Vec<bool> = vec![false; n_edges];
    for i in 0..n_edges {
        let c = edge_canon[i];
        if c == i {
            continue;
        }
        let ei = &all_edges[i];
        let ee = &all_edges[c];
        dup_edge_flip[i] = ei.start != ee.start;
    }

    // Remap wire edge indices in all faces
    for face in &mut all_faces {
        let remap = |w: &mut Wire| {
            for we in &mut w.edges {
                let old_idx = we.idx;
                let c = edge_canon[old_idx];
                we.idx = c;
                if old_idx != c && dup_edge_flip[old_idx] {
                    we.forward = !we.forward;
                }
            }
        };
        remap(&mut face.outer_wire);
        for iw in &mut face.inner_wires {
            remap(iw);
        }
    }

    // ── Step 4: build result BRep using topods builder ─────────────────────

    // Compact edges: keep only canonical ones
    let kept_edges: Vec<Edge> = (0..n_edges)
        .filter(|&i| edge_canon[i] == i)
        .map(|i| all_edges[i])
        .collect();

    // Remap edge indices: old canonical → compacted index
    let mut compact_edge_idx: Vec<usize> = vec![0; n_edges];
    {
        let mut ci = 0;
        for i in 0..n_edges {
            if edge_canon[i] == i {
                compact_edge_idx[i] = ci;
                ci += 1;
            }
        }
    }

    // Remap wires to compacted edge indices
    for face in &mut all_faces {
        let remap2 = |w: &mut Wire| {
            for we in &mut w.edges {
                we.idx = compact_edge_idx[edge_canon[we.idx]];
            }
        };
        remap2(&mut face.outer_wire);
        for iw in &mut face.inner_wires {
            remap2(iw);
        }
    }

    // Compact vertices (keep canonical representatives, renumber)
    let mut seen_verts: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &c in &canon {
        seen_verts.insert(c);
    }
    let mut canon_to_compact: Vec<usize> = vec![0; n_verts];
    let mut compact_vertices: Vec<Vertex> = Vec::new();
    for i in 0..n_verts {
        if canon[i] == i {
            canon_to_compact[i] = compact_vertices.len();
            compact_vertices.push(all_vertices[i]);
        }
    }

    // Remap edge endpoints to compacted vertices
    let mut compact_edges = kept_edges;
    for e in &mut compact_edges {
        e.start = canon_to_compact[e.start];
        e.end = canon_to_compact[e.end];
    }

    // Remap face triangles to compacted vertices (triangles are empty in the
    // new pipeline but the field exists on topology::Face).
    for face in &mut all_faces {
        for tri in &mut face.triangles {
            tri[0] = canon_to_compact[canon[tri[0]]];
            tri[1] = canon_to_compact[canon[tri[1]]];
            tri[2] = canon_to_compact[canon[tri[2]]];
        }
    }

    // ── Step 5: find free edges ───────────────────────────────────────────

    let n_compact_edges = compact_edges.len();
    let mut edge_face_count: Vec<usize> = vec![0; n_compact_edges];
    for face in &all_faces {
        for we in &face.outer_wire.edges {
            if we.idx < n_compact_edges {
                edge_face_count[we.idx] += 1;
            }
        }
        for iw in &face.inner_wires {
            for we in &iw.edges {
                if we.idx < n_compact_edges {
                    edge_face_count[we.idx] += 1;
                }
            }
        }
    }
    let free_edges: Vec<usize> = (0..n_compact_edges)
        .filter(|&i| edge_face_count[i] == 1)
        .collect();

    // ── Step 6: build topods::BRep result ─────────────────────────────────

    let mut result = topods::BRep::new();

    // Add all vertices, build flat-index → Shape map
    let mut vert_refs: Vec<topods::Shape> = Vec::with_capacity(compact_vertices.len());
    for v in &compact_vertices {
        vert_refs.push(result.add_tvertex(v.point));
    }

    // Add all edges
    let mut edge_refs: Vec<topods::Shape> = Vec::with_capacity(compact_edges.len());
    for e in &compact_edges {
        let first = vert_refs[e.start].clone();
        let last = vert_refs[e.end].clone();
        edge_refs.push(result.add_tedge(None, first, last, [0.0, 1.0]));
    }

    // Add all faces (build wires, then faces)
    let mut face_refs = Vec::with_capacity(all_faces.len());
    for f in &all_faces {
        let mut build_wire = |w: &Wire| -> topods::Shape {
            let edge_srs: Vec<topods::Shape> = w
                .edges
                .iter()
                .map(|we| {
                    let orient = if we.forward {
                        Orientation::Forward
                    } else {
                        Orientation::Reversed
                    };
                    topods::Shape::synthetic(edge_refs[we.idx].index, orient)
                })
                .collect();
            result.add_twire(edge_srs)
        };

        let outer_wire_sr = build_wire(&f.outer_wire);
        let inner_wire_refs: Vec<topods::Shape> = f.inner_wires.iter().map(build_wire).collect();

        let face_sr = result.add_tface(
            None,
            outer_wire_sr,
            inner_wire_refs,
            f.sample_point,
            None,
            vec![],
            false,
        );
        face_refs.push(face_sr);
    }

    // Add shell and solid
    let shell_sr = result.add_tshell(face_refs);
    result.add_tsolid(vec![shell_sr]);

    SewingResult {
        brep: result,
        stitched_pairs,
        free_edges,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sew_single_brep_is_identity() {
        let (brep, _root) = topods::BRep::build_unit_cube();
        let result = sew_shells(std::slice::from_ref(&brep), 1e-6);
        assert_eq!(result.stitched_pairs, 0);
        assert_eq!(
            result.free_edges.len(),
            0,
            "closed box should have no free edges, got {:?}",
            result.free_edges
        );
    }

    #[test]
    fn sew_two_boxes_identifies_shared_face() {
        let (a, _) = topods::BRep::build_unit_cube();
        let (b, _) = topods::BRep::build_unit_cube();
        let result = sew_shells(&[a, b], 1e-6);
        assert!(
            result.stitched_pairs > 0,
            "expected stitched pairs for coincident boxes, got 0"
        );
        println!(
            "sew two coincident boxes: stitched={}, free={:?}",
            result.stitched_pairs, result.free_edges
        );
    }

    #[test]
    fn sew_empty_input() {
        let result = sew_shells(&[], 1e-6);
        assert_eq!(result.stitched_pairs, 0);
        assert!(!result.brep.has_solids());
    }
}
