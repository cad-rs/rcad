//! Graph-topology wrapper for `BRep`, analogous to OCCT's `BRepGraph` module.
//!
//! `BRepGraph` pre-computes and caches all adjacency relations so that repeated
//! topology queries run in O(1) instead of the O(n) scans performed by the
//! free functions in [`crate::topo_query`].
//!
//! # What is cached
//!
//! | Query | OCCT equivalent |
//! |---|---|
//! | edge → adjacent face indices | `TopExp::MapShapesAndAncestors` |
//! | face → edge indices | `TopExp::MapShapes` |
//! | vertex → adjacent edge indices | `TopExp::MapShapesAndAncestors` |
//! | vertex → adjacent face indices | derived from above |
//! | edge → (start_vertex, end_vertex) | `BRep_Tool::Vertices` |
//!
//! # Mutation tracking
//!
//! Each vertex, edge, and face has a corresponding dirty bit.  Algorithms that
//! mutate a BRep entity can set these bits via the `mark_*_modified` methods,
//! and downstream steps can query them to decide what needs re-evaluation.
//! This is a lightweight analogue to OCCT's `BRepCheck_Analyzer` invalidation
//! mechanism and the history-event bus.
//!
//! # Traversal
//!
//! `BRepGraph` provides DFS and BFS iterators over face/edge/vertex connectivity:
//! - [`BRepGraph::dfs_faces`] — visits all faces reachable (via shared edges) from a seed face.
//! - [`BRepGraph::bfs_faces`] — same graph, breadth-first order.
//! - [`BRepGraph::dfs_edges_from_vertex`] — walks edge-adjacency from a seed vertex DFS.
//!
//! # Examples
//!
//! ```rust
//! use rcad_kernel::{BRep, BRepGraph, PrimitiveSolid};
//!
//! let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
//! let graph = BRepGraph::from_brep(&brep);
//!
//! // O(1) adjacency
//! let adj_faces = graph.edge_adjacent_faces(0);
//! assert_eq!(adj_faces.len(), 2);
//!
//! // Traversal: all 6 faces of a box are connected
//! let visited = graph.bfs_faces(0).collect::<Vec<_>>();
//! assert_eq!(visited.len(), 6);
//!
//! // Manifold / closed inspection
//! assert!(graph.is_manifold());
//! assert!(graph.is_closed());
//! ```

use crate::BRep;
use std::collections::{HashSet, VecDeque};

// ─────────────────────────────────────────────────────────────────────────────
// BRepGraph
// ─────────────────────────────────────────────────────────────────────────────

/// Cached graph-topology view of a [`BRep`].
///
/// Build with [`BRepGraph::from_brep`].  The graph is a **snapshot** — it does
/// not update automatically when the source `BRep` is mutated.  Call
/// `from_brep` again to refresh.
#[derive(Debug, Clone)]
pub struct BRepGraph {
    // ── Counts ────────────────────────────────────────────────────────────────
    pub vertex_count: usize,
    pub edge_count: usize,
    pub face_count: usize,

    // ── Adjacency tables (indexed by the corresponding entity index) ──────────
    /// `edge_to_faces[edge_idx]` = list of flat face indices that include this edge
    /// in their outer wire.
    edge_to_faces: Vec<Vec<usize>>,

    /// `face_to_edges[face_idx]` = list of edge indices in the outer wire of that face.
    face_to_edges: Vec<Vec<usize>>,

    /// `vertex_to_edges[vertex_idx]` = list of edge indices where start == vertex_idx
    /// or end == vertex_idx.
    vertex_to_edges: Vec<Vec<usize>>,

    /// `vertex_to_faces[vertex_idx]` = list of flat face indices whose outer wire
    /// references a wire-edge that touches this vertex.
    vertex_to_faces: Vec<Vec<usize>>,

    /// `edge_endpoints[edge_idx]` = (start_vertex_idx, end_vertex_idx).
    edge_endpoints: Vec<(usize, usize)>,

    // ── Dirty / modified bits ─────────────────────────────────────────────────
    /// Dirty flags for vertices (by vertex index).
    vertex_dirty: Vec<bool>,
    /// Dirty flags for edges (by edge index).
    edge_dirty: Vec<bool>,
    /// Dirty flags for faces (by flat face index).
    face_dirty: Vec<bool>,
}

/// Non-manifold topology summary derived from edge-face adjacency.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct NonManifoldSummary {
    /// Edge indices with exactly one adjacent face.
    pub boundary_edges: Vec<usize>,
    /// Edge indices with zero adjacent faces.
    pub orphan_edges: Vec<usize>,
    /// Edge indices with more than two adjacent faces.
    pub multi_face_edges: Vec<usize>,
    /// Vertex indices touched by at least one multi-face edge.
    pub non_manifold_vertices: Vec<usize>,
}

impl NonManifoldSummary {
    pub fn is_clean(&self) -> bool {
        self.boundary_edges.is_empty() && self.orphan_edges.is_empty() && self.multi_face_edges.is_empty()
    }
}

/// Categorised repair hints for a non-manifold BRep.
///
/// Each hint describes one problem and suggests a concrete remedy action.
/// Analogous to OCCT `BRepCheck_Analyzer` diagnostic entries.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RepairHint {
    /// Two boundary edges share both endpoints — they can be stitched together.
    StitchablePair {
        edge_a: usize,
        edge_b: usize,
        face_a: usize,
        face_b: usize,
    },
    /// A boundary edge has no stitch candidate — the hole must be capped.
    UnmatchedBoundaryEdge { edge_idx: usize, face_idx: usize },
    /// An orphan edge is attached to no face — it should be removed.
    OrphanEdge { edge_idx: usize },
    /// An edge is shared by more than two faces — it must be split into
    /// separate copies, one per face pair.
    MultiManifoldEdge { edge_idx: usize, face_count: usize },
    /// A vertex lies on a multi-face edge — it may require duplication when
    /// the surrounding edge is split.
    NonManifoldVertex { vertex_idx: usize, connected_multi_edges: Vec<usize> },
}

/// Detailed actionable repair hints derived from a `NonManifoldSummary`.
///
/// Build with [`BRepGraph::repair_hints`].
#[derive(Debug, Clone, Default)]
pub struct ManifoldRepairHints {
    pub hints: Vec<RepairHint>,
}

impl ManifoldRepairHints {
    /// Returns `true` when there are no repair items.
    pub fn is_empty(&self) -> bool {
        self.hints.is_empty()
    }

    /// All hints of the `StitchablePair` variant.
    pub fn stitchable_pairs(&self) -> impl Iterator<Item = &RepairHint> {
        self.hints.iter().filter(|h| matches!(h, RepairHint::StitchablePair { .. }))
    }

    /// All hints of the `OrphanEdge` variant.
    pub fn orphan_edges(&self) -> impl Iterator<Item = &RepairHint> {
        self.hints.iter().filter(|h| matches!(h, RepairHint::OrphanEdge { .. }))
    }
}

impl BRepGraph {
    // ── Construction ──────────────────────────────────────────────────────────

    /// Build a `BRepGraph` from the current state of `brep`.
    ///
    /// Iterates all topology entities once (O(V + E + F + wire_edges)) to
    /// populate the adjacency tables.  Subsequent adjacency queries are O(1).
    pub fn from_brep(brep: &BRep) -> Self {
        let vc = brep.vertices.len();
        let ec = brep.edges.len();

        // ── Pre-fill edge endpoints ───────────────────────────────────────────
        let edge_endpoints: Vec<(usize, usize)> = brep
            .edges
            .iter()
            .map(|e| (e.start, e.end))
            .collect();

        // ── vertex → edges ────────────────────────────────────────────────────
        let mut vertex_to_edges: Vec<Vec<usize>> = vec![Vec::new(); vc];
        for (ei, ep) in edge_endpoints.iter().enumerate() {
            if ep.0 < vc {
                vertex_to_edges[ep.0].push(ei);
            }
            if ep.1 < vc && ep.1 != ep.0 {
                vertex_to_edges[ep.1].push(ei);
            }
        }

        // ── Build flat face list and face/edge/vertex adjacency ───────────────
        // Count faces first.
        let fc: usize = brep
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();

        let mut edge_to_faces: Vec<Vec<usize>> = vec![Vec::new(); ec];
        let mut face_to_edges: Vec<Vec<usize>> = vec![Vec::new(); fc];
        let mut vertex_to_faces: Vec<Vec<usize>> = vec![Vec::new(); vc];

        let mut flat_fi = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let mut face_edge_set = Vec::new();
                    let mut face_vert_set = HashSet::new();
                    for we in &face.outer_wire.edges {
                        if we.idx < ec {
                            // edge → face
                            if !edge_to_faces[we.idx].contains(&flat_fi) {
                                edge_to_faces[we.idx].push(flat_fi);
                            }
                            // face → edge
                            if !face_edge_set.contains(&we.idx) {
                                face_edge_set.push(we.idx);
                            }
                            // vertices touched by this wire edge
                            let (vs, ve) = edge_endpoints[we.idx];
                            if vs < vc {
                                face_vert_set.insert(vs);
                            }
                            if ve < vc {
                                face_vert_set.insert(ve);
                            }
                        }
                    }
                    // Also walk inner wires.
                    for inner in &face.inner_wires {
                        for we in &inner.edges {
                            if we.idx < ec {
                                if !edge_to_faces[we.idx].contains(&flat_fi) {
                                    edge_to_faces[we.idx].push(flat_fi);
                                }
                                if !face_edge_set.contains(&we.idx) {
                                    face_edge_set.push(we.idx);
                                }
                                let (vs, ve) = edge_endpoints[we.idx];
                                if vs < vc {
                                    face_vert_set.insert(vs);
                                }
                                if ve < vc {
                                    face_vert_set.insert(ve);
                                }
                            }
                        }
                    }
                    face_to_edges[flat_fi] = face_edge_set;
                    // vertex → faces
                    for vi in face_vert_set {
                        if !vertex_to_faces[vi].contains(&flat_fi) {
                            vertex_to_faces[vi].push(flat_fi);
                        }
                    }
                    flat_fi += 1;
                }
            }
        }

        BRepGraph {
            vertex_count: vc,
            edge_count: ec,
            face_count: fc,
            edge_to_faces,
            face_to_edges,
            vertex_to_edges,
            vertex_to_faces,
            edge_endpoints,
            vertex_dirty: vec![false; vc],
            edge_dirty: vec![false; ec],
            face_dirty: vec![false; fc],
        }
    }

    // ── O(1) adjacency queries ────────────────────────────────────────────────

    /// Faces that share `edge_idx` in their wire.  Expected length 2 for a
    /// manifold edge, 1 for a boundary edge, 0 if the index is out of range.
    #[inline]
    pub fn edge_adjacent_faces(&self, edge_idx: usize) -> &[usize] {
        self.edge_to_faces
            .get(edge_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Edge indices in the outer (+inner) wires of `face_idx`.
    #[inline]
    pub fn face_edges(&self, face_idx: usize) -> &[usize] {
        self.face_to_edges
            .get(face_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Edge indices incident to `vertex_idx`.
    #[inline]
    pub fn vertex_adjacent_edges(&self, vertex_idx: usize) -> &[usize] {
        self.vertex_to_edges
            .get(vertex_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Face indices adjacent to `vertex_idx` (faces whose wire touches it).
    #[inline]
    pub fn vertex_adjacent_faces(&self, vertex_idx: usize) -> &[usize] {
        self.vertex_to_faces
            .get(vertex_idx)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// (start_vertex, end_vertex) for `edge_idx`.
    #[inline]
    pub fn edge_endpoints(&self, edge_idx: usize) -> Option<(usize, usize)> {
        self.edge_endpoints.get(edge_idx).copied()
    }

    // ── Manifold / closed inspection ──────────────────────────────────────────

    /// Returns `true` if every edge is shared by exactly 2 faces (no boundary
    /// or non-manifold edges).
    ///
    /// Analogous to `BRepCheck_Shell::Closed()` manifold test in OCCT.
    pub fn is_manifold(&self) -> bool {
        self.edge_to_faces.iter().all(|adj| adj.len() == 2)
    }

    /// Returns `true` if no edge has fewer than 2 adjacent faces (i.e., no
    /// boundary / free edges exist).
    pub fn is_closed(&self) -> bool {
        self.edge_to_faces.iter().all(|adj| adj.len() >= 2)
    }

    /// Returns edge indices where the number of adjacent faces is not exactly 2.
    ///
    /// - `adj.len() == 1` → free (boundary) edge
    /// - `adj.len() > 2`  → non-manifold edge
    /// - `adj.len() == 0` → orphan edge (not used in any face wire)
    pub fn non_manifold_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() != 2)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Edge indices with exactly one adjacent face.
    pub fn boundary_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() == 1)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Edge indices with zero adjacent faces.
    pub fn orphan_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.is_empty())
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Edge indices with more than two adjacent faces.
    pub fn multi_face_edges(&self) -> Vec<usize> {
        self.edge_to_faces
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() > 2)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// Vertex indices touched by multi-face edges (>2 adjacent faces).
    pub fn non_manifold_vertices(&self) -> Vec<usize> {
        let mut verts = HashSet::new();
        for ei in self.multi_face_edges() {
            if let Some((vs, ve)) = self.edge_endpoints(ei) {
                verts.insert(vs);
                verts.insert(ve);
            }
        }
        let mut out: Vec<usize> = verts.into_iter().collect();
        out.sort_unstable();
        out
    }

    /// Combined non-manifold summary report.
    pub fn non_manifold_summary(&self) -> NonManifoldSummary {
        NonManifoldSummary {
            boundary_edges: self.boundary_edges(),
            orphan_edges: self.orphan_edges(),
            multi_face_edges: self.multi_face_edges(),
            non_manifold_vertices: self.non_manifold_vertices(),
        }
    }

    /// The number of faces sharing `edge_idx` (the edge "valence").
    ///
    /// - 0 → orphan edge
    /// - 1 → boundary / free edge
    /// - 2 → manifold edge (expected)
    /// - >2 → non-manifold / T-junction edge
    #[inline]
    pub fn edge_valence(&self, edge_idx: usize) -> usize {
        self.edge_to_faces
            .get(edge_idx)
            .map(|v| v.len())
            .unwrap_or(0)
    }

    /// The number of edges incident to `vertex_idx` (vertex degree).
    #[inline]
    pub fn vertex_degree(&self, vertex_idx: usize) -> usize {
        self.vertex_to_edges
            .get(vertex_idx)
            .map(|v| v.len())
            .unwrap_or(0)
    }

    /// Generate actionable [`ManifoldRepairHints`] for this graph.
    ///
    /// The algorithm classifies each problem edge and groups compatible
    /// boundary edges as stitchable pairs (same orientation·length within 1e-6).
    ///
    /// Analogous to OCCT `BRepCheck_Analyzer` hint emission.
    pub fn repair_hints(&self, brep: &crate::BRep) -> ManifoldRepairHints {
        let mut hints = Vec::new();

        // ── Multi-face (non-manifold) edges ───────────────────────────────────
        let multi = self.multi_face_edges();
        for ei in &multi {
            hints.push(RepairHint::MultiManifoldEdge {
                edge_idx: *ei,
                face_count: self.edge_valence(*ei),
            });
        }

        // ── Non-manifold vertices ─────────────────────────────────────────────
        for vi in self.non_manifold_vertices() {
            let connected: Vec<usize> = self
                .vertex_adjacent_edges(vi)
                .iter()
                .copied()
                .filter(|&ei| self.edge_valence(ei) > 2)
                .collect();
            if !connected.is_empty() {
                hints.push(RepairHint::NonManifoldVertex {
                    vertex_idx: vi,
                    connected_multi_edges: connected,
                });
            }
        }

        // ── Orphan edges ──────────────────────────────────────────────────────
        for ei in self.orphan_edges() {
            hints.push(RepairHint::OrphanEdge { edge_idx: ei });
        }

        // ── Boundary edges: find stitchable pairs ─────────────────────────────
        let boundary = self.boundary_edges();
        let mut paired: HashSet<usize> = HashSet::new();

        for i in 0..boundary.len() {
            let ei = boundary[i];
            if paired.contains(&ei) {
                continue;
            }
            let face_i = self.edge_adjacent_faces(ei).first().copied().unwrap_or(0);
            let (va, vb) = match self.edge_endpoints(ei) {
                Some(ep) => ep,
                None => continue,
            };
            let pa = brep.vertices.get(va).map(|v| v.point).unwrap_or_default();
            let pb = brep.vertices.get(vb).map(|v| v.point).unwrap_or_default();
            let len_i = (pb - pa).length();

            let mut matched = false;
            for j in (i + 1)..boundary.len() {
                let ej = boundary[j];
                if paired.contains(&ej) {
                    continue;
                }
                let (vc2, vd) = match self.edge_endpoints(ej) {
                    Some(ep) => ep,
                    None => continue,
                };
                let pc = brep.vertices.get(vc2).map(|v| v.point).unwrap_or_default();
                let pd = brep.vertices.get(vd).map(|v| v.point).unwrap_or_default();
                let len_j = (pd - pc).length();

                // Stitchable if both endpoints are coincident (within 1e-6),
                // possibly with reversed orientation.
                let direct  = (pa - pc).length() < 1e-6 && (pb - pd).length() < 1e-6;
                let reverse = (pa - pd).length() < 1e-6 && (pb - pc).length() < 1e-6;
                let same_len = (len_i - len_j).abs() < 1e-6;

                if same_len && (direct || reverse) {
                    let face_j = self.edge_adjacent_faces(ej).first().copied().unwrap_or(0);
                    hints.push(RepairHint::StitchablePair {
                        edge_a: ei,
                        edge_b: ej,
                        face_a: face_i,
                        face_b: face_j,
                    });
                    paired.insert(ei);
                    paired.insert(ej);
                    matched = true;
                    break;
                }
            }
            if !matched {
                hints.push(RepairHint::UnmatchedBoundaryEdge {
                    edge_idx: ei,
                    face_idx: face_i,
                });
            }
        }

        ManifoldRepairHints { hints }
    }

    /// Returns vertex indices that are referenced by a number of edges other
    /// than the expected valence.  For a closed manifold each vertex has
    /// valence ≥ 2 (connected by at least 2 edges).
    pub fn low_valence_vertices(&self) -> Vec<usize> {
        self.vertex_to_edges
            .iter()
            .enumerate()
            .filter(|(_, adj)| adj.len() < 2)
            .map(|(vi, _)| vi)
            .collect()
    }

    // ── Dirty / modification tracking ─────────────────────────────────────────

    /// Mark vertex `vertex_idx` as having been modified.
    pub fn mark_vertex_modified(&mut self, vertex_idx: usize) {
        if let Some(d) = self.vertex_dirty.get_mut(vertex_idx) {
            *d = true;
        }
    }

    /// Mark edge `edge_idx` as having been modified.
    pub fn mark_edge_modified(&mut self, edge_idx: usize) {
        if let Some(d) = self.edge_dirty.get_mut(edge_idx) {
            *d = true;
        }
    }

    /// Mark face `face_idx` (flat index) as having been modified.
    pub fn mark_face_modified(&mut self, face_idx: usize) {
        if let Some(d) = self.face_dirty.get_mut(face_idx) {
            *d = true;
        }
    }

    /// Clear all dirty flags.
    pub fn clear_dirty(&mut self) {
        self.vertex_dirty.iter_mut().for_each(|d| *d = false);
        self.edge_dirty.iter_mut().for_each(|d| *d = false);
        self.face_dirty.iter_mut().for_each(|d| *d = false);
    }

    /// All vertex indices currently marked dirty.
    pub fn modified_vertices(&self) -> Vec<usize> {
        self.vertex_dirty
            .iter()
            .enumerate()
            .filter(|(_, d)| **d)
            .map(|(vi, _)| vi)
            .collect()
    }

    /// All edge indices currently marked dirty.
    pub fn modified_edges(&self) -> Vec<usize> {
        self.edge_dirty
            .iter()
            .enumerate()
            .filter(|(_, d)| **d)
            .map(|(ei, _)| ei)
            .collect()
    }

    /// All face indices (flat) currently marked dirty.
    pub fn modified_faces(&self) -> Vec<usize> {
        self.face_dirty
            .iter()
            .enumerate()
            .filter(|(_, d)| **d)
            .map(|(fi, _)| fi)
            .collect()
    }

    /// Returns `true` if any entity is dirty.
    pub fn has_modifications(&self) -> bool {
        self.vertex_dirty.iter().any(|&d| d)
            || self.edge_dirty.iter().any(|&d| d)
            || self.face_dirty.iter().any(|&d| d)
    }

    // ── Graph traversal ───────────────────────────────────────────────────────

    /// Depth-first traversal of faces, starting from `seed_face_idx`.
    ///
    /// Faces are considered adjacent when they share at least one edge.
    /// Returns a `DfsFaces` iterator that yields flat face indices in DFS
    /// discovery order.
    ///
    /// Analogous to `BRepTools_WireExplorer` / `TopExp_Explorer` over `TopAbs_FACE`.
    pub fn dfs_faces(&self, seed_face_idx: usize) -> DfsFaces<'_> {
        DfsFaces::new(self, seed_face_idx)
    }

    /// Breadth-first traversal of faces, starting from `seed_face_idx`.
    ///
    /// Returns a `BfsFaces` iterator that yields flat face indices in BFS
    /// discovery order.
    pub fn bfs_faces(&self, seed_face_idx: usize) -> BfsFaces<'_> {
        BfsFaces::new(self, seed_face_idx)
    }

    /// Depth-first traversal of edges, starting from `seed_vertex_idx`.
    ///
    /// Visits each edge reachable by following edge-adjacency from the seed
    /// vertex.  Yields edge indices in DFS discovery order.
    pub fn dfs_edges_from_vertex(&self, seed_vertex_idx: usize) -> DfsEdgesFromVertex<'_> {
        DfsEdgesFromVertex::new(self, seed_vertex_idx)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Iterators
// ─────────────────────────────────────────────────────────────────────────────

/// DFS iterator over faces reachable from a seed face via shared edges.
pub struct DfsFaces<'g> {
    graph: &'g BRepGraph,
    stack: Vec<usize>,
    visited: HashSet<usize>,
}

impl<'g> DfsFaces<'g> {
    fn new(graph: &'g BRepGraph, seed: usize) -> Self {
        let mut visited = HashSet::new();
        let mut stack = Vec::new();
        if seed < graph.face_count {
            stack.push(seed);
            visited.insert(seed);
        }
        DfsFaces { graph, stack, visited }
    }
}

impl Iterator for DfsFaces<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let fi = self.stack.pop()?;
        // Push unvisited neighbors: faces reachable via edges of `fi`.
        for &ei in self.graph.face_edges(fi) {
            for &adj_fi in self.graph.edge_adjacent_faces(ei) {
                if self.visited.insert(adj_fi) {
                    self.stack.push(adj_fi);
                }
            }
        }
        Some(fi)
    }
}

/// BFS iterator over faces reachable from a seed face via shared edges.
pub struct BfsFaces<'g> {
    graph: &'g BRepGraph,
    queue: VecDeque<usize>,
    visited: HashSet<usize>,
}

impl<'g> BfsFaces<'g> {
    fn new(graph: &'g BRepGraph, seed: usize) -> Self {
        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();
        if seed < graph.face_count {
            queue.push_back(seed);
            visited.insert(seed);
        }
        BfsFaces { graph, queue, visited }
    }
}

impl Iterator for BfsFaces<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let fi = self.queue.pop_front()?;
        for &ei in self.graph.face_edges(fi) {
            for &adj_fi in self.graph.edge_adjacent_faces(ei) {
                if self.visited.insert(adj_fi) {
                    self.queue.push_back(adj_fi);
                }
            }
        }
        Some(fi)
    }
}

/// DFS iterator over edges reachable from a seed vertex via edge-vertex adjacency.
pub struct DfsEdgesFromVertex<'g> {
    graph: &'g BRepGraph,
    edge_stack: Vec<usize>,
    vertex_stack: Vec<usize>,
    visited_edges: HashSet<usize>,
    visited_verts: HashSet<usize>,
}

impl<'g> DfsEdgesFromVertex<'g> {
    fn new(graph: &'g BRepGraph, seed_vertex: usize) -> Self {
        let mut visited_verts = HashSet::new();
        let mut vertex_stack = Vec::new();
        if seed_vertex < graph.vertex_count {
            vertex_stack.push(seed_vertex);
            visited_verts.insert(seed_vertex);
        }
        DfsEdgesFromVertex {
            graph,
            edge_stack: Vec::new(),
            vertex_stack,
            visited_edges: HashSet::new(),
            visited_verts,
        }
    }
}

impl Iterator for DfsEdgesFromVertex<'_> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        // Expand vertices on the vertex stack into unvisited edges.
        while let Some(vi) = self.vertex_stack.pop() {
            for &ei in self.graph.vertex_adjacent_edges(vi) {
                if self.visited_edges.insert(ei) {
                    self.edge_stack.push(ei);
                }
            }
        }
        let ei = self.edge_stack.pop()?;
        // Push far-end vertex.
        if let Some((vs, ve)) = self.graph.edge_endpoints(ei) {
            for &vn in &[vs, ve] {
                if self.visited_verts.insert(vn) {
                    self.vertex_stack.push(vn);
                }
            }
        }
        // Drain new vertices immediately.
        while let Some(vi) = self.vertex_stack.pop() {
            for &nei in self.graph.vertex_adjacent_edges(vi) {
                if self.visited_edges.insert(nei) {
                    self.edge_stack.push(nei);
                }
            }
        }
        Some(ei)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Checkpoint / Mutation Guard
// ─────────────────────────────────────────────────────────────────────────────

/// A point-in-time snapshot of all adjacency tables and dirty flags in a
/// [`BRepGraph`].
///
/// Created by [`BRepGraph::checkpoint`]; restored by
/// [`BRepGraph::restore_from_checkpoint`] or automatically on drop of an
/// uncommitted [`BRepGraphMutGuard`].
///
/// The snapshot owns cloned copies of all internal adjacency data, so it is
/// independent of the live `BRepGraph`.  For large models this can be
/// memory-heavy; prefer scoped guards for short mutations rather than
/// storing long-lived checkpoints.
///
/// Analogous to `BRepGraph_Compact` / rollback semantics in OCCT 8.0.
#[derive(Debug, Clone)]
pub struct BRepGraphCheckpoint {
    vertex_count: usize,
    edge_count: usize,
    face_count: usize,
    edge_to_faces: Vec<Vec<usize>>,
    face_to_edges: Vec<Vec<usize>>,
    vertex_to_edges: Vec<Vec<usize>>,
    vertex_to_faces: Vec<Vec<usize>>,
    edge_endpoints: Vec<(usize, usize)>,
    vertex_dirty: Vec<bool>,
    edge_dirty: Vec<bool>,
    face_dirty: Vec<bool>,
}

impl BRepGraph {
    // ── Checkpoint / restore ──────────────────────────────────────────────────

    /// Snapshot the current graph state into a [`BRepGraphCheckpoint`].
    ///
    /// The checkpoint can later be passed to [`BRepGraph::restore_from_checkpoint`]
    /// to undo all mutations applied to the graph since the checkpoint was taken.
    pub fn checkpoint(&self) -> BRepGraphCheckpoint {
        BRepGraphCheckpoint {
            vertex_count: self.vertex_count,
            edge_count: self.edge_count,
            face_count: self.face_count,
            edge_to_faces: self.edge_to_faces.clone(),
            face_to_edges: self.face_to_edges.clone(),
            vertex_to_edges: self.vertex_to_edges.clone(),
            vertex_to_faces: self.vertex_to_faces.clone(),
            edge_endpoints: self.edge_endpoints.clone(),
            vertex_dirty: self.vertex_dirty.clone(),
            edge_dirty: self.edge_dirty.clone(),
            face_dirty: self.face_dirty.clone(),
        }
    }

    /// Restore the graph to the state captured by `cp`.
    ///
    /// After this call the graph's adjacency tables, counts, and dirty flags
    /// are exactly as they were when `cp` was created.  Any mutations applied
    /// after `cp` was taken are undone.
    pub fn restore_from_checkpoint(&mut self, cp: &BRepGraphCheckpoint) {
        self.vertex_count = cp.vertex_count;
        self.edge_count = cp.edge_count;
        self.face_count = cp.face_count;
        self.edge_to_faces = cp.edge_to_faces.clone();
        self.face_to_edges = cp.face_to_edges.clone();
        self.vertex_to_edges = cp.vertex_to_edges.clone();
        self.vertex_to_faces = cp.vertex_to_faces.clone();
        self.edge_endpoints = cp.edge_endpoints.clone();
        self.vertex_dirty = cp.vertex_dirty.clone();
        self.edge_dirty = cp.edge_dirty.clone();
        self.face_dirty = cp.face_dirty.clone();
    }

    /// Open a scoped mutation guard over this graph.
    ///
    /// While the guard is alive, the caller is free to mutate the `BRepGraph`
    /// (dirty-marking, adjacency updates, etc.) via [`BRepGraphMutGuard::graph`].
    /// On drop the guard automatically rolls back all changes unless
    /// [`BRepGraphMutGuard::commit`] or [`BRepGraphMutGuard::commit_unchecked`]
    /// was called first.
    ///
    /// Analogous to `BRepGraph_MutGuard` in OCCT 8.0.
    pub fn begin_mutation(&mut self) -> BRepGraphMutGuard<'_> {
        BRepGraphMutGuard::new(self)
    }

    // ── Invariant validation ──────────────────────────────────────────────────

    /// Validate internal topology invariants and return a list of error strings.
    ///
    /// An empty list means the graph is consistent.  Checks include:
    /// - Table lengths match the declared counts.
    /// - Every (start, end) vertex index in `edge_endpoints` is in bounds.
    /// - Every face index referenced by `edge_to_faces` is in bounds.
    /// - Every edge index referenced by `face_to_edges` is in bounds.
    ///
    /// Analogous to `BRepGraph_Validate` in OCCT 8.0.
    pub fn validate_invariants(&self) -> Vec<String> {
        let mut errors: Vec<String> = Vec::new();

        if self.edge_to_faces.len() != self.edge_count {
            errors.push(format!(
                "edge_to_faces length {} ≠ edge_count {}",
                self.edge_to_faces.len(),
                self.edge_count
            ));
        }
        if self.face_to_edges.len() != self.face_count {
            errors.push(format!(
                "face_to_edges length {} ≠ face_count {}",
                self.face_to_edges.len(),
                self.face_count
            ));
        }
        if self.vertex_to_edges.len() != self.vertex_count {
            errors.push(format!(
                "vertex_to_edges length {} ≠ vertex_count {}",
                self.vertex_to_edges.len(),
                self.vertex_count
            ));
        }
        if self.vertex_to_faces.len() != self.vertex_count {
            errors.push(format!(
                "vertex_to_faces length {} ≠ vertex_count {}",
                self.vertex_to_faces.len(),
                self.vertex_count
            ));
        }
        if self.edge_endpoints.len() != self.edge_count {
            errors.push(format!(
                "edge_endpoints length {} ≠ edge_count {}",
                self.edge_endpoints.len(),
                self.edge_count
            ));
        }
        if self.vertex_dirty.len() != self.vertex_count {
            errors.push(format!(
                "vertex_dirty length {} ≠ vertex_count {}",
                self.vertex_dirty.len(),
                self.vertex_count
            ));
        }
        if self.edge_dirty.len() != self.edge_count {
            errors.push(format!(
                "edge_dirty length {} ≠ edge_count {}",
                self.edge_dirty.len(),
                self.edge_count
            ));
        }
        if self.face_dirty.len() != self.face_count {
            errors.push(format!(
                "face_dirty length {} ≠ face_count {}",
                self.face_dirty.len(),
                self.face_count
            ));
        }

        // Check edge_endpoints vertex indices are in bounds.
        for (ei, &(vs, ve)) in self.edge_endpoints.iter().enumerate() {
            if vs >= self.vertex_count {
                errors.push(format!(
                    "edge {} start vertex {} out of range (vertex_count={})",
                    ei, vs, self.vertex_count
                ));
            }
            if ve >= self.vertex_count {
                errors.push(format!(
                    "edge {} end vertex {} out of range (vertex_count={})",
                    ei, ve, self.vertex_count
                ));
            }
        }

        // Check face indices referenced by edge_to_faces are in bounds.
        for (ei, faces) in self.edge_to_faces.iter().enumerate() {
            for &fi in faces {
                if fi >= self.face_count {
                    errors.push(format!(
                        "edge {} references face {} out of range (face_count={})",
                        ei, fi, self.face_count
                    ));
                }
            }
        }

        // Check edge indices referenced by face_to_edges are in bounds.
        for (fi, edges) in self.face_to_edges.iter().enumerate() {
            for &ei in edges {
                if ei >= self.edge_count {
                    errors.push(format!(
                        "face {} references edge {} out of range (edge_count={})",
                        fi, ei, self.edge_count
                    ));
                }
            }
        }

        errors
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRepGraphMutGuard
// ─────────────────────────────────────────────────────────────────────────────

/// RAII-scoped mutation guard for a [`BRepGraph`].
///
/// On construction a checkpoint is taken of the entire graph state.  The
/// caller mutates the graph through [`BRepGraphMutGuard::graph`].  When the
/// guard goes out of scope:
///
/// - If [`BRepGraphMutGuard::commit`] or
///   [`BRepGraphMutGuard::commit_unchecked`] was called first, the new state
///   is kept.
/// - Otherwise the graph is **rolled back** to the state at guard creation.
///
/// # Example
///
/// ```rust
/// use rcad_kernel::{BRep, BRepGraph};
/// use rcad_kernel::geom::PrimitiveSolid;
///
/// let brep = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let mut graph = BRepGraph::from_brep(&brep);
///
/// {
///     let mut guard = graph.begin_mutation();
///     guard.graph().mark_face_modified(0);
///     guard.graph().mark_edge_modified(1);
///     // Validate before committing.
///     guard.commit().expect("invariants should hold");
/// }
/// assert!(graph.modified_faces().contains(&0));
/// ```
pub struct BRepGraphMutGuard<'g> {
    graph: &'g mut BRepGraph,
    checkpoint: BRepGraphCheckpoint,
    committed: bool,
}

impl<'g> BRepGraphMutGuard<'g> {
    fn new(graph: &'g mut BRepGraph) -> Self {
        let checkpoint = graph.checkpoint();
        Self { graph, checkpoint, committed: false }
    }

    /// Access the wrapped graph for mutation.
    pub fn graph(&mut self) -> &mut BRepGraph {
        self.graph
    }

    /// Read-only view of the graph (useful for inspect-before-commit checks).
    pub fn graph_ref(&self) -> &BRepGraph {
        self.graph
    }

    /// Validate topology invariants and commit if they all pass.
    ///
    /// Returns `Ok(())` on success; `Err(errors)` listing every violated
    /// invariant if validation fails — in which case the graph is **not**
    /// rolled back, allowing callers to inspect the invalid intermediate state
    /// before deciding to call [`BRepGraphMutGuard::rollback`] explicitly.
    pub fn commit(mut self) -> Result<(), Vec<String>> {
        let errors = self.graph.validate_invariants();
        if errors.is_empty() {
            self.committed = true;
            Ok(())
        } else {
            Err(errors)
        }
    }

    /// Commit the current state without running invariant validation.
    ///
    /// Use when you have guaranteed the invariants externally, or when
    /// performance is critical and you have already validated.
    pub fn commit_unchecked(mut self) {
        self.committed = true;
    }

    /// Explicitly roll back all mutations made since the guard was created.
    ///
    /// This is equivalent to simply dropping the guard without committing,
    /// but makes the intent explicit.
    pub fn rollback(self) {
        // `committed` remains false; Drop handles the actual restore.
    }

    /// Consume the guard and return the snapshot taken at guard creation,
    /// keeping the current (mutated) graph state.
    ///
    /// This is a lower-level escape hatch for callers that want to commit
    /// unconditionally but retain the checkpoint for a manual restore later.
    pub fn into_checkpoint(mut self) -> BRepGraphCheckpoint {
        self.committed = true;
        self.checkpoint.clone()
    }
}

impl Drop for BRepGraphMutGuard<'_> {
    fn drop(&mut self) {
        if !self.committed {
            self.graph.restore_from_checkpoint(&self.checkpoint);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        BRep, Edge, Face, Shell, Solid, Vertex, Wire, WireEdge,
        geom::PrimitiveSolid,
    };
    use glam::DVec3;

    fn unit_box() -> BRep {
        BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        })
    }

    /// Build a minimal non-manifold BRep where edge 0 is shared by 3 faces.
    fn non_manifold_tripod() -> BRep {
        let vertices = vec![
            Vertex { point: DVec3::new(0.0, 0.0, 0.0) }, // 0
            Vertex { point: DVec3::new(1.0, 0.0, 0.0) }, // 1
            Vertex { point: DVec3::new(0.0, 1.0, 0.0) }, // 2
            Vertex { point: DVec3::new(0.0, 0.0, 1.0) }, // 3
            Vertex { point: DVec3::new(0.0, -1.0, 0.0) }, // 4
        ];
        // Edge 0 is the shared spine. Other edges are unique per triangle face.
        let edges = vec![
            Edge { start: 0, end: 1 }, // shared by 3 faces
            Edge { start: 1, end: 2 },
            Edge { start: 2, end: 0 },
            Edge { start: 1, end: 3 },
            Edge { start: 3, end: 0 },
            Edge { start: 1, end: 4 },
            Edge { start: 4, end: 0 },
        ];

        let f0 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(1, true),
                    WireEdge::new(2, true),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Z,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f1 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(3, true),
                    WireEdge::new(4, true),
                ],
            },
            inner_wires: vec![],
            normal: DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };
        let f2 = Face {
            outer_wire: Wire {
                edges: vec![
                    WireEdge::new(0, true),
                    WireEdge::new(5, true),
                    WireEdge::new(6, true),
                ],
            },
            inner_wires: vec![],
            normal: -DVec3::Y,
            triangles: vec![],
            mesh_dirty: true,
        };

        BRep {
            vertices,
            edges,
            solids: vec![Solid {
                shells: vec![Shell {
                    faces: vec![f0, f1, f2],
                }],
            }],
            geom: Default::default(),
        }
    }

    // ── Construction & counts ─────────────────────────────────────────────────

    #[test]
    fn box_entity_counts_correct() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        assert_eq!(g.vertex_count, 8, "box: 8 vertices");
        assert_eq!(g.edge_count, 12, "box: 12 edges");
        assert_eq!(g.face_count, 6, "box: 6 faces");
    }

    #[test]
    fn empty_brep_yields_zero_counts() {
        let brep = BRep::new();
        let g = BRepGraph::from_brep(&brep);
        assert_eq!(g.vertex_count, 0);
        assert_eq!(g.edge_count, 0);
        assert_eq!(g.face_count, 0);
        assert!(g.edge_adjacent_faces(0).is_empty());
        assert!(!g.has_modifications());
    }

    // ── O(1) adjacency ────────────────────────────────────────────────────────

    #[test]
    fn box_every_edge_has_two_adjacent_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for ei in 0..g.edge_count {
            let adj = g.edge_adjacent_faces(ei);
            assert_eq!(adj.len(), 2, "edge {ei}: expected 2 adjacent faces, got {adj:?}");
        }
    }

    #[test]
    fn box_every_face_has_four_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for fi in 0..g.face_count {
            let edges = g.face_edges(fi);
            assert_eq!(edges.len(), 4, "face {fi}: expected 4 edges, got {edges:?}");
        }
    }

    #[test]
    fn box_every_vertex_has_three_adjacent_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for vi in 0..g.vertex_count {
            let adj = g.vertex_adjacent_edges(vi);
            assert_eq!(adj.len(), 3, "vertex {vi}: expected 3 adjacent edges, got {adj:?}");
        }
    }

    #[test]
    fn box_every_vertex_has_three_adjacent_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for vi in 0..g.vertex_count {
            let adj = g.vertex_adjacent_faces(vi);
            assert_eq!(adj.len(), 3, "vertex {vi}: expected 3 adjacent faces, got {adj:?}");
        }
    }

    #[test]
    fn edge_endpoints_match_brep_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for (ei, edge) in brep.edges.iter().enumerate() {
            let ep = g.edge_endpoints(ei).expect("endpoints should exist");
            assert_eq!(ep.0, edge.start);
            assert_eq!(ep.1, edge.end);
        }
    }

    // ── Manifold / closed inspection ──────────────────────────────────────────

    #[test]
    fn box_is_manifold_and_closed() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        assert!(g.is_manifold(), "unit box should be manifold");
        assert!(g.is_closed(), "unit box should be closed");
        assert!(g.non_manifold_edges().is_empty(), "no non-manifold edges");
        assert!(g.low_valence_vertices().is_empty(), "no low-valence vertices");

        let summary = g.non_manifold_summary();
        assert!(summary.is_clean());
        assert!(summary.boundary_edges.is_empty());
        assert!(summary.orphan_edges.is_empty());
        assert!(summary.multi_face_edges.is_empty());
        assert!(summary.non_manifold_vertices.is_empty());
    }

    #[test]
    fn detects_multi_face_non_manifold_edge() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);

        assert!(!g.is_manifold(), "tripod should be non-manifold");
        assert!(!g.is_closed(), "tripod has boundary edges on each side face");

        let multi = g.multi_face_edges();
        assert_eq!(multi, vec![0], "edge 0 should be the only multi-face edge");
        assert_eq!(g.boundary_edges(), vec![1, 2, 3, 4, 5, 6]);
        assert!(g.orphan_edges().is_empty(), "no orphan edges expected");

        let verts = g.non_manifold_vertices();
        assert_eq!(verts, vec![0, 1], "shared edge endpoints should be non-manifold vertices");

        let summary = g.non_manifold_summary();
        assert!(!summary.is_clean());
        assert_eq!(summary.multi_face_edges, vec![0]);
        assert_eq!(summary.non_manifold_vertices, vec![0, 1]);
    }

    // ── Dirty / modification tracking ─────────────────────────────────────────

    #[test]
    fn dirty_marking_and_query() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);
        assert!(!g.has_modifications());

        g.mark_vertex_modified(0);
        g.mark_edge_modified(3);
        g.mark_face_modified(5);

        assert!(g.has_modifications());
        assert_eq!(g.modified_vertices(), vec![0]);
        assert_eq!(g.modified_edges(), vec![3]);
        assert_eq!(g.modified_faces(), vec![5]);

        g.clear_dirty();
        assert!(!g.has_modifications());
        assert!(g.modified_vertices().is_empty());
        assert!(g.modified_edges().is_empty());
        assert!(g.modified_faces().is_empty());
    }

    #[test]
    fn out_of_range_mark_is_safe() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);
        // Should not panic.
        g.mark_vertex_modified(9999);
        g.mark_edge_modified(9999);
        g.mark_face_modified(9999);
        assert!(!g.has_modifications());
    }

    // ── DFS / BFS face traversal ──────────────────────────────────────────────

    #[test]
    fn bfs_faces_visits_all_box_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let mut visited: Vec<usize> = g.bfs_faces(0).collect();
        visited.sort_unstable();
        assert_eq!(visited, vec![0, 1, 2, 3, 4, 5], "BFS should visit all 6 faces");
    }

    #[test]
    fn dfs_faces_visits_all_box_faces() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let mut visited: Vec<usize> = g.dfs_faces(0).collect();
        visited.sort_unstable();
        assert_eq!(visited, vec![0, 1, 2, 3, 4, 5], "DFS should visit all 6 faces");
    }

    #[test]
    fn dfs_bfs_face_traversal_no_duplicates() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let bfs: Vec<usize> = g.bfs_faces(2).collect();
        let dfs: Vec<usize> = g.dfs_faces(2).collect();
        // No duplicates.
        let bfs_set: HashSet<usize> = bfs.iter().copied().collect();
        let dfs_set: HashSet<usize> = dfs.iter().copied().collect();
        assert_eq!(bfs.len(), bfs_set.len(), "BFS duplicates");
        assert_eq!(dfs.len(), dfs_set.len(), "DFS duplicates");
    }

    #[test]
    fn out_of_range_seed_returns_empty_traversal() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        assert!(g.bfs_faces(999).next().is_none());
        assert!(g.dfs_faces(999).next().is_none());
    }

    // ── DFS edge traversal ────────────────────────────────────────────────────

    #[test]
    fn dfs_edges_from_vertex_visits_all_box_edges() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let mut visited: Vec<usize> = g.dfs_edges_from_vertex(0).collect();
        visited.sort_unstable();
        assert_eq!(visited.len(), 12, "DFS from any vertex should reach all 12 edges");
        assert_eq!(visited, (0..12).collect::<Vec<_>>());
    }

    // ── Edge valence and vertex degree ────────────────────────────────────────

    #[test]
    fn edge_valence_manifold_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        for ei in 0..brep.edges.len() {
            assert_eq!(g.edge_valence(ei), 2, "every box edge should have valence 2");
        }
        assert_eq!(g.edge_valence(9999), 0, "out-of-range returns 0");
    }

    #[test]
    fn edge_valence_tripod_multi_face_edge() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        assert_eq!(g.edge_valence(0), 3, "shared spine should have valence 3");
        assert_eq!(g.edge_valence(1), 1, "side edge should have valence 1 (boundary)");
    }

    #[test]
    fn vertex_degree_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        // Every vertex of a box is incident to exactly 3 edges.
        for vi in 0..brep.vertices.len() {
            assert_eq!(g.vertex_degree(vi), 3, "box vertex {vi} should have degree 3");
        }
    }

    // ── repair_hints ──────────────────────────────────────────────────────────

    #[test]
    fn repair_hints_empty_for_manifold_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        assert!(
            hints.is_empty(),
            "a closed manifold box should have no repair hints"
        );
    }

    #[test]
    fn repair_hints_detect_multi_manifold_edge_in_tripod() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        // Must include at least one MultiManifoldEdge hint for edge 0.
        let has_multi = hints.hints.iter().any(|h| {
            matches!(h, RepairHint::MultiManifoldEdge { edge_idx: 0, .. })
        });
        assert!(has_multi, "tripod must generate a MultiManifoldEdge hint for edge 0");
    }

    #[test]
    fn repair_hints_detect_non_manifold_vertex_in_tripod() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        let has_nm_vert = hints.hints.iter().any(|h| {
            matches!(h, RepairHint::NonManifoldVertex { .. })
        });
        assert!(has_nm_vert, "tripod must generate NonManifoldVertex hints");
    }

    #[test]
    fn repair_hints_detect_unmatched_boundary_in_tripod() {
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let hints = g.repair_hints(&brep);
        // Tripod side edges (1-6) are all boundary edges;
        // they cannot be stitched to each other because there are no matching
        // partner edges with coincident endpoints.
        let has_unmatched = hints.hints.iter().any(|h| {
            matches!(h, RepairHint::UnmatchedBoundaryEdge { .. })
        });
        assert!(has_unmatched, "tripod must have unmatched boundary edges");
    }

    // ── Checkpoint / restore ──────────────────────────────────────────────────

    #[test]
    fn checkpoint_captures_clean_state() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);
        let cp = g.checkpoint();
        // Checkpoint counts should match the live graph.
        assert_eq!(cp.vertex_count, g.vertex_count);
        assert_eq!(cp.edge_count, g.edge_count);
        assert_eq!(cp.face_count, g.face_count);
        // Dirty flags are all false in a freshly built graph.
        assert!(cp.vertex_dirty.iter().all(|&d| !d));
        assert!(cp.edge_dirty.iter().all(|&d| !d));
        assert!(cp.face_dirty.iter().all(|&d| !d));
        // Mark something dirty after the checkpoint.
        g.mark_face_modified(0);
        // Restoring should clear the flag.
        g.restore_from_checkpoint(&cp);
        assert!(!g.modified_faces().contains(&0), "restore should clear face 0 dirty bit");
    }

    #[test]
    fn restore_from_checkpoint_undoes_dirty_marks() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);
        // Nothing dirty initially.
        assert!(!g.has_modifications());

        let cp = g.checkpoint();

        // Mark several entities dirty.
        g.mark_vertex_modified(2);
        g.mark_edge_modified(7);
        g.mark_face_modified(4);
        assert!(g.has_modifications());

        // Restore.
        g.restore_from_checkpoint(&cp);
        assert!(!g.has_modifications(), "all dirty bits should be cleared after restore");
        assert!(g.modified_vertices().is_empty());
        assert!(g.modified_edges().is_empty());
        assert!(g.modified_faces().is_empty());
    }

    // ── BRepGraphMutGuard ─────────────────────────────────────────────────────

    #[test]
    fn mut_guard_commit_keeps_changes() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);

        {
            let mut guard = g.begin_mutation();
            guard.graph().mark_face_modified(1);
            guard.graph().mark_edge_modified(3);
            guard.commit().expect("box invariants should hold");
        }

        // Changes persisted after commit.
        assert!(g.modified_faces().contains(&1), "face 1 should be dirty after commit");
        assert!(g.modified_edges().contains(&3), "edge 3 should be dirty after commit");
    }

    #[test]
    fn mut_guard_drop_without_commit_rolls_back() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);

        {
            let mut guard = g.begin_mutation();
            guard.graph().mark_face_modified(0);
            guard.graph().mark_vertex_modified(5);
            // Drop without commit → automatic rollback.
        }

        assert!(!g.has_modifications(), "uncommitted guard drop should roll back dirty bits");
    }

    #[test]
    fn mut_guard_explicit_rollback_reverts_changes() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);

        {
            let mut guard = g.begin_mutation();
            guard.graph().mark_face_modified(2);
            guard.rollback(); // explicit rollback
        }

        assert!(!g.modified_faces().contains(&2), "explicit rollback should clear face 2");
    }

    #[test]
    fn mut_guard_commit_unchecked_keeps_changes() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);

        {
            let mut guard = g.begin_mutation();
            guard.graph().mark_edge_modified(0);
            guard.commit_unchecked();
        }

        assert!(g.modified_edges().contains(&0), "commit_unchecked should keep edge 0 dirty");
    }

    #[test]
    fn mut_guard_into_checkpoint_keeps_state_and_returns_snapshot() {
        let brep = unit_box();
        let mut g = BRepGraph::from_brep(&brep);

        let cp = {
            let mut guard = g.begin_mutation();
            guard.graph().mark_face_modified(3);
            guard.into_checkpoint()
        };

        // State should be kept (committed via into_checkpoint).
        assert!(g.modified_faces().contains(&3));
        // The returned checkpoint should represent the pre-mutation state (face 3 not dirty).
        assert!(!cp.face_dirty.get(3).copied().unwrap_or(true),
            "checkpoint should capture pre-mutation state (face 3 was clean)");
    }

    // ── validate_invariants ───────────────────────────────────────────────────

    #[test]
    fn validate_invariants_passes_for_unit_box() {
        let brep = unit_box();
        let g = BRepGraph::from_brep(&brep);
        let errors = g.validate_invariants();
        assert!(errors.is_empty(), "unit box graph should have no invariant violations: {errors:?}");
    }

    #[test]
    fn validate_invariants_passes_for_non_manifold_tripod() {
        // The tripod is non-manifold but structurally consistent (table lengths match,
        // indices are in bounds).
        let brep = non_manifold_tripod();
        let g = BRepGraph::from_brep(&brep);
        let errors = g.validate_invariants();
        assert!(errors.is_empty(), "tripod graph should have no index invariant violations: {errors:?}");
    }
}
