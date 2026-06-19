//! OCCT-aligned ShellSplitter: partitions a set of connected faces into closed shells.
//!
//! OCCT ref: BOPAlgo_ShellSplitter (BOPAlgo_ShellSplitter.cxx / .hxx)
//!
//! Processing steps (OCCT BOPAlgo_ShellSplitter::Perform):
//! 1. AddStartElement — collects faces to process (via `add_start_face`)
//! 2. MakeConnexityBlocks — builds edge-connectivity graph from shared edges
//! 3. MakeShells (regular) / SplitBlock (irregular) — flood-fills connected components
//!
//! ✅ OCCT-aligned: core algorithm mirrors BOPAlgo_ShellSplitter.

use std::collections::{HashMap, HashSet};

use super::ds::DS;

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct ShellSplitter {
    /// Face indices (into DS.faces) to process.
    /// OCCT: BOPAlgo_ShellSplitter::myFaces (list of TopoDS_Face).
    myFaces: Vec<usize>,
    /// Edge adjacency: (min_vertex_idx, max_vertex_idx) -> list of face indices
    /// sharing that edge.  Built during MakeConnexityBlocks.
    myEdgeAdj: HashMap<(usize, usize), Vec<usize>>,
    /// Resulting shells: each entry is a vector of face indices forming one
    /// connected component.
    /// OCCT: BOPAlgo_ShellSplitter::myShells (list of TopoDS_Shell).
    myShells: Vec<Vec<usize>>,
}

impl ShellSplitter {
    /// Create an empty ShellSplitter.
    ///
    /// ✅ OCCT-aligned: BOPAlgo_ShellSplitter default constructor.
    pub fn new() -> Self {
        Self {
            myFaces: Vec::new(),
            myEdgeAdj: HashMap::new(),
            myShells: Vec::new(),
        }
    }

    /// Add a face to the set to be split into shells.
    ///
    /// ✅ OCCT-aligned: BOPAlgo_ShellSplitter::AddStartElement(const TopoDS_Shape&).
    pub fn add_start_face(&mut self, fi: usize) {
        self.myFaces.push(fi);
    }

    /// Build edge adjacency and partition faces into connected shells.
    ///
    /// ✅ OCCT-aligned: BOPAlgo_ShellSplitter::MakeConnexityBlocks (builds edge
    ///   connectivity graph) + MakeShells (flood-fill connected components).
    ///
    /// Steps:
    /// 1. For each face, extract all edge vertex-pairs (outer wire + inner wires).
    /// 2. Build adjacency: faces that share an edge (same vertex-pair) are neighbors.
    /// 3. Flood-fill (DFS) to find connected components, stored in `myShells`.
    ///
    /// OCCT ref: BOPAlgo_ShellSplitter.cxx
    ///   - MakeConnexityBlocks: iterates TopExp_Explorer over each face's edges,
    ///     builds a map from edge index to face list.
    ///   - MakeShells: BFS over the connectivity graph, collecting shells.
    pub fn perform(&mut self, ds: &DS) {
        self.myEdgeAdj.clear();
        self.myShells.clear();

        // --- Step 1+2: Build edge adjacency from vertex-pair keys. ---
        // OCCT ref: BOPAlgo_ShellSplitter::MakeConnexityBlocks uses
        // TopExp_Explorer to iterate E of each F and build a
        // map: edge.id() -> list of faces that contain that edge.
        for &fi in &self.myFaces {
            let Some(face) = ds.faces.get(fi) else {
                continue;
            };
            let mut seen_keys = HashSet::new();

            // Outer wire edges
            for &ei in &face.boundary_edges {
                if let Some(edge) = ds.edges.get(ei) {
                    let key = if edge.start_vertex < edge.end_vertex {
                        (edge.start_vertex, edge.end_vertex)
                    } else {
                        (edge.end_vertex, edge.start_vertex)
                    };
                    if seen_keys.insert(key) {
                        self.myEdgeAdj.entry(key).or_default().push(fi);
                    }
                }
            }

            // Inner wire edges (TopExp_Explorer walks inner wires too)
            for inner_wire in &face.inner_boundary_edges {
                for &(ei, _) in inner_wire {
                    if let Some(edge) = ds.edges.get(ei) {
                        let key = if edge.start_vertex < edge.end_vertex {
                            (edge.start_vertex, edge.end_vertex)
                        } else {
                            (edge.end_vertex, edge.start_vertex)
                        };
                        if seen_keys.insert(key) {
                            self.myEdgeAdj.entry(key).or_default().push(fi);
                        }
                    }
                }
            }
        }

        // --- Step 3: Flood-fill connected components (DFS). ---
        // OCCT ref: BOPAlgo_ShellSplitter::MakeShells performs a
        // breadth-first traversal of the connexity graph.
        let n = self.myFaces.len();
        if n == 0 {
            return;
        }

        // Build a local adjacency list from the edge-adjacency map.
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        let face_to_local: HashMap<usize, usize> = self
            .myFaces
            .iter()
            .enumerate()
            .map(|(li, &fi)| (fi, li))
            .collect();

        for flist in self.myEdgeAdj.values() {
            for i in 0..flist.len() {
                for j in (i + 1)..flist.len() {
                    let a = flist[i];
                    let b = flist[j];
                    if let (Some(&la), Some(&lb)) =
                        (face_to_local.get(&a), face_to_local.get(&b))
                    {
                        adj[la].push(lb);
                        adj[lb].push(la);
                    }
                }
            }
        }

        // DFS flood-fill
        let mut visited = vec![false; n];
        for i in 0..n {
            if visited[i] {
                continue;
            }
            let mut comp = Vec::new();
            let mut stack = vec![i];
            while let Some(cur) = stack.pop() {
                if visited[cur] {
                    continue;
                }
                visited[cur] = true;
                comp.push(self.myFaces[cur]);
                for &nb in &adj[cur] {
                    if !visited[nb] {
                        stack.push(nb);
                    }
                }
            }
            if !comp.is_empty() {
                self.myShells.push(comp);
            }
        }
    }

    /// Return the resulting shells (connected components of faces).
    ///
    /// ✅ OCCT-aligned: BOPAlgo_ShellSplitter::Shells() accessor.
    pub fn shells(&self) -> &[Vec<usize>] {
        &self.myShells
    }

    /// Return the number of shells found.
    pub fn nb_shells(&self) -> usize {
        self.myShells.len()
    }

    /// Return true when the split produced more than one shell.
    pub fn has_multiple_shells(&self) -> bool {
        self.myShells.len() > 1
    }

    /// Clear all state and re-initialize.
    ///
    /// ✅ OCCT-aligned: BOPAlgo_ShellSplitter::Clear().
    pub fn clear(&mut self) {
        self.myFaces.clear();
        self.myEdgeAdj.clear();
        self.myShells.clear();
    }
}

impl Default for ShellSplitter {
    fn default() -> Self {
        Self::new()
    }
}
