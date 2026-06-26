//! OCCT-aligned ShellSplitter: partitions a set of connected faces into closed shells.
//!
//! OCCT ref: BOPAlgo_ShellSplitter (BOPAlgo_ShellSplitter.cxx / .hxx)
//!
//! ✅ OCCT-aligned: full algorithm mirrors BOPAlgo_ShellSplitter.
//!   - Perform → MakeConnexityBlocks + MakeShells (L137-149)
//!   - MakeShells: regular blocks → MakeShell; irregular → SplitBlock via
//!     parallel CBK (L621-679)
//!   - SplitBlock: free-edge removal + flood-fill + RefineShell (L153-421)
//!   - RefineShell: multi-connected edge splitting (L443-617)
//!
//! Architecture note: OCCT uses TopoDS_Shape; rcad uses DS face/edge indices.

use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;

use super::ds::{ConnexityBlock, DS};
use crate::boptools;

// ============================================================================
// BOPAlgo_CBK — OCCT L48-88: parallel SplitBlock wrapper
// ============================================================================

/// OCCT-aligned: BOPAlgo_CBK (L48-88).
/// Wraps a ConnexityBlock for parallel SplitBlock execution.
#[derive(Debug, Clone)]
#[allow(non_snake_case)]
struct SplitBlockCbk {
    /// OCCT: BOPAlgo_CBK::myPCB — pointer to the connexity block.
    myFaces: Vec<usize>,
    /// Loops produced by SplitBlock (stored here before collecting).
    myLoops: Vec<Vec<usize>>,
}

impl SplitBlockCbk {
    fn new() -> Self {
        Self {
            myFaces: Vec::new(),
            myLoops: Vec::new(),
        }
    }

    /// OCCT: SetConnexityBlock.
    fn set_faces(&mut self, faces: Vec<usize>) {
        self.myFaces = faces;
    }

    /// OCCT: ConnexityBlock accessor.
    fn loops(&self) -> &[Vec<usize>] {
        &self.myLoops
    }

    /// OCCT: Perform — calls SplitBlock on the stored block.
    fn perform(&mut self, ds: &DS) {
        self.myLoops.clear();
        split_block(&self.myFaces, ds, &mut self.myLoops);
    }
}

/// OCCT-aligned: BOPAlgo_VectorOfCBK — dynamic array of CBK.
type VectorOfCbk = Vec<SplitBlockCbk>;

/// OCCT-aligned: GetEdgeOff (BOPTools_AlgoTools.cxx L1099-1130).
/// Given edge `theE1` and face `theF2`, find whether the edge appears
/// in the face and return the same edge index.  Orientation matching
/// is deferred to GetFaceOff; rcad returns the raw edge index.
fn get_edge_off(ei: usize, fi: usize, ds: &DS) -> Option<usize> {
    let face = &ds.faces[fi];
    if face.boundary_edges.contains(&ei) {
        return Some(ei);
    }
    for wire in &face.inner_boundary_edges {
        for &(e, _) in wire {
            if e == ei {
                return Some(ei);
            }
        }
    }
    None
}

/// Build edge-face adjacency map for a set of face indices.
/// OCCT: TopExp::MapShapesAndAncestors(shell/face, EDGE, FACE, map).
fn build_ef_map(faces: &[usize], ds: &DS) -> BTreeMap<usize, Vec<usize>> {
    let mut ef_map: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &fi in faces {
        if let Some(face) = ds.faces.get(fi) {
            for &ei in &face.boundary_edges {
                ef_map.entry(ei).or_default().push(fi);
            }
            for wire in &face.inner_boundary_edges {
                for &(ei, _) in wire {
                    ef_map.entry(ei).or_default().push(fi);
                }
            }
        }
    }
    ef_map
}

/// Check whether a shell (set of faces) is closed.
/// OCCT: BRep_Tool::IsClosed(shell).  A shell is closed when every edge
/// is referenced exactly twice (once from each adjacent face).
fn is_shell_closed(shell_faces: &[usize], ds: &DS) -> bool {
    let ef = build_ef_map(shell_faces, ds);
    ef.values().all(|flist| flist.len() == 2)
}

/// OCCT-aligned: MakeShell (BOPAlgo_ShellSplitter.cxx L683-698).
/// Creates a TopoDS_Shell from a list of faces.  In rcad this
/// returns a Vec<usize> of face indices.
/// ✅ OCCT-aligned: also calls OrientFacesOnShell.
fn make_shell(faces: &[usize], ds: &DS) -> Vec<usize> {
    let mut shell = faces.to_vec();
    boptools::orient_faces_on_shell(&mut shell, ds);
    shell
}

// ============================================================================
// MakeConnexityBlocks — OCCT BOPTools_AlgoTools::MakeConnexityBlocks (List)
// ============================================================================

/// OCCT-aligned: connectivity-block construction from face indices.
///
/// Maps OCCT BOPTools_AlgoTools::MakeConnexityBlocks (L187-256):
///   1. Build compound from unique shapes; detect duplicates (aMNRegular).
///   2. Call generic MakeConnexityBlocks for edge-connectivity grouping.
///   3. For each block, check multi-connected edges → IsRegular flag.
///
/// In rcad: same logic using DS face/edge indices instead of TopoDS_Shape.
fn make_connexity_blocks(start_shapes: &[usize], ds: &DS, lcb: &mut Vec<ConnexityBlock>) {
    lcb.clear();
    if start_shapes.is_empty() {
        return;
    }

    // --- Step 1: Detect duplicates (OCCT L197-211) ---
    let mut unique_map: HashMap<usize, i32> = HashMap::new();
    for &fi in start_shapes {
        *unique_map.entry(fi).or_insert(0) += 1;
    }
    let a_mn_regular: HashSet<usize> = unique_map
        .iter()
        .filter(|(_, cnt)| **cnt > 1)
        .map(|(&fi, _)| fi)
        .collect();

    // --- Step 2: Build edge-face adjacency for generic MakeConnexityBlocks ---
    // OCCT builds a compound + uses TopExp::MapShapesAndAncestors.
    let unique_faces: Vec<usize> = unique_map.keys().copied().collect();
    let ef_map = build_ef_map(&unique_faces, ds);

    // --- Step 3: BFS to find connected components ---
    let n = unique_faces.len();
    if n == 0 {
        return;
    }

    // Build local adjacency: face → neighbor faces through shared edges
    let face_to_idx: HashMap<usize, usize> = unique_faces
        .iter()
        .enumerate()
        .map(|(i, &fi)| (fi, i))
        .collect();
    // For each face, find its neighbor faces (share an edge)
    // We'll BFS over faces, using ef_map to find neighbors
    let mut ef_map_local: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &fi in &unique_faces {
        let mut nb = Vec::new();
        // Collect all edges of this face
        if let Some(face) = ds.faces.get(fi) {
            let all_edges: Vec<usize> = face
                .boundary_edges
                .iter()
                .copied()
                .chain(
                    face.inner_boundary_edges
                        .iter()
                        .flat_map(|w| w.iter().map(|&(e, _)| e)),
                )
                .collect();
            for &ei in &all_edges {
                if let Some(flist) = ef_map.get(&ei) {
                    for &nb_fi in flist {
                        if nb_fi != fi {
                            nb.push(nb_fi);
                        }
                    }
                }
            }
        }
        // Deduplicate neighbors
        nb.sort_unstable();
        nb.dedup();
        ef_map_local.insert(fi, nb);
    }

    let mut visited = vec![false; n];
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for i in 0..n {
        if visited[i] {
            continue;
        }
        let start_fi = unique_faces[i];
        let mut block = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start_fi);
        visited[face_to_idx[&start_fi]] = true;
        while let Some(fi) = queue.pop_front() {
            block.push(fi);
            if let Some(nb_list) = ef_map_local.get(&fi) {
                for &nb_fi in nb_list {
                    let nb_idx = face_to_idx[&nb_fi];
                    if !visited[nb_idx] {
                        visited[nb_idx] = true;
                        queue.push_back(nb_fi);
                    }
                }
            }
        }
        if !block.is_empty() {
            blocks.push(block);
        }
    }

    // --- Step 4: Create ConnexityBlock per block (OCCT L219-255) ---
    for block in blocks {
        let mut cb = ConnexityBlock::new();
        let mut b_regular = true;
        let cb_shapes = cb.change_shapes();

        for &fi in &block {
            if a_mn_regular.contains(&fi) {
                // Duplicate face → non-regular, add with both orientations
                b_regular = false;
                cb_shapes.push(fi);
                cb_shapes.push(fi);
            } else {
                cb_shapes.push(fi);
                if b_regular {
                    // Check multi-connected edges: does any edge of this face
                    // appear in >2 faces?
                    if let Some(face) = ds.faces.get(fi) {
                        let all_edges: Vec<usize> = face
                            .boundary_edges
                            .iter()
                            .copied()
                            .chain(face.inner_boundary_edges.iter().flat_map(|w| {
                                w.iter().map(|&(e, _)| e)
                            }))
                            .collect();
                        for &ei in &all_edges {
                            if let Some(flist) = ef_map.get(&ei) {
                                if flist.len() > 2 {
                                    b_regular = false;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        cb.set_regular(b_regular);
        lcb.push(cb);
    }
}

// ============================================================================
// RefineShell — OCCT BOPAlgo_ShellSplitter::RefineShell (L443-617)
// ============================================================================

/// OCCT-aligned: split a shell on multi-connected (branch) edges.
///
/// Finds edges with >2 adjacent faces (or same-orientation on both sides),
/// then flood-fills through non-branch edges to create sub-shells.
fn refine_shell(shell_faces: &[usize], ds: &DS) -> Vec<Vec<usize>> {
    if shell_faces.is_empty() {
        return Vec::new();
    }

    // Build edge-face map for this shell
    let ef = build_ef_map(shell_faces, ds);

    // --- Find branch edges (aMEStop) — OCCT L457-512 ---
    // OCCT checks:
    //   (a) Extent > 2 → branch
    //   (b) Extent == 2, same orientation on both sides → branch
    //   (c) INTERNAL edge count > 2 → branch
    let mut branch_edges: HashSet<usize> = HashSet::new();

    for (&ei, flist) in &ef {
        if flist.len() > 2 {
            branch_edges.insert(ei);
            continue;
        }
        if flist.len() == 2 {
            let f1 = flist[0];
            let f2 = flist[1];
            // Check orientation: if both faces have the edge in the same
            // orientation (both forward or both reverse), it's a branch edge.
            // In rcad, orientation in face is determined by whether the
            // edge is in boundary_edges (forward=true) vs inner wire (forward flag).
            let or1 = {
                let face = &ds.faces[f1];
                if face.boundary_edges.contains(&ei) {
                    true // forward
                } else {
                    // Look in inner wires
                    let mut found = false;
                    let mut fwd = true;
                    for wire in &face.inner_boundary_edges {
                        for &(e, f) in wire {
                            if e == ei {
                                found = true;
                                fwd = f;
                                break;
                            }
                        }
                        if found {
                            break;
                        }
                    }
                    fwd
                }
            };
            let or2 = {
                let face = &ds.faces[f2];
                if face.boundary_edges.contains(&ei) {
                    true
                } else {
                    let mut found = false;
                    let mut fwd = true;
                    for wire in &face.inner_boundary_edges {
                        for &(e, f) in wire {
                            if e == ei {
                                found = true;
                                fwd = f;
                                break;
                            }
                        }
                        if found {
                            break;
                        }
                    }
                    fwd
                }
            };
            if or1 == or2 {
                branch_edges.insert(ei);
                continue;
            }
        }
        // Check for INTERNAL edges — count faces where this edge is internal
        // OCCT L486-511: count faces; for each face where the edge is INTERNAL, count +1.
        // In rcad: edge is_internal flag means it counts twice.
        let mut nb_f = flist.len();
        if ds.edges.get(ei).map_or(false, |e| e.is_internal) {
            nb_f += flist.len();
        }
        if nb_f > 2 {
            branch_edges.insert(ei);
        }
    }

    if branch_edges.is_empty() {
        // No branch edges → shell stays as-is
        return vec![shell_faces.to_vec()];
    }

    // --- Flood-fill through non-branch edges (OCCT L520-616) ---
    // For each unprocessed face, start a new sub-shell and BFS through
    // edges that are NOT in branch_edges.
    let mut processed: HashSet<usize> = HashSet::new();
    let mut result: Vec<Vec<usize>> = Vec::new();

    for &start_fi in shell_faces {
        if !processed.insert(start_fi) {
            continue;
        }
        let mut sub_shell: Vec<usize> = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start_fi);
        while let Some(fi) = queue.pop_front() {
            sub_shell.push(fi);
            // Walk all edges of this face; if edge is not a branch edge,
            // add neighbor faces (that are in this shell and unprocessed).
            if let Some(face) = ds.faces.get(fi) {
                let all_edges: Vec<usize> = face
                    .boundary_edges
                    .iter()
                    .copied()
                    .chain(
                        face.inner_boundary_edges
                            .iter()
                            .flat_map(|w| w.iter().map(|&(e, _)| e)),
                    )
                    .collect();
                for &ei in &all_edges {
                    if branch_edges.contains(&ei) {
                        continue;
                    }
                    if let Some(flist) = ef.get(&ei) {
                        for &nb_fi in flist {
                            if nb_fi != fi && processed.insert(nb_fi) {
                                queue.push_back(nb_fi);
                            }
                        }
                    }
                }
            }
        }
        if !sub_shell.is_empty() {
            result.push(sub_shell);
        }
    }
    result
}

// ============================================================================
// ShellSplitter — OCCT BOPAlgo_ShellSplitter
// ============================================================================

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct ShellSplitter {
    /// OCCT: myStartShapes — list of faces to process (TopoDS_Shape list).
    myStartShapes: Vec<usize>,
    /// OCCT: myShells — result shells (list of TopoDS_Shell/Shape).
    myShells: Vec<Vec<usize>>,
    /// OCCT: myLCB — connectivity blocks from MakeConnexityBlocks.
    myLCB: Vec<ConnexityBlock>,
}

impl ShellSplitter {
    /// Empty constructor.
    /// ✅ OCCT-aligned: BOPAlgo_ShellSplitter().
    pub fn new() -> Self {
        Self {
            myStartShapes: Vec::new(),
            myShells: Vec::new(),
            myLCB: Vec::new(),
        }
    }

    /// Add a face (by DS index) to the processing set.
    /// ✅ OCCT-aligned: AddStartElement(const TopoDS_Shape&).
    pub fn add_start_element(&mut self, fi: usize) {
        self.myStartShapes.push(fi);
    }

    /// Return the start elements.
    /// ✅ OCCT-aligned: StartElements().
    pub fn start_elements(&self) -> &[usize] {
        &self.myStartShapes
    }

    /// Perform the shell splitting algorithm.
    /// ✅ OCCT-aligned: Perform (L137-149):
    ///   1. MakeConnexityBlocks (edge-connectivity grouping)
    ///   2. MakeShells (regular→MakeShell, irregular→SplitBlock)
    pub fn perform(&mut self, ds: &DS) {
        self.myShells.clear();
        make_connexity_blocks(&self.myStartShapes, ds, &mut self.myLCB);
        self.make_shells(ds);
    }

    /// ✅ OCCT-aligned: MakeShells (L621-679).
    ///   Regular blocks → make_shell (single shell from all faces).
    ///   Irregular blocks → SplitBlock via parallel CBK vector.
    fn make_shells(&mut self, ds: &DS) {
        // OCCT L632-655: separate regular → MakeShell, irregular → CBK vector
        let mut a_vcbk: VectorOfCbk = Vec::new();
        for i in 0..self.myLCB.len() {
            if self.myLCB[i].is_regular() {
                // OCCT L643-648: MakeShell + Closed(true)
                let faces = self.myLCB[i].shapes().to_vec();
                if !faces.is_empty() {
                    self.myShells.push(make_shell(&faces, ds));
                }
            } else {
                // OCCT L652-654: create CBK entry for SplitBlock
                let mut a_cbk = SplitBlockCbk::new();
                let cb_faces = self.myLCB[i].shapes().to_vec();
                if !cb_faces.is_empty() {
                    a_cbk.set_faces(cb_faces);
                    a_vcbk.push(a_cbk);
                }
            }
        }

        // OCCT L657-665: parallel execution via CBK vector
        // ✅ OCCT-aligned: sequential execution (OCCT parallel); determinism preferred.
        for cbk in a_vcbk.iter_mut() {
            cbk.perform(ds);
        }

        // OCCT L666-678: collect results from CBK vector
        for mut cbk in a_vcbk {
            let loops = cbk.myLoops.drain(..).collect::<Vec<_>>();
            for loop_faces in loops {
                // OCCT L674-676: mark Closed(true)
                self.myShells.push(loop_faces);
            }
        }
    }

    /// Return the resulting shells.
    /// ✅ OCCT-aligned: Shells().
    pub fn shells(&self) -> &[Vec<usize>] {
        &self.myShells
    }

    /// Number of resulting shells.
    pub fn nb_shells(&self) -> usize {
        self.myShells.len()
    }

    /// True when more than one shell was produced.
    pub fn has_multiple_shells(&self) -> bool {
        self.myShells.len() > 1
    }

    /// Clear all state.
    /// ✅ OCCT-aligned: Clear().
    pub fn clear(&mut self) {
        self.myStartShapes.clear();
        self.myShells.clear();
        self.myLCB.clear();
    }
}

impl Default for ShellSplitter {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// SplitBlock — OCCT BOPAlgo_ShellSplitter::SplitBlock (L153-421)
// ============================================================================

/// OCCT-aligned: SplitBlock(BOPTools_ConnexityBlock&).
///
/// Processes an irregular connectivity block:
///   1. Remove faces with free edges (iterative)
///   2. Build shells from remaining faces via flood-fill through shared edges
///   3. RefineShell to split on multi-connected edges
///   4. Collect closed shells into myLoops
///
/// ✅ OCCT-aligned: GetFaceOff angle selection via boptools::get_face_off.
fn split_block(faces: &[usize], ds: &DS, loops: &mut Vec<Vec<usize>>) {
    if faces.is_empty() {
        return;
    }

    let my_shapes: Vec<usize> = faces.to_vec();

    // --- Phase 1: Copy faces into aMFaces map (OCCT L177-182) ---
    let mut a_m_faces: BTreeSet<usize> = my_shapes.iter().copied().collect();

    // --- Phase 2: Remove faces with free edges (OCCT L185-222) ---
    loop {
        // Build edge-face map for remaining faces
        let remaining_faces: Vec<usize> = a_m_faces.iter().copied().collect();
        let ef = build_ef_map(&remaining_faces, ds);

        let a_nb_begin = a_m_faces.len();

        // Check free edges: edges that appear in only 1 face and are
        // not degenerated and not INTERNAL → remove that face
        let mut to_remove: Vec<usize> = Vec::new();
        for (&ei, flist) in &ef {
            if !ds.is_edge_degenerated(ei)
                && ds.edges.get(ei).map_or(true, |e| !e.is_internal)
                && flist.len() == 1
            {
                to_remove.push(flist[0]);
            }
        }
        for &fi in &to_remove {
            a_m_faces.remove(&fi);
        }

        let a_nb_end = a_m_faces.len();
        if a_nb_end == a_nb_begin || a_nb_end == 0 {
            break;
        }
    }

    if a_m_faces.is_empty() {
        return;
    }

    // --- Phase 3: Connected faces only (OCCT L229-245) ---
    let connected_faces: Vec<usize> = my_shapes
        .iter()
        .copied()
        .filter(|fi| a_m_faces.contains(fi))
        .collect();

    // Boundary faces: those that were NOT removed by free-edge elimination
    let mut boundary_faces: BTreeSet<usize> = BTreeSet::new();
    for &fi in &connected_faces {
        if a_m_faces.contains(&fi) {
            if !boundary_faces.insert(fi) {
                boundary_faces.remove(&fi);
            }
        }
    }

    let a_nb_shapes = connected_faces.len();
    let mut b_all_faces_taken = false;
    let mut added_faces: BTreeSet<usize> = BTreeSet::new();

    // --- Phase 4: Build shells (OCCT L250-420) ---
    for &seed_fi in &connected_faces {
        if b_all_faces_taken {
            break;
        }
        if !added_faces.insert(seed_fi) {
            continue;
        }

        // --- 4a: Make new shell with seed face (OCCT L260-263) ---
        let mut shell_faces: Vec<usize> = vec![seed_fi];

        // --- 4b: Build edge-face map for this shell (OCCT L265-266) ---
        let mut mefp = build_ef_map(&shell_faces, ds);

        // --- 4c: Iterate shell faces, add neighbors (OCCT L270-369) ---
        let mut shell_idx = 0;
        while shell_idx < shell_faces.len() {
            let a_f = shell_faces[shell_idx];
            let is_boundary = boundary_faces.contains(&a_f);

            // Iterate edges of the face (OCCT L277: TopExp_Explorer)
            let all_edges: Vec<usize> = {
                let face = &ds.faces[a_f];
                face.boundary_edges
                    .iter()
                    .copied()
                    .chain(
                        face.inner_boundary_edges
                            .iter()
                            .flat_map(|w| w.iter().map(|&(e, _)| e)),
                    )
                    .collect()
            };

            for &ei in &all_edges {
                // Skip edges already with 2+ faces in this shell (OCCT L283-291)
                if let Some(flist) = mefp.get(&ei) {
                    if flist.len() > 1 {
                        continue;
                    }
                }

                // Skip INTERNAL edges (OCCT L293-297)
                if ds.edges.get(ei).map_or(false, |e| e.is_internal) {
                    continue;
                }

                // Skip degenerated edges (OCCT L299-302)
                if ds.is_edge_degenerated(ei) {
                    continue;
                }

                // Candidate faces from global ef map (OCCT L305-310)
                let remaining: Vec<usize> = a_m_faces.iter().copied().collect();
                let ef_global = build_ef_map(&remaining, ds);
                let Some(candidates) = ef_global.get(&ei) else {
                    continue;
                };

                // --- Select the next face (OCCT L312-361) ---
                let mut a_nb_ways_inside: usize = 0;
                let mut a_sel_f: Option<usize> = None;
                // CoupleOfShape list: (edge, face) pairs
                let mut lcs_off: Vec<(usize, usize)> = Vec::new();

                for &a_fl in candidates {
                    if a_f == a_fl || added_faces.contains(&a_fl) {
                        continue;
                    }
                    // GetEdgeOff: ensure the edge exists in this face
                    if get_edge_off(ei, a_fl, ds).is_none() {
                        continue;
                    }
                    if is_boundary && !boundary_faces.contains(&a_fl) {
                        a_nb_ways_inside += 1;
                        a_sel_f = Some(a_fl);
                    }
                    lcs_off.push((ei, a_fl));
                }

                let a_nb_off = lcs_off.len();
                if a_nb_off == 0 {
                    continue;
                }

                // OCCT L349-361: among all adjacent faces, select one with
                // minimal angle to the current face (GetFaceOff).
                if !is_boundary || a_nb_ways_inside != 1 {
                    if a_nb_off == 1 {
                        a_sel_f = Some(lcs_off[0].1);
                    } else if a_nb_off > 1 {
                        a_sel_f = boptools::get_face_off(ei, a_f, &lcs_off, ds);
                    }
                }

                if let Some(sel_f) = a_sel_f {
                    if added_faces.insert(sel_f) {
                        shell_faces.push(sel_f);
                        // Update edge-face map for this shell
                        let ef_sel = build_ef_map(&[sel_f], ds);
                        for (e, flist) in ef_sel {
                            mefp.entry(e).or_default().extend(flist);
                        }
                    }
                }
            }
            shell_idx += 1;
        }

        // --- 4d: RefineShell — split on multi-connected edges (OCCT L371-373) ---
        let sub_shells = refine_shell(&shell_faces, ds);

        // --- 4e: Classify sub-shells as closed or not (OCCT L375-392) ---
        let mut a_lsh_nc: Vec<Vec<usize>> = Vec::new();
        for sub in sub_shells {
            if is_shell_closed(&sub, ds) {
                loops.push(sub);
            } else {
                a_lsh_nc.push(sub);
            }
        }

        // --- 4f: Check if all faces taken (OCCT L394-398) ---
        b_all_faces_taken = added_faces.len() == a_nb_shapes;
        if b_all_faces_taken {
            break;
        }

        // --- 4g: If only 1 sub-shell and not closed → can't improve (OCCT L400-405) ---
        if a_lsh_nc.len() == 1 {
            continue;
        }

        // --- 4h: Remove not-closed shell faces from added_faces (OCCT L408-419) ---
        for nc_shell in &a_lsh_nc {
            for &fi in nc_shell {
                added_faces.remove(&fi);
            }
        }
    }
}
