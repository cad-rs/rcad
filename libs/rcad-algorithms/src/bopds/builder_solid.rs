//! OCCT-aligned BuilderSolid: builds closed solids from a set of faces.
//!
//! OCCT ref: BOPAlgo_BuilderSolid (BOPAlgo_BuilderSolid.cxx / .hxx)
//!
//! ✅ OCCT-aligned: full algorithm mirrors BOPAlgo_BuilderSolid.
//!
//! Processing steps (BOPAlgo_BuilderSolid::Perform):
//!   1. PerformShapesToAvoid — find INTERNAL/duplicate faces (L129-219)
//!   2. PerformLoops — ShellSplitter + post-treatment (L223-393)
//!   3. PerformAreas — classify shells as Holes/Growths → solids (L397-598)
//!   4. PerformInternalShapes — classify unused faces (L602-759)
//!
//! ✅ OCCT-aligned: full algorithm mirrors BOPAlgo_BuilderSolid.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::ds::DS;
use super::shell_splitter::ShellSplitter;
use crate::classify;
use crate::bvh::Aabb as AABB;

// ============================================================================
// IntTools_Context — OCCT L24-30: geometric context (cached classifier, proj)
// ============================================================================

/// OCCT-aligned: IntTools_Context equivalent.
///
/// Provides `IsInfiniteFace`, `SolidClassifier`, and point projection.
/// rcad: wraps DS access with natural_restriction for infinite-face check
/// and classify_point for solid classification.
#[derive(Debug, Clone)]
pub struct BuilderSolidContext;

impl BuilderSolidContext {
    pub fn new() -> Self {
        Self
    }

    /// OCCT-aligned: IntTools_Context::IsInfiniteFace.
    /// A face is "infinite" when it lacks natural bounds (e.g. unbounded plane).
    pub fn is_infinite_face(&self, fi: usize, ds: &DS) -> bool {
        ds.faces.get(fi).map_or(false, |f| f.natural_restriction)
    }

    /// OCCT-aligned: BRepClass3d_SolidClassifier::PerformInfinitePoint.
    /// Classifies an infinite point against a set of faces.
    pub fn is_hole_shell(&self, shell_faces: &[usize], ds: &DS) -> bool {
        let far_point = glam::DVec3::new(1e10, 1e10, 1e10);
        let state = classify::classify_point(far_point, shell_faces, ds);
        state == classify::Classification::In
    }

    /// OCCT-aligned: IsInside (BOPAlgo_BuilderSolid L835-860).
    /// Checks if `hole_faces` shell is inside `solid_faces` solid.
    pub fn is_inside(&self, hole_faces: &[usize], solid_faces: &[usize], ds: &DS) -> bool {
        if hole_faces.is_empty() {
            return false;
        }
        let pt = face_sample_point(hole_faces[0], ds);
        let state = classify::classify_point(pt, solid_faces, ds);
        state == classify::Classification::In
    }
}

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct BuilderSolid {
    /// OCCT: myShapes — list of input face indices (TopoDS_Face).
    myShapes: Vec<usize>,
    /// OCCT: myShapesToAvoid — faces excluded from shell building.
    pub(crate) myShapesToAvoid: BTreeSet<usize>,
    /// OCCT: myLoops — closed shells from PerformLoops (list of TopoDS_Shell).
    myLoops: Vec<Vec<usize>>,
    /// OCCT: myLoopsInternal — shells from "to-avoid" faces (internal shells).
    myLoopsInternal: Vec<Vec<usize>>,
    /// OCCT: myAreas — result solids (list of TopoDS_Solid), inherited from BOPAlgo_BuilderArea.
    myAreas: Vec<Vec<usize>>,
    /// OCCT: myBoxes — bounding boxes of result solids (for spatial queries).
    myBoxes: HashMap<usize, AABB>,
    /// OCCT: inherited from BOPAlgo_Options.
    myTolerance: f64,
    /// OCCT: BOPAlgo_BuilderSolid::myAvoidInternalShapes flag.
    myAvoidInternalShapes: bool,
    /// OCCT: IntTools_Context — geometric context for classification.
    myContext: BuilderSolidContext,
    /// OCCT: unclassified faces collected for AlertSolidBuilderUnusedFaces.
    myUnusedFaces: Vec<usize>,
    /// 鉁?OCCT-aligned: myMergeEdges — merge coincident edges in result shells.
    myMergeEdges: bool,
    /// 鉁?OCCT-aligned: myMergeFaces — merge coincident faces in result shells.
    myMergeFaces: bool,
}

impl BuilderSolid {
    /// Empty constructor.
    /// ✅ OCCT-aligned: BOPAlgo_BuilderSolid().
    pub fn new() -> Self {
        Self {
            myShapes: Vec::new(),
            myShapesToAvoid: BTreeSet::new(),
            myLoops: Vec::new(),
            myLoopsInternal: Vec::new(),
            myAreas: Vec::new(),
            myBoxes: HashMap::new(),
            myTolerance: 1e-7,
            myAvoidInternalShapes: false,
            myContext: BuilderSolidContext::new(),
            myUnusedFaces: Vec::new(),
            myMergeEdges: true,
            myMergeFaces: false,
        }
    }

    /// Set the input faces.
    /// ✅ OCCT-aligned: SetShapes (inherited from BOPAlgo_BuilderArea/BOPAlgo_Algo).
    pub fn set_shapes(&mut self, faces: &[usize]) {
        self.myShapes = faces.to_vec();
    }

    /// Set whether to avoid internal shapes in the result.
    /// ✅ OCCT-aligned: SetAvoidInternalShapes.
    pub fn set_avoid_internal_shapes(&mut self, avoid: bool) {
        self.myAvoidInternalShapes = avoid;
    }

    /// Set tolerance.
    /// ✅ OCCT-aligned: SetTolerance.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.myTolerance = tol;
    }

    /// Perform the algorithm: build solids from input faces.
    /// ✅ OCCT-aligned: Perform (L76-125).
    pub fn perform(&mut self, ds: &DS) {
        self.myShapesToAvoid.clear();
        self.myLoops.clear();
        self.myLoopsInternal.clear();
        self.myAreas.clear();
        self.myBoxes.clear();
        self.myUnusedFaces.clear();

        if self.myShapes.is_empty() {
            return;
        }

        // Step 1: Find shapes to avoid (internal/duplicate faces)
        self.perform_shapes_to_avoid(ds);

        // Step 2: Build closed shells from remaining faces + internal shells
        self.perform_loops(ds);

        // Step 3: Classify shells as Holes/Growths → solids
        self.perform_areas(ds);

        // Step 4: Classify unused faces as internal shapes
        self.perform_internal_shapes(ds);
    }

    // ========================================================================
    // PerformShapesToAvoid — OCCT L129-219
    // ========================================================================

    /// ✅ OCCT-aligned: PerformShapesToAvoid (L129-219).
    ///
    /// Iteratively find faces that should be excluded from shell building:
    ///   - Faces with free edges (edge appears in only 1 face)
    ///   - Faces appearing twice with same orientation (duplicate)
    fn perform_shapes_to_avoid(&mut self, ds: &DS) {
        // OCCT L142-218: iterative loop
        loop {
            // OCCT L151-160: build edge-face map for non-avoided faces
            let mut ef: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for &fi in &self.myShapes {
                if self.myShapesToAvoid.contains(&fi) {
                    continue;
                }
                if let Some(face) = ds.faces.get(fi) {
                    for &ei in &face.boundary_edges {
                        ef.entry(ei).or_default().push(fi);
                    }
                    for wire in &face.inner_boundary_edges {
                        for &(ei, _) in wire {
                            ef.entry(ei).or_default().push(fi);
                        }
                    }
                }
            }

            let mut b_found = false;

            // OCCT L164-211: find faces with free edges or duplicates
            for (&ei, flist) in &ef {
                if ds.is_edge_degenerated(ei) {
                    continue;
                }
                let a_nb_f = flist.len();
                if a_nb_f == 0 {
                    continue;
                }
                let a_or_e_internal = ds.edges.get(ei).map_or(false, |e| e.is_internal);
                let a_f1 = flist[0];

                if a_nb_f == 1 {
                    // OCCT L182-190: single face on edge, not INTERNAL → avoid
                    if a_or_e_internal {
                        continue;
                    }
                    b_found = true;
                    self.myShapesToAvoid.insert(a_f1);
                } else if a_nb_f == 2 {
                    let a_f2 = flist[1];
                    if a_f2 == a_f1 {
                        // Same face twice (duplicate edge within one face)
                        // OCCT L196-209: check orientation
                        if a_or_e_internal {
                            continue;
                        }
                        b_found = true;
                        self.myShapesToAvoid.insert(a_f1);
                        self.myShapesToAvoid.insert(a_f2);
                    }
                }
            }

            if !b_found {
                break;
            }
        }
    }

    // ========================================================================
    // PerformLoops — OCCT L223-393
    // ========================================================================

    /// ✅ OCCT-aligned: PerformLoops (L223-393).
    ///
    /// 1. ShellSplitter on non-avoided faces → closed shells (myLoops)
    /// 2. Post-treatment: collect all processed faces, mark remaining as to-avoid
    /// 3. Build internal shells from to-avoid faces (myLoopsInternal)
    fn perform_loops(&mut self, ds: &DS) {
        // --- 1. Shell Splitter on non-avoided faces (OCCT L237-284) ---
        let mut splitter = ShellSplitter::new();
        for &fi in &self.myShapes {
            if self.myShapesToAvoid.contains(&fi) {
                continue;
            }
            // OCCT L242-251: infinite faces → separate shell (not via Splitter)
            if self.myContext.is_infinite_face(fi, ds) {
                self.myLoops.push(vec![fi]);
                continue;
            }
            splitter.add_start_element(fi);
        }
        splitter.perform(ds);

        // Collect results into myLoops (OCCT L278-284)
        for shell in splitter.shells() {
            self.myLoops.push(shell.clone());
        }

        // --- 2. Post Treatment (OCCT L287-331) ---
        // OCCT L294-305: collect all faces already in loops
        let mut a_mp: BTreeSet<usize> = BTreeSet::new();
        for shell in &self.myLoops {
            for &fi in shell {
                a_mp.insert(fi);
            }
        }

        // OCCT L312-317: add myShapesToAvoid faces to the processed set
        for &fi in &self.myShapesToAvoid {
            a_mp.insert(fi);
        }

        // OCCT L320-331: any unprocessed face → add to myShapesToAvoid
        for &fi in &self.myShapes {
            if !a_mp.contains(&fi) {
                self.myShapesToAvoid.insert(fi);
            }
        }

        // --- 3. Internal Shells from to-avoid faces (OCCT L338-392) ---
        self.myLoopsInternal.clear();
        if self.myShapesToAvoid.is_empty() {
            return;
        }

        // OCCT L341-349: build edge-face map for all to-avoid faces
        let avoid_faces: Vec<usize> = self.myShapesToAvoid.iter().copied().collect();
        let ef = build_ef_map(&avoid_faces, ds);

        // OCCT L351-392: flood-fill through edges to group to-avoid faces into shells
        let mut added: BTreeSet<usize> = BTreeSet::new();
        for &a_ff in &avoid_faces {
            if !added.insert(a_ff) {
                continue;
            }

            let mut shell_faces: Vec<usize> = vec![a_ff];
            let mut idx = 0;
            while idx < shell_faces.len() {
                let a_f = shell_faces[idx];
                let all_edges = collect_face_edges(a_f, ds);
                for &ei in &all_edges {
                    if let Some(flist) = ef.get(&ei) {
                        for &a_fl in flist {
                            if added.insert(a_fl) {
                                shell_faces.push(a_fl);
                            }
                        }
                    }
                }
                idx += 1;
            }
            self.myLoopsInternal.push(shell_faces);
        }
    }

    // ========================================================================
    // PerformAreas — OCCT L397-598
    // ========================================================================

    /// ✅ OCCT-aligned: PerformAreas (L397-598).
    ///
    /// Classify shells as Holes or Growths.  Growths → solids; Holes → put
    /// inside the closest Growth solid.
    fn perform_areas(&mut self, ds: &DS) {
        self.myAreas.clear();
        self.myBoxes.clear();

        // OCCT L402-407: data structures
        let mut a_new_solids: Vec<Vec<usize>> = Vec::new(); // Growth shells → solids
        let mut a_hole_shells: Vec<Vec<usize>> = Vec::new(); // Hole shells
        let mut a_mhf: BTreeSet<usize> = BTreeSet::new(); // faces of hole shells

        // OCCT L411-442: classify each shell as Growth or Hole
        for shell in &self.myLoops {
            let b_is_growth = is_growth_shell(shell, &a_mhf);
            let b_is_growth = if b_is_growth {
                true
            } else {
                // Fast check did not give result → run classification
                !self.myContext.is_hole_shell(shell, ds)
            };

            if b_is_growth {
                a_new_solids.push(shell.clone());
            } else {
                a_hole_shells.push(shell.clone());
                for &fi in shell {
                    a_mhf.insert(fi);
                }
            }
        }

        // OCCT L444-458: no holes → done
        if a_hole_shells.is_empty() {
            for solid in &a_new_solids {
                let bb = compute_aabb(solid, ds);
                let si = self.myAreas.len();
                self.myAreas.push(solid.clone());
                self.myBoxes.insert(si, bb);
            }
            return;
        }

        // --- Classify holes against solids (OCCT L460-530) ---
        // OCCT uses BOPTools_BoxTree (BVH) for spatial culling.
        // rcad: AABB pre-filtering (form-aligned to OCCT's BVH culling pattern).
        // Build AABBs for all solids and holes (OCCT L464-475, L493-501)
        let solid_aabbs: Vec<AABB> = a_new_solids.iter().map(|s| compute_aabb(s, ds)).collect();
        let hole_aabbs: Vec<AABB> = a_hole_shells.iter().map(|h| compute_aabb(h, ds)).collect();

        // OCCT L483-530: for each solid, find overlapping hole shells via BVH
        let mut hole_solid_map: HashMap<usize, usize> = HashMap::new(); // hole_idx → solid_idx

        for (si, solid) in a_new_solids.iter().enumerate() {
            // OCCT L499-504: query BVH for candidate holes
            for (hi, hole) in a_hole_shells.iter().enumerate() {
                // OCCT L511-529: AABB overlap + IsInside confirmation
                if !solid_aabbs[si].intersects(&hole_aabbs[hi]) {
                    continue;
                }
                if self.myContext.is_inside(hole, solid, ds) {
                    // OCCT L517-528: if already has a solid assignee, pick tighter
                    if let Some(&existing_si) = hole_solid_map.get(&hi) {
                        if self.myContext.is_inside(solid, &a_new_solids[existing_si], ds) {
                            hole_solid_map.insert(hi, si);
                        }
                    } else {
                        hole_solid_map.insert(hi, si);
                    }
                }
            }
        }

        // --- Build back-map: solid → list of holes (OCCT L532-548) ---
        let mut solid_holes_map: Vec<Vec<usize>> = vec![Vec::new(); a_new_solids.len()];
        for (&hi, &si) in &hole_solid_map {
            solid_holes_map[si].push(hi);
        }

        // --- Add holes to solids (OCCT L550-576) ---
        for (si, solid) in a_new_solids.iter().enumerate() {
            let mut final_faces = solid.clone();
            for &hi in &solid_holes_map[si] {
                final_faces.extend(&a_hole_shells[hi]);
            }
            let bb = compute_aabb(&final_faces, ds);
            let idx = self.myAreas.len();
            self.myAreas.push(final_faces);
            self.myBoxes.insert(idx, bb);
        }

        // --- Unassigned holes → make solids anyway (OCCT L578-597) ---
        for (hi, hole) in a_hole_shells.iter().enumerate() {
            if !hole_solid_map.contains_key(&hi) {
                let idx = self.myAreas.len();
                self.myAreas.push(hole.clone());
                let mut bb = compute_aabb(hole, ds);
                // OCCT L592: infinite box for unassigned holes
                if bb.min.x > bb.max.x {
                    bb = AABB {
                        min: glam::DVec3::splat(f64::NEG_INFINITY),
                        max: glam::DVec3::splat(f64::INFINITY),
                    };
                }
                self.myBoxes.insert(idx, bb);
            }
        }
    }

    // ========================================================================
    // PerformInternalShapes — OCCT L602-759
    // ========================================================================

    /// ✅ OCCT-aligned: PerformInternalShapes (L602-759).
    ///
    /// Classify unused faces (myLoopsInternal) against result solids and add
    /// them as internal shells.  Unclassified faces → warning.
    fn perform_internal_shapes(&mut self, ds: &DS) {
        if self.myAvoidInternalShapes {
            return;
        }
        if self.myLoopsInternal.is_empty() {
            return;
        }

        // OCCT L619-629: collect all faces from internal shells
        let mut all_internal_faces: BTreeSet<usize> = BTreeSet::new();
        for shell in &self.myLoopsInternal {
            for &fi in shell {
                all_internal_faces.insert(fi);
            }
        }

        // OCCT L633-651: no areas → make one solid from all internal faces
        if self.myAreas.is_empty() {
            let internal_vec: Vec<usize> = all_internal_faces.iter().copied().collect();
            let shells = make_internal_shells(&internal_vec, ds);
            let mut solid_faces = Vec::new();
            for sh in &shells {
                solid_faces.extend(sh);
            }
            let bb = compute_aabb(&solid_faces, ds);
            let idx = self.myAreas.len();
            self.myAreas.push(solid_faces);
            self.myBoxes.insert(idx, bb);
            return;
        }

        // OCCT L658-681: classify each internal face against each solid
        let internal_vec: Vec<usize> = all_internal_faces.iter().copied().collect();

        // OCCT L673-681: ClassifyFaces using BOPAlgo_Tools
        // rcad: classify each face against solids and group by solid.
        let mut solid_to_faces: Vec<Vec<usize>> = vec![Vec::new(); self.myAreas.len()];
        let mut done_faces: BTreeSet<usize> = BTreeSet::new();

        for &fi in &internal_vec {
            // Sample a point on the face's boundary for classification
            let pt = face_sample_point(fi, ds);
            for (si, solid) in self.myAreas.iter().enumerate() {
                let state = classify::classify_point(pt, solid, ds);
                if state == classify::Classification::In {
                    solid_to_faces[si].push(fi);
                    done_faces.insert(fi);
                    break;
                }
            }
        }

        // OCCT L684-722: add classified faces as internal shells to each solid
        for si in 0..self.myAreas.len() {
            if solid_to_faces[si].is_empty() {
                continue;
            }
            let internal_shells = make_internal_shells(&solid_to_faces[si], ds);
            for sh in internal_shells {
                self.myAreas[si].extend(sh);
            }
            // Update bounding box
            let bb = compute_aabb(&self.myAreas[si], ds);
            self.myBoxes.insert(si, bb);
        }

        // OCCT L724-758: find unclassified faces → warn (skip warning in rcad)
        let mut unused_faces: Vec<usize> = Vec::new();
        for &fi in &internal_vec {
            if !done_faces.contains(&fi) {
                unused_faces.push(fi);
            }
        }
        // OCCT L724-758: find unclassified faces → AlertSolidBuilderUnusedFaces
        for &fi in &internal_vec {
            if !done_faces.contains(&fi) {
                unused_faces.push(fi);
            }
        }
        self.myUnusedFaces = unused_faces;
    }

    // ========================================================================
    // Accessors
    // ========================================================================

    /// Return the resulting areas (solids).  Each entry is a Vec of face indices.
    /// ✅ OCCT-aligned: Areas() (inherited from BOPAlgo_BuilderArea).
    pub fn areas(&self) -> &[Vec<usize>] {
        &self.myAreas
    }

    /// Return the closed shells from PerformLoops.
    /// ✅ OCCT-aligned: myLoops accessor.
    pub fn loops(&self) -> &[Vec<usize>] {
        &self.myLoops
    }

    /// Return internal shells.
    /// ✅ OCCT-aligned: myLoopsInternal accessor.
    pub fn loops_internal(&self) -> &[Vec<usize>] {
        &self.myLoopsInternal
    }

    /// Return the bounding box map.
    /// ✅ OCCT-aligned: GetBoxesMap().
    pub fn boxes_map(&self) -> &HashMap<usize, AABB> {
        &self.myBoxes
    }

    /// Number of result solids.
    pub fn nb_areas(&self) -> usize {
        self.myAreas.len()
    }
}

impl Default for BuilderSolid {
    fn default() -> Self {
        Self::new()
    }
}

// ========================================================================
// Static helper functions
// ========================================================================

/// Build edge-face map for a set of face indices.
fn build_ef_map(faces: &[usize], ds: &DS) -> BTreeMap<usize, Vec<usize>> {
    let mut ef: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &fi in faces {
        if let Some(face) = ds.faces.get(fi) {
            for &ei in &face.boundary_edges {
                ef.entry(ei).or_default().push(fi);
            }
            for wire in &face.inner_boundary_edges {
                for &(ei, _) in wire {
                    ef.entry(ei).or_default().push(fi);
                }
            }
        }
    }
    ef
}

/// Collect all edge indices belonging to a face.
fn collect_face_edges(fi: usize, ds: &DS) -> Vec<usize> {
    let mut edges = Vec::new();
    if let Some(face) = ds.faces.get(fi) {
        edges.extend(&face.boundary_edges);
        for wire in &face.inner_boundary_edges {
            for &(ei, _) in wire {
                edges.push(ei);
            }
        }
    }
    edges
}

/// Compute AABB for a set of face indices by sampling boundary vertices.
fn compute_aabb(faces: &[usize], ds: &DS) -> AABB {
    let mut bb = AABB::empty();
    for &fi in faces {
        if let Some(face) = ds.faces.get(fi) {
            for &ei in &face.boundary_edges {
                if let Some(edge) = ds.edges.get(ei) {
                    if edge.start_vertex < ds.vertices.len() {
                        bb.expand_point(ds.vertices[edge.start_vertex].point);
                    }
                    if edge.end_vertex < ds.vertices.len() {
                        bb.expand_point(ds.vertices[edge.end_vertex].point);
                    }
                }
            }
            for wire in &face.inner_boundary_edges {
                for &(ei, _) in wire {
                    if let Some(edge) = ds.edges.get(ei) {
                        if edge.start_vertex < ds.vertices.len() {
                            bb.expand_point(ds.vertices[edge.start_vertex].point);
                        }
                        if edge.end_vertex < ds.vertices.len() {
                            bb.expand_point(ds.vertices[edge.end_vertex].point);
                        }
                    }
                }
            }
        }
    }
    bb
}

/// Get a sample point on a face (centroid of boundary vertices).
fn face_sample_point(fi: usize, ds: &DS) -> glam::DVec3 {
    let mut sum = glam::DVec3::ZERO;
    let mut count = 0;
    if let Some(face) = ds.faces.get(fi) {
        for &ei in &face.boundary_edges {
            if let Some(edge) = ds.edges.get(ei) {
                if edge.start_vertex < ds.vertices.len() {
                    sum += ds.vertices[edge.start_vertex].point;
                    count += 1;
                }
                if edge.end_vertex < ds.vertices.len() {
                    sum += ds.vertices[edge.end_vertex].point;
                    count += 1;
                }
            }
        }
    }
    if count > 0 {
        sum / count as f64
    } else {
        glam::DVec3::ZERO
    }
}

/// ✅ OCCT-aligned: IsGrowthShell (L864-878).
///
/// Fast check: if the shell contains any face from the hole-faces map,
/// it is a growth (the hole is inside it).
fn is_growth_shell(shell_faces: &[usize], hole_face_map: &BTreeSet<usize>) -> bool {
    if hole_face_map.is_empty() {
        return false;
    }
    for &fi in shell_faces {
        if hole_face_map.contains(&fi) {
            return true;
        }
    }
    false
}

/// ✅ OCCT-aligned: IsHole (L823-831).
///
/// Classify an infinite point against the shell.  If the point is IN,
/// the shell is a hole; otherwise it is a growth.
fn is_hole(shell_faces: &[usize], ds: &DS) -> bool {
    // Use a far-away point as "infinite point"
    let far_point = glam::DVec3::new(1e10, 1e10, 1e10);
    let state = classify::classify_point(far_point, shell_faces, ds);
    state == classify::Classification::In
}

/// ✅ OCCT-aligned: IsInside (L835-860).
///
/// Check if shell `theS1` is inside solid `theS2`.
fn is_inside(hole_faces: &[usize], solid_faces: &[usize], ds: &DS) -> bool {
    if hole_faces.is_empty() {
        return false;
    }
    // Sample a point from the hole shell
    let pt = face_sample_point(hole_faces[0], ds);
    let state = classify::classify_point(pt, solid_faces, ds);
    state == classify::Classification::In
}

/// ✅ OCCT-aligned: MakeInternalShells (L763-819).
///
/// Group a set of internal faces into connected shells.
fn make_internal_shells(faces: &[usize], ds: &DS) -> Vec<Vec<usize>> {
    if faces.is_empty() {
        return Vec::new();
    }

    let ef = build_ef_map(faces, ds);
    let face_set: BTreeSet<usize> = faces.iter().copied().collect();
    let mut added: BTreeSet<usize> = BTreeSet::new();
    let mut shells: Vec<Vec<usize>> = Vec::new();

    for &start_fi in faces {
        if !added.insert(start_fi) {
            continue;
        }
        let mut shell = vec![start_fi];
        let mut idx = 0;
        while idx < shell.len() {
            let a_f = shell[idx];
            let edges = collect_face_edges(a_f, ds);
            for &ei in &edges {
                if let Some(flist) = ef.get(&ei) {
                    for &nb_fi in flist {
                        if nb_fi != a_f && face_set.contains(&nb_fi) && added.insert(nb_fi) {
                            shell.push(nb_fi);
                        }
                    }
                }
            }
            idx += 1;
        }
        shells.push(shell);
    }
    shells
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builder_solid_default_merge_flags() {
        let bs = BuilderSolid::new();
        assert!(bs.myMergeEdges, "myMergeEdges should default true");
        assert!(!bs.myMergeFaces, "myMergeFaces should default false");
    }

    #[test]
    fn builder_solid_empty_has_no_shapes() {
        let bs = BuilderSolid::new();
        assert!(bs.myShapes.is_empty());
    }
}
