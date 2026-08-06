// OCCT BOPAlgo_ShellSplitter — shell partitioning.
//
// OCCT BOPAlgo_ShellSplitter.cxx
// Performs: MakeConnexityBlocks -> MakeShells (regular blocks directly,
// non-regular blocks via SplitBlock + RefineShell).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, Orientation, TShape};
use std::collections::{HashMap, HashSet, VecDeque};
use glam::DVec3;

/// OCCT BOPTools_ConnexityBlock — a connected block of faces with regularity.
#[derive(Debug, Clone)]
pub struct ConnexityBlock {
    pub shapes: Vec<Shape>,
    pub regular: bool,
}

/// OCCT BOPTools_AlgoTools::MakeConnexityBlocks (BOPTools_AlgoTools.cxx L187-256).
/// Builds connected blocks of faces; multi-connected faces (appearing more than
/// once in the input) are expanded into both orientations and make the block
/// non-regular.
pub fn make_connexity_blocks(shapes: &[Shape]) -> Vec<ConnexityBlock> {
    // OCCT L197-211: aMFence/aMNRegular — faces appearing more than once.
    let mut a_mfence: HashSet<u64> = HashSet::new();
    let mut a_mn_regular: HashSet<u64> = HashSet::new();
    let mut a_start: Vec<Shape> = Vec::new();
    for s in shapes {
        if !a_mfence.insert(s.ptr_id()) {
            a_mn_regular.insert(s.ptr_id());
        } else {
            a_start.push(s.clone());
        }
    }
    // Edge -> [face indices] map (OCCT TopExp::MapShapesAndAncestors).
    let mut a_cmap: HashMap<u64, Vec<usize>> = HashMap::new();
    for (i, s) in a_start.iter().enumerate() {
        for eptr in face_edge_ptrs(s) {
            a_cmap.entry(eptr).or_default().push(i);
        }
    }
    // BFS grouping by shared edges (OCCT L213-216).
    let n = a_start.len();
    let mut visited = vec![false; n];
    let mut blocks: Vec<Vec<usize>> = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut block: Vec<usize> = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;
        while let Some(i) = queue.pop_front() {
            block.push(i);
            for eptr in face_edge_ptrs(&a_start[i]) {
                if let Some(nei) = a_cmap.get(&eptr) {
                    for &ni in nei {
                        if !visited[ni] {
                            visited[ni] = true;
                            queue.push_back(ni);
                        }
                    }
                }
            }
        }
        if !block.is_empty() {
            blocks.push(block);
        }
    }
    // OCCT L219-255: save the blocks and check their regularity.
    let mut result = Vec::new();
    for block in &blocks {
        let mut a_lcs: Vec<Shape> = Vec::new();
        let mut b_regular = true;
        for &idx in block {
            let mut a_s = a_start[idx].clone();
            if a_mn_regular.contains(&a_s.ptr_id()) {
                // OCCT L231-238: multi-connected shape — both orientations.
                b_regular = false;
                a_s.orientation = Orientation::Forward;
                a_lcs.push(a_s.clone());
                a_s.orientation = Orientation::Reversed;
                a_lcs.push(a_s);
            } else {
                a_lcs.push(a_s);
                if b_regular {
                    // OCCT L243-247: no multi-connected shapes — check every
                    // connection edge is used by exactly 2 elements.
                    for eptr in face_edge_ptrs(&a_start[idx]) {
                        if a_cmap.get(&eptr).map_or(0, |v| v.len()) != 2 {
                            b_regular = false;
                            break;
                        }
                    }
                }
            }
        }
        result.push(ConnexityBlock {
            shapes: a_lcs,
            regular: b_regular,
        });
    }
    result
}

/// OCCT BOPAlgo_ShellSplitter::MakeShells (BOPAlgo_ShellSplitter.cxx L621-679).
/// Regular blocks become shells directly (via MakeShell, which orients the
/// faces); non-regular blocks are split via SplitBlock.
pub fn make_shells(blocks: &[ConnexityBlock]) -> Vec<Vec<Shape>> {
    let mut shells: Vec<Vec<Shape>> = Vec::new();
    for cb in blocks {
        if cb.regular {
            // OCCT L643-648: MakeShell(aLF, aShell) — adds the faces and
            // calls BOPTools_AlgoTools::OrientFacesOnShell; Closed(true).
            let mut a_shell = cb.shapes.clone();
            orient_faces_on_shell(&mut a_shell);
            shells.push(a_shell);
        } else {
            // OCCT L651-678: SplitBlock.
            let loops = split_block(&cb.shapes);
            for lp in loops {
                shells.push(lp);
            }
        }
    }
    shells
}

/// OCCT BOPTools_AlgoTools::OrientFacesOnShell (BOPTools_AlgoTools.cxx L363-503).
/// Reorients the faces of a closed shell so that every shared edge is used by
/// the two adjacent faces with opposite orientations (consistent outward
/// orientation). Seam edges (BRep_Tool::IsClosed(aE, aF)) are not flipped.
fn orient_faces_on_shell(faces: &mut Vec<Shape>) {
    // OCCT L385-387: TopExp::MapShapesAndAncestors(aShell, EDGE, FACE, aEFMap).
    let mut a_ef_map: HashMap<u64, Vec<usize>> = HashMap::new();
    let mut a_edge_map: HashMap<u64, Shape> = HashMap::new();
    for (fi, f) in faces.iter().enumerate() {
        for e in face_edges(f) {
            a_ef_map.entry(e.ptr_id()).or_default().push(fi);
            a_edge_map.entry(e.ptr_id()).or_insert_with(|| e.clone());
        }
    }
    // OCCT L381-403: dedup equivalent faces per edge. A seam edge of a
    // periodic face (cylinder/cone/sphere lateral wire stores the seam edge
    // twice) maps the same face twice; aFM (TopTools_ShapeMapHasher,
    // orientation-insensitive) keeps the first occurrence only.
    for a_lf in a_ef_map.values_mut() {
        if a_lf.len() > 1 {
            let mut a_fm: HashSet<u64> = HashSet::new();
            a_lf.retain(|&fi| a_fm.insert(faces[fi].ptr_id()));
        }
    }
    //
    // aProcessedFaces — IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>
    // (orientation-insensitive); stores the face as currently oriented.
    let mut a_processed_keys: HashSet<u64> = HashSet::new();
    // New shell contents (face indices, in the order faces get added).
    let mut a_shell_new: Vec<usize> = Vec::new();
    //
    // OCCT L410-459: process the edges with exactly 2 faces.
    let keys: Vec<u64> = a_ef_map.keys().copied().collect();
    for eptr in keys {
        let a_e = a_edge_map[&eptr].clone();
        // OCCT L417: skip degenerated edges.
        let degen = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
        if degen {
            continue;
        }
        let a_lf = &a_ef_map[&eptr];
        let a_nb_f = a_lf.len();
        if a_nb_f != 2 {
            continue;
        }
        let f1_idx = a_lf[0];
        let f2_idx = a_lf[1];
        //
        let mut b_is_processed1 = a_processed_keys.contains(&faces[f1_idx].ptr_id());
        let mut b_is_processed2 = a_processed_keys.contains(&faces[f2_idx].ptr_id());
        if b_is_processed1 && b_is_processed2 {
            continue;
        }
        //
        if !b_is_processed1 && !b_is_processed2 {
            a_processed_keys.insert(faces[f1_idx].ptr_id());
            a_shell_new.push(f1_idx);
            b_is_processed1 = true;
        }
        //
        let an_or_e1 = edge_orientation_in_face(eptr, &faces[f1_idx]);
        let an_or_e2 = edge_orientation_in_face(eptr, &faces[f2_idx]);
        //
        if b_is_processed1 && !b_is_processed2 {
            if an_or_e1 == an_or_e2 {
                if !edge_closed_on_face(&a_e, faces[f1_idx].index)
                    && !edge_closed_on_face(&a_e, faces[f2_idx].index)
                {
                    faces[f2_idx].orientation = flip_orientation(faces[f2_idx].orientation);
                }
            }
            a_processed_keys.insert(faces[f2_idx].ptr_id());
            a_shell_new.push(f2_idx);
        } else if !b_is_processed1 && b_is_processed2 {
            if an_or_e1 == an_or_e2 {
                if !edge_closed_on_face(&a_e, faces[f1_idx].index)
                    && !edge_closed_on_face(&a_e, faces[f2_idx].index)
                {
                    faces[f1_idx].orientation = flip_orientation(faces[f1_idx].orientation);
                }
            }
            a_processed_keys.insert(faces[f1_idx].ptr_id());
            a_shell_new.push(f1_idx);
        }
    }
    // OCCT L460-497: add the unprocessed faces of the other edges (free edges
    // and multi-connected edges) to the new shell.
    for eptr in a_ef_map.keys() {
        let a_e = a_edge_map[eptr].clone();
        let degen = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
        if degen {
            continue;
        }
        let a_lf = &a_ef_map[eptr];
        let a_nb_f = a_lf.len();
        if a_nb_f != 2 {
            for &fi in a_lf {
                if !a_processed_keys.contains(&faces[fi].ptr_id()) {
                    a_processed_keys.insert(faces[fi].ptr_id());
                    a_shell_new.push(fi);
                }
            }
        }
    }
    // OCCT L502: aShell = aShellNew — rebuild the shell in the new order.
    let new_faces: Vec<Shape> = a_shell_new.iter().map(|&i| faces[i].clone()).collect();
    *faces = new_faces;
}

/// OCCT BRep_Tool::IsClosed(aE, aF) — the edge has two pcurves on the closed
/// surface (seam edge). BOPAlgo_Builder_2.cxx L397. `f_index` is the face's
/// BRep-slot index (`Shape::index`), matching the `face` field stored in the
/// edge's CurveOnClosedSurface representation.
fn edge_closed_on_face(a_e: &Shape, f_index: usize) -> bool {
    a_e.as_edge()
        .map(|ed| {
            ed.representations.iter().any(|r| {
                matches!(
                    r,
                    topods::CurveRepresentation::CurveOnClosedSurface { face, .. }
                        if *face == f_index
                )
            })
        })
        .unwrap_or(false)
}

fn flip_orientation(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// Key of a Shape in an OCCT NCollection_Map<TopoDS_Shape> — includes the
/// orientation, because TopoDS_Shape::IsEqual compares myOrient
/// (TopoDS_Shape.hxx L276-280). NCollection_Map with the default hasher
/// (e.g. aMFaces, AddedFacesMap) is orientation-sensitive; maps using
/// TopTools_ShapeMapHasher (aBoundaryFaces) ignore it.
fn shape_map_key(s: &Shape) -> (u64, Orientation) {
    (s.ptr_id(), s.orientation)
}

/// OCCT BOPAlgo_ShellSplitter::SplitBlock (BOPAlgo_ShellSplitter.cxx L153-421).
fn split_block(shapes: &[Shape]) -> Vec<Vec<Shape>> {
    // OCCT L176-182: aMFaces — all faces of the block (orientation-sensitive).
    let mut a_mfaces: HashSet<(u64, Orientation)> =
        shapes.iter().map(|f| shape_map_key(f)).collect();
    // Edge -> (edge shape, [faces]) map for the faces still in aMFaces.
    let mut a_ef_map: HashMap<u64, (Shape, Vec<Shape>)> = HashMap::new();
    // OCCT L184-222: remove the faces with free edges, iteratively.
    loop {
        a_ef_map.clear();
        for f in shapes {
            if !a_mfaces.contains(&shape_map_key(f)) {
                continue;
            }
            for e in face_edges(f) {
                let eptr = e.ptr_id();
                let entry = a_ef_map.entry(eptr).or_insert_with(|| (e.clone(), Vec::new()));
                entry.1.push(f.clone());
            }
        }
        let a_nb_begin = a_mfaces.len();
        for (eptr, (a_e, a_lf)) in &a_ef_map {
            // OCCT L205: skip degenerated and INTERNAL edges.
            let degen = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
            if degen {
                continue;
            }
            if a_e.orientation == Orientation::Internal {
                continue;
            }
            if a_lf.len() == 1 {
                a_mfaces.remove(&shape_map_key(&a_lf[0]));
            }
        }
        let a_nb_end = a_mfaces.len();
        if a_nb_end == a_nb_begin || a_nb_end == 0 {
            break;
        }
    }
    if a_mfaces.is_empty() {
        return Vec::new();
    }
    // OCCT L230-245: connected faces + boundary faces (used exactly once).
    // aBoundaryFaces uses TopTools_ShapeMapHasher (orientation-insensitive).
    let mut a_lf_connected: Vec<Shape> = Vec::new();
    let mut a_boundary_faces: HashSet<u64> = HashSet::new();
    let mut a_bf_seen: HashSet<u64> = HashSet::new();
    for f in shapes {
        if !a_mfaces.contains(&shape_map_key(f)) {
            continue;
        }
        a_lf_connected.push(f.clone());
        if !a_bf_seen.insert(f.ptr_id()) {
            a_boundary_faces.remove(&f.ptr_id());
        } else {
            a_boundary_faces.insert(f.ptr_id());
        }
    }
    let a_nb_shapes = a_lf_connected.len();
    let mut b_all_faces_taken = false;
    let mut a_added: HashSet<(u64, Orientation)> = HashSet::new();
    let mut loops: Vec<Vec<Shape>> = Vec::new();
    for a_ff in &a_lf_connected {
        if b_all_faces_taken {
            break;
        }
        // OCCT L255: if (!AddedFacesMap.Add(aFF)) continue;
        if !a_added.insert(shape_map_key(a_ff)) {
            continue;
        }
        // OCCT L260-266: MakeShell; Add(aShell, aFF);
        // MapShapesAndAncestors(aShell, EDGE, FACE, aMEFP).
        let mut a_shell: Vec<Shape> = vec![a_ff.clone()];
        // Edge -> [face indices in the shell] (OCCT aMEFP).
        let mut a_mefp: HashMap<u64, Vec<usize>> = HashMap::new();
        for e in face_edges(a_ff) {
            a_mefp.entry(e.ptr_id()).or_default().push(0);
        }
        let mut i = 0;
        while i < a_shell.len() {
            let a_f = a_shell[i].clone();
            let is_boundary = a_boundary_faces.contains(&a_f.ptr_id());
            // OCCT L277-368: expand the shell along the free edges of a_f.
            for e in face_edges(&a_f) {
                let eptr = e.ptr_id();
                // OCCT L283-291: proceed only free edges in this shell.
                if let Some(users) = a_mefp.get(&eptr) {
                    if users.len() > 1 {
                        continue;
                    }
                }
                // OCCT L293-297: avoid INTERNAL edges.
                if e.orientation == Orientation::Internal {
                    continue;
                }
                // OCCT L299-302: avoid degenerated edges.
                let degen = e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                if degen {
                    continue;
                }
                // OCCT L305-310: candidate faces using this edge.
                let a_lf = match a_ef_map.get(&eptr) {
                    Some((_, lf)) => lf.clone(),
                    None => continue,
                };
                if a_lf.is_empty() {
                    continue;
                }
                // OCCT L314-341: prepare the candidate list.
                let mut a_lcs_off: Vec<Shape> = Vec::new();
                let mut a_nb_ways_inside = 0usize;
                let mut a_sel_f: Option<Shape> = None;
                for a_fl in &a_lf {
                    // OCCT L322: aF.IsSame(aFL) — same TShape, any orientation.
                    if a_f.ptr_id() == a_fl.ptr_id() {
                        continue;
                    }
                    // OCCT L322: AddedFacesMap.Contains(aFL) — orientation-sensitive.
                    if a_added.contains(&shape_map_key(a_fl)) {
                        continue;
                    }
                    // OCCT L328: GetEdgeOff(aE, aFL, aEL).
                    if !get_edge_off(&e, a_fl) {
                        continue;
                    }
                    if is_boundary && !a_boundary_faces.contains(&a_fl.ptr_id()) {
                        a_nb_ways_inside += 1;
                        a_sel_f = Some(a_fl.clone());
                    }
                    a_lcs_off.push(a_fl.clone());
                }
                let a_nb_off = a_lcs_off.len();
                if a_nb_off == 0 {
                    continue;
                }
                // OCCT L349-361: select the next face.
                if !is_boundary || a_nb_ways_inside != 1 {
                    if a_nb_off == 1 {
                        a_sel_f = Some(a_lcs_off[0].clone());
                    } else if a_nb_off > 1 {
                        a_sel_f = get_face_off(&e, &a_f, &a_lcs_off);
                    }
                }
                if let Some(a_sel) = a_sel_f {
                    // OCCT L363-367: if (!aSelF.IsNull() && AddedFacesMap.Add(aSelF)).
                    if a_added.insert(shape_map_key(&a_sel)) {
                        let new_idx = a_shell.len();
                        a_shell.push(a_sel);
                        // OCCT L366: MapShapesAndAncestors(aSelF, EDGE, FACE, aMEFP).
                        for se in face_edges(&a_shell[new_idx]) {
                            a_mefp.entry(se.ptr_id()).or_default().push(new_idx);
                        }
                    }
                }
            }
            i += 1;
        }
        // OCCT L371-392: RefineShell on multi-connected edges; collect the
        // closed sub-shells into the result, the not-closed ones for later
        // re-use of their faces.
        let a_lsh_sp = refine_shell(&a_shell, &a_mefp);
        let mut a_lsh_nc: Vec<Vec<Shape>> = Vec::new();
        for sh_sp in &a_lsh_sp {
            if is_shell_closed(sh_sp) {
                loops.push(sh_sp.clone());
            } else {
                a_lsh_nc.push(sh_sp.clone());
            }
        }
        b_all_faces_taken = a_added.len() == a_nb_shapes;
        if b_all_faces_taken {
            break;
        }
        // OCCT L400-405: a single sub-shell — not further processing of
        // not-closed shells is needed, it will not bring any new results.
        if a_lsh_sp.len() == 1 {
            continue;
        }
        // OCCT L408-419: remove the faces of the not-closed shells from the map
        // of processed faces and try to rebuild the shells using all not
        // processed faces, because faces of one shell might be needed for
        // building the other.
        for sh_nc in &a_lsh_nc {
            for f in sh_nc {
                a_added.remove(&shape_map_key(f));
            }
        }
    }
    loops
}

/// OCCT BOPTools_AlgoTools::GetEdgeOff (BOPTools_AlgoTools.cxx L1107-1135).
/// Finds the edge in the face theF2 that is the same shape as theE1 and has
/// the reversed orientation.
fn get_edge_off(a_e1: &Shape, a_f2: &Shape) -> bool {
    let a_or1 = a_e1.orientation;
    let a_or1c = reverse_orientation(a_or1);
    for e in face_edges(a_f2) {
        if e.ptr_id() == a_e1.ptr_id() {
            if e.orientation == a_or1c {
                return true;
            }
        }
    }
    false
}

/// OCCT BOPTools_AlgoTools::GetFaceOff (BOPTools_AlgoTools.cxx L994-1103) —
/// selects the candidate face with the minimal signed angle around the shared
/// edge, measured between the face bi-normal directions (aDB = aDN ^ aDTgt).
///
/// rcad: semantic translation for straight/analytic edges. aT = edge tangent
/// at the intermediate point; aP/aT only feed the MinStep3D/FindPointInFace
/// refinement of aDB, which is a no-op for planar faces and is omitted here.
fn get_face_off(a_e1: &Shape, a_f1: &Shape, candidates: &[Shape]) -> Option<Shape> {
    // OCCT L1014-1019: aT = IntermediatePoint; aPx = aC3D->D0(aT);
    // aVTgt = EdgeTangent(aE1, aT) (reversed for a REVERSED edge).
    let a_t = edge_tangent(a_e1)?;
    let a_or = a_e1.orientation;
    // OCCT GetFaceDir(aE1, aF1) L2118-2160: aDN1 = normal to the face at the
    // edge point (flipped by the face orientation); aDBF = aDN1 ^ aDTgt.
    let a_dn1 = face_normal_at(a_e1, a_f1)?;
    let a_dbf = a_dn1.cross(a_t);
    // OCCT L1038: aDTF = aDN1 ^ aDBF.
    let a_dtf = a_dn1.cross(a_dbf);
    let mut a_angle_min = 100.0f64;
    let mut a_sel: Option<Shape> = None;
    for a_f2 in candidates {
        // OCCT L1052: aDTgt2 = (aE2.Orientation() == aOr) ? aDTgt :
        // aDTgt.Reversed(). GetEdgeOff guarantees the candidate edge aE2 has
        // the reversed orientation, so the tangent is flipped.
        let a_dtgt2 = -a_t;
        // OCCT GetFaceDir(aE2, aF2): aDN2; aDBF2 = aDN2 ^ aDTgt2.
        let a_dn2 = face_normal_at(a_e1, a_f2)?;
        let a_dbf2 = a_dn2.cross(a_dtgt2);
        // OCCT L1063: aAngle = AngleWithRef(aDBF, aDBF2, aDTF).
        let mut a_angle = angle_with_ref_bop(a_dbf, a_dbf2, a_dtf);
        // OCCT L1065-1083: near-zero angle special cases set the angle to PI or
        // 2*PI when aF2 equals (or IsSame as) theF1, or when bIsComputed is
        // false. In split_block the current face is filtered out as a candidate,
        // so aF2 never equals/IsSame theF1, and the planar/analytic faces
        // supported by face_normal_at always compute the bi-normal (bIsComputed
        // true). Thus the generic branch below applies in practice.
        // OCCT L1091-1094: map the negative angles into the [0, 2*PI) range.
        if a_angle < 0.0 {
            a_angle = std::f64::consts::TAU + a_angle;
        }
        // OCCT L1096-1100: keep the minimal angle.
        if a_angle < a_angle_min {
            a_angle_min = a_angle;
            a_sel = Some(a_f2.clone());
        }
    }
    a_sel
}

/// OCCT BOPTools_AlgoTools::AngleWithRef (BOPTools_AlgoTools.cxx L1946-1975).
/// Signed angle aBeta = PI/2*(1 - d1.d2), sign from (d1 x d2).dRef.
/// (The |d1 x d2| >= 0 branch is always taken, so the else-branch is dead.)
fn angle_with_ref_bop(d1: DVec3, d2: DVec3, d_ref: DVec3) -> f64 {
    let a_xyz = d1.cross(d2);
    let a_cosinus = d1.dot(d2);
    let mut a_beta = 0.5 * std::f64::consts::PI * (1.0 - a_cosinus);
    if a_xyz.dot(d_ref) < 0.0 {
        a_beta = -a_beta;
    }
    a_beta
}

/// OCCT BOPAlgo_ShellSplitter::RefineShell (BOPAlgo_ShellSplitter.cxx L443-617).
/// Splits a shell on the edges shared by more than two faces (branch edges).
fn refine_shell(a_shell: &[Shape], a_mef: &HashMap<u64, Vec<usize>>) -> Vec<Vec<Shape>> {
    if a_shell.is_empty() {
        return Vec::new();
    }
    // OCCT L455-512: find the branch edges (edges with >2 adjacent faces, or
    // with 2 faces traversing the edge with the same orientation, or internal
    // edges counted twice).
    let mut a_me_stop: HashSet<u64> = HashSet::new();
    for (eptr, a_lf) in a_mef {
        if a_lf.len() > 2 {
            a_me_stop.insert(*eptr);
            continue;
        }
        if a_lf.len() == 2 {
            let f1 = &a_shell[a_lf[0]];
            let f2 = &a_shell[a_lf[1]];
            // OCCT L475-481: FindShape(aE, aF1/2) — same edge orientation check.
            let e1_or = edge_orientation_in_face(*eptr, f1);
            let e2_or = edge_orientation_in_face(*eptr, f2);
            if let (Some(o1), Some(o2)) = (e1_or, e2_or) {
                if o1 == o2 {
                    a_me_stop.insert(*eptr);
                    continue;
                }
            }
        }
        // OCCT L486-511: count faces, counting INTERNAL edges twice.
        let mut a_nb_f = 0usize;
        for &fi in a_lf {
            a_nb_f += 1;
            if edge_internal_in_face(*eptr, &a_shell[fi]) {
                a_nb_f += 1;
            }
            if a_nb_f > 2 {
                break;
            }
        }
        if a_nb_f > 2 {
            a_me_stop.insert(*eptr);
        }
    }
    if a_me_stop.is_empty() {
        // OCCT L514-518: no branch edges — the whole shell is the result.
        return vec![a_shell.to_vec()];
    }
    // OCCT L520-617: split the shell into sub-shells grown from each face
    // without crossing the branch edges.
    let mut a_mf_processed: HashSet<u64> = HashSet::new();
    let mut shells: Vec<Vec<Shape>> = Vec::new();
    for a_f1 in a_shell {
        if !a_mf_processed.insert(a_f1.ptr_id()) {
            continue;
        }
        // OCCT L536-601: BFS from aF1 avoiding the branch edges.
        let mut a_mfb: Vec<Shape> = vec![a_f1.clone()];
        let mut a_lfp: Vec<Shape> = vec![a_f1.clone()];
        loop {
            let mut a_lfp1: Vec<Shape> = Vec::new();
            for a_fp in &a_lfp {
                for e in face_edges(a_fp) {
                    let eptr = e.ptr_id();
                    if a_me_stop.contains(&eptr) {
                        continue;
                    }
                    if e.orientation == Orientation::Internal {
                        continue;
                    }
                    let degen = e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                    if degen {
                        continue;
                    }
                    let a_lf = match a_mef.get(&eptr) {
                        Some(lf) => lf.clone(),
                        None => continue,
                    };
                    for &fi in &a_lf {
                        let a_fp1 = &a_shell[fi];
                        if a_fp1.ptr_id() == a_fp.ptr_id() {
                            continue;
                        }
                        if a_mfb.iter().any(|f| f.ptr_id() == a_fp1.ptr_id()) {
                            continue;
                        }
                        if a_mf_processed.insert(a_fp1.ptr_id()) {
                            a_mfb.push(a_fp1.clone());
                            a_lfp1.push(a_fp1.clone());
                        }
                    }
                }
            }
            if a_lfp1.is_empty() {
                break;
            }
            a_lfp = a_lfp1;
        }
        if !a_mfb.is_empty() {
            shells.push(a_mfb);
        }
    }
    shells
}

/// A shell is closed when every non-degenerate, non-INTERNAL/EXTERNAL edge is
/// used by an even number of faces (OCCT BRep_Tool::IsClosed, BRep_Tool.cxx
/// L1707-1728: toggle each edge in a map; closed iff hasBound && map empty).
fn is_shell_closed(faces: &[Shape]) -> bool {
    let mut a_map: HashSet<u64> = HashSet::new();
    let mut has_bound = false;
    for f in faces {
        for e in face_edges(f) {
            let degen = e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
            if degen
                || e.orientation == Orientation::Internal
                || e.orientation == Orientation::External
            {
                continue;
            }
            has_bound = true;
            if !a_map.insert(e.ptr_id()) {
                a_map.remove(&e.ptr_id());
            }
        }
    }
    has_bound && a_map.is_empty()
}

/// Orientation of the edge shape in a face (OCCT FindShape).
fn edge_orientation_in_face(eptr: u64, face: &Shape) -> Option<Orientation> {
    for e in face_edges(face) {
        if e.ptr_id() == eptr {
            return Some(e.orientation);
        }
    }
    None
}

/// True when the edge appears in the face with INTERNAL orientation.
fn edge_internal_in_face(eptr: u64, face: &Shape) -> bool {
    for e in face_edges(face) {
        if e.ptr_id() == eptr {
            return e.orientation == Orientation::Internal;
        }
    }
    false
}

fn reverse_orientation(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// Midpoint of an edge from its two vertices.
fn edge_midpoint(e: &Shape) -> Option<DVec3> {
    match &*e.data {
        TShape::Edge(ed) => {
            let p1 = vertex_point(&ed.first)?;
            let p2 = vertex_point(&ed.last)?;
            Some((p1 + p2) * 0.5)
        }
        _ => None,
    }
}

/// Edge tangent at the edge parameterization direction (straight edges).
/// OCCT BOPTools_AlgoTools2D::EdgeTangent (BOPTools_AlgoTools2D.cxx L74-103):
/// the curve derivative at aT, normalized, reversed for a REVERSED edge.
fn edge_tangent(e: &Shape) -> Option<DVec3> {
    match &*e.data {
        TShape::Edge(ed) => {
            let d = vertex_point(&ed.last)? - vertex_point(&ed.first)?;
            if d.length_squared() < 1e-24 {
                None
            } else {
                let mut t = d.normalize();
                if e.orientation == Orientation::Reversed {
                    t = -t;
                }
                Some(t)
            }
        }
        _ => None,
    }
}

fn vertex_point(v: &Shape) -> Option<DVec3> {
    match &*v.data {
        TShape::Vertex(vd) => Some(vd.point),
        _ => None,
    }
}

/// Outward normal of a face at the shared edge's midpoint.
/// OCCT BOPTools_AlgoTools3D::GetNormalToFaceOnEdge + GetFaceDir L2133-2137:
/// the surface normal at the edge point, flipped by the face orientation.
fn face_normal_at(a_e: &Shape, a_f: &Shape) -> Option<DVec3> {
    use rcad_kernel::geom::SurfaceEval;
    let p = edge_midpoint(a_e)?;
    match &*a_f.data {
        TShape::Face(fd) => {
            let surf = fd.surface.as_ref()?;
            let n = match surf {
                rcad_kernel::geom::Surface3::Plane(pl) => pl.normal,
                _ => {
                    // Project the edge point onto the surface and evaluate the
                    // normal there (analytic surfaces dispatch in
                    // closest_point_on_surface).
                    let (uv, _) = crate::bop::closest_point_on_surface(surf, p);
                    surf.normal_at(uv.x, uv.y)
                }
            };
            Some(if a_f.orientation == Orientation::Reversed { -n } else { n })
        }
        _ => None,
    }
}

/// Extract edge ptr_ids from a Face Shape (outer + inner wires).
fn face_edge_ptrs(face: &Shape) -> Vec<u64> {
    face_edges(face).iter().map(|e| e.ptr_id()).collect()
}

/// Extract edge Shapes from a Face Shape (outer + inner wires), composing the
/// face and wire orientations into each edge (OCCT TopExp_Explorer composes
/// the parent orientation at every level: TopExp_Explorer.cxx L152, L110-170;
/// TopoDS_Iterator.cxx L35-37, L72-80).
fn face_edges(face: &Shape) -> Vec<Shape> {
    let mut edges = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                let w_or = fd.outer_wire.orientation;
                for e in &wd.edges {
                    let mut e2 = e.clone();
                    e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
                    edges.push(e2);
                }
            }
            for iw in &fd.inner_wires {
                if let TShape::Wire(wd) = &*iw.data {
                    let w_or = iw.orientation;
                    for e in &wd.edges {
                        let mut e2 = e.clone();
                        e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
                        edges.push(e2);
                    }
                }
            }
        }
        _ => {}
    }
    edges
}
