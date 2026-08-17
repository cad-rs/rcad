// OCCT BOPAlgo_ShellSplitter — shell partitioning.
//
// OCCT BOPAlgo_ShellSplitter.cxx
// Performs: MakeConnexityBlocks -> MakeShells (regular blocks directly,
// non-regular blocks via SplitBlock + RefineShell).

use crate::bop::ds::DS;
use indexmap::IndexMap;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, Orientation, TShape};
use std::collections::{HashMap, HashSet, VecDeque};

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
    // TopTools_ShapeMapHasher (TShape + Location).
    let mut a_mfence: HashSet<(u64, u32)> = HashSet::new();
    let mut a_mn_regular: HashSet<(u64, u32)> = HashSet::new();
    let mut a_start: Vec<Shape> = Vec::new();
    for s in shapes {
        let key = (s.ptr_id(), s.location);
        if !a_mfence.insert(key) {
            a_mn_regular.insert(key);
        } else {
            a_start.push(s.clone());
        }
    }
    // Edge -> [face indices] map (OCCT TopExp::MapShapesAndAncestors,
    // TopTools_ShapeMapHasher).
    let mut a_cmap: HashMap<(u64, u32), Vec<usize>> = HashMap::new();
    for (i, s) in a_start.iter().enumerate() {
        for ekey in face_edge_keys(s) {
            a_cmap.entry(ekey).or_default().push(i);
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
            for ekey in face_edge_keys(&a_start[i]) {
                if let Some(nei) = a_cmap.get(&ekey) {
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
            if a_mn_regular.contains(&(a_s.ptr_id(), a_s.location)) {
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
                    for ekey in face_edge_keys(&a_start[idx]) {
                        if a_cmap.get(&ekey).map_or(0, |v| v.len()) != 2 {
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
pub fn make_shells(blocks: &[ConnexityBlock], ds: &DS) -> Vec<Vec<Shape>> {
    let mut shells: Vec<Vec<Shape>> = Vec::new();
    for cb in blocks {
        if std::env::var("RCAD_BS_DEBUG").is_ok() {
            let desc: Vec<String> = cb.shapes.iter().map(|f| {
                let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                    rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                    _ => "O".into(),
                }).unwrap_or_else(|| "?".into());
                format!("{}:{}:{}", f.ptr_id(), n, if f.orientation == Orientation::Reversed { "R" } else { "F" })
            }).collect();
            eprintln!("[BS-BLOCK] regular={} n={} faces=[{}]", cb.regular, cb.shapes.len(), desc.join(" "));
        }
        if cb.regular {
            // OCCT L643-648: MakeShell(aLF, aShell) — adds the faces and
            // calls BOPTools_AlgoTools::OrientFacesOnShell; Closed(true).
            let mut a_shell = cb.shapes.clone();
            orient_faces_on_shell(&mut a_shell);
            shells.push(a_shell);
        } else {
            // OCCT L651-678: SplitBlock.
            let loops = split_block(&cb.shapes, ds);
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
pub(crate) fn orient_faces_on_shell(faces: &mut Vec<Shape>) {
    // OCCT L385-387: TopExp::MapShapesAndAncestors(aShell, EDGE, FACE, aEFMap).
    // OCCT aEFMap is NCollection_IndexedDataMap with TopTools_ShapeMapHasher —
    // iteration is in insertion order; IndexMap reproduces that order (a
    // HashMap would be random).
    let mut a_ef_map: IndexMap<(u64, u32), Vec<usize>> = IndexMap::new();
    let mut a_edge_map: IndexMap<(u64, u32), Shape> = IndexMap::new();
    for (fi, f) in faces.iter().enumerate() {
        for e in face_edges(f) {
            let ekey = (e.ptr_id(), e.location);
            a_ef_map.entry(ekey).or_default().push(fi);
            a_edge_map.entry(ekey).or_insert_with(|| e.clone());
        }
    }
    // OCCT L381-403: dedup equivalent faces per edge. A seam edge of a
    // periodic face (cylinder/cone/sphere lateral wire stores the seam edge
    // twice) maps the same face twice; aFM (TopTools_ShapeMapHasher,
    // orientation-insensitive) keeps the first occurrence only.
    for a_lf in a_ef_map.values_mut() {
        if a_lf.len() > 1 {
            let mut a_fm: HashSet<(u64, u32)> = HashSet::new();
            a_lf.retain(|&fi| a_fm.insert((faces[fi].ptr_id(), faces[fi].location)));
        }
    }
    //
    // aProcessedFaces — IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>
    // (orientation-insensitive); stores the face as currently oriented.
    let mut a_processed_keys: HashSet<(u64, u32)> = HashSet::new();
    // New shell contents (face indices, in the order faces get added).
    let mut a_shell_new: Vec<usize> = Vec::new();
    //
    // OCCT L410-459: process the edges with exactly 2 faces.
    let keys: Vec<(u64, u32)> = a_ef_map.keys().copied().collect();
    for ekey in keys {
        let a_e = a_edge_map[&ekey].clone();
        // OCCT L417: skip degenerated edges.
        let degen = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
        if degen {
            continue;
        }
        let a_lf = &a_ef_map[&ekey];
        let a_nb_f = a_lf.len();
        if a_nb_f != 2 {
            continue;
        }
        let f1_idx = a_lf[0];
        let f2_idx = a_lf[1];
        //
        let mut b_is_processed1 = a_processed_keys.contains(&(faces[f1_idx].ptr_id(), faces[f1_idx].location));
        let mut b_is_processed2 = a_processed_keys.contains(&(faces[f2_idx].ptr_id(), faces[f2_idx].location));
        if b_is_processed1 && b_is_processed2 {
            continue;
        }
        //
        if !b_is_processed1 && !b_is_processed2 {
            a_processed_keys.insert((faces[f1_idx].ptr_id(), faces[f1_idx].location));
            a_shell_new.push(f1_idx);
            b_is_processed1 = true;
        }
        //
        let an_or_e1 = edge_orientation_in_face(ekey, &faces[f1_idx]);
        let an_or_e2 = edge_orientation_in_face(ekey, &faces[f2_idx]);
        //
        if b_is_processed1 && !b_is_processed2 {
            if an_or_e1 == an_or_e2 {
                if !edge_closed_on_face(&a_e, &faces[f1_idx])
                    && !edge_closed_on_face(&a_e, &faces[f2_idx])
                {
                    faces[f2_idx].orientation = flip_orientation(faces[f2_idx].orientation);
                }
            }
            a_processed_keys.insert((faces[f2_idx].ptr_id(), faces[f2_idx].location));
            a_shell_new.push(f2_idx);
        } else if !b_is_processed1 && b_is_processed2 {
            if an_or_e1 == an_or_e2 {
                if !edge_closed_on_face(&a_e, &faces[f1_idx])
                    && !edge_closed_on_face(&a_e, &faces[f2_idx])
                {
                    faces[f1_idx].orientation = flip_orientation(faces[f1_idx].orientation);
                }
            }
            a_processed_keys.insert((faces[f1_idx].ptr_id(), faces[f1_idx].location));
            a_shell_new.push(f1_idx);
        }
    }
    // OCCT L460-497: add the unprocessed faces of the other edges (free edges
    // and multi-connected edges) to the new shell.
    for ekey in a_ef_map.keys() {
        let a_e = a_edge_map[ekey].clone();
        let degen = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
        if degen {
            continue;
        }
        let a_lf = &a_ef_map[ekey];
        let a_nb_f = a_lf.len();
        if a_nb_f != 2 {
            for &fi in a_lf {
                if !a_processed_keys.contains(&(faces[fi].ptr_id(), faces[fi].location)) {
                    a_processed_keys.insert((faces[fi].ptr_id(), faces[fi].location));
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
fn edge_closed_on_face(a_e: &Shape, a_f: &Shape) -> bool {
    let f_key = (a_f.ptr_id(), a_f.location);
    a_e.as_edge()
        .map(|ed| {
            ed.representations.iter().any(|r| {
                matches!(
                    r,
                    topods::CurveRepresentation::CurveOnClosedSurface { face, .. }
                        if *face == f_key
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
fn shape_map_key(s: &Shape) -> (u64, u32, Orientation) {
    (s.ptr_id(), s.location, s.orientation)
}

/// OCCT BOPAlgo_ShellSplitter::SplitBlock (BOPAlgo_ShellSplitter.cxx L153-421).
fn split_block(shapes: &[Shape], ds: &DS) -> Vec<Vec<Shape>> {
    // OCCT L176-182: aMFaces — all faces of the block (orientation-sensitive).
    let mut a_mfaces: HashSet<(u64, u32, Orientation)> =
        shapes.iter().map(|f| shape_map_key(f)).collect();
    // Edge -> (edge shape, [faces]) map for the faces still in aMFaces.
    // OCCT aEFMap is NCollection_IndexedDataMap with TopTools_ShapeMapHasher —
    // insertion order; IndexMap reproduces it (a HashMap would iterate in random order).
    let mut a_ef_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
    // OCCT L184-222: remove the faces with free edges, iteratively.
    loop {
        a_ef_map.clear();
        for f in shapes {
            if !a_mfaces.contains(&shape_map_key(f)) {
                continue;
            }
            for e in face_edges(f) {
                let ekey = (e.ptr_id(), e.location);
                let entry = a_ef_map.entry(ekey).or_insert_with(|| (e.clone(), Vec::new()));
                entry.1.push(f.clone());
            }
        }
        let a_nb_begin = a_mfaces.len();
        for (ekey, (a_e, a_lf)) in &a_ef_map {
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
    let mut a_boundary_faces: HashSet<(u64, u32)> = HashSet::new();
    let mut a_bf_seen: HashSet<(u64, u32)> = HashSet::new();
    for f in shapes {
        if !a_mfaces.contains(&shape_map_key(f)) {
            continue;
        }
        a_lf_connected.push(f.clone());
        let fkey = (f.ptr_id(), f.location);
        if !a_bf_seen.insert(fkey) {
            a_boundary_faces.remove(&fkey);
        } else {
            a_boundary_faces.insert(fkey);
        }
    }
    let a_nb_shapes = a_lf_connected.len();
    let mut b_all_faces_taken = false;
    let mut a_added: HashSet<(u64, u32, Orientation)> = HashSet::new();
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
        // Edge -> [face indices in the shell] (OCCT aMEFP — IndexedDataMap
        // with TopTools_ShapeMapHasher).
        let mut a_mefp: IndexMap<(u64, u32), Vec<usize>> = IndexMap::new();
        for e in face_edges(a_ff) {
            a_mefp.entry((e.ptr_id(), e.location)).or_default().push(0);
        }
        let mut i = 0;
        while i < a_shell.len() {
            let a_f = a_shell[i].clone();
            let is_boundary = a_boundary_faces.contains(&(a_f.ptr_id(), a_f.location));
            // OCCT L277-368: expand the shell along the free edges of a_f.
            for e in face_edges(&a_f) {
                let ekey = (e.ptr_id(), e.location);
                // OCCT L283-291: proceed only free edges in this shell.
                if let Some(users) = a_mefp.get(&ekey) {
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
                let a_lf = match a_ef_map.get(&ekey) {
                    Some((_, lf)) => lf.clone(),
                    None => continue,
                };
                if a_lf.is_empty() {
                    continue;
                }
                // OCCT L314-341: prepare the candidate list. Each candidate is
                // the couple (aEL, aFL): aEL is the edge in aFL that is the same
                // shape as aE with the reversed orientation (GetEdgeOff).
                let mut a_lcs_off: Vec<(Shape, Shape)> = Vec::new();
                let mut a_nb_ways_inside = 0usize;
                let mut a_sel_f: Option<Shape> = None;
                for a_fl in &a_lf {
                    // OCCT L322: aF.IsSame(aFL) — same TShape, any orientation.
                    if a_f.ptr_id() == a_fl.ptr_id() && a_f.location == a_fl.location {
                        continue;
                    }
                    // OCCT L322: AddedFacesMap.Contains(aFL) — orientation-sensitive.
                    if a_added.contains(&shape_map_key(a_fl)) {
                        continue;
                    }
                    // OCCT L328: GetEdgeOff(aE, aFL, aEL).
                    let a_el = match get_edge_off(&e, a_fl) {
                        Some(el) => el,
                        None => continue,
                    };
                    if is_boundary && !a_boundary_faces.contains(&(a_fl.ptr_id(), a_fl.location)) {
                        a_nb_ways_inside += 1;
                        a_sel_f = Some(a_fl.clone());
                    }
                    a_lcs_off.push((a_el, a_fl.clone()));
                }
                let a_nb_off = a_lcs_off.len();
                if std::env::var("RCAD_BS_DEBUG").is_ok() {
                    let cur = format!("{}", a_f.ptr_id() % 100000);
                    let cands: Vec<String> = a_lcs_off.iter().map(|(el, fl)| format!("{}:{}", fl.ptr_id() % 100000, if el.orientation == Orientation::Reversed { "R" } else { "F" })).collect();
                    eprintln!("[BS-GROW-EDGE] shell_start={} edge={} or={} is_boundary={} n_off={} n_inside={} sel={} cands=[{}]", cur, ekey.0 % 100000, if e.orientation == Orientation::Reversed { "R" } else { "F" }, is_boundary, a_nb_off, a_nb_ways_inside, a_sel_f.as_ref().map(|s| format!("{}", s.ptr_id() % 100000)).unwrap_or_else(|| "-".into()), cands.join(" "));
                    if a_nb_off == 0 {
                        let why: Vec<String> = a_lf.iter().map(|fl| {
                            let same = a_f.ptr_id() == fl.ptr_id() && a_f.location == fl.location;
                            let added = a_added.contains(&shape_map_key(fl));
                            let off = get_edge_off(&e, fl).map(|el| format!("{}", if el.orientation == Orientation::Reversed { "R" } else { "F" })).unwrap_or_else(|| "no-off".into());
                            format!("{}:same={}:added={}:off={}", fl.ptr_id() % 100000, same, added, off)
                        }).collect();
                        eprintln!("[BS-GROW-NONE] edge={} raw=[{}]", ekey.0 % 100000, why.join(" "));
                    }
                }
                if a_nb_off == 0 {
                    continue;
                }
                // OCCT L349-361: select the next face.
                if !is_boundary || a_nb_ways_inside != 1 {
                    if a_nb_off == 1 {
                        a_sel_f = Some(a_lcs_off[0].1.clone());
                    } else if a_nb_off > 1 {
                        // OCCT L359: GetFaceOff(aE, aF, aLCSOff, aSelF, aContext).
                        a_sel_f = super::builder::get_face_off(&e, &a_f, &a_lcs_off, ds).0;
                    }
                }
                if let Some(a_sel) = a_sel_f {
                    // OCCT L363-367: if (!aSelF.IsNull() && AddedFacesMap.Add(aSelF)).
                    if a_added.insert(shape_map_key(&a_sel)) {
                        let new_idx = a_shell.len();
                        a_shell.push(a_sel);
                        // OCCT L366: MapShapesAndAncestors(aSelF, EDGE, FACE, aMEFP).
                        for se in face_edges(&a_shell[new_idx]) {
                            a_mefp.entry((se.ptr_id(), se.location)).or_default().push(new_idx);
                        }
                    }
                }
            }
            i += 1;
        }
        // OCCT L371-392: RefineShell on multi-connected edges; collect the
        // closed sub-shells into the result, the not-closed ones for later
        // re-use of their faces.
        if std::env::var("RCAD_BS_DEBUG").is_ok() {
            let desc: Vec<String> = a_shell.iter().map(|f| {
                let n = f.as_face().and_then(|fd| fd.surface.clone()).map(|s| match s {
                    rcad_kernel::geom::Surface3::Plane(p) => format!("({:.2},{:.2},{:.2})", p.normal.x, p.normal.y, p.normal.z),
                    _ => "O".into(),
                }).unwrap_or_else(|| "?".into());
                format!("{}:{}", f.ptr_id(), n)
            }).collect();
            eprintln!("[BS-GROW] start={} grown={} faces=[{}]", a_ff.ptr_id(), a_shell.len(), desc.join(" "));
        }
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
fn get_edge_off(a_e1: &Shape, a_f2: &Shape) -> Option<Shape> {
    let a_or1 = a_e1.orientation;
    let a_or1c = reverse_orientation(a_or1);
    for e in face_edges(a_f2) {
        // OCCT L1123: aEF2.IsSame(theE1) — TShape + Location.
        if e.ptr_id() == a_e1.ptr_id() && e.location == a_e1.location {
            if e.orientation == a_or1c {
                return Some(e);
            }
        }
    }
    None
}

/// OCCT BOPAlgo_ShellSplitter::RefineShell (BOPAlgo_ShellSplitter.cxx L443-617).
/// Splits a shell on the edges shared by more than two faces (branch edges).
fn refine_shell(a_shell: &[Shape], a_mef: &IndexMap<(u64, u32), Vec<usize>>) -> Vec<Vec<Shape>> {
    if a_shell.is_empty() {
        return Vec::new();
    }
    // OCCT L455-512: find the branch edges (edges with >2 adjacent faces, or
    // with 2 faces traversing the edge with the same orientation, or internal
    // edges counted twice).
    let mut a_me_stop: HashSet<(u64, u32)> = HashSet::new();
    for (ekey, a_lf) in a_mef {
        if a_lf.len() > 2 {
            a_me_stop.insert(*ekey);
            continue;
        }
        if a_lf.len() == 2 {
            let f1 = &a_shell[a_lf[0]];
            let f2 = &a_shell[a_lf[1]];
            // OCCT L475-481: FindShape(aE, aF1/2) — same edge orientation check.
            let e1_or = edge_orientation_in_face(*ekey, f1);
            let e2_or = edge_orientation_in_face(*ekey, f2);
            if let (Some(o1), Some(o2)) = (e1_or, e2_or) {
                if o1 == o2 {
                    a_me_stop.insert(*ekey);
                    continue;
                }
            }
        }
        // OCCT L486-511: count faces, counting INTERNAL edges twice.
        let mut a_nb_f = 0usize;
        for &fi in a_lf {
            a_nb_f += 1;
            if edge_internal_in_face(*ekey, &a_shell[fi]) {
                a_nb_f += 1;
            }
            if a_nb_f > 2 {
                break;
            }
        }
        if a_nb_f > 2 {
            a_me_stop.insert(*ekey);
        }
    }
    if a_me_stop.is_empty() {
        // OCCT L514-518: no branch edges — the whole shell is the result.
        return vec![a_shell.to_vec()];
    }
    // OCCT L520-617: split the shell into sub-shells grown from each face
    // without crossing the branch edges. aMFProcessed is NCollection_Map
    // (orientation-sensitive), aMFB is an IndexedMap with
    // TopTools_ShapeMapHasher (orientation-insensitive).
    let mut a_mf_processed: HashSet<(u64, u32, Orientation)> = HashSet::new();
    let mut shells: Vec<Vec<Shape>> = Vec::new();
    for a_f1 in a_shell {
        if !a_mf_processed.insert(shape_map_key(a_f1)) {
            continue;
        }
        // OCCT L536-601: BFS from aF1 avoiding the branch edges.
        // OCCT aLFP1 (L529) is declared once; aLFP.Append(aLFP1) (L604-605)
        // CLEARS aLFP1 (NCollection_List::Append(List&) moves the nodes), so
        // aLFP takes exactly the faces discovered in the previous iteration —
        // a layer-by-layer BFS, reproduced by a per-iteration a_lfp1.
        let mut a_mfb: Vec<Shape> = vec![a_f1.clone()];
        let mut a_lfp: Vec<Shape> = vec![a_f1.clone()];
        loop {
            let mut a_lfp1: Vec<Shape> = Vec::new();
            for a_fp in &a_lfp {
                for e in face_edges(a_fp) {
                    let ekey = (e.ptr_id(), e.location);
                    if a_me_stop.contains(&ekey) {
                        continue;
                    }
                    if e.orientation == Orientation::Internal {
                        continue;
                    }
                    let degen = e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                    if degen {
                        continue;
                    }
                    let a_lf = match a_mef.get(&ekey) {
                        Some(lf) => lf.clone(),
                        None => continue,
                    };
                    for &fi in &a_lf {
                        let a_fp1 = &a_shell[fi];
                        if a_fp1.ptr_id() == a_fp.ptr_id() && a_fp1.location == a_fp.location {
                            continue;
                        }
                        if a_mfb.iter().any(|f| f.ptr_id() == a_fp1.ptr_id() && f.location == a_fp1.location) {
                            continue;
                        }
                        if a_mf_processed.insert(shape_map_key(a_fp1)) {
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
    // OCCT BRep_Tool::IsClosed(Shell) (BRep_Tool.cxx L1707-1728) toggles each
    // edge in an orientation-insensitive map; rcad keys by (TShape, Location).
    let mut a_map: HashSet<(u64, u32)> = HashSet::new();
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
            let ekey = (e.ptr_id(), e.location);
            if !a_map.insert(ekey) {
                a_map.remove(&ekey);
            }
        }
    }
    has_bound && a_map.is_empty()
}

/// Orientation of the edge shape in a face (OCCT FindShape).
fn edge_orientation_in_face(ekey: (u64, u32), face: &Shape) -> Option<Orientation> {
    for e in face_edges(face) {
        if e.ptr_id() == ekey.0 && e.location == ekey.1 {
            return Some(e.orientation);
        }
    }
    None
}

/// True when the edge appears in the face with INTERNAL orientation.
fn edge_internal_in_face(ekey: (u64, u32), face: &Shape) -> bool {
    for e in face_edges(face) {
        if e.ptr_id() == ekey.0 && e.location == ekey.1 {
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

/// Extract edge (ptr_id, location) pairs from a Face Shape (outer + inner wires).
/// OCCT TopExp::MapShapesAndAncestors uses TopTools_ShapeMapHasher — key
/// identity TShape + Location, orientation ignored.
fn face_edge_keys(face: &Shape) -> Vec<(u64, u32)> {
    face_edges(face).iter().map(|e| (e.ptr_id(), e.location)).collect()
}

/// Extract edge Shapes from a Face Shape (outer + inner wires), composing the
/// face and wire orientations into each edge (OCCT TopExp_Explorer composes
/// the parent orientation at every level: TopExp_Explorer.cxx L152, L110-170;
/// TopoDS_Iterator.cxx L35-37, L72-80). The face Location is composed into the
/// edges too (TopoDS_Iterator cumLoc, L76-78): an edge of a located face (the
/// revolve end cap at the rotation Location L1) is keyed by its effective
/// location, matching the located boundary edges of the neighboring faces.
pub(crate) fn face_edges(face: &Shape) -> Vec<Shape> {
    let mut edges = Vec::new();
    match &*face.data {
        TShape::Face(fd) => {            let face_loc = face.location;
            let compose = |eloc: u32| -> u32 {
                if face_loc == 0 {
                    eloc
                } else if eloc == 0 {
                    face_loc
                } else {
                    // Both non-identity: the composed transform is one of the
                    // two operands in the single-fold cases; fall back to the
                    // face location index.
                    face_loc
                }
            };
            if let TShape::Wire(wd) = &*fd.outer_wire.data {
                let w_or = fd.outer_wire.orientation;
                for e in &wd.edges {
                    let mut e2 = e.clone();
                    e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
                    e2.location = compose(e2.location);
                    edges.push(e2);
                }
            }
            for iw in &fd.inner_wires {
                if let TShape::Wire(wd) = &*iw.data {
                    let w_or = iw.orientation;
                    for e in &wd.edges {
                        let mut e2 = e.clone();
                        e2.orientation = face.orientation.compose(w_or).compose(e.orientation);
                        e2.location = compose(e2.location);
                        edges.push(e2);
                    }
                }
            }
        }
        _ => {}
    }
    edges
}
