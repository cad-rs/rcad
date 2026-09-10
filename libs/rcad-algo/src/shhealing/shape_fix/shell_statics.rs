//! 1:1 translation of OCCT `ShapeFix_Shell` — the file statics
//! (`ShapeFix_Shell.cxx` L178-1421): `GetFreeEdges` (L178-197),
//! `GetConnectedFaceGroups` (L206-297), `GetShells` (L308-648),
//! `AddMultiConexityFaces` (L656-913), `BoxIn` (L919-942),
//! `GetClosedShells` (L950-986), `GlueClosedCandidate` (L992-1164),
//! `CreateNonManifoldShells` (L1171-1320), `CreateClosedShell` (L1327-1421).
//!
//! GAP carrier (the iron rule): [`brep_bnd_lib_add_close_gap`] — OCCT
//! `BRepBndLib::AddClose` (TKTopAlgo/BRepBndLib.cxx) is not translated; the
//! box stays void, so `BoxIn` reports no nesting and `GetClosedShells` keeps
//! all candidate shells (the OCCT "no candidate nested" path).
//!
//! Bridges:
//! - `FaceEdgesMap` (IndexedDataMap<Face, Array1<Edge>>) ->
//!   [`IndexedShapeListMap`] (the ordered (key, list) map);
//! - `EdgeFacesMap` / EdgeOrientedMap / TempProcessedEdges
//!   (DataMap<Edge, ...>) -> `HashMap<(u64, u32), ...>`;
//! - `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` ->
//!   [`ShapeSet`] (the ordered shape set over the IsSame identity);
//! - the DFS `std::stack` -> `Vec`.

use std::collections::{HashMap, HashSet};

use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType, TShape};

use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes, topexp_explorer};

/// OCCT `NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>` — the
/// ordered shape set over the IsSame identity (the map stores the shapes).
#[derive(Default)]
pub(crate) struct ShapeSet {
    items: Vec<Shape>,
}

impl ShapeSet {
    pub fn new() -> Self {
        ShapeSet { items: Vec::new() }
    }

    fn key(s: &Shape) -> (u64, u32) {
        (s.ptr_id(), s.location)
    }

    /// OCCT Add — True when the shape was not yet in the set.
    pub fn add(&mut self, s: &Shape) -> bool {
        if self.contains(Self::key(s)) {
            false
        } else {
            self.items.push(s.clone());
            true
        }
    }

    /// OCCT Remove.
    pub fn remove(&mut self, s: &Shape) {
        let k = Self::key(s);
        self.items.retain(|it| Self::key(it) != k);
    }

    /// OCCT Contains.
    pub fn contains(&self, k: (u64, u32)) -> bool {
        self.items.iter().any(|it| Self::key(it) == k)
    }

    pub fn contains_shape(&self, s: &Shape) -> bool {
        self.contains(Self::key(s))
    }

    pub fn extent(&self) -> usize {
        self.items.len()
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    pub fn clear(&mut self) {
        self.items.clear();
    }

    pub fn iter(&self) -> impl Iterator<Item = &Shape> {
        self.items.iter()
    }
}

/// OCCT `FaceEdgesMap = NCollection_IndexedDataMap<TopoDS_Face,
/// NCollection_Array1<TopoDS_Edge>>` — the ordered (face, edges) map.
pub(crate) struct IndexedShapeListMap {
    keys: Vec<Shape>,
    values: Vec<Vec<Shape>>,
    index: HashMap<(u64, u32), usize>,
}

impl IndexedShapeListMap {
    pub fn new() -> Self {
        IndexedShapeListMap {
            keys: Vec::new(),
            values: Vec::new(),
            index: HashMap::new(),
        }
    }

    /// OCCT Add(key, value) — binds (the key must be fresh, as in OCCT).
    pub fn add(&mut self, key: Shape, value: Vec<Shape>) {
        self.index.insert((key.ptr_id(), key.location), self.keys.len());
        self.keys.push(key);
        self.values.push(value);
    }

    pub fn extent(&self) -> usize {
        self.keys.len()
    }

    /// OCCT FindKey(i).
    pub fn find_key(&self, i: usize) -> &Shape {
        &self.keys[i - 1]
    }

    /// OCCT FindFromIndex(i).
    pub fn find_from_index(&self, i: usize) -> &Vec<Shape> {
        &self.values[i - 1]
    }

    /// OCCT FindFromKey(key).
    pub fn find_from_key(&self, s: &Shape) -> Option<&Vec<Shape>> {
        self.index.get(&ShapeSet::key(s)).map(|i| &self.values[*i])
    }

    /// OCCT Seek(key) — the entry pointer (None when not bound).
    pub fn seek(&self, s: &Shape) -> Option<&Vec<Shape>> {
        self.find_from_key(s)
    }

    /// The rcad entry index of the key (the FindIndex service).
    pub fn index_of(&self, key: &(u64, u32)) -> usize {
        *self.index.get(key).unwrap()
    }

    /// Append the ancestor to the entry (the OCCT M(index).Append(anc)).
    pub fn append(&mut self, idx: usize, anc: Shape) {
        self.values[idx].push(anc);
    }
}

impl Default for IndexedShapeListMap {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L178-197 — GetFreeEdges (static).
// ---------------------------------------------------------------------------

/// OCCT static GetFreeEdges (cxx L178-197).
pub(crate) fn get_free_edges(brep: &mut BRep, a_shape: &Shape, map_edges: &mut ShapeSet) -> bool {
    for a_exp_f in topexp_explorer(brep, a_shape, ShapeType::Face) {
        for a_exp_e in topexp_explorer(brep, &a_exp_f, ShapeType::Edge) {
            let edge = a_exp_e;
            if !map_edges.contains_shape(&edge) {
                map_edges.add(&edge);
            } else {
                map_edges.remove(&edge);
            }
        }
    }
    !map_edges.is_empty()
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L206-297 — GetConnectedFaceGroups (static).
// ---------------------------------------------------------------------------

/// OCCT static GetConnectedFaceGroups (cxx L206-297): groups the connected
/// faces into separate sequences using the existing connectivity data (the
/// depth-first search through shared edges; the groups are sorted by size
/// with the largest group first).
pub(crate) fn get_connected_face_groups(
    the_face_edges: &IndexedShapeListMap,
    the_edge_faces: &HashMap<(u64, u32), Vec<Shape>>,
) -> Vec<Vec<Shape>> {
    let mut a_connected_groups: Vec<Vec<Shape>> = Vec::new();

    if the_face_edges.extent() == 0 {
        return a_connected_groups;
    }

    let mut a_visited_faces: HashSet<(u64, u32)> = HashSet::new();

    for a_face_idx in 1..=(the_face_edges.extent()) {
        let a_start_face = the_face_edges.find_key(a_face_idx).clone();
        let a_start_face_key = (a_start_face.ptr_id(), a_start_face.location);

        if a_visited_faces.contains(&a_start_face_key) {
            continue;
        }

        // Start new connected group
        let mut a_connected_group: Vec<Shape> = Vec::new();

        // DFS traversal (the std::stack of the OCCT body).
        let mut a_stack: Vec<Shape> = Vec::new();
        a_stack.push(a_start_face);
        a_visited_faces.insert(a_start_face_key);

        while let Some(a_current_face) = a_stack.pop() {
            a_connected_group.push(a_current_face.clone());

            // Find connected faces through shared edges
            let a_face_edges_iter = the_face_edges.seek(&a_current_face);
            if let Some(a_face_edges_array) = a_face_edges_iter {
                for an_edge in a_face_edges_array.iter() {
                    let an_edge_faces_ptr = the_edge_faces.get(&(an_edge.ptr_id(), an_edge.location));
                    if let Some(a_connected_faces) = an_edge_faces_ptr {
                        for a_neighbor_face in a_connected_faces.iter() {
                            let nkey = (a_neighbor_face.ptr_id(), a_neighbor_face.location);
                            if !a_visited_faces.contains(&nkey) {
                                a_visited_faces.insert(nkey);
                                a_stack.push(a_neighbor_face.clone());
                            }
                        }
                    }
                }
            }
        }

        // Insert in sorted order (largest groups first)
        let mut an_is_inserted = false;
        for an_iter in 0..a_connected_groups.len() {
            if a_connected_group.len() > a_connected_groups[an_iter].len() {
                a_connected_groups.insert(an_iter, a_connected_group.clone());
                an_is_inserted = true;
                break;
            }
        }

        if !an_is_inserted {
            a_connected_groups.push(a_connected_group);
        }
    }

    a_connected_groups
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L308-648 — GetShells (static).
// ---------------------------------------------------------------------------

/// OCCT static GetShells (cxx L308-648): creates shells from the connected
/// face groups using the connectivity analysis; processes only the largest
/// connected group for shell construction.  `the_lfaces` returns the
/// unprocessed faces.
pub(crate) fn get_shells(
    brep: &mut BRep,
    the_lfaces: &mut Vec<Shape>,
    the_map_multi_connect_edges: &ShapeSet,
    the_seq_shells: &mut Vec<Shape>,
    the_map_face_shells: &mut HashMap<(u64, u32), Shape>,
    the_err_faces: &mut Vec<Shape>,
) -> bool {
    let mut a_done = false;
    if the_lfaces.is_empty() {
        return false;
    }
    let mut nshell = brep.add_tshell(Vec::new());
    let an_is_multi_connex = !the_map_multi_connect_edges.is_empty();
    let mut a_face_idx = 1i32;
    let mut a_faces_in_shell_count = 1i32;
    let mut a_seq_unconnect_faces: Vec<Shape> = Vec::new();

    // OCCT L327-334: the EdgeOrientedMap / TempProcessedEdges.
    let mut a_processed_edges: HashMap<(u64, u32), (bool, bool)> = HashMap::new();

    // OCCT L334-354: the face-edges connectivity.
    let mut a_face_edges = IndexedShapeListMap::new();
    let mut a_number_of_edges = 0usize;
    for an_face_iter in the_lfaces.iter() {
        let a_face = an_face_iter.clone();
        let mut a_temp_edges: Vec<Shape> = Vec::new();
        for an_edge_exp in topexp_explorer(brep, &a_face, ShapeType::Edge) {
            a_temp_edges.push(an_edge_exp);
            a_number_of_edges += 1;
        }
        a_face_edges.add(a_face.clone(), a_temp_edges);
    }

    // OCCT L356-391: the edge-faces connectivity.
    let mut a_edge_faces: HashMap<(u64, u32), Vec<Shape>> = HashMap::new();

    for a_face_ind in 1..=(a_face_edges.extent()) {
        let a_face = a_face_edges.find_key(a_face_ind).clone();
        let a_face_edges_array = a_face_edges.find_from_index(a_face_ind).clone();

        for an_edge in a_face_edges_array.iter() {
            let ekey = (an_edge.ptr_id(), an_edge.location);
            let a_faces_array = a_edge_faces.entry(ekey).or_default();

            // Check if the face already exists in the array
            let mut a_face_exists = false;
            for a_face_check_idx in 0..a_faces_array.len() {
                if a_faces_array[a_face_check_idx].is_same(&a_face) {
                    a_face_exists = true;
                    break;
                }
            }

            if !a_face_exists {
                a_faces_array.push(a_face.clone());
            }
        }
    }

    // OCCT L393-404: get the connected groups of faces.
    let a_connected_groups = get_connected_face_groups(&a_face_edges, &a_edge_faces);

    if a_connected_groups.is_empty() {
        return false;
    }

    // Some assumption that each edge can be in two orientations
    a_number_of_edges = (a_number_of_edges / 2) + 1;

    let mut a_processing_faces: Vec<Shape> = a_connected_groups[0].clone();

    let mut a_temp_processed_edges: HashMap<(u64, u32), (bool, bool)> = HashMap::new();
    loop {
        if !(a_face_idx <= a_processing_faces.len() as i32) {
            break;
        }
        a_temp_processed_edges.clear();

        let mut a_bad_orientation_count = 0i32;
        let mut a_good_orientation_count = 0i32;
        let mut f1 = a_processing_faces[(a_face_idx - 1) as usize].clone();
        // Get edges of the face
        let a_face_edges_array = a_face_edges
            .seek(&f1)
            .cloned()
            .unwrap_or_default();

        for an_edge_ind in 0..a_face_edges_array.len() {
            let edge = a_face_edges_array[an_edge_ind].clone();
            let ekey = (edge.ptr_id(), edge.location);

            // if multiconnexity mode is equal to true faces contains the same
            // multiconnexity edges are not added to one shell.
            if an_is_multi_connex && the_map_multi_connect_edges.contains_shape(&edge) {
                continue;
            }

            let a_processed_edge_it = a_processed_edges.get_mut(&ekey);

            match a_processed_edge_it {
                None => {
                    let a_temp_processed_edge_it = a_temp_processed_edges.get_mut(&ekey);
                    match a_temp_processed_edge_it {
                        None => {
                            let an_edge_orientation_pair =
                                (edge.orientation == Orientation::Forward, edge.orientation == Orientation::Reversed);
                            a_temp_processed_edges.insert(ekey, an_edge_orientation_pair);
                        }
                        Some(a_temp_processed_edge_it) => {
                            a_temp_processed_edge_it.0 =
                                a_temp_processed_edge_it.0 || (edge.orientation == Orientation::Forward);
                            a_temp_processed_edge_it.1 =
                                a_temp_processed_edge_it.1 || (edge.orientation == Orientation::Reversed);
                        }
                    }
                    continue;
                }
                Some(a_pair) => {
                    let is_direct = a_pair.0;
                    let is_reversed = a_pair.1;

                    if (edge.orientation == Orientation::Forward && is_direct)
                        || (edge.orientation == Orientation::Reversed && is_reversed)
                    {
                        a_bad_orientation_count += 1;
                    } else if (edge.orientation == Orientation::Forward && is_reversed)
                        || (edge.orientation == Orientation::Reversed && is_direct)
                    {
                        a_good_orientation_count += 1;
                    }

                    if is_direct {
                        a_pair.0 = false;
                    } else if is_reversed {
                        a_pair.1 = false;
                    }

                    if !a_pair.0 && !a_pair.1 {
                        // if edge is processed in this face it is removed from
                        // the map of processed edges
                        a_processed_edges.remove(&ekey);
                    }
                }
            }
        }

        if a_bad_orientation_count == 0 && a_good_orientation_count == 0 && a_temp_processed_edges.is_empty() {
            a_face_idx += 1;
            continue;
        }

        // if face can not be added to shell it added to sequence of error faces.

        if a_good_orientation_count != 0 && a_bad_orientation_count != 0 {
            the_err_faces.push(f1.clone());
            a_processing_faces.remove((a_face_idx - 1) as usize);
            a_faces_in_shell_count += 1;
            continue;
        }

        // Addition of face to shell. In the dependance of orientation faces in
        // the shell added face can be reversed.

        if (a_good_orientation_count != 0 || a_bad_orientation_count != 0) || a_faces_in_shell_count == 1 {
            if a_bad_orientation_count != 0 {
                // OCCT L507: F1.Reverse().
                f1.orientation = match f1.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };

                let a_temp_iter: Vec<((u64, u32), (bool, bool))> =
                    a_temp_processed_edges.iter().map(|(k, v)| (*k, *v)).collect();
                for (edge_key, an_edge_orientation_pair) in a_temp_iter {
                    let a_reverted_pair = (!an_edge_orientation_pair.0, !an_edge_orientation_pair.1);
                    a_processed_edges.insert(edge_key, a_reverted_pair);
                }
                a_done = true;
            } else {
                let a_temp_iter: Vec<((u64, u32), (bool, bool))> =
                    a_temp_processed_edges.iter().map(|(k, v)| (*k, *v)).collect();
                for (edge_key, an_edge_orientation_pair) in a_temp_iter {
                    a_processed_edges.insert(edge_key, an_edge_orientation_pair);
                }
            }
            a_faces_in_shell_count += 1;
            builder_add(brep, &nshell, &f1);
            the_map_face_shells.insert((f1.ptr_id(), f1.location), nshell.clone());
            a_processing_faces.remove((a_face_idx - 1) as usize);

            // check if closed shell is obtained in multi connex mode and add
            // to sequence of shells and new shell begin to construct.
            // (check is n*2)
            if an_is_multi_connex && crate::shhealing::shape_build::brep_tool::brep_tool_is_closed(brep, &nshell) {
                set_shell_closed(brep, &nshell, true);
                the_seq_shells.push(nshell.clone());
                let nshellnext = brep.add_tshell(Vec::new());
                nshell = nshellnext;
                a_faces_in_shell_count = 1;
            }

            a_face_idx = 0;
        }
        // if shell contains of one face. This face is added to sequence of
        // faces. This shell is removed.
        if !a_processing_faces.is_empty()
            && a_face_idx == a_processing_faces.len() as i32
            && a_faces_in_shell_count <= 2
        {
            let mut a_itf = iter_subshapes(brep, &nshell, true, true).into_iter();
            if let Some(first) = a_itf.next() {
                a_seq_unconnect_faces.push(first.clone());
                the_map_face_shells.remove(&(first.ptr_id(), first.location));
            }
            let nshellnext = brep.add_tshell(Vec::new());
            nshell = nshellnext;
            a_face_idx = 0;
            a_faces_in_shell_count = 1;
        }
        a_face_idx += 1;
    }
    let mut is_contains = false;
    for k in 1..=(the_seq_shells.len()) {
        if is_contains {
            break;
        }
        is_contains = nshell.is_same(&the_seq_shells[k - 1]);
    }
    if !is_contains {
        let mut num_face = 0i32;
        let mut a_face = Shape::null();
        for a_itf in iter_subshapes(brep, &nshell, true, true) {
            a_face = a_itf;
            num_face += 1;
        }
        if num_face > 1 {
            // close all closed shells in no multi connex mode
            if !an_is_multi_connex {
                set_shell_closed(brep, &nshell, crate::shhealing::shape_build::brep_tool::brep_tool_is_closed(brep, &nshell));
            }
            the_seq_shells.push(nshell.clone());
        } else if num_face == 1 {
            the_map_face_shells.remove(&(a_face.ptr_id(), a_face.location));
            a_processing_faces.push(a_face);
        }
    }

    // Add all unprocessed connected groups (second group and after) to
    // unconnected faces
    for (a_group_index, a_unprocessed_group) in a_connected_groups.iter().enumerate() {
        if a_group_index == 0 {
            continue; // Skip first group (already processed)
        }
        for an_unproc_face_idx in 1..=(a_unprocessed_group.len()) {
            a_seq_unconnect_faces.push(a_unprocessed_group[an_unproc_face_idx - 1].clone());
        }
    }

    *the_lfaces = a_processing_faces;

    // Add unconnected faces from the largest group that couldn't be added to
    // shells
    for j1 in 1..=(a_seq_unconnect_faces.len()) {
        the_lfaces.push(a_seq_unconnect_faces[j1 - 1].clone());
    }

    a_done
}

/// OCCT `TopoDS_Shape::Closed(flag)` — the shell Closed-flag write.
pub(crate) fn set_shell_closed(brep: &mut BRep, shell: &Shape, flag: bool) {
    use rcad_kernel::topods::tshape_flags;
    // SAFETY: the in-place TShape mutation through the shared pool handle
    // (the TopoDS_TShape flag write of OCCT).
    let ptr = std::sync::Arc::as_ptr(&brep.tshapes[shell.index]) as *mut TShape;
    let ts = unsafe { &mut *ptr };
    if let TShape::Shell(sd) = ts {
        if flag {
            sd.flags |= tshape_flags::CLOSED;
        } else {
            sd.flags &= !tshape_flags::CLOSED;
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L656-913 — AddMultiConexityFaces (static).
// ---------------------------------------------------------------------------

/// OCCT static AddMultiConexityFaces (cxx L656-913): the faces having only a
/// multiconnexity boundary are added to the shells having a free boundary
/// containing the same multiconnexity edges.
pub(crate) fn add_multi_conexity_faces(
    brep: &mut BRep,
    lface: &mut Vec<Shape>,
    a_map_multi_connect_edges: &ShapeSet,
    seq_shells: &mut Vec<Shape>,
    a_map_face_shells: &HashMap<(u64, u32), Shape>,
    a_map_edge_faces: &IndexedShapeListMap,
    err_faces: &mut Vec<Shape>,
    non_manifold: bool,
) -> bool {
    let mut done = false;
    //  BRep_Builder aB;
    let mut ll_posible_shells: Vec<Shape> = Vec::new();
    let mut add_shapes: Vec<Shape> = Vec::new();
    for i1 in 1..=(lface.len()) {
        let a_shape = lface[i1 - 1].clone();

        let mut a_nb_mult_edges = 0i32;

        // Finds faces having only multiconnexity boundary.
        for a_it_wires in iter_subshapes(brep, &a_shape, false, true) {
            let mut a_nb_edges = 0i32;
            for a_it_edges in iter_subshapes(brep, &a_it_wires, false, true) {
                a_nb_edges += 1;
                let edge = a_it_edges;
                if !a_map_multi_connect_edges.contains_shape(&edge) {
                    continue;
                }
                a_nb_mult_edges += 1;
            }
            if a_nb_mult_edges == 0 {
                continue;
            }

            if a_nb_mult_edges == a_nb_edges {
                add_shapes.push(a_shape.clone());
            } else {
                ll_posible_shells.push(a_shape.clone());
            }
        }
    }

    // Attempt to create a shell from the unconnected which have not only
    // multiconnexity boundary.
    let mut a_tmp_shells: Vec<Shape> = Vec::new();
    if !ll_posible_shells.is_empty() {
        let mut a_map = ShapeSet::new();
        let mut a_tmp: Vec<Shape> = Vec::new();
        let mut a_tmp_face_shell: HashMap<(u64, u32), Shape> = HashMap::new();
        let mut ll = ll_posible_shells.clone();
        if get_shells(brep, &mut ll, &a_map, &mut a_tmp_shells, &mut a_tmp_face_shell, &mut a_tmp) {
            for kk in 1..=(a_tmp_shells.len()) {
                let a_sh = a_tmp_shells[kk - 1].clone();
                let mut map_edges = ShapeSet::new();
                if get_free_edges(brep, &a_sh, &mut map_edges) {
                    let mut nbedge = 0i32;
                    for amap_iter in map_edges.iter() {
                        if a_map_multi_connect_edges.contains_shape(amap_iter) {
                            nbedge += 1;
                        }
                    }
                    if nbedge != 0 && nbedge == map_edges.extent() as i32 {
                        add_shapes.push(a_sh);
                    }
                }
            }
        }
    }

    // Add chosen faces to shells.
    for k1 in 1..=(add_shapes.len()) {
        let mut map_other_shells: HashMap<(u64, u32), i32> = HashMap::new();
        let mut dire = ShapeSet::new();
        let mut reve = ShapeSet::new();
        let a_sh = add_shapes[k1 - 1].clone();
        let mut map_edges = ShapeSet::new();
        if !get_free_edges(brep, &a_sh, &mut map_edges) {
            continue;
        }
        let mut lfaces: Vec<Shape> = Vec::new();

        // Fill MapOtherShells which will contain shells with orientation in
        // which the selected shape aSh will be added.
        for amap_iter in map_edges.iter() {
            if !a_map_multi_connect_edges.contains_shape(amap_iter) {
                continue;
            }
            let edge = amap_iter;
            if edge.orientation == Orientation::Forward {
                dire.add(edge);
            } else {
                reve.add(edge);
            }
            if let Some(lf) = a_map_edge_faces.find_from_key(edge) {
                for f in lf.iter() {
                    lfaces.push(f.clone());
                }
            }
        }
        for a_itl in lfaces.iter() {
            let a_f = a_itl;
            if !a_map_face_shells.contains_key(&(a_f.ptr_id(), a_f.location)) {
                continue;
            }

            let a_othershell = a_map_face_shells
                .get(&(a_f.ptr_id(), a_f.location))
                .cloned()
                .unwrap();
            if map_other_shells.contains_key(&(a_othershell.ptr_id(), a_othershell.location)) {
                continue;
            }
            if !non_manifold && crate::shhealing::shape_build::brep_tool::brep_tool_is_closed(brep, &a_othershell) {
                continue;
            }

            let mut map_shell_edges = ShapeSet::new();
            get_free_edges(brep, &a_othershell, &mut map_shell_edges);
            let mut is_add = true;
            for amap_iter1 in map_edges.iter() {
                if !is_add {
                    break;
                }
                is_add = map_shell_edges.contains_shape(amap_iter1);
            }

            if !is_add {
                continue;
            }
            let mut nbdir = 0i32;
            let mut nbrev = 0i32;

            // add only free face whome all edges contains in the shell as open
            // boundary.
            for a_ite in map_shell_edges.iter() {
                let edge_s = a_ite;
                if !a_map_multi_connect_edges.contains_shape(edge_s) {
                    continue;
                }
                if (edge_s.orientation == Orientation::Forward && dire.contains_shape(edge_s))
                    || (edge_s.orientation == Orientation::Reversed && reve.contains_shape(edge_s))
                {
                    nbrev += 1;
                } else if (edge_s.orientation == Orientation::Forward && reve.contains_shape(edge_s))
                    || (edge_s.orientation == Orientation::Reversed && dire.contains_shape(edge_s))
                {
                    nbdir += 1;
                }
            }
            if nbdir != 0 && nbrev != 0 {
                err_faces.push(a_sh.clone());
                continue;
            }
            if nbdir != 0 || nbrev != 0 {
                let is_reverse = if nbrev != 0 { 1 } else { 0 };
                map_other_shells.insert((a_othershell.ptr_id(), a_othershell.location), is_reverse);
            }
        }
        if map_other_shells.is_empty() {
            //      i1++;
            continue;
        }

        // Adds face to open shells containing the same multishared edges.
        // For nonmanifold mode creation of one shell from the face and the
        // shells containing the same multishared edges.
        //  If one face can be added to a few shells (case of compsolid) face
        //  will be added to each shell.
        done = true;
        let mut first_rev = 0i32;
        let mut first_ind = 0usize;
        let mut ind = 0i32;
        let mut l = 0usize;
        while l < seq_shells.len() {
            let shell_key = (seq_shells[l].ptr_id(), seq_shells[l].location);
            if !map_other_shells.contains_key(&shell_key) {
                l += 1;
                continue;
            }
            ind += 1;
            let is_rev = map_other_shells.get(&shell_key).copied().unwrap_or(0);
            let mut anew_shape = a_sh.clone();
            if is_rev != 0 {
                // OCCT L868: aSh.Reversed().
                anew_shape.orientation = match anew_shape.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    o => o,
                };
            }

            if ind == 1 || !non_manifold {
                if ind == 1 {
                    first_rev = is_rev;
                    first_ind = l;
                }
                for a_e in topexp_explorer(brep, &anew_shape, ShapeType::Face) {
                    builder_add(brep, &seq_shells[l], &a_e);
                }
                let shell = seq_shells[l].clone();
                let _ = shell;
            } else if non_manifold {
                let is_reversed = ((is_rev != 0) || (first_rev != 0))
                    && (!(is_rev != 0) || !(first_rev != 0));
                // OCCT L888-895: the shells are merged into the first one.
                for a_it_f in iter_subshapes(brep, &seq_shells[l], false, true) {
                    let mut n_f = a_it_f;
                    if is_reversed {
                        n_f.orientation = match n_f.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            o => o,
                        };
                    }
                    builder_add(brep, &seq_shells[first_ind], &n_f);
                }
                seq_shells.remove(l);
                continue;
            }
            l += 1;
        }

        dire.clear();
        reve.clear();
        for a_et in topexp_explorer(brep, &a_sh, ShapeType::Face) {
            let mut kk = 0usize;
            while kk < lface.len() {
                if a_et.is_same(&lface[kk]) {
                    lface.remove(kk);
                } else {
                    kk += 1;
                }
            }
        }
    }
    done
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L919-942 — BoxIn (static).
// ---------------------------------------------------------------------------

/// OCCT static BoxIn (cxx L919-942): checks if one face contains inside
/// other.
pub(crate) fn box_in(the_box1: &BndBox, the_box2: &BndBox) -> i32 {
    let mut a_num_in = 0i32;
    let corners1 = the_box1.get();
    let corners2 = the_box2.get();
    let (a_xmin1, a_ymin1, a_zmin1, a_xmax1, a_ymax1, a_zmax1) = match corners1 {
        Some(v) => v,
        None => return a_num_in,
    };
    let (a_xmin2, a_ymin2, a_zmin2, a_xmax2, a_ymax2, a_zmax2) = match corners2 {
        Some(v) => v,
        None => return a_num_in,
    };
    if a_xmin1 == a_xmin2
        && a_xmax1 == a_xmax2
        && a_ymin1 == a_ymin2
        && a_ymax1 == a_ymax2
        && a_zmin1 == a_zmin2
        && a_zmax1 == a_zmax2
    {
        a_num_in = 0;
    } else if a_xmin1 >= a_xmin2
        && a_xmax1 <= a_xmax2
        && a_ymin1 >= a_ymin2
        && a_ymax1 <= a_ymax2
        && a_zmin1 >= a_zmin2
        && a_zmax1 <= a_zmax2
    {
        a_num_in = 1;
    } else if a_xmin1 <= a_xmin2
        && a_xmax1 >= a_xmax2
        && a_ymin1 <= a_ymin2
        && a_ymax1 >= a_ymax2
        && a_zmin1 <= a_zmin2
        && a_zmax1 >= a_zmax2
    {
        a_num_in = 2;
    }
    a_num_in
}

// ---------------------------------------------------------------------------
// OCCT BRepBndLib::AddClose — GAP carrier (module doc).
// ---------------------------------------------------------------------------

/// OCCT `BRepBndLib::AddClose(S, B)` (TKTopAlgo/BRepBndLib.cxx) — GAP
/// carrier: the BRepBndLib bounding machinery is not translated; the box
/// stays void, so `BoxIn` reports no nesting and `GetClosedShells` keeps all
/// candidate shells (the OCCT "no candidate nested" path).
pub(crate) fn brep_bnd_lib_add_close_gap(_brep: &mut BRep, _s: &Shape, _b: &mut BndBox) {}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L950-986 — GetClosedShells (static).
// ---------------------------------------------------------------------------

/// OCCT static GetClosedShells (cxx L950-986): checks if one shell is a part
/// of another shell (the compsolid case when a few shells are created from
/// the same set of faces).
pub(crate) fn get_closed_shells(brep: &mut BRep, shells: &mut Vec<Shape>, a_remain_shells: &mut Vec<Shape>) {
    let mut a_boxes: Vec<BndBox> = Vec::new();
    for i in 1..=(shells.len()) {
        let mut box3d = BndBox::new();
        brep_bnd_lib_add_close_gap(brep, &shells[i - 1], &mut box3d);
        a_boxes.push(box3d);
    }
    let mut a_map_num: HashSet<i32> = HashSet::new();
    for j in 1..=(a_boxes.len()) {
        for k in (j + 1)..=(a_boxes.len()) {
            let num_in = box_in(&a_boxes[j - 1], &a_boxes[k - 1]);
            match num_in {
                1 => {
                    a_map_num.insert(k as i32);
                }
                2 => {
                    a_map_num.insert(j as i32);
                }
                _ => {}
            }
        }
    }
    for i1 in 1..=(shells.len()) {
        if !a_map_num.contains(&(i1 as i32)) {
            a_remain_shells.push(shells[i1 - 1].clone());
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L992-1164 — GlueClosedCandidate (static).
// ---------------------------------------------------------------------------

/// OCCT static GlueClosedCandidate (cxx L992-1164): first, attempts to
/// create closed shells from the sequence of open shells.
pub(crate) fn glue_closed_candidate(
    brep: &mut BRep,
    open_shells: &mut Vec<Shape>,
    a_map_multi_connect_edges: &ShapeSet,
    a_seq_new_shells: &mut Vec<Shape>,
) {
    // Creating new shells if some open shells contain the same free boundary.
    let mut i = 1usize;
    while i < open_shells.len() {
        let a_shell = open_shells[i - 1].clone();
        let mut map_edges1 = ShapeSet::new();
        let mut dire = ShapeSet::new();
        let mut reve = ShapeSet::new();
        if !get_free_edges(brep, &a_shell, &mut map_edges1) {
            // OCCT L1162: OpenShells.Remove(i--) — kept for the step parity.
            i += 1;
            continue;
        }

        for a_ite in map_edges1.iter() {
            let edge1 = a_ite;
            if !a_map_multi_connect_edges.contains_shape(edge1) {
                break;
            }
            if edge1.orientation == Orientation::Forward {
                dire.add(edge1);
            } else if edge1.orientation == Orientation::Reversed {
                reve.add(edge1);
            }
        }
        if map_edges1.extent() > (dire.extent() + reve.extent()) {
            i += 1;
            continue;
        }

        // Filling map MapOtherShells which contains candidates to creation of
        // a closed shell with aShell.
        let mut map_other_shells: HashMap<(u64, u32), bool> = HashMap::new();
        for j in (i + 1)..=(open_shells.len()) {
            let mut is_add_shell = true;
            let mut is_reversed = false;
            let mut map_edges2 = ShapeSet::new();
            let a_shell2 = open_shells[j - 1].clone();
            if !get_free_edges(brep, &a_shell2, &mut map_edges2) {
                continue;
            }
            for a_ite2 in map_edges2.iter() {
                if !is_add_shell {
                    break;
                }
                let edge2 = a_ite2;
                if !a_map_multi_connect_edges.contains_shape(edge2) {
                    is_add_shell = false;
                    break;
                    // continue;
                }
                is_add_shell = dire.contains_shape(edge2) || reve.contains_shape(edge2);
                if (edge2.orientation == Orientation::Forward && dire.contains_shape(edge2))
                    || (edge2.orientation == Orientation::Reversed && reve.contains_shape(edge2))
                {
                    is_reversed = true;
                }
            }

            if !is_add_shell {
                continue;
            }
            map_other_shells.insert(
                (open_shells[j - 1].ptr_id(), open_shells[j - 1].location),
                is_reversed,
            );
        }
        if map_other_shells.is_empty() {
            i += 1;
            continue;
        }

        if !map_other_shells.is_empty() {
            // Case of compsolid when more than two shells have the same free
            // boundary.
            let mut a_seq_candidate: Vec<Shape> = Vec::new();
            a_seq_candidate.push(open_shells[i - 1].clone());

            for (key, _v) in map_other_shells.iter() {
                if let Some(sh) = open_shells.iter().find(|s| (s.ptr_id(), s.location) == *key) {
                    a_seq_candidate.push(sh.clone());
                }
            }

            // Creation of all possible shells from the chosen candidates.
            //  And the addition of them to the temporary sequence.
            let mut a_tmp_seq: Vec<Shape> = Vec::new();
            for k in 1..=(a_seq_candidate.len()) {
                for l in (k + 1)..=(a_seq_candidate.len()) {
                    let a_new_sh = brep.add_tshell(Vec::new());
                    for a_it1 in iter_subshapes(brep, &a_seq_candidate[k - 1], false, true) {
                        builder_add(brep, &a_new_sh, &a_it1);
                    }
                    let mut is_rev = *map_other_shells
                        .get(&(a_seq_candidate[l - 1].ptr_id(), a_seq_candidate[l - 1].location))
                        .unwrap_or(&false);
                    if k != 1 {
                        is_rev = is_rev
                            == *map_other_shells
                                .get(&(a_seq_candidate[k - 1].ptr_id(), a_seq_candidate[k - 1].location))
                                .unwrap_or(&false);
                    }
                    for a_exp in topexp_explorer(brep, &a_seq_candidate[l - 1], ShapeType::Face) {
                        let mut a_face = a_exp;
                        if is_rev {
                            // OCCT L1111: aExp.Current().Reversed().
                            a_face.orientation = match a_face.orientation {
                                Orientation::Forward => Orientation::Reversed,
                                Orientation::Reversed => Orientation::Forward,
                                o => o,
                            };
                        }
                        builder_add(brep, &a_new_sh, &a_face);
                    }
                    a_tmp_seq.push(a_new_sh);
                }
            }

            // Choice from the temporary sequence of the shells containing
            // different sets of faces (case of compsolid)
            let mut a_remain_shells: Vec<Shape> = Vec::new();
            get_closed_shells(brep, &mut a_tmp_seq, &mut a_remain_shells);
            a_seq_new_shells.append(&mut a_remain_shells);

            let mut j1 = i + 1;
            while j1 <= open_shells.len() {
                let key = (open_shells[j1 - 1].ptr_id(), open_shells[j1 - 1].location);
                if !map_other_shells.contains_key(&key) {
                    j1 += 1;
                    continue;
                }
                open_shells.remove(j1 - 1);
            }
        }
        // OCCT L1162: OpenShells.Remove(i--) — the scan restarts at the same
        // index (the i-- and the loop ++ cancel).
        open_shells.remove(i - 1);
        if i > open_shells.len() {
            break;
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L1171-1320 — CreateNonManifoldShells (static).
// ---------------------------------------------------------------------------

/// OCCT static CreateNonManifoldShells (cxx L1171-1320): attempts to create
/// the max possible shells from the open shells.
pub(crate) fn create_non_manifold_shells(
    brep: &mut BRep,
    seq_shells: &mut Vec<Shape>,
    a_map_multi_connect_edges: &ShapeSet,
) {
    // OCCT L1175-1205.
    let mut a_map = IndexedShapeListMap::new();
    for i in 1..=(seq_shells.len()) {
        let a_shell = seq_shells[i - 1].clone();
        let mut medeg = ShapeSet::new();
        // OCCT L1181: TopExp::MapShapes(aShell, TopAbs_EDGE, medeg).
        for e in topexp_explorer(brep, &a_shell, ShapeType::Edge) {
            medeg.add(&e);
        }
        for mit in a_map_multi_connect_edges.iter() {
            let ae = mit;
            // if( aMapMultiConnectEdges.Contains(aExp.Current())) {
            if medeg.contains_shape(ae) {
                if let Some(list) = a_map.seek(ae) {
                    let mut list = list.clone();
                    list.push(a_shell.clone());
                    // (replace the entry value in place)
                    let idx = *a_map.index.get(&ShapeSet::key(ae)).unwrap();
                    a_map.values[idx] = list;
                } else {
                    a_map.add(ae.clone(), vec![a_shell.clone()]);
                }
            }
        }
    }
    // OCCT L1206-1295.
    let mut a_map_shells: HashMap<(u64, u32), Shape> = HashMap::new();
    for j in 1..=(a_map.extent()) {
        let l_shells = a_map.find_from_index(j).clone();
        let a_new_shell = brep.add_tshell(Vec::new());
        let mut mapmerge = ShapeSet::new();
        let mut ismerged = false;
        let mut num = 1i32;
        for alit in l_shells.iter() {
            let alit_key = (alit.ptr_id(), alit.location);
            if !a_map_shells.contains_key(&alit_key) {
                for a_ef in topexp_explorer(brep, alit, ShapeType::Face) {
                    builder_add(brep, &a_new_shell, &a_ef);
                }
                ismerged = true;
                mapmerge.add(alit);
            } else if ismerged {
                let mut arshell = a_map_shells.get(&alit_key).cloned().unwrap();
                loop {
                    if !a_map_shells.contains_key(&(arshell.ptr_id(), arshell.location)) {
                        break;
                    }
                    let ss = a_map_shells
                        .get(&(arshell.ptr_id(), arshell.location))
                        .cloned()
                        .unwrap();
                    if ss.is_same(&arshell) {
                        break;
                    }
                    arshell = ss;
                }

                if !mapmerge.contains_shape(&arshell) {
                    for a_ef in topexp_explorer(brep, &arshell, ShapeType::Face) {
                        builder_add(brep, &a_new_shell, &a_ef);
                    }
                    mapmerge.add(&arshell);
                }
            } else {
                let mut arshell = a_map_shells.get(&alit_key).cloned().unwrap();
                loop {
                    if !a_map_shells.contains_key(&(arshell.ptr_id(), arshell.location)) {
                        break;
                    }
                    let ss = a_map_shells
                        .get(&(arshell.ptr_id(), arshell.location))
                        .cloned()
                        .unwrap();
                    if ss.is_same(&arshell) {
                        break;
                    }
                    arshell = ss;
                }
                if num == 1 {
                    for a_ef in topexp_explorer(brep, &arshell, ShapeType::Face) {
                        builder_add(brep, &a_new_shell, &a_ef);
                    }

                    mapmerge.add(&arshell);
                } else if !mapmerge.contains_shape(&arshell) {
                    for a_ef in topexp_explorer(brep, &arshell, ShapeType::Face) {
                        builder_add(brep, &a_new_shell, &a_ef);
                    }
                    mapmerge.add(&arshell);
                }
            }
            num += 1;
        }
        if mapmerge.extent() > 1 || ismerged {
            for alit1 in mapmerge.iter() {
                a_map_shells.insert((alit1.ptr_id(), alit1.location), a_new_shell.clone());
            }
        }
    }
    // OCCT L1296-1319.
    let mut map_new_shells = ShapeSet::new();
    let mut nn = 0usize;
    while nn < seq_shells.len() {
        let key = (seq_shells[nn].ptr_id(), seq_shells[nn].location);
        if a_map_shells.contains_key(&key) {
            let mut a_new_shell = a_map_shells.get(&key).cloned().unwrap();
            loop {
                if !a_map_shells.contains_key(&(a_new_shell.ptr_id(), a_new_shell.location)) {
                    break;
                }
                let ss = a_map_shells
                    .get(&(a_new_shell.ptr_id(), a_new_shell.location))
                    .cloned()
                    .unwrap();
                if ss.is_same(&a_new_shell) {
                    break;
                }
                a_new_shell = ss;
            }
            map_new_shells.add(&a_new_shell);

            seq_shells.remove(nn);
        } else {
            nn += 1;
        }
    }
    for ii in map_new_shells.iter() {
        seq_shells.push(ii.clone());
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell.cxx L1327-1421 — CreateClosedShell (static).
// ---------------------------------------------------------------------------

/// OCCT static CreateClosedShell (cxx L1327-1421): attempts to create the
/// max possible shells from the open shells.
pub(crate) fn create_closed_shell(
    brep: &mut BRep,
    open_shells: &mut Vec<Shape>,
    a_map_multi_connect_edges: &ShapeSet,
) {
    let mut a_new_shells: Vec<Shape> = Vec::new();
    // First, attempt to create closed shells.
    glue_closed_candidate(brep, open_shells, a_map_multi_connect_edges, &mut a_new_shells);

    // Creating new shells if some open shells contain the multishared same
    // edges.
    let mut i = 1usize;
    while i < open_shells.len() {
        let mut is_add_shell = false;
        let mut a_shell = open_shells[i - 1].clone();
        let mut is_reversed = false;
        let mut j = i + 1usize;
        while j <= open_shells.len() {
            let mut map_edges1 = ShapeSet::new();
            let mut dire = ShapeSet::new();
            let mut reve = ShapeSet::new();
            if !get_free_edges(brep, &a_shell, &mut map_edges1) {
                break;
            }
            for a_ite in map_edges1.iter() {
                let edge1 = a_ite;
                if !a_map_multi_connect_edges.contains_shape(edge1) {
                    continue;
                }
                if edge1.orientation == Orientation::Forward {
                    dire.add(edge1);
                } else if edge1.orientation == Orientation::Reversed {
                    reve.add(edge1);
                }
            }
            if dire.is_empty() && reve.is_empty() {
                break;
            }
            let mut map_edges2 = ShapeSet::new();
            let a_shell2 = open_shells[j - 1].clone();
            if !get_free_edges(brep, &a_shell2, &mut map_edges2) {
                j += 1;
                continue;
            }
            for a_ite2 in map_edges2.iter() {
                let edge2 = a_ite2;
                if !a_map_multi_connect_edges.contains_shape(edge2) {
                    continue;
                }
                if !dire.contains_shape(edge2) && !reve.contains_shape(edge2) {
                    continue;
                }
                is_add_shell = true;
                if (edge2.orientation == Orientation::Forward && dire.contains_shape(edge2))
                    || (edge2.orientation == Orientation::Reversed && reve.contains_shape(edge2))
                {
                    is_reversed = true;
                }
            }

            if !is_add_shell {
                j += 1;
                continue;
            }

            for a_exp_f21 in topexp_explorer(brep, &open_shells[j - 1], ShapeType::Face) {
                let mut a_face = a_exp_f21;
                if is_reversed {
                    // OCCT L1409: aFace.Reverse().
                    a_face.orientation = match a_face.orientation {
                        Orientation::Forward => Orientation::Reversed,
                        Orientation::Reversed => Orientation::Forward,
                        o => o,
                    };
                }
                builder_add(brep, &a_shell, &a_face);
            }

            open_shells[i - 1] = a_shell.clone();
            open_shells.remove(j - 1);
        }
        i += 1;
    }

    open_shells.append(&mut a_new_shells);
}
