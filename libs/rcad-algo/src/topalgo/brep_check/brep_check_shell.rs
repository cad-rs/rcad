//! OCCT BRepCheck_Shell (TKTopAlgo/BRepCheck).
//!
//! Source: `$OCCT_SRC/src/ModelingAlgorithms/TKTopAlgo/BRepCheck/BRepCheck_Shell.cxx`
//! (L39-1007) and `BRepCheck_Shell.hxx` (L30-77).

use std::collections::HashMap;
use std::sync::Arc;

use rcad_kernel::topo::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepTool, Orientation, ShapeType, TShape, TShellData};

use super::brep_check_result::{
    brep_check_add, explorer, ShapeKey, BRepCheckResultBase, BRepCheckStatus,
};
use super::brep_check_wire::{IndexedShapeMap, ShapeSet};

/// OCCT Shell.cxx L107-110: `IsOriented(S)` (a local copy in the OCCT file).
pub fn is_oriented_shell(s: &Shape) -> bool {
    s.orientation == Orientation::Forward || s.orientation == Orientation::Reversed
}

/// OCCT BRepCheck_Shell (Shell.hxx L30-77).
#[derive(Debug)]
pub struct BRepCheckShell {
    /// OCCT protected base (Result.hxx L82-90).
    pub base: BRepCheckResultBase,
    /// OCCT myNbori.
    pub my_nbori: i32,
    /// OCCT myCdone.
    pub my_cdone: bool,
    /// OCCT myCstat.
    pub my_cstat: BRepCheckStatus,
    /// OCCT myOdone.
    pub my_odone: bool,
    /// OCCT myOstat.
    pub my_ostat: BRepCheckStatus,
    /// OCCT myMapEF.
    pub my_map_ef: IndexedShapeMap,
}

impl BRepCheckShell {
    /// OCCT BRepCheck_Shell::BRepCheck_Shell(const TopoDS_Shell& S)
    /// (Shell.cxx L114-122).
    pub fn new(brep: &BRep, s: &Shape) -> Self {
        let mut r = BRepCheckShell {
            base: BRepCheckResultBase::new(),
            my_nbori: 0,
            my_cdone: false,
            my_cstat: BRepCheckStatus::NoError,
            my_odone: false,
            my_ostat: BRepCheckStatus::NoError,
            my_map_ef: IndexedShapeMap::new(),
        };
        r.base.init(s);
        r.minimum(brep);
        r
    }

    /// OCCT BRepCheck_Shell::Minimum (Shell.cxx L126-184).
    pub fn minimum(&mut self, brep: &BRep) {
        // OCCT L127-128.
        self.my_cdone = false;
        self.my_odone = false;

        if !self.base.my_min {
            // OCCT L133-135.
            self.base.my_map.bound(&self.base.my_shape);
            let my_shape = self.base.my_shape.clone();

            // OCCT L138-157: it is checked if the shell is "connected".
            let faces = explorer(brep, &my_shape, ShapeType::Face);
            let mut nbface = 0i32;
            self.my_map_ef.clear();
            for f in &faces {
                nbface += 1;
                for edg in explorer(brep, f, ShapeType::Edge) {
                    let index = self.my_map_ef.find_index_or_add(&edg);
                    self.my_map_ef.value_mut(index).push(f.clone());
                }
            }

            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");
            if nbface == 0 {
                // OCCT L159-162.
                brep_check_add(lst, BRepCheckStatus::EmptyShell);
            } else if nbface >= 2 {
                // OCCT L163-174.
                let mut map_f = ShapeSet::new();
                let first = faces[0].clone();
                shell_propagate(brep, &self.my_map_ef, &first, &mut map_f);
                if map_f.extent() as i32 != nbface {
                    brep_check_add(lst, BRepCheckStatus::NotConnected);
                }
            }

            // OCCT L176-183.
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Minimum: myShape must be bound");
            if lst.is_empty() {
                lst.push(BRepCheckStatus::NoError);
            }
            self.my_map_ef.clear();
            self.base.my_min = true;
        }
    }

    /// OCCT BRepCheck_Shell::InContext (Shell.cxx L188-248).
    pub fn in_context(&mut self, brep: &BRep, s: &Shape) {
        // OCCT L190-204: bound check under the (parallel) lock.
        if self.base.my_map.is_bound(s) {
            return;
        }
        self.base.my_map.bind(s.clone(), Vec::new());
        let lst = self
            .base
            .my_map
            .find_mut(s)
            .expect("InContext: the context must be bound");

        // OCCT L208-220.
        let mut found = false;
        let mut more = false;
        for cur in explorer(brep, s, ShapeType::Shell) {
            more = true;
            if cur.is_same(&self.base.my_shape) {
                found = true;
                break;
            }
        }
        if !more || !found {
            brep_check_add(lst, BRepCheckStatus::SubshapeNotInShape);
            return;
        }

        // OCCT L222-242. (The OCCT code keeps `lst` alive across the Closed/
        // Orientation calls; they write the myShape entry, a different map
        // slot — the rcad borrows are re-taken at each write.)
        match s.shape_type() {
            ShapeType::Solid => {
                let fst = self.closed(brep, false);
                if (fst == BRepCheckStatus::NotClosed && brep.is_closed(s))
                    || (fst != BRepCheckStatus::NoError)
                {
                    let lst = self
                        .base
                        .my_map
                        .find_mut(s)
                        .expect("InContext: the context must be bound");
                    brep_check_add(lst, fst);
                } else if !self.is_unorientable() {
                    let fst = self.orientation(brep, false);
                    let lst = self
                        .base
                        .my_map
                        .find_mut(s)
                        .expect("InContext: the context must be bound");
                    brep_check_add(lst, fst);
                }
            }
            _ => {
                // OCCT L240-242: default — nothing.
            }
        }

        // OCCT L244-247.
        let lst = self
            .base
            .my_map
            .find_mut(s)
            .expect("InContext: the context must be bound");
        if lst.is_empty() {
            lst.push(BRepCheckStatus::NoError);
        }
    }

    /// OCCT BRepCheck_Shell::Blind (Shell.cxx L252-259).
    pub fn blind(&mut self) {
        if !self.base.my_blind {
            self.base.my_blind = true;
        }
    }

    /// OCCT BRepCheck_Shell::Closed (Shell.cxx L263-445).
    pub fn closed(&mut self, brep: &BRep, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        if self.my_cdone {
            // OCCT L276-284.
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed: myShape must be bound");
                brep_check_add(lst, self.my_cstat);
            }
            return self.my_cstat;
        }

        // OCCT L286.
        self.my_cdone = true; // it will be done...

        // OCCT L288-293: an already stored error is returned as-is.
        {
            let list = self
                .base
                .my_map
                .find(&my_shape)
                .expect("Closed: myShape must be bound");
            if let Some(&first) = list.first() {
                if first != BRepCheckStatus::NoError {
                    self.my_cstat = first;
                    return self.my_cstat; // already saved
                }
            }
        }

        self.my_cstat = BRepCheckStatus::NoError;

        // OCCT L297-300.
        let mut map_s = ShapeSet::new();
        let mut a_me_to_avoid = ShapeSet::new();
        self.my_map_ef.clear();

        // OCCT L307-323: the non-oriented edges of the oriented faces.
        for a_f in explorer(brep, &my_shape, ShapeType::Face) {
            if is_oriented_shell(&a_f) {
                for a_e in explorer(brep, &a_f, ShapeType::Edge) {
                    if !is_oriented_shell(&a_e) {
                        a_me_to_avoid.add(&a_e);
                    }
                }
            }
        }

        // OCCT L326-366: the oriented face/edge connectivity.
        for a_f in explorer(brep, &my_shape, ShapeType::Face) {
            if is_oriented_shell(&a_f) {
                if !map_s.add(&a_f) {
                    // OCCT L332-342.
                    self.my_cstat = BRepCheckStatus::RedundantFace;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("Closed: myShape must be bound");
                        brep_check_add(lst, self.my_cstat);
                    }
                    return self.my_cstat;
                }

                for a_e in explorer(brep, &a_f, ShapeType::Edge) {
                    // OCCT L351: if (!aMEToAvoid.Contains(aE)).
                    if !a_me_to_avoid.contains(&a_e) {
                        let index = self.my_map_ef.find_index_or_add(&a_e);
                        self.my_map_ef.value_mut(index).push(a_f.clone());
                    }
                }
            }
        }

        // OCCT L369-386.
        self.my_nbori = map_s.extent() as i32;
        if self.my_nbori >= 2 {
            map_s.clear();
            // Search for the first oriented face.
            let mut a_f: Option<Shape> = None;
            for f in explorer(brep, &my_shape, ShapeType::Face) {
                a_f = Some(f.clone());
                if is_oriented_shell(&f) {
                    break;
                }
            }
            if let Some(f) = a_f {
                shell_propagate(brep, &self.my_map_ef, &f, &mut map_s);
            }
        }

        // OCCT L390-399.
        let a_nb_f = map_s.extent() as i32;
        if self.my_nbori != a_nb_f {
            self.my_cstat = BRepCheckStatus::NotConnected;
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Closed: myShape must be bound");
                brep_check_add(lst, self.my_cstat);
            }
            return self.my_cstat;
        }

        // OCCT L401-438: the per-edge connection counts.
        let nbedges = self.my_map_ef.extent();
        for i in 0..nbedges {
            let nboc = self.my_map_ef.value(i).len() as i32;
            if nboc == 0 || nboc >= 3 {
                let mut the_set: Vec<Shape> = Vec::new();
                let nb_set = self.nb_connected_set(brep, &mut the_set);
                // If there is more than one closed cavity the shell is
                // considered invalid — this corresponds to the criteria of a
                // solid (not those of a shell).
                if nb_set > 1 {
                    self.my_cstat = BRepCheckStatus::InvalidMultiConnexity;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("Closed: myShape must be bound");
                        brep_check_add(lst, self.my_cstat);
                    }
                    return self.my_cstat;
                }
            } else if nboc == 1 {
                let edge = self.my_map_ef.find_key(i).clone();
                if !brep.is_edge_degenerated(&edge) {
                    self.my_cstat = BRepCheckStatus::NotClosed;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("Closed: myShape must be bound");
                        brep_check_add(lst, self.my_cstat);
                    }
                    return self.my_cstat;
                }
            }
        }

        // OCCT L440-444.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Closed: myShape must be bound");
            brep_check_add(lst, self.my_cstat);
        }
        self.my_cstat
    }

    /// OCCT BRepCheck_Shell::Orientation (Shell.cxx L449-853).
    pub fn orientation(&mut self, brep: &BRep, update: bool) -> BRepCheckStatus {
        let my_shape = self.base.my_shape.clone();
        if self.my_odone {
            // OCCT L462-469.
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Orientation: myShape must be bound");
                brep_check_add(lst, self.my_ostat);
            }
            return self.my_ostat;
        }
        // OCCT L470.
        self.my_odone = true;

        // OCCT L472-480.
        self.my_ostat = self.closed(brep, false);
        if self.my_ostat != BRepCheckStatus::NotClosed
            && self.my_ostat != BRepCheckStatus::NoError
        {
            if update {
                let lst = self
                    .base
                    .my_map
                    .find_mut(&my_shape)
                    .expect("Orientation: myShape must be bound");
                brep_check_add(lst, self.my_ostat);
            }
            return self.my_ostat;
        }

        self.my_ostat = BRepCheckStatus::NoError;

        // OCCT L484-504: the orientation of each face in relation to the shell
        // (used to check BRepCheck_RedundantFace).
        let mut map_of_shape_orientation: HashMap<ShapeKey, Orientation> = HashMap::new();
        for f in explorer(brep, &my_shape, ShapeType::Face) {
            let key = ShapeKey::of(&f);
            if map_of_shape_orientation.contains_key(&key) {
                // OCCT L492-503.
                self.my_ostat = BRepCheckStatus::RedundantFace;
                if update {
                    let lst = self
                        .base
                        .my_map
                        .find_mut(&my_shape)
                        .expect("Orientation: myShape must be bound");
                    brep_check_add(lst, self.my_ostat);
                    break;
                } else {
                    return self.my_ostat;
                }
            } else {
                map_of_shape_orientation.insert(key, f.orientation);
            }
        }

        // OCCT L520-699: the orientation of faces by their connectivity
        // (BRepCheck_BadOrientationOfSubshape / BRepCheck_SubshapeNotInShape).
        let nbedges = self.my_map_ef.extent();
        let mut fref: Option<Shape> = None;

        for i in 0..nbedges {
            let edg = self.my_map_ef.find_key(i).clone();
            if brep.is_edge_degenerated(&edg) {
                continue;
            }
            let lface = self.my_map_ef.value(i).clone();

            if lface.len() <= 2 {
                // OCCT L539-631.
                if lface.is_empty() {
                    continue;
                }
                let fref_raw = lface[0].clone();
                fref = Some(fref_raw.clone());

                let key_ref = ShapeKey::of(&fref_raw);
                let Some(&orf_ref) = map_of_shape_orientation.get(&key_ref) else {
                    // OCCT L544-553.
                    self.my_ostat = BRepCheckStatus::SubshapeNotInShape;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("Orientation: myShape must be bound");
                        brep_check_add(lst, self.my_ostat);
                    }
                    // quit because no workaround for the incoherence is possible
                    return self.my_ostat;
                };
                let fref_s = oriented_local(&fref_raw, orf_ref);

                if lface.len() > 1 {
                    // Edge of connectivity.
                    let fcur_raw = lface[1].clone();
                    // OCCT L565-574: the edge on Fref.
                    let edges_of_fref = explorer(brep, &fref_s, ShapeType::Edge);
                    let orient = edges_of_fref
                        .iter()
                        .find(|e| e.is_same(&edg))
                        .map(|e| e.orientation);
                    let Some(orient) = orient else {
                        continue;
                    };

                    let key_cur = ShapeKey::of(&fcur_raw);
                    let Some(&orf_cur) = map_of_shape_orientation.get(&key_cur) else {
                        // OCCT L576-585.
                        self.my_ostat = BRepCheckStatus::SubshapeNotInShape;
                        if update {
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("Orientation: myShape must be bound");
                            brep_check_add(lst, self.my_ostat);
                        }
                        return self.my_ostat;
                    };
                    let fcur = oriented_local(&fcur_raw, orf_cur);

                    // OCCT L593-600.
                    let edges_of_fcur = explorer(brep, &fcur, ShapeType::Edge);
                    let found_pos = edges_of_fcur.iter().position(|e| e.is_same(&edg));
                    let cur_orient = found_pos.map(|p| edges_of_fcur[p].orientation);
                    if let Some(cur_orient) = cur_orient {
                        if cur_orient == orient {
                            // OCCT L602-628: the loop continues while the same
                            // edge appears again in the wire.
                            let start = found_pos.unwrap() + 1;
                            let bfound = edges_of_fcur[start..].iter().any(|e| e.is_same(&edg));
                            // OCCT L619.
                            if !bfound || cur_orient == orient {
                                self.my_ostat = BRepCheckStatus::BadOrientationOfSubshape;
                                if update {
                                    let lst = self
                                        .base
                                        .my_map
                                        .find_mut(&my_shape)
                                        .expect("Orientation: myShape must be bound");
                                    brep_check_add(lst, self.my_ostat);
                                    break;
                                }
                                return self.my_ostat;
                            }
                        }
                    }
                }
            } else {
                // OCCT L633-698: more than two faces.
                let mut num_f = 0i32;
                let mut num_r = 0i32;
                let mut fmap = ShapeSet::new();

                for fcur_raw in &lface {
                    let key_cur = ShapeKey::of(fcur_raw);
                    let Some(&orf_cur) = map_of_shape_orientation.get(&key_cur) else {
                        // OCCT L641-650.
                        self.my_ostat = BRepCheckStatus::SubshapeNotInShape;
                        if update {
                            let lst = self
                                .base
                                .my_map
                                .find_mut(&my_shape)
                                .expect("Orientation: myShape must be bound");
                            brep_check_add(lst, self.my_ostat);
                        }
                        return self.my_ostat;
                    };
                    let fcur = oriented_local(fcur_raw, orf_cur);

                    // OCCT L657-663: the first occurrence of the edge on Fcur.
                    let edges_of_fcur = explorer(brep, &fcur, ShapeType::Edge);
                    let pos = edges_of_fcur
                        .iter()
                        .position(|e| e.is_same(&edg))
                        .map(|p| p + 1)
                        .unwrap_or(edges_of_fcur.len());
                    // OCCT L664-674: the "closed" edge on Fcur — the second
                    // occurrence is taken when Fcur is met twice.
                    let pos2 = if fmap.contains(fcur_raw) {
                        edges_of_fcur[pos..]
                            .iter()
                            .position(|e| e.is_same(&edg))
                            .map(|p| p + pos + 1)
                            .unwrap_or(edges_of_fcur.len())
                    } else {
                        pos
                    };
                    let orient = edges_of_fcur
                        .get(pos2.min(edges_of_fcur.len().saturating_sub(1)))
                        .map(|e| e.orientation)
                        .unwrap_or(Orientation::Forward);
                    if orient == Orientation::Forward {
                        num_f += 1;
                    } else {
                        num_r += 1;
                    }

                    fmap.add(fcur_raw);
                }

                if num_f != num_r {
                    // OCCT L688-697.
                    self.my_ostat = BRepCheckStatus::BadOrientationOfSubshape;
                    if update {
                        let lst = self
                            .base
                            .my_map
                            .find_mut(&my_shape)
                            .expect("Orientation: myShape must be bound");
                        brep_check_add(lst, self.my_ostat);
                        break;
                    }
                    return self.my_ostat;
                }
            }
        }

        // OCCT L701-846: if at least one incorrectly oriented face has been
        // found, check whether the shell can be oriented (the Moebius walk).
        if self.my_ostat == BRepCheckStatus::BadOrientationOfSubshape {
            if let Some(mut fref_s) = fref.clone() {
                if nbedges > 0 {
                    let mut alre = ShapeSet::new();
                    let mut voisin: Vec<Shape> = vec![fref_s.clone()];
                    alre.clear();
                    while !voisin.is_empty() {
                        fref_s = voisin.remove(0);
                        let key_ref = ShapeKey::of(&fref_s);
                        let Some(&orf_ref) = map_of_shape_orientation.get(&key_ref) else {
                            // OCCT L721-730.
                            self.my_ostat = BRepCheckStatus::SubshapeNotInShape;
                            if update {
                                let lst = self
                                    .base
                                    .my_map
                                    .find_mut(&my_shape)
                                    .expect("Orientation: myShape must be bound");
                                brep_check_add(lst, self.my_ostat);
                            }
                            return self.my_ostat;
                        };
                        fref_s = oriented_local(&fref_s, orf_ref);

                        alre.add(&fref_s);

                        for ede in explorer(brep, &fref_s, ShapeType::Edge) {
                            let edg = ede.clone();
                            let orient = edg.orientation;
                            let Some(lface) = self.my_map_ef.seek(&edg) else {
                                continue;
                            };
                            if lface.is_empty() {
                                continue;
                            }
                            let mut fcur_raw = lface[0].clone();
                            if fcur_raw.is_same(&fref_s) {
                                if lface.len() > 1 {
                                    fcur_raw = lface[1].clone();
                                } else {
                                    // from the free border one goes to the
                                    // next edge
                                    continue;
                                }
                            }

                            let key_cur = ShapeKey::of(&fcur_raw);
                            let Some(&orf_cur) = map_of_shape_orientation.get(&key_cur) else {
                                // OCCT L770-779.
                                self.my_ostat = BRepCheckStatus::SubshapeNotInShape;
                                if update {
                                    let lst = self
                                        .base
                                        .my_map
                                        .find_mut(&my_shape)
                                        .expect("Orientation: myShape must be bound");
                                    brep_check_add(lst, self.my_ostat);
                                }
                                return self.my_ostat;
                            };
                            let fcur = oriented_local(&fcur_raw, orf_cur);

                            // OCCT L794-801.
                            let cur_orient = explorer(brep, &fcur, ShapeType::Edge)
                                .into_iter()
                                .find(|e| e.is_same(&edg))
                                .map(|e| e.orientation);
                            if let Some(cur_orient) = cur_orient {
                                if cur_orient == orient {
                                    if alre.contains(&fcur) {
                                        // OCCT L803-823.
                                        self.my_ostat = BRepCheckStatus::UnorientableShape;
                                        if update {
                                            let lst = self
                                                .base
                                                .my_map
                                                .find_mut(&my_shape)
                                                .expect(
                                                    "Orientation: myShape must be bound",
                                                );
                                            brep_check_add(lst, self.my_ostat);
                                        }
                                        // quit, otherwise there is a risk of
                                        // taking too much time.
                                        return self.my_ostat;
                                    }
                                    // OCCT L825-826: MapOfShapeOrientation(Fcur)
                                    // = TopAbs::Reverse(orf).
                                    let reversed = match orf_cur {
                                        Orientation::Forward => Orientation::Reversed,
                                        Orientation::Reversed => Orientation::Forward,
                                        o => o,
                                    };
                                    map_of_shape_orientation.insert(key_cur, reversed);
                                }
                            }
                            if alre.add(&fcur) {
                                voisin.push(fcur);
                            }
                        }
                    }
                }
            }
        }

        // OCCT L848-852.
        if update {
            let lst = self
                .base
                .my_map
                .find_mut(&my_shape)
                .expect("Orientation: myShape must be bound");
            brep_check_add(lst, self.my_ostat);
        }
        self.my_ostat
    }

    /// OCCT BRepCheck_Shell::SetUnorientable (Shell.cxx L857-865).
    pub fn set_unorientable(&mut self) {
        let lst = self
            .base
            .my_map
            .find_mut(&self.base.my_shape)
            .expect("SetUnorientable: myShape must be bound");
        brep_check_add(lst, BRepCheckStatus::UnorientableShape);
    }

    /// OCCT BRepCheck_Shell::IsUnorientable (Shell.cxx L869-895).
    pub fn is_unorientable(&self) -> bool {
        if self.my_odone {
            return self.my_ostat != BRepCheckStatus::NoError;
        }
        if let Some(list) = self.base.my_map.find(&self.base.my_shape) {
            for itl in list {
                if *itl == BRepCheckStatus::UnorientableShape {
                    return true;
                }
            }
        }
        false
    }

    /// OCCT BRepCheck_Shell::NbConnectedSet (Shell.cxx L899-1007) — the
    /// number of face groups separated by multi-connected edges.
    ///
    /// Architecture note: the OCCT body accumulates the current group into a
    /// BRep_Builder-built TopoDS_Shell (CurShell) and stores the CLOSED flag
    /// on it; rcad collects the group faces in a Vec wrapped in an ephemeral
    /// Shell TShape (the flag store has no observable effect on the count).
    #[allow(unused_assignments)]
    pub fn nb_connected_set(&self, brep: &BRep, the_sets: &mut Vec<Shape>) -> usize {
        let my_shape = self.base.my_shape.clone();
        // OCCT L901-904: parents = MapShapesAndAncestors(S, EDGE, FACE).
        let mut parents: IndexedShapeMap = IndexedShapeMap::new();
        for f in explorer(brep, &my_shape, ShapeType::Face) {
            for e in explorer(brep, &f, ShapeType::Edge) {
                let i = parents.find_index_or_add(&e);
                parents.value_mut(i).push(f.clone());
            }
        }
        // OCCT L906-911: all the faces.
        let mut the_faces = ShapeSet::new();
        for exsh in explorer(brep, &my_shape, ShapeType::Face) {
            the_faces.add(&exsh);
        }
        // OCCT L913-927: the multi-connected and non-oriented edges.
        let mut the_multi_ed = ShapeSet::new();
        let mut the_un_ori_ed = ShapeSet::new();
        for i_cur in 0..parents.extent() {
            let ed = parents.find_key(i_cur).clone();
            if parents.value(i_cur).len() > 2 {
                the_multi_ed.add(&ed);
            }
            if ed.orientation != Orientation::Reversed && ed.orientation != Orientation::Forward {
                the_un_ori_ed.add(&ed);
            }
        }
        // OCCT L929-1005: the group walk.
        let mut cur_shell: Vec<Shape> = Vec::new();
        let mut new_cur = true;
        'outer: for itmsh in the_multi_ed.items().to_vec() {
            let ed = itmsh;
            if !the_un_ori_ed.contains(&ed) {
                let faces_of_ed = parents.seek(&ed).cloned().unwrap_or_default();
                for lconx1 in faces_of_ed {
                    if the_faces.contains(&lconx1) {
                        let mut ad_fac = lconx1.clone();
                        cur_shell.push(ad_fac.clone());
                        the_faces.remove(&ad_fac);
                        new_cur = false;
                        if the_faces.extent() == 0 {
                            break;
                        }
                        let mut les_cur: Vec<Shape> = vec![ad_fac.clone()];
                        while !les_cur.is_empty() {
                            ad_fac = les_cur.remove(0);
                            for exsh in explorer(brep, &ad_fac, ShapeType::Edge) {
                                let ced = exsh;
                                if !the_multi_ed.contains(&ced) {
                                    let faces_of_ced =
                                        parents.seek(&ced).cloned().unwrap_or_default();
                                    for lconx2 in faces_of_ced {
                                        if the_faces.contains(&lconx2) {
                                            ad_fac = lconx2.clone();
                                            cur_shell.push(ad_fac.clone());
                                            the_faces.remove(&ad_fac);
                                            new_cur = false;
                                            if the_faces.extent() == 0 {
                                                break;
                                            }
                                            les_cur.push(ad_fac.clone());
                                        }
                                    }
                                }
                                if the_faces.extent() == 0 {
                                    break;
                                }
                            }
                        }
                        if !new_cur {
                            // OCCT L988-992: CurShell.Closed(BRep_Tool::
                            // IsClosed(CurShell)) — the flag store on the
                            // ephemeral shell has no observable effect — then
                            // theSets.Append(CurShell).
                            let group = TShape::Shell(TShellData {
                                my_shapes: cur_shell.clone(),
                                flags: rcad_kernel::topods::tshape_flags::DEFAULT,
                                faces: cur_shell.clone(),
                            });
                            the_sets.push(Shape {
                                data: Arc::new(group),
                                index: usize::MAX,
                                location: 0,
                                orientation: Orientation::Forward,
                            });
                            cur_shell.clear();
                            new_cur = true;
                        }
                    }
                    if the_faces.extent() == 0 {
                        break 'outer;
                    }
                }
            }
            if the_faces.extent() == 0 {
                break;
            }
        }
        the_sets.len()
    }
}

/// OCCT TopoDS_Shape::Orientation(theOrientation) — the local re-orientation
/// of a copy.
fn oriented_local(s: &Shape, o: Orientation) -> Shape {
    let mut c = s.clone();
    c.orientation = o;
    c
}

/// OCCT Shell.cxx L43-77: `Propagate` — the width-first traverse of the faces
/// connected through the edges of <mapEF>.
pub fn shell_propagate(
    brep: &BRep,
    map_ef: &IndexedShapeMap,
    the_face: &Shape,
    the_map_f: &mut ShapeSet,
) {
    // OCCT L50: base for the traverse procedure.
    the_map_f.add(the_face);

    // OCCT L53-76: the breadth-first traverse (the map grows at the end).
    let mut idx = 0usize;
    while idx < the_map_f.extent() {
        let a_face = the_map_f.items()[idx].clone();
        for ex in explorer(brep, &a_face, ShapeType::Edge) {
            let edg = ex;
            // OCCT L61-66: the edge list in the map (only oriented edges are
            // present).
            let Some(a_list) = map_ef.seek(&edg) else {
                continue;
            };
            for itl in a_list.clone() {
                // OCCT L68-74: existing objects are not added; new objects go
                // to the end.
                the_map_f.add(&itl);
            }
        }
        idx += 1;
    }
}
