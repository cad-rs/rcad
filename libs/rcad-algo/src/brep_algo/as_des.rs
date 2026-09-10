//! OCCT BRepAlgo_AsDes — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_AsDes.cxx
//!         (L25-298) + BRepAlgo_AsDes.hxx (L30-84).
//! OCCT inheritance chain: BRepAlgo_AsDes : Standard_Transient.
//!
//! Architecture differences:
//! 1. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//!    TopTools_ShapeMapHasher> -> HashMap<ShapeKey, Vec<Shape>> (the key
//!    shape is never read back in this class).
//! 2. Standard_Transient RTTI has no Rust equivalent (no ref-count probe
//!    needed here).

use crate::brep_algo::tool::{shape_key, ShapeKey};
use rcad_kernel::topo_shape::Shape;
use std::collections::HashMap;

/// OCCT BRepAlgo_AsDes (BRepAlgo_AsDes.hxx L30-84) — SD to store descendants
/// and ascendants of Shapes.
pub struct BRepAlgoAsDes {
    /// OCCT: up.
    up: HashMap<ShapeKey, Vec<Shape>>,
    /// OCCT: down.
    down: HashMap<ShapeKey, Vec<Shape>>,
}

impl Default for BRepAlgoAsDes {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepAlgoAsDes {
    /// OCCT BRepAlgo_AsDes::BRepAlgo_AsDes() (cxx L27) — creates an empty
    /// AsDes.
    pub fn new() -> Self {
        BRepAlgoAsDes {
            up: HashMap::new(),
            down: HashMap::new(),
        }
    }

    /// OCCT BRepAlgo_AsDes::Clear() (cxx L61-65).
    pub fn clear(&mut self) {
        self.up.clear();
        self.down.clear();
    }

    /// OCCT BRepAlgo_AsDes::Add(S, SS) (cxx L31-46) — stores SS as a futur
    /// subshape of S.
    pub fn add(&mut self, s: &Shape, ss: &Shape) {
        // OCCT L33-38.
        let k = shape_key(s);
        if !self.down.contains_key(&k) {
            self.down.insert(k, Vec::new());
        }
        self.down.get_mut(&k).expect("down").push(ss.clone());

        // OCCT L40-45.
        let ks = shape_key(ss);
        if !self.up.contains_key(&ks) {
            self.up.insert(ks, Vec::new());
        }
        self.up.get_mut(&ks).expect("up").push(s.clone());
    }

    /// OCCT BRepAlgo_AsDes::Add(S, SS list) (cxx L50-57) — stores SS as
    /// futurs SubShapes of S.
    pub fn add_list(&mut self, s: &Shape, ss: &[Shape]) {
        for it in ss {
            self.add(s, it);
        }
    }

    /// OCCT BRepAlgo_AsDes::HasAscendant(S) (cxx L69-72).
    pub fn has_ascendant(&self, s: &Shape) -> bool {
        self.up.contains_key(&shape_key(s))
    }

    /// OCCT BRepAlgo_AsDes::HasDescendant(S) (cxx L76-79).
    pub fn has_descendant(&self, s: &Shape) -> bool {
        self.down.contains_key(&shape_key(s))
    }

    /// OCCT BRepAlgo_AsDes::Ascendant(S) (cxx L83-91) — returns the Shape
    /// containing S (the static empty list maps to a shared empty slice).
    pub fn ascendant(&self, s: &Shape) -> &[Shape] {
        match self.up.get(&shape_key(s)) {
            Some(l) => l,
            None => &[],
        }
    }

    /// OCCT BRepAlgo_AsDes::Descendant(S) (cxx L95-103) — returns futur
    /// subhapes of S.
    pub fn descendant(&self, s: &Shape) -> &[Shape] {
        match self.down.get(&shape_key(s)) {
            Some(l) => l,
            None => &[],
        }
    }

    /// OCCT BRepAlgo_AsDes::ChangeDescendant(S) (cxx L107-115) — returns
    /// futur subhapes of S.
    ///
    /// OCCT returns a reference to the static empty list when S is unbound
    /// (mutation through it is lost); the rcad call sites re-bind on the
    /// bound case only, so the unbound case returns None here.
    pub fn change_descendant(&mut self, s: &Shape) -> Option<&mut Vec<Shape>> {
        self.down.get_mut(&shape_key(s))
    }

    /// OCCT BRepAlgo_AsDes::Replace(OldS, NewS) (cxx L221-278) — replaces
    /// theOldS by theNewS; theOldS disappears from this.
    pub fn replace(&mut self, old_s: &Shape, new_s: &Shape) {
        let old_key = shape_key(old_s);
        let new_key = shape_key(new_s);
        for i in 0..2 {
            // OCCT L225-226: aMap = !i ? up : down; the BackReplace target is
            // the other map (InUp = i != 0 selects `up`) — disjoint field
            // borrows keep the OCCT pointer-aliasing shape.
            let (a_map, other_map) =
                if i == 0 { (&mut self.up, &mut self.down) } else { (&mut self.down, &mut self.up) };
            // OCCT L227-231: pLSOld = aMap.ChangeSeek(OldS); skip when
            // unbound.
            if !a_map.contains_key(&old_key) {
                continue;
            }
            // OCCT L233-234: InUp = i != 0; BackReplace(OldS, NewS, *pLSOld,
            // InUp).  The const-ref pass is an owned snapshot of the old
            // list (the map itself is what BackReplace mutates).
            let ls_old_snapshot = a_map.get(&old_key).cloned().expect("pLSOld");
            back_replace(old_s, new_s, &ls_old_snapshot, i != 0, other_map);

            // OCCT L236-238: pLSNew = aMap.ChangeSeek(NewS).
            if !a_map.contains_key(&new_key) {
                // OCCT L239-253: filter the list, then Bind(NewS, *pLSOld).
                let mut a_ms: std::collections::HashSet<ShapeKey> =
                    std::collections::HashSet::new();
                let mut ls_old = a_map.get(&old_key).cloned().expect("pLSOld");
                ls_old.retain(|it| a_ms.insert(shape_key(it)));
                a_map.insert(new_key, ls_old);
            } else {
                // OCCT L257-274: avoid duplicates — collect aMS from the
                // new list, then append the old-list entries not in aMS.
                let mut a_ms: std::collections::HashSet<ShapeKey> =
                    std::collections::HashSet::new();
                if let Some(ls_new) = a_map.get(&new_key) {
                    for it in ls_new {
                        a_ms.insert(shape_key(it));
                    }
                }
                let ls_old = a_map.get(&old_key).cloned().expect("pLSOld");
                if let Some(ls_new) = a_map.get_mut(&new_key) {
                    for a_s in &ls_old {
                        if a_ms.insert(shape_key(a_s)) {
                            ls_new.push(a_s.clone());
                        }
                    }
                }
            }
            // OCCT L276: aMap.UnBind(OldS).
            a_map.remove(&old_key);
        }
    }

    /// OCCT BRepAlgo_AsDes::Remove(SS) (cxx L282-298) — removes SS from me.
    pub fn remove(&mut self, ss: &Shape) {
        let k = shape_key(ss);
        if self.down.contains_key(&k) {
            panic!(" BRepAlgo_AsDes::Remove");
        }
        if !self.up.contains_key(&k) {
            panic!(" BRepAlgo_AsDes::Remove");
        }
        let ups = self.up.get(&k).cloned().expect("up(SS)");
        for it in &ups {
            if let Some(l) = self.down.get_mut(&shape_key(it)) {
                remove_in_list(ss, l);
            }
        }
        self.up.remove(&k);
    }

    /// OCCT BRepAlgo_AsDes::HasCommonDescendant(S1, S2, LC) (cxx L166-189) —
    /// returns True if S1 and S2 have common Descendants; stores in LC the
    /// commons Descendants.
    pub fn has_common_descendant(&self, s1: &Shape, s2: &Shape, lc: &mut Vec<Shape>) -> bool {
        lc.clear();
        if self.has_descendant(s1) && self.has_descendant(s2) {
            let ds_list = self.descendant(s1).to_vec();
            for ds1 in &ds_list {
                let ads_list = self.ascendant(ds1).to_vec();
                for ads1 in &ads_list {
                    if ads1.is_same(s2) {
                        lc.push(ds1.clone());
                    }
                }
            }
        }
        !lc.is_empty()
    }
}

/// OCCT BRepAlgo_AsDes::BackReplace(OldS, NewS, L, InUp) (cxx L193-217) —
/// replaces theOldS by theNewS; theOldS disappears from this.  The map
/// selected by theInUp (`up` when true, `down` otherwise) is passed
/// explicitly (disjoint field borrow of the OCCT this-pointer aliasing).
fn back_replace(
    old_s: &Shape,
    new_s: &Shape,
    l: &[Shape],
    in_up: bool,
    target: &mut HashMap<ShapeKey, Vec<Shape>>,
) {
    let _ = in_up;
    for it in l {
        let s = it;
        if let Some(found) = target.get_mut(&shape_key(s)) {
            replace_in_list(old_s, new_s, found);
        }
    }
}

/// OCCT static ReplaceInList (cxx L119-146).
fn replace_in_list(old_s: &Shape, new_s: &Shape, l: &mut Vec<Shape>) {
    // OCCT L123-128: aMS collects the list identities.
    let mut a_ms: std::collections::HashSet<ShapeKey> = std::collections::HashSet::new();
    for it in l.iter() {
        a_ms.insert(shape_key(it));
    }
    // OCCT L129-145: walk the list; every OldS entry is replaced by
    // NewS.Oriented(O) (deduplicated through aMS) and removed.
    let mut i = 0usize;
    while i < l.len() {
        if l[i].is_same(old_s) {
            let o = l[i].orientation;
            let mut new_oriented = new_s.clone();
            new_oriented.orientation = o;
            if a_ms.insert(shape_key(&new_oriented)) {
                l.insert(i, new_oriented);
                i += 1;
            }
            l.remove(i);
        } else {
            i += 1;
        }
    }
}

/// OCCT static RemoveInList (cxx L150-162).
fn remove_in_list(s: &Shape, l: &mut Vec<Shape>) {
    for i in 0..l.len() {
        if l[i].is_same(s) {
            l.remove(i);
            break;
        }
    }
}
