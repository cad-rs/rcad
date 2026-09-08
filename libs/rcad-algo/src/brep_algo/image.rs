//! OCCT BRepAlgo_Image — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_Image.cxx
//!         (L28-353) + BRepAlgo_Image.hxx (L34-100).
//!
//! Architecture differences:
//! 1. NCollection_DataMap keyed by TopTools_ShapeMapHasher ->
//!    HashMap<ShapeKey, _>.  The `up` map keeps the key shape (its
//!    orientation is read in Filter); `down` values keep nothing extra.
//! 2. The OCCT function-local `static` list in Image(S) (cxx L207-209) maps
//!    to the owned list returned by image() (see the accessor comment).

use crate::brep_algo::tool::{shape_key, ShapeKey};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
use std::collections::HashMap;

/// OCCT BRepAlgo_Image (BRepAlgo_Image.hxx L34-100) — stores link between a
/// shape S and a shape NewS obtained from S; NewS is an image of S.
pub struct BRepAlgoImage {
    /// OCCT: roots.
    roots: Vec<Shape>,
    /// OCCT: up (image shape -> its generator).
    up: HashMap<ShapeKey, (Shape, Shape)>,
    /// OCCT: down (shape -> its images).
    down: HashMap<ShapeKey, Vec<Shape>>,
}

impl Default for BRepAlgoImage {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepAlgoImage {
    /// OCCT BRepAlgo_Image::BRepAlgo_Image() (cxx L28).
    pub fn new() -> Self {
        BRepAlgoImage {
            roots: Vec::new(),
            up: HashMap::new(),
            down: HashMap::new(),
        }
    }

    /// OCCT BRepAlgo_Image::SetRoot(S) (cxx L32-35).
    pub fn set_root(&mut self, s: &Shape) {
        self.roots.push(s.clone());
    }

    /// OCCT BRepAlgo_Image::Bind(OldS, NewS) (cxx L39-50) — links NewS as
    /// image of OldS.
    pub fn bind(&mut self, old_s: &Shape, new_s: &Shape) {
        if self.down.contains_key(&shape_key(old_s)) {
            panic!(" BRepAlgo_Image::Bind");
        }
        let ko = shape_key(old_s);
        self.down.insert(ko, vec![new_s.clone()]);
        let kn = shape_key(new_s);
        self.up.insert(kn, (new_s.clone(), old_s.clone()));
    }

    /// OCCT BRepAlgo_Image::Bind(OldS, L) (cxx L54-73) — links NewS as image
    /// of OldS.
    pub fn bind_list(&mut self, old_s: &Shape, l: &[Shape]) {
        if self.has_image(old_s) {
            panic!(" BRepAlgo_Image::Bind");
        }
        for it in l {
            if !self.has_image(old_s) {
                self.bind(old_s, it);
            } else {
                self.add(old_s, it);
            }
        }
    }

    /// OCCT BRepAlgo_Image::Clear() (cxx L77-82).
    pub fn clear(&mut self) {
        self.roots.clear();
        self.up.clear();
        self.down.clear();
    }

    /// OCCT BRepAlgo_Image::Add(OldS, NewS) (cxx L86-94) — adds NewS to the
    /// image of OldS.
    pub fn add(&mut self, old_s: &Shape, new_s: &Shape) {
        if !self.has_image(old_s) {
            panic!(" BRepAlgo_Image::Add");
        }
        let ko = shape_key(old_s);
        if let Some(l) = self.down.get_mut(&ko) {
            l.push(new_s.clone());
        }
        let kn = shape_key(new_s);
        self.up.insert(kn, (new_s.clone(), old_s.clone()));
    }

    /// OCCT BRepAlgo_Image::Add(OldS, L) (cxx L98-105) — adds NewS to the
    /// image of OldS.
    pub fn add_list(&mut self, old_s: &Shape, l: &[Shape]) {
        for it in l {
            self.add(old_s, it);
        }
    }

    /// OCCT BRepAlgo_Image::Remove(S) (cxx L109-132) — removes S from the
    /// set of images.
    pub fn remove(&mut self, s: &Shape) {
        let ks = shape_key(s);
        let old_s = match self.up.get(&ks) {
            Some((_, generator)) => generator.clone(),
            None => panic!(" BRepAlgo_Image::Remove"),
        };
        let ko = shape_key(&old_s);
        if let Some(l) = self.down.get_mut(&ko) {
            // OCCT L117-126: remove the first entry IsSame(S), then break.
            for i in 0..l.len() {
                if l[i].is_same(s) {
                    l.remove(i);
                    break;
                }
            }
        }
        // OCCT L127-131.
        let now_empty = self.down.get(&ko).map_or(true, |l| l.is_empty());
        if now_empty {
            self.down.remove(&ko);
        }
        self.up.remove(&ks);
    }

    /// OCCT BRepAlgo_Image::Roots() (cxx L139-142).
    pub fn roots(&self) -> &[Shape] {
        &self.roots
    }

    /// OCCT BRepAlgo_Image::IsImage(S) (cxx L146-149).
    pub fn is_image(&self, s: &Shape) -> bool {
        self.up.contains_key(&shape_key(s))
    }

    /// OCCT BRepAlgo_Image::ImageFrom(S) (cxx L153-160) — returns the
    /// generator of S.
    pub fn image_from(&self, s: &Shape) -> &Shape {
        match self.up.get(&shape_key(s)) {
            Some((_, generator)) => generator,
            None => panic!(" BRepAlgo_Image::ImageFrom"),
        }
    }

    /// OCCT BRepAlgo_Image::Root(S) (cxx L164-189) — returns the upper
    /// generator of S.
    pub fn root(&self, s: &Shape) -> &Shape {
        if !self.up.contains_key(&shape_key(s)) {
            panic!(" BRepAlgo_Image::FirstImageFrom");
        }

        // OCCT L171-172.
        let mut s1 = self.up.get(&shape_key(s)).expect("up(S)").1.clone();
        let mut s2 = s.clone();

        if s1.is_same(&s2) {
            return &self.up.get(&shape_key(s)).expect("up(S)").1;
        }

        while self.up.contains_key(&shape_key(&s1)) {
            s2 = s1.clone();
            s1 = self.up.get(&shape_key(&s1)).expect("up(S1)").1.clone();
            if s1.is_same(&s2) {
                break;
            }
        }
        &self.up.get(&shape_key(&s2)).expect("up(S2)").1
    }

    /// OCCT BRepAlgo_Image::HasImage(S) (cxx L193-196).
    pub fn has_image(&self, s: &Shape) -> bool {
        self.down.contains_key(&shape_key(s))
    }

    /// OCCT BRepAlgo_Image::Image(S) (cxx L203-212) — returns the Image of S;
    /// returns S in the list if HasImage(S) is false.  The OCCT function-local
    /// `static` list (cxx L207-209) maps to the returned owned list (the
    /// documented contract is "S in the list when HasImage(S) is false"; the
    /// OCCT static's never-cleared accumulation is an implementation artifact
    /// and is not carried over).
    pub fn image(&self, s: &Shape) -> Vec<Shape> {
        if !self.has_image(s) {
            let mut l: Vec<Shape> = Vec::new();
            l.push(s.clone());
            return l;
        }
        self.down.get(&shape_key(s)).cloned().unwrap_or_default()
    }

    /// OCCT BRepAlgo_Image::LastImage(S, L) (cxx L218-239) — stores in L the
    /// images of images of...images of S; L contains only S if HasImage(S)
    /// is false.
    pub fn last_image(&self, s: &Shape, l: &mut Vec<Shape>) {
        if !self.down.contains_key(&shape_key(s)) {
            l.push(s.clone());
        } else {
            let down_list = self.down.get(&shape_key(s)).cloned().unwrap_or_default();
            for it in &down_list {
                if it.is_same(s) {
                    l.push(s.clone());
                } else {
                    self.last_image(it, l);
                }
            }
        }
    }

    /// OCCT BRepAlgo_Image::Compact() (cxx L243-266) — keeps only the link
    /// between roots and lastimage.
    pub fn compact(&mut self) {
        // OCCT L245: M — the new down-map built from the roots.
        let mut m: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();
        let roots_snapshot = self.roots.clone();
        for s in &roots_snapshot {
            let mut li: Vec<Shape> = Vec::new();
            if self.has_image(s) {
                self.last_image(s, &mut li);
            }
            m.insert(shape_key(s), li);
        }
        // OCCT L257-258.
        self.up.clear();
        self.down.clear();
        // OCCT L259-265.
        for it in &roots_snapshot {
            if m.contains_key(&shape_key(it)) {
                let li = m.get(&shape_key(it)).cloned().unwrap_or_default();
                // OCCT Bind(S, LI): an empty LI binds nothing (the loop body
                // of Bind(OldS, L) never runs) — bind_list reproduces that.
                self.bind_list(it, &li);
            }
        }
    }

    /// OCCT BRepAlgo_Image::Filter(S, T) (cxx L270-295) — deletes in the
    /// images the shape of type T which are not in S.  Warning: Compact()
    /// must be called before.
    pub fn filter(&mut self, s: &Shape, t: ShapeType) {
        // OCCT L273-278: M = the shapes of type T in S.
        let mut m: std::collections::HashSet<ShapeKey> = std::collections::HashSet::new();
        for e in crate::brep_algo::tool::explorer(s, t, ShapeType::Shape) {
            m.insert(shape_key(&e));
        }
        // OCCT L279-294: remove every up-keyed image of type T not in M,
        // restarting the scan after each removal.
        let mut change = true;
        while change {
            change = false;
            let keys: Vec<ShapeKey> = self.up.keys().copied().collect();
            for k in keys {
                let (a_s, _) = match self.up.get(&k) {
                    Some(entry) => entry.clone(),
                    None => continue,
                };
                if a_s.shape_type() == t && !m.contains(&k) {
                    self.remove(&a_s);
                    change = true;
                    break;
                }
            }
        }
    }

    /// OCCT BRepAlgo_Image::RemoveRoot(Root) (cxx L299-330) — removes the
    /// root from the list of roots and up and down maps.
    pub fn remove_root(&mut self, root: &Shape) {
        // OCCT L301-311.
        let mut is_removed = false;
        for i in 0..self.roots.len() {
            if root.is_same(&self.roots[i]) {
                self.roots.remove(i);
                is_removed = true;
                break;
            }
        }
        if !is_removed {
            return;
        }

        // OCCT L317-329.
        let kr = shape_key(root);
        let new_s = match self.down.get(&kr).cloned() {
            Some(l) => l,
            None => return,
        };
        for it in &new_s {
            if let Some((_, old_s)) = self.up.get(&shape_key(it)) {
                if old_s.is_same(root) {
                    self.up.remove(&shape_key(it));
                }
            }
        }
        self.down.remove(&kr);
    }

    /// OCCT BRepAlgo_Image::ReplaceRoot(OldRoot, NewRoot) (cxx L334-353) —
    /// replaces the OldRoot with the NewRoot so all images of the OldRoot
    /// become images of the NewRoot; the OldRoot is removed.
    pub fn replace_root(&mut self, old_root: &Shape, new_root: &Shape) {
        if !self.has_image(old_root) {
            return;
        }

        let a_l_image = self.image(old_root);
        if self.has_image(new_root) {
            self.add_list(new_root, &a_l_image);
        } else {
            self.bind_list(new_root, &a_l_image);
        }

        self.set_root(new_root);
        self.remove_root(old_root);
    }
}
