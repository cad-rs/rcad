//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_Shell`
//! (`ShapeAnalysis_Shell.hxx` L17-110 + `.cxx` L1-321).
//!
//! Analyzes the shell(s): detects the edges encountered more than once in
//! the same orientation (bad edges) and the edges encountered only once
//! (free edges).
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — `BRep_Builder` / `TopExp_Explorer` /
//!    `TopoDS_Iterator` read and mutate the TShape graph through
//!    `rcad_kernel::BRep` (the edge.rs bridge #1).
//! 2. `NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>` ->
//!    the local [`IndexedMapShape`] keyed by the IsSame identity (TShape
//!    pointer + location index) with an insertion-order key list.
//! 3. `TopoDS_Iterator` / `TopExp_Explorer` -> `brep_tool::iter_subshapes`
//!    (CumOri = true) / `brep_tool::topexp_explorer`.
//! 4. `TopoDS_Compound` construction -> `brep.add_tcompound(Vec<Shape>)`.

use std::collections::HashMap;

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, ShapeType, TShape};

use crate::shhealing::shape_build::brep_tool::{iter_subshapes, topexp_explorer};

/// OCCT BRep_Tool::Degenerated(E).
fn brep_tool_degenerated(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT TopoDS_Shape::IsSame(S) — the TopTools_ShapeMapHasher key: same
/// TShape and same location (the orientation is ignored).
fn shape_same_key(s: &Shape) -> (u64, u32) {
    (std::sync::Arc::as_ptr(&s.data) as u64, s.location)
}

/// OCCT NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> —
/// the insertion-order map with the IsSame key.
#[derive(Default)]
pub(crate) struct IndexedMapShape {
    keys: Vec<Shape>,
    index: HashMap<(u64, u32), usize>,
}

impl IndexedMapShape {
    /// OCCT Clear().
    pub(crate) fn clear(&mut self) {
        self.keys.clear();
        self.index.clear();
    }

    /// OCCT Add(S) — returns true when the shape was added (not present).
    pub(crate) fn add(&mut self, s: &Shape) -> bool {
        let key = shape_same_key(s);
        if self.index.contains_key(&key) {
            return false;
        }
        self.keys.push(s.clone());
        self.index.insert(key, self.keys.len());
        true
    }

    /// OCCT FindIndex(S) — the 1-based index, 0 when absent.
    pub(crate) fn find_index(&self, s: &Shape) -> i32 {
        self.index
            .get(&shape_same_key(s))
            .map(|&i| i as i32)
            .unwrap_or(0)
    }

    /// OCCT Contains(S).
    pub(crate) fn contains(&self, s: &Shape) -> bool {
        self.index.contains_key(&shape_same_key(s))
    }

    /// OCCT Extent().
    pub(crate) fn extent(&self) -> i32 {
        self.keys.len() as i32
    }

    /// OCCT FindKey(i) — the 1-based key access.
    pub(crate) fn find_key(&self, i: i32) -> Shape {
        self.keys[(i - 1) as usize].clone()
    }
}

/// OCCT static CheckEdges (cxx L72-135): fills the bad/direct/reversed/
/// internal edge maps; returns true when a bad edge was found.
fn check_edges(
    brep: &mut BRep,
    shape: &Shape,
    bads: &mut IndexedMapShape,
    dirs: &mut IndexedMapShape,
    revs: &mut IndexedMapShape,
    ints: &mut IndexedMapShape,
) -> bool {
    let mut res = false;

    if shape.shape_type() != ShapeType::Edge {
        // OCCT L82: for (TopoDS_Iterator it(shape); it.More(); it.Next())
        // — CumOri = true (the composed orientations).
        for child in iter_subshapes(brep, shape, true, false) {
            if check_edges(brep, &child, bads, dirs, revs, ints) {
                res = true;
            }
        }
    } else {
        let e = shape;
        if brep_tool_degenerated(e) {
            return false;
        }

        if shape.orientation == Orientation::Forward {
            // szv#4:S4163:12Mar99 optimized
            if dirs.find_index(shape) == 0 {
                dirs.add(shape);
            } else {
                bads.add(shape);
                res = true;
            }
        }
        if shape.orientation == Orientation::Reversed {
            // szv#4:S4163:12Mar99 optimized
            if revs.find_index(shape) == 0 {
                revs.add(shape);
            } else {
                bads.add(shape);
                res = true;
            }
        }
        if shape.orientation == Orientation::Internal {
            if ints.find_index(shape) == 0 {
                ints.add(shape);
            }
            // else { bads.Add (shape); res = true; }
        }
    }

    res
}

/// OCCT ShapeAnalysis_Shell (hxx L29-110).
#[derive(Default)]
pub struct ShapeAnalysisShell {
    /// OCCT `myShells`.
    my_shells: IndexedMapShape,
    /// OCCT `myBad`.
    my_bad: IndexedMapShape,
    /// OCCT `myFree`.
    my_free: IndexedMapShape,
    /// OCCT `myConex`.
    my_conex: bool,
}

impl ShapeAnalysisShell {
    /// OCCT ShapeAnalysis_Shell() (cxx L28-31).
    pub fn new() -> Self {
        ShapeAnalysisShell {
            my_shells: IndexedMapShape::default(),
            my_bad: IndexedMapShape::default(),
            my_free: IndexedMapShape::default(),
            my_conex: false,
        }
    }

    /// OCCT Clear() (cxx L35-41).
    pub fn clear(&mut self) {
        self.my_shells.clear();
        self.my_bad.clear();
        self.my_free.clear();
        self.my_conex = false;
    }

    /// OCCT LoadShells(shape) (cxx L45-64).
    pub fn load_shells(&mut self, brep: &mut BRep, shape: &Shape) {
        if shape.is_null() {
            return;
        }

        if shape.shape_type() == ShapeType::Shell {
            self.my_shells.add(shape); // szv#4:S4163:12Mar99 i =
        } else {
            for sh in topexp_explorer(brep, shape, ShapeType::Shell) {
                self.my_shells.add(&sh); // szv#4:S4163:12Mar99 i =
            }
        }
    }

    /// OCCT CheckOrientedShells(shape, alsofree = true, checkinternaledges
    /// = false) (cxx L139-244): feeds BadEdges and FreeEdges.
    pub fn check_oriented_shells(
        &mut self,
        brep: &mut BRep,
        shape: &Shape,
        alsofree: bool,
        checkinternaledges: bool,
    ) -> bool {
        self.my_conex = false;
        if shape.is_null() {
            return false;
        }
        let mut res = false;

        let mut dirs = IndexedMapShape::default();
        let mut revs = IndexedMapShape::default();
        let mut ints = IndexedMapShape::default();
        for sh in topexp_explorer(brep, shape, ShapeType::Shell) {
            // szv#4:S4163:12Mar99 optimized
            if check_edges(brep, &sh, &mut self.my_bad, &mut dirs, &mut revs, &mut ints) {
                if self.my_shells.add(&sh) {
                    res = true;
                }
            }
        }

        //  Resteraient a faire les FreeEdges
        if !alsofree {
            return res;
        }

        //  Free Edges. Ce sont les edges d une map pas dans l autre
        //  et lycee de Versailles  (les maps dirs et revs)
        let mut nb = dirs.extent();
        let mut i = 1;
        while i <= nb {
            let sh = dirs.find_key(i);
            if !self.my_bad.contains(&sh) {
                if !revs.contains(&sh) {
                    if checkinternaledges {
                        if !ints.contains(&sh) {
                            self.my_free.add(&sh);
                        } else {
                            self.my_conex = true;
                        }
                    } else {
                        self.my_free.add(&sh);
                    }
                } else {
                    self.my_conex = true;
                }
            } else {
                self.my_conex = true;
            }
            i += 1;
        }

        nb = revs.extent();
        i = 1;
        while i <= nb {
            let sh = revs.find_key(i);
            if !self.my_bad.contains(&sh) {
                if !dirs.contains(&sh) {
                    if checkinternaledges {
                        if !ints.contains(&sh) {
                            self.my_free.add(&sh);
                        } else {
                            self.my_conex = true;
                        }
                    } else {
                        self.my_free.add(&sh);
                    }
                } else {
                    self.my_conex = true;
                }
            } else {
                self.my_conex = true;
            }
            i += 1;
        }

        res
    }

    /// OCCT IsLoaded(shape) (cxx L248-255).
    pub fn is_loaded(&self, shape: &Shape) -> bool {
        if shape.is_null() {
            return false;
        }
        self.my_shells.contains(shape)
    }

    /// OCCT NbLoaded() (cxx L259-262).
    pub fn nb_loaded(&self) -> i32 {
        self.my_shells.extent()
    }

    /// OCCT Loaded(num) (cxx L266-269).
    pub fn loaded(&self, num: i32) -> Shape {
        self.my_shells.find_key(num)
    }

    /// OCCT HasBadEdges() (cxx L273-276).
    pub fn has_bad_edges(&self) -> bool {
        self.my_bad.extent() > 0
    }

    /// OCCT BadEdges() (cxx L280-291).
    pub fn bad_edges(&self, brep: &mut BRep) -> Shape {
        // OCCT: TopoDS_Compound C; BRep_Builder B; B.MakeCompound(C).
        let mut edges = Vec::new();
        let n = self.my_bad.extent();
        for i in 1..=n {
            edges.push(self.my_bad.find_key(i));
        }
        brep.add_tcompound(edges)
    }

    /// OCCT HasFreeEdges() (cxx L295-298).
    pub fn has_free_edges(&self) -> bool {
        self.my_free.extent() > 0
    }

    /// OCCT FreeEdges() (cxx L302-313).
    pub fn free_edges(&self, brep: &mut BRep) -> Shape {
        let mut edges = Vec::new();
        let n = self.my_free.extent();
        for i in 1..=n {
            edges.push(self.my_free.find_key(i));
        }
        brep.add_tcompound(edges)
    }

    /// OCCT HasConnectedEdges() (cxx L317-320).
    pub fn has_connected_edges(&self) -> bool {
        self.my_conex
    }
}
