//! OCCT ShapeExtend_Explorer (TKShHealing): `.hxx` L37-107 and `.cxx`
//! L28-405 — explores shapes and converts different representations (list,
//! sequence, compound) of complex shapes.
//!
//! Architecture mapping:
//! - `handle(NCollection_HSequence(TopoDS_Shape))` ->
//!   `Option<Vec<Shape>>` (the null handle maps to `None`).
//! - `NCollection_List(TopoDS_Shape)` -> `Vec<Shape>`.
//! - `TopoDS_Iterator` -> [`crate::shhealing::shape_build::brep_tool::iter_subshapes`]
//!   and `TopExp_Explorer` ->
//!   [`crate::shhealing::shape_build::brep_tool::topexp_explorer`] — the
//!   W1-3 dedup: the Explorer's walks delegate to the single TopExp/TopoDS
//!   primitive re-hosts in `shape_build/brep_tool.rs` instead of duplicating
//!   them.
//! - `BRep_Builder::MakeCompound/Add` -> the kernel pool +
//!   [`crate::shhealing::shape_build::brep_tool::builder_add`], so every
//!   building method carries the `&mut BRep` pool parameter.

use crate::shhealing::shape_build::brep_tool::{
    builder_add, iter_subshapes, topexp_explorer,
};
use rcad_kernel::topo::topods::{BRep, Shape, ShapeType};

/// OCCT ShapeExtend_Explorer (ShapeExtend_Explorer.hxx L37-107).
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeExtendExplorer;

impl ShapeExtendExplorer {
    /// OCCT ShapeExtend_Explorer() (cxx L28): creates an object Explorer.
    pub fn new() -> Self {
        ShapeExtendExplorer
    }

    /// OCCT CompoundFromSeq(seqval) (cxx L32-44): converts a sequence of
    /// Shapes to a Compound.
    pub fn compound_from_seq(&self, brep: &mut BRep, seqval: &[Shape]) -> Shape {
        // OCCT L35-37: BRep_Builder B; TopoDS_Compound C; B.MakeCompound(C).
        let c = brep.add_tcompound(Vec::new());
        let n = seqval.len() as i32;
        for i in 1..=n {
            // OCCT L41: B.Add(C, seqval->Value(i)).
            builder_add(brep, &c, &seqval[(i - 1) as usize]);
        }
        c
    }

    /// OCCT SeqFromCompound(comp, expcomp) (cxx L48-86): converts a Compound
    /// to a sequence of Shapes; when `comp` is not a compound the sequence
    /// contains only `comp`; when `comp` is Null the sequence is empty; when
    /// `expcomp` is True nested compounds are flattened recursively.
    pub fn seq_from_compound(&self, brep: &mut BRep, comp: &Shape, expcomp: bool) -> Vec<Shape> {
        // OCCT L74-78: list = new HSequence; if comp is Null return it.
        let mut list: Vec<Shape> = Vec::new();
        if comp.is_null() {
            return list;
        }
        if comp.shape_type() != ShapeType::Compound {
            list.push(comp.clone());
            return list;
        }
        fill_list(brep, comp, expcomp, &mut list);
        list
    }

    /// OCCT ListFromSeq(seqval, lisval, clear) (cxx L90-108): converts a
    /// Sequence of Shapes to a List of Shapes; `clear` (True default)
    /// commands the list to start from scratch, else the list is cumulated.
    pub fn list_from_seq(&self, seqval: Option<&[Shape]>, lisval: &mut Vec<Shape>, clear: bool) {
        if clear {
            lisval.clear();
        }
        let Some(seqval) = seqval else {
            return;
        };
        let nb = seqval.len() as i32;
        for i in 1..=nb {
            lisval.push(seqval[(i - 1) as usize].clone());
        }
    }

    /// OCCT SeqFromList(lisval) (cxx L112-123): converts a List of Shapes to
    /// a Sequence of Shapes.
    pub fn seq_from_list(&self, lisval: &[Shape]) -> Vec<Shape> {
        let mut seqval: Vec<Shape> = Vec::new();
        for it in lisval {
            seqval.push(it.clone());
        }
        seqval
    }

    /// OCCT ShapeType(shape, compound) (cxx L127-179): returns the type of a
    /// Shape — the true type when `compound` is False.  When `compound` is
    /// True and `shape` is a Compound, iterates on its items; if all are of
    /// the same type returns this type, else returns COMPOUND (an empty or
    /// Null shape returns SHAPE).  Equality tolerates the EDGE/WIRE and
    /// FACE/SHELL pseudo-equalities.
    pub fn shape_type(&self, brep: &mut BRep, shape: &Shape, compound: bool) -> ShapeType {
        if shape.is_null() {
            return ShapeType::Shape;
        }
        let mut res = shape.shape_type();
        if !compound || res != ShapeType::Compound {
            return res;
        }
        res = ShapeType::Shape;
        for sh in iter_subshapes(brep, shape, true, true) {
            if sh.is_null() {
                continue;
            }
            let mut typ = sh.shape_type();
            if typ == ShapeType::Compound {
                typ = self.shape_type(brep, &sh, compound);
            }
            if res == ShapeType::Shape {
                res = typ;
                // Egalite : OK;  Pseudo-Egalite : EDGE/WIRE ou FACE/SHELL
            } else if res == ShapeType::Edge && typ == ShapeType::Wire {
                res = typ;
            } else if res == ShapeType::Wire && typ == ShapeType::Edge {
                continue;
            } else if res == ShapeType::Face && typ == ShapeType::Shell {
                res = typ;
            } else if res == ShapeType::Shell && typ == ShapeType::Face {
                continue;
            } else if res != typ {
                return ShapeType::Compound;
            }
        }
        res
    }

    /// OCCT SortedCompound(shape, type, explore, compound) (cxx L183-315):
    /// builds a COMPOUND from the given shape, exploring it level by level
    /// according to `explore`; free edges become wires when `type` is WIRE
    /// and free faces become shells when `type` is SHELL; when `compound` is
    /// True items are gathered in compounds corresponding to the starting
    /// COMPOUND/SOLID/SHELL containers.
    pub fn sorted_compound(
        &self,
        brep: &mut BRep,
        shape: &Shape,
        typ: ShapeType,
        explore: bool,
        compound: bool,
    ) -> Shape {
        if shape.is_null() {
            return shape.clone();
        }
        let mut typ_sh = shape.shape_type();
        let mut sh = Shape::null();
        let mut sh0;
        let mut nb = 0;

        // Compound : on le prend, soit tel quel, soit son contenu
        // (cxx L197-235).
        if typ_sh == ShapeType::Compound || typ_sh == ShapeType::CompSolid {
        // OCCT L199-201: TopoDS_Compound C; BRep_Builder B; MakeCompound.
        let c = brep.add_tcompound(Vec::new());
            for it in iter_subshapes(brep, shape, true, true) {
                sh0 = self.sorted_compound(brep, &it, typ, explore, compound);
                if sh0.is_null() {
                    continue;
                }
                sh = sh0;
                typ_sh = sh.shape_type();
                if typ_sh == ShapeType::Compound && !compound {
                    for it2 in iter_subshapes(brep, &sh, true, true) {
                        nb += 1;
                        sh = it2;
                        builder_add(brep, &c, &sh);
                    }
                } else {
                    nb += 1;
                    builder_add(brep, &c, &sh);
                }
            }
            if nb == 0 {
                // OCCT L226-229: C.Nullify().
                return Shape::null();
            } else if nb == 1 {
                return sh;
            }
            return c;
        }

        // Egalite : OK;  Pseudo-Egalite : EDGE/WIRE ou FACE/SHELL
        // (cxx L237-257).
        if typ_sh == typ {
            return shape.clone();
        }
        if typ_sh == ShapeType::Edge && typ == ShapeType::Wire {
            let w = brep.add_twire(Vec::new());
            builder_add(brep, &w, shape);
            return w;
        }
        if typ_sh == ShapeType::Face && typ == ShapeType::Shell {
            let s = brep.add_tshell(Vec::new());
            builder_add(brep, &s, shape);
            return s;
        }
        // Le reste : selon exploration (cxx L259-264).
        if !explore {
            return Shape::null();
        }

        // Ici, on doit explorer
        // SOLID + mode COMPOUND : reconduire les SHELLs (cxx L267-293).
        if typ_sh == ShapeType::Solid && compound {
            let c = brep.add_tcompound(Vec::new());
            for it in iter_subshapes(brep, shape, true, true) {
                sh0 = self.sorted_compound(brep, &it, typ, explore, compound);
                if sh0.is_null() {
                    continue;
                }
                sh = sh0;
                nb += 1;
                builder_add(brep, &c, &sh);
            }
            if nb == 0 {
                return Shape::null();
            } else if nb == 1 {
                return sh;
            }
            return c;
        }

        // Exploration classique (cxx L295-315): TopExp_Explorer walk that
        // delegates to the brep_tool primitive (W1-3 dedup).
        let cc = brep.add_tcompound(Vec::new());
        for a_exp in topexp_explorer(brep, shape, typ) {
            nb += 1;
            sh = a_exp;
            builder_add(brep, &cc, &sh);
        }
        if nb == 0 {
            return Shape::null();
        } else if nb == 1 {
            return sh;
        }
        cc
    }

    /// OCCT DispatchList(list, vertices, edges, wires, faces, shells, solids,
    /// compsols, compounds) (cxx L319-405): dispatches the starting list of
    /// shapes according to their type to the appropriate resulting lists; a
    /// null result list is firstly created, else new items are appended.
    pub fn dispatch_list(
        &self,
        list: Option<&[Shape]>,
        vertices: &mut Option<Vec<Shape>>,
        edges: &mut Option<Vec<Shape>>,
        wires: &mut Option<Vec<Shape>>,
        faces: &mut Option<Vec<Shape>>,
        shells: &mut Option<Vec<Shape>>,
        solids: &mut Option<Vec<Shape>>,
        compsols: &mut Option<Vec<Shape>>,
        compounds: &mut Option<Vec<Shape>>,
    ) {
        let Some(list) = list else {
            return;
        };
        if vertices.is_none() {
            *vertices = Some(Vec::new());
        }
        if edges.is_none() {
            *edges = Some(Vec::new());
        }
        if wires.is_none() {
            *wires = Some(Vec::new());
        }
        if faces.is_none() {
            *faces = Some(Vec::new());
        }
        if shells.is_none() {
            *shells = Some(Vec::new());
        }
        if solids.is_none() {
            *solids = Some(Vec::new());
        }
        if compsols.is_none() {
            *compsols = Some(Vec::new());
        }
        if compounds.is_none() {
            *compounds = Some(Vec::new());
        }

        let nb = list.len() as i32;
        for i in 1..=nb {
            let sh = list[(i - 1) as usize].clone();
            if sh.is_null() {
                continue;
            }
            match sh.shape_type() {
                ShapeType::Vertex => vertices.as_mut().unwrap().push(sh),
                ShapeType::Edge => edges.as_mut().unwrap().push(sh),
                ShapeType::Wire => wires.as_mut().unwrap().push(sh),
                ShapeType::Face => faces.as_mut().unwrap().push(sh),
                ShapeType::Shell => shells.as_mut().unwrap().push(sh),
                ShapeType::Solid => solids.as_mut().unwrap().push(sh),
                ShapeType::CompSolid => compsols.as_mut().unwrap().push(sh),
                ShapeType::Compound => compounds.as_mut().unwrap().push(sh),
                _ => {}
            }
        }
    }
}

/// OCCT FillList (cxx L48-68): appends the compound items to the list,
/// recursing into nested compounds when `expcomp` is True.
fn fill_list(brep: &mut BRep, comp: &Shape, expcomp: bool, list: &mut Vec<Shape>) {
    for sub in iter_subshapes(brep, comp, true, true) {
        if sub.shape_type() != ShapeType::Compound {
            list.push(sub);
        } else if !expcomp {
            list.push(sub);
        } else {
            fill_list(brep, &sub, expcomp, list);
        }
    }
}
