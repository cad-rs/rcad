// OCCT LocOpe_LinearForm.hxx L32-83 + LocOpe_LinearForm.lxx L21-53 +
// LocOpe_LinearForm.cxx L39-244 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_LinearForm.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_LinearForm.lxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_LinearForm.cxx
//
// OCCT inheritance chain (LocOpe_LinearForm.hxx L32): none —
// LocOpe_LinearForm is a standalone class in the OCCT 8.0 sources at
// $OCCT_SRC (no LocOpe_GeneratedShape virtuals to translate; see
// loc_ope_prism.rs).
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//    TopTools_ShapeMapHasher> (myMap) — HashMap keyed by (TShape ptr,
//    Location); never iterated.
// 2. BRepSweep_Prism — the sweep engine is the crate::brep_sweep
//    translation (BRepSweepPrism); the former loc_ope_prism.rs GAP carrier
//    is removed.
// 3. BRepTools_Modifier + BRepTools_TrsfModification — re-hosted in
//    loc_ope_prism.rs (pub(crate), same GAP).
// 4. gp_Trsf::SetTranslation maps to Trsf::identity() +
//    set_translation_part.
// 5. gp_Pnt -> glam::DVec3 (myPnt1/myPnt2 are stored but never read by the
//    OCCT IntPerf — the fields are part of the hxx form).
//
// first consumer: BRepFeat_MakeLinearForm (3b).

use crate::brep_sweep::BRepSweepPrism;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_build_shape::LocOpeBuildShape;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::feat::loc_ope_prism::{ BRepToolsModifier, BRepToolsTrsfModification };
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::math::gp::Trsf;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
use std::collections::HashMap;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT LocOpe_LinearForm (LocOpe_LinearForm.hxx L32-83).
pub struct LocOpeLinearForm {
    my_base: Shape,        // OCCT: myBase
    my_vec: DVec3,         // OCCT: myVec (gp_Vec)
    my_tra: DVec3,         // OCCT: myTra (gp_Vec)
    my_done: bool,         // OCCT: myDone
    my_is_trans: bool,     // OCCT: myIsTrans
    my_res: Shape,         // OCCT: myRes
    my_first_shape: Shape, // OCCT: myFirstShape
    my_last_shape: Shape,  // OCCT: myLastShape
    // OCCT: myMap (NCollection_DataMap) — arch. diff. #1
    my_map: HashMap<(u64, u32), Vec<Shape>>,
    my_pnt1: DVec3, // OCCT: myPnt1 (gp_Pnt)
    my_pnt2: DVec3, // OCCT: myPnt2 (gp_Pnt)
}

impl Default for LocOpeLinearForm {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeLinearForm {
    /// OCCT LocOpe_LinearForm::LocOpe_LinearForm() (lxx L21-30).
    pub fn new() -> Self {
        LocOpeLinearForm {
            my_base: Shape::null(),
            my_vec: DVec3::ZERO,
            my_tra: DVec3::ZERO,
            my_done: false,
            my_is_trans: false,
            my_res: Shape::null(),
            my_first_shape: Shape::null(),
            my_last_shape: Shape::null(),
            my_map: HashMap::new(),
            my_pnt1: DVec3::ZERO,
            my_pnt2: DVec3::ZERO,
        }
    }

    /// OCCT LocOpe_LinearForm::LocOpe_LinearForm(Base, V, Pnt1, Pnt2)
    /// (lxx L32-42).
    pub fn new_full(the_base: &Shape, v: DVec3, pnt1: DVec3, pnt2: DVec3) -> Self {
        let mut res = LocOpeLinearForm::new();
        res.perform(the_base, v, pnt1, pnt2);
        res
    }

    /// OCCT LocOpe_LinearForm::LocOpe_LinearForm(Base, V, Vectra, Pnt1,
    /// Pnt2) (lxx L44-53).
    pub fn new_trans(
        the_base: &Shape,
        v: DVec3,
        vectra: DVec3,
        pnt1: DVec3,
        pnt2: DVec3,
    ) -> Self {
        let mut res = LocOpeLinearForm::new();
        res.perform_trans(the_base, v, vectra, pnt1, pnt2);
        res
    }

    /// OCCT LocOpe_LinearForm::Perform(Base, V, Pnt1, Pnt2) (cxx L39-60).
    pub fn perform(&mut self, the_base: &Shape, v: DVec3, pnt1: DVec3, pnt2: DVec3) {
        self.my_is_trans = false;
        self.my_map.clear();
        self.my_first_shape = Shape::null();
        self.my_last_shape = Shape::null();
        self.my_base = Shape::null();
        self.my_res = Shape::null();

        self.my_base = the_base.clone();
        self.my_vec = v;

        // myEdge = E;
        self.my_pnt1 = pnt1;
        self.my_pnt2 = pnt2;

        self.int_perf();
    }

    /// OCCT LocOpe_LinearForm::Perform(Base, V, Vectra, Pnt1, Pnt2)
    /// (cxx L64-87) — Rust has no overloading: the `_trans` suffix.
    pub fn perform_trans(
        &mut self,
        the_base: &Shape,
        v: DVec3,
        vectra: DVec3,
        pnt1: DVec3,
        pnt2: DVec3,
    ) {
        self.my_is_trans = true;
        self.my_tra = vectra;
        self.my_map.clear();
        self.my_first_shape = Shape::null();
        self.my_last_shape = Shape::null();
        self.my_base = Shape::null();
        self.my_res = Shape::null();

        self.my_base = the_base.clone();
        self.my_vec = v;

        // myEdge = E;
        self.my_pnt1 = pnt1;
        self.my_pnt2 = pnt2;

        self.int_perf();
    }

    /// OCCT LocOpe_LinearForm::IntPerf() (cxx L91-212).
    fn int_perf(&mut self) {
        // OCCT cxx L93-94.
        let mut the_base = self.my_base.clone();
        let mut modif = BRepToolsModifier::new();

        if self.my_is_trans {
            // OCCT cxx L98-99 (arch. diff. #4).
            let mut t = Trsf::identity();
            t.set_translation_part(self.my_tra);
            // OCCT cxx L100-103 (arch. diff. #3).
            let modbase = BRepToolsTrsfModification::new(t);
            modif.init(&the_base);
            modif.perform(&modbase);
            the_base = modif.modified_shape(&the_base);
        }

        // OCCT cxx L106 (arch. diff. #2): BRepSweep_Prism thePrism(theBase,
        // myVec) — the OCCT default arguments are C=false, Canonize=true.
        let mut my_prism = BRepSweepPrism::with_vec(&the_base, self.my_vec, false, true);

        // OCCT cxx L108-109.
        self.my_first_shape = my_prism.first_shape();
        self.my_last_shape = my_prism.last_shape();

        // OCCT cxx L111-129: the base-FACE branch.
        if the_base.shape_type() == ShapeType::Face {
            for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L117.
                if !self.my_map.contains_key(&shape_key(&edg)) {
                    // OCCT cxx L119-120.
                    self.my_map.insert(shape_key(&edg), Vec::new());
                    // OCCT cxx L121.
                    let desc = my_prism.shape_of(&edg);
                    // OCCT cxx L122-125: if (!desc.IsNull()) — the engine
                    // null result is the Vertex-typed dummy of Shape::null();
                    // the kernel is_null() (index == usize::MAX) also fires
                    // for pool-built real shapes, so the emptiness test is
                    // the type test (a real generated shape is non-Vertex).
                    if desc.shape_type() != ShapeType::Vertex {
                        self.my_map
                            .get_mut(&shape_key(&edg))
                            .expect("myMap(edg)")
                            .push(desc);
                    }
                }
            }
            // OCCT cxx L128.
            self.my_res = my_prism.shape();
        } else {
            // Cas base != FACE
            //
            // OCCT cxx L134-138: theEFMap.
            let mut the_ef_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            map_shapes_and_ancestors(
                &the_base,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_ef_map,
            );
            // OCCT cxx L139-140.
            let mut lfaces: Vec<Shape> = Vec::new();
            let mut toremove = false;
            // OCCT cxx L141-159.
            for i in 1..=the_ef_map.len() {
                // OCCT cxx L143: edg = theEFMap.FindKey(i).
                let (_, entry) = the_ef_map.get_index(i - 1).expect("theEFMap entry");
                let edg = entry.0.clone();
                // OCCT cxx L144-145.
                self.my_map.insert(shape_key(&edg), Vec::new());
                // OCCT cxx L146.
                let desc = my_prism.shape_of(&edg);
                // OCCT cxx L147-158: if (!desc.IsNull()) — the engine null
                // result is the Vertex-typed dummy of Shape::null(); the
                // kernel is_null() (index == usize::MAX) also fires for
                // pool-built real shapes, so the emptiness test is the type
                // test (a real generated shape is non-Vertex).
                if desc.shape_type() != ShapeType::Vertex {
                    if entry.1.len() >= 2 {
                        toremove = true;
                    } else {
                        self.my_map
                            .get_mut(&shape_key(&edg))
                            .expect("myMap(edg)")
                            .push(desc.clone());
                        lfaces.push(desc);
                    }
                }
            }
            // OCCT cxx L160-174.
            if toremove {
                // Rajouter les faces de FirstShape et LastShape
                for f in explorer(&self.my_first_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }
                for f in explorer(&self.my_last_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }

                // OCCT cxx L172-173.
                let bs = LocOpeBuildShape::with_faces(&lfaces);
                self.my_res = bs.shape().cloned().unwrap_or_else(Shape::null);
            } else {
                // OCCT cxx L177-191.
                for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT cxx L180.
                    if !self.my_map.contains_key(&shape_key(&edg)) {
                        // OCCT cxx L182-183.
                        self.my_map.insert(shape_key(&edg), Vec::new());
                        // OCCT cxx L184.
                        let desc = my_prism.shape_of(&edg);
                        // OCCT cxx L185-188: if (!desc.IsNull()) — the engine
                        // null result is the Vertex-typed dummy of
                        // Shape::null(); the kernel is_null()
                        // (index == usize::MAX) also fires for pool-built
                        // real shapes, so the emptiness test is the type
                        // test (a real generated shape is non-Vertex).
                        if desc.shape_type() != ShapeType::Vertex {
                            self.my_map
                                .get_mut(&shape_key(&edg))
                                .expect("myMap(edg)")
                                .push(desc);
                        }
                    }
                }
                // OCCT cxx L191.
                self.my_res = my_prism.shape();
            }
        }

        // OCCT cxx L195-209: m-a-j des descendants.
        if self.my_is_trans {
            for an_edg in explorer(&self.my_base.clone(), ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L201-202.
                let edg = an_edg.clone();
                let edgbis = modif.modified_shape(&edg);
                // OCCT cxx L203-207.
                if !edgbis.is_same(&edg) && self.my_map.contains_key(&shape_key(&edgbis)) {
                    let list = self
                        .my_map
                        .get(&shape_key(&edgbis))
                        .expect("myMap(edgbis)")
                        .clone();
                    self.my_map.insert(shape_key(&edg), list);
                    self.my_map.remove(&shape_key(&edgbis));
                }
            }
        }

        // OCCT cxx L211.
        self.my_done = true;
    }

    /// OCCT LocOpe_LinearForm::Shape() (cxx L216-223).
    pub fn shape(&self) -> Shape {
        if !self.my_done {
            // OCCT cxx L218-221: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_res.clone()
    }

    /// OCCT LocOpe_LinearForm::FirstShape() (cxx L227-230).
    pub fn first_shape(&self) -> Shape {
        self.my_first_shape.clone()
    }

    /// OCCT LocOpe_LinearForm::LastShape() (cxx L234-237).
    pub fn last_shape(&self) -> Shape {
        self.my_last_shape.clone()
    }

    /// OCCT LocOpe_LinearForm::Shapes(S) (cxx L241-244).
    pub fn shapes(&self, s: &Shape) -> &Vec<Shape> {
        // OCCT cxx L243: return myMap(S) — operator() asserts IsBound.
        self.my_map.get(&shape_key(s)).expect("Standard_NoSuchObject")
    }
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
