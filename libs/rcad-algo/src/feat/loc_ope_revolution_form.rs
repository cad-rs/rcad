// OCCT LocOpe_RevolutionForm.hxx L33-71 + LocOpe_RevolutionForm.cxx
// L38-218 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_RevolutionForm.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_RevolutionForm.cxx
//
// OCCT inheritance chain (LocOpe_RevolutionForm.hxx L33): none —
// LocOpe_RevolutionForm is a standalone class in the OCCT 8.0 sources at
// $OCCT_SRC (no LocOpe_GeneratedShape virtuals to translate; see
// loc_ope_prism.rs).
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//    TopTools_ShapeMapHasher> (myMap) — HashMap keyed by (TShape ptr,
//    Location); never iterated.
// 2. BRepSweep_Revol — the sweep engine has no rcad equivalent yet; the
//    BRepSweepRevol carrier of loc_ope_revol.rs is reused (same GAP).
// 3. BRepTools_Modifier + BRepTools_TrsfModification — re-hosted in
//    loc_ope_prism.rs (pub(crate), same GAP).
// 4. gp_Trsf::SetRotation(Ax, Ang) — the trsf_set_rotation GAP helper of
//    loc_ope_revol.rs is reused (arch. diff. #4 there).
// 5. gp_Pnt -> glam::DVec3 (myPnt1/myPnt2/myVec/myTra are hxx-form fields
//    never read by the OCCT IntPerf).
//
// The value constructor declared in the OCCT header (hxx L40-42) has no
// definition in LocOpe_RevolutionForm.cxx — nothing to translate.
//
// first consumer: BRepFeat_MakeRevolutionForm (3b).

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_build_shape::LocOpeBuildShape;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::feat::loc_ope_prism::{ BRepToolsModifier, BRepToolsTrsfModification };
use crate::feat::loc_ope_revol::{ trsf_set_rotation, BRepSweepRevol };
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::math::gp::{ Ax1, Trsf };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
use std::collections::HashMap;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT LocOpe_RevolutionForm (LocOpe_RevolutionForm.hxx L33-71).
pub struct LocOpeRevolutionForm {
    my_base: Shape,        // OCCT: myBase
    my_vec: DVec3,         // OCCT: myVec (gp_Vec)
    my_tra: DVec3,         // OCCT: myTra (gp_Vec)
    my_angle: f64,         // OCCT: myAngle
    my_axis: Ax1,          // OCCT: myAxis (gp_Ax1)
    my_ang_tra: f64,       // OCCT: myAngTra
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

impl Default for LocOpeRevolutionForm {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpeRevolutionForm {
    /// OCCT LocOpe_RevolutionForm::LocOpe_RevolutionForm() (cxx L38-44).
    pub fn new() -> Self {
        LocOpeRevolutionForm {
            my_base: Shape::null(),
            my_vec: DVec3::ZERO,
            my_tra: DVec3::ZERO,
            my_angle: 0.0,
            my_axis: Ax1::new(DVec3::ZERO, DVec3::Z),
            my_ang_tra: 0.0,
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

    /// OCCT LocOpe_RevolutionForm::Perform(Base, Axis, Angle)
    /// (cxx L48-63).
    pub fn perform(&mut self, the_base: &Shape, the_axis: &Ax1, the_angle: f64) {
        self.my_map.clear();
        self.my_first_shape = Shape::null();
        self.my_last_shape = Shape::null();
        self.my_base = Shape::null();
        self.my_res = Shape::null();
        self.my_base = the_base.clone();
        self.my_angle = the_angle;
        self.my_axis = *the_axis;
        self.my_ang_tra = 0.;
        self.my_is_trans = false;
        self.int_perf();
    }

    /// OCCT LocOpe_RevolutionForm::IntPerf() (cxx L67-186).
    fn int_perf(&mut self) {
        // OCCT cxx L69-70.
        let mut the_base = self.my_base.clone();
        let mut modif = BRepToolsModifier::new();
        if self.my_is_trans {
            // OCCT cxx L73-74 (arch. diff. #4).
            let mut t = Trsf::identity();
            trsf_set_rotation(&mut t, &self.my_axis, self.my_ang_tra);
            // OCCT cxx L75-78 (arch. diff. #3).
            let modbase = BRepToolsTrsfModification::new(t);
            modif.init(&the_base);
            modif.perform(&modbase);
            the_base = modif.modified_shape(&the_base);
        }

        // OCCT cxx L81 (arch. diff. #2).
        let the_revol = BRepSweepRevol::new(&the_base, &self.my_axis, self.my_angle);

        // OCCT cxx L83-84.
        self.my_first_shape = the_revol.first_shape();
        self.my_last_shape = the_revol.last_shape();

        // OCCT cxx L86-104: the base-FACE branch.
        if the_base.shape_type() == ShapeType::Face {
            for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L92.
                if !self.my_map.contains_key(&shape_key(&edg)) {
                    // OCCT cxx L94-95.
                    self.my_map.insert(shape_key(&edg), Vec::new());
                    // OCCT cxx L96.
                    let desc = the_revol.shape_of_edge(&edg);
                    // OCCT cxx L97-100.
                    if !desc.is_null() {
                        self.my_map
                            .get_mut(&shape_key(&edg))
                            .expect("myMap(edg)")
                            .push(desc);
                    }
                }
            }
            // OCCT cxx L103.
            self.my_res = the_revol.shape();
        } else {
            // Cas base != FACE
            //
            // OCCT cxx L109-113: theEFMap.
            let mut the_ef_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            map_shapes_and_ancestors(
                &the_base,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_ef_map,
            );
            // OCCT cxx L114-115.
            let mut lfaces: Vec<Shape> = Vec::new();
            let mut toremove = false;
            // OCCT cxx L116-134.
            for i in 1..=the_ef_map.len() {
                // OCCT cxx L118: edg = theEFMap.FindKey(i).
                let (_, entry) = the_ef_map.get_index(i - 1).expect("theEFMap entry");
                let edg = entry.0.clone();
                // OCCT cxx L119-120.
                self.my_map.insert(shape_key(&edg), Vec::new());
                // OCCT cxx L121.
                let desc = the_revol.shape_of_edge(&edg);
                // OCCT cxx L122-133.
                if !desc.is_null() {
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
            // OCCT cxx L135-149.
            if toremove {
                // Rajouter les faces de FirstShape et LastShape
                for f in explorer(&self.my_first_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }
                for f in explorer(&self.my_last_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }

                // OCCT cxx L147-148.
                let bs = LocOpeBuildShape::with_faces(&lfaces);
                self.my_res = bs.shape().cloned().unwrap_or_else(Shape::null);
            } else {
                // OCCT cxx L152-166.
                for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT cxx L155.
                    if !self.my_map.contains_key(&shape_key(&edg)) {
                        // OCCT cxx L157-158.
                        self.my_map.insert(shape_key(&edg), Vec::new());
                        // OCCT cxx L159.
                        let desc = the_revol.shape_of_edge(&edg);
                        // OCCT cxx L160-163.
                        if !desc.is_null() {
                            self.my_map
                                .get_mut(&shape_key(&edg))
                                .expect("myMap(edg)")
                                .push(desc);
                        }
                    }
                }
                // OCCT cxx L166.
                self.my_res = the_revol.shape();
            }
        }

        // OCCT cxx L170-184: m-a-j des descendants.
        if self.my_is_trans {
            for an_edg in explorer(&self.my_base.clone(), ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L175-176.
                let edg = an_edg.clone();
                let edgbis = modif.modified_shape(&edg);
                // OCCT cxx L177-181.
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
        // OCCT cxx L185.
        self.my_done = true;
    }

    /// OCCT LocOpe_RevolutionForm::Shape() (cxx L190-197).
    pub fn shape(&self) -> Shape {
        if !self.my_done {
            // OCCT cxx L192-195: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_res.clone()
    }

    /// OCCT LocOpe_RevolutionForm::FirstShape() (cxx L201-204).
    pub fn first_shape(&self) -> Shape {
        self.my_first_shape.clone()
    }

    /// OCCT LocOpe_RevolutionForm::LastShape() (cxx L208-211).
    pub fn last_shape(&self) -> Shape {
        self.my_last_shape.clone()
    }

    /// OCCT LocOpe_RevolutionForm::Shapes(S) (cxx L215-218).
    pub fn shapes(&self, s: &Shape) -> &Vec<Shape> {
        // OCCT cxx L217: return myMap(S) — operator() asserts IsBound.
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
