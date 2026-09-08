// OCCT LocOpe_Prism.hxx L35-74 + LocOpe_Prism.cxx L44-297 — 1:1
// translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Prism.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/LocOpe/LocOpe_Prism.cxx
//
// OCCT inheritance chain (LocOpe_Prism.hxx L35): none — LocOpe_Prism is a
// standalone class in the OCCT 8.0 sources at $OCCT_SRC. The pre-8.0
// LocOpe_GeneratedShape base has no counterpart in this source line, so no
// LocOpeGeneratedShape impl exists here to translate.
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>,
//    TopTools_ShapeMapHasher> (myMap) — HashMap keyed by (TShape ptr,
//    Location); the map is never iterated, only Bound/UnBound/operator()
//    accessed.
// 2. BRepSweep_Prism (TKTopAlgo/BRepSweep) — the sweep engine has no rcad
//    equivalent yet; the BRepSweepPrism carrier below carries the OCCT
//    constructor/accessor surface with a GAP panic (port plan §0.6 gap
//    annotation). Everything downstream of the construction (IntPerf) is
//    translated 1:1 against the carrier.
// 3. BRepTools_Modifier + BRepTools_TrsfModification (TKTopAlgo/BRepTools)
//    — the myIsTrans branch vehicle; the BRepToolsModifier /
//    BRepToolsTrsfModification carriers below carry the same GAP (they are
//    shared with loc_ope_revol / loc_ope_linear_form /
//    loc_ope_revolution_form).
// 4. gp_Trsf::SetTranslation maps to Trsf::identity() +
//    set_translation_part (rcad_kernel::math::gp).
// 5. NCollection_Sequence<gp_Pnt> is Vec<DVec3>; NCollection_Sequence<
//    Handle(Geom_Curve)> is Vec<Curve3> — the out-parameter form of
//    Curves() is kept.
// 6. Geom_Line -> geom::Line3 (Curve3::Line), Geom_TrimmedCurve ->
//    TrimmedCurve3 (the OCCT Sense=true argument has no rcad counterpart:
//    TrimmedCurve3 keeps the direction).
//
// first consumer: BRepFeat_MakePrism (3b) — Perform/IntPerf drive the prism
// engine and read Shape/FirstShape/LastShape/Shapes/Curves/BarycCurve.

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::feat::loc_ope_build_shape::LocOpeBuildShape;
use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::geom::{ Curve3, Line3, TrimmedCurve3 };
use rcad_kernel::math::gp::{ Ax1, Trsf };
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
use std::collections::HashMap;

/// OCCT BRepTools_TrsfModification (BRepTools_TrsfModification.hxx) — the
/// gp_Trsf carrier (architecture difference #3; GAP: no rcad translation of
/// the BRepTools modification framework yet).
pub(crate) struct BRepToolsTrsfModification {
    my_trsf: Trsf,
}

impl BRepToolsTrsfModification {
    /// OCCT BRepTools_TrsfModification(T).
    pub(crate) fn new(t: Trsf) -> Self {
        BRepToolsTrsfModification { my_trsf: t }
    }

    /// OCCT BRepTools_TrsfModification::Trsf() — the carried transformation.
    pub(crate) fn trsf(&self) -> &Trsf {
        &self.my_trsf
    }
}

/// OCCT BRepTools_Modifier (BRepTools_Modifier.hxx) — the shape rebuild
/// vehicle of the myIsTrans branches (architecture difference #3; GAP: the
/// Perform/ModifiedShape engine has no rcad translation yet — the GAP panic
/// is the §0.6 annotation; Init keeps the OCCT storage form).
pub(crate) struct BRepToolsModifier {
    my_shape: Shape, // OCCT: myShape
}

impl BRepToolsModifier {
    /// OCCT BRepTools_Modifier::BRepTools_Modifier().
    pub(crate) fn new() -> Self {
        BRepToolsModifier {
            my_shape: Shape::null(),
        }
    }

    /// OCCT BRepTools_Modifier::Init(S).
    pub(crate) fn init(&mut self, the_shape: &Shape) {
        self.my_shape = the_shape.clone();
    }

    /// OCCT BRepTools_Modifier::Perform(M).
    pub(crate) fn perform(&mut self, m: &BRepToolsTrsfModification) {
        let _ = m;
        panic!("GAP: BRepTools_Modifier::Perform (TKTopAlgo/BRepTools not translated)");
    }

    /// OCCT BRepTools_Modifier::ModifiedShape(S).
    pub(crate) fn modified_shape(&self, the_shape: &Shape) -> Shape {
        let _ = the_shape;
        panic!("GAP: BRepTools_Modifier::ModifiedShape (TKTopAlgo/BRepTools not translated)");
    }
}

/// OCCT BRepSweep_Prism (BRepSweep_Prism.hxx) — the prism sweep engine
/// (architecture difference #2; GAP: BRepSweep has no rcad translation yet
/// — the GAP panic is the §0.6 annotation; the LocOpe_Prism::IntPerf body
/// below is translated 1:1 against this surface).
pub(crate) struct BRepSweepPrism;

impl BRepSweepPrism {
    /// OCCT BRepSweep_Prism::BRepSweep_Prism(S, V).
    pub(crate) fn new(_the_base: &Shape, _the_vec: DVec3) -> Self {
        panic!("GAP: BRepSweep_Prism (TKTopAlgo/BRepSweep not translated)");
    }

    /// OCCT BRepSweep_Prism::FirstShape().
    pub(crate) fn first_shape(&self) -> Shape {
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)");
    }

    /// OCCT BRepSweep_Prism::LastShape().
    pub(crate) fn last_shape(&self) -> Shape {
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)");
    }

    /// OCCT BRepSweep_Prism::Shape() — the whole sweep result.
    pub(crate) fn shape(&self) -> Shape {
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)");
    }

    /// OCCT BRepSweep_Prism::Shape(E) — the generated shape of an ancestor.
    pub(crate) fn shape_of_edge(&self, the_e: &Shape) -> Shape {
        let _ = the_e;
        panic!("GAP: BRepSweep_Prism (unreachable while the engine is a GAP)");
    }
}

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT LocOpe_Prism (LocOpe_Prism.hxx L35-74).
pub struct LocOpePrism {
    my_base: Shape,        // OCCT: myBase
    my_vec: DVec3,         // OCCT: myVec (gp_Vec)
    my_tra: DVec3,         // OCCT: myTra (gp_Vec)
    my_is_trans: bool,     // OCCT: myIsTrans
    my_done: bool,         // OCCT: myDone
    my_res: Shape,         // OCCT: myRes
    my_first_shape: Shape, // OCCT: myFirstShape
    my_last_shape: Shape,  // OCCT: myLastShape
    // OCCT: myMap (NCollection_DataMap) — arch. diff. #1
    my_map: HashMap<(u64, u32), Vec<Shape>>,
}

impl Default for LocOpePrism {
    fn default() -> Self {
        Self::new()
    }
}

impl LocOpePrism {
    /// OCCT LocOpe_Prism::LocOpe_Prism() (cxx L44-48).
    pub fn new() -> Self {
        LocOpePrism {
            my_base: Shape::null(),
            my_vec: DVec3::ZERO,
            my_tra: DVec3::ZERO,
            my_is_trans: false,
            my_done: false,
            my_res: Shape::null(),
            my_first_shape: Shape::null(),
            my_last_shape: Shape::null(),
            my_map: HashMap::new(),
        }
    }

    /// OCCT LocOpe_Prism::LocOpe_Prism(Base, V) (cxx L52-59) — Rust has no
    /// overloading: the `_vec` suffix.
    pub fn with_vec(the_base: &Shape, v: DVec3) -> Self {
        let mut res = LocOpePrism {
            my_base: the_base.clone(),
            my_vec: v,
            my_tra: DVec3::ZERO,
            my_is_trans: false,
            my_done: false,
            my_res: Shape::null(),
            my_first_shape: Shape::null(),
            my_last_shape: Shape::null(),
            my_map: HashMap::new(),
        };
        res.int_perf();
        res
    }

    /// OCCT LocOpe_Prism::LocOpe_Prism(Base, V, Vtra) (cxx L63-71) — Rust
    /// has no overloading: the `_vec_trans` suffix.
    pub fn with_vec_trans(the_base: &Shape, v: DVec3, v_tra: DVec3) -> Self {
        let mut res = LocOpePrism {
            my_base: the_base.clone(),
            my_vec: v,
            my_tra: v_tra,
            my_is_trans: true,
            my_done: false,
            my_res: Shape::null(),
            my_first_shape: Shape::null(),
            my_last_shape: Shape::null(),
            my_map: HashMap::new(),
        };
        res.int_perf();
        res
    }

    /// OCCT LocOpe_Prism::Perform(Base, V) (cxx L75-87).
    pub fn perform(&mut self, the_base: &Shape, v: DVec3) {
        self.my_map.clear();
        self.my_first_shape = Shape::null();
        self.my_last_shape = Shape::null();
        self.my_base = Shape::null();
        self.my_res = Shape::null();

        self.my_base = the_base.clone();
        self.my_vec = v;
        self.my_is_trans = false;
        self.int_perf();
    }

    /// OCCT LocOpe_Prism::Perform(Base, V, Vtra) (cxx L91-104).
    pub fn perform_trans(&mut self, the_base: &Shape, v: DVec3, v_tra: DVec3) {
        self.my_map.clear();
        self.my_first_shape = Shape::null();
        self.my_last_shape = Shape::null();
        self.my_base = Shape::null();
        self.my_res = Shape::null();

        self.my_base = the_base.clone();
        self.my_vec = v;
        self.my_tra = v_tra;
        self.my_is_trans = true;
        self.int_perf();
    }

    /// OCCT LocOpe_Prism::IntPerf() (cxx L108-227).
    fn int_perf(&mut self) {
        // OCCT cxx L110-111.
        let mut the_base = self.my_base.clone();
        let mut modif = BRepToolsModifier::new();
        if self.my_is_trans {
            // OCCT cxx L114-115 (arch. diff. #4).
            let mut t = Trsf::identity();
            t.set_translation_part(self.my_tra);
            // OCCT cxx L116-119 (arch. diff. #3).
            let modbase = BRepToolsTrsfModification::new(t);
            modif.init(&the_base);
            modif.perform(&modbase);
            the_base = modif.modified_shape(&the_base);
        }

        // OCCT cxx L122 (arch. diff. #2).
        let the_prism = BRepSweepPrism::new(&the_base, self.my_vec);

        // OCCT cxx L124-125.
        self.my_first_shape = the_prism.first_shape();
        self.my_last_shape = the_prism.last_shape();

        // OCCT cxx L127-145: the base-FACE branch.
        if the_base.shape_type() == ShapeType::Face {
            for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L133: if (!myMap.IsBound(edg)).
                if !self.my_map.contains_key(&shape_key(&edg)) {
                    // OCCT cxx L135-136.
                    self.my_map.insert(shape_key(&edg), Vec::new());
                    // OCCT cxx L137.
                    let desc = the_prism.shape_of_edge(&edg);
                    // OCCT cxx L138-141.
                    if !desc.is_null() {
                        self.my_map
                            .get_mut(&shape_key(&edg))
                            .expect("myMap(edg)")
                            .push(desc);
                    }
                }
            }
            // OCCT cxx L144.
            self.my_res = the_prism.shape();
        } else {
            // Cas base != FACE
            //
            // OCCT cxx L150-154: theEFMap (arch. diff. #1 vehicle: the
            // pub(crate) TopExp::MapShapesAndAncestors re-host).
            let mut the_ef_map: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
            map_shapes_and_ancestors(
                &the_base,
                ShapeType::Edge,
                ShapeType::Face,
                &mut the_ef_map,
            );
            // OCCT cxx L155-156.
            let mut lfaces: Vec<Shape> = Vec::new();
            let mut toremove = false;
            // OCCT cxx L157-175.
            for i in 1..=the_ef_map.len() {
                // OCCT cxx L159: edg = theEFMap.FindKey(i).
                let (_, entry) = the_ef_map.get_index(i - 1).expect("theEFMap entry");
                let edg = entry.0.clone();
                // OCCT cxx L160-161.
                self.my_map.insert(shape_key(&edg), Vec::new());
                // OCCT cxx L162.
                let desc = the_prism.shape_of_edge(&edg);
                // OCCT cxx L163-174.
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
            // OCCT cxx L176-190.
            if toremove {
                // Rajouter les faces de FirstShape et LastShape
                for f in explorer(&self.my_first_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }
                for f in explorer(&self.my_last_shape.clone(), ShapeType::Face, ShapeType::Shape) {
                    lfaces.push(f);
                }

                // OCCT cxx L188-189.
                let bs = LocOpeBuildShape::with_faces(&lfaces);
                self.my_res = bs.shape().cloned().unwrap_or_else(Shape::null);
            } else {
                // OCCT cxx L193-206.
                for edg in explorer(&the_base, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT cxx L196.
                    if !self.my_map.contains_key(&shape_key(&edg)) {
                        // OCCT cxx L198-199.
                        self.my_map.insert(shape_key(&edg), Vec::new());
                        // OCCT cxx L200.
                        let desc = the_prism.shape_of_edge(&edg);
                        // OCCT cxx L201-204.
                        if !desc.is_null() {
                            self.my_map
                                .get_mut(&shape_key(&edg))
                                .expect("myMap(edg)")
                                .push(desc);
                        }
                    }
                }
                // OCCT cxx L207.
                self.my_res = the_prism.shape();
            }
        }

        // OCCT cxx L211-225: m-a-j des descendants.
        if self.my_is_trans {
            for an_edg in explorer(&self.my_base.clone(), ShapeType::Edge, ShapeType::Shape) {
                // OCCT cxx L217-218.
                let edg = an_edg.clone();
                let edgbis = modif.modified_shape(&edg);
                // OCCT cxx L219-223.
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
        // OCCT cxx L226.
        self.my_done = true;
    }

    /// OCCT LocOpe_Prism::Shape() (cxx L231-238).
    pub fn shape(&self) -> Shape {
        if !self.my_done {
            // OCCT cxx L233-236: throw StdFail_NotDone().
            panic!("StdFail_NotDone");
        }
        self.my_res.clone()
    }

    /// OCCT LocOpe_Prism::FirstShape() (cxx L242-245).
    pub fn first_shape(&self) -> Shape {
        self.my_first_shape.clone()
    }

    /// OCCT LocOpe_Prism::LastShape() (cxx L249-252).
    pub fn last_shape(&self) -> Shape {
        self.my_last_shape.clone()
    }

    /// OCCT LocOpe_Prism::Shapes(S) (cxx L256-259).
    pub fn shapes(&self, s: &Shape) -> &Vec<Shape> {
        // OCCT cxx L258: return myMap(S) — operator() asserts IsBound.
        self.my_map.get(&shape_key(s)).expect("Standard_NoSuchObject")
    }

    /// OCCT LocOpe_Prism::Curves(SCurves) (cxx L263-279).
    pub fn curves(&self, s_curves: &mut Vec<Curve3>) {
        // OCCT cxx L265.
        s_curves.clear();
        // OCCT cxx L266-267.
        let mut spt: Vec<glam::DVec3> = Vec::new();
        crate::feat::loc_ope::sample_edges(&self.my_first_shape, &mut spt);
        // OCCT cxx L268-270.
        let height = (self.my_vec.x * self.my_vec.x
            + self.my_vec.y * self.my_vec.y
            + self.my_vec.z * self.my_vec.z)
            .sqrt();
        let u1 = -2. * height;
        let u2 = 2. * height;

        // OCCT cxx L272-278.
        for pvt in &spt {
            let the_ax = Ax1::new(*pvt, self.my_vec);
            let the_lin = Curve3::Line(Line3::new(the_ax.location, the_ax.direction));
            // OCCT cxx L276: Geom_TrimmedCurve(theLin, u1, u2, true) — the
            // Sense argument has no rcad counterpart (arch. diff. #6).
            let trlin = Curve3::Trimmed(TrimmedCurve3::new(the_lin, u1, u2));
            s_curves.push(trlin);
        }
    }

    /// OCCT LocOpe_Prism::BarycCurve() (cxx L283-297).
    pub fn baryc_curve(&self) -> Curve3 {
        // OCCT cxx L285-287.
        let mut bar = glam::DVec3::ZERO;
        let mut spt: Vec<glam::DVec3> = Vec::new();
        crate::feat::loc_ope::sample_edges(&self.my_first_shape, &mut spt);
        for pvt in &spt {
            // OCCT cxx L290-291: bar.ChangeCoord() += pvt.XYZ().
            bar += *pvt;
        }
        // OCCT cxx L293.
        bar /= spt.len() as f64;
        // OCCT cxx L294-296.
        let new_ax = Ax1::new(bar, self.my_vec);
        let the_lin = Curve3::Line(Line3::new(new_ax.location, new_ax.direction));
        the_lin
    }
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
