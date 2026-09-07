// OCCT BRepFeat_Form.hxx L17-181 + BRepFeat_Form.cxx L17-1579 +
// BRepFeat_Form.lxx L17-85 + BRepFeat.cxx L202-333 / L467-520 / L524-638 /
// L642-712 (the consumed BRepFeat statics) — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Form.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Form.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Form.lxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKBool/BRepAlgo/BRepAlgo_1.cxx (IsValid)
//
// OCCT inheritance chain (BRepFeat_Form.hxx L66):
//   BRepFeat_Form : BRepBuilderAPI_MakeShape
//     BRepBuilderAPI_MakeShape : BRepBuilderAPI_Command
// Rust has no inheritance -> composition + delegation. The MakeShape base
// sub-object is carried by the BRepFeatForm fields listed below (the same
// member-by-member mapping as brep_feat_builder.rs / brep_feat_gluer.rs).
//
// Field mapping (hxx L151-176 protected + L169-176 private):
//   BRepBuilderAPI_Command::myDone          -> my_done
//   BRepBuilderAPI_MakeShape::myShape       -> my_shape (None = null)
//   BRepBuilderAPI_MakeShape::myGenerated   -> my_generated
//   BRepFeat_Form::myFuse                   -> my_fuse
//   BRepFeat_Form::myModify                 -> my_modify
//   BRepFeat_Form::myMap                    -> my_map (DataMap<Shape,
//                                              List<Shape>>; the key shape is
//                                              carried as the tuple head — the
//                                              SkipFace test of
//                                              UpdateDescendants reads the key
//                                              ShapeType)
//   BRepFeat_Form::myFShape                 -> my_f_shape
//   BRepFeat_Form::myLShape                 -> my_l_shape
//   BRepFeat_Form::myNewEdges               -> my_new_edges
//   BRepFeat_Form::myTgtEdges               -> my_tgt_edges
//   BRepFeat_Form::myPerfSelection          -> my_perf_selection
//   BRepFeat_Form::myJustGluer              -> my_just_gluer
//   BRepFeat_Form::myJustFeat               -> my_just_feat
//   BRepFeat_Form::mySbase                  -> my_sbase
//   BRepFeat_Form::mySkface                 -> my_skface
//   BRepFeat_Form::myGShape                 -> my_gshape
//   BRepFeat_Form::mySFrom                  -> my_sfrom
//   BRepFeat_Form::mySUntil                 -> my_suntil
//   BRepFeat_Form::myGluedF                 -> my_glued_f
//   BRepFeat_Form::mySbOK..myPSOK           -> my_sb_ok..my_ps_ok
//   BRepFeat_Form::myStatusError            -> my_status_error
//
// Architecture differences (referenced from the affected functions):
// 1. Curves(S) / BarycCurve() are pure virtuals of BRepFeat_Form (hxx
//    L130-132), overridden by the form-feature subclasses (MakePrism /
//    MakeDPrism / MakeRevol / MakePipe / MakeRevolutionForm — stage 3c).
//    The two virtual slots are carried as methods that stop with a pure
//    virtual marker (the C++ base-class slot is not callable either); the
//    Rust dispatch mechanism (trait object or equivalent) lands with the
//    first subclass translation.
// 2. BRepAlgoAPI_Cut is driven on the rcad (PaveFiller, Builder) vehicle,
//    re-hosted below as CutVehicle (same vehicle as
//    bop/brep_algo_api::run_build and BRepFeatBuilder::perform_bop); the
//    myFillHistory knob is set so trP.Modified/IsDeleted (=
//    BRepAlgoAPI_BuilderAlgo::Modified / IsDeleted, BuilderAlgo.cxx
//    L203-235) read the BRepTools_History carried by Builder::my_history.
// 3. BRepFeat::ParametricMinMax / FaceUntil / Tool are re-hosted below from
//    BRepFeat.cxx (the package statics have no dedicated rcad module yet —
//    same re-hosting model as BopToolsSet in brep_feat_builder.rs).
//    ParametricMinMax uses LocOpe_CSIntersector + Extrema_ExtPC (both
//    translated); the edge 3D curve is read from the rcad TShape (curve +
//    range) and the TopLoc_Location transformation of OCCT L267/L142 is the
//    identity-location reduction of the loc_ope modules.
// 4. BRepFeat::IsInside (BRepFeat.cxx L467-520) needs BRepTopAdaptor_FClass2d
//    and GCPnts_QuasiUniformDeflection (TKTopAlgo/TKMath — not translated);
//    the re-host stops at the gap with the OCCT structure documented.
// 5. BRepAlgo::IsValid(S) = BRepCheck_Analyzer(S).IsValid()
//    (BRepAlgo_1.cxx L39-43); rcad's brep_check works on a BRep pool, the
//    standalone-Shape analyzer is pending — the re-host stops at the gap.
// 6. BRepFeat_Builder::Modified(fdsc) (the Descendants static, cxx L1501) is
//    BOPAlgo_BuilderShape::Modified = myHistory->Modified(S)
//    (BOPAlgo_BuilderShape.hxx L52-56). rcad's BRepFeatBuilder does not
//    carry the myHistory member (reported gap); the re-host returns the
//    empty list — the OCCT myFillHistory=false fallback of the same accessor
//    (BuilderShape.hxx L57) — so the Descendants map stays empty and the
//    OCCT "if (!mapFuntil.IsEmpty())" guards skip exactly as in that
//    fallback.
// 7. LocOpe_Gluer::Perform (loc_ope_gluer.rs deferral note) is pending
//    (LocOpe_WiresOnShape / LocOpe_Spliter / LocOpe_Generator); the
//    theGlue.Perform() call of the still-gluer branch is carried at its spot
//    and the not-done branch of OCCT is the one taken.
// 8. NCollection_DataMap maps to HashMap keyed by (TShape ptr, Location) —
//    the TopTools_ShapeMapHasher identity; the OCCT bucket iteration order
//    is not reproduced (same reduction as brep_feat_builder.rs myImages).
// 9. BRepFeat::Tool builds the result solid through rcad's BRepBuilder pool
//    (the ResultBuilder pattern); the escaped shapes keep their TShape
//    identity, the pool locations table is not carried (identity-location
//    reduction of the loc_ope modules).

use super::brep_feat_form_2::{
    brep_algo_is_valid, brep_feat_face_until, brep_feat_is_inside, brep_feat_tool,
    brep_feat_parametric_min_max, builder_result_shape, CutVehicle,
};
use crate::feat::brep_feat_builder::{explorer, BRepFeatBuilder, OcctShapeMap};
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use crate::feat::loc_ope_find_edges::LocOpeFindEdges;
use crate::feat::loc_ope_gluer::LocOpeGluer;
use crate::feat::loc_ope_operation::LocOpeOperation;
use rcad_kernel::geom::{Curve3, Surface3};
use rcad_kernel::math::el::in_period;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use std::collections::HashMap;


/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::IsSame(S).
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT TopoDS_Shape::IsEqual(S) — IsSame + same orientation.
fn shape_is_equal(a: &Shape, b: &Shape) -> bool {
    shape_is_same(a, b) && a.orientation == b.orientation
}

/// OCCT TopAbs::Reverse (TopAbs.hxx) — FORWARD<->REVERSED, INTERNAL/
/// EXTERNAL unchanged.
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// OCCT NCollection_Map::Add — returns true when newly added.
fn map_add(m: &mut OcctShapeMap, s: &Shape) -> bool {
    m.add(shape_key(s), s.clone())
}

/// OCCT Geom_Curve::IsPeriodic over the rcad Curve3 variants (the rcad
/// CurveEval surface does not carry the periodicity query; the variant map
/// reproduces the OCCT per-class semantics: conics and helices are periodic,
/// trimmed/offset curves delegate to the basis).
fn geom_curve_is_periodic(the_c: &Curve3) -> bool {
    match the_c {
        Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::CircularHelix(_) => true,
        Curve3::Offset(o) => geom_curve_is_periodic(&o.basis),
        Curve3::Trimmed(t) => geom_curve_is_periodic(&t.curve),
        _ => false,
    }
}

/// OCCT Geom_Curve::Period over the rcad Curve3 variants (2*pi for the
/// periodic conics; only consumed behind geom_curve_is_periodic).
fn geom_curve_period(the_c: &Curve3) -> f64 {
    match the_c {
        Curve3::Circle(_) | Curve3::Ellipse(_) | Curve3::CircularHelix(_) => {
            2.0 * std::f64::consts::PI
        }
        Curve3::Offset(o) => geom_curve_period(&o.basis),
        Curve3::Trimmed(t) => geom_curve_period(&t.curve),
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Surface(fac) (architecture difference #3).
fn brep_tool_surface(fac: &Shape) -> Option<Surface3> {
    match fac.data.as_ref() {
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::NaturalRestriction(fac).
fn brep_tool_natural_restriction(fac: &Shape) -> bool {
    match fac.data.as_ref() {
        TShape::Face(fd) => fd.natural_restriction,
        _ => false,
    }
}

fn descendants(the_s: &Shape, map_f: &mut OcctShapeMap) {
    map_f.clear();
    for _fdsc in explorer(the_s, ShapeType::Face, ShapeType::Shape) {
        // OCCT L1501: theFB.Modified(fdsc) — see the gap note above.
        let a_lm: Vec<Shape> = Vec::new();
        for it in &a_lm {
            map_add(map_f, it);
        }
    }
}

/// OCCT BRepFeat_Form — provides general functions to build form features
/// (BRepFeat_Form.hxx L39-66).
pub struct BRepFeatForm {
    my_done: bool,            // BRepBuilderAPI_Command::myDone
    my_shape: Option<Shape>,  // BRepBuilderAPI_MakeShape::myShape (None = null)
    my_generated: Vec<Shape>, // BRepBuilderAPI_MakeShape::myGenerated
    // --- BRepFeat_Form protected members (hxx L151-166) ---
    pub(crate) my_fuse: bool,
    // myModify is written by the form-feature subclasses (stage 3c); the base
    // class only carries the member.
    #[allow(dead_code)]
    pub(crate) my_modify: bool,
    // OCCT myMap: DataMap<Shape, List<Shape>>; the key shape is carried as
    // the tuple head (the UpdateDescendants SkipFace test reads the key
    // ShapeType).
    my_map: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    my_f_shape: Shape,
    my_l_shape: Shape,
    my_new_edges: Vec<Shape>,
    my_tgt_edges: Vec<Shape>,
    my_perf_selection: BRepFeatPerfSelection,
    my_just_gluer: bool,
    my_just_feat: bool,
    pub(crate) my_sbase: Shape,
    // mySkface is written/read by the sketch-feature subclasses (stage 3c);
    // the base class only carries the member.
    #[allow(dead_code)]
    pub(crate) my_skface: Shape,
    pub(crate) my_gshape: Shape,
    pub(crate) my_sfrom: Shape,
    pub(crate) my_suntil: Shape,
    // OCCT myGluedF: DataMap<Shape, Shape>; the key shape is carried as the
    // tuple head (the GlobalPerform iterations read the Key and the Value
    // shapes, cxx L475/L488/L518-519).
    pub(crate) my_glued_f: HashMap<(u64, u32), (Shape, Shape)>,
    // --- BRepFeat_Form private members (hxx L169-176) ---
    my_sb_ok: bool,
    my_sk_ok: bool,
    my_gs_ok: bool,
    my_sf_ok: bool,
    my_su_ok: bool,
    my_gf_ok: bool,
    my_ps_ok: bool,
    my_status_error: BRepFeatStatusError,
}

impl BRepFeatForm {
    /// OCCT BRepFeat_Form::BRepFeat_Form() (Form.lxx L19-35).
    pub fn new() -> Self {
        BRepFeatForm {
            my_done: false,
            my_shape: None,
            my_generated: Vec::new(),
            my_fuse: false,
            my_modify: false,
            my_map: HashMap::new(),
            my_f_shape: Shape::null(),
            my_l_shape: Shape::null(),
            my_new_edges: Vec::new(),
            my_tgt_edges: Vec::new(),
            my_perf_selection: BRepFeatPerfSelection::NoSelection,
            my_just_gluer: false,
            my_just_feat: false,
            my_sbase: Shape::null(),
            my_skface: Shape::null(),
            my_gshape: Shape::null(),
            my_sfrom: Shape::null(),
            my_suntil: Shape::null(),
            my_glued_f: HashMap::new(),
            my_sb_ok: false,
            my_sk_ok: false,
            my_gs_ok: false,
            my_sf_ok: false,
            my_su_ok: false,
            my_gf_ok: false,
            my_ps_ok: false,
            my_status_error: BRepFeatStatusError::NotInitialized,
        }
    }

    // --- BRepBuilderAPI_Command / MakeShape base ---

    /// OCCT BRepBuilderAPI_Command::Done().
    fn done(&mut self) {
        self.my_done = true;
    }

    /// OCCT BRepBuilderAPI_Command::NotDone().
    fn not_done(&mut self) {
        self.my_done = false;
    }

    /// OCCT BRepBuilderAPI_Command::IsDone().
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    // --- pure virtual slots (architecture difference #1) ---

    /// OCCT BRepFeat_Form::Curves(S) = 0 (hxx L130) — pure virtual; the
    /// override lands with the form-feature subclasses (stage 3c).
    pub fn curves(&mut self, _s: &mut Vec<Option<Curve3>>) {
        panic!("pure virtual BRepFeat_Form::Curves called — overridden by the form-feature subclasses (MakePrism/MakeDPrism/MakeRevol/MakePipe/MakeRevolutionForm, stage 3c)");
    }

    /// OCCT BRepFeat_Form::BarycCurve() = 0 (hxx L132) — pure virtual; the
    /// override lands with the form-feature subclasses (stage 3c). None
    /// carries the OCCT null handle.
    pub fn baryc_curve(&mut self) -> Option<Curve3> {
        panic!("pure virtual BRepFeat_Form::BarycCurve called — overridden by the form-feature subclasses (MakePrism/MakeDPrism/MakeRevol/MakePipe/MakeRevolutionForm, stage 3c)");
    }

    // --- Form.lxx inline members (L39-85) ---

    /// OCCT BRepFeat_Form::BasisShapeValid (lxx L39-42).
    pub fn basis_shape_valid(&mut self) {
        self.my_sb_ok = true;
    }

    /// OCCT BRepFeat_Form::PerfSelectionValid (lxx L46-49).
    pub fn perf_selection_valid(&mut self) {
        self.my_ps_ok = true;
    }

    /// OCCT BRepFeat_Form::GeneratedShapeValid (lxx L53-56).
    pub fn generated_shape_valid(&mut self) {
        self.my_gs_ok = true;
    }

    /// OCCT BRepFeat_Form::ShapeFromValid (lxx L60-63).
    pub fn shape_from_valid(&mut self) {
        self.my_sf_ok = true;
    }

    /// OCCT BRepFeat_Form::ShapeUntilValid (lxx L67-70).
    pub fn shape_until_valid(&mut self) {
        self.my_su_ok = true;
    }

    /// OCCT BRepFeat_Form::GluedFacesValid (lxx L74-77).
    pub fn glued_faces_valid(&mut self) {
        self.my_gf_ok = true;
    }

    /// OCCT BRepFeat_Form::SketchFaceValid (lxx L81-84).
    pub fn sketch_face_valid(&mut self) {
        self.my_sk_ok = true;
    }

    // --- accessors ---

    /// OCCT BRepFeat_Form::CurrentStatusError (cxx L1483-1486).
    pub fn current_status_error(&self) -> BRepFeatStatusError {
        self.my_status_error
    }

    /// OCCT BRepFeat_Form::FirstShape (cxx L1339-1346). The absent-map case
    /// of the OCCT myMap(myFShape) Find (Standard_NoSuchObject) panics; the
    /// null myFShape returns myGenerated.
    pub fn first_shape(&self) -> &Vec<Shape> {
        if !self.my_f_shape.is_null() {
            return match self.my_map.get(&shape_key(&self.my_f_shape)) {
                Some((_, l)) => l,
                None => panic!("Standard_NoSuchObject"),
            };
        }
        &self.my_generated // empty list
    }

    /// OCCT BRepFeat_Form::LastShape (cxx L1350-1357).
    pub fn last_shape(&self) -> &Vec<Shape> {
        if !self.my_l_shape.is_null() {
            return match self.my_map.get(&shape_key(&self.my_l_shape)) {
                Some((_, l)) => l,
                None => panic!("Standard_NoSuchObject"),
            };
        }
        &self.my_generated // empty list
    }

    /// OCCT BRepFeat_Form::NewEdges (cxx L1361-1364).
    pub fn new_edges(&self) -> &Vec<Shape> {
        &self.my_new_edges
    }

    /// OCCT BRepFeat_Form::TgtEdges (cxx L1368-1371).
    pub fn tgt_edges(&self) -> &Vec<Shape> {
        &self.my_tgt_edges
    }

    /// OCCT BRepFeat_Form::IsDeleted(F) (cxx L1243-1250).
    pub fn is_deleted(&self, the_f: &Shape) -> bool {
        if let Some((_, l)) = self.my_map.get(&shape_key(the_f)) {
            return l.is_empty();
        }
        false
    }

    /// OCCT BRepFeat_Form::Modified(F) (cxx L1254-1281) — returns the list
    /// of generated faces (the myGenerated carrier).
    pub fn modified(&mut self, the_f: &Shape) -> &Vec<Shape> {
        // OCCT L1256: myGenerated.Clear().
        self.my_generated.clear();
        // OCCT L1257-1260.
        if !self.is_done() {
            return &self.my_generated;
        }
        // OCCT L1262-1266.
        if shape_is_equal(&self.my_sbase, the_f) {
            if let Some(my_shape) = self.my_shape.as_ref() {
                self.my_generated.push(my_shape.clone());
            }
            return &self.my_generated;
        }
        // OCCT L1268-1279.
        if let Some((_, l)) = self.my_map.get(&shape_key(the_f)) {
            for sh in l.clone() {
                if !shape_is_same(&sh, the_f) && sh.shape_type() == the_f.shape_type() {
                    self.my_generated.push(sh);
                }
            }
        }
        &self.my_generated // empty list
    }

    /// OCCT BRepFeat_Form::Generated(S) (cxx L1285-1306).
    pub fn generated(&mut self, the_s: &Shape) -> &Vec<Shape> {
        // OCCT L1287: myGenerated.Clear().
        self.my_generated.clear();
        // OCCT L1288-1291.
        if !self.is_done() {
            return &self.my_generated;
        }
        // OCCT L1292-1304: check if filter on face or not.
        if self.my_map.contains_key(&shape_key(the_s)) && the_s.shape_type() != ShapeType::Face {
            if let Some((_, l)) = self.my_map.get(&shape_key(the_s)) {
                for sh in l.clone() {
                    if !shape_is_same(&sh, the_s) {
                        self.my_generated.push(sh);
                    }
                }
            }
            return &self.my_generated;
        }
        &self.my_generated
    }

    /// OCCT BRepFeat_Form::GlobalPerform (cxx L59-1239) — topological
    /// reconstruction of the result.
    pub fn global_perform(&mut self) {
        // OCCT L68-77.
        if !self.my_sb_ok
            || !self.my_gs_ok
            || !self.my_sf_ok
            || !self.my_su_ok
            || !self.my_gf_ok
            || !self.my_sk_ok
            || !self.my_ps_ok
        {
            self.my_status_error = BRepFeatStatusError::NotInitialized;
            self.not_done();
            return;
        }
        //
        // OCCT L80-82: exp, exp2; theOpe = 2; itm.
        let mut the_ope = 2i32;
        //
        // OCCT L84-104.
        if self.my_just_feat && !self.my_fuse {
            self.my_status_error = BRepFeatStatusError::InvOption;
            self.not_done();
            return;
        } else if self.my_just_feat {
            the_ope = 2;
        } else if !self.my_glued_f.is_empty() {
            the_ope = 1;
        }
        let mut change_ope = false;
        //
        let mut from_in_shape = false;
        let mut until_in_shape = false;
        //
        // OCCT L110-133.
        if !self.my_sfrom.is_null() {
            from_in_shape = true;
            for ffrom in explorer(&self.my_sfrom, ShapeType::Face, ShapeType::Shape) {
                // OCCT L116-122: the inner explorer over mySbase.
                let mut found = false;
                for cur in explorer(&self.my_sbase, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&cur, &ffrom) {
                        found = true;
                        break;
                    }
                }
                // OCCT L123-131: if (!exp.More()).
                if !found {
                    from_in_shape = false;
                    break;
                }
            }
        }
        //
        // OCCT L135-158.
        if !self.my_suntil.is_null() {
            until_in_shape = true;
            for funtil in explorer(&self.my_suntil, ShapeType::Face, ShapeType::Shape) {
                let mut found = false;
                for cur in explorer(&self.my_sbase, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&cur, &funtil) {
                        found = true;
                        break;
                    }
                }
                if !found {
                    until_in_shape = false;
                    break;
                }
            }
        }
        //
        // OCCT L160-161: it, it2; sens = 0.
        let mut sens = 0i32;
        //
        // OCCT L163-164: NCollection_Sequence<Geom_Curve> scur; Curves(scur).
        let mut scur: Vec<Option<Curve3>> = Vec::new();
        self.curves(&mut scur);
        //
        // OCCT L166: mf, Mf, mu, Mu (carried in the loop below).
        // OCCT L168-170.
        let mut orifuntil = Orientation::Internal;
        let mut oriffrom = Orientation::Internal;
        let mut f_from = Shape::null();
        let mut f_until = Shape::null();
        //
        // OCCT L172-173: LocOpe_CSIntersector ASI1, ASI2.
        let mut a_si1 = LocOpeCSIntersector::new();
        let mut a_si2 = LocOpeCSIntersector::new();
        //
        // OCCT L175-176: NCollection_List IntList; IntList.Clear() — the dead
        // local of the OCCT source is not translated (arch. note, cf.
        // loc_ope_build_shape.rs diff #6).
        //
        // --- 1) by intersection ---
        //
        // OCCT L181-185: Intersection Tool Shape From.
        if !self.my_sfrom.is_null() {
            a_si1.init(&self.my_sfrom);
            a_si1.perform_cur(&scur);
        }
        //
        // OCCT L188-192: Intersection Tool Shape Until.
        if !self.my_suntil.is_null() {
            a_si2.init(&self.my_suntil);
            a_si2.perform_cur(&scur);
        }
        //
        // OCCT L194-327: find sens, FFrom, FUntil.
        'find_sens: for jj in 1..=(scur.len() as i32) {
            if a_si1.is_done() && a_si2.is_done() {
                // OCCT L200-203.
                if a_si1.nb_points(jj) <= 0 {
                    continue;
                }
                let mf = a_si1.point(jj, 1).parameter();
                let big_mf = a_si1.point(jj, a_si1.nb_points(jj)).parameter();
                // OCCT L206-209.
                if a_si2.nb_points(jj) <= 0 {
                    continue;
                }
                let mu = a_si2.point(jj, 1).parameter();
                let big_mu = a_si2.point(jj, a_si2.nb_points(jj)).parameter();
                // OCCT L212: if (!scur(jj)->IsPeriodic()).
                let periodic = match &scur[(jj - 1) as usize] {
                    Some(c) => geom_curve_is_periodic(c),
                    None => false, // the null handle has no OCCT counterpart
                };
                if !periodic {
                    let kf: i32;
                    let ku: i32;
                    if mu <= big_mf && mf <= big_mu {
                        // overlapping intervals
                        sens = 1;
                        kf = 1;
                        ku = a_si2.nb_points(jj);
                    } else if mu > big_mf {
                        if sens == -1 {
                            self.my_status_error = BRepFeatStatusError::IntervalOverlap;
                            self.not_done();
                            return;
                        }
                        sens = 1;
                        kf = 1;
                        ku = a_si2.nb_points(jj);
                    } else {
                        if sens == 1 {
                            self.my_status_error = BRepFeatStatusError::IntervalOverlap;
                            self.not_done();
                            return;
                        }
                        sens = -1;
                        kf = a_si1.nb_points(jj);
                        ku = 1;
                    }
                    // OCCT L245-257.
                    if oriffrom == Orientation::Internal {
                        let mut oript = a_si1.point(jj, kf).orientation();
                        if oript == Orientation::Forward || oript == Orientation::Reversed {
                            if sens == -1 {
                                oript = top_abs_reverse(oript);
                            }
                            oriffrom = top_abs_reverse(oript);
                            if let Some(f) = a_si1.point(jj, kf).face() {
                                f_from = f.clone();
                            }
                        }
                    }
                    // OCCT L258-270.
                    if orifuntil == Orientation::Internal {
                        let mut oript = a_si2.point(jj, ku).orientation();
                        if oript == Orientation::Forward || oript == Orientation::Reversed {
                            if sens == -1 {
                                oript = top_abs_reverse(oript);
                            }
                            orifuntil = oript;
                            if let Some(f) = a_si2.point(jj, ku).face() {
                                f_until = f.clone();
                            }
                        }
                    }
                }
            } else if a_si2.is_done() {
                // OCCT L275-278.
                if a_si2.nb_points(jj) <= 0 {
                    continue;
                }
                // OCCT L280-296: for base case prism on mySUntil ->
                // ambivalent direction -> preferable direction = 1.
                if sens != 1 {
                    let p1 = a_si2.point(jj, 1).parameter();
                    let p2 = a_si2.point(jj, a_si2.nb_points(jj)).parameter();
                    if p1 * p2 <= 0.0 {
                        sens = 1;
                    } else if p1 < 0.0 {
                        sens = -1;
                    } else {
                        sens = 1;
                    }
                }
                // OCCT L298-306.
                let ku = if sens == -1 {
                    1
                } else {
                    a_si2.nb_points(jj)
                };
                // OCCT L307-319.
                if orifuntil == Orientation::Internal && sens != 0 {
                    let mut oript = a_si2.point(jj, ku).orientation();
                    if oript == Orientation::Forward || oript == Orientation::Reversed {
                        if sens == -1 {
                            oript = top_abs_reverse(oript);
                        }
                        orifuntil = oript;
                        if let Some(f) = a_si2.point(jj, ku).face() {
                            f_until = f.clone();
                        }
                    }
                }
            } else {
                // OCCT L321-325.
                sens = 1;
                break 'find_sens;
            }
        }
        //
        // OCCT L329: LocOpe_Gluer theGlue.
        let mut the_glue = LocOpeGluer::new();
        //
        // --- case of gluing (OCCT L331-600) ---
        if the_ope == 1 {
            let mut collage = true;
            // OCCT L340-343: cut by FFrom && FUntil — the compound Comp.
            let mut pool = BRep::new();
            let mut b = BRepBuilder::new();
            let comp = b.make_compound(&mut pool, Vec::new());
            // OCCT L344-351.
            if !self.my_sfrom.is_null() {
                if let Some(s) = brep_feat_tool(&self.my_sfrom, &f_from, oriffrom) {
                    b.add_to_compound(&mut pool, comp.clone(), s);
                }
            }
            // OCCT L352-359.
            if !self.my_suntil.is_null() {
                if let Some(s) = brep_feat_tool(&self.my_suntil, &f_until, orifuntil) {
                    b.add_to_compound(&mut pool, comp.clone(), s);
                }
            }
            //
            // OCCT L361-364.
            let mut the_fe = LocOpeFindEdges::new();
            let mut locmap: HashMap<(u64, u32), (Shape, Vec<Shape>)> = HashMap::new();
            let comp_solids = explorer(&comp, ShapeType::Solid, ShapeType::Shape);
            // OCCT L365: if (expp.More() && !Comp.IsNull() && !myGShape.IsNull()).
            if !comp_solids.is_empty() && !comp.is_null() && !self.my_gshape.is_null() {
                // OCCT L367: BRepAlgoAPI_Cut trP(myGShape, Comp).
                let tr_p = CutVehicle::new(&self.my_gshape, &comp);
                // OCCT L368-374: exp over the SOLIDs of trP.Shape();
                // exp.Current().IsNull() == the explorer found nothing.
                let res_solids = tr_p
                    .shape()
                    .map(|s| explorer(s, ShapeType::Solid, ShapeType::Shape))
                    .unwrap_or_default();
                if res_solids.is_empty() {
                    the_ope = 2;
                    change_ope = true;
                    collage = false;
                } else {
                    // else X0 (OCCT L376-511).
                    // OCCT L377-383: only solids are preserved — theGShape.
                    let mut pool_g = BRep::new();
                    let mut b_g = BRepBuilder::new();
                    let the_gshape = b_g.make_compound(&mut pool_g, res_solids.clone());
                    // OCCT L384-389.
                    if !brep_algo_is_valid(&the_gshape) {
                        the_ope = 2;
                        change_ope = true;
                        collage = false;
                    } else {
                        // else X1 (OCCT L390-510).
                        // OCCT L392-429.
                        if !self.my_sfrom.is_null() {
                            for fac in explorer(&self.my_sfrom, ShapeType::Face, ShapeType::Shape)
                            {
                                if !from_in_shape {
                                    // OCCT L401-402.
                                    self.my_map
                                        .insert(shape_key(&fac), (fac.clone(), Vec::new()));
                                } else {
                                    // OCCT L406-407.
                                    locmap.insert(shape_key(&fac), (fac.clone(), Vec::new()));
                                }
                                // OCCT L409-427.
                                if tr_p.is_deleted(&fac) {
                                    // OCCT L409-410: empty block.
                                } else if !from_in_shape {
                                    let e = self.my_map
                                        .entry(shape_key(&fac))
                                        .or_insert_with(|| (fac.clone(), Vec::new()));
                                    e.1 = tr_p.modified(&fac);
                                    if e.1.is_empty() {
                                        e.1.push(fac.clone());
                                    }
                                } else {
                                    let e = locmap
                                        .entry(shape_key(&fac))
                                        .or_insert_with(|| (fac.clone(), Vec::new()));
                                    e.1 = tr_p.modified(&fac);
                                    if e.1.is_empty() {
                                        e.1.push(fac.clone());
                                    }
                                }
                            }
                        } // if(!mySFrom.IsNull())
                        //
                        // OCCT L431-468.
                        if !self.my_suntil.is_null() {
                            for fac in
                                explorer(&self.my_suntil, ShapeType::Face, ShapeType::Shape)
                            {
                                if !until_in_shape {
                                    self.my_map
                                        .insert(shape_key(&fac), (fac.clone(), Vec::new()));
                                } else {
                                    locmap.insert(shape_key(&fac), (fac.clone(), Vec::new()));
                                }
                                if tr_p.is_deleted(&fac) {
                                } else if !until_in_shape {
                                    let e = self.my_map
                                        .entry(shape_key(&fac))
                                        .or_insert_with(|| (fac.clone(), Vec::new()));
                                    e.1 = tr_p.modified(&fac);
                                    if e.1.is_empty() {
                                        e.1.push(fac.clone());
                                    }
                                } else {
                                    let e = locmap
                                        .entry(shape_key(&fac))
                                        .or_insert_with(|| (fac.clone(), Vec::new()));
                                    e.1 = tr_p.modified(&fac);
                                    if e.1.is_empty() {
                                        e.1.push(fac.clone());
                                    }
                                }
                            }
                        } // if(!mySUntil.IsNull())
                        //
                        // OCCT L470: UpdateDescendants(trP, theGShape, true).
                        self.update_descendants_bop(&tr_p, &the_gshape, true);
                        //
                        // OCCT L472-509.
                        the_glue.init(&self.my_sbase, &the_gshape);
                        // (Key, Value) pairs of myGluedF (the key shape is
                        // the tuple head).
                        let glued_items: Vec<(Shape, Shape)> = self
                            .my_glued_f
                            .values()
                            .map(|(k, v)| (k.clone(), v.clone()))
                            .collect();
                        for (gl, glface) in glued_items {
                            // OCCT L475-487: gl = Face(itm.Key()); ldsc.
                            let mut ldsc: Vec<Shape> = Vec::new();
                            if tr_p.is_deleted(&gl) {
                                // OCCT L477-478: empty block.
                            } else {
                                ldsc = tr_p.modified(&gl);
                                if ldsc.is_empty() {
                                    ldsc.push(gl.clone());
                                }
                            }
                            // OCCT L488-508.
                            for it in ldsc.clone() {
                                let fac = it;
                                collage = brep_feat_is_inside(&fac, &glface);
                                if !collage {
                                    the_ope = 2;
                                    change_ope = true;
                                    break;
                                } else {
                                    the_glue.bind_face(&fac, &glface);
                                    the_fe.set(&fac, &glface);
                                    the_fe.init_iterator();
                                    while the_fe.more() {
                                        let ef = the_fe.edge_from();
                                        let et = the_fe.edge_to();
                                        the_glue.bind_edge(&ef, &et);
                                        the_fe.next();
                                    }
                                }
                            }
                        }
                    } // else X1
                } // else X0
            } // if(expp.More() && !Comp.IsNull() && !myGShape.IsNull())
            else {
                // OCCT L513-547.
                the_glue.init(&self.my_sbase, &self.my_gshape);
                // (Key, Value) pairs of myGluedF: glface = Key, fac = Value
                // (cxx L518-519).
                let glued_pairs: Vec<(Shape, Shape)> = self
                    .my_glued_f
                    .values()
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect();
                for (glface, fac) in glued_pairs {
                    // OCCT L520-526: find glface among the faces of myGShape.
                    let mut found = false;
                    for cur in explorer(&self.my_gshape, ShapeType::Face, ShapeType::Shape) {
                        if shape_is_same(&cur, &glface) {
                            found = true;
                            break;
                        }
                    }
                    // OCCT L527-545: if (exp.More()).
                    if found {
                        collage = brep_feat_is_inside(&glface, &fac);
                        if !collage {
                            the_ope = 2;
                            change_ope = true;
                            break;
                        } else {
                            the_glue.bind_face(&glface, &fac);
                            the_fe.set(&glface, &fac);
                            the_fe.init_iterator();
                            while the_fe.more() {
                                let ef = the_fe.edge_from();
                                let et = the_fe.edge_to();
                                the_glue.bind_edge(&ef, &et);
                                the_fe.next();
                            }
                        }
                    }
                }
            }
            //
            // OCCT L549-569: add gluing on start and end face if necessary.
            if from_in_shape && collage {
                for fac2 in explorer(&self.my_sfrom, ShapeType::Face, ShapeType::Shape) {
                    // OCCT L557: for (it.Initialize(locmap(fac2))).
                    if let Some((_, l)) = locmap.get(&shape_key(&fac2)) {
                        for it in l.clone() {
                            let fac1 = it;
                            the_fe.set(&fac1, &fac2);
                            the_glue.bind_face(&fac1, &fac2);
                            the_fe.init_iterator();
                            while the_fe.more() {
                                let ef = the_fe.edge_from();
                                let et = the_fe.edge_to();
                                the_glue.bind_edge(&ef, &et);
                                the_fe.next();
                            }
                        }
                    }
                }
            }
            //
            // OCCT L571-591.
            if until_in_shape && collage {
                for fac2 in explorer(&self.my_suntil, ShapeType::Face, ShapeType::Shape) {
                    if let Some((_, l)) = locmap.get(&shape_key(&fac2)) {
                        for it in l.clone() {
                            let fac1 = it;
                            the_glue.bind_face(&fac1, &fac2);
                            the_fe.set(&fac1, &fac2);
                            the_fe.init_iterator();
                            while the_fe.more() {
                                let ef = the_fe.edge_from();
                                let et = the_fe.edge_to();
                                the_glue.bind_edge(&ef, &et);
                                the_fe.next();
                            }
                        }
                    }
                }
            }
            //
            // OCCT L593-599.
            let ope = the_glue.ope_type();
            if ope == LocOpeOperation::Invalid
                || (self.my_fuse && ope != LocOpeOperation::Fuse)
                || (!self.my_fuse && ope != LocOpeOperation::Cut)
                || !collage
            {
                the_ope = 2;
                change_ope = true;
            }
        }
        //
        // --- if the gluing is always applicable (OCCT L602-638) ---
        if the_ope == 1 {
            // OCCT L610: theGlue.Perform() — deferred dependency
            // (architecture difference #7): the call is carried at its spot;
            // with the deferred body myGluer stays not-done and the OCCT
            // else branch below is the one taken.
            if the_glue.is_done() {
                // OCCT L613: shshs = theGlue.ResultingShape().
                if let Some(shshs) = the_glue.resulting_shape().cloned() {
                    // OCCT L615.
                    if brep_algo_is_valid(&shshs) {
                        self.update_descendants_gluer(&the_glue);
                        // OCCT L618-619.
                        self.my_new_edges = the_glue.edges().clone();
                        self.my_tgt_edges = the_glue.tgt_edges().clone();
                        self.done();
                        self.my_shape = Some(shshs);
                    } else {
                        the_ope = 2;
                        change_ope = true;
                    }
                } else {
                    the_ope = 2;
                    change_ope = true;
                }
            } else {
                the_ope = 2;
                change_ope = true;
            }
        }
        //
        // --- case without gluing + Tool with proper dimensions
        //     (OCCT L640-652) ---
        if the_ope == 2 && change_ope && self.my_just_gluer {
            self.my_just_gluer = false;
            the_ope = 0;
        }
        //
        // --- case without gluing (OCCT L654-1236) ---
        if the_ope == 2 {
            // OCCT L662: theGShape = myGShape.
            let mut the_gshape = self.my_gshape.clone();
            // OCCT L663-669: if (ChangeOpe) — the debug trace only.
            //
            // OCCT L671-673: the compound Comp.
            let mut pool = BRep::new();
            let mut b = BRepBuilder::new();
            let comp = b.make_compound(&mut pool, Vec::new());
            // OCCT L674-706.
            if !self.my_sfrom.is_null() || !self.my_suntil.is_null() {
                if !self.my_sfrom.is_null() && !from_in_shape {
                    if let Some(s) = brep_feat_tool(&self.my_sfrom, &f_from, oriffrom) {
                        b.add_to_compound(&mut pool, comp.clone(), s);
                    }
                }
                if !self.my_suntil.is_null() && !until_in_shape {
                    if !self.my_sfrom.is_null() {
                        if !shape_is_same(&self.my_sfrom, &self.my_suntil) {
                            if let Some(s) = brep_feat_tool(&self.my_suntil, &f_until, orifuntil) {
                                b.add_to_compound(&mut pool, comp.clone(), s);
                            }
                        }
                    } else {
                        if let Some(s) = brep_feat_tool(&self.my_suntil, &f_until, orifuntil) {
                            b.add_to_compound(&mut pool, comp.clone(), s);
                        }
                    }
                }
            }
            //
            // OCCT L708-723: update type of selection.
            if self.my_perf_selection == BRepFeatPerfSelection::SelectionU && !until_in_shape {
                self.my_perf_selection = BRepFeatPerfSelection::NoSelection;
            } else if self.my_perf_selection == BRepFeatPerfSelection::SelectionFU
                && !from_in_shape
                && !until_in_shape
            {
                self.my_perf_selection = BRepFeatPerfSelection::NoSelection;
            } else if self.my_perf_selection == BRepFeatPerfSelection::SelectionShU
                && !until_in_shape
            {
                self.my_perf_selection = BRepFeatPerfSelection::NoSelection;
            }
            //
            // OCCT L725-799.
            let comp_solids = explorer(&comp, ShapeType::Solid, ShapeType::Shape);
            if !comp_solids.is_empty() && !comp.is_null() && !self.my_gshape.is_null() {
                // OCCT L728: BRepAlgoAPI_Cut trP(myGShape, Comp).
                let tr_p = CutVehicle::new(&self.my_gshape, &comp);
                // OCCT L729-736: the result is necessarily a compound.
                let res_solids = tr_p
                    .shape()
                    .map(|s| explorer(s, ShapeType::Solid, ShapeType::Shape))
                    .unwrap_or_default();
                if res_solids.is_empty() {
                    self.my_status_error = BRepFeatStatusError::EmptyCutResult;
                    self.not_done();
                    return;
                }
                // OCCT L737-743: only solids are preserved.
                the_gshape = Shape::null();
                let mut pool_g = BRep::new();
                let mut b_g = BRepBuilder::new();
                the_gshape = b_g.make_compound(&mut pool_g, res_solids.clone());
                if !brep_algo_is_valid(&the_gshape) {
                    self.my_status_error = BRepFeatStatusError::InvShape;
                    self.not_done();
                    return;
                }
                // OCCT L750-773.
                if !self.my_sfrom.is_null() {
                    if !from_in_shape {
                        for fac in explorer(&self.my_sfrom, ShapeType::Face, ShapeType::Shape) {
                            // OCCT L758-759.
                            self.my_map.insert(shape_key(&fac), (fac.clone(), Vec::new()));
                            // OCCT L760-770.
                            if tr_p.is_deleted(&fac) {
                            } else {
                                let e = self.my_map
                                        .entry(shape_key(&fac))
                                        .or_insert_with(|| (fac.clone(), Vec::new()));
                                e.1 = tr_p.modified(&fac);
                                if e.1.is_empty() {
                                    e.1.push(fac.clone());
                                }
                            }
                        }
                    }
                }
                // OCCT L774-797.
                if !self.my_suntil.is_null() {
                    if !until_in_shape {
                        for fac in explorer(&self.my_suntil, ShapeType::Face, ShapeType::Shape) {
                            self.my_map.insert(shape_key(&fac), (fac.clone(), Vec::new()));
                            if tr_p.is_deleted(&fac) {
                            } else {
                                let a_modified = tr_p.modified(&fac);
                                let e = self.my_map
                                        .entry(shape_key(&fac))
                                        .or_insert_with(|| (fac.clone(), Vec::new()));
                                e.1 = a_modified;
                                // OCCT L790: if (myMap.IsEmpty()) — the OCCT
                                // source tests the MAP (not the list) after
                                // the assignment, so the bound entry makes
                                // the test false and the append is dead;
                                // kept with the same evaluation order.
                                let b_map_empty = self.my_map.is_empty();
                                if b_map_empty {
                                    let e = self
                                        .my_map
                                        .get_mut(&shape_key(&fac))
                                        .expect("map entry");
                                    e.1.push(fac.clone());
                                }
                            }
                        }
                    }
                }
                // OCCT L798: UpdateDescendants(trP, theGShape, true).
                self.update_descendants_bop(&tr_p, &the_gshape, true);
            } // if(expp.More() && !Comp.IsNull() && !myGShape.IsNull())
            //
            // OCCT L802: generation of "just feature" for assembly = Parts of
            // tool.
            let b_flag = self.my_perf_selection != BRepFeatPerfSelection::NoSelection;
            let mut the_builder = BRepFeatBuilder::new();
            the_builder.init_with_tool(&self.my_sbase, &the_gshape);
            the_builder.set_operation_with_flag(self.my_fuse as i32, b_flag);
            the_builder.perform_bop();
            //
            // OCCT L809-810: NCollection_List lshape; PartsOfTool(lshape).
            let lshape = the_builder.parts_of_tool();
            //
            // OCCT L812-816.
            let mut pbmin = f64::MAX;
            let mut pbmax = f64::MIN;
            let mut prmin = f64::MAX - 2.0 * CONFUSION;
            let mut prmax = f64::MIN + 2.0 * CONFUSION;
            let mut flag1 = false;
            //
            // --- Selection of pieces of tool to be preserved
            //     (OCCT L818-1212) ---
            if !lshape.is_empty() && self.my_perf_selection != BRepFeatPerfSelection::NoSelection
            {
                // OCCT L823-829.
                let c = self.baryc_curve();
                let Some(c) = c else {
                    self.my_status_error = BRepFeatStatusError::EmptyBaryCurve;
                    self.not_done();
                    return;
                };
                //
                // OCCT L831-927.
                if self.my_perf_selection == BRepFeatPerfSelection::SelectionSh {
                    let (nprmin, nprmax, npbmin, npbmax, nflag1) =
                        brep_feat_parametric_min_max(&self.my_sbase, &c, false);
                    prmin = nprmin;
                    prmax = nprmax;
                    pbmin = npbmin;
                    pbmax = npbmax;
                    flag1 = nflag1;
                } else if self.my_perf_selection == BRepFeatPerfSelection::SelectionFU {
                    // OCCT L837-841 (flag1 is the in/out parameter of the two
                    // calls).
                    let (prmin1, prmax1, prbmin1, prbmax1, flag1a) =
                        brep_feat_parametric_min_max(&self.my_sfrom, &c, false);
                    flag1 = flag1a;
                    let (prmin2, prmax2, prbmin2, prbmax2, flag1b) =
                        brep_feat_parametric_min_max(&self.my_suntil, &c, false);
                    flag1 = flag1b;
                    // OCCT L843-865: case of revolutions.
                    if geom_curve_is_periodic(&c) {
                        let period = geom_curve_period(&c);
                        prmax = prmax2;
                        if flag1 {
                            prmin = in_period(prmin1, prmax - period, prmax);
                        } else {
                            prmin = prmin1.min(prmin2);
                        }
                        pbmax = prbmax2;
                        pbmin = in_period(prbmin1, pbmax - period, pbmax);
                    } else {
                        prmin = prmin1.min(prmin2);
                        prmax = prmax1.max(prmax2);
                        pbmin = prbmin1.min(prbmin2);
                        pbmax = prbmax1.max(prbmax2);
                    }
                } else if self.my_perf_selection == BRepFeatPerfSelection::SelectionShU {
                    // OCCT L869-899.
                    if !self.my_just_feat && sens == 0 {
                        sens = 1;
                    }
                    if sens == 0 {
                        self.my_status_error = BRepFeatStatusError::IncDirection;
                        self.not_done();
                        return;
                    }
                    let (prmin1, prmax1, prbmin1, prbmax1, flag1a) =
                        brep_feat_parametric_min_max(&self.my_suntil, &c, false);
                    flag1 = flag1a;
                    let (prmin2, prmax2, prbmin2, prbmax2, flag1b) =
                        brep_feat_parametric_min_max(&self.my_sbase, &c, false);
                    flag1 = flag1b;
                    if sens == 1 {
                        prmin = prmin2;
                        prmax = prmax1;
                        pbmin = prbmin2;
                        pbmax = prbmax1;
                    } else if sens == -1 {
                        prmin = prmin1;
                        prmax = prmax2;
                        pbmin = prbmin1;
                        pbmax = prbmax2;
                    }
                } else if self.my_perf_selection == BRepFeatPerfSelection::SelectionU {
                    // OCCT L901-927.
                    if sens == 0 {
                        self.my_status_error = BRepFeatStatusError::IncDirection;
                        self.not_done();
                        return;
                    }
                    // Find parts of the tool containing descendants of Shape
                    // Until (OCCT L912 passes flag1).
                    let (prmin1, prmax1, prbmin1, prbmax1, flag1a) =
                        brep_feat_parametric_min_max(&self.my_suntil, &c, false);
                    flag1 = flag1a;
                    if sens == 1 {
                        prmin = f64::MIN;
                        prmax = prmax1;
                        pbmin = f64::MIN;
                        pbmax = prbmax1;
                    } else if sens == -1 {
                        prmin = prmin1;
                        prmax = f64::MAX;
                        pbmin = prbmin1;
                        pbmax = f64::MAX;
                    }
                }
                //
                // OCCT L929-1113: finer choice of ParametricMinMax in case
                // when the tool intersects Shapes From and Until.
                let delta = CONFUSION;
                if self.my_perf_selection != BRepFeatPerfSelection::NoSelection {
                    // OCCT L939-1025.
                    if !self.my_suntil.is_null() {
                        let mut map_funtil = OcctShapeMap::new();
                        descendants(&self.my_suntil, &mut map_funtil);
                        if !map_funtil.is_empty() {
                            for it in &lshape {
                                for expf in explorer(it, ShapeType::Face, ShapeType::Shape) {
                                    if map_funtil.contains(shape_key(&expf)) {
                                        let (prmin1, prmax1, prbmin1, prbmax1, flag3) =
                                            brep_feat_parametric_min_max(&expf, &c, false);
                                        let (prmin2, prmax2, prbmin2, prbmax2, flag2) =
                                            brep_feat_parametric_min_max(it, &c, false);
                                        if sens == 1 {
                                            let mut test_ok = !flag2;
                                            if flag2 {
                                                test_ok = !flag1;
                                                if flag1 && prmax2 > prmin + delta {
                                                    test_ok = !flag3;
                                                    if flag3 && prmax1 == prmax2 {
                                                        test_ok = true;
                                                    }
                                                }
                                            }
                                            if prbmin1 < pbmax && test_ok {
                                                if flag2 {
                                                    flag1 = flag2;
                                                    prmax = prmax2;
                                                }
                                                pbmax = prbmin1;
                                            }
                                        } else if sens == -1 {
                                            let mut test_ok = !flag2;
                                            if flag2 {
                                                test_ok = !flag1;
                                                if flag1 && prmin2 < prmax - delta {
                                                    test_ok = !flag3;
                                                    if flag3 && prmin1 == prmin2 {
                                                        test_ok = true;
                                                    }
                                                }
                                            }
                                            if prbmax1 > pbmin && test_ok {
                                                if flag2 {
                                                    flag1 = flag2;
                                                    prmin = prmin2;
                                                }
                                                pbmin = prbmax1;
                                            }
                                        }
                                        let _ = (prmin1, prmax1, prbmin1, prbmax1);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    // OCCT L1026-1112.
                    if !self.my_sfrom.is_null() {
                        let mut map_ffrom = OcctShapeMap::new();
                        descendants(&self.my_sfrom, &mut map_ffrom);
                        if !map_ffrom.is_empty() {
                            for it in &lshape {
                                for expf in explorer(it, ShapeType::Face, ShapeType::Shape) {
                                    if map_ffrom.contains(shape_key(&expf)) {
                                        let (prmin1, prmax1, prbmin1, prbmax1, flag3) =
                                            brep_feat_parametric_min_max(&expf, &c, false);
                                        let (prmin2, prmax2, prbmin2, prbmax2, flag2) =
                                            brep_feat_parametric_min_max(it, &c, false);
                                        if sens == 1 {
                                            let mut test_ok = !flag2;
                                            if flag2 {
                                                test_ok = !flag1;
                                                if flag1 && prmin2 < prmax - delta {
                                                    test_ok = !flag3;
                                                    if flag3 && prmin1 == prmin2 {
                                                        test_ok = true;
                                                    }
                                                }
                                            }
                                            if prbmax1 > pbmin && test_ok {
                                                if flag2 {
                                                    flag1 = flag2;
                                                    prmin = prmin2;
                                                }
                                                pbmin = prbmax1;
                                            }
                                        } else if sens == -1 {
                                            let mut test_ok = !flag2;
                                            if flag2 {
                                                test_ok = !flag1;
                                                if flag1 && prmax2 > prmin + delta {
                                                    test_ok = !flag3;
                                                    if flag3 && prmax1 == prmax2 {
                                                        test_ok = true;
                                                    }
                                                }
                                            }
                                            if prbmin1 < pbmax && test_ok {
                                                if flag2 {
                                                    flag1 = flag2;
                                                    prmax = prmax2;
                                                }
                                                pbmax = prbmin1;
                                            }
                                        }
                                        let _ = (prmin1, prmax1, prbmin1, prbmax1);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                //
                // OCCT L1115-1176: parse PartsOfTool to preserve or not
                // depending on ParametricMinMax.
                if !self.my_just_feat {
                    let mut keep_parts = false;
                    for it in &lshape {
                        // OCCT L1124-1142.
                        let (prmin1, prmax1, prbmin1, prbmax1, flag2);
                        if geom_curve_is_periodic(&c) {
                            let period = geom_curve_period(&c);
                            let (pr, prmax1l, prb, prbmax1l, flag2l) =
                                brep_feat_parametric_min_max(it, &c, true);
                            prmax1 = prmax1l;
                            prbmax1 = prbmax1l;
                            flag2 = flag2l;
                            if flag2 {
                                prmin1 = in_period(pr, prmax1 - period, prmax1);
                            } else {
                                prmin1 = pr;
                            }
                            prbmin1 = in_period(prb, prbmax1 - period, prbmax1);
                        } else {
                            let (a, bpx, bmin, bmax, f2) =
                                brep_feat_parametric_min_max(it, &c, false);
                            prmin1 = a;
                            prmax1 = bpx;
                            prbmin1 = bmin;
                            prbmax1 = bmax;
                            flag2 = f2;
                        }
                        // OCCT L1143-1156.
                        let (pmin, pmax, min, max) = if !flag2 || !flag1 {
                            (pbmin, pbmax, prbmin1, prbmax1)
                        } else {
                            (prmin, prmax, prmin1, prmax1)
                        };
                        // OCCT L1157-1162.
                        if (min <= pmax - delta) && (max >= pmin + delta) {
                            keep_parts = true;
                            the_builder.keep_part(it);
                        }
                    }
                    // OCCT L1165-1175: case when no part of the tool is
                    // preserved.
                    if !keep_parts {
                        self.my_status_error = BRepFeatStatusError::NoParts;
                        self.not_done();
                        return;
                    }
                } else {
                    // OCCT L1177-1211: case JustFeature -> all PartsOfTool
                    // are preserved.
                    let mut pool_c = BRep::new();
                    let mut b_c = BRepBuilder::new();
                    let mut compo_children: Vec<Shape> = Vec::new();
                    for it in &lshape {
                        let (prmin1, prmax1, prbmin1, prbmax1, flag2) =
                            brep_feat_parametric_min_max(it, &c, false);
                        let (pmin, pmax, min, max) = if !flag2 || !flag1 {
                            (pbmin, pbmax, prbmin1, prbmax1)
                        } else {
                            (prmin, prmax, prmin1, prmax1)
                        };
                        if (min < pmax - delta) && (max > pmin + delta) {
                            if !it.is_null() {
                                compo_children.push(it.clone());
                            }
                        }
                    }
                    let compo = b_c.make_compound(&mut pool_c, compo_children);
                    self.my_shape = Some(compo);
                }
            }
            //
            // --- Generation of result myShape (OCCT L1214-1235) ---
            if !self.my_just_feat {
                // OCCT L1220-1228.
                if b_flag {
                    the_builder.perform_result();
                    self.my_shape = builder_result_shape(&mut the_builder);
                } else {
                    self.my_shape = builder_result_shape(&mut the_builder);
                }
                self.done();
            } else {
                // OCCT L1233: all is already done.
                self.done();
            }
        }
        //
        // OCCT L1238.
        self.my_status_error = BRepFeatStatusError::OK;
    }

    /// OCCT BRepFeat_Form::UpdateDescendants(const LocOpe_Gluer& G)
    /// (cxx L1310-1335).
    fn update_descendants_gluer(&mut self, the_g: &LocOpeGluer) {
        let keys: Vec<(u64, u32)> = self.my_map.keys().copied().collect();
        for orig in keys {
            // OCCT L1319-1320: orig = itdm.Key(); newdsc map.
            let mut newdsc = OcctShapeMap::new();
            let items = match self.my_map.get(&orig) {
                Some((_, l)) => l.clone(),
                None => Vec::new(),
            };
            // OCCT L1321-1328.
            for it in items {
                let fdsc = it;
                for it2 in the_g.descendant_faces(&fdsc) {
                    map_add(&mut newdsc, it2);
                }
            }
            // OCCT L1329: myMap.ChangeFind(orig).Clear().
            let entry = match self.my_map.get_mut(&orig) {
                Some((_, l)) => l,
                None => continue,
            };
            entry.clear();
            // OCCT L1330-1333.
            for (k, itm) in newdsc.iter() {
                let _ = k;
                entry.push(itm.clone());
            }
        }
    }

    /// OCCT BRepFeat_Form::UpdateDescendants(const BRepAlgoAPI_BooleanOperation&
    /// aBOP, const TopoDS_Shape& S, const bool SkipFace) (cxx L1512-1579).
    /// The aBOP argument carries the CutVehicle of the performed cut
    /// (architecture difference #2).
    fn update_descendants_bop(
        &mut self,
        a_bop: &CutVehicle,
        the_s: &Shape,
        skip_face: bool,
    ) {
        let keys: Vec<(u64, u32)> = self.my_map.keys().copied().collect();
        for orig_key in keys {
            // OCCT L1524: orig = itdm.Key().
            let orig = match self.my_map.get(&orig_key) {
                Some((sh, _)) => sh.clone(),
                None => continue,
            };
            // OCCT L1525-1528.
            if skip_face && orig.shape_type() == ShapeType::Face {
                continue;
            }
            let mut newdsc = OcctShapeMap::new();
            let items = match self.my_map.get(&orig_key) {
                Some((_, l)) => l.clone(),
                None => Vec::new(),
            };
            // OCCT L1531-1534.
            if items.is_empty() {
                if let Some((_, l)) = self.my_map.get_mut(&orig_key) {
                    l.push(orig.clone());
                }
            }
            // OCCT L1536-1562.
            for it in items {
                let sh = it;
                if sh.shape_type() != ShapeType::Face {
                    continue;
                }
                let fdsc = sh;
                // OCCT L1544-1551: preserved in S?
                let mut preserved = false;
                for cur in explorer(the_s, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&cur, &fdsc) {
                        preserved = true;
                        map_add(&mut newdsc, &fdsc);
                        break;
                    }
                }
                // OCCT L1552-1562: if (!exp.More()).
                if !preserved {
                    let a_lm = a_bop.modified(&fdsc);
                    for a_it in a_lm {
                        map_add(&mut newdsc, &a_it);
                    }
                }
            }
            // OCCT L1564: myMap.ChangeFind(orig).Clear().
            let entry = match self.my_map.get_mut(&orig_key) {
                Some((_, l)) => l,
                None => continue,
            };
            entry.clear();
            // OCCT L1565-1577: check the appartenance to the shape.
            for (_k, itm) in newdsc.iter() {
                for cur in explorer(the_s, ShapeType::Face, ShapeType::Shape) {
                    if shape_is_same(&cur, itm) {
                        entry.push(itm.clone());
                        break;
                    }
                }
            }
        }
    }

    /// OCCT BRepFeat_Form::TransformShapeFU(const int flag) (cxx L1378-1479)
    /// — limitation of the shape until the case of infinite faces.
    pub fn transform_shape_fu(&mut self, flag: i32) -> bool {
        let mut trf = false;
        //
        // OCCT L1385-1397.
        let shapefu = if flag == 0 {
            self.my_sfrom.clone()
        } else if flag == 1 {
            self.my_suntil.clone()
        } else {
            return trf;
        };
        //
        // OCCT L1399-1407: no faces — return an error.
        let exp = explorer(&shapefu, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            return trf;
        }
        //
        // OCCT L1409-1458: the only face. Is it infinite?
        if exp.len() == 1 {
            // OCCT L1412-1413: exp.ReInit(); fac = the single face.
            let mut fac = exp[0].clone();
            //
            // OCCT L1415-1421: S = BRep_Tool::Surface(fac); unwrap the
            // rectangular trimmed surface once.
            let Some(mut s) = brep_tool_surface(&fac) else {
                return trf;
            };
            if let Surface3::Trimmed(t) = &s {
                s = (*t.basis).clone();
            }
            // OCCT L1423-1435.
            match &s {
                Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) => {
                    let exp1 = explorer(&fac, ShapeType::Wire, ShapeType::Shape);
                    if exp1.is_empty() {
                        trf = true;
                    } else {
                        trf = brep_tool_natural_restriction(&fac);
                    }
                }
                _ => {}
            }
            if trf {
                // OCCT L1438: BRepFeat::FaceUntil(mySbase, fac).
                let my_sbase = self.my_sbase.clone();
                brep_feat_face_until(&my_sbase, &mut fac);
            }
            //
            // OCCT L1441-1457.
            if flag == 0 {
                self.my_map
                    .insert(shape_key(&self.my_sfrom), (self.my_sfrom.clone(), Vec::new()));
                if let Some((_, l)) = self.my_map.get_mut(&shape_key(&self.my_sfrom)) {
                    l.push(fac.clone());
                }
                self.my_sfrom = fac;
            } else if flag == 1 {
                self.my_map
                    .insert(shape_key(&self.my_suntil), (self.my_suntil.clone(), Vec::new()));
                if let Some((_, l)) = self.my_map.get_mut(&shape_key(&self.my_suntil)) {
                    l.push(fac.clone());
                }
                self.my_suntil = fac;
            }
        } else {
            // OCCT L1459-1468.
            for fac in exp {
                self.my_map.insert(shape_key(&fac), (fac.clone(), Vec::new()));
                if let Some((_, l)) = self.my_map.get_mut(&shape_key(&fac)) {
                    l.push(fac.clone());
                }
            }
        }
        trf
    }
}

impl Default for BRepFeatForm {
    fn default() -> Self {
        Self::new()
    }
}
