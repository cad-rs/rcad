// OCCT BRepFeat_MakeRevol.hxx L17-105 + BRepFeat_MakeRevol.cxx L17-990 +
// BRepFeat_MakeRevol.lxx L17-55 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeRevol.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeRevol.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeRevol.lxx
//
// OCCT inheritance chain (BRepFeat_MakeRevol.hxx L):
//   BRepFeat_MakeRevol : BRepFeat_Form : BRepBuilderAPI_MakeShape
// Rust has no inheritance -> composition + slots trait (architecture
// decision 2026-09-08): the BRepFeat_Form base sub-object is the `form`
// field; the pure virtual slots Curves/BarycCurve are the
// BRepFeatFormSlots impl below; GlobalPerform is entered through
// brep_feat_form::global_perform(self).
//
// Architecture differences (referenced from the affected functions):
// 1. The derived class declares its own private myStatusError (hxx L100);
//    both members are carried (see brep_feat_make_prism.rs note 1).
// 2. NCollection_DataMap (mySlface) maps to HashMap keyed by (TShape ptr,
//    Location); the key shape is carried as the tuple head.
// 3. BRepAlgoAPI_Fuse/Cut are driven on the CutVehicle
//    (with_operation / with_args_tools).
// 4. BRep_Tool::Surface is the identity-location reduction (the
//    "to apply the location" re-read of ToFuse, cxx L978-979, collapses to
//    the same identity-location surface).
// 5. The dead local `sl` of Perform(Angle) (cxx L293-294) is not
//    translated (the MakePrism Perform(Length) note).
// 6. LocOpe_Revol is consumed through the loc_ope_revol.rs translation
//    (its BRepSweep_Revol engine is the crate::brep_sweep translation);
//    BarycCurve() carries the null handle as None.

use crate::bop::algo::builder::BooleanOpType;
use crate::feat::brep_feat_builder::{explorer, BRepFeatBuilder};
use crate::feat::brep_feat_form::{global_perform, BRepFeatForm, BRepFeatFormSlots};
use crate::feat::brep_feat_form_2::{
    brep_feat_is_inside, brep_feat_parametric_barycenter, brep_feat_tool, CutVehicle,
};
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use crate::feat::loc_ope_revol::LocOpeRevol;
use rcad_kernel::geom::{Curve3, Surface3};
use rcad_kernel::math::el::in_period;
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::precision::ANGULAR;
use rcad_kernel::topo::topods::{BRep, BRepBuilder};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Shape identity key (TopTools_ShapeMapHasher: TShape + Location,
/// orientation ignored).
fn shape_key(s: &Shape) -> (u64, u32) {
    (s.ptr_id(), s.location)
}

/// OCCT TopoDS_Shape::IsSame(S).
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT TopAbs::Reverse (TopAbs.hxx).
fn top_abs_reverse(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        Orientation::Internal => Orientation::Internal,
        Orientation::External => Orientation::External,
    }
}

/// OCCT BRep_Tool::Surface(fac) (the identity-location reduction).
fn brep_tool_surface(fac: &Shape) -> Option<Surface3> {
    match fac.data.as_ref() {
        rcad_kernel::topo::topods::TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx L142-146).
fn gp_vec_is_parallel(dir1: glam::DVec3, dir2: glam::DVec3, angular_tolerance: f64) -> bool {
    let an_ang = glam::DVec3::angle_between(dir1, dir2);
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT gp_Ax3::IsCoplanar (gp_Ax3.hxx L603-620).
fn gp_ax3_is_coplanar(
    dir1: glam::DVec3,
    loc1: glam::DVec3,
    dir2: glam::DVec3,
    loc2: glam::DVec3,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> bool {
    let a_vec = loc2 - loc1;
    let mut a_d1 = dir1.dot(a_vec);
    if a_d1 < 0.0 {
        a_d1 = -a_d1;
    }
    let mut a_d2 = dir2.dot(a_vec);
    if a_d2 < 0.0 {
        a_d2 = -a_d2;
    }
    a_d1 <= linear_tolerance
        && a_d2 <= linear_tolerance
        && gp_vec_is_parallel(dir1, dir2, angular_tolerance)
}

/// OCCT BRepFeat_MakeRevol — describes functions to build revolved shell
/// features (BRepFeat_MakeRevol.hxx L36-40).
pub struct BRepFeatMakeRevol {
    /// The BRepFeat_Form base sub-object (the composition carrier).
    pub form: BRepFeatForm,
    my_pbase: Shape, // OCCT: myPbase
    // OCCT: mySlface (architecture difference #2).
    my_slface: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    my_axis: Ax1, // OCCT: myAxis (gp_Ax1)
    my_curves: Vec<Option<Curve3>>, // OCCT: myCurves
    my_b_curve: Option<Curve3>,     // OCCT: myBCurve
    #[allow(dead_code)]
    my_status_error: BRepFeatStatusError, // OCCT: myStatusError (the derived member)
}

impl BRepFeatMakeRevol {
    /// OCCT BRepFeat_MakeRevol::BRepFeat_MakeRevol() (lxx: the empty
    /// constructor with myStatusError(BRepFeat_OK)).
    pub fn new() -> Self {
        BRepFeatMakeRevol {
            form: BRepFeatForm::new(),
            my_pbase: Shape::null(),
            my_slface: HashMap::new(),
            my_axis: Ax1::new(glam::DVec3::ZERO, glam::DVec3::Z),
            my_curves: Vec::new(),
            my_b_curve: None,
            my_status_error: BRepFeatStatusError::OK,
        }
    }

    /// OCCT BRepFeat_MakeRevol::BRepFeat_MakeRevol(Sbase, Pbase, Skface,
    /// Axis, Fuse, Modify) (hxx L54-59) — Rust has no overloading: the
    /// `with_init` suffix.
    pub fn with_init(
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        axis: Ax1,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut res = BRepFeatMakeRevol::new();
        res.init(sbase, pbase, skface, axis, fuse, modify);
        res
    }

    /// OCCT BRepFeat_MakeRevol::Init (cxx L71-140).
    pub fn init(
        &mut self,
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        axis: Ax1,
        mode: i32,
        modify: bool,
    ) {
        // OCCT L83-85.
        self.my_axis = axis;
        self.my_pbase = pbase.clone();
        self.form.my_sbase = sbase.clone();
        self.form.basis_shape_valid();
        // OCCT L87-90.
        self.form.my_skface = skface.clone();
        self.form.sketch_face_valid();
        self.my_pbase = pbase.clone();
        self.my_slface.clear();
        // OCCT L91-108.
        if mode == 0 {
            self.form.my_fuse = false;
            self.form.my_just_feat = false;
        } else if mode == 1 {
            self.form.my_fuse = true;
            self.form.my_just_feat = false;
        } else if mode == 2 {
            self.form.my_fuse = true;
            self.form.my_just_feat = true;
        } else {
        }
        // OCCT L109-110.
        self.form.my_modify = modify;
        self.form.my_just_gluer = false;
        //
        // OCCT L116-119.
        self.form.my_shape = None;
        self.form.my_map.clear();
        self.form.my_f_shape = Shape::null();
        self.form.my_l_shape = Shape::null();
        // OCCT L120-126.
        for exp in explorer(&self.form.my_sbase, ShapeType::Face, ShapeType::Shape) {
            self.form
                .my_map
                .insert(shape_key(&exp), (exp.clone(), Vec::new()));
            self.form
                .my_map
                .get_mut(&shape_key(&exp))
                .expect("myMap(exp.Current())")
                .1
                .push(exp.clone());
        }
    }

    /// OCCT BRepFeat_MakeRevol::Add (cxx L147-196).
    pub fn add(&mut self, the_e: &Shape, the_f: &Shape) {
        // OCCT L154-165.
        let mut found = false;
        for exp in explorer(&self.form.my_sbase, ShapeType::Face, ShapeType::Shape) {
            if shape_is_same(&exp, the_f) {
                found = true;
                break;
            }
        }
        if !found {
            panic!("Standard_ConstructionError");
        }
        // OCCT L167-177.
        found = false;
        for exp in explorer(&self.my_pbase, ShapeType::Edge, ShapeType::Shape) {
            if shape_is_same(&exp, the_e) {
                found = true;
                break;
            }
        }
        if !found {
            panic!("Standard_ConstructionError");
        }
        // OCCT L179-183.
        if !self.my_slface.contains_key(&shape_key(the_f)) {
            self.my_slface
                .insert(shape_key(the_f), (the_f.clone(), Vec::new()));
        }
        // OCCT L184-191.
        let list = self
            .my_slface
            .get(&shape_key(the_f))
            .expect("mySlface(F)")
            .1
            .clone();
        let mut found_e = false;
        for itl in &list {
            if shape_is_same(itl, the_e) {
                found_e = true;
                break;
            }
        }
        // OCCT L192-195.
        if !found_e {
            self.my_slface
                .get_mut(&shape_key(the_f))
                .expect("mySlface(F)")
                .1
                .push(the_e.clone());
        }
    }

    /// OCCT BRepFeat_MakeRevol::Perform(const double Angle)
    /// (cxx L200-336).
    pub fn perform(&mut self, angle: f64) {
        // OCCT L207-213.
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = Shape::null();
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        self.form.perf_selection_valid();
        // OCCT L214.
        let revol_comp = 2.0 * PI - angle.abs() <= ANGULAR;
        // OCCT L215-216.
        let mut the_revol = LocOpeRevol::new();
        let angledec = 0.0f64;
        // OCCT L218-232 (the commented skface block of the OCCT source is
        // not translated).
        if revol_comp {
            self.form.my_skface = Shape::null();
        }
        // OCCT L233-240.
        if angledec == 0.0 {
            the_revol.perform(&self.my_pbase, &self.my_axis, angle);
        } else {
            the_revol.perform_angle_dec(&self.my_pbase, &self.my_axis, angle, angledec);
        }

        // OCCT L242.
        let mut vrai_revol = the_revol.shape();

        // OCCT L244.
        maj_map(
            &self.my_pbase,
            &the_revol,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // OCCT L246-257.
        self.form.my_gshape = vrai_revol.clone();
        self.form.generated_shape_valid();
        let base = the_revol.first_shape();
        let mut expf = explorer(&base, ShapeType::Face, ShapeType::Shape).into_iter();
        let the_base = expf.next().unwrap_or_else(Shape::null); // OCCT: exp.Current()
        if expf.next().is_some() {
            self.form.not_done();
            self.my_status_error = BRepFeatStatusError::InvFirstShape;
            return;
        }

        // OCCT L259-261.
        let mut f_face = Shape::null();
        let mut found = false;

        // OCCT L263-301.
        if !self.form.my_skface.is_null() || !self.my_slface.is_empty() {
            // OCCT L265-285.
            if self.form.my_l_shape.shape_type() == ShapeType::Wire {
                'ex1: for ex1 in explorer(&vrai_revol, ShapeType::Face, ShapeType::Shape) {
                    for ex2 in explorer(&ex1, ShapeType::Wire, ShapeType::Shape) {
                        if shape_is_same(&ex2, &self.form.my_l_shape) {
                            f_face = ex1.clone();
                            found = true;
                            break;
                        }
                    }
                    if found {
                        break 'ex1;
                    }
                }
            }
            // OCCT L287-300.
            for an_exp in explorer(&self.form.my_sbase, ShapeType::Face, ShapeType::Shape) {
                let ff = an_exp;
                if to_fuse(&ff, &f_face) {
                    // OCCT L293-294: the dead local sl is not translated
                    // (architecture difference #5).
                    if !shape_is_same(&f_face, &self.my_pbase)
                        && brep_feat_is_inside(&ff, &f_face)
                    {
                        break;
                    }
                }
            }
        }
        // OCCT L302.
        self.form.glued_faces_valid();
        // OCCT L303-306.
        if !self.form.my_skface.is_null() {
            verif_glued_faces(
                &self.form.my_skface,
                &the_base,
                &mut self.my_b_curve,
                &mut self.my_curves,
                &the_revol,
                &mut self.form.my_glued_f,
            );
        }

        // OCCT L308-335.
        if self.form.my_glued_f.is_empty() {
            if self.form.my_fuse {
                // OCCT L312-315.
                let f = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &self.form.my_gshape,
                    BooleanOpType::Union,
                );
                self.form.my_shape = f.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&f, &f_shape, false);
                self.form.done();
            } else if !self.form.my_fuse {
                // OCCT L319-322.
                let c = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &self.form.my_gshape,
                    BooleanOpType::Cut,
                );
                self.form.my_shape = c.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&c, &f_shape, false);
                self.form.done();
            } else {
                // OCCT L326-327.
                self.form.my_shape = Some(self.form.my_gshape.clone());
                self.form.done();
            }
        } else {
            // OCCT L332-334.
            let mut seq: Vec<Curve3> = Vec::new();
            the_revol.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = the_revol.baryc_curve();
            global_perform(self);
        }
    }

    /// OCCT BRepFeat_MakeRevol::Perform(const TopoDS_Shape& Until)
    /// (cxx L340-469).
    pub fn perform_until(&mut self, until: &Shape) {
        // OCCT L347-348.
        let mut angle = 0.0f64;
        let mut tour_complet = false;

        // OCCT L350-358.
        if until.is_null() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L359-363.
        if !self.form.my_skface.is_null() && shape_is_same(until, &self.form.my_skface) {
            angle = 2.0 * PI;
            tour_complet = true;
        }
        // OCCT L364-371.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // OCCT L375-384.
        let mut the_revol = LocOpeRevol::new();
        if !tour_complet {
            angle = 2.0 * PI - 3.0 * PI / 180.0;
        }
        the_revol.perform(&self.my_pbase, &self.my_axis, angle);
        let mut vrai_revol = the_revol.shape();
        // OCCT L386.
        maj_map(
            &self.my_pbase,
            &the_revol,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // OCCT L388-410.
        if !trf {
            self.form.my_gshape = vrai_revol.clone();
            self.form.generated_shape_valid();
            // OCCT L394-403.
            let base = the_revol.first_shape();
            let mut expf = explorer(&base, ShapeType::Face, ShapeType::Shape).into_iter();
            let _the_base = expf.next().unwrap_or_else(Shape::null); // OCCT: exp.Current()
            if expf.next().is_some() {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::InvFirstShape;
                return;
            }
            // OCCT L404.
            self.form.glued_faces_valid();
            // OCCT L405: the VerifGluedFaces call is commented out in the
            // OCCT source.
            // OCCT L407-409.
            let mut seq: Vec<Curve3> = Vec::new();
            the_revol.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = the_revol.baryc_curve();
            global_perform(self);
        } else {
            // OCCT L413-419.
            let mut seq: Vec<Curve3> = Vec::new();
            the_revol.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = the_revol.baryc_curve();
            let scur = vec![self.my_b_curve.clone()];
            let mut a_si = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            a_si.perform_cur(&scur);
            // OCCT L420-467.
            if a_si.is_done() && a_si.nb_points(1) >= 1 {
                let or_ = a_si.point(1, 1).orientation();
                let f_until = a_si
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L424-431.
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                let comp = b.make_compound(&mut pool, Vec::new());
                if let Some(s) = brep_feat_tool(&self.form.my_suntil, &f_until, or_) {
                    b.add_to_compound(&mut pool, comp.clone(), s);
                }
                // OCCT L432-433.
                let tr_p =
                    CutVehicle::with_operation(&vrai_revol, &comp, BooleanOpType::Cut);
                let cutsh = tr_p.shape().cloned().unwrap_or_else(Shape::null);
                // OCCT L434-447: the last matching solid wins (no outer
                // break in the OCCT source).
                for ex in explorer(&cutsh, ShapeType::Solid, ShapeType::Shape) {
                    for ex1 in explorer(&ex, ShapeType::Face, ShapeType::Shape) {
                        let fac = ex1;
                        if shape_is_same(&fac, &self.my_pbase) {
                            vrai_revol = ex.clone();
                            break;
                        }
                    }
                }
                // OCCT L448-466.
                if self.form.my_fuse {
                    let f = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &vrai_revol,
                        BooleanOpType::Union,
                    );
                    self.form.my_shape = f.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&f, &f_shape, false);
                    self.form.done();
                } else if !self.form.my_fuse {
                    let c = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &vrai_revol,
                        BooleanOpType::Cut,
                    );
                    self.form.my_shape = c.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&c, &f_shape, false);
                    self.form.done();
                } else {
                    self.form.my_shape = Some(vrai_revol);
                    self.form.done();
                }
            }
        }
    }

    /// OCCT BRepFeat_MakeRevol::Perform(const TopoDS_Shape& From, const
    /// TopoDS_Shape& Until) (cxx L473-648).
    pub fn perform_from_until(&mut self, from: &Shape, until: &Shape) {
        // OCCT L480-483.
        if from.is_null() || until.is_null() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L484-505.
        if !self.form.my_skface.is_null() {
            if shape_is_same(from, &self.form.my_skface) {
                self.form.my_just_gluer = true;
                self.perform_until(until);
                if self.form.my_just_gluer {
                    return;
                }
            } else if shape_is_same(until, &self.form.my_skface) {
                self.form.my_just_gluer = true;
                // OCCT L498: myAxis.Reverse().
                self.my_axis.direction = -self.my_axis.direction;
                self.perform_until(from);
                if self.form.my_just_gluer {
                    return;
                }
            }
        }
        // OCCT L507-509.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionFU;
        self.form.perf_selection_valid();
        // OCCT L511-520.
        let exp = explorer(from, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L522-534.
        self.form.my_sfrom = from.clone();
        let trff = self.form.transform_shape_fu(0);
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trfu = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        if trfu != trff {
            self.form.not_done();
            self.my_status_error = BRepFeatStatusError::IncTypes;
            return;
        }
        // OCCT L536-538.
        let mut the_revol = LocOpeRevol::new();
        the_revol.perform(&self.my_pbase, &self.my_axis, 2.0 * PI);
        let mut vrai_revol = the_revol.shape();
        // OCCT L540.
        maj_map(
            &self.my_pbase,
            &the_revol,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // OCCT L542-552.
        if !trff {
            self.form.my_gshape = vrai_revol.clone();
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            // OCCT L547: the VerifGluedFaces call is commented out in the
            // OCCT source.
            // OCCT L549-551.
            let mut seq: Vec<Curve3> = Vec::new();
            the_revol.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = the_revol.baryc_curve();
            global_perform(self);
        } else {
            // OCCT L555-563.
            let mut seq: Vec<Curve3> = Vec::new();
            the_revol.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = the_revol.baryc_curve();
            let scur = vec![self.my_b_curve.clone()];
            let mut a_si1 = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            let mut a_si2 = LocOpeCSIntersector::with_shape(&self.form.my_sfrom);
            a_si1.perform_cur(&scur);
            a_si2.perform_cur(&scur);
            // OCCT L564-566.
            let or_u: Orientation;
            let or_f: Orientation;
            let mut f_from = Shape::null();
            let mut f_until = Shape::null();
            let pr_f: f64;
            let pr_u: f64;
            // OCCT L567-578.
            if a_si1.is_done() && a_si1.nb_points(1) >= 1 {
                or_u = a_si1.point(1, 1).orientation();
                f_until = a_si1
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                pr_u = a_si1.point(1, 1).parameter();
            } else {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::NoIntersectU;
                return;
            }
            // OCCT L579-595.
            if a_si2.is_done() && a_si2.nb_points(1) >= 1 {
                let mut pr1 = a_si2.point(1, 1).parameter();
                pr1 = in_period(pr1, pr_u - 2.0 * PI, pr_u);
                let mut pr2 = a_si2.point(1, a_si2.nb_points(1)).parameter();
                pr2 = in_period(pr2, pr_u - 2.0 * PI, pr_u);
                // OCCT L585-588.
                or_f = top_abs_reverse(or_u);
                f_from = a_si2
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                pr_f = pr1.max(pr2);
            } else {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::NoIntersectF;
                return;
            }
            // OCCT L596-601.
            if !(pr_u > pr_f) {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::IncParameter;
                return;
            }
            // OCCT L602-614.
            let mut pool = BRep::new();
            let mut b = BRepBuilder::new();
            let comp = b.make_compound(&mut pool, Vec::new());
            // OCCT L605-609.
            if let Some(sf) = brep_feat_tool(&self.form.my_sfrom, &f_from, or_f) {
                b.add_to_compound(&mut pool, comp.clone(), sf);
            }
            // OCCT L610-614.
            if let Some(su) = brep_feat_tool(&self.form.my_suntil, &f_until, or_u) {
                b.add_to_compound(&mut pool, comp.clone(), su);
            }
            // OCCT L615-616.
            let tr_p =
                CutVehicle::with_operation(&vrai_revol, &comp, BooleanOpType::Cut);
            let cutsh = tr_p.shape().cloned().unwrap_or_else(Shape::null);
            // OCCT L617-627.
            vrai_revol = explorer(&cutsh, ShapeType::Solid, ShapeType::Shape)
                .into_iter()
                .next()
                .unwrap_or_else(Shape::null);
            let my_b_curve = self.my_b_curve.as_ref().expect("myBCurve").clone();
            for ex in explorer(&cutsh, ShapeType::Solid, ShapeType::Shape) {
                let pr_cur = brep_feat_parametric_barycenter(&ex, &my_b_curve);
                if pr_f <= pr_cur && pr_u >= pr_cur {
                    vrai_revol = ex.clone();
                    break;
                }
            }
            // OCCT L628-646.
            if self.form.my_fuse && !self.form.my_just_feat {
                let f = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &vrai_revol,
                    BooleanOpType::Union,
                );
                self.form.my_shape = f.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&f, &f_shape, false);
                self.form.done();
            } else if !self.form.my_fuse && !self.form.my_just_feat {
                let c = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &vrai_revol,
                    BooleanOpType::Cut,
                );
                self.form.my_shape = c.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&c, &f_shape, false);
                self.form.done();
            } else {
                self.form.my_shape = Some(vrai_revol);
                self.form.done();
            }
        }
    }

    /// OCCT BRepFeat_MakeRevol::PerformThruAll (cxx L655-663) — feature
    /// throughout the initial shape.
    pub fn perform_thru_all(&mut self) {
        // OCCT L662.
        self.perform(2.0 * PI);
    }

    /// OCCT BRepFeat_MakeRevol::PerformUntilAngle (cxx L670-792) — feature
    /// till shape Until defined with the angle.
    pub fn perform_until_angle(&mut self, until: &Shape, angle: f64) {
        // OCCT L677-684 (no return after the nested Performs — source form
        // kept).
        if until.is_null() {
            self.perform(angle);
        }
        if angle == 0.0 {
            self.perform_until(until);
        }
        // OCCT L685-689.
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L690-694.
        if !self.form.my_skface.is_null() && shape_is_same(until, &self.form.my_skface) {
            self.perform(angle);
            return;
        }
        // OCCT L695-702.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();

        // OCCT L706-708.
        let mut the_revol = LocOpeRevol::new();
        the_revol.perform(&self.my_pbase, &self.my_axis, angle);
        let mut vrai_revol = the_revol.shape();

        // OCCT L710.
        maj_map(
            &self.my_pbase,
            &the_revol,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // OCCT L712-733 (the Trf test is INVERTED w.r.t. Perform(Until) —
        // source form kept).
        if trf {
            self.form.my_gshape = vrai_revol.clone();
            self.form.generated_shape_valid();
            // OCCT L717-726.
            let base = the_revol.first_shape();
            let mut expf = explorer(&base, ShapeType::Face, ShapeType::Shape).into_iter();
            let _the_base = expf.next().unwrap_or_else(Shape::null); // OCCT: exp.Current()
            if expf.next().is_some() {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::InvFirstShape;
                return;
            }
            // OCCT L727.
            self.form.glued_faces_valid();
            // OCCT L728: the VerifGluedFaces call is commented out in the
            // OCCT source.
            // OCCT L730-732.
            let mut seq: Vec<Curve3> = Vec::new();
            the_revol.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = the_revol.baryc_curve();
            global_perform(self);
        } else {
            // OCCT L736-742.
            let mut seq: Vec<Curve3> = Vec::new();
            the_revol.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = the_revol.baryc_curve();
            let scur = vec![self.my_b_curve.clone()];
            let mut a_si = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            a_si.perform_cur(&scur);
            // OCCT L743-790.
            if a_si.is_done() && a_si.nb_points(1) >= 1 {
                let or_ = a_si.point(1, 1).orientation();
                let f_until = a_si
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L747-754.
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                let comp = b.make_compound(&mut pool, Vec::new());
                if let Some(s) = brep_feat_tool(&self.form.my_suntil, &f_until, or_) {
                    b.add_to_compound(&mut pool, comp.clone(), s);
                }
                // OCCT L755-756.
                let tr_p =
                    CutVehicle::with_operation(&vrai_revol, &comp, BooleanOpType::Cut);
                let cutsh = tr_p.shape().cloned().unwrap_or_else(Shape::null);
                // OCCT L757-770: the last matching solid wins (no outer
                // break in the OCCT source).
                for ex in explorer(&cutsh, ShapeType::Solid, ShapeType::Shape) {
                    for ex1 in explorer(&ex, ShapeType::Face, ShapeType::Shape) {
                        let fac = ex1;
                        if shape_is_same(&fac, &self.my_pbase) {
                            vrai_revol = ex.clone();
                            break;
                        }
                    }
                }
                // OCCT L771-789.
                if self.form.my_fuse {
                    let f = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &vrai_revol,
                        BooleanOpType::Union,
                    );
                    self.form.my_shape = f.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&f, &f_shape, false);
                    self.form.done();
                } else if !self.form.my_fuse {
                    let c = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &vrai_revol,
                        BooleanOpType::Cut,
                    );
                    self.form.my_shape = c.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&c, &f_shape, false);
                    self.form.done();
                } else {
                    self.form.my_shape = Some(vrai_revol);
                    self.form.done();
                }
            }
        }
    }
}

/// OCCT BRepFeat_Form::Curves override (cxx L799-802) — circles parallel to
/// the generating edge of revolution.
impl BRepFeatFormSlots for BRepFeatMakeRevol {
    fn curves(&mut self, s: &mut Vec<Option<Curve3>>) {
        // OCCT L801: scur = myCurves.
        *s = self.my_curves.clone();
    }

    /// OCCT BRepFeat_Form::BarycCurve override (cxx L809-812) — pass
    /// through the center of mass of the primitive.
    fn baryc_curve(&mut self) -> Option<Curve3> {
        // OCCT L811: return myBCurve.
        self.my_b_curve.clone()
    }

    /// The BRepFeat_Form base sub-object.
    fn form(&mut self) -> &mut BRepFeatForm {
        &mut self.form
    }
}

/// OCCT static VerifGluedFaces (cxx L820-890) — check intersection
/// Tool/theSkface = thePbase; if yes -> OK, otherwise -> case without
/// gluing.
fn verif_glued_faces(
    the_skface: &Shape,
    the_pbase: &Shape,
    the_b_curve: &mut Option<Curve3>,
    the_curves: &mut Vec<Option<Curve3>>,
    the_revol: &LocOpeRevol,
    the_map: &mut HashMap<(u64, u32), (Shape, Shape)>,
) {
    // OCCT L828-829.
    let mut glued_faces = true;
    let vrai_revol = the_revol.shape();

    // OCCT L831-835.
    let mut seq: Vec<Curve3> = Vec::new();
    the_revol.curves(&mut seq);
    *the_curves = seq.into_iter().map(Some).collect();
    *the_b_curve = the_revol.baryc_curve();
    let scur = vec![the_b_curve.clone()];
    // OCCT L836-837.
    let mut a_si = LocOpeCSIntersector::with_shape(the_skface);
    a_si.perform_cur(&scur);
    // OCCT L838-849.
    if a_si.is_done() && a_si.nb_points(1) >= 1 {
        let or_ = a_si.point(1, 1).orientation();
        let f_sk = a_si
            .point(1, 1)
            .face()
            .cloned()
            .unwrap_or_else(Shape::null);
        let mut pool = BRep::new();
        let mut b = BRepBuilder::new();
        let comp = b.make_compound(&mut pool, Vec::new());
        if let Some(s) = brep_feat_tool(the_skface, &f_sk, or_) {
            b.add_to_compound(&mut pool, comp.clone(), s);
        }
        // OCCT L850-851.
        let tr_p = CutVehicle::with_operation(&vrai_revol, &comp, BooleanOpType::Cut);
        let cutsh = tr_p.shape().cloned().unwrap_or_else(Shape::null);
        // OCCT L852-879.
        'ex: for ex in explorer(&cutsh, ShapeType::Solid, ShapeType::Shape) {
            let mut ex1_more = false;
            for ex1 in explorer(&ex, ShapeType::Face, ShapeType::Shape) {
                let fac1 = ex1;
                for ex2 in explorer(the_pbase, ShapeType::Face, ShapeType::Shape) {
                    let fac2 = ex2;
                    if shape_is_same(&fac1, &fac2) {
                        ex1_more = true;
                        break;
                    }
                }
                if ex1_more {
                    break;
                }
            }
            // OCCT L873-878.
            if ex1_more {
                continue 'ex;
            }
            glued_faces = false;
            break 'ex;
        }
        // OCCT L880-888.
        if !glued_faces {
            the_map.clear();
        }
    }
}

/// OCCT static MajMap (cxx L894-935) — the Revol variant.
fn maj_map(
    the_b: &Shape,
    the_p: &LocOpeRevol,
    the_map: &mut HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_f_shape: &mut Shape,
    the_l_shape: &mut Shape,
) {
    // OCCT L902-912.
    let exp = explorer(&the_p.first_shape(), ShapeType::Wire, ShapeType::Shape);
    if let Some(cur) = exp.first() {
        *the_f_shape = cur.clone();
        the_map.insert(shape_key(the_f_shape), (the_f_shape.clone(), Vec::new()));
        for exp in explorer(&the_p.first_shape(), ShapeType::Face, ShapeType::Shape) {
            the_map
                .get_mut(&shape_key(the_f_shape))
                .expect("theMap(theFShape)")
                .1
                .push(exp);
        }
    }
    // OCCT L914-924.
    let exp = explorer(&the_p.last_shape(), ShapeType::Wire, ShapeType::Shape);
    if let Some(cur) = exp.first() {
        *the_l_shape = cur.clone();
        the_map.insert(shape_key(the_l_shape), (the_l_shape.clone(), Vec::new()));
        for exp in explorer(&the_p.last_shape(), ShapeType::Face, ShapeType::Shape) {
            the_map
                .get_mut(&shape_key(the_l_shape))
                .expect("theMap(theLShape)")
                .1
                .push(exp);
        }
    }
    // OCCT L926-934.
    for exp in explorer(the_b, ShapeType::Edge, ShapeType::Shape) {
        if !the_map.contains_key(&shape_key(&exp)) {
            let shapes = the_p.shapes(&exp).clone();
            the_map.insert(shape_key(&exp), (exp.clone(), shapes));
        }
    }
}

/// OCCT static ToFuse (cxx L939-990) — the Revol variant (same body as the
/// MakePrism static).
fn to_fuse(the_f1: &Shape, the_f2: &Shape) -> bool {
    if the_f1.is_null() || the_f2.is_null() {
        return false;
    }
    // OCCT L949-950.
    let tollin = rcad_kernel::precision::CONFUSION;
    let tolang = ANGULAR;
    // OCCT L952-953 (the identity-location reduction; the null-surface case
    // has no OCCT counterpart).
    let Some(mut s1) = brep_tool_surface(the_f1) else {
        return false;
    };
    let Some(mut s2) = brep_tool_surface(the_f2) else {
        return false;
    };
    // OCCT L958-968.
    if let Surface3::Trimmed(t) = &s1 {
        s1 = (*t.basis).clone();
    }
    if let Surface3::Trimmed(t) = &s2 {
        s2 = (*t.basis).clone();
    }
    // OCCT L970-973.
    if std::mem::discriminant(&s1) != std::mem::discriminant(&s2) {
        return false;
    }
    // OCCT L975-987.
    let mut val_ret = false;
    if let Surface3::Plane(pl1) = &s1 {
        if let Surface3::Plane(pl2) = &s2 {
            // OCCT L983 (the location re-read of L978-979 is the
            // identity-location reduction).
            if gp_ax3_is_coplanar(
                pl1.normal,
                pl1.origin,
                pl2.normal,
                pl2.origin,
                tollin,
                tolang,
            ) {
                val_ret = true;
            }
        }
    }
    val_ret
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
