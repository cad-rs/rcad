// OCCT BRepFeat_MakePrism.hxx L17-149 + BRepFeat_MakePrism.cxx L17-1257 +
// BRepFeat_MakePrism.lxx L17-55 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePrism.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePrism.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePrism.lxx
//
// OCCT inheritance chain (BRepFeat_MakePrism.hxx L54):
//   BRepFeat_MakePrism : BRepFeat_Form : BRepBuilderAPI_MakeShape
// Rust has no inheritance -> composition + slots trait (architecture
// decision 2026-09-08): the BRepFeat_Form base sub-object is the `form`
// field; the pure virtual slots Curves/BarycCurve are the
// BRepFeatFormSlots impl below; GlobalPerform is entered through
// brep_feat_form::global_perform(self).
//
// Architecture differences (referenced from the affected functions):
// 1. The derived class declares its own private myStatusError (hxx L144)
//    on top of the base private member; both are carried (the base member
//    lives in `form`, the derived one here) — the OCCT CurrentStatusError()
//    reads the base member only.
// 2. NCollection_DataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>>
//    (mySlface) maps to HashMap keyed by (TShape ptr, Location); the key
//    shape is carried as the tuple head (arch. note of brep_feat_form.rs).
// 3. BRepAlgoAPI_Fuse/Cut are driven on the CutVehicle (the shared
//    PaveFiller+Builder vehicle of brep_feat_form_2.rs): with_operation for
//    the two-shape constructors, with_args_tools for the
//    SetArguments/SetTools/Build form.
// 4. BRep_Tool::Surface(F, loc) is the identity-location reduction of the
//    loc_ope modules; the loc1/loc2 transformations of ToFuse (cxx
//    L1243-1248) are accordingly the identity.
// 5. The dead local `sl` of Perform(Length) (cxx L274-275) is not
//    translated (same note as the IntList of BRepFeat_Form.cxx L175-176).
// 6. LocOpe_Prism is consumed through the loc_ope_prism.rs translation
//    (its BRepSweep_Prism engine is the crate::brep_sweep translation);
//    Curves() stores the sequence of handles, BarycCurve() the single
//    handle.

use crate::bop::algo::builder::BooleanOpType;
use crate::feat::brep_feat_builder::{explorer, BRepFeatBuilder};
use crate::feat::brep_feat_form::{global_perform, BRepFeatForm, BRepFeatFormSlots};
use crate::feat::brep_feat_form_2::{
    brep_feat_is_inside, brep_feat_parametric_barycenter, brep_feat_tool, CutVehicle,
};
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use crate::feat::loc_ope_prism::LocOpePrism;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, BRepBuilder};
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

/// OCCT BRep_Tool::Surface(fac) (the identity-location reduction).
fn brep_tool_surface(fac: &Shape) -> Option<Surface3> {
    match fac.data.as_ref() {
        rcad_kernel::topo::topods::TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx L142-146): the angle between the two
/// directions within theAngularTolerance, or its complement to PI.
fn gp_vec_is_parallel(dir1: DVec3, dir2: DVec3, angular_tolerance: f64) -> bool {
    let an_ang = DVec3::angle_between(dir1, dir2);
    an_ang <= angular_tolerance || std::f64::consts::PI - an_ang <= angular_tolerance
}

/// OCCT gp_Ax3::IsCoplanar (gp_Ax3.hxx L603-620) over the axis directions
/// and locations (the rcad Plane carries the normal + origin of its
/// position).
fn gp_ax3_is_coplanar(
    dir1: DVec3,
    loc1: DVec3,
    dir2: DVec3,
    loc2: DVec3,
    linear_tolerance: f64,
    angular_tolerance: f64,
) -> bool {
    // OCCT L607: gp_Vec aVec(axis.Location(), theOther.axis.Location()).
    let a_vec = loc2 - loc1;
    // OCCT L608-612.
    let mut a_d1 = dir1.dot(a_vec);
    if a_d1 < 0.0 {
        a_d1 = -a_d1;
    }
    // OCCT L613-617.
    let mut a_d2 = dir2.dot(a_vec);
    if a_d2 < 0.0 {
        a_d2 = -a_d2;
    }
    // OCCT L618-619.
    a_d1 <= linear_tolerance
        && a_d2 <= linear_tolerance
        && gp_vec_is_parallel(dir1, dir2, angular_tolerance)
}

/// OCCT BRepFeat_MakePrism — describes functions to build prism features
/// (BRepFeat_MakePrism.hxx L36-54).
pub struct BRepFeatMakePrism {
    /// The BRepFeat_Form base sub-object (the composition carrier).
    pub form: BRepFeatForm,
    my_pbase: Shape, // OCCT: myPbase
    // OCCT: mySlface (DataMap<Shape, List<Shape>>; the key shape is the
    // tuple head — architecture difference #2).
    my_slface: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    my_dir: DVec3, // OCCT: myDir (gp_Dir)
    my_curves: Vec<Option<Curve3>>, // OCCT: myCurves
    my_b_curve: Option<Curve3>,     // OCCT: myBCurve
    #[allow(dead_code)]
    my_status_error: BRepFeatStatusError, // OCCT: myStatusError (the derived member)
}

impl BRepFeatMakePrism {
    /// OCCT BRepFeat_MakePrism::BRepFeat_MakePrism() (lxx L21-24).
    pub fn new() -> Self {
        BRepFeatMakePrism {
            form: BRepFeatForm::new(),
            my_pbase: Shape::null(),
            my_slface: HashMap::new(),
            my_dir: DVec3::ZERO,
            my_curves: Vec::new(),
            my_b_curve: None,
            my_status_error: BRepFeatStatusError::OK,
        }
    }

    /// OCCT BRepFeat_MakePrism::BRepFeat_MakePrism(Sbase, Pbase, Skface,
    /// Direc, Fuse, Modify) (lxx L27-35) — Rust has no overloading: the
    /// `with_init` suffix.
    pub fn with_init(
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        direc: DVec3,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut res = BRepFeatMakePrism::new();
        res.init(sbase, pbase, skface, direc, fuse, modify);
        res
    }

    /// OCCT BRepFeat_MakePrism::Init (cxx L76-146).
    pub fn init(
        &mut self,
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        direc: DVec3,
        mode: i32,
        modify: bool,
    ) {
        // OCCT L88-89.
        self.form.my_skface = skface.clone();
        self.form.sketch_face_valid();
        // OCCT L90-91.
        self.form.my_sbase = sbase.clone();
        self.form.basis_shape_valid();
        // OCCT L92-94.
        self.my_pbase = pbase.clone();
        self.my_slface.clear();
        self.my_dir = direc;
        // OCCT L95-112.
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
        // OCCT L113-114.
        self.form.my_modify = modify;
        self.form.my_just_gluer = false;
        //
        // OCCT L120-125.
        self.form.my_shape = None;
        self.form.my_new_edges.clear();
        self.form.my_tgt_edges.clear();
        self.form.my_map.clear();
        self.form.my_f_shape = Shape::null();
        self.form.my_l_shape = Shape::null();
        // OCCT L126-132: myMap.Bind(face, [face]) for every face of mySbase.
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

    /// OCCT BRepFeat_MakePrism::Add (cxx L153-202) — add elements of sliding
    /// (edge on face).
    pub fn add(&mut self, the_e: &Shape, the_f: &Shape) {
        // OCCT L160-171: find F among the faces of mySbase.
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
        // OCCT L173-183: find E among the edges of myPbase.
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
        // OCCT L185-189.
        if !self.my_slface.contains_key(&shape_key(the_f)) {
            self.my_slface
                .insert(shape_key(the_f), (the_f.clone(), Vec::new()));
        }
        // OCCT L190-197.
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
        // OCCT L198-201.
        if !found_e {
            self.my_slface
                .get_mut(&shape_key(the_f))
                .expect("mySlface(F)")
                .1
                .push(the_e.clone());
        }
    }

    /// OCCT BRepFeat_MakePrism::Perform(const double Length) (cxx L210-319)
    /// — construction of prism of length Length and call of reconstruction
    /// topo. Rust has no overloading: Perform(Length) keeps the plain name.
    pub fn perform(&mut self, length: f64) {
        // OCCT L217-223.
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = Shape::null();
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        self.form.perf_selection_valid();
        // OCCT L224: gp_Vec V(Length * myDir).
        let v = self.my_dir * length;

        // construction of prism of height Length

        // OCCT L228-229.
        let the_prism = LocOpePrism::with_vec(&self.my_pbase, v);
        let vrai_prism = the_prism.shape();

        // management of descendants

        // OCCT L232.
        maj_map(
            &self.my_pbase,
            &the_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // OCCT L234-235.
        self.form.my_gshape = vrai_prism.clone(); // the primitive
        self.form.generated_shape_valid();

        // OCCT L237-239.
        let mut f_face = Shape::null();
        let mut found = false;

        // try to detect the faces of gluing
        // in case if the top of the prism is tangent to the initial shape

        // OCCT L244-282.
        if !self.form.my_skface.is_null() || !self.my_slface.is_empty() {
            // OCCT L246-266.
            if self.form.my_l_shape.shape_type() == ShapeType::Wire {
                'ex1: for ex1 in explorer(&vrai_prism, ShapeType::Face, ShapeType::Shape) {
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
            // OCCT L268-281.
            for exp in explorer(&self.form.my_sbase, ShapeType::Face, ShapeType::Shape) {
                let ff = exp;
                if to_fuse(&ff, &f_face) {
                    // OCCT L274-275: the dead local sl is not translated
                    // (architecture difference #5).
                    if !shape_is_same(&f_face, &self.my_pbase)
                        && brep_feat_is_inside(&ff, &f_face)
                    {
                        break;
                    }
                }
            }
        }

        // management of faces of gluing given by the user

        // OCCT L286.
        self.form.glued_faces_valid();

        // OCCT L288-294: case gluing.
        if !self.form.my_glued_f.is_empty() {
            self.form.my_just_gluer = true;
            // OCCT L291-292.
            let mut seq: Vec<Curve3> = Vec::new();
            the_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_prism.baryc_curve());
            // OCCT L293: GlobalPerform() — topological reconstruction.
            global_perform(self);
        }

        // if there is no gluing -> call of ope topo

        // OCCT L297-318.
        if !self.form.my_just_gluer {
            if self.form.my_fuse && !self.form.my_just_feat {
                // OCCT L301-304.
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
                // OCCT L308-311.
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
                // OCCT L315-316.
                self.form.my_shape = Some(self.form.my_gshape.clone());
                self.form.done();
            }
        }
    }

    /// OCCT BRepFeat_MakePrism::Perform(const TopoDS_Shape& Until)
    /// (cxx L327-432) — construction of prism oriented at the face Until,
    /// sufficiently long; call of topological reconstruction.
    pub fn perform_until(&mut self, until: &Shape) {
        // OCCT L334-342.
        if until.is_null() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L343-350.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L351-354.
        let c = test_curve(&self.my_pbase, self.my_dir);
        let sens = sens_of_prism(&c, &self.form.my_suntil);
        let height = height_max(
            &self.form.my_sbase,
            &self.form.my_skface,
            &self.form.my_sfrom,
            &self.form.my_suntil,
        );
        let v = self.my_dir * (2.0 * sens as f64 * height);

        // construction of long prism

        // OCCT L357-358.
        let the_prism = LocOpePrism::with_vec(&self.my_pbase, v);
        let vrai_prism = the_prism.shape();

        // OCCT L361-370: in case of support of face Until.
        if !trf {
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = vrai_prism.clone();
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            // OCCT L367-368.
            let mut seq: Vec<Curve3> = Vec::new();
            the_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_prism.baryc_curve());
            // OCCT L369.
            global_perform(self);
        } else {
            // until support -> passage to topological operations

            // OCCT L373-376.
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let scur = vec![Some(c.clone())];

            // direction of the prism depending on Until

            // OCCT L380-381.
            let mut a_si = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            a_si.perform_cur(&scur);
            // OCCT L382-430.
            if a_si.is_done() && a_si.nb_points(1) >= 1 {
                let mut or_: Orientation;
                if self.form.my_fuse {
                    or_ = a_si.point(1, 1).orientation();
                } else {
                    or_ = a_si.point(1, a_si.nb_points(1)).orientation();
                }
                if sens == -1 {
                    or_ = top_abs_reverse(or_);
                }
                let f_until = a_si
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L398-400: Comp compound.
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                let comp = b.make_compound(&mut pool, Vec::new());
                // OCCT L401-405.
                if let Some(s) = brep_feat_tool(&self.form.my_suntil, &f_until, or_) {
                    b.add_to_compound(&mut pool, comp.clone(), s);
                }
                // OCCT L406.
                let tr_p =
                    CutVehicle::with_operation(&vrai_prism, &comp, BooleanOpType::Cut);
                // OCCT L407.
                let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&tr_p, &trp_shape, false);
                //
                // OCCT L409-410: Cutsh = the first SOLID of trP.Shape().
                let cutsh = explorer(&trp_shape, ShapeType::Solid, ShapeType::Shape)
                    .into_iter()
                    .next()
                    .unwrap_or_else(Shape::null);
                // OCCT L411-429.
                if self.form.my_fuse && !self.form.my_just_feat {
                    // OCCT L413-416.
                    let f = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &cutsh,
                        BooleanOpType::Union,
                    );
                    self.form.my_shape = f.shape().cloned();
                    let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&f, &f_shape, false);
                    self.form.done();
                } else if !self.form.my_fuse {
                    // OCCT L420-423.
                    let c2 = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &cutsh,
                        BooleanOpType::Cut,
                    );
                    self.form.my_shape = c2.shape().cloned();
                    let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&c2, &f_shape, false);
                    self.form.done();
                } else {
                    // OCCT L427-428.
                    self.form.my_shape = Some(cutsh);
                    self.form.done();
                }
            }
        }
    }

    /// OCCT BRepFeat_MakePrism::Perform(const TopoDS_Shape& From, const
    /// TopoDS_Shape& Until) (cxx L440-652) — construction of a sufficiently
    /// long and properly oriented prism; call of topological reconstruction.
    pub fn perform_from_until(&mut self, from: &Shape, until: &Shape) {
        // OCCT L447-450.
        if from.is_null() || until.is_null() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L452-472.
        if !self.form.my_skface.is_null() {
            if shape_is_same(from, &self.form.my_skface) {
                self.form.my_just_gluer = true;
                self.perform_until(until);
                if self.form.my_just_gluer {
                    return;
                }
            } else if shape_is_same(until, &self.form.my_skface) {
                self.form.my_just_gluer = true;
                self.perform_until(from);
                if self.form.my_just_gluer {
                    return;
                }
            }
        }
        // OCCT L474-476.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionFU;
        self.form.perf_selection_valid();
        // OCCT L478-487.
        let exp = explorer(from, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L488-493.
        self.form.my_sfrom = from.clone();
        let trff = self.form.transform_shape_fu(0);
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trfu = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L494-499.
        if trfu != trff {
            self.form.not_done();
            self.my_status_error = BRepFeatStatusError::IncTypes;
            return;
        }

        // length depending on bounding boxes

        // OCCT L503-506.
        let height = height_max(
            &self.form.my_sbase,
            &self.form.my_skface,
            &self.form.my_sfrom,
            &self.form.my_suntil,
        );
        let c = test_curve(&self.my_pbase, self.my_dir);
        // OCCT L505-516: sens = direction of prism; tran = transfer of prism.
        let sens: i32;
        let tran: i32;
        if shape_is_same(from, until) {
            sens = 1;
            tran = -1;
        } else {
            sens = sens_of_prism(&c, &self.form.my_suntil);
            tran = sens * sens_of_prism(&c, &self.form.my_sfrom);
        }
        // OCCT L517-527.
        let mut the_prism = LocOpePrism::new();
        if tran < 0 {
            let v_tra = self.my_dir * (-3.0 * height * sens as f64 / 2.0);
            let v = self.my_dir * (3.0 * sens as f64 * height);
            the_prism.perform_trans(&self.my_pbase, v, v_tra);
        } else {
            let v = self.my_dir * (2.0 * sens as f64 * height);
            the_prism.perform(&self.my_pbase, v);
        }
        let vrai_prism = the_prism.shape();

        // OCCT L529-539.
        if !trff {
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = vrai_prism.clone();
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            // OCCT L536-537.
            let mut seq: Vec<Curve3> = Vec::new();
            the_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_prism.baryc_curve());
            // OCCT L538.
            global_perform(self);
        } else {
            // case until support -> topological operation

            // OCCT L542-549.
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let scur = vec![Some(c.clone())];
            let mut a_si1 = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            let mut a_si2 = LocOpeCSIntersector::with_shape(&self.form.my_sfrom);
            a_si1.perform_cur(&scur);
            a_si2.perform_cur(&scur);
            // OCCT L550-551.
            let mut or_u: Orientation;
            let mut or_f: Orientation;
            let mut f_from = Shape::null();
            let mut f_until = Shape::null();
            let par_f: f64;
            let par_u: f64;
            // OCCT L553-575.
            if a_si1.is_done() && a_si1.nb_points(1) >= 1 {
                if self.form.my_fuse {
                    or_u = a_si1.point(1, 1).orientation();
                } else {
                    or_u = a_si1.point(1, a_si1.nb_points(1)).orientation();
                }
                if sens == -1 {
                    or_u = top_abs_reverse(or_u);
                }
                f_until = a_si1
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                par_u = a_si1.point(1, 1).parameter();
            } else {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::NoIntersectU;
                return;
            }
            // OCCT L576-591.
            if a_si2.is_done() && a_si2.nb_points(1) >= 1 {
                or_f = a_si2.point(1, 1).orientation();
                if sens == 1 {
                    or_f = top_abs_reverse(or_f);
                }
                f_from = a_si2
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                par_f = a_si2.point(1, 1).parameter();
            } else {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::NoIntersectF;
                return;
            }
            // OCCT L592-598.
            if tran > 0 && par_u.abs() < par_f.abs() {
                let or_tmp = or_u;
                or_u = or_f;
                or_f = or_tmp;
            }
            //
            // OCCT L600-622.
            let mut a_l_tools: Vec<Shape> = Vec::new();
            // OCCT L601-611.
            match brep_feat_tool(&self.form.my_suntil, &f_until, or_u) {
                Some(s) => a_l_tools.push(s),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolU;
                    return;
                }
            }
            // OCCT L612-622.
            match brep_feat_tool(&self.form.my_sfrom, &f_from, or_f) {
                Some(s) => a_l_tools.push(s),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolF;
                    return;
                }
            }
            //
            // OCCT L624-630.
            let a_l_obj = vec![vrai_prism.clone()];
            let tr_p = CutVehicle::with_args_tools(a_l_obj, a_l_tools, BooleanOpType::Cut);
            // OCCT L631.
            let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
            self.form.update_descendants_bop(&tr_p, &trp_shape, false);
            // OCCT L632-651.
            if self.form.my_fuse && !self.form.my_just_feat {
                // OCCT L634-637.
                let f = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &trp_shape,
                    BooleanOpType::Union,
                );
                self.form.my_shape = f.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&f, &f_shape, false);
                self.form.done();
            } else if !self.form.my_fuse {
                // OCCT L641-644.
                let c2 = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &trp_shape,
                    BooleanOpType::Cut,
                );
                self.form.my_shape = c2.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&c2, &f_shape, false);
                self.form.done();
            } else {
                // OCCT L648-649.
                self.form.my_shape = Some(trp_shape);
                self.form.done();
            }
        }
    }

    /// OCCT BRepFeat_MakePrism::PerformUntilEnd (cxx L659-701) —
    /// construction of a prism and reconstruction.
    pub fn perform_until_end(&mut self) {
        // OCCT L666-672.
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionSh;
        self.form.perf_selection_valid();
        self.form.my_glued_f.clear();
        self.form.my_suntil = Shape::null();
        self.form.shape_until_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        // OCCT L673-674.
        let height = height_max(
            &self.form.my_sbase,
            &self.form.my_skface,
            &self.form.my_sfrom,
            &self.form.my_suntil,
        );
        let v = self.my_dir * (2.0 * height);

        // OCCT L676-677.
        let the_prism = LocOpePrism::with_vec(&self.my_pbase, v);
        let vrai_prism = the_prism.shape();

        // OCCT L679.
        maj_map(
            &self.my_pbase,
            &the_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // OCCT L681-683.
        self.form.my_gshape = vrai_prism.clone();
        self.form.generated_shape_valid();
        self.form.glued_faces_valid();

        // OCCT L685-700.
        if !self.form.my_fuse {
            // OCCT L687-694.
            let c = CutVehicle::with_operation(
                &self.form.my_sbase,
                &self.form.my_gshape,
                BooleanOpType::Cut,
            );
            if c.shape().is_some() {
                // OCCT: c.IsDone().
                self.form.my_shape = c.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&c, &f_shape, false);
                self.form.done();
            }
        } else {
            // OCCT L697-699.
            let mut seq: Vec<Curve3> = Vec::new();
            the_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_prism.baryc_curve());
            global_perform(self);
        }
    }

    /// OCCT BRepFeat_MakePrism::PerformFromEnd (cxx L705-845).
    pub fn perform_from_end(&mut self, until: &Shape) {
        // OCCT L712-715.
        if until.is_null() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L716-721.
        if !self.form.my_skface.is_null() && shape_is_same(until, &self.form.my_skface) {
            // OCCT L718: myDir.Reverse().
            self.my_dir = -self.my_dir;
            self.perform_until_end();
            return;
        }
        // OCCT L722-727.
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L728-734.
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionShU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let mut trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L735-741.
        let c = test_curve(&self.my_pbase, self.my_dir);
        let sens = sens_of_prism(&c, &self.form.my_suntil);
        let height = height_max(
            &self.form.my_sbase,
            &self.form.my_skface,
            &self.form.my_sfrom,
            &self.form.my_suntil,
        );
        let v_tra = self.my_dir * (-3.0 * height * sens as f64 / 2.0);
        let vect = self.my_dir * (3.0 * sens as f64 * height);
        let the_prism = LocOpePrism::with_vec_trans(&self.my_pbase, vect, v_tra);
        let vrai_prism = the_prism.shape();

        // OCCT L743-753: case face until.
        if !trf {
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = vrai_prism.clone();
            self.form.generated_shape_valid();
            self.form.my_glued_f.clear();
            self.form.glued_faces_valid();
            // OCCT L750-751.
            let mut seq: Vec<Curve3> = Vec::new();
            the_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_prism.baryc_curve());
            // OCCT L752.
            global_perform(self);
        } else {
            // case support

            // OCCT L756-763.
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let scur = vec![Some(c.clone())];
            let mut a_si1 = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            let mut a_si2 = LocOpeCSIntersector::with_shape(&self.form.my_sbase);
            a_si1.perform_cur(&scur);
            a_si2.perform_cur(&scur);
            // OCCT L764-765.
            let mut or_u = Orientation::Forward;
            let mut or_f = Orientation::Forward;
            let mut f_until = Shape::null();
            let mut f_from = Shape::null();
            // OCCT L766-774.
            if a_si1.is_done() && a_si1.nb_points(1) >= 1 {
                or_u = a_si1.point(1, 1).orientation();
                if sens == -1 {
                    or_u = top_abs_reverse(or_u);
                }
                f_until = a_si1
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
            }
            // OCCT L775-789.
            if a_si2.is_done() && a_si2.nb_points(1) >= 1 {
                or_f = a_si2.point(1, 1).orientation();
                // OCCT L778: if(sens==1) OrF = TopAbs::Reverse(OrF);
                // (commented out in the OCCT source).
                f_from = a_si2
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L780-784: the surface of FFrom, basis of a trimmed
                // surface (the identity-location reduction).
                let Some(mut s) = brep_tool_surface(&f_from) else {
                    // OCCT: the surface of an intersection face is never
                    // null; no OCCT counterpart.
                    panic!("Geom_Surface null");
                };
                if let Surface3::Trimmed(t) = &s {
                    s = (*t.basis).clone();
                }
                // OCCT L785-786: BRepLib_MakeFace fac(S, Precision::Confusion()).
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                self.form.my_sfrom = b.make_face(&mut pool, Some(s), Shape::null());
                // OCCT L787-788.
                trf = self.form.transform_shape_fu(0);
                f_from = self.form.my_sfrom.clone();
            }

            // OCCT L791-814.
            let mut a_l_tools: Vec<Shape> = Vec::new();
            // OCCT L792-802.
            match brep_feat_tool(&self.form.my_suntil, &f_until, or_u) {
                Some(sol) => a_l_tools.push(sol),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolU;
                    return;
                }
            }
            // OCCT L804-814.
            match brep_feat_tool(&self.form.my_sfrom, &f_from, or_f) {
                Some(sol1) => a_l_tools.push(sol1),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolF;
                    return;
                }
            }
            //
            // OCCT L816-822.
            let a_l_obj = vec![vrai_prism.clone()];
            let tr_p = CutVehicle::with_args_tools(a_l_obj, a_l_tools, BooleanOpType::Cut);
            //
            // OCCT L824.
            let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
            self.form.update_descendants_bop(&tr_p, &trp_shape, false);
            // OCCT L825-843.
            if self.form.my_fuse && !self.form.my_just_feat {
                // OCCT L827-830.
                let f = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &trp_shape,
                    BooleanOpType::Union,
                );
                self.form.my_shape = f.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&f, &f_shape, false);
                self.form.done();
            } else if !self.form.my_fuse {
                // OCCT L834-837.
                let c2 = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &trp_shape,
                    BooleanOpType::Cut,
                );
                self.form.my_shape = c2.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&c2, &f_shape, false);
                self.form.done();
            } else {
                // OCCT L841-842.
                self.form.my_shape = Some(trp_shape);
                self.form.done();
            }
        }
    }

    /// OCCT BRepFeat_MakePrism::PerformThruAll (cxx L849-898) — builds an
    /// infinite prism; the infinite descendants will not be kept in the
    /// result.
    pub fn perform_thru_all(&mut self) {
        // OCCT L856-870.
        self.form.my_suntil = Shape::null();
        self.form.shape_until_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        if !self.form.my_fuse {
            self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        } else {
            self.form.my_perf_selection = BRepFeatPerfSelection::SelectionSh;
        }
        self.form.perf_selection_valid();
        self.form.my_glued_f.clear();
        self.form.glued_faces_valid();

        // OCCT L872-876.
        let height = height_max(
            &self.form.my_sbase,
            &self.form.my_skface,
            &self.form.my_sfrom,
            &self.form.my_suntil,
        );
        let v = self.my_dir * (3.0 * height);
        let v_tra = self.my_dir * (-3.0 * height / 2.0);
        let the_prism = LocOpePrism::with_vec_trans(&self.my_pbase, v, v_tra);
        let vrai_prism = the_prism.shape();
        // OCCT L877.
        maj_map(
            &self.my_pbase,
            &the_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );

        // OCCT L879-880.
        self.form.my_gshape = vrai_prism.clone();
        self.form.generated_shape_valid();

        // OCCT L882-897.
        if !self.form.my_fuse {
            // OCCT L884-890.
            let c = CutVehicle::with_operation(
                &self.form.my_sbase,
                &self.form.my_gshape,
                BooleanOpType::Cut,
            );
            if c.shape().is_some() {
                // OCCT: c.IsDone().
                self.form.my_shape = c.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&c, &f_shape, false);
                self.form.done();
            }
        } else {
            // OCCT L894-896.
            let mut seq: Vec<Curve3> = Vec::new();
            the_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_prism.baryc_curve());
            global_perform(self);
        }
    }

    /// OCCT BRepFeat_MakePrism::PerformUntilHeight (cxx L902-1003).
    pub fn perform_until_height(&mut self, until: &Shape, length: f64) {
        // OCCT L909-916 (no return after the nested Performs — source form
        // kept).
        if until.is_null() {
            self.perform(length);
        }
        if length == 0.0 {
            self.perform_until(until);
        }
        // OCCT L917-921.
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L922-929.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L930-934.
        let c = test_curve(&self.my_pbase, self.my_dir);
        let sens = sens_of_prism(&c, &self.form.my_suntil);
        let v = self.my_dir * (sens as f64 * length);
        let the_prism = LocOpePrism::with_vec(&self.my_pbase, v);
        let vrai_prism = the_prism.shape();

        // OCCT L936-946.
        if !trf {
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = vrai_prism.clone();
            self.form.generated_shape_valid();
            self.form.glued_faces_valid();
            // OCCT L943-944.
            let mut seq: Vec<Curve3> = Vec::new();
            the_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_prism.baryc_curve());
            // OCCT L945.
            global_perform(self);
        } else {
            // OCCT L949-954.
            maj_map(
                &self.my_pbase,
                &the_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let scur = vec![Some(c.clone())];
            let mut a_si = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            a_si.perform_cur(&scur);
            // OCCT L955-1001.
            if a_si.is_done() && a_si.nb_points(1) >= 1 {
                let mut or_: Orientation;
                if self.form.my_fuse {
                    or_ = a_si.point(1, 1).orientation();
                } else {
                    or_ = a_si.point(1, a_si.nb_points(1)).orientation();
                }
                if sens == -1 {
                    or_ = top_abs_reverse(or_);
                }
                let f_until = a_si
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L971-978: Comp compound.
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                let comp = b.make_compound(&mut pool, Vec::new());
                // OCCT L974-978.
                if let Some(s) = brep_feat_tool(&self.form.my_suntil, &f_until, or_) {
                    b.add_to_compound(&mut pool, comp.clone(), s);
                }
                // OCCT L980.
                let tr_p =
                    CutVehicle::with_operation(&vrai_prism, &comp, BooleanOpType::Cut);
                // OCCT L981.
                let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&tr_p, &trp_shape, false);
                // OCCT L982-1000 (the result is trP.Shape(), not a solid
                // extraction — source form kept).
                if self.form.my_fuse && !self.form.my_just_feat {
                    // OCCT L984-987.
                    let f = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &trp_shape,
                        BooleanOpType::Union,
                    );
                    self.form.my_shape = f.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&f, &f_shape, false);
                    self.form.done();
                } else if !self.form.my_fuse {
                    // OCCT L991-994.
                    let c2 = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &trp_shape,
                        BooleanOpType::Cut,
                    );
                    self.form.my_shape = c2.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&c2, &f_shape, false);
                    self.form.done();
                } else {
                    // OCCT L998-999.
                    self.form.my_shape = Some(trp_shape);
                    self.form.done();
                }
            }
        }
    }
}

/// OCCT BRepFeat_Form::Curves override (cxx L1010-1013) — sequence of curves
/// parallel to the axis of prism.
impl BRepFeatFormSlots for BRepFeatMakePrism {
    fn curves(&mut self, s: &mut Vec<Option<Curve3>>) {
        // OCCT L1012: scur = myCurves.
        *s = self.my_curves.clone();
    }

    /// OCCT BRepFeat_Form::BarycCurve override (cxx L1021-1024) — curve
    /// parallel to the axis of the prism passing through the center of
    /// masses.
    fn baryc_curve(&mut self) -> Option<Curve3> {
        // OCCT L1023: return myBCurve.
        self.my_b_curve.clone()
    }

    /// The BRepFeat_Form base sub-object.
    fn form(&mut self) -> &mut BRepFeatForm {
        &mut self.form
    }
}

/// OCCT static HeightMax (cxx L1032-1102) — calculate the height of the
/// prism following the parameters of bounding box.
fn height_max(
    the_sbase: &Shape,
    the_skface: &Shape,
    the_sfrom: &Shape,
    the_suntil: &Shape,
) -> f64 {
    // OCCT L1037-1039: Bnd_Box Box; BRepBndLib::Add(theSbase); Add(theSkface).
    // The rcad shape_box evaluates the per-shape box; the accumulation is
    // the union (the identity-location reduction).
    let mut box_min = glam::DVec3::splat(f64::MAX);
    let mut box_max = glam::DVec3::splat(f64::MIN);
    add_to_box(the_sbase, &mut box_min, &mut box_max);
    add_to_box(the_skface, &mut box_min, &mut box_max);
    // OCCT L1040-1059.
    if !the_sfrom.is_null() {
        let mut fac_revol_infini = false;
        for exp in explorer(the_sfrom, ShapeType::Edge, ShapeType::Shape) {
            if explorer(&exp, ShapeType::Vertex, ShapeType::Shape).is_empty() {
                fac_revol_infini = true;
                break;
            }
        }
        if !fac_revol_infini {
            add_to_box(the_sfrom, &mut box_min, &mut box_max);
        }
    }
    // OCCT L1060-1079.
    if !the_suntil.is_null() {
        let mut fac_revol_infini = false;
        for exp in explorer(the_suntil, ShapeType::Edge, ShapeType::Shape) {
            if explorer(&exp, ShapeType::Vertex, ShapeType::Shape).is_empty() {
                fac_revol_infini = true;
                break;
            }
        }
        if !fac_revol_infini {
            add_to_box(the_suntil, &mut box_min, &mut box_max);
        }
    }
    // OCCT L1081-1095: Box.Get(c[0], c[2], c[4], c[1], c[3], c[5]) —
    // c = (xmin, xmax, ymin, ymax, zmin, zmax).
    let c = [
        box_min.x,
        box_max.x,
        box_min.y,
        box_max.y,
        box_min.z,
        box_max.z,
    ];
    let mut parmin = c[0];
    let mut parmax = c[0];
    for i in 0..6 {
        if c[i] > parmax {
            parmax = c[i];
        }
        if c[i] < parmin {
            parmin = c[i];
        }
    }
    // OCCT L1097.
    let height = (2.0 * (parmax - parmin)).abs();
    height
}

/// The BRepBndLib::Add carrier of height_max (the rcad shape_box union).
fn add_to_box(the_s: &Shape, box_min: &mut DVec3, box_max: &mut DVec3) {
    if let Some((mn, mx)) = BRepFeatBuilder::shape_box(the_s, &[]) {
        *box_min = box_min.min(mn);
        *box_max = box_max.max(mx);
    }
}

/// OCCT static SensOfPrism (cxx L1108-1131) — direction of the prism
/// depending on the shape Until.
fn sens_of_prism(the_c: &Curve3, until: &Shape) -> i32 {
    // OCCT L1110-1113.
    let mut a_si1 = LocOpeCSIntersector::with_shape(until);
    let scur = vec![Some(the_c.clone())];
    a_si1.perform_cur(&scur);
    let mut sens = 1;
    // OCCT L1115-1122.
    if a_si1.is_done() && a_si1.nb_points(1) >= 1 {
        if a_si1.point(1, 1).parameter() + CONFUSION < 0.0
            && a_si1.point(1, a_si1.nb_points(1)).parameter() + CONFUSION < 0.0
        {
            sens = -1;
        }
    } else if brep_feat_parametric_barycenter(until, the_c) < 0.0 {
        // OCCT L1123-1126.
        sens = -1;
    } else {
        // OCCT L1127-1129.
    }
    sens
}

/// OCCT static MajMap (cxx L1135-1176) — the management of myMap, myFShape
/// and myLShape.
fn maj_map(
    the_b: &Shape,
    the_p: &LocOpePrism,
    the_map: &mut HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_f_shape: &mut Shape,
    the_l_shape: &mut Shape,
) {
    // OCCT L1143-1153.
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
    // OCCT L1155-1165.
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
    // OCCT L1167-1175.
    for exp in explorer(the_b, ShapeType::Edge, ShapeType::Shape) {
        if !the_map.contains_key(&shape_key(&exp)) {
            // OCCT L1172-1173: theMap.Bind(cur, thelist); theMap(cur) =
            // theP.Shapes(cur).
            let shapes = the_p.shapes(&exp).clone();
            the_map.insert(shape_key(&exp), (exp.clone(), shapes));
        }
    }
}

/// OCCT static TestCurve (cxx L1180-1194) — the line through the
/// barycentre of the sampled edges, along V.
fn test_curve(the_base: &Shape, v: DVec3) -> Curve3 {
    // OCCT L1182-1184.
    let mut bar = DVec3::ZERO;
    let mut spt: Vec<DVec3> = Vec::new();
    crate::feat::loc_ope::sample_edges(the_base, &mut spt);
    // OCCT L1185-1189.
    for pvt in &spt {
        bar += *pvt;
    }
    // OCCT L1190.
    bar /= spt.len() as f64;
    // OCCT L1191-1193: gp_Ax1 newAx(bar, V); Geom_Line.
    Curve3::Line(Line3::new(bar, v))
}

/// OCCT static ToFuse (cxx L1198-1257) — the two planar faces are coplanar
/// (same surface type, plane positions coplanar within Precision).
fn to_fuse(the_f1: &Shape, the_f2: &Shape) -> bool {
    if the_f1.is_null() || the_f2.is_null() {
        return false;
    }
    // OCCT L1206-1209: tollin/tolang.
    let tollin = CONFUSION;
    let tolang = rcad_kernel::precision::ANGULAR;
    // OCCT L1211-1212: S1 = BRep_Tool::Surface(F1, loc1); S2 = ... — the
    // identity-location reduction (architecture difference #4); the OCCT
    // null-surface case has no counterpart (the handles are never null for
    // the faces consumed here), so the rcad None returns false.
    let Some(mut s1) = brep_tool_surface(the_f1) else {
        return false;
    };
    let Some(mut s2) = brep_tool_surface(the_f2) else {
        return false;
    };
    // OCCT L1217-1227: unwrap the rectangular trimmed surfaces.
    if let Surface3::Trimmed(t) = &s1 {
        s1 = (*t.basis).clone();
    }
    if let Surface3::Trimmed(t) = &s2 {
        s2 = (*t.basis).clone();
    }
    // OCCT L1229-1232: typS1 != typS2.
    if std::mem::discriminant(&s1) != std::mem::discriminant(&s2) {
        return false;
    }
    // OCCT L1234-1254.
    let mut val_ret = false;
    if let Surface3::Plane(pl1) = &s1 {
        if let Surface3::Plane(pl2) = &s2 {
            // OCCT L1250: pl1.Position().IsCoplanar(pl2.Position(), tollin,
            // tolang) — the location transformations of OCCT L1241-1248 are
            // the identity (architecture difference #4).
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
