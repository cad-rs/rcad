// OCCT BRepFeat_MakePipe.hxx L17-122 + BRepFeat_MakePipe.cxx L17-397 +
// BRepFeat_MakePipe.lxx L17-55 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePipe.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePipe.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakePipe.lxx
//
// OCCT inheritance chain (BRepFeat_MakePipe.hxx L51):
//   BRepFeat_MakePipe : BRepFeat_Form : BRepBuilderAPI_MakeShape
// Rust has no inheritance -> composition + slots trait (architecture
// decision 2026-09-08): the BRepFeat_Form base sub-object is the `form`
// field; the pure virtual slots Curves/BarycCurve are the
// BRepFeatFormSlots impl below; GlobalPerform is entered through
// brep_feat_form::global_perform(self).
//
// Architecture differences (referenced from the affected functions):
// 1. The derived class declares its own private myStatusError (hxx L118);
//    both members are carried (see brep_feat_make_prism.rs note 1).
// 2. NCollection_DataMap (mySlface) maps to HashMap keyed by (TShape ptr,
//    Location); the key shape is carried as the tuple head.
// 3. BRepAlgoAPI_Fuse/Cut are driven on the CutVehicle (with_operation).
// 4. LocOpe_Pipe is consumed through the loc_ope_pipe.rs translation (its
//    BRepFill_Pipe engine is the GAP carrier there — the result shapes stay
//    null until the TKTopAlgo/BRepFill stage lands); the Option carriers
//    reproduce the OCCT null-handle cases.

use crate::bop::algo::builder::BooleanOpType;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::brep_feat_form::{global_perform, BRepFeatForm, BRepFeatFormSlots};
use crate::feat::brep_feat_form_2::CutVehicle;
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_pipe::LocOpePipe;
use rcad_kernel::geom::Curve3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::ShapeType;
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

/// OCCT BRepFeat_MakePipe — describes functions to build pipe features
/// (BRepFeat_MakePipe.hxx L36-50).
pub struct BRepFeatMakePipe {
    /// The BRepFeat_Form base sub-object (the composition carrier).
    pub form: BRepFeatForm,
    my_pbase: Shape, // OCCT: myPbase
    // OCCT: mySlface (architecture difference #2).
    my_slface: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    my_spine: Shape, // OCCT: mySpine
    my_curves: Vec<Option<Curve3>>, // OCCT: myCurves
    my_b_curve: Option<Curve3>,     // OCCT: myBCurve
    #[allow(dead_code)]
    my_status_error: BRepFeatStatusError, // OCCT: myStatusError (the derived member)
}

impl BRepFeatMakePipe {
    /// OCCT BRepFeat_MakePipe::BRepFeat_MakePipe() (lxx: the empty
    /// constructor with myStatusError(BRepFeat_OK)).
    pub fn new() -> Self {
        BRepFeatMakePipe {
            form: BRepFeatForm::new(),
            my_pbase: Shape::null(),
            my_slface: HashMap::new(),
            my_spine: Shape::null(),
            my_curves: Vec::new(),
            my_b_curve: None,
            my_status_error: BRepFeatStatusError::OK,
        }
    }

    /// OCCT BRepFeat_MakePipe::BRepFeat_MakePipe(Sbase, Pbase, Skface,
    /// Spine, Fuse, Modify) (hxx) — Rust has no overloading: the
    /// `with_init` suffix.
    pub fn with_init(
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        spine: &Shape,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut res = BRepFeatMakePipe::new();
        res.init(sbase, pbase, skface, spine, fuse, modify);
        res
    }

    /// OCCT BRepFeat_MakePipe::Init (cxx L49-117).
    pub fn init(
        &mut self,
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        spine: &Shape,
        mode: i32,
        modify: bool,
    ) {
        // OCCT L61-64.
        self.form.my_sbase = sbase.clone();
        self.form.basis_shape_valid();
        self.form.my_skface = skface.clone();
        self.form.sketch_face_valid();
        // OCCT L65-67.
        self.my_pbase = pbase.clone();
        self.my_slface.clear();
        self.my_spine = spine.clone();
        // OCCT L68-85.
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
        // OCCT L86-87.
        self.form.my_modify = modify;
        self.form.my_just_gluer = false;
        //
        // OCCT L93-96.
        self.form.my_shape = None;
        self.form.my_map.clear();
        self.form.my_f_shape = Shape::null();
        self.form.my_l_shape = Shape::null();
        // OCCT L97-103.
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

    /// OCCT BRepFeat_MakePipe::Add (cxx L121-170).
    pub fn add(&mut self, the_e: &Shape, the_f: &Shape) {
        // OCCT L128-139.
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
        // OCCT L141-151.
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
        // OCCT L153-157.
        if !self.my_slface.contains_key(&shape_key(the_f)) {
            self.my_slface
                .insert(shape_key(the_f), (the_f.clone(), Vec::new()));
        }
        // OCCT L158-165.
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
        // OCCT L166-169.
        if !found_e {
            self.my_slface
                .get_mut(&shape_key(the_f))
                .expect("mySlface(F)")
                .1
                .push(the_e.clone());
        }
    }

    /// OCCT BRepFeat_MakePipe::Perform() (cxx L174-228).
    pub fn perform(&mut self) {
        // OCCT L181-187.
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = Shape::null();
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        self.form.perf_selection_valid();
        // OCCT L188-190.
        let the_base = self.my_pbase.clone();
        let mut the_pipe = LocOpePipe::new(&self.my_spine, &the_base);
        let vrai_pipe = the_pipe.shape().cloned().unwrap_or_else(Shape::null);
        // OCCT L191.
        maj_map(
            &self.my_pbase,
            &mut the_pipe,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        // OCCT L192-193.
        self.form.my_gshape = vrai_pipe;
        self.form.generated_shape_valid();

        // OCCT L195.
        self.form.glued_faces_valid();

        // OCCT L197-227.
        if self.form.my_glued_f.is_empty() {
            if self.form.my_fuse {
                // OCCT L201-204.
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
                // OCCT L208-211.
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
                // OCCT L215-216.
                self.form.my_shape = Some(self.form.my_gshape.clone());
                self.form.done();
            }
        } else {
            // OCCT L221-226.
            self.form.my_f_shape = the_pipe.first_shape().unwrap_or_else(Shape::null);
            let mut spt: Vec<glam::DVec3> = Vec::new();
            crate::feat::loc_ope::sample_edges(&self.form.my_f_shape, &mut spt);
            self.my_curves = the_pipe.curves(&spt).to_vec();
            self.my_b_curve = the_pipe.baryc_curve();
            global_perform(self);
        }
    }

    /// OCCT BRepFeat_MakePipe::Perform(const TopoDS_Shape& Until)
    /// (cxx L232-270).
    pub fn perform_until(&mut self, until: &Shape) {
        // OCCT L239-247.
        if until.is_null() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L248-255.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L256-269.
        let mut the_pipe = LocOpePipe::new(&self.my_spine, &self.my_pbase);
        let vrai_tuyau = the_pipe.shape().cloned().unwrap_or_else(Shape::null);
        maj_map(
            &self.my_pbase,
            &mut the_pipe,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = vrai_tuyau;
        self.form.generated_shape_valid();

        self.form.glued_faces_valid();

        self.form.my_f_shape = the_pipe.first_shape().unwrap_or_else(Shape::null);
        let mut spt: Vec<glam::DVec3> = Vec::new();
        crate::feat::loc_ope::sample_edges(&self.form.my_f_shape, &mut spt);
        self.my_curves = the_pipe.curves(&spt).to_vec();
        self.my_b_curve = the_pipe.baryc_curve();
        global_perform(self);
    }

    /// OCCT BRepFeat_MakePipe::Perform(const TopoDS_Shape& From, const
    /// TopoDS_Shape& Until) (cxx L274-332).
    pub fn perform_from_until(&mut self, from: &Shape, until: &Shape) {
        // OCCT L281-284.
        if from.is_null() || until.is_null() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L285-297.
        if !self.form.my_skface.is_null() {
            if shape_is_same(from, &self.form.my_skface) {
                self.perform_until(until);
                return;
            } else if shape_is_same(until, &self.form.my_skface) {
                self.perform_until(from);
                return;
            }
        }
        // OCCT L298-300.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionFU;
        self.form.perf_selection_valid();
        // OCCT L301-310.
        let exp = explorer(from, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L311-316.
        self.form.my_sfrom = from.clone();
        self.form.transform_shape_fu(0);
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L317-331.
        let mut the_pipe = LocOpePipe::new(&self.my_spine, &self.my_pbase);
        let vrai_tuyau = the_pipe.shape().cloned().unwrap_or_else(Shape::null);
        maj_map(
            &self.my_pbase,
            &mut the_pipe,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = vrai_tuyau;
        self.form.generated_shape_valid();

        // OCCT L324.
        self.form.glued_faces_valid();

        self.form.my_f_shape = the_pipe.first_shape().unwrap_or_else(Shape::null);
        let mut spt: Vec<glam::DVec3> = Vec::new();
        crate::feat::loc_ope::sample_edges(&self.form.my_f_shape, &mut spt);
        self.my_curves = the_pipe.curves(&spt).to_vec();
        self.my_b_curve = the_pipe.baryc_curve();
        global_perform(self);
    }
}

/// OCCT BRepFeat_Form::Curves override (cxx L339-342) — curves parallel to
/// the generating wire of the pipe.
impl BRepFeatFormSlots for BRepFeatMakePipe {
    fn curves(&mut self, s: &mut Vec<Option<Curve3>>) {
        // OCCT L341: scur = myCurves.
        *s = self.my_curves.clone();
    }

    /// OCCT BRepFeat_Form::BarycCurve override (cxx L349-352) — pass
    /// through the center of mass.
    fn baryc_curve(&mut self) -> Option<Curve3> {
        // OCCT L351: return myBCurve.
        self.my_b_curve.clone()
    }

    /// The BRepFeat_Form base sub-object.
    fn form(&mut self) -> &mut BRepFeatForm {
        &mut self.form
    }
}

/// OCCT static MajMap (cxx L356-397) — the Pipe variant.
fn maj_map(
    the_b: &Shape,
    the_p: &mut LocOpePipe,
    the_map: &mut HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_f_shape: &mut Shape,
    the_l_shape: &mut Shape,
) {
    // OCCT L364-374.
    let first_shape = the_p.first_shape().unwrap_or_else(Shape::null);
    let exp = explorer(&first_shape, ShapeType::Wire, ShapeType::Shape);
    if let Some(cur) = exp.first() {
        *the_f_shape = cur.clone();
        the_map.insert(shape_key(the_f_shape), (the_f_shape.clone(), Vec::new()));
        for exp in explorer(&first_shape, ShapeType::Face, ShapeType::Shape) {
            the_map
                .get_mut(&shape_key(the_f_shape))
                .expect("theMap(theFShape)")
                .1
                .push(exp);
        }
    }
    // OCCT L376-386.
    let last_shape = the_p.last_shape().unwrap_or_else(Shape::null);
    let exp = explorer(&last_shape, ShapeType::Wire, ShapeType::Shape);
    if let Some(cur) = exp.first() {
        *the_l_shape = cur.clone();
        the_map.insert(shape_key(the_l_shape), (the_l_shape.clone(), Vec::new()));
        for exp in explorer(&last_shape, ShapeType::Face, ShapeType::Shape) {
            the_map
                .get_mut(&shape_key(the_l_shape))
                .expect("theMap(theLShape)")
                .1
                .push(exp);
        }
    }
    // OCCT L388-396.
    for exp in explorer(the_b, ShapeType::Edge, ShapeType::Shape) {
        if !the_map.contains_key(&shape_key(&exp)) {
            let shapes = the_p.shapes(&exp).to_vec();
            the_map.insert(shape_key(&exp), (exp.clone(), shapes));
        }
    }
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
