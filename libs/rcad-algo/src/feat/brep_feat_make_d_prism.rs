// OCCT BRepFeat_MakeDPrism.hxx L17-156 + BRepFeat_MakeDPrism.cxx L17-1313 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeDPrism.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeDPrism.cxx
//
// OCCT inheritance chain (BRepFeat_MakeDPrism.hxx L52):
//   BRepFeat_MakeDPrism : BRepFeat_Form : BRepBuilderAPI_MakeShape
// Rust has no inheritance -> composition + slots trait (architecture
// decision 2026-09-08): the BRepFeat_Form base sub-object is the `form`
// field; the pure virtual slots Curves/BarycCurve are the
// BRepFeatFormSlots impl below; GlobalPerform is entered through
// brep_feat_form::global_perform(self).
//
// Architecture differences (referenced from the affected functions):
// 1. The derived class declares its own private myStatusError (hxx L153);
//    both members are carried (see brep_feat_make_prism.rs note 1).
// 2. NCollection_DataMap (mySlface) maps to HashMap keyed by (TShape ptr,
//    Location); the key shape is carried as the tuple head.
// 3. BRepAlgoAPI_Fuse/Cut are driven on the CutVehicle
//    (with_operation / with_args_tools).
// 4. BRep_Tool::Surface is the identity-location reduction; the trimmed
//    surface unwrap of TestCurve / PerformFromEnd (cxx L758-762, L1297-1301)
//    is the Surface3::Trimmed carrier.
// 5. NCollection_Map (the MapE of BossEdges, cxx L1060) maps to a
//    HashMap<(TShape ptr, Location), Shape>; the OCCT bucket iteration
//    order is not reproduced (same reduction as brep_feat_form.rs myMap).
// 6. Geom_Curve::Reversed is the per-variant carrier geom_curve_reversed
//    below (the DPrism test curve is a Geom_Line; Geom_Line::Reversed
//    reverses the direction).
// 7. The myNewEdges fixing pass tail of Perform(Until) (cxx L403-413) calls
//    BRepAlgo::IsValid (the brep_algo_is_valid GAP marker), then
//    BRep_Builder::SameRange/SameParameter (the rcad kernel edge flags —
//    the Arc-shared TShape is not mutated through the handle) and
//    BRepLib::SameParameter (the topalgo BRepLib stub).

use crate::bop::algo::builder::BooleanOpType;
use crate::feat::brep_feat_builder::{explorer, BRepFeatBuilder};
use crate::feat::brep_feat_form::{global_perform, BRepFeatForm, BRepFeatFormSlots};
use crate::feat::brep_feat_form_2::{
    brep_algo_is_valid, brep_feat_parametric_barycenter, brep_feat_tool, CutVehicle,
};
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use crate::feat::loc_ope_d_prism::LocOpeDPrism;
use crate::topalgo::brep_lib::brep_lib::BRepLib;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3, Surface3};
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
        TShape::Face(fd) => fd.surface.clone(),
        _ => None,
    }
}

/// OCCT BRep_Tool::Degenerated(edg).
fn brep_tool_degenerated(edg: &Shape) -> bool {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::Tolerance(edg).
fn brep_tool_tolerance(edg: &Shape) -> f64 {
    match edg.data.as_ref() {
        TShape::Edge(ed) => ed.tolerance,
        _ => 0.0,
    }
}

/// OCCT TopExp::Vertices(E, V1, V2) — the two vertices of the edge.
fn top_exp_vertices(the_e: &Shape) -> (Shape, Shape) {
    match the_e.data.as_ref() {
        TShape::Edge(ed) => (ed.first.clone(), ed.last.clone()),
        _ => (Shape::null(), Shape::null()),
    }
}

/// OCCT Geom_Curve::Reversed over the rcad Curve3 variants — the reverse
/// parameterization carrier (architecture difference #6: the DPrism test
/// curve is a Geom_Line and Geom_Line::Reversed reverses the direction).
fn geom_curve_reversed(the_c: &Curve3) -> Curve3 {
    match the_c {
        Curve3::Line(l) => Curve3::Line(Line3::new(l.origin, -l.direction)),
        _ => panic!("unreached: the DPrism test curve is a Geom_Line"),
    }
}

/// OCCT BRepFeat_MakeDPrism — describes functions to build draft prism
/// features (BRepFeat_MakeDPrism.hxx L36-52).
pub struct BRepFeatMakeDPrism {
    /// The BRepFeat_Form base sub-object (the composition carrier).
    pub form: BRepFeatForm,
    my_pbase: Shape, // OCCT: myPbase
    // OCCT: mySlface (architecture difference #2).
    my_slface: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    my_angle: f64, // OCCT: myAngle
    my_curves: Vec<Option<Curve3>>, // OCCT: myCurves
    my_b_curve: Option<Curve3>,     // OCCT: myBCurve
    my_top_edges: Vec<Shape>,       // OCCT: myTopEdges
    my_lat_edges: Vec<Shape>,       // OCCT: myLatEdges
    #[allow(dead_code)]
    my_status_error: BRepFeatStatusError, // OCCT: myStatusError (the derived member)
}

impl BRepFeatMakeDPrism {
    /// OCCT BRepFeat_MakeDPrism::BRepFeat_MakeDPrism() (hxx L76-80).
    pub fn new() -> Self {
        BRepFeatMakeDPrism {
            form: BRepFeatForm::new(),
            my_pbase: Shape::null(),
            my_slface: HashMap::new(),
            my_angle: f64::MAX, // OCCT: RealLast()
            my_curves: Vec::new(),
            my_b_curve: None,
            my_top_edges: Vec::new(),
            my_lat_edges: Vec::new(),
            my_status_error: BRepFeatStatusError::OK,
        }
    }

    /// OCCT BRepFeat_MakeDPrism::BRepFeat_MakeDPrism(Sbase, Pbase, Skface,
    /// Angle, Fuse, Modify) (hxx L66-74) — Rust has no overloading: the
    /// `with_init` suffix.
    pub fn with_init(
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        angle: f64,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut res = BRepFeatMakeDPrism::new();
        res.init(sbase, pbase, skface, angle, fuse, modify);
        res
    }

    /// OCCT BRepFeat_MakeDPrism::Init (cxx L76-148).
    pub fn init(
        &mut self,
        sbase: &Shape,
        pbase: &Shape,
        skface: &Shape,
        angle: f64,
        mode: i32,
        modify: bool,
    ) {
        // OCCT L89-90.
        self.form.my_skface = skface.clone();
        self.form.sketch_face_valid();
        // OCCT L91-92.
        self.form.my_sbase = sbase.clone();
        self.form.basis_shape_valid();
        // OCCT L93-94.
        self.my_pbase = pbase.clone();
        self.my_slface.clear();
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
        self.form.my_map.clear();
        self.form.my_f_shape = Shape::null();
        self.form.my_l_shape = Shape::null();
        self.my_top_edges.clear();
        self.my_lat_edges.clear();
        // OCCT L126-132.
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
        // OCCT L133.
        self.my_angle = angle;
    }

    /// OCCT BRepFeat_MakeDPrism::Add (cxx L152-201).
    pub fn add(&mut self, the_e: &Shape, the_f: &Shape) {
        // OCCT L159-170.
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
        // OCCT L172-182.
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
        // OCCT L184-188.
        if !self.my_slface.contains_key(&shape_key(the_f)) {
            self.my_slface
                .insert(shape_key(the_f), (the_f.clone(), Vec::new()));
        }
        // OCCT L189-196.
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
        // OCCT L197-200.
        if !found_e {
            self.my_slface
                .get_mut(&shape_key(the_f))
                .expect("mySlface(F)")
                .1
                .push(the_e.clone());
        }
    }

    /// OCCT BRepFeat_MakeDPrism::Perform(const double Height)
    /// (cxx L205-276).
    pub fn perform(&mut self, height: f64) {
        // OCCT L212-218.
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = Shape::null();
        self.form.shape_until_valid();
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        self.form.perf_selection_valid();
        // OCCT L220.
        let theheight = height / self.my_angle.cos();
        // OCCT L223-224.
        let the_d_prism = LocOpeDPrism::with_height_angle(&self.my_pbase, theheight, self.my_angle);
        let vrai_d_prism = the_d_prism.shape();
        // OCCT L226.
        maj_map(
            &self.my_pbase,
            &the_d_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        // OCCT L228-239.
        self.form.my_gshape = vrai_d_prism.clone();
        self.form.generated_shape_valid();
        // OCCT L230-232: Base = theDPrism.FirstShape(); theBase = the face.
        let base = the_d_prism.first_shape();
        let mut expf = explorer(&base, ShapeType::Face, ShapeType::Shape).into_iter();
        let _the_base = expf.next().unwrap_or_else(Shape::null); // OCCT: exp.Current()
        // OCCT L233-239.
        if expf.next().is_some() {
            self.form.not_done();
            self.my_status_error = BRepFeatStatusError::InvFirstShape;
            return;
        }

        // management of gluing faces

        // OCCT L243.
        self.form.glued_faces_valid();

        // OCCT L245-251: case gluing.
        if !self.form.my_glued_f.is_empty() {
            self.form.my_just_gluer = true;
            // OCCT L248-249.
            let mut seq: Vec<Curve3> = Vec::new();
            the_d_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_d_prism.baryc_curve());
            // OCCT L250.
            global_perform(self);
        }

        // if there is no gluing -> call topological operations

        // OCCT L254-275.
        if !self.form.my_just_gluer {
            if self.form.my_fuse {
                // OCCT L258-261.
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
                // OCCT L265-268.
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
                // OCCT L272-273.
                self.form.my_shape = Some(self.form.my_gshape.clone());
                self.form.done();
            }
        }
    }

    /// OCCT BRepFeat_MakeDPrism::Perform(const TopoDS_Shape& Until)
    /// (cxx L283-414) — feature limited by the shape Until.
    pub fn perform_until(&mut self, until: &Shape) {
        // OCCT L290-298.
        if until.is_null() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L301-308.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L309-313.
        let c = test_curve(&self.my_pbase);
        let sens = sens_of_prism(&c, &self.form.my_suntil);
        let height = sens as f64
            * height_max(
                &self.form.my_sbase,
                &self.form.my_skface,
                &self.form.my_sfrom,
                &self.form.my_suntil,
            );
        // OCCT L314-315.
        let the_d_prism = LocOpeDPrism::with_height_angle(&self.my_pbase, height, self.my_angle);
        let vrai_d_prism = the_d_prism.shape();
        // OCCT L316-336.
        if !trf {
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = vrai_d_prism.clone();
            self.form.generated_shape_valid();
            // OCCT L321-330.
            let base = the_d_prism.first_shape();
            let mut expf = explorer(&base, ShapeType::Face, ShapeType::Shape).into_iter();
            let _the_base = expf.next().unwrap_or_else(Shape::null); // OCCT: exp.Current()
            if expf.next().is_some() {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::InvFirstShape;
                return;
            }
            // OCCT L332-335.
            self.form.glued_faces_valid();
            let mut seq: Vec<Curve3> = Vec::new();
            the_d_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_d_prism.baryc_curve());
            global_perform(self);
        } else {
            // OCCT L339-354.
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let c1 = if sens == -1 {
                Some(geom_curve_reversed(c.as_ref().expect("test curve")))
            } else {
                c.clone()
            };
            let scur = vec![c1];
            let mut a_si = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            a_si.perform_cur(&scur);
            // OCCT L355-401.
            if a_si.is_done() && a_si.nb_points(1) >= 1 {
                let mut or_: Orientation;
                if self.form.my_fuse {
                    or_ = a_si.point(1, 1).orientation();
                } else {
                    or_ = a_si.point(1, a_si.nb_points(1)).orientation();
                }
                let f_until = a_si
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L369-375: Comp compound.
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                let comp = b.make_compound(&mut pool, Vec::new());
                // OCCT L371-375.
                if let Some(s) = brep_feat_tool(&self.form.my_suntil, &f_until, or_) {
                    pool_add(&mut pool, &mut b, &comp, &s);
                }
                // OCCT L377.
                let tr_p =
                    CutVehicle::with_operation(&vrai_d_prism, &comp, BooleanOpType::Cut);
                // OCCT L378.
                let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&tr_p, &trp_shape, false);
                // OCCT L380-381.
                let cutsh = explorer(&trp_shape, ShapeType::Solid, ShapeType::Shape)
                    .into_iter()
                    .next()
                    .unwrap_or_else(Shape::null);
                // OCCT L382-400.
                if self.form.my_fuse {
                    let f = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &cutsh,
                        BooleanOpType::Union,
                    );
                    self.form.my_shape = f.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&f, &f_shape, false);
                    self.form.done();
                } else if !self.form.my_fuse {
                    let c2 = CutVehicle::with_operation(
                        &self.form.my_sbase,
                        &cutsh,
                        BooleanOpType::Cut,
                    );
                    self.form.my_shape = c2.shape().cloned();
                    let f_shape =
                        self.form.my_shape.clone().unwrap_or_else(Shape::null);
                    self.form.update_descendants_bop(&c2, &f_shape, false);
                    self.form.done();
                } else {
                    self.form.my_shape = Some(cutsh);
                    self.form.done();
                }
            }
        }
        // OCCT L403-413: the myNewEdges fixing pass (architecture
        // difference #7).
        for ited in self.form.my_new_edges.clone() {
            let ledg = ited;
            if !brep_algo_is_valid(&ledg) {
                // OCCT L409-411: bB.SameRange(ledg, false);
                // bB.SameParameter(ledg, false); BRepLib::SameParameter(ledg,
                // BRep_Tool::Tolerance(ledg)).
                let tolerance = brep_tool_tolerance(&ledg);
                BRepLib::same_parameter(&ledg, tolerance);
            }
        }
    }

    /// OCCT BRepFeat_MakeDPrism::Perform(const TopoDS_Shape& From, const
    /// TopoDS_Shape& Until) (cxx L418-616).
    pub fn perform_from_until(&mut self, from: &Shape, until: &Shape) {
        // OCCT L425-428.
        if from.is_null() || until.is_null() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L430-450.
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
        // OCCT L453-455.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionFU;
        self.form.perf_selection_valid();
        // OCCT L457-466.
        let exp = explorer(from, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L467-478.
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
        // OCCT L479-488.
        let c = test_curve(&self.my_pbase);
        let sens: i32;
        if shape_is_same(from, until) {
            sens = 1;
        } else {
            sens = sens_of_prism(&c, &self.form.my_suntil);
        }
        // OCCT L490-492 (the OCCT source passes myPbase in the Skface
        // parameter slot of HeightMax — source form kept).
        let height = sens as f64
            * height_max(
                &self.form.my_sbase,
                &self.my_pbase,
                &self.form.my_sfrom,
                &self.form.my_suntil,
            );
        let the_d_prism =
            LocOpeDPrism::with_heights_angle(&self.my_pbase, height, height, self.my_angle);
        let vrai_d_prism = the_d_prism.shape();

        // OCCT L494-510.
        if !trff {
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            // Make systematically dprism
            self.form.my_gshape = vrai_d_prism.clone();
            self.form.generated_shape_valid();
            // management of gluing faces
            self.form.glued_faces_valid();
            // OCCT L505-506.
            let mut seq: Vec<Curve3> = Vec::new();
            the_d_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_d_prism.baryc_curve());
            // topologic reconstruction
            global_perform(self);
        } else {
            // management of descendants

            // OCCT L514-530.
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let c1 = if sens == -1 {
                Some(geom_curve_reversed(c.as_ref().expect("test curve")))
            } else {
                c.clone()
            };
            let scur = vec![c1];
            let mut a_si1 = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            let mut a_si2 = LocOpeCSIntersector::with_shape(&self.form.my_sfrom);
            a_si1.perform_cur(&scur);
            a_si2.perform_cur(&scur);
            // OCCT L531-532.
            let mut or_u: Orientation;
            let mut or_f: Orientation;
            let mut f_from = Shape::null();
            let mut f_until = Shape::null();
            // direction of dprism
            // OCCT L534-553.
            if a_si1.is_done() && a_si1.nb_points(1) >= 1 {
                if self.form.my_fuse {
                    or_u = a_si1.point(1, 1).orientation();
                } else {
                    or_u = a_si1.point(1, a_si1.nb_points(1)).orientation();
                }
                f_until = a_si1
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
            } else {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::NoIntersectU;
                return;
            }
            // OCCT L554-566.
            if a_si2.is_done() && a_si2.nb_points(1) >= 1 {
                or_f = a_si2.point(1, 1).orientation();
                or_f = top_abs_reverse(or_f);
                f_from = a_si2
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
            } else {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::NoIntersectF;
                return;
            }
            // OCCT L567-591.
            let mut pool = BRep::new();
            let mut b = BRepBuilder::new();
            let comp = b.make_compound(&mut pool, Vec::new());
            // OCCT L570-580.
            match brep_feat_tool(&self.form.my_suntil, &f_until, or_u) {
                Some(s) => pool_add(&mut pool, &mut b, &comp, &s),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolU;
                    return;
                }
            }
            // OCCT L581-591.
            match brep_feat_tool(&self.form.my_sfrom, &f_from, or_f) {
                Some(ss) => pool_add(&mut pool, &mut b, &comp, &ss),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolF;
                    return;
                }
            }
            // OCCT L593.
            let tr_p =
                CutVehicle::with_operation(&vrai_d_prism, &comp, BooleanOpType::Cut);
            let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
            // OCCT L595-614.
            if self.form.my_fuse {
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
                self.form.my_shape = Some(trp_shape);
                self.form.done();
            }
        }
    }

    /// OCCT BRepFeat_MakeDPrism::PerformUntilEnd (cxx L620-648).
    pub fn perform_until_end(&mut self) {
        // OCCT L627-633.
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionSh;
        self.form.perf_selection_valid();
        self.form.my_glued_f.clear();
        self.form.my_suntil = Shape::null();
        self.form.shape_until_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        // OCCT L634-638.
        let height = height_max(
            &self.form.my_sbase,
            &self.form.my_skface,
            &self.form.my_sfrom,
            &self.form.my_suntil,
        );
        let the_d_prism = LocOpeDPrism::with_height_angle(&self.my_pbase, height, self.my_angle);
        let vrai_d_prism = the_d_prism.shape();
        // OCCT L640-647.
        maj_map(
            &self.my_pbase,
            &the_d_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        self.form.my_gshape = vrai_d_prism;
        self.form.generated_shape_valid();
        self.form.glued_faces_valid();
        let mut seq: Vec<Curve3> = Vec::new();
        the_d_prism.curves(&mut seq);
        self.my_curves = seq.into_iter().map(Some).collect();
        self.my_b_curve = Some(the_d_prism.baryc_curve());
        global_perform(self);
    }

    /// OCCT BRepFeat_MakeDPrism::PerformFromEnd (cxx L655-817) — feature
    /// semiinfinite limited by the shape Until from the other side.
    pub fn perform_from_end(&mut self, until: &Shape) {
        // OCCT L662-670.
        if until.is_null() {
            panic!("Standard_ConstructionError");
        }
        if !self.form.my_skface.is_null() && shape_is_same(until, &self.form.my_skface) {
            self.perform_until_end();
            return;
        }
        // OCCT L671-676.
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L678-684.
        self.form.my_perf_selection = BRepFeatPerfSelection::SelectionShU;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let mut trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L685-687.
        let c = test_curve(&self.my_pbase);
        let sens = sens_of_prism(&c, &self.form.my_suntil);
        let height = sens as f64
            * height_max(
                &self.form.my_sbase,
                &self.form.my_skface,
                &self.form.my_sfrom,
                &self.form.my_suntil,
            );
        // OCCT L689-696.
        let the_d_prism =
            LocOpeDPrism::with_heights_angle(&self.my_pbase, height, height, self.my_angle);
        let vrai_d_prism = the_d_prism.shape();
        if vrai_d_prism.is_null() {
            self.form.not_done();
            self.my_status_error = BRepFeatStatusError::NullRealTool;
            return;
        }

        // OCCT L698-708: case finite face.
        if !trf {
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = vrai_d_prism.clone();
            self.form.generated_shape_valid();
            self.form.my_glued_f.clear();
            self.form.glued_faces_valid();
            // OCCT L705-707.
            let mut seq: Vec<Curve3> = Vec::new();
            the_d_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_d_prism.baryc_curve());
            global_perform(self);
        } else {
            // case support

            // OCCT L711-727.
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let c2 = if sens == -1 {
                Some(geom_curve_reversed(c.as_ref().expect("test curve")))
            } else {
                c.clone()
            };
            let scur = vec![c2];
            let mut a_si1 = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            let mut a_si2 = LocOpeCSIntersector::with_shape(&self.form.my_sbase);
            a_si1.perform_cur(&scur);
            a_si2.perform_cur(&scur);
            // OCCT L728-729.
            let mut or_u = Orientation::Forward;
            let mut or_f = Orientation::Forward;
            let mut f_until = Shape::null();
            let mut f_from = Shape::null();
            // OCCT L730-739.
            if a_si1.is_done() && a_si1.nb_points(1) >= 1 {
                or_u = a_si1.point(1, 1).orientation();
                let prm = a_si1.point(1, 1).parameter();
                if prm < 0.0 {
                    or_u = top_abs_reverse(or_u);
                }
                f_until = a_si1
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
            }
            // OCCT L741-767.
            if a_si2.is_done() && a_si2.nb_points(1) >= 1 {
                let jj = a_si2.nb_points(1);
                let mut prm = a_si2.point(1, 1).parameter();
                f_from = a_si2
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                or_f = a_si2.point(1, 1).orientation();
                or_f = top_abs_reverse(or_f);
                // OCCT L748-757.
                for iii in 1..=jj {
                    if a_si2.point(1, iii).parameter() < prm {
                        prm = a_si2.point(1, iii).parameter();
                        f_from = a_si2
                            .point(1, iii)
                            .face()
                            .cloned()
                            .unwrap_or_else(Shape::null);
                        or_f = a_si2.point(1, iii).orientation();
                        or_f = top_abs_reverse(or_f);
                    }
                }
                // OCCT L758-764: the surface of FFrom, basis; MakeFace (the
                // identity-location reduction).
                let Some(mut s) = brep_tool_surface(&f_from) else {
                    panic!("Geom_Surface null");
                };
                if let Surface3::Trimmed(t) = &s {
                    s = (*t.basis).clone();
                }
                // OCCT L763-764: BRepLib_MakeFace fac(S, Precision::Confusion()).
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                self.form.my_sfrom = b.make_face(&mut pool, Some(s), Shape::null());
                // OCCT L765 (FFrom stays as the intersection face — the
                // assignment is commented out in the OCCT source).
                trf = self.form.transform_shape_fu(0);
            }
            // OCCT L769-794.
            let mut pool = BRep::new();
            let mut b = BRepBuilder::new();
            let comp = b.make_compound(&mut pool, Vec::new());
            // OCCT L772-782.
            match brep_feat_tool(&self.form.my_suntil, &f_until, or_u) {
                Some(sol) => pool_add(&mut pool, &mut b, &comp, &sol),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolU;
                    return;
                }
            }
            // OCCT L784-794.
            match brep_feat_tool(&self.form.my_sfrom, &f_from, or_f) {
                Some(sol1) => pool_add(&mut pool, &mut b, &comp, &sol1),
                None => {
                    self.form.not_done();
                    self.my_status_error = BRepFeatStatusError::NullToolF;
                    return;
                }
            }
            // OCCT L796-815.
            let tr_p =
                CutVehicle::with_operation(&vrai_d_prism, &comp, BooleanOpType::Cut);
            let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
            if self.form.my_fuse {
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
                let c3 = CutVehicle::with_operation(
                    &self.form.my_sbase,
                    &trp_shape,
                    BooleanOpType::Cut,
                );
                self.form.my_shape = c3.shape().cloned();
                let f_shape = self.form.my_shape.clone().unwrap_or_else(Shape::null);
                self.form.update_descendants_bop(&c3, &f_shape, false);
                self.form.done();
            } else {
                self.form.my_shape = Some(trp_shape);
                self.form.done();
            }
        }
    }

    /// OCCT BRepFeat_MakeDPrism::PerformThruAll (cxx L824-873) — feature
    /// throughout the entire initial shape.
    pub fn perform_thru_all(&mut self) {
        // OCCT L831-847.
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
        // OCCT L849-852.
        let height = height_max(
            &self.form.my_sbase,
            &self.form.my_skface,
            &self.form.my_sfrom,
            &self.form.my_suntil,
        );
        let the_d_prism =
            LocOpeDPrism::with_heights_angle(&self.my_pbase, height, height, self.my_angle);
        let vrai_d_prism = the_d_prism.shape();
        maj_map(
            &self.my_pbase,
            &the_d_prism,
            &mut self.form.my_map,
            &mut self.form.my_f_shape,
            &mut self.form.my_l_shape,
        );
        // OCCT L854-855.
        self.form.my_gshape = vrai_d_prism;
        self.form.generated_shape_valid();
        // OCCT L857-872.
        if !self.form.my_fuse {
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
            let mut seq: Vec<Curve3> = Vec::new();
            the_d_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_d_prism.baryc_curve());
            global_perform(self);
        }
    }

    /// OCCT BRepFeat_MakeDPrism::PerformUntilHeight (cxx L880-997).
    pub fn perform_until_height(&mut self, until: &Shape, height: f64) {
        // OCCT L887-894 (no return after the nested Performs — source form
        // kept).
        if until.is_null() {
            self.perform(height);
        }
        if height == 0.0 {
            self.perform_until(until);
        }
        // OCCT L895-899.
        let exp = explorer(until, ShapeType::Face, ShapeType::Shape);
        if exp.is_empty() {
            panic!("Standard_ConstructionError");
        }
        // OCCT L901-908.
        self.form.my_glued_f.clear();
        self.form.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        self.form.perf_selection_valid();
        self.form.my_sfrom = Shape::null();
        self.form.shape_from_valid();
        self.form.my_suntil = until.clone();
        let trf = self.form.transform_shape_fu(1);
        self.form.shape_until_valid();
        // OCCT L909-913.
        let c = test_curve(&self.my_pbase);
        let sens = sens_of_prism(&c, &self.form.my_suntil);
        let the_d_prism =
            LocOpeDPrism::with_height_angle(&self.my_pbase, sens as f64 * height, self.my_angle);
        let vrai_d_prism = the_d_prism.shape();

        // OCCT L915-935: case face finished.
        if !trf {
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            self.form.my_gshape = vrai_d_prism.clone();
            self.form.generated_shape_valid();
            // OCCT L920-929.
            let base = the_d_prism.first_shape();
            let mut expf = explorer(&base, ShapeType::Face, ShapeType::Shape).into_iter();
            let _the_base = expf.next().unwrap_or_else(Shape::null); // OCCT: exp.Current()
            if expf.next().is_some() {
                self.form.not_done();
                self.my_status_error = BRepFeatStatusError::InvFirstShape;
                return;
            }
            // OCCT L931-934.
            self.form.glued_faces_valid();
            let mut seq: Vec<Curve3> = Vec::new();
            the_d_prism.curves(&mut seq);
            self.my_curves = seq.into_iter().map(Some).collect();
            self.my_b_curve = Some(the_d_prism.baryc_curve());
            global_perform(self);
        } else {
            // case support

            // OCCT L938-952.
            maj_map(
                &self.my_pbase,
                &the_d_prism,
                &mut self.form.my_map,
                &mut self.form.my_f_shape,
                &mut self.form.my_l_shape,
            );
            let c1 = if sens == -1 {
                Some(geom_curve_reversed(c.as_ref().expect("test curve")))
            } else {
                c.clone()
            };
            let scur = vec![c1];
            let mut a_si = LocOpeCSIntersector::with_shape(&self.form.my_suntil);
            a_si.perform_cur(&scur);
            // OCCT L953-995.
            if a_si.is_done() && a_si.nb_points(1) >= 1 {
                let mut or_: Orientation;
                if self.form.my_fuse {
                    or_ = a_si.point(1, 1).orientation();
                } else {
                    or_ = a_si.point(1, a_si.nb_points(1)).orientation();
                }
                let f_until = a_si
                    .point(1, 1)
                    .face()
                    .cloned()
                    .unwrap_or_else(Shape::null);
                // OCCT L967-974.
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                let comp = b.make_compound(&mut pool, Vec::new());
                if let Some(s) = brep_feat_tool(&self.form.my_suntil, &f_until, or_) {
                    pool_add(&mut pool, &mut b, &comp, &s);
                }
                // OCCT L975.
                let tr_p =
                    CutVehicle::with_operation(&vrai_d_prism, &comp, BooleanOpType::Cut);
                let trp_shape = tr_p.shape().cloned().unwrap_or_else(Shape::null);
                // OCCT L976-994.
                if self.form.my_fuse {
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
                    self.form.my_shape = Some(trp_shape);
                    self.form.done();
                }
            }
        }
    }

    /// OCCT BRepFeat_MakeDPrism::BossEdges (cxx L1015-1127) — determination
    /// of TopEdges and LatEdges (sig = 1 -> TopEdges = FirstShape of the
    /// DPrism; sig = 2 -> TopEdges = LastShape of the DPrism).
    pub fn boss_edges(&mut self, signature: i32) {
        // OCCT L1022-1035.
        let the_last_shape: Vec<Shape>;
        if signature == 1 || signature == -1 {
            the_last_shape = self.form.first_shape().clone();
        } else if signature == 2 || signature == -2 {
            the_last_shape = self.form.last_shape().clone();
        } else {
            return;
        }

        // Edges Top

        // OCCT L1038-1048.
        for it_ls in &the_last_shape {
            let ff = it_ls;
            for exp_e in explorer(ff, ShapeType::Edge, ShapeType::Shape) {
                let ee = exp_e;
                self.my_top_edges.push(ee);
            }
        }

        // Edges Bottom

        // OCCT L1051-1055.
        if signature < 0 {
            // Attention check if TgtEdges is important
            self.my_lat_edges = self.form.new_edges().clone();
        } else if signature > 0 {
            // OCCT L1056-1126.
            let Some(my_shape) = self.form.my_shape.clone() else {
                return;
            };
            if !my_shape.is_null() {
                // OCCT L1060: NCollection_Map MapE (architecture difference
                // #5: the bucket iteration order is not reproduced).
                let mut map_e: HashMap<(u64, u32), Shape> = HashMap::new();

                // OCCT L1063-1106.
                for exp_f in explorer(&my_shape, ShapeType::Face, ShapeType::Shape) {
                    let mut found = false;
                    let ff = exp_f;
                    for it_ls in &the_last_shape {
                        let top_face = it_ls;
                        if !shape_is_same(&ff, top_face) {
                            let exp_edges =
                                explorer(&ff, ShapeType::Edge, ShapeType::Shape);
                            for exp_e in exp_edges {
                                if found {
                                    break;
                                }
                                let e1 = exp_e;
                                let (v1, v2) = top_exp_vertices(&e1);
                                for it in self.my_top_edges.clone() {
                                    if found {
                                        break;
                                    }
                                    let e2 = it;
                                    let (vt1, vt2) = top_exp_vertices(&e2);
                                    if shape_is_same(&v1, &vt1)
                                        || shape_is_same(&v1, &vt2)
                                        || shape_is_same(&v2, &vt1)
                                        || shape_is_same(&v2, &vt2)
                                    {
                                        found = true;
                                        for exp_e2 in
                                            explorer(&ff, ShapeType::Edge, ShapeType::Shape)
                                        {
                                            let e3 = exp_e2;
                                            if map_e.contains_key(&shape_key(&e3)) {
                                                map_e.remove(&shape_key(&e3));
                                            } else {
                                                map_e.insert(shape_key(&e3), e3);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L1108-1115.
                for it in self.my_top_edges.clone() {
                    if map_e.contains_key(&shape_key(&it)) {
                        map_e.remove(&shape_key(&it));
                    }
                }

                // OCCT L1117-1124.
                for (_k, itm) in map_e.iter() {
                    if !brep_tool_degenerated(itm) {
                        self.my_lat_edges.push(itm.clone());
                    }
                }
            }
        }
    }

    /// OCCT BRepFeat_MakeDPrism::TopEdges (cxx L1133-1136) — the list of
    /// edges of the top of the boss.
    pub fn top_edges(&self) -> &Vec<Shape> {
        &self.my_top_edges
    }

    /// OCCT BRepFeat_MakeDPrism::LatEdges (cxx L1142-1145) — the list of
    /// edges of the bottom of the boss.
    pub fn lat_edges(&self) -> &Vec<Shape> {
        &self.my_lat_edges
    }
}

/// OCCT BRep_Builder::Add on a compound over the pool (the local carrier
/// shared by the Comp constructions).
fn pool_add(
    pool: &mut BRep,
    b: &mut BRepBuilder,
    comp: &Shape,
    s: &Shape,
) {
    b.add_to_compound(pool, comp.clone(), s.clone());
}

/// OCCT BRepFeat_Form::Curves override (cxx L1004-1007).
impl BRepFeatFormSlots for BRepFeatMakeDPrism {
    fn curves(&mut self, s: &mut Vec<Option<Curve3>>) {
        // OCCT L1006: scur = myCurves.
        *s = self.my_curves.clone();
    }

    /// OCCT BRepFeat_Form::BarycCurve override (cxx L1152-1155) — passe par
    /// le centre de masses de la primitive.
    fn baryc_curve(&mut self) -> Option<Curve3> {
        // OCCT L1154: return myBCurve.
        self.my_b_curve.clone()
    }

    /// The BRepFeat_Form base sub-object.
    fn form(&mut self) -> &mut BRepFeatForm {
        &mut self.form
    }
}

/// OCCT static HeightMax (cxx L1162-1198) — calculate the height of the
/// prism following the parameters of the bounding box (the DPrism variant:
/// the largest span, no FacRevolInfini guards).
fn height_max(
    the_sbase: &Shape,
    the_skface: &Shape,
    the_sfrom: &Shape,
    the_suntil: &Shape,
) -> f64 {
    // OCCT L1167-1177 (the identity-location reduction of the per-shape
    // box; the accumulation is the union).
    let mut box_min = glam::DVec3::splat(f64::MAX);
    let mut box_max = glam::DVec3::splat(f64::MIN);
    add_to_box(the_sbase, &mut box_min, &mut box_max);
    add_to_box(the_skface, &mut box_min, &mut box_max);
    if !the_sfrom.is_null() {
        add_to_box(the_sfrom, &mut box_min, &mut box_max);
    }
    if !the_suntil.is_null() {
        add_to_box(the_suntil, &mut box_min, &mut box_max);
    }
    // OCCT L1178-1180: c = (xmin, xmax, ymin, ymax, zmin, zmax).
    let c = [
        box_min.x,
        box_max.x,
        box_min.y,
        box_max.y,
        box_min.z,
        box_max.z,
    ];
    // OCCT L1189.
    let par = (c[1] - c[0])
        .abs()
        .max((c[3] - c[2]).abs())
        .max((c[5] - c[4]).abs());
    // OCCT L1197.
    par
}

/// The BRepBndLib::Add carrier of height_max (the rcad shape_box union).
fn add_to_box(the_s: &Shape, box_min: &mut DVec3, box_max: &mut DVec3) {
    if let Some((mn, mx)) = BRepFeatBuilder::shape_box(the_s, &[]) {
        *box_min = box_min.min(mn);
        *box_max = box_max.max(mx);
    }
}

/// OCCT static SensOfPrism (cxx L1204-1229) — determine the direction of
/// prism generation (the DPrism variant: no Precision::Confusion()).
fn sens_of_prism(the_c: &Option<Curve3>, until: &Shape) -> i32 {
    // OCCT L1206-1209.
    let mut a_si1 = LocOpeCSIntersector::with_shape(until);
    let scur = vec![the_c.clone()];
    a_si1.perform_cur(&scur);
    let mut sens = 1;
    // OCCT L1211-1220.
    if a_si1.is_done() && a_si1.nb_points(1) >= 1 {
        let nb = a_si1.nb_points(1);
        let prm1 = a_si1.point(1, 1).parameter();
        let prm2 = a_si1.point(1, nb).parameter();
        if prm1 < 0.0 && prm2 < 0.0 {
            sens = -1;
        }
    } else {
        // OCCT L1221-1224: BRepFeat::ParametricBarycenter(Until, C) — C is
        // the TestCurve result (the null-handle case has no OCCT
        // counterpart for the non-planar base).
        let Some(c) = the_c else {
            panic!("null test curve (no OCCT counterpart)");
        };
        if brep_feat_parametric_barycenter(until, c) < 0.0 {
            sens = -1;
        }
    }
    sens
}

/// OCCT static MajMap (cxx L1233-1282) — the DPrism variant (the null
/// FirstShape/LastShape guards).
fn maj_map(
    the_b: &Shape,
    the_p: &LocOpeDPrism,
    the_map: &mut HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_f_shape: &mut Shape,
    the_l_shape: &mut Shape,
) {
    // OCCT L1241-1255.
    if !the_p.first_shape().is_null() {
        let exp = explorer(&the_p.first_shape(), ShapeType::Wire, ShapeType::Shape);
        if let Some(cur) = exp.first() {
            *the_f_shape = cur.clone();
            the_map.insert(shape_key(the_f_shape), (the_f_shape.clone(), Vec::new()));
            for exp in
                explorer(&the_p.first_shape(), ShapeType::Face, ShapeType::Shape)
            {
                the_map
                    .get_mut(&shape_key(the_f_shape))
                    .expect("theMap(theFShape)")
                    .1
                    .push(exp);
            }
        }
    }
    // OCCT L1257-1270.
    if !the_p.last_shape().is_null() {
        let exp = explorer(&the_p.last_shape(), ShapeType::Wire, ShapeType::Shape);
        if let Some(cur) = exp.first() {
            *the_l_shape = cur.clone();
            the_map.insert(shape_key(the_l_shape), (the_l_shape.clone(), Vec::new()));
            for exp in
                explorer(&the_p.last_shape(), ShapeType::Face, ShapeType::Shape)
            {
                the_map
                    .get_mut(&shape_key(the_l_shape))
                    .expect("theMap(theLShape)")
                    .1
                    .push(exp);
            }
        }
    }
    // OCCT L1272-1281.
    for exp in explorer(the_b, ShapeType::Edge, ShapeType::Shape) {
        if !the_map.contains_key(&shape_key(&exp)) {
            let shapes = the_p.shapes(&exp);
            the_map.insert(shape_key(&exp), (exp.clone(), shapes));
        }
    }
}

/// OCCT static TestCurve(const TopoDS_Face&) (cxx L1286-1313) — the line
/// through the barycentre of the sampled edges along the plane normal; the
/// OCCT null handle (non-planar base) maps to None.
fn test_curve(the_base: &Shape) -> Option<Curve3> {
    // OCCT L1288-1290.
    let mut bar = DVec3::ZERO;
    let mut spt: Vec<DVec3> = Vec::new();
    crate::feat::loc_ope::sample_edges(the_base, &mut spt);
    // OCCT L1291-1296.
    for pvt in &spt {
        bar += *pvt;
    }
    bar /= spt.len() as f64;
    // OCCT L1297-1301: the surface of the face, basis of a trimmed surface.
    let Some(mut s) = brep_tool_surface(the_base) else {
        return None;
    };
    if let Surface3::Trimmed(t) = &s {
        s = (*t.basis).clone();
    }
    // OCCT L1302-1307: the plane extraction; the null handle otherwise.
    let Surface3::Plane(p) = &s else {
        return None;
    };
    // OCCT L1308-1311: Normale = XDirection ^ YDirection; Geom_Line(bar,
    // Normale). The rcad Plane carries the normal invariant.
    let normale = p.normal;
    Some(Curve3::Line(Line3::new(bar, normale)))
}

#[cfg(test)]
mod tests {
    //! Translation-period placeholder: anchor tests are a stage-2 asset
    //! (acceptance = cargo check + formal alignment, no test runs).

    #[test]
    fn placeholder() {}
}
