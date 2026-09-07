// OCCT BRepFeat_MakeCylindricalHole.cxx L1-811 + BRepFeat_MakeCylindricalHole.hxx
// L1-117 + BRepFeat_MakeCylindricalHole.lxx L1-52 + BRepFeat_Status.hxx L1-28 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeCylindricalHole.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeCylindricalHole.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeCylindricalHole.lxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_Status.hxx
//
// OCCT inheritance chain (BRepFeat_MakeCylindricalHole.hxx L34):
//   BRepFeat_MakeCylindricalHole : BRepFeat_Builder : BOPAlgo_BOP : ...
// Rust has no inheritance -> composition + delegation: the BRepFeat_Builder
// base sub-object is rcad's BRepFeatBuilder
// (crate::feat::brep_feat_builder::BRepFeatBuilder), reached through
// `self.base`; every inherited member/method keeps its OCCT name shape.
//
// Architecture differences (referenced from the affected functions):
// 1. LocOpe_CurveShapeIntersector / LocOpe_PntFace (TKFeat/LocOpe) are not
//    translated yet (the loc_ope_* module is a later stage of the port plan).
//    Their API surface is carried below exactly as OCCT declares it
//    (LocOpe_CurveShapeIntersector.hxx: IsDone/NbPoints/LocalizeAfter/
//    LocalizeBefore/Point; LocOpe_PntFace.hxx: Pnt/Face/Parameter/
//    UParameter/VParameter); the intersection body lands with the LocOpe
//    stage. Until then IsDone() reports false, so the OCCT
//    "if (!theASI.IsDone())" guards take the BRepFeat_InvalidPlacement early
//    return — the same guard structure, driven by the missing dependency.
// 2. BRepPrim_Cylinder (TKPrim/BRepPrim) is not translated yet; the cylinder
//    shell/cap-face production stays null until the TKPrim stage lands. The
//    tool solid of every Perform* therefore cannot be built yet and the
//    BOPAlgo_BOP::Perform step is skipped (guarded) rather than run with an
//    empty tool.
// 3. GetOffset (L749-776) needs BRepAdaptor_Surface::D1 + CSLib::Normal
//    (surface derivative evaluation); until that support lands it reports
//    "not defined", and the OCCT own fallback (offF/offL = Radius) applies.
// 4. BRepFeat_Status / BRepFeat_StatusError (BRepFeat_Status.hxx L20-25 /
//    BRepFeat_StatusError.hxx L21-51) are carried locally below with their
//    OCCT anchors (no dedicated rcad module exists for them yet).
// 5. Standard_ConstructionError throws (L61, L111, L262, L365, L481) map to
//    panic! at the same conditions (rcad has no exception machinery on this
//    facade).
// 6. OCCT out-parameters (PartsOfTool(aLT), GetOffset(..., outOff), CreateCyl
//    (... Cyl, CylTopF, CylBottF)) map to return values / Option<f64>.

use crate::feat::brep_feat_builder::{pool_faces, BRepFeatBuilder};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods;

/// OCCT BRepFeat_Status (BRepFeat_Status.hxx L20-25) — carried locally
/// (architecture difference #4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepFeatStatus {
    /// BRepFeat_NoError
    NoError,
    /// BRepFeat_InvalidPlacement
    InvalidPlacement,
    /// BRepFeat_HoleTooLong
    HoleTooLong,
}

/// OCCT BRepFeat_StatusError (BRepFeat_StatusError.hxx L21-51) — carried
/// locally (architecture difference #4); used by the BRepFeat_Form family
/// translated in a later stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum BRepFeatStatusError {
    OK,
    BadDirect,
    BadIntersect,
    EmptyBaryCurve,
    EmptyCutResult,
    FalseSide,
    IncDirection,
    IncSlidFace,
    IncParameter,
    IncTypes,
    IntervalOverlap,
    InvFirstShape,
    InvOption,
    InvShape,
    LocOpeNotDone,
    LocOpeInvNotDone,
    NoExtFace,
    NoFaceProf,
    NoGluer,
    NoIntersectF,
    NoIntersectU,
    NoParts,
    NoProjPt,
    NotInitialized,
    NotYetImplemented,
    NullRealTool,
    NullToolF,
    NullToolU,
}

// OCCT gp_Ax1 -> rcad_kernel::math::gp::Ax1;  OCCT LocOpe_PntFace and
// LocOpe_CurveShapeIntersector -> the loc_ope_* translations (Stage 3c
// front half).  The Stage 3a local stubs were retired when the real
// modules landed; the int-overload quirk (LocalizeBefore(int) implicitly
// converts to the Real overload, cxx L294 -> L381) is replicated inside
// `LocOpeCurveShapeIntersector::localize_before_index`.
use rcad_kernel::math::gp::Ax1;
use super::loc_ope_curve_shape_intersector::LocOpeCurveShapeIntersector;
use super::loc_ope_pnt_face::LocOpePntFace;

/// OCCT BRepPrim_Cylinder (TKPrim/BRepPrim/BRepPrim_Cylinder.cxx) — a
/// cylinder solid primitive. API surface carried per OCCT (Shell/TopFace/
/// BottomFace); the shell production lands with the TKPrim stage
/// (architecture difference #2).
#[allow(dead_code)]
pub(crate) struct BRepPrimCylinder {
    // OCCT BRepPrim_OneAxis fields (axis, radius, height).
    pub(crate) my_location: glam::DVec3,
    pub(crate) my_direction: glam::DVec3,
    pub(crate) my_radius: f64,
    pub(crate) my_height: f64,
}

impl BRepPrimCylinder {
    /// OCCT BRepPrim_Cylinder::BRepPrim_Cylinder(Axis, Radius, Height).
    pub fn new(a_location: glam::DVec3, a_direction: glam::DVec3, radius: f64, height: f64) -> Self {
        BRepPrimCylinder {
            my_location: a_location,
            my_direction: a_direction,
            my_radius: radius,
            my_height: height,
        }
    }
    /// OCCT BRepPrim_OneAxis::Shell() — null until the TKPrim stage lands.
    pub fn shell(&self) -> Option<Shape> {
        None
    }
    /// OCCT BRepPrim_OneAxis::TopFace() — null until the TKPrim stage lands.
    pub fn top_face(&self) -> Option<Shape> {
        None
    }
    /// OCCT BRepPrim_OneAxis::BottomFace() — null until the TKPrim stage
    /// lands.
    pub fn bottom_face(&self) -> Option<Shape> {
        None
    }
}

/// OCCT static Baryc (MakeCylindricalHole.cxx L690-717) — the barycentre of
/// the non-degenerated edges of S, sampled with 11 points per edge. OCCT
/// fills the gp_Pnt out-parameter; rcad returns the point.
fn baryc(s: &Shape, brep: Option<&topods::BRep>) -> glam::DVec3 {
    use rcad_kernel::CurveEval;
    //
    let mut bar = glam::DVec3::ZERO;
    let mut nbp = 0i32;
    for a_e in crate::feat::brep_feat_builder::explorer(s, topods::ShapeType::Edge, topods::ShapeType::Shape) {
        // OCCT L702: const TopoDS_Edge& E = TopoDS::Edge(exp.Current());
        // Calculate points by non-degenerated edges (L701-703).
        let Some(ed) = a_e.as_edge() else { continue };
        if ed.degenerated {
            continue;
        }
        // OCCT L705: C = BRep_Tool::Curve(E, L, First, Last);
        let Some(c) = ed.curve.clone() else { continue };
        let [first, last] = ed.range;
        // OCCT L706: C = C->Transformed(L.Transformation());
        let loc = brep
            .map(|b| b.get_location(a_e.location))
            .unwrap_or(glam::DAffine3::IDENTITY);
        let c = rcad_kernel::geom::transform_curve(&c, &loc);
        // OCCT L707-712: 11 sample points per edge.
        for i in 1..=11 {
            let prm = ((11 - i) as f64 * first + (i - 1) as f64 * last) / 10.;
            bar += c.point_at(prm);
            nbp += 1;
        }
    }
    // OCCT L715-716: Bar.Divide((double)nbp); B.SetXYZ(Bar);
    bar /= nbp as f64;
    bar
}

/// OCCT static BoxParameters (MakeCylindricalHole.cxx L719-747) — the
/// parameters of a bounding box in the direction of the axis of the hole.
/// OCCT fills parmin/parmax; rcad returns the pair.
fn box_parameters(s: &Shape, axis: &Ax1) -> (f64, f64) {
    //
    // OCCT L723-726: Bnd_Box B; BRepBndLib::Add(S, B);
    // B.Get(c[0], c[2], c[4], c[1], c[3], c[5]);
    // Architecture gap: the argument shape's own TopLoc_Location table is not
    // carried by the facade, so the box is evaluated at the identity
    // location (the Init argument is expected unlocated until the location
    // table plumbing lands).
    let (c_xmin, c_ymin, c_zmin, c_xmax, c_ymax, c_zmax) =
        match crate::feat::brep_feat_builder::BRepFeatBuilder::shape_box(s, &[]) {
            Some((mn, mx)) => (mn.x, mn.y, mn.z, mx.x, mx.y, mx.z),
            None => return (f64::MAX, f64::MIN),
        };
    // OCCT stores c[0..2] = (xmin, xmax), c[2..4] = (ymin, ymax),
    // c[4..6] = (zmin, zmax).
    let c = [c_xmin, c_xmax, c_ymin, c_ymax, c_zmin, c_zmax];
    let mut p = glam::DVec3::ZERO;
    let mut parmin = f64::MAX; // OCCT: RealLast()
    let mut parmax = f64::MIN; // OCCT: RealFirst()
    let mut param: f64;
    for i in 0..=1usize {
        p.x = c[i];
        for j in 2..=3usize {
            p.y = c[j];
            for k in 4..=5usize {
                p.z = c[k];
                // OCCT L741: param = ElCLib::LineParameter(Axis, P);
                param =
                    crate::geomalgo::int_patch::elclib::line_parameter_of_axis(
                        axis.location,
                        axis.direction,
                        p,
                    );
                // OCCT L742-743.
                parmin = parmin.min(param);
                parmax = parmax.max(param);
            }
        }
    }
    (parmin, parmax)
}

/// OCCT static GetOffset (MakeCylindricalHole.cxx L749-776) — the offset of
/// the cylinder cap along the axis from the intersection point on the face.
/// OCCT returns bool with the outOff out-parameter; rcad returns Option<f64>
/// (None == OCCT false).
fn get_offset(pnt_info: &LocOpePntFace, radius: f64, axis: &Ax1) -> Option<f64> {
    // OCCT L754-755: FF = PntInfo.Face(); BRepAdaptor_Surface FFA(FF);
    let _ff = pnt_info.face();
    let _radius = radius;
    let _axis = axis;
    // OCCT L757-764: Up/Vp = PntInfo.U/VParameter(); FFA.D1(Up, Vp, PP, D1U,
    // D1V); CSLib::Normal(D1U, D1V, Precision::Angular(), stat, NormF);
    // Architecture difference #3: BRepAdaptor_Surface::D1 + CSLib::Normal
    // (surface derivative evaluation) are not available yet.
    // OCCT L765-775: stat != CSLib_Defined -> false; angle near PI/2 ->
    // false; outOff = Radius * |tan(angle)|.
    None
}

/// Result of OCCT static CreateCyl (MakeCylindricalHole.cxx L778-810): the
/// cylinder shell and its cap faces (OCCT out-parameters Cyl, CylTopF,
/// CylBottF).
pub(crate) struct CreatedCylinder {
    pub cyl: Option<Shape>,      // OCCT: Cyl (TopoDS_Shell)
    pub cyl_top_f: Option<Shape>, // OCCT: CylTopF
    pub cyl_bott_f: Option<Shape>, // OCCT: CylBottF
}

/// OCCT static CreateCyl (MakeCylindricalHole.cxx L778-810) — creates the
/// cylinder along the axis from 'First - offF' to 'Last + offL' params.
#[allow(unused_assignments)]
fn create_cyl(
    pnt_info_first: &LocOpePntFace,
    pnt_info_last: &LocOpePntFace,
    radius: f64,
    axis: &Ax1,
) -> CreatedCylinder {
    let mut off_f = 0.0f64; // OCCT L787: offF = 0.
    let mut off_l = 0.0f64; // OCCT L787: offL = 0.
    // OCCT L788-790: Last/First/Heigth.
    let last = pnt_info_last.parameter();
    let first = pnt_info_first.parameter();
    let heigth = last - first;
    //
    // OCCT L792-799: if (!GetOffset(...)) offF/offL = Radius;
    if let Some(off) = get_offset(pnt_info_first, radius, axis) {
        off_f = off;
    } else {
        off_f = radius;
    }
    if let Some(off) = get_offset(pnt_info_last, radius, axis) {
        off_l = off;
    } else {
        off_l = radius;
    }
    //
    // OCCT L803: theOrig = PntInfoFirst.Pnt().XYZ() - offF *
    // Axis.Direction().XYZ();
    let the_orig = pnt_info_first.pnt() - off_f * axis.direction;
    // OCCT L804-806: gp_Pnt p2_ao1(theOrig); gp_Ax2 a2_ao1(p2_ao1,
    // Axis.Direction()); BRepPrim_Cylinder theCylinder(a2_ao1, Radius,
    // Heigth + offF + offL);
    let the_cylinder = BRepPrimCylinder::new(the_orig, axis.direction, radius, heigth + off_f + off_l);
    // OCCT L807-809: Cyl = theCylinder.Shell(); CylTopF =
    // theCylinder.TopFace(); CylBottF = theCylinder.BottomFace();
    // Architecture difference #2: the shell production lands with TKPrim.
    CreatedCylinder {
        cyl: the_cylinder.shell(),
        cyl_top_f: the_cylinder.top_face(),
        cyl_bott_f: the_cylinder.bottom_face(),
    }
}

/// OCCT BRepFeat_MakeCylindricalHole — provides a tool to make cylindrical
/// holes on a shape (BRepFeat_MakeCylindricalHole.hxx L33).
pub struct BRepFeatMakeCylindricalHole {
    /// The BRepFeat_Builder base sub-object (composition + delegation).
    pub(crate) base: BRepFeatBuilder,
    my_axis: Ax1,                     // OCCT L105: myAxis
    my_ax_def: bool,                    // OCCT L106: myAxDef
    my_status: BRepFeatStatus,          // OCCT L107: myStatus
    my_is_blind: bool,                  // OCCT L108: myIsBlind
    my_validate: bool,                  // OCCT L109: myValidate
    my_top_face: Option<Shape>,         // OCCT L110: myTopFace
    my_bot_face: Option<Shape>,         // OCCT L111: myBotFace
    /// OCCT myShape = Shape() in Build() — the result of the operation. rcad
    /// keeps the result pool (BRepFeatBuilder::myShape).
    my_result: Option<topods::BRep>,
}

impl BRepFeatMakeCylindricalHole {
    /// OCCT BRepFeat_MakeCylindricalHole::BRepFeat_MakeCylindricalHole()
    /// (lxx L21-27): myAxDef(false), myStatus(BRepFeat_NoError),
    /// myIsBlind(false), myValidate(false).
    pub fn new() -> Self {
        BRepFeatMakeCylindricalHole {
            base: BRepFeatBuilder::new(),
            my_axis: Ax1 {
                location: glam::DVec3::ZERO,
                direction: glam::DVec3::Z,
            },
            my_ax_def: false,
            my_status: BRepFeatStatus::NoError,
            my_is_blind: false,
            my_validate: false,
            my_top_face: None,
            my_bot_face: None,
            my_result: None,
        }
    }

    /// OCCT BRepFeat_MakeCylindricalHole::Init(const gp_Ax1& Axis) (lxx
    /// L31-35) — sets the axis of the hole(s) (the Init overload without a
    /// shape; Rust has no overloading).
    pub fn init_axis(&mut self, axis: Ax1) {
        self.my_axis = axis;
        self.my_ax_def = true;
    }

    /// OCCT BRepFeat_MakeCylindricalHole::Init(const TopoDS_Shape& S, const
    /// gp_Ax1& Axis) (lxx L39-44) — sets the shape and axis on which hole(s)
    /// will be performed.
    pub fn init(&mut self, s: &Shape, axis: Ax1) {
        self.base.init(s);
        self.my_axis = axis;
        self.my_ax_def = true;
    }

    /// OCCT BRepFeat_MakeCylindricalHole::Status() (lxx L48-51) — returns the
    /// status after a hole is performed.
    pub fn status(&self) -> BRepFeatStatus {
        self.my_status
    }

    /// OCCT L58-62 / L108-112 / L259-263 / L362-366 / L479-483: the argument
    /// fetch and the Standard_ConstructionError guard
    /// (aObject.IsNull() || !myAxDef).
    fn object_arg(&self) -> Shape {
        match self.base.my_arguments.first() {
            Some(s) if self.my_ax_def => s.clone(),
            _ => panic!("Standard_ConstructionError"), // OCCT throw (difference #5)
        }
    }

    /// OCCT BRepFeat_MakeCylindricalHole::Perform(const double Radius)
    /// (L56-101) — performs every hole of radius Radius (a cut with an
    /// infinite cylinder defined by the axis and Radius).
    pub fn perform(&mut self, radius: f64) {
        // OCCT L58-62.
        let a_object = self.object_arg();
        //
        self.my_is_blind = false; // OCCT L64
        self.my_status = BRepFeatStatus::NoError; // OCCT L65
        //
        // OCCT L67-72: LocOpe_CurveShapeIntersector theASI(myAxis, aObject);
        // if (!theASI.IsDone() || theASI.NbPoints() <= 0) { myStatus =
        // BRepFeat_InvalidPlacement; return; }
        let the_asi = LocOpeCurveShapeIntersector::new_axis(&self.my_axis, &a_object);
        if !the_asi.is_done() || the_asi.nb_points() <= 0 {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // It is not possible to use infinite cylinder for topological
        // operations. (OCCT L74-82)
        let (p_min, p_max) = box_parameters(&a_object, &self.my_axis);
        let heigth = 2. * (p_max - p_min);
        let the_orig =
            self.my_axis.location + ((3. * p_min - p_max) / 2.) * self.my_axis.direction;
        // OCCT L80-82: gp_Pnt p1_ao1(theOrig); gp_Ax2 a1_ao1(p1_ao1,
        // myAxis.Direction()); BRepPrim_Cylinder theCylinder(a1_ao1, Radius,
        // Heigth);
        let the_cylinder = BRepPrimCylinder::new(the_orig, self.my_axis.direction, radius, heigth);
        //
        // Probably it is better to make cut directly (OCCT L84).
        // OCCT L86-89: BRep_Builder B; TopoDS_Solid theTool;
        // B.MakeSolid(theTool); B.Add(theTool, theCylinder.Shell());
        // Architecture difference #2: the tool solid stays null until TKPrim
        // produces the cylinder shell.
        let the_tool: Option<Shape> = the_cylinder.shell();
        //
        self.my_top_face = the_cylinder.top_face(); // OCCT L91
        self.my_bot_face = the_cylinder.bottom_face(); // OCCT L92
        self.my_validate = false; // OCCT L93
        //
        // OCCT L96-100: bool Fuse = false; AddTool(theTool);
        // SetOperation(Fuse); BOPAlgo_BOP::Perform();
        let fuse = 0i32;
        if let Some(the_tool) = the_tool {
            self.base.add_tool(the_tool);
            self.base.set_operation(fuse);
            self.base.perform_bop();
        }
    }

    /// OCCT BRepFeat_MakeCylindricalHole::PerformThruNext(const double
    /// Radius, const bool Cont) (L105-252) — performs the first hole of
    /// radius Radius in the direction of the defined axis.
    pub fn perform_thru_next(&mut self, radius: f64, cont: bool) {
        //
        // OCCT L108-112.
        let a_object = self.object_arg();
        //
        self.my_is_blind = false; // OCCT L114
        self.my_validate = cont; // OCCT L115
        self.my_status = BRepFeatStatus::NoError; // OCCT L116
        //
        // OCCT L118-123.
        let the_asi = LocOpeCurveShapeIntersector::new_axis(&self.my_axis, &a_object);
        if !the_asi.is_done() {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L125-163: the localization of the first two crossings.
        let mut ind_from: i32 = 0;
        let mut ind_to: i32 = 0;
        let mut the_or = topods::Orientation::Forward;
        let mut pnt_info_first = LocOpePntFace::new();
        let mut pnt_info_last = LocOpePntFace::new();
        let mut ok = the_asi.localize_after(0., &mut the_or, &mut ind_from, &mut ind_to);
        if ok {
            if the_or == topods::Orientation::Forward {
                pnt_info_first = the_asi.point(ind_from);
                // OCCT L134: LocalizeAfter(IndTo, ...) — the index overload.
                ok = the_asi.localize_after_index(ind_to, &mut the_or, &mut ind_from, &mut ind_to);
                if ok {
                    if the_or != topods::Orientation::Reversed {
                        ok = false;
                    } else {
                        pnt_info_last = the_asi.point(ind_to);
                    }
                }
            } else {
                // TopAbs_REVERSED
                pnt_info_last = the_asi.point(ind_to);
                ok = the_asi.localize_before_index(ind_from, &mut the_or, &mut ind_from, &mut ind_to);
                if ok {
                    if the_or != topods::Orientation::Forward {
                        ok = false;
                    } else {
                        pnt_info_first = the_asi.point(ind_from);
                    }
                }
            }
        }
        if !ok {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L170-176: TopoDS_Shell Cyl; CreateCyl(...); BRep_Builder B;
        // TopoDS_Solid theTool; B.MakeSolid(theTool); B.Add(theTool, Cyl);
        let the_cyl = create_cyl(&pnt_info_first, &pnt_info_last, radius, &self.my_axis);
        self.my_top_face = the_cyl.cyl_top_f.clone();
        self.my_bot_face = the_cyl.cyl_bott_f.clone();
        // Architecture difference #2: the tool solid stays null until TKPrim
        // produces the cylinder shell.
        let the_tool: Option<Shape> = the_cyl.cyl.clone();
        //
        // OCCT L178-181: bool Fuse = false; AddTool(theTool);
        // SetOperation(Fuse); BOPAlgo_BOP::Perform();
        let fuse = 0i32;
        if let Some(the_tool) = the_tool {
            self.base.add_tool(the_tool);
            self.base.set_operation(fuse);
            self.base.perform_bop();
            // OCCT L182-251: PartsOfTool(parts) and the part selection.
            let parts = self.base.parts_of_tool();
            //
            let mut nbparts = 0i32;
            for _its in &parts {
                nbparts += 1;
            }
            if nbparts == 0 {
                self.my_status = BRepFeatStatus::InvalidPlacement;
                return;
            }
            //
            if nbparts >= 2 {
                // preserve the smallest as parameter along the axis (OCCT
                // L197-251).
                let first = pnt_info_first.parameter();
                let last = pnt_info_last.parameter();
                let mut tokeep: Option<Shape> = None;
                let mut parmin = last;
                for its in &parts {
                    // OCCT L208-209: Baryc(...); parbar =
                    // ElCLib::LineParameter(myAxis, Barycentre);
                    let barycentre = baryc(its, self.base.my_shape.as_ref());
                    let parbar = crate::geomalgo::int_patch::elclib::line_parameter_of_axis(
                        self.my_axis.location,
                        self.my_axis.direction,
                        barycentre,
                    );
                    if parbar >= first && parbar <= last && parbar <= parmin {
                        parmin = parbar;
                        tokeep = Some(its.clone());
                    }
                }
                //
                if tokeep.is_none() {
                    // preserve the closest interval (OCCT L218-242).
                    let mut dmin = f64::MAX; // OCCT: RealLast()
                    for its in &parts {
                        let barycentre = baryc(its, self.base.my_shape.as_ref());
                        let parbar = crate::geomalgo::int_patch::elclib::line_parameter_of_axis(
                            self.my_axis.location,
                            self.my_axis.direction,
                            barycentre,
                        );
                        if parbar < first {
                            if first - parbar < dmin {
                                dmin = first - parbar;
                                tokeep = Some(its.clone());
                            }
                        } else {
                            // parbar > Last (OCCT L233)
                            if parbar - last < dmin {
                                dmin = parbar - last;
                                tokeep = Some(its.clone());
                            }
                        }
                    }
                }
                for its in &parts {
                    let is_same = match &tokeep {
                        Some(tk) => shape_is_same(tk, its),
                        None => false,
                    };
                    if is_same {
                        self.base.keep_part(its);
                        break;
                    }
                }
            }
        }
    }

    /// OCCT BRepFeat_MakeCylindricalHole::PerformUntilEnd(const double
    /// Radius, const bool Cont) (L256-352) — performs every hole of radius
    /// Radius located after the origin of the given axis.
    pub fn perform_until_end(&mut self, radius: f64, cont: bool) {
        //
        // OCCT L259-263.
        let a_object = self.object_arg();
        //
        self.my_is_blind = false; // OCCT L265
        self.my_validate = cont; // OCCT L266
        self.my_status = BRepFeatStatus::NoError; // OCCT L267
        //
        // OCCT L269-274.
        let the_asi = LocOpeCurveShapeIntersector::new_axis(&self.my_axis, &a_object);
        if !the_asi.is_done() {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L276-304: the localization.
        let mut ind_from: i32 = 0;
        let mut ind_to: i32 = 0;
        let mut the_or = topods::Orientation::Forward;
        let mut pnt_info_first = LocOpePntFace::new();
        let mut pnt_info_last = LocOpePntFace::new();
        let mut ok = the_asi.localize_after(0., &mut the_or, &mut ind_from, &mut ind_to);
        if ok {
            if the_or == topods::Orientation::Reversed {
                // on reset (OCCT L285)
                ok = the_asi.localize_before_index(ind_from, &mut the_or, &mut ind_from, &mut ind_to);
                // It is possible to search for the next.
            }
            if ok && the_or == topods::Orientation::Forward {
                pnt_info_first = the_asi.point(ind_from);
                ok = the_asi.localize_before_index(
                    the_asi.nb_points() + 1,
                    &mut the_or,
                    &mut ind_from,
                    &mut ind_to,
                );
                if ok {
                    if the_or != topods::Orientation::Reversed {
                        ok = false;
                    } else {
                        pnt_info_last = the_asi.point(ind_to);
                    }
                }
            }
        }
        if !ok {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L311-317.
        let the_cyl = create_cyl(&pnt_info_first, &pnt_info_last, radius, &self.my_axis);
        self.my_top_face = the_cyl.cyl_top_f.clone();
        self.my_bot_face = the_cyl.cyl_bott_f.clone();
        let the_tool: Option<Shape> = the_cyl.cyl.clone();
        //
        // OCCT L319-322.
        let fuse = 0i32;
        if let Some(the_tool) = the_tool {
            self.base.add_tool(the_tool);
            self.base.set_operation(fuse);
            self.base.perform_bop();
            // OCCT L323-324.
            let parts = self.base.parts_of_tool();
            //
            let mut nbparts = 0i32;
            for _its in &parts {
                nbparts += 1;
            }
            if nbparts == 0 {
                self.my_status = BRepFeatStatus::InvalidPlacement;
                return;
            }
            //
            if nbparts >= 2 {
                // preserve everything above the First (OCCT L339).
                for its in &parts {
                    let barycentre = baryc(its, self.base.my_shape.as_ref());
                    let parbar = crate::geomalgo::int_patch::elclib::line_parameter_of_axis(
                        self.my_axis.location,
                        self.my_axis.direction,
                        barycentre,
                    );
                    if parbar > pnt_info_first.parameter() {
                        self.base.keep_part(its);
                    }
                }
            }
        }
    }

    /// OCCT BRepFeat_MakeCylindricalHole::Perform(const double Radius, const
    /// double PFrom, const double PTo, const bool Cont) (L356-470) —
    /// performs every hole of radius Radius located between PFrom and PTo on
    /// the given axis (the second Perform overload; Rust has no overloading).
    pub fn perform_between(&mut self, radius: f64, p_from: f64, p_to: f64, cont: bool) {
        //
        // OCCT L362-366.
        let a_object = self.object_arg();
        //
        self.my_is_blind = false; // OCCT L368
        self.my_validate = cont; // OCCT L369
        self.my_status = BRepFeatStatus::NoError; // OCCT L370
        //
        // OCCT L372-377.
        let the_asi = LocOpeCurveShapeIntersector::new_axis(&self.my_axis, &a_object);
        if !the_asi.is_done() {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L379-389: thePFrom/thePTo ordering.
        let (the_p_from, the_p_to) = if p_from < p_to {
            (p_from, p_to)
        } else {
            (p_to, p_from)
        };
        //
        // OCCT L392-419: the localization.
        let mut pnt_info_first = LocOpePntFace::new();
        let mut pnt_info_last = LocOpePntFace::new();
        let mut ind_from: i32 = 0;
        let mut ind_to: i32 = 0;
        let mut the_or = topods::Orientation::Forward;
        let mut ok = the_asi.localize_after(the_p_from, &mut the_or, &mut ind_from, &mut ind_to);
        if ok {
            if the_or == topods::Orientation::Reversed {
                // reset (OCCT L400)
                ok = the_asi.localize_before_index(ind_from, &mut the_or, &mut ind_from, &mut ind_to);
                // It is possible to find the next.
            }
            if ok && the_or == topods::Orientation::Forward {
                pnt_info_first = the_asi.point(ind_from);
                // OCCT L406: LocalizeBefore(thePTo, ...) — the parameter
                // overload.
                ok = the_asi.localize_before(the_p_to, &mut the_or, &mut ind_from, &mut ind_to);
                if ok {
                    if the_or == topods::Orientation::Forward {
                        // OCCT L411: LocalizeAfter(IndTo, ...) — the index
                        // overload.
                        ok = the_asi.localize_after_index(ind_to, &mut the_or, &mut ind_from, &mut ind_to);
                    }
                    if ok && the_or == topods::Orientation::Reversed {
                        pnt_info_last = the_asi.point(ind_to);
                    }
                }
            }
        }
        //
        if !ok {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L427-433.
        let the_cyl = create_cyl(&pnt_info_first, &pnt_info_last, radius, &self.my_axis);
        self.my_top_face = the_cyl.cyl_top_f.clone();
        self.my_bot_face = the_cyl.cyl_bott_f.clone();
        let the_tool: Option<Shape> = the_cyl.cyl.clone();
        //
        // OCCT L435-438.
        let fuse = 0i32;
        if let Some(the_tool) = the_tool {
            self.base.add_tool(the_tool);
            self.base.set_operation(fuse);
            self.base.perform_bop();
            // OCCT L439-440.
            let parts = self.base.parts_of_tool();
            //
            let mut nbparts = 0i32;
            for _its in &parts {
                nbparts += 1;
            }
            if nbparts == 0 {
                self.my_status = BRepFeatStatus::InvalidPlacement;
                return;
            }
            //
            if nbparts >= 2 {
                // preserve parts between First and Last (OCCT L455).
                for its in &parts {
                    let barycentre = baryc(its, self.base.my_shape.as_ref());
                    let parbar = crate::geomalgo::int_patch::elclib::line_parameter_of_axis(
                        self.my_axis.location,
                        self.my_axis.direction,
                        barycentre,
                    );
                    if parbar >= pnt_info_first.parameter()
                        && parbar <= pnt_info_last.parameter()
                    {
                        self.base.keep_part(its);
                    }
                }
            }
        }
    }

    /// OCCT BRepFeat_MakeCylindricalHole::PerformBlind(const double Radius,
    /// const double Length, const bool Cont) (L474-626) — performs a blind
    /// hole of radius Radius and length Length.
    pub fn perform_blind(&mut self, radius: f64, length: f64, cont: bool) {
        //
        // OCCT L479-483: aObject.IsNull() || !myAxDef || Length <= 0.
        let a_object = self.object_arg();
        if length <= 0. {
            panic!("Standard_ConstructionError"); // OCCT throw (difference #5)
        }
        //
        self.my_is_blind = true; // OCCT L485
        self.my_validate = cont; // OCCT L486
        self.my_status = BRepFeatStatus::NoError; // OCCT L487
        //
        // OCCT L489-494.
        let the_asi = LocOpeCurveShapeIntersector::new_axis(&self.my_axis, &a_object);
        if !the_asi.is_done() {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L496-514: the localization of the hole origin.
        let first: f64;
        let mut ind_from: i32 = 0;
        let mut ind_to: i32 = 0;
        let mut the_or = topods::Orientation::Forward;
        let mut ok = the_asi.localize_after(0., &mut the_or, &mut ind_from, &mut ind_to);
        //
        if ok {
            if the_or == topods::Orientation::Reversed {
                // reset (OCCT L505)
                ok = the_asi.localize_before_index(ind_from, &mut the_or, &mut ind_from, &mut ind_to);
                // it is possible to find the next
            }
            ok = ok && the_or == topods::Orientation::Forward;
        }
        if !ok {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // check a priori the length of the hole (OCCT L516-528).
        let mut if_next: i32 = 0;
        let mut it_next: i32 = 0;
        // OCCT L518: LocalizeAfter(IndTo, ...) — the index overload.
        ok = the_asi.localize_after_index(ind_to, &mut the_or, &mut if_next, &mut it_next);
        if !ok {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        if the_asi.point(if_next).parameter() <= length {
            self.my_status = BRepFeatStatus::HoleTooLong;
            return;
        }
        //
        // OCCT L530-536: NCollection_List<TopoDS_Shape> theList; the version
        // for advanced control fills the faces from IndFrom to ITNext. (The
        // list stays unused: the myBuilder.Perform(theTool, theList, Fuse)
        // call is commented out at OCCT L566-567.)
        let mut the_list: Vec<Shape> = Vec::new();
        for i in ind_from..=it_next {
            if let Some(f) = the_asi.point(i).face() {
                the_list.push(f.clone());
            }
        }
        //
        // OCCT L538: First = theASI.Point(IndFrom).Parameter();
        first = the_asi.point(ind_from).parameter();
        //
        // It is not possible to use infinite cylinder for topological
        // operations. (OCCT L540-547)
        let (p_min, _p_max) = box_parameters(&a_object, &self.my_axis);
        if p_min > length {
            self.my_status = BRepFeatStatus::InvalidPlacement;
            return;
        }
        //
        // OCCT L549-554: Heigth = 3.*(Length - PMin)/2.; theOrig = ...;
        // BRepPrim_Cylinder theCylinder(a5_ao1, Radius, Heigth);
        let heigth = 3. * (length - p_min) / 2.;
        let the_orig = self.my_axis.location + ((3. * p_min - length) / 2.) * self.my_axis.direction;
        let the_cylinder = BRepPrimCylinder::new(the_orig, self.my_axis.direction, radius, heigth);
        //
        // OCCT L556-559: the tool solid.
        // Architecture difference #2: the tool solid stays null until TKPrim
        // produces the cylinder shell.
        let the_tool: Option<Shape> = the_cylinder.shell();
        //
        self.my_top_face = the_cylinder.top_face(); // OCCT L561
        self.my_bot_face = None; // OCCT L562: myBotFace.Nullify();
        //
        // OCCT L565-570.
        let fuse = 0i32;
        if let Some(the_tool) = the_tool {
            self.base.add_tool(the_tool);
            self.base.set_operation(fuse);
            self.base.perform_bop();
            // OCCT L571-572.
            let parts = self.base.parts_of_tool();
            //
            let mut nbparts = 0i32;
            for _its in &parts {
                nbparts += 1;
            }
            if nbparts == 0 {
                self.my_status = BRepFeatStatus::InvalidPlacement;
                return;
            }
            //
            if nbparts >= 2 {
                // preserve the smallest as parameter along the axis (OCCT
                // L587).
                let mut tokeep: Option<Shape> = None;
                let mut parmin = f64::MAX; // OCCT: RealLast()
                for its in &parts {
                    let barycentre = baryc(its, self.base.my_shape.as_ref());
                    let parbar = crate::geomalgo::int_patch::elclib::line_parameter_of_axis(
                        self.my_axis.location,
                        self.my_axis.direction,
                        barycentre,
                    );
                    if parbar >= first && parbar <= parmin {
                        parmin = parbar;
                        tokeep = Some(its.clone());
                    }
                }
                //
                if tokeep.is_none() {
                    // preserve the closest interval (OCCT L603).
                    let mut dmin = f64::MAX; // OCCT: RealLast()
                    for its in &parts {
                        let barycentre = baryc(its, self.base.my_shape.as_ref());
                        let parbar = crate::geomalgo::int_patch::elclib::line_parameter_of_axis(
                            self.my_axis.location,
                            self.my_axis.direction,
                            barycentre,
                        );
                        if (first - parbar).abs() < dmin {
                            dmin = (first - parbar).abs();
                            tokeep = Some(its.clone());
                        }
                    }
                }
                for its in &parts {
                    let is_same = match &tokeep {
                        Some(tk) => shape_is_same(tk, its),
                        None => false,
                    };
                    if is_same {
                        self.base.keep_part(its);
                        break;
                    }
                }
            }
        }
        let _ = &mut the_list;
    }

    /// OCCT BRepFeat_MakeCylindricalHole::Build() (L630-648) — builds the
    /// resulting shape; invalidates the given parts of tools if any, and
    /// performs the result of the local operation.
    pub fn build(&mut self) {
        if self.my_status == BRepFeatStatus::NoError {
            // OCCT L634: PerformResult();
            self.base.perform_result();
            if !self.base.has_errors() {
                // OCCT L637: myStatus = (myValidate) ? Validate() :
                // BRepFeat_NoError;
                self.my_status = if self.my_validate {
                    self.validate()
                } else {
                    BRepFeatStatus::NoError
                };
                if self.my_status == BRepFeatStatus::NoError {
                    // OCCT L640: myShape = Shape();
                    self.my_result = self.base.shape_brep();
                }
            } else {
                // OCCT L645: myStatus = BRepFeat_InvalidPlacement; // why not
                self.my_status = BRepFeatStatus::InvalidPlacement; // why not
            }
        }
    }

    /// OCCT BRepFeat_MakeCylindricalHole::Validate() (L652-688).
    fn validate(&self) -> BRepFeatStatus {
        let mut the_stat = BRepFeatStatus::NoError;
        // OCCT L655: TopExp_Explorer ex(Shape(), TopAbs_FACE);
        let faces: Vec<Shape> = self
            .base
            .my_shape
            .as_ref()
            .map(pool_faces)
            .unwrap_or_default();
        if self.my_is_blind {
            // limit of the hole (OCCT L657)
            let mut found = false;
            for f in &faces {
                if face_is_same(f, &self.my_top_face) {
                    found = true;
                    break;
                }
            }
            // OCCT L665-668: if (!ex.More()) thestat = BRepFeat_HoleTooLong;
            if !found {
                the_stat = BRepFeatStatus::HoleTooLong;
            }
        } else {
            // OCCT L672-678.
            for f in &faces {
                if face_is_same(f, &self.my_top_face) {
                    return BRepFeatStatus::InvalidPlacement;
                }
            }
            // OCCT L679: ex.ReInit() — re-explore from the beginning.
            for f in &faces {
                if face_is_same(f, &self.my_bot_face) {
                    return BRepFeatStatus::InvalidPlacement;
                }
            }
        }
        the_stat
    }
}

/// OCCT TopoDS_Shape::IsSame — same TShape and same Location (orientation
/// ignored); the null face never matches.
fn face_is_same(f: &Shape, other: &Option<Shape>) -> bool {
    match other {
        Some(o) => f.ptr_id() == o.ptr_id() && f.location == o.location,
        None => false,
    }
}

/// OCCT TopoDS_Shape::IsSame between two non-null shapes.
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

// OCCT BRepFeat_MakeCylindricalHole::Perform / Status / Init — see the impl
// block above. The Perform out-argument style of OCCT (PartsOfTool(aLT),
// GetOffset(..., outOff), CreateCyl(..., Cyl, CylTopF, CylBottF)) is mapped
// to return values (see the header, difference #6).
