// OCCT BRepFeat_MakeRevolutionForm.hxx L17-125 + BRepFeat_MakeRevolutionForm.cxx
// L17-2030 + BRepFeat_MakeRevolutionForm.lxx L17-55 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeRevolutionForm.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeRevolutionForm.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeRevolutionForm.lxx
//
// OCCT inheritance chain (BRepFeat_MakeRevolutionForm.hxx L52):
//   BRepFeat_MakeRevolutionForm : BRepFeat_RibSlot : BRepBuilderAPI_MakeShape
// NOTE (OCCT 8.0): RibSlot does NOT derive from BRepFeat_Form, so this
// subclass has no Curves/BarycCurve overrides and never enters
// BRepFeat_Form::GlobalPerform — the topological reconstruction is
// BRepFeat_RibSlot::LFPerform (rib_slot.lf_perform()). Rust composition: the
// RibSlot base sub-object is the `rib_slot` field; the OCCT protected
// members are pub(crate) on BRepFeatRibSlot and are driven exactly as the
// C++ subclass does.
//
// Architecture differences (referenced from the affected functions):
// 1. NCollection_DataMap (mySlface) maps to HashMap keyed by (TShape ptr,
//    Location); the key shape is carried as the tuple head.
// 2. BRepAlgoAPI_Cut / Common are the CutVehicle of brep_feat_form_2
//    (with_operation); BRepAlgoAPI_Section is the SectionOp of
//    bop/brep_algo_api (the true BOPAlgo_Section vehicle).
// 3. Geom2dAPI_ExtremaCurveCurve (Init, cxx L154-158) is not translated
//    (TKGeomBase/Extrema_ExtCC on 2d curves); the carrier below reports
//    NbExtrema() = 0, so the OCCT "if (NbExtrema() >= 1)" guard keeps
//    L = RealLast — the same guard structure driven by the missing
//    dependency (the MakeCylindricalHole arch.-difference #1 model).
// 4. BRepExtrema_ExtCF (Init, cxx L435-447) is not translated; the carrier
//    reports NbExt() = 0, so the OCCT guard drives Sliding = false into the
//    NoSlidingProfile branch.
// 5. BRepPrimAPI_MakeBox (Init, cxx L268-269) is not translated (TKPrim —
//    the same gap as the BRepPrim_Cylinder carrier of
//    brep_feat_make_cylindrical_hole.rs); the carrier returns the null
//    solid.
// 6. BRepBuilderAPI_Transform (Perform, cxx L1204-1206) — the myAngle2 != 0
//    rotation branch; the shape-transform engine is pending (the same gap
//    as the BRepTools_Modifier carrier of loc_ope_prism.rs) and stops at
//    the GAP panic at its spot.
// 7. gp_Circ / gp_Ax2 / gp_Pnt::Rotated / gp_Pln::Translated / Geom_Plane::D0
//    are small analytic carriers below over the rcad geom types.
// 8. BRepTools_WireExplorer maps to the wire edge list (the same reduction
//    as the wire_explorer_edges carrier of brep_feat_rib_slot_b.rs).
// 9. BRep_Tool::IsClosed / TopExp::FirstVertex/LastVertex / IntPar /
//    Normal / CheckPoint / ExtremeFaces / SlidingProfile / NoSlidingProfile
//    / HeightMax are the RibSlot re-hosts consumed through the rib_slot
//    sub-object (brep_feat_rib_slot.rs / brep_feat_rib_slot_b.rs).

use crate::bop::algo::builder::BooleanOpType;
use crate::bop::brep_algo_api::SectionOp;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::brep_feat_form_2::{
    brep_algo_is_valid, brep_feat_face_until, brep_feat_tool, CutVehicle,
};
use crate::feat::brep_feat_rib_slot::{
    brep_tool_curve, brep_tool_is_closed, brep_tool_pnt, brep_tool_surface,
    brep_tool_tolerance, geom_api_to_2d, make_edge_c_v_v, make_vertex,
    top_exp_first_vertex, top_exp_last_vertex, BRepFeatRibSlot, Geom2dAPIInterCurveCurve,
    GeomAPIProjectPointOnCurve,
};
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_cs_intersector::LocOpeCSIntersector;
use crate::feat::loc_ope_revolution_form::LocOpeRevolutionForm;
use glam::{DVec2, DVec3};
use indexmap::IndexMap;
use rcad_kernel::base::extrema::ExtPC;
use rcad_kernel::geom::{
    Circle3, Curve2d, Curve3, Line3, Plane, Surface3, TrimmedCurve3, TrimmedSurface,
};
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo::topods::{BRep, BRepBuilder, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, State};
use rcad_kernel::CurveEval;
use std::collections::HashMap;
use std::f64::consts::PI;

/// OCCT BRepTools_WireExplorer — the wire edge list (architecture
/// difference #8: the same reduction as the wire_explorer_edges carrier of
/// brep_feat_rib_slot_b.rs).
fn wire_edges(the_wire: &Shape) -> Vec<Shape> {
    match the_wire.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// The stable X direction of a normal (the gp_Ax3 construction rule of
/// gp_Ax3.cxx L29-80 reduced to the rcad Plane::new carrier).
fn stable_x_dir(normal: DVec3) -> DVec3 {
    let n = normal;
    let (a, b, c) = (n.x.abs(), n.y.abs(), n.z.abs());
    let d = if b <= a && b <= c {
        if a > c {
            DVec3::new(-c, 0.0, a)
        } else {
            DVec3::new(c, 0.0, -a)
        }
    } else if a <= b && a <= c {
        if b > c {
            DVec3::new(0.0, -c, b)
        } else {
            DVec3::new(0.0, c, -b)
        }
    } else if a > b {
        DVec3::new(-b, a, 0.0)
    } else {
        DVec3::new(b, -a, 0.0)
    };
    d.normalize_or_zero()
}

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

/// OCCT gp_Pnt::Rotated(A1, Ang) — the rotation of a point about the axis
/// (the Rodrigues form of gp_Pnt::Rotated/gp_Trsf::SetRotation).
fn gp_pnt_rotated(the_p: DVec3, the_axe: &Ax1, the_ang: f64) -> DVec3 {
    let d = the_axe.direction;
    let (sin_a, cos_a) = the_ang.sin_cos();
    let v = the_p - the_axe.location;
    let cross = d.cross(v);
    let dot = d.dot(v);
    the_axe.location + v * cos_a + cross * sin_a + d * (dot * (1.0 - cos_a))
}

/// OCCT gp_Ax1::IsCoaxial(theOther, theAngTol, theLinTol) — parallel
/// directions and the other location on the axis line (gp_Ax1.hxx
/// IsCoaxial).
fn gp_ax1_is_coaxial(the_axe: &Ax1, the_other: &Ax1, ang_tol: f64, lin_tol: f64) -> bool {
    let d1 = the_axe.direction;
    let d2 = the_other.direction;
    let an_ang = DVec3::angle_between(d1, d2);
    let parallel = an_ang <= ang_tol || std::f64::consts::PI - an_ang <= ang_tol;
    let v = the_other.location - the_axe.location;
    let cross = d1.cross(v).length();
    parallel && cross <= lin_tol
}

/// OCCT Geom_Plane::D0(U, V, P) — the point of the plane at (u, v).
fn plane_d0(the_pln: &Plane, u: f64, v: f64) -> DVec3 {
    the_pln.origin + the_pln.u_dir * u + the_pln.v_dir * v
}

/// OCCT Geom2dAPI_ExtremaCurveCurve carrier (architecture difference #3) —
/// NbExtrema() = 0 keeps the OCCT guard structure.
#[allow(dead_code)]
struct Geom2dAPIExtremaCurveCurve;

impl Geom2dAPIExtremaCurveCurve {
    #[allow(dead_code)]
    fn new(_c1: &Curve2d, _c2: &Curve2d, _f1: f64, _l1: f64, _f2: f64, _l2: f64) -> Self {
        Geom2dAPIExtremaCurveCurve
    }
    #[allow(dead_code)]
    fn nb_extrema(&self) -> i32 {
        0
    }
    #[allow(dead_code)]
    fn lower_distance(&self) -> f64 {
        f64::MAX
    }
}

/// OCCT BRepExtrema_ExtCF carrier (architecture difference #4) — NbExt() = 0
/// drives the OCCT Sliding = false guard.
struct BRepExtremaExtCF;

impl BRepExtremaExtCF {
    fn new(_e: &Shape, _f: &Shape) -> Self {
        BRepExtremaExtCF
    }
    fn nb_ext(&self) -> i32 {
        0
    }
    #[allow(dead_code)]
    fn square_distance(&self, _n: i32) -> f64 {
        f64::MAX
    }
}

/// OCCT BRepPrimAPI_MakeBox(P1, P2).Solid() carrier (architecture
/// difference #5) — the null solid until the TKPrim stage lands.
fn brep_prim_make_box(_p1: DVec3, _p2: DVec3) -> Shape {
    Shape::null()
}

/// OCCT BRepLib_MakeFace(Pln, U1, U2, V1, V2) — the face on the plane with
/// the uv bounds (the rectangular trimmed surface carrier of
/// brep_feat_face_until).
fn make_face_pln_bounds(
    the_b: &mut BRepBuilder,
    the_pool: &mut BRep,
    the_pln: &Plane,
    u1: f64,
    u2: f64,
    v1: f64,
    v2: f64,
) -> Shape {
    let str_surf = TrimmedSurface::new(Surface3::Plane(*the_pln), u1, u2, v1, v2);
    the_b.make_face(the_pool, Some(Surface3::Trimmed(str_surf)), Shape::null())
}

/// OCCT BRepLib_MakeFace(Pln, wire, Inside) — the face on the plane limited
/// by the wire.
fn make_face_plane_wire(
    the_b: &mut BRepBuilder,
    the_pool: &mut BRep,
    the_pln: &Plane,
    the_wire: &Shape,
) -> Shape {
    the_b.make_face(the_pool, Some(Surface3::Plane(*the_pln)), the_wire.clone())
}

/// OCCT BRepBuilderAPI_Transform(T).Perform(S, false) carrier (architecture
/// difference #6).
fn brep_builder_api_transform(the_shape: &Shape, _the_trsf_rot: (DVec3, DVec3, f64)) -> Shape {
    let _ = the_shape;
    panic!("GAP(BRepFeat_MakeRevolutionForm): BRepBuilderAPI_Transform (the shape rotation engine is pending translation)");
}

/// OCCT BRepFeat_MakeRevolutionForm — describes functions to build
/// revolved form features (BRepFeat_MakeRevolutionForm.hxx L36-51).
pub struct BRepFeatMakeRevolutionForm {
    /// The BRepFeat_RibSlot base sub-object (the composition carrier).
    pub rib_slot: BRepFeatRibSlot,
    my_axe: Ax1, // OCCT: myAxe (gp_Ax1)
    my_height1: f64, // OCCT: myHeight1
    my_height2: f64, // OCCT: myHeight2
    #[allow(dead_code)]
    my_sliding: bool, // OCCT: mySliding
    my_pln: Plane,   // OCCT: myPln (Handle(Geom_Plane))
    #[allow(dead_code)]
    my_bnd: f64, // OCCT: myBnd
    // OCCT: mySlface (architecture difference #1).
    my_slface: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    my_list_of_edges: Vec<Shape>, // OCCT: myListOfEdges
    #[allow(dead_code)]
    my_tol: f64, // OCCT: myTol
    #[allow(dead_code)]
    my_angle1: f64, // OCCT: myAngle1
    #[allow(dead_code)]
    my_angle2: f64, // OCCT: myAngle2
}

impl BRepFeatMakeRevolutionForm {
    /// OCCT BRepFeat_MakeRevolutionForm::BRepFeat_MakeRevolutionForm()
    /// (lxx: the empty constructor).
    pub fn new() -> Self {
        BRepFeatMakeRevolutionForm {
            rib_slot: BRepFeatRibSlot::new(),
            my_axe: Ax1::new(DVec3::ZERO, DVec3::Z),
            my_height1: 0.0,
            my_height2: 0.0,
            my_sliding: false,
            my_pln: Plane::new(DVec3::ZERO, DVec3::Z),
            my_bnd: 0.0,
            my_slface: HashMap::new(),
            my_list_of_edges: Vec::new(),
            my_tol: 0.0,
            my_angle1: 0.0,
            my_angle2: 0.0,
        }
    }

    /// OCCT BRepFeat_MakeRevolutionForm::BRepFeat_MakeRevolutionForm(Sbase,
    /// W, Plane, Axis, Height1, Height2, Fuse, Sliding) (hxx L67-74) — Rust
    /// has no overloading: the `with_init` suffix.
    #[allow(clippy::too_many_arguments)]
    pub fn with_init(
        sbase: &Shape,
        w: &Shape,
        plane: &Plane,
        axis: Ax1,
        height1: f64,
        height2: f64,
        fuse: i32,
        sliding: &mut bool,
    ) -> Self {
        let mut res = BRepFeatMakeRevolutionForm::new();
        res.init(sbase, w, plane, axis, height1, height2, fuse, sliding);
        res
    }

    /// OCCT BRepFeat_MakeRevolutionForm::Init (cxx L101-1121).
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn init(
        &mut self,
        sbase: &Shape,
        w: &Shape,
        plane: &Plane,
        axis: Ax1,
        h1: f64,
        h2: f64,
        mode: i32,
        modify: &mut bool,
    ) {
        // OCCT L115-116.
        let revol_rib = true;
        self.rib_slot.done();

        // modify = 0 if it is not required to make sliding
        //        = 1 if it is intended to try to make sliding
        // OCCT L120.
        let mut sliding = *modify;

        // OCCT L122-123.
        self.my_axe = axis;
        let line = Curve3::Line(Line3::new(self.my_axe.location, self.my_axe.direction));

        // OCCT L126-140.
        let mut a_si = LocOpeCSIntersector::with_shape(sbase);
        let scur = vec![Some(line.clone())];
        a_si.perform_cur(&scur);
        let line_first: f64;
        let line_last: f64;
        if a_si.is_done() && a_si.nb_points(1) >= 2 {
            line_last = a_si.point(1, a_si.nb_points(1)).parameter();
            line_first = a_si.point(1, 1).parameter();
        } else {
            line_first = f64::MIN; // OCCT: RealFirst()
            line_last = f64::MAX; // OCCT: RealLast()
        }

        // OCCT L142.
        let Some(ln2d) = geom_api_to_2d(&line, plane) else {
            panic!("GeomAPI::To2d failed (no OCCT counterpart)");
        };

        // OCCT L144-185.
        let mut rad = f64::MAX; // OCCT: RealLast()
        for exx in explorer(w, ShapeType::Edge, ShapeType::Shape) {
            let e = exx;
            // OCCT L152.
            let Some((c, f, l)) = brep_tool_curve(&e) else {
                continue;
            };
            // OCCT L153.
            let Some(c2d) = geom_api_to_2d(&c, plane) else {
                continue;
            };
            // OCCT L154-159 (architecture difference #3).
            let extr = Geom2dAPIExtremaCurveCurve::new(&ln2d, &c2d, line_first, line_last, f, l);
            let mut big_l = f64::MAX; // OCCT: RealLast()
            if extr.nb_extrema() >= 1 {
                big_l = extr.lower_distance();
            }
            // OCCT L160-163.
            let p1 = c.point_at(f);
            let p2 = c.point_at(l);
            let proj1 = GeomAPIProjectPointOnCurve::new(p1, &line);
            let proj2 = GeomAPIProjectPointOnCurve::new(p2, &line);
            // OCCT L164-173.
            if proj1.nb_points() < 1 || proj2.nb_points() < 1 {
                self.rib_slot.my_status_error = BRepFeatStatusError::NoProjPt;
                self.rib_slot.not_done();
                return;
            }
            // OCCT L174-180.
            let par1 = proj1.distance(1);
            let par2 = proj2.distance(1);
            let par = par1.min(par2);
            if par < big_l {
                big_l = par;
            }
            // OCCT L181-184.
            if big_l < rad && big_l > 0.0 {
                rad = big_l;
            }
        }

        // OCCT L187-204.
        let height = h1.min(h2);
        if rad <= height {
            rad = height + 0.01 * height;
        }
        self.my_angle1 = (h1 / rad).asin() + PI / 10.0;
        self.my_angle2 = (h2 / rad).asin() + PI / 10.0;
        if (self.my_angle1 - PI / 2.0) > CONFUSION {
            self.my_angle1 = PI / 2.0;
        }
        if (self.my_angle2 - PI / 2.0) > CONFUSION {
            self.my_angle2 = PI / 2.0;
        }

        // OCCT L206-209.
        self.rib_slot.my_skface = Shape::null();
        self.rib_slot.my_pbase = Shape::null();
        self.rib_slot.my_fuse = mode != 0;

        // ---Determination Tolerance : tolerance max on parameters

        // OCCT L221.
        self.my_tol = CONFUSION;
        // OCCT L223-231.
        for exx in explorer(w, ShapeType::Vertex, ShapeType::Shape) {
            let tol = brep_tool_tolerance(&exx);
            if tol > self.my_tol {
                self.my_tol = tol;
            }
        }
        // OCCT L233-241.
        for exx in explorer(sbase, ShapeType::Vertex, ShapeType::Shape) {
            let tol = brep_tool_tolerance(&exx);
            if tol > self.my_tol {
                self.my_tol = tol;
            }
        }

        // OCCT L243-248.
        self.rib_slot.my_wire = w.clone();
        self.my_pln = *plane;
        self.my_height1 = h1;
        self.my_height2 = h2;

        // OCCT L250-255.
        self.rib_slot.my_sbase = sbase.clone();
        self.my_slface.clear();
        self.rib_slot.my_shape = None;
        self.rib_slot.my_map.clear();
        self.rib_slot.my_f_shape = Shape::null();
        self.rib_slot.my_l_shape = Shape::null();

        // ---Calculate bounding box

        // OCCT L260-266.
        let mut the_list: Vec<Shape> = Vec::new();
        let mut u = Shape::null(); // OCCT: U.Nullify()
        let mut first_corner = DVec3::ZERO;
        let mut last_corner = DVec3::ZERO;
        let bnd = self.rib_slot.height_max(&self.rib_slot.my_sbase, &mut u, &mut first_corner, &mut last_corner);
        self.my_bnd = bnd;

        // OCCT L268-269 (architecture difference #5).
        let bnd_box = brep_prim_make_box(first_corner, last_corner);

        // ---Construction of the working plane face (section bounding box)

        // OCCT L272-273.
        let mut pool = BRep::new();
        let mut bb = BRepBuilder::new();
        let plane_face = make_face_pln_bounds(
            &mut bb,
            &mut pool,
            &self.my_pln,
            -6.0 * self.my_bnd,
            6.0 * self.my_bnd,
            -6.0 * self.my_bnd,
            6.0 * self.my_bnd,
        );

        // OCCT L275-281.
        let plane_s =
            CutVehicle::with_operation(&bnd_box, &plane_face, BooleanOpType::Intersection);
        let plane_sect = plane_s.shape().cloned().unwrap_or_else(Shape::null);
        let www = explorer(&plane_sect, ShapeType::Wire, ShapeType::Shape)
            .into_iter()
            .next()
            .unwrap_or_else(Shape::null);
        let mut bndface_pool = BRep::new();
        let mut bndface_b = BRepBuilder::new();
        let bnd_face = make_face_plane_wire(&mut bndface_b, &mut bndface_pool, &self.my_pln, &www);

        // ---Find base faces of the rib

        // OCCT L284-294.
        let mut first_edge = Shape::null();
        let mut last_edge = Shape::null();
        let mut first_face = Shape::null();
        let mut last_face = Shape::null();
        let mut first_vertex = Shape::null();
        let mut last_vertex = Shape::null();
        let mut on_first_face = false;
        let mut on_last_face = false;
        let mut pt_on_first_edge = false;
        let mut pt_on_last_edge = false;
        let mut on_first_edge = Shape::null();
        let mut on_last_edge = Shape::null();

        // OCCT L296-321.
        let data = self.rib_slot.extreme_faces(
            revol_rib,
            self.my_bnd,
            &self.my_pln,
            &mut first_edge,
            &mut last_edge,
            &mut first_face,
            &mut last_face,
            &mut first_vertex,
            &mut last_vertex,
            &mut on_first_face,
            &mut on_last_face,
            &mut pt_on_first_edge,
            &mut pt_on_last_edge,
            &mut on_first_edge,
            &mut on_last_edge,
        );
        if !data {
            self.set_derived_status(BRepFeatStatusError::NoExtFace);
            self.rib_slot.not_done();
            return;
        }

        // ---Proofing Point for the side of the wire to be filled - material
        // side

        // OCCT L324.
        let check_pnt = self.rib_slot.check_point(&first_edge, bnd / 10.0, &self.my_pln);

        // ---Control sliding valid
        // Many cases when the sliding is abandoned

        // OCCT L330.
        let mut concavite = 3i32; // a priori the profile is not concave

        // OCCT L332-333.
        self.rib_slot.my_first_pnt = brep_tool_pnt(&first_vertex);
        self.rib_slot.my_last_pnt = brep_tool_pnt(&last_vertex);

        // SliList : list of faces concerned by the rib

        // OCCT L336-337.
        let mut sli_list: Vec<Shape> = Vec::new();
        sli_list.push(first_face.clone());

        // OCCT L339-357.
        if sliding {
            let s = match brep_tool_surface(&first_face) {
                Some(mut s) => {
                    if let Surface3::Trimmed(t) = &s {
                        s = (*t.basis).clone();
                    }
                    s
                }
                None => Surface3::Plane(self.my_pln),
            };
            match &s {
                Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) => {}
                _ => {
                    sliding = false;
                }
            }
        }

        // OCCT L359-373.
        if sliding {
            let ss = match brep_tool_surface(&last_face) {
                Some(mut ss) => {
                    if let Surface3::Trimmed(t) = &ss {
                        ss = (*t.basis).clone();
                    }
                    ss
                }
                None => Surface3::Plane(self.my_pln),
            };
            match &ss {
                Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) => {}
                _ => {
                    sliding = false;
                }
            }
        }

        // OCCT L379-380.
        let mut first_circle: Option<Circle3> = None;
        let mut last_circle: Option<Circle3> = None;

        // OCCT L382-514.
        if sliding {
            // OCCT L384-394.
            let proj = GeomAPIProjectPointOnCurve::new(self.rib_slot.my_first_pnt, &line);
            if proj.nb_points() < 1 {
                self.set_derived_status(BRepFeatStatusError::NoProjPt);
                self.rib_slot.not_done();
                return;
            }
            let first_rayon = proj.distance(1);
            let first_center = line.point_at(
                self.line_project_parameter(&line, self.rib_slot.my_first_pnt),
            );
            // OCCT L398-410.
            let proj1 = GeomAPIProjectPointOnCurve::new(self.rib_slot.my_last_pnt, &line);
            if proj.nb_points() < 1 {
                self.set_derived_status(BRepFeatStatusError::NoProjPt);
                self.rib_slot.not_done();
                return;
            }
            let last_rayon = proj1.distance(1);
            let last_center = line.point_at(
                self.line_project_parameter(&line, self.rib_slot.my_last_pnt),
            );

            // OCCT L412-421: the circles of the sliding check.
            let axv = self.my_axe.direction;
            let the_fc = Circle3 {
                center: first_center,
                normal: axv,
                x_dir: stable_x_dir(axv),
                y_dir: axv.cross(stable_x_dir(axv)),
                radius: first_rayon,
            };
            let the_lc = Circle3 {
                center: last_center,
                normal: axv,
                x_dir: stable_x_dir(axv),
                y_dir: axv.cross(stable_x_dir(axv)),
                radius: last_rayon,
            };
            let r_first_pnt1 = gp_pnt_rotated(self.rib_slot.my_first_pnt, &self.my_axe, self.my_angle1);
            let r_last_pnt1 = gp_pnt_rotated(self.rib_slot.my_last_pnt, &self.my_axe, self.my_angle1);
            let r_first_pnt2 = gp_pnt_rotated(self.rib_slot.my_first_pnt, &self.my_axe, -self.my_angle2);
            let r_last_pnt2 = gp_pnt_rotated(self.rib_slot.my_last_pnt, &self.my_axe, -self.my_angle2);

            // OCCT L423-428.
            let mut bu_pool = BRep::new();
            let v1 = make_vertex(&mut bu_pool, r_first_pnt2);
            let v2 = make_vertex(&mut bu_pool, r_first_pnt1);
            let v3 = make_vertex(&mut bu_pool, r_last_pnt2);
            let v4 = make_vertex(&mut bu_pool, r_last_pnt1);

            // OCCT L430-431: BRepLib_MakeEdge(theFC, v1, v2) — the circle
            // edges (the curve+vertices carrier).
            let mut ee_pool = BRep::new();
            let ee1 = make_edge_c_v_v(
                &mut ee_pool,
                &Curve3::Circle(the_fc),
                &v1,
                &v2,
            );
            let ee2 = make_edge_c_v_v(
                &mut ee_pool,
                &Curve3::Circle(the_lc),
                &v3,
                &v4,
            );

            // OCCT L433-448 (architecture difference #4).
            if sliding && !pt_on_first_edge {
                let ext1 = BRepExtremaExtCF::new(&ee1, &first_face);
                if ext1.nb_ext() < 1 || ext1.square_distance(1) > 1e-14 {
                    // OCCT: Precision::SquareConfusion().
                    sliding = false;
                }
            }
            if sliding && !pt_on_last_edge {
                let ext2 = BRepExtremaExtCF::new(&ee2, &last_face);
                if ext2.nb_ext() < 1 || ext2.square_distance(1) > 1e-14 {
                    sliding = false;
                }
            }
            // OCCT L449-480.
            if sliding && pt_on_first_edge {
                match brep_tool_curve(&on_first_edge) {
                    Some((first_crv, _f, _l)) => {
                        if let Curve3::Circle(cc) = &first_crv {
                            let circ = *cc;
                            first_circle = Some(circ);
                            let circax = Ax1::new(circ.center, circ.normal);
                            if !gp_ax1_is_coaxial(&circax, &self.my_axe, CONFUSION, CONFUSION)
                            {
                                sliding = false;
                            } else if (circ.radius - first_rayon).abs() >= CONFUSION {
                                sliding = false;
                            }
                        } else {
                            sliding = false;
                        }
                    }
                    None => {
                        sliding = false;
                    }
                }
            }
            // OCCT L482-513.
            if sliding && pt_on_last_edge {
                match brep_tool_curve(&on_last_edge) {
                    Some((last_crv, _f, _l)) => {
                        if let Curve3::Circle(cc) = &last_crv {
                            let circ = *cc;
                            last_circle = Some(circ);
                            let rad_c = circ.radius;
                            let circax = Ax1::new(circ.center, circ.normal);
                            if !gp_ax1_is_coaxial(&circax, &self.my_axe, CONFUSION, CONFUSION)
                            {
                                sliding = false;
                            } else if (rad_c - last_rayon).abs() >= CONFUSION {
                                sliding = false;
                            }
                        } else {
                            sliding = false;
                        }
                    }
                    None => {
                        sliding = false;
                    }
                }
            }
        }
        let _ = (&first_circle, &last_circle);

        // ---case of sliding : construction of the face profile

        // OCCT L522-572.
        if sliding {
            let mut prof = Shape::null();
            let profile_ok = self.rib_slot.sliding_profile(
                &mut prof,
                revol_rib,
                self.my_tol,
                &mut concavite,
                &self.my_pln,
                &bnd_face,
                check_pnt,
                &first_face,
                &last_face,
                &first_vertex,
                &last_vertex,
                &first_edge,
                &last_edge,
            );
            if !profile_ok {
                self.set_derived_status(BRepFeatStatusError::NoFaceProf);
                self.rib_slot.not_done();
                return;
            }
            // ---Propagation on faces of the initial shape
            // to find the faces concerned by the rib
            // OCCT L560-571.
            let mut falseside = true;
            sliding = self.propagate(
                &mut sli_list,
                &prof,
                self.rib_slot.my_first_pnt,
                self.rib_slot.my_last_pnt,
                &mut falseside,
            );
            if !falseside {
                self.set_derived_status(BRepFeatStatusError::FalseSide);
                self.rib_slot.not_done();
                return;
            }
        }

        // ---Generation of the base profile of the rib

        // OCCT L576-583.
        let mut bb_pool = BRep::new();
        let mut bb = BRepBuilder::new();
        let w = bb.make_wire(&mut bb_pool);
        let mut the_previous_edge = Shape::null();
        let mut the_fv = Shape::null();

        // counter of the number of edges to fill the map
        // OCCT L583.
        let mut counter = 1i32;

        // ---case of sliding

        // OCCT L586-908 (the sliding wire reconstruction + the SlidMap
        // re-bind).
        if sliding && !self.my_list_of_edges.is_empty() {
            // OCCT L588: BRepTools_WireExplorer EX1(myWire) (architecture
            // difference #8).
            let wire_edges_list = wire_edges(&self.rib_slot.my_wire);
            let mut idx = 0usize;
            while idx < wire_edges_list.len() {
                let e = &wire_edges_list[idx];
                idx += 1;
                // OCCT L592-596.
                if !self.rib_slot.my_lfmap.contains_key(&shape_key(e)) {
                    self.rib_slot
                        .my_lfmap
                        .insert(shape_key(e), (e.clone(), Vec::new()));
                }
                if shape_is_same(e, &first_edge) {
                    // OCCT L599-615.
                    let Some((mut cc, _f, _l)) = brep_tool_curve(e) else {
                        break;
                    };
                    let pt = if !shape_is_same(&first_edge, &last_edge) {
                        brep_tool_pnt(&top_exp_last_vertex(e, true))
                    } else {
                        let pt = self.rib_slot.my_last_pnt;
                        let fpar = BRepFeatRibSlot::int_par(&cc, self.rib_slot.my_first_pnt);
                        let lpar = BRepFeatRibSlot::int_par(&cc, pt);
                        if fpar > lpar {
                            cc = crate::feat::brep_feat_rib_slot::geom_curve_reversed(&cc);
                        }
                        pt
                    };
                    // OCCT L616-633.
                    let ee1 = if the_previous_edge.is_null() {
                        let v1 = make_vertex(&mut bb_pool, self.rib_slot.my_first_pnt);
                        let v2 = make_vertex(&mut bb_pool, pt);
                        make_edge_c_v_v(&mut bb_pool, &cc, &v1, &v2)
                    } else {
                        let v1 = top_exp_last_vertex(&the_previous_edge, true);
                        let v2 = make_vertex(&mut bb_pool, pt);
                        make_edge_c_v_v(&mut bb_pool, &cc, &v1, &v2)
                    };
                    let mut ee1 = ee1;
                    ee1.orientation = e.orientation; // OCCT L632: ee1.Oriented(E.Orientation())
                    // OCCT L635-644.
                    if counter == 1 {
                        the_fv = top_exp_first_vertex(&ee1, true);
                    }
                    self.rib_slot
                        .my_lfmap
                        .entry(shape_key(e))
                        .or_insert_with(|| (e.clone(), Vec::new()))
                        .1
                        .push(ee1.clone());
                    bb.add_to_wire(&mut bb_pool, w.clone(), ee1.clone());
                    the_previous_edge = ee1;
                    counter += 1;
                    idx += 1;
                    break;
                }
            }

            // Case of several edges

            // OCCT L649-722.
            if !shape_is_same(&first_edge, &last_edge) {
                while idx < wire_edges_list.len() {
                    let e = &wire_edges_list[idx];
                    idx += 1;
                    // OCCT L654-658.
                    if !self.rib_slot.my_lfmap.contains_key(&shape_key(e)) {
                        self.rib_slot
                            .my_lfmap
                            .insert(shape_key(e), (e.clone(), Vec::new()));
                    }
                    the_list.push(e.clone());
                    // OCCT L660-689.
                    if !shape_is_same(e, &last_edge) {
                        let Some((ccc, _f, _l)) = brep_tool_curve(e) else {
                            break;
                        };
                        let (v1, v2);
                        if !the_previous_edge.is_null() {
                            v1 = top_exp_last_vertex(&the_previous_edge, true);
                            v2 = top_exp_last_vertex(e, true);
                        } else {
                            v1 = top_exp_first_vertex(e, true);
                            v2 = top_exp_last_vertex(e, true);
                        }
                        let e11 = make_edge_c_v_v(&mut bb_pool, &ccc, &v1, &v2);
                        let mut e11 = e11;
                        e11.orientation = e.orientation; // OCCT L678
                        the_previous_edge = e11.clone();
                        self.rib_slot
                            .my_lfmap
                            .entry(shape_key(e))
                            .or_insert_with(|| (e.clone(), Vec::new()))
                            .1
                            .push(e11.clone());
                        bb.add_to_wire(&mut bb_pool, w.clone(), e11.clone());
                        if counter == 1 {
                            the_fv = top_exp_first_vertex(&e11, true);
                        }
                        counter += 1;
                    } else {
                        // OCCT L690-720.
                        let Some((cc, _f, _l)) = brep_tool_curve(e) else {
                            break;
                        };
                        let pf = brep_tool_pnt(&top_exp_first_vertex(e, true));
                        let pl = self.rib_slot.my_last_pnt;
                        let ee = if the_previous_edge.is_null() {
                            let v1 = make_vertex(&mut bb_pool, pf);
                            let v2 = make_vertex(&mut bb_pool, pl);
                            make_edge_c_v_v(&mut bb_pool, &cc, &v1, &v2)
                        } else {
                            let v1 = top_exp_last_vertex(&the_previous_edge, true);
                            let v2 = make_vertex(&mut bb_pool, pl);
                            make_edge_c_v_v(&mut bb_pool, &cc, &v1, &v2)
                        };
                        let mut ee = ee;
                        ee.orientation = e.orientation; // OCCT L708
                        bb.add_to_wire(&mut bb_pool, w.clone(), ee.clone());
                        self.rib_slot
                            .my_lfmap
                            .entry(shape_key(e))
                            .or_insert_with(|| (e.clone(), Vec::new()))
                            .1
                            .push(ee.clone());
                        if counter == 1 {
                            the_fv = top_exp_first_vertex(&ee, true);
                        }
                        the_previous_edge = ee;
                        counter += 1;
                        break;
                    }
                }
            }

            // OCCT L724-853: the while(!FirstOK) pass.
            let mut it_idx = 0usize;
            let mut first_ok = false;
            let mut last_ok = false;
            let mut the_last_pnt = self.rib_slot.my_last_pnt;
            let mut sens = 0i32;
            let mut the_edge = Shape::null();
            let mut the_l_edge = Shape::null();
            let mut counter1 = counter;
            let mut new_list_of_edges: Vec<Shape> = Vec::new();
            while !first_ok {
                if it_idx >= self.my_list_of_edges.len() {
                    // OCCT L843-851: the it.More() exhaustion -> Sliding =
                    // false.
                    sliding = false;
                    break;
                }
                let edg = self.my_list_of_edges[it_idx].clone();
                // OCCT L736-747.
                let Some((ccc, f, l)) = brep_tool_curve(&edg) else {
                    break;
                };
                let mut cc = Curve3::Trimmed(TrimmedCurve3::new(ccc, f, l));
                if edg.orientation == rcad_kernel::topods::Orientation::Reversed {
                    // OCCT L743: cc->Reverse().
                    cc = Curve3::Trimmed(TrimmedCurve3::new(
                        crate::feat::brep_feat_rib_slot::geom_curve_reversed(
                            match &cc {
                                Curve3::Trimmed(t) => &t.curve,
                                _ => unreachable!(),
                            },
                        ),
                        f,
                        l,
                    ));
                }
                let dom = cc.default_domain();
                let mut fp = cc.point_at(dom[0]);
                let mut lp = cc.point_at(dom[1]);
                let mut dist = fp.distance(the_last_pnt);
                // OCCT L748-763.
                if dist <= self.my_tol {
                    sens = 1;
                    last_ok = true;
                } else {
                    dist = lp.distance(the_last_pnt);
                    if dist <= self.my_tol {
                        sens = 2;
                        last_ok = true;
                        // OCCT L761: cc->Reverse().
                        let rev = crate::feat::brep_feat_rib_slot::geom_curve_reversed(&cc);
                        cc = rev;
                        let rdom = cc.default_domain();
                        fp = cc.point_at(rdom[0]);
                        lp = cc.point_at(rdom[1]);
                    }
                }
                // OCCT L764-774.
                let mut first_flag = 0i32;
                if sens == 1 && lp.distance(self.rib_slot.my_first_pnt) <= self.my_tol {
                    first_ok = true;
                    first_flag = 1;
                } else if sens == 2 && fp.distance(self.rib_slot.my_first_pnt) <= self.my_tol {
                    first_ok = true;
                    first_flag = 2;
                }

                // OCCT L776-827.
                if last_ok {
                    let fpar = cc.default_domain()[0];
                    let lpar = cc.default_domain()[1];
                    let eeee = if !first_ok {
                        if the_previous_edge.is_null() {
                            // OCCT L785: BRepLib_MakeEdge e(cc, fpar, lpar).
                            let mut e_pool = BRep::new();
                            crate::feat::brep_feat_rib_slot::make_edge_cl(
                                &mut e_pool,
                                &cc,
                                fpar,
                                lpar,
                            )
                        } else {
                            let v1 = top_exp_last_vertex(&the_previous_edge, true);
                            let v2 = make_vertex(&mut bb_pool, cc.point_at(lpar));
                            make_edge_c_v_v(&mut bb_pool, &cc, &v1, &v2)
                        }
                    } else if the_previous_edge.is_null() {
                        let v1 = make_vertex(&mut bb_pool, cc.point_at(fpar));
                        make_edge_c_v_v(&mut bb_pool, &cc, &v1, &the_fv)
                    } else {
                        let v1 = top_exp_last_vertex(&the_previous_edge, true);
                        make_edge_c_v_v(&mut bb_pool, &cc, &v1, &the_fv)
                    };

                    the_previous_edge = eeee.clone();
                    bb.add_to_wire(&mut bb_pool, w.clone(), eeee.clone());
                    if counter == 1 {
                        the_fv = top_exp_first_vertex(&eeee, true);
                    }
                    counter1 += 1;
                    new_list_of_edges.push(edg.clone());
                    the_edge = eeee;

                    if dist <= self.my_tol {
                        the_l_edge = edg.clone();
                        let _ = the_l_edge;
                    }
                    the_last_pnt = brep_tool_pnt(&top_exp_last_vertex(&the_edge, true));
                }

                // OCCT L829-836.
                if first_flag == 1 {
                    the_l_edge = edg.clone();
                } else if first_flag == 2 {
                    the_l_edge = the_edge.clone();
                }
                let _ = &the_l_edge;

                // OCCT L838-852.
                if last_ok {
                    it_idx = 0;
                    last_ok = false;
                } else {
                    it_idx += 1;
                }
                sens = 0;
            }
            let _ = (&first_ok, &last_ok, &sens);

            // OCCT L855-908: the SlidMap re-bind.
            let mut slid_map: HashMap<(u64, u32), (Shape, Vec<Shape>)> = HashMap::new();
            if sliding && counter1 > counter {
                // OCCT L863: the edges of w with their positions.
                let w_edges = wire_edges(&w);
                for (ii, e) in w_edges.iter().enumerate() {
                    let ii = ii as i32 + 1;
                    if ii >= counter && ii <= counter1 {
                        for (jj, e2) in new_list_of_edges.iter().enumerate() {
                            let jj = jj as i32 + 1;
                            if jj == (ii - counter + 1) {
                                let items: Vec<(Shape, Vec<Shape>)> =
                                    self.my_slface.values().cloned().collect();
                                for (fac, ledg) in items {
                                    for e1 in &ledg {
                                        if shape_is_same(e1, e2) {
                                            let entry = slid_map
                                                .entry(shape_key(&fac))
                                                .or_insert_with(|| (fac.clone(), Vec::new()));
                                            entry.1.push(e.clone());
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // OCCT L906-907.
            self.my_slface = slid_map;
        }

        // ---Arguments of LocOpe_LinearForm : arguments of the prism
        // sliding

        // OCCT L912-921.
        if sliding {
            let mut f_pool = BRep::new();
            let mut fb = BRepBuilder::new();
            let f = fb.make_face(
                &mut f_pool,
                Some(Surface3::Plane(self.my_pln)),
                Shape::null(),
            );
            // OCCT L916: w.Closed(BRep_Tool::IsClosed(w)).
            let closed = brep_tool_is_closed(&w);
            let _ = closed;
            fb.add_to_face(&mut f_pool, f.clone(), w.clone());
            self.rib_slot.my_skface = f.clone();
            self.rib_slot.my_pbase = self.rib_slot.my_skface.clone();
            self.rib_slot.my_suntil = Shape::null();
        }

        // ---Case without sliding : construction of the face profile

        // OCCT L924-1110.
        if !sliding {
            // OCCT L934-943.
            let explo1 = explorer(&bnd_face, ShapeType::Wire, ShapeType::Shape);
            let www2 = explo1.first().cloned().unwrap_or_else(Shape::null);
            let wi_edges = wire_edges(&www2);
            let mut bu_pool = BRep::new();
            let mut bu = BRepBuilder::new();
            let wiwiwi = bu.make_wire(&mut bu_pool);
            let mut new_v1 = Shape::null();
            let mut new_v2 = Shape::null();
            let mut last_v = Shape::null();

            // OCCT L945-1012.
            for e in &wi_edges {
                let v1 = top_exp_first_vertex(e, true);
                let v2 = top_exp_last_vertex(e, true);
                // OCCT L952.
                let Some((ln, f, l)) = brep_tool_curve(e) else {
                    continue;
                };
                // OCCT L963-964.
                let Some(l2d) = geom_api_to_2d(&ln, plane) else {
                    continue;
                };
                let intcc = Geom2dAPIInterCurveCurve::new(&l2d, &ln2d, CONFUSION);
                let mut vv = Shape::null();

                // OCCT L968-978.
                if intcc.nb_points() > 0 {
                    let p = intcc.point(1);
                    let point = plane_d0(&self.my_pln, p.x, p.y);
                    let par = BRepFeatRibSlot::int_par(&ln, point);
                    if f <= par && l >= par {
                        vv = make_vertex(&mut bu_pool, point);
                    }
                }

                // OCCT L980-1011.
                if vv.is_null() && new_v1.is_null() {
                    continue;
                }
                if !vv.is_null() && new_v1.is_null() {
                    new_v1 = vv;
                    last_v = v2;
                    let ee1 = make_edge_c_v_v(&mut bu_pool, &ln, &new_v1, &last_v);
                    bu.add_to_wire(&mut bu_pool, wiwiwi.clone(), ee1);
                    continue;
                }
                if vv.is_null() && !new_v1.is_null() {
                    let _ee1 = make_edge_c_v_v(&mut bu_pool, &ln, &last_v, &v2);
                    last_v = v2;
                    bu.add_to_wire(&mut bu_pool, wiwiwi.clone(), e.clone());
                    continue;
                }
                if !vv.is_null() && !new_v1.is_null() {
                    new_v2 = vv;
                    let ee1 = make_edge_c_v_v(&mut bu_pool, &ln, &last_v, &new_v2);
                    last_v = new_v2;
                    bu.add_to_wire(&mut bu_pool, wiwiwi.clone(), ee1);
                    let ee2 = make_edge_c_v_v(&mut bu_pool, &ln, &last_v, &new_v1);
                    bu.add_to_wire(&mut bu_pool, wiwiwi.clone(), ee2);
                    break;
                }
            }

            // OCCT L1013-1016.
            let mut nb_pool = BRep::new();
            let mut nb = BRepBuilder::new();
            let mut new_bnd_face =
                make_face_plane_wire(&mut nb, &mut nb_pool, &self.my_pln, &wiwiwi);

            // OCCT L1018-1029: the classification of CheckPnt (architecture
            // note: the rib_slot fclass2d carrier).
            let (paru, parv) = rcad_kernel::math::el::elslib_plane_parameters(
                check_pnt,
                self.my_pln.origin,
                self.my_pln.u_dir,
                self.my_pln.v_dir,
            );
            if crate::feat::brep_feat_rib_slot::brep_top_adaptor_fclass2d_perform(
                &new_bnd_face,
                DVec2::new(paru, parv),
            ) == State::Out
            {
                let c = CutVehicle::with_operation(&bnd_face, &new_bnd_face, BooleanOpType::Cut);
                let c_shape = c.shape().cloned().unwrap_or_else(Shape::null);
                let a_cur_wire = explorer(&c_shape, ShapeType::Wire, ShapeType::Shape)
                    .into_iter()
                    .next()
                    .unwrap_or_else(Shape::null);
                let mut ff_pool = BRep::new();
                let mut ff = BRepBuilder::new();
                new_bnd_face =
                    make_face_plane_wire(&mut ff, &mut ff_pool, &self.my_pln, &a_cur_wire);
            }

            // OCCT L1031-1041.
            if !brep_algo_is_valid(&new_bnd_face) {
                self.set_derived_status(BRepFeatStatusError::InvShape);
                self.rib_slot.not_done();
                return;
            }
            let bnd_face = new_bnd_face;

            // OCCT L1043-1074.
            let mut prof = Shape::null();
            let profile_ok = self.rib_slot.no_sliding_profile(
                &mut prof,
                revol_rib,
                self.my_tol,
                &mut concavite,
                &self.my_pln,
                bnd,
                &bnd_face,
                check_pnt,
                &first_face,
                &last_face,
                &first_vertex,
                &last_vertex,
                &first_edge,
                &last_edge,
                on_first_face,
                on_last_face,
            );
            if !profile_ok {
                self.set_derived_status(BRepFeatStatusError::NoFaceProf);
                self.rib_slot.not_done();
                return;
            }

            // ---Propagation on the faces of the initial shape
            // to find the faces concerned by the rib
            // OCCT L1078-1089.
            let mut falseside = true;
            self.propagate(
                &mut sli_list,
                &prof,
                self.rib_slot.my_first_pnt,
                self.rib_slot.my_last_pnt,
                &mut falseside,
            );
            if !falseside {
                self.set_derived_status(BRepFeatStatusError::FalseSide);
                self.rib_slot.not_done();
                return;
            }

            // OCCT L1091-1109.
            self.my_slface.clear();
            let mut pool_c = BRep::new();
            let comp = bb.make_shell(&mut pool_c);
            for it in &sli_list {
                bb.add_to_shell(&mut pool_c, comp.clone(), it.clone());
            }
            self.rib_slot.my_suntil = comp.clone();
            self.rib_slot.my_skface = prof.clone();
            self.rib_slot.my_pbase = prof;
        }

        // OCCT L1112.
        self.my_sliding = sliding;

        // OCCT L1114-1120.
        for exp in explorer(&self.rib_slot.my_sbase, ShapeType::Face, ShapeType::Shape) {
            self.rib_slot
                .my_map
                .insert(shape_key(&exp), (exp.clone(), Vec::new()));
            self.rib_slot
                .my_map
                .get_mut(&shape_key(&exp))
                .expect("myMap(exp.Current())")
                .1
                .push(exp.clone());
        }
    }

    /// The projected parameter of a point on a Geom_Line (the
    /// GeomAPI_ProjectPointOnCurve.Point(1) parameter carrier of the Init
    /// circle construction).
    fn line_project_parameter(&self, the_line: &Curve3, the_p: DVec3) -> f64 {
        let dom = the_line.default_domain();
        rcad_kernel::base::geom_api::project::closest_point_on_curve_range(
            the_line,
            the_p,
            dom[0],
            dom[1],
            64,
        )
        .param
    }

    /// The derived myStatusError assignments go through the RibSlot member
    /// (the OCCT derived class owns a single myStatusError via RibSlot).
    fn set_derived_status(&mut self, e: BRepFeatStatusError) {
        self.rib_slot.my_status_error = e;
    }

    /// OCCT BRepFeat_MakeRevolutionForm::Add (cxx L1125-1165).
    pub fn add(&mut self, the_e: &Shape, the_f: &Shape) {
        // OCCT L1132-1145.
        if self.my_slface.is_empty() {
            let mut found = false;
            for exp in explorer(&self.rib_slot.my_sbase, ShapeType::Face, ShapeType::Shape) {
                if shape_is_same(&exp, the_f) {
                    found = true;
                    break;
                }
            }
            if !found {
                panic!("Standard_ConstructionError");
            }
            // OCCT L1147-1163.
            if !self.my_slface.contains_key(&shape_key(the_f)) {
                self.my_slface
                    .insert(shape_key(the_f), (the_f.clone(), Vec::new()));
            }
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
            if !found_e {
                self.my_slface
                    .get_mut(&shape_key(the_f))
                    .expect("mySlface(F)")
                    .1
                    .push(the_e.clone());
            }
        }
    }

    /// OCCT BRepFeat_MakeRevolutionForm::Perform (cxx L1169-1438) —
    /// performs a prism from the wire to the plane along the basis shape S.
    pub fn perform(&mut self) {
        // OCCT L1176-1185.
        if self.rib_slot.my_sbase.is_null()
            || self.rib_slot.my_skface.is_null()
            || self.rib_slot.my_pbase.is_null()
        {
            self.rib_slot.my_status_error = BRepFeatStatusError::NotInitialized;
            self.rib_slot.not_done();
            return;
        }

        // OCCT L1187-1198.
        let mut pt = DVec3::ZERO;
        for exx in explorer(&self.rib_slot.my_pbase, ShapeType::Vertex, ShapeType::Shape) {
            pt = brep_tool_pnt(&exx);
            break;
        }

        // OCCT L1200-1249 (architecture difference #6: the rotation branch
        // carries the transform GAP).
        if self.my_angle2 != 0.0 {
            let pbase_trsf = brep_builder_api_transform(
                &self.rib_slot.my_pbase,
                (self.my_axe.location, self.my_axe.direction, -self.my_angle2),
            );
            // OCCT L1207-1243: the myLFMap/mySlface re-attribution.
            let lf_keys: Vec<(u64, u32)> = self.rib_slot.my_lfmap.keys().copied().collect();
            for key in lf_keys {
                let Some((_, list)) = self.rib_slot.my_lfmap.get(&key) else {
                    continue;
                };
                let Some(e1) = list.first().cloned() else {
                    continue;
                };
                let pbase_edges = explorer(&self.rib_slot.my_pbase, ShapeType::Edge, ShapeType::Shape);
                let pbase_trsf_edges = explorer(&pbase_trsf, ShapeType::Edge, ShapeType::Shape);
                for (i, cur) in pbase_edges.iter().enumerate() {
                    if shape_is_same(cur, &e1) {
                        let entry =
                            self.rib_slot.my_lfmap.get_mut(&key).expect("myLFMap");
                        entry.1.clear();
                        if let Some(ex2c) = pbase_trsf_edges.get(i) {
                            entry.1.push(ex2c.clone());
                        }
                        break;
                    }
                }
            }
            let sl_keys: Vec<(u64, u32)> = self.my_slface.keys().copied().collect();
            for key in sl_keys {
                let Some((_, list)) = self.my_slface.get(&key) else {
                    continue;
                };
                let Some(e1) = list.first().cloned() else {
                    continue;
                };
                let pbase_edges = explorer(&self.rib_slot.my_pbase, ShapeType::Edge, ShapeType::Shape);
                let pbase_trsf_edges = explorer(&pbase_trsf, ShapeType::Edge, ShapeType::Shape);
                for (i, cur) in pbase_edges.iter().enumerate() {
                    if shape_is_same(cur, &e1) {
                        let entry = self.my_slface.get_mut(&key).expect("mySlface");
                        entry.1.clear();
                        if let Some(ex2c) = pbase_trsf_edges.get(i) {
                            entry.1.push(ex2c.clone());
                        }
                        break;
                    }
                }
            }
            self.rib_slot.my_pbase = pbase_trsf.clone();
            // OCCT L1245-1248.
            let skface_trsf = brep_builder_api_transform(
                &self.rib_slot.my_skface,
                (self.my_axe.location, self.my_axe.direction, -self.my_angle2),
            );
            self.rib_slot.my_skface = skface_trsf;
        }

        // OCCT L1251-1253.
        let mut the_form = LocOpeRevolutionForm::new();
        the_form.perform(
            &self.rib_slot.my_pbase,
            &self.my_axe,
            self.my_angle1 + self.my_angle2,
        );
        let mut vrai_form = the_form.shape(); // uncut primitive

        // management of descendants

        // OCCT L1256.
        maj_map(
            &self.rib_slot.my_pbase,
            &the_form,
            &mut self.rib_slot.my_map,
            &mut self.rib_slot.my_f_shape,
            &mut self.rib_slot.my_l_shape,
        );

        // OCCT L1258.
        self.rib_slot.my_glued_f.clear();

        // OCCT L1260-1274.
        let pln0 = self.my_pln;
        let normale_dir = self.rib_slot.normal(
            &{
                let mut pool = BRep::new();
                let mut b = BRepBuilder::new();
                make_face_pln_bounds(&mut b, &mut pool, &pln0, -1.0, 1.0, -1.0, 1.0)
            },
            pt,
        );
        let vec1 = self.my_height1 * normale_dir;
        let vec2 = -(self.my_height2) * normale_dir;

        let pln1 = Plane::new(pln0.origin + vec1, pln0.normal);
        let pln2 = Plane::new(pln0.origin + vec2, pln0.normal);

        let mut f_pool = BRep::new();
        let mut fb = BRepBuilder::new();
        let mut f1 = make_face_pln_bounds(&mut fb, &mut f_pool, &pln1, -1.0, 1.0, -1.0, 1.0);
        let mut f2 = make_face_pln_bounds(&mut fb, &mut f_pool, &pln2, -1.0, 1.0, -1.0, 1.0);
        let my_sbase = self.rib_slot.my_sbase.clone();
        brep_feat_face_until(&my_sbase, &mut f1);
        brep_feat_face_until(&my_sbase, &mut f2);

        // OCCT L1276-1284.
        let mut a_si1 = LocOpeCSIntersector::with_shape(&f1);
        let mut a_si2 = LocOpeCSIntersector::with_shape(&f2);
        let normale = Curve3::Line(Line3::new(pt, vec1));
        let scur = vec![Some(normale)];
        a_si1.perform_cur(&scur);
        a_si2.perform_cur(&scur);

        // OCCT L1286-1295.
        if !a_si1.is_done() || !a_si2.is_done() || a_si1.nb_points(1) != 1 || a_si2.nb_points(1) != 1
        {
            self.rib_slot.my_status_error = BRepFeatStatusError::BadIntersect;
            self.rib_slot.not_done();
            return;
        }

        // OCCT L1297-1300.
        let ori1 = a_si1.point(1, 1).orientation();
        let ori2 = top_abs_reverse(a_si2.point(1, 1).orientation());
        let ff1 = a_si1.point(1, 1).face().cloned().unwrap_or_else(Shape::null);
        let ff2 = a_si2.point(1, 1).face().cloned().unwrap_or_else(Shape::null);

        // OCCT L1302-1314.
        let mut pool = BRep::new();
        let mut b = BRepBuilder::new();
        let comp = b.make_compound(&mut pool, Vec::new());
        if let Some(s1) = brep_feat_tool(&f1, &ff1, ori1) {
            b.add_to_compound(&mut pool, comp.clone(), s1);
        }
        if let Some(s2) = brep_feat_tool(&f2, &ff2, ori2) {
            b.add_to_compound(&mut pool, comp.clone(), s2);
        }

        // OCCT L1316-1319: coupe de la nervure par deux plans parallels.
        let tr_p = CutVehicle::with_operation(&vrai_form, &comp, BooleanOpType::Cut);
        let mut sliding_map: HashMap<(u64, u32), (Shape, Vec<Shape>)> = HashMap::new();

        // management of descendants

        // OCCT L1323-1366.
        let map_items: Vec<(Shape, Vec<Shape>)> = self.rib_slot.my_map.values().cloned().collect();
        for (orig, list) in &map_items {
            if list.is_empty() {
                continue;
            }
            let sh = &list[0];
            for fac in explorer(&vrai_form, ShapeType::Face, ShapeType::Shape) {
                let thew = explorer(&fac, ShapeType::Wire, ShapeType::Shape)
                    .into_iter()
                    .next()
                    .unwrap_or_else(Shape::null);
                if shape_is_same(&thew, &self.rib_slot.my_f_shape) {
                    let desfaces = tr_p.modified(&f2);
                    self.rib_slot
                        .my_map
                        .insert(shape_key(&self.rib_slot.my_f_shape), (self.rib_slot.my_f_shape.clone(), desfaces));
                    continue;
                } else if shape_is_same(&thew, &self.rib_slot.my_l_shape) {
                    let desfaces = tr_p.modified(&f1);
                    self.rib_slot
                        .my_map
                        .insert(shape_key(&self.rib_slot.my_l_shape), (self.rib_slot.my_l_shape.clone(), desfaces));
                    continue;
                }
                if shape_is_same(&fac, sh) {
                    if !tr_p.is_deleted(&fac) {
                        let desfaces = tr_p.modified(&fac);
                        if !desfaces.is_empty() {
                            self.rib_slot
                                .my_map
                                .insert(shape_key(orig), (orig.clone(), tr_p.modified(&fac)));
                            break;
                        }
                    }
                }
            }
        }

        // OCCT L1368-1389.
        for fac in explorer(&vrai_form, ShapeType::Face, ShapeType::Shape) {
            sliding_map.insert(shape_key(&fac), (fac.clone(), Vec::new()));
            if tr_p.is_deleted(&fac) {
            } else {
                let desfaces = tr_p.modified(&fac);
                if !desfaces.is_empty() {
                    let e = sliding_map.get_mut(&shape_key(&fac)).expect("SlidingMap(fac)");
                    e.1 = desfaces;
                } else {
                    let e = sliding_map.get_mut(&shape_key(&fac)).expect("SlidingMap(fac)");
                    e.1.push(fac.clone());
                }
            }
        }

        // gestion of faces of sliding

        // OCCT L1392.
        set_glued_faces(&self.my_slface, &the_form, &sliding_map, &mut self.rib_slot.my_glued_f);

        // OCCT L1394.
        vrai_form = tr_p.shape().cloned().unwrap_or_else(Shape::null); // primitive cut

        // OCCT L1396-1403.
        if !self.rib_slot.my_glued_f.is_empty() {
            self.rib_slot.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        } else {
            self.rib_slot.my_perf_selection = BRepFeatPerfSelection::SelectionSh;
        }

        // OCCT L1405-1419.
        for exx in explorer(&self.rib_slot.my_pbase, ShapeType::Edge, ShapeType::Shape) {
            if !self.rib_slot.my_map.contains_key(&shape_key(&exx)) {
                self.rib_slot.my_status_error = BRepFeatStatusError::IncSlidFace;
                self.rib_slot.not_done();
                return;
            }
        }

        // OCCT L1421.
        self.rib_slot.my_gshape = vrai_form;

        // OCCT L1423-1435.
        if !self.rib_slot.my_glued_f.is_empty() && !self.rib_slot.my_suntil.is_null() {
            self.rib_slot.my_status_error = BRepFeatStatusError::InvShape;
            self.rib_slot.not_done();
            return;
        }

        // OCCT L1437: LFPerform() — topological reconstruction.
        self.rib_slot.lf_perform();
    }

    /// OCCT BRepFeat_MakeRevolutionForm::Propagate (cxx L1446-1933) —
    /// propagation on the faces of the initial shape, find faces concerned
    /// by the rib.
    pub fn propagate(
        &mut self,
        sli_list: &mut Vec<Shape>,
        fac: &Shape,
        firstpnt: DVec3,
        lastpnt: DVec3,
        falseside: &mut bool,
    ) -> bool {
        // OCCT L1457-1458.
        let mut firstpoint = firstpnt;
        let mut lastpoint = lastpnt;

        // OCCT L1460-1465.
        let result = true;
        let mut current_face = sli_list.first().cloned().unwrap_or_else(Shape::null);
        let save_face = current_face.clone();
        let mut last_ok = false;
        let mut first_ok = false;
        let mut v1 = Shape::null();
        let mut v2 = Shape::null();
        // OCCT L1467-1469: BRepAlgoAPI_Section sect(fac, CurrentFace, false)
        // + Approximation(true) + Build().
        let mut sect = SectionOp::from_shapes(fac.clone(), current_face.clone());
        sect.approximation(true);
        sect.build();
        let mut e = Shape::null();
        let mut e1 = Shape::null();
        // OCCT L1474-1486.
        let sect_shape = sect.algo.bs.result.clone().unwrap_or_else(Shape::null);
        let mut ii = 0i32;
        for ex in explorer(&sect_shape, ShapeType::Edge, ShapeType::Shape) {
            ii += 1;
            if ii == 1 {
                e = ex;
            } else {
                e1 = ex;
                break;
            }
        }
        // OCCT L1487-1491.
        if e.is_null() {
            *falseside = false;
            return false;
        }
        // OCCT L1493-1530.
        if !e1.is_null() {
            self.my_list_of_edges.clear();
            self.my_slface
                .insert(shape_key(&current_face), (current_face.clone(), Vec::new()));
            self.my_slface
                .get_mut(&shape_key(&current_face))
                .expect("mySlface(CurrentFace)")
                .1
                .push(e1.clone());
            self.my_list_of_edges.push(e1.clone());

            v1 = top_exp_first_vertex(&e1, true);
            v2 = top_exp_last_vertex(&e1, true);

            let fp = brep_tool_pnt(&v1);
            let lp = brep_tool_pnt(&v2);

            let a_tol_v1 = brep_tool_tolerance(&v1);
            let a_tol_v2 = brep_tool_tolerance(&v2);

            if fp.distance(firstpoint) <= a_tol_v1 || fp.distance(lastpoint) <= a_tol_v1 {
                first_ok = true;
            }
            if lp.distance(firstpoint) <= a_tol_v2 || lp.distance(lastpoint) <= a_tol_v2 {
                last_ok = true;
            }
            if last_ok && first_ok {
                return result;
            }
            self.my_list_of_edges.clear();
        }
        // OCCT L1532-1629.
        if !e1.is_null() {
            self.my_list_of_edges.clear();
            self.my_slface
                .insert(shape_key(&current_face), (current_face.clone(), Vec::new()));
            {
                let entry = self
                    .my_slface
                    .get_mut(&shape_key(&current_face))
                    .expect("mySlface(CurrentFace)");
                entry.1.push(e.clone());
                entry.1.push(e1.clone());
            }
            self.my_list_of_edges.push(e.clone());

            v1 = top_exp_first_vertex(&e, true);
            v2 = top_exp_last_vertex(&e, true);
            let v3 = top_exp_first_vertex(&e1, true);
            let v4 = top_exp_last_vertex(&e1, true);
            let p1 = brep_tool_pnt(&v1);
            let mut fp = p1;
            let p2 = brep_tool_pnt(&v2);
            let mut lp = p2;
            let p3 = brep_tool_pnt(&v3);
            let p4 = brep_tool_pnt(&v4);
            // OCCT L1556-1628.
            if p1.distance(firstpoint) <= brep_tool_tolerance(&v1) {
                if p3.distance(lastpoint) <= brep_tool_tolerance(&v3) {
                    first_ok = true;
                    lastpoint = p4;
                } else if p4.distance(lastpoint) <= brep_tool_tolerance(&v4) {
                    first_ok = true;
                    lastpoint = p3;
                } else {
                    e1 = Shape::null();
                }
            } else if p1.distance(lastpoint) <= brep_tool_tolerance(&v1) {
                if p3.distance(firstpoint) <= brep_tool_tolerance(&v3) {
                    first_ok = true;
                    firstpoint = p4;
                } else if p4.distance(firstpoint) <= brep_tool_tolerance(&v4) {
                    first_ok = true;
                    firstpoint = p3;
                } else {
                    e1 = Shape::null();
                }
            } else if p2.distance(firstpoint) <= brep_tool_tolerance(&v2) {
                if p3.distance(lastpoint) <= brep_tool_tolerance(&v3) {
                    last_ok = true;
                    lastpoint = p4;
                } else if p4.distance(lastpoint) <= brep_tool_tolerance(&v4) {
                    last_ok = true;
                    lastpoint = p3;
                } else {
                    e1 = Shape::null();
                }
            } else if p2.distance(lastpoint) <= brep_tool_tolerance(&v2) {
                if p3.distance(firstpoint) <= brep_tool_tolerance(&v3) {
                    last_ok = true;
                    firstpoint = p4;
                } else if p4.distance(firstpoint) <= brep_tool_tolerance(&v4) {
                    last_ok = true;
                    firstpoint = p3;
                } else {
                    e1 = Shape::null();
                }
            } else {
                e = e1.clone();
                e1 = Shape::null();
            }
            let _ = (fp, lp);
        }
        // OCCT L1630-1660.
        if e1.is_null() {
            self.my_list_of_edges.clear();
            self.my_slface
                .insert(shape_key(&current_face), (current_face.clone(), Vec::new()));
            self.my_slface
                .get_mut(&shape_key(&current_face))
                .expect("mySlface(CurrentFace)")
                .1
                .push(e.clone());
            self.my_list_of_edges.push(e.clone());

            v1 = top_exp_first_vertex(&e, true);
            v2 = top_exp_last_vertex(&e, true);

            let fp = brep_tool_pnt(&v1);
            let lp = brep_tool_pnt(&v2);

            if fp.distance(firstpoint) <= brep_tool_tolerance(&v1)
                || fp.distance(lastpoint) <= brep_tool_tolerance(&v1)
            {
                first_ok = true;
            }
            if lp.distance(firstpoint) <= brep_tool_tolerance(&v2)
                || lp.distance(lastpoint) <= brep_tool_tolerance(&v2)
            {
                last_ok = true;
            }
            if last_ok && first_ok {
                return result;
            }
        }

        // OCCT L1662-1671.
        let mut mapedges: IndexMap<(u64, u32), (Shape, Vec<Shape>)> = IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &self.rib_slot.my_sbase,
            ShapeType::Edge,
            ShapeType::Face,
            &mut mapedges,
        );
        let mut first_edge = Shape::null();
        let mut v_previous = Shape::null();
        let mut v_preprevious = Shape::null();

        // OCCT L1673-1802: the while(!FirstOK) pass.
        while !first_ok {
            // find edge connected to v1:
            let pt = if !v1.is_null() {
                brep_tool_pnt(&v1)
            } else {
                DVec3::ZERO
            };
            if (!v_previous.is_null() && brep_tool_pnt(&v_previous).distance(pt) <= self.my_tol)
                || (!v_preprevious.is_null()
                    && brep_tool_pnt(&v_preprevious).distance(pt) <= self.my_tol)
            {
                *falseside = false;
                return false;
            }

            // OCCT L1699-1726.
            for ex in explorer(&current_face, ShapeType::Edge, ShapeType::Shape) {
                let a_cur_edge = ex;
                // OCCT L1703: BRepExtrema_ExtPC projF(v1, aCurEdge).
                let Some((c, f, l)) = brep_tool_curve(&a_cur_edge) else {
                    continue;
                };
                let proj_f = ExtPC::new(brep_tool_pnt(&v1), &c, CONFUSION, f, l);
                if proj_f.is_done() && proj_f.nb_ext() >= 1 {
                    let mut dist2min = f64::MAX; // OCCT: RealLast()
                    let mut index = 0usize;
                    for sol in 1..=proj_f.nb_ext() {
                        if proj_f.square_distance(sol) <= dist2min {
                            index = sol;
                            dist2min = proj_f.square_distance(sol);
                        }
                    }
                    if index != 0
                        && dist2min
                            <= brep_tool_tolerance(&a_cur_edge) * brep_tool_tolerance(&a_cur_edge)
                    {
                        first_edge = a_cur_edge;
                        break;
                    }
                }
            }

            // OCCT L1728-1739.
            if let Some((_, l)) = mapedges.get(&shape_key(&first_edge)) {
                for ff in l {
                    if !shape_is_same(ff, &current_face) {
                        // current_face = FF (the copy below).
                        break;
                    }
                }
            }
            let mut current_face_new = current_face.clone();
            if let Some((_, l)) = mapedges.get(&shape_key(&first_edge)) {
                for ff in l {
                    if !shape_is_same(ff, &current_face) {
                        current_face_new = ff.clone();
                        break;
                    }
                }
            }
            current_face = current_face_new;

            // OCCT L1741-1756.
            let mut sectf = SectionOp::from_shapes(fac.clone(), current_face.clone());
            sectf.approximation(true);
            sectf.build();
            let mut edg1 = Shape::null();
            let sectf_shape = sectf.algo.bs.result.clone().unwrap_or_else(Shape::null);
            for ex in explorer(&sectf_shape, ShapeType::Edge, ShapeType::Shape) {
                edg1 = ex;
                let ppp1 = brep_tool_pnt(&top_exp_first_vertex(&edg1, true));
                let ppp2 = brep_tool_pnt(&top_exp_last_vertex(&edg1, true));
                if ppp1.distance(brep_tool_pnt(&v1)) <= brep_tool_tolerance(&v1)
                    || ppp2.distance(brep_tool_pnt(&v1)) <= brep_tool_tolerance(&v1)
                {
                    break;
                }
            }

            // OCCT L1758-1770.
            self.my_slface
                .insert(shape_key(&current_face), (current_face.clone(), Vec::new()));
            self.my_slface
                .get_mut(&shape_key(&current_face))
                .expect("mySlface(CurrentFace)")
                .1
                .push(edg1.clone());
            self.my_list_of_edges.push(edg1.clone());

            if !edg1.is_null() {
                sli_list.insert(0, current_face.clone()); // Prepend
            } else {
                return false;
            }

            // OCCT L1772-1793.
            let vert = top_exp_first_vertex(&edg1, true);
            let pp = brep_tool_pnt(&vert);
            let fp = brep_tool_pnt(&v1);
            let mut tol = brep_tool_tolerance(&edg1);
            let tol1 = brep_tool_tolerance(&v1);
            if tol1 > tol {
                tol = tol1;
            }
            let dist = pp.distance(fp);
            let v_preprevious_new = v_previous.clone();
            let v_previous_new = v1.clone();
            if dist <= tol {
                v1 = top_exp_last_vertex(&edg1, true);
            } else {
                v1 = vert;
            }
            v_preprevious = v_preprevious_new;
            v_previous = v_previous_new;

            // OCCT L1795-1801.
            let fp = brep_tool_pnt(&v1);
            if fp.distance(firstpoint) <= brep_tool_tolerance(&v1)
                || fp.distance(lastpoint) <= brep_tool_tolerance(&v1)
            {
                first_ok = true;
            }
        }

        // OCCT L1804-1806.
        let mut current_face = save_face;
        v_previous = Shape::null();
        v_preprevious = Shape::null();

        // OCCT L1808-1927: the while(!LastOK) pass.
        while !last_ok {
            // find edge connected to v2:
            let pt = if !v2.is_null() {
                brep_tool_pnt(&v2)
            } else {
                DVec3::ZERO
            };
            if (!v_previous.is_null() && brep_tool_pnt(&v_previous).distance(pt) <= self.my_tol)
                || (!v_preprevious.is_null()
                    && brep_tool_pnt(&v_preprevious).distance(pt) <= self.my_tol)
            {
                *falseside = false;
                return false;
            }

            // OCCT L1834-1860.
            for ex in explorer(&current_face, ShapeType::Edge, ShapeType::Shape) {
                let a_cur_edge = ex;
                let Some((c, f, l)) = brep_tool_curve(&a_cur_edge) else {
                    continue;
                };
                let proj_f = ExtPC::new(brep_tool_pnt(&v2), &c, CONFUSION, f, l);
                if proj_f.is_done() && proj_f.nb_ext() >= 1 {
                    let mut dist2min = f64::MAX;
                    let mut index = 0usize;
                    for sol in 1..=proj_f.nb_ext() {
                        if proj_f.square_distance(sol) <= dist2min {
                            index = sol;
                            dist2min = proj_f.square_distance(sol);
                        }
                    }
                    if index != 0
                        && dist2min
                            <= brep_tool_tolerance(&a_cur_edge) * brep_tool_tolerance(&a_cur_edge)
                    {
                        first_edge = a_cur_edge;
                        break;
                    }
                }
            }

            // OCCT L1862-1873.
            let mut current_face_new = current_face.clone();
            if let Some((_, l)) = mapedges.get(&shape_key(&first_edge)) {
                for ff in l {
                    if !shape_is_same(ff, &current_face) {
                        current_face_new = ff.clone();
                        break;
                    }
                }
            }
            current_face = current_face_new;

            // OCCT L1875-1892.
            let ii = 0i32;
            let _ = ii;
            let mut sectf = SectionOp::from_shapes(fac.clone(), current_face.clone());
            sectf.approximation(true);
            sectf.build();
            let mut edg2 = Shape::null();
            let sectf_shape = sectf.algo.bs.result.clone().unwrap_or_else(Shape::null);
            for ex in explorer(&sectf_shape, ShapeType::Edge, ShapeType::Shape) {
                edg2 = ex;
                let ppp1 = brep_tool_pnt(&top_exp_first_vertex(&edg2, true));
                let ppp2 = brep_tool_pnt(&top_exp_last_vertex(&edg2, true));
                if ppp1.distance(brep_tool_pnt(&v2)) <= brep_tool_tolerance(&v2)
                    || ppp2.distance(brep_tool_pnt(&v2)) <= brep_tool_tolerance(&v2)
                {
                    break;
                }
            }

            // OCCT L1893-1905.
            self.my_slface
                .insert(shape_key(&current_face), (current_face.clone(), Vec::new()));
            self.my_slface
                .get_mut(&shape_key(&current_face))
                .expect("mySlface(CurrentFace)")
                .1
                .push(edg2.clone());
            self.my_list_of_edges.push(edg2.clone());

            if !edg2.is_null() {
                sli_list.push(current_face.clone());
            } else {
                return false;
            }

            // OCCT L1907-1920.
            let vert = top_exp_first_vertex(&edg2, true);
            let pp = brep_tool_pnt(&vert);
            let fp = brep_tool_pnt(&v2);
            if pp.distance(fp) <= brep_tool_tolerance(&v2) {
                v_preprevious = v_previous.clone();
                v_previous = v2.clone();
                v2 = top_exp_last_vertex(&edg2, true);
            } else {
                v2 = vert;
            }

            // OCCT L1922-1926.
            let fp = brep_tool_pnt(&v2);
            if fp.distance(firstpoint) <= brep_tool_tolerance(&v2)
                || fp.distance(lastpoint) <= brep_tool_tolerance(&v2)
            {
                last_ok = true;
            }
        }
        // OCCT L1928-1932.
        if !e1.is_null() {
            self.my_list_of_edges.push(e1);
        }
        result
    }
}

/// OCCT static MajMap (cxx L1937-1980) — the RevolutionForm variant.
fn maj_map(
    the_b: &Shape,
    the_p: &LocOpeRevolutionForm,
    the_map: &mut HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_f_shape: &mut Shape,
    the_l_shape: &mut Shape,
) {
    // OCCT L1945-1956.
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
    // OCCT L1958-1969.
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
    // OCCT L1971-1979.
    for exp in explorer(the_b, ShapeType::Edge, ShapeType::Shape) {
        if !the_map.contains_key(&shape_key(&exp)) {
            let shapes = the_p.shapes(&exp).clone();
            the_map.insert(shape_key(&exp), (exp.clone(), shapes));
        }
    }
}

/// OCCT static SetGluedFaces (cxx L1984-2030) — the gestion of faces of
/// sliding.
fn set_glued_faces(
    the_slmap: &HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_prism: &LocOpeRevolutionForm,
    sliding_map: &HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_map: &mut HashMap<(u64, u32), (Shape, Shape)>,
) {
    // Slidings
    // OCCT L1993-2029.
    if !the_slmap.is_empty() {
        for (fac_key, (fac, ledg)) in the_slmap.iter() {
            let fac = fac;
            let _ = fac_key;
            for it in ledg {
                let gfac = the_prism.shapes(it);
                if gfac.len() != 1 {
                    // OCCT L2007: the debug print only.
                }
                for (ff_key, (ff, lfaces)) in sliding_map.iter() {
                    let ff = ff;
                    let _ = ff_key;
                    if lfaces.is_empty() {
                        continue;
                    }
                    let fff = &lfaces[0];
                    if let Some(g_first) = gfac.first() {
                        if shape_is_same(g_first, ff) {
                            the_map.insert(shape_key(fff), (fff.clone(), fac.clone()));
                        }
                    }
                }
            }
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
