// OCCT BRepFeat_MakeLinearForm.hxx L17-136 + BRepFeat_MakeLinearForm.cxx
// L17-1366 + BRepFeat_MakeLinearForm.lxx L17-55 — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeLinearForm.hxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeLinearForm.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKFeat/BRepFeat/BRepFeat_MakeLinearForm.lxx
//
// OCCT inheritance chain (BRepFeat_MakeLinearForm.hxx L52):
//   BRepFeat_MakeLinearForm : BRepFeat_RibSlot : BRepBuilderAPI_MakeShape
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
// 2. BRepAlgoAPI_Common is the CutVehicle of brep_feat_form_2 (with
//    BooleanOpType::Intersection); BRepAlgoAPI_Section is the SectionOp of
//    bop/brep_algo_api (the true BOPAlgo_Section vehicle).
// 3. BRepExtrema_ExtCF (Init, cxx L307-317 / L333-345) and BRepExtrema_ExtPF
//    (Propagate, cxx L1211-1227) are not translated; the carriers report
//    NbExt() = 0 / not-done, so the OCCT guards drive Sliding = false /
//    result = false into the no-sliding branches (the MakeCylindricalHole
//    arch.-difference #1 model).
// 4. BRepPrimAPI_MakeBox (Init, cxx L208-209) is the leaf re-host below
//    (brep_prim_make_box): the modeling MakeBox pool builder + the root
//    Solid read (the TKPrim BRepPrimAPI batch can re-home it).
// 5. BRepBuilderAPI_Transform (Init, cxx L189-194, the centred-rib
//    translation branch) is the topalgo/brep_builderapi_transform port; the
//    rigid translation takes the OCCT Perform "Moved" branch, materialised
//    on the wire TShapes because the feat BRep_Tool re-hosts read geometry
//    flat (arch. difference #1 of brep_algo/tool.rs).
// 6. gp_Vec::IsEqual / gp_Vec::Angle are the small analytic carriers below.
// 7. BRepTools_WireExplorer maps to the wire edge list (the same reduction
//    as in brep_feat_make_revolution_form.rs).
// 8. ExtremeFaces / CheckPoint / SlidingProfile / NoSlidingProfile /
//    HeightMax / IntPar / Normal are the RibSlot re-hosts consumed through
//    the rib_slot sub-object; BB.UpdateVertex (cxx L630 / L1247-1280) is
//    the rcad vertex-tolerance carrier (the Arc-shared TShape flag is not
//    mutated through the handle — the kernel edge flags keep the value).
// 9. BRepFeat::FaceUntil is the brep_feat_form_2 re-host.

use crate::bop::algo::builder::BooleanOpType;
use crate::bop::brep_algo_api::SectionOp;
use crate::feat::brep_feat_builder::explorer;
use crate::feat::brep_feat_form_2::{brep_feat_tool, CutVehicle};
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::feat::brep_feat_rib_slot::{
    brep_tool_curve, brep_tool_is_closed, brep_tool_pnt, brep_tool_surface,
    make_edge_c_v_v, make_edge_p_p, make_vertex,
    top_exp_first_vertex, top_exp_last_vertex, BRepFeatRibSlot,
    Geom2dAPIInterCurveCurve, GeomAPIProjectPointOnCurve,
};
use crate::feat::brep_feat_status::{BRepFeatPerfSelection, BRepFeatStatusError};
use crate::feat::loc_ope_linear_form::LocOpeLinearForm;
use glam::DVec3;
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::geom::{Plane, Surface3, TrimmedSurface};
use rcad_kernel::math::gp::Trsf;
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, TShape};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType};
use rcad_kernel::CurveEval;
use std::collections::HashMap;

/// OCCT BRepTools_WireExplorer — the wire edge list (architecture
/// difference #7).
fn wire_edges(the_wire: &Shape) -> Vec<Shape> {
    match the_wire.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT BRepExtrema_ExtPF(ve, fac) carrier (architecture difference #3) —
/// IsDone() = false keeps the OCCT guard structure.
pub(crate) struct BRepExtremaExtPFCarrier;

impl BRepExtremaExtPFCarrier {
    fn new(_the_ve: &Shape, _the_fac: &Shape) -> Self {
        BRepExtremaExtPFCarrier
    }
    pub(crate) fn is_done(&self) -> bool {
        false
    }
    #[allow(dead_code)]
    fn point(&self, _n: i32) -> DVec3 {
        DVec3::ZERO
    }
    pub(crate) fn perform(&mut self, _the_ve: &Shape, _the_fac: &Shape) {}
}

impl BRepExtremaExtPFCarrier {
    /// OCCT BRepExtrema_ExtPF() (BRepExtrema_ExtPF.hxx L44) — the default
    /// constructor of the Initialize(TheFace) + Perform(Vertex, Face) call
    /// pattern (the LocOpe_Gluer::AddEdges site).
    pub(crate) fn new_default() -> Self {
        BRepExtremaExtPFCarrier
    }

    /// OCCT BRepExtrema_ExtPF::Initialize(TheFace) (cxx L42-64).
    pub(crate) fn initialize(&mut self, _the_face: &Shape) {}

    /// OCCT BRepExtrema_ExtPF::NbExt() (hxx L57) — the carrier reports the
    /// not-done count (mySqDist stays cleared, the OCCT early-out of cxx
    /// L74-77 / L82).
    pub(crate) fn nb_ext(&self) -> i32 {
        0
    }

    /// OCCT BRepExtrema_ExtPF::SquareDistance(N) (hxx L53).
    pub(crate) fn square_distance(&self, _n: i32) -> f64 {
        0.0
    }
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

/// OCCT gp_Vec::IsEqual(theOther, theLinTol, theAngTol) — the magnitude and
/// direction equality (gp_Vec.hxx L105-125).
fn gp_vec_is_equal(v1: DVec3, v2: DVec3, lin_tol: f64, ang_tol: f64) -> bool {
    let dist = (v1 - v2).length();
    if dist <= lin_tol {
        return true;
    }
    let m1 = v1.length();
    let m2 = v2.length();
    if m1 <= lin_tol || m2 <= lin_tol {
        return false;
    }
    let an_ang = DVec3::angle_between(v1, v2);
    let _ = ang_tol;
    (an_ang <= ang_tol && dist <= lin_tol * m1) || dist <= lin_tol
}

/// OCCT BRepExtrema_ExtCF carrier (architecture difference #3) — NbExt() = 0
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

/// OCCT BRepPrimAPI_MakeBox(P1, P2).Solid() (BRepPrimAPI_MakeBox.cxx L78-86):
/// myWedge = BRepPrim_Wedge(gp_Ax2(pmin(P1, P2), gp_Dir(Z), gp_Dir(X)),
/// Abs(P2.X()-P1.X()), Abs(P2.Y()-P1.Y()), Abs(P2.Z()-P1.Z())) — the
/// axis-aligned box solid between the two corners.  The rcad leaf re-host
/// consumes the modeling MakeBox pool builder and reads the root Solid from
/// the pool (the TKPrim BRepPrimAPI batch can re-home it).
fn brep_prim_make_box(p1: DVec3, p2: DVec3) -> Shape {
    // OCCT pmin(P1, P2).
    let pmin = DVec3::new(p1.x.min(p2.x), p1.y.min(p2.y), p1.z.min(p2.z));
    let brep = rcad_modeling::make_box_brep(
        pmin,
        glam::DVec3::X,
        glam::DVec3::Y,
        (p2.x - p1.x).abs(),
        (p2.y - p1.y).abs(),
        (p2.z - p1.z).abs(),
    )
    .unwrap_or_else(|_| BRep::new());
    // OCCT Solid().
    brep.tshapes
        .iter()
        .enumerate()
        .rev()
        .find_map(|(i, ts)| match ts.as_ref() {
            TShape::Solid(_) => {
                Some(Shape::from_parts(ts.clone(), i, 0, Orientation::Forward))
            }
            _ => None,
        })
        .unwrap_or_else(Shape::null)
}

/// OCCT BRepLib_MakeFace(Pln, U1, U2, V1, V2) — the face on the plane with
/// the uv bounds (the rectangular trimmed surface carrier, the
/// brep_feat_face_until model).
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
fn make_face_of_plane_wire(
    the_b: &mut BRepBuilder,
    the_pool: &mut BRep,
    the_pln: &Plane,
    the_wire: &Shape,
) -> Shape {
    the_b.make_face(the_pool, Some(Surface3::Plane(*the_pln)), the_wire.clone())
}

/// OCCT BB.MakeWire() — the empty wire of the pool.
fn make_wire_local(the_pool: &mut BRep) -> Shape {
    let mut b = BRepBuilder::new();
    b.make_wire(the_pool)
}

/// OCCT BB.Add(w, edge).
fn add_to_wire_local(the_pool: &mut BRep, w: &Shape, e: &Shape) {
    let mut b = BRepBuilder::new();
    b.add_to_wire(the_pool, w.clone(), e.clone())
}

/// OCCT BB.MakeFace(F, myPln, myTol) — the plane face.
fn make_face_of_plane_local(the_pool: &mut BRep, pln: &Plane) -> Shape {
    let mut b = BRepBuilder::new();
    b.make_face(the_pool, Some(Surface3::Plane(*pln)), Shape::null())
}

/// OCCT BB.Add(F, w).
fn add_face_wire_local(the_pool: &mut BRep, f: &Shape, w: &Shape) {
    let mut b = BRepBuilder::new();
    b.add_to_face(the_pool, f.clone(), w.clone())
}

/// OCCT BB.MakeShell(comp).
fn make_shell_local(the_pool: &mut BRep) -> Shape {
    let mut b = BRepBuilder::new();
    b.make_shell(the_pool)
}

/// OCCT BB.Add(comp, face).
fn add_to_shell_local(the_pool: &mut BRep, comp: &Shape, face: &Shape) {
    let mut b = BRepBuilder::new();
    b.add_to_shell(the_pool, comp.clone(), face.clone())
}

/// OCCT BRepFeat_MakeLinearForm — describes functions to build linear form
/// features (BRepFeat_MakeLinearForm.hxx L36-51).
pub struct BRepFeatMakeLinearForm {
    /// The BRepFeat_RibSlot base sub-object (the composition carrier).
    pub rib_slot: BRepFeatRibSlot,
    my_dir: DVec3, // OCCT: myDir (gp_Vec)
    my_dir1: DVec3, // OCCT: myDir1 (gp_Vec)
    my_pln: Plane, // OCCT: myPln (Handle(Geom_Plane))
    #[allow(dead_code)]
    my_bnd: f64, // OCCT: myBnd
    // OCCT: mySlface (architecture difference #1).
    my_slface: HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    #[allow(dead_code)]
    my_list_of_edges: Vec<Shape>, // OCCT: myListOfEdges
    #[allow(dead_code)]
    my_tol: f64, // OCCT: myTol
}

impl BRepFeatMakeLinearForm {
    /// OCCT BRepFeat_MakeLinearForm::BRepFeat_MakeLinearForm() (lxx: the
    /// empty constructor).
    pub fn new() -> Self {
        BRepFeatMakeLinearForm {
            rib_slot: BRepFeatRibSlot::new(),
            my_dir: DVec3::ZERO,
            my_dir1: DVec3::ZERO,
            my_pln: Plane::new(DVec3::ZERO, DVec3::Z),
            my_bnd: 0.0,
            my_slface: HashMap::new(),
            my_list_of_edges: Vec::new(),
            my_tol: 0.0,
        }
    }

    /// OCCT BRepFeat_MakeLinearForm::BRepFeat_MakeLinearForm(Sbase, W, P,
    /// Direction, Direction1, Fuse, Modify) (hxx) — Rust has no overloading:
    /// the `with_init` suffix.
    #[allow(clippy::too_many_arguments)]
    pub fn with_init(
        sbase: &Shape,
        w: &Shape,
        plane: &Plane,
        direc: DVec3,
        direc1: DVec3,
        fuse: i32,
        modify: bool,
    ) -> Self {
        let mut res = BRepFeatMakeLinearForm::new();
        res.init(sbase, w, plane, direc, direc1, fuse, modify);
        res
    }

    /// OCCT BRepFeat_MakeLinearForm::Init (cxx L84-852).
    #[allow(clippy::too_many_arguments)]
    #[allow(unused_assignments)]
    pub fn init(
        &mut self,
        sbase: &Shape,
        w: &Shape,
        plane: &Plane,
        direc: DVec3,
        direc1: DVec3,
        mode: i32,
        modify: bool,
    ) {
        // OCCT L97-99.
        let revol_rib = false;
        self.rib_slot.done();
        self.rib_slot.my_generated.clear();

        // modify = 0 if there is no intention to make sliding
        //        = 1 if one tries to make sliding
        // OCCT L103.
        let mut sliding = modify;
        self.rib_slot.my_lfmap.clear();

        // OCCT L106-117.
        self.rib_slot.my_shape = None;
        self.rib_slot.my_map.clear();
        self.rib_slot.my_f_shape = Shape::null();
        self.rib_slot.my_l_shape = Shape::null();
        self.rib_slot.my_sbase = sbase.clone();
        self.rib_slot.my_skface = Shape::null();
        self.rib_slot.my_pbase = Shape::null();
        self.rib_slot.my_gshape = Shape::null();
        self.rib_slot.my_suntil = Shape::null();
        self.my_list_of_edges.clear();
        self.my_slface.clear();

        // OCCT L119-124.
        self.rib_slot.my_wire = w.clone();
        self.my_dir = direc;
        self.my_dir1 = direc1;
        self.my_pln = *plane;

        // OCCT L126.
        self.rib_slot.my_fuse = mode != 0;

        // ---Determine Tolerance : max tolerance on parameters

        // OCCT L138.
        self.my_tol = CONFUSION;
        // OCCT L140-149.
        for exx in explorer(&self.rib_slot.my_wire, ShapeType::Vertex, ShapeType::Shape) {
            let tol = brep_tool_tolerance(&exx);
            if tol > self.my_tol {
                self.my_tol = tol;
            }
        }
        // OCCT L151-159.
        for exx in explorer(sbase, ShapeType::Vertex, ShapeType::Shape) {
            let tol = brep_tool_tolerance(&exx);
            if tol > self.my_tol {
                self.my_tol = tol;
            }
        }

        // ---Control of directions
        //    the wire should be in the rib

        // OCCT L163-177.
        let nulldir = DVec3::ZERO;
        if !gp_vec_is_equal(self.my_dir1, nulldir, self.my_tol, self.my_tol) {
            let ang = DVec3::angle_between(self.my_dir1, self.my_dir);
            if ang != std::f64::consts::PI {
                self.rib_slot.my_status_error = BRepFeatStatusError::BadDirect;
                self.rib_slot.not_done();
                return;
            }
        } else {
            // Rib is centre in the middle of translation
            // OCCT L186: DirTranslation = (Direc + Direc1) * 0.5.
            let dir_translation = (direc + direc1) * 0.5;
            // OCCT L187-188: gp_Trsf T; T.SetTranslation(DirTranslation).
            let mut t = Trsf::identity();
            t.set_translation(dir_translation);
            // OCCT L189-191: BRepBuilderAPI_Transform trf(T);
            // trf.Perform(myWire); myWire = TopoDS::Wire(trf.Shape()).
            // A rigid translation takes the Perform "Moved" branch (cxx L58-63:
            // the transform rides on the result location).  rcad arch.
            // difference #1 (the feat BRep_Tool re-hosts read TShape geometry
            // flat, see the module header #5) makes the located-but-
            // untransformed wire unobservable, so the transformation is
            // materialised on the wire's TShapes with TShape identity
            // preserved (in-place mutation = the Moved same-handle model).
            crate::topalgo::brep_builderapi_transform::perform_shape(
                &self.rib_slot.my_wire,
                &t,
            );
            // OCCT L192-193: myDir = Direc - DirTranslation;
            //                 myDir1 = Direc1 - DirTranslation.
            self.my_dir = direc - dir_translation;
            self.my_dir1 = direc1 - dir_translation;
            // OCCT L194: myPln->Transform(T) — gp_Pln::Transform.
            self.my_pln.transform(&t);
        }

        // ---Calculate bounding box

        // OCCT L200-206.
        let mut the_list: Vec<Shape> = Vec::new();
        let mut u = Shape::null(); // OCCT: U.Nullify()
        let mut first_corner = DVec3::ZERO;
        let mut last_corner = DVec3::ZERO;
        let bnd = self.rib_slot.height_max(&self.rib_slot.my_sbase, &mut u, &mut first_corner, &mut last_corner);
        self.my_bnd = bnd;

        // OCCT L208-209 (architecture difference #4).
        let bnd_box = brep_prim_make_box(first_corner, last_corner);

        // ---Construction of the face workplane (section bounding box)

        // OCCT L212-213.
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

        // OCCT L215-221.
        let plane_s =
            CutVehicle::with_operation(&bnd_box, &plane_face, BooleanOpType::Intersection);
        let plane_sect = plane_s.shape().cloned().unwrap_or_else(Shape::null);
        let www = explorer(&plane_sect, ShapeType::Wire, ShapeType::Shape)
            .into_iter()
            .next()
            .unwrap_or_else(Shape::null);
        let bnd_face = make_face_of_plane_wire(&mut bb, &mut pool, &self.my_pln, &www);

        // ---Find support faces of the rib

        // OCCT L224-234.
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

        // OCCT L236-261.
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
            self.rib_slot.my_status_error = BRepFeatStatusError::NoExtFace;
            self.rib_slot.not_done();
            return;
        }

        // ---Proofing Point for the side of the wire to be filled - side
        // material

        // OCCT L264.
        let check_pnt = self.rib_slot.check_point(&first_edge, bnd / 10.0, &self.my_pln);

        // ---Control sliding valuable
        // Many cases when the sliding is abandoned

        // OCCT L270.
        let mut concavite = 3i32; // a priori the profile is not concave

        // OCCT L272-273.
        self.rib_slot.my_first_pnt = brep_tool_pnt(&first_vertex);
        self.rib_slot.my_last_pnt = brep_tool_pnt(&last_vertex);

        // SliList : list of faces concerned by the rib

        // OCCT L276-277.
        let mut sli_list: Vec<Shape> = Vec::new();
        sli_list.push(first_face.clone());

        // OCCT L279-297.
        if sliding {
            sliding = false;
            let s = match brep_tool_surface(&first_face) {
                Some(mut s) => {
                    if let Surface3::Trimmed(t) = &s {
                        s = (*t.basis).clone();
                    }
                    s
                }
                None => Surface3::Plane(self.my_pln),
            };
            // if plane or cylinder : sliding is possible
            match &s {
                Surface3::Plane(_) | Surface3::Cylinder(_) => {
                    sliding = true;
                }
                _ => {}
            }
        }

        // OCCT L303-323 (architecture difference #3).
        if sliding {
            let p1 = self.rib_slot.my_first_pnt + self.my_dir;
            let mut pool = BRep::new();
            let ee1 = make_edge_p_p(&mut pool, self.rib_slot.my_first_pnt, p1);
            let ext1 = BRepExtremaExtCF::new(&ee1, &first_face);
            if ext1.nb_ext() == 1
                && ext1.square_distance(1)
                    <= brep_tool_tolerance(&first_face) * brep_tool_tolerance(&first_face)
            {
                let p2 = self.rib_slot.my_last_pnt + self.my_dir;
                let ee2 = make_edge_p_p(&mut pool, self.rib_slot.my_last_pnt, p2);
                let ext2 = BRepExtremaExtCF::new(&ee2, &last_face);
                sliding = ext2.nb_ext() == 1
                    && ext2.square_distance(1)
                        <= brep_tool_tolerance(&last_face) * brep_tool_tolerance(&last_face);
            } else {
                sliding = false;
            }
        }

        // OCCT L325-352 (architecture difference #3).
        if !gp_vec_is_equal(self.my_dir1, nulldir, CONFUSION, CONFUSION) {
            if sliding {
                let p1 = self.rib_slot.my_first_pnt + self.my_dir1;
                let mut pool = BRep::new();
                let ee1 = make_edge_p_p(&mut pool, self.rib_slot.my_first_pnt, p1);
                let ext1 = BRepExtremaExtCF::new(&ee1, &first_face);
                if ext1.nb_ext() == 1
                    && ext1.square_distance(1)
                        <= brep_tool_tolerance(&first_face) * brep_tool_tolerance(&first_face)
                {
                    let p2 = self.rib_slot.my_last_pnt + self.my_dir1;
                    let ee2 = make_edge_p_p(&mut pool, self.rib_slot.my_last_pnt, p2);
                    let ext2 = BRepExtremaExtCF::new(&ee2, &last_face);
                    sliding = ext2.nb_ext() == 1
                        && ext2.square_distance(1)
                            <= brep_tool_tolerance(&last_face)
                                * brep_tool_tolerance(&last_face);
                } else {
                    sliding = false;
                }
            }
        }

        // ---case of sliding : construction of the profile face

        // OCCT L360-410.
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
                self.rib_slot.my_status_error = BRepFeatStatusError::NoFaceProf;
                self.rib_slot.not_done();
                return;
            }
            // ---Propagation on faces of the initial shape
            // to find the faces concerned by the rib
            // OCCT L398-409.
            let mut falseside = true;
            sliding = self.propagate(
                &mut sli_list,
                &prof,
                self.rib_slot.my_first_pnt,
                self.rib_slot.my_last_pnt,
                &mut falseside,
            );
            if !falseside {
                self.rib_slot.my_status_error = BRepFeatStatusError::FalseSide;
                self.rib_slot.not_done();
                return;
            }
        }

        // ---Generation of the base of the rib profile

        // OCCT L414-421.
        let mut bb_pool = BRep::new();
        let w = make_wire_local(&mut bb_pool);
        let mut the_previous_edge = Shape::null();
        let mut the_fv = Shape::null();
        let mut counter = 1i32;

        // ---case of sliding

        // OCCT L424-750.
        if sliding && !self.my_list_of_edges.is_empty() {
            // OCCT L426: BRepTools_WireExplorer EX1(myWire) (architecture
            // difference #7).
            let wire_edges_list = wire_edges(&self.rib_slot.my_wire);
            let mut idx = 0usize;
            while idx < wire_edges_list.len() {
                let e = &wire_edges_list[idx];
                idx += 1;
                if !self.rib_slot.my_lfmap.contains_key(&shape_key(e)) {
                    self.rib_slot
                        .my_lfmap
                        .insert(shape_key(e), (e.clone(), Vec::new()));
                }
                if shape_is_same(e, &first_edge) {
                    // OCCT L437-454.
                    let Some((c_raw, f, l)) = brep_tool_curve(e) else {
                        break;
                    };
                    let mut cc = rcad_kernel::geom::Curve3::Trimmed(
                        rcad_kernel::geom::TrimmedCurve3::new(c_raw, f, l),
                    );
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
                    // OCCT L455-470.
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
                    ee1.orientation = e.orientation; // OCCT L471
                    // OCCT L474-483.
                    if counter == 1 {
                        the_fv = top_exp_first_vertex(&ee1, true);
                    }
                    self.rib_slot
                        .my_lfmap
                        .entry(shape_key(e))
                        .or_insert_with(|| (e.clone(), Vec::new()))
                        .1
                        .push(ee1.clone());
                    add_to_wire_local(&mut bb_pool, &w, &ee1);
                    the_previous_edge = ee1;
                    counter += 1;
                    idx += 1;
                    break;
                }
            }

            // Case of several edges

            // OCCT L488-561.
            if !shape_is_same(&first_edge, &last_edge) {
                while idx < wire_edges_list.len() {
                    let e = &wire_edges_list[idx];
                    idx += 1;
                    if !self.rib_slot.my_lfmap.contains_key(&shape_key(e)) {
                        self.rib_slot
                            .my_lfmap
                            .insert(shape_key(e), (e.clone(), Vec::new()));
                    }
                    the_list.push(e.clone());
                    if !shape_is_same(e, &last_edge) {
                        // OCCT L502-527.
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
                        e11.orientation = e.orientation; // OCCT L517
                        the_previous_edge = e11.clone();
                        self.rib_slot
                            .my_lfmap
                            .entry(shape_key(e))
                            .or_insert_with(|| (e.clone(), Vec::new()))
                            .1
                            .push(e11.clone());
                        add_to_wire_local(&mut bb_pool, &w, &e11);
                        if counter == 1 {
                            the_fv = top_exp_first_vertex(&e11, true);
                        }
                        counter += 1;
                    } else {
                        // OCCT L531-558.
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
                        ee.orientation = e.orientation; // OCCT L547
                        add_to_wire_local(&mut bb_pool, &w, &ee);
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

            // OCCT L563-695: the while(!FirstOK) pass.
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
                    sliding = false;
                    break;
                }
                let edg = self.my_list_of_edges[it_idx].clone();
                let Some((ccc, f, l)) = brep_tool_curve(&edg) else {
                    break;
                };
                let mut cc = rcad_kernel::geom::Curve3::Trimmed(
                    rcad_kernel::geom::TrimmedCurve3::new(ccc, f, l),
                );
                if edg.orientation == rcad_kernel::topods::Orientation::Reversed {
                    // OCCT L582: cc->Reverse().
                    let rev = crate::feat::brep_feat_rib_slot::geom_curve_reversed(&cc);
                    cc = rev;
                }
                let dom = cc.default_domain();
                let mut fp = cc.point_at(dom[0]);
                let mut lp = cc.point_at(dom[1]);
                let mut dist = fp.distance(the_last_pnt);
                if dist <= self.my_tol {
                    sens = 1;
                    last_ok = true;
                } else {
                    dist = lp.distance(the_last_pnt);
                    if dist <= self.my_tol {
                        sens = 2;
                        last_ok = true;
                        let rev = crate::feat::brep_feat_rib_slot::geom_curve_reversed(&cc);
                        cc = rev;
                        let rdom = cc.default_domain();
                        fp = cc.point_at(rdom[0]);
                        lp = cc.point_at(rdom[1]);
                    }
                }
                let mut first_flag = 0i32;
                if sens == 1 && lp.distance(self.rib_slot.my_first_pnt) <= self.my_tol {
                    first_ok = true;
                    first_flag = 1;
                } else if sens == 2 && fp.distance(self.rib_slot.my_first_pnt) <= self.my_tol {
                    first_ok = true;
                    first_flag = 2;
                }

                // OCCT L615-668.
                if last_ok {
                    let fpar = cc.default_domain()[0];
                    let lpar = cc.default_domain()[1];
                    let eeee = if !first_ok {
                        if the_previous_edge.is_null() {
                            let mut e_pool = BRep::new();
                            crate::feat::brep_feat_rib_slot::make_edge_cl(
                                &mut e_pool,
                                &cc,
                                fpar,
                                lpar,
                            )
                        } else {
                            let v1 = top_exp_last_vertex(&the_previous_edge, true);
                            // OCCT L630: BB.UpdateVertex(v1, dist).
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
                    add_to_wire_local(&mut bb_pool, &w, &eeee);
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

                // OCCT L670-677.
                if first_flag == 1 {
                    the_l_edge = edg.clone();
                } else if first_flag == 2 {
                    the_l_edge = the_edge.clone();
                }
                let _ = &the_l_edge;

                // OCCT L679-694: the list Remove + re-Initialize.
                if last_ok {
                    if it_idx < self.my_list_of_edges.len() {
                        self.my_list_of_edges.remove(it_idx);
                    }
                    it_idx = 0;
                    last_ok = false;
                } else {
                    it_idx += 1;
                }
                sens = 0;
            }
            let _ = (&first_ok, &last_ok, &sens);

            // OCCT L697-749: the SlidMap re-bind.
            let mut slid_map: HashMap<(u64, u32), (Shape, Vec<Shape>)> = HashMap::new();
            if sliding && counter1 > counter {
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

            // OCCT L748-749.
            self.my_slface = slid_map;
        }

        // ---Arguments of LocOpe_LinearForm : arguments of the prism sliding

        // OCCT L753-761.
        if sliding {
            let mut f_pool = BRep::new();
            let f = make_face_of_plane_local(&mut f_pool, &self.my_pln);
            let closed = brep_tool_is_closed(&w);
            let _ = closed;
            add_face_wire_local(&mut f_pool, &f, &w);
            self.rib_slot.my_skface = f.clone();
            self.rib_slot.my_pbase = self.rib_slot.my_skface.clone();
            self.rib_slot.my_suntil = Shape::null();
        }
        let _ = brep_tool_is_closed;

        // ---Case without sliding : construction of the profile face

        // OCCT L764-846.
        if !sliding {
            // OCCT L778-780.
            let explo1 = explorer(&bnd_face, ShapeType::Wire, ShapeType::Shape);
            let www2 = explo1.first().cloned().unwrap_or_else(Shape::null);
            let _ = www2;
            // OCCT L781-841.
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
                self.rib_slot.my_status_error = BRepFeatStatusError::NoFaceProf;
                self.rib_slot.not_done();
                return;
            }
            // OCCT L820-832.
            let mut falseside = true;
            self.propagate(
                &mut sli_list,
                &prof,
                self.rib_slot.my_first_pnt,
                self.rib_slot.my_last_pnt,
                &mut falseside,
            );
            if !falseside {
                self.rib_slot.my_status_error = BRepFeatStatusError::FalseSide;
                self.rib_slot.not_done();
                return;
            }

            // OCCT L834-845.
            self.my_slface.clear();
            let mut pool_c = BRep::new();
            let comp = make_shell_local(&mut pool_c);
            for it in &sli_list {
                add_to_shell_local(&mut pool_c, &comp, it);
            }
            self.rib_slot.my_suntil = comp.clone();
            self.rib_slot.my_skface = prof.clone();
            self.rib_slot.my_pbase = prof;
        }

        // OCCT L848.
        self.rib_slot.my_sliding = sliding;

        // OCCT L850-856.
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

    /// OCCT BRepFeat_MakeLinearForm::Add (cxx L858-902).
    pub fn add(&mut self, the_e: &Shape, the_f: &Shape) {
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

    /// OCCT BRepFeat_MakeLinearForm::Perform (cxx L905-1032) — construction
    /// of rib from a profile and the initial shape.
    pub fn perform(&mut self) {
        // OCCT L911-919.
        if self.rib_slot.my_sbase.is_null()
            || self.rib_slot.my_skface.is_null()
            || self.rib_slot.my_pbase.is_null()
        {
            self.rib_slot.my_status_error = BRepFeatStatusError::NotInitialized;
            self.rib_slot.not_done();
            return;
        }

        // OCCT L921-929.
        let nulldir = DVec3::ZERO;
        let length = self.my_dir.length() + self.my_dir1.length();
        self.rib_slot.my_glued_f.clear();
        if !self.rib_slot.my_suntil.is_null() {
            self.rib_slot.my_perf_selection = BRepFeatPerfSelection::SelectionU;
        } else {
            self.rib_slot.my_perf_selection = BRepFeatPerfSelection::NoSelection;
        }

        // OCCT L931-932.
        let dir = self.my_dir.normalize_or_zero();
        let v = dir * length;

        // OCCT L934-946.
        let mut the_form = LocOpeLinearForm::new();
        if gp_vec_is_equal(self.my_dir1, nulldir, CONFUSION, CONFUSION) {
            the_form.perform(
                &self.rib_slot.my_pbase,
                v,
                self.rib_slot.my_first_pnt,
                self.rib_slot.my_last_pnt,
            );
        } else {
            the_form.perform_trans(
                &self.rib_slot.my_pbase,
                v,
                self.my_dir1,
                self.rib_slot.my_first_pnt,
                self.rib_slot.my_last_pnt,
            );
        }

        // OCCT L948.
        let vrai_form = the_form.shape(); // primitive of the rib

        // OCCT L950-952.
        self.rib_slot.my_faces_for_draft.push(the_form.first_shape());
        self.rib_slot.my_faces_for_draft.push(the_form.last_shape());
        // OCCT L952: MajMap — management of descendants.
        maj_map(
            &self.rib_slot.my_pbase,
            &the_form,
            &mut self.rib_slot.my_map,
            &mut self.rib_slot.my_f_shape,
            &mut self.rib_slot.my_l_shape,
        );

        // OCCT L954-968.
        for exx in explorer(&self.rib_slot.my_pbase, ShapeType::Edge, ShapeType::Shape) {
            if !self.rib_slot.my_map.contains_key(&shape_key(&exx)) {
                self.rib_slot.my_status_error = BRepFeatStatusError::IncSlidFace;
                self.rib_slot.not_done();
                return;
            }
        }

        // OCCT L970.
        self.rib_slot.my_gshape = vrai_form;
        // OCCT L971: SetGluedFaces — management of sliding faces.
        set_glued_faces(&self.my_slface, &the_form, &mut self.rib_slot.my_glued_f);

        // OCCT L973-986.
        if !self.rib_slot.my_glued_f.is_empty() && !self.rib_slot.my_suntil.is_null() {
            self.rib_slot.my_status_error = BRepFeatStatusError::InvShape;
            self.rib_slot.not_done();
            return;
        }

        // OCCT L988: LFPerform().
        self.rib_slot.lf_perform();

        // OCCT L991-1029: the commented DBRep debug block is not translated.
    }

    /// OCCT BRepFeat_MakeLinearForm::Propagate (cxx L1035-1288) —
    /// propagation on faces of the initial shape, find faces concerned by
    /// the rib.
    pub fn propagate(
        &mut self,
        sli_list: &mut Vec<Shape>,
        fac: &Shape,
        firstpnt: DVec3,
        lastpnt: DVec3,
        falseside: &mut bool,
    ) -> bool {
        // OCCT L1049-1050.
        let mut firstpoint = firstpnt;
        let mut lastpoint = lastpnt;

        // OCCT L1052-1058.
        let mut result = true;
        let mut current_face = sli_list.first().cloned().unwrap_or_else(Shape::null);
        let _save_face = current_face.clone();

        // OCCT L1060-1064.
        let mut last_ok = false;
        let mut first_ok = false;
        let mut v1_ok = false;
        let mut v2_ok = false;
        let mut v1 = Shape::null();
        let mut v2 = Shape::null();

        // OCCT L1066-1068.
        let mut sect = SectionOp::from_shapes(fac.clone(), current_face.clone());
        sect.approximation(true);
        sect.build();

        // OCCT L1070-1101.
        let mut eb = Shape::null();
        let sect_shape = sect.algo.bs.result.clone().unwrap_or_else(Shape::null);
        for ex in explorer(&sect_shape, ShapeType::Edge, ShapeType::Shape) {
            let ec = ex;
            v1 = top_exp_first_vertex(&ec, true);
            v2 = top_exp_last_vertex(&ec, true);
            let p1 = brep_tool_pnt(&v1);
            let p2 = brep_tool_pnt(&v2);
            let t1 = brep_tool_tolerance(&v1);
            let t2 = brep_tool_tolerance(&v2);
            let c1f = p1.distance(firstpoint) <= t1;
            let c2f = p2.distance(firstpoint) <= t2;
            let c1l = p1.distance(lastpoint) <= t1;
            let c2l = p2.distance(lastpoint) <= t2;
            if c1f || c2f || c1l || c2l {
                eb = ec;
                if c1f || c1l {
                    v1_ok = true;
                }
                if c2f || c2l {
                    v2_ok = true;
                }
                if c1f || c2f {
                    first_ok = true;
                }
                if c1l || c2l {
                    last_ok = true;
                }
                break;
            }
        }

        // OCCT L1103-1107.
        if eb.is_null() {
            *falseside = false;
            return false;
        }
        // OCCT L1108-1115.
        self.my_slface
            .insert(shape_key(&current_face), (current_face.clone(), Vec::new()));
        self.my_slface
            .get_mut(&shape_key(&current_face))
            .expect("mySlface(CurrentFace)")
            .1
            .push(eb.clone());
        self.my_list_of_edges.clear();
        self.my_list_of_edges.push(eb.clone());

        // two points are on the same face.
        // OCCT L1118-1121.
        if last_ok && first_ok {
            return result;
        }

        // OCCT L1123-1131.
        let mut mapedges: indexmap::IndexMap<(u64, u32), (Shape, Vec<Shape>)> =
            indexmap::IndexMap::new();
        crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors(
            &self.rib_slot.my_sbase,
            ShapeType::Edge,
            ShapeType::Face,
            &mut mapedges,
        );
        let mut first_edge = Shape::null();

        let mut v_previous = Shape::null();
        let mut ptprev = DVec3::ZERO;

        // OCCT L1136-1283: the while(!(LastOK && FirstOK)) pass.
        while !(last_ok && first_ok) {
            // OCCT L1137-1146.
            if v1_ok {
                v_previous = v2.clone();
                ptprev = brep_tool_pnt(&v2);
            } else {
                v_previous = v1.clone();
                ptprev = brep_tool_pnt(&v1);
            }

            // find edge connected to v1 or v2:
            // OCCT L1149-1204.
            for ex in explorer(&current_face, ShapeType::Edge, ShapeType::Shape) {
                let rfe = ex;
                // OCCT L1152: BRepExtrema_ExtPC projF(Vprevious, rfe) — the
                // wrapper ctor is Initialize(E) + Perform(V)
                // (BRepExtrema_ExtPC.cxx L30-33); a non-geometric edge leaves
                // the wrapper not-done (the rcad continue skips the same
                // block the OCCT !IsDone() gate skips).
                let Some((c, f, l)) = brep_tool_curve(&rfe) else {
                    continue;
                };
                // OCCT BRepExtrema_ExtPC::Initialize (cxx L35-47): myHC = new
                // BRepAdaptor_Curve(E); Tol = min(BRep_Tool::Tolerance(E),
                // Precision::Confusion()); Tol = max(myHC->Resolution(Tol),
                // Precision::PConfusion()); BRep_Tool::Range(E, U1, U2).
                let a_adaptor = GeomCurveAdaptor::new(c.clone());
                let mut a_tol = brep_tool_tolerance(&rfe).min(CONFUSION);
                a_tol = a_adaptor.resolution(a_tol).max(PCONFUSION);
                let a_tool = CurveToolHandle::for_curve3(&c, &a_adaptor, &a_adaptor);
                // OCCT cxx L50-56: BRep_Tool::Pnt(V); myExtPC.Perform(P).
                let proj_f = ExtremaExtPC::new_point_curve_ranged(
                    brep_tool_pnt(&v_previous),
                    &a_tool,
                    f,
                    l,
                    a_tol,
                );
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
                        && dist2min <= brep_tool_tolerance(&rfe) * brep_tool_tolerance(&rfe)
                    {
                        first_edge = rfe.clone();
                        // If the edge is not perpendicular to the plane of
                        // the rib it is required to set Sliding(result) to
                        // false. (architecture difference #3: the
                        // BRepExtrema_ExtPF carrier reports not-done, so the
                        // OCCT guard keeps result = false.)
                        if result {
                            result = false;
                            let ve1 = top_exp_first_vertex(&rfe, true);
                            let ve2 = top_exp_last_vertex(&rfe, true);
                            let mut perp = BRepExtremaExtPFCarrier::new(&ve1, fac);
                            if perp.is_done() {
                                let pe1 = perp.point(1);
                                perp.perform(&ve2, fac);
                                if perp.is_done() {
                                    let pe2 = perp.point(1);
                                    if pe1.distance(pe2) <= brep_tool_tolerance(&rfe) {
                                        result = true;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }

            // OCCT L1206-1217.
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

            // OCCT L1219-1243.
            let mut sectf = SectionOp::from_shapes(fac.clone(), current_face.clone());
            sectf.approximation(true);
            sectf.build();
            let mut edg1 = Shape::null();
            let sectf_shape = sectf.algo.bs.result.clone().unwrap_or_else(Shape::null);
            for ex in explorer(&sectf_shape, ShapeType::Edge, ShapeType::Shape) {
                edg1 = ex;
                v1 = top_exp_first_vertex(&edg1, true);
                v2 = top_exp_last_vertex(&edg1, true);
                let t1 = brep_tool_tolerance(&v1);
                let t2 = brep_tool_tolerance(&v2);
                let p1 = brep_tool_pnt(&v1);
                let p2 = brep_tool_pnt(&v2);
                v1_ok = p1.distance(ptprev) <= t1;
                v2_ok = p2.distance(ptprev) <= t2;
                if v1_ok || v2_ok {
                    break;
                }
            }

            // OCCT L1245-1281.
            if v1_ok {
                if !first_ok {
                    let dp = brep_tool_pnt(&v2).distance(firstpoint);
                    if dp <= 2.0 * brep_tool_tolerance(&v2) {
                        first_ok = true;
                        // OCCT: BB.UpdateVertex(v2, dp).
                    }
                }
                if !last_ok {
                    let dp = brep_tool_pnt(&v2).distance(lastpoint);
                    if dp <= 2.0 * brep_tool_tolerance(&v2) {
                        last_ok = true;
                        // OCCT: BB.UpdateVertex(v2, dp).
                    }
                }
            } else if v2_ok {
                if !first_ok {
                    let dp = brep_tool_pnt(&v1).distance(firstpoint);
                    if dp <= 2.0 * brep_tool_tolerance(&v1) {
                        first_ok = true;
                        // OCCT: BB.UpdateVertex(v1, dp).
                    }
                }
                if !last_ok {
                    let dp = brep_tool_pnt(&v1).distance(lastpoint);
                    if dp <= 2.0 * brep_tool_tolerance(&v1) {
                        last_ok = true;
                        // OCCT: BB.UpdateVertex(v1, dp).
                    }
                }
            } else {
                // end by chaining the section
                return false;
            }

            // OCCT L1282-1286 (hmm: L1275-1281 in the OCCT source — the
            // mySlface/myListOfEdges appends).
            self.my_slface
                .insert(shape_key(&current_face), (current_face.clone(), Vec::new()));
            self.my_slface
                .get_mut(&shape_key(&current_face))
                .expect("mySlface(CurrentFace)")
                .1
                .push(edg1.clone());
            self.my_list_of_edges.push(edg1);
        }

        // OCCT L1285-1287.
        result
    }
}
/// OCCT static MajMap (cxx L1292-1335) — the LinearForm variant.
fn maj_map(
    the_b: &Shape,
    the_p: &LocOpeLinearForm,
    the_map: &mut HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_f_shape: &mut Shape,
    the_l_shape: &mut Shape,
) {
    // OCCT L1305-1314.
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
    // OCCT L1316-1325.
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
    // OCCT L1327-1334.
    for exp in explorer(the_b, ShapeType::Edge, ShapeType::Shape) {
        if !the_map.contains_key(&shape_key(&exp)) {
            let shapes = the_p.shapes(&exp).clone();
            the_map.insert(shape_key(&exp), (exp.clone(), shapes));
        }
    }
}

/// OCCT static SetGluedFaces (cxx L1337-1366).
fn set_glued_faces(
    the_slmap: &HashMap<(u64, u32), (Shape, Vec<Shape>)>,
    the_prism: &LocOpeLinearForm,
    the_map: &mut HashMap<(u64, u32), (Shape, Shape)>,
) {
    // Slidings
    // OCCT L1345-1365.
    if !the_slmap.is_empty() {
        for (_key, (fac, ledg)) in the_slmap.iter() {
            let fac = fac;
            for it in ledg {
                let gfac = the_prism.shapes(it);
                if gfac.len() != 1 {
                    // OCCT L1356: the debug print only.
                }
                if let Some(g_first) = gfac.first() {
                    the_map.insert(shape_key(g_first), (g_first.clone(), fac.clone()));
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
