//! OCCT ChFi3d_Builder_2.cxx — 1:1 translation (Stage 1f).
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKFillet/ChFi3d/
//!         ChFi3d_Builder_2.cxx (3900 lines).
//!
//! Coverage of this file: OCCT L74-1627 (anonymous-namespace helpers,
//! CallPerformSurf, StripeOrientations, ConexFaces, StartSol overloads,
//! ChFi3d_BuildPlane, SearchFace).  The rest of the .cxx is split by size:
//! `chfi3d_builder_2b` (PerformSetOfSurfOnElSpine, ChFi3d_BoxDiag,
//! PerformSetOfKGen) and `chfi3d_builder_2c` (ChFi3d_SingularExtremity,
//! IsFree, ChFi3d_MakeExtremities, ChFi3d_Purge, InsertAfter, RemoveSD,
//! InsertBefore).
//!
//! Pre-existing translations (name collisions on `impl ChFi3dBuilder`,
//! they live in `chfi3d.rs` and are NOT redefined here):
//!   - StripeOrientations (OCCT L826-855)  -> chfi3d.rs stripe_orientations
//!   - ConexFaces         (OCCT L859-889)  -> chfi3d.rs conex_faces
//!   - PerformSetOfKPart  (OCCT L3004-3280)-> chfi3d.rs perform_set_of_k_part
//!   - PerformSetOfSurf   (OCCT L3882-3900)-> chfi3d.rs perform_set_of_surf
//!     (the chfi3d.rs PerformSetOfSurf keeps its PerformSetOfKGen /
//!     ChFi3d_MakeExtremities calls pending; rewiring onto this file's
//!     full translations needs an edit of chfi3d.rs and is left to the
//!     acceptance pass).
//!
//! rcad architecture notes (chfi3d_builder_0 convention): the owning BRep
//! is passed to the free functions as an extra first argument where OCCT
//! reads geometry from the TopoDS handle graph.

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema::{extrema_locate_ext_pc, ExtPC, ExtPS};
use rcad_kernel::core::precision::CONFUSION;
use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::topo::topods::BRepTool as _;
use rcad_kernel::topo::topods::{Orientation, Shape};
use rcad_kernel::topods;

use super::chfi3d::{is_tangent_faces, next_side, topabs_reverse};
use super::brep_blend_func_consrad::BlendFuncConstRad;
use super::brep_blend_func_evolrad::BlendFuncEvolRad;
use super::brep_blend_walking::BRepBlendWalking;
use super::chfi3d_builder_0::{
    brep_tool_parameter, chfi3d_compute_curves, topexp_face_edges, topexp_vertices,
    vec_is_parallel, BRepAdaptorSurface, GeomAdaptorSurface, P_CONFUSION,
};
use super::chfi3d_ds::TopOpeBRepDSHDataStructure;
use super::chfi_ds::{
    ChFiDSElSpine, ChFiDS_ErrorStatus, ChFiDSMap, ChFiDS_CommonPoint, ChFiDSSpineHandle,
    ChFiDSSurfData, SharedStripe, SharedSurfData,
};
use crate::geomalgo::gtests_stubs::GeomAbsShape;

// =========================================================================
// OCCT BRepAdaptor_Curve2d (TKTopAlgo) — the 2d curve of an edge in a face
// ("representation of the obstacle" in StartSol).  rcad form carrier: the
// pcurve is resolved through the owning BRep on each Value query.
// =========================================================================
#[derive(Debug, Clone)]
pub struct BRepAdaptorCurve2d {
    pub brep: topods::BRep,
    pub edge: Shape,
    pub face: Shape,
}

impl Default for BRepAdaptorCurve2d {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepAdaptorCurve2d {
    /// OCCT BRepAdaptor_Curve2d default construction (null handles).
    pub fn new() -> Self {
        BRepAdaptorCurve2d {
            brep: topods::BRep::default(),
            edge: Shape::null(),
            face: Shape::null(),
        }
    }

    /// OCCT BRepAdaptor_Curve2d::Initialize(E, F).
    pub fn initialize(&mut self, brep: &topods::BRep, e: &Shape, f: &Shape) {
        self.brep = brep.clone();
        self.edge = e.clone();
        self.face = f.clone();
    }

    /// OCCT BRepAdaptor_Curve2d::Value(U).
    pub fn value(&self, u: f64) -> DVec2 {
        self.brep
            .curve_on_surface(&self.edge, &self.face)
            .map(|(pc, _, _)| pc.point_at(u))
            .unwrap_or(DVec2::ZERO)
    }

    /// OCCT BRepAdaptor_Curve2d::Edge().
    pub fn edge(&self) -> &Shape {
        &self.edge
    }

    /// OCCT occ::handle nullity (an empty carrier encodes the null handle).
    pub fn is_null(&self) -> bool {
        self.edge.is_null()
    }
}

// =========================================================================
// OCCT BRepTopAdaptor_TopolTool / Adaptor3d_TopolTool (TKTopAlgo) — pending
// translation.  The carrier keeps the loaded surface so the
// Initialize/Classify call sites keep their form; Classify falls back to
// the surface UV box until the real wire classifier lands (Adaptor3d_
// TopolTool::Classify over the face restrictions).
// =========================================================================

/// OCCT TopAbs_State (TopAbs_ShapeState enumeration values used here).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopAbsState {
    In,
    Out,
    On,
    Unknown,
}

#[derive(Debug, Clone, Default)]
pub struct BRepTopAdaptorTopolTool {
    pub surface: Option<BRepAdaptorSurface>,
}

impl BRepTopAdaptorTopolTool {
    /// OCCT Adaptor3d_TopolTool::Initialize(theSurface) (the
    /// BRepTopAdaptor_TopolTool(BRepAdaptor_Surface) overload).
    pub fn initialize_brep(&mut self, hs: &BRepAdaptorSurface) {
        self.surface = Some(hs.clone());
    }

    /// OCCT Adaptor3d_TopolTool::Classify(P, Tol, WithBound) — pending the
    /// real classifier; the UV-box fallback reports IN inside the loaded
    /// adaptor bounds and OUT outside.
    pub fn classify(&self, p2d: DVec2, _tol: f64, _with_bound: bool) -> TopAbsState {
        match &self.surface {
            Some(s) => {
                if p2d.x >= s.ufirst && p2d.x <= s.ulast && p2d.y >= s.vfirst && p2d.y <= s.vlast
                {
                    TopAbsState::In
                } else {
                    TopAbsState::Out
                }
            }
            None => TopAbsState::Unknown,
        }
    }
}

impl BRepAdaptorSurface {
    /// OCCT BRepAdaptor_Surface default construction (`new
    /// BRepAdaptor_Surface()` — the empty adaptor before Initialize).  rcad
    /// carries a null face until Initialize.
    pub fn empty() -> Self {
        BRepAdaptorSurface {
            brep: topods::BRep::default(),
            face: Shape::null(),
            surface: Surface3::Plane(rcad_kernel::geom::Plane {
                origin: DVec3::ZERO,
                normal: DVec3::Z,
                u_dir: DVec3::X,
                v_dir: DVec3::Y,
            }),
            ufirst: 0.0,
            ulast: 0.0,
            vfirst: 0.0,
            vlast: 0.0,
            bounds_set: false,
        }
    }

    /// OCCT BRepAdaptor_Surface::Initialize(F) (the in-place overload).
    pub fn initialize_face(&mut self, brep: &topods::BRep, f: &Shape) {
        *self = BRepAdaptorSurface::initialize(brep, f);
    }
}

// =========================================================================
// OCCT ChFiDS_ElSpine — accessors and the curve queries needed by Builder_2
// (the ElSpine composite curve machinery / BRepAdaptor_Curve myCurve is the
// pending boundary; the linear extension model through the stored endpoint
// tangents stands in, see the stage report).
// =========================================================================

impl ChFiDSElSpine {
    // The ChFiDS_ElSpine method translations (SetPeriodic, Period,
    // SaveFirst/LastParameter, SetOrigin, Resolution, SetCurve, the
    // curve-backed Value/D1/D2 and the vertex list) live in chfi_ds_spine
    // per the ChFiDS_ElSpine.hxx field completion — not redefined here.
}

// =========================================================================
// OCCT BRepAdaptor_Curve stand-ins over a spine edge (orientation-adjusted
// Value/D1, the BRepAdaptor_Curve conventions for REVERSED edges).
// =========================================================================

/// OCCT BRepAdaptor_Curve::Value(U) over a topological edge.
pub(crate) fn brep_adaptor_curve_value(e: &Shape, u: f64) -> DVec3 {
    let ed = e.as_edge().expect("not an edge");
    let c = ed.curve.as_ref().expect("edge curve");
    if e.orientation == Orientation::Reversed {
        c.point_at(ed.range[0] + ed.range[1] - u)
    } else {
        c.point_at(u)
    }
}

/// OCCT BRepAdaptor_Curve::D1(U, P, V) over a topological edge.
pub(crate) fn brep_adaptor_curve_d1(e: &Shape, u: f64) -> (DVec3, DVec3) {
    let ed = e.as_edge().expect("not an edge");
    let c = ed.curve.as_ref().expect("edge curve");
    if e.orientation == Orientation::Reversed {
        let uu = ed.range[0] + ed.range[1] - u;
        (c.point_at(uu), -c.derivative_at(uu))
    } else {
        (c.point_at(u), c.derivative_at(u))
    }
}

/// OCCT BRep_Tool::Parameters(V, F) — the UV of a vertex on a face surface
/// (GeomAPI_ProjectPointOnSurf over the natural domain).
pub(crate) fn brep_tool_parameters_vf(brep: &topods::BRep, v: &Shape, f: &Shape) -> DVec2 {
    let Some(surf) = f.as_face().and_then(|fd| fd.surface.clone()) else {
        return DVec2::ZERO;
    };
    let Some(pv) = brep
        .tshapes
        .get(v.index)
        .and_then(|ts| match ts.as_ref() {
            topods::TShape::Vertex(vd) => Some(vd.point),
            _ => None,
        })
    else {
        return DVec2::ZERO;
    };
    let ext = ExtPS::new(pv, &surf, 1.0e-6, 1.0e-6);
    if ext.nb_ext() >= 1 {
        let p = ext.point(1);
        DVec2::new(p.u, p.v)
    } else {
        DVec2::ZERO
    }
}

/// OCCT gp_Vec::AngleWithRef(Other, Dir) — the signed angle in [-PI, PI],
/// positive when (this ^ Other) agrees with Dir.
pub(crate) fn vec_angle_with_ref(a: DVec3, b: DVec3, dir: DVec3) -> f64 {
    let tri_prod = a.cross(b);
    let sin_tri_prod = tri_prod.length();
    let cos_ang = a.dot(b);
    let mut ang = sin_tri_prod.atan2(cos_ang);
    let dir_dot_tri_prod = dir.dot(tri_prod);
    if dir_dot_tri_prod < 0.0 {
        ang = -ang;
    }
    ang
}

/// OCCT gp_Trsf::SetTransformation(gp_Ax3) + Invert — the global->local
/// (resp. local->global after Invert) change of frame for a point.
pub(crate) fn trsf_global_to_local(origin: DVec3, xdir: DVec3, ydir: DVec3, zdir: DVec3, p: DVec3) -> DVec3 {
    let d = p - origin;
    DVec3::new(d.dot(xdir), d.dot(ydir), d.dot(zdir))
}

pub(crate) fn trsf_local_to_global(origin: DVec3, xdir: DVec3, ydir: DVec3, zdir: DVec3, c: DVec3) -> DVec3 {
    origin + xdir * c.x + ydir * c.y + zdir * c.z
}

/// OCCT BRepTools::IsReallyClosed(E, F) — the edge occurs twice among the
/// face's wires (local copy; the builder_0 twin is private).
pub(crate) fn brep_tools_is_really_closed(brep: &topods::BRep, e: &Shape, f: &Shape) -> bool {
    let mut n = 0usize;
    for we in topexp_face_edges(brep, f) {
        if we.is_same(e) {
            n += 1;
        }
    }
    n > 1
}

/// OCCT anonymous-namespace getCurveOnSurface (ChFi3d_Builder_2.cxx
/// L77-85) — the 2d curve of an edge on a face without the parameters.
pub(crate) fn get_curve_on_surface(
    brep: &topods::BRep,
    the_edge: &Shape,
    the_face: &Shape,
) -> Option<rcad_kernel::geom::Curve2d> {
    let (pc, _first_param, _last_param) = brep.curve_on_surface(the_edge, the_face)?;
    Some(pc)
}

/// OCCT ChFi3d_Builder_2.cxx L87-188 — ChFi3d_CoupeParPlan.
#[allow(clippy::too_many_arguments)]
pub(crate) fn chfi3d_coupe_par_plan(
    _brep: &topods::BRep,
    compoint1: &ChFiDS_CommonPoint,
    compoint2: &ChFiDS_CommonPoint,
    hconge: &GeomAdaptorSurface,
    uv1: DVec2,
    uv2: DVec2,
    tol3d: f64,
    tol2d: f64,
    c3d: &mut Option<rcad_kernel::geom::Curve3>,
    pcurve: &mut Option<rcad_kernel::geom::Curve2d>,
    tolreached: &mut f64,
    pardeb: &mut f64,
    parfin: &mut f64,
    plane: &mut bool,
) {
    *plane = true;
    if compoint1.is_on_arc() && compoint2.is_on_arc() {
        let bcurv1 = compoint1.arc().clone();
        let bcurv2 = compoint2.arc().clone();
        let par_e1 = compoint1.parameter_on_arc();
        let par_e2 = compoint2.parameter_on_arc();
        let (p1, t1) = brep_adaptor_curve_d1(&bcurv1, par_e1);
        let (p2, t2) = brep_adaptor_curve_d1(&bcurv2, par_e2);
        // OCCT: gp_Dir tgt1(t1) / tgt2(t2) — normalized on construction.
        let tgt1 = t1.normalize();
        let tgt2 = t2.normalize();
        let v12 = p2 - p1;
        let d12 = v12.normalize();
        let nor = tgt1.cross(d12);
        let plan = Surface3::Plane(rcad_kernel::geom::Plane::new(p1, nor));
        let scal = nor.dot(tgt2).abs();
        if scal < 0.01 {
            let hplan = GeomAdaptorSurface::new(plan.clone());
            let as_surf = plan.clone();
            let mut an_ext_ps = ExtPS::with_domain(
                p1,
                &as_surf,
                hplan.first_u_parameter(),
                hplan.last_u_parameter(),
                hplan.first_v_parameter(),
                hplan.last_v_parameter(),
                1.0e-3,
                1.0e-3,
            );
            let mut u1;
            let mut v1;
            {
                let p = an_ext_ps.point(1);
                u1 = p.u;
                v1 = p.v;
            }
            let mut pdeb = [0.0f64; 4];
            let mut pfin = [0.0f64; 4];
            pdeb[0] = uv1.x;
            pdeb[1] = uv1.y;
            pdeb[2] = u1;
            pdeb[3] = v1;
            an_ext_ps.perform(
                p2,
                &as_surf,
                hplan.first_u_parameter(),
                hplan.last_u_parameter(),
                hplan.first_v_parameter(),
                hplan.last_v_parameter(),
                1.0e-3,
                1.0e-3,
            );
            {
                let p = an_ext_ps.point(1);
                u1 = p.u;
                v1 = p.v;
            }
            pfin[0] = uv2.x;
            pfin[1] = uv2.y;
            pfin[2] = u1;
            pfin[3] = v1;
            // OCCT: C2dint2 (the pcurve on the cut plane) is not consumed.
            if let Some(cc) = chfi3d_compute_curves(hconge, &hplan, pdeb, pfin, tol3d, tol2d, tolreached)
            {
                *c3d = Some(cc.c3d.clone());
                *pcurve = Some(cc.pc1);
                // OCCT L171-172: Pardeb = C3d->FirstParameter().
                (*pardeb, *parfin) = curve_first_last(&cc.c3d);
            } else {
                *plane = false;
            }
        } else {
            *plane = false;
        }
    } else {
        *plane = false;
    }
}

/// OCCT Geom_Curve::FirstParameter()/LastParameter() over the rcad curve
/// kinds produced by ChFi3d_ComputeCurves (trimmed curves carry their range).
pub(crate) fn curve_first_last(c: &rcad_kernel::geom::Curve3) -> (f64, f64) {
    match c {
        rcad_kernel::geom::Curve3::Trimmed(t) => (t.first, t.last),
        other => {
            let d = other.default_domain();
            (d[0], d[1])
        }
    }
}

/// OCCT ChFi3d_Builder_2.cxx L190-209 — isTangentToArc.
pub(crate) fn is_tangent_to_arc(a_common_point: &ChFiDS_CommonPoint, an_angular_tolerance: f64) -> bool {
    if !a_common_point.has_vector() {
        return false;
    }
    let (_a_dummy_point, an_arc_tangent) =
        brep_adaptor_curve_d1(a_common_point.arc(), a_common_point.parameter_on_arc());
    let a_common_point_vector = a_common_point.vector;
    vec_is_parallel(a_common_point_vector, an_arc_tangent, an_angular_tolerance)
}

/// OCCT ChFi3d_Builder_2.cxx L718-765 (Builder_0.cxx) — ChFi3d_InterPlaneEdge.
/// The IntCurveSurface_HInter curve/plane intersection is pending
/// (TKGeomAlgo); the OCCT not-done path (no intersection point) is reported.
pub(crate) fn chfi3d_inter_plane_edge(
    _plan: &GeomAdaptorSurface,
    _c: &Shape,
    _w: &mut f64,
    _sens: bool,
    _tolc: f64,
) -> bool {
    // Pending: IntCurveSurface_HInter::Perform(C, Plan) translation.
    false
}

/// OCCT ChFi3d_Builder_2.cxx L211-316 — BonVoisin.
#[allow(clippy::too_many_arguments)]
pub(crate) fn bon_voisin(
    brep: &topods::BRep,
    point: DVec3,
    hs: &mut BRepAdaptorSurface,
    f: &mut Shape,
    plane: &GeomAdaptorSurface,
    cured: &Shape,
    sol_dep: &mut [f64; 4],
    x_dep: usize,
    y_dep: usize,
    ef_map: &ChFiDSMap,
    tol3d: f64,
) -> bool {
    let mut bonvoisin = true;
    let mut winter = 0.0f64;
    let papp = hs.surface.point_at(sol_dep[x_dep], sol_dep[y_dep]);
    let mut dist = f64::INFINITY; // OCCT: RealLast()
    let mut pc: Option<rcad_kernel::geom::Curve2d> = None;
    let mut found = false;

    for ecur in topexp_face_edges(brep, f) {
        if !ecur.is_same(cured) {
            let tolc = {
                let ed = ecur.as_edge().expect("not an edge");
                // OCCT: hc->Resolution(tol3d) over BRepAdaptor_Curve.
                let Some(c) = &ed.curve else { continue };
                curve_resolution(c, tol3d)
            };
            if chfi3d_inter_plane_edge(plane, &ecur, &mut winter, true, tolc) {
                let np = brep_adaptor_curve_value(&ecur, winter);
                let ndist = np.distance_squared(papp);
                if ndist < dist {
                    let mut ff = Shape::null();
                    let isclosed = brep.is_edge_closed_on_face(&ecur, f);
                    let isreallyclosed = brep_tools_is_really_closed(brep, &ecur, f);
                    for it in ef_map.find(&ecur).clone() {
                        ff = it;
                        let issame = ff.is_same(f);
                        let istg = is_tangent_faces(brep, &ecur, &ff, f, GeomAbsShape::G1);
                        if (!issame || (issame && isreallyclosed)) && istg {
                            found = true;
                            let mut newe = ecur.clone();
                            newe.orientation = Orientation::Forward;
                            dist = ndist;
                            hs.initialize_face(brep, &ff);
                            if isclosed && !isreallyclosed {
                                let mut fff = ff.clone();
                                fff.orientation = Orientation::Forward;
                                for ex2 in topexp_face_edges(brep, &fff) {
                                    if newe.is_same(&ex2) {
                                        newe = ex2.clone();
                                        pc = get_curve_on_surface(brep, &newe, &fff);
                                        break;
                                    }
                                }
                            } else {
                                pc = get_curve_on_surface(brep, &newe, &ff);
                            }
                            if let Some(pc) = &pc {
                                let coord = pc.point_at(winter);
                                sol_dep[x_dep] = coord.x;
                                sol_dep[y_dep] = coord.y;
                            }
                            if issame {
                                let (spt, sdu, sdv) =
                                    hs.surface.derivatives(sol_dep[x_dep], sol_dep[y_dep]);
                                let nors = sdu.cross(sdv);
                                let (cpt, cd) = brep_adaptor_curve_d1(&ecur, winter);
                                let _ = spt;
                                let vref = cpt - point;
                                let mut fff = ff.clone();
                                fff.orientation = Orientation::Forward;
                                if vref.dot(nors.cross(cd)) < 0.0 {
                                    newe.orientation = Orientation::Reversed;
                                }
                                pc = get_curve_on_surface(brep, &newe, &fff);
                                if let Some(pc) = &pc {
                                    let coord = pc.point_at(winter);
                                    sol_dep[x_dep] = coord.x;
                                    sol_dep[y_dep] = coord.y;
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    if !found {
        bonvoisin = false;
    }
    bonvoisin
}

/// OCCT BRepAdaptor_Curve::Resolution(R3d) — the parameter step for a 3d
/// tolerance (BRepAdaptor_Curve.cxx Resolution over the curve kind; the
/// rcad line/circle forms use the geometric resolution).
fn curve_resolution(c: &rcad_kernel::geom::Curve3, tol3d: f64) -> f64 {
    use rcad_kernel::geom::CurveEval as _;
    match c {
        rcad_kernel::geom::Curve3::Line(l) => {
            // OCCT: Res = 3d / |tangent| (a line has a unit-speed arc length
            // parameterization when the direction is normalized).
            let mag = l.direction.length();
            if mag > 0.0 {
                tol3d / mag
            } else {
                tol3d
            }
        }
        _ => {
            // OCCT L58-77: sample-based resolution (Icc = point/derivative).
            let d = c.default_domain();
            let t = 0.5 * (d[0] + d[1]);
            let dv = c.derivative_at(t);
            let mag = dv.length();
            if mag > 0.0 {
                tol3d / mag
            } else {
                tol3d
            }
        }
    }
}

/// OCCT ChFi3d_Builder_2.cxx L318-369 — Projection.
pub(crate) fn projection(
    pext: &mut Option<ExtPC>,
    p: DVec3,
    c: &rcad_kernel::geom::Curve3,
    w: &mut f64,
    tol: f64,
) -> bool {
    use rcad_kernel::geom::CurveEval as _;
    let d = c.default_domain();
    let mut dist2 = c.point_at(*w).distance_squared(p);

    // It is checked if it is not already a solution
    if dist2 < tol * tol {
        return true;
    }

    let mut ok = false;

    // On essai une resolution initialise
    if let Some(poc) = extrema_locate_ext_pc(p, c, *w, d[0], d[1], tol / 10.0) {
        let daux2 = c.point_at(poc.param).distance_squared(p);
        if daux2 < dist2 {
            *w = poc.param;
            dist2 = daux2;
            ok = true;
            if dist2 < tol * tol {
                return true;
            }
        }
    }

    // Global resolution
    if let Some(pext) = pext {
        pext.perform(p, c, d[0], d[1]);
        if pext.is_done() {
            for ii in 1..=pext.nb_ext() {
                if pext.square_distance(ii) < dist2 {
                    dist2 = pext.square_distance(ii);
                    *w = pext.point(ii).param;
                    ok = true;
                }
            }
        }
    }
    ok
}

/// OCCT ChFi3d_Builder_2.cxx L373-392 — TgtKP.
#[allow(dead_code)]
pub(crate) fn tgt_kp(
    cd: &ChFiDSSurfData,
    spine: &ChFiDSSpineHandle,
    iedge: usize,
    isfirst: bool,
) -> (DVec3, DVec3) {
    let wtg = cd.interference_on_s1().parameter(isfirst);
    let bc = spine.base().edges(iedge).clone();
    let bc_first = bc.as_edge().expect("not an edge").range[0];
    let bc_last = bc.as_edge().expect("not an edge").range[1];
    if spine.base().edges(iedge).orientation == Orientation::Forward {
        brep_adaptor_curve_d1(&bc, wtg + bc_first)
    } else {
        let (ped, mut ded) = brep_adaptor_curve_d1(&bc, -wtg + bc_last);
        ded = -ded;
        (ped, ded.normalize())
    }
}

/// OCCT ChFi3d_Builder_2.cxx L394-471 — IsInput.
pub(crate) fn is_input(brep: &topods::BRep, vec: DVec3, ve: &Shape, fa: &Shape) -> bool {
    let mut trouve = 0usize;
    let mut vec3d = [DVec3::ZERO; 2];
    let mut point = DVec3::ZERO;

    // Find edges and compute 3D vectors (TopExp_Explorer(Fa, WIRE) x
    // BRepTools_WireExplorer — the rcad face-edge walk keeps the order).
    for e in topexp_face_edges(brep, fa) {
        if trouve >= 2 {
            break;
        }
        let e = e;
        let (vf, vl) = topexp_vertices(&e);
        if vf.is_same(ve) {
            let par = brep_tool_parameter(brep, ve, &e);
            let (pt, tg) = brep_adaptor_curve_d1(&e, par);
            point = pt;
            vec3d[trouve] = tg;
            trouve += 1;
        } else if vl.is_same(ve) {
            let par = brep_tool_parameter(brep, ve, &e);
            let (pt, tg) = brep_adaptor_curve_d1(&e, par);
            point = pt;
            vec3d[trouve] = -tg;
            trouve += 1;
        }
        let _ = e;
    }
    if trouve < 2 {
        return false;
    }
    // Calculate the normal and the angles in the associated vector plane
    let normal = vec3d[0].cross(vec3d[1]);
    if normal.length_squared() < CONFUSION {
        // Colinear case
        return vec_is_parallel(vec, vec3d[0], CONFUSION);
    }

    let amin;
    let mut amax = vec_angle_with_ref(vec3d[1], vec3d[0], normal);
    if amax < 0.0 {
        amin = amax;
        amax = 0.0;
    } else {
        amin = 0.0;
    }

    // Projection of the vector (gp_Ax3 Axe(Point, Normal, Vec3d[0]) +
    // gp_Trsf::SetTransformation — global->local, zero z, then back).
    let xdir = vec3d[0].normalize();
    let zdir = normal.normalize();
    let ydir = zdir.cross(xdir);
    let mut coord = vec;
    coord = trsf_global_to_local(point, xdir, ydir, zdir, coord);
    coord.z = 0.0;
    coord = trsf_local_to_global(point, xdir, ydir, zdir, coord);
    let the_proj = coord;

    // and finally...
    let angle = vec_angle_with_ref(vec3d[0], the_proj, normal);
    (angle >= amin) && (angle <= amax)
}

/// OCCT ChFi3d_Builder_2.cxx L473-522 — IsG1.
pub(crate) fn is_g1(
    brep: &topods::BRep,
    the_map: &ChFiDSMap,
    e: &Shape,
    f_ref: &Shape,
    f_voi: &mut Shape,
) -> bool {
    // Find a neighbor of E different from FRef (general case).
    for it in the_map.find(e).clone() {
        let f = it;
        if !f.is_same(f_ref) {
            *f_voi = f;
            if is_tangent_faces(brep, e, f_ref, f_voi, GeomAbsShape::G1) {
                return true;
            }
        }
    }
    // If is was not found it is checked if E is a cutting edge,
    // in which case FVoi = FRef is returned (less frequent case).
    let mut orset = false;
    let mut orient = Orientation::Forward;
    for ex in topexp_face_edges(brep, f_ref) {
        let ed = ex;
        if ed.is_same(e) {
            if !orset {
                orient = ed.orientation;
                orset = true;
            } else if ed.orientation == topabs_reverse(orient) {
                *f_voi = f_ref.clone();
                return is_tangent_faces(brep, e, f_ref, f_ref, GeomAbsShape::G1);
            }
        }
    }
    false
}

/// OCCT ChFi3d_Builder_2.cxx L524-589 — SearchFaceOnV.  Returns the number
/// of candidate faces; F1/F2 are the out faces.
pub(crate) fn search_face_on_v(
    brep: &topods::BRep,
    pc: &ChFiDS_CommonPoint,
    f_ref: &Shape,
    ve_map: &ChFiDSMap,
    ef_map: &ChFiDSMap,
    f1: &mut Shape,
    f2: &mut Shape,
) -> i32 {
    // it is checked that it leaves the current face.
    let mut find_face = is_input(brep, pc.vector, pc.vertex(), f_ref);
    if find_face {
        find_face = is_input(brep, -pc.vector, pc.vertex(), f_ref);
    }
    // If it does not leave, it is finished
    if find_face {
        *f1 = f_ref.clone();
        return 1;
    }
    let mut num = 0i32;
    let mut trouve;
    let mut f_voi = Shape::null();

    for it_e in ve_map.find(pc.vertex()).clone() {
        if num >= 2 {
            break;
        }
        let e = it_e;
        trouve = false;
        for it_f in ef_map.find(&e).clone() {
            if it_f.is_same(f_ref) {
                trouve = true;
                break;
            }
        }
        if trouve {
            trouve = is_g1(brep, ef_map, &e, f_ref, &mut f_voi);
        }
        if trouve {
            trouve = is_input(brep, pc.vector, pc.vertex(), &f_voi);
        }
        if trouve {
            if num == 0 {
                *f1 = f_voi.clone();
            } else {
                *f2 = f_voi.clone();
            }
            num += 1;
        }
    }
    num
}

/// OCCT ChFi3d_Builder_2.cxx L591-628 — ChangeTransition.
pub(crate) fn change_transition(
    brep: &topods::BRep,
    precedant: &ChFiDS_CommonPoint,
    courant: &mut ChFiDS_CommonPoint,
    face_index: i32,
    ds: &TopOpeBRepDSHDataStructure,
) {
    let mut tochange = true;
    let f = ds.shape(face_index).clone();
    let arc = precedant.arc().clone();
    let pcurve1 = get_curve_on_surface(brep, &arc, &f);
    let arc_rev = {
        let mut x = arc.clone();
        x.orientation = topabs_reverse(x.orientation);
        x
    };
    let pcurve2 = get_curve_on_surface(brep, &arc_rev, &f);
    // OCCT compares the Geom2d_Curve handles by identity; both queries
    // address the same stored pcurve slot of (Arc, F), so the OCCT handles
    // are equal whenever the rcad payloads evaluate identically.
    if pcurve_payloads_differ(&pcurve1, &pcurve2) {
        // This is a cutting edge, it is necessary to make a small Geometric test
        let (_p, tgarc) = brep_adaptor_curve_d1(&arc, precedant.parameter_on_arc());
        tochange = vec_is_parallel(tgarc, precedant.vector, CONFUSION);
    }

    if tochange {
        courant.set_arc(
            CONFUSION,
            arc,
            precedant.parameter_on_arc(),
            topabs_reverse(precedant.transition_on_arc()),
        );
    }
}

/// The rcad encoding of the OCCT pcurve handle identity test (both handles
/// address the same stored slot -> the payloads are compared).
fn pcurve_payloads_differ(a: &Option<rcad_kernel::geom::Curve2d>, b: &Option<rcad_kernel::geom::Curve2d>) -> bool {
    match (a, b) {
        (None, None) => false,
        (Some(_), None) | (None, Some(_)) => true,
        (Some(ca), Some(cb)) => {
            let da = ca.default_domain();
            let db = cb.default_domain();
            let mid_a = 0.5 * (da[0] + da[1]);
            let mid_b = 0.5 * (db[0] + db[1]);
            ca.point_at(mid_a).distance(cb.point_at(mid_b)) > P_CONFUSION
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder::CallPerformSurf — L630-818.
// =========================================================================

/// Pending-leaf stand-in of ChFi3d_Builder::SimulSurf (the 2-surface
/// overload, SD out) — the owning translation (ChFi3d_Builder_C2/CnCrn
/// family) has not landed; the OCCT failure path (IsDone() == false) is
/// reported.
#[allow(clippy::too_many_arguments)]
fn simul_surf_2faces_pending(
    _sd: &SharedSurfData,
    _hguide: &ChFiDSElSpine,
    _spine: &ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_on_s1: bool,
    _rec_on_s2: bool,
    _soldep: &[f64; 4],
    _intf: &mut i32,
    _intl: &mut i32,
) -> bool {
    false
}

/// Pending-leaf stand-in of ChFi3d_Builder::PerformSurf (the 2-surface
/// overload, SeqSD out) — see simul_surf_2faces_pending.
#[allow(clippy::too_many_arguments)]
fn perform_surf_2faces_pending(
    _seqsd: &mut Vec<SharedSurfData>,
    _hguide: &ChFiDSElSpine,
    _spine: &ChFiDSSpineHandle,
    _choix: i32,
    _hs1: &BRepAdaptorSurface,
    _it1: &BRepTopAdaptorTopolTool,
    _hs2: &BRepAdaptorSurface,
    _it2: &BRepTopAdaptorTopolTool,
    _max_step: f64,
    _fleche: f64,
    _tolesp: f64,
    _first: &mut f64,
    _last: &mut f64,
    _inside_f: bool,
    _inside_l: bool,
    _forward: bool,
    _rec_on_s1: bool,
    _rec_on_s2: bool,
    _soldep: &[f64; 4],
    _intf: &mut i32,
    _intl: &mut i32,
) -> bool {
    false
}

impl super::chfi3d::ChFi3dBuilder {
    /// OCCT ChFi3d_Builder_2.cxx L630-818 — CallPerformSurf (encapsulates
    /// the call to PerformSurf/SimulSurf).
    /// OCCT's unnamed parameters (TolGuide, Appro) are dropped.
    #[allow(clippy::too_many_arguments)]
    pub fn call_perform_surf(
        &mut self,
        stripe: &SharedStripe,
        simul: bool,
        seqsd: &mut Vec<SharedSurfData>,
        sd: &SharedSurfData,
        hguide: &ChFiDSElSpine,
        spine: &ChFiDSSpineHandle,
        hs1: &mut BRepAdaptorSurface,
        hs3: &mut Option<BRepAdaptorSurface>,
        pp1: DVec2,
        pp3: DVec2,
        it1: &mut BRepTopAdaptorTopolTool,
        hs2: &mut BRepAdaptorSurface,
        hs4: &mut Option<BRepAdaptorSurface>,
        pp2: DVec2,
        pp4: DVec2,
        it2: &mut BRepTopAdaptorTopolTool,
        max_step: f64,
        fleche: f64,
        first: &mut f64,
        last: &mut f64,
        inside: bool,
        forward: bool,
        rec_on_s1: bool,
        rec_on_s2: bool,
        soldep: &mut [f64; 4],
        intf: &mut i32,
        intl: &mut i32,
        surf1: &mut Option<BRepAdaptorSurface>,
        surf2: &mut Option<BRepAdaptorSurface>,
    ) {
        let mut hson1 = hs1.clone();
        let mut hson2 = hs2.clone();
        // Definition of the domain of path It1, It2
        it1.initialize_brep(&hson1);
        it2.initialize_brep(&hson2);

        let mut or1 = hs1.face.orientation;
        let mut or2 = hs2.face.orientation;
        let (stripe_or1, stripe_or2, stripe_choix) = {
            let st = stripe.read().expect("stripe lock");
            (
                st.orientation_on_face1(),
                st.orientation_on_face2(),
                st.choix(),
            )
        };
        let mut choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, stripe_choix);
        soldep[0] = pp1.x;
        soldep[1] = pp1.y;
        soldep[2] = pp2.x;
        soldep[3] = pp2.y;

        let thef = *first;
        let thel = *last;
        let mut isdone;

        if simul {
            isdone = simul_surf_2faces_pending(
                sd, hguide, spine, choix, hs1, it1, hs2, it2, self.tolesp, first, last, inside,
                inside, forward, rec_on_s1, rec_on_s2, soldep, intf, intl,
            );
        } else {
            isdone = perform_surf_2faces_pending(
                seqsd, hguide, spine, choix, hs1, it1, hs2, it2, max_step, fleche, self.tolesp,
                first, last, inside, inside, forward, rec_on_s1, rec_on_s2, soldep, intf, intl,
            );
        }

        // Case of error
        if !isdone {
            *first = thef;
            *last = thel;
            let mut reprise = false;
            if let Some(hs3v) = hs3 {
                hson1 = hs3v.clone();
                it1.initialize_brep(&hson1);
                or1 = hs3v.face.orientation;
                soldep[0] = pp3.x;
                soldep[1] = pp3.y;
                reprise = true;
            } else if let Some(hs4v) = hs4 {
                hson2 = hs4v.clone();
                it2.initialize_brep(&hson2);
                or2 = hs4v.face.orientation;
                soldep[2] = pp4.x;
                soldep[3] = pp4.y;
                reprise = true;
            }

            if reprise {
                choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, stripe_choix);
                if simul {
                    isdone = simul_surf_2faces_pending(
                        sd, hguide, spine, choix, &hson1, it1, &hson2, it2, self.tolesp, first,
                        last, inside, inside, forward, rec_on_s1, rec_on_s2, soldep, intf, intl,
                    );
                } else {
                    isdone = perform_surf_2faces_pending(
                        seqsd, hguide, spine, choix, &hson1, it1, &hson2, it2, max_step, fleche,
                        self.tolesp, first, last, inside, inside, forward, rec_on_s1, rec_on_s2,
                        soldep, intf, intl,
                    );
                }
            }
        }
        *surf1 = Some(hson1);
        *surf2 = Some(hson2);
    }

    /// OCCT ChFi3d_Builder_2.cxx L891-1109 — StartSol (the Stripe overload;
    /// calculates a starting solution by sampling the spine).
    #[allow(clippy::too_many_arguments)]
    pub fn start_sol_on_stripe(
        &mut self,
        stripe: &SharedStripe,
        hguide: &mut ChFiDSElSpine,
        hs1: &mut BRepAdaptorSurface,
        hs2: &mut BRepAdaptorSurface,
        i1: &mut BRepTopAdaptorTopolTool,
        i2: &mut BRepTopAdaptorTopolTool,
        p1: &mut DVec2,
        p2: &mut DVec2,
        first: &mut f64,
    ) {
        // OCCT: occ::handle<ChFiDS_Spine>& Spine = Stripe->ChangeSpine();
        let mut spine = {
            let mut st = stripe.write().expect("stripe lock");
            st.my_spine.take()
        }
        .expect("null spine");
        let els = hguide;
        let nbed = spine.base().nb_edges();
        let mut nbessaimax = 3 * nbed;
        if nbessaimax < 10 {
            nbessaimax = 10;
        }
        let unsurnbessaimax = 1.0 / nbessaimax as f64;
        let wf = 0.9981 * spine.base().first_parameter_of(1)
            + 0.0019 * spine.base().last_parameter_of(1);
        let wl = 0.9973 * spine.base().last_parameter_of(nbed)
            + 0.0027 * spine.base().first_parameter_of(nbed);

        let mut tol_e = 1.0e-7;

        let mut nbessai;
        let mut iedge = 0usize;
        let rc = {
            let st = stripe.read().expect("stripe lock");
            st.choix()
        };
        // OCCT: gp_Vec2d derive; gp_Pnt2d P2d; TopoDS_Edge cured;
        // TopoDS_Face f1, f2 — the uninitialized declarations.
        let mut derive = DVec2::ZERO;
        let mut p2d = DVec2::ZERO;
        let mut cured = Shape::null();
        let (mut f1, mut f2) = (Shape::null(), Shape::null());
        let (mut or1, mut or2) = (Orientation::Forward, Orientation::Forward);
        let mut choix = 0i32;
        let mut sol_dep = [0.0f64; 4];
        let mut pc: Option<rcad_kernel::geom::Curve2d>;
        // OCCT L936-940: Extrema_ExtPC PExt; PExt.Initialize(els, First,
        // Last, Confusion) — deferred construction; the pending ElSpine
        // curve (adaptor_curve() == None) keeps PExt empty.
        let mut pext: Option<ExtPC> = None;
        let (mut pos1, mut pos2) = (TopAbsState::Unknown, TopAbsState::Unknown);
        nbessai = 0;
        while nbessai <= nbessaimax {
            let t = nbessai as f64 * unsurnbessaimax;
            let mut w = wf * (1.0 - t) + wl * t;
            let ie = spine.base().index_of_param(w, true);
            if iedge != ie {
                iedge = ie;
                cured = spine.base().edges(iedge).clone();
                tol_e = cured.as_edge().expect("not an edge").tolerance;
                let (ff1, ff2) = self.conex_faces(&spine, iedge);
                // OCCT L952: ConexFaces(Spine, iedge, HS1, HS2).
                *hs1 = BRepAdaptorSurface::initialize(&self.my_brep, &ff1);
                *hs2 = BRepAdaptorSurface::initialize(&self.my_brep, &ff2);
                f1 = hs1.face.clone();
                f2 = hs2.face.clone();
                or1 = f1.orientation;
                or2 = f2.orientation;
                let (stripe_or1, stripe_or2) = {
                    let st = stripe.read().expect("stripe lock");
                    (st.orientation_on_face1(), st.orientation_on_face2())
                };
                choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, rc);
            }

            let woned = spine.base_mut().parameter_on(w, iedge, true);
            cured.orientation = Orientation::Forward;
            let mut f1forward = f1.clone();
            let mut f2forward = f2.clone();
            f1forward.orientation = Orientation::Forward;
            f2forward.orientation = Orientation::Forward;
            pc = get_curve_on_surface(&self.my_brep, &cured, &f1forward);
            i1.initialize_brep(hs1);
            {
                let (pt, dv) = pc
                    .as_ref()
                    .map(|pc| {
                        let pt = pc.point_at(woned);
                        let dv = pc.derivative_at(woned);
                        (pt, dv)
                    })
                    .expect("pcurve");
                *p1 = pt;
                derive = dv;
            }
            // There are points on the border, and internal points are found
            if derive.length() > P_CONFUSION {
                derive = derive.normalize();
                // OCCT: derive.Rotate(M_PI / 2) — (x, y) -> (-y, x).
                derive = DVec2::new(-derive.y, derive.x);
                let as_surf = f1
                    .as_face()
                    .and_then(|fd| fd.surface.clone())
                    .expect("face surface");
                let res_u = GeomAdaptorSurface::new(as_surf.clone()).u_resolution(tol_e);
                let res_v = GeomAdaptorSurface::new(as_surf).v_resolution(tol_e);
                derive *= 2.0 * (derive.x.abs() * res_u + derive.y.abs() * res_v);
                p2d = *p1 + derive;
                if i1.classify(p2d, res_u.min(res_v), false) == TopAbsState::In {
                    *p1 = p2d;
                } else {
                    p2d = *p1 - derive;
                    if i1.classify(p2d, res_u.min(res_v), false) == TopAbsState::In {
                        *p1 = p2d;
                    }
                }
            }
            if f1.is_same(&f2) {
                cured.orientation = Orientation::Reversed;
            }
            pc = get_curve_on_surface(&self.my_brep, &cured, &f2forward);
            {
                *p2 = pc.map(|pc| pc.point_at(woned)).unwrap_or(DVec2::ZERO);
            }
            i2.initialize_brep(hs2);

            sol_dep[0] = p1.x;
            sol_dep[1] = p1.y;
            sol_dep[2] = p2.x;
            sol_dep[3] = p2.y;
            let pnt = brep_adaptor_curve_value(spine.base().edges(iedge), woned);

            // OCCT L1009-1010: Projection over the pending ElSpine curve.
            let projection_ok = els
                .adaptor_curve()
                .map_or(false, |c| projection(&mut pext, pnt, &c, &mut w, self.tolapp3d));
            if projection_ok
                && self.perform_first_section(
                    &spine, els, choix, hs1, hs2, i1, i2, w, &mut sol_dep, &mut pos1, &mut pos2,
                )
            {
                *p1 = DVec2::new(sol_dep[0], sol_dep[1]);
                *p2 = DVec2::new(sol_dep[2], sol_dep[3]);
                *first = w;
                {
                    let mut st = stripe.write().expect("stripe lock");
                    st.my_spine = Some(spine);
                }
                return;
            }
            nbessai += 1;
        }
        // No solution was found for the faces adjacent to the trajectory.
        // Now one tries the neighbor faces.
        iedge = 0;
        nbessai = 0;
        while nbessai <= nbessaimax {
            let t = nbessai as f64 * unsurnbessaimax;
            let mut w = wf * (1.0 - t) + wl * t;
            iedge = spine.base().index_of_param(w, true);
            cured = spine.base().edges(iedge).clone();
            let (ff1, ff2) = self.conex_faces(&spine, iedge);
            *hs1 = BRepAdaptorSurface::initialize(&self.my_brep, &ff1);
            *hs2 = BRepAdaptorSurface::initialize(&self.my_brep, &ff2);
            f1 = hs1.face.clone();
            f2 = hs2.face.clone();
            or1 = f1.orientation;
            or2 = f2.orientation;
            let (stripe_or1, stripe_or2) = {
                let st = stripe.read().expect("stripe lock");
                (st.orientation_on_face1(), st.orientation_on_face2())
            };
            choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, rc);
            let woned = spine.base_mut().parameter_on(w, iedge, true);
            let mut f1forward = f1.clone();
            let mut f2forward = f2.clone();
            f1forward.orientation = Orientation::Forward;
            f2forward.orientation = Orientation::Forward;
            pc = get_curve_on_surface(&self.my_brep, &cured, &f1forward);
            {
                *p1 = pc.map(|pc| pc.point_at(woned)).unwrap_or(DVec2::ZERO);
            }
            pc = get_curve_on_surface(&self.my_brep, &cured, &f2forward);
            {
                *p2 = pc.map(|pc| pc.point_at(woned)).unwrap_or(DVec2::ZERO);
            }
            i1.initialize_brep(hs1);
            i2.initialize_brep(hs2);
            sol_dep[0] = p1.x;
            sol_dep[1] = p1.y;
            sol_dep[2] = p2.x;
            sol_dep[3] = p2.y;
            let pnt = brep_adaptor_curve_value(spine.base().edges(iedge), woned);
            // OCCT L1056: Projection over the pending ElSpine curve.
            if els
                .adaptor_curve()
                .map_or(false, |c| projection(&mut pext, pnt, &c, &mut w, self.tolapp3d))
            {
                self.perform_first_section(
                    &spine, els, choix, hs1, hs2, i1, i2, w, &mut sol_dep, &mut pos1, &mut pos2,
                );
                let (p, v) = els.d1(w);
                let pl = Surface3::Plane(rcad_kernel::geom::Plane::new(p, v));
                let plane = GeomAdaptorSurface::new(pl);

                let mut bonvoisin = true;
                let mut found = false;
                let mut nb_changement = 1;
                while bonvoisin && (!found) && (nb_changement < 5) {
                    if pos1 != TopAbsState::In {
                        bonvoisin = bon_voisin(
                            &self.my_brep,
                            p,
                            hs1,
                            &mut f1,
                            &plane,
                            &cured,
                            &mut sol_dep,
                            0,
                            1,
                            &self.my_ef_map,
                            self.tolapp3d,
                        );
                    }
                    if pos2 != TopAbsState::In && bonvoisin {
                        bonvoisin = bon_voisin(
                            &self.my_brep,
                            p,
                            hs2,
                            &mut f2,
                            &plane,
                            &cured,
                            &mut sol_dep,
                            2,
                            3,
                            &self.my_ef_map,
                            self.tolapp3d,
                        );
                    }
                    if bonvoisin {
                        f1 = hs1.face.clone();
                        f2 = hs2.face.clone();
                        or1 = f1.orientation;
                        or2 = f2.orientation;
                        let (stripe_or1, stripe_or2) = {
                            let st = stripe.read().expect("stripe lock");
                            (st.orientation_on_face1(), st.orientation_on_face2())
                        };
                        choix = next_side(&mut or1, &mut or2, stripe_or1, stripe_or2, rc);
                        let hson1new = hs1.clone();
                        let hson2new = hs2.clone();
                        i1.initialize_brep(&hson1new);
                        i2.initialize_brep(&hson2new);
                        if self.perform_first_section(
                            &spine, els, choix, hs1, hs2, i1, i2, w, &mut sol_dep, &mut pos1,
                            &mut pos2,
                        ) {
                            *p1 = DVec2::new(sol_dep[0], sol_dep[1]);
                            *p2 = DVec2::new(sol_dep[2], sol_dep[3]);
                            *first = w;
                            found = true;
                        }
                    }
                    nb_changement += 1;
                }
                if found {
                    {
                        let mut st = stripe.write().expect("stripe lock");
                        st.my_spine = Some(spine);
                    }
                    return;
                }
            }
            nbessai += 1;
        }
        {
            spine.base_mut().set_error_status(ChFiDS_ErrorStatus::StartsolFailure);
            let mut st = stripe.write().expect("stripe lock");
            st.my_spine = Some(spine);
        }
        panic!("Standard_Failure: StartSol echec");
    }

    /// OCCT ChFi3d_FilBuilder.cxx L1500-1534 — PerformFirstSection(Spine,
    /// HGuide, Choix, S1, S2, I1, I2, Par, SolDep, Pos1, Pos2).  The
    /// FilBuilder override called from StartSol (Builder_2.cxx L1058/L1092);
    /// rcad models the override on the shared ChFi3dBuilder (the blend
    /// builder configuration), which is where the StartSol calls land.
    #[allow(clippy::too_many_arguments)]
    fn perform_first_section(
        &mut self,
        spine: &ChFiDSSpineHandle,
        hguide: &ChFiDSElSpine,
        choix: i32,
        s1: &BRepAdaptorSurface,
        s2: &BRepAdaptorSurface,
        i1: &BRepTopAdaptorTopolTool,
        i2: &BRepTopAdaptorTopolTool,
        par: f64,
        soldep: &mut [f64; 4],
        pos1: &mut TopAbsState,
        pos2: &mut TopAbsState,
    ) -> bool {
        // OCCT L1512-1516: down_cast<ChFiDS_FilSpine>(Spine).
        let Some(fsp) = spine.down_cast_fil() else {
            panic!(
                "Standard_ConstructionError: PerformSurf : this is not the spine of a fillet"
            );
        };
        // OCCT L1517: TolGuide = HGuide->Resolution(tolapp3d).
        let tol_guide = hguide.resolution(self.tolapp3d);
        if fsp.is_constant() {
            // OCCT L1520-1524: BRepBlend_ConstRad Func(S1, S2, HGuide);
            // Func.Set(fsp->Radius(), Choix); Func.Set(myShape);
            // BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide);
            // return TheWalk.PerformFirstSection(Func, Par, SolDep,
            //     tolapp3d, TolGuide, Pos1, Pos2);
            // (the second Set is the BlendFunc_SectionShape overload — the
            // builder's myShape / BlendFunc_SectionShape member,
            // ChFi3d_FilBuilder.hxx L377.)
            let guide = hguide
                .curve
                .as_ref()
                .expect("ElSpine curve (ChFiDS_ElSpine carries a loaded curve)");
            let mut func = BlendFuncConstRad::new(&s1.surface, &s2.surface, guide);
            func.set(fsp.radius(), choix);
            func.set_section_shape(self.my_blend_shape);
            let mut the_walk = BRepBlendWalking::new(&s1.surface, &s2.surface, i1, i2, guide);
            the_walk.perform_first_section(
                &mut func,
                par,
                soldep,
                self.tolapp3d,
                tol_guide,
                pos1,
                pos2,
            )
        } else {
            // OCCT L1526-1533: BRepBlend_EvolRad Func(S1, S2, HGuide,
            // fsp->Law(HGuide)); Func.Set(Choix); Func.Set(myShape);
            // BRepBlend_Walking TheWalk(S1, S2, I1, I2, HGuide);
            // return TheWalk.PerformFirstSection(Func, Par, SolDep,
            //     tolapp3d, TolGuide, Pos1, Pos2);
            let guide = hguide
                .curve
                .as_ref()
                .expect("ElSpine curve (ChFiDS_ElSpine carries a loaded curve)");
            // OCCT: fsp->Law(HGuide) — a null law cannot exist in OCCT
            // (AppendElSpine has appended the ComputeLaw composite); the
            // expect keeps the OCCT null-handle failure path explicit.
            let law = fsp
                .law_of(hguide)
                .expect("ChFiDS_FilSpine::Law: no law for this elspine");
            let mut func = BlendFuncEvolRad::new(&s1.surface, &s2.surface, guide, law);
            func.set(choix);
            func.set_section_shape(self.my_blend_shape);
            let mut the_walk = BRepBlendWalking::new(&s1.surface, &s2.surface, i1, i2, guide);
            the_walk.perform_first_section(
                &mut func,
                par,
                soldep,
                self.tolapp3d,
                tol_guide,
                pos1,
                pos2,
            )
        }
    }

    /// OCCT ChFi3d_Builder_2.cxx L1146-1519 — StartSol (the Spine overload;
    /// chains the obstacle over the common point).
    #[allow(clippy::too_many_arguments)]
    pub fn start_sol(
        &self,
        spine: &mut ChFiDSSpineHandle,
        hs: &mut BRepAdaptorSurface, // New face
        pons: &mut DVec2,            // Localization
        hc: &mut Option<BRepAdaptorCurve2d>, // Representation of the obstacle
        w: &mut f64,
        sd: &ChFiDSSurfData,
        isfirst: bool,
        ons: i32,
        hsref: &mut BRepAdaptorSurface, // The other representation
        hcref: &mut Option<BRepAdaptorCurve2d>, // of the obstacle
        rec_p: &mut bool,
        rec_s: &mut bool,
        rec_rst: &mut bool,
        c1obstacle: &mut bool,
        hsbis: &mut Option<BRepAdaptorSurface>, // Face of support
        pbis: &mut DVec2,                       // and its point
        decroch: bool,
        vref: &Shape,
    ) -> bool {
        *rec_p = false;
        *rec_s = false;
        *rec_rst = false;
        *c1obstacle = false;
        let mut fv = Shape::null();
        let mut a_pcurve: Option<rcad_kernel::geom::Curve2d> = None;

        let fref = hsref.face.clone();
        let f = self
            .my_ds
            .as_ref()
            .expect("DS")
            .shape(sd.index_of(ons))
            .clone();

        let a_common_point = sd.vertex(isfirst, ons).clone();
        *hsbis = None;

        if a_common_point.is_on_arc() {
            let notons = if ons == 1 { 2 } else { 1 };
            let cpbis = sd.vertex(isfirst, notons).clone();
            if cpbis.is_on_arc() {
                // It is checked if it is not the extension zone
                let ts = sd.interference(ons).parameter(isfirst);
                let tns = sd.interference(notons).parameter(isfirst);
                // Arbitrary test (to precise)
                let is_extend = if isfirst {
                    ts - tns > 100.0 * self.tolesp
                } else {
                    tns - ts > 100.0 * self.tolesp
                };
                if is_extend && a_common_point.point.distance(cpbis.point) > 0.0 {
                    //  the state is preserved and False is returned
                    //  (extension by the expected plane).
                    *hs = BRepAdaptorSurface::initialize(&self.my_brep, &f);
                    let pc = sd.interference(ons).pcurve_on_face().cloned();
                    // The 2nd point is given by its trace on the support surface
                    *rec_s = false;
                    *pons = pc
                        .map(|pc| pc.point_at(tns))
                        .expect("PCurveOnFace");
                    return false;
                }
            }
        }

        if a_common_point.is_vertex() && hc.is_some() && !decroch {
            // The edge is changed, the parameter is updated and
            // eventually the support face and(or) the reference face.
            let vcp = a_common_point.vertex().clone();
            let ehc = hc.as_ref().expect("HC").edge().clone();
            // One starts by searching in Fref another edge referencing VCP.
            let mut newedge = Shape::null();
            let mut edgereg = Shape::null();
            let mut facereg = Shape::null();
            let mut bidface = fref.clone();
            bidface.orientation = Orientation::Forward;
            for cured in topexp_face_edges(&self.my_brep, &bidface) {
                let mut found = false;
                if !cured.is_same(&ehc) {
                    let (v1, v2) = topexp_vertices(&cured);
                    for vx in [&v1, &v2] {
                        if vx.is_same(&vcp) && !found {
                            if is_g1(&self.my_brep, &self.my_ef_map, &cured, &fref, &mut fv) {
                                edgereg = cured.clone();
                                facereg = fv.clone();
                            } else {
                                found = true;
                            }
                        }
                    }
                }
                if found {
                    newedge = cured;
                    break;
                }
            }
            if newedge.is_null() {
                // It is checked if EHC is not a closed edge.
                let (v1, v2) = topexp_vertices(&ehc);
                if v1.is_same(&v2) {
                    newedge = ehc.clone();
                    let w1 = brep_tool_parameter(&self.my_brep, &v1, &ehc);
                    let w2 = brep_tool_parameter(&self.my_brep, &v2, &ehc);
                    let fi = sd.interference(ons);
                    let pcf = fi.pcurve_on_face().cloned();
                    let ww = fi.parameter(isfirst);
                    let pww = pcf
                        .map(|pc| pc.point_at(ww))
                        .unwrap_or_else(|| sd.get_2d_point(isfirst, ons));
                    let hc_v = hc.as_ref().expect("HC");
                    let p1 = hc_v.value(w1);
                    let p2 = hc_v.value(w2);

                    if p1.distance(pww) > p2.distance(pww) {
                        *w = w1;
                        *pons = p1;
                    } else {
                        *w = w2;
                        *pons = p2;
                    }
                    *rec_p = true;
                    *c1obstacle = true;
                    return true;
                } else if !edgereg.is_null() {
                    // the reference edge and face are changed.
                    let fref2 = facereg.clone();
                    *hsref = BRepAdaptorSurface::initialize(&self.my_brep, &fref2);
                    let facereg2 = facereg.clone();
                    for cured in topexp_face_edges(&self.my_brep, &facereg2) {
                        if newedge.is_null() {
                            if !cured.is_same(&edgereg) {
                                let (v1, v2) = topexp_vertices(&cured);
                                for vx in [&v1, &v2] {
                                    if vx.is_same(&vcp)
                                        && !is_g1(
                                            &self.my_brep,
                                            &self.my_ef_map,
                                            &cured,
                                            &fref2,
                                            &mut fv,
                                        )
                                    {
                                        newedge = cured.clone();
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // it is necessary to find the new support face of the fillet :
            // connected to FRef along the newedge.
            if newedge.is_null() {
                panic!(
                    "Standard_Failure: StartSol : chain is not possible, new obstacle not found"
                );
            }
            if is_g1(&self.my_brep, &self.my_ef_map, &newedge, &fref, &mut fv) {
                panic!(
                    "Standard_Failure: StartSol : chain is not possible, config non processed"
                );
            } else if fv.is_null() {
                panic!(
                    "Standard_Failure: StartSol : chain is not possible, new obstacle not found"
                );
            } else {
                *hs = BRepAdaptorSurface::initialize(&self.my_brep, &fv);
                *w = brep_tool_parameter(&self.my_brep, &vcp, &newedge);
                if let Some(hcref) = hcref.as_mut() {
                    hcref.initialize(&self.my_brep, &newedge, &fref);
                }
                let mut newface = fv.clone();
                newface.orientation = Orientation::Forward;
                for ex in topexp_face_edges(&self.my_brep, &newface) {
                    if ex.is_same(&newedge) {
                        newedge = ex;
                        break;
                    }
                }
                if let Some(hc) = hc.as_mut() {
                    hc.initialize(&self.my_brep, &newedge, &fv);
                }
                *pons = hc.as_ref().expect("HC").value(*w);
            }
            *rec_p = true;
            *c1obstacle = true;
            return true;
        } // End of Case Vertex && Obstacle
        if a_common_point.is_on_arc() && hc.is_some() && !decroch {
            // Nothing is changed, the parameter is only updated.
            *w = a_common_point.parameter_on_arc();
            *c1obstacle = true;
            return true;
        }

        *hc = None;

        if a_common_point.is_on_arc() {
            let an_arc_edge = a_common_point.arc().clone();

            // Lambda to avoid code duplication.
            // Sets HS, W, aPCurve and pons to default return values.
            let prepare_default_return = |hs: &mut BRepAdaptorSurface,
                                              w: &mut f64,
                                              a_pcurve: &mut Option<rcad_kernel::geom::Curve2d>,
                                              pons: &mut DVec2| {
                *hs = BRepAdaptorSurface::initialize(&self.my_brep, &f);
                *w = a_common_point.parameter_on_arc();
                *a_pcurve = get_curve_on_surface(&self.my_brep, &an_arc_edge, &f);
                *pons = a_pcurve
                    .as_ref()
                    .map(|pc| pc.point_at(*w))
                    .expect("pcurve");
            };

            if decroch {
                *hs = BRepAdaptorSurface::initialize(&self.my_brep, &fref);
                *w = a_common_point.parameter_on_arc();
                a_pcurve = get_curve_on_surface(&self.my_brep, &an_arc_edge, &fref);
                *pons = a_pcurve
                    .as_ref()
                    .map(|pc| pc.point_at(*w))
                    .expect("pcurve");
                *rec_s = true;
                return true;
            }

            if self.search_face(spine, &a_common_point, &f, &mut fv) {
                *hs = BRepAdaptorSurface::initialize(&self.my_brep, &fv);
                *rec_s = true;
                if a_common_point.is_vertex() {
                    // One goes directly by the Vertex
                    // And it is checked that there are no other candidates
                    let mut aux = Shape::null();
                    let nb = search_face_on_v(
                        &self.my_brep,
                        &a_common_point,
                        &f,
                        &self.my_ve_map,
                        &self.my_ef_map,
                        &mut fv,
                        &mut aux,
                    );

                    *pons = brep_tool_parameters_vf(&self.my_brep, a_common_point.vertex(), &fv);
                    *hs = BRepAdaptorSurface::initialize(&self.my_brep, &fv);
                    if nb >= 2 {
                        *hsbis = Some(BRepAdaptorSurface::initialize(&self.my_brep, &aux));
                        *pbis = brep_tool_parameters_vf(&self.my_brep, a_common_point.vertex(), &aux);
                    }
                    return true;
                }
                // otherwise one passes by the arc...
                if !fv.is_same(&f) {
                    fv.orientation = Orientation::Forward;
                    let mut newedge = Shape::null();
                    for ex in topexp_face_edges(&self.my_brep, &fv) {
                        if ex.is_same(&an_arc_edge) {
                            newedge = ex;
                            break;
                        }
                    }
                    //  In cas of Tangent output, the current face becomes the support face
                    if is_tangent_to_arc(&a_common_point, 0.1) {
                        a_pcurve = get_curve_on_surface(&self.my_brep, a_common_point.arc(), &f);
                        *hsbis = Some(BRepAdaptorSurface::initialize(&self.my_brep, &f));
                        *pbis = a_pcurve
                            .as_ref()
                            .map(|pc| pc.point_at(a_common_point.parameter_on_arc()))
                            .expect("pcurve");
                    }

                    a_pcurve = get_curve_on_surface(&self.my_brep, &newedge, &fv);
                } else {
                    let mut newedge = an_arc_edge.clone();
                    newedge.orientation = topabs_reverse(newedge.orientation);
                    fv.orientation = Orientation::Forward;
                    a_pcurve = get_curve_on_surface(&self.my_brep, &newedge, &fv);
                }
                *pons = a_pcurve
                    .as_ref()
                    .map(|pc| pc.point_at(a_common_point.parameter_on_arc()))
                    .expect("pcurve");
                return true;
            } else if !fv.is_null() {
                *c1obstacle = true;
                if !vref.is_null() {
                    let (v1, v2) = topexp_vertices(&an_arc_edge);
                    for vx in [&v1, &v2] {
                        if vx.is_same(vref) {
                            *c1obstacle = false;
                            break;
                        }
                    }
                }
                if *c1obstacle {
                    *hs = BRepAdaptorSurface::initialize(&self.my_brep, &fv);
                    *hsref = BRepAdaptorSurface::initialize(&self.my_brep, &f);
                    *w = a_common_point.parameter_on_arc();
                    *hc = Some(BRepAdaptorCurve2d::new());
                    let mut newedge = Shape::null();
                    let mut newface = fv.clone();
                    newface.orientation = Orientation::Forward;
                    for ex in topexp_face_edges(&self.my_brep, &newface) {
                        if ex.is_same(&an_arc_edge) {
                            newedge = ex;
                            break;
                        }
                    }

                    if !newedge.is_null() {
                        if let Some(hc) = hc.as_mut() {
                            hc.initialize(&self.my_brep, &newedge, &fv);
                        }
                        *pons = hc.as_ref().expect("HC").value(*w);
                        if let Some(hcref) = hcref.as_mut() {
                            hcref.initialize(&self.my_brep, &an_arc_edge, &f);
                        }
                        if a_common_point.is_vertex() {
                            *rec_p = true;
                        } else {
                            *rec_rst = true;
                        }
                        return true;
                    } else {
                        prepare_default_return(hs, w, &mut a_pcurve, pons);
                        return false;
                    }
                } else {
                    prepare_default_return(hs, w, &mut a_pcurve, pons);
                    return false;
                }
            } else {
                // there is no neighbor face, the state is preserved and False is returned.
                prepare_default_return(hs, w, &mut a_pcurve, pons);
                return false;
            }
        } else {
            *hs = BRepAdaptorSurface::initialize(&self.my_brep, &f);
            let fi = sd.interference(ons);
            if fi.pcurve_on_face().is_none() {
                *pons = sd.get_2d_point(isfirst, ons);
            } else {
                *pons = fi
                    .pcurve_on_face()
                    .map(|pc| pc.point_at(fi.parameter(isfirst)))
                    .expect("PCurveOnFace");
            }
        }
        true
    }

    /// OCCT ChFi3d_Builder_2.cxx L1521-1627 — SearchFace.
    pub fn search_face(
        &self,
        spine: &mut ChFiDSSpineHandle,
        pc: &ChFiDS_CommonPoint,
        f_ref: &Shape,
        f_voi: &mut Shape,
    ) -> bool {
        let mut trouve;
        if !pc.is_on_arc() {
            return false;
        }
        *f_voi = Shape::null();
        let mut e;
        if pc.is_vertex() {
            // attention it is necessary to analyze all faces that turn around of the vertex
            if pc.has_vector() {
                // General processing
                let mut fbis = Shape::null();
                let nb_faces = search_face_on_v(
                    &self.my_brep,
                    pc,
                    f_ref,
                    &self.my_ve_map,
                    &self.my_ef_map,
                    f_voi,
                    &mut fbis,
                );
                return nb_faces > 0;
            } else {
                // Processing using the spine
                let mut find_face;
                let (point, mut vec_spine) = spine.base_mut().d1(pc.parameter());

                // It is checked if one leaves from the current face.
                find_face = is_input(&self.my_brep, vec_spine, pc.vertex(), f_ref);
                if find_face {
                    vec_spine = -vec_spine;
                    find_face = is_input(&self.my_brep, vec_spine, pc.vertex(), f_ref);
                }
                // If one does not leave, it is ended
                if find_face {
                    *f_voi = f_ref.clone();
                    return true;
                }
                let _ = point;

                // Otherwise one finds the next among shared Faces
                // by a common edge G1
                for it_e in self.my_ve_map.find(pc.vertex()).clone() {
                    if find_face {
                        break;
                    }
                    e = it_e;
                    trouve = false;
                    for it_f in self.my_ef_map.find(&e).clone() {
                        if it_f.is_same(f_ref) {
                            trouve = true;
                            break;
                        }
                    }
                    if trouve {
                        find_face = is_g1(&self.my_brep, &self.my_ef_map, &e, f_ref, f_voi);
                    }
                    if find_face {
                        find_face = false;
                        // OCCT L1589-1593: La Spine peut etre nulle (ThreeCorner).
                        // (the rcad handle is non-null by construction.)

                        // It is checked if the selected face actually possesses edges of the spine
                        // containing the vertex on its front
                        // This processing should go only if the Vertex belongs to the spine
                        // This is a single case, for other vertexes it is required to do other things
                        trouve = false;
                        let mut ie = 1usize;
                        while ie <= spine.base().nb_edges() && !trouve {
                            e = spine.base().edges(ie).clone();
                            let (ef, el) = topexp_vertices(&e);
                            if ef.is_same(pc.vertex()) || el.is_same(pc.vertex()) {
                                for it_f in self.my_ef_map.find(&e).clone() {
                                    if it_f.is_same(f_voi) {
                                        trouve = true;
                                        break;
                                    }
                                }
                            }
                            ie += 1;
                        }
                        find_face = trouve;
                    }
                }
            }
        } else {
            return is_g1(&self.my_brep, &self.my_ef_map, pc.arc(), f_ref, f_voi);
        }
        false
    }
}

/// OCCT BRep_Tool::CurveOnSurface(E, F) without the range (BuildPlane's
/// local use).
fn brep_curve_on_surface_raw(
    brep: &topods::BRep,
    e: &Shape,
    f: &Shape,
) -> Option<rcad_kernel::geom::Curve2d> {
    brep.curve_on_surface(e, f).map(|(pc, _, _)| pc)
}

/// OCCT ChFi3d_Builder_2.cxx L1111-1144 — ChFi3d_BuildPlane.
pub(crate) fn chfi3d_build_plane(
    dstr: &mut TopOpeBRepDSHDataStructure,
    hs: &mut BRepAdaptorSurface,
    pons: &mut DVec2,
    sd: &ChFiDSSurfData,
    isfirst: bool,
    ons: i32,
) {
    let f = dstr.shape(sd.index_of(ons)).clone();
    if f.is_null() {
        return;
    }

    if sd.vertex(isfirst, ons).is_on_arc() {
        let arc = sd.vertex(isfirst, ons).arc().clone();
        let hc = brep_curve_on_surface_raw(&hs.brep, &arc, &f);
        if let Some(hc) = &hc {
            let uv = hc.point_at(sd.vertex(isfirst, ons).parameter_on_arc());
            let (u, v) = (uv.x, uv.y);
            // OCCT L1132: BRepLProp_SLProps theProp(*HS, u, v, 1, 1.e-12).
            let (_p, du, dv) = hs.surface.derivatives(u, v);
            let n = du.cross(dv);
            if n.length() > CONFUSION {
                // OCCT L1135-1139: a bare face over the plane with the
                // orientation of F; rcad folds the bare-face construction
                // into initialize_surface and keeps F's identity and
                // orientation on the adaptor (architecture note,
                // chfi3d_builder_0).
                let pln = Surface3::Plane(rcad_kernel::geom::Plane::new(
                    hs.surface.point_at(u, v),
                    n.normalize(),
                ));
                let mut new_f = BRepAdaptorSurface::initialize_surface(pln);
                new_f.face = f.clone();
                new_f.brep = hs.brep.clone();
                *pons = DVec2::new(0.0, 0.0);
                *hs = new_f;
                return; // everything is good !
            }
        }
    }
    panic!("Standard_Failure: ChFi3d_BuildPlane : echec .");
}

