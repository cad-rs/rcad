//! OCCT ChFi3d_Builder_6.cxx — 1:1 translation (Stage 1f).
//!
//! Source: OCCT src/ModelingAlgorithms/TKFillet/ChFi3d/ChFi3d_Builder_6.cxx
//! (2,703 lines).  This file carries the file-static helpers (SearchIndex,
//! nbedconnex, IsVois, IsObst, CompParam, CompBlendPoint, UpdateLine), the
//! ChFi3d_Builder_0 helpers consumed here but not yet translated elsewhere
//! (ChFi3d_TrsfTrans, ChFi3d_FilCommonPoint — annotated with their Builder_0
//! line numbers), the four CompleteData overloads and StoreData.  The three
//! ComputeData and three SimulData overloads live in `chfi3d_builder_6b`
//! (2000-line guideline).
//!
//! Architecture mappings:
//!   - `occ::handle<Adaptor3d_Surface>` (face/DS-surface adaptor) maps to
//!     `&BRepAdaptorSurface` (rcad single adaptor kind, chfi3d_builder_0);
//!     the nullable S2 handle maps to `Option<&BRepAdaptorSurface>`.
//!   - `occ::handle<BRepBlend_Line>` maps to a plain `BRepBlendLine` value
//!     (Arc + RwLock, SharedStripe precedent).
//!   - `occ::handle<Geom2d_Curve>` maps to `Option<Curve2d>` / `Curve2d`.
//!   - `Geom2dInt_GInter` is bridged by the analytic AnaIntersection2d for
//!     conic pcurves; non-analytic pcurves fall to the OCCT !found branch
//!     (projection) inside CompParam.
//!
//! Cross-file references (Stage 1e second batch, translated in parallel;
//! types named per the OCCT class, signatures per the OCCT .hxx):
//!   BRepBlend_Line           -> super::brep_blend_line::BRepBlendLine
//!   BRepBlend_Extremity      -> super::brep_blend_extremity::BRepBlendExtremity
//!   BRepBlend_PointOnRst     -> super::brep_blend_point_on_rst::BRepBlendPointOnRst
//!   AppBlend_Approx          -> super::brep_blend_app_surf::AppBlendApprox
//!   BRepBlend_AppFunc        -> super::brep_blend_app_surf::BRepBlendAppFunc
//!   BRepBlend_AppFuncRst     -> super::brep_blend_app_surf::BRepBlendAppFuncRst
//!   BRepBlend_AppFuncRstRst  -> super::brep_blend_app_surf::BRepBlendAppFuncRstRst
//!   BRepBlend_AppSurface     -> super::brep_blend_app_surf::BRepBlendAppSurface
//!   Blend_SurfRstFunction    -> super::brep_blend_surf_rst_function::BlendSurfRstFunction
//!   Blend_RstRstFunction     -> super::brep_blend_rst_rst_function::BlendRstRstFunction

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Curve3, CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};
use rcad_kernel::topods;

use super::brep_blend_app_surf::{
    AppBlendApprox, BRepBlendAppFunc, BRepBlendAppFuncRst, BRepBlendAppFuncRstRst, BRepBlendAppSurface,
};
use super::brep_blend_extremity::BRepBlendExtremity;
use super::brep_blend_function::BlendFunction;
use super::brep_blend_line::BRepBlendLine;
use super::brep_blend_point::BlendPoint;
use super::brep_blend_point_on_rst::BRepBlendPointOnRst;
use super::brep_blend_rst_rst_function::BlendRstRstFunction;
use super::brep_blend_surf_rst_function::BlendSurfRstFunction;
use super::chfi3d::{topabs_compose, topabs_reverse, ChFi3dBuilder};
use super::chfi3d_builder_0::{brep_tool_parameter, chfi3d_same_parameter, topexp_vertices, BRepAdaptorSurface};
use super::chfi_ds::{ChFiDS_CommonPoint, ChFiDSMap, ChFiDSSurfData};
use super::chfi3d_ds::{TopOpeBRepDSCurve, TopOpeBRepDSSurface};
use crate::topalgo::adaptor3d::hvertex::HVertex;

/// OCCT IntSurf_TypeTrans (IntSurf_TypeTrans.hxx) — rcad carrier.
pub use crate::geomalgo::int_patch::transitions::TypeTrans as IntSurfTypeTrans;

// =========================================================================
// OCCT ChFi3d_Builder_6.cxx L94-103 — SearchIndex.
// =========================================================================
pub fn search_index(value: f64, lin: &BRepBlendLine) -> i32 {
    let nb_pnt = lin.nb_points();
    let mut ind = 1;
    while ind < nb_pnt && lin.point(ind).parameter() < value {
        ind += 1;
    }
    ind
}

// =========================================================================
// OCCT ChFi3d_Builder_6.cxx L107-130 — nbedconnex.
// =========================================================================
pub fn nbedconnex(l: &[Shape]) -> i32 {
    let mut nb = 0;
    for (i, curs) in l.iter().enumerate() {
        let mut dejavu = false;
        for j in 0..i {
            if curs.is_same(&l[j]) {
                dejavu = true;
                break;
            }
        }
        if !dejavu {
            nb += 1;
        }
    }
    nb
}

// =========================================================================
// OCCT ChFi3d_Builder_6.cxx L132-193 — IsVois.
// OCCT NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> DONE maps to
// a Vec with IsSame membership (shape identity, chfi3d_builder_0 note).
// =========================================================================
pub fn is_vois(
    e: &Shape,
    vref: &Shape,
    ve_map: &ChFiDSMap,
    done: &mut Vec<Shape>,
    prof: i32,
    profmax: i32,
) -> bool {
    if prof > profmax {
        return false;
    }
    if done.iter().any(|d| d.is_same(e)) {
        return false;
    }
    let (v1, v2) = topexp_vertices(e);
    if vref.is_same(&v1) || vref.is_same(&v2) {
        return true;
    }
    done.push(e.clone());
    let l1 = if ve_map.contains(&v1) { ve_map.find(&v1).clone() } else { Vec::new() };
    let i1 = nbedconnex(&l1);
    for cur in &l1 {
        let cur_e = cur;
        if i1 <= 2 {
            if is_vois(cur_e, vref, ve_map, done, prof, profmax) {
                return true;
            }
        } else if is_vois(cur_e, vref, ve_map, done, prof + 1, profmax) {
            return true;
        }
    }
    let l2 = if ve_map.contains(&v2) { ve_map.find(&v2).clone() } else { Vec::new() };
    for cur in &l2 {
        let cur_e = cur;
        if i1 <= 2 {
            if is_vois(cur_e, vref, ve_map, done, prof, profmax) {
                return true;
            }
        } else if is_vois(cur_e, vref, ve_map, done, prof + 1, profmax) {
            return true;
        }
    }
    false
}

// =========================================================================
// OCCT ChFi3d_Builder_6.cxx L195-205 — IsObst.
// =========================================================================
pub fn is_obst(cp: &ChFiDS_CommonPoint, vref: &Shape, ve_map: &ChFiDSMap) -> bool {
    if !cp.is_on_arc() {
        return false;
    }
    let e = cp.arc();
    let mut done: Vec<Shape> = Vec::new();
    let prof = 4;
    !is_vois(&e, vref, ve_map, &mut done, 0, prof)
}

// =========================================================================
// OCCT ChFi3d_Builder_6.cxx L209-302 — CompParam.
// =========================================================================
pub fn comp_param(carc: &Curve2d, ctg: &Curve2d, parc: &mut f64, ptg: &mut f64, prefarc: f64, preftg: f64) {
    let mut found = false;
    // (1) It is checked if the provided parameters are good
    //     if pcurves have the same parameters as the spine.
    let point = carc.point_at(prefarc);
    let distini = point.distance(ctg.point_at(preftg));
    if distini <= rcad_kernel::core::precision::PCONFUSION {
        *parc = prefarc;
        *ptg = preftg;
        found = true;
    } else {
        // (2) Intersection — Geom2dInt_GInter bridged by AnaIntersection2d
        // (architecture note in the file header).
        let intersection = geom2d_int_g_inter(ctg, carc);
        let mut dist = rcad_kernel::core::precision::INFINITE_VALUE;
        if intersection.is_done() && !intersection.is_empty() {
            // OCCT reads NbSegments() here; the rcad analytic bridge
            // carries intersection points only.
            let nbpt = intersection.nb_points();
            for i in 0..nbpt {
                let int2d = intersection.point(i);
                let p1 = int2d.param_on_first();
                let p2 = int2d.param_on_second();
                if (prefarc - p2).abs() < dist {
                    *ptg = p1;
                    *parc = p2;
                    dist = (prefarc - p2).abs();
                    found = true;
                }
            }
        }
    }

    if !found {
        // (3) Projection...
        *parc = prefarc;
        // OCCT: Geom2dAPI_ProjectPointOnCurve projector(point, Ctg).
        let (nbpoints, lower_distance, lower_distance_parameter) = geom2d_api_project_point_on_curve(point, ctg);
        if nbpoints == 0 {
            // This happens in some cases when there is a vertex
            // at the end of spine...
            *ptg = preftg;
        } else {
            // It is checked if everything was calculated correctly (EDC402 C2)
            if lower_distance < distini {
                *ptg = lower_distance_parameter;
            } else {
                *ptg = preftg;
            }
        }
    }
}

/// OCCT Geom2dInt_GInter::Perform(C1, C2, tol1, tol2) bridged to the rcad
/// analytic intersector (boundary note in the file header).  C1 keeps the
/// OCCT first-curve position (ParamOnFirst); non-analytic inputs yield the
/// not-done state (the caller falls to the OCCT !found branch).
fn geom2d_int_g_inter(c1: &Curve2d, c2: &Curve2d) -> rcad_kernel::base::int_ana2d::AnaIntersection2d {
    use rcad_kernel::base::int_ana2d::{AnaIntersection2d, Conic2d};
    fn conic_of(c: &Curve2d) -> Option<Conic2d> {
        match c {
            Curve2d::Line(l) => Some(Conic2d::from_line(l)),
            Curve2d::Circle(ci) => Some(Conic2d::from_circle(ci)),
            Curve2d::Ellipse(e) => Some(Conic2d::from_ellipse(e)),
            Curve2d::Parabola(p) => Some(Conic2d::from_parabola(p)),
            Curve2d::Hyperbola(h) => Some(Conic2d::from_hyperbola(h)),
            _ => None,
        }
    }
    let mut inter = AnaIntersection2d::new();
    match (c1, c2) {
        (Curve2d::Line(l1), Curve2d::Line(l2)) => inter.perform_lin_lin(l1, l2),
        (Curve2d::Line(l), Curve2d::Circle(ci)) => inter.perform_lin_circ(l, ci),
        (Curve2d::Line(l), other) => {
            if let Some(k) = conic_of(other) {
                inter.perform_lin_conic(l, &k);
            }
        }
        (Curve2d::Circle(ci), Curve2d::Line(l)) => {
            inter.perform_circ_conic(ci, &Conic2d::from_line(l));
        }
        (Curve2d::Circle(c1), Curve2d::Circle(c2)) => inter.perform_circ_circ(c1, c2),
        (Curve2d::Circle(ci), other) => {
            if let Some(k) = conic_of(other) {
                inter.perform_circ_conic(ci, &k);
            }
        }
        (Curve2d::Ellipse(e), other) => {
            if let Some(k) = conic_of(other) {
                inter.perform_ellipse_conic(e, &k);
            }
        }
        (Curve2d::Parabola(p), other) => {
            if let Some(k) = conic_of(other) {
                inter.perform_parabola_conic(p, &k);
            }
        }
        (Curve2d::Hyperbola(h), other) => {
            if let Some(k) = conic_of(other) {
                inter.perform_hyperbola_conic(h, &k);
            }
        }
        _ => {
            // OCCT handles BSpline pcurves here; the rcad analytic bridge
            // does not — stays not-done (falls to the projection branch).
        }
    }
    inter
}

/// OCCT Geom2dAPI_ProjectPointOnCurve(point, C2d) — (NbPoints(),
/// LowerDistance(), LowerDistanceParameter()).  rcad boundary: sampled
/// projection over the curve domain (pending AppDef/ExtremaPKG port).
fn geom2d_api_project_point_on_curve(p: DVec2, c: &Curve2d) -> (usize, f64, f64) {
    let [f0, l0] = c.default_domain();
    let n = 64usize;
    let mut best = f0;
    let mut bestd = f64::MAX;
    for i in 0..=n {
        let t = f0 + (l0 - f0) * (i as f64) / (n as f64);
        let d = p.distance(c.point_at(t));
        if d < bestd {
            bestd = d;
            best = t;
        }
    }
    (1, bestd, best)
}

// =========================================================================
// OCCT ChFi3d_Builder_6.cxx L304-339 — CompBlendPoint (create BlendPoint
// corresponding to a tangency on Vertex; pmn 15/10/1997: returns false
// if there is no pcurve).
// =========================================================================
pub fn comp_blend_point(
    brep: &topods::BRep,
    v: &Shape,
    e: &Shape,
    w: f64,
    f1: &Shape,
    f2: &Shape,
    bp: &mut BlendPoint,
) -> bool {
    let p3d = brep.vertex_position(v);
    // OCCT: param = BRep_Tool::Parameter(V, E, F1) — the vertex parameter
    // on the edge pcurve of F1 (rcad fallback: builder_0 resolver).
    let mut param = match brep.parameter_on_edge(v, e, f1) {
        Some(p) => p,
        None => brep_tool_parameter(brep, v, e),
    };
    let Some((pc1, _, _)) = brep.curve_on_surface(e, f1) else {
        return false;
    };
    let p1 = pc1.point_at(param);
    param = match brep.parameter_on_edge(v, e, f2) {
        Some(p) => p,
        None => brep_tool_parameter(brep, v, e),
    };
    let Some((pc2, _, _)) = brep.curve_on_surface(e, f2) else {
        return false;
    };
    let p2 = pc2.point_at(param);
    bp.set_value_on_2_surfaces(p3d, p3d, w, p1.x, p1.y, p2.x, p2.y);
    true
}

// =========================================================================
// OCCT ChFi3d_Builder_6.cxx L341-387 — UpdateLine (updates extremities
// after a partial invalidation).
// =========================================================================
pub fn update_line(line: &mut BRepBlendLine, isfirst: bool) {
    if isfirst {
        let bp = line.point(1).clone();
        let tguide = bp.parameter();
        if line.start_point_on_first().parameter_on_guide() < tguide {
            let mut be = BRepBlendExtremity::new();
            let (u, v) = bp.parameters_on_s1();
            be.set_value(
                bp.point_on_s1(),
                u,
                v,
                bp.parameter(),
                rcad_kernel::core::precision::CONFUSION,
            );
            let second = line.start_point_on_second().clone();
            line.set_start_points(&be, &second);
        }
        if line.start_point_on_second().parameter_on_guide() < tguide {
            let mut be = BRepBlendExtremity::new();
            let (u, v) = bp.parameters_on_s2();
            be.set_value(
                bp.point_on_s2(),
                u,
                v,
                bp.parameter(),
                rcad_kernel::core::precision::CONFUSION,
            );
            let first = line.start_point_on_first().clone();
            line.set_start_points(&first, &be);
        }
    } else {
        let nbp = line.nb_points();
        let bp = line.point(nbp).clone();
        let tguide = bp.parameter();
        if line.end_point_on_first().parameter_on_guide() > tguide {
            let mut be = BRepBlendExtremity::new();
            let (u, v) = bp.parameters_on_s1();
            be.set_value(
                bp.point_on_s1(),
                u,
                v,
                bp.parameter(),
                rcad_kernel::core::precision::CONFUSION,
            );
            let second = line.end_point_on_second().clone();
            line.set_end_points(&be, &second);
        }
        if line.end_point_on_second().parameter_on_guide() > tguide {
            let mut be = BRepBlendExtremity::new();
            let (u, v) = bp.parameters_on_s2();
            be.set_value(
                bp.point_on_s2(),
                u,
                v,
                bp.parameter(),
                rcad_kernel::core::precision::CONFUSION,
            );
            let first = line.end_point_on_first().clone();
            line.set_end_points(&first, &be);
        }
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2174-2191 — ChFi3d_TrsfTrans (consumed by
// CompleteData / StoreData / ComputeData / ChFi3d_FilCommonPoint).
// =========================================================================
pub fn chfi3d_trsf_trans(t1: IntSurfTypeTrans) -> Orientation {
    match t1 {
        IntSurfTypeTrans::In => Orientation::Forward,
        IntSurfTypeTrans::Out => Orientation::Reversed,
        _ => Orientation::Internal,
    }
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L2193-2320 — ChFi3d_FilCommonPoint (loading
// of the common point; management of the case when it happens on an
// already existing vertex).
// =========================================================================
pub fn chfi3d_fil_common_point(
    brep: &topods::BRep,
    sp: &BRepBlendExtremity,
    trans_line: IntSurfTypeTrans,
    start: bool,
    cp: &mut ChFiDS_CommonPoint,
    tol: f64,
) {
    let mut maxtol = tol.max(cp.tolerance());

    cp.set_point(sp.value()); // One starts with the point and the vector
    if sp.has_tangent() {
        if start {
            cp.set_vector(-sp.tangent()); // The tangent is oriented to the exit
        } else {
            cp.set_vector(sp.tangent());
        }
    }

    cp.set_parameter(sp.parameter_on_guide()); // and the parameter of the spine

    if sp.is_vertex() {
        // the Vertex is loaded if required (inside of a face)
        // OCCT: V = down_cast<BRepTopAdaptor_HVertex>(SP.Vertex())->Vertex()
        // — pending boundary: the rcad HVertex carries no TopoDS_Vertex
        // payload (topalgo/adaptor3d/hvertex.rs), the load is deferred.
        if let Some(v) = hvertex_vertex(sp.vertex()) {
            cp.set_vertex(v.clone());
            let dist = sp.value().distance(brep.vertex_position(&v));
            // modified by jgv, 18.09.02 for OCC571
            maxtol = dist.max(maxtol);
            cp.set_point(brep.vertex_position(&v));
        }
        // the sequence of arcs the information is known by the vertex
        // (ancestor); in this case the transitions are not computed, it is
        // done by this program
    }

    if sp.nb_point_on_rst() != 0 {
        // An arc, and/or a vertex is loaded
        let pr: &BRepBlendPointOnRst = sp.point_on_rst(1);
        // OCCT: Harc = down_cast<BRepAdaptor_Curve2d>(PR.Arc());
        //       E = Harc->Edge() — pending boundary: the rcad PointOnRst
        // carries the pcurve only, the TopoDS_Edge is deferred.
        if let Some(e) = point_on_rst_edge(pr) {
            let distf;
            let distl;
            let le_param_amoi;
            let v = topexp_vertices(&e);
            let v0 = &v.0;
            let v1 = &v.1;

            distf = sp.value().distance(brep.vertex_position(v0));
            distl = sp.value().distance(brep.vertex_position(v1));
            let index_min;
            let dist;
            if distf < distl {
                index_min = 0;
                dist = distf;
            } else {
                index_min = 1;
                dist = distl;
            }
            let v_index_min = if index_min == 0 { v0 } else { v1 };

            if dist <= maxtol + brep.vertex_tolerance(v_index_min) {
                // a preexisting vertex has been met
                cp.set_vertex(v_index_min.clone()); // the old vertex is loaded
                cp.set_point(brep.vertex_position(v_index_min));
                maxtol = brep.vertex_tolerance(v_index_min).max(maxtol);
                // modified by jgv, 18.09.02 for OCC571
                maxtol = dist.max(maxtol);
                le_param_amoi = brep_tool_parameter(brep, v_index_min, &e);
            } else {
                // Creation of an arc only
                maxtol = brep.tolerance(&e).max(maxtol);
                maxtol = sp.tolerance().max(maxtol);
                le_param_amoi = pr.parameter_on_arc();
            }

            // Definition of the arc
            let tr;
            let or = e.orientation;
            if start {
                tr = topabs_reverse(topabs_compose(chfi3d_trsf_trans(trans_line), or));
            } else {
                tr = topabs_compose(chfi3d_trsf_trans(trans_line), or);
            }
            cp.set_arc(maxtol, e, le_param_amoi, tr);
        }
    }
    cp.set_tolerance(maxtol); // Finally, the tolerance.
}

/// OCCT `SP.Vertex()->Vertex()` — pending boundary: the rcad HVertex
/// (topalgo/adaptor3d/hvertex.rs) carries no TopoDS_Vertex payload yet.
fn hvertex_vertex(_hv: &HVertex) -> Option<Shape> {
    None
}

/// OCCT `Harc->Edge()` — pending boundary: the rcad BRepBlendPointOnRst
/// carries the restriction pcurve only; the owning TopoDS_Edge is deferred.
fn point_on_rst_edge(_pr: &BRepBlendPointOnRst) -> Option<Shape> {
    None
}

// =========================================================================
// OCCT ChFi3d_Builder_0.cxx L1565-1606 — ChFi3d_CheckSameParameter
// (pending boundary: the identity keeps the input tolerance, same
// convention as chfi3d_same_parameter in chfi3d_builder_0).  Returns true
// (the reached tolerance is within tol3d); `tolcheck_out` receives the
// reached tolerance.
// =========================================================================
pub fn chfi3d_check_same_parameter(
    _c3d: Option<&Curve3>,
    _pcurv: &Curve2d,
    _s: &BRepAdaptorSurface,
    tol3d: f64,
    tolcheck_out: &mut f64,
) -> bool {
    *tolcheck_out = tol3d;
    true
}

// =========================================================================
// Local geometry helpers (OCCT Geom_* equivalents on rcad Surface3).
// =========================================================================

/// OCCT Geom_Surface::Bounds — the natural (underlying) bounds.
fn geom_surface_bounds(s: &rcad_kernel::geom::Surface3) -> (f64, f64, f64, f64) {
    let b = s.default_domain();
    (b[0], b[1], b[2], b[3])
}

/// OCCT Geom_Surface::VIso(V) for the surfaces reaching the Surfcoin path
/// (plane; cylinder anchor kept from chfi3d_builder_0::cylinder_v_iso).
fn surface_v_iso(s: &rcad_kernel::geom::Surface3, v: f64) -> Option<Curve3> {
    use rcad_kernel::geom::Surface3;
    match s {
        Surface3::Plane(p) => {
            // OCCT Geom_Plane::VIso — the u-line at v.
            let (pt, du, _) = p.derivatives(0.0, v);
            Some(Curve3::Line(rcad_kernel::geom::Line3::new(pt, du.normalize())))
        }
        Surface3::Cylinder(_) => {
            // OCCT Geom_CylindricalSurface::VIso — circle at height v.
            let circle = super::chfi3d_builder_0::cylinder_v_iso(
                cyl_origin(s),
                cyl_xdir(s),
                cyl_axis(s),
                cyl_radius(s),
                v,
            );
            Some(Curve3::Circle(circle))
        }
        _ => None,
    }
}

fn cyl_origin(s: &rcad_kernel::geom::Surface3) -> DVec3 {
    match s {
        rcad_kernel::geom::Surface3::Cylinder(c) => c.origin,
        _ => DVec3::ZERO,
    }
}

fn cyl_xdir(s: &rcad_kernel::geom::Surface3) -> DVec3 {
    match s {
        rcad_kernel::geom::Surface3::Cylinder(c) => c.ref_dir,
        _ => DVec3::ZERO,
    }
}

fn cyl_axis(s: &rcad_kernel::geom::Surface3) -> DVec3 {
    match s {
        rcad_kernel::geom::Surface3::Cylinder(c) => c.axis,
        _ => DVec3::ZERO,
    }
}

fn cyl_radius(s: &rcad_kernel::geom::Surface3) -> f64 {
    match s {
        rcad_kernel::geom::Surface3::Cylinder(c) => c.radius,
        _ => 0.0,
    }
}

/// OCCT Geom_BSplineSurface::UIso(U) — pending the AppSurf batch (the
/// BSpline approximation surfaces are produced there); returns None which
/// the DS curve keeps null.
fn surface_u_iso_bspline(_s: &rcad_kernel::geom::BSplineSurface, _u: f64) -> Option<Curve3> {
    None
}

/// OCCT Geom_BSplineSurface ctor from (poles, weights, UKnots, VKnots,
/// UMults, VMults, UDeg, VDeg) mapped onto the rcad expanded-knots
/// BSplineSurface (knot expansion per BSplCLib::BuildKnots).
#[allow(clippy::too_many_arguments)]
fn make_bspline_surface(
    poles: &[Vec<DVec3>],
    weights: &[Vec<f64>],
    u_knots: &[f64],
    v_knots: &[f64],
    u_mults: &[i32],
    v_mults: &[i32],
    u_degree: i32,
    v_degree: i32,
) -> rcad_kernel::geom::BSplineSurface {
    fn expand(knots: &[f64], mults: &[i32]) -> Vec<f64> {
        let mut out = Vec::new();
        for (k, m) in knots.iter().zip(mults.iter()) {
            for _ in 0..*m {
                out.push(*k);
            }
        }
        out
    }
    rcad_kernel::geom::BSplineSurface {
        degree_u: u_degree.max(0) as usize,
        degree_v: v_degree.max(0) as usize,
        knots_u: expand(u_knots, u_mults),
        knots_v: expand(v_knots, v_mults),
        control_points: poles.to_vec(),
        weights: weights.to_vec(),
    }
}

/// OCCT Geom2d_BSplineCurve ctor from (poles, knots, mults, degree) on the
/// rcad expanded-knots BSplineCurve2.
fn make_bspline_curve2d(poles: &[DVec2], knots: &[f64], mults: &[i32], degree: i32) -> Curve2d {
    let mut expanded = Vec::new();
    for (k, m) in knots.iter().zip(mults.iter()) {
        for _ in 0..*m {
            expanded.push(*k);
        }
    }
    Curve2d::BSpline(rcad_kernel::geom::BSplineCurve2 {
        degree: degree.max(0) as usize,
        knots: expanded,
        control_points: poles.to_vec(),
        weights: vec![1.0; poles.len()],
    })
}

/// OCCT Geom2d_BSplineCurve::FirstParameter / LastParameter — the first /
/// last distinct knot.
fn bspline2d_parameter_bounds(c: &Curve2d) -> (f64, f64) {
    match c {
        Curve2d::BSpline(b) => (
            b.knots[b.degree],
            b.knots[b.knots.len() - 1 - b.degree],
        ),
        other => {
            let d = other.default_domain();
            (d[0], d[1])
        }
    }
}

impl ChFi3dBuilder {
    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L395-505 — CompleteData (Surfcoin variant:
    // calculates curves and CommonPoints from the data calculated by
    // filling).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn complete_data_surfcoin(
        &mut self,
        data: &mut ChFiDSSurfData,
        surfcoin: &rcad_kernel::geom::Surface3,
        s1: &BRepAdaptorSurface,
        pc1: Option<&Curve2d>,
        s2: &BRepAdaptorSurface,
        pc2: Option<&Curve2d>,
        or: Orientation,
        on1: bool,
        gd1: bool,
        gd2: bool,
        gf1: bool,
        gf2: bool,
    ) -> bool {
        let dstr = self.my_ds.as_mut().expect("DS");
        let surf_index = dstr.add_surface(TopOpeBRepDSSurface::new(surfcoin.clone(), self.tolapp3d));
        data.change_surf(surf_index);

        let (u_first, u_last, v_first, v_last) = geom_surface_bounds(surfcoin);
        if !gd1 {
            data.change_vertex_first_on_s1()
                .set_point(surfcoin.point_at(u_first, v_first));
        }
        if !gd2 {
            data.change_vertex_first_on_s2()
                .set_point(surfcoin.point_at(u_first, v_last));
        }
        if !gf1 {
            data.change_vertex_last_on_s1()
                .set_point(surfcoin.point_at(u_last, v_first));
        }
        if !gf2 {
            data.change_vertex_last_on_s2()
                .set_point(surfcoin.point_at(u_last, v_last));
        }

        // calculate curves side S1
        let mut crv3d1: Option<Curve3> = None;
        if pc1.is_some() {
            crv3d1 = surface_v_iso(surfcoin, v_first);
        }
        let pd1 = DVec2::new(u_first, v_first);
        let pf1 = DVec2::new(u_last, v_first);
        let lfil1 = rcad_kernel::geom::Line2d {
            origin: pd1,
            direction: (pf1 - pd1).normalize(),
        };
        let pcurve_on_surf = Curve2d::Line(lfil1);
        let mut tra1 = Orientation::Forward;
        let mut orsurf = or;
        let w = 0.5 * (u_first + u_last);
        let mut tolreached = 1.0e-5;
        let mut c2dtrim: Option<Curve2d> = None;
        if let Some(pc1) = pc1 {
            // OCCT: c2dtrim = new Geom2d_TrimmedCurve(PC1, UFirst, ULast)
            let mut c2d = Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                curve: Box::new(pc1.clone()),
                t_min: u_first,
                t_max: u_last,
            });
            if let Some(c3d) = &crv3d1 {
                chfi3d_same_parameter(c3d, &mut c2d, &s1.surface, self.tolapp3d, &mut tolreached);
            }
            let (x, y) = {
                let p = c2d.point_at(w);
                (p.x, p.y)
            };
            let (_, du1, dv1) = s1.surface.derivatives(x, y);
            let nf = du1.cross(dv1);
            let (_, du2, dv2) = surfcoin.derivatives(w, v_first);
            let ns = du2.cross(dv2);
            if nf.dot(ns) > 0.0 {
                tra1 = Orientation::Reversed;
            } else if on1 {
                orsurf = topabs_reverse(orsurf);
            }
            c2dtrim = Some(c2d);
        }
        let index1_of_curve = dstr.add_curve(TopOpeBRepDSCurve::new(crv3d1.clone(), tolreached));
        {
            let fint1 = data.change_interference_on_s1();
            fint1.set_first_parameter(u_first);
            fint1.set_last_parameter(u_last);
            fint1.set_interference(index1_of_curve, tra1, c2dtrim.clone(), Some(pcurve_on_surf.clone()));
        }
        // calculate curves side S2
        let mut crv3d2: Option<Curve3> = None;
        if pc2.is_some() {
            crv3d2 = surface_v_iso(surfcoin, v_last);
        }
        let pd2 = DVec2::new(u_first, v_last);
        let pf2 = DVec2::new(u_last, v_last);
        let lfil2 = rcad_kernel::geom::Line2d {
            origin: pd2,
            direction: (pf2 - pd2).normalize(),
        };
        let pcurve_on_surf2 = Curve2d::Line(lfil2);
        let mut tra2 = Orientation::Forward;
        if let Some(pc2) = pc2 {
            let mut c2d = Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 {
                curve: Box::new(pc2.clone()),
                t_min: u_first,
                t_max: u_last,
            });
            if let Some(c3d) = &crv3d2 {
                chfi3d_same_parameter(c3d, &mut c2d, &s2.surface, self.tolapp3d, &mut tolreached);
            }
            let (x, y) = {
                let p = c2d.point_at(w);
                (p.x, p.y)
            };
            let (_, du1, dv1) = s2.surface.derivatives(x, y);
            let np = du1.cross(dv1);
            let (_, du2, dv2) = surfcoin.derivatives(w, v_last);
            let ns = du2.cross(dv2);
            if np.dot(ns) < 0.0 {
                tra2 = Orientation::Reversed;
                if !on1 {
                    orsurf = topabs_reverse(orsurf);
                }
            }
            c2dtrim = Some(c2d);
        }
        let index2_of_curve = dstr.add_curve(TopOpeBRepDSCurve::new(crv3d2, tolreached));
        {
            let fint2 = data.change_interference_on_s2();
            fint2.set_first_parameter(u_first);
            fint2.set_last_parameter(u_last);
            fint2.set_interference(index2_of_curve, tra2, c2dtrim, Some(pcurve_on_surf2));
        }
        *data.change_orientation() = orsurf;
        true
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L516-551 — CompleteData (Blend_Function
    // overload: calculates the surface of curves and eventually
    // CommonPoints from the data calculated in ComputeData; use of F(t)).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn complete_data_function(
        &mut self,
        data: &mut ChFiDSSurfData,
        func: &mut dyn BlendFunction,
        lin: &BRepBlendLine,
        s1: &BRepAdaptorSurface,
        s2: Option<&BRepAdaptorSurface>,
        or1: Orientation,
        gd1: bool,
        gd2: bool,
        gf1: bool,
        gf2: bool,
        reversed: bool,
    ) -> bool {
        // OCCT: TheFunc = new BRepBlend_AppFunc(lin, Func, tolapp3d, 1.e-5)
        let the_func = BRepBlendAppFunc::new(lin, func, self.tolapp3d, 1.0e-5);

        let degmax = 20;
        let segmax = 5000;
        let approx = BRepBlendAppSurface::new(
            &the_func,
            lin.point(1).parameter(),
            lin.point(lin.nb_points()).parameter(),
            self.tolapp3d,
            1.0e-5,           // tolapp2d, tolerance max
            self.tolappangle, // Contact G1
            self.my_conti,
            degmax,
            segmax,
        );
        if !approx.is_done() {
            return false;
        }
        self.store_data(
            data, &approx, lin, s1, s2, or1, gd1, gd2, gf1, gf2, reversed,
        )
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L559-588 — CompleteData (Blend_SurfRstFunction
    // overload; jlr le 28/07/97 branchement F(t)).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn complete_data_surf_rst(
        &mut self,
        data: &mut ChFiDSSurfData,
        func: &mut dyn BlendSurfRstFunction,
        lin: &BRepBlendLine,
        s1: &BRepAdaptorSurface,
        s2: Option<&BRepAdaptorSurface>,
        or: Orientation,
        reversed: bool,
    ) -> bool {
        // OCCT: TheFunc = new BRepBlend_AppFuncRst(lin, Func, tolapp3d, 1.e-5)
        let the_func = BRepBlendAppFuncRst::new(lin, func, self.tolapp3d, 1.0e-5);
        let approx = BRepBlendAppSurface::new_defaulted(
            &the_func,
            lin.point(1).parameter(),
            lin.point(lin.nb_points()).parameter(),
            self.tolapp3d,
            1.0e-5,           // tolapp2d, tolerance max
            self.tolappangle, // Contact G1
            self.my_conti,
        );
        if !approx.is_done() {
            return false;
        }
        self.store_data(data, &approx, lin, s1, s2, or, false, false, false, false, reversed)
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L596-624 — CompleteData (Blend_RstRstFunction
    // overload; jlr le 28/07/97 branchement F(t)).
    // =====================================================================
    pub fn complete_data_rst_rst(
        &mut self,
        data: &mut ChFiDSSurfData,
        func: &mut dyn BlendRstRstFunction,
        lin: &BRepBlendLine,
        s1: &BRepAdaptorSurface,
        s2: Option<&BRepAdaptorSurface>,
        or: Orientation,
    ) -> bool {
        // OCCT: TheFunc = new BRepBlend_AppFuncRstRst(lin, Func, tolapp3d, 1.e-5)
        let the_func = BRepBlendAppFuncRstRst::new(lin, func, self.tolapp3d, 1.0e-5);
        let approx = BRepBlendAppSurface::new_defaulted(
            &the_func,
            lin.point(1).parameter(),
            lin.point(lin.nb_points()).parameter(),
            self.tolapp3d,
            1.0e-5,           // tolapp2d, tolerance max
            self.tolappangle, // Contact G1
            self.my_conti,
        );
        if !approx.is_done() {
            return false;
        }
        self.store_data(data, &approx, lin, s1, s2, or, false, false, false, false, false)
    }

    // =====================================================================
    // OCCT ChFi3d_Builder_6.cxx L631-1004 — StoreData (copy of an
    // approximation result in SurfData).
    // =====================================================================
    #[allow(clippy::too_many_arguments)]
    pub fn store_data(
        &mut self,
        data: &mut ChFiDSSurfData,
        approx: &dyn AppBlendApprox,
        lin: &BRepBlendLine,
        s1: &BRepAdaptorSurface,
        s2: Option<&BRepAdaptorSurface>,
        or1: Orientation,
        gd1: bool,
        gd2: bool,
        gf1: bool,
        gf2: bool,
        reversed: bool,
    ) -> bool {
        let mut tolget3d: f64 = 0.0;
        let mut tolget2d: f64 = 0.0;
        approx.tol_reached(&mut tolget3d, &mut tolget2d);
        let tolaux = approx.tol_curve_on_surf(1);
        let mut tol_c1 = tolget3d + tolaux;
        let mut tol_c2 = 0.0;
        if s2.is_some() {
            let tolaux = approx.tol_curve_on_surf(2);
            tol_c2 = tolget3d + tolaux;
        }

        // By default parametric space is created using a square surface
        // to be able to parameterize in U by # R*teta // a revoir lbo 29/08/97
        let ku = approx.surf_u_knots();
        let kv = approx.surf_v_knots();
        let larg = kv[kv.len() - 1] - kv[0];
        let mut kku = ku.to_vec();
        rcad_kernel::math::bspl_lib::reparametrize(0.0, larg, &mut kku);
        let mut surf = make_bspline_surface(
            approx.surf_poles(),
            approx.surf_weights(),
            &kku,
            kv,
            approx.surf_u_mults(),
            approx.surf_v_mults(),
            approx.u_degree(),
            approx.v_degree(),
        );
        // extension of the surface

        let length1 = data.first_extension_value();
        let length2 = data.last_extension_value();

        let mut ext1 = false;
        let mut ext2 = false;
        let eps = tolget3d.max(2.0 * rcad_kernel::core::precision::CONFUSION);
        if length1 > eps {
            let p11 = surf.control_points[0][0];
            let p21 = surf.control_points[surf.control_points.len() - 1][0];
            if p11.distance(p21) > eps {
                // to avoid extending surface with singular boundary
                // OCCT: GeomLib::ExtendSurfByLength(aBndSurf, length1, 1, false, false)
                // — pending GeomLib surface extension (chfi3d_builder_c1 note).
                ext1 = true;
            }
        }
        if length2 > eps {
            let nv = surf.control_points[0].len();
            let p12 = surf.control_points[0][nv - 1];
            let p22 = surf.control_points[surf.control_points.len() - 1][nv - 1];
            if p12.distance(p22) > eps {
                // to avoid extending surface with singular boundary
                // OCCT: GeomLib::ExtendSurfByLength(aBndSurf, length2, 1, false, true)
                // — pending GeomLib surface extension (chfi3d_builder_c1 note).
                ext2 = true;
            }
        }
        // OCCT L710: Surf = down_cast<Geom_BSplineSurface>(aBndSurf) — the
        // extension keeps the BSpline surface; the pending GeomLib branch
        // above leaves `surf` unchanged.

        // Correction of surface on extremities
        if !ext1 {
            let p11 = lin.start_point_on_first().value();
            let p21 = lin.start_point_on_second().value();
            surf.control_points[0][0] = p11;
            let nu = surf.control_points.len();
            surf.control_points[nu - 1][0] = p21;
        }
        if !ext2 {
            let p12 = lin.end_point_on_first().value();
            let p22 = lin.end_point_on_second().value();
            let nv = surf.control_points[0].len();
            surf.control_points[0][nv - 1] = p12;
            let nu = surf.control_points.len();
            surf.control_points[nu - 1][nv - 1] = p22;
        }

        let surf3 = rcad_kernel::geom::Surface3::BSpline(surf.clone());
        let surf_index = self.my_ds.as_mut().expect("DS").add_surface(TopOpeBRepDSSurface::new(surf3.clone(), tolget3d));
        data.change_surf(surf_index);

        // OCCT L732-733: Surf->Bounds(...)
        let nu = surf.knots_u.len() - 1 - surf.degree_u;
        let nv = surf.knots_v.len() - 1 - surf.degree_v;
        let (u_first, u_last, v_first, v_last) = (
            surf.knots_u[surf.degree_u],
            surf.knots_u[nu],
            surf.knots_v[surf.degree_v],
            surf.knots_v[nv],
        );

        let mut uon1 = u_first;
        let mut uon2 = u_last;
        let mut ion1 = 1;
        let mut ion2 = 2;
        if reversed {
            uon1 = u_last;
            uon2 = u_first;
            ion1 = 2;
            ion2 = 1;
        }

        // The SurfData is filled in what concerns S1,
        // OCCT: Crv3d1 = Surf->UIso(Uon1)
        let crv3d1 = surface_u_iso_bspline(&surf, uon1);
        let pori1 = DVec2::new(uon1, 0.0);
        let lfil1 = rcad_kernel::geom::Line2d {
            origin: pori1,
            direction: DVec2::new(0.0, 1.0), // gp::DY2d()
        };
        let pcurve_on_surf1 = Curve2d::Line(lfil1);
        // OCCT: PCurveOnFace = new Geom2d_BSplineCurve(approx.Curve2dPoles(ion1), ...)
        let pcurve_on_face1 = make_bspline_curve2d(
            approx.curve2d_poles(ion1),
            approx.curves2d_knots(),
            approx.curves2d_mults(),
            approx.curves2d_degree(),
        );
        let (par1, par2) = bspline2d_parameter_bounds(&pcurve_on_face1);

        // OCCT: chc.Load(Crv3d1, par1, par2) + ChFi3d_CheckSameParameter(...)
        // (checkcurve pending boundary above)
        // OCCT leaves TolChk uninitialized; seeded 0.0 for definite-init
        // (the callee assigns on the true path).
        let mut tolcheck = 0.0;
        if !chfi3d_check_same_parameter(crv3d1.as_ref(), &pcurve_on_face1, s1, tol_c1, &mut tolcheck) {
            tol_c1 = tolcheck;
        }
        let index1_of_curve = self.my_ds.as_mut().expect("DS").add_curve(TopOpeBRepDSCurve::new(crv3d1, tol_c1));

        let mut uarc = 0.0;
        let mut utg = 0.0;
        let mut pppdeb;
        let mut pppfin;
        if gd1 {
            // OCCT: forwfac = BS1->Face(); forwfac.Orientation(TopAbs_FORWARD)
            let forwfac = self.forward_face(&s1.face);
            let brc_arc = data.vertex_first_on_s1().arc().clone();
            let carg = self.arc_pcurve_on_face(&brc_arc, &forwfac);
            let ca_arc = self.arc_3d_curve(&brc_arc);
            let v = data.change_vertex_first_on_s1();
            comp_param(
                &carg,
                &pcurve_on_face1,
                &mut uarc,
                &mut utg,
                v.parameter_on_arc(),
                v.parameter(),
            );
            let brep = self.my_brep.clone();
            tolcheck = ca_arc.point_at(uarc).distance(v.point());
            let arc = v.arc().clone();
            let tarc = v.transition_on_arc();
            v.set_arc(tol_c1 + tolcheck, arc, uarc, tarc);
            pppdeb = utg;
        } else {
            pppdeb = v_first;
        }
        if gf1 {
            let forwfac = self.forward_face(&s1.face);
            let brc_arc = data.vertex_last_on_s1().arc().clone();
            let carg = self.arc_pcurve_on_face(&brc_arc, &forwfac);
            let ca_arc = self.arc_3d_curve(&brc_arc);
            let v = data.change_vertex_last_on_s1();
            comp_param(
                &carg,
                &pcurve_on_face1,
                &mut uarc,
                &mut utg,
                v.parameter_on_arc(),
                v.parameter(),
            );
            let brep = self.my_brep.clone();
            tolcheck = ca_arc.point_at(uarc).distance(v.point());
            let arc = v.arc().clone();
            let tarc = v.transition_on_arc();
            v.set_arc(tol_c1 + tolcheck, arc, uarc, tarc);
            pppfin = utg;
        } else {
            pppfin = v_last;
        }
        let tr_on_1 = if reversed {
            chfi3d_trsf_trans(lin.transition_on_s2())
        } else {
            chfi3d_trsf_trans(lin.transition_on_s1())
        };
        {
            let fint1 = data.change_interference_on_s1();
            fint1.set_first_parameter(pppdeb);
            fint1.set_last_parameter(pppfin);
            fint1.set_interference(
                index1_of_curve,
                tr_on_1,
                Some(pcurve_on_face1.clone()),
                Some(pcurve_on_surf1),
            );
        }

        // SurfData is filled in what concerns S2,
        // OCCT: Crv3d2 = Surf->UIso(Uon2)
        let crv3d2 = surface_u_iso_bspline(&surf, uon2);
        let pori2 = DVec2::new(uon2, 0.0);
        let lfil2 = rcad_kernel::geom::Line2d {
            origin: pori2,
            direction: DVec2::new(0.0, 1.0),
        };
        let pcurve_on_surf2 = Curve2d::Line(lfil2);
        let pcurve_on_face2: Option<Curve2d> = if s2.is_some() {
            let pc = make_bspline_curve2d(
                approx.curve2d_poles(ion2),
                approx.curves2d_knots(),
                approx.curves2d_mults(),
                approx.curves2d_degree(),
            );
            if !chfi3d_check_same_parameter(crv3d2.as_ref(), &pc, s2.unwrap(), tol_c2, &mut tolcheck) {
                tol_c2 = tolcheck;
            }
            Some(pc)
        } else {
            None
        };
        let index2_of_curve = self.my_ds.as_mut().expect("DS").add_curve(TopOpeBRepDSCurve::new(crv3d2, tol_c2));
        if gd2 {
            let forwfac = self.forward_face(&s2.as_ref().expect("BS2").face);
            let brc_arc = data.vertex_first_on_s2().arc().clone();
            let carg = self.arc_pcurve_on_face(&brc_arc, &forwfac);
            let ca_arc = self.arc_3d_curve(&brc_arc);
            let v = data.change_vertex_first_on_s2();
            let pc_ref = pcurve_on_face2.clone().expect("PCurveOnFace");
            comp_param(
                &carg,
                &pc_ref,
                &mut uarc,
                &mut utg,
                v.parameter_on_arc(),
                v.parameter(),
            );
            let brep = self.my_brep.clone();
            tolcheck = ca_arc.point_at(uarc).distance(v.point());
            let arc = v.arc().clone();
            let tarc = v.transition_on_arc();
            v.set_arc(tol_c2 + tolcheck, arc, uarc, tarc);
            pppdeb = utg;
        } else {
            pppdeb = v_first;
        }
        if gf2 {
            let forwfac = self.forward_face(&s2.as_ref().expect("BS2").face);
            let brc_arc = data.vertex_last_on_s2().arc().clone();
            let carg = self.arc_pcurve_on_face(&brc_arc, &forwfac);
            let ca_arc = self.arc_3d_curve(&brc_arc);
            let v = data.change_vertex_last_on_s2();
            let pc_ref = pcurve_on_face2.clone().expect("PCurveOnFace");
            comp_param(
                &carg,
                &pc_ref,
                &mut uarc,
                &mut utg,
                v.parameter_on_arc(),
                v.parameter(),
            );
            let brep = self.my_brep.clone();
            tolcheck = ca_arc.point_at(uarc).distance(v.point());
            let arc = v.arc().clone();
            let tarc = v.transition_on_arc();
            v.set_arc(tol_c2 + tolcheck, arc, uarc, tarc);
            pppfin = utg;
        } else {
            pppfin = v_last;
        }
        if s2.is_some() {
            let tr_on_2 = if reversed {
                chfi3d_trsf_trans(lin.transition_on_s1())
            } else {
                chfi3d_trsf_trans(lin.transition_on_s2())
            };
            let fint2 = data.change_interference_on_s2();
            fint2.set_first_parameter(pppdeb);
            fint2.set_last_parameter(pppfin);
            fint2.set_interference(index2_of_curve, tr_on_2, pcurve_on_face2, Some(pcurve_on_surf2));
        } else {
            // OCCT: bidpc (null pcurve on face)
            let fint2 = data.change_interference_on_s2();
            fint2.set_first_parameter(pppdeb);
            fint2.set_last_parameter(pppfin);
            fint2.set_interference(index2_of_curve, Orientation::Forward, None, Some(pcurve_on_surf2));
        }

        // the orientation of the fillet in relation to the faces is evaluated,
        let sref: &BRepAdaptorSurface = if reversed {
            s2.expect("Sref = S2")
        } else {
            s1
        };
        let pcurve_on_face_ref: Curve2d = if reversed {
            data.interference_on_s2()
                .pcurve_on_face()
                .cloned()
                .expect("Fint2.PCurveOnFace()")
        } else {
            data.interference_on_s1()
                .pcurve_on_face()
                .cloned()
                .expect("Fint1.PCurveOnFace()")
        };

        //  Modified by skv - Wed Jun  9 17:16:26 2004 OCC5898 Begin
        let a_delta = v_last - v_first;
        let mut a_denom = 2.0f64;

        loop {
            let a_deltav = a_delta / a_denom;
            let a_param = v_first + a_deltav;
            let puv = pcurve_on_face_ref.point_at(a_param);
            let (p, du1_v, dv1_v) = sref.surface.derivatives(puv.x, puv.y);
            let mut du1 = du1_v.cross(dv1_v);

            if or1 == Orientation::Reversed {
                du1 = -du1;
            }

            let (p2, du2_v, dv2_v) = surf3.derivatives(u_first, a_param);
            let du2 = du2_v.cross(dv2_v);
            let _ = (p, p2);

            if du1.length() <= tolget3d || du2.length() <= tolget3d {
                a_denom += 1.0;

                if a_deltav.abs() <= tolget2d {
                    return false;
                }

                continue;
            }

            if du1.dot(du2) > 0.0 {
                *data.change_orientation() = Orientation::Forward;
            } else {
                *data.change_orientation() = Orientation::Reversed;
            }

            break;
        }
        //  Modified by skv - Wed Jun  9 17:16:26 2004 OCC5898 End

        if !gd1 {
            let brep = self.my_brep.clone();
            chfi3d_fil_common_point(
                &brep,
                lin.start_point_on_first(),
                lin.transition_on_s1(),
                true,
                data.change_vertex(true, ion1),
                tol_c1,
            );
        }
        if !gf1 {
            let brep = self.my_brep.clone();
            chfi3d_fil_common_point(
                &brep,
                lin.end_point_on_first(),
                lin.transition_on_s1(),
                false,
                data.change_vertex(false, ion1),
                tol_c1,
            );
        }
        if !gd2 && s2.is_some() {
            let brep = self.my_brep.clone();
            chfi3d_fil_common_point(
                &brep,
                lin.start_point_on_second(),
                lin.transition_on_s2(),
                true,
                data.change_vertex(true, ion2),
                tol_c2,
            );
        }
        if !gf2 && s2.is_some() {
            let brep = self.my_brep.clone();
            chfi3d_fil_common_point(
                &brep,
                lin.end_point_on_second(),
                lin.transition_on_s2(),
                false,
                data.change_vertex(false, ion2),
                tol_c2,
            );
        }
        // Parameters on ElSpine
        let nbp = lin.nb_points();
        data.set_first_spine_param(lin.point(1).parameter());
        data.set_last_spine_param(lin.point(nbp).parameter());
        true
    }

    /// OCCT `TopoDS_Face forwfac = BSx->Face(); forwfac.Orientation(TopAbs_FORWARD)`
    /// — the face with a forced FORWARD orientation.
    fn forward_face(&self, face: &Shape) -> Shape {
        let mut f = face.clone();
        f.orientation = Orientation::Forward;
        f
    }

    /// OCCT `BRepAdaptor_Curve2d brc; brc.Initialize(Arc, forwfac)` — the
    /// arc's pcurve on the face (BRep_Tool::CurveOnSurface via the pool).
    fn arc_pcurve_on_face(&self, arc: &Shape, forwfac: &Shape) -> Curve2d {
        self.my_brep
            .curve_on_surface(arc, forwfac)
            .map(|(pc, _, _)| pc)
            .expect("BRepAdaptor_Curve2d: no pcurve of arc on face")
    }

    /// OCCT `BRepAdaptor_Curve CArc; CArc.Initialize(V.Arc())` — the arc's
    /// 3d curve.
    fn arc_3d_curve(&self, arc: &Shape) -> Curve3 {
        self.my_brep
            .edge_curve_world(arc)
            .map(|(c, _)| c)
            .expect("BRepAdaptor_Curve: no 3d curve on arc")
    }
}
