//! OCCT `ShapeFix_Wire` — the file statics of `ShapeFix_Wire.cxx` and
//! `ShapeFix_Wire_1.cxx`, and the GAP re-hosts for the not-yet-landed
//! out-of-tranche dependencies.
//!
//! Statics (Wire.cxx): `UpdateEdgeUVPoints` (L1652-1657) / `TryNewPCurve`
//! (L2215-2252) / `howMuchPCurves` (L2262-2279) / `RemoveLoop(6-arg)`
//! (L2283-2523) / `RemoveLoop(E1,E2)` (L2527-2704) / `ComputeLocalDeviation`
//! (L2917-2962) / `TryBendingPCurve` (L3532-3613) / `CopyReversePcurves`
//! (L4108-4174).  Statics (Wire_1.cxx): `AdjustOnPeriodic3d` (L126-152) /
//! `AdjustOnPeriodic2d` (L865-891).
//!
//! GAP re-hosts (the iron rule: dependency + anchor + OCCT failure path):
//! - [`GeomConvertCompCurveToBSplineGap`] — OCCT
//!   `GeomConvert_CompCurveToBSplineCurve` (TKGeomBase, kernel-geom scope)
//!   reduced to `Add`/`BSplineCurve`; keeps OCCT's `!Add` failure branch
//!   (`Add` -> false), which makes the enclosing OCCT body return false.
//! - [`geom_api_to3d`] / [`geom_api_to2d`] — OCCT `GeomAPI::To3d`/`To2d`
//!   (the planar 2d<->3d conversion, TKGeomBase) reduced to the null-result
//!   path; the consumers keep OCCT's `IsNull` guards.
//! - `param_on_first` / `param_on_second` — OCCT
//!   `IntRes2d_IntersectionPoint::ParamOnFirst/Second`: the landed
//!   `ShapeAnalysis_Wire` self-intersection sampler keeps the OCCT `!IsDone`
//!   path (the output sequences stay empty — wire_checks.rs bridge #5), so
//!   the accessors are unreachable; they read the (x, y) parameter slots of
//!   the sampler element type for form.
//!
//! Retired by the W3 tranche 3 (the real tools landed 1:1 in
//! `shape_fix/split_tool.rs` and `shape_fix/intersection_tool.rs`): the
//! former `ShapeFixIntersectionToolGap` (consumed by `fix_api.rs` L958) and
//! `ShapeFixSplitToolGap` (consumed by `fix_adv.rs` and `fix_intersect.rs`
//! L226) carriers — deleted, Rule 4; every consumer now calls the real
//! `ShapeFix_IntersectionTool` / `ShapeFix_SplitTool`.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::geom::Plane;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, CurveRepresentation, Orientation, TShape};

use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_analysis::transfer_parameters_proj::ShapeAnalysisTransferParametersProj;
use crate::shhealing::shape_build::edge::{builder_range_on_face, ShapeBuildEdge};
use crate::shhealing::shape_build::reshape::ShapeBuildReShape;
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::edge::ShapeFixEdge;
use crate::shhealing::shape_fix::wire::fix_api::brep_face_surface_loc;
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::shhealing::shape_fix::wire::REAL_LAST;

/// OCCT `Geom_Curve::Period()` (Geom_Curve.cxx) — the period is the span of
/// the natural parameter domain of the periodic basis curve.
pub(crate) fn curve_period(c: &Curve3) -> f64 {
    let d = c.default_domain();
    d[1] - d[0]
}

/// OCCT `Geom2d_Curve::Period()` — the 2d form.
pub(crate) fn curve_period_2d(c: &Curve2d) -> f64 {
    let d = c.default_domain();
    d[1] - d[0]
}

// ---------------------------------------------------------------------------
// GAP re-hosts (module doc).
// ---------------------------------------------------------------------------

/// OCCT `GeomConvert_CompCurveToBSplineCurve` GAP re-host (module doc).
pub struct GeomConvertCompCurveToBSplineGap;

impl GeomConvertCompCurveToBSplineGap {
    /// OCCT GeomConvert_CompCurveToBSplineCurve(BSplineCurve)
    /// (GeomConvert_CompCurveToBSplineCurve.cxx L52-70).
    pub fn new(_curve: &Curve3) -> Self {
        GeomConvertCompCurveToBSplineGap
    }

    /// OCCT GeomConvert_CompCurveToBSplineCurve::Add(NewCurve, Tolerance,
    /// After, WithRatio) — GAP: keeps OCCT's `!Add` failure branch.
    pub fn add(
        &mut self,
        _new_curve: &Curve3,
        _tolerance: f64,
        _after: bool,
        _with_ratio: bool,
    ) -> bool {
        false
    }

    /// OCCT GeomConvert_CompCurveToBSplineCurve::BSplineCurve() — GAP: the
    /// no-result path of the failed concatenation.
    pub fn bspline_curve(&self) -> Option<Curve3> {
        None
    }
}

/// OCCT `GeomAPI::To3d(2dCurve, plane)` (GeomAPI_XX.cxx) — GAP: the
/// TKGeomBase planar conversion is kernel-geom scope; keeps the null-result
/// path consumed by the `IsNull` guards.
pub fn geom_api_to3d(_the_2d_curve: &Curve2d, _the_plane: &Plane) -> Option<Curve3> {
    None
}

/// OCCT `GeomAPI::To2d(3dCurve, plane)` — GAP: the TKGeomBase planar
/// conversion is kernel-geom scope; keeps the null-result path consumed by
/// the `IsNull` guards.
pub fn geom_api_to2d(_the_3d_curve: &Curve3, _the_plane: &Plane) -> Option<Curve2d> {
    None
}

/// OCCT `IntRes2d_IntersectionPoint::ParamOnFirst` (module doc — GAP note).
pub fn param_on_first(p: &DVec2) -> f64 {
    p.x
}

/// OCCT `IntRes2d_IntersectionPoint::ParamOnSecond` (module doc — GAP note).
pub fn param_on_second(p: &DVec2) -> f64 {
    p.y
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L1652-1657 — UpdateEdgeUVPoints (static).
// ---------------------------------------------------------------------------

/// OCCT static UpdateEdgeUVPoints (cxx L1652-1657).
pub(crate) fn update_edge_uv_points(brep: &mut BRep, e: &Shape, f: &Shape) {
    // OCCT L1655: BRep_Tool::Range(E, F, first, last).
    let (_c2d, first, last) = pcurve_range_on_face(brep, e, f);
    // OCCT L1656: BRep_Builder().Range(E, F, first, last).
    builder_range_on_face(brep, e, f, first, last);
}

/// OCCT BRep_Tool::Range(edge, face, first, last) — the pcurve range on the
/// face (the face-matched representation row).
pub(crate) fn pcurve_range_on_face(brep: &BRep, edge: &Shape, face: &Shape) -> (Option<Curve2d>, f64, f64) {
    let fptr = face.ptr_id();
    let reps = match brep.tshapes[edge.index].as_ref() {
        TShape::Edge(ed) => &ed.representations,
        _ => return (None, 0.0, 0.0),
    };
    for r in reps {
        let (key, pc, range) = match r {
            CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve, range),
            CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                range,
                ..
            } => (*face, pcurve1, range),
            _ => continue,
        };
        if key.0 == fptr {
            return (Some(pc.clone()), range[0], range[1]);
        }
    }
    (None, 0.0, 0.0)
}

/// OCCT BRep_Builder::UpdateEdge(edge, curve, tol) — sets the 3d curve of
/// the edge (the TEdgeData 3d-curve row) and raises the tolerance.
pub(crate) fn builder_update_edge_curve3d_value(brep: &mut BRep, edge: &Shape, curve: Curve3, tol: f64) {
    let ed = brep.edge_mut_inplace(edge.clone());
    ed.curve = Some(curve);
    ed.tolerance = ed.tolerance.max(tol);
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L2215-2252 — TryNewPCurve (static).
// ---------------------------------------------------------------------------

/// OCCT static TryNewPCurve (cxx L2215-2252): creates an edge on the basis
/// of E with the new pcurve and calls FixSameParameter; returns the
/// resulting tolerance and the modified pcurve.
pub(crate) fn try_new_pcurve(
    brep: &mut BRep,
    e: &Shape,
    face: &Shape,
    c2d: &mut Option<Curve2d>,
    first: &mut f64,
    last: &mut f64,
    tol: &mut f64,
) -> bool {
    // OCCT L2222-2223: crv = BRep_Tool::Curve(E, f, l).
    let (crv, mut f, mut l) = {
        match brep.tshapes[e.index].as_ref() {
            TShape::Edge(ed) => (ed.curve.clone(), ed.range[0], ed.range[1]),
            _ => (None, 0.0, 0.0),
        }
    };
    let crv = match crv {
        Some(c) => c,
        None => return false, // OCCT L2224-2227: return false.
    };

    // OCCT L2230: BRepBuilderAPI_MakeEdge mkedge(crv, f, l) — the rcad edge
    // constructor builds the end vertices at the curve ends (the
    // BRepBuilderAPI_MakeEdge(crv, f, l) contract).
    let mut builder = BRepBuilder::new();
    let p_f = rcad_kernel::geom::CurveEval::point_at(&crv, f);
    let p_l = rcad_kernel::geom::CurveEval::point_at(&crv, l);
    let v1 = builder.add_vertex(brep, p_f, CONFUSION);
    let v2 = builder.add_vertex(brep, p_l, CONFUSION);
    let edge = builder.add_edge(brep, Some(crv), v1, v2, [f, l]);

    // OCCT L2232-2233 (skl 17.07.2001): SBE.SetRange3d(mkedge, f, l).
    let sbe = ShapeBuildEdge;
    sbe.set_range3d(brep, &edge, f, l);

    // OCCT L2235-2238: mkedge.IsDone() — the rcad constructor has no failure
    // mode (annotated, bridge #5).
    let _edge2 = edge.clone();

    // OCCT L2242-2245.
    let c2d_val = c2d.clone().unwrap();
    builder.update_edge_pcurve(brep, edge.clone(), c2d_val.clone(), face.clone(), 0.0);
    builder_range_on_face(brep, &edge, face, *first, *last);
    builder.set_edge_same_range(brep, edge.clone(), false);
    // no call to BRepLib:  B.SameParameter ( edge, false );

    // OCCT L2247-2248.
    let mut sfe = ShapeFixEdge::new();
    sfe.fix_same_parameter_face(brep, &edge, face, 0.0);
    // OCCT L2249-2251.
    let (pc, fi, la) = brep_tool_curve_on_surface(brep, &edge, face);
    if let Some(pc) = pc {
        *c2d = Some(pc);
    }
    *first = fi;
    *last = la;
    *tol = brep_tool_tolerance(&edge);
    f += 0.0;
    l += 0.0;
    true
}

/// OCCT BRep_Tool::CurveOnSurface(edge, face) — the (pcurve, first, last)
/// of the face representation.
pub(crate) fn brep_tool_curve_on_surface(
    brep: &BRep,
    edge: &Shape,
    face: &Shape,
) -> (Option<Curve2d>, f64, f64) {
    pcurve_range_on_face(brep, edge, face)
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L2262-2279 — howMuchPCurves (static).
// ---------------------------------------------------------------------------

/// OCCT static howMuchPCurves (cxx L2262-2279): counts the pcurve
/// representations of the edge.
pub(crate) fn how_much_pcurves(e: &Shape) -> i32 {
    let mut count = 0;
    if let TShape::Edge(ed) = e.data.as_ref() {
        for cr in &ed.representations {
            match cr {
                CurveRepresentation::CurveOnSurface { .. }
                | CurveRepresentation::CurveOnClosedSurface { .. } => count += 1,
                _ => {}
            }
        }
    }
    count
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L2283-2523 — RemoveLoop (6-arg static).
// ---------------------------------------------------------------------------

/// OCCT static RemoveLoop(E, face, IP, tolfact, prec, RemoveLoop3d)
/// (cxx L2283-2523): cuts out the self-intersection loop on the pcurve of
/// the edge.
#[allow(clippy::too_many_arguments)]
pub(crate) fn remove_loop(
    brep: &mut BRep,
    e: &mut Shape,
    face: &Shape,
    ip: &DVec2,
    tolfact: f64,
    prec: f64,
    remove_loop3d: bool,
) -> bool {
    // OCCT L2291-2294.
    if brep.is_edge_closed_on_face(e, face) {
        return false;
    }

    // OCCT L2296-2297: crv = BRep_Tool::Curve(E, f, l).
    let (crv, f, l) = match brep.tshapes[e.index].as_ref() {
        TShape::Edge(ed) => (ed.curve.clone(), ed.range[0], ed.range[1]),
        _ => return false,
    };
    let crv = match crv {
        Some(c) => c,
        None => return false,
    };

    // OCCT L2299-2310.
    let mut t1 = param_on_first(ip);
    let mut t2 = param_on_second(ip);
    if t1 > t2 {
        std::mem::swap(&mut t1, &mut t2);
    }

    let sae = ShapeAnalysisEdge::new();
    let mut a = 0.0f64;
    let mut b = 0.0f64;
    let mut c2d: Option<Curve2d> = None;
    if !sae.pcurve_face(brep, e, face, &mut c2d, &mut a, &mut b, false) {
        return false;
    }

    // OCCT L2319-2322: GAC.Resolution(prec).
    let dt = tolfact * resolution_of_curve(&crv, prec);
    t1 -= dt; // 1e-3;//::Precision::PConfusion();
    t2 += dt; // 1e-3;//::Precision::PConfusion();

    // OCCT L2327-2330.
    if t1 <= a || t2 >= b {
        // should not be so, but to be sure ..
        return false;
    }

    // OCCT L2339-2348: PTV OCC884.
    let (a_plane_surf, _) = brep_face_surface_loc(brep, face);
    let pln = plane_of_surface(&a_plane_surf);
    let mut pcurve3d: Option<Curve3> = Some(crv.clone());
    if a_plane_surf.is_none() {
        // OCCT L2347: pcurve3d = GeomAPI::To3d(c2d, Pln).
        pcurve3d = c2d.as_ref().and_then(|c| geom_api_to3d(c, &pln));
    }

    // OCCT L2352: the first segment.
    let c2d_val = match &c2d {
        Some(c) => c.clone(),
        None => return false,
    };
    let _ = &c2d_val;
    let pcurve3d_val = match &pcurve3d {
        Some(c) => c.clone(),
        None => return false,
    };
    let trim = rcad_kernel::geom::TrimmedCurve3::new(pcurve3d_val.clone(), a, t1);
    let mut connect = GeomConvertCompCurveToBSplineGap::new(&Curve3::Trimmed(trim));

    // OCCT L2356-2369: the null-length segment patch instead of the loop.
    let trim = rcad_kernel::geom::TrimmedCurve3::new(pcurve3d_val.clone(), t2, b);
    let patch = rcad_kernel::geom::Curve3::Trimmed(trim);

    // OCCT L2370-2373.
    if !connect.add(&patch, PCONFUSION, true, false) {
        return false;
    }

    // OCCT L2376-2379: the last segment.
    if !connect.add(&patch, PCONFUSION, true, false) {
        return false;
    }
    // OCCT L2382: PTV OCC884 — keep the created 3d curve.
    let a_new_3d_crv = connect.bspline_curve();

    // OCCT L2384-2388.
    let bs = a_new_3d_crv.as_ref().and_then(|c| geom_api_to2d(c, &pln));
    let bs = match bs {
        Some(c) => c,
        None => return false,
    };

    // OCCT L2390-2436.
    let mut builder = BRepBuilder::new();
    if !remove_loop3d {
        // old variant (not remove loop 3d)
        let mut newtol = 0.0f64;
        // OCC901
        let nb_c2d = how_much_pcurves(e);
        if nb_c2d <= 1 && a_plane_surf.is_some() {
            builder.update_edge_pcurve(brep, e.clone(), bs.clone(), face.clone(), 0.0);
        } else {
            let mut bs_var: Option<Curve2d> = Some(bs.clone());
            if !try_new_pcurve(brep, e, face, &mut bs_var, &mut a, &mut b, &mut newtol) {
                return false;
            }
        }

        let tol = brep_tool_tolerance(e);
        if newtol > prec.max(tol) {
            return false;
        }
        //: s2  bs = BRep_Tool::CurveOnSurface ( edge, face, a, b );
        if (a - f).abs() > PCONFUSION || (b - l).abs() > PCONFUSION {
            // smth strange, cancel
            return false;
        }
        // PTV OCC884
        if a_plane_surf.is_some() {
            if let Some(c3d) = &a_new_3d_crv {
                builder_update_edge_curve3d_value(brep, e, c3d.clone(), newtol.max(tol));
            }
            // OCC901
            let mut bs_var: Option<Curve2d> = Some(bs.clone());
            if !try_new_pcurve(brep, e, face, &mut bs_var, &mut a, &mut b, &mut newtol) {
                return false;
            }
        }
        builder.update_edge_pcurve(brep, e.clone(), bs, face.clone(), newtol);
        let fv = sae.first_vertex(brep, e);
        let lv = sae.last_vertex(brep, e);
        builder.update_vertex_tolerance(brep, fv, newtol);
        builder.update_vertex_tolerance(brep, lv, newtol);
        return true;
    }

    // OCCT L2438-2444: :q1 — the Adaptor3d_CurveOnSurface values re-hosted
    // through the landed ShapeAnalysisSurface evaluation.
    let p1 = c2d
        .as_ref()
        .map(|c| {
            let uv = rcad_kernel::geom::Curve2dEval::point_at(c, t1);
            surface_value_at(brep, face, uv)
        })
        .unwrap_or(DVec3::ZERO);
    let p2 = c2d
        .as_ref()
        .map(|c| {
            let uv = rcad_kernel::geom::Curve2dEval::point_at(c, t2);
            surface_value_at(brep, face, uv)
        })
        .unwrap_or(DVec3::ZERO);
    let pcur_pnt = DVec3::new(
        (p1.x + p2.x) / 2.0,
        (p1.y + p2.y) / 2.0,
        (p1.z + p2.z) / 2.0,
    );

    let mut sftp = ShapeAnalysisTransferParametersProj::new_edge_face(brep, e, face);
    let seq2d = [t1, t2, (t1 + t2) / 2.0];
    let seq3d = sftp.perform_params(brep, &seq2d, false);

    // OCCT L2457-2462: correcting Seq3d already project.
    let dist1 = pcur_pnt.distance(value_at_param(&crv, seq3d[0]));
    let dist2 = pcur_pnt.distance(value_at_param(&crv, seq3d[1]));
    let dist3 = pcur_pnt.distance(value_at_param(&crv, seq3d[2]));
    let loop_removed3d = dist3 <= dist1.max(dist2);

    let bs1: Option<Curve3>;

    if !loop_removed3d {
        // OCCT L2466-2508: create the new 3d curve.
        let mut ftrim = seq3d[0];
        let mut ltrim = seq3d[1];
        ftrim -= dt;
        ltrim += dt;
        let trim1 = rcad_kernel::geom::TrimmedCurve3::new(crv.clone(), f, ftrim);
        let mut connect1 = GeomConvertCompCurveToBSplineGap::new(&Curve3::Trimmed(trim1));
        let trim1 = rcad_kernel::geom::TrimmedCurve3::new(crv.clone(), ltrim, l);
        let patch1 = rcad_kernel::geom::Curve3::Trimmed(trim1);

        // OCCT L2491-2494.
        if !connect1.add(&patch1, PCONFUSION, true, false) {
            return false;
        }
        // OCCT L2496-2500: the last segment.
        if !connect1.add(&patch1, PCONFUSION, true, false) {
            return false;
        }
        // OCCT L2502.
        bs1 = connect1.bspline_curve();
        if bs1.is_none() {
            return false;
        }
    } else {
        bs1 = None;
    }
    // OCCT L2509: double oldtol = BRep_Tool::Tolerance ( E );

    // OCCT L2511-2517.
    if !loop_removed3d {
        if let Some(c) = bs1 {
            builder_update_edge_curve3d_value(brep, e, c, 0.0);
        }
    }
    builder.update_edge_pcurve(brep, e.clone(), bs, face.clone(), 0.0);
    builder_range_on_face(brep, e, face, f, l);
    builder.set_edge_same_range(brep, e.clone(), false);

    // OCCT L2519-2520.
    let mut sfe = ShapeFixEdge::new();
    sfe.fix_same_parameter(brep, e, 0.0);

    true
}

/// OCCT GeomAdaptor_Curve::Resolution(prec) re-host — the parameter-space
/// resolution from the curve tolerance (the landed
/// `BSplineCurve3::bsplclib_resolution` covers the bspline arm; the generic
/// curve falls back to the tolerance itself).
fn resolution_of_curve(c: &Curve3, prec: f64) -> f64 {
    match c {
        Curve3::BSpline(bs) => bs.bsplclib_resolution(prec),
        _ => prec,
    }
}

/// OCCT `gp_Pln` of the face surface (the plane arm of RemoveLoop L2339).
fn plane_of_surface(s: &Option<Surface3>) -> Plane {
    match s {
        Some(Surface3::Plane(p)) => *p,
        // OCCT L2333-2335: the direct construction on the Z plane.
        None => Plane::new(DVec3::ZERO, DVec3::Z),
        _ => Plane::new(DVec3::ZERO, DVec3::Z),
    }
}

/// OCCT Adaptor3d_CurveOnSurface::Value — the surface 3d value of the 2d
/// point.
fn surface_value_at(brep: &BRep, face: &Shape, uv: DVec2) -> DVec3 {
    let (s, l) = brep_face_surface_loc(brep, face);
    match s {
        Some(s) => {
            let world = if l != 0 {
                let trsf = brep.get_location(l);
                rcad_kernel::geom::transform_surface(&s, &trsf)
            } else {
                s
            };
            surface_value(&world, uv)
        }
        None => DVec3::ZERO,
    }
}

/// The surface evaluation (the ShapeAnalysis_Surface::Value re-host).
fn surface_value(s: &Surface3, uv: DVec2) -> DVec3 {
    ShapeAnalysisSurface::new(s.clone()).value(uv)
}

/// OCCT Geom_Curve::Value.
fn value_at_param(c: &Curve3, u: f64) -> DVec3 {
    rcad_kernel::geom::CurveEval::point_at(c, u)
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L2527-2704 — RemoveLoop (E1,E2 static).
// ---------------------------------------------------------------------------

/// OCCT static RemoveLoop(E, face, IP2d, E1, E2) (cxx L2527-2704): tries to
/// insert a vertex at the self-intersection point and split the edge.
pub(crate) fn remove_loop_split(
    brep: &mut BRep,
    e: &Shape,
    face: &Shape,
    ip2d: &DVec2,
    e1: &mut Option<Shape>,
    e2: &mut Option<Shape>,
) -> bool {
    // OCCT L2537-2540.
    if brep.is_edge_closed_on_face(e, face) {
        return false;
    }

    // OCCT L2542-2543.
    let (crv, f, l) = match brep.tshapes[e.index].as_ref() {
        TShape::Edge(ed) => (ed.curve.clone(), ed.range[0], ed.range[1]),
        _ => return false,
    };
    let crv = match crv {
        Some(c) => c,
        None => return false,
    };

    // OCCT L2545-2553.
    let mut t1 = param_on_first(ip2d);
    let mut t2 = param_on_second(ip2d);
    if t1 > t2 {
        std::mem::swap(&mut t1, &mut t2);
    }

    let sae = ShapeAnalysisEdge::new();

    // OCCT L2557-2561: define Vfirst, Vlast, Vmid.
    let vfirst = sae.first_vertex(brep, e);
    let vlast = sae.last_vertex(brep, e);

    // OCCT L2563-2569.
    let mut a = 0.0f64;
    let mut b = 0.0f64;
    let mut c2d: Option<Curve2d> = None;
    if !sae.pcurve_face(brep, e, face, &mut c2d, &mut a, &mut b, false) {
        return false;
    }
    let c2d = match c2d {
        Some(c) => c,
        None => return false,
    };

    // OCCT L2571-2584: the first and second segments for the 2d curve.
    let mut trim1: Option<Curve2d> = None;
    if (t1 - a) > PCONFUSION {
        trim1 = Some(Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 { curve: Box::new(c2d.clone()), t_min: a, t_max: t1 }));
    }
    // OCCT L2578: the second segment for the 2d curve.
    let trim2 = Curve2d::Trimmed(rcad_kernel::geom::TrimmedCurve2 { curve: Box::new(c2d.clone()), t_min: t2, t_max: b });

    // OCCT L2580-2584.
    if trim1.is_none() {
        // OCCT: `if (trim1.IsNull() || trim2.IsNull())` — trim2 is never
        // null in the rcad value model (bridge #5 note); the trim1 arm is
        // the live check.
    }

    // OCCT L2586-2594: the Adaptor3d_CurveOnSurface values (the landed
    // ShapeAnalysisSurface evaluation).
    let p1 = {
        let uv = rcad_kernel::geom::Curve2dEval::point_at(&c2d, t1);
        surface_value_at(brep, face, uv)
    };
    let p2 = {
        let uv = rcad_kernel::geom::Curve2dEval::point_at(&c2d, t2);
        surface_value_at(brep, face, uv)
    };
    let pcur_pnt = DVec3::new(
        (p1.x + p2.x) / 2.0,
        (p1.y + p2.y) / 2.0,
        (p1.z + p2.z) / 2.0,
    );

    // OCCT L2596-2602: transfer the parameters from the pcurve to the 3d
    // curve.
    let mut sftp = ShapeAnalysisTransferParametersProj::new_edge_face(brep, e, face);
    let seq2d = [t1, t2, (t1 + t2) / 2.0];
    let seq3d = sftp.perform_params(brep, &seq2d, false);

    // OCCT L2605-2619.
    let dist1 = pcur_pnt.distance(value_at_param(&crv, seq3d[0]));
    let dist2 = pcur_pnt.distance(value_at_param(&crv, seq3d[1]));
    let dist3 = pcur_pnt.distance(value_at_param(&crv, seq3d[2]));
    let (ftrim, ltrim);
    if dist3 > dist1.max(dist2) {
        // is loop in 3d
        ftrim = seq3d[0];
        ltrim = seq3d[1];
    } else {
        // not loop in 3d
        ftrim = seq3d[2];
        ltrim = seq3d[2];
    }

    // OCCT L2624-2637: the trims for the 3d curve.
    let mut trim3: Option<Curve3> = None;
    if trim1.is_some() {
        trim3 = Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
            crv.clone(),
            f,
            ftrim,
        )));
    }
    // OCCT L2631: the second segment for the 3d curve.
    let trim4 = Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(crv.clone(), ltrim, l));

    // OCCT L2633-2637.
    if trim4_instantiable(&trim4).is_none() {
        return false;
    }

    // OCCT L2639-2642: the point for the middle vertex.
    let pnt1 = value_at_param(&crv, ftrim);
    let pnt2 = value_at_param(&crv, ltrim);
    let pmid = DVec3::new(
        (pnt1.x + pnt2.x) / 2.0,
        (pnt1.y + pnt2.y) / 2.0,
        (pnt1.z + pnt2.z) / 2.0,
    );

    let mut b = BRepBuilder::new();

    // OCCT L2645-2651: create the new copies for E1 and E2.
    let mut e1v: Option<Shape> = None;
    if trim1.is_some() {
        let ec = brep.empty_copied(e);
        e1v = Some(ec);
    }
    let e2v = brep.empty_copied(e);

    // OCCT L2653-2661: initialize Vmid.
    let vmid = if trim1.is_none() {
        b.add_vertex(brep, pnt2, 0.0)
    } else {
        b.add_vertex(brep, pmid, 0.0)
    };

    let sbe = ShapeBuildEdge;

    // OCCT L2664-2681: replace the vertices for the new edges E1 and E2.
    if e.orientation == Orientation::Forward {
        if e1v.is_some() {
            let e1s = e1v.as_ref().unwrap();
            let replaced = sbe.copy_replace_vertices(brep, e1s, &vfirst, &vmid);
            e1v = Some(replaced);
        }
        let t = sbe.copy_replace_vertices(brep, &e2v, &vmid, &vlast);
        *e2 = Some(t);
    } else {
        if e1v.is_some() {
            let e1s = e1v.as_ref().unwrap();
            let replaced = sbe.copy_replace_vertices(brep, e1s, &vmid, &vlast);
            e1v = Some(replaced);
        }
        let t = sbe.copy_replace_vertices(brep, &e2v, &vfirst, &vmid);
        *e2 = Some(t);
    }

    // OCCT L2683-2701: update the edges by the 2d and 3d curves.
    let mut my_sfe = ShapeFixEdge::new();
    if let Some(e1s) = &mut e1v {
        if let Some(t1c) = &trim1 {
            b.update_edge_pcurve(brep, e1s.clone(), t1c.clone(), face.clone(), 0.0);
        }
        if let Some(t3c) = &trim3 {
            builder_update_edge_curve3d_value(brep, e1s, t3c.clone(), 0.0);
        }
        b.set_edge_range(brep, e1s.clone(), f, ftrim);
        b.set_edge_same_range(brep, e1s.clone(), false);
        //    B.SameParameter(E1,false);
        my_sfe.fix_same_parameter(brep, e1s, 0.0);
        my_sfe.fix_vertex_tolerance(brep, e1s);
    }
    if e2.is_some() {
        let e2s = e2.as_mut().unwrap();
        b.update_edge_pcurve(brep, e2s.clone(), trim2.clone(), face.clone(), 0.0);
        builder_update_edge_curve3d_value(brep, e2s, trim4.clone(), 0.0);
        b.set_edge_range(brep, e2s.clone(), ltrim, l);
        b.set_edge_same_range(brep, e2s.clone(), false);
        //  B.SameParameter(E2, false);
        my_sfe.fix_same_parameter(brep, e2s, 0.0);
        my_sfe.fix_vertex_tolerance(brep, e2s);
    }
    *e1 = e1v;

    true
}

/// The always-instantiable trim4 guard (the OCCT `trim4.IsNull()` check is
/// dead in the rcad value model — bridge #5 note).
fn trim4_instantiable(t: &Curve3) -> Option<&Curve3> {
    Some(t)
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L2917-2962 — ComputeLocalDeviation (static).
// ---------------------------------------------------------------------------

/// OCCT static ComputeLocalDeviation (cxx L2917-2962).
pub(crate) fn compute_local_deviation(
    brep: &BRep,
    edge: &Shape,
    pint: DVec3,
    pnt: DVec3,
    mut f: f64,
    mut l: f64,
    face: &Shape,
) -> f64 {
    let sae = ShapeAnalysisEdge::new();
    let mut c3d: Option<Curve3> = None;
    let mut a = 0.0f64;
    let mut b = 0.0f64;
    if !sae.curve3d(brep, edge, &mut c3d, &mut a, &mut b, false) {
        return REAL_LAST;
    }
    let c3d = match c3d {
        Some(c) => c,
        None => return REAL_LAST,
    };

    // OCCT L2932: gp_Lin line(pint, gp_Vec(pint, pnt)).
    let line_dir = (pnt - pint).normalize_or_zero();

    // OCCT L2934-2947.
    let mut crv: Option<Curve2d> = None;
    let mut fp = 0.0f64;
    let mut lp = 0.0f64;
    if sae.pcurve_face(brep, edge, face, &mut crv, &mut fp, &mut lp, false) {
        if let Some(c) = &crv {
            if let Curve2d::Trimmed(tc) = c {
                if matches!(tc.curve.as_ref(), Curve2d::Line(_)) {
                    f = a + (f - fp) * (b - a) / (lp - fp);
                    l = a + (l - fp) * (b - a) / (lp - fp);
                }
            }
        }
    }

    // OCCT L2949-2961.
    const NSEG: i32 = 10;
    let step = (l - f) / NSEG as f64;
    let mut dev = 0.0f64;
    for i in 1..NSEG {
        let p = rcad_kernel::geom::CurveEval::point_at(&c3d, f + i as f64 * step);
        // OCCT gp_Lin::Distance(p) — the point-to-line distance.
        let d = (p - pint).cross(line_dir).length();
        if dev < d {
            dev = d;
        }
    }
    dev
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L3532-3613 — TryBendingPCurve (static).
// ---------------------------------------------------------------------------

/// OCCT static TryBendingPCurve (cxx L3532-3613): :s2 abv 21 Apr 99 — tries
/// to bend the pcurve towards the given 2d point.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_bending_pcurve(
    brep: &mut BRep,
    e: &Shape,
    face: &Shape,
    p2d: DVec2,
    end: bool,
    c2d: &mut Option<Curve2d>,
    first: &mut f64,
    last: &mut f64,
    tol: &mut f64,
) -> bool {
    let sae = ShapeAnalysisEdge::new();
    if !sae.pcurve_face(brep, e, face, c2d, first, last, false) {
        return false;
    }

    // OCCT L3547-3610: the try/catch arm — the failure arm returns false
    // (bridge #5); the bspline conversion arm below keeps the OCCT flow.
    let bs: Option<Curve2d> = match c2d.as_ref() {
        Some(Curve2d::BSpline(_)) => c2d.clone(),
        // OCCT L3556-3559: Geom2dConvert::CurveToBSplineCurve(trim) — the
        // TKGeomBase conversion GAP keeps the catch path.
        _ => None,
    };
    let mut bs = match bs {
        Some(Curve2d::BSpline(bs)) => bs,
        _ => return false,
    };

    let par = if end { *last } else { *first };
    let nb_knots = bspline2_nb_knots(&bs);
    let mult_first = bspline2_multiplicity(&bs, 1);
    let mult_last = bspline2_multiplicity(&bs, nb_knots);
    if (bs.knots[0] - par).abs() < PCONFUSION && mult_first > bs.degree as i32 {
        // OCCT L3570: bs->SetPole(1, p2d).
        bs.control_points[0] = p2d;
    } else if (bs.knots[bs.knots.len() - 1] - par).abs() < PCONFUSION && mult_last > bs.degree as i32 {
        // OCCT L3575: bs->SetPole(bs->NbPoles(), p2d).
        let n = bs.control_points.len();
        bs.control_points[n - 1] = p2d;
    } else {
        // OCCT L3579: bs->Segment(first, last) — the TKGeomBase segment
        // reduction GAP keeps the catch path.
        return false;
    }
    *c2d = Some(Curve2d::BSpline(bs));

    // OCCT L3597-3600.
    if !try_new_pcurve(brep, e, face, c2d, first, last, tol) {
        return false;
    }

    true
}

/// The number of distinct knots (OCCT Geom2d_BSplineCurve::NbKnots — the
/// knots vector is stored expanded).
fn bspline2_nb_knots(bs: &rcad_kernel::geom::BSplineCurve2) -> usize {
    let mut n = 0usize;
    let mut prev = f64::NAN;
    for k in &bs.knots {
        if *k != prev {
            n += 1;
            prev = *k;
        }
    }
    n
}

/// The multiplicity of the num-th knot (OCCT
/// Geom2d_BSplineCurve::Multiplicity — the run length in the expanded
/// vector).
fn bspline2_multiplicity(bs: &rcad_kernel::geom::BSplineCurve2, num: usize) -> i32 {
    let mut n = 0usize;
    let mut prev = f64::NAN;
    let mut run = 0i32;
    for k in &bs.knots {
        if *k != prev {
            n += 1;
            run = 1;
            prev = *k;
        } else {
            run += 1;
        }
        if n == num {
            break;
        }
    }
    run
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire.cxx L4108-4174 — CopyReversePcurves (static).
// ---------------------------------------------------------------------------

/// OCCT static CopyReversePcurves (cxx L4108-4174).
pub(crate) fn copy_reverse_pcurves(brep: &mut BRep, toedge: &Shape, fromedge: &Shape, reverse: bool) {
    // OCCT L4112-4113: the edge locations.
    let from_loc = fromedge.location;
    let to_loc = toedge.location;

    // OCCT L4114-4117: iterate the pcurve representations of fromedge.
    let from_reps: Vec<CurveRepresentation> = match fromedge.data.as_ref() {
        TShape::Edge(ed) => ed.representations.clone(),
        _ => return,
    };
    for from_gc in &from_reps {
        let (surface_key, pcurve, range, is_closed_rep, pcurve2) = match from_gc {
            CurveRepresentation::CurveOnSurface {
                face,
                pcurve,
                range,
            } => (*face, pcurve.clone(), *range, false, None),
            CurveRepresentation::CurveOnClosedSurface {
                face,
                pcurve1,
                pcurve2,
                range,
                ..
            } => (*face, pcurve1.clone(), *range, true, Some(pcurve2.clone())),
            _ => continue,
        };

        // OCCT L4129-4144: search the matching representation in toedge.
        let mut found = false;
        if let TShape::Edge(ed) = toedge.data.as_ref() {
            for to_gc in &ed.representations {
                let (tkey, _) = match to_gc {
                    CurveRepresentation::CurveOnSurface { face, .. } => (*face, 0),
                    CurveRepresentation::CurveOnClosedSurface { face, .. } => (*face, 0),
                    _ => continue,
                };
                if tkey == surface_key {
                    found = true;
                    break;
                }
            }
        }
        if found {
            continue;
        }

        // OCCT L4145-4170: append the copied representation.
        let mut fp = range[0];
        let mut lp = range[1];
        let mut pcurve_new = pcurve.clone();
        if reverse {
            // OCCT L4154-4159: pcurve->ReversedParameter(U) = -U; Reverse().
            fp = -fp;
            lp = -lp;
            pcurve_new = rcad_kernel::geom::reverse_curve2d(&pcurve_new);
            let tmp = fp;
            fp = lp;
            lp = tmp;
        }
        // OCCT L4162: bug OCC209 — newLoc = (fromLoc * L).Predivided(toLoc).
        let l_loc = surface_key.1;
        let new_loc = predivided_location(brep, from_loc, l_loc, to_loc);
        let new_key = (surface_key.0, new_loc);
        let rep = if is_closed_rep {
            CurveRepresentation::CurveOnClosedSurface {
                face: new_key,
                pcurve1: pcurve_new.clone(),
                pcurve2: pcurve2.clone().unwrap(),
                range: [fp, lp],
            }
        } else {
            CurveRepresentation::CurveOnSurface {
                face: new_key,
                pcurve: pcurve_new,
                range: [fp, lp],
            }
        };
        brep.edge_mut_inplace(toedge.clone())
            .representations
            .push(rep);
    }
}

/// OCCT `(fromLoc * L).Predivided(toLoc)` — the location composition over
/// the BRep location table.
fn predivided_location(brep: &mut BRep, from_loc: u32, l_loc: u32, to_loc: u32) -> u32 {
    let from_m = brep.get_location(from_loc);
    let l_m = brep.get_location(l_loc);
    let to_m = brep.get_location(to_loc);
    let new_m = from_m * l_m * to_m.inverse();
    brep.add_location(new_m)
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire_1.cxx L126-152 — AdjustOnPeriodic3d (static).
// ---------------------------------------------------------------------------

/// OCCT static AdjustOnPeriodic3d (Wire_1.cxx L126-152).
pub(crate) fn adjust_on_periodic3d(
    c: &Curve3,
    takefirst: bool,
    first: f64,
    last: f64,
    param: f64,
) -> f64 {
    // 15.11.2002 PTV OCC966
    let sa = crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
    if sa.is_periodic(c) {
        let t = curve_period(c);
        let mut shift = -(first / t).trunc() * t;
        if first < 0.0 {
            shift += t;
        }
        let sfirst = first + shift;
        let slast = last + shift;
        if takefirst && (param > slast) && (param > sfirst) {
            return param - t - shift;
        }
        if !takefirst && (param < slast) && (param < sfirst) {
            return param + t - shift;
        }
    }
    param
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire_1.cxx L865-891 — AdjustOnPeriodic2d (static).
// ---------------------------------------------------------------------------

/// OCCT static AdjustOnPeriodic2d (Wire_1.cxx L865-891).
pub(crate) fn adjust_on_periodic2d(
    pc: &Curve2d,
    takefirst: bool,
    first: f64,
    last: f64,
    param: f64,
) -> f64 {
    // 15.11.2002 PTV OCC966
    let sa = crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
    if sa.is_periodic_2d(pc) {
        let t = curve_period_2d(pc);
        let mut shift = -(first / t).trunc() * t;
        if first < 0.0 {
            shift += t;
        }
        let sfirst = first + shift;
        let slast = last + shift;
        if takefirst && (param > slast) && (param > sfirst) {
            return param - t - shift;
        }
        if !takefirst && (param < slast) && (param < sfirst) {
            return param + t - shift;
        }
    }
    param
}
