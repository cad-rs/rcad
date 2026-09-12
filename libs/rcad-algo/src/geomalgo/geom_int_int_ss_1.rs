// OCCT GeomInt_IntSS — the intersection engine (TKGeomAlgo/GeomInt), 1:1 Rust
// translation.  This file carries the part of the batch that lives in
// `GeomInt_IntSS_1.cxx` plus the two `IntPatch_Intersection` statics the
// GeomInt_IntSS pipeline calls.
//
// OCCT sources:
//   GeomInt_IntSS_1.cxx L55-102   file-static AdjustUPeriodic
//   GeomInt_IntSS_1.cxx L106-147  file-static GetQuadric / Parameters
//   GeomInt_IntSS_1.cxx L151-175  file-static ParametersOfNearestPointOnSurface
//   GeomInt_IntSS_1.cxx L179-241  file-static GetSegmentBoundary /
//                                 IntersectCurveAndBoundary
//   GeomInt_IntSS_1.cxx L247-271  file-static isDegenerated
//   GeomInt_IntSS_1.cxx L275-1094 GeomInt_IntSS::MakeCurve
//   GeomInt_IntSS_1.cxx L1098-1168 GeomInt_IntSS::TreatRLine
//   GeomInt_IntSS_1.cxx L1172-1304 GeomInt_IntSS::BuildPCurves (full bounds)
//   GeomInt_IntSS_1.cxx L1308-1325 GeomInt_IntSS::BuildPCurves (surface bounds)
//   GeomInt_IntSS_1.cxx L1333-1448 GeomInt_IntSS::TrimILineOnSurfBoundaries
//   GeomInt_IntSS_1.cxx L1452-1469 GeomInt_IntSS::MakeBSpline
//   GeomInt_IntSS_1.cxx L1473-1502 GeomInt_IntSS::MakeBSpline2d
//   IntPatch_Intersection.cxx L2286-2368 CheckSingularPoints
//   IntPatch_Intersection.cxx L2372-2408 DefineUVMaxStep
//   IntPatch_Intersection.cxx L2412-2477 splitCone + PrepareSurfaces
//
// Architecture differences (Rust vs C++):
//   - `Adaptor3d_Surface` -> [`GeomSurfaceAdapter`]; `Adaptor3d_TopolTool` ->
//     the landed `TopolTool`.
//   - `occ::handle<Geom_Curve>` / `handle<Geom2d_Curve>` -> `Option<Curve3>` /
//     `Option<Curve2d>` (OCCT appends NULL handles for the "not approximated"
//     arms).
//   - `Bnd_Box2d` -> `rcad_kernel::math::bnd::BndBox2d`.
//   - `GeomProjLib::Curve2d` -> `rcad_kernel::base::geom_proj_lib::curve2d`.
//   - `Geom2dAdaptor::MakeCurve(Adaptor2d_Curve2d)` is the identity here: the
//     rcad `IntPatchLine` carries the restriction arcs as `Curve2d` values
//     already (`arc_on_s1` / `arc_on_s2`), so the adaptor-to-curve conversion of
//     OCCT has no counterpart.

use glam::{DVec2, DVec3};

use rcad_kernel::geom::{
    BSplineCurve3, Curve2d, Curve2dEval, Curve3, CurveEval, Plane, Surface3, TrimmedCurve2,
    TrimmedCurve3,
};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::precision::{
    is_infinite_value, is_negative_infinite_value, is_positive_infinite_value, CONFUSION,
    INFINITE_VALUE, PCONFUSION,
};
use rcad_kernel::topods::State;

use super::geom_int_int_ss::{uv_rect, GeomIntIntSS};
// The MakeBSpline/BuildPCurves/TrimILineOnSurfBoundaries/TreatRLine group
// lives in the sibling file (2000-line split); imported back verbatim.
use super::geom_int_int_ss_1_curves::{
    bspline2_from_poles, bspline3_from_poles, build_p_curves_bounds, make_b_spline,
    make_b_spline_2d, treat_r_line, trim_i_line_on_surf_boundaries,
};
use crate::geomalgo::approx_int::{define_par_type, WLineAccess, WLineApprox};
use crate::geomalgo::geom2d_int::{Curve2dAdaptor, GInter};
use crate::geomalgo::geom_int_line_constructor::GeomIntLineConstructor;
use crate::geomalgo::int_patch::{GeomAbsSurfaceType, IntPatchIType, IntPatchLine};
use crate::geomalgo::int_surf::Quadric;
use crate::hlr::contap::geom_tool::GeomTool;
use crate::hlr::contap::surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter};
use crate::topalgo::adaptor3d::topol_tool::TopolTool;

/// OCCT `RealEpsilon()` (Standard_Real.hxx) = DBL_EPSILON.
const REAL_EPSILON: f64 = f64::EPSILON;

/// OCCT compares `occ::handle<Geom_Surface>` identity (`S1 == S2`,
/// `myHS1 != myHS2`, `theS1 != theS2`).  rcad carries the `Geom_Surface` by
/// value and `Surface3` has no `PartialEq`, so the structural identity of the
/// value (the derived `Debug` rendering, which covers every variant payload in
/// field order) is the counterpart.  Called O(1) times per `Perform`.
pub(crate) fn surface3_identity(a: &Surface3, b: &Surface3) -> bool {
    format!("{a:?}") == format!("{b:?}")
}

/// The adaptor identity (`aHS1 == aHS2` compares the OCCT adaptor handles):
/// the carried surface AND the restricted UV window must match, because a
/// `VTrim`/`UTrim` product is a distinct OCCT object.
pub(crate) fn adapter_identity(a: &GeomSurfaceAdapter, b: &GeomSurfaceAdapter) -> bool {
    surface3_identity(a.surface3(), b.surface3()) && a.uv_window() == b.uv_window()
}

// ============================================================================
// OCCT IntPatch_Intersection.cxx L2372-2477 — DefineUVMaxStep /
// CheckSingularPoints / splitCone / PrepareSurfaces.
//
// OCCT homes these in IntPatch_Intersection; the rcad home is this batch's
// `geomalgo` file because their only caller is GeomInt_IntSS::InternalPerform
// (the FF pipeline passes its own corrected domains to the intersector).  They
// are translated verbatim from the OCCT statements, statics included.
// ============================================================================

/// OCCT `static void splitCone(theS, theD, theTol, theVecHS)`
/// (IntPatch_Intersection.cxx L2412-2443).
fn split_cone(
    the_s: &GeomSurfaceAdapter,
    the_d: &TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
    the_tol: f64,
    the_vec_hs: &mut Vec<GeomSurfaceAdapter>,
) {
    if the_s.get_type() != GeomAbsSurfaceType::Cone {
        panic!("Standard_NoSuchObject: IntPatch_Intersection : Surface is not Cone");
    }

    let a_cone = the_s.cone();

    let mut a_u0 = 0.0;
    let mut a_v0 = 0.0;
    crate::topalgo::adaptor3d::topol_tool::get_cone_apex_param(&a_cone, &mut a_u0, &mut a_v0);

    let a_state = the_d.classify(DVec2::new(a_u0, a_v0), the_tol, true);

    if a_state == State::In || a_state == State::On {
        // OCCT L2431-2434:
        //   theS->VTrim(theS->FirstVParameter(), aV0, Precision::PConfusion())
        //   theS->VTrim(aV0, theS->LastVParameter(), Precision::PConfusion())
        let a_hs_dn = the_s.clone_with_window([
            the_s.first_u_parameter(),
            the_s.last_u_parameter(),
            the_s.first_v_parameter(),
            a_v0,
        ]);
        let a_hs_up = the_s.clone_with_window([
            the_s.first_u_parameter(),
            the_s.last_u_parameter(),
            a_v0,
            the_s.last_v_parameter(),
        ]);

        the_vec_hs.push(a_hs_dn);
        the_vec_hs.push(a_hs_up);
    } else {
        the_vec_hs.push(the_s.clone());
    }
}

/// OCCT IntPatch_Intersection::PrepareSurfaces (IntPatch_Intersection.cxx
/// L2449-2477).
pub(crate) fn int_patch_intersection_prepare_surfaces(
    the_s1: &GeomSurfaceAdapter,
    the_d1: &TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
    the_s2: &GeomSurfaceAdapter,
    the_d2: &TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
    the_tol: f64,
    the_vec_hs1: &mut Vec<GeomSurfaceAdapter>,
    the_vec_hs2: &mut Vec<GeomSurfaceAdapter>,
) {
    // OCCT L2458-2466.
    if the_s1.get_type() == GeomAbsSurfaceType::Cone
        && (std::f64::consts::FRAC_PI_2 - the_s1.cone().half_angle_rad.abs()).abs() < the_tol
    {
        split_cone(the_s1, the_d1, the_tol, the_vec_hs1);
    } else {
        the_vec_hs1.push(the_s1.clone());
    }

    // OCCT L2468-2476.
    if the_s2.get_type() == GeomAbsSurfaceType::Cone
        && (std::f64::consts::FRAC_PI_2 - the_s2.cone().half_angle_rad.abs()).abs() < the_tol
    {
        split_cone(the_s2, the_d2, the_tol, the_vec_hs2);
    } else {
        the_vec_hs2.push(the_s2.clone());
    }
}

/// OCCT IntPatch_Intersection::CheckSingularPoints (IntPatch_Intersection.cxx
/// L2286-2368).
fn check_singular_points(
    the_s1: &GeomSurfaceAdapter,
    the_d1: &mut TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
    the_s2: &GeomSurfaceAdapter,
    the_dist: &mut f64,
) -> bool {
    *the_dist = INFINITE_VALUE;
    let mut is_singular = false;
    // OCCT L2293-2296: if (theS1 == theS2) return isSingular.
    if adapter_identity(the_s1, the_s2) {
        return is_singular;
    }
    // OCCT L2298-2299.
    const A_NB_BND_PNTS: i32 = 5;
    const A_TOL: f64 = CONFUSION;
    // OCCT L2301-2303: theD1->Init(); for (; theD1->More(); theD1->Next()).
    the_d1.init();
    let mut is_u = true;
    while the_d1.more() {
        let a_bnd = *the_d1.value();
        let pinf = a_bnd.first_parameter();
        let psup = a_bnd.last_parameter();
        if is_negative_infinite_value(pinf) || is_positive_infinite_value(psup) {
            the_d1.next();
            continue;
        }
        let dt = (psup - pinf) / (A_NB_BND_PNTS - 1) as f64;
        // OCCT L2314: aBnd->D1((pinf + psup) / 2., aP1, aDir).
        let (a_p1_mid, a_dir) = a_bnd.d1((pinf + psup) / 2.);
        let _ = a_p1_mid;
        is_u = a_dir.x.abs() > a_dir.y.abs();
        let mut a_d1_norm_max = 0.0f64;
        let mut a_pmid = DVec3::ZERO;
        let mut a_nb = 0i32;
        // OCCT L2321: for (t = pinf; t <= psup; t += dt).
        let mut t = pinf;
        let mut a_p1;
        loop {
            if !(t <= psup) {
                break;
            }
            a_p1 = a_bnd.value(t);
            // OCCT L2324: theS1->D1(aP1.X(), aP1.Y(), aPP1, aDU, aDV).
            let (a_pp1, a_du, a_dv) = the_s1.d1(a_p1.x, a_p1.y);
            if is_u {
                a_d1_norm_max = a_d1_norm_max.max(a_du.length());
            } else {
                a_d1_norm_max = a_d1_norm_max.max(a_dv.length());
            }

            a_pmid += a_pp1;
            a_nb += 1;

            if a_d1_norm_max > A_TOL {
                break;
            }
            t += dt;
        }

        if a_d1_norm_max <= A_TOL {
            // Singular point aPP1.
            let mut a_pp1 = DVec3::ZERO;
            a_pp1 += a_pmid / a_nb as f64;
            let a_tol_u = PCONFUSION;
            let a_tol_v = PCONFUSION;
            // OCCT L2349: Extrema_ExtPS aProj(aPP1, *theS2.get(), aTolU,
            // aTolV, Extrema_ExtFlag_MIN).
            let a_proj = rcad_kernel::base::extrema::ExtPS::new(
                a_pp1,
                the_s2.surface3(),
                a_tol_u,
                a_tol_v,
            );

            if a_proj.is_done() {
                let a_nb_ext = a_proj.nb_ext();
                let mut i = 1usize;
                while i <= a_nb_ext {
                    *the_dist = the_dist.min(a_proj.square_distance(i));
                    i += 1;
                }
            }
        }
        the_d1.next();
    }
    if !is_infinite_value(*the_dist) {
        *the_dist = the_dist.sqrt();
        is_singular = true;
    }

    is_singular
}

/// OCCT IntPatch_Intersection::DefineUVMaxStep (IntPatch_Intersection.cxx
/// L2372-2408).
pub(crate) fn define_uv_max_step(
    the_s1: &GeomSurfaceAdapter,
    the_d1: &mut TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
    the_s2: &GeomSurfaceAdapter,
    the_d2: &mut TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
) -> f64 {
    let mut an_uv_max_step = 0.001;
    let mut a_dist_to_sing1 = INFINITE_VALUE;
    let mut a_dist_to_sing2 = INFINITE_VALUE;
    let a_tol_min = CONFUSION;
    let a_tol_max = 1.0e-5;
    // OCCT L2381: if (theS1 != theS2) { ... }.
    if !adapter_identity(the_s1, the_s2) {
        let mut is_sing1 = check_singular_points(the_s1, the_d1, the_s2, &mut a_dist_to_sing1);
        if is_sing1 {
            if a_dist_to_sing1 > a_tol_min && a_dist_to_sing1 < a_tol_max {
                an_uv_max_step = 0.0001;
            } else {
                is_sing1 = false;
            }
        }
        if !is_sing1 {
            let is_sing2 = check_singular_points(the_s2, the_d2, the_s1, &mut a_dist_to_sing2);
            if is_sing2 && a_dist_to_sing2 > a_tol_min && a_dist_to_sing2 < a_tol_max {
                an_uv_max_step = 0.0001;
            }
        }
    }
    an_uv_max_step
}

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L55-102: file-static AdjustUPeriodic
// ============================================================================

/// OCCT `static void AdjustUPeriodic(aS, aC2D)` (cxx L55-102).
fn adjust_u_periodic(a_s: &GeomSurfaceAdapter, a_c2d: &mut Option<Curve2d>) {
    let Some(c2d) = a_c2d.as_ref() else {
        return;
    };
    if !a_s.is_u_periodic() {
        return;
    }
    // constexpr double aEps = Precision::PConfusion(); // 1.e-9
    let a_eps = PCONFUSION;
    // const double aEpsilon = Epsilon(10.); // 1.77e-15
    let a_epsilon = epsilon_of(10.);
    //
    let (umin, umax) = (a_s.first_u_parameter(), a_s.last_u_parameter());
    let a_period = a_s.u_period();

    let a_t1 = curve2d_first_parameter(c2d);
    let a_t2 = curve2d_last_parameter(c2d);
    let a_tx = a_t1 + 0.467 * (a_t2 - a_t1);
    let a_px = c2d.point_at(a_tx);
    //
    let mut a_ux = a_px.x;
    if a_ux.abs() < a_epsilon {
        a_ux = 0.;
    }
    if (a_ux - a_period).abs() < a_epsilon {
        a_ux = a_period;
    }
    //
    let mut d_u = 0.;
    while a_ux < (umin - a_eps) {
        a_ux += a_period;
        d_u += a_period;
    }
    while a_ux > (umax + a_eps) {
        a_ux -= a_period;
        d_u -= a_period;
    }
    //
    if d_u != 0. {
        // gp_Vec2d aV2D(dU, 0.); aC2D->Translate(aV2D);
        let translated = rcad_kernel::geom::translate_curve2d(c2d, DVec2::new(d_u, 0.));
        *a_c2d = Some(translated);
    }
}

/// OCCT `Epsilon(theValue)` (Standard_Real.hxx L242-246).
pub(super) fn epsilon_of(x: f64) -> f64 {
    if x >= 0.0 {
        next_after(x, f64::INFINITY) - x
    } else {
        x - next_after(x, f64::NEG_INFINITY)
    }
}

/// OCCT `std::nextafter` (bit-level construction; see the kernel twin).
fn next_after(x: f64, to: f64) -> f64 {
    if x.is_nan() || to.is_nan() || x == to {
        return x;
    }
    if x == 0.0 {
        return if to > 0.0 { f64::from_bits(1) } else { -f64::from_bits(1) };
    }
    if (to > x) == (x > 0.0) {
        f64::from_bits(x.to_bits() + 1)
    } else {
        f64::from_bits(x.to_bits() - 1)
    }
}

/// OCCT Geom2d_Curve::FirstParameter.
pub(super) fn curve2d_first_parameter(c: &Curve2d) -> f64 {
    c.default_domain()[0]
}

/// OCCT Geom2d_Curve::LastParameter.
pub(super) fn curve2d_last_parameter(c: &Curve2d) -> f64 {
    c.default_domain()[1]
}

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L106-147: file-static GetQuadric / Parameters
// ============================================================================

/// OCCT `static void GetQuadric(HS1, quad1)` (cxx L106-128).
fn get_quadric(hs1: &GeomSurfaceAdapter) -> Quadric {
    match hs1.get_type() {
        GeomAbsSurfaceType::Plane => Quadric::from_plane(&hs1.plane()),
        GeomAbsSurfaceType::Cylinder => Quadric::from_cylinder(&hs1.cylinder()),
        GeomAbsSurfaceType::Cone => Quadric::from_cone(&hs1.cone()),
        GeomAbsSurfaceType::Sphere => Quadric::from_sphere(&hs1.sphere()),
        GeomAbsSurfaceType::Torus => Quadric::from_torus(&hs1.torus()),
        // OCCT L126: throw Standard_ConstructionError("GeomInt_IntSS::MakeCurve").
        _ => panic!("Standard_ConstructionError: GeomInt_IntSS::MakeCurve"),
    }
}

/// OCCT `static void Parameters(HS1, HS2, Ptref, U1, V1, U2, V2)` (cxx
/// L132-147).
fn parameters(
    hs1: &GeomSurfaceAdapter,
    hs2: &GeomSurfaceAdapter,
    ptref: DVec3,
) -> (f64, f64, f64, f64) {
    let quad1 = get_quadric(hs1);
    let quad2 = get_quadric(hs2);
    let (u1, v1) = quad1.parameters(ptref);
    let (u2, v2) = quad2.parameters(ptref);
    (u1, v1, u2, v2)
}

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L151-175:
// file-static ParametersOfNearestPointOnSurface
// ============================================================================

/// OCCT `static bool ParametersOfNearestPointOnSurface(theExtr, theU, theV)`
/// (cxx L151-175).
pub(super) fn parameters_of_nearest_point_on_surface(
    the_extr: &rcad_kernel::base::extrema::ExtPS,
    the_u: &mut f64,
    the_v: &mut f64,
) -> bool {
    if !the_extr.is_done() || the_extr.nb_ext() == 0 {
        return false;
    }

    let mut an_index = 1usize;
    let mut a_min_sq_dist = the_extr.square_distance(an_index);
    let mut i = 2usize;
    while i <= the_extr.nb_ext() {
        let a_sqd = the_extr.square_distance(i);
        if a_sqd < a_min_sq_dist {
            a_min_sq_dist = a_sqd;
            an_index = i;
        }
        i += 1;
    }

    // theExtr.Point(anIndex).Parameter(theU, theV)
    let p = the_extr.point(an_index);
    *the_u = p.u;
    *the_v = p.v;

    true
}

/// OCCT Extrema_ExtPS::Point(Index) is 1-based; the rcad kernel accessor is
/// 1-based as well (`square_distance(n)` / `point(n)`).
const _EXT_REM_1BASED: () = ();

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L179-241: file-static GetSegmentBoundary /
// IntersectCurveAndBoundary
// ============================================================================

/// OCCT `static void GetSegmentBoundary(theSegm, theCurve,
/// theArrayOfParameters)` (cxx L179-199).
fn get_segment_boundary(
    the_segm: &crate::geomalgo::int_res2d::IntersectionSegment,
    the_curve: &Curve2d,
    the_array_of_parameters: &mut Vec<f64>,
) {
    let mut a_u1 = curve2d_first_parameter(the_curve);
    let mut a_u2 = curve2d_last_parameter(the_curve);

    if the_segm.has_first_point() {
        let an_ipf = the_segm.first_point();
        a_u1 = an_ipf.param_on_first();
    }

    if the_segm.has_last_point() {
        let an_ipl = the_segm.last_point();
        a_u2 = an_ipl.param_on_first();
    }

    the_array_of_parameters.push(a_u1);
    the_array_of_parameters.push(a_u2);
}

/// OCCT `static void IntersectCurveAndBoundary(theC2d, theArrBounds,
/// theNumberOfCurves, theTol, theArrayOfParameters)` (cxx L203-241).
pub(super) fn intersect_curve_and_boundary(
    the_c2d: &Option<Curve2d>,
    the_arr_bounds: &[Option<Curve2d>],
    the_number_of_curves: usize,
    the_tol: f64,
    the_array_of_parameters: &mut Vec<f64>,
) {
    let Some(c2d) = the_c2d.as_ref() else {
        return;
    };

    // Geom2dAdaptor_Curve anAC1(theC2d);
    // rcad: `impl Curve2dAdaptor for Curve2d` makes the Geom2d_Curve itself the
    // adaptor (the OCCT Geom2dAdaptor_Curve wrapper has no separate state).
    let an_ac1: &dyn Curve2dAdaptor = c2d;
    let mut a_cur_id = 0usize;
    while a_cur_id < the_number_of_curves {
        let Some(bound) = the_arr_bounds[a_cur_id].as_ref() else {
            a_cur_id += 1;
            continue;
        };

        // Geom2dAdaptor_Curve anAC2(theArrBounds[aCurID]);
        let an_ac2: &dyn Curve2dAdaptor = bound;
        // Geom2dInt_GInter anIntCC2d(anAC1, anAC2, theTol, theTol);
        let mut an_int_cc2d = GInter::new();
        an_int_cc2d.perform_cc(an_ac1, an_ac2, the_tol, the_tol);

        // OCCT IntRes2d_Intersection::IsEmpty() == no points and no segments;
        // the rcad IntCurveCurveGen exposes the two counts but not the base
        // predicate.
        if !an_int_cc2d.is_done()
            || (an_int_cc2d.nb_points() == 0 && an_int_cc2d.nb_segments() == 0)
        {
            a_cur_id += 1;
            continue;
        }

        let mut a_pnt_id = 1usize;
        while a_pnt_id <= an_int_cc2d.nb_points() {
            let a_param = an_int_cc2d.point(a_pnt_id).param_on_first();
            the_array_of_parameters.push(a_param);
            a_pnt_id += 1;
        }

        let mut a_segm_id = 1usize;
        while a_segm_id <= an_int_cc2d.nb_segments() {
            get_segment_boundary(
                an_int_cc2d.segment(a_segm_id),
                c2d,
                the_array_of_parameters,
            );
            a_segm_id += 1;
        }
        a_cur_id += 1;
    }
}

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L247-271: file-static isDegenerated
// ============================================================================

/// OCCT `static bool isDegenerated(theGAHS, theAHC2d, theFirstPar, theLastPar)`
/// (cxx L247-271).
pub(super) fn is_degenerated(
    the_gahs: &GeomSurfaceAdapter,
    the_ahc2d: &Curve2d,
    the_first_par: f64,
    the_last_par: f64,
) -> bool {
    // constexpr double aSqTol = Precision::Confusion() * Precision::Confusion();
    let a_sq_tol = CONFUSION * CONFUSION;

    let a_p2d = the_ahc2d.point_at(the_first_par);
    let a_p1 = the_gahs.value(a_p2d.x, a_p2d.y);
    let a_p2d = the_ahc2d.point_at(the_last_par);
    let mut a_p2 = the_gahs.value(a_p2d.x, a_p2d.y);

    if a_p1.distance_squared(a_p2) > a_sq_tol {
        return false;
    }

    let a_p2d = the_ahc2d.point_at(0.5 * (the_first_par + the_last_par));
    a_p2 = the_gahs.value(a_p2d.x, a_p2d.y);

    a_p1.distance_squared(a_p2) <= a_sq_tol
}

/// The effective [tf, tl] of the trimmed arc (OCCT's `tf`/`tl` after the
/// min/max clamping against the arc's own range).

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L275-1094: MakeCurve
// ============================================================================

impl GeomIntIntSS {
    /// OCCT GeomInt_IntSS::MakeCurve(Ind, D1, D2, Tol, Approx, Approx1,
    /// Approx2) (cxx L275-1094).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn make_curve(
        &mut self,
        ind: usize,
        dom1: &TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
        dom2: &TopolTool<'_, GeomSurfaceAdapter, GeomTool>,
        tol: f64,
        approx: bool,
        approx_s1: bool,
        approx_s2: bool,
    ) {
        // OCCT L284-285.
        let my_approx;
        let my_approx1;
        let my_approx2;
        let mut my_tol_approx;
        //
        // OCCT L290-294.
        let mut tolpc = tol;
        my_approx = approx;
        my_approx1 = approx_s1;
        my_approx2 = approx_s2;
        my_tol_approx = 0.0000001;
        //
        // OCCT L296-297.
        let a_s1 = self.my_hs1().clone();
        let a_s2 = self.my_hs2().clone();
        //
        // OCCT L299-300.
        let l0 = self.my_intersector().line(ind - 1).clone();
        let mut typl = l0.line_type;
        //
        // OCCT L302-310.
        if typl == IntPatchIType::Walking {
            if !l0.is_wline() {
                return;
            }
            typl = IntPatchIType::Walking;
        }
        // OCCT L299: L = myIntersector.Line(Index); OCCT L304 reassigns L to the
        // down_cast WLine — in rcad both are the same `IntPatchLine`.
        let l = l0;
        //
        // Line Constructor
        // OCCT L313: myLConstruct.Perform(L);
        let a_l = l.clone();
        self.my_l_construct_mut().perform(&a_l);
        if !self.my_l_construct().is_done() || self.my_l_construct().nb_parts() == 0 {
            return;
        }
        // OCCT L287: occ::handle<Geom2d_BSplineCurve> H1 (a NULL handle).
        let h1: Option<Curve2d> = None;

        // Running copies of the OCCT result state (written back at the end;
        // the append order of sline / slineS1 / slineS2 is preserved because
        // every arm appends to all three in lockstep).
        let mut sline: Vec<Option<Curve3>> = Vec::new();
        let mut sline_s1: Vec<Option<Curve2d>> = Vec::new();
        let mut sline_s2: Vec<Option<Curve2d>> = Vec::new();
        let mut my_tol_reached_2d = *self.my_tol_reached_2d_mut();

        match typl {
            // ########################################
            //  Line, Parabola, Hyperbola
            // ########################################
            IntPatchIType::Line | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
                // OCCT L332-343.
                let newc = match typl {
                    IntPatchIType::Line => match &l.curve {
                        Curve3::Line(ln) => Curve3::Line(ln.clone()),
                        _ => panic!("IntPatch_GLine down-cast: Line"),
                    },
                    IntPatchIType::Parabola => match &l.curve {
                        Curve3::Parabola(p) => Curve3::Parabola(p.clone()),
                        _ => panic!("IntPatch_GLine down-cast: Parabola"),
                    },
                    IntPatchIType::Hyperbola => match &l.curve {
                        Curve3::Hyperbola(h) => Curve3::Hyperbola(h.clone()),
                        _ => panic!("IntPatch_GLine down-cast: Hyperbola"),
                    },
                    _ => unreachable!(),
                };
                //
                let a_nb_parts = self.my_l_construct().nb_parts();
                let mut i = 1usize;
                while i <= a_nb_parts {
                    let (fprm, lprm) = self.my_l_construct().part(i);

                    // OCCT L350.
                    if !is_negative_infinite_value(fprm) && !is_positive_infinite_value(lprm) {
                        // OCCT L352-353.
                        let a_ct3d = Curve3::Trimmed(TrimmedCurve3::new(
                            newc.clone(),
                            fprm,
                            lprm,
                        ));
                        sline.push(Some(a_ct3d.clone()));
                        //
                        // OCCT L355-368.
                        if my_approx1 {
                            let mut c2d: Option<Curve2d> = None;
                            build_p_curves_bounds(
                                fprm,
                                lprm,
                                &mut tolpc,
                                a_s1.surface3(),
                                &newc,
                                &mut c2d,
                            );
                            if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                my_tol_reached_2d = tolpc;
                            }
                            sline_s1.push(c2d.map(|c| {
                                Curve2d::Trimmed(TrimmedCurve2 {
                                    curve: Box::new(c),
                                    t_min: fprm,
                                    t_max: lprm,
                                })
                            }));
                        } else {
                            sline_s1.push(h1.clone());
                        }
                        //
                        // OCCT L370-384.
                        if my_approx2 {
                            let mut c2d: Option<Curve2d> = None;
                            build_p_curves_bounds(
                                fprm,
                                lprm,
                                &mut tolpc,
                                a_s2.surface3(),
                                &newc,
                                &mut c2d,
                            );
                            if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                my_tol_reached_2d = tolpc;
                            }
                            //
                            sline_s2.push(c2d.map(|c| {
                                Curve2d::Trimmed(TrimmedCurve2 {
                                    curve: Box::new(c),
                                    t_min: fprm,
                                    t_max: lprm,
                                })
                            }));
                        } else {
                            sline_s2.push(h1.clone());
                        }
                    } else {
                        // OCCT L387-432.
                        let typ_s1 = a_s1.get_type();
                        let typ_s2 = a_s2.get_type();
                        if typ_s1 == GeomAbsSurfaceType::SurfaceOfExtrusion
                            || typ_s1 == GeomAbsSurfaceType::OffsetSurface
                            || typ_s1 == GeomAbsSurfaceType::SurfaceOfRevolution
                            || typ_s2 == GeomAbsSurfaceType::SurfaceOfExtrusion
                            || typ_s2 == GeomAbsSurfaceType::OffsetSurface
                            || typ_s2 == GeomAbsSurfaceType::SurfaceOfRevolution
                        {
                            sline.push(Some(newc.clone()));
                            sline_s1.push(h1.clone());
                            sline_s2.push(h1.clone());
                            i += 1;
                            continue;
                        }
                        // OCCT L400-405.
                        let d_t = 100.;
                        let b_fn_it = is_negative_infinite_value(fprm);
                        let b_lp_it = is_positive_infinite_value(lprm);

                        let mut a_test_prm = 0.;
                        // OCCT L409-416.
                        if b_fn_it && !b_lp_it {
                            a_test_prm = lprm - d_t;
                        } else if !b_fn_it && b_lp_it {
                            a_test_prm = fprm + d_t;
                        }
                        //
                        // gp_Pnt ptref(newc->Value(aTestPrm));
                        let ptref = newc.point_at(a_test_prm);
                        //
                        // OCCT L420-432.
                        let tol_x = CONFUSION;
                        let (u1, v1, u2, v2) = parameters(&a_s1, &a_s2, ptref);
                        let mut ok = dom1.classify(DVec2::new(u1, v1), tol_x, true) != State::Out;
                        if ok {
                            ok = dom2.classify(DVec2::new(u2, v2), tol_x, true) != State::Out;
                        }
                        if ok {
                            sline.push(Some(newc.clone()));
                            sline_s1.push(h1.clone());
                            sline_s2.push(h1.clone());
                        }
                    }
                    i += 1;
                }
            }

            // ########################################
            //  Circle and Ellipse
            // ########################################
            IntPatchIType::Circle | IntPatchIType::Ellipse => {
                // OCCT L444-451.
                let newc = if typl == IntPatchIType::Circle {
                    match &l.curve {
                        Curve3::Circle(c) => Curve3::Circle(c.clone()),
                        _ => panic!("IntPatch_GLine down-cast: Circle"),
                    }
                } else {
                    match &l.curve {
                        Curve3::Ellipse(e) => Curve3::Ellipse(e.clone()),
                        _ => panic!("IntPatch_GLine down-cast: Ellipse"),
                    }
                };
                //
                // OCCT L453-456.
                let a_real_epsilon = REAL_EPSILON;
                let a_period = std::f64::consts::PI + std::f64::consts::PI;
                //
                let a_nb_parts = self.my_l_construct().nb_parts();
                //
                let mut i = 1usize;
                while i <= a_nb_parts {
                    let (fprm, lprm) = self.my_l_construct().part(i);
                    //
                    // OCCT L464.
                    if fprm.abs() > a_real_epsilon
                        || (lprm - a_period).abs() > a_real_epsilon
                    {
                        // ==============================================
                        let a_tc3d =
                            Curve3::Trimmed(TrimmedCurve3::new(newc.clone(), fprm, lprm));
                        //
                        sline.push(Some(a_tc3d.clone()));
                        //
                        // OCCT L471-472: fprm = aTC3D->FirstParameter();
                        // lprm = aTC3D->LastParameter().  Geom_TrimmedCurve
                        // stores the bounds given to Geom_Curve::Trimmed, so
                        // the reassignment is the identity in rcad's
                        // TrimmedCurve3 (same [first, last] pair).
                        ////
                        if my_approx1 {
                            let mut c2d: Option<Curve2d> = None;
                            build_p_curves_bounds(
                                fprm,
                                lprm,
                                &mut tolpc,
                                a_s1.surface3(),
                                &newc,
                                &mut c2d,
                            );
                            if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                my_tol_reached_2d = tolpc;
                            }
                            sline_s1.push(c2d);
                        } else {
                            ////
                            sline_s1.push(h1.clone());
                        }
                        //
                        if my_approx2 {
                            let mut c2d: Option<Curve2d> = None;
                            build_p_curves_bounds(
                                fprm,
                                lprm,
                                &mut tolpc,
                                a_s2.surface3(),
                                &newc,
                                &mut c2d,
                            );
                            if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                my_tol_reached_2d = tolpc;
                            }
                            sline_s2.push(c2d);
                        } else {
                            sline_s2.push(h1.clone());
                        }
                        // ==============================================
                    } else {
                        // on regarde si on garde
                        //
                        // OCCT L509-550.
                        if a_nb_parts == 1 {
                            if fprm.abs() < a_real_epsilon
                                && (lprm - 2. * std::f64::consts::PI).abs() < a_real_epsilon
                            {
                                let a_tc3d =
                                    Curve3::Trimmed(TrimmedCurve3::new(newc.clone(), fprm, lprm));
                                //
                                sline.push(Some(a_tc3d.clone()));

                                if my_approx1 {
                                    let mut c2d: Option<Curve2d> = None;
                                    build_p_curves_bounds(
                                        fprm,
                                        lprm,
                                        &mut tolpc,
                                        a_s1.surface3(),
                                        &newc,
                                        &mut c2d,
                                    );
                                    if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                        my_tol_reached_2d = tolpc;
                                    }
                                    sline_s1.push(c2d);
                                } else {
                                    ////
                                    sline_s1.push(h1.clone());
                                }

                                if my_approx2 {
                                    let mut c2d: Option<Curve2d> = None;
                                    build_p_curves_bounds(
                                        fprm,
                                        lprm,
                                        &mut tolpc,
                                        a_s2.surface3(),
                                        &newc,
                                        &mut c2d,
                                    );
                                    if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                        my_tol_reached_2d = tolpc;
                                    }
                                    sline_s2.push(c2d);
                                } else {
                                    sline_s2.push(h1.clone());
                                }
                                break;
                            }
                        }
                        //
                        // OCCT L552-602.
                        let a_two_pi_div17 = 2. * std::f64::consts::PI / 17.;
                        //
                        let mut j = 0i32;
                        while j <= 17 {
                            let ptref = newc.point_at(j as f64 * a_two_pi_div17);
                            let tol_x = CONFUSION;

                            let (u1, v1, u2, v2) = parameters(&a_s1, &a_s2, ptref);
                            let mut ok =
                                dom1.classify(DVec2::new(u1, v1), tol_x, true) != State::Out;
                            if ok {
                                ok = dom2.classify(DVec2::new(u2, v2), tol_x, true) != State::Out;
                            }
                            if ok {
                                sline.push(Some(newc.clone()));
                                // ==============================================
                                if my_approx1 {
                                    let mut c2d: Option<Curve2d> = None;
                                    build_p_curves_bounds(
                                        fprm,
                                        lprm,
                                        &mut tolpc,
                                        a_s1.surface3(),
                                        &newc,
                                        &mut c2d,
                                    );
                                    if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                        my_tol_reached_2d = tolpc;
                                    }
                                    sline_s1.push(c2d);
                                } else {
                                    sline_s1.push(h1.clone());
                                }

                                if my_approx2 {
                                    let mut c2d: Option<Curve2d> = None;
                                    build_p_curves_bounds(
                                        fprm,
                                        lprm,
                                        &mut tolpc,
                                        a_s2.surface3(),
                                        &newc,
                                        &mut c2d,
                                    );
                                    if tolpc > my_tol_reached_2d || my_tol_reached_2d == 0. {
                                        my_tol_reached_2d = tolpc;
                                    }
                                    sline_s2.push(c2d);
                                } else {
                                    sline_s2.push(h1.clone());
                                }
                                break;
                            }
                            j += 1;
                        }
                    }
                    i += 1;
                }
            }

            // ########################################
            //  Analytic
            // ########################################
            IntPatchIType::Analytic => {
                // This case was processed earlier (in IntPatch_Intersection).
            }

            // ########################################
            //  Walking
            // ########################################
            IntPatchIType::Walking => {
                let mut wl = l.clone();

                //
                let mut ifprm: i32;
                let mut ilprm: i32;
                //
                if !my_approx {
                    // OCCT L630-654.
                    let a_nb_parts = self.my_l_construct().nb_parts();
                    let mut i = 1usize;
                    while i <= a_nb_parts {
                        let (fprm, lprm) = self.my_l_construct().part(i);
                        ifprm = fprm as i32;
                        ilprm = lprm as i32;
                        //
                        let mut a_h1: Option<Curve2d> = None;
                        let mut a_h2: Option<Curve2d> = None;

                        if my_approx1 {
                            a_h1 = make_b_spline_2d(&wl, ifprm, ilprm, true);
                        }
                        if my_approx2 {
                            a_h2 = make_b_spline_2d(&wl, ifprm, ilprm, false);
                        }
                        //
                        let a_b_sp = make_b_spline(&wl, ifprm, ilprm);
                        //
                        sline.push(a_b_sp);
                        sline_s1.push(a_h1);
                        sline_s2.push(a_h2);
                        i += 1;
                    }
                } else {
                    // OCCT L657-1002.
                    let mut nbiter: usize;
                    let a_nb_seq_of_l: usize;
                    let b_is_decomposited: bool;
                    let tol2d: f64;
                    let a_tol_ss = 2.0e-7;
                    //
                    tol2d = my_tol_approx;
                    // theapp3d.SetParameters(myTolApprox, tol2d, 4, 8, 0, 30,
                    // myHS1 != myHS2);
                    let mut theapp3d = WLineApprox::new();
                    theapp3d.set_parameters(
                        my_tol_approx,
                        tol2d,
                        4,
                        8,
                        0,
                        30,
                        !adapter_identity(&a_s1, &a_s2),
                        crate::geomalgo::approx_int::ApproxParamType::ChordLength,
                    );
                    //
                    // OCCT L669-670: bIsDecomposited =
                    // GeomInt_LineTool::DecompositionOfWLine(WL, myHS1, myHS2,
                    // aTolSS, myLConstruct, aSeqOfL);
                    let mut a_seq_of_l: Vec<IntPatchLine> = Vec::new();
                    b_is_decomposited = decomposition_of_w_line(
                        &wl,
                        &a_s1,
                        &a_s2,
                        a_tol_ss,
                        self.my_l_construct(),
                        &mut a_seq_of_l,
                    );
                    //
                    // OCCT L672-675.
                    let a_nb_parts = self.my_l_construct().nb_parts();
                    a_nb_seq_of_l = a_seq_of_l.len();
                    //
                    nbiter = if b_is_decomposited {
                        a_nb_seq_of_l
                    } else {
                        a_nb_parts
                    };
                    //
                    let mut i = 1usize;
                    while i <= nbiter {
                        if b_is_decomposited {
                            // OCCT L681-683.
                            wl = a_seq_of_l[i - 1].clone();
                            ifprm = 1;
                            ilprm = wl.nb_points() as i32;
                        } else {
                            let (fprm, lprm) = self.my_l_construct().part(i);
                            ifprm = fprm as i32;
                            ilprm = lprm as i32;
                        }

                        // OCCT L692-708.
                        let mut an_approx = my_approx;
                        let mut an_approx1 = my_approx1;
                        let mut an_approx2 = my_approx2;
                        let typs1 = a_s1.get_type();
                        let typs2 = a_s2.get_type();

                        if typs1 == GeomAbsSurfaceType::Plane {
                            an_approx = false;
                            an_approx1 = true;
                        } else if typs2 == GeomAbsSurfaceType::Plane {
                            an_approx = false;
                            an_approx2 = true;
                        }

                        // OCCT L710-713.
                        let a_par_type = define_par_type(
                            &WLineAccess {
                                line: &wl,
                                indicemin: ifprm as usize,
                                indicemax: ilprm as usize,
                                nbp3d: 1,
                                nbp2d: 2,
                                approx_u1v1: an_approx1,
                                approx_u2v2: an_approx2,
                                p2d_on_first: true,
                                xo: 0.0,
                                yo: 0.0,
                                zo: 0.0,
                                u1o: 0.0,
                                v1o: 0.0,
                                u2o: 0.0,
                                v2o: 0.0,
                                s1: a_s1.surface3(),
                                s2: a_s2.surface3(),
                                uv1: uv_rect(&a_s1),
                                uv2: uv_rect(&a_s2),
                            },
                            ifprm as usize,
                            ilprm as usize,
                            an_approx,
                            an_approx1,
                            an_approx2,
                        );

                        theapp3d.set_parameters(
                            my_tol_approx,
                            tol2d,
                            4,
                            8,
                            0,
                            30,
                            !adapter_identity(&a_s1, &a_s2),
                            a_par_type,
                        );

                        // -- lbr :
                        // -- Si une des surfaces est un plan, on approxime en 2d
                        // -- sur cette surface et on remonte les points 2d en 3d.
                        //
                        // OCCT L719-741.
                        if typs1 == GeomAbsSurfaceType::Plane {
                            theapp3d.perform(
                                &WLineAccess {
                                    line: &wl,
                                    indicemin: ifprm as usize,
                                    indicemax: ilprm as usize,
                                    nbp3d: 1,
                                    nbp2d: 2,
                                    approx_u1v1: an_approx1,
                                    approx_u2v2: an_approx2,
                                    p2d_on_first: true,
                                    xo: 0.0,
                                    yo: 0.0,
                                    zo: 0.0,
                                    u1o: 0.0,
                                    v1o: 0.0,
                                    u2o: 0.0,
                                    v2o: 0.0,
                                    s1: a_s1.surface3(),
                                    s2: a_s2.surface3(),
                                    uv1: uv_rect(&a_s1),
                                    uv2: uv_rect(&a_s2),
                                },
                                false,
                                true,
                                my_approx2,
                                ifprm as usize,
                                ilprm as usize,
                            );
                        } else if typs2 == GeomAbsSurfaceType::Plane {
                            theapp3d.perform(
                                &WLineAccess {
                                    line: &wl,
                                    indicemin: ifprm as usize,
                                    indicemax: ilprm as usize,
                                    nbp3d: 1,
                                    nbp2d: 2,
                                    approx_u1v1: an_approx1,
                                    approx_u2v2: an_approx2,
                                    p2d_on_first: true,
                                    xo: 0.0,
                                    yo: 0.0,
                                    zo: 0.0,
                                    u1o: 0.0,
                                    v1o: 0.0,
                                    u2o: 0.0,
                                    v2o: 0.0,
                                    s1: a_s1.surface3(),
                                    s2: a_s2.surface3(),
                                    uv1: uv_rect(&a_s1),
                                    uv2: uv_rect(&a_s2),
                                },
                                false,
                                my_approx1,
                                true,
                                ifprm as usize,
                                ilprm as usize,
                            );
                        } else {
                            //
                            // OCCT L730-738.
                            if !adapter_identity(&a_s1, &a_s2) {
                                if (typs1 == GeomAbsSurfaceType::BezierSurface
                                    || typs1 == GeomAbsSurfaceType::BSplineSurface)
                                    && (typs2 == GeomAbsSurfaceType::BezierSurface
                                        || typs2 == GeomAbsSurfaceType::BSplineSurface)
                                {
                                    theapp3d.set_parameters(
                                        my_tol_approx,
                                        tol2d,
                                        4,
                                        8,
                                        0,
                                        30,
                                        true,
                                        a_par_type,
                                    );
                                }
                            }
                            //
                            theapp3d.perform(
                                &WLineAccess {
                                    line: &wl,
                                    indicemin: ifprm as usize,
                                    indicemax: ilprm as usize,
                                    nbp3d: 1,
                                    nbp2d: 2,
                                    approx_u1v1: an_approx1,
                                    approx_u2v2: an_approx2,
                                    p2d_on_first: true,
                                    xo: 0.0,
                                    yo: 0.0,
                                    zo: 0.0,
                                    u1o: 0.0,
                                    v1o: 0.0,
                                    u2o: 0.0,
                                    v2o: 0.0,
                                    s1: a_s1.surface3(),
                                    s2: a_s2.surface3(),
                                    uv1: uv_rect(&a_s1),
                                    uv2: uv_rect(&a_s2),
                                },
                                true,
                                my_approx1,
                                my_approx2,
                                ifprm as usize,
                                ilprm as usize,
                            );
                        }

                        if !theapp3d.is_done() {
                            // OCCT L743-761.
                            let mut a_h1: Option<Curve2d> = None;
                            let mut a_h2: Option<Curve2d> = None;
                            //
                            let a_b_sp = make_b_spline(&wl, ifprm, ilprm);
                            if my_approx1 {
                                a_h1 = make_b_spline_2d(&wl, ifprm, ilprm, true);
                            }
                            if my_approx2 {
                                a_h2 = make_b_spline_2d(&wl, ifprm, ilprm, false);
                            }
                            //
                            sline.push(a_b_sp);
                            sline_s1.push(a_h1);
                            sline_s2.push(a_h2);
                        } else {
                            // OCCT L763-1001.
                            if my_approx1
                                || my_approx2
                                || (typs1 == GeomAbsSurfaceType::Plane
                                    || typs2 == GeomAbsSurfaceType::Plane)
                            {
                                if theapp3d.my_tol_reached2d > my_tol_reached_2d
                                    || my_tol_reached_2d == 0.
                                {
                                    my_tol_reached_2d = theapp3d.my_tol_reached2d;
                                }
                            }
                            if typs1 == GeomAbsSurfaceType::Plane
                                || typs2 == GeomAbsSurfaceType::Plane
                            {
                                *self.my_tol_reached_3d_mut() = my_tol_reached_2d;
                            } else if theapp3d.my_tol_reached3d
                                > *self.my_tol_reached_3d_mut()
                                || *self.my_tol_reached_3d_mut() == 0.
                            {
                                *self.my_tol_reached_3d_mut() = theapp3d.my_tol_reached3d;
                            }

                            let a_nb_multi_curves = theapp3d.nb_multi_curves();
                            //
                            let mut j = 1usize;
                            while j <= a_nb_multi_curves {
                                if typs1 == GeomAbsSurfaceType::Plane {
                                    // OCCT L786-853.
                                    let mbspc = theapp3d.value();
                                    let nbpoles = mbspc.nb_poles();

                                    let mut tpoles2d: Vec<DVec2> = Vec::new();
                                    let mut tpoles: Vec<DVec3> = Vec::new();

                                    mbspc.curve2d(1, &mut tpoles2d);
                                    let pln = a_s1.plane();
                                    //
                                    for p2d in tpoles2d.iter() {
                                        tpoles.push(el_s_lib_plane_value(p2d.x, p2d.y, &pln));
                                    }
                                    //
                                    let bs = bspline3_from_poles(
                                        tpoles,
                                        &mbspc.knots,
                                        &mbspc.mults,
                                        mbspc.degree,
                                    );
                                    // GeomLib_CheckBSplineCurve Check(BS,
                                    // myTolCheck, myTolAngCheck);
                                    // Check.FixTangent(true, true);
                                    let mut bs = bs;
                                    geom_lib_check_bspline_curve_fix_tangent(
                                        &mut bs,
                                        self.my_tol_check(),
                                        self.my_tol_ang_check(),
                                    );
                                    //
                                    sline.push(Some(bs));
                                    //
                                    if my_approx1 {
                                        let mut bs1 = Some(bspline2_from_poles(
                                            tpoles2d.clone(),
                                            &mbspc.knots,
                                            &mbspc.mults,
                                            mbspc.degree,
                                        ));
                                        geom_lib_check2d_bspline_curve_fix_tangent(
                                            &mut bs1,
                                            self.my_tol_check(),
                                            self.my_tol_ang_check(),
                                        );
                                        //
                                        adjust_u_periodic(&a_s1, &mut bs1);
                                        //
                                        sline_s1.push(bs1);
                                    } else {
                                        sline_s1.push(h1.clone());
                                    }

                                    if my_approx2 {
                                        mbspc.curve2d(2, &mut tpoles2d);

                                        let mut bs2 = Some(bspline2_from_poles(
                                            tpoles2d.clone(),
                                            &mbspc.knots,
                                            &mbspc.mults,
                                            mbspc.degree,
                                        ));
                                        geom_lib_check2d_bspline_curve_fix_tangent(
                                            &mut bs2,
                                            self.my_tol_check(),
                                            self.my_tol_ang_check(),
                                        );
                                        //
                                        adjust_u_periodic(&a_s2, &mut bs2);
                                        //
                                        sline_s2.push(bs2);
                                    } else {
                                        sline_s2.push(h1.clone());
                                    }
                                    let _ = nbpoles;
                                } else if typs2 == GeomAbsSurfaceType::Plane {
                                    // OCCT L855-922.
                                    let mbspc = theapp3d.value();
                                    let nbpoles = mbspc.nb_poles();

                                    let mut tpoles2d: Vec<DVec2> = Vec::new();
                                    let mut tpoles: Vec<DVec3> = Vec::new();
                                    mbspc.curve2d(if my_approx1 { 2 } else { 1 }, &mut tpoles2d);
                                    let pln = a_s2.plane();
                                    //
                                    for p2d in tpoles2d.iter() {
                                        tpoles.push(el_s_lib_plane_value(p2d.x, p2d.y, &pln));
                                    }
                                    //
                                    let bs = bspline3_from_poles(
                                        tpoles,
                                        &mbspc.knots,
                                        &mbspc.mults,
                                        mbspc.degree,
                                    );
                                    let mut bs = bs;
                                    geom_lib_check_bspline_curve_fix_tangent(
                                        &mut bs,
                                        self.my_tol_check(),
                                        self.my_tol_ang_check(),
                                    );
                                    //
                                    sline.push(Some(bs));
                                    //
                                    if my_approx2 {
                                        let mut bs1 = Some(bspline2_from_poles(
                                            tpoles2d.clone(),
                                            &mbspc.knots,
                                            &mbspc.mults,
                                            mbspc.degree,
                                        ));
                                        geom_lib_check2d_bspline_curve_fix_tangent(
                                            &mut bs1,
                                            self.my_tol_check(),
                                            self.my_tol_ang_check(),
                                        );
                                        //
                                        adjust_u_periodic(&a_s2, &mut bs1);
                                        //
                                        sline_s2.push(bs1);
                                    } else {
                                        sline_s2.push(h1.clone());
                                    }

                                    if my_approx1 {
                                        mbspc.curve2d(1, &mut tpoles2d);
                                        let mut bs2 = Some(bspline2_from_poles(
                                            tpoles2d.clone(),
                                            &mbspc.knots,
                                            &mbspc.mults,
                                            mbspc.degree,
                                        ));
                                        geom_lib_check2d_bspline_curve_fix_tangent(
                                            &mut bs2,
                                            self.my_tol_check(),
                                            self.my_tol_ang_check(),
                                        );
                                        //
                                        adjust_u_periodic(&a_s1, &mut bs2);
                                        //
                                        sline_s1.push(bs2);
                                    } else {
                                        sline_s1.push(h1.clone());
                                    }
                                    let _ = nbpoles;
                                } else {
                                    // typs1 != GeomAbs_Plane && typs2 != GeomAbs_Plane
                                    // OCCT L924-999.
                                    let mbspc = theapp3d.value();
                                    let nbpoles = mbspc.nb_poles();
                                    let mut tpoles: Vec<DVec3> = Vec::new();
                                    mbspc.curve(1, &mut tpoles);
                                    let bs = bspline3_from_poles(
                                        tpoles,
                                        &mbspc.knots,
                                        &mbspc.mults,
                                        mbspc.degree,
                                    );
                                    let mut bs = match bs {
                                        Curve3::BSpline(b) => b,
                                        _ => unreachable!(),
                                    };
                                    geom_lib_check_bspline_curve_fix_tangent(
                                        &mut Curve3::BSpline(bs.clone()),
                                        self.my_tol_check(),
                                        self.my_tol_ang_check(),
                                    );
                                    //
                                    // Check IsClosed() (OCCT L937-956).
                                    let dom = bs.default_domain();
                                    let a_dist = bs
                                        .point_at(dom[0])
                                        .length_squared()
                                        .max(bs.point_at(dom[1]).length_squared());
                                    let eps = epsilon_of(a_dist);
                                    if bs.point_at(dom[0]).distance_squared(bs.point_at(dom[1]))
                                        < 2. * eps
                                    {
                                        // Avoid creating B-splines containing two
                                        // coincident poles only.
                                        if mbspc.degree == 1 && nbpoles == 2 {
                                            j += 1;
                                            continue;
                                        }

                                        if !curve3_bspline_is_closed(&bs) && !bs.is_periodic {
                                            // force Closed()
                                            let a_pm = (bs.control_points[0]
                                                + bs.control_points[nbpoles - 1])
                                                / 2.;
                                            bs.control_points[0] = a_pm;
                                            bs.control_points[nbpoles - 1] = a_pm;
                                        }
                                    }
                                    sline.push(Some(Curve3::BSpline(bs)));

                                    if my_approx1 {
                                        let mut tpoles2d: Vec<DVec2> = Vec::new();
                                        mbspc.curve2d(2, &mut tpoles2d);
                                        let mut bs1 = Some(bspline2_from_poles(
                                            tpoles2d,
                                            &mbspc.knots,
                                            &mbspc.mults,
                                            mbspc.degree,
                                        ));
                                        geom_lib_check2d_bspline_curve_fix_tangent(
                                            &mut bs1,
                                            self.my_tol_check(),
                                            self.my_tol_ang_check(),
                                        );
                                        //
                                        adjust_u_periodic(&a_s1, &mut bs1);
                                        //
                                        sline_s1.push(bs1);
                                    } else {
                                        sline_s1.push(h1.clone());
                                    }
                                    if my_approx2 {
                                        let mut tpoles2d: Vec<DVec2> = Vec::new();
                                        mbspc.curve2d(if my_approx1 { 3 } else { 2 }, &mut tpoles2d);
                                        let mut bs2 = Some(bspline2_from_poles(
                                            tpoles2d,
                                            &mbspc.knots,
                                            &mbspc.mults,
                                            mbspc.degree,
                                        ));
                                        geom_lib_check2d_bspline_curve_fix_tangent(
                                            &mut bs2,
                                            self.my_tol_check(),
                                            self.my_tol_ang_check(),
                                        );
                                        //
                                        adjust_u_periodic(&a_s2, &mut bs2);
                                        //
                                        sline_s2.push(bs2);
                                    } else {
                                        sline_s2.push(h1.clone());
                                    }
                                }
                                j += 1;
                            }
                        }
                        i += 1;
                    }
                }
            }

            IntPatchIType::Restriction => {
                // OCCT L1007-1092.
                let rl = l.clone();
                let mut a_c3d: Option<Curve3> = None;
                let mut a_c2d1: Option<Curve2d> = None;
                let mut a_c2d2: Option<Curve2d> = None;
                let mut a_tol_reached = 0.0;
                treat_r_line(
                    &rl,
                    &a_s1,
                    &a_s2,
                    &mut a_c3d,
                    &mut a_c2d1,
                    &mut a_c2d2,
                    &mut a_tol_reached,
                );

                if a_c3d.is_none() {
                    return;
                }
                let a_c3d = a_c3d.unwrap();

                let mut a_box1 = BndBox2d::new();
                let mut a_box2 = BndBox2d::new();

                // OCCT L1021-1029.
                let a_u1f = a_s1.first_u_parameter();
                let a_v1f = a_s1.first_v_parameter();
                let a_u1l = a_s1.last_u_parameter();
                let a_v1l = a_s1.last_v_parameter();
                let a_u2f = a_s2.first_u_parameter();
                let a_v2f = a_s2.first_v_parameter();
                let a_u2l = a_s2.last_u_parameter();
                let a_v2l = a_s2.last_v_parameter();

                a_box1.add_point(DVec2::new(a_u1f, a_v1f));
                a_box1.add_point(DVec2::new(a_u1l, a_v1l));
                a_box2.add_point(DVec2::new(a_u2f, a_v2f));
                a_box2.add_point(DVec2::new(a_u2l, a_v2l));

                let mut an_array_of_parameters: Vec<f64> = Vec::new();

                // We consider here that the intersection line is
                // same-parameter-line.
                an_array_of_parameters.push(curve3_first_parameter(&a_c3d));
                an_array_of_parameters.push(curve3_last_parameter(&a_c3d));

                trim_i_line_on_surf_boundaries(
                    &a_c2d1,
                    &a_c2d2,
                    &a_box1,
                    &a_box2,
                    &mut an_array_of_parameters,
                );

                let a_nb_inters_solutions_m1 = an_array_of_parameters.len() as i64 - 1;

                // Trim RLine found.  OCCT L1042-1090.
                let mut an_ind = 0i64;
                while an_ind < a_nb_inters_solutions_m1 {
                    let a_par_f = an_array_of_parameters[an_ind as usize];
                    let a_par_l = an_array_of_parameters[an_ind as usize + 1];

                    if (a_par_l - a_par_f) <= PCONFUSION {
                        an_ind += 1;
                        continue;
                    }

                    let a_par = 0.5 * (a_par_f + a_par_l);

                    let mut a_curv2d1: Option<Curve2d> = None;
                    let mut a_curv2d2: Option<Curve2d> = None;
                    if let Some(c2d1) = a_c2d1.as_ref() {
                        let a_pt = c2d1.point_at(a_par);

                        if a_box1.is_out_point(a_pt) {
                            an_ind += 1;
                            continue;
                        }

                        if my_approx1 {
                            a_curv2d1 = Some(Curve2d::Trimmed(TrimmedCurve2 {
                                curve: Box::new(c2d1.clone()),
                                t_min: a_par_f,
                                t_max: a_par_l,
                            }));
                        }
                    }

                    if let Some(c2d2) = a_c2d2.as_ref() {
                        let a_pt = c2d2.point_at(a_par);

                        if a_box2.is_out_point(a_pt) {
                            an_ind += 1;
                            continue;
                        }

                        if my_approx2 {
                            a_curv2d2 = Some(Curve2d::Trimmed(TrimmedCurve2 {
                                curve: Box::new(c2d2.clone()),
                                t_min: a_par_f,
                                t_max: a_par_l,
                            }));
                        }
                    }

                    let a_curv3d = Curve3::Trimmed(TrimmedCurve3::new(
                        a_c3d.clone(),
                        a_par_f,
                        a_par_l,
                    ));

                    sline.push(Some(a_curv3d));
                    sline_s1.push(a_curv2d1);
                    sline_s2.push(a_curv2d2);
                    an_ind += 1;
                }
            }

            IntPatchIType::Unknown => {}
        }

        // Write the collected results back into the OCCT member sequences.
        self.sline_mut().extend(sline);
        self.sline_s1_mut().extend(sline_s1);
        self.sline_s2_mut().extend(sline_s2);
        *self.my_tol_reached_2d_mut() = my_tol_reached_2d;
    }
}

/// OCCT ElSLib::PlaneValue(U, V, Pln) = Pln.Location() + U * XDirection +
/// V * YDirection.
fn el_s_lib_plane_value(u: f64, v: f64, pln: &Plane) -> DVec3 {
    pln.origin + pln.u_dir * u + pln.v_dir * v
}

/// OCCT Geom_BSplineCurve::IsClosed().
fn curve3_bspline_is_closed(b: &BSplineCurve3) -> bool {
    if b.control_points.len() < 2 {
        return false;
    }
    b.control_points[0].distance_squared(b.control_points[b.control_points.len() - 1]) < CONFUSION
}

/// OCCT Geom_Curve::FirstParameter of the 3D curve.
fn curve3_first_parameter(c: &Curve3) -> f64 {
    c.default_domain()[0]
}

/// OCCT Geom_Curve::LastParameter of the 3D curve.
fn curve3_last_parameter(c: &Curve3) -> f64 {
    c.default_domain()[1]
}

// ============================================================================
// GAP carriers for the dependencies of this batch that are not translated yet.
// Each carrier is placed exactly where the OCCT statement sits, and follows the
// OCCT branch semantics; the batch report lists them as untranslated.
// ============================================================================

/// OCCT GeomInt_LineTool::DecompositionOfWLine (GeomInt_LineTool.cxx L42-...)
/// — NOT translated in this batch.
///
/// The carrier returns `false`, which is OCCT's own "not decomposited" result:
/// the caller (`MakeCurve`) then iterates `myLConstruct` parts exactly as it
/// does for a WLine that needs no decomposition.
fn decomposition_of_w_line(
    _wl: &IntPatchLine,
    _s1: &GeomSurfaceAdapter,
    _s2: &GeomSurfaceAdapter,
    _tol_ss: f64,
    _l_construct: &GeomIntLineConstructor,
    _a_seq_of_l: &mut Vec<IntPatchLine>,
) -> bool {
    false
}

/// OCCT GeomLib_CheckBSplineCurve::FixTangent(theFirst, theLast)
/// (GeomLib_CheckBSplineCurve.cxx) — NOT translated in this batch.
///
/// The OCCT step has no failure path: it only adjusts the end poles/knots of
/// the approximated B-spline.  The carrier leaves the curve untouched, which is
/// exactly the OCCT result for a curve whose ends need no fixing; the batch
/// report lists this class as untranslated.
fn geom_lib_check_bspline_curve_fix_tangent(
    _bs: &mut Curve3,
    _tol_check: f64,
    _tol_ang_check: f64,
) {
}

/// OCCT GeomLib_Check2dBSplineCurve::FixTangent(theFirst, theLast)
/// (GeomLib_Check2dBSplineCurve.cxx) — NOT translated in this batch; see
/// [`geom_lib_check_bspline_curve_fix_tangent`].
fn geom_lib_check2d_bspline_curve_fix_tangent(
    _bs: &mut Option<Curve2d>,
    _tol_check: f64,
    _tol_ang_check: f64,
) {
}
