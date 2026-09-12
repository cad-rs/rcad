//! The curve-construction half of the `GeomInt_IntSS` engine: MakeBSpline /
//! MakeBSpline2d, BuildPCurves (both overloads), TrimILineOnSurfBoundaries,
//! TreatRLine and the RLine parameter helpers.
//!
//! Split out of `geom_int_int_ss_1.rs` verbatim (the 2000-line file limit,
//! AGENTS.md Rule 5) and registered as a sibling module in `geomalgo/mod.rs`.
//! The seven file-local helpers this half calls back into
//! (`epsilon_of`, `curve2d_first_parameter`, `curve2d_last_parameter`,
//! `parameters_of_nearest_point_on_surface`, `intersect_curve_and_boundary`,
//! `is_degenerated`) stay in `geom_int_int_ss_1.rs` and are widened to
//! `pub(super)` for this import only — bodies unchanged.

use glam::{DVec2, DVec3};

use rcad_kernel::base::geom_proj_lib;
use rcad_kernel::geom::{
    BSplineCurve2, BSplineCurve3, Curve2d, Curve2dEval, Curve3, CurveEval, SurfaceEval,
    TrimmedCurve2,
};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::math::bspl_lib::reparametrize;
use rcad_kernel::precision::{is_infinite_value, CONFUSION, PCONFUSION};

use super::geom_int_int_ss_1::{
    curve2d_first_parameter, curve2d_last_parameter, epsilon_of, intersect_curve_and_boundary,
    is_degenerated, parameters_of_nearest_point_on_surface,
};
use crate::geomalgo::int_patch::IntPatchLine;
use crate::hlr::contap::surface_adaptor::GeomSurfaceAdapter;

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L1452-1502: MakeBSpline / MakeBSpline2d
// ============================================================================

/// OCCT GeomInt_IntSS::MakeBSpline(WL, ideb, ifin) (cxx L1452-1469).
pub(crate) fn make_b_spline(wl: &IntPatchLine, ideb: i32, ifin: i32) -> Option<Curve3> {
    let nbpnt = (ifin - ideb + 1) as usize;
    let mut poles: Vec<DVec3> = Vec::with_capacity(nbpnt);
    let mut knots: Vec<f64> = Vec::with_capacity(nbpnt);
    let mut mults: Vec<usize> = Vec::with_capacity(nbpnt);
    // int i = 1, ipidebm1 = ideb;
    // for (; i <= nbpnt; ipidebm1++, i++) { poles(i) = WL->Point(ipidebm1).Value(); ... }
    let mut i = 1i32;
    let mut ipidebm1 = ideb;
    while i <= nbpnt as i32 {
        // WL->Point(ipidebm1) — OCCT 1-based, rcad 0-based.
        let p = &wl.wline_pnts[(ipidebm1 - 1) as usize];
        poles.push(p.p3d);
        mults.push(1);
        knots.push(i as f64 - 1.);
        ipidebm1 += 1;
        i += 1;
    }
    // mults(1) = mults(nbpnt) = 2;
    mults[0] = 2;
    mults[nbpnt - 1] = 2;
    Some(bspline3_from_poles(poles, &knots, &mults, 1))
}

/// OCCT GeomInt_IntSS::MakeBSpline2d(theWLine, ideb, ifin, onFirst) (cxx
/// L1473-1502).
pub(crate) fn make_b_spline_2d(
    the_w_line: &IntPatchLine,
    ideb: i32,
    ifin: i32,
    on_first: bool,
) -> Option<Curve2d> {
    let nbpnt = (ifin - ideb + 1) as usize;
    let mut poles: Vec<DVec2> = Vec::with_capacity(nbpnt);
    let mut knots: Vec<f64> = Vec::with_capacity(nbpnt);
    let mut mults: Vec<usize> = Vec::with_capacity(nbpnt);
    let mut i = 1i32;
    let mut ipidebm1 = ideb;
    while i <= nbpnt as i32 {
        let p = &the_w_line.wline_pnts[(ipidebm1 - 1) as usize];
        let (u, v) = if on_first {
            (p.u1, p.v1)
        } else {
            (p.u2, p.v2)
        };
        poles.push(DVec2::new(u, v));
        mults.push(1);
        knots.push(i as f64 - 1.);
        ipidebm1 += 1;
        i += 1;
    }

    mults[0] = 2;
    mults[nbpnt - 1] = 2;
    Some(bspline2_from_poles(poles, &knots, &mults, 1))
}

/// `new Geom_BSplineCurve(poles, knots, mults, degree)` — the NCollection
/// arrays carry the COMPRESSED (knots, multiplicities) pair, while the rcad
/// `BSplineCurve3` stores the expanded knot vector.
pub(crate) fn bspline3_from_poles(
    poles: Vec<DVec3>,
    knots: &[f64],
    mults: &[usize],
    degree: usize,
) -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree,
        knots: expand_knots(knots, mults),
        control_points: poles,
        weights: vec![],
        is_periodic: false,
    })
}

/// `new Geom2d_BSplineCurve(poles, knots, mults, degree)`.
pub(crate) fn bspline2_from_poles(
    poles: Vec<DVec2>,
    knots: &[f64],
    mults: &[usize],
    degree: usize,
) -> Curve2d {
    Curve2d::BSpline(BSplineCurve2 {
        degree,
        knots: expand_knots(knots, mults),
        control_points: poles,
        weights: vec![],
    })
}

/// Expand the compressed (knots, multiplicities) representation into the full
/// knot vector.
fn expand_knots(knots: &[f64], mults: &[usize]) -> Vec<f64> {
    let mut out = Vec::new();
    for (i, &k) in knots.iter().enumerate() {
        let m = if i < mults.len() { mults[i] } else { 1 };
        for _ in 0..m {
            out.push(k);
        }
    }
    out
}

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L1172-1325: BuildPCurves
// ============================================================================

/// OCCT GeomInt_IntSS::BuildPCurves(f, l, umin, umax, vmin, vmax, tol, S, C,
/// C2d) (cxx L1172-1304) — the full-bounds overload.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_p_curves(
    the_first: f64,
    the_last: f64,
    the_umin: f64,
    the_umax: f64,
    the_vmin: f64,
    the_vmax: f64,
    the_tol: &mut f64,
    the_surface: &rcad_kernel::geom::Surface3,
    the_curve: &Curve3,
    the_curve2d: &mut Option<Curve2d>,
) {
    // OCCT L1183-1186: if (!theCurve2d.IsNull() || theSurface.IsNull()) return;
    if the_curve2d.is_some() {
        return;
    }
    //
    // in class ProjLib_Function the range of parameters is shrank by 1.e-09
    if (the_last - the_first) > 2.0e-09 {
        // OCCT L1191-1199: theCurve2d = GeomProjLib::Curve2d(theCurve,
        // theFirst, theLast, theSurface, theUmin, theUmax, theVmin, theVmax,
        // theTol);  (the TolReached goes through the reference parameter)
        let mut a_tol = *the_tol;
        let projected = geom_proj_lib::curve2d(
            the_curve,
            the_first,
            the_last,
            the_surface,
            the_umin,
            the_umax,
            the_vmin,
            the_vmax,
        );
        // rcad note: `geom_proj_lib::curve2d` does not report a TolReached
        // (the kernel translation returns Option<Curve2d>); the OCCT
        // `theTol` in/out value is therefore unchanged.  See the batch report.
        let _ = &mut a_tol;
        *the_curve2d = projected;
        if the_curve2d.is_none() {
            // proj. a circle that goes through the pole on a sphere to the
            // sphere: theTol += Precision::Confusion();
            *the_tol += CONFUSION;
            *the_curve2d =
                geom_proj_lib::curve2d_simple(the_curve, the_first, the_last, the_surface);
        }
        // OCCT L1206-1221: if the result is a Geom2d_BSplineCurve, re-check the
        // first/last knots.
        if let Some(c2d) = the_curve2d.as_mut() {
            if matches!(c2d, Curve2d::BSpline(_)) {
                if (curve2d_first_parameter(c2d) - the_first > PCONFUSION)
                    || (the_last - curve2d_last_parameter(c2d) > PCONFUSION)
                {
                    // NCollection_Array1<double> aKnots(aBspl->Knots());
                    // BSplCLib::Reparametrize(theFirst, theLast, aKnots);
                    // aBspl->SetKnots(aKnots);
                    if let Curve2d::BSpline(bspl) = c2d {
                        let mut a_knots = bspl.knots.clone();
                        reparametrize(the_first, the_last, &mut a_knots);
                        bspl.knots = a_knots;
                    }
                }
            }
        }
    } else {
        // OCCT L1223-1281.
        if (the_last - the_first) > epsilon_of(the_first.abs()) {
            // The domain of C2d is [Epsilon(std::abs(f)), 2.e-09]; on this
            // small range C2d can be considered as segment of line.
            let mut a_u = 0.;
            let mut a_v = 0.;
            let a_p3d1 = the_curve.point_at(the_first);
            let a_p3d2 = the_curve.point_at(the_last);

            // OCCT L1232-1245: GeomAdaptor_Surface anAS; anAS.Load(theSurface);
            // Extrema_ExtPS anExtr; anExtr.SetAlgo(Extrema_ExtAlgo_Grad);
            // anExtr.Initialize(anAS, theUmin, theUmax, theVmin, theVmax,
            // Precision::Confusion(), Precision::Confusion());
            // anExtr.Perform(aP3d1);
            //
            // rcad: the kernel Extrema ExtPS performs the Grad Initialize +
            // Perform pair (kernel base/extrema.rs ExtPS::perform), so the
            // adaptor load and the two-phase init collapse into one call.
            let mut an_extr = rcad_kernel::base::extrema::ExtPS::new(
                a_p3d1,
                the_surface,
                CONFUSION,
                CONFUSION,
            );
            let mut an_extr = {
                an_extr.perform(
                    a_p3d1,
                    the_surface,
                    the_umin,
                    the_umax,
                    the_vmin,
                    the_vmax,
                    CONFUSION,
                    CONFUSION,
                );
                an_extr
            };

            if parameters_of_nearest_point_on_surface(&an_extr, &mut a_u, &mut a_v) {
                let a_p2d1 = DVec2::new(a_u, a_v);

                an_extr.perform(
                    a_p3d2,
                    the_surface,
                    the_umin,
                    the_umax,
                    the_vmin,
                    the_vmax,
                    CONFUSION,
                    CONFUSION,
                );

                if parameters_of_nearest_point_on_surface(&an_extr, &mut a_u, &mut a_v) {
                    let a_p2d2 = DVec2::new(a_u, a_v);

                    if a_p2d1.distance(a_p2d2) > rcad_kernel::precision::CONFUSION {
                        let poles = vec![a_p2d1, a_p2d2];
                        let knots = vec![the_first, the_last];
                        let mults = vec![2usize, 2usize];

                        *the_curve2d = Some(bspline2_from_poles(poles, &knots, &mults, 1));

                        // Check same parameter in middle point .begin
                        let p_mid =
                            the_curve.point_at(0.5 * (the_first + the_last));
                        let pmid_curve2d = DVec2::new(
                            0.5 * (a_p2d1.x + a_p2d2.x),
                            0.5 * (a_p2d1.y + a_p2d2.y),
                        );
                        let a_pc = the_surface.point_at(pmid_curve2d.x, pmid_curve2d.y);
                        let a_dist = p_mid.distance(a_pc);
                        *the_tol = a_dist.max(*the_tol);
                        // Check same parameter in middle point .end
                    }
                }
            }
        }
    }
    //
    // OCCT L1284-1303: recadre dans le domaine UV de la face.
    if the_surface.is_u_periodic() {
        if let Some(c2d) = the_curve2d.as_mut() {
            let a_eps = PCONFUSION;
            let period = period_of_u(the_surface);
            //
            let a_tm = 0.5 * (the_first + the_last);
            let pm = c2d.point_at(a_tm);
            let u0 = pm.x;
            //
            let mut u0x = 0.0;
            let mut du = 0.0;
            let b_adjust = crate::geomalgo::geom_int_line_constructor::geom_int_adjust_periodic(
                u0, the_umin, the_umax, period, &mut u0x, &mut du, a_eps,
            );
            if b_adjust {
                // gp_Vec2d T1(du, 0.); theCurve2d->Translate(T1);
                let translated = rcad_kernel::geom::translate_curve2d(c2d, DVec2::new(du, 0.));
                *the_curve2d = Some(translated);
            }
        }
    }
}

/// OCCT GeomInt_IntSS::BuildPCurves(f, l, Tol, S, C, C2d) (cxx L1308-1325) —
/// the surface-bounds overload.
pub(crate) fn build_p_curves_bounds(
    f: f64,
    l: f64,
    tol: &mut f64,
    s: &rcad_kernel::geom::Surface3,
    c: &Curve3,
    c2d: &mut Option<Curve2d>,
) {
    // OCCT L1315-1318.
    if c2d.is_some() {
        return;
    }
    // OCCT L1320-1322: S->Bounds(umin, umax, vmin, vmax).
    let (umin, umax, vmin, vmax) = surface_bounds(s);

    build_p_curves(f, l, umin, umax, vmin, vmax, tol, s, c, c2d);
}

/// OCCT Geom_Surface::Bounds (the adaptor-free bounds of the carried surface).
fn surface_bounds(s: &rcad_kernel::geom::Surface3) -> (f64, f64, f64, f64) {
    let d = rcad_kernel::geom::SurfaceEval::default_domain(s);
    (d[0], d[1], d[2], d[3])
}

/// OCCT Geom_Surface::UPeriod().
fn period_of_u(s: &rcad_kernel::geom::Surface3) -> f64 {
    if rcad_kernel::geom::SurfaceEval::is_u_periodic(s) {
        std::f64::consts::TAU
    } else {
        0.0
    }
}

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L1333-1448: TrimILineOnSurfBoundaries
// ============================================================================

/// OCCT GeomInt_IntSS::TrimILineOnSurfBoundaries(theC2d1, theC2d2, theBound1,
/// theBound2, theArrayOfParameters) (cxx L1333-1448).
pub(crate) fn trim_i_line_on_surf_boundaries(
    the_c2d1: &Option<Curve2d>,
    the_c2d2: &Option<Curve2d>,
    the_bound1: &BndBox2d,
    the_bound2: &BndBox2d,
    the_array_of_parameters: &mut Vec<f64>,
) {
    // Rectangular boundaries of two surfaces:
    // [0]:U=Ufirst, [1]:U=Ulast, [2]:V=Vfirst, [3]:V=Vlast
    const A_NUMBER_OF_CURVES: usize = 4;
    let mut a_cur_s1_bounds: [Option<Curve2d>; A_NUMBER_OF_CURVES] = [None, None, None, None];
    let mut a_cur_s2_bounds: [Option<Curve2d>; A_NUMBER_OF_CURVES] = [None, None, None, None];

    // OCCT L1346-1350.
    let (a_u1f, a_v1f, a_u1l, a_v1l) = the_bound1.get().expect("Bnd_Box2d::Get");
    let (a_u2f, a_v2f, a_u2l, a_v2l) = the_bound2.get().expect("Bnd_Box2d::Get");

    // OCCT L1352-1373.
    let mut a_delta = a_v1l - a_v1f;
    if a_delta.abs() > real_small() {
        if !is_infinite_value(a_u1f) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u1f, a_v1f),
                DVec2::new(0., 1.),
            ));

            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s1_bounds[0] = Some(c);
        }

        if !is_infinite_value(a_u1l) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u1l, a_v1f),
                DVec2::new(0., 1.),
            ));
            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s1_bounds[1] = Some(c);
        }
    }

    // OCCT L1375-1395.
    a_delta = a_u1l - a_u1f;
    if a_delta.abs() > real_small() {
        if !is_infinite_value(a_v1f) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u1f, a_v1f),
                DVec2::new(1., 0.),
            ));
            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s1_bounds[2] = Some(c);
        }

        if !is_infinite_value(a_v1l) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u1f, a_v1l),
                DVec2::new(1., 0.),
            ));
            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s1_bounds[3] = Some(c);
        }
    }

    // OCCT L1397-1417.
    a_delta = a_v2l - a_v2f;
    if a_delta.abs() > real_small() {
        if !is_infinite_value(a_u2f) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u2f, a_v2f),
                DVec2::new(0., 1.),
            ));
            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s2_bounds[0] = Some(c);
        }

        if !is_infinite_value(a_u2l) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u2l, a_v2f),
                DVec2::new(0., 1.),
            ));
            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s2_bounds[1] = Some(c);
        }
    }

    // OCCT L1419-1439.
    a_delta = a_u2l - a_u2f;
    if a_delta.abs() > real_small() {
        if !is_infinite_value(a_v2f) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u2f, a_v2f),
                DVec2::new(1., 0.),
            ));
            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s2_bounds[2] = Some(c);
        }

        if !is_infinite_value(a_v2l) {
            let mut c = Curve2d::Line(rcad_kernel::geom::Line2d::new(
                DVec2::new(a_u2f, a_v2l),
                DVec2::new(1., 0.),
            ));
            if !is_infinite_value(a_delta) {
                c = Curve2d::Trimmed(TrimmedCurve2 {
                    curve: Box::new(c),
                    t_min: 0.,
                    t_max: a_delta,
                });
            }
            a_cur_s2_bounds[3] = Some(c);
        }
    }

    // OCCT L1441-1447.
    let an_int_tol = 10.0 * CONFUSION;

    intersect_curve_and_boundary(
        the_c2d1,
        &a_cur_s1_bounds,
        A_NUMBER_OF_CURVES,
        an_int_tol,
        the_array_of_parameters,
    );

    intersect_curve_and_boundary(
        the_c2d2,
        &a_cur_s2_bounds,
        A_NUMBER_OF_CURVES,
        an_int_tol,
        the_array_of_parameters,
    );

    the_array_of_parameters.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
}

/// OCCT RealSmall() = 1.0e-10 (Standard_Real.hxx L87-91).
fn real_small() -> f64 {
    1.0e-10
}

// ============================================================================
// OCCT GeomInt_IntSS_1.cxx L1098-1168: TreatRLine
// ============================================================================

/// OCCT GeomInt_IntSS::TreatRLine(theRL, theHS1, theHS2, theC3d, theC2d1,
/// theC2d2, theTolReached) (cxx L1098-1168).
pub(crate) fn treat_r_line(
    the_rl: &IntPatchLine,
    the_hs1: &GeomSurfaceAdapter,
    the_hs2: &GeomSurfaceAdapter,
    the_c3d: &mut Option<Curve3>,
    the_c2d1: &mut Option<Curve2d>,
    the_c2d2: &mut Option<Curve2d>,
    the_tol_reached: &mut f64,
) {
    // OCCT L1106-1108.
    let a_gahs: GeomSurfaceAdapter;
    let an_ahc2d: Curve2d;
    let tf: f64;
    let tl: f64;

    // It is assumed that 2d curve is 2d line (rectangular surface domain)
    if the_rl.is_arc_on_s1() {
        // OCCT L1112-1118.
        a_gahs = the_hs1.clone();
        an_ahc2d = the_rl.arc_on_s1().expect("RLine::ArcOnS1").clone();
        let (a_tf, a_tl) = param_on_s1(the_rl);
        // theC2d1 = Geom2dAdaptor::MakeCurve(*anAHC2d); — rcad: the RLine
        // carries the arc as a Geom2d_Curve already (architecture note in the
        // file header), hence the identity.
        tf = a_tf.max(curve2d_first_parameter(&an_ahc2d));
        tl = a_tl.min(curve2d_last_parameter(&an_ahc2d));
        // theC2d1 = new Geom2d_TrimmedCurve(theC2d1, tf, tl);
        *the_c2d1 = Some(Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(an_ahc2d.clone()),
            t_min: tf,
            t_max: tl,
        }));
    } else if the_rl.is_arc_on_s2() {
        // OCCT L1120-1129.
        a_gahs = the_hs2.clone();
        an_ahc2d = the_rl.arc_on_s2().expect("RLine::ArcOnS2").clone();
        let (a_tf, a_tl) = param_on_s2(the_rl);
        tf = a_tf.max(curve2d_first_parameter(&an_ahc2d));
        tl = a_tl.min(curve2d_last_parameter(&an_ahc2d));
        *the_c2d2 = Some(Curve2d::Trimmed(TrimmedCurve2 {
            curve: Box::new(an_ahc2d.clone()),
            t_min: tf,
            t_max: tl,
        }));
    } else {
        // OCCT L1130-1133: return;
        return;
    }

    // Restriction line can correspond to a degenerated edge.  In this case we
    // return null-curve.
    if is_degenerated(&a_gahs, &an_ahc2d, tf, tl) {
        return;
    }

    //
    // To provide sameparameter it is necessary to get 3d curve as
    // approximation of curve on surface.
    let a_max_deg = 8;
    let a_max_seg = 1000;
    // OCCT L1147-1148:
    //   Approx_CurveOnSurface anApp(anAHC2d, aGAHS, tf, tl, Precision::Confusion());
    //   anApp.Perform(aMaxSeg, aMaxDeg, GeomAbs_C1, true, false);
    let c2d_handle = std::sync::Arc::new(
        rcad_kernel::base::proj_lib::adaptor::Geom2dCurveAdaptor::new(an_ahc2d.clone()),
    );
    let surf_handle = std::sync::Arc::new(
        rcad_kernel::base::proj_lib::geom_adaptor_surface::GeomSurfaceAdaptor::new(
            a_gahs.surface3().clone(),
        ),
    );
    let mut an_app = crate::geomalgo::approx_curve_on_surface::ApproxCurveOnSurface::new(
        c2d_handle,
        surf_handle,
        tf,
        tl,
        CONFUSION,
    );
    an_app.perform(a_max_seg, a_max_deg, rcad_kernel::math::GeomAbsShape::C1, true, false);
    if !an_app.has_result() {
        // OCCT L1149-1152: if (!anApp.HasResult()) return;
        return;
    }

    // OCCT L1154-1155.
    *the_c3d = an_app.curve3d().map(Curve3::BSpline);
    *the_tol_reached = an_app.max_error3d();
    let mut a_tol = CONFUSION;
    // OCCT L1157-1161.
    if the_rl.is_arc_on_s1() {
        // GeomAdaptor::MakeSurface(*theHS2)
        let a_s = the_hs2.surface3().clone();
        let c3d = the_c3d.clone().expect("TreatRLine theC3d");
        build_p_curves_bounds(tf, tl, &mut a_tol, &a_s, &c3d, the_c2d2);
    }
    // OCCT L1162-1166.
    if the_rl.is_arc_on_s2() {
        let a_s = the_hs1.surface3().clone();
        let c3d = the_c3d.clone().expect("TreatRLine theC3d");
        build_p_curves_bounds(tf, tl, &mut a_tol, &a_s, &c3d, the_c2d1);
    }
    // OCCT L1167: theTolReached = std::max(theTolReached, aTol);
    *the_tol_reached = the_tol_reached.max(a_tol);
}

/// OCCT IntPatch_RLine::ParamOnS1(tf, tl) — the arc range on surface 1 stored
/// on the RLine.
fn param_on_s1(rl: &IntPatchLine) -> (f64, f64) {
    if let Some(v) = rl.vertices.first() {
        if let Some(l) = rl.vertices.last() {
            return (v.param_on_arc1, l.param_on_arc1);
        }
    }
    (f64::NEG_INFINITY, f64::INFINITY)
}

/// OCCT IntPatch_RLine::ParamOnS2(tf, tl).
fn param_on_s2(rl: &IntPatchLine) -> (f64, f64) {
    if let Some(v) = rl.vertices.first() {
        if let Some(l) = rl.vertices.last() {
            return (v.param_on_arc2, l.param_on_arc2);
        }
    }
    (f64::NEG_INFINITY, f64::INFINITY)
}
