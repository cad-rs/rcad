//! OCCT GeomLib statics (TKGeomBase/GeomLib — GeomLib.cxx), the isoline
//! group:
//! - `GeomLib::isIsoLine` (GeomLib.cxx L2991-3078) — the Line /
//!   degree-1-2-pole-BSpline / degree-1-2-pole-Bezier classification plus the
//!   `IsParallel(gp::DX2d(), Precision::Angular())` /
//!   `IsParallel(gp::DY2d(), ...)` return with `theIsU`, `theParam`,
//!   `theIsForward`.
//! - `GeomLib::buildC3dOnIsoLine` (GeomLib.cxx L3079-3235) — the UIso/VIso
//!   construction with the infinite-bound trimming rules, the
//!   `CurveToBSplineCurve(..., Convert_QuasiAngular)` conversion, the
//!   `Reverse()` when `!theIsForward`, the
//!   `BSplCLib::Reparametrize(theC2D->FirstParameter(), theC2D->LastParameter(),
//!   aKnots)` + `SetKnots`, and the 23-sample error evaluation that returns a
//!   null handle when the error exceeds the tolerance.
//!
//! This is GeomLib's own copy of the two statics: OCCT itself carries them
//! twice (once in `GeomLib`, once in `Approx_CurveOnSurface`), and the
//! `Approx_CurveOnSurface` copies live in
//! `crate::geomalgo::approx_curve_on_surface` — they are NOT shared with this
//! module, matching the OCCT duplication.
//!
//! Architecture difference: the `Geom_Surface::UIso` / `VIso` virtual
//! dispatch and the `Geom_RectangularTrimmedSurface` construction are the
//! kernel `Surface3` evaluation set; the rcad `Surface3::Trimmed` payload is
//! the both-trimmed `Geom_RectangularTrimmedSurface`
//! (Geom_RectangularTrimmedSurface.cxx L67-112 sets `isutrimmed = true` and
//! `isvtrimmed = true`), so the UIso/VIso on it apply the complementary trim
//! (Geom_RectangularTrimmedSurface.cxx L444-478).

use glam::DVec2;

use rcad_kernel::base::proj_lib::adaptor::{Adaptor2dCurve2d, Adaptor3dSurface};
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::base::convert::{ConvertParameterisation, geom_convert_curve_to_bspline_curve};
use rcad_kernel::core::precision::{p_confusion, ANGULAR, CONFUSION};
use rcad_kernel::geom::{
    BSplineCurve3, Curve3, CurveEval, Surface3, SurfaceEval, TrimmedCurve3, TrimmedSurface,
};
use rcad_kernel::math::bspl_lib;

/// OCCT GeomLib::isIsoLine(theC2D, theIsU, theParam, theIsForward)
/// (GeomLib.cxx L2991-3078).
pub fn is_iso_line(
    the_c2d: &dyn Adaptor2dCurve2d,
    the_is_u: &mut bool,
    the_param: &mut f64,
    the_is_forward: &mut bool,
) -> bool {
    // These variables are used to check line state (vertical or horizontal).
    let mut is_appropriate_type = false;
    let mut a_loc2d = DVec2::ZERO;
    let mut a_dir2d = DVec2::ZERO;

    // Test type (cxx L3000-3048).
    let a_type = the_c2d.get_type();
    if a_type == CurveType::Line {
        let a_lin2d = the_c2d.line();
        a_loc2d = a_lin2d.origin;
        a_dir2d = a_lin2d.direction;
        is_appropriate_type = true;
    } else if a_type == CurveType::BSpline {
        // OCCT: aBSpline2d = theC2D->BSpline() — the downcast; the rcad
        // trait models the null handle as None.
        let Some(a_bspline2d) = the_c2d.bspline() else {
            return false;
        };
        if a_bspline2d.degree != 1 || a_bspline2d.control_points.len() != 2 {
            return false; // Not a line or uneven parameterization.
        }

        a_loc2d = a_bspline2d.control_points[0];

        // Vector should be non-degenerated.
        let a_vec2d = a_bspline2d.control_points[1] - a_bspline2d.control_points[0];
        if a_vec2d.length_squared() < CONFUSION {
            return false; // Degenerated spline.
        }
        a_dir2d = a_vec2d;

        is_appropriate_type = true;
    } else if a_type == CurveType::Bezier {
        let Some(a_bezier2d) = the_c2d.bezier() else {
            return false;
        };
        // OCCT: aBezier2d->Degree() != 1 — a 2-pole Bezier is degree 1.
        if a_bezier2d.control_points.len() != 2 {
            return false; // Not a line or uneven parameterization.
        }

        a_loc2d = a_bezier2d.control_points[0];

        // Vector should be non-degenerated.
        let a_vec2d = a_bezier2d.control_points[1] - a_bezier2d.control_points[0];
        if a_vec2d.length_squared() < CONFUSION {
            return false; // Degenerated spline.
        }
        a_dir2d = a_vec2d;

        is_appropriate_type = true;
    }

    if !is_appropriate_type {
        return false;
    }

    // Check line to be vertical or horizontal (cxx L3054-3076).
    if dir2d_is_parallel(a_dir2d, DVec2::X, ANGULAR) {
        // Horizontal line. V = const.
        *the_is_u = false;
        *the_param = a_loc2d.y;
        *the_is_forward = a_dir2d.dot(DVec2::X) > 0.0;
        return true;
    } else if dir2d_is_parallel(a_dir2d, DVec2::Y, ANGULAR) {
        // Vertical line. U = const.
        *the_is_u = true;
        *the_param = a_loc2d.x;
        *the_is_forward = a_dir2d.dot(DVec2::Y) > 0.0;
        return true;
    }

    false
}

/// OCCT GeomLib::buildC3dOnIsoLine(theC2D, theSurf, theFirst, theLast,
/// theTolerance, theIsU, theParam, theIsForward) (GeomLib.cxx L3079-3235).
///
/// Returns the OCCT `NewCurvePtr`; the OCCT null-handle outcomes (the failed
/// `GeomAdaptor_Surface` down-cast, the sphere guard, the out-of-bounds /
/// degenerate isoline guards and the error-over-tolerance tail) are `None`.
#[allow(clippy::too_many_arguments)]
pub fn build_c3d_on_iso_line(
    the_c2d: &dyn Adaptor2dCurve2d,
    the_surf: &dyn Adaptor3dSurface,
    the_first: f64,
    the_last: f64,
    the_tolerance: f64,
    the_is_u: bool,
    the_param: f64,
    the_is_forward: bool,
) -> Option<Curve3> {
    // Convert adapter to the appropriate type (cxx L3082-3086): the OCCT
    // down_cast to GeomAdaptor_Surface maps to the kernel-surface bridge.
    let kernel_surf = the_surf.kernel_surface()?;
    let mut a_surf: Surface3 = kernel_surf.clone();

    if the_surf.get_type() == rcad_kernel::base::proj_lib::GeomAbsSurfaceType::Sphere {
        return None;
    }

    // Extract isoline (cxx L3089-3092).
    let a_f2d = the_c2d.value(the_c2d.first_parameter());
    let a_l2d = the_c2d.value(the_c2d.last_parameter());

    let mut is_to_trim = true;
    // OCCT L3095-3096: aSurf->Bounds(U1, U2, V1, V2).
    let domain = a_surf.default_domain();
    let (u1, u2, v1, v2) = (domain[0], domain[1], domain[2], domain[3]);

    let a_c3d: Curve3;
    if the_is_u {
        // cxx L3098-3137.
        let mut a_v1_param = a_f2d.y.min(a_l2d.y);
        let mut a_v2_param = a_f2d.y.max(a_l2d.y);
        if a_v2_param < v1 - the_tolerance || a_v1_param > v2 + the_tolerance {
            return None;
        } else if rcad_kernel::precision::is_infinite_value(v1)
            || rcad_kernel::precision::is_infinite_value(v2)
        {
            if (a_v2_param - a_v1_param).abs() < p_confusion() {
                return None;
            }
            // OCCT L3113: aSurf = new Geom_RectangularTrimmedSurface(aSurf,
            // U1, U2, aV1Param, aV2Param); isToTrim = false.
            a_surf = surface_rectangular_trimmed(&a_surf, u1, u2, a_v1_param, a_v2_param);
            is_to_trim = false;
        } else {
            a_v1_param = a_v1_param.max(v1);
            a_v2_param = a_v2_param.min(v2);
            if (a_v2_param - a_v1_param).abs() < p_confusion() {
                return None;
            }
        }
        // OCCT L3126-3130: aC3d = aSurf->UIso(theParam); if (isToTrim) aC3d =
        // new Geom_TrimmedCurve(aC3d, aV1Param, aV2Param).
        let iso = surface_u_iso(&a_surf, the_param);
        a_c3d = if is_to_trim {
            Curve3::Trimmed(TrimmedCurve3::new(iso, a_v1_param, a_v2_param))
        } else {
            iso
        };
    } else {
        // cxx L3138-3177.
        let mut a_u1_param = a_f2d.x.min(a_l2d.x);
        let mut a_u2_param = a_f2d.x.max(a_l2d.x);
        if a_u2_param < u1 - the_tolerance || a_u1_param > u2 + the_tolerance {
            return None;
        } else if rcad_kernel::precision::is_infinite_value(u1)
            || rcad_kernel::precision::is_infinite_value(u2)
        {
            if (a_u2_param - a_u1_param).abs() < p_confusion() {
                return None;
            }
            // OCCT L3153: aSurf = new Geom_RectangularTrimmedSurface(aSurf,
            // aU1Param, aU2Param, V1, V2); isToTrim = false.
            a_surf = surface_rectangular_trimmed(&a_surf, a_u1_param, a_u2_param, v1, v2);
            is_to_trim = false;
        } else {
            a_u1_param = a_u1_param.max(u1);
            a_u2_param = a_u2_param.min(u2);
            if (a_u2_param - a_u1_param).abs() < p_confusion() {
                return None;
            }
        }
        // OCCT L3166-3170: aC3d = aSurf->VIso(theParam); if (isToTrim) aC3d =
        // new Geom_TrimmedCurve(aC3d, aU1Param, aU2Param).
        let iso = surface_v_iso(&a_surf, the_param);
        a_c3d = if is_to_trim {
            Curve3::Trimmed(TrimmedCurve3::new(iso, a_u1_param, a_u2_param))
        } else {
            iso
        };
    }

    // Convert arbitrary curve type to the b-spline (cxx L3180-3185).
    let mut a_curve3d = geom_convert_curve_to_bspline(&a_c3d, ConvertParameterisation::QuasiAngular);
    if !the_is_forward {
        // OCCT: aCurve3d->Reverse().
        a_curve3d = a_curve3d.reversed();
    }

    // Rebuild parameterization for the 3d curve to have the same
    // parameterization with a two-dimensional curve (cxx L3189-3191):
    // BSplCLib::Reparametrize(theC2D->FirstParameter(),
    // theC2D->LastParameter(), aKnots); aCurve3d->SetKnots(aKnots).
    let mut a_knots = a_curve3d.knots.clone();
    bspl_lib::reparametrize(
        the_c2d.first_parameter(),
        the_c2d.last_parameter(),
        &mut a_knots,
    );
    a_curve3d.knots = a_knots;

    // Evaluate error (cxx L3193-3214).
    let mut an_error3d = 0.0f64;

    let a_par_f = the_first;
    let a_par_l = the_last;
    let a_nb_pnt = 23usize;
    for an_idx in 0..=a_nb_pnt {
        let a_par = a_par_f + ((a_par_l - a_par_f) * an_idx as f64) / a_nb_pnt as f64;

        let a_pnt2d = the_c2d.value(a_par);

        let a_pnt_c3d = a_curve3d.point_at(a_par);
        let a_pnt_c2d = the_surf.value(a_pnt2d.x, a_pnt2d.y);

        let a_sq_deviation = a_pnt_c3d.distance_squared(a_pnt_c2d);
        an_error3d = a_sq_deviation.max(an_error3d);
    }

    an_error3d = an_error3d.sqrt();

    // Target tolerance is not obtained. This situation happens for isolines on
    // the sphere. OCCT is unable to convert it keeping original
    // parameterization, while the geometric form of the result is entirely
    // identical. In that case, it is better to utilize a general-purpose
    // approach (cxx L3218-3232).
    if an_error3d > the_tolerance {
        return None;
    }

    Some(Curve3::BSpline(a_curve3d))
}

// ---------------------------------------------------------------------------
// Geom_Surface virtual dispatch support (the UIso / VIso / rectangular-trim
// constructions the GeomLib static calls).
// ---------------------------------------------------------------------------

/// OCCT Geom_RectangularTrimmedSurface(S, U1, U2, V1, V2, USense, VSense)
/// (Geom_RectangularTrimmedSurface.cxx L67-112) — the nested trimmed basis is
/// killed and the resulting surface carries `isutrimmed = isvtrimmed = true`.
fn surface_rectangular_trimmed(surf: &Surface3, u1: f64, u2: f64, v1: f64, v2: f64) -> Surface3 {
    // OCCT: kill trimmed basis surfaces.
    let basis = match surf {
        Surface3::Trimmed(t) => (*t.basis).clone(),
        other => other.clone(),
    };
    Surface3::Trimmed(TrimmedSurface::new(basis, u1, u2, v1, v2))
}

/// OCCT Geom_Surface::UIso (the per-type overrides) + the
/// Geom_RectangularTrimmedSurface::UIso complementary-V-trim
/// (Geom_RectangularTrimmedSurface.cxx L444-459).
fn surface_u_iso(surf: &Surface3, param: f64) -> Curve3 {
    use rcad_kernel::base::proj_lib::elslib_iso as el;
    match surf {
        Surface3::Trimmed(t) => {
            // OCCT L447-458: C = basisSurf->UIso(U); if (isvtrimmed) return
            // new Geom_TrimmedCurve(C, vtrim1, vtrim2, true).
            let c = surface_u_iso(&t.basis, param);
            Curve3::Trimmed(TrimmedCurve3::new(c, t.trim[2], t.trim[3]))
        }
        Surface3::Plane(p) => Curve3::Line(el::elslib_plane_u_iso(
            &el::Ax3View::from_axes(p.origin, p.normal, p.u_dir),
            param,
        )),
        Surface3::Cylinder(c) => Curve3::Line(el::elslib_cylinder_u_iso(
            &el::Ax3View::from_axes(c.origin, c.axis, c.ref_dir),
            c.radius,
            param,
        )),
        Surface3::Cone(c) => Curve3::Line(el::elslib_cone_u_iso(
            &el::Ax3View::from_axes(c.apex, c.axis, c.ref_dir),
            c.radius,
            c.half_angle_rad,
            param,
        )),
        Surface3::Sphere(s) => Curve3::Circle(el::elslib_sphere_u_iso(
            &el::Ax3View::from_axes(s.center, s.axis, s.ref_dir),
            s.radius,
            param,
        )),
        Surface3::Torus(t) => Curve3::Circle(el::elslib_torus_u_iso(
            &el::Ax3View::from_axes(t.center, t.axis, t.ref_dir),
            t.major_radius,
            t.minor_radius,
            param,
        )),
        // OCCT Geom_BSplineSurface::UIso (Geom_BSplineSurface_1.cxx L598-635):
        // BSplSLib::Iso on the U direction; the result curve's knots and
        // degree are the V direction's (and its periodicity myVPeriodic).
        Surface3::BSpline(b) => {
            let weights = if b.is_rational_u() || b.is_rational_v() {
                Some(b.weights.as_slice())
            } else {
                None
            };
            let (cpoles, cweights) = bspl_lib::bspl_slib_iso(
                param,
                true,
                b.degree_u,
                &b.knots_u,
                &b.control_points,
                weights,
                b.is_periodic_u,
            );
            Curve3::BSpline(BSplineCurve3 {
                degree: b.degree_v,
                knots: b.knots_v.clone(),
                control_points: cpoles,
                weights: cweights,
                is_periodic: b.is_periodic_v,
            })
        }
        _ => panic!("GAP: Geom_Surface::UIso not translated for this surface type"),
    }
}

/// OCCT Geom_Surface::VIso (the per-type overrides) + the
/// Geom_RectangularTrimmedSurface::VIso complementary-U-trim
/// (Geom_RectangularTrimmedSurface.cxx L463-478).
fn surface_v_iso(surf: &Surface3, param: f64) -> Curve3 {
    use rcad_kernel::base::proj_lib::elslib_iso as el;
    match surf {
        Surface3::Trimmed(t) => {
            // OCCT L466-477: C = basisSurf->VIso(V); if (isutrimmed) return
            // new Geom_TrimmedCurve(C, utrim1, utrim2, true).
            let c = surface_v_iso(&t.basis, param);
            Curve3::Trimmed(TrimmedCurve3::new(c, t.trim[0], t.trim[1]))
        }
        Surface3::Plane(p) => Curve3::Line(el::elslib_plane_v_iso(
            &el::Ax3View::from_axes(p.origin, p.normal, p.u_dir),
            param,
        )),
        Surface3::Cylinder(c) => Curve3::Circle(el::elslib_cylinder_v_iso(
            &el::Ax3View::from_axes(c.origin, c.axis, c.ref_dir),
            c.radius,
            param,
        )),
        Surface3::Cone(c) => Curve3::Circle(el::elslib_cone_v_iso(
            &el::Ax3View::from_axes(c.apex, c.axis, c.ref_dir),
            c.radius,
            c.half_angle_rad,
            param,
        )),
        Surface3::Sphere(s) => Curve3::Circle(el::elslib_sphere_v_iso(
            &el::Ax3View::from_axes(s.center, s.axis, s.ref_dir),
            s.radius,
            param,
        )),
        Surface3::Torus(t) => Curve3::Circle(el::elslib_torus_v_iso(
            &el::Ax3View::from_axes(t.center, t.axis, t.ref_dir),
            t.major_radius,
            t.minor_radius,
            param,
        )),
        // OCCT Geom_BSplineSurface::VIso (Geom_BSplineSurface_1.cxx L775-812).
        Surface3::BSpline(b) => {
            let weights = if b.is_rational_u() || b.is_rational_v() {
                Some(b.weights.as_slice())
            } else {
                None
            };
            let (cpoles, cweights) = bspl_lib::bspl_slib_iso(
                param,
                false,
                b.degree_v,
                &b.knots_v,
                &b.control_points,
                weights,
                b.is_periodic_v,
            );
            Curve3::BSpline(BSplineCurve3 {
                degree: b.degree_u,
                knots: b.knots_u.clone(),
                control_points: cpoles,
                weights: cweights,
                is_periodic: b.is_periodic_u,
            })
        }
        _ => panic!("GAP: Geom_Surface::VIso not translated for this surface type"),
    }
}

/// OCCT GeomConvert::CurveToBSplineCurve(C, Convert_QuasiAngular)
/// (GeomConvert.cxx L157-209).
///
/// The exact kernel conversion has landed the Trimmed(Line) and
/// Trimmed(Circle) arms (GeomConvert.cxx L200-282).  This isoline path wraps
/// the `Geom_Surface::UIso/VIso` result in a trimmed curve or hands the raw
/// iso curve over, so the trimmed line / circle arms are the reachable ones;
/// inputs outside the landed arms keep the kernel's sampling conversion (the
/// pre-existing stand-in for the staged arms).
fn geom_convert_curve_to_bspline(c: &Curve3, parameterisation: ConvertParameterisation) -> BSplineCurve3 {
    if let Curve3::Trimmed(tc) = c {
        let mut basis = tc.basis_curve();
        let mut depth = 0;
        while let Curve3::Trimmed(inner) = basis {
            basis = inner.basis_curve();
            depth += 1;
            if depth > 8 {
                break;
            }
        }
        if matches!(basis, Curve3::Line(_) | Curve3::Circle(_)) {
            return geom_convert_curve_to_bspline_curve(
                &Curve3::Trimmed(TrimmedCurve3::new(basis.clone(), tc.first, tc.last)),
                parameterisation,
            );
        }
    }
    rcad_kernel::base::convert::curve_to_bspline(c, 23)
}

/// OCCT gp_Dir2d::IsParallel(Other, AngularTolerance) — true when the angular
/// distance between the two directions is within the tolerance (mod pi).
fn dir2d_is_parallel(a: DVec2, b: DVec2, angular_tolerance: f64) -> bool {
    let ang_a = a.y.atan2(a.x);
    let ang_b = b.y.atan2(b.x);
    let diff = (ang_a - ang_b).abs();
    let diff = diff.min(std::f64::consts::PI * 2.0 - diff);
    diff <= angular_tolerance || (diff - std::f64::consts::PI).abs() <= angular_tolerance
}
