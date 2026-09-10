//! OCCT CPnts_AbscissaPoint + GCPnts_AbscissaPoint (TKGeomBase/CPnts,
//! TKGeomBase/GCPnts) — the arc-length machinery of the Length statics.
//!
//! Delivered 1:1 (the members consumed by the GeomPlate curve path —
//! GeomPlate_CurveConstraint::Length and GeomPlate_BuildPlateSurface::Discretise):
//! - CPnts_AbscissaPoint.cxx L41-99: the f3d/f2d integrands and the two
//!   `order` dispatchers,
//! - CPnts_AbscissaPoint.cxx L131-207: the four `Length` statics,
//! - math_GaussSingleIntegration.cxx L55-155: the constructor integration
//!   drivers (the Tol overload is delivered for the Length+Tol forms),
//! - GCPnts_AbscissaPoint.cxx L26-65: `computeType`,
//! - GCPnts_AbscissaPoint.cxx L372-423: the private `length` template,
//! - GCPnts_AbscissaPoint.cxx L305-367: the public `Length` statics.
//!
//! Not delivered (unconsumed by the plate pipeline; they require
//! CPnts_MyRootFunction + math_FunctionRoot): the CPnts_AbscissaPoint
//! instance machinery (Init/Perform/AdvPerform) and the GCPnts
//! constructors (compute/advCompute drivers).
//!
//! The template parameter `TheCurve` is [`GCPntsCurve`] (see
//! [`super::gcpnts_curve`]); the 2D flavor is [`GCPntsCurve2d`].  In OCCT
//! the 3d/2d flavor is resolved by C++ overload resolution on the adaptor
//! static type; here it is carried by the `*_3d` / `*_2d` function pairs.

use glam::DVec3;

use rcad_kernel::math::gauss_points::{gauss_points, gauss_points_max, gauss_weights};
use rcad_kernel::math::VecD;

use super::gcpnts_curve::GCPntsCurve;

/// OCCT GCPnts_AbscissaType.hxx — compute/length classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GCPntsAbscissaType {
    AbsComposite,
    LengthParametrized,
    Parametrized,
}

/// OCCT math_GaussSingleIntegration::Perform
/// (math_GaussSingleIntegration.cxx L100-155) — the Gauss quadrature of
/// `value` over [lower, upper] with `order` points; returns None for the
/// OCCT `Done = false` path (a failed F.Value callback cannot happen for the
/// GCPnts integrands, so None is only the initialization state).
fn gauss_single_integration_perform(
    value: &mut dyn FnMut(f64) -> f64,
    lower: f64,
    upper: f64,
    order: usize,
) -> Option<f64> {
    let the_order = order.min(gauss_points_max());
    let mut gauss_p = VecD::new(the_order + 1);
    let mut gauss_w = VecD::new(the_order + 1);
    // math::GaussPoints(Order, GaussP); math::GaussWeights(Order, GaussW).
    gauss_points(the_order, &mut gauss_p);
    gauss_weights(the_order, &mut gauss_w);

    // Changement de variable pour la mise a l'echelle [Lower, Upper].
    let xm = 0.5 * (upper + lower);
    let xr = 0.5 * (upper - lower);
    let mut val = 0.0f64;

    let ind = the_order / 2;
    let ind1 = (the_order + 1) / 2;
    if ind1 > ind {
        // odder case
        val = value(xm);
        val *= gauss_w.get(ind1);
    }
    // Sommation sur tous les points de Gauss: avec utilisation de la symetrie.
    for j in 1..=ind {
        let dx = xr * gauss_p.get(j);
        let f1 = value(xm - dx);
        let f2 = value(xm + dx);
        // Multiplication par les poids de Gauss.
        let ft = f1 + f2;
        val += gauss_w.get(j) * ft;
    }
    // Mise a l'echelle de l'intervalle [Lower, Upper]
    val *= xr;
    Some(val)
}

/// OCCT math_GaussSingleIntegration(F, Lower, Upper, Order, Tol)
/// (math_GaussSingleIntegration.cxx L64-98) — the adaptive subdivision loop.
fn gauss_single_integration_tol(
    value: &mut dyn FnMut(f64) -> f64,
    lower: f64,
    upper: f64,
    order: usize,
    tol: f64,
) -> Option<f64> {
    let the_order = order.min(gauss_points_max());
    let iter_max = 13; // Max number of iteration
    let mut n_iter = 1; // current number of iteration
    let mut nb_interval = 1; // current number of subintervals

    let mut len = match gauss_single_integration_perform(value, lower, upper, the_order) {
        Some(v) => v,
        None => return None,
    };
    loop {
        let old_len = len;
        len = 0.0;
        nb_interval *= 2;
        let du = (upper - lower) / nb_interval as f64;
        for i in 1..=nb_interval {
            let v = gauss_single_integration_perform(
                value,
                lower + (i - 1) as f64 * du,
                lower + i as f64 * du,
                the_order,
            )?;
            len += v;
        }
        n_iter += 1;
        if !((old_len - len).abs() > tol && n_iter <= iter_max) {
            break;
        }
    }
    Some(len)
}

/// OCCT CPnts_AbscissaPoint.cxx L41-47 — the 3D integrand |C'(X)|.
fn f3d(c: &dyn GCPntsCurve, x: f64) -> f64 {
    let (_p, v): (DVec3, DVec3) = c.d1(x);
    v.length()
}

/// OCCT CPnts_AbscissaPoint.cxx L49-55 — the 2D integrand |C'(X)| (lifted
/// to the 3D slot by the trait shim, the magnitude is unaffected).
fn f2d(c: &dyn GCPntsCurve, x: f64) -> f64 {
    let (_p, v): (DVec3, DVec3) = c.d1(x);
    v.length()
}

/// OCCT CPnts_AbscissaPoint.cxx L57-77 — the 3D Gauss order.
fn order_3d(c: &dyn GCPntsCurve) -> i32 {
    use rcad_kernel::base::proj_lib::CurveType;
    match c.get_type() {
        CurveType::Line => 2,
        CurveType::Parabola => 5,
        CurveType::Bezier => (2 * c.curve_degree()).min(24),
        CurveType::BSpline => (2 * c.nb_poles() - 1).min(24),
        _ => 10,
    }
}

/// OCCT CPnts_AbscissaPoint.cxx L79-99 — the 2D Gauss order.
fn order_2d(c: &dyn GCPntsCurve) -> i32 {
    use rcad_kernel::base::proj_lib::CurveType;
    match c.get_type() {
        CurveType::Line => 2,
        CurveType::Parabola => 5,
        CurveType::Bezier => (2 * c.curve_degree()).min(24),
        CurveType::BSpline => (2 * c.nb_poles() - 1).min(24),
        _ => 10,
    }
}

/// OCCT computeType (GCPnts_AbscissaPoint.cxx L26-65) — computes the type
/// and the length ratio if GCPnts_LengthParametrized.
fn compute_type(c: &dyn GCPntsCurve, the_ratio: &mut f64) -> GCPntsAbscissaType {
    use rcad_kernel::base::proj_lib::CurveType;
    if c.nb_intervals_cn() > 1 {
        return GCPntsAbscissaType::AbsComposite;
    }

    match c.get_type() {
        CurveType::Line => {
            *the_ratio = 1.0;
            GCPntsAbscissaType::LengthParametrized
        }
        CurveType::Circle => {
            *the_ratio = c.circle_radius();
            GCPntsAbscissaType::LengthParametrized
        }
        CurveType::Bezier => {
            // Handle(BezierCurve) aBz = theC.Bezier();
            // if (aBz->NbPoles() == 2 && !aBz->IsRational()).
            if c.nb_poles() == 2 && !c.is_rational() {
                *the_ratio = c.dn1(0.0).length();
                GCPntsAbscissaType::LengthParametrized
            } else {
                GCPntsAbscissaType::Parametrized
            }
        }
        CurveType::BSpline => {
            // Handle(BSplineCurve) aBs = theC.BSpline();
            // if (aBs->NbPoles() == 2 && !aBs->IsRational()).
            if c.nb_poles() == 2 && !c.is_rational() {
                *the_ratio = c.dn1(c.first_parameter()).length();
                GCPntsAbscissaType::LengthParametrized
            } else {
                GCPntsAbscissaType::Parametrized
            }
        }
        _ => GCPntsAbscissaType::Parametrized,
    }
}

/// OCCT CPnts_AbscissaPoint::Length(C, U1, U2) — the `f3d` flavor
/// (CPnts_AbscissaPoint.cxx L131-144).
pub fn cpnts_length_3d(c: &dyn GCPntsCurve, u1: f64, u2: f64) -> f64 {
    let mut integrand = |x: f64| f3d(c, x);
    let the_length = gauss_single_integration_perform(&mut integrand, u1, u2, order_3d(c) as usize);
    match the_length {
        Some(v) => v.abs(),
        // throw Standard_ConstructionError().
        None => panic!("Standard_ConstructionError"),
    }
}

/// OCCT CPnts_AbscissaPoint::Length(C, U1, U2) — the `f2d` flavor
/// (CPnts_AbscissaPoint.cxx L148-161).
pub fn cpnts_length_2d(c: &dyn GCPntsCurve, u1: f64, u2: f64) -> f64 {
    let mut integrand = |x: f64| f2d(c, x);
    let the_length = gauss_single_integration_perform(&mut integrand, u1, u2, order_2d(c) as usize);
    match the_length {
        Some(v) => v.abs(),
        // throw Standard_ConstructionError().
        None => panic!("Standard_ConstructionError"),
    }
}

/// OCCT CPnts_AbscissaPoint::Length(C, U1, U2, Tol) — the `f3d` flavor
/// (CPnts_AbscissaPoint.cxx L168-184).
pub fn cpnts_length_3d_tol(c: &dyn GCPntsCurve, u1: f64, u2: f64, tol: f64) -> f64 {
    let mut integrand = |x: f64| f3d(c, x);
    let the_length =
        gauss_single_integration_tol(&mut integrand, u1, u2, order_3d(c) as usize, tol);
    match the_length {
        Some(v) => v.abs(),
        // throw Standard_ConstructionError().
        None => panic!("Standard_ConstructionError"),
    }
}

/// OCCT CPnts_AbscissaPoint::Length(C, U1, U2, Tol) — the `f2d` flavor
/// (CPnts_AbscissaPoint.cxx L191-207).
pub fn cpnts_length_2d_tol(c: &dyn GCPntsCurve, u1: f64, u2: f64, tol: f64) -> f64 {
    let mut integrand = |x: f64| f2d(c, x);
    let the_length =
        gauss_single_integration_tol(&mut integrand, u1, u2, order_2d(c) as usize, tol);
    match the_length {
        Some(v) => v.abs(),
        // throw Standard_ConstructionError().
        None => panic!("Standard_ConstructionError"),
    }
}

/// OCCT GCPnts_AbscissaPoint::length — the private template instantiated
/// for `Adaptor3d_Curve` (GCPnts_AbscissaPoint.cxx L372-423).
fn length_3d_over(c: &dyn GCPntsCurve, the_u1: f64, the_u2: f64, the_tol: Option<f64>) -> f64 {
    let mut a_ratio = 1.0f64;
    let a_type = compute_type(c, &mut a_ratio);
    match a_type {
        GCPntsAbscissaType::LengthParametrized => (the_u2 - the_u1).abs() * a_ratio,
        GCPntsAbscissaType::Parametrized => match the_tol {
            Some(tol) => cpnts_length_3d_tol(c, the_u1, the_u2, tol),
            None => cpnts_length_3d(c, the_u1, the_u2),
        },
        GCPntsAbscissaType::AbsComposite => {
            let a_nb_intervals = c.nb_intervals_cn();
            // NCollection_Array1<double> aTI(1, aNbIntervals + 1);
            // theC.Intervals(aTI, GeomAbs_CN);
            let a_ti = c.intervals_cn();
            let a_uu1 = the_u1.min(the_u2);
            let a_uu2 = the_u1.max(the_u2);
            let mut a_l = 0.0f64;
            for an_index in 1..=a_nb_intervals {
                let lo = a_ti[an_index - 1];
                let hi = a_ti[an_index];
                if lo > a_uu2 {
                    break;
                }
                if hi < a_uu1 {
                    continue;
                }
                match the_tol {
                    Some(tol) => {
                        a_l += cpnts_length_3d_tol(c, lo.max(a_uu1), hi.min(a_uu2), tol);
                    }
                    None => {
                        a_l += cpnts_length_3d(c, lo.max(a_uu1), hi.min(a_uu2));
                    }
                }
            }
            a_l
        }
    }
}

/// OCCT GCPnts_AbscissaPoint::length — the private template instantiated
/// for `Adaptor2d_Curve2d` (GCPnts_AbscissaPoint.cxx L372-423).
fn length_2d_over(c: &dyn GCPntsCurve, the_u1: f64, the_u2: f64, the_tol: Option<f64>) -> f64 {
    let mut a_ratio = 1.0f64;
    let a_type = compute_type(c, &mut a_ratio);
    match a_type {
        GCPntsAbscissaType::LengthParametrized => (the_u2 - the_u1).abs() * a_ratio,
        GCPntsAbscissaType::Parametrized => match the_tol {
            Some(tol) => cpnts_length_2d_tol(c, the_u1, the_u2, tol),
            None => cpnts_length_2d(c, the_u1, the_u2),
        },
        GCPntsAbscissaType::AbsComposite => {
            let a_nb_intervals = c.nb_intervals_cn();
            let a_ti = c.intervals_cn();
            let a_uu1 = the_u1.min(the_u2);
            let a_uu2 = the_u1.max(the_u2);
            let mut a_l = 0.0f64;
            for an_index in 1..=a_nb_intervals {
                let lo = a_ti[an_index - 1];
                let hi = a_ti[an_index];
                if lo > a_uu2 {
                    break;
                }
                if hi < a_uu1 {
                    continue;
                }
                match the_tol {
                    Some(tol) => {
                        a_l += cpnts_length_2d_tol(c, lo.max(a_uu1), hi.min(a_uu2), tol);
                    }
                    None => {
                        a_l += cpnts_length_2d(c, lo.max(a_uu1), hi.min(a_uu2));
                    }
                }
            }
            a_l
        }
    }
}

/// OCCT GCPnts_AbscissaPoint::Length(const Adaptor3d_Curve& theC)
/// (GCPnts_AbscissaPoint.cxx L305-308).
pub fn gcpnts_length_3d(
    c: &dyn GCPntsCurve,
    ) -> f64 {
    length_3d_over(c, c.first_parameter(), c.last_parameter(), None)
}

/// OCCT GCPnts_AbscissaPoint::Length(const Adaptor2d_Curve2d& theC)
/// (GCPnts_AbscissaPoint.cxx L312-315).
pub fn gcpnts_length_2d(
    c: &dyn GCPntsCurve,
    ) -> f64 {
    length_2d_over(c, c.first_parameter(), c.last_parameter(), None)
}

/// OCCT GCPnts_AbscissaPoint::Length(const Adaptor3d_Curve& theC,
/// theU1, theU2) (GCPnts_AbscissaPoint.cxx L333-338).
pub fn gcpnts_length_3d_range(c: &dyn GCPntsCurve, the_u1: f64, the_u2: f64) -> f64 {
    length_3d_over(c, the_u1, the_u2, None)
}

/// OCCT GCPnts_AbscissaPoint::Length(const Adaptor2d_Curve2d& theC,
/// theU1, theU2) (GCPnts_AbscissaPoint.cxx L342-347).  Used by the
/// Discretise ACR law construction
/// (GeomPlate_BuildPlateSurface.cxx L2207/L2219).
pub fn gcpnts_length_2d_range(c: &dyn GCPntsCurve, the_u1: f64, the_u2: f64) -> f64 {
    length_2d_over(c, the_u1, the_u2, None)
}

/// OCCT GCPnts_AbscissaPoint::Length(C, U1, U2, Tol) — Adaptor3d_Curve
/// flavor (GCPnts_AbscissaPoint.cxx L351-357).
pub fn gcpnts_length_3d_range_tol(
    c: &dyn GCPntsCurve,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> f64 {
    length_3d_over(c, the_u1, the_u2, Some(the_tol))
}

/// OCCT GCPnts_AbscissaPoint::Length(C, U1, U2, Tol) — Adaptor2d_Curve2d
/// flavor (GCPnts_AbscissaPoint.cxx L361-367).
pub fn gcpnts_length_2d_range_tol(
    c: &dyn GCPntsCurve,
    the_u1: f64,
    the_u2: f64,
    the_tol: f64,
) -> f64 {
    length_2d_over(c, the_u1, the_u2, Some(the_tol))
}
