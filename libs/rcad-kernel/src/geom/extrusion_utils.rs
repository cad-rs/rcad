//! OCCT Geom_ExtrusionUtils (TKG3d/Geom/Geom_ExtrusionUtils.pxx) — the
//! extrusion-surface evaluation helpers, plus the rcad re-hosts of
//! `Geom_SurfaceOfLinearExtrusion::EvalD1/EvalD2/EvalDN`
//! (Geom_SurfaceOfLinearExtrusion.cxx L166-184 / L188-210 / L245-271) which
//! consume them.
//!
//! Scope: the bodies needed by the `Geom_OffsetSurface::UIso` extrusion arm
//! (Geom_OffsetSurface.cxx L609-623 — `basisSurf->D1(UU, 0., aP, aD1U, aD1V)`)
//! and by the `Surface3::dn` dispatch (`GeomAdaptor_Surface::EvalDN` L1754:
//! `Geom_ExtrusionUtils::DN`), i.e. `CalculateD0/D1/D2/DN` (pxx L40-71 /
//! L81-99 / L154-168) with the `D0` / `D1` / `D2` / `DN` template wrappers
//! (pxx L179-207 / L228-252 / L318-332).  The `CalculateD3` / `D3` bodies
//! (pxx L119-146 / L272-306) are NOT translated here — rcad has no consumer of
//! `Geom_SurfaceOfLinearExtrusion::EvalD3`.
//!
//! Extrusion surface: `P(U,V) = C(U) + V * Direction`.
//!
//! Architecture differences: `occ::handle(Geom_Curve) basisCurve` maps to the
//! rcad [`Curve3`] value, the `gp_Ax2 loc` axis of the OCCT class does not
//! exist in the rcad [`LinearExtrusionSurface`] payload (the pxx bodies do not
//! consume it either), and `theBasis.EvalDN(U, N)` (the `Geom_Curve` virtual)
//! is the rcad [`crate::geom::curve_dn::curve_dn`] union.

use glam::DVec3;

use super::curve_dn::curve_dn;
use super::offset_surface_utils::{ResD1, ResD2};
use crate::geom::{Curve3, CurveEval, LinearExtrusionSurface, Surface3};

/// OCCT `Geom_ExtrusionUtils::CalculateD0` (pxx L40-50).
pub fn calculate_d0(the_curve_pt: DVec3, the_v: f64, the_dir: DVec3) -> DVec3 {
    the_curve_pt + the_v * the_dir
}

/// OCCT `Geom_ExtrusionUtils::CalculateD1` (pxx L56-71).
pub fn calculate_d1(
    the_curve_pt: DVec3,
    the_curve_d1: DVec3,
    the_v: f64,
    the_dir: DVec3,
) -> ResD1 {
    ResD1 {
        point: the_curve_pt + the_v * the_dir,
        d1u: the_curve_d1,
        d1v: the_dir,
    }
}

/// OCCT `Geom_ExtrusionUtils::CalculateD2` (pxx L81-99).
pub fn calculate_d2(
    the_curve_pt: DVec3,
    the_curve_d1: DVec3,
    the_curve_d2: DVec3,
    the_v: f64,
    the_dir: DVec3,
) -> ResD2 {
    ResD2 {
        point: the_curve_pt + the_v * the_dir,
        d1u: the_curve_d1,
        d1v: the_dir,
        d2u: the_curve_d2,
        d2v: DVec3::ZERO,
        d2uv: DVec3::ZERO,
    }
}

/// OCCT `Geom_ExtrusionUtils::CalculateDN` (pxx L154-168).
pub fn calculate_dn(the_curve_dn: DVec3, the_dir: DVec3, the_der_u: i32, the_der_v: i32) -> DVec3 {
    if the_der_v == 0 {
        the_curve_dn
    } else if the_der_u == 0 && the_der_v == 1 {
        the_dir
    } else {
        DVec3::ZERO
    }
}

/// OCCT `Geom_ExtrusionUtils::D0` (pxx L179-199) — the basis-curve
/// `EvalD0(U)` + `CalculateD0`.
pub fn d0(the_u: f64, the_v: f64, the_basis: &Curve3, the_dir: DVec3) -> DVec3 {
    calculate_d0(the_basis.point_at(the_u), the_v, the_dir)
}

/// OCCT `Geom_ExtrusionUtils::D1` (pxx L201-220) — the basis-curve
/// `EvalD1(U)` + `CalculateD1`.
pub fn d1(the_u: f64, the_v: f64, the_basis: &Curve3, the_dir: DVec3) -> ResD1 {
    calculate_d1(
        the_basis.point_at(the_u),
        the_basis.derivative_at(the_u),
        the_v,
        the_dir,
    )
}

/// OCCT `Geom_ExtrusionUtils::D2` (pxx L228-252) — the basis-curve
/// `EvalD2(U)` + `CalculateD2`.
pub fn d2(the_u: f64, the_v: f64, the_basis: &Curve3, the_dir: DVec3) -> ResD2 {
    calculate_d2(
        the_basis.point_at(the_u),
        the_basis.derivative_at(the_u),
        the_basis.derivative2_at(the_u),
        the_v,
        the_dir,
    )
}

/// OCCT `Geom_ExtrusionUtils::DN` (pxx L318-332) — the basis-curve
/// `EvalDN(U, theDerU)` + `CalculateDN` (the curve derivative is only
/// requested for `theDerV == 0`, as in the OCCT body).
pub fn dn(
    the_u: f64,
    the_basis: &Curve3,
    the_dir: DVec3,
    the_der_u: i32,
    the_der_v: i32,
) -> DVec3 {
    let mut a_curve_dn = DVec3::ZERO;
    if the_der_v == 0 {
        a_curve_dn = curve_dn(the_basis, the_u, the_der_u);
    }
    calculate_dn(a_curve_dn, the_dir, the_der_u, the_der_v)
}

/// OCCT `Geom_SurfaceOfLinearExtrusion::EvalD0` (cxx L150-163) — the
/// `GeomEval_RepUtils::TryEvalSurfaceD0(myEvalRep, ...)` short circuit is the
/// rcad `Surface3` payload's own concern, not reproduced here.
pub fn linear_extrusion_eval_d0(le: &LinearExtrusionSurface, u: f64, v: f64) -> DVec3 {
    d0(u, v, &le.profile, le.direction)
}

/// OCCT `Geom_SurfaceOfLinearExtrusion::EvalD1` (cxx L166-184).
pub fn linear_extrusion_eval_d1(le: &LinearExtrusionSurface, u: f64, v: f64) -> ResD1 {
    d1(u, v, &le.profile, le.direction)
}

/// OCCT `Geom_SurfaceOfLinearExtrusion::EvalD2` (cxx L188-210).
pub fn linear_extrusion_eval_d2(le: &LinearExtrusionSurface, u: f64, v: f64) -> ResD2 {
    d2(u, v, &le.profile, le.direction)
}

/// OCCT `Geom_SurfaceOfLinearExtrusion::EvalDN` (cxx L245-271).
///
/// Architecture difference: the OCCT `V` parameter only feeds the
/// `GeomEval_RepUtils::TryEvalSurfaceDN(myEvalRep, U, V, Nu, Nv)` short
/// circuit (the occluded-surface representation of a translated basis), which
/// is the rcad payload's own concern — the `CalculateDN` body itself is
/// V-independent (pxx L154-168).
pub fn linear_extrusion_eval_dn(
    le: &LinearExtrusionSurface,
    u: f64,
    v: f64,
    nu: i32,
    nv: i32,
) -> DVec3 {
    let _ = v;
    // OCCT: if (Nu + Nv < 1 || Nu < 0 || Nv < 0) throw Geom_UndefinedDerivative.
    assert!(
        nu + nv >= 1 && nu >= 0 && nv >= 0,
        "Geom_UndefinedDerivative: Geom_SurfaceOfLinearExtrusion::EvalDN"
    );

    if nv == 0 {
        let a_dn = curve_dn(&le.profile, u, nu);
        return calculate_dn(a_dn, le.direction, nu, nv);
    } else if nu == 0 && nv == 1 {
        return le.direction;
    }
    DVec3::ZERO
}

/// OCCT `Geom_Surface::EvalD1` for the rcad extrusion surface value — the
/// arm of the offset UIso extrusion branch (`aGAsurf.GetType() ==
/// GeomAbs_SurfaceOfExtrusion`).
pub fn surface_eval_d1(the_s: &Surface3, u: f64, v: f64) -> ResD1 {
    match the_s {
        Surface3::LinearExtrusion(le) => linear_extrusion_eval_d1(le, u, v),
        _ => panic!(
            "GAP: Geom_Surface::EvalD1 (TKG3d/Geom) is not translated for this surface type \
             (rcad carries the Geom_SurfaceOfLinearExtrusion arm only) — \
             Geom_OffsetSurface::UIso extrusion branch"
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::{Circle3, Curve3};

    /// A circle of radius `R` in the XY plane extruded along `+Z`: the surface
    /// is the cylinder `P(u, v) = (R cos u, R sin u, v)`, whose exact
    /// derivatives are
    ///   `Pu = R(-sin u, cos u, 0)`, `Pv = (0, 0, 1)`,
    ///   `Puu = -R(cos u, sin u, 0)`, `Puuu = R(sin u, -cos u, 0)`,
    /// with every derivative of order >= 2 in V equal to zero.
    fn extruded_circle() -> LinearExtrusionSurface {
        LinearExtrusionSurface {
            profile: Box::new(Curve3::Circle(Circle3 {
                center: DVec3::ZERO,
                normal: DVec3::Z,
                radius: 1.5,
                x_dir: DVec3::X,
                y_dir: DVec3::Y,
            })),
            direction: DVec3::Z,
        }
    }

    #[test]
    fn extrusion_eval_d1_d2_dn_match_the_closed_form() {
        let le = extruded_circle();
        let r = 1.5;
        let (u, v) = (0.7f64, 1.3f64);
        let (su, cu) = u.sin_cos();
        let want_p = DVec3::new(r * cu, r * su, v);
        let want_d1u = DVec3::new(-r * su, r * cu, 0.0);
        let want_d1v = DVec3::Z;
        let want_d2u = DVec3::new(-r * cu, -r * su, 0.0);
        let want_d3u = DVec3::new(r * su, -r * cu, 0.0);
        let want_d4u = DVec3::new(r * cu, r * su, 0.0);

        let d1 = linear_extrusion_eval_d1(&le, u, v);
        assert!((d1.point - want_p).length() < 1e-15, "p={:?}", d1.point);
        assert!((d1.d1u - want_d1u).length() < 1e-15, "d1u={:?}", d1.d1u);
        assert!((d1.d1v - want_d1v).length() < 1e-15);
        let d2 = linear_extrusion_eval_d2(&le, u, v);
        assert!((d2.point - want_p).length() < 1e-15);
        assert!((d2.d1u - want_d1u).length() < 1e-15);
        assert!((d2.d1v - want_d1v).length() < 1e-15);
        assert!((d2.d2u - want_d2u).length() < 1e-15, "d2u={:?}", d2.d2u);
        assert!(d2.d2v.length() < 1e-15);
        assert!(d2.d2uv.length() < 1e-15);
        assert!((linear_extrusion_eval_dn(&le, u, v, 1, 0) - want_d1u).length() < 1e-15);
        assert!((linear_extrusion_eval_dn(&le, u, v, 0, 1) - want_d1v).length() < 1e-15);
        assert!((linear_extrusion_eval_dn(&le, u, v, 2, 0) - want_d2u).length() < 1e-14);
        assert!((linear_extrusion_eval_dn(&le, u, v, 3, 0) - want_d3u).length() < 1e-14);
        assert!((linear_extrusion_eval_dn(&le, u, v, 4, 0) - want_d4u).length() < 1e-13);
        assert!(linear_extrusion_eval_dn(&le, u, v, 0, 2).length() < 1e-15);
        assert!(linear_extrusion_eval_dn(&le, u, v, 1, 1).length() < 1e-15);
    }
}
