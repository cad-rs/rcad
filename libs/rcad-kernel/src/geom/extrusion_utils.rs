//! OCCT Geom_ExtrusionUtils (TKG3d/Geom/Geom_ExtrusionUtils.pxx) — the
//! extrusion-surface evaluation helpers, plus the rcad re-host of
//! `Geom_SurfaceOfLinearExtrusion::EvalD1` (Geom_SurfaceOfLinearExtrusion.cxx
//! L166-184) which consumes them.
//!
//! Scope: the bodies needed by the `Geom_OffsetSurface::UIso` extrusion arm
//! (Geom_OffsetSurface.cxx L609-623 — `basisSurf->D1(UU, 0., aP, aD1U, aD1V)`),
//! i.e. `Geom_ExtrusionUtils::CalculateD0` (pxx L40-50) and `CalculateD1`
//! (pxx L56-71) with the `D0` / `D1` template wrappers (pxx L179-207).  The
//! `CalculateD2/D3/DN` bodies (pxx L81-340) belong to
//! `Geom_SurfaceOfLinearExtrusion::EvalD2/EvalD3/DN` and are NOT translated
//! here — no consumer in rcad today.
//!
//! Extrusion surface: `P(U,V) = C(U) + V * Direction`.
//!
//! Architecture differences: `occ::handle(Geom_Curve) basisCurve` maps to the
//! rcad [`Curve3`] value, and the `gp_Ax2 loc` axis of the OCCT class does not
//! exist in the rcad [`LinearExtrusionSurface`] payload (the pxx bodies do not
//! consume it either).

use glam::DVec3;

use super::offset_surface_utils::ResD1;
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
