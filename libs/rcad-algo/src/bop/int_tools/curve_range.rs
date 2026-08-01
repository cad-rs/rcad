// OCCT Adaptor3d_Curve::Resolution — curve parameter step helpers.
// Ported from the legacy rcad-algorithms implementation (curve_range.rs).
use rcad_kernel::geom::{Curve3, CurveEval};

/// Curve parameter step: parameter increment needed to move `tol` distance along curve.
///
/// OCCT: Adaptor3d_Curve::Resolution(tol) = tol / |dP/dt|
/// (BRepLib_1.cxx L61, IntTools_ShrunkRange.cxx L162)
///
/// When the tangent speed is nearly zero (singularity), the resolution is clamped to `tol`.
pub fn curve_resolution(curve: &Curve3, t: f64, tol: f64) -> f64 {
    let speed = curve.tangent_at(t).length();
    if speed < 1e-15 {
        tol
    } else {
        tol / speed
    }
}
