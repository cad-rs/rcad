//! OCCT Extrema_Curve2dTool (TKGeomBase/Extrema/Extrema_Curve2dTool.hxx
//! L38-157) — the static facade over `Adaptor2d_Curve2d` used by the 2D half
//! of the Extrema package (`Extrema_ExtCC2d`, `Extrema_ECC2d`,
//! `Extrema_GGenExtCC<..., Extrema_Curve2dTool, ...>`).
//!
//! rcad encoding (architecture difference, annotated): in OCCT every
//! `Extrema_Curve2dTool` member is a one-liner forwarding to the adaptor
//! (`FirstParameter(C)` -> `C.FirstParameter()`, `Value(C, U)` ->
//! `C.Value(U)`, ...).  The rcad `base::proj_lib::adaptor::Adaptor2dCurve2d`
//! trait carries those members verbatim (evaluation, interval machinery,
//! periodicity, resolution and the analytic downcasts), so the facade
//! collapses onto the trait object type: the Extrema 2D algorithms are written
//! against [`ExtremaCurve2dTool`] exactly as the OCCT templates are written
//! against `TheCurveTool2 = Extrema_Curve2dTool`, and every adaptor
//! implementing `Adaptor2dCurve2d` (e.g. `Geom2dAdaptor_Curve`) is a tool.
//! The same encoding is used by `extrema_glob_opt_func_cc`'s
//! `GlobOptCurves::Curves2d` arm.
//!
//! The one member the rcad adaptor trait does not carry is the arc length
//! `Extrema_GGenExtCC::Perform` requests through
//! `GCPnts_AbscissaPoint::Length(C)` (Extrema_GGenExtCC.hxx L476-491): the
//! kernel `base::gcpnts::abscissa_point` engine is `Curve3`-based, so the
//! facade answers the `-1` sentinel — the state OCCT itself reaches in the
//! same place whenever a curve parameter is infinite, which keeps
//! `indmax == -1` and skips the interval-count optimization.

use crate::base::proj_lib::adaptor::Adaptor2dCurve2d;

/// OCCT Extrema_Curve2dTool (Extrema_Curve2dTool.hxx L38-157) — the static
/// facade over `Adaptor2d_Curve2d`.
pub trait ExtremaCurve2dTool: Adaptor2dCurve2d {
    /// OCCT `GCPnts_AbscissaPoint::Length(C)` (CPnts_AbscissaPoint.cxx
    /// L103-106) — the arc length of the adaptor over its whole domain, as
    /// `Extrema_GGenExtCC::Perform` consumes it (hxx L476-491).
    ///
    /// rcad encoding: the arc-length engine is `Curve3`-based and there is no
    /// 2D engine yet, so the facade answers `-1`; `Perform` then keeps
    /// `indmax == -1`, exactly the OCCT state when a curve parameter is
    /// infinite (the `Length` result is only consumed when all four boundary
    /// parameters are finite).
    fn abscissa_length(&self) -> f64 {
        -1.0
    }
}

impl<T: Adaptor2dCurve2d + ?Sized> ExtremaCurve2dTool for T {}
