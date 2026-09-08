//! GAP carrier: OCCT Approx_CurvlinFunc (TKGeomAlgo/Approx) — the
//! curvilinear-abscissa function of a curve.
//!
//! This class is NOT part of the GeomFill_Pipe + GeomFill_Sweep closure
//! batch (plan §9 E0 R4 staged); rcad does not host a translation yet, so
//! the carrier preserves the OCCT failure path (panic at first use) until
//! the Approx package batch lands.  Consumers in this module:
//! GeomFill_GuideTrihedronAC / GeomFill_GuideTrihedronPlan-era guide laws
//! and GeomFill_LocationGuide (GeomFill_LocationGuide.cxx L103+).
//!
//! OCCT anchor: Approx_CurvlinFunc.hxx + Approx_CurvlinFunc.cxx
//! (GetLength / GetSParameter / GetUParameter / NbIntervals / Intervals /
//! Trim).

use rcad_kernel::geom::Curve3;
use rcad_kernel::math::GeomAbsShape;

/// OCCT handle(Approx_CurvlinFunc) — GAP carrier.
#[derive(Debug, Clone)]
pub struct ApproxCurvlinFunc;

impl ApproxCurvlinFunc {
    /// OCCT Approx_CurvlinFunc::Approx_CurvlinFunc(Curve, Tol).
    pub fn new(_curve: &Curve3, _tol: f64) -> Self {
        panic!(
            "GAP: Approx_CurvlinFunc (TKGeomAlgo/Approx) is not translated — see file header"
        )
    }

    /// OCCT Approx_CurvlinFunc::GetLength().
    pub fn get_length(&self) -> f64 {
        panic!(
            "GAP: Approx_CurvlinFunc (TKGeomAlgo/Approx) is not translated — see file header"
        )
    }

    /// OCCT Approx_CurvlinFunc::GetSParameter(U).
    pub fn get_s_parameter(&self, _u: f64) -> f64 {
        panic!(
            "GAP: Approx_CurvlinFunc (TKGeomAlgo/Approx) is not translated — see file header"
        )
    }

    /// OCCT Approx_CurvlinFunc::GetUParameter(C, S, Order).
    pub fn get_u_parameter(&self, _c: &Curve3, _s: f64, _order: i32) -> f64 {
        panic!(
            "GAP: Approx_CurvlinFunc (TKGeomAlgo/Approx) is not translated — see file header"
        )
    }

    /// OCCT Approx_CurvlinFunc::NbIntervals(S).
    pub fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        panic!(
            "GAP: Approx_CurvlinFunc (TKGeomAlgo/Approx) is not translated — see file header"
        )
    }

    /// OCCT Approx_CurvlinFunc::Intervals(T, S).
    pub fn intervals(&self, _t: &mut Vec<f64>, _s: GeomAbsShape) {
        panic!(
            "GAP: Approx_CurvlinFunc (TKGeomAlgo/Approx) is not translated — see file header"
        )
    }
}
