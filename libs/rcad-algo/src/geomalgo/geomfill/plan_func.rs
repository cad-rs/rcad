//! OCCT GeomFill_PlanFunc (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_PlanFunc.hxx (members) + GeomFill_PlanFunc.cxx (whole file
//! L22-91).  The class derives `math_FunctionWithDerivative`; the rcad form
//! implements the 1-variable/1-equation [`FunctionSetWithDerivatives`]
//! (the OCCT internal `math_MyFunctionSetWithDerivatives` wrapper protocol).

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

/// OCCT GeomFill_PlanFunc (GeomFill_PlanFunc.hxx): myPnt, myVec, G, V
/// (gp_XYZ) + the myCurve handle.
#[derive(Debug, Clone)]
pub struct PlanFunc {
    /// OCCT gp_XYZ myPnt.
    my_pnt: DVec3,
    /// OCCT gp_XYZ myVec.
    my_vec: DVec3,
    /// OCCT gp_XYZ G (scratch).
    g: DVec3,
    /// OCCT gp_XYZ V (scratch).
    v: DVec3,
    /// OCCT handle(Adaptor3d_Curve) myCurve.
    my_curve: Option<Curve3>,
}

impl PlanFunc {
    /// OCCT GeomFill_PlanFunc::GeomFill_PlanFunc (L22-29).
    pub fn new(the_p: DVec3, the_v: DVec3, the_c: &Curve3) -> Self {
        PlanFunc {
            my_pnt: the_p,
            my_vec: the_v,
            g: DVec3::ZERO,
            v: DVec3::ZERO,
            my_curve: Some(the_c.clone()),
        }
    }

    /// OCCT GeomFill_PlanFunc::Value (L31-37).
    pub fn value(&mut self, x: f64, f: &mut f64) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        self.g = my_curve.point_at(x);
        // OCCT: V.SetLinearForm(-1, myPnt, G.XYZ()).
        self.v = self.g - self.my_pnt;
        *f = self.my_vec.dot(self.v);
        true
    }

    /// OCCT GeomFill_PlanFunc::Derivative (L39-46).
    pub fn derivative(&mut self, x: f64, d: &mut f64) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        self.g = my_curve.point_at(x);
        let dg = my_curve.derivative_at(x);
        *d = self.my_vec.dot(dg);
        true
    }

    /// OCCT GeomFill_PlanFunc::Values (L48-59).
    pub fn values(&mut self, x: f64, f: &mut f64, d: &mut f64) -> bool {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        self.g = my_curve.point_at(x);
        let dg = my_curve.derivative_at(x);
        self.v = self.g - self.my_pnt;
        *f = self.my_vec.dot(self.v);
        *d = self.my_vec.dot(dg);
        true
    }

    /// OCCT GeomFill_PlanFunc::D2 (L61-65) — the literal OCCT empty body.
    pub fn d2(&mut self, _x: f64, _f: &mut f64, _d1: &mut f64, _d2: &mut f64) {}

    /// OCCT GeomFill_PlanFunc::DEDT (L67-73).
    pub fn dedt(&mut self, x: f64, d_pnt: DVec3, d_vec: DVec3, dfdt: &mut f64) {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        self.g = my_curve.point_at(x);
        self.v = self.g - self.my_pnt;
        *dfdt = d_vec.dot(self.v) - self.my_vec.dot(d_pnt);
    }

    /// OCCT GeomFill_PlanFunc::D2E (L75-91).
    #[allow(clippy::too_many_arguments)]
    pub fn d2e(&mut self, x: f64, dp: DVec3, dv: DVec3, dfdt: &mut f64) {
        let my_curve = self.my_curve.as_ref().expect("null myCurve");
        self.g = my_curve.point_at(x);
        let dg = my_curve.derivative_at(x);
        self.v = self.g - self.my_pnt;
        // OCCT: DVDT.SetLinearForm(-1, DP.XYZ(), G.XYZ()).
        let dvdt = self.g - dp;
        *dfdt = dv.dot(self.v) + self.my_vec.dot(dvdt);
    }
}

impl FunctionSetWithDerivatives for PlanFunc {
    /// OCCT math_MyFunctionSetWithDerivatives::NbVariables — 1.
    fn nb_variables(&self) -> usize {
        1
    }

    /// OCCT math_MyFunctionSetWithDerivatives::NbEquations — 1.
    fn nb_equations(&self) -> usize {
        1
    }

    /// OCCT math_MyFunctionSetWithDerivatives::Value.
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let mut fv = 0.0;
        let ok = self.value(x[0], &mut fv);
        f[0] = fv;
        ok
    }

    /// OCCT math_MyFunctionSetWithDerivatives::Derivatives.
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        let mut dv = 0.0;
        let ok = self.derivative(x[0], &mut dv);
        df[0][0] = dv;
        ok
    }

    /// OCCT math_MyFunctionSetWithDerivatives::Values.
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        let (mut fv, mut dv) = (0.0, 0.0);
        let ok = self.values(x[0], &mut fv, &mut dv);
        f[0] = fv;
        df[0][0] = dv;
        ok
    }
}
