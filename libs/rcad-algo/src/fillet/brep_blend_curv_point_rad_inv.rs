//! OCCT BRepBlend_CurvPointRadInv (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_CurvPointRadInv.hxx (L41-95) and BRepBlend_CurvPointRadInv.cxx
//! (whole file L24-130).
//!
//! Architecture mappings: `class BRepBlend_CurvPointRadInv :
//! public Blend_CurvPointFuncInv` is expressed by implementing the
//! [`BlendCurvPointFuncInv`] trait and the `math_FunctionSetWithDerivatives`
//! base ([`FunctionSetWithDerivatives`]); `occ::handle<Adaptor3d_Curve>`
//! maps to a `&Curve3` reference; `math_Vector` maps to `&[f64]` slices and
//! `math_Matrix` to the row-major `&mut [Vec<f64>]` Jacobian of the rcad
//! solver (OCCT D(i, j) -> d[i - 1][j - 1]).

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::brep_blend_func_inv::BlendCurvPointFuncInv;

/// OCCT BRepBlend_CurvPointRadInv — function of reframing between a point
/// and a curve, valid in cases of constant and progressive radius
/// (BRepBlend_CurvPointRadInv.hxx L41).  The vector X is w, U.
pub struct BRepBlendCurvPointRadInv<'a> {
    curv1: &'a Curve3,
    curv2: &'a Curve3,
    point: DVec3,
    choix: i32,
}

impl<'a> BRepBlendCurvPointRadInv<'a> {
    /// OCCT BRepBlend_CurvPointRadInv(C1, C2) (BRepBlend_CurvPointRadInv.cxx
    /// L24-30).
    pub fn new(c1: &'a Curve3, c2: &'a Curve3) -> Self {
        BRepBlendCurvPointRadInv {
            curv1: c1,
            curv2: c2,
            point: DVec3::ZERO,
            choix: 0,
        }
    }

    /// OCCT Set(Choix) (BRepBlend_CurvPointRadInv.cxx L34-37).
    pub fn set_choix(&mut self, choix: i32) {
        self.choix = choix;
    }

    /// OCCT NbEquations() (BRepBlend_CurvPointRadInv.cxx L41-44) — returns 2.
    pub fn nb_equations(&self) -> usize {
        2
    }

    /// OCCT Value(X, F) (BRepBlend_CurvPointRadInv.cxx L48-61).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: curv1->D1(X(1), ptcur1, d1cur1);
        let ptcur1 = self.curv1.point_at(x[0]);
        let d1cur1 = self.curv1.derivative_at(x[0]);
        let nplan = d1cur1.normalize_or_zero();
        let the_d = -nplan.dot(ptcur1);
        // OCCT: curv2->D1(X(2), ptcur2, d1cur2);
        let ptcur2 = self.curv2.point_at(x[1]);
        let d1cur2 = self.curv2.derivative_at(x[1]);
        f[0] = nplan.dot(self.point) + the_d;
        f[1] = nplan.dot(ptcur2) + the_d;
        let _ = d1cur2; // OCCT: d1cur2 is computed but not used in F.
        true
    }

    /// OCCT Derivatives(X, D) (BRepBlend_CurvPointRadInv.cxx L65-86).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        let (ptcur1, d1cur1, d2cur1) = {
            let p = self.curv1.point_at(x[0]);
            let d1 = self.curv1.derivative_at(x[0]);
            let d2 = self.curv1.derivative2_at(x[0]);
            (p, d1, d2)
        };

        let normd1cur1 = d1cur1.length();
        let unsurnormd1cur1 = 1.0 / normd1cur1;
        let nplan = unsurnormd1cur1 * d1cur1;
        // OCCT: dnplan.SetLinearForm(-nplan.Dot(d2cur1), nplan, d2cur1);
        let mut dnplan = -(nplan.dot(d2cur1)) * nplan + d2cur1;
        dnplan *= unsurnormd1cur1;
        let dthe_d = -nplan.dot(d1cur1) - dnplan.dot(ptcur1);
        d[0][0] = dnplan.dot(self.point) + dthe_d;
        d[0][1] = 0.0;
        let ptcur2 = self.curv2.point_at(x[1]);
        let d1cur2 = self.curv2.derivative_at(x[1]);
        d[1][0] = dnplan.dot(ptcur2) + dthe_d;
        d[1][1] = nplan.dot(d1cur2);

        true
    }

    /// OCCT Values(X, F, D) (BRepBlend_CurvPointRadInv.cxx L90-96).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        self.value(x, f);
        self.derivatives(x, d);
        true
    }

    /// OCCT Set(P) (BRepBlend_CurvPointRadInv.cxx L100-103).
    pub fn set_point(&mut self, p: DVec3) {
        self.point = p;
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BRepBlend_CurvPointRadInv.cxx
    /// L107-111) — Tolerance(1) = curv1->Resolution(Tol); Tolerance(2) =
    /// curv2->Resolution(Tol).  Pending: rcad Curve3 has no Resolution
    /// (Adaptor3d_Curve::Resolution / per-type Geom resolution is an
    /// untranslated adaptor-layer subsystem).
    pub fn get_tolerance(&self, _tolerance: &mut [f64], _tol: f64) {
        unimplemented!(
            "BRepBlend_CurvPointRadInv::GetTolerance: pending Adaptor3d_Curve::Resolution"
        );
    }

    /// OCCT GetBounds(InfBound, SupBound) (BRepBlend_CurvPointRadInv.cxx
    /// L115-121).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.curv1.default_domain()[0]; // OCCT: curv1->FirstParameter()
        sup_bound[0] = self.curv1.default_domain()[1]; // OCCT: curv1->LastParameter()
        inf_bound[1] = self.curv2.default_domain()[0]; // OCCT: curv2->FirstParameter()
        sup_bound[1] = self.curv2.default_domain()[1]; // OCCT: curv2->LastParameter()
    }

    /// OCCT IsSolution(Sol, Tol) (BRepBlend_CurvPointRadInv.cxx L125-129).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector valsol(1, 2); Value(Sol, valsol);
        let mut valsol = [0.0f64; 2];
        self.value(sol, &mut valsol);
        valsol[0].abs() <= tol && valsol[1].abs() <= tol
    }
}

impl<'a> FunctionSetWithDerivatives for BRepBlendCurvPointRadInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: NbVariables is implemented once at Blend_CurvPointFuncInv
        // level and inherited.
        BlendCurvPointFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        self.nb_equations()
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BRepBlendCurvPointRadInv::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BRepBlendCurvPointRadInv::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BRepBlendCurvPointRadInv::values(self, x, f, df)
    }
}

impl<'a> BlendCurvPointFuncInv for BRepBlendCurvPointRadInv<'a> {
    fn set_point(&mut self, p: DVec3) {
        BRepBlendCurvPointRadInv::set_point(self, p)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BRepBlendCurvPointRadInv::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BRepBlendCurvPointRadInv::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BRepBlendCurvPointRadInv::is_solution(self, sol, tol)
    }
}
