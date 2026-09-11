//! OCCT BRepBlend_CurvPointRadInv (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_CurvPointRadInv.hxx + BRepBlend_CurvPointRadInv.cxx
//! (whole file L24-130), bound at the ChFi3d rst/rst call sites to an
//! `Adaptor3d_CurveOnSurface` restriction curve (ChFi3d_FilBuilder.cxx
//! L1308 / L1310 / L2119 / L2121: `BRepBlend_CurvPointRadInv finvp1(HGuide,
//! HC2)` with `HC2->Load(PC2, HS2)`).
//!
//! This is the same OCCT class as
//! [`super::brep_blend_curv_point_rad_inv::BRepBlendCurvPointRadInv`]; the
//! plain-curve port models `curv2` as `&Curve3`, while the ChFi3d rst/rst
//! arms bind an Adaptor3d_CurveOnSurface handle.  Following the established
//! HC-payload architecture mapping (pattern of
//! [`super::brep_blend_surf_curv_const_rad_inv`]), this port carries that
//! payload as the (curv2_pcurve, curv2_surf) pair and transcribes the
//! consumed Adaptor3d_CurveOnSurface operations (D1/EvalD1 generic branch,
//! FirstParameter/LastParameter L977-987, Resolution L1364-1370) in the
//! `curv2_*` helpers.
//!
//! Pending kernel dependency (marked GAP, plan 0.6):
//! Adaptor3d_Curve::Resolution for curv1 (same pending as the plain-curve
//! port) and Adaptor2d_Curve2d::Resolution for the curv2 chain.

use glam::DVec3;

use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending;
use super::brep_blend_func_inv::BlendCurvPointFuncInv;

/// OCCT BRepBlend_CurvPointRadInv bound to an Adaptor3d_CurveOnSurface
/// restriction — function of reframing between a point and a curve, valid
/// in cases of constant and progressive radius (hxx L41).  The vector X is
/// w, U.
pub struct BRepBlendCurvPointRadInvHc<'a> {
    // OCCT: occ::handle<Adaptor3d_Curve> curv1 — the guide (elspine) curve.
    pub(crate) curv1: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> curv2 — bound to the
    // Adaptor3d_CurveOnSurface(PC, HS); carried as the wrapped pair.
    pub(crate) curv2_pcurve: &'a Curve2d,
    pub(crate) curv2_surf: &'a Surface3,
    pub(crate) point: DVec3,
    pub(crate) choix: i32,
}

impl<'a> BRepBlendCurvPointRadInvHc<'a> {
    /// OCCT Adaptor3d_CurveOnSurface::D1 bound to curv2
    /// (Adaptor3d_CurveOnSurface.cxx EvalD1 L1200-1241, the generic branch:
    /// `D1.SetLinearForm(Duv.X(), D1U, Duv.Y(), D1V)`).
    fn curv2_d1(&self, t: f64) -> (DVec3, DVec3) {
        let puv = self.curv2_pcurve.point_at(t);
        let duv = self.curv2_pcurve.derivative_at(t);
        let (p, d1u, d1v) = self.curv2_surf.derivatives(puv.x, puv.y);
        (p, duv.x * d1u + duv.y * d1v)
    }

    /// OCCT Adaptor3d_CurveOnSurface::Resolution bound to curv2
    /// (Adaptor3d_CurveOnSurface.cxx L1364-1370).  GAP (plan 0.6): the
    /// final Adaptor2d_Curve2d::Resolution step has no rcad equivalent; the
    /// established pending marker preserves the OCCT failure path.  The
    /// helper is reached once the GetTolerance pending Adaptor3d_Curve::
    /// Resolution gap closes (see below); kept with its OCCT anchor.
    #[allow(dead_code)]
    fn curv2_resolution(&self, r3d: f64) -> f64 {
        let ru = self.curv2_surf.u_resolution(r3d);
        let rv = self.curv2_surf.v_resolution(r3d);
        let _ = ru.min(rv);
        adaptor2d_curve2d_resolution_pending()
    }

    /// OCCT BRepBlend_CurvPointRadInv(C1, C2)
    /// (BRepBlend_CurvPointRadInv.cxx L24-30) — choix(0); the OCCT `point`
    /// member is uninitialized until Set(P); the rcad port zero-initializes
    /// it.  The rcad constructor receives the Adaptor3d_CurveOnSurface
    /// payload (curv2_pcurve, curv2_surf) in place of the C2 handle.
    pub fn new(c1: &'a Curve3, curv2_pcurve: &'a Curve2d, curv2_surf: &'a Surface3) -> Self {
        BRepBlendCurvPointRadInvHc {
            curv1: c1,
            curv2_pcurve,
            curv2_surf,
            point: DVec3::ZERO,
            choix: 0,
        }
    }

    /// OCCT Set(Choix) (BRepBlend_CurvPointRadInv.cxx L34-37).
    pub fn set_choix(&mut self, choix: i32) {
        self.choix = choix;
    }

    /// OCCT NbEquations() (BRepBlend_CurvPointRadInv.cxx L41-44) —
    /// returns 2.
    pub fn nb_equations(&self) -> usize {
        2
    }

    /// OCCT Value(X, F) (BRepBlend_CurvPointRadInv.cxx L48-61).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L54-56: curv1->D1(X(1), ptcur1, d1cur1); nplan =
        // d1cur1.Normalized().XYZ(); theD = -(nplan.Dot(ptcur1.XYZ())).
        let ptcur1 = self.curv1.point_at(x[0]);
        let d1cur1 = self.curv1.derivative_at(x[0]);
        let nplan = d1cur1.normalize_or_zero();
        let the_d = -(nplan.dot(ptcur1));
        // OCCT L57: curv2->D1(X(2), ptcur2, d1cur2).
        let (ptcur2, _d1cur2) = self.curv2_d1(x[1]);
        // OCCT L58-59.
        f[0] = nplan.dot(self.point) + the_d;
        f[1] = nplan.dot(ptcur2) + the_d;
        true
    }

    /// OCCT Derivatives(X, D) (BRepBlend_CurvPointRadInv.cxx L65-86).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L71: curv1->D2(X(1), ptcur1, d1cur1, d2cur1).
        let ptcur1 = self.curv1.point_at(x[0]);
        let d1cur1 = self.curv1.derivative_at(x[0]);
        let d2cur1 = self.curv1.derivative2_at(x[0]);

        // OCCT L73-77.
        let normd1cur1 = d1cur1.length();
        let unsurnormd1cur1 = 1.0 / normd1cur1;
        let nplan = unsurnormd1cur1 * d1cur1;
        // OCCT L76-77: dnplan.SetLinearForm(-nplan.Dot(d2cur1), nplan, d2cur1);
        // dnplan.Multiply(unsurnormd1cur1).
        let mut dnplan = -(nplan.dot(d2cur1)) * nplan + d2cur1;
        dnplan *= unsurnormd1cur1;
        // OCCT L78.
        let dthe_d = -nplan.dot(d1cur1) - dnplan.dot(ptcur1);
        // OCCT L79-80.
        d[0][0] = dnplan.dot(self.point) + dthe_d;
        d[0][1] = 0.0;
        // OCCT L81-83.
        let (ptcur2, d1cur2) = self.curv2_d1(x[1]);
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

    /// OCCT GetTolerance(Tolerance, Tol)
    /// (BRepBlend_CurvPointRadInv.cxx L107-111) — Tolerance(1) =
    /// curv1->Resolution(Tol); Tolerance(2) = curv2->Resolution(Tol).
    /// Pending: rcad Curve3 has no Resolution (Adaptor3d_Curve::Resolution /
    /// per-type Geom resolution is an untranslated adaptor-layer subsystem;
    /// same pending as the plain-curve port), and the curv2 chain ends at
    /// Adaptor2d_Curve2d::Resolution (pending marker).
    pub fn get_tolerance(&self, _tolerance: &mut [f64], _tol: f64) {
        unimplemented!(
            "BRepBlend_CurvPointRadInv::GetTolerance: pending Adaptor3d_Curve::Resolution"
        );
    }

    /// OCCT GetBounds(InfBound, SupBound)
    /// (BRepBlend_CurvPointRadInv.cxx L115-121) — the curv2 FirstParameter /
    /// LastParameter are the pcurve's range (Adaptor3d_CurveOnSurface
    /// L977-987).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.curv1.default_domain()[0]; // OCCT: curv1->FirstParameter()
        sup_bound[0] = self.curv1.default_domain()[1]; // OCCT: curv1->LastParameter()
        inf_bound[1] = self.curv2_pcurve.default_domain()[0]; // OCCT: curv2->FirstParameter()
        sup_bound[1] = self.curv2_pcurve.default_domain()[1]; // OCCT: curv2->LastParameter()
    }

    /// OCCT IsSolution(Sol, Tol) (BRepBlend_CurvPointRadInv.cxx L125-129).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector valsol(1, 2); Value(Sol, valsol).
        let mut valsol = [0.0f64; 2];
        self.value(sol, &mut valsol);
        valsol[0].abs() <= tol && valsol[1].abs() <= tol
    }
}

impl<'a> FunctionSetWithDerivatives for BRepBlendCurvPointRadInvHc<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: NbVariables is implemented once at Blend_CurvPointFuncInv
        // level and inherited.
        BlendCurvPointFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        self.nb_equations()
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BRepBlendCurvPointRadInvHc::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BRepBlendCurvPointRadInvHc::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BRepBlendCurvPointRadInvHc::values(self, x, f, df)
    }
}

impl<'a> BlendCurvPointFuncInv for BRepBlendCurvPointRadInvHc<'a> {
    fn set_point(&mut self, p: DVec3) {
        BRepBlendCurvPointRadInvHc::set_point(self, p)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BRepBlendCurvPointRadInvHc::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BRepBlendCurvPointRadInvHc::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BRepBlendCurvPointRadInvHc::is_solution(self, sol, tol)
    }
}
