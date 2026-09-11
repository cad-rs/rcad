//! OCCT BRepBlend_SurfCurvConstRadInv (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_SurfCurvConstRadInv.hxx + BRepBlend_SurfCurvConstRadInv.cxx
//! (whole file L23-284).  Inverse function to find a solution on a curve on
//! surface for a constant-radius blend with a curve; the SimulSurf /
//! PerformSurf face/rst constant-radius arms instantiate it as `finvc`
//! (ChFi3d_FilBuilder.cxx L833 / L1763).
//!
//! Architecture mapping: `class BRepBlend_SurfCurvConstRadInv : public
//! Blend_SurfCurvFuncInv` is expressed by implementing the
//! [`BlendSurfCurvFuncInv`] trait over the `math_FunctionSetWithDerivatives`
//! base; `math_Vector` / `math_Matrix` map to `[f64; 3]` / `Vec<Vec<f64>>`
//! (OCCT X(i) -> x[i - 1], D(i, j) -> d[i - 1][j - 1]).
//!
//! The OCCT `occ::handle<Adaptor3d_Curve> curv` member is bound at every
//! ChFi3d construction site to an `Adaptor3d_CurveOnSurface` HC loaded with
//! (PC1, HS1) (ChFi3d_FilBuilder.cxx L831-833); the rcad port carries that
//! payload as the (curv_pcurve, curv_surf) pair and transcribes the consumed
//! Adaptor3d_Curve operations from Adaptor3d_CurveOnSurface.cxx (Value /
//! EvalD0, D1 / EvalD1 generic branch, FirstParameter / LastParameter
//! L977-987, Resolution L1364-1370) in the `curv_*` helpers below.
//!
//! Pending kernel dependency (marked GAP, plan 0.6):
//! Adaptor2d_Curve2d::Resolution.

use glam::DVec3;

use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending;
use super::brep_blend_func_inv::BlendSurfCurvFuncInv;

/// OCCT BRepBlend_SurfCurvConstRadInv.
pub struct BlendSurfCurvConstRadInv<'a> {
    // OCCT BRepBlend_SurfCurvConstRadInv.hxx fields.
    pub(crate) surf: &'a Surface3,
    // OCCT: occ::handle<Adaptor3d_Curve> curv — bound to the
    // Adaptor3d_CurveOnSurface(PC1, HS1); carried as the wrapped pair.
    pub(crate) curv_pcurve: &'a Curve2d,
    pub(crate) curv_surf: &'a Surface3,
    pub(crate) guide: &'a Curve3,
    pub(crate) rst: Option<&'a Curve2d>,
    pub(crate) ray: f64,
    pub(crate) choix: i32,
}

impl<'a> BlendSurfCurvConstRadInv<'a> {
    /// OCCT Adaptor3d_CurveOnSurface::Value bound to curv (the EvalD0
    /// generic branch of Adaptor3d_CurveOnSurface.cxx — the pcurve point
    /// carried on the wrapped surface).
    fn curv_value(&self, t: f64) -> DVec3 {
        let puv = self.curv_pcurve.point_at(t);
        self.curv_surf.point_at(puv.x, puv.y)
    }

    /// OCCT Adaptor3d_CurveOnSurface::D1 bound to curv
    /// (Adaptor3d_CurveOnSurface.cxx EvalD1 L1200-1241, the generic branch:
    /// `D1.SetLinearForm(Duv.X(), D1U, Duv.Y(), D1V)`).
    fn curv_d1(&self, t: f64) -> (DVec3, DVec3) {
        let puv = self.curv_pcurve.point_at(t);
        let duv = self.curv_pcurve.derivative_at(t);
        let (p, d1u, d1v) = self.curv_surf.derivatives(puv.x, puv.y);
        (p, duv.x * d1u + duv.y * d1v)
    }

    /// OCCT Adaptor3d_CurveOnSurface::Resolution bound to curv
    /// (Adaptor3d_CurveOnSurface.cxx L1364-1370): the pcurve resolution at
    /// min(u, v)-resolution of the wrapped surface.  GAP (plan 0.6): the
    /// final Adaptor2d_Curve2d::Resolution step has no rcad equivalent; the
    /// established pending marker preserves the OCCT failure path.
    fn curv_resolution(&self, r3d: f64) -> f64 {
        let ru = self.curv_surf.u_resolution(r3d);
        let rv = self.curv_surf.v_resolution(r3d);
        let _ = ru.min(rv);
        adaptor2d_curve2d_resolution_pending()
    }

    /// OCCT BRepBlend_SurfCurvConstRadInv(S, C, Cg) (…SurfCurvConstRadInv.cxx
    /// L23-33) — ray(0.0), choix(0); rst is null until Set(Rst).  The rcad
    /// constructor receives the Adaptor3d_CurveOnSurface payload
    /// (curv_pcurve, curv_surf) in place of the C handle.
    pub fn new(
        s: &'a Surface3,
        curv_pcurve: &'a Curve2d,
        curv_surf: &'a Surface3,
        cg: &'a Curve3,
    ) -> Self {
        BlendSurfCurvConstRadInv {
            surf: s,
            curv_pcurve,
            curv_surf,
            guide: cg,
            rst: None,
            ray: 0.0,
            choix: 0,
        }
    }

    /// OCCT Set(R, Choix) (…SurfCurvConstRadInv.cxx L37-56).
    pub fn set(&mut self, r: f64, choix: i32) {
        self.choix = choix;
        match self.choix {
            1 | 2 => {
                self.ray = -r.abs();
            }
            3 | 4 => {
                self.ray = r.abs();
            }
            _ => {
                self.ray = -r.abs();
            }
        }
    }

    /// OCCT NbEquations() (…SurfCurvConstRadInv.cxx L60-63) — returns 3.
    pub fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(X, F) (…SurfCurvConstRadInv.cxx L67-90).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L71: guide->D1(X(1), ptgui, d1gui);
        let ptgui = self.guide.point_at(x[0]);
        let d1gui = self.guide.derivative_at(x[0]);
        // OCCT L72-73.
        let nplan = d1gui.normalize_or_zero();
        let the_d = -nplan.dot(ptgui);
        // OCCT L74-75: ptcur = curv->Value(X(2)); F(1) = nplan . ptcur + theD;
        let ptcur = self.curv_value(x[1]);
        f[0] = nplan.dot(ptcur) + the_d;
        // OCCT L76-80: p2drst = rst->Value(X(3)); surf->D1(p2drst.X(), p2drst.Y(), ...);
        let rst = self.rst.expect("rst");
        let p2drst = rst.point_at(x[2]);
        let (pts, du, dv) = self.surf.derivatives(p2drst.x, p2drst.y);
        f[1] = nplan.dot(pts) + the_d;
        // OCCT L81-88.
        let mut ns = du.cross(dv);
        let norm = nplan.cross(ns).length();
        let unsurnorm = 1.0 / norm;
        ns = (nplan.dot(ns)) * nplan + (-1.0) * ns;
        ns *= unsurnorm;
        // OCCT L86-88: gp_Vec ref(ptcur, pts) = pts - ptcur;
        let refv = self.ray * ns + (pts - ptcur);
        f[2] = refv.length_squared() - self.ray * self.ray;
        true
    }

    /// OCCT Derivatives(X, D) (…SurfCurvConstRadInv.cxx L94-165).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L98: guide->D2(X(1), ptgui, d1gui, d2gui);
        let ptgui = self.guide.point_at(x[0]);
        let d1gui = self.guide.derivative_at(x[0]);
        let d2gui = self.guide.derivative2_at(x[0]);
        // OCCT L99-104.
        let normd1gui = d1gui.length();
        let unsurnormd1gui = 1.0 / normd1gui;
        let nplan = unsurnormd1gui * d1gui;
        let mut dnplan = (-nplan.dot(d2gui)) * nplan + d2gui;
        dnplan *= unsurnormd1gui;
        // OCCT L105.
        let dthe_d = -nplan.dot(d1gui) - dnplan.dot(ptgui);
        // OCCT L106-108: curv->D1(X(2), ptcur, d1cur);
        let ptcur = self.curv_value(x[1]);
        let d1cur = self.curv_d1(x[1]).1;
        // OCCT L109-111.
        d[0][0] = dnplan.dot(ptcur) + dthe_d;
        d[0][1] = nplan.dot(d1cur);
        d[0][2] = 0.0;

        // OCCT L113-118: rst->D1(X(3), p2drst, d1rst); surf->D2(...);
        let rst = self.rst.expect("rst");
        let p2drst = rst.point_at(x[2]);
        let d1rst = rst.derivative_at(x[2]);
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.derivatives2(p2drst.x, p2drst.y);
        // OCCT L119-123.
        d[1][0] = dnplan.dot(pts) + dthe_d;
        d[1][1] = 0.0;
        let dwrstpts = d1rst.x * d1u + d1rst.y * d1v;
        d[1][2] = nplan.dot(dwrstpts);

        // OCCT L125-129.
        let nsurf = d1u.cross(d1v);
        let dunsurf = d2u.cross(d1v) + d1u.cross(duv);
        let dvnsurf = d1u.cross(d2v) + duv.cross(d1v);
        let dwrstnsurf = d1rst.x * dunsurf + d1rst.y * dvnsurf;

        // OCCT L131-133.
        let nplancrosnsurf = nplan.cross(nsurf);
        let dwguinplancrosnsurf = dnplan.cross(nsurf);
        let dwrstnplancrosnsurf = nplan.cross(dwrstnsurf);

        // OCCT L135-142.
        let norm2 = nplancrosnsurf.length_squared();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwguinorm = unsurnorm * nplancrosnsurf.dot(dwguinplancrosnsurf);
        let dwrstnorm = unsurnorm * nplancrosnsurf.dot(dwrstnplancrosnsurf);

        // OCCT L144-146.
        let nplandotnsurf = nplan.dot(nsurf);
        let dwguinplandotnsurf = dnplan.dot(nsurf);
        let dwrstnplandotnsurf = nplan.dot(dwrstnsurf);

        // OCCT L148-151.
        let temp = nplandotnsurf * nplan + (-1.0) * nsurf;
        let dwguitemp = nplandotnsurf * dnplan + dwguinplandotnsurf * nplan;
        let dwrsttemp = dwrstnplandotnsurf * nplan + (-1.0) * dwrstnsurf;

        // OCCT L153-157: corde(ptcur, pts) = pts - ptcur.
        let corde = pts - ptcur;
        let mut refv = raysurnorm * temp + corde;
        let dwguiref = raysurnorm * dwguitemp + (-raysurnorm2 * dwguinorm) * temp;
        let dwrstref = raysurnorm * dwrsttemp + (-raysurnorm2 * dwrstnorm) * temp + dwrstpts;

        // OCCT L159-162.
        refv += refv;
        d[2][0] = refv.dot(dwguiref);
        d[2][1] = -refv.dot(d1cur);
        d[2][2] = refv.dot(dwrstref);

        true
    }

    /// OCCT Values(X, F, D) (…SurfCurvConstRadInv.cxx L169-243).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L173: guide->D2(X(1), ptgui, d1gui, d2gui);
        let ptgui = self.guide.point_at(x[0]);
        let d1gui = self.guide.derivative_at(x[0]);
        let d2gui = self.guide.derivative2_at(x[0]);
        // OCCT L174-181.
        let normd1gui = d1gui.length();
        let unsurnormd1gui = 1.0 / normd1gui;
        let nplan = unsurnormd1gui * d1gui;
        let the_d = -nplan.dot(ptgui);
        let mut dnplan = (-nplan.dot(d2gui)) * nplan + d2gui;
        dnplan *= unsurnormd1gui;
        let dthe_d = -nplan.dot(d1gui) - dnplan.dot(ptgui);
        // OCCT L182-184.
        let ptcur = self.curv_value(x[1]);
        let d1cur = self.curv_d1(x[1]).1;
        // OCCT L185-188.
        f[0] = nplan.dot(ptcur) + the_d;
        d[0][0] = dnplan.dot(ptcur) + dthe_d;
        d[0][1] = nplan.dot(d1cur);
        d[0][2] = 0.0;

        // OCCT L190-195.
        let rst = self.rst.expect("rst");
        let p2drst = rst.point_at(x[2]);
        let d1rst = rst.derivative_at(x[2]);
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.derivatives2(p2drst.x, p2drst.y);
        // OCCT L196-201.
        f[1] = nplan.dot(pts) + the_d;
        d[1][0] = dnplan.dot(pts) + dthe_d;
        d[1][1] = 0.0;
        let dwrstpts = d1rst.x * d1u + d1rst.y * d1v;
        d[1][2] = nplan.dot(dwrstpts);

        // OCCT L203-207.
        let nsurf = d1u.cross(d1v);
        let dunsurf = d2u.cross(d1v) + d1u.cross(duv);
        let dvnsurf = d1u.cross(d2v) + duv.cross(d1v);
        let dwrstnsurf = d1rst.x * dunsurf + d1rst.y * dvnsurf;

        // OCCT L209-211.
        let nplancrosnsurf = nplan.cross(nsurf);
        let dwguinplancrosnsurf = dnplan.cross(nsurf);
        let dwrstnplancrosnsurf = nplan.cross(dwrstnsurf);

        // OCCT L213-220.
        let norm2 = nplancrosnsurf.length_squared();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwguinorm = unsurnorm * nplancrosnsurf.dot(dwguinplancrosnsurf);
        let dwrstnorm = unsurnorm * nplancrosnsurf.dot(dwrstnplancrosnsurf);

        // OCCT L222-224.
        let nplandotnsurf = nplan.dot(nsurf);
        let dwguinplandotnsurf = dnplan.dot(nsurf);
        let dwrstnplandotnsurf = nplan.dot(dwrstnsurf);

        // OCCT L226-229.
        let temp = nplandotnsurf * nplan + (-1.0) * nsurf;
        let dwguitemp = nplandotnsurf * dnplan + dwguinplandotnsurf * nplan;
        let dwrsttemp = dwrstnplandotnsurf * nplan + (-1.0) * dwrstnsurf;

        // OCCT L231-236: corde(ptcur, pts) = pts - ptcur.
        let corde = pts - ptcur;
        let mut refv = raysurnorm * temp + corde;
        f[2] = refv.length_squared() - self.ray * self.ray;
        let dwguiref = raysurnorm * dwguitemp + (-raysurnorm2 * dwguinorm) * temp;
        let dwrstref = raysurnorm * dwrsttemp + (-raysurnorm2 * dwrstnorm) * temp + dwrstpts;

        // OCCT L238-242.
        refv += refv;
        d[2][0] = refv.dot(dwguiref);
        d[2][1] = -refv.dot(d1cur);
        d[2][2] = refv.dot(dwrstref);
        true
    }

    /// OCCT Set(Rst) (…SurfCurvConstRadInv.cxx L247-250).
    pub fn set_rst(&mut self, rst: &'a Curve2d) {
        self.rst = Some(rst);
    }

    /// OCCT GetTolerance(Tolerance, Tol) (…SurfCurvConstRadInv.cxx L254-262).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L256: Tolerance(1) = guide->Resolution(Tol);
        tolerance[0] = self.guide.resolution(tol);
        // OCCT L257: Tolerance(2) = curv->Resolution(Tol) — the
        // Adaptor3d_CurveOnSurface resolution bound to curv.
        tolerance[1] = self.curv_resolution(tol);
        // OCCT L258-261: Tolerance(3) = rst->Resolution(std::min(ru, rv)).
        // GAP (plan 0.6): the final Adaptor2d_Curve2d::Resolution step has no
        // rcad equivalent; the established pending marker preserves the OCCT
        // failure path.
        let ru = self.surf.u_resolution(tol);
        let rv = self.surf.v_resolution(tol);
        let _ = ru.min(rv);
        tolerance[2] = adaptor2d_curve2d_resolution_pending();
    }

    /// OCCT GetBounds(InfBound, SupBound) (…SurfCurvConstRadInv.cxx L266-274).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.guide.default_domain()[0]; // FirstParameter
        sup_bound[0] = self.guide.default_domain()[1]; // LastParameter
        // OCCT L270-271: curv->FirstParameter / LastParameter — the pcurve's
        // range (Adaptor3d_CurveOnSurface L977-987).
        inf_bound[1] = self.curv_pcurve.default_domain()[0]; // FirstParameter
        sup_bound[1] = self.curv_pcurve.default_domain()[1]; // LastParameter
        let rst = self.rst.expect("rst");
        inf_bound[2] = rst.default_domain()[0]; // FirstParameter
        sup_bound[2] = rst.default_domain()[1]; // LastParameter
    }

    /// OCCT IsSolution(Sol, Tol) (…SurfCurvConstRadInv.cxx L278-284).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 3];
        self.value(sol, &mut valsol);
        valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= 2.0 * tol * self.ray.abs()
    }
}

impl<'a> FunctionSetWithDerivatives for BlendSurfCurvConstRadInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: Blend_SurfCurvFuncInv::NbVariables (Blend_SurfCurvFuncInv.cxx
        // L19-27) — returns 3; dispatches to the trait override below.
        BlendSurfCurvFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendSurfCurvConstRadInv::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendSurfCurvConstRadInv::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfCurvConstRadInv::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfCurvConstRadInv::values(self, x, f, df)
    }
}

impl<'a> BlendSurfCurvFuncInv for BlendSurfCurvConstRadInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: Blend_SurfCurvFuncInv::NbVariables (Blend_SurfCurvFuncInv.cxx
        // L19-27) — returns 3.
        3
    }

    fn nb_equations(&self) -> usize {
        BlendSurfCurvConstRadInv::nb_equations(self)
    }

    fn set_rst(&mut self, rst: &Curve2d) {
        // OCCT SurfCurvConstRadInv.cxx L247-250; the rcad port stores the
        // curve reference (OCCT copies the handle).  SAFETY: the caller owns
        // the Curve2d for the lifetime 'a of this function object — same
        // invariant as the OCCT handle (pattern of BlendFunc_ConstRadInv).
        let r: &'a Curve2d = unsafe { &*(rst as *const Curve2d) };
        BlendSurfCurvConstRadInv::set_rst(self, r)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendSurfCurvConstRadInv::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendSurfCurvConstRadInv::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendSurfCurvConstRadInv::is_solution(self, sol, tol)
    }
}
