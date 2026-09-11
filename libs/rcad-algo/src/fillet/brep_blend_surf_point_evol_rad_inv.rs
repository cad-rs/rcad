//! OCCT BRepBlend_SurfPointEvolRadInv (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_SurfPointEvolRadInv.hxx (L36-95) + BRepBlend_SurfPointEvolRadInv.cxx
//! (whole file L24-269).  Inverse function to find a solution on the surface
//! for a variable-radius blend with a point of a curve; the SimulSurf /
//! PerformSurf rst/face and rst/rst variable-radius arms instantiate it as
//! `finvp` (ChFi3d_FilBuilder.cxx L1150 / L1394 / L2003 / L2197).
//!
//! Architecture mapping: `class BRepBlend_SurfPointEvolRadInv : public
//! Blend_SurfPointFuncInv` is expressed by implementing the
//! [`BlendSurfPointFuncInv`] trait over the `math_FunctionSetWithDerivatives`
//! base; `math_Vector` / `math_Matrix` map to `[f64; 3]` / `Vec<Vec<f64>>`
//! (OCCT X(i) -> x[i - 1], D(i, j) -> d[i - 1][j - 1]);
//! `occ::handle<Law_Function> tevol` maps to the shared
//! [`LawFunctionHandle`].

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use crate::geomalgo::law::law_function::LawFunctionHandle;

use super::brep_blend_func_inv::BlendSurfPointFuncInv;

/// OCCT BRepBlend_SurfPointEvolRadInv.
pub struct BlendSurfPointEvolRadInv<'a> {
    // OCCT BRepBlend_SurfPointEvolRadInv.hxx fields (L88-94).
    pub(crate) surf: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) ray: f64,
    pub(crate) choix: i32,
    pub(crate) tevol: LawFunctionHandle,
    pub(crate) sg1: f64,
    /// OCCT: gp_Pnt point — the point set by Set(P).
    pub(crate) point: DVec3,
}

impl<'a> BlendSurfPointEvolRadInv<'a> {
    /// OCCT BRepBlend_SurfPointEvolRadInv(S, C, Evol)
    /// (BRepBlend_SurfPointEvolRadInv.cxx L24-32).  The OCCT `point` member
    /// is uninitialized until Set(P); the rcad port zero-initializes it.
    pub fn new(s: &'a Surface3, c: &'a Curve3, evol: LawFunctionHandle) -> Self {
        BlendSurfPointEvolRadInv {
            surf: s,
            curv: c,
            ray: 0.0,
            choix: 0,
            tevol: evol, // OCCT L31: tevol = Evol.
            sg1: 0.0,
            point: DVec3::ZERO,
        }
    }

    /// OCCT Set(Choix) (BRepBlend_SurfPointEvolRadInv.cxx L36-53).
    pub fn set(&mut self, choix: i32) {
        self.choix = choix;
        match self.choix {
            1 | 2 => {
                self.sg1 = -1.0;
            }
            3 | 4 => {
                self.sg1 = 1.0;
            }
            _ => {
                self.sg1 = -1.0;
            }
        }
    }

    /// OCCT NbEquations() (BRepBlend_SurfPointEvolRadInv.cxx L57-60) —
    /// returns 3.
    pub fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(X, F) (BRepBlend_SurfPointEvolRadInv.cxx L64-86).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L70: curv->D1(X(1), ptcur, d1cur).
        let ptcur = self.curv.point_at(x[0]);
        let d1cur = self.curv.derivative_at(x[0]);
        // OCCT L71.
        self.ray = self.sg1 * self.tevol.borrow_mut().value(x[0]);
        // OCCT L72.
        let nplan = d1cur.normalize_or_zero();
        // OCCT L73.
        let the_d = -(nplan.dot(ptcur));
        // OCCT L74.
        let (pts, d1u, d1v) = self.surf.derivatives(x[1], x[2]);
        // OCCT L75-76.
        f[0] = nplan.dot(self.point) + the_d;
        f[1] = nplan.dot(pts) + the_d;
        // OCCT L77-81.
        let mut ns = d1u.cross(d1v);
        let norm = nplan.cross(ns).length();
        let unsurnorm = 1.0 / norm;
        ns = (nplan.dot(ns)) * nplan + (-1.0) * ns;
        ns *= unsurnorm;
        // OCCT L82-84: ref = pts - point; ref.SetLinearForm(ray, ns, ref).
        let refv = self.ray * ns + (pts - self.point);
        f[2] = refv.length_squared() - self.ray * self.ray;
        true
    }

    /// OCCT Derivatives(X, D) (BRepBlend_SurfPointEvolRadInv.cxx L90-156).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L96: curv->D2(X(1), ptcur, d1cur, d2cur).
        let ptcur = self.curv.point_at(x[0]);
        let d1cur = self.curv.derivative_at(x[0]);
        let d2cur = self.curv.derivative2_at(x[0]);
        // OCCT L97-100.
        let (mut ray, mut dray) = (0.0f64, 0.0f64);
        self.tevol.borrow_mut().d1(x[0], &mut ray, &mut dray);
        ray = self.sg1 * ray;
        dray = self.sg1 * dray;
        self.ray = ray;
        // OCCT L101-103.
        let normd1cur = d1cur.length();
        let unsurnormd1cur = 1.0 / normd1cur;
        let nplan = unsurnormd1cur * d1cur;
        // OCCT L103-105.
        let mut dnplan = (-nplan.dot(d2cur)) * nplan + d2cur;
        dnplan *= unsurnormd1cur;
        // OCCT L105.
        let dthe_d = -nplan.dot(d1cur) - dnplan.dot(ptcur);
        // OCCT L106-107.
        d[0][0] = dnplan.dot(self.point) + dthe_d;
        d[0][1] = 0.0;
        d[0][2] = 0.0;
        // OCCT L108-111.
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.derivatives2(x[1], x[2]);
        d[1][0] = dnplan.dot(pts) + dthe_d;
        d[1][1] = nplan.dot(d1u);
        d[1][2] = nplan.dot(d1v);

        // OCCT L113-116.
        let nsurf = d1u.cross(d1v);
        let dunsurf = d2u.cross(d1v) + d1u.cross(duv);
        let dvnsurf = d1u.cross(d2v) + duv.cross(d1v);

        // OCCT L117-120.
        let nplancrosnsurf = nplan.cross(nsurf);
        let dwnplancrosnsurf = dnplan.cross(nsurf);
        let dunplancrosnsurf = nplan.cross(dunsurf);
        let dvnplancrosnsurf = nplan.cross(dvnsurf);

        // OCCT L122-130.
        let norm2 = nplancrosnsurf.length_squared();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = ray * unsurnorm2;
        let dwnorm = unsurnorm * nplancrosnsurf.dot(dwnplancrosnsurf);
        let dunorm = unsurnorm * nplancrosnsurf.dot(dunplancrosnsurf);
        let dvnorm = unsurnorm * nplancrosnsurf.dot(dvnplancrosnsurf);

        // OCCT L132-135.
        let nplandotnsurf = nplan.dot(nsurf);
        let dwnplandotnsurf = dnplan.dot(nsurf);
        let dunplandotnsurf = nplan.dot(dunsurf);
        let dvnplandotnsurf = nplan.dot(dvnsurf);

        // OCCT L137-141.
        let temp = nplandotnsurf * nplan + (-1.0) * nsurf;
        let mut dwtemp = nplandotnsurf * dnplan + dwnplandotnsurf * nplan;
        let dutemp = dunplandotnsurf * nplan + (-1.0) * dunsurf;
        let dvtemp = dvnplandotnsurf * nplan + (-1.0) * dvnsurf;

        // OCCT L143-148: corde(point, pts) = pts - point.
        let corde = pts - self.point;
        let mut refv = raysurnorm * temp + corde;
        dwtemp = raysurnorm * dwtemp + (-raysurnorm2 * dwnorm) * temp;
        // OCCT L146: dwref.SetLinearForm(1., dwref, dray * unsurnorm, temp).
        dwtemp = dwtemp + dray * unsurnorm * temp;
        let duref = raysurnorm * dutemp + (-raysurnorm2 * dunorm) * temp + d1u;
        let dvref = raysurnorm * dvtemp + (-raysurnorm2 * dvnorm) * temp + d1v;

        // OCCT L150-153: ref.Add(ref).
        refv += refv;
        d[2][0] = refv.dot(dwtemp) - 2.0 * dray * ray;
        d[2][1] = refv.dot(duref);
        d[2][2] = refv.dot(dvref);

        true
    }

    /// OCCT Values(X, F, D) (BRepBlend_SurfPointEvolRadInv.cxx L160-231).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L166: curv->D2(X(1), ptcur, d1cur, d2cur).
        let ptcur = self.curv.point_at(x[0]);
        let d1cur = self.curv.derivative_at(x[0]);
        let d2cur = self.curv.derivative2_at(x[0]);
        // OCCT L167-169.
        let (mut ray, mut dray) = (0.0f64, 0.0f64);
        self.tevol.borrow_mut().d1(x[0], &mut ray, &mut dray);
        ray = self.sg1 * ray;
        dray = self.sg1 * dray;
        self.ray = ray;
        // OCCT L170.
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.derivatives2(x[1], x[2]);
        // OCCT L171-175.
        let normd1cur = d1cur.length();
        let unsurnormd1cur = 1.0 / normd1cur;
        let nplan = unsurnormd1cur * d1cur;
        let the_d = -(nplan.dot(ptcur));
        // OCCT L175-176.
        f[0] = nplan.dot(self.point) + the_d;
        f[1] = nplan.dot(pts) + the_d;

        // OCCT L178-185.
        let mut dnplan = (-nplan.dot(d2cur)) * nplan + d2cur;
        dnplan *= unsurnormd1cur;
        let dthe_d = -nplan.dot(d1cur) - dnplan.dot(ptcur);
        d[0][0] = dnplan.dot(self.point) + dthe_d;
        d[0][1] = 0.0;
        d[0][2] = 0.0;
        d[1][0] = dnplan.dot(pts) + dthe_d;
        d[1][1] = nplan.dot(d1u);
        d[1][2] = nplan.dot(d1v);

        // OCCT L187-190.
        let nsurf = d1u.cross(d1v);
        let dunsurf = d2u.cross(d1v) + d1u.cross(duv);
        let dvnsurf = d1u.cross(d2v) + duv.cross(d1v);

        // OCCT L191-194.
        let nplancrosnsurf = nplan.cross(nsurf);
        let dwnplancrosnsurf = dnplan.cross(nsurf);
        let dunplancrosnsurf = nplan.cross(dunsurf);
        let dvnplancrosnsurf = nplan.cross(dvnsurf);

        // OCCT L196-204.
        let norm2 = nplancrosnsurf.length_squared();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = ray * unsurnorm2;
        let dwnorm = unsurnorm * nplancrosnsurf.dot(dwnplancrosnsurf);
        let dunorm = unsurnorm * nplancrosnsurf.dot(dunplancrosnsurf);
        let dvnorm = unsurnorm * nplancrosnsurf.dot(dvnplancrosnsurf);

        // OCCT L206-209.
        let nplandotnsurf = nplan.dot(nsurf);
        let dwnplandotnsurf = dnplan.dot(nsurf);
        let dunplandotnsurf = nplan.dot(dunsurf);
        let dvnplandotnsurf = nplan.dot(dvnsurf);

        // OCCT L211-215.
        let temp = nplandotnsurf * nplan + (-1.0) * nsurf;
        let mut dwtemp = nplandotnsurf * dnplan + dwnplandotnsurf * nplan;
        let dutemp = dunplandotnsurf * nplan + (-1.0) * dunsurf;
        let dvtemp = dvnplandotnsurf * nplan + (-1.0) * dvnsurf;

        // OCCT L217-223: corde(point, pts) = pts - point.
        let corde = pts - self.point;
        let mut refv = raysurnorm * temp + corde;
        f[2] = refv.length_squared() - ray * ray;
        dwtemp = raysurnorm * dwtemp + (-raysurnorm2 * dwnorm) * temp;
        // OCCT L221: dwref.SetLinearForm(1., dwref, dray * unsurnorm, temp).
        dwtemp = dwtemp + dray * unsurnorm * temp;
        let duref = raysurnorm * dutemp + (-raysurnorm2 * dunorm) * temp + d1u;
        let dvref = raysurnorm * dvtemp + (-raysurnorm2 * dvnorm) * temp + d1v;

        // OCCT L225-228: ref.Add(ref).
        refv += refv;
        d[2][0] = refv.dot(dwtemp) - 2.0 * dray * ray;
        d[2][1] = refv.dot(duref);
        d[2][2] = refv.dot(dvref);

        true
    }

    /// OCCT Set(P) (BRepBlend_SurfPointEvolRadInv.cxx L235-238).
    pub fn set_point(&mut self, p: DVec3) {
        self.point = p;
    }

    /// OCCT GetTolerance(Tolerance, Tol)
    /// (BRepBlend_SurfPointEvolRadInv.cxx L242-247).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        tolerance[0] = self.curv.resolution(tol);
        tolerance[1] = self.surf.u_resolution(tol);
        tolerance[2] = self.surf.v_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound)
    /// (BRepBlend_SurfPointEvolRadInv.cxx L251-259).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.curv.default_domain()[0]; // FirstParameter
        sup_bound[0] = self.curv.default_domain()[1]; // LastParameter
        inf_bound[1] = self.surf.default_domain()[0]; // FirstUParameter
        sup_bound[1] = self.surf.default_domain()[1]; // LastUParameter
        inf_bound[2] = self.surf.default_domain()[2]; // FirstVParameter
        sup_bound[2] = self.surf.default_domain()[3]; // LastVParameter
    }

    /// OCCT IsSolution(Sol, Tol)
    /// (BRepBlend_SurfPointEvolRadInv.cxx L263-269).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 3];
        self.value(sol, &mut valsol);
        valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= 2.0 * tol * self.ray.abs()
    }
}

impl<'a> FunctionSetWithDerivatives for BlendSurfPointEvolRadInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: Blend_SurfPointFuncInv::NbVariables (Blend_SurfPointFuncInv.cxx
        // L19-22) — returns 3 for this derived function; dispatches to the
        // trait override below.
        BlendSurfPointFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendSurfPointEvolRadInv::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendSurfPointEvolRadInv::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfPointEvolRadInv::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfPointEvolRadInv::values(self, x, f, df)
    }
}

impl<'a> BlendSurfPointFuncInv for BlendSurfPointEvolRadInv<'a> {
    /// OCCT Blend_SurfPointFuncInv::NbVariables — the OCCT base returns 2
    /// but the const-rad sibling overrides it to 3 for the (U, V) surface
    /// unknowns; the same override is carried here (pattern of
    /// BlendSurfPointConstRadInv).
    fn nb_variables(&self) -> usize {
        3
    }

    fn nb_equations(&self) -> usize {
        BlendSurfPointEvolRadInv::nb_equations(self)
    }

    fn set_point(&mut self, p: DVec3) {
        BlendSurfPointEvolRadInv::set_point(self, p)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendSurfPointEvolRadInv::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendSurfPointEvolRadInv::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendSurfPointEvolRadInv::is_solution(self, sol, tol)
    }
}
