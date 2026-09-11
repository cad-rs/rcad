//! OCCT BRepBlend_SurfPointConstRadInv (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_SurfPointConstRadInv.hxx + BRepBlend_SurfPointConstRadInv.cxx
//! (whole file L23-276).  Inverse function to find a solution on the surface
//! for a constant-radius blend with a point of a curve; the SimulSurf /
//! PerformSurf face/rst constant-radius arms instantiate it as `finvp`
//! (ChFi3d_FilBuilder.cxx L838 / L1768).
//!
//! Architecture mapping: `class BRepBlend_SurfPointConstRadInv : public
//! Blend_SurfPointFuncInv` is expressed by implementing the
//! [`BlendSurfPointFuncInv`] trait over the `math_FunctionSetWithDerivatives`
//! base; `math_Vector` / `math_Matrix` map to `[f64; 3]` / `Vec<Vec<f64>>`
//! (OCCT X(i) -> x[i - 1], D(i, j) -> d[i - 1][j - 1]).

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::brep_blend_func_inv::BlendSurfPointFuncInv;

/// OCCT BRepBlend_SurfPointConstRadInv.
pub struct BlendSurfPointConstRadInv<'a> {
    // OCCT BRepBlend_SurfPointConstRadInv.hxx fields.
    pub(crate) surf: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) ray: f64,
    pub(crate) choix: i32,
    /// OCCT: gp_Pnt point — the point set by Set(P).
    pub(crate) point: DVec3,
}

impl<'a> BlendSurfPointConstRadInv<'a> {
    /// OCCT BRepBlend_SurfPointConstRadInv(S, C) (…SurfPointConstRadInv.cxx
    /// L23-31) — ray(0.0), choix(0); the OCCT `point` member is
    /// uninitialized until Set(P); the rcad port zero-initializes it.
    pub fn new(s: &'a Surface3, c: &'a Curve3) -> Self {
        BlendSurfPointConstRadInv {
            surf: s,
            curv: c,
            ray: 0.0,
            choix: 0,
            point: DVec3::ZERO,
        }
    }

    /// OCCT Set(R, Choix) (…SurfPointConstRadInv.cxx L35-54).
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

    /// OCCT NbEquations() (…SurfPointConstRadInv.cxx L58-61) — returns 3.
    pub fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(X, F) (…SurfPointConstRadInv.cxx L65-90).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L71: curv->D1(X(1), ptcur, d1cur);
        let ptcur = self.curv.point_at(x[0]);
        let d1cur = self.curv.derivative_at(x[0]);
        // OCCT L72: nplan = d1cur.Normalized().XYZ();
        let nplan = d1cur.normalize_or_zero();
        // OCCT L74-76: theD = nplan.Dot(ptcur); theD = theD * (-1.);
        let the_d = nplan.dot(ptcur) * -1.0;

        // OCCT L78: surf->D1(X(2), X(3), pts, d1u, d1v);
        let (pts, d1u, d1v) = self.surf.derivatives(x[1], x[2]);
        // OCCT L79-80.
        f[0] = nplan.dot(self.point) + the_d;
        f[1] = nplan.dot(pts) + the_d;
        // OCCT L81-85.
        let mut ns = d1u.cross(d1v);
        let norm = nplan.cross(ns).length();
        let unsurnorm = 1.0 / norm;
        ns = (nplan.dot(ns)) * nplan + (-1.0) * ns;
        ns *= unsurnorm;
        // OCCT L86-88: ref = pts - point; ref.SetLinearForm(ray, ns, ref);
        let refv = self.ray * ns + (pts - self.point);
        f[2] = refv.length_squared() - self.ray * self.ray;
        true
    }

    /// OCCT Derivatives(X, D) (…SurfPointConstRadInv.cxx L94-162).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L100: curv->D2(X(1), ptcur, d1cur, d2cur);
        let ptcur = self.curv.point_at(x[0]);
        let d1cur = self.curv.derivative_at(x[0]);
        let d2cur = self.curv.derivative2_at(x[0]);
        // OCCT L101-103.
        let normd1cur = d1cur.length();
        let unsurnormd1cur = 1.0 / normd1cur;
        let nplan = unsurnormd1cur * d1cur;
        // OCCT L105-108.  (theD is a dead store in the OCCT Derivatives body
        // — computed at L107-108, never read; kept with the underscore for
        // statement parity.)
        let _the_d = nplan.dot(ptcur) * -1.0;

        // OCCT L110-111: dnplan.SetLinearForm(-nplan.Dot(d2cur), nplan, d2cur);
        // dnplan.Multiply(unsurnormd1cur);
        let mut dnplan = (-nplan.dot(d2cur)) * nplan + d2cur;
        dnplan *= unsurnormd1cur;
        // OCCT L112.
        let dthe_d = -nplan.dot(d1cur) - dnplan.dot(ptcur);
        // OCCT L113-114.
        d[0][0] = dnplan.dot(self.point) + dthe_d;
        d[0][1] = 0.0;
        d[0][2] = 0.0;
        // OCCT L115-118.
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.derivatives2(x[1], x[2]);
        d[1][0] = dnplan.dot(pts) + dthe_d;
        d[1][1] = nplan.dot(d1u);
        d[1][2] = nplan.dot(d1v);

        // OCCT L120-127.
        let nsurf = d1u.cross(d1v);
        let dunsurf = d2u.cross(d1v) + d1u.cross(duv);
        let dvnsurf = d1u.cross(d2v) + duv.cross(d1v);

        let nplancrosnsurf = nplan.cross(nsurf);
        let dwnplancrosnsurf = dnplan.cross(nsurf);
        let dunplancrosnsurf = nplan.cross(dunsurf);
        let dvnplancrosnsurf = nplan.cross(dvnsurf);

        // OCCT L129-137.
        let norm2 = nplancrosnsurf.length_squared();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwnorm = unsurnorm * nplancrosnsurf.dot(dwnplancrosnsurf);
        let dunorm = unsurnorm * nplancrosnsurf.dot(dunplancrosnsurf);
        let dvnorm = unsurnorm * nplancrosnsurf.dot(dvnplancrosnsurf);

        // OCCT L139-142.
        let nplandotnsurf = nplan.dot(nsurf);
        let dwnplandotnsurf = dnplan.dot(nsurf);
        let dunplandotnsurf = nplan.dot(dunsurf);
        let dvnplandotnsurf = nplan.dot(dvnsurf);

        // OCCT L144-148.
        let temp = nplandotnsurf * nplan + (-1.0) * nsurf;
        let dwtemp = nplandotnsurf * dnplan + dwnplandotnsurf * nplan;
        let dutemp = dunplandotnsurf * nplan + (-1.0) * dunsurf;
        let dvtemp = dvnplandotnsurf * nplan + (-1.0) * dvnsurf;

        // OCCT L150-154: corde(point, pts) = pts - point.
        let corde = pts - self.point;
        let mut refv = raysurnorm * temp + corde;
        let dwref = raysurnorm * dwtemp + (-raysurnorm2 * dwnorm) * temp;
        let duref = raysurnorm * dutemp + (-raysurnorm2 * dunorm) * temp + d1u;
        let dvref = raysurnorm * dvtemp + (-raysurnorm2 * dvnorm) * temp + d1v;

        // OCCT L156-159: ref.Add(ref);
        refv += refv;
        d[2][0] = refv.dot(dwref);
        d[2][1] = refv.dot(duref);
        d[2][2] = refv.dot(dvref);

        true
    }

    /// OCCT Values(X, F, D) (…SurfPointConstRadInv.cxx L166-238).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L172-173.
        let ptcur = self.curv.point_at(x[0]);
        let d1cur = self.curv.derivative_at(x[0]);
        let d2cur = self.curv.derivative2_at(x[0]);
        let (pts, d1u, d1v, d2u, d2v, duv) = self.surf.derivatives2(x[1], x[2]);
        // OCCT L174-176.
        let normd1cur = d1cur.length();
        let unsurnormd1cur = 1.0 / normd1cur;
        let nplan = unsurnormd1cur * d1cur;
        // OCCT L178-181.
        let the_d = nplan.dot(ptcur) * -1.0;

        // OCCT L183-184.
        f[0] = nplan.dot(self.point) + the_d;
        f[1] = nplan.dot(pts) + the_d;

        // OCCT L186-190.
        let mut dnplan = (-nplan.dot(d2cur)) * nplan + d2cur;
        dnplan *= unsurnormd1cur;
        let dthe_d = -nplan.dot(d1cur) - dnplan.dot(ptcur);
        d[0][0] = dnplan.dot(self.point) + dthe_d;
        d[0][1] = 0.0;
        d[0][2] = 0.0;
        // OCCT L191-193.
        d[1][0] = dnplan.dot(pts) + dthe_d;
        d[1][1] = nplan.dot(d1u);
        d[1][2] = nplan.dot(d1v);

        // OCCT L195-202.
        let nsurf = d1u.cross(d1v);
        let dunsurf = d2u.cross(d1v) + d1u.cross(duv);
        let dvnsurf = d1u.cross(d2v) + duv.cross(d1v);

        let nplancrosnsurf = nplan.cross(nsurf);
        let dwnplancrosnsurf = dnplan.cross(nsurf);
        let dunplancrosnsurf = nplan.cross(dunsurf);
        let dvnplancrosnsurf = nplan.cross(dvnsurf);

        // OCCT L204-212.
        let norm2 = nplancrosnsurf.length_squared();
        let norm = norm2.sqrt();
        let unsurnorm = 1.0 / norm;
        let raysurnorm = self.ray * unsurnorm;
        let unsurnorm2 = unsurnorm * unsurnorm;
        let raysurnorm2 = self.ray * unsurnorm2;
        let dwnorm = unsurnorm * nplancrosnsurf.dot(dwnplancrosnsurf);
        let dunorm = unsurnorm * nplancrosnsurf.dot(dunplancrosnsurf);
        let dvnorm = unsurnorm * nplancrosnsurf.dot(dvnplancrosnsurf);

        // OCCT L214-217.
        let nplandotnsurf = nplan.dot(nsurf);
        let dwnplandotnsurf = dnplan.dot(nsurf);
        let dunplandotnsurf = nplan.dot(dunsurf);
        let dvnplandotnsurf = nplan.dot(dvnsurf);

        // OCCT L219-223.
        let temp = nplandotnsurf * nplan + (-1.0) * nsurf;
        let dwtemp = nplandotnsurf * dnplan + dwnplandotnsurf * nplan;
        let dutemp = dunplandotnsurf * nplan + (-1.0) * dunsurf;
        let dvtemp = dvnplandotnsurf * nplan + (-1.0) * dvnsurf;

        // OCCT L225-230: corde(point, pts) = pts - point.
        let corde = pts - self.point;
        let mut refv = raysurnorm * temp + corde;
        f[2] = refv.length_squared() - self.ray * self.ray;
        let dwref = raysurnorm * dwtemp + (-raysurnorm2 * dwnorm) * temp;
        let duref = raysurnorm * dutemp + (-raysurnorm2 * dunorm) * temp + d1u;
        let dvref = raysurnorm * dvtemp + (-raysurnorm2 * dvnorm) * temp + d1v;

        // OCCT L232-235.
        refv += refv;
        d[2][0] = refv.dot(dwref);
        d[2][1] = refv.dot(duref);
        d[2][2] = refv.dot(dvref);

        true
    }

    /// OCCT Set(P) (…SurfPointConstRadInv.cxx L242-245).
    pub fn set_point(&mut self, p: DVec3) {
        self.point = p;
    }

    /// OCCT GetTolerance(Tolerance, Tol) (…SurfPointConstRadInv.cxx L249-254).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        tolerance[0] = self.curv.resolution(tol);
        tolerance[1] = self.surf.u_resolution(tol);
        tolerance[2] = self.surf.v_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (…SurfPointConstRadInv.cxx
    /// L258-266).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.curv.default_domain()[0]; // FirstParameter
        sup_bound[0] = self.curv.default_domain()[1]; // LastParameter
        inf_bound[1] = self.surf.default_domain()[0]; // FirstUParameter
        sup_bound[1] = self.surf.default_domain()[1]; // LastUParameter
        inf_bound[2] = self.surf.default_domain()[2]; // FirstVParameter
        sup_bound[2] = self.surf.default_domain()[3]; // LastVParameter
    }

    /// OCCT IsSolution(Sol, Tol) (…SurfPointConstRadInv.cxx L270-276).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 3];
        self.value(sol, &mut valsol);
        valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= 2.0 * tol * self.ray.abs()
    }
}

impl<'a> FunctionSetWithDerivatives for BlendSurfPointConstRadInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: Blend_SurfPointFuncInv::NbVariables (Blend_SurfPointFuncInv.cxx
        // L19-27) — returns 3; dispatches to the trait override below.
        BlendSurfPointFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendSurfPointConstRadInv::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendSurfPointConstRadInv::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfPointConstRadInv::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfPointConstRadInv::values(self, x, f, df)
    }
}

impl<'a> BlendSurfPointFuncInv for BlendSurfPointConstRadInv<'a> {
    /// OCCT Blend_SurfPointFuncInv::NbVariables — returns 3 (see the
    /// FunctionSetWithDerivatives note above; the rcad trait default of 2
    /// predates this anchor and is overridden here).
    fn nb_variables(&self) -> usize {
        3
    }

    fn nb_equations(&self) -> usize {
        BlendSurfPointConstRadInv::nb_equations(self)
    }

    fn set_point(&mut self, p: DVec3) {
        BlendSurfPointConstRadInv::set_point(self, p)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendSurfPointConstRadInv::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendSurfPointConstRadInv::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendSurfPointConstRadInv::is_solution(self, sol, tol)
    }
}
