//! OCCT BlendFunc_EvolRadInv (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_EvolRadInv.hxx (L27-81) + BlendFunc_EvolRadInv.cxx (whole file
//! L26-470).  `BRepBlend_EvolRadInv` is a pure typedef
//! (`typedef BlendFunc_EvolRadInv BRepBlend_EvolRadInv;`,
//! BRepBlend_EvolRadInv.hxx L21) so this single translation serves both
//! names — the SimulSurf (ChFi3d_FilBuilder.cxx L666) and PerformSurf
//! (L1650) EvolRad arms instantiate this class.
//!
//! Architecture mappings (mirroring [`super::brep_blend_func_consrad_c`]):
//! `class BlendFunc_EvolRadInv : public Blend_FuncInv` is expressed by
//! implementing the [`BlendFuncInv`] trait over the
//! `math_FunctionSetWithDerivatives` base; `occ::handle<Adaptor2d_Curve2d>
//! csurf` (null until Set(OnFirst, C)) maps to `Option<&'a Curve2d>`;
//! `occ::handle<Law_Function> fevol` maps to the shared
//! [`LawFunctionHandle`]; `math_Vector` / `math_Matrix` map to `[f64; 4]` /
//! `Vec<Vec<f64>>` (OCCT D(i, j) -> d[i - 1][j - 1], X(i) -> x[i - 1]).
//! Pending kernel dependency (marked GAP, plan 0.6):
//! Adaptor2d_Curve2d::Resolution.

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use crate::geomalgo::law::LawFunctionHandle;

use super::brep_blend_func::blend_func_compute_normal;
use super::brep_blend_func_evolrad::EPS;
use super::brep_blend_func_inv::BlendFuncInv;

/// OCCT BlendFunc_EvolRadInv — inverse of the variable-radius function used
/// to find a solution on a restriction of one of the surfaces
/// (BlendFunc_EvolRadInv.hxx L27).  The vector X is t, w, U, V.
pub struct BlendFuncEvolRadInv<'a> {
    // OCCT BlendFunc_EvolRadInv.hxx fields (L70-78).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) csurf: Option<&'a Curve2d>,
    pub(crate) fevol: LawFunctionHandle,
    pub(crate) sg1: f64,
    pub(crate) sg2: f64,
    pub(crate) choix: i32,
    pub(crate) first: bool,
}

use glam::DVec3;

impl<'a> BlendFuncEvolRadInv<'a> {
    /// OCCT BlendFunc_EvolRadInv(S1, S2, C, Law) (BlendFunc_EvolRadInv.cxx
    /// L26-35) — surf1(S1), surf2(S2), curv(C), fevol = Law.  The remaining
    /// members (sg1, sg2, choix, first) are uninitialized in OCCT; the rcad
    /// port zero-initializes them like the ConstRadInv constructor.
    pub fn new(
        s1: &'a Surface3,
        s2: &'a Surface3,
        c: &'a Curve3,
        law: LawFunctionHandle,
    ) -> Self {
        BlendFuncEvolRadInv {
            surf1: s1,
            surf2: s2,
            curv: c,
            csurf: None,
            fevol: law,
            sg1: 0.0,
            sg2: 0.0,
            choix: 0,
            first: false,
        }
    }

    /// OCCT Set(Choix) (BlendFunc_EvolRadInv.cxx L37-70) — inits the
    /// "quadrant" signs.
    pub fn set(&mut self, choix: i32) {
        self.choix = choix;
        match self.choix {
            1 | 2 => {
                self.sg1 = -1.0;
                self.sg2 = -1.0;
            }
            3 | 4 => {
                self.sg1 = 1.0;
                self.sg2 = -1.0;
            }
            5 | 6 => {
                self.sg1 = 1.0;
                self.sg2 = 1.0;
            }
            7 | 8 => {
                self.sg1 = -1.0;
                self.sg2 = 1.0;
            }
            _ => {
                self.sg1 = -1.0;
                self.sg2 = -1.0;
            }
        }
    }

    /// OCCT Set(OnFirst, C) (BlendFunc_EvolRadInv.cxx L72-76).
    pub fn set_curve_on_surface(&mut self, on_first: bool, c: &'a Curve2d) {
        self.first = on_first;
        self.csurf = Some(c);
    }

    /// OCCT NbEquations() (BlendFunc_EvolRadInv.cxx L78-81) — returns 4.
    pub fn nb_equations(&self) -> usize {
        4
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BlendFunc_EvolRadInv.cxx L83-97).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L85: Tolerance(1) = csurf->Resolution(Tol).
        // GAP (plan 0.6): rcad-kernel has no Adaptor2d_Curve2d::Resolution
        // equivalent yet; the call panics with the pending marker until the
        // kernel exposes it (same carrier as BlendFunc_ConstRadInv).
        tolerance[0] = super::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending();
        // OCCT L86: Tolerance(2) = curv->Resolution(Tol).
        tolerance[1] = self.curv.resolution(tol);
        if self.first {
            // OCCT L89-90: surf2->UResolution / VResolution.
            tolerance[2] = self.surf2.u_resolution(tol);
            tolerance[3] = self.surf2.v_resolution(tol);
        } else {
            // OCCT L94-95: surf1->UResolution / VResolution.
            tolerance[2] = self.surf1.u_resolution(tol);
            tolerance[3] = self.surf1.v_resolution(tol);
        }
    }

    /// OCCT GetBounds(InfBound, SupBound) (BlendFunc_EvolRadInv.cxx L99-144).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        let csurf = self.csurf.expect("csurf");
        // OCCT L101-104: the csurf / curv parameter ranges.
        inf_bound[0] = csurf.default_domain()[0]; // FirstParameter
        inf_bound[1] = self.curv.default_domain()[0];
        sup_bound[0] = csurf.default_domain()[1]; // LastParameter
        sup_bound[1] = self.curv.default_domain()[1];

        if self.first {
            // OCCT L108-111: the surf2 parameter ranges.
            inf_bound[2] = self.surf2.default_domain()[0]; // FirstUParameter
            inf_bound[3] = self.surf2.default_domain()[2]; // FirstVParameter
            sup_bound[2] = self.surf2.default_domain()[1]; // LastUParameter
            sup_bound[3] = self.surf2.default_domain()[3]; // LastVParameter
            // OCCT L112-117: the periodic-u range widening.
            if !is_infinite_value(inf_bound[2]) && !is_infinite_value(sup_bound[2]) {
                let range = sup_bound[2] - inf_bound[2];
                inf_bound[2] -= range;
                sup_bound[2] += range;
            }
            // OCCT L118-123: the periodic-v range widening.
            if !is_infinite_value(inf_bound[3]) && !is_infinite_value(sup_bound[3]) {
                let range = sup_bound[3] - inf_bound[3];
                inf_bound[3] -= range;
                sup_bound[3] += range;
            }
        } else {
            // OCCT L127-130: the surf1 parameter ranges.
            inf_bound[2] = self.surf1.default_domain()[0];
            inf_bound[3] = self.surf1.default_domain()[2];
            sup_bound[2] = self.surf1.default_domain()[1];
            sup_bound[3] = self.surf1.default_domain()[3];
            // OCCT L131-136.
            if !is_infinite_value(inf_bound[2]) && !is_infinite_value(sup_bound[2]) {
                let range = sup_bound[2] - inf_bound[2];
                inf_bound[2] -= range;
                sup_bound[2] += range;
            }
            // OCCT L137-142.
            if !is_infinite_value(inf_bound[3]) && !is_infinite_value(sup_bound[3]) {
                let range = sup_bound[3] - inf_bound[3];
                inf_bound[3] -= range;
                sup_bound[3] += range;
            }
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_EvolRadInv.cxx L146-152).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 4];
        self.value(sol, &mut valsol);
        valsol[0].abs() <= tol
            && valsol[1] * valsol[1] + valsol[2] * valsol[2] + valsol[3] * valsol[3] <= tol * tol
    }

    /// OCCT Value(X, F) (BlendFunc_EvolRadInv.cxx L154-238).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L156: const double ray = fevol->Value(X(2));
        let ray = self.fevol.borrow_mut().value(x[1]);

        // OCCT L158-160: gp_Pnt ptcur; gp_Vec d1cur; curv->D1(X(2), ptcur, d1cur);
        let ptcur = self.curv.point_at(x[1]);
        let d1cur = self.curv.derivative_at(x[1]);

        // OCCT L162-163: nplan = d1cur.Normalized(); theD = -(nplan . ptcur).
        let nplan = d1cur.normalize_or_zero();
        let the_d = -nplan.dot(ptcur);

        // OCCT L165: const gp_Pnt2d pt2d(csurf->Value(X(1)));
        let csurf = self.csurf.expect("csurf");
        let pt2d = csurf.point_at(x[0]);

        // OCCT L167-178: the surface evaluations.
        let pts1: DVec3;
        let pts2: DVec3;
        let d1u1: DVec3;
        let d1v1: DVec3;
        let d1u2: DVec3;
        let d1v2: DVec3;
        if self.first {
            // OCCT L171: surf1->D1(pt2d.X(), pt2d.Y(), pts1, d1u1, d1v1);
            let (p, du, dv) = self.surf1.derivatives(pt2d.x, pt2d.y);
            pts1 = p;
            d1u1 = du;
            d1v1 = dv;
            // OCCT L172: surf2->D1(X(3), X(4), pts2, d1u2, d1v2);
            let (p, du, dv) = self.surf2.derivatives(x[2], x[3]);
            pts2 = p;
            d1u2 = du;
            d1v2 = dv;
        } else {
            // OCCT L176: surf1->D1(X(3), X(4), pts1, d1u1, d1v1);
            let (p, du, dv) = self.surf1.derivatives(x[2], x[3]);
            pts1 = p;
            d1u1 = du;
            d1v1 = dv;
            // OCCT L177: surf2->D1(pt2d.X(), pt2d.Y(), pts2, d1u2, d1v2);
            let (p, du, dv) = self.surf2.derivatives(pt2d.x, pt2d.y);
            pts2 = p;
            d1u2 = du;
            d1v2 = dv;
        }

        // OCCT L180-183: F(1) — the plane equation.
        f[0] = (nplan.x * (pts1.x + pts2.x) + nplan.y * (pts1.y + pts2.y)
            + nplan.z * (pts1.z + pts2.z))
            / 2.0
            + the_d;

        // OCCT L185-197: ns1 with the degenerate-normal fallback.
        let mut ns1 = d1u1.cross(d1v1);
        if ns1.length() < EPS {
            if self.first {
                // OCCT L190: BlendFunc::ComputeNormal(surf1, pt2d, ns1);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (pt2d.x, pt2d.y), &mut n);
                ns1 = n;
            } else {
                // OCCT L194-195: gp_Pnt2d P(X(3), X(4));
                //                BlendFunc::ComputeNormal(surf1, P, ns1);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (x[2], x[3]), &mut n);
                ns1 = n;
            }
        }

        // OCCT L199-211: ns2 with the degenerate-normal fallback (the OCCT
        // `.XYZ()` on L199 is a gp_XYZ round-trip of the same vector).
        let mut ns2 = d1u2.cross(d1v2);
        if ns2.length() < EPS {
            if !self.first {
                // OCCT L204: BlendFunc::ComputeNormal(surf2, pt2d, ns2);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (pt2d.x, pt2d.y), &mut n);
                ns2 = n;
            } else {
                // OCCT L208-209: gp_Pnt2d P(X(3), X(4));
                //                BlendFunc::ComputeNormal(surf2, P, ns2);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (x[2], x[3]), &mut n);
                ns2 = n;
            }
        }

        // OCCT L213-216: the cross norms.
        let ncrossns1 = nplan.cross(ns1);
        let ncrossns2 = nplan.cross(ns2);
        let mut norm1 = ncrossns1.length();
        let mut norm2 = ncrossns2.length();

        // OCCT L218-225.
        if norm1 < EPS {
            norm1 = 1.0;
        }
        if norm2 < EPS {
            norm2 = 1.0;
        }

        // OCCT L227-231: the projection of the normals onto the plane.
        let ndotns1 = nplan.dot(ns1);
        let ndotns2 = nplan.dot(ns2);
        ns1 = (ndotns1 / norm1) * nplan + (-1.0 / norm1) * ns1;
        ns2 = (ndotns2 / norm2) * nplan + (-1.0 / norm2) * ns2;
        // OCCT L232:
        // resul.SetLinearForm(sg1 * ray, ns1, -1., pts2.XYZ(), -sg2 * ray, ns2, pts1.XYZ());
        let resul = self.sg1 * ray * ns1 + (-1.0) * pts2 + (-self.sg2 * ray) * ns2 + pts1;
        // OCCT L233-235.
        f[1] = resul.x;
        f[2] = resul.y;
        f[3] = resul.z;

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_EvolRadInv.cxx L240-244) —
    /// `math_Vector F(1, 4); return Values(X, F, D);`.
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        let mut f = [0.0f64; 4];
        self.values(x, &mut f, d)
    }

    /// OCCT Values(X, F, D) (BlendFunc_EvolRadInv.cxx L246-470).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L248-249: double ray, dray; fevol->D1(X(2), ray, dray);
        let mut ray = 0.0f64;
        let mut dray = 0.0f64;
        self.fevol.borrow_mut().d1(x[1], &mut ray, &mut dray);

        // OCCT L251-253: curv->D2(X(2), ptcur, d1cur, d2cur);
        let ptcur = self.curv.point_at(x[1]);
        let d1cur = self.curv.derivative_at(x[1]);
        let d2cur = self.curv.derivative2_at(x[1]);

        // OCCT L255-257.
        let normtgcur = d1cur.length();
        let nplan = d1cur.normalize_or_zero();
        let the_d = -nplan.dot(ptcur);

        // OCCT L259-261:
        // dnplan.SetLinearForm(-nplan.Dot(d2cur), nplan, d2cur); dnplan /= normtgcur;
        let mut dnplan = (-nplan.dot(d2cur)) * nplan + d2cur;
        dnplan /= normtgcur;

        // OCCT L263-265: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let v2d = csurf.derivative_at(x[0]);

        // OCCT L267-294: the surface D2 evaluations and the D(1, *) row.
        let pts1: DVec3;
        let pts2: DVec3;
        let d1u1: DVec3;
        let d1v1: DVec3;
        let d1u2: DVec3;
        let d1v2: DVec3;
        let d2u1: DVec3;
        let d2v1: DVec3;
        let d2uv1: DVec3;
        let d2u2: DVec3;
        let d2v2: DVec3;
        let d2uv2: DVec3;
        if self.first {
            // OCCT L273: surf1->D2(p2d.X(), p2d.Y(), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1);
            let (p1, du1, dv1, d2u_1, d2uv_1, d2v_1) = self.surf1.derivatives2(p2d.x, p2d.y);
            pts1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = d2u_1;
            d2uv1 = d2uv_1;
            d2v1 = d2v_1;
            // OCCT L274: surf2->D2(X(3), X(4), pts2, d1u2, d1v2, d2u2, d2v2, d2uv2);
            let (p2, du2, dv2, d2u_2, d2uv_2, d2v_2) = self.surf2.derivatives2(x[2], x[3]);
            pts2 = p2;
            d1u2 = du2;
            d1v2 = dv2;
            d2u2 = d2u_2;
            d2uv2 = d2uv_2;
            d2v2 = d2v_2;
            // OCCT L276-281.
            let temp = v2d.x * d1u1 + v2d.y * d1v1;
            d[0][0] = nplan.dot(temp) / 2.0;
            let temp = 0.5 * (pts1 + pts2) - ptcur;
            d[0][1] = dnplan.dot(temp) - normtgcur;
            d[0][2] = nplan.dot(d1u2) / 2.0;
            d[0][3] = nplan.dot(d1v2) / 2.0;
        } else {
            // OCCT L285: surf1->D2(X(3), X(4), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1);
            let (p1, du1, dv1, d2u_1, d2uv_1, d2v_1) = self.surf1.derivatives2(x[2], x[3]);
            pts1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = d2u_1;
            d2uv1 = d2uv_1;
            d2v1 = d2v_1;
            // OCCT L286: surf2->D2(p2d.X(), p2d.Y(), pts2, d1u2, d1v2, d2u2, d2v2, d2uv2);
            let (p2, du2, dv2, d2u_2, d2uv_2, d2v_2) = self.surf2.derivatives2(p2d.x, p2d.y);
            pts2 = p2;
            d1u2 = du2;
            d1v2 = dv2;
            d2u2 = d2u_2;
            d2uv2 = d2uv_2;
            d2v2 = d2v_2;
            // OCCT L288-293.
            let temp = v2d.x * d1u2 + v2d.y * d1v2;
            d[0][0] = nplan.dot(temp) / 2.0;
            let temp = 0.5 * (pts1 + pts2) - ptcur;
            d[0][1] = dnplan.dot(temp) - normtgcur;
            d[0][2] = nplan.dot(d1u1) / 2.0;
            d[0][3] = nplan.dot(d1v1) / 2.0;
        }

        // OCCT L296-299: F(1) — the plane equation.
        f[0] = (nplan.x * (pts1.x + pts2.x) + nplan.y * (pts1.y + pts2.y)
            + nplan.z * (pts1.z + pts2.z))
            / 2.0
            + the_d;

        // OCCT L301-313: ns1 with the degenerate-normal fallback.
        let mut ns1 = d1u1.cross(d1v1);
        if ns1.length() < EPS {
            if self.first {
                // OCCT L306: BlendFunc::ComputeNormal(surf1, p2d, ns1);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (p2d.x, p2d.y), &mut n);
                ns1 = n;
            } else {
                // OCCT L310-311: gp_Pnt2d P(X(3), X(4)); ComputeNormal(surf1, P, ns1);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (x[2], x[3]), &mut n);
                ns1 = n;
            }
        }

        // OCCT L315-327: ns2 with the degenerate-normal fallback.
        let mut ns2 = d1u2.cross(d1v2);
        if ns2.length() < EPS {
            if !self.first {
                // OCCT L320: BlendFunc::ComputeNormal(surf2, p2d, ns2);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (p2d.x, p2d.y), &mut n);
                ns2 = n;
            } else {
                // OCCT L324-325: gp_Pnt2d P(X(3), X(4)); ComputeNormal(surf2, P, ns2);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (x[2], x[3]), &mut n);
                ns2 = n;
            }
        }

        // OCCT L329-340: the cross norms and their Eps guards.
        let ncrossns1 = nplan.cross(ns1);
        let ncrossns2 = nplan.cross(ns2);
        let mut norm1 = ncrossns1.length();
        let mut norm2 = ncrossns2.length();
        if norm1 < EPS {
            norm1 = 1.0;
        }
        if norm2 < EPS {
            norm2 = 1.0;
        }

        // OCCT L342-350: the function value rows F(2)..F(4).
        let mut resul1: DVec3;
        let mut resul2: DVec3;
        let mut grosterme: f64;

        let ndotns1 = nplan.dot(ns1);
        let ndotns2 = nplan.dot(ns2);
        // OCCT L347: temp.SetLinearForm(ndotns1 / norm1, nplan, -1. / norm1, ns1);
        let temp = (ndotns1 / norm1) * nplan + (-1.0 / norm1) * ns1;
        // OCCT L348: resul1.SetLinearForm(sg1 * ray, temp, gp_Vec(pts2, pts1));
        resul1 = self.sg1 * ray * temp + (pts1 - pts2);
        // OCCT L349: temp.SetLinearForm(ndotns2 / norm2, nplan, -1. / norm2, ns2);
        let temp = (ndotns2 / norm2) * nplan + (-1.0 / norm2) * ns2;
        // OCCT L350: resul1.Subtract(sg2 * ray * temp);
        resul1 -= self.sg2 * ray * temp;

        // OCCT L352-354.
        f[1] = resul1.x;
        f[2] = resul1.y;
        f[3] = resul1.z;

        // OCCT L356-366: Derivee par rapport a u1.
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        grosterme = ncrossns1.dot(nplan.cross(temp)) / norm1 / norm1;
        resul1 = (-self.sg1 * ray / norm1 * (grosterme * ndotns1 - nplan.dot(temp))) * nplan
            + (self.sg1 * ray * grosterme / norm1) * ns1
            + (-self.sg1 * ray / norm1) * temp
            + d1u1;

        // OCCT L368-378: Derivee par rapport a v1.
        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        grosterme = ncrossns1.dot(nplan.cross(temp)) / norm1 / norm1;
        resul2 = (-self.sg1 * ray / norm1 * (grosterme * ndotns1 - nplan.dot(temp))) * nplan
            + (self.sg1 * ray * grosterme / norm1) * ns1
            + (-self.sg1 * ray / norm1) * temp
            + d1v1;

        // OCCT L380-395.
        if self.first {
            d[1][0] = v2d.x * resul1.x + v2d.y * resul2.x;
            d[2][0] = v2d.x * resul1.y + v2d.y * resul2.y;
            d[3][0] = v2d.x * resul1.z + v2d.y * resul2.z;
        } else {
            d[1][2] = resul1.x;
            d[2][2] = resul1.y;
            d[3][2] = resul1.z;

            d[1][3] = resul2.x;
            d[2][3] = resul2.y;
            d[3][3] = resul2.z;
        }

        // OCCT L397-405: derivee par rapport a w (parametre sur ligne guide).
        grosterme = ncrossns1.dot(dnplan.cross(ns1)) / norm1 / norm1;
        resul1 = (-self.sg1 / norm1 * (grosterme * ndotns1 - dnplan.dot(ns1))) * nplan
            + (self.sg1 * ndotns1 / norm1) * dnplan
            + (self.sg1 * grosterme / norm1) * ns1;

        // OCCT L407-413.
        grosterme = ncrossns2.dot(dnplan.cross(ns2)) / norm2 / norm2;
        resul2 = (self.sg2 / norm2 * (grosterme * ndotns2 - dnplan.dot(ns2))) * nplan
            + (-self.sg2 * ndotns2 / norm2) * dnplan
            + (-self.sg2 * grosterme / norm2) * ns2;

        // OCCT L415-417.
        d[1][1] = ray * (resul1.x + resul2.x);
        d[2][1] = ray * (resul1.y + resul2.y);
        d[3][1] = ray * (resul1.z + resul2.z);

        // OCCT L419-424: the dray term direction.
        let temp = (self.sg1 * ndotns1 / norm1 - self.sg2 * ndotns2 / norm2) * nplan
            + (-self.sg1 / norm1) * ns1
            + (self.sg2 / norm2) * ns2;

        // OCCT L426-428.
        d[1][1] += dray * temp.x;
        d[2][1] += dray * temp.y;
        d[3][1] += dray * temp.z;

        // OCCT L430-439: Derivee par rapport a u2.
        let temp = d2u2.cross(d1v2) + d1u2.cross(d2uv2);
        grosterme = ncrossns2.dot(nplan.cross(temp)) / norm2 / norm2;
        resul1 = (self.sg2 * ray / norm2 * (grosterme * ndotns2 - nplan.dot(temp))) * nplan
            + (-self.sg2 * ray * grosterme / norm2) * ns2
            + (self.sg2 * ray / norm2) * temp;
        // OCCT L439: resul1.Subtract(d1u2);
        resul1 -= d1u2;

        // OCCT L441-450: Derivee par rapport a v2.
        let temp = d2uv2.cross(d1v2) + d1u2.cross(d2v2);
        grosterme = ncrossns2.dot(nplan.cross(temp)) / norm2 / norm2;
        resul2 = (self.sg2 * ray / norm2 * (grosterme * ndotns2 - nplan.dot(temp))) * nplan
            + (-self.sg2 * ray * grosterme / norm2) * ns2
            + (self.sg2 * ray / norm2) * temp;
        // OCCT L450: resul2.Subtract(d1v2);
        resul2 -= d1v2;

        // OCCT L452-467.
        if !self.first {
            d[1][0] = v2d.x * resul1.x + v2d.y * resul2.x;
            d[2][0] = v2d.x * resul1.y + v2d.y * resul2.y;
            d[3][0] = v2d.x * resul1.z + v2d.y * resul2.z;
        } else {
            d[1][2] = resul1.x;
            d[2][2] = resul1.y;
            d[3][2] = resul1.z;

            d[1][3] = resul2.x;
            d[2][3] = resul2.y;
            d[3][3] = resul2.z;
        }

        true
    }
}

impl<'a> FunctionSetWithDerivatives for BlendFuncEvolRadInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: inherited from Blend_FuncInv::NbVariables (returns 4).
        BlendFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendFuncEvolRadInv::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendFuncEvolRadInv::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncEvolRadInv::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncEvolRadInv::values(self, x, f, df)
    }
}

impl<'a> BlendFuncInv for BlendFuncEvolRadInv<'a> {
    fn set_curve_on_surface(&mut self, on_first: bool, c_on_surf: &Curve2d) {
        // OCCT EvolRadInv.cxx L72-76; the rcad port stores the curve
        // reference (OCCT copies the handle).  SAFETY: the caller owns the
        // Curve2d for the lifetime 'a of this function object — same
        // invariant as the OCCT handle (pattern of BlendFunc_ConstRadInv).
        let c: &'a Curve2d = unsafe { &*(c_on_surf as *const Curve2d) };
        BlendFuncEvolRadInv::set_curve_on_surface(self, on_first, c)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendFuncEvolRadInv::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendFuncEvolRadInv::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendFuncEvolRadInv::is_solution(self, sol, tol)
    }
}
