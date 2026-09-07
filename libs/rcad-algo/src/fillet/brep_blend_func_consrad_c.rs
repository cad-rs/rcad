//! OCCT BlendFunc_ConstRadInv (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_ConstRadInv.hxx (L26-78) + BlendFunc_ConstRadInv.cxx (whole
//! file L25-646).  The struct is defined in
//! [`super::brep_blend_func_consrad`]; this file carries its method bodies
//! and trait wiring.
//!
//! Architecture mappings: `class BlendFunc_ConstRadInv : public
//! Blend_FuncInv` is expressed by implementing the [`BlendFuncInv`] trait
//! and the `math_FunctionSetWithDerivatives` base; `math_Vector` /
//! `math_Matrix` map to `[f64; 4]` / `Vec<Vec<f64>>` (OCCT D(i, j) ->
//! d[i - 1][j - 1]).  Pending kernel dependency (marked GAP, plan 0.6):
//! Adaptor2d_Curve2d::Resolution.

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::brep_blend_func::blend_func_compute_normal;
use super::brep_blend_func_consrad::EPS;
use super::brep_blend_func_consrad::BlendFuncConstRadInv;
use super::brep_blend_func_inv::BlendFuncInv;

impl<'a> BlendFuncConstRadInv<'a> {
    /// OCCT BlendFunc_ConstRadInv(S1, S2, C) (BlendFunc_ConstRadInv.cxx
    /// L25-36).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncConstRadInv {
            surf1: s1,
            surf2: s2,
            curv: c,
            csurf: None,
            ray1: 0.0,
            ray2: 0.0,
            choix: 0,
            first: false,
        }
    }

    /// OCCT Set(R, Choix) (BlendFunc_ConstRadInv.cxx L38-71) — inits the
    /// value of radius, and the "quadrant".
    pub fn set(&mut self, r: f64, choix: i32) {
        self.choix = choix;
        match self.choix {
            1 | 2 => {
                self.ray1 = -r;
                self.ray2 = -r;
            }
            3 | 4 => {
                self.ray1 = r;
                self.ray2 = -r;
            }
            5 | 6 => {
                self.ray1 = r;
                self.ray2 = r;
            }
            7 | 8 => {
                self.ray1 = -r;
                self.ray2 = r;
            }
            _ => {
                self.ray1 = -r;
                self.ray2 = -r;
            }
        }
    }

    /// OCCT Set(OnFirst, C) (BlendFunc_ConstRadInv.cxx L73-77).
    pub fn set_curve_on_surface(&mut self, on_first: bool, c: &'a Curve2d) {
        self.first = on_first;
        self.csurf = Some(c);
    }

    /// OCCT NbEquations() (BlendFunc_ConstRadInv.cxx L79-82) — returns 4.
    pub fn nb_equations(&self) -> usize {
        4
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BlendFunc_ConstRadInv.cxx L84-98).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L86: Tolerance(1) = csurf->Resolution(Tol).
        // GAP (plan 0.6): rcad-kernel has no Adaptor2d_Curve2d::Resolution
        // equivalent yet; the call panics with the pending marker until the
        // kernel exposes it.
        tolerance[0] = super::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending();
        // OCCT L87: Tolerance(2) = curv->Resolution(Tol).
        tolerance[1] = self.curv.resolution(tol);
        if self.first {
            tolerance[2] = self.surf2.u_resolution(tol);
            tolerance[3] = self.surf2.v_resolution(tol);
        } else {
            tolerance[2] = self.surf1.u_resolution(tol);
            tolerance[3] = self.surf1.v_resolution(tol);
        }
    }

    /// OCCT GetBounds(InfBound, SupBound) (BlendFunc_ConstRadInv.cxx
    /// L100-145).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        let csurf = self.csurf.expect("csurf");
        inf_bound[0] = csurf.default_domain()[0]; // FirstParameter
        inf_bound[1] = self.curv.default_domain()[0];
        sup_bound[0] = csurf.default_domain()[1]; // LastParameter
        sup_bound[1] = self.curv.default_domain()[1];

        if self.first {
            inf_bound[2] = self.surf2.default_domain()[0];
            inf_bound[3] = self.surf2.default_domain()[2];
            sup_bound[2] = self.surf2.default_domain()[1];
            sup_bound[3] = self.surf2.default_domain()[3];
            if !is_infinite_value(inf_bound[2]) && !is_infinite_value(sup_bound[2]) {
                let range = sup_bound[2] - inf_bound[2];
                inf_bound[2] -= range;
                sup_bound[2] += range;
            }
            if !is_infinite_value(inf_bound[3]) && !is_infinite_value(sup_bound[3]) {
                let range = sup_bound[3] - inf_bound[3];
                inf_bound[3] -= range;
                sup_bound[3] += range;
            }
        } else {
            inf_bound[2] = self.surf1.default_domain()[0];
            inf_bound[3] = self.surf1.default_domain()[2];
            sup_bound[2] = self.surf1.default_domain()[1];
            sup_bound[3] = self.surf1.default_domain()[3];
            if !is_infinite_value(inf_bound[2]) && !is_infinite_value(sup_bound[2]) {
                let range = sup_bound[2] - inf_bound[2];
                inf_bound[2] -= range;
                sup_bound[2] += range;
            }
            if !is_infinite_value(inf_bound[3]) && !is_infinite_value(sup_bound[3]) {
                let range = sup_bound[3] - inf_bound[3];
                inf_bound[3] -= range;
                sup_bound[3] += range;
            }
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_ConstRadInv.cxx L147-153).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 4];
        self.value(sol, &mut valsol);
        valsol[0].abs() <= tol
            && valsol[1] * valsol[1] + valsol[2] * valsol[2] + valsol[3] * valsol[3] <= tol * tol
    }

    /// OCCT Value(X, F) (BlendFunc_ConstRadInv.cxx L155-232).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: curv->D1(X(2), ptcur, d1cur);
        let ptcur = self.curv.point_at(x[1]);
        let d1cur = self.curv.derivative_at(x[1]);

        let nplan = d1cur.normalize_or_zero();
        let the_d = -nplan.dot(ptcur);

        // OCCT: const gp_Pnt2d pt2d(csurf->Value(X(1)));
        let csurf = self.csurf.expect("csurf");
        let pt2d = csurf.point_at(x[0]);

        let pts1: DVec3;
        let pts2: DVec3;
        let d1u1: DVec3;
        let d1v1: DVec3;
        let d1u2: DVec3;
        let d1v2: DVec3;
        if self.first {
            // OCCT: surf1->D1(pt2d.X(), pt2d.Y(), pts1, d1u1, d1v1);
            let (p, du, dv) = self.surf1.derivatives(pt2d.x, pt2d.y);
            pts1 = p;
            d1u1 = du;
            d1v1 = dv;
            // OCCT: surf2->D1(X(3), X(4), pts2, d1u2, d1v2);
            let (p, du, dv) = self.surf2.derivatives(x[2], x[3]);
            pts2 = p;
            d1u2 = du;
            d1v2 = dv;
        } else {
            // OCCT: surf1->D1(X(3), X(4), pts1, d1u1, d1v1);
            let (p, du, dv) = self.surf1.derivatives(x[2], x[3]);
            pts1 = p;
            d1u1 = du;
            d1v1 = dv;
            // OCCT: surf2->D1(pt2d.X(), pt2d.Y(), pts2, d1u2, d1v2);
            let (p, du, dv) = self.surf2.derivatives(pt2d.x, pt2d.y);
            pts2 = p;
            d1u2 = du;
            d1v2 = dv;
        }

        f[0] = (nplan.x * (pts1.x + pts2.x) + nplan.y * (pts1.y + pts2.y)
            + nplan.z * (pts1.z + pts2.z))
            / 2.0
            + the_d;

        let mut ns1 = d1u1.cross(d1v1);
        if ns1.length() < EPS {
            if self.first {
                // OCCT: BlendFunc::ComputeNormal(surf1, pt2d, ns1);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (pt2d.x, pt2d.y), &mut n);
                ns1 = n;
            } else {
                // OCCT: gp_Pnt2d P(X(3), X(4)); BlendFunc::ComputeNormal(surf1, P, ns1);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (x[2], x[3]), &mut n);
                ns1 = n;
            }
        }

        let mut ns2 = d1u2.cross(d1v2);
        if ns2.length() < EPS {
            if !self.first {
                // OCCT: BlendFunc::ComputeNormal(surf2, pt2d, ns2);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (pt2d.x, pt2d.y), &mut n);
                ns2 = n;
            } else {
                // OCCT: gp_Pnt2d P(X(3), X(4)); BlendFunc::ComputeNormal(surf2, P, ns2);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (x[2], x[3]), &mut n);
                ns2 = n;
            }
        }

        let mut norm1 = nplan.cross(ns1).length();
        let mut norm2 = nplan.cross(ns2).length();
        if norm1 < EPS {
            norm1 = 1.0;
        }
        if norm2 < EPS {
            norm2 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }

        // OCCT: ns1.SetLinearForm(nplan.Dot(ns1) / norm1, nplan, -1. / norm1, ns1);
        ns1 = (nplan.dot(ns1) / norm1) * nplan + (-1.0 / norm1) * ns1;
        ns2 = (nplan.dot(ns2) / norm2) * nplan + (-1.0 / norm2) * ns2;
        // OCCT: resul.SetLinearForm(ray1, ns1, -1., pts2.XYZ(), -ray2, ns2, pts1.XYZ());
        let resul = self.ray1 * ns1 + (-1.0) * pts2 + (-self.ray2) * ns2 + pts1;
        f[1] = resul.x;
        f[2] = resul.y;
        f[3] = resul.z;

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_ConstRadInv.cxx L234-432).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        let ns1: DVec3;
        let ns2: DVec3;
        let nplan: DVec3;
        let mut dnplan: DVec3;
        let pts1: DVec3;
        let pts2: DVec3;
        let p2d: DVec2;
        let v2d: DVec2;
        let d1u1: DVec3;
        let d1v1: DVec3;
        let d2u1: DVec3;
        let d2v1: DVec3;
        let d2uv1: DVec3;
        let d1u2: DVec3;
        let d1v2: DVec3;
        let d2u2: DVec3;
        let d2v2: DVec3;
        let d2uv2: DVec3;
        let mut norm1: f64;
        let mut norm2: f64;
        let ndotns1: f64;
        let ndotns2: f64;
        let normtgcur: f64;
        let mut grosterme: f64;
        let the_d: f64;

        // OCCT: curv->D2(X(2), ptcur, d1cur, d2cur);
        let ptcur = self.curv.point_at(x[1]);
        let d1cur = self.curv.derivative_at(x[1]);
        let d2cur = self.curv.derivative2_at(x[1]);
        normtgcur = d1cur.length();
        nplan = d1cur.normalize_or_zero();
        the_d = -nplan.dot(ptcur);

        // OCCT: dnplan.SetLinearForm(theD, nplan, d2cur); dnplan /= normtgcur;
        dnplan = the_d * nplan + d2cur;
        dnplan /= normtgcur;

        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        p2d = csurf.point_at(x[0]);
        v2d = csurf.derivative_at(x[0]);

        if self.first {
            // OCCT: surf1->D2(p2d.X(), p2d.Y(), pts1, ...); surf2->D2(X(3), X(4), pts2, ...);
            let (p1, du1, dv1, d2u_1, d2uv_1, d2v_1) = self.surf1.derivatives2(p2d.x, p2d.y);
            pts1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = d2u_1;
            d2uv1 = d2uv_1;
            d2v1 = d2v_1;
            let (p2, du2, dv2, d2u_2, d2uv_2, d2v_2) = self.surf2.derivatives2(x[2], x[3]);
            pts2 = p2;
            d1u2 = du2;
            d1v2 = dv2;
            d2u2 = d2u_2;
            d2uv2 = d2uv_2;
            d2v2 = d2v_2;
            // OCCT: temp.SetLinearForm(v2d.X(), d1u1, v2d.Y(), d1v1);
            let temp = v2d.x * d1u1 + v2d.y * d1v1;
            d[0][0] = nplan.dot(temp) / 2.0;
            // OCCT: temp.SetXYZ(0.5 * (pts1.XYZ() + pts2.XYZ()) - ptcur.XYZ());
            let temp = 0.5 * (pts1 + pts2) - ptcur;
            d[0][1] = dnplan.dot(temp) - normtgcur;
            d[0][2] = nplan.dot(d1u2) / 2.0;
            d[0][3] = nplan.dot(d1v2) / 2.0;
        } else {
            // OCCT: surf1->D2(X(3), X(4), pts1, ...); surf2->D2(p2d.X(), p2d.Y(), pts2, ...);
            let (p1, du1, dv1, d2u_1, d2uv_1, d2v_1) = self.surf1.derivatives2(x[2], x[3]);
            pts1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = d2u_1;
            d2uv1 = d2uv_1;
            d2v1 = d2v_1;
            let (p2, du2, dv2, d2u_2, d2uv_2, d2v_2) = self.surf2.derivatives2(p2d.x, p2d.y);
            pts2 = p2;
            d1u2 = du2;
            d1v2 = dv2;
            d2u2 = d2u_2;
            d2uv2 = d2uv_2;
            d2v2 = d2v_2;
            // OCCT: temp.SetLinearForm(v2d.X(), d1u2, v2d.Y(), d1v2);
            let temp = v2d.x * d1u2 + v2d.y * d1v2;
            d[0][0] = nplan.dot(temp) / 2.0;
            let temp = 0.5 * (pts1 + pts2) - ptcur;
            d[0][1] = dnplan.dot(temp) - normtgcur;
            d[0][2] = nplan.dot(d1u1) / 2.0;
            d[0][3] = nplan.dot(d1v1) / 2.0;
        }

        let mut ns1v = d1u1.cross(d1v1);
        if ns1v.length() < EPS {
            if self.first {
                // OCCT: BlendFunc::ComputeNormal(surf1, p2d, ns1);
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (p2d.x, p2d.y), &mut n);
                ns1v = n;
            } else {
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (x[2], x[3]), &mut n);
                ns1v = n;
            }
        }

        let mut ns2v = d1u2.cross(d1v2);
        if ns2v.length() < EPS {
            if !self.first {
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (p2d.x, p2d.y), &mut n);
                ns2v = n;
            } else {
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (x[2], x[3]), &mut n);
                ns2v = n;
            }
        }

        let ncrossns1 = nplan.cross(ns1v);
        let ncrossns2 = nplan.cross(ns2v);
        norm1 = ncrossns1.length();
        norm2 = ncrossns2.length();
        if norm1 < EPS {
            norm1 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }
        if norm2 < EPS {
            norm2 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }

        ndotns1 = nplan.dot(ns1v);
        ndotns2 = nplan.dot(ns2v);

        ns1 = ns1v;
        ns2 = ns2v;
        let _ = (ns1, ns2);

        // Derived compared to u1

        // OCCT: temp = d2u1.Crossed(d1v1).Added(d1u1.Crossed(d2uv1));
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        grosterme = ncrossns1.dot(nplan.cross(temp)) / norm1 / norm1;
        // OCCT: resul1.SetLinearForm(-ray1 / norm1 * (grosterme * ndotns1 - nplan.Dot(temp)),
        //                           nplan, ray1 * grosterme / norm1, ns1,
        //                           -ray1 / norm1, temp, d1u1);
        let resul1 = (-self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(temp))) * nplan
            + (self.ray1 * grosterme / norm1) * ns1v
            + (-self.ray1 / norm1) * temp
            + d1u1;

        // Derived compared to v1

        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        grosterme = ncrossns1.dot(nplan.cross(temp)) / norm1 / norm1;
        let resul2 = (-self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(temp))) * nplan
            + (self.ray1 * grosterme / norm1) * ns1v
            + (-self.ray1 / norm1) * temp
            + d1v1;

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

        // derived compared to w (parameter on guideline)
        // It is assumed that the radius is constant

        grosterme = ncrossns1.dot(dnplan.cross(ns1v)) / norm1 / norm1;
        let resul1 = (-self.ray1 / norm1 * (grosterme * ndotns1 - dnplan.dot(ns1v))) * nplan
            + (self.ray1 * ndotns1 / norm1) * dnplan
            + (self.ray1 * grosterme / norm1) * ns1v;

        grosterme = ncrossns2.dot(dnplan.cross(ns2v)) / norm2 / norm2;
        let resul2 = (self.ray2 / norm2 * (grosterme * ndotns2 - dnplan.dot(ns2v))) * nplan
            + (-self.ray2 * ndotns2 / norm2) * dnplan
            + (-self.ray2 * grosterme / norm2) * ns2v;

        d[1][1] = resul1.x + resul2.x;
        d[2][1] = resul1.y + resul2.y;
        d[3][1] = resul1.z + resul2.z;

        // Derived compared to u2
        let mut temp = d2u2.cross(d1v2) + d1u2.cross(d2uv2);
        grosterme = ncrossns2.dot(nplan.cross(temp)) / norm2 / norm2;
        let mut resul1 = (self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(temp))) * nplan
            + (-self.ray2 * grosterme / norm2) * ns2v
            + (self.ray2 / norm2) * temp;
        // OCCT: resul1.Subtract(d1u2);
        resul1 -= d1u2;

        // Derived compared to v2
        temp = d2uv2.cross(d1v2) + d1u2.cross(d2v2);
        grosterme = ncrossns2.dot(nplan.cross(temp)) / norm2 / norm2;
        let mut resul2 = (self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(temp))) * nplan
            + (-self.ray2 * grosterme / norm2) * ns2v
            + (self.ray2 / norm2) * temp;
        // OCCT: resul2.Subtract(d1v2);
        resul2 -= d1v2;

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

    /// OCCT Values(X, F, D) (BlendFunc_ConstRadInv.cxx L434-646) — the
    /// combined evaluation; the OCCT body recomputes everything with the
    /// same formulas as Value + Derivatives (the dnplan differs by its
    /// SetLinearForm(-nplan.Dot(d2cur), nplan, d2cur) first coefficient).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        let nplan: DVec3;
        let mut dnplan: DVec3;
        let mut norm1: f64;
        let mut norm2: f64;
        let ndotns1: f64;
        let ndotns2: f64;
        let normtgcur: f64;
        let mut grosterme: f64;
        let the_d: f64;

        // OCCT: curv->D2(X(2), ptcur, d1cur, d2cur);
        let ptcur = self.curv.point_at(x[1]);
        let d1cur = self.curv.derivative_at(x[1]);
        let d2cur = self.curv.derivative2_at(x[1]);
        normtgcur = d1cur.length();
        nplan = d1cur.normalize_or_zero();
        the_d = -nplan.dot(ptcur);

        // NOTE: the Values body computes dnplan with the coefficient
        // -nplan.Dot(d2cur) (L449), unlike Derivatives which uses theD
        // (L251) — translated literally.
        dnplan = (-nplan.dot(d2cur)) * nplan + d2cur;
        dnplan /= normtgcur;

        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let v2d = csurf.derivative_at(x[0]);

        let (pts1, d1u1, d1v1, d2u1, d2uv1, d2v1, pts2, d1u2, d1v2, d2u2, d2uv2, d2v2);
        if self.first {
            let (p1, du1, dv1, d2u_1, d2uv_1, d2v_1) = self.surf1.derivatives2(p2d.x, p2d.y);
            pts1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = d2u_1;
            d2uv1 = d2uv_1;
            d2v1 = d2v_1;
            let (p2, du2, dv2, d2u_2, d2uv_2, d2v_2) = self.surf2.derivatives2(x[2], x[3]);
            pts2 = p2;
            d1u2 = du2;
            d1v2 = dv2;
            d2u2 = d2u_2;
            d2uv2 = d2uv_2;
            d2v2 = d2v_2;
            let temp = v2d.x * d1u1 + v2d.y * d1v1;
            d[0][0] = nplan.dot(temp) / 2.0;
            let temp = 0.5 * (pts1 + pts2) - ptcur;
            d[0][1] = dnplan.dot(temp) - normtgcur;
            d[0][2] = nplan.dot(d1u2) / 2.0;
            d[0][3] = nplan.dot(d1v2) / 2.0;
        } else {
            let (p1, du1, dv1, d2u_1, d2uv_1, d2v_1) = self.surf1.derivatives2(x[2], x[3]);
            pts1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = d2u_1;
            d2uv1 = d2uv_1;
            d2v1 = d2v_1;
            let (p2, du2, dv2, d2u_2, d2uv_2, d2v_2) = self.surf2.derivatives2(p2d.x, p2d.y);
            pts2 = p2;
            d1u2 = du2;
            d1v2 = dv2;
            d2u2 = d2u_2;
            d2uv2 = d2uv_2;
            d2v2 = d2v_2;
            let temp = v2d.x * d1u2 + v2d.y * d1v2;
            d[0][0] = nplan.dot(temp) / 2.0;
            let temp = 0.5 * (pts1 + pts2) - ptcur;
            d[0][1] = dnplan.dot(temp) - normtgcur;
            d[0][2] = nplan.dot(d1u1) / 2.0;
            d[0][3] = nplan.dot(d1v1) / 2.0;
        }

        f[0] = (nplan.x * (pts1.x + pts2.x) + nplan.y * (pts1.y + pts2.y)
            + nplan.z * (pts1.z + pts2.z))
            / 2.0
            + the_d;

        let mut ns1v = d1u1.cross(d1v1);
        if ns1v.length() < EPS {
            if self.first {
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (p2d.x, p2d.y), &mut n);
                ns1v = n;
            } else {
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf1, (x[2], x[3]), &mut n);
                ns1v = n;
            }
        }

        let mut ns2v = d1u2.cross(d1v2);
        if ns2v.length() < EPS {
            if !self.first {
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (p2d.x, p2d.y), &mut n);
                ns2v = n;
            } else {
                let mut n = DVec3::ZERO;
                blend_func_compute_normal(self.surf2, (x[2], x[3]), &mut n);
                ns2v = n;
            }
        }

        let ncrossns1 = nplan.cross(ns1v);
        let ncrossns2 = nplan.cross(ns2v);
        norm1 = ncrossns1.length();
        norm2 = ncrossns2.length();
        if norm1 < EPS {
            norm1 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }
        if norm2 < EPS {
            norm2 = 1.0; // Unsatisfactory, but it is not necessary to stop
        }

        ndotns1 = nplan.dot(ns1v);
        ndotns2 = nplan.dot(ns2v);

        // OCCT L534-541:
        // temp.SetLinearForm(ndotns1 / norm1, nplan, -1. / norm1, ns1);
        // resul1.SetLinearForm(ray1, temp, gp_Vec(pts2, pts1));
        // temp.SetLinearForm(ndotns2 / norm2, nplan, -1. / norm2, ns2);
        // resul1.Subtract(ray2 * temp);
        let temp = (ndotns1 / norm1) * nplan + (-1.0 / norm1) * ns1v;
        let mut resul1 = self.ray1 * temp + (pts1 - pts2);
        let temp = (ndotns2 / norm2) * nplan + (-1.0 / norm2) * ns2v;
        resul1 -= self.ray2 * temp;

        f[1] = resul1.x;
        f[2] = resul1.y;
        f[3] = resul1.z;

        // Derived compared to u1

        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        grosterme = ncrossns1.dot(nplan.cross(temp)) / norm1 / norm1;
        let resul1 = (-self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(temp))) * nplan
            + (self.ray1 * grosterme / norm1) * ns1v
            + (-self.ray1 / norm1) * temp
            + d1u1;

        // Derived compared to v1

        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        grosterme = ncrossns1.dot(nplan.cross(temp)) / norm1 / norm1;
        let resul2 = (-self.ray1 / norm1 * (grosterme * ndotns1 - nplan.dot(temp))) * nplan
            + (self.ray1 * grosterme / norm1) * ns1v
            + (-self.ray1 / norm1) * temp
            + d1v1;

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

        // derived compared to w (parameter on guideline)
        // It is assumed that the raduis is constant

        grosterme = ncrossns1.dot(dnplan.cross(ns1v)) / norm1 / norm1;
        let resul1 = (-self.ray1 / norm1 * (grosterme * ndotns1 - dnplan.dot(ns1v))) * nplan
            + (self.ray1 * ndotns1 / norm1) * dnplan
            + (self.ray1 * grosterme / norm1) * ns1v;

        grosterme = ncrossns2.dot(dnplan.cross(ns2v)) / norm2 / norm2;
        let resul2 = (self.ray2 / norm2 * (grosterme * ndotns2 - dnplan.dot(ns2v))) * nplan
            + (-self.ray2 * ndotns2 / norm2) * dnplan
            + (-self.ray2 * grosterme / norm2) * ns2v;

        d[1][1] = resul1.x + resul2.x;
        d[2][1] = resul1.y + resul2.y;
        d[3][1] = resul1.z + resul2.z;

        // Derived compared to u2
        let mut temp = d2u2.cross(d1v2) + d1u2.cross(d2uv2);
        grosterme = ncrossns2.dot(nplan.cross(temp)) / norm2 / norm2;
        let mut resul1 = (self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(temp))) * nplan
            + (-self.ray2 * grosterme / norm2) * ns2v
            + (self.ray2 / norm2) * temp;
        resul1 -= d1u2;

        // Derived compared to v2
        temp = d2uv2.cross(d1v2) + d1u2.cross(d2v2);
        grosterme = ncrossns2.dot(nplan.cross(temp)) / norm2 / norm2;
        let mut resul2 = (self.ray2 / norm2 * (grosterme * ndotns2 - nplan.dot(temp))) * nplan
            + (-self.ray2 * grosterme / norm2) * ns2v
            + (self.ray2 / norm2) * temp;
        resul2 -= d1v2;

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

impl<'a> FunctionSetWithDerivatives for BlendFuncConstRadInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: inherited from Blend_FuncInv::NbVariables (returns 4).
        BlendFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendFuncConstRadInv::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendFuncConstRadInv::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncConstRadInv::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncConstRadInv::values(self, x, f, df)
    }
}

impl<'a> BlendFuncInv for BlendFuncConstRadInv<'a> {
    fn set_curve_on_surface(&mut self, on_first: bool, c_on_surf: &Curve2d) {
        // OCCT ConstRadInv.cxx L73-77; the rcad port stores the curve
        // reference (OCCT copies the handle).  SAFETY: the caller owns the
        // Curve2d for the lifetime 'a of this function object — same
        // invariant as the OCCT handle.
        let c: &'a Curve2d = unsafe { &*(c_on_surf as *const Curve2d) };
        BlendFuncConstRadInv::set_curve_on_surface(self, on_first, c)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendFuncConstRadInv::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendFuncConstRadInv::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendFuncConstRadInv::is_solution(self, sol, tol)
    }
}
