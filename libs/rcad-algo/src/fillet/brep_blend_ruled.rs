//! OCCT BlendFunc_Ruled (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_Ruled.hxx (fields L186-200) and BlendFunc_Ruled.cxx (whole file
//! L31-744).  OCCT `typedef BlendFunc_Ruled BRepBlend_Ruled`
//! (BRepBlend_Ruled.hxx L22) — the rcad name is [`BlendFuncRuled`].
//!
//! Architecture mappings: `class BlendFunc_Ruled : public Blend_Function`
//! is expressed by implementing the [`BlendFunction`] trait and the
//! `math_FunctionSetWithDerivatives` base ([`FunctionSetWithDerivatives`]);
//! `occ::handle<Adaptor3d_Surface>` / `occ::handle<Adaptor3d_Curve>` map to
//! `&Surface3` / `&Curve3` references; `math_Vector` maps to `&[f64]` /
//! `Vec<f64>` and `math_Matrix` to the row-major `Vec<Vec<f64>>` Jacobian of
//! the rcad solver (OCCT D(i, j) -> d[i - 1][j - 1]); `math_Gauss` maps to
//! [`MathGauss`]; `gp_Vec` / `gp_Vec2d` map to `DVec3` / `DVec2` (OCCT
//! `SetLinearForm(A, V1, B, V2)` becomes `A * V1 + B * V2`).

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::gp::Ax1;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{GeomAbsShape, MatD, VecD};

use super::brep_blend_func::blend_func_next_shape;
use super::brep_blend_function::{BlendAppFunction, BlendFunction};
use super::brep_blend_point::BlendPoint;

/// OCCT BlendFunc_Ruled — function for a ruled blending surface between two
/// surfaces, using a guide line.
pub struct BlendFuncRuled<'a> {
    // OCCT BlendFunc_Ruled.hxx fields (L186-200).
    surf1: &'a Surface3,
    surf2: &'a Surface3,
    curv: &'a Curve3,
    pts1: DVec3,
    pts2: DVec3,
    istangent: bool,
    tg1: DVec3,
    tg12d: DVec2,
    tg2: DVec3,
    tg22d: DVec2,
    ptgui: DVec3,
    d1gui: DVec3,
    d2gui: DVec3,
    nplan: DVec3,
    normtg: f64,
    the_d: f64,
    distmin: f64,
}

impl<'a> BlendFuncRuled<'a> {
    /// OCCT BlendFunc_Ruled(S1, S2, C) (BlendFunc_Ruled.cxx L31-44).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncRuled {
            surf1: s1,
            surf2: s2,
            curv: c,
            pts1: DVec3::ZERO,
            pts2: DVec3::ZERO,
            istangent: true,
            tg1: DVec3::ZERO,
            tg12d: DVec2::ZERO,
            tg2: DVec3::ZERO,
            tg22d: DVec2::ZERO,
            ptgui: DVec3::ZERO,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            normtg: 0.0,
            the_d: 0.0,
            distmin: f64::MAX, // OCCT: RealLast()
        }
    }

    /// OCCT NbEquations() (BlendFunc_Ruled.cxx L46-49) — returns 4.
    pub fn nb_equations(&self) -> usize {
        4
    }

    /// OCCT Set(Param) (BlendFunc_Ruled.cxx L51-58).
    pub fn set_param(&mut self, param: f64) {
        // OCCT: curv->D2(Param, ptgui, d1gui, d2gui);
        self.ptgui = self.curv.point_at(param);
        self.d1gui = self.curv.derivative_at(param);
        self.d2gui = self.curv.derivative2_at(param);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -(self.nplan.dot(self.ptgui));
        self.istangent = true;
    }

    /// OCCT Set(First, Last) (BlendFunc_Ruled.cxx L60-63) — throws
    /// Standard_NotImplemented in OCCT.
    pub fn set_interval(&mut self, _first: f64, _last: f64) {
        panic!("Standard_NotImplemented: BlendFunc_Ruled::Set");
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BlendFunc_Ruled.cxx L65-71) —
    /// UResolution / VResolution of the two surfaces.
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        tolerance[0] = self.surf1.u_resolution(tol);
        tolerance[1] = self.surf1.v_resolution(tol);
        tolerance[2] = self.surf2.u_resolution(tol);
        tolerance[3] = self.surf2.v_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BlendFunc_Ruled.cxx L73-93).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.surf1.default_domain()[0]; // FirstUParameter
        inf_bound[1] = self.surf1.default_domain()[2]; // FirstVParameter
        inf_bound[2] = self.surf2.default_domain()[0];
        inf_bound[3] = self.surf2.default_domain()[2];
        sup_bound[0] = self.surf1.default_domain()[1]; // LastUParameter
        sup_bound[1] = self.surf1.default_domain()[3]; // LastVParameter
        sup_bound[2] = self.surf2.default_domain()[1];
        sup_bound[3] = self.surf2.default_domain()[3];

        for i in 0..4 {
            if !is_infinite_value(inf_bound[i]) && !is_infinite_value(sup_bound[i]) {
                let range = sup_bound[i] - inf_bound[i];
                inf_bound[i] -= range;
                sup_bound[i] += range;
            }
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_Ruled.cxx L95-169).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector valsol(1, 4), secmember(1, 4);
        //       math_Matrix gradsol(1, 4, 1, 4);
        let mut valsol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];
        let mut gradsol = vec![vec![0.0f64; 4]; 4];

        self.values(sol, &mut valsol, &mut gradsol);
        if valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= tol
            && valsol[3].abs() <= tol
        {
            // Calcul des tangentes

            let (_, d1u1, d1v1) = self.surf1.derivatives(sol[0], sol[1]);
            self.pts1 = self.surf1.point_at(sol[0], sol[1]);
            let (_, d1u2, d1v2) = self.surf2.derivatives(sol[2], sol[3]);
            self.pts2 = self.surf2.point_at(sol[2], sol[3]);
            // OCCT: dnplan.SetLinearForm(1./normtg, d2gui,
            //                           -1./normtg * (nplan.Dot(d2gui)), nplan);
            let dnplan =
                (1.0 / self.normtg) * self.d2gui + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            // OCCT: secmember(1) = normtg - dnplan.Dot(gp_Vec(ptgui, pts1));
            secmember[0] = self.normtg - dnplan.dot(self.pts1 - self.ptgui);
            secmember[1] = self.normtg - dnplan.dot(self.pts2 - self.ptgui);

            let ns = d1u1.cross(d1v1);
            let ncrossns = self.nplan.cross(ns);
            let ndotns = self.nplan.dot(ns);
            let norm = ncrossns.length();

            // Derivee de nor1 par rapport au parametre sur la ligne guide
            let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
            // OCCT: temp.SetLinearForm((dnplan.Dot(ns) - grosterme * ndotns) / norm,
            //                          nplan, ndotns / norm, dnplan,
            //                          grosterme / norm, ns);
            let temp = ((dnplan.dot(ns) - grosterme * ndotns) / norm) * self.nplan
                + (ndotns / norm) * dnplan
                + (grosterme / norm) * ns;

            secmember[2] = -(temp.dot(self.pts2 - self.pts1));

            let ns = d1u2.cross(d1v2);
            let ncrossns = self.nplan.cross(ns);
            let ndotns = self.nplan.dot(ns);
            let norm = ncrossns.length();

            // Derivee de nor2 par rapport au parametre sur la ligne guide
            let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
            let temp = ((dnplan.dot(ns) - grosterme * ndotns) / norm) * self.nplan
                + (ndotns / norm) * dnplan
                + (grosterme / norm) * ns;

            secmember[3] = -(temp.dot(self.pts2 - self.pts1));

            // OCCT: math_Gauss Resol(gradsol); — the Jacobian is copied into
            // the rcad 1-based MatD for the solver.
            let mut a = MatD::new(4, 4);
            for r in 1..=4 {
                for c in 1..=4 {
                    a.set(r, c, gradsol[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                let mut x = VecD::new(4);
                for i in 1..=4 {
                    x.set(i, secmember[i - 1]);
                }
                resol.solve(&mut x);

                // OCCT: tg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1);
                self.tg1 = secmember[0] * d1u1 + secmember[1] * d1v1;
                self.tg2 = secmember[2] * d1u2 + secmember[3] * d1v2;
                self.tg12d = DVec2::new(secmember[0], secmember[1]);
                self.tg22d = DVec2::new(secmember[2], secmember[3]);
                self.istangent = false;
            } else {
                self.istangent = true;
            }
            return true;
        }
        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BlendFunc_Ruled.cxx L173-176).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT Value(X, F) (BlendFunc_Ruled.cxx L178-201).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let (_, d1u1, d1v1) = self.surf1.derivatives(x[0], x[1]);
        self.pts1 = self.surf1.point_at(x[0], x[1]);
        let (_, d1u2, d1v2) = self.surf2.derivatives(x[2], x[3]);
        self.pts2 = self.surf2.point_at(x[2], x[3]);

        let temp = self.pts2 - self.pts1; // OCCT: gp_Vec temp(pts1, pts2)

        let ns1 = d1u1.cross(d1v1);
        let ns2 = d1u2.cross(d1v2);

        let norm1 = self.nplan.cross(ns1).length();
        let norm2 = self.nplan.cross(ns2).length();
        // OCCT: ns1.SetLinearForm(nplan.Dot(ns1) / norm1, nplan, -1. / norm1, ns1);
        let ns1 = (self.nplan.dot(ns1) / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        let ns2 = (self.nplan.dot(ns2) / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        f[0] = self.nplan.dot(self.pts1) + self.the_d;
        f[1] = self.nplan.dot(self.pts2) + self.the_d;

        f[2] = temp.dot(ns1);
        f[3] = temp.dot(ns2);

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_Ruled.cxx L203-295).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: surf1->D2(X(1), X(2), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1);
        // (SurfaceEval::derivatives2 returns P, dPu, dPv, dPu2, dPuv, dPv2.)
        let (pts1, d1u1, d1v1, d2u1, d2uv1, d2v1) = self.surf1.derivatives2(x[0], x[1]);
        self.pts1 = pts1;
        let (pts2, d1u2, d1v2, d2u2, d2uv2, d2v2) = self.surf2.derivatives2(x[2], x[3]);
        self.pts2 = pts2;

        d[0][0] = self.nplan.dot(d1u1);
        d[0][1] = self.nplan.dot(d1v1);
        d[0][2] = 0.0;
        d[0][3] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(d1u2);
        d[1][3] = self.nplan.dot(d1v2);

        let ns1 = d1u1.cross(d1v1);
        let ns2 = d1u2.cross(d1v2);
        let ncrossns1 = self.nplan.cross(ns1);
        let ncrossns2 = self.nplan.cross(ns2);
        let norm1 = ncrossns1.length();
        let norm2 = ncrossns2.length();

        let ndotns1 = self.nplan.dot(ns1);
        let ndotns2 = self.nplan.dot(ns2);

        // OCCT: nor1.SetLinearForm(ndotns1 / norm1, nplan, -1. / norm1, ns1);
        let nor1 = (ndotns1 / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        let nor2 = (ndotns2 / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        let p1p2 = self.pts2 - self.pts1; // OCCT: gp_Vec p1p2(pts1, pts2)

        // Derivee de nor1 par rapport a u1
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        let grosterme = ncrossns1.dot(self.nplan.cross(temp)) / norm1 / norm1;
        // OCCT: resul.SetLinearForm(-(grosterme * ndotns1 - nplan.Dot(temp)) / norm1,
        //                           nplan, grosterme / norm1, ns1, -1. / norm1, temp);
        let resul = (-(grosterme * ndotns1 - self.nplan.dot(temp)) / norm1) * self.nplan
            + (grosterme / norm1) * ns1
            + (-1.0 / norm1) * temp;

        d[2][0] = -(d1u1.dot(nor1)) + p1p2.dot(resul);

        // Derivee par rapport a v1
        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        let grosterme = ncrossns1.dot(self.nplan.cross(temp)) / norm1 / norm1;
        let resul = (-(grosterme * ndotns1 - self.nplan.dot(temp)) / norm1) * self.nplan
            + (grosterme / norm1) * ns1
            + (-1.0 / norm1) * temp;

        d[2][1] = -(d1v1.dot(nor1)) + p1p2.dot(resul);

        d[2][2] = d1u2.dot(nor1);
        d[2][3] = d1v2.dot(nor1);

        d[3][0] = -(d1u2.dot(nor1));
        d[3][1] = -(d1v2.dot(nor1));

        // Derivee de nor2 par rapport a u2
        let temp = d2u2.cross(d1v2) + d1u2.cross(d2uv2);
        let grosterme = ncrossns2.dot(self.nplan.cross(temp)) / norm2 / norm2;
        let resul = (-(grosterme * ndotns2 - self.nplan.dot(temp)) / norm2) * self.nplan
            + (grosterme / norm2) * ns2
            + (-1.0 / norm2) * temp;

        d[3][2] = d1u2.dot(nor2) + p1p2.dot(resul);

        // Derivee par rapport a v2
        let temp = d2uv2.cross(d1v2) + d1u2.cross(d2v2);
        let grosterme = ncrossns2.dot(self.nplan.cross(temp)) / norm2 / norm2;
        let resul = (-(grosterme * ndotns2 - self.nplan.dot(temp)) / norm2) * self.nplan
            + (grosterme / norm2) * ns2
            + (-1.0 / norm2) * temp;

        d[3][3] = d1v2.dot(nor2) + p1p2.dot(resul);

        true
    }

    /// OCCT Values(X, F, D) (BlendFunc_Ruled.cxx L297-394).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: surf1->D2(X(1), X(2), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1);
        let (pts1, d1u1, d1v1, d2u1, d2uv1, d2v1) = self.surf1.derivatives2(x[0], x[1]);
        self.pts1 = pts1;
        let (pts2, d1u2, d1v2, d2u2, d2uv2, d2v2) = self.surf2.derivatives2(x[2], x[3]);
        self.pts2 = pts2;

        let p1p2 = self.pts2 - self.pts1;

        let ns1 = d1u1.cross(d1v1);
        let ns2 = d1u2.cross(d1v2);
        let ncrossns1 = self.nplan.cross(ns1);
        let ncrossns2 = self.nplan.cross(ns2);
        let norm1 = ncrossns1.length();
        let norm2 = ncrossns2.length();

        let ndotns1 = self.nplan.dot(ns1);
        let ndotns2 = self.nplan.dot(ns2);

        let nor1 = (ndotns1 / norm1) * self.nplan + (-1.0 / norm1) * ns1;
        let nor2 = (ndotns2 / norm2) * self.nplan + (-1.0 / norm2) * ns2;

        f[0] = self.nplan.dot(self.pts1) + self.the_d;
        f[1] = self.nplan.dot(self.pts2) + self.the_d;
        f[2] = p1p2.dot(nor1);
        f[3] = p1p2.dot(nor2);

        d[0][0] = self.nplan.dot(d1u1);
        d[0][1] = self.nplan.dot(d1v1);
        d[0][2] = 0.0;
        d[0][3] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(d1u2);
        d[1][3] = self.nplan.dot(d1v2);

        // Derivee de nor1 par rapport a u1
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        let grosterme = ncrossns1.dot(self.nplan.cross(temp)) / norm1 / norm1;
        let resul = (-(grosterme * ndotns1 - self.nplan.dot(temp)) / norm1) * self.nplan
            + (grosterme / norm1) * ns1
            + (-1.0 / norm1) * temp;

        d[2][0] = -(d1u1.dot(nor1)) + p1p2.dot(resul);

        // Derivee par rapport a v1
        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        let grosterme = ncrossns1.dot(self.nplan.cross(temp)) / norm1 / norm1;
        let resul = (-(grosterme * ndotns1 - self.nplan.dot(temp)) / norm1) * self.nplan
            + (grosterme / norm1) * ns1
            + (-1.0 / norm1) * temp;

        d[2][1] = -(d1v1.dot(nor1)) + p1p2.dot(resul);

        d[2][2] = d1u2.dot(nor1);
        d[2][3] = d1v2.dot(nor1);

        d[3][0] = -(d1u2.dot(nor1));
        d[3][1] = -(d1v2.dot(nor1));

        // Derivee de nor2 par rapport a u2
        let temp = d2u2.cross(d1v2) + d1u2.cross(d2uv2);
        let grosterme = ncrossns2.dot(self.nplan.cross(temp)) / norm2 / norm2;
        let resul = (-(grosterme * ndotns2 - self.nplan.dot(temp)) / norm2) * self.nplan
            + (grosterme / norm2) * ns2
            + (-1.0 / norm2) * temp;

        d[3][2] = d1u2.dot(nor2) + p1p2.dot(resul);

        // Derivee par rapport a v2
        let temp = d2uv2.cross(d1v2) + d1u2.cross(d2v2);
        let grosterme = ncrossns2.dot(self.nplan.cross(temp)) / norm2 / norm2;
        let resul = (-(grosterme * ndotns2 - self.nplan.dot(temp)) / norm2) * self.nplan
            + (grosterme / norm2) * ns2
            + (-1.0 / norm2) * temp;

        d[3][3] = d1v2.dot(nor2) + p1p2.dot(resul);

        true
    }

    /// OCCT Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)
    /// (BlendFunc_Ruled.cxx L396-416).
    #[allow(clippy::too_many_arguments)]
    pub fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_f: &mut DVec3,
        tg_l: &mut DVec3,
        nm_f: &mut DVec3,
        nm_l: &mut DVec3,
    ) {
        // OCCT: surf2->D1(U2, V2, bid, d1u, d1v); NmL = d1u.Crossed(d1v);
        let (_, d1u2, d1v2) = self.surf2.derivatives(u2, v2);
        *nm_l = d1u2.cross(d1v2);

        // OCCT: surf1->D1(U1, V1, bid, d1u, d1v); NmF = ns1 = d1u.Crossed(d1v);
        let (_, d1u1, d1v1) = self.surf1.derivatives(u1, v1);
        let ns1 = d1u1.cross(d1v1);
        *nm_f = ns1;

        *tg_f = self.pts2 - self.pts1; // OCCT: TgF = TgL = gp_Vec(pts1, pts2)
        *tg_l = *tg_f;
    }

    /// OCCT PointOnS1() (BlendFunc_Ruled.cxx L418-421).
    pub fn point_on_s1(&self) -> DVec3 {
        self.pts1
    }

    /// OCCT PointOnS2() (BlendFunc_Ruled.cxx L423-426).
    pub fn point_on_s2(&self) -> DVec3 {
        self.pts2
    }

    /// OCCT IsTangencyPoint() (BlendFunc_Ruled.cxx L428-431).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnS1() (BlendFunc_Ruled.cxx L433-440).
    pub fn tangent_on_s1(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_Ruled::TangentOnS1");
        }
        self.tg1
    }

    /// OCCT TangentOnS2() (BlendFunc_Ruled.cxx L442-449).
    pub fn tangent_on_s2(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_Ruled::TangentOnS2");
        }
        self.tg2
    }

    /// OCCT Tangent2dOnS1() (BlendFunc_Ruled.cxx L451-458).
    pub fn tangent_2d_on_s1(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_Ruled::Tangent2dOnS1");
        }
        self.tg12d
    }

    /// OCCT Tangent2dOnS2() (BlendFunc_Ruled.cxx L460-467).
    pub fn tangent_2d_on_s2(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_Ruled::Tangent2dOnS2");
        }
        self.tg22d
    }

    /// OCCT GetSection(Param, U1, V1, U2, V2, tabP, tabV)
    /// (BlendFunc_Ruled.cxx L469-568).
    pub fn get_section(
        &mut self,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tab_p: &mut [DVec3],
        tab_v: &mut [DVec3],
    ) -> bool {
        let nb_point = tab_p.len();
        if nb_point != tab_v.len() || nb_point < 2 {
            panic!("Standard_RangeError");
        }
        let lowp = 0usize; // OCCT: tabP.Lower()
        let lowv = 0usize; // OCCT: tabV.Lower()

        // OCCT: curv->D2(Param, ptgui, d1gui, d2gui);
        self.ptgui = self.curv.point_at(param);
        self.d1gui = self.curv.derivative_at(param);
        self.d2gui = self.curv.derivative2_at(param);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -(self.nplan.dot(self.ptgui));

        let sol = [u1, v1, u2, v2];
        let mut valsol = [0.0f64; 4];
        let mut gradsol = vec![vec![0.0f64; 4]; 4];
        self.values(&sol, &mut valsol, &mut gradsol);

        let (_, d1u1, d1v1) = self.surf1.derivatives(sol[0], sol[1]);
        self.pts1 = self.surf1.point_at(sol[0], sol[1]);
        let (_, d1u2, d1v2) = self.surf2.derivatives(sol[2], sol[3]);
        self.pts2 = self.surf2.point_at(sol[2], sol[3]);
        let dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        let mut secmember = [0.0f64; 4];
        secmember[0] = self.normtg - dnplan.dot(self.pts1 - self.ptgui);
        secmember[1] = self.normtg - dnplan.dot(self.pts2 - self.ptgui);

        let ns = d1u1.cross(d1v1);
        let ncrossns = self.nplan.cross(ns);
        let ndotns = self.nplan.dot(ns);
        let norm = ncrossns.length();

        // Derivee de nor1 par rapport au parametre sur la ligne guide
        let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
        let temp = ((dnplan.dot(ns) - grosterme * ndotns) / norm) * self.nplan
            + (ndotns / norm) * dnplan
            + (grosterme / norm) * ns;

        secmember[2] = -(temp.dot(self.pts2 - self.pts1));

        let ns = d1u2.cross(d1v2);
        let ncrossns = self.nplan.cross(ns);
        let ndotns = self.nplan.dot(ns);
        let norm = ncrossns.length();

        // Derivee de nor2 par rapport au parametre sur la ligne guide
        let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
        let temp = ((dnplan.dot(ns) - grosterme * ndotns) / norm) * self.nplan
            + (ndotns / norm) * dnplan
            + (grosterme / norm) * ns;

        secmember[3] = -(temp.dot(self.pts2 - self.pts1));

        let mut a = MatD::new(4, 4);
        for r in 1..=4 {
            for c in 1..=4 {
                a.set(r, c, gradsol[r - 1][c - 1]);
            }
        }
        let resol = MathGauss::new(&a);
        if resol.is_done() {
            let mut x = VecD::new(4);
            for i in 1..=4 {
                x.set(i, secmember[i - 1]);
            }
            resol.solve(&mut x);
            for i in 1..=4 {
                secmember[i - 1] = x.get(i);
            }

            self.tg1 = secmember[0] * d1u1 + secmember[1] * d1v1;
            self.tg2 = secmember[2] * d1u2 + secmember[3] * d1v2;

            tab_p[lowp] = self.pts1;
            tab_p[lowp + nb_point - 1] = self.pts2;

            tab_v[lowv] = self.tg1;
            tab_v[lowv + nb_point - 1] = self.tg2;

            for i in 2..=nb_point - 1 {
                let lambda = (i - 1) as f64 / (nb_point - 1) as f64;
                tab_p[lowp + i - 1] = (1.0 - lambda) * self.pts1 + lambda * self.pts2;
                tab_v[lowv + i - 1] = (1.0 - lambda) * self.tg1 + lambda * self.tg2;
            }
            return true;
        }
        false
    }

    /// OCCT IsRational() (BlendFunc_Ruled.cxx L572-575).
    pub fn is_rational(&self) -> bool {
        false
    }

    /// OCCT GetSectionSize() (BlendFunc_Ruled.cxx L579-583) — throws
    /// Standard_NotImplemented in OCCT.
    pub fn get_section_size(&self) -> f64 {
        panic!("Standard_NotImplemented: BlendFunc_Ruled::GetSectionSize()");
    }

    /// OCCT GetMinimalWeight(Weigths) (BlendFunc_Ruled.cxx L586-590).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        for w in weigths.iter_mut() {
            *w = 1.0;
        }
    }

    /// OCCT NbIntervals(S) (BlendFunc_Ruled.cxx L593-599) —
    /// `curv->NbIntervals(BlendFunc::NextShape(S))`.
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.curv.nb_intervals(blend_func_next_shape(s))
    }

    /// OCCT Intervals(T, S) (BlendFunc_Ruled.cxx L603-606) —
    /// `curv->Intervals(T, BlendFunc::NextShape(S))`.
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let mut intervals = Vec::new();
        self.curv.intervals(&mut intervals, blend_func_next_shape(s));
        for (dst, src) in t.iter_mut().zip(intervals) {
            *dst = src;
        }
    }

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
    /// (BlendFunc_Ruled.cxx L607-616).
    pub fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        *nb_poles = 2;
        *nb_knots = 2;
        *degree = 1;
        *nb_poles_2d = 2;
    }

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)
    /// (BlendFunc_Ruled.cxx L619-626).
    pub fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        _surf_tol: f64,
        _angle_tol: f64,
        tol3d: &mut [f64],
        _tol1d: &mut [f64],
    ) {
        for v in tol3d.iter_mut() {
            *v = bound_tol;
        }
    }

    /// OCCT Knots(TKnots) (BlendFunc_Ruled.cxx L628-632).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        tknots[0] = 0.0;
        tknots[tknots.len() - 1] = 1.0;
    }

    /// OCCT Mults(TMults) (BlendFunc_Ruled.cxx L634-637).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        tmults[0] = 2;
        tmults[tmults.len() - 1] = 2;
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)
    /// (BlendFunc_Ruled.cxx L653-691) — used for the first and last section.
    pub fn section_d1(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool {
        let lowp = 0usize; // OCCT: Poles.Lower()
        let lowp_2d = 0usize; // OCCT: Poles2d.Lower()

        poles[lowp] = p.point_on_s1();
        poles[lowp + 1] = p.point_on_s2();

        let (u, v) = p.parameters_on_s1();
        poles_2d[lowp_2d] = DVec2::new(u, v);
        let (u, v) = p.parameters_on_s2();
        poles_2d[lowp_2d + 1] = DVec2::new(u, v);

        weigths[lowp] = 1.0;
        weigths[lowp + 1] = 1.0;

        if !p.is_tangency_point() {
            d_poles[lowp] = p.tangent_on_s1();
            d_poles[lowp + 1] = p.tangent_on_s2();

            d_poles_2d[lowp_2d] = p.tangent_2d_on_s1();
            d_poles_2d[lowp_2d + 1] = p.tangent_2d_on_s2();

            d_weigths[lowp] = 0.0;
            d_weigths[lowp + 1] = 0.0;

            return true;
        }

        false
    }

    /// OCCT Section(P, Poles, Poles2d, Weights) (BlendFunc_Ruled.cxx
    /// L693-711).
    pub fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        let lowp = 0usize; // OCCT: Poles.Lower()
        let lowp_2d = 0usize; // OCCT: Poles2d.Lower()

        poles[lowp] = p.point_on_s1();
        poles[lowp + 1] = p.point_on_s2();

        let (u, v) = p.parameters_on_s1();
        poles_2d[lowp_2d] = DVec2::new(u, v);
        let (u, v) = p.parameters_on_s2();
        poles_2d[lowp_2d + 1] = DVec2::new(u, v);

        weigths[lowp] = 1.0;
        weigths[lowp + 1] = 1.0;
    }

    /// OCCT AxeRot(Prm) (BlendFunc_Ruled.cxx L713-730).
    pub fn axe_rot(&mut self, prm: f64) -> Ax1 {
        let mut axrot = Ax1::new(DVec3::ZERO, DVec3::Z);

        // OCCT: curv->D2(Prm, ptgui, d1gui, d2gui);
        self.ptgui = self.curv.point_at(prm);
        self.d1gui = self.curv.derivative_at(prm);
        self.d2gui = self.curv.derivative2_at(prm);

        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        let dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        let dirax = self.nplan.cross(dnplan);
        axrot.set_direction(dirax);
        let oriax = self.ptgui + (self.normtg / dnplan.length()) * dnplan.normalize_or_zero();
        axrot.set_location(oriax);
        axrot
    }

    /// OCCT Resolution(IC2d, Tol, TolU, TolV) (BlendFunc_Ruled.cxx L732-744).
    pub fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic_2d == 1 {
            *tol_u = self.surf1.u_resolution(tol);
            *tol_v = self.surf1.v_resolution(tol);
        } else {
            *tol_u = self.surf2.u_resolution(tol);
            *tol_v = self.surf2.v_resolution(tol);
        }
    }
}

impl<'a> FunctionSetWithDerivatives for BlendFuncRuled<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: inherited from Blend_Function::NbVariables (returns 4).
        BlendFunction::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendFuncRuled::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendFuncRuled::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncRuled::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncRuled::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendFuncRuled<'a> {
    fn set_param(&mut self, param: f64) {
        BlendFuncRuled::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendFuncRuled::set_interval(self, first, last)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendFuncRuled::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendFuncRuled::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendFuncRuled::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendFuncRuled::get_minimal_distance(self)
    }

    fn pnt1(&self) -> DVec3 {
        BlendFuncRuled::point_on_s1(self)
    }

    fn pnt2(&self) -> DVec3 {
        BlendFuncRuled::point_on_s2(self)
    }

    fn is_rational(&self) -> bool {
        BlendFuncRuled::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendFuncRuled::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendFuncRuled::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendFuncRuled::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendFuncRuled::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendFuncRuled::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendFuncRuled::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendFuncRuled::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendFuncRuled::mults(self, tmults)
    }

    fn section_d1(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool {
        BlendFuncRuled::section_d1(
            self, p, poles, d_poles, poles_2d, d_poles_2d, weigths, d_weigths,
        )
    }

    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        BlendFuncRuled::section(self, p, poles, poles_2d, weigths)
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weights, DWeights, D2Weights) (BlendFunc_Ruled.cxx L639-651) returns
    /// false; the body is identical to the Blend_Function implementation
    /// (Blend_Function.cxx L44-56), which the Rust subtrait cannot inherit
    /// for this struct, so the body is restated here.
    fn section_d2(
        &mut self,
        _p: &BlendPoint,
        _poles: &mut [DVec3],
        _d_poles: &mut [DVec3],
        _d2_poles: &mut [DVec3],
        _poles_2d: &mut [DVec2],
        _d_poles_2d: &mut [DVec2],
        _d2_poles_2d: &mut [DVec2],
        _weigths: &mut [f64],
        _d_weigths: &mut [f64],
        _d2_weigths: &mut [f64],
    ) -> bool {
        false
    }

    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        BlendFuncRuled::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendFunction for BlendFuncRuled<'a> {
    fn point_on_s1(&self) -> DVec3 {
        BlendFuncRuled::point_on_s1(self)
    }

    fn point_on_s2(&self) -> DVec3 {
        BlendFuncRuled::point_on_s2(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendFuncRuled::is_tangency_point(self)
    }

    fn tangent_on_s1(&self) -> DVec3 {
        BlendFuncRuled::tangent_on_s1(self)
    }

    fn tangent_2d_on_s1(&self) -> DVec2 {
        BlendFuncRuled::tangent_2d_on_s1(self)
    }

    fn tangent_on_s2(&self) -> DVec3 {
        BlendFuncRuled::tangent_on_s2(self)
    }

    fn tangent_2d_on_s2(&self) -> DVec2 {
        BlendFuncRuled::tangent_2d_on_s2(self)
    }

    fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut DVec3,
        tg_last: &mut DVec3,
        norm_first: &mut DVec3,
        norm_last: &mut DVec3,
    ) {
        BlendFuncRuled::tangent(self, u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last)
    }
}
