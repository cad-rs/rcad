//! OCCT BlendFunc_ConstThroat family (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_ConstThroat.hxx (fields L85-112) + BlendFunc_ConstThroat.cxx
//! (whole file L27-280); BlendFunc_ConstThroatWithPenetration.hxx +
//! BlendFunc_ConstThroatWithPenetration.cxx (whole file L27-185);
//! BlendFunc_ConstThroatInv.hxx + BlendFunc_ConstThroatInv.cxx (whole file
//! L22-237); BlendFunc_ConstThroatWithPenetrationInv.hxx +
//! BlendFunc_ConstThroatWithPenetrationInv.cxx (whole file L25-199).  The
//! structs are defined in [`super::brep_blend_func_chamfer`]; this file
//! carries their method bodies and trait wiring.
//!
//! Architecture mappings: the GenChamfer / GenChamfInv inherited bodies are
//! stamped by the macros of [`super::brep_blend_func_chamfer`]; the bodies
//! ConstThroatWithPenetration / ConstThroatWithPenetrationInv inherit from
//! ConstThroat / ConstThroatInv (identical OCCT text) are stamped by the
//! `const_throat_common!` / `const_throat_inv_common!` macros below.

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{
    Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _,
};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::gp::Lin;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{GeomAbsShape, MatD, VecD};

use super::brep_blend_func::blend_func_next_shape;
use super::brep_blend_func_chamfer::{
    gen_chamf_inv_common, gen_chamfer_common, impl_gen_chamf_inv_traits, impl_gen_chamfer_traits,
    BlendFuncConstThroat, BlendFuncConstThroatInv, BlendFuncConstThroatWithPenetration,
    BlendFuncConstThroatWithPenetrationInv, BlendFuncGenChamfInv, BlendFuncGenChamfer,
};
use super::brep_blend_func_inv::BlendFuncInv;
use super::brep_blend_function::{BlendAppFunction, BlendFunction};
use super::brep_blend_point::BlendPoint;

impl<'a> BlendFuncConstThroat<'a> {
    /// OCCT BlendFunc_ConstThroat(S1, S2, C) (BlendFunc_ConstThroat.cxx
    /// L27-37) and the inherited BlendFunc_GenChamfer ctor
    /// (GenChamfer.cxx L29-39).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncConstThroat {
            surf1: s1,
            surf2: s2,
            curv: c,
            choix: 0,
            tol: 0.0,
            distmin: f64::MAX, // OCCT: RealLast()
            pts1: DVec3::ZERO,
            pts2: DVec3::ZERO,
            d1u1: DVec3::ZERO,
            d1v1: DVec3::ZERO,
            d1u2: DVec3::ZERO,
            d1v2: DVec3::ZERO,
            istangent: false,
            tg1: DVec3::ZERO,
            tg12d: DVec2::ZERO,
            tg2: DVec3::ZERO,
            tg22d: DVec2::ZERO,
            param: 0.0,
            throat: 0.0,
            ptgui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            normtg: 0.0,
            the_d: 0.0,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
        }
    }
}

impl<'a> BlendFuncConstThroatWithPenetration<'a> {
    /// OCCT BlendFunc_ConstThroatWithPenetration(S1, S2, C)
    /// (BlendFunc_ConstThroatWithPenetration.cxx L27-33) — derives from
    /// BlendFunc_ConstThroat.
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncConstThroatWithPenetration {
            surf1: s1,
            surf2: s2,
            curv: c,
            choix: 0,
            tol: 0.0,
            distmin: f64::MAX, // OCCT: RealLast()
            pts1: DVec3::ZERO,
            pts2: DVec3::ZERO,
            d1u1: DVec3::ZERO,
            d1v1: DVec3::ZERO,
            d1u2: DVec3::ZERO,
            d1v2: DVec3::ZERO,
            istangent: false,
            tg1: DVec3::ZERO,
            tg12d: DVec2::ZERO,
            tg2: DVec3::ZERO,
            tg22d: DVec2::ZERO,
            param: 0.0,
            throat: 0.0,
            ptgui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            normtg: 0.0,
            the_d: 0.0,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
        }
    }
}

/// Stamps the bodies BlendFunc_ConstThroatWithPenetration inherits from
/// BlendFunc_ConstThroat (Set(aThroat, ...), Set(Param), the point/tangent
/// accessors, Tangent and GetSectionSize — identical OCCT text).
#[allow(clippy::too_many_arguments)]
macro_rules! const_throat_common {
    ($t:ident) => {
        impl<'a> $t<'a> {
            /// OCCT Set(aThroat, ..., Choix) (BlendFunc_ConstThroat.cxx
            /// L41-45) — sets the throat and the "quadrant".
            pub fn set(&mut self, a_throat: f64, _unused: f64, choix: i32) {
                self.throat = a_throat;
                self.choix = choix;
            }

            /// OCCT Set(Param) (BlendFunc_ConstThroat.cxx L49-56).
            pub fn set_param(&mut self, param: f64) {
                self.param = param;
                // OCCT: curv->D2(param, ptgui, d1gui, d2gui);
                self.ptgui = self.curv.point_at(param);
                self.d1gui = self.curv.derivative_at(param);
                self.d2gui = self.curv.derivative2_at(param);
                self.normtg = self.d1gui.length();
                self.nplan = self.d1gui.normalize_or_zero();
                self.the_d = -self.nplan.dot(self.ptgui);
            }

            /// OCCT PointOnS1() (BlendFunc_ConstThroat.cxx L163-166).
            pub fn point_on_s1(&self) -> DVec3 {
                self.pts1
            }

            /// OCCT PointOnS2() (BlendFunc_ConstThroat.cxx L170-173).
            pub fn point_on_s2(&self) -> DVec3 {
                self.pts2
            }

            /// OCCT IsTangencyPoint() (BlendFunc_ConstThroat.cxx L177-180).
            pub fn is_tangency_point(&self) -> bool {
                self.istangent
            }

            /// OCCT TangentOnS1() (BlendFunc_ConstThroat.cxx L184-192).
            pub fn tangent_on_s1(&self) -> DVec3 {
                if self.istangent {
                    panic!("Standard_DomainError: BlendFunc_ConstThroat::TangentOnS1");
                }
                self.tg1
            }

            /// OCCT TangentOnS2() (BlendFunc_ConstThroat.cxx L195-202).
            pub fn tangent_on_s2(&self) -> DVec3 {
                if self.istangent {
                    panic!("Standard_DomainError: BlendFunc_ConstThroat::TangentOnS2");
                }
                self.tg2
            }

            /// OCCT Tangent2dOnS1() (BlendFunc_ConstThroat.cxx L206-213).
            pub fn tangent_2d_on_s1(&self) -> DVec2 {
                if self.istangent {
                    panic!("Standard_DomainError: BlendFunc_ConstThroat::Tangent2dOnS1");
                }
                self.tg12d
            }

            /// OCCT Tangent2dOnS2() (BlendFunc_ConstThroat.cxx L217-224).
            pub fn tangent_2d_on_s2(&self) -> DVec2 {
                if self.istangent {
                    panic!("Standard_DomainError: BlendFunc_ConstThroat::Tangent2dOnS2");
                }
                self.tg22d
            }

            /// OCCT Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)
            /// (BlendFunc_ConstThroat.cxx L228-273).
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
                let mut rev_f = false;
                let mut rev_l = false;

                // OCCT: surf1->D1(U1, V1, pt, d1u, d1v); NmF = d1u.Crossed(d1v);
                let (_, d1u1, d1v1) = self.surf1.derivatives(u1, v1);
                *nm_f = d1u1.cross(d1v1);

                // OCCT: surf2->D1(U2, V2, pt, d1u, d1v); NmL = d1u.Crossed(d1v);
                let (_, d1u2, d1v2) = self.surf2.derivatives(u2, v2);
                *nm_l = d1u2.cross(d1v2);

                *tg_f = self.nplan.cross(*nm_f).normalize_or_zero();
                *tg_l = self.nplan.cross(*nm_l).normalize_or_zero();

                if (self.choix == 2) || (self.choix == 5) {
                    rev_f = true;
                    rev_l = true;
                }
                if (self.choix == 4) || (self.choix == 7) {
                    rev_l = true;
                }
                if (self.choix == 3) || (self.choix == 8) {
                    rev_f = true;
                }

                if rev_f {
                    *tg_f = -*tg_f;
                }
                if rev_l {
                    *tg_l = -*tg_l;
                }
            }

            /// OCCT GetSectionSize() (BlendFunc_ConstThroat.cxx L277-280).
            pub fn get_section_size(&self) -> f64 {
                panic!("Standard_NotImplemented: BlendFunc_ConstThroat::GetSectionSize()");
            }
        }
    };
}

const_throat_common!(BlendFuncConstThroat);
const_throat_common!(BlendFuncConstThroatWithPenetration);

impl<'a> BlendFuncConstThroat<'a> {
    /// OCCT IsSolution(Sol, Tol) (BlendFunc_ConstThroat.cxx L60-109).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector secmember(1, 4), valsol(1, 4);
        //       math_Matrix gradsol(1, 4, 1, 4);
        let mut valsol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];
        let mut gradsol = vec![vec![0.0f64; 4]; 4];

        self.value(sol, &mut valsol);
        self.derivatives(sol, &mut gradsol);

        self.tol = tol;

        if valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= tol * tol
            && valsol[3].abs() <= tol * tol
        {
            // OCCT: dnplan.SetLinearForm(1./normtg, d2gui,
            //                           -1./normtg * (nplan.Dot(d2gui)), nplan);
            let dnplan = (1.0 / self.normtg) * self.d2gui
                + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            let temp1 = self.pts1 - self.ptgui;
            let temp2 = self.pts2 - self.ptgui;
            let tempmid = (self.pts1 + self.pts2) / 2.0 - self.ptgui;
            // OCCT: surf1->D1(Sol(1), Sol(2), pts1, d1u1, d1v1);
            let (_, du1, dv1) = self.surf1.derivatives(sol[0], sol[1]);
            self.d1u1 = du1;
            self.d1v1 = dv1;
            // OCCT: surf2->D1(Sol(3), Sol(4), pts2, d1u2, d1v2);
            let (_, du2, dv2) = self.surf2.derivatives(sol[2], sol[3]);
            self.d1u2 = du2;
            self.d1v2 = dv2;

            secmember[0] = self.nplan.dot(self.d1gui) - dnplan.dot(temp1);
            secmember[1] = self.nplan.dot(self.d1gui) - dnplan.dot(temp2);
            secmember[2] = 2.0 * self.d1gui.dot(tempmid);
            secmember[3] = 2.0 * self.d1gui.dot(temp2) - 2.0 * self.d1gui.dot(temp1);

            // OCCT: math_Gauss Resol(gradsol);
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
                // OCCT: tg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1);
                self.tg1 = secmember[0] * self.d1u1 + secmember[1] * self.d1v1;
                self.tg2 = secmember[2] * self.d1u2 + secmember[3] * self.d1v2;
                self.tg12d = DVec2::new(secmember[0], secmember[1]);
                self.tg22d = DVec2::new(secmember[2], secmember[3]);
                self.istangent = false;
            } else {
                self.istangent = true;
            }

            self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

            return true;
        }

        false
    }

    /// OCCT Value(X, F) (BlendFunc_ConstThroat.cxx L113-132).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: surf1->D0(X(1), X(2), pts1);
        self.pts1 = self.surf1.point_at(x[0], x[1]);
        // OCCT: surf2->D0(X(3), X(4), pts2);
        self.pts2 = self.surf2.point_at(x[2], x[3]);

        f[0] = self.nplan.dot(self.pts1) + self.the_d;
        f[1] = self.nplan.dot(self.pts2) + self.the_d;

        // OCCT: const gp_Pnt ptmid((pts1.XYZ() + pts2.XYZ()) / 2);
        //       const gp_Vec vmid(ptgui, ptmid);
        let ptmid = (self.pts1 + self.pts2) / 2.0;
        let vmid = ptmid - self.ptgui;

        f[2] = vmid.length_squared() - self.throat * self.throat;

        // OCCT: const gp_Vec vref1(ptgui, pts1); const gp_Vec vref2(ptgui, pts2);
        let vref1 = self.pts1 - self.ptgui;
        let vref2 = self.pts2 - self.ptgui;

        f[3] = vref1.length_squared() - vref2.length_squared();

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_ConstThroat.cxx L136-159).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: surf1->D1(X(1), X(2), pts1, d1u1, d1v1);
        let (_, du1, dv1) = self.surf1.derivatives(x[0], x[1]);
        self.d1u1 = du1;
        self.d1v1 = dv1;
        // OCCT: surf2->D1(X(3), X(4), pts2, d1u2, d1v2);
        let (_, du2, dv2) = self.surf2.derivatives(x[2], x[3]);
        self.d1u2 = du2;
        self.d1v2 = dv2;

        // OCCT: gp_Vec((pts1.XYZ() + pts2.XYZ()) / 2 - ptgui.XYZ())
        let pmid_ptgui = (self.pts1 + self.pts2) / 2.0 - self.ptgui;
        let ptp1 = self.pts1 - self.ptgui; // OCCT: gp_Vec(ptgui, pts1)
        let ptp2 = self.pts2 - self.ptgui; // OCCT: gp_Vec(ptgui, pts2)

        d[0][0] = self.nplan.dot(self.d1u1);
        d[0][1] = self.nplan.dot(self.d1v1);
        d[0][2] = 0.0;
        d[0][3] = 0.0;
        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(self.d1u2);
        d[1][3] = self.nplan.dot(self.d1v2);
        d[2][0] = pmid_ptgui.dot(self.d1u1);
        d[2][1] = pmid_ptgui.dot(self.d1v1);
        d[2][2] = pmid_ptgui.dot(self.d1u2);
        d[2][3] = pmid_ptgui.dot(self.d1v2);
        d[3][0] = 2.0 * ptp1.dot(self.d1u1);
        d[3][1] = 2.0 * ptp1.dot(self.d1v1);
        d[3][2] = -2.0 * ptp2.dot(self.d1u2);
        d[3][3] = -2.0 * ptp2.dot(self.d1v2);

        true
    }
}

gen_chamfer_common!(BlendFuncConstThroat);
impl_gen_chamfer_traits!(BlendFuncConstThroat);

impl<'a> BlendFuncConstThroatWithPenetration<'a> {
    /// OCCT IsSolution(Sol, Tol)
    /// (BlendFunc_ConstThroatWithPenetration.cxx L37-86).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector secmember(1, 4), valsol(1, 4);
        //       math_Matrix gradsol(1, 4, 1, 4);
        let mut valsol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];
        let mut gradsol = vec![vec![0.0f64; 4]; 4];

        self.value(sol, &mut valsol);
        self.derivatives(sol, &mut gradsol);

        self.tol = tol;

        if valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= tol * tol
            && valsol[3].abs() <= tol
        {
            // OCCT: dnplan.SetLinearForm(1./normtg, d2gui,
            //                           -1./normtg * (nplan.Dot(d2gui)), nplan);
            let dnplan = (1.0 / self.normtg) * self.d2gui
                + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            let temp1 = self.pts1 - self.ptgui;
            let temp2 = self.pts2 - self.ptgui;
            let temp3 = self.pts2 - self.pts1;
            // OCCT: surf1->D1(Sol(1), Sol(2), pts1, d1u1, d1v1);
            let (_, du1, dv1) = self.surf1.derivatives(sol[0], sol[1]);
            self.d1u1 = du1;
            self.d1v1 = dv1;
            // OCCT: surf2->D1(Sol(3), Sol(4), pts2, d1u2, d1v2);
            let (_, du2, dv2) = self.surf2.derivatives(sol[2], sol[3]);
            self.d1u2 = du2;
            self.d1v2 = dv2;

            secmember[0] = self.nplan.dot(self.d1gui) - dnplan.dot(temp1);
            secmember[1] = self.nplan.dot(self.d1gui) - dnplan.dot(temp2);
            secmember[2] = 2.0 * self.d1gui.dot(temp1);
            secmember[3] = self.d1gui.dot(temp3);

            // OCCT: math_Gauss Resol(gradsol);
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
                // OCCT: tg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1);
                self.tg1 = secmember[0] * self.d1u1 + secmember[1] * self.d1v1;
                self.tg2 = secmember[2] * self.d1u2 + secmember[3] * self.d1v2;
                self.tg12d = DVec2::new(secmember[0], secmember[1]);
                self.tg22d = DVec2::new(secmember[2], secmember[3]);
                self.istangent = false;
            } else {
                self.istangent = true;
            }

            self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

            return true;
        }

        false
    }

    /// OCCT Value(X, F) (BlendFunc_ConstThroatWithPenetration.cxx L90-107).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: surf1->D0(X(1), X(2), pts1);
        self.pts1 = self.surf1.point_at(x[0], x[1]);
        // OCCT: surf2->D0(X(3), X(4), pts2);
        self.pts2 = self.surf2.point_at(x[2], x[3]);

        f[0] = self.nplan.dot(self.pts1) + self.the_d;
        f[1] = self.nplan.dot(self.pts2) + self.the_d;

        // OCCT: const gp_Vec vref(ptgui, pts1);
        let vref = self.pts1 - self.ptgui;

        f[2] = vref.length_squared() - self.throat * self.throat;

        // OCCT: const gp_Vec vec12(pts1, pts2);
        let vec12 = self.pts2 - self.pts1;

        f[3] = vref.dot(vec12);

        true
    }

    /// OCCT Derivatives(X, D)
    /// (BlendFunc_ConstThroatWithPenetration.cxx L111-134).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: surf1->D1(X(1), X(2), pts1, d1u1, d1v1);
        let (_, du1, dv1) = self.surf1.derivatives(x[0], x[1]);
        self.d1u1 = du1;
        self.d1v1 = dv1;
        // OCCT: surf2->D1(X(3), X(4), pts2, d1u2, d1v2);
        let (_, du2, dv2) = self.surf2.derivatives(x[2], x[3]);
        self.d1u2 = du2;
        self.d1v2 = dv2;

        let ptp1 = self.pts1 - self.ptgui; // OCCT: gp_Vec(ptgui, pts1)
        let pts1_pts2 = self.pts2 - self.pts1; // OCCT: gp_Vec(pts1, pts2)

        d[0][0] = self.nplan.dot(self.d1u1);
        d[0][1] = self.nplan.dot(self.d1v1);
        d[0][2] = 0.0;
        d[0][3] = 0.0;
        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(self.d1u2);
        d[1][3] = self.nplan.dot(self.d1v2);
        d[2][0] = 2.0 * ptp1.dot(self.d1u1);
        d[2][1] = 2.0 * ptp1.dot(self.d1v1);
        d[2][2] = 0.0;
        d[2][3] = 0.0;
        // OCCT: D(4, 1) = d1u1.Dot(gp_Vec(pts1, pts2)) - gp_Vec(ptgui, pts1).Dot(d1u1);
        d[3][0] = self.d1u1.dot(pts1_pts2) - ptp1.dot(self.d1u1);
        d[3][1] = self.d1v1.dot(pts1_pts2) - ptp1.dot(self.d1v1);
        d[3][2] = ptp1.dot(self.d1u2);
        d[3][3] = ptp1.dot(self.d1v2);

        true
    }
}

gen_chamfer_common!(BlendFuncConstThroatWithPenetration);
impl_gen_chamfer_traits!(BlendFuncConstThroatWithPenetration);

impl<'a> BlendFuncConstThroatInv<'a> {
    /// OCCT BlendFunc_ConstThroatInv(S1, S2, C)
    /// (BlendFunc_ConstThroatInv.cxx L22-33) and the inherited
    /// BlendFunc_GenChamfInv ctor (GenChamfInv.cxx L23-32).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncConstThroatInv {
            surf1: s1,
            surf2: s2,
            curv: c,
            csurf: None,
            choix: 0,
            first: false,
            throat: 0.0,
            param: 0.0,
            sign1: 0.0,
            sign2: 0.0,
            ptgui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            normtg: 0.0,
            the_d: 0.0,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            pts1: DVec3::ZERO,
            pts2: DVec3::ZERO,
            d1u1: DVec3::ZERO,
            d1v1: DVec3::ZERO,
            d1u2: DVec3::ZERO,
            d1v2: DVec3::ZERO,
        }
    }
}

impl<'a> BlendFuncConstThroatWithPenetrationInv<'a> {
    /// OCCT BlendFunc_ConstThroatWithPenetrationInv(S1, S2, C)
    /// (BlendFunc_ConstThroatWithPenetrationInv.cxx L27-40) — derives from
    /// BlendFunc_ConstThroatInv.
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncConstThroatWithPenetrationInv {
            surf1: s1,
            surf2: s2,
            curv: c,
            csurf: None,
            choix: 0,
            first: false,
            throat: 0.0,
            param: 0.0,
            sign1: 0.0,
            sign2: 0.0,
            ptgui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            normtg: 0.0,
            the_d: 0.0,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            pts1: DVec3::ZERO,
            pts2: DVec3::ZERO,
            d1u1: DVec3::ZERO,
            d1v1: DVec3::ZERO,
            d1u2: DVec3::ZERO,
            d1v2: DVec3::ZERO,
        }
    }
}

/// Stamps the body BlendFunc_ConstThroatWithPenetrationInv inherits from
/// BlendFunc_ConstThroatInv (Set(theThroat, ..., Choix) — identical OCCT
/// text).
macro_rules! const_throat_inv_common {
    ($t:ident) => {
        impl<'a> $t<'a> {
            /// OCCT Set(theThroat, ..., Choix)
            /// (BlendFunc_ConstThroatInv.cxx L37-74).
            pub fn set(&mut self, the_throat: f64, _unused: f64, choix: i32) {
                self.throat = the_throat;

                self.choix = choix;
                match self.choix {
                    1 | 2 => {
                        self.sign1 = -1.0;
                        self.sign2 = -1.0;
                    }
                    3 | 4 => {
                        self.sign1 = 1.0;
                        self.sign2 = -1.0;
                    }
                    5 | 6 => {
                        self.sign1 = 1.0;
                        self.sign2 = 1.0;
                    }
                    7 | 8 => {
                        self.sign1 = -1.0;
                        self.sign2 = 1.0;
                    }
                    _ => {
                        self.sign1 = -1.0;
                        self.sign2 = -1.0;
                    }
                }
            }
        }
    };
}

const_throat_inv_common!(BlendFuncConstThroatInv);
const_throat_inv_common!(BlendFuncConstThroatWithPenetrationInv);

impl<'a> BlendFuncConstThroatInv<'a> {
    /// OCCT IsSolution(Sol, Tol) (BlendFunc_ConstThroatInv.cxx L78-85).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 4];
        self.value(sol, &mut valsol);

        valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= tol * tol
            && valsol[3].abs() <= tol * tol
    }

    /// OCCT Value(X, F) (BlendFunc_ConstThroatInv.cxx L89-135).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let _v2d = csurf.derivative_at(x[0]);
        self.param = x[1];
        // OCCT: curv->D2(param, ptgui, d1gui, d2gui);
        self.ptgui = self.curv.point_at(self.param);
        self.d1gui = self.curv.derivative_at(self.param);
        self.d2gui = self.curv.derivative2_at(self.param);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -self.nplan.dot(self.ptgui);

        // OCCT: math_Vector XX(1, 4);
        let xx: [f64; 4];

        if self.first {
            xx = [p2d.x, p2d.y, x[2], x[3]];
        } else {
            xx = [x[2], x[3], p2d.x, p2d.y];
        }

        // OCCT: surf1->D0(XX(1), XX(2), pts1);
        self.pts1 = self.surf1.point_at(xx[0], xx[1]);
        // OCCT: surf2->D0(XX(3), XX(4), pts2);
        self.pts2 = self.surf2.point_at(xx[2], xx[3]);

        f[0] = self.nplan.dot(self.pts1) + self.the_d;
        f[1] = self.nplan.dot(self.pts2) + self.the_d;

        // OCCT: const gp_Pnt ptmid((pts1.XYZ() + pts2.XYZ()) / 2);
        //       const gp_Vec vmid(ptgui, ptmid);
        let ptmid = (self.pts1 + self.pts2) / 2.0;
        let vmid = ptmid - self.ptgui;

        f[2] = vmid.length_squared() - self.throat * self.throat;

        let vref1 = self.pts1 - self.ptgui; // OCCT: gp_Vec vref1(ptgui, pts1)
        let vref2 = self.pts2 - self.ptgui; // OCCT: gp_Vec vref2(ptgui, pts2)

        f[3] = vref1.length_squared() - vref2.length_squared();

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_ConstThroatInv.cxx L139-237).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let v2d = csurf.derivative_at(x[0]);
        self.param = x[1];
        // OCCT: curv->D2(param, ptgui, d1gui, d2gui);
        self.ptgui = self.curv.point_at(self.param);
        self.d1gui = self.curv.derivative_at(self.param);
        self.d2gui = self.curv.derivative2_at(self.param);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -self.nplan.dot(self.ptgui);

        // OCCT: dnplan.SetLinearForm(1./normtg, d2gui,
        //                           -1./normtg * (nplan.Dot(d2gui)), nplan);
        let dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        let temp1 = self.pts1 - self.ptgui;
        let temp2 = self.pts2 - self.ptgui;
        let tempmid = (self.pts1 + self.pts2) / 2.0 - self.ptgui;

        let xx: [f64; 4];
        if self.first {
            xx = [p2d.x, p2d.y, x[2], x[3]];
        } else {
            xx = [x[2], x[3], p2d.x, p2d.y];
        }

        // OCCT: surf1->D1(XX(1), XX(2), pts1, d1u1, d1v1);
        let (_, du1, dv1) = self.surf1.derivatives(xx[0], xx[1]);
        self.d1u1 = du1;
        self.d1v1 = dv1;
        // OCCT: surf2->D1(XX(3), XX(4), pts2, d1u2, d1v2);
        let (_, du2, dv2) = self.surf2.derivatives(xx[2], xx[3]);
        self.d1u2 = du2;
        self.d1v2 = dv2;

        if self.first {
            // p2d = pts est sur surf1
            // OCCT: temp.SetLinearForm(v2d.X(), d1u1, v2d.Y(), d1v1);
            let temp = v2d.x * self.d1u1 + v2d.y * self.d1v1;

            d[0][0] = self.nplan.dot(temp);
            d[1][0] = 0.0;
            d[2][0] = (self.pts1 - self.ptgui).dot(temp); // OCCT: gp_Vec(ptgui, pts1).Dot(temp)
            d[3][0] = 2.0 * (self.pts1 - self.ptgui).dot(temp);

            d[0][2] = 0.0;
            d[0][3] = 0.0;
            d[1][2] = self.nplan.dot(self.d1u2);
            d[1][3] = self.nplan.dot(self.d1v2);
            d[2][2] = ((self.pts1 + self.pts2) / 2.0 - self.ptgui).dot(self.d1u2);
            d[2][3] = ((self.pts1 + self.pts2) / 2.0 - self.ptgui).dot(self.d1v2);
            d[3][2] = -2.0 * (self.pts2 - self.ptgui).dot(self.d1u2);
            d[3][3] = -2.0 * (self.pts2 - self.ptgui).dot(self.d1v2);
        } else {
            //  p2d = pts est sur surf2
            // OCCT: temp.SetLinearForm(v2d.X(), d1u2, v2d.Y(), d1v2);
            let temp = v2d.x * self.d1u2 + v2d.y * self.d1v2;

            d[0][0] = 0.0;
            d[1][0] = self.nplan.dot(temp);
            d[2][0] = (self.pts2 - self.ptgui).dot(temp); // OCCT: gp_Vec(ptgui, pts2).Dot(temp)
            d[3][0] = -2.0 * (self.pts2 - self.ptgui).dot(temp);

            d[0][2] = self.nplan.dot(self.d1u1);
            d[0][3] = self.nplan.dot(self.d1v1);
            d[1][2] = 0.0;
            d[1][3] = 0.0;
            d[2][2] = ((self.pts1 + self.pts2) / 2.0 - self.ptgui).dot(self.d1u1);
            d[2][3] = ((self.pts1 + self.pts2) / 2.0 - self.ptgui).dot(self.d1v1);
            d[3][2] = 2.0 * (self.pts1 - self.ptgui).dot(self.d1u1);
            d[3][3] = 2.0 * (self.pts1 - self.ptgui).dot(self.d1v1);
        }

        d[0][1] = dnplan.dot(temp1) - self.nplan.dot(self.d1gui);
        d[1][1] = dnplan.dot(temp2) - self.nplan.dot(self.d1gui);
        d[2][1] = -2.0 * self.d1gui.dot(tempmid);
        d[3][1] = 2.0 * self.d1gui.dot(temp1) - 2.0 * self.d1gui.dot(temp2);

        true
    }
}

gen_chamf_inv_common!(BlendFuncConstThroatInv);
impl_gen_chamf_inv_traits!(BlendFuncConstThroatInv);

impl<'a> BlendFuncConstThroatWithPenetrationInv<'a> {
    /// OCCT IsSolution(Sol, Tol)
    /// (BlendFunc_ConstThroatWithPenetrationInv.cxx L31-42).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 4];
        self.value(sol, &mut valsol);

        valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= tol * tol
            && valsol[3].abs() <= tol
    }

    /// OCCT Value(X, F)
    /// (BlendFunc_ConstThroatWithPenetrationInv.cxx L46-97).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let _v2d = csurf.derivative_at(x[0]);
        self.param = x[1];
        // OCCT: curv->D2(param, ptgui, d1gui, d2gui);
        self.ptgui = self.curv.point_at(self.param);
        self.d1gui = self.curv.derivative_at(self.param);
        self.d2gui = self.curv.derivative2_at(self.param);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -self.nplan.dot(self.ptgui);

        let xx: [f64; 4];
        if self.first {
            xx = [p2d.x, p2d.y, x[2], x[3]];
        } else {
            xx = [x[2], x[3], p2d.x, p2d.y];
        }

        // OCCT: surf1->D0(XX(1), XX(2), pts1);
        self.pts1 = self.surf1.point_at(xx[0], xx[1]);
        // OCCT: surf2->D0(XX(3), XX(4), pts2);
        self.pts2 = self.surf2.point_at(xx[2], xx[3]);

        f[0] = self.nplan.dot(self.pts1) + self.the_d;
        f[1] = self.nplan.dot(self.pts2) + self.the_d;

        let vref = self.pts1 - self.ptgui; // OCCT: gp_Vec vref(ptgui, pts1)

        f[2] = vref.length_squared() - self.throat * self.throat;

        let vec12 = self.pts2 - self.pts1; // OCCT: gp_Vec vec12(pts1, pts2)

        f[3] = vref.dot(vec12);

        true
    }

    /// OCCT Derivatives(X, D)
    /// (BlendFunc_ConstThroatWithPenetrationInv.cxx L101-199).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let v2d = csurf.derivative_at(x[0]);
        self.param = x[1];
        // OCCT: curv->D2(param, ptgui, d1gui, d2gui);
        self.ptgui = self.curv.point_at(self.param);
        self.d1gui = self.curv.derivative_at(self.param);
        self.d2gui = self.curv.derivative2_at(self.param);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -self.nplan.dot(self.ptgui);

        // OCCT: dnplan.SetLinearForm(1./normtg, d2gui,
        //                           -1./normtg * (nplan.Dot(d2gui)), nplan);
        let dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        let temp1 = self.pts1 - self.ptgui;
        let temp2 = self.pts2 - self.ptgui;
        let temp3 = self.pts2 - self.pts1;

        let xx: [f64; 4];
        if self.first {
            xx = [p2d.x, p2d.y, x[2], x[3]];
        } else {
            xx = [x[2], x[3], p2d.x, p2d.y];
        }

        // OCCT: surf1->D1(XX(1), XX(2), pts1, d1u1, d1v1);
        let (_, du1, dv1) = self.surf1.derivatives(xx[0], xx[1]);
        self.d1u1 = du1;
        self.d1v1 = dv1;
        // OCCT: surf2->D1(XX(3), XX(4), pts2, d1u2, d1v2);
        let (_, du2, dv2) = self.surf2.derivatives(xx[2], xx[3]);
        self.d1u2 = du2;
        self.d1v2 = dv2;

        if self.first {
            // p2d = pts est sur surf1
            // OCCT: temp.SetLinearForm(v2d.X(), d1u1, v2d.Y(), d1v1);
            let temp = v2d.x * self.d1u1 + v2d.y * self.d1v1;

            d[0][0] = self.nplan.dot(temp);
            d[1][0] = 0.0;
            // OCCT: D(3, 1) = 2 * temp1.Dot(temp);
            d[2][0] = 2.0 * temp1.dot(temp);
            // OCCT: D(4, 1) = temp.Dot(temp3) - temp.Dot(temp1);
            d[3][0] = temp.dot(temp3) - temp.dot(temp1);

            d[0][2] = 0.0;
            d[0][3] = 0.0;
            d[1][2] = self.nplan.dot(self.d1u2);
            d[1][3] = self.nplan.dot(self.d1v2);
            d[2][2] = 0.0;
            d[2][3] = 0.0;
            // OCCT: D(4, 3) = temp1.Dot(d1u2);
            d[3][2] = temp1.dot(self.d1u2);
            // OCCT: D(4, 4) = temp1.Dot(d1v2);
            d[3][3] = temp1.dot(self.d1v2);
        } else {
            //  p2d = pts est sur surf2
            // OCCT: temp.SetLinearForm(v2d.X(), d1u2, v2d.Y(), d1v2);
            let temp = v2d.x * self.d1u2 + v2d.y * self.d1v2;

            d[0][0] = 0.0;
            d[1][0] = self.nplan.dot(temp);
            d[2][0] = 0.0;
            // OCCT: D(4, 1) = temp1.Dot(temp);
            d[3][0] = temp1.dot(temp);

            d[0][2] = self.nplan.dot(self.d1u1);
            d[0][3] = self.nplan.dot(self.d1v1);
            d[1][2] = 0.0;
            d[1][3] = 0.0;
            // OCCT: D(3, 3) = 2. * temp1.Dot(d1u1);
            d[2][2] = 2.0 * temp1.dot(self.d1u1);
            // OCCT: D(3, 4) = 2. * temp1.Dot(d1v1);
            d[2][3] = 2.0 * temp1.dot(self.d1v1);
            // OCCT: D(4, 3) = d1u1.Dot(temp3) - d1u1.Dot(temp1);
            d[3][2] = self.d1u1.dot(temp3) - self.d1u1.dot(temp1);
            d[3][3] = self.d1v1.dot(temp3) - self.d1v1.dot(temp1);
        }

        // OCCT: D(1, 2) = dnplan.Dot(temp1) - nplan.Dot(d1gui);
        d[0][1] = dnplan.dot(temp1) - self.nplan.dot(self.d1gui);
        d[1][1] = dnplan.dot(temp2) - self.nplan.dot(self.d1gui);
        // OCCT: D(3, 2) = -2. * d1gui.Dot(temp1);
        d[2][1] = -2.0 * self.d1gui.dot(temp1);
        // OCCT: D(4, 2) = -d1gui.Dot(temp3);
        d[3][1] = -self.d1gui.dot(temp3);

        true
    }
}

gen_chamf_inv_common!(BlendFuncConstThroatWithPenetrationInv);
impl_gen_chamf_inv_traits!(BlendFuncConstThroatWithPenetrationInv);
