//! OCCT BlendFunc_Chamfer family (TKFillet/BlendFunc) — the chamfer section
//! functions and their restriction inverses:
//! GenChamfer (BlendFunc_GenChamfer.hxx L24-162, BlendFunc_GenChamfer.cxx
//! L29-300), Corde (BlendFunc_Corde.hxx L34-100, BlendFunc_Corde.cxx
//! L30-200), Chamfer (BlendFunc_Chamfer.hxx L33-94, BlendFunc_Chamfer.cxx
//! L29-245), ChamfInv (BlendFunc_ChamfInv.hxx, BlendFunc_ChamfInv.cxx
//! L24-230) and the abstract inverse GenChamfInv (BlendFunc_GenChamfInv.hxx
//! L25-59, BlendFunc_GenChamfInv.cxx L23-123).  ChAsym lives in
//! [`super::brep_blend_func_chamfer_b`], ConstThroat / ConstThroatWithPenetration /
//! ConstThroatInv / ConstThroatWithPenetrationInv in
//! [`super::brep_blend_func_chamfer_c`].
//!
//! Architecture mappings: `class BlendFunc_GenChamfer : public Blend_Function`
//! maps to the Rust subtrait [`BlendFuncGenChamfer`] over
//! [`BlendFunction`]; `class BlendFunc_GenChamfInv : public Blend_FuncInv`
//! maps to [`BlendFuncGenChamfInv`] over [`BlendFuncInv`].  OCCT
//! implementation inheritance (the concrete method bodies of GenChamfer /
//! GenChamfInv shared by their subclasses) is expressed by the
//! `gen_chamfer_common!` / `impl_gen_chamfer_traits!` /
//! `gen_chamf_inv_common!` / `impl_gen_chamf_inv_traits!` macros which stamp
//! the single OCCT body verbatim onto every concrete class.  The four OCCT
//! `Section` overloads are named `section` (the obsolete
//! `Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)`), `section_d1`,
//! `section_d2` and `section_simple` (Section(P, Poles, Poles2d, Weigths)).

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Curve3, CurveEval as _, Curve2d, Curve2dEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::gp::Lin;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{GeomAbsShape, MatD, VecD};

use super::brep_blend_func::blend_func_next_shape;
use super::brep_blend_func_inv::BlendFuncInv;
use super::brep_blend_function::{BlendAppFunction, BlendFunction};
use super::brep_blend_point::BlendPoint;

/// OCCT BlendFunc_GenChamfer — abstract base class for the chamfer section
/// functions (BlendFunc_GenChamfer.hxx L24).
pub trait BlendFuncGenChamfer<'a>: BlendFunction {
    /// OCCT Set(Dist1, Dist2, Choix) (BlendFunc_GenChamfer.hxx L72) — sets
    /// the distances and the "quadrant"; pure virtual in OCCT.
    fn set(&mut self, dist1: f64, dist2: f64, choix: i32);

    /// OCCT Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)
    /// (BlendFunc_GenChamfer.hxx L113) — obsolete method kept for the
    /// Simul sections; implemented once in GenChamfer.cxx L103-121.
    #[allow(clippy::too_many_arguments)]
    fn section(
        &mut self,
        param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pdeb: &mut f64,
        pfin: &mut f64,
        c: &mut Lin,
    );

    /// OCCT Set(First, Last) (BlendFunc_GenChamfer.cxx L50) — the GenChamfer
    /// implementation is an empty body, shared by every subclass.
    fn set_interval(&mut self, _first: f64, _last: f64) {}
}

/// OCCT BlendFunc_GenChamfInv — abstract base class for the chamfer inverse
/// functions (BlendFunc_GenChamfInv.hxx L25).
pub trait BlendFuncGenChamfInv<'a>: BlendFuncInv {
    /// OCCT Set(Dist1, Dist2, Choix) (BlendFunc_GenChamfInv.hxx L50) — pure
    /// virtual in OCCT.
    fn set(&mut self, dist1: f64, dist2: f64, choix: i32);
}

/// OCCT BlendFunc_Corde — function computing a point on a surface at a given
/// distance from the guide curve, in the normal plane
/// (BlendFunc_Corde.hxx L29).  X(1), X(2) are the parameters U, V of pts on
/// surf.
pub struct BlendFuncCorde<'a> {
    // OCCT BlendFunc_Corde.hxx fields (L85-99).
    surf: &'a Surface3,
    guide: &'a Curve3,
    pts: DVec3,
    // OCCT BlendFunc_Corde.hxx L89 — declared in OCCT but never accessed by
    // the Corde bodies (vestigial field, kept for structure parity).
    #[allow(dead_code)]
    pt2d: DVec2,
    dis: f64,
    normtg: f64,
    the_d: f64,
    ptgui: DVec3,
    nplan: DVec3,
    d1gui: DVec3,
    d2gui: DVec3,
    tgs: DVec3,
    tg2d: DVec2,
    istangent: bool,
}

impl<'a> BlendFuncCorde<'a> {
    /// OCCT BlendFunc_Corde(S, CG) (BlendFunc_Corde.cxx L30-39).
    pub fn new(s: &'a Surface3, cg: &'a Curve3) -> Self {
        BlendFuncCorde {
            surf: s,
            guide: cg,
            pts: DVec3::ZERO,
            pt2d: DVec2::ZERO,
            dis: 0.0,
            normtg: 0.0,
            the_d: 0.0,
            ptgui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            tgs: DVec3::ZERO,
            tg2d: DVec2::ZERO,
            istangent: false,
        }
    }

    /// OCCT SetDist(Dist) (BlendFunc_Corde.cxx L43-46).
    pub fn set_dist(&mut self, dist: f64) {
        self.dis = dist;
    }

    /// OCCT SetParam(Param) (BlendFunc_Corde.cxx L50-56).
    pub fn set_param(&mut self, param: f64) {
        // OCCT: guide->D2(Param, ptgui, d1gui, d2gui);
        self.ptgui = self.guide.point_at(param);
        self.d1gui = self.guide.derivative_at(param);
        self.d2gui = self.guide.derivative2_at(param);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -self.nplan.dot(self.ptgui);
    }

    /// OCCT Value(X, F) (BlendFunc_Corde.cxx L63-73) — returns F(U, V).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: surf->D1(X(1), X(2), pts, d1u, d1v);
        // (d1u / d1v are computed by OCCT but not used in F.)
        let (_, _d1u, _d1v) = self.surf.derivatives(x[0], x[1]);
        self.pts = self.surf.point_at(x[0], x[1]);

        f[0] = self.nplan.dot(self.pts) + self.the_d;
        // OCCT: const gp_Vec vref(ptgui, pts);
        let vref = self.pts - self.ptgui;
        f[1] = vref.length_squared() - self.dis * self.dis;

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_Corde.cxx L80-91) — D = grad F(U, V).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: surf->D1(X(1), X(2), pts, d1u, d1v);
        let (_, d1u, d1v) = self.surf.derivatives(x[0], x[1]);
        self.pts = self.surf.point_at(x[0], x[1]);

        // OCCT: D(2, 1) = 2. * gp_Vec(ptgui, pts).Dot(d1u);
        let ptp = self.pts - self.ptgui;
        d[0][0] = self.nplan.dot(d1u);
        d[0][1] = self.nplan.dot(d1v);
        d[1][0] = 2.0 * ptp.dot(d1u);
        d[1][1] = 2.0 * ptp.dot(d1v);

        true
    }

    /// OCCT PointOnS() (BlendFunc_Corde.cxx L95-98).
    pub fn point_on_s(&self) -> DVec3 {
        self.pts
    }

    /// OCCT PointOnGuide() (BlendFunc_Corde.cxx L102-105) — the point of
    /// parameter Param on CGuide.
    pub fn point_on_guide(&self) -> DVec3 {
        self.ptgui
    }

    /// OCCT NPlan() (BlendFunc_Corde.cxx L109-112) — the normal to CGuide at
    /// Ptgui.
    pub fn n_plan(&self) -> DVec3 {
        self.nplan
    }

    /// OCCT IsTangencyPoint() (BlendFunc_Corde.cxx L116-119).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnS() (BlendFunc_Corde.cxx L123-131).
    pub fn tangent_on_s(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_Corde::TangentOnS");
        }
        self.tgs
    }

    /// OCCT Tangent2dOnS() (BlendFunc_Corde.cxx L134-141).
    pub fn tangent_2d_on_s(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_Corde::Tangent2dOnS");
        }
        self.tg2d
    }

    /// OCCT DerFguide(Sol, DerF) (BlendFunc_Corde.cxx L145-157) — derivative
    /// of the function compared to the parameter of the guideline.
    pub fn der_fguide(&mut self, sol: &[f64], der_f: &mut DVec2) {
        // OCCT: surf->D1(Sol(1), Sol(2), pts, d1u, d1v);
        // (d1u / d1v are computed by OCCT but not used in DerFguide.)
        let (_, _d1u, _d1v) = self.surf.derivatives(sol[0], sol[1]);
        self.pts = self.surf.point_at(sol[0], sol[1]);

        // OCCT: dnplan.SetLinearForm(1./normtg, d2gui,
        //                           -1./normtg * (nplan.Dot(d2gui)), nplan);
        let dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        let temp = self.pts - self.ptgui;

        der_f.x = dnplan.dot(temp) - self.nplan.dot(self.d1gui);
        der_f.y = -2.0 * self.d1gui.dot(temp);
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_Corde.cxx L161-200) — returns
    /// false if Sol is not solution else returns true and updates the fields
    /// tgs and tg2d.
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector secmember(1, 2), valsol(1, 2);
        //       math_Matrix gradsol(1, 2, 1, 2);
        let mut valsol = [0.0f64; 2];
        let mut secmember = [0.0f64; 2];
        let mut gradsol = vec![vec![0.0f64; 2]; 2];

        self.value(sol, &mut valsol);
        self.derivatives(sol, &mut gradsol);
        if valsol[0].abs() <= tol && valsol[1].abs() <= tol * tol {
            // OCCT: surf->D1(Sol(1), Sol(2), pts, d1u, d1v);
            let (_, d1u, d1v) = self.surf.derivatives(sol[0], sol[1]);
            self.pts = self.surf.point_at(sol[0], sol[1]);
            // OCCT: dnplan.SetLinearForm(1./normtg, d2gui,
            //                           -1./normtg * (nplan.Dot(d2gui)), nplan);
            let dnplan = (1.0 / self.normtg) * self.d2gui
                + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            let temp = self.pts - self.ptgui;

            secmember[0] = self.nplan.dot(self.d1gui) - dnplan.dot(temp);
            secmember[1] = 2.0 * self.d1gui.dot(temp);

            //  gradsol*der = secmember
            //  with  der(1) = dU/dW, der(2) = dU/dW, W is the guide parameter

            // OCCT: math_Gauss Resol(gradsol);
            let mut a = MatD::new(2, 2);
            for r in 1..=2 {
                for c in 1..=2 {
                    a.set(r, c, gradsol[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                let mut x = VecD::new(2);
                for i in 1..=2 {
                    x.set(i, secmember[i - 1]);
                }
                resol.solve(&mut x);
                for i in 1..=2 {
                    secmember[i - 1] = x.get(i);
                }
                // OCCT: tgs.SetLinearForm(secmember(1), d1u, secmember(2), d1v);
                self.tgs = secmember[0] * d1u + secmember[1] * d1v;
                self.tg2d = DVec2::new(secmember[0], secmember[1]);
                self.istangent = false;
            } else {
                self.istangent = true;
            }
            return true;
        }

        false
    }
}

/// OCCT BlendFunc_Chamfer — function for a symmetric chamfer: the distances
/// from spine to surfaces are constant (BlendFunc_Chamfer.hxx L31).
pub struct BlendFuncChamfer<'a> {
    // OCCT BlendFunc_GenChamfer.hxx fields (L154-159).
    surf1: &'a Surface3,
    surf2: &'a Surface3,
    curv: &'a Curve3,
    choix: i32,
    tol: f64,
    distmin: f64,
    // OCCT BlendFunc_Chamfer.hxx fields (L92-93).
    corde1: BlendFuncCorde<'a>,
    corde2: BlendFuncCorde<'a>,
}

/// OCCT BlendFunc_ChAsym — function for an asymmetric chamfer
/// (BlendFunc_ChAsym.hxx L34, fields L198-220).  Methods live in
/// [`super::brep_blend_func_chamfer_b`].
pub struct BlendFuncChAsym<'a> {
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> tcurv — a handle aliasing curv until
    // Set(First, Last) trims it.  The rcad port owns the trimmed copy
    // (architecture mapping: Handle copy -> owned enum + accessor).
    pub(crate) tcurv: Option<Curve3>,
    pub(crate) param: f64,
    pub(crate) dist1: f64,
    pub(crate) angle: f64,
    pub(crate) tgang: f64,
    pub(crate) nplan: DVec3,
    pub(crate) pt1: DVec3,
    pub(crate) tsurf1: DVec3,
    pub(crate) pt2: DVec3,
    pub(crate) fx: [f64; 4],
    pub(crate) dx: Vec<Vec<f64>>,
    pub(crate) istangent: bool,
    pub(crate) tg1: DVec3,
    pub(crate) tg12d: DVec2,
    pub(crate) tg2: DVec3,
    pub(crate) tg22d: DVec2,
    pub(crate) choix: i32,
    pub(crate) distmin: f64,
}

/// OCCT BlendFunc_ConstThroat — symmetric chamfer with constant throat
/// (BlendFunc_ConstThroat.hxx, fields L85-112).  Methods live in
/// [`super::brep_blend_func_chamfer_c`].
pub struct BlendFuncConstThroat<'a> {
    // OCCT BlendFunc_GenChamfer.hxx fields (L154-159).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) choix: i32,
    pub(crate) tol: f64,
    pub(crate) distmin: f64,
    // OCCT BlendFunc_ConstThroat.hxx fields (L85-112).
    pub(crate) pts1: DVec3,
    pub(crate) pts2: DVec3,
    pub(crate) d1u1: DVec3,
    pub(crate) d1v1: DVec3,
    pub(crate) d1u2: DVec3,
    pub(crate) d1v2: DVec3,
    pub(crate) istangent: bool,
    pub(crate) tg1: DVec3,
    pub(crate) tg12d: DVec2,
    pub(crate) tg2: DVec3,
    pub(crate) tg22d: DVec2,
    pub(crate) param: f64,
    pub(crate) throat: f64,
    pub(crate) ptgui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
}

/// OCCT BlendFunc_ConstThroatWithPenetration — chamfer with constant throat
/// and penetration (BlendFunc_ConstThroatWithPenetration.hxx).  It derives
/// from BlendFunc_ConstThroat in OCCT, so it carries the same fields.
/// Methods live in [`super::brep_blend_func_chamfer_c`].
pub struct BlendFuncConstThroatWithPenetration<'a> {
    // OCCT BlendFunc_GenChamfer.hxx fields (L154-159).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) choix: i32,
    pub(crate) tol: f64,
    pub(crate) distmin: f64,
    // OCCT BlendFunc_ConstThroat.hxx fields (L85-112), inherited.
    pub(crate) pts1: DVec3,
    pub(crate) pts2: DVec3,
    pub(crate) d1u1: DVec3,
    pub(crate) d1v1: DVec3,
    pub(crate) d1u2: DVec3,
    pub(crate) d1v2: DVec3,
    pub(crate) istangent: bool,
    pub(crate) tg1: DVec3,
    pub(crate) tg12d: DVec2,
    pub(crate) tg2: DVec3,
    pub(crate) tg22d: DVec2,
    pub(crate) param: f64,
    pub(crate) throat: f64,
    pub(crate) ptgui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
}

/// OCCT BlendFunc_ChamfInv — inverse of the symmetric chamfer function
/// (BlendFunc_ChamfInv.hxx L24).  The vector X is t, w, U, V.
pub struct BlendFuncChamfInv<'a> {
    // OCCT BlendFunc_GenChamfInv.hxx fields (L53-58).
    surf1: &'a Surface3,
    surf2: &'a Surface3,
    curv: &'a Curve3,
    csurf: Option<&'a Curve2d>,
    choix: i32,
    first: bool,
    // OCCT BlendFunc_ChamfInv.hxx fields: corde1, corde2.
    corde1: BlendFuncCorde<'a>,
    corde2: BlendFuncCorde<'a>,
}

/// OCCT BlendFunc_ChAsymInv — inverse of the asymmetric chamfer function
/// (BlendFunc_ChAsymInv.hxx L26, fields L74-86).  Methods live in
/// [`super::brep_blend_func_chamfer_b`].
pub struct BlendFuncChAsymInv<'a> {
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) dist1: f64,
    pub(crate) angle: f64,
    pub(crate) tgang: f64,
    pub(crate) curv: &'a Curve3,
    pub(crate) csurf: Option<&'a Curve2d>,
    pub(crate) choix: i32,
    pub(crate) first: bool,
    pub(crate) fx: [f64; 4],
    pub(crate) dx: Vec<Vec<f64>>,
}

/// OCCT BlendFunc_ConstThroatInv — inverse of the constant throat function
/// (BlendFunc_ConstThroatInv.hxx, fields L52-88).  Methods live in
/// [`super::brep_blend_func_chamfer_c`].
pub struct BlendFuncConstThroatInv<'a> {
    // OCCT BlendFunc_GenChamfInv.hxx fields (L53-58).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) csurf: Option<&'a Curve2d>,
    pub(crate) choix: i32,
    pub(crate) first: bool,
    // OCCT BlendFunc_ConstThroatInv.hxx protected fields (L52-88).
    pub(crate) throat: f64,
    pub(crate) param: f64,
    pub(crate) sign1: f64,
    pub(crate) sign2: f64,
    pub(crate) ptgui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) pts1: DVec3,
    pub(crate) pts2: DVec3,
    pub(crate) d1u1: DVec3,
    pub(crate) d1v1: DVec3,
    pub(crate) d1u2: DVec3,
    pub(crate) d1v2: DVec3,
}

/// OCCT BlendFunc_ConstThroatWithPenetrationInv — inverse of the constant
/// throat with penetration function
/// (BlendFunc_ConstThroatWithPenetrationInv.hxx).  It derives from
/// BlendFunc_ConstThroatInv in OCCT, so it carries the same fields.
/// Methods live in [`super::brep_blend_func_chamfer_c`].
pub struct BlendFuncConstThroatWithPenetrationInv<'a> {
    // OCCT BlendFunc_GenChamfInv.hxx fields (L53-58).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    pub(crate) csurf: Option<&'a Curve2d>,
    pub(crate) choix: i32,
    pub(crate) first: bool,
    // OCCT BlendFunc_ConstThroatInv.hxx protected fields, inherited.
    pub(crate) throat: f64,
    pub(crate) param: f64,
    pub(crate) sign1: f64,
    pub(crate) sign2: f64,
    pub(crate) ptgui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) pts1: DVec3,
    pub(crate) pts2: DVec3,
    pub(crate) d1u1: DVec3,
    pub(crate) d1v1: DVec3,
    pub(crate) d1u2: DVec3,
    pub(crate) d1v2: DVec3,
}

impl<'a> BlendFuncChamfer<'a> {
    /// OCCT BlendFunc_Chamfer(S1, S2, CG) (BlendFunc_Chamfer.cxx L29-36) and
    /// the inherited BlendFunc_GenChamfer ctor (GenChamfer.cxx L29-39).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, cg: &'a Curve3) -> Self {
        BlendFuncChamfer {
            surf1: s1,
            surf2: s2,
            curv: cg,
            choix: 0,
            tol: 0.0,
            distmin: f64::MAX, // OCCT: RealLast()
            corde1: BlendFuncCorde::new(s1, cg),
            corde2: BlendFuncCorde::new(s2, cg),
        }
    }
}

impl<'a> BlendFuncChamfInv<'a> {
    /// OCCT BlendFunc_ChamfInv(S1, S2, C) (BlendFunc_ChamfInv.cxx L24-31) and
    /// the inherited BlendFunc_GenChamfInv ctor (GenChamfInv.cxx L23-32).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncChamfInv {
            surf1: s1,
            surf2: s2,
            curv: c,
            csurf: None,
            choix: 0,
            first: false,
            corde1: BlendFuncCorde::new(s1, c),
            corde2: BlendFuncCorde::new(s2, c),
        }
    }
}

/// Stamps the concrete methods inherited from BlendFunc_GenChamfer
/// (BlendFunc_GenChamfer.cxx L43-300) onto a subclass.  OCCT expresses the
/// sharing through C++ inheritance; the Rust port has no implementation
/// inheritance, so the single OCCT body is stamped verbatim per class.
#[allow(clippy::too_many_arguments)]
macro_rules! gen_chamfer_common {
    ($t:ident) => {
        impl<'a> $t<'a> {
            /// OCCT NbEquations() (BlendFunc_GenChamfer.cxx L43-46) —
            /// returns 4.
            pub fn nb_equations(&self) -> usize {
                4
            }

            /// OCCT GetTolerance(Tolerance, Tol)
            /// (BlendFunc_GenChamfer.cxx L54-60).
            pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
                tolerance[0] = self.surf1.u_resolution(tol);
                tolerance[1] = self.surf1.v_resolution(tol);
                tolerance[2] = self.surf2.u_resolution(tol);
                tolerance[3] = self.surf2.v_resolution(tol);
            }

            /// OCCT GetBounds(InfBound, SupBound)
            /// (BlendFunc_GenChamfer.cxx L64-84).
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

            /// OCCT GetMinimalDistance()
            /// (BlendFunc_GenChamfer.cxx L88-91).
            pub fn get_minimal_distance(&self) -> f64 {
                self.distmin
            }

            /// OCCT Values(X, F, D) (BlendFunc_GenChamfer.cxx L95-99).
            pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
                let val = self.value(x, f);
                val && self.derivatives(x, d)
            }

            /// OCCT Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)
            /// (BlendFunc_GenChamfer.cxx L103-121) — obsolete method.
            #[allow(clippy::too_many_arguments)]
            pub fn section(
                &mut self,
                _param: f64,
                u1: f64,
                v1: f64,
                u2: f64,
                v2: f64,
                pdeb: &mut f64,
                pfin: &mut f64,
                c: &mut Lin,
            ) {
                let pts1 = self.surf1.point_at(u1, v1);
                let pts2 = self.surf2.point_at(u2, v2);
                // OCCT: const gp_Dir dir(gp_Vec(pts1, pts2));
                let dir = (pts2 - pts1).normalize_or_zero();

                c.pos = pts1;
                c.dir = dir;

                *pdeb = 0.0;
                // OCCT: Pfin = ElCLib::Parameter(C, pts2);
                *pfin = (pts2 - c.pos).dot(c.dir);
            }

            /// OCCT IsRational() (BlendFunc_GenChamfer.cxx L125-128).
            pub fn is_rational(&self) -> bool {
                false
            }

            /// OCCT GetMinimalWeight(Weights)
            /// (BlendFunc_GenChamfer.cxx L132-135).
            pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
                for w in weigths.iter_mut() {
                    *w = 1.0;
                }
            }

            /// OCCT NbIntervals(S) (BlendFunc_GenChamfer.cxx L139-142).
            pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
                self.curv.nb_intervals(blend_func_next_shape(s))
            }

            /// OCCT Intervals(T, S) (BlendFunc_GenChamfer.cxx L146-149).
            pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
                let mut intervals = Vec::new();
                self.curv.intervals(&mut intervals, blend_func_next_shape(s));
                for (dst, src) in t.iter_mut().zip(intervals) {
                    *dst = src;
                }
            }

            /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
            /// (BlendFunc_GenChamfer.cxx L153-159).
            pub fn get_shape(
                &mut self,
                nb_poles: &mut i32,
                nb_knots: &mut i32,
                degree: &mut i32,
                nb_poles_2d: &mut i32,
            ) {
                *nb_poles = 2;
                *nb_poles_2d = 2;
                *nb_knots = 2;
                *degree = 1;
            }

            /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1D)
            /// (BlendFunc_GenChamfer.cxx L165-172).
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

            /// OCCT Knots(TKnots) (BlendFunc_GenChamfer.cxx L176-180).
            pub fn knots(&mut self, tknots: &mut [f64]) {
                tknots[0] = 0.0;
                tknots[1] = 1.0;
            }

            /// OCCT Mults(TMults) (BlendFunc_GenChamfer.cxx L184-188).
            pub fn mults(&mut self, tmults: &mut [i32]) {
                tmults[0] = 2;
                tmults[1] = 2;
            }

            /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d,
            /// D2Poles2d, Weights, DWeights, D2Weights)
            /// (BlendFunc_GenChamfer.cxx L192-204).
            #[allow(clippy::too_many_arguments)]
            pub fn section_d2(
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

            /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights,
            /// DWeights) (BlendFunc_GenChamfer.cxx L208-254) — used for the
            /// first and last section.
            #[allow(clippy::too_many_arguments)]
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
                // OCCT: math_Vector sol(1, 4), valsol(1, 4), secmember(1, 4);
                //       math_Matrix gradsol(1, 4, 1, 4);
                let mut sol = [0.0f64; 4];
                let mut valsol = [0.0f64; 4];
                let mut gradsol = vec![vec![0.0f64; 4]; 4];

                let prm = p.parameter();
                let low = 0usize; // OCCT: Poles.Lower()
                let upp = poles.len() - 1; // OCCT: Poles.Upper()

                let (u1, v1) = p.parameters_on_s1();
                sol[0] = u1;
                sol[1] = v1;
                let (u2, v2) = p.parameters_on_s2();
                sol[2] = u2;
                sol[3] = v2;

                self.set_param(prm);

                self.values(&sol, &mut valsol, &mut gradsol);
                self.is_solution(&sol, self.tol);

                let istgt = self.is_tangency_point();

                poles_2d[0] = DVec2::new(sol[0], sol[1]);
                poles_2d[poles_2d.len() - 1] = DVec2::new(sol[2], sol[3]);
                if !istgt {
                    let t2d1 = self.tangent_2d_on_s1();
                    d_poles_2d[0] = DVec2::new(t2d1.x, t2d1.y);
                    let t2d2 = self.tangent_2d_on_s2();
                    d_poles_2d[poles_2d.len() - 1] = DVec2::new(t2d2.x, t2d2.y);
                }
                poles[low] = self.point_on_s1();
                poles[upp] = self.point_on_s2();
                weigths[low] = 1.0;
                weigths[upp] = 1.0;
                if !istgt {
                    d_poles[low] = self.tangent_on_s1();
                    d_poles[upp] = self.tangent_on_s2();
                    d_weigths[low] = 0.0;
                    d_weigths[upp] = 0.0;
                }

                !istgt
            }

            /// OCCT Section(P, Poles, Poles2d, Weights)
            /// (BlendFunc_GenChamfer.cxx L258-283) — the Rust name carries the
            /// `_simple` suffix because the obsolete 8-argument overload owns
            /// the plain `section` name.
            pub fn section_simple(
                &mut self,
                p: &BlendPoint,
                poles: &mut [DVec3],
                poles_2d: &mut [DVec2],
                weigths: &mut [f64],
            ) {
                // OCCT: double u1, v1, u2, v2, prm = P.Parameter();
                let prm = p.parameter();
                let low = 0usize; // OCCT: Poles.Lower()
                let upp = poles.len() - 1; // OCCT: Poles.Upper()
                let mut x = [0.0f64; 4];
                let mut f = [0.0f64; 4];

                let (u1, v1) = p.parameters_on_s1();
                let (u2, v2) = p.parameters_on_s2();
                x[0] = u1;
                x[1] = v1;
                x[2] = u2;
                x[3] = v2;
                poles_2d[0] = DVec2::new(u1, v1);
                poles_2d[poles_2d.len() - 1] = DVec2::new(u2, v2);

                self.set_param(prm);
                self.value(&x, &mut f);
                poles[low] = self.point_on_s1();
                poles[upp] = self.point_on_s2();
                weigths[low] = 1.0;
                weigths[upp] = 1.0;
            }

            /// OCCT Resolution(IC2d, Tol, TolU, TolV)
            /// (BlendFunc_GenChamfer.cxx L285-300).
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
    };
}

/// Stamps the trait implementations shared by every BlendFunc_GenChamfer
/// subclass: `math_FunctionSetWithDerivatives`, `Blend_AppFunction`,
/// `Blend_Function` and `BlendFunc_GenChamfer`.  The concrete class must
/// provide the OCCT deferred bodies (Set(Param), Set(Dist1, Dist2, Choix),
/// Value, Derivatives, IsSolution, PointOnS1/S2, IsTangencyPoint,
/// TangentOnS1/S2, Tangent2dOnS1/S2, Tangent, GetSectionSize).
macro_rules! impl_gen_chamfer_traits {
    ($t:ident) => {
        impl<'a> FunctionSetWithDerivatives for $t<'a> {
            fn nb_variables(&self) -> usize {
                // OCCT: inherited from Blend_Function::NbVariables (returns 4).
                BlendFunction::nb_variables(self)
            }

            fn nb_equations(&self) -> usize {
                $t::nb_equations(self)
            }

            fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
                $t::value(self, x, f)
            }

            fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
                $t::derivatives(self, x, df)
            }

            fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
                $t::values(self, x, f, df)
            }
        }

        impl<'a> BlendAppFunction for $t<'a> {
            fn set_param(&mut self, param: f64) {
                $t::set_param(self, param)
            }

            // OCCT Blend_Function.cxx L24-32 — Pnt1/Pnt2 delegate to
            // PointOnS1/PointOnS2 (Blend_Function override of the AppFunction
            // deferred methods).
            fn pnt1(&self) -> DVec3 {
                BlendFunction::pnt1(self)
            }

            fn pnt2(&self) -> DVec3 {
                BlendFunction::pnt2(self)
            }

            fn set_interval(&mut self, first: f64, last: f64) {
                BlendFuncGenChamfer::set_interval(self, first, last)
            }

            fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
                $t::get_tolerance(self, tolerance, tol)
            }

            fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
                $t::get_bounds(self, inf_bound, sup_bound)
            }

            fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
                $t::is_solution(self, sol, tol)
            }

            fn get_minimal_distance(&self) -> f64 {
                $t::get_minimal_distance(self)
            }

            fn is_rational(&self) -> bool {
                $t::is_rational(self)
            }

            fn get_section_size(&self) -> f64 {
                $t::get_section_size(self)
            }

            fn get_minimal_weight(&self, weigths: &mut [f64]) {
                $t::get_minimal_weight(self, weigths)
            }

            fn nb_intervals(&self, s: GeomAbsShape) -> usize {
                $t::nb_intervals(self, s)
            }

            fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
                $t::intervals(self, t, s)
            }

            fn get_shape(
                &mut self,
                nb_poles: &mut i32,
                nb_knots: &mut i32,
                degree: &mut i32,
                nb_poles_2d: &mut i32,
            ) {
                $t::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
            }

            fn get_approx_tolerance(
                &self,
                bound_tol: f64,
                surf_tol: f64,
                angle_tol: f64,
                tol3d: &mut [f64],
                tol1d: &mut [f64],
            ) {
                $t::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
            }

            fn knots(&mut self, tknots: &mut [f64]) {
                $t::knots(self, tknots)
            }

            fn mults(&mut self, tmults: &mut [i32]) {
                $t::mults(self, tmults)
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
                $t::section_d1(
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
                $t::section_simple(self, p, poles, poles_2d, weigths)
            }

            fn section_d2(
                &mut self,
                p: &BlendPoint,
                poles: &mut [DVec3],
                d_poles: &mut [DVec3],
                d2_poles: &mut [DVec3],
                poles_2d: &mut [DVec2],
                d_poles_2d: &mut [DVec2],
                d2_poles_2d: &mut [DVec2],
                weigths: &mut [f64],
                d_weigths: &mut [f64],
                d2_weigths: &mut [f64],
            ) -> bool {
                $t::section_d2(
                    self,
                    p,
                    poles,
                    d_poles,
                    d2_poles,
                    poles_2d,
                    d_poles_2d,
                    d2_poles_2d,
                    weigths,
                    d_weigths,
                    d2_weigths,
                )
            }

            fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
                $t::resolution(self, ic_2d, tol, tol_u, tol_v)
            }
        }

        impl<'a> BlendFunction for $t<'a> {
            fn point_on_s1(&self) -> DVec3 {
                $t::point_on_s1(self)
            }

            fn point_on_s2(&self) -> DVec3 {
                $t::point_on_s2(self)
            }

            fn is_tangency_point(&self) -> bool {
                $t::is_tangency_point(self)
            }

            fn tangent_on_s1(&self) -> DVec3 {
                $t::tangent_on_s1(self)
            }

            fn tangent_2d_on_s1(&self) -> DVec2 {
                $t::tangent_2d_on_s1(self)
            }

            fn tangent_on_s2(&self) -> DVec3 {
                $t::tangent_on_s2(self)
            }

            fn tangent_2d_on_s2(&self) -> DVec2 {
                $t::tangent_2d_on_s2(self)
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
                $t::tangent(self, u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last)
            }
        }

        impl<'a> BlendFuncGenChamfer<'a> for $t<'a> {
            fn set(&mut self, dist1: f64, dist2: f64, choix: i32) {
                $t::set(self, dist1, dist2, choix)
            }

            fn section(
                &mut self,
                param: f64,
                u1: f64,
                v1: f64,
                u2: f64,
                v2: f64,
                pdeb: &mut f64,
                pfin: &mut f64,
                c: &mut Lin,
            ) {
                $t::section(self, param, u1, v1, u2, v2, pdeb, pfin, c)
            }
        }
    };
}

impl<'a> BlendFuncChamfer<'a> {
    /// OCCT Set(Dist1, Dist2, Choix) (BlendFunc_Chamfer.cxx L40-45).
    pub fn set(&mut self, dist1: f64, dist2: f64, choix: i32) {
        self.corde1.set_dist(dist1);
        self.corde2.set_dist(dist2);
        self.choix = choix;
    }

    /// OCCT Set(Param) (BlendFunc_Chamfer.cxx L49-53).
    pub fn set_param(&mut self, param: f64) {
        self.corde1.set_param(param);
        self.corde2.set_param(param);
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_Chamfer.cxx L57-75).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector Sol1(1, 2), Sol2(1, 2);
        let sol1 = [sol[0], sol[1]];
        let sol2 = [sol[2], sol[3]];

        let mut issol = self.corde1.is_solution(&sol1, tol);
        issol = issol && self.corde2.is_solution(&sol2, tol);
        self.tol = tol;
        if issol {
            self.distmin = self.distmin.min(self.corde1.point_on_s().distance(self.corde2.point_on_s()));
        }

        issol
    }

    /// OCCT Value(X, F) (BlendFunc_Chamfer.cxx L79-96).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: math_Vector x(1, 2), f(1, 2);
        let mut cx = [0.0f64; 2];
        let mut cf = [0.0f64; 2];

        cx[0] = x[0];
        cx[1] = x[1];
        self.corde1.value(&cx, &mut cf);
        f[0] = cf[0];
        f[1] = cf[1];

        cx[0] = x[2];
        cx[1] = x[3];
        self.corde2.value(&cx, &mut cf);
        f[2] = cf[0];
        f[3] = cf[1];

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_Chamfer.cxx L100-131).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: math_Vector x(1, 2); math_Matrix d(1, 2, 1, 2);
        let mut cx = [0.0f64; 2];
        let mut cd = vec![vec![0.0f64; 2]; 2];

        cx[0] = x[0];
        cx[1] = x[1];
        self.corde1.derivatives(&cx, &mut cd);
        for i in 1..3 {
            for j in 1..3 {
                d[i - 1][j - 1] = cd[i - 1][j - 1];
                d[i - 1][j + 1] = 0.0;
            }
        }

        cx[0] = x[2];
        cx[1] = x[3];
        self.corde2.derivatives(&cx, &mut cd);
        for i in 1..3 {
            for j in 1..3 {
                d[i + 1][j + 1] = cd[i - 1][j - 1];
                d[i + 1][j - 1] = 0.0;
            }
        }

        true
    }

    /// OCCT PointOnS1() (BlendFunc_Chamfer.cxx L135-138).
    pub fn point_on_s1(&self) -> DVec3 {
        self.corde1.point_on_s()
    }

    /// OCCT PointOnS2() (BlendFunc_Chamfer.cxx L142-145).
    pub fn point_on_s2(&self) -> DVec3 {
        self.corde2.point_on_s()
    }

    /// OCCT IsTangencyPoint() (BlendFunc_Chamfer.cxx L149-152).
    pub fn is_tangency_point(&self) -> bool {
        self.corde1.is_tangency_point() && self.corde2.is_tangency_point()
    }

    /// OCCT TangentOnS1() (BlendFunc_Chamfer.cxx L156-159).
    pub fn tangent_on_s1(&self) -> DVec3 {
        self.corde1.tangent_on_s()
    }

    /// OCCT TangentOnS2() (BlendFunc_Chamfer.cxx L163-166).
    pub fn tangent_on_s2(&self) -> DVec3 {
        self.corde2.tangent_on_s()
    }

    /// OCCT Tangent2dOnS1() (BlendFunc_Chamfer.cxx L170-173).
    pub fn tangent_2d_on_s1(&self) -> DVec2 {
        self.corde1.tangent_2d_on_s()
    }

    /// OCCT Tangent2dOnS2() (BlendFunc_Chamfer.cxx L177-180).
    pub fn tangent_2d_on_s2(&self) -> DVec2 {
        self.corde2.tangent_2d_on_s()
    }

    /// OCCT Tangent(U1, V1, U2, V2, TgF, TgL, NmF, NmL)
    /// (BlendFunc_Chamfer.cxx L188-236) — TgF, NmF et TgL, NmL les tangentes
    /// et normales respectives aux surfaces S1 et S2.
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

        let ptgui = self.corde1.point_on_guide();
        let nplan = self.corde1.n_plan();
        // OCCT: surf1->D1(U1, V1, pt1, d1u1, d1v1); NmF = d1u1.Crossed(d1v1);
        let (_, d1u1, d1v1) = self.surf1.derivatives(u1, v1);
        *nm_f = d1u1.cross(d1v1);

        // OCCT: surf2->D1(U2, V2, pt2, d1u2, d1v2); NmL = d1u2.Crossed(d1v2);
        let (_, d1u2, d1v2) = self.surf2.derivatives(u2, v2);
        *nm_l = d1u2.cross(d1v2);
        let _ = ptgui;

        *tg_f = nplan.cross(*nm_f).normalize_or_zero();
        *tg_l = nplan.cross(*nm_l).normalize_or_zero();

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

    /// OCCT GetSectionSize() (BlendFunc_Chamfer.cxx L242-245) — non
    /// implementee (non necessaire car non rationel).
    pub fn get_section_size(&self) -> f64 {
        panic!("Standard_NotImplemented: BlendFunc_Chamfer::GetSectionSize()");
    }
}

gen_chamfer_common!(BlendFuncChamfer);
impl_gen_chamfer_traits!(BlendFuncChamfer);

impl<'a> BlendFuncChamfInv<'a> {
    /// OCCT Set(Dist1, Dist2, Choix) (BlendFunc_ChamfInv.cxx L35-72).
    pub fn set(&mut self, dist1: f64, dist2: f64, choix: i32) {
        let dis1: f64;
        let dis2: f64;

        self.choix = choix;
        match self.choix {
            1 | 2 => {
                dis1 = -dist1;
                dis2 = -dist2;
            }
            3 | 4 => {
                dis1 = dist1;
                dis2 = -dist2;
            }
            5 | 6 => {
                dis1 = dist1;
                dis2 = dist2;
            }
            7 | 8 => {
                dis1 = -dist1;
                dis2 = dist2;
            }
            _ => {
                dis1 = -dist1;
                dis2 = -dist2;
            }
        }
        self.corde1.set_dist(dis1);
        self.corde2.set_dist(dis2);
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_ChamfInv.cxx L76-103).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: csurf->D1(Sol(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(sol[0]);
        let _v2d = csurf.derivative_at(sol[0]);

        // OCCT: math_Vector Sol1(1, 2), Sol2(1, 2);
        let sol1 = [p2d.x, p2d.y];
        let sol2 = [sol[2], sol[3]];

        let issol;
        if self.first {
            let mut r = self.corde1.is_solution(&sol1, tol);
            r = r && self.corde2.is_solution(&sol2, tol);
            issol = r;
        } else {
            let mut r = self.corde1.is_solution(&sol2, tol);
            r = r && self.corde2.is_solution(&sol1, tol);
            issol = r;
        }

        issol
    }

    /// OCCT Value(X, F) (BlendFunc_ChamfInv.cxx L107-138).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let _v2d = csurf.derivative_at(x[0]);
        self.corde1.set_param(x[1]);
        self.corde2.set_param(x[1]);

        // OCCT: math_Vector x1(1, 2), f1(1, 2), x2(1, 2), f2(1, 2);
        let x1 = [p2d.x, p2d.y];
        let mut f1 = [0.0f64; 2];
        let x2 = [x[2], x[3]];
        let mut f2 = [0.0f64; 2];

        if self.first {
            self.corde1.value(&x1, &mut f1);
            self.corde2.value(&x2, &mut f2);
        } else {
            self.corde1.value(&x2, &mut f1);
            self.corde2.value(&x1, &mut f2);
        }
        f[0] = f1[0];
        f[1] = f1[1];
        f[2] = f2[0];
        f[3] = f2[1];

        true
    }

    /// OCCT Derivatives(X, D) (BlendFunc_ChamfInv.cxx L142-230).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT: math_Matrix d1(1, 2, 1, 2), d2(1, 2, 1, 2);
        let mut d1 = vec![vec![0.0f64; 2]; 2];
        let mut d2 = vec![vec![0.0f64; 2]; 2];
        let mut df1 = DVec2::ZERO;
        let mut df2 = DVec2::ZERO;

        // OCCT: csurf->D1(X(1), p2d, v2d);
        let csurf = self.csurf.expect("csurf");
        let p2d = csurf.point_at(x[0]);
        let v2d = csurf.derivative_at(x[0]);
        self.corde1.set_param(x[1]);
        self.corde2.set_param(x[1]);

        let x1 = [p2d.x, p2d.y];
        let x2 = [x[2], x[3]];

        let ptgui: DVec3;
        let nplan: DVec3;
        let pts: DVec3;
        let d1u: DVec3;
        let d1v: DVec3;
        if self.first {
            // p2d = pts est sur surf1
            ptgui = self.corde1.point_on_guide();
            nplan = self.corde1.n_plan();
            self.corde2.derivatives(&x2, &mut d2);
            self.corde1.der_fguide(&x1, &mut df1);
            self.corde2.der_fguide(&x2, &mut df2);
            let (p, du, dv) = self.surf1.derivatives(x1[0], x1[1]);
            pts = p;
            d1u = du;
            d1v = dv;
        } else {
            //  p2d = pts est sur surf2
            ptgui = self.corde2.point_on_guide();
            nplan = self.corde2.n_plan();
            self.corde1.derivatives(&x2, &mut d1);
            self.corde1.der_fguide(&x2, &mut df1);
            self.corde2.der_fguide(&x1, &mut df2);
            let (p, du, dv) = self.surf2.derivatives(x1[0], x1[1]);
            pts = p;
            d1u = du;
            d1v = dv;
        }

        // derivees par rapport a T
        // OCCT: temp.SetLinearForm(v2d.X(), d1u, v2d.Y(), d1v);
        let temp = v2d.x * d1u + v2d.y * d1v;
        if self.first {
            d[0][0] = nplan.dot(temp);
            d[1][0] = 2.0 * (pts - ptgui).dot(temp);
            d[2][0] = 0.0;
            d[3][0] = 0.0;
        } else {
            d[0][0] = 0.0;
            d[1][0] = 0.0;
            d[2][0] = nplan.dot(temp);
            d[3][0] = 2.0 * (pts - ptgui).dot(temp);
        }

        // derivees par rapport a W
        d[0][1] = df1.x;
        d[1][1] = df1.y;
        d[2][1] = df2.x;
        d[3][1] = df2.y;

        // derivees par rapport a U et V
        if self.first {
            for i in 1..3 {
                for j in 3..5 {
                    d[i - 1][j - 1] = 0.0;
                    d[i + 1][j - 1] = d2[i - 1][j - 3];
                }
            }
        } else {
            for i in 1..3 {
                for j in 3..5 {
                    d[i - 1][j - 1] = d1[i - 1][j - 3];
                    d[i + 1][j - 1] = 0.0;
                }
            }
        }

        true
    }
}


/// GAP (plan 0.6): OCCT csurf->Resolution(Tol)
/// (Adaptor2d_Curve2d::Resolution) — pending in rcad-kernel.  The helper
/// wraps the pending marker so that the following OCCT statements stay
/// reachable for the compiler (no unreachable-code noise at the call sites).
pub(crate) fn adaptor2d_curve2d_resolution_pending() -> f64 {
    unimplemented!("Adaptor2d_Curve2d::Resolution pending in rcad-kernel")
}

/// GAP (plan 0.6): GeomFill::GetCircle etc. share the marker convention.

/// Stamps the concrete methods inherited from BlendFunc_GenChamfInv
/// (BlendFunc_GenChamfInv.cxx L36-123) onto a subclass.  OCCT expresses the
/// sharing through C++ inheritance; the single OCCT body is stamped verbatim
/// per class.
macro_rules! gen_chamf_inv_common {
    ($t:ident) => {
        impl<'a> $t<'a> {
            /// OCCT NbEquations() (BlendFunc_GenChamfInv.cxx L36-39).
            pub fn nb_equations(&self) -> usize {
                4
            }

            /// OCCT Set(OnFirst, C) (BlendFunc_GenChamfInv.cxx L43-47).
            pub fn set_curve_on_surface(&mut self, on_first: bool, c: &'a Curve2d) {
                self.first = on_first;
                self.csurf = Some(c);
            }

            /// OCCT GetTolerance(Tolerance, Tol)
            /// (BlendFunc_GenChamfInv.cxx L51-65).
            pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
                // OCCT L53: Tolerance(1) = csurf->Resolution(Tol).
                // GAP (plan 0.6): rcad-kernel has no Adaptor2d_Curve2d::
                // Resolution equivalent yet; the call panics with the pending
                // marker until the kernel exposes it (same treatment as
                // BRepBlend_CurvPointRadInv::GetTolerance).
                tolerance[0] = crate::fillet::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending();
                // OCCT L54: Tolerance(2) = curv->Resolution(Tol).
                tolerance[1] = self.curv.resolution(tol);
                if self.first {
                    tolerance[2] = self.surf2.u_resolution(tol);
                    tolerance[3] = self.surf2.v_resolution(tol);
                } else {
                    tolerance[2] = self.surf1.u_resolution(tol);
                    tolerance[3] = self.surf1.v_resolution(tol);
                }
            }

            /// OCCT GetBounds(InfBound, SupBound)
            /// (BlendFunc_GenChamfInv.cxx L69-114).
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

            /// OCCT Values(X, F, D) (BlendFunc_GenChamfInv.cxx L118-123).
            pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
                self.value(x, f);
                self.derivatives(x, d);
                true
            }
        }
    };
}

/// Stamps the trait implementations shared by every BlendFunc_GenChamfInv
/// subclass: `math_FunctionSetWithDerivatives`, `Blend_FuncInv` and
/// `BlendFunc_GenChamfInv`.
macro_rules! impl_gen_chamf_inv_traits {
    ($t:ident) => {
        impl<'a> FunctionSetWithDerivatives for $t<'a> {
            fn nb_variables(&self) -> usize {
                // OCCT: inherited from Blend_FuncInv::NbVariables (returns 4).
                BlendFuncInv::nb_variables(self)
            }

            fn nb_equations(&self) -> usize {
                $t::nb_equations(self)
            }

            fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
                $t::value(self, x, f)
            }

            fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
                $t::derivatives(self, x, df)
            }

            fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
                $t::values(self, x, f, df)
            }
        }

        impl<'a> BlendFuncInv for $t<'a> {
            fn set_curve_on_surface(&mut self, on_first: bool, c_on_surf: &Curve2d) {
                // OCCT GenChamfInv.cxx L43-47; the rcad port stores the curve
                // reference (OCCT copies the handle).  The unsized cast keeps
                // the trait signature while recovering the 'a lifetime of the
                // stored restriction — the callers own the Curve2d for 'a.
                let c: &'a Curve2d = unsafe { &*(c_on_surf as *const Curve2d) };
                $t::set_curve_on_surface(self, on_first, c)
            }

            fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
                $t::get_tolerance(self, tolerance, tol)
            }

            fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
                $t::get_bounds(self, inf_bound, sup_bound)
            }

            fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
                $t::is_solution(self, sol, tol)
            }
        }

        impl<'a> BlendFuncGenChamfInv<'a> for $t<'a> {
            fn set(&mut self, dist1: f64, dist2: f64, choix: i32) {
                $t::set(self, dist1, dist2, choix)
            }
        }
    };
}

gen_chamf_inv_common!(BlendFuncChamfInv);
impl_gen_chamf_inv_traits!(BlendFuncChamfInv);

// Macro exports for the sibling translation files (chamfer_b: ChAsym /
// ChAsymInv; chamfer_c: ConstThroat / ConstThroatWithPenetration /
// ConstThroatInv / ConstThroatWithPenetrationInv).
pub(crate) use gen_chamf_inv_common;
pub(crate) use gen_chamfer_common;
pub(crate) use impl_gen_chamf_inv_traits;
pub(crate) use impl_gen_chamfer_traits;
