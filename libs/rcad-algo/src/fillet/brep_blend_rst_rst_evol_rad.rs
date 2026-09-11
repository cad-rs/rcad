//! OCCT BRepBlend_RstRstEvolRad (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_RstRstEvolRad.hxx (L40-268) + BRepBlend_RstRstEvolRad.cxx
//! (whole file L40-1092).  Function to approximate by AppSurface for
//! Edge/Edge with evolutif radius (hxx L40-41); the SimulSurf / PerformSurf
//! rst/rst variable-radius arms instantiate it as `func`
//! (ChFi3d_FilBuilder.cxx L1386 / L2189).
//!
//! Architecture mappings (mirroring [`super::brep_blend_rst_rst_const_rad`]):
//! `class BRepBlend_RstRstEvolRad : public Blend_RstRstFunction` is
//! expressed by implementing the [`BlendRstRstFunction`] and
//! [`BlendAppFunction`] traits over the `math_FunctionSetWithDerivatives`
//! base; `occ::handle<Adaptor3d_Curve> tguide` (aliasing `guide` until
//! Set(First, Last) trims it) maps to an owned trimmed [`Curve3`] copy plus
//! an accessor, and likewise `tevol` (aliasing `fevol` until Set(First,
//! Last) trims it); `math_Vector` / `math_Matrix` map to `[f64; 2]` /
//! `Vec<Vec<f64>>` (OCCT X(i) -> x[i - 1], D(i, j) -> d[i - 1][j - 1]);
//! `gp_Circ` maps to the kernel [`Circle3`].
//!
//! The OCCT `Adaptor3d_CurveOnSurface cons1 / cons2` members (hxx L233-234)
//! have no rcad adaptor equivalent; the rcad port stores the (rst, surf)
//! pairs the adaptors wrap and transcribes the consumed operations from
//! Adaptor3d_CurveOnSurface.cxx (Value/EvalD0, D1/EvalD1 generic branch,
//! FirstParameter/LastParameter L977-987, Resolution L1364-1370) in the
//! `cons1_*` / `cons2_*` helpers below.
//!
//! Pending kernel dependency (marked GAP, plan 0.6): math_SVD (the
//! second-chance solver of IsSolution / Section-d1 — the OCCT !IsDone()
//! route is preserved).

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::p_confusion;
use rcad_kernel::geom::{Circle3, Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::math_svd::MathSvd;
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::{MatD, VecD};

use crate::geomalgo::geomfill::geom_fill::{get_circle, get_circle_d1};
use crate::geomalgo::law::law_function::LawFunctionHandle;

use super::brep_blend::BlendDecrochStatus;
use super::brep_blend_func::{
    blend_func_get_minimal_weights, blend_func_get_shape, blend_func_next_shape,
    BlendFuncSectionShape, ConvertParameterisationType,
};
use super::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending;
use super::brep_blend_func_consrad::{
    elclib_circle_parameter, geomfill_get_tolerance, geomfill_knots, geomfill_mults,
};
use super::brep_blend_func_evolrad::fusionne_intervalles;
use super::brep_blend_function::BlendAppFunction;
use super::brep_blend_point::BlendPoint;
use super::brep_blend_rst_rst_function::BlendRstRstFunction;

/// OCCT Eps constant (BRepBlend_RstRstEvolRad.cxx L40) — defined but unused
/// in the original file; kept for parity.
#[allow(dead_code)]
const EPS: f64 = 1.0e-15;

/// OCCT static t3dto2d (BRepBlend_RstRstEvolRad.cxx L42-52).
fn t3dto2d(a: &mut f64, b: &mut f64, av: DVec3, bv: DVec3, cv: DVec3) {
    let ab = av.dot(bv);
    let ac = av.dot(cv);
    let bc = bv.dot(cv);
    let bb = bv.dot(bv);
    let cc = cv.dot(cv);
    let deno = bb * cc - bc * bc;
    *a = (ab * cc - ac * bc) / deno;
    *b = (ac * bb - ab * bc) / deno;
}

/// Architecture mapping: OCCT consumes GeomFill's
/// `Convert_ParameterisationType` through both BlendFunc and GeomFill; the
/// identity mapping between the two rcad spellings of the OCCT enum
/// (pattern of brep_blend_surf_rst_const_rad::tconv).
fn tconv(t_conv: ConvertParameterisationType) -> rcad_kernel::base::convert::ConvertParameterisation {
    use rcad_kernel::base::convert::ConvertParameterisation;
    match t_conv {
        ConvertParameterisationType::TgtThetaOver2 => ConvertParameterisation::TgtThetaOver2,
        ConvertParameterisationType::TgtThetaOver2_1 => ConvertParameterisation::TgtThetaOver2_1,
        ConvertParameterisationType::TgtThetaOver2_2 => ConvertParameterisation::TgtThetaOver2_2,
        ConvertParameterisationType::TgtThetaOver2_3 => ConvertParameterisation::TgtThetaOver2_3,
        ConvertParameterisationType::TgtThetaOver2_4 => ConvertParameterisation::TgtThetaOver2_4,
        ConvertParameterisationType::QuasiAngular => ConvertParameterisation::QuasiAngular,
        ConvertParameterisationType::RationalC1 => ConvertParameterisation::RationalC1,
        ConvertParameterisationType::Polynomial => ConvertParameterisation::Polynomial,
    }
}

/// OCCT BRepBlend_RstRstEvolRad.
pub struct BlendRstRstEvolRad<'a> {
    // OCCT BRepBlend_RstRstEvolRad.hxx fields (L229-267).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) rst1: &'a Curve2d,
    pub(crate) rst2: &'a Curve2d,
    // OCCT: Adaptor3d_CurveOnSurface cons1(Rst1, Surf1) / cons2(Rst2, Surf2)
    // — carried as the wrapped pairs (see the module header architecture
    // note).
    pub(crate) guide: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> tguide — a handle aliasing guide
    // until Set(First, Last) trims it (cxx L230).  The rcad port owns the
    // trimmed copy.
    pub(crate) tguide: Option<Curve3>,
    pub(crate) ptrst1: DVec3,
    pub(crate) ptrst2: DVec3,
    pub(crate) pt2drst1: DVec2,
    pub(crate) pt2drst2: DVec2,
    pub(crate) prmrst1: f64,
    pub(crate) prmrst2: f64,
    pub(crate) istangent: bool,
    pub(crate) tgrst1: DVec3,
    pub(crate) tg2drst1: DVec2,
    pub(crate) tgrst2: DVec3,
    pub(crate) tg2drst2: DVec2,
    pub(crate) ray: f64,
    pub(crate) dray: f64,
    pub(crate) choix: i32,
    pub(crate) ptgui: DVec3,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    // OCCT: occ::handle<Adaptor3d_Surface> surfref1 / rstref1 / surfref2 /
    // rstref2 — null until Set(SurfRef1, RstRef1, SurfRef2, RstRef2)
    // (cxx L198-207).
    pub(crate) surfref1: Option<&'a Surface3>,
    pub(crate) rstref1: Option<&'a Curve2d>,
    pub(crate) surfref2: Option<&'a Surface3>,
    pub(crate) rstref2: Option<&'a Curve2d>,
    pub(crate) maxang: f64,
    pub(crate) minang: f64,
    pub(crate) distmin: f64,
    pub(crate) my_s_shape: BlendFuncSectionShape,
    pub(crate) my_t_conv: ConvertParameterisationType,
    // OCCT: occ::handle<Law_Function> tevol / fevol.
    pub(crate) tevol: Option<LawFunctionHandle>,
    pub(crate) fevol: LawFunctionHandle,
}

impl<'a> BlendRstRstEvolRad<'a> {
    /// Architecture mapping: OCCT `tguide` is a handle aliasing `guide`
    /// until Set(First, Last) replaces it with a trimmed copy (cxx L230).
    #[inline]
    pub(crate) fn tguide(&self) -> &Curve3 {
        self.tguide.as_ref().unwrap_or(self.guide)
    }

    /// Architecture mapping: OCCT `tevol` is a handle aliasing `fevol`
    /// until Set(First, Last) replaces it with a trimmed law (cxx L231).
    #[inline]
    pub(crate) fn tevol(&self) -> &LawFunctionHandle {
        self.tevol.as_ref().unwrap_or(&self.fevol)
    }

    /// OCCT Adaptor3d_CurveOnSurface cons1::Value (the EvalD0 generic branch
    /// of Adaptor3d_CurveOnSurface.cxx — the pcurve point carried on the
    /// wrapped surface).
    fn cons1_value(&self, w: f64) -> DVec3 {
        let puv = self.rst1.point_at(w);
        self.surf1.point_at(puv.x, puv.y)
    }

    /// OCCT Adaptor3d_CurveOnSurface cons2::Value (the EvalD0 generic branch
    /// of Adaptor3d_CurveOnSurface.cxx).
    fn cons2_value(&self, w: f64) -> DVec3 {
        let puv = self.rst2.point_at(w);
        self.surf2.point_at(puv.x, puv.y)
    }

    /// OCCT Adaptor3d_CurveOnSurface cons1::D1 (Adaptor3d_CurveOnSurface.cxx
    /// EvalD1 L1200-1241, the generic branch: the pcurve D1 composed with
    /// the carrier surface D1 — `D1.SetLinearForm(Duv.X(), D1U, Duv.Y(),
    /// D1V)`).
    fn cons1_d1(&self, w: f64) -> (DVec3, DVec3) {
        let puv = self.rst1.point_at(w);
        let duv = self.rst1.derivative_at(w);
        let (p, d1u, d1v) = self.surf1.derivatives(puv.x, puv.y);
        (p, duv.x * d1u + duv.y * d1v)
    }

    /// OCCT Adaptor3d_CurveOnSurface cons2::D1 (Adaptor3d_CurveOnSurface.cxx
    /// EvalD1 L1200-1241, the generic branch).
    fn cons2_d1(&self, w: f64) -> (DVec3, DVec3) {
        let puv = self.rst2.point_at(w);
        let duv = self.rst2.derivative_at(w);
        let (p, d1u, d1v) = self.surf2.derivatives(puv.x, puv.y);
        (p, duv.x * d1u + duv.y * d1v)
    }

    /// OCCT Adaptor3d_CurveOnSurface::Resolution (Adaptor3d_CurveOnSurface.cxx
    /// L1364-1370): `ru = mySurface->UResolution(R3d); rv = ...; return
    /// myCurve->Resolution(std::min(ru, rv));`  GAP (plan 0.6): the final
    /// Adaptor2d_Curve2d::Resolution step has no rcad equivalent; the
    /// established pending marker preserves the OCCT failure path.
    fn cons_resolution(&self, surf: &Surface3, r3d: f64) -> f64 {
        let ru = surf.u_resolution(r3d);
        let rv = surf.v_resolution(r3d);
        let _ = ru.min(rv);
        adaptor2d_curve2d_resolution_pending()
    }

    /// OCCT BRepBlend_RstRstEvolRad(Surf1, Rst1, Surf2, Rst2, CGuide, Evol)
    /// (BRepBlend_RstRstEvolRad.cxx L117-139).
    pub fn new(
        surf1: &'a Surface3,
        rst1: &'a Curve2d,
        surf2: &'a Surface3,
        rst2: &'a Curve2d,
        cguide: &'a Curve3,
        evol: LawFunctionHandle,
    ) -> Self {
        BlendRstRstEvolRad {
            surf1,
            surf2,
            rst1,
            rst2,
            guide: cguide,
            tguide: None, // OCCT: tguide(CGuide) — alias, see tguide().
            ptrst1: DVec3::ZERO,
            ptrst2: DVec3::ZERO,
            pt2drst1: DVec2::ZERO,
            pt2drst2: DVec2::ZERO,
            prmrst1: 0.0,
            prmrst2: 0.0,
            istangent: true,
            tgrst1: DVec3::ZERO,
            tg2drst1: DVec2::ZERO,
            tgrst2: DVec3::ZERO,
            tg2drst2: DVec2::ZERO,
            ray: 0.0,
            dray: 0.0,
            choix: 0,
            ptgui: DVec3::ZERO,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            normtg: 0.0,
            the_d: 0.0,
            surfref1: None,
            rstref1: None,
            surfref2: None,
            rstref2: None,
            maxang: -f64::MAX, // OCCT: RealFirst()
            minang: f64::MAX,  // OCCT: RealLast()
            distmin: f64::MAX, // OCCT: RealLast()
            my_s_shape: BlendFuncSectionShape::Rational,
            my_t_conv: ConvertParameterisationType::TgtThetaOver2,
            tevol: None, // OCCT: tevol(Evol) — alias, see tevol().
            fevol: evol, // OCCT L137-138: tevol = Evol; fevol = Evol.
        }
    }

    /// OCCT NbVariables() (BRepBlend_RstRstEvolRad.cxx L143-146) —
    /// returns 2.
    pub fn nb_variables(&self) -> usize {
        2
    }

    /// OCCT NbEquations() (BRepBlend_RstRstEvolRad.cxx L150-153) —
    /// returns 2.
    pub fn nb_equations(&self) -> usize {
        2
    }

    /// OCCT Value(X, F) (BRepBlend_RstRstEvolRad.cxx L157-166).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L159-160.
        self.ptrst1 = self.cons1_value(x[0]);
        self.ptrst2 = self.cons2_value(x[1]);

        // OCCT L162-163.
        f[0] = self.nplan.dot(self.ptrst1) + self.the_d;
        f[1] = self.nplan.dot(self.ptrst2) + self.the_d;

        true
    }

    /// OCCT Derivatives(X, D) (BRepBlend_RstRstEvolRad.cxx L170-184).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L174-175.
        let (ptrst1_v, d11) = self.cons1_d1(x[0]);
        self.ptrst1 = ptrst1_v;
        let (ptrst2_v, d21) = self.cons2_d1(x[1]);
        self.ptrst2 = ptrst2_v;

        // OCCT L177-181.
        d[0][0] = self.nplan.dot(d11);
        d[0][1] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = self.nplan.dot(d21);

        true
    }

    /// OCCT Values(X, F, D) (BRepBlend_RstRstEvolRad.cxx L188-194).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        self.value(x, f);
        self.derivatives(x, d);

        true
    }

    /// OCCT Set(SurfRef1, RstRef1, SurfRef2, RstRef2)
    /// (BRepBlend_RstRstEvolRad.cxx L198-207).
    pub fn set_ref(
        &mut self,
        surf_ref1: &'a Surface3,
        rst_ref1: &'a Curve2d,
        surf_ref2: &'a Surface3,
        rst_ref2: &'a Curve2d,
    ) {
        self.surfref1 = Some(surf_ref1);
        self.surfref2 = Some(surf_ref2);
        self.rstref1 = Some(rst_ref1);
        self.rstref2 = Some(rst_ref2);
    }

    /// OCCT Set(Param) (BRepBlend_RstRstEvolRad.cxx L211-224).
    pub fn set_param(&mut self, param: f64) {
        // OCCT L213-214.
        self.d1gui = DVec3::ZERO;
        self.nplan = DVec3::ZERO;
        // OCCT L215: tguide->D2(Param, ptgui, d1gui, d2gui).
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        self.d2gui = self.tguide().derivative2_at(param);
        // OCCT L216-221: theD is carried through the two-step gp_XYZ form
        // (the commented-out theD = -(nplan.Dot(ptgui)) is kept in the
        // original as a comment).
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        let the_d = self.nplan.dot(self.ptgui);
        self.the_d = the_d * (-1.0);
        // OCCT L223: tevol->D1(Param, ray, dray).
        let (mut ray, mut dray) = (0.0f64, 0.0f64);
        self.tevol().borrow_mut().d1(param, &mut ray, &mut dray);
        self.ray = ray;
        self.dray = dray;
    }

    /// OCCT Set(First, Last) (BRepBlend_RstRstEvolRad.cxx L228-232).
    pub fn set_interval(&mut self, first: f64, last: f64) {
        self.tguide = Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(self.guide.clone()),
            first,
            last,
        }));
        // OCCT L231: tevol = fevol->Trim(First, Last, 1.e-12).
        self.tevol = Some(self.fevol.borrow().trim(first, last, 1.0e-12));
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BRepBlend_RstRstEvolRad.cxx
    /// L236-240).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L238: Tolerance(1) = cons1.Resolution(Tol).
        tolerance[0] = self.cons_resolution(self.surf1, tol);
        // OCCT L239: Tolerance(2) = cons2.Resolution(Tol).
        tolerance[1] = self.cons_resolution(self.surf2, tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BRepBlend_RstRstEvolRad.cxx
    /// L244-250) — the cons FirstParameter / LastParameter are the pcurves'
    /// ranges (Adaptor3d_CurveOnSurface L977-987).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.rst1.default_domain()[0]; // cons1.FirstParameter
        inf_bound[1] = self.rst2.default_domain()[0]; // cons2.FirstParameter
        sup_bound[0] = self.rst1.default_domain()[1]; // cons1.LastParameter
        sup_bound[1] = self.rst2.default_domain()[1]; // cons2.LastParameter
    }

    /// OCCT IsSolution(Sol, Tol) (BRepBlend_RstRstEvolRad.cxx L254-368).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 2];
        let mut secmember = [0.0f64; 2];
        let mut gradsol = vec![vec![0.0f64; 2]; 2];

        // OCCT L265: Values(Sol, valsol, gradsol).
        self.values(sol, &mut valsol, &mut gradsol);

        // OCCT L267.
        if valsol[0].abs() <= tol && valsol[1].abs() <= tol {
            // Calculation of tangents
            // OCCT L271-274.
            self.prmrst1 = sol[0];
            self.pt2drst1 = self.rst1.point_at(self.prmrst1);
            self.prmrst2 = sol[1];
            self.pt2drst2 = self.rst2.point_at(self.prmrst2);

            // OCCT L276-277.
            let (ptrst1_v, d11) = self.cons1_d1(sol[0]);
            self.ptrst1 = ptrst1_v;
            let (ptrst2_v, d21) = self.cons2_d1(sol[1]);
            self.ptrst2 = ptrst2_v;

            // OCCT L279: dnplan.SetLinearForm(1. / normtg, d2gui,
            // -1. / normtg * (nplan.Dot(d2gui)), nplan).
            let dnplan = (1.0 / self.normtg) * self.d2gui
                + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            // OCCT L281-282.
            let temp = self.ptrst1 - self.ptgui;
            secmember[0] = self.normtg - dnplan.dot(temp);

            // OCCT L284-285.
            let temp = self.ptrst2 - self.ptgui;
            secmember[1] = self.normtg - dnplan.dot(temp);

            // OCCT L287-308: math_Gauss Resol(gradsol) and the SVD fallback.
            let mut a = MatD::new(2, 2);
            for r in 1..=2 {
                for c in 1..=2 {
                    a.set(r, c, gradsol[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                // OCCT L291: Resol.Solve(secmember).
                let mut x = VecD::new(2);
                for i in 1..=2 {
                    x.set(i, secmember[i - 1]);
                }
                resol.solve(&mut x);
                for i in 1..=2 {
                    secmember[i - 1] = x.get(i);
                }
                self.istangent = false;
            } else {
                // OCCT L296-307: math_SVD SingRS(gradsol);
                // if (SingRS.IsDone()) { math_Vector DEDT = secmember;
                // SingRS.Solve(DEDT, secmember, 1.e-6); istangent = false; }
                // else { istangent = true; }
                let mut sing_rs = MathSvd::new(&a);
                if sing_rs.is_done() {
                    let mut dedt = VecD::new(2);
                    let mut sol = VecD::new(2);
                    for i in 1..=2 {
                        dedt.set(i, secmember[i - 1]);
                        sol.set(i, secmember[i - 1]);
                    }
                    sing_rs.solve(&dedt, &mut sol, 1.0e-6);
                    for i in 1..=2 {
                        secmember[i - 1] = sol.get(i);
                    }
                    self.istangent = false;
                } else {
                    self.istangent = true;
                }
            }

            // OCCT L310-322.
            if !self.istangent {
                // OCCT L312-313.
                self.tgrst1 = secmember[0] * d11;
                self.tgrst2 = secmember[1] * d21;

                // OCCT L315-321 (gp_Pnt bid unused in the rcad form).
                let (_, d1urst1, d1vrst1) = self.surf1.derivatives(self.pt2drst1.x, self.pt2drst1.y);
                let mut a2d = 0.0f64;
                let mut b2d = 0.0f64;
                t3dto2d(&mut a2d, &mut b2d, self.tgrst1, d1urst1, d1vrst1);
                self.tg2drst1 = DVec2::new(a2d, b2d);
                let (_, d1urst2, d1vrst2) = self.surf2.derivatives(self.pt2drst2.x, self.pt2drst2.y);
                // OCCT L320 passes tgrst1 (not tgrst2) — translated literally.
                t3dto2d(&mut a2d, &mut b2d, self.tgrst1, d1urst2, d1vrst2);
                self.tg2drst2 = DVec2::new(a2d, b2d);
            }

            // OCCT L324-328.
            let mut center = DVec3::ZERO;
            let mut not_used = DVec3::ZERO;
            let is_center = self.center_circle_rst1_rst2(
                self.ptrst1,
                self.ptrst2,
                self.nplan,
                &mut center,
                &mut not_used,
            );

            // OCCT L330-333.
            if !is_center {
                return false;
            }

            // OCCT L335-338.
            let mut n1 = self.ptrst1 - center;
            let mut n2 = self.ptrst2 - center;
            n1 = n1.normalize_or_zero();
            n2 = n2.normalize_or_zero();

            // OCCT L340-341.
            let cosa = n1.dot(n2);
            let mut sina = self.nplan.dot(n1.cross(n2));

            // OCCT L343-346.
            if self.choix % 2 != 0 {
                sina = -sina; // nplan is changed into -nplan
            }

            // OCCT L348-352.
            let mut angle = cosa.acos();
            if sina < 0.0 {
                angle = 2.0 * std::f64::consts::PI - angle;
            }

            // OCCT L354-361: update of maxang / minang.
            if angle > self.maxang {
                self.maxang = angle;
            }
            if angle < self.minang {
                self.minang = angle;
            }
            // OCCT L362.
            self.distmin = self.distmin.min(self.ptrst1.distance(self.ptrst2));

            return true;
        }
        // OCCT L366-367.
        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BRepBlend_RstRstEvolRad.cxx L372-375).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT PointOnRst1() (BRepBlend_RstRstEvolRad.cxx L379-382).
    pub fn point_on_rst1(&self) -> DVec3 {
        self.ptrst1
    }

    /// OCCT PointOnRst2() (BRepBlend_RstRstEvolRad.cxx L386-389).
    pub fn point_on_rst2(&self) -> DVec3 {
        self.ptrst2
    }

    /// OCCT Pnt2dOnRst1() (BRepBlend_RstRstEvolRad.cxx L393-396).
    pub fn pnt2d_on_rst1(&self) -> DVec2 {
        self.pt2drst1
    }

    /// OCCT Pnt2dOnRst2() (BRepBlend_RstRstEvolRad.cxx L400-403).
    pub fn pnt2d_on_rst2(&self) -> DVec2 {
        self.pt2drst2
    }

    /// OCCT ParameterOnRst1() (BRepBlend_RstRstEvolRad.cxx L407-410).
    pub fn parameter_on_rst1(&self) -> f64 {
        self.prmrst1
    }

    /// OCCT ParameterOnRst2() (BRepBlend_RstRstEvolRad.cxx L414-417).
    pub fn parameter_on_rst2(&self) -> f64 {
        self.prmrst2
    }

    /// OCCT IsTangencyPoint() (BRepBlend_RstRstEvolRad.cxx L421-424).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnRst1() (BRepBlend_RstRstEvolRad.cxx L428-435).
    pub fn tangent_on_rst1(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstEvolRad::TangentOnRst1");
        }
        self.tgrst1
    }

    /// OCCT Tangent2dOnRst1() (BRepBlend_RstRstEvolRad.cxx L439-446).
    pub fn tangent_2d_on_rst1(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstEvolRad::Tangent2dOnRst1");
        }
        self.tg2drst1
    }

    /// OCCT TangentOnRst2() (BRepBlend_RstRstEvolRad.cxx L450-457).
    pub fn tangent_on_rst2(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstEvolRad::TangentOnRst2");
        }
        self.tgrst2
    }

    /// OCCT Tangent2dOnRst2() (BRepBlend_RstRstEvolRad.cxx L461-468).
    pub fn tangent_2d_on_rst2(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstEvolRad::Tangent2dOnRst2");
        }
        self.tg2drst2
    }

    /// OCCT Decroch(Sol, NRst1, TgRst1, NRst2, TgRst2)
    /// (BRepBlend_RstRstEvolRad.cxx L472-550).
    ///
    /// Parameter-name note: the OCCT header names the out-parameters
    /// (NRst1, TgRst1, NRst2, TgRst2) while the consumer call site
    /// (BRepBlend_RstRstLineBuilder.cxx L1961) passes (tgrst1, norst1,
    /// tgrst2, norst2); the rcad trait parameter names follow the call
    /// site, so the OCCT `NRst1` corresponds to the rcad `tgrst1` slot and
    /// the OCCT `TgRst1` to the rcad `nrrst1` slot.
    pub fn decroch(
        &self,
        sol: &[f64],
        tgrst1: &mut DVec3,
        nrrst1: &mut DVec3,
        tgrst2: &mut DVec3,
        nrrst2: &mut DVec3,
    ) -> BlendDecrochStatus {
        // OCCT L484-485: rstref1->Value(Sol(1)).Coord(u, v);
        // surfref1->D1(u, v, PtTmp1, d1u, d1v).
        let (u, v) = {
            let p = self.rstref1.expect("rstref1").point_at(sol[0]);
            (p.x, p.y)
        };
        let (pt_tmp1, d1u, d1v) = self.surfref1.expect("surfref1").derivatives(u, v);
        // Normal to the reference surface 1
        // OCCT L487: NRst1 = d1u.Crossed(d1v).
        let n_rst1 = d1u.cross(d1v);
        *tgrst1 = n_rst1;
        // OCCT L488-489.
        let (u, v) = {
            let p = self.rstref2.expect("rstref2").point_at(sol[1]);
            (p.x, p.y)
        };
        let (pt_tmp2, d1u, d1v) = self.surfref2.expect("surfref2").derivatives(u, v);
        // Normal to the reference surface 2
        // OCCT L491: NRst2 = d1u.Crossed(d1v).
        let n_rst2 = d1u.cross(d1v);
        *tgrst2 = n_rst2;

        // OCCT L493.
        let mut center = DVec3::ZERO;
        let mut not_used = DVec3::ZERO;
        self.center_circle_rst1_rst2(pt_tmp1, pt_tmp2, self.nplan, &mut center, &mut not_used);

        // OCCT L495-498.
        let norm = self.nplan.cross(n_rst1).length();
        let unsurnorm = 1.0 / norm;
        let mut n_rst1_in_plane =
            (self.nplan.dot(n_rst1) * unsurnorm) * self.nplan + (-unsurnorm) * n_rst1;

        // OCCT L500: centptrst.SetXYZ(PtTmp1.XYZ() - Center.XYZ()).
        let mut centptrst = pt_tmp1 - center;

        // OCCT L502-505.
        if centptrst.dot(n_rst1_in_plane) < 0.0 {
            n_rst1_in_plane = -n_rst1_in_plane;
        }

        // OCCT L507: TgRst1 = nplan.Crossed(centptrst).
        *nrrst1 = self.nplan.cross(centptrst);

        // OCCT L509-519.
        let norm = self.nplan.cross(n_rst2).length();
        let unsurnorm = 1.0 / norm;
        let mut n_rst2_in_plane =
            (self.nplan.dot(n_rst2) * unsurnorm) * self.nplan + (-unsurnorm) * n_rst2;
        centptrst = pt_tmp2 - center;

        // OCCT L514-517.
        if centptrst.dot(n_rst2_in_plane) < 0.0 {
            n_rst2_in_plane = -n_rst2_in_plane;
        }

        // OCCT L519: TgRst2 = nplan.Crossed(centptrst).
        *nrrst2 = self.nplan.cross(centptrst);

        // OCCT L521-525.
        if self.choix % 2 != 0 {
            *nrrst1 = -*nrrst1;
            *nrrst2 = -*nrrst2;
        }

        // Vectors are returned
        // OCCT L528-549.
        if n_rst1_in_plane.dot(*nrrst1) > -1.0e-10 {
            if n_rst2_in_plane.dot(*nrrst2) < 1.0e-10 {
                BlendDecrochStatus::DecrochBoth
            } else {
                BlendDecrochStatus::DecrochRst1
            }
        } else {
            if n_rst2_in_plane.dot(*nrrst2) < 1.0e-10 {
                BlendDecrochStatus::DecrochRst2
            } else {
                BlendDecrochStatus::NoDecroch
            }
        }
    }

    /// OCCT Set(Choix) (BRepBlend_RstRstEvolRad.cxx L554-557).
    pub fn set(&mut self, choix: i32) {
        self.choix = choix;
    }

    /// OCCT Set(BlendFunc_SectionShape) (BRepBlend_RstRstEvolRad.cxx
    /// L561-564).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_s_shape = type_section;
    }

    /// OCCT CenterCircleRst1Rst2(PtRst1, PtRst2, np, Center, VdMed)
    /// (BRepBlend_RstRstEvolRad.cxx L570-609) — calculates the center of
    /// the circle passing by two points of restrictions.
    pub fn center_circle_rst1_rst2(
        &self,
        pt_rst1: DVec3,
        pt_rst2: DVec3,
        np: DVec3,
        center: &mut DVec3,
        vd_med: &mut DVec3,
    ) -> bool {
        // OCCT L577: gp_Vec rst1rst2(PtRst1, PtRst2).
        let rst1rst2 = pt_rst2 - pt_rst1;

        // Calculate the center of the circle
        // OCCT L583-585.
        *vd_med = rst1rst2.cross(np);
        let norm2 = rst1rst2.length_squared();
        let mut dist = self.ray * self.ray - 0.25 * norm2;

        // OCCT L587-590.
        if self.choix > 2 {
            *vd_med = -*vd_med;
        }

        // OCCT L592-595.
        if dist < -1.0e-07 {
            return false;
        }

        // OCCT L597-606.
        if dist > 1.0e-07 {
            dist = dist.sqrt();
            let vdmed_nor = vd_med.normalize_or_zero();
            *center = 0.5 * rst1rst2 + pt_rst1 + dist * vdmed_nor;
        } else {
            *center = 0.5 * rst1rst2 + pt_rst1;
        }

        true
    }

    /// OCCT Section(Param, U, V, Pdeb, Pfin, C)
    /// (BRepBlend_RstRstEvolRad.cxx L613-654).
    #[allow(clippy::too_many_arguments)]
    pub fn section(
        &mut self,
        param: f64,
        u: f64,
        v: f64,
        pdeb: &mut f64,
        pfin: &mut f64,
        c: &mut Circle3,
    ) {
        // OCCT L623-625: tguide->D1(Param, ptgui, d1gui);
        // ray = tevol->Value(Param); np = d1gui.Normalized().
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        // OCCT L624: ray = tevol->Value(Param) — the law read and the
        // self.ray write are split across statements (Rc<RefCell> borrow).
        let law_ray = self.tevol().borrow_mut().value(param);
        self.ray = law_ray;
        let mut np = self.d1gui.normalize_or_zero();        // OCCT L626-627.
        self.ptrst1 = self.cons1_value(u);
        self.ptrst2 = self.cons2_value(v);

        // OCCT L629 (the IsCenter flag is discarded in the original).
        let mut center = DVec3::ZERO;
        let mut not_used = DVec3::ZERO;
        self.center_circle_rst1_rst2(self.ptrst1, self.ptrst2, np, &mut center, &mut not_used);

        // OCCT L631-632.
        c.radius = self.ray.abs();
        let ns = (self.ptrst1 - center).normalize_or_zero();

        // OCCT L634-637.
        if self.choix % 2 != 0 {
            np = -np;
        }

        // OCCT L639: C.SetPosition(gp_Ax2(Center, np, ns))
        // (gp_Ax2(P, N, Vx): the Y direction is N ^ Vx).
        c.center = center;
        c.normal = np;
        c.x_dir = ns;
        c.y_dir = np.cross(ns);
        // OCCT L640-641: Pdeb = 0; Pfin = ElCLib::Parameter(C, ptrst2).
        *pdeb = 0.0;
        *pfin = elclib_circle_parameter(c, self.ptrst2);

        // Test negative and quasi null angles: Special case
        // OCCT L644-649.
        if *pfin > 1.5 * std::f64::consts::PI {
            np = -np;
            c.normal = np;
            c.y_dir = np.cross(ns);
            *pfin = elclib_circle_parameter(c, self.ptrst2);
        }
        // OCCT L650-653.
        if *pfin < p_confusion() {
            *pfin += p_confusion();
        }
    }

    /// OCCT IsRational() (BRepBlend_RstRstEvolRad.cxx L658-661).
    pub fn is_rational(&self) -> bool {
        self.my_s_shape == BlendFuncSectionShape::Rational
            || self.my_s_shape == BlendFuncSectionShape::QuasiAngular
    }

    /// OCCT GetSectionSize() (BRepBlend_RstRstEvolRad.cxx L665-668).
    pub fn get_section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// OCCT GetMinimalWeight(Weights) (BRepBlend_RstRstEvolRad.cxx
    /// L672-676).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        blend_func_get_minimal_weights(self.my_s_shape, self.my_t_conv, self.minang, self.maxang, weigths);
        // It is supposed that it does not depend on the Radius!
    }

    /// OCCT NbIntervals(S) (BRepBlend_RstRstEvolRad.cxx L680-699).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let nb_int_courbe = self.guide.nb_intervals(blend_func_next_shape(s));
        let nb_int_loi = self.fevol.borrow().nb_intervals(s);

        if nb_int_loi == 1 {
            return nb_int_courbe;
        }

        let mut int_c = vec![0.0f64; nb_int_courbe + 1];
        let mut int_l = vec![0.0f64; nb_int_loi + 1];
        let mut inter: Vec<f64> = Vec::new();
        self.guide.intervals(&mut int_c, blend_func_next_shape(s));
        self.fevol.borrow().intervals(&mut int_l, s);

        fusionne_intervalles(&int_c, &int_l, &mut inter);
        inter.len() - 1
    }

    /// OCCT Intervals(T, S) (BRepBlend_RstRstEvolRad.cxx L703-727).
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let nb_int_courbe = self.guide.nb_intervals(blend_func_next_shape(s));
        let nb_int_loi = self.fevol.borrow().nb_intervals(s);

        if nb_int_loi == 1 {
            let mut intervals = Vec::new();
            self.guide.intervals(&mut intervals, blend_func_next_shape(s));
            for (dst, src) in t.iter_mut().zip(intervals) {
                *dst = src;
            }
        } else {
            let mut int_c = vec![0.0f64; nb_int_courbe + 1];
            let mut int_l = vec![0.0f64; nb_int_loi + 1];
            let mut inter: Vec<f64> = Vec::new();
            self.guide.intervals(&mut int_c, blend_func_next_shape(s));
            self.fevol.borrow().intervals(&mut int_l, s);

            fusionne_intervalles(&int_c, &int_l, &mut inter);
            for (dst, src) in t.iter_mut().zip(inter) {
                *dst = src;
            }
        }
    }

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
    /// (BRepBlend_RstRstEvolRad.cxx L731-735).
    pub fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        *nb_poles_2d = 2;
        blend_func_get_shape(
            self.my_s_shape,
            self.maxang,
            nb_poles,
            nb_knots,
            degree,
            &mut self.my_t_conv,
        );
    }

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)
    /// (BRepBlend_RstRstEvolRad.cxx L742-755) — tolerances used for
    /// approximations.
    pub fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        let low = 0usize; // OCCT: Tol3d.Lower()
        let up = tol3d.len() - 1; // OCCT: Tol3d.Upper()
        let tol = geomfill_get_tolerance(self.my_t_conv, self.minang, self.ray.abs(), angle_tol, surf_tol);
        for v in tol1d.iter_mut() {
            *v = surf_tol;
        }
        for v in tol3d.iter_mut() {
            *v = surf_tol;
        }
        tol3d[low + 1] = tol.min(surf_tol);
        tol3d[up - 1] = tol.min(surf_tol);
        tol3d[low] = tol.min(bound_tol);
        tol3d[up] = tol.min(bound_tol);
    }

    /// OCCT Knots(TKnots) (BRepBlend_RstRstEvolRad.cxx L759-762).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        geomfill_knots(self.my_t_conv, tknots);
    }

    /// OCCT Mults(TMults) (BRepBlend_RstRstEvolRad.cxx L766-769).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        geomfill_mults(self.my_t_conv, tmults);
    }

    /// OCCT Section(P, Poles, Poles2d, Weights)
    /// (BRepBlend_RstRstEvolRad.cxx L773-835).
    pub fn section_simple(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        let prm = p.parameter();
        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()

        // OCCT L786-788: tguide->D1(prm, ptgui, d1gui);
        // ray = tevol->Value(prm); nplan = d1gui.Normalized().
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        // OCCT L787: ray = tevol->Value(prm) — the law read and the
        // self.ray write are split across statements (Rc<RefCell> borrow).
        let law_ray = self.tevol().borrow_mut().value(prm);
        self.ray = law_ray;
        self.nplan = self.d1gui.normalize_or_zero();

        // OCCT L790-791.
        let u = p.parameter_on_c1();
        let v = p.parameter_on_c2();

        // OCCT L793-794.
        let pt2d1 = self.rst1.point_at(u);
        let pt2d2 = self.rst2.point_at(v);

        // OCCT L796-798.
        self.ptrst1 = self.cons1_value(u);
        self.ptrst2 = self.cons2_value(v);
        self.distmin = self.distmin.min(self.ptrst1.distance(self.ptrst2));

        // OCCT L800-801.
        poles_2d[0] = DVec2::new(pt2d1.x, pt2d1.y);
        poles_2d[poles_2d.len() - 1] = DVec2::new(pt2d2.x, pt2d2.y);

        // Linear Case
        // OCCT L804-811.
        if self.my_s_shape == BlendFuncSectionShape::Linear {
            poles[low] = self.ptrst1;
            poles[upp] = self.ptrst2;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            return;
        }

        // Calculate the center of the circle
        // OCCT L814 (the IsCenter flag is discarded in the original).
        let mut center = DVec3::ZERO;
        let mut not_used = DVec3::ZERO;
        self.center_circle_rst1_rst2(self.ptrst1, self.ptrst2, self.nplan, &mut center, &mut not_used);

        // normals to the section with points
        // OCCT L817-818.
        let n1 = (self.ptrst1 - center).normalize_or_zero();
        let n2 = (self.ptrst2 - center).normalize_or_zero();

        // OCCT L820-823.
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
        }

        // OCCT L825-834: GeomFill::GetCircle(myTConv, n1, n2, nplan, ptrst1,
        // ptrst2, std::abs(ray), Center, Poles, Weights).
        get_circle(
            tconv(self.my_t_conv),
            n1,
            n2,
            self.nplan,
            self.ptrst1,
            self.ptrst2,
            self.ray.abs(),
            center,
            poles,
            weigths,
        );
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)
    /// (BRepBlend_RstRstEvolRad.cxx L839-1059) — used for the first and
    /// last section.
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
        let prm = p.parameter();
        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()

        // OCCT L864-868.
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        self.d2gui = self.tguide().derivative2_at(prm);
        let (mut ray, mut dray) = (0.0f64, 0.0f64);
        self.tevol().borrow_mut().d1(prm, &mut ray, &mut dray);
        self.ray = ray;
        self.dray = dray;
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        // OCCT L868: dnplan.SetLinearForm(1. / normtg, d2gui,
        // -1. / normtg * (nplan.Dot(d2gui)), nplan).
        let mut dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        // OCCT L870-873.
        let mut sol = [0.0f64; 2];
        self.prmrst1 = p.parameter_on_c1();
        sol[0] = self.prmrst1;
        self.prmrst2 = p.parameter_on_c2();
        sol[1] = self.prmrst2;
        self.pt2drst1 = self.rst1.point_at(self.prmrst1);
        self.pt2drst2 = self.rst2.point_at(self.prmrst2);

        // OCCT L875: Values(sol, valsol, gradsol).
        let mut valsol = [0.0f64; 2];
        let mut gradsol = vec![vec![0.0f64; 2]; 2];
        self.values(&sol, &mut valsol, &mut gradsol);

        // OCCT L877-878.
        let (ptrst1_v, d11) = self.cons1_d1(sol[0]);
        self.ptrst1 = ptrst1_v;
        let (ptrst2_v, d21) = self.cons2_d1(sol[1]);
        self.ptrst2 = ptrst2_v;

        // OCCT L880-884.
        let mut secmember = [0.0f64; 2];
        let temp = self.ptrst1 - self.ptgui;
        secmember[0] = self.normtg - dnplan.dot(temp);

        let temp = self.ptrst2 - self.ptgui;
        secmember[1] = self.normtg - dnplan.dot(temp);

        // OCCT L886-907: math_Gauss Resol(gradsol, 1.e-9) and the SVD
        // fallback (the rcad MathGauss carries no tolerance argument).
        let mut a = MatD::new(2, 2);
        for r in 1..=2 {
            for c in 1..=2 {
                a.set(r, c, gradsol[r - 1][c - 1]);
            }
        }
        let resol = MathGauss::new(&a);
        let istgt: bool;
        if resol.is_done() {
            istgt = false;
            // OCCT L891: Resol.Solve(secmember).
            let mut x = VecD::new(2);
            for i in 1..=2 {
                x.set(i, secmember[i - 1]);
            }
            resol.solve(&mut x);
            for i in 1..=2 {
                secmember[i - 1] = x.get(i);
            }
        } else {
            // OCCT L895-906: math_SVD SingRS(gradsol);
            // if (SingRS.IsDone()) { math_Vector DEDT(1, 2); DEDT = secmember;
            // SingRS.Solve(DEDT, secmember, 1.e-6); istgt = false; }
            // else { istgt = true; }
            let mut sing_rs = MathSvd::new(&a);
            if sing_rs.is_done() {
                let mut dedt = VecD::new(2);
                let mut sol = VecD::new(2);
                for i in 1..=2 {
                    dedt.set(i, secmember[i - 1]);
                    sol.set(i, secmember[i - 1]);
                }
                sing_rs.solve(&dedt, &mut sol, 1.0e-6);
                for i in 1..=2 {
                    secmember[i - 1] = sol.get(i);
                }
                istgt = false;
            } else {
                istgt = true;
            }
        }

        // OCCT L909-921.
        let mut med = DVec3::ZERO;
        let rst1rst2 = self.ptrst2 - self.ptrst1;
        let mut center = DVec3::ZERO;
        let is_center =
            self.center_circle_rst1_rst2(self.ptrst1, self.ptrst2, self.nplan, &mut center, &mut med);
        if !is_center {
            return false;
        }

        let normmed = med.length();
        med = med.normalize_or_zero();
        // OCCT L921: gp_Vec n1(Center, ptrst1), n2(Center, ptrst2) — these
        // are normalized only later, at L974-975, after the !istgt block
        // (the !istgt branch consumes the unnormalized vectors).
        let mut n1 = self.ptrst1 - center;
        let mut n2 = self.ptrst2 - center;

        // OCCT L923-972: secmember contains derivatives of parameters on
        // curves corresponding to t.
        let mut d1n1 = DVec3::ZERO; // only read on the !istgt path
        let mut d1n2 = DVec3::ZERO; // only read on the !istgt path
        if !istgt {
            // OCCT L927-928.
            self.tgrst1 = secmember[0] * d11;
            self.tgrst2 = secmember[1] * d21;

            // OCCT L930-935.
            let mut d1rst1rst2;
            let norm2 = rst1rst2.length_squared();
            d1rst1rst2 = self.tgrst2 - self.tgrst1;
            let mut dist = ray * ray - 0.25 * norm2;
            let invdray = dray / ray;

            // OCCT L937-962.
            if dist > 1.0e-07 {
                // OCCT L939-944.
                let d1p1p2crosnp = d1rst1rst2.cross(self.nplan) + rst1rst2.cross(dnplan);
                // derivative of the bisector
                let mut dmed = d1p1p2crosnp - med.dot(d1p1p2crosnp) * med;
                dmed /= normmed;
                dist = dist.sqrt();
                // OCCT L947.
                let d1dist = (ray * dray - 0.25 * rst1rst2.dot(d1rst1rst2)) / dist;

                // OCCT L949-952.
                if self.choix > 2 {
                    dmed = -dmed;
                }

                // derivative of the coefficient Dist is located in dmed
                // OCCT L955-961.
                dmed = dist * dmed + d1dist * med;
                d1rst1rst2 *= 0.5;
                // derivative of the Normal to the curve in P1
                d1n1 = -(d1rst1rst2 + dmed + invdray * n1) / ray;

                // derivative of the Normal to the curve in P2
                d1n2 = (d1rst1rst2 - dmed - invdray * n2) / ray;
            } else {
                // OCCT L964-971.
                d1rst1rst2 *= 0.5;
                // Normal to the curve in P1
                d1n1 = -(d1rst1rst2 + invdray * n1) / ray;

                // Normal to the curve in P2
                d1n2 = (d1rst1rst2 - invdray * n2) / ray;
            }
        }

        // OCCT L974-975.
        n1 = n1.normalize_or_zero();
        n2 = n2.normalize_or_zero();

        // Tops 2D
        // OCCT L979-991.
        poles_2d[0] = DVec2::new(self.pt2drst1.x, self.pt2drst1.y);
        poles_2d[poles_2d.len() - 1] = DVec2::new(self.pt2drst2.x, self.pt2drst2.y);
        if !istgt {
            let (_, d1urst, d1vrst) = self.surf1.derivatives(self.pt2drst1.x, self.pt2drst1.y);
            let mut a2d = 0.0f64;
            let mut b2d = 0.0f64;
            t3dto2d(&mut a2d, &mut b2d, self.tgrst1, d1urst, d1vrst);
            d_poles_2d[0] = DVec2::new(a2d, b2d);

            let (_, d1urst, d1vrst) = self.surf2.derivatives(self.pt2drst2.x, self.pt2drst2.y);
            t3dto2d(&mut a2d, &mut b2d, self.tgrst2, d1urst, d1vrst);
            d_poles_2d[d_poles_2d.len() - 1] = DVec2::new(a2d, b2d);
        }

        // Linear Case
        // OCCT L994-1008.
        if self.my_s_shape == BlendFuncSectionShape::Linear {
            poles[low] = self.ptrst1;
            poles[upp] = self.ptrst2;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            if !istgt {
                d_poles[low] = self.tgrst1;
                d_poles[upp] = self.tgrst2;
                d_weigths[low] = 0.0;
                d_weigths[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle
        // tangent to the center of the circle
        // OCCT L1012-1015: tgct.SetLinearForm(-ray, d1n1, -dray, n1, tgrst1).
        let tgct: DVec3;
        if !istgt {
            tgct = -ray * d1n1 + (-dray) * n1 + self.tgrst1;
        } else {
            tgct = DVec3::ZERO; // the OCCT tgct is only read on the !istgt path
        }

        // OCCT L1017-1021.
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
            dnplan = -dnplan;
        }

        // OCCT L1023-1058.
        if !istgt {
            get_circle_d1(
                tconv(self.my_t_conv),
                n1,
                n2,
                d1n1,
                d1n2,
                self.nplan,
                dnplan,
                self.ptrst1,
                self.ptrst2,
                self.tgrst1,
                self.tgrst2,
                ray.abs(),
                dray,
                center,
                tgct,
                poles,
                d_poles,
                weigths,
                d_weigths,
            )
        } else {
            // OCCT L1047-1056.
            get_circle(
                tconv(self.my_t_conv),
                n1,
                n2,
                self.nplan,
                self.ptrst1,
                self.ptrst2,
                ray.abs(),
                center,
                poles,
                weigths,
            );
            false
        }
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weights, DWeights, D2Weights) (BRepBlend_RstRstEvolRad.cxx
    /// L1063-1075) — returns false.
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

    /// OCCT Resolution(IC2d, Tol, TolU, TolV)
    /// (BRepBlend_RstRstEvolRad.cxx L1077-1092).
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

impl<'a> FunctionSetWithDerivatives for BlendRstRstEvolRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendRstRstEvolRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendRstRstEvolRad::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendRstRstEvolRad::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendRstRstEvolRad::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendRstRstEvolRad::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendRstRstEvolRad<'a> {
    fn set_param(&mut self, param: f64) {
        BlendRstRstEvolRad::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendRstRstEvolRad::set_interval(self, first, last)
    }

    fn pnt1(&self) -> DVec3 {
        BlendRstRstFunction::pnt1(self)
    }

    fn pnt2(&self) -> DVec3 {
        BlendRstRstFunction::pnt2(self)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendRstRstEvolRad::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendRstRstEvolRad::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendRstRstEvolRad::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendRstRstEvolRad::get_minimal_distance(self)
    }

    fn is_rational(&self) -> bool {
        BlendRstRstEvolRad::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendRstRstEvolRad::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendRstRstEvolRad::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendRstRstEvolRad::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendRstRstEvolRad::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendRstRstEvolRad::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendRstRstEvolRad::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendRstRstEvolRad::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendRstRstEvolRad::mults(self, tmults)
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
        BlendRstRstEvolRad::section_d1(
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
        BlendRstRstEvolRad::section_simple(self, p, poles, poles_2d, weigths)
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
        BlendRstRstEvolRad::section_d2(
            self, p, poles, d_poles, d2_poles, poles_2d, d_poles_2d, d2_poles_2d, weigths,
            d_weigths, d2_weigths,
        )
    }

    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        BlendRstRstEvolRad::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendRstRstFunction for BlendRstRstEvolRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendRstRstEvolRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendRstRstEvolRad::nb_equations(self)
    }

    fn point_on_rst1(&self) -> DVec3 {
        BlendRstRstEvolRad::point_on_rst1(self)
    }

    fn point_on_rst2(&self) -> DVec3 {
        BlendRstRstEvolRad::point_on_rst2(self)
    }

    fn pnt2d_on_rst1(&self) -> DVec2 {
        BlendRstRstEvolRad::pnt2d_on_rst1(self)
    }

    fn pnt2d_on_rst2(&self) -> DVec2 {
        BlendRstRstEvolRad::pnt2d_on_rst2(self)
    }

    fn parameter_on_rst1(&self) -> f64 {
        BlendRstRstEvolRad::parameter_on_rst1(self)
    }

    fn parameter_on_rst2(&self) -> f64 {
        BlendRstRstEvolRad::parameter_on_rst2(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendRstRstEvolRad::is_tangency_point(self)
    }

    fn tangent_on_rst1(&self) -> DVec3 {
        BlendRstRstEvolRad::tangent_on_rst1(self)
    }

    fn tangent_2d_on_rst1(&self) -> DVec2 {
        BlendRstRstEvolRad::tangent_2d_on_rst1(self)
    }

    fn tangent_on_rst2(&self) -> DVec3 {
        BlendRstRstEvolRad::tangent_on_rst2(self)
    }

    fn tangent_2d_on_rst2(&self) -> DVec2 {
        BlendRstRstEvolRad::tangent_2d_on_rst2(self)
    }

    fn decroch(
        &self,
        sol: &[f64],
        tgrst1: &mut DVec3,
        nrrst1: &mut DVec3,
        tgrst2: &mut DVec3,
        nrrst2: &mut DVec3,
    ) -> BlendDecrochStatus {
        BlendRstRstEvolRad::decroch(self, sol, tgrst1, nrrst1, tgrst2, nrrst2)
    }
}
