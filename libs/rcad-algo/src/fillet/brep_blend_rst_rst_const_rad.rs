//! OCCT BRepBlend_RstRstConstRad (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_RstRstConstRad.hxx (L40-263) + BRepBlend_RstRstConstRad.cxx
//! (whole file L38-981).  Copy of CSConstRad with pcurves on surfaces as
//! supports (hxx L39-40); the SimulSurf / PerformSurf rst/rst constant-radius
//! arms instantiate it as `func` (ChFi3d_FilBuilder.cxx L1301 / L2112).
//!
//! Architecture mappings (mirroring [`super::brep_blend_surf_rst_const_rad`]):
//! `class BRepBlend_RstRstConstRad : public Blend_RstRstFunction` is
//! expressed by implementing the [`BlendRstRstFunction`] and
//! [`BlendAppFunction`] traits over the `math_FunctionSetWithDerivatives`
//! base; `occ::handle<Adaptor3d_Curve> tguide` (aliasing `guide` until
//! Set(First, Last) trims it) maps to an owned trimmed [`Curve3`] copy plus
//! an accessor; `math_Vector` / `math_Matrix` map to `[f64; 2]` /
//! `Vec<Vec<f64>>` (OCCT X(i) -> x[i - 1], D(i, j) -> d[i - 1][j - 1]);
//! `gp_Circ` maps to the kernel [`Circle3`].
//!
//! The OCCT `Adaptor3d_CurveOnSurface cons1 / cons2` members (hxx L231-232)
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

use super::brep_blend::BlendDecrochStatus;
use super::brep_blend_func::{
    blend_func_get_minimal_weights, blend_func_get_shape, blend_func_next_shape,
    BlendFuncSectionShape, ConvertParameterisationType,
};
use super::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending;
use super::brep_blend_func_consrad::{
    elclib_circle_parameter, geomfill_get_tolerance, geomfill_knots, geomfill_mults,
};
use super::brep_blend_function::BlendAppFunction;
use super::brep_blend_point::BlendPoint;
use super::brep_blend_rst_rst_function::BlendRstRstFunction;

/// OCCT Eps constant (BRepBlend_RstRstConstRad.cxx L38) — defined but
/// unused in the original file; kept for parity.
#[allow(dead_code)]
const EPS: f64 = 1.0e-15;

/// OCCT static t3dto2d (BRepBlend_RstRstConstRad.cxx L40-50).
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

/// OCCT BRepBlend_RstRstConstRad — function to build a surfaces blend
/// between two pcurves (restrictions), constant radius rolling ball.
pub struct BlendRstRstConstRad<'a> {
    // OCCT BRepBlend_RstRstConstRad.hxx fields (L227-262).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) rst1: &'a Curve2d,
    pub(crate) rst2: &'a Curve2d,
    // OCCT: Adaptor3d_CurveOnSurface cons1(Rst1, Surf1) — carried as the
    // wrapped pair (see the module header architecture note).
    // OCCT: Adaptor3d_CurveOnSurface cons2(Rst2, Surf2) — likewise.
    pub(crate) guide: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> tguide — a handle aliasing guide
    // until Set(First, Last) trims it (cxx L66).  The rcad port owns the
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
    pub(crate) choix: i32,
    pub(crate) ptgui: DVec3,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    // OCCT: occ::handle<Adaptor3d_Surface> surfref1 / rstref1 / surfref2 /
    // rstref2 — null until Set(SurfRef1, RstRef1, SurfRef2, RstRef2)
    // (cxx L138-147).
    pub(crate) surfref1: Option<&'a Surface3>,
    pub(crate) rstref1: Option<&'a Curve2d>,
    pub(crate) surfref2: Option<&'a Surface3>,
    pub(crate) rstref2: Option<&'a Curve2d>,
    pub(crate) maxang: f64,
    pub(crate) minang: f64,
    pub(crate) distmin: f64,
    pub(crate) my_s_shape: BlendFuncSectionShape,
    pub(crate) my_t_conv: ConvertParameterisationType,
}

impl<'a> BlendRstRstConstRad<'a> {
    /// Architecture mapping: OCCT `tguide` is a handle aliasing `guide`
    /// until Set(First, Last) replaces it with a trimmed copy (cxx L163).
    #[inline]
    pub(crate) fn tguide(&self) -> &Curve3 {
        self.tguide.as_ref().unwrap_or(self.guide)
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

    /// OCCT BRepBlend_RstRstConstRad(Surf1, Rst1, Surf2, Rst2, CGuide)
    /// (BRepBlend_RstRstConstRad.cxx L54-79).
    pub fn new(
        surf1: &'a Surface3,
        rst1: &'a Curve2d,
        surf2: &'a Surface3,
        rst2: &'a Curve2d,
        cguide: &'a Curve3,
    ) -> Self {
        BlendRstRstConstRad {
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
        }
    }

    /// OCCT NbVariables() (BRepBlend_RstRstConstRad.cxx L83-86) — returns 2.
    pub fn nb_variables(&self) -> usize {
        2
    }

    /// OCCT NbEquations() (BRepBlend_RstRstConstRad.cxx L90-93) — returns 2.
    pub fn nb_equations(&self) -> usize {
        2
    }

    /// OCCT Value(X, F) (BRepBlend_RstRstConstRad.cxx L97-106).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L99-100.
        self.ptrst1 = self.cons1_value(x[0]);
        self.ptrst2 = self.cons2_value(x[1]);

        // OCCT L102-103.
        f[0] = self.nplan.dot(self.ptrst1) + self.the_d;
        f[1] = self.nplan.dot(self.ptrst2) + self.the_d;

        true
    }

    /// OCCT Derivatives(X, D) (BRepBlend_RstRstConstRad.cxx L110-124).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L114-115.
        let (ptrst1_v, d11) = self.cons1_d1(x[0]);
        self.ptrst1 = ptrst1_v;
        let (ptrst2_v, d21) = self.cons2_d1(x[1]);
        self.ptrst2 = ptrst2_v;

        // OCCT L117-121.
        d[0][0] = self.nplan.dot(d11);
        d[0][1] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = self.nplan.dot(d21);

        true
    }

    /// OCCT Values(X, F, D) (BRepBlend_RstRstConstRad.cxx L128-134).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        self.value(x, f);
        self.derivatives(x, d);

        true
    }

    /// OCCT Set(SurfRef1, RstRef1, SurfRef2, RstRef2)
    /// (BRepBlend_RstRstConstRad.cxx L138-147).
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

    /// OCCT Set(Param) (BRepBlend_RstRstConstRad.cxx L151-157).
    pub fn set_param(&mut self, param: f64) {
        // OCCT L153: tguide->D2(Param, ptgui, d1gui, d2gui).
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        self.d2gui = self.tguide().derivative2_at(param);
        // OCCT L154-156.
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        self.the_d = -(self.nplan.dot(self.ptgui));
    }

    /// OCCT Set(First, Last) (BRepBlend_RstRstConstRad.cxx L161-164) —
    /// `tguide = guide->Trim(First, Last, 1.e-12);`.
    pub fn set_interval(&mut self, first: f64, last: f64) {
        self.tguide = Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(self.guide.clone()),
            first,
            last,
        }));
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BRepBlend_RstRstConstRad.cxx
    /// L168-172).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L170: Tolerance(1) = cons1.Resolution(Tol).
        tolerance[0] = self.cons_resolution(self.surf1, tol);
        // OCCT L171: Tolerance(2) = cons2.Resolution(Tol).
        tolerance[1] = self.cons_resolution(self.surf2, tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BRepBlend_RstRstConstRad.cxx
    /// L176-182) — the cons FirstParameter / LastParameter are the pcurves'
    /// ranges (Adaptor3d_CurveOnSurface L977-987).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.rst1.default_domain()[0]; // cons1.FirstParameter
        inf_bound[1] = self.rst2.default_domain()[0]; // cons2.FirstParameter
        sup_bound[0] = self.rst1.default_domain()[1]; // cons1.LastParameter
        sup_bound[1] = self.rst2.default_domain()[1]; // cons2.LastParameter
    }

    /// OCCT IsSolution(Sol, Tol) (BRepBlend_RstRstConstRad.cxx L186-300).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 2];
        let mut secmember = [0.0f64; 2];
        let mut gradsol = vec![vec![0.0f64; 2]; 2];

        // OCCT L197: Values(Sol, valsol, gradsol).
        self.values(sol, &mut valsol, &mut gradsol);

        // OCCT L199.
        if valsol[0].abs() <= tol && valsol[1].abs() <= tol {
            // Calculation of tangents
            // OCCT L203-206.
            self.prmrst1 = sol[0];
            self.pt2drst1 = self.rst1.point_at(self.prmrst1);
            self.prmrst2 = sol[1];
            self.pt2drst2 = self.rst2.point_at(self.prmrst2);

            // OCCT L208-209.
            let (ptrst1_v, d11) = self.cons1_d1(sol[0]);
            self.ptrst1 = ptrst1_v;
            let (ptrst2_v, d21) = self.cons2_d1(sol[1]);
            self.ptrst2 = ptrst2_v;

            // OCCT L211: dnplan.SetLinearForm(1. / normtg, d2gui,
            // -1. / normtg * (nplan.Dot(d2gui)), nplan).
            let dnplan = (1.0 / self.normtg) * self.d2gui
                + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            // OCCT L213-214.
            let temp = self.ptrst1 - self.ptgui;
            secmember[0] = self.normtg - dnplan.dot(temp);

            // OCCT L216-217.
            let temp = self.ptrst2 - self.ptgui;
            secmember[1] = self.normtg - dnplan.dot(temp);

            // OCCT L219-240: math_Gauss Resol(gradsol) and the SVD fallback.
            let mut a = MatD::new(2, 2);
            for r in 1..=2 {
                for c in 1..=2 {
                    a.set(r, c, gradsol[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                // OCCT L223: Resol.Solve(secmember).
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
                // OCCT L228-239: math_SVD SingRS(gradsol);
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

            // OCCT L242-254.
            if !self.istangent {
                // OCCT L244-245.
                self.tgrst1 = secmember[0] * d11;
                self.tgrst2 = secmember[1] * d21;

                // OCCT L247-253 (gp_Pnt bid unused in the rcad form).
                let (_, d1urst1, d1vrst1) = self.surf1.derivatives(self.pt2drst1.x, self.pt2drst1.y);
                let mut a2d = 0.0f64;
                let mut b2d = 0.0f64;
                t3dto2d(&mut a2d, &mut b2d, self.tgrst1, d1urst1, d1vrst1);
                self.tg2drst1 = DVec2::new(a2d, b2d);
                let (_, d1urst2, d1vrst2) = self.surf2.derivatives(self.pt2drst2.x, self.pt2drst2.y);
                // OCCT L252 passes tgrst1 (not tgrst2) — translated literally.
                t3dto2d(&mut a2d, &mut b2d, self.tgrst1, d1urst2, d1vrst2);
                self.tg2drst2 = DVec2::new(a2d, b2d);
            }

            // OCCT L256-260.
            let mut center = DVec3::ZERO;
            let mut not_used = DVec3::ZERO;
            let is_center = self.center_circle_rst1_rst2(
                self.ptrst1,
                self.ptrst2,
                self.nplan,
                &mut center,
                &mut not_used,
            );

            // OCCT L262-265.
            if !is_center {
                return false;
            }

            // OCCT L267-270.
            let mut n1 = self.ptrst1 - center;
            let mut n2 = self.ptrst2 - center;
            n1 = n1.normalize_or_zero();
            n2 = n2.normalize_or_zero();

            // OCCT L272-273.
            let cosa = n1.dot(n2);
            let mut sina = self.nplan.dot(n1.cross(n2));

            // OCCT L275-278.
            if self.choix % 2 != 0 {
                sina = -sina; // nplan is changed into -nplan
            }

            // OCCT L280-284.
            let mut angle = cosa.acos();
            if sina < 0.0 {
                angle = 2.0 * std::f64::consts::PI - angle;
            }

            // OCCT L286-293: update of maxang / minang.
            if angle > self.maxang {
                self.maxang = angle;
            }
            if angle < self.minang {
                self.minang = angle;
            }
            // OCCT L294.
            self.distmin = self.distmin.min(self.ptrst1.distance(self.ptrst2));

            return true;
        }
        // OCCT L298-299.
        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BRepBlend_RstRstConstRad.cxx L304-307).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT PointOnRst1() (BRepBlend_RstRstConstRad.cxx L311-314).
    pub fn point_on_rst1(&self) -> DVec3 {
        self.ptrst1
    }

    /// OCCT PointOnRst2() (BRepBlend_RstRstConstRad.cxx L318-321).
    pub fn point_on_rst2(&self) -> DVec3 {
        self.ptrst2
    }

    /// OCCT Pnt2dOnRst1() (BRepBlend_RstRstConstRad.cxx L325-328).
    pub fn pnt2d_on_rst1(&self) -> DVec2 {
        self.pt2drst1
    }

    /// OCCT Pnt2dOnRst2() (BRepBlend_RstRstConstRad.cxx L332-335).
    pub fn pnt2d_on_rst2(&self) -> DVec2 {
        self.pt2drst2
    }

    /// OCCT ParameterOnRst1() (BRepBlend_RstRstConstRad.cxx L339-342).
    pub fn parameter_on_rst1(&self) -> f64 {
        self.prmrst1
    }

    /// OCCT ParameterOnRst2() (BRepBlend_RstRstConstRad.cxx L346-349).
    pub fn parameter_on_rst2(&self) -> f64 {
        self.prmrst2
    }

    /// OCCT IsTangencyPoint() (BRepBlend_RstRstConstRad.cxx L353-356).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnRst1() (BRepBlend_RstRstConstRad.cxx L360-367).
    pub fn tangent_on_rst1(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstConstRad::TangentOnRst1");
        }
        self.tgrst1
    }

    /// OCCT Tangent2dOnRst1() (BRepBlend_RstRstConstRad.cxx L371-378).
    pub fn tangent_2d_on_rst1(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstConstRad::Tangent2dOnRst1");
        }
        self.tg2drst1
    }

    /// OCCT TangentOnRst2() (BRepBlend_RstRstConstRad.cxx L382-389).
    pub fn tangent_on_rst2(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstConstRad::TangentOnRst2");
        }
        self.tgrst2
    }

    /// OCCT Tangent2dOnRst2() (BRepBlend_RstRstConstRad.cxx L393-400).
    pub fn tangent_2d_on_rst2(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_RstRstConstRad::Tangent2dOnRst2");
        }
        self.tg2drst2
    }

    /// OCCT Decroch(Sol, NRst1, TgRst1, NRst2, TgRst2)
    /// (BRepBlend_RstRstConstRad.cxx L404-482).
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
        // OCCT L416-417: rstref1->Value(Sol(1)).Coord(u, v);
        // surfref1->D1(u, v, PtTmp1, d1u, d1v).
        let (u, v) = {
            let p = self.rstref1.expect("rstref1").point_at(sol[0]);
            (p.x, p.y)
        };
        let (pt_tmp1, d1u, d1v) = self.surfref1.expect("surfref1").derivatives(u, v);
        // Normal to the reference surface 1
        // OCCT L419: NRst1 = d1u.Crossed(d1v).
        let n_rst1 = d1u.cross(d1v);
        *tgrst1 = n_rst1;
        // OCCT L420-421.
        let (u, v) = {
            let p = self.rstref2.expect("rstref2").point_at(sol[1]);
            (p.x, p.y)
        };
        let (pt_tmp2, d1u, d1v) = self.surfref2.expect("surfref2").derivatives(u, v);
        // Normal to the reference surface 2
        // OCCT L423: NRst2 = d1u.Crossed(d1v).
        let n_rst2 = d1u.cross(d1v);
        *tgrst2 = n_rst2;

        // OCCT L425.
        let mut center = DVec3::ZERO;
        let mut not_used = DVec3::ZERO;
        self.center_circle_rst1_rst2(pt_tmp1, pt_tmp2, self.nplan, &mut center, &mut not_used);

        // OCCT L427-430.
        let norm = self.nplan.cross(n_rst1).length();
        let unsurnorm = 1.0 / norm;
        let mut n_rst1_in_plane =
            (self.nplan.dot(n_rst1) * unsurnorm) * self.nplan + (-unsurnorm) * n_rst1;

        // OCCT L432: centptrst.SetXYZ(PtTmp1.XYZ() - Center.XYZ()).
        let mut centptrst = pt_tmp1 - center;

        // OCCT L434-437.
        if centptrst.dot(n_rst1_in_plane) < 0.0 {
            n_rst1_in_plane = -n_rst1_in_plane;
        }

        // OCCT L439: TgRst1 = nplan.Crossed(centptrst).
        *nrrst1 = self.nplan.cross(centptrst);

        // OCCT L441-451.
        let norm = self.nplan.cross(n_rst2).length();
        let unsurnorm = 1.0 / norm;
        let mut n_rst2_in_plane =
            (self.nplan.dot(n_rst2) * unsurnorm) * self.nplan + (-unsurnorm) * n_rst2;
        centptrst = pt_tmp2 - center;

        // OCCT L446-449.
        if centptrst.dot(n_rst2_in_plane) < 0.0 {
            n_rst2_in_plane = -n_rst2_in_plane;
        }

        // OCCT L451: TgRst2 = nplan.Crossed(centptrst).
        *nrrst2 = self.nplan.cross(centptrst);

        // OCCT L453-457.
        if self.choix % 2 != 0 {
            *nrrst1 = -*nrrst1;
            *nrrst2 = -*nrrst2;
        }

        // The vectors are returned
        // OCCT L460-481.
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

    /// OCCT Set(Radius, Choix) (BRepBlend_RstRstConstRad.cxx L486-490).
    pub fn set(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        self.ray = radius.abs();
    }

    /// OCCT Set(BlendFunc_SectionShape) (BRepBlend_RstRstConstRad.cxx
    /// L494-497).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_s_shape = type_section;
    }

    /// OCCT CenterCircleRst1Rst2(PtRst1, PtRst2, np, Center, VdMed)
    /// (BRepBlend_RstRstConstRad.cxx L503-542) — calculates the center of
    /// the circle passing by two points of restrictions.
    pub fn center_circle_rst1_rst2(
        &self,
        pt_rst1: DVec3,
        pt_rst2: DVec3,
        np: DVec3,
        center: &mut DVec3,
        vd_med: &mut DVec3,
    ) -> bool {
        // OCCT L510: gp_Vec rst1rst2(PtRst1, PtRst2).
        let rst1rst2 = pt_rst2 - pt_rst1;

        // Calculate the center of the circle
        // OCCT L516-518.
        *vd_med = rst1rst2.cross(np);
        let norm2 = rst1rst2.length_squared();
        let mut dist = self.ray * self.ray - 0.25 * norm2;

        // OCCT L520-523.
        if self.choix > 2 {
            *vd_med = -*vd_med;
        }

        // OCCT L525-528.
        if dist < -1.0e-07 {
            return false;
        }

        // OCCT L530-539.
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
    /// (BRepBlend_RstRstConstRad.cxx L546-586).
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
        // OCCT L556-557: tguide->D1(Param, ptgui, d1gui); np = d1gui.Normalized().
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        let mut np = self.d1gui.normalize_or_zero();
        // OCCT L558-559.
        self.ptrst1 = self.cons1_value(u);
        self.ptrst2 = self.cons2_value(v);

        // OCCT L561 (the IsCenter flag is discarded in the original).
        let mut center = DVec3::ZERO;
        let mut not_used = DVec3::ZERO;
        self.center_circle_rst1_rst2(self.ptrst1, self.ptrst2, np, &mut center, &mut not_used);

        // OCCT L563-564.
        c.radius = self.ray.abs();
        let ns = (self.ptrst1 - center).normalize_or_zero();

        // OCCT L566-569.
        if self.choix % 2 != 0 {
            np = -np;
        }

        // OCCT L571: C.SetPosition(gp_Ax2(Center, np, ns))
        // (gp_Ax2(P, N, Vx): the Y direction is N ^ Vx).
        c.center = center;
        c.normal = np;
        c.x_dir = ns;
        c.y_dir = np.cross(ns);
        // OCCT L572-573: Pdeb = 0; Pfin = ElCLib::Parameter(C, ptrst2).
        *pdeb = 0.0;
        *pfin = elclib_circle_parameter(c, self.ptrst2);

        // Test of angles negative and almost null : Special Case
        // OCCT L576-581.
        if *pfin > 1.5 * std::f64::consts::PI {
            np = -np;
            c.normal = np;
            c.y_dir = np.cross(ns);
            *pfin = elclib_circle_parameter(c, self.ptrst2);
        }
        // OCCT L582-585.
        if *pfin < p_confusion() {
            *pfin += p_confusion();
        }
    }

    /// OCCT IsRational() (BRepBlend_RstRstConstRad.cxx L590-593).
    pub fn is_rational(&self) -> bool {
        self.my_s_shape == BlendFuncSectionShape::Rational
            || self.my_s_shape == BlendFuncSectionShape::QuasiAngular
    }

    /// OCCT GetSectionSize() (BRepBlend_RstRstConstRad.cxx L597-600).
    pub fn get_section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// OCCT GetMinimalWeight(Weights) (BRepBlend_RstRstConstRad.cxx
    /// L604-608).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        blend_func_get_minimal_weights(self.my_s_shape, self.my_t_conv, self.minang, self.maxang, weigths);
        // It is supposed that it does not depend on the Radius!
    }

    /// OCCT NbIntervals(S) (BRepBlend_RstRstConstRad.cxx L612-615).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.guide.nb_intervals(blend_func_next_shape(s))
    }

    /// OCCT Intervals(T, S) (BRepBlend_RstRstConstRad.cxx L619-622).
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let mut intervals = Vec::new();
        self.guide.intervals(&mut intervals, blend_func_next_shape(s));
        for (dst, src) in t.iter_mut().zip(intervals) {
            *dst = src;
        }
    }

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
    /// (BRepBlend_RstRstConstRad.cxx L626-630).
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
    /// (BRepBlend_RstRstConstRad.cxx L637-650) — tolerances used for
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

    /// OCCT Knots(TKnots) (BRepBlend_RstRstConstRad.cxx L654-657).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        geomfill_knots(self.my_t_conv, tknots);
    }

    /// OCCT Mults(TMults) (BRepBlend_RstRstConstRad.cxx L661-664).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        geomfill_mults(self.my_t_conv, tmults);
    }

    /// OCCT Section(P, Poles, Poles2d, Weights)
    /// (BRepBlend_RstRstConstRad.cxx L668-729).
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

        // OCCT L681-682: tguide->D1(prm, ptgui, d1gui); nplan = d1gui.Normalized().
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        self.nplan = self.d1gui.normalize_or_zero();

        // OCCT L684-685.
        let u = p.parameter_on_c1();
        let v = p.parameter_on_c2();

        // OCCT L687-688.
        let pt2d1 = self.rst1.point_at(u);
        let pt2d2 = self.rst2.point_at(v);

        // OCCT L690-692.
        self.ptrst1 = self.cons1_value(u);
        self.ptrst2 = self.cons2_value(v);
        self.distmin = self.distmin.min(self.ptrst1.distance(self.ptrst2));

        // OCCT L694-695.
        poles_2d[0] = DVec2::new(pt2d1.x, pt2d1.y);
        poles_2d[poles_2d.len() - 1] = DVec2::new(pt2d2.x, pt2d2.y);

        // Linear case
        // OCCT L698-705.
        if self.my_s_shape == BlendFuncSectionShape::Linear {
            poles[low] = self.ptrst1;
            poles[upp] = self.ptrst2;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            return;
        }

        // Calculate the center of the circle
        // OCCT L708 (the IsCenter flag is discarded in the original).
        let mut center = DVec3::ZERO;
        let mut not_used = DVec3::ZERO;
        self.center_circle_rst1_rst2(self.ptrst1, self.ptrst2, self.nplan, &mut center, &mut not_used);

        // normals to the section with points
        // OCCT L711-712.
        let ns = (self.ptrst1 - center).normalize_or_zero();
        let ns2 = (self.ptrst2 - center).normalize_or_zero();

        // OCCT L714-717.
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
        }

        // OCCT L719-728: GeomFill::GetCircle(myTConv, ns, ns2, nplan, ptrst1,
        // ptrst2, std::abs(ray), Center, Poles, Weights).
        get_circle(
            tconv(self.my_t_conv),
            ns,
            ns2,
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
    /// (BRepBlend_RstRstConstRad.cxx L733-948) — used for the first and
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

        // OCCT L758-761.
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        self.d2gui = self.tguide().derivative2_at(prm);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        // OCCT L761: dnplan.SetLinearForm(1. / normtg, d2gui,
        // -1. / normtg * (nplan.Dot(d2gui)), nplan).
        let mut dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        // OCCT L763-766.
        let mut sol = [0.0f64; 2];
        self.prmrst1 = p.parameter_on_c1();
        sol[0] = self.prmrst1;
        self.prmrst2 = p.parameter_on_c2();
        sol[1] = self.prmrst2;
        self.pt2drst1 = self.rst1.point_at(self.prmrst1);
        self.pt2drst2 = self.rst2.point_at(self.prmrst2);

        // OCCT L768: Values(sol, valsol, gradsol).
        let mut valsol = [0.0f64; 2];
        let mut gradsol = vec![vec![0.0f64; 2]; 2];
        self.values(&sol, &mut valsol, &mut gradsol);

        // OCCT L770-771.
        let (ptrst1_v, d11) = self.cons1_d1(sol[0]);
        self.ptrst1 = ptrst1_v;
        let (ptrst2_v, d21) = self.cons2_d1(sol[1]);
        self.ptrst2 = ptrst2_v;

        // OCCT L773-777.
        let mut secmember = [0.0f64; 2];
        let temp = self.ptrst1 - self.ptgui;
        secmember[0] = self.normtg - dnplan.dot(temp);

        let temp = self.ptrst2 - self.ptgui;
        secmember[1] = self.normtg - dnplan.dot(temp);

        // OCCT L779-800: math_Gauss Resol(gradsol, 1.e-9) and the SVD
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
            // OCCT L784: Resol.Solve(secmember).
            let mut x = VecD::new(2);
            for i in 1..=2 {
                x.set(i, secmember[i - 1]);
            }
            resol.solve(&mut x);
            for i in 1..=2 {
                secmember[i - 1] = x.get(i);
            }
        } else {
            // OCCT L788-799: math_SVD SingRS(gradsol);
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

        // OCCT L802-816.
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
        let mut n1 = self.ptrst1 - center;
        let mut n2 = self.ptrst2 - center;
        n1 = n1.normalize_or_zero();
        n2 = n2.normalize_or_zero();

        // OCCT L818-864: secmember contains derivatives of parameters on
        // curves compared to t.
        let mut d1n1 = DVec3::ZERO; // only read on the !istgt path
        let mut d1n2 = DVec3::ZERO; // only read on the !istgt path
        if !istgt {
            // OCCT L822-823.
            self.tgrst1 = secmember[0] * d11;
            self.tgrst2 = secmember[1] * d21;

            // OCCT L825-829.
            let mut d1rst1rst2;
            let norm2 = rst1rst2.length_squared();
            d1rst1rst2 = self.tgrst2 - self.tgrst1;
            let mut dist = self.ray * self.ray - 0.25 * norm2;

            // OCCT L831-854.
            if dist > 1.0e-07 {
                // OCCT L833-838.
                let d1p1p2crosnp = d1rst1rst2.cross(self.nplan) + rst1rst2.cross(dnplan);
                // derivative of the perpendicular bisector
                let mut dmed = d1p1p2crosnp - med.dot(d1p1p2crosnp) * med;
                dmed /= normmed;
                dist = dist.sqrt();
                let d1dist = -(0.25 / dist) * rst1rst2.dot(d1rst1rst2);

                // OCCT L841-844.
                if self.choix > 2 {
                    dmed = -dmed;
                }

                // the derivative of coefficient Dist is located in dmed
                // OCCT L847-853.
                dmed = dist * dmed + d1dist * med;
                d1rst1rst2 *= 0.5;
                // derivative of the Normal to the curve in P1
                d1n1 = -(dmed + d1rst1rst2) / self.ray;

                // derivative of the Normal to the curve in P2
                d1n2 = (d1rst1rst2 - dmed) / self.ray;
            } else {
                // OCCT L857-863.
                d1rst1rst2 *= 0.5;
                // Normal to the curve in P1
                d1n1 = -d1rst1rst2 / self.ray;

                // Normal to the curve in P2
                d1n2 = d1rst1rst2 / self.ray;
            }
        }

        // Tops 2d
        // OCCT L868-880.
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

        // Linear case
        // OCCT L883-897.
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
        // OCCT L901-904.
        let tgct: DVec3;
        if !istgt {
            tgct = -self.ray * d1n1 + self.tgrst1;
        } else {
            tgct = DVec3::ZERO; // the OCCT tgct is only read on the !istgt path
        }

        // OCCT L906-910.
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
            dnplan = -dnplan;
        }

        // OCCT L912-947.
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
                self.ray.abs(),
                0.0,
                center,
                tgct,
                poles,
                d_poles,
                weigths,
                d_weigths,
            )
        } else {
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
            false
        }
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weights, DWeights, D2Weights) (BRepBlend_RstRstConstRad.cxx
    /// L952-964) — returns false.
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
    /// (BRepBlend_RstRstConstRad.cxx L966-981).
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

impl<'a> FunctionSetWithDerivatives for BlendRstRstConstRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendRstRstConstRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendRstRstConstRad::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendRstRstConstRad::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendRstRstConstRad::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendRstRstConstRad::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendRstRstConstRad<'a> {
    fn set_param(&mut self, param: f64) {
        BlendRstRstConstRad::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendRstRstConstRad::set_interval(self, first, last)
    }

    fn pnt1(&self) -> DVec3 {
        BlendRstRstFunction::pnt1(self)
    }

    fn pnt2(&self) -> DVec3 {
        BlendRstRstFunction::pnt2(self)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendRstRstConstRad::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendRstRstConstRad::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendRstRstConstRad::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendRstRstConstRad::get_minimal_distance(self)
    }

    fn is_rational(&self) -> bool {
        BlendRstRstConstRad::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendRstRstConstRad::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendRstRstConstRad::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendRstRstConstRad::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendRstRstConstRad::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendRstRstConstRad::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendRstRstConstRad::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendRstRstConstRad::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendRstRstConstRad::mults(self, tmults)
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
        BlendRstRstConstRad::section_d1(
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
        BlendRstRstConstRad::section_simple(self, p, poles, poles_2d, weigths)
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
        BlendRstRstConstRad::section_d2(
            self, p, poles, d_poles, d2_poles, poles_2d, d_poles_2d, d2_poles_2d, weigths,
            d_weigths, d2_weigths,
        )
    }

    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        BlendRstRstConstRad::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendRstRstFunction for BlendRstRstConstRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendRstRstConstRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendRstRstConstRad::nb_equations(self)
    }

    fn point_on_rst1(&self) -> DVec3 {
        BlendRstRstConstRad::point_on_rst1(self)
    }

    fn point_on_rst2(&self) -> DVec3 {
        BlendRstRstConstRad::point_on_rst2(self)
    }

    fn pnt2d_on_rst1(&self) -> DVec2 {
        BlendRstRstConstRad::pnt2d_on_rst1(self)
    }

    fn pnt2d_on_rst2(&self) -> DVec2 {
        BlendRstRstConstRad::pnt2d_on_rst2(self)
    }

    fn parameter_on_rst1(&self) -> f64 {
        BlendRstRstConstRad::parameter_on_rst1(self)
    }

    fn parameter_on_rst2(&self) -> f64 {
        BlendRstRstConstRad::parameter_on_rst2(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendRstRstConstRad::is_tangency_point(self)
    }

    fn tangent_on_rst1(&self) -> DVec3 {
        BlendRstRstConstRad::tangent_on_rst1(self)
    }

    fn tangent_2d_on_rst1(&self) -> DVec2 {
        BlendRstRstConstRad::tangent_2d_on_rst1(self)
    }

    fn tangent_on_rst2(&self) -> DVec3 {
        BlendRstRstConstRad::tangent_on_rst2(self)
    }

    fn tangent_2d_on_rst2(&self) -> DVec2 {
        BlendRstRstConstRad::tangent_2d_on_rst2(self)
    }

    fn decroch(
        &self,
        sol: &[f64],
        tgrst1: &mut DVec3,
        nrrst1: &mut DVec3,
        tgrst2: &mut DVec3,
        nrrst2: &mut DVec3,
    ) -> BlendDecrochStatus {
        BlendRstRstConstRad::decroch(self, sol, tgrst1, nrrst1, tgrst2, nrrst2)
    }
}
