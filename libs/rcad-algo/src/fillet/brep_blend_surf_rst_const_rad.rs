//! OCCT BRepBlend_SurfRstConstRad (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_SurfRstConstRad.hxx (L40-244) + BRepBlend_SurfRstConstRad.cxx
//! (whole file L37-1048).  Copy of CSConstRad with pcurve on surface as
//! support (hxx L38-39); the SimulSurf / PerformSurf face/rst constant-radius
//! arms instantiate it as `func` (ChFi3d_FilBuilder.cxx L833 / L1763).
//!
//! Architecture mappings (mirroring [`super::brep_blend_func_consrad`]):
//! `class BRepBlend_SurfRstConstRad : public Blend_SurfRstFunction` is
//! expressed by implementing the [`BlendSurfRstFunction`] and
//! [`BlendAppFunction`] traits over the `math_FunctionSetWithDerivatives`
//! base; `occ::handle<Adaptor3d_Curve> tguide` (aliasing `guide` until
//! Set(First, Last) trims it) maps to an owned trimmed [`Curve3`] copy plus
//! an accessor; `math_Vector` / `math_Matrix` map to `[f64; 3]` /
//! `Vec<Vec<f64>>` (OCCT X(i) -> x[i - 1], D(i, j) -> d[i - 1][j - 1]);
//! `gp_Circ` maps to the kernel [`Circle3`].
//!
//! The OCCT `Adaptor3d_CurveOnSurface cons` member (hxx L214) has no rcad
//! adaptor equivalent; the rcad port stores the (rst, surfrst) pair the
//! adaptor wraps and transcribes the consumed operations from
//! Adaptor3d_CurveOnSurface.cxx (Value/EvalD0, D1/EvalD1 generic branch,
//! FirstParameter/LastParameter L977-987, Resolution L1364-1370) in the
//! `cons_*` helpers below.
//!
//! Pending kernel dependency (marked GAP, plan 0.6): Adaptor2d_Curve2d::
//! Resolution (the last step of the cons Resolution chain). math_SVD is
//! available (`rcad_kernel::math::math_svd::MathSvd`) and used by the
//! IsSolution / Section-d1 second-chance solvers.

use glam::{DVec2, DVec3};

use rcad_kernel::base::convert::ConvertParameterisation;
use rcad_kernel::core::precision::{is_infinite_value, p_confusion};
use rcad_kernel::geom::{Circle3, Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::math_svd::MathSvd;
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::{MatD, VecD};

use crate::geomalgo::geomfill::geom_fill::{get_circle, get_circle_d1};

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
use super::brep_blend_surf_rst_function::BlendSurfRstFunction;

/// OCCT Eps constant (BRepBlend_SurfRstConstRad.cxx L35).
const EPS: f64 = 1.0e-15;

/// OCCT static t3dto2d (BRepBlend_SurfRstConstRad.cxx L37-47).
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
/// (pattern of brep_blend_func_evolrad_b::tconv).
fn tconv(t_conv: ConvertParameterisationType) -> ConvertParameterisation {
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

/// OCCT BRepBlend_SurfRstConstRad — function to build a variable evolutive
/// section : constant radius rolling ball on a surface and a pcurve on
/// another surface.
pub struct BlendSurfRstConstRad<'a> {
    // OCCT BRepBlend_SurfRstConstRad.hxx fields (L211-241).
    pub(crate) surf: &'a Surface3,
    pub(crate) surfrst: &'a Surface3,
    pub(crate) rst: &'a Curve2d,
    // OCCT: Adaptor3d_CurveOnSurface cons(Rst, SurfRst) — carried as the
    // wrapped pair (see the module header architecture note).
    pub(crate) guide: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> tguide — a handle aliasing guide
    // until Set(First, Last) trims it (cxx L272).  The rcad port owns the
    // trimmed copy.
    pub(crate) tguide: Option<Curve3>,
    pub(crate) pts: DVec3,
    pub(crate) ptrst: DVec3,
    pub(crate) pt2ds: DVec2,
    pub(crate) pt2drst: DVec2,
    pub(crate) prmrst: f64,
    pub(crate) istangent: bool,
    pub(crate) tgs: DVec3,
    pub(crate) tg2ds: DVec2,
    pub(crate) tgrst: DVec3,
    pub(crate) tg2drst: DVec2,
    pub(crate) ray: f64,
    pub(crate) choix: i32,
    pub(crate) ptgui: DVec3,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    // OCCT: occ::handle<Adaptor3d_Surface> surfref — null until Set(SurfRef,
    // RstRef) (cxx L246-251).
    pub(crate) surfref: Option<&'a Surface3>,
    pub(crate) rstref: Option<&'a Curve2d>,
    pub(crate) maxang: f64,
    pub(crate) minang: f64,
    pub(crate) distmin: f64,
    pub(crate) my_s_shape: BlendFuncSectionShape,
    pub(crate) my_t_conv: ConvertParameterisationType,
}

impl<'a> BlendSurfRstConstRad<'a> {
    /// Architecture mapping: OCCT `tguide` is a handle aliasing `guide`
    /// until Set(First, Last) replaces it with a trimmed copy (cxx L272).
    #[inline]
    pub(crate) fn tguide(&self) -> &Curve3 {
        self.tguide.as_ref().unwrap_or(self.guide)
    }

    /// OCCT Adaptor3d_CurveOnSurface::Value (the EvalD0 generic branch of
    /// Adaptor3d_CurveOnSurface.cxx — the pcurve point carried on the
    /// wrapped surface).
    fn cons_value(&self, w: f64) -> DVec3 {
        let puv = self.rst.point_at(w);
        self.surfrst.point_at(puv.x, puv.y)
    }

    /// OCCT Adaptor3d_CurveOnSurface::D1 (Adaptor3d_CurveOnSurface.cxx
    /// EvalD1 L1200-1241, the generic branch: the pcurve D1 composed with
    /// the carrier surface D1 — `D1.SetLinearForm(Duv.X(), D1U, Duv.Y(),
    /// D1V)`).
    fn cons_d1(&self, w: f64) -> (DVec3, DVec3) {
        let puv = self.rst.point_at(w);
        let duv = self.rst.derivative_at(w);
        let (p, d1u, d1v) = self.surfrst.derivatives(puv.x, puv.y);
        (p, duv.x * d1u + duv.y * d1v)
    }

    /// OCCT Adaptor3d_CurveOnSurface::Resolution (Adaptor3d_CurveOnSurface.cxx
    /// L1364-1370): `ru = mySurface->UResolution(R3d); rv = ...; return
    /// myCurve->Resolution(std::min(ru, rv));`  GAP (plan 0.6): the final
    /// Adaptor2d_Curve2d::Resolution step has no rcad equivalent; the
    /// established pending marker preserves the OCCT failure path.
    fn cons_resolution(&self, r3d: f64) -> f64 {
        let ru = self.surfrst.u_resolution(r3d);
        let rv = self.surfrst.v_resolution(r3d);
        let _ = ru.min(rv);
        adaptor2d_curve2d_resolution_pending()
    }

    /// OCCT BRepBlend_SurfRstConstRad(Surf, SurfRst, Rst, CGuide)
    /// (BRepBlend_SurfRstConstRad.cxx L51-72).
    pub fn new(
        surf: &'a Surface3,
        surfrst: &'a Surface3,
        rst: &'a Curve2d,
        cguide: &'a Curve3,
    ) -> Self {
        BlendSurfRstConstRad {
            surf,
            surfrst,
            rst,
            guide: cguide,
            tguide: None, // OCCT: tguide(CGuide) — alias, see tguide().
            pts: DVec3::ZERO,
            ptrst: DVec3::ZERO,
            pt2ds: DVec2::ZERO,
            pt2drst: DVec2::ZERO,
            prmrst: 0.0,
            istangent: true,
            tgs: DVec3::ZERO,
            tg2ds: DVec2::ZERO,
            tgrst: DVec3::ZERO,
            tg2drst: DVec2::ZERO,
            ray: 0.0,
            choix: 0,
            ptgui: DVec3::ZERO,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            nplan: DVec3::ZERO,
            normtg: 0.0,
            the_d: 0.0,
            surfref: None,
            rstref: None,
            maxang: -f64::MAX, // OCCT: RealFirst()
            minang: f64::MAX,  // OCCT: RealLast()
            distmin: f64::MAX, // OCCT: RealLast()
            my_s_shape: BlendFuncSectionShape::Rational,
            my_t_conv: ConvertParameterisationType::TgtThetaOver2,
        }
    }

    /// OCCT NbVariables() (BRepBlend_SurfRstConstRad.cxx L76-79) — returns 3.
    pub fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations() (BRepBlend_SurfRstConstRad.cxx L83-86) — returns 3.
    pub fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(X, F) (BRepBlend_SurfRstConstRad.cxx L90-108).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L95: surf->D1(X(1), X(2), pts, d1u1, d1v1);
        let (p, d1u1, d1v1) = self.surf.derivatives(x[0], x[1]);
        self.pts = p;
        // OCCT L96: ptrst = cons.Value(X(3));
        self.ptrst = self.cons_value(x[2]);

        // OCCT L98-99.
        f[0] = self.nplan.dot(self.pts) + self.the_d;
        f[1] = self.nplan.dot(self.ptrst) + self.the_d;

        // OCCT L101-103: ns and its projection on the plane.
        let mut ns = d1u1.cross(d1v1);
        let norm = self.nplan.cross(ns).length();
        ns = (self.nplan.dot(ns) / norm) * self.nplan + (-1.0 / norm) * ns;
        // OCCT L104-105: vref.SetLinearForm(ray, ns, gp_Vec(ptrst, pts));
        // vref /= ray;
        let mut vref = self.ray * ns + (self.pts - self.ptrst);
        vref /= self.ray;
        // OCCT L106.
        f[2] = (vref.length_squared() - 1.0) * self.ray * self.ray;
        true
    }

    /// OCCT Derivatives(X, D) (BRepBlend_SurfRstConstRad.cxx L112-171).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L119: surf->D2(X(1), X(2), pts, d1u1, d1v1, d2u1, d2v1, d2uv1);
        let (p, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.derivatives2(x[0], x[1]);
        self.pts = p;
        // OCCT L120: cons.D1(X(3), ptrst, d1);
        let (ptrst_v, d1) = self.cons_d1(x[2]);
        self.ptrst = ptrst_v;

        // OCCT L122-128.
        d[0][0] = self.nplan.dot(d1u1);
        d[0][1] = self.nplan.dot(d1v1);
        d[0][2] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(d1);

        // OCCT L130-137.
        let ns = d1u1.cross(d1v1);
        let ncrossns = self.nplan.cross(ns);
        let norm = ncrossns.length();
        let ndotns = self.nplan.dot(ns);

        // OCCT L135-137: vref.SetLinearForm(ndotns, nplan, -1., ns);
        // vref.Divide(norm); vref.SetLinearForm(ray, vref, gp_Vec(ptrst, pts));
        let vref = (ndotns * self.nplan + (-1.0) * ns) / norm;
        let vref = self.ray * vref + (self.pts - self.ptrst);

        // Derivative by u1
        // OCCT L140-148.
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1u1;

        // OCCT L150-151.
        d[2][0] = resul.dot(vref);
        d[2][0] *= 2.0;

        // Derivative by v1
        // OCCT L154-162.
        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1v1;

        // OCCT L164-165.
        d[2][1] = resul.dot(vref);
        d[2][1] *= 2.0;

        // OCCT L167-168.
        d[2][2] = d1.dot(vref);
        d[2][2] *= -2.0;

        true
    }

    /// OCCT Values(X, F, D) (BRepBlend_SurfRstConstRad.cxx L175-242).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L183: surf->D2(X(1), X(2), pts, d1u1, d1v1, d2u1, d2v1, d2uv1);
        let (p, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.derivatives2(x[0], x[1]);
        self.pts = p;
        // OCCT L184: cons.D1(X(3), ptrst, d1);
        let (ptrst_v, d1) = self.cons_d1(x[2]);
        self.ptrst = ptrst_v;

        // OCCT L186-187.
        f[0] = self.nplan.dot(self.pts) + self.the_d;
        f[1] = self.nplan.dot(self.ptrst) + self.the_d;

        // OCCT L189-195.
        d[0][0] = self.nplan.dot(d1u1);
        d[0][1] = self.nplan.dot(d1v1);
        d[0][2] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(d1);

        // OCCT L197-204.
        let ns = d1u1.cross(d1v1);
        let ncrossns = self.nplan.cross(ns);
        let norm = ncrossns.length();
        let ndotns = self.nplan.dot(ns);

        // OCCT L202-204.
        let vref = (ndotns * self.nplan + (-1.0) * ns) / norm;
        let vref = self.ray * vref + (self.pts - self.ptrst);

        // OCCT L206-208:
        // temp = vref / ray;
        // F(3) = (temp.SquareMagnitude() - 1) * ray * ray; // more stable numerically
        let temp = vref / self.ray;
        f[2] = (temp.length_squared() - 1.0) * self.ray * self.ray;

        // Derivative by u1
        // OCCT L211-219.
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1u1;

        // OCCT L221-222.
        d[2][0] = resul.dot(vref);
        d[2][0] *= 2.0;

        // Derivative by v1
        // OCCT L225-233.
        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1v1;

        // OCCT L235-236.
        d[2][1] = resul.dot(vref);
        d[2][1] *= 2.0;

        // OCCT L238-239.
        d[2][2] = d1.dot(vref);
        d[2][2] *= -2.0;

        true
    }

    /// OCCT Set(SurfRef, RstRef) (BRepBlend_SurfRstConstRad.cxx L246-251).
    pub fn set_ref(&mut self, surf_ref: &'a Surface3, rst_ref: &'a Curve2d) {
        self.surfref = Some(surf_ref);
        self.rstref = Some(rst_ref);
    }

    /// OCCT Set(Param) (BRepBlend_SurfRstConstRad.cxx L255-266).
    pub fn set_param(&mut self, param: f64) {
        self.d1gui = DVec3::ZERO;
        self.nplan = DVec3::ZERO;
        // OCCT L259: tguide->D2(Param, ptgui, d1gui, d2gui);
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        self.d2gui = self.tguide().derivative2_at(param);
        // OCCT L260-261.
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        // OCCT L262-265: theD = nplan.Dot(ptgui); theD = theD * (-1.).
        self.the_d = self.nplan.dot(self.ptgui);
        self.the_d *= -1.0;
    }

    /// OCCT Set(First, Last) (BRepBlend_SurfRstConstRad.cxx L270-273) —
    /// `tguide = guide->Trim(First, Last, 1.e-12);`.
    pub fn set_interval(&mut self, first: f64, last: f64) {
        self.tguide = Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(self.guide.clone()),
            first,
            last,
        }));
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BRepBlend_SurfRstConstRad.cxx
    /// L277-282).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L279-280.
        tolerance[0] = self.surf.u_resolution(tol);
        tolerance[1] = self.surf.v_resolution(tol);
        // OCCT L281: Tolerance(3) = cons.Resolution(Tol).
        tolerance[2] = self.cons_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BRepBlend_SurfRstConstRad.cxx
    /// L286-307).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        // OCCT L288-293 (cons.FirstParameter / LastParameter = the pcurve's).
        inf_bound[0] = self.surf.default_domain()[0]; // FirstUParameter
        inf_bound[1] = self.surf.default_domain()[2]; // FirstVParameter
        inf_bound[2] = self.rst.default_domain()[0]; // cons.FirstParameter
        sup_bound[0] = self.surf.default_domain()[1]; // LastUParameter
        sup_bound[1] = self.surf.default_domain()[3]; // LastVParameter
        sup_bound[2] = self.rst.default_domain()[1]; // cons.LastParameter

        // OCCT L295-300.
        if !is_infinite_value(inf_bound[0]) && !is_infinite_value(sup_bound[0]) {
            let range = sup_bound[0] - inf_bound[0];
            inf_bound[0] -= range;
            sup_bound[0] += range;
        }
        // OCCT L301-306.
        if !is_infinite_value(inf_bound[1]) && !is_infinite_value(sup_bound[1]) {
            let range = sup_bound[1] - inf_bound[1];
            inf_bound[1] -= range;
            sup_bound[1] += range;
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BRepBlend_SurfRstConstRad.cxx L311-426).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT L322: Values(Sol, valsol, gradsol);
        let mut valsol = [0.0f64; 3];
        let mut gradsol = vec![vec![0.0f64; 3]; 3];
        self.values(sol, &mut valsol, &mut gradsol);
        // OCCT L323-325.
        if valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= 2.0 * tol * self.ray.abs()
        {
            // Calculation of tangents

            // OCCT L329-334.
            self.pt2ds = DVec2::new(sol[0], sol[1]);
            self.prmrst = sol[2];
            self.pt2drst = self.rst.point_at(self.prmrst);
            let (p, d1u1, d1v1) = self.surf.derivatives(sol[0], sol[1]);
            self.pts = p;
            let (ptrst_v, d1) = self.cons_d1(sol[2]);
            self.ptrst = ptrst_v;
            // OCCT L334: dnplan.SetLinearForm(1. / normtg, d2gui,
            // -1. / normtg * (nplan.Dot(d2gui)), nplan);
            let dnplan = (1.0 / self.normtg) * self.d2gui
                + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            // OCCT L336-337.
            let temp = self.pts - self.ptgui;
            let mut secmember = [0.0f64; 3];
            secmember[0] = self.normtg - dnplan.dot(temp);

            // OCCT L339-340.
            let temp = self.ptrst - self.ptgui;
            secmember[1] = self.normtg - dnplan.dot(temp);

            // OCCT L342-353.
            let mut ns = d1u1.cross(d1v1);
            let ncrossns = self.nplan.cross(ns);
            let ndotns = self.nplan.dot(ns);
            let norm = ncrossns.length();

            let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
            let temp = (self.ray / norm * (dnplan.dot(ns) - grosterme * ndotns)) * self.nplan
                + (self.ray * ndotns / norm) * dnplan
                + (self.ray * grosterme / norm) * ns;

            // OCCT L355-357.
            ns = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;
            let resul = self.ray * ns + (self.pts - self.ptrst);
            secmember[2] = -2.0 * temp.dot(resul);

            // OCCT L359-364: math_Gauss Resol(gradsol); if (Resol.IsDone())
            // { Resol.Solve(secmember); istangent = false; }
            let mut a = MatD::new(3, 3);
            for r in 1..=3 {
                for c in 1..=3 {
                    a.set(r, c, gradsol[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                let mut x = VecD::new(3);
                for i in 1..=3 {
                    x.set(i, secmember[i - 1]);
                }
                resol.solve(&mut x);
                for i in 1..=3 {
                    secmember[i - 1] = x.get(i);
                }
                self.istangent = false;
            } else {
                // OCCT L366-378: math_SVD SingRS(gradsol);
                // if (SingRS.IsDone()) { math_Vector DEDT(1,3);
                // DEDT = secmember; SingRS.Solve(DEDT, secmember, 1.e-6);
                // istangent = false; } else { istangent = true; }
                let mut sing_rs = MathSvd::new(&a);
                if sing_rs.is_done() {
                    let mut dedt = VecD::new(3);
                    let mut sol = VecD::new(3);
                    for i in 1..=3 {
                        dedt.set(i, secmember[i - 1]);
                        sol.set(i, secmember[i - 1]);
                    }
                    sing_rs.solve(&dedt, &mut sol, 1.0e-6);
                    for i in 1..=3 {
                        secmember[i - 1] = sol.get(i);
                    }
                    self.istangent = false;
                } else {
                    self.istangent = true;
                }
            }

            // OCCT L381-390.
            if !self.istangent {
                // OCCT L383: tgs.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1);
                self.tgs = secmember[0] * d1u1 + secmember[1] * d1v1;
                // OCCT L384: tgrst = secmember(3) * d1;
                self.tgrst = secmember[2] * d1;
                // OCCT L385.
                self.tg2ds = DVec2::new(secmember[0], secmember[1]);
                // OCCT L386-389.
                let (_, d1urst, d1vrst) = self.surfrst.derivatives(self.pt2drst.x, self.pt2drst.y);
                let mut a2d = 0.0f64;
                let mut b2d = 0.0f64;
                t3dto2d(&mut a2d, &mut b2d, self.tgrst, d1urst, d1vrst);
                self.tg2drst = DVec2::new(a2d, b2d);
            }

            // update of maxang
            // OCCT L393-396.
            if self.ray > 0.0 {
                ns = -ns;
            }
            // OCCT L397: ns2 = -resul.Normalized();
            let ns2 = -resul.normalize_or_zero();

            // OCCT L399-404.
            let cosa = ns.dot(ns2);
            let mut sina = self.nplan.dot(ns.cross(ns2));
            if self.choix % 2 != 0 {
                sina = -sina; // nplan is changed to -nplan
            }

            // OCCT L406-410.
            let angle = cosa.acos();
            let angle = if sina < 0.0 {
                2.0 * std::f64::consts::PI - angle
            } else {
                angle
            };

            // OCCT L412-419.
            if angle > self.maxang {
                self.maxang = angle;
            }
            if angle < self.minang {
                self.minang = angle;
            }
            // OCCT L420.
            self.distmin = self.distmin.min(self.pts.distance(self.ptrst));

            return true;
        }
        // OCCT L424-425.
        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BRepBlend_SurfRstConstRad.cxx L430-433).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT PointOnS() (BRepBlend_SurfRstConstRad.cxx L437-440).
    pub fn point_on_s(&self) -> DVec3 {
        self.pts
    }

    /// OCCT PointOnRst() (BRepBlend_SurfRstConstRad.cxx L444-447).
    pub fn point_on_rst(&self) -> DVec3 {
        self.ptrst
    }

    /// OCCT Pnt2dOnS() (BRepBlend_SurfRstConstRad.cxx L451-454).
    pub fn pnt2d_on_s(&self) -> DVec2 {
        self.pt2ds
    }

    /// OCCT Pnt2dOnRst() (BRepBlend_SurfRstConstRad.cxx L458-461).
    pub fn pnt2d_on_rst(&self) -> DVec2 {
        self.pt2drst
    }

    /// OCCT ParameterOnRst() (BRepBlend_SurfRstConstRad.cxx L465-468).
    pub fn parameter_on_rst(&self) -> f64 {
        self.prmrst
    }

    /// OCCT IsTangencyPoint() (BRepBlend_SurfRstConstRad.cxx L472-475).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnS() (BRepBlend_SurfRstConstRad.cxx L479-486).
    pub fn tangent_on_s(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstConstRad::TangentOnS");
        }
        self.tgs
    }

    /// OCCT Tangent2dOnS() (BRepBlend_SurfRstConstRad.cxx L490-497).
    pub fn tangent_2d_on_s(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstConstRad::Tangent2dOnS");
        }
        self.tg2ds
    }

    /// OCCT TangentOnRst() (BRepBlend_SurfRstConstRad.cxx L500-508).
    pub fn tangent_on_rst(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstConstRad::TangentOnRst");
        }
        self.tgrst
    }

    /// OCCT Tangent2dOnRst() (BRepBlend_SurfRstConstRad.cxx L512-519).
    pub fn tangent_2d_on_rst(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstConstRad::Tangent2dOnRst");
        }
        self.tg2drst
    }

    /// OCCT Decroch(Sol, NS, TgS) (BRepBlend_SurfRstConstRad.cxx L523-575).
    // The OCCT body assigns NSInPlane at L531 and reverses it at L540
    // without any further read (a dead store in the original); the rcad
    // port keeps the statements for parity.
    #[allow(unused_assignments)]
    pub fn decroch(&self, sol: &[f64], ns_out: &mut DVec3, tg_s: &mut DVec3) -> bool {
        // OCCT L530-531: surf->D1(Sol(1), Sol(2), bid, d1u, d1v);
        // NS = NSInPlane = d1u.Crossed(d1v);
        let (bid, d1u, d1v) = self.surf.derivatives(sol[0], sol[1]);
        let ns = d1u.cross(d1v);
        *ns_out = ns;
        let mut ns_in_plane = ns;

        // OCCT L533-535.
        let norm = self.nplan.cross(ns).length();
        let unsurnorm = 1.0 / norm;
        ns_in_plane = (self.nplan.dot(ns) * unsurnorm) * self.nplan + (-unsurnorm) * ns;

        // OCCT L537-541: Center.SetXYZ(bid.XYZ() + ray * NSInPlane.XYZ());
        let center = bid + self.ray * ns_in_plane;
        if self.choix > 2 {
            ns_in_plane = -ns_in_plane;
        }
        // OCCT L542-546: TgS = nplan.Crossed(gp_Vec(Center, bid));
        *tg_s = self.nplan.cross(bid - center);
        if self.choix % 2 == 1 {
            *tg_s = -*tg_s;
        }

        // OCCT L547-549: rstref->Value(Sol(3)).Coord(u, v);
        // surfref->D1(u, v, bid, d1u, d1v);
        let (u, v) = {
            let p = self.rstref.expect("rstref").point_at(sol[2]);
            (p.x, p.y)
        };
        let (bid, d1u, d1v) = self.surfref.expect("surfref").derivatives(u, v);

        // OCCT L550-553.
        let n_rst = d1u.cross(d1v);
        let norm = self.nplan.cross(n_rst).length();
        let unsurnorm = 1.0 / norm;
        let mut n_rst_in_plane =
            (self.nplan.dot(n_rst) * unsurnorm) * self.nplan + (-unsurnorm) * n_rst;
        // OCCT L554-558: gp_Vec centptrst(Center, bid);
        let centptrst = bid - center;
        if centptrst.dot(n_rst_in_plane) < 0.0 {
            n_rst_in_plane = -n_rst_in_plane;
        }
        // OCCT L559-563.
        let mut tg_rst = self.nplan.cross(centptrst);
        if self.choix % 2 == 1 {
            tg_rst = -tg_rst;
        }

        // OCCT L565-574.
        let mut nt = n_rst_in_plane.length();
        nt *= tg_rst.length();
        if nt.abs() < 1.0e-7 {
            return false; // Singularity or Incoherence.
        }
        let mut dot = n_rst_in_plane.dot(tg_rst);
        dot /= nt;

        dot < 1.0e-10
    }

    /// OCCT Set(Radius, Choix) (BRepBlend_SurfRstConstRad.cxx L579-596).
    pub fn set(&mut self, radius: f64, choix: i32) {
        self.choix = choix;
        match self.choix {
            1 | 2 => {
                self.ray = -radius.abs();
            }
            3 | 4 => {
                self.ray = radius.abs();
            }
            _ => {
                self.ray = -radius.abs();
            }
        }
    }

    /// OCCT Set(BlendFunc_SectionShape) (BRepBlend_SurfRstConstRad.cxx
    /// L600-603).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_s_shape = type_section;
    }

    /// OCCT Section(Param, U, V, W, Pdeb, Pfin, C)
    /// (BRepBlend_SurfRstConstRad.cxx L607-658).
    #[allow(clippy::too_many_arguments)]
    pub fn section(
        &mut self,
        param: f64,
        u: f64,
        v: f64,
        w: f64,
        pdeb: &mut f64,
        pfin: &mut f64,
        c: &mut Circle3,
    ) {
        // OCCT L620-621: tguide->D1(Param, ptgui, d1gui); np = d1gui.Normalized();
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        let mut np = self.d1gui.normalize_or_zero();

        // OCCT L623-624: surf->D1(U, V, pts, d1u1, d1v1); ptrst = cons.Value(W);
        let (p, d1u1, d1v1) = self.surf.derivatives(u, v);
        self.pts = p;
        self.ptrst = self.cons_value(w);

        // OCCT L626-630.
        let mut ns = d1u1.cross(d1v1);
        let norm = self.nplan.cross(ns).length();
        ns = (self.nplan.dot(ns) / norm) * self.nplan + (-1.0 / norm) * ns;
        let center = self.pts + self.ray * ns;
        c.radius = self.ray.abs();

        // OCCT L633-641.
        if self.ray > 0.0 {
            ns = -ns;
        }
        if self.choix % 2 != 0 {
            np = -np;
        }

        // OCCT L643: C.SetPosition(gp_Ax2(Center, np, ns))
        // (gp_Ax2(P, N, Vx): the Y direction is N ^ Vx).
        c.center = center;
        c.normal = np;
        c.x_dir = ns;
        c.y_dir = np.cross(ns);
        // OCCT L644-645: Pdeb = 0; Pfin = ElCLib::Parameter(C, ptrst);
        *pdeb = 0.0;
        *pfin = elclib_circle_parameter(c, self.ptrst);

        // Test negative and almost null angles : Special case
        // OCCT L648-653.
        if *pfin > 1.5 * std::f64::consts::PI {
            np = -np;
            c.normal = np;
            c.y_dir = np.cross(ns);
            *pfin = elclib_circle_parameter(c, self.ptrst);
        }
        // OCCT L654-657.
        if *pfin < p_confusion() {
            *pfin += p_confusion();
        }
    }

    /// OCCT IsRational() (BRepBlend_SurfRstConstRad.cxx L662-665).
    pub fn is_rational(&self) -> bool {
        self.my_s_shape == BlendFuncSectionShape::Rational
            || self.my_s_shape == BlendFuncSectionShape::QuasiAngular
    }

    /// OCCT GetSectionSize() (BRepBlend_SurfRstConstRad.cxx L669-672).
    pub fn get_section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// OCCT GetMinimalWeight(Weights) (BRepBlend_SurfRstConstRad.cxx
    /// L676-680).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        blend_func_get_minimal_weights(self.my_s_shape, self.my_t_conv, self.minang, self.maxang, weigths);
        // It is supposed that it does not depend on the Radius!
    }

    /// OCCT NbIntervals(S) (BRepBlend_SurfRstConstRad.cxx L684-687).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.guide.nb_intervals(blend_func_next_shape(s))
    }

    /// OCCT Intervals(T, S) (BRepBlend_SurfRstConstRad.cxx L691-695).
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let mut intervals = Vec::new();
        self.guide.intervals(&mut intervals, blend_func_next_shape(s));
        for (dst, src) in t.iter_mut().zip(intervals) {
            *dst = src;
        }
    }

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
    /// (BRepBlend_SurfRstConstRad.cxx L699-703).
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
    /// (BRepBlend_SurfRstConstRad.cxx L710-723) — tolerances used for
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

    /// OCCT Knots(TKnots) (BRepBlend_SurfRstConstRad.cxx L727-730).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        geomfill_knots(self.my_t_conv, tknots);
    }

    /// OCCT Mults(TMults) (BRepBlend_SurfRstConstRad.cxx L734-737).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        geomfill_mults(self.my_t_conv, tmults);
    }

    /// OCCT Section(P, Poles, Poles2d, Weights)
    /// (BRepBlend_SurfRstConstRad.cxx L741-798).
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

        // OCCT L756-757: tguide->D1(prm, ptgui, d1gui); nplan = d1gui.Normalized();
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        self.nplan = self.d1gui.normalize_or_zero();

        // OCCT L759-761.
        let (u1, v1) = p.parameters_on_s();
        let w = p.parameter_on_c(); // jlr : point on curve not on surface
        let pt2d = self.rst.point_at(w);

        // OCCT L763-765.
        let (ptsv, d1u1, d1v1) = self.surf.derivatives(u1, v1);
        self.pts = ptsv;
        self.ptrst = self.cons_value(w);
        self.distmin = self.distmin.min(self.pts.distance(self.ptrst));

        // OCCT L767-768.
        poles_2d[0] = DVec2::new(u1, v1);
        poles_2d[poles_2d.len() - 1] = DVec2::new(pt2d.x, pt2d.y);

        // Linear Case
        // OCCT L770-778.
        if self.my_s_shape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptrst;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            return;
        }

        // OCCT L780-787.
        let mut ns = d1u1.cross(d1v1);
        let norm = self.nplan.cross(ns).length();
        ns = (self.nplan.dot(ns) / norm) * self.nplan + (-1.0 / norm) * ns;
        let center = self.pts + self.ray * ns;
        // OCCT L787: ns2 = gp_Vec(Center, ptrst).Normalized();
        let ns2 = (self.ptrst - center).normalize_or_zero();

        // OCCT L788-795.
        if self.ray > 0.0 {
            ns = -ns;
        }
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
        }

        // OCCT L797: GeomFill::GetCircle(myTConv, ns, ns2, nplan, pts, ptrst,
        // std::abs(ray), Center, Poles, Weights).
        get_circle(
            tconv(self.my_t_conv),
            ns,
            ns2,
            self.nplan,
            self.pts,
            self.ptrst,
            self.ray.abs(),
            center,
            poles,
            weigths,
        );
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)
    /// (BRepBlend_SurfRstConstRad.cxx L802-1015) — used for the first and
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

        // OCCT L828-831.
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        self.d2gui = self.tguide().derivative2_at(prm);
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        // OCCT L831: dnplan.SetLinearForm(1. / normtg, d2gui,
        // -1. / normtg * (nplan.Dot(d2gui)), nplan);
        let mut dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        // OCCT L833-835.
        let (u1, v1) = p.parameters_on_s();
        let mut sol = [0.0f64; 3];
        sol[0] = u1;
        sol[1] = v1;
        self.prmrst = p.parameter_on_c();
        sol[2] = self.prmrst;
        self.pt2drst = self.rst.point_at(self.prmrst);

        // OCCT L837: Values(sol, valsol, gradsol);
        let mut valsol = [0.0f64; 3];
        let mut gradsol = vec![vec![0.0f64; 3]; 3];
        self.values(&sol, &mut valsol, &mut gradsol);

        // OCCT L839-840.
        let (ptsv, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.derivatives2(sol[0], sol[1]);
        self.pts = ptsv;
        let (ptrst_v, d1) = self.cons_d1(sol[2]);
        self.ptrst = ptrst_v;

        // OCCT L842-846.
        let mut secmember = [0.0f64; 3];
        let temp = self.pts - self.ptgui;
        secmember[0] = self.normtg - dnplan.dot(temp);

        let temp = self.ptrst - self.ptgui;
        secmember[1] = self.normtg - dnplan.dot(temp);

        // OCCT L848-858.
        let mut ns = d1u1.cross(d1v1);
        let ncrossns = self.nplan.cross(ns);
        let ndotns = self.nplan.dot(ns);
        let mut norm = ncrossns.length();
        if norm < EPS {
            norm = 1.0; // Not enough, but it is not necessary to stop
        }

        // Derivative of n1 corresponding to w
        // OCCT L862-868.
        let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
        let mut dnw = ((dnplan.dot(ns) - grosterme * ndotns) / norm) * self.nplan
            + (ndotns / norm) * dnplan
            + (grosterme / norm) * ns;

        // OCCT L870-873.
        let temp = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;
        let resul = self.ray * temp + (self.pts - self.ptrst);
        secmember[2] = dnw.dot(resul);
        secmember[2] = -2.0 * self.ray * secmember[2];

        // OCCT L875-896: math_Gauss Resol(gradsol, 1.e-9) and the SVD
        // fallback.
        let mut a = MatD::new(3, 3);
        for r in 1..=3 {
            for c in 1..=3 {
                a.set(r, c, gradsol[r - 1][c - 1]);
            }
        }
        let resol = MathGauss::new(&a);
        let istgt: bool;
        if resol.is_done() {
            istgt = false;
            // OCCT L880: Resol.Solve(secmember);
            let mut x = VecD::new(3);
            for i in 1..=3 {
                x.set(i, secmember[i - 1]);
            }
            resol.solve(&mut x);
            for i in 1..=3 {
                secmember[i - 1] = x.get(i);
            }
        } else {
            // OCCT L884-891: math_SVD SingRS(gradsol);
            // if (SingRS.IsDone()) { math_Vector DEDT(1,3); DEDT = secmember;
            // SingRS.Solve(DEDT, secmember, 1.e-6); istgt = false; }
            // else { istgt = true; }
            let mut sing_rs = MathSvd::new(&a);
            if sing_rs.is_done() {
                let mut dedt = VecD::new(3);
                let mut sol = VecD::new(3);
                for i in 1..=3 {
                    dedt.set(i, secmember[i - 1]);
                    sol.set(i, secmember[i - 1]);
                }
                sing_rs.solve(&dedt, &mut sol, 1.0e-6);
                for i in 1..=3 {
                    secmember[i - 1] = sol.get(i);
                }
                istgt = false;
            } else {
                istgt = true;
            }
        }

        // OCCT L898-936.
        let ns2: DVec3;
        let mut dn2w = DVec3::ZERO;
        if !istgt {
            // OCCT L900-901.
            self.tgs = secmember[0] * d1u1 + secmember[1] * d1v1;
            self.tgrst = secmember[2] * d1;

            // Derivative of n1 corresponding to u1
            // OCCT L904-911.
            let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
            let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
            let resulu = (-(grosterme * ndotns - self.nplan.dot(temp)) / norm) * self.nplan
                + (grosterme / norm) * ns
                + (-1.0 / norm) * temp;

            // Derivative of n1 corresponding to v1
            // OCCT L914-921.
            let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
            let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
            let resulv = (-(grosterme * ndotns - self.nplan.dot(temp)) / norm) * self.nplan
                + (grosterme / norm) * ns
                + (-1.0 / norm) * temp;

            // OCCT L923: dnw.SetLinearForm(secmember(1), resulu, secmember(2), resulv, dnw);
            dnw = secmember[0] * resulu + secmember[1] * resulv + dnw;
            // OCCT L924.
            ns = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;

            // OCCT L926-930.
            dn2w = self.ray * dnw + (-1.0) * self.tgrst + self.tgs;
            norm = resul.length();
            dn2w /= norm;
            ns2 = -resul.normalize_or_zero();
            dn2w = ns2.dot(dn2w) * ns2 + (-1.0) * dn2w;
        } else {
            // OCCT L934-935.
            ns = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;
            ns2 = -resul.normalize_or_zero();
        }

        // Tops 2D
        // OCCT L940-949.
        poles_2d[0] = DVec2::new(sol[0], sol[1]);
        poles_2d[poles_2d.len() - 1] = DVec2::new(self.pt2drst.x, self.pt2drst.y);
        if !istgt {
            d_poles_2d[0] = DVec2::new(secmember[0], secmember[1]);
            let (_, d1urst, d1vrst) = self.surfrst.derivatives(self.pt2drst.x, self.pt2drst.y);
            let mut a2d = 0.0f64;
            let mut b2d = 0.0f64;
            t3dto2d(&mut a2d, &mut b2d, self.tgrst, d1urst, d1vrst);
            d_poles_2d[d_poles_2d.len() - 1] = DVec2::new(a2d, b2d);
        }

        // Linear Case
        // OCCT L951-966.
        if self.my_s_shape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptrst;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            if !istgt {
                d_poles[low] = self.tgs;
                d_poles[upp] = self.tgrst;
                d_weigths[low] = 0.0;
                d_weigths[upp] = 0.0;
            }
            return !istgt;
        }

        // Case of the circle
        // OCCT L969-972.
        let center = self.pts + self.ray * ns;
        let tgct: DVec3;
        if !istgt {
            tgct = self.tgs + self.ray * dnw;
        } else {
            tgct = DVec3::ZERO; // the OCCT tgct is only read on the !istgt path
        }

        // OCCT L975-987.
        if self.ray > 0.0 {
            ns = -ns;
            if !istgt {
                dnw = -dnw;
            }
        }
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
            dnplan = -dnplan;
        }

        // OCCT L988-1014.
        if !istgt {
            get_circle_d1(
                tconv(self.my_t_conv),
                ns,
                ns2,
                dnw,
                dn2w,
                self.nplan,
                dnplan,
                self.pts,
                self.ptrst,
                self.tgs,
                self.tgrst,
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
                ns,
                ns2,
                self.nplan,
                self.pts,
                self.ptrst,
                self.ray.abs(),
                center,
                poles,
                weigths,
            );
            false
        }
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weigths, DWeigths, D2Weigths) (BRepBlend_SurfRstConstRad.cxx
    /// L1019-1031) — returns false.
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
    /// (BRepBlend_SurfRstConstRad.cxx L1033-1048).
    pub fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic_2d == 1 {
            *tol_u = self.surf.u_resolution(tol);
            *tol_v = self.surf.v_resolution(tol);
        } else {
            *tol_u = self.surfrst.u_resolution(tol);
            *tol_v = self.surfrst.v_resolution(tol);
        }
    }
}

impl<'a> FunctionSetWithDerivatives for BlendSurfRstConstRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendSurfRstConstRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendSurfRstConstRad::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendSurfRstConstRad::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfRstConstRad::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfRstConstRad::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendSurfRstConstRad<'a> {
    fn set_param(&mut self, param: f64) {
        BlendSurfRstConstRad::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendSurfRstConstRad::set_interval(self, first, last)
    }

    fn pnt1(&self) -> DVec3 {
        BlendSurfRstFunction::pnt1(self)
    }

    fn pnt2(&self) -> DVec3 {
        BlendSurfRstFunction::pnt2(self)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendSurfRstConstRad::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendSurfRstConstRad::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendSurfRstConstRad::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendSurfRstConstRad::get_minimal_distance(self)
    }

    fn is_rational(&self) -> bool {
        BlendSurfRstConstRad::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendSurfRstConstRad::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendSurfRstConstRad::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendSurfRstConstRad::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendSurfRstConstRad::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendSurfRstConstRad::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendSurfRstConstRad::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendSurfRstConstRad::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendSurfRstConstRad::mults(self, tmults)
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
        BlendSurfRstConstRad::section_d1(
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
        BlendSurfRstConstRad::section_simple(self, p, poles, poles_2d, weigths)
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
        BlendSurfRstConstRad::section_d2(
            self, p, poles, d_poles, d2_poles, poles_2d, d_poles_2d, d2_poles_2d, weigths,
            d_weigths, d2_weigths,
        )
    }

    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        BlendSurfRstConstRad::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendSurfRstFunction for BlendSurfRstConstRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendSurfRstConstRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendSurfRstConstRad::nb_equations(self)
    }

    fn point_on_s(&self) -> DVec3 {
        BlendSurfRstConstRad::point_on_s(self)
    }

    fn point_on_rst(&self) -> DVec3 {
        BlendSurfRstConstRad::point_on_rst(self)
    }

    fn pnt2d_on_s(&self) -> DVec2 {
        BlendSurfRstConstRad::pnt2d_on_s(self)
    }

    fn pnt2d_on_rst(&self) -> DVec2 {
        BlendSurfRstConstRad::pnt2d_on_rst(self)
    }

    fn parameter_on_rst(&self) -> f64 {
        BlendSurfRstConstRad::parameter_on_rst(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendSurfRstConstRad::is_tangency_point(self)
    }

    fn tangent_on_s(&self) -> DVec3 {
        BlendSurfRstConstRad::tangent_on_s(self)
    }

    fn tangent_2d_on_s(&self) -> DVec2 {
        BlendSurfRstConstRad::tangent_2d_on_s(self)
    }

    fn tangent_on_rst(&self) -> DVec3 {
        BlendSurfRstConstRad::tangent_on_rst(self)
    }

    fn tangent_2d_on_rst(&self) -> DVec2 {
        BlendSurfRstConstRad::tangent_2d_on_rst(self)
    }

    fn decroch(&self, sol: &[f64], ns: &mut DVec3, tg_s: &mut DVec3) -> bool {
        BlendSurfRstConstRad::decroch(self, sol, ns, tg_s)
    }
}
