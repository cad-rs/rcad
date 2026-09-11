//! OCCT BRepBlend_SurfRstEvolRad (TKFillet/BRepBlend) — 1:1 port of
//! BRepBlend_SurfRstEvolRad.hxx (L40-244) + BRepBlend_SurfRstEvolRad.cxx
//! (whole file L37-1166).  Variable evolutive section : variable radius
//! rolling ball on a surface and a pcurve on another surface; the
//! SimulSurf / PerformSurf face/rst and rst/face variable-radius arms
//! instantiate it as `func` (ChFi3d_FilBuilder.cxx L1146 / L1998).
//!
//! Architecture mappings (mirroring [`super::brep_blend_surf_rst_const_rad`]):
//! `class BRepBlend_SurfRstEvolRad : public Blend_SurfRstFunction` is
//! expressed by implementing the [`BlendSurfRstFunction`] and
//! [`BlendAppFunction`] traits over the `math_FunctionSetWithDerivatives`
//! base; `occ::handle<Adaptor3d_Curve> tguide` (aliasing `guide` until
//! Set(First, Last) trims it) maps to an owned trimmed [`Curve3`] copy plus
//! an accessor, and likewise `tevol` (aliasing `fevol` until Set(First,
//! Last) trims it); `math_Vector` / `math_Matrix` map to `[f64; 3]` /
//! `Vec<Vec<f64>>` (OCCT X(i) -> x[i - 1], D(i, j) -> d[i - 1][j - 1]);
//! `gp_Circ` maps to the kernel [`Circle3`].
//!
//! The OCCT `Adaptor3d_CurveOnSurface cons` member (hxx) has no rcad
//! adaptor equivalent; the rcad port stores the (rst, surfrst) pair the
//! adaptor wraps and transcribes the consumed operations from
//! Adaptor3d_CurveOnSurface.cxx (Value/EvalD0, D1/EvalD1 generic branch,
//! FirstParameter/LastParameter L977-987, Resolution L1364-1370) in the
//! `cons_*` helpers below.
//!
//! Pending kernel dependency (marked GAP, plan 0.6): math_SVD (the
//! second-chance solver of IsSolution / Section-d1 — the OCCT !IsDone()
//! route is preserved).

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::{is_infinite_value, p_confusion};
use rcad_kernel::geom::{Circle3, Curve2d, Curve2dEval as _, Curve3, CurveEval as _, Surface3, SurfaceEval as _};
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::math_svd::MathSvd;
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::{MatD, VecD};

use crate::geomalgo::geomfill::geom_fill::{get_circle, get_circle_d1};
use crate::geomalgo::law::law_function::LawFunctionHandle;

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
use super::brep_blend_surf_rst_function::BlendSurfRstFunction;

/// OCCT Eps constant (BRepBlend_SurfRstEvolRad.cxx L37).
const EPS: f64 = 1.0e-15;

/// OCCT static t3dto2d (BRepBlend_SurfRstEvolRad.cxx L39-49).
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

/// OCCT BRepBlend_SurfRstEvolRad.
pub struct BlendSurfRstEvolRad<'a> {
    // OCCT BRepBlend_SurfRstEvolRad.hxx fields.
    pub(crate) surf: &'a Surface3,
    pub(crate) surfrst: &'a Surface3,
    pub(crate) rst: &'a Curve2d,
    // OCCT: Adaptor3d_CurveOnSurface cons(Rst, SurfRst) — carried as the
    // wrapped pair (see the module header architecture note).
    pub(crate) guide: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> tguide — a handle aliasing guide
    // until Set(First, Last) trims it (cxx L332).  The rcad port owns the
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
    pub(crate) dray: f64,
    pub(crate) choix: i32,
    pub(crate) sg1: f64,
    pub(crate) ptgui: DVec3,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) normtg: f64,
    pub(crate) the_d: f64,
    // OCCT: occ::handle<Adaptor3d_Surface> surfref — null until Set(SurfRef,
    // RstRef) (cxx L300-305).
    pub(crate) surfref: Option<&'a Surface3>,
    pub(crate) rstref: Option<&'a Curve2d>,
    pub(crate) maxang: f64,
    pub(crate) minang: f64,
    pub(crate) distmin: f64,
    pub(crate) my_s_shape: BlendFuncSectionShape,
    pub(crate) my_t_conv: ConvertParameterisationType,
    // OCCT: occ::handle<Law_Function> tevol / fevol.
    pub(crate) tevol: Option<LawFunctionHandle>,
    pub(crate) fevol: LawFunctionHandle,
}

impl<'a> BlendSurfRstEvolRad<'a> {
    /// Architecture mapping: OCCT `tguide` is a handle aliasing `guide`
    /// until Set(First, Last) replaces it with a trimmed copy (cxx L332).
    #[inline]
    pub(crate) fn tguide(&self) -> &Curve3 {
        self.tguide.as_ref().unwrap_or(self.guide)
    }

    /// Architecture mapping: OCCT `tevol` is a handle aliasing `fevol`
    /// until Set(First, Last) replaces it with a trimmed law (cxx L333).
    #[inline]
    pub(crate) fn tevol(&self) -> &LawFunctionHandle {
        self.tevol.as_ref().unwrap_or(&self.fevol)
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

    /// OCCT BRepBlend_SurfRstEvolRad(Surf, SurfRst, Rst, CGuide, Evol)
    /// (BRepBlend_SurfRstEvolRad.cxx L114-134).
    pub fn new(
        surf: &'a Surface3,
        surfrst: &'a Surface3,
        rst: &'a Curve2d,
        cguide: &'a Curve3,
        evol: LawFunctionHandle,
    ) -> Self {
        BlendSurfRstEvolRad {
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
            dray: 0.0,
            choix: 0,
            sg1: 0.0,
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
            tevol: None, // OCCT: tevol(Evol) — alias, see tevol().
            fevol: evol, // OCCT L132-133: tevol = Evol; fevol = Evol.
        }
    }

    /// OCCT NbVariables() (BRepBlend_SurfRstEvolRad.cxx L138-141) —
    /// returns 3.
    pub fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations() (BRepBlend_SurfRstEvolRad.cxx L145-148) —
    /// returns 3.
    pub fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value(X, F) (BRepBlend_SurfRstEvolRad.cxx L152-170).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // OCCT L157: surf->D1(X(1), X(2), pts, d1u1, d1v1).
        let (p, d1u1, d1v1) = self.surf.derivatives(x[0], x[1]);
        self.pts = p;
        // OCCT L158: ptrst = cons.Value(X(3)).
        self.ptrst = self.cons_value(x[2]);

        // OCCT L160-162.
        f[0] = self.nplan.dot(self.pts) + self.the_d;
        f[1] = self.nplan.dot(self.ptrst) + self.the_d;

        // OCCT L164-166.
        let ns = d1u1.cross(d1v1);
        let norm = self.nplan.cross(ns).length();
        let ns = (self.nplan.dot(ns) / norm) * self.nplan + (-1.0 / norm) * ns;
        // OCCT L167-168.
        let vref = self.ray * ns + (self.pts - self.ptrst);
        f[2] = vref.length_squared() - self.ray * self.ray;
        true
    }

    /// OCCT Derivatives(X, D) (BRepBlend_SurfRstEvolRad.cxx L174-230).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L181: surf->D2(X(1), X(2), pts, d1u1, d1v1, d2u1, d2v1, d2uv1).
        let (p, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.derivatives2(x[0], x[1]);
        self.pts = p;
        // OCCT L182: cons.D1(X(3), ptrst, d1).
        let (ptrst_v, d1) = self.cons_d1(x[2]);
        self.ptrst = ptrst_v;

        // OCCT L184-190.
        d[0][0] = self.nplan.dot(d1u1);
        d[0][1] = self.nplan.dot(d1v1);
        d[0][2] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(d1);

        // OCCT L192-199.
        let ns = d1u1.cross(d1v1);
        let ncrossns = self.nplan.cross(ns);
        let norm = ncrossns.length();
        let ndotns = self.nplan.dot(ns);

        // OCCT L197-199: vref.SetLinearForm(ndotns, nplan, -1., ns);
        // vref.Divide(norm); vref.SetLinearForm(ray, vref, gp_Vec(ptrst, pts)).
        let vref = (ndotns * self.nplan + (-1.0) * ns) / norm;
        let vref = self.ray * vref + (self.pts - self.ptrst);

        // Derivative corresponding to u1
        // OCCT L202-210.
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1u1;

        // OCCT L212.
        d[2][0] = 2.0 * resul.dot(vref);

        // Derivative corresponding to v1
        // OCCT L215-223.
        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1v1;

        // OCCT L225.
        d[2][1] = 2.0 * resul.dot(vref);

        // OCCT L227.
        d[2][2] = -2.0 * d1.dot(vref);

        true
    }

    /// OCCT Values(X, F, D) (BRepBlend_SurfRstEvolRad.cxx L234-296).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // OCCT L242: surf->D2(X(1), X(2), pts, d1u1, d1v1, d2u1, d2v1, d2uv1).
        let (p, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.derivatives2(x[0], x[1]);
        self.pts = p;
        // OCCT L243: cons.D1(X(3), ptrst, d1).
        let (ptrst_v, d1) = self.cons_d1(x[2]);
        self.ptrst = ptrst_v;

        // OCCT L245-246.
        f[0] = self.nplan.dot(self.pts) + self.the_d;
        f[1] = self.nplan.dot(self.ptrst) + self.the_d;

        // OCCT L248-254.
        d[0][0] = self.nplan.dot(d1u1);
        d[0][1] = self.nplan.dot(d1v1);
        d[0][2] = 0.0;

        d[1][0] = 0.0;
        d[1][1] = 0.0;
        d[1][2] = self.nplan.dot(d1);

        // OCCT L256-263.
        let ns = d1u1.cross(d1v1);
        let ncrossns = self.nplan.cross(ns);
        let norm = ncrossns.length();
        let ndotns = self.nplan.dot(ns);

        // OCCT L261-263.
        let vref = (ndotns * self.nplan + (-1.0) * ns) / norm;
        let vref = self.ray * vref + (self.pts - self.ptrst);

        // OCCT L265.
        f[2] = vref.length_squared() - self.ray * self.ray;

        // Derivative corresponding to u1
        // OCCT L268-276.
        let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1u1;

        // OCCT L278.
        d[2][0] = 2.0 * resul.dot(vref);

        // Derivative corresponding to v1
        // OCCT L281-289.
        let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
        let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
        let resul = (-self.ray / norm * (grosterme * ndotns - self.nplan.dot(temp))) * self.nplan
            + (self.ray * grosterme / norm) * ns
            + (-self.ray / norm) * temp
            + d1v1;

        // OCCT L291.
        d[2][1] = 2.0 * resul.dot(vref);

        // OCCT L293.
        d[2][2] = -2.0 * d1.dot(vref);

        true
    }

    /// OCCT Set(SurfRef, RstRef) (BRepBlend_SurfRstEvolRad.cxx L300-305).
    pub fn set_ref(&mut self, surf_ref: &'a Surface3, rst_ref: &'a Curve2d) {
        self.surfref = Some(surf_ref);
        self.rstref = Some(rst_ref);
    }

    /// OCCT Set(Param) (BRepBlend_SurfRstEvolRad.cxx L309-323).
    pub fn set_param(&mut self, param: f64) {
        self.d1gui = DVec3::ZERO;
        self.nplan = DVec3::ZERO;
        // OCCT L313: tguide->D2(Param, ptgui, d1gui, d2gui).
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        self.d2gui = self.tguide().derivative2_at(param);
        // OCCT L314-319: theD is carried through the two-step gp_XYZ form.
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        let the_d = self.nplan.dot(self.ptgui);
        self.the_d = the_d * (-1.0);
        // OCCT L320-322: tevol->D1(Param, ray, dray); ray = sg1 * ray;
        // dray = sg1 * dray.
        let (mut ray, mut dray) = (0.0f64, 0.0f64);
        self.tevol().borrow_mut().d1(param, &mut ray, &mut dray);
        self.ray = self.sg1 * ray;
        self.dray = self.sg1 * dray;
    }

    /// OCCT Set(First, Last) (BRepBlend_SurfRstEvolRad.cxx L330-334) —
    /// segments the curve and the law in their useful part.
    pub fn set_interval(&mut self, first: f64, last: f64) {
        self.tguide = Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(self.guide.clone()),
            first,
            last,
        }));
        // OCCT L333: tevol = fevol->Trim(First, Last, 1.e-12).
        self.tevol = Some(self.fevol.borrow().trim(first, last, 1.0e-12));
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BRepBlend_SurfRstEvolRad.cxx
    /// L338-343).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L340-341.
        tolerance[0] = self.surf.u_resolution(tol);
        tolerance[1] = self.surf.v_resolution(tol);
        // OCCT L342: Tolerance(3) = cons.Resolution(Tol).
        tolerance[2] = self.cons_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BRepBlend_SurfRstEvolRad.cxx
    /// L347-368).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        // OCCT L349-354 (cons.FirstParameter / LastParameter = the pcurve's).
        inf_bound[0] = self.surf.default_domain()[0]; // FirstUParameter
        inf_bound[1] = self.surf.default_domain()[2]; // FirstVParameter
        inf_bound[2] = self.rst.default_domain()[0]; // cons.FirstParameter
        sup_bound[0] = self.surf.default_domain()[1]; // LastUParameter
        sup_bound[1] = self.surf.default_domain()[3]; // LastVParameter
        sup_bound[2] = self.rst.default_domain()[1]; // cons.LastParameter

        // OCCT L356-361.
        if !is_infinite_value(inf_bound[0]) && !is_infinite_value(sup_bound[0]) {
            let range = sup_bound[0] - inf_bound[0];
            inf_bound[0] -= range;
            sup_bound[0] += range;
        }
        // OCCT L362-367.
        if !is_infinite_value(inf_bound[1]) && !is_infinite_value(sup_bound[1]) {
            let range = sup_bound[1] - inf_bound[1];
            inf_bound[1] -= range;
            sup_bound[1] += range;
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BRepBlend_SurfRstEvolRad.cxx L372-491).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let mut valsol = [0.0f64; 3];
        let mut secmember = [0.0f64; 3];
        let mut gradsol = vec![vec![0.0f64; 3]; 3];

        // OCCT L382: Values(Sol, valsol, gradsol).
        self.values(sol, &mut valsol, &mut gradsol);
        // OCCT L383-385.
        if valsol[0].abs() <= tol
            && valsol[1].abs() <= tol
            && valsol[2].abs() <= 2.0 * tol * self.ray.abs()
        {
            // Calculation of tangents

            // OCCT L389-391.
            self.pt2ds = DVec2::new(sol[0], sol[1]);
            self.prmrst = sol[2];
            self.pt2drst = self.rst.point_at(self.prmrst);
            // OCCT L392-393.
            let (p, d1u1, d1v1) = self.surf.derivatives(sol[0], sol[1]);
            self.pts = p;
            let (ptrst_v, d1) = self.cons_d1(sol[2]);
            self.ptrst = ptrst_v;
            // OCCT L394: dnplan.SetLinearForm(1. / normtg, d2gui,
            // -1. / normtg * (nplan.Dot(d2gui)), nplan).
            let dnplan = (1.0 / self.normtg) * self.d2gui
                + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

            // OCCT L396-397.
            let temp = self.pts - self.ptgui;
            secmember[0] = self.normtg - dnplan.dot(temp);

            // OCCT L399-400.
            let temp = self.ptrst - self.ptgui;
            secmember[1] = self.normtg - dnplan.dot(temp);

            // OCCT L402-414.
            let mut ns = d1u1.cross(d1v1);
            let ncrossns = self.nplan.cross(ns);
            let ndotns = self.nplan.dot(ns);
            let norm = ncrossns.length();

            let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
            let dnw = ((dnplan.dot(ns) - grosterme * ndotns) / norm) * self.nplan
                + (ndotns / norm) * dnplan
                + (grosterme / norm) * ns;

            // OCCT L416-417.
            ns = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;
            // OCCT L417: resul.SetLinearForm(ray, ns, gp_Vec(ptrst, pts))
            // (gp_Vec(ptrst, pts) = pts - ptrst).
            let resul = self.ray * ns + (self.pts - self.ptrst);
            // OCCT L419.
            secmember[2] =
                -2.0 * self.ray * dnw.dot(resul) - 2.0 * self.dray * ns.dot(resul) + 2.0 * self.ray * self.dray;

            // OCCT L420-440: math_Gauss Resol(gradsol) and the SVD fallback.
            let mut a = MatD::new(3, 3);
            for r in 1..=3 {
                for c in 1..=3 {
                    a.set(r, c, gradsol[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);
            if resol.is_done() {
                // OCCT L423: Resol.Solve(secmember).
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
                // OCCT L428-439: math_SVD SingRS(gradsol);
                // if (SingRS.IsDone()) { math_Vector DEDT = secmember;
                // SingRS.Solve(DEDT, secmember, 1.e-6); istangent = false; }
                // else { istangent = true; }
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

            // OCCT L442-456 (the trailing istangent = false / = true
            // re-assignments are dead in the original; the rcad port keeps
            // the state as set by the solver branches).
            if !self.istangent {
                // OCCT L444-450.
                self.tgs = secmember[0] * d1u1 + secmember[1] * d1v1;
                self.tgrst = secmember[2] * d1;
                self.tg2ds = DVec2::new(secmember[0], secmember[1]);
                let (_, d1urst, d1vrst) = self.surfrst.derivatives(self.pt2drst.x, self.pt2drst.y);
                let mut a2d = 0.0f64;
                let mut b2d = 0.0f64;
                t3dto2d(&mut a2d, &mut b2d, self.tgrst, d1urst, d1vrst);
                self.tg2drst = DVec2::new(a2d, b2d);
            }

            // update of maxang
            // OCCT L458-461.
            if self.ray > 0.0 {
                ns = -ns;
            }
            // OCCT L462: ns2 = -resul.Normalized().
            let ns2 = -resul.normalize_or_zero();

            // OCCT L464-469.
            let cosa = ns.dot(ns2);
            let mut sina = self.nplan.dot(ns.cross(ns2));
            if self.choix % 2 != 0 {
                sina = -sina; // nplan is changed into -nplan
            }

            // OCCT L471-475.
            let mut angle = cosa.acos();
            if sina < 0.0 {
                angle = 2.0 * std::f64::consts::PI - angle;
            }

            // OCCT L477-484.
            if angle > self.maxang {
                self.maxang = angle;
            }
            if angle < self.minang {
                self.minang = angle;
            }
            // OCCT L485.
            self.distmin = self.distmin.min(self.pts.distance(self.ptrst));

            return true;
        }
        // OCCT L489-490.
        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BRepBlend_SurfRstEvolRad.cxx L495-498).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT PointOnS() (BRepBlend_SurfRstEvolRad.cxx L502-505).
    pub fn point_on_s(&self) -> DVec3 {
        self.pts
    }

    /// OCCT PointOnRst() (BRepBlend_SurfRstEvolRad.cxx L509-512).
    pub fn point_on_rst(&self) -> DVec3 {
        self.ptrst
    }

    /// OCCT Pnt2dOnS() (BRepBlend_SurfRstEvolRad.cxx L516-519).
    pub fn pnt2d_on_s(&self) -> DVec2 {
        self.pt2ds
    }

    /// OCCT Pnt2dOnRst() (BRepBlend_SurfRstEvolRad.cxx L523-526).
    pub fn pnt2d_on_rst(&self) -> DVec2 {
        self.pt2drst
    }

    /// OCCT ParameterOnRst() (BRepBlend_SurfRstEvolRad.cxx L530-533).
    pub fn parameter_on_rst(&self) -> f64 {
        self.prmrst
    }

    /// OCCT IsTangencyPoint() (BRepBlend_SurfRstEvolRad.cxx L537-540).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnS() (BRepBlend_SurfRstEvolRad.cxx L544-551).
    pub fn tangent_on_s(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstEvolRad::TangentOnS");
        }
        self.tgs
    }

    /// OCCT Tangent2dOnS() (BRepBlend_SurfRstEvolRad.cxx L555-562).
    pub fn tangent_2d_on_s(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstEvolRad::Tangent2dOnS");
        }
        self.tg2ds
    }

    /// OCCT TangentOnRst() (BRepBlend_SurfRstEvolRad.cxx L566-573).
    pub fn tangent_on_rst(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstEvolRad::TangentOnRst");
        }
        self.tgrst
    }

    /// OCCT Tangent2dOnRst() (BRepBlend_SurfRstEvolRad.cxx L577-584).
    pub fn tangent_2d_on_rst(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BRepBlend_SurfRstEvolRad::Tangent2dOnRst");
        }
        self.tg2drst
    }

    /// OCCT Decroch(Sol, NS, TgS) (BRepBlend_SurfRstEvolRad.cxx L588-640).
    // The OCCT body assigns NSInPlane at L596 and reverses it at L605
    // without any further read (a dead store in the original); the rcad
    // port keeps the statements for parity.
    #[allow(unused_assignments)]
    pub fn decroch(&self, sol: &[f64], ns_out: &mut DVec3, tg_s: &mut DVec3) -> bool {
        // OCCT L595-596: surf->D1(Sol(1), Sol(2), bid, d1u, d1v);
        // NS = NSInPlane = d1u.Crossed(d1v).
        let (bid, d1u, d1v) = self.surf.derivatives(sol[0], sol[1]);
        let ns = d1u.cross(d1v);
        *ns_out = ns;
        let mut ns_in_plane = ns;

        // OCCT L598-600.
        let norm = self.nplan.cross(ns).length();
        let unsurnorm = 1.0 / norm;
        ns_in_plane = (self.nplan.dot(ns) * unsurnorm) * self.nplan + (-unsurnorm) * ns;

        // OCCT L602: Center.SetXYZ(bid.XYZ() + ray * NSInPlane.XYZ()).
        let center = bid + self.ray * ns_in_plane;
        // OCCT L603-606.
        if self.choix > 2 {
            ns_in_plane = -ns_in_plane;
        }
        // OCCT L607-611: TgS = nplan.Crossed(gp_Vec(Center, bid)).
        *tg_s = self.nplan.cross(bid - center);
        if self.choix % 2 == 1 {
            *tg_s = -*tg_s;
        }
        // OCCT L612-614: rstref->Value(Sol(3)).Coord(u, v);
        // surfref->D1(u, v, bid, d1u, d1v).
        let (u, v) = {
            let p = self.rstref.expect("rstref").point_at(sol[2]);
            (p.x, p.y)
        };
        let (bid, d1u, d1v) = self.surfref.expect("surfref").derivatives(u, v);

        // OCCT L615-618.
        let n_rst = d1u.cross(d1v);
        let norm = self.nplan.cross(n_rst).length();
        let unsurnorm = 1.0 / norm;
        let mut n_rst_in_plane =
            (self.nplan.dot(n_rst) * unsurnorm) * self.nplan + (-unsurnorm) * n_rst;
        // OCCT L619: gp_Vec centptrst(Center, bid).
        let centptrst = bid - center;
        // OCCT L620-623.
        if centptrst.dot(n_rst_in_plane) < 0.0 {
            n_rst_in_plane = -n_rst_in_plane;
        }
        // OCCT L624-628.
        let mut tg_rst = self.nplan.cross(centptrst);
        if self.choix % 2 == 1 {
            tg_rst = -tg_rst;
        }

        // OCCT L630-639.
        let mut nt = n_rst_in_plane.length();
        nt *= tg_rst.length();
        if nt.abs() < 1.0e-7 {
            return false; // Singularity or Incoherence.
        }
        let mut dot = n_rst_in_plane.dot(tg_rst);
        dot /= nt;

        dot < 1.0e-10
    }

    /// OCCT Set(Choix) (BRepBlend_SurfRstEvolRad.cxx L644-661).
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

    /// OCCT Set(BlendFunc_SectionShape) (BRepBlend_SurfRstEvolRad.cxx
    /// L665-668).
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_s_shape = type_section;
    }

    /// OCCT Section(Param, U, V, W, Pdeb, Pfin, C)
    /// (BRepBlend_SurfRstEvolRad.cxx L672-723).
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
        // OCCT L685-687: tguide->D1(Param, ptgui, d1gui); np = d1gui.Normalized();
        // ray = sg1 * tevol->Value(Param).
        self.ptgui = self.tguide().point_at(param);
        self.d1gui = self.tguide().derivative_at(param);
        let mut np = self.d1gui.normalize_or_zero();
        // OCCT L687: ray = sg1 * tevol->Value(Param) — the law read and the
        // self.ray write are split across statements (Rc<RefCell> borrow).
        let law_ray = self.tevol().borrow_mut().value(param);
        self.ray = self.sg1 * law_ray;

        // OCCT L689-690.
        let (p, d1u1, d1v1) = self.surf.derivatives(u, v);
        self.pts = p;
        self.ptrst = self.cons_value(w);

        // OCCT L692-697.
        let ns = d1u1.cross(d1v1);
        let norm = self.nplan.cross(ns).length();
        let ns = (self.nplan.dot(ns) / norm) * self.nplan + (-1.0 / norm) * ns;
        let center = self.pts + self.ray * ns;
        c.radius = self.ray.abs();

        // OCCT L699-706.
        let ns = if self.ray > 0.0 { -ns } else { ns };
        if self.choix % 2 != 0 {
            np = -np;
        }
        // OCCT L707: C.SetPosition(gp_Ax2(Center, np, ns))
        // (gp_Ax2(P, N, Vx): the Y direction is N ^ Vx).
        c.center = center;
        c.normal = np;
        c.x_dir = ns;
        c.y_dir = np.cross(ns);

        // OCCT L709-710: Pdeb = 0.; Pfin = ElCLib::Parameter(C, ptrst).
        *pdeb = 0.0;
        *pfin = elclib_circle_parameter(c, self.ptrst);

        // Test negative and almost null angles : Single Case
        // OCCT L713-718.
        if *pfin > 1.5 * std::f64::consts::PI {
            np = -np;
            c.normal = np;
            c.y_dir = np.cross(ns);
            *pfin = elclib_circle_parameter(c, self.ptrst);
        }
        // OCCT L719-722.
        if *pfin < p_confusion() {
            *pfin += p_confusion();
        }
    }

    /// OCCT IsRational() (BRepBlend_SurfRstEvolRad.cxx L727-730).
    pub fn is_rational(&self) -> bool {
        self.my_s_shape == BlendFuncSectionShape::Rational
            || self.my_s_shape == BlendFuncSectionShape::QuasiAngular
    }

    /// OCCT GetSectionSize() (BRepBlend_SurfRstEvolRad.cxx L734-737).
    pub fn get_section_size(&self) -> f64 {
        self.maxang * self.ray.abs()
    }

    /// OCCT GetMinimalWeight(Weights) (BRepBlend_SurfRstEvolRad.cxx
    /// L741-745).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        blend_func_get_minimal_weights(self.my_s_shape, self.my_t_conv, self.minang, self.maxang, weigths);
        // It is supposed that it does not depend on the Radius!
    }

    /// OCCT NbIntervals(S) (BRepBlend_SurfRstEvolRad.cxx L749-768).
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

    /// OCCT Intervals(T, S) (BRepBlend_SurfRstEvolRad.cxx L772-796).
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
    /// (BRepBlend_SurfRstEvolRad.cxx L800-804).
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
    /// (BRepBlend_SurfRstEvolRad.cxx L808-821) — tolerances used for
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

    /// OCCT Knots(TKnots) (BRepBlend_SurfRstEvolRad.cxx L825-828).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        geomfill_knots(self.my_t_conv, tknots);
    }

    /// OCCT Mults(TMults) (BRepBlend_SurfRstEvolRad.cxx L832-835).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        geomfill_mults(self.my_t_conv, tmults);
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)
    /// (BRepBlend_SurfRstEvolRad.cxx L839-1069) — used for the first and
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

        // OCCT L865-872.
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        self.d2gui = self.tguide().derivative2_at(prm);
        let (mut ray, mut a_dray) = (0.0f64, 0.0f64);
        self.tevol().borrow_mut().d1(prm, &mut ray, &mut a_dray);
        ray = self.sg1 * ray;
        a_dray = self.sg1 * a_dray;
        self.ray = ray;
        self.dray = a_dray;
        self.normtg = self.d1gui.length();
        self.nplan = self.d1gui.normalize_or_zero();
        // OCCT L872: dnplan.SetLinearForm(1. / normtg, d2gui,
        // -1. / normtg * (nplan.Dot(d2gui)), nplan).
        let mut dnplan = (1.0 / self.normtg) * self.d2gui
            + (-(1.0 / self.normtg) * self.nplan.dot(self.d2gui)) * self.nplan;

        // OCCT L874-876.
        let mut sol = [0.0f64; 3];
        let (u1, v1) = p.parameters_on_s();
        sol[0] = u1;
        sol[1] = v1;
        self.prmrst = p.parameter_on_c();
        sol[2] = self.prmrst;
        self.pt2drst = self.rst.point_at(self.prmrst);

        // OCCT L878: Values(sol, valsol, gradsol).
        let mut valsol = [0.0f64; 3];
        let mut gradsol = vec![vec![0.0f64; 3]; 3];
        self.values(&sol, &mut valsol, &mut gradsol);

        // OCCT L880-881.
        let (ptsv, d1u1, d1v1, d2u1, d2v1, d2uv1) = self.surf.derivatives2(sol[0], sol[1]);
        self.pts = ptsv;
        let (ptrst_v, d1) = self.cons_d1(sol[2]);
        self.ptrst = ptrst_v;

        // OCCT L883-887.
        let mut secmember = [0.0f64; 3];
        let temp = self.pts - self.ptgui;
        secmember[0] = self.normtg - dnplan.dot(temp);

        let temp = self.ptrst - self.ptgui;
        secmember[1] = self.normtg - dnplan.dot(temp);

        // OCCT L889-899.
        let mut ns = d1u1.cross(d1v1);
        let ncrossns = self.nplan.cross(ns);
        let ndotns = self.nplan.dot(ns);
        let mut norm = ncrossns.length();
        if norm < EPS {
            norm = 1.0; // Not enough, but it is not necessary to stop
        }

        // Derivative of n1 corresponding to w
        // OCCT L903-909.
        let grosterme = ncrossns.dot(dnplan.cross(ns)) / norm / norm;
        let mut dnw = ((dnplan.dot(ns) - grosterme * ndotns) / norm) * self.nplan
            + (ndotns / norm) * dnplan
            + (grosterme / norm) * ns;

        // OCCT L911-915.
        let temp = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;
        let resul = ray * temp + (self.pts - self.ptrst);
        // secmember(3) = -2.*ray*(dnw.Dot(resul)); // jag 950105 il manquait ray
        secmember[2] = -2.0 * ray * dnw.dot(resul) - 2.0 * a_dray * temp.dot(resul)
            + 2.0 * ray * a_dray;

        // OCCT L916-937: math_Gauss Resol(gradsol) and the SVD fallback.
        let mut a = MatD::new(3, 3);
        for r in 1..=3 {
            for c in 1..=3 {
                a.set(r, c, gradsol[r - 1][c - 1]);
            }
        }
        let resol = MathGauss::new(&a);
        let mut istgt: bool;
        if resol.is_done() {
            // OCCT L920: Resol.Solve(secmember).
            let mut x = VecD::new(3);
            for i in 1..=3 {
                x.set(i, secmember[i - 1]);
            }
            resol.solve(&mut x);
            for i in 1..=3 {
                secmember[i - 1] = x.get(i);
            }
            istgt = false;
        } else {
            // OCCT L925-936: math_SVD SingRS(gradsol);
            // if (SingRS.IsDone()) { math_Vector DEDT = secmember;
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

        // OCCT L939-981.
        let ns2: DVec3;
        let mut dn2w = DVec3::ZERO;
        if !istgt {
            // OCCT L942-943.
            self.tgs = secmember[0] * d1u1 + secmember[1] * d1v1;
            self.tgrst = secmember[2] * d1;

            // Derivative of n1 corresponding to u1
            // OCCT L945-952.
            let temp = d2u1.cross(d1v1) + d1u1.cross(d2uv1);
            let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
            let resulu = (-(grosterme * ndotns - self.nplan.dot(temp)) / norm) * self.nplan
                + (grosterme / norm) * ns
                + (-1.0 / norm) * temp;

            // Derivative of n1 corresponding to v1
            // OCCT L954-962.
            let temp = d2uv1.cross(d1v1) + d1u1.cross(d2v1);
            let grosterme = ncrossns.dot(self.nplan.cross(temp)) / norm / norm;
            let resulv = (-(grosterme * ndotns - self.nplan.dot(temp)) / norm) * self.nplan
                + (grosterme / norm) * ns
                + (-1.0 / norm) * temp;

            // OCCT L964: dnw.SetLinearForm(secmember(1), resulu, secmember(2),
            // resulv, dnw).
            dnw = secmember[0] * resulu + secmember[1] * resulv + dnw;
            // OCCT L965.
            ns = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;

            // OCCT L967-968: dn2w.SetLinearForm(ray, dnw, -1., tgrst, tgs);
            // dn2w.SetLinearForm(aDray, ns, dn2w).
            dn2w = ray * dnw + (-1.0) * self.tgrst + self.tgs;
            dn2w = a_dray * ns + dn2w;
            norm = resul.length();
            dn2w /= norm;
            ns2 = -resul.normalize_or_zero();
            dn2w = ns2.dot(dn2w) * ns2 + (-1.0) * dn2w;

            // OCCT L974: istgt = false — a dead re-assignment in the
            // original (istgt is already false here); noted for parity.
        } else {
            // OCCT L978-980.
            ns = (ndotns / norm) * self.nplan + (-1.0 / norm) * ns;
            ns2 = -resul.normalize_or_zero();
            istgt = true;
        }

        // Tops 2D
        // OCCT L985-994.
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
        // OCCT L997-1011.
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
        // OCCT L1014: Center.SetXYZ(pts.XYZ() + ray * ns.XYZ()).
        let center = self.pts + ray * ns;
        // OCCT L1015-1018: tgct.SetLinearForm(ray, dnw, aDray, ns, tgs).
        let tgct: DVec3;
        if !istgt {
            tgct = ray * dnw + a_dray * ns + self.tgs;
        } else {
            tgct = DVec3::ZERO; // the OCCT tgct is only read on the !istgt path
        }

        // OCCT L1020-1027.
        if ray > 0.0 {
            ns = -ns;
            if !istgt {
                dnw = -dnw;
            }
        }
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
            dnplan = -dnplan;
        }

        // OCCT L1033-1063.
        if !istgt {
            // OCCT L1035-1042: to avoid std::abs(dray) some lines below.
            let rayprim = if ray < 0.0 { -a_dray } else { a_dray };

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
                ray.abs(),
                rayprim,
                center,
                tgct,
                poles,
                d_poles,
                weigths,
                d_weigths,
            )
        } else {
            // OCCT L1066.
            get_circle(
                tconv(self.my_t_conv),
                ns,
                ns2,
                self.nplan,
                self.pts,
                self.ptrst,
                ray.abs(),
                center,
                poles,
                weigths,
            );
            false
        }
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weigths, DWeigths, D2Weigths) (BRepBlend_SurfRstEvolRad.cxx
    /// L1073-1085) — returns false.
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

    /// OCCT Section(P, Poles, Poles2d, Weigths)
    /// (BRepBlend_SurfRstEvolRad.cxx L1089-1149).
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

        // OCCT L1104-1107: tguide->D1(prm, ptgui, d1gui); ray =
        // tevol->Value(prm); ray = sg1 * ray; nplan = d1gui.Normalized().
        self.ptgui = self.tguide().point_at(prm);
        self.d1gui = self.tguide().derivative_at(prm);
        // OCCT L1105-1106: ray = tevol->Value(prm); ray = sg1 * ray — the
        // law read and the self.ray write are split across statements
        // (Rc<RefCell> borrow).
        let law_ray = self.tevol().borrow_mut().value(prm);
        self.ray = self.sg1 * law_ray;
        self.nplan = self.d1gui.normalize_or_zero();

        // OCCT L1109-1111.
        let (u1, v1) = p.parameters_on_s();
        let w = p.parameter_on_c(); // jlr : point on curve not on surface
        let pt2d = self.rst.point_at(w);

        // OCCT L1113-1114.
        let (ptsv, d1u1, d1v1) = self.surf.derivatives(u1, v1);
        self.pts = ptsv;
        self.ptrst = self.cons_value(w);

        // OCCT L1116.
        self.distmin = self.distmin.min(self.pts.distance(self.ptrst));

        // OCCT L1118-1119.
        poles_2d[0] = DVec2::new(u1, v1);
        poles_2d[poles_2d.len() - 1] = DVec2::new(pt2d.x, pt2d.y);

        // Linear case
        // OCCT L1122-1129.
        if self.my_s_shape == BlendFuncSectionShape::Linear {
            poles[low] = self.pts;
            poles[upp] = self.ptrst;
            weigths[low] = 1.0;
            weigths[upp] = 1.0;
            return;
        }

        // OCCT L1131-1136.
        let mut ns = d1u1.cross(d1v1);
        let norm = self.nplan.cross(ns).length();
        ns = (self.nplan.dot(ns) / norm) * self.nplan + (-1.0 / norm) * ns;
        let center = self.pts + self.ray * ns;

        // OCCT L1138-1146.
        let ns2 = (self.ptrst - center).normalize_or_zero();
        if self.ray > 0.0 {
            ns = -ns;
        }
        if self.choix % 2 != 0 {
            self.nplan = -self.nplan;
        }

        // OCCT L1148: GeomFill::GetCircle(myTConv, ns, ns2, nplan, pts,
        // ptrst, std::abs(ray), Center, Poles, Weigths).
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

    /// OCCT Resolution(IC2d, Tol, TolU, TolV)
    /// (BRepBlend_SurfRstEvolRad.cxx L1151-1166).
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

impl<'a> FunctionSetWithDerivatives for BlendSurfRstEvolRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendSurfRstEvolRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendSurfRstEvolRad::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendSurfRstEvolRad::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfRstEvolRad::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendSurfRstEvolRad::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendSurfRstEvolRad<'a> {
    fn set_param(&mut self, param: f64) {
        BlendSurfRstEvolRad::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendSurfRstEvolRad::set_interval(self, first, last)
    }

    fn pnt1(&self) -> DVec3 {
        BlendSurfRstFunction::pnt1(self)
    }

    fn pnt2(&self) -> DVec3 {
        BlendSurfRstFunction::pnt2(self)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendSurfRstEvolRad::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendSurfRstEvolRad::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendSurfRstEvolRad::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendSurfRstEvolRad::get_minimal_distance(self)
    }

    fn is_rational(&self) -> bool {
        BlendSurfRstEvolRad::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendSurfRstEvolRad::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendSurfRstEvolRad::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendSurfRstEvolRad::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendSurfRstEvolRad::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendSurfRstEvolRad::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendSurfRstEvolRad::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendSurfRstEvolRad::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendSurfRstEvolRad::mults(self, tmults)
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
        BlendSurfRstEvolRad::section_d1(
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
        BlendSurfRstEvolRad::section_simple(self, p, poles, poles_2d, weigths)
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
        BlendSurfRstEvolRad::section_d2(
            self, p, poles, d_poles, d2_poles, poles_2d, d_poles_2d, d2_poles_2d, weigths,
            d_weigths, d2_weigths,
        )
    }

    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        BlendSurfRstEvolRad::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendSurfRstFunction for BlendSurfRstEvolRad<'a> {
    fn nb_variables(&self) -> usize {
        BlendSurfRstEvolRad::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendSurfRstEvolRad::nb_equations(self)
    }

    fn point_on_s(&self) -> DVec3 {
        BlendSurfRstEvolRad::point_on_s(self)
    }

    fn point_on_rst(&self) -> DVec3 {
        BlendSurfRstEvolRad::point_on_rst(self)
    }

    fn pnt2d_on_s(&self) -> DVec2 {
        BlendSurfRstEvolRad::pnt2d_on_s(self)
    }

    fn pnt2d_on_rst(&self) -> DVec2 {
        BlendSurfRstEvolRad::pnt2d_on_rst(self)
    }

    fn parameter_on_rst(&self) -> f64 {
        BlendSurfRstEvolRad::parameter_on_rst(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendSurfRstEvolRad::is_tangency_point(self)
    }

    fn tangent_on_s(&self) -> DVec3 {
        BlendSurfRstEvolRad::tangent_on_s(self)
    }

    fn tangent_2d_on_s(&self) -> DVec2 {
        BlendSurfRstEvolRad::tangent_2d_on_s(self)
    }

    fn tangent_on_rst(&self) -> DVec3 {
        BlendSurfRstEvolRad::tangent_on_rst(self)
    }

    fn tangent_2d_on_rst(&self) -> DVec2 {
        BlendSurfRstEvolRad::tangent_2d_on_rst(self)
    }

    fn decroch(&self, sol: &[f64], ns: &mut DVec3, tg_s: &mut DVec3) -> bool {
        BlendSurfRstEvolRad::decroch(self, sol, ns, tg_s)
    }
}
