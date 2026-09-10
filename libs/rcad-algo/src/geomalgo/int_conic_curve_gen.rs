//! OCCT IntCurve_IntConicCurveGen (TKGeomAlgo IntCurve package) —
//! intersection of a conic with a parametric curve.
//!
//! 1:1 translation of `IntCurve_IntConicCurveGen.gxx` (L1-92, the C/E/P/H
//! constructors) + `IntCurve_IntConicCurveGen.lxx` (L23-130, the default
//! constructor, the Lin constructor, the five Perform overloads and the
//! IConicTool Perform).  The template parameters map to the curve type `C`
//! plus the [`ParTool`]/[`ProjectOnPCurveTool`] traits of
//! [`super::int_imp_par_gen`]:
//! - Geom2dInt_TheIntConicCurveOfGInter — ThePCurve = Adaptor2d_Curve2d
//!   (`C = dyn Curve2dAdaptor`), TheIntersector =
//!   Geom2dInt_TheIntersectorOfTheIntConicCurveOfGInter;
//! - HLRBRep_TheIntConicCurveOfCInter — ThePCurve = HLRBRep_CurvePtr,
//!   TheIntersector = HLRBRep_TheIntersectorOfTheIntConicCurveOfCInter
//!   (Stage 3a supplies the tool impls over the same traits).

use rcad_kernel::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};

use super::geom2d_int::IConicTool;
use super::int_imp_par_gen::{Intersector, ParTool, ProjectOnPCurveTool};
use super::int_res2d::{Domain as Res2dDomain, IntersectionBase};

/// 2·π (OCCT M_PI + M_PI).
const PI2: f64 = std::f64::consts::TAU;

/// OCCT IntCurve_IntConicCurveGen — the conic x parametric-curve
/// intersection, generic over the ThePCurve/TheIntersector template
/// parameters.
pub struct IntConicCurveGen<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> {
    pub base: IntersectionBase,
    _tools: std::marker::PhantomData<fn(&C, &PT, &JT)>,
}

impl<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> Clone for IntConicCurveGen<C, PT, JT> {
    fn clone(&self) -> Self {
        IntConicCurveGen {
            base: self.base.clone(),
            _tools: std::marker::PhantomData,
        }
    }
}

impl<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> std::fmt::Debug
    for IntConicCurveGen<C, PT, JT>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("IntConicCurveGen").field("base", &self.base).finish()
    }
}

impl<C: ?Sized, PT: ParTool<C>, JT: ProjectOnPCurveTool<C>> IntConicCurveGen<C, PT, JT> {
    /// OCCT IntCurve_IntConicCurveGen() (IntCurve_IntConicCurveGen.lxx L30).
    pub fn bare() -> Self {
        IntConicCurveGen {
            base: IntersectionBase::new(),
            _tools: std::marker::PhantomData,
        }
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Lin2d& L, D1, PCurve, D2,
    /// TolConf, Tol) (IntCurve_IntConicCurveGen.lxx L33-42).
    #[allow(clippy::too_many_arguments)]
    pub fn new_line(
        line: &Line2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntConicCurveGen::bare();
        r.perform_line(line, d1, pcurve, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Circ2d& C, D1, PCurve, D2,
    /// TolConf, Tol) (IntCurve_IntConicCurveGen.gxx L24-42).
    #[allow(clippy::too_many_arguments)]
    pub fn new_circle(
        c: &Circle2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntConicCurveGen::bare();
        let tool = IConicTool::new_circle(c);
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            r.perform_imp(&tool, &d, pcurve, d2, tol_conf, tol);
        } else {
            r.perform_imp(&tool, d1, pcurve, d2, tol_conf, tol);
        }
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Elips2d& E, D1, PCurve, D2,
    /// TolConf, Tol) (IntCurve_IntConicCurveGen.gxx L46-64).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ellipse(
        e: &Ellipse2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntConicCurveGen::bare();
        let tool = IConicTool::new_ellipse(e);
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            r.perform_imp(&tool, &d, pcurve, d2, tol_conf, tol);
        } else {
            r.perform_imp(&tool, d1, pcurve, d2, tol_conf, tol);
        }
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Parab2d& Prb, D1, PCurve, D2,
    /// TolConf, Tol) (IntCurve_IntConicCurveGen.gxx L68-77).
    #[allow(clippy::too_many_arguments)]
    pub fn new_parabola(
        p: &Parabola2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntConicCurveGen::bare();
        r.perform_imp(&IConicTool::new_parabola(p), d1, pcurve, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntConicCurveGen(const gp_Hypr2d& H, D1, PCurve, D2,
    /// TolConf, Tol) (IntCurve_IntConicCurveGen.gxx L81-90).
    #[allow(clippy::too_many_arguments)]
    pub fn new_hyperbola(
        h: &Hyperbola2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntConicCurveGen::bare();
        r.perform_imp(&IConicTool::new_hyperbola(h), d1, pcurve, d2, tol_conf, tol);
        r
    }

    /// OCCT Perform(const gp_Lin2d& L, D1, PCurve, D2, TolConf, Tol)
    /// (IntCurve_IntConicCurveGen.lxx L44-53).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_line(
        &mut self,
        line: &Line2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.perform_imp(&IConicTool::new_line(line), d1, pcurve, d2, tol_conf, tol);
    }

    /// OCCT Perform(const gp_Circ2d& C, D1, PCurve, D2, TolConf, Tol)
    /// (IntCurve_IntConicCurveGen.lxx L56-74).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_circle(
        &mut self,
        c: &Circle2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            self.perform_imp(&IConicTool::new_circle(c), &d, pcurve, d2, tol_conf, tol);
        } else {
            self.perform_imp(&IConicTool::new_circle(c), d1, pcurve, d2, tol_conf, tol);
        }
    }

    /// OCCT Perform(const gp_Elips2d& E, D1, PCurve, D2, TolConf, Tol)
    /// (IntCurve_IntConicCurveGen.lxx L77-94).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_ellipse(
        &mut self,
        e: &Ellipse2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        if !d1.is_closed() {
            let mut d = d1.clone();
            d.set_equivalent_parameters(d1.first_parameter(), d1.first_parameter() + PI2);
            self.perform_imp(&IConicTool::new_ellipse(e), &d, pcurve, d2, tol_conf, tol);
        } else {
            self.perform_imp(&IConicTool::new_ellipse(e), d1, pcurve, d2, tol_conf, tol);
        }
    }

    /// OCCT Perform(const gp_Parab2d& Prb, D1, PCurve, D2, TolConf, Tol)
    /// (IntCurve_IntConicCurveGen.lxx L97-105).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_parabola(
        &mut self,
        p: &Parabola2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.perform_imp(&IConicTool::new_parabola(p), d1, pcurve, d2, tol_conf, tol);
    }

    /// OCCT Perform(const gp_Hypr2d& H, D1, PCurve, D2, TolConf, Tol)
    /// (IntCurve_IntConicCurveGen.lxx L108-116).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_hyperbola(
        &mut self,
        h: &Hyperbola2d,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.perform_imp(&IConicTool::new_hyperbola(h), d1, pcurve, d2, tol_conf, tol);
    }

    /// OCCT Perform(const IntCurve_IConicTool& ICurve, D1, PCurve, D2,
    /// TolConf, Tol) (IntCurve_IntConicCurveGen.lxx L119-130).
    #[allow(clippy::too_many_arguments)]
    fn perform_imp(
        &mut self,
        imp_tool: &IConicTool,
        d1: &Res2dDomain,
        pcurve: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let mut myintersection = Intersector::<C, PT, JT>::new();
        myintersection
            .base
            .set_reversed_parameters(self.base.reversed_parameters());
        myintersection.perform(imp_tool, d1, pcurve, d2, tol_conf, tol);
        self.base.set_values(&myintersection.base);
    }
}
