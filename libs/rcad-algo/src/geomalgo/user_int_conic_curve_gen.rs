//! OCCT IntCurve_UserIntConicCurveGen (TKGeomAlgo IntCurve package) —
//! intersection of a conic with a parametric curve, dispatching on the
//! parametric curve kind (conic kinds go conic x conic, other kinds go
//! through IntConicCurveGen / IntImpParGen_Intersector).
//!
//! 1:1 translation of `IntCurve_UserIntConicCurveGen.gxx` (L1-889): five
//! constructors, five Perform overloads (each with the NbIntervals>1
//! per-interval composite loop) and five InternalPerform overloads (each
//! dispatching on ThePCurveTool::GetType).  The template parameters map to
//! the [`PCurveTool`] trait (ThePCurveTool / ThePCurve) plus the
//! [`ParTool`]/[`ProjectOnPCurveTool`] traits consumed by the inner
//! IntConicCurveGen:
//! - HLRBRep_IntConicCurveOfCInter — ThePCurve = HLRBRep_CurvePtr,
//!   ThePCurveTool = HLRBRep_CurveTool, intconicurv =
//!   HLRBRep_TheIntConicCurveOfCInter (Stage 3a supplies the HLRBRep tool
//!   impls over the same traits).

use glam::DVec2;
use rcad_kernel::geom::{Circle2d, Ellipse2d, Hyperbola2d, Line2d, Parabola2d};

use super::geom2d_int::Curve2dType;
use super::int_conic_conic::IntConicConic;
use super::int_conic_curve_gen::IntConicCurveGen;
use super::int_imp_par_gen::{ParTool, ProjectOnPCurveTool};
use super::int_res2d::{Domain as Res2dDomain, IntersectionBase};

/// OCCT Precision::Infinite() (Precision.hxx).
const PRECISION_INFINITE: f64 = 2.0e100;
/// OCCT Standard_Real RealEpsilon() (DBL_EPSILON).
const REAL_EPSILON: f64 = f64::EPSILON;

/// OCCT `ThePCurveTool` template parameter of IntCurve_UserIntConicCurveGen
/// (HLRBRep_CurveTool for the CInter instantiation) — the static accessors
/// over the parametric curve, including the conic kind extractors.
pub trait PCurveTool<C: ?Sized> {
    /// OCCT ThePCurveTool::NbIntervals(C).
    fn nb_intervals(c: &C) -> i32;
    /// OCCT ThePCurveTool::Intervals(C, Tab) — fills Tab(1, NbIntervals+1).
    fn intervals(c: &C, tab: &mut [f64]);
    /// OCCT ThePCurveTool::GetInterval(C, Index, Tab, U1, U2).
    fn get_interval(c: &C, index: usize, tab: &[f64]) -> (f64, f64);
    /// OCCT ThePCurveTool::FirstParameter(C).
    fn first_parameter(c: &C) -> f64;
    /// OCCT ThePCurveTool::LastParameter(C).
    fn last_parameter(c: &C) -> f64;
    /// OCCT ThePCurveTool::Value(C, U).
    fn value(c: &C, u: f64) -> DVec2;
    /// OCCT ThePCurveTool::GetType(C).
    fn get_type(c: &C) -> Curve2dType;
    /// OCCT ThePCurveTool::Line(C).
    fn line(c: &C) -> Line2d;
    /// OCCT ThePCurveTool::Circle(C).
    fn circle(c: &C) -> Circle2d;
    /// OCCT ThePCurveTool::Ellipse(C).
    fn ellipse(c: &C) -> Ellipse2d;
    /// OCCT ThePCurveTool::Parabola(C).
    fn parabola(c: &C) -> Parabola2d;
    /// OCCT ThePCurveTool::Hyperbola(C).
    fn hyperbola(c: &C) -> Hyperbola2d;
}

/// OCCT IntCurve_UserIntConicCurveGen — conic x parametric-curve
/// intersection with kind dispatch.
pub struct UserIntConicCurveGen<C: ?Sized, PT: ParTool<C> + PCurveTool<C>, JT: ProjectOnPCurveTool<C>>
{
    pub base: IntersectionBase,
    /// OCCT param1inf.
    param1inf: f64,
    /// OCCT param1sup.
    param1sup: f64,
    /// OCCT param2inf.
    param2inf: f64,
    /// OCCT param2sup.
    param2sup: f64,
    /// OCCT intconiconi.
    intconiconi: IntConicConic,
    /// OCCT intconicurv.
    intconicurv: IntConicCurveGen<C, PT, JT>,
    _tools: std::marker::PhantomData<fn(&C, &PT, &JT)>,
}

impl<C: ?Sized, PT: ParTool<C> + PCurveTool<C>, JT: ProjectOnPCurveTool<C>> Clone for UserIntConicCurveGen<C, PT, JT>
{
    fn clone(&self) -> Self {
        UserIntConicCurveGen {
            base: self.base.clone(),
            param1inf: self.param1inf,
            param1sup: self.param1sup,
            param2inf: self.param2inf,
            param2sup: self.param2sup,
            intconiconi: self.intconiconi.clone(),
            intconicurv: self.intconicurv.clone(),
            _tools: std::marker::PhantomData,
        }
    }
}

impl<C: ?Sized, PT: ParTool<C> + PCurveTool<C>, JT: ProjectOnPCurveTool<C>> std::fmt::Debug for UserIntConicCurveGen<C, PT, JT>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("UserIntConicCurveGen")
            .field("base", &self.base)
            .finish()
    }
}

impl<C: ?Sized, PT: ParTool<C> + PCurveTool<C>, JT: ProjectOnPCurveTool<C>> UserIntConicCurveGen<C, PT, JT>
{
    /// OCCT IntCurve_UserIntConicCurveGen() (gxx L25-28).
    pub fn new() -> Self {
        UserIntConicCurveGen {
            base: IntersectionBase::new(),
            param1inf: 0.0,
            param1sup: 0.0,
            param2inf: 0.0,
            param2sup: 0.0,
            intconiconi: IntConicConic::new(),
            intconicurv: IntConicCurveGen::bare(),
            _tools: std::marker::PhantomData,
        }
    }

    /// OCCT IntCurve_UserIntConicCurveGen(const gp_Lin2d& Lin1, D1, C2, D2,
    /// TolConf, Tol) (gxx L32-40).
    #[allow(clippy::too_many_arguments)]
    pub fn new_line(
        lin1: &Line2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = UserIntConicCurveGen::new();
        r.perform_line(lin1, d1, c2, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_UserIntConicCurveGen(const gp_Circ2d& Circ1, D1, C2, D2,
    /// TolConf, Tol) (gxx L44-52).
    #[allow(clippy::too_many_arguments)]
    pub fn new_circle(
        circ1: &Circle2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = UserIntConicCurveGen::new();
        r.perform_circle(circ1, d1, c2, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_UserIntConicCurveGen(const gp_Parab2d& Parab1, D1, C2,
    /// D2, TolConf, Tol) (gxx L56-64).
    #[allow(clippy::too_many_arguments)]
    pub fn new_parabola(
        parab1: &Parabola2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = UserIntConicCurveGen::new();
        r.perform_parabola(parab1, d1, c2, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_UserIntConicCurveGen(const gp_Elips2d& Elips1, D1, C2,
    /// D2, TolConf, Tol) (gxx L68-76).
    #[allow(clippy::too_many_arguments)]
    pub fn new_ellipse(
        elips1: &Ellipse2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = UserIntConicCurveGen::new();
        r.perform_ellipse(elips1, d1, c2, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_UserIntConicCurveGen(const gp_Hypr2d& Hyper1, D1, C2,
    /// D2, TolConf, Tol) (gxx L80-88).
    #[allow(clippy::too_many_arguments)]
    pub fn new_hyperbola(
        hyper1: &Hyperbola2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = UserIntConicCurveGen::new();
        r.perform_hyperbola(hyper1, d1, c2, d2, tol_conf, tol);
        r
    }

    /// OCCT Perform(const gp_Lin2d& Lin1, D1, C2, D2, TolConf, Tol)
    /// (gxx L94-146).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_line(
        &mut self,
        lin1: &Line2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.reset_fields();
        let nb_inter_c2 = <PT as PCurveTool<C>>::nb_intervals(c2);
        if nb_inter_c2 > 1 {
            let mut param_inf;
            let mut param_sup;
            let d2_first_param = d2.first_parameter();
            let d2_last_param = d2.last_parameter();
            let mut ok = true;
            self.param1inf = if d1.has_first_point() {
                d1.first_parameter()
            } else {
                -PRECISION_INFINITE
            };
            self.param1sup = if d1.has_last_point() {
                d1.last_parameter()
            } else {
                PRECISION_INFINITE
            };
            self.param2inf = <PT as PCurveTool<C>>::first_parameter(c2);
            self.param2sup = <PT as PCurveTool<C>>::last_parameter(c2);
            let mut domain_c2_num_inter = Res2dDomain::infinite();

            // NCollection_Array1<double> Tab2(1, NbInterC2 + 1) — 0-based
            // storage of the 1-based array (logical index i at tab[i - 1]).
            let mut tab2 = vec![0.0; (nb_inter_c2 + 2) as usize];
            <PT as PCurveTool<C>>::intervals(c2, &mut tab2);

            let mut num_inter_c2 = 1;
            while ok && num_inter_c2 <= nb_inter_c2 {
                (param_inf, param_sup) = <PT as PCurveTool<C>>::get_interval(c2, num_inter_c2 as usize, &tab2);
                if (param_inf > d2_last_param) || (param_sup < d2_first_param) {
                    ok = false;
                } else {
                    if param_inf < d2_first_param {
                        param_inf = d2_first_param;
                    }
                    if param_sup > d2_last_param {
                        param_sup = d2_last_param;
                    }
                    if (param_sup - param_inf) > REAL_EPSILON {
                        domain_c2_num_inter.set_values_bounded(
                            <PT as PCurveTool<C>>::value(c2, param_inf),
                            param_inf,
                            d2.first_tolerance(),
                            <PT as PCurveTool<C>>::value(c2, param_sup),
                            param_sup,
                            d2.last_tolerance(),
                        );
                        self.internal_perform_line(lin1, d1, c2, &domain_c2_num_inter, tol_conf, tol, true);
                    }
                }
                num_inter_c2 += 1;
            }
        } else {
            self.internal_perform_line(lin1, d1, c2, d2, tol_conf, tol, false);
        }
    }

    /// OCCT Perform(const gp_Circ2d& Circ1, D1, C2, D2, TolConf, Tol)
    /// (gxx L150-203).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_circle(
        &mut self,
        circ1: &Circle2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.reset_fields();
        let nb_inter_c2 = <PT as PCurveTool<C>>::nb_intervals(c2);
        if nb_inter_c2 > 1 {
            let mut param_inf;
            let mut param_sup;
            let d2_first_param = d2.first_parameter();
            let d2_last_param = d2.last_parameter();
            let mut ok = true;
            self.param1inf = if d1.has_first_point() {
                d1.first_parameter()
            } else {
                -PRECISION_INFINITE
            };
            self.param1sup = if d1.has_last_point() {
                d1.last_parameter()
            } else {
                PRECISION_INFINITE
            };
            self.param2inf = <PT as PCurveTool<C>>::first_parameter(c2);
            self.param2sup = <PT as PCurveTool<C>>::last_parameter(c2);
            let mut domain_c2_num_inter = Res2dDomain::infinite();

            let mut tab2 = vec![0.0; (nb_inter_c2 + 2) as usize];
            <PT as PCurveTool<C>>::intervals(c2, &mut tab2);

            let mut num_inter_c2 = 1;
            while ok && num_inter_c2 <= nb_inter_c2 {
                (param_inf, param_sup) = <PT as PCurveTool<C>>::get_interval(c2, num_inter_c2 as usize, &tab2);
                if (param_inf > d2_last_param) || (param_sup < d2_first_param) {
                    ok = false;
                } else {
                    if param_inf < d2_first_param {
                        param_inf = d2_first_param;
                    }
                    if param_sup > d2_last_param {
                        param_sup = d2_last_param;
                    }
                    if (param_sup - param_inf) > REAL_EPSILON {
                        domain_c2_num_inter.set_values_bounded(
                            <PT as PCurveTool<C>>::value(c2, param_inf),
                            param_inf,
                            d2.first_tolerance(),
                            <PT as PCurveTool<C>>::value(c2, param_sup),
                            param_sup,
                            d2.last_tolerance(),
                        );
                        self.internal_perform_circle(circ1, d1, c2, &domain_c2_num_inter, tol_conf, tol, true);
                    }
                }
                num_inter_c2 += 1;
            }
        } else {
            self.internal_perform_circle(circ1, d1, c2, d2, tol_conf, tol, false);
        }
    }

    /// OCCT Perform(const gp_Parab2d& Parab1, D1, C2, D2, TolConf, Tol)
    /// (gxx L207-259).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_parabola(
        &mut self,
        parab1: &Parabola2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.reset_fields();
        let nb_inter_c2 = <PT as PCurveTool<C>>::nb_intervals(c2);
        if nb_inter_c2 > 1 {
            let mut param_inf;
            let mut param_sup;
            let d2_first_param = d2.first_parameter();
            let d2_last_param = d2.last_parameter();
            let mut ok = true;
            self.param1inf = if d1.has_first_point() {
                d1.first_parameter()
            } else {
                -PRECISION_INFINITE
            };
            self.param1sup = if d1.has_last_point() {
                d1.last_parameter()
            } else {
                PRECISION_INFINITE
            };
            self.param2inf = <PT as PCurveTool<C>>::first_parameter(c2);
            self.param2sup = <PT as PCurveTool<C>>::last_parameter(c2);
            let mut domain_c2_num_inter = Res2dDomain::infinite();

            let mut tab2 = vec![0.0; (nb_inter_c2 + 2) as usize];
            <PT as PCurveTool<C>>::intervals(c2, &mut tab2);

            let mut num_inter_c2 = 1;
            while ok && num_inter_c2 <= nb_inter_c2 {
                (param_inf, param_sup) = <PT as PCurveTool<C>>::get_interval(c2, num_inter_c2 as usize, &tab2);
                if (param_inf > d2_last_param) || (param_sup < d2_first_param) {
                    ok = false;
                } else {
                    if param_inf < d2_first_param {
                        param_inf = d2_first_param;
                    }
                    if param_sup > d2_last_param {
                        param_sup = d2_last_param;
                    }
                    if (param_sup - param_inf) > REAL_EPSILON {
                        domain_c2_num_inter.set_values_bounded(
                            <PT as PCurveTool<C>>::value(c2, param_inf),
                            param_inf,
                            d2.first_tolerance(),
                            <PT as PCurveTool<C>>::value(c2, param_sup),
                            param_sup,
                            d2.last_tolerance(),
                        );
                        self.internal_perform_parabola(parab1, d1, c2, &domain_c2_num_inter, tol_conf, tol, true);
                    }
                }
                num_inter_c2 += 1;
            }
        } else {
            self.internal_perform_parabola(parab1, d1, c2, d2, tol_conf, tol, false);
        }
    }

    /// OCCT Perform(const gp_Elips2d& Elips1, D1, C2, D2, TolConf, Tol)
    /// (gxx L263-315).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_ellipse(
        &mut self,
        elips1: &Ellipse2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.reset_fields();
        let nb_inter_c2 = <PT as PCurveTool<C>>::nb_intervals(c2);
        if nb_inter_c2 > 1 {
            let mut param_inf;
            let mut param_sup;
            let d2_first_param = d2.first_parameter();
            let d2_last_param = d2.last_parameter();
            let mut ok = true;
            self.param1inf = if d1.has_first_point() {
                d1.first_parameter()
            } else {
                -PRECISION_INFINITE
            };
            self.param1sup = if d1.has_last_point() {
                d1.last_parameter()
            } else {
                PRECISION_INFINITE
            };
            self.param2inf = <PT as PCurveTool<C>>::first_parameter(c2);
            self.param2sup = <PT as PCurveTool<C>>::last_parameter(c2);
            let mut domain_c2_num_inter = Res2dDomain::infinite();

            let mut tab2 = vec![0.0; (nb_inter_c2 + 2) as usize];
            <PT as PCurveTool<C>>::intervals(c2, &mut tab2);

            let mut num_inter_c2 = 1;
            while ok && num_inter_c2 <= nb_inter_c2 {
                (param_inf, param_sup) = <PT as PCurveTool<C>>::get_interval(c2, num_inter_c2 as usize, &tab2);
                if (param_inf > d2_last_param) || (param_sup < d2_first_param) {
                    ok = false;
                } else {
                    if param_inf < d2_first_param {
                        param_inf = d2_first_param;
                    }
                    if param_sup > d2_last_param {
                        param_sup = d2_last_param;
                    }
                    if (param_sup - param_inf) > REAL_EPSILON {
                        domain_c2_num_inter.set_values_bounded(
                            <PT as PCurveTool<C>>::value(c2, param_inf),
                            param_inf,
                            d2.first_tolerance(),
                            <PT as PCurveTool<C>>::value(c2, param_sup),
                            param_sup,
                            d2.last_tolerance(),
                        );
                        self.internal_perform_ellipse(elips1, d1, c2, &domain_c2_num_inter, tol_conf, tol, true);
                    }
                }
                num_inter_c2 += 1;
            }
        } else {
            self.internal_perform_ellipse(elips1, d1, c2, d2, tol_conf, tol, false);
        }
    }

    /// OCCT Perform(const gp_Hypr2d& Hyper1, D1, C2, D2, TolConf, Tol)
    /// (gxx L319-371).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_hyperbola(
        &mut self,
        hyper1: &Hyperbola2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.reset_fields();
        let nb_inter_c2 = <PT as PCurveTool<C>>::nb_intervals(c2);
        if nb_inter_c2 > 1 {
            let mut param_inf;
            let mut param_sup;
            let d2_first_param = d2.first_parameter();
            let d2_last_param = d2.last_parameter();
            let mut ok = true;
            self.param1inf = if d1.has_first_point() {
                d1.first_parameter()
            } else {
                -PRECISION_INFINITE
            };
            self.param1sup = if d1.has_last_point() {
                d1.last_parameter()
            } else {
                PRECISION_INFINITE
            };
            self.param2inf = <PT as PCurveTool<C>>::first_parameter(c2);
            self.param2sup = <PT as PCurveTool<C>>::last_parameter(c2);
            let mut domain_c2_num_inter = Res2dDomain::infinite();

            let mut tab2 = vec![0.0; (nb_inter_c2 + 2) as usize];
            <PT as PCurveTool<C>>::intervals(c2, &mut tab2);

            let mut num_inter_c2 = 1;
            while ok && num_inter_c2 <= nb_inter_c2 {
                (param_inf, param_sup) = <PT as PCurveTool<C>>::get_interval(c2, num_inter_c2 as usize, &tab2);
                if (param_inf > d2_last_param) || (param_sup < d2_first_param) {
                    ok = false;
                } else {
                    if param_inf < d2_first_param {
                        param_inf = d2_first_param;
                    }
                    if param_sup > d2_last_param {
                        param_sup = d2_last_param;
                    }
                    if (param_sup - param_inf) > REAL_EPSILON {
                        domain_c2_num_inter.set_values_bounded(
                            <PT as PCurveTool<C>>::value(c2, param_inf),
                            param_inf,
                            d2.first_tolerance(),
                            <PT as PCurveTool<C>>::value(c2, param_sup),
                            param_sup,
                            d2.last_tolerance(),
                        );
                        self.internal_perform_hyperbola(hyper1, d1, c2, &domain_c2_num_inter, tol_conf, tol, true);
                    }
                }
                num_inter_c2 += 1;
            }
        } else {
            self.internal_perform_hyperbola(hyper1, d1, c2, d2, tol_conf, tol, false);
        }
    }

    /// OCCT InternalPerform(const gp_Lin2d& Lin1, ...) (gxx L381-479).
    /// Si Composite == True les resultats sont Ajoutes, sinon ils sont
    /// Copies.
    #[allow(clippy::too_many_arguments)]
    fn internal_perform_line(
        &mut self,
        lin1: &Line2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        composite: bool,
    ) {
        let typ2 = <PT as PCurveTool<C>>::get_type(c2);

        // Shared tail of every InternalPerform arm (gxx L398-406 and
        // repeats): Append in composite mode, SetValues otherwise.
        macro_rules! finish {
            ($inter:expr) => {
                if composite {
                    self.base.append_intersector(
                        &$inter.base,
                        self.param1inf,
                        self.param1sup,
                        self.param2inf,
                        self.param2sup,
                    );
                } else {
                    self.base.set_values(&$inter.base);
                }
            };
        }

        match typ2 {
            Curve2dType::Line => {
                self.intconiconi.base.set_reversed_parameters(false);
                let l2 = <PT as PCurveTool<C>>::line(c2);
                self.intconiconi.perform_line_line(lin1, d1, &l2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Circle => {
                self.intconiconi.base.set_reversed_parameters(false);
                let c2c = <PT as PCurveTool<C>>::circle(c2);
                self.intconiconi.perform_line_circle(lin1, d1, &c2c, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Ellipse => {
                self.intconiconi.base.set_reversed_parameters(false);
                let e2 = <PT as PCurveTool<C>>::ellipse(c2);
                self.intconiconi.perform_line_ellipse(lin1, d1, &e2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Parabola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let p2 = <PT as PCurveTool<C>>::parabola(c2);
                self.intconiconi.perform_line_parabola(lin1, d1, &p2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Hyperbola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let h2 = <PT as PCurveTool<C>>::hyperbola(c2);
                self.intconiconi.perform_line_hyperbola(lin1, d1, &h2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            _ => {
                self.intconicurv.base.set_reversed_parameters(false);
                self.intconicurv.perform_line(lin1, d1, c2, d2, tol_conf, tol);
                finish!(self.intconicurv);
            }
        }
    }

    /// OCCT InternalPerform(const gp_Circ2d& Circ1, ...) (gxx L483-581).
    #[allow(clippy::too_many_arguments)]
    fn internal_perform_circle(
        &mut self,
        circ1: &Circle2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        composite: bool,
    ) {
        let typ2 = <PT as PCurveTool<C>>::get_type(c2);

        macro_rules! finish {
            ($inter:expr) => {
                if composite {
                    self.base.append_intersector(
                        &$inter.base,
                        self.param1inf,
                        self.param1sup,
                        self.param2inf,
                        self.param2sup,
                    );
                } else {
                    self.base.set_values(&$inter.base);
                }
            };
        }

        match typ2 {
            Curve2dType::Line => {
                self.intconiconi.base.set_reversed_parameters(true);
                let l2 = <PT as PCurveTool<C>>::line(c2);
                self.intconiconi.perform_line_circle(&l2, d2, circ1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Circle => {
                self.intconiconi.base.set_reversed_parameters(false);
                let c2c = <PT as PCurveTool<C>>::circle(c2);
                self.intconiconi.perform_circle_circle(circ1, d1, &c2c, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Ellipse => {
                self.intconiconi.base.set_reversed_parameters(false);
                let e2 = <PT as PCurveTool<C>>::ellipse(c2);
                self.intconiconi.perform_circle_ellipse(circ1, d1, &e2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Parabola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let p2 = <PT as PCurveTool<C>>::parabola(c2);
                self.intconiconi.perform_circle_parabola(circ1, d1, &p2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Hyperbola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let h2 = <PT as PCurveTool<C>>::hyperbola(c2);
                self.intconiconi.perform_circle_hyperbola(circ1, d1, &h2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            _ => {
                self.intconicurv.base.set_reversed_parameters(false);
                self.intconicurv.perform_circle(circ1, d1, c2, d2, tol_conf, tol);
                finish!(self.intconicurv);
            }
        }
    }

    /// OCCT InternalPerform(const gp_Elips2d& Elips1, ...) (gxx L585-683).
    #[allow(clippy::too_many_arguments)]
    fn internal_perform_ellipse(
        &mut self,
        elips1: &Ellipse2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        composite: bool,
    ) {
        let typ2 = <PT as PCurveTool<C>>::get_type(c2);

        macro_rules! finish {
            ($inter:expr) => {
                if composite {
                    self.base.append_intersector(
                        &$inter.base,
                        self.param1inf,
                        self.param1sup,
                        self.param2inf,
                        self.param2sup,
                    );
                } else {
                    self.base.set_values(&$inter.base);
                }
            };
        }

        match typ2 {
            Curve2dType::Line => {
                self.intconiconi.base.set_reversed_parameters(true);
                let l2 = <PT as PCurveTool<C>>::line(c2);
                self.intconiconi.perform_line_ellipse(&l2, d2, elips1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Circle => {
                self.intconiconi.base.set_reversed_parameters(true);
                let c2c = <PT as PCurveTool<C>>::circle(c2);
                self.intconiconi.perform_circle_ellipse(&c2c, d2, elips1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Ellipse => {
                self.intconiconi.base.set_reversed_parameters(false);
                let e2 = <PT as PCurveTool<C>>::ellipse(c2);
                self.intconiconi.perform_ellipse_ellipse(elips1, d1, &e2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Parabola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let p2 = <PT as PCurveTool<C>>::parabola(c2);
                self.intconiconi.perform_ellipse_parabola(elips1, d1, &p2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Hyperbola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let h2 = <PT as PCurveTool<C>>::hyperbola(c2);
                self.intconiconi.perform_ellipse_hyperbola(elips1, d1, &h2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            _ => {
                self.intconicurv.base.set_reversed_parameters(false);
                self.intconicurv.perform_ellipse(elips1, d1, c2, d2, tol_conf, tol);
                finish!(self.intconicurv);
            }
        }
    }

    /// OCCT InternalPerform(const gp_Parab2d& Parab1, ...) (gxx L687-785).
    #[allow(clippy::too_many_arguments)]
    fn internal_perform_parabola(
        &mut self,
        parab1: &Parabola2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        composite: bool,
    ) {
        let typ2 = <PT as PCurveTool<C>>::get_type(c2);

        macro_rules! finish {
            ($inter:expr) => {
                if composite {
                    self.base.append_intersector(
                        &$inter.base,
                        self.param1inf,
                        self.param1sup,
                        self.param2inf,
                        self.param2sup,
                    );
                } else {
                    self.base.set_values(&$inter.base);
                }
            };
        }

        match typ2 {
            Curve2dType::Line => {
                self.intconiconi.base.set_reversed_parameters(true);
                let l2 = <PT as PCurveTool<C>>::line(c2);
                self.intconiconi.perform_line_parabola(&l2, d2, parab1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Circle => {
                self.intconiconi.base.set_reversed_parameters(true);
                let c2c = <PT as PCurveTool<C>>::circle(c2);
                self.intconiconi.perform_circle_parabola(&c2c, d2, parab1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Ellipse => {
                self.intconiconi.base.set_reversed_parameters(true);
                let e2 = <PT as PCurveTool<C>>::ellipse(c2);
                self.intconiconi.perform_ellipse_parabola(&e2, d2, parab1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Parabola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let p2 = <PT as PCurveTool<C>>::parabola(c2);
                self.intconiconi.perform_parabola_parabola(parab1, d1, &p2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Hyperbola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let h2 = <PT as PCurveTool<C>>::hyperbola(c2);
                self.intconiconi.perform_parabola_hyperbola(parab1, d1, &h2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            _ => {
                self.intconicurv.base.set_reversed_parameters(false);
                self.intconicurv.perform_parabola(parab1, d1, c2, d2, tol_conf, tol);
                finish!(self.intconicurv);
            }
        }
    }

    /// OCCT InternalPerform(const gp_Hypr2d& Hyper1, ...) (gxx L789-887).
    #[allow(clippy::too_many_arguments)]
    fn internal_perform_hyperbola(
        &mut self,
        hyper1: &Hyperbola2d,
        d1: &Res2dDomain,
        c2: &C,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        composite: bool,
    ) {
        let typ2 = <PT as PCurveTool<C>>::get_type(c2);

        macro_rules! finish {
            ($inter:expr) => {
                if composite {
                    self.base.append_intersector(
                        &$inter.base,
                        self.param1inf,
                        self.param1sup,
                        self.param2inf,
                        self.param2sup,
                    );
                } else {
                    self.base.set_values(&$inter.base);
                }
            };
        }

        match typ2 {
            Curve2dType::Line => {
                self.intconiconi.base.set_reversed_parameters(true);
                let l2 = <PT as PCurveTool<C>>::line(c2);
                self.intconiconi.perform_line_hyperbola(&l2, d2, hyper1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Circle => {
                self.intconiconi.base.set_reversed_parameters(true);
                let c2c = <PT as PCurveTool<C>>::circle(c2);
                self.intconiconi.perform_circle_hyperbola(&c2c, d2, hyper1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Ellipse => {
                self.intconiconi.base.set_reversed_parameters(true);
                let e2 = <PT as PCurveTool<C>>::ellipse(c2);
                self.intconiconi.perform_ellipse_hyperbola(&e2, d2, hyper1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Parabola => {
                self.intconiconi.base.set_reversed_parameters(true);
                let p2 = <PT as PCurveTool<C>>::parabola(c2);
                self.intconiconi.perform_parabola_hyperbola(&p2, d2, hyper1, d1, tol_conf, tol);
                finish!(self.intconiconi);
            }
            Curve2dType::Hyperbola => {
                self.intconiconi.base.set_reversed_parameters(false);
                let h2 = <PT as PCurveTool<C>>::hyperbola(c2);
                self.intconiconi.perform_hyperbola_hyperbola(hyper1, d1, &h2, d2, tol_conf, tol);
                finish!(self.intconiconi);
            }
            _ => {
                self.intconicurv.base.set_reversed_parameters(false);
                self.intconicurv.perform_hyperbola(hyper1, d1, c2, d2, tol_conf, tol);
                finish!(self.intconicurv);
            }
        }
    }
}

impl<C: ?Sized, PT: ParTool<C> + PCurveTool<C>, JT: ProjectOnPCurveTool<C>> Default
    for UserIntConicCurveGen<C, PT, JT>
{
    fn default() -> Self {
        UserIntConicCurveGen::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::geom2d_int::{Curve2dAdaptor, Geom2dCurveTool, TheProjPCurOfGInter};
    use rcad_kernel::geom::{BezierCurve2, Circle2d, Curve2d, Line2d};

    const TAU: f64 = std::f64::consts::TAU;
    const TOL: f64 = 1.0e-6;

    /// Unit circle at the origin, full closed domain [0, 2pi].
    fn unit_circle() -> (Circle2d, Res2dDomain) {
        let c = Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::X,
            y_dir: DVec2::Y,
            radius: 1.0,
        };
        let mut d = Res2dDomain::bounded(DVec2::X, 0.0, TOL, DVec2::X, TAU, TOL);
        d.set_equivalent_parameters(0.0, TAU);
        (c, d)
    }

    /// A "straight" Bezier: colinear poles P0=(-2,0.5), P1=(0,0), P2=(2,-0.5)
    /// give the affine point map P(u) = (-2+4u, 0.5-u) on [0,1].
    fn straight_bezier() -> (Curve2d, Res2dDomain) {
        let bez = BezierCurve2 {
            control_points: vec![DVec2::new(-2.0, 0.5), DVec2::ZERO, DVec2::new(2.0, -0.5)],
            weights: vec![1.0, 1.0, 1.0],
        };
        let c = Curve2d::Bezier(bez);
        let d = Res2dDomain::bounded(c.value(0.0), 0.0, TOL, c.value(1.0), 1.0, TOL);
        (c, d)
    }

    fn make_user<'a>(
    ) -> UserIntConicCurveGen<dyn Curve2dAdaptor + 'a, Geom2dCurveTool, TheProjPCurOfGInter> {
        UserIntConicCurveGen::new()
    }

    /// Exact crossings of the unit circle with P(u) = (-2+4u, 0.5-u):
    /// 17u^2 - 17u + 3.25 = 0 -> u = (17 +- 2*sqrt(17))/34.
    fn expected_crossings() -> [(DVec2, f64, f64); 2] {
        let sq = 2.0 * 17.0f64.sqrt();
        let mut out = [(DVec2::ZERO, 0.0, 0.0); 2];
        for (k, u) in [(17.0 - sq) / 34.0, (17.0 + sq) / 34.0].into_iter().enumerate() {
            let p = DVec2::new(-2.0 + 4.0 * u, 0.5 - u);
            let theta = p.y.atan2(p.x);
            let theta = if theta < 0.0 { theta + TAU } else { theta };
            out[k] = (p, u, theta);
        }
        out
    }

    fn assert_crossings(base: &IntersectionBase) {
        assert!(base.is_done());
        assert_eq!(base.nb_points(), 2, "two transversal crossings expected");
        assert_eq!(base.nb_segments(), 0);
        let exp = expected_crossings();
        let mut seen = [false; 2];
        for i in 1..=base.nb_points() {
            let pt = base.point(i);
            let v = pt.value();
            assert!((v.length() - 1.0).abs() < 1.0e-7, "point on circle: {v:?}");
            let mut matched = false;
            for k in 0..2 {
                let (p, u, theta) = exp[k];
                if v.distance(p) < 1.0e-6 {
                    matched = true;
                    seen[k] = true;
                    assert!(
                        (pt.param_on_second() - u).abs() < 1.0e-6,
                        "pcurve param: {} vs {u}",
                        pt.param_on_second()
                    );
                    let got = pt.param_on_first();
                    let got = if got < 0.0 { got + TAU } else { got };
                    assert!((got - theta).abs() < 1.0e-6, "circle param: {got} vs {theta}");
                }
            }
            assert!(matched, "unexpected crossing point {v:?}");
        }
        assert!(seen[0] && seen[1]);
    }

    /// OCCT Geom2dInt_TheIntConicCurveOfGInter(C, D1, PCurve, D2, TolConf,
    /// Tol) (IntCurve_IntConicCurveGen.gxx L24-42) driving the
    /// IntImpParGen_Intersector engine through a non-conic pcurve.
    #[test]
    fn shell_circle_x_straight_bezier() {
        let (circle, d1) = unit_circle();
        let (c2, d2) = straight_bezier();
        let intp = crate::geomalgo::geom2d_int::TheIntConicCurveOfGInter::new_circle(
            &circle,
            &d1,
            &c2,
            &d2,
            TOL,
            TOL,
        );
        assert_crossings(&intp.base);
    }

    /// OCCT IntCurve_UserIntConicCurveGen(Circ1, D1, C2, D2, TolConf, Tol)
    /// (gxx L44-52) -> Perform(Circ1, ...) single-interval path ->
    /// InternalPerform default arm -> intconicurv (the IntConicCurveGen
    /// shell).
    #[test]
    fn user_circle_x_straight_bezier() {
        let (circle, d1) = unit_circle();
        let (c2, d2) = straight_bezier();
        let mut intp = make_user();
        intp.perform_circle(&circle, &d1, &c2, &d2, TOL, TOL);
        assert_crossings(&intp.base);
    }

    /// The IntConicCurveGen circle-ctor/Perform quirk (lxx L63-69): a
    /// non-closed D1 gains SetEquivalentParameters(F, F+2pi).  The domain
    /// bounds themselves stay untouched: an arc domain away from the
    /// crossings still yields no results, while a near-full domain keeps
    /// both.
    #[test]
    fn circle_ctor_quirk_non_closed_domain() {
        let (circle, _) = unit_circle();
        let (c2, d2) = straight_bezier();
        let arc = |a: f64, b: f64| {
            Res2dDomain::bounded(
                DVec2::new(a.cos(), a.sin()),
                a,
                TOL,
                DVec2::new(b.cos(), b.sin()),
                b,
                TOL,
            )
        };
        // Arc [0.1, 1.5]: both crossing angles (2.897, 6.038) are outside
        // the bounds and get dropped (OCCT behaves the same).
        let intp = crate::geomalgo::geom2d_int::TheIntConicCurveOfGInter::new_circle(
            &circle,
            &arc(0.1, 1.5),
            &c2,
            &d2,
            TOL,
            TOL,
        );
        assert!(intp.base.is_done());
        assert_eq!(intp.base.nb_points(), 0);
        assert_eq!(intp.base.nb_segments(), 0);
        // Near-full domain [0.1, 6.3]: both crossings survive.
        let intp = crate::geomalgo::geom2d_int::TheIntConicCurveOfGInter::new_circle(
            &circle,
            &arc(0.1, 6.3),
            &c2,
            &d2,
            TOL,
            TOL,
        );
        assert_crossings(&intp.base);
    }

    /// OCCT UserIntConicCurveGen(Lin1, ...) default arm with an implicit
    /// line: x = -0.3 meets P(u) = (-2+4u, 0.5-u) at u = 0.425 only.
    #[test]
    fn user_line_x_straight_bezier() {
        let line = Line2d::new(DVec2::new(-0.3, 0.0), DVec2::Y);
        let d1 = Res2dDomain::infinite();
        let (c2, d2) = straight_bezier();
        let mut intp = make_user();
        intp.perform_line(&line, &d1, &c2, &d2, TOL, TOL);
        assert!(intp.base.is_done());
        assert_eq!(intp.base.nb_points(), 1, "single crossing expected");
        let pt = intp.base.point(1);
        let v = pt.value();
        assert!((v.x - -0.3).abs() < 1.0e-7 && (v.y - 0.075).abs() < 1.0e-7, "got {v:?}");
        assert!((pt.param_on_second() - 0.425).abs() < 1.0e-7);
        assert!((pt.param_on_first() - 0.075).abs() < 1.0e-7);
    }

    /// IntImpParGen statics anchors (IntImpParGen.cxx).
    #[test]
    fn imp_par_gen_statics() {
        // NormalizeOnDomain on the closed circle domain.
        let (_, d1) = unit_circle();
        let got = crate::geomalgo::int_imp_par_gen::normalize_on_domain(-0.5, &d1);
        assert!((got - (TAU - 0.5)).abs() < 1.0e-12);
        let got = crate::geomalgo::int_imp_par_gen::normalize_on_domain(TAU + 0.25, &d1);
        assert!((got - 0.25).abs() < 1.0e-12);
        // DeterminePosition: middle / head / end.
        use crate::geomalgo::int_res2d::Position;
        let mut pos = Position::Middle;
        crate::geomalgo::int_imp_par_gen::determine_position(&mut pos, &d1, DVec2::new(5.0, 5.0), 1.0);
        assert_eq!(pos, Position::Middle);
        crate::geomalgo::int_imp_par_gen::determine_position(&mut pos, &d1, DVec2::X, 0.0);
        assert_eq!(pos, Position::Head);
        crate::geomalgo::int_imp_par_gen::determine_position(&mut pos, &d1, DVec2::X, TAU);
        assert_eq!(pos, Position::End);
    }


}
