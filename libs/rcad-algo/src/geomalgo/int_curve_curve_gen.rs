// OCCT Geom2dInt_GInter == IntCurve_IntCurveCurveGen — 1:1 Rust translation
// of IntCurve_IntCurveCurveGen.hxx/.lxx/.gxx as instantiated by
// Geom2dInt_GInter_0.cxx (TheCurve=Adaptor2d_Curve2d,
// TheCurveTool=Geom2dInt_Geom2dCurveTool, IntCurve_TheIntConicCurve=
// Geom2dInt_TheIntConicCurveOfGInter, IntCurve_IntConicConic,
// IntCurve_TheIntPCurvePCurve=Geom2dInt_TheIntPCurvePCurveOfGInter).
//
// Members:
//   - intconiconi: IntConicConic (conic x conic)              — ported
//     (Ellipse-Ellipse overload; see int_conic_conic.rs).
//   - intconicurv: TheIntConicCurveOfGInter (conic x pcurve)  — ported
//     (geom2d_int.rs).
//   - intcurvcurv: TheIntPCurvePCurveOfGInter (pcurve x pcurve,
//     IntCurve_IntPolyPolyGen.gxx, 1797 lines)                — NOT yet
//     ported; the dispatch arms that would call it carry an
//     unimplemented!() (ThruSections skeleton precedent).

use super::geom2d_int::{
    geom2d_curve_tool, TheIntConicCurveOfGInter, Curve2dAdaptor, Curve2dType,
};
use super::int_conic_conic::IntConicConic;
use super::int_res2d::{Domain as Res2dDomain, IntersectionBase};

/// OCCT Precision::Infinite() (Precision.hxx — 1e100).
const PRECISION_INFINITE: f64 = 1.0e100;

/// OCCT Precision::IsInfinite(V).
fn precision_is_infinite(v: f64) -> bool {
    v >= PRECISION_INFINITE || v <= -PRECISION_INFINITE
}

/// OCCT IntCurve_IntCurveCurveGen (Geom2dInt_GInter).
#[derive(Debug, Clone)]
pub struct IntCurveCurveGen {
    pub base: IntersectionBase,
    /// OCCT intconiconi.
    intconiconi: IntConicConic,
    /// OCCT intconicurv.
    intconicurv: TheIntConicCurveOfGInter,
    /// OCCT param1inf / param1sup / param2inf / param2sup.
    param1inf: f64,
    param1sup: f64,
    param2inf: f64,
    param2sup: f64,
}

impl IntCurveCurveGen {
    /// OCCT IntCurve_IntCurveCurveGen() (lxx L20-26).
    pub fn new() -> Self {
        IntCurveCurveGen {
            base: IntersectionBase::new(),
            intconiconi: IntConicConic::new(),
            intconicurv: TheIntConicCurveOfGInter::bare(),
            param1inf: -PRECISION_INFINITE,
            param1sup: PRECISION_INFINITE,
            param2inf: -PRECISION_INFINITE,
            param2sup: PRECISION_INFINITE,
        }
    }

    /// OCCT IntCurve_IntCurveCurveGen(C, TolConf, Tol) (lxx L29-38).
    pub fn new_c(c: &dyn Curve2dAdaptor, tol_conf: f64, tol: f64) -> Self {
        let mut r = IntCurveCurveGen::new();
        r.perform_c(c, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntCurveCurveGen(C, D, TolConf, Tol) (lxx L41-51).
    pub fn new_cd(c: &dyn Curve2dAdaptor, d: &Res2dDomain, tol_conf: f64, tol: f64) -> Self {
        let mut r = IntCurveCurveGen::new();
        r.perform_cd(c, d, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntCurveCurveGen(C1, C2, TolConf, Tol) (lxx L54-64).
    pub fn new_cc(
        c1: &dyn Curve2dAdaptor,
        c2: &dyn Curve2dAdaptor,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntCurveCurveGen::new();
        r.perform_cc(c1, c2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntCurveCurveGen(C1, D1, C2, TolConf, Tol) (lxx L67-78).
    pub fn new_cd_c(
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        c2: &dyn Curve2dAdaptor,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntCurveCurveGen::new();
        r.perform_cd_c(c1, d1, c2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntCurveCurveGen(C1, C2, D2, TolConf, Tol) (lxx L81-92).
    pub fn new_cc_d(
        c1: &dyn Curve2dAdaptor,
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntCurveCurveGen::new();
        r.perform_cc_d(c1, c2, d2, tol_conf, tol);
        r
    }

    /// OCCT IntCurve_IntCurveCurveGen(C1, D1, C2, D2, TolConf, Tol)
    /// (lxx L95-107).
    #[allow(clippy::too_many_arguments)]
    pub fn new_cd_cd(
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) -> Self {
        let mut r = IntCurveCurveGen::new();
        r.perform_cd_cd(c1, d1, c2, d2, tol_conf, tol);
        r
    }

    // -- IntRes2d_Intersection result accessors (base delegates) ------------

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.base.is_done()
    }

    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.base.nb_points()
    }

    /// OCCT Point(N) — 1-based.
    pub fn point(&self, n: usize) -> &super::int_res2d::IntersectionPoint {
        self.base.point(n)
    }

    /// OCCT NbSegments().
    pub fn nb_segments(&self) -> usize {
        self.base.nb_segments()
    }

    /// OCCT Segment(N) — 1-based.
    pub fn segment(&self, n: usize) -> &super::int_res2d::IntersectionSegment {
        self.base.segment(n)
    }

    // -- Perform overloads --------------------------------------------------

    /// OCCT Perform(C, TolConf, Tol) — self-intersection
    /// (IntCurve_IntCurveCurveGen.gxx L29-89).
    pub fn perform_c(&mut self, c: &dyn Curve2dAdaptor, tol_conf: f64, tol: f64) {
        let mut d1 = Res2dDomain::infinite();
        let tol_domain = if tol < tol_conf { tol_conf } else { tol };
        let typ = geom2d_curve_tool::get_type(c);
        match typ {
            Curve2dType::Ellipse
            | Curve2dType::Circle
            | Curve2dType::Parabola
            | Curve2dType::Hyperbola
            | Curve2dType::Line => {
                self.base.reset_fields();
                self.base.done = true;
                return;
            }
            _ => {
                let paraminf = geom2d_curve_tool::first_parameter(c);
                let paramsup = geom2d_curve_tool::last_parameter(c);
                if precision_is_infinite(paraminf) && precision_is_infinite(paramsup) {
                    self.base.done = false;
                    return;
                }
                //
                if paraminf > -PRECISION_INFINITE {
                    if paramsup < PRECISION_INFINITE {
                        //--         paraminf-----------paramsup
                        d1.set_values_bounded(
                            geom2d_curve_tool::value(c, paraminf),
                            paraminf,
                            tol_domain,
                            geom2d_curve_tool::value(c, paramsup),
                            paramsup,
                            tol_domain,
                        );
                    } else {
                        //--        paraminf------------...
                        d1.set_values_semi(
                            geom2d_curve_tool::value(c, paraminf),
                            paraminf,
                            tol_domain,
                            true,
                        );
                    }
                } else if paramsup < PRECISION_INFINITE {
                    //--    ...-----------------paramsup
                    d1.set_values_semi(
                        geom2d_curve_tool::value(c, paramsup),
                        paramsup,
                        tol_domain,
                        false,
                    );
                }
                self.base.reset_fields();
                self.perform_cd(c, &d1, tol_conf, tol);
            }
        }
    }

    /// OCCT Perform(C, D, TolConf, Tol) — self-intersection with domain
    /// (IntCurve_IntCurveCurveGen.gxx L91-116).
    pub fn perform_cd(&mut self, c: &dyn Curve2dAdaptor, d: &Res2dDomain, tol_conf: f64, tol: f64) {
        let typ = geom2d_curve_tool::get_type(c);
        match typ {
            Curve2dType::Ellipse
            | Curve2dType::Circle
            | Curve2dType::Parabola
            | Curve2dType::Hyperbola
            | Curve2dType::Line => {
                self.base.reset_fields();
                self.base.done = true;
                return;
            }
            _ => {
                self.base.reset_fields();
                // intcurvcurv.SetReversedParameters(false);
                // intcurvcurv.Perform(C, D1, TolConf, Tol);
                // -- Geom2dInt_TheIntPCurvePCurveOfGInter (IntCurve_IntPolyPolyGen)
                //    is not yet ported (see file header).
                unimplemented!(
                    "IntCurve_IntCurveCurveGen: pcurve x same-pcurve self-intersection \
                     (Geom2dInt_TheIntPCurvePCurveOfGInter / IntCurve_IntPolyPolyGen) is not \
                     ported yet"
                );
            }
        }
    }

    /// OCCT Perform(C1, C2, TolConf, Tol) (IntCurve_IntCurveCurveGen.lxx L109-121).
    pub fn perform_cc(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        c2: &dyn Curve2dAdaptor,
        tol_conf: f64,
        tol: f64,
    ) {
        let tol_domain = tol;
        let tol_domain = if tol_conf > tol_domain { tol_conf } else { tol_domain };
        let d1 = self.compute_domain(c1, tol_domain);
        let d2 = self.compute_domain(c2, tol_domain);
        self.perform_cd_cd(c1, &d1, c2, &d2, tol_conf, tol);
    }

    /// OCCT Perform(C1, D1, C2, TolConf, Tol) (IntCurve_IntCurveCurveGen.lxx L124-136).
    pub fn perform_cd_c(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        c2: &dyn Curve2dAdaptor,
        tol_conf: f64,
        tol: f64,
    ) {
        let tol_domain = tol;
        let tol_domain = if tol_conf > tol_domain { tol_conf } else { tol_domain };
        let d2 = self.compute_domain(c2, tol_domain);
        self.perform_cd_cd(c1, d1, c2, &d2, tol_conf, tol);
    }

    /// OCCT Perform(C1, C2, D2, TolConf, Tol) (IntCurve_IntCurveCurveGen.lxx L139-151).
    pub fn perform_cc_d(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        let tol_domain = tol;
        let tol_domain = if tol_conf > tol_domain { tol_conf } else { tol_domain };
        let d1 = self.compute_domain(c1, tol_domain);
        self.perform_cd_cd(c1, &d1, c2, d2, tol_conf, tol);
    }

    /// OCCT ComputeDomain(C1, TolDomain) (IntCurve_IntCurveCurveGen.gxx L120-177).
    pub fn compute_domain(&self, c1: &dyn Curve2dAdaptor, tol_domain: f64) -> Res2dDomain {
        let mut d1 = Res2dDomain::infinite();

        let typ = geom2d_curve_tool::get_type(c1);
        match typ {
            Curve2dType::Ellipse | Curve2dType::Circle => {
                //---------------------------------------------------------------
                //-- if the curve is a trimmed curve, first and last parameters
                //-- will be the parameters used to build the domain
                //--
                let firstparameter = geom2d_curve_tool::first_parameter(c1);
                let lastparameter = geom2d_curve_tool::last_parameter(c1);

                let p1 = geom2d_curve_tool::value(c1, firstparameter);
                let p2 = geom2d_curve_tool::value(c1, lastparameter);
                d1.set_values_bounded(p1, firstparameter, tol_domain, p2, lastparameter, tol_domain);
                d1.set_equivalent_parameters(firstparameter, firstparameter + std::f64::consts::PI + std::f64::consts::PI);
            }
            _ => {
                let paraminf = geom2d_curve_tool::first_parameter(c1);
                let paramsup = geom2d_curve_tool::last_parameter(c1);
                if paraminf > -PRECISION_INFINITE {
                    if paramsup < PRECISION_INFINITE {
                        //--         paraminf-----------paramsup
                        d1.set_values_bounded(
                            geom2d_curve_tool::value(c1, paraminf),
                            paraminf,
                            tol_domain,
                            geom2d_curve_tool::value(c1, paramsup),
                            paramsup,
                            tol_domain,
                        );
                    } else {
                        //--        paraminf------------...
                        d1.set_values_semi(
                            geom2d_curve_tool::value(c1, paraminf),
                            paraminf,
                            tol_domain,
                            true,
                        );
                    }
                } else if paramsup < PRECISION_INFINITE {
                    //--    ...-----------------paramsup
                    d1.set_values_semi(
                        geom2d_curve_tool::value(c1, paramsup),
                        paramsup,
                        tol_domain,
                        false,
                    );
                }
            }
        }
        d1
    }

    /// OCCT Perform(C1, D1, C2, D2, TolConf, Tol)
    /// (IntCurve_IntCurveCurveGen.gxx L182-225).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_cd_cd(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        self.base.reset_fields();
        let nbi1 = geom2d_curve_tool::nb_intervals(c1);
        if nbi1 > 1 {
            self.param1inf = geom2d_curve_tool::first_parameter(c1);
            self.param1sup = geom2d_curve_tool::last_parameter(c1);
        } else {
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
        }
        let nbi2 = geom2d_curve_tool::nb_intervals(c2);
        if nbi2 > 1 {
            self.param2inf = geom2d_curve_tool::first_parameter(c2);
            self.param2sup = geom2d_curve_tool::last_parameter(c2);
        } else {
            self.param2inf = if d2.has_first_point() {
                d2.first_parameter()
            } else {
                -PRECISION_INFINITE
            };
            self.param2sup = if d2.has_last_point() {
                d2.last_parameter()
            } else {
                PRECISION_INFINITE
            };
        }
        if nbi1 > 1 || nbi2 > 1 {
            // NCollection_Array1<double> Tab(1, nbi+1) — 0-based storage of a
            // 1-based array of length nbi+2 (logical indices 1..=nbi+1).
            let mut tab1 = vec![0.0; (nbi1 + 2) as usize];
            let mut tab2 = vec![0.0; (nbi2 + 2) as usize];
            geom2d_curve_tool::intervals(c1, &mut tab1);
            geom2d_curve_tool::intervals(c2, &mut tab2);
            self.internal_composite_perform(
                c1, d1, 1, nbi1, &tab1, c2, d2, 1, nbi2, &tab2, tol_conf, tol, true,
            );
            return;
        } else {
            self.internal_perform(c1, d1, c2, d2, tol_conf, tol, false);
        }
    }

    /// OCCT InternalPerform — Suppose des Courbes Lin...Other
    /// (IntCurve_IntCurveCurveGen.gxx L235-816). Si Composite == True les
    /// resultats sont Ajoutes, sinon ils sont Copies.
    #[allow(clippy::too_many_arguments)]
    fn internal_perform(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
        composite: bool,
    ) {
        let typ1 = geom2d_curve_tool::get_type(c1);
        let typ2 = geom2d_curve_tool::get_type(c2);

        // The per-arm tail shared by every switch arm: Append in composite
        // mode, SetValues otherwise (gxx L255-262 and all repeats).
        macro_rules! finish_coniconiconi {
            () => {
                if composite {
                    self.base.append_intersector(
                        &self.intconiconi.base,
                        self.param1inf,
                        self.param1sup,
                        self.param2inf,
                        self.param2sup,
                    );
                } else {
                    self.base.set_values(&self.intconiconi.base);
                }
            };
        }
        macro_rules! finish_conicurv {
            () => {
                if composite {
                    self.base.append_intersector(
                        &self.intconicurv.base,
                        self.param1inf,
                        self.param1sup,
                        self.param2inf,
                        self.param2sup,
                    );
                } else {
                    self.base.set_values(&self.intconicurv.base);
                }
            };
        }

        match typ1 {
            Curve2dType::Line => match typ2 {
                Curve2dType::Line => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_line_line(
                        &geom2d_curve_tool::line(c1),
                        d1,
                        &geom2d_curve_tool::line(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Circle => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_line_circle(
                        &geom2d_curve_tool::line(c1),
                        d1,
                        &geom2d_curve_tool::circle(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Ellipse => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_line_ellipse(
                        &geom2d_curve_tool::line(c1),
                        d1,
                        &geom2d_curve_tool::ellipse(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Parabola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_line_parabola(
                        &geom2d_curve_tool::line(c1),
                        d1,
                        &geom2d_curve_tool::parabola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Hyperbola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_line_hyperbola(
                        &geom2d_curve_tool::line(c1),
                        d1,
                        &geom2d_curve_tool::hyperbola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                _ => {
                    self.intconicurv.base.set_reversed_parameters(false);
                    self.intconicurv.perform_line(
                        &geom2d_curve_tool::line(c1),
                        d1,
                        c2,
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
            },
            Curve2dType::Circle => match typ2 {
                Curve2dType::Line => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_line_circle(
                        &geom2d_curve_tool::line(c2),
                        d2,
                        &geom2d_curve_tool::circle(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Circle => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_circle_circle(
                        &geom2d_curve_tool::circle(c1),
                        d1,
                        &geom2d_curve_tool::circle(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Ellipse => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_circle_ellipse(
                        &geom2d_curve_tool::circle(c1),
                        d1,
                        &geom2d_curve_tool::ellipse(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Parabola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_circle_parabola(
                        &geom2d_curve_tool::circle(c1),
                        d1,
                        &geom2d_curve_tool::parabola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Hyperbola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_circle_hyperbola(
                        &geom2d_curve_tool::circle(c1),
                        d1,
                        &geom2d_curve_tool::hyperbola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                _ => {
                    self.intconicurv.base.set_reversed_parameters(false);
                    self.intconicurv.perform_circle(
                        &geom2d_curve_tool::circle(c1),
                        d1,
                        c2,
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
            },
            Curve2dType::Ellipse => match typ2 {
                Curve2dType::Line => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_line_ellipse(
                        &geom2d_curve_tool::line(c2),
                        d2,
                        &geom2d_curve_tool::ellipse(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Circle => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_circle_ellipse(
                        &geom2d_curve_tool::circle(c2),
                        d2,
                        &geom2d_curve_tool::ellipse(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Ellipse => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_ellipse_ellipse(
                        &geom2d_curve_tool::ellipse(c1),
                        d1,
                        &geom2d_curve_tool::ellipse(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Parabola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_ellipse_parabola(
                        &geom2d_curve_tool::ellipse(c1),
                        d1,
                        &geom2d_curve_tool::parabola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Hyperbola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_ellipse_hyperbola(
                        &geom2d_curve_tool::ellipse(c1),
                        d1,
                        &geom2d_curve_tool::hyperbola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                _ => {
                    self.intconicurv.base.set_reversed_parameters(false);
                    self.intconicurv.perform_ellipse(
                        &geom2d_curve_tool::ellipse(c1),
                        d1,
                        c2,
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
            },
            Curve2dType::Parabola => match typ2 {
                Curve2dType::Line => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_line_parabola(
                        &geom2d_curve_tool::line(c2),
                        d2,
                        &geom2d_curve_tool::parabola(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Circle => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_circle_parabola(
                        &geom2d_curve_tool::circle(c2),
                        d2,
                        &geom2d_curve_tool::parabola(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Ellipse => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_ellipse_parabola(
                        &geom2d_curve_tool::ellipse(c2),
                        d2,
                        &geom2d_curve_tool::parabola(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Parabola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_parabola_parabola(
                        &geom2d_curve_tool::parabola(c1),
                        d1,
                        &geom2d_curve_tool::parabola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Hyperbola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_parabola_hyperbola(
                        &geom2d_curve_tool::parabola(c1),
                        d1,
                        &geom2d_curve_tool::hyperbola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                _ => {
                    self.intconicurv.base.set_reversed_parameters(false);
                    self.intconicurv.perform_parabola(
                        &geom2d_curve_tool::parabola(c1),
                        d1,
                        c2,
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
            },
            Curve2dType::Hyperbola => match typ2 {
                Curve2dType::Line => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_line_hyperbola(
                        &geom2d_curve_tool::line(c2),
                        d2,
                        &geom2d_curve_tool::hyperbola(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Circle => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_circle_hyperbola(
                        &geom2d_curve_tool::circle(c2),
                        d2,
                        &geom2d_curve_tool::hyperbola(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Ellipse => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_ellipse_hyperbola(
                        &geom2d_curve_tool::ellipse(c2),
                        d2,
                        &geom2d_curve_tool::hyperbola(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Parabola => {
                    self.intconiconi.base.set_reversed_parameters(true);
                    self.intconiconi.perform_parabola_hyperbola(
                        &geom2d_curve_tool::parabola(c2),
                        d2,
                        &geom2d_curve_tool::hyperbola(c1),
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                Curve2dType::Hyperbola => {
                    self.intconiconi.base.set_reversed_parameters(false);
                    self.intconiconi.perform_hyperbola_hyperbola(
                        &geom2d_curve_tool::hyperbola(c1),
                        d1,
                        &geom2d_curve_tool::hyperbola(c2),
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_coniconiconi!();
                }
                _ => {
                    self.intconicurv.base.set_reversed_parameters(false);
                    self.intconicurv.perform_hyperbola(
                        &geom2d_curve_tool::hyperbola(c1),
                        d1,
                        c2,
                        d2,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
            },
            _ => match typ2 {
                Curve2dType::Line => {
                    self.intconicurv.base.set_reversed_parameters(true);
                    self.intconicurv.perform_line(
                        &geom2d_curve_tool::line(c2),
                        d2,
                        c1,
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
                Curve2dType::Circle => {
                    self.intconicurv.base.set_reversed_parameters(true);
                    self.intconicurv.perform_circle(
                        &geom2d_curve_tool::circle(c2),
                        d2,
                        c1,
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
                Curve2dType::Ellipse => {
                    self.intconicurv.base.set_reversed_parameters(true);
                    self.intconicurv.perform_ellipse(
                        &geom2d_curve_tool::ellipse(c2),
                        d2,
                        c1,
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
                Curve2dType::Parabola => {
                    self.intconicurv.base.set_reversed_parameters(true);
                    self.intconicurv.perform_parabola(
                        &geom2d_curve_tool::parabola(c2),
                        d2,
                        c1,
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
                Curve2dType::Hyperbola => {
                    self.intconicurv.base.set_reversed_parameters(true);
                    self.intconicurv.perform_hyperbola(
                        &geom2d_curve_tool::hyperbola(c2),
                        d2,
                        c1,
                        d1,
                        tol_conf,
                        tol,
                    );
                    finish_conicurv!();
                }
                _ => {
                    // intcurvcurv.SetReversedParameters(false);
                    // intcurvcurv.Perform(C1, D1, C2, D2, TolConf, Tol);
                    // -- Geom2dInt_TheIntPCurvePCurveOfGInter
                    //    (IntCurve_IntPolyPolyGen.gxx) is not yet ported
                    //    (see file header).
                    let _ = composite;
                    unimplemented!(
                        "IntCurve_IntCurveCurveGen: pcurve x pcurve intersection \
                         (Geom2dInt_TheIntPCurvePCurveOfGInter / IntCurve_IntPolyPolyGen) is \
                         not ported yet"
                    );
                }
            },
        }
    }

    /// OCCT InternalCompositePerform_noRecurs
    /// (IntCurve_IntCurveCurveGen.gxx L818-932).
    #[allow(clippy::too_many_arguments)]
    fn internal_composite_perform_no_recurs(
        &mut self,
        nb_inter_c1: i32,
        c1: &dyn Curve2dAdaptor,
        num_inter_c1: i32,
        tab1: &[f64],
        d1: &Res2dDomain,
        nb_inter_c2: i32,
        c2: &dyn Curve2dAdaptor,
        num_inter_c2: i32,
        tab2: &[f64],
        d2: &Res2dDomain,
        tol_conf: f64,
        tol: f64,
    ) {
        if num_inter_c2 > nb_inter_c2 {
            return;
        }

        let mut domain_c1_num_inter;
        let mut domain_c2_num_inter;

        //----------------------------------------------------------------------
        //-- Creation du domaine associe a la portion de C1
        //----------------------------------------------------------------------
        let mut domain_is_ok = true;
        let (mut param_inf, mut param_sup);

        if nb_inter_c1 > 1 {
            let (pi, ps) = geom2d_curve_tool::get_interval(tab1, num_inter_c1 as usize);
            param_inf = pi;
            param_sup = ps;
            //--------------------------------------------------------------
            //-- Verification : Domaine Inclu dans Intervalle de Definition
            //--------------------------------------------------------------

            let mut u;

            u = d1.first_parameter();
            if param_inf < u {
                param_inf = u;
            }

            u = d1.last_parameter();
            if param_sup > u {
                param_sup = u;
            }

            if (param_sup - param_inf) > 1e-10 {
                domain_c1_num_inter = Res2dDomain::infinite();
                domain_c1_num_inter.set_values_bounded(
                    geom2d_curve_tool::value(c1, param_inf),
                    param_inf,
                    d1.first_tolerance(),
                    geom2d_curve_tool::value(c1, param_sup),
                    param_sup,
                    d1.last_tolerance(),
                );
            } else {
                domain_is_ok = false;
                domain_c1_num_inter = Res2dDomain::infinite();
            }
        } else {
            domain_c1_num_inter = d1.clone();
            param_inf = 0.0;
            param_sup = 0.0;
            let _ = (param_inf, param_sup);
        }

        //----------------------------------------------------------------------
        //-- Creation du domaine associe a la portion de C2
        //----------------------------------------------------------------------
        if nb_inter_c2 > 1 {
            let (pi, ps) = geom2d_curve_tool::get_interval(tab2, num_inter_c2 as usize);
            param_inf = pi;
            param_sup = ps;
            //--------------------------------------------------------------
            //-- Verification : Domaine Inclu dans Intervalle de Definition
            //--------------------------------------------------------------

            let mut u;

            u = d2.first_parameter();
            if param_inf < u {
                param_inf = u;
            }
            u = d2.last_parameter();

            if param_sup > u {
                param_sup = u;
            }

            if (param_sup - param_inf) > 1e-10 {
                domain_c2_num_inter = Res2dDomain::infinite();
                domain_c2_num_inter.set_values_bounded(
                    geom2d_curve_tool::value(c2, param_inf),
                    param_inf,
                    d2.first_tolerance(),
                    geom2d_curve_tool::value(c2, param_sup),
                    param_sup,
                    d2.last_tolerance(),
                );
            } else {
                domain_is_ok = false;
                domain_c2_num_inter = Res2dDomain::infinite();
            }
        } else {
            domain_c2_num_inter = d2.clone();
        }

        if domain_is_ok {
            // OCCT swaps the two curves here (C2 first).
            self.internal_perform(
                c2,
                &domain_c2_num_inter,
                c1,
                &domain_c1_num_inter,
                tol_conf,
                tol,
                true,
            );
        }
    }

    //-- C1 ou C2 sont des courbes composites
    //--
    /// OCCT InternalCompositePerform (IntCurve_IntCurveCurveGen.gxx L937-1019).
    #[allow(clippy::too_many_arguments)]
    fn internal_composite_perform(
        &mut self,
        c1: &dyn Curve2dAdaptor,
        d1: &Res2dDomain,
        xxx_num_inter_c1: i32,
        nb_inter_c1: i32,
        tab1: &[f64],
        c2: &dyn Curve2dAdaptor,
        d2: &Res2dDomain,
        xxx_num_inter_c2: i32,
        nb_inter_c2: i32,
        tab2: &[f64],
        tol_conf: f64,
        tol: f64,
        recurs_on_c2: bool,
    ) {
        let mut num_inter_c2 = xxx_num_inter_c2;
        let mut num_inter_c1 = xxx_num_inter_c1;

        if num_inter_c2 > nb_inter_c2 {
            return;
        }

        if !recurs_on_c2 {
            self.internal_composite_perform_no_recurs(
                nb_inter_c1, c1, num_inter_c1, tab1, d1, nb_inter_c2, c2, num_inter_c2, tab2, d2,
                tol_conf, tol,
            );
            return;
        }

        for i in num_inter_c1..=nb_inter_c1 {
            num_inter_c1 = i;

            self.internal_composite_perform_no_recurs(
                nb_inter_c2, c2, num_inter_c2, tab2, d2, nb_inter_c1, c1, num_inter_c1, tab1, d1,
                tol_conf, tol,
            );
        }

        if num_inter_c2 < nb_inter_c2 {
            num_inter_c2 += 1;
            num_inter_c1 = 1;

            self.internal_composite_perform(
                c1, d1, num_inter_c1, nb_inter_c1, tab1, c2, d2, num_inter_c2, nb_inter_c2, tab2,
                tol_conf, tol, true,
            );
        }
    }
}

impl Default for IntCurveCurveGen {
    fn default() -> Self {
        IntCurveCurveGen::new()
    }
}
