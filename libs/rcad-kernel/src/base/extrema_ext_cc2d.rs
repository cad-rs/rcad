//! OCCT Extrema_ExtCC2d (TKGeomBase/Extrema/Extrema_ExtCC2d.hxx L36-143 and
//! Extrema_ExtCC2d.cxx L34-685) — all the distances between two 2D curves.
//!
//! The elementary branches run [`ExtremaExtElC2d`] (Extrema_ExtElC2d, see
//! [`super::extrema_ext_elc2d`]); the general branches run `myECC`
//! (= `Extrema_ECC2d`, see [`super::extrema_gen_ext_cc2d`]).
//!
//! GAP (module-level): the `Extrema_ECC2d` instantiation carries one
//! untranslated dependency — the `Extrema_ExtPC2d` projection refinement of
//! the parallel-verification block (Extrema_GGenExtCC.hxx L771-783); on that
//! path the solver abandons with the OCCT not-done equivalent and
//! `Extrema_ExtCC2d::Results(Extrema_ECC2d&)` propagates `myDone == false`
//! exactly as the OCCT body does when `AlgExt.IsDone()` is false (cxx
//! L621-622).  All elementary (analytic) branches are complete.

use glam::DVec2;

use crate::base::extrema::POnCurve2d;
use crate::base::extrema_curve2d_tool::ExtremaCurve2dTool;
use crate::base::extrema_ext_elc::elclib_in_period;
use crate::base::extrema_ext_elc2d::ExtremaExtElC2d;
use crate::base::extrema_gen_ext_cc2d::GenExtCC2d;
use crate::base::proj_lib::CurveType;
use crate::core::precision::PCONFUSION;

/// OCCT 2 * M_PI — the period literal the Perform switch passes to Results.
const TWO_PI: f64 = std::f64::consts::PI + std::f64::consts::PI;

/// OCCT Extrema_ExtCC2d (hxx L36-143).
pub struct ExtremaExtCC2d<'a> {
    /// hxx L121: bool myIsFindSingleSolution — default value is false.
    my_is_find_single_solution: bool,
    /// hxx L122: bool myDone.
    my_done: bool,
    /// hxx L123: bool myIsPar.
    my_is_par: bool,
    /// hxx L124: NCollection_Sequence<Extrema_POnCurv2d> mypoints.
    my_points: Vec<POnCurve2d>,
    /// hxx L125: NCollection_Sequence<double> mySqDist.
    my_sq_dist: Vec<f64>,
    /// hxx L126: int mynbext.
    my_nb_ext: usize,
    /// hxx L127: bool inverse.
    inverse: bool,
    /// hxx L128: const Adaptor2d_Curve2d* myC.
    my_c: Option<&'a dyn ExtremaCurve2dTool>,
    /// hxx L129-130: double myv1, myv2.
    my_v1: f64,
    my_v2: f64,
    /// hxx L131-132: double mytolc1, mytolc2.
    my_tolc1: f64,
    my_tolc2: f64,
    /// hxx L133-136: gp_Pnt2d P1f, P1l, P2f, P2l.
    p1f: DVec2,
    p1l: DVec2,
    p2f: DVec2,
    p2l: DVec2,
    /// hxx L137-140: double mydist11, mydist12, mydist21, mydist22.
    my_dist11: f64,
    my_dist12: f64,
    my_dist21: f64,
    my_dist22: f64,
}

impl<'a> ExtremaExtCC2d<'a> {
    /// OCCT Extrema_ExtCC2d() (cxx L34-50).
    pub fn new() -> Self {
        ExtremaExtCC2d {
            my_is_find_single_solution: false,
            my_done: false,
            my_is_par: false,
            my_points: Vec::new(),
            my_sq_dist: Vec::new(),
            my_nb_ext: 0,
            inverse: false,
            my_c: None,
            my_v1: 0.0,
            my_v2: 0.0,
            my_tolc1: 0.0,
            my_tolc2: 0.0,
            p1f: DVec2::ZERO,
            p1l: DVec2::ZERO,
            p2f: DVec2::ZERO,
            p2l: DVec2::ZERO,
            my_dist11: 0.0,
            my_dist12: 0.0,
            my_dist21: 0.0,
            my_dist22: 0.0,
        }
    }

    /// OCCT Extrema_ExtCC2d(C1, C2, TolC1 = 1.0e-10, TolC2 = 1.0e-10)
    /// (cxx L54-66).
    pub fn new_curves(
        c1: &'a dyn ExtremaCurve2dTool,
        c2: &'a dyn ExtremaCurve2dTool,
        tol_c1: f64,
        tol_c2: f64,
    ) -> Self {
        let mut this = ExtremaExtCC2d::new();
        this.initialize(c2, c2.first_parameter(), c2.last_parameter(), tol_c1, tol_c2);
        this.perform(c1, c1.first_parameter(), c1.last_parameter());
        this
    }

    /// OCCT Extrema_ExtCC2d(C1, C2, U1, U2, V1, V2, TolC1 = 1.0e-10,
    /// TolC2 = 1.0e-10) (cxx L70-82).
    #[allow(clippy::too_many_arguments)]
    pub fn new_curves_ranged(
        c1: &'a dyn ExtremaCurve2dTool,
        c2: &'a dyn ExtremaCurve2dTool,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        tol_c1: f64,
        tol_c2: f64,
    ) -> Self {
        let mut this = ExtremaExtCC2d::new();
        this.initialize(c2, v1, v2, tol_c1, tol_c2);
        this.perform(c1, u1, u2);
        this
    }

    /// OCCT Initialize(C2, V1, V2, TolC1 = 1.0e-10, TolC2 = 1.0e-10)
    /// (cxx L86-97).
    pub fn initialize(
        &mut self,
        c2: &'a dyn ExtremaCurve2dTool,
        v1: f64,
        v2: f64,
        tol_c1: f64,
        tol_c2: f64,
    ) {
        self.my_c = Some(c2);
        self.my_v1 = v1;
        self.my_v2 = v2;
        self.my_tolc1 = tol_c1;
        self.my_tolc2 = tol_c2;
    }

    /// OCCT Perform(C1, U1, U2) (cxx L101-453).
    pub fn perform(&mut self, c1: &'a dyn ExtremaCurve2dTool, u1: f64, u2: f64) {
        // cxx L103-104.
        self.my_points.clear();
        self.my_sq_dist.clear();

        // cxx L105-109.
        let type1 = c1.get_type();
        let type2 = self.my_c.expect("Extrema_ExtCC2d: myC").get_type();
        let a_tol = self.my_tolc1.min(self.my_tolc2);
        self.my_nb_ext = 0;
        self.inverse = false;
        self.my_is_par = false;

        // cxx L113-120.
        let a_u11 = u1;
        let a_u12 = u2;
        let a_u21 = self.my_v1;
        let a_u22 = self.my_v2;
        self.p1f = c1.value(a_u11);
        self.p1l = c1.value(a_u12);
        let a_my_c = self.my_c.expect("Extrema_ExtCC2d: myC");
        self.p2f = a_my_c.value(a_u21);
        self.p2l = a_my_c.value(a_u22);

        // cxx L122-452: the type dispatch; `a_xtream` / `a_param_solver` are
        // the OCCT local shared_ptr slots.
        match type1 {
            // cxx L130-181: the first curve is a circle.
            CurveType::Circle => match type2 {
                CurveType::Line => {
                    self.inverse = true;
                    let a_c2 = a_my_c.line();
                    let a_c1 = c1.circle();
                    let an_xtream = ExtremaExtElC2d::line_circle(&a_c2, &a_c1, a_tol);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, TWO_PI, 0.0);
                }
                CurveType::Circle => {
                    let a_c1 = c1.circle();
                    let a_c2 = a_my_c.circle();
                    let an_xtream = ExtremaExtElC2d::circle_circle(&a_c1, &a_c2);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, TWO_PI, TWO_PI);
                }
                CurveType::Ellipse => {
                    let a_c1 = c1.circle();
                    let a_c2 = a_my_c.ellipse();
                    let an_xtream = ExtremaExtElC2d::circle_ellipse(&a_c1, &a_c2);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, TWO_PI, TWO_PI);
                }
                CurveType::Parabola => {
                    let a_c1 = c1.circle();
                    let a_c2 = a_my_c.parabola();
                    let an_xtream = ExtremaExtElC2d::circle_parabola(&a_c1, &a_c2);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, TWO_PI, 0.0);
                }
                CurveType::Hyperbola => {
                    let a_c1 = c1.circle();
                    let a_c2 = a_my_c.hyperbola();
                    let an_xtream = ExtremaExtElC2d::circle_hyperbola(&a_c1, &a_c2);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, TWO_PI, 0.0);
                }
                _ => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    let mut period2 = 0.0;
                    if a_my_c.is_periodic() {
                        period2 = a_my_c.period();
                    }
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, TWO_PI, period2);
                }
            },
            // cxx L186-243: the first curve is an ellipse.
            CurveType::Ellipse => match type2 {
                CurveType::Line => {
                    self.inverse = true;
                    let a_c2 = a_my_c.line();
                    let a_c1 = c1.ellipse();
                    let an_xtream = ExtremaExtElC2d::line_ellipse(&a_c2, &a_c1);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, TWO_PI, 0.0);
                }
                CurveType::Circle => {
                    self.inverse = true;
                    let a_c2 = a_my_c.circle();
                    let a_c1 = c1.ellipse();
                    let an_xtream = ExtremaExtElC2d::circle_ellipse(&a_c2, &a_c1);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, TWO_PI, TWO_PI);
                }
                CurveType::Ellipse => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, TWO_PI, TWO_PI);
                }
                CurveType::Parabola | CurveType::Hyperbola => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, TWO_PI, 0.0);
                }
                _ => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    let mut period2 = 0.0;
                    if a_my_c.is_periodic() {
                        period2 = a_my_c.period();
                    }
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, TWO_PI, period2);
                }
            },
            // cxx L248-309: the first curve is a parabola.
            CurveType::Parabola => match type2 {
                CurveType::Line => {
                    self.inverse = true;
                    let a_c2 = a_my_c.line();
                    let a_c1 = c1.parabola();
                    let an_xtream = ExtremaExtElC2d::line_parabola(&a_c2, &a_c1);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, 0.0);
                }
                CurveType::Circle => {
                    self.inverse = true;
                    let a_c2 = a_my_c.circle();
                    let a_c1 = c1.parabola();
                    let an_xtream = ExtremaExtElC2d::circle_parabola(&a_c2, &a_c1);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, TWO_PI);
                }
                CurveType::Ellipse => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, 0.0, TWO_PI);
                }
                CurveType::Parabola | CurveType::Hyperbola => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, 0.0, 0.0);
                }
                _ => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    let mut period2 = 0.0;
                    if a_my_c.is_periodic() {
                        period2 = a_my_c.period();
                    }
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, 0.0, period2);
                }
            },
            // cxx L314-374: the first curve is a hyperbola.
            CurveType::Hyperbola => match type2 {
                CurveType::Line => {
                    self.inverse = true;
                    let a_c2 = a_my_c.line();
                    let a_c1 = c1.hyperbola();
                    let an_xtream = ExtremaExtElC2d::line_hyperbola(&a_c2, &a_c1);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, 0.0);
                }
                CurveType::Circle => {
                    self.inverse = true;
                    let a_c2 = a_my_c.circle();
                    let a_c1 = c1.hyperbola();
                    let an_xtream = ExtremaExtElC2d::circle_hyperbola(&a_c2, &a_c1);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, TWO_PI);
                }
                CurveType::Ellipse => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, 0.0, TWO_PI);
                }
                CurveType::Parabola | CurveType::Hyperbola => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, 0.0, 0.0);
                }
                _ => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    let mut period2 = 0.0;
                    if a_my_c.is_periodic() {
                        period2 = a_my_c.period();
                    }
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, 0.0, period2);
                }
            },
            // cxx L379-430: the first curve is a line.
            CurveType::Line => match type2 {
                CurveType::Line => {
                    let a_c1 = c1.line();
                    let a_c2 = a_my_c.line();
                    let an_xtream = ExtremaExtElC2d::line_line(&a_c1, &a_c2, a_tol);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, 0.0);
                }
                CurveType::Circle => {
                    let a_c1 = c1.line();
                    let a_c2 = a_my_c.circle();
                    let an_xtream = ExtremaExtElC2d::line_circle(&a_c1, &a_c2, a_tol);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, TWO_PI);
                }
                CurveType::Ellipse => {
                    let a_c1 = c1.line();
                    let a_c2 = a_my_c.ellipse();
                    let an_xtream = ExtremaExtElC2d::line_ellipse(&a_c1, &a_c2);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, TWO_PI);
                }
                CurveType::Parabola => {
                    let a_c1 = c1.line();
                    let a_c2 = a_my_c.parabola();
                    let an_xtream = ExtremaExtElC2d::line_parabola(&a_c1, &a_c2);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, 0.0);
                }
                CurveType::Hyperbola => {
                    let a_c1 = c1.line();
                    let a_c2 = a_my_c.hyperbola();
                    let an_xtream = ExtremaExtElC2d::line_hyperbola(&a_c1, &a_c2);
                    self.results_elc(&an_xtream, a_u11, a_u12, a_u21, a_u22, 0.0, 0.0);
                }
                _ => {
                    let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                    a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                    a_param_solver.perform();
                    let mut period2 = 0.0;
                    if a_my_c.is_periodic() {
                        period2 = a_my_c.period();
                    }
                    self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, 0.0, period2);
                }
            },
            // cxx L435-452: the first curve is a BezierCurve or a BSplineCurve.
            _ => {
                let mut a_param_solver = GenExtCC2d::new(c1, a_my_c);
                a_param_solver.set_single_solution_flag(self.get_single_solution_flag());
                a_param_solver.perform();
                let mut period1 = 0.0;
                if c1.is_periodic() {
                    period1 = c1.period();
                }
                let mut period2 = 0.0;
                if a_my_c.is_periodic() {
                    period2 = a_my_c.period();
                }
                self.results_ecc(&a_param_solver, a_u11, a_u12, a_u21, a_u22, period1, period2);
            }
        }
    }

    /// OCCT IsDone() (cxx L457-460).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT SquareDistance(N) (cxx L464-475) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if !self.my_done {
            panic!("StdFail_NotDone");
        }
        if n == 0 || n > self.my_nb_ext {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT NbExt() (cxx L479-486).
    pub fn nb_ext(&self) -> usize {
        if !self.my_done {
            panic!("StdFail_NotDone");
        }
        self.my_nb_ext
    }

    /// OCCT Points(N, P1, P2) (cxx L490-502) — 1-based.
    pub fn points(&self, n: usize, p1: &mut POnCurve2d, p2: &mut POnCurve2d) {
        if !self.my_done {
            panic!("StdFail_NotDone");
        }
        if n == 0 || n > self.my_nb_ext {
            panic!("Standard_OutOfRange");
        }
        *p1 = self.my_points[2 * n - 2].clone();
        *p2 = self.my_points[2 * n - 1].clone();
    }

    /// OCCT TrimmedSquareDistances(dist11, dist12, dist21, dist22, P11, P12,
    /// P21, P22) (cxx L506-523).
    #[allow(clippy::type_complexity)]
    pub fn trimmed_square_distances(
        &self,
    ) -> (f64, f64, f64, f64, DVec2, DVec2, DVec2, DVec2) {
        (
            self.my_dist11,
            self.my_dist12,
            self.my_dist21,
            self.my_dist22,
            self.p1f,
            self.p1l,
            self.p2f,
            self.p2l,
        )
    }

    /// OCCT Results(const Extrema_ExtElC2d& AlgExt, Ut11, Ut12, Ut21, Ut22,
    /// Period1 = 0.0, Period2 = 0.0) (cxx L527-605).
    pub(crate) fn results_elc(
        &mut self,
        alg_ext: &ExtremaExtElC2d,
        ut11: f64,
        ut12: f64,
        ut21: f64,
        ut22: f64,
        period1: f64,
        period2: f64,
    ) {
        // cxx L539-540.
        self.my_done = alg_ext.is_done();
        self.my_is_par = alg_ext.is_parallel();
        if self.my_done {
            if !self.my_is_par {
                // cxx L545-597: verification of the validity of parameters
                // for the trimmed case.
                let nb_ext = alg_ext.nb_ext();
                for an_i in 1..=nb_ext {
                    let mut a_p1 = POnCurve2d {
                        param: 0.0,
                        point: DVec2::ZERO,
                    };
                    let mut a_p2 = POnCurve2d {
                        param: 0.0,
                        point: DVec2::ZERO,
                    };
                    alg_ext.points(an_i, &mut a_p1, &mut a_p2);
                    let mut a_u;
                    let mut a_u2;
                    if !self.inverse {
                        a_u = a_p1.param;
                        if period1 != 0.0 {
                            a_u = elclib_in_period(a_u, ut11, ut11 + period1);
                        }
                        a_u2 = a_p2.param;
                        if period2 != 0.0 {
                            a_u2 = elclib_in_period(a_u2, ut21, ut21 + period2);
                        }
                    } else {
                        a_u2 = a_p1.param;
                        if period2 != 0.0 {
                            a_u2 = elclib_in_period(a_u2, ut21, ut21 + period2);
                        }
                        a_u = a_p2.param;
                        if period1 != 0.0 {
                            a_u = elclib_in_period(a_u, ut11, ut11 + period1);
                        }
                    }
                    if (a_u >= ut11 - PCONFUSION)
                        && (a_u <= ut12 + PCONFUSION)
                        && (a_u2 >= ut21 - PCONFUSION)
                        && (a_u2 <= ut22 + PCONFUSION)
                    {
                        self.my_nb_ext += 1;
                        let a_val = alg_ext.square_distance(an_i);
                        self.my_sq_dist.push(a_val);
                        if !self.inverse {
                            a_p1.param = a_u;
                            a_p2.param = a_u2;
                            self.my_points.push(a_p1);
                            self.my_points.push(a_p2);
                        } else {
                            a_p1.param = a_u2;
                            a_p2.param = a_u;
                            self.my_points.push(a_p2);
                            self.my_points.push(a_p1);
                        }
                    }
                }
            }

            // cxx L600-603.
            self.my_dist11 = self.p1f.distance_squared(self.p2f);
            self.my_dist12 = self.p1f.distance_squared(self.p2l);
            self.my_dist21 = self.p1l.distance_squared(self.p2f);
            self.my_dist22 = self.p1l.distance_squared(self.p2l);
        }
    }

    /// OCCT Results(const Extrema_ECC2d& AlgExt, Ut11, Ut12, Ut21, Ut22,
    /// Period1 = 0.0, Period2 = 0.0) (cxx L609-659).
    pub(crate) fn results_ecc(
        &mut self,
        alg_ext: &GenExtCC2d,
        ut11: f64,
        ut12: f64,
        ut21: f64,
        ut22: f64,
        period1: f64,
        period2: f64,
    ) {
        // cxx L621.
        self.my_done = alg_ext.is_done();
        if self.my_done {
            // cxx L624-626.
            self.my_is_par = alg_ext.is_parallel();
            let nb_ext = alg_ext.nb_ext();
            for an_i in 1..=nb_ext {
                // cxx L629-639: verification of parameter validity for the
                // trimmed case.
                let mut a_p1 = POnCurve2d {
                    param: 0.0,
                    point: DVec2::ZERO,
                };
                let mut a_p2 = POnCurve2d {
                    param: 0.0,
                    point: DVec2::ZERO,
                };
                alg_ext.points(an_i, &mut a_p1, &mut a_p2);
                let mut a_u = a_p1.param;
                if period1 != 0.0 {
                    a_u = elclib_in_period(a_u, ut11, ut11 + period1);
                }
                let mut a_u2 = a_p2.param;
                if period2 != 0.0 {
                    a_u2 = elclib_in_period(a_u2, ut21, ut21 + period2);
                }

                if (a_u >= ut11 - PCONFUSION)
                    && (a_u <= ut12 + PCONFUSION)
                    && (a_u2 >= ut21 - PCONFUSION)
                    && (a_u2 <= ut22 + PCONFUSION)
                {
                    self.my_nb_ext += 1;
                    let a_val = alg_ext.square_distance(an_i);
                    a_p1.param = a_u;
                    a_p2.param = a_u2;
                    self.my_sq_dist.push(a_val);
                    self.my_points.push(a_p1);
                    self.my_points.push(a_p2);
                }
            }

            // cxx L654-657.
            self.my_dist11 = self.p1f.distance_squared(self.p2f);
            self.my_dist12 = self.p1f.distance_squared(self.p2l);
            self.my_dist21 = self.p1l.distance_squared(self.p2f);
            self.my_dist22 = self.p1l.distance_squared(self.p2l);
        }
    }

    /// OCCT IsParallel() (cxx L663-670).
    pub fn is_parallel(&self) -> bool {
        if !self.my_done {
            panic!("StdFail_NotDone");
        }
        self.my_is_par
    }

    /// OCCT SetSingleSolutionFlag(theFlag) (cxx L674-677).
    pub fn set_single_solution_flag(&mut self, the_flag: bool) {
        self.my_is_find_single_solution = the_flag;
    }

    /// OCCT GetSingleSolutionFlag() (cxx L681-684).
    pub fn get_single_solution_flag(&self) -> bool {
        self.my_is_find_single_solution
    }
}

impl Default for ExtremaExtCC2d<'_> {
    fn default() -> Self {
        ExtremaExtCC2d::new()
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::proj_lib::adaptor::Geom2dCurveAdaptor;
    use crate::geom::{Circle2d, Curve2d, Line2d};

    /// The minimum over the extrema the class reports — the OCCT accessor
    /// sequence `IsDone() / NbExt() / SquareDistance(i)`.
    fn min_sq_dist(the_ext: &ExtremaExtCC2d) -> f64 {
        let mut a_min = f64::INFINITY;
        for an_idx in 1..=the_ext.nb_ext() {
            a_min = a_min.min(the_ext.square_distance(an_idx));
        }
        a_min
    }

    /// Extrema_ExtCC2d(gp_Lin2d, gp_Lin2d): two lines crossing at the origin
    /// -> distance 0 at parameters 0 / 0.  The accessors are called in the
    /// OCCT order Perform -> IsDone -> NbExt -> SquareDistance -> Points
    /// (the regression guard for the consumer API shape).
    #[test]
    fn two_lines_crossing_at_origin() {
        let a_line1 = Curve2d::Line(Line2d::new(DVec2::ZERO, DVec2::X));
        let a_line2 = Curve2d::Line(Line2d::new(DVec2::ZERO, DVec2::Y));
        let an_ad1 = Geom2dCurveAdaptor::with_range(a_line1, -10.0, 10.0);
        let an_ad2 = Geom2dCurveAdaptor::with_range(a_line2, -10.0, 10.0);

        // OCCT Extrema_ExtCC2d(C1, C2, TolC1, TolC2) — constructor runs
        // Initialize + Perform.
        let an_ext = ExtremaExtCC2d::new_curves(&an_ad1, &an_ad2, 1.0e-10, 1.0e-10);

        assert!(an_ext.is_done());
        assert!(!an_ext.is_parallel());
        assert_eq!(an_ext.nb_ext(), 1);
        let a_d = an_ext.square_distance(1);
        assert!(a_d < 1.0e-20, "expected 0 distance, got {a_d}");

        let mut a_p1 = POnCurve2d {
            param: 0.0,
            point: DVec2::ZERO,
        };
        let mut a_p2 = POnCurve2d {
            param: 0.0,
            point: DVec2::ZERO,
        };
        an_ext.points(1, &mut a_p1, &mut a_p2);
        assert!(a_p1.point.distance(DVec2::ZERO) < 1.0e-12);
        assert!(a_p2.point.distance(DVec2::ZERO) < 1.0e-12);
        assert!(a_p1.param.abs() < 1.0e-12);
        assert!(a_p2.param.abs() < 1.0e-12);
    }

    /// Extrema_ExtCC2d(gp_Lin2d, gp_Circ2d): the line y = 3 against the
    /// circle R = 1 at the origin -> the nearest point of the circle is
    /// (0, 1) at parameter PI/2, the minimum distance is 2, and the far
    /// branch at 3*PI/2 reports 4.
    #[test]
    fn line_circle_known_minimum() {
        let a_line = Curve2d::Line(Line2d::new(DVec2::new(0.0, 3.0), DVec2::X));
        let a_circle = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 1.0));
        let an_ad1 = Geom2dCurveAdaptor::with_range(a_line, -10.0, 10.0);
        let an_ad2 = Geom2dCurveAdaptor::with_range(a_circle, 0.0, TWO_PI);

        let an_ext = ExtremaExtCC2d::new_curves(&an_ad1, &an_ad2, 1.0e-10, 1.0e-10);

        assert!(an_ext.is_done());
        assert!(!an_ext.is_parallel());
        assert_eq!(an_ext.nb_ext(), 2);
        // The nearest branch: distance 3 - 1 = 2 at the circle parameter
        // PI/2 (0.5 * PI within the [0, 2*PI] window).
        let a_d1 = an_ext.square_distance(1).sqrt();
        assert!((a_d1 - 2.0).abs() < 1.0e-9, "expected 2.0, got {a_d1}");
        // The far branch: distance 3 + 1 = 4 at PI/2 + PI.
        let a_d2 = an_ext.square_distance(2).sqrt();
        assert!((a_d2 - 4.0).abs() < 1.0e-9, "expected 4.0, got {a_d2}");

        let mut a_p1 = POnCurve2d {
            param: 0.0,
            point: DVec2::ZERO,
        };
        let mut a_p2 = POnCurve2d {
            param: 0.0,
            point: DVec2::ZERO,
        };
        an_ext.points(1, &mut a_p1, &mut a_p2);
        assert!(a_p2.point.distance(DVec2::new(0.0, 1.0)) < 1.0e-9);
        assert!(a_p1.point.distance(DVec2::new(0.0, 3.0)) < 1.0e-9);
        assert!(a_p1.param.abs() < 1.0e-9);
        assert!((a_p2.param - std::f64::consts::FRAC_PI_2).abs() < 1.0e-9);
    }

    /// Extrema_ExtCC2d(gp_Lin2d, gp_Lin2d): parallel lines report
    /// IsParallel() with no extrema stored and the trimmed corner distances
    /// equal to the offset squared (the OCCT Results arm skips the extrema
    /// loop when myIsPar, cxx L543/L598).
    #[test]
    fn parallel_lines_report_parallel_and_trimmed_distances() {
        let a_line1 = Curve2d::Line(Line2d::new(DVec2::ZERO, DVec2::X));
        let a_line2 = Curve2d::Line(Line2d::new(DVec2::new(0.0, 2.0), DVec2::X));
        let an_ad1 = Geom2dCurveAdaptor::with_range(a_line1, -10.0, 10.0);
        let an_ad2 = Geom2dCurveAdaptor::with_range(a_line2, -10.0, 10.0);

        let an_ext = ExtremaExtCC2d::new_curves(&an_ad1, &an_ad2, 1.0e-10, 1.0e-10);

        assert!(an_ext.is_done());
        assert!(an_ext.is_parallel());
        assert_eq!(an_ext.nb_ext(), 0);
        let (d11, d12, d21, d22, _p11, _p12, _p21, _p22) = an_ext.trimmed_square_distances();
        assert!((d11 - 4.0).abs() < 1.0e-12);
        assert!((d12 - 404.0).abs() < 1.0e-9);
        assert!((d21 - 404.0).abs() < 1.0e-9);
        assert!((d22 - 4.0).abs() < 1.0e-12);
        let _ = min_sq_dist(&an_ext);
    }

    /// Extrema_ExtCC2d(gp_Lin2d, gp_Circ2d) over a line window that misses
    /// the circle: both elementary extrema fall outside the trimmed ranges
    /// and the Results filter drops them (cxx L576-597) — nbext stays 0.
    #[test]
    fn trimmed_window_filter_drops_out_of_range_extrema() {
        // The line spans x in [20, 30]; the extrema sit at x = 0.
        let a_line = Curve2d::Line(Line2d::new(DVec2::new(20.0, 3.0), DVec2::X));
        let a_circle = Curve2d::Circle(Circle2d::new(DVec2::ZERO, 1.0));
        let an_ad1 = Geom2dCurveAdaptor::with_range(a_line, 0.0, 10.0);
        let an_ad2 = Geom2dCurveAdaptor::with_range(a_circle, 0.0, TWO_PI);

        let an_ext = ExtremaExtCC2d::new_curves_ranged(
            &an_ad1, &an_ad2, 0.0, 10.0, 0.0, TWO_PI, 1.0e-10, 1.0e-10,
        );

        assert!(an_ext.is_done());
        assert_eq!(an_ext.nb_ext(), 0);
    }
}

/// The `BRepExtrema_ExtCF` regression guard (BRepFeat_MakeLinearForm sliding
/// check, BRepFeat_MakeLinearForm.cxx L307-345).
///
/// Survey finding (architecture difference #3 of the feat carrier): OCCT has
/// no `Extrema_ExtCF` class — `BRepExtrema_ExtCF` (TKTopAlgo) is the
/// topological wrapper whose kernel engine is `Extrema_ExtCS`
/// (TKGeomBase/Extrema/Extrema_ExtCS.cxx, `Initialize` L90-107 + `Perform`
/// L109-...), already translated in [`crate::base::extrema_ext_cs`] with real
/// analytic bodies.  These tests exercise that kernel class in the exact
/// OCCT call order the wrapper uses (`Initialize` -> `Perform` -> `IsDone` ->
/// `NbExt` -> `SquareDistance`), on the sliding scenarios the feat consumer
/// needs: a line inside / above its support plane.
#[cfg(test)]
mod extcf_regression {
    use crate::base::extrema_curve_tool::CurveToolHandle;
    use crate::base::extrema_ext_cs::ExtremaExtCS;
    use crate::base::proj_lib::geom_adaptor_curve::GeomCurveAdaptor;
    use crate::base::proj_lib::geom_adaptor_surface::GeomSurfaceAdaptor;
    use glam::DVec3;
    use crate::geom::{Curve3, Line3, Plane, Surface3};

    /// BRepFeat_MakeLinearForm.cxx L307-311: the sliding probe of the segment
    /// (myFirstPnt -> myFirstPnt + myDir) against the (planar) FirstFace.
    /// The segment lies IN the plane -> the extrema are parallel with the
    /// single square distance 0 <= tol^2.
    #[test]
    fn sliding_line_inside_plane() {
        let a_curve3 = Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, 0.0), DVec3::X));
        let an_adaptor = GeomCurveAdaptor::with_range(a_curve3.clone(), 0.0, 1.0);
        let a_tool = CurveToolHandle::for_curve3(&a_curve3, &an_adaptor, &an_adaptor);

        let a_surface = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let a_surf_adaptor = GeomSurfaceAdaptor::new(a_surface);

        // OCCT BRepExtrema_ExtCF::Initialize: myExtCS.Initialize(S, U1, U2,
        // V1, V2, aTolC, aTolS) over the face UV bounds.
        let mut an_ext_cs = ExtremaExtCS::new();
        an_ext_cs.initialize_range(&a_surf_adaptor, -10.0, 10.0, -10.0, 10.0, 1.0e-7, 1.0e-7);
        // OCCT BRepExtrema_ExtCF::Perform: myExtCS.Perform(C, U1, U2).
        an_ext_cs.perform(&a_tool, 0.0, 1.0);

        // The consumer read-back, exactly BRepFeat_MakeLinearForm.cxx
        // L309-312: NbExt() == 1 && SquareDistance(1) <= tol^2.
        assert!(an_ext_cs.is_done());
        assert!(an_ext_cs.is_parallel());
        assert_eq!(an_ext_cs.nb_ext(), 1);
        let a_sq = an_ext_cs.square_distance(1);
        assert!(a_sq <= 1.0e-7 * 1.0e-7, "expected in-plane, got {a_sq}");
    }

    /// The same probe with the segment lifted by 2 above the plane: still one
    /// parallel extremum, the square distance is exactly 4.
    #[test]
    fn sliding_line_above_plane() {
        let a_curve3 = Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, 2.0), DVec3::X));
        let an_adaptor = GeomCurveAdaptor::with_range(a_curve3.clone(), 0.0, 1.0);
        let a_tool = CurveToolHandle::for_curve3(&a_curve3, &an_adaptor, &an_adaptor);

        let a_surface = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let a_surf_adaptor = GeomSurfaceAdaptor::new(a_surface);

        let mut an_ext_cs = ExtremaExtCS::new();
        an_ext_cs.initialize_range(&a_surf_adaptor, -10.0, 10.0, -10.0, 10.0, 1.0e-7, 1.0e-7);
        an_ext_cs.perform(&a_tool, 0.0, 1.0);

        assert!(an_ext_cs.is_done());
        assert!(an_ext_cs.is_parallel());
        assert_eq!(an_ext_cs.nb_ext(), 1);
        let a_sq = an_ext_cs.square_distance(1);
        assert!((a_sq - 4.0).abs() < 1.0e-12, "expected 4.0, got {a_sq}");
    }

    /// A line crossing the plane (non-parallel): OCCT routes it to the
    /// generic search `Extrema_GenExtCS` (Extrema_ExtCS.cxx L141-144), which
    /// is the documented GAP carrier of the [`crate::base::extrema_ext_cs`]
    /// port — the tool answers `IsDone() == false` and the wrapper's
    /// `Perform` keeps the empty result (the preserved OCCT failure path).
    /// This test pins that contract: the sliding consumer falls back to
    /// `Sliding = false` exactly as it does when the generic search finds
    /// nothing.
    #[test]
    fn crossing_line_takes_generic_gap_path() {
        let a_curve3 = Curve3::Line(Line3::new(DVec3::new(0.0, 0.0, 1.0), DVec3::Z));
        let an_adaptor = GeomCurveAdaptor::with_range(a_curve3.clone(), -5.0, 5.0);
        let a_tool = CurveToolHandle::for_curve3(&a_curve3, &an_adaptor, &an_adaptor);

        let a_surface = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));
        let a_surf_adaptor = GeomSurfaceAdaptor::new(a_surface);

        let mut an_ext_cs = ExtremaExtCS::new();
        an_ext_cs.initialize_range(&a_surf_adaptor, -10.0, 10.0, -10.0, 10.0, 1.0e-7, 1.0e-7);
        an_ext_cs.perform(&a_tool, -5.0, 5.0);

        // The GAP carrier: the generic search dependency is not translated,
        // the OCCT `myDone = Ext.IsDone()` propagation (cxx L435-436) keeps
        // the not-done state.  The wrapper answers exactly as
        // BRepExtrema_ExtCF.cxx L80-83: `if (!IsDone()) return;` — NbExt()
        // is never reached (OCCT raises StdFail_NotDone there).
        assert!(!an_ext_cs.is_done());
    }
}
