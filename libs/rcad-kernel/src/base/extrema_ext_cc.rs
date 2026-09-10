//! OCCT Extrema_ExtCC (TKGeomBase/Extrema/Extrema_ExtCC.hxx L34-176 and
//! Extrema_ExtCC.cxx L46-969) — all the distances between two curves.
//!
//! The general (non-elementary) branch runs `myECC`
//! (= `Extrema_GGenExtCC`, see [`super::extrema_gen_ext_cc`]); the elementary
//! branches run `Extrema_ExtElC` and, in `PrepareParallelResult`, the
//! point-elementary-curve extrema `Extrema_ExtPElC`.

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::base::extrema_curve_tool::ExtremaCurveTool;
use crate::base::extrema_ext_elc::{
    elclib_circle_parameter, elclib_circle_value, elclib_in_period, elclib_line_parameter,
    elclib_line_value, epsilon_of, ExtremaExtElC, REAL_FIRST, REAL_LAST,
};
use crate::base::extrema_ext_p_elc::ExtremaExtPElC;
use crate::base::extrema_gen_ext_cc::GenExtCC;
use crate::base::proj_lib::CurveType;
use crate::core::precision::{is_infinite_value, ANGULAR, CONFUSION, SQUARE_CONFUSION};
use crate::geom::{Circle3, Line3};
use crate::math::bnd::range::Range;
use crate::math::GeomAbsShape;

/// OCCT Extrema_ExtCC (hxx L34-176).
pub struct ExtremaExtCC<'a> {
    /// hxx L158: bool myIsFindSingleSolution.
    my_is_find_single_solution: bool,
    /// hxx L161: bool myDone.
    my_done: bool,
    /// hxx L162: bool myIsParallel.
    my_is_parallel: bool,
    /// hxx L163: NCollection_Sequence<Extrema_POnCurv> mypoints.
    my_points: Vec<POnCurve>,
    /// hxx L164: NCollection_Sequence<double> mySqDist.
    my_sq_dist: Vec<f64>,
    /// hxx L165: const Adaptor3d_Curve* myC[2].
    my_c1: Option<&'a dyn ExtremaCurveTool>,
    my_c2: Option<&'a dyn ExtremaCurveTool>,
    /// hxx L166: double myInf[2].
    my_inf: [f64; 2],
    /// hxx L167: double mySup[2].
    my_sup: [f64; 2],
    /// hxx L168: double myTol[2].
    my_tol: [f64; 2],
    /// hxx L169: gp_Pnt myP1f.
    my_p1f: DVec3,
    /// hxx L170: gp_Pnt myP1l.
    my_p1l: DVec3,
    /// hxx L171: gp_Pnt myP2f.
    my_p2f: DVec3,
    /// hxx L172: gp_Pnt myP2l.
    my_p2l: DVec3,
    /// hxx L173: double mydist11.
    my_dist11: f64,
    /// hxx L174: double mydist12.
    my_dist12: f64,
    /// hxx L175: double mydist21.
    my_dist21: f64,
    /// hxx L176: double mydist22.
    my_dist22: f64,
}

impl<'a> ExtremaExtCC<'a> {
    /// OCCT Extrema_ExtCC(TolC1 = 1.0e-10, TolC2 = 1.0e-10) (cxx L46-58).
    pub fn new(tol_c1: f64, tol_c2: f64) -> Self {
        ExtremaExtCC {
            my_is_find_single_solution: false,
            my_done: false,
            my_is_parallel: false,
            my_points: Vec::new(),
            my_sq_dist: Vec::new(),
            my_c1: None,
            my_c2: None,
            my_inf: [-crate::core::precision::INFINITE_VALUE; 2],
            my_sup: [crate::core::precision::INFINITE_VALUE; 2],
            my_tol: [tol_c1, tol_c2],
            my_p1f: DVec3::ZERO,
            my_p1l: DVec3::ZERO,
            my_p2f: DVec3::ZERO,
            my_p2l: DVec3::ZERO,
            my_dist11: REAL_FIRST,
            my_dist12: REAL_FIRST,
            my_dist21: REAL_FIRST,
            my_dist22: REAL_FIRST,
        }
    }

    /// OCCT Extrema_ExtCC(C1, C2, TolC1, TolC2) (cxx L84-98).
    pub fn new_curves(
        c1: &'a dyn ExtremaCurveTool,
        c2: &'a dyn ExtremaCurveTool,
        tol_c1: f64,
        tol_c2: f64,
    ) -> Self {
        let mut this = ExtremaExtCC::new(tol_c1, tol_c2);
        this.set_curve3(
            1,
            c1,
            c1.first_parameter(),
            c1.last_parameter(),
        );
        this.set_curve3(2, c2, c2.first_parameter(), c2.last_parameter());
        this.set_tolerance(1, tol_c1);
        this.set_tolerance(2, tol_c2);
        this.my_dist11 = REAL_FIRST;
        this.my_dist12 = REAL_FIRST;
        this.my_dist21 = REAL_FIRST;
        this.my_dist22 = REAL_FIRST;
        this.perform();
        this
    }

    /// OCCT Extrema_ExtCC(C1, C2, U1, U2, V1, V2, TolC1, TolC2)
    /// (cxx L62-80).
    pub fn new_curves_ranged(
        c1: &'a dyn ExtremaCurveTool,
        c2: &'a dyn ExtremaCurveTool,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        tol_c1: f64,
        tol_c2: f64,
    ) -> Self {
        let mut this = ExtremaExtCC::new(tol_c1, tol_c2);
        this.set_curve3(1, c1, u1, u2);
        this.set_curve3(2, c2, v1, v2);
        this.set_tolerance(1, tol_c1);
        this.set_tolerance(2, tol_c2);
        this.my_dist11 = REAL_FIRST;
        this.my_dist12 = REAL_FIRST;
        this.my_dist21 = REAL_FIRST;
        this.my_dist22 = REAL_FIRST;
        this.perform();
        this
    }

    /// OCCT Initialize(C1, C2, TolC1, TolC2) (cxx L102-114).
    pub fn initialize(
        &mut self,
        c1: &'a dyn ExtremaCurveTool,
        c2: &'a dyn ExtremaCurveTool,
        tol_c1: f64,
        tol_c2: f64,
    ) {
        self.my_done = false;
        self.set_curve3(1, c1, c1.first_parameter(), c1.last_parameter());
        self.set_curve3(2, c2, c2.first_parameter(), c2.last_parameter());
        self.set_tolerance(1, tol_c1);
        self.set_tolerance(2, tol_c2);
        self.my_dist11 = REAL_FIRST;
        self.my_dist12 = REAL_FIRST;
        self.my_dist21 = REAL_FIRST;
        self.my_dist22 = REAL_FIRST;
    }

    /// OCCT Initialize(C1, C2, U1, U2, V1, V2, TolC1, TolC2) (cxx L118-134).
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_ranged(
        &mut self,
        c1: &'a dyn ExtremaCurveTool,
        c2: &'a dyn ExtremaCurveTool,
        u1: f64,
        u2: f64,
        v1: f64,
        v2: f64,
        tol_c1: f64,
        tol_c2: f64,
    ) {
        self.my_done = false;
        self.set_curve3(1, c1, u1, u2);
        self.set_curve3(2, c2, v1, v2);
        self.set_tolerance(1, tol_c1);
        self.set_tolerance(2, tol_c2);
        self.my_dist11 = REAL_FIRST;
        self.my_dist12 = REAL_FIRST;
        self.my_dist21 = REAL_FIRST;
        self.my_dist22 = REAL_FIRST;
    }

    /// OCCT SetCurve(theRank, C) (cxx L138-143).
    pub fn set_curve(&mut self, the_rank: i32, c: &'a dyn ExtremaCurveTool) {
        if the_rank < 1 || the_rank > 2 {
            panic!("Standard_OutOfRange: Extrema_ExtCC::SetCurve()");
        }
        let an_ind = (the_rank - 1) as usize;
        if an_ind == 0 {
            self.my_c1 = Some(c);
        } else {
            self.my_c2 = Some(c);
        }
    }

    /// OCCT SetCurve(theRank, C, Uinf, Usup) (cxx L147-154).
    pub fn set_curve3(
        &mut self,
        the_rank: i32,
        c: &'a dyn ExtremaCurveTool,
        uinf: f64,
        usup: f64,
    ) {
        self.set_curve(the_rank, c);
        self.set_range(the_rank, uinf, usup);
    }

    /// OCCT SetRange(theRank, Uinf, Usup) (cxx L158-164).
    pub fn set_range(&mut self, the_rank: i32, uinf: f64, usup: f64) {
        if the_rank < 1 || the_rank > 2 {
            panic!("Standard_OutOfRange: Extrema_ExtCC::SetRange()");
        }
        let an_ind = (the_rank - 1) as usize;
        self.my_inf[an_ind] = uinf;
        self.my_sup[an_ind] = usup;
    }

    /// OCCT SetTolerance(theRank, theTol) (cxx L168-173).
    pub fn set_tolerance(&mut self, the_rank: i32, the_tol: f64) {
        if the_rank < 1 || the_rank > 2 {
            panic!("Standard_OutOfRange: Extrema_ExtCC::SetTolerance()");
        }
        let an_ind = (the_rank - 1) as usize;
        self.my_tol[an_ind] = the_tol;
    }

    /// OCCT SetSingleSolutionFlag(theFlag) (cxx L959-962).
    pub fn set_single_solution_flag(&mut self, the_flag: bool) {
        self.my_is_find_single_solution = the_flag;
    }

    /// OCCT GetSingleSolutionFlag() (cxx L966-969).
    pub fn get_single_solution_flag(&self) -> bool {
        self.my_is_find_single_solution
    }

    /// OCCT Perform() (cxx L177-317).
    pub fn perform(&mut self) {
        // cxx L179-180: Standard_NullObject_Raise_if(!myC[0] || !myC[1], ...).
        let c1 = self.my_c1.expect("Standard_NullObject: Extrema_ExtCC::Perform()");
        let c2 = self.my_c2.expect("Standard_NullObject: Extrema_ExtCC::Perform()");

        // cxx L186-191.
        self.my_done = false;
        self.my_points.clear();
        self.my_sq_dist.clear();
        self.my_is_parallel = false;

        // cxx L188-195.
        let type1 = c1.get_type();
        let type2 = c2.get_type();
        let tol = self.my_tol[0].min(self.my_tol[1]);

        let u11 = self.my_inf[0];
        let u12 = self.my_sup[0];
        let u21 = self.my_inf[1];
        let u22 = self.my_sup[1];

        // cxx L197-212.
        if !is_infinite_value(u11) {
            self.my_p1f = c1.value(u11);
        }
        if !is_infinite_value(u12) {
            self.my_p1l = c1.value(u12);
        }
        if !is_infinite_value(u21) {
            self.my_p2f = c2.value(u21);
        }
        if !is_infinite_value(u22) {
            self.my_p2l = c2.value(u22);
        }

        // cxx L214-245.
        self.my_dist11 = if is_infinite_value(u11) || is_infinite_value(u21) {
            REAL_LAST
        } else {
            (self.my_p1f - self.my_p2f).length_squared()
        };
        self.my_dist12 = if is_infinite_value(u11) || is_infinite_value(u22) {
            REAL_LAST
        } else {
            (self.my_p1f - self.my_p2l).length_squared()
        };
        self.my_dist21 = if is_infinite_value(u12) || is_infinite_value(u21) {
            REAL_LAST
        } else {
            (self.my_p1l - self.my_p2f).length_squared()
        };
        self.my_dist22 = if is_infinite_value(u12) || is_infinite_value(u22) {
            REAL_LAST
        } else {
            (self.my_p1l - self.my_p2l).length_squared()
        };

        // cxx L251-295: the analytical case (one curve is always a line).
        if (type1 == CurveType::Line && is_elementary(type2))
            || (type2 == CurveType::Line && is_elementary(type1))
        {
            let mut an_ind1 = 0usize;
            let mut an_ind2 = 1usize;
            let mut a_type2 = type2;
            let is_inverse = curve_rank(type1) > curve_rank(type2);
            if is_inverse {
                an_ind1 = 1;
                an_ind2 = 0;
                a_type2 = type1;
            }
            let (ca, cb) = if an_ind1 == 0 { (c1, c2) } else { (c2, c1) };
            match a_type2 {
                CurveType::Line => {
                    // Extrema_ExtElC Xtrem(myC[anInd1]->Line(), myC[anInd2]->Line(), Tol).
                    let a_l1 = ca.line();
                    let a_l2 = cb.line();
                    let xtrem = ExtremaExtElC::line_line(&a_l1, &a_l2, tol);
                    self.prepare_results_ext_el_c(&xtrem, is_inverse, u11, u12, u21, u22);
                }
                CurveType::Circle => {
                    // Extrema_ExtElC Xtrem(myC[anInd1]->Line(), myC[anInd2]->Circle(), Tol).
                    let a_l1 = ca.line();
                    let a_c2 = cb.circle();
                    let xtrem = ExtremaExtElC::line_circle(&a_l1, &a_c2, tol);
                    self.prepare_results_ext_el_c(&xtrem, is_inverse, u11, u12, u21, u22);
                }
                CurveType::Ellipse => {
                    let a_l1 = ca.line();
                    let an_e2 = cb.ellipse();
                    let xtrem = ExtremaExtElC::line_ellipse(&a_l1, &an_e2);
                    self.prepare_results_ext_el_c(&xtrem, is_inverse, u11, u12, u21, u22);
                }
                CurveType::Hyperbola => {
                    let a_l1 = ca.line();
                    let a_h2 = cb.hyperbola();
                    let xtrem = ExtremaExtElC::line_hyperbola(&a_l1, &a_h2);
                    self.prepare_results_ext_el_c(&xtrem, is_inverse, u11, u12, u21, u22);
                }
                CurveType::Parabola => {
                    let a_l1 = ca.line();
                    let a_p2 = cb.parabola();
                    let xtrem = ExtremaExtElC::line_parabola(&a_l1, &a_p2);
                    self.prepare_results_ext_el_c(&xtrem, is_inverse, u11, u12, u21, u22);
                }
                _ => {}
            }
        } else if type1 == CurveType::Circle && type2 == CurveType::Circle {
            // cxx L296-311: analytical case — two circles.
            let cc_xtrem = ExtremaExtElC::circle_circle(&c1.circle(), &c2.circle());
            let b_is_done = cc_xtrem.is_done();
            if b_is_done {
                self.prepare_results_ext_el_c(&cc_xtrem, false, u11, u12, u21, u22);
            } else {
                let mut my_ecc = GenExtCC::new(c1, c2);
                my_ecc.set_params(c1, c2, self.my_inf[0], self.my_sup[0], self.my_inf[1], self.my_sup[1]);
                my_ecc.set_tolerance(tol);
                my_ecc.set_single_solution_flag(self.get_single_solution_flag());
                my_ecc.perform();
                let is_parallel = my_ecc.is_parallel();
                if my_ecc.is_done() {
                    self.prepare_results_ecc(&my_ecc, is_parallel, u11, u12, u21, u22);
                }
            }
        } else {
            // cxx L312-316.
            let mut my_ecc = GenExtCC::new(c1, c2);
            my_ecc.set_params(
                c1,
                c2,
                self.my_inf[0],
                self.my_sup[0],
                self.my_inf[1],
                self.my_sup[1],
            );
            my_ecc.set_tolerance(tol);
            my_ecc.set_single_solution_flag(self.get_single_solution_flag());
            my_ecc.perform();
            let is_parallel = if my_ecc.is_done() {
                my_ecc.is_parallel()
            } else {
                false
            };
            if my_ecc.is_done() {
                self.prepare_results_ecc(&my_ecc, is_parallel, u11, u12, u21, u22);
            }
        }
    }

    /// OCCT IsDone() (cxx L321-324).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT IsParallel() (cxx L328-336).
    pub fn is_parallel(&self) -> bool {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_is_parallel
    }

    /// OCCT SquareDistance(N) (cxx L340-347) — 1-based.
    pub fn square_distance(&self, n: usize) -> f64 {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.my_sq_dist[n - 1]
    }

    /// OCCT NbExt() (cxx L351-358).
    pub fn nb_ext(&self) -> usize {
        if !self.my_done {
            panic!("StdFail_NotDone");
        }
        self.my_sq_dist.len()
    }

    /// OCCT Points(N, P1, P2) (cxx L362-371) — 1-based.
    pub fn points(&self, n: usize, p1: &mut POnCurve, p2: &mut POnCurve) {
        if n < 1 || n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        *p1 = self.my_points[2 * n - 2].clone();
        *p2 = self.my_points[2 * n - 1].clone();
    }

    /// OCCT TrimmedSquareDistances (cxx L375-393) — the four corner distances.
    #[allow(clippy::type_complexity)]
    pub fn trimmed_square_distances(&self) -> (f64, f64, f64, f64, DVec3, DVec3, DVec3, DVec3) {
        (
            self.my_dist11,
            self.my_dist12,
            self.my_dist21,
            self.my_dist22,
            self.my_p1f,
            self.my_p1l,
            self.my_p2f,
            self.my_p2l,
        )
    }

    /// OCCT PrepareParallelResult(theUt11, theUt12, theUt21, theUt22, theSqDist)
    /// (cxx L397-828).
    fn prepare_parallel_result(
        &mut self,
        the_ut11: f64,
        the_ut12: f64,
        the_ut21: f64,
        the_ut22: f64,
        the_sq_dist: f64,
    ) {
        if !self.my_is_parallel {
            return;
        }

        let c1 = self.my_c1.expect("PrepareParallelResult");
        let c2 = self.my_c2.expect("PrepareParallelResult");

        // cxx L408-417.
        let a_type1 = c1.get_type();
        let a_type2 = c2.get_type();
        if ((a_type1 != CurveType::Line) && (a_type1 != CurveType::Circle))
            || ((a_type2 != CurveType::Line) && (a_type2 != CurveType::Circle))
        {
            self.my_sq_dist.push(the_sq_dist);
            self.my_done = true;
            self.my_is_parallel = true;
            return;
        }

        // cxx L421-445: line-circle / circle-line.
        if a_type1 != a_type2 {
            let is_reversed = a_type1 != CurveType::Circle;
            let a_p_on_c = if !is_reversed {
                c1.value(the_ut11)
            } else {
                c2.value(the_ut21)
            };

            let a_l = if !is_reversed { c2.line() } else { c1.line() };
            let ext_p_lin = ExtremaExtPElC::point_line(
                a_p_on_c,
                &a_l,
                CONFUSION,
                if !is_reversed { the_ut21 } else { the_ut11 },
                if !is_reversed { the_ut22 } else { the_ut12 },
            );

            if ext_p_lin.is_done() {
                self.my_sq_dist.push(the_sq_dist);
            } else {
                self.my_is_parallel = false;
            }
            return;
        }

        if a_type1 == CurveType::Line {
            // cxx L447-590: Line - Line.
            let is_first_infinite =
                is_infinite_value(the_ut11) && is_infinite_value(the_ut12);
            let is_last_infinite = is_infinite_value(the_ut21) && is_infinite_value(the_ut22);

            if is_first_infinite || is_last_infinite {
                // Infinite number of solutions.
                self.my_sq_dist.push(the_sq_dist);
            } else {
                // cxx L468-588.
                self.my_is_parallel = false;

                let a_lin1 = c1.line();
                let a_lin2 = c2.line();
                let is_opposite = a_lin1.direction.dot(a_lin2.direction) < 0.0;

                let mut a_range2 = Range::from_bounds(the_ut21, the_ut22);
                let mut a_proj_rng12 = Range::new();

                if is_infinite_value(the_ut11) {
                    if is_opposite {
                        a_proj_rng12.add_parameter(crate::core::precision::INFINITE_VALUE);
                    } else {
                        a_proj_rng12.add_parameter(-crate::core::precision::INFINITE_VALUE);
                    }
                } else {
                    let a_p_on_c1 = elclib_line_value(the_ut11, &a_lin1);
                    let a_par = elclib_line_parameter(&a_lin2, a_p_on_c1);
                    a_proj_rng12.add_parameter(a_par);
                }

                if is_infinite_value(the_ut12) {
                    if is_opposite {
                        a_proj_rng12.add_parameter(-crate::core::precision::INFINITE_VALUE);
                    } else {
                        a_proj_rng12.add_parameter(crate::core::precision::INFINITE_VALUE);
                    }
                } else {
                    let a_p_on_c1 = elclib_line_value(the_ut12, &a_lin1);
                    let a_par = elclib_line_parameter(&a_lin2, a_p_on_c1);
                    a_proj_rng12.add_parameter(a_par);
                }

                a_range2.common(&a_proj_rng12);
                if a_range2.delta() > CONFUSION {
                    self.clear_solutions();
                    self.my_sq_dist.push(the_sq_dist);
                    self.my_is_parallel = true;
                } else if !a_range2.is_void() {
                    // cxx L520-541.
                    self.clear_solutions();
                    let mut a_par1 = 0.0f64;
                    let mut a_par2 = 0.0f64;
                    if let Some((lo, hi)) = a_range2.get_bounds() {
                        a_par1 = lo;
                        a_par2 = hi;
                    }
                    a_par2 = 0.5 * (a_par1 + a_par2);
                    let mut a_p = elclib_line_value(a_par2, &a_lin2);
                    let a_p2 = POnCurve {
                        param: a_par2,
                        point: a_p,
                    };
                    a_par1 = elclib_line_parameter(&a_lin1, a_p);
                    a_p = elclib_line_value(a_par1, &a_lin1);
                    let a_p1 = POnCurve {
                        param: a_par1,
                        point: a_p,
                    };
                    self.my_points.push(a_p1);
                    self.my_points.push(a_p2);
                    self.my_sq_dist.push(the_sq_dist);
                } else {
                    // cxx L542-588.
                    let a_dists = [self.my_dist11, self.my_dist12, self.my_dist21, self.my_dist22];
                    let mut a_dmin = a_dists[0];
                    let mut imin = 0usize;
                    for i in 1..4 {
                        if a_dmin > a_dists[i] {
                            a_dmin = a_dists[i];
                            imin = i;
                        }
                    }
                    let a_p1;
                    let a_p2;
                    if imin == 0 {
                        a_p1 = POnCurve {
                            param: self.my_inf[0],
                            point: self.my_p1f,
                        };
                        a_p2 = POnCurve {
                            param: self.my_inf[1],
                            point: self.my_p2f,
                        };
                    } else if imin == 1 {
                        a_p1 = POnCurve {
                            param: self.my_inf[0],
                            point: self.my_p1f,
                        };
                        a_p2 = POnCurve {
                            param: self.my_sup[1],
                            point: self.my_p2l,
                        };
                    } else if imin == 2 {
                        a_p1 = POnCurve {
                            param: self.my_sup[0],
                            point: self.my_p1l,
                        };
                        a_p2 = POnCurve {
                            param: self.my_inf[1],
                            point: self.my_p2f,
                        };
                    } else {
                        a_p1 = POnCurve {
                            param: self.my_sup[0],
                            point: self.my_p1l,
                        };
                        a_p2 = POnCurve {
                            param: self.my_sup[1],
                            point: self.my_p2l,
                        };
                    }
                    self.clear_solutions();
                    self.my_points.push(a_p1);
                    self.my_points.push(a_p2);
                    self.my_sq_dist.push(a_dmin);
                }
            }
        } else {
            // cxx L591-827: Circle - Circle.
            self.my_is_parallel = false;

            let a_work_circ = c2.circle();
            let a_period = std::f64::consts::PI + std::f64::consts::PI;
            let a_p12 = c1.value(the_ut12);
            let (a_p11, a_v_tg1) = c1.d1(the_ut11);

            let mut a_range = Range::from_bounds(the_ut21, the_ut22);
            let mut a_proj_rng1 = Range::new();

            // cxx L614-615: aPrecision = Max(Epsilon(R1), Epsilon(R2)).
            let a_precision =
                epsilon_of(c1.circle().radius).max(epsilon_of(c2.circle().radius));

            // cxx L620-621.
            let mut a_par1 = elclib_in_period(
                elclib_circle_parameter(&a_work_circ, a_p11),
                the_ut21,
                the_ut21 + a_period,
            );
            let a_v_tg2 = c2.dn(a_par1, 1);

            // cxx L624: same/opposite directions.
            let is_opposite = a_v_tg1.dot(a_v_tg2) < 0.0;

            let mut a_par2 = elclib_in_period(
                elclib_circle_parameter(&a_work_circ, a_p12),
                the_ut21,
                the_ut21 + a_period,
            );

            // cxx L628-643.
            if is_opposite {
                if a_range.delta() > ANGULAR && (a_par1 - a_par2) < ANGULAR {
                    a_par2 -= a_period;
                }
            } else if a_range.delta() > ANGULAR && (a_par2 - a_par1) < ANGULAR {
                a_par2 += a_period;
            }

            // cxx L648-826.
            let mut a_min_square_dist = REAL_LAST;

            a_proj_rng1.add_parameter(a_par1 - a_period);
            a_proj_rng1.add_parameter(a_par2 - a_period);
            for _i in 0..3 {
                let mut a_rng = a_proj_rng1.clone();
                a_rng.common(&a_range);

                if a_rng.delta() > CONFUSION {
                    // cxx L675-719.
                    let a_par = a_rng.get_intermediate_point(0.5).unwrap_or(0.0);
                    let a_p_circ2 = elclib_circle_value(a_par, &a_work_circ);
                    let ext_p_cir = ExtremaExtPElC::point_circle(
                        a_p_circ2,
                        &c1.circle(),
                        CONFUSION,
                        the_ut11,
                        the_ut12,
                    );
                    if ext_p_cir.nb_ext() < 1 {
                        a_par1_reset(&mut a_par1, a_par);
                        a_proj_rng1.shift(std::f64::consts::PI);
                        continue;
                    }
                    let mut a_min_sq_d = ext_p_cir.square_distance(1);
                    for an_ext_id in 2..=ext_p_cir.nb_ext() {
                        a_min_sq_d = a_min_sq_d.min(ext_p_cir.square_distance(an_ext_id));
                    }

                    if a_min_sq_d <= a_min_square_dist + (1.0 + a_min_sq_d) * a_precision {
                        self.clear_solutions();
                        self.my_sq_dist.push(a_min_sq_d);
                        self.my_is_parallel = true;

                        let a_delta_sq_dist = a_min_sq_d - the_sq_dist;
                        let a_sq_d = a_min_sq_d.max(the_sq_dist);

                        if a_delta_sq_dist * a_delta_sq_dist
                            < 4.0 * a_sq_d * SQUARE_CONFUSION
                        {
                            break;
                        }
                    }
                } else if !a_rng.is_void() {
                    // cxx L721-772.
                    let a_par = a_rng.get_intermediate_point(0.5).unwrap_or(0.0);
                    let a_p_circ2 = elclib_circle_value(a_par, &a_work_circ);
                    let a_p2 = POnCurve {
                        param: a_par,
                        point: a_p_circ2,
                    };

                    let ext_p_cir = ExtremaExtPElC::point_circle(
                        a_p_circ2,
                        &c1.circle(),
                        CONFUSION,
                        the_ut11,
                        the_ut12,
                    );

                    let mut is_found = !self.my_is_parallel;

                    if !is_found {
                        for an_ext_id in 1..=ext_p_cir.nb_ext() {
                            if ext_p_cir.square_distance(an_ext_id) < a_min_square_dist {
                                is_found = true;
                                break;
                            }
                        }
                    }

                    if is_found {
                        self.clear_solutions();
                        self.my_is_parallel = false;
                        for an_ext_id in 1..=ext_p_cir.nb_ext() {
                            self.my_points.push(ext_p_cir.point(an_ext_id).clone());
                            self.my_points.push(a_p2.clone());
                            self.my_sq_dist.push(ext_p_cir.square_distance(an_ext_id));
                            a_min_square_dist =
                                a_min_square_dist.min(ext_p_cir.square_distance(an_ext_id));
                        }
                    }
                } else {
                    // cxx L773-824.
                    self.my_is_parallel = false;
                    let a_dists =
                        [self.my_dist11, self.my_dist12, self.my_dist21, self.my_dist22];
                    let mut a_dmin = a_dists[0];
                    let mut imin = 0usize;
                    for k in 1..4 {
                        if a_dmin > a_dists[k] {
                            a_dmin = a_dists[k];
                            imin = k;
                        }
                    }
                    if a_dmin <= a_min_square_dist + (1.0 + a_dmin) * a_precision {
                        let a_p1;
                        let a_p2;
                        if imin == 0 {
                            a_p1 = POnCurve {
                                param: self.my_inf[0],
                                point: self.my_p1f,
                            };
                            a_p2 = POnCurve {
                                param: self.my_inf[1],
                                point: self.my_p2f,
                            };
                        } else if imin == 1 {
                            a_p1 = POnCurve {
                                param: self.my_inf[0],
                                point: self.my_p1f,
                            };
                            a_p2 = POnCurve {
                                param: self.my_sup[1],
                                point: self.my_p2l,
                            };
                        } else if imin == 2 {
                            a_p1 = POnCurve {
                                param: self.my_sup[0],
                                point: self.my_p1l,
                            };
                            a_p2 = POnCurve {
                                param: self.my_inf[1],
                                point: self.my_p2f,
                            };
                        } else {
                            a_p1 = POnCurve {
                                param: self.my_sup[0],
                                point: self.my_p1l,
                            };
                            a_p2 = POnCurve {
                                param: self.my_sup[1],
                                point: self.my_p2l,
                            };
                        }
                        self.clear_solutions();
                        self.my_points.push(a_p1);
                        self.my_points.push(a_p2);
                        self.my_sq_dist.push(a_dmin);
                        a_min_square_dist = a_min_square_dist.min(a_dmin);
                    }
                }
                a_proj_rng1.shift(std::f64::consts::PI);
            }
        }
    }

    /// OCCT PrepareResults(const Extrema_ExtElC& AlgExt, theIsInverse, ...)
    /// (cxx L832-901).
    fn prepare_results_ext_el_c(
        &mut self,
        alg_ext: &ExtremaExtElC,
        the_is_inverse: bool,
        ut11: f64,
        ut12: f64,
        ut21: f64,
        ut22: f64,
    ) {
        let c1 = self.my_c1.expect("PrepareResults");
        let c2 = self.my_c2.expect("PrepareResults");

        self.my_done = alg_ext.is_done();
        if self.my_done {
            self.my_is_parallel = alg_ext.is_parallel();
            if self.my_is_parallel {
                self.prepare_parallel_result(ut11, ut12, ut21, ut22, alg_ext.square_distance(1));
            } else {
                let nb_ext = alg_ext.nb_ext();
                for i in 1..=nb_ext {
                    // cxx L857-867: verification of the validity of parameters.
                    let mut p1 = POnCurve {
                        param: 0.0,
                        point: DVec3::ZERO,
                    };
                    let mut p2 = POnCurve {
                        param: 0.0,
                        point: DVec3::ZERO,
                    };
                    alg_ext.points(i, &mut p1, &mut p2);
                    let mut u;
                    let mut u2;
                    if !the_is_inverse {
                        u = p1.param;
                        u2 = p2.param;
                    } else {
                        u2 = p1.param;
                        u = p2.param;
                    }

                    // cxx L869-876.
                    if c1.is_periodic() {
                        u = elclib_in_period(u, ut11, ut11 + c1.period());
                    }
                    if c2.is_periodic() {
                        u2 = elclib_in_period(u2, ut21, ut21 + c2.period());
                    }

                    // cxx L878-897.
                    if (u >= ut11 - f64::EPSILON)
                        && (u <= ut12 + f64::EPSILON)
                        && (u2 >= ut21 - f64::EPSILON)
                        && (u2 <= ut22 + f64::EPSILON)
                    {
                        let val = alg_ext.square_distance(i);
                        self.my_sq_dist.push(val);
                        if !the_is_inverse {
                            p1.param = u;
                            p2.param = u2;
                            self.my_points.push(p1);
                            self.my_points.push(p2);
                        } else {
                            p1.param = u2;
                            p2.param = u;
                            self.my_points.push(p2);
                            self.my_points.push(p1);
                        }
                    }
                }
            }
        }
    }

    /// OCCT PrepareResults(const Extrema_ECC& AlgExt, ...) (cxx L905-955).
    fn prepare_results_ecc(
        &mut self,
        alg_ext: &GenExtCC<'_>,
        is_parallel: bool,
        ut11: f64,
        ut12: f64,
        ut21: f64,
        ut22: f64,
    ) {
        let c1 = self.my_c1.expect("PrepareResults");
        let c2 = self.my_c2.expect("PrepareResults");

        self.my_done = alg_ext.is_done();
        if !self.my_done {
            return;
        }
        self.my_is_parallel = is_parallel;
        if self.my_is_parallel {
            self.prepare_parallel_result(ut11, ut12, ut21, ut22, alg_ext.square_distance(1));
        } else {
            let nb_ext = alg_ext.nb_ext();
            for i in 1..=nb_ext {
                let mut p1 = POnCurve {
                    param: 0.0,
                    point: DVec3::ZERO,
                };
                let mut p2 = POnCurve {
                    param: 0.0,
                    point: DVec3::ZERO,
                };
                alg_ext.points(i, &mut p1, &mut p2);
                let mut u = p1.param;
                let mut u2 = p2.param;

                // cxx L933-940.
                if c1.is_periodic() {
                    u = elclib_in_period(u, ut11, ut11 + c1.period());
                }
                if c2.is_periodic() {
                    u2 = elclib_in_period(u2, ut21, ut21 + c2.period());
                }

                // cxx L942-951.
                if (u >= ut11 - f64::EPSILON)
                    && (u <= ut12 + f64::EPSILON)
                    && (u2 >= ut21 - f64::EPSILON)
                    && (u2 <= ut22 + f64::EPSILON)
                {
                    let val = alg_ext.square_distance(i);
                    self.my_sq_dist.push(val);
                    p1.param = u;
                    p2.param = u2;
                    self.my_points.push(p1);
                    self.my_points.push(p2);
                }
            }
        }
    }

    /// OCCT ClearSolutions() (hxx L146-150) — clears the results without
    /// touching the Done / IsParallel flags.
    fn clear_solutions(&mut self) {
        self.my_sq_dist.clear();
        self.my_points.clear();
    }
}

/// The OCCT GeomAbs_CurveType ordering used by the analytical dispatch
/// (cxx L251-252): `type2 <= GeomAbs_Parabola`.
fn is_elementary(t: CurveType) -> bool {
    matches!(
        t,
        CurveType::Line
            | CurveType::Circle
            | CurveType::Ellipse
            | CurveType::Hyperbola
            | CurveType::Parabola
    )
}

/// The OCCT `GeomAbs_CurveType` ordinal of the elementary kinds — the rcad
/// `CurveType` keeps the same relative order for Line..Parabola.
fn curve_rank(t: CurveType) -> i32 {
    match t {
        CurveType::Line => 0,
        CurveType::Circle => 1,
        CurveType::Ellipse => 2,
        CurveType::Hyperbola => 3,
        CurveType::Parabola => 4,
        CurveType::Bezier => 5,
        CurveType::BSpline => 6,
        CurveType::Other => 8,
    }
}

/// The no-op helper keeping the `continue` arm of the circle-circle loop shaped
/// like the OCCT `continue` (the OCCT body re-derives `aPar` on the next pass,
/// so the binding is unused).
fn a_par1_reset(_a_par1: &mut f64, _a_par: f64) {}

/// The `Circle3`/`Line3` reuse markers keep the unused-import surface honest.
#[allow(dead_code)]
fn _unused_types(_l: Line3, _c: Circle3, _s: GeomAbsShape) {}
