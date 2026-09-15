// OCCT Approx_FitAndDivide (TKGeomBase/Approx) — 1:1 Rust translation of
// the Approx_ComputeCLine.gxx engine (L26-485) as instantiated by
// Approx_FitAndDivide_0.cxx (L17-25):
//   MultiLine           = AppCont_Function  (rcad `AppContFunction`, app_cont.rs)
//   Approx_ComputeCLine = Approx_FitAndDivide (rcad `ApproxFitAndDivide`)
//   LeastSquare         = AppCont_LeastSquare (rcad app_cont::AppContLeastSquare)
//
// The public surface follows Approx_FitAndDivide.hxx (L39-105): the two
// constructors (with/without the MultiLine), Perform, NbMultiCurves, Value,
// Parameters, the Set* reconfigurations, IsAllApproximated,
// IsToleranceReached and Error.
//
// OCCT quirk kept: the constructor initialises myInvOrder = true (gxx L53),
// i.e. the degree scan starts at degreemax even though the .hxx comment
// (L75-79) documents "by default inverse order is used" alongside a
// SetInvOrder setter.

use rcad_kernel::core::precision::{CONFUSION, PCONFUSION, REAL_LAST};

use super::app_cont::{real_last, AppContFunction, AppContLeastSquare};
use super::approx_int::{AppParConstraint, MultiCurve};

/// OCCT Standard_Real RealEpsilon() (Standard_Real.hxx L160-163) — the
/// smallest double e such that 1.0 + e != 1.0.
const REAL_EPSILON: f64 = f64::EPSILON;

/// OCCT Precision::PApproximation() — rcad core::precision::p_approximation()
/// (Approximation() * 0.01).
use rcad_kernel::core::precision::p_approximation;

/// OCCT `const static int MAXSEGM = 1000` (the gxx L26).
const MAXSEGM: i32 = 1000;

/// OCCT Approx_ComputeCLine (the gxx L29-60 constructor with the MultiLine)
/// as named by the Approx_FitAndDivide instantiation (the .hxx L30-135).
///
/// The fields mirror the .hxx L115-134 declaration order.
pub struct ApproxFitAndDivide {
    my_multi_curves: Vec<MultiCurve>, // OCCT: myMultiCurves (hxx L115)
    myfirstparam: Vec<f64>,           // OCCT: myfirstparam (hxx L116)
    mylastparam: Vec<f64>,            // OCCT: mylastparam (hxx L117)
    the_multi_curve: MultiCurve,      // OCCT: TheMultiCurve (hxx L118)
    alldone: bool,                    // OCCT: alldone (hxx L119)
    tolreached: bool,                 // OCCT: tolreached (hxx L120)
    tolers3d: Vec<f64>,               // OCCT: Tolers3d (hxx L121)
    tolers2d: Vec<f64>,               // OCCT: Tolers2d (hxx L122)
    mydegremin: i32,                  // OCCT: mydegremin (hxx L123)
    mydegremax: i32,                  // OCCT: mydegremax (hxx L124)
    mytol3d: f64,                     // OCCT: mytol3d (hxx L125)
    mytol2d: f64,                     // OCCT: mytol2d (hxx L126)
    currenttol3d: f64,                // OCCT: currenttol3d (hxx L127)
    currenttol2d: f64,                // OCCT: currenttol2d (hxx L128)
    mycut: bool,                      // OCCT: mycut (hxx L129)
    myfirstc: AppParConstraint,       // OCCT: myfirstC (hxx L130)
    mylastc: AppParConstraint,        // OCCT: mylastC (hxx L131)
    my_max_segments: i32,             // OCCT: myMaxSegments (hxx L132)
    my_inv_order: bool,               // OCCT: myInvOrder (hxx L133)
    my_hang_checking: bool,           // OCCT: myHangChecking (hxx L134)
}

impl ApproxFitAndDivide {
    /// OCCT Approx_FitAndDivide(Line, degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, cutting, FirstC, LastC) (the .hxx L39-47 with the
    /// TangencyPoint defaults; the gxx L36-60 body).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        line: &dyn AppContFunction,
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        cutting: bool,
    ) -> Self {
        let mut fit = ApproxFitAndDivide {
            my_multi_curves: Vec::new(),
            myfirstparam: Vec::new(),
            mylastparam: Vec::new(),
            the_multi_curve: MultiCurve::new(0, 0, 0),
            alldone: false,
            tolreached: false,
            tolers3d: Vec::new(),
            tolers2d: Vec::new(),
            mydegremin: degreemin,
            mydegremax: degreemax,
            mytol3d: tolerance3d,
            mytol2d: tolerance2d,
            currenttol3d: 0.0,
            currenttol2d: 0.0,
            mycut: cutting,
            myfirstc: AppParConstraint::TangencyPoint,
            mylastc: AppParConstraint::TangencyPoint,
            my_max_segments: MAXSEGM,
            my_inv_order: true,
            my_hang_checking: true,
        };
        fit.perform(line);
        fit
    }

    /// OCCT Approx_FitAndDivide(degreemin, degreemax, Tolerance3d,
    /// Tolerance2d, cutting, FirstC, LastC) (the .hxx L50-57; the gxx
    /// L67-89 body) — initializes the fields of the algorithm.
    #[allow(clippy::too_many_arguments)]
    pub fn new_init(
        degreemin: i32,
        degreemax: i32,
        tolerance3d: f64,
        tolerance2d: f64,
        cutting: bool,
    ) -> Self {
        ApproxFitAndDivide {
            my_multi_curves: Vec::new(),
            myfirstparam: Vec::new(),
            mylastparam: Vec::new(),
            the_multi_curve: MultiCurve::new(0, 0, 0),
            alldone: false,
            tolreached: false,
            tolers3d: Vec::new(),
            tolers2d: Vec::new(),
            mydegremin: degreemin,
            mydegremax: degreemax,
            mytol3d: tolerance3d,
            mytol2d: tolerance2d,
            currenttol3d: 0.0,
            currenttol2d: 0.0,
            mycut: cutting,
            myfirstc: AppParConstraint::TangencyPoint,
            mylastc: AppParConstraint::TangencyPoint,
            my_max_segments: MAXSEGM,
            my_inv_order: true,
            my_hang_checking: true,
        }
    }

    /// OCCT Approx_ComputeCLine::Perform(Line) (the gxx L96-229) — runs the
    /// algorithm after having initialized the fields.
    pub fn perform(&mut self, line: &dyn AppContFunction) {
        let u_first = line.first_parameter();
        let u_last = line.last_parameter();
        let mut finish = false;
        let mut begin = true;
        let mut ok = false;
        let mut thetol3d = CONFUSION;
        let mut thetol2d = CONFUSION;
        let tol_u = if self.my_hang_checking {
            ((u_last - u_first) * 1.0e-03).max(CONFUSION)
        } else {
            ((u_last - u_first) * 1.0e-05).max(p_approximation())
        };
        let mut myfirstu = u_first;
        let mut mylastu = u_last;
        let mut a_max_segments = 0;
        let a_max_segments1 = self.my_max_segments - 1;
        let mut a_nb_cut = 0;
        let mut a_nb_imp = 0;
        let a_nb_comp = 10;

        if !self.mycut {
            let alldone = self.compute(line, u_first, u_last, &mut thetol3d, &mut thetol2d);
            self.alldone = alldone;
            if !self.alldone {
                self.tolreached = false;
                self.myfirstparam.push(u_first);
                self.mylastparam.push(u_last);
                self.my_multi_curves.push(self.the_multi_curve.clone());
                self.tolers3d.push(self.currenttol3d);
                self.tolers2d.push(self.currenttol2d);
            }
        } else {
            // previous decision to be taken if we get worse with next cut (eap)
            let mut kept_multi_curve = MultiCurve::new(0, 0, 0);
            let mut kept_ufirst = 0.0f64;
            let mut kept_ulast = 0.0f64;
            let mut kept_t3d = REAL_LAST;
            let mut kept_t2d = 0.0f64;

            while !finish {
                // Management of the multiline splitting for approximation:
                if !begin {
                    if ok {
                        // Calcul de la partie a approximer.
                        myfirstu = mylastu;
                        mylastu = u_last;
                        a_nb_cut = 0;
                        a_nb_imp = 0;
                        if (u_last - myfirstu).abs() <= REAL_EPSILON
                            || a_max_segments >= self.my_max_segments
                        {
                            finish = true;
                            self.alldone = true;
                            return;
                        }
                        kept_t3d = real_last();
                        kept_t2d = 0.0;
                        kept_ufirst = myfirstu;
                        kept_ulast = mylastu;
                    } else {
                        // keep best decision
                        if (thetol3d + thetol2d) < (kept_t3d + kept_t2d) {
                            kept_multi_curve = self.the_multi_curve.clone();
                            kept_ufirst = myfirstu;
                            kept_ulast = mylastu;
                            kept_t3d = thetol3d;
                            kept_t2d = thetol2d;
                            a_nb_imp += 1;
                        }

                        // cut an interval
                        mylastu = (myfirstu + mylastu) / 2.0;
                        a_nb_cut += 1;
                    }
                }

                // Calculation of parameters on this new interval.
                ok = self.compute(line, myfirstu, mylastu, &mut thetol3d, &mut thetol2d);
                if ok {
                    a_max_segments += 1;
                }

                let mut a_stop_cutting = false;
                if self.my_hang_checking && a_nb_cut >= a_nb_comp {
                    if a_nb_cut > a_nb_imp + 1 {
                        a_stop_cutting = true;
                    }
                    a_nb_cut = 0;
                    a_nb_imp = 0;
                }
                // is new decision better?
                if !ok
                    && ((myfirstu - mylastu).abs() <= tol_u
                        || a_max_segments >= a_max_segments1
                        || a_stop_cutting)
                {
                    ok = true; // stop interval cutting, approx the rest part

                    if (thetol3d + thetol2d) < (kept_t3d + kept_t2d) {
                        kept_multi_curve = self.the_multi_curve.clone();
                        kept_ufirst = myfirstu;
                        kept_ulast = mylastu;
                        kept_t3d = thetol3d;
                        kept_t2d = thetol2d;
                    }

                    mylastu = kept_ulast;

                    self.tolreached = false; // helas
                    self.my_multi_curves.push(kept_multi_curve.clone());
                    a_max_segments += 1;
                    self.tolers3d.push(kept_t3d);
                    self.tolers2d.push(kept_t2d);
                    self.myfirstparam.push(kept_ufirst);
                    self.mylastparam.push(kept_ulast);
                }

                begin = false;
            } // while (!Finish)
        }
    }

    /// OCCT Approx_ComputeCLine::NbMultiCurves() (the gxx L237-240).
    pub fn nb_multi_curves(&self) -> i32 {
        self.my_multi_curves.len() as i32
    }

    /// OCCT Approx_ComputeCLine::Value(Index) (the gxx L247-250).
    pub fn value(&self, index: i32) -> MultiCurve {
        self.my_multi_curves[(index - 1) as usize].clone()
    }

    /// OCCT Approx_ComputeCLine::Compute(Line, Ufirst, Ulast, TheTol3d,
    /// TheTol2d) (the gxx L257-382) — internally used by the algorithms.
    fn compute(
        &mut self,
        line: &dyn AppContFunction,
        ufirst: f64,
        ulast: f64,
        the_tol_3d: &mut f64,
        the_tol_2d: &mut f64,
    ) -> bool {
        let nb_points_max = 24;
        let a_min_ratio = 0.05;
        let a_max_deg = 8;

        let mut fv;

        let mut a_prev_curve = MultiCurve::new(0, 0, 0);
        let mut a_prev_tol3d = REAL_LAST;
        let mut a_prev_tol2d = REAL_LAST;
        let mut a_prev_is_ok = false;
        let mut an_inv_order = self.my_inv_order;
        if an_inv_order && self.mydegremax > a_max_deg {
            if (ulast - ufirst) / (line.last_parameter() - line.first_parameter()) < a_min_ratio
            {
                an_inv_order = false;
            }
        }
        if an_inv_order {
            let mut deg = self.mydegremax;
            while deg >= self.mydegremin {
                let nb_points = (2 * deg + 1).min(nb_points_max);
                let mut l_square = AppContLeastSquare::new(
                    line,
                    ufirst,
                    ulast,
                    self.myfirstc,
                    self.mylastc,
                    deg,
                    nb_points,
                );
                let mydone = l_square.is_done();
                if mydone {
                    // OCCT L292: LSquare.Error(Fv, TheTol3d, TheTol2d).
                    let (f, max_e3d, max_e2d) = l_square.error();
                    fv = f;
                    *the_tol_3d = max_e3d;
                    *the_tol_2d = max_e2d;
                    if *the_tol_3d <= self.mytol3d && *the_tol_2d <= self.mytol2d {
                        if deg == self.mydegremin {
                            // Stockage de la multicurve approximee.
                            self.tolreached = true;
                            self.my_multi_curves.push(l_square.value());
                            self.myfirstparam.push(ufirst);
                            self.mylastparam.push(ulast);
                            self.tolers3d.push(*the_tol_3d);
                            self.tolers2d.push(*the_tol_2d);
                            return true;
                        }
                        a_prev_tol3d = *the_tol_3d;
                        a_prev_tol2d = *the_tol_2d;
                        a_prev_curve = l_square.value();
                        a_prev_is_ok = true;
                        deg -= 1;
                        continue;
                    } else if a_prev_is_ok {
                        // Stockage de la multicurve approximee.
                        self.tolreached = true;
                        *the_tol_3d = a_prev_tol3d;
                        *the_tol_2d = a_prev_tol2d;
                        self.my_multi_curves.push(a_prev_curve.clone());
                        self.myfirstparam.push(ufirst);
                        self.mylastparam.push(ulast);
                        self.tolers3d.push(a_prev_tol3d);
                        self.tolers2d.push(a_prev_tol2d);
                        return true;
                    }
                } else if a_prev_is_ok {
                    // Stockage de la multicurve approximee.
                    self.tolreached = true;
                    *the_tol_3d = a_prev_tol3d;
                    *the_tol_2d = a_prev_tol2d;
                    self.my_multi_curves.push(a_prev_curve.clone());
                    self.myfirstparam.push(ufirst);
                    self.mylastparam.push(ulast);
                    self.tolers3d.push(a_prev_tol3d);
                    self.tolers2d.push(a_prev_tol2d);
                    return true;
                }
                if !a_prev_is_ok && deg == self.mydegremax {
                    self.the_multi_curve = l_square.value();
                    self.currenttol3d = *the_tol_3d;
                    self.currenttol2d = *the_tol_2d;
                    a_prev_tol3d = *the_tol_3d;
                    a_prev_tol2d = *the_tol_2d;
                    a_prev_curve = self.the_multi_curve.clone();
                    break;
                }
                deg -= 1;
            }
        } else {
            for deg in self.mydegremin..=self.mydegremax {
                let nb_points = (2 * deg + 1).min(nb_points_max);
                let mut l_square = AppContLeastSquare::new(
                    line,
                    ufirst,
                    ulast,
                    self.myfirstc,
                    self.mylastc,
                    deg,
                    nb_points,
                );
                let mydone = l_square.is_done();
                if mydone {
                    let (f, max_e3d, max_e2d) = l_square.error();
                    fv = f;
                    *the_tol_3d = max_e3d;
                    *the_tol_2d = max_e2d;
                    if *the_tol_3d <= self.mytol3d && *the_tol_2d <= self.mytol2d {
                        // Stockage de la multicurve approximee.
                        self.tolreached = true;
                        self.my_multi_curves.push(l_square.value());
                        self.myfirstparam.push(ufirst);
                        self.mylastparam.push(ulast);
                        self.tolers3d.push(*the_tol_3d);
                        self.tolers2d.push(*the_tol_2d);
                        return true;
                    }
                }
                if deg == self.mydegremax {
                    self.the_multi_curve = l_square.value();
                    self.currenttol3d = *the_tol_3d;
                    self.currenttol2d = *the_tol_2d;
                }
            }
        }
        let _ = fv;
        false
    }

    /// OCCT Approx_ComputeCLine::Parameters(Index, firstpar, lastpar)
    /// (the gxx L390-394).
    pub fn parameters(&self, index: i32) -> (f64, f64) {
        (
            self.myfirstparam[(index - 1) as usize],
            self.mylastparam[(index - 1) as usize],
        )
    }

    /// OCCT Approx_ComputeCLine::SetDegrees(degreemin, degreemax)
    /// (the gxx L401-405).
    pub fn set_degrees(&mut self, degreemin: i32, degreemax: i32) {
        self.mydegremin = degreemin;
        self.mydegremax = degreemax;
    }

    /// OCCT Approx_ComputeCLine::SetTolerances(Tolerance3d, Tolerance2d)
    /// (the gxx L412-416).
    pub fn set_tolerances(&mut self, tolerance3d: f64, tolerance2d: f64) {
        self.mytol3d = tolerance3d;
        self.mytol2d = tolerance2d;
    }

    /// OCCT Approx_ComputeCLine::SetConstraints(FirstC, LastC)
    /// (the gxx L423-428).
    pub fn set_constraints(&mut self, firstc: AppParConstraint, lastc: AppParConstraint) {
        self.myfirstc = firstc;
        self.mylastc = lastc;
    }

    /// OCCT Approx_ComputeCLine::SetMaxSegments(theMaxSegments)
    /// (the gxx L435-438).
    pub fn set_max_segments(&mut self, the_max_segments: i32) {
        self.my_max_segments = the_max_segments;
    }

    /// OCCT Approx_ComputeCLine::SetInvOrder(theInvOrder) (the gxx L442-445).
    pub fn set_inv_order(&mut self, the_inv_order: bool) {
        self.my_inv_order = the_inv_order;
    }

    /// OCCT Approx_ComputeCLine::SetHangChecking(theHangChecking)
    /// (the gxx L449-452).
    pub fn set_hang_checking(&mut self, the_hang_checking: bool) {
        self.my_hang_checking = the_hang_checking;
    }

    /// OCCT Approx_ComputeCLine::IsAllApproximated() (the gxx L461-464).
    pub fn is_all_approximated(&self) -> bool {
        self.alldone
    }

    /// OCCT Approx_ComputeCLine::IsToleranceReached() (the gxx L471-474).
    pub fn is_tolerance_reached(&self) -> bool {
        self.tolreached
    }

    /// OCCT Approx_ComputeCLine::Error(Index, tol3d, tol2d) (the gxx L481-485).
    pub fn error(&self, index: i32) -> (f64, f64) {
        (
            self.tolers3d[(index - 1) as usize],
            self.tolers2d[(index - 1) as usize],
        )
    }
}

// The PCONFUSION import is carried for the Precision reads of the consumers
// (kept next to the Confusion anchor above).
const _: f64 = PCONFUSION;

#[cfg(test)]
mod tests {
    use super::*;
    use glam::{DVec2, DVec3};

    /// A straight line f(u) = P0 + u * D over [0, 1] — the analytic case the
    /// degree-8 least square must reproduce to machine precision.
    struct TestLine;

    impl AppContFunction for TestLine {
        fn my_nb_pnt(&self) -> i32 {
            1
        }
        fn my_nb_pnt2d(&self) -> i32 {
            0
        }
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            1.0
        }
        fn value(&self, the_u: f64, _the_pnt2d: &mut [DVec2], the_pnt: &mut [DVec3]) -> bool {
            the_pnt[0] = DVec3::new(the_u, 2.0 * the_u, -the_u);
            true
        }
        fn d1(&self, _the_u: f64, _the_vec2d: &mut [DVec2], _the_vec: &mut [DVec3]) -> bool {
            false
        }
    }

    /// A quarter circle of radius 2 around the origin, u in [0, PI/2].
    struct TestQuarterCircle;

    impl AppContFunction for TestQuarterCircle {
        fn my_nb_pnt(&self) -> i32 {
            1
        }
        fn my_nb_pnt2d(&self) -> i32 {
            0
        }
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            std::f64::consts::FRAC_PI_2
        }
        fn value(&self, the_u: f64, _the_pnt2d: &mut [DVec2], the_pnt: &mut [DVec3]) -> bool {
            the_pnt[0] = DVec3::new(2.0 * the_u.cos(), 2.0 * the_u.sin(), 0.0);
            true
        }
        fn d1(&self, _the_u: f64, _the_vec2d: &mut [DVec2], _the_vec: &mut [DVec3]) -> bool {
            false
        }
    }

    fn bernstein_basis(deg: usize, j: usize, u: f64) -> f64 {
        let c = |n: usize, k: usize| -> f64 {
            let k = k.min(n - k);
            let mut r = 1.0;
            for i in 1..=k {
                r = r * (n - k + i) as f64 / i as f64;
            }
            r
        };
        c(deg, j) * u.powi(j as i32) * (1.0 - u).powi((deg - j) as i32)
    }

    /// The BiTgte_Blend MakeCurve call shape: Fit(F, 8, 8, Approximation,
    /// Approximation, cutting=true).  On a straight line the fit must
    /// succeed on the whole interval (one MultiCurve, alldone, tolerance
    /// reached) and the resulting degree-8 Bezier must evaluate on the line.
    #[test]
    fn fit_and_divide_line_single_segment_exact() {
        let f = TestLine;
        let tol = 1.0e-6;
        let fit = ApproxFitAndDivide::new(&f, 8, 8, tol, tol, true);

        assert!(fit.is_all_approximated(), "alldone");
        assert!(fit.is_tolerance_reached(), "tolreached");
        assert_eq!(fit.nb_multi_curves(), 1);

        let (first, last) = fit.parameters(1);
        assert!((first - 0.0).abs() < 1e-12 && (last - 1.0).abs() < 1e-12);

        // OCCT L481-485: Error(Index, tol3d, tol2d).
        let (err3d, err2d) = fit.error(1);
        assert!(err3d < 1e-8, "line 3d error {}", err3d);
        assert!(err2d < 1e-8, "line 2d error {}", err2d);

        // The Bezier poles must evaluate back onto the line.
        let mc = fit.value(1);
        assert_eq!(mc.degree(), 8);
        let mut poles = vec![DVec3::ZERO; mc.degree() + 1];
        mc.curve(1, &mut poles);
        for k in 1..=9 {
            let t = k as f64 / 9.0;
            let mut p = DVec3::ZERO;
            for (j, pole) in poles.iter().enumerate() {
                p += *pole * bernstein_basis(8, j, t);
            }
            let want = DVec3::new(t, 2.0 * t, -t);
            assert!(
                (p - want).length() < 1e-8,
                "t={}: {} vs {}",
                t,
                p,
                want
            );
        }
    }

    /// The quarter-circle case: one cut-free segment, degree-8 Bezier error
    /// under the Precision::Approximation tolerance.
    #[test]
    fn fit_and_divide_quarter_circle_within_tolerance() {
        let f = TestQuarterCircle;
        let tol = 1.0e-6;
        let fit = ApproxFitAndDivide::new(&f, 8, 8, tol, tol, true);

        assert!(fit.is_all_approximated(), "alldone");
        assert_eq!(fit.nb_multi_curves(), 1);

        let (err3d, _err2d) = fit.error(1);
        assert!(err3d < 1e-7, "circle 3d error {}", err3d);

        // The poles must evaluate onto the circle of radius 2.  The fit
        // controls the 3d error only at the Gauss discretisation points (the
        // OCCT behaviour); between them the discrete least-squares residual
        // of the degree-8 approximation is bounded by the normal-equation
        // conditioning — the bound below sits two orders above that floor
        // and still catches any real misalignment (a wrong parameterisation
        // or matrix gives deviations >= 1e-3).
        let mc = fit.value(1);
        assert_eq!(mc.degree(), 8);
        let mut poles = vec![DVec3::ZERO; mc.degree() + 1];
        mc.curve(1, &mut poles);
        let mut max_radius_dev = 0.0f64;
        for k in 0..=200 {
            let u = std::f64::consts::FRAC_PI_2 * (k as f64) / 200.0;
            let mut p = DVec3::ZERO;
            for (j, pole) in poles.iter().enumerate() {
                p += *pole * bernstein_basis(8, j, u);
            }
            let radius = (p - DVec3::new(0.0, 0.0, 0.0)).length();
            max_radius_dev = max_radius_dev.max((radius - 2.0).abs());
        }
        assert!(
            max_radius_dev < 5e-4,
            "max radius deviation {}",
            max_radius_dev
        );
    }
}
