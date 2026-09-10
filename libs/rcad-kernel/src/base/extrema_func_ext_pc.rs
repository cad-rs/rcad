//! OCCT Extrema_GFuncExtPC (TKGeomBase/Extrema/Extrema_GFuncExtPC.hxx L32-480)
//! — the point-curve extremum distance function solved by math_FunctionRoots.
//!
//! OCCT template instantiation (Extrema_PCFOfEPCOfExtPC.hxx L28-34):
//! `Extrema_PCFOfEPCOfExtPC = Extrema_GFuncExtPC<Adaptor3d_Curve,
//! Extrema_CurveTool, Extrema_POnCurv, gp_Pnt, gp_Vec,
//! NCollection_Sequence<Extrema_POnCurv>>`, so this single translation unit is
//! both the `Extrema_GFuncExtPC` body and the `Extrema_PCFOfEPCOfExtPC` alias.
//!
//! If D1c and D2c are the first and second derivatives:
//! F(u) = (C(u)-P).D1c(u) / ||D1c||
//! DF(u) = ||D1c|| + (C(u)-P).D2c(u)/||D1c|| - F(u)*D2c.D1c/||D1c||^2
//!
//! Template-parameter mapping: `TheCurve` + `TheCurveTool` map to the
//! [`ExtPCurveTool`] trait (the rcad `ExtremaCurveTool` facade extended with
//! the Extrema_CurveTool.hxx L128-141 Bezier/BSpline statics, defined in
//! `extrema_ext_pc`); `ThePOnC` / `TheSeqPOnC` map to [`POnCurve`]; `ThePoint`
//! and `TheVector` map to `glam::DVec3` (the rcad gp_Pnt / gp_Vec).
//!
//! OCCT `math_FunctionWithDerivative` (the Value / Derivative / Values /
//! GetStateNumber virtuals) is represented by the rcad
//! `math::root::FunctionWithDerivative` trait; the OCCT bool return + out
//! parameter signature maps to `Option` (None = the OCCT `false` return).

use glam::DVec3;

use crate::base::extrema::POnCurve;
use crate::base::extrema_ext_elc::{REAL_FIRST, REAL_LAST};
use crate::base::extrema_ext_pc::ExtPCurveTool;
use crate::base::proj_lib::CurveType;
use crate::core::precision::{is_infinite_value, INFINITE_VALUE};
use crate::math::root::{FunctionValue, FunctionWithDerivative};

/// hxx L460: static constexpr double TolFactor = 1.e-12.
const TOL_FACTOR: f64 = 1.0e-12;
/// hxx L461: static constexpr double MinTol = 1.e-20.
const MIN_TOL: f64 = 1.0e-20;
/// hxx L462: static constexpr double MinStep = 1.e-7.
const MIN_STEP: f64 = 1.0e-7;
/// hxx L463: static constexpr int MaxOrder = 3.
const MAX_ORDER: i32 = 3;

/// OCCT Extrema_GFuncExtPC (hxx L47-480); the `Extrema_PCFOfEPCOfExtPC` alias.
pub struct FuncExtPC<'a> {
    /// hxx L465: ThePoint myP.
    my_p: DVec3,
    /// hxx L466: TheCurve* myC.
    my_c: Option<&'a dyn ExtPCurveTool>,
    /// hxx L467: double myU.
    my_u: f64,
    /// hxx L468: ThePoint myPc.
    my_pc: DVec3,
    /// hxx L469: double myD1f.
    my_d1f: f64,
    /// hxx L470: NCollection_Sequence<double> mySqDist.
    my_sq_dist: Vec<f64>,
    /// hxx L471: NCollection_Sequence<int> myIsMin.
    my_is_min: Vec<i32>,
    /// hxx L472: TheSeqPOnC myPoint.
    my_point: Vec<POnCurve>,
    /// hxx L473: bool myPinit.
    my_pinit: bool,
    /// hxx L474: bool myCinit.
    my_cinit: bool,
    /// hxx L475: bool myD1Init — written by the OCCT body, never read there.
    #[allow(dead_code)]
    my_d1init: bool,
    /// hxx L476: double myTol.
    my_tol: f64,
    /// hxx L477: int myMaxDerivOrder.
    my_max_deriv_order: i32,
    /// hxx L478: double myUinfium.
    my_u_inf_ium: f64,
    /// hxx L479: double myUsupremum.
    my_u_sup_remum: f64,
}

impl<'a> FuncExtPC<'a> {
    /// OCCT Extrema_GFuncExtPC() default constructor (hxx L59-71).
    pub fn new() -> Self {
        FuncExtPC {
            my_p: DVec3::ZERO,
            my_c: None,
            my_u: 0.0,
            my_pc: DVec3::ZERO,
            my_d1f: 0.0,
            my_sq_dist: Vec::new(),
            my_is_min: Vec::new(),
            my_point: Vec::new(),
            my_pinit: false,
            my_cinit: false,
            my_d1init: false,
            my_tol: MIN_TOL,
            my_max_deriv_order: 0,
            my_u_inf_ium: 0.0,
            my_u_sup_remum: 0.0,
        }
    }

    /// OCCT Extrema_GFuncExtPC(const ThePoint& theP, const TheCurve& theC)
    /// (hxx L76-101).
    pub fn new_point_curve(the_p: DVec3, the_c: &'a dyn ExtPCurveTool) -> Self {
        let mut this = FuncExtPC {
            my_p: the_p,
            my_c: Some(the_c),
            my_u: 0.0,
            my_pc: DVec3::ZERO,
            my_d1f: 0.0,
            my_sq_dist: Vec::new(),
            my_is_min: Vec::new(),
            my_point: Vec::new(),
            my_pinit: true,
            my_cinit: true,
            my_d1init: false,
            my_tol: MIN_TOL,
            my_max_deriv_order: 0,
            my_u_inf_ium: 0.0,
            my_u_sup_remum: 0.0,
        };
        // hxx L85.
        this.sub_interval_initialize(the_c.first_parameter(), the_c.last_parameter());

        // hxx L87-100.  OCCT lists GeomAbs_BezierCurve, GeomAbs_BSplineCurve,
        // GeomAbs_OffsetCurve, GeomAbs_OtherCurve; the rcad CurveType enum
        // (proj_lib/mod.rs L44-53) carries no offset value — an OCCT offset
        // curve reports Other at the rcad adaptor layer.
        match the_c.get_type() {
            CurveType::Bezier | CurveType::BSpline | CurveType::Other => {
                this.my_max_deriv_order = MAX_ORDER;
                this.my_tol = this.search_of_tolerance();
            }
            _ => {
                this.my_max_deriv_order = 0;
                this.my_tol = MIN_TOL;
            }
        }
        this
    }

    /// OCCT Initialize(const TheCurve& theC) (hxx L105-129).
    pub fn initialize(&mut self, the_c: &'a dyn ExtPCurveTool) {
        self.my_c = Some(the_c);
        self.my_cinit = true;
        self.my_point.clear();
        self.my_sq_dist.clear();
        self.my_is_min.clear();

        // hxx L113.
        self.sub_interval_initialize(the_c.first_parameter(), the_c.last_parameter());

        // hxx L115-128 (same type switch as the constructor).
        match the_c.get_type() {
            CurveType::Bezier | CurveType::BSpline | CurveType::Other => {
                self.my_max_deriv_order = MAX_ORDER;
                self.my_tol = self.search_of_tolerance();
            }
            _ => {
                self.my_max_deriv_order = 0;
                self.my_tol = MIN_TOL;
            }
        }
    }

    /// OCCT SetPoint(const ThePoint& theP) (hxx L133-140).
    pub fn set_point(&mut self, the_p: DVec3) {
        self.my_p = the_p;
        self.my_pinit = true;
        self.my_point.clear();
        self.my_sq_dist.clear();
        self.my_is_min.clear();
    }

    /// OCCT NbExt() (hxx L382).
    pub fn nb_ext(&self) -> usize {
        self.my_sq_dist.len()
    }

    /// OCCT SquareDistance(const int theN) (hxx L386-393) — 1-based.
    pub fn square_distance(&self, the_n: usize) -> f64 {
        if !self.my_pinit || !self.my_cinit {
            panic!("Standard_TypeMismatch");
        }
        self.my_sq_dist[the_n - 1]
    }

    /// OCCT IsMin(const int theN) (hxx L397-404) — 1-based.
    pub fn is_min(&self, the_n: usize) -> bool {
        if !self.my_pinit || !self.my_cinit {
            panic!("Standard_TypeMismatch");
        }
        self.my_is_min[the_n - 1] == 1
    }

    /// OCCT Point(const int theN) (hxx L408-415) — 1-based.
    pub fn point(&self, the_n: usize) -> &POnCurve {
        if !self.my_pinit || !self.my_cinit {
            panic!("Standard_TypeMismatch");
        }
        &self.my_point[the_n - 1]
    }

    /// OCCT SubIntervalInitialize(theUfirst, theUlast) (hxx L420-424).
    pub fn sub_interval_initialize(&mut self, the_ufirst: f64, the_ulast: f64) {
        self.my_u_inf_ium = the_ufirst;
        self.my_u_sup_remum = the_ulast;
    }

    /// OCCT SearchOfTolerance() (hxx L428-457) — if 1st derivative of curve
    /// |D1| < Tol, it is considered D1=0.
    pub fn search_of_tolerance(&mut self) -> f64 {
        let c = self.my_c.expect("Extrema_GFuncExtPC: no curve");
        let n_point = 10i32;
        let a_step = (self.my_u_sup_remum - self.my_u_inf_ium) / n_point as f64;

        let mut a_num = 0i32;
        let mut a_max = -INFINITE_VALUE;

        // OCCT do { ... } while (++aNum < NPoint + 1); the `continue` inside
        // the OCCT body reaches the increment, so the increment runs on every
        // path before the condition is re-tested.
        loop {
            // hxx L438-440.
            let mut u = self.my_u_inf_ium + a_num as f64 * a_step;
            if u > self.my_u_sup_remum {
                u = self.my_u_sup_remum;
            }

            // hxx L442-444.
            let (_ptemp, v_der) = c.d1(u);

            // hxx L446-449: the infinite check skips the sample.
            if !(is_infinite_value(v_der.x) || is_infinite_value(v_der.y)) {
                let vm = v_der.length();
                if vm > a_max {
                    a_max = vm;
                }
            }

            a_num += 1;
            if a_num < n_point + 1 {
                continue;
            }
            break;
        }

        // hxx L456.
        (a_max * TOL_FACTOR).max(MIN_TOL)
    }
}

impl FunctionValue for FuncExtPC<'_> {
    /// OCCT Value(const double theU, double& theF) (hxx L146-253).
    fn value(&mut self, the_u: f64) -> Option<f64> {
        if !self.my_pinit || !self.my_cinit {
            panic!("Standard_TypeMismatch: No init");
        }

        self.my_u = the_u;
        let c = self.my_c.expect("Extrema_GFuncExtPC: no curve");
        let (pc, mut d1c) = c.d1(self.my_u);
        self.my_pc = pc;

        let mut the_f = 0.0;
        // hxx L157: OCCT tests only the X and Y components.
        if is_infinite_value(d1c.x) || is_infinite_value(d1c.y) {
            the_f = INFINITE_VALUE;
            return None;
        }

        let mut ndu = d1c.length();

        if self.my_max_deriv_order != 0 {
            if ndu <= self.my_tol {
                // Singular case (hxx L167-241).
                let division_factor = 1.0e-3;
                let du = if (self.my_u_sup_remum >= REAL_LAST)
                    || (self.my_u_inf_ium <= REAL_FIRST)
                {
                    0.0
                } else {
                    self.my_u_sup_remum - self.my_u_inf_ium
                };

                let a_delta = (du * division_factor).max(MIN_STEP);
                // Derivative is approximated by Taylor-series.
                // hxx L179-188: do { V = DN(*myC, myU, ++n); ... } while (...).
                let mut n = 1i32; // Derivative order
                let mut v = DVec3::ZERO; // OCCT TheVector V; (default gp_Vec)
                let mut is_derive_found;
                loop {
                    n += 1;
                    v = c.dn(self.my_u, n);
                    ndu = v.length();
                    is_derive_found = ndu > self.my_tol;
                    if !(!is_derive_found && n < self.my_max_deriv_order) {
                        break;
                    }
                }

                if is_derive_found {
                    // hxx L192-197.
                    let u = if self.my_u - self.my_u_inf_ium < a_delta {
                        self.my_u + a_delta
                    } else {
                        self.my_u - a_delta
                    };

                    // hxx L199-201.
                    let p1 = c.d0(self.my_u.min(u));
                    let p2 = c.d0(self.my_u.max(u));

                    // hxx L203-204: TheVector V1(P1, P2) = P2 - P1.
                    let v1 = p2 - p1;
                    let a_dir_factor = v.dot(v1);

                    if a_dir_factor < 0.0 {
                        d1c = -v;
                    } else {
                        d1c = v;
                    }
                } else {
                    // Derivative is approximated by three points
                    // (hxx L211-239).  OCCT measures from the default
                    // (zero) point Ptemp.
                    let ptemp = DVec3::ZERO;
                    let p1;
                    let p2;
                    let p3;
                    let is_parameter_grown;

                    if self.my_u - self.my_u_inf_ium < 2.0 * a_delta {
                        p1 = c.d0(self.my_u);
                        p2 = c.d0(self.my_u + a_delta);
                        p3 = c.d0(self.my_u + 2.0 * a_delta);
                        is_parameter_grown = true;
                    } else {
                        p1 = c.d0(self.my_u - 2.0 * a_delta);
                        p2 = c.d0(self.my_u - a_delta);
                        p3 = c.d0(self.my_u);
                        is_parameter_grown = false;
                    }

                    // hxx L233: TheVector V1(Ptemp, P1), V2(Ptemp, P2),
                    // V3(Ptemp, P3).
                    let v1 = p1 - ptemp;
                    let v2 = p2 - ptemp;
                    let v3 = p3 - ptemp;

                    if is_parameter_grown {
                        d1c = -3.0 * v1 + 4.0 * v2 - v3;
                    } else {
                        d1c = v1 - 4.0 * v2 + 3.0 * v3;
                    }
                }
                ndu = d1c.length();
            }
        }

        if ndu <= MIN_TOL {
            // Warning: 1st derivative is equal to zero!
            return None;
        }

        // hxx L250-251: TheVector PPc(myP, myPc) = myPc - myP.
        let ppc = self.my_pc - self.my_p;
        the_f = ppc.dot(d1c) / ndu;
        Some(the_f)
    }
}

impl FunctionWithDerivative for FuncExtPC<'_> {
    /// OCCT Derivative(const double theU, double& theDF) (hxx L259-267).
    fn derivative(&mut self, the_u: f64) -> Option<f64> {
        if !self.my_pinit || !self.my_cinit {
            panic!("Standard_TypeMismatch");
        }
        // hxx L265-266: double F; return Values(theU, F, theDF);
        match self.values(the_u) {
            Some((_f, the_df)) => Some(the_df),
            None => None,
        }
    }

    /// OCCT Values(const double theU, double& theF, double& theDF)
    /// (hxx L274-353).
    fn values(&mut self, the_u: f64) -> Option<(f64, f64)> {
        if !self.my_pinit || !self.my_cinit {
            panic!("Standard_TypeMismatch: No init");
        }

        // hxx L281.
        let my_pc_old = self.my_pc;
        let my_p_old = self.my_p;

        let mut the_f = 0.0;
        // hxx L283-287.
        match self.value(the_u) {
            Some(f) => the_f = f,
            None => {
                self.my_d1init = false;
                return None;
            }
        }

        // hxx L289-291.
        self.my_u = the_u;
        self.my_pc = my_pc_old;
        self.my_p = my_p_old;

        // hxx L293-294.
        let c = self.my_c.expect("Extrema_GFuncExtPC: no curve");
        let (pc, d1c, d2c) = c.d2(self.my_u);
        self.my_pc = pc;

        let ndu = d1c.length();
        let the_df;
        if ndu <= self.my_tol {
            // Singular case (hxx L297-342): derivative is approximated by
            // three points.
            let division_factor = 0.01;
            let du = if (self.my_u_sup_remum >= REAL_LAST)
                || (self.my_u_inf_ium <= REAL_FIRST)
            {
                0.0
            } else {
                self.my_u_sup_remum - self.my_u_inf_ium
            };

            let a_delta = (du * division_factor).max(MIN_STEP);

            let f1: f64;
            let f2: f64;
            let f3: f64;

            if self.my_u - self.my_u_inf_ium < 2.0 * a_delta {
                f1 = the_f;
                let u2 = self.my_u + a_delta;
                let u3 = self.my_u + a_delta * 2.0;

                // hxx L317: if (!((Value(U2, F2)) && (Value(U3, F3)))).
                let f2v = match self.value(u2) {
                    Some(v) => v,
                    None => {
                        self.my_d1init = false;
                        return None;
                    }
                };
                let f3v = match self.value(u3) {
                    Some(v) => v,
                    None => {
                        self.my_d1init = false;
                        return None;
                    }
                };
                f2 = f2v;
                f3 = f3v;

                the_df = (-3.0 * f1 + 4.0 * f2 - f3) / (2.0 * a_delta);
            } else {
                f3 = the_f;
                let u1 = self.my_u - a_delta * 2.0;
                let u2 = self.my_u - a_delta;

                // hxx L331: if (!((Value(U2, F2)) && (Value(U1, F1)))).
                let f2v = match self.value(u2) {
                    Some(v) => v,
                    None => {
                        self.my_d1init = false;
                        return None;
                    }
                };
                let f1v = match self.value(u1) {
                    Some(v) => v,
                    None => {
                        self.my_d1init = false;
                        return None;
                    }
                };
                f2 = f2v;
                f1 = f1v;

                the_df = (f1 - 4.0 * f2 + 3.0 * f3) / (2.0 * a_delta);
            }
            // hxx L339-341.
            self.my_u = the_u;
            self.my_pc = my_pc_old;
            self.my_p = my_p_old;
        } else {
            // hxx L345-346: TheVector PPc(myP, myPc).
            let ppc = self.my_pc - self.my_p;
            the_df = ndu + (ppc.dot(d2c) / ndu) - the_f * (d1c.dot(d2c)) / (ndu * ndu);
        }

        self.my_d1f = the_df;

        self.my_d1init = true;
        Some((the_f, the_df))
    }

    /// OCCT GetStateNumber() (hxx L357-379) — saves the found extremum.
    fn get_state_number(&mut self) -> i32 {
        if !self.my_pinit || !self.my_cinit {
            panic!("Standard_TypeMismatch");
        }
        self.my_sq_dist.push((self.my_pc - self.my_p).length_squared());

        // It is necessary to always compute myD1f.
        self.my_d1init = true;
        let _ = self.values(self.my_u);

        let mut int_val = 0;
        if self.my_d1f > 0.0 {
            int_val = 1;
        }

        self.my_is_min.push(int_val);
        self.my_point.push(POnCurve {
            param: self.my_u,
            point: self.my_pc,
        });
        0
    }
}
