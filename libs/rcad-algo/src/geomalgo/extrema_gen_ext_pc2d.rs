// OCCT Extrema_EPCOfExtPC2d / Extrema_PCFOfEPCOfExtPC2d (TKGeomBase) —
// the point-to-2D-curve extremum engine.
//
// OCCT 8.0 defines both as template aliases over the generic engines:
//   Extrema_PCFOfEPCOfExtPC2d = Extrema_GFuncExtPC<Adaptor2d_Curve2d,
//       Extrema_Curve2dTool, Extrema_POnCurv2d, gp_Pnt2d, gp_Vec2d,
//       NCollection_Sequence<Extrema_POnCurv2d>>           (Extrema_GFuncExtPC.hxx)
//   Extrema_EPCOfExtPC2d      = Extrema_GGenExtPC<Adaptor2d_Curve2d,
//       Extrema_Curve2dTool, Extrema_POnCurv2d, gp_Pnt2d,
//       Extrema_PCFOfEPCOfExtPC2d>                         (Extrema_GGenExtPC.hxx)
//
// rcad mapping (along the established template-argument-to-trait pattern):
//   TheCurve     = dyn Curve2dAdaptor   (Adaptor2d_Curve2d)
//   TheCurveTool = module fns over Curve2dAdaptor (Extrema_Curve2dTool —
//                  thin statics delegating to the adaptor, mirroring the
//                  OCCT hxx one-liners)
//   ThePOnC      = POnCurv2d (param + point)
//   ThePoint     = DVec2 (gp_Pnt2d), TheVector = DVec2 (gp_Vec2d)
//
// Consumer: Contap_HContTool::Project (TKHLR).
//
// Note: geom2d_int.rs carries a private GFuncExtPC copy used by the
// GenLocateExtPC projection (Extrema_GenLocateExtPC); this module is the
// full Extrema_GGenExtPC engine driven by math_FunctionRoots.

use glam::DVec2;
use rcad_kernel::math::root::{FunctionValue, FunctionWithDerivative, FunctionRoots};

use super::geom2d_int::{Curve2dAdaptor, Curve2dType, geom2d_curve_tool};

/// OCCT Extrema_POnCurv2d — a point on a 2D curve with its parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct POnCurv2d {
    param: f64,
    pt: DVec2,
}

impl POnCurv2d {
    /// OCCT Extrema_POnCurv2d(U, P).
    pub fn new(u: f64, p: DVec2) -> Self {
        POnCurv2d { param: u, pt: p }
    }
    /// OCCT Parameter().
    pub fn parameter(&self) -> f64 {
        self.param
    }
    /// OCCT Value().
    pub fn value(&self) -> DVec2 {
        self.pt
    }
}

// OCCT Extrema_Curve2dTool statics used by the engines — thin delegates to
// the adaptor (Extrema_Curve2dTool.hxx L43-128).
mod curve2d_tool {
    use super::*;

    pub fn first_parameter(c: &dyn Curve2dAdaptor) -> f64 {
        c.first_parameter()
    }
    pub fn last_parameter(c: &dyn Curve2dAdaptor) -> f64 {
        c.last_parameter()
    }
    pub fn d0(c: &dyn Curve2dAdaptor, u: f64) -> DVec2 {
        c.value(u)
    }
    pub fn d1(c: &dyn Curve2dAdaptor, u: f64) -> (DVec2, DVec2) {
        c.d1(u)
    }
    pub fn d2(c: &dyn Curve2dAdaptor, u: f64) -> (DVec2, DVec2, DVec2) {
        c.d2(u)
    }
    pub fn dn(c: &dyn Curve2dAdaptor, u: f64, n: i32) -> DVec2 {
        c.dn(u, n)
    }
    pub fn get_type(c: &dyn Curve2dAdaptor) -> Curve2dType {
        c.get_type()
    }
}

// OCCT Extrema_GFuncExtPC private constants (hxx L460-463).
const TOL_FACTOR: f64 = 1.0e-12;
const MIN_TOL: f64 = 1.0e-20;
const MIN_STEP: f64 = 1.0e-7;
const MAX_ORDER: i32 = 3;

/// OCCT Extrema_GFuncExtPC (Extrema_GFuncExtPC.hxx L47-480) — the 2D
/// instantiation `Extrema_PCFOfEPCOfExtPC2d`: F(u) = (C(u)-P)·D1c/|D1c|
/// with derivative, consumed by math_FunctionRoots.
pub struct PCFOfEPCOfExtPC2d<'a> {
    p: DVec2,
    c: Option<&'a dyn Curve2dAdaptor>,
    u: f64,
    pc: DVec2,
    d1f: f64,
    sq_dist: Vec<f64>,
    is_min: Vec<i32>,
    point: Vec<POnCurv2d>,
    pinit: bool,
    cinit: bool,
    d1_init: bool,
    tol: f64,
    max_deriv_order: i32,
    uinfium: f64,
    usupremum: f64,
}

impl<'a> PCFOfEPCOfExtPC2d<'a> {
    /// OCCT Extrema_GFuncExtPC() default constructor (L59-71).
    pub fn new() -> Self {
        PCFOfEPCOfExtPC2d {
            p: DVec2::ZERO,
            c: None,
            u: 0.0,
            pc: DVec2::ZERO,
            d1f: 0.0,
            sq_dist: Vec::new(),
            is_min: Vec::new(),
            point: Vec::new(),
            pinit: false,
            cinit: false,
            d1_init: false,
            tol: MIN_TOL,
            max_deriv_order: 0,
            uinfium: 0.0,
            usupremum: 0.0,
        }
    }

    /// OCCT Extrema_GFuncExtPC(theP, theC) constructor (L76-101).
    pub fn with_point_curve(the_p: DVec2, the_c: &'a dyn Curve2dAdaptor) -> Self {
        let mut f = PCFOfEPCOfExtPC2d {
            p: the_p,
            c: Some(the_c),
            u: 0.0,
            pc: DVec2::ZERO,
            d1f: 0.0,
            sq_dist: Vec::new(),
            is_min: Vec::new(),
            point: Vec::new(),
            pinit: true,
            cinit: true,
            d1_init: false,
            tol: MIN_TOL,
            max_deriv_order: 0,
            uinfium: 0.0,
            usupremum: 0.0,
        };
        f.sub_interval_initialize(
            curve2d_tool::first_parameter(the_c),
            curve2d_tool::last_parameter(the_c),
        );
        match curve2d_tool::get_type(the_c) {
            Curve2dType::BezierCurve
            | Curve2dType::BSplineCurve
            | Curve2dType::OffsetCurve
            | Curve2dType::OtherCurve => {
                f.max_deriv_order = MAX_ORDER;
                f.tol = f.search_of_tolerance();
            }
            _ => {
                f.max_deriv_order = 0;
                f.tol = MIN_TOL;
            }
        }
        f
    }

    /// OCCT Initialize(theC) (L105-129).
    pub fn initialize(&mut self, the_c: &'a dyn Curve2dAdaptor) {
        self.c = Some(the_c);
        self.cinit = true;
        self.point.clear();
        self.sq_dist.clear();
        self.is_min.clear();

        self.sub_interval_initialize(
            curve2d_tool::first_parameter(the_c),
            curve2d_tool::last_parameter(the_c),
        );

        match curve2d_tool::get_type(the_c) {
            Curve2dType::BezierCurve
            | Curve2dType::BSplineCurve
            | Curve2dType::OffsetCurve
            | Curve2dType::OtherCurve => {
                self.max_deriv_order = MAX_ORDER;
                self.tol = self.search_of_tolerance();
            }
            _ => {
                self.max_deriv_order = 0;
                self.tol = MIN_TOL;
            }
        }
    }

    /// OCCT SetPoint(theP) (L133-140).
    pub fn set_point(&mut self, the_p: DVec2) {
        self.p = the_p;
        self.pinit = true;
        self.point.clear();
        self.sq_dist.clear();
        self.is_min.clear();
    }

    /// OCCT Value(theU, theF) (L146-253).
    fn value_f(&mut self, the_u: f64, the_f: &mut f64) -> bool {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch: No init");
        }

        self.u = the_u;
        let (pc, mut d1c) = curve2d_tool::d1(self.c.unwrap(), self.u);
        self.pc = pc;

        if d1c.x.is_infinite() || d1c.y.is_infinite() {
            *the_f = f64::INFINITY;
            return false;
        }

        let mut ndu = d1c.length();

        if self.max_deriv_order != 0 {
            if ndu <= self.tol {
                // Singular case
                let division_factor = 1.0e-3;
                let du = if self.usupremum >= f64::MAX || self.uinfium <= f64::MIN {
                    0.0
                } else {
                    self.usupremum - self.uinfium
                };

                let a_delta = (du * division_factor).max(MIN_STEP);
                // Derivative is approximated by Taylor-series

                let mut n = 1; // Derivative order
                let mut v = DVec2::ZERO;
                let mut is_deriv_found;

                loop {
                    n += 1;
                    v = curve2d_tool::dn(self.c.unwrap(), self.u, n);
                    ndu = v.length();
                    is_deriv_found = ndu > self.tol;
                    if is_deriv_found || n >= self.max_deriv_order {
                        break;
                    }
                }

                if is_deriv_found {
                    let u2;
                    if self.u - self.uinfium < a_delta {
                        u2 = self.u + a_delta;
                    } else {
                        u2 = self.u - a_delta;
                    }

                    let p1 = curve2d_tool::d0(self.c.unwrap(), self.u.min(u2));
                    let p2 = curve2d_tool::d0(self.c.unwrap(), self.u.max(u2));

                    let v1 = p2 - p1;
                    let a_dir_factor = v.dot(v1);

                    if a_dir_factor < 0.0 {
                        d1c = -v;
                    } else {
                        d1c = v;
                    }
                } else {
                    // Derivative is approximated by three points
                    let is_parameter_grown;
                    let (p1, p2, p3);
                    if self.u - self.uinfium < 2.0 * a_delta {
                        p1 = curve2d_tool::d0(self.c.unwrap(), self.u);
                        p2 = curve2d_tool::d0(self.c.unwrap(), self.u + a_delta);
                        p3 = curve2d_tool::d0(self.c.unwrap(), self.u + 2.0 * a_delta);
                        is_parameter_grown = true;
                    } else {
                        p1 = curve2d_tool::d0(self.c.unwrap(), self.u - 2.0 * a_delta);
                        p2 = curve2d_tool::d0(self.c.unwrap(), self.u - a_delta);
                        p3 = curve2d_tool::d0(self.c.unwrap(), self.u);
                        is_parameter_grown = false;
                    }

                    let v1 = p1 - self.pc;
                    let v2 = p2 - self.pc;
                    let v3 = p3 - self.pc;

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
            return false;
        }

        let ppc = self.pc - self.p; // OCCT TheVector PPc(myP, myPc) = myPc - myP
        *the_f = ppc.dot(d1c) / ndu;
        true
    }

    /// OCCT Derivative(theU, theDF) (L259-267).
    fn derivative_f(&mut self, the_u: f64, the_df: &mut f64) -> bool {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        let mut f = 0.0;
        self.values_f(the_u, &mut f, the_df)
    }

    /// OCCT Values(theU, theF, theDF) (L274-353).
    fn values_f(&mut self, the_u: f64, the_f: &mut f64, the_df: &mut f64) -> bool {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch: No init");
        }

        let my_pc_old = self.pc;
        let my_p_old = self.p;

        if !self.value_f(the_u, the_f) {
            self.d1_init = false;
            return false;
        }

        self.u = the_u;
        self.pc = my_pc_old;
        self.p = my_p_old;

        let (pc, d1c, d2c) = curve2d_tool::d2(self.c.unwrap(), self.u);
        self.pc = pc;

        let ndu = d1c.length();
        if ndu <= self.tol {
            // Singular case — derivative approximated by three points
            let division_factor = 0.01;
            let du = if self.usupremum >= f64::MAX || self.uinfium <= f64::MIN {
                0.0
            } else {
                self.usupremum - self.uinfium
            };

            let a_delta = (du * division_factor).max(MIN_STEP);

            let the_df_v;
            if self.u - self.uinfium < 2.0 * a_delta {
                let f1 = *the_f;
                let u2 = self.u + a_delta;
                let u3 = self.u + a_delta * 2.0;

                let mut f2 = 0.0;
                let mut f3 = 0.0;
                if !(self.value_f(u2, &mut f2) && self.value_f(u3, &mut f3)) {
                    self.d1_init = false;
                    return false;
                }

                the_df_v = (-3.0 * f1 + 4.0 * f2 - f3) / (2.0 * a_delta);
            } else {
                let f3 = *the_f;
                let u1 = self.u - a_delta * 2.0;
                let u2 = self.u - a_delta;

                let mut f1 = 0.0;
                let mut f2 = 0.0;
                if !(self.value_f(u2, &mut f2) && self.value_f(u1, &mut f1)) {
                    self.d1_init = false;
                    return false;
                }

                the_df_v = (f1 - 4.0 * f2 + 3.0 * f3) / (2.0 * a_delta);
            }
            *the_df = the_df_v;
            self.u = the_u;
            self.pc = my_pc_old;
            self.p = my_p_old;
        } else {
            let ppc = self.pc - self.p; // OCCT TheVector PPc(myP, myPc) = myPc - myP
            *the_df = ndu + (ppc.dot(d2c) / ndu) - *the_f * (d1c.dot(d2c)) / (ndu * ndu);
        }

        self.d1f = *the_df;

        self.d1_init = true;
        true
    }

    /// OCCT GetStateNumber (L357-379).
    fn get_state_number_f(&mut self) -> i32 {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        self.sq_dist.push(self.pc.distance_squared(self.p));

        // It is necessary to always compute myD1f.
        self.d1_init = true;
        let mut ff = 0.0;
        let mut dd = 0.0;
        self.values_f(self.u, &mut ff, &mut dd);

        let int_val = if self.d1f > 0.0 { 1 } else { 0 };

        self.is_min.push(int_val);
        self.point.push(POnCurv2d::new(self.u, self.pc));
        0
    }

    /// OCCT NbExt (L382).
    pub fn nb_ext(&self) -> usize {
        self.sq_dist.len()
    }

    /// OCCT SquareDistance(theN) (L386-393).
    pub fn square_distance(&self, the_n: usize) -> f64 {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        self.sq_dist[the_n - 1]
    }

    /// OCCT IsMin(theN) (L397-404).
    pub fn is_min(&self, the_n: usize) -> bool {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        self.is_min[the_n - 1] == 1
    }

    /// OCCT Point(theN) (L408-415).
    pub fn point(&self, the_n: usize) -> &POnCurv2d {
        if !self.pinit || !self.cinit {
            panic!("Standard_TypeMismatch");
        }
        &self.point[the_n - 1]
    }

    /// OCCT SubIntervalInitialize(theUfirst, theUlast) (L420-424).
    pub fn sub_interval_initialize(&mut self, the_ufirst: f64, the_ulast: f64) {
        self.uinfium = the_ufirst;
        self.usupremum = the_ulast;
    }

    /// OCCT SearchOfTolerance (L428-457).
    fn search_of_tolerance(&mut self) -> f64 {
        let n_point = 10;
        let a_step = (self.usupremum - self.uinfium) / n_point as f64;

        let mut a_num = 0;
        let mut a_max = f64::NEG_INFINITY;

        loop {
            let mut u = self.uinfium + a_num as f64 * a_step;
            if u > self.usupremum {
                u = self.usupremum;
            }

            let (_, v_der) = curve2d_tool::d1(self.c.unwrap(), u);

            if !(v_der.x.is_infinite() || v_der.y.is_infinite()) {
                let vm = v_der.length();
                if vm > a_max {
                    a_max = vm;
                }
            }

            a_num += 1;
            if a_num >= n_point + 1 {
                break;
            }
        }

        (a_max * TOL_FACTOR).max(MIN_TOL)
    }
}

impl Default for PCFOfEPCOfExtPC2d<'_> {
    fn default() -> Self {
        Self::new()
    }
}

impl FunctionValue for PCFOfEPCOfExtPC2d<'_> {
    fn value(&mut self, x: f64) -> Option<f64> {
        let mut f = 0.0;
        if self.value_f(x, &mut f) {
            Some(f)
        } else {
            None
        }
    }
}

impl FunctionWithDerivative for PCFOfEPCOfExtPC2d<'_> {
    fn derivative(&mut self, x: f64) -> Option<f64> {
        let mut d = 0.0;
        if self.derivative_f(x, &mut d) {
            Some(d)
        } else {
            None
        }
    }
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let mut f = 0.0;
        let mut d = 0.0;
        if self.values_f(x, &mut f, &mut d) {
            Some((f, d))
        } else {
            None
        }
    }
    fn get_state_number(&mut self) -> i32 {
        self.get_state_number_f()
    }
}

/// OCCT Extrema_GGenExtPC (Extrema_GGenExtPC.hxx L37-234) — the 2D
/// instantiation `Extrema_EPCOfExtPC2d`: all extrem distances between a
/// point and a curve, driven by math_FunctionRoots.
pub struct EPCOfExtPC2d<'a> {
    done: bool,
    init: bool,
    nbsample: i32,
    umin: f64,
    usup: f64,
    tolu: f64,
    tolf: f64,
    f: PCFOfEPCOfExtPC2d<'a>,
}

impl<'a> EPCOfExtPC2d<'a> {
    /// OCCT Extrema_GGenExtPC() default constructor (L44-53).
    pub fn new() -> Self {
        EPCOfExtPC2d {
            done: false,
            init: false,
            nbsample: 0,
            umin: 0.0,
            usup: 0.0,
            tolu: 0.0,
            tolf: 0.0,
            f: PCFOfEPCOfExtPC2d::new(),
        }
    }

    /// OCCT Extrema_GGenExtPC(theP, theC, theNbSample, theTolU, theTolF)
    /// (L61-70).
    pub fn with_full_domain(
        the_p: DVec2,
        the_c: &'a dyn Curve2dAdaptor,
        the_nb_sample: i32,
        the_tol_u: f64,
        the_tol_f: f64,
    ) -> Self {
        let mut e = EPCOfExtPC2d {
            done: false,
            init: false,
            nbsample: 0,
            umin: 0.0,
            usup: 0.0,
            tolu: 0.0,
            tolf: 0.0,
            f: PCFOfEPCOfExtPC2d::with_point_curve(the_p, the_c),
        };
        e.initialize_c(the_c, the_nb_sample, the_tol_u, the_tol_f);
        e.perform(the_p);
        e
    }

    /// OCCT Extrema_GGenExtPC(theP, theC, theNbSample, theUmin, theUsup,
    /// theTolU, theTolF) (L80-91).
    pub fn with_range(
        the_p: DVec2,
        the_c: &'a dyn Curve2dAdaptor,
        the_nb_sample: i32,
        the_umin: f64,
        the_usup: f64,
        the_tol_u: f64,
        the_tol_f: f64,
    ) -> Self {
        let mut e = EPCOfExtPC2d {
            done: false,
            init: false,
            nbsample: 0,
            umin: 0.0,
            usup: 0.0,
            tolu: 0.0,
            tolf: 0.0,
            f: PCFOfEPCOfExtPC2d::with_point_curve(the_p, the_c),
        };
        e.initialize_range(the_c, the_nb_sample, the_umin, the_usup, the_tol_u, the_tol_f);
        e.perform(the_p);
        e
    }

    /// OCCT Initialize(theC, theNbU, theTolU, theTolF) (L98-110).
    pub fn initialize_c(
        &mut self,
        the_c: &'a dyn Curve2dAdaptor,
        the_nb_u: i32,
        the_tol_u: f64,
        the_tol_f: f64,
    ) {
        self.init = true;
        self.nbsample = the_nb_u;
        self.tolu = the_tol_u;
        self.tolf = the_tol_f;
        self.f.initialize(the_c);
        self.umin = curve2d_tool::first_parameter(the_c);
        self.usup = curve2d_tool::last_parameter(the_c);
    }

    /// OCCT Initialize(theC, theNbU, theUmin, theUsup, theTolU, theTolF)
    /// (L119-133).
    pub fn initialize_range(
        &mut self,
        the_c: &'a dyn Curve2dAdaptor,
        the_nb_u: i32,
        the_umin: f64,
        the_usup: f64,
        the_tol_u: f64,
        the_tol_f: f64,
    ) {
        self.init = true;
        self.nbsample = the_nb_u;
        self.tolu = the_tol_u;
        self.tolf = the_tol_f;
        self.f.initialize(the_c);
        self.umin = the_umin;
        self.usup = the_usup;
    }

    /// OCCT Initialize(theNbU, theUmin, theUsup, theTolU, theTolF) (L141-152).
    pub fn initialize_params(
        &mut self,
        the_nb_u: i32,
        the_umin: f64,
        the_usup: f64,
        the_tol_u: f64,
        the_tol_f: f64,
    ) {
        self.nbsample = the_nb_u;
        self.tolu = the_tol_u;
        self.tolf = the_tol_f;
        self.umin = the_umin;
        self.usup = the_usup;
    }

    /// OCCT Initialize(theC) (L156).
    pub fn initialize_curve(&mut self, the_c: &'a dyn Curve2dAdaptor) {
        self.f.initialize(the_c);
    }

    /// OCCT Perform(theP) (L160-173).
    pub fn perform(&mut self, the_p: DVec2) {
        self.f.set_point(the_p);
        self.f.sub_interval_initialize(self.umin, self.usup);
        self.done = false;

        // math_FunctionRoots S(myF, myumin, myusup, mynbsample,
        //                      mytolu, mytolF, mytolF) — K defaults to 0.
        let s = FunctionRoots::new(
            &mut self.f,
            self.umin,
            self.usup,
            self.nbsample,
            self.tolu,
            self.tolf,
            self.tolf,
            0.0,
        );
        if !s.is_done() || s.is_all_null() {
            return;
        }

        self.done = true;
    }

    /// OCCT IsDone (L176).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT NbExt (L180-187).
    pub fn nb_ext(&self) -> usize {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.f.nb_ext()
    }

    /// OCCT SquareDistance(theN) (L192-199).
    pub fn square_distance(&self, the_n: usize) -> f64 {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.f.square_distance(the_n)
    }

    /// OCCT IsMin(theN) (L204-211).
    pub fn is_min(&self, the_n: usize) -> bool {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.f.is_min(the_n)
    }

    /// OCCT Point(theN) (L216-223).
    pub fn point(&self, the_n: usize) -> &POnCurv2d {
        if the_n < 1 || the_n > self.nb_ext() {
            panic!("Standard_OutOfRange");
        }
        self.f.point(the_n)
    }
}

impl Default for EPCOfExtPC2d<'_> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::Curve2d;

            /// Anchor: extrema of a circle of radius 2 around the origin from the
    /// point (5, 0) — one minimum at u=0 (dist 3, min) and one maximum at
    /// u=pi (dist 7). OCCT reference computed analytically.
    #[test]
    fn epc_of_ext_pc2d_circle_min_max() {
        // gp_Circ2d with axis at origin, X/Y directions, radius 2.
        let c = Curve2d::Circle(rcad_kernel::geom::Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: 2.0,
        });
        let p = DVec2::new(5.0, 0.0);
        // Contap_HContTool::Project parameters: epsX = 1e-8, Nbu = 20, Tol = 1e-5.
        let extrema = EPCOfExtPC2d::with_full_domain(p, &c, 20, 1.0e-8, 1.0e-5);

        assert!(extrema.is_done());
        assert!(extrema.nb_ext() >= 2);
        let mut dmin = f64::MAX;
        let mut dmax = 0.0f64;
        for i in 1..=extrema.nb_ext() {
            let d = extrema.square_distance(i);
            if extrema.is_min(i) {
                dmin = dmin.min(d);
            } else {
                dmax = dmax.max(d);
            }
        }
        // Minimum square distance = 3^2 = 9 (at u=0), maximum = 7^2 = 49.
        assert!((dmin - 9.0).abs() < 1e-6, "dmin={}", dmin);
        assert!((dmax - 49.0).abs() < 1e-6, "dmax={}", dmax);
        // The minimum point is the nearest extremum; check its parameter.
        for i in 1..=extrema.nb_ext() {
            if extrema.is_min(i) && extrema.square_distance(i) < 9.0 + 1e-6 {
                let pt = extrema.point(i);
                assert!((pt.parameter().rem_euclid(std::f64::consts::TAU) - 0.0).abs() < 1e-4
                    || (pt.parameter().rem_euclid(std::f64::consts::TAU)
                        - std::f64::consts::TAU)
                    .abs()
                        < 1e-4);
                assert!((pt.value().x - 2.0).abs() < 1e-4);
                assert!(pt.value().y.abs() < 1e-4);
            }
        }
    }

    /// Anchor: point projected onto a bounded line — the single extremum is
    /// the orthogonal projection (minimum). Like the OCCT Contap usage, the
    /// arc has a finite parameter range (face boundary pcurves are trimmed).
    #[test]
    fn epc_of_ext_pc2d_line_projection() {
        // Line y = 0 parameterized as u -> (u, 0), restricted to [0, 6].
        let c = Curve2d::Line(rcad_kernel::geom::Line2d {
            origin: DVec2::ZERO,
            direction: DVec2::new(1.0, 0.0),
        });
        let p = DVec2::new(3.0, 4.0);
        let extrema = EPCOfExtPC2d::with_range(p, &c, 20, 0.0, 6.0, 1.0e-8, 1.0e-5);
        assert!(extrema.is_done());
        assert!(extrema.nb_ext() >= 1);
        // Closest extremum is at u=3, square distance 16.
        let mut best = (f64::MAX, 0.0);
        for i in 1..=extrema.nb_ext() {
            let d = extrema.square_distance(i);
            if d < best.0 {
                best = (d, extrema.point(i).parameter());
            }
        }
        assert!((best.0 - 16.0).abs() < 1e-6, "d={}", best.0);
        assert!((best.1 - 3.0).abs() < 1e-6, "u={}", best.1);
    }

    /// Anchor: restricted parameter range — extrema are searched only in
    /// [umin, usup] (Initialize(theC, NbU, Umin, Usup, TolU, TolF)). On the
    /// circle from the previous anchor, restricting the range to [pi/2,
    /// 3pi/2] leaves only the maximum extremum at u=pi (square distance 49).
    #[test]
    fn epc_of_ext_pc2d_restricted_range() {
        let c = Curve2d::Circle(rcad_kernel::geom::Circle2d {
            center: DVec2::ZERO,
            x_dir: DVec2::new(1.0, 0.0),
            y_dir: DVec2::new(0.0, 1.0),
            radius: 2.0,
        });
        let p = DVec2::new(5.0, 0.0);
        let half_pi = std::f64::consts::FRAC_PI_2;
        let extrema = EPCOfExtPC2d::with_range(p, &c, 20, half_pi, 3.0 * half_pi, 1.0e-8, 1.0e-5);
        assert!(extrema.is_done());
        assert_eq!(extrema.nb_ext(), 1);
        assert!(!extrema.is_min(1));
        assert!((extrema.square_distance(1) - 49.0).abs() < 1e-6);
        let pt = extrema.point(1);
        assert!((pt.parameter() - std::f64::consts::PI).abs() < 1e-4);
        assert!((pt.value().x + 2.0).abs() < 1e-4);
        assert!(pt.value().y.abs() < 1e-4);
    }
}
