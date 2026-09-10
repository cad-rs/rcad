// OCCT TKMath ports for the Gradient minimization chain:
//   math_MultipleVarFunctionWithGradient (hxx) — the trait,
//   math_BracketMinimum (hxx/.cxx L17-305 + lxx) — the bracketing class,
//   math_BFGS (hxx/.cxx L17-510 + lxx) — the BFGS minimizer class with the
//     DirFunction helper and the ComputeInitScale / ComputeMinMaxScale /
//     MinimizeDirection statics.
//
// math_BrentMinimum is NOT duplicated here: the existing
// math::opt::BrentMinimum is the OCCT math_BrentMinimum port and is reused
// by MinimizeDirection exactly as math_BFGS.cxx uses math_BrentMinimum.
//
// Rust note: OCCT's DirFunction aliases F through a pointer while Perform
// also calls F directly; the borrow checker forbids that aliasing, so the
// DirFunction here carries only the line scratch vectors and every method
// takes `&mut F` as a parameter (identical call graph, no aliasing).

use crate::math::math_matrix::{Matrix, Vector};
use crate::math::opt::BrentMinimum;
use crate::math::root::FunctionValue;

// ---------------------------------------------------------------------------
// math_MultipleVarFunctionWithGradient
// ---------------------------------------------------------------------------

/// OCCT math_MultipleVarFunction (hxx) — the base trait.
pub trait MultipleVarFunction {
    /// OCCT NbVariables().
    fn nb_variables(&self) -> i32;
    /// OCCT Value(X, F).
    fn value(&mut self, x: &Vector, f: &mut f64) -> bool;
}

/// OCCT math_MultipleVarFunctionWithGradient
/// (math_MultipleVarFunctionWithGradient.hxx L8-32).
pub trait MultipleVarFunctionWithGradient: MultipleVarFunction {
    /// OCCT Gradient(X, G) — the gradient G of the function at X.
    fn gradient(&mut self, x: &Vector, g: &mut Vector) -> bool;
    /// OCCT Values(X, F, G) — the value F and the gradient G at X.
    fn values(&mut self, x: &Vector, f: &mut f64, g: &mut Vector) -> bool;
}

// ---------------------------------------------------------------------------
// math_BracketMinimum
// ---------------------------------------------------------------------------

/// OCCT GOLD (math_BracketMinimum.cxx L21).
const GOLD: f64 = 1.618034;
/// OCCT GLIMIT (math_BracketMinimum.cxx L23).
const GLIMIT: f64 = 100.0;
/// OCCT TINY (math_BracketMinimum.cxx L24).
const TINY: f64 = 1.0e-20;

/// OCCT scalar math_Function adapter — a closure over (x, &mut f) -> bool.
pub trait ScalarFunction {
    fn value(&mut self, x: f64, f: &mut f64) -> bool;
}

impl<T: FnMut(f64, &mut f64) -> bool> ScalarFunction for T {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        self(x, f)
    }
}

/// OCCT math_BracketMinimum (math_BracketMinimum.hxx L46-133) — computes a
/// bracketing triplet of abscissae Ax, Bx, Cx (such that Bx is between Ax
/// and Cx and F(Bx) is less than both F(Ax) and F(Cx)); the Brent
/// minimization is then done on the function F.
#[derive(Debug, Clone)]
pub struct BracketMinimum {
    done: bool,
    ax: f64,
    bx: f64,
    cx: f64,
    fax: f64,
    fbx: f64,
    fcx: f64,
    my_left: f64,
    my_right: f64,
    my_is_limited: bool,
    my_fa: bool,
    my_fb: bool,
}

impl BracketMinimum {
    /// OCCT math_BracketMinimum(A, B) (lxx L8-25) — prepares A and B only;
    /// Perform must be called explicitly.
    pub fn new(a: f64, b: f64) -> Self {
        BracketMinimum {
            done: false,
            ax: a,
            bx: b,
            cx: 0.0,
            fax: 0.0,
            fbx: 0.0,
            fcx: 0.0,
            my_left: -f64::INFINITY,
            my_right: f64::INFINITY,
            my_is_limited: false,
            my_fa: false,
            my_fb: false,
        }
    }

    /// OCCT SetLimits(theLeft, theRight) (lxx L33-38).
    pub fn set_limits(&mut self, the_left: f64, the_right: f64) {
        self.my_left = the_left;
        self.my_right = the_right;
        self.my_is_limited = true;
    }

    /// OCCT SetFA(theValue) (lxx L40-44).
    pub fn set_fa(&mut self, the_value: f64) {
        self.fax = the_value;
        self.my_fa = true;
    }

    /// OCCT SetFB(theValue) (lxx L46-50).
    pub fn set_fb(&mut self, the_value: f64) {
        self.fbx = the_value;
        self.my_fb = true;
    }

    /// OCCT Limited(theValue) (lxx L65-68).
    fn limited(&self, the_value: f64) -> f64 {
        if the_value < self.my_left {
            self.my_left
        } else if the_value > self.my_right {
            self.my_right
        } else {
            the_value
        }
    }

    /// OCCT LimitAndMayBeSwap(F, theA, theB, theFB, theC, theFC)
    /// (math_BracketMinimum.cxx L34-58). OCCT Precision::PConfusion().
    fn limit_and_may_be_swap(
        &self,
        f: &mut dyn ScalarFunction,
        the_a: f64,
        the_b: &mut f64,
        the_fb: &mut f64,
        the_c: &mut f64,
        the_fc: &mut f64,
    ) -> bool {
        const P_CONFUSION: f64 = 1.0e-12;
        *the_c = self.limited(*the_c);
        if (*the_b - *the_c).abs() < P_CONFUSION {
            return false;
        }
        let ok = f.value(*the_c, the_fc);
        if !ok {
            return false;
        }
        // check that B is between A and C
        if (the_a - *the_b) * (*the_b - *the_c) < 0.0 {
            // swap B and C — OCCT SHFT(dum, theB, theC, dum) x2.
            std::mem::swap(the_b, the_c);
            std::mem::swap(the_fb, the_fc);
        }
        true
    }

    /// OCCT Perform(F) (math_BracketMinimum.cxx L60-222).
    pub fn perform(&mut self, f: &mut dyn ScalarFunction) {
        self.done = false;
        let lambda = GOLD;
        if !self.my_fa {
            let ok = f.value(self.ax, &mut self.fax);
            if !ok {
                return;
            }
        }
        if !self.my_fb {
            let ok = f.value(self.bx, &mut self.fbx);
            if !ok {
                return;
            }
        }
        // OCCT: if (FBx > FAx) { SHFT(dum, Ax, Bx, dum); SHFT(dum, FBx, FAx, dum); }
        // — swap Ax/Bx and FAx/FBx.
        if self.fbx > self.fax {
            std::mem::swap(&mut self.ax, &mut self.bx);
            std::mem::swap(&mut self.fax, &mut self.fbx);
        }
        // get next prob after (A, B)
        self.cx = self.bx + lambda * (self.bx - self.ax);
        if self.my_is_limited {
            let (mut bx, mut fbx, mut cx, mut fcx) = (self.bx, self.fbx, self.cx, self.fcx);
            let ok = self.limit_and_may_be_swap(f, self.ax, &mut bx, &mut fbx, &mut cx, &mut fcx);
            self.bx = bx;
            self.fbx = fbx;
            self.cx = cx;
            self.fcx = fcx;
            if !ok {
                return;
            }
        } else {
            let ok = f.value(self.cx, &mut self.fcx);
            if !ok {
                return;
            }
        }
        while self.fbx > self.fcx {
            let r = (self.bx - self.ax) * (self.fbx - self.fcx);
            let q = (self.bx - self.cx) * (self.fbx - self.fax);
            // OCCT SIGN(MAX(fabs(q - r), TINY), q - r).
            let m = (q - r).abs().max(TINY);
            let denom = 2.0 * if q - r > 0.0 { m } else { -m };
            let mut u = self.bx - ((self.bx - self.cx) * q - (self.bx - self.ax) * r) / denom;
            let mut ulim = self.bx + GLIMIT * (self.cx - self.bx);
            if self.my_is_limited {
                ulim = self.limited(ulim);
            }
            let mut fu = 0.0;
            if (self.bx - u) * (u - self.cx) > 0.0 {
                // u is between B and C
                let ok = f.value(u, &mut fu);
                if !ok {
                    return;
                }
                if fu < self.fcx {
                    // solution is found (B, u, c)
                    self.ax = self.bx;
                    self.bx = u;
                    self.fax = self.fbx;
                    self.fbx = fu;
                    self.done = true;
                    return;
                } else if fu > self.fbx {
                    // solution is found (A, B, u)
                    self.cx = u;
                    self.fcx = fu;
                    self.done = true;
                    return;
                }
                // get next prob after (B, C)
                u = self.cx + lambda * (self.cx - self.bx);
                if self.my_is_limited {
                    let (mut bx, mut cx, mut fcx, mut uo, mut fuo) =
                        (self.bx, self.cx, self.fcx, u, fu);
                    let ok =
                        self.limit_and_may_be_swap(f, bx, &mut cx, &mut fcx, &mut uo, &mut fuo);
                    self.bx = bx;
                    self.cx = cx;
                    self.fcx = fcx;
                    u = uo;
                    fu = fuo;
                    if !ok {
                        return;
                    }
                } else {
                    let ok = f.value(u, &mut fu);
                    if !ok {
                        return;
                    }
                }
            } else if (self.cx - u) * (u - ulim) > 0.0 {
                // u is beyond C but between C and limit
                let ok = f.value(u, &mut fu);
                if !ok {
                    return;
                }
            } else if (u - ulim) * (ulim - self.cx) >= 0.0 {
                // u is beyond limit
                u = ulim;
                let ok = f.value(u, &mut fu);
                if !ok {
                    return;
                }
            } else {
                // u tends to approach to the side of A,
                // so reset it to the next prob after (B, C)
                u = self.cx + GOLD * (self.cx - self.bx);
                if self.my_is_limited {
                    let (mut bx, mut cx, mut fcx, mut uo, mut fuo) =
                        (self.bx, self.cx, self.fcx, u, fu);
                    let ok =
                        self.limit_and_may_be_swap(f, bx, &mut cx, &mut fcx, &mut uo, &mut fuo);
                    self.bx = bx;
                    self.cx = cx;
                    self.fcx = fcx;
                    u = uo;
                    fu = fuo;
                    if !ok {
                        return;
                    }
                } else {
                    let ok = f.value(u, &mut fu);
                    if !ok {
                        return;
                    }
                }
            }
            // SHFT(Ax, Bx, Cx, u); SHFT(FAx, FBx, FCx, fu);
            self.ax = self.bx;
            self.bx = self.cx;
            self.cx = u;
            self.fax = self.fbx;
            self.fbx = self.fcx;
            self.fcx = fu;
        }
        self.done = true;
    }

    /// OCCT IsDone() (lxx L27-30).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Values(A, B, C) (cxx L283-291).
    pub fn values(&self) -> (f64, f64, f64) {
        assert!(self.done, "StdFail_NotDone: math_BracketMinimum::Values");
        (self.ax, self.bx, self.cx)
    }

    /// OCCT FunctionValues(FA, FB, FC) (cxx L293-301).
    pub fn function_values(&self) -> (f64, f64, f64) {
        assert!(
            self.done,
            "StdFail_NotDone: math_BracketMinimum::FunctionValues"
        );
        (self.fax, self.fbx, self.fcx)
    }
}

// ---------------------------------------------------------------------------
// math_BFGS
// ---------------------------------------------------------------------------

/// OCCT math_Status values used by math_BFGS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MathStatus {
    Ok,
    SingularMatrix,
    ArgumentError,
    NoConvergence,
    FunctionError,
    DirectionSearchError,
    TooManyIterations,
    UserAborted,
}

/// OCCT RealSmall().
const REAL_SMALL: f64 = 2.225_073_858_507_2014e-308;
/// OCCT Precision::PConfusion().
const P_CONFUSION: f64 = 1.0e-12;
/// OCCT Precision::Infinite().
const INFINITE: f64 = 1.0e100;

/// OCCT DirFunction (math_BFGS.cxx L20-80) — the scratch vectors of the
/// line function F(P0 + x*Dir). The OCCT class also aliases F; here F is
/// passed per-call (no aliasing under the borrow checker).
struct DirFunction<'a> {
    p0: &'a mut Vector,
    dir: &'a mut Vector,
    p: &'a mut Vector,
    g: &'a mut Vector,
}

impl DirFunction<'_> {
    /// OCCT Initialize(p0, dir) (cxx L46-50).
    fn initialize(&mut self, p0: &Vector, dir: &Vector) {
        for i in self.p0.lower..=self.p0.upper() {
            let v = p0.get(i);
            self.p0.set(i, v);
        }
        for i in self.dir.lower..=self.dir.upper() {
            let v = dir.get(i);
            self.dir.set(i, v);
        }
    }

    /// OCCT Value(x, fval) (cxx L56-63): *P = *Dir; P->Multiply(x);
    /// P->Add(*P0).
    fn line_value<F: MultipleVarFunctionWithGradient>(
        &mut self,
        f: &mut F,
        x: f64,
        fval: &mut f64,
    ) -> bool {
        for i in self.p.lower..=self.p.upper() {
            let v = self.dir.get(i) * x + self.p0.get(i);
            self.p.set(i, v);
        }
        *fval = 0.0;
        f.value(self.p, fval)
    }

    /// OCCT Values(x, fval, D) (cxx L65-76).
    fn line_values<F: MultipleVarFunctionWithGradient>(
        &mut self,
        f: &mut F,
        x: f64,
        fval: &mut f64,
        d: &mut f64,
    ) -> bool {
        for i in self.p.lower..=self.p.upper() {
            let v = self.dir.get(i) * x + self.p0.get(i);
            self.p.set(i, v);
        }
        *fval = 0.0;
        *d = 0.0;
        if f.values(self.p, fval, self.g) {
            // D = (*G).Multiplied(*Dir).
            let mut dot = 0.0;
            for i in self.g.lower..=self.g.upper() {
                dot += self.g.get(i) * self.dir.get(i);
            }
            *d = dot;
            return true;
        }
        false
    }
}

/// ScalarFunction adapter over DirFunction + F (OCCT passes the DirFunction
/// itself; math_FunctionWithDerivative IS-A math_Function).
struct DirScalar<'a, 'b, F: MultipleVarFunctionWithGradient> {
    d: &'a mut DirFunction<'b>,
    f: &'a mut F,
}

impl<F: MultipleVarFunctionWithGradient> ScalarFunction for DirScalar<'_, '_, F> {
    fn value(&mut self, x: f64, f: &mut f64) -> bool {
        self.d.line_value(self.f, x, f)
    }
}

/// FunctionValue adapter for math_BrentMinimum::Perform.
struct DirFunctionValue<'a, 'b, F: MultipleVarFunctionWithGradient> {
    d: &'a mut DirFunction<'b>,
    f: &'a mut F,
}

impl<F: MultipleVarFunctionWithGradient> FunctionValue for DirFunctionValue<'_, '_, F> {
    fn value(&mut self, x: f64) -> Option<f64> {
        let mut fval = 0.0;
        if self.d.line_value(self.f, x, &mut fval) {
            Some(fval)
        } else {
            None
        }
    }
}

/// OCCT static ComputeInitScale (math_BFGS.cxx L82-96). OCCT RealSmall().
fn compute_init_scale(the_f0: f64, the_dir: &Vector, the_gr: &Vector, the_scale: &mut f64) -> bool {
    let mut dy1 = 0.0;
    for i in the_dir.lower..=the_dir.upper() {
        dy1 += the_gr.get(i) * the_dir.get(i);
    }
    if dy1.abs() < REAL_SMALL {
        return false;
    }
    let mut hnr1 = 0.0;
    for i in the_dir.lower..=the_dir.upper() {
        hnr1 += the_dir.get(i) * the_dir.get(i);
    }
    let alfa = 0.7 * (-the_f0) / dy1;
    *the_scale = 0.015 / hnr1.sqrt();
    if *the_scale > alfa {
        *the_scale = alfa;
    }
    true
}

/// OCCT static ComputeMinMaxScale (math_BFGS.cxx L98-176).
fn compute_min_max_scale(
    the_point: &Vector,
    the_dir: &Vector,
    the_left: &Vector,
    the_right: &Vector,
    the_min_scale: &mut f64,
    the_max_scale: &mut f64,
) -> bool {
    for an_idx in the_left.lower..=the_left.upper() {
        let a_left = the_left.get(an_idx) - the_point.get(an_idx);
        let a_right = the_right.get(an_idx) - the_point.get(an_idx);
        if the_dir.get(an_idx).abs() > REAL_SMALL {
            // Use PConfusion to get off a little from the bounds to prevent
            // possible refuse in Value function.
            let a_lscale = (a_left + P_CONFUSION) / the_dir.get(an_idx);
            let a_rscale = (a_right - P_CONFUSION) / the_dir.get(an_idx);
            if a_left.abs() < P_CONFUSION {
                // Point is on the left border.
                *the_max_scale = (*the_max_scale).min(0.0f64.max(a_rscale));
                *the_min_scale = (*the_min_scale).max(0.0f64.min(a_rscale));
            } else if a_right.abs() < P_CONFUSION {
                // Point is on the right border.
                *the_max_scale = (*the_max_scale).min(0.0f64.max(a_lscale));
                *the_min_scale = (*the_min_scale).max(0.0f64.min(a_lscale));
            } else if a_left * a_right < 0.0 {
                // Point is inside allowed range.
                *the_max_scale = (*the_max_scale).min(a_lscale.max(a_rscale));
                *the_min_scale = (*the_min_scale).max(a_lscale.min(a_rscale));
            } else {
                // point is out of bounds
                return false;
            }
        } else {
            // Direction is parallel to the border.
            // Check that the point is not out of bounds
            if a_left > P_CONFUSION || a_right < -P_CONFUSION {
                return false;
            }
        }
    }
    true
}

/// OCCT static MinimizeDirection (math_BFGS.cxx L178-321).
#[allow(clippy::too_many_arguments)]
fn minimize_direction<F: MultipleVarFunctionWithGradient>(
    p: &mut Vector,
    f0: f64,
    gr: &Vector,
    dir: &mut Vector,
    result: &mut f64,
    f_dir: &mut DirFunction,
    f: &mut F,
    is_bounds: bool,
    the_left: &Vector,
    the_right: &Vector,
) -> bool {
    let mut lambda = 0.0;
    if !compute_init_scale(f0, dir, gr, &mut lambda) {
        return false;
    }
    // by default the scaling range is unlimited
    let mut a_min_lambda = -INFINITE;
    let mut a_max_lambda = INFINITE;
    if is_bounds {
        // limit the scaling range taking into account the bounds
        if !compute_min_max_scale(p, dir, the_left, the_right, &mut a_min_lambda, &mut a_max_lambda)
        {
            return false;
        }
        if a_min_lambda > -P_CONFUSION && a_max_lambda < P_CONFUSION {
            // Point is on the border and the direction shows outside.
            // Make direction to go along the border
            for an_idx in the_left.lower..=the_left.upper() {
                if ((p.get(an_idx) - the_right.get(an_idx)).abs() < P_CONFUSION
                    && dir.get(an_idx) > 0.0)
                    || ((p.get(an_idx) - the_left.get(an_idx)).abs() < P_CONFUSION
                        && dir.get(an_idx) < 0.0)
                {
                    dir.set(an_idx, 0.0);
                }
            }
            // re-compute scale values with new direction
            if !compute_init_scale(f0, dir, gr, &mut lambda) {
                return false;
            }
            if !compute_min_max_scale(
                p,
                dir,
                the_left,
                the_right,
                &mut a_min_lambda,
                &mut a_max_lambda,
            ) {
                return false;
            }
        }
        lambda = lambda.min(a_max_lambda);
        lambda = lambda.max(a_min_lambda);
    }
    f_dir.initialize(p, dir);
    let mut f1 = 0.0;
    if !f_dir.line_value(f, lambda, &mut f1) {
        return false;
    }
    let mut bracket = BracketMinimum::new(0.0, lambda);
    if is_bounds {
        bracket.set_limits(a_min_lambda, a_max_lambda);
    }
    bracket.set_fa(f0);
    bracket.set_fb(f1);
    let mut scalar = DirScalar { d: f_dir, f };
    bracket.perform(&mut scalar);
    if bracket.is_done() {
        // find minimum inside the bracket
        let (ax, xx, bx) = bracket.values();
        let (_fax, fxx, _fbx) = bracket.function_values();
        let niter = 100;
        let tol = 1.0e-03;
        // math_BrentMinimum Sol(tol, Fxx, niter, 1.e-08) — the existing
        // math::opt::BrentMinimum is the math_BrentMinimum port.
        let mut sol = BrentMinimum::new_with_fbx(tol, fxx, niter, 1.0e-08);
        let mut sol_f = DirFunctionValue { d: f_dir, f };
        sol.perform(&mut sol_f, ax, xx, bx);
        if sol.is_done() {
            let scale = sol.location();
            *result = sol.minimum();
            for i in dir.lower..=dir.upper() {
                let v = dir.get(i) * scale;
                dir.set(i, v);
            }
            for i in p.lower..=p.upper() {
                let v = p.get(i) + dir.get(i);
                p.set(i, v);
            }
            return true;
        }
    } else if is_bounds {
        // Bracket definition is failure. If the bounds are defined then
        // set current point to intersection with bounds
        let mut a_fmin = 0.0;
        let mut a_fmax = 0.0;
        if !f_dir.line_value(f, a_min_lambda, &mut a_fmin) {
            return false;
        }
        if !f_dir.line_value(f, a_max_lambda, &mut a_fmax) {
            return false;
        }
        let a_best_lambda;
        if a_fmin < a_fmax {
            a_best_lambda = a_min_lambda;
            *result = a_fmin;
        } else {
            a_best_lambda = a_max_lambda;
            *result = a_fmax;
        }
        for i in dir.lower..=dir.upper() {
            let v = dir.get(i) * a_best_lambda;
            dir.set(i, v);
        }
        for i in p.lower..=p.upper() {
            let v = p.get(i) + dir.get(i);
            p.set(i, v);
        }
        return true;
    }
    false
}

/// OCCT math_BFGS (math_BFGS.hxx L38-115).
#[derive(Debug, Clone)]
pub struct Bfgs {
    /// OCCT TheStatus.
    the_status: MathStatus,
    /// OCCT TheLocation.
    the_location: Vector,
    /// OCCT TheGradient.
    the_gradient: Vector,
    /// OCCT PreviousMinimum / TheMinimum / XTol / EPSZ.
    previous_minimum: f64,
    the_minimum: f64,
    xtol: f64,
    epsz: f64,
    /// OCCT nbiter / Itermax / Done.
    nbiter: i32,
    itermax: i32,
    done: bool,
    /// OCCT myIsBoundsDefined / myLeft / myRight.
    my_is_bounds_defined: bool,
    my_left: Vector,
    my_right: Vector,
}

impl Bfgs {
    /// OCCT math_BFGS(NbVariables, Tolerance = 1.0e-8, NbIterations = 200,
    /// ZEPS = 1.0e-12) (cxx L408-427).
    pub fn new(nb_variables: i32, tolerance: f64, nb_iterations: i32, zeps: f64) -> Self {
        Bfgs {
            the_status: MathStatus::Ok,
            the_location: Vector::new(1, nb_variables),
            the_gradient: Vector::new(1, nb_variables),
            previous_minimum: 0.0,
            the_minimum: 0.0,
            xtol: tolerance,
            epsz: zeps,
            nbiter: 0,
            itermax: nb_iterations,
            done: false,
            my_is_bounds_defined: false,
            my_left: Vector::new_init(1, nb_variables, 0.0),
            my_right: Vector::new_init(1, nb_variables, 0.0),
        }
    }

    /// OCCT SetBoundary(theLeftBorder, theRightBorder) (cxx L500-506).
    pub fn set_boundary(&mut self, the_left_border: &Vector, the_right_border: &Vector) {
        for i in self.my_left.lower..=self.my_left.upper() {
            let l = the_left_border.get(i);
            let r = the_right_border.get(i);
            self.my_left.set(i, l);
            self.my_right.set(i, r);
        }
        self.my_is_bounds_defined = true;
    }

    /// OCCT Perform(F, StartingPoint) (math_BFGS.cxx L323-445) — uses the
    /// base IsSolutionReached test (cxx L447-452).
    pub fn perform<F: MultipleVarFunctionWithGradient>(
        &mut self,
        f: &mut F,
        starting_point: &Vector,
    ) {
        // The OCCT virtual IsSolutionReached dispatch maps to a checker
        // closure capturing the base-class tolerances.
        let xtol = self.xtol;
        let epsz = self.epsz;
        self.perform_with_checker(f, starting_point, &mut |the_minimum, previous_minimum, _f| {
            2.0 * (the_minimum - previous_minimum).abs()
                <= xtol * (the_minimum.abs() + previous_minimum.abs() + epsz)
        });
    }

    /// OCCT Perform body with the virtual IsSolutionReached dispatched
    /// through `is_solution(the_minimum, previous_minimum, F)` — the
    /// Gradient_BFGS sub-class override plugs its own test here (the
    /// sub-class test also reads F's MaxError3d/MaxError2d, hence the F
    /// argument).
    pub fn perform_with_checker<F: MultipleVarFunctionWithGradient>(
        &mut self,
        f: &mut F,
        starting_point: &Vector,
        is_solution: &mut dyn FnMut(f64, f64, &mut F) -> bool,
    ) {
        let n = self.the_location.length();
        let mut xi = Vector::new(1, n);
        let mut dg = Vector::new(1, n);
        let mut hdg = Vector::new(1, n);
        let mut hessin = Matrix::new_init(1, n, 1, n, 0.0);
        let mut temp1 = Vector::new(1, n);
        let mut temp2 = Vector::new(1, n);
        let mut temp3 = Vector::new(1, n);
        let mut temp4 = Vector::new(1, n);
        let mut f_dir = DirFunction {
            p0: &mut temp1,
            dir: &mut temp2,
            p: &mut temp3,
            g: &mut temp4,
        };
        for i in self.the_location.lower..=self.the_location.upper() {
            let v = starting_point.get(i);
            self.the_location.set(i, v);
        }
        let mut previous_minimum = 0.0;
        let mut gradient = Vector::new(1, n);
        let good = f.values(&self.the_location, &mut previous_minimum, &mut gradient);
        for i in self.the_gradient.lower..=self.the_gradient.upper() {
            let v = gradient.get(i);
            self.the_gradient.set(i, v);
        }
        self.previous_minimum = previous_minimum;
        if !good {
            self.done = false;
            self.the_status = MathStatus::FunctionError;
            return;
        }
        for i in 1..=n {
            hessin.set(i, i, 1.0);
            let v = -self.the_gradient.get(i);
            xi.set(i, v);
        }
        self.nbiter = 1;
        while self.nbiter <= self.itermax {
            self.the_minimum = self.previous_minimum;
            let is_good = minimize_direction(
                &mut self.the_location,
                self.the_minimum,
                &self.the_gradient,
                &mut xi,
                &mut self.the_minimum,
                &mut f_dir,
                f,
                self.my_is_bounds_defined,
                &self.my_left,
                &self.my_right,
            );
            if is_solution(self.the_minimum, self.previous_minimum, f) {
                self.done = true;
                self.the_status = MathStatus::Ok;
                return;
            }
            if !is_good {
                self.done = false;
                self.the_status = MathStatus::DirectionSearchError;
                return;
            }
            self.previous_minimum = self.the_minimum;
            // dg = TheGradient (the values before the line minimization).
            for i in 1..=n {
                let v = self.the_gradient.get(i);
                dg.set(i, v);
            }
            let mut the_minimum = 0.0;
            let mut gradient = Vector::new(1, n);
            let good = f.values(&self.the_location, &mut the_minimum, &mut gradient);
            for i in self.the_gradient.lower..=self.the_gradient.upper() {
                let v = gradient.get(i);
                self.the_gradient.set(i, v);
            }
            self.the_minimum = the_minimum;
            if !good {
                self.done = false;
                self.the_status = MathStatus::FunctionError;
                return;
            }
            for i in 1..=n {
                let v = self.the_gradient.get(i) - dg.get(i);
                dg.set(i, v);
            }
            for i in 1..=n {
                hdg.set(i, 0.0);
                for j in 1..=n {
                    let v = hdg.get(i) + hessin.get(i, j) * dg.get(j);
                    hdg.set(i, v);
                }
            }
            let mut fac = 0.0;
            let mut fae = 0.0;
            for i in 1..=n {
                fac += dg.get(i) * xi.get(i);
                fae += dg.get(i) * hdg.get(i);
            }
            fac = 1.0 / fac;
            let fad = 1.0 / fae;
            for i in 1..=n {
                let v = fac * xi.get(i) - fad * hdg.get(i);
                dg.set(i, v);
            }
            for i in 1..=n {
                for j in 1..=n {
                    let v = hessin.get(i, j)
                        + fac * xi.get(i) * xi.get(j)
                        - fad * hdg.get(i) * hdg.get(j)
                        + fae * dg.get(i) * dg.get(j);
                    hessin.set(i, j, v);
                }
            }
            for i in 1..=n {
                xi.set(i, 0.0);
                for j in 1..=n {
                    let v = xi.get(i) - hessin.get(i, j) * self.the_gradient.get(j);
                    xi.set(i, v);
                }
            }
            self.nbiter += 1;
        }
        self.done = false;
        self.the_status = MathStatus::TooManyIterations;
    }

    /// OCCT IsSolutionReached(F) (cxx L447-452) — called at the end of each
    /// iteration; redefined by the Gradient_BFGS shell for its specific
    /// test.
    pub fn is_solution_reached<F: MultipleVarFunctionWithGradient>(&self, _f: &F) -> bool {
        2.0 * (self.the_minimum - self.previous_minimum).abs()
            <= self.xtol * (self.the_minimum.abs() + self.previous_minimum.abs() + self.epsz)
    }

    /// OCCT IsDone() (lxx).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Location() (lxx).
    pub fn location(&self) -> &Vector {
        assert!(self.done, "StdFail_NotDone: math_BFGS::Location");
        &self.the_location
    }

    /// OCCT Location(Loc) (lxx).
    pub fn location_into(&self, loc: &mut Vector) {
        assert!(self.done, "StdFail_NotDone: math_BFGS::Location");
        for i in loc.lower..=loc.upper() {
            let v = self.the_location.get(i);
            loc.set(i, v);
        }
    }

    /// OCCT Minimum() (lxx).
    pub fn minimum(&self) -> f64 {
        assert!(self.done, "StdFail_NotDone: math_BFGS::Minimum");
        self.the_minimum
    }

    /// OCCT Gradient() (lxx).
    pub fn gradient(&self) -> &Vector {
        assert!(self.done, "StdFail_NotDone: math_BFGS::Gradient");
        &self.the_gradient
    }

    /// OCCT Gradient(Grad) (lxx).
    pub fn gradient_into(&self, grad: &mut Vector) {
        assert!(self.done, "StdFail_NotDone: math_BFGS::Gradient");
        for i in grad.lower..=grad.upper() {
            let v = self.the_gradient.get(i);
            grad.set(i, v);
        }
    }

    /// OCCT NbIterations() (lxx).
    pub fn nb_iterations(&self) -> i32 {
        assert!(self.done, "StdFail_NotDone: math_BFGS::NbIterations");
        self.nbiter
    }

    /// OCCT protected TheMinimum read access (used by the Gradient_BFGS
    /// IsSolutionReached override).
    pub fn the_minimum_value(&self) -> f64 {
        self.the_minimum
    }

    /// OCCT protected PreviousMinimum read access.
    pub fn previous_minimum_value(&self) -> f64 {
        self.previous_minimum
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Quadratic bowl F = (x-1)^2 + 10*(y-2)^2 with analytic gradient —
    /// math_BFGS must drive both components to (1, 2).
    #[derive(Debug, Clone)]
    struct Bowl;
    impl MultipleVarFunction for Bowl {
        fn nb_variables(&self) -> i32 {
            2
        }
        fn value(&mut self, x: &Vector, f: &mut f64) -> bool {
            let dx = x.get(1) - 1.0;
            let dy = x.get(2) - 2.0;
            *f = dx * dx + 10.0 * dy * dy;
            true
        }
    }
    impl MultipleVarFunctionWithGradient for Bowl {
        fn gradient(&mut self, x: &Vector, g: &mut Vector) -> bool {
            g.set(1, 2.0 * (x.get(1) - 1.0));
            g.set(2, 20.0 * (x.get(2) - 2.0));
            true
        }
        fn values(&mut self, x: &Vector, f: &mut f64, g: &mut Vector) -> bool {
            self.value(x, f);
            self.gradient(x, g);
            true
        }
    }

    #[test]
    fn bfgs_quadratic_bowl() {
        let mut b = Bfgs::new(2, 1.0e-8, 200, 1.0e-12);
        let start = Vector::new_init(1, 2, 0.0);
        let mut f = Bowl;
        b.perform(&mut f, &start);
        assert!(b.is_done(), "BFGS must converge");
        assert!((b.location().get(1) - 1.0).abs() < 1.0e-6);
        assert!((b.location().get(2) - 2.0).abs() < 1.0e-6);
        assert!(b.minimum() < 1.0e-12);
    }
}
