// OCCT math_NewtonFunctionRoot (math_NewtonFunctionRoot.cxx) 1:1 Rust
// translation — Newton-Raphson root finding with bounds clamping and best-X
// tracking.
//
// OCCT ref: math_NewtonFunctionRoot.cxx (whole file); the "full range"
// constructor (F, Guess, EpsX, EpsF, NbIterations) used by
// math_TrigonometricFunctionRoots sets Binf=RealFirst, Bsup=RealLast.
//
// The file also carries the bounded math_FunctionRoot constructor
// (math_FunctionRoot.cxx L95-119) that wraps math_FunctionSetRoot — the
// form ChFi3d_FilBuilder::SplitSurf instantiates as
// `math_FunctionRoot Resol(Fonc, (a+c)/2, tol2d, a, c, 50)`.

use super::root::FunctionWithDerivative;
use crate::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};

/// OCCT math_MyFunctionSetWithDerivatives (math_FunctionRoot.cxx L27-71) —
/// the adapter of a 1-variable math_FunctionWithDerivative into the
/// math_FunctionSetWithDerivatives interface consumed by math_FunctionSetRoot.
struct MyFunctionSetWithDerivatives<'a> {
    /// OCCT math_FunctionWithDerivative* Ff.
    ff: &'a mut dyn FunctionWithDerivative,
}

impl<'a> MyFunctionSetWithDerivatives<'a> {
    /// OCCT math_MyFunctionSetWithDerivatives(F) (L43-46).
    fn new(f: &'a mut dyn FunctionWithDerivative) -> Self {
        MyFunctionSetWithDerivatives { ff: f }
    }
}

impl<'a> FunctionSetWithDerivatives for MyFunctionSetWithDerivatives<'a> {
    /// OCCT NbVariables() (L48-51) — returns 1.
    fn nb_variables(&self) -> usize {
        1
    }

    /// OCCT NbEquations() (L53-56) — returns 1.
    fn nb_equations(&self) -> usize {
        1
    }

    /// OCCT Value(X, Fs) (L58-61) — Ff->Value(X(1), Fs(1)).
    fn value(&mut self, x: &[f64], fs: &mut [f64]) -> bool {
        match self.ff.value(x[0]) {
            Some(v) => {
                fs[0] = v;
                true
            }
            None => false,
        }
    }

    /// OCCT Derivatives(X, D) (L63-66) — Ff->Derivative(X(1), D(1, 1)).
    fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        match self.ff.derivative(x[0]) {
            Some(v) => {
                d[0][0] = v;
                true
            }
            None => false,
        }
    }

    /// OCCT Values(X, F, D) (L68-71) — Ff->Values(X(1), F(1), D(1, 1)).
    fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        match self.ff.values(x[0]) {
            Some((fv, dv)) => {
                f[0] = fv;
                d[0][0] = dv;
                true
            }
            None => false,
        }
    }
}

/// OCCT math_NewtonFunctionRoot.
pub struct NewtonFunctionRoot {
    /// OCCT Binf / Bsup — clamping bounds.
    b_inf: f64,
    b_sup: f64,
    /// OCCT EpsilonX / EpsilonF.
    eps_x: f64,
    eps_f: f64,
    /// OCCT Itermax.
    itermax: i32,
    /// OCCT Done.
    done: bool,
    /// OCCT X — the found root.
    x: f64,
    /// OCCT Fx / DFx.
    _fx: f64,
    _dfx: f64,
    /// OCCT It.
    it: i32,
}

impl NewtonFunctionRoot {
    /// OCCT math_NewtonFunctionRoot(F, Guess, EpsX, EpsF, NbIterations) — the
    /// full-range constructor used by math_TrigonometricFunctionRoots.
    pub fn new_full_range(
        f: &mut dyn FunctionWithDerivative,
        guess: f64,
        eps_x: f64,
        eps_f: f64,
        nb_iterations: i32,
    ) -> Self {
        let mut s = NewtonFunctionRoot {
            b_inf: -f64::MAX,
            b_sup: f64::MAX,
            eps_x,
            eps_f,
            itermax: nb_iterations,
            done: false,
            x: f64::MAX,
            _fx: f64::MAX,
            _dfx: 0.0,
            it: 0,
        };
        s.perform(f, guess);
        s
    }

    /// OCCT Perform(F, Guess).
    fn perform(&mut self, f: &mut dyn FunctionWithDerivative, guess: f64) {
        let (aa, bb) = if self.b_inf < self.b_sup {
            (self.b_inf, self.b_sup)
        } else {
            (self.b_sup, self.b_inf)
        };

        let mut dx = f64::MAX;
        let mut fx = f64::MAX;
        self.x = guess;
        self.it = 1;

        // OCCT: the best estimate is tracked and returned even when the
        // iteration diverges.
        let mut best_x = self.x;
        let mut best_fx = f64::MAX;

        while self.it <= self.itermax && (dx.abs() > self.eps_x || fx.abs() > self.eps_f) {
            let ok = f.values(self.x);

            if let Some((fxv, dfxv)) = ok {
                // OCCT: F.Values(X, Fx, DFx) writes Fx before the step, so the
                // loop condition sees the value at the current X.
                fx = fxv;
                self._dfx = dfxv;
                let abs_fx = fxv.abs();
                if abs_fx < best_fx {
                    best_fx = abs_fx;
                    best_x = self.x;
                }

                if dfxv == 0.0 {
                    self.done = false;
                    self.it = self.itermax + 1;
                } else {
                    dx = fxv / dfxv;
                    self.x -= dx;
                    // Limit the variations of X.
                    if self.x <= aa {
                        self.x = aa;
                    }
                    if self.x >= bb {
                        self.x = bb;
                    }
                    self.it += 1;
                }
            } else {
                self.done = false;
                self.it = self.itermax + 1;
            }
        }

        self.x = best_x;
        self._fx = fx;
        self.done = self.it <= self.itermax;
    }

    /// OCCT math_FunctionRoot::math_FunctionRoot(F, Guess, Tolerance, A, B,
    /// NbIterations) (math_FunctionRoot.cxx L95-119) — the bounded form that
    /// wraps math_FunctionSetRoot with the (V, Aa, Bb, Tol) vectors.
    pub fn new_bounded(
        f: &mut dyn FunctionWithDerivative,
        guess: f64,
        tolerance: f64,
        a: f64,
        b: f64,
        nb_iterations: i32,
    ) -> Self {
        // OCCT L102: math_Vector V(1, 1), Aa(1, 1), Bb(1, 1), Tol(1, 1);
        let v = [guess];
        let aa = [a];
        let bb = [b];
        let tol = [tolerance];
        // OCCT L103: math_MyFunctionSetWithDerivatives Ff(F);
        // OCCT L104-107: V(1) = Guess; Tol(1) = Tolerance; Aa(1) = A;
        //                Bb(1) = B;
        // OCCT L108-109: math_FunctionSetRoot Sol(Ff, Tol, NbIterations);
        //                Sol.Perform(Ff, V, Aa, Bb);
        // (the OCCT Perform default theStopOnDivergent = Standard_False)
        let sol = {
            let mut ff = MyFunctionSetWithDerivatives::new(f);
            let mut sol = FunctionSetRoot::new(&ff, &tol, nb_iterations);
            sol.perform(&mut ff, &v, &aa, &bb, false);
            sol
        };
        let mut s = NewtonFunctionRoot {
            // The rcad struct merges the math_NewtonFunctionRoot member set
            // with this math_FunctionRoot constructor; the fields the bounded
            // OCCT class has no equivalent of mirror the constructor inputs.
            b_inf: a,
            b_sup: b,
            eps_x: tolerance,
            eps_f: 0.0,
            itermax: nb_iterations,
            done: false,
            x: f64::MAX,
            _fx: f64::MAX,
            _dfx: 0.0,
            it: 0,
        };
        // OCCT L110: Done = Sol.IsDone();
        s.done = sol.is_done();
        if s.done {
            // OCCT L113: NbIter = Sol.NbIterations();
            s.it = sol.nb_iterations();
            // OCCT L114: F.GetStateNumber();
            f.get_state_number();
            // OCCT L115: TheRoot = Sol.Root()(1);
            s.x = sol.root()[0];
            // OCCT L116: TheDerivative = Sol.Derivative()(1, 1);
            s._dfx = sol.derivative();
            // OCCT L117: F.Value(TheRoot, TheError);
            if let Some(err) = f.value(s.x) {
                s._fx = err;
            }
        }
        s
    }

    /// OCCT math_FunctionRoot::IsDone() (math_FunctionRoot.hxx L63).
    /// OCCT math_NewtonFunctionRoot::IsDone() — same accessor.
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT Root() — TheRoot.
    pub fn root(&self) -> f64 {
        self.x
    }
    /// OCCT math_FunctionRoot::Derivative() (math_FunctionRoot.hxx L71) —
    /// TheDerivative.
    pub fn derivative(&self) -> f64 {
        self._dfx
    }
    /// OCCT math_FunctionRoot::Value() (math_FunctionRoot.hxx L75) —
    /// TheError (the value of the function at the root).
    pub fn value(&self) -> f64 {
        self._fx
    }
    /// OCCT math_FunctionRoot::NbIterations() (math_FunctionRoot.hxx L80).
    pub fn nb_iterations(&self) -> i32 {
        self.it
    }
}
