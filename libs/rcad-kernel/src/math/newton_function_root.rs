// OCCT math_NewtonFunctionRoot (math_NewtonFunctionRoot.cxx) 1:1 Rust
// translation — Newton-Raphson root finding with bounds clamping and best-X
// tracking.
//
// OCCT ref: math_NewtonFunctionRoot.cxx (whole file); the "full range"
// constructor (F, Guess, EpsX, EpsF, NbIterations) used by
// math_TrigonometricFunctionRoots sets Binf=RealFirst, Bsup=RealLast.

use super::root::FunctionWithDerivative;

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

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT Root().
    pub fn root(&self) -> f64 {
        self.x
    }
}
