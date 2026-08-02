// OCCT TKMath ports used by IntStart_SearchOnBoundaries:
//   math_FunctionSample, math_FunctionRoots (NEWCODE), math_FunctionAllRoots,
//   math_BracketedRoot.
//
// 1:1 Rust translation.  math_FunctionWithDerivative / math_Function are
// represented as Rust traits; exceptions (Standard_NumericError_Raise_if) are
// mapped to panic! with the same control flow.

/// OCCT math_Function — abstract function of a single variable (Value only).
pub trait FunctionValue {
    /// Computes the value of the function <F> for variable <X>.
    fn value(&mut self, x: f64) -> Option<f64>;
}

/// OCCT math_FunctionWithDerivative — function with first derivative.
pub trait FunctionWithDerivative: FunctionValue {
    /// Computes the derivative <D> of the function for variable <X>.
    fn derivative(&mut self, x: f64) -> Option<f64>;
    /// Computes the value <F> and derivative <D> for variable <X>.
    fn values(&mut self, x: f64) -> Option<(f64, f64)>;
    /// Returns the state of the function after the latest call.
    fn get_state_number(&mut self) -> i32 {
        0
    }
}

/// OCCT math_FunctionSample (math_FunctionSample.hxx/.cxx).
/// Default sample with constant parameter step between A and B.
pub struct FunctionSample {
    a: f64,
    b: f64,
    n: i32,
}

impl FunctionSample {
    pub fn new(a: f64, b: f64, n: i32) -> Self {
        FunctionSample { a, b, n }
    }
    /// OCCT Bounds(A, B).
    pub fn bounds(&self) -> (f64, f64) {
        (self.a, self.b)
    }
    /// OCCT NbPoints().
    pub fn nb_points(&self) -> i32 {
        self.n
    }
    /// OCCT GetParameter(Index): A + ((Index-1)/(NbPoints-1))*B.
    /// 1-based Index; panics if Index <= 0 or Index > n.
    pub fn get_parameter(&self, index: i32) -> f64 {
        assert!(index > 0 && index <= self.n, "FunctionSample::GetParameter out of range");
        ((self.n - index) as f64 * self.a + (index - 1) as f64 * self.b) / ((self.n - 1) as f64)
    }
}

/// OCCT math_FunctionRoots (math_FunctionRoots.cxx, NEWCODE path).
/// Finds all real roots of F-K within [A,B] using a sample.
pub struct FunctionRoots {
    done: bool,
    all_null: bool,
    sol: Vec<f64>,
    nb_state_sol: Vec<i32>,
}

impl FunctionRoots {
    /// OCCT math_FunctionRoots(F, A, B, NbSample, EpsX, EpsF, EpsNull, K).
    pub fn new(
        f: &mut dyn FunctionWithDerivative,
        a: f64,
        b: f64,
        nb_sample: i32,
        eps_x: f64,
        eps_f: f64,
        eps_null: f64,
        k: f64,
    ) -> Self {
        let itmax = 100;
        let epseps = 2e-14;

        let mut sol: Vec<f64> = Vec::new();
        let mut nb_state_sol: Vec<i32> = Vec::new();

        let mut done = true;
        let (mut x0, mut xn) = (a, b);
        let mut n = nb_sample;
        if b < a {
            x0 = b;
            xn = a;
        }
        n *= 2;
        if n < 20 {
            n = 20;
        }
        let mut eps_x = eps_x;
        let delta_u = x0.abs() + xn.abs();
        let neps_x = 0.0000000001 * delta_u;
        if eps_x < neps_x {
            eps_x = neps_x;
        }

        let mut x = x0;
        let dx = (xn - x0) / n as f64;
        let mut ptrval: Vec<f64> = Vec::with_capacity((n + 1) as usize);
        let mut nvalid = -1;
        let mut aux = 0.0;
        let mut i = 0;
        while i <= n {
            if x > xn {
                x = xn;
            }
            if let Some(v) = f.value(x) {
                aux = v;
                nvalid += 1;
                ptrval.push(aux - k);
            }
            x += dx;
            i += 1;
        }
        if nvalid < n {
            done = false;
            return FunctionRoots {
                done,
                all_null: false,
                sol,
                nb_state_sol,
            };
        }

        let mut all_null = true;
        for i in 0..=n as usize {
            if ptrval[i] > eps_null || ptrval[i] < -eps_null {
                all_null = false;
            }
        }

        let mut solve = |f: &mut dyn FunctionWithDerivative,
                         k: f64,
                         x1: f64,
                         y1: f64,
                         x2: f64,
                         y2: f64,
                         tol: f64,
                         dx: f64,
                         sol: &mut Vec<f64>,
                         nb_state_sol: &mut Vec<i32>| {
            let tols2 = 0.5 * tol;
            let (mut a_, mut b_, mut c_) = (x1, x2, x2);
            let (mut fa, mut fb, mut fc) = (y1, y2, y2);
            let mut d = 0.0;
            let mut e = 0.0;
            let mut iter = 0;
            while iter < itmax {
                if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
                    c_ = a_;
                    fc = fa;
                    e = b_ - a_;
                    d = e;
                }
                if fc.abs() < fb.abs() {
                    a_ = b_;
                    b_ = c_;
                    c_ = a_;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }
                let tol1 = epseps * b_.abs() + tols2;
                let xm = 0.5 * (c_ - b_);
                if xm.abs() < tol1 || fb == 0.0 {
                    // Try a Newton iteration.
                    let mut xp = b_;
                    let mut itern = 5;
                    let mut ok = true;
                    while ok && itern >= 0 {
                        if let Some((yp, dp)) = f.values(xp) {
                            let mut ok2 = false;
                            if dp > 1e-10 || dp < -1e-10 {
                                xp = xp - (yp - k) / dp;
                            }
                            if xp <= x2 && xp >= x1 {
                                if let Some(y2v) = f.value(xp) {
                                    let y2c = y2v - k;
                                    if y2c.abs() < fb.abs() {
                                        b_ = xp;
                                        fb = y2c;
                                        ok2 = true;
                                    }
                                }
                            }
                            ok = ok2;
                        } else {
                            ok = false;
                        }
                        itern -= 1;
                    }
                    append_root(sol, nb_state_sol, b_, f, dx);
                    return;
                }
                if e.abs() >= tol1 && fa.abs() > fb.abs() {
                    let s = fb / fa;
                    let (mut p, mut q);
                    if a_ == c_ {
                        p = xm * s;
                        p += p;
                        q = 1.0 - s;
                    } else {
                        q = fa / fc;
                        let r = fb / fc;
                        p = s * ((xm + xm) * q * (q - r) - (b_ - a_) * (r - 1.0));
                        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                    }
                    if p > 0.0 {
                        q = -q;
                    }
                    let p_abs = p.abs();
                    let min1 = 3.0 * xm * q - (tol1 * q).abs();
                    let min2 = (e * q).abs();
                    if p_abs + p_abs < if min1 < min2 { min1 } else { min2 } {
                        e = d;
                        d = p / q;
                    } else {
                        d = xm;
                        e = d;
                    }
                } else {
                    d = xm;
                    e = d;
                }
                a_ = b_;
                fa = fb;
                if d.abs() > tol1 {
                    b_ += d;
                } else if xm >= 0.0 {
                    b_ += tol1.abs();
                } else {
                    b_ += -tol1.abs();
                }
                if let Some(v) = f.value(b_) {
                    fb = v - k;
                } else {
                    done = false;
                    return;
                }
                iter += 1;
            }
            // Non-convergence: falls through (OCCT prints debug warning only).
        };

        let append_root_sig = |sol: &mut Vec<f64>,
                               nb_state_sol: &mut Vec<i32>,
                               x: f64,
                               f: &mut dyn FunctionWithDerivative,
                               dx: f64| {
            append_root(sol, nb_state_sol, x, f, dx);
        };

        if !all_null {
            let mut ip1: i32 = 1;
            let mut x = x0;
            let tol = eps_x;
            let mut i = 0;
            while i < n {
                let mut x2 = x + dx;
                if x2 > xn {
                    x2 = xn;
                }
                if ptrval[i as usize] < 0.0 {
                    if ptrval[ip1 as usize] > 0.0 {
                        solve(f, k, x, ptrval[i as usize], x2, ptrval[ip1 as usize], tol, neps_x, &mut sol, &mut nb_state_sol);
                    }
                } else if ptrval[ip1 as usize] < 0.0 {
                    solve(f, k, x, ptrval[i as usize], x2, ptrval[ip1 as usize], tol, neps_x, &mut sol, &mut nb_state_sol);
                }
                i += 1;
                ip1 += 1;
                x += dx;
            }
            // Zeros on sample points.
            let mut i = 0;
            while i <= n {
                if ptrval[i as usize] == 0.0 {
                    let mut x = x0 + i as f64 * dx;
                    if x > xn {
                        x = xn;
                    }
                    let mut u0 = dx * 0.5;
                    let mut u1 = x + u0;
                    u0 += x;
                    if u0 < x0 {
                        u0 = x0;
                    }
                    if u0 > xn {
                        u0 = xn;
                    }
                    if u1 < x0 {
                        u1 = x0;
                    }
                    if u1 > xn {
                        u1 = xn;
                    }
                    let (mut y0, mut y1) = (0.0, 0.0);
                    if let Some(v0) = f.value(u0) {
                        y0 = v0 - k;
                    }
                    if let Some(v1) = f.value(u1) {
                        y1 = v1 - k;
                    }
                    if y0 * y1 < 0.0 {
                        solve(f, k, u0, y0, u1, y1, tol, neps_x, &mut sol, &mut nb_state_sol);
                    } else if y0 != 0.0 || y1 != 0.0 {
                        append_root_sig(&mut sol, &mut nb_state_sol, x, f, neps_x);
                    }
                }
                i += 1;
            }
            // Endpoints.
            if ptrval[0] <= eps_f && ptrval[0] >= -eps_f {
                append_root_sig(&mut sol, &mut nb_state_sol, x0, f, neps_x);
            }
            if ptrval[n as usize] <= eps_f && ptrval[n as usize] >= -eps_f {
                append_root_sig(&mut sol, &mut nb_state_sol, xn, f, neps_x);
            }
            // Extrema re-discretization (positive minima / negative maxima).
            let majdx = 5.0 * dx;
            let mut im1: i32 = 0;
            let mut ip1: i32 = 2;
            let mut i = 1;
            let mut xm = x0 + dx;
            while i < n {
                let mut rediscr = false;
                if xm > xn {
                    xm = xn;
                }
                if ptrval[i as usize] > 0.0 {
                    if (ptrval[im1 as usize] > ptrval[i as usize]) && (ptrval[ip1 as usize] > ptrval[i as usize]) {
                        // Estimate from X_{i-1}.
                        let mut xm1 = xm - dx;
                        if xm1 < x0 {
                            xm1 = x0;
                        }
                        if let Some((mut ym, dym)) = f.values(xm1) {
                            ym -= k;
                            if dym < -1e-10 || dym > 1e-10 {
                                let t = ym / dym;
                                if t < majdx && t > -majdx {
                                    rediscr = true;
                                }
                            }
                        }
                        // Estimate from X_{i+1}.
                        if !rediscr {
                            let mut xp1 = xm + dx;
                            if xp1 > xn {
                                xp1 = xn;
                            }
                            if let Some((mut ym, dym)) = f.values(xp1) {
                                ym -= k;
                                if dym < -1e-10 || dym > 1e-10 {
                                    let t = ym / dym;
                                    if t < majdx && t > -majdx {
                                        rediscr = true;
                                    }
                                }
                            }
                        }
                    }
                } else if ptrval[i as usize] < 0.0 {
                    if (ptrval[im1 as usize] < ptrval[i as usize]) && (ptrval[ip1 as usize] < ptrval[i as usize]) {
                        let mut xm1 = xm - dx;
                        if xm1 < x0 {
                            xm1 = x0;
                        }
                        if let Some((mut ym, dym)) = f.values(xm1) {
                            ym -= k;
                            if dym > 1e-10 || dym < -1e-10 {
                                let t = ym / dym;
                                if t < majdx && t > -majdx {
                                    rediscr = true;
                                }
                            }
                        }
                        if !rediscr {
                            let mut xm2 = xm - dx;
                            if xm2 < x0 {
                                xm2 = x0;
                            }
                            if let Some((mut ym, dym)) = f.values(xm2) {
                                ym -= k;
                                if dym > 1e-10 || dym < -1e-10 {
                                    let t = ym / dym;
                                    if t < majdx && t > -majdx {
                                        rediscr = true;
                                    }
                                }
                            }
                        }
                    }
                }
                if rediscr {
                    let mut x0r = xm - dx;
                    let mut x3r = xm + dx;
                    if x0r < x0 {
                        x0r = x0;
                    }
                    if x3r > xn {
                        x3r = xn;
                    }
                    let mut a_sol_x1 = 0.0;
                    let mut a_sol_x2 = 0.0;
                    let mut a_val1 = 0.0;
                    let mut a_val2 = 0.0;
                    let mut a_der1 = 0.0;
                    let mut a_der2 = 0.0;
                    let mut is_sol1 = false;
                    let mut is_sol2 = false;
                    // Find minimum of |F| between x0r and x3r via derivative zero.
                    {
                        let mut a_der_f = DerivFunction(f);
                        let a_br = BracketedRoot::new(&mut a_der_f, x0r, x3r, eps_x, 100, 1.0e-12);
                        if a_br.is_done() {
                            a_sol_x1 = a_br.root();
                            if let Some(v) = f.value(a_sol_x1) {
                                a_val1 = v.abs();
                                if a_val1 < eps_f {
                                    is_sol1 = true;
                                    a_der1 = a_br.value();
                                }
                            }
                        }
                    }
                    // Golden-section search for the extrema between x0r and x3r.
                    let (mut x1, mut x2);
                    let (mut f0, mut f3);
                    let r = 0.61803399;
                    let c = 1.0 - r;
                    let tol_cr = neps_x * 10.0;
                    f0 = ptrval[im1 as usize];
                    f3 = ptrval[ip1 as usize];
                    let recherche_minimum = f0 > 0.0;
                    if (x3r - xm).abs() > (x0r - xm).abs() {
                        x1 = xm;
                        x2 = xm + c * (x3r - xm);
                    } else {
                        x2 = xm;
                        x1 = xm - c * (xm - x0r);
                    }
                    let (mut f1, mut f2) = (0.0, 0.0);
                    if let Some(v) = f.value(x1) {
                        f1 = v - k;
                    }
                    if let Some(v) = f.value(x2) {
                        f2 = v - k;
                    }
                    let tol_x = 0.001 * neps_x;
                    while (x3r - x0r).abs() > tol_cr * (x1.abs() + x2.abs())
                        && (x1 - x2).abs() > tol_x
                    {
                        if recherche_minimum {
                            if f2 < f1 {
                                x0r = x1;
                                x1 = x2;
                                x2 = r * x1 + c * x3r;
                                f0 = f1;
                                f1 = f2;
                                if let Some(v) = f.value(x2) {
                                    f2 = v - k;
                                }
                            } else {
                                x3r = x2;
                                x2 = x1;
                                x1 = r * x2 + c * x0r;
                                f3 = f2;
                                f2 = f1;
                                if let Some(v) = f.value(x1) {
                                    f1 = v - k;
                                }
                            }
                        } else if f2 > f1 {
                            x0r = x1;
                            x1 = x2;
                            x2 = r * x1 + c * x3r;
                            f0 = f1;
                            f1 = f2;
                            if let Some(v) = f.value(x2) {
                                f2 = v - k;
                            }
                        } else {
                            x3r = x2;
                            x2 = x1;
                            x1 = r * x2 + c * x0r;
                            f3 = f2;
                            f2 = f1;
                            if let Some(v) = f.value(x1) {
                                f1 = v - k;
                            }
                        }
                        if f1 * f0 < 0.0 {
                            solve(f, k, x0r, f0, x1, f1, tol, neps_x, &mut sol, &mut nb_state_sol);
                        }
                        if f2 * f3 < 0.0 {
                            solve(f, k, x2, f2, x3r, f3, tol, neps_x, &mut sol, &mut nb_state_sol);
                        }
                    }
                    if (recherche_minimum && f1 < f2) || (!recherche_minimum && f1 > f2) {
                        if f1.abs() < eps_f {
                            is_sol2 = true;
                            a_sol_x2 = x1;
                            a_val2 = f1.abs();
                        }
                    } else if f2.abs() < eps_f {
                        is_sol2 = true;
                        a_sol_x2 = x2;
                        a_val2 = f2.abs();
                    }
                    // Choose the best solution between aSolX1, aSolX2.
                    if is_sol1 && is_sol2 {
                        if a_val2 - a_val1 > eps_f {
                            append_root_sig(&mut sol, &mut nb_state_sol, a_sol_x1, f, neps_x);
                        } else if a_val1 - a_val2 > eps_f {
                            append_root_sig(&mut sol, &mut nb_state_sol, a_sol_x2, f, neps_x);
                        } else {
                            a_der1 = a_der1.abs();
                            if let Some(v) = f.derivative(a_sol_x2) {
                                a_der2 = v.abs();
                            }
                            if a_der1 < a_der2 {
                                append_root_sig(&mut sol, &mut nb_state_sol, a_sol_x1, f, neps_x);
                            } else {
                                append_root_sig(&mut sol, &mut nb_state_sol, a_sol_x2, f, neps_x);
                            }
                        }
                    } else if is_sol1 {
                        append_root_sig(&mut sol, &mut nb_state_sol, a_sol_x1, f, neps_x);
                    } else if is_sol2 {
                        append_root_sig(&mut sol, &mut nb_state_sol, a_sol_x2, f, neps_x);
                    }
                }
                i += 1;
                im1 += 1;
                ip1 += 1;
                xm += dx;
            }
        }

        FunctionRoots {
            done,
            all_null,
            sol,
            nb_state_sol,
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT IsAllNull().
    pub fn is_all_null(&self) -> bool {
        self.all_null
    }
    /// OCCT NbSolutions().
    pub fn nb_solutions(&self) -> usize {
        self.sol.len()
    }
    /// OCCT Value(Nieme) — 1-based.
    pub fn value(&self, n: usize) -> f64 {
        self.sol[n - 1]
    }
    /// OCCT StateNumber(Nieme) — 1-based.
    pub fn state_number(&self, n: usize) -> i32 {
        self.nb_state_sol[n - 1]
    }
}

/// OCCT static AppendRoot (math_FunctionRoots.cxx L51-108).
fn append_root(
    sol: &mut Vec<f64>,
    nb_state_sol: &mut Vec<i32>,
    x: f64,
    f: &mut dyn FunctionWithDerivative,
    dx: f64,
) {
    let n = sol.len();
    if n == 0 {
        sol.push(x);
        if let Some(t) = f.value(x) {
            let _ = t;
        }
        nb_state_sol.push(f.get_state_number());
    } else {
        let mut i = 1usize;
        let mut pl = n + 1;
        while i <= n {
            let t = sol[i - 1];
            if t >= x {
                pl = i;
                i = n;
            }
            if (x - t).abs() <= dx {
                pl = 0;
                i = n;
            }
            i += 1;
        }
        if pl > n {
            sol.push(x);
            if let Some(_t) = f.value(x) {}
            nb_state_sol.push(f.get_state_number());
        } else if pl > 0 {
            sol.insert(pl - 1, x);
            if let Some(_t) = f.value(x) {}
            nb_state_sol.insert(pl - 1, f.get_state_number());
        }
    }
}

/// OCCT class DerivFunction (math_FunctionRoots.cxx L38-49) — the derivative of
/// a FunctionWithDerivative viewed as a plain FunctionValue.
struct DerivFunction<'a>(&'a mut dyn FunctionWithDerivative);

impl FunctionValue for DerivFunction<'_> {
    fn value(&mut self, x: f64) -> Option<f64> {
        self.0.derivative(x)
    }
}

/// OCCT math_BracketedRoot (math_BracketedRoot.cxx) — Brent's method root
/// finding on a bracketed interval (Numerical Recipes p.269).
pub struct BracketedRoot {
    the_root: f64,
    the_error: f64,
    nb_iter: i32,
    done: bool,
}

impl BracketedRoot {
    /// OCCT math_BracketedRoot(F, Bound1, Bound2, Tolerance, NbIterations, ZEPS).
    pub fn new(
        f: &mut dyn FunctionValue,
        bound1: f64,
        bound2: f64,
        tolerance: f64,
        nb_iterations: i32,
        zeps: f64,
    ) -> Self {
        let mut a = bound1;
        let mut the_root = bound2;
        let mut fa = f.value(a).unwrap_or(f64::MAX);
        let mut the_error = f.value(the_root).unwrap_or(f64::MAX);
        if fa * the_error > 0.0 {
            return BracketedRoot {
                the_root,
                the_error,
                nb_iter: 0,
                done: false,
            };
        }
        let mut fc = the_error;
        let mut c = 0.0;
        let mut d = 0.0;
        let mut e = 0.0;
        let mut nb_iter = 0;
        while nb_iter < nb_iterations {
            nb_iter += 1;
            if the_error * fc > 0.0 {
                c = a;
                fc = fa;
                d = the_root - a;
                e = d;
            }
            if fc.abs() < fa.abs() {
                a = the_root;
                the_root = c;
                c = a;
                fa = the_error;
                the_error = fc;
                fc = fa;
            }
            let tol1 = 2.0 * zeps * the_root.abs() + 0.5 * tolerance;
            let xm = 0.5 * (c - the_root);
            if xm.abs() <= tol1 || the_error == 0.0 {
                return BracketedRoot {
                    the_root,
                    the_error,
                    nb_iter,
                    done: true,
                };
            }
            if e.abs() >= tol1 && fa.abs() > the_error.abs() {
                let s = the_error / fa;
                let (mut p, mut q);
                if a == c {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    let r = the_error / fc;
                    p = s * (2.0 * xm * q * (q - r) - (the_root - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if p > 0.0 {
                    q = -q;
                }
                let p = p.abs();
                let min1 = 3.0 * xm * q - (tol1 * q).abs();
                let min2 = (e * q).abs();
                if 2.0 * p < if min1 < min2 { min1 } else { min2 } {
                    e = d;
                    d = p / q;
                } else {
                    d = xm;
                    e = d;
                }
            } else {
                d = xm;
                e = d;
            }
            a = the_root;
            fa = the_error;
            if d.abs() > tol1 {
                the_root += d;
            } else if xm > 0.0 {
                the_root += tol1.abs();
            } else {
                the_root += -tol1.abs();
            }
            the_error = f.value(the_root).unwrap_or(f64::MAX);
        }
        BracketedRoot {
            the_root,
            the_error,
            nb_iter,
            done: false,
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT Root().
    pub fn root(&self) -> f64 {
        self.the_root
    }
    /// OCCT Value() — the function value at the root.
    pub fn value(&self) -> f64 {
        self.the_error
    }
    /// OCCT NbIterations().
    pub fn nb_iterations(&self) -> i32 {
        self.nb_iter
    }
}

/// OCCT math_FunctionAllRoots (math_FunctionAllRoots.cxx).
/// Finds all null intervals and all isolated roots of F on a sample.
pub struct FunctionAllRoots {
    done: bool,
    pdeb: Vec<f64>,
    pfin: Vec<f64>,
    piso: Vec<f64>,
    ideb: Vec<i32>,
    ifin: Vec<i32>,
    iiso: Vec<i32>,
}

impl FunctionAllRoots {
    /// OCCT math_FunctionAllRoots(F, S, EpsX, EpsF, EpsNul).
    pub fn new(
        f: &mut dyn FunctionWithDerivative,
        s: &FunctionSample,
        eps_x: f64,
        eps_f: f64,
        eps_nul: f64,
    ) -> Self {
        let mut pdeb: Vec<f64> = Vec::new();
        let mut pfin: Vec<f64> = Vec::new();
        let mut piso: Vec<f64> = Vec::new();
        let mut ideb: Vec<i32> = Vec::new();
        let mut ifin: Vec<i32> = Vec::new();
        let mut iiso: Vec<i32> = Vec::new();

        let nbp = s.nb_points();
        let mut pnul;
        let mut inter_nul = false;
        let mut nuld = false;
        let mut nulf = false;
        let mut deb_nul = 0.0;
        let mut fin_nul = 0.0;
        let mut indd = 0i32;
        let mut indf = 0i32;

        let mut val = f.value(s.get_parameter(1)).unwrap_or(f64::MAX);
        pnul = val.abs() <= eps_nul;
        let mut valsav = if pnul { 0.0 } else { val };

        let mut i = 2i32;
        let mut fini = i > nbp;

        while !fini {
            val = f.value(s.get_parameter(i)).unwrap_or(f64::MAX);
            let nul = val.abs() <= eps_nul;
            if !nul {
                valsav = val;
            }
            if inter_nul && !nul {
                inter_nul = false;
                pdeb.push(deb_nul);
                ideb.push(indd);
                let mut cst = if val > 0.0 { eps_nul } else { -eps_nul };
                let mut res1 = FunctionRoots::new(
                    f,
                    s.get_parameter(i - 1),
                    s.get_parameter(i),
                    10,
                    eps_x,
                    eps_f,
                    0.0,
                    cst,
                );
                assert!(
                    res1.is_done() && !res1.is_all_null() && res1.nb_solutions() != 0,
                    "math_FunctionAllRoots: Res1 failed"
                );
                fin_nul = res1.value(1);
                indf = res1.state_number(1);

                cst = -cst;
                let mut res2 = FunctionRoots::new(
                    f,
                    s.get_parameter(i - 1),
                    s.get_parameter(i),
                    10,
                    eps_x,
                    eps_f,
                    0.0,
                    cst,
                );
                assert!(res2.is_done() && !res2.is_all_null(), "math_FunctionAllRoots: Res2 failed");
                if res2.nb_solutions() != 0 {
                    if res2.value(1) < fin_nul {
                        fin_nul = res2.value(1);
                        indf = res2.state_number(1);
                    }
                }
                pfin.push(fin_nul);
                ifin.push(indf);
            } else if !inter_nul && pnul && nul {
                inter_nul = true;
                if i == 2 {
                    deb_nul = s.get_parameter(1);
                    f.value(deb_nul);
                    indd = f.get_state_number();
                    nuld = true;
                } else {
                    let mut cst = if valsav > 0.0 { eps_nul } else { -eps_nul };
                    let mut res1 = FunctionRoots::new(
                        f,
                        s.get_parameter(i - 2),
                        s.get_parameter(i - 1),
                        10,
                        eps_x,
                        eps_f,
                        0.0,
                        cst,
                    );
                    assert!(
                        res1.is_done() && !res1.is_all_null() && res1.nb_solutions() != 0,
                        "math_FunctionAllRoots: Res1b failed"
                    );
                    deb_nul = res1.value(res1.nb_solutions());
                    indd = res1.state_number(res1.nb_solutions());

                    cst = -cst;
                    let mut res3 = FunctionRoots::new(
                        f,
                        s.get_parameter(i - 2),
                        s.get_parameter(i - 1),
                        10,
                        eps_x,
                        eps_f,
                        0.0,
                        cst,
                    );
                    assert!(res3.is_done() && !res3.is_all_null(), "math_FunctionAllRoots: Res3 failed");
                    if res3.nb_solutions() != 0 {
                        if res3.value(res3.nb_solutions()) > deb_nul {
                            deb_nul = res3.value(res3.nb_solutions());
                            indd = res3.state_number(res3.nb_solutions());
                        }
                    }
                }
            }
            i += 1;
            pnul = nul;
            fini = i > nbp;
        }

        if inter_nul {
            // Add the interval ending at the last point.
            pdeb.push(deb_nul);
            ideb.push(indd);
            fin_nul = s.get_parameter(nbp);
            f.value(fin_nul);
            indf = f.get_state_number();
            pfin.push(fin_nul);
            ifin.push(indf);
            nulf = true;
        }

        if pdeb.is_empty() {
            // No null interval.
            let mut res = FunctionRoots::new(
                f,
                s.get_parameter(1),
                s.get_parameter(nbp),
                nbp,
                eps_x,
                eps_f,
                0.0,
                0.0,
            );
            assert!(res.is_done() && !res.is_all_null(), "math_FunctionAllRoots: Res failed");
            for j in 1..=res.nb_solutions() {
                piso.push(res.value(j));
                iiso.push(res.state_number(j));
            }
        } else {
            let nbp_min = 3;
            if !nuld {
                // Roots between the first sample point and the start of the 1st null interval.
                let nbrpt = (((pdeb[0] - s.get_parameter(1)).abs()
                    / (s.get_parameter(nbp) - s.get_parameter(1)))
                    * nbp as f64).trunc() as i32;
                let mut res = FunctionRoots::new(
                    f,
                    s.get_parameter(1),
                    pdeb[0],
                    if nbrpt > nbp_min { nbrpt } else { nbp_min },
                    eps_x,
                    eps_f,
                    0.0,
                    0.0,
                );
                assert!(res.is_done() && !res.is_all_null(), "math_FunctionAllRoots: Res (pre) failed");
                for j in 1..=res.nb_solutions() {
                    piso.push(res.value(j));
                    iiso.push(res.state_number(j));
                }
            }
            let npdeb = pdeb.len();
            for k in 2..=npdeb {
                let nbrpt = (((pdeb[k - 1] - pfin[k - 2]).abs()
                    / (s.get_parameter(nbp) - s.get_parameter(1)))
                    * nbp as f64).trunc() as i32;
                let mut res = FunctionRoots::new(
                    f,
                    pfin[k - 2],
                    pdeb[k - 1],
                    if nbrpt > nbp_min { nbrpt } else { nbp_min },
                    eps_x,
                    eps_f,
                    0.0,
                    0.0,
                );
                assert!(res.is_done() && !res.is_all_null(), "math_FunctionAllRoots: Res (mid) failed");
                for j in 1..=res.nb_solutions() {
                    piso.push(res.value(j));
                    iiso.push(res.state_number(j));
                }
            }
            if !nulf {
                // Roots between the end of the last null interval and the last sample point.
                let nbrpt = (((s.get_parameter(nbp) - pfin[pdeb.len() - 1]).abs()
                    / (s.get_parameter(nbp) - s.get_parameter(1)))
                    * nbp as f64).trunc() as i32;
                let mut res = FunctionRoots::new(
                    f,
                    pfin[pdeb.len() - 1],
                    s.get_parameter(nbp),
                    if nbrpt > nbp_min { nbrpt } else { nbp_min },
                    eps_x,
                    eps_f,
                    0.0,
                    0.0,
                );
                assert!(res.is_done() && !res.is_all_null(), "math_FunctionAllRoots: Res (post) failed");
                for j in 1..=res.nb_solutions() {
                    piso.push(res.value(j));
                    iiso.push(res.state_number(j));
                }
            }
        }

        FunctionAllRoots {
            done: true,
            pdeb,
            pfin,
            piso,
            ideb,
            ifin,
            iiso,
        }
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT NbIntervals().
    pub fn nb_intervals(&self) -> usize {
        self.pdeb.len()
    }
    /// OCCT GetInterval(Index, A, B) — 1-based.
    pub fn get_interval(&self, index: usize) -> (f64, f64) {
        (self.pdeb[index - 1], self.pfin[index - 1])
    }
    /// OCCT GetIntervalState(Index, IFirst, ILast) — 1-based.
    pub fn get_interval_state(&self, index: usize) -> (i32, i32) {
        (self.ideb[index - 1], self.ifin[index - 1])
    }
    /// OCCT NbPoints().
    pub fn nb_points(&self) -> usize {
        self.piso.len()
    }
    /// OCCT GetPoint(Index) — 1-based.
    pub fn get_point(&self, index: usize) -> f64 {
        self.piso[index - 1]
    }
    /// OCCT GetPointState(Index) — 1-based.
    pub fn get_point_state(&self, index: usize) -> i32 {
        self.iiso[index - 1]
    }
}
