// OCCT math_DirectPolynomialRoots (math_DirectPolynomialRoots.cxx) 1:1 Rust
// translation — robust real-root finding for polynomials up to degree 4
// (Ferrari / Cardano / stable quadratic, with Newton-Raphson refinement).
//
// OCCT refs:
//   - math_DirectPolynomialRoots.cxx (whole file)
//   - Standard_Real.hxx L132-246 (RealSmall/RealEpsilon/Epsilon)

/// OCCT Standard_Real.hxx Epsilon(Value) — the gap to the next representable
/// double in the direction of the value's sign.
pub fn epsilon(value: f64) -> f64 {
    if value >= 0.0 {
        value.next_up() - value
    } else {
        value - value.next_down()
    }
}

/// OCCT math_DirectPolynomialRoots.hxx — the direct (analytic) polynomial
/// root solver.  `value(i)` is 1-based like OCCT's Value(Index).
pub struct DirectPolynomialRoots {
    done: bool,
    infinite: bool,
    nb_sol: usize,
    roots: [f64; 4],
}

impl DirectPolynomialRoots {
    /// OCCT `math_DirectPolynomialRoots(A, B, C, D, E)` (degree 4).
    pub fn new_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> Self {
        let mut s = DirectPolynomialRoots {
            done: true,
            infinite: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        s.solve_quartic(a, b, c, d, e);
        s
    }

    /// OCCT `math_DirectPolynomialRoots(A, B, C, D)` (degree 3).
    pub fn new_cubic(a: f64, b: f64, c: f64, d: f64) -> Self {
        let mut s = DirectPolynomialRoots {
            done: true,
            infinite: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        s.solve_cubic(a, b, c, d);
        s
    }

    /// OCCT `math_DirectPolynomialRoots(A, B, C)` (degree 2).
    pub fn new_quadratic(a: f64, b: f64, c: f64) -> Self {
        let mut s = DirectPolynomialRoots {
            done: true,
            infinite: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        s.solve_quadratic(a, b, c);
        s
    }

    /// OCCT `math_DirectPolynomialRoots(A, B)` (degree 1).
    pub fn new_linear(a: f64, b: f64) -> Self {
        let mut s = DirectPolynomialRoots {
            done: true,
            infinite: false,
            nb_sol: 0,
            roots: [0.0; 4],
        };
        s.solve_linear(a, b);
        s
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT InfiniteRoots().
    pub fn infinite_roots(&self) -> bool {
        self.infinite
    }
    /// OCCT NbSolutions().
    pub fn nb_solutions(&self) -> usize {
        self.nb_sol
    }
    /// OCCT Value(Index) — 1-based.
    pub fn value(&self, index: usize) -> f64 {
        self.roots[index - 1]
    }

    // =====================================================================
    // OCCT Solve(A, B, C, D, E) — quartic
    // =====================================================================
    fn solve_quartic(&mut self, a: f64, b: f64, c: f64, d: f64, e: f64) {
        const ZERO_THRESHOLD: f64 = 1.0e-30;
        const MACHINE_EPSILON: f64 = f64::EPSILON;
        const OVERFLOW_LIMIT: f64 = 1.0e+80;
        let _ = (ZERO_THRESHOLD, OVERFLOW_LIMIT);

        // OCCT L562-567: degree reduction.
        if should_reduce_degree_quartic(a, b, c, d, e) {
            self.solve_cubic(b, c, d, e);
            return;
        }

        // OCCT L569-573: normalize coefficients.
        let an = b / a;
        let bn = c / a;
        let cn = d / a;
        let dn = e / a;

        // OCCT L575-577: scale to avoid overflow/underflow.
        let mut scaled = ScaledCoefficients::default();
        scaled.scale_quartic(an, bn, cn, dn, dn);

        // OCCT L580-587: Ferrari's resolvent cubic.
        let mut a_success = false;
        let y0 = solve_ferrari_resolvent(scaled.a, scaled.b, scaled.c, scaled.d, &mut a_success);
        if std::env::var("RCAD_DPR_DEBUG").is_ok() {
            eprintln!("[DPR] quartic a={} b={} c={} d={} e={} scaled=({},{},{},{}) y0={} success={}",
                a, b, c, d, e, scaled.a, scaled.b, scaled.c, scaled.d, y0, a_success);
        }
        if !a_success {
            self.done = false;
            return;
        }

        // OCCT L589-591: factor into two quadratics.
        let factors = factor_quartic_via_ferrari(scaled.a, scaled.b, scaled.c, scaled.d, y0);

        // OCCT L593-599: solve first quadratic.
        let mut q1 = DirectPolynomialRoots::new_quadratic(1.0, factors.p1, factors.q1);
        if std::env::var("RCAD_DPR_DEBUG").is_ok() {
            eprintln!("[DPR] ferrari factors p1={} q1={} p2={} q2={} | q1 done={}", factors.p1, factors.q1, factors.p2, factors.q2, q1.is_done());
        }
        if !q1.is_done() {
            self.done = false;
            return;
        }

        // OCCT L601-607: solve second quadratic.
        let mut q2 = DirectPolynomialRoots::new_quadratic(1.0, factors.p2, factors.q2);
        if !q2.is_done() {
            self.done = false;
            return;
        }

        // OCCT L609-621: collect all roots.
        self.nb_sol = q1.nb_sol + q2.nb_sol;
        let mut a_index = 0;
        for i in 0..q1.nb_sol {
            self.roots[a_index] = q1.roots[i];
            a_index += 1;
        }
        for i in 0..q2.nb_sol {
            self.roots[a_index] = q2.roots[i];
            a_index += 1;
        }

        // OCCT L623-624: inverse scaling + Newton refinement.
        let coeffs = [a, b, c, d, e];
        for i in 0..self.nb_sol {
            self.roots[i] *= scaled.scale_factor;
            self.roots[i] = refine_polynomial_root(self.roots[i], &coeffs);
        }
    }

    // =====================================================================
    // OCCT Solve(A, B, C, D) — cubic
    // =====================================================================
    fn solve_cubic(&mut self, a: f64, b: f64, c: f64, d: f64) {
        const ZERO_THRESHOLD: f64 = 1.0e-30;
        const MACHINE_EPSILON: f64 = f64::EPSILON;
        const OVERFLOW_LIMIT: f64 = 1.0e+80;

        // OCCT L634-639: degree reduction.
        if a.abs() <= ZERO_THRESHOLD {
            self.solve_quadratic(b, c, d);
            return;
        }

        // OCCT L641-644: normalize coefficients.
        let a_beta = b / a;
        let a_gamma = c / a;
        let a_del = d / a;

        // OCCT L646-648: scale.
        let mut a_scaled = ScaledCoefficients::default();
        a_scaled.scale_cubic(a_beta, a_gamma, a_del, a_del);

        // OCCT L650-658: depressed cubic t^3 + P t + Q = 0.
        let a_p1 = a_scaled.b;
        let a_p2 = -(a_scaled.a * a_scaled.a) / 3.0;
        let mut a_p = a_p1 + a_p2;
        let a_ep = 5.0 * MACHINE_EPSILON * (a_p1.abs() + a_p2.abs());
        if a_p.abs() <= a_ep {
            a_p = 0.0;
        }

        let a_q1 = a_scaled.c;
        let a_q2 = -a_scaled.a * a_scaled.b / 3.0;
        let a_q3 = 2.0 * (a_scaled.a * a_scaled.a * a_scaled.a) / 27.0;
        let mut a_q = a_q1 + a_q2 + a_q3;
        let a_eq = 10.0 * MACHINE_EPSILON * (a_q1.abs() + a_q2.abs() + a_q3.abs());
        if a_q.abs() <= a_eq {
            a_q = 0.0;
        }

        // OCCT L670-675: overflow check.
        if a_p.abs() > OVERFLOW_LIMIT {
            self.done = false;
            return;
        }

        // OCCT L677-686: discriminant.
        let a_a1 = (a_p * a_p * a_p) / 27.0;
        let a_a2 = (a_q * a_q) / 4.0;
        let mut a_discr = a_a1 + a_a2;
        if a_p < 0.0 {
            a_discr = compute_special_discriminant(a_scaled.a, a_scaled.b, a_scaled.c, a_a1);
        }

        // OCCT L689-705: solve by discriminant.
        if a_discr < 0.0 {
            self.nb_sol = 3;
            solve_cubic_three_real_roots(
                a_scaled.a, a_scaled.b, a_scaled.c, a_p, a_q, a_discr, &mut self.roots,
            );
        } else if a_discr > 0.0 {
            self.nb_sol = 1;
            solve_cubic_one_real_root(a_scaled.a, a_scaled.c, a_p, a_q, a_discr, &mut self.roots);
        } else {
            solve_cubic_multiple_roots(
                a_scaled.a, a_scaled.b, a_scaled.c, a_p, a_q, &mut self.roots, &mut self.nb_sol,
            );
        }

        // OCCT L707-708: inverse scaling + refinement.
        let coeffs = [a, b, c, d];
        for i in 0..self.nb_sol {
            self.roots[i] *= a_scaled.scale_factor;
            self.roots[i] = refine_polynomial_root(self.roots[i], &coeffs);
        }
    }

    // =====================================================================
    // OCCT Solve(A, B, C) — quadratic
    // =====================================================================
    fn solve_quadratic(&mut self, a: f64, b: f64, c: f64) {
        const ZERO_THRESHOLD: f64 = 1.0e-30;
        const MACHINE_EPSILON: f64 = f64::EPSILON;

        // OCCT L716-720: degree reduction.
        if a.abs() <= ZERO_THRESHOLD {
            self.solve_linear(b, c);
            return;
        }

        // OCCT L722-724: normalize x^2 + P x + Q = 0.
        let p = b / a;
        let q = c / a;

        // OCCT L726-733: discriminant with error bounds.
        let a_eps_d = 3.0 * MACHINE_EPSILON * (p * p + (4.0 * q).abs());
        let mut a_discrim = p * p - 4.0 * q;
        if a_discrim.abs() <= a_eps_d {
            a_discrim = 0.0;
        }

        if a_discrim < 0.0 {
            // OCCT L735-739: no real roots.
            self.nb_sol = 0;
        } else if a_discrim == 0.0 {
            // OCCT L740-747: double root.
            self.nb_sol = 2;
            self.roots[0] = -0.5 * p;
            self.roots[0] = refine_polynomial_root(self.roots[0], &[1.0, p, q]);
            self.roots[1] = self.roots[0];
        } else {
            // OCCT L748-763: two distinct real roots, numerically stable.
            self.nb_sol = 2;
            if p > 0.0 {
                self.roots[0] = -(p + a_discrim.sqrt()) / 2.0;
            } else {
                self.roots[0] = -(p - a_discrim.sqrt()) / 2.0;
            }
            self.roots[0] = refine_polynomial_root(self.roots[0], &[1.0, p, q]);
            self.roots[1] = q / self.roots[0];
            self.roots[1] = refine_polynomial_root(self.roots[1], &[1.0, p, q]);
        }
    }

    // =====================================================================
    // OCCT Solve(A, B) — linear
    // =====================================================================
    fn solve_linear(&mut self, a: f64, b: f64) {
        const ZERO_THRESHOLD: f64 = 1.0e-30;

        if a.abs() <= ZERO_THRESHOLD {
            if b.abs() <= ZERO_THRESHOLD {
                // OCCT L774-776: 0 = 0: infinite solutions.
                self.infinite = true;
                return;
            }
            // OCCT L777-780: no solution.
            self.nb_sol = 0;
            return;
        }

        // OCCT L782-786: unique solution.
        self.nb_sol = 1;
        self.roots[0] = -b / a;
    }
}

/// OCCT anonymous-namespace ScaledCoefficients (L161-193).
#[derive(Default)]
struct ScaledCoefficients {
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    scale_factor: f64,
}

impl ScaledCoefficients {
    fn scale_quartic(&mut self, a: f64, b: f64, c: f64, d: f64, e: f64) {
        const FLOATING_RADIX: f64 = 2.0;
        let a_exp = compute_base_exponent(e) / 4;
        self.scale_factor = FLOATING_RADIX.powi(a_exp);
        let a_scale_factor2 = self.scale_factor * self.scale_factor;
        self.a = a / self.scale_factor;
        self.b = b / a_scale_factor2;
        self.c = c / (a_scale_factor2 * self.scale_factor);
        self.d = d / (a_scale_factor2 * a_scale_factor2);
        // E = theE / (aScaleFactor2 * aScaleFactor2) — unused downstream.
    }

    fn scale_cubic(&mut self, a: f64, b: f64, c: f64, d: f64) {
        const FLOATING_RADIX: f64 = 2.0;
        let a_exp = compute_base_exponent(d) / 3;
        self.scale_factor = FLOATING_RADIX.powi(a_exp);
        let a_scale_factor2 = self.scale_factor * self.scale_factor;
        self.a = a / self.scale_factor;
        self.b = b / a_scale_factor2;
        self.c = c / (a_scale_factor2 * self.scale_factor);
        // D = theD — not used for cubic.
    }
}

/// OCCT ComputeBaseExponent (L147-158).
fn compute_base_exponent(value: f64) -> i32 {
    const INV_LOG_RADIX: f64 = 1.0 / std::f64::consts::LN_2;
    if value > 1.0 {
        (value.ln() * INV_LOG_RADIX) as i32
    } else if value < -1.0 {
        (-(-value).ln() * INV_LOG_RADIX) as i32
    } else {
        0
    }
}

/// OCCT EvaluatePolynomial (L64-72) — Horner.
fn evaluate_polynomial(n: usize, poly: &[f64], x: f64) -> f64 {
    let mut result = poly[0];
    for i in 1..n {
        result = result * x + poly[i];
    }
    result
}

/// OCCT EvaluatePolynomialWithDerivative (L74-89).
fn evaluate_polynomial_with_derivative(n: usize, poly: &[f64], x: f64) -> (f64, f64) {
    let mut value = poly[0] * x + poly[1];
    let mut derivative = poly[0];
    for i in 2..n {
        derivative = derivative * x + value;
        value = value * x + poly[i];
    }
    (value, derivative)
}

/// OCCT RefineRoot (L92-120) — Newton-Raphson refinement.
fn refine_root(n: usize, poly: &[f64], initial_guess: f64) -> f64 {
    const ZERO_THRESHOLD: f64 = 1.0e-30;
    const MACHINE_EPSILON: f64 = f64::EPSILON;
    const MAX_NEWTON_ITERATIONS: i32 = 10;

    let mut a_value = 0.0;
    let mut a_derivative = 0.0;
    let mut a_solution = initial_guess;
    let a_initial_value = evaluate_polynomial(n, poly, initial_guess);

    let mut iter = 1;
    while iter < MAX_NEWTON_ITERATIONS {
        let (v, d) = evaluate_polynomial_with_derivative(n, poly, a_solution);
        a_value = v;
        a_derivative = d;

        if a_derivative.abs() <= ZERO_THRESHOLD {
            break;
        }

        let a_delta = -a_value / a_derivative;
        if a_delta.abs() <= MACHINE_EPSILON * a_solution.abs() {
            break;
        }

        a_solution += a_delta;
        iter += 1;
    }

    // OCCT L119: return the improved solution only if it is better.
    if a_value.abs() <= a_initial_value.abs() {
        a_solution
    } else {
        initial_guess
    }
}

/// OCCT RefinePolynomialRoot (L125-130) — variadic template; Rust passes the
/// coefficient slice.
fn refine_polynomial_root(initial_guess: f64, coeffs: &[f64]) -> f64 {
    refine_root(coeffs.len(), coeffs, initial_guess)
}

/// OCCT ComputeSpecialDiscriminant (L196-223).
fn compute_special_discriminant(beta: f64, gamma: f64, del: f64, a1: f64) -> f64 {
    const MACHINE_EPSILON: f64 = f64::EPSILON;
    let a_sigma = beta * gamma / 3.0 - 2.0 * beta * beta * beta / 27.0;
    let a_psi = gamma * gamma * (4.0 * gamma - beta * beta) / 27.0;

    let a_d1 = if a_sigma >= 0.0 {
        a_sigma + 2.0 * (-a1).sqrt()
    } else {
        a_sigma - 2.0 * (-a1).sqrt()
    };

    let a_d2 = a_psi / a_d1;

    if (del - a_d1).abs() >= 18.0 * MACHINE_EPSILON * (del.abs() + a_d1.abs())
        && (del - a_d2).abs() >= 24.0 * MACHINE_EPSILON * (del.abs() + a_d2.abs())
    {
        return (del - a_d1) * (del - a_d2) / 4.0;
    }
    0.0
}

/// OCCT SolveCubicThreeRealRoots (L226-268).
fn solve_cubic_three_real_roots(
    beta: f64, gamma: f64, del: f64, p: f64, q: f64, discr: f64, roots: &mut [f64; 4],
) {
    if beta == 0.0 && q == 0.0 {
        // x^3 + P x = 0.
        roots[0] = (-p).sqrt();
        roots[1] = -roots[0];
        roots[2] = 0.0;
    } else {
        let a_sb = if beta >= 0.0 { 1.0 } else { -1.0 };
        let a_omega = (0.5 * q / (-discr).sqrt()).atan();
        let a_sp3 = (-p / 3.0).sqrt();
        let a_y1 = -2.0 * a_sb * a_sp3 * (std::f64::consts::PI / 6.0 - a_sb * a_omega / 3.0).cos();

        roots[0] = -beta / 3.0 + a_y1;

        if beta * q <= 0.0 {
            roots[1] = -beta / 3.0 + 2.0 * a_sp3 * (a_omega / 3.0).sin();
        } else {
            let a_dbg = del - beta * gamma;
            let a_sdbg = if a_dbg >= 0.0 { 1.0 } else { -1.0 };
            let a_den1 = 8.0 * beta * beta / 9.0 - 4.0 * beta * a_y1 / 3.0 - 2.0 * q / a_y1;
            let a_den2 = 2.0 * a_y1 * a_y1 - q / a_y1;
            roots[1] = a_dbg / a_den1 + a_sdbg * (-27.0 * discr).sqrt() / a_den2;
        }

        roots[2] = -del / (roots[0] * roots[1]);
    }
}

/// OCCT SolveCubicOneRealRoot (L271-306).
fn solve_cubic_one_real_root(
    beta: f64, del: f64, p: f64, q: f64, discr: f64, roots: &mut [f64; 4],
) {
    const REAL_SMALL: f64 = f64::MIN_POSITIVE;

    // aU = sqrt(discr) + |q/2| is always >= 0, so the OCCT sign ternary takes
    // the positive branch (std::pow(aU, 1/3) = cube root).
    let a_u = (discr.sqrt() + (q / 2.0).abs()).cbrt();

    let a_h = if p >= 0.0 {
        a_u * a_u + p / 3.0 + (p / a_u) * (p / a_u) / 9.0
    } else {
        a_u * q.abs() / (a_u * a_u - p / 3.0)
    };

    if beta * q >= 0.0 {
        if a_h.abs() <= REAL_SMALL && q.abs() <= REAL_SMALL {
            roots[0] = -beta / 3.0 - a_u + p / (3.0 * a_u);
        } else {
            roots[0] = -beta / 3.0 - q / a_h;
        }
    } else {
        roots[0] = -del / (beta * beta / 9.0 + a_h - beta * q / (3.0 * a_h));
    }
}

/// OCCT SolveCubicMultipleRoots (L309-341).
fn solve_cubic_multiple_roots(
    beta: f64, gamma: f64, del: f64, p: f64, q: f64, roots: &mut [f64; 4], nb_roots: &mut usize,
) {
    *nb_roots = 3;
    let a_sq = if q >= 0.0 { 1.0 } else { -1.0 };
    let a_sp3 = (-p / 3.0).sqrt();

    if beta * q <= 0.0 {
        roots[0] = -beta / 3.0 + a_sq * a_sp3;
        roots[1] = roots[0];
        if beta * q == 0.0 {
            roots[2] = -beta / 3.0 - 2.0 * a_sq * a_sp3;
        } else {
            roots[2] = -del / (roots[0] * roots[1]);
        }
    } else {
        roots[0] = -gamma / (beta + 3.0 * a_sq * a_sp3);
        roots[1] = roots[0];
        roots[2] = -beta / 3.0 - 2.0 * a_sq * a_sp3;
    }
}

/// OCCT ShouldReduceDegreeQuartic (L344-393).
fn should_reduce_degree_quartic(a: f64, b: f64, c: f64, d: f64, e: f64) -> bool {
    const ZERO_THRESHOLD: f64 = 1.0e-30;
    if a.abs() <= ZERO_THRESHOLD {
        return true;
    }

    let mut a_max_coeff = ZERO_THRESHOLD;
    a_max_coeff = a_max_coeff.max(b.abs());
    a_max_coeff = a_max_coeff.max(c.abs());
    a_max_coeff = a_max_coeff.max(d.abs());
    a_max_coeff = a_max_coeff.max(e.abs());

    if a_max_coeff > ZERO_THRESHOLD {
        a_max_coeff = epsilon(100.0 * a_max_coeff);
    }

    if a.abs() <= a_max_coeff {
        let a_max_coeff1000 = 1000.0 * a_max_coeff;
        let mut a_with_a = false;
        if b.abs() > ZERO_THRESHOLD && b.abs() <= a_max_coeff1000 {
            a_with_a = true;
        }
        if c.abs() > ZERO_THRESHOLD && c.abs() <= a_max_coeff1000 {
            a_with_a = true;
        }
        if d.abs() > ZERO_THRESHOLD && d.abs() <= a_max_coeff1000 {
            a_with_a = true;
        }
        if e.abs() > ZERO_THRESHOLD && e.abs() <= a_max_coeff1000 {
            a_with_a = true;
        }
        return !a_with_a;
    }
    false
}

/// OCCT SolveFerrariResolvent (L396-428).
fn solve_ferrari_resolvent(a: f64, b: f64, c: f64, d: f64, success: &mut bool) -> f64 {
    let a_r3 = -b;
    let a_s3 = a * c - 4.0 * d;
    let a_t3 = d * (4.0 * b - a * a) - c * c;

    let a_cubic_solver = DirectPolynomialRoots::new_cubic(1.0, a_r3, a_s3, a_t3);
    if !a_cubic_solver.is_done() {
        *success = false;
        return 0.0;
    }

    *success = true;
    let mut a_y0 = a_cubic_solver.value(1);
    for i in 2..=a_cubic_solver.nb_solutions() {
        if a_cubic_solver.value(i) > a_y0 {
            a_y0 = a_cubic_solver.value(i);
        }
    }
    a_y0
}

/// OCCT FactorQuarticViaFerrari (L445-505).
struct QuarticFactorization {
    p1: f64,
    q1: f64,
    p2: f64,
    q2: f64,
}

fn factor_quartic_via_ferrari(a: f64, b: f64, c: f64, d: f64, y0: f64) -> QuarticFactorization {
    const MACHINE_EPSILON: f64 = f64::EPSILON;

    let a_discr = a * y0 * 0.5 - c;
    let a_sdiscr = if a_discr >= 0.0 { 1.0 } else { -1.0 };

    let mut a_p0 = a * a * 0.25 - b + y0;
    a_p0 = if a_p0 < 0.0 { 0.0 } else { a_p0.sqrt() };

    let mut a_q0 = y0 * y0 * 0.25 - d;
    if a_q0.abs() < 10.0 * MACHINE_EPSILON {
        a_q0 = 0.0;
    } else {
        a_q0 = if a_q0 < 0.0 { 0.0 } else { a_q0.sqrt() };
    }

    let a_ademi = a * 0.5;
    let a_ydemi = y0 * 0.5;
    let a_sdiscr_q0 = a_sdiscr * a_q0;

    let mut p1 = a_ademi + a_p0;
    let mut q1 = a_ydemi + a_sdiscr_q0;
    let mut p2 = a_ademi - a_p0;
    let mut q2 = a_ydemi - a_sdiscr_q0;

    let an_eps = 100.0 * MACHINE_EPSILON;
    if p1.abs() <= an_eps {
        p1 = 0.0;
    }
    if p2.abs() <= an_eps {
        p2 = 0.0;
    }
    if q1.abs() <= an_eps {
        q1 = 0.0;
    }
    if q2.abs() <= an_eps {
        q2 = 0.0;
    }

    QuarticFactorization { p1, q1, p2, q2 }
}
