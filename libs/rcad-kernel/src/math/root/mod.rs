//! OCCT MathRoot: root-finding algorithms.
//!
//! Corresponds to OCCT `math_FunctionRoot`, `math_FunctionSetRoot`,
//! `math_TrigonometricFunctionRoots`, `math_BissecNewton`,
//! `math_FunctionSample`, `math_FunctionRoots`, `math_FunctionAllRoots`,
//! `math_BracketedRoot`.
//!
//! - `math_FunctionRoot` — newton_raphson, bisection, secant
//! - `math_BissecNewton` — biss_newton
//! - `math_TrigonometricFunctionRoots` — trig_roots, find_roots_in, bracket_root
//! - `math_FunctionSample/Roots/AllRoots/BracketedRoot` — function_all_roots

pub mod function_all_roots;

pub use function_all_roots::{
    BracketedRoot, FunctionAllRoots, FunctionRoots, FunctionSample, FunctionValue,
    FunctionWithDerivative,
};

use std::f64::consts::TAU;

const TOL_FLOAT_DEDUP: f64 = 1e-15;

// =============================================================================
// math_FunctionRoot — 1D root finding
// =============================================================================

/// Newton-Raphson method for finding roots of a function.
///
/// Returns the root if found within tolerance and iteration limit.
pub fn newton_raphson(
    f: fn(f64) -> f64,
    df: fn(f64) -> f64,
    x0: f64,
    tol: f64,
    max_iter: usize,
) -> Option<f64> {
    let mut x = x0;
    for _ in 0..max_iter {
        let fx = f(x);
        if fx.abs() < tol {
            return Some(x);
        }
        let dfx = df(x);
        if dfx.abs() < TOL_FLOAT_DEDUP {
            return None;
        }
        let x_new = x - fx / dfx;
        if (x_new - x).abs() < tol {
            return Some(x_new);
        }
        x = x_new;
    }
    if f(x).abs() < tol { Some(x) } else { None }
}

/// Bisection method for finding roots of a function.
///
/// Returns the root if found within interval and tolerance.
pub fn bisection(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> Option<f64> {
    let mut lo = a;
    let mut hi = b;
    let f_lo = f(lo);
    let f_hi = f(hi);
    if f_lo * f_hi > 0.0 {
        return None;
    }
    if f_lo.abs() < tol {
        return Some(lo);
    }
    if f_hi.abs() < tol {
        return Some(hi);
    }
    let max_iter = ((hi - lo) / tol).ceil() as usize;
    for _ in 0..max_iter {
        let mid = (lo + hi) / 2.0;
        let f_mid = f(mid);
        if f_mid.abs() < tol || (hi - lo) / 2.0 < tol {
            return Some(mid);
        }
        if f_lo * f_mid < 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some((lo + hi) / 2.0)
}

/// Secant method for finding roots of a function.
///
/// Returns the root if found within tolerance.
pub fn secant(f: fn(f64) -> f64, x0: f64, x1: f64, tol: f64) -> Option<f64> {
    let mut x_prev = x0;
    let mut x_curr = x1;
    let max_iter = 100;
    for _ in 0..max_iter {
        let f_prev = f(x_prev);
        let f_curr = f(x_curr);
        if f_curr.abs() < tol {
            return Some(x_curr);
        }
        let denom = f_curr - f_prev;
        if denom.abs() < TOL_FLOAT_DEDUP {
            return None;
        }
        let x_new = x_curr - f_curr * (x_curr - x_prev) / denom;
        if (x_new - x_curr).abs() < tol {
            return Some(x_new);
        }
        x_prev = x_curr;
        x_curr = x_new;
    }
    if f(x_curr).abs() < tol { Some(x_curr) } else { None }
}

// =============================================================================
// math_BissecNewton — hybrid bisection + Newton root finding
// =============================================================================

/// Hybrid bisection-Newton method (safe + fast). OCCT `math_BissecNewton`.
pub fn biss_newton(
    f: impl Fn(f64) -> f64,
    df: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    tol: f64,
) -> Option<f64> {
    let mut lo = a.min(b);
    let mut hi = a.max(b);
    let mut x = (lo + hi) / 2.0;
    for _ in 0..100 {
        let fx = f(x);
        if fx.abs() < tol {
            return Some(x);
        }
        let dfx = df(x);
        let xn = if dfx.abs() > TOL_FLOAT_DEDUP {
            let nv = x - fx / dfx;
            if nv > lo && nv < hi { nv } else { (lo + hi) / 2.0 }
        } else {
            (lo + hi) / 2.0
        };
        if (xn - x).abs() < tol {
            return Some(xn);
        }
        if f(lo) * f(xn) <= 0.0 {
            hi = xn;
        } else {
            lo = xn;
        }
        x = xn;
    }
    if f(x).abs() < tol { Some(x) } else { None }
}

// =============================================================================
// math_TrigonometricFunctionRoots — trigonometric equation root finding
// =============================================================================

fn trig_f(a: f64, b: f64, c: f64, d: f64, e: f64) -> impl Fn(f64) -> f64 {
    move |x: f64| {
        let (s, c_) = x.sin_cos();
        a * c_ * c_ + 2.0 * b * c_ * s + c * c_ + d * s + e
    }
}

/// Find roots of trigonometric equation in [x_min, x_max].
/// Solves: a·cos²(x) + 2b·cos(x)sin(x) + c·cos(x) + d·sin(x) + e = 0
pub fn trig_roots(a: f64, b: f64, c: f64, d: f64, e: f64, x_min: f64, x_max: f64) -> Vec<f64> {
    let f = trig_f(a, b, c, d, e);
    find_roots_in(f, x_min, x_max, 200)
}

/// Solve d·sin(x) + e = 0 in [x_min, x_max].
pub fn trig_roots_sin_only(d: f64, e: f64, x_min: f64, x_max: f64) -> Vec<f64> {
    trig_roots(0.0, 0.0, 0.0, d, e, x_min, x_max)
}

/// Solve c·cos(x) + d·sin(x) + e = 0 in [x_min, x_max].
pub fn trig_roots_cos_sin(c: f64, d: f64, e: f64, x_min: f64, x_max: f64) -> Vec<f64> {
    trig_roots(0.0, 0.0, c, d, e, x_min, x_max)
}

/// Find all roots of f(x) = 0 in [a, b] by scanning for sign changes.
pub fn find_roots_in(f: impl Fn(f64) -> f64, a: f64, b: f64, n_intervals: usize) -> Vec<f64> {
    let step = (b - a) / n_intervals.max(1) as f64;
    let mut roots = Vec::new();
    for i in 0..n_intervals {
        let x1 = a + i as f64 * step;
        let x2 = x1 + step;
        let f1 = f(x1);
        let f2 = f(x2);
        if f1 * f2 < 0.0 {
            if let Some(r) = bisection(&f, x1, x2, 1e-10) {
                roots.push(r);
            }
        } else if f1.abs() < 1e-10 && !roots.iter().any(|r| (r - x1).abs() < 1e-8) {
            roots.push(x1);
        }
    }
    if f(b).abs() < 1e-10 && !roots.iter().any(|r| (r - b).abs() < 1e-8) {
        roots.push(b);
    }
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

/// Find a bracket [a,b] where f(a) and f(b) have opposite signs.
pub fn bracket_root(
    f: impl Fn(f64) -> f64,
    x0: f64,
    step: f64,
    max_steps: usize,
) -> Option<(f64, f64)> {
    let mut a = x0;
    let mut fa = f(a);
    for _ in 0..max_steps {
        let b = a + step;
        let fb = f(b);
        if fa * fb <= 0.0 {
            return Some((a, b));
        }
        a = b;
        fa = fb;
    }
    None
}

// =============================================================================
// math_TrigonometricFunctionRoots — trigonometric polynomial roots (OCCT 1:1)
// =============================================================================

/// Result of [`trig_function_roots`] (OCCT math_TrigonometricFunctionRoots).
pub struct TrigRoots {
    pub done: bool,
    pub infinite: bool,
    /// Sorted roots within [inf_bound, sup_bound].
    pub roots: Vec<f64>,
}

/// OCCT math_TrigonometricFunctionRoots::Perform (math_TrigonometricFunctionRoots.cxx L75-551).
///
/// Solves `a·cos²(x) + 2b·sin(x)cos(x) + c·cos(x) + d·sin(x) + e = 0` in
/// `[inf_bound, sup_bound]`.
pub fn trig_function_roots(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    e: f64,
    inf_bound: f64,
    sup_bound: f64,
) -> TrigRoots {
    use std::f64::consts::PI;

    let depi = TAU;
    let eps = 1.5e-12;
    let tol1 = 1e-15;
    let nit = 10;

    let mut infinite = false;
    let mut done = true;
    let mut sol: Vec<f64> = Vec::new(); // up to 4, sorted on insert
    let mut zer = [0.0f64; 4];
    let mut n_zer = 0usize;
    let mut ko = [0.0f64; 5];

    // Bounds setup (OCCT L95-123)
    let (my_borne_inf, delta, mod0) = if inf_bound <= f64::MIN && sup_bound >= f64::MAX {
        (0.0, depi, 0.0)
    } else if sup_bound >= f64::MAX {
        (inf_bound, depi, (inf_bound / depi).trunc())
    } else if inf_bound <= f64::MIN {
        (sup_bound - depi, depi, ((sup_bound - depi) / depi).trunc())
    } else {
        let mut delta = sup_bound - inf_bound;
        if delta > depi {
            delta = depi;
        }
        (inf_bound, delta, (inf_bound / depi).trunc())
    };

    // OCCT L125-237: A==0 && B==0
    if a.abs() <= eps && b.abs() <= eps {
        if c.abs() <= eps {
            if d.abs() <= eps {
                if e.abs() <= eps {
                    return TrigRoots { done: true, infinite: true, roots: sol };
                } else {
                    return TrigRoots { done: true, infinite: false, roots: sol };
                }
            }
            // d*sin(x) + e = 0  (OCCT L144-172)
            let aa = -e / d;
            if aa.abs() > 1.0 {
                return TrigRoots { done: true, infinite: false, roots: sol };
            }
            zer[0] = aa.asin();
            zer[1] = PI - zer[0];
            n_zer = 2;
            for i in 0..n_zer {
                if zer[i] <= -eps {
                    zer[i] = depi - zer[i].abs();
                }
                zer[i] += mod0 * depi;
                let x = zer[i] - my_borne_inf;
                if x > -1e-12 && x < delta + 1e-12 {
                    sol.push(zer[i]);
                }
            }
            return TrigRoots { done: true, infinite: false, roots: sol };
        } else if d.abs() <= eps {
            // c*cos(x) + e = 0  (OCCT L175-207)
            let aa = -e / c;
            if aa.abs() > 1.0 {
                return TrigRoots { done: true, infinite: false, roots: sol };
            }
            zer[0] = aa.acos();
            zer[1] = -zer[0];
            n_zer = 2;
            for i in 0..n_zer {
                if zer[i] <= -eps {
                    zer[i] = depi - zer[i].abs();
                }
                zer[i] += mod0 * TAU;
                let x = zer[i] - my_borne_inf;
                if x >= -1e-12 && x <= delta + 1e-12 {
                    sol.push(zer[i]);
                }
            }
            return TrigRoots { done: true, infinite: false, roots: sol };
        } else {
            // quadratic (OCCT L211-236): AA*t² + BB*t + CC = 0, t = tan(x/2)
            let aa = e - c;
            let bb = 2.0 * d;
            let cc = e + c;
            let roots_t = crate::math::math_poly::solve_quadratic(aa, bb, cc);
            n_zer = roots_t.len();
            for (i, t) in roots_t.iter().enumerate() {
                zer[i] = *t;
            }
            if roots_t.is_empty() {
                return TrigRoots { done: true, infinite: false, roots: sol };
            }
        }
    } else {
        // OCCT L240-353: two additional analytical cases (A==0 && E==0)
        if a.abs() <= eps && e.abs() <= eps {
            if c.abs() <= eps {
                // 2*B*sin*cos + D*sin = 0  (OCCT L243-296)
                n_zer = 2;
                zer[0] = 0.0;
                zer[1] = PI;
                let aa = -d / (b * 2.0);
                if aa.abs() <= 1.0 + 1e-9 {
                    if aa >= 1.0 {
                        zer[2] = 0.0;
                        zer[3] = 0.0;
                    } else if aa <= -1.0 {
                        zer[2] = PI;
                        zer[3] = PI;
                    } else {
                        zer[2] = aa.acos();
                        zer[3] = depi - zer[2];
                    }
                    n_zer = 4;
                }
                for i in 0..n_zer {
                    if zer[i] <= my_borne_inf - eps {
                        zer[i] += depi;
                    }
                    zer[i] += mod0 * TAU;
                    let x = zer[i] - my_borne_inf;
                    if x >= -1e-9 && x <= delta + 1e-9 {
                        if zer[i] < inf_bound { zer[i] = inf_bound; }
                        if zer[i] > sup_bound { zer[i] = sup_bound; }
                        sol.push(zer[i]);
                    }
                }
                return TrigRoots { done: true, infinite: false, roots: sol };
            }
            if d.abs() <= eps {
                // 2*B*sin*cos + C*cos = 0  (OCCT L298-353)
                n_zer = 2;
                zer[0] = PI / 2.0;
                zer[1] = PI * 3.0 / 2.0;
                let aa = -c / (b * 2.0);
                if aa.abs() <= 1.0 + 1e-9 {
                    if aa >= 1.0 {
                        zer[2] = PI / 2.0;
                        zer[3] = PI / 2.0;
                    } else if aa <= -1.0 {
                        zer[2] = PI * 3.0 / 2.0;
                        zer[3] = PI * 3.0 / 2.0;
                    } else {
                        zer[2] = aa.asin();
                        zer[3] = PI - zer[2];
                    }
                    n_zer = 4;
                }
                for i in 0..n_zer {
                    if zer[i] <= my_borne_inf - eps {
                        zer[i] += depi;
                    }
                    zer[i] += mod0 * TAU;
                    let x = zer[i] - my_borne_inf;
                    if x >= -1e-9 && x <= delta + 1e-9 {
                        if zer[i] < inf_bound { zer[i] = inf_bound; }
                        if zer[i] > sup_bound { zer[i] = sup_bound; }
                        sol.push(zer[i]);
                    }
                }
                return TrigRoots { done: true, infinite: false, roots: sol };
            }
        }

        // General quartic (OCCT L356-435)
        ko[0] = a - c + e;
        ko[1] = 2.0 * d - 4.0 * b;
        ko[2] = 2.0 * e - 2.0 * a;
        ko[3] = 4.0 * b + 2.0 * d;
        ko[4] = a + c + e;
        let mut bko;
        let mut iterations = 0;
        loop {
            bko = false;
            let roots_t = crate::math::math_poly::solve_quartic(ko[0], ko[1], ko[2], ko[3], ko[4]);
            n_zer = roots_t.len();
            for (i, t) in roots_t.iter().enumerate() {
                zer[i] = *t;
            }
            if roots_t.is_empty() {
                break;
            }
            // sort (OCCT L387-400)
            zer[..n_zer].sort_by(|a, b| a.partial_cmp(b).unwrap());
            // dedup double roots via derivative check (OCCT L402-433)
            for i in 0..n_zer.saturating_sub(1) {
                if (zer[i + 1] - zer[i]).abs() < eps {
                    let qw = zer[i + 1];
                    let va = ko[3] + qw * (2.0 * ko[2] + qw * (3.0 * ko[1] + qw * (4.0 * ko[0])));
                    if va.abs() > eps {
                        bko = true;
                        break;
                    }
                }
            }
            if bko {
                // scale coefficients down (OCCT L427-433)
                for v in ko.iter_mut() {
                    *v *= 0.0001;
                }
                iterations += 1;
                if iterations > 4 {
                    break;
                }
            } else {
                break;
            }
        }
    }

    // Verification against bounds + Newton refinement (OCCT L437-504)
    let sup_min_inf_100 = (sup_bound - inf_bound) * 0.01;
    let trig_f = |x: f64| {
        let (sn, cs) = x.sin_cos();
        cs * (a * cs + 2.0 * b * sn + c) + d * sn + e
    };
    let trig_df = |x: f64| {
        let (sn, cs) = x.sin_cos();
        -2.0 * a * sn * cs + 2.0 * b * (cs * cs - sn * sn) - c * sn + d * cs
    };
    let mut n_sol = 0usize;
    let mut tmp = [0.0f64; 4];
    for i in 0..n_zer {
        let mut teta = 2.0 * zer[i].atan();
        if zer[i] <= -eps {
            teta = depi - teta.abs();
        }
        teta += mod0 * depi;
        if teta - my_borne_inf < 0.0 {
            teta += depi;
        }
        let x = teta - my_borne_inf;
        if x >= -1e-12 && x <= delta + 1e-12 {
            // Newton refinement (OCCT L460-478: math_NewtonFunctionRoot)
            let mut teta_newton = teta;
            let mut x_n = teta;
            for _ in 0..nit {
                let fv = trig_f(x_n);
                let dfv = trig_df(x_n);
                if dfv.abs() < 1e-30 {
                    break;
                }
                let dx = fv / dfv;
                x_n -= dx;
                if dx.abs() < tol1 {
                    break;
                }
            }
            teta_newton = x_n;
            let d_newton = teta_newton - teta;
            if d_newton <= sup_min_inf_100 && d_newton >= -sup_min_inf_100 {
                teta = teta_newton;
            }
            // insert sorted (OCCT L480-502)
            let mut inserted = false;
            for k in 0..n_sol {
                if teta < tmp[k] {
                    for l in (k..n_sol).rev() {
                        tmp[l + 1] = tmp[l];
                    }
                    tmp[k] = teta;
                    n_sol += 1;
                    inserted = true;
                    break;
                }
            }
            if !inserted {
                if n_sol < 4 {
                    tmp[n_sol] = teta;
                    n_sol += 1;
                }
            }
        }
    }
    // Special case x = PI (OCCT L506-550)
    if n_sol < 4 && (a - c + e).abs() <= eps {
        let teta = PI + mod0 * TAU;
        let x = teta - my_borne_inf;
        if x >= -1e-12 && x <= delta + 1e-12 {
            let mut j = 0usize;
            let mut found = false;
            for k in 0..n_sol {
                if teta < tmp[k] {
                    found = true;
                    break;
                }
                if (teta - tmp[k]).abs() <= eps {
                    return TrigRoots { done: true, infinite: false, roots: tmp[..n_sol].to_vec() };
                }
                j = k + 1;
            }
            if !found {
                if n_sol < 4 {
                    tmp[n_sol] = teta;
                    n_sol += 1;
                }
            } else {
                for k in (j..n_sol).rev() {
                    tmp[k + 1] = tmp[k];
                }
                tmp[j] = teta;
                n_sol += 1;
            }
        }
    }

    TrigRoots { done, infinite, roots: tmp[..n_sol].to_vec() }
}
