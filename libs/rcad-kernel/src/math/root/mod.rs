//! OCCT MathRoot: root-finding algorithms.
//!
//! Corresponds to OCCT `math_FunctionRoot`, `math_FunctionSetRoot`,
//! `math_TrigonometricFunctionRoots`, `math_BissecNewton`.
//!
//! - `math_FunctionRoot` — newton_raphson, bisection, secant
//! - `math_BissecNewton` — biss_newton
//! - `math_TrigonometricFunctionRoots` — trig_roots, find_roots_in, bracket_root

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
