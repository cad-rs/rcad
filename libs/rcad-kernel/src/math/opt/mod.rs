//! OCCT MathOpt: optimization and minimization algorithms.
//!
//! Corresponds to OCCT `math_BFGS`, `math_FRPR`, `math_NewtonMinimum`,
//! `math_Powell`, `math_BrentMinimum`, `math_BracketMinimum`,
//! `math_GlobOptMin`, `math_PSO`, `math_LevenbergMarquardt`.

const TOL_FLOAT_DEDUP: f64 = 1e-15;
const PHI: f64 = 1.618033988749895;
const RESPHI: f64 = 0.3819660112501051; // 1/phi^2

// =============================================================================
// math_BFGS — Broyden-Fletcher-Goldfarb-Shanno quasi-Newton optimization
// =============================================================================

fn line_search_backtracking(
    x: &[f64],
    p: &[f64],
    _grad: &[f64],
    f_current: f64,
    f_grad: impl Fn(&[f64], &mut [f64]) -> f64,
) -> f64 {
    let c = 1e-4;
    let rho = 0.5;
    let mut alpha = 1.0;
    let mut trial = vec![0.0; x.len()];
    let mut g = vec![0.0; x.len()];
    for _ in 0..20 {
        for i in 0..x.len() {
            trial[i] = x[i] + alpha * p[i];
        }
        let f_trial = f_grad(&trial, &mut g);
        if f_trial
            <= f_current + c * alpha * p.iter().zip(_grad.iter()).map(|(p, g)| p * g).sum::<f64>()
        {
            return alpha;
        }
        alpha *= rho;
    }
    alpha
}

/// Minimize using BFGS quasi-Newton method. Returns minimizer or None.
pub fn bfgs_minimize(
    x0: &[f64],
    f_grad: impl Fn(&[f64], &mut [f64]) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut grad = vec![0.0; n];
    let mut f = f_grad(&x, &mut grad);
    let mut h_inv = vec![0.0; n * n];
    for i in 0..n {
        h_inv[i * n + i] = 1.0;
    }
    for _ in 0..max_iter {
        let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gn < tol {
            return Some(x);
        }
        let mut p = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                p[i] -= h_inv[i * n + j] * grad[j];
            }
        }
        let alpha = line_search_backtracking(&x, &p, &grad, f, &f_grad);
        let mut s = vec![0.0; n];
        for i in 0..n {
            s[i] = alpha * p[i];
            x[i] += s[i];
        }
        let mut new_grad = vec![0.0; n];
        let new_f = f_grad(&x, &mut new_grad);
        let mut y = vec![0.0; n];
        for i in 0..n {
            y[i] = new_grad[i] - grad[i];
        }
        grad = new_grad;
        f = new_f;
        let sy = s.iter().zip(y.iter()).map(|(s, y)| s * y).sum::<f64>();
        if sy.abs() < TOL_FLOAT_DEDUP {
            continue;
        }
        let rho = 1.0 / sy;
        let mut hy = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                hy[i] += h_inv[i * n + j] * y[j];
            }
        }
        let ythy = y.iter().zip(hy.iter()).map(|(y, hy)| y * hy).sum::<f64>();
        let factor = 1.0 + rho * ythy;
        let mut h_new = h_inv.clone();
        for i in 0..n {
            for j in 0..n {
                h_new[i * n + j] += rho * (factor * s[i] * s[j] - s[i] * hy[j] - hy[i] * s[j]);
            }
        }
        h_inv = h_new;
    }
    let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
    if gn < tol { Some(x) } else { None }
}

// =============================================================================
// math_FRPR — Fletcher-Reeves Polak-Ribiere conjugate gradient optimization
// =============================================================================

/// Minimize using FRPR conjugate gradient method.
pub fn frpr_minimize(
    x0: &[f64],
    f_grad: impl Fn(&[f64], &mut [f64]) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut grad = vec![0.0; n];
    let _f = f_grad(&x, &mut grad);
    let mut d = grad.iter().map(|g| -g).collect::<Vec<_>>();
    let mut prev_norm2 = grad.iter().map(|g| g * g).sum::<f64>();

    for _ in 0..max_iter {
        if prev_norm2.sqrt() < tol {
            return Some(x);
        }
        let alpha = line_search_backtracking(&x, &d, &grad, _f, &f_grad);
        for i in 0..n {
            x[i] += alpha * d[i];
        }
        let mut new_grad = vec![0.0; n];
        let _f = f_grad(&x, &mut new_grad);
        let new_norm2 = new_grad.iter().map(|g| g * g).sum::<f64>();
        let gg_diff: f64 = new_grad
            .iter()
            .zip(grad.iter())
            .map(|(ng, g)| ng * (ng - g))
            .sum();
        let beta = if prev_norm2 > TOL_FLOAT_DEDUP {
            (gg_diff / prev_norm2).max(0.0)
        } else {
            0.0
        };
        for i in 0..n {
            d[i] = -new_grad[i] + beta * d[i];
        }
        grad = new_grad;
        prev_norm2 = new_norm2;
    }
    if prev_norm2.sqrt() < tol { Some(x) } else { None }
}

// =============================================================================
// math_NewtonMinimum — Newton's method with Hessian
// =============================================================================

fn solve_linear_system_gauss(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut aug = vec![0.0; n * (n + 1)];
    for i in 0..n {
        for j in 0..n {
            aug[i * (n + 1) + j] = a[i * n + j];
        }
        aug[i * (n + 1) + n] = b[i];
    }
    for col in 0..n {
        let mut max_row = col;
        for row in (col + 1)..n {
            if aug[row * (n + 1) + col].abs() > aug[max_row * (n + 1) + col].abs() {
                max_row = row;
            }
        }
        if aug[max_row * (n + 1) + col].abs() < TOL_FLOAT_DEDUP {
            return None;
        }
        if max_row != col {
            for j in col..=n {
                aug.swap(col * (n + 1) + j, max_row * (n + 1) + j);
            }
        }
        for row in (col + 1)..n {
            let factor = aug[row * (n + 1) + col] / aug[col * (n + 1) + col];
            for j in col..=n {
                aug[row * (n + 1) + j] -= factor * aug[col * (n + 1) + j];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i * (n + 1) + n];
        for j in (i + 1)..n {
            sum -= aug[i * (n + 1) + j] * x[j];
        }
        x[i] = sum / aug[i * (n + 1) + i];
    }
    Some(x)
}

/// Minimize using Newton's method with gradient and Hessian.
pub fn newton_minimize(
    x0: &[f64],
    f_grad_hess: impl Fn(&[f64], &mut [f64], &mut [f64]) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut grad = vec![0.0; n];
    let mut hess = vec![0.0; n * n];
    let mut f = f_grad_hess(&x, &mut grad, &mut hess);

    for _ in 0..max_iter {
        let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gn < tol {
            return Some(x);
        }
        if let Some(p) = solve_linear_system_gauss(&hess, &grad.iter().map(|g| -g).collect::<Vec<_>>(), n) {
            let mut alpha = 1.0;
            let mut trial = x.clone();
            let mut g_trial = vec![0.0; n];
            for _ in 0..20 {
                for i in 0..n {
                    trial[i] = x[i] + alpha * p[i];
                }
                let f_trial = f_grad_hess(&trial, &mut g_trial, &mut hess);
                let directional = p.iter().zip(grad.iter()).map(|(p, g)| p * g).sum::<f64>();
                if f_trial <= f + 1e-4 * alpha * directional {
                    break;
                }
                alpha *= 0.5;
            }
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            f = f_grad_hess(&x, &mut grad, &mut hess);
        } else {
            return None;
        }
    }
    let gn = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
    if gn < tol { Some(x) } else { None }
}

// =============================================================================
// math_Powell — derivative-free direction-set optimization
// =============================================================================

fn minimize_1d(x: &[f64], dir: &[f64], f: impl Fn(&[f64]) -> f64) -> (Vec<f64>, f64) {
    let n = x.len();
    let at = |alpha: f64| -> f64 {
        let mut t = x.to_vec();
        for i in 0..n {
            t[i] += alpha * dir[i];
        }
        f(&t)
    };
    let alpha = golden_section_min(at, 0.0, 10.0, 1e-10);
    let mut xn = x.to_vec();
    for i in 0..n {
        xn[i] += alpha * dir[i];
    }
    let fxn = f(&xn);
    (xn, fxn)
}

/// Minimize using Powell's derivative-free conjugate direction method.
pub fn powell_minimize(
    x0: &[f64],
    f: impl Fn(&[f64]) -> f64,
    tol: f64,
    max_iter: usize,
) -> Option<Vec<f64>> {
    let n = x0.len();
    let mut x = x0.to_vec();
    let mut dirs: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut d = vec![0.0; n];
            d[i] = 1.0;
            d
        })
        .collect();
    let mut prev_f = f(&x);
    for _ in 0..max_iter {
        let x_start = x.clone();
        let mut delta = 0.0;
        let mut best_dir = 0;
        for i in 0..n {
            let (xn, fn_) = minimize_1d(&x, &dirs[i], &f);
            let dec = prev_f - fn_;
            if dec > delta {
                delta = dec;
                best_dir = i;
            }
            x = xn;
            prev_f = fn_;
        }
        if x.iter()
            .zip(x_start.iter())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>() < tol
        {
            return Some(x);
        }
        let mut nd = Vec::with_capacity(n);
        for i in 0..n {
            nd.push(x[i] - x_start[i]);
        }
        let nn = nd.iter().map(|d| d * d).sum::<f64>().sqrt();
        if nn > TOL_FLOAT_DEDUP {
            for i in 0..n {
                nd[i] /= nn;
            }
            for i in best_dir..n - 1 {
                dirs.swap(i, i + 1);
            }
            *dirs.last_mut().unwrap() = nd;
        }
    }
    None
}

// =============================================================================
// math_BracketMinimum — golden section search
// =============================================================================

/// Golden section search for finding minimum. OCCT `math_BracketMinimum`.
pub fn golden_section_min<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> f64 {
    let mut lo = a;
    let mut hi = b;
    let mut c = lo + RESPHI * (hi - lo);
    let mut d = hi - RESPHI * (hi - lo);
    let mut fc = f(c);
    let mut fd = f(d);
    while (hi - lo).abs() > tol {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = lo + RESPHI * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = hi - RESPHI * (hi - lo);
            fd = f(d);
        }
    }
    (lo + hi) / 2.0
}

/// Golden section search for finding maximum.
pub fn golden_section_max<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64) -> f64 {
    golden_section_min(|x| -f(x), a, b, tol)
}

// =============================================================================
// math_BrentMinimum — Brent's method for 1D minimization
// =============================================================================

/// Minimize a 1D function using Brent's method. OCCT `math_BrentMinimum`.
///
/// `a` and `b` are the bracket endpoints; the initial guess is their midpoint.
/// OCCT-aligned: `math_BrentMinimum::Perform(F, Ax, Bx, Cx)` with the bracketing
/// triplet `(a, (a+b)/2, b)`.
pub fn brent_minimize(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> f64 {
    // OCCT math_BrentMinimum.cxx: CGOLD = 0.5 * (3 - sqrt(5)).
    const CGOLD: f64 = 0.3819660;
    // OCCT: EPSZ = 1e-12, NbIterations = 100.
    const EPSZ: f64 = 1e-12;
    const ITERMAX: usize = 100;

    let (ax, cx) = if a < b { (a, b) } else { (b, a) };
    let (mut lo, mut hi) = if ax < cx { (ax, cx) } else { (cx, ax) };
    let bx = (ax + cx) * 0.5;
    let (mut x, mut w, mut v) = (bx, bx, bx);
    let fx0 = f(x);
    let (mut fx, mut fw, mut fv) = (fx0, fx0, fx0);
    let (mut e, mut d) = (0.0f64, f64::MAX); // OCCT: d = RealLast()

    for _ in 1..=ITERMAX {
        let xm = 0.5 * (lo + hi);
        let tol1 = tol * x.abs() + EPSZ;
        let tol2 = 2.0 * tol1;
        // OCCT IsSolutionReached: x <= 2*tol1 + a && x >= b - 2*tol1.
        if x <= tol2 + lo && x >= hi - tol2 {
            return x;
        }
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let q = (x - v) * (fx - fw);
            let p = (x - v) * q - (x - w) * r;
            let mut q2 = 2.0 * (q - r);
            let mut p2 = p;
            if q2 > 0.0 {
                p2 = -p2;
            }
            q2 = q2.abs();
            let etemp = e;
            e = d;
            if p2.abs() >= (0.5 * q2 * etemp).abs() || p2 <= q2 * (lo - x) || p2 >= q2 * (hi - x)
            {
                e = if x >= xm { lo - x } else { hi - x };
                d = CGOLD * e;
            } else {
                d = p2 / q2;
                let u = x + d;
                if u - lo < tol2 || hi - u < tol2 {
                    d = tol1.copysign(xm - x);
                }
            }
        } else {
            e = if x >= xm { lo - x } else { hi - x };
            d = CGOLD * e;
        }
        let u = if d.abs() >= tol1 { x + d } else { x + tol1.copysign(d) };
        let fu = f(u);
        if fu <= fx {
            if u >= x {
                lo = x;
            } else {
                hi = x;
            }
            // SHFT(v, w, x, u) and SHFT(fv, fw, fx, fu).
            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if u < x {
                lo = u;
            } else {
                hi = u;
            }
            if fu <= fw || w == x {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    }
    x
}

// =============================================================================
// math_GlobOptMin — global optimization via grid + local refinement
// =============================================================================

/// Global optimizer: coarse grid evaluation + local coordinate descent.
pub fn glob_opt_min(
    f: impl Fn(&[f64]) -> f64,
    lower: &[f64],
    upper: &[f64],
    n_cells: usize,
    n_local: usize,
) -> Vec<f64> {
    let n = lower.len();
    if n == 0 { return vec![]; }
    let mut candidates: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut current = vec![0.0; n];
    fn grid_eval(
        f: &impl Fn(&[f64]) -> f64,
        lower: &[f64], upper: &[f64],
        nc: usize, cand: &mut Vec<(f64, Vec<f64>)>,
        cur: &mut Vec<f64>, dim: usize,
    ) {
        if dim == cur.len() {
            cand.push((f(cur), cur.clone()));
            return;
        }
        let step = (upper[dim] - lower[dim]) / nc as f64;
        for i in 0..=nc {
            cur[dim] = lower[dim] + i as f64 * step;
            grid_eval(f, lower, upper, nc, cand, cur, dim + 1);
        }
    }
    grid_eval(&f, lower, upper, n_cells, &mut candidates, &mut current, 0);
    candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    candidates.truncate(n_local.max(1));
    let mut best_x = candidates[0].1.clone();
    let mut best_f = candidates[0].0;
    for (_, mut x) in candidates {
        for _ in 0..5 {
            for i in 0..n {
                let step = (upper[i] - lower[i]) * 0.01;
                let fx = f(&x);
                let mut xt = x.clone();
                xt[i] += step;
                if xt[i] <= upper[i] && f(&xt) < fx { x[i] = xt[i]; continue; }
                xt[i] = x[i] - step;
                if xt[i] >= lower[i] && f(&xt) < fx { x[i] = xt[i]; }
            }
        }
        let fv = f(&x);
        if fv < best_f { best_f = fv; best_x = x; }
    }
    best_x
}

// =============================================================================
// math_PSO — Particle Swarm Optimization
// =============================================================================

/// Minimize using Particle Swarm Optimization.
pub fn pso_minimize(
    f: impl Fn(&[f64]) -> f64,
    lower: &[f64], upper: &[f64],
    n_particles: usize, max_iter: usize, tol: f64,
) -> Vec<f64> {
    let n = lower.len();
    // OCCT math_PSO uses math_BullardGenerator with default seed 1 (deterministic).
    let mut rng = fastrand::Rng::with_seed(1);
    let mut positions = Vec::with_capacity(n_particles);
    let mut velocities = Vec::with_capacity(n_particles);
    let mut pbest = Vec::with_capacity(n_particles);
    let mut pbest_val = Vec::with_capacity(n_particles);
    for _ in 0..n_particles {
        let mut pos = Vec::with_capacity(n);
        for i in 0..n {
            pos.push(lower[i] + rng.f64() * (upper[i] - lower[i]));
        }
        let fv = f(&pos);
        positions.push(pos.clone());
        velocities.push(vec![0.0; n]);
        pbest.push(pos);
        pbest_val.push(fv);
    }
    let mut gbest = pbest[0].clone();
    let mut gbest_val = pbest_val[0];
    for i in 1..n_particles {
        if pbest_val[i] < gbest_val { gbest = pbest[i].clone(); gbest_val = pbest_val[i]; }
    }
    for _ in 0..max_iter {
        let prev = gbest_val;
        for i in 0..n_particles {
            for j in 0..n {
                velocities[i][j] = 0.72 * velocities[i][j]
                    + 1.49 * rng.f64() * (pbest[i][j] - positions[i][j])
                    + 1.49 * rng.f64() * (gbest[j] - positions[i][j]);
                let vmax = (upper[j] - lower[j]) * 0.2;
                velocities[i][j] = velocities[i][j].clamp(-vmax, vmax);
                positions[i][j] = (positions[i][j] + velocities[i][j]).clamp(lower[j], upper[j]);
            }
            let fv = f(&positions[i]);
            if fv < pbest_val[i] { pbest_val[i] = fv; pbest[i] = positions[i].clone(); }
            if fv < gbest_val { gbest_val = fv; gbest = positions[i].clone(); }
        }
        if (prev - gbest_val).abs() < tol && gbest_val.abs() > 1e-15 {
            break;
        }
    }
    gbest
}

// =============================================================================
// math_LevenbergMarquardt — nonlinear least squares
// =============================================================================

/// Levenberg-Marquardt solver for nonlinear least squares.
pub fn lm_solve(
    x0: &[f64],
    mut func: impl FnMut(&[f64], &mut [f64], &mut [f64]) -> f64,
    n_eq: usize, max_iter: usize, tol: f64,
) -> Option<Vec<f64>> {
    let n = x0.len();
    if n == 0 || n_eq == 0 { return None; }
    let mut x = x0.to_vec();
    let mut f = vec![0.0; n_eq];
    let mut jac = vec![0.0; n_eq * n];
    let mut lambda = 1.0;
    let mut cost = func(&x, &mut f, &mut jac);
    for _ in 0..max_iter {
        let mut jtj = vec![0.0; n * n];
        let mut jtf_vec = vec![0.0; n];
        for i in 0..n {
            for k in 0..n {
                let mut s = 0.0;
                for r in 0..n_eq {
                    s += jac[r * n + i] * jac[r * n + k];
                }
                jtj[i * n + k] = s;
            }
        }
        for i in 0..n {
            let mut s = 0.0;
            for r in 0..n_eq {
                s += jac[r * n + i] * f[r];
            }
            jtf_vec[i] = s;
        }
        if jtf_vec.iter().map(|v| v * v).sum::<f64>().sqrt() < tol {
            return Some(x);
        }
        let mut h = jtj.clone();
        for i in 0..n {
            h[i * n + i] += lambda;
        }
        let rhs: Vec<f64> = jtf_vec.iter().map(|v| -v).collect();
        let delta = solve_linear_system_gauss(&h, &rhs, n);
        let (cost_new, gain_ratio) = match delta {
            Some(ref d) => {
                let mut xn = vec![0.0; n];
                for i in 0..n {
                    xn[i] = x[i] + d[i];
                }
                let mut fn_ = vec![0.0; n_eq];
                let mut jn = vec![0.0; n_eq * n];
                let cn = func(&xn, &mut fn_, &mut jn);
                let pred: f64 = 0.5 * d.iter().zip(rhs.iter()).map(|(d, r)| d * r).sum::<f64>();
                let gr = if pred.abs() > 1e-15 { (cost - cn) / pred } else { 0.0 };
                (cn, gr)
            }
            None => (cost, -1.0),
        };
        if gain_ratio > 0.0 {
            let d = delta.as_ref().unwrap();
            for i in 0..n {
                x[i] += d[i];
            }
            let _ = func(&x, &mut f, &mut jac);
            cost = cost_new;
            lambda *= if gain_ratio > 0.75 { 0.5 } else if gain_ratio < 0.25 { 2.0 } else { 1.0 };
            lambda = lambda.max(1e-12).min(1e12);
            let dn: f64 = d.iter().map(|d| d * d).sum::<f64>().sqrt();
            if dn < tol { return Some(x); }
        } else {
            lambda *= 2.0;
            if lambda > 1e12 { break; }
        }
    }
    if cost < tol { Some(x) } else { None }
}

// =============================================================================
// math_BrentMinimum — OCCT 1:1 three-parameter Brent minimization
// =============================================================================

/// OCCT math_BrentMinimum (math_BrentMinimum.hxx/.cxx/.lxx) — Brent's method to
/// find the minimum of a function of a single variable.  No knowledge of the
/// derivative is required.
///
/// Unlike the simplified [`brent_minimize`] above (which always uses the
/// bracket midpoint as the initial guess), `Perform(F, Ax, Bx, Cx)` takes the
/// bracketing triplet (Ax, Bx, Cx) explicitly: Bx is the initial guess, not the
/// midpoint.  Used by IntStart_SearchOnBoundaries tangent refinement.
pub struct BrentMinimum {
    a: f64,
    b: f64,
    x: f64,
    fx: f64,
    fv: f64,
    fw: f64,
    x_tol: f64,
    eps_z: f64,
    done: bool,
    iter: i32,
    iter_max: i32,
    my_f: bool,
}

impl BrentMinimum {
    /// OCCT math_BrentMinimum(TolX, NbIterations = 100, ZEPS = 1.0e-12).
    pub fn new(tol_x: f64, nb_iterations: i32, zeps: f64) -> Self {
        BrentMinimum {
            a: 0.0,
            b: 0.0,
            x: 0.0,
            fx: 0.0,
            fv: 0.0,
            fw: 0.0,
            x_tol: tol_x,
            eps_z: zeps,
            done: false,
            iter: 0,
            iter_max: nb_iterations,
            my_f: false,
        }
    }

    /// OCCT math_BrentMinimum(TolX, Fbx, NbIterations = 100, ZEPS = 1.0e-12).
    pub fn new_with_fbx(tol_x: f64, f_bx: f64, nb_iterations: i32, zeps: f64) -> Self {
        BrentMinimum {
            a: 0.0,
            b: 0.0,
            x: 0.0,
            fx: f_bx,
            fv: 0.0,
            fw: 0.0,
            x_tol: tol_x,
            eps_z: zeps,
            done: false,
            iter: 0,
            iter_max: nb_iterations,
            my_f: true,
        }
    }

    /// OCCT Perform(F, Ax, Bx, Cx) — Brent minimization on function F from a
    /// bracketing triplet (Ax, Bx, Cx) with Bx between Ax and Cx.
    /// F is a math_Function (Value only), like OCCT math_BrentMinimum::Perform.
    pub fn perform(&mut self, f: &mut dyn crate::math::root::FunctionValue, ax: f64, bx: f64, cx: f64) {
        let cgold = 0.3819660; // 0.5*(3 - sqrt(5))
        let mut ok: bool;
        let (mut etemp, mut fu, mut p, mut q, mut r): (f64, f64, f64, f64, f64) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut tol1, mut tol2, mut u, mut v, mut w, mut xm): (f64, f64, f64, f64, f64, f64) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut e: f64 = 0.0;
        let mut d: f64 = f64::MAX; // OCCT RealLast()

        self.a = if ax < cx { ax } else { cx };
        self.b = if ax > cx { ax } else { cx };
        self.x = bx;
        w = bx;
        v = bx;
        if !self.my_f {
            ok = match f.value(self.x) {
                Some(fv) => {
                    self.fx = fv;
                    true
                }
                None => false,
            };
            if !ok {
                return;
            }
        }
        self.fw = self.fx;
        self.fv = self.fx;
        self.iter = 1;
        while self.iter <= self.iter_max {
            xm = 0.5 * (self.a + self.b);
            tol1 = self.x_tol * self.x.abs() + self.eps_z;
            tol2 = 2.0 * tol1;
            if self.is_solution_reached() {
                self.done = true;
                return;
            }
            if e.abs() > tol1 {
                r = (self.x - w) * (self.fx - self.fv);
                q = (self.x - v) * (self.fx - self.fw);
                p = (self.x - v) * q - (self.x - w) * r;
                q = 2.0 * (q - r);
                if q > 0.0 {
                    p = -p;
                }
                q = q.abs();
                etemp = e;
                e = d;
                if p.abs() >= (0.5 * q * etemp).abs()
                    || p <= q * (self.a - self.x)
                    || p >= q * (self.b - self.x)
                {
                    e = if self.x >= xm { self.a - self.x } else { self.b - self.x };
                    d = cgold * e;
                } else {
                    d = p / q;
                    u = self.x + d;
                    if u - self.a < tol2 || self.b - u < tol2 {
                        d = tol1.copysign(xm - self.x);
                    }
                }
            } else {
                e = if self.x >= xm { self.a - self.x } else { self.b - self.x };
                d = cgold * e;
            }
            u = if d.abs() >= tol1 {
                self.x + d
            } else {
                self.x + tol1.copysign(d)
            };
            ok = match f.value(u) {
                Some(fu2) => {
                    fu = fu2;
                    true
                }
                None => false,
            };
            if !ok {
                return;
            }
            if fu <= self.fx {
                if u >= self.x {
                    self.a = self.x;
                } else {
                    self.b = self.x;
                }
                // SHFT(v, w, x, u); SHFT(fv, fw, fx, fu);
                v = w;
                w = self.x;
                self.x = u;
                self.fv = self.fw;
                self.fw = self.fx;
                self.fx = fu;
            } else {
                if u < self.x {
                    self.a = u;
                } else {
                    self.b = u;
                }
                if fu <= self.fw || w == self.x {
                    v = w;
                    w = u;
                    self.fv = self.fw;
                    self.fw = fu;
                } else if fu <= self.fv || v == self.x || v == w {
                    v = u;
                    self.fv = fu;
                }
            }
            self.iter += 1;
        }
        self.done = false;
    }

    /// OCCT IsSolutionReached (math_BrentMinimum.lxx L17-21).
    fn is_solution_reached(&self) -> bool {
        let two_tol = 2.0 * (self.x_tol * self.x.abs() + self.eps_z);
        (self.x <= two_tol + self.a) && (self.x >= self.b - two_tol)
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
    /// OCCT Location().
    pub fn location(&self) -> f64 {
        assert!(self.done, "BrentMinimum::Location on not-done");
        self.x
    }
    /// OCCT Minimum().
    pub fn minimum(&self) -> f64 {
        assert!(self.done, "BrentMinimum::Minimum on not-done");
        self.fx
    }
    /// OCCT NbIterations().
    pub fn nb_iterations(&self) -> i32 {
        assert!(self.done, "BrentMinimum::NbIterations on not-done");
        self.iter
    }
}

// =============================================================================
// math_BracketMinimum — 1:1 of math_BracketMinimum::Perform
// (math_BracketMinimum.cxx L62-211).
// =============================================================================

fn sign_macro(a: f64, b: f64) -> f64 {
    if b > 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

/// OCCT math_BracketMinimum::Perform (math_BracketMinimum.cxx L62-211).
/// Bracket a minimum of `f` starting from [a0, b0] with known fa/fb when
/// provided.  Returns Some((ax, bx, cx, fax, fbx, fcx)) on success.
pub fn bracket_minimum(
    f: &mut dyn FnMut(f64) -> Option<f64>,
    a0: f64,
    b0: f64,
    fa0: Option<f64>,
    fb0: Option<f64>,
) -> Option<(f64, f64, f64, f64, f64, f64)> {
    let gold = 1.618034;
    let glimit = 100.0;
    let tiny = 1.0e-20;
    let mut ax = a0;
    let mut bx = b0;
    let mut fax = 0.0;
    let mut fbx = 0.0;
    if let Some(fa) = fa0 {
        fax = fa;
    } else {
        fax = f(ax)?;
    }
    if let Some(fb) = fb0 {
        fbx = fb;
    } else {
        fbx = f(bx)?;
    }
    if fbx > fax {
        std::mem::swap(&mut ax, &mut bx);
        std::mem::swap(&mut fax, &mut fbx);
    }
    let lambda = gold;
    let mut cx = bx + lambda * (bx - ax);
    let mut fcx = f(cx)?;
    while fbx > fcx {
        let r = (bx - ax) * (fbx - fcx);
        let q = (bx - cx) * (fbx - fax);
        let mut u = bx - ((bx - cx) * q - (bx - ax) * r)
            / (2.0 * sign_macro((q - r).abs().max(tiny), q - r));
        let ulim = bx + glimit * (cx - bx);
        let mut fu;
        if (bx - u) * (u - cx) > 0.0 {
            fu = f(u)?;
            if fu < fcx {
                ax = bx;
                bx = u;
                fax = fbx;
                fbx = fu;
                return Some((ax, bx, cx, fax, fbx, fcx));
            } else if fu > fbx {
                cx = u;
                fcx = fu;
                return Some((ax, bx, cx, fax, fbx, fcx));
            }
            // Get the next probe after (B, C).
            u = cx + lambda * (cx - bx);
            fu = f(u)?;
        } else if (cx - u) * (u - ulim) > 0.0 {
            fu = f(u)?;
        } else if (u - ulim) * (ulim - cx) >= 0.0 {
            u = ulim;
            fu = f(u)?;
        } else {
            u = cx + gold * (cx - bx);
            fu = f(u)?;
        }
        ax = bx;
        bx = cx;
        cx = u;
        fax = fbx;
        fbx = fcx;
        fcx = fu;
    }
    Some((ax, bx, cx, fax, fbx, fcx))
}

// =============================================================================
// math_BFGS — 1:1 of math_BFGS::Perform (math_BFGS.cxx L110-454) with the
// DirFunction 1D line search (MinimizeDirection L203-321: ComputeInitScale +
// BracketMinimum + BrentMinimum).
// =============================================================================

/// OCCT math_BFGS::ComputeInitScale (math_BFGS.cxx L115-135).
fn bfgs_compute_init_scale(f0: f64, dir: &[f64], gr: &[f64], scale: &mut f64) -> bool {
    let dy1: f64 = gr.iter().zip(dir.iter()).map(|(g, d)| g * d).sum();
    if dy1.abs() < 1.0e-12 {
        return false;
    }
    let a_hnr1: f64 = dir.iter().map(|x| x * x).sum();
    let alfa = 0.7 * (-f0) / dy1;
    *scale = 0.015 / a_hnr1.sqrt();
    if *scale > alfa {
        *scale = alfa;
    }
    true
}

/// Adapter from a closure to the OCCT math_Function (FunctionValue) trait.
struct FnValueAdapter<'a>(&'a mut dyn FnMut(f64) -> Option<f64>);

impl<'a> crate::math::root::FunctionValue for FnValueAdapter<'a> {
    fn value(&mut self, x: f64) -> Option<f64> {
        (self.0)(x)
    }
}

/// OCCT math_BFGS::Perform (math_BFGS.cxx L327-443).  `f_val` is the
/// n-variable function returning (F, gradient).  The solution vector is
/// returned.
pub fn bfgs_minimize_occt(
    n: usize,
    starting: &[f64],
    tolerance: f64,
    itermax: i32,
    zeps: f64,
    f_val: &mut dyn FnMut(&[f64]) -> Option<(f64, Vec<f64>)>,
) -> Option<Vec<f64>> {
    let mut location = starting.to_vec();
    let (mut prev_min, mut grad) = f_val(&location)?;
    let x_tol = tolerance;
    let epsz = zeps;
    let mut hessin = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        hessin[i][i] = 1.0;
    }
    let mut xi = vec![0.0f64; n];
    for i in 0..n {
        xi[i] = -grad[i];
    }
    for _ in 0..itermax {
        let the_min0 = prev_min;
        // MinimizeDirection (L203-321).
        let mut scale = 0.0;
        if !bfgs_compute_init_scale(the_min0, &xi, &grad, &mut scale) {
            return None;
        }
        // DirFunction f1d(alpha) = F(P0 + alpha*Dir).
        let p0 = location.clone();
        let dir = xi.clone();
        let mut dir_fn = |alpha: f64| -> Option<f64> {
            let mut p = vec![0.0f64; n];
            for i in 0..n {
                p[i] = p0[i] + alpha * dir[i];
            }
            f_val(&p).map(|(f, _g)| f)
        };
        // math_BracketMinimum Bracket(0.0, lambda) with FA = F0 (L264-270).
        let (ax, xx, bx2, _fax, fxx, _fbx) =
            bracket_minimum(&mut dir_fn, 0.0, scale, Some(the_min0), None)?;
        // math_BrentMinimum Sol(tol, Fxx, niter, 1.e-08) (L281-282).
        let mut brent_f = dir_fn;
        let mut brent = BrentMinimum::new_with_fbx(1.0e-3, fxx, 100, 1.0e-8);
        let mut brent_fn = FnValueAdapter(&mut brent_f);
        brent.perform(&mut brent_fn, ax, xx, bx2);
        if !brent.is_done() {
            return None;
        }
        let loc_scale = brent.location();
        let min_val = brent.minimum();
        // P += Dir * Scale (L287-288).
        for i in 0..n {
            location[i] += xi[i] * loc_scale;
        }
        // OCCT IsSolutionReached (L449-454).
        if 2.0 * (min_val - prev_min).abs()
            <= x_tol * (min_val.abs() + prev_min.abs() + epsz)
        {
            return Some(location);
        }
        prev_min = min_val;
        // dg = grad_new - grad_old (L386-399).
        let (_f_new, grad_new) = f_val(&location)?;
        let mut dg = vec![0.0f64; n];
        for i in 0..n {
            dg[i] = grad_new[i] - grad[i];
        }
        let mut hdg = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..n {
                hdg[i] += hessin[i][j] * dg[j];
            }
        }
        let mut fac = 0.0;
        let mut fae = 0.0;
        for i in 0..n {
            fac += dg[i] * xi[i];
            fae += dg[i] * hdg[i];
        }
        if fac.abs() < 1.0e-300 || fae.abs() < 1.0e-300 {
            return Some(location);
        }
        fac = 1.0 / fac;
        let fad = 1.0 / fae;
        for i in 0..n {
            dg[i] = fac * xi[i] - fad * hdg[i];
        }
        for i in 0..n {
            for j in 0..n {
                hessin[i][j] += fac * xi[i] * xi[j] - fad * hdg[i] * hdg[j] + fae * dg[i] * dg[j];
            }
        }
        for i in 0..n {
            xi[i] = 0.0;
            for j in 0..n {
                xi[i] -= hessin[i][j] * grad_new[j];
            }
        }
        grad = grad_new;
    }
    Some(location)
}
