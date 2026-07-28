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
pub fn brent_minimize(f: impl Fn(f64) -> f64, a: f64, b: f64, tol: f64) -> f64 {
    let (mut lo, mut hi) = if a < b { (a, b) } else { (b, a) };
    let (mut x, mut w, mut v) = ((a + b) / 2.0, (a + b) / 2.0, (a + b) / 2.0);
    let (mut fx, mut fw, mut fv) = (f(x), f(x), f(x));
    let (mut d, mut e) = (0.0f64, 0.0f64);
    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        let tol1 = tol * x.abs() + 1e-12;
        let tol2 = 2.0 * tol1;
        if (x - mid).abs() <= tol2 - (hi - lo) / 2.0 {
            return x;
        }
        let mut use_para = false;
        let mut u = 0.0;
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let qq = (x - v) * (fx - fw);
            let p = (x - v) * qq - (x - w) * r;
            let q = 2.0 * (qq - r);
            if q.abs() > tol1 {
                u = x - p / q;
                if u > lo + tol1 && u < hi - tol1 && (u - x).abs() < e {
                    use_para = true;
                }
            }
        }
        if !use_para {
            u = if x >= mid { x - PHI * (x - lo) } else { x + PHI * (hi - x) };
            e = d;
            d = u - x;
        }
        let fu = f(u);
        if fu <= fx {
            if u >= x { lo = x; } else { hi = x; }
            v = w; fv = fw;
            w = x; fw = fx;
            x = u; fx = fu;
        } else {
            if u >= x { hi = u; } else { lo = u; }
            if fu <= fw || w == x {
                v = w; fv = fw;
                w = u; fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u; fv = fu;
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
    let mut rng = fastrand::Rng::new();
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
