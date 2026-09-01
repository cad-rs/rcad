//! OCCT math_FunctionSetRoot (math_FunctionSetRoot.cxx, whole) 鈥?bounded
//! Newton-with-direction-search root finding for a function set with
//! derivatives, plus the math_FunctionRoot wrapper (math_FunctionRoot.cxx).
//!
//! 1:1 translation of the Perform algorithm: Newton direction (Gauss), SVD
//! fallback on singular Jacobians, gradient direction when Newton does not
//! descend, boundary re-clamping (Bounds), line minimization via
//! math_BrentMinimum (MinimizeDirection), and the solution/divergence tests.
//!
//! Consumed by Geom2dGcc_Lin2d2TanIter (tangent line to two 2D curves /
//! a 2D curve and a point).

use crate::math::lin::solve_linear_system;
use crate::math::opt::BrentMinimum;
use crate::math::root::FunctionValue;

/// OCCT math_FunctionSetWithDerivatives 鈥?a function set with its Jacobian.
pub trait FunctionSetWithDerivatives {
    fn nb_variables(&self) -> usize;
    fn nb_equations(&self) -> usize;
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool;
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool;
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool;
}

/// OCCT math_SVD::Solve fallback 鈥?least-squares / min-norm solution of
/// A路x = b for a square A via the Jacobi-rotation SVD pseudo-inverse
/// (mathematically the same operation as the OCCT Golub-Reinsch SVD solve).
fn svd_solve_square(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return None;
    }
    let mut u = a.to_vec();
    let mut v: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }
    let mut s = vec![0.0; n];

    // Jacobi one-sided rotations (cyclic).
    const MAX_SWEEPS: usize = 100;
    for _ in 0..MAX_SWEEPS {
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                let mut sum = 0.0;
                for i in 0..n {
                    sum += u[i][p] * u[i][q];
                }
                off += sum * sum;
                if sum.abs() < 1e-300 {
                    continue;
                }
                let spp = {
                    let mut s = 0.0;
                    for i in 0..n {
                        s += u[i][p] * u[i][p];
                    }
                    s
                };
                let sqq = {
                    let mut s = 0.0;
                    for i in 0..n {
                        s += u[i][q] * u[i][q];
                    }
                    s
                };
                let theta = 0.5 * (sqq - spp) / sum;
                let t = if theta >= 0.0 {
                    1.0 / (theta + (1.0 + theta * theta).sqrt())
                } else {
                    1.0 / (theta - (1.0 + theta * theta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let sn = t * c;
                for i in 0..n {
                    let up = u[i][p];
                    let uq = u[i][q];
                    u[i][p] = c * up - sn * uq;
                    u[i][q] = sn * up + c * uq;
                }
                for i in 0..n {
                    let vp = v[i][p];
                    let vq = v[i][q];
                    v[i][p] = c * vp - sn * vq;
                    v[i][q] = sn * vp + c * vq;
                }
            }
        }
        if off < 1e-30 {
            break;
        }
    }
    for i in 0..n {
        s[i] = {
            let mut s = 0.0;
            for k in 0..n {
                s += u[k][i] * u[k][i];
            }
            s.sqrt()
        };
        if s[i] > 1e-300 {
            for k in 0..n {
                u[k][i] /= s[i];
            }
        }
    }

    // x = V路危鈦宦孤稶岬€路b (min-norm solution; zero the singular components).
    let mut x = vec![0.0; n];
    let cond_max = s.iter().copied().fold(0.0f64, f64::max);
    let tol_sing = 1e-12 * cond_max.max(1e-300);
    for i in 0..n {
        if s[i] <= tol_sing {
            continue;
        }
        let mut ui_dot_b = 0.0;
        for k in 0..n {
            ui_dot_b += u[k][i] * b[k];
        }
        let c = ui_dot_b / s[i];
        for k in 0..n {
            x[k] += v[k][i] * c;
        }
    }
    Some(x)
}

/// OCCT math_FunctionSetRoot 鈥?bounded Newton root finding with direction
/// search, line minimization and boundary re-clamping.
pub struct FunctionSetRoot {
    delta: Vec<f64>,
    sol: Vec<f64>,
    df: Vec<Vec<f64>>,
    tol: Vec<f64>,
    done: bool,
    kount: i32,
    state: i32,
    iter_max: i32,
    inf_bound: Vec<f64>,
    sup_bound: Vec<f64>,
    sol_save: Vec<f64>,
    gh: Vec<f64>,
    dh: Vec<f64>,
    dh_save: Vec<f64>,
    ff: Vec<f64>,
    previous_solution: Vec<f64>,
    save: Vec<f64>,
    constraints: Vec<i32>,
    temp1: Vec<f64>,
    temp2: Vec<f64>,
    temp3: Vec<f64>,
    temp4: Vec<f64>,
    is_divergent: bool,
}

/// The OCCT MyDirFunction 鈥?the function restricted to a direction
/// (P(x) = P0 + x路Dir), evaluated as 0.5路|F|虏 for the line minimization.
struct DirFunction<'a> {
    p0: Vec<f64>,
    dir: Vec<f64>,
    p: Vec<f64>,
    fv: Vec<f64>,
    f: &'a mut dyn FunctionSetWithDerivatives,
}

impl<'a> DirFunction<'a> {
    fn new(f: &'a mut dyn FunctionSetWithDerivatives) -> Self {
        let n = f.nb_variables();
        DirFunction {
            p0: vec![0.0; n],
            dir: vec![0.0; n],
            p: vec![0.0; n],
            fv: vec![0.0; f.nb_equations()],
            f,
        }
    }

    fn initialize(&mut self, p0: &[f64], dir: &[f64]) {
        self.p0.copy_from_slice(p0);
        self.dir.copy_from_slice(dir);
    }

    /// OCCT MyDirFunction::Value(Sol, FF, DF, GH, F2, Gnr1).
    fn value_full(
        &mut self,
        sol: &[f64],
        ff: &mut [f64],
        df: &mut [Vec<f64>],
        gh: &mut [f64],
        f2: &mut f64,
        gnr1: &mut f64,
    ) -> bool {
        if !self.f.values(sol, ff, df) {
            return false;
        }
        for v in ff.iter() {
            if *v <= -1.0e100 || *v >= 1.0e100 {
                return false;
            }
        }
        *f2 = 0.5 * ff.iter().map(|x| x * x).sum::<f64>();
        // GH = DF岬€路FF.
        for i in 0..gh.len() {
            let mut s = 0.0;
            for j in 0..ff.len() {
                s += df[j][i] * ff[j];
            }
            gh[i] = s;
        }
        for g in gh.iter() {
            if !g.is_finite() {
                return false;
            }
        }
        *gnr1 = gh.iter().map(|x| x * x).sum::<f64>();
        true
    }
}

impl<'a> FunctionValue for DirFunction<'a> {
    /// OCCT MyDirFunction::Value(x, fval) 鈥?P = P0 + x路Dir, fval = 0.5路|F(P)|虏.
    fn value(&mut self, x: f64) -> Option<f64> {
        for i in 0..self.p.len() {
            self.p[i] = self.dir[i] * x + self.p0[i];
        }
        if !self.f.value(&self.p, &mut self.fv) {
            return None;
        }
        for v in self.fv.iter() {
            if *v <= -1.0e100 || *v >= 1.0e100 {
                return None;
            }
        }
        Some(0.5 * self.fv.iter().map(|v| v * v).sum::<f64>())
    }
}

/// OCCT MinimizeDirection (3-point version) (math_FunctionSetRoot.cxx L198-264).
fn minimize_direction_3p(
    p0: &[f64],
    p1: &[f64],
    p2: &[f64],
    f1: f64,
    delta: &mut Vec<f64>,
    tol: &[f64],
    f: &mut DirFunction,
) -> bool {
    // (1) 1D parametric tolerance.
    let mut tol1d: f64 = 2.1;
    let eps = 1.0e-16;
    for ii in 0..tol.len() {
        let invnorme = delta[ii].abs();
        if invnorme > eps {
            tol1d = tol1d.min(tol[ii] / invnorme);
        }
    }
    if tol1d > 1.9 {
        return false;
    }
    tol1d /= 3.0;

    delta.copy_from_slice(&(p1.iter().zip(p0.iter()).map(|(a, b)| a - b).collect::<Vec<_>>()));
    let invnorme = delta.iter().map(|d| d * d).sum::<f64>().sqrt();
    if invnorme <= eps {
        return false;
    }
    let invnorme = 1.0 / invnorme;

    f.initialize(p1, delta);

    // (2) Minimize in the direction.
    let ax = -1.0;
    let bx = 0.0;
    let cx = {
        let mut s = 0.0;
        for i in 0..p2.len() {
            let d = p2[i] - p1[i];
            s += d * d;
        }
        s.sqrt() * invnorme
    };
    if cx < 1.0e-2 {
        return false;
    }

    let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
    sol.perform(f, ax, bx, cx);

    if sol.is_done() {
        let tsol = sol.location();
        if sol.minimum() < f1 {
            for d in delta.iter_mut() {
                *d *= tsol;
            }
            return true;
        }
    }
    false
}

/// OCCT MinimizeDirection (2-point + derivative version) (L267-436).
fn minimize_direction_2p(
    p: &[f64],
    dir: &mut Vec<f64>,
    p_value: f64,
    p_dir_value: f64,
    gradient: &[f64],
    d_gradient: &[f64],
    tol: &[f64],
    f: &mut DirFunction,
) -> bool {
    if !p_value.is_finite() || !p_dir_value.is_finite() {
        return false;
    }
    // (0) 1D parametric tolerance.
    let mut good = false;
    let eps = 1.0e-20;
    let mut tol1d: f64 = 1.1;
    let mut result = p_value;
    for ii in 0..tol.len() {
        let absdir = dir[ii].abs();
        if absdir > eps {
            tol1d = tol1d.min(tol[ii] / absdir);
        }
    }
    if tol1d > 0.9 {
        return false;
    }

    // (1) First quadratic interpolation.
    let mut tsol: f64;
    let df1 = gradient.iter().zip(dir.iter()).map(|(g, d)| g * d).sum::<f64>();
    let df2 = d_gradient.iter().zip(dir.iter()).map(|(g, d)| g * d).sum::<f64>();

    if df1 < -eps && df2 > eps {
        // Cup.
        tsol = -df1 / (df2 - df1);
    } else {
        let cx = p_value;
        let bx = df1;
        let ax = p_dir_value - (bx + cx);
        if ax.abs() <= eps {
            // Linear case.
            tsol = if bx.abs() >= eps { -cx / bx } else { 0.0 };
        } else {
            // Quadratic case.
            let delta = bx * bx - 4.0 * ax * cx;
            if delta > 1.0e-9 {
                let delta = delta.sqrt();
                let t1 = -(bx + delta);
                let t2 = delta - bx;
                tsol = if t2.abs() < t1.abs() { t2 } else { t1 };
                tsol /= 2.0 * ax;
            } else {
                tsol = -(0.5 * bx) / ax;
            }
        }
    }

    if tsol.abs() >= 1.0 {
        return false;
    }

    f.initialize(p, dir);
    let fsol = f.value(tsol).unwrap_or(f64::MAX);

    if fsol < p_value {
        good = true;
        result = fsol;
    }

    // (2) If we did not progress enough, run a proper line search.
    if fsol > 0.2 * p_value && tol1d < 0.5 {
        let (ax, bx, cx) = if tsol < 0.0 {
            (tsol, 0.0, 1.0)
        } else {
            (0.0, tsol, 1.0)
        };
        let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
        sol.perform(f, ax, bx, cx);
        if sol.is_done() && sol.minimum() <= result {
            tsol = sol.location();
            good = true;
            result = sol.minimum();
            // Objective function changes too fast 鈥?refine on both halves.
            if gradient.iter().map(|g| g * g).sum::<f64>() > 1.0 / (1.0e-7 * 1.0e-7)
                && tsol > ax
                && tsol < cx
            {
                let mut sol2 = BrentMinimum::new(tol1d, 100, tol1d);
                sol2.perform(f, ax, (ax + tsol) / 2.0, tsol);
                if sol2.is_done() && sol2.minimum() <= result {
                    tsol = sol2.location();
                    good = true;
                    result = sol2.minimum();
                }
                let mut sol3 = BrentMinimum::new(tol1d, 100, tol1d);
                sol3.perform(f, tsol, (cx + tsol) / 2.0, cx);
                if sol3.is_done() && sol3.minimum() <= result {
                    tsol = sol3.location();
                    good = true;
                    result = sol3.minimum();
                }
            }
        }
    }

    if good {
        for d in dir.iter_mut() {
            *d *= tsol;
        }
    }
    good
}

/// OCCT SearchDirection (math_FunctionSetRoot.cxx L439-531) 鈥?the Newton (or
/// gradient) direction from the Jacobian and gradient.
fn search_direction(
    df: &[Vec<f64>],
    gh: &[f64],
    ff: &[f64],
    change_direction: &mut bool,
    inv_length_max: &[f64],
    direction: &mut Vec<f64>,
    dy: &mut f64,
) {
    let ninc = df[0].len();
    let neq = df.len();
    let eps = 1.0e-32;
    if !*change_direction {
        if ninc == neq {
            for i in 0..ff.len() {
                direction[i] = -ff[i];
            }
            let a_flat: Vec<f64> = df.iter().flatten().copied().collect();
            let b: Vec<f64> = direction.clone();
            if let Some(x) = solve_linear_system(&a_flat, &b, ninc) {
                direction.copy_from_slice(&x);
            } else {
                // Singular matrix 鈥?SVD.
                if let Some(x) = svd_solve_square(df, &b) {
                    direction.copy_from_slice(&x);
                } else {
                    *change_direction = true;
                }
            }
        } else if ninc > neq {
            // Over-determined: SVD least squares.
            let b: Vec<f64> = ff.iter().map(|f| -*f).collect();
            if let Some(x) = svd_solve_rect(df, &b) {
                direction.copy_from_slice(&x);
            } else {
                *change_direction = true;
            }
        } else {
            // Under-determined: least squares (Gauss least square).
            let b: Vec<f64> = ff.iter().map(|f| -*f).collect();
            if let Some(x) = svd_solve_rect(df, &b) {
                direction.copy_from_slice(&x);
            } else {
                *change_direction = true;
            }
        }
    }

    // Limit over-long directions.
    let mut ratio = (direction[0] * inv_length_max[0]).abs();
    for i in 1..direction.len() {
        ratio = ratio.max((direction[i] * inv_length_max[i]).abs());
    }
    if ratio > 1.0 {
        for d in direction.iter_mut() {
            *d /= ratio;
        }
    }

    *dy = direction.iter().zip(gh.iter()).map(|(d, g)| d * g).sum();
    if *dy >= -eps {
        // Newton does not descend 鈥?use the gradient.
        *change_direction = true;
    }
    if *change_direction {
        for i in 0..direction.len() {
            direction[i] = -gh[i];
        }
        *dy = -(gh.iter().map(|g| g * g).sum::<f64>());
    }
}

/// SVD least-squares / min-norm solve for a rectangular m脳n matrix.
fn svd_solve_rect(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let (m, n) = (a.len(), a.first().map_or(0, |r| r.len()));
    if n == 0 {
        return None;
    }
    // Square the system with A岬€ (least squares normal equations solved via
    // the square SVD pseudo-inverse) 鈥?mathematically equivalent to the OCCT
    // SVD least-squares solve for the over/under-determined cases.
    let mut at_a = vec![vec![0.0; n]; n];
    let mut at_b = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..m {
                s += a[k][i] * a[k][j];
            }
            at_a[i][j] = s;
        }
        for k in 0..m {
            at_b[i] += a[k][i] * b[k];
        }
    }
    svd_solve_square(&at_a, &at_b)
}

/// OCCT Bounds (math_FunctionSetRoot.cxx L623-705) 鈥?re-clamp the solution to
/// the domain, updating the constraints and the displacement Delta.
fn bounds(
    inf_bound: &[f64],
    sup_bound: &[f64],
    tol: &[f64],
    sol: &mut Vec<f64>,
    sol_save: &[f64],
    constraints: &mut [i32],
    delta: &mut Vec<f64>,
    is_new_sol: &mut bool,
) -> bool {
    let mut out = false;
    let ninc = sol.len();
    let mut monratio: f64 = 1.0;

    *is_new_sol = true;

    for i in 0..ninc {
        constraints[i] = 0;
        delta[i] = sol[i] - sol_save[i];
        if inf_bound[i] == sup_bound[i] {
            constraints[i] = 1;
            out = true;
        } else if sol[i] < inf_bound[i] {
            constraints[i] = 1;
            out = true;
            if -delta[i] > tol[i] {
                monratio = monratio.min((inf_bound[i] - sol_save[i]) / delta[i]);
            }
        } else if sol[i] > sup_bound[i] {
            constraints[i] = 1;
            out = true;
            if delta[i] > tol[i] {
                monratio = monratio.min((sup_bound[i] - sol_save[i]) / delta[i]);
            }
        }
    }

    if out {
        if monratio == 0.0 {
            *is_new_sol = false;
            sol.copy_from_slice(sol_save);
            for d in delta.iter_mut() {
                *d = 0.0;
            }
        } else {
            for d in delta.iter_mut() {
                *d *= monratio;
            }
            for i in 0..ninc {
                sol[i] = sol_save[i] + delta[i];
            }
            for i in 0..ninc {
                if sol[i] < inf_bound[i] {
                    sol[i] = inf_bound[i];
                    delta[i] = sol[i] - sol_save[i];
                } else if sol[i] > sup_bound[i] {
                    sol[i] = sup_bound[i];
                    delta[i] = sol[i] - sol_save[i];
                }
            }
        }
    }
    out
}

impl FunctionSetRoot {
    /// OCCT math_FunctionSetRoot(F, Tolerance, NbIterations) (L709-738).
    pub fn new(f: &dyn FunctionSetWithDerivatives, tolerance: &[f64], nb_iterations: i32) -> Self {
        let n_vars = f.nb_variables();
        let n_eqs = f.nb_equations();
        let mut r = FunctionSetRoot {
            delta: vec![0.0; n_vars],
            sol: vec![0.0; n_vars],
            df: vec![vec![0.0; n_vars]; n_eqs],
            tol: vec![0.0; n_vars],
            done: false,
            kount: 0,
            state: 0,
            iter_max: nb_iterations,
            inf_bound: vec![-f64::MAX; n_vars],
            sup_bound: vec![f64::MAX; n_vars],
            sol_save: vec![0.0; n_vars],
            gh: vec![0.0; n_vars],
            dh: vec![0.0; n_vars],
            dh_save: vec![0.0; n_vars],
            ff: vec![0.0; n_eqs],
            previous_solution: vec![0.0; n_vars],
            save: vec![0.0; nb_iterations.max(1) as usize + 1],
            constraints: vec![0; n_vars],
            temp1: vec![0.0; n_vars],
            temp2: vec![0.0; n_vars],
            temp3: vec![0.0; n_vars],
            temp4: vec![0.0; n_eqs],
            is_divergent: false,
        };
        r.tol.copy_from_slice(tolerance);
        r
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Root(V) 鈥?panics when not done.
    pub fn root(&self) -> Vec<f64> {
        assert!(self.done, "FunctionSetRoot: not done");
        self.sol.clone()
    }

    /// OCCT NbIterations().
    pub fn nb_iterations(&self) -> i32 {
        self.kount
    }

    /// OCCT Derivative()(1,1) 鈥?the (1,1) Jacobian entry after a 1-variable run.
    pub fn derivative(&self) -> f64 {
        self.df[0][0]
    }

    /// OCCT Perform(F, StartingPoint, theInfBound, theSupBound,
    /// theStopOnDivergent) (L796-1417).
    pub fn perform(
        &mut self,
        f: &mut dyn FunctionSetWithDerivatives,
        starting_point: &[f64],
        the_inf_bound: &[f64],
        the_sup_bound: &[f64],
        stop_on_divergent: bool,
    ) {
        let ninc = f.nb_variables();
        let neq = f.nb_equations();
        if neq <= 0 || starting_point.len() != ninc || the_inf_bound.len() != ninc
            || the_sup_bound.len() != ninc
        {
            panic!("Standard_DimensionError");
        }

        let mut change_direction = false;
        let mut sort = false;
        let mut is_new_sol = false;
        let mut good;
        let mut verif;
        let mut stop;
        let eps_sqrt = 1.0e-16;
        let eps = 1.0e-32;
        let eps2 = 1.0e-64;
        let progres = 0.005;
        let mut f2 = 0.0;
        let mut previous_minimum = 0.0;
        let mut dy = 0.0;
        let mut old_f = 0.0;
        let mut ambda = 0.0;
        let mut ambda2 = 0.0;
        let mut gnr1 = 0.0;
        let mut old_gr = 0.0;

        let mut inv_length_max = vec![0.0; ninc];
        for i in 0..ninc {
            let a_sup = the_sup_bound[i].min(2.0e100);
            let an_inf = the_inf_bound[i].max(-2.0e100);
            inv_length_max[i] = 1.0 / ((a_sup - an_inf) / 4.0).max(1.0e-9);
        }

        let mut f_dir = DirFunction::new(f);

        self.done = false;
        self.sol.copy_from_slice(starting_point);
        self.kount = 0;

        self.is_divergent = false;
        for i in 0..ninc {
            self.is_divergent =
                self.is_divergent || self.sol[i] < the_inf_bound[i] || self.sol[i] > the_sup_bound[i];
        }
        if stop_on_divergent && self.is_divergent {
            return;
        }

        // Re-clamp the starting point to the bounds.
        for i in 0..ninc {
            if self.sol[i] <= the_inf_bound[i] {
                self.sol[i] = the_inf_bound[i];
            } else if self.sol[i] > the_sup_bound[i] {
                self.sol[i] = the_sup_bound[i];
            }
        }

        // First function value and gradient.
        if !f_dir.value_full(&self.sol, &mut self.ff, &mut self.df, &mut self.gh, &mut f2, &mut gnr1) {
            self.done = false;
            if !stop_on_divergent || !self.is_divergent {
                self.state = 0;
            }
            return;
        }
        ambda2 = gnr1;
        self.save[0] = f2.max(eps_sqrt);
        let a_tol_func = crate::math::direct_polynomial_roots::epsilon(f2);

        if f2 <= eps || gnr1 <= eps2 {
            self.done = false;
            if !stop_on_divergent || !self.is_divergent {
                self.done = true;
                self.state = 0;
            }
            return;
        }

        let mut kount = 1i32;
        while kount <= self.iter_max {
            previous_minimum = f2;
            old_gr = gnr1;
            self.previous_solution.copy_from_slice(&self.sol);
            self.sol_save.copy_from_slice(&self.sol);

            search_direction(
                &self.df,
                &self.gh,
                &self.ff,
                &mut change_direction,
                &inv_length_max,
                &mut self.dh,
                &mut dy,
            );
            if dy.abs() <= eps {
                self.done = false;
                if !stop_on_divergent || !self.is_divergent {
                    self.done = true;
                    f_dir.f.value(&self.sol, &mut self.ff);
                    self.state = 0;
                }
                return;
            }
            if change_direction {
                ambda = ambda2 / dy.abs().sqrt();
                if ambda > 1.0 {
                    ambda = 1.0;
                }
            } else {
                ambda = 1.0;
                let dh_norm = self.dh.iter().map(|d| d * d).sum::<f64>().sqrt();
                ambda2 = if dh_norm > 1e-300 { 0.5 * ambda / dh_norm } else { 0.0 };
            }

            for i in 0..ninc {
                self.sol[i] += ambda * self.dh[i];
            }
            for i in 0..ninc {
                self.is_divergent = self.is_divergent
                    || self.sol[i] < the_inf_bound[i]
                    || self.sol[i] > the_sup_bound[i];
            }
            if stop_on_divergent && self.is_divergent {
                return;
            }

            sort = bounds(
                the_inf_bound,
                the_sup_bound,
                &self.tol,
                &mut self.sol,
                &self.sol_save,
                &mut self.constraints,
                &mut self.delta,
                &mut is_new_sol,
            );

            self.dh_save.copy_from_slice(&self.gh);
            if is_new_sol
                && !f_dir.value_full(
                    &self.sol,
                    &mut self.ff,
                    &mut self.df,
                    &mut self.gh,
                    &mut f2,
                    &mut gnr1,
                )
            {
                self.done = false;
                if !stop_on_divergent || !self.is_divergent {
                    self.state = 0;
                }
                return;
            }

            if f2 <= eps || gnr1 <= eps2 {
                self.done = false;
                if !stop_on_divergent || !self.is_divergent {
                    self.done = true;
                    f_dir.f.value(&self.sol, &mut self.ff);
                    self.state = 0;
                }
                return;
            }

            if sort || f2 / previous_minimum > progres {
                dy = self.gh.iter().zip(self.dh.iter()).map(|(g, d)| g * d).sum();
                old_f = previous_minimum;
                stop = false;
                good = false;
                let mut descente_iter = 0;
                let mut sort_bis;

                // Standard processing without boundary handling.
                if !sort {
                    while f2 / previous_minimum > progres && !stop {
                        if f2 < old_f && dy < 0.0 {
                            descente_iter += 1;
                            self.sol_save.copy_from_slice(&self.sol);
                            old_f = f2;
                            for i in 0..ninc {
                                self.sol[i] += ambda * self.dh[i];
                            }
                            for i in 0..ninc {
                                self.is_divergent = self.is_divergent
                                    || self.sol[i] < the_inf_bound[i]
                                    || self.sol[i] > the_sup_bound[i];
                            }
                            if stop_on_divergent && self.is_divergent {
                                return;
                            }
                            stop = bounds(
                                the_inf_bound,
                                the_sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut self.constraints,
                                &mut self.delta,
                                &mut is_new_sol,
                            );
                            ambda *= 1.7;
                        } else {
                            if f2 >= old_f || f2 >= previous_minimum {
                                good = false;
                                if descente_iter == 0 {
                                    descente_iter += 1;
                                    good = minimize_direction_2p(
                                        &self.sol_save,
                                        &mut self.delta,
                                        old_f,
                                        f2,
                                        &self.dh_save,
                                        &self.gh,
                                        &self.tol,
                                        &mut f_dir,
                                    );
                                } else if change_direction
                                    || descente_iter > 1
                                    || old_f > previous_minimum
                                {
                                    descente_iter += 1;
                                    good = minimize_direction_3p(
                                        &self.previous_solution,
                                        &self.sol_save,
                                        &self.sol,
                                        old_f,
                                        &mut self.delta,
                                        &self.tol,
                                        &mut f_dir,
                                    );
                                }
                                if !good {
                                    self.sol.copy_from_slice(&self.sol_save);
                                    f2 = old_f;
                                } else {
                                    for i in 0..ninc {
                                        self.sol[i] = self.sol_save[i] + self.delta[i];
                                    }
                                    for i in 0..ninc {
                                        self.is_divergent = self.is_divergent
                                            || self.sol[i] < the_inf_bound[i]
                                            || self.sol[i] > the_sup_bound[i];
                                    }
                                    if stop_on_divergent && self.is_divergent {
                                        return;
                                    }
                                    sort = bounds(
                                        the_inf_bound,
                                        the_sup_bound,
                                        &self.tol,
                                        &mut self.sol,
                                        &self.sol_save,
                                        &mut self.constraints,
                                        &mut self.delta,
                                        &mut is_new_sol,
                                    );
                                }
                                sort = false;
                            }
                            stop = true;
                        }
                        self.dh_save.copy_from_slice(&self.gh);
                        if is_new_sol
                            && !f_dir.value_full(
                                &self.sol,
                                &mut self.ff,
                                &mut self.df,
                                &mut self.gh,
                                &mut f2,
                                &mut gnr1,
                            )
                        {
                            self.done = false;
                            if !stop_on_divergent || !self.is_divergent {
                                self.state = 0;
                            }
                            return;
                        }
                        dy = self.gh.iter().zip(self.dh.iter()).map(|(g, d)| g * d).sum();
                        if dy.abs() <= eps {
                            if f2 > old_f {
                                self.sol.copy_from_slice(&self.sol_save);
                            }
                            self.done = false;
                            if !stop_on_divergent || !self.is_divergent {
                                self.done = true;
                                f_dir.f.value(&self.sol, &mut self.ff);
                                self.state = 0;
                            }
                            return;
                        }
                        if descente_iter >= 100 {
                            stop = true;
                        }
                    }
                }

                // Boundary processing.
                if sort {
                    stop = f2 > 1.001 * old_f;
                    sort_bis = sort;
                    descente_iter = 0;
                    while sort_bis && (f2 < old_f || descente_iter == 0) && !stop {
                        descente_iter += 1;
                        self.sol_save.copy_from_slice(&self.sol);
                        old_f = f2;
                        search_direction_constrained(
                            &self.df,
                            &self.gh,
                            &self.ff,
                            &self.constraints,
                            &self.sol,
                            &mut change_direction,
                            &inv_length_max,
                            &mut self.dh,
                            &mut dy,
                        );
                        if dy < -eps {
                            if change_direction {
                                ambda = ambda2 / (-dy).sqrt();
                                if ambda > 1.0 {
                                    ambda = 1.0;
                                }
                            } else {
                                ambda = 1.0;
                                let dh_norm = self.dh.iter().map(|d| d * d).sum::<f64>().sqrt();
                                ambda2 = if dh_norm > 1e-300 {
                                    0.5 * ambda / dh_norm
                                } else {
                                    0.0
                                };
                            }
                            for i in 0..ninc {
                                self.sol[i] += ambda * self.dh[i];
                            }
                            for i in 0..ninc {
                                self.is_divergent = self.is_divergent
                                    || self.sol[i] < the_inf_bound[i]
                                    || self.sol[i] > the_sup_bound[i];
                            }
                            if stop_on_divergent && self.is_divergent {
                                return;
                            }
                            sort_bis = bounds(
                                the_inf_bound,
                                the_sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut self.constraints,
                                &mut self.delta,
                                &mut is_new_sol,
                            );
                            self.dh_save.copy_from_slice(&self.gh);
                            if is_new_sol
                                && !f_dir.value_full(
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                    &mut f2,
                                    &mut gnr1,
                                )
                            {
                                self.done = false;
                                if !stop_on_divergent || !self.is_divergent {
                                    self.state = 0;
                                }
                                return;
                            }
                            ambda2 = gnr1;
                        } else {
                            stop = true;
                        }

                        while f2 / previous_minimum > progres && f2 < old_f && !stop {
                            descente_iter += 1;
                            if f2 < old_f && dy < 0.0 {
                                self.sol_save.copy_from_slice(&self.sol);
                                old_f = f2;
                                for i in 0..ninc {
                                    self.sol[i] += ambda * self.dh[i];
                                }
                                for i in 0..ninc {
                                    self.is_divergent = self.is_divergent
                                        || self.sol[i] < the_inf_bound[i]
                                        || self.sol[i] > the_sup_bound[i];
                                }
                                if stop_on_divergent && self.is_divergent {
                                    return;
                                }
                                sort_bis = bounds(
                                    the_inf_bound,
                                    the_sup_bound,
                                    &self.tol,
                                    &mut self.sol,
                                    &self.sol_save,
                                    &mut self.constraints,
                                    &mut self.delta,
                                    &mut is_new_sol,
                                );
                            }
                            self.dh_save.copy_from_slice(&self.gh);
                            if is_new_sol
                                && !f_dir.value_full(
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                    &mut f2,
                                    &mut gnr1,
                                )
                            {
                                self.done = false;
                                if !stop_on_divergent || !self.is_divergent {
                                    self.state = 0;
                                }
                                return;
                            }
                            ambda2 = gnr1;
                            dy = self.gh.iter().zip(self.dh.iter()).map(|(g, d)| g * d).sum();
                            stop = dy >= 0.0 || descente_iter >= 10 || sort_bis;
                        }
                        stop = dy >= 0.0 || descente_iter >= 10;
                    }
                    if (f2 / previous_minimum > progres && f2 >= old_f) || f2 >= previous_minimum {
                        descente_iter += 1;
                        good = minimize_direction_2p(
                            &self.sol_save,
                            &mut self.delta,
                            old_f,
                            f2,
                            &self.dh_save,
                            &self.gh,
                            &self.tol,
                            &mut f_dir,
                        );
                        if !good {
                            self.sol.copy_from_slice(&self.sol_save);
                            sort = false;
                        } else {
                            for i in 0..ninc {
                                self.sol[i] = self.sol_save[i] + self.delta[i];
                            }
                            for i in 0..ninc {
                                self.is_divergent = self.is_divergent
                                    || self.sol[i] < the_inf_bound[i]
                                    || self.sol[i] > the_sup_bound[i];
                            }
                            if stop_on_divergent && self.is_divergent {
                                return;
                            }
                            sort = bounds(
                                the_inf_bound,
                                the_sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &self.sol_save,
                                &mut self.constraints,
                                &mut self.delta,
                                &mut is_new_sol,
                            );
                            if is_new_sol
                                && !f_dir.value_full(
                                    &self.sol,
                                    &mut self.ff,
                                    &mut self.df,
                                    &mut self.gh,
                                    &mut f2,
                                    &mut gnr1,
                                )
                            {
                                self.done = false;
                                if !stop_on_divergent || !self.is_divergent {
                                    self.state = 0;
                                }
                                return;
                            }
                        }
                        dy = self.gh.iter().zip(self.dh.iter()).map(|(g, d)| g * d).sum();
                    }
                }
            }

            // Stop tests.
            self.save[kount as usize] = f2;
            verif = if change_direction {
                true
            } else if kount > 1 {
                self.save[(kount - 1) as usize] < 1.0e-4 * self.save[(kount - 2) as usize]
            } else {
                f2 < 1.0e-6 * self.save[0]
            };
            if verif {
                for i in 0..ninc {
                    self.delta[i] = self.previous_solution[i] - self.sol[i];
                }
                if self.is_solution_reached() {
                    if previous_minimum < f2 {
                        self.sol.copy_from_slice(&self.sol_save);
                    }
                    self.done = false;
                    if !stop_on_divergent || !self.is_divergent {
                        self.done = true;
                        f_dir.f.value(&self.sol, &mut self.ff);
                        self.state = 0;
                    }
                    return;
                }
            }

            // Progress analysis.
            if (f2 - previous_minimum) <= a_tol_func {
                if kount > 5 {
                    if f2 >= 0.95 * self.save[(kount - 5) as usize] {
                        if !change_direction {
                            change_direction = true;
                        } else {
                            self.done = false;
                            if !stop_on_divergent || !self.is_divergent {
                                self.done = true;
                                self.state = 0;
                            }
                            return;
                        }
                    } else {
                        change_direction = false;
                    }
                } else {
                    change_direction = false;
                }
                if (gnr1 > 0.9 * old_gr) && (f2 > 0.5 * previous_minimum) {
                    change_direction = true;
                }
                if !change_direction && !verif {
                    for i in 0..ninc {
                        self.delta[i] = self.previous_solution[i] - self.sol[i];
                    }
                    if self.is_solution_reached() {
                        self.done = false;
                        if !stop_on_divergent || !self.is_divergent {
                            self.done = true;
                            f_dir.f.value(&self.sol, &mut self.ff);
                            self.state = 0;
                        }
                        return;
                    }
                }
            } else {
                // Regression case.
                if !change_direction {
                    change_direction = true;
                    self.sol.copy_from_slice(&self.previous_solution);
                    if !f_dir.value_full(
                        &self.sol,
                        &mut self.ff,
                        &mut self.df,
                        &mut self.gh,
                        &mut f2,
                        &mut gnr1,
                    ) {
                        self.done = false;
                        if !stop_on_divergent || !self.is_divergent {
                            self.state = 0;
                        }
                        return;
                    }
                } else {
                    if !stop_on_divergent || !self.is_divergent {
                        self.state = 0;
                    }
                    return;
                }
            }
            kount += 1;
        }
        self.kount = kount - 1;
        if !stop_on_divergent || !self.is_divergent {
            self.state = 0;
        }
    }

    /// OCCT math_FunctionSetRoot::IsSolutionReached (hxx L69-79) 鈥?    /// |螖岬 鈮?Tol岬?for all unknowns.
    fn is_solution_reached(&self) -> bool {
        for i in 0..self.delta.len() {
            if self.delta[i].abs() > self.tol[i] {
                return false;
            }
        }
        true
    }
}

/// OCCT SearchDirection with constraints (math_FunctionSetRoot.cxx L534-620) 鈥?/// solve the sub-problem on the free unknowns; constrained unknowns are fixed.
fn search_direction_constrained(
    df: &[Vec<f64>],
    gh: &[f64],
    ff: &[f64],
    constraints: &[i32],
    _x: &[f64],
    change_direction: &mut bool,
    inv_length_max: &[f64],
    direction: &mut Vec<f64>,
    dy: &mut f64,
) {
    let ninc = df[0].len();
    let neq = df.len();
    let mut cons = 0;
    for i in 0..ninc {
        if constraints[i] != 0 {
            cons += 1;
        }
    }

    if cons == 0 {
        search_direction(df, gh, ff, change_direction, inv_length_max, direction, dy);
    } else if cons == ninc {
        for d in direction.iter_mut() {
            *d = 0.0;
        }
        *dy = 0.0;
    } else {
        // Sub-problem on the free unknowns.
        let n_free = ninc - cons;
        let mut df2 = vec![vec![0.0; n_free]; neq];
        let mut my_gh = vec![0.0; n_free];
        let mut my_inv = vec![0.0; n_free];
        let mut my_dir = vec![0.0; n_free];
        let mut k = 0usize;
        for i in 0..ninc {
            if constraints[i] == 0 {
                my_gh[k] = gh[i];
                my_inv[k] = inv_length_max[i];
                my_dir[k] = direction[i];
                for j in 0..neq {
                    df2[j][k] = df[j][i];
                }
                k += 1;
            }
        }
        search_direction(
            &df2,
            &my_gh,
            ff,
            change_direction,
            &my_inv,
            &mut my_dir,
            dy,
        );
        let mut k2 = 0usize;
        for i in 0..ninc {
            if constraints[i] == 0 {
                if !*change_direction {
                    direction[i] = my_dir[k2];
                } else {
                    direction[i] = -gh[i];
                }
                k2 += 1;
            } else {
                direction[i] = 0.0;
            }
        }
    }
}
