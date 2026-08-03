// OCCT math_FunctionSetRoot (math_FunctionSetRoot.cxx L796-1100, L439-705)
// Rust translation, specialized to the 2-variable / 1-equation case used by
// IntStart_SearchInside and IntWalk_IWalking (the algebraic function
// F(u,v) = Q(P(u,v))).
//
// The OCCT class is fully generic (N variables, M equations); here Ninc = 2,
// Neq = 1.  The SearchDirection SVD branch (Ninc > Neq) collapses to the
// minimum-norm Gauss-Newton step  Direction = -FF * Gradient / |Gradient|^2,
// which is exactly what OCCT's math_SVD::Solve computes for the under-determined
// system DF . Direction = -FF.

use rcad_kernel::math::opt::BrentMinimum;
use rcad_kernel::math::root::FunctionValue;

use super::surf_function::SurfFunction;

const EPS: f64 = 1e-32;
const EPS2: f64 = 1e-64;
const EPS_SQRT: f64 = 1e-16;
const PROGRES: f64 = 0.005;

/// A 1-D restriction of the function along a direction, used by the line
/// search / minimization (OCCT MyDirFunction, math_FunctionSetRoot.cxx L70-195).
/// `f` is a raw pointer, mirroring the OCCT `void* F` (the function is owned
/// by the caller, alive for the whole Perform).
struct DirFunction {
    p0: [f64; 2],
    dir: [f64; 2],
    p: [f64; 2],
    fv: [f64; 1],
    f: *mut SurfFunction,
}

impl DirFunction {
    fn new(f: &mut SurfFunction) -> Self {
        DirFunction {
            p0: [0.0; 2],
            dir: [0.0; 2],
            p: [0.0; 2],
            fv: [0.0],
            f: f as *mut SurfFunction,
        }
    }

    fn initialize(&mut self, p0: [f64; 2], dir: [f64; 2]) {
        self.p0 = p0;
        self.dir = dir;
    }

    /// OCCT MyDirFunction::Value(Sol, FF, DF, GH, F2, Gnr1) (L152-195).
    fn value_vec(
        &mut self,
        sol: [f64; 2],
        ff: &mut [f64; 1],
        df: &mut [[f64; 2]; 1],
        gh: &mut [f64; 2],
        f2: &mut f64,
        gnr1: &mut f64,
    ) -> bool {
        let func = unsafe { &mut *self.f };
        let Some((val, d)) = func.values(&sol) else {
            return false;
        };
        if val.abs() >= 1e+100 {
            return false;
        }
        ff[0] = val;
        df[0][0] = d[0];
        df[0][1] = d[1];
        *f2 = 0.5 * (ff[0] * ff[0]);
        // GH = DF^T * FF.
        gh[0] = df[0][0] * ff[0];
        gh[1] = df[0][1] * ff[0];
        *gnr1 = gh[0] * gh[0] + gh[1] * gh[1];
        true
    }
}

impl FunctionValue for DirFunction {
    /// OCCT MyDirFunction::Value(x, fval) (L120-150): F along the direction.
    fn value(&mut self, x: f64) -> Option<f64> {
        for i in 0..2 {
            self.p[i] = self.dir[i] * x + self.p0[i];
        }
        let func = unsafe { &mut *self.f };
        let Some(val) = func.value(&self.p) else {
            return None;
        };
        if val.abs() >= 1e+100 {
            return None;
        }
        Some(0.5 * (val * val))
    }
}

/// OCCT MinimizeDirection (math_FunctionSetRoot.cxx L198-264) — minimization
/// from three points P0, P1, P2.  `delta` is updated to `tsol * (P1 - P0)`.
fn minimize_direction_3(
    p0: &[f64; 2],
    p1: &[f64; 2],
    p2: &[f64; 2],
    f1: f64,
    delta: &mut [f64; 2],
    tol: &[f64; 2],
    f: &mut DirFunction,
) -> bool {
    // (1) Evaluation d'une tolerance parametrique 1D.
    let mut tol1d = 2.1f64;
    let eps = 1e-16;
    for ii in 0..tol.len() {
        let invnorme = delta[ii].abs();
        if invnorme > eps {
            tol1d = tol1d.min(tol[ii] / invnorme);
        }
    }
    if tol1d > 1.9 {
        return false; // Pas la peine de se fatiguer
    }
    tol1d /= 3.0;

    // Delta = P1 - P0.
    delta[0] = p1[0] - p0[0];
    delta[1] = p1[1] - p0[1];
    let mut invnorme = (delta[0] * delta[0] + delta[1] * delta[1]).sqrt();
    if invnorme <= eps {
        return false;
    }
    invnorme = 1.0 / invnorme;

    f.initialize(*p1, *delta);

    // (2) On minimise.
    let ax = -1.0;
    let bx = 0.0;
    let cx = ((p2[0] - p1[0]).powi(2) + (p2[1] - p1[1]).powi(2)).sqrt() * invnorme;
    if cx < 1e-2 {
        return false;
    }

    let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
    sol.perform(f, ax, bx, cx);

    if sol.is_done() {
        let tsol = sol.location();
        if sol.minimum() < f1 {
            delta[0] *= tsol;
            delta[1] *= tsol;
            return true;
        }
    }
    false
}

/// OCCT MinimizeDirection (math_FunctionSetRoot.cxx L266-436) — minimization
/// from two points and a derivative.  `dir` is updated to `tsol * dir`.
#[allow(clippy::too_many_arguments)]
fn minimize_direction_2(
    p: &[f64; 2],
    dir: &mut [f64; 2],
    p_value: f64,
    p_dir_value: f64,
    gradient: &[f64; 2],
    d_gradient: &[f64; 2],
    tol: &[f64; 2],
    f: &mut DirFunction,
) -> bool {
    if !p_value.is_finite() || !p_dir_value.is_finite() {
        return false;
    }
    // (0) Evaluation d'une tolerance parametrique 1D.
    let mut good = false;
    let eps = 1e-20;
    let mut tol1d = 1.1f64;
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

    // (1) On realise une premiere interpolation quadratique.
    let mut ax;
    let mut bx;
    let mut cx;
    let mut df1;
    let mut df2;
    let mut tsol;
    let mut fsol;
    let mut tsolbis;
    let mut delta;

    df1 = gradient[0] * dir[0] + gradient[1] * dir[1];
    df2 = d_gradient[0] * dir[0] + d_gradient[1] * dir[1];

    if df1 < -eps && df2 > eps {
        // cuvette
        tsol = -df1 / (df2 - df1);
    } else {
        cx = p_value;
        bx = df1;
        ax = p_dir_value - (bx + cx);

        if ax.abs() <= eps {
            // cas lineaire
            if bx.abs() >= eps {
                tsol = -cx / bx;
            } else {
                tsol = 0.0;
            }
        } else {
            // cas quadratique
            delta = bx * bx - 4.0 * ax * cx;
            if delta > 1e-9 {
                // il y a des racines, on prend la plus proche de 0.
                delta = delta.sqrt();
                tsol = -(bx + delta);
                tsolbis = delta - bx;
                if tsolbis.abs() < tsol.abs() {
                    tsol = tsolbis;
                }
                tsol /= 2.0 * ax;
            } else {
                // pas ou peu de racine : on "extremise".
                tsol = -(0.5 * bx) / ax;
            }
        }
    }

    if tsol.abs() >= 1.0 {
        return false; // resultat sans interet
    }

    f.initialize(*p, *dir);
    fsol = match f.value(tsol) {
        Some(v) => v,
        None => return false,
    };

    if fsol < p_value {
        good = true;
        result = fsol;
    }

    // (2) Si l'on a pas assez progresse on realise une recherche en bonne et
    //     due forme, a partir des inits precedents.
    if (fsol > 0.2 * p_value) && (tol1d < 0.5) {
        if tsol < 0.0 {
            ax = tsol;
            bx = 0.0;
            cx = 1.0;
        } else {
            ax = 0.0;
            bx = tsol;
            cx = 1.0;
        }

        let mut sol = BrentMinimum::new(tol1d, 100, tol1d);
        sol.perform(f, ax, bx, cx);
        if sol.is_done() {
            if sol.minimum() <= result {
                tsol = sol.location();
                good = true;
                result = sol.minimum();

                // Objective function changes too fast -> perform additional
                // computations.
                if (gradient[0] * gradient[0] + gradient[1] * gradient[1])
                    > 1.0 / rcad_kernel::precision::SQUARE_CONFUSION
                    && tsol > ax
                    && tsol < cx
                {
                    // First and second part invocation.
                    let mut sol2 = BrentMinimum::new(tol1d, 100, tol1d);
                    sol2.perform(f, ax, (ax + tsol) / 2.0, tsol);
                    if sol2.is_done() {
                        if sol2.minimum() <= result {
                            tsol = sol2.location();
                            good = true;
                            result = sol2.minimum();
                        }
                    }

                    let mut sol3 = BrentMinimum::new(tol1d, 100, tol1d);
                    sol3.perform(f, tsol, (cx + tsol) / 2.0, cx);
                    if sol3.is_done() {
                        if sol3.minimum() <= result {
                            tsol = sol3.location();
                            good = true;
                            result = sol3.minimum();
                        }
                    }
                }
            }
        }
    }

    if good {
        // mise a jour du Delta.
        dir[0] *= tsol;
        dir[1] *= tsol;
    }
    good
}

/// OCCT SearchDirection (math_FunctionSetRoot.cxx L439-531) — for Ninc > Neq
/// the SVD minimum-norm solve reduces to Direction = -FF * Gradient/|G|^2.
fn search_direction(
    df: &[[f64; 2]; 1],
    gh: &[f64; 2],
    ff: &[f64; 1],
    change_direction: bool,
    inv_length_max: &[f64; 2],
    direction: &mut [f64; 2],
    dy: &mut f64,
) {
    let ninc = 2;
    let _neq = 1;
    let eps = 1e-32;
    let mut change = change_direction;
    if !change {
        // Ninc > Neq: minimum-norm Gauss-Newton step for DF.Direction = -FF.
        let grad2 = df[0][0] * df[0][0] + df[0][1] * df[0][1];
        if grad2 > 1e-64 {
            let inv = -ff[0] / grad2;
            direction[0] = inv * df[0][0];
            direction[1] = inv * df[0][1];
        } else {
            direction[0] = 0.0;
            direction[1] = 0.0;
            change = true;
        }
    }

    // Il vaut mieux interdire des directions trop longues.
    let mut ratio = (direction[0] * inv_length_max[0]).abs();
    ratio = ratio.max((direction[1] * inv_length_max[1]).abs());
    if ratio > 1.0 {
        direction[0] /= ratio;
        direction[1] /= ratio;
    }

    *dy = direction[0] * gh[0] + direction[1] * gh[1];
    if *dy >= -eps {
        change = true;
    }
    if change {
        // On va faire un gradient!
        direction[0] = -gh[0];
        direction[1] = -gh[1];
        *dy = -(gh[0] * gh[0] + gh[1] * gh[1]);
    }
}

/// OCCT Bounds (math_FunctionSetRoot.cxx L623-705).
#[allow(clippy::too_many_arguments)]
fn bounds(
    inf_bound: &[f64; 2],
    sup_bound: &[f64; 2],
    tol: &[f64; 2],
    sol: &mut [f64; 2],
    sol_save: &[f64; 2],
    constraints: &mut [i32; 2],
    delta: &mut [f64; 2],
    the_is_new_sol: &mut bool,
) -> bool {
    let mut out = false;
    let mut monratio = 1.0f64;

    *the_is_new_sol = true;
    for i in 0..2 {
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
            *the_is_new_sol = false;
            sol.copy_from_slice(sol_save);
            delta[0] = 0.0;
            delta[1] = 0.0;
        } else {
            delta[0] *= monratio;
            delta[1] *= monratio;
            sol[0] = sol_save[0] + delta[0];
            sol[1] = sol_save[1] + delta[1];
            for i in 0..2 {
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

/// OCCT math_FunctionSetRoot — Newton solver for F(u,v) = 0 with bounds.
pub struct FunctionSetRoot {
    done: bool,
    sol: [f64; 2],
    state: i32,
    itermax: i32,
    tol: [f64; 2],
}

impl FunctionSetRoot {
    /// OCCT math_FunctionSetRoot(F, Tolerance, NbIterations = 100).
    pub fn new(_f: &mut SurfFunction, tol: [f64; 2]) -> Self {
        FunctionSetRoot {
            done: false,
            sol: [0.0; 2],
            state: 0,
            itermax: 100,
            tol,
        }
    }

    /// OCCT SetTolerance.
    pub fn set_tolerance(&mut self, tol: [f64; 2]) {
        self.tol = tol;
    }

    /// OCCT Perform(F, StartingPoint, InfBound, SupBound) (L796-1100).
    pub fn perform(
        &mut self,
        f: &mut SurfFunction,
        starting_point: [f64; 2],
        inf_bound: [f64; 2],
        sup_bound: [f64; 2],
    ) {
        let ninc = 2;

        let mut inv_length_max = [0.0; 2];
        for i in 0..ninc {
            let a_sup_bound = sup_bound[i].min(rcad_kernel::precision::INFINITE_VALUE);
            let an_inf_bound = inf_bound[i].max(-rcad_kernel::precision::INFINITE_VALUE);
            inv_length_max[i] = 1.0 / ((a_sup_bound - an_inf_bound) / 4.0).max(1e-9);
        }

        let mut f_dir = DirFunction::new(f);
        let mut descente_iter: i32;

        self.done = false;
        self.sol = starting_point;

        // Verification de la validite des inconnues par rapport aux bornes.
        for i in 0..ninc {
            if self.sol[i] <= inf_bound[i] {
                self.sol[i] = inf_bound[i];
            } else if self.sol[i] > sup_bound[i] {
                self.sol[i] = sup_bound[i];
            }
        }

        // Calcul de la premiere valeur de F et de son gradient.
        let mut ff = [0.0; 1];
        let mut df = [[0.0; 2]; 1];
        let mut gh = [0.0; 2];
        let mut dh = [0.0; 2];
        let mut dh_save = [0.0; 2];
        let mut delta = [0.0; 2];
        let mut f2 = 0.0f64;
        let mut gnr1 = 0.0f64;
        if !f_dir.value_vec(self.sol, &mut ff, &mut df, &mut gh, &mut f2, &mut gnr1) {
            self.done = false;
            self.state = f.get_state_number();
            return;
        }
        let mut ambda2 = gnr1;
        let mut save0 = f2.max(EPS_SQRT);
        let _a_tol_func = f.tolerance();

        if f2 <= EPS || gnr1 <= EPS2 {
            self.done = false;
            let _ = &mut save0;
            self.done = true;
            self.state = f.get_state_number();
            return;
        }

        let mut ambda: f64;
        let mut previous_minimum: f64;
        let mut old_f: f64;
        let mut sol_save = [0.0; 2];
        let mut previous_solution = [0.0; 2];
        let mut constraints = [0i32; 2];
        let mut change_direction = false;
        let mut sort = false;
        let mut is_new_sol = false;
        let mut stop = false;
        let mut good = false;
        let mut dy = 0.0f64;

        let mut kount = 0;
        while kount < self.itermax {
            kount += 1;
            previous_minimum = f2;
            old_f = gnr1;
            previous_solution = self.sol;
            sol_save = self.sol;

            change_direction = false;
            search_direction(&df, &gh, &ff, change_direction, &inv_length_max, &mut dh, &mut dy);
            if dy.abs() <= EPS {
                self.done = false;
                let _ = f.value(&self.sol);
                self.done = true;
                self.state = f.get_state_number();
                return;
            }
            if change_direction {
                ambda = ambda2 / dy.abs().sqrt();
                if ambda > 1.0 {
                    ambda = 1.0;
                }
            } else {
                ambda = 1.0;
                let n = (dh[0] * dh[0] + dh[1] * dh[1]).sqrt();
                ambda2 = if n > 0.0 { 0.5 * ambda / n } else { 0.0 };
            }

            for i in 0..ninc {
                self.sol[i] = self.sol[i] + ambda * dh[i];
            }

            sort = bounds(
                &inf_bound,
                &sup_bound,
                &self.tol,
                &mut self.sol,
                &sol_save,
                &mut constraints,
                &mut delta,
                &mut is_new_sol,
            );

            dh_save = gh;
            if is_new_sol {
                if !f_dir.value_vec(self.sol, &mut ff, &mut df, &mut gh, &mut f2, &mut gnr1) {
                    self.done = false;
                    self.state = f.get_state_number();
                    return;
                }
            }

            if f2 <= EPS || gnr1 <= EPS2 {
                self.done = false;
                let _ = f.value(&self.sol);
                self.done = true;
                self.state = f.get_state_number();
                return;
            }

            if sort || (f2 / previous_minimum > PROGRES) {
                dy = gh[0] * dh[0] + gh[1] * dh[1];
                old_f = previous_minimum;
                stop = false;
                good = false;
                descente_iter = 0;
                let mut sortbis;

                // -------------------------------------------
                // Standard processing without boundary handling
                // -------------------------------------------
                if !sort {
                    // if we haven't exited, we try to progress forward.
                    while (f2 / previous_minimum > PROGRES) && !stop {
                        if f2 < old_f && dy < 0.0 {
                            // We try to progress in this direction.
                            descente_iter += 1;
                            sol_save = self.sol;
                            old_f = f2;
                            for i in 0..ninc {
                                self.sol[i] = self.sol[i] + ambda * dh[i];
                            }
                            stop = bounds(
                                &inf_bound,
                                &sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &sol_save,
                                &mut constraints,
                                &mut delta,
                                &mut is_new_sol,
                            );
                            ambda *= 1.7;
                        } else {
                            if f2 >= old_f || f2 >= previous_minimum {
                                good = false;
                                if descente_iter == 0 {
                                    // C'est le premier pas qui flanche, on fait
                                    // une interpolation.
                                    descente_iter += 1;
                                    good = minimize_direction_3(
                                        &previous_solution,
                                        &sol_save,
                                        &self.sol,
                                        old_f,
                                        &mut delta,
                                        &self.tol,
                                        &mut f_dir,
                                    );
                                } else if change_direction
                                    || descente_iter > 1
                                    || old_f > previous_minimum
                                {
                                    // La progression a ete utile, on minimise.
                                    descente_iter += 1;
                                    good = minimize_direction_2(
                                        &sol_save,
                                        &mut delta,
                                        old_f,
                                        f2,
                                        &dh_save,
                                        &gh,
                                        &self.tol,
                                        &mut f_dir,
                                    );
                                }
                                if !good {
                                    self.sol = sol_save;
                                    f2 = old_f;
                                } else {
                                    self.sol[0] = sol_save[0] + delta[0];
                                    self.sol[1] = sol_save[1] + delta[1];
                                    sort = bounds(
                                        &inf_bound,
                                        &sup_bound,
                                        &self.tol,
                                        &mut self.sol,
                                        &sol_save,
                                        &mut constraints,
                                        &mut delta,
                                        &mut is_new_sol,
                                    );
                                }
                                sort = false; // On a rejete le point sur la frontiere
                            }
                            stop = true; // et on sort dans tous les cas...
                        }
                        dh_save = gh;
                        if is_new_sol {
                            if !f_dir.value_vec(
                                self.sol, &mut ff, &mut df, &mut gh, &mut f2, &mut gnr1,
                            ) {
                                self.done = false;
                                self.state = f.get_state_number();
                                return;
                            }
                        }
                        dy = gh[0] * dh[0] + gh[1] * dh[1];
                        if dy.abs() <= EPS {
                            if f2 > old_f {
                                self.sol = sol_save;
                            }
                            self.done = false;
                            let _ = f.value(&self.sol);
                            self.done = true;
                            self.state = f.get_state_number();
                            return;
                        }
                        if descente_iter >= 100 {
                            stop = true;
                        }
                    }
                }

                // --------------------------------------
                //  on passe au traitement des bords
                // --------------------------------------
                if sort {
                    stop = f2 > 1.001 * old_f; // Pour ne pas progresser sur le bord
                    sortbis = sort;
                    descente_iter = 0;
                    while sortbis && ((f2 < old_f) || (descente_iter == 0)) && !stop {
                        descente_iter += 1;
                        // On essaye de progresser sur le bord.
                        sol_save = self.sol;
                        old_f = f2;
                        // Conditional SearchDirection uses constraints; for the
                        // 1-eq case with a fixed unknown the step is zero on it.
                        let mut cond_dir = change_direction;
                        search_direction(
                            &df,
                            &gh,
                            &ff,
                            cond_dir,
                            &inv_length_max,
                            &mut dh,
                            &mut dy,
                        );
                        if dy < -EPS {
                            if cond_dir {
                                ambda = ambda2 / (-dy).sqrt();
                                if ambda > 1.0 {
                                    ambda = 1.0;
                                }
                            } else {
                                ambda = 1.0;
                                let n = (dh[0] * dh[0] + dh[1] * dh[1]).sqrt();
                                ambda2 = if n > 0.0 { 0.5 * ambda / n } else { 0.0 };
                            }
                            for i in 0..ninc {
                                self.sol[i] = self.sol[i] + ambda * dh[i];
                            }
                            sortbis = bounds(
                                &inf_bound,
                                &sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &sol_save,
                                &mut constraints,
                                &mut delta,
                                &mut is_new_sol,
                            );
                            dh_save = gh;
                            if is_new_sol {
                                if !f_dir.value_vec(
                                    self.sol, &mut ff, &mut df, &mut gh, &mut f2, &mut gnr1,
                                ) {
                                    self.done = false;
                                    self.state = f.get_state_number();
                                    return;
                                }
                            }
                            ambda2 = gnr1;
                        } else {
                            stop = true;
                        }

                        while (f2 / previous_minimum > PROGRES) && (f2 < old_f) && !stop {
                            descente_iter += 1;
                            if f2 < old_f && dy < 0.0 {
                                // On essaye de progresser dans cette direction.
                                sol_save = self.sol;
                                old_f = f2;
                                for i in 0..ninc {
                                    self.sol[i] = self.sol[i] + ambda * dh[i];
                                }
                                sortbis = bounds(
                                    &inf_bound,
                                    &sup_bound,
                                    &self.tol,
                                    &mut self.sol,
                                    &sol_save,
                                    &mut constraints,
                                    &mut delta,
                                    &mut is_new_sol,
                                );
                            }
                            dh_save = gh;
                            if is_new_sol {
                                if !f_dir.value_vec(
                                    self.sol, &mut ff, &mut df, &mut gh, &mut f2, &mut gnr1,
                                ) {
                                    self.done = false;
                                    self.state = f.get_state_number();
                                    return;
                                }
                            }
                            ambda2 = gnr1;
                            dy = gh[0] * dh[0] + gh[1] * dh[1];
                            stop = (dy >= 0.0) || (descente_iter >= 10) || sortbis;
                        }
                        stop = (dy >= 0.0) || (descente_iter >= 10);
                    }
                    if ((f2 / previous_minimum > PROGRES) && (f2 >= old_f)) || (f2 >= previous_minimum)
                    {
                        // On minimise par Brent.
                        descente_iter += 1;
                        good = minimize_direction_2(
                            &sol_save,
                            &mut delta,
                            old_f,
                            f2,
                            &dh_save,
                            &gh,
                            &self.tol,
                            &mut f_dir,
                        );
                        if !good {
                            self.sol = sol_save;
                            sort = false;
                        } else {
                            self.sol[0] = sol_save[0] + delta[0];
                            self.sol[1] = sol_save[1] + delta[1];
                            sort = bounds(
                                &inf_bound,
                                &sup_bound,
                                &self.tol,
                                &mut self.sol,
                                &sol_save,
                                &mut constraints,
                                &mut delta,
                                &mut is_new_sol,
                            );
                            if is_new_sol {
                                if !f_dir.value_vec(
                                    self.sol, &mut ff, &mut df, &mut gh, &mut f2, &mut gnr1,
                                ) {
                                    self.done = false;
                                    self.state = f.get_state_number();
                                    return;
                                }
                            }
                        }
                    }
                }
            }
        }
        self.done = false;
        self.state = f.get_state_number();
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Root(Sol) — the found root.
    pub fn root(&self) -> [f64; 2] {
        self.sol
    }
}
