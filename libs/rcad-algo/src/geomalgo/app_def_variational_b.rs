//! OCCT AppDef_Variational - private engine methods, part 2
//! (AppDef_Variational.cxx L1160-2126): TheMotor, Optimization, Project,
//! ACR, the statics NearIndex / GettingKnots and SplitCurve.
//!
//! Split off from `app_def_variational.rs` via `#[path]`; the parent's
//! private fields and helpers are visible here (child module).

use std::cell::RefCell;
use std::rc::Rc;

use rcad_kernel::math::fem_tool::{Assembly, Curve};
use rcad_kernel::math::math_matrix::{Matrix, Vector};

use super::{stable_sort_real_array, AppDefVariational, CurveHandle, RealArray1, CriterionHandle};

impl AppDefVariational {
    /// OCCT AppDef_Variational::TheMotor (cxx L1160-1590).  The first
    /// floating argument is unnamed in OCCT (the commented-out WQuadratic).
    /// `allow(unused_assignments)`: the OCCT keeps dead stores (L1253
    /// `lconst = true`, L1194-ish `u_new` reset) that Rust would flag.
    #[allow(unused_assignments)]
    pub(super) fn the_motor(
        &mut self,
        j: &CriterionHandle,
        _w_quadratic: f64, // OCCT: unnamed `const double` (WQuadratic, unused)
        w_quality: f64,
        the_curve: &mut CurveHandle,
        ecarts: &mut RealArray1,
    ) {
        // ...

        const BIG_VALUE: f64 = 1.0e37;
        const SMALL_VALUE: f64 = 1.0e-6;
        const SMALLEST_VALUE: f64 = 1.0e-9;

        // OCCT handles CurrentTi / NewTi / OldTi - mirrored by values with
        // explicit clones at the OCCT handle assignments (no OCCT statement
        // mutates through an alias of these arrays).
        let mut current_ti: RealArray1;
        let mut new_ti: RealArray1;
        let mut old_ti: RealArray1;
        let dependence;
        let mut lestim;
        let mut to_optim;
        let mut iscut = false;
        let mut isnear = false;
        let mut again = true;
        let mut nb_est;
        let mut icdana;
        let mut num_pnt = 0; // OCCT: uninitialized `int NumPnt` (set by Project)
        let mut iter;
        let max_nb_est = 5;
        let mut vocri = [BIG_VALUE, BIG_VALUE, BIG_VALUE];
        let mut erold = BIG_VALUE;
        let mut valcri = [0.0f64; 3]; // OCCT: uninitialized VALCRI[3]
        let mut errmax = BIG_VALUE;
        let mut errmoy = 0.0; // OCCT: uninitialized double
        let mut errqua = 0.0; // OCCT: uninitialized double
        let mut cblong;
        let mut lnold;
        let nbr_pnt = self.my_last_point - self.my_first_point + 1;
        let nbr_constraint = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;
        let mut c_current: CurveHandle;
        let mut c_old: CurveHandle;
        // OCCT: null handle until (1.1); initialized for Rust flow analysis
        // (the OCCT handle is always assigned before any read).
        let mut c_new: CurveHandle = the_curve.clone();
        let mut eps_length = SMALL_VALUE;
        let mut eps_deg;
        // OCCT: uninitialized doubles e1 / e2 / e3.
        let mut e1 = 0.0;
        let mut e2 = 0.0;
        let mut e3 = 0.0;
        let j1min;
        let j2min;
        let j3min;
        let mut iprog;

        // (0) Init

        j.borrow().get_estimation(&mut e1, &mut e2, &mut e3);
        j1min = 1.0e-8;
        j2min = (e1 + 1.0e-8) * 1.0e-6;
        j3min = j2min;

        if e1 < j1min {
            e1 = j1min; // Like in
        }
        if e2 < j2min {
            e2 = j2min; // MOTLIS
        }
        if e3 < j3min {
            e3 = j3min;
        }

        j.borrow_mut().set_estimation(e1, e2, e3);

        c_current = the_curve.clone();
        // CurrentTi = new HArray1(1, myParameters->Length());
        // CurrentTi->ChangeArray1() = myParameters->Array1(); (bounds + data)
        current_ti = self.my_parameters.clone();
        // OldTi = new HArray1(1, CurrentTi->Length());
        // OldTi->ChangeArray1() = CurrentTi->Array1();
        old_ti = current_ti.clone();
        c_old = c_current.clone();
        cblong = *j.borrow_mut().est_length();
        lnold = cblong;
        dependence = j.borrow().dependence_table();

        j.borrow_mut().set_curve(&c_current);
        // OCCT: FEmTool_Assembly* TheAssembly = new FEmTool_Assembly(...);
        // (the OCCT `delete TheAssembly` at cxx L1544/L1589 is the Rust drop.)
        let mut the_assembly =
            Assembly::new(&dependence, j.borrow().assembly_table());

        //============        Optimization      ============================
        //  int inagain = 0;
        'great_loop: while again {

            // (1) Loop  Optimization / Estimation
            lestim = true;
            let mut lconst = true;
            nb_est = 0;

            j.borrow_mut().set_curve(&c_current);

            while lestim {

                //     (1.1) Curve's Optimization.
                eps_length = SMALL_VALUE * cblong / nbr_pnt as f64;
                c_new = {
                    let cc = c_current.borrow();
                    Rc::new(RefCell::new(Curve::new(
                        cc.dimension(),
                        cc.nb_elements(),
                        cc.base(),
                        eps_length,
                    )))
                };
                // CNew->Knots() = CCurrent->Knots(); (same element count)
                {
                    let src = c_current.borrow_mut().knots().to_vec();
                    c_new.borrow_mut().knots().copy_from_slice(&src);
                }

                j.borrow_mut().set_parameters(&current_ti);
                eps_deg = (w_quality * 0.1).min(cblong * 0.001);

                self.optimization(j, &mut the_assembly, lconst, eps_deg, &c_new, &current_ti);

                lconst = false;

                //        (1.2) calculation of quality criteria and improvement
                //              of estimation.
                // OCCT passes &VALCRI[0..2]; the Rust array needs
                // split_at_mut for the three disjoint &mut cells.
                let (vc0, vc12) = valcri.split_at_mut(1);
                let (vc1, vc2) = vc12.split_at_mut(1);
                icdana = j.borrow_mut().quality_values(
                    j1min,
                    j2min,
                    j3min,
                    &mut vc0[0],
                    &mut vc1[0],
                    &mut vc2[0],
                );

                if icdana > 0 {
                    lconst = true;
                }

                j.borrow_mut().error_values(&mut errmax, &mut errqua, &mut errmoy);

                isnear = ((errqua / nbr_pnt as f64).sqrt() < 2.0 * w_quality)
                    && (self.my_nb_iterations > 1);

                //       (1.3) Optimization of ti by orthogonal projection
                //             and calculation of the error at points.

                if isnear {
                    // NewTi = new HArray1(1, CurrentTi->Length());
                    new_ti = RealArray1::new(1, current_ti.length());
                    self.project(
                        &c_new,
                        &current_ti,
                        &mut new_ti,
                        ecarts,
                        &mut num_pnt,
                        &mut errmax,
                        &mut errqua,
                        &mut errmoy,
                        2,
                    );
                } else {
                    new_ti = current_ti.clone(); // OCCT: NewTi = CurrentTi (handle alias)
                }

                //        (1.4) Progression's test
                iprog = 0;
                if (erold > w_quality) && (errmax < 0.95 * erold) {
                    iprog += 1;
                }
                if (erold > w_quality) && (errmax < 0.8 * erold) {
                    iprog += 1;
                }
                if (erold > w_quality) && (errmax < w_quality) {
                    iprog += 1;
                }
                if (erold > w_quality) && (errmax < 0.99 * erold) && (errmax < 1.1 * w_quality) {
                    iprog += 1;
                }
                if valcri[0] < 0.975 * vocri[0] {
                    iprog += 1;
                }
                if valcri[0] < 0.9 * vocri[0] {
                    iprog += 1;
                }
                if valcri[1] < 0.95 * vocri[1] {
                    iprog += 1;
                }
                if valcri[1] < 0.8 * vocri[1] {
                    iprog += 1;
                }
                if valcri[2] < 0.95 * vocri[2] {
                    iprog += 1;
                }
                if valcri[2] < 0.8 * vocri[2] {
                    iprog += 1;
                }
                if (vocri[1] > SMALLEST_VALUE) && (vocri[2] > SMALLEST_VALUE) {
                    if (valcri[1] / vocri[1] + 2.0 * valcri[2] / vocri[2]) < 2.8 {
                        iprog += 1;
                    }
                }

                if iprog < 2 && nb_est == 0 {
                    //             (1.5) Invalidation of new knots.
                    valcri[0] = vocri[0];
                    valcri[1] = vocri[1];
                    valcri[2] = vocri[2];
                    errmax = erold;
                    cblong = lnold;
                    c_current = c_old.clone();
                    current_ti = old_ti.clone();

                    break 'great_loop; // OCCT goto L8000 (exit)
                }

                vocri[0] = valcri[0];
                vocri[1] = valcri[1];
                vocri[2] = valcri[2];
                lnold = cblong;
                erold = errmax;

                c_current = c_new.clone();
                current_ti = new_ti.clone();

                //       (1.6) Test if the Estimations seems  OK, else repeat
                nb_est += 1;
                lestim = (nb_est < max_nb_est) && (icdana == 2) && (iprog > 0);

                if lestim && isnear {
                    //           (1.7) Optimization of ti by ACR.

                    stable_sort_real_array(&mut current_ti);

                    let decima = 4;

                    c_current.borrow_mut().length(0.0, 1.0, &mut cblong);
                    *j.borrow_mut().est_length() = cblong;

                    self.acr(&c_current, &mut current_ti, decima);
                    lconst = true;
                }
            }

            //     (2) loop of parametric / geometric optimization

            iter = 1;
            to_optim = (iter < self.my_nb_iterations) && (isnear);

            while to_optim {
                iter += 1;
                //     (2.1) Save current results
                vocri[0] = valcri[0];
                vocri[1] = valcri[1];
                vocri[2] = valcri[2];
                erold = errmax;
                lnold = cblong;
                c_old = c_current.clone();
                // OldTi->ChangeArray1() = CurrentTi->Array1();
                old_ti = current_ti.clone();

                //     (2.2) Optimization of ti by ACR.

                stable_sort_real_array(&mut current_ti);

                let decima = 4;

                c_current.borrow_mut().length(0.0, 1.0, &mut cblong);
                *j.borrow_mut().est_length() = cblong;

                self.acr(&c_current, &mut current_ti, decima);
                lconst = true;

                //      (2.3) Optimization of curves
                eps_length = SMALL_VALUE * cblong / nbr_pnt as f64;

                c_new = {
                    let cc = c_current.borrow();
                    Rc::new(RefCell::new(Curve::new(
                        cc.dimension(),
                        cc.nb_elements(),
                        cc.base(),
                        eps_length,
                    )))
                };
                {
                    let src = c_current.borrow_mut().knots().to_vec();
                    c_new.borrow_mut().knots().copy_from_slice(&src);
                }

                j.borrow_mut().set_parameters(&current_ti);

                eps_deg = (w_quality * 0.1).min(cblong * 0.001);
                self.optimization(j, &mut the_assembly, lconst, eps_deg, &c_new, &current_ti);

                c_current = c_new.clone();

                //      (2.4) calculation of quality criteria and improvement
                //             of estimation.
                // OCCT passes &VALCRI[0..2]; the Rust array needs
                // split_at_mut for the three disjoint &mut cells.
                let (vc0, vc12) = valcri.split_at_mut(1);
                let (vc1, vc2) = vc12.split_at_mut(1);
                icdana = j.borrow_mut().quality_values(
                    j1min,
                    j2min,
                    j3min,
                    &mut vc0[0],
                    &mut vc1[0],
                    &mut vc2[0],
                );
                if icdana > 0 {
                    lconst = true;
                }

                j.borrow().get_estimation(&mut e1, &mut e2, &mut e3);
                //       (2.5) Optimization of ti by orthogonal projection

                new_ti = RealArray1::new(1, current_ti.length());
                self.project(
                    &c_current,
                    &current_ti,
                    &mut new_ti,
                    ecarts,
                    &mut num_pnt,
                    &mut errmax,
                    &mut errqua,
                    &mut errmoy,
                    2,
                );

                //       (2.6)  Test de non regression

                let mut iregre = 0;
                if nbr_constraint < nbr_pnt {
                    if (errmax > w_quality) && (errmax > 1.05 * erold) {
                        iregre += 1;
                    }
                    if (errmax > w_quality) && (errmax > 2.0 * erold) {
                        iregre += 1;
                    }
                    if (erold > w_quality) && (errmax <= 0.5 * erold) {
                        iregre -= 1;
                    }
                }
                let (mut big_e1, mut big_e2, mut big_e3) = (0.0f64, 0.0f64, 0.0f64);
                j.borrow().get_estimation(&mut big_e1, &mut big_e2, &mut big_e3);
                if (valcri[0] > big_e1) && (valcri[0] > 1.1 * vocri[0]) {
                    iregre += 1;
                }
                if (valcri[1] > big_e2) && (valcri[1] > 1.1 * vocri[1]) {
                    iregre += 1;
                }
                if (valcri[2] > big_e3) && (valcri[2] > 1.1 * vocri[2]) {
                    iregre += 1;
                }

                if iregre >= 2 {
                    //      if (iregre >= 1) {
                    // (2.7) restore the previous iteration
                    valcri[0] = vocri[0];
                    valcri[1] = vocri[1];
                    valcri[2] = vocri[2];
                    errmax = erold;
                    cblong = lnold;
                    c_current = c_old.clone();
                    // CurrentTi->ChangeArray1() = OldTi->Array1();
                    current_ti = old_ti.clone();
                    to_optim = false;
                } else {
                    // Iteration is Ok.
                    c_current = c_new.clone();
                    current_ti = new_ti.clone();
                }
                if iter >= self.my_nb_iterations {
                    to_optim = false;
                }
            }

            // (3) Optional splitting

            if (c_current.borrow().nb_elements() < self.my_max_segment) && self.my_with_cutting {

                //    (3.1) Sauvgarde de l'etat precedent
                vocri[0] = valcri[0];
                vocri[1] = valcri[1];
                vocri[2] = valcri[2];
                erold = errmax;
                c_old = c_current.clone();
                // OldTi->ChangeArray1() = CurrentTi->Array1();
                old_ti = current_ti.clone();

                //       (3.2) Arrange the ti: Sort + rescale to (0,1)
                //         ---> Sort to ensure proper ordering afterwards.

                stable_sort_real_array(&mut current_ti);

                if (current_ti.value(1) != 0.0) || (current_ti.value(nbr_pnt) != 1.0) {
                    let mut t;
                    let delat_t = 1.0 / (current_ti.value(nbr_pnt) - current_ti.value(1));
                    for ii in 2..nbr_pnt {
                        t = (current_ti.value(ii) - current_ti.value(1)) * delat_t;
                        current_ti.set_value(ii, t);
                    }
                    current_ti.set_value(1, 0.0);
                    current_ti.set_value(nbr_pnt, 1.0);
                }

                //       (3.3) Insert new Knots

                self.split_curve(&c_current, &current_ti, eps_length, &mut c_new, &mut iscut);
                if !iscut {
                    again = false;
                } else {
                    c_current = c_new.clone();
                    // New Knots => New Assembly.
                    j.borrow_mut().set_curve(&c_new);
                    // delete TheAssembly;
                    the_assembly = Assembly::new(&dependence, j.borrow().assembly_table());
                }
            } else {
                again = false;
            }
        }

        //    ================   Great loop end   ===================

        // L8000:

        // (4) Compute the best Error.
        new_ti = RealArray1::new(1, current_ti.length());
        self.project(
            &c_current,
            &current_ti,
            &mut new_ti,
            ecarts,
            &mut num_pnt,
            &mut errmax,
            &mut errqua,
            &mut errmoy,
            10,
        );

        // (5) field's update

        *the_curve = c_current.clone();
        *j.borrow_mut().est_length() = cblong;
        // myParameters->ChangeArray1() = NewTi->Array1();
        self.my_parameters = new_ti.clone();
        self.my_criterium[0] = errqua;
        self.my_criterium[1] = valcri[0].sqrt();
        self.my_criterium[2] = valcri[1].sqrt();
        self.my_criterium[3] = valcri[2].sqrt();
        self.my_max_error = errmax;
        self.my_max_error_index = num_pnt;
        if nbr_pnt > nbr_constraint {
            self.my_average_error = errmoy / (nbr_pnt - nbr_constraint) as f64;
        } else {
            self.my_average_error = errmoy / nbr_constraint as f64;
        }

        // delete TheAssembly; (drop)
    }

    /// OCCT AppDef_Variational::Optimization (cxx L1593-1689).
    fn optimization(
        &self,
        j: &CriterionHandle,
        a: &mut Assembly,
        to_assemble: bool,
        eps_deg: f64,
        curve: &CurveHandle,
        parameters: &RealArray1,
    ) {
        let (mx_deg, nb_elm, nb_dim) = {
            let c = curve.borrow();
            (c.base().work_degree(), c.nb_elements(), c.dimension())
        };

        // OCCT: math_Matrix H(0, MxDeg, 0, MxDeg);
        //       math_Vector G(0, MxDeg), Sol(1, A.NbGlobVar());
        let mut h = Matrix::new(0, mx_deg, 0, mx_deg);
        let mut g = Vector::new(0, mx_deg);
        let mut sol = Vector::new(1, a.nb_glob_var());

        let nb_constr = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;

        let cblong = *j.borrow_mut().est_length();

        // Updating Assembly
        if to_assemble {
            a.nullify_matrix();
        }
        a.nullify_vector();

        for el in 1..=nb_elm {
            if to_assemble {
                j.borrow_mut().hessian(el, 1, 1, &mut h);
                for dim in 1..=nb_dim {
                    a.add_matrix(el, dim, dim, &h);
                }
            }

            for dim in 1..=nb_dim {
                j.borrow_mut().gradient(el, dim, &mut g);
                a.add_vector(el, dim, &g);
            }
        }

        // Solution of system
        if to_assemble {
            if nb_constr != 0 {
                // Treatment of constraints
                self.assembling_constraints(curve, parameters, cblong, a);
            }
            a.solve();
        }
        a.solution(&mut sol);

        // OCCT L1609: AssTable = A.AssemblyTable(); (an immutable reference;
        // re-bound here because the Rust &mut-A calls above have ended - the
        // table is not modified in between).
        let ass_table = a.assembly_table();

        // Updating J
        j.borrow_mut().set_curve(curve);
        j.borrow_mut().input_vector(&sol, ass_table);

        // Updating Curve and reduction of degree

        let mut newdeg = 0;
        let mut max_error = 0.0;

        if nb_constr == 0 {
            for el in 1..=nb_elm {
                curve.borrow_mut().reduce_degree(el, eps_deg, &mut newdeg, &mut max_error);
            }
        } else {

            // OCCT: NCollection_Array1<double>& TabInt = Curve->Knots(); - a
            // read-only view in this scope (ReduceDegree / SetDegree do not
            // modify the knots); the Rust translation snapshots it.
            let tab_int: Vec<f64> = curve.borrow_mut().knots().to_vec();
            let mut icnt = 1;
            let p0 = parameters.lower() - self.my_first_point;
            for el in 1..=nb_elm {
                while (icnt < nb_constr)
                    && (parameters.value(p0 + self.my_typ_constraints.value(2 * icnt - 1))
                        <= tab_int[(el - 1) as usize])
                {
                    icnt += 1;
                }
                let point = p0 + self.my_typ_constraints.value(2 * icnt - 1);
                if parameters.value(point) <= tab_int[(el - 1) as usize]
                    || parameters.value(point) >= tab_int[el as usize]
                {
                    curve.borrow_mut().reduce_degree(el, eps_deg, &mut newdeg, &mut max_error);
                } else if curve.borrow().degree(el) < mx_deg {
                    curve.borrow_mut().set_degree(el, mx_deg);
                }
            }
        }
    }

    /// OCCT AppDef_Variational::Project (cxx L1691-1827).
    #[allow(clippy::too_many_arguments)]
    fn project(
        &self,
        c: &CurveHandle,
        ti: &RealArray1,
        proj_ti: &mut RealArray1,
        distance: &mut RealArray1,
        num_points: &mut i32,
        max_err: &mut f64,
        qua_err: &mut f64,
        ave_err: &mut f64,
        nb_iterations: i32,
    ) {
        // Initialisation

        let seuil = 1.0e-9;
        let eps = 1.0e-12;

        *max_err = 0.0;
        *qua_err = 0.0;
        *ave_err = 0.0;

        // int i0 = -myDimension, d0 = Distance.Lower() - 1;
        let mut i0 = -self.my_dimension;
        let d0 = distance.lower() - 1;

        // OCCT: NCollection_Array1<double> ValOfC(1, myDimension), ...
        let mut val_of_c = vec![0.0f64; self.my_dimension as usize];
        let mut first_der_of_c = vec![0.0f64; self.my_dimension as usize];
        let mut secnd_der_of_c = vec![0.0f64; self.my_dimension as usize];

        for ipnt in 1..=proj_ti.length() {

            i0 += self.my_dimension;

            let mut t_new = ti.value(ipnt);

            let mut en_cour = true;
            let mut nit_cv = 0;
            let mut iter = 0;
            c.borrow_mut().d0(t_new, &mut val_of_c);

            let mut dist = 0.0;
            for i in 1..=self.my_dimension {
                let aux = val_of_c[(i - 1) as usize] - self.my_tab_points.value(i0 + i);
                dist += aux * aux;
            }
            dist = dist.sqrt();

            // ------- Newton's method for solving (C'(t),C(t) - P) = 0

            while en_cour {

                iter += 1;
                let t0 = t_new;
                let dist0 = dist;

                c.borrow_mut().d2(t_new, &mut secnd_der_of_c);
                c.borrow_mut().d1(t_new, &mut first_der_of_c);

                let mut f1 = 0.0;
                let mut f2 = 0.0;
                for i in 1..=self.my_dimension {
                    let aux = val_of_c[(i - 1) as usize] - self.my_tab_points.value(i0 + i);
                    let df = first_der_of_c[(i - 1) as usize];
                    f1 += aux * df; // (C'(t),C(t) - P)
                    f2 += df * df + aux * secnd_der_of_c[(i - 1) as usize]; // ((C'(t),C(t) - P))'
                }

                if f2.abs() < eps {
                    en_cour = false;
                } else {
                    // Formula of Newton x(k+1) = x(k) - F(x(k))/F'(x(k))
                    t_new -= f1 / f2;
                    if t_new < 0.0 {
                        t_new = 0.0;
                    }
                    if t_new > 1.0 {
                        t_new = 1.0;
                    }

                    // Analysis of result

                    c.borrow_mut().d0(t_new, &mut val_of_c);

                    dist = 0.0;
                    for i in 1..=self.my_dimension {
                        let aux = val_of_c[(i - 1) as usize] - self.my_tab_points.value(i0 + i);
                        dist += aux * aux;
                    }
                    dist = dist.sqrt();

                    let ecart = dist0 - dist;

                    if ecart <= -seuil {
                        // No improvement, stop here
                        en_cour = false;
                        t_new = t0;
                        dist = dist0;
                    } else if ecart <= seuil {
                        // Convergence
                        nit_cv += 1;
                    } else {
                        nit_cv = 0;
                    }

                    if (nit_cv >= 2) || (iter >= nb_iterations) {
                        en_cour = false;
                    }
                }
            }

            proj_ti.set_value(ipnt, t_new);
            distance.set_value(d0 + ipnt, dist);
            if dist > *max_err {
                *max_err = dist;
                *num_points = ipnt;
            }
            *qua_err += dist * dist;
            *ave_err += dist;
        }

        // Setting NumPoints to interval [myFirstPoint, myLastPoint]
        *num_points = *num_points + self.my_first_point - 1;
    }

    /// OCCT AppDef_Variational::ACR (cxx L1829-1959).
    fn acr(&self, curve: &CurveHandle, ti: &mut RealArray1, decima: i32) {

        let eps = 1.0e-8;

        // OCCT: NCollection_Array1<double>& Knots = Curve->Knots(); - an
        // in-place reference; the Rust translation re-borrows the handle at
        // each knot statement (no Curve evaluation intervenes between a knot
        // read and its write, so the aliasing semantics are identical).

        let nbr_pnt = ti.length();
        let ti_first = ti.lower();
        let ti_last = ti.upper();
        // KFirst = Knots.Lower(), KLast = Knots.Upper() (1-based bounds).
        let k_first = 1;
        let k_last = curve.borrow().nb_elements() + 1;

        let mut cb_long = 0.0;
        let delta_t;
        let mut v_test;
        let mut u_new = 0.0;
        let mut u_old;
        let mut du;
        let mut t_para;
        let mut t_old;
        let mut ratio;
        let mut ii;
        let mut i_elm;
        let mut i_old;
        let mut p_old;
        let mut p_cnt;
        let mut i_cnt = 0;
        let nb_cntr = self.my_nb_pass_points + self.my_nb_tang_points + self.my_nb_curv_points;

        //     (1) Compute the curve length

        curve.borrow_mut().length(ti.value(ti_first), ti.value(ti_last), &mut cb_long);

        //     (2)  Place the arc-length parameterization in Ti

        if nbr_pnt >= 2 {

            //     (2.0) Initialisation
            delta_t = (ti.value(ti_last) - ti.value(ti_first)) / decima as f64;
            v_test = ti.value(ti_first) + delta_t;

            if nb_cntr > 0 {
                p_cnt = self.my_typ_constraints.value(1) - self.my_first_point + ti_first;
                i_cnt = 1;
            } else {
                p_cnt = ti_last + 1;
            }

            u_old = 0.0;

            t_old = ti.value(ti_first);
            p_old = ti_first;

            i_elm = k_first;
            i_old = i_elm;

            ti.set_value(ti_first, 0.0);

            let mut ipnt = ti_first + 1;
            while ipnt <= ti_last {

                while (i_cnt <= nb_cntr) && (p_cnt < ipnt) {
                    i_cnt += 1;
                    p_cnt = self.my_typ_constraints.value(2 * i_cnt - 1)
                        - self.my_first_point
                        + ti_first;
                }

                t_para = ti.value(ipnt);

                if t_para >= v_test || p_cnt == ipnt {

                    if ti.value(ti_last) - t_para <= 1.0e-2 * delta_t {
                        ipnt = ti_last;
                        t_para = ti.value(ipnt);
                    }
                    //        (2.2), (2.3) Compute the curve length
                    curve.borrow_mut().length(ti.value(ti_first), t_para, &mut u_new);

                    u_new /= cb_long;

                    // while (Knots(IElm + 1) < TPara && IElm < KLast - 1) IElm++;
                    while curve.borrow_mut().knots()[((i_elm + 1) - 1) as usize] < t_para
                        && i_elm < k_last - 1
                    {
                        i_elm += 1;
                    }

                    //         (2.4) Update the splitting parameters
                    let dt_inv = 1.0 / (t_para - t_old);
                    du = u_new - u_old;

                    for ii in (i_old + 1)..=i_elm {
                        // Ratio = (Knots(ii) - TOld) * DTInv;
                        ratio = (curve.borrow_mut().knots()[(ii - 1) as usize] - t_old) * dt_inv;
                        // Knots(ii) = UOld + Ratio * DU;
                        curve.borrow_mut().knots()[(ii - 1) as usize] = u_old + ratio * du;
                    }

                    //           (2.5) Update the point parameters.

                    // Very strange loop, because it never works (POld+1 > ipnt-1)
                    for ii in (p_old + 1)..=(ipnt - 1) {
                        ratio = (ti.value(ii) - t_old) * dt_inv;
                        ti.set_value(ii, u_old + ratio * du);
                    }

                    ti.set_value(ipnt, u_new);

                    u_old = u_new;
                    i_old = i_elm;
                    t_old = t_para;
                    p_old = ipnt;
                }
                //        --> New parametric threshold for decimation

                if t_para >= v_test {
                    // ii = RealToInt(...); VTest += (ii + 1) * DeltaT; (disabled)
                    v_test += ((t_para - v_test + eps) / delta_t).ceil() * delta_t;
                    if v_test > 1.0 - eps {
                        v_test = 1.0;
                    }
                }
                ipnt += 1;
            }
        }

        //     --- Adjust the extreme values

        ti.set_value(ti_first, 0.0);
        ti.set_value(ti_last, 1.0);
        ii = ti_last - 1;
        while ti.value(ii) > curve.borrow_mut().knots()[(k_last - 1) as usize] {
            ti.set_value(ii, 1.0);
            ii -= 1;
        }
        curve.borrow_mut().knots()[(k_first - 1) as usize] = 0.0;
        curve.borrow_mut().knots()[(k_last - 1) as usize] = 1.0;
    }

    /// OCCT AppDef_Variational::SplitCurve (cxx L2074-2126).
    fn split_curve(
        &self,
        in_curve: &CurveHandle,
        ti: &RealArray1,
        curve_tol: f64,
        out_curve: &mut CurveHandle,
        iscut: &mut bool,
    ) {
        let nb_elm_old = in_curve.borrow().nb_elements();

        if nb_elm_old >= self.my_max_segment {
            *iscut = false;
            return;
        }
        // InCurve->Base().WorkDegree() (the OCCT_DEBUG naming kept out).
        let work_degree = in_curve.borrow().base().work_degree();
        let mut nb_elm = nb_elm_old;
        let mut new_knots = RealArray1::new(nb_elm + 1, self.my_max_segment);
        getting_knots(ti, in_curve, work_degree, &mut nb_elm, &mut new_knots);
        getting_knots(ti, in_curve, work_degree - 1, &mut nb_elm, &mut new_knots);

        if nb_elm > nb_elm_old {

            *iscut = true;

            let out = Rc::new(RefCell::new(Curve::new(
                in_curve.borrow().dimension(),
                nb_elm,
                in_curve.borrow().base(),
                curve_tol,
            )));
            // NCollection_Array1<double>& OutKnots = OutCurve->Knots();
            // NCollection_Array1<double>& InKnots  = InCurve->Knots();
            let in_knots: Vec<f64> = in_curve.borrow_mut().knots().to_vec();
            {
                let mut om = out.borrow_mut();
                let out_knots = om.knots();
                let i0 = 1; // OutKnots.Lower()
                for i in 1..=(nb_elm_old + 1) {
                    out_knots[(i - 1) as usize] = in_knots[(i - 1) as usize];
                }
                for i in (nb_elm_old + 1)..=nb_elm {
                    out_knots[((i + i0) - 1) as usize] = new_knots.value(i);
                }
            }

            // std::sort(OutKnots.begin(), OutKnots.end());
            {
                let mut om = out.borrow_mut();
                let mut v: Vec<f64> = om.knots().to_vec();
                v.sort_by(|x, y| x.partial_cmp(y).unwrap());
                om.knots().copy_from_slice(&v);
            }
            *out_curve = out;
        } else {
            *iscut = false;
        }
    }
}

/// OCCT static NearIndex (AppDef_Variational.cxx L1963-2004).  `tab_par` is
/// the 0-based rcad mirror of the OCCT 1-based Array1 (Lower() == 1); the
/// OCCT 1-based indices are kept and shifted on the slice access.
/// `pub(super)` so the sibling module `variational_c` (Adjusting) can call
/// it - the OCCT static is file-scope in AppDef_Variational.cxx.
pub(super) fn near_index(t: f64, tab_par: &[f64], eps: f64, flag: &mut i32) -> i32 {
    let loi = 1; // TabPar.Lower()
    let upi = tab_par.len() as i32; // TabPar.Upper()

    *flag = 0;

    if t < tab_par[(loi - 1) as usize] {
        *flag = -1;
        return loi;
    }
    if t > tab_par[(upi - 1) as usize] {
        *flag = 1;
        return upi;
    }

    let mut ibeg = loi;
    let mut ifin = upi;
    let mut imidl;

    while ibeg + 1 != ifin {
        imidl = (ibeg + ifin) / 2;
        if (t >= tab_par[(ibeg - 1) as usize]) && (t <= tab_par[(imidl - 1) as usize]) {
            ifin = imidl;
        } else {
            ibeg = imidl;
        }
    }

    if (t - tab_par[(ifin - 1) as usize]).abs() < eps {
        return ifin;
    }

    ibeg
}

/// OCCT static GettingKnots (AppDef_Variational.cxx L2006-2072) - calculating
/// values of new Knots for elements with degree that is equal Deg.
fn getting_knots(
    tab_par: &RealArray1,
    in_curve: &CurveHandle,
    deg: i32,
    nb_elm: &mut i32,
    new_knots: &mut RealArray1,
) {
    let eps = 1.0e-12;

    // OCCT: NCollection_Array1<double>& OldKnots = InCurve->Knots(); -
    // read-only in this routine; snapshot (identical semantics).
    let old_knots: Vec<f64> = in_curve.borrow_mut().knots().to_vec();
    let nb_max_old = in_curve.borrow().nb_elements();
    let nb_max = new_knots.upper();
    let mut flag = 0i32;
    // OCCT: int el = 0, i1 = OldKnots.Lower(), i0 = i1 - 1;
    let mut el = 0;
    let mut i1 = 1;
    let mut i0 = i1 - 1;
    let mut t_par;

    while (*nb_elm < nb_max) && (el < nb_max_old) {

        el += 1;
        i0 += 1;
        i1 += 1; // i0, i1 are indexes of left and right knots of element el

        if in_curve.borrow().degree(el) == deg {

            *nb_elm += 1;

            // TabPar is flat-copied for the 0-based NearIndex mirror.
            let tab_par_flat: Vec<f64> =
                (tab_par.lower()..=tab_par.upper()).map(|i| tab_par.value(i)).collect();
            let mut ipt1 = near_index(old_knots[(i0 - 1) as usize], &tab_par_flat, eps, &mut flag);
            if flag != 0 {
                ipt1 = tab_par.lower(); // TabPar.Lower()
            }
            let mut ipt2 = near_index(old_knots[(i1 - 1) as usize], &tab_par_flat, eps, &mut flag);
            if flag != 0 {
                ipt2 = tab_par.upper(); // TabPar.Upper()
            }

            if ipt2 - ipt1 >= 1 {

                let ipt = (ipt1 + ipt2) / 2;
                if 2 * ipt == ipt1 + ipt2 {
                    t_par = 2.0 * tab_par.value(ipt);
                } else {
                    t_par = tab_par.value(ipt) + tab_par.value(ipt + 1);
                }

                new_knots.set_value(
                    *nb_elm,
                    (old_knots[(i0 - 1) as usize] + old_knots[(i1 - 1) as usize] + t_par) / 4.0,
                );
            } else {
                new_knots.set_value(
                    *nb_elm,
                    (old_knots[(i0 - 1) as usize] + old_knots[(i1 - 1) as usize]) / 2.0,
                );
            }
        }
    }
}
