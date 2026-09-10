//! Math support classes required by the Extrema translations:
//! - `math_MultipleVarFunctionWithHessian` (TKMath/math) — the trait,
//! - `math_NewtonMinimum` (math_NewtonMinimum.cxx L31-278 + .lxx),
//! - `math_Powell` (math_Powell.cxx L36-247 + .lxx) with the file-local
//!   `DirFunctionBis` (cxx L36-78) and the static `MinimizeDirection`
//!   (cxx L80-105),
//! - `NCollection_CellFilter` (NCollection_CellFilter.hxx L112-506) reduced to
//!   the member subset the Extrema objective uses (Reset/Add/Inspect), keeping
//!   the hash-map-of-cells structure and the newest-first list order.
//!
//! rcad placement note: the E3-Q file domain is `base/extrema*.rs`, so these
//! TKMath/NCollection classes live here instead of `math/`; the bodies are the
//! OCCT ones.  Interface request: relocate to `math/opt` and
//! `math/collection` once the domain opens.
//!
//! rcad encodings of the OCCT vector arithmetic
//! (`math_Vector::Add/Subtract/Multiply/Divide/Norm/Norm2`):
//! the kernel `math::math_matrix::Vector` carries storage only, so the
//! operators are free functions here ([`vec_add`], [`vec_sub`],
//! [`vec_scale`], [`vec_norm`], [`vec_norm2`]).

use crate::math::math_bfgs::{BracketMinimum, MultipleVarFunction, MultipleVarFunctionWithGradient, ScalarFunction};
use crate::math::math_matrix::{Matrix, Vector};
use crate::math::opt::BrentMinimum;
use crate::math::root::FunctionValue;

// =============================================================================
// math_Vector arithmetic (math_Vector.cxx / .lxx)
// =============================================================================

/// OCCT `math_Vector::Add(theRight)` — component-wise sum.
pub(crate) fn vec_add(a: &Vector, b: &Vector) -> Vector {
    let mut r = Vector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set(i, a.get(i) + b.get(i));
    }
    r
}

/// OCCT `math_Vector::Subtract(theRight)`.
pub(crate) fn vec_sub(a: &Vector, b: &Vector) -> Vector {
    let mut r = Vector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set(i, a.get(i) - b.get(i));
    }
    r
}

/// OCCT `math_Vector::Multiply(theRight)` — scalar multiply.
pub(crate) fn vec_scale(a: &Vector, s: f64) -> Vector {
    let mut r = Vector::new(a.lower(), a.upper());
    for i in a.lower()..=a.upper() {
        r.set(i, a.get(i) * s);
    }
    r
}

/// OCCT `math_Vector::Norm2()`.
pub(crate) fn vec_norm2(a: &Vector) -> f64 {
    let mut s = 0.0;
    for i in a.lower()..=a.upper() {
        s += a.get(i) * a.get(i);
    }
    s
}

/// OCCT `math_Vector::Norm()`.
pub(crate) fn vec_norm(a: &Vector) -> f64 {
    vec_norm2(a).sqrt()
}

/// OCCT `math_Vector::Dump` / the tuple-copy the OCCT assignments perform.
pub(crate) fn vec_copy_assign(dst: &mut Vector, src: &Vector) {
    for i in dst.lower()..=dst.upper() {
        dst.set(i, src.get(i));
    }
}

// =============================================================================
// math_MultipleVarFunctionWithHessian (hxx L8-64)
// =============================================================================

/// OCCT math_MultipleVarFunctionWithHessian (math_MultipleVarFunctionWithHessian.hxx).
pub trait MultipleVarFunctionWithHessian: MultipleVarFunctionWithGradient {
    /// OCCT Values(X, F, G, H) — the value, the gradient and the Hessian.
    fn values_hessian(&mut self, x: &Vector, f: &mut f64, g: &mut Vector, h: &mut Matrix) -> bool;
}

// =============================================================================
// math_NewtonMinimum (math_NewtonMinimum.hxx L31-129, .cxx L31-278, .lxx)
// =============================================================================

/// OCCT math_NewtonMinimum (hxx L31-129).
pub struct MathNewtonMinimum {
    /// hxx L109: math_Status TheStatus.
    the_status: MathStatus,
    /// hxx L110: math_Vector TheLocation.
    the_location: Vector,
    /// hxx L111: math_Vector TheGradient.
    the_gradient: Vector,
    /// hxx L112: math_Vector TheStep.
    the_step: Vector,
    /// hxx L113: math_Matrix TheHessian.
    the_hessian: Matrix,
    /// hxx L114: double PreviousMinimum.
    previous_minimum: f64,
    /// hxx L115: double TheMinimum.
    the_minimum: f64,
    /// hxx L116: double MinEigenValue.
    min_eigen_value: f64,
    /// hxx L117: double XTol.
    x_tol: f64,
    /// hxx L118: double CTol.
    c_tol: f64,
    /// hxx L119: int nbiter.
    nbiter: i32,
    /// hxx L120: bool NoConvexTreatement.
    no_convex_treatement: bool,
    /// hxx L121: bool Convex.
    convex: bool,
    /// hxx L122: bool myIsBoundsDefined.
    my_is_bounds_defined: bool,
    /// hxx L123: math_Vector myLeft.
    my_left: Vector,
    /// hxx L124: math_Vector myRight.
    my_right: Vector,
    /// hxx L127: bool Done.
    done: bool,
    /// hxx L128: int Itermax.
    itermax: i32,
}

/// OCCT math_Status (math_Status.hxx) — the consumed enumerators.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MathStatus {
    /// math_OK.
    Ok,
    /// math_TooManyIterations.
    TooManyIterations,
    /// math_FunctionError.
    FunctionError,
    /// math_NotBracketed.
    NotBracketed,
    /// math_DirectionSearchError.
    DirectionSearchError,
}

impl MathNewtonMinimum {
    /// OCCT math_NewtonMinimum(theFunction, theTolerance, theNbIterations,
    /// theConvexity, theWithSingularity) (cxx L33-57).
    pub fn new(
        the_function: &dyn MultipleVarFunctionWithHessian,
        the_tolerance: f64,
        the_nb_iterations: i32,
        the_convexity: f64,
        the_with_singularity: bool,
    ) -> Self {
        let n = the_function.nb_variables();
        MathNewtonMinimum {
            the_status: MathStatus::NotBracketed,
            the_location: Vector::new(1, n),
            the_gradient: Vector::new(1, n),
            the_step: Vector::new_init(1, n, 10.0 * the_tolerance),
            the_hessian: Matrix::new(1, n, 1, n),
            previous_minimum: 0.0,
            the_minimum: 0.0,
            min_eigen_value: 0.0,
            x_tol: the_tolerance,
            c_tol: the_convexity,
            nbiter: 0,
            no_convex_treatement: the_with_singularity,
            convex: true,
            my_is_bounds_defined: false,
            my_left: Vector::new_init(1, n, 0.0),
            my_right: Vector::new_init(1, n, 0.0),
            done: false,
            itermax: the_nb_iterations,
        }
    }

    /// OCCT SetBoundary(theLeftBorder, theRightBorder) (cxx L67-73).
    pub fn set_boundary(&mut self, the_left_border: &Vector, the_right_border: &Vector) {
        vec_copy_assign(&mut self.my_left, the_left_border);
        vec_copy_assign(&mut self.my_right, the_right_border);
        self.my_is_bounds_defined = true;
    }

    /// OCCT Perform(F, StartingPoint) (cxx L77-263).
    pub fn perform(&mut self, f: &mut dyn MultipleVarFunctionWithHessian, starting_point: &Vector) {
        let n = f.nb_variables();
        let mut point1 = Vector::new(1, n);
        vec_copy_assign(&mut point1, starting_point);
        let mut point2 = Vector::new(1, n);

        // cxx L86-92.
        let mut ok = true;
        let mut nb_conv = 0;
        let mut v_precedent;
        let mut v_ithere;

        self.done = true;
        self.the_status = MathStatus::Ok;
        self.nbiter = 0;

        // The C++ swaps two pointers; the Rust body mirrors that with an
        // explicit flag selecting which of the two vectors is `precedent`.
        let mut precedent_is_point1 = true;

        // cxx L94: while (Ok && (NbConv < 2))
        while ok && nb_conv < 2 {
            self.nbiter += 1;

            // cxx L100: Ok = F.Values(*precedent, VPrecedent, TheGradient, TheHessian).
            let v_prec_val;
            {
                let precedent = if precedent_is_point1 { &point1 } else { &point2 };
                let mut val = 0.0;
                let mut grad = Vector::new(1, n);
                let mut hess = Matrix::new(1, n, 1, n);
                let vals = f.values_hessian(precedent, &mut val, &mut grad, &mut hess);
                if !vals {
                    self.done = false;
                    self.the_status = MathStatus::FunctionError;
                    return;
                }
                v_prec_val = val;
                self.the_gradient = grad;
                self.the_hessian = hess;
            }
            v_precedent = v_prec_val;
            if self.nbiter == 1 {
                self.previous_minimum = v_precedent;
                self.the_minimum = v_precedent;
            }

            // cxx L115-142: treatment of the non-convexity.
            let jacobi = crate::math::math_jacobi::MathJacobi::new(&self.the_hessian.data);
            if !jacobi.is_done() {
                self.done = false;
                self.the_status = MathStatus::FunctionError;
                return;
            }
            // OCCT: MinEigenValue = Values()(Values().Min()) — the smallest
            // eigenvalue by index.
            let eig = jacobi.values();
            let mut min_idx = 1usize;
            for idx in 2..=eig.len() {
                if eig.get(idx) < eig.get(min_idx) {
                    min_idx = idx;
                }
            }
            self.min_eigen_value = eig.get(min_idx);
            if self.min_eigen_value < self.c_tol {
                self.convex = false;
                if self.no_convex_treatement && self.min_eigen_value.abs() > self.c_tol {
                    // cxx L130-134.
                    let delta = self.c_tol + 0.1 * self.min_eigen_value.abs()
                        - self.min_eigen_value;
                    for ii in 1..=self.the_gradient.length() {
                        let v = self.the_hessian.get(ii, ii) + delta;
                        self.the_hessian.set(ii, ii, v);
                    }
                } else {
                    self.done = false;
                    self.the_status = MathStatus::FunctionError;
                    return;
                }
            }

            // cxx L146-153: Newton scheme (math_Gauss LU(TheHessian, CTol/100)).
            let gauss = crate::math::math_gauss::MathGauss::new(&self.the_hessian.data);
            if !gauss.is_done() {
                self.done = false;
                self.the_status = MathStatus::DirectionSearchError;
                return;
            }
            let mut step = Vector::new(1, n);
            {
                let grad_data = self.the_gradient.data.clone();
                let mut sol = grad_data;
                gauss.solve(&mut sol);
                for i in 1..=n {
                    step.set(i, sol.get(i as usize));
                }
            }
            self.the_step = step;

            // cxx L155-205: project the point on the bounds.
            if self.my_is_bounds_defined {
                let mut mult = f64::MAX; // RealLast()
                let precedent = if precedent_is_point1 { &point1 } else { &point2 };
                for an_idx in 1..=self.my_left.upper() {
                    let an_abs_step = self.the_step.get(an_idx).abs();
                    if an_abs_step < crate::base::extrema_ext_elc::GP_RESOLUTION {
                        continue;
                    }
                    let next_val = precedent.get(an_idx) - self.the_step.get(an_idx);
                    if next_val < self.my_left.get(an_idx) {
                        let a_value =
                            (precedent.get(an_idx) - self.my_left.get(an_idx)).abs() / an_abs_step;
                        mult = mult.min(a_value);
                    }
                    if next_val > self.my_right.get(an_idx) {
                        let a_value =
                            (precedent.get(an_idx) - self.my_right.get(an_idx)).abs() / an_abs_step;
                        mult = mult.min(a_value);
                    }
                }
                if mult != f64::MAX {
                    if mult > crate::core::precision::PCONFUSION {
                        // Project the point into the parameter space.
                        let new_step = vec_scale(&self.the_step, mult);
                        self.the_step = new_step;
                    } else {
                        // cxx L192-202: nullify the step on the boundary axes.
                        let precedent = if precedent_is_point1 { &point1 } else { &point2 };
                        for an_idx in 1..=self.my_left.upper() {
                            if ((precedent.get(an_idx) - self.my_right.get(an_idx)).abs()
                                < crate::core::precision::PCONFUSION
                                && self.the_step.get(an_idx) < 0.0)
                                || ((precedent.get(an_idx) - self.my_left.get(an_idx)).abs()
                                    < crate::core::precision::PCONFUSION
                                    && self.the_step.get(an_idx) > 0.0)
                            {
                                self.the_step.set(an_idx, 0.0);
                            }
                        }
                    }
                }
            }

            // cxx L207-219: convergence guard.
            let mut the_minimum = 0.0;
            loop {
                let next = {
                    let precedent = if precedent_is_point1 { &point1 } else { &point2 };
                    vec_sub(precedent, &self.the_step)
                };
                if precedent_is_point1 {
                    vec_copy_assign(&mut point2, &next);
                } else {
                    vec_copy_assign(&mut point1, &next);
                }
                let has_problem = {
                    let suivant = if precedent_is_point1 { &point2 } else { &point1 };
                    let mut val = 0.0;
                    let ok_val = f.value(suivant, &mut val);
                    the_minimum = val;
                    !ok_val
                };
                if has_problem {
                    let new_step = vec_scale(&self.the_step, 1.0 / 2.0);
                    self.the_step = new_step;
                } else {
                    break;
                }
            }
            self.the_minimum = the_minimum;

            // cxx L221-228.
            if self.is_converged() {
                nb_conv += 1;
            } else {
                nb_conv = 0;
            }

            // cxx L230-241.
            v_ithere = self.the_minimum;
            self.the_minimum = self.previous_minimum;
            let mut n_reduction = 0;
            while v_ithere > v_precedent && n_reduction < 10 {
                let new_step = vec_scale(&self.the_step, 0.4);
                self.the_step = new_step;
                let suivant = {
                    let precedent = if precedent_is_point1 { &point1 } else { &point2 };
                    vec_sub(precedent, &self.the_step)
                };
                if precedent_is_point1 {
                    vec_copy_assign(&mut point2, &suivant);
                } else {
                    vec_copy_assign(&mut point1, &suivant);
                }
                let suivant = if precedent_is_point1 { &point2 } else { &point1 };
                let mut val = 0.0;
                f.value(suivant, &mut val);
                v_ithere = val;
                n_reduction += 1;
            }

            // cxx L243-260.
            if v_ithere <= v_precedent {
                precedent_is_point1 = !precedent_is_point1;
                self.previous_minimum = v_precedent;
                self.the_minimum = v_ithere;
                ok = self.nbiter < self.itermax;
                if !ok && nb_conv < 2 {
                    self.the_status = MathStatus::TooManyIterations;
                }
            } else {
                ok = false;
                self.the_status = MathStatus::DirectionSearchError;
            }
        }
        // cxx L262: TheLocation = *precedent.
        let precedent = if precedent_is_point1 { &point1 } else { &point2 };
        let loc = precedent.clone();
        self.the_location = loc;
    }

    /// OCCT IsConverged() (math_NewtonMinimum.lxx L22-26).
    pub fn is_converged(&self) -> bool {
        vec_norm(&self.the_step) <= self.x_tol
            || (self.the_minimum - self.previous_minimum).abs()
                <= self.x_tol * self.previous_minimum.abs()
    }

    /// OCCT IsDone() (lxx L29-32).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT IsConvex() (hxx L63).
    pub fn is_convex(&self) -> bool {
        self.convex
    }

    /// OCCT Location() (lxx L35-39).
    pub fn location(&self) -> &Vector {
        if !self.done {
            panic!("math_NewtonMinimum: NotDone");
        }
        &self.the_location
    }

    /// OCCT Minimum() (lxx L55-59).
    pub fn minimum(&self) -> f64 {
        if !self.done {
            panic!("math_NewtonMinimum: NotDone");
        }
        self.the_minimum
    }

    /// OCCT NbIterations() (lxx L70-74).
    pub fn nb_iterations(&self) -> i32 {
        if !self.done {
            panic!("math_NewtonMinimum: NotDone");
        }
        self.nbiter
    }

    /// OCCT GetStatus() (lxx L76-79).
    pub fn get_status(&self) -> MathStatus {
        self.the_status
    }
}

// =============================================================================
// math_Powell (math_Powell.hxx L30-105, .cxx L36-247, .lxx)
// =============================================================================

/// OCCT DirFunctionBis (math_Powell.cxx L36-78) — the function restricted to a
/// direction, `P = P0 + x * Dir`.
struct DirFunctionBis<'a> {
    p0: Vector,
    dir: Vector,
    p: Vector,
    f: &'a mut dyn MultipleVarFunction,
}

impl<'a> DirFunctionBis<'a> {
    /// OCCT DirFunctionBis(V1, V2, V3, f) (cxx L52-61).
    fn new(v1: Vector, v2: Vector, v3: Vector, f: &'a mut dyn MultipleVarFunction) -> Self {
        DirFunctionBis {
            p0: v1,
            dir: v2,
            p: v3,
            f,
        }
    }

    /// OCCT Initialize(p0, dir) (cxx L63-68).
    fn initialize(&mut self, p0: &Vector, dir: &Vector) {
        vec_copy_assign(&mut self.p0, p0);
        vec_copy_assign(&mut self.dir, dir);
    }
}

impl ScalarFunction for DirFunctionBis<'_> {
    /// OCCT DirFunctionBis::Value(x, fval) (cxx L70-78).
    fn value(&mut self, x: f64, fval: &mut f64) -> bool {
        let p = vec_add(&vec_scale(&self.dir, x), &self.p0);
        self.p = p;
        *fval = 0.0;
        self.f.value(&self.p, fval)
    }
}

/// The `math_Function` view of [`DirFunctionBis`] required by the kernel
/// `BrentMinimum` (OCCT math_BrentMinimum takes a math_Function, i.e. the
/// `Value(x, f)` -> bool interface the ScalarFunction trait already carries;
/// the kernel BrentMinimum consumes `FunctionValue`, so this adapter re-exposes
/// the same body — no behavioural difference, a pure interface bridge).
struct DirFunctionValueAdapter<'a, 'b>(&'a mut DirFunctionBis<'b>);

impl FunctionValue for DirFunctionValueAdapter<'_, '_> {
    fn value(&mut self, x: f64) -> Option<f64> {
        let mut fval = 0.0;
        if self.0.value(x, &mut fval) {
            Some(fval)
        } else {
            None
        }
    }
}

/// OCCT static MinimizeDirection(P, Dir, Result, F) (cxx L80-105).
fn minimize_direction(
    p: &mut Vector,
    dir: &mut Vector,
    result: &mut f64,
    f_dir: &mut DirFunctionBis<'_>,
) -> bool {
    f_dir.initialize(p, dir);

    // cxx L89-90.
    let mut bracket = BracketMinimum::new(0.0, 1.0);
    if !bracket.is_done() {
        return false;
    }
    let (ax, xx, bx) = bracket.values();
    let mut solver = BrentMinimum::new(1.0e-10, 100, 1.0e-12);
    {
        let mut adapter = DirFunctionValueAdapter(f_dir);
        solver.perform(&mut adapter, ax, xx, bx);
    }
    if solver.is_done() {
        let scale = solver.location();
        *result = solver.minimum();
        *dir = vec_scale(dir, scale);
        *p = vec_add(p, dir);
        return true;
    }
    false
}

/// OCCT math_Powell (hxx L30-105).
pub struct MathPowell {
    /// hxx L78: math_Vector TheLocation.
    the_location: Vector,
    /// hxx L79: double TheMinimum.
    the_minimum: f64,
    /// hxx L81: double PreviousMinimum.
    previous_minimum: f64,
    /// hxx L83: double XTol.
    x_tol: f64,
    /// hxx L84: double EPSZ.
    epsz: f64,
    /// hxx L87: bool Done.
    done: bool,
    /// hxx L88: int Iter.
    iter: i32,
    /// hxx L89: math_Status TheStatus.
    the_status: MathStatus,
    /// hxx L90: math_Matrix TheDirections.
    the_directions: Matrix,
    /// hxx L91: int State.
    state: i32,
    /// hxx L92: int Itermax.
    itermax: i32,
}

impl MathPowell {
    /// OCCT math_Powell(theFunction, theTolerance, theNbIterations, theZEPS)
    /// (cxx L109-126).
    pub fn new(
        the_function: &dyn MultipleVarFunction,
        the_tolerance: f64,
        the_nb_iterations: i32,
        the_zeps: f64,
    ) -> Self {
        let n = the_function.nb_variables();
        MathPowell {
            the_location: Vector::new(1, n),
            the_minimum: f64::MAX,      // RealLast()
            previous_minimum: f64::MAX, // RealLast()
            x_tol: the_tolerance,
            epsz: the_zeps,
            done: false,
            iter: 0,
            the_status: MathStatus::NotBracketed,
            the_directions: Matrix::new(1, n, 1, n),
            state: 0,
            itermax: the_nb_iterations,
        }
    }

    /// OCCT Perform(F, StartingPoint, StartingDirections) (cxx L134-229).
    pub fn perform(
        &mut self,
        f: &mut dyn MultipleVarFunction,
        starting_point: &Vector,
        starting_directions: &Matrix,
    ) {
        self.done = false;
        let n = self.the_location.length();
        let mut pt = Vector::new(1, n);
        let mut ptt = Vector::new(1, n);
        let mut xit = Vector::new(1, n);

        vec_copy_assign(&mut self.the_location, starting_point);
        self.the_directions = starting_directions.clone();
        vec_copy_assign(&mut pt, &self.the_location);

        let mut temp1 = Vector::new(1, n);
        let mut temp2 = Vector::new(1, n);
        let mut temp3 = Vector::new(1, n);
        // cxx L140: double t, fptt, del.
        let mut fptt = 0.0f64;

        for iter in 1..=self.itermax {
            self.iter = iter;
            let mut previous_minimum = 0.0;
            f.value(&self.the_location, &mut previous_minimum);
            self.previous_minimum = previous_minimum;
            let mut ibig = 0i32;
            let mut del = 0.0f64;
            for i in 1..=n {
                for j in 1..=n {
                    xit.set(j, self.the_directions.get(j, i));
                }
                f.value(&self.the_location, &mut fptt);

                let mut the_minimum = 0.0;
                let mut f_dir = DirFunctionBis::new(temp1, temp2, temp3, f);
                let is_good =
                    minimize_direction(&mut self.the_location, &mut xit, &mut the_minimum, &mut f_dir);
                // The DirFunctionBis borrows `f`; rebind the temporaries the OCCT
                // DirFunctionBis owns (V1/V2/V3 are pass-by-reference scratch in
                // OCCT; here the scratch is copied back for the next iteration).
                temp1 = f_dir.p0.clone();
                temp2 = f_dir.dir.clone();
                temp3 = f_dir.p.clone();
                self.the_minimum = the_minimum;

                if !is_good {
                    self.done = false;
                    self.the_status = MathStatus::DirectionSearchError;
                    return;
                }
                if (fptt - the_minimum).abs() > del {
                    del = (fptt - the_minimum).abs();
                    ibig = i as i32;
                }
            }

            if self.is_solution_reached() {
                // OCCT: State = F.GetStateNumber() — math_MultipleVarFunction's
                // state hook; the kernel trait carries no such query, so the
                // recorded state stays at its constructed value.
                self.state = 0;
                self.done = true;
                self.the_status = MathStatus::Ok;
                return;
            }

            if iter == self.itermax {
                self.done = false;
                self.the_status = MathStatus::TooManyIterations;
                return;
            }

            ptt = vec_sub(&vec_scale(&self.the_location, 2.0), &pt);
            xit = vec_sub(&self.the_location, &pt);
            vec_copy_assign(&mut pt, &self.the_location);

            f.value(&ptt, &mut fptt);

            if fptt < self.previous_minimum {
                let t = 2.0
                    * (self.previous_minimum - 2.0 * self.the_minimum + fptt)
                    * (self.previous_minimum - self.the_minimum - del).powi(2)
                    - del * (self.previous_minimum - fptt).powi(2);
                if t < 0.0 {
                    let mut the_minimum = 0.0;
                    let mut f_dir = DirFunctionBis::new(temp1, temp2, temp3, f);
                    let is_good = minimize_direction(
                        &mut self.the_location,
                        &mut xit,
                        &mut the_minimum,
                        &mut f_dir,
                    );
                    temp1 = f_dir.p0.clone();
                    temp2 = f_dir.dir.clone();
                    temp3 = f_dir.p.clone();
                    self.the_minimum = the_minimum;
                    if !is_good {
                        self.done = false;
                        self.the_status = MathStatus::FunctionError;
                        return;
                    }
                    for j in 1..=n {
                        self.the_directions.set(j, ibig, xit.get(j));
                    }
                }
            }
        }
    }

    /// OCCT IsSolutionReached(F) (lxx L26-30).
    pub fn is_solution_reached(&self) -> bool {
        2.0 * (self.previous_minimum - self.the_minimum).abs()
            <= self.x_tol * (self.previous_minimum.abs() + self.the_minimum.abs() + self.epsz)
    }

    /// OCCT IsDone() (lxx L32-35).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT Location() (lxx L45-49).
    pub fn location(&self) -> &Vector {
        if !self.done {
            panic!("math_Powell: NotDone");
        }
        &self.the_location
    }

    /// OCCT Minimum() (lxx L57-61).
    pub fn minimum(&self) -> f64 {
        if !self.done {
            panic!("math_Powell: NotDone");
        }
        self.the_minimum
    }

    /// OCCT NbIterations() (lxx L63-67).
    pub fn nb_iterations(&self) -> i32 {
        if !self.done {
            panic!("math_Powell: NotDone");
        }
        self.iter
    }
}

// =============================================================================
// NCollection_CellFilter (NCollection_CellFilter.hxx L112-506) — used subset
// =============================================================================

/// OCCT NCollection_CellFilter_Action (hxx L26-30).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CellFilterAction {
    /// CellFilter_Keep.
    Keep,
    /// CellFilter_Purge.
    Purge,
}

/// OCCT NCollection_LocalArray<int, 10> cell index — a fixed-size index vector.
pub(crate) type CellIndex = Vec<i64>;

/// OCCT NCollection_CellFilter (hxx L112-506) with the member subset the
/// Extrema objective uses: `Reset`, `Add(target, point)`,
/// `Inspect(min, max, inspector)`.
///
/// The cell map is a hash map keyed by the cell index; every cell holds the
/// list of targets, most recently added first (OCCT prepends to a singly
/// linked list, hxx L355-359), and `Inspect` walks the list in that order
/// (hxx L452-469).
pub(crate) struct CellFilter<T, I> {
    /// hxx L500: int myDim.
    my_dim: usize,
    /// hxx L503: NCollection_Array1<double> myCellSize.
    my_cell_size: Vec<f64>,
    /// hxx L502: CellMap myCells.
    my_cells: std::collections::HashMap<CellIndex, Vec<T>>,
    _inspector: std::marker::PhantomData<I>,
}

impl<T: Clone, I> CellFilter<T, I> {
    /// OCCT NCollection_CellFilter(theCellSize) (hxx L140-146) for the
    /// compile-time dimension form.
    pub(crate) fn new(dim: usize, the_cell_size: f64) -> Self {
        let mut this = CellFilter {
            my_dim: dim,
            my_cell_size: vec![0.0; dim],
            my_cells: std::collections::HashMap::new(),
            _inspector: std::marker::PhantomData,
        };
        this.reset_scalar(the_cell_size);
        this
    }

    /// OCCT Reset(theCellSize) (hxx L149-154).
    pub(crate) fn reset_scalar(&mut self, the_cell_size: f64) {
        for i in 0..self.my_dim {
            self.my_cell_size[i] = the_cell_size;
        }
        // OCCT resetAllocator(theAlloc) (hxx L340-347): clear the cells.
        self.my_cells.clear();
    }

    /// OCCT Reset(theCellSize Array1) (hxx L157-162).
    pub(crate) fn reset_array(&mut self, the_cell_size: &[f64]) {
        self.my_cell_size = the_cell_size.to_vec();
        self.my_cells.clear();
    }

    /// OCCT Cell(thePnt, theCellSize) (hxx L252-267) — the cell index of a
    /// point.
    fn cell_index<C: Fn(usize, &T) -> f64>(&self, pnt: &T, coord: &C) -> CellIndex {
        const INT_MAX: f64 = i32::MAX as f64;
        const INT_MIN: f64 = i32::MIN as f64;
        let mut index = vec![0i64; self.my_dim];
        for i in 0..self.my_dim {
            let a_val = coord(i, pnt) / self.my_cell_size[i];
            index[i] = if a_val > INT_MAX - 1.0 {
                a_val % INT_MAX
            } else if a_val < INT_MIN + 1.0 {
                a_val % INT_MIN
            } else {
                a_val
            } as i64;
        }
        index
    }

    /// OCCT Add(theTarget, thePnt) (hxx L165-169) — into the single cell that
    /// contains the point.
    pub(crate) fn add<C: Fn(usize, &T) -> f64>(&mut self, the_target: &T, the_pnt: &T, coord: &C) {
        let index = self.cell_index(the_pnt, coord);
        let cell = self.my_cells.entry(index).or_default();
        // hxx L355-359: the new node is prepended to the cell's list.
        cell.insert(0, the_target.clone());
    }

    /// OCCT Add(theTarget, thePntMin, thePntMax) (hxx L174-182) — into all the
    /// cells covered by the range.
    pub(crate) fn add_range<C: Fn(usize, &T) -> f64>(
        &mut self,
        the_target: &T,
        the_pnt_min: &T,
        the_pnt_max: &T,
        coord: &C,
    ) {
        let cell_min = self.cell_index(the_pnt_min, coord);
        let cell_max = self.cell_index(the_pnt_max, coord);
        let mut index = cell_min.clone();
        self.iterate_add(self.my_dim - 1, &mut index, &cell_min, &cell_max, the_target);
    }

    /// OCCT iterateAdd (hxx L364-384).
    fn iterate_add(
        &mut self,
        idim: usize,
        the_index: &mut CellIndex,
        the_min_index: &CellIndex,
        the_max_index: &CellIndex,
        the_target: &T,
    ) {
        let a_start = the_min_index[idim];
        let an_end = the_max_index[idim];
        let mut i = a_start;
        while i <= an_end {
            the_index[idim] = i;
            if idim > 0 {
                self.iterate_add(idim - 1, the_index, the_min_index, the_max_index, the_target);
            } else {
                let cell = self.my_cells.entry(the_index.clone()).or_default();
                cell.insert(0, the_target.clone());
            }
            i += 1;
        }
    }

    /// OCCT Inspect(thePntMin, thePntMax, theInspector) (hxx L217-225) — the
    /// inspector's `Inspect` is invoked on every target of every cell in the
    /// range, in the OCCT (newest-first) order.
    pub(crate) fn inspect_range<C, F>(
        &mut self,
        the_pnt_min: &T,
        the_pnt_max: &T,
        coord: &C,
        inspect: &mut F,
    ) where
        C: Fn(usize, &T) -> f64,
        F: FnMut(&T) -> CellFilterAction,
    {
        let cell_min = self.cell_index(the_pnt_min, coord);
        let cell_max = self.cell_index(the_pnt_max, coord);
        let mut cell = cell_min.clone();
        self.iterate_inspect(self.my_dim - 1, &mut cell, &cell_min, &cell_max, inspect);
    }

    /// OCCT iterateInspect (hxx L477-497).
    fn iterate_inspect<F>(
        &mut self,
        idim: usize,
        the_cell: &mut CellIndex,
        the_cell_min: &CellIndex,
        the_cell_max: &CellIndex,
        inspect: &mut F,
    ) where
        F: FnMut(&T) -> CellFilterAction,
    {
        let a_start = the_cell_min[idim];
        let an_end = the_cell_max[idim];
        let mut i = a_start;
        while i <= an_end {
            the_cell[idim] = i;
            if idim > 0 {
                self.iterate_inspect(idim - 1, the_cell, the_cell_min, the_cell_max, inspect);
            } else {
                // OCCT inspect(Cell, Inspector) (hxx L443-474): walk the list and
                // drop the entries the inspector purges.
                if let Some(list) = self.my_cells.get_mut(the_cell) {
                    let mut kept: Vec<T> = Vec::with_capacity(list.len());
                    for obj in list.iter() {
                        if inspect(obj) == CellFilterAction::Keep {
                            kept.push(obj.clone());
                        }
                    }
                    *list = kept;
                }
                if let Some(list) = self.my_cells.get(the_cell) {
                    if list.is_empty() {
                        self.my_cells.remove(the_cell);
                    }
                }
            }
            i += 1;
        }
    }
}
