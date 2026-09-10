//! OCCT math_GlobOptMin (FoundationClasses/TKMath/math/math_GlobOptMin.hxx
//! L50-267 and math_GlobOptMin.cxx L32-669) — Evtushenko's global optimization
//! on a non-uniform mesh.  This is the solver `Extrema_GGenExtCC::Perform`
//! drives for the curve-curve extremum search.
//!
//! Also carries the file-local static `DistanceToBorder` (cxx L32-47) and the
//! nested `NCollection_CellFilter_Inspector` (hxx L131-188).
//!
//! rcad placement note: the body is TKMath; it lives under `base/extrema*.rs`
//! because the E3-Q file domain is restricted (interface request: relocate to
//! `math/opt` when the domain opens).
//!
//! rcad encoding of the C++ `dynamic_cast` chain in `computeLocalExtremum`
//! (cxx L273-336): the objective is held as [`MathMultiVarFunc`], whose
//! variants are the OCCT interfaces the pointer could be cast to, so the three
//! `dynamic_cast` tests become variant tests with the same outcomes.

use crate::base::extrema_math_opt::{
    vec_add, vec_norm, vec_norm2, vec_scale, vec_sub, CellFilter, CellFilterAction,
    MathNewtonMinimum, MathPowell, MultipleVarFunctionWithHessian,
};
use crate::core::precision::{is_infinite_value, PCONFUSION};
use crate::math::math_bfgs::{Bfgs, MultipleVarFunction, MultipleVarFunctionWithGradient};
use crate::math::math_matrix::{Matrix, Vector};

/// OCCT M_SQRT2 (math_GlobOptMin.cxx L26-28).
const M_SQRT2: f64 = 1.41421356237309504880168872420969808;

/// OCCT RealLast().
const REAL_LAST: f64 = f64::MAX;

/// OCCT RealFirst().
const REAL_FIRST: f64 = -f64::MAX;

/// OCCT Precision::Infinite().
const PRECISION_INFINITE: f64 = crate::core::precision::INFINITE_VALUE;

/// The rcad encoding of the `math_MultipleVarFunction*` the GlobOptMin solver
/// receives: the variant names the OCCT interface the pointer refers to, so the
/// `dynamic_cast` tests of `computeLocalExtremum` (cxx L273, L294, L315) become
/// variant tests.
pub enum MathMultiVarFunc<'a> {
    /// A plain `math_MultipleVarFunction`.
    Plain(&'a mut dyn MultipleVarFunction),
    /// A `math_MultipleVarFunctionWithGradient`.
    WithGradient(&'a mut dyn MultipleVarFunctionWithGradient),
    /// A `math_MultipleVarFunctionWithHessian`.
    WithHessian(&'a mut dyn MultipleVarFunctionWithHessian),
}

impl MathMultiVarFunc<'_> {
    /// OCCT `myFunc->NbVariables()`.
    pub(crate) fn nb_variables(&self) -> i32 {
        match self {
            MathMultiVarFunc::Plain(f) => f.nb_variables(),
            MathMultiVarFunc::WithGradient(f) => f.nb_variables(),
            MathMultiVarFunc::WithHessian(f) => f.nb_variables(),
        }
    }

    /// OCCT `myFunc->Value(X, F)`.
    pub(crate) fn value(&mut self, x: &Vector, f: &mut f64) -> bool {
        match self {
            MathMultiVarFunc::Plain(p) => p.value(x, f),
            MathMultiVarFunc::WithGradient(p) => p.value(x, f),
            MathMultiVarFunc::WithHessian(p) => p.value(x, f),
        }
    }
}

/// The `MultipleVarFunction` view of a gradient function — the rcad bridge that
/// lets the kernel `math_BFGS` (generic over `F: Sized`) run on the trait
/// object the OCCT code passes by pointer.  Pure delegation, no behaviour.
struct GradFnRef<'a>(&'a mut dyn MultipleVarFunctionWithGradient);

impl MultipleVarFunction for GradFnRef<'_> {
    fn nb_variables(&self) -> i32 {
        self.0.nb_variables()
    }
    fn value(&mut self, x: &Vector, f: &mut f64) -> bool {
        self.0.value(x, f)
    }
}

impl MultipleVarFunctionWithGradient for GradFnRef<'_> {
    fn gradient(&mut self, x: &Vector, g: &mut Vector) -> bool {
        self.0.gradient(x, g)
    }
    fn values(&mut self, x: &Vector, f: &mut f64, g: &mut Vector) -> bool {
        self.0.values(x, f, g)
    }
}

/// OCCT static DistanceToBorder(theX, theMin, theMax) (cxx L32-47).
fn distance_to_border(the_x: &Vector, the_min: &Vector, the_max: &Vector) -> f64 {
    let mut a_dist = REAL_LAST;
    for an_idx in the_min.lower()..=the_min.upper() {
        let a_dist1 = (the_x.get(an_idx) - the_min.get(an_idx)).abs();
        let a_dist2 = (the_x.get(an_idx) - the_max.get(an_idx)).abs();
        a_dist = a_dist.min(a_dist1.min(a_dist2));
    }
    a_dist
}

/// OCCT `isInside(thePnt)` body (cxx L494-507) as a free function so it can be
/// called while the objective (a field of the solver) is mutably borrowed.
fn is_inside_bounds(the_pnt: &Vector, my_n: i32, glob_a: &Vector, glob_b: &Vector) -> bool {
    for i in 1..=my_n {
        if the_pnt.get(i) < glob_a.get(i) || the_pnt.get(i) > glob_b.get(i) {
            return false;
        }
    }
    true
}

/// OCCT NCollection_CellFilter_Inspector (hxx L131-188) of math_GlobOptMin.
struct GlobOptInspector {
    /// hxx L184: double myTol.
    my_tol: f64,
    /// hxx L185: math_Vector myCurrent.
    my_current: Vector,
    /// hxx L186: bool myIsFind.
    my_is_find: bool,
    /// hxx L187: int Dimension.
    dimension: usize,
}

impl GlobOptInspector {
    /// OCCT NCollection_CellFilter_Inspector(theDim, theTol) (hxx L138-144).
    fn new(the_dim: usize, the_tol: f64) -> Self {
        GlobOptInspector {
            my_tol: the_tol * the_tol,
            my_current: Vector::new(1, the_dim as i32),
            my_is_find: false,
            dimension: the_dim,
        }
    }

    /// OCCT Coord(i, thePnt) (hxx L147).
    fn coord(i: usize, the_pnt: &Vector) -> f64 {
        the_pnt.get((i + 1) as i32)
    }

    /// OCCT Shift(thePnt, theTol, theLowPnt, theUppPnt) (hxx L151-161).
    fn shift(&self, the_pnt: &Vector, the_tol: &[f64]) -> (Vector, Vector) {
        let mut low = Vector::new(1, self.dimension as i32);
        let mut upp = Vector::new(1, self.dimension as i32);
        for an_idx in 1..=self.dimension as i32 {
            low.set(an_idx, the_pnt.get(an_idx) - the_tol[(an_idx - 1) as usize]);
            upp.set(an_idx, the_pnt.get(an_idx) + the_tol[(an_idx - 1) as usize]);
        }
        (low, upp)
    }

    /// OCCT ClearFind() (hxx L163).
    fn clear_find(&mut self) {
        self.my_is_find = false;
    }

    /// OCCT isFind() (hxx L165).
    fn is_find(&self) -> bool {
        self.my_is_find
    }

    /// OCCT SetCurrent(theCurPnt) (hxx L168).
    fn set_current(&mut self, the_cur_pnt: &Vector) {
        self.my_current = the_cur_pnt.clone();
    }

    /// OCCT Inspect(theObject) (hxx L171-181).
    fn inspect(&mut self, the_object: &Vector) -> CellFilterAction {
        let diff = vec_sub(&self.my_current, the_object);
        let a_sq_dist = vec_norm2(&diff);
        if a_sq_dist < self.my_tol {
            self.my_is_find = true;
        }
        CellFilterAction::Keep
    }
}

/// OCCT math_GlobOptMin (hxx L50-267).
pub struct MathGlobOptMin<'a> {
    /// hxx L225: math_MultipleVarFunction* myFunc.
    my_func: MathMultiVarFunc<'a>,
    /// hxx L226: int myN.
    my_n: i32,
    /// hxx L227: math_Vector myA.
    my_a: Vector,
    /// hxx L228: math_Vector myB.
    my_b: Vector,
    /// hxx L229: math_Vector myGlobA.
    my_glob_a: Vector,
    /// hxx L230: math_Vector myGlobB.
    my_glob_b: Vector,
    /// hxx L231: double myTol.
    my_tol: f64,
    /// hxx L232: double mySameTol.
    my_same_tol: f64,
    /// hxx L235: double myC.
    my_c: f64,
    /// hxx L236: double myInitC.
    my_init_c: f64,
    /// hxx L237: bool myIsFindSingleSolution.
    my_is_find_single_solution: bool,
    /// hxx L238: double myFunctionalMinimalValue.
    my_functional_minimal_value: f64,
    /// hxx L239: bool myIsConstLocked.
    my_is_const_locked: bool,
    /// hxx L242: bool myDone.
    my_done: bool,
    /// hxx L243: NCollection_Sequence<double> myY — the solution components,
    /// stored flat as `sol * myN + dim` in OCCT (hxx L526, L551).
    my_y: Vec<f64>,
    /// hxx L244: int mySolCount.
    my_sol_count: usize,
    /// hxx L247: double myZ.
    my_z: f64,
    /// hxx L248: double myE1.
    my_e1: f64,
    /// hxx L249: double myE2.
    my_e2: f64,
    /// hxx L250: double myE3.
    my_e3: f64,
    /// hxx L252: math_Vector myX.
    my_x: Vector,
    /// hxx L253: math_Vector myTmp.
    my_tmp: Vector,
    /// hxx L254: math_Vector myV.
    my_v: Vector,
    /// hxx L255: math_Vector myMaxV.
    my_max_v: Vector,
    /// hxx L256: double myLastStep.
    my_last_step: f64,
    /// hxx L258: NCollection_Array1<double> myCellSize.
    my_cell_size: Vec<f64>,
    /// hxx L259: int myMinCellFilterSol.
    my_min_cell_filter_sol: usize,
    /// hxx L260: bool isFirstCellFilterInvoke.
    is_first_cell_filter_invoke: bool,
    /// hxx L261: NCollection_CellFilter<...> myFilter.
    my_filter: CellFilter<Vector, GlobOptInspector>,
    /// hxx L264: int myCont.
    my_cont: i32,
    /// hxx L266: double myF.
    my_f: f64,
}

impl<'a> MathGlobOptMin<'a> {
    /// OCCT math_GlobOptMin(theFunc, theA, theB, theC, theDiscretizationTol,
    /// theSameTol) (cxx L51-106).
    pub fn new(
        the_func: MathMultiVarFunc<'a>,
        the_a: &Vector,
        the_b: &Vector,
        the_c: f64,
        the_discretization_tol: f64,
        the_same_tol: f64,
    ) -> Self {
        let my_n = the_func.nb_variables();
        let mut this = MathGlobOptMin {
            my_func: the_func,
            my_n,
            my_a: Vector::new(1, my_n),
            my_b: Vector::new(1, my_n),
            my_glob_a: Vector::new(1, my_n),
            my_glob_b: Vector::new(1, my_n),
            my_tol: the_discretization_tol,
            my_same_tol: the_same_tol,
            my_c: the_c,
            my_init_c: the_c,
            my_is_find_single_solution: false,
            my_functional_minimal_value: -PRECISION_INFINITE,
            my_is_const_locked: false,
            my_done: false,
            my_y: Vec::new(),
            my_sol_count: 0,
            my_z: -1.0,
            my_e1: 0.0,
            my_e2: 0.0,
            my_e3: 0.0,
            my_x: Vector::new(1, my_n),
            my_tmp: Vector::new(1, my_n),
            my_v: Vector::new(1, my_n),
            my_max_v: Vector::new(1, my_n),
            my_last_step: 0.0,
            my_cell_size: vec![0.0; my_n as usize],
            my_min_cell_filter_sol: 0,
            is_first_cell_filter_invoke: false,
            my_filter: CellFilter::new(my_n as usize, 0.0),
            my_cont: 2,
            my_f: PRECISION_INFINITE,
        };

        // cxx L72-89.
        this.my_done = false;
        for i in 1..=my_n {
            this.my_glob_a.set(i, the_a.get(i));
            this.my_glob_b.set(i, the_b.get(i));
            this.my_a.set(i, the_a.get(i));
            this.my_b.set(i, the_b.get(i));
        }
        for i in 1..=my_n {
            let v = (this.my_b.get(i) - this.my_a.get(i)) / 3.0;
            this.my_max_v.set(i, v);
        }
        this.my_tol = the_discretization_tol;
        this.my_same_tol = the_same_tol;

        // cxx L99-103.
        let a_max_square_search_sol = 200usize;
        let a_sol_nb = 3.0f64.powi(my_n) as usize;
        this.my_min_cell_filter_sol = (2 * a_sol_nb).max(a_max_square_search_sol);
        this.init_cell_size();
        this.compute_init_sol();
        this.my_done = false;
        this
    }

    /// OCCT SetGlobalParams (cxx L112-148).
    pub fn set_global_params(
        &mut self,
        the_func: MathMultiVarFunc<'a>,
        the_a: &Vector,
        the_b: &Vector,
        the_c: f64,
        the_discretization_tol: f64,
        the_same_tol: f64,
    ) {
        self.my_func = the_func;
        self.my_c = the_c;
        self.my_init_c = the_c;
        self.my_z = -1.0;
        self.my_sol_count = 0;
        for i in 1..=self.my_n {
            self.my_glob_a.set(i, the_a.get(i));
            self.my_glob_b.set(i, the_b.get(i));
            self.my_a.set(i, the_a.get(i));
            self.my_b.set(i, the_b.get(i));
        }
        for i in 1..=self.my_n {
            let v = (self.my_b.get(i) - self.my_a.get(i)) / 3.0;
            self.my_max_v.set(i, v);
        }
        self.my_tol = the_discretization_tol;
        self.my_same_tol = the_same_tol;
        self.init_cell_size();
        self.compute_init_sol();
        self.my_done = false;
    }

    /// OCCT SetLocalParams(theLocalA, theLocalB) (cxx L154-171).
    pub fn set_local_params(&mut self, the_local_a: &Vector, the_local_b: &Vector) {
        self.my_z = -1.0;
        for i in 1..=self.my_n {
            self.my_a.set(i, the_local_a.get(i));
            self.my_b.set(i, the_local_b.get(i));
        }
        for i in 1..=self.my_n {
            let v = (self.my_b.get(i) - self.my_a.get(i)) / 3.0;
            self.my_max_v.set(i, v);
        }
        self.my_done = false;
    }

    /// OCCT SetTol(theDiscretizationTol, theSameTol) (cxx L175-179).
    pub fn set_tol(&mut self, the_discretization_tol: f64, the_same_tol: f64) {
        self.my_tol = the_discretization_tol;
        self.my_same_tol = the_same_tol;
    }

    /// OCCT GetTol(theDiscretizationTol, theSameTol) (cxx L183-187).
    pub fn get_tol(&self) -> (f64, f64) {
        (self.my_tol, self.my_same_tol)
    }

    /// OCCT Perform(isFindSingleSolution) (cxx L192-262).
    pub fn perform(&mut self, is_find_single_solution: bool) {
        self.my_done = false;

        // cxx L197-212: compute the parameters range.
        let mut min_length = REAL_LAST;
        let mut max_length = REAL_FIRST;
        for i in 1..=self.my_n {
            let current_length = self.my_b.get(i) - self.my_a.get(i);
            if current_length < min_length {
                min_length = current_length;
            }
            if current_length > max_length {
                max_length = current_length;
            }
            self.my_v.set(i, 0.0);
        }

        // cxx L214-221.
        if min_length < PCONFUSION {
            return;
        }

        // cxx L223-227.
        if !self.my_is_const_locked {
            self.compute_initial_values();
        }

        // cxx L229-248.
        self.my_e1 = min_length * self.my_tol;
        self.my_e2 = max_length * self.my_tol;
        self.my_is_find_single_solution = is_find_single_solution;
        if is_find_single_solution {
            self.my_e3 = 0.0;
        } else if self.my_c > 1.0 {
            self.my_e3 = -max_length * self.my_tol / 4.0;
        } else {
            self.my_e3 = -max_length * self.my_tol * self.my_c / 4.0;
        }

        // cxx L251-255.
        if self.check_functional_stop_criteria() {
            self.my_done = true;
            return;
        }

        // cxx L257-261.
        self.my_last_step = 0.0;
        self.is_first_cell_filter_invoke = true;
        let n = self.my_n;
        self.compute_global_extremum(n as usize);

        self.my_done = true;
    }

    /// OCCT computeLocalExtremum(thePnt, theVal, theOutPnt) (cxx L266-339).
    fn compute_local_extremum(
        &mut self,
        the_pnt: &Vector,
        the_val: &mut f64,
        the_out_pnt: &mut Vector,
    ) -> bool {
        let my_n = self.my_n;
        let glob_a = self.my_glob_a.clone();
        let glob_b = self.my_glob_b.clone();

        // Newton method (cxx L272-291).
        if self.my_cont >= 2 {
            let mut is_hessian = false;
            if let MathMultiVarFunc::WithHessian(_) = self.my_func {
                is_hessian = true;
            }
            if is_hessian {
                // OCCT math_NewtonMinimum newtonMinimum(*aTmp) — the defaults
                // are Tolerance = Precision::Confusion(), NbIterations = 40,
                // Convexity = 1.0e-6, WithSingularity = true.
                let (mut newton_minimum, f_ref) = match &mut self.my_func {
                    MathMultiVarFunc::WithHessian(f) => (
                        MathNewtonMinimum::new(
                            &**f,
                            crate::core::precision::CONFUSION,
                            40,
                            1.0e-6,
                            true,
                        ),
                        &mut **f,
                    ),
                    _ => unreachable!(),
                };
                newton_minimum.set_boundary(&glob_a, &glob_b);
                newton_minimum.perform(f_ref, the_pnt);

                if newton_minimum.is_done() {
                    *the_out_pnt = newton_minimum.location().clone();
                    *the_val = newton_minimum.minimum();

                    if is_inside_bounds(the_out_pnt, my_n, &glob_a, &glob_b) {
                        return true;
                    }
                }
            }
        }

        // BFGS method used (cxx L293-312).
        if self.my_cont >= 1 {
            let mut is_gradient = false;
            match self.my_func {
                MathMultiVarFunc::WithGradient(_) | MathMultiVarFunc::WithHessian(_) => {
                    is_gradient = true
                }
                _ => {}
            }
            if is_gradient {
                // OCCT math_BFGS bfgs(aTmp->NbVariables()) — defaults
                // Tolerance = 1.0e-8, NbIterations = 200, ZEPS = 1.0e-12.
                let mut bfgs = Bfgs::new(my_n, 1.0e-8, 200, 1.0e-12);
                bfgs.set_boundary(&glob_a, &glob_b);
                let done = match &mut self.my_func {
                    MathMultiVarFunc::WithGradient(p) => {
                        bfgs.perform(&mut GradFnRef(&mut **p), the_pnt);
                        true
                    }
                    MathMultiVarFunc::WithHessian(p) => {
                        bfgs.perform(&mut GradFnRef(&mut **p), the_pnt);
                        true
                    }
                    MathMultiVarFunc::Plain(_) => false,
                };
                if done && bfgs.is_done() {
                    *the_out_pnt = bfgs.location().clone();
                    *the_val = bfgs.minimum();

                    if is_inside_bounds(the_out_pnt, my_n, &glob_a, &glob_b) {
                        return true;
                    }
                }
            }
        }

        // Powell method used (cxx L314-336).
        {
            let mut m = Matrix::new(1, my_n, 1, my_n);
            for i in 1..=my_n {
                m.set(i, i, 1.0);
            }

            // OCCT math_Powell powell(*myFunc, 1e-10) — defaults
            // NbIterations = 200, ZEPS = 1.0e-12.
            let mut powell = {
                let f_ref: &dyn MultipleVarFunction = match &self.my_func {
                    MathMultiVarFunc::Plain(p) => &**p,
                    MathMultiVarFunc::WithGradient(p) => &**p,
                    MathMultiVarFunc::WithHessian(p) => &**p,
                };
                MathPowell::new(f_ref, 1e-10, 200, 1.0e-12)
            };
            let pnt_clone = the_pnt.clone();
            {
                let f_mut: &mut dyn MultipleVarFunction = match &mut self.my_func {
                    MathMultiVarFunc::Plain(p) => &mut **p,
                    MathMultiVarFunc::WithGradient(p) => &mut **p,
                    MathMultiVarFunc::WithHessian(p) => &mut **p,
                };
                powell.perform(f_mut, &pnt_clone, &m);
            }
            if powell.is_done() {
                *the_out_pnt = powell.location().clone();
                *the_val = powell.minimum();

                if is_inside_bounds(the_out_pnt, my_n, &glob_a, &glob_b) {
                    return true;
                }
            }
        }

        false
    }

    /// OCCT computeInitialValues() (cxx L343-388).
    fn compute_initial_values(&mut self) {
        let a_min_lc = 0.01;
        let a_max_lc = 1000.0;
        let a_min_eps = 0.1;
        let a_max_eps = 100.0;

        // Lipchitz const approximation (cxx L355-376).
        let mut a_lip_const = 0.0;
        let a_pnt_nb = 13usize;
        let a_ref = self.my_a.clone();
        let mut a_prev_val_diag = 0.0;
        self.my_func.value(&a_ref, &mut a_prev_val_diag);
        let mut a_prev_val_proj = a_prev_val_diag;
        let diff = vec_sub(&self.my_b, &self.my_a);
        let a_step = vec_norm(&diff) / a_pnt_nb as f64;
        let a_param_step = vec_scale(&diff, 1.0 / a_pnt_nb as f64);
        for i in 1..=a_pnt_nb {
            let mut a_curr_pnt = vec_add(&self.my_a, &vec_scale(&a_param_step, i as f64));

            // Walk over the diagonal.
            let mut a_curr_val = 0.0;
            self.my_func.value(&a_curr_pnt, &mut a_curr_val);
            a_lip_const = (a_curr_val - a_prev_val_diag).abs().max(a_lip_const);
            a_prev_val_diag = a_curr_val;

            // Walk over the diagonal in the projected space (aPnt(1) = myA(1)).
            a_curr_pnt.set(1, self.my_a.get(1));
            self.my_func.value(&a_curr_pnt, &mut a_curr_val);
            a_lip_const = (a_curr_val - a_prev_val_proj).abs().max(a_lip_const);
            a_prev_val_proj = a_curr_val;
        }

        // cxx L378-387.
        self.my_c = self.my_init_c;
        let a_lip_const = a_lip_const * (self.my_n as f64).sqrt() / a_step;
        if a_lip_const < self.my_c * a_min_eps {
            self.my_c = (a_lip_const * a_min_eps).max(a_min_lc);
        } else if a_lip_const > self.my_c * a_max_eps {
            self.my_c = (self.my_c * a_max_eps).min(a_max_lc);
        }
    }

    /// OCCT computeGlobalExtremum(j) (cxx L392-490).
    fn compute_global_extremum(&mut self, j: usize) {
        let mut d = REAL_LAST;
        let mut a_prev_val = 0.0;
        let mut val = REAL_LAST;
        let mut a_step_best_value = REAL_LAST;
        let mut a_step_best_point = Vector::new(1, self.my_n);
        let mut is_inside = false;
        let mut is_reached = false;

        // cxx L403: for (myX(j) = myA(j) + myE1; !isReached; myX(j) += myV(j))
        let jj = j as i32;
        let start = self.my_a.get(jj) + self.my_e1;
        self.my_x.set(jj, start);
        while !is_reached {
            if self.my_x.get(jj) > self.my_b.get(jj) {
                self.my_x.set(jj, self.my_b.get(jj));
                is_reached = true;
            }

            // cxx L411-414.
            if self.check_functional_stop_criteria() {
                return;
            }

            if j == 1 {
                // cxx L416-462.
                is_inside = false;
                a_prev_val = d;
                {
                    let x_ref = self.my_x.clone();
                    self.my_func.value(&x_ref, &mut d);
                }
                let r1 = (d + self.my_z * self.my_c * self.my_last_step - self.my_f) * self.my_z;
                let r2 =
                    ((d + a_prev_val - self.my_c * self.my_last_step) * 0.5 - self.my_f) * self.my_z;
                let r = r1.min(r2);
                if r > self.my_e3 {
                    let a_save_param = self.my_x.get(1);

                    // Piyavsky midpoint estimation (cxx L430-434).
                    let mut a_param = (2.0 * self.my_x.get(1) - self.my_v.get(1)) * 0.5
                        + (a_prev_val - d) * 0.5 / self.my_c;
                    if is_infinite_value(a_prev_val) {
                        a_param = self.my_x.get(1) - self.my_v.get(1) * 0.5;
                    }

                    self.my_x.set(1, a_param);
                    let mut a_val = 0.0;
                    {
                        let x_ref = self.my_x.clone();
                        self.my_func.value(&x_ref, &mut a_val);
                    }
                    self.my_x.set(1, a_save_param);

                    // cxx L442-447.
                    if (a_val < d && a_val < a_prev_val)
                        || distance_to_border(&self.my_x, &self.my_a, &self.my_b) < self.my_e1
                    {
                        let x_clone = self.my_x.clone();
                        let mut tmp = Vector::new(1, self.my_n);
                        is_inside = self.compute_local_extremum(&x_clone, &mut val, &mut tmp);
                        if is_inside {
                            self.my_tmp = tmp;
                        }
                    }
                }
                a_step_best_value = if is_inside && (val < d) { val } else { d };
                a_step_best_point = if is_inside && (val < d) {
                    self.my_tmp.clone()
                } else {
                    self.my_x.clone()
                };

                // cxx L453.
                self.check_add_candidate(&a_step_best_point, a_step_best_value);

                if self.check_functional_stop_criteria() {
                    return;
                }

                // cxx L460-461.
                let v1 = (self.my_e2 + (self.my_f - d).abs() / self.my_c).min(self.my_max_v.get(1));
                self.my_v.set(1, v1);
                self.my_last_step = v1;
            } else {
                // cxx L464-473.
                self.my_v.set(jj, REAL_LAST / 2.0);
                self.compute_global_extremum(j - 1);

                // Nullify the steps on the lower dimensions.
                for i in 1..j {
                    self.my_v.set(i as i32, 0.0);
                }
            }

            // cxx L474-488.
            if j < self.my_n as usize {
                let a_upper_dim_step = self.my_v.get(jj).max(self.my_e2);
                if self.my_v.get(jj + 1) > a_upper_dim_step {
                    if a_upper_dim_step > self.my_max_v.get(jj + 1) {
                        self.my_v.set(jj + 1, self.my_max_v.get(jj + 1));
                    } else {
                        self.my_v.set(jj + 1, a_upper_dim_step);
                    }
                }
            }

            // cxx L403: myX(j) += myV(j).
            let new_val = self.my_x.get(jj) + self.my_v.get(jj);
            self.my_x.set(jj, new_val);
        }
    }

    /// OCCT isStored(thePnt) (cxx L511-573).
    fn is_stored(&mut self, the_pnt: &Vector) -> bool {
        let mut a_tol = Vector::new(1, self.my_n);
        let diff = vec_sub(&self.my_b, &self.my_a);
        for j in 1..=self.my_n {
            a_tol.set(j, diff.get(j) * self.my_same_tol);
        }

        // cxx L519-537.
        if self.my_sol_count < self.my_min_cell_filter_sol {
            for i in 0..self.my_sol_count {
                let mut is_same = true;
                for j in 1..=self.my_n {
                    let idx = i * (self.my_n as usize) + (j - 1) as usize;
                    if (the_pnt.get(j) - self.my_y[idx]).abs() > a_tol.get(j) {
                        is_same = false;
                        break;
                    }
                }
                if is_same {
                    return true;
                }
            }
        } else {
            // cxx L540-570.
            let mut inspector = GlobOptInspector::new(self.my_n as usize, PCONFUSION);
            if self.is_first_cell_filter_invoke {
                self.my_filter.reset_array(&self.my_cell_size);

                // Copy the initial data into the cell filter.
                for a_sol_idx in 0..self.my_sol_count {
                    let mut a_vec = Vector::new(1, self.my_n);
                    for a_sol_dim in 1..=self.my_n {
                        let idx = a_sol_idx * (self.my_n as usize) + (a_sol_dim - 1) as usize;
                        a_vec.set(a_sol_dim, self.my_y[idx]);
                    }
                    let v_clone = a_vec.clone();
                    self.my_filter.add(&v_clone, &a_vec, &GlobOptInspector::coord);
                }
            }
            self.is_first_cell_filter_invoke = false;

            let cell_size = self.my_cell_size.clone();
            let (a_low, an_up) = inspector.shift(the_pnt, &cell_size);

            inspector.clear_find();
            inspector.set_current(the_pnt);
            let mut found = false;
            {
                let mut inspector_ref = &mut inspector;
                self.my_filter.inspect_range(
                    &a_low,
                    &an_up,
                    &GlobOptInspector::coord,
                    &mut |obj: &Vector| {
                        let r = inspector_ref.inspect(obj);
                        if inspector_ref.is_find() {
                            found = true;
                        }
                        r
                    },
                );
            }
            if !found {
                // Point is out of close cells, add a new one.
                let p_clone = the_pnt.clone();
                self.my_filter.add(&p_clone, the_pnt, &GlobOptInspector::coord);
            }
        }
        false
    }

    /// OCCT Points(theIndex, theSol) (cxx L577-585) — 1-based.
    pub fn points(&self, the_index: usize, the_sol: &mut Vector) {
        for j in 1..=self.my_n {
            let idx = (the_index - 1) * (self.my_n as usize) + (j - 1) as usize;
            the_sol.set(j, self.my_y[idx]);
        }
    }

    /// OCCT initCellSize() (cxx L589-596).
    fn init_cell_size(&mut self) {
        for an_idx in 1..=self.my_n {
            self.my_cell_size[(an_idx - 1) as usize] = (self.my_glob_b.get(an_idx)
                - self.my_glob_a.get(an_idx))
                * PCONFUSION
                / (2.0 * M_SQRT2);
        }
    }

    /// OCCT CheckFunctionalStopCriteria() (cxx L600-604).
    fn check_functional_stop_criteria(&self) -> bool {
        self.my_is_find_single_solution
            && (self.my_f - self.my_functional_minimal_value).abs() < self.my_same_tol * 0.01
    }

    /// OCCT ComputeInitSol() (cxx L608-630).
    fn compute_init_sol(&mut self) {
        // cxx L610-618.
        let a_pnt = vec_scale(&vec_add(&self.my_glob_a, &self.my_glob_b), 0.5);
        let mut a_val = 0.0;
        self.my_func.value(&a_pnt, &mut a_val);
        self.check_add_candidate(&a_pnt, a_val);

        // cxx L620-629.
        for i in 1..=3i32 {
            let mut a_pnt = vec_add(
                &self.my_a,
                &vec_scale(&vec_sub(&self.my_b, &self.my_a), (i - 1) as f64 / 2.0),
            );
            let mut a_val = 0.0;
            let pnt_clone = a_pnt.clone();
            if self.compute_local_extremum(&pnt_clone, &mut a_val, &mut a_pnt) {
                self.check_add_candidate(&a_pnt, a_val);
            }
        }
    }

    /// OCCT checkAddCandidate(thePnt, theValue) (cxx L634-669).
    fn check_add_candidate(&mut self, the_pnt: &Vector, the_value: f64) {
        // cxx L636-651.
        if (the_value - self.my_f).abs() < self.my_same_tol * 0.01 && !self.my_is_find_single_solution
        {
            if !self.is_stored(the_pnt) {
                if (the_value - self.my_f) * self.my_z > 0.0 {
                    self.my_f = the_value;
                }
                for j in 1..=self.my_n {
                    self.my_y.push(the_pnt.get(j));
                }
                self.my_sol_count += 1;
            }
        }

        // cxx L653-668.
        let a_delta = (the_value - self.my_f) * self.my_z;
        if a_delta > self.my_same_tol * 0.01
            || (a_delta > 0.0 && self.my_is_find_single_solution)
        {
            self.my_f = the_value;
            self.my_y.clear();
            for j in 1..=self.my_n {
                self.my_y.push(the_pnt.get(j));
            }
            self.my_sol_count = 1;
            self.is_first_cell_filter_invoke = true;
        }
    }

    /// OCCT SetContinuity (hxx L102).
    pub fn set_continuity(&mut self, the_cont: i32) {
        self.my_cont = the_cont;
    }

    /// OCCT GetContinuity (hxx L104).
    pub fn get_continuity(&self) -> i32 {
        self.my_cont
    }

    /// OCCT SetFunctionalMinimalValue (hxx L107-110).
    pub fn set_functional_minimal_value(&mut self, the_minimal_value: f64) {
        self.my_functional_minimal_value = the_minimal_value;
    }

    /// OCCT GetFunctionalMinimalValue (hxx L112).
    pub fn get_functional_minimal_value(&self) -> f64 {
        self.my_functional_minimal_value
    }

    /// OCCT SetLipConstState (hxx L116).
    pub fn set_lip_const_state(&mut self, the_flag: bool) {
        self.my_is_const_locked = the_flag;
    }

    /// OCCT GetLipConstState (hxx L118).
    pub fn get_lip_const_state(&self) -> bool {
        self.my_is_const_locked
    }

    /// OCCT isDone() (hxx L121).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT GetF() (hxx L124).
    pub fn get_f(&self) -> f64 {
        self.my_f
    }

    /// OCCT NbExtrema() (hxx L127).
    pub fn nb_extrema(&self) -> usize {
        self.my_sol_count
    }
}
