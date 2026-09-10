//! OCCT ProjLib_PrjFunc + ProjLib_PrjResolve (TKGeomBase/ProjLib).
//!
//! 1:1 translation of:
//! - ProjLib_PrjFunc.hxx L17-76 + ProjLib_PrjFunc.cxx L18-158: the E(u,v) = 0
//!   function set solved by the root finders (E1 = (S(t)-C) . DS1_u,
//!   E2 = (S(t)-C) . DS1_v, scaled by myNorm),
//! - ProjLib_PrjResolve.hxx L17-64 + ProjLib_PrjResolve.cxx L26-170: the
//!   Newton solver driver with the FunctionSetRoot fallback.
//!
//! GAP carrier (staged): OCCT solves first with `math_NewtonFunctionSetRoot`
//! (TKMath) and falls back to `math_FunctionSetRoot`.  rcad has the
//! `math_FunctionSetRoot` translation; the Newton wrapper below models the
//! not-done output so the driver always takes the OCCT fallback path.

use glam::DVec2;

use crate::base::proj_lib::adaptor::{Adaptor3dCurve, Adaptor3dSurface};
use crate::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};

// ---------------------------------------------------------------------------
// math_NewtonFunctionSetRoot — GAP carrier (staged)
// ---------------------------------------------------------------------------

/// OCCT math_NewtonFunctionSetRoot (TKMath/math_NewtonFunctionSetRoot.hxx) —
/// pure-Newton root finder with a functional tolerance.  GAP (staged): the
/// class is not translated; Perform reports IsDone() == false which drives
/// the ProjLib_PrjResolve driver into the OCCT math_FunctionSetRoot fallback.
pub struct NewtonFunctionSetRoot {
    done: bool,
}

impl NewtonFunctionSetRoot {
    /// OCCT math_NewtonFunctionSetRoot(Function, Tolerance, FunctionTol).
    pub fn new(
        _f: &mut dyn FunctionSetWithDerivatives,
        _tolerance: &[f64],
        _function_tol: f64,
    ) -> Self {
        NewtonFunctionSetRoot { done: false }
    }

    /// OCCT Perform(Function, StartingPoint, InfBound, SupBound).
    pub fn perform(
        &mut self,
        _f: &mut dyn FunctionSetWithDerivatives,
        _starting_point: &[f64],
        _inf_bound: &[f64],
        _sup_bound: &[f64],
    ) {
        // GAP: deferred until the math_NewtonFunctionSetRoot body lands.
        self.done = false;
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }
}

// ---------------------------------------------------------------------------
// ProjLib_PrjFunc
// ---------------------------------------------------------------------------

/// OCCT ProjLib_PrjFunc (PrjFunc.hxx L27-76) — the function E(u, v) = 0 whose
/// root is the normal projection of the fixed curve parameter on the surface.
pub struct PrjFunc<'a, 'b> {
    /// OCCT: const Adaptor3d_Curve* myCurve.
    my_curve: &'a dyn Adaptor3dCurve,
    /// OCCT: const Adaptor3d_Surface* mySurface.
    my_surface: &'b dyn Adaptor3dSurface,
    myt: f64,
    my_u: f64,
    my_v: f64,
    /// OCCT: Standard_Integer myFix (1: t fixed; 2: u fixed; 3: v fixed).
    my_fix: i32,
    /// OCCT: Standard_Real myNorm.
    my_norm: f64,
}

impl<'a, 'b> PrjFunc<'a, 'b> {
    /// OCCT ProjLib_PrjFunc::ProjLib_PrjFunc(C, FixVal, S, Fix)
    /// (PrjFunc.cxx L20-39).
    pub fn new(
        c: &'a dyn Adaptor3dCurve,
        fix_val: f64,
        s: &'b dyn Adaptor3dSurface,
        fix: i32,
    ) -> Self {
        let mut f = PrjFunc {
            my_curve: c,
            my_surface: s,
            myt: 0.0,
            my_u: 0.0,
            my_v: 0.0,
            my_fix: fix,
            my_norm: 1.0f64.min(s.u_resolution(1.0)).min(s.v_resolution(1.0)),
        };
        match f.my_fix {
            1 => f.myt = fix_val,
            2 => f.my_u = fix_val,
            3 => f.my_v = fix_val,
            _ => panic!("Standard_ConstructionError"),
        }
        f
    }

    /// OCCT ProjLib_PrjFunc::Solution() (PrjFunc.cxx L145-158).
    pub fn solution(&self) -> DVec2 {
        match self.my_fix {
            1 => DVec2::new(self.my_u, self.my_v),
            2 => DVec2::new(self.myt, self.my_v),
            3 => DVec2::new(self.myt, self.my_u),
            // For NT, even if we never reach this point.
            _ => DVec2::new(0.0, 0.0),
        }
    }
}

impl FunctionSetWithDerivatives for PrjFunc<'_, '_> {
    /// OCCT ProjLib_PrjFunc::NbVariables() (PrjFunc.cxx L41-44).
    fn nb_variables(&self) -> usize {
        2
    }

    /// OCCT ProjLib_PrjFunc::NbEquations() (L46-49).
    fn nb_equations(&self) -> usize {
        2
    }

    /// OCCT ProjLib_PrjFunc::Value(X, F) (L51-57).
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let mut d = vec![vec![0.0f64; 2]; 2];
        self.values(x, f, &mut d)
    }

    /// OCCT ProjLib_PrjFunc::Derivatives(X, D) (L59-64).
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        let mut f = vec![0.0f64; 2];
        self.values(x, &mut f, df)
    }

    /// OCCT ProjLib_PrjFunc::Values(X, F, D) (L66-143).
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        let (t, u, v) = match self.my_fix {
            1 => (self.myt, x[0], x[1]),
            2 => (x[0], self.my_u, x[1]),
            _ => (x[0], x[1], self.my_v),
        };

        let (c, dc1_t) = self.my_curve.d1(t);
        let (s, ds1_u, ds1_v, ds2_u, ds2_v, ds2_uv) = self.my_surface.d2(u, v);

        let big_v = s - c; // OCCT gp_Vec V(C, S).

        f[0] = big_v.dot(ds1_u) * self.my_norm;
        f[1] = big_v.dot(ds1_v) * self.my_norm;

        match self.my_fix {
            1 => {
                df[0][0] = (ds1_u.length_squared() + big_v.dot(ds2_u)) * self.my_norm; // dE1/du
                df[0][1] = (ds1_v.dot(ds1_u) + big_v.dot(ds2_uv)) * self.my_norm; // dE1/dv
                df[1][0] = df[0][1]; // dE2/du
                df[1][1] = (ds1_v.length_squared() + big_v.dot(ds2_v)) * self.my_norm; // dE2/dv
            }
            2 => {
                df[0][0] = (-dc1_t.dot(ds1_u)) * self.my_norm; // dE1/dt
                df[0][1] = (ds1_v.dot(ds1_u) + big_v.dot(ds2_uv)) * self.my_norm; // dE1/dv
                df[1][0] = (-dc1_t.dot(ds1_v)) * self.my_norm; // dE2/dt
                df[1][1] = (ds1_v.length_squared() + big_v.dot(ds2_v)) * self.my_norm; // dE2/dv
            }
            _ => {
                df[0][0] = -dc1_t.dot(ds1_u) * self.my_norm; // dE1/dt
                df[0][1] = (ds1_u.length_squared() + big_v.dot(ds2_u)) * self.my_norm; // dE1/du
                df[1][0] = -dc1_t.dot(ds1_v) * self.my_norm; // dE2/dt
                df[1][1] = (ds1_v.dot(ds1_u) + big_v.dot(ds2_uv)) * self.my_norm; // dE2/du
            }
        }

        self.my_u = u;
        self.my_v = v;
        self.myt = t;

        true
    }
}

// ---------------------------------------------------------------------------
// ProjLib_PrjResolve
// ---------------------------------------------------------------------------

/// OCCT ProjLib_PrjResolve (PrjResolve.hxx L22-64) — search of the projection
/// of a point of parameter t of the curve on the surface, with one of the
/// three (t, u, v) coordinates fixed.
pub struct PrjResolve<'a, 'b> {
    /// OCCT: const Adaptor3d_Curve* myCurve.
    my_curve: &'a dyn Adaptor3dCurve,
    /// OCCT: const Adaptor3d_Surface* mySurface.
    my_surface: &'b dyn Adaptor3dSurface,
    my_done: bool,
    /// OCCT: Standard_Integer myFix.
    my_fix: i32,
    /// OCCT: gp_Pnt2d mySolution.
    my_solution: DVec2,
}

impl<'a, 'b> PrjResolve<'a, 'b> {
    /// OCCT ProjLib_PrjResolve::ProjLib_PrjResolve(C, S, Fix)
    /// (PrjResolve.cxx L26-39).
    pub fn new(c: &'a dyn Adaptor3dCurve, s: &'b dyn Adaptor3dSurface, fix: i32) -> Self {
        if fix > 3 || fix < 1 {
            panic!("Standard_ConstructionError");
        }
        PrjResolve {
            my_curve: c,
            my_surface: s,
            my_done: false,
            my_fix: fix,
            my_solution: DVec2::new(0.0, 0.0),
        }
    }

    /// OCCT ProjLib_PrjResolve::Perform(t, U, V, Tol2d, Inf, Sup, FuncTol,
    /// StrictInside) (PrjResolve.cxx L44-156).  The StrictInside flag is
    /// unused by the OCCT body (kept for signature parity).
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        t: f64,
        u: f64,
        v: f64,
        tol2d: DVec2,
        inf: DVec2,
        sup: DVec2,
        func_tol: f64,
        _strict_inside: bool,
    ) {
        self.my_done = false;
        let mut fix_val = 0.0f64;
        let ext_u = 10.0 * tol2d.x;
        let ext_v = 10.0 * tol2d.y;
        let ext_inf = DVec2::new(inf.x - ext_u, inf.y - ext_v);
        let ext_sup = DVec2::new(sup.x + ext_u, sup.y + ext_v);
        let b_inf = [ext_inf.x, ext_inf.y];
        let b_sup = [ext_sup.x, ext_sup.y];
        let tol = [tol2d.x, tol2d.y];

        let mut start = [0.0f64; 2];
        match self.my_fix {
            1 => {
                start[0] = u;
                start[1] = v;
                fix_val = t;
            }
            2 => {
                start[0] = t;
                start[1] = v;
                fix_val = u;
            }
            3 => {
                start[0] = t;
                start[1] = u;
                fix_val = v;
            }
            _ => {}
        }

        // OCCT L95-96: math_NewtonFunctionSetRoot SR(F, Tol, FuncTol);
        //               SR.Perform(F, Start, BInf, BSup);
        let mut f = PrjFunc::new(self.my_curve, fix_val, self.my_surface, self.my_fix);
        let mut sr = NewtonFunctionSetRoot::new(&mut f, &tol, func_tol);
        sr.perform(&mut f, &start, &b_inf, &b_sup);

        // OCCT L98-107: if (!SR.IsDone()) { math_FunctionSetRoot S1(F, Tol);
        //               S1.Perform(F, Start, BInf, BSup); ... }
        if !sr.is_done() {
            let mut s1 = FunctionSetRoot::new(&mut f, &tol, 100);
            s1.perform(&mut f, &start, &b_inf, &b_sup, false);

            if !s1.is_done() {
                return;
            }
        }

        let sol = f.solution();
        self.my_solution = DVec2::new(sol.x, sol.y);

        // computation of myDone (L111-112).
        self.my_done = true;

        let extra_u = 2.0 * tol2d.x;
        let extra_v = 2.0 * tol2d.y;
        if self.my_solution.x > inf.x - tol2d.x && self.my_solution.x < inf.x {
            self.my_solution.x = inf.x;
        }
        if self.my_solution.x > sup.x && self.my_solution.x < sup.x + tol2d.x {
            self.my_solution.x = sup.x;
        }
        if self.my_solution.y > inf.y - tol2d.y && self.my_solution.y < inf.y {
            self.my_solution.y = inf.y;
        }
        if self.my_solution.y > sup.y && self.my_solution.y < sup.y + tol2d.y {
            self.my_solution.y = sup.y;
        }
        if self.my_solution.x < inf.x - extra_u
            || self.my_solution.x > sup.x + extra_u
            || self.my_solution.y < inf.y - extra_v
            || self.my_solution.y > sup.y + extra_v
        {
            self.my_done = false;
        } else if func_tol > 0.0 {
            // OCCT L142-154: X(1..2) = mySolution; F.Value(X, FVal);
            // if (!SR.IsDone() && |FVal|^2 > FuncTol) myDone = false.
            let x = [self.my_solution.x, self.my_solution.y];
            let mut f_val = vec![0.0f64; 2];
            let mut f2 = PrjFunc::new(self.my_curve, fix_val, self.my_surface, self.my_fix);
            f2.value(&x, &mut f_val);

            if !sr.is_done() && (f_val[0] * f_val[0] + f_val[1] * f_val[1]) > func_tol {
                self.my_done = false;
            }
        }
    }

    /// OCCT ProjLib_PrjResolve::IsDone() (L158-161).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT ProjLib_PrjResolve::Solution() (L163-170) — throws
    /// StdFail_NotDone when not done (mirrored as a panic).
    pub fn solution(&self) -> DVec2 {
        if !self.is_done() {
            panic!("StdFail_NotDone");
        }
        self.my_solution
    }
}
