//! OCCT Adaptor3d_InterFunc (TKG3d/Adaptor3d) — the function U(t) = U0 or
//! V(t) = V0 used to find the points where a 2D pcurve crosses a surface
//! isoparametric discontinuity, in order to determine the Cn discontinuities
//! of an Adaptor3d_CurveOnSurface relative to the discontinuities of the
//! surface.
//!
//! 1:1 translation of Adaptor3d_InterFunc.hxx / Adaptor3d_InterFunc.cxx.
//! `math_FunctionWithDerivative` is the rcad trait of the same shape in
//! `crate::math::root::function_all_roots`; `occ::handle<Adaptor2d_Curve2d>`
//! maps to [`Curve2dHandle`] (shared ownership).

use super::adaptor::Curve2dHandle;
use crate::math::root::function_all_roots::{FunctionValue, FunctionWithDerivative};

/// OCCT Adaptor3d_InterFunc (Adaptor3d_InterFunc.hxx L29-61).
pub struct Adaptor3dInterFunc {
    /// OCCT: handle(Adaptor2d_Curve2d) myCurve2d.
    pub my_curve2d: Curve2dHandle,
    /// OCCT: Standard_Real myFixVal.
    pub my_fix_val: f64,
    /// OCCT: Standard_Integer myFix (1 = fix X, 2 = fix Y).
    pub my_fix: i32,
}

impl Adaptor3dInterFunc {
    /// OCCT Adaptor3d_InterFunc::Adaptor3d_InterFunc(C, FixVal, Fix)
    /// (Adaptor3d_InterFunc.cxx L24-35): builds the function U(t) = FixVal if
    /// Fix = 1 or V(t) = FixVal if Fix = 2; raises Standard_ConstructionError
    /// for any other Fix.
    pub fn new(c: Curve2dHandle, fix_val: f64, fix: i32) -> Self {
        if fix != 1 && fix != 2 {
            panic!("Standard_ConstructionError: Adaptor3d_InterFunc");
        }
        Adaptor3dInterFunc {
            my_curve2d: c,
            my_fix_val: fix_val,
            my_fix: fix,
        }
    }
}

impl FunctionValue for Adaptor3dInterFunc {
    /// OCCT Adaptor3d_InterFunc::Value(X, F)
    /// (Adaptor3d_InterFunc.cxx L37-51): F = C.X() - myFixVal when myFix == 1,
    /// F = C.Y() - myFixVal otherwise; always succeeds.
    fn value(&mut self, x: f64) -> Option<f64> {
        let c = self.my_curve2d.d0(x);
        if self.my_fix == 1 {
            Some(c.x - self.my_fix_val)
        } else {
            Some(c.y - self.my_fix_val)
        }
    }
}

impl FunctionWithDerivative for Adaptor3dInterFunc {
    /// OCCT Adaptor3d_InterFunc::Derivative(X, D)
    /// (Adaptor3d_InterFunc.cxx L53-57): delegates to Values.
    fn derivative(&mut self, x: f64) -> Option<f64> {
        // OCCT: double F; return Values(X, F, D);
        self.values(x).map(|(_f, d)| d)
    }

    /// OCCT Adaptor3d_InterFunc::Values(X, F, D)
    /// (Adaptor3d_InterFunc.cxx L59-75): F and D are the coordinate
    /// (myFix == 1: X, else Y) of the curve point and derivative; always
    /// succeeds.
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let (c, dc) = self.my_curve2d.d1(x);
        if self.my_fix == 1 {
            Some((c.x - self.my_fix_val, dc.x))
        } else {
            Some((c.y - self.my_fix_val, dc.y))
        }
    }
}

// OCCT `occ::handle<Adaptor3d_InterFunc>` has no consumer in the rcad
// encoding; the function is passed by reference to math_FunctionRoots
// (Adaptor3d_CurveOnSurface.cxx L1088, L1099), which maps to
// `&mut dyn FunctionWithDerivative`.
