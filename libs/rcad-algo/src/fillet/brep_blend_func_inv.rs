//! OCCT Blend_FuncInv (TKFillet/Blend) — 1:1 port of Blend_FuncInv.hxx
//! (L42-92) + Blend_FuncInv.cxx (L18-21); and Blend_CurvPointFuncInv —
//! 1:1 port of Blend_CurvPointFuncInv.hxx (L39-85) +
//! Blend_CurvPointFuncInv.cxx (L20-24).
//!
//! Architecture mapping: both deferred classes derive from
//! `math_FunctionSetWithDerivatives`, expressed as Rust subtraits of
//! [`FunctionSetWithDerivatives`].  `occ::handle<Adaptor2d_Curve2d>` maps to
//! a `&Curve2d` reference (rcad `Curve2d` enum).

use glam::DVec3;

use rcad_kernel::geom::Curve2d;
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

/// OCCT Blend_FuncInv — deferred class for a function used to compute a
/// blending surface between two surfaces, using a guide line; this function
/// is used to find a solution on a restriction of one of the surfaces
/// (Blend_FuncInv.hxx L42).  The vector X is t, w, U, V.
pub trait BlendFuncInv: FunctionSetWithDerivatives {
    /// OCCT NbVariables() (Blend_FuncInv.cxx L18-21) — returns 4.
    fn nb_variables(&self) -> usize {
        4
    }

    /// OCCT Set(OnFirst, COnSurf) — sets the CurveOnSurface on which a
    /// solution has to be found; if OnFirst is true the curve is on the
    /// first surface, otherwise on the second one.
    fn set_curve_on_surface(&mut self, on_first: bool, c_on_surf: &Curve2d);

    /// OCCT GetTolerance(Tolerance, Tol) — the parametric tolerance for each
    /// of the 4 variables; Tol is the tolerance used in 3d space.
    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64);

    /// OCCT GetBounds(InfBound, SupBound) — the lowest / greatest values
    /// allowed for each of the 4 variables.
    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]);

    /// OCCT IsSolution(Sol, Tol) — true if Sol is a zero of the function.
    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool;
}

/// OCCT Blend_CurvPointFuncInv — deferred class for a function used to
/// compute a blending surface between a surface and a curve, using a guide
/// line; this function is used to find a solution on a done point of the
/// curve (Blend_CurvPointFuncInv.hxx L39).  The vector X is w, U, V.
pub trait BlendCurvPointFuncInv: FunctionSetWithDerivatives {
    /// OCCT NbVariables() (Blend_CurvPointFuncInv.cxx L21-24) — returns 2
    /// (sic: the OCCT body returns 2 although the header comment says 3;
    /// translated literally).
    fn nb_variables(&self) -> usize {
        2
    }

    /// OCCT Set(P) — sets the Point on which a solution has to be found.
    fn set_point(&mut self, p: DVec3);

    /// OCCT GetTolerance(Tolerance, Tol) — the parametric tolerance for the
    /// variables; Tol is the tolerance used in 3d space.
    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64);

    /// OCCT GetBounds(InfBound, SupBound) — the lowest / greatest values
    /// allowed for the variables.
    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]);

    /// OCCT IsSolution(Sol, Tol) — true if Sol is a zero of the function.
    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool;
}

/// OCCT Blend_SurfPointFuncInv — deferred class for a function used to
/// compute a blending surface between a surface and a point on a curve,
/// finding a solution on the surface (Blend_SurfPointFuncInv.hxx L25-60).
/// The vector X is U, V.
pub trait BlendSurfPointFuncInv: FunctionSetWithDerivatives {
    /// OCCT NbVariables() (Blend_SurfPointFuncInv.cxx L19-22) — returns 2.
    fn nb_variables(&self) -> usize {
        2
    }

    /// OCCT NbEquations() — the number of equations of the function.
    fn nb_equations(&self) -> usize;

    /// OCCT Set(P) — sets the Point on which a solution has to be found.
    fn set_point(&mut self, p: DVec3);

    /// OCCT GetTolerance(Tolerance, Tol) — the parametric tolerance for the
    /// variables; Tol is the tolerance used in 3d space.
    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64);

    /// OCCT GetBounds(InfBound, SupBound) — the lowest / greatest values
    /// allowed for the variables.
    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]);

    /// OCCT IsSolution(Sol, Tol) — true if Sol is a zero of the function.
    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool;
}

/// OCCT Blend_SurfCurvFuncInv — deferred class for a function used to compute
/// a blending surface between a surface and a curve, finding a solution on a
/// restriction of the surface (Blend_SurfCurvFuncInv.hxx L26-70).
/// The vector X is t, U, V.
pub trait BlendSurfCurvFuncInv: FunctionSetWithDerivatives {
    /// OCCT NbVariables() (Blend_SurfCurvFuncInv.cxx L19-22) — returns 3.
    fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations() — the number of equations of the function.
    fn nb_equations(&self) -> usize;

    /// OCCT Set(Rst) — sets the restriction on which a solution has to be
    /// found.
    fn set_rst(&mut self, rst: &Curve2d);

    /// OCCT GetTolerance(Tolerance, Tol) — the parametric tolerance for the
    /// variables; Tol is the tolerance used in 3d space.
    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64);

    /// OCCT GetBounds(InfBound, SupBound) — the lowest / greatest values
    /// allowed for the variables.
    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]);

    /// OCCT IsSolution(Sol, Tol) — true if Sol is a zero of the function.
    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool;
}
