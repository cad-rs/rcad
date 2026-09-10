//! OCCT AdvApprox package (ModelingData/TKG3d) — function approximation.
//!
//! 1:1 translation of:
//! - `AdvApprox_EvaluatorFunction.hxx` -> [`EvaluatorFunction`] trait
//! - `AdvApprox_Cutting.hxx` + `AdvApprox_DichoCutting` -> [`Cutting`] /
//!   [`DichoCutting`]
//! - `AdvApprox_SimpleApprox.cxx` -> [`SimpleApprox`]
//! - `AdvApprox_ApproxAFunction.cxx` -> [`ApproxAFunction`] (including the
//!   static `PrepareConvert` and `Approximation` routines)
//!
//! The engine approximates a function by piecewise polynomials in a Jacobi
//! basis (`PLib_JacobiPolynomial`), cuts intervals with a dichotomy until the
//! tolerance is met, and converts the result to BSpline poles via
//! `Convert_CompPolynomialToPoles`.

mod approx_a_function;
mod cutting;
mod evaluator_function;
mod simple_approx;

pub use approx_a_function::{approximation, ApproxAFunction};
pub use cutting::{Cutting, DichoCutting};
pub use evaluator_function::EvaluatorFunction;
pub use simple_approx::SimpleApprox;

use super::GeomAbsShape;

/// OCCT GeomAbs continuity -> ContinuityOrder used by the engine.
pub(crate) fn continuity_order(continuity: GeomAbsShape) -> i32 {
    match continuity {
        GeomAbsShape::C0 => 0,
        GeomAbsShape::C1 => 1,
        GeomAbsShape::C2 => 2,
        _ => panic!("AdvApprox: invalid Continuity"),
    }
}
