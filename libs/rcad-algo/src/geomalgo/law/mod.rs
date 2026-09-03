//! OCCT Law package (TKGeomAlgo/Law) — 1:1 port (in progress).
//!
//! Complete: Law_Function (base), Law_Constant, Law_Composite.
//! Remaining (needed by GeomFill_CorrectedFrenet / GeomFill_Frenet unit):
//! Law_BSpline (Law_BSpline.cxx L26-1768), Law_BSpFunc, Law_BSplineKnotSplitting,
//! Law_Interpolate (plus the standalone Law_Linear / Law_S).

pub mod law_bsp_func;
pub mod law_bspline;
pub mod law_bspline_knot_splitting;
pub mod law_composite;
pub mod law_constant;
pub mod law_function;
pub mod law_interpolate;

pub use law_bsp_func::LawBSpFunc;
pub use law_bspline::LawBSpline;
pub use law_bspline_knot_splitting::LawBSplineKnotSplitting;
pub use law_composite::LawComposite;
pub use law_constant::LawConstant;
pub use law_function::{LawFunction, LawFunctionHandle};
pub use law_interpolate::LawInterpolate;
