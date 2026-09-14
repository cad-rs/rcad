//! OCCT Hermit package (TKGeomBase/Hermit) — 1D Hermite reparameterization
//! splines for rational BSpline C1 concatenation, plus the shared BSplCLib
//! evaluation kernels (scalar/point D0/D1, Bohm, RationalDerivative) that the
//! package and `Geom2d_BSplineCurve::EvalD0/EvalD1` consume.  See the child
//! module docs for the OCCT line anchors and the hosting rationale.

pub mod bspl_eval_kernels;
pub mod hermit_solution;

pub use hermit_solution::solution as hermit_solution_2d;
