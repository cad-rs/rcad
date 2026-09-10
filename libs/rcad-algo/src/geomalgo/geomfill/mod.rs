//! OCCT GeomFill package (TKGeomAlgo/GeomFill) — 1:1 port (in progress).
//!
//! Complete: GeomFill_Filling base + GeomFill_Stretch / GeomFill_Coons /
//! GeomFill_Curved + GeomFill_BSplineCurves (anchor:
//! GeomFill_BSplineCurves_Test.cxx / OCC28131 boundary setup); the
//! trihedron-law family (Fixed / ConstantBiNormal / Darboux /
//! TrihedronWithGuide / GuideTrihedronAC / GuideTrihedronPlan); the
//! location-law family (LocationLaw / CurveAndTrihedron / PlanFunc); the
//! section laws (UniformSection); the sweep machinery (SweepSectionGenerator
//! / CircularBlendFunc / Approx_SweepFunction).
//! Remaining (later units): DiscreteTrihedron / NSections consumers of the
//! sweep / Gordon.

pub mod approx_curvlin_func;
pub mod approx_sweep_function;
pub mod bspline_curves;
pub mod circular_blend_func;
pub mod coons;
pub mod constant_bi_normal;
pub mod corrected_frenet;
pub mod curved;
pub mod curve_and_trihedron;
pub mod darboux;
pub mod filling;
pub mod fixed;
pub mod frenet;
pub mod function_guide;
pub mod geom_fill;
pub mod gp_mat;
pub mod guide_trihedron_ac;
pub mod guide_trihedron_plan;
pub mod int_curve_surface_h_inter;
pub mod line;
pub mod location_guide;
pub mod location_law;
pub mod plan_func;
pub mod pipe;
pub mod polynomial_convertor;
pub mod profiler;
pub mod quasi_angular_convertor;
pub mod nsections;
pub mod section_generator;
pub mod section_law;
pub mod section_placement;
pub mod stretch;
pub mod sweep;
pub mod sweep_function;
pub mod sweep_section_generator;
pub mod sngrl_func;
pub mod trihedron_law;
pub mod trihedron_with_guide;
pub mod uniform_section;

pub use bspline_curves::{BSplineCurves, FillingStyle};
pub use trihedron_law::{PipeError, TrihedronLaw, TrihedronLawBase};
pub use coons::Coons;
pub use constant_bi_normal::ConstantBiNormal;
pub use corrected_frenet::{CorrectedFrenet, Trihedron};
pub use curved::Curved;
pub use filling::FillingBase;
pub use fixed::Fixed;
pub use frenet::Frenet;
pub use nsections::NSections;
pub use section_law::SectionLaw;
pub use sngrl_func::SnglrFunc;
pub use stretch::Stretch;
