//! OCCT GeomFill package (TKGeomAlgo/GeomFill) — 1:1 port (in progress).
//!
//! Complete: GeomFill_Filling base + GeomFill_Stretch / GeomFill_Coons /
//! GeomFill_Curved + GeomFill_BSplineCurves (anchor:
//! GeomFill_BSplineCurves_Test.cxx / OCC28131 boundary setup).
//! Remaining (later units): CorrectedFrenet / NSections / Gordon /
//! GuideTrihedronAC and the sweep machinery.

pub mod bspline_curves;
pub mod coons;
pub mod corrected_frenet;
pub mod curved;
pub mod filling;
pub mod frenet;
pub mod nsections;
pub mod section_law;
pub mod stretch;
pub mod sngrl_func;
pub mod trihedron_law;

pub use bspline_curves::{BSplineCurves, FillingStyle};
pub use trihedron_law::{PipeError, TrihedronLaw, TrihedronLawBase};
pub use coons::Coons;
pub use corrected_frenet::{CorrectedFrenet, Trihedron};
pub use curved::Curved;
pub use filling::FillingBase;
pub use frenet::Frenet;
pub use nsections::NSections;
pub use section_law::SectionLaw;
pub use sngrl_func::SnglrFunc;
pub use stretch::Stretch;
