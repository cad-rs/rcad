//! OCCT GeomFill package (TKGeomAlgo/GeomFill) — 1:1 port (in progress).
//!
//! Complete: GeomFill_Filling base + GeomFill_Stretch / GeomFill_Coons /
//! GeomFill_Curved + GeomFill_BSplineCurves (anchor:
//! GeomFill_BSplineCurves_Test.cxx / OCC28131 boundary setup).
//! Remaining (later units): CorrectedFrenet / NSections / Gordon /
//! GuideTrihedronAC and the sweep machinery.

pub mod bspline_curves;
pub mod coons;
pub mod curved;
pub mod filling;
pub mod stretch;

pub use bspline_curves::{BSplineCurves, FillingStyle};
pub use coons::Coons;
pub use curved::Curved;
pub use filling::FillingBase;
pub use stretch::Stretch;
