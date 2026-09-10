//! OCCT Plate package (TKGeomAlgo/Plate) — 1:1 port.
//!
//! Variational spline algorithm defining a two-variable function satisfying
//! constraints while minimizing an energy-like criterion
//! (Plate_Plate.hxx L48-50).

pub mod constraints;
pub mod d123;
pub mod pinpoint_constraint;
pub mod plate;

pub use constraints::{
    FreeGtoCConstraint, GlobalTranslationConstraint, GtoCConstraint, LineConstraint,
    LinearScalarConstraint, LinearXYZConstraint, PlaneConstraint, SampledCurveConstraint,
};
pub use d123::{PlateD1, PlateD2, PlateD3};
pub use pinpoint_constraint::PinpointConstraint;


pub use plate::Plate;
