//! OCCT GCPnts: points on curves.
//!
//! Sub-modules:
//! - abscissa_point: GCPnts_AbscissaPoint (arc length computation)

pub mod abscissa_point;

pub use abscissa_point::arc_length;
