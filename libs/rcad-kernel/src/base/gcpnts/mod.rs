//! OCCT GCPnts: points on curves.
//!
//! Sub-modules:
//! - abscissa_point: GCPnts_AbscissaPoint (arc length computation)
//! - gcpnts_curve: the GCPnts adaptor-interface projection
//!   (Adaptor3d_Curve / Adaptor2d_Curve2d shims)
//! - gcpnts_curve_bridge: the Extrema_CurveTool-facade adapter onto
//!   `GCPntsCurve` (DeflCurvIntervals cxx L83)
//! - gcpnts_tangential_deflection: GCPnts_TangentialDeflection +
//!   GCPnts_DistFunction

pub mod abscissa_point;
pub mod gcpnts_curve;
pub mod gcpnts_curve_bridge;
pub mod gcpnts_tangential_deflection;

pub use abscissa_point::arc_length;
