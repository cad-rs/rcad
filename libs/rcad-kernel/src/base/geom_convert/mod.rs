//! 3D curve conversion utilities (GeomConvert).
//!
//! OCCT TKGeomBase GeomConvert package
//! (`src/ModelingData/TKGeomBase/GeomConvert/`).
//!
//! Parallel to Geom2dConvert but for 3D curves (Geom_*).  Provides the
//! analytic-curve-to-BSpline approximation framework.

pub mod approx_curve;

pub use approx_curve::GeomConvertApproxCurve;
