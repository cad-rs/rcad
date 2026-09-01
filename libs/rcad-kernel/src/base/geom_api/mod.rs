//! OCCT GeomAPI: high-level geometry algorithms.
//!
//! Sub-modules:
//! - project: GeomAPI_ProjectPointOnCurve + GeomAPI_ProjectPointOnSurf
//! - interpolate: GeomAPI_Interpolate + GeomAPI_PointsToBSpline
//! - extrema: GeomAPI_ExtremaCurveCurve

pub mod interpolate;
pub mod project;
pub mod extrema;
pub mod int_cs;
pub mod int_ss;
pub mod geom2d_interpolate;

pub use geom2d_interpolate::Geom2dInterpolate;
pub use interpolate::{approximate_points, interpolate_points, interpolate_points_2d};
pub use project::{
    closest_point_on_curve, closest_point_on_curve_range, closest_point_on_surface,
    closest_point_on_surface_near, make_pcurve_on_surface, CurveProjection, SurfaceProjection,
};
pub use extrema::{CurveCurveExtrema, ExtremaPair, extrema_curve_curve};
