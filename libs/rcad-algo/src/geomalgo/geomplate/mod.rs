//! OCCT GeomPlate package (TKGeomAlgo/GeomPlate) — 1:1 port (in progress).
//!
//! Point-constraint path complete (anchor: GeomPlate_BuildPlateSurface_Test).
//! Curve-constraint machinery (ProjectCurve/Approx_CurveOnSurface/
//! Geom2dInt_GInter/Discretise/LoadCurve/Intersect) is anchor-out-of-scope
//! and follows the ThruSections precedent: API skeleton, backfill later.

pub mod point_constraint;
pub mod surface;

pub use point_constraint::PointConstraint;
pub use surface::GeomPlateSurface;
