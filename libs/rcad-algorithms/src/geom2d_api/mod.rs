//! 2D curve intersection, tangent circles, BSpline fitting, and distance/angle queries.

use rcad_kernel::geom::*;
use std::f64::consts::PI;
use glam::DVec2;
use crate::tolerance::*;

// Circle/curve tangent geometry (~850 lines)
include!("tangent_inc.rs");
// 2D curve intersection (~75 lines)
include!("intersect_inc.rs");
// BSpline2D creation (~55 lines)
include!("bspline_inc.rs");
// Projection, distance, angle queries (~160 lines)
include!("query_inc.rs");
// Internal helpers (~355 lines)
include!("helpers_inc.rs");

