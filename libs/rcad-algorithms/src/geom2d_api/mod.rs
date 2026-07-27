//! 2D curve intersection, tangent circles, BSpline fitting, and distance/angle queries.

use crate::tolerance::*;
use glam::DVec2;
use rcad_kernel::geom::*;
use std::f64::consts::PI;

// Circle/curve tangent geometry (~850 lines)
pub mod tangent;
pub use self::tangent::*;
// 2D curve intersection (~75 lines)
pub mod intersect;
pub use self::intersect::*;
// BSpline2D creation (~55 lines)
pub mod bspline;
pub use self::bspline::*;
// Projection, distance, angle queries (~160 lines)
pub mod query;
pub use self::query::*;
// 3D curve to 2D on plane projection (~40 lines)
mod project;
pub use self::project::*;
// Internal helpers (~355 lines)
pub(crate) mod helpers;
