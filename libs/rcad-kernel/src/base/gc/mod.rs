//! Geometric construction algorithms (GC package).
//!
//! Provides validated construction of geometric objects (curves, surfaces,
//! transformations) from various input combinations.  Each constructor returns
//! `Result<T, GceError>` where `GceError` mirrors OCCT `gce_ErrorType`.
//!
//! OCCT TKGeomBase GC package: GC_MakeCircle, GC_MakeLine, GC_MakePlane, etc.

#![allow(clippy::manual_clamp)]

use glam::{DVec2, DVec3};
use std::fmt;

// Re-export all public constructors
mod curves;
mod curves2d;
mod surfaces;
mod transforms;

pub use curves::*;
pub use curves2d::*;
pub use surfaces::*;
pub use transforms::*;

// ============================================================================
// GceError — mirrors OCCT gce_ErrorType
// ============================================================================

/// Error status for geometric construction operations.
///
/// Mirrors OCCT `gce_ErrorType`:
///
/// | Variant | OCCT | Meaning |
/// |---------|------|---------|
/// | `NegativeRadius` | `gce_NegativeRadius` | Radius < 0 |
/// | `NullAxis` | `gce_NullAxis` | Null direction vector for axis |
/// | `ConfusedPoints` | `gce_ConfusedPoints` | Two points are coincident |
/// | `InvertAxis` | `gce_InvertAxis` | Invalid axis orientation |
/// | `BadEquation` | `gce_BadEquation` | Degenerate plane equation |
/// | `ColinearPoints` | `gce_ColinearPoints` | Three points are collinear |
/// | `NegativeLength` | `gce_NegativeLength` | Negative length parameter |
/// | `NullAngle` | `gce_NullAngle` | Zero angle |
/// | `NullRadius` | `gce_NullRadius` | Zero radius outside tolerance |
/// | `SameParameters` | — | Two parameters are the same |
/// | `ZeroDistance` | — | Distance is zero (returned from parallel offset 0) |
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GceError {
    NegativeRadius,
    NullAxis,
    ConfusedPoints,
    InvertAxis,
    BadEquation,
    ColinearPoints,
    NegativeLength,
    NullAngle,
    NullRadius,
    SameParameters,
    NullLength,
    ZeroDistance,
}

impl fmt::Display for GceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GceError::NegativeRadius => write!(f, "negative radius"),
            GceError::NullAxis => write!(f, "null axis"),
            GceError::ConfusedPoints => write!(f, "confused points"),
            GceError::InvertAxis => write!(f, "inverted axis"),
            GceError::BadEquation => write!(f, "bad equation"),
            GceError::ColinearPoints => write!(f, "colinear points"),
            GceError::NegativeLength => write!(f, "negative length"),
            GceError::NullAngle => write!(f, "null angle"),
            GceError::NullRadius => write!(f, "null radius"),
            GceError::SameParameters => write!(f, "same parameters"),
            GceError::NullLength => write!(f, "null length"),
            GceError::ZeroDistance => write!(f, "zero distance"),
        }
    }
}

impl std::error::Error for GceError {}

// ============================================================================
// Internal helpers — shared tolerances
// ============================================================================

/// Angular tolerance used for direction comparison (1e-12 rad).
const TOL_ANG: f64 = 1e-12;
/// Confusion tolerance for point coincidence (1e-7).
const TOL_CONF: f64 = 1e-7;
/// Float deduplication tolerance (1e-15).
const TOL_FLOAT: f64 = 1e-15;

/// Returns `true` when two vectors are parallel (cross product near zero).
fn vectors_parallel(a: DVec3, b: DVec3) -> bool {
    a.cross(b).length_squared() < TOL_ANG * TOL_ANG
}

/// Returns `true` when two 2D vectors are parallel (cross product near zero).
fn vectors_parallel_2d(a: DVec2, b: DVec2) -> bool {
    (a.x * b.y - a.y * b.x).abs() < TOL_ANG
}

/// Returns `true` when two points are coincident.
fn points_coincident(a: DVec3, b: DVec3) -> bool {
    (a - b).length_squared() < TOL_CONF * TOL_CONF
}

/// Returns `true` when two 2D points are coincident.
fn points_coincident_2d(a: DVec2, b: DVec2) -> bool {
    (a - b).length_squared() < TOL_CONF * TOL_CONF
}

/// Stable orthonormal frame from axis and reference direction.
/// Returns (x_dir, y_dir) where y_dir = axis × x_dir.
fn orthonormal_from_ref(axis: DVec3, ref_dir: DVec3) -> (DVec3, DVec3) {
    let axis = axis.normalize_or_zero();
    let mut x_dir = ref_dir - axis * ref_dir.dot(axis);
    if x_dir.length_squared() <= 1e-24 {
        // Fallback: use any_perpendicular
        let alt = if axis.x.abs() > 1.0 - 1e-12 {
            DVec3::Z
        } else {
            DVec3::X
        };
        x_dir = (alt - axis * alt.dot(axis)).normalize_or_zero();
    } else {
        x_dir = x_dir.normalize();
    }
    let y_dir = axis.cross(x_dir).normalize_or_zero();
    (x_dir, y_dir)
}
