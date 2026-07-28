//! OCCT GProp: global shape properties (surface area, volume, inertia).
//!
//! Sub-modules:
//! - plate: thin-plate spline surface reconstruction
//! - surface, volume, inertia, tri: (in progress — currently forwarded from properties.rs)
//!
//! TODO: Split properties.rs into these sub-modules for OCCT alignment.

pub mod plate;

/// Re-export from the existing properties.rs (migration path).
/// All surface_area, volume, centroid, inertia_tensor functions are here.
pub use crate::math::properties::*;
