//! OCCT GProp: global shape properties (surface area, volume, inertia).
//!
//! Sub-modules:
//! - tri: triangulation helpers (earcut, UV grid, winding numbers)
//! - surface: BRepGProp::SurfaceProperties (surface area)
//! - volume: BRepGProp::VolumeProperties (volume, centroid)
//! - inertia: GProp_PrincipalProps (inertia tensor)
//! - plate: thin-plate spline surface reconstruction

pub mod tri;
pub mod surface;
pub mod volume;
pub mod inertia;
pub mod plate;

// Re-export public API
pub use tri::{
    face_flat_iter, face_triangles_pub, point_in_spherical_polygon_3d_pub,
    sample_wire_polyline_3d, sample_wire_polyline_3d_with_n,
    trim_almost_closed_polyline,
};
pub use surface::{
    face_surface_area, surface_area, try_analytic_face_surface_area_pub,
    try_spherical_uv_masked_raster,
};
pub use volume::{centroid, signed_volume, volume};
pub use inertia::{InertiaTensor, inertia_tensor};
pub use plate::*;
