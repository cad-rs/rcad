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
pub mod linear;
pub mod plate;

// Re-export public API
pub use tri::{
    face_flat_iter, face_triangles_pub, point_in_spherical_polygon_3d_pub,
    sample_wire_polyline_3d, sample_wire_polyline_3d_with_n,
    trim_almost_closed_polyline,
};
pub use surface::{face_surface_area, surface_area};
pub use volume::{
    centroid, face_volume_gauss_domain, face_volume_gauss_domain_full,
    face_volume_gauss_natural, face_volume_gauss_natural_full, shape_vinert, signed_volume,
    VinertFace, volume,
};
pub use inertia::{InertiaTensor, PrincipalProps, inertia_tensor, principal_properties};
pub use linear::linear_properties;
pub use plate::*;
