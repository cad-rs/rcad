pub mod coplanar;
pub mod curve_surface;
pub mod edge_face;
pub mod intss;
pub mod marching;
pub mod pcurve_derive;
pub mod plane_cone;
pub mod plane_cylinder;
pub mod plane_plane;
pub mod plane_sphere;
pub mod vertex_ops;

pub use intss::{
    SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection, intersect_surfaces,
};
