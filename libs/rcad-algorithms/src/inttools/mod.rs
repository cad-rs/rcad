pub mod coplanar;
pub mod curve_surface;
pub mod cylinder_cylinder;
pub mod edge_face;
pub mod intss;
pub mod marching;
pub mod pcurve_derive;
pub mod plane_cone;
pub mod plane_cylinder;
pub mod plane_plane;
pub mod plane_sphere;
pub mod sphere_cylinder;
pub mod vertex_ops;

pub use intss::{
    SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection, intersect_surfaces,
    intersect_surfaces_with_density,
};
