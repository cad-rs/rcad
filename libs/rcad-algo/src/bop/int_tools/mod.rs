//! OCCT IntTools — geometric intersection algorithms.

pub const CHORD_TOLERANCE: f64 = 1e-7;
pub const CHORD_REFINE_DEPTH: usize = 8;

pub mod common_prt;
pub mod curve;
pub mod range;
pub mod root;
pub mod pnt_on_face;
pub mod pnt_on_2_faces;
pub mod edge_edge;
pub mod edge_face;
pub mod face_face;
pub mod curve_range;
pub mod context;
pub mod fclass2d;
pub mod intss;
pub mod marching;
pub mod pcurve_derive;
pub mod plane_plane;
pub mod plane_cylinder;
pub mod plane_cone;
pub mod plane_sphere;
pub mod plane_torus;
pub mod cylinder_cylinder;
pub mod cylinder_cone;
pub mod cylinder_torus;
pub mod sphere_cylinder;
pub mod sphere_cone;
pub mod sphere_torus;
pub mod cone_cone;
pub mod torus_cone;
pub mod torus_torus;

pub use cylinder_torus::CylinderTorusResult;
pub use plane_torus::PlaneTorusResult;
