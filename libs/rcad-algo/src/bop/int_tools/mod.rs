//! OCCT IntTools — geometric intersection algorithms.

pub const CHORD_TOLERANCE: f64 = 1e-7;
/// Max refinement depth for chord-based adaptive sampling.
pub const CHORD_REFINE_DEPTH: usize = 8;

pub mod common_prt;
pub mod curve;
pub mod range;
pub mod root;
pub mod pnt_on_face;
pub mod pnt_on_2_faces;
pub mod int_surf_quadric;
pub mod geom_abs_surface_type;
pub mod int_patch_line;
pub mod int_patch_point;
pub mod int_patch_type;
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
pub mod ellipse_intersection;
pub mod hyperbola_intersection;
pub mod parabola_intersection;
pub mod pcurve_derive;
pub mod edge_edge;
pub mod edge_face;
pub mod face_face;

pub use cylinder_torus::CylinderTorusResult;
pub use plane_torus::PlaneTorusResult;
