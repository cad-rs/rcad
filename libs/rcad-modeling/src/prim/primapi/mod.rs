// OCCT BRepPrimAPI — high-level API for building primitive shapes.
//
// OCCT src: ModelingAlgorithms/TKPrim/BRepPrimAPI
//
// Provides MakeBox, MakeCylinder, MakeCone, MakeSphere, MakeTorus, etc.
// Each type wraps a low-level BRepPrim_* builder and exposes the
// constructed TopoDS_Shape.

pub mod make_box;
pub mod make_cylinder;
pub mod make_cone;
pub mod make_sphere;
pub mod make_torus;
pub mod make_prism;

pub use make_box::{MakeBox, box_brep, make_box_brep};
pub use make_cylinder::{MakeCylinder, cylinder_brep, make_cylinder_brep};
pub use make_cone::{MakeCone, cone_brep, make_cone_brep};
pub use make_sphere::{MakeSphere, sphere_brep, make_sphere_brep};
pub use make_torus::{MakeTorus, torus_brep, make_torus_brep};
pub use make_prism::{make_prism_brep, prism_brep};
