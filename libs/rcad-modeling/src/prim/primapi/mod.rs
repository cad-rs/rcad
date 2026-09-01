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
pub mod make_planar;
pub mod make_half_space;
pub mod make_edge;
pub mod make_wire;
pub mod make_face;
pub mod transform;
pub mod make_polygon;

pub use make_box::{MakeBox, box_brep, make_box_brep};
pub use make_cylinder::{MakeCylinder, cylinder_brep, make_cylinder_brep, prism_face_solid_brep};
pub use make_cone::{MakeCone, cone_brep, make_cone_brep};
pub use make_sphere::{MakeSphere, sphere_brep, make_sphere_brep};
pub use make_torus::{MakeTorus, torus_brep, make_torus_brep};
pub use make_planar::{make_planar_polygon_brep, make_planar_rect_brep};
pub use make_half_space::make_half_space_brep;
pub use make_prism::{make_prism_brep, make_prism_from_face_brep, prism_brep};
pub use make_edge::{
    edge_count_brep, make_edge_brep, make_edge_circle_brep, make_edge_circle_range_brep,
    make_edge_line_brep, make_edge_line_range_brep,
};
pub use make_wire::MakeWire;
pub use make_face::{
    make_face_cylinder_bounds_brep, make_face_from_wire_brep, make_face_plane_brep,
    make_face_plane_bounds_brep,
};
pub use transform::transform_brep;
pub use make_polygon::MakePolygon;
