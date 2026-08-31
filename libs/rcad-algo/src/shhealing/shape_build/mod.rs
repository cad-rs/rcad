//! OCCT ShapeBuild package (TKShHealing) — foundation classes for the
//! ShapeFix stack: the ReShape substitution engine and edge construction
//! helpers, plus the kernel-side TopoDS/TopExp/BRep_Tool machinery they call.

pub mod brep_tool;
pub mod edge;
pub mod reshape;

pub use brep_tool::{
    brep_tool_is_closed, builder_add, iter_subshapes, occt_is_partner, occt_is_same,
    set_flag_inplace, shape_is_null, topexp_explorer,
};
pub use edge::ShapeBuildEdge;
pub use reshape::ShapeBuildReShape;
