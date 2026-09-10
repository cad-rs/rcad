//! OCCT ShapeBuild package (TKShHealing) — 1:1 translation of the four
//! package classes (`ShapeBuild.hxx` L17-40, `ShapeBuild_Edge.hxx` L17-156,
//! `ShapeBuild_ReShape.hxx` L17-120, `ShapeBuild_Vertex.hxx` L17-53):
//! the ReShape substitution engine, the edge construction helpers, the
//! vertex combination tool and the package namespace class, plus the
//! kernel-side TopoDS/TopExp/BRep_Tool machinery they call.

pub mod brep_tool;
pub mod edge;
pub mod reshape;
pub mod vertex;

pub use brep_tool::{
    brep_tool_is_closed, builder_add, iter_subshapes, occt_is_partner, occt_is_same,
    set_flag_inplace, shape_is_null, topexp_explorer,
};
pub use edge::ShapeBuildEdge;
pub use reshape::ShapeBuildReShape;
pub use vertex::ShapeBuildVertex;

/// OCCT ShapeBuild package namespace class (`ShapeBuild.hxx` L28-38).
#[derive(Debug, Clone, Copy, Default)]
pub struct ShapeBuild;

impl ShapeBuild {
    /// OCCT ShapeBuild::PlaneXOY (ShapeBuild.cxx L23-31): returns a
    /// Geom_Surface which is the Plane XOY (Z positive). This allows to
    /// consider an UV space homologous to a 3D space, with this support
    /// surface.
    ///
    /// Architecture note: OCCT caches the plane in a function-local static
    /// handle (`xoy.IsNull()` lazy init); rcad's `Plane` is a plain value, so
    /// the cache has no identity to preserve and the construction is returned
    /// directly. `Geom_Plane(0, 0, 1, 0)` is the plane z = 0 with the +Z
    /// normal; `Plane::new` derives the same gp_Ax3-convention u/v axes.
    pub fn plane_xoy() -> rcad_kernel::geom::Plane {
        rcad_kernel::geom::Plane::new(glam::DVec3::ZERO, glam::DVec3::Z)
    }
}
