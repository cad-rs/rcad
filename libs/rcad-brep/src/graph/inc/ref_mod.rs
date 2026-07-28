//! Reference entry structs (typed incidence edges).
//!
//! OCCT BRepGraphInc: BRepGraphInc_Reference.hxx
//!
//! Each reference entry connects a parent entity to a child definition,
//! carrying orientation and (for some kinds) local placement.

use crate::graph::inc::id::*;

/// Common generation counter shared by all reference types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OwnGen(pub u32);

/// Shell → Solid membership.
#[derive(Debug, Clone)]
pub struct ShellRef {
    pub own_gen: OwnGen,
    pub shell_id: ShellId,
    pub solid_id: SolidId,
}

/// Face → Shell membership.
#[derive(Debug, Clone)]
pub struct FaceRef {
    pub own_gen: OwnGen,
    pub face_id: FaceId,
    pub shell_id: ShellId,
    pub orientation: rcad_kernel::topods::Orientation,
}

/// Wire → Face membership.
#[derive(Debug, Clone)]
pub struct WireRef {
    pub own_gen: OwnGen,
    pub wire_id: WireId,
    pub face_id: FaceId,
    pub orientation: rcad_kernel::topods::Orientation,
}

/// Vertex → Edge endpoint.
#[derive(Debug, Clone)]
pub struct VertexRef {
    pub own_gen: OwnGen,
    pub vertex_id: VertexId,
    pub edge_id: EdgeId,
    pub orientation: rcad_kernel::topods::Orientation,
}

/// Solid → CompSolid / Compound membership.
#[derive(Debug, Clone)]
pub struct SolidRef {
    pub own_gen: OwnGen,
    pub solid_id: SolidId,
    /// Parent id (CompSolidId or CompoundId depending on context).
    pub parent_id: u32,
    pub orientation: rcad_kernel::topods::Orientation,
}

/// Child (compound) → Compound membership.
#[derive(Debug, Clone)]
pub struct ChildRef {
    pub own_gen: OwnGen,
    pub child_id: u32,
    pub parent_id: CompoundId,
    pub orientation: rcad_kernel::topods::Orientation,
}

/// Occurrence → Product placement.
#[derive(Debug, Clone)]
pub struct OccurrenceRef {
    pub own_gen: OwnGen,
    pub parent_product_id: ProductId,
    pub child_occurrence_id: OccurrenceId,
}
