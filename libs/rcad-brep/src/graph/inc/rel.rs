//! Relation tables: ordered parent→child lists stored independently of entities.
//!
//! OCCT BRepGraphInc: BRepGraphInc_Relations.hxx
//!
//! Each struct stores the ordered list of child references for one entity kind.
//! Incoming (reverse) relations are derived via the ref endpoint IDs.

use crate::graph::inc::id::*;

/// Solid → Shells.
#[derive(Debug, Clone, Default)]
pub struct SolidRelations {
    pub shell_ref_ids: Vec<ShellRefId>,
}

/// Shell → Faces.
#[derive(Debug, Clone, Default)]
pub struct ShellRelations {
    pub face_ref_ids: Vec<FaceRefId>,
}

/// Face → Wires.
#[derive(Debug, Clone, Default)]
pub struct FaceRelations {
    pub wire_ref_ids: Vec<WireRefId>,
}

/// Wire → CoEdges.
#[derive(Debug, Clone, Default)]
pub struct WireRelations {
    pub co_edge_ids: Vec<CoEdgeId>,
}

/// Edge → CoEdges (for reverse traversal).
#[derive(Debug, Clone, Default)]
pub struct EdgeRelations {
    pub co_edge_ids: Vec<CoEdgeId>,
}

/// Vertex → Edges (incident edges).
#[derive(Debug, Clone, Default)]
pub struct VertexRelations {
    pub edge_ids: Vec<EdgeId>,
}

/// Compound → Children (solids, shells, faces, or sub-compounds).
#[derive(Debug, Clone, Default)]
pub struct CompoundRelations {
    pub child_ref_ids: Vec<ChildRefId>,
}

/// CompSolid → Solids.
#[derive(Debug, Clone, Default)]
pub struct CompSolidRelations {
    pub solid_ref_ids: Vec<SolidRefId>,
}

/// Product → Occurrences (assembly children).
#[derive(Debug, Clone, Default)]
pub struct ProductRelations {
    pub occurrence_ref_ids: Vec<OccurrenceRefId>,
}

/// Occurrence → Parent product references (reverse).
#[derive(Debug, Clone, Default)]
pub struct OccurrenceRelations {
    pub parent_occurrence_ref_ids: Vec<OccurrenceRefId>,
}

/// All relation tables bundled for easy access.
#[derive(Debug, Clone)]
pub struct RelationTables {
    pub solid: Vec<SolidRelations>,
    pub shell: Vec<ShellRelations>,
    pub face: Vec<FaceRelations>,
    pub wire: Vec<WireRelations>,
    pub edge: Vec<EdgeRelations>,
    pub vertex: Vec<VertexRelations>,
    pub compound: Vec<CompoundRelations>,
    pub comp_solid: Vec<CompSolidRelations>,
    pub product: Vec<ProductRelations>,
    pub occurrence: Vec<OccurrenceRelations>,
}

impl RelationTables {
    pub fn new() -> Self {
        RelationTables {
            solid: Vec::new(), shell: Vec::new(), face: Vec::new(),
            wire: Vec::new(), edge: Vec::new(), vertex: Vec::new(),
            compound: Vec::new(), comp_solid: Vec::new(),
            product: Vec::new(), occurrence: Vec::new(),
        }
    }

    /// Ensure relation slots exist for N solids.
    pub fn prepare_solid(&mut self, n: usize) { self.solid.resize_with(n.max(self.solid.len()), Default::default); }
    pub fn prepare_shell(&mut self, n: usize) { self.shell.resize_with(n.max(self.shell.len()), Default::default); }
    pub fn prepare_face(&mut self, n: usize) { self.face.resize_with(n.max(self.face.len()), Default::default); }
    pub fn prepare_wire(&mut self, n: usize) { self.wire.resize_with(n.max(self.wire.len()), Default::default); }
    pub fn prepare_edge(&mut self, n: usize) { self.edge.resize_with(n.max(self.edge.len()), Default::default); }
    pub fn prepare_vertex(&mut self, n: usize) { self.vertex.resize_with(n.max(self.vertex.len()), Default::default); }
    pub fn prepare_compound(&mut self, n: usize) { self.compound.resize_with(n.max(self.compound.len()), Default::default); }
    pub fn prepare_comp_solid(&mut self, n: usize) { self.comp_solid.resize_with(n.max(self.comp_solid.len()), Default::default); }
}
