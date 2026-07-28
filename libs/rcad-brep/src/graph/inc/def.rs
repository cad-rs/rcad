//! Entity definition structs.
//!
//! OCCT BRepGraphInc: BRepGraphInc_Definition.hxx
//!
//! Each struct holds the intrinsic data for one topological entity kind.
//! Relation tables (parent→child lists) are stored separately in `super::rel`.

use glam::DVec3;

use crate::graph::inc::id::*;

// ── Topology entities ───────────────────────────────────────────────────────

/// Vertex: a point in 3D with tolerance.
#[derive(Debug, Clone)]
pub struct VertexDef {
    /// 3D position (definition frame, no Location applied).
    pub point: DVec3,
    /// Vertex tolerance (coincidence radius).
    pub tolerance: f64,
}

/// Edge: a 3D curve segment between two vertices.
#[derive(Debug, Clone)]
pub struct EdgeDef {
    /// 3D curve representation (optional).
    pub curve_rep_id: Option<Curve3DRepId>,
    /// 3D polygon representation (optional).
    pub polygon_rep_id: Option<Polygon3DRepId>,
    /// Reference to start vertex.
    pub start_vertex_ref_id: VertexRefId,
    /// Reference to end vertex.
    pub end_vertex_ref_id: VertexRefId,
    /// Edge tolerance.
    pub tolerance: f64,
    /// Whether the edge is degenerate (zero-length in 3D).
    pub degenerated: bool,
}

/// CoEdge: oriented usage of an edge within a wire, with PCurve.
///
/// OCCT's Weiler half-edge pattern: each CoEdge binds an edge to a face
/// with a specific orientation and PCurve.
#[derive(Debug, Clone)]
pub struct CoEdgeDef {
    /// Index of the parent edge.
    pub child_edge_id: EdgeId,
    /// Face that this co-edge belongs to (for PCurve lookup).
    pub face_id: FaceId,
    /// Orientation relative to the parent edge.
    pub orientation: rcad_kernel::topods::Orientation,
    /// PCurve representation (optional).
    pub curve_2d_rep_id: Option<Curve2DRepId>,
    /// 2D polygon representation (optional).
    pub polygon_2d_rep_id: Option<Polygon2DRepId>,
    /// Polygon-on-triangulation representation (optional).
    pub polygon_on_tri_rep_id: Option<PolygonOnTriRepId>,
}

/// Wire: an ordered sequence of co-edges forming a closed or open chain.
#[derive(Debug, Clone)]
pub struct WireDef;

/// Face: a surface bounded by wires (outer + holes).
#[derive(Debug, Clone)]
pub struct FaceDef {
    /// Surface representation.
    pub surface_rep_id: Option<SurfaceRepId>,
    /// Triangulation representation (optional).
    pub triangulation_rep_id: Option<TriangulationRepId>,
    /// Face tolerance.
    pub tolerance: f64,
    /// Whether this face has natural restriction (no explicit wire).
    pub natural_restriction: bool,
}

/// Shell: a set of faces forming a connected manifold boundary.
#[derive(Debug, Clone)]
pub struct ShellDef;

/// Solid: a volume bounded by shells.
#[derive(Debug, Clone)]
pub struct SolidDef;

/// Compound: a non-topological grouping of entities.
#[derive(Debug, Clone)]
pub struct CompoundDef;

/// CompSolid: a set of solids sharing faces.
#[derive(Debug, Clone)]
pub struct CompSolidDef;

// ── Assembly entities ───────────────────────────────────────────────────────

/// Product: an assembly node (root or sub-assembly).
#[derive(Debug, Clone)]
pub struct ProductDef;

/// Occurrence: a placement of a child node (product or topology root) within
/// a parent product.
#[derive(Debug, Clone)]
pub struct OccurrenceDef {
    /// The child node index (ProductId for sub-assemblies,
    /// SolidId/CompoundId/etc. for topology roots).
    pub child_node_id: u32,
    /// Kind of the child node (0=Product, 1=Solid, 2=Compound, etc.)
    pub child_kind: u8,
}
