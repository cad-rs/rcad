use rcad_kernel::geom::{Curve2d, Curve3};

/// A parametric point along an edge's curve (OCCT: BOPDS_Pave).
#[derive(Debug, Clone, Copy)]
pub struct Pave {
    /// Index of the vertex at this parametric point (in DS.vertices).
    pub vertex_idx: usize,
    /// Parametric value on the edge's curve.
    pub param: f64,
}

/// A segment of an edge between two paves (OCCT: BOPDS_PaveBlock).
/// When an edge is split by intersections, it becomes multiple PaveBlocks.
#[derive(Debug, Clone)]
pub struct PaveBlock {
    /// Original edge index in DS.edges.
    pub original_edge: usize,
    pub pave1: Pave,
    pub pave2: Pave,
    /// New edge index assigned during result building.
    pub new_edge: Option<usize>,
    /// 3D curve of this edge segment (trimmed to [pave1.param, pave2.param]).
    pub curve: Option<Curve3>,
    /// 2D pcurve on face A.
    pub pcurve_on_a: Option<Curve2d>,
    /// 2D pcurve on face B.
    pub pcurve_on_b: Option<Curve2d>,
    /// OCCT-aligned: shrunk range from IntTools_ShrunkRange.
    pub shrunk_range: Option<[f64; 2]>,
    /// OCCT-aligned: whether this PaveBlock can be split.
    pub is_splittable: bool,
}

impl PaveBlock {
    pub fn new(original_edge: usize, pave1: Pave, pave2: Pave) -> Self {
        Self {
            original_edge,
            pave1,
            pave2,
            new_edge: None,
            curve: None,
            pcurve_on_a: None,
            pcurve_on_b: None,
            shrunk_range: None,
            is_splittable: false,
        }
    }
}
