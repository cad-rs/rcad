use std::collections::BTreeSet;

/// Per-face intersection bookkeeping (OCCT: BOPDS_FaceInfo).
///
/// We use `BTreeSet` (not `HashSet`) so iteration order is **deterministic** for the same
/// indices — `HashSet` iteration is insertion-order dependent and can vary between boolean
/// runs, which breaks stable face splitting and `total_surface_area`.
#[derive(Debug, Clone, Default)]
pub struct FaceInfo {
    /// Indices of PaveBlocks that lie ON this face (from E-F intersection).
    pub pave_blocks_on: BTreeSet<usize>,
    /// Indices of intersection curves that lie IN this face (from F-F intersection).
    pub curves_in: BTreeSet<usize>,
    /// Vertex indices that lie ON this face.
    pub vertices_on: BTreeSet<usize>,
    /// Vertex indices that lie IN this face (from F-F intersection).
    pub vertices_in: BTreeSet<usize>,
}
