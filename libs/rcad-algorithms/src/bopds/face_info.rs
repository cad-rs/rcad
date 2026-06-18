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
    /// Indices of PaveBlocks that lie IN this face (from E-F intersection).
    pub pave_blocks_in: BTreeSet<usize>,
    /// OCCT PaveBlocksSc: IntersectionCurve indices from FF intersection.
    pub curves_sc: BTreeSet<usize>,
    /// Vertex indices that lie ON this face.
    pub vertices_on: BTreeSet<usize>,
    /// Vertex indices that lie IN this face (from F-F intersection).
    pub vertices_in: BTreeSet<usize>,
}

impl FaceInfo {
    /// True if this face has any interference data (PaveBlocksIn/On/Sc).
    /// OCCT equivalent: `bHasFaceInfo = myDS->HasFaceInfo(i)` which checks for
    /// any PaveBlocksIn/On/Sc or alone vertices.
    pub fn has_any_interference(&self) -> bool {
        !self.pave_blocks_in.is_empty()
            || !self.pave_blocks_on.is_empty()
            || !self.curves_sc.is_empty()
    }

    /// OCCT PaveBlocksSc: section curve indices from FF intersection only.
    pub fn curves_sc_only(&self) -> Vec<usize> {
        self.curves_sc.iter().copied().collect()
    }
}
