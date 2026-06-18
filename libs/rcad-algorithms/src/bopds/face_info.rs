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
    /// Indices of intersection curves that lie IN this face (both endpoints on face,
    /// traversal through interior). From F-F intersection.
    pub curves_in: BTreeSet<usize>,
    /// Indices of intersection curves that lie ON the face boundary (coincident with
    /// existing boundary edges). From E-F intersection at boundary.
    pub curves_on: BTreeSet<usize>,
    /// Indices of section curves from F-F intersection (section edges crossing the face).
    pub curves_sc: BTreeSet<usize>,
    /// Vertex indices that lie ON this face.
    pub vertices_on: BTreeSet<usize>,
    /// Vertex indices that lie IN this face (from F-F intersection).
    pub vertices_in: BTreeSet<usize>,
}

impl FaceInfo {
    /// True if this face has any intersection curve data (In/On/Sc).
    /// OCCT equivalent: `bHasFaceInfo = myDS->HasFaceInfo(i)` which checks for
    /// any PaveBlocksIn/On/Sc or alone vertices.
    pub fn has_any_curves(&self) -> bool {
        !self.curves_in.is_empty() || !self.curves_on.is_empty() || !self.curves_sc.is_empty()
    }

    /// All intersection curve indices (In + On + Sc combined).
    /// Used for iteration over all curves on this face.
    pub fn all_curves(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self.curves_in.iter().chain(&self.curves_on).chain(&self.curves_sc).copied().collect();
        v.sort();
        v
    }
}
