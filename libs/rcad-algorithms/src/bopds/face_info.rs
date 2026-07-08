use indexmap::IndexSet;

/// Per-face intersection bookkeeping (OCCT: BOPDS_FaceInfo).
///
/// ✅ OCCT-aligned: uses `IndexSet` which preserves insertion order,
/// matching OCCT's `NCollection_IndexedMap` (not `NCollection_Map`/`HashSet`).
#[derive(Debug, Clone, Default)]
pub struct FaceInfo {
    /// Indices of PaveBlocks that lie ON this face (from E-F intersection).
    pub pave_blocks_on: IndexSet<usize>,
    /// Indices of PaveBlocks that lie IN this face (from E-F intersection).
    pub pave_blocks_in: IndexSet<usize>,
    /// ✅ OCCT-aligned: BOPDS_FaceInfo::PaveBlocksSc (hxx:115-117).
    ///   Section curve PaveBlock indices (sub-segments of intersection curves).
    ///   Populated by post_treat_ff from each curve's pave_blocks via update().
    pub pave_blocks_sc: IndexSet<usize>,
    /// IntersectionCurve indices from FF intersection.
    pub curves_sc: IndexSet<usize>,
    /// Vertex indices that lie ON this face.
    pub vertices_on: IndexSet<usize>,
    /// Vertex indices that lie IN this face (from F-F intersection).
    pub vertices_in: IndexSet<usize>,
    /// ✅ OCCT-aligned: BOPDS_FaceInfo::VerticesSc (hxx:122-123).
    ///   Vertex indices from section curves (FF intersection).
    pub vertices_sc: IndexSet<usize>,
}

impl FaceInfo {
    /// True if this face has any interference data (PaveBlocksIn/On/Sc).
    /// OCCT equivalent: `bHasFaceInfo = myDS->HasFaceInfo(i)` which checks for
    /// any PaveBlocksIn/On/Sc or alone vertices.
    pub fn has_any_interference(&self) -> bool {
        !self.pave_blocks_in.is_empty()
            || !self.pave_blocks_on.is_empty()
            || !self.pave_blocks_sc.is_empty()
            || !self.curves_sc.is_empty()
    }

    /// OCCT PaveBlocksSc: section curve indices from FF intersection only.
    pub fn curves_sc_only(&self) -> Vec<usize> {
        self.curves_sc.iter().copied().collect()
    }
}


