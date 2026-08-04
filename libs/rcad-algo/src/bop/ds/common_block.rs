use crate::bop::ds::pave::SharedPB;
/// OCCT BOPDS_CommonBlock: groups geometrically coincident PaveBlocks
/// (edges that lie on the same geometry but belong to different faces).
///
/// This is used by:
/// - `PaveFiller::force_interf_ee()` to mark collinear/coincident edges
/// - `BooleanBuilder::fill_same_domain_faces()` to identify shared edges
///   between coplanar faces
///
/// mirrors `BOPDS_CommonBlock` (BOPDS_CommonBlock.cxx / .hxx).
/// Field names follow OCCT convention (`myPaveBlocks`, `myFaces`, `myEdge`).
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct CommonBlock {
    /// PaveBlock handles paired with their face indices.
    /// Each entry stores `(pave_block, face_index)`, matching OCCT's
    /// `myPaveBlocks` list of `Handle(BOPDS_PaveBlock)`. A pool index alone
    /// cannot identify the PaveBlock when an edge is split into several
    /// blocks, so the specific handle is stored (OCCT BOPDS_CommonBlock.hxx).
    myPaveBlocks: Vec<(SharedPB, usize)>,
    /// Face indices that share this common block.
    /// OCCT: `myFaces` (`TColStd_ListOfInteger`).
    myFaces: Vec<usize>,
    /// The real edge index (if this common block sits on a real edge).
    /// OCCT: `myEdge` ?set when the common block corresponds to an edge
    /// from the original shape (not a section edge).
    myEdge: Option<usize>,
    /// Tolerance of this CommonBlock (max deviation between merged PaveBlocks).
    /// OCCT: `myTolerance` ?computed by `BOPAlgo_Tools::ComputeToleranceOfCB`.
    myTolerance: f64,
}
impl CommonBlock {
    /// Create an empty `CommonBlock`.
    /// default constructor.
    pub fn new() -> Self {
        Self {
            myPaveBlocks: Vec::new(),
            myFaces: Vec::new(),
            myEdge: None,
            myTolerance: 0.0,
        }
    }
    /// Add a `PaveBlock` with its associated face index.
    /// `AddPaveBlock(const Handle(BOPDS_PaveBlock)&, const int)`.
    ///
    /// `pb` ?the specific PaveBlock handle.
    /// `face_idx` ?face index in the DS that this PaveBlock belongs to.
    ///
    /// OCCT BOPDS_CommonBlock::AddPaveBlock inserts at the front of the list so
    /// that `PaveBlock1()` returns the first-added block (the one on the real edge).
    /// rcad: defers ordering to `sort_by_edge()` which is called after all blocks
    /// are added, ensuring the real-edge PaveBlock comes first.
    pub fn add_pave_block(&mut self, pb: SharedPB, face_idx: usize) {
        self.myPaveBlocks.push((pb, face_idx));
    }
    /// Replace all `PaveBlock`/face pairs at once.
    /// `SetPaveBlocks(const BOPDS_ListOfPaveBlock&)`.
    ///
    /// Each entry is `(pave_block, face_index)`.
    pub fn set_pave_blocks(&mut self, blocks: Vec<(SharedPB, usize)>) {
        self.myPaveBlocks = blocks;
    }
    /// Add a face index to this common block.
    /// `AddFace(const int)`.
    /// Duplicate faces are silently ignored.
    pub fn add_face(&mut self, face_idx: usize) {
        if !self.myFaces.contains(&face_idx) {
            self.myFaces.push(face_idx);
        }
    }
    /// Replace all face indices at once.
    /// `SetFaces(const TColStd_ListOfInteger&)`.
    pub fn set_faces(&mut self, faces: Vec<usize>) {
        self.myFaces = faces;
    }
    /// Append additional face indices to the existing set.
    /// `AppendFaces(const TColStd_ListOfInteger&)`.
    /// Duplicates are silently ignored.
    pub fn append_faces(&mut self, faces: &[usize]) {
        for &f in faces {
            self.add_face(f);
        }
    }
    /// Get the `PaveBlock`/face pair list.
    /// `PaveBlocks()`.
    pub fn pave_blocks(&self) -> &[(SharedPB, usize)] {
        &self.myPaveBlocks
    }
    /// Get the face indices.
    /// `Faces()`.
    pub fn faces(&self) -> &[usize] {
        &self.myFaces
    }
    /// Get the first `PaveBlock` handle, if any.
    /// `PaveBlock1()`.
    pub fn pave_block1(&self) -> Option<SharedPB> {
        self.myPaveBlocks.first().map(|(pb, _)| pb.clone())
    }
    /// sort PaveBlocks so that the block on the real edge
    ///   (original_edge == myEdge) is first, matching OCCT's insertion order
    ///   invariant where `AddPaveBlock` inserts the real-edge block first.
    ///   `is_real_edge` is a closure `(&SharedPB) -> bool` that returns true when
    ///   the PaveBlock has `original_edge == self.myEdge`.
    ///   Call this after all PaveBlocks have been added, before using `PaveBlock1()`.
    pub fn sort_by_edge(&mut self, is_real_edge: impl Fn(&SharedPB) -> bool) {
        if self.myEdge.is_none() || self.myPaveBlocks.len() < 2 {
            return;
        }
        // Find the first PaveBlock on the real edge and swap it to front
        if let Some(pos) = self
            .myPaveBlocks
            .iter()
            .position(|(pb, _)| is_real_edge(pb))
        {
            if pos != 0 {
                self.myPaveBlocks.swap(0, pos);
            }
        }
    }
    /// Get the `PaveBlock` on the real edge.
    /// `PaveBlockOnEdge()`.
    ///
    /// Returns the first `PaveBlock` handle as a proxy (the first PaveBlock is
    /// usually the one on the real edge).
    pub fn pave_block_on_edge(&self) -> Option<SharedPB> {
        self.myPaveBlocks.first().map(|(pb, _)| pb.clone())
    }
    /// Check if the given face index has a `PaveBlock` in this common block.
    /// `IsPaveBlockOnFace(const int)`.
    pub fn is_pave_block_on_face(&self, face_idx: usize) -> bool {
        self.myPaveBlocks.iter().any(|(_, fi)| *fi == face_idx)
    }
    /// Check if this common block has a `PaveBlock` on a real edge
    /// (i.e., `myEdge` is set).
    /// `IsPaveBlockOnEdge()`.
    pub fn is_pave_block_on_edge(&self) -> bool {
        self.myEdge.is_some()
    }
    /// Check if the given `PaveBlock` is contained in this common block.
    /// `Contains(const Handle(BOPDS_PaveBlock)&)`.
    pub fn contains(&self, pb: &SharedPB) -> bool {
        let ptr = std::sync::Arc::as_ptr(&pb.0) as u64;
        self.myPaveBlocks
            .iter()
            .any(|(p, _)| std::sync::Arc::as_ptr(&p.0) as u64 == ptr)
    }
    /// Set the edge of all PaveBlocks in the common block.
    /// `SetEdge(const int)` (BOPDS_CommonBlock.cxx L195-207): OCCT propagates
    /// the edge to every PaveBlock in `myPaveBlocks`.
    pub fn set_edge(&mut self, edge_idx: usize) {
        self.myEdge = Some(edge_idx);
        for (pb, _) in &self.myPaveBlocks {
            pb.write().edge = edge_idx;
        }
    }
    /// Get the real edge index, if set.
    /// `Edge()`.
    pub fn edge(&self) -> Option<usize> {
        self.myEdge
    }
    /// Set the tolerance of this CommonBlock.
    /// `SetTolerance()`.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.myTolerance = tol;
    }
    /// Get the tolerance of this CommonBlock.
    /// `Tolerance()`.
    pub fn tolerance(&self) -> f64 {
        self.myTolerance
    }
    /// Debug dump.
    /// `Dump()`.
    pub fn dump(&self) -> String {
        let edge_str = match self.myEdge {
            Some(e) => format!("{}", e),
            None => "N/A".to_string(),
        };
        format!(
            "CommonBlock: {} PaveBlocks on edge={}, faces={:?}",
            self.myPaveBlocks.len(),
            edge_str,
            self.myFaces
        )
    }
}
impl Default for CommonBlock {
    fn default() -> Self {
        Self::new()
    }
}