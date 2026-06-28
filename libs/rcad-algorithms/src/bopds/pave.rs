use std::collections::HashSet;

use rcad_kernel::geom::{Curve2d, Curve3};

/// A parametric point along an edge's curve (OCCT: BOPDS_Pave).
/// ✅ OCCT-aligned: BOPDS_Pave (hxx:27-78).
#[derive(Debug, Clone, Copy)]
pub struct Pave {
    /// Index of the vertex at this parametric point (in DS.vertices).
    pub vertex_idx: usize,
    /// Parametric value on the edge's curve.
    pub param: f64,
}

impl Pave {
    /// OCCT: BOPDS_Pave::Index()
    pub fn index(&self) -> usize {
        self.vertex_idx
    }

    /// OCCT: BOPDS_Pave::Parameter()
    pub fn parameter(&self) -> f64 {
        self.param
    }
}

impl PartialEq for Pave {
    fn eq(&self, other: &Self) -> bool {
        self.vertex_idx == other.vertex_idx && self.param == other.param
    }
}

impl PartialOrd for Pave {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.param.partial_cmp(&other.param)
    }
}

/// Sentinel value meaning "no original edge" (section edge PaveBlock).
pub const NO_EDGE: usize = usize::MAX;

/// A segment of an edge between two paves (OCCT: BOPDS_PaveBlock).
/// When an edge is split by intersections, it becomes multiple PaveBlocks.
/// ✅ OCCT-aligned: BOPDS_PaveBlock (hxx:30-194, cxx:25-358).
#[derive(Debug, Clone)]
pub struct PaveBlock {
    /// ✅ OCCT-aligned: BOPDS_PaveBlock::myOriginalEdge (cxx:32).
    ///   Index of the original edge in DS.edges (NO_EDGE for section edge PaveBlocks).
    pub original_edge: usize,
    /// ✅ OCCT-aligned: BOPDS_PaveBlock::myPave1 (hxx:186).
    pub pave1: Pave,
    /// ✅ OCCT-aligned: BOPDS_PaveBlock::myPave2 (hxx:187).
    pub pave2: Pave,
    /// ✅ OCCT-aligned: BOPDS_PaveBlock::myEdge (cxx:31).
    ///   New edge index assigned during result building.
    pub new_edge: Option<usize>,
    /// ✅ OCCT-aligned: BOPDS_PaveBlock::myExtPaves (hxx:188, cxx:29-31).
    ///   Extra paves used to split this block via Update().
    pub ext_paves: Vec<Pave>,
    /// ✅ OCCT-aligned: BOPDS_PaveBlock::myMFence (hxx:192, cxx:43).
    ///   Dedup fence for AppendExtPave (maps vertex index → seen flag).
    pub ext_paves_fence: HashSet<usize>,
    /// 3D curve of this edge segment (trimmed to [pave1.param, pave2.param]).
    pub curve: Option<Curve3>,
    /// 2D pcurve on face A.
    pub pcurve_on_a: Option<Curve2d>,
    /// 2D pcurve on face B.
    pub pcurve_on_b: Option<Curve2d>,
    /// ✅ OCCT-aligned: shrunk range from IntTools_ShrunkRange (myTS1/myTS2, cxx:33-34).
    pub shrunk_range: Option<[f64; 2]>,
    /// ✅ OCCT-aligned: myIsSplittable (hxx:193, cxx:35).
    pub is_splittable: bool,
    /// ✅ OCCT-aligned: BOPDS_PaveBlock::myShrunkBox (hxx:191).
    ///   Bounding box of the shrunk range, used by GetPBBox.
    pub my_shrunk_box: Option<(glam::DVec3, glam::DVec3)>,
    /// ✅ OCCT-aligned: index of the CommonBlock this PaveBlock belongs to,
    ///   or None if not on any CommonBlock (BOPDS_PaveBlock::myCommonBlock).
    pub common_block_idx: Option<usize>,
}

impl PaveBlock {
    /// OCCT: BOPDS_PaveBlock default constructor (cxx:27-36).
    pub fn new(original_edge: usize, pave1: Pave, pave2: Pave) -> Self {
        Self {
            original_edge,
            pave1,
            pave2,
            new_edge: None,
            ext_paves: Vec::new(),
            ext_paves_fence: HashSet::new(),
            curve: None,
            pcurve_on_a: None,
            pcurve_on_b: None,
            shrunk_range: None,
            my_shrunk_box: None,
            is_splittable: false,
            common_block_idx: None,
        }
    }

    /// OCCT: BOPDS_PaveBlock default constructor (cxx:27-36) for curve PaveBlocks.
    ///   myEdge=-1, myOriginalEdge=-1, myTS1=myTS2=-99., myIsSplittable=false,
    ///   then InitPaveBlock1 overrides myIsSplittable=true.
    pub fn new_curve_block() -> Self {
        Self {
            original_edge: NO_EDGE,
            pave1: Pave { vertex_idx: NO_EDGE, param: 0.0 },
            pave2: Pave { vertex_idx: NO_EDGE, param: 0.0 },
            new_edge: None,
            ext_paves: Vec::new(),
            ext_paves_fence: HashSet::new(),
            curve: None,
            pcurve_on_a: None,
            pcurve_on_b: None,
            shrunk_range: None,
            my_shrunk_box: None,
            is_splittable: true,
            common_block_idx: None,
        }
    }

    // ----- OCCT-aligned: Edge accessors (cxx:54-100) -----

    /// OCCT: BOPDS_PaveBlock::SetEdge (cxx:54-57).
    pub fn set_edge(&mut self, edge_idx: usize) {
        self.new_edge = Some(edge_idx);
    }

    /// OCCT: BOPDS_PaveBlock::Edge (cxx:61-63).
    pub fn edge(&self) -> Option<usize> {
        self.new_edge
    }

    /// OCCT: BOPDS_PaveBlock::HasEdge (cxx:68-71).
    pub fn has_edge(&self) -> bool {
        self.new_edge.is_some()
    }

    /// OCCT: BOPDS_PaveBlock::HasEdge(int&) (cxx:75-79).
    pub fn has_edge_value(&self) -> (bool, Option<usize>) {
        (self.new_edge.is_some(), self.new_edge)
    }

    /// OCCT: BOPDS_PaveBlock::SetOriginalEdge (cxx:83-86).
    pub fn set_original_edge(&mut self, edge_idx: usize) {
        self.original_edge = edge_idx;
    }

    /// OCCT: BOPDS_PaveBlock::OriginalEdge (cxx:90-93).
    pub fn original_edge(&self) -> usize {
        self.original_edge
    }

    /// OCCT: BOPDS_PaveBlock::IsSplitEdge (cxx:97-100).
    pub fn is_split_edge(&self) -> bool {
        self.new_edge.map_or(false, |e| e != self.original_edge)
    }

    // ----- Pave accessors (cxx:104-145) -----

    /// OCCT: BOPDS_PaveBlock::SetPave1 (cxx:104-107).
    pub fn set_pave1(&mut self, pave: Pave) {
        self.pave1 = pave;
    }

    /// OCCT: BOPDS_PaveBlock::SetPave2 (cxx:118-121).
    pub fn set_pave2(&mut self, pave: Pave) {
        self.pave2 = pave;
    }

    /// OCCT: BOPDS_PaveBlock::Range (cxx:132-136).
    pub fn range(&self) -> (f64, f64) {
        (self.pave1.param, self.pave2.param)
    }

    /// OCCT: BOPDS_PaveBlock::Indices (cxx:140-144).
    pub fn indices(&self) -> (usize, usize) {
        (self.pave1.vertex_idx, self.pave2.vertex_idx)
    }

    /// OCCT: BOPDS_PaveBlock::HasSameBounds (cxx:148-160).
    pub fn has_same_bounds(&self, other: &Self) -> bool {
        let (n11, n12) = self.indices();
        let (n21, n22) = other.indices();
        (n11 == n21 && n12 == n22) || (n11 == n22 && n12 == n21)
    }

    // ----- ExtPave methods (cxx:167-312) -----

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::AppendExtPave (cxx:167-173).
    ///   Adds an extra pave; deduplicates via myMFence (vertex index fence).
    pub fn append_ext_pave(&mut self, pave: Pave) {
        if self.ext_paves_fence.insert(pave.vertex_idx) {
            self.ext_paves.push(pave);
        }
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::AppendExtPave1 (cxx:177-180).
    ///   Adds an extra pave without dedup fence check.
    pub fn append_ext_pave1(&mut self, pave: Pave) {
        self.ext_paves.push(pave);
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::RemoveExtPave (cxx:184-202).
    ///   Removes the extra pave with the given vertex index.
    pub fn remove_ext_pave(&mut self, vertex_idx: usize) {
        if self.ext_paves_fence.remove(&vertex_idx) {
            self.ext_paves.retain(|p| p.vertex_idx != vertex_idx);
        }
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::ExtPaves (cxx:206-209).
    pub fn ext_paves(&self) -> &[Pave] {
        &self.ext_paves
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::ChangeExtPaves (cxx:213-216).
    pub fn change_ext_paves(&mut self) -> &mut Vec<Pave> {
        &mut self.ext_paves
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::IsToUpdate (cxx:220-223).
    pub fn is_to_update(&self) -> bool {
        !self.ext_paves.is_empty()
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::ContainsParameter (cxx:227-245).
    ///   Returns true if an ext pave exists with |param - thePrm| < theTol.
    ///   Sets theIndex to the matching pave's vertex index.
    pub fn contains_parameter(&self, the_prm: f64, the_tol: f64, the_index: &mut usize) -> bool {
        for pave in &self.ext_paves {
            if (pave.param - the_prm).abs() < the_tol {
                *the_index = pave.vertex_idx;
                return true;
            }
        }
        false
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::Update (cxx:249-312).
    ///   Produces new PaveBlocks from the current block's ext_paves.
    ///   When theFlag=true, includes myPave1/myPave2 as boundaries.
    ///   Clears ext_paves and fence after splitting.
    ///   Returns the list of new PaveBlocks (OCCT: output parameter theLPB).
    pub fn update(&mut self, the_flag: bool) -> Vec<PaveBlock> {
        let mut a_nb = self.ext_paves.len();
        if the_flag {
            a_nb += 2; // include myPave1, myPave2
        }

        if a_nb <= 1 {
            self.ext_paves.clear();
            self.ext_paves_fence.clear();
            return Vec::new();
        }

        // Collect all paves into a sorted array
        let mut p_paves: Vec<Pave> = Vec::with_capacity(a_nb);
        if the_flag {
            p_paves.push(self.pave1);
            p_paves.push(self.pave2);
        }
        p_paves.extend(self.ext_paves.drain(..));
        self.ext_paves_fence.clear();

        // Sort by parameter (OCCT: std::sort with BOPDS_Pave::operator< by parameter)
        p_paves.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));

        // Create new PaveBlocks from adjacent pairs
        let mut result = Vec::with_capacity(p_paves.len() - 1);
        let mut a_pave1 = p_paves[0];
        for i in 1..p_paves.len() {
            let a_pave2 = p_paves[i];
            let mut pb = PaveBlock::new(self.original_edge, a_pave1, a_pave2);
            pb.set_original_edge(self.original_edge);
            result.push(pb);
            a_pave1 = a_pave2;
        }

        result
    }

    // ----- ShrunkData (cxx:317-346) -----

    /// OCCT: BOPDS_PaveBlock::HasShrunkData (cxx:317-320).
    pub fn has_shrunk_data(&self) -> bool {
        self.shrunk_range.is_some()
    }

    /// OCCT: BOPDS_PaveBlock::SetShrunkData (cxx:324-333).
    pub fn set_shrunk_data(&mut self, ts1: f64, ts2: f64, the_is_splittable: bool) {
        self.shrunk_range = Some([ts1, ts2]);
        self.is_splittable = the_is_splittable;
    }

    /// OCCT: BOPDS_PaveBlock::SetShrunkData with bounding box (full sig).
    pub fn set_shrunk_data_with_box(&mut self, ts1: f64, ts2: f64,
        box_min: glam::DVec3, box_max: glam::DVec3, the_is_splittable: bool) {
        self.shrunk_range = Some([ts1, ts2]);
        self.my_shrunk_box = Some((box_min, box_max));
        self.is_splittable = the_is_splittable;
    }

    /// OCCT: BOPDS_PaveBlock::ShrunkData (cxx:337-346).
    pub fn shrunk_data(&self) -> (f64, f64, bool) {
        match self.shrunk_range {
            Some([ts1, ts2]) => (ts1, ts2, self.is_splittable),
            None => (0.0, 0.0, self.is_splittable),
        }
    }
}
