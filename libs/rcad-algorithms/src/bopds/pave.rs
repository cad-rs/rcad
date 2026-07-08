use std::collections::HashSet;
use std::sync::{Arc, RwLock};

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

    /// OCCT: BOPDS_PaveBlock::IsSplitEdge (cxx:97-100).
    pub fn is_split_edge(&self) -> bool {
        self.new_edge.map_or(false, |e| e != self.original_edge)
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
    pub fn append_ext_pave(&mut self, pave: Pave) {
        if self.ext_paves_fence.insert(pave.vertex_idx) {
            self.ext_paves.push(pave);
        }
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::AppendExtPave1 (cxx:177-180).
    pub fn append_ext_pave1(&mut self, pave: Pave) {
        self.ext_paves.push(pave);
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::RemoveExtPave (cxx:184-202).
    pub fn remove_ext_pave(&mut self, vertex_idx: usize) {
        if self.ext_paves_fence.remove(&vertex_idx) {
            self.ext_paves.retain(|p| p.vertex_idx != vertex_idx);
        }
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::IsToUpdate (cxx:220-223).
    pub fn is_to_update(&self) -> bool {
        !self.ext_paves.is_empty()
    }

    /// ✅ OCCT-aligned: BOPDS_PaveBlock::ContainsParameter (cxx:227-245).
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
    pub fn update(&mut self, the_flag: bool) -> Vec<PaveBlock> {
        let mut a_nb = self.ext_paves.len();
        if the_flag {
            a_nb += 2;
        }
        if std::env::var("RCAD_DEBUG_MB").is_ok() {
            eprintln!("[MB_update] a_nb={} the_flag={} ext={}", a_nb, the_flag, self.ext_paves.len());
        }

        if a_nb <= 1 {
            self.ext_paves.clear();
            self.ext_paves_fence.clear();
            return Vec::new();
        }

        let mut p_paves: Vec<Pave> = Vec::with_capacity(a_nb);
        if the_flag {
            p_paves.push(self.pave1);
            p_paves.push(self.pave2);
        }
        p_paves.extend(self.ext_paves.drain(..));
        self.ext_paves_fence.clear();

        p_paves.sort_by(|a, b| a.param.partial_cmp(&b.param).unwrap_or(std::cmp::Ordering::Equal));

        let mut result = Vec::with_capacity(p_paves.len() - 1);
        let mut a_pave1 = p_paves[0];
        for i in 1..p_paves.len() {
            let a_pave2 = p_paves[i];
            let mut pb = PaveBlock::new(self.original_edge, a_pave1, a_pave2);
            pb.original_edge = self.original_edge;
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



/// OCCT-aligned: shared PaveBlock via `Arc<RwLock<PaveBlock>>`.
#[derive(Debug, Clone)]
pub struct SharedPB(pub Arc<RwLock<PaveBlock>>);

impl SharedPB {
    pub fn new(pb: PaveBlock) -> Self { SharedPB(Arc::new(RwLock::new(pb))) }
}
