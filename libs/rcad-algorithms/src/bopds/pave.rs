use std::collections::HashSet;
use std::sync::{Arc, RwLock};

use rcad_kernel::geom::{Curve2d, Curve3};

/// A parametric point along an edge's curve (OCCT: BOPDS_Pave).
/// BOPDS_Pave (hxx:27-78).
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
/// BOPDS_PaveBlock (hxx:30-194, cxx:25-358).
#[derive(Debug, Clone)]
pub struct PaveBlock {
    /// BOPDS_PaveBlock::myOriginalEdge (cxx:32).
    ///   Index of the original edge in DS.edges (NO_EDGE for section edge PaveBlocks).
    pub original_edge: usize,
    /// BOPDS_PaveBlock::myPave1 (hxx:186).
    pub pave1: Pave,
    /// BOPDS_PaveBlock::myPave2 (hxx:187).
    pub pave2: Pave,
    /// BOPDS_PaveBlock::myEdge (cxx:31).
    ///   New edge index assigned during result building.
    pub new_edge: Option<usize>,
    /// BOPDS_PaveBlock::myExtPaves (hxx:188, cxx:29-31).
    ///   Extra paves used to split this block via Update().
    pub ext_paves: Vec<Pave>,
    /// BOPDS_PaveBlock::myMFence (hxx:192, cxx:43).
    ///   Dedup fence for AppendExtPave (maps vertex index → seen flag).
    pub ext_paves_fence: HashSet<usize>,
    /// 3D curve of this edge segment (trimmed to [pave1.param, pave2.param]).
    pub curve: Option<Curve3>,
    /// 2D pcurve on face A.
    pub pcurve_on_a: Option<Curve2d>,
    /// 2D pcurve on face B.
    pub pcurve_on_b: Option<Curve2d>,
    /// shrunk range from IntTools_ShrunkRange (myTS1/myTS2, cxx:33-34).
    pub shrunk_range: Option<[f64; 2]>,
    /// myIsSplittable (hxx:193, cxx:35).
    pub is_splittable: bool,
    /// BOPDS_PaveBlock::myShrunkBox (hxx:191).
    ///   Bounding box of the shrunk range, used by GetPBBox.
    pub my_shrunk_box: Option<(glam::DVec3, glam::DVec3)>,
    /// index of the CommonBlock this PaveBlock belongs to,
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

    // ----- Edge accessors (cxx:54-100) -----

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

    /// BOPDS_PaveBlock::AppendExtPave (cxx:167-173).
    pub fn append_ext_pave(&mut self, pave: Pave) {
        if self.ext_paves_fence.insert(pave.vertex_idx) {
            self.ext_paves.push(pave);
        }
    }

    /// BOPDS_PaveBlock::AppendExtPave1 (cxx:177-180).
    pub fn append_ext_pave1(&mut self, pave: Pave) {
        if self.ext_paves_fence.insert(pave.vertex_idx) {
            self.ext_paves.push(pave);
        }
    }

    /// BOPDS_PaveBlock::RemoveExtPave (cxx:184-202).
    pub fn remove_ext_pave(&mut self, vertex_idx: usize) {
        self.ext_paves.retain(|p| p.vertex_idx != vertex_idx);
    }

    /// BOPDS_PaveBlock::IsToUpdate (cxx:220-223).
    pub fn is_to_update(&self) -> bool {
        !self.ext_paves.is_empty()
    }

    /// BOPDS_PaveBlock::ContainsParameter (cxx:227-245).
    pub fn contains_parameter(&self, the_prm: f64, the_tol: f64, the_index: &mut usize) -> bool {
        for pave in &self.ext_paves {
            if (pave.param - the_prm).abs() < the_tol {
                *the_index = pave.vertex_idx;
                return true;
            }
        }
        false
    }

    /// BOPDS_PaveBlock::Update (cxx:249-312).
    pub fn update(&mut self, the_flag: bool) -> Vec<PaveBlock> {
        let mut a_nb = self.ext_paves.len();
        if the_flag {
            a_nb += 2;
        }
        if std::env::var("RCAD_DEBUG_MB").is_ok() {
            eprintln!("[MB_update] a_nb={} the_flag={} ext={}", a_nb, the_flag, self.ext_paves.len());
        }

        // OCCT L288: aNb <= 1 → return (no split possible).
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
            // OCCT L312: skip identical paves (zero-length sub-block)
            if a_pave2.param == a_pave1.param {
                continue;
            }
            let mut pb = PaveBlock::new(self.original_edge, a_pave1, a_pave2);
            pb.original_edge = self.original_edge;
            // OCCT: sub-PBs inherit curve from parent (needed for PostTreatFF edge creation)
            pb.curve = self.curve.clone();
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



/// shared PaveBlock via `Arc<RwLock<PaveBlock>>`.
#[derive(Debug, Clone)]
pub struct SharedPB(pub Arc<RwLock<PaveBlock>>);

impl SharedPB {
    pub fn new(pb: PaveBlock) -> Self { SharedPB(Arc::new(RwLock::new(pb))) }
}

/// BOPDS_CoupleOfPaveBlocks (hxx:28-108).
/// Stores two PaveBlocks and satellite data for PB splitting during
/// PerformNewVertices.  For EF intersections both PBs are the same.
#[derive(Debug, Clone)]
pub struct CoupleOfPaveBlocks {
    /// Index of the EF/EE interference (myIndexInterf).
    pub interf_idx: usize,
    /// The new vertex index (myIndex). Set after vertex creation/fusion.
    pub vertex_index: usize,
    /// First PaveBlock (always the same as second for EF).
    pub pb1: SharedPB,
    /// Second PaveBlock.
    pub pb2: SharedPB,
    /// Tolerance of the new vertex.
    pub tolerance: f64,
}

#[cfg(test)]
mod pave_block_tests {
    use super::*;
    use crate::bopds::common_block::CommonBlock;
    use glam::DVec3;

    #[test]
    fn pave_block_new() {
        let pv1 = Pave { vertex_idx: 0, param: 0.0 };
        let pv2 = Pave { vertex_idx: 1, param: 1.0 };
        let pb = PaveBlock::new(5, pv1, pv2);

        assert_eq!(pb.original_edge, 5);
        assert_eq!(pb.pave1.vertex_idx, 0);
        assert_eq!(pb.pave1.param, 0.0);
        assert_eq!(pb.pave2.vertex_idx, 1);
        assert_eq!(pb.pave2.param, 1.0);
        assert!(pb.ext_paves.is_empty());
        assert!(!pb.is_splittable);
        assert!(pb.new_edge.is_none());
    }

    #[test]
    fn pave_block_range_and_indices() {
        let pb = PaveBlock::new(3,
            Pave { vertex_idx: 2, param: 0.5 },
            Pave { vertex_idx: 5, param: 2.5 });
        let (t1, t2) = pb.range();
        assert!((t1 - 0.5).abs() < 1e-15);
        assert!((t2 - 2.5).abs() < 1e-15);
        let (v1, v2) = pb.indices();
        assert_eq!(v1, 2);
        assert_eq!(v2, 5);
    }

    #[test]
    fn pave_block_has_same_bounds() {
        let a = PaveBlock::new(1,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        let b = PaveBlock::new(2,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        let c = PaveBlock::new(3,
            Pave { vertex_idx: 1, param: 0.0 },
            Pave { vertex_idx: 0, param: 1.0 });
        let d = PaveBlock::new(4,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 2, param: 1.0 });

        assert!(a.has_same_bounds(&b), "Same bounds (forward)");
        assert!(a.has_same_bounds(&c), "Same bounds (reversed vertex order)");
        assert!(!a.has_same_bounds(&d), "Different vertex should not match");
    }

    #[test]
    fn pave_block_remove_ext_pave() {
        let mut pb = PaveBlock::new(1,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 2, param: 0.3 });
        pb.append_ext_pave(Pave { vertex_idx: 3, param: 0.7 });
        assert_eq!(pb.ext_paves.len(), 2);

        pb.remove_ext_pave(2);
        assert_eq!(pb.ext_paves.len(), 1);
        assert_eq!(pb.ext_paves[0].vertex_idx, 3);

        // Remove non-existent → no-op
        pb.remove_ext_pave(99);
        assert_eq!(pb.ext_paves.len(), 1);
    }

    #[test]
    fn pave_block_is_to_update() {
        let mut pb = PaveBlock::new(1,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        assert!(!pb.is_to_update(), "No ext_paves → not to update");

        pb.append_ext_pave(Pave { vertex_idx: 2, param: 0.5 });
        assert!(pb.is_to_update(), "Has ext_paves → to update");
    }

    #[test]
    fn pave_block_contains_parameter() {
        let mut pb = PaveBlock::new(1,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 2, param: 0.5 });

        let mut idx = 0usize;
        assert!(pb.contains_parameter(0.5, 1e-10, &mut idx));
        assert_eq!(idx, 2);

        assert!(!pb.contains_parameter(0.8, 1e-10, &mut idx));
    }

    #[test]
    fn pave_block_update_splits() {
        // BOPDS_PaveBlock::Update splits the block at ext_pave params
        let mut pb = PaveBlock::new(1,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 3, param: 3.0 });
        pb.append_ext_pave(Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 2, param: 2.0 });

        let result = pb.update(true); // the_flag=true includes endpoint paves
        assert_eq!(result.len(), 3, "Three sub-blocks from two ext_paves + endpoints");
        // Sorted endpoints + ext_paves: [0.0, 1.0, 2.0, 3.0] → blocks [0,1], [1,2], [2,3]
        assert!((result[0].pave1.param - 0.0).abs() < 1e-15);
        assert!((result[0].pave2.param - 1.0).abs() < 1e-15);
        assert!((result[1].pave1.param - 1.0).abs() < 1e-15);
        assert!((result[1].pave2.param - 2.0).abs() < 1e-15);
        assert!((result[2].pave1.param - 2.0).abs() < 1e-15);
        assert!((result[2].pave2.param - 3.0).abs() < 1e-15);
        assert!(pb.ext_paves.is_empty(), "Ext_paves drained after update");
    }

    #[test]
    fn pave_block_update_includes_endpoints() {
        let mut pb = PaveBlock::new(1,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 3, param: 3.0 });
        pb.append_ext_pave(Pave { vertex_idx: 2, param: 2.0 }); // One interior split

        let result = pb.update(true); // the_flag=true (include endpoint paves)
        // Endpoints (0.0, 3.0) + ext (2.0) = sorted [0.0, 2.0, 3.0] → 2 blocks
        assert_eq!(result.len(), 2, "Two sub-blocks: [0,2] and [2,3]");
    }

    #[test]
    fn pave_block_new_curve_block() {
        let pb = PaveBlock::new_curve_block();
        assert_eq!(pb.original_edge, NO_EDGE);
        assert!(pb.is_splittable);
        assert!(pb.new_edge.is_none());
    }

    #[test]
    fn pave_block_has_same_bounds_reversed() {
        let a = PaveBlock::new(1,
            Pave { vertex_idx: 5, param: 0.0 },
            Pave { vertex_idx: 9, param: 1.0 });
        let b = PaveBlock::new(2,
            Pave { vertex_idx: 9, param: 1.0 },
            Pave { vertex_idx: 5, param: 0.0 });

        assert!(a.has_same_bounds(&b), "Reversed endpoints should be considered same bounds");
        assert!(b.has_same_bounds(&a), "Symmetric check");
    }

    #[test]
    fn common_block_new() {
        let cb = CommonBlock::new();
        assert!(cb.pave_blocks().is_empty());
        assert!(cb.faces().is_empty());
        assert_eq!(cb.tolerance(), 0.0);
        assert!(cb.edge().is_none());
    }

    #[test]
    fn common_block_add_pave_block() {
        let mut cb = CommonBlock::new();
        cb.add_pave_block(0, 1);
        cb.add_pave_block(2, 1);
        assert_eq!(cb.pave_blocks().len(), 2);
        assert_eq!(cb.pave_blocks()[0], (0, 1));
        assert_eq!(cb.pave_blocks()[1], (2, 1));
    }

    #[test]
    fn common_block_add_face_dedup() {
        let mut cb = CommonBlock::new();
        cb.add_face(1);
        cb.add_face(2);
        cb.add_face(1); // duplicate
        assert_eq!(cb.faces().len(), 2, "Duplicate faces should be ignored");
    }

    #[test]
    fn common_block_set_pave_blocks() {
        let mut cb = CommonBlock::new();
        cb.set_pave_blocks(vec![(0, 1), (1, 2), (2, 3)]);
        assert_eq!(cb.pave_blocks().len(), 3);
        cb.set_pave_blocks(vec![]);
        assert!(cb.pave_blocks().is_empty(), "set_pave_blocks should replace all");
    }

    #[test]
    fn common_block_set_and_get_edge() {
        let mut cb = CommonBlock::new();
        assert!(cb.edge().is_none());
        cb.set_edge(42);
        assert!(cb.edge().is_some());
        assert_eq!(cb.edge(), Some(42));
    }

    #[test]
    fn common_block_tolerance() {
        let mut cb = CommonBlock::new();
        cb.set_tolerance(1.5);
        assert!((cb.tolerance() - 1.5).abs() < 1e-15);
    }

    #[test]
    fn common_block_set_faces() {
        let mut cb = CommonBlock::new();
        cb.set_faces(vec![10, 20, 30]);
        assert_eq!(cb.faces().len(), 3);
        assert_eq!(cb.faces()[0], 10);
        assert_eq!(cb.faces()[2], 30);
    }

    #[test]
    fn common_block_append_faces() {
        let mut cb = CommonBlock::new();
        cb.add_face(1);
        cb.append_faces(&[2, 3, 1]); // 1 is duplicate
        assert_eq!(cb.faces().len(), 3, "append_faces with duplicate should not add dups");
        assert_eq!(cb.faces()[0], 1);
        assert_eq!(cb.faces()[1], 2);
        assert_eq!(cb.faces()[2], 3);
    }

    #[test]
    fn shared_pb_wraps_pave_block() {
        let pb = PaveBlock::new(7,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        let spb = SharedPB::new(pb);
        let read = spb.0.read().unwrap();
        assert_eq!(read.original_edge, 7);
    }

    #[test]
    fn common_block_pave_block1() {
        let mut cb = CommonBlock::new();
        assert!(cb.pave_block1().is_none());
        cb.add_pave_block(5, 1);
        assert_eq!(cb.pave_block1(), Some(5));
    }

    #[test]
    fn common_block_is_pave_block_on_face() {
        let mut cb = CommonBlock::new();
        cb.add_pave_block(0, 1);
        cb.add_pave_block(1, 2);
        assert!(cb.is_pave_block_on_face(1));
        assert!(cb.is_pave_block_on_face(2));
        assert!(!cb.is_pave_block_on_face(3));
    }

    #[test]
    fn common_block_is_pave_block_on_edge() {
        let mut cb = CommonBlock::new();
        assert!(!cb.is_pave_block_on_edge());
        cb.set_edge(5);
        assert!(cb.is_pave_block_on_edge());
    }

    #[test]
    fn common_block_contains() {
        let mut cb = CommonBlock::new();
        cb.add_pave_block(0, 1);
        cb.add_pave_block(2, 1);
        assert!(cb.contains(0));
        assert!(cb.contains(2));
        assert!(!cb.contains(1));
    }

    #[test]
    fn common_block_sort_by_edge() {
        let mut cb = CommonBlock::new();
        cb.add_pave_block(3, 1); // not on real edge
        cb.add_pave_block(1, 2); // on real edge
        cb.add_pave_block(2, 3); // not on real edge
        cb.set_edge(5);
        // sort: PBs where original_edge == 5 → this is checked by the closure
        cb.sort_by_edge(|pb_idx| pb_idx == 1);
        assert_eq!(cb.pave_block1(), Some(1), "Real-edge PB should be first");
    }

    #[test]
    fn common_block_pave_block_on_edge() {
        let mut cb = CommonBlock::new();
        assert!(cb.pave_block_on_edge().is_none());
        cb.add_pave_block(7, 1);
        assert_eq!(cb.pave_block_on_edge(), Some(7));
    }

    #[test]
    fn common_block_dump_contains_info() {
        let mut cb = CommonBlock::new();
        cb.add_pave_block(0, 1);
        cb.set_edge(3);
        let dump = cb.dump();
        assert!(dump.contains("CommonBlock"));
        assert!(dump.contains("3"));
        assert!(dump.contains("1"));
    }
}
