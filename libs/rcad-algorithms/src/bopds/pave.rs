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

#[cfg(test)]
mod tests {
    use super::*;

    // ===== Pave: Eq semantics (BOPDS_Pave::IsEqual) =====

    #[test]
    fn pave_eq_matches_both_vertex_and_param() {
        let a = Pave { vertex_idx: 1, param: 0.5 };
        let b = Pave { vertex_idx: 1, param: 0.5 };
        let c = Pave { vertex_idx: 1, param: 0.7 };
        let d = Pave { vertex_idx: 2, param: 0.5 };
        assert_eq!(a, b);
        assert_ne!(a, c);
        assert_ne!(a, d);
    }

    #[test]
    fn pave_partial_cmp_orders_by_param_only() {
        let early = Pave { vertex_idx: 9, param: 0.2 };
        let mid   = Pave { vertex_idx: 3, param: 0.5 };
        let late  = Pave { vertex_idx: 7, param: 0.8 };
        assert!(early < mid);
        assert!(mid < late);
        assert!(early < late);
        // Same param → not less (PartialOrd returns Equal, not Less)
        let same = Pave { vertex_idx: 99, param: 0.5 };
        assert!(!(mid < same));
        assert!(!(same < mid));
    }

    // ===== PaveBlock: constructor invariants =====

    #[test]
    fn new_block_not_splittable_and_no_ext_paves() {
        let p1 = Pave { vertex_idx: 0, param: 0.0 };
        let p2 = Pave { vertex_idx: 1, param: 1.0 };
        let pb = PaveBlock::new(5, p1, p2);
        assert_eq!(pb.original_edge, 5);
        assert_eq!(pb.pave1, p1);
        assert_eq!(pb.pave2, p2);
        assert!(!pb.is_splittable);
        assert!(pb.ext_paves.is_empty());
        assert!(pb.ext_paves_fence.is_empty());
        assert!(pb.new_edge.is_none());
    }

    #[test]
    fn new_curve_block_uses_no_edge_and_splittable() {
        let pb = PaveBlock::new_curve_block();
        assert_eq!(pb.original_edge, NO_EDGE);
        assert_eq!(pb.pave1.vertex_idx, NO_EDGE);
        assert_eq!(pb.pave2.vertex_idx, NO_EDGE);
        assert!(pb.is_splittable);
        assert!(pb.ext_paves.is_empty());
    }

    // ===== Accessors =====

    #[test]
    fn range_and_indices_from_paves() {
        let p1 = Pave { vertex_idx: 2, param: 0.3 };
        let p2 = Pave { vertex_idx: 5, param: 1.7 };
        let pb = PaveBlock::new(0, p1, p2);
        assert_eq!(pb.range(), (0.3, 1.7));
        assert_eq!(pb.indices(), (2, 5));
    }

    #[test]
    fn has_same_bounds_matches_both_orientations() {
        let a = PaveBlock::new(0,
            Pave { vertex_idx: 1, param: 0.0 },
            Pave { vertex_idx: 2, param: 1.0 });
        // exact same orientation
        let b = PaveBlock::new(0,
            Pave { vertex_idx: 1, param: 0.0 },
            Pave { vertex_idx: 2, param: 1.0 });
        // reversed orientation (edge is undirected)
        let c = PaveBlock::new(0,
            Pave { vertex_idx: 2, param: 1.0 },
            Pave { vertex_idx: 1, param: 0.0 });
        assert!(a.has_same_bounds(&b));
        assert!(a.has_same_bounds(&c));
        assert!(b.has_same_bounds(&c));
    }

    #[test]
    fn has_same_bounds_rejects_different_vertices() {
        let a = PaveBlock::new(0,
            Pave { vertex_idx: 1, param: 0.0 },
            Pave { vertex_idx: 2, param: 1.0 });
        let d = PaveBlock::new(0,
            Pave { vertex_idx: 1, param: 0.0 },
            Pave { vertex_idx: 3, param: 1.1 });
        assert!(!a.has_same_bounds(&d));
    }

    #[test]
    fn is_split_edge_cases() {
        let pb = PaveBlock::new(7,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        // None → not split
        assert!(!pb.is_split_edge());
        // Some(original) → not split
        let mut e2 = pb.clone();
        e2.new_edge = Some(7);
        assert!(!e2.is_split_edge());
        // Some(different) → split
        e2.new_edge = Some(42);
        assert!(e2.is_split_edge());
    }

    // ===== ExtPave fence management =====

    #[test]
    fn append_ext_pave_dedup_by_vertex_idx() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 3, param: 0.5 });
        assert_eq!(pb.ext_paves.len(), 1);
        // Same vertex_idx, different param → still deduped (OCCT dedup by vertex_idx)
        pb.append_ext_pave(Pave { vertex_idx: 3, param: 0.51 });
        assert_eq!(pb.ext_paves.len(), 1, "second add with same vertex_idx must be deduped");
        // Different vertex_idx → allowed
        pb.append_ext_pave(Pave { vertex_idx: 4, param: 0.75 });
        assert_eq!(pb.ext_paves.len(), 2);
    }

    #[test]
    fn append_ext_pave1_no_dedup() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave1(Pave { vertex_idx: 3, param: 0.5 });
        pb.append_ext_pave1(Pave { vertex_idx: 3, param: 0.5 });
        assert_eq!(pb.ext_paves.len(), 2, "append_ext_pave1 must skip dedup fence");
    }

    #[test]
    fn remove_ext_pave_removes_from_vec_and_fence() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 3, param: 0.3 });
        pb.append_ext_pave(Pave { vertex_idx: 4, param: 0.6 });
        assert!(pb.ext_paves_fence.contains(&3));
        assert_eq!(pb.ext_paves.len(), 2);

        pb.remove_ext_pave(3);
        assert!(!pb.ext_paves_fence.contains(&3));
        assert_eq!(pb.ext_paves.len(), 1);
        assert_eq!(pb.ext_paves[0].vertex_idx, 4);
    }

    #[test]
    fn remove_ext_pave_unknown_is_noop() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 3, param: 0.5 });
        pb.remove_ext_pave(42); // not in fence
        assert_eq!(pb.ext_paves.len(), 1);
    }

    #[test]
    fn is_to_update_reflects_ext_paves() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        assert!(!pb.is_to_update());
        pb.append_ext_pave(Pave { vertex_idx: 3, param: 0.5 });
        assert!(pb.is_to_update());

        // fence removed but vec still has the entry — OCCT's IsToUpdate checks
        // myExtPaves->Extent(), which is the LIST size, not fence.
        // rcad uses ext_paves.len(), which matches the vec size.
        pb.ext_paves_fence.remove(&3);
        assert!(pb.is_to_update(), "is_to_update must use ext_paves.len(), not fence");
    }

    // ===== ContainsParameter =====

    #[test]
    fn contains_parameter_matches_within_tolerance() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.5 });
        let mut idx = 0;
        assert!(pb.contains_parameter(0.5001, 1e-3, &mut idx));
        assert_eq!(idx, 5);
    }

    #[test]
    fn contains_parameter_rejects_outside_tolerance() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.5 });
        let mut idx = 0;
        assert!(!pb.contains_parameter(0.6, 0.01, &mut idx));
        assert_eq!(idx, 0, "index must remain unchanged on miss");
    }

    #[test]
    fn contains_parameter_only_checks_ext_paves() {
        // does NOT check pave1 or pave2
        let pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        let mut idx = 0;
        assert!(!pb.contains_parameter(0.0, 1e-15, &mut idx), "must not match pave1");
        assert!(!pb.contains_parameter(1.0, 1e-15, &mut idx), "must not match pave2");
    }

    // ===== Update (the critical split method, cxx:249-312) =====

    /// OCCT: a_nb=0, the_flag=false → no split.
    #[test]
    fn update_no_ext_paves_no_flag_returns_empty() {
        let mut pb = PaveBlock::new(3,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        assert!(pb.update(false).is_empty());
    }

    /// OCCT: a_nb=1 (single ext pave) → a_nb<=1 → no split, clears ext_paves.
    #[test]
    fn update_single_ext_pave_no_flag_clears_and_returns_empty() {
        let mut pb = PaveBlock::new(3,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.5 });
        assert!(pb.update(false).is_empty(), "a_nb=1 must produce no sub-blocks");
        assert!(pb.ext_paves.is_empty(), "ext_paves must be cleared");
        assert!(pb.ext_paves_fence.is_empty(), "fence must be cleared");
    }

    /// OCCT: a_nb=2, the_flag=false → two ext paves → one sub-block.
    #[test]
    fn update_two_ext_paves_no_flag_produces_one_sub_block() {
        let mut pb = PaveBlock::new(3,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.3 });
        pb.append_ext_pave(Pave { vertex_idx: 6, param: 0.7 });
        let result = pb.update(false);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].indices(), (5, 6));
        assert_eq!(result[0].range(), (0.3, 0.7));
        assert_eq!(result[0].original_edge, 3, "sub-block must inherit original_edge");
        assert!(pb.ext_paves.is_empty(), "original ext_paves cleared");
    }

    /// OCCT: a_nb=3, the_flag=false → three ext paves → two sub-blocks.
    #[test]
    fn update_three_ext_paves_no_flag_produces_two_sub_blocks() {
        let mut pb = PaveBlock::new(3,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.2 });
        pb.append_ext_pave(Pave { vertex_idx: 6, param: 0.5 });
        pb.append_ext_pave(Pave { vertex_idx: 7, param: 0.8 });
        let result = pb.update(false);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].indices(), (5, 6));
        assert_eq!(result[0].range(), (0.2, 0.5));
        assert_eq!(result[1].indices(), (6, 7));
        assert_eq!(result[1].range(), (0.5, 0.8));
    }

    /// OCCT: the_flag=true includes pave1/pave2 as boundaries.
    /// a_nb = 0 + 2 = 2 → one sub-block = [pave1, pave2] (the full range).
    #[test]
    fn update_flag_true_without_ext_paves_covers_full_range() {
        let mut pb = PaveBlock::new(3,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        let result = pb.update(true);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].indices(), (0, 1));
        assert_eq!(result[0].range(), (0.0, 1.0));
    }

    /// OCCT: the_flag=true + two ext paves → a_nb=4 → three sub-blocks.
    /// total sorted: [pave1(0.0), ext1(0.3), ext2(0.7), pave2(1.0)]
    #[test]
    fn update_flag_true_with_ext_paves_includes_range_bounds_as_endpoints() {
        let mut pb = PaveBlock::new(3,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.3 });
        pb.append_ext_pave(Pave { vertex_idx: 6, param: 0.7 });
        let result = pb.update(true);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].indices(), (0, 5));
        assert_eq!(result[0].range(), (0.0, 0.3));
        assert_eq!(result[1].indices(), (5, 6));
        assert_eq!(result[1].range(), (0.3, 0.7));
        assert_eq!(result[2].indices(), (6, 1));
        assert_eq!(result[2].range(), (0.7, 1.0));
    }

    /// OCCT: update clears ext_paves and fence regardless of split.
    #[test]
    fn update_clears_ext_paves_and_fence() {
        let mut pb = PaveBlock::new(3,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.3 });
        pb.append_ext_pave(Pave { vertex_idx: 6, param: 0.7 });
        let _ = pb.update(false);
        assert!(pb.ext_paves.is_empty(), "ext_paves must be drained");
        assert!(pb.ext_paves_fence.is_empty(), "fence must be cleared");
    }

    /// OCCT: each sub-block from update inherits the original_edge.
    #[test]
    fn update_sub_blocks_inherit_original_edge() {
        let mut pb = PaveBlock::new(42,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.append_ext_pave(Pave { vertex_idx: 5, param: 0.3 });
        pb.append_ext_pave(Pave { vertex_idx: 6, param: 0.7 });
        let result = pb.update(false);
        for sb in &result {
            assert_eq!(sb.original_edge, 42);
        }
    }

    // ===== ShrunkData =====

    #[test]
    fn shrunk_data_default_state() {
        let pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        assert!(!pb.has_shrunk_data());
        let (ts1, ts2, splittable) = pb.shrunk_data();
        assert_eq!(ts1, 0.0);
        assert_eq!(ts2, 0.0);
        assert!(!splittable);
    }

    #[test]
    fn set_shrunk_data_roundtrip() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.set_shrunk_data(0.2, 0.8, true);
        assert!(pb.has_shrunk_data());
        let (ts1, ts2, splittable) = pb.shrunk_data();
        assert!((ts1 - 0.2).abs() < 1e-15);
        assert!((ts2 - 0.8).abs() < 1e-15);
        assert!(splittable);
    }

    #[test]
    fn set_shrunk_data_with_box_roundtrip() {
        let mut pb = PaveBlock::new(0,
            Pave { vertex_idx: 0, param: 0.0 },
            Pave { vertex_idx: 1, param: 1.0 });
        pb.set_shrunk_data_with_box(
            0.1, 0.9,
            glam::DVec3::new(0.0, 0.0, 0.0),
            glam::DVec3::new(1.0, 1.0, 0.0),
            true,
        );
        assert!(pb.has_shrunk_data());
        assert!(pb.my_shrunk_box.is_some());
        let (ts1, ts2, splittable) = pb.shrunk_data();
        assert!((ts1 - 0.1).abs() < 1e-15);
        assert!((ts2 - 0.9).abs() < 1e-15);
        assert!(splittable);
    }
}
