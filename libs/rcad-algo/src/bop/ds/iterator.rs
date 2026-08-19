// OCCT BOPDS_Iterator — BVH-based pair enumeration.
// OCCT ref: BOPDS_Iterator.hxx / BOPDS_Iterator.cxx
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::tools::box_tree::{Aabb, BoxTree, PairSelector};
use rcad_kernel::topods::ShapeType;
use std::collections::HashSet;

/// BOPDS_Iterator — iterates pairs of interfering shapes.
pub struct BOPDS_Iterator {
    fuzzy_tol: f64,
    my_lists: Vec<Vec<(usize, usize)>>,
    my_ext_lists: Vec<Vec<(usize, usize)>>,
    my_length: usize,
    my_run_parallel: bool,
    my_use_ext: bool,
    current_list: Vec<(usize, usize)>,
    current_pos: usize,
}

fn type_to_integer_single(t: ShapeType) -> isize {
    match t { ShapeType::Vertex => 7, ShapeType::Edge => 6, ShapeType::Face => 4, _ => 9, }
}

fn type_to_integer(t1: ShapeType, t2: ShapeType) -> usize {
    let i_t1 = type_to_integer_single(t1); let i_t2 = type_to_integer_single(t2);
    match i_t2 * 10 + i_t1 {
        77 => 0, 76 | 67 => 1, 66 => 2, 74 | 47 => 3, 64 | 46 => 4, 44 => 5, _ => 6,
    }
}

pub(crate) fn boxes_overlap(si1: &crate::bop::ds::ShapeInfo, si2: &crate::bop::ds::ShapeInfo, gap: f64) -> bool {
    if si1.bbox.is_void() || si2.bbox.is_void() { return true; }
    let mut b1 = si1.bbox.clone(); let mut b2 = si2.bbox.clone();
    b1.set_gap(b1.get_gap() + gap); b2.set_gap(b2.get_gap() + gap);
    !b1.is_out_box(&b2)
}

impl BOPDS_Iterator {
    pub fn new(fuzzy_tol: f64) -> Self {
        let n = 10;
        let mut my_lists = Vec::with_capacity(n);
        for _ in 0..n { my_lists.push(Vec::new()); }
        let n_ext = Self::nb_ext_interfs();
        let mut my_ext_lists = Vec::with_capacity(n_ext);
        for _ in 0..n_ext { my_ext_lists.push(Vec::new()); }
        BOPDS_Iterator { fuzzy_tol, my_lists, my_ext_lists, my_length: 0, my_run_parallel: false, my_use_ext: false, current_list: Vec::new(), current_pos: 0 }
    }

    pub fn expected_length(&self) -> usize { self.my_length }
    pub fn block_length(&self) -> usize {
        let a_nb_iis = self.my_length; let a_cf_predict = 0.5;
        if a_nb_iis <= 1 { return 1; } (a_cf_predict * a_nb_iis as f64) as usize
    }
    pub fn set_run_parallel(&mut self, flag: bool) { self.my_run_parallel = flag; }
    pub fn run_parallel(&self) -> bool { self.my_run_parallel }
    pub fn nb_ext_interfs() -> usize { 4 }

    // OCCT BOPDS_Iterator::Prepare (L247-265).
    pub fn prepare(&mut self, ds: &DS, _ctx: Option<&IntToolsContext>, _check_obb: bool, _fuzzy: f64) {
        let a_nb_interf_types = DS::nb_interf_types();
        self.my_length = 0; self.my_use_ext = false;
        for i in 0..a_nb_interf_types { self.my_lists[i].clear(); }
        if ds.nb_shapes() == 0 { return; }
        self.intersect(ds);
    }

    // OCCT BOPDS_Iterator::Intersect (L270-359).
    fn intersect(&mut self, ds: &DS) {
        // L273: const int aNb = myDS->NbSourceShapes();
        let a_nb = ds.nb_source_shapes();
        // L276: BOPTools_BoxTree aBoxTree;
        let mut a_box_tree = BoxTree::new();
        // L277: aBoxTree.SetSize(aNb);
        a_box_tree.set_size(a_nb);
        // L279-289: for (int i = 0; i < aNb; ++i)
        for i in 0..a_nb {
            // L282: const BOPDS_ShapeInfo& aSI = myDS->ShapeInfo(i);
            let a_si = &ds.shapes[i];
            // L283: if (!aSI.HasBRep()) continue;
            if !a_si.shape_type.has_brep() { continue; }
            // L284: aBoxTree.Add(i, Bnd_Tools::Bnd2BVH(aSI.Box()));
            if let Some((xmn, ymn, zmn, xmx, ymx, zmx)) = a_si.bbox.get() {
                a_box_tree.add(i, Aabb { min: glam::DVec3::new(xmn, ymn, zmn), max: glam::DVec3::new(xmx, ymx, zmx), gap: 0.0 });
            }
        }
        // L291: aBoxTree.Build();
        a_box_tree.build();
        // L294: BOPTools_BoxPairSelector aPairSelector;
        let mut a_ps = PairSelector::new();
        // L295: aPairSelector.SetBVHSets(&aBoxTree, &aBoxTree);
        a_ps.set_bvh_sets(&a_box_tree);
        // L296: aPairSelector.SetSame(true);
        a_ps.set_same(true);
        // L297: aPairSelector.Select();
        a_ps.select();
        // L298: aPairSelector.Sort();
        a_ps.sort();
        // L302: const auto& aPairs = aPairSelector.Pairs();
        let a_pairs = a_ps.pairs();
        // L303: const int aNbPairs = (int)aPairs.Size();
        let a_nb_pairs = a_pairs.len();
        // L306: int iPair = 0;
        let mut i_pair = 0;
        // L307: const int aNbR = myDS->NbRanges();
        let a_nb_r = ds.nb_ranges();
        // L308: for (int iR = 0; iR < aNbR; ++iR)
        for i_r in 0..a_nb_r {
            // L309: const BOPDS_IndexRange& aRange = myDS->Range(iR);
            let a_range = ds.range(i_r);
            // L310-354: for (; iPair < aNbPairs; ++iPair)
            while i_pair < a_nb_pairs {
                // L312: const auto& aPair = aPairs[iPair];
                let (id1, id2) = a_pairs[i_pair];
                // L313: if (!aRange.Contains(aPair.ID1)) { break; }
                if !a_range.contains(id1) { break; }
                // L317: if (aRange.Contains(aPair.ID2)) { continue; }
                if a_range.contains(id2) { i_pair += 1; continue; }
                // L326-329: ShapeInfo + ShapeType
                let a_si1 = &ds.shapes[id1]; let a_si2 = &ds.shapes[id2];
                let a_type1 = a_si1.shape_type; let a_type2 = a_si2.shape_type;
                let i_type1 = type_to_integer_single(a_type1); let i_type2 = type_to_integer_single(a_type2);
                // L336-340: sub-shape self-interference check
                if ((i_type1 < i_type2) && a_si1.has_sub_shape(id2))
                    || ((i_type1 > i_type2) && a_si2.has_sub_shape(id1)) { i_pair += 1; continue; }
                // L352: int iX = BOPDS_Tools::TypeToInteger(aType1, aType2);
                let i_x = type_to_integer(a_type1, a_type2);
                // L353: myLists(iX).Append(BOPDS_Pair(min, max)).
                // OCCT BOPDS_Iterator::Value() (BOPDS_Iterator.cxx L226-243) then
                // returns the pair with the HIGHER shape type first (vertex before
                // edge before face). Store in that Value() order so that callers
                // (FillShrunkData, PerformVF/EF) see (type1, type2) ordering,
                // exactly as OCCT's myIterator->Value(nS[0], nS[1]) does.
                let (i_1, i_2) = if i_type1 == i_type2 {
                    (id1.min(id2), id1.max(id2))
                } else if i_type1 > i_type2 {
                    (id1, id2)
                } else {
                    (id2, id1)
                };
                self.my_lists[i_x].push((i_1, i_2));
                i_pair += 1;
            }
        }
    }

    // OCCT BOPDS_Iterator::IntersectExt (L363-463).
    pub fn intersect_ext(&mut self, ds: &DS, the_indices: &HashSet<usize>) {
        if ds.nb_shapes() == 0 { return; }
        let a_nb = ds.nb_source_shapes();
        let mut a_mp_fence: HashSet<(usize, usize)> = HashSet::new();
        for &i in the_indices {
            let si = &ds.shapes[i];
            if si.shape_type != ShapeType::Vertex { continue; }
            let i_rank_i = ds.rank(i);
            let mut n_vsd = i; ds.has_shape_sd(i, &mut n_vsd);
            let si_sd = &ds.shapes[n_vsd];
            let i_ti = type_to_integer_single(si_sd.shape_type);
            for j in 0..a_nb {
                if j == i { continue; }
                let sj = &ds.shapes[j]; let i_rank_j = ds.rank(j);
                if i_rank_i >= 0 && i_rank_i == i_rank_j { continue; }
                let i_tj = type_to_integer_single(sj.shape_type);
                if ((i_ti < i_tj) && si_sd.has_sub_shape(j))
                    || ((i_ti > i_tj) && sj.has_sub_shape(i)) { continue; }
                if boxes_overlap(si_sd, sj, self.fuzzy_tol) {
                    // Same Value() ordering as Intersect: higher type first.
                    let a_pair = if i_ti == i_tj {
                        if i < j { (i, j) } else { (j, i) }
                    } else if i_ti > i_tj {
                        (i, j)
                    } else {
                        (j, i)
                    };
                    if a_mp_fence.insert(a_pair) {
                        let i_x = type_to_integer(si_sd.shape_type, sj.shape_type);
                        if i_x < Self::nb_ext_interfs() { self.my_ext_lists[i_x].push(a_pair); }
                    }
                }
            }
        }
        self.my_use_ext = true;
    }

    // OCCT BOPDS_Iterator::Initialize (L192-208).
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        self.my_length = 0; let i_x = type_to_integer(t1, t2);
        if i_x < self.my_lists.len() {
            let pairs = if self.my_use_ext && i_x < Self::nb_ext_interfs() { &self.my_ext_lists[i_x] } else { &self.my_lists[i_x] };
            let mut sorted = pairs.clone(); sorted.sort();
            self.current_list = sorted; self.my_length = self.current_list.len();
        } else { self.current_list = Vec::new(); }
        self.current_pos = 0;
    }
    pub fn more(&self) -> bool { self.current_pos < self.current_list.len() }
    pub fn next(&mut self) { self.current_pos += 1; }
    pub fn value(&self) -> (usize, usize) { self.current_list[self.current_pos] }
    pub fn pairs(&self, t1: ShapeType, t2: ShapeType) -> &[(usize, usize)] {
        let idx = type_to_integer(t1, t2);
        if self.my_use_ext && idx < self.my_ext_lists.len() { &self.my_ext_lists[idx] }
        else if idx < self.my_lists.len() { &self.my_lists[idx] } else { &[] }
    }
}
