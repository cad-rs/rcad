// OCCT BOPDS_Iterator — BVH-based pair enumeration.
// OCCT ref: BOPDS_Iterator.hxx / BOPDS_Iterator.cxx
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use crate::bop::tools::box_tree::{Aabb, BoxTree};
use rcad_kernel::topods::ShapeType;
use std::collections::HashSet;

/// BOPDS_Iterator — iterates pairs of interfering shapes.
/// OCCT ref: BOPDS_Iterator.hxx
pub struct BOPDS_Iterator {
    fuzzy_tol: f64,
    my_lists: Vec<Vec<(usize, usize)>>,
    my_ext_lists: Vec<Vec<(usize, usize)>>, // OCCT: myExtLists
    my_length: usize, // OCCT: myLength
    my_run_parallel: bool, // OCCT: myRunParallel
    my_use_ext: bool, // OCCT: myUseExt
    current_list: Vec<(usize, usize)>,
    current_pos: usize,
}

// OCCT BOPDS_Tools::TypeToInteger (BOPDS_Tools.lxx L86-123) — single type.
fn type_to_integer_single(t: ShapeType) -> isize {
    match t {
        ShapeType::Vertex => 7,
        ShapeType::Edge => 6,
        ShapeType::Face => 4,
        _ => 9,
    }
}

// OCCT BOPDS_Tools::TypeToInteger (BOPDS_Tools.lxx L31-82) — type pair → list index.
// Uses iT2*10+iT1 formula to handle both orderings symmetrically.
fn type_to_integer(t1: ShapeType, t2: ShapeType) -> usize {
    let i_t1 = type_to_integer_single(t1);
    let i_t2 = type_to_integer_single(t2);
    let i_x = i_t2 * 10 + i_t1;
    match i_x {
        77 => 0,       // VV: 7*10+7
        76 | 67 => 1,  // VE: 7*10+6, EV: 6*10+7
        66 => 2,       // EE: 6*10+6
        74 | 47 => 3,  // VF: 7*10+4, FV: 4*10+7
        64 | 46 => 4,  // EF: 6*10+4, FE: 4*10+6
        44 => 5,       // FF: 4*10+4
        _ => 0,
    }
}

/// Returns true if bounding boxes of two shapes overlap (with gap tolerance).
/// OCCT: Bnd_Box::IsOut (Bnd_Box.cxx L889-966).
fn boxes_overlap(si1: &crate::bop::ds::ShapeInfo, si2: &crate::bop::ds::ShapeInfo, gap: f64) -> bool {
    if si1.bbox.is_void() || si2.bbox.is_void() { return true; }
    let mut b1 = si1.bbox.clone();
    let mut b2 = si2.bbox.clone();
    b1.set_gap(b1.get_gap() + gap);
    b2.set_gap(b2.get_gap() + gap);
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
        BOPDS_Iterator {
            fuzzy_tol,
            my_lists,
            my_ext_lists,
            my_length: 0,
            my_run_parallel: false,
            my_use_ext: false,
            current_list: Vec::new(),
            current_pos: 0,
        }
    }

    // OCCT BOPDS_Iterator.cxx L167-188: ExpectedLength, BlockLength
    pub fn expected_length(&self) -> usize { self.my_length }
    pub fn block_length(&self) -> usize {
        let a_nb_iis = self.my_length;
        let a_cf_predict = 0.5;
        if a_nb_iis <= 1 { return 1; }
        (a_cf_predict * a_nb_iis as f64) as usize
    }

    // OCCT: SetRunParallel / RunParallel
    pub fn set_run_parallel(&mut self, flag: bool) { self.my_run_parallel = flag; }
    pub fn run_parallel(&self) -> bool { self.my_run_parallel }

    // OCCT BOPDS_Iterator.hxx L107: NbExtInterfs
    pub fn nb_ext_interfs() -> usize { 4 }

    // OCCT BOPDS_Iterator::Prepare (BOPDS_Iterator.cxx L247-265).
    pub fn prepare(&mut self, ds: &DS, _ctx: Option<&IntToolsContext>, _check_obb: bool, _fuzzy: f64) {
        let a_nb_interf_types = DS::nb_interf_types();
        self.my_length = 0;
        self.my_use_ext = false;
        for i in 0..a_nb_interf_types {
            self.my_lists[i].clear();
        }
        if ds.nb_shapes() == 0 { return; }
        self.intersect(ds);
    }

    // OCCT BOPDS_Iterator::Intersect (BOPDS_Iterator.cxx L270-359).
    fn intersect(&mut self, ds: &DS) {
        let a_nb = ds.nb_source_shapes();

        // OCCT L279-284: Prepare BVH — add shapes with BRep
        let mut a_box_indices: Vec<usize> = Vec::new();
        let mut a_box_aabbs: Vec<Aabb> = Vec::new();
        for i in 0..a_nb {
            // OCCT L282: const BOPDS_ShapeInfo& aSI = myDS->ShapeInfo(i);
            let a_si = &ds.shapes[i];
            // OCCT L283: if (!aSI.HasBRep()) continue;
            if !a_si.shape_type.has_brep() { continue; }
            // OCCT L284: aBoxTree.Add(i, Bnd_Tools::Bnd2BVH(aSI.Box()));
            if let Some((x_min, y_min, z_min, x_max, y_max, z_max)) = a_si.bbox.get() {
                a_box_indices.push(i);
                a_box_aabbs.push(Aabb {
                    min: glam::DVec3::new(x_min, y_min, z_min),
                    max: glam::DVec3::new(x_max, y_max, z_max),
                    gap: a_si.bbox.get_gap(),
                });
            }
        }

        // OCCT L291: aBoxTree.Build();
        let a_box_tree = BoxTree::build(a_box_indices, a_box_aabbs);

        // OCCT L294-298: BOPTools_BoxPairSelector with same BVH sets
        // OCCT: aPairSelector.SetBVHSets(&aBoxTree, &aBoxTree); SetSame(true);
        // OCCT: aPairSelector.Select(); aPairSelector.Sort();
        let a_pairs = a_box_tree.self_pairs();
        let a_nb_pairs = a_pairs.len();

        // OCCT L306-358: iterate selected pairs by ranges
        let mut i_pair = 0;
        let a_nb_r = ds.nb_ranges();
        for i_r in 0..a_nb_r {
            // OCCT L309: const BOPDS_IndexRange& aRange = myDS->Range(iR);
            let a_range = ds.range(i_r);
            // OCCT L310: for (; iPair < aNbPairs; ++iPair)
            while i_pair < a_nb_pairs {
                let (id1, id2) = a_pairs[i_pair];
                // OCCT L313: if (!aRange.Contains(aPair.ID1)) break;
                if !a_range.contains(id1) { break; }
                // OCCT L317: if (aRange.Contains(aPair.ID2)) continue;
                if a_range.contains(id2) { i_pair += 1; continue; }

                // OCCT L326: const BOPDS_ShapeInfo& aSI1 = myDS->ShapeInfo(aPair.ID1);
                let a_si1 = &ds.shapes[id1];
                let a_si2 = &ds.shapes[id2];
                let a_type1 = a_si1.shape_type;
                let a_type2 = a_si2.shape_type;
                let i_type1 = type_to_integer_single(a_type1);
                let i_type2 = type_to_integer_single(a_type2);

                // OCCT L336-340: sub-shape self-interference check
                if ((i_type1 < i_type2) && a_si1.has_sub_shape(id2))
                    || ((i_type1 > i_type2) && a_si2.has_sub_shape(id1))
                {
                    i_pair += 1;
                    continue;
                }

                // OCCT L346-348: OBB check — skipped (theCheckOBB=false)

                // OCCT L352: int iX = BOPDS_Tools::TypeToInteger(aType1, aType2);
                let i_x = type_to_integer(a_type1, a_type2);
                // OCCT L353: myLists(iX).Append(BOPDS_Pair(min, max));
                self.my_lists[i_x].push((id1.min(id2), id1.max(id2)));
                i_pair += 1;
            }
        }
    }

    // OCCT BOPDS_Iterator::IntersectExt (BOPDS_Iterator.cxx L363-463).
    pub fn intersect_ext(&mut self, ds: &DS, the_indices: &HashSet<usize>) {
        if ds.nb_shapes() == 0 { return; }
        let a_nb = ds.nb_source_shapes();

        // rcad: brute-force equivalent of BVH with increased boxes for map indices.
        // OCCT builds a BVH tree + TSR selectors; we iterate pairs directly.

        // Fence map to avoid duplicating pairs (OCCT L412: NCollection_Map<BOPDS_Pair>)
        let mut a_mp_fence: HashSet<(usize, usize)> = HashSet::new();

        // For each vertex in the map, pair with all other shapes
        for &i in the_indices {
            let si = &ds.shapes[i];
            if si.shape_type != ShapeType::Vertex { continue; } // OCCT only adds VERTEX to map
            let i_rank_i = ds.rank(i);
            // OCCT L386-397: get SD vertex box (with increased tolerance)
            let mut n_vsd = i;
            ds.has_shape_sd(i, &mut n_vsd);
            let si_sd = &ds.shapes[n_vsd];
            let i_ti = type_to_integer_single(si_sd.shape_type);

            for j in 0..a_nb {
                if j == i { continue; }
                let sj = &ds.shapes[j];
                // OCCT L435: if (iRankI == iRankJ) continue;
                let i_rank_j = ds.rank(j);
                if i_rank_i >= 0 && i_rank_i == i_rank_j { continue; }

                // OCCT L440-448: sub-shape check
                let i_tj = type_to_integer_single(sj.shape_type);
                if ((i_ti < i_tj) && si_sd.has_sub_shape(j))
                    || ((i_ti > i_tj) && sj.has_sub_shape(i))
                {
                    continue;
                }

                if boxes_overlap(si_sd, sj, self.fuzzy_tol) {
                    // OCCT L450-458: fence + add to ext lists
                    let a_pair = if i < j { (i, j) } else { (j, i) };
                    if a_mp_fence.insert(a_pair) {
                        let i_x = type_to_integer(si_sd.shape_type, sj.shape_type);
                        if i_x < Self::nb_ext_interfs() {
                            self.my_ext_lists[i_x].push(a_pair);
                        }
                    }
                }
            }
        }

        self.my_use_ext = true; // OCCT L462
    }

    // OCCT BOPDS_Iterator::Initialize (BOPDS_Iterator.cxx L192-208).
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        self.my_length = 0;
        let i_x = type_to_integer(t1, t2);
        if i_x < self.my_lists.len() {
            // OCCT L202: stable_sort for deterministic iteration order
            let pairs = if self.my_use_ext && i_x < Self::nb_ext_interfs() {
                &self.my_ext_lists[i_x]
            } else {
                &self.my_lists[i_x]
            };
            let mut sorted = pairs.clone();
            sorted.sort();
            self.current_list = sorted;
            self.my_length = self.current_list.len();
        } else {
            self.current_list = Vec::new();
        }
        self.current_pos = 0;
    }
    pub fn more(&self) -> bool { self.current_pos < self.current_list.len() }
    pub fn next(&mut self) { self.current_pos += 1; }
    pub fn value(&self) -> (usize, usize) { self.current_list[self.current_pos] }
    pub fn pairs(&self, t1: ShapeType, t2: ShapeType) -> &[(usize, usize)] {
        let idx = type_to_integer(t1, t2);
        if self.my_use_ext && idx < self.my_ext_lists.len() {
            &self.my_ext_lists[idx]
        } else if idx < self.my_lists.len() {
            &self.my_lists[idx]
        } else {
            &[]
        }
    }
}