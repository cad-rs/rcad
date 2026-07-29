// OCCT BOPDS_Iterator — BVH-based pair enumeration.
// OCCT ref: BOPDS_Iterator.hxx / BOPDS_Iterator.cxx
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
use rcad_kernel::topods::ShapeType;
use std::collections::HashSet;

/// BOPDS_Iterator — iterates pairs of interfering shapes.
/// OCCT ref: BOPDS_Iterator.hxx
pub struct BOPDS_Iterator<'a> {
    ds: &'a DS,
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
fn boxes_overlap(si1: &crate::bop::ds::ShapeInfo, si2: &crate::bop::ds::ShapeInfo, gap: f64) -> bool {
    use rcad_kernel::math::bnd::BndBox;
    let b1 = match (si1.box_min, si1.box_max) {
        (Some(min), Some(max)) => {
            let mut b = BndBox::from_corners(min.x, min.y, min.z, max.x, max.y, max.z);
            b.set_gap(si1.box_gap + gap);
            b
        }
        _ => return true,
    };
    let b2 = match (si2.box_min, si2.box_max) {
        (Some(min), Some(max)) => {
            let mut b = BndBox::from_corners(min.x, min.y, min.z, max.x, max.y, max.z);
            b.set_gap(si2.box_gap + gap);
            b
        }
        _ => return true,
    };
    !b1.is_out_box(&b2)
}

impl<'a> BOPDS_Iterator<'a> {
    pub fn new(ds: &'a DS, fuzzy_tol: f64) -> Self {
        let n = 10;
        let mut my_lists = Vec::with_capacity(n);
        for _ in 0..n { my_lists.push(Vec::new()); }
        let n_ext = Self::nb_ext_interfs();
        let mut my_ext_lists = Vec::with_capacity(n_ext);
        for _ in 0..n_ext { my_ext_lists.push(Vec::new()); }
        BOPDS_Iterator {
            ds, fuzzy_tol,
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
    pub fn prepare(&mut self, _ctx: Option<&IntToolsContext>, _check_obb: bool, _fuzzy: f64) {
        let a_nb_interf_types = DS::nb_interf_types();
        self.my_length = 0;
        self.my_use_ext = false;
        for i in 0..a_nb_interf_types {
            self.my_lists[i].clear();
        }
        if self.ds.nb_shapes() == 0 { return; }
        self.intersect();
    }

    // OCCT BOPDS_Iterator::Intersect (BOPDS_Iterator.cxx L270-359).
    fn intersect(&mut self) {
        let a_nb = self.ds.nb_source_shapes();
        // OCCT L279-284: Add shapes with BRep to BVH
        // OCCT L291: Build BVH
        // OCCT L294-298: Select pairs via BVH + Sort
        // rcad: brute-force AABB + range-filter (same semantics)

        let has_brep = |st: ShapeType| -> bool {
            matches!(st, ShapeType::Vertex | ShapeType::Edge | ShapeType::Face)
        };

        // OCCT L306-358: iterate all (i,j) pairs, skip same-range.
        // OCCT's BVH produces all pairs; then each range consumes those where
        // ID1 is in this range and ID2 is in a different range.
        // rcad: iterate i<j globally, skip if same range.
        for i in 0..a_nb {
            let si1 = &self.ds.shapes[i];
            if !has_brep(si1.shape_type) { continue; }
            let a_type1 = si1.shape_type;
            let i_type1 = type_to_integer_single(a_type1);
            let r1 = self.ds.rank(i);

            for j in (i + 1)..a_nb {
                let si2 = &self.ds.shapes[j];
                if !has_brep(si2.shape_type) { continue; }
                let a_type2 = si2.shape_type;
                let i_type2 = type_to_integer_single(a_type2);

                // OCCT L320-324: skip if both are from the same argument range
                if r1 >= 0 && r1 == self.ds.rank(j) { continue; }

                // OCCT L336-340: avoid interfering shape with its sub-shapes
                if ((i_type1 < i_type2) && si1.has_sub_shape(j))
                    || ((i_type1 > i_type2) && si2.has_sub_shape(i))
                {
                    continue;
                }

                if boxes_overlap(si1, si2, self.fuzzy_tol) {
                    let i_x = type_to_integer(a_type1, a_type2);
                    self.my_lists[i_x].push((i.min(j), i.max(j)));
                }
            }
        }
    }

    // OCCT BOPDS_Iterator::IntersectExt (BOPDS_Iterator.cxx L363-463).
    pub fn intersect_ext(&mut self, the_indices: &HashSet<usize>) {
        if self.ds.nb_shapes() == 0 { return; }
        let a_nb = self.ds.nb_source_shapes();

        // rcad: brute-force equivalent of BVH with increased boxes for map indices.
        // OCCT builds a BVH tree + TSR selectors; we iterate pairs directly.

        // Fence map to avoid duplicating pairs (OCCT L412: NCollection_Map<BOPDS_Pair>)
        let mut a_mp_fence: HashSet<(usize, usize)> = HashSet::new();

        // For each vertex in the map, pair with all other shapes
        for &i in the_indices {
            let si = &self.ds.shapes[i];
            if si.shape_type != ShapeType::Vertex { continue; } // OCCT only adds VERTEX to map
            let i_rank_i = self.ds.rank(i);
            // OCCT L386-397: get SD vertex box (with increased tolerance)
            let mut n_vsd = i;
            self.ds.has_shape_sd(i, &mut n_vsd);
            let si_sd = &self.ds.shapes[n_vsd];
            let i_ti = type_to_integer_single(si_sd.shape_type);

            for j in 0..a_nb {
                if j == i { continue; }
                let sj = &self.ds.shapes[j];
                // OCCT L435: if (iRankI == iRankJ) continue;
                let i_rank_j = self.ds.rank(j);
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