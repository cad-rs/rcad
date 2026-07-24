//! BOPDS_Iterator — 1:1 translation of OCCT BOPDS_Iterator
//!
//! OCCT BOPDS_Iterator.hxx / BOPDS_Iterator.cxx
//!
//! Computes intersections between BRep sub-shapes of arguments of an
//! operation (see the class BOPDS_DS) in terms of their bounding boxes,
//! and provides an interface to iterate the pairs of intersected sub-shapes
//! of given type.

use crate::bnd_box::Aabb;
use crate::bopds::ds::DS;
use crate::bvh::BoxTree;
use rcad_kernel::topods::ShapeType;

/// Number of interference types (VV, VE, EE, VF, EF, FF, VZ, EZ, FZ, ZZ)
const NB_INTERF_TYPES: usize = 10;

/// Number of extra interfering types (V/V, V/E, V/F, E/E-unused).
/// Extra lists contain only V/V, V/E, V/F interfering pairs.
/// Although E/E is also initialized (but never filled) for code simplicity.
const NB_EXT_INTERFS: usize = 4;

/// The class BOPDS_Iterator:
/// 1. Computes intersections between BRep sub-shapes of arguments of an
///    operation (see the class BOPDS_DS) in terms of their bounding boxes.
/// 2. Provides an interface to iterate the pairs of intersected sub-shapes
///    of given type.
pub struct BOPDS_Iterator<'a> {
    ds: &'a DS,
    /// Length of the intersection vector of particular intersection type
    my_length: usize,
    /// Pairs with interfering bounding boxes
    my_lists: Vec<Vec<(usize, usize)>>,
    /// Current bucket index for iteration
    my_current_bucket: usize,
    /// Current position within the current bucket
    my_current_pos: usize,
    /// Flag for parallel processing
    my_run_parallel: bool,
    /// Extra pairs of sub-shapes found after intersection of increased sub-shapes
    my_ext_lists: Vec<Vec<(usize, usize)>>,
    /// Information flag for using the extra lists
    my_use_ext: bool,
}

impl<'a> BOPDS_Iterator<'a> {
    /// Empty constructor
    pub fn new(ds: &'a DS) -> Self {
        let mut my_lists = Vec::with_capacity(NB_INTERF_TYPES);
        for _ in 0..NB_INTERF_TYPES {
            my_lists.push(Vec::new());
        }

        let mut my_ext_lists = Vec::with_capacity(NB_EXT_INTERFS);
        for _ in 0..NB_EXT_INTERFS {
            my_ext_lists.push(Vec::new());
        }

        BOPDS_Iterator {
            ds,
            my_length: 0,
            my_lists,
            my_current_bucket: usize::MAX,
            my_current_pos: 0,
            my_run_parallel: false,
            my_ext_lists,
            my_use_ext: false,
        }
    }

    // ===== Modifier / Selector =====

    /// Sets the data structure to process
    pub fn set_ds(&mut self, ds: &'a DS) {
        self.ds = ds;
    }

    /// Returns the data structure
    pub fn ds(&self) -> &DS {
        self.ds
    }

    // ===== Parallel processing =====

    /// Set the flag of parallel processing
    pub fn set_run_parallel(&mut self, flag: bool) {
        self.my_run_parallel = flag;
    }

    /// Returns the flag of parallel processing
    pub fn run_parallel(&self) -> bool {
        self.my_run_parallel
    }

    // ===== ExpectedLength / BlockLength =====

    /// Returns the number of intersections found
    pub fn expected_length(&self) -> usize {
        self.my_length
    }

    /// Returns the block length
    pub fn block_length(&self) -> usize {
        let a_nb_iis = self.expected_length();
        let a_cf_predict = 0.5;
        if a_nb_iis <= 1 {
            return 1;
        }
        (a_cf_predict * a_nb_iis as f64) as usize
    }

    /// Returns a reference to the pre-computed pair list for (t1, t2).
    /// Must be called after `prepare()`.
    pub fn pairs(&self, t1: ShapeType, t2: ShapeType) -> &[(usize, usize)] {
        let i_x = type_to_bucket(t1, t2);
        if i_x >= 0 {
            let i_x = i_x as usize;
            if self.my_use_ext && i_x < NB_EXT_INTERFS {
                &self.my_ext_lists[i_x]
            } else {
                &self.my_lists[i_x]
            }
        } else {
            &[]
        }
    }

    // ===== NbExtInterfs =====

    /// Number of extra interfering types
    pub fn nb_ext_interfs() -> usize {
        NB_EXT_INTERFS
    }

    // ===== Initialize / More / Next / Value =====

    /// Initializes the iterator.
    /// theType1 - the first type of shape
    /// theType2 - the second type of shape
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        // OCCT BOPDS_Iterator::Initialize L192-208
        self.my_length = 0;
        let i_x = type_to_bucket(t1, t2);
        if i_x >= 0 {
            let i_x = i_x as usize;
            let a_pairs = if self.my_use_ext && i_x < NB_EXT_INTERFS {
                &mut self.my_ext_lists[i_x]
            } else {
                &mut self.my_lists[i_x]
            };
            // Sort interfering pairs for constant order of intersection
            a_pairs.sort();
            // Initialize iterator to access the pairs
            self.my_current_bucket = i_x;
            self.my_current_pos = 0;
            self.my_length = a_pairs.len();
        } else {
            self.my_current_bucket = usize::MAX;
            self.my_current_pos = 0;
        }
    }

    /// Returns true if there are still pairs of intersected shapes
    pub fn more(&self) -> bool {
        // OCCT BOPDS_Iterator::More L213-215
        self.my_current_pos < self.current_list_len()
    }

    /// Moves iteration ahead
    pub fn next(&mut self) {
        // OCCT BOPDS_Iterator::Next L219-222
        self.my_current_pos += 1;
    }

    /// Returns indices (DS) of intersected shapes.
    /// theIndex1 - the index of the first shape
    /// theIndex2 - the index of the second shape
    pub fn value(&self) -> (usize, usize) {
        // OCCT BOPDS_Iterator::Value L226-243
        let (n1, n2) = self.current_pair();
        let si1 = self.ds.shape_info_at(n1);
        let si2 = self.ds.shape_info_at(n2);
        let i_t1 = occt_shape_type_index(si1.shape_type);
        let i_t2 = occt_shape_type_index(si2.shape_type);
        if i_t1 < i_t2 {
            (n2, n1)
        } else {
            (n1, n2)
        }
    }

    // ===== Prepare / Intersect =====

    /// Perform the intersection algorithm and prepare the results to be used.
    ///
    /// theCtx        - context (optional, for OBB checks) [not yet used in rcad]
    /// theCheckOBB   - check oriented bounding boxes [not yet used in rcad]
    /// theFuzzyValue - fuzzy tolerance [not yet used in rcad]
    pub fn prepare(&mut self) {
        // OCCT BOPDS_Iterator::Prepare L247-265
        let a_nb_interf_types = NB_INTERF_TYPES;
        for i in 0..a_nb_interf_types {
            self.my_lists[i].clear();
        }

        if self.ds.nb_source_shapes() == 0 {
            return;
        }
        self.intersect();

        // stable_sort each bucket (matching old rcad behavior for pairs() API)
        for list in &mut self.my_lists {
            list.sort();
        }
    }

    /// Intersects the bounding boxes of sub-shapes of the arguments with
    /// the tree and saves the interfering pairs for further geometrical
    /// intersection.
    fn intersect(&mut self) {
        // OCCT BOPDS_Iterator::Intersect L270-359

        let a_nb = self.ds.nb_source_shapes();

        // ---- Prepare BVH ----
        // Collect shapes with BRep and their bounding boxes
        let mut tree_indices: Vec<usize> = Vec::new();
        let mut tree_aabbs: Vec<Aabb> = Vec::new();
        for i in 0..a_nb {
            let si = self.ds.shape_info_at(i);
            if !si.has_brep() {
                continue;
            }
            let (Some(bmin), Some(bmax)) = (si.box_min, si.box_max) else {
                continue;
            };
            tree_indices.push(i);
            tree_aabbs.push(Aabb {
                min: bmin,
                max: bmax,
                gap: si.box_gap,
            });
        }

        if tree_indices.len() < 2 {
            return;
        }

        // Build BVH (OCCT: aBoxTree.Build())
        let a_box_tree = BoxTree::build(tree_indices, tree_aabbs);

        // Select pairs of shapes with interfering bounding boxes
        // (OCCT: BOPTools_BoxPairSelector with SetSame(true), Select(), Sort())
        let mut pairs = a_box_tree.self_pairs();
        pairs.sort();

        // ---- Treat the selected pairs ----
        // Determine ranges for operand A and B
        // In OCCT, NbRanges = number of operands = 2
        // Range(0) = [0, nA_shapes), Range(1) = [nA_shapes, nTotal)
        // Shapes from operand A come first in the flat index.
        let a_end_a = self
            .ds
            .shape_info
            .iter()
            .position(|s| s.rank == 1)
            .unwrap_or(a_nb);

        let ranges: [(usize, usize); 2] = [(0, a_end_a), (a_end_a, a_nb)];

        let a_nb_pairs = pairs.len();
        let mut i_pair: usize = 0;

        for &(r_start, r_end) in &ranges {
            loop {
                if i_pair >= a_nb_pairs { break; }
                let (id1, id2) = pairs[i_pair];
                i_pair += 1;

                // If ID1 is not in the current range, go to the next range
                if id1 < r_start || id1 >= r_end {
                    i_pair -= 1; // undo increment so outer loop continues from this pair
                    break;
                }

                // If ID2 is also in the same range, skip (same operand)
                if id2 >= r_start && id2 < r_end {
                    continue;
                }

                let si1 = self.ds.shape_info_at(id1);
                let si2 = self.ds.shape_info_at(id2);

                let t1 = si1.shape_type;
                let t2 = si2.shape_type;

                let i_type1 = occt_shape_type_index(t1);
                let i_type2 = occt_shape_type_index(t2);

                // Avoid interfering of the shape with its sub-shapes
                if (i_type1 < i_type2 && si1.has_sub_shape(id2))
                    || (i_type1 > i_type2 && si2.has_sub_shape(id1))
                {
                    continue;
                }

                let i_x = type_to_bucket(t1, t2);
                if i_x >= 0 {
                    let i_x = i_x as usize;
                    let (min_id, max_id) = if id1 < id2 { (id1, id2) } else { (id2, id1) };
                    self.my_lists[i_x].push((min_id, max_id));
                }
            }
        }
    }

    /// Updates the tree of Bounding Boxes with increased boxes and
    /// intersects such elements with the tree.
    ///
    /// theIndices - indices of shapes whose boxes have been increased
    pub fn intersect_ext(&mut self, the_indices: &std::collections::HashSet<usize>) {
        // OCCT BOPDS_Iterator::IntersectExt L363-463
        if self.ds.nb_source_shapes() == 0 {
            return;
        }

        let a_nb = self.ds.nb_source_shapes();

        // Build BVH tree with increased boxes
        let mut tree_indices: Vec<usize> = Vec::new();
        let mut tree_aabbs: Vec<Aabb> = Vec::new();

        // First pass: collect all shapes and build the tree
        for i in 0..a_nb {
            let si = self.ds.shape_info_at(i);
            if !si.has_brep() || (si.shape_type == ShapeType::Solid) {
                continue;
            }

            let (bmin, bmax) = match (si.box_min, si.box_max) {
                (Some(min), Some(max)) => (min, max),
                _ => continue,
            };

            if the_indices.contains(&i) {
                // For increased shapes, use the shape's SD (same-domain) partner's box
                // In rcad: find the SD partner via shape_sd
                let n_vsd = self.ds.shape_sd.find_sd_partner(i).unwrap_or(i);
                let si_sd = self.ds.shape_info_at(n_vsd);
                let (sd_bmin, sd_bmax) = match (si_sd.box_min, si_sd.box_max) {
                    (Some(min), Some(max)) => (min, max),
                    _ => (bmin, bmax),
                };
                tree_indices.push(i);
                tree_aabbs.push(Aabb {
                    min: sd_bmin,
                    max: sd_bmax,
                    gap: si_sd.box_gap,
                });
            } else {
                tree_indices.push(i);
                tree_aabbs.push(Aabb {
                    min: bmin,
                    max: bmax,
                    gap: si.box_gap,
                });
            }
        }

        if tree_indices.len() < 2 {
            return;
        }

        let box_tree = BoxTree::build(tree_indices, tree_aabbs);

        // Select overlapping pairs for increased shapes only
        let mut pairs = box_tree.self_pairs();
        pairs.sort();

        // Fence map to avoid duplicating pairs
        let mut a_mp_fence: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();

        for &(id1, id2) in &pairs {
            // Only process pairs involving at least one increased shape
            if !the_indices.contains(&id1) && !the_indices.contains(&id2) {
                continue;
            }

            let si1 = self.ds.shape_info_at(id1);
            let si2 = self.ds.shape_info_at(id2);

            let i_rank_i = si1.rank;
            let i_rank_j = si2.rank;

            // Same operand -> skip
            if i_rank_i == i_rank_j {
                continue;
            }

            let t1 = si1.shape_type;
            let t2 = si2.shape_type;

            let i_type1 = occt_shape_type_index(t1);
            let i_type2 = occt_shape_type_index(t2);

            // Avoid interfering of the shape with its sub-shapes
            if (i_type1 < i_type2 && si1.has_sub_shape(id2))
                || (i_type1 > i_type2 && si2.has_sub_shape(id1))
            {
                continue;
            }

            let (min_id, max_id) = if id1 < id2 { (id1, id2) } else { (id2, id1) };
            let a_pair = (min_id, max_id);

            if a_mp_fence.insert(a_pair) {
                let i_x = type_to_bucket(t1, t2);
                if i_x >= 0 && (i_x as usize) < NB_EXT_INTERFS {
                    self.my_ext_lists[i_x as usize].push(a_pair);
                }
            }
        }

        self.my_use_ext = true;
    }

    // ===== Private helpers =====

    /// Returns the length of the currently selected list
    fn current_list_len(&self) -> usize {
        if self.my_current_bucket == usize::MAX {
            return 0;
        }
        if self.my_use_ext && self.my_current_bucket < NB_EXT_INTERFS {
            self.my_ext_lists[self.my_current_bucket].len()
        } else if self.my_current_bucket < self.my_lists.len() {
            self.my_lists[self.my_current_bucket].len()
        } else {
            0
        }
    }

    /// Returns the current pair from the selected list
    fn current_pair(&self) -> (usize, usize) {
        debug_assert!(
            self.my_current_pos < self.current_list_len(),
            "BOPDS_Iterator::current_pair: position out of range"
        );
        if self.my_use_ext && self.my_current_bucket < NB_EXT_INTERFS {
            self.my_ext_lists[self.my_current_bucket][self.my_current_pos]
        } else {
            self.my_lists[self.my_current_bucket][self.my_current_pos]
        }
    }
}

// ===== OCCT BOPDS_Tools::TypeToInteger equivalents =====

/// OCCT TopAbs_ShapeEnum integer indices:
///   VERTEX=0, EDGE=1, FACE=2, SHELL=3, SOLID=4, WIRE=5, COMPOUND=6, COMPSOLID=7
fn occt_shape_type_index(t: ShapeType) -> isize {
    match t {
        ShapeType::Vertex => 0,
        ShapeType::Edge => 1,
        ShapeType::Face => 2,
        ShapeType::Shell => 3,
        ShapeType::Solid => 4,
        ShapeType::Wire => 5,
        ShapeType::Compound => 6,
        ShapeType::CompSolid => 7,
        _ => -1,
    }
}

/// OCCT BOPDS_Tools::TypeToInteger(theType1, theType2)
/// Returns the bucket index for a pair of shape types:
///   0 = VV, 1 = VE, 2 = EE, 3 = VF, 4 = EF, 5 = FF,
///   6 = VZ, 7 = EZ, 8 = FZ, 9 = ZZ, -1 = invalid
fn type_to_bucket(t1: ShapeType, t2: ShapeType) -> isize {
    let i1 = occt_shape_type_index(t1);
    let i2 = occt_shape_type_index(t2);
    if i1 < 0 || i2 < 0 {
        return -1;
    }
    let i_type = i2 * 10 + i1;
    match i_type {
        0 => 0,   // VV
        1 | 10 => 1, // VE | EV
        11 => 2,  // EE
        2 | 20 => 3, // VF | FV
        12 | 21 => 4, // EF | FE
        22 => 5,  // FF
        3 | 30 => 6, // VZ | ZV
        13 | 31 => 7, // EZ | ZE
        23 | 32 => 8, // FZ | ZF
        33 => 9,  // ZZ
        _ => -1,
    }
}
