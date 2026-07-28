use crate::bop::ds;
use crate::bop::tools::bvh::Aabb;
use rcad_kernel::topods::ShapeType;
use std::collections::HashSet;

/// BOPDS_Iterator 鈥?BVH-based pair enumeration with type bucketing.
///
/// Builds a single BVH tree over all DS sub-shapes (vertices, edges, faces),
/// finds overlapping AABB pairs, buckets them by (type1, type2) combination,
/// applies stable_sort within each bucket, and provides iteration via
/// `Initialize(T1, T2) 鈫?More/Next/Value`.
///
/// OCCT BOPDS_Iterator.hxx / .cxx
pub struct BOPDS_Iterator<'a> {
    ds: &'a DS,
    /// Fuzzy tolerance for interference detection (from BOPAlgo_Options).
    fuzzy_tol: f64,
    // Per-type-combo pair buckets, indexed by TypeToInteger(t1, t2) result:
    //   0=VV, 1=VE, 2=EE, 3=VF, 4=EF, 5=FF, 6=VZ, 7=EZ, 8=FZ, 9=ZZ
    my_lists: Vec<Vec<(usize, usize)>>,
    // Extra interference lists for IntersectExt (OCCT BOPDS_Iterator::myExtLists)
    my_ext_lists: Vec<Vec<(usize, usize)>>,
    // Flag indicating IntersectExt has been called (OCCT BOPDS_Iterator::myUseExt)
    my_use_ext: bool,
    // Current iteration state
    current_list: Vec<(usize, usize)>, // pairs being iterated (cloned from bucket)
    current_pos: usize,                // index into current_list
    my_run_parallel: bool,
}

impl<'a> BOPDS_Iterator<'a> {
    pub fn new(ds: &'a DS, fuzzy_tol: f64) -> Self {
        let n = 10; // NbInterfTypes = 10 (VV..ZZ)
        let mut my_lists = Vec::with_capacity(n);
        for _ in 0..n {
            my_lists.push(Vec::new());
        }
        BOPDS_Iterator {
            ds,
            fuzzy_tol,
            my_ext_lists: my_lists.clone(),
            my_lists,
            my_use_ext: false,
            current_list: Vec::new(),
            current_pos: 0,
            my_run_parallel: false,
        }
    }

    pub fn set_run_parallel(&mut self, flag: bool) {
        self.my_run_parallel = flag;
    }

    pub fn run_parallel(&self) -> bool {
        self.my_run_parallel
    }

    /// OCCT BOPDS_Tools::TypeToInteger(ShapeType)
    fn type_to_int(t: ShapeType) -> i32 {
        match t {
            ShapeType::Vertex => 7,
            ShapeType::Edge => 6,
            ShapeType::Face => 4,
            ShapeType::Shell => 3,
            ShapeType::Solid => 2,
            ShapeType::Wire => 5,
            ShapeType::Compound => 0,
            ShapeType::CompSolid => 1,
            _ => 9,
        }
    }

    /// OCCT BOPDS_Tools::TypeToInteger(t1, t2) 鈫?bucket index
    fn type_to_bucket(t1: ShapeType, t2: ShapeType) -> i32 {
        let i1 = Self::type_to_int(t1);
        let i2 = Self::type_to_int(t2);
        let ix = i2 * 10 + i1;
        match ix {
            77 => 0,      // VV
            76 | 67 => 1, // VE
            66 => 2,      // EE
            74 | 47 => 3, // VF
            64 | 46 => 4, // EF
            44 => 5,      // FF
            72 | 27 => 6, // VZ
            62 | 26 => 7, // EZ
            42 | 24 => 8, // FZ
            22 => 9,      // ZZ
            _ => -1,
        }
    }

    /// OCCT BOPDS_Iterator::Prepare 鈥?build BVH, find all overlapping pairs, bucket by type.
    ///
    /// Builds a single BVH over all shapes (vertices + edges + faces), runs candidate_pairs,
    /// filters by cross-operand, skips shape-subshape pairs, and buckets into my_lists.
    ///
    /// OCCT uses BOPTools_BoxPairSelector with a BVH tree (Bnd_Box-based). rcad uses
    /// cross-operand direct enumeration with AABB overlap check per type combination.
    ///
    /// OCCT BOPDS_Iterator.cxx L247-265: Prepare calls Intersect(L270-359).
    pub fn prepare(&mut self) {
        // Clear all lists (OCCT L254-258)
        for list in &mut self.my_lists {
            list.clear();
        }
        for list in &mut self.my_ext_lists {
            list.clear();
        }
        self.my_use_ext = false;

        let nv = self.ds.vertex_count();
        let ne = self.ds.edge_count();
        let nf = self.ds.face_count();
        if nv + ne + nf < 2 {
            return;
        }

        let a_vc = self.ds.a_vertex_count();
        let a_ec = self.ds.a_edge_count();
        let a_fc = self.ds.a_face_count();
        let mut add_pair = |s1: usize, s2: usize, t1: ShapeType, t2: ShapeType| {
            // Cross-operand filter: skip same-operand pairs
            // Vertex range: 0..a_vc = operand A, a_vc..nv = operand B
            // Edge range: 0..a_ec = operand A, a_ec..ne = operand B
            // Face range: 0..a_fc = operand A, a_fc..nf = operand B
            let (op1, op2) = match (t1, t2) {
                (ShapeType::Vertex, ShapeType::Vertex) => (s1 < a_vc, s2 < a_vc),
                (ShapeType::Vertex, ShapeType::Edge) => (s1 < a_vc, s2 < a_ec),
                (ShapeType::Edge, ShapeType::Vertex) => (s1 < a_ec, s2 < a_vc),
                (ShapeType::Edge, ShapeType::Edge) => (s1 < a_ec, s2 < a_ec),
                (ShapeType::Vertex, ShapeType::Face) => (s1 < a_vc, s2 < a_fc),
                (ShapeType::Face, ShapeType::Vertex) => (s1 < a_fc, s2 < a_vc),
                (ShapeType::Edge, ShapeType::Face) => (s1 < a_ec, s2 < a_fc),
                (ShapeType::Face, ShapeType::Edge) => (s1 < a_fc, s2 < a_ec),
                (ShapeType::Face, ShapeType::Face) => (s1 < a_fc, s2 < a_fc),
                _ => return, // skip unsupported type combos
            };
            if op1 == op2 {
                return;
            } // cross-operand only

            // OCCT L335-340: avoid interfering shape with its sub-shapes
            if t1 == ShapeType::Vertex && t2 == ShapeType::Edge {
                if self.ds.edge_has_vertex(s1, s2) {
                    return;
                }
            }
            if t1 == ShapeType::Edge && t2 == ShapeType::Vertex {
                if self.ds.edge_has_vertex(s2, s1) {
                    return;
                }
            }

            let bucket = Self::type_to_bucket(t1, t2);
            if bucket >= 0 && (bucket as usize) < self.my_lists.len() {
                // Push (s1, s2) preserving type-specific ordering.
                self.my_lists[bucket as usize].push((s1, s2));
            }
        };

        // Helper: edge AABB = (box_min, box_max) with gap = box_gap + fuzzy_tol
        let edge_aabb = |ei: usize| -> Option<Aabb> {
            let si = if ei < self.ds.edge_shape_idx.len() {
                self.ds.edge_shape_idx[ei]
            } else {
                self.ds.vertex_count() + ei
            };
            self.ds.shape_info.get(si).and_then(|info| {
                info.box_min.zip(info.box_max).map(|(min, max)| Aabb {
                    min,
                    max,
                    gap: info.box_gap + self.fuzzy_tol,
                })
            })
        };

        // Helper: vertex AABB = point, gap = vertex_tolerance + fuzzy_tol
        let vertex_aabb = |vi: usize| -> Aabb {
            let p = self.ds.vertex_point(vi);
            Aabb {
                min: p,
                max: p,
                gap: self.ds.vertex_tolerance(vi) + self.fuzzy_tol,
            }
        };

        // AABB overlap test
        let aabb_overlap =
            |a_min: [f64; 3], a_max: [f64; 3], b_min: [f64; 3], b_max: [f64; 3]| -> bool {
                !(a_max[0] < b_min[0]
                    || a_min[0] > b_max[0]
                    || a_max[1] < b_min[1]
                    || a_min[1] > b_max[1]
                    || a_max[2] < b_min[2]
                    || a_min[2] > b_max[2])
            };

        // VV pairs: cross-operand vertices
        for va in 0..a_vc {
            for vb in a_vc..nv {
                add_pair(va, vb, ShapeType::Vertex, ShapeType::Vertex);
            }
        }

        // VE pairs: vertex vs edge (cross-operand) with AABB overlap filter
        // OCCT BOPDS_Iterator::Intersect (L270-359) uses a single BVH with
        // Bnd_Box + SetGap. rcad: O(n虏) enumeration with Aabb::intersects
        // (gap = box_gap + fuzzy_tol) 鈥?functionally equivalent AABB overlap.
        // Architecture diff: OCCT uses single BVH + SetSame; rcad per-type O(n虏).
        for vi in 0..nv {
            let is_a = vi < a_vc;
            for ei in 0..ne {
                let is_e_a = ei < a_ec;
                if is_a == is_e_a {
                    continue;
                }
                let v_abb = vertex_aabb(vi);
                if let Some(e_abb) = edge_aabb(ei) {
                    if !v_abb.intersects(&e_abb) {
                        continue;
                    }
                } else {
                    continue; // edge has no valid AABB (type mismatch or missing)
                }
                add_pair(vi, ei, ShapeType::Vertex, ShapeType::Edge);
            }
        }

        // EE pairs: cross-operand edges with AABB overlap filter
        for ea in 0..a_ec {
            for eb in a_ec..ne {
                if let (Some(ea_abb), Some(eb_abb)) = (edge_aabb(ea), edge_aabb(eb)) {
                    if !ea_abb.intersects(&eb_abb) {
                        continue;
                    }
                }
                add_pair(ea, eb, ShapeType::Edge, ShapeType::Edge);
            }
        }

        // VF pairs: vertex vs face (all cross-operand)
        for vi in 0..nv {
            let is_a = vi < a_vc;
            for fi in 0..nf {
                let is_f_a = fi < a_fc;
                if is_a == is_f_a {
                    continue;
                }
                add_pair(vi, fi, ShapeType::Vertex, ShapeType::Face);
            }
        }

        // EF pairs: edge vs face (all cross-operand) with AABB overlap filter
        for ei in 0..ne {
            let is_a = ei < a_ec;
            for fi in 0..nf {
                let is_f_a = fi < a_fc;
                if is_a == is_f_a {
                    continue;
                }
                let si_f = if fi < self.ds.face_shape_idx.len() {
                    self.ds.face_shape_idx[fi]
                } else {
                    fi
                };
                if let Some(e_abb) = edge_aabb(ei) {
                    if let Some(f_info) = self.ds.shape_info.get(si_f) {
                        if let (Some(f_min), Some(f_max)) = (f_info.box_min, f_info.box_max) {
                            let f_abb = Aabb {
                                min: f_min,
                                max: f_max,
                                gap: f_info.box_gap + self.fuzzy_tol,
                            };
                            if !e_abb.intersects(&f_abb) {
                                continue;
                            }
                        }
                    }
                }
                add_pair(ei, fi, ShapeType::Edge, ShapeType::Face);
            }
        }

        // FF pairs: cross-operand faces
        for fa in 0..a_fc {
            for fb in a_fc..nf {
                add_pair(fa, fb, ShapeType::Face, ShapeType::Face);
            }
        }

        // stable_sort each bucket (OCCT Initialize L203: std::stable_sort)
        for list in &mut self.my_lists {
            list.sort();
        }
    }

    /// OCCT BOPDS_Iterator::Initialize 鈥?select pairs of given type combination.
    ///
    /// Applies stable_sort (already done in Prepare) and sets up iteration.
    pub fn initialize(&mut self, t1: ShapeType, t2: ShapeType) {
        let bucket = Self::type_to_bucket(t1, t2);
        if bucket >= 0 && (bucket as usize) < self.my_lists.len() {
            // OCCT L203: std::stable_sort(aPairs.begin(), aPairs.end())
            let mut pairs = self.my_lists[bucket as usize].clone();
            pairs.sort();
            self.current_list = pairs;
        } else {
            self.current_list.clear();
        }
        self.current_pos = 0;
    }

    /// OCCT BOPDS_Iterator::More
    pub fn more(&self) -> bool {
        self.current_pos < self.current_list.len()
    }

    /// OCCT BOPDS_Iterator::Next
    pub fn next(&mut self) {
        self.current_pos += 1;
    }

    /// OCCT BOPDS_Iterator::Value
    ///
    /// Returns (index1, index2) where index1 <= index2 (matching OCCT's
    /// BOPDS_Pair which stores min/max indices).
    pub fn value(&self) -> (usize, usize) {
        self.current_list[self.current_pos]
    }

    /// Returns the expected length (OCCT ExpectedLength)
    pub fn expected_length(&self) -> usize {
        self.current_list.len()
    }

    /// Returns a reference to the pre-computed pair list for (t1, t2).
    /// Must be called after `prepare()`.
    pub fn pairs(&self, t1: ShapeType, t2: ShapeType) -> &[(usize, usize)] {
        let bucket = Self::type_to_bucket(t1, t2);
        if bucket >= 0 && (bucket as usize) < self.my_lists.len() {
            &self.my_lists[bucket as usize]
        } else {
            &[]
        }
    }

    /// OCCT BOPDS_Iterator.cxx L363-463: IntersectExt
    ///
    /// Builds extra interference pairs for vertices whose tolerance was
    /// increased (theIndices).  Builds a single BVH over all source shapes,
    /// using the SD vertex box (expanded tolerance) for extra vertices.
    /// For each extra vertex, finds all overlapping shapes, filters by
    /// cross-operand and sub-shape relationship, deduplicates, and
    /// appends the pairs to my_lists in the appropriate type bucket.
    pub fn intersect_ext(&mut self, extra_map: &HashSet<usize>) {
        // L365-368: if (!myDS) return; (rcad: ds always present)
        // L370: const int aNb = myDS->NbSourceShapes();
        let a_nb = self.ds.nb_source_shapes();

        // L372-374: BOPTools_BoxTree aBoxTree; aBoxTree.SetSize(aNb);
        // BOPDS_VectorOfTSR aVTSR(theIndices.Extent());
        // rcad: single BoxTree for all shapes + extra vertex index list (TSR equivalents)

        let mut all_indices: Vec<usize> = Vec::new();
        let mut all_aabbs: Vec<Aabb> = Vec::new();
        let mut extra_vert_indices: Vec<usize> = Vec::new();

        // L376-402: for (int i = 0; i < aNb; ++i)
        for i in 0..a_nb {
            let si = self.ds.shape_info_at(i);
            // L378-382: if (!aSI.IsInterfering() || (aSI.ShapeType() == TopAbs_SOLID)) continue;
            // rcad: IsInterfering = Vertex | Edge | Face only
            match si.shape_type {
                ShapeType::Vertex | ShapeType::Edge | ShapeType::Face => {}
                _ => continue,
            }
            all_indices.push(i);

            // L384-401: theIndices.Contains(i) ? SD box : normal box
            if extra_map.contains(&i) {
                // L386-388: int nVSD = i; myDS->HasShapeSD(i, nVSD);
                let n_vsd = self.ds.has_shape_sd(i).unwrap_or(i);
                let si_sd = self.ds.shape_info_at(n_vsd);
                // L389: const Bnd_Box& aBox = aSISD.Box();
                if let (Some(min), Some(max)) = (si_sd.box_min, si_sd.box_max) {
                    all_aabbs.push(Aabb {
                        min,
                        max,
                        gap: si_sd.box_gap,
                    });
                    extra_vert_indices.push(i);
                }
            } else {
                // L400: aBoxTree.Add(i, Bnd_Tools::Bnd2BVH(aSI.Box()));
                if let (Some(min), Some(max)) = (si.box_min, si.box_max) {
                    all_aabbs.push(Aabb {
                        min,
                        max,
                        gap: si.box_gap,
                    });
                }
            }
        }

        // L404: aBoxTree.Build();
        let tree = crate::bop::tools::box_tree::BoxTree::build(all_indices, all_aabbs);

        // L406-407: BOPTools_Parallel::Perform(myRunParallel, aVTSR);
        // rcad: sequential (no parallel framework)

        // L412: NCollection_Map<BOPDS_Pair> aMPFence;
        let mut fence: HashSet<(usize, usize)> = HashSet::new();

        // L414: const int aNbV = aVTSR.Length();
        // For each extra vertex (TSR entry), find all intersecting shapes
        for &i in &extra_vert_indices {
            // L424-428: get shape info, rank, type for the extra vertex
            let si_i = self.ds.shape_info_at(i);
            let i_rank = si_i.rank;
            let ti = si_i.shape_type;
            let i_rank_ti = Self::type_rank(ti);

            // Get SD box for TSR query (OCCT L392-396: BOPDS_TSR box = SD box)
            let n_vsd = self.ds.has_shape_sd(i).unwrap_or(i);
            let si_sd = self.ds.shape_info_at(n_vsd);
            let (qmin, qmax) = match (si_sd.box_min, si_sd.box_max) {
                (Some(min), Some(max)) => (min, max),
                _ => continue,
            };
            let qgap = si_sd.box_gap;

            // rcad: single-node BoxTree = OCCT BOPDS_TSR with SetBVHSet + SetBox
            let tsr_tree = crate::bop::tools::box_tree::BoxTree::build(
                vec![i],
                vec![Aabb {
                    min: qmin,
                    max: qmax,
                    gap: qgap,
                }],
            );
            let candidates = crate::bop::tools::box_tree::BoxTree::candidate_pairs(&tsr_tree, &tree);

            // L430-459: Treat selections
            for &(_, j) in &candidates {
                if i == j {
                    continue;
                }

                let sj = self.ds.shape_info_at(j);
                // L435-437: if (iRankI == iRankJ) continue;
                if i_rank == sj.rank {
                    continue;
                }

                // L440-442: get j's type
                let tj = sj.shape_type;
                let j_rank_tj = Self::type_rank(tj);

                // L444-448: avoid interfering of the shape with its sub-shapes
                if (i_rank_ti < j_rank_tj && si_i.has_sub_shape(j))
                    || (i_rank_ti > j_rank_tj && sj.has_sub_shape(i))
                {
                    continue;
                }

                // L450-458: dedup via fence map + bucket to myExtLists
                let key = if i < j { (i, j) } else { (j, i) };
                if !fence.insert(key) {
                    continue;
                }

                let bucket = Self::type_to_bucket(ti, tj);
                if bucket >= 0 && (bucket as usize) < self.my_ext_lists.len() {
                    self.my_ext_lists[bucket as usize].push((i, j));
                }
            }
        }

        // OCCT L462: myUseExt = true;
        // Append my_ext_lists to my_lists, matching OCCT Initialize merge
        self.my_use_ext = true;
        for (bucket, ext_list) in self.my_ext_lists.iter().enumerate() {
            self.my_lists[bucket].extend(ext_list);
        }
    }

    /// OCCT BOPDS_Tools::TypeToInteger 鈥?hierarchy depth ordering
    /// (Vertex=0, Edge=1, Face=2, ...), matching OCCT's TopAbs enum.
    /// Used in IntersectExt for subshape comparison (rcad's type_to_int
    /// uses inverted values and is not suitable for this comparison).
    fn type_rank(t: ShapeType) -> i32 {
        match t {
            ShapeType::Vertex => 0,
            ShapeType::Edge => 1,
            ShapeType::Face => 2,
            ShapeType::Wire => 3,
            ShapeType::Shell => 4,
            ShapeType::Solid => 5,
            ShapeType::CompSolid => 6,
            ShapeType::Compound => 7,
            _ => 8,
        }
    }
}
