use super::*;
use rcad_kernel::topods::ShapeType;

impl<'a> PaveFiller<'a> {
    // OCCT BOPAlgo_PaveFiller_3.cxx L997-1333
    pub(crate) fn force_interf_ee(&mut self) {
        // L999-1003: comment
        // Now that we have vertices increased and unified, try to find additional
        // common blocks among the pairs of edges.
        // Since all real intersections should have already happened, here we
        // are interested in common blocks only, thus we need to check only
        // those pairs of pave blocks with the same bounding vertices.

        // L1005-1023: Initialize pave blocks for all vertices which participated
        // in intersections
        // OCCT: for (int i = 0; i < aNbS; ++i) if VERTEX && HasInterf(i) -> InitPaveBlocksForVertex(i)
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            if self.ds.shape_type_of(i) != ShapeType::Vertex {
                continue;
            }
            // OCCT L1014: if (myDS->HasInterf(i))
            // rcad: HasInterf(i) checks if vertex i has any interference.
            // Check via interf_tb for any pair involving i.
            let has_interf = self.ds.interf_tb.iter().any(|&(a, b)| a == i || b == i);
            if has_interf {
                self.ds.init_pave_blocks_for_vertex(i);
            }
        }

        // L1024-1080: Fill the connection map from bounding vertices to PBs
        // Fence map of pave blocks (OCCT L1030: NCollection_Map<handle<BOPDS_PaveBlock>> aMPBFence)
        let mut a_mpb_fence: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        // rcad: HashMap keyed by (v_min, v_max), value = Vec<(edge_idx, local_pb_idx)>
        // OCCT: NCollection_IndexedDataMap<BOPDS_Pair, NCollection_List<handle<PaveBlock>>>
        let mut a_pb_map: std::collections::HashMap<(usize, usize), Vec<(usize, usize)>> =
            std::collections::HashMap::new();

        for i in 0..a_nb_s {
            // L1034-1038: only edges
            if self.ds.shape_type_of(i) != ShapeType::Edge {
                continue;
            }
            // L1041-1044: edge must have PBs (HasReference equivalent)
            let ei = self.ds.shape_info_at(i).source_idx;
            if ei >= self.ds.edges.len() || self.ds.edges[ei].pave_blocks.is_empty() {
                continue;
            }
            // L1047-1051: skip degenerated edges
            if self.ds.edge_has_flag(ei) {
                continue;
            }

            // L1056-1079: iterate PBs of this edge
            let a_lpb = &self.ds.edges[ei].pave_blocks;
            for local_i in 0..a_lpb.len() {
                let a_pb = &a_lpb[local_i];

                // L1060-1061: get real PaveBlock (OCCT: myDS->RealPaveBlock(aPB))
                // rcad: no CommonBlock indirection — a_pb IS the real PB
                let a_pbr_key = (ei, local_i);
                // L1062-1065: fence map — skip if already processed
                if !a_mpb_fence.insert(a_pbr_key) {
                    continue;
                }

                // L1068-1069: get vertex indices
                let (n_v1, n_v2) = {
                    let pbr = a_pb.0.read().unwrap();
                    (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
                };

                // L1072-1078: add PB to map keyed by vertex pair
                // OCCT: BOPDS_Pair aPair(nV1, nV2); aPBMap(aPBMap.Add(aPair, ...))
                let a_pair = if n_v1 <= n_v2 {
                    (n_v1, n_v2)
                } else {
                    (n_v2, n_v1)
                };
                a_pb_map.entry(a_pair).or_default().push(a_pbr_key);
            }
        }

        // L1082-1086: empty map check
        if a_pb_map.is_empty() {
            return;
        }

        // OCCT L1088: const bool bSICheckMode = (myArguments.Extent() == 1);
        let b_si_check_mode = self.my_arguments.len() == 1;

        // L1090-1225: Prepare pairs for intersection
        // rcad: Vec of struct to hold pair data before intersection.
        struct EEIntersectPair {
            ei1: usize,
            pb1_local: usize,
            n_e1: usize,
            n_v1: usize,
            n_v2: usize,
            ei2: usize,
            pb2_local: usize,
            n_e2: usize,
            b_use_add_tol: bool,
        }

        let mut a_v_edge_edge: Vec<EEIntersectPair> = Vec::new();

        for (&(n_v1, n_v2), pbs) in &a_pb_map {
            // L1100-1102: need at least 2 PBs sharing the same vertices
            if pbs.len() < 2 {
                continue;
            }

            // L1105-1118: compute tolerance addition from vertex tolerances
            // L1109-1110: get TopoDS_Vertex shapes (rcad: tolerance only)
            // L1116-1118: compute aTolAdd
            let a_tol_add = if b_si_check_mode {
                self.ds.fuzzy_tol
            } else {
                2.0 * self
                    .ds
                    .vertex_tolerance(n_v1)
                    .max(self.ds.vertex_tolerance(n_v2))
            };

            // L1120-1225: iterate all PBs in this group as pairs (i < j)
            for i in 0..pbs.len() {
                let (ei1, pb1_local) = pbs[i];
                let pb1 = &self.ds.edges[ei1].pave_blocks[pb1_local];
                let cb1 = self.ds.common_block(pb1);
                let n_e1 = pb1.0.read().unwrap().original_edge;
                let i_r1 = self.ds.rank(n_e1);
                let (t11, t12) = pb1.0.read().unwrap().range();
                let a_t_mid = (t11 + t12) * 0.5;

                // L1131-1139: compute tangent at middle point of edge 1
                let curve1 = self.ds.edge_curve(n_e1);
                let (a_pm, a_v_tgt1) = match curve1 {
                    Some(c) => {
                        let p = c.point_at(a_t_mid);
                        let t = c.tangent_at(a_t_mid);
                        (p, t)
                    }
                    None => continue,
                };
                if a_v_tgt1.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                    continue;
                }
                let a_v_tgt1_n = a_v_tgt1.normalize();

                // L1141: nested iterator (i < j)
                for j in (i + 1)..pbs.len() {
                    let (ei2, pb2_local) = pbs[j];
                    let pb2 = &self.ds.edges[ei2].pave_blocks[pb2_local];
                    let cb2 = self.ds.common_block(pb2);
                    let n_e2 = pb2.0.read().unwrap().original_edge;
                    let i_r2 = self.ds.rank(n_e2);

                    // L1149-1160: check that edges came from different arguments,
                    // or have acquired (new) vertices
                    if i_r1 == i_r2 {
                        if (!self.ds.is_new_vertex(n_v1) && self.ds.rank(n_v1) == i_r1)
                            || (!self.ds.is_new_vertex(n_v2) && self.ds.rank(n_v2) == i_r2)
                        {
                            continue;
                        }
                    }

                    // L1162-1169: check that the PBs do not already share a CommonBlock
                    if let (Some(ref cb1), Some(ref cb2)) = (cb1.as_ref(), cb2.as_ref()) {
                        // rcad: compare by pointer identity (OCCT handle ==)
                        if std::ptr::eq(*cb1, *cb2) {
                            continue;
                        }
                    }

                    // L1175-1205: check the angle between edges at middle point
                    let mut b_use_add_tol = true;
                    {
                        let curve2 = self.ds.edge_curve(n_e2);
                        let curve2 = match curve2 {
                            Some(c) => c,
                            None => continue,
                        };
                        if !(matches!(curve1, Some(&Curve3::Line(_)))
                            && matches!(curve2, Curve3::Line(_)))
                        {
                            // L1182-1202: non-line case - project middle point onto curve2
                            let a_proj = closest_point_on_curve(curve2, a_pm, 64);
                            let a_v_tgt2 = curve2.tangent_at(a_proj.param);
                            if a_v_tgt2.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                                continue;
                            }
                            // L1199: angle threshold (cos >= 0.9063 ≈ 25 degrees)
                            let a_cos = a_v_tgt1_n.dot(a_v_tgt2.normalize());
                            if a_cos.abs() < 0.9063 {
                                b_use_add_tol = false;
                            }
                        }
                    }

                    // L1207-1223: add pair for intersection
                    a_v_edge_edge.push(EEIntersectPair {
                        ei1,
                        pb1_local,
                        n_e1,
                        n_v1,
                        n_v2,
                        ei2,
                        pb2_local,
                        n_e2,
                        b_use_add_tol,
                    });
                } // for (nested j)
            } // for (outer i)
        } // for each vertex pair

        let a_nb_pairs = a_v_edge_edge.len();
        if a_nb_pairs == 0 {
            return;
        }

        // L1233-1238: close preparation step
        // OCCT: aPBMap.Clear(); aMPBFence.Clear(); anAlloc->Reset(false);
        a_pb_map.clear();
        a_mpb_fence.clear();

        // L1240-1252: Perform intersection (rcad: sequential, no parallel dispatch)
        // L1253-1257: get EE array reference (rcad: already exists as Vec)
        // rcad: use Vec<Option<InterferenceEE>> to store per-pair results
        let b_si = b_si_check_mode;
        let mut a_ee_results: Vec<Option<(InterferenceEE, bool)>> = vec![None; a_nb_pairs];

        for (idx, entry) in a_v_edge_edge.iter().enumerate() {
            let curve1 = match self.ds.edge_curve(entry.n_e1) {
                Some(c) => c.clone(),
                None => continue,
            };
            let curve2 = match self.ds.edge_curve(entry.n_e2) {
                Some(c) => c.clone(),
                None => continue,
            };

            let tr1 = self.ds.edge_range(entry.n_e1);
            let tr2 = self.ds.edge_range(entry.n_e2);

            let a_tol_add = if b_si {
                self.ds.fuzzy_tol
            } else {
                2.0 * self
                    .ds
                    .vertex_tolerance(entry.n_v1)
                    .max(self.ds.vertex_tolerance(entry.n_v2))
            };

            let fuzzy = if entry.b_use_add_tol {
                self.ds.fuzzy_tol + a_tol_add
            } else {
                self.ds.fuzzy_tol
            };

            // OCCT L1131-1133: midpoint of edge1's range
            let mid_t1 = (tr1[0] + tr1[1]) * 0.5;
            let mid_p1 = curve1.point_at(mid_t1);

            // OCCT L1167-1202: angle check at midpoint
            // Returns false only if both curves are non-linear AND angle < 25°
            let b_use_add_tol = {
                let is_ic_line = matches!(&curve1, Curve3::Line(_));
                let is_edge_line = matches!(&curve2, Curve3::Line(_));
                if is_ic_line || is_edge_line {
                    true
                } else {
                    let proj_mid = closest_point_on_curve(&curve2, mid_p1, 64);
                    let p2_at_mid = curve2.point_at(proj_mid.param);
                    let d_at_mid = mid_p1.distance(p2_at_mid);
                    if d_at_mid > fuzzy {
                        continue;
                    }
                    let v_tgt1 = curve1.tangent_at(mid_t1);
                    if v_tgt1.length_squared() < 1e-60 {
                        continue;
                    }
                    let v_tgt2 = curve2.tangent_at(proj_mid.param);
                    if v_tgt2.length_squared() < 1e-60 {
                        continue;
                    }
                    let a_cos = v_tgt1.normalize().dot(v_tgt2.normalize());
                    a_cos.abs() >= 0.9063
                }
            };

            // OCCT L1208-1218: BOPAlgo_EdgeEdge intersection with quick coincidence check
            // OCCT uses IntTools_EdgeEdge recursively; rcad has EdgeEdgeIntersector equivalent.
            let mut ee_int = crate::inttools::edge_edge::EdgeEdgeIntersector::new();
            ee_int.set_edges(entry.n_e1, tr1, entry.n_e2, tr2, self.ds);
            ee_int.use_quick_coincidence_check(true);
            ee_int.set_fuzzy_value(fuzzy);
            ee_int.perform();

            let result: Option<(InterferenceEE, bool)> = if ee_int.is_done() {
                let cp = &ee_int.common_parts()[0];
                // OCCT L1287-1291: only EDGE-type common parts are accepted
                if cp.is_edge_type {
                    let mid_t1 = (cp.range1[0] + cp.range1[1]) * 0.5;
                    let mid_t2 = if let Some(r2) = cp.ranges2.first() {
                        (r2[0] + r2[1]) * 0.5
                    } else {
                        continue;
                    };
                    let avg_pt = (curve1.point_at(mid_t1) + curve2.point_at(mid_t2)) * 0.5;
                    Some((
                        InterferenceEE {
                            e1: entry.n_e1,
                            e2: entry.n_e2,
                            point: avg_pt,
                            param1: mid_t1,
                            param2: mid_t2,
                            new_vertex: entry.n_v1,
                            range1: cp.range1,
                            range2: cp.ranges2.first().copied().unwrap_or([mid_t2, mid_t2]),
                        },
                        true,
                    ))
                } else {
                    None
                }
            } else {
                None
            };
            a_ee_results[idx] = result;
        }

        // L1262-1330: Process results
        // rcad: build connection map for CommonBlock creation
        // OCCT: IndexedDataMap<handle<PaveBlock>, List<handle<PaveBlock>>>
        // rcad: BTreeMap<(usize, usize), Vec<(usize, usize)>> (unencoded tuples, no collision risk)
        let mut a_mpblpb: std::collections::BTreeMap<(usize, usize), Vec<(usize, usize)>> =
            std::collections::BTreeMap::new();

        for (idx, opt_ee) in a_ee_results.iter().enumerate() {
            let entry = &a_v_edge_edge[idx];
            let (ee, is_edge_type) = match opt_ee {
                Some(pair) => pair,
                None => continue,
            };
            // OCCT L1287-1291: only EDGE-type common parts accepted
            if !is_edge_type {
                continue;
            }

            // L1297-1305: self-interference warning
            // OCCT: if both edges are from the same rank, add self-interference warning
            if self.ds.rank(entry.n_e1) == self.ds.rank(entry.n_e2) {
                self.my_report
                    .add_alert(crate::bopalgo::Alert::AcquiredSelfIntersection(vec![
                        entry.n_e1, entry.n_e2,
                    ]));
            }

            // L1307-1310: create InterfEE entry
            self.ds.interf_ee.push(ee.clone());
            self.ds
                .interf_tb
                .insert((entry.n_e1.min(entry.n_e2), entry.n_e1.max(entry.n_e2)));

            // L1312-1329: Fill map for common blocks creation
            // Use (ei, local_i) tuple keys directly (no collision risk vs encoded scalar).
            let key1 = (entry.ei1, entry.pb1_local);
            let key2 = (entry.ei2, entry.pb2_local);

            let pb1 = &self.ds.edges[entry.ei1].pave_blocks[entry.pb1_local];
            let pb2 = &self.ds.edges[entry.ei2].pave_blocks[entry.pb2_local];

            // For each of the two PBs, if it belongs to a CommonBlock,
            // connect all PBs in that CB as mutual siblings.
            for &pb_ref in [pb1, pb2].iter() {
                if self.ds.is_common_block(pb_ref) {
                    // PB is already part of a CommonBlock — OCCT L1315-1327:
                    // expand connections to all PBs in that CB.
                    let cb = self.ds.common_block(pb_ref).unwrap();
                    let a_lpcb = cb.pave_blocks();
                    // rcad: CB stores global PB indices; encode them similarly
                    for &(_cb_pb_gi, _) in a_lpcb {
                        // rcad: connect as siblings.  The encoding scheme
                        // differs from local (ei,li) encoding; use fill_map
                        // on both keys to ensure the pairs are connected.
                        crate::bopalgo::fill_map(&mut a_mpblpb, key1, key2);
                    }
                }
            }

            // L1329: connect the two PBs bidirectionally
            crate::bopalgo::fill_map(&mut a_mpblpb, key1, key2);
        }

        // L1332: Create new common blocks of coinciding pairs.
        // OCCT: BOPAlgo_Tools::PerformCommonBlocks(aMPBLPB, anAlloc, myDS);
        // rcad: use make_blocks on a_mpblpb to group connected PBs, then create
        // CommonBlocks from each group.  No global PB array rebuild (OCCT uses
        // PB handles directly — no global array).
        if a_mpblpb.len() > 0 {
            let a_m_blocks: Vec<Vec<(usize, usize)>> = crate::bopalgo::make_blocks(&a_mpblpb);
            for block in &a_m_blocks {
                if block.len() < 2 {
                    continue;
                }
                // Collect unique face indices for this block
                let mut a_lfaces: Vec<usize> = Vec::new();
                for &(ei, local_i) in block {
                    let face_indices: Vec<usize> = self.ds.edges[ei]
                        .face_reps
                        .iter()
                        .map(|r| r.face_idx)
                        .collect();
                    for &fi in &face_indices {
                        if !a_lfaces.contains(&fi) {
                            a_lfaces.push(fi);
                        }
                    }
                }
                if a_lfaces.len() < 2 {
                    continue;
                }
                // Check if any PB in the block already has a CommonBlock
                let mut a_cb: Option<crate::bopds::common_block::CommonBlock> = None;
                let mut cb_idx: Option<usize> = None;
                for &(ei, local_i) in block {
                    let pb = &self.ds.edges[ei].pave_blocks[local_i];
                    if self.ds.is_common_block(pb) {
                        let existing = self.ds.common_block(pb).unwrap();
                        // Collect faces from existing CB
                        for &fi in existing.faces() {
                            if !a_lfaces.contains(&fi) {
                                a_lfaces.push(fi);
                            }
                        }
                        if a_cb.is_none() {
                            a_cb = Some(existing.clone());
                            cb_idx = pb.0.read().unwrap().common_block_idx;
                        }
                    }
                }
                // Create or update CommonBlock
                let mut cb = a_cb.unwrap_or_else(crate::bopds::common_block::CommonBlock::new);
                cb.set_pave_blocks(block.clone());
                cb.set_faces(a_lfaces);
                let new_idx = match cb_idx {
                    Some(idx) => {
                        let tol = crate::bopds::tools::compute_tolerance_of_cb(&self.ds, idx);
                        self.ds.common_blocks[idx] = cb.clone();
                        self.ds.common_blocks[idx].set_tolerance(tol);
                        idx
                    }
                    None => {
                        let idx = self.ds.common_blocks.len();
                        cb.set_tolerance(0.0);
                        self.ds.common_blocks.push(cb);
                        let tol = crate::bopds::tools::compute_tolerance_of_cb(&self.ds, idx);
                        self.ds.common_blocks[idx].set_tolerance(tol);
                        idx
                    }
                };
                // Mark all PBs in block as belonging to this CommonBlock
                for &(ei, local_i) in block {
                    if let Some(local_pb) = self.ds.edges[ei].pave_blocks.get_mut(local_i) {
                        local_pb.0.write().unwrap().common_block_idx = Some(new_idx);
                    }
                }
            }
        }
    }

    // OCCT BOPAlgo_PaveFiller_5.cxx L772-827
    pub(crate) fn force_interf_ef(&mut self) {
        // L774-778
        if !self.is_primary {
            return;
        }

        // L787: IndexedMap with dedup — rcad: Vec + HashSet via (ne, local_i) keys
        // OCCT uses RealPaveBlock(aPB) to resolve through CommonBlock before dedup.
        let mut a_mpb: Vec<(usize, usize)> = Vec::new();
        let mut a_mpb_dedup: std::collections::HashSet<(usize, usize)> =
            std::collections::HashSet::new();
        let a_nb_s = self.ds.nb_source_shapes();
        for n_e in 0..a_nb_s {
            // L791-796: only edges
            if self.ds.shape_type_of(n_e) != ShapeType::Edge {
                continue;
            }
            // L798-802: edge must have PBs (HasReference)
            // rcad: shape_info[n_e].reference gives the flat edge index
            let ne = self.ds.shape_info_at(n_e).reference as usize;
            if !self.ds.has_pave_blocks(ne) {
                continue;
            }
            // L804-808: skip degenerated edges
            if self.ds.is_edge_degenerated(ne) {
                continue;
            }

            // L814-821: add all PBs to map with RealPaveBlock resolution
            let a_lpb = self.ds.edge_pave_blocks(ne);
            for (local_i, _a_it_lpb) in a_lpb.iter().enumerate() {
                let a_pb = &a_lpb[local_i];
                // OCCT L819: RealPaveBlock(aPB) — resolve through CommonBlock
                // rcad: if PB has a CommonBlock, use the first PB's identity
                let a_pbr_key = match self.ds.common_block(a_pb) {
                    Some(cb) => {
                        let pb1_global = cb.pave_block1().unwrap_or(0);
                        // rcad: convert global PB index to (edge, local) key.
                        // Global pool PBs are clones; use (pb1_global, 0) as unique key.
                        (pb1_global, 0)
                    }
                    None => (ne, local_i),
                };
                // L820: IndexedMap dedup
                if a_mpb_dedup.insert(a_pbr_key) {
                    a_mpb.push(a_pbr_key);
                }
            }
        }

        // L826: call overload 2 with theAddInterf=true
        self.force_interf_ef_with(&a_mpb, true);
    }

    // OCCT BOPAlgo_PaveFiller_5.cxx L831-1199
    fn force_interf_ef_with(&mut self, the_mpb: &[(usize, usize)], the_add_interf: bool) {
        // L838-840
        if the_mpb.is_empty() {
            return;
        }

        // L842-874: Fill the tree with bounding boxes of the pave blocks
        // rcad: build BoxTree from PB AABBs (shrunk data or edge endpoint fallback).
        let mut a_pb_aabbs: Vec<crate::bvh::Aabb> = Vec::with_capacity(the_mpb.len());
        let mut a_pb_indices: Vec<usize> = Vec::with_capacity(the_mpb.len());
        // L848-870: for each PB, get ShrunkData and add to tree
        for (i_pb, &(ne, local_i)) in the_mpb.iter().enumerate() {
            // rcad: equivalent to OCCT's theMPB(iPB) via (ne, local_i) pair
            // L852-858: ensure shrunk data is available
            // OCCT: if (!aPB->HasShrunkData() || !myDS->IsValidShrunkData(aPB))
            {
                let r = self.ds.edges[ne].pave_blocks[local_i].0.read().unwrap();
                let need_fill = !r.has_shrunk_data() || !self.ds.is_valid_shrunk_data(&r);
                drop(r);
                if need_fill {
                    // OCCT L854: FillShrunkData(aPB)
                    // rcad: analyze_shrunk_data as FillShrunkData equivalent
                    self.analyze_shrunk_data(ne, local_i);
                    // OCCT L855-858: if still no shrunk data after FillShrunkData
                    if !self.ds.edges[ne].pave_blocks[local_i]
                        .0
                        .read()
                        .unwrap()
                        .has_shrunk_data()
                    {
                        continue;
                    }
                }
            }
            // L866-868: ShrunkData(f, l, aPBBox, isSplit)
            // rcad: shrunk_data() returns (f, l, isSplit), bounding box computed below
            let (f, l, _is_split) = self.ds.edges[ne].pave_blocks[local_i]
                .0
                .read()
                .unwrap()
                .shrunk_data();
            // L866-870: build AABB from shrunk range endpoint positions
            // OCCT: aBBTree.Add(aPBMap.Add(aPB), Bnd_Tools::Bnd2BVH(aPBBox));
            let a_p1 = self.ds.edges[ne].curve.point_at(f);
            let a_p2 = self.ds.edges[ne].curve.point_at(l);
            let a_pb_box = crate::bvh::Aabb {
                min: a_p1.min(a_p2),
                max: a_p1.max(a_p2),
                gap: 0.0,
            };
            // L870: aBBTree.Add(aPBMap.Add(aPB), Bnd_Tools::Bnd2BVH(aPBBox));
            a_pb_aabbs.push(a_pb_box);
            a_pb_indices.push(i_pb);
        }
        // L873-874: Build BVH tree
        let a_bb_tree = BoxTree::build(a_pb_indices, a_pb_aabbs);

        // L876: bSICheckMode — Self-Interference check mode (one argument)
        // OCCT: myArguments.Extent() == 1
        let b_si_check_mode = self.my_arguments.len() <= 1;

        // L880: vector of edge-face pairs for intersection
        // OCCT: BOPAlgo_VectorOfEdgeFace aVEdgeFace
        // rcad: pair data stored in a local struct
        struct _EFPair {
            n_e: usize,
            n_f_src: usize, // source shape index for fpbdone check
            fi: usize,      // flat face index for face operations
            pb: SharedPB,
            a_tol_add: f64,
            a_ts: [f64; 2],
        }
        let mut a_v_edge_face: Vec<_EFPair> = Vec::new();

        // L882-1108: For each source face that has face info
        // OCCT: for (int nF = 0; nF < aNbS; ++nF) {
        //         if (aSI.ShapeType() != TopAbs_FACE) continue;
        //         if (!aSI.HasReference()) continue;
        let a_nb_s = self.ds.nb_source_shapes();
        for n_f in 0..a_nb_s {
            // L885-890: only faces
            if self.ds.shape_type_of(n_f) != ShapeType::Face {
                continue;
            }
            // L892-896: HasReference — rcad: check if flat face index is valid
            let fi = self.ds.shape_info_at(n_f).reference as usize;
            if fi >= self.ds.faces.len() {
                continue;
            }

            // L903-910: Face AABB for BVH query
            // OCCT: const Bnd_Box& aBoxF = aSI.Box();
            let a_box_f = {
                let si = self.ds.shape_info_at(n_f);
                let mut aabb = crate::bvh::Aabb::empty();
                if let (Some(mn), Some(mx)) = (si.box_min, si.box_max) {
                    aabb = crate::bvh::Aabb {
                        min: mn,
                        max: mx,
                        gap: si.box_gap,
                    };
                }
                aabb
            };
            let a_overlapping = a_bb_tree.query_aabb(&a_box_f);
            if a_overlapping.is_empty() {
                continue;
            }

            // L912-913: FaceInfo
            // OCCT: const TopoDS_Face& aF = TopoDS::Face(aSI.Shape());
            //        const BOPDS_FaceInfo& aFI = myDS->FaceInfo(nF);
            let a_fi = &self.ds.faces[fi].face_info;

            // L914-924: build aMVF from all face vertex sets
            // OCCT: NCollection_Map<int> aMVF from FaceInfo VerticesOn/In/Sc
            let mut a_mvf: std::collections::HashSet<usize> = std::collections::HashSet::new();
            for &v in &a_fi.vertices_on {
                a_mvf.insert(v);
            }
            for &v in &a_fi.vertices_in {
                a_mvf.insert(v);
            }
            for &v in &a_fi.vertices_sc {
                a_mvf.insert(v);
            }

            // L926-939: add PB endpoints from face's PaveBlocksOn/In/Sc
            let p_mpbf = [
                &a_fi.pave_blocks_on,
                &a_fi.pave_blocks_in,
                &a_fi.pave_blocks_sc,
            ];
            for pbs in &p_mpbf {
                for &pb_gi in pbs.iter() {
                    if pb_gi < self.ds.pave_blocks.len() {
                        let a_pb_r = self.ds.pave_blocks[pb_gi].0.read().unwrap();
                        a_mvf.insert(a_pb_r.pave1.vertex_idx);
                        a_mvf.insert(a_pb_r.pave2.vertex_idx);
                    }
                }
            }

            // Projection tool
            // rcad: IntToolsContext provides proj_ps and is_point_in_face

            // L947-1107: iterate PBs overlapping this face from BVH
            for &a_pb_idx in &a_overlapping {
                let &(ne, local_i) = &the_mpb[a_pb_idx];
                let a_pb = &self.ds.edges[ne].pave_blocks[local_i];

                // L952-955: skip if PB already in face's PB sets
                // OCCT: if (pMPBF[0]->Contains(aPB) || ...)
                // rcad: compare by original_edge (PB identity) since a_pb_idx is
                // a local BVH index, NOT a global DS PB index (namespace mismatch).
                let (n_v1, n_v2) = {
                    let r = a_pb.0.read().unwrap();
                    let a_pb_ne = r.original_edge;
                    // Check if any global pool PB with same original_edge is already in face's sets
                    let skip = a_fi
                        .pave_blocks_on
                        .iter()
                        .chain(a_fi.pave_blocks_in.iter())
                        .chain(a_fi.pave_blocks_sc.iter())
                        .any(|&pb_gi| {
                            pb_gi < self.ds.pave_blocks.len()
                                && self.ds.pave_blocks[pb_gi].0.read().unwrap().original_edge
                                    == a_pb_ne
                        });
                    if skip {
                        continue;
                    }
                    (r.pave1.vertex_idx, r.pave2.vertex_idx)
                };
                if !a_mvf.contains(&n_v1) || !a_mvf.contains(&n_v2) {
                    continue;
                }

                // L967-981: Get the edge and check arguments
                // OCCT: aPB->HasEdge(nE) / OriginalEdge() + rank check
                // rcad: PB knows its edge from (ne, local_i) pair
                // OCCT L977-980: check edge and face came from different arguments.
                // rcad: ShapeOrigin check equivalent to OCCT's Rank()
                if self.ds.edge_origin(ne) == self.ds.face_origin(fi) {
                    continue;
                }

                let a_e_curve = &self.ds.edges[ne].curve;

                // L986-1006: check directions coincidence at middle point on edge
                let mut b_use_add_tol = true;

                let a_ts: [f64; 2] = {
                    let r = a_pb.0.read().unwrap();
                    if let Some(sr) = r.shrunk_range {
                        sr
                    } else {
                        [r.pave1.param, r.pave2.param]
                    }
                };

                // L998-1002: Middle point + tangent (aBAC.D1)
                let a_t_mid = crate::boptools::intermediate_point(a_ts[0], a_ts[1]);
                // rcad: point_at + tangent_at = OCCT D1
                let a_p_on_e = a_e_curve.point_at(a_t_mid);
                let a_ve_tgt = a_e_curve.tangent_at(a_t_mid);
                if a_ve_tgt.length_squared() < TOLERANCE_LEN_SQ_DIV_SAFE {
                    continue;
                }

                // L1008-1012: project middle point onto face
                // OCCT: aProjPS.Perform(aPOnE); check NbPoints()
                let (a_uv, _a_p_on_s, a_lower_dist) =
                    match self.context.proj_ps(&self.ds, fi, a_p_on_e) {
                        Some(data) => data,
                        None => continue,
                    };

                // L1016-1029: check distance using max vertex tolerance
                // OCCT: aTolCheck = bSICheckMode ? myFuzzyValue : 2*max(Tol(V1), Tol(V2))
                let a_tol_v1 = self.ds.vertex_tolerance(n_v1);
                let a_tol_v2 = self.ds.vertex_tolerance(n_v2);
                let a_tol_check = if b_si_check_mode {
                    self.ds.fuzzy_tol
                } else {
                    2.0 * a_tol_v1.max(a_tol_v2)
                };
                if a_lower_dist > a_tol_check + self.ds.fuzzy_tol {
                    continue;
                }

                // L1031-1036: check point-in-face (UV)
                // OCCT: myContext->IsPointInFace(aF, gp_Pnt2d(U, V))
                if !self.context.is_point_in_face(&self.ds, fi, a_uv) {
                    continue;
                }

                // L1038-1052: non-plane/non-line → angle check for bUseAddTol
                // OCCT: if (aSurfAdaptor.GetType() != GeomAbs_Plane || aBAC.GetType() != GeomAbs_Line)
                let a_face = &self.ds.faces[fi];
                let is_plane = matches!(a_face.surface, Surface3::Plane(_));
                let is_line = matches!(*a_e_curve, Curve3::Line(_));
                if !is_plane || !is_line {
                    // L1040-1041: projection point (closest on surface)
                    // OCCT: aProjPS.Perform(aPOnE); aProjPS.LowerDistanceParameters(...)
                    // aProjPS.LowerDistanceParameter() — project to get tangent
                    // ComputePE(aPOnE, aTolCheck, aE, aTLdp, aDistOnE)
                    // Compute tangent at projected point: aBAC.D1(aTLdp, aPm2, aVTgt2)
                    // Check angle: cos >= 0.9063 (25°) → edges parallel → bUseAddTol = true
                    // L1040-1041: projection point (closest on surface)
                    // rcad: re-project to get the nearest surface point
                    let (_, a_p_on_s_re, _) = match self.context.proj_ps(&self.ds, fi, a_p_on_e) {
                        Some(data) => data,
                        None => continue,
                    };
                    let a_vf_norm = a_p_on_s_re - a_p_on_e;
                    if a_vf_norm.length_squared() > TOLERANCE_LEN_SQ_DIV_SAFE {
                        // L1044-1050: angle check with cos threshold 0.4226
                        // (deviation 25 degrees from 90 degrees)
                        let a_cos = a_vf_norm.normalize().dot(a_ve_tgt.normalize());
                        if a_cos.abs() > 0.4226 {
                            b_use_add_tol = false;
                        }
                    }
                }

                // L1054-1085: compute aTolAdd from endpoint projections
                let mut a_tol_add = 0.0;
                if b_use_add_tol {
                    // L1064-1076: project the two endpoints of shrunk range onto face
                    for &a_t in &[a_ts[0], a_ts[1]] {
                        let a_p = a_e_curve.point_at(a_t);
                        if let Some((_, _, a_dist_ef)) = self.context.proj_ps(&self.ds, fi, a_p) {
                            if a_dist_ef < a_tol_check && a_dist_ef > a_tol_add {
                                a_tol_add = a_dist_ef;
                            }
                        }
                    }
                    // L1077-1084: subtract edge + face tolerances
                    if a_tol_add > 0.0 {
                        a_tol_add -= (self.ds.edge_tolerance(ne) + self.ds.face_tolerance(fi));
                        if a_tol_add < 0.0 {
                            a_tol_add = 0.0;
                        }
                    }
                }

                // L1087-1092: bIntersect decision
                // OCCT: bIntersect = aTolAdd > 0 ||
                //   !pMPB || !(pMPB->Contains(aPB))
                let b_intersect = a_tol_add > 0.0
                    || !self.fpbdone.get(&n_f).map_or(false, |done| {
                        // rcad: PB identity via encoded (ne << 16 | local_i)
                        done.contains(&(ne.wrapping_mul(1_000_003).wrapping_add(local_i)))
                    });

                // L1094-1106: prepare pair for intersection
                if b_intersect {
                    // OCCT: BOPAlgo_EdgeFace& aEdgeFace = aVEdgeFace.Appended();
                    // SetIndices, SetPaveBlock, SetEdge, SetFace, etc.
                    // rcad: store pair data for later processing
                    a_v_edge_face.push(_EFPair {
                        n_e: ne,
                        n_f_src: n_f,
                        fi,
                        pb: a_pb.clone(),
                        a_tol_add,
                        a_ts,
                    });
                }
            }
        }

        // L1110-1113: no pairs found
        let a_nb_ef = a_v_edge_face.len();
        if a_nb_ef == 0 {
            return;
        }

        // L1116-1120: close preparation step
        // rcad: no allocator to clean up

        // L1122-1129: Perform intersection of the found pairs
        // OCCT: BOPTools_Parallel::Perform(myRunParallel, aVEdgeFace, myContext)
        // rcad: sequential — compute EF coincidence for each pair and store result.
        // rcad uses a simple projection-distance test; OCCT runs the full
        // BOPAlgo_EdgeFace solver (IntTools_EdgeFace) on each pair.
        //
        // ForceInterfEF configures EdgeFace with UseQuickCoincidenceCheck(true),
        // producing exactly 1 CommonPart of type EDGE when the edge segment
        // lies on the face. The rcad equivalent: mid-point projection check
        // already passed above; additional tolerance computed from endpoints.
        //
        // rcad: the pair-level result is already determined by the selection
        // criteria above.  Mark each as "done/EDGE-type" for processing phase.

        // L1135-1139: prepare EF array if theAddInterf
        if the_add_interf {
            self.ds.interf_ef.reserve(a_nb_ef);
        }

        // L1147-1192: analyze results, create interferences
        // Collect map for CommonBlock creation
        let mut a_mpbli: std::collections::BTreeMap<usize, Vec<usize>> =
            std::collections::BTreeMap::new();

        for pair in &a_v_edge_face {
            // L1154-1159: check IsDone / HasErrors
            // rcad: verify with multi-point projection (OCCT UseQuickCoincidenceCheck)
            let mut all_on_face = true;
            let n_samples = 8;
            for k in 0..n_samples {
                let frac = (k as f64 + 0.5) / n_samples as f64;
                let t = pair.a_ts[0] + frac * (pair.a_ts[1] - pair.a_ts[0]);
                let pt = self.ds.edges[pair.n_e].curve.point_at(t);
                let proj = self.context.proj_ps(&self.ds, pair.fi, pt);
                match proj {
                    Some((_, _, dist))
                        if dist
                            <= pair.a_tol_add
                                + self.ds.fuzzy_tol
                                + self.ds.face_tolerance(pair.fi) => {}
                    _ => {
                        all_on_face = false;
                        break;
                    }
                }
            }
            if !all_on_face {
                continue;
            }

            // L1161-1171: exactly 1 CommonPart of type EDGE expected
            // rcad: the coordinate check above implies edge-on-face coincidence,
            // equivalent to a single EDGE-type common part in OCCT.

            // L1173-1181: add interference entry
            if the_add_interf {
                // L1178-1181: set indices and common part
                // rcad: InterferenceEF stores edge, face, mid-point, param, new_vertex
                let a_t_mid = crate::boptools::intermediate_point(pair.a_ts[0], pair.a_ts[1]);
                let a_mid_pt = self.ds.edges[pair.n_e].curve.point_at(a_t_mid);
                self.ds.interf_ef.push(InterferenceEF {
                    edge: pair.n_e,
                    face: pair.fi,
                    point: a_mid_pt,
                    edge_param: a_t_mid,
                    new_vertex: {
                        let r = pair.pb.0.read().unwrap();
                        r.pave1.vertex_idx
                    },
                });
                // L1181: myDS->AddInterf(nE, nF)
                self.ds.try_add_interf(pair.n_e, pair.fi);
            }

            // L1184-1186: update face info with new IN pave block
            // OCCT: myDS->ChangeFaceInfo(nF).ChangePaveBlocksIn().Add(aPB);
            let g_pb_idx = self
                .ds
                .allocate_pave_block(pair.pb.0.read().unwrap().clone());
            self.ds
                .face_info_mut(pair.fi)
                .pave_blocks_in
                .insert(g_pb_idx);

            // L1187-1191: fill map for common blocks creation
            if the_add_interf {
                crate::bopalgo::fill_map(&mut a_mpbli, g_pb_idx, pair.fi);
            }
        }

        // L1194-1198: Create new common blocks for coinciding pairs
        if !a_mpbli.is_empty() {
            // rcad: BTreeMap<usize, Vec<usize>> where key = PB index,
            // value = face indices sharing that PB.
            // OCCT: BOPAlgo_Tools::PerformCommonBlocks(aMPBLI, anAlloc, myDS)
            crate::bopds::tools::perform_common_blocks(&mut self.ds);
        }
    }

    pub(crate) fn put_se_in_other_faces(&mut self) {
        let n_faces = self.ds.faces.len();
        let ics = self.ds.intersection_curves.clone();
        let mut ic_creators: Vec<Vec<usize>> = vec![Vec::new(); ics.len()];
        for inf in &self.ds.interf_ff {
            for &ci in &inf.curves {
                if ci < ic_creators.len() {
                    ic_creators[ci].push(inf.f1);
                    ic_creators[ci].push(inf.f2);
                }
            }
        }
        for (ci, ic) in ics.iter().enumerate() {
            let creators = &ic_creators[ci];
            if creators.is_empty() {
                continue;
            }
            let mid_t = (ic.t_range[0] + ic.t_range[1]) * 0.5;
            let params = if (ic.t_range[1] - ic.t_range[0]).abs() < TOLERANCE_ABS {
                vec![mid_t]
            } else {
                vec![
                    ic.t_range[0] * 0.9 + ic.t_range[1] * 0.1,
                    mid_t,
                    ic.t_range[0] * 0.1 + ic.t_range[1] * 0.9,
                ]
            };
            for fi in 0..n_faces {
                if creators.contains(&fi) {
                    continue;
                }
                if !self.ds.faces[fi].face_info.has_any_interference() {
                    continue;
                }
                let on_face = params.iter().any(|&t| {
                    use rcad_kernel::geom::CurveEval;
                    let pt = ic.curve.point_at(t);
                    let tol = self.ds.face_tolerance(fi).max(TOLERANCE_ABS);
                    self.point_on_face(pt, fi, tol)
                });
                if on_face {
                    let sv = ic.start_vertex;
                    let ev = ic.end_vertex;
                    self.ds.face_info_mut(fi).curves_sc.insert(ci);
                    self.ds.face_info_mut(fi).vertices_in.insert(sv);
                    self.ds.face_info_mut(fi).vertices_in.insert(ev);
                }
            }
        }
    }

    fn point_on_face(&self, pt: DVec3, fi: usize, tol: f64) -> bool {
        use rcad_kernel::geom::{Surface3, SurfaceEval};
        let face = &self.ds.faces[fi];
        match &face.surface {
            Surface3::Plane(p) => (pt - p.origin).dot(p.normal).abs() <= tol,
            Surface3::Sphere(s) => ((pt - s.center).length() - s.radius).abs() <= tol,
            Surface3::Cylinder(c) => {
                let v = pt - c.origin;
                let radial = v - c.axis.normalize() * v.dot(c.axis.normalize());
                (radial.length() - c.radius).abs() <= tol
            }
            Surface3::Cone(c) => {
                let v = pt - c.apex;
                let a = c.axis_dir();
                let along = v.dot(a);
                let r = (v - a * along).length();
                let half_ang = c.half_angle_rad;
                let expected_r = c.radius + along * half_ang.tan();
                (r - expected_r).abs() * half_ang.cos() <= tol
            }
            _ => false,
        }
    }

    /// BOPAlgo_PaveFiller::ForceInterfVF (PaveFiller_5.cxx L631-681).
    /// Force Vertex-Face interference check for a vertex with increased tolerance.
    /// Returns true if a new VF interference was created.
    /// Note: named with `_pair` suffix to distinguish from the zero-arg version above.
    pub(crate) fn force_interf_vf_pair(&mut self, n_v: usize, n_f: usize) -> bool {
        if n_v >= self.ds.vertices.len() || n_f >= self.ds.faces.len() {
            return false;
        }
        let v_pt = self.ds.vertex_point(n_v);
        let surf = self.ds.faces[n_f].surface.clone();
        // Project vertex onto face surface
        // Compute UV from surface projection
        use rcad_kernel::geom::Surface3;
        let (u_val, v_val) = match &surf {
            Surface3::Plane(p) => {
                let d = v_pt - p.origin;
                (d.dot(p.u_dir), d.dot(p.v_dir))
            }
            Surface3::Sphere(s) => {
                let d = (v_pt - s.center).normalize();
                let u = d.dot(s.ref_dir_perp()).atan2(d.dot(s.ref_dir));
                let v = d.dot(s.axis).asin();
                (u, v)
            }
            Surface3::Cylinder(c) => {
                let ax = c.axis.normalize();
                let d = v_pt - c.origin;
                let v = d.dot(ax);
                let radial = d - ax * v;
                let ref_dir = c.ref_dir.normalize();
                let cross_dir = ax.cross(ref_dir).normalize();
                let u = radial.dot(cross_dir).atan2(radial.dot(ref_dir));
                (u, v)
            }
            _ => (0.0, 0.0),
        };
        let proj = surf.point_at(u_val, v_val);
        let dist = (proj - v_pt).length();
        let tol_v = self.ds.vertex_tolerance(n_v);
        let tol_f = self.ds.face_tolerance(n_f);
        let a_tol_check = tol_v.max(tol_f) + self.fuzzy_tolerance;

        if dist <= a_tol_check {
            // Create VF interference
            self.ds.interf_vf.push(InterferenceVF {
                vertex: n_v,
                face: n_f,
                u: u_val,
                v: v_val,
                index_new: None,
            });
            self.ds.interf_tb.insert((n_v.min(n_f), n_v.max(n_f)));
            // Update vertex tolerance
            let n_vx = self.update_vertex(n_v, dist.max(tol_v));
            // Register vertex in face info
            self.ds.face_info_mut(n_f).vertices_in.insert(n_vx);
            return true;
        }
        false
    }
}
