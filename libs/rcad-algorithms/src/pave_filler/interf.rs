use super::*;

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
        // rcad: PBs already initialized in earlier pipeline steps; skip.

        // L1024-1080: Fill the connection map from bounding vertices to PBs
        // OCCT: iterates over NbSourceShapes, filters to EDGE type.
        // rcad: source edges are 0..a_edge_count; skip intersection-created edges.
        let a_nb_s = self.ds.a_edge_count;
        // rcad: HashMap keyed by (v_min, v_max), value = Vec<(edge_idx, local_pb_idx)>
        // OCCT: NCollection_IndexedDataMap<BOPDS_Pair, NCollection_List<handle<PaveBlock>>>
        let mut a_pb_map: std::collections::HashMap<(usize, usize), Vec<(usize, usize)>> =
            std::collections::HashMap::new();

        for i in 0..a_nb_s {
            // L1034-1038: only edges (rcad: edges vec contains only edges)

            // L1040-1044: edge must have PBs (HasReference equivalent)
            if self.ds.edges[i].pave_blocks.is_empty() {
                // L1042-1044: No pave blocks
                continue;
            }

            // L1046-1050: skip degenerated edges
            if self.ds.is_edge_degenerated(i) {
                continue;
            }

            // L1056-1079: iterate PBs of this edge
            let a_lpb = &self.ds.edges[i].pave_blocks;
            for local_i in 0..a_lpb.len() {
                let a_pb = &a_lpb[local_i];

                // L1060-1065: get real PaveBlock + fence
                // rcad: PBs are unique per (edge, local). No RealPaveBlock indirection.
                // rcad: PBs have no CommonBlock indirection - skip fence check.

                // L1067-1069: get vertex indices
                let (n_v1, n_v2) = {
                    let pbr = a_pb.0.read().unwrap();
                    (pbr.pave1.vertex_idx, pbr.pave2.vertex_idx)
                };

                // L1071-1078: add PB to map keyed by vertex pair
                let a_pair = if n_v1 <= n_v2 { (n_v1, n_v2) } else { (n_v2, n_v1) };
                a_pb_map.entry(a_pair).or_default().push((i, local_i));
            }
        }

        // L1082-1086: empty map check
        if a_pb_map.is_empty() {
            return;
        }

        // L1088: Self-Interference check mode (single argument)
        let b_si_check_mode = self.my_arguments.len() <= 1;

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
                2.0 * self.ds.vertex_tolerance(n_v1).max(self.ds.vertex_tolerance(n_v2))
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

        // L1233-1238: close preparation step (rcad: no allocator to cleanup)

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
                2.0 * self.ds.vertex_tolerance(entry.n_v1)
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
                    if d_at_mid > fuzzy { continue; }
                    let v_tgt1 = curve1.tangent_at(mid_t1);
                    if v_tgt1.length_squared() < 1e-60 { continue; }
                    let v_tgt2 = curve2.tangent_at(proj_mid.param);
                    if v_tgt2.length_squared() < 1e-60 { continue; }
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
                    Some((InterferenceEE {
                        e1: entry.n_e1,
                        e2: entry.n_e2,
                        point: avg_pt,
                        param1: mid_t1,
                        param2: mid_t2,
                        new_vertex: entry.n_v1,
                    }, true))
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
                self.my_report.add_alert(crate::bopalgo::Alert::AcquiredSelfIntersection(
                    vec![entry.n_e1, entry.n_e2],
                ));
            }

            // L1307-1310: create InterfEE entry
            self.ds.interf_ee.push(ee.clone());
            self.ds.interf_tb.insert(
                (entry.n_e1.min(entry.n_e2), entry.n_e1.max(entry.n_e2)),
            );

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

        // L1332: Create new common blocks of coinciding pairs
        if a_mpblpb.len() > 0 {
            crate::bopds::tools::perform_common_blocks(self.ds);
        }
    }

    // OCCT BOPAlgo_PaveFiller_5.cxx L772-827
    pub(crate) fn force_interf_ef(&mut self) {
        // L774-778
        if !self.is_primary {
            return;
        }

        // L787-822: Collect all pave blocks from all source edges
        // rcad: Vec<(edge_idx, local_pb_idx)> — OCCT uses IndexedMap<handle<PaveBlock>>
        let mut a_mpb: Vec<(usize, usize)> = Vec::new();
        let a_nb_s = self.ds.nb_source_shapes();
        for n_e in 0..a_nb_s {
            // L791-796: only edges
            if self.ds.shape_type_of(n_e) != rcad_kernel::topods::ShapeType::Edge {
                continue;
            }
            // rcad: shape_info[n_e].reference gives the flat edge index
            let ne = self.ds.shape_info[n_e].reference as usize;

            // L798-802: edge must have PBs (HasReference)
            if !self.ds.has_pave_blocks(ne) {
                continue;
            }

            // L804-808: skip degenerated edges
            if self.ds.is_edge_degenerated(ne) {
                continue;
            }

            // L814-821: add all PBs to map
            // OCCT uses RealPaveBlock(aPB) to resolve the real edge's PB
            let a_lpb = self.ds.edge_pave_blocks(ne);
            for (local_i, _a_it_lpb) in a_lpb.iter().enumerate() {
                // rcad: real_pave_block_edge resolves through CommonBlock to
                // the original edge index. OCCT uses RealPaveBlock to get the
                // handle of the real pave block for IndexedMap dedup.
                // In rcad, we store the index pair and resolve per-pair later.
                a_mpb.push((ne, local_i));
            }
        }

        // L826: call overload 2 with theAddInterf=true
        self.force_interf_ef_with(&a_mpb, true);
    }

    // OCCT BOPAlgo_PaveFiller_5.cxx L831-1199
    fn force_interf_ef_with(
        &mut self,
        the_mpb: &[(usize, usize)],
        the_add_interf: bool,
    ) {
        // L838-840
        if the_mpb.is_empty() {
            return;
        }

        // L842-874: Fill the tree with bounding boxes of the pave blocks
        // rcad: build DsBvh from PB AABBs (shrunk data or edge endpoint fallback).
        let mut a_pb_aabbs: Vec<crate::bvh::Aabb> = Vec::with_capacity(the_mpb.len());
        let mut a_pb_indices: Vec<usize> = Vec::with_capacity(the_mpb.len());
        // L848-870: for each PB, get ShrunkData and add to tree
        for (i_pb, &(ne, local_i)) in the_mpb.iter().enumerate() {
            let a_pb = &self.ds.edges[ne].pave_blocks[local_i];
            // L852-858: ensure shrunk data is available
            // rcad: PB.shrunk_range holds the shrunk data (f, l) and splits
            let a_pb_r = a_pb.0.read().unwrap();
            let (shrunk, _is_split) = match a_pb_r.shrunk_range {
                Some(sr) => (sr, a_pb_r.is_splittable),
                None => continue,
            };
            drop(a_pb_r);
            // L866-870: build AABB from shrunk range endpoint positions
            let a_p1 = self.ds.edges[ne].curve.point_at(shrunk[0]);
            let a_p2 = self.ds.edges[ne].curve.point_at(shrunk[1]);
            let a_pb_box = crate::bvh::Aabb {
                min: a_p1.min(a_p2),
                max: a_p1.max(a_p2), gap: 0.0 };
            // L870: aBBTree.Add(aPBMap.Add(aPB), Bnd_Tools::Bnd2BVH(aPBBox));
            a_pb_aabbs.push(a_pb_box);
            a_pb_indices.push(i_pb);
        }
        // L873-874: Build BVH tree
        let a_bb_tree = DsBvh::build(a_pb_indices, a_pb_aabbs);

        // L876: bSICheckMode — Self-Interference check mode (one argument)
        // OCCT: myArguments.Extent() == 1
        let b_si_check_mode = self.my_arguments.len() <= 1;

        // L880: vector of edge-face pairs for intersection
        // OCCT: BOPAlgo_VectorOfEdgeFace aVEdgeFace
        // rcad: pair data stored in a local struct
        struct _EFPair {
            n_e: usize,
            n_f: usize,
            pb: SharedPB,
            a_tol_add: f64,
            a_ts: [f64; 2],
        }
        let mut a_v_edge_face: Vec<_EFPair> = Vec::new();

        // L882-1108: For each face that has face info
        // OCCT iterates over all source shapes with type FACE + HasReference check.
        // rcad: iterate over flat face array and map source idx for rank/origin checks.
        let a_nb_f = self.ds.faces.len();
        for n_f in 0..a_nb_f {
            let a_fi = &self.ds.faces[n_f].face_info;

            // L885-896: only faces with face info (HasReference equivalent)
            // rcad: check that face info has any data
            if a_fi.pave_blocks_on.is_empty()
                && a_fi.pave_blocks_in.is_empty()
                && a_fi.pave_blocks_sc.is_empty()
                && a_fi.curves_sc.is_empty()
                && a_fi.vertices_on.is_empty()
                && a_fi.vertices_in.is_empty()
                && a_fi.vertices_sc.is_empty()
            {
                continue;
            }

            // L903-910: Face AABB for BVH query
            // OCCT: BOPTools_BoxTreeSelector with face box from ShapeInfo
            let a_box_f = face_aabb::face_aabb(&self.ds, n_f);
            let a_overlapping = a_bb_tree.query_aabb(&a_box_f);
            if a_overlapping.is_empty() {
                continue;
            }

            // L914-924: build aMVF from all face vertex sets
            // OCCT: NCollection_Map<int> aMVF from FaceInfo VerticesOn/In/Sc
            let mut a_mvf: std::collections::HashSet<usize> =
                std::collections::HashSet::new();
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
                    let skip = a_fi.pave_blocks_on.iter()
                        .chain(a_fi.pave_blocks_in.iter())
                        .chain(a_fi.pave_blocks_sc.iter())
                        .any(|&pb_gi| {
                            pb_gi < self.ds.pave_blocks.len()
                                && self.ds.pave_blocks[pb_gi].0.read().unwrap().original_edge == a_pb_ne
                        });
                    if skip { continue; }
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
                if self.ds.edge_origin(ne) == self.ds.face_origin(n_f) {
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
                    match self.context.proj_ps(&self.ds, n_f, a_p_on_e) {
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
                if !self.context.is_point_in_face(&self.ds, n_f, a_uv) {
                    continue;
                }

                // L1038-1052: non-plane/non-line → angle check for bUseAddTol
                // OCCT: if (aSurfAdaptor.GetType() != GeomAbs_Plane || aBAC.GetType() != GeomAbs_Line)
                let a_face = &self.ds.faces[n_f];
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
                    let (_, a_p_on_s_re, _) =
                        match self.context.proj_ps(&self.ds, n_f, a_p_on_e) {
                            Some(data) => data,
                            None => continue,
                        };
                    // L1041: aVFNorm = aPOnS - aPOnE
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
                        if let Some((_, _, a_dist_ef)) =
                            self.context.proj_ps(&self.ds, n_f, a_p)
                        {
                            if a_dist_ef < a_tol_check && a_dist_ef > a_tol_add {
                                a_tol_add = a_dist_ef;
                            }
                        }
                    }
                    // L1077-1084: subtract edge + face tolerances
                    if a_tol_add > 0.0 {
                        a_tol_add -= (self.ds.edge_tolerance(ne)
                            + self.ds.face_tolerance(n_f));
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
                        n_f,
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
                let proj = self.context.proj_ps(&self.ds, pair.n_f, pt);
                match proj {
                    Some((_, _, dist)) if dist <= pair.a_tol_add + self.ds.fuzzy_tol + self.ds.face_tolerance(pair.n_f) => {}
                    _ => { all_on_face = false; break; }
                }
            }
            if !all_on_face { continue; }

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
                    face: pair.n_f,
                    point: a_mid_pt,
                    edge_param: a_t_mid,
                    new_vertex: {
                        let r = pair.pb.0.read().unwrap();
                        r.pave1.vertex_idx
                    },
                });
                // L1181: myDS->AddInterf(nE, nF)
                self.ds.try_add_interf(pair.n_e, pair.n_f);
            }

            // L1184-1186: update face info with new IN pave block
            // OCCT: myDS->ChangeFaceInfo(nF).ChangePaveBlocksIn().Add(aPB);
            let g_pb_idx = self.ds.allocate_pave_block(
                pair.pb.0.read().unwrap().clone()
            );
            self.ds.face_info_mut(pair.n_f).pave_blocks_in.insert(g_pb_idx);

            // L1187-1191: fill map for common blocks creation
            if the_add_interf {
                crate::bopalgo::fill_map(&mut a_mpbli, g_pb_idx, pair.n_f);
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
            for &ci in &inf.curves { if ci < ic_creators.len() { ic_creators[ci].push(inf.f1); ic_creators[ci].push(inf.f2); } }
        }
        for (ci, ic) in ics.iter().enumerate() {
            let creators = &ic_creators[ci];
            if creators.is_empty() { continue; }
            let mid_t = (ic.t_range[0] + ic.t_range[1]) * 0.5;
            let params = if (ic.t_range[1] - ic.t_range[0]).abs() < TOLERANCE_ABS { vec![mid_t] }
            else { vec![ic.t_range[0]*0.9+ic.t_range[1]*0.1, mid_t, ic.t_range[0]*0.1+ic.t_range[1]*0.9] };
            for fi in 0..n_faces {
                if creators.contains(&fi) { continue; }
                if !self.ds.faces[fi].face_info.has_any_interference() { continue; }
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
        use rcad_kernel::geom::{SurfaceEval, Surface3};
        let face = &self.ds.faces[fi];
        match &face.surface {
            Surface3::Plane(p) => (pt - p.origin).dot(p.normal).abs() <= tol,
            Surface3::Sphere(s) => ((pt - s.center).length() - s.radius).abs() <= tol,
            Surface3::Cylinder(c) => { let v = pt - c.origin; let radial = v - c.axis.normalize() * v.dot(c.axis.normalize()); (radial.length() - c.radius).abs() <= tol }
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
        if n_v >= self.ds.vertices.len() || n_f >= self.ds.faces.len() { return false; }
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
