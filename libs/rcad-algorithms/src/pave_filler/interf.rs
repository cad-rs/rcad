use super::*;

impl<'a> PaveFiller<'a> {
    pub(crate) fn force_interf_ee(&mut self) {
        // OCCT L999-1003: ForceInterfEE — find additional common blocks among
        // pairs of edges.  Since all real intersections have already happened,
        // here we are interested in common blocks only, thus we check only
        // pairs of pave blocks with the same bounding vertices.

        // OCCT L1008-1023: InitPaveBlocksForVertex for interfered vertices.
        // rcad: PBs already initialized in earlier pipeline steps; skip.

        // OCCT L1026-1080: build (v1,v2) → Vec<(edge, local_pb)> map from ALL edge PBs
        // OCCT L1029-1030: aMPBFence — in OCCT, a fence avoids double-counting the same
        //   PB (shared via CommonBlock). rcad PBs are per-edge and unique; skip fence.
        let mut pb_map: std::collections::HashMap<(usize, usize), Vec<(usize, usize)>> =
            std::collections::HashMap::new();

        for (ei, edge) in self.ds.edges.iter().enumerate() {
            // OCCT L1034-1039: ShapeType != TopAbs_EDGE → skip (all DS entries are edges)
            // OCCT L1041-1045: HasReference (non-empty pave_blocks)
            if edge.pave_blocks.is_empty() { continue; }
            // OCCT L1047-1051: HasFlag → skip degenerated edges
            if self.ds.is_edge_degenerated(ei) { continue; }

            for (local_pbi, pb) in edge.pave_blocks.iter().enumerate() {
                // OCCT L1060-1065: RealPaveBlock + fence dedup (skip — see above)
                // OCCT L1068-1069: aPBR->Indices(nV1, nV2)
                let (nV1, nV2) = pb.0.write().unwrap().indices();
                let key = if nV1 < nV2 { (nV1, nV2) } else { (nV2, nV1) };
                // OCCT L1073-1078: append to map
                pb_map.entry(key).or_default().push((ei, local_pbi));
            }
        }

        // OCCT L1082-1086: early return if no entries
        if pb_map.is_empty() { return; }

        // OCCT L1090: aVEdgeEdge — instead of a parallel vector, process pairs inline
        // OCCT L1093-1225: iterate map entries with ≥2 PBs
        for (&(nV1, nV2), pbs) in &pb_map {
            // OCCT L1100-1103: skip entries with < 2 PBs
            if pbs.len() < 2 { continue; }

            // OCCT L1105-1118: vertex tolerances for aTolAdd
            // OCCT L1116: aTolAdd = 2 * max(Tol(V1), Tol(V2))
            let a_tol_add = 2.0 * self.ds.vertices[nV1].geom_tol
                                .max(self.ds.vertices[nV2].geom_tol);

            // OCCT L1120-1223: iterate all (pb1,pb2) pairs from the list
            for i in 0..pbs.len() {
                let (ei1, pb1_local) = pbs[i];
                let pb1 = &self.ds.edges[ei1].pave_blocks[pb1_local];
                let nE1 = pb1.0.read().unwrap().original_edge;                  // OCCT L1127
                let r1 = self.ds.edges[nE1].origin;
                let (t11, t12) = pb1.0.write().unwrap().range();                 // OCCT L1130
                let mid_t1 = (t11 + t12) * 0.5;              // OCCT L1134
                let c1 = &self.ds.edges[nE1].curve;           // OCCT L1131: BRepAdaptor_Curve
                // OCCT L1132-1134: tangent at midpoint
                let tgt1 = c1.tangent_at(mid_t1);
                // OCCT L1135-1139: skip if |tangent| < Resolution
                if tgt1.length_squared() < 1e-60 { continue; }

                for j in (i + 1)..pbs.len() {
                    let (ei2, pb2_local) = pbs[j];
                    let pb2 = &self.ds.edges[ei2].pave_blocks[pb2_local];
                    let nE2 = pb2.0.read().unwrap().original_edge;              // OCCT L1145
                    let r2 = self.ds.edges[nE2].origin;       // OCCT L1147: iR2
                    let (t21, t22) = pb2.0.write().unwrap().range();             // OCCT L1173

                    // OCCT L1149-1160: skip edges from the same argument
                    //   if the bounding vertices are original (not acquired during operation)
                    if r1 == r2 {
                        let o_rank = if r1 == ShapeOrigin::ShapeA { 0 } else { 1 };
                        // OCCT L1155-1158: IsNewShape + Rank check
                        if (!self.ds.is_new_vertex(nV1) && self.ds.rank(nV1) == o_rank)
                            || (!self.ds.is_new_vertex(nV2) && self.ds.rank(nV2) == o_rank)
                        {
                            continue;
                        }
                    }

                    // OCCT L1162-1169: skip if PBs already form a CommonBlock
                    // rcad: has_interf_ee checks for existing EE interference
                    if self.ds.has_interf_ee(nE1, nE2) { continue; }

                    // OCCT L1178: use angle between edges to decide bUseAddTol
                    let c2 = &self.ds.edges[nE2].curve;       // OCCT L1180: BRepAdaptor_Curve aBAC2
                    let tgt1_n = tgt1.normalize();            // OCCT L1139
                    let b_use_add_tol = match (c1, c2) {
                        // OCCT L1181: both lines → direction dot check
                        (Curve3::Line(l1), Curve3::Line(l2)) => {
                            let cos_angle = l1.direction.dot(l2.direction).abs();
                            cos_angle >= 0.9063  // OCCT L1200: 25°
                        }
                        _ => {
                            // OCCT L1183-1191: project midpoint of curve1 onto curve2
                            let mid_pt = c1.point_at(mid_t1);
                            // OCCT L1185: GeomAPI_ProjectPointOnCurve aProjPC(aE2);
                            let proj = closest_point_on_curve(c2, mid_pt, 64);
                            if !proj.param.is_finite() { false } else {
                                // OCCT L1192: aBAC2.D1(projPC.LowerDistanceParameter(), ...)
                                let tgt2 = c2.tangent_at(proj.param);
                                if tgt2.length_squared() < 1e-60 { false } else {
                                    let cos_angle = tgt1_n.dot(tgt2.normalize()).abs();
                                    cos_angle >= 0.9063  // OCCT L1200
                                }
                            }
                        }
                    };

                    // OCCT L1207-1222: create BOPAlgo_EdgeEdge
                    let fuzzy = if b_use_add_tol {
                        self.ds.fuzzy_tol + a_tol_add  // OCCT L1217: myFuzzyValue + aTolAdd
                    } else {
                        self.ds.fuzzy_tol              // OCCT L1221
                    };

                    // ── BEGIN existing intersection code (preserved verbatim) ──
                    match (c1, c2) {
                        (Curve3::Line(l1), Curve3::Line(l2)) => {
                            // OCCT L1097-1102: midpoint direction vector, check angle
                            // intersect_line_line returns Option<(f64,f64,DVec3)>
                            if let Some((t1, t2, pt)) = intersect_line_line(
                                l1, self.ds.edges[nE1].t_range,
                                l2, self.ds.edges[nE2].t_range, fuzzy)
                            {
                                self.ds.interf_ee.push(InterferenceEE{
                                    e1: nE1, e2: nE2, point: pt, param1: t1, param2: t2, new_vertex: nV1,
                                });
                            }
                        }
                        (Curve3::Circle(circ), Curve3::Circle(_)) => {
                            // intersect_circle_circle returns Vec<DVec3>
                            let cp_hits = intersect_circle_circle(circ, circ, fuzzy);
                            if let Some(&pt) = cp_hits.first() {
                                self.ds.interf_ee.push(InterferenceEE{
                                    e1: nE1, e2: nE2, point: pt, param1: 0.0, param2: 0.0, new_vertex: nV1,
                                });
                            }
                        }
                        _ => {
                            // OCCT IntTools_EdgeEdge: coarse + adaptive + Newton
                            // (1) Coarse 21x21 grid -> find best (t1,t2)
                            // (2) Recursive subdivision around best: 2x denser per level
                            // (3) Converge when distance < fuzzy OR subrange < 1e-6
                            let tr1 = self.ds.edges[nE1].t_range;
                            let tr2 = self.ds.edges[nE2].t_range;
                            let mid_t1 = (tr1[0] + tr1[1]) * 0.5;
                            let mid_t2 = (tr2[0] + tr2[1]) * 0.5;
                            let tgt1 = c1.tangent_at(mid_t1);
                            let tgt2 = c2.tangent_at(mid_t2);
                            let cos_angle = if tgt1.length_squared() > 1e-30 && tgt2.length_squared() > 1e-30 {
                                tgt1.normalize().dot(tgt2.normalize()).abs()
                            } else { 0.0 };
                            let fuzzy = if cos_angle >= 0.9063 {
                                self.ds.fuzzy_tol + a_tol_add
                            } else {
                                self.ds.fuzzy_tol
                            };
                            let mut best_t1 = mid_t1;
                            let mut best_t2 = mid_t2;
                            let mut best_d = f64::MAX;
                            // OCCT N=20 -> 21 samples per curve
                            for si in 0..21 {
                                let t1 = tr1[0] + (tr1[1] - tr1[0]) * (si as f64 / 20.0);
                                let p1 = c1.point_at(t1);
                                for sj in 0..21 {
                                    let t2 = tr2[0] + (tr2[1] - tr2[0]) * (sj as f64 / 20.0);
                                    let d = p1.distance(c2.point_at(t2));
                                    if d < best_d { best_d = d; best_t1 = t1; best_t2 = t2; }
                                }
                            }
                            // (2) Adaptive refinement: subdivide around min point
                            let mut r1_lo = (best_t1 - (tr1[1] - tr1[0]) / 20.0).max(tr1[0]);
                            let mut r1_hi = (best_t1 + (tr1[1] - tr1[0]) / 20.0).min(tr1[1]);
                            let mut r2_lo = (best_t2 - (tr2[1] - tr2[0]) / 20.0).max(tr2[0]);
                            let mut r2_hi = (best_t2 + (tr2[1] - tr2[0]) / 20.0).min(tr2[1]);
                            for _ in 0..4 {
                                let mid1 = (r1_lo + r1_hi) * 0.5;
                                let mid2 = (r2_lo + r2_hi) * 0.5;
                                let test_t1 = [r1_lo, mid1, r1_hi];
                                let test_t2 = [r2_lo, mid2, r2_hi];
                                for &t1 in &test_t1 {
                                    let pt1 = c1.point_at(t1);
                                    for &t2 in &test_t2 {
                                        let d = pt1.distance(c2.point_at(t2));
                                        if d < best_d { best_d = d; best_t1 = t1; best_t2 = t2; }
                                    }
                                }
                                let span = (r1_hi - r1_lo) * 0.5;
                                r1_lo = (best_t1 - span).max(tr1[0]);
                                r1_hi = (best_t1 + span).min(tr1[1]);
                                r2_lo = (best_t2 - span).max(tr2[0]);
                                r2_hi = (best_t2 + span).min(tr2[1]);
                            }
                            // (3) OCCT IntTools_CurveRange L230-260: Newton-Raphson iteration
                            // Minimize F(t1,t2) = ||C1(t1)-C2(t2)||^2 using gradient+Hessian.
                            let mut nr_t1 = best_t1;
                            let mut nr_t2 = best_t2;
                            for _ in 0..8 {
                                let p1 = c1.point_at(nr_t1);
                                let p2 = c2.point_at(nr_t2);
                                let diff = p1 - p2;
                                if diff.length_squared() < 1e-30 { break; }
                                let t1 = c1.tangent_at(nr_t1);
                                let t2 = c2.tangent_at(nr_t2);
                                if t1.length_squared() < 1e-30 || t2.length_squared() < 1e-30 { break; }
                                let d1 = t1.normalize();
                                let d2 = t2.normalize();
                                // Hessian H and gradient grad of F(t1,t2) = ||C1-C2||^2
                                let h00 = 2.0;
                                let h01 = -2.0 * d1.dot(d2);
                                let h10 = h01;
                                let h11 = 2.0;
                                let g0 = 2.0 * diff.dot(d1);
                                let g1 = 2.0 * diff.dot(d2);
                                let det = h00 * h11 - h01 * h01;
                                if det.abs() < 1e-30 { break; }
                                let dt1 = (-g0 * h11 - g1 * h01) / det;
                                let dt2 = (g1 * h00 + g0 * h10) / det;
                                let new_t1 = (nr_t1 + dt1).clamp(tr1[0], tr1[1]);
                                let new_t2 = (nr_t2 + dt2).clamp(tr2[0], tr2[1]);
                                if (new_t1 - nr_t1).abs() < 1e-12 && (new_t2 - nr_t2).abs() < 1e-12 { break; }
                                nr_t1 = new_t1; nr_t2 = new_t2;
                            }
                            let nr_d = c1.point_at(nr_t1).distance(c2.point_at(nr_t2));
                            if nr_d < best_d { best_d = nr_d; best_t1 = nr_t1; best_t2 = nr_t2; }
                            if best_d <= fuzzy {
                                let best_pt = c1.point_at(best_t1);
                                self.ds.interf_ee.push(InterferenceEE{
                                    e1: nE1, e2: nE2, point: best_pt,
                                    param1: best_t1, param2: best_t2, new_vertex: nV1,
                                });
                            }
                        }
                    }
                    // ── END existing intersection code ──
                }
            }
        }
        // OCCT L1332: PerformCommonBlocks — in rcad, the EE interferences above
        // are consumed downstream by the builder.
    }

    /// OCCT L772-827: ForceInterfEF (overload 1) — collect all edge PBs
    pub(crate) fn force_interf_ef(&mut self) {
        // L774-778
        if !self.is_primary {
            return;
        }

        // L787-822: collect all PBs from all edges (skip no-PB, degenerated)
        // rcad: Vec<(edge_idx, local_pb_idx)> — OCCT uses IndexedMap<handle<PaveBlock>>
        let mut a_mpb: Vec<(usize, usize)> = Vec::new();
        let a_nb_s = self.ds.edges.len();
        for ne in 0..a_nb_s {
            // L791-796: only edges (rcad: edges Vec contains only edges)

            // L798-802: edge must have PBs
            if self.ds.edges[ne].pave_blocks.is_empty() {
                continue;
            }

            // L804-808: skip degenerated edges
            if self.ds.is_edge_degenerated(ne) {
                continue;
            }

            // L814-821: add all PBs to map
            for local_i in 0..self.ds.edges[ne].pave_blocks.len() {
                a_mpb.push((ne, local_i));
            }
        }

        // L826: call overload 2
        self.force_interf_ef_with(&a_mpb);
    }

    /// OCCT L831-1199: ForceInterfEF (overload 2, with theMPB)
    fn force_interf_ef_with(&mut self, the_mpb: &[(usize, usize)]) {
        // L838-840
        if the_mpb.is_empty() {
            return;
        }

        // L842-874: Build BVH tree of PB shrunk-data boxes
        // OCCT: BOPTools_BoxTree over PB AABBs for spatial pair filtering.
        // rcad: build DsBvh from PB shrunk-range AABBs (or edge endpoint fallback).
        let mut pb_aabbs: Vec<crate::bvh::Aabb> = Vec::with_capacity(the_mpb.len());
        let mut pb_indices: Vec<usize> = Vec::with_capacity(the_mpb.len());
        for (pi, &(ne, local_i)) in the_mpb.iter().enumerate() {
            if ne >= self.ds.edges.len() { continue; }
            let pb = &self.ds.edges[ne].pave_blocks[local_i];
            let aabb = if let Some((mn, mx)) = pb.0.read().unwrap().my_shrunk_box {
                crate::bvh::Aabb { min: mn, max: mx }
            } else {
                // Fallback: AABB from edge endpoint vertices
                let v1 = if pb.0.read().unwrap().pave1.vertex_idx < self.ds.vertices.len() {
                    self.ds.vertices[pb.0.read().unwrap().pave1.vertex_idx].point
                } else { continue; };
                let v2 = if pb.0.read().unwrap().pave2.vertex_idx < self.ds.vertices.len() {
                    self.ds.vertices[pb.0.read().unwrap().pave2.vertex_idx].point
                } else { continue; };
                crate::bvh::Aabb { min: v1.min(v2), max: v1.max(v2) }
            };
            pb_aabbs.push(aabb);
            pb_indices.push(pi);
        }
        let pb_tree = crate::bvh::DsBvh::build(pb_indices, pb_aabbs);

        // L876: bSICheckMode — Self-Interference check mode (one argument)
        let b_si_check_mode = self.my_arguments.len() <= 1;

        // L880: EdgeFace pairs for intersection
        // rcad: Vec<(edge_idx, face_idx, local_pb_idx)>
        let mut a_v_edge_face: Vec<(usize, usize, usize)> = Vec::new();

        // L882-1108: For each face with face info
        // G3: parallel face processing with per-thread IntToolsContext
        // (OCCT BOPAlgo_Parallel equivalent — each thread has its own context).
        use rayon::prelude::*;
        let n_faces = self.ds.faces.len();
        let ds = &self.ds;
        let pb_tree = &pb_tree;
        let fpbdone = &self.fpbdone;
        let context_tol = self.context.tol_uv();
        let b_si_check_mode = b_si_check_mode;

        let per_face_pairs: Vec<Vec<(usize, usize, usize)>> = (0..n_faces)
            .into_iter()
            .map(|nf| {
                let mut ctx = crate::inttools::context::Context::new(n_faces, context_tol);
                let fi = &ds.faces[nf].face_info;

                // L885-896: skip faces with no face info
                if fi.pave_blocks_on.is_empty()
                    && fi.pave_blocks_in.is_empty()
                    && fi.pave_blocks_sc.is_empty()
                    && fi.curves_sc.is_empty()
                    && fi.vertices_on.is_empty()
                    && fi.vertices_in.is_empty()
                    && fi.vertices_sc.is_empty()
                {
                    return Vec::new();
                }

                // L914-924: build aMVF from all face vertex sets
                let mut a_mvf: std::collections::HashSet<usize> = std::collections::HashSet::new();
                for &v in &fi.vertices_on { a_mvf.insert(v); }
                for &v in &fi.vertices_in { a_mvf.insert(v); }
                for &v in &fi.vertices_sc { a_mvf.insert(v); }

                // L926-938: add PB endpoints from face's PaveBlocksOn/In/Sc
                for &pb_gi in &fi.pave_blocks_on {
                    if pb_gi < ds.pave_blocks.len() {
                        a_mvf.insert(ds.pave_blocks[pb_gi].0.read().unwrap().pave1.vertex_idx);
                        a_mvf.insert(ds.pave_blocks[pb_gi].0.read().unwrap().pave2.vertex_idx);
                    }
                }
                for &pb_gi in &fi.pave_blocks_in {
                    if pb_gi < ds.pave_blocks.len() {
                        a_mvf.insert(ds.pave_blocks[pb_gi].0.read().unwrap().pave1.vertex_idx);
                        a_mvf.insert(ds.pave_blocks[pb_gi].0.read().unwrap().pave2.vertex_idx);
                    }
                }
                for &pb_gi in &fi.pave_blocks_sc {
                    if pb_gi < ds.pave_blocks.len() {
                        a_mvf.insert(ds.pave_blocks[pb_gi].0.read().unwrap().pave1.vertex_idx);
                        a_mvf.insert(ds.pave_blocks[pb_gi].0.read().unwrap().pave2.vertex_idx);
                    }
                }

                // L947-949: compute face AABB
                let a_face = &ds.faces[nf];
                let mut face_aabb = crate::bvh::Aabb::empty();
                for &vi in &a_face.boundary_verts {
                    if vi < ds.vertices.len() {
                        face_aabb.expand_point(ds.vertices[vi].point);
                    }
                }

                // L950-952: query BVH for PBs overlapping this face
                let overlapping = pb_tree.query_aabb(&face_aabb);

                // L954-1107: iterate overlapping PBs
                let mut local_pairs: Vec<(usize, usize, usize)> = Vec::new();
                for &pb_idx in &overlapping {
                    let &(ne, local_i) = &the_mpb[pb_idx];
                    let pb = &ds.edges[ne].pave_blocks[local_i];

                    // skip if PB already in face's EF set
                    if ds.interf_ef.iter().any(|inf| inf.edge == ne && inf.face == nf) {
                        continue;
                    }

                    let n_v1 = pb.0.read().unwrap().pave1.vertex_idx;
                    let n_v2 = pb.0.read().unwrap().pave2.vertex_idx;
                    if !a_mvf.contains(&n_v1) || !a_mvf.contains(&n_v2) { continue; }

                    if ds.edges[ne].origin == ds.faces[nf].origin { continue; }

                    let a_e_curve = &ds.edges[ne].curve;
                    let mut b_use_add_tol = true;
                    let a_ts = pb.0.read().unwrap().shrunk_range.unwrap_or([pb.0.read().unwrap().pave1.param, pb.0.read().unwrap().pave2.param]);
                    let a_t_mid = crate::boptools::intermediate_point(a_ts[0], a_ts[1]);
                    let a_p_on_e = a_e_curve.point_at(a_t_mid);
                    let a_ve_tgt = a_e_curve.tangent_at(a_t_mid);
                    if a_ve_tgt.length_squared() < 1e-30 { continue; }

                    let (a_uv, a_p_on_s, a_lower_dist) =
                        match ctx.proj_ps(ds, nf, a_p_on_e) {
                            Some(data) => data,
                            None => continue,
                        };

                    let a_tol_v1 = ds.vertices.get(n_v1).map(|v| v.geom_tol).unwrap_or(0.0);
                    let a_tol_v2 = ds.vertices.get(n_v2).map(|v| v.geom_tol).unwrap_or(0.0);
                    let a_tol_check = if b_si_check_mode {
                        ds.fuzzy_tol
                    } else {
                        2.0 * a_tol_v1.max(a_tol_v2)
                    };
                    if a_lower_dist > a_tol_check + ds.fuzzy_tol { continue; }

                    if !ctx.is_point_in_face(ds, nf, a_uv) { continue; }

                    if !matches!(a_face.surface, Surface3::Plane(_))
                        || !matches!(*a_e_curve, Curve3::Line(_))
                    {
                        let a_vf_norm = a_p_on_s - a_p_on_e;
                        if a_vf_norm.length_squared() > 1e-30 {
                            let a_cos = a_vf_norm.normalize().dot(a_ve_tgt.normalize());
                            if a_cos.abs() > 0.4226 { b_use_add_tol = false; }
                        }
                    }

                    let mut a_tol_add = 0.0;
                    if b_use_add_tol {
                        for &a_t in &[a_ts[0], a_ts[1]] {
                            let a_p = a_e_curve.point_at(a_t);
                            if let Some((_, _, a_dist_ef)) = ctx.proj_ps(ds, nf, a_p) {
                                if a_dist_ef < a_tol_check && a_dist_ef > a_tol_add {
                                    a_tol_add = a_dist_ef;
                                }
                            }
                        }
                        if a_tol_add > 0.0 {
                            a_tol_add -= (ds.edges[ne].geom_tol + a_face.geom_tol);
                            if a_tol_add < 0.0 { a_tol_add = 0.0; }
                        }
                    }

                    let b_intersect = a_tol_add > 0.0
                        || !fpbdone.get(&nf).map_or(false, |done| {
                            done.contains(&(ne << 16 | local_i))
                        });
                    if b_intersect {
                        local_pairs.push((ne, nf, local_i));
                    }
                }
                local_pairs
            })
            .collect();

        // Merge per-face results into a_v_edge_face
        for pairs in &per_face_pairs {
            a_v_edge_face.extend(pairs);
        }

        // L1110-1113: no pairs
        let a_nb_ef = a_v_edge_face.len();
        if a_nb_ef == 0 {
            return;
        }

        // L1129: process EF pairs sequentially.
        // OCCT uses BOPAlgo_EdgeFace parallel solver with thread-local contexts.
        // rcad: sequential — IntTools_Context (proj_ps) requires &mut self,
        // preventing parallel closure capture.  Thread-local contexts could
        // enable parallel execution in a future refactor.

        // L1147-1192: process results
        // rcad: create Interference::EdgeFace + update FaceInfo
        // Collect map for CommonBlock creation: PB index -> face list
        // rcad: use BTreeMap for ordered iteration
        let mut a_mpbli: std::collections::BTreeMap<usize, Vec<usize>> =
            std::collections::BTreeMap::new();

        for &(ne, nf, local_i) in &a_v_edge_face {
            let pb = &self.ds.edges[ne].pave_blocks[local_i];

            // L1154-1158: (rcad: skip error check — no intersection engine)
            // L1162-1170: results check (rcad: all pairs produce a result)

            // L1174-1183: Add EdgeFace interference (rcad equivalent)
            let a_t_mid = crate::boptools::intermediate_point(pb.0.read().unwrap().pave1.param, pb.0.read().unwrap().pave2.param);
            let a_mid_pt = self.ds.edges[ne].curve.point_at(a_t_mid);

            // OCCT L1176-1182: create BOPDS_InterfEF + AddInterf
            self.ds.interf_ef.push(InterferenceEF{
                edge: ne,
                face: nf,
                point: a_mid_pt,
                edge_param: a_t_mid,
                // Use first vertex index as new_vertex placeholder
                new_vertex: pb.0.read().unwrap().pave1.vertex_idx,
            });

            // L1184-1186: Update face info with new IN pave block
            let g_pb_idx = self.ds.pave_blocks.len();
            self.ds.pave_blocks.push(pb.clone());
            self.ds.faces[nf].face_info.pave_blocks_in.insert(g_pb_idx);

            // L1188-1192: Fill map for common blocks creation
            crate::bopalgo::fill_map(&mut a_mpbli, g_pb_idx, nf);
        }

        // L1194-1198: PerformCommonBlocks for coinciding pairs
        if a_mpbli.len() > 0 {
            crate::bopds::tools::perform_common_blocks(self.ds);
        }
    }

    fn force_interf_ve(&mut self) {
        // Build set of existing VE interferences for dedup
        let mut ve_done: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for inf in &self.ds.interf_ve {
            ve_done.insert((inf.vertex, inf.edge));
        }

        for fi in 0..self.ds.faces.len() {
            // Collect all boundary edges of this face (outer + inner wires)
            let face = &self.ds.faces[fi];
            let face_vertices: Vec<usize> = face.face_info.vertices_on
                .iter()
                .chain(face.face_info.vertices_in.iter())
                .copied()
                .collect();
            if face_vertices.is_empty() {
                continue;
            }

            let boundary_edges: Vec<usize> = {
                let f = &self.ds.faces[fi];
                let mut edges = f.boundary_edges.clone();
                for inner in &f.inner_boundary_edges {
                    for &(ei, _) in inner {
                        edges.push(ei);
                    }
                }
                edges
            };

            for &vi in &face_vertices {
                let v_origin = self.ds.vertices[vi].origin;
                if v_origin.is_none() { continue; }
                for &ei in &boundary_edges {
                    let e_origin = self.ds.edges[ei].origin;
                    if e_origin == v_origin.unwrap() { continue; }
                    if self.ds.is_edge_degenerated(ei) { continue; }
                    if ve_done.contains(&(vi, ei)) { continue; }

                    self.check_vertex_edge(vi, ei);
                }
            }
        }
    }

    fn force_interf_vf(&mut self) {
        // Build set of existing VF interferences for dedup
        let mut vf_done: std::collections::HashSet<(usize, usize)> = std::collections::HashSet::new();
        for inf in &self.ds.interf_vf {
            vf_done.insert((inf.vertex, inf.face));
        }

        for vi in 0..self.ds.vertices.len() {
            let v_origin = self.ds.vertices[vi].origin;
            let opposite_faces: Vec<usize> = match v_origin {
                Some(ShapeOrigin::ShapeA) => self.faces_of(ShapeOrigin::ShapeB),
                Some(ShapeOrigin::ShapeB) => self.faces_of(ShapeOrigin::ShapeA),
                _ => continue,
            };

            for &fi in &opposite_faces {
                if vf_done.contains(&(vi, fi)) { continue; }
                self.check_vertex_face(vi, fi);
            }
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
                    let tol = self.ds.faces[fi].geom_tol.max(TOLERANCE_ABS);
                    self.point_on_face(pt, fi, tol)
                });
                if on_face {
                    self.ds.faces[fi].face_info.curves_sc.insert(ci);
                    self.ds.faces[fi].face_info.vertices_in.insert(ic.start_vertex);
                    self.ds.faces[fi].face_info.vertices_in.insert(ic.end_vertex);
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
            Surface3::Cone(c) => { let v = pt - c.apex; let a = c.axis_dir(); let pj = v.dot(a); (v - a*pj).length(); false }
            _ => false,
        }
    }
}
