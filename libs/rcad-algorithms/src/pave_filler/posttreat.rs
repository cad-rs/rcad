//! PostTreatFF family of functions.
//!
//! OCCT ref: BOPAlgo_PaveFiller_6.cxx L1197-1701 (PostTreatFF),
//! L2340-2400 (PutBoundPaveOnCurve), L3105-3310 (ProcessExistingPaveBlocks),
//! L3311-3530 (UpdateExistingPaveBlocks), L3642-3710 (PreparePostTreatFF).

use std::collections::{HashMap, HashSet};
use super::*;
use crate::bvh::{Aabb, BoxTree};
use crate::bopds::pave::*;
use crate::bopds::ds::BOPDS_Iterator;
use rcad_kernel::topods::ShapeType;

impl<'a> PaveFiller<'a> {
    /// BOPAlgo_PaveFiller::PostTreatFF (PaveFiller_6.cxx L1197-1701).
    /// Performs post-treatment of section edges: fuses coincident vertices via
    /// a nested PaveFiller VV pass, creates SD mappings, updates section edge PBs.
    #[allow(non_snake_case)]
    pub(super) fn post_treat_ff(
        &mut self,
        a_mscpb: &mut HashMap<usize, (usize, usize)>,     // theMSCPB: edge/vertex -> (ff_idx, curve_idx)
        a_dm_ex_edges: &mut HashMap<usize, Vec<usize>>,   // aDMExEdges: existing PB -> replacement PBs
        a_dm_new_sd: &mut HashMap<usize, usize>,           // aDMNewSD: old_vertex -> new_vertex
        a_micro_pb: &[PaveBlock],                          // theMicroPB: micro section edges
        a_verts_on_rejected_pb: &HashSet<usize>,           // vertices on rejected PBs
    ) {
        // ===== OCCT L1235-1263: Find unused vertices =====
        let a_nb_s = a_mscpb.len();
        let a_nb_me = a_micro_pb.len();
        let a_nb_v_on_rpb = a_verts_on_rejected_pb.len();
        
        let mut verts_unused: HashSet<usize> = HashSet::new();
        let mut a_ind_map: HashSet<usize> = HashSet::new();
        for inf in &self.ds.interf_ff {
            let (n_f1, n_f2) = (inf.f1, inf.f2);
            let mut a_mv: HashSet<usize> = HashSet::new();
            let mut a_mv_ef: HashSet<usize> = HashSet::new();
            let mut a_mi = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
            let a_mi_b = crate::pave_filler::build_face_shape_map(self.ds, n_f2);
            for v in &a_mi_b { a_mi.insert(*v); }
            // Collect stick vertices
            for inf_ve in &self.ds.interf_ve {
                if !a_mi.contains(&inf_ve.vertex) { continue; }
                a_mv.insert(inf_ve.vertex);
            }
            for inf_vf in &self.ds.interf_vf {
                if !a_mi.contains(&inf_vf.vertex) { continue; }
                a_mv.insert(inf_vf.vertex);
            }
            for inf_ee in &self.ds.interf_ee {
                if inf_ee.new_vertex == usize::MAX { continue; }
                if !a_mi.contains(&inf_ee.e1) || !a_mi.contains(&inf_ee.e2) { continue; }
                let n_v = self.ds.has_shape_sd(inf_ee.new_vertex).unwrap_or(inf_ee.new_vertex);
                a_mv.insert(n_v);
            }
            for inf_vv in &self.ds.interf_vv {
                if inf_vv.merged_vertex == usize::MAX { continue; }
                if !a_mi.contains(&inf_vv.v1) || !a_mi.contains(&inf_vv.v2) { continue; }
                let n_v = self.ds.has_shape_sd(inf_vv.merged_vertex).unwrap_or(inf_vv.merged_vertex);
                a_mv.insert(n_v);
            }
            for inf_ef in &self.ds.interf_ef {
                if inf_ef.new_vertex == usize::MAX { continue; }
                if !a_mi.contains(&inf_ef.edge) || !a_mi.contains(&inf_ef.face) { continue; }
                let n_v = self.ds.has_shape_sd(inf_ef.new_vertex).unwrap_or(inf_ef.new_vertex);
                a_mv.insert(n_v);
                a_mv_ef.insert(n_v);
            }
            // Remove used vertices (those that appear as PB endpoints on intersection curves)
            for &ci in &inf.curves {
                if ci >= self.ds.intersection_curves.len() { continue; }
                let ic = &self.ds.intersection_curves[ci];
                for spb in &ic.pave_blocks {
                    let pb = spb.0.read().unwrap();
                    a_mv.remove(&pb.pave1.vertex_idx);
                    a_mv.remove(&pb.pave2.vertex_idx);
                    for ep in &pb.ext_paves {
                        a_mv.remove(&ep.vertex_idx);
                    }
                }
            }
            // OCCT: IndMap fence — vertices appearing once → VertsUnused, twice → removed
            for &vi in &a_mv {
                if !a_ind_map.insert(vi) {
                    verts_unused.remove(&vi);
                } else {
                    verts_unused.insert(vi);
                }
            }
        }

        // ===== OCCT L1266-1308: Early return for single-entry case =====
        if a_nb_s == 1 && a_nb_me == 0 && a_nb_v_on_rpb == 0 && verts_unused.is_empty() {
            let (&key, &(cur_ind, j)) = a_mscpb.iter().next().unwrap();
            if key < self.ds.edges.len() {
                // It's an edge
                let pb_key = key;
                let (nv1, nv2, a_t1, a_t2, has_edge) = if pb_key < self.ds.pave_blocks.len() {
                    let pb = &self.ds.pave_blocks[pb_key];
                    let r = pb.0.read().unwrap();
                    (r.pave1.vertex_idx, r.pave2.vertex_idx, r.pave1.param, r.pave2.param, r.new_edge.is_some())
                } else if pb_key < self.ds.edges.len() {
                    let pb = self.ds.edges[pb_key].pave_blocks.first()
                        .map(|spb| spb.0.read().unwrap());
                    let r = pb.as_ref().unwrap();
                    (r.pave1.vertex_idx, r.pave2.vertex_idx, r.pave1.param, r.pave2.param, r.new_edge.is_some())
                } else { return; };
                if has_edge {
                    let a_lpbx = vec![pb_key];
                    a_dm_ex_edges.insert(pb_key, a_lpbx);
                } else {
                    let new_ei = self.split_edge(usize::MAX, nv1, a_t1, nv2, a_t2);
                    self.ds.pave_blocks[pb_key].0.write().unwrap().new_edge = Some(new_ei);
                }
            } else {
                // It's a vertex: create a DS copy and set vertex_index on the FFPoint
                let vp = self.ds.vertex_point(key);
                let i_v = self.ds.add_vertex(vp);
                if cur_ind < self.ds.interf_ff.len() {
                    let mut ffp = crate::bopds::ds::types::FFPoint::new(vp, glam::DVec2::ZERO, glam::DVec2::ZERO);
                    ffp.vertex_index = i_v;
                    self.ds.interf_ff[cur_ind].points.push(ffp);
                }
            }
            return;
        }

        // ===== OCCT L1310-1348: Prepare arguments for nested PaveFiller =====
        let mut a_ls: Vec<(usize, usize)> = Vec::new(); // (shape_idx, type: 0=vertex, 1=edge)
        let mut an_existing_edges: HashSet<usize> = HashSet::new();
        let mut an_added_sd: HashSet<usize> = HashSet::new();

        for (&key, &(cur_ind, _j)) in a_mscpb.iter() {
            if key < self.ds.edges.len() {
                let ei = key;
                // Check if PB has an edge (existing edge from source shapes)
                if ei < self.ds.pave_blocks.len() {
                    let pb = &self.ds.pave_blocks[ei];
                    if pb.0.read().unwrap().new_edge.is_some() {
                        // Existing edge
                        an_existing_edges.insert(ei);
                    }
                }
            }
            // Add candidate vertices from aDMNewSD
            let (nv1, nv2) = if key < self.ds.pave_blocks.len() {
                let pb = &self.ds.pave_blocks[key];
                let r = pb.0.read().unwrap();
                (r.pave1.vertex_idx, r.pave2.vertex_idx)
            } else {
                (key, usize::MAX)
            };
            if nv1 < self.ds.vertices.len() {
                let vi = nv1;
                if let Some(&sd_v) = a_dm_new_sd.get(&vi) {
                    if an_added_sd.insert(sd_v) {
                        a_ls.push((sd_v, 0));
                    }
                }
            }
            if nv2 < self.ds.vertices.len() {
                let vi = nv2;
                if let Some(&sd_v) = a_dm_new_sd.get(&vi) {
                    if an_added_sd.insert(sd_v) {
                        a_ls.push((sd_v, 0));
                    }
                }
            }
        }

        // Add existing edges compound
        if !an_existing_edges.is_empty() {
            for &ei in &an_existing_edges {
                a_ls.push((ei, 1));
            }
        }

        // ===== OCCT L1350-1391: Handle micro section edges =====
        for pb in a_micro_pb {
            let (nv1, nv2) = pb.indices();
            let verts = [nv1, nv2];
            for &nv in &verts {
                let sd_vi = a_dm_new_sd.get(&nv).copied().unwrap_or(nv);
                if an_added_sd.insert(sd_vi) {
                    a_ls.push((sd_vi, 0));
                }
            }
            // Increase tolerance to ensure fusion (OCCT: compare points, increase if needed)
            if nv1 < self.ds.vertices.len() && nv2 < self.ds.vertices.len() {
                let p1 = self.ds.vertex_point(nv1);
                let p2 = self.ds.vertex_point(nv2);
                let a_dist = p1.distance(p2);
                let a_tol_v1 = self.ds.vertex_tolerance(nv1);
                let a_tol_v2 = self.ds.vertex_tolerance(nv2);
                let a_dist_adj = a_dist - (a_tol_v1 + a_tol_v2);
                if a_dist_adj > 0.0 {
                    let a_inc = a_dist_adj / 2.0;
                    self.ds.vertex_data_mut(nv1).tolerance = a_tol_v1 + a_inc;
                    self.ds.vertex_data_mut(nv2).tolerance = a_tol_v2 + a_inc;
                }
            }
        }

        // ===== OCCT L1393-1417: Add vertices from rejected PBs and unused =====
        for &vi in a_verts_on_rejected_pb {
            let sd_vi = a_dm_new_sd.get(&vi).copied().unwrap_or(vi);
            if an_added_sd.insert(sd_vi) {
                a_ls.push((sd_vi, 0));
            }
        }
        for &vi in &verts_unused {
            let sd_vi = a_dm_new_sd.get(&vi).copied().unwrap_or(vi);
            if an_added_sd.insert(sd_vi) {
                a_ls.push((sd_vi, 0));
            }
        }

        // ===== OCCT L1387-1398: Nested PaveFiller → aPF.Perform() =====
        // OCCT: BOPAlgo_PaveFiller aPF(theAllocator);
        // OCCT: aPF.SetIsPrimary(false);
        // OCCT: aPF.SetNonDestructive(myNonDestructive);
        // OCCT: aPF.SetRunParallel(myRunParallel);
        // OCCT: aPF.SetArguments(aLS);
        // OCCT: aPF.Perform(aPS.Next());
        // OCCT: aPDS = aPF.PDS();
        if a_ls.is_empty() { return; }

        // Create nested DS (aPDS) and populate with vertices from aLS
        let mut a_pds = DS::new_empty();
        let mut a_map_nested_to_main: Vec<usize> = Vec::new();
        for &(shape_idx, st) in &a_ls {
            if st == 0 && shape_idx < self.ds.vertices.len() {
                let pt = self.ds.vertex_point(shape_idx);
                let tol = self.ds.vertex_tolerance(shape_idx);
                let nested_vi = a_pds.add_vertex(pt);
                if nested_vi < a_pds.vertices.len() {
                    a_pds.vertex_data_mut(nested_vi).tolerance = tol;
                }
                a_map_nested_to_main.push(shape_idx);
            }
        }
        if a_map_nested_to_main.len() < 2 { return; }

        // Set up DS for PaveFiller pipeline
        a_pds.a_vertex_count = a_map_nested_to_main.len() / 2;
        a_pds.nb_source_shapes = a_map_nested_to_main.len();
        for i in 0..a_map_nested_to_main.len() {
            let mut si = crate::bopds::ds::types::ShapeInfo::new(rcad_kernel::topods::ShapeType::Vertex);
            si.is_new = (i >= a_pds.a_vertex_count);
            if i < a_pds.shape_info.len() {
                a_pds.shape_info[i] = si;
            } else {
                a_pds.shape_info.push(si);
            }
        }
        a_pds.build_map_ve();

        // Create nested PaveFiller and run full pipeline
        let mut a_pf = super::PaveFiller::new(&mut a_pds);
        a_pf.is_primary = false;
        a_pf.non_destructive = self.non_destructive;
        a_pf.run_parallel = self.run_parallel;
        a_pf.context.resize(a_pf.ds.faces.len());
        a_pf.perform_body();

        // ===== OCCT L1398-1466: Process results from nested PaveFiller =====
        // OCCT: aItLS.Initialize(aLS); for (; aItLS.More(); aItLS.Next())
        // For each shape in aLS: get type, if VERTEX handle SD or FF point update
        let sd_pairs: Vec<(usize, usize)> = a_pds.shape_sd.sd_vertices_iter()
            .map(|&(a, b)| {
                let main_a = if a < a_map_nested_to_main.len() { a_map_nested_to_main[a] } else { a };
                let main_b = if b < a_map_nested_to_main.len() { a_map_nested_to_main[b] } else { b };
                (main_a, main_b)
            })
            .collect();
        let mut sd_map: std::collections::HashMap<usize, usize> = HashMap::new();
        for &(a, b) in &sd_pairs {
            sd_map.insert(a, b);
        }

        // OCCT L1407-1600: iterate aLS shapes
        for &(shape_idx, st) in &a_ls {
            // OCCT L1411-1418: flatten COMPOUND — not needed in rcad (no compounds in a_ls)

            // OCCT L1419-1422: nSx = aPDS->Index(aSx); aType = aSIx.ShapeType()
            // In rcad, shape_idx is already the main DS index; st is 0=VERTEX, 1=EDGE

            if st == 0 {
                // ===== VERTEX (OCCT L1424-1466) =====
                // OCCT L1426: bIntersectionPoint = theMSCPB.Contains(aSx)
                let b_intersection_point = a_mscpb.contains_key(&shape_idx);
                // OCCT L1428-1435: if (aPDS->HasShapeSD(nSx, nVSD)) aV = aPDS->Shape(nVSD); else aV = aSx
                // rcad: check nested VV pass result — was this vertex fused?
                let n_vx = if let Some(&sd_target) = sd_map.get(&shape_idx) {
                    sd_target
                } else {
                    shape_idx
                };
                // OCCT L1437-1443: iV = myDS->Index(aV); if (iV<0) { append to myDS }
                // rcad: n_vx is already the main DS index (vertices have fixed indices)
                let i_v = n_vx;
                // OCCT L1445-1465: if (!bIntersectionPoint) save SD; else update FF point
                if !b_intersection_point {
                    // OCCT L1448-1453: save SD connection
                    let n_sx = shape_idx;
                    if n_sx != i_v {
                        a_dm_new_sd.insert(n_sx, i_v);
                        self.ds.add_shape_sd(n_sx, i_v);
                    }
                } else {
                    // OCCT L1457-1464: update FF interference point
                    if let Some(&(cur_ind, j)) = a_mscpb.get(&shape_idx) {
                        if cur_ind < self.ds.interf_ff.len() {
                            let a_nb_c = self.ds.interf_ff[cur_ind].curves.len();
                            if j >= a_nb_c {
                                let pt_idx = j - a_nb_c;
                                if pt_idx < self.ds.interf_ff[cur_ind].points.len() {
                                    // OCCT: aNP.SetIndex(iV)
                                    self.ds.interf_ff[cur_ind].points[pt_idx].vertex_index = i_v;
                                }
                            }
                        }
                    }
                }
            } else if st == 1 {
                // ===== EDGE (OCCT L1468-1600) =====
                if shape_idx >= self.ds.edges.len() { continue; }
                // OCCT L1470-1474: get CPB from theMSCPB, iX, iC, aPB1
                let cpb_entry = match a_mscpb.get(&shape_idx) {
                    Some(&val) => val,
                    None => continue,
                };
                let (_i_x, _i_c) = (cpb_entry.0, cpb_entry.1);
                // OCCT L1474: aPB1 = theCPB.PaveBlock1()
                let a_pb1_idx = shape_idx;
                // OCCT L1476: bOld = aPB1->HasEdge()
                let b_old = a_pb1_idx < self.ds.pave_blocks.len()
                    && self.ds.pave_blocks[a_pb1_idx].0.read().unwrap().new_edge.is_some();
                // OCCT L1477-1481: if (bOld) aDMExEdges.Bind(aPB1, aLPBx)
                if b_old {
                    if !a_dm_ex_edges.contains_key(&a_pb1_idx) {
                        a_dm_ex_edges.insert(a_pb1_idx, Vec::new());
                    }
                }
                // OCCT L1470: bHasPaveBlocks = aPDS->HasPaveBlocks(nSx)
                // rcad nested VV pass does not split edges, so bHasPaveBlocks is always false
                let b_has_pave_blocks = false;
                if !b_has_pave_blocks {
                    // OCCT L1483-1497: if (!bHasPaveBlocks)
                    if b_old {
                        // OCCT L1485-1488: aDMExEdges.ChangeFind(aPB1).Append(aPB1)
                        a_dm_ex_edges.get_mut(&a_pb1_idx).unwrap().push(a_pb1_idx);
                    } else {
                        // OCCT L1489-1496: create new edge in myDS
                        // aSI.SetShapeType(aType); aSI.SetShape(aSx);
                        // iE = myDS->Append(aSI); aPB1->SetEdge(iE);
                        // rcad: edge already created by MakeEdge in make_blocks()
                    }
                }
            }
        }

        // ===== OCCT L1690-1700: Follow transitive SD chains =====
        let sd_keys: Vec<usize> = a_dm_new_sd.keys().copied().collect();
        for &k in &sd_keys {
            if let Some(&v) = a_dm_new_sd.get(&k) {
                let mut chain = v;
                let mut depth = 0;
                while let Some(&next) = a_dm_new_sd.get(&chain) {
                    if next == chain || depth > 100 { break; }
                    chain = next;
                    depth += 1;
                }
                if chain != v {
                    a_dm_new_sd.insert(k, chain);
                    self.ds.add_shape_sd(k, chain);
                }
            }
        }
    } // end post_treat_ff

    /// Helper: compute VV tolerance between two vertices in a given DS (possibly nested).
    /// OCCT PerformVV: aTolSum = Tol(V1) + Tol(V2) + myFuzzyValue
    fn vv_pair_tol_ds(&self, ds: &DS, n1: usize, n2: usize) -> f64 {
        let tol1 = if n1 < ds.vertices.len() { ds.vertex_tolerance(n1) } else { TOLERANCE_ABS };
        let tol2 = if n2 < ds.vertices.len() { ds.vertex_tolerance(n2) } else { TOLERANCE_ABS };
        tol1 + tol2 + self.fuzzy_tolerance
    }

    /// BOPAlgo_PaveFiller::PutBoundPaveOnCurve (PaveFiller_6.cxx L2340-2400).
    /// OCCT BOPAlgo_PaveFiller::PutBoundPaveOnCurve (PaveFiller_6.cxx L2340-2399).
    /// Creates new vertices at curve endpoints (if none exist) and adds them as ext_paves.
    pub(super) fn put_bound_pave_on_curve(&mut self, n_f1: usize, n_f2: usize, ci: usize, a_lbv: &mut Vec<usize>) {
        if ci >= self.ds.intersection_curves.len() { return; }
        let (a_t, a_curve, a_geom_tol, a_tang_tol) = {
            let ic = &self.ds.intersection_curves[ci];
            (ic.t_range, ic.curve.clone(), ic.geom_tol, ic.curve_extra.tangential_tol)
        };
        let a_p = [a_curve.point_at(a_t[0]), a_curve.point_at(a_t[1])];
        let a_tol_r3d = a_geom_tol.max(a_tang_tol);

        // getBoundPaves — find extreme ext_paves from PaveBlock1
        let mut a_bnd_nv = [usize::MAX, usize::MAX];
        if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first() {
            let r = pb.0.read().unwrap();
            if !r.ext_paves.is_empty() {
                let min_pave = r.ext_paves.iter().min_by(|a, b| a.param.partial_cmp(&b.param).unwrap());
                let max_pave = r.ext_paves.iter().max_by(|a, b| a.param.partial_cmp(&b.param).unwrap());
                if let Some(p) = min_pave { a_bnd_nv[0] = p.vertex_idx; }
                if let Some(p) = max_pave { a_bnd_nv[1] = p.vertex_idx; }
            }
        }

        // Check if curve is closed (endpoints within tolerance)
        let is_closed = a_p[0].distance_squared(a_p[1]) < TOLERANCE_ABS_SQ;
        if is_closed && (a_bnd_nv[0] != usize::MAX || a_bnd_nv[1] != usize::MAX) {
            return;
        }

        for j in 0..2 {
            if a_bnd_nv[j] != usize::MAX {
                // OCCT L2326-2335: verify bound vertex is near the endpoint
                let n_v = a_bnd_nv[j];
                if n_v < self.ds.vertices.len() {
                    let v = &self.ds.vertices[n_v];
                    let i_flag = crate::boptools::compute_vv_p(v, a_p[j], a_tol_r3d + TOLERANCE_ABS);
                    if i_flag != 0 {
                        // bound vertex NOT at this endpoint → treat as no bound, create new one
                        a_bnd_nv[j] = usize::MAX;
                    } else {
                        continue; // bound vertex exists and matches
                    }
                }
            }

            // OCCT L2372-2397: No bound vertex, create new one if endpoint is valid for both faces
            // Ensure surface cache is populated (OCCT BRepAdaptor_Surface is lazy).
            self.context.surface_adaptor(&self.ds, n_f1);
            self.context.surface_adaptor(&self.ds, n_f2);
            let b_vf = self.context.is_valid_point_for_faces(self.ds, a_p[j], n_f1, n_f2, a_tol_r3d);
            if !b_vf { continue; }

            // Create new vertex at curve endpoint
            // OCCT: myDS->Append(aV) uses TopoDS identity (not position dedup), but the
            // new vertex at the curve endpoint may later be SD-merged with coincident
            // vertices. rcad: use add_vertex (position-based dedup) to avoid creating
            // duplicate vertices at the same position as existing intersection vertices.
            let n_vn = self.ds.add_vertex(a_p[j]);

            // Append ext_pave to PaveBlock1
            if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first_mut() {
                pb.0.write().unwrap().append_ext_pave(Pave { vertex_idx: n_vn, param: a_t[j] });
            }

            a_lbv.push(n_vn);
        }
    }

    /// ProcessExistingPaveBlocks (PaveFiller_6.cxx L3105-3203).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn process_existing_pave_blocks(
        &mut self, ci: usize, j: usize, n_f1: usize, n_f2: usize,
        a_es: usize, a_mpb_on_in: &HashSet<usize>, a_pb_tree: &Option<BoxTree>,
        a_mscpb: &mut HashMap<usize, (usize, usize)>,
        a_mvi: &mut HashMap<usize, usize>, _a_lpb: &mut Vec<PaveBlock>,
        a_pb_faces_map: &mut HashMap<usize, Vec<usize>>,
        a_mpb_add: &mut HashSet<usize>,
    ) {
        // L3087-3088: build Bnd_Box from theES (the new edge a_es)
        if a_es >= self.ds.edges.len() { return; }
        // Build AABB from the edge's curve (equivalent of BRepBndLib::Add)
        let a_box_es = {
            let e = &self.ds.edges[a_es];
            let tol = e.geom_tol;
            let [t0, t1] = e.t_range;
            let n_samples = 8.max(((t1 - t0).abs() / 0.1).ceil() as usize).min(32);
            let mut mn = DVec3::splat(f64::MAX);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for k in 0..=n_samples {
                let t = t0 + (t1 - t0) * (k as f64) / (n_samples as f64);
                let pt = e.curve.point_at(t);
                mn = mn.min(pt);
                mx = mx.max(pt);
            }
            Aabb { min: mn - DVec3::splat(tol), max: mx + DVec3::splat(tol), gap: 0.0 }
        };

        // L3090-3096: BOPTools_BoxTreeSelector query
        let candidates: Vec<usize> = if let Some(pb_tree) = a_pb_tree.as_ref() {
            pb_tree.query_aabb(&a_box_es)
        } else { Vec::new() };
        if candidates.is_empty() { return; }

        // L3098: aTolES = BRep_Tool::Tolerance(theES)
        let a_tol_es = self.ds.edge_tolerance(a_es);

        // L3100-3101: face infos — queried inline to avoid borrow conflicts
        // L3103-3166: for each matching PB from BVH query
        for &pb_idx in &candidates {
            if pb_idx >= self.ds.pave_blocks.len() { continue; }
            // L3106-3109: skip if already in theMPB (a_mpb_add)
            if a_mpb_add.contains(&pb_idx) { continue; }

            // L3111-3112: check face membership
            let b_in_f1 = self.ds.face_info(n_f1).pave_blocks_on.contains(&pb_idx)
                || self.ds.face_info(n_f1).pave_blocks_in.contains(&pb_idx);
            let b_in_f2 = self.ds.face_info(n_f2).pave_blocks_on.contains(&pb_idx)
                || self.ds.face_info(n_f2).pave_blocks_in.contains(&pb_idx);

            if b_in_f1 && b_in_f2 {
                // L3113-3119: both faces — add for post treatment
                a_mpb_add.insert(pb_idx);
                self.prepare_post_treat_ff(ci, j, pb_idx, a_mscpb, a_mvi, ci);
                continue;
            }

            // L3121: one face only
            let n_f = if b_in_f1 { n_f2 } else { n_f1 };
            // L3122-3123: get distance list from myDistances
            // rcad: self.distances is HashMap<(edge_idx, face_idx), Vec<EdgeRangeDistance>>
            let pb_orig_edge = {
                let pb_r = self.ds.pave_blocks[pb_idx].0.read().unwrap();
                pb_r.original_edge
            };
            let p_list = if pb_orig_edge != NO_EDGE {
                self.distances.get(&(pb_orig_edge, n_f))
            } else { None };

            let p_list = match p_list {
                Some(list) => list,
                None => continue,
            };

            // L3129-3130: aPBF->Range(aT1, aT2)
            let (a_t1, a_t2) = {
                let pb_r = self.ds.pave_blocks[pb_idx].0.read().unwrap();
                (pb_r.pave1.param, pb_r.pave2.param)
            };

            // L3132-3144: find distance with range overlap
            let mut a_dist = f64::MAX;
            for range_dist in p_list {
                if (a_t1 <= range_dist.first && range_dist.first <= a_t2)
                    || (a_t1 <= range_dist.last && range_dist.last <= a_t2)
                    || (range_dist.first <= a_t1 && a_t1 <= range_dist.last)
                    || (range_dist.first <= a_t2 && a_t2 <= range_dist.last)
                {
                    a_dist = range_dist.distance;
                    break;
                }
            }

            // L3145-3163: if distance found and within tolerance
            if a_dist < f64::MAX {
                // L3147-3148: aTolSum = aTolES + BRep_Tool::Tolerance(aEF)
                let a_ef_tol = if pb_orig_edge < self.ds.edges.len() {
                    self.ds.edge_tolerance(pb_orig_edge)
                } else { a_tol_es };
                let a_tol_sum = a_tol_es + a_ef_tol;

                if a_dist <= a_tol_sum {
                    // L3152-3153: theMPB.Add + PreparePostTreatFF
                    a_mpb_add.insert(pb_idx);
                    self.prepare_post_treat_ff(ci, j, pb_idx, a_mscpb, a_mvi, ci);

                    // L3155-3163: update aPBFacesMap
                    a_pb_faces_map.entry(pb_idx).or_default().push(n_f);
                }
            }
        }
    }

    /// ProcessExistingPaveBlocks second overload (PaveFiller_6.cxx L3204-3310).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn process_existing_pave_blocks_after(
        &mut self, _a_cur_ind: usize, n_f1: usize, n_f2: usize,
        a_mpb_on_in: &HashSet<usize>, _a_pb_tree: &Option<BoxTree>,
        _a_dmbv: &mut HashMap<usize, Vec<usize>>,
        a_mscpb: &mut HashMap<usize, (usize, usize)>,
        a_mvi: &mut HashMap<usize, usize>,
        a_pb_faces_map: &mut HashMap<usize, Vec<usize>>,
        a_mpb_add: &mut HashSet<usize>,
    ) {
        let se_refs: Vec<Vec<usize>> = (0..self.ds.intersection_curves.len())
            .map(|ci| self.ds.section_edge_refs[ci].clone()).collect();
        for &pb_idx in a_mpb_on_in {
            if pb_idx >= self.ds.pave_blocks.len() || a_mpb_add.contains(&pb_idx) { continue; }
            let (pbsv, pbev) = { let r = self.ds.pave_blocks[pb_idx].0.read().unwrap();
                (r.pave1.vertex_idx, r.pave2.vertex_idx) };
            let mut shares = false;
            for refs in &se_refs {
                for &sei in refs {
                    if sei >= self.ds.edges.len() { continue; }
                    let e = &self.ds.edges[sei];
                    if e.start_vertex == pbsv || e.start_vertex == pbev
                        || e.end_vertex == pbsv || e.end_vertex == pbev { shares = true; break; }
                }
                if shares { break; }
            }
            if !shares { continue; }
            a_mpb_add.insert(pb_idx);
            let b_in_f1 = self.ds.face_info(n_f1).pave_blocks_on.contains(&pb_idx)
                || self.ds.face_info(n_f1).pave_blocks_in.contains(&pb_idx);
            let b_in_f2 = self.ds.face_info(n_f2).pave_blocks_on.contains(&pb_idx)
                || self.ds.face_info(n_f2).pave_blocks_in.contains(&pb_idx);
            if b_in_f1 && b_in_f2 {
                for ci in 0..self.ds.intersection_curves.len() {
                    self.ds.intersection_curves[ci].pave_blocks.push(self.ds.pave_blocks[pb_idx].clone());
                }
                self.ds.face_info_mut(n_f1).pave_blocks_sc.insert(pb_idx);
                self.ds.face_info_mut(n_f2).pave_blocks_sc.insert(pb_idx);
            } else {
                let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                a_pb_faces_map.entry(pb_idx).or_default().push(n_f);
                for ci in 0..self.ds.intersection_curves.len() {
                    self.ds.intersection_curves[ci].pave_blocks.push(self.ds.pave_blocks[pb_idx].clone());
                }
                self.ds.face_info_mut(n_f1).pave_blocks_sc.insert(pb_idx);
                self.ds.face_info_mut(n_f2).pave_blocks_sc.insert(pb_idx);
            }
        }
    }
}