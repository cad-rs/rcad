//! OCCT-aligned: PostTreatFF family of functions.
//!
//! OCCT ref: BOPAlgo_PaveFiller_6.cxx L1197-1701 (PostTreatFF),
//! L2340-2400 (PutBoundPaveOnCurve), L3105-3310 (ProcessExistingPaveBlocks),
//! L3311-3530 (UpdateExistingPaveBlocks), L3642-3710 (PreparePostTreatFF).

use std::collections::{HashMap, HashSet};
use super::*;
use crate::bopds::pave::*;
use crate::bopds::ds::BOPDS_Iterator;
use rcad_kernel::topods::ShapeType;

impl<'a> PaveFiller<'a> {
    /// OCCT-aligned: BOPAlgo_PaveFiller::PostTreatFF (PaveFiller_6.cxx L1197-1701).
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
                // It's a vertex
                let i_v = self.ds.add_vertex(self.ds.vertex_point(key));
                if cur_ind < self.ds.interf_ff.len() {
                    self.ds.interf_ff[cur_ind].points.push(i_v);
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

        // ===== OCCT L1419-1430: Nested PaveFiller → VV pass for vertex fusion =====
        if a_ls.is_empty() { return; }

        // Create nested DS and populate it with the argument vertices
        let mut nested_ds = DS::new_empty();
        // Populate vertices from the argument list
        let mut nested_vertex_map: Vec<usize> = Vec::new(); // nested_idx → main_ds_idx
        for &(shape_idx, st) in &a_ls {
            if st == 0 && shape_idx < self.ds.vertices.len() {
                let pt = self.ds.vertex_point(shape_idx);
                let tol = self.ds.vertex_tolerance(shape_idx);
                let nested_vi = nested_ds.add_vertex(pt);
                if nested_vi < nested_ds.vertices.len() {
                    nested_ds.vertex_data_mut(nested_vi).tolerance = tol;
                }
                nested_vertex_map.push(shape_idx);
            }
        }

        if nested_vertex_map.len() < 2 {
            // No VV pairs possible (need at least 2 vertices)
            return;
        }

        // Build iterator for cross-group VV pairs
        nested_ds.a_vertex_count = nested_vertex_map.len() / 2;
        nested_ds.nb_source_shapes = nested_vertex_map.len();
        // Add shape_info entries for each vertex
        for i in 0..nested_vertex_map.len() {
            let mut si = crate::bopds::ds::types::ShapeInfo::new(rcad_kernel::topods::ShapeType::Vertex);
            si.is_new = (i >= nested_ds.a_vertex_count);
            if i < nested_ds.shape_info.len() {
                nested_ds.shape_info[i] = si;
            } else {
                nested_ds.shape_info.push(si);
            }
        }
        nested_ds.build_map_ve();

        // Run VV pass via BOPDS_Iterator
        let vv_pairs: Vec<(usize, usize)> = {
            use crate::bopds::ds::types::PairIterator;
            let mut pairs = Vec::new();
            // Generate cross-group VV pairs
            let n_a = nested_ds.a_vertex_count;
            let n_b = nested_ds.vertices.len();
            for i in 0..n_a {
                for j in n_a..n_b {
                    pairs.push((i, j));
                }
            }
            pairs
        };

        // Perform VV: check distances and create SD mappings in nested DS
        for &(n1, n2) in &vv_pairs {
            let tol = if n1 < nested_ds.vertices.len() && n2 < nested_ds.vertices.len() {
                self.vv_pair_tol_ds(&nested_ds, n1, n2)
            } else { TOLERANCE_ABS };
            let dist = if n1 < nested_ds.vertices.len() && n2 < nested_ds.vertices.len() {
                (nested_ds.vertex_point(n1) - nested_ds.vertex_point(n2)).length()
            } else { f64::MAX };
            if dist <= tol {
                // Create SD mapping in nested DS
                nested_ds.add_shape_sd(n1, n2);
            }
        }

        // ===== Translate nested DS SD mappings back to main DS =====
        let sd_pairs: Vec<(usize, usize)> = nested_ds.shape_sd.sd_vertices_iter()
            .map(|&(a, b)| {
                let main_a = if a < nested_vertex_map.len() { nested_vertex_map[a] } else { a };
                let main_b = if b < nested_vertex_map.len() { nested_vertex_map[b] } else { b };
                (main_a, main_b)
            })
            .collect();

        // Register SD mappings in main DS and update aDMNewSD
        for &(main_a, main_b) in &sd_pairs {
            if main_a != main_b && main_a < self.ds.vertices.len() && main_b < self.ds.vertices.len() {
                let n_vx = main_b;
                // Register in SD map
                self.ds.add_shape_sd(main_a, n_vx);
                a_dm_new_sd.insert(main_a, n_vx);
                
                // Update FF interference point indices if this vertex was an intersection point
                let b_intersection_point = a_mscpb.contains_key(&main_a);
                if b_intersection_point {
                    if let Some(&(cur_ind, j)) = a_mscpb.get(&main_a) {
                        if cur_ind < self.ds.interf_ff.len() {
                            // Check if j is a point index (not curve index)
                            let a_nb_c = self.ds.interf_ff[cur_ind].curves.len();
                            if j >= a_nb_c {
                                let pt_idx = j - a_nb_c; // point index within this FF
                                if pt_idx < self.ds.interf_ff[cur_ind].points.len() {
                                    self.ds.interf_ff[cur_ind].points[pt_idx] = n_vx;
                                }
                            }
                        }
                    }
                }
            }
        }

        // ===== OCCT L1690-1700: Follow transitive SD chains =====
        let sd_keys: Vec<usize> = a_dm_new_sd.keys().copied().collect();
        for &k in &sd_keys {
            if let Some(&v) = a_dm_new_sd.get(&k) {
                // Follow chain
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
    }

    /// Helper: compute VV tolerance between two vertices in a given DS (possibly nested).
    /// OCCT-aligned: equivalent to the VV pair tolerance computation.
    fn vv_pair_tol_ds(&self, ds: &DS, n1: usize, n2: usize) -> f64 {
        let tol1 = if n1 < ds.vertices.len() { ds.vertex_tolerance(n1) } else { TOLERANCE_ABS };
        let tol2 = if n2 < ds.vertices.len() { ds.vertex_tolerance(n2) } else { TOLERANCE_ABS };
        let a_tol_vv = tol1.max(tol2);
        a_tol_vv + self.fuzzy_tolerance
    }

    /// OCCT-aligned: BOPAlgo_PaveFiller::PutBoundPaveOnCurve (PaveFiller_6.cxx L2340-2400).
    pub(super) fn put_bound_pave_on_curve(&mut self, n_f1: usize, n_f2: usize, ci: usize) {
        if ci >= self.ds.intersection_curves.len() { return; }
        let ic_data = {
            let ic = &self.ds.intersection_curves[ci];
            (ic.curve.clone(), ic.t_range, ic.geom_tol,
             ic.pcurve_on_a.clone(), ic.pcurve_on_b.clone())
        };
        let (a_curve, a_t_range, a_geom_tol, pcurve_on_a, pcurve_on_b) = ic_data;
        let a_tol_r3d = a_geom_tol.max(TOLERANCE_ABS);
        for (k, &fi) in ([n_f1, n_f2]).iter().enumerate() {
            if fi >= self.ds.faces.len() { continue; }
            let pc = if k == 0 { pcurve_on_a.as_ref() } else { pcurve_on_b.as_ref() };
            let Some(pc) = pc else { continue };
            let tt0 = a_t_range[0]; let tt1 = a_t_range[1];
            let span = tt1 - tt0;
            if span <= TOLERANCE_CLAMP_MIN { continue; }
            let n_samp = 129;
            let first_state = self.context.is_point_in_on_face(self.ds, fi, pc.point_at(tt0));
            let mut prev_t = tt0; let mut prev_state = first_state;
            for i in 1..=n_samp {
                let t = tt0 + span * i as f64 / n_samp as f64;
                let state = self.context.is_point_in_on_face(self.ds, fi, pc.point_at(t));
                if state != prev_state {
                    let mut lo = prev_t; let mut hi = t;
                    for _ in 0..20 {
                        let mid = (lo + hi) * 0.5;
                        let mid_state = self.context.is_point_in_on_face(self.ds, fi, pc.point_at(mid));
                        if mid_state == prev_state { lo = mid; } else { hi = mid; }
                    }
                    let ct = (lo + hi) * 0.5;
                    let cp = a_curve.point_at(ct);
                    let nv = self.ds.add_vertex(cp);
                    self.ds.vertex_data_mut(nv).tolerance = a_tol_r3d;
                    if ci < self.ds.intersection_curves.len() {
                        if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first_mut() {
                            pb.0.write().unwrap().append_ext_pave(Pave { vertex_idx: nv, param: ct });
                        }
                    }
                    prev_state = state;
                }
                prev_t = t;
            }
        }
    }

    /// OCCT-aligned: ProcessExistingPaveBlocks (PaveFiller_6.cxx L3105-3203).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn process_existing_pave_blocks(
        &mut self, a_cur_ind: usize, _j: usize, n_f1: usize, n_f2: usize,
        a_es: usize, a_mpb_on_in: &HashSet<usize>, _a_pb_tree: &Option<DsBvh>,
        a_mscpb: &mut HashMap<usize, (usize, usize)>,
        a_mvi: &mut HashMap<usize, usize>, _a_lpb: &mut Vec<PaveBlock>,
        a_pb_faces_map: &mut HashMap<usize, Vec<usize>>,
        a_mpb_add: &mut HashSet<usize>,
    ) {
        if a_es >= self.ds.edges.len() { return; }
        let (sv, ev) = { let e = &self.ds.edges[a_es]; (e.start_vertex, e.end_vertex) };
        for &pb_idx in a_mpb_on_in {
            if pb_idx >= self.ds.pave_blocks.len() || a_mpb_add.contains(&pb_idx) { continue; }
            let (pbsv, pbev) = { let r = self.ds.pave_blocks[pb_idx].0.read().unwrap();
                (r.pave1.vertex_idx, r.pave2.vertex_idx) };
            if pbsv != sv && pbsv != ev && pbev != sv && pbev != ev { continue; }
            a_mpb_add.insert(pb_idx);
            let b_in_f1 = self.ds.face_info(n_f1).pave_blocks_on.contains(&pb_idx)
                || self.ds.face_info(n_f1).pave_blocks_in.contains(&pb_idx);
            let b_in_f2 = self.ds.face_info(n_f2).pave_blocks_on.contains(&pb_idx)
                || self.ds.face_info(n_f2).pave_blocks_in.contains(&pb_idx);
            if b_in_f1 && b_in_f2 {
                if let Some(ic) = self.ds.intersection_curves.get_mut(a_cur_ind) {
                    ic.pave_blocks.push(self.ds.pave_blocks[pb_idx].clone());
                }
                self.ds.face_info_mut(n_f1).pave_blocks_sc.insert(pb_idx);
                self.ds.face_info_mut(n_f2).pave_blocks_sc.insert(pb_idx);
            } else {
                let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                a_pb_faces_map.entry(pb_idx).or_default().push(n_f);
                if let Some(ic) = self.ds.intersection_curves.get_mut(a_cur_ind) {
                    ic.pave_blocks.push(self.ds.pave_blocks[pb_idx].clone());
                }
                self.ds.face_info_mut(n_f1).pave_blocks_sc.insert(pb_idx);
                self.ds.face_info_mut(n_f2).pave_blocks_sc.insert(pb_idx);
            }
        }
    }

    /// OCCT-aligned: ProcessExistingPaveBlocks second overload (PaveFiller_6.cxx L3204-3310).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn process_existing_pave_blocks_after(
        &mut self, _a_cur_ind: usize, n_f1: usize, n_f2: usize,
        a_mpb_on_in: &HashSet<usize>, _a_pb_tree: &Option<DsBvh>,
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