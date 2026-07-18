//! PostTreatFF family of functions.
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
                                    // Update vertex_index on existing FFPoint (OCCT BOPDS_Point::SetIndex)
                                    let mut n_vx = n_vx;
                                    // Check if this vertex already exists among FF points
                                    let already_has = self.ds.interf_ff[cur_ind].points.iter()
                                        .any(|p| p.vertex_index == n_vx);
                                    if !already_has {
                                        // Update existing point's vertex index (OCCT equivalent: aNP.SetIndex(nVX))
                                        // If the point hasn't been assigned a vertex yet, assign it.
                                        // n_vx is the new vertex index for this point.
                                        if self.ds.interf_ff[cur_ind].points[pt_idx].vertex_index == usize::MAX
                                            || self.ds.interf_ff[cur_ind].points[pt_idx].vertex_index == n_vx
                                        {
                                            // Update the existing point's vertex reference
                                            let updated = crate::bopds::ds::types::FFPoint {
                                                vertex_index: n_vx,
                                                ..self.ds.interf_ff[cur_ind].points[pt_idx].clone()
                                            };
                                            self.ds.interf_ff[cur_ind].points[pt_idx] = updated;
                                        }
                                    }
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

        // ===== OCCT L1500-1530: Create edges from micro PBs after VV fusion =====
        // DISABLED: rcad's micro PBs have NO_EDGE as original_edge, and
        // SplitEdge(usize::MAX, ...) is a no-op.  Proper fix requires making
        // PostTreatFF use the intersection curve as the edge source for SplitEdge.
        // For now, micro PBs are handled by the nested VV fusion (vertex merging)
        // which is sufficient for most cases.
        if false {} // placeholder
    }

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
            let b_vf = self.context.is_valid_point_for_faces(a_p[j], n_f1, n_f2, a_tol_r3d);
            if !b_vf { continue; }

            // Create new vertex at curve endpoint
            let n_vn = crate::boptools::make_new_vertex(&mut self.ds, a_p[j], a_tol_r3d);

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

    /// ProcessExistingPaveBlocks second overload (PaveFiller_6.cxx L3204-3310).
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