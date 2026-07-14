// OCCT-aligned: BOPAlgo_PaveFiller::MakeBlocks (PaveFiller_6.cxx L650-1169).
// Creates section edges from FF intersection curves and handles post-treatment:
// vertex fusion (PostTreatFF), tolerance correction, face info updates.
//
// Variable naming follows OCCT convention: aXxx = local var matching OCCT aXxx.

use std::collections::{HashMap, HashSet};
use glam::DVec3;
use rcad_kernel::geom::*;
use rcad_kernel::PCurve;
use crate::bvh::{Aabb, DsBvh};
use crate::bopds::ds::{
    DS, DSEdge, DSCurveRepOnFace, Interference, IntersectionCurve, ShapeOrigin,
};
use crate::bopds::pave::*;
use crate::tolerance::*;
use super::helpers::*;
use super::*;

impl<'a> super::PaveFiller<'a> {
    /// OCCT-aligned: BOPAlgo_PaveFiller::CorrectToleranceOfSE (PaveFiller_6.cxx L4105-4306).
    pub(super) fn correct_tolerance_of_se(&mut self) {
        for ci in 0..self.ds.intersection_curves.len() {
            let refs = self.ds.section_edge_refs[ci].clone();
            for &sei in &refs {
                if sei < self.ds.edges.len() {
                    let edge_tol = self.ds.edge_tolerance(sei);
                    let curve_tol = if ci < self.ds.intersection_curves.len() {
                        self.ds.intersection_curves[ci].geom_tol
                    } else { edge_tol };
                    self.ds.edge_data_mut(sei).tolerance = edge_tol.min(curve_tol).max(TOLERANCE_ABS);
                }
            }
        }
    }

    /// OCCT-aligned: GetStickVertices (PaveFiller_6.cxx L2879-2937).
    pub(super) fn get_stick_vertices_ff(
        &self,
        n_f1: usize, n_f2: usize,
        a_mv_stick: &mut HashSet<usize>,
        a_mv_ef: &mut HashSet<usize>,
        a_mi: &mut HashSet<usize>,
    ) {
        a_mi.clear();
        let a_mi_1 = crate::pave_filler::build_face_shape_map(self.ds, n_f1);
        let a_mi_2 = crate::pave_filler::build_face_shape_map(self.ds, n_f2);
        for &v in &a_mi_1 { a_mi.insert(v); }
        for &v in &a_mi_2 { a_mi.insert(v); }
        for inf in &self.ds.interf_ve {
            if !a_mi.contains(&inf.vertex) { continue; }
            a_mv_stick.insert(inf.vertex); a_mi.insert(inf.vertex);
        }
        for inf in &self.ds.interf_vf {
            if !a_mi.contains(&inf.vertex) { continue; }
            a_mv_stick.insert(inf.vertex); a_mi.insert(inf.vertex);
        }
        for inf in &self.ds.interf_ee {
            if inf.new_vertex == usize::MAX { continue; }
            if !a_mi.contains(&inf.e1) || !a_mi.contains(&inf.e2) { continue; }
            let n_v_new = self.ds.has_shape_sd(inf.new_vertex).unwrap_or(inf.new_vertex);
            a_mv_stick.insert(n_v_new); a_mi.insert(n_v_new);
        }
        for inf in &self.ds.interf_vv {
            if inf.merged_vertex == usize::MAX { continue; }
            if !a_mi.contains(&inf.v1) || !a_mi.contains(&inf.v2) { continue; }
            let n_v_new = self.ds.has_shape_sd(inf.merged_vertex).unwrap_or(inf.merged_vertex);
            a_mv_stick.insert(n_v_new); a_mi.insert(n_v_new);
        }
        for inf in &self.ds.interf_ef {
            if inf.new_vertex == usize::MAX { continue; }
            if !a_mi.contains(&inf.edge) || !a_mi.contains(&inf.face) { continue; }
            let n_v_new = self.ds.has_shape_sd(inf.new_vertex).unwrap_or(inf.new_vertex);
            a_mv_stick.insert(n_v_new); a_mv_ef.insert(n_v_new); a_mi.insert(n_v_new);
        }
    }

    /// OCCT-aligned: MakeBlocks (PaveFiller_6.cxx L650-1169).
    #[allow(non_snake_case)]
    pub(super) fn make_blocks(&mut self) {
        // OCCT L654-657: skip if Glue mode
        if self.use_glue() { return; }

        let a_nb_ff_prev = self.ds.interf_ff.len();
        if a_nb_ff_prev == 0 { return; }

        // Cross-iteration collections (persist through entire loop)
        let mut a_mpb_add: HashSet<usize> = HashSet::new();
        let mut a_lpb: Vec<PaveBlock> = Vec::new();
        let mut a_mscpb: HashMap<usize, (usize, usize)> = HashMap::new();  // edge/vertex -> (ff_idx, curve_idx)
        let mut a_mvi: HashMap<usize, usize> = HashMap::new();  // vertex shape -> vertex index
        let mut a_dm_ex_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut a_dm_new_sd: HashMap<usize, usize> = HashMap::new();
        let mut a_dm_vlv: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut a_micro_pb: Vec<PaveBlock> = Vec::new();
        let mut a_verts_on_rejected_pb: HashSet<usize> = HashSet::new();
        let mut a_pb_faces_map: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut a_ff_to_recheck: Vec<usize> = Vec::new();
        let mut a_nb_ff = a_nb_ff_prev;

        // Ensure section_edge_refs is populated
        self.ds.section_edge_refs = vec![Vec::new(); self.ds.intersection_curves.len()];

        // Pre-collect FF data: (f1, f2, curves, points) to avoid borrow conflicts
        let ff_data: Vec<(usize, usize, Vec<usize>, Vec<usize>)> = self.ds.interf_ff.iter()
            .map(|ff| (ff.f1, ff.f2, ff.curves.clone(), ff.points.clone()))
            .collect();

        for i in 0..a_nb_ff {
            let a_cur_ind = if i < a_nb_ff_prev { i } else { a_ff_to_recheck[i - a_nb_ff_prev] };
            let (n_f1, n_f2, curves_of_ff, points_of_ff) = &ff_data[a_cur_ind];
            let (n_f1, n_f2) = (*n_f1, *n_f2);
            let a_nb_c = curves_of_ff.len();
            let a_nb_p = points_of_ff.len();
            if a_nb_p == 0 && a_nb_c == 0 { continue; }

            let a_tol_ff = self.ff_tol(n_f1, n_f2);

            // Per-iteration collections
            let mut a_mv_on_in: HashSet<usize> = HashSet::new();
            let mut a_mv_common: HashSet<usize> = HashSet::new();
            let mut a_mpb_on_in: HashSet<usize> = HashSet::new();
            let mut a_mpb_common: HashSet<usize> = HashSet::new();
            let mut a_dmbv: HashMap<usize, Vec<usize>> = HashMap::new();
            let mut a_mv_tol: Vec<(usize, f64)> = Vec::new();
            let mut a_lse: Vec<usize> = Vec::new();
            let mut a_lbv: Vec<usize> = Vec::new();
            let mut a_mv_stick: HashSet<usize> = HashSet::new();
            let mut a_mv_ef: HashSet<usize> = HashSet::new();
            let mut a_mv_bounds: HashSet<usize> = HashSet::new();
            let mut a_mi: HashSet<usize> = HashSet::new();

            // OCCT L772-773: SubShapesOnIn + SharedEdges
            self.ds.sub_shapes_on_in(n_f1, n_f2, &mut a_mv_on_in, &mut a_mv_common,
                                     &mut a_mpb_on_in, &mut a_mpb_common);
            self.ds.shared_edges(n_f1, n_f2, &mut a_lse);

            // OCCT L775-793: 1. Treat Points (FF point contacts)
            for &pi in points_of_ff {
                if pi >= self.ds.ff_points.len() { continue; }
                let a_p = self.ds.ff_points[pi];
                let b_exist = self.is_existing_vertex_at_point(a_p, a_tol_ff, &a_mv_on_in);
                if !b_exist {
                    let n_v = self.ds.add_vertex(a_p);
                    self.ds.vertex_data_mut(n_v).tolerance = a_tol_ff;
                    a_mscpb.insert(n_v, (a_cur_ind, a_nb_c + pi));  // point index after curves
                }
            }

            // OCCT L796: GetStickVertices
            self.get_stick_vertices_ff(n_f1, n_f2, &mut a_mv_stick, &mut a_mv_ef, &mut a_mi);

            // OCCT L798-829: For each curve — PutPavesOnCurve
            for &ci in curves_of_ff {
                if ci >= self.ds.intersection_curves.len() { continue; }
                // InitPaveBlock1 equivalent
                if self.ds.intersection_curves[ci].pave_blocks.is_empty() {
                    let pb = PaveBlock::new_curve_block();
                    self.ds.intersection_curves[ci].pave_blocks.push(SharedPB::new(pb));
                }
                self.put_paves_on_curve(&a_mv_on_in, &a_mv_common, ci, &a_mi, &a_mv_ef);
            }

            // OCCT L834: FilterPavesOnCurves
            self.filter_paves_on_curves(curves_of_ff);

            // OCCT L836-864: For each curve — stick/EF/boundary paves
            for (j, &ci) in curves_of_ff.iter().enumerate() {
                if ci >= self.ds.intersection_curves.len() { continue; }
                self.put_stick_paves_on_curve(ci, &a_mi, &a_mv_stick);
                if a_nb_c == 1 {
                    self.put_ef_paves_on_curve(ci, &a_mi, &a_mv_ef);
                }
                // PutBoundPaveOnCurve
                a_lbv.clear();
                self.put_bound_pave_on_curve(n_f1, n_f2, ci);
                // Collect bound vertices
                if ci < self.ds.intersection_curves.len() {
                    if let Some(pb) = self.ds.intersection_curves[ci].pave_blocks.first() {
                        let r = pb.0.read().unwrap();
                        for ep in &r.ext_paves {
                            if a_mv_bounds.insert(ep.vertex_idx) {
                                a_lbv.push(ep.vertex_idx);
                            }
                        }
                    }
                }
                if !a_lbv.is_empty() {
                    a_dmbv.insert(j, a_lbv.clone());
                }
            }

            // OCCT L866-871: PutClosingPaveOnCurve for each curve
            for &ci in curves_of_ff {
                if ci >= self.ds.intersection_curves.len() { continue; }
                self.put_closing_pave_on_curve(ci);
            }

            // OCCT L874-895: Build PB tree from aMPBOnIn
            let a_pb_tree = {
                let mut a_pb_indices: Vec<usize> = Vec::new();
                let mut a_pb_aabbs: Vec<Aabb> = Vec::new();
                for &pb_idx in &a_mpb_on_in {
                    if pb_idx >= self.ds.pave_blocks.len() { continue; }
                    let pb = &self.ds.pave_blocks[pb_idx];
                    let r = pb.0.read().unwrap();
                    if r.new_edge.is_none() && r.original_edge == NO_EDGE { continue; }
                    let ei = r.new_edge.unwrap_or(r.original_edge);
                    if ei >= self.ds.edges.len() { continue; }
                    if self.ds.is_edge_degenerated(ei) { continue; }
                    let sv = self.ds.edge_start_vertex_ds(ei);
                    let ev = self.ds.edge_end_vertex_ds(ei);
                    if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
                        let tol = self.ds.edge_tolerance(ei);
                        let mn = self.ds.vertex_point(sv).min(self.ds.vertex_point(ev)) - DVec3::splat(tol);
                        let mx = self.ds.vertex_point(sv).max(self.ds.vertex_point(ev)) + DVec3::splat(tol);
                        a_pb_indices.push(pb_idx);
                        a_pb_aabbs.push(Aabb { min: mn, max: mx });
                    }
                }
                if !a_pb_indices.is_empty() {
                    Some(DsBvh::build(a_pb_indices, a_pb_aabbs))
                } else { None }
            };

            // OCCT L899: isToRecheck flag
            let mut is_to_recheck = a_nb_c > 0 && i < a_nb_ff_prev;

            // OCCT L902-1098: 3. Make section edges
            for (j, &ci) in curves_of_ff.iter().enumerate() {
                if ci >= self.ds.intersection_curves.len() { continue; }

                let a_tol_r3d = {
                    let ic = &self.ds.intersection_curves[ci];
                    ic.geom_tol.max(ic.curve_extra.tangential_tol)
                };

                // aLPB.Clear(); aPB1->Update(aLPB, false);
                a_lpb.clear();
                {
                    let ic = &mut self.ds.intersection_curves[ci];
                    if let Some(mut pb1) = ic.pave_blocks.first().map(|spb| spb.0.write().unwrap()) {
                        // Sync PB endpoints with IC endpoints
                        if pb1.original_edge == NO_EDGE
                            && ic.start_vertex < self.ds.vertices.len()
                            && ic.end_vertex < self.ds.vertices.len()
                        {
                            pb1.pave1 = Pave { vertex_idx: ic.start_vertex, param: ic.t_range[0] };
                            pb1.pave2 = Pave { vertex_idx: ic.end_vertex, param: ic.t_range[1] };
                        }
                        let sub_pbs = pb1.update(false);
                        a_lpb = sub_pbs;
                    }
                }

                if !a_lpb.is_empty() {
                    is_to_recheck = false;
                }

                // Process each sub-PB
                for a_pb in &a_lpb {
                    let (n_v1, n_v2) = a_pb.indices();
                    let (a_t1, a_t2) = a_pb.range();
                    let a_t_range = a_t2 - a_t1;
                    if a_t_range.abs() < 1e-9 { continue; }

                    // OCCT L946-950: IsValidBlockForFaces (midpoint check)
                    let ic = &self.ds.intersection_curves[ci];
                    let mid_t = 0.56786082 * a_t1 + 0.43213918 * a_t2;
                    let mut b_valid_2d = true;
                    for (k, &fi) in [n_f1, n_f2].iter().enumerate() {
                        if fi == usize::MAX { continue; }
                        let pc_opt = if k == 0 { ic.pcurve_on_a.as_ref() } else { ic.pcurve_on_b.as_ref() };
                        let ok = if let Some(pc) = pc_opt {
                            let uv = pc.point_at(mid_t);
                            self.context.is_point_in_on_face(self.ds, fi, uv)
                        } else {
                            let mid_pt = ic.curve.point_at(mid_t);
                            self.context.is_valid_point_for_face(mid_pt, fi, a_tol_r3d)
                        };
                        if !ok { b_valid_2d = false; break; }
                    }
                    if !b_valid_2d { continue; }

                    // OCCT L952-962: IsExistingPaveBlock via LSE (shared edges)
                    let mut n_e_out: usize = usize::MAX;
                    let mut a_tol_new: f64 = -1.0;
                    let b_exist_lse = self.is_existing_pb_via_lse(&a_lse, a_pb, ci, &mut n_e_out, &mut a_tol_new);
                    if b_exist_lse {
                        if a_tol_new > 0.0 && n_e_out < self.ds.edges.len() {
                            self.update_edge_tolerance(n_e_out, a_tol_new);
                        }
                        continue;
                    }

                    // OCCT L967-992: FindValidRange check (micro edge detection)
                    let has_valid_range = {
                        if n_v1 < self.ds.vertices.len() && n_v2 < self.ds.vertices.len() {
                            let v1_pt = self.ds.vertex_point(n_v1);
                            let v2_pt = self.ds.vertex_point(n_v2);
                            let v1_tol = a_tol_r3d.max(self.ds.vertex_tolerance(n_v1));
                            let v2_tol = a_tol_r3d.max(self.ds.vertex_tolerance(n_v2));
                            let ic = &self.ds.intersection_curves[ci];
                            find_valid_range(&ic.curve, a_t1, a_t2, a_tol_r3d,
                                             v1_pt, v1_tol, v2_pt, v2_tol).is_some()
                        } else { false }
                    };
                    if !has_valid_range {
                        if !a_mv_bounds.contains(&n_v1) && !a_mv_bounds.contains(&n_v2) {
                            a_micro_pb.push(a_pb.clone());
                            a_mvi.insert(n_v1, n_v1);
                            a_mvi.insert(n_v2, n_v2);
                        }
                        continue;
                    }

                    // OCCT L994-1052: IsExistingPaveBlock via ON/IN (BVH)
                    let mut a_pb_out: usize = usize::MAX;
                    let mut a_tol_new2: f64 = -1.0;
                    let b_exist_on_in = self.is_existing_pb_via_bvh(
                        a_pb, ci, a_tol_r3d, &a_mpb_on_in, &a_pb_tree,
                        &a_mpb_common, &mut a_pb_out, &mut a_tol_new2,
                    );
                    if b_exist_on_in {
                        if a_pb_out < self.ds.pave_blocks.len() {
                            let edge_of_pb = {
                                let pb_r = self.ds.pave_blocks[a_pb_out].0.read().unwrap();
                                pb_r.new_edge.unwrap_or(pb_r.original_edge)
                            };
                            let b_in_f1 = self.ds.face_info(n_f1).pave_blocks_on.contains(&a_pb_out)
                                || self.ds.face_info(n_f1).pave_blocks_in.contains(&a_pb_out);
                            let b_in_f2 = self.ds.face_info(n_f2).pave_blocks_on.contains(&a_pb_out)
                                || self.ds.face_info(n_f2).pave_blocks_in.contains(&a_pb_out);
                            if !b_in_f1 || !b_in_f2 {
                                // Update edge tolerance
                                if a_tol_new2 >= 0.0 && a_tol_new2 < a_tol_r3d {
                                    a_tol_new2 = a_tol_r3d;
                                }
                                if a_tol_new2 > 0.0 {
                                    self.update_edge_tolerance(edge_of_pb, a_tol_new2);
                                }
                                // PBFacesMap
                                let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                                a_pb_faces_map.entry(a_pb_out).or_default().push(n_f);
                                // Vertices on rejected PB
                                let (n_v_out1, n_v_out2) = {
                                    let pb_r = self.ds.pave_blocks[a_pb_out].0.read().unwrap();
                                    (pb_r.pave1.vertex_idx, pb_r.pave2.vertex_idx)
                                };
                                if n_v1 != n_v_out1 && n_v1 != n_v_out2 && !a_mv_bounds.contains(&n_v1) {
                                    a_verts_on_rejected_pb.insert(n_v1);
                                }
                                if n_v2 != n_v_out1 && n_v2 != n_v_out2 && !a_mv_bounds.contains(&n_v2) {
                                    a_verts_on_rejected_pb.insert(n_v2);
                                }
                                // PreparePostTreatFF
                                if a_mpb_add.insert(a_pb_out) {
                                    if ci < self.ds.intersection_curves.len() {
                                        self.ds.intersection_curves[ci].pave_blocks.push(
                                            self.ds.pave_blocks[a_pb_out].clone()
                                        );
                                    }
                                    a_mscpb.insert(a_pb_out, (a_cur_ind, j));
                                }
                            }
                        }
                        continue;
                    }

                    // OCCT L1055-1094: Make section edge
                    let a_curve = &self.ds.intersection_curves[ci].curve;
                    let pca = self.ds.intersection_curves[ci].pcurve_on_a.clone();
                    let pcb = self.ds.intersection_curves[ci].pcurve_on_b.clone();
                    let new_ei = crate::boptools::make_edge(self.ds, ci, n_v1, n_v2, a_t1, a_t2, a_tol_r3d);
                    crate::boptools::make_pcurve(
                        self.ds, new_ei, n_f1, n_f2, ci,
                        self.section_attribute.pcurve_on_s1,
                        self.section_attribute.pcurve_on_s2,
                        pca.as_ref(), pcb.as_ref(),
                        Some([a_t1, a_t2]), Some([a_t1, a_t2]),
                    );
                    if new_ei < self.ds.edges.len() {
                        if let Some(epb) = self.ds.edge_pave_blocks_mut(new_ei).first_mut() {
                            epb.0.write().unwrap().new_edge = Some(new_ei);
                        }
                    }
                    self.ds.section_edge_refs[ci].push(new_ei);

                    let mut sub_pb = a_pb.clone();
                    sub_pb.new_edge = Some(new_ei);
                    a_mscpb.insert(new_ei, (a_cur_ind, j));
                    a_mvi.insert(n_v1, n_v1);
                    a_mvi.insert(n_v2, n_v2);

                    // Allocate global PB and register on faces
                    let g_pb_idx = self.ds.allocate_pave_block(sub_pb.clone());
                    for &fi in &[n_f1, n_f2] {
                        if fi != usize::MAX {
                            self.ds.face_info_mut(fi).pave_blocks_sc.insert(g_pb_idx);
                        }
                    }

                    // ProcessExistingPaveBlocks (first overload) — register existing PBs
                    let mut tmp_lpb: Vec<PaveBlock> = Vec::new();
                    self.process_existing_pave_blocks(
                        a_cur_ind, j, n_f1, n_f2, new_ei,
                        &a_mpb_on_in, &a_pb_tree, &mut a_mscpb, &mut a_mvi,
                        &mut tmp_lpb, &mut a_pb_faces_map, &mut a_mpb_add,
                    );

                    // Clear MV tolerance for used vertices
                    a_mv_tol.retain(|&(v, _)| v != n_v1 && v != n_v2);
                } // for a_pb in &a_lpb

                // OCCT L1097: aLPBC.RemoveFirst() — remove InitPaveBlock1
                if !self.ds.intersection_curves[ci].pave_blocks.is_empty() {
                    self.ds.intersection_curves[ci].pave_blocks.remove(0);
                }
            } // for (j, &ci)

            // OCCT L1099-1103: isToRecheck
            if is_to_recheck {
                a_ff_to_recheck.push(a_cur_ind);
                a_nb_ff += 1;
            }

            // OCCT L1105-1127: Restore vertex tolerances
            for &(n_v, saved_tol) in &a_mv_tol {
                if n_v < self.ds.vertices.len() {
                    self.ds.vertex_data_mut(n_v).tolerance = saved_tol;
                }
            }
            for &(n_v, _) in &a_mv_tol {
                a_dm_vlv.remove(&n_v);
            }

            // OCCT L1129-1138: ProcessExistingPaveBlocks (second overload)
            self.process_existing_pave_blocks_after(
                a_cur_ind, n_f1, n_f2, &a_mpb_on_in, &a_pb_tree,
                &mut a_dmbv, &mut a_mscpb, &mut a_mvi,
                &mut a_pb_faces_map, &mut a_mpb_add,
            );
        } // for i in 0..a_nb_ff

        // ===== Post-loop phase =====
        // OCCT L1141-1142: RemoveMicroSectionEdges
        self.remove_micro_section_edges(&mut a_mscpb, &mut a_micro_pb);

        // OCCT L1145: MakeSDVerticesFF
        self.make_sd_vertices_ff();

        // OCCT L1146-1152: PostTreatFF (vertex fusion via nested PaveFiller)
        let a_verts_set: HashSet<usize> = a_verts_on_rejected_pb.iter().copied().collect();
        self.post_treat_ff(&mut a_mscpb, &mut a_dm_ex_edges, &mut a_dm_new_sd,
                           &a_micro_pb, &a_verts_set);

        // OCCT L1158: CorrectToleranceOfSE
        self.correct_tolerance_of_se();

        // OCCT L1161: UpdateFaceInfo
        for (&pb_idx, faces) in &a_pb_faces_map {
            if pb_idx < self.ds.pave_blocks.len() {
                for &fi in faces {
                    if fi < self.ds.faces.len() {
                        self.ds.face_info_mut(fi).pave_blocks_sc.insert(pb_idx);
                    }
                }
            }
        }
        for fi in 0..self.ds.faces.len() {
            let curves_sc: Vec<usize> = self.ds.face_info(fi).curves_sc_only().iter().copied().collect();
            for &ci in &curves_sc {
                if ci < self.ds.intersection_curves.len() {
                    let (sv, ev) = {
                        let ic = &self.ds.intersection_curves[ci];
                        (ic.start_vertex, ic.end_vertex)
                    };
                    self.ds.face_info_mut(fi).vertices_in.insert(sv);
                    self.ds.face_info_mut(fi).vertices_in.insert(ev);
                }
            }
        }

        // OCCT L1163: UpdatePaveBlocks
        for (old_v, new_v) in &a_dm_new_sd {
            for ei in 0..self.ds.edges.len() {
                for spb in &mut self.ds.edges[ei].pave_blocks {
                    let mut pb = spb.0.write().unwrap();
                    if pb.pave1.vertex_idx == *old_v { pb.pave1.vertex_idx = *new_v; }
                    if pb.pave2.vertex_idx == *old_v { pb.pave2.vertex_idx = *new_v; }
                }
            }
            for fi in 0..self.ds.faces.len() {
                if self.ds.face_info(fi).vertices_in.contains(old_v) {
                    self.ds.face_info_mut(fi).vertices_in.remove(old_v);
                    self.ds.face_info_mut(fi).vertices_in.insert(*new_v);
                }
            }
        }

        // Remove micro edges
        let micro_edge_set: HashSet<usize> = a_micro_pb.iter()
            .filter_map(|pb| pb.new_edge)
            .collect();
        if !micro_edge_set.is_empty() {
            self.remove_pave_blocks(&micro_edge_set);
        }

        // OCCT L1168: PutSEInOtherFaces
        self.put_se_in_other_faces();

        // Build edge images
        self.ds.build_edge_images();
    }

    /// OCCT-aligned: IsExistingPaveBlock via LSE (shared edges).
    /// OCCT PaveFiller_6.cxx L952-961.
    fn is_existing_pb_via_lse(
        &self, a_lse: &[usize], a_pb: &PaveBlock, ci: usize,
        n_e_out: &mut usize, a_tol_new: &mut f64,
    ) -> bool {
        if a_lse.is_empty() { return false; }
        let (n_v1, n_v2) = a_pb.indices();
        let (a_t1, a_t2) = a_pb.range();
        let a_tm = 0.56786082 * a_t1 + 0.43213918 * a_t2;
        let a_pm = if ci < self.ds.intersection_curves.len() {
            self.ds.intersection_curves[ci].curve.point_at(a_tm)
        } else { return false };
        let a_tol = {
            let v1_tol = if n_v1 < self.ds.vertices.len() { self.ds.vertex_tolerance(n_v1) } else { TOLERANCE_ABS };
            let v2_tol = if n_v2 < self.ds.vertices.len() { self.ds.vertex_tolerance(n_v2) } else { TOLERANCE_ABS };
            v1_tol.max(v2_tol)
        };
        let mut found = false;
        let mut best_dist = f64::MAX;
        for &sei in a_lse {
            if sei >= self.ds.edges.len() { continue; }
            let se = &self.ds.edges[sei];
            let a_tol_check = se.geom_tol.max(a_tol);
            let (_t, a_proj) = crate::extrema::closest_point_on_curve(&se.curve, a_pm);
            let dist = (a_proj - a_pm).length();
            if dist <= a_tol_check && dist < best_dist {
                found = true;
                *n_e_out = sei;
                *a_tol_new = dist;
                best_dist = dist;
            }
        }
        found
    }

    /// OCCT-aligned: IsExistingPaveBlock via ON/IN + BVH.
    /// OCCT PaveFiller_6.cxx L994-1052.
    #[allow(clippy::too_many_arguments)]
    fn is_existing_pb_via_bvh(
        &self, a_pb: &PaveBlock, ci: usize, a_tol_r3d: f64,
        a_mpb_on_in: &HashSet<usize>,
        a_pb_tree: &Option<DsBvh>,
        a_mpb_common: &HashSet<usize>,
        a_pb_out: &mut usize,
        a_tol_new: &mut f64,
    ) -> bool {
        let (n_v1, n_v2) = a_pb.indices();
        let (a_t1, a_t2) = a_pb.range();
        let a_tm = 0.56786082 * a_t1 + 0.43213918 * a_t2;
        let a_pm = if ci < self.ds.intersection_curves.len() {
            self.ds.intersection_curves[ci].curve.point_at(a_tm)
        } else { return false };
        let a_p1 = self.ds.intersection_curves[ci].curve.point_at(a_t1);
        let a_p2 = self.ds.intersection_curves[ci].curve.point_at(a_t2);

        let a_tol_v11 = if n_v1 < self.ds.vertices.len() { self.ds.vertex_tolerance(n_v1) } else { a_tol_r3d };
        let a_tol_v12 = if n_v2 < self.ds.vertices.len() { self.ds.vertex_tolerance(n_v2) } else { a_tol_r3d };
        let a_tol_v1 = a_tol_v11.max(a_tol_v12);
        let a_tol_check = a_tol_r3d;

        // Query BVH
        let candidates: Vec<usize> = if let Some(pb_tree) = a_pb_tree.as_ref() {
            let query_box = Aabb {
                min: a_pm - DVec3::splat(a_tol_v1 + a_tol_check),
                max: a_pm + DVec3::splat(a_tol_v1 + a_tol_check),
            };
            pb_tree.query_aabb(&query_box)
        } else { Vec::new() };

        if candidates.is_empty() { return false; }

        let mut found_pb = usize::MAX;
        let mut best_dist = f64::MAX;
        let mut best_tol = -1.0;

        for &pb_idx in &candidates {
            if pb_idx >= self.ds.pave_blocks.len() { continue; }
            let existing_pb = &self.ds.pave_blocks[pb_idx];
            let (n_v21, n_v22) = existing_pb.0.read().unwrap().indices();
            let a_tol_v21 = if n_v21 < self.ds.vertices.len() { self.ds.vertex_tolerance(n_v21) } else { a_tol_r3d };
            let a_tol_v22 = if n_v22 < self.ds.vertices.len() { self.ds.vertex_tolerance(n_v22) } else { a_tol_r3d };
            let a_tol_v2 = a_tol_v21.max(a_tol_v22);
            let i_flag1 = n_v1 == n_v21 || n_v1 == n_v22;

            let edge_ei = existing_pb.0.read().unwrap().new_edge
                .unwrap_or(existing_pb.0.read().unwrap().original_edge);
            let i_flag2 = if n_v2 == n_v21 || n_v2 == n_v22 {
                true
            } else if edge_ei < self.ds.edges.len() {
                let sv = self.ds.edge_start_vertex_ds(edge_ei);
                let ev = self.ds.edge_end_vertex_ds(edge_ei);
                let e_min = if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
                    self.ds.vertex_point(sv).min(self.ds.vertex_point(ev))
                } else { a_p2 };
                let e_max = if sv < self.ds.vertices.len() && ev < self.ds.vertices.len() {
                    self.ds.vertex_point(sv).max(self.ds.vertex_point(ev))
                } else { a_p2 };
                let e_tol = a_tol_v21.max(a_tol_v22);
                let sp_min = e_min - DVec3::splat(e_tol);
                let sp_max = e_max + DVec3::splat(e_tol);
                let p2_min = a_p2 - DVec3::splat(a_tol_v12);
                let p2_max = a_p2 + DVec3::splat(a_tol_v12);
                !(sp_max.x < p2_min.x || sp_min.x > p2_max.x
                    || sp_max.y < p2_min.y || sp_min.y > p2_max.y
                    || sp_max.z < p2_min.z || sp_min.z > p2_max.z)
            } else { false };
            if !i_flag2 { continue; }

            if edge_ei >= self.ds.edges.len() { continue; }
            let existing_edge = &self.ds.edges[edge_ei];
            let mut a_real_tol = a_tol_check;
            if a_mpb_common.contains(&pb_idx) {
                a_real_tol = a_real_tol.max(a_tol_v1.max(a_tol_v2));
                a_real_tol *= 2.0;
            }
            let (_t, proj) = crate::extrema::closest_point_on_curve(&existing_edge.curve, a_pm);
            let dist_to_sp = (proj - a_pm).length();
            if dist_to_sp > a_real_tol { continue; }

            let mut dist_p1 = f64::MAX;
            if !i_flag1 {
                let (_t1, p1_proj) = crate::extrema::closest_point_on_curve(&existing_edge.curve, a_p1);
                dist_p1 = (p1_proj - a_p1).length();
            }
            let mut dist_to_use = dist_to_sp;
            if n_v2 != n_v21 && n_v2 != n_v22 {
                let (_t2, p2_proj) = crate::extrema::closest_point_on_curve(&existing_edge.curve, a_p2);
                let dist_p2 = (p2_proj - a_p2).length();
                if dist_to_use < dist_p2 { dist_to_use = dist_p2; }
            }
            let i_flag1_ok = i_flag1 || dist_p1 <= a_real_tol;
            if i_flag1_ok && dist_to_use < best_dist {
                found_pb = pb_idx;
                best_tol = dist_to_use;
                best_dist = dist_to_use;
            }
        }

        if found_pb != usize::MAX {
            *a_pb_out = found_pb;
            *a_tol_new = best_tol;
            true
        } else { false }
    }
}
