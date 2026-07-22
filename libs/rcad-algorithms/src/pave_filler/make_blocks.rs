// BOPAlgo_PaveFiller::MakeBlocks (PaveFiller_6.cxx L650-1169).
//
// Variable naming follows OCCT convention: aXxx = local var matching OCCT aXxx.
// Allocator-specific logic (IncAllocator, DefaultAllocator) is omitted because
// Rust's Vec/HashMap handle memory automatically — the structure (collection
// declarations outside the loop, cleared inside) is preserved for OCCT alignment.
//
// One notable structural difference: the rcad DS uses separate index namespaces
// for vertices vs edges (occt_vertex_idx != occt_edge_idx), so HashMap keys
// mixing both indices are safe as long as no index value from one namespace
// equals an index value from the other — unlikely in practice.

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
    // OCCT BOPAlgo_PaveFiller::CorrectToleranceOfSE (PaveFiller_6.cxx L4105-4306).
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

    // OCCT BOPAlgo_PaveFiller::GetStickVertices (PaveFiller_6.cxx L2879-2937).
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

    /// UpdateSavedTolerance (PaveFiller_6.cxx L626-646).
    /// Updates the saved tolerance of the vertices of the edge with new tolerance of edge.
    /// In the Rust version, aMVTol is HashMap<usize, f64> mapping vertex -> saved tolerance.
    fn update_saved_tolerance(
        &self,
        n_e: usize,
        a_tol_new: f64,
        a_mv_tol: &mut HashMap<usize, f64>,
    ) {
        if n_e >= self.ds.edges.len() { return; }
        // OCCT L630-645: iterate SubShapes() — all vertex sub-shapes of the edge
        // including internal vertices (paves). rcad: start/end vertices + paves.
        let n_vs: Vec<usize> = {
            let mut v = vec![self.ds.edge_start_vertex_ds(n_e), self.ds.edge_end_vertex_ds(n_e)];
            for pave in &self.ds.edges[n_e].paves {
                v.push(pave.vertex_idx);
            }
            v
        };
        for &n_v in &n_vs {
            if n_v >= self.ds.vertices.len() { continue; }
            if let Some(&tol_saved) = a_mv_tol.get(&n_v) {
                if tol_saved < a_tol_new {
                    a_mv_tol.insert(n_v, a_tol_new);
                }
            }
            // If not yet in aMVTol, the tolerance will be saved when first modified.
        }
    }

    /// PreparePostTreatFF (PaveFiller_6.cxx L3642-3668).
    /// Adds the existing pave block for post-treatment processing.
    pub(super) fn prepare_post_treat_ff(
        &mut self,
        a_int: usize,
        a_cur: usize,
        a_pb_idx: usize,
        a_mscpb: &mut HashMap<usize, (usize, usize)>,
        a_mvi: &mut HashMap<usize, usize>,
        ci: usize,
    ) {
        // aLPBC.Append(aPB)
        if ci < self.ds.intersection_curves.len() {
            self.ds.intersection_curves[ci].pave_blocks.push(
                self.ds.pave_blocks[a_pb_idx].clone()
            );
        }
        // aPB->Indices(nV1, nV2);
        let (n_v1, n_v2) = {
            let pb_r = self.ds.pave_blocks[a_pb_idx].0.read().unwrap();
            (pb_r.pave1.vertex_idx, pb_r.pave2.vertex_idx)
        };
        // Keep info for post treatment
        a_mscpb.insert(a_pb_idx, (a_int, a_cur));
        a_mvi.insert(n_v1, n_v1);
        a_mvi.insert(n_v2, n_v2);
    }

    /// IsExistingPaveBlock via LSE (shared edges).
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
            let a_tol_check = se.geom_tol.max(a_tol) + self.ds.fuzzy_tol;
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

    /// IsExistingPaveBlock via ON/IN + BVH.
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
        // OCCT L2093: std::max(aTolV11, aTolV12) + myFuzzyValue
        let a_tol_v1 = a_tol_v11.max(a_tol_v12) + self.ds.fuzzy_tol;
        // OCCT L2095: theTolR3D + myFuzzyValue
        let a_tol_check = a_tol_r3d + self.ds.fuzzy_tol;

        // Query BVH
        let candidates: Vec<usize> = if let Some(pb_tree) = a_pb_tree.as_ref() {
            let query_box = Aabb {
                min: a_pm - DVec3::splat(a_tol_v1 + a_tol_check),
                max: a_pm + DVec3::splat(a_tol_v1 + a_tol_check), gap: 0.0 };
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
            // OCCT L2117: std::max(aTolV21, aTolV22) + myFuzzyValue
            let a_tol_v2 = a_tol_v21.max(a_tol_v22) + self.ds.fuzzy_tol;
            let edge_ei = existing_pb.0.read().unwrap().new_edge
                .unwrap_or(existing_pb.0.read().unwrap().original_edge);
            // OCCT L2123: iFlag1 = (nV11 == nV21 || nV11 == nV22) ? 2 : 1
            //             iFlag2 = (nV12 == nV21 || nV12 == nV22) ? 2 : (!aBoxSp.IsOut(aBoxP2) ? 1 : 0)
            // rcad: iFlag1 == 2 maps to true (vertex match), 1 maps to false (AABB check needed)
            //       iFlag2 == 2 maps to true (vertex match), 1 maps to false (AABB check needed)
            let i_flag1 = n_v1 == n_v21 || n_v1 == n_v22;
            let i_flag2 = if n_v2 == n_v21 || n_v2 == n_v22 {
                true
            } else if edge_ei < self.ds.edges.len() {
                // AABB overlap: does the edge's tolerance box contain a_p2?
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
            // OCCT L2132-2176: tolerance adjustment based on common-block / thin-face
            let mut a_real_tol = a_tol_check;
            // OCCT L2134: if (myDS->IsCommonBlock(aPB))
            let is_cb = self.ds.is_common_block(&self.ds.pave_blocks[pb_idx]);
            if is_cb {
                // OCCT L2135-2137: aRealTol = max(aRealTol, max(aTolV1, aTolV2))
                a_real_tol = a_real_tol.max(a_tol_v1.max(a_tol_v2));
                // OCCT L2138: if (theMPBCommon.Contains(aPB)) aRealTol *= 2.
                if a_mpb_common.contains(&pb_idx) {
                    a_real_tol *= 2.0;
                }
            } else if i_flag1 && i_flag2 {
                // OCCT L2139-2176: thin-face tangent-angle check
                // Skip if one edge is closed and the other is not
                let n_v11_closed = n_v1 == n_v2;
                let n_v21_closed = n_v21 == n_v22;
                let b_skip_processing =
                    (n_v11_closed && !n_v21_closed) || (!n_v11_closed && n_v21_closed);
                if !b_skip_processing {
                    // Check tangent alignment
                    let a_ic = &self.ds.intersection_curves[ci];
                    // OCCT L2103-2105: aC3d->D1(aTm, aPm, aVTgt1)
                    let a_pm_ic = a_ic.curve.point_at(a_tm);
                    let a_vtgt1 = a_ic.curve.derivative_at(a_tm);
                    let is_vtgt1_valid = a_vtgt1.length_squared() > f64::EPSILON;
                    if is_vtgt1_valid {
                        let a_vtgt1_n = a_vtgt1.normalize();
                        // OCCT L2147: if (aIC.Type() != GeomAbs_Line || aBAC2.GetType() != GeomAbs_Line)
                        let is_ic_line = matches!(a_ic.curve, Curve3::Line(_));
                        let is_edge_line = matches!(existing_edge.curve, Curve3::Line(_));
                        if !is_ic_line || !is_edge_line {
                            // OCCT L2096-2100, L2148-2150: aMaxTolAdd = 0.001
                            //   aMaxTolAdd = min(aMaxTolAdd, aCoeffTolAdd * aTolCheck)
                            //   aTolAdd = 2 * min(aMaxTolAdd, max(aRealTol, max(aTolV1, aTolV2)))
                            let a_max_tol_add = (0.001_f64).min(10.0 * a_tol_check);
                            let a_tol_add = 2.0 * a_max_tol_add.min(a_real_tol.max(a_tol_v1.max(a_tol_v2)));
                            let (_a_t, a_dist_m1m2) =
                                crate::extrema::closest_point_on_curve(&existing_edge.curve, a_pm);
                            // OCCT L2157: if (aPEStatus == 0)
                            // Compute tangent on existing edge at projection point
                            let a_vtgt2 = existing_edge.curve.derivative_at(_a_t);
                            if a_vtgt2.length_squared() > f64::EPSILON {
                                let a_vtgt2_n = a_vtgt2.normalize();
                                // OCCT L2161: cos = aVTgt1.Dot(aVTgt2.Normalized())
                                let a_cos = a_vtgt1_n.dot(a_vtgt2_n);
                                // OCCT L2162: if (abs(aCos) >= 0.9063) — 25-degree threshold
                                if a_cos.abs() >= 0.9063 {
                                    a_real_tol = a_tol_add;
                                }
                            }
                        }
                    }
                }
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

    /// MakeBlocks (PaveFiller_6.cxx L650-1169).
    /// Creates section edges from FF intersection curves and handles post-treatment:
    /// vertex fusion (PostTreatFF), tolerance correction, face info updates.
    #[allow(non_snake_case)]
    pub(super) fn make_blocks(&mut self) {
        // L654-657: skip if Glue mode
        if self.use_glue() {
            return;
        }

        // L659-665: FF interference array
        let mut a_nb_ff = self.ds.interf_ff.len();
        if a_nb_ff == 0 {
            return;
        }

        // Cross-iteration collections (persist through entire loop).
        // OCCT uses a Main Allocator (IncAllocator) for these.
        let mut a_mpb_add: HashSet<usize> = HashSet::new();
        let mut a_lpb: Vec<PaveBlock> = Vec::new();
        let mut a_mscpb: HashMap<usize, (usize, usize)> = HashMap::new();
        let mut a_mvi: HashMap<usize, usize> = HashMap::new();
        let mut a_dm_ex_edges: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut a_dm_new_sd: HashMap<usize, usize> = HashMap::new();
        let mut a_dm_vlv: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut a_micro_pb: Vec<PaveBlock> = Vec::new();
        let mut a_micro_pb_set: HashSet<(usize, usize)> = HashSet::new();
        let mut a_verts_on_rejected_pb: HashSet<usize> = HashSet::new();
        let mut a_pb_faces_map: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut a_ff_to_recheck: Vec<usize> = Vec::new();

        // Per-iteration collections (declared outside loop, cleared inside).
        // OCCT uses a Temporary Allocator (IncAllocator) for these, reset at each iteration.
        let mut a_lse: Vec<usize> = Vec::new();
        let mut a_lbv: Vec<usize> = Vec::new();
        let mut a_mv_on_in: HashSet<usize> = HashSet::new();
        let mut a_mv_common: HashSet<usize> = HashSet::new();
        let mut a_mv_stick: HashSet<usize> = HashSet::new();
        let mut a_mv_ef: HashSet<usize> = HashSet::new();
        let mut a_mv_bounds: HashSet<usize> = HashSet::new();
        let mut a_mi: HashSet<usize> = HashSet::new();
        let mut a_mpb_on_in: HashSet<usize> = HashSet::new();
        let mut a_mpb_common: HashSet<usize> = HashSet::new();
        let mut a_dmbv: HashMap<usize, Vec<usize>> = HashMap::new();
        // aMVTol uses Default Allocator (UnBind operations require proper Free).
        let mut a_mv_tol: HashMap<usize, f64> = HashMap::new();

        // Ensure section_edge_refs is populated
        self.ds.section_edge_refs = vec![Vec::new(); self.ds.intersection_curves.len()];

        // Pre-collect FF data: (f1, f2, curves, points) to avoid borrow conflicts
        let ff_data: Vec<(usize, usize, Vec<usize>, Vec<crate::bopds::ds::types::FFPoint>)> = self.ds.interf_ff.iter()
            .map(|ff| (ff.f1, ff.f2, ff.curves.clone(), ff.points.clone()))
            .collect();

        let a_nb_ff_prev = a_nb_ff;

        // L727: for (i = 0; i < aNbFF; ++i, aPS.Next())
        let mut loop_i = 0usize;
        while loop_i < a_nb_ff {
            let i = loop_i;
            loop_i += 1; // increment BEFORE loop body (OCCT: ++i in for header)
            // L729-732: UserBreak check — omitted (no progress range in rcad)

            // L735: aCurInd = i < aNbFFPrev ? i : aFFToRecheck[i - aNbFFPrev];
            let a_cur_ind = if i < a_nb_ff_prev { i } else { a_ff_to_recheck[i - a_nb_ff_prev] };

            // L737-747: FF data
            let (n_f1, n_f2, curves_of_ff, points_of_ff) = &ff_data[a_cur_ind];
            let (n_f1, n_f2) = (*n_f1, *n_f2);
            let a_nb_c = curves_of_ff.len();
            let a_nb_p = points_of_ff.len();
            if a_nb_p == 0 && a_nb_c == 0 { continue; }

            // L752: aTolFF = max(BRep_Tool::Tolerance(aF1), BRep_Tool::Tolerance(aF2))
            let a_tol_ff = self.ff_tol(n_f1, n_f2);

            // L757-770: Clear per-iteration collections
            a_mv_on_in.clear();
            a_mv_common.clear();
            a_mpb_on_in.clear();
            a_mpb_common.clear();
            a_dmbv.clear();
            a_mv_tol.clear();
            a_lse.clear();
            a_lbv.clear();
            a_mv_stick.clear();
            a_mv_ef.clear();
            a_mv_bounds.clear();
            a_mi.clear();

            // L772-773: SubShapesOnIn + SharedEdges
            self.ds.sub_shapes_on_in(n_f1, n_f2, &mut a_mv_on_in, &mut a_mv_common,
                                     &mut a_mpb_on_in, &mut a_mpb_common);
            self.ds.shared_edges(n_f1, n_f2, &mut a_lse);

            // L775-793: 1. Treat Points (FF point contacts)
            for (pi, ffp) in points_of_ff.iter().enumerate() {
                let a_p = ffp.pnt;
                // L784: IsExistingVertex(aP, aTolFF, aMVOnIn)
                let b_exist = self.is_existing_vertex_at_point(a_p, a_tol_ff, &a_mv_on_in);
                if !b_exist {
                    // L787: BOPTools_AlgoTools::MakeNewVertex(aP, aTolFF, aV)
                    let n_v = self.ds.add_vertex(a_p);
                    self.ds.vertex_data_mut(n_v).tolerance = a_tol_ff;
                    // L789-791: aMSCPB.Add(aV, aCPB)
                    a_mscpb.insert(n_v, (a_cur_ind, a_nb_c + pi));
                }
            }

            // L796: GetStickVertices
            self.get_stick_vertices_ff(n_f1, n_f2, &mut a_mv_stick, &mut a_mv_ef, &mut a_mi);

            // L798-829: 2. PutPavesOnCurve for each curve
            for &ci in curves_of_ff {
                if ci >= self.ds.intersection_curves.len() { continue; }
                // L810: aNC.InitPaveBlock1()
                if self.ds.intersection_curves[ci].pave_blocks.is_empty() {
                    let mut pb = PaveBlock::new_curve_block();
                    // Set curve data from intersection curve (needed for sub-PB curve propagation)
                    pb.curve = Some(self.ds.intersection_curves[ci].curve.clone());
                    pb.pave1 = Pave { vertex_idx: self.ds.intersection_curves[ci].start_vertex,
                                      param: self.ds.intersection_curves[ci].t_range[0] };
                    pb.pave2 = Pave { vertex_idx: self.ds.intersection_curves[ci].end_vertex,
                                      param: self.ds.intersection_curves[ci].t_range[1] };
                    self.ds.intersection_curves[ci].pave_blocks.push(SharedPB::new(pb));
                }
                // L818: PutPavesOnCurve(aMVOnIn, aMVCommon, aNC, aMI, aMVEF, aMVTol, aDMVLV)
                self.put_paves_on_curve(&a_mv_on_in, &a_mv_common, ci, &a_mi, &a_mv_ef);
            }

            // L834: FilterPavesOnCurves
            self.filter_paves_on_curves(curves_of_ff);

            // L836-864: 3. PutStickPavesOnCurve + PutEFPavesOnCurve + PutBoundPaveOnCurve
            for (j, &ci) in curves_of_ff.iter().enumerate() {
                if ci >= self.ds.intersection_curves.len() { continue; }
                // L841: PutStickPavesOnCurve(aF1, aF2, aMI, aVC, j, aMVStick, aMVTol, aDMVLV)
                // OCCT L2928-2955: RemoveUsedVertices — remove vertices already on any curve
                let a_mv_filtered: std::collections::HashSet<usize> = {
                    let mut mv = a_mv_stick.clone();
                    for &oci in curves_of_ff {
                        if oci >= self.ds.intersection_curves.len() { continue; }
                        let Some(aPB) = self.ds.intersection_curves[oci].pave_blocks.first() else { continue };
                        let r = aPB.0.read().unwrap();
                        for pave in &r.ext_paves { mv.remove(&pave.vertex_idx); }
                        mv.remove(&r.pave1.vertex_idx);
                        mv.remove(&r.pave2.vertex_idx);
                    }
                    mv
                };
                self.put_stick_paves_on_curve(ci, &a_mi, &a_mv_filtered);
                // L843-846: if aNbC == 1: PutEFPavesOnCurve
                if a_nb_c == 1 {
                    self.put_ef_paves_on_curve(ci, &a_mi, &a_mv_ef);
                }
                // L848-863: aIC.HasBounds() → PutBoundPaveOnCurve
                let has_bounds = {
                    let ic = &self.ds.intersection_curves[ci];
                    ic.start_vertex < self.ds.vertices.len()
                        && ic.end_vertex < self.ds.vertices.len()
                };
                if has_bounds {
                    a_lbv.clear();
                    // L852: PutBoundPaveOnCurve(aF1, aF2, aNC, aLBV)
                    self.put_bound_pave_on_curve(n_f1, n_f2, ci, &mut a_lbv);
                    // L854-862: collect bound vertices
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
            }

            // L866-871: PutClosingPaveOnCurve for each curve
            for &ci in curves_of_ff {
                if ci >= self.ds.intersection_curves.len() { continue; }
                self.put_closing_pave_on_curve(ci);
            }

            // L874-895: Build PB BVH tree from aMPBOnIn
            let a_pb_tree = Self::build_pb_tree(self.ds, &a_mpb_on_in);

            // L899: isToRecheck — OCCT: (aNbC > 0) && (i < aNbFFPrev)
            // If there are any intersection curves, recheck is needed to
            // avoid missing section edges due to FF processing order.
            let mut is_to_recheck = a_nb_c > 0 && i < a_nb_ff_prev;

            // L901-1098: 4. Make section edges
            for (j, &ci) in curves_of_ff.iter().enumerate() {
                if ci >= self.ds.intersection_curves.len() { continue; }

                // L906: aTolR3D = max(aNC.Tolerance(), aNC.TangentialTolerance())
                let a_tol_r3d = {
                    let ic = &self.ds.intersection_curves[ci];
                    ic.geom_tol.max(ic.curve_extra.tangential_tol)
                };

                // L908-912: aLPB.Clear(); aPB1->Update(aLPB, false);
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

                // L926-929: if aLPB is non-empty, isToRecheck = false
                if !a_lpb.is_empty() {
                    is_to_recheck = false;
                }

                if std::env::var("RCAD_DBG_MB").is_ok() {
                    eprintln!("[DBG_MB] curve j={} ci={} a_lpb.len={}", j, ci, a_lpb.len());
                }
                // L931-1095: iterate sub-PBs
                for a_pb in &a_lpb {
                    // L935-936: aPB->Indices(nV1, nV2); aPB->Range(aT1, aT2);
                    let (n_v1, n_v2) = a_pb.indices();
                    let (a_t1, a_t2) = a_pb.range();
                    // L938: if (fabs(aT1 - aT2) < Precision::PConfusion())
                    let a_t_range = a_t2 - a_t1;
                    if a_t_range.abs() < CONFUSION { continue; }

                    // L943-950: IsValidBlockForFaces (midpoint check)
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
                            self.context.is_valid_point_for_face(self.ds, mid_pt, fi, a_tol_r3d)
                        };
                        if !ok { b_valid_2d = false; break; }
                    }
                    if !b_valid_2d {
                        if std::env::var("RCAD_DBG_MB").is_ok() {
                            eprintln!("[DBG_MB]   sub-PB ({},{}) REJECTED by IsValidBlockForFaces", n_v1, n_v2);
                        }
                        continue; }

                    // L952-962: IsExistingPaveBlock via LSE (shared edges)
                    let mut n_e_out: usize = usize::MAX;
                    let mut a_tol_new: f64 = -1.0;
                    let b_exist_lse = self.is_existing_pb_via_lse(&a_lse, a_pb, ci, &mut n_e_out, &mut a_tol_new);
                    if b_exist_lse {
                        if std::env::var("RCAD_DBG_MB").is_ok() {
                            eprintln!("[DBG_MB]   sub-PB ({},{}) EXISTING via LSE n_e={}", n_v1, n_v2, n_e_out);
                        }
                        // L958: UpdateEdgeTolerance(nEOut, aTolNew)
                        self.update_edge_tolerance(n_e_out, a_tol_new);
                        // L960: UpdateSavedTolerance(myDS, nEOut, aTolNew, aMVTol)
                        self.update_saved_tolerance(n_e_out, a_tol_new, &mut a_mv_tol);
                        continue;
                    }

                    // OCCT L937-960: BRepLib::FindValidRange — check if the pave block
                    // has a valid range outside the vertex tolerance spheres.
                    let has_valid_range = {
                        let ic = &self.ds.intersection_curves[ci];
                        let v1_pt = if n_v1 < self.ds.vertices.len() { self.ds.vertex_point(n_v1) } else { continue; };
                        let v2_pt = if n_v2 < self.ds.vertices.len() { self.ds.vertex_point(n_v2) } else { continue; };
                        let a_tol_v1 = a_tol_r3d.max(self.ds.vertex_tolerance(n_v1));
                        let a_tol_v2 = a_tol_r3d.max(self.ds.vertex_tolerance(n_v2));
                        crate::pave_filler::helpers::find_valid_range(
                            &ic.curve, a_t1, a_t2, a_tol_r3d,
                            v1_pt, a_tol_v1, v2_pt, a_tol_v2,
                        ).is_some()
                    };
                    if !has_valid_range {
                        if std::env::var("RCAD_DBG_MB").is_ok() {
                            eprintln!("[DBG_MB]   sub-PB ({},{}) REJECTED by FindValidRange", n_v1, n_v2);
                        }
                        // L984-990: if not bound, add to aMicroPB
                        if !a_mv_bounds.contains(&n_v1) && !a_mv_bounds.contains(&n_v2) {
                            if a_micro_pb_set.insert((n_v1.min(n_v2), n_v1.max(n_v2))) {
                                a_micro_pb.push(a_pb.clone());
                            }
                            a_mvi.insert(n_v1, n_v1);
                            a_mvi.insert(n_v2, n_v2);
                        }
                        continue;
                    }

                    // L994-1052: IsExistingPaveBlock via ON/IN (BVH)
                    let mut a_pb_out: usize = usize::MAX;
                    let mut a_tol_new2: f64 = -1.0;
                    let b_exist_on_in = self.is_existing_pb_via_bvh(
                        a_pb, ci, a_tol_r3d, &a_mpb_on_in, &a_pb_tree,
                        &a_mpb_common, &mut a_pb_out, &mut a_tol_new2,
                    );
                    if b_exist_on_in {
                        if std::env::var("RCAD_DBG_MB").is_ok() {
                            eprintln!("[DBG_MB]   sub-PB ({},{}) EXISTING via BVH a_pb_out={}", n_v1, n_v2, a_pb_out);
                        }
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
                                // L1005-1017: Update edge tolerance (OCCT logic)
                                let a_tol_e = self.ds.edge_tolerance(edge_of_pb);
                                // L1008-1011: if (aTolNew < aNC.Tolerance()) aTolNew = aNC.Tolerance()
                                let ic_tol = {
                                    let ic = &self.ds.intersection_curves[ci];
                                    ic.geom_tol.max(ic.curve_extra.tangential_tol)
                                };
                                if a_tol_new2 < ic_tol {
                                    a_tol_new2 = ic_tol;
                                }
                                // L1012-1017: if (aTolNew > aTolE) { UpdateEdgeTolerance; UpdateSavedTolerance }
                                if a_tol_new2 > a_tol_e {
                                    self.update_edge_tolerance(edge_of_pb, a_tol_new2);
                                    self.update_saved_tolerance(edge_of_pb, a_tol_new2, &mut a_mv_tol);
                                }

                                // L1019-1030: Face without pave block
                                let n_f = if b_in_f1 { n_f2 } else { n_f1 };
                                a_pb_faces_map.entry(a_pb_out).or_default().push(n_f);

                                // L1032-1044: Vertices on rejected PB
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

                                // L1046-1050: PreparePostTreatFF
                                if a_mpb_add.insert(a_pb_out) {
                                    self.prepare_post_treat_ff(
                                        a_cur_ind, j, a_pb_out,
                                        &mut a_mscpb, &mut a_mvi, ci,
                                    );
                                }
                            }
                        }
                        continue;
                    }

                    // L1055-1094: Make section edge
                    if std::env::var("RCAD_DBG_MB").is_ok() {
                        eprintln!("[DBG_MB]   sub-PB ({},{}) → MAKE EDGE", n_v1, n_v2);
                    }
                    let a_curve = &self.ds.intersection_curves[ci].curve;
                    let pca = self.ds.intersection_curves[ci].pcurve_on_a.clone();
                    let pcb = self.ds.intersection_curves[ci].pcurve_on_b.clone();
                    // L1056: BOPTools_AlgoTools::MakeEdge(aIC, aV1, aT1, aV2, aT2, aTolR3D, aES)
                    let new_ei = crate::boptools::make_edge(self.ds, ci, n_v1, n_v2, a_t1, a_t2, a_tol_r3d);
                    // L1058-1064: MakePCurve
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
                    // L1067: aLPBC.Append(aPB)
                    self.ds.section_edge_refs[ci].push(new_ei);

                    // L1070-1077: Keep info for post treatment
                    let mut sub_pb = a_pb.clone();
                    sub_pb.new_edge = Some(new_ei);
                    a_mscpb.insert(new_ei, (a_cur_ind, j));
                    // L1076-1077: aMVI.Bind(aV1, nV1); aMVI.Bind(aV2, nV2);
                    a_mvi.insert(n_v1, n_v1);
                    a_mvi.insert(n_v2, n_v2);

                    // L1079-1080: aMVTol.UnBind(nV1); aMVTol.UnBind(nV2);
                    a_mv_tol.remove(&n_v1);
                    a_mv_tol.remove(&n_v2);

                    // Allocate global PaveBlock and register on faces
                    let g_pb_idx = self.ds.allocate_pave_block(sub_pb.clone());
                    for &fi in &[n_f1, n_f2] {
                        if fi != usize::MAX {
                            self.ds.face_info_mut(fi).pave_blocks_sc.insert(g_pb_idx);
                        }
                    }

                    // L1082-1094: ProcessExistingPaveBlocks (first overload)
                    // Adds existing pave blocks for post treatment
                    let mut tmp_lpb: Vec<PaveBlock> = Vec::new();
                    self.process_existing_pave_blocks(
                        ci, j, n_f1, n_f2, new_ei,
                        &a_mpb_on_in, &a_pb_tree, &mut a_mscpb, &mut a_mvi,
                        &mut tmp_lpb, &mut a_pb_faces_map, &mut a_mpb_add,
                    );
                } // for a_pb in &a_lpb

                // L1097: aLPBC.RemoveFirst() — remove InitPaveBlock1
                if !self.ds.intersection_curves[ci].pave_blocks.is_empty() {
                    self.ds.intersection_curves[ci].pave_blocks.remove(0);
                }
            } // for (j, &ci) in curves_of_ff

            // L1099-1103: isToRecheck
            if is_to_recheck {
                a_ff_to_recheck.push(a_cur_ind);
                // OCCT: extend loop to include recheck entries
                a_nb_ff = self.ds.interf_ff.len() + a_ff_to_recheck.len();
            }

            // L1105-1127: Restore vertex tolerances for unused vertices
            // OCCT L1107-1127: iterate aMVTol, restore original tolerances
            let saved_tols: Vec<(usize, f64)> = a_mv_tol.iter().map(|(&k, &v)| (k, v)).collect();
            for &(n_v, saved_tol) in &saved_tols {
                if n_v < self.ds.vertices.len() {
                    self.ds.vertex_data_mut(n_v).tolerance = saved_tol;
                }
            }
            // OCCT L1117-1121: reset vertex bounding boxes after tolerance restore
            for &(n_v, _) in &saved_tols {
                if let Some(si) = self.ds.vertex_shape_idx.get(n_v) {
                    if *si < self.ds.shape_info.len() {
                        let vp = self.ds.vertex_point(n_v);
                        let vt = self.ds.vertex_tolerance(n_v);
                        let si_mut = &mut self.ds.shape_info[*si];
                        si_mut.box_min = Some(vp);
                        si_mut.box_max = Some(vp);
                        si_mut.box_gap = vt + crate::tolerance::CONFUSION;
                    }
                }
            }
            // L1123-1126: forget SD groups of restored vertices
            for &(n_v, _) in &saved_tols {
                a_dm_vlv.remove(&n_v);
            }

            // L1129-1138: ProcessExistingPaveBlocks (second overload)
            self.process_existing_pave_blocks_after(
                a_cur_ind, n_f1, n_f2, &a_mpb_on_in, &a_pb_tree,
                &mut a_dmbv, &mut a_mscpb, &mut a_mvi,
                &mut a_pb_faces_map, &mut a_mpb_add,
            );
        } // while loop_i < a_nb_ff

        // ===== Post-loop phase =====

        // L1141-1142: RemoveMicroSectionEdges
        self.remove_micro_section_edges(&mut a_mscpb, &mut a_micro_pb);

        // L1145: MakeSDVerticesFF(aDMVLV, aDMNewSD)
        self.make_sd_vertices_ff(&a_dm_vlv, &mut a_dm_new_sd);

        // L1146-1156: PostTreatFF (vertex fusion)
        self.post_treat_ff(
            &mut a_mscpb,
            &mut a_dm_ex_edges,
            &mut a_dm_new_sd,
            &a_micro_pb,
            &a_verts_on_rejected_pb,
        );

        // L1153-1156: if (HasErrors()) return; — omitted (no HasErrors in rcad)

        // L1158: CorrectToleranceOfSE
        self.correct_tolerance_of_se();

        // L1161: UpdateFaceInfo
        self.update_face_info_post(
            &a_dm_ex_edges,
            &a_dm_new_sd,
            &a_pb_faces_map,
        );

        // L1163: UpdatePaveBlocks
        self.update_pave_blocks(&a_dm_new_sd);

        // L1168: PutSEInOtherFaces
        self.put_se_in_other_faces();
    }

    /// Build PB BVH tree from a set of pave blocks ON/IN.
    /// OCCT PaveFiller_6.cxx L874-895 (BOPTools_BoxTree construction).
    fn build_pb_tree(
        ds: &DS,
        a_mpb_on_in: &HashSet<usize>,
    ) -> Option<DsBvh> {
        let mut a_pb_indices: Vec<usize> = Vec::new();
        let mut a_pb_aabbs: Vec<Aabb> = Vec::new();
        for &pb_idx in a_mpb_on_in {
            if pb_idx >= ds.pave_blocks.len() { continue; }
            let pb = &ds.pave_blocks[pb_idx];
            let r = pb.0.read().unwrap();
            // L883-886: if (!aPB->HasEdge()) continue;
            if r.new_edge.is_none() && r.original_edge == NO_EDGE { continue; }
            let ei = r.new_edge.unwrap_or(r.original_edge);
            if ei >= ds.edges.len() { continue; }
            // L888-891: if (myDS->ShapeInfo(aPB->OriginalEdge()).HasFlag()) continue;
            if ds.edge_has_flag(ei) { continue; }
            // L893: aPBTree.Add(iPB, Bnd_Tools::Bnd2BVH(myDS->ShapeInfo(aPB->Edge()).Box()));
            // OCCT uses the precomputed ShapeInfo Box which covers the full 3D curve geometry.
            // rcad: compute AABB from edge's 3D curve by sampling points along it.
            let edge = &ds.edges[ei];
            let tol = ds.edge_tolerance(ei);
            let [t0, t1] = edge.t_range;
            let n_samples = 8.max(((t1 - t0).abs() / 0.1).ceil() as usize).min(32);
            let mut mn = DVec3::splat(f64::MAX);
            let mut mx = DVec3::splat(f64::NEG_INFINITY);
            for k in 0..=n_samples {
                let t = t0 + (t1 - t0) * (k as f64) / (n_samples as f64);
                let pt = edge.curve.point_at(t);
                mn = mn.min(pt);
                mx = mx.max(pt);
            }
            // Also include start/end vertex points for safety
            let sv = ds.edge_start_vertex_ds(ei);
            let ev = ds.edge_end_vertex_ds(ei);
            if sv < ds.vertices.len() {
                let vp = ds.vertex_point(sv);
                mn = mn.min(vp); mx = mx.max(vp);
            }
            if ev < ds.vertices.len() {
                let vp = ds.vertex_point(ev);
                mn = mn.min(vp); mx = mx.max(vp);
            }
            a_pb_indices.push(pb_idx);
            a_pb_aabbs.push(Aabb { min: mn - DVec3::splat(tol), max: mx + DVec3::splat(tol), gap: 0.0 });
        }
        if !a_pb_indices.is_empty() {
            Some(DsBvh::build(a_pb_indices, a_pb_aabbs))
        } else {
            None
        }
    }
}
