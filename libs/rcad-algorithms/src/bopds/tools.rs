///  ?BOPAlgo_Tools  ?CommonBlock merging and tolerance computation.
///
/// OCCT references:
/// - `BOPAlgo_Tools::PerformCommonBlocks` (BOPAlgo_Tools.cxx L107-243)
/// - `BOPAlgo_Tools::ComputeToleranceOfCB` (BOPAlgo_Tools.cxx L248-356)
///
/// These functions process CommonBlocks after PaveBlocks have been created
/// by the PaveFiller.  CommonBlocks group PaveBlocks from different edges
/// (and different faces) that lie on the same geometry.
use std::collections::HashMap;

use super::common_block::CommonBlock;
use super::ds::DS;
use crate::tolerance::*;

use crate::bopds::pave::NO_EDGE;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};

//  € € PerformCommonBlocks  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

///  ?`BOPAlgo_Tools::PerformCommonBlocks` (overload 1  ?PB B
/// connection map).
///
/// Scans all edges in the DS, groups PaveBlocks that share the same
/// (start_vertex, end_vertex) pair (geometrically coincident), and creates
/// CommonBlocks for each group with at least two participating faces.
///
/// Additionally populates the global `ds.pave_blocks` array from the
/// edge-local PaveBlocks so that CommonBlock indices are valid for later
/// consumption by `FillSameDomainFaces` in the Builder.
///
/// Calling this after all paves have been placed (post-`build_split_edges`)
/// ensures the vertex indices reflect the final shared topology.
pub fn perform_common_blocks(ds: &mut DS) {
    //  € € Phase 1: Populate global PaveBlock array  € € € € € € € € € € € € € € € € € € € € € € € €
    // Map (edge_idx, local_pb_idx)  ?global PaveBlock index.
    let mut edge_local_to_global: HashMap<(usize, usize), usize> = HashMap::new();
    // Reverse map: global_pb_idx  ?(edge_idx, local_pb_idx).
    let mut global_to_edge_local: HashMap<usize, (usize, usize)> = HashMap::new();

    for (ei, edge) in ds.edges.iter().enumerate() {
        for local_i in 0..edge.pave_blocks.len() {
            let global_i = ds.pave_blocks.len();
            ds.pave_blocks.push(edge.pave_blocks[local_i].clone());
            edge_local_to_global.insert((ei, local_i), global_i);
            global_to_edge_local.insert(global_i, (ei, local_i));
        }
    }

    //  € € Phase 2: Group coincident PaveBlocks  € € € € € € € € € € € € € € € € € € € € € € € € € € € €
    // Two PaveBlocks from different edges are geometrically coincident when
    // they share the same ordered (v_min, v_max) vertex pair.  After OCCT's
    // PutPaveOnCurve / PutBoundPaveOnCurve all intersection vertices are
    // merged into the DS pool, so the same 3D position on coincident edges
    // maps to the same DS vertex index.
    //
    // Key: (v_min, v_max) with v_min  ?v_max.
    // Value: Vec<(global_pb_idx, face_idx)>  ?the PaveBlocks sharing this
    // vertex pair, along with the face they belong to.
    let mut vertex_groups: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();

    for (ei, edge) in ds.edges.iter().enumerate() {
        for local_i in 0..edge.pave_blocks.len() {
            let pb = &edge.pave_blocks[local_i];
            // OCCT: skip section edges (no original edge)
            if pb.0.read().unwrap().original_edge == NO_EDGE {
                continue;
            }
            let v1 = pb.0.read().unwrap().pave1.vertex_idx;
            let v2 = pb.0.read().unwrap().pave2.vertex_idx;
            let key = if v1 <= v2 { (v1, v2) } else { (v2, v1) };

            let global_pb = match edge_local_to_global.get(&(ei, local_i)) {
                Some(&idx) => idx,
                None => continue,
            };

            // Derive face indices from edge's face_reps (OCCT records this at EF intersection).
            let face_indices: Vec<usize> = edge.face_reps.iter().map(|r| r.face_idx).collect();

            for &fi in &face_indices {
                vertex_groups.entry(key).or_default().push((global_pb, fi));
            }
        }
    }

    //  € € Phase 3: Create CommonBlocks  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
    for ((_v_min, _v_max), entries) in vertex_groups.iter() {
        if entries.len() < 2 {
            continue;
        }

        // Require at least two distinct faces.
        let mut unique_faces: Vec<usize> = entries.iter().map(|&(_, fi)| fi).collect();
        unique_faces.sort();
        unique_faces.dedup();
        if unique_faces.len() < 2 {
            continue;
        }

        let a_original_edge = entries
            .first()
            .and_then(|(gp, _)| ds.pave_blocks.get(*gp))
            .map_or(usize::MAX, |spb| spb.0.read().unwrap().original_edge);
        let mut cb = CommonBlock::new();
        // OCCT: set the real edge index for sort_by_edge
        cb.set_edge(a_original_edge);
        for &(global_pb, fi) in entries {
            cb.add_pave_block(global_pb, fi);
        }
        cb.set_faces(unique_faces);
        // OCCT 1:1: ensure the real-edge PB comes first
        cb.sort_by_edge(|pbi| {
            ds.pave_blocks
                .get(pbi)
                .map_or(false, |spb| a_original_edge == a_original_edge)
        });

        // Compute tolerance for this CommonBlock.
        let cb_idx = ds.common_blocks.len();
        ds.common_blocks.push(cb);
        let tol = compute_tolerance_of_cb(ds, cb_idx);
        ds.common_blocks[cb_idx].set_tolerance(tol);

        // Mark local PaveBlocks as belonging to this CommonBlock.
        //  ?BOPDS_PaveBlock::myCommonBlock (L103-108 in ds.cxx).
        for &(global_pb, _) in entries {
            if let Some(&(ei, local_i)) = global_to_edge_local.get(&global_pb) {
                if let Some(local_pb) = ds
                    .edges
                    .get_mut(ei)
                    .and_then(|e| e.pave_blocks.get_mut(local_i))
                {
                    local_pb.0.write().unwrap().common_block_idx = Some(cb_idx);
                }
            }
        }
    }
}

//  € € ComputeToleranceOfCB  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

///  ?`BOPAlgo_Tools::ComputeToleranceOfCB`.
///
/// Computes the maximum tolerance for a CommonBlock by sampling points along
/// the reference PaveBlock's curve and measuring the deviation to:
/// - Each other PaveBlock's curve (via numeric projection)
/// - Each face's surface (via closest-point-on-surface)
///
/// The result starts from `geom_tol` of the reference edge and accumulates
/// the projection distance of each sample.
///
/// OCCT reference: BOPAlgo_Tools.cxx L248-356.
pub fn compute_tolerance_of_cb(ds: &DS, cb_idx: usize) -> f64 {
    let cb = &ds.common_blocks[cb_idx];
    let pb_entries = cb.pave_blocks();
    if pb_entries.is_empty() {
        return TOLERANCE_ABS;
    }

    // Reference PaveBlock (the first one).
    let (ref_pb_idx, _) = pb_entries[0];
    let ref_pb = &ds.pave_blocks[ref_pb_idx];
    let ref_edge = &ds.edges[ref_pb.0.read().unwrap().original_edge];
    let ref_curve = &ref_edge.curve;
    let t1 = ref_pb
        .0
        .read()
        .unwrap()
        .pave1
        .param
        .min(ref_pb.0.read().unwrap().pave2.param);
    let t2 = ref_pb
        .0
        .read()
        .unwrap()
        .pave1
        .param
        .max(ref_pb.0.read().unwrap().pave2.param);

    // Start tolerance from the reference edge's model tolerance (OCCT aTolMax = BRep_Tool::Tolerance(aEOr)).
    let mut tol_max = ref_edge.geom_tol;

    // OCCT uses 11 interior sample points: aDt = (aT2 - aT1) / (aNbPnt + 1).
    let n_samples = 11usize;
    if (t2 - t1).abs() < TOLERANCE_ABS * 0.01 {
        return tol_max;
    }
    let dt = (t2 - t1) / (n_samples as f64 + 1.0);

    //  € € 1. Other PaveBlocks: sample reference curve, project onto other curves  € €
    if pb_entries.len() > 1 {
        for &(pb_idx, _face_idx) in &pb_entries[1..] {
            let pb = &ds.pave_blocks[pb_idx];
            let other_edge = &ds.edges[pb.0.read().unwrap().original_edge];
            let other_curve = &other_edge.curve;
            let other_t_range = [
                pb.0.read()
                    .unwrap()
                    .pave1
                    .param
                    .min(pb.0.read().unwrap().pave2.param),
                pb.0.read()
                    .unwrap()
                    .pave1
                    .param
                    .max(pb.0.read().unwrap().pave2.param),
            ];
            let other_tol = other_edge.geom_tol;

            let mut t = t1;
            for _ in 0..n_samples {
                t += dt;
                let pt = ref_curve.point_at(t);
                let dist = min_distance_to_curve(pt, other_curve, other_t_range);
                let tol_new = other_tol + dist;
                if tol_new > tol_max {
                    tol_max = tol_new;
                }
            }
        }
    }

    //  € € 2. Faces: project sample points onto each face's surface  € € € € € € € € € € € € €
    {
        let faces_list = cb.faces();
        if !faces_list.is_empty() {
            for &fi in faces_list {
                let face = &ds.faces[fi];
                let surf = &face.surface;
                let face_tol = face.geom_tol;

                let mut t = t1;
                for _ in 0..n_samples {
                    t += dt;
                    let pt = ref_curve.point_at(t);
                    let (_uv, proj_pt) = crate::extrema::closest_point_on_surface(surf, pt);
                    let dist = proj_pt.distance(pt);
                    let tol_new = face_tol + dist;
                    if tol_new > tol_max {
                        tol_max = tol_new;
                    }
                }
            }
        }
    }

    tol_max
}

//  € € Helpers  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €

/// Compute the minimum distance from a point to a curve by sampling.
///
/// OCCT uses `GeomAPI_ProjectPointOnCurve` which internally uses a
/// Newton-based minimiser.  This is a simplified numeric fallback:
/// 64 uniform samples over the curve's parameter range, followed by
/// a coarse Newton refinement for the nearest candidate.
fn min_distance_to_curve(pt: DVec3, curve: &Curve3, t_range: [f64; 2]) -> f64 {
    let t0 = t_range[0];
    let t1 = t_range[1];
    let span = t1 - t0;
    if span.abs() < TOLERANCE_ABS * 0.01 {
        // Degenerate range: just measure distance to midpoint.
        let mid = curve.point_at(t0);
        return mid.distance(pt);
    }

    // 64 uniform samples to find the nearest candidate.
    let n_coarse = 64usize;
    let mut best_t = t0;
    let mut best_dist_sq = f64::MAX;
    for i in 0..=n_coarse {
        let t = t0 + span * (i as f64 / n_coarse as f64);
        let cp = curve.point_at(t);
        let d_sq = cp.distance_squared(pt);
        if d_sq < best_dist_sq {
            best_dist_sq = d_sq;
            best_t = t;
        }
    }

    // Coarse Newton-style refinement: sample closer around the best candidate.
    let mut search_span = span / n_coarse as f64;
    for _ in 0..4 {
        search_span *= 0.5;
        if search_span < TOLERANCE_ABS {
            break;
        }
        for offset in [-1.0, 0.0, 1.0] {
            let t_candidate = (best_t + offset * search_span).clamp(t0, t1);
            let cp = curve.point_at(t_candidate);
            let d_sq = cp.distance_squared(pt);
            if d_sq < best_dist_sq {
                best_dist_sq = d_sq;
                best_t = t_candidate;
            }
        }
    }

    best_dist_sq.sqrt()
}

//  € € Tests  € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € € €
