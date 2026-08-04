// BOPAlgo_Tools ?CommonBlock merging and tolerance computation.
//
// OCCT: BOPAlgo_Tools.cxx L107-356
// Operates on BOPDS_DS to group geometrically coincident PaveBlocks.
use std::collections::HashMap;
use glam::DVec3;
use rcad_kernel::geom::{Curve3, CurveEval};
use crate::bop::ds::common_block::CommonBlock;
use crate::bop::ds::pave::{NO_EDGE, Pave, PaveBlock, SharedPB};
use crate::bop::ds::DS;

/// BOPAlgo_Tools::PerformCommonBlocks ?group PBs with same (v1,v2) into CommonBlocks.
pub fn perform_common_blocks(ds: &mut DS) {
    // Collect all PBs from the per-edge pool.
    let mut all_pbs: Vec<(usize, usize, SharedPB)> = Vec::new();
    for ei in 0..ds.edge_count() {
        let pbs = ds.pave_blocks(ei);
        for local_i in 0..pbs.len() {
            all_pbs.push((ei, local_i, pbs[local_i].clone()));
        }
    }
    // Build map: (v1, v2) -> Vec<(global_idx, ei, local_i, SharedPB)>
    let mut groups: HashMap<(usize, usize), Vec<(usize, usize, usize, SharedPB)>> = HashMap::new();
    for &(ei, local_i, ref pb) in &all_pbs {
        let r = pb.0.read().unwrap();
        let v1 = r.pave1.vertex_idx;
        let v2 = r.pave2.vertex_idx;
        let key = if v1 <= v2 { (v1, v2) } else { (v2, v1) };
        groups.entry(key).or_default().push((all_pbs.len(), ei, local_i, pb.clone()));
    }
    // Create CommonBlocks for groups that span multiple edges.
    for (_key, entries) in groups.iter() {
        if entries.len() < 2 { continue; }
        let mut cb = CommonBlock::new();
        let mut faces_seen: Vec<usize> = Vec::new();
        for &(_gi, ei, _li, ref pb) in entries {
            let oe = pb.0.read().unwrap().original_edge;
            if oe < ds.nb_shapes() {
                // Find faces that reference this edge
                for fi in 0..ds.face_count() {
                    let fi_shape = ds.face_shape_idx(fi);
                    if fi_shape < ds.nb_shapes() {
                        let si = &ds.shapes[fi_shape];
                        if si.sub_shapes.contains(&oe) {
                            if !faces_seen.contains(&fi) {
                                faces_seen.push(fi);
                            }
                        }
                    }
                }
            }
        }
        // Add PBs and faces to CommonBlock.
        for &(_gi, ei, _li, ref pb) in entries {
            let oe = pb.0.read().unwrap().original_edge;
            for &fi in &faces_seen {
                if oe < ds.nb_shapes() {
                    cb.add_pave_block(pb.clone(), fi);
                    cb.add_face(fi);
                }
            }
        }
        if !cb.pave_blocks().is_empty() {
            let tol = compute_tolerance_of_cb(ds, cb.pave_blocks());
            cb.set_tolerance(tol);
            ds.common_blocks.push(cb);
        }
    }
}
/// BOPAlgo_Tools::ComputeToleranceOfCB ?compute combined tolerance for a CommonBlock.
fn compute_tolerance_of_cb(ds: &DS, pb_list: &[(SharedPB, usize)]) -> f64 {
    let mut tol_max = rcad_kernel::CONFUSION;
    for (pb, _) in pb_list {
        let pb_idx = pb.0.read().unwrap().original_edge;
        if pb_idx >= ds.common_blocks.len() { continue; }
        // Use PaveBlock data to estimate tolerance from edge geometry.
        if let Some(curve) = ds.edge_curve(pb_idx) {
            let is_line = matches!(curve, Curve3::Line(_));
            let base = if is_line { rcad_kernel::CONFUSION } else { rcad_kernel::CONFUSION * 10.0 };
            tol_max = tol_max.max(base);
        }
    }
    // Sample extra tolerance from the first face's surface (if available).
    let _ = ds; // placeholder for future surface-based tolerance computation
    tol_max
}