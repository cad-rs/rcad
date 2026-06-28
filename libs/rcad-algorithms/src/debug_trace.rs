/// Debug tracing for Boolean alignment debugging.
///
/// Usage: enable by setting env var `RCAD_TRACE=1`.
/// All trace output goes to stderr (eprintln!) so it is visible even when
/// test harness captures stdout.
use std::collections::BTreeSet;

/// Returns true if RCAD_TRACE is set (any non-empty value).
pub fn is_enabled() -> bool {
    std::env::var("RCAD_TRACE").is_ok()
}

/// Conditional eprintln! — only prints when RCAD_TRACE=1.
#[macro_export]
macro_rules! trace {
    ($($arg:tt)*) => {
        if $crate::debug_trace::is_enabled() {
            eprintln!($($arg)*);
        }
    }
}

/// Count unique DS vertex indices used by result edges in the final BRep.
pub fn count_used_vertices(result: &rcad_kernel::BRep) -> usize {
    let mut used = BTreeSet::new();
    for s in &result.solids {
        for sh in &s.shells {
            for f in &sh.faces {
                for we in &f.outer_wire.edges {
                    if let Some(e) = result.edges.get(we.idx) {
                        used.insert(e.start);
                        used.insert(e.end);
                    }
                }
                for w in &f.inner_wires {
                    for we in &w.edges {
                        if let Some(e) = result.edges.get(we.idx) {
                            used.insert(e.start);
                            used.insert(e.end);
                        }
                    }
                }
            }
        }
    }
    for fiv in &result.geom.face_internal_vertices {
        used.extend(fiv);
    }
    used.len()
}

/// Print result topology summary.
pub fn dump_result_topo(result: &rcad_kernel::BRep, label: &str) {
    let n_verts_used = count_used_vertices(result);
    let n_verts_total = result.vertices.len();
    let n_edges = result.edges.len();
    let n_faces: usize = result.solids.iter()
        .flat_map(|s| &s.shells)
        .map(|sh| sh.faces.len())
        .sum();
    let n_shells: usize = result.solids.iter()
        .flat_map(|s| &s.shells)
        .count();
    let n_solids = result.solids.len();
    eprintln!("[TRACE] {}: verts_total={} verts_used={} edges={} faces={} shells={} solids={}",
        label, n_verts_total, n_verts_used, n_edges, n_faces, n_shells, n_solids);
}

/// Print DS vertex/edge/face/interference counts.
pub fn dump_ds(ds: &crate::bopds::ds::DS, label: &str) {
    eprintln!("[TRACE] {}: vertices={} edges={} faces={} intersect_curves={} interferences={}",
        label, ds.vertices.len(), ds.edges.len(), ds.faces.len(),
        ds.intersection_curves.len(), ds.interferences.len());
}
