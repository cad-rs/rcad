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
    for ts in &result.tshapes {
        if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let rcad_kernel::topods::TShape::Shell(shd) = &*result.tshapes[shell_sr.index] {
                    for face_sr in &shd.faces {
                        if let rcad_kernel::topods::TShape::Face(fd) = &*result.tshapes[face_sr.index] {
                            // outer wire
                            if let rcad_kernel::topods::TShape::Wire(wd) = &*result.tshapes[fd.outer_wire.index] {
                                for e_sr in &wd.edges {
                                    if let rcad_kernel::topods::TShape::Edge(ed) = &*result.tshapes[e_sr.index] {
                                        used.insert(ed.first.index);
                                        used.insert(ed.last.index);
                                    }
                                }
                            }
                            // inner wires
                            for iw_sr in &fd.inner_wires {
                                if let rcad_kernel::topods::TShape::Wire(wd) = &*result.tshapes[iw_sr.index] {
                                    for e_sr in &wd.edges {
                                        if let rcad_kernel::topods::TShape::Edge(ed) = &*result.tshapes[e_sr.index] {
                                            used.insert(ed.first.index);
                                            used.insert(ed.last.index);
                                        }
                                    }
                                }
                            }
                            // internal vertices
                            for iv_sr in &fd.internal_vertices {
                                used.insert(iv_sr.index);
                            }
                        }
                    }
                }
            }
        }
    }
    used.len()
}

/// Print result topology summary.
pub fn dump_result_topo(result: &rcad_kernel::BRep, label: &str) {
    let n_verts_used = count_used_vertices(result);
    let n_verts_total = result.vertex_count();
    let n_edges = result.edge_count();
    let mut n_faces = 0usize;
    let mut n_shells = 0usize;
    for ts in &result.tshapes {
        if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
            for shell_sr in &sd.shells {
                if let rcad_kernel::topods::TShape::Shell(shd) = &*result.tshapes[shell_sr.index] {
                    n_shells += 1;
                    n_faces += shd.faces.len();
                }
            }
        }
    }
    let n_solids = result.solid_count();
    eprintln!("[TRACE] {}: verts_total={} verts_used={} edges={} faces={} shells={} solids={}",
        label, n_verts_total, n_verts_used, n_edges, n_faces, n_shells, n_solids);
}

/// Print DS vertex/edge/face/interference counts.
pub fn dump_ds(ds: &crate::bopds::ds::DS, label: &str) {
    let n_interf = ds.interf_vv.len() + ds.interf_ve.len() + ds.interf_vf.len()
        + ds.interf_ee.len() + ds.interf_ef.len() + ds.interf_ff.len();
    eprintln!("[TRACE] {}: vertices={} edges={} faces={} intersect_curves={} interferences={}",
        label, ds.vertices.len(), ds.edges.len(), ds.faces.len(),
        ds.intersection_curves.len(), n_interf);
}

/// Count unique TShape::Edge indices referenced by all Face TShapes in a topods BRep.
pub fn count_edges_in_faces(t: &rcad_kernel::topods::BRep) -> usize {
    let mut used = std::collections::BTreeSet::new();
    for ts in &t.tshapes {
        if let rcad_kernel::topods::TShape::Face(fd) = &**ts {
            for wire_sr in std::iter::once(&fd.outer_wire).chain(&fd.inner_wires) {
                let wd = t.wire(*wire_sr);
                for e_sr in &wd.edges {
                    used.insert(e_sr.index);
                }
            }
        }
    }
    used.len()
}

/// Dump edge usage stats: total tshape edges vs edges referenced by faces.
pub fn dump_edge_stats(t: &rcad_kernel::topods::BRep, label: &str) {
    let mut n_tshape_edges = 0usize;
    for ts in &t.tshapes {
        if matches!(ts.as_ref(), rcad_kernel::topods::TShape::Edge(_)) { n_tshape_edges += 1; }
    }
    let n_face_refs = count_edges_in_faces(t);
    eprintln!("[TRACE] {}: tshape_edges={} edges_in_faces={}",
        label, n_tshape_edges, n_face_refs);
}
