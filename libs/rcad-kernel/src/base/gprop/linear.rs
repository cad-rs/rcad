//! OCCT BRepGProp::LinearProperties — total edge length of a shape.
//!
//! OCCT BRepGProp.cxx L55-110: LinearProperties(shape, props, skipShared)
//! computes the total mass (length) of the boundary edges.  With
//! skipShared=false the edges are counted once per face occurrence (a shared
//! edge belongs to two faces, so it is counted twice); with skipShared=true
//! each distinct edge is counted once.

use std::collections::HashSet;

use crate::base::gcpnts::arc_length;
use crate::topods::{self, TShape};

/// Length of the edge whose TShape pointer equals `ptr`.
fn edge_length(brep: &topods::BRep, ptr: u64) -> f64 {
    for ts in &brep.tshapes {
        if std::sync::Arc::as_ptr(ts) as u64 == ptr {
            if let TShape::Edge(ed) = &**ts {
                if let Some(c) = &ed.curve {
                    return arc_length(c, ed.range[0], ed.range[1]).abs();
                }
            }
        }
    }
    0.0
}

/// OCCT BRepGProp::LinearProperties(shape, props, skipShared).
pub fn linear_properties(brep: &topods::BRep, skip_shared: bool) -> f64 {
    let mut total = 0.0;
    if skip_shared {
        // Each distinct edge counted once.
        let mut seen: HashSet<u64> = HashSet::new();
        for ts in &brep.tshapes {
            if let TShape::Edge(ed) = &**ts {
                if seen.insert(std::sync::Arc::as_ptr(ts) as u64) {
                    if let Some(c) = &ed.curve {
                        total += arc_length(c, ed.range[0], ed.range[1]).abs();
                    }
                }
            }
        }
    } else {
        // Each edge counted once per face occurrence (per wire).
        for ts in &brep.tshapes {
            if let TShape::Face(fd) = &**ts {
                for ws in std::iter::once(&fd.outer_wire).chain(fd.inner_wires.iter()) {
                    for wt in &brep.tshapes {
                        if std::sync::Arc::as_ptr(wt) as u64 == ws.ptr_id() {
                            if let TShape::Wire(wd) = &**wt {
                                for e in &wd.edges {
                                    total += edge_length(brep, e.ptr_id());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    total
}
