//! BRepAlgo-style helpers (migrated / re-implemented for the OCCT grid tests).
//!
//! The legacy `rcad-algorithms` crate exposed `brep_algo::total_edge_length`;
//! it was never implemented there, so this is a fresh implementation over the
//! current `topods::BRep` pool.

use rcad_kernel::topods::{self, TShape};

/// Total arc length of every edge in the BRep.
///
/// Mirrors the `checkprops -l` style "curve length" used by the generated
/// OCCT boolean tests. Each `TShape::Edge` contributes its 3D curve arc
/// length over its parameter range.
pub fn total_edge_length(brep: &topods::BRep) -> f64 {
    let mut total = 0.0;
    for ts in &brep.tshapes {
        if let TShape::Edge(ed) = ts.as_ref() {
            if let Some(curve) = &ed.curve {
                let (t1, t2) = (ed.range[0], ed.range[1]);
                total += rcad_kernel::base::gcpnts::abscissa_point::arc_length(curve, t1, t2);
            }
        }
    }
    total
}
