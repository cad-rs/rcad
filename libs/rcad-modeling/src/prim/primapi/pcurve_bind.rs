//! The shared writer-migration dual-write helper for the primapi pcurve
//! sites.  Every pcurve write in the primitives is the OCCT
//! `BRep_Builder::UpdateEdge` statement chain, whose representation arm is
//! the static UpdateCurves (BRep_Builder.cxx L104-167): L133-146 removes any
//! existing curve-on-surface representation of the same (S, L), then L149-167
//! appends the new `BRep_CurveOnSurface(C, S, L)` representation to the edge
//! TShape's single curve-representation list (BRep_TEdge myCurves).  The
//! `pcurves` map insert at the call sites is retained for the map-era readers
//! during the writer migration; the representation is the authority.

use rcad_kernel::geom::Curve2d;
use rcad_kernel::topods::{CurveRepresentation, TEdgeData};

/// The UpdateCurves single-pcurve representation arm (BRep_Builder.cxx
/// L133-146 removal + L149-167 append).
pub(crate) fn bind_pcurve_representation(
    ed: &mut TEdgeData,
    key: (u64, u32),
    pc: Curve2d,
    range: [f64; 2],
) {
    ed.representations.retain(|a_cr| match a_cr {
        CurveRepresentation::CurveOnSurface { face, .. }
        | CurveRepresentation::CurveOnClosedSurface { face, .. } => *face != key,
        _ => true,
    });
    ed.representations
        .push(CurveRepresentation::CurveOnSurface {
            face: key,
            pcurve: pc,
            range,
        });
}
