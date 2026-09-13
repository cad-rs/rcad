//! OCCT `BRep_Tool::Tolerance` — the kernel-local canonical reader.
//!
//! OCCT hosts the three `Tolerance` overloads in `TKBRep/BRep/BRep_Tool.cxx`;
//! the rcad equivalents live on the `BRepTool` trait in
//! [`crate::topo::topods`].  The shared body — which floors the stored
//! tolerance at `Precision::Confusion()` — is factored out here so the kernel
//! has exactly one copy, and re-exported from `topods` (the BRep_Tool home),
//! so the family stays discoverable from the trait that fronts it.
//!
//! rcad-algo cannot be the host: `rcad-kernel` must not depend on
//! `rcad-algo`.  rcad-algo keeps its own canonical
//! (`rcad_algo::brep_algo::tool::brep_tool_tolerance`).

use crate::core::precision::CONFUSION;
use crate::topo::topods::TShape;

/// OCCT `BRep_Tool::Tolerance` — the three shape-kind overloads, each of which
/// floors the stored tolerance at `Precision::Confusion()`:
///   - `BRep_Tool::Tolerance(const TopoDS_Vertex& V)`  BRep_Tool.cxx L1314-1333
///   - `BRep_Tool::Tolerance(const TopoDS_Edge& E)`    BRep_Tool.cxx L881-895
///   - `BRep_Tool::Tolerance(const TopoDS_Face& F)`    BRep_Tool.cxx L137-149
/// All three bodies are identical: `p = TE->Tolerance(); pMin =
/// Precision::Confusion(); if (p > pMin) return p; else return pMin;`
/// OCCT has no generic `BRep_Tool::Tolerance(const TopoDS_Shape&)`
/// dispatcher (BRep_Tool.hxx declares only the three overloads), so the match
/// over the shape kind stands in for the C++ overload resolution.  The OCCT
/// vertex body throws Standard_NullObject when the TVertex is null; the rcad
/// TShape::Vertex payload is non-nullable, so that branch has no counterpart.
/// For a shape kind without such an overload OCCT has no answer at all; rcad
/// returns 0.0, which is what the pre-canonical readers returned.
///
/// The argument is the TShape payload rather than a `Shape`, because the two
/// rcad call-site families reach the TShape differently and the OCCT
/// equivalence holds for both: `BRepTool` is implemented on `BRep` and reads
/// the pool slot `BRep::tshapes[Shape::index]` (the store every rcad mutator
/// edits — `BRep::edge_mut` even copy-on-writes the slot, so a caller's
/// `Shape::data` Arc can be a stale snapshot of it), while the
/// `BRepAdaptor_*::Tolerance` forwarders hold only the `Shape` handle and read
/// `Shape::data`, exactly as OCCT's adaptor reads through `myFace` / `myEdge`.
///
/// OCCT `BRepAdaptor_Surface::Tolerance` (BRepAdaptor_Surface.cxx L92-95) and
/// `BRepAdaptor_Curve::Tolerance` (BRepAdaptor_Curve.cxx L146-149) are both
/// literally `return BRep_Tool::Tolerance(myFace / myEdge);` — the same OCCT
/// function, not a separate one.
pub fn brep_tool_tolerance(ts: &TShape) -> f64 {
    // OCCT `constexpr double pMin = Precision::Confusion();`
    const P_MIN: f64 = CONFUSION;
    match ts {
        // OCCT BRep_Tool.cxx L1314-1333.
        TShape::Vertex(vd) => {
            let p = vd.tolerance;
            if p > P_MIN { p } else { P_MIN }
        }
        // OCCT BRep_Tool.cxx L881-895.
        TShape::Edge(ed) => {
            let p = ed.tolerance;
            if p > P_MIN { p } else { P_MIN }
        }
        // OCCT BRep_Tool.cxx L137-149.
        TShape::Face(fd) => {
            let p = fd.tolerance;
            if p > P_MIN { p } else { P_MIN }
        }
        _ => 0.0,
    }
}
