//! OCCT BRepLib::EncodeRegularity (TKTopAlgo/BRepLib — BRepLib.cxx) — GAP
//! carriers (plan §0.6) for the two forms consumed by
//! BRepFill_Sweep::BuildShell (part B, cxx L3063 / L3132-3135):
//! - EncodeRegularity(S, Tol) — BRepLib.cxx EncodeRegularity over a shape
//!   (TopExp::MapShapesAndAncestors EDGE/FACE walk + regularity encoding);
//! - EncodeRegularity(E, F1, F2, Tol) — the single-edge form.
//!
//! The real bodies compute the tangency between the faces sharing an edge
//! (BRepLib_UpdateEdgeTol / BRep_Tool regularity machinery); they are out of
//! this dispatch and the OCCT failure path is preserved.

use rcad_kernel::topo::topods::{BRep, GeomAbsShape, Shape};

/// OCCT BRepLib::EncodeRegularity(const TopoDS_Shape& S, const double Tol).
pub fn encode_regularity(_brep: &mut BRep, _s: &Shape, _tol: f64) {
    panic!(
        "GAP: BRepLib::EncodeRegularity(S, Tol) (TKTopAlgo/BRepLib.cxx) is \
         not translated — see file header (plan section 0.6)"
    );
}

/// OCCT BRepLib::EncodeRegularity(const TopoDS_Edge& E, const TopoDS_Face& F1,
/// const TopoDS_Face& F2, const double Tol).
pub fn encode_regularity_edge(
    _brep: &mut BRep,
    _e: &Shape,
    _f1: &Shape,
    _f2: &Shape,
    _tol: f64,
) {
    panic!(
        "GAP: BRepLib::EncodeRegularity(E, F1, F2, Tol) (TKTopAlgo/BRepLib.cxx) \
         is not translated — see file header (plan section 0.6)"
    );
}

/// OCCT GeomAbs_G1 — the continuity order requested by the sweep callers.
pub const GEOM_ABS_G1: GeomAbsShape = GeomAbsShape::G1;
