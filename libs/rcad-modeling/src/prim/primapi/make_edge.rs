//! OCCT BRepBuilderAPI_MakeEdge (TKTopAlgo) — straight edge between two points.
//!
//! OCCT BRepBuilderAPI_MakeEdge(gp_Pnt, gp_Pnt): builds a single straight
//! edge from P1 to P2.  rcad equivalent: a line curve edge with two vertices.

use glam::DVec3;
use rcad_kernel::geom::{Curve3, Line3};
use rcad_kernel::topods::{self, BRep};

/// OCCT BRepBuilderAPI_MakeEdge(P1, P2) — straight edge from `p1` to `p2`.
pub fn make_edge_brep(p1: DVec3, p2: DVec3) -> Result<BRep, crate::BuildError> {
    let mut brep = BRep::new();
    let v1 = brep.add_tvertex(p1);
    let v2 = brep.add_tvertex(p2);
    let dir = (p2 - p1).normalize_or_zero();
    if dir.length_squared() < 1e-30 {
        return Err(crate::BuildError::DegenerateGeometry("degenerate edge"));
    }
    let range = [0.0, (p2 - p1).length()];
    let line = Line3::new(p1, dir);
    brep.add_tedge(Some(Curve3::Line(line)), v1, v2, range);
    Ok(brep)
}

/// OCCT BRepBuilderAPI_MakeEdge(Line) — edge from a 3D line.
pub fn make_edge_line_brep(line: &Line3) -> Result<BRep, crate::BuildError> {
    use rcad_kernel::geom::CurveEval;
    let p1 = line.point_at(0.0);
    let p2 = line.point_at(1.0);
    make_edge_brep(p1, p2)
}

/// Number of edges in a BRep (OCCT TopExp_Explorer(shape, TopAbs_EDGE) count
/// with distinct TShapes).
pub fn edge_count_brep(brep: &BRep) -> usize {
    let mut seen = std::collections::HashSet::new();
    for ts in &brep.tshapes {
        if let topods::TShape::Edge(_) = &**ts {
            seen.insert(std::sync::Arc::as_ptr(ts) as u64);
        }
    }
    seen.len()
}
