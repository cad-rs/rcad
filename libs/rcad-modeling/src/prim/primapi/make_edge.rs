//! OCCT BRepBuilderAPI_MakeEdge (TKTopAlgo) — straight/circular edge builders.
//!
//! OCCT BRepBuilderAPI_MakeEdge(gp_Pnt, gp_Pnt): builds a single straight
//! edge from P1 to P2.  BRepBuilderAPI_MakeEdge(Geom_Circle[, P1, P2]) builds
//! a circular edge (full circle when the range is omitted, trimmed otherwise).
//! BRepBuilderAPI_MakeEdge(Geom_Line, P1, P2) builds a line edge over the
//! parameter range [P1, P2].  rcad equivalents: line curve edges with two
//! vertices (BRepLib_MakeEdge, TKTopAlgo).

use glam::DVec3;
use rcad_kernel::geom::{Circle3, Curve3, Line3};
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

/// OCCT BRepBuilderAPI_MakeEdge(Line) — edge from a 3D line over [0, 1].
pub fn make_edge_line_brep(line: &Line3) -> Result<BRep, crate::BuildError> {
    make_edge_line_range_brep(line, 0.0, 1.0)
}

/// OCCT BRepBuilderAPI_MakeEdge(Line, P1, P2) — line edge over [p1, p2].
pub fn make_edge_line_range_brep(
    line: &Line3,
    t1: f64,
    t2: f64,
) -> Result<BRep, crate::BuildError> {
    use rcad_kernel::geom::CurveEval;
    let p1 = line.point_at(t1);
    let p2 = line.point_at(t2);
    if (p2 - p1).length_squared() < 1e-30 {
        return Err(crate::BuildError::DegenerateGeometry("degenerate edge"));
    }
    let mut brep = BRep::new();
    let v1 = brep.add_tvertex(p1);
    let v2 = brep.add_tvertex(p2);
    brep.add_tedge(Some(Curve3::Line(*line)), v1, v2, [t1, t2]);
    Ok(brep)
}

/// OCCT BRepBuilderAPI_MakeEdge(Geom_Circle) — full circle edge over the
/// natural parameter range [0, 2*PI].  Both endpoints coincide; the vertex
/// dedup (identity sharing by quantized position) yields a single seam
/// vertex, mirroring BRepLib_MakeEdge::Init for a closed curve.
pub fn make_edge_circle_brep(circle: &Circle3) -> Result<BRep, crate::BuildError> {
    make_edge_circle_range_brep(circle, 0.0, 2.0 * std::f64::consts::PI)
}

/// OCCT BRepBuilderAPI_MakeEdge(Geom_Circle, P1, P2) — trimmed circular edge
/// over [t1, t2].
pub fn make_edge_circle_range_brep(
    circle: &Circle3,
    t1: f64,
    t2: f64,
) -> Result<BRep, crate::BuildError> {
    use rcad_kernel::geom::CurveEval;
    let mut brep = BRep::new();
    let p1 = circle.point_at(t1);
    let p2 = circle.point_at(t2);
    let v1 = brep.add_tvertex(p1);
    let v2 = brep.add_tvertex(p2);
    brep.add_tedge(Some(Curve3::Circle(*circle)), v1, v2, [t1, t2]);
    Ok(brep)
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
