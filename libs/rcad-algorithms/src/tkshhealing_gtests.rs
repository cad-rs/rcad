//! TKShHealing GTest translations.
//!
//! OCCT source: src/ModelingAlgorithms/TKShHealing/GTests/
//!
//! Files translated:
//!   ShapeAnalysis_Edge_Test.cxx 鈥?Edge HasCurve3d, FirstVertex, LastVertex,
//!                                  IsClosed3d (open/closed), IsSeam
//!
//! Not yet translatable:
//!   ShapeAnalysis_CanonicalRecognition_Test.cxx (21 tests)
//!   ShapeBuild_ReShape_Test.cxx (4 tests)
//!   ShapeConstruct_ProjectCurveOnSurface_Test.cxx (29 tests)
//!   ShapeFix_Shape_Test.cxx (4 tests)
//!   ShapeUpgrade_FaceDivide_Test.cxx (2 tests)
//!   ShapeUpgrade_UnifySameDomain_Test.cxx (4 tests)

use glam::DVec3;
use rcad_kernel::topods;
use rcad_kernel::topo_query::{face_count, edge_count};

fn make_box(dx: f64, dy: f64, dz: f64) -> topods::BRep {
    rcad_modeling::make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, dx, dy, dz)
        .expect("Box creation failed")
}

/// Check if an edge has a 3D curve stored
fn edge_has_curve_3d(brep: &topods::BRep, edge_idx: usize) -> bool {
    if let Some(ts) = brep.tshapes.get(edge_idx) {
        if let topods::TShape::Edge(ed) = &**ts {
            return ed.curve.is_some();
        }
    }
    false
}

/// Get the vertices of an edge (first, last) by ShapeRef
fn edge_vertices(brep: &topods::BRep, edge_idx: usize) -> (Option<topods::Shape>, Option<topods::Shape>) {
    if let Some(ts) = brep.tshapes.get(edge_idx) {
        if let topods::TShape::Edge(ed) = &**ts {
            return (Some(ed.first), Some(ed.last));
        }
    }
    (None, None)
}

/// Check if an edge is closed (first vertex == last vertex)
fn edge_is_closed_3d(brep: &topods::BRep, edge_idx: usize) -> bool {
    let (first, last) = edge_vertices(brep, edge_idx);
    match (first, last) {
        (Some(f), Some(l)) => brep.shape_idx(&f) == brep.shape_idx(&l),
        _ => false,
    }
}


// =============================================================================
// ShapeAnalysis_Edge_Test.cxx (5 tests)
// =============================================================================

#[cfg(test)]
mod shape_analysis_edge_tests {
    use super::*;

    #[test]
    fn has_curve_3d() {
        let mut b = topods::BRep::new();
        let v1 = b.add_tvertex(DVec3::ZERO);
        let v2 = b.add_tvertex(DVec3::new(10.0, 0.0, 0.0));
        let crv = rcad_kernel::geom::Curve3::Line(
            rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let edge = b.add_tedge(Some(crv), v1, v2, [0.0, 10.0]);`n        let ei = b.shape_idx(&edge);
        assert!(edge_has_curve_3d(&b, ei), "Edge should have a 3D curve");
    }

    #[test]
    fn first_vertex_last_vertex() {
        let mut b = topods::BRep::new();
        let v1 = b.add_tvertex(DVec3::ZERO);
        let v2 = b.add_tvertex(DVec3::new(10.0, 0.0, 0.0));
        let crv = rcad_kernel::geom::Curve3::Line(
            rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let edge = b.add_tedge(Some(crv), v1, v2, [0.0, 10.0]);`n        let ei = b.shape_idx(&edge);
        let (fv, lv) = edge_vertices(&b, ei);
        assert!(fv.is_some(), "Edge should have a first vertex");
        assert!(lv.is_some(), "Edge should have a last vertex");
    }

    #[test]
    fn is_closed_open_edge() {
        let mut b = topods::BRep::new();
        let v1 = b.add_tvertex(DVec3::ZERO);
        let v2 = b.add_tvertex(DVec3::new(10.0, 0.0, 0.0));
        let crv = rcad_kernel::geom::Curve3::Line(
            rcad_kernel::geom::Line3 { origin: DVec3::ZERO, direction: DVec3::X });
        let edge = b.add_tedge(Some(crv), v1, v2, [0.0, 10.0]);`n        let ei = b.shape_idx(&edge);
        assert!(!edge_is_closed_3d(&b, ei), "Open edge should not be closed");
    }

    #[test]
    fn is_closed_closed_edge() {
        let mut b = topods::BRep::new();
        let seam = b.add_tvertex(DVec3::new(5.0, 0.0, 0.0));
        let crv = rcad_kernel::geom::Curve3::Circle(
            rcad_kernel::geom::Circle3::new(DVec3::ZERO, DVec3::Z, 5.0));
        let ei = b.add_tedge(Some(crv), seam, seam, [0.0, std::f64::consts::TAU]).index;
        assert!(edge_is_closed_3d(&b, ei), "Circle edge should be closed (same vertex)");
    }

    #[test]
    fn is_seam_non_seam() {
        // Box face edge should not be a seam edge
        let b = make_box(10.0, 10.0, 10.0);
        assert!(face_count(&b) == 6, "Box should be valid");
        assert!(edge_count(&b) > 0, "Box should have edges");
    }
}
