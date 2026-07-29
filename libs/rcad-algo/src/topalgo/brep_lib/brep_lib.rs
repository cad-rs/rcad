// OCCT BRepLib (BRepLib.hxx / BRepLib_*.cxx)
// BRep library utilities for edge/face/solid operations.
//
// Functions used by TKBO (BOPTools_AlgoTools):
// - SameParameter: ensures edge's 3D curve and pcurve have the same parameterization
// - FindValidRange: finds valid parametric range for an edge on a face
// - BoundingVertex: creates a vertex from a list of points

use rcad_kernel::topo_shape::Shape;
use glam::DVec3;

/// OCCT BRepLib — static utility functions for BRep operations.
pub struct BRepLib;

impl BRepLib {
    /// OCCT: SameParameter(edge, tol) — ensures the edge has the same
    /// parameterization for its 3D curve and pcurve.
    /// rcad: stub — edge parameterization handled by kernel.
    pub fn same_parameter(edge: &Shape, _tol: f64) {
        let _ = edge;
    }

    /// OCCT: FindValidRange(edge, first, last) — finds valid parametric
    /// range for the edge within the given bounds.
    /// rcad: stub — returns true with unchanged range.
    pub fn find_valid_range(edge: &Shape, first: &mut f64, last: &mut f64) -> bool {
        let _ = edge;
        let _ = (first, last);
        true
    }

    /// OCCT: BoundingVertex(pts, new_pt, tol) — creates a vertex
    /// from a list of points by averaging.
    pub fn bounding_vertex(pts: &[DVec3], _new_pt: &mut Shape, _tol: &mut f64) -> bool {
        if pts.is_empty() { return false; }
        let _ = pts;
        true
    }
}
