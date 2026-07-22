//! Precision constants and per-entity tolerance query helpers.
//!
//! Analogous to OCCT's `Precision` class and `BRep_Builder::UpdateVertex` /
//! `BRep_Builder::UpdateEdge` tolerance API.
//!
//! Tolerances are stored on individual TShapes (TVertexData.tolerance,
//! TEdgeData.tolerance, TFaceData.tolerance), matching OCCT's per-entity model.
//! When absent or zero the functions fall back to the `CONFUSION` constant.

use crate::topods;

// ── Precision constants ───────────────────────────────────────────────────────

/// Point-coincidence tolerance.
/// Analogous to `Precision::Confusion()` = 1e-7 in OCCT.
pub const CONFUSION: f64 = 1e-7;

/// Angular tolerance (radians).
/// Analogous to `Precision::Angular()` = 1e-12 in OCCT.
pub const ANGULAR: f64 = 1e-12;

/// Tessellation / approximation tolerance.
/// Analogous to `Precision::Approximation()` = 1e-6 in OCCT.
pub const APPROXIMATION: f64 = 1e-6;

/// Intersection tolerance.  Used by intersection algorithms to decide
/// when a solution is reached.
/// Analogous to `Precision::Intersection()` = Confusion / 100 = 1e-9.
pub const INTERSECTION: f64 = CONFUSION * 0.01;

/// Machine epsilon for double-precision floating-point.
/// Analogous to `Precision::Computational()` = DBL_EPSILON ≈ 2.22e-16.
pub const COMPUTATIONAL: f64 = f64::EPSILON;
pub const SQUARE_COMPUTATIONAL: f64 = COMPUTATIONAL * COMPUTATIONAL;

/// Parametric tolerance for small-range checks.
/// Analogous to `Precision::PConfusion()` = 1e-12 in OCCT.
pub const PCONFUSION: f64 = 1e-12;

/// Square of CONFUSION — used for squared-distance comparisons.
pub const SQUARE_CONFUSION: f64 = CONFUSION * CONFUSION;

/// Analogous to `Precision::Infinite()` = 2e+100 in OCCT.
pub const INFINITE_VALUE: f64 = 2e100;

/// OCCT-aligned: Precision::IsInfinite (Precision.hxx L350-353).
/// OCCT: std::abs(R) >= 0.5 * Precision::Infinite() where Precision::Infinite() = 2e100.
pub fn is_infinite_value(r: f64) -> bool { r.abs() >= 0.5 * INFINITE_VALUE }

/// OCCT-aligned: Precision::Epsilon (Precision.hxx L336-341).
/// Returns `|thePar| * 1e-12` (or `1e-12` if thePar is zero).
/// Used by BRepLib::FindValidRange for parametric convergence threshold.
pub fn parametric_epsilon(the_par: f64) -> f64 {
    if the_par == 0.0 { 1e-12 } else { the_par.abs() * 1e-12 }
}

/// OCCT-aligned: Precision::IsPositiveInfinite (Precision.hxx L357-360).
pub fn is_positive_infinite_value(r: f64) -> bool { r >= 0.5 * INFINITE_VALUE }

/// OCCT-aligned: Precision::IsNegativeInfinite (Precision.hxx L364-367).
pub fn is_negative_infinite_value(r: f64) -> bool { r <= -0.5 * INFINITE_VALUE }

// ── Per-shape tolerance helpers (topods::BRep) ─────────────────────────────

/// Vertex tolerance from a TShape vertex.
fn vtol(vd: &topods::TVertexData) -> f64 {
    if vd.tolerance > 0.0 { vd.tolerance } else { CONFUSION }
}

/// Edge tolerance from a TShape edge.
fn etol(ed: &topods::TEdgeData) -> f64 {
    if ed.tolerance > 0.0 { ed.tolerance } else { CONFUSION }
}

/// Face tolerance from a TShape face.
fn ftol(fd: &topods::TFaceData) -> f64 {
    if fd.tolerance > 0.0 { fd.tolerance } else { CONFUSION }
}

/// Returns the tolerance for the vertex at `tshape_idx`.
pub fn vertex_tolerance(brep: &topods::BRep, tshape_idx: usize) -> f64 {
    brep.tshapes.get(tshape_idx)
        .and_then(|ts| {
            if let topods::TShape::Vertex(vd) = &**ts { Some(vtol(vd)) } else { None }
        })
        .unwrap_or(CONFUSION)
}

/// Set the tolerance for the vertex at `tshape_idx`.
pub fn set_vertex_tolerance(brep: &mut topods::BRep, tshape_idx: usize, tol: f64) {
    if let Some(ts) = brep.tshapes.get(tshape_idx) {
        if let topods::TShape::Vertex(vd) = &**ts {
            let mut new_vd = vd.clone();
            new_vd.tolerance = tol.max(CONFUSION);
            brep.tshapes[tshape_idx] = std::sync::Arc::new(topods::TShape::Vertex(new_vd));
        }
    }
}

/// Update the tolerance for the vertex at `tshape_idx` — sets to max of existing and new.
pub fn update_vertex_tolerance(brep: &mut topods::BRep, tshape_idx: usize, tol: f64) {
    let cur = vertex_tolerance(brep, tshape_idx);
    if tol > cur {
        set_vertex_tolerance(brep, tshape_idx, tol);
    }
}

/// Returns the tolerance for the edge at `tshape_idx`.
pub fn edge_tolerance(brep: &topods::BRep, tshape_idx: usize) -> f64 {
    brep.tshapes.get(tshape_idx)
        .and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { Some(etol(ed)) } else { None }
        })
        .unwrap_or(CONFUSION)
}

/// Set the tolerance for the edge at `tshape_idx`.
pub fn set_edge_tolerance(brep: &mut topods::BRep, tshape_idx: usize, tol: f64) {
    // We can't directly mutate through Arc, so rebuild the edge with new tolerance
    // This requires replacing the TShape at the given index
    if let Some(ts) = brep.tshapes.get(tshape_idx) {
        if let topods::TShape::Edge(ed) = &**ts {
            let mut new_ed = ed.clone();
            new_ed.tolerance = tol.max(CONFUSION);
            brep.tshapes[tshape_idx] = std::sync::Arc::new(topods::TShape::Edge(new_ed));
        }
    }
}

/// Update the tolerance for the edge at `tshape_idx` — sets to max of existing and new.
pub fn update_edge_tolerance(brep: &mut topods::BRep, tshape_idx: usize, tol: f64) {
    let cur = edge_tolerance(brep, tshape_idx);
    if tol > cur {
        set_edge_tolerance(brep, tshape_idx, tol);
    }
}

/// Returns the tolerance for the face at `tshape_idx`.
pub fn face_tolerance(brep: &topods::BRep, tshape_idx: usize) -> f64 {
    brep.tshapes.get(tshape_idx)
        .and_then(|ts| {
            if let topods::TShape::Face(fd) = &**ts { Some(ftol(fd)) } else { None }
        })
        .unwrap_or(CONFUSION)
}

/// Set the tolerance for the face at `tshape_idx`.
pub fn set_face_tolerance(brep: &mut topods::BRep, tshape_idx: usize, tol: f64) {
    if let Some(ts) = brep.tshapes.get(tshape_idx) {
        if let topods::TShape::Face(fd) = &**ts {
            let mut new_fd = fd.clone();
            new_fd.tolerance = tol.max(CONFUSION);
            brep.tshapes[tshape_idx] = std::sync::Arc::new(topods::TShape::Face(new_fd));
        }
    }
}

/// Update the tolerance for the face at `tshape_idx` — sets to max of existing and new.
pub fn update_face_tolerance(brep: &mut topods::BRep, tshape_idx: usize, tol: f64) {
    let cur = face_tolerance(brep, tshape_idx);
    if tol > cur {
        set_face_tolerance(brep, tshape_idx, tol);
    }
}

/// Initialize tolerance arrays for a new topods BRep.
/// Sets all vertex/edge/face tolerances to CONFUSION.
pub fn finalize_tolerance_hierarchy(brep: &mut topods::BRep) {
    for ts in &mut brep.tshapes {
        let ts_mut = std::sync::Arc::make_mut(ts);
        match ts_mut {
            topods::TShape::Vertex(vd) => {
                if vd.tolerance <= 0.0 { vd.tolerance = CONFUSION; }
            }
            topods::TShape::Edge(ed) => {
                if ed.tolerance <= 0.0 { ed.tolerance = CONFUSION; }
            }
            topods::TShape::Face(fd) => {
                if fd.tolerance <= 0.0 { fd.tolerance = CONFUSION; }
            }
            _ => {}
        }
    }
}

/// Compute per-edge SameParameter flag.
/// In OCCT's model this is stored on the edge itself.
pub fn edge_same_parameter(brep: &topods::BRep, tshape_idx: usize) -> bool {
    brep.tshapes.get(tshape_idx)
        .and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { Some(ed.same_parameter) } else { None }
        })
        .unwrap_or(true)
}

/// Compute per-edge SameRange flag.
pub fn edge_same_range(brep: &topods::BRep, tshape_idx: usize) -> bool {
    brep.tshapes.get(tshape_idx)
        .and_then(|ts| {
            if let topods::TShape::Edge(ed) = &**ts { Some(ed.same_range) } else { None }
        })
        .unwrap_or(true)
}

/// Compute per-face natural domain flag.
pub fn face_domain(brep: &topods::BRep, tshape_idx: usize) -> Option<[f64; 4]> {
    brep.tshapes.get(tshape_idx)
        .and_then(|ts| {
            if let topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None }
        })
}

/// Model-level tolerance: max of all vertex/edge tolerances.
pub fn model_tolerance(brep: &topods::BRep) -> f64 {
    let mut max_tol = CONFUSION;
    for ts in &brep.tshapes {
        match &**ts {
            topods::TShape::Vertex(vd) => {
                if vd.tolerance > max_tol { max_tol = vd.tolerance; }
            }
            topods::TShape::Edge(ed) => {
                if ed.tolerance > max_tol { max_tol = ed.tolerance; }
            }
            topods::TShape::Face(fd) => {
                if fd.tolerance > max_tol { max_tol = fd.tolerance; }
            }
            _ => {}
        }
    }
    max_tol
}

/// Compute brep SameParameter by checking all edges.
pub fn brep_same_parameter(brep: &topods::BRep) -> bool {
    brep.tshapes.iter().all(|ts| {
        if let topods::TShape::Edge(ed) = &**ts { ed.same_parameter } else { true }
    })
}

/// STEP export uncertainty: max tolerance of all entities, times 10.
pub fn step_export_uncertainty(brep: &topods::BRep) -> f64 {
    let tol = model_tolerance(brep);
    if tol > 0.0 { tol * 10.0 } else { 1e-6 }
}

// ── Legacy helpers (for backward compat with callers that use flat indices) ──
// These map flat vertex/edge/face tshape indices in the tshapes array.

/// Resize tolerance arrays (legacy compat — no-op, tolerances are per-entity).
/// Retained for API compatibility with callers that expect a flat resize operation.
pub fn resize_tolerance_arrays(_brep: &mut topods::BRep) {
    // Tolerances are stored per-entity on each TShape — no flat arrays to resize.
}
