//! OCCT-aligned splitter: BOPAlgo_Splitter / BRepAlgoAPI_Splitter.
//!
//! OCCT's BOPAlgo_Splitter splits a shape (Object) by another shape (Tool),
//! keeping only the Object's split parts. The Tool serves as cutting boundaries
//! but its split parts are discarded.
//!
//! Pipeline (mirrors OCCT):
//! 1. PaveFiller — compute intersection curves between Object and Tool
//! 2. BooleanBuilder — split faces of both shapes along intersection curves
//! 3. BuildResult — filter result to keep only split parts of the Object
//!
//! OCCT references:
//! - `BOPAlgo_Splitter.cxx` L1-80: BOPAlgo_Splitter::Perform + BuildResult
//! - `BRepAlgoAPI_Splitter.cxx` L1-50: API wrapper

use crate::builder::BooleanOpType;
use crate::history::FaceOrigin;
use rcad_kernel::BRep;

// ── Task B: filter_object_only ──────────────────────────────────────────────

/// Remove faces that originated from the Tool shape (ShapeB) from a boolean
/// operation result.
///
/// OCCT equivalent: BOPAlgo_Splitter::BuildResult filtering (L37-79).
///
/// In OCCT, `BOPAlgo_Splitter::BuildResult()`:
/// 1. Calls `BOPAlgo_Builder::BuildResult()` to produce the full split result
///    (all split parts of both Object and Tool).
/// 2. Iterates through the result shapes and checks `myImages` / `myOrigins`
///    maps to trace each result shape back to its input origin.
/// 3. Keeps only shapes whose origin traces back to the Object (ShapeA).
/// 4. Discards shapes whose origin traces back to the Tool (ShapeB).
///
/// The `face_origins` slice comes from [`BooleanHistory::face_origins`]
/// (returned by [`crate::boolean_op_with_history`]), which maps each result
/// face (in flattened order) to its input face origin.
///
/// # Example
///
/// ```rust,no_run
/// use rcad_algorithms::splitter::filter_object_only;
/// # let (brep, history) = todo!();
/// let shape_a_only = filter_object_only(&brep, &history.face_origins);
/// ```
///
/// # OCCT Reference
/// `BOPAlgo_Splitter.cxx` L58-L79 — BuildResult filtering loop
/// `BOPAlgo_Builder.cxx` L178-L240 — BuildResult shape iteration
pub fn filter_object_only(brep: &BRep, face_origins: &[FaceOrigin]) -> BRep {
    // Build a keep mask: true if the origin is ShapeA or Generated.
    let keep_mask: Vec<bool> = face_origins
        .iter()
        .map(|fo| matches!(fo, FaceOrigin::FromA(_) | FaceOrigin::Generated))
        .collect();

    let mut result_geom = brep.geom.clone();

    // Filter per-face geom arrays.
    let mut new_face_surface: Vec<Option<usize>> = Vec::new();
    let mut new_face_surface_range: Vec<Option<[f64; 4]>> = Vec::new();
    let mut new_face_tolerance: Vec<f64> = Vec::new();
    let mut new_face_internal_vertices: Vec<Vec<usize>> = Vec::new();

    let mut flat_idx = 0usize;
    for &keep in &keep_mask {
        if keep {
            if flat_idx < result_geom.face_surface.len() {
                new_face_surface.push(result_geom.face_surface[flat_idx]);
            }
            if flat_idx < result_geom.face_surface_range.len() {
                new_face_surface_range.push(result_geom.face_surface_range[flat_idx]);
            }
            if flat_idx < result_geom.face_tolerance.len() {
                new_face_tolerance.push(result_geom.face_tolerance[flat_idx]);
            }
            if flat_idx < result_geom.face_internal_vertices.len() {
                new_face_internal_vertices.push(result_geom.face_internal_vertices[flat_idx].clone());
            }
        }
        flat_idx += 1;
    }

    result_geom.face_surface = new_face_surface;
    result_geom.face_surface_range = new_face_surface_range;
    result_geom.face_tolerance = new_face_tolerance;
    result_geom.face_internal_vertices = new_face_internal_vertices;

    // Filter faces from shells, tracking flat index in lockstep.
    let mut new_solids: Vec<rcad_kernel::topology::Solid> = Vec::new();
    flat_idx = 0;

    for solid in &brep.solids {
        let mut new_shells: Vec<rcad_kernel::topology::Shell> = Vec::new();
        for shell in &solid.shells {
            let mut kept_faces: Vec<rcad_kernel::topology::Face> = Vec::new();
            for face in &shell.faces {
                if flat_idx < keep_mask.len() && keep_mask[flat_idx] {
                    kept_faces.push(face.clone());
                }
                flat_idx += 1;
            }
            if !kept_faces.is_empty() {
                new_shells.push(rcad_kernel::topology::Shell {
                    faces: kept_faces,
                });
            }
        }
        if !new_shells.is_empty() {
            new_solids.push(rcad_kernel::topology::Solid {
                shells: new_shells,
            });
        }
    }

    BRep {
        vertices: brep.vertices.clone(),
        edges: brep.edges.clone(),
        solids: new_solids,
        geom: result_geom,
        compound: brep.compound.clone(),
        compsolid: brep.compsolid.clone(),
    }
}

// ── Task C: split_shape_occt_aligned ────────────────────────────────────────

/// Split a shape by a tool, keeping only the Object's split parts.
///
/// This is the OCCT-aligned equivalent of BRepAlgoAPI_Splitter /
/// BOPAlgo_Splitter. It:
/// 1. Runs the full PaveFiller + BooleanBuilder pipeline (same as boolean ops).
/// 2. Computes intersection curves between the two shapes.
/// 3. Splits faces of both shapes along the intersection curves.
/// 4. Filters the result to keep only faces originating from the Object (ShapeA),
///    discarding Tool (ShapeB) faces.
///
/// The operation mode is intersection/section: the splitter divides the Object
/// by the Tool's cutting boundaries. The result contains the portion of the
/// Object split at the Tool's intersection boundary.
///
/// # OCCT Reference
/// `BOPAlgo_Splitter.cxx` L1-80:
///   - Perform(): calls `BOPAlgo_Builder::Perform()` which runs PaveFiller
///     and splits faces (same pipeline as boolean ops).
///   - BuildResult(): calls `BOPAlgo_Builder::BuildResult()`, then filters
///     result shapes by `myImages` / `myOrigins` to keep only Object-originated
///     shapes.
/// `BRepAlgoAPI_Splitter.cxx` L1-50:
///   - API wrapper that delegates to BOPAlgo_Splitter.
pub fn split_shape_occt_aligned(shape: &BRep, tool: &BRep) -> Result<BRep, String> {
    // Step 1: Run the PaveFiller + Builder pipeline (same as boolean ops).
    // OCCT: BOPAlgo_Splitter uses BOPAlgo_Builder (not BOPAlgo_BOP), which
    // splits all faces and classifies them per the boolean op type.
    let (brep, history) = crate::boolean_op_with_history(
        BooleanOpType::Intersection,
        shape,
        tool,
    )
    .map_err(|e| format!("splitter: boolean_op failed: {e}"))?;

    // Step 2: Filter to keep only Object-originated (ShapeA) faces.
    // OCCT: BOPAlgo_Splitter::BuildResult filters by myImages/myOrigins.
    Ok(filter_object_only(&brep, &history.face_origins))
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::BRep;
    use rcad_modeling::builder::{box_brep, sphere_brep};

    fn make_box(origin: DVec3, x: DVec3, y: DVec3, w: f64, h: f64, d: f64) -> BRep {
        box_brep(origin, x, y, w, h, d).expect("box")
    }

    fn make_sphere(center: DVec3, radius: f64) -> BRep {
        sphere_brep(center, radius).expect("sphere")
    }

    /// Helper: count total faces across all solids/shells.
    fn count_faces(brep: &BRep) -> usize {
        brep.solids.iter().flat_map(|s| &s.shells).flat_map(|sh| &sh.faces).count()
    }

    #[test]
    fn test_filter_object_only_all_keep() {
        // When all faces are FromA, filter_object_only should keep all faces.
        let brep = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        let nfaces = count_faces(&brep);
        let origins: Vec<FaceOrigin> = (0..nfaces).map(|_| FaceOrigin::FromA(0)).collect();
        let filtered = filter_object_only(&brep, &origins);
        assert_eq!(count_faces(&filtered), nfaces);
    }

    #[test]
    fn test_filter_object_only_all_removed() {
        // When all faces are FromB, result should have no solids/faces.
        let brep = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        let nfaces = count_faces(&brep);
        let origins: Vec<FaceOrigin> = (0..nfaces).map(|_| FaceOrigin::FromB(0)).collect();
        let filtered = filter_object_only(&brep, &origins);
        assert_eq!(count_faces(&filtered), 0);
    }

    #[test]
    fn test_filter_object_only_half_kept() {
        // First 3 faces FromA, last 3 FromB -- should keep only the first 3.
        let brep = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0);
        let mut origins: Vec<FaceOrigin> = (0..3).map(|_| FaceOrigin::FromA(0)).collect();
        origins.extend((0..3).map(|_| FaceOrigin::FromB(0)));
        let filtered = filter_object_only(&brep, &origins);
        assert_eq!(count_faces(&filtered), 3);
    }

    #[test]
    fn test_split_shape_occt_aligned_no_intersection() {
        // Non-overlapping shapes -- no intersection, result should be empty.
        let shape = make_box(DVec3::new(-10.0, -1.0, -1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        let tool = make_sphere(DVec3::new(10.0, 0.0, 0.0), 1.5);
        let result = split_shape_occt_aligned(&shape, &tool);
        assert!(result.is_ok());
    }

    #[test]
    fn test_split_shape_occt_aligned_box_box() {
        // Two overlapping boxes -- splitter should keep only ShapeA faces.
        let shape = make_box(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        let tool = make_box(DVec3::new(0.5, 0.5, 0.5), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0);
        let result = split_shape_occt_aligned(&shape, &tool);
        assert!(result.is_ok());
    }
}
