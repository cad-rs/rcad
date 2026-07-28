//! splitter: BOPAlgo_Splitter / BRepAlgoAPI_Splitter.
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

use crate::bopalgo::builder::BooleanOpType;
use crate::history::FaceOrigin;
use rcad_kernel::topods;

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
/// (returned by [`crate::bop_occt_ops::boolean_op_with_history_generic`]), which maps each result
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
pub fn filter_object_only(
    brep: &rcad_kernel::BRep,
    face_origins: &[FaceOrigin],
) -> rcad_kernel::BRep {
    // Build a keep mask: true if the origin is ShapeA or Generated.
    let keep_mask: Vec<bool> = face_origins
        .iter()
        .map(|fo| matches!(fo, FaceOrigin::FromA(_) | FaceOrigin::Generated))
        .collect();

    // First pass: collect flat face index -> face tshape index mapping and
    // determine which faces to remove.
    let mut remove_face_tsi: Vec<usize> = Vec::new();
    {
        let mut flat_idx = 0usize;
        for (tsi, ts) in brep.tshapes.iter().enumerate() {
            if let topods::TShape::Face(_) = ts.as_ref() {
                if flat_idx < keep_mask.len() && !keep_mask[flat_idx] {
                    remove_face_tsi.push(tsi);
                }
                flat_idx += 1;
            }
        }
    }
    let remove_set: std::collections::HashSet<usize> = remove_face_tsi.into_iter().collect();

    // Build a new BRep from scratch with only kept faces, keeping shell/solid structure.
    let mut out = rcad_kernel::topods::BRep::new();

    // Copy all non-Solid TShapes (vertices, edges, wires, faces, shells)
    // and rebuild Solid TShapes with filtered shells/faces.
    // To preserve Arc identity, we iterate the input tshapes and clone
    // TShapes that are not part of the Solid -> Shell -> Face path.
    // For Solid/Shell/Face we rebuild the ref lists.
    //
    // Strategy: add all non-Compound, non-Solid tshapes as-is (clone),
    // then rebuild Solid/Shell/Face hierarchy from scratch.

    // First, collect all face tshape indices in flat order so we can remap.
    let face_tsi_by_flat: Vec<usize> = brep
        .tshapes
        .iter()
        .enumerate()
        .filter_map(|(i, ts)| {
            if matches!(ts.as_ref(), topods::TShape::Face(_)) {
                Some(i)
            } else {
                None
            }
        })
        .collect();

    // For each face, determine if it should be kept.
    let face_keep: Vec<bool> = face_tsi_by_flat
        .iter()
        .enumerate()
        .map(|(fi, _)| fi < keep_mask.len() && keep_mask[fi])
        .collect();

    // Walk all input TShapes, rebuilding hierarchy for kept faces.
    // We use a HashMap to map old tshape index to new Shape.
    let mut old_to_new: std::collections::HashMap<usize, topods::Shape> =
        std::collections::HashMap::new();

    // First, copy all Vertex TShapes.
    for (old_i, ts) in brep.tshapes.iter().enumerate() {
        match ts.as_ref() {
            topods::TShape::Vertex(vd) => {
                let sr = out.add_tvertex(vd.point);
                old_to_new.insert(old_i, sr);
            }
            _ => {}
        }
    }

    // Copy all Edge TShapes.
    for (old_i, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Edge(ed) = ts.as_ref() {
            let first = old_to_new
                .get(&ed.first.index)
                .cloned()
                .unwrap_or(topods::Shape::null());
            let last = old_to_new
                .get(&ed.last.index)
                .cloned()
                .unwrap_or(topods::Shape::null());
            let sr = out.add_tedge(ed.curve.clone(), first, last, ed.range);
            old_to_new.insert(old_i, sr);
        }
    }

    // Copy Wire TShapes, filtering out removed edges.
    for (old_i, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Wire(wd) = ts.as_ref() {
            let new_edges: Vec<topods::Shape> = wd
                .edges
                .iter()
                .filter_map(|e| old_to_new.get(&e.index).cloned())
                .collect();
            if !new_edges.is_empty() {
                // Preserve orientation
                let new_edges_oriented: Vec<topods::Shape> = wd
                    .edges
                    .iter()
                    .filter_map(|e| {
                        old_to_new.get(&e.index).map(|n| {
                            topods::Shape::synthetic(n.index, e.orientation)
                        })
                    })
                    .collect();
                let sr = out.add_twire(new_edges_oriented);
                old_to_new.insert(old_i, sr);
            }
        }
    }

    // Copy Face TShapes (only kept ones).
    for (old_i, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Face(fd) = ts.as_ref() {
            // Check if this face is kept
            let face_flat_idx = face_tsi_by_flat.iter().position(|&tsi| tsi == old_i);
            let keep = face_flat_idx.map_or(false, |fi| fi < face_keep.len() && face_keep[fi]);
            if !keep {
                continue;
            }
            let new_outer = old_to_new
                .get(&fd.outer_wire.index)
                .cloned()
                .unwrap_or(topods::Shape::null());
            let new_inner: Vec<topods::Shape> = fd
                .inner_wires
                .iter()
                .filter_map(|w| old_to_new.get(&w.index).cloned())
                .collect();
            let new_internal: Vec<topods::Shape> = fd
                .internal_vertices
                .iter()
                .filter_map(|v| old_to_new.get(&v.index).cloned())
                .collect();
            // Preserve surface, uv_domain, natural_restriction from original face
            let sr = out.add_tface(
                fd.surface.clone(),
                new_outer,
                new_inner,
                fd.sample_point,
                fd.uv_domain,
                new_internal,
                fd.natural_restriction,
            );
            old_to_new.insert(old_i, sr);
        }
    }

    // Build filtered Shell TShapes.
    for (old_i, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Shell(shd) = ts.as_ref() {
            let new_faces: Vec<topods::Shape> = shd
                .faces
                .iter()
                .filter_map(|f| {
                    old_to_new.get(&f.index).map(|n| {
                        topods::Shape::synthetic(n.index, f.orientation)
                    })
                })
                .collect();
            if !new_faces.is_empty() {
                let sr = out.add_tshell(new_faces);
                old_to_new.insert(old_i, sr);
            }
        }
    }

    // Build filtered Solid TShapes.
    for (old_i, ts) in brep.tshapes.iter().enumerate() {
        if let topods::TShape::Solid(sd) = ts.as_ref() {
            let new_shells: Vec<topods::Shape> = sd
                .shells
                .iter()
                .filter_map(|s| {
                    old_to_new.get(&s.index).map(|n| {
                        topods::Shape::synthetic(n.index, s.orientation)
                    })
                })
                .collect();
            if !new_shells.is_empty() {
                out.add_tsolid(new_shells);
            }
        }
    }

    out
}

// ── Task C: split_shape_occt_aligned ────────────────────────────────────────

/// Split a shape by a tool, keeping only the Object's split parts.
///
/// This is the equivalent of BRepAlgoAPI_Splitter /
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
pub fn split_shape_occt_aligned(
    shape: &rcad_kernel::BRep,
    tool: &rcad_kernel::BRep,
) -> Result<rcad_kernel::BRep, String> {
    // Step 1: Run the PaveFiller + Builder pipeline (same as boolean ops).
    // OCCT: BOPAlgo_Splitter uses BOPAlgo_Builder (not BOPAlgo_BOP), which
    // splits all faces and classifies them per the boolean op type.
    let (brep, history) = crate::bop_occt_ops::boolean_op_with_history_generic(
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
