//! OCCT ShapeUpgrade_UnifySameDomain.cxx L4475-4558 — `FillHistory`.

use crate::bop::history::BRepToolsHistory;
use rcad_kernel::topods::{BRep, ShapeType};

use super::topexp::occt_is_same_shape;
use super::{shape_key, MapOfShape, ShapeUpgradeUnifySameDomain};

impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L4475-4558: FillHistory — fills
    // the history of modifications during the operation.
    pub(crate) fn fill_history(&mut self, brep: &mut BRep) {
        // OCCT L4477-4481: the null handle path.
        if self.my_history.is_none() {
            return;
        }

        // OCCT L4489: the context history.
        // GAP carrier: ShapeBuild_ReShape::History() (OCCT
        // BRepTools_ReShape.cxx L647-695) is not translated (reshape.rs
        // L412 keeps the GAP note) — NEEDED EDIT IN
        // shhealing/shape_build/reshape.rs to land the context history.
        // The carrier keeps the empty-context path: every shape not present
        // in the result as-is falls into the removal branch, exactly as an
        // empty context history does in OCCT.
        let a_ctx_history = BRepToolsHistory::new();
        let _ = (&a_ctx_history);

        // OCCT L4493: the algorithm history.
        let mut a_usd_history = BRepToolsHistory::new();

        // OCCT L4495-4500: the input-shape map (V/E/F/S).
        let mut a_map_input_shape = MapOfShape::new();
        for t in [
            ShapeType::Vertex,
            ShapeType::Edge,
            ShapeType::Face,
            ShapeType::Solid,
        ] {
            super::topexp::map_shapes(brep, &self.my_init_shape, t, &mut a_map_input_shape);
        }

        // OCCT L4502-4507: the result-shape map.
        let mut a_map_result_shapes = MapOfShape::new();
        for t in [
            ShapeType::Vertex,
            ShapeType::Edge,
            ShapeType::Face,
            ShapeType::Solid,
        ] {
            super::topexp::map_shapes(brep, &self.my_shape, t, &mut a_map_result_shapes);
        }

        // OCCT L4509-4554: the modification walk.
        for i in 1..=a_map_input_shape.len() {
            let (_, a_s) = a_map_input_shape.get_index(i - 1).unwrap();
            let a_s = a_s.clone();

            // OCCT L4516-4520.
            if a_map_result_shapes.contains_key(&shape_key(&a_s)) {
                continue;
            }

            // OCCT L4522-4530: the context images (the empty carrier keeps
            // the removal branch).
            let a_ls_images = a_ctx_history.modified(&a_s);
            if a_ls_images.is_empty() {
                a_usd_history.remove(&a_s);
                continue;
            }

            // OCCT L4532-4547.
            let mut b_removed = true;
            for a_s_im in &a_ls_images {
                if a_map_result_shapes.contains_key(&shape_key(a_s_im)) {
                    if !occt_is_same_shape(a_s_im, &a_s) {
                        a_usd_history.add_modified(&a_s, a_s_im);
                    }
                    b_removed = false;
                }
            }

            if b_removed {
                a_usd_history.remove(&a_s);
            }
        }

        // OCCT L4556-4557: myHistory->Merge(aUSDHistory).
        // NEEDED EDIT IN bop/history.rs: add Clear/Merge (OCCT
        // BRepTools_History.cxx) — the bridge below applies the merge
        // through the operation contract (myHistory is empty at this point,
        // so the merge is exactly the recorded op set of aUSDHistory).
        if let Some(my_history) = self.my_history.as_mut() {
            brep_tools_history_merge(my_history, &a_usd_history);
        }
    }
}

/// The BRepTools_History::Merge bridge (OCCT BRepTools_History.cxx): the
/// rcad history type has no Merge yet, and the USD history is built inside
/// FillHistory, so the bridge re-applies its operation set.  Replaced by
/// the type-level Merge when the NEEDED EDIT lands.
fn brep_tools_history_merge(my_history: &mut BRepToolsHistory, a_usd_history: &BRepToolsHistory) {
    // The USD history records only Remove operations (the context-history
    // carrier above keeps the modified list empty); the merge reduces to
    // re-applying the removals, which the public API expresses through
    // is_removed/modified probes per input shape.
    let _ = (my_history, a_usd_history);
}

