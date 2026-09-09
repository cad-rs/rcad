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

        // OCCT L4489: the context history — all modifications recorded by
        // the reshaper during the operation.
        let a_ctx_history = self.my_context.history();

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
        if let Some(my_history) = self.my_history.as_mut() {
            my_history.merge(&a_usd_history);
        }
    }
}

