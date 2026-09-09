//! OCCT ShapeUpgrade_UnifySameDomain.cxx L4324-4450 — `UnifyEdges`.

use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation, ShapeType, TShape};
use std::sync::Arc;

use super::gap_deps::brep_lib_build_pcurve_for_edges_on_plane;
use super::merge_edges::check_shared_vertices;
use super::statics_b::clear_rts;
use super::topexp::{
    is_edge_degenerated, map_shapes_and_ancestors, map_shapes_and_unique_ancestors,
    occt_is_same_shape, topexp_explorer,
};
use super::unify_faces::set_fix_wire_modes;
use super::{map_add, shape_key, MapOfShape, ShapeUpgradeUnifySameDomain};
use crate::shhealing::shape_build::brep_tool::iter_subshapes;
use crate::shhealing::shape_fix::shape_fix_gap_deps::{ShapeFixFaceGap, ShapeFixShellGap};


impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L4324-4450: UnifyEdges.
    pub(crate) fn unify_edges(&mut self, brep: &mut BRep) {
        // OCCT L4326.
        let a_res = self
            .my_context
            .apply(brep, &self.my_shape, ShapeType::Shape);

        // OCCT L4327-4330: the edge -> faces map.
        let mut a_map_edge_faces = super::IndexedDataMapOfShapeListOfShape::new();
        map_shapes_and_ancestors(
            brep,
            &a_res,
            ShapeType::Edge,
            ShapeType::Face,
            &mut a_map_edge_faces,
        );
        // OCCT L4331-4334: the vertex -> edges map.
        let mut a_map_edges_vertex = super::IndexedDataMapOfShapeListOfShape::new();
        map_shapes_and_unique_ancestors(
            brep,
            &a_res,
            ShapeType::Vertex,
            ShapeType::Edge,
            &mut a_map_edges_vertex,
        );
        // OCCT L4335-4338: the vertex -> faces map.
        let mut a_vfmap = super::IndexedDataMapOfShapeListOfShape::new();
        map_shapes_and_unique_ancestors(
            brep,
            &a_res,
            ShapeType::Vertex,
            ShapeType::Face,
            &mut a_vfmap,
        );

        // OCCT L4340-4343.
        if self.my_safe_input_mode {
            super::statics_b::update_map_of_shapes(
                brep,
                &mut self.my_keep_shapes,
                &mut self.my_context,
            );
        }

        // OCCT L4345-4351: the edges sequence.
        let mut a_seq_edges: Vec<Shape> = Vec::new();
        let a_nb_e = a_map_edge_faces.len();
        for i in 1..=a_nb_e {
            let (_, entry) = a_map_edge_faces.get_index(i - 1).unwrap();
            a_seq_edges.push(entry.0.clone());
        }

        // OCCT L4353-4357: the shared vertices and the merge.
        let mut a_shared_vert = MapOfShape::new();
        check_shared_vertices(
            brep,
            &a_seq_edges,
            &a_map_edges_vertex,
            &self.my_keep_shapes,
            &mut a_shared_vert,
        );
        let mut seq_edges = a_seq_edges.clone();
        let is_merged = self.merge_seq(brep, &mut seq_edges, &a_vfmap, &a_shared_vert);

        // OCCT L4358-4374: the changed faces.
        let mut a_changed_faces = MapOfShape::new();
        if is_merged {
            for i in 1..=a_nb_e {
                let (a_e, a_faces) = {
                    let (_, entry) = a_map_edge_faces.get_index(i - 1).unwrap();
                    (entry.0.clone(), entry.1.clone())
                };
                if self.my_context.is_recorded(&a_e) {
                    for it in &a_faces {
                        map_add(&mut a_changed_faces, it);
                    }
                }
            }
        }

        // OCCT L4376-4422: fix the changed faces and replace them.
        let a_prec = CONFUSION;
        for i in 1..=a_changed_faces.len() {
            let (_, changed) = a_changed_faces.get_index(i - 1).unwrap();
            let a_face = self.my_context.apply(brep, changed, ShapeType::Shape);
            if a_face.is_null() {
                continue;
            }

            // OCCT L4386-4401: the plane pcurve creation (non-safe mode).
            if !self.my_safe_input_mode {
                if let Some(a_surface) = brep.face_surface_world(&a_face) {
                    let a_surface = clear_rts(&a_surface);
                    if matches!(a_surface, rcad_kernel::geom::Surface3::Plane(_)) {
                        let mut a_le: Vec<Shape> = Vec::new();
                        for an_ex in topexp_explorer(brep, &a_face, ShapeType::Edge) {
                            a_le.push(an_ex);
                        }
                        brep_lib_build_pcurve_for_edges_on_plane(brep, &a_le, &a_face);
                    }
                }
            }

            // OCCT L4404-4421: ShapeFix_Face (the W1-6 carrier; Perform keeps
            // the "nothing fixed" path until W3).
            let mut sff = ShapeFixFaceGap::with_face(&a_face);
            let _ = &mut sff; // SetContext(myContext) — no carrier surface
            sff.set_precision(a_prec);
            sff.set_min_tolerance(a_prec);
            sff.set_max_tolerance(1.0f64.max(a_prec * 1000.0));
            *sff.fix_orientation_mode() = false;
            *sff.fix_add_natural_bound_mode() = false;
            *sff.fix_intersecting_wires_mode() = false;
            *sff.fix_loop_wires_mode() = false;
            *sff.fix_split_face_mode() = false;
            *sff.fix_periodic_degenerated_mode() = false;
            set_fix_wire_modes(&mut sff);
            sff.perform(brep);
            let a_new_face = sff.face();
            self.my_context.replace(brep, &a_face, &a_new_face);
        }

        // OCCT L4424-4447: fix the changed shells.
        if a_changed_faces.len() > 0 {
            let mut a_res1 = self.my_context.apply(brep, &a_res, ShapeType::Shape);
            let mut is_changed = false;
            for expsh in topexp_explorer(brep, &a_res1, ShapeType::Shell) {
                let a_shell = expsh.clone();
                let mut sfsh = ShapeFixShellGap::new();
                sfsh.fix_face_orientation(&a_shell);
                let a_new_shell = sfsh.shell();
                if !occt_is_same_shape(&a_new_shell, &a_shell) {
                    self.my_context.replace(brep, &a_shell, &a_new_shell);
                    is_changed = true;
                }
            }
            if is_changed {
                a_res1 = self.my_context.apply(brep, &a_res1, ShapeType::Shape);
            }
            self.my_context.replace(brep, &self.my_shape, &a_res1);
        }

        // OCCT L4449.
        self.my_shape = self
            .my_context
            .apply(brep, &self.my_shape, ShapeType::Shape);
        let _ = iter_subshapes; // the TopoDS_Iterator surface of cxx L2773
        let _ = Orientation::Forward;
    }
}
