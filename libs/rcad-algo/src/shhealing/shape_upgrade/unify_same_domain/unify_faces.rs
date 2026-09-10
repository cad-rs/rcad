//! OCCT ShapeUpgrade_UnifySameDomain.cxx L3047-3181 — `UnifyFaces`
//! (L3047-3129), static `SetFixWireModes` (L3133-3144), static
//! `isSameSets` (L3146-3181).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType};

use super::topexp::{
    map_shapes_and_ancestors, topexp_explorer, topexp_explorer_avoid,
};
use super::{
    map_add, shape_key, DataMapOfShapeMapOfShape, IndexedDataMapOfShapeListOfShape, MapOfShape,
    ShapeUpgradeUnifySameDomain,
};
use crate::shhealing::shape_build::brep_tool::iter_subshapes;

impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L3047-3129: UnifyFaces.
    pub(crate) fn unify_faces(&mut self, brep: &mut BRep) {
        // OCCT L3049-3058: the global edge -> faces map.
        let mut a_gmap_edge_faces = IndexedDataMapOfShapeListOfShape::new();
        let mut a_face_map = MapOfShape::new();
        super::topexp::map_shapes(brep, &self.my_shape, ShapeType::Face, &mut a_face_map);
        for i in 1..=a_face_map.len() {
            let (_, a_face) = a_face_map.get_index(i - 1).unwrap();
            map_shapes_and_ancestors(
                brep,
                a_face,
                ShapeType::Edge,
                ShapeType::Face,
                &mut a_gmap_edge_faces,
            );
        }

        // OCCT L3060-3080: the face -> shells map (the unification of faces
        // from different shells is avoided).
        let mut a_gmap_face_shells = DataMapOfShapeMapOfShape::new();
        for an_exp in topexp_explorer(brep, &self.my_shape, ShapeType::Shell) {
            let a_shell = an_exp.clone();
            for an_it_f in iter_subshapes(brep, &a_shell, true, true) {
                let a_f = an_it_f;
                let key = shape_key(&a_f);
                match a_gmap_face_shells.get_mut(&key) {
                    Some((_, shells)) => {
                        map_add(shells, &a_shell);
                    }
                    None => {
                        let mut shells = MapOfShape::new();
                        map_add(&mut shells, &a_shell);
                        a_gmap_face_shells.insert(key, (a_f, shells));
                    }
                }
            }
        }

        // OCCT L3082-3103: the free-boundary map (only shells not belonging
        // to solids).
        let mut a_free_bound_map = MapOfShape::new();
        for a_shell in
            topexp_explorer_avoid(brep, &self.my_shape, ShapeType::Shell, ShapeType::Solid)
        {
            let mut a_efmap = IndexedDataMapOfShapeListOfShape::new();
            map_shapes_and_ancestors(
                brep,
                &a_shell,
                ShapeType::Edge,
                ShapeType::Face,
                &mut a_efmap,
            );
            for ii in 1..=a_efmap.len() {
                let (an_edge, a_face_list) = {
                    let (_, (s, l)) = a_efmap.get_index(ii - 1).unwrap();
                    (s.clone(), l.len())
                };
                if !super::topexp::is_edge_degenerated(brep, &an_edge) && a_face_list == 1 {
                    map_add(&mut a_free_bound_map, &an_edge);
                }
            }
        }

        // OCCT L3105-3110: unify faces in each shell separately.
        for exps in topexp_explorer(brep, &self.my_shape, ShapeType::Shell) {
            self.int_unify_faces(
                brep,
                &exps,
                &a_gmap_edge_faces,
                &a_gmap_face_shells,
                &a_free_bound_map,
            );
        }

        // OCCT L3112-3126: gather all faces out of shells in one compound
        // and unify them at once.
        let mut a_bb = BRepBuilder::new();
        let a_cmp = a_bb.make_compound(brep, vec![]);
        let mut nbf = 0i32;
        for exps in topexp_explorer_avoid(brep, &self.my_shape, ShapeType::Face, ShapeType::Shell) {
            a_bb.add_to_compound(brep, a_cmp.clone(), exps.clone());
            nbf += 1;
        }

        if nbf > 0 {
            // OCCT L3124-3125: no connection to shells — the empty
            // face-shell map.
            self.int_unify_faces(
                brep,
                &a_cmp,
                &a_gmap_edge_faces,
                &DataMapOfShapeMapOfShape::new(),
                &a_free_bound_map,
            );
        }

        // OCCT L3128.
        self.my_shape = self
            .my_context
            .apply(brep, &self.my_shape, ShapeType::Shape);
    }
}

/// OCCT static SetFixWireModes (cxx L3133-3144): the ShapeFix_Wire tool
/// flags of the ShapeFix_Face tool (the real W3 tranche 2 `ShapeFixWire`).
pub fn set_fix_wire_modes(
    the_sff: &mut crate::shhealing::shape_fix::face_a::ShapeFixFace,
) {
    let a_fix_wire = the_sff.fix_wire_tool();
    *a_fix_wire.fix_self_intersection_mode() = 0;
    *a_fix_wire.fix_non_adjacent_intersecting_edges_mode() = 0;
    *a_fix_wire.fix_lacking_mode() = 0;
    *a_fix_wire.fix_notched_edges_mode() = 0;
    *a_fix_wire.modify_topology_mode() = false;
    *a_fix_wire.modify_remove_loop_mode() = 0;
    *a_fix_wire.fix_gaps_by_ranges_mode() = false;
    *a_fix_wire.fix_small_mode() = 0;
}

/// OCCT static isSameSets (cxx L3146-3181): compares two sets of shapes
/// (the template maps over the face-shell map values; the None case is the
/// OCCT null pointer).
pub fn is_same_sets(
    the_fshells1: Option<&(Shape, MapOfShape)>,
    the_fshells2: Option<&(Shape, MapOfShape)>,
) -> bool {
    // OCCT L3156-3159: both null — no problem.
    let (Some((_, s1)), Some((_, s2))) = (the_fshells1, the_fshells2) else {
        return the_fshells1.is_none() && the_fshells2.is_none();
    };
    // OCCT L3166-3169.
    if s1.len() != s2.len() {
        return false;
    }
    // OCCT L3170-3180: the mutual containment (the sets are small).
    for (_, v1) in s1.iter() {
        if !s2.contains_key(&shape_key(v1)) {
            return false;
        }
    }
    for (_, v2) in s2.iter() {
        if !s1.contains_key(&shape_key(v2)) {
            return false;
        }
    }
    true
}
