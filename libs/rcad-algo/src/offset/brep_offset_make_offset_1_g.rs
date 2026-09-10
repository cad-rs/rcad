// OCCT BRepOffset_MakeOffset_1.cxx L5924-6974 — module g of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L5924-6974:
//   - IntersectFaces, the main form (cxx L5924-6570)
//   - PrepareFacesForIntersection (cxx L6571-6679)
//   - FindVerticesToAvoid (cxx L6680-6787)
//   - FindFacesForIntersection (cxx L6788-6974)
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#49).

use std::collections::HashMap;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    add_to_container_shape, append_to_list, block_shape, empty_compound, find_common_parts,
    map_shapes_and_ancestors_map, map_shapes_indexed, shape_key_of, BRepOffsetBuildOffsetFaces,
    };
use super::brep_offset_make_offset_1_f::map_shapes_static;
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, shape_indexed_data_map, OcctIndexedShapeMap,
    OcctShapeSet, ShapeDataMap, ShapeIndexedDataMap,
};
use crate::feat::loc_ope_wires_on_shape_b::ShapeKey;

// ---------------------------------------------------------------------------
// OCCT IntersectFaces, the main form (cxx L5924-6570).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::IntersectFaces (cxx L5924-6570) —
/// intersection of the faces that should be rebuilt to resolve all
/// invalidities.
pub(crate) fn intersect_faces_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_verts_to_avoid: &mut OcctShapeSet,
    _the_range: &ProgressScope,
) {
    // OCCT L5926-5931.
    let a_nb_fr = a_bf.my_faces_to_rebuild.len();
    if a_nb_fr == 0 {
        return;
    }

    let _a_ps_outer = ProgressScope::new(&NoopProgress, "Rebuilding invalid faces", 10);

    // OCCT L5936-5948: get the vertices from the invalid edges.
    let mut a_mv_inv: OcctShapeSet = HashMap::new();
    let mut a_mv_inv_all: OcctShapeSet = HashMap::new();
    let a_nb_inv = a_bf.my_invalid_edges.extent();
    for i in 1..=a_nb_inv {
        let a_e_inv = a_bf.my_invalid_edges.find_key_1(i).clone();
        let b_valid = a_bf.my_valid_edges.contains(&a_e_inv);
        for a_v in explorer(&a_e_inv, ShapeType::Vertex, ShapeType::Shape) {
            set_add(&mut a_mv_inv_all, &a_v);
            if !b_valid {
                set_add(&mut a_mv_inv, &a_v);
            }
        }
    }

    // OCCT L5951.
    let b_look_vert_to_avoid = !a_mv_inv.is_empty();

    // OCCT L5953-5958.
    let mut a_dmsf: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let mut a_m_done: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let mut a_meinf_etrim: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let mut a_dmvefull: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let mut a_fle: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let mut a_dmefinv: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();

    // OCCT L5964-5971: PrepareFacesForIntersection.
    a_bf.prepare_faces_for_intersection(
        b_look_vert_to_avoid,
        &mut a_fle,
        &mut a_m_done,
        &mut a_dmsf,
        &mut a_meinf_etrim,
        &mut a_dmvefull,
        &mut a_dmefinv,
    );

    // OCCT L5978-5984: FindVerticesToAvoid.
    let mut a_mvr_inv: OcctShapeSet = the_verts_to_avoid.clone();
    a_bf.find_vertices_to_avoid(&a_dmefinv, &a_dmvefull, &mut a_mvr_inv);
    let _ = ProgressScope::new(&NoopProgress, "", 1);

    // OCCT L5993-6009: find the blocks of the artificially invalid faces.
    let mut a_dmf_im_f: ShapeDataMap<Shape> = HashMap::new();
    let mut a_cf_art = empty_compound();
    let art_keys: Vec<ShapeKey> = a_bf.my_art_invalid_faces.keys().cloned().collect();
    for a_key in art_keys.iter() {
        let a_f = match a_bf.my_art_invalid_faces.get(a_key) {
            Some(v) => v.0.clone(),
            None => continue,
        };
        let a_lf_inv = match a_bf.my_invalid_faces.get(&shape_key_of(&a_f)) {
            Some(v) => v.1.clone(),
            None => continue,
        };
        for a_f_im in a_lf_inv.iter() {
            add_to_container_shape(a_f_im, &mut a_cf_art);
            shape_data_map::bind(&mut a_dmf_im_f, a_f_im, a_f.clone());
        }
    }

    // OCCT L6011-6013: make the connexity blocks — the VERTEX/FACE form of
    // BOPTools_AlgoTools::MakeConnexityBlocks (cxx L6016) has no rcad
    // translation (architecture difference #44) — the GAP carrier takes
    // the OCCT empty-list path; the block list carries the face lists of
    // the blocks (the aCB list carriers).
    let a_lcb_art: Vec<Vec<Shape>> = Vec::new();

    // OCCT L6015-6016: alone edges.
    let mut a_me_alone: OcctShapeSet = HashMap::new();
    let mut a_me_inv_on_art: OcctShapeSet = HashMap::new();

    let my_oe_origins = a_bf.my_oe_origins.clone();
    let my_oe_images = a_bf.my_oe_images.clone();

    // OCCT L6018-6120: the alone-edge analysis over the artificial blocks.
    for a_cb_faces in a_lcb_art.iter() {
        let _a_ps_art = ProgressScope::new(&NoopProgress, "", 1);
        // OCCT L6022: aCB — the block (the face-list carrier).
        let a_cb = block_shape(a_cb_faces);

        // OCCT L6025-6031: check if aCB contains the splits of only one
        // offset face.
        let mut a_mf_art: OcctShapeSet = HashMap::new();
        for a_f in explorer(&a_cb, ShapeType::Face, ShapeType::Shape) {
            let a_f_origin = shape_data_map::find(&a_dmf_im_f, &a_f);
            set_add(&mut a_mf_art, &a_f_origin);
        }
        let b_alone = a_mf_art.len() == 1;

        // OCCT L6033-6043: the vertices on the invalid edges.
        let mut a_mve_inv: OcctShapeSet = HashMap::new();
        let mut a_m_fence: OcctShapeSet = HashMap::new();
        // OCCT L6038-6039: the edges that should not be marked as alone.
        let mut a_me_avoid: OcctShapeSet = HashMap::new();
        // OCCT L6040-6043: the map to find the alone edges by looking for
        // the free vertices.
        let mut a_dmve_val: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();

        for a_e in explorer(&a_cb, ShapeType::Edge, ShapeType::Shape) {
            if a_bf.my_invalid_edges.contains(&a_e) {
                set_add(&mut a_me_inv_on_art, &a_e);
                for a_v in bat::sub_shapes(&a_e) {
                    set_add(&mut a_mve_inv, &a_v);
                }
                if b_alone {
                    if let Some(p_le_or) = shape_data_map::seek(&my_oe_origins, &a_e) {
                        let p_le_or = p_le_or.clone();
                        for a_e_or in p_le_or.iter() {
                            let a_le_im = shape_data_map::find(&my_oe_images, a_e_or);
                            for a_e_im in a_le_im.iter() {
                                set_add(&mut a_me_avoid, a_e_im);
                            }
                        }
                    }
                }
                continue;
            }
            if set_add(&mut a_m_fence, &a_e) {
                map_shapes_and_ancestors_map(
                    &a_e,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                    &mut a_dmve_val,
                );
            }
        }

        // OCCT L6069-6119: find the edges with free vertices.
        let a_nb_v = a_dmve_val.len();
        for i in 1..=a_nb_v {
            let a_v = super::brep_offset_make_offset_1::find_key_1_local(&a_dmve_val, i).clone();
            if !set_contains(&a_mve_inv, &a_v) {
                continue;
            }
            let a_lev = super::brep_offset_make_offset_1::value_1_local(&a_dmve_val, i);
            if a_lev.len() > 1 {
                continue;
            }
            let a_e = a_lev[0].clone();
            if set_contains(&a_me_avoid, &a_e) {
                continue;
            }
            set_add(&mut a_me_alone, &a_e);
            // OCCT L6085-6088: if the alone edge adds nothing to the
            // intersection list the origin has been split.
            let a_e_in_dmsf = shape_data_map::find(&a_dmsf, &a_e);
            if a_e_in_dmsf.len() > 1 {
                continue;
            }
            // OCCT L6090-6100: check also its vertices.
            let mut b_broke = false;
            for a_ve in bat::sub_shapes(&a_e) {
                let a_ve_in_dmsf = shape_data_map::find(&a_dmsf, &a_ve);
                if a_ve_in_dmsf.len() > 2 {
                    b_broke = true;
                    break;
                }
            }
            if b_broke {
                continue;
            }
            // OCCT L6102-6119: the edge is useless — look for the other
            // images.
            let p_le_or = match shape_data_map::seek(&my_oe_origins, &a_e) {
                Some(v) => v.clone(),
                None => continue,
            };
            for a_e_or in p_le_or.iter() {
                let a_le_im = shape_data_map::find(&my_oe_images, a_e_or);
                for a_e_im in a_le_im.iter() {
                    if set_contains(&a_m_fence, a_e_im) {
                        set_add(&mut a_me_alone, a_e_im);
                    }
                }
            }
        }
    }

    // OCCT L6123-6141: get all the invalidities from all the faces.
    let mut a_m_all_invs: OcctShapeSet = HashMap::new();
    let a_nb_inv = a_bf.my_invalid_faces.len();
    for k in 1..=a_nb_inv {
        let a_lf = {
            let (_, v) = a_bf.my_invalid_faces.get_index(k - 1).expect("index");
            v.1.clone()
        };
        for a_f in a_lf.iter() {
            for a_e in explorer(a_f, ShapeType::Edge, ShapeType::Shape) {
                if a_bf.my_invalid_edges.contains(&a_e) || set_contains(&a_me_alone, &a_e) {
                    set_add(&mut a_m_all_invs, &a_e);
                    for a_v in bat::sub_shapes(&a_e) {
                        set_add(&mut a_m_all_invs, &a_v);
                    }
                }
            }
        }
    }

    // OCCT L6143-6156.
    let mut a_mv_bounds: OcctShapeSet = HashMap::new();
    let mut a_me_check_ext: OcctShapeSet = HashMap::new();
    let mut a_dmeetrim: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let mut a_e_images: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let mut a_dmoen_edges: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();

    // OCCT L6158-6546: the invalid-face walk.
    let a_nb_inv = a_bf.my_invalid_faces.len();
    for k in 1..=a_nb_inv {
        let (a_f_inv, a_lf_inv) = {
            let (_k2, v) = a_bf.my_invalid_faces.get_index(k - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        let b_self_reb_avoid = set_contains(&a_bf.my_f_self_reb_avoid, &a_f_inv);

        // OCCT L6167-6184: the connexity blocks of the invalid face splits.
        let mut a_lcb: Vec<Shape> = Vec::new();
        if a_lf_inv.len() > 1 {
            let mut a_cf_inv = empty_compound();
            for a_f_im in a_lf_inv.iter() {
                add_to_container_shape(a_f_im, &mut a_cf_inv);
            }
            // OCCT L6182: MakeConnexityBlocks(aCFInv, EDGE, FACE, aLCB).
            let a_cf_faces = explorer(&a_cf_inv, ShapeType::Face, ShapeType::Shape);
            let a_blocks = crate::bop::algo::shell_splitter::make_connexity_blocks(&a_cf_faces);
            for a_block in a_blocks.iter() {
                a_lcb.push(block_shape(&a_block.shapes));
            }
        } else {
            a_lcb = a_lf_inv.clone();
        }

        // OCCT L6186-6187.
        let b_artificial = shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f_inv);

        for a_cb_inv in a_lcb.iter() {
            let mut a_me_fence: OcctShapeSet = HashMap::new();
            let mut a_cbe = empty_compound();

            // OCCT L6193-6213: remember the inside edges and vertices.
            let mut an_inside_edges: OcctShapeSet = HashMap::new();
            let mut an_inside_vertices: OcctShapeSet = HashMap::new();
            for a_e in explorer(a_cb_inv, ShapeType::Edge, ShapeType::Shape) {
                if a_bf.my_invalid_edges.contains(&a_e)
                    || (b_artificial && set_contains(&a_me_alone, &a_e))
                {
                    if set_add(&mut a_me_fence, &a_e) {
                        add_to_container_shape(&a_e, &mut a_cbe);
                        if !a_bf.my_edges_to_avoid.contains(&a_e)
                            && a_bf.my_invalid_edges.contains(&a_e)
                        {
                            set_add(&mut an_inside_edges, &a_e);
                            for a_v in bat::sub_shapes(&a_e) {
                                set_add(&mut an_inside_vertices, &a_v);
                            }
                        }
                    }
                }
            }

            // OCCT L6215-6217: make the connexity blocks of the edges (the
            // VERTEX/EDGE re-host form).
            let a_cbe_edges = explorer(&a_cbe, ShapeType::Edge, ShapeType::Shape);
            let a_locations = [glam::DAffine3::IDENTITY];
            let a_lcbe = crate::bop::algo::wire_splitter::make_connexity_blocks(
                &a_cbe_edges,
                &a_locations,
            );

            for a_block in a_lcbe.iter() {
                let a_cbe_loc = block_shape(&a_block.shapes);

                // OCCT L6221-6224: the maps of the processing invalidity.
                let mut a_me = OcctIndexedShapeMap::new();
                let mut a_mecv = OcctIndexedShapeMap::new();
                map_shapes_indexed(&a_cbe_loc, ShapeType::Edge, &mut a_me);
                map_shapes_indexed(&a_cbe_loc, ShapeType::Edge, &mut a_mecv);
                map_shapes_indexed(&a_cbe_loc, ShapeType::Vertex, &mut a_me);
                let mut a_mf_int = OcctIndexedShapeMap::new();
                let mut a_mf_int_ext = OcctIndexedShapeMap::new();
                let mut a_lf_int: Vec<Shape> = Vec::new();
                let mut a_mf_avoid = OcctIndexedShapeMap::new();

                // OCCT L6229-6239: FindFacesForIntersection.
                a_bf.find_faces_for_intersection(
                    &a_f_inv,
                    &a_me,
                    &a_dmsf,
                    &a_mv_inv_all,
                    b_artificial,
                    &mut a_mf_avoid,
                    &mut a_mf_int,
                    &mut a_mf_int_ext,
                    &mut a_lf_int,
                );
                // OCCT L6240-6244: nothing to intersect.
                if a_mf_int.extent() < 3 {
                    continue;
                }

                // OCCT L6246-6249.
                let p_mf_inter = shape_data_map::seek(&a_bf.my_intersection_pairs, &a_f_inv).cloned();

                // OCCT L6251-6253.
                let mut a_me_to_int = OcctIndexedShapeMap::new();
                let a_nb = a_mf_int.extent();
                for i in 1..=a_nb {
                    let a_fi = a_mf_int.find_key_1(i).clone();
                    if b_self_reb_avoid && a_fi.is_same(&a_f_inv) {
                        continue;
                    }
                    // OCCT L6259-6262.
                    let a_lf_imi = match super::brep_offset_make_offset_1::idm_seek(&a_bf.my_of_images, &a_fi) {
                        Some(v) => v.clone(),
                        None => continue,
                    };
                    // OCCT L6264-6267.
                    let has_lfei = a_fle.contains_key(&shape_key_of(&a_fi));
                    if !has_lfei {
                        continue;
                    }
                    // OCCT L6269: aLFDone = aMDone.ChangeFind(aFi) — the
                    // snapshot form (the live list of the treated pairs).
                    let mut a_lf_done: Vec<Shape> = shape_data_map::find(&a_m_done, &a_fi);

                    // OCCT L6271-6275.
                    let p_inter_fi = match &p_mf_inter {
                        Some(m) => m.get(&shape_key_of(&a_fi)).map(|e| e.1.clone()),
                        None => None,
                    };
                    if p_mf_inter.is_some() && p_inter_fi.is_none() {
                        continue;
                    }

                    // OCCT L6277-6280: the map of the edges and vertices of
                    // aLFImi.
                    let mut a_mev_im: OcctShapeSet = HashMap::new();
                    map_shapes_static(&a_lf_imi, ShapeType::Edge, &mut a_mev_im);
                    map_shapes_static(&a_lf_imi, ShapeType::Vertex, &mut a_mev_im);

                    // OCCT L6282-6283: NCollection_MapAlgo::HasIntersection.
                    let is_i_contains_e = a_mev_im.keys().any(|k| an_inside_edges.contains_key(k));
                    let is_i_contains_v =
                        a_mev_im.keys().any(|k| an_inside_vertices.contains_key(k));

                    for j in (i + 1)..=a_nb {
                        let a_fj = a_mf_int.find_key_1(j).clone();
                        if b_self_reb_avoid && a_fj.is_same(&a_f_inv) {
                            continue;
                        }
                        // OCCT L6288-6290.
                        if let Some(p_inter_fi) = &p_inter_fi {
                            if !set_contains(p_inter_fi, &a_fj) {
                                continue;
                            }
                        }

                        // OCCT L6292-6300.
                        let a_lf_imj = match super::brep_offset_make_offset_1::idm_seek(&a_bf.my_of_images, &a_fj) {
                            Some(v) => v,
                            None => continue,
                        };
                        let has_lfej = a_fle.contains_key(&shape_key_of(&a_fj));
                        if !has_lfej {
                            continue;
                        }

                        // OCCT L6302-6306.
                        let mut a_mev_im: OcctShapeSet = HashMap::new();
                        map_shapes_static(&a_lf_imj, ShapeType::Edge, &mut a_mev_im);
                        map_shapes_static(&a_lf_imj, ShapeType::Vertex, &mut a_mev_im);
                        let is_j_contains_e =
                            a_mev_im.keys().any(|k| an_inside_edges.contains_key(k));
                        let is_j_contains_v =
                            a_mev_im.keys().any(|k| an_inside_vertices.contains_key(k));

                        // OCCT L6308-6322: check the inside-edge connection
                        // symmetry.
                        if (is_i_contains_e && !is_j_contains_v)
                            || (is_j_contains_e && !is_i_contains_v)
                        {
                            let mut a_lvc: Vec<Shape> = Vec::new();
                            find_common_parts(&a_lf_imi, &a_lf_imj, &mut a_lvc, ShapeType::Vertex);
                            if a_lvc.is_empty() {
                                continue;
                            }
                        }

                        // OCCT L6324-6327: the common edges.
                        let mut a_lec: Vec<Shape> = Vec::new();
                        find_common_parts(&a_lf_imi, &a_lf_imj, &mut a_lec, ShapeType::Edge);

                        if !a_lec.is_empty() {
                            // OCCT L6329-6349: process the common edges.
                            let b_force_use = a_mf_int_ext.contains(&a_fi)
                                || a_mf_int_ext.contains(&a_fj);
                            let mut a_lfei: Vec<Shape> =
                                super::brep_offset_make_offset_1::idm_find(&a_fle, &a_fi);
                            let mut a_lfej: Vec<Shape> =
                                super::brep_offset_make_offset_1::idm_find(&a_fle, &a_fj);
                            a_bf.process_common_edges(
                                &a_lec,
                                &a_me,
                                &a_meinf_etrim,
                                &a_m_all_invs,
                                b_force_use,
                                &mut a_mecv,
                                &mut a_me_check_ext,
                                &mut a_dmeetrim,
                                &mut a_lfei,
                                &mut a_lfej,
                                &mut a_me_to_int,
                            );
                            super::brep_offset_make_offset_1::idm_bind(&mut a_fle, &a_fi, a_lfei);
                            super::brep_offset_make_offset_1::idm_bind(&mut a_fle, &a_fj, a_lfej);

                            // OCCT L6352-6366: add the common vertices not
                            // belonging to the common edges.
                            let mut a_mv_on_ce = OcctIndexedShapeMap::new();
                            for a_e in a_lec.iter() {
                                map_shapes_indexed(a_e, ShapeType::Vertex, &mut a_mv_on_ce);
                            }
                            let mut a_lev: Vec<Shape> = Vec::new();
                            find_common_parts(
                                &a_lf_imi,
                                &a_lf_imj,
                                &mut a_lev,
                                ShapeType::Vertex,
                            );
                            for a_v in a_lev.iter() {
                                if !a_mv_on_ce.contains(a_v) {
                                    a_mecv.add(a_v);
                                }
                            }
                            continue;
                        }

                        // OCCT L6368-6374: both faces invalid and sharing
                        // edges.
                        if a_bf.my_invalid_faces.contains_key(&shape_key_of(&a_fi))
                            && a_bf.my_invalid_faces.contains_key(&shape_key_of(&a_fj))
                            && !shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_fi)
                            && !shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_fj)
                        {
                            continue;
                        }

                        // OCCT L6376-6390: check if the two faces have
                        // already been treated.
                        let mut b_treated = false;
                        for a_f in a_lf_done.iter() {
                            if a_f.is_same(&a_fj) {
                                b_treated = true;
                                break;
                            }
                        }
                        if b_treated {
                            // OCCT L6392-6401: UpdateIntersectedFaces.
                            let a_lfei = super::brep_offset_make_offset_1::idm_find(&a_fle, &a_fi);
                            let a_lfej = super::brep_offset_make_offset_1::idm_find(&a_fle, &a_fj);
                            a_bf.update_intersected_faces(
                                &a_f_inv,
                                &a_fi,
                                &a_fj,
                                &a_lf_inv,
                                &a_lf_imi,
                                &a_lf_imj,
                                &a_lfei,
                                &a_lfej,
                                &mut a_me_to_int,
                            );
                            continue;
                        }

                        // OCCT L6403-6406.
                        if a_mf_avoid.contains(&a_fi) || a_mf_avoid.contains(&a_fj) {
                            continue;
                        }

                        // OCCT L6408-6409.
                        a_lf_done.push(a_fj.clone());
                        {
                            let entry = a_m_done
                                .entry(shape_key_of(&a_fj))
                                .or_insert((a_fj.clone(), Vec::new()));
                            entry.1.push(a_fi.clone());
                        }

                        // OCCT L6411-6421: IntersectFaces (the pair form).
                        let mut a_lfei = super::brep_offset_make_offset_1::idm_find(&a_fle, &a_fi);
                        let mut a_lfej = super::brep_offset_make_offset_1::idm_find(&a_fle, &a_fj);
                        a_bf.intersect_faces_pair(
                            &a_f_inv,
                            &a_fi,
                            &a_fj,
                            &a_lf_inv,
                            &a_lf_imi,
                            &a_lf_imj,
                            &mut a_lfei,
                            &mut a_lfej,
                            &mut a_mecv,
                            &mut a_me_to_int,
                        );
                        super::brep_offset_make_offset_1::idm_bind(&mut a_fle, &a_fi, a_lfei);
                        super::brep_offset_make_offset_1::idm_bind(&mut a_fle, &a_fj, a_lfej);
                    }

                    // OCCT L6423: write back the treated list of aFi.
                    shape_data_map::bind(&mut a_m_done, &a_fi, a_lf_done);
                }

                // OCCT L6425-6438: intersect and trim the edges of this
                // chain.
                let a_ss_interfs_art: Option<ShapeDataMap<Vec<Shape>>> = if b_artificial {
                    Some(a_bf.my_ss_interfs_art.clone())
                } else {
                    None
                };
                a_bf.intersect_and_trim_edges(
                    &a_mf_int,
                    &a_me_to_int,
                    &a_dmeetrim,
                    &a_me,
                    &a_mecv,
                    &a_mv_inv,
                    &a_mvr_inv,
                    &a_me_check_ext,
                    a_ss_interfs_art.as_ref(),
                    &mut a_mv_bounds,
                    &mut a_e_images,
                );

                // OCCT L6440-6455.
                let a_nb_e_to_int = a_me_to_int.extent();
                for i_e in 1..=a_nb_e_to_int {
                    let a_e_int = a_me_to_int.find_key_1(i_e).clone();
                    for a_e in explorer(&a_cbe_loc, ShapeType::Edge, ShapeType::Shape) {
                        let entry = a_dmoen_edges
                            .entry(shape_key_of(&a_e))
                            .or_insert((a_e.clone(), Vec::new()));
                        append_to_list(&mut entry.1, &a_e_int);
                    }
                }
            }
        }
    }

    // OCCT L6548-6557: filter the obtained edges.
    a_bf.update_valid_edges(
        &a_fle,
        &a_dmoen_edges,
        &a_mv_bounds,
        &a_me_inv_on_art,
        &mut a_me_check_ext,
        the_verts_to_avoid,
        &mut a_e_images,
        &mut a_dmeetrim,
        _the_range,
    );
    let _ = a_dmvefull;
}

// ---------------------------------------------------------------------------
// OCCT PrepareFacesForIntersection (cxx L6571-6679).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::PrepareFacesForIntersection (cxx
/// L6571-6679) — preparation of the maps for analyzing the intersections
/// of the faces.
pub(crate) fn prepare_faces_for_intersection_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_look_vert_to_avoid: bool,
    the_fle: &mut ShapeIndexedDataMap<Vec<Shape>>,
    the_mdone: &mut ShapeDataMap<Vec<Shape>>,
    the_dmsf: &mut ShapeDataMap<Vec<Shape>>,
    the_meinf_etrim: &mut ShapeDataMap<Vec<Shape>>,
    the_dmvefull: &mut ShapeDataMap<Vec<Shape>>,
    the_dmefinv: &mut ShapeIndexedDataMap<Vec<Shape>>,
) {
    let my_e_trim_e_inf = a_bf.my_e_trim_e_inf.clone().unwrap_or_default();

    // OCCT L6581-6679.
    let a_nb = a_bf.my_faces_to_rebuild.len();
    for i in 1..=a_nb {
        let a_f = super::brep_offset_make_offset_1::find_key_1_local(&a_bf.my_faces_to_rebuild, i)
            .clone();

        // OCCT L6585-6586: theFLE.Add(aF, aLE); theMDone.Bind(aF, aLE).
        let a_key = shape_key_of(&a_f);
        if !the_fle.contains_key(&a_key) {
            the_fle.insert(a_key, (a_f.clone(), Vec::new()));
        }
        shape_data_map::bind(the_mdone, &a_f, Vec::new());

        // OCCT L6588-6622.
        let a_lf_im = match a_bf.my_of_images.get(&shape_key_of(&a_f)) {
            Some(v) => v.1.clone(),
            None => continue,
        };
        for a_f_im in a_lf_im.iter() {
            for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                // OCCT L6594-6598: save the connection to the untrimmed
                // face.
                {
                    let entry = the_dmsf
                        .entry(shape_key_of(&a_e))
                        .or_insert((a_e.clone(), Vec::new()));
                    append_to_list(&mut entry.1, &a_f);
                }
                // OCCT L6600-6605: save the connection to the untrimmed
                // edge.
                let a_e_inf = shape_data_map::find(&my_e_trim_e_inf, &a_e);
                {
                    let entry = the_meinf_etrim
                        .entry(shape_key_of(&a_e_inf))
                        .or_insert((a_e_inf.clone(), Vec::new()));
                    append_to_list(&mut entry.1, &a_e);
                }
                // OCCT L6607-6622.
                for a_v in explorer(&a_e, ShapeType::Vertex, ShapeType::Shape) {
                    {
                        let entry = the_dmsf
                            .entry(shape_key_of(&a_v))
                            .or_insert((a_v.clone(), Vec::new()));
                        append_to_list(&mut entry.1, &a_f);
                    }
                    if the_look_vert_to_avoid {
                        let entry = the_dmvefull
                            .entry(shape_key_of(&a_v))
                            .or_insert((a_v.clone(), Vec::new()));
                        append_to_list(&mut entry.1, &a_e);
                    }
                }
            }
        }

        // OCCT L6624-6652: get the edges of the invalid faces (from the
        // invalid splits only).
        if the_look_vert_to_avoid {
            let p_lf_inv = super::brep_offset_make_offset_1::idm_seek(&a_bf.my_invalid_faces, &a_f);
            let p_lf_inv = match p_lf_inv {
                Some(v) if !shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f) => v,
                _ => continue,
            };
            for a_f_inv in p_lf_inv.iter() {
                for a_e in explorer(a_f_inv, ShapeType::Edge, ShapeType::Shape) {
                    let entry = the_dmefinv
                        .entry(shape_key_of(&a_e))
                        .or_insert((a_e.clone(), Vec::new()));
                    append_to_list(&mut entry.1, &a_f);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT FindVerticesToAvoid (cxx L6680-6787).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FindVerticesToAvoid (cxx L6680-6787)
/// — looking for the invalid vertices.
pub(crate) fn find_vertices_to_avoid_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_dmefinv: &ShapeIndexedDataMap<Vec<Shape>>,
    the_dmvefull: &ShapeDataMap<Vec<Shape>>,
    the_mvr_inv: &mut OcctShapeSet,
) {
    let my_oe_origins = a_bf.my_oe_origins.clone();
    let my_oe_images = a_bf.my_oe_images.clone();

    // OCCT L6686.
    let mut a_m_fence: OcctShapeSet = HashMap::new();
    let a_nb = the_dmefinv.len();
    for i in 1..=a_nb {
        let (a_e, a_lf_inv) = {
            let (_k, v) = the_dmefinv.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        // OCCT L6688-6691.
        if a_lf_inv.len() == 1 {
            continue;
        }
        // OCCT L6693-6695.
        if !a_bf.my_invalid_edges.contains(&a_e) || a_bf.my_valid_edges.contains(&a_e) {
            continue;
        }
        // OCCT L6697-6699.
        if !set_add(&mut a_m_fence, &a_e) {
            continue;
        }

        // OCCT L6701-6721: the ending vertices of the images (do not check
        // the splitting vertices).
        let mut a_mve_edges: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
        match shape_data_map::seek(&my_oe_origins, &a_e) {
            Some(p_le_or) => {
                let p_le_or = p_le_or.clone();
                for a_e_or in p_le_or.iter() {
                    let a_le_im = shape_data_map::find(&my_oe_images, a_e_or);
                    for a_e_im in a_le_im.iter() {
                        set_add(&mut a_m_fence, a_e_im);
                        map_shapes_and_ancestors_map(
                            a_e_im,
                            ShapeType::Vertex,
                            ShapeType::Edge,
                            &mut a_mve_edges,
                        );
                    }
                }
            }
            None => {
                map_shapes_and_ancestors_map(
                    &a_e,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                    &mut a_mve_edges,
                );
            }
        }

        // OCCT L6723-6762.
        let a_nb_v = a_mve_edges.len();
        for j in 1..=a_nb_v {
            if shape_indexed_data_map::value_1(&a_mve_edges, j).len() != 1 {
                continue;
            }
            let a_v = super::brep_offset_make_offset_1::find_key_1_local(&a_mve_edges, j).clone();
            if !set_add(&mut a_m_fence, &a_v) {
                continue;
            }
            let p_le = match shape_data_map::seek(the_dmvefull, &a_v) {
                Some(v) => v.clone(),
                None => {
                    // OCCT L6741-6744: isolated vertex.
                    set_add(the_mvr_inv, &a_v);
                    continue;
                }
            };
            // OCCT L6746-6760.
            let mut i_nb_e_inverted = 0usize;
            let mut b_all_edges_inv = true;
            for a_ev in p_le.iter() {
                if a_bf.my_inverted_edges.contains(a_ev) {
                    i_nb_e_inverted += 1;
                }
                if b_all_edges_inv {
                    b_all_edges_inv = a_bf.my_invalid_edges.contains(a_ev);
                }
            }
            if i_nb_e_inverted > 1 || b_all_edges_inv {
                set_add(the_mvr_inv, &a_v);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT FindFacesForIntersection (cxx L6788-6974).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FindFacesForIntersection (cxx
/// L6788-6974) — looking for the faces around each invalidity for the
/// intersection.
#[allow(unused_variables)]
pub(crate) fn find_faces_for_intersection_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_f_inv: &Shape,
    the_me: &OcctIndexedShapeMap,
    the_dmsf: &ShapeDataMap<Vec<Shape>>,
    the_mv_inv_all: &OcctShapeSet,
    the_art_case: bool,
    the_mf_avoid: &mut OcctIndexedShapeMap,
    the_mf_int: &mut OcctIndexedShapeMap,
    the_mf_int_ext: &mut OcctIndexedShapeMap,
    the_lf_im_int: &mut Vec<Shape>,
) {
    let my_oe_origins = ();

    // OCCT L6789-6813.
    let a_nb_e = the_me.extent();
    let mut a_m_shapes = OcctIndexedShapeMap::new();
    for i in 1..=a_nb_e {
        let a_s = the_me.find_key_1(i).clone();
        if !shape_data_map::is_bound(the_dmsf, &a_s) {
            continue;
        }
        // OCCT L6798-6799: in the artificial case intersect the faces
        // which are close to invalidity.
        let b_avoid = if the_art_case {
            a_s.shape_type() == ShapeType::Vertex && !set_contains(the_mv_inv_all, &a_s)
        } else {
            false
        };
        let a_lf = shape_data_map::find(the_dmsf, &a_s);
        for a_f in a_lf.iter() {
            if the_mf_int.contains(a_f) {
                continue;
            }
            if b_avoid && shape_data_map::is_bound(&a_bf.my_art_invalid_faces, a_f) {
                the_mf_avoid.add(a_f);
            }
            the_mf_int.add(a_f);
            let b_use = !a_f.is_same(the_f_inv);
            let a_lf_im = match a_bf.my_of_images.get(&shape_key_of(a_f)) {
                Some(v) => v.1.clone(),
                None => continue,
            };
            for a_f_im in a_lf_im.iter() {
                the_lf_im_int.push(a_f_im.clone());
                if b_use {
                    map_shapes_indexed(a_f_im, ShapeType::Edge, &mut a_m_shapes);
                }
            }
        }
    }

    // OCCT L6815-6820.
    let a_ss_interfs_map: &ShapeDataMap<Vec<Shape>> = if the_art_case {
        &a_bf.my_ss_interfs_art
    } else {
        &a_bf.my_ss_interfs
    };
    let p_lf_inv = match shape_data_map::seek(a_ss_interfs_map, the_f_inv) {
        Some(v) => v.clone(),
        None => return,
    };

    // OCCT L6822-6830.
    let mut a_mf: OcctShapeSet = HashMap::new();
    for a_f in p_lf_inv.iter() {
        set_add(&mut a_mf, a_f);
    }

    // OCCT L6832-6846: the faces should be unique in each place.
    let mut a_cf = empty_compound();
    let mut a_mf_to_add = OcctIndexedShapeMap::new();
    let mut a_dmfor: ShapeDataMap<Shape> = HashMap::new();

    for i in 1..=a_nb_e {
        let a_s = the_me.find_key_1(i).clone();
        let p_lf = match shape_data_map::seek(a_ss_interfs_map, &a_s) {
            Some(v) => v.clone(),
            None => continue,
        };
        for a_f in p_lf.iter() {
            if the_mf_int.contains(a_f) || a_mf_to_add.contains(a_f) || !set_contains(&a_mf, a_f)
            {
                continue;
            }
            // OCCT L6856-6876: check if the face has some connection to the
            // already added for-intersection faces.
            let a_lf_im = match a_bf.my_of_images.get(&shape_key_of(a_f)) {
                Some(v) => v.1.clone(),
                None => continue,
            };
            if !the_art_case {
                let mut b_found = false;
                'outer: for a_f_im in a_lf_im.iter() {
                    for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                        if a_m_shapes.contains(&a_e) {
                            b_found = true;
                            break 'outer;
                        }
                    }
                }
                if !b_found {
                    continue;
                }
            }
            // OCCT L6878-6884.
            a_mf_to_add.add(a_f);
            for a_f_im in a_lf_im.iter() {
                shape_data_map::bind(&mut a_dmfor, a_f_im, a_f.clone());
                add_to_container_shape(a_f_im, &mut a_cf);
            }
        }
    }

    // OCCT L6886-6889.
    if a_mf_to_add.is_empty() {
        return;
    }

    // OCCT L6891-6897: MakeConnexityBlocks(aCF, EDGE, FACE, aLCB).
    let a_cf_faces = explorer(&a_cf, ShapeType::Face, ShapeType::Shape);
    let a_lcb = crate::bop::algo::shell_splitter::make_connexity_blocks(&a_cf_faces);
    if a_lcb.len() == 1 && a_mf_to_add.extent() > 1 {
        return;
    }

    // OCCT L6899-6925.
    for a_block in a_lcb.iter() {
        let a_cb = block_shape(&a_block.shapes);
        a_mf_to_add.clear();
        for a_f_im in explorer(&a_cb, ShapeType::Face, ShapeType::Shape) {
            let a_f = shape_data_map::find(&a_dmfor, &a_f_im);
            a_mf_to_add.add(&a_f);
        }
        if a_mf_to_add.extent() == 1 {
            let a_f = a_mf_to_add.find_key_1(1).clone();
            the_mf_int.add(&a_f);
            the_mf_int_ext.add(&a_f);
            let a_lf_im = match a_bf.my_of_images.get(&shape_key_of(&a_f)) {
                Some(v) => v.1.clone(),
                None => continue,
            };
            for a_f_im in a_lf_im.iter() {
                the_lf_im_int.push(a_f_im.clone());
            }
        }
    }
}
