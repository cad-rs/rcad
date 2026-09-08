// OCCT BRepOffset_MakeOffset_1.cxx L7607-8550 — module i of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L7607-8550:
//   - GetInvalidEdges (cxx L7607-7678)
//   - UpdateValidEdges (cxx L7679-8238)
//   - TrimNewIntersectionEdges (cxx L8239-8407)
//   - IntersectEdges (cxx L8408-8550)
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#49).

use std::collections::HashMap;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::builder::Builder;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    add_to_container_map_list, add_to_container_shape, append_to_list, builder_modified,
    empty_compound, map_shapes_and_ancestors_map, map_shapes_indexed, shape_key_of,
    update_origins, update_images, BRepOffsetBuildOffsetFaces, BuilderRef,
};
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, OcctIndexedShapeMap, OcctShapeSet, ShapeDataMap,
    ShapeIndexedDataMap, top_exp_vertices,
};

// ---------------------------------------------------------------------------
// OCCT GetInvalidEdges (cxx L7607-7678).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::GetInvalidEdges (cxx L7607-7678) —
/// looking for the invalid edges by intersecting with the invalid
/// vertices.
pub(crate) fn get_invalid_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_verts_to_avoid: &OcctShapeSet,
    the_mv_bounds: &OcctShapeSet,
    the_gf: &BuilderRef,
    the_me_inv: &mut OcctShapeSet,
) {
    let _ = a_bf;
    // OCCT L7614-7617.
    if the_verts_to_avoid.is_empty() {
        return;
    }
    // OCCT L7619-7624: get the vertices created with the intersection
    // edges — aRes = theGF.Shape().
    let mut a_dmve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    if let Some(a_res) = the_gf.shape() {
        map_shapes_and_ancestors_map(&a_res, ShapeType::Vertex, ShapeType::Edge, &mut a_dmve);
    }
    // OCCT L7626: pDS = theGF.PDS() — the DS access for the IsNewShape
    // probes (the Builder form; the MakerVolume form takes the OCCT
    // nV < 0 path — the shape is not in the DS).
    let p_ds_index = |a_v: &Shape| -> isize {
        match the_gf {
            BuilderRef::Builder(b) => b.ds.index(a_v),
            BuilderRef::MakerVolume(_) => -1,
        }
    };
    let p_ds_is_new_shape = |n_v: isize| -> bool {
        match the_gf {
            BuilderRef::Builder(b) => n_v >= 0 && b.ds.is_new_shape(n_v as usize),
            BuilderRef::MakerVolume(_) => false,
        }
    };

    // OCCT L7634-7672: find the invalid splits of the edges.
    let mut a_mv_inv: OcctShapeSet = HashMap::new();
    let a_nb = a_dmve.len();
    for i in 1..=a_nb {
        let (a_v, a_lve) = {
            let (_k, v) = a_dmve.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        if set_contains(the_mv_bounds, &a_v) {
            continue;
        }
        // OCCT L7638-7641.
        let n_v = p_ds_index(&a_v);
        if n_v >= 0 && !p_ds_is_new_shape(n_v) {
            continue;
        }
        // OCCT L7643-7653: check the vertex SD with the vertices to avoid
        // (BOPTools_AlgoTools::ComputeVV).
        let mut b_found = false;
        for a_v_inv in the_verts_to_avoid.values() {
            let a_tol = bat::brep_tool_tolerance(&a_v);
            let a_tol_inv = bat::brep_tool_tolerance(a_v_inv);
            let a_p = bat::brep_tool_pnt(&a_v).unwrap_or(glam::DVec3::ZERO);
            let a_p_inv = bat::brep_tool_pnt(a_v_inv).unwrap_or(glam::DVec3::ZERO);
            let i_flag =
                crate::bop::tools::algo_tools::compute_vv(a_tol, a_p, a_tol_inv, a_p_inv, 0.0);
            if i_flag == 0 {
                set_add(&mut a_mv_inv, &a_v);
                b_found = true;
                break;
            }
        }
        // OCCT L7656-7666.
        if b_found {
            for a_e in a_lve.iter() {
                set_add(the_me_inv, a_e);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT UpdateValidEdges (cxx L7679-8238).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::UpdateValidEdges (cxx L7679-8238) —
/// making the new splits and updating the maps.
#[allow(clippy::too_many_arguments)]
pub(crate) fn update_valid_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_fle: &ShapeIndexedDataMap<Vec<Shape>>,
    the_oen_edges: &ShapeIndexedDataMap<Vec<Shape>>,
    the_mv_bounds: &OcctShapeSet,
    the_me_inv_on_art: &OcctShapeSet,
    the_me_check_ext: &mut OcctShapeSet,
    the_verts_to_avoid: &mut OcctShapeSet,
    the_e_images: &mut ShapeDataMap<Vec<Shape>>,
    the_eetrim: &mut ShapeDataMap<Vec<Shape>>,
    _the_range: &ProgressScope,
) {
    let _a_ps_outer = ProgressScope::new(&NoopProgress, "Updating edges", 10);

    // OCCT L7687-7689: the new edges and the back connection from the
    // edges to the faces.
    let mut a_le: Vec<Shape> = Vec::new();
    let mut a_melf: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let mut a_me_tmp: OcctShapeSet = HashMap::new();

    // OCCT L7691-7713.
    let a_nb = the_fle.len();
    for i in 1..=a_nb {
        let (a_f, a_le_int) = {
            let (_k, v) = the_fle.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        for a_e in a_le_int.iter() {
            // OCCT L7697-7702.
            if (set_contains(the_me_check_ext, a_e) || set_contains(&a_me_tmp, a_e))
                && !shape_data_map::is_bound(the_e_images, a_e)
            {
                the_me_check_ext.remove(&shape_key_of(a_e));
                set_add(&mut a_me_tmp, a_e);
                continue;
            }
            // OCCT L7703-7711.
            let is_new = !shape_data_map::is_bound(&a_melf, a_e);
            add_to_container_map_list(&mut a_melf, a_e, &a_f);
            if is_new {
                a_le.push(a_e.clone());
            }
        }
    }

    // OCCT L7715-7718.
    if a_le.is_empty() {
        return;
    }

    // OCCT L7720-7730.
    let mut a_meb: OcctShapeSet = HashMap::new();
    let mut a_me_new: OcctShapeSet = HashMap::new();
    let mut a_mv_old: OcctShapeSet = HashMap::new();
    let mut a_dme_or: ShapeDataMap<Vec<Shape>> = HashMap::new();

    // OCCT L7732-7742: trim the new intersection edges.
    a_bf.trim_new_intersection_edges(
        &a_le,
        the_eetrim,
        the_mv_bounds,
        the_me_check_ext,
        the_e_images,
        &mut a_meb,
        &mut a_mv_old,
        &mut a_me_new,
        &mut a_dme_or,
        &mut a_melf,
    );

    // OCCT L7744-7750.
    if the_e_images.is_empty() {
        // OCCT L7747-7749: no new splits preserved — update the
        // intersection edges and exit.
        a_bf.update_new_intersection_edges(&a_le, &a_melf, the_e_images, the_eetrim);
        return;
    }

    let _ = ProgressScope::new(&NoopProgress, "", 1);

    // OCCT L7752-7770: the compound of all the invalid edges — aCEAll.
    let mut a_ce_all = empty_compound();
    let a_nb_e = the_oen_edges.len();
    for i in 1..=a_nb_e {
        let a_e = super::brep_offset_make_offset_1::find_key_1_local(the_oen_edges, i).clone();
        add_to_container_shape(&a_e, &mut a_ce_all);
    }

    // OCCT L7772-7774: separate the edges into blocks (the VERTEX/EDGE
    // re-host form).
    let a_ce_all_edges = explorer(&a_ce_all, ShapeType::Edge, ShapeType::Shape);
    let a_locations = [glam::DAffine3::IDENTITY];
    let a_l_blocks =
        crate::bop::algo::wire_splitter::make_connexity_blocks(&a_ce_all_edges, &a_locations);

    // OCCT L7776-7779: the intersected splits.
    let mut a_m_blocks_sp: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();

    // OCCT L7781-7866: the block walk.
    for a_block in a_l_blocks.iter() {
        let _a_psb = ProgressScope::new(&NoopProgress, "", 1);
        let a_block = block_shape_of(&a_block.shapes);

        // OCCT L7785-7801: the fence and the block's new edges.
        let mut a_block_le_new: Vec<Shape> = Vec::new();
        {
            let mut a_me_fence: OcctShapeSet = HashMap::new();
            for a_e in explorer(&a_block, ShapeType::Edge, ShapeType::Shape) {
                let a_le_int = match super::brep_offset_make_offset_1::idm_seek(the_oen_edges, &a_e) {
                    Some(v) => v,
                    None => continue,
                };
                for a_e_int in a_le_int.iter() {
                    if set_add(&mut a_me_fence, a_e_int) {
                        a_block_le_new.push(a_e_int.clone());
                    }
                }
            }
        }

        // OCCT L7803-7806.
        if a_block_le_new.is_empty() {
            continue;
        }

        // OCCT L7808-7828: get the splits of the new edges to intersect.
        let mut a_lsplits: Vec<Shape> = Vec::new();
        for a_e in a_block_le_new.iter() {
            let p_le_im = shape_data_map::seek(the_e_images, a_e).cloned();
            let p_le_im = match p_le_im {
                Some(v) if !v.is_empty() => v,
                _ => continue,
            };
            for a_e_im in p_le_im.iter() {
                a_lsplits.push(a_e_im.clone());
            }
        }
        if a_lsplits.is_empty() {
            continue;
        }

        // OCCT L7830-7849.
        let a_ce: Shape = if a_lsplits.len() > 1 {
            // OCCT L7834-7843: intersect the new splits among themselves.
            let mut a_ce = Shape::null();
            a_bf.intersect_edges(
                &a_lsplits,
                &a_block_le_new,
                the_mv_bounds,
                the_verts_to_avoid,
                &mut a_me_new,
                the_me_check_ext,
                the_e_images,
                &mut a_dme_or,
                &mut a_melf,
                &mut a_ce,
            );
            a_ce
        } else {
            a_lsplits[0].clone()
        };

        // OCCT L7851: aMBlocksSp.Add(aCE, aBlockLENew).
        let key = shape_key_of(&a_ce);
        if let Some(entry) = a_m_blocks_sp.get_mut(&key) {
            entry.1 = a_block_le_new.clone();
        } else {
            a_m_blocks_sp.insert(key, (a_ce.clone(), a_block_le_new.clone()));
        }
    }

    // OCCT L7854-7860: the first stage — the separate treatment of the
    // blocks.
    let mut a_me_val: OcctShapeSet = HashMap::new();
    let mut a_l_val_blocks: Vec<Shape> = Vec::new();

    let a_nb_b = a_m_blocks_sp.len();
    for i in 1..=a_nb_b {
        let _a_psb_sp = ProgressScope::new(&NoopProgress, "", 1);
        let (a_ce, a_block_le_new) = {
            let (_k, v) = a_m_blocks_sp.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };

        // OCCT L7867-7882: get all the participating faces for the bounds.
        let mut a_lfaces: Vec<Shape> = Vec::new();
        for a_e in a_block_le_new.iter() {
            let p_lf = match shape_data_map::seek(&a_melf, a_e) {
                Some(v) => v.clone(),
                None => continue,
            };
            for a_f in p_lf.iter() {
                append_to_list(&mut a_lfaces, a_f);
            }
        }

        // OCCT L7884-7886: the localized bounds.
        let mut a_filter_bounds = Shape::null();
        a_bf.get_bounds(&a_lfaces, &a_meb, &mut a_filter_bounds);

        // OCCT L7888-7897: filter the splits by the bounds.
        let mut a_me_inv_loc: OcctShapeSet = HashMap::new();
        a_bf.get_invalid_edges_by_bounds(
            &a_ce,
            &a_filter_bounds,
            &a_mv_old,
            &a_me_new,
            &a_dme_or,
            &a_melf,
            the_e_images,
            the_me_check_ext,
            the_me_inv_on_art,
            the_verts_to_avoid,
            &mut a_me_inv_loc,
        );

        // OCCT L7899-7912: keep only the valid edges of the block.
        let mut a_ce_val = empty_compound();
        let mut b_kept = false;
        for a_esp in explorer(&a_ce, ShapeType::Edge, ShapeType::Shape) {
            if !set_contains(&a_me_inv_loc, &a_esp) && set_add(&mut a_me_val, &a_esp) {
                add_to_container_shape(&a_esp, &mut a_ce_val);
                b_kept = true;
            }
        }
        if b_kept {
            a_l_val_blocks.push(a_ce_val);
        }
    }

    // OCCT L7914-7916: filter the images of the edges after the first
    // stage.
    let mut a_splits1 = Shape::null();
    a_bf.filter_splits(&a_le, &a_me_val, false, the_e_images, &mut a_splits1);

    // OCCT L7918-7924.
    if a_l_val_blocks.is_empty() {
        a_bf.update_new_intersection_edges(&a_le, &a_melf, the_e_images, the_eetrim);
        return;
    }

    let _ = ProgressScope::new(&NoopProgress, "", 1);

    // OCCT L7926-7940: the second stage — add the already removed new
    // edges as the markers.
    let a_nb_b = a_m_blocks_sp.len();
    for i in 1..=a_nb_b {
        let a_ce = {
            let (_, v) = a_m_blocks_sp.get_index(i - 1).expect("index");
            v.0.clone()
        };
        for a_e_im in explorer(&a_ce, ShapeType::Edge, ShapeType::Shape) {
            if set_contains(&a_me_new, &a_e_im) && !set_contains(&a_me_val, &a_e_im) {
                a_l_val_blocks.push(a_e_im.clone());
            }
        }
    }

    // OCCT L7942-7956.
    if a_l_val_blocks.len() > 1 {
        // OCCT L7946-7954: intersect the new splits among themselves.
        a_bf.intersect_edges(
            &a_l_val_blocks,
            &a_le,
            the_mv_bounds,
            the_verts_to_avoid,
            &mut a_me_new,
            the_me_check_ext,
            the_e_images,
            &mut a_dme_or,
            &mut a_melf,
            &mut a_splits1,
        );
    } else {
        a_splits1 = a_l_val_blocks[0].clone();
    }

    let _ = ProgressScope::new(&NoopProgress, "", 1);

    // OCCT L7958-7966: get all the faces for the bounds.
    let mut a_lfaces: Vec<Shape> = Vec::new();
    let a_nb = a_bf.my_of_images.len();
    for i in 1..=a_nb {
        let a_f = super::brep_offset_make_offset_1::find_key_1_local(&a_bf.my_of_images, i)
            .clone();
        a_lfaces.push(a_f);
    }

    // OCCT L7968-7970.
    let mut a_filter_bounds = Shape::null();
    a_bf.get_bounds(&a_lfaces, &a_meb, &mut a_filter_bounds);

    // OCCT L7972-7984.
    let mut a_me_inv: OcctShapeSet = HashMap::new();
    a_bf.get_invalid_edges_by_bounds(
        &a_splits1,
        &a_filter_bounds,
        &a_mv_old,
        &a_me_new,
        &a_dme_or,
        &a_melf,
        the_e_images,
        the_me_check_ext,
        the_me_inv_on_art,
        the_verts_to_avoid,
        &mut a_me_inv,
    );

    // OCCT L7986-7990.
    let mut a_splits = Shape::null();
    a_bf.filter_splits(&a_le, &a_me_inv, true, the_e_images, &mut a_splits);

    let _ = ProgressScope::new(&NoopProgress, "", 1);

    // OCCT L7992-8026: get the bounds to update.
    let mut a_lf: Vec<Shape> = Vec::new();
    let mut a_mv_sp = OcctIndexedShapeMap::new();
    map_shapes_indexed(&a_splits, ShapeType::Vertex, &mut a_mv_sp);
    let a_nb_f = a_bf.my_of_images.len();
    for i in 1..=a_nb_f {
        let (a_f, a_lf_im) = {
            let (_k, v) = a_bf.my_of_images.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        if the_fle.contains_key(&shape_key_of(&a_f)) {
            a_lf.push(a_f);
            continue;
        }
        // OCCT L8007-8022: check the splits of the faces to have the
        // vertices from the splits.
        let mut b_found = false;
        'outer: for a_f_im in a_lf_im.iter() {
            for a_v in explorer(a_f_im, ShapeType::Vertex, ShapeType::Shape) {
                if a_mv_sp.contains(&a_v) {
                    b_found = true;
                    break 'outer;
                }
            }
        }
        if b_found {
            a_lf.push(a_f);
        }
    }

    // OCCT L8028-8030.
    let mut a_bounds = Shape::null();
    let mut a_la_valid: Vec<Shape> = Vec::new();
    let mut a_la_bounds: Vec<Shape> = Vec::new();
    a_bf.get_bounds_to_update(&a_lf, &a_meb, &mut a_la_bounds, &mut a_la_valid, &mut a_bounds);

    // OCCT L8032-8037: intersect the valid splits with the bounds.
    let mut a_gf_args: Vec<Shape> = Vec::new();
    a_gf_args.push(a_bounds.clone());
    a_gf_args.push(a_splits.clone());
    let mut a_gf_filler = PaveFiller::new();
    a_gf_filler.set_arguments(a_gf_args);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_Builder", 1);
    a_gf_filler.perform(&a_ps);
    let mut a_gf = Builder::new(
        a_gf_filler.ds(),
        crate::bop::algo::builder::BooleanOpType::Union,
        a_gf_filler.fuzzy_value(),
    );
    a_gf.my_arguments = a_gf_filler.ds().arguments.clone();
    let _a_root = a_gf
        .build()
        .ok()
        .map(|brep| super::brep_offset_make_offset_1::brep_root_shape(&brep));

    // OCCT L8039: update the splits — UpdateImages(aLE, theEImages, aGF,
    // myModifiedEdges).
    {
        let mut my_modified_edges = a_bf.my_modified_edges.clone();
        let mut e_images = the_e_images.clone();
        update_images(&a_le, &mut e_images, &a_gf, &mut my_modified_edges);
        *the_e_images = e_images;
        a_bf.my_modified_edges = my_modified_edges;
    }

    // OCCT L8041: update the new intersection edges.
    a_bf.update_new_intersection_edges(&a_le, &a_melf, the_e_images, the_eetrim);

    // OCCT L8043-8047: update the bounds.
    {
        let mut my_modified_edges = a_bf.my_modified_edges.clone();
        let mut my_oe_images = a_bf.my_oe_images.clone();
        let mut my_oe_origins = a_bf.my_oe_origins.clone();
        let mut my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();
        update_images(&a_la_valid, &mut my_oe_images, &a_gf, &mut my_modified_edges);
        update_origins(&a_la_bounds, &mut my_oe_origins, &a_gf);
        update_origins(&a_la_bounds, &mut my_edges_origins, &a_gf);
        super::brep_offset_make_offset_1_b::update_intersected_edges(
            a_bf,
            &a_la_bounds,
            &a_gf,
        );
        a_bf.my_modified_edges = my_modified_edges;
        a_bf.my_oe_images = my_oe_images;
        a_bf.my_oe_origins = my_oe_origins;
        a_bf.my_edges_origins = Some(my_edges_origins);
    }

    // OCCT L8049-8062: update the EdgesToAvoid with the splits.
    let mut a_new_edges = OcctIndexedShapeMap::new();
    {
        let p_splits_im = a_gf
            .my_images
            .get((a_splits.ptr_id(), a_splits.location))
            .cloned();
        if let Some(p_splits_im) = p_splits_im {
            for a_sp_im in p_splits_im.iter() {
                map_shapes_indexed(a_sp_im, ShapeType::Edge, &mut a_new_edges);
            }
        }
    }

    // OCCT L8064-8068.
    let mut an_inside_edges = empty_compound();
    for i_e in 1..=a_bf.my_inside_edges.extent() {
        let a_e = a_bf.my_inside_edges.find_key_1(i_e).clone();
        add_to_container_shape(&a_e, &mut an_inside_edges);
    }

    // OCCT L8070-8074: rebuild the map of the edges to avoid.
    let mut a_me_avoid = OcctIndexedShapeMap::new();
    let mut a_ce_avoid = empty_compound();
    let p_ds = a_gf.ds;

    let a_nb_e = a_bf.my_edges_to_avoid.extent();
    for i in 1..=a_nb_e {
        let a_e = a_bf.my_edges_to_avoid.find_key_1(i).clone();
        let a_le_im = builder_modified(&a_gf, &a_e);

        // OCCT L8082-8105: only the untouched and fully coinciding edges
        // should be kept in the avoid map.
        let mut b_keep = a_le_im.is_empty();
        if a_le_im.len() == 1 && a_e.is_same(&a_le_im[0]) {
            let n_e = p_ds.index(&a_e);
            if n_e >= 0 {
                let a_lpb = p_ds.pave_blocks(n_e as usize);
                if a_lpb.len() == 1 {
                    let a_pb = &a_lpb[0];
                    if let Some(cb_idx) = p_ds.common_block(a_pb) {
                        let a_lpbc = p_ds.common_blocks[cb_idx].pave_blocks();
                        let mut b_broke = false;
                        for (a_pb_cb, _) in a_lpbc.iter() {
                            let n_e_orig: usize = a_pb_cb.read().original_edge();
                            if p_ds.pave_blocks(n_e_orig).len() > 1 {
                                b_broke = true;
                                break;
                            }
                        }
                        b_keep = !b_broke;
                    }
                }
            }
        }

        if b_keep {
            // OCCT L8107-8110: keep the original edge.
            a_me_avoid.add(&a_e);
            continue;
        }

        // OCCT L8112-8119.
        for a_e_im in a_le_im.iter() {
            if !a_new_edges.contains(a_e_im) {
                add_to_container_shape(a_e_im, &mut a_ce_avoid);
            }
        }
    }

    // OCCT L8121-8140: the BOPAlgo_BOP CUT form (the rcad Builder with the
    // Cut operation; IsDeleted takes the OCCT false path — the
    // architecture difference #49 query surface).
    let mut is_cut = false;
    if !bat::sub_shapes(&a_ce_avoid).is_empty() {
        // OCCT L8122-8128: BOPAlgo_BOP aBOP; AddArgument(aCEAvoid);
        // AddTool(anInsideEdges); SetOperation(BOPAlgo_CUT); Perform().
        let mut a_bop_args: Vec<Shape> = Vec::new();
        a_bop_args.push(a_ce_avoid.clone());
        a_bop_args.push(an_inside_edges.clone());
        let mut a_bop_filler = PaveFiller::new();
        a_bop_filler.set_arguments(a_bop_args);
        let a_prog = NoopProgress;
        let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_BOP", 1);
        a_bop_filler.perform(&a_ps);
        let mut a_bop = Builder::new(
            a_bop_filler.ds(),
            crate::bop::algo::builder::BooleanOpType::Cut,
            a_bop_filler.fuzzy_value(),
        );
        a_bop.my_arguments = a_bop_filler.ds().arguments.clone();
        a_bop.my_tools = a_bop_filler.ds().arguments[1..].to_vec();
        let a_bop_root = a_bop
            .build()
            .ok()
            .map(|brep| super::brep_offset_make_offset_1::brep_root_shape(&brep));
        is_cut = a_bop_root.is_some() && !a_bop.has_errors();
        if is_cut {
            for a_e in bat::sub_shapes(&a_ce_avoid) {
                // OCCT L8129-8132: if (!aBOP.IsDeleted(itCE.Value())).
                if !a_bop_is_deleted_gap(&a_e) {
                    a_me_avoid.add(&a_e);
                }
            }
        }
    }

    if !is_cut {
        // OCCT L8136-8138.
        map_shapes_indexed(&a_ce_avoid, ShapeType::Edge, &mut a_me_avoid);
    }

    // OCCT L8140: myEdgesToAvoid = aMEAvoid.
    a_bf.my_edges_to_avoid = a_me_avoid;
}

/// The block-shape carrier of the connexity-block lists.
fn block_shape_of(the_list: &[Shape]) -> Shape {
    super::brep_offset_make_offset_1::block_shape(the_list)
}

/// OCCT BOPAlgo_BOP::IsDeleted(S) — the GAP form (the rcad Builder keeps
/// no deleted-shapes table) — the OCCT false path.
fn a_bop_is_deleted_gap(_the_s: &Shape) -> bool {
    false
}

// ---------------------------------------------------------------------------
// OCCT TrimNewIntersectionEdges (cxx L8239-8407).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::TrimNewIntersectionEdges (cxx
/// L8239-8407) — trims the intersection edges.
#[allow(clippy::too_many_arguments)]
pub(crate) fn trim_new_intersection_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_le: &[Shape],
    the_eetrim: &ShapeDataMap<Vec<Shape>>,
    the_mv_bounds: &OcctShapeSet,
    the_mecheckext: &mut OcctShapeSet,
    the_e_images: &mut ShapeDataMap<Vec<Shape>>,
    the_meb: &mut OcctShapeSet,
    the_mv_old: &mut OcctShapeSet,
    the_me_new: &mut OcctShapeSet,
    the_dme_or: &mut ShapeDataMap<Vec<Shape>>,
    the_melf: &mut ShapeDataMap<Vec<Shape>>,
) {
    // OCCT L8245-8404.
    for a_e in the_le.iter() {
        // OCCT L8250: bCheckExt = theMECheckExt.Remove(aE).
        let b_check_ext = set_contains(the_mecheckext, a_e);
        if b_check_ext {
            the_mecheckext.remove(&shape_key_of(a_e));
        }

        // OCCT L8252-8267: bOld = theEETrim.IsBound(aE).
        let b_old = shape_data_map::is_bound(the_eetrim, a_e);
        if b_old {
            let a_let = shape_data_map::find(the_eetrim, a_e);
            for a_et in a_let.iter() {
                set_add(the_meb, a_et);
                for a_v in explorer(a_et, ShapeType::Vertex, ShapeType::Shape) {
                    set_add(the_mv_old, &a_v);
                }
            }
        }

        // OCCT L8269-8275.
        if !shape_data_map::is_bound(the_e_images, a_e) {
            continue;
        }
        let a_le_im = shape_data_map::find(the_e_images, a_e);
        if a_le_im.is_empty() {
            shape_data_map::un_bind(the_e_images, a_e);
            continue;
        }

        let mut a_ce_im = Shape::null();
        let mut a_mev_bounds: OcctShapeSet = HashMap::new();

        if a_le_im.len() > 1 {
            // OCCT L8283-8306: fuse these parts.
            let mut a_mv = OcctIndexedShapeMap::new();
            let mut a_gfe_args: Vec<Shape> = Vec::new();
            for a_e_im in a_le_im.iter() {
                a_gfe_args.push(a_e_im.clone());
                map_shapes_indexed(a_e_im, ShapeType::Vertex, &mut a_mv);
            }
            // OCCT L8296-8301: add the two bounding vertices of this edge.
            let (a_v1, a_v2) = top_exp_vertices(a_e);
            a_gfe_args.push(a_v1.clone());
            a_gfe_args.push(a_v2.clone());
            a_mv.add(&a_v1);
            a_mv.add(&a_v2);

            // OCCT L8303-8305: aGFE.Perform().
            let mut a_gfe_filler = PaveFiller::new();
            a_gfe_filler.set_arguments(a_gfe_args);
            let a_prog = NoopProgress;
            let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_Builder", 1);
            a_gfe_filler.perform(&a_ps);
            let mut a_gfe = Builder::new(
                a_gfe_filler.ds(),
                crate::bop::algo::builder::BooleanOpType::Union,
                a_gfe_filler.fuzzy_value(),
            );
            a_gfe.my_arguments = a_gfe_filler.ds().arguments.clone();
            let _a_root = a_gfe
                .build()
                .ok()
                .map(|brep| super::brep_offset_make_offset_1::brep_root_shape(&brep));
            let b_ok = _a_root.is_some() && !a_gfe.has_errors();
            if b_ok {
                // OCCT L8306-8320: get the images of the bounding vertices.
                let a_nb_v = a_mv.extent();
                for i_v in 1..=a_nb_v {
                    let a_v = a_mv.find_key_1(i_v).clone();
                    if the_mv_bounds.contains_key(&shape_key_of(&a_v))
                        || a_v.is_same(&a_v1)
                        || a_v.is_same(&a_v2)
                    {
                        let a_lv_im = builder_modified(&a_gfe, &a_v);
                        if a_lv_im.is_empty() {
                            set_add(&mut a_mev_bounds, &a_v);
                        } else {
                            set_add(&mut a_mev_bounds, &a_lv_im[0]);
                        }
                    }
                }
                // OCCT L8322: aCEIm = aGFE.Shape().
                a_ce_im = a_gfe
                    .my_shape
                    .as_ref()
                    .map(super::brep_offset_make_offset_1::brep_root_shape)
                    .unwrap_or_else(Shape::null);
            }
        } else {
            // OCCT L8324-8326: aCEIm = aLEIm.First().
            a_ce_im = a_le_im[0].clone();
        }

        // OCCT L8328: aLEIm.Clear().
        let mut a_le_im: Vec<Shape> = Vec::new();

        // OCCT L8330-8348: explore the split edges.
        for a_e_im in explorer(&a_ce_im, ShapeType::Edge, ShapeType::Shape) {
            // OCCT L8334-8343: check the split not to contain the bounding
            // vertices.
            let mut b_broke = false;
            for a_v in bat::sub_shapes(&a_e_im) {
                if set_contains(&a_mev_bounds, &a_v) || set_contains(the_mv_bounds, &a_v) {
                    b_broke = true;
                    break;
                }
            }
            if !b_broke {
                a_le_im.push(a_e_im.clone());
                // OCCT L8346: theDMEOr.Bound(aEIm, ..)->Append(aE).
                add_to_container_map_list(the_dme_or, &a_e_im, a_e);
            }
        }

        // OCCT L8350-8383.
        if a_le_im.is_empty() {
            shape_data_map::un_bind(the_e_images, a_e);
        } else {
            let a_lfe = shape_data_map::find(the_melf, a_e);
            for a_e_im in a_le_im.iter() {
                let entry = the_melf
                    .entry(shape_key_of(a_e_im))
                    .or_insert((a_e_im.clone(), Vec::new()));
                for a_f in a_lfe.iter() {
                    append_to_list(&mut entry.1, a_f);
                }
                if b_check_ext {
                    set_add(the_mecheckext, a_e_im);
                } else if !b_old {
                    set_add(the_me_new, a_e_im);
                }
            }
            // OCCT L8269-8275 write-back: the images list of aE.
            if let Some(entry) = the_e_images.get_mut(&shape_key_of(a_e)) {
                entry.1 = a_le_im;
            }
        }
    }
    let _ = a_bf;
}

// ---------------------------------------------------------------------------
// OCCT IntersectEdges (cxx L8408-8550).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::IntersectEdges (cxx L8408-8550) —
/// intersecting the trimmed edges to avoid the self-intersections.
#[allow(clippy::too_many_arguments)]
#[allow(unused_variables)]
pub(crate) fn intersect_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_la: &[Shape],
    the_le: &[Shape],
    the_mv_bounds: &OcctShapeSet,
    the_verts_to_avoid: &OcctShapeSet,
    the_me_new: &mut OcctShapeSet,
    the_mecheckext: &mut OcctShapeSet,
    the_e_images: &mut ShapeDataMap<Vec<Shape>>,
    the_dme_or: &mut ShapeDataMap<Vec<Shape>>,
    the_melf: &mut ShapeDataMap<Vec<Shape>>,
    the_splits: &mut Shape,
) {
    // OCCT L8414-8431: BOPAlgo_Builder aGFA; SetArguments(theLA);
    // Perform(); on errors the input is copied into the result compound.
    let a_gfa_args = the_la.to_vec();
    let mut a_gfa_filler = PaveFiller::new();
    a_gfa_filler.set_arguments(a_gfa_args);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_Builder", 1);
    a_gfa_filler.perform(&a_ps);
    let mut a_gfa = Builder::new(
        a_gfa_filler.ds(),
        crate::bop::algo::builder::BooleanOpType::Union,
        a_gfa_filler.fuzzy_value(),
    );
    a_gfa.my_arguments = a_gfa_filler.ds().arguments.clone();
    let a_root = a_gfa
        .build()
        .ok()
        .map(|brep| super::brep_offset_make_offset_1::brep_root_shape(&brep));
    let b_ok = a_root.is_some() && !a_gfa.has_errors();
    if !b_ok {
        // OCCT L8419-8429: just copy the input into the result.
        let mut a_sp = empty_compound();
        for a_e in the_la.iter() {
            add_to_container_shape(a_e, &mut a_sp);
        }
        *the_splits = a_sp;
        return;
    }

    // OCCT L8433: UpdateImages(theLE, theEImages, aGFA, myModifiedEdges).
    {
        let mut my_modified_edges = a_bf.my_modified_edges.clone();
        update_images(the_le, the_e_images, &a_gfa, &mut my_modified_edges);
        a_bf.my_modified_edges = my_modified_edges;
    }

    // OCCT L8435: theSplits = aGFA.Shape().
    *the_splits = a_root.unwrap_or_else(Shape::null);

    // OCCT L8439-8447: prepare the list of the edges to update.
    let mut a_le_input: Vec<Shape> = Vec::new();
    for a_s in the_la.iter() {
        for a_e in explorer(a_s, ShapeType::Edge, ShapeType::Shape) {
            a_le_input.push(a_e);
        }
    }

    // OCCT L8449-8465: update the new edges.
    for a_e in a_le_input.iter() {
        if !set_contains(the_me_new, a_e) {
            continue;
        }
        let a_le_im = builder_modified(&a_gfa, a_e);
        if a_le_im.is_empty() {
            continue;
        }
        the_me_new.remove(&shape_key_of(a_e));
        for a_e_im in a_le_im.iter() {
            set_add(the_me_new, a_e_im);
        }
    }

    // OCCT L8467-8492: update the edges after the intersection for the
    // extended checking.
    for a_e in a_le_input.iter() {
        let a_le_im = builder_modified(&a_gfa, a_e);
        if a_le_im.is_empty() {
            continue;
        }
        if set_contains(the_mecheckext, a_e) {
            for a_e_im in a_le_im.iter() {
                set_add(the_mecheckext, a_e_im);
            }
            the_mecheckext.remove(&shape_key_of(a_e));
        }
        let a_lfe = shape_data_map::find(the_melf, a_e);
        for a_e_im in a_le_im.iter() {
            let entry = the_melf
                .entry(shape_key_of(a_e_im))
                .or_insert((a_e_im.clone(), Vec::new()));
            for a_f in a_lfe.iter() {
                append_to_list(&mut entry.1, a_f);
            }
        }
    }

    // OCCT L8494-8497: GetInvalidEdges.
    let mut a_me_inv: OcctShapeSet = HashMap::new();
    a_bf.get_invalid_edges(
        the_verts_to_avoid,
        the_mv_bounds,
        &super::brep_offset_make_offset_1::BuilderRef::Builder(&a_gfa),
        &mut a_me_inv,
    );

    // OCCT L8498-8512: update the shape.
    if !a_me_inv.is_empty() {
        let mut a_sp = empty_compound();
        for a_e in explorer(the_splits, ShapeType::Edge, ShapeType::Shape) {
            if !set_contains(&a_me_inv, &a_e) {
                add_to_container_shape(&a_e, &mut a_sp);
            }
        }
        *the_splits = a_sp;
    }

    // OCCT L8514: update the origins.
    {
        let mut my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();
        update_origins(&a_le_input, &mut my_edges_origins, &a_gfa);
        a_bf.my_edges_origins = Some(my_edges_origins);
    }
}
