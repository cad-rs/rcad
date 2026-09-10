// OCCT BRepOffset_MakeOffset_1.cxx L6975-7606 — module h of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L6975-7606:
//   - ProcessCommonEdges (cxx L6975-7112)
//   - FindOrigins (cxx L7113-7158, file static)
//   - UpdateIntersectedFaces (cxx L7159-7232)
//   - IntersectFaces, the pair form (cxx L7233-7318)
//   - IntersectAndTrimEdges (cxx L7319-7606)
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#49).

use std::collections::HashMap;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo::topods::{ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::builder::Builder;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    append_to_list, builder_modified, empty_compound,
    find_common_parts, map_shapes_and_ancestors_map, map_shapes_indexed, shape_key_of,
    BRepOffsetBuildOffsetFaces,
};
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, OcctIndexedShapeMap, OcctShapeSet,
    ShapeDataMap, ShapeIndexedDataMap,
};

// ---------------------------------------------------------------------------
// OCCT ProcessCommonEdges (cxx L6975-7112).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::ProcessCommonEdges (cxx L6975-7112) —
/// analyzing the common edges between the splits of the offset faces.
#[allow(clippy::too_many_arguments)]
pub(crate) fn process_common_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_lec: &[Shape],
    the_me: &OcctIndexedShapeMap,
    the_meinf_etrim: &ShapeDataMap<Vec<Shape>>,
    the_all_invs: &OcctShapeSet,
    the_force_use: bool,
    the_mecv: &mut OcctIndexedShapeMap,
    the_mecheckext: &mut OcctShapeSet,
    the_dmeetrim: &mut ShapeDataMap<Vec<Shape>>,
    the_lfei: &mut Vec<Shape>,
    the_lfej: &mut Vec<Shape>,
    the_me_to_int: &mut OcctIndexedShapeMap,
) {
    let my_oe_origins = a_bf.my_oe_origins.clone();
    let my_e_trim_e_inf = a_bf.my_e_trim_e_inf.clone().unwrap_or_default();

    // OCCT L7000: aLEC.
    let mut a_lec: Vec<Shape> = Vec::new();

    // OCCT L7001-7030: process the common edges.
    for a_ec in the_lec.iter() {
        // OCCT L7004-7007: check first if the common edges are valid.
        if a_bf.my_invalid_edges.contains(a_ec) && !a_bf.my_valid_edges.contains(a_ec) {
            continue;
        }
        // OCCT L7009-7013: the common edge should have a connection to the
        // current invalidity.
        if the_me.contains(a_ec) {
            a_lec.push(a_ec.clone());
            continue;
        }
        // OCCT L7015-7029.
        let mut b_broke = false;
        for a_ve in bat::sub_shapes(a_ec) {
            if the_me.contains(&a_ve) {
                a_lec.push(a_ec.clone());
                b_broke = true;
                break;
            }
        }
        let _ = b_broke;
    }

    // OCCT L7031-7055.
    let b_use_only_inf = a_lec.is_empty();
    if b_use_only_inf {
        if the_force_use {
            a_lec = the_lec.to_vec();
        } else {
            for a_ec in the_lec.iter() {
                // OCCT L7042-7051: check if all the images of the origin of
                // this edge are not connected to any invalidity.
                let a_e_int = shape_data_map::find(&my_e_trim_e_inf, a_ec);
                let a_lve = shape_data_map::find(the_meinf_etrim, &a_e_int);
                let mut b_return = false;
                for a_ecx in a_lve.iter() {
                    if set_contains(the_all_invs, a_ecx) || a_bf.my_invalid_edges.contains(a_ecx)
                    {
                        b_return = true;
                        break;
                    }
                    for a_v in bat::sub_shapes(a_ecx) {
                        if set_contains(the_all_invs, &a_v) {
                            b_return = true;
                            break;
                        }
                    }
                    if b_return {
                        break;
                    }
                }
                if b_return {
                    return;
                }
                // OCCT L7053-7054: use only one element.
                if a_lec.is_empty() {
                    a_lec.push(a_ec.clone());
                }
            }
        }
    }

    // OCCT L7057-7090.
    for a_ec in a_lec.iter() {
        let a_e_int = shape_data_map::find(&my_e_trim_e_inf, a_ec);
        if !b_use_only_inf {
            // OCCT L7062-7072: find the edges of the same original edge and
            // take their vertices as well.
            let a_lve = shape_data_map::find(the_meinf_etrim, &a_e_int);
            for a_ecx in a_lve.iter() {
                let alone_origin = match shape_data_map::seek(&my_oe_origins, a_ecx)
                {
                    Some(p) => p.len() == 1,
                    None => true,
                };
                if !shape_data_map::is_bound(&my_oe_origins, a_ecx)
                    || alone_origin
                {
                    map_shapes_indexed(a_ecx, ShapeType::Vertex, the_mecv);
                }
            }
            // OCCT L7074-7081: bind the unlimited edge to its trimmed part.
            {
                let entry = the_dmeetrim
                    .entry(shape_key_of(&a_e_int))
                    .or_insert((a_e_int.clone(), Vec::new()));
                append_to_list(&mut entry.1, a_ec);
            }
        } else if !the_force_use {
            // OCCT L7083-7085.
            set_add(the_mecheckext, &a_e_int);
        }
        // OCCT L7087-7089.
        append_to_list(the_lfei, &a_e_int);
        append_to_list(the_lfej, &a_e_int);
        the_me_to_int.add(&a_e_int);
    }
}

// ---------------------------------------------------------------------------
// OCCT FindOrigins (cxx L7113-7158, file static).
// ---------------------------------------------------------------------------

/// OCCT static FindOrigins (cxx L7113-7158) — finds the origin edges of
/// the given images.
pub(crate) fn find_origins(
    the_lf_im1: &[Shape],
    the_lf_im2: &[Shape],
    the_me: &OcctIndexedShapeMap,
    the_origins: &ShapeDataMap<Vec<Shape>>,
    the_le_or: &mut Vec<Shape>,
) {
    // OCCT L7117.
    let mut a_m_fence: OcctShapeSet = HashMap::new();
    for i in 0..2 {
        let a_lf = if i == 0 { the_lf_im1 } else { the_lf_im2 };
        for a_f in a_lf.iter() {
            for a_e in explorer(a_f, ShapeType::Edge, ShapeType::Shape) {
                // OCCT L7128-7138.
                if the_me.contains(&a_e) && shape_data_map::is_bound(the_origins, &a_e) {
                    let a_le_or = shape_data_map::find(the_origins, &a_e);
                    for a_e_or in a_le_or.iter() {
                        if set_add(&mut a_m_fence, a_e_or)
                            && a_e_or.shape_type() == ShapeType::Edge
                        {
                            the_le_or.push(a_e_or.clone());
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT UpdateIntersectedFaces (cxx L7159-7232).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::UpdateIntersectedFaces (cxx
/// L7159-7232) — updating the already interfered faces.
#[allow(clippy::too_many_arguments)]
pub(crate) fn update_intersected_faces_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_f_inv: &Shape,
    the_fi: &Shape,
    the_fj: &Shape,
    the_lf_inv: &[Shape],
    the_lf_imi: &[Shape],
    the_lf_imj: &[Shape],
    the_lfei: &[Shape],
    the_lfej: &[Shape],
    the_me_to_int: &mut OcctIndexedShapeMap,
) {
    let mut my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();

    // OCCT L7164-7171: find the common edges in these two lists.
    let mut a_mei: OcctShapeSet = HashMap::new();
    for a_e in the_lfei.iter() {
        set_add(&mut a_mei, a_e);
    }

    // OCCT L7173-7183: find the origins.
    let mut a_me_to_find_origins = OcctIndexedShapeMap::new();
    let mut a_le_to_find_origins: Vec<Shape> = Vec::new();
    if !the_fi.is_same(the_f_inv) {
        find_common_parts(the_lf_imi, the_lf_inv, &mut a_le_to_find_origins, ShapeType::Edge);
    }
    if !the_fj.is_same(the_f_inv) {
        find_common_parts(the_lf_imj, the_lf_inv, &mut a_le_to_find_origins, ShapeType::Edge);
    }
    for a_ec in a_le_to_find_origins.iter() {
        a_me_to_find_origins.add(a_ec);
    }

    // OCCT L7185.
    let mut a_le_or_init: Vec<Shape> = Vec::new();
    find_origins(
        the_lf_imi,
        the_lf_imj,
        &a_me_to_find_origins,
        &my_edges_origins,
        &mut a_le_or_init,
    );

    // OCCT L7187-7214.
    for a_e in the_lfej.iter() {
        if set_contains(&a_mei, a_e) {
            the_me_to_int.add(a_e);
            if !a_le_or_init.is_empty() {
                if shape_data_map::is_bound(&my_edges_origins, a_e) {
                    let a_le_or = shape_data_map::change_find(&mut my_edges_origins, a_e);
                    for a_e_or in a_le_or_init.iter() {
                        append_to_list(a_le_or, a_e_or);
                    }
                } else {
                    shape_data_map::bind(&mut my_edges_origins, a_e, a_le_or_init.clone());
                }
            }
        }
    }

    a_bf.my_edges_origins = Some(my_edges_origins);
}

// ---------------------------------------------------------------------------
// OCCT IntersectFaces, the pair form (cxx L7233-7318).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::IntersectFaces (cxx L7233-7318) —
/// intersection of the pair of faces.
#[allow(clippy::too_many_arguments)]
pub(crate) fn intersect_faces_pair_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_f_inv: &Shape,
    the_fi: &Shape,
    the_fj: &Shape,
    the_lf_inv: &[Shape],
    the_lf_imi: &[Shape],
    the_lf_imj: &[Shape],
    the_lfei: &mut Vec<Shape>,
    the_lfej: &mut Vec<Shape>,
    the_mecv: &mut OcctIndexedShapeMap,
    the_me_to_int: &mut OcctIndexedShapeMap,
) {
    let mut my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();

    // OCCT L7239-7250: intersect the faces — BRepOffset_Tool::Inter3D
    // (the brep_offset_tool_b.rs translation).
    let a_side = State::Out;
    let mut a_l_int1: Vec<Shape> = Vec::new();
    let mut a_l_int2: Vec<Shape> = Vec::new();
    let a_null_edge = Shape::null();
    let a_null_face = Shape::null();
    super::brep_offset_tool_b::inter3d(
        the_fi,
        the_fj,
        &mut a_l_int1,
        &mut a_l_int2,
        a_side,
        &a_null_edge,
        &a_null_face,
        &a_null_face,
    );

    // OCCT L7252-7255.
    if a_l_int1.is_empty() {
        return;
    }

    // OCCT L7257-7265: find the common vertices for trimming the edges.
    let mut a_lcv: Vec<Shape> = Vec::new();
    find_common_parts(the_lf_imi, the_lf_imj, &mut a_lcv, ShapeType::Vertex);
    if a_lcv.len() > 1 {
        for a_cv in a_lcv.iter() {
            the_mecv.add(a_cv);
        }
    }

    // OCCT L7267-7284: find the origins.
    let mut a_me_to_find_origins = OcctIndexedShapeMap::new();
    let mut a_le_to_find_origins: Vec<Shape> = Vec::new();
    if !the_fi.is_same(the_f_inv) {
        find_common_parts(the_lf_imi, the_lf_inv, &mut a_le_to_find_origins, ShapeType::Edge);
    }
    if !the_fj.is_same(the_f_inv) {
        find_common_parts(the_lf_imj, the_lf_inv, &mut a_le_to_find_origins, ShapeType::Edge);
    }
    for a_ec in a_le_to_find_origins.iter() {
        a_me_to_find_origins.add(a_ec);
    }
    let mut a_le_or_init: Vec<Shape> = Vec::new();
    find_origins(
        the_lf_imi,
        the_lf_imj,
        &a_me_to_find_origins,
        &my_edges_origins,
        &mut a_le_or_init,
    );

    // OCCT L7286-7303.
    for a_e_int in a_l_int1.iter() {
        the_lfei.push(a_e_int.clone());
        the_lfej.push(a_e_int.clone());
        if !a_le_or_init.is_empty() {
            shape_data_map::bind(&mut my_edges_origins, a_e_int, a_le_or_init.clone());
        }
        the_me_to_int.add(a_e_int);
    }

    a_bf.my_edges_origins = Some(my_edges_origins);
}

// ---------------------------------------------------------------------------
// OCCT IntersectAndTrimEdges (cxx L7319-7606).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::IntersectAndTrimEdges (cxx
/// L7319-7606) — intersection of the new intersection edges among
/// themselves.
#[allow(clippy::too_many_arguments, unused_mut)]
pub(crate) fn intersect_and_trim_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_mf_int: &OcctIndexedShapeMap,
    the_me_int: &OcctIndexedShapeMap,
    the_dmeetrim: &ShapeDataMap<Vec<Shape>>,
    the_ms_inv: &OcctIndexedShapeMap,
    the_mve: &OcctIndexedShapeMap,
    the_verts_to_avoid: &OcctShapeSet,
    the_new_verts_to_avoid: &OcctShapeSet,
    the_mecheckext: &OcctShapeSet,
    the_ss_interfs: Option<&ShapeDataMap<Vec<Shape>>>,
    the_mv_bounds: &mut OcctShapeSet,
    the_e_images: &mut ShapeDataMap<Vec<Shape>>,
) {
    // OCCT L7333-7336.
    let a_nb = the_me_int.extent();
    if a_nb == 0 {
        return;
    }

    // OCCT L7338-7342.
    let mut a_largs: Vec<Shape> = Vec::new();
    let mut a_m_fence: OcctShapeSet = HashMap::new();

    // OCCT L7344-7368: get the vertices from the splits of the intersected
    // faces (the edges close to invalidity).
    let mut a_dmve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let a_nb = the_mf_int.extent();
    for i in 1..=a_nb {
        let a_f = the_mf_int.find_key_1(i).clone();
        let a_le = match a_bf.my_faces_to_rebuild.get(&shape_key_of(&a_f)) {
            Some(v) => v.1.clone(),
            None => continue,
        };
        for a_e in a_le.iter() {
            map_shapes_and_ancestors_map(
                a_e,
                ShapeType::Vertex,
                ShapeType::Edge,
                &mut a_dmve,
            );
            for a_v1 in explorer(a_e, ShapeType::Vertex, ShapeType::Shape) {
                if !set_contains(the_verts_to_avoid, &a_v1)
                    && the_mve.contains(&a_v1)
                    && set_add(&mut a_m_fence, &a_v1)
                {
                    a_largs.push(a_v1.clone());
                }
            }
        }
    }

    // OCCT L7370-7405.
    let a_nb = the_ms_inv.extent();
    for i in 1..=a_nb {
        let a_s = the_ms_inv.find_key_1(i).clone();
        // OCCT L7374-7386: the edge case.
        if let Some(the_ss_interfs) = the_ss_interfs {
            if let Some(p_lv) = shape_data_map::seek(the_ss_interfs, &a_s) {
                let p_lv = p_lv.clone();
                for a_v in p_lv.iter() {
                    if a_v.shape_type() == ShapeType::Vertex {
                        a_largs.push(a_v.clone());
                    }
                }
            }
        }
        // OCCT L7388-7405: the vertex case.
        if let Some(p_lve) = a_dmve.get(&shape_key_of(&a_s)) {
            let p_lve = p_lve.1.clone();
            for a_e in p_lve.iter() {
                for a_v1 in explorer(a_e, ShapeType::Vertex, ShapeType::Shape) {
                    if !set_contains(the_verts_to_avoid, &a_v1)
                        && set_add(&mut a_m_fence, &a_v1)
                    {
                        a_largs.push(a_v1.clone());
                    }
                }
            }
        }
    }

    // OCCT L7407-7432: the bounding vertices of the untrimmed edges, the
    // new intersection edges, the edges to intersect, the common edges.
    let mut a_lv_bounds: Vec<Shape> = Vec::new();
    let mut a_le_new: Vec<Shape> = Vec::new();
    let mut a_le_int: Vec<Shape> = Vec::new();
    let mut a_lce: Vec<Shape> = Vec::new();
    let a_nb = the_me_int.extent();
    for i in 1..=a_nb {
        let a_e = the_me_int.find_key_1(i).clone();
        if set_contains(the_mecheckext, &a_e) {
            // OCCT L7416-7420: avoid the trimming of the intersection edges
            // by the additional common edges.
            a_lce.push(a_e);
            continue;
        }
        if !shape_data_map::is_bound(the_dmeetrim, &a_e) {
            a_le_new.push(a_e.clone());
        }
        a_le_int.push(a_e.clone());
        a_largs.push(a_e.clone());
        for a_v in explorer(&a_e, ShapeType::Vertex, ShapeType::Shape) {
            a_lv_bounds.push(a_v.clone());
        }
    }
    let _ = &a_le_new;

    // OCCT L7434-7440: BOPAlgo_Builder aGF; SetArguments(aLArgs);
    // Perform(); if (aGF.HasErrors()) return.
    let mut a_gf_filler = PaveFiller::new();
    a_gf_filler.set_arguments(a_largs);
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
    if _a_root.is_none() || a_gf.has_errors() {
        return;
    }

    // OCCT L7442-7453: update the vertices to avoid with the SD vertices.
    for a_v in a_lv_bounds.iter() {
        let a_lv_im = builder_modified(&a_gf, a_v);
        if a_lv_im.is_empty() {
            set_add(the_mv_bounds, a_v);
        } else {
            let a_v_im = &a_lv_im[0];
            set_add(the_mv_bounds, a_v_im);
        }
    }

    // OCCT L7455-7456: find the invalid splits of the edges.
    let mut a_me_inv: OcctShapeSet = HashMap::new();
    a_bf.get_invalid_edges(
        the_new_verts_to_avoid,
        the_mv_bounds,
        &super::brep_offset_make_offset_1::BuilderRef::Builder(&a_gf),
        &mut a_me_inv,
    );

    // OCCT L7458-7462: get the valid splits to intersect with the commons.
    let mut a_ce_im = empty_compound();

    // OCCT L7464-7488: remove the splits containing the vertices from the
    // invalid edges.
    for a_e in a_le_int.iter() {
        let mut a_le_im = builder_modified(&a_gf, a_e);
        if a_le_im.is_empty() {
            continue;
        }
        // OCCT L7474-7484: the iterator-removal walk.
        let mut j = 0usize;
        while j < a_le_im.len() {
            let a_e_im = a_le_im[j].clone();
            if set_contains(&a_me_inv, &a_e_im) {
                a_le_im.remove(j);
            } else {
                bat::builder_add_compound_shape(&mut a_ce_im, &a_e_im);
                j += 1;
            }
        }
        // OCCT L7486-7488.
        if !a_le_im.is_empty() {
            let entry = the_e_images
                .entry(shape_key_of(a_e))
                .or_insert((a_e.clone(), Vec::new()));
            entry.1.extend(a_le_im.iter().cloned());
        }
    }

    // OCCT L7490-7493.
    if a_lce.is_empty() {
        return;
    }

    // OCCT L7495-7505: trim the common edges by the other intersection
    // edges — BOPAlgo_Builder aGFCE; SetArguments(aLCE);
    // AddArgument(aCEIm); Perform(); HasErrors.
    let mut a_gfce_args = a_lce.clone();
    a_gfce_args.push(a_ce_im.clone());
    let mut a_gfce_filler = PaveFiller::new();
    a_gfce_filler.set_arguments(a_gfce_args);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_Builder", 1);
    a_gfce_filler.perform(&a_ps);
    let mut a_gfce = Builder::new(
        a_gfce_filler.ds(),
        crate::bop::algo::builder::BooleanOpType::Union,
        a_gfce_filler.fuzzy_value(),
    );
    a_gfce.my_arguments = a_gfce_filler.ds().arguments.clone();
    let _a_root = a_gfce
        .build()
        .ok()
        .map(|brep| super::brep_offset_make_offset_1::brep_root_shape(&brep));
    if _a_root.is_none() || a_gfce.has_errors() {
        return;
    }

    // OCCT L7507-7605: the common-edge walk over the DS (the
    // pDS->PaveBlocks / IsCommonBlock / CommonBlock->PaveBlocks /
    // OriginalEdge re-hosts).
    let p_ds = a_gfce.ds;
    for a_e in a_lce.iter() {
        let mut a_le_im = builder_modified(&a_gfce, a_e);
        if a_le_im.is_empty() {
            continue;
        }
        // OCCT L7518-7545: check if it does not coincide with some
        // intersection edge.
        let n_e = p_ds.index(a_e);
        let mut b_broke_lp = false;
        let mut b_broke_lpbc = false;
        if n_e >= 0 {
            let a_lpb = p_ds.pave_blocks(n_e as usize);
            for a_pb in a_lpb.iter() {
                if p_ds.is_common_block(a_pb) {
                    // OCCT L7527-7541: find with what it is a common.
                    if let Some(cb_idx) = p_ds.common_block(a_pb) {
                        let a_lpbc = p_ds.common_blocks[cb_idx].pave_blocks();
                        for (a_pb_c, _) in a_lpbc.iter() {
                            let n_e_orig: usize = a_pb_c.read().original_edge();
                            let a_ec = p_ds.shape(n_e_orig);
                            if !set_contains(the_mecheckext, a_ec) {
                                b_broke_lpbc = true;
                                break;
                            }
                        }
                    }
                    if b_broke_lpbc {
                        b_broke_lp = true;
                        break;
                    }
                }
            }
        }
        if b_broke_lp {
            // OCCT L7547-7549: avoid the creation of the unnecessary splits
            // from the commons which coincide with the intersection edges.
            continue;
        }

        // OCCT L7551-7557: save the images.
        {
            let entry = the_e_images
                .entry(shape_key_of(a_e))
                .or_insert((a_e.clone(), Vec::new()));
            entry.1.extend(a_le_im.iter().cloned());
        }
        // OCCT L7559-7565: save the bounding vertices.
        for a_v in bat::sub_shapes(a_e) {
            let a_lv_im = builder_modified(&a_gfce, &a_v);
            if a_lv_im.is_empty() {
                set_add(the_mv_bounds, &a_v);
            } else {
                set_add(the_mv_bounds, &a_lv_im[0]);
            }
        }
    }
}

