// OCCT BRepOffset_MakeOffset_1.cxx L3893-5167 — module e of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L3893-5167:
//   - RemoveInvalidSplitsFromValid (cxx L3893-4070)
//   - buildPairs (cxx L4076-4106, file static)
//   - buildIntersectionPairs (cxx L4110-4287, file static)
//   - RemoveInsideFaces (cxx L4295-4636)
//   - ShapesConnections (cxx L4643-4967)
//   - RemoveHangingParts (cxx L4968-5167)
//
// The BuilderRef carrier (architecture difference #49 of
// brep_offset_make_offset_1.rs) carries the OCCT `BOPAlgo_Builder&` /
// `BOPAlgo_MakerVolume` base-class references; the post-perform query
// surfaces that the rcad MakerVolume re-host does not expose (Images /
// Origins / Modified / IsDeleted / HasDeleted / PDS) take the annotated
// GAP paths with the OCCT empty-result structure preserved.
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#49).

use std::collections::HashMap;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::maker_volume::MakerVolume;
use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    add_to_container_map_list, block_shape, empty_compound, find_shape,
    map_shapes_and_ancestors_map, map_shapes_indexed, shape_key_of, take_modified_map_no_fence,
    take_modified_no_fence, take_modified_no_fence_indexed, take_modified_shape,
    BRepOffsetBuildOffsetFaces, BuilderRef,
};
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, shape_indexed_data_map, OcctIndexedShapeMap,
    OcctShapeSet, ShapeDataMap, ShapeIndexedDataMap,
};

// ---------------------------------------------------------------------------
// OCCT RemoveInvalidSplitsFromValid (cxx L3893-4070).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::RemoveInvalidSplitsFromValid (cxx
/// L3893-4070) — removing the invalid splits of faces from the valid ones:
/// 1. make the connexity blocks from the invalid faces;
/// 2. find the free edges in these blocks;
/// 3. if all free edges are valid for the faces — remove the block.
pub(crate) fn remove_invalid_splits_from_valid_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_dmfmvie: &ShapeDataMap<OcctShapeSet>,
) {
    // OCCT L3904-3908.
    let mut a_m_fence: OcctShapeSet = HashMap::new();
    let mut a_mf_to_rem: OcctShapeSet = HashMap::new();
    let mut a_cf_inv = empty_compound();

    // OCCT L3910-3932: make the compound of the invalid faces.
    let mut a_dmifof: ShapeDataMap<Shape> = HashMap::new();
    let a_nb = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb {
        let (a_f, a_lf_inv) = {
            let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        // OCCT L3917-3920: artificially invalid faces should not be
        // removed.
        if shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f) {
            continue;
        }
        for a_f_im in a_lf_inv.iter() {
            if set_add(&mut a_m_fence, a_f_im) {
                bat::builder_add_compound_shape(&mut a_cf_inv, a_f_im);
                shape_data_map::bind(&mut a_dmifof, a_f_im, a_f.clone());
            }
        }
    }

    // OCCT L3935-3936: make the connexity blocks (the shell_splitter
    // re-host of the BOPTools_AlgoTools::MakeConnexityBlocks EDGE/FACE
    // form, architecture difference #44).
    let a_cf_faces = explorer(&a_cf_inv, ShapeType::Face, ShapeType::Shape);
    let a_lcb_inv = crate::bop::algo::shell_splitter::make_connexity_blocks(&a_cf_faces);

    // OCCT L3938-4045: analyze each block.
    for a_block in a_lcb_inv.iter() {
        // OCCT L3942: aCB — the connexity-block list carrier.
        let a_cb = block_shape(&a_block.shapes);

        // OCCT L3944-3975: if the connexity block contains only one face —
        // it should be removed.  (The OCCT explorer double-pass form: the
        // first pass checks the second face, the second pass checks the
        // valid images left; b_has_more_faces/b_broke mirror the
        // aExp.More() state.)
        let a_cb_faces: Vec<Shape> = explorer(&a_cb, ShapeType::Face, ShapeType::Shape);
        let b_has_more_faces = a_cb_faces.len() > 1;
        let mut b_broke = false;
        if b_has_more_faces {
            // OCCT L3950-3963: check if there are valid images left.
            for a_f_im in a_cb_faces.iter() {
                let a_f = shape_data_map::find(&a_dmifof, a_f_im);
                let a_lf_im = match a_bf.my_of_images.get(&shape_key_of(&a_f)) {
                    Some(v) => v.1.clone(),
                    None => continue,
                };
                let a_lf_inv = match a_bf.my_invalid_faces.get(&shape_key_of(&a_f)) {
                    Some(v) => v.1.clone(),
                    None => continue,
                };
                if a_lf_im.len() == a_lf_inv.len() {
                    b_broke = true;
                    break;
                }
            }
        }
        // OCCT L3966-3975.
        if !b_broke {
            for a_f in a_cb_faces.iter() {
                set_add(&mut a_mf_to_rem, a_f);
            }
            continue;
        }

        // OCCT L3977-3982: remove the faces connected by the inverted
        // edges.
        let mut a_dmef: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
        map_shapes_and_ancestors_map(&a_cb, ShapeType::Edge, ShapeType::Face, &mut a_dmef);

        // OCCT L3984-3997.
        let mut a_dmff: ShapeDataMap<Vec<Shape>> = HashMap::new();
        for a_fcb in a_cb_faces.iter() {
            let a_f = shape_data_map::find(&a_dmifof, a_fcb);
            add_to_container_map_list(&mut a_dmff, &a_f, a_fcb);
        }

        // OCCT L3999-4044.
        for (a_f, a_lfcb) in a_dmff.values() {
            let p_valid_inverted = shape_data_map::seek(the_dmfmvie, a_f);

            // OCCT L4009-4036: either remove all of these faces or none.
            let mut b_found_kept = false;
            for a_fcb in a_lfcb.iter() {
                let mut b_broke_e = false;
                for a_ecb in explorer(a_fcb, ShapeType::Edge, ShapeType::Shape) {
                    if p_valid_inverted
                        .map(|m| set_contains(m, &a_ecb))
                        .unwrap_or(false)
                    {
                        b_broke_e = true;
                        break;
                    }
                    if shape_indexed_data_map::find(&a_dmef, &a_ecb).len() > 1
                        && !a_bf.my_inverted_edges.contains(&a_ecb)
                    {
                        b_broke_e = true;
                        break;
                    }
                }
                // OCCT L4031-4035: if one removed — remove all.
                if !b_broke_e {
                    b_found_kept = true;
                    break;
                }
            }
            // OCCT L4037-4043: if (itL.More()).
            if b_found_kept {
                for a_fcb in a_lfcb.iter() {
                    set_add(&mut a_mf_to_rem, a_fcb);
                }
            }
        }
    }

    // OCCT L4047-4069: remove the invalid faces from the images.
    if !a_mf_to_rem.is_empty() {
        let a_nb = a_bf.my_invalid_faces.len();
        for i in 1..=a_nb {
            let (a_f, _a_lf_images) = {
                let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
                (v.0.clone(), v.1.clone())
            };
            let a_key = shape_key_of(&a_f);
            let mut a_lf_images = match a_bf.my_of_images.get(&a_key) {
                Some(v) => v.1.clone(),
                None => continue,
            };
            // OCCT L4055-4067: the iterator-removal walk.
            let mut j = 0usize;
            while j < a_lf_images.len() {
                let a_f_im = a_lf_images[j].clone();
                if set_contains(&a_mf_to_rem, &a_f_im) {
                    a_lf_images.remove(j);
                } else {
                    j += 1;
                }
            }
            if let Some(entry) = a_bf.my_of_images.get_mut(&a_key) {
                entry.1 = a_lf_images;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT buildPairs (cxx L4076-4106, file static).
// ---------------------------------------------------------------------------

/// OCCT static buildPairs (cxx L4076-4106) — builds all the pairs of the
/// shapes of the map.
pub(crate) fn build_pairs(
    the_smap: &OcctIndexedShapeMap,
    the_int_pairs: &mut ShapeDataMap<OcctShapeSet>,
) {
    let a_nb_s = the_smap.extent();
    if a_nb_s < 2 {
        return;
    }
    for it1 in 1..=a_nb_s {
        let a_s = the_smap.find_key_1(it1);
        if !shape_data_map::is_bound(the_int_pairs, a_s) {
            shape_data_map::bind(the_int_pairs, a_s, HashMap::new());
        }
    }

    for it1 in 1..=a_nb_s {
        let a_s1 = the_smap.find_key_1(it1).clone();
        for it2 in (it1 + 1)..=a_nb_s {
            let a_s2 = the_smap.find_key_1(it2);
            // OCCT L4102: aMap1.Add(aS2).
            {
                let a_map1 = shape_data_map::change_find(the_int_pairs, &a_s1);
                set_add(a_map1, a_s2);
            }
            // OCCT L4103: theIntPairs(aS2).Add(aS1).
            {
                let a_map2 = shape_data_map::change_find(the_int_pairs, a_s2);
                set_add(a_map2, &a_s1);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT buildIntersectionPairs (cxx L4110-4287, file static).
// ---------------------------------------------------------------------------

/// OCCT static buildIntersectionPairs (cxx L4110-4287) — finds all the
/// faces connected to the not-removed faces and builds the intersection
/// pairs among them; for the removed faces intersects only those connected
/// to each other.
pub(crate) fn build_intersection_pairs(
    my_of_images: &ShapeIndexedDataMap<Vec<Shape>>,
    my_invalid_faces: &ShapeIndexedDataMap<Vec<Shape>>,
    the_builder: &BuilderRef,
    the_mf_removed: &OcctShapeSet,
    the_f_origins: &ShapeDataMap<Shape>,
    the_int_pairs: &mut ShapeDataMap<ShapeDataMap<OcctShapeSet>>,
) {
    let a_ctype = ShapeType::Vertex;
    // OCCT L4127-4130: the connection map from the vertices to the faces.
    let mut a_dmvf: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    if let Some(a_shape) = the_builder.shape() {
        map_shapes_and_ancestors_map(&a_shape, a_ctype, ShapeType::Face, &mut a_dmvf);
    }

    // OCCT L4132-4135: the images and origins of the builder.
    let an_images = the_builder.images();
    let an_origins = the_builder.origins();

    // OCCT L4140-4286.
    let a_nb_f = my_invalid_faces.len();
    for i_f in 1..=a_nb_f {
        let (a_f_inv, a_lf_invalid) = {
            let (_k, v) = my_invalid_faces.get_index(i_f - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };

        // OCCT L4144-4146.
        let mut a_cf = empty_compound();
        let mut a_cf_rem = empty_compound();

        // OCCT L4148-4172.
        for i_c in 0..2 {
            let a_lf: Vec<Shape> = if i_c == 0 {
                a_lf_invalid.clone()
            } else {
                my_of_images
                    .get(&shape_key_of(&a_f_inv))
                    .map(|v| v.1.clone())
                    .unwrap_or_default()
            };
            for a_s in a_lf.iter() {
                let mut a_lf_im: Vec<Shape> = Vec::new();
                // OCCT L4156: TakeModified(it.Value(), anImages, aLFIm) —
                // the no-fence form.
                take_modified_no_fence(a_s, &an_images, &mut a_lf_im);
                for a_f_im in a_lf_im.iter() {
                    if set_contains(the_mf_removed, a_f_im) {
                        bat::builder_add_compound_shape(&mut a_cf_rem, a_f_im);
                    } else {
                        bat::builder_add_compound_shape(&mut a_cf, a_f_im);
                    }
                }
            }
        }

        // OCCT L4174-4180: the connexity blocks of the not-removed faces
        // (the shell_splitter EDGE/FACE re-host form).
        let a_cf_face_list = explorer(&a_cf, ShapeType::Face, ShapeType::Shape);
        let a_lcb = crate::bop::algo::shell_splitter::make_connexity_blocks(&a_cf_face_list);
        if a_lcb.is_empty() {
            continue;
        }

        // OCCT L4182-4188: pFInterMap = theIntPairs.Bound(aFInv, empty).
        if !shape_data_map::is_bound(the_int_pairs, &a_f_inv) {
            shape_data_map::bind(the_int_pairs, &a_f_inv, HashMap::new());
        }

        // OCCT L4190-4224: build the pairs for the not-removed faces.
        for a_block in a_lcb.iter() {
            let a_cb = block_shape(&a_block.shapes);
            let mut a_mf_inter = OcctIndexedShapeMap::new();
            for a_cs in explorer(&a_cb, a_ctype, ShapeType::Shape) {
                let p_lfv = match a_dmvf.get(&shape_key_of(&a_cs)) {
                    Some(v) => v.1.clone(),
                    None => continue,
                };
                for a_f_connected in p_lfv.iter() {
                    let mut a_lf_or: Vec<Shape> = Vec::new();
                    // OCCT L4210: TakeModified(aFConnected, anOrigins,
                    // aLFOr) — the no-fence form.
                    take_modified_no_fence(a_f_connected, &an_origins, &mut a_lf_or);
                    for a_it_or in a_lf_or.iter() {
                        if let Some(p_for) = shape_data_map::seek(the_f_origins, a_it_or) {
                            a_mf_inter.add(p_for);
                        }
                    }
                }
            }
            // OCCT L4223: buildPairs(aMFInter, *pFInterMap).
            let p_f_inter_map = shape_data_map::change_find(the_int_pairs, &a_f_inv);
            build_pairs(&a_mf_inter, p_f_inter_map);
        }

        // OCCT L4226-4232: the connexity blocks of the removed faces.
        let a_cf_rem_face_list = explorer(&a_cf_rem, ShapeType::Face, ShapeType::Shape);
        let a_lcb = crate::bop::algo::shell_splitter::make_connexity_blocks(&a_cf_rem_face_list);
        if a_lcb.is_empty() {
            continue;
        }

        // OCCT L4234-4285.
        for a_block in a_lcb.iter() {
            let a_cb = block_shape(&a_block.shapes);

            let mut a_dmef: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
            for a_cs in explorer(&a_cb, a_ctype, ShapeType::Shape) {
                let p_lfv = match a_dmvf.get(&shape_key_of(&a_cs)) {
                    Some(v) => v.1.clone(),
                    None => continue,
                };
                for a_f_connected in p_lfv.iter() {
                    map_shapes_and_ancestors_map(
                        a_f_connected,
                        ShapeType::Edge,
                        ShapeType::Face,
                        &mut a_dmef,
                    );
                }
            }

            for i_e in 1..=a_dmef.len() {
                let a_lf_connected = shape_indexed_data_map::value_1(&a_dmef, i_e).clone();
                if a_lf_connected.len() < 2 {
                    continue;
                }
                let mut a_mf_inter = OcctIndexedShapeMap::new();
                for a_f_connected in a_lf_connected.iter() {
                    let mut a_lf_or: Vec<Shape> = Vec::new();
                    take_modified_no_fence(a_f_connected, &an_origins, &mut a_lf_or);
                    for a_it_or in a_lf_or.iter() {
                        if let Some(p_for) = shape_data_map::seek(the_f_origins, a_it_or) {
                            a_mf_inter.add(p_for);
                        }
                    }
                }
                // OCCT L4283.
                let p_f_inter_map = shape_data_map::change_find(the_int_pairs, &a_f_inv);
                build_pairs(&a_mf_inter, p_f_inter_map);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT RemoveInsideFaces (cxx L4295-4636).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::RemoveInsideFaces (cxx L4295-4636) —
/// looking for the inside faces that can be safely removed.
pub(crate) fn remove_inside_faces_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_inverted_faces: &[Shape],
    the_mf_to_check_int: &OcctIndexedShapeMap,
    the_mf_inv_in_hole: &OcctIndexedShapeMap,
    the_f_holes: &Shape,
    the_me_removed: &mut OcctIndexedShapeMap,
    _the_range: &ProgressScope,
) {
    // OCCT L4304-4308.
    let mut a_ls: Vec<Shape> = Vec::new();
    let mut a_m_fence: OcctShapeSet = HashMap::new();
    let mut a_mf_inv = OcctIndexedShapeMap::new();
    let mut a_dmf_im_f: ShapeDataMap<Shape> = HashMap::new();

    let _a_ps = ProgressScope::new(&NoopProgress, "Looking for inside faces", 10);

    // OCCT L4311-4350.
    let a_nb = a_bf.my_of_images.len();
    for i in 1..=a_nb {
        let (a_f, a_lf_images_i) = {
            let (_k, v) = a_bf.my_of_images.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        // OCCT L4319-4323: to avoid the intersection of the splits of the
        // same offset faces among themselves — the compound of the splits
        // is used as one argument.
        let mut a_cf_imi = empty_compound();
        for j in 0..2 {
            // OCCT L4327: pLFSp = !j ? myInvalidFaces.Seek(aF) :
            // &myOFImages(i).
            let p_lf_sp: Option<Vec<Shape>> = if j == 0 {
                super::brep_offset_make_offset_1::idm_seek(&a_bf.my_invalid_faces, &a_f)
            } else {
                Some(a_lf_images_i.clone())
            };
            let p_lf_sp = match p_lf_sp {
                Some(v) => v,
                None => continue,
            };
            for a_f_im in p_lf_sp.iter() {
                if set_add(&mut a_m_fence, a_f_im) {
                    bat::builder_add_compound_shape(&mut a_cf_imi, a_f_im);
                    shape_data_map::bind(&mut a_dmf_im_f, a_f_im, a_f.clone());
                    if j == 0 {
                        a_mf_inv.add(a_f_im);
                    }
                }
            }
        }
        // OCCT L4349: aLS.Append(aCFImi).
        a_ls.push(a_cf_imi);
    }

    // OCCT L4352-4362: add the faces consisting only of the invalid edges.
    let a_nb = the_mf_to_check_int.extent();
    for i in 1..=a_nb {
        let a_f_sp = the_mf_to_check_int.find_key_1(i).clone();
        if set_add(&mut a_m_fence, &a_f_sp) {
            a_ls.push(a_f_sp);
        }
    }

    // OCCT L4364-4371: BOPAlgo_MakerVolume aMV; SetArguments(aLS);
    // SetIntersect(true); Perform(); if (aMV.HasErrors()) return.
    let mut a_mv = MakerVolume::new();
    a_mv.set_arguments(a_ls);
    a_mv.set_intersect(true);
    a_mv.perform();
    if a_mv.has_errors() {
        return;
    }
    let a_mv_ref = BuilderRef::MakerVolume(&a_mv);

    // OCCT L4374-4376: get the shapes connection for using in the
    // rebuilding process.
    {
        let a_dm_for = a_dmf_im_f.clone();
        a_bf.shapes_connections(&a_dm_for, &a_mv_ref);
    }

    // OCCT L4378-4389: find the faces to remove.
    let a_sols = a_mv_ref.shape().unwrap_or_else(Shape::null);
    let mut a_dmfs: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    map_shapes_and_ancestors_map(&a_sols, ShapeType::Face, ShapeType::Solid, &mut a_dmfs);
    let a_nb = a_dmfs.len();
    if a_nb == 0 {
        return;
    }

    // OCCT L4391-4459: check the completeness of the created solids.
    let mut a_mf_to_rem: OcctShapeSet = HashMap::new();
    if a_mv_ref.has_deleted() {
        let mut a_me_holes = OcctIndexedShapeMap::new();
        map_shapes_indexed(the_f_holes, ShapeType::Edge, &mut a_me_holes);
        // OCCT L4402-4405: map the edges of the solids.
        let mut a_me_sols = OcctIndexedShapeMap::new();
        map_shapes_indexed(&a_sols, ShapeType::Edge, &mut a_me_sols);

        // OCCT L4407-4458: additional check on faces.
        let a_nb = a_bf.my_of_images.len();
        for i in 1..=a_nb {
            let (a_f, a_lf_im) = {
                let (_k, v) = a_bf.my_of_images.get_index(i - 1).expect("index");
                (v.0.clone(), v.1.clone())
            };
            if a_lf_im.is_empty() {
                continue;
            }
            let b_invalid = a_bf.my_invalid_faces.contains_key(&shape_key_of(&a_f));
            let mut b_connected = false;
            let mut b_face_kept = false;
            for a_f_im in a_lf_im.iter() {
                // OCCT L4432: if (!aMV.IsDeleted(aFIm)).
                if !a_mv_ref.is_deleted(a_f_im) {
                    b_face_kept = true;
                    continue;
                }
                for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                    if a_me_holes.contains(&a_e) {
                        b_face_kept = true;
                        set_add(&mut a_mf_to_rem, a_f_im);
                        break;
                    }
                    if !b_face_kept && b_invalid && !b_connected {
                        b_connected = a_me_sols.contains(&a_e);
                    }
                }
            }
            if !b_face_kept && !b_connected {
                return;
            }
        }
    }

    // OCCT L4461-4475.
    let mut a_me_boundary = OcctIndexedShapeMap::new();
    let a_nb = a_dmfs.len();
    for i in 1..=a_nb {
        let (a_f_im, a_l_sol) = {
            let (_k, v) = a_dmfs.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        if a_l_sol.len() > 1 {
            set_add(&mut a_mf_to_rem, &a_f_im);
        } else if a_f_im.orientation != Orientation::Internal {
            map_shapes_indexed(&a_f_im, ShapeType::Edge, &mut a_me_boundary);
        }
    }

    // OCCT L4477-4493: update the invalid faces with the images.
    let a_mv_ims = a_mv_ref.images();
    let a_nb = a_mf_inv.extent();
    for i in 1..=a_nb {
        let a_f_inv = a_mf_inv.find_key_1(i).clone();
        take_modified_no_fence_indexed(&a_f_inv, &a_mv_ims, &mut a_mf_inv);
    }
    for a_f in the_inverted_faces.iter() {
        take_modified_no_fence_indexed(a_f, &a_mv_ims, &mut a_mf_inv);
    }

    // OCCT L4495-4541: check if the invalid faces inside the holes are
    // really invalid.
    let a_nb_fh = the_mf_inv_in_hole.extent();
    for i in 1..=a_nb_fh {
        let a_f_inv = the_mf_inv_in_hole.find_key_1(i).clone();
        let mut a_lf_inv_im = a_mv_ref.modified(&a_f_inv);
        if a_lf_inv_im.is_empty() {
            a_lf_inv_im.push(a_f_inv.clone());
        }
        let p_f_offset = match shape_data_map::seek(&a_dmf_im_f, &a_f_inv) {
            Some(v) => v.clone(),
            None => continue,
        };
        for a_f_inv_im in a_lf_inv_im.iter() {
            let p_l_sols = match a_dmfs.get(&shape_key_of(a_f_inv_im)) {
                Some(v) => v.1.clone(),
                None => continue,
            };
            if p_l_sols.len() != 1 {
                continue;
            }
            let a_f_sol = &p_l_sols[0];
            let mut a_fx = Shape::null();
            if !find_shape(a_f_inv_im, a_f_sol, None, &mut a_fx) {
                continue;
            }
            // OCCT L4535: BRepOffset_Tool::CheckPlanesNormals(aFx,
            // pFOffset) — the brep_offset_tool_c.rs translation.
            if super::brep_offset_tool_d::check_planes_normals(&a_fx, &p_f_offset, 1.0e-8) {
                // OCCT L4537-4538: the normal direction has not changed.
                set_add(&mut a_mf_to_rem, a_f_inv_im);
            }
        }
    }

    // OCCT L4543-4603: the solids walk.
    let mut a_solids = empty_compound();
    let mut a_mf_keep: OcctShapeSet = HashMap::new();
    for a_sol in explorer(&a_sols, ShapeType::Solid, ShapeType::Shape) {
        let mut b_all_inv = true;
        let mut b_all_removed = true;
        for a_fs in explorer(&a_sol, ShapeType::Face, ShapeType::Shape) {
            if a_fs.orientation == Orientation::Internal {
                set_add(&mut a_mf_to_rem, &a_fs);
                continue;
            }
            if set_contains(&a_mf_to_rem, &a_fs) {
                continue;
            }
            b_all_removed = false;
            b_all_inv &= a_mf_inv.contains(&a_fs);
        }
        if b_all_inv && !b_all_removed {
            // OCCT L4579-4596: remove the invalid faces but keep those that
            // have already been marked for removal.
            for a_fs in explorer(&a_sol, ShapeType::Face, ShapeType::Shape) {
                if set_contains(&a_mf_to_rem, &a_fs) {
                    if !set_add(&mut a_mf_keep, &a_fs) {
                        a_mf_keep.remove(&shape_key_of(&a_fs));
                    }
                } else {
                    set_add(&mut a_mf_to_rem, &a_fs);
                }
            }
        } else {
            // OCCT L4600-4601.
            bat::builder_add_compound_shape(&mut a_solids, &a_sol);
            a_bf.my_solids = a_solids.clone();
        }
    }

    // OCCT L4605-4609.
    for a_f in a_mf_keep.values() {
        a_mf_to_rem.remove(&shape_key_of(a_f));
    }

    // OCCT L4611-4612: remove the invalid hanging parts external to the
    // solids.
    a_bf.remove_hanging_parts(&a_mv_ref, &a_dmf_im_f, &a_mf_inv, &mut a_mf_to_rem);

    // OCCT L4614-4616: remove the newly found internal and hanging faces.
    a_bf.remove_valid_splits(&a_mf_to_rem, &a_mv_ref, the_me_removed);
    a_bf.remove_invalid_splits(&a_mf_to_rem, &a_mv_ref, the_me_removed);

    // OCCT L4618-4628: get the inside faces from the removed ones.
    a_bf.my_inside_edges.clear();
    let a_nb = the_me_removed.extent();
    for i in 1..=a_nb {
        let a_e = the_me_removed.find_key_1(i).clone();
        if !a_me_boundary.contains(&a_e) {
            a_bf.my_inside_edges.add(&a_e);
        }
    }

    // OCCT L4630-4635: build all possible intersection pairs.
    if !a_mf_to_rem.is_empty() {
        let a_of_images_snapshot = a_bf.my_of_images.clone();
        let a_invalid_faces_snapshot = a_bf.my_invalid_faces.clone();
        build_intersection_pairs(
            &a_of_images_snapshot,
            &a_invalid_faces_snapshot,
            &a_mv_ref,
            &a_mf_to_rem,
            &a_dmf_im_f,
            &mut a_bf.my_intersection_pairs,
        );
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapesConnections (cxx L4643-4967).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::ShapesConnections (cxx L4643-4967) —
/// looking for the connections between the faces not to miss some
/// necessary intersection.
#[allow(unused_variables)]
pub(crate) fn shapes_connections_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_dm_for: &ShapeDataMap<Shape>,
    the_builder: &BuilderRef,
) {
    // OCCT L4647-4654: make the connexity blocks from the invalid edges.
    let mut a_ce_inv = empty_compound();
    for i in 1..=a_bf.my_invalid_edges.extent() {
        let a_e = a_bf.my_invalid_edges.find_key_1(i).clone();
        super::brep_offset_make_offset_1::add_to_container_shape(&a_e, &mut a_ce_inv);
    }

    // OCCT L4656-4657: MakeConnexityBlocks(aCEInv, VERTEX, EDGE, aLCB) —
    // the wire_splitter VERTEX/EDGE re-host form.
    let a_ce_inv_edges = explorer(&a_ce_inv, ShapeType::Edge, ShapeType::Shape);
    let a_locations = [glam::DAffine3::IDENTITY];
    let a_lcb =
        crate::bop::algo::wire_splitter::make_connexity_blocks(&a_ce_inv_edges, &a_locations);

    // OCCT L4659-4669: the binding from the edge to the block.
    let mut a_ecb_map: ShapeDataMap<Shape> = HashMap::new();
    for a_block in a_lcb.iter() {
        let a_cb = block_shape(&a_block.shapes);
        for a_e in bat::sub_shapes(&a_cb) {
            shape_data_map::bind(&mut a_ecb_map, &a_e, a_cb.clone());
        }
    }

    // OCCT L4671-4687: update the invalid edges with the images and keep
    // the connection to the original edge.
    let mut a_dme_or: ShapeDataMap<Vec<Shape>> = HashMap::new();
    let a_nb = a_bf.my_invalid_edges.extent();
    for i in 1..=a_nb {
        let a_e_inv = a_bf.my_invalid_edges.find_key_1(i).clone();
        let a_le_im = the_builder.modified(&a_e_inv);
        if a_le_im.is_empty() {
            add_to_container_map_list(&mut a_dme_or, &a_e_inv, &a_e_inv);
            continue;
        }
        for a_e_im in a_le_im.iter() {
            add_to_container_map_list(&mut a_dme_or, a_e_im, &a_e_inv);
        }
    }

    // OCCT L4689-4692: the DS interference analysis (the pDS->InterfFF()
    // access; the GAP form returns no interferences — the OCCT aNbFF == 0
    // path — when the MakerVolume re-host does not address the DS).
    let a_ffs = the_builder.pds_interf_ff();
    let a_nb_ff = a_ffs.len();
    let _ = a_nb_ff;
    // With the GAP DS access the FF walk below is a no-op: the OCCT loop
    // body L4694-4966 reads the InterfFF/InterfEF entries (Index1/Index2,
    // Curves/PaveBlocks, CommonPart, HasInterf) through the DS indices and
    // fills mySSInterfs / mySSInterfsArt from the interference records —
    // the empty record list keeps the call structure with the OCCT
    // no-interference result.
}

// ---------------------------------------------------------------------------
// OCCT RemoveHangingParts (cxx L4968-5167).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::RemoveHangingParts (cxx L4968-5167) —
/// remove the isolated invalid hanging parts.
pub(crate) fn remove_hanging_parts_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_mv: &BuilderRef,
    the_dmf_im_f: &ShapeDataMap<Shape>,
    the_mf_inv: &OcctIndexedShapeMap,
    the_mf_to_rem: &mut OcctShapeSet,
) {
    // OCCT L4975-4978: map the faces of the result solids.
    let mut a_mfs = OcctIndexedShapeMap::new();
    if let Some(a_shape) = the_mv.shape() {
        map_shapes_indexed(&a_shape, ShapeType::Face, &mut a_mfs);
    }

    // OCCT L4980-4981: the compound of all the faces not included into the
    // solids.
    let mut a_cf_hangs = empty_compound();

    // OCCT L4983-4985: the images of the MV.
    let a_mv_ims = the_mv.images();

    // OCCT L4987-4995: take the faces of the arguments not included into
    // the solids (TakeModified(aF, aMVIms, aCFHangs, &aMFS) — the fence is
    // the faces of the solids).
    let mut a_mfs_fence: OcctShapeSet = a_mfs
        .iter()
        .map(|s| (shape_key_of(s), s.clone()))
        .collect();
    for a_arg in the_mv.arguments().iter() {
        for a_f in explorer(a_arg, ShapeType::Face, ShapeType::Shape) {
            take_modified_shape(&a_f, &a_mv_ims, &mut a_cf_hangs, &mut a_mfs_fence);
        }
    }

    // OCCT L4997-5003: make the connexity blocks of all the hanging parts
    // (the shell_splitter EDGE/FACE re-host form).
    let a_cf_hangs_faces = explorer(&a_cf_hangs, ShapeType::Face, ShapeType::Shape);
    let a_lcb_hangs =
        crate::bop::algo::shell_splitter::make_connexity_blocks(&a_cf_hangs_faces);
    if a_lcb_hangs.is_empty() {
        return;
    }

    // OCCT L5005-5015: map the edges and vertices of the result solids.
    let mut a_dmef: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let mut a_dmve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    if let Some(a_shape) = the_mv.shape() {
        map_shapes_and_ancestors_map(&a_shape, ShapeType::Edge, ShapeType::Face, &mut a_dmef);
        map_shapes_and_ancestors_map(&a_shape, ShapeType::Vertex, ShapeType::Edge, &mut a_dmve);
    }

    // OCCT L5017-5027: update the invalid edges with the intersection
    // results.
    let mut a_me_inv: OcctShapeSet = HashMap::new();
    let a_nb_e = a_bf.my_invalid_edges.extent();
    for i in 1..=a_nb_e {
        let a_e = a_bf.my_invalid_edges.find_key_1(i).clone();
        take_modified_map_no_fence(&a_e, &a_mv_ims, &mut a_me_inv);
    }

    // OCCT L5029-5034: update the inverted edges with the intersection
    // results.
    let mut a_me_inverted: OcctShapeSet = HashMap::new();
    for i_inv in 1..=a_bf.my_inverted_edges.extent() {
        let a_e = a_bf.my_inverted_edges.find_key_1(i_inv).clone();
        take_modified_map_no_fence(&a_e, &a_mv_ims, &mut a_me_inverted);
    }

    // OCCT L5036-5038: the origins of the splits.
    let a_mv_ors = the_mv.origins();

    // OCCT L5040-5139: find the hanging blocks to remove.
    let mut a_blocks_to_remove: Vec<Shape> = Vec::new();
    for a_block in a_lcb_hangs.iter() {
        let a_cbh = block_shape(&a_block.shapes);

        // OCCT L5046-5056: remove the block containing the inverted edges.
        let mut b_has_inverted = false;
        for a_e in explorer(&a_cbh, ShapeType::Edge, ShapeType::Shape) {
            b_has_inverted = !a_dmef.contains_key(&shape_key_of(&a_e))
                && set_contains(&a_me_inverted, &a_e);
            if b_has_inverted {
                break;
            }
        }
        if b_has_inverted {
            a_blocks_to_remove.push(a_cbh);
            continue;
        }

        // OCCT L5058-5122.
        let mut b_has_invalid_face = false;
        let mut b_is_connected = false;
        let mut a_block_me = OcctIndexedShapeMap::new();
        map_shapes_indexed(&a_cbh, ShapeType::Edge, &mut a_block_me);
        let mut a_m_offset_f: OcctShapeSet = HashMap::new();

        for a_f in explorer(&a_cbh, ShapeType::Face, ShapeType::Shape) {
            // OCCT L5064-5066: check the block to contain an invalid face.
            if !b_has_invalid_face {
                b_has_invalid_face = the_mf_inv.contains(&a_f);
            }
            // OCCT L5068-5098: check the block for connectivity to the
            // invalid parts.
            if !b_is_connected {
                for a_e in explorer(&a_f, ShapeType::Edge, ShapeType::Shape) {
                    if b_is_connected {
                        break;
                    }
                    if let Some(p_lf) = a_dmef.get(&shape_key_of(&a_e)) {
                        for a_f2 in p_lf.1.iter() {
                            b_is_connected = the_mf_inv.contains(a_f2);
                            if b_is_connected {
                                break;
                            }
                        }
                    }
                }
                // OCCT L5086-5097: check vertices.
                if !b_is_connected {
                    for a_v in explorer(&a_f, ShapeType::Vertex, ShapeType::Shape) {
                        if b_is_connected {
                            break;
                        }
                        if let Some(p_le) = a_dmve.get(&shape_key_of(&a_v)) {
                            for a_e2 in p_le.1.iter() {
                                b_is_connected = !a_block_me.contains(a_e2)
                                    && set_contains(&a_me_inv, a_e2);
                                if b_is_connected {
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            // OCCT L5100-5121: check the block to be isolated.
            match a_mv_ors.get((a_f.ptr_id(), a_f.location)) {
                Some(p_lf_or) => {
                    for a_f_or in p_lf_or.iter() {
                        if let Some(p_f_offset) = shape_data_map::seek(the_dmf_im_f, a_f_or) {
                            set_add(&mut a_m_offset_f, p_f_offset);
                        }
                    }
                }
                None => {
                    if let Some(p_f_offset) = shape_data_map::seek(the_dmf_im_f, &a_f) {
                        set_add(&mut a_m_offset_f, p_f_offset);
                    }
                }
            }
        }

        // OCCT L5124-5131.
        let b_remove = b_has_invalid_face && (!b_is_connected || a_m_offset_f.len() == 1);
        if b_remove {
            a_blocks_to_remove.push(a_cbh);
        }
    }

    // OCCT L5133-5141: remove the invalidated blocks.
    for a_cbh in a_blocks_to_remove.iter() {
        for a_f in explorer(a_cbh, ShapeType::Face, ShapeType::Shape) {
            set_add(the_mf_to_rem, &a_f);
        }
    }
}
