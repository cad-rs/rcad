// OCCT BRepOffset_MakeOffset_1.cxx L5168-5923 — module f of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L5168-5923:
//   - RemoveValidSplits (cxx L5168-5227)
//   - RemoveInvalidSplits (cxx L5228-5319)
//   - FilterEdgesImages (cxx L5320-5368)
//   - FilterInvalidFaces (cxx L5369-5554)
//   - CheckEdgesCreatedByVertex (cxx L5555-5605)
//   - FilterInvalidEdges (cxx L5606-5757)
//   - FindFacesToRebuild (cxx L5758-5905)
//   - mapShapes (cxx L5906-5921, file static template)
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#49).

use std::collections::HashMap;

use rcad_kernel::topo::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    append_to_list, block_shape, empty_compound, map_shapes_indexed,
    shape_key_of, BRepOffsetBuildOffsetFaces, BuilderRef, DataMapOfShapeIndexedMapOfShape,
};
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, shape_indexed_data_map, OcctIndexedShapeMap,
    OcctShapeSet, ShapeDataMap, ShapeIndexedDataMap,
};
use crate::feat::loc_ope_wires_on_shape_b::ShapeKey;

// ---------------------------------------------------------------------------
// OCCT RemoveValidSplits (cxx L5168-5227).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::RemoveValidSplits (cxx L5168-5227) —
/// removing the valid splits according to the results of the intersection.
pub(crate) fn remove_valid_splits_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_sp_rem: &OcctShapeSet,
    the_gf: &BuilderRef,
    the_me_removed: &mut OcctIndexedShapeMap,
) {
    // OCCT L5170-5174.
    let a_nb = a_bf.my_of_images.len();
    if a_nb == 0 {
        return;
    }

    // OCCT L5176-5226: the iterator-removal walk over myOFImages.
    for i in 1..=a_nb {
        let (a_key, mut a_ls_im) = {
            let (_, v) = a_bf.my_of_images.get_index(i - 1).expect("index");
            (shape_key_of(&v.0), v.1.clone())
        };
        let mut j = 0usize;
        while j < a_ls_im.len() {
            let a_s_im = a_ls_im[j].clone();
            if set_contains(the_sp_rem, &a_s_im) {
                map_shapes_indexed(&a_s_im, ShapeType::Edge, the_me_removed);
                a_ls_im.remove(j);
                continue;
            }
            // OCCT L5187-5208: check if all its images have to be removed.
            let a_ls_im_im = the_gf.modified(&a_s_im);
            if !a_ls_im_im.is_empty() {
                let mut b_all_rem = true;
                for a_s_im_im in a_ls_im_im.iter() {
                    if set_contains(the_sp_rem, a_s_im_im) {
                        map_shapes_indexed(a_s_im_im, ShapeType::Edge, the_me_removed);
                    } else {
                        b_all_rem = false;
                    }
                }
                if b_all_rem {
                    map_shapes_indexed(&a_s_im, ShapeType::Edge, the_me_removed);
                    a_ls_im.remove(j);
                    continue;
                }
            }
            // OCCT L5211: aIt.Next().
            j += 1;
        }
        if let Some(entry) = a_bf.my_of_images.get_mut(&a_key) {
            entry.1 = a_ls_im;
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT RemoveInvalidSplits (cxx L5228-5319).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::RemoveInvalidSplits (cxx L5228-5319) —
/// removing the invalid splits according to the results of the
/// intersection.
pub(crate) fn remove_invalid_splits_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_sp_rem: &OcctShapeSet,
    the_gf: &BuilderRef,
    the_me_removed: &mut OcctIndexedShapeMap,
) {
    // OCCT L5230-5234.
    let a_nb = a_bf.my_invalid_faces.len();
    if a_nb == 0 {
        return;
    }

    // OCCT L5236-5318.
    for i in 1..=a_nb {
        let (a_s, mut a_ls_im) = {
            let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        // OCCT L5237: bArt = myArtInvalidFaces.IsBound(aS).
        let b_art = shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_s);

        let mut j = 0usize;
        while j < a_ls_im.len() {
            let a_s_im = a_ls_im[j].clone();
            if set_contains(the_sp_rem, &a_s_im) {
                map_shapes_indexed(&a_s_im, ShapeType::Edge, the_me_removed);
                a_ls_im.remove(j);
                continue;
            }
            // OCCT L5250-5257.
            let a_ls_im_im = the_gf.modified(&a_s_im);
            if a_ls_im_im.is_empty() {
                j += 1;
                continue;
            }
            // OCCT L5259-5271.
            let mut b_all_rem = true;
            let mut a_me_removed_loc = OcctIndexedShapeMap::new();
            for a_s_im_im in a_ls_im_im.iter() {
                if set_contains(the_sp_rem, a_s_im_im) {
                    map_shapes_indexed(a_s_im_im, ShapeType::Edge, &mut a_me_removed_loc);
                } else {
                    b_all_rem = false;
                }
            }
            if b_all_rem {
                a_ls_im.remove(j);
                continue;
            }
            // OCCT L5273-5276.
            if b_art {
                j += 1;
                continue;
            }
            // OCCT L5278-5291: remove the face from the invalid ones if all
            // the invalid edges of this face have been marked for removal.
            let mut b_broke = false;
            for a_e_inv in explorer(&a_s_im, ShapeType::Edge, ShapeType::Shape) {
                if a_bf.my_invalid_edges.contains(&a_e_inv)
                    && !a_me_removed_loc.contains(&a_e_inv)
                {
                    b_broke = true;
                    break;
                }
            }
            if !b_broke {
                map_shapes_indexed(&a_s_im, ShapeType::Edge, the_me_removed);
                a_ls_im.remove(j);
            } else {
                j += 1;
            }
        }
        if let Some(entry) = a_bf.my_invalid_faces.get_mut(&shape_key_of(&a_s)) {
            entry.1 = a_ls_im;
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT FilterEdgesImages (cxx L5320-5368).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FilterEdgesImages (cxx L5320-5368) —
/// updating the maps of images and origins of the offset edges.
pub(crate) fn filter_edges_images_impl(a_bf: &mut BRepOffsetBuildOffsetFaces, the_s: &Shape) {
    // OCCT L5323-5325: map edges.
    let mut a_me = OcctIndexedShapeMap::new();
    map_shapes_indexed(the_s, ShapeType::Edge, &mut a_me);

    // OCCT L5327: myOEOrigins.Clear().
    a_bf.my_oe_origins.clear();

    // OCCT L5328-5358: walk the myOEImages entries.
    let a_keys: Vec<ShapeKey> = a_bf.my_oe_images.keys().cloned().collect();
    for a_key in a_keys.iter() {
        let (a_e, a_le_im) = match a_bf.my_oe_images.get(a_key) {
            Some(v) => (v.0.clone(), v.1.clone()),
            None => continue,
        };
        // OCCT L5335-5352: the iterator-removal walk.
        let mut kept: Vec<Shape> = Vec::new();
        for a_e_im in a_le_im.iter() {
            // OCCT L5340-5346: filter images — the edges with no images
            // left should be kept in the map to avoid their usage when
            // building the splits of faces.
            if !a_me.contains(a_e_im) {
                continue;
            }
            kept.push(a_e_im.clone());
            // OCCT L5349-5355: save origins.
            if shape_data_map::is_bound(&a_bf.my_oe_origins, a_e_im) {
                append_to_list(
                    shape_data_map::change_find(&mut a_bf.my_oe_origins, a_e_im),
                    &a_e,
                );
            } else {
                let a_l_or = vec![a_e.clone()];
                shape_data_map::bind(&mut a_bf.my_oe_origins, a_e_im, a_l_or);
            }
        }
        if let Some(entry) = a_bf.my_oe_images.get_mut(a_key) {
            entry.1 = kept;
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT FilterInvalidFaces (cxx L5369-5554).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FilterInvalidFaces (cxx L5369-5554) —
/// filtering of the invalid faces: faces having only valid images left with
/// non-free edges are considered valid; the invalid faces that would create
/// free edges are not removed.
pub(crate) fn filter_invalid_faces_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_dmef: &ShapeIndexedDataMap<Vec<Shape>>,
    the_me_removed: &OcctIndexedShapeMap,
) {
    // OCCT L5376-5380.
    let mut a_really_inv_faces: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    // OCCT L5377-5379: Edge-Face connexity map of all the splits, both
    // invalid and valid.
    let mut a_dmef_all: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();

    let a_nb = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb {
        let (a_f, a_lf_inv) = {
            let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };

        // OCCT L5385-5397.
        if shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f) {
            if a_lf_inv.is_empty() {
                shape_data_map::un_bind(&mut a_bf.my_art_invalid_faces, &a_f);
            } else {
                indexed_add(&mut a_really_inv_faces, &a_f, a_lf_inv.clone());
            }
            continue;
        }
        // OCCT L5399-5402.
        if a_lf_inv.is_empty() {
            continue;
        }

        // OCCT L5404-5406: aLFIm = myOFImages.ChangeFromKey(aF);
        // bInvalid = aLFIm.IsEmpty().
        let a_key = shape_key_of(&a_f);
        let mut a_lf_im: Vec<Shape> = match a_bf.my_of_images.get(&a_key) {
            Some(v) => v.1.clone(),
            None => Vec::new(),
        };
        let mut b_invalid = a_lf_im.is_empty();

        // OCCT L5408-5430: check the two lists on common splits.
        if !b_invalid {
            let mut b_broke_outer = false;
            for a_f_inv in a_lf_inv.iter() {
                let mut b_broke_inner = false;
                for a_f_im in a_lf_im.iter() {
                    if a_f_inv.is_same(a_f_im) {
                        b_broke_inner = true;
                        break;
                    }
                }
                if b_broke_inner {
                    b_broke_outer = true;
                    break;
                }
            }
            b_invalid = b_broke_outer;
        }

        // OCCT L5432-5456: check for free edges.
        if !b_invalid {
            for j in 0..2 {
                if b_invalid {
                    break;
                }
                let a_li = if j == 0 { &a_lf_im } else { &a_lf_inv };
                let mut b_broke_outer = false;
                for a_f_im in a_li.iter() {
                    let mut b_broke_inner = false;
                    for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                        if !the_me_removed.contains(&a_e) {
                            let p_lef = the_dmef
                                .get(&shape_key_of(&a_e))
                                .map(|e| &e.1);
                            if let Some(p_lef) = p_lef {
                                if p_lef.len() == 1 {
                                    b_broke_inner = true;
                                    break;
                                }
                            }
                        }
                    }
                    if b_broke_inner {
                        b_broke_outer = true;
                        break;
                    }
                }
                b_invalid = b_broke_outer;
            }
        }

        // OCCT L5458-5502.
        if b_invalid {
            // OCCT L5460-5470: fill the aDMEFAll copy once.
            if a_dmef_all.is_empty() {
                a_dmef_all = the_dmef.clone();
                for i_f in 1..=a_nb {
                    let lf_inv = {
                        let (_, v) = a_bf
                            .my_invalid_faces
                            .get_index(i_f - 1)
                            .expect("index");
                        v.1.clone()
                    };
                    for a_f_inv in lf_inv.iter() {
                        super::brep_offset_make_offset_1::map_shapes_and_ancestors_map(
                            a_f_inv,
                            ShapeType::Edge,
                            ShapeType::Face,
                            &mut a_dmef_all,
                        );
                    }
                }
            }

            // OCCT L5472-5479: the local splits.
            let mut a_local_splits: OcctShapeSet = HashMap::new();
            for j in 0..2 {
                let a_li = if j == 0 { &a_lf_im } else { &a_lf_inv };
                for a_f in a_li.iter() {
                    set_add(&mut a_local_splits, a_f);
                }
            }

            // OCCT L5481-5500: check if all the invalid edges are located
            // inside the split and do not touch any other faces.
            let mut b_broke_outer = false;
            for a_f_im in a_lf_inv.iter() {
                let mut b_broke_inner = false;
                for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                    if a_bf.my_invalid_edges.contains(&a_e) && !the_me_removed.contains(&a_e) {
                        let a_lf = shape_indexed_data_map::find(&a_dmef_all, &a_e);
                        let mut b_broke_lf = false;
                        for a_f2 in a_lf.iter() {
                            if !set_contains(&a_local_splits, a_f2) {
                                b_broke_lf = true;
                                break;
                            }
                        }
                        if b_broke_lf {
                            b_broke_inner = true;
                            break;
                        }
                    }
                }
                if b_broke_inner {
                    b_broke_outer = true;
                    break;
                }
            }
            b_invalid = b_broke_outer;

            // OCCT L5501-5504: if (!bInvalid) — keep the images.
            if !b_invalid {
                for a_f_inv in a_lf_inv.iter() {
                    a_lf_im.push(a_f_inv.clone());
                }
                add_of_images_f(&mut a_bf.my_of_images, &a_key, a_lf_im);
            }
        }

        // OCCT L5506-5509.
        if b_invalid {
            indexed_add(&mut a_really_inv_faces, &a_f, a_lf_inv);
        }
    }

    // OCCT L5512: myInvalidFaces = aReallyInvFaces.
    a_bf.my_invalid_faces = a_really_inv_faces;
}

/// The IndexedDataMap::Add replace-on-Add form (the module-local form of
/// the OCCT aReallyInvFaces.Add).
fn indexed_add(m: &mut ShapeIndexedDataMap<Vec<Shape>>, k: &Shape, v: Vec<Shape>) {
    let key = shape_key_of(k);
    if let Some(entry) = m.get_mut(&key) {
        entry.1 = v;
    } else {
        m.insert(key, (k.clone(), v));
    }
}

/// The myOFImages write-back of the cxx L5404 ChangeFromKey form.
fn add_of_images_f(
    m: &mut ShapeIndexedDataMap<Vec<Shape>>,
    key: &ShapeKey,
    v: Vec<Shape>,
) {
    if let Some(entry) = m.get_mut(key) {
        entry.1 = v;
    }
}

// ---------------------------------------------------------------------------
// OCCT CheckEdgesCreatedByVertex (cxx L5555-5605).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::CheckEdgesCreatedByVertex (cxx
/// L5555-5605) — checks additionally the unchecked edges originated from
/// the vertices.
pub(crate) fn check_edges_created_by_vertex_impl(a_bf: &mut BRepOffsetBuildOffsetFaces) {
    let my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();

    // OCCT L5559-5592.
    let a_nb_f = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb_f {
        let (a_f, a_lf_im) = {
            let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        if shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f) {
            continue;
        }
        for a_f_im in a_lf_im.iter() {
            for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                if a_bf.my_invalid_edges.contains(&a_e) || a_bf.my_valid_edges.contains(&a_e) {
                    continue;
                }
                // OCCT L5575-5581: check if this edge is not created from a
                // vertex and mark it as invalid.
                let p_le_or = match shape_data_map::seek(&my_edges_origins, &a_e) {
                    Some(v) => v.clone(),
                    None => continue,
                };
                let mut b_broke = false;
                for a_s in p_le_or.iter() {
                    if a_s.shape_type() != ShapeType::Vertex {
                        b_broke = true;
                        break;
                    }
                }
                // OCCT L5585-5588: if (!itLEO.More()) myInvalidEdges.Add.
                if !b_broke {
                    a_bf.my_invalid_edges.add(&a_e);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT FilterInvalidEdges (cxx L5606-5757).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FilterInvalidEdges (cxx L5606-5757) —
/// filtering the invalid edges according to the currently invalid faces.
pub(crate) fn filter_invalid_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_dmfmie: &DataMapOfShapeIndexedMapOfShape,
    the_me_removed: &OcctIndexedShapeMap,
    the_me_inside: &OcctIndexedShapeMap,
    the_me_use_in_rebuild: &mut OcctShapeSet,
) {
    let my_oe_origins = a_bf.my_oe_origins.clone();
    let my_oe_images = a_bf.my_oe_images.clone();

    // OCCT L5561-5566 (cxx numbering L5612-5616).
    let mut a_ce_inv = empty_compound();
    let mut a_me_inv = OcctIndexedShapeMap::new();

    // OCCT L5619-5636.
    let a_nb = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb {
        let a_lf_inv = {
            let (_, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            v.1.clone()
        };
        for a_f_im in a_lf_inv.iter() {
            map_shapes_indexed(a_f_im, ShapeType::Edge, &mut a_me_inv);
            for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                if a_bf.my_invalid_edges.contains(&a_e) {
                    bat::builder_add_compound_shape(&mut a_ce_inv, &a_e);
                }
            }
        }
    }

    // OCCT L5638-5658: remove the edges which have been marked for removal.
    let mut a_me_inv_to_avoid = OcctIndexedShapeMap::new();
    let a_ce_inv_edges = explorer(&a_ce_inv, ShapeType::Edge, ShapeType::Shape);
    let a_locations = [glam::DAffine3::IDENTITY];
    let a_lcbe = crate::bop::algo::wire_splitter::make_connexity_blocks(
        &a_ce_inv_edges,
        &a_locations,
    );
    for a_block in a_lcbe.iter() {
        let a_cbe = block_shape(&a_block.shapes);
        let mut b_broke = false;
        for a_e in explorer(&a_cbe, ShapeType::Edge, ShapeType::Shape) {
            if !the_me_removed.contains(&a_e) {
                b_broke = true;
                break;
            }
        }
        if !b_broke {
            map_shapes_indexed(&a_cbe, ShapeType::Edge, &mut a_me_inv_to_avoid);
        }
    }

    // OCCT L5660-5700: the really invalid edges.
    let mut a_really_inv_edges = OcctIndexedShapeMap::new();
    let a_nb = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb {
        let (a_f, a_lf_inv) = {
            let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        if shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f) {
            if let Some(a_mie) = shape_data_map::seek(the_dmfmie, &a_f) {
                let a_nb_ie = a_mie.extent();
                for i_e in 1..=a_nb_ie {
                    let a_e = a_mie.find_key_1(i_e).clone();
                    if a_me_inv.contains(&a_e) && !a_me_inv_to_avoid.contains(&a_e) {
                        a_really_inv_edges.add(&a_e);
                    }
                }
            }
        } else {
            for a_f_im in a_lf_inv.iter() {
                for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                    if a_bf.my_invalid_edges.contains(&a_e)
                        && !a_me_inv_to_avoid.contains(&a_e)
                    {
                        a_really_inv_edges.add(&a_e);
                    }
                }
            }
        }
    }

    // OCCT L5702: myInvalidEdges = aReallyInvEdges.
    a_bf.my_invalid_edges = a_really_inv_edges;

    // OCCT L5704-5756: check if any of the currently invalid edges may be
    // used for rebuilding the splits of the invalid faces.
    let a_nb = a_bf.my_invalid_edges.extent();
    for i in 1..=a_nb {
        let a_e = a_bf.my_invalid_edges.find_key_1(i).clone();
        if !the_me_inside.contains(&a_e) || !a_bf.my_valid_edges.contains(&a_e) {
            continue;
        }
        let p_e_origins = match shape_data_map::seek(&my_oe_origins, &a_e) {
            Some(v) => v.clone(),
            None => continue,
        };
        let mut b_has_inv_outside = false;
        for a_e_or in p_e_origins.iter() {
            if b_has_inv_outside {
                break;
            }
            let p_e_ims = shape_data_map::seek(&my_oe_images, a_e_or).cloned();
            if let Some(p_e_ims) = p_e_ims {
                for a_e_im in p_e_ims.iter() {
                    b_has_inv_outside = a_bf.my_invalid_edges.contains(a_e_im)
                        && !the_me_inside.contains(a_e_im);
                    if b_has_inv_outside {
                        break;
                    }
                }
            }
        }
        if !b_has_inv_outside {
            set_add(the_me_use_in_rebuild, &a_e);
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT FindFacesToRebuild (cxx L5758-5905).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FindFacesToRebuild (cxx L5758-5905) —
/// looking for the faces that have to be rebuilt: 1. faces close to
/// invalidity; 2. faces containing some invalid parts.
pub(crate) fn find_faces_to_rebuild_impl(a_bf: &mut BRepOffsetBuildOffsetFaces) {
    // OCCT L5760-5762.
    let a_nb = a_bf.my_of_images.len();
    if a_nb == 0 {
        return;
    }

    let mut b_rebuild: bool;
    let mut a_le_valid: Vec<Shape> = Vec::new();
    let mut a_m_fence: OcctShapeSet = HashMap::new();
    let mut a_me_reb: OcctShapeSet = HashMap::new();
    let mut a_mf_reb: OcctShapeSet = HashMap::new();

    // OCCT L5769: aDMFLV.
    let mut a_dmflv: ShapeDataMap<Vec<Shape>> = HashMap::new();

    // OCCT L5771-5802: get the edges from the invalid faces.
    let a_nb = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb {
        let (a_f, a_lf_im) = {
            let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        a_m_fence.clear();
        let mut a_lv_avoid: Vec<Shape> = Vec::new();
        for a_f_im in a_lf_im.iter() {
            for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                set_add(&mut a_me_reb, &a_e);
                if a_bf.my_invalid_edges.contains(&a_e) {
                    for a_v in explorer(&a_e, ShapeType::Vertex, ShapeType::Shape) {
                        if set_add(&mut a_m_fence, &a_v) {
                            a_lv_avoid.push(a_v.clone());
                            set_add(&mut a_me_reb, &a_v);
                        }
                    }
                }
            }
        }
        // OCCT L5790-5793.
        if !a_lv_avoid.is_empty() {
            shape_data_map::bind(&mut a_dmflv, &a_f, a_lv_avoid);
        }
        // OCCT L5795-5802.
        let p_lf = if !shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f) {
            shape_data_map::seek(&a_bf.my_ss_interfs, &a_f)
        } else {
            shape_data_map::seek(&a_bf.my_ss_interfs_art, &a_f)
        };
        if let Some(p_lf) = p_lf {
            let p_lf = p_lf.clone();
            for a_fe in p_lf.iter() {
                set_add(&mut a_mf_reb, a_fe);
            }
        }
    }

    // OCCT L5804-5879: get the faces to rebuild.
    let a_nb = a_bf.my_of_images.len();
    for i in 1..=a_nb {
        let (a_f, a_lf_im) = {
            let (_k, v) = a_bf.my_of_images.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        let mut a_mv_avoid: OcctShapeSet = HashMap::new();
        if shape_data_map::is_bound(&a_dmflv, &a_f) {
            let a_lv_avoid = shape_data_map::find(&a_dmflv, &a_f);
            for a_v in a_lv_avoid.iter() {
                set_add(&mut a_mv_avoid, a_v);
            }
        }

        b_rebuild = set_contains(&a_mf_reb, &a_f);
        a_le_valid.clear();
        a_m_fence.clear();

        for a_f_im in a_lf_im.iter() {
            for an_e_im in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                if !a_bf.my_invalid_edges.contains(&an_e_im) {
                    if set_add(&mut a_m_fence, &an_e_im) {
                        a_le_valid.push(an_e_im.clone());
                    }
                }
                // OCCT L5833-5836.
                if !b_rebuild {
                    b_rebuild = set_contains(&a_me_reb, &an_e_im);
                }
                // OCCT L5838-5847: check vertices.
                if !b_rebuild {
                    for a_v in explorer(&an_e_im, ShapeType::Vertex, ShapeType::Shape) {
                        if b_rebuild {
                            break;
                        }
                        if !set_contains(&a_mv_avoid, &a_v) {
                            b_rebuild = set_contains(&a_me_reb, &a_v);
                        }
                    }
                }
            }
        }

        // OCCT L5850-5858.
        if !b_rebuild {
            b_rebuild = !a_lf_im.is_empty() && a_bf.my_invalid_faces.contains_key(&shape_key_of(&a_f));
            if b_rebuild {
                set_add(&mut a_bf.my_f_self_reb_avoid, &a_f);
            }
        }

        // OCCT L5860-5863.
        if b_rebuild {
            // OCCT L5862: myFacesToRebuild.Add(aF, aLEValid) — the
            // IndexedDataMap::Add replace-on-Add form.
            let key = shape_key_of(&a_f);
            if let Some(entry) = a_bf.my_faces_to_rebuild.get_mut(&key) {
                entry.1 = a_le_valid.clone();
            } else {
                a_bf.my_faces_to_rebuild.insert(key, (a_f.clone(), a_le_valid.clone()));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT mapShapes (cxx L5906-5921, file static template).
// ---------------------------------------------------------------------------

/// OCCT static mapShapes (cxx L5906-5921) — collect theVecShapes into
/// theMap with theType set (the template is taken in the shape-list form
/// used by the cxx call sites).
pub(crate) fn map_shapes_static(the_vec_shapes: &[Shape], the_type: ShapeType, the_map: &mut OcctShapeSet) {
    for a_shape in the_vec_shapes.iter() {
        for an_exp in explorer(a_shape, the_type, ShapeType::Shape) {
            set_add(the_map, &an_exp);
        }
    }
}
