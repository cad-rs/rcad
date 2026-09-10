// OCCT BRepOffset_MakeOffset_1.cxx L8551-9533 — module j of the 1:1
// translation (the final split of brep_offset_make_offset_1.rs).
//
// This module carries cxx L8551-9533:
//   - GetBounds (cxx L8551-8593)
//   - GetBoundsToUpdate (cxx L8594-8674)
//   - GetInvalidEdgesByBounds (cxx L8675-8996)
//   - FilterSplits (cxx L8997-9047)
//   - UpdateNewIntersectionEdges (cxx L9048-9296)
//   - FillGaps (cxx L9297-9412)
//   - FillHistory (cxx L9413-9491)
//   - BRepOffset_MakeOffset::BuildSplitsOfTrimmedFaces /
//     BuildSplitsOfExtendedFaces dispatchers (cxx L9497-9533)
//
// Local dependency carriers of this module:
//   - IntTools_Context::ComputeVF / IsPointInOnFace (TKBO IntTools) — the
//     GAP carriers with the OCCT error paths (iStatus != 0 skips the
//     classification; the classification result keeps the OCCT structure);
//     the rcad IntToolsContext is DS-bound and the offset pipeline works
//     on the BRep shapes.
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#49).

use std::collections::HashMap;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::geom::{CurveEval, Curve3};
use rcad_kernel::topo::topods::{Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::algo::section::BOPAlgoSection;
use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    add_to_container_shape, append_to_list, build_splits_of_face, empty_compound,
    map_shapes_and_ancestors_map, map_shapes_indexed, shape_key_of,
    BRepOffsetBuildOffsetFaces,
};
use crate::feat::loc_ope_wires_on_shape_b::ShapeKey;
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, OcctIndexedShapeMap, OcctShapeSet,
    ShapeDataMap, ShapeIndexedDataMap,
};

// ---------------------------------------------------------------------------
// Local dependency carriers.
// ---------------------------------------------------------------------------

/// OCCT IntTools_Context::ComputeVF(V, F, U, V, aTol) (IntTools_Context.cxx)
/// — GAP carrier: the rcad IntToolsContext re-host is DS-bound and exposes
/// no BRep-shape ComputeVF; the OCCT error path is taken (iStatus != 0 —
/// the cxx classification branch is skipped).
fn compute_vf_gap(_the_v: &Shape, _the_f: &Shape) -> i32 {
    1
}

/// OCCT IntTools_Context::IsPointInOnFace(F, P2d) (IntTools_Context.cxx) —
/// GAP carrier with the OCCT false path (the point is not ON the face).
fn is_point_in_on_face_gap(_the_f: &Shape, _the_p2d: glam::DVec2) -> bool {
    false
}

// ---------------------------------------------------------------------------
// OCCT GetBounds (cxx L8551-8593).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::GetBounds (cxx L8551-8593) — getting
/// the edges from the splits of the offset faces.
pub(crate) fn get_bounds_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_lfaces: &[Shape],
    the_meb: &OcctShapeSet,
    the_bounds: &mut Shape,
) {
    // OCCT L8556-8559: the compound of the edges contained in the splits of
    // the faces.
    let mut a_bounds = empty_compound();
    // OCCT L8560: fence map.
    let mut a_m_fence: OcctShapeSet = HashMap::new();

    for a_f in the_lfaces.iter() {
        let p_lf_im = match super::brep_offset_make_offset_1::idm_seek(&a_bf.my_of_images, a_f) {
            Some(v) => v,
            None => continue,
        };
        for a_f_im in p_lf_im.iter() {
            for a_e_im in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                if !set_contains(the_meb, &a_e_im) && set_add(&mut a_m_fence, &a_e_im) {
                    add_to_container_shape(&a_e_im, &mut a_bounds);
                }
            }
        }
    }
    // OCCT L8585: theBounds = aBounds.
    *the_bounds = a_bounds;
}

// ---------------------------------------------------------------------------
// OCCT GetBoundsToUpdate (cxx L8594-8674).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::GetBoundsToUpdate (cxx L8594-8674) —
/// get the bounding edges that should be updated.
pub(crate) fn get_bounds_to_update_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_lf: &[Shape],
    the_meb: &OcctShapeSet,
    the_la_bounds: &mut Vec<Shape>,
    the_la_valid: &mut Vec<Shape>,
    the_bounds: &mut Shape,
) {
    let my_oe_images = a_bf.my_oe_images.clone();
    let my_oe_origins = a_bf.my_oe_origins.clone();

    // OCCT L8599-8601: get all the edges.
    let mut a_bounds = empty_compound();
    // OCCT L8603.
    let mut a_ma_valid: OcctShapeSet = HashMap::new();
    let mut a_m_fence: OcctShapeSet = HashMap::new();

    for a_f in the_lf.iter() {
        // OCCT L8606-8622: the descendants of the face with their images.
        let mut a_mde = OcctIndexedShapeMap::new();
        let a_lf_des: Vec<Shape> = a_bf
            .my_as_des
            .as_ref()
            .expect("myAsDes")
            .borrow()
            .descendant(a_f)
            .to_vec();
        for a_ed in a_lf_des.iter() {
            match shape_data_map::seek(&my_oe_images, a_ed) {
                None => {
                    a_mde.add(a_ed);
                }
                Some(p_led_im) => {
                    let p_led_im = p_led_im.clone();
                    for a_ed_im in p_led_im.iter() {
                        a_mde.add(a_ed_im);
                    }
                }
            }
        }
        // OCCT L8624-8657.
        let a_nb_e = a_mde.extent();
        for j in 1..=a_nb_e {
            let a_e_im = a_mde.find_key_1(j).clone();
            if !set_contains(the_meb, &a_e_im) && set_add(&mut a_m_fence, &a_e_im) {
                add_to_container_shape(&a_e_im, &mut a_bounds);
                the_la_bounds.push(a_e_im.clone());
            }
            // OCCT L8638-8656.
            match shape_data_map::seek(&my_oe_origins, &a_e_im) {
                Some(p_lo) => {
                    let p_lo = p_lo.clone();
                    for a_eo in p_lo.iter() {
                        if set_add(&mut a_ma_valid, a_eo) {
                            the_la_valid.push(a_eo.clone());
                        }
                    }
                }
                None => {
                    if set_add(&mut a_ma_valid, &a_e_im) {
                        the_la_valid.push(a_e_im.clone());
                    }
                }
            }
        }
    }
    // OCCT L8658: theBounds = aBounds.
    *the_bounds = a_bounds;
}

// ---------------------------------------------------------------------------
// OCCT GetInvalidEdgesByBounds (cxx L8675-8996).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::GetInvalidEdgesByBounds (cxx
/// L8675-8996) — filter the new splits by the intersection with the
/// bounds.
#[allow(clippy::too_many_arguments)]
pub(crate) fn get_invalid_edges_by_bounds_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_splits: &Shape,
    the_bounds: &Shape,
    the_mv_old: &OcctShapeSet,
    the_me_new: &OcctShapeSet,
    the_dme_or: &ShapeDataMap<Vec<Shape>>,
    the_melf: &ShapeDataMap<Vec<Shape>>,
    the_e_images: &ShapeDataMap<Vec<Shape>>,
    the_mecheckext: &OcctShapeSet,
    the_me_inv_on_art: &OcctShapeSet,
    the_verts_to_avoid: &mut OcctShapeSet,
    the_me_inv: &mut OcctShapeSet,
) {
    // OCCT L8684-8686: map the splits to check the vertices of the edges.
    let mut a_dmve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    map_shapes_and_ancestors_map(the_splits, ShapeType::Vertex, ShapeType::Edge, &mut a_dmve);

    // OCCT L8688-8692: BOPAlgo_Section aSec; AddArgument(theSplits);
    // AddArgument(theBounds); Perform().
    let mut a_sec_args: Vec<Shape> = Vec::new();
    a_sec_args.push(the_splits.clone());
    a_sec_args.push(the_bounds.clone());
    let mut a_sec_filler = PaveFiller::new();
    a_sec_filler.set_arguments(a_sec_args);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_Section", 1);
    a_sec_filler.perform(&a_ps);
    let mut a_sec = BOPAlgoSection::new(a_sec_filler.ds(), a_sec_filler.fuzzy_value());
    a_sec.set_arguments(a_sec_filler.ds().arguments.clone());
    a_sec.perform();
    let a_sec_has_errors = a_sec.has_errors();

    // OCCT L8695-8698: the invalid vertices and the vertices to check
    // additionally.
    let mut a_mv_inv = OcctIndexedShapeMap::new();
    let mut a_mv_check_add: OcctShapeSet = HashMap::new();

    // OCCT L8700-8701: pDS = aSec.PDS() — the DS of the section builder.
    let p_ds = a_sec.builder.ds;

    // OCCT L8703-8753: check the edge/edge intersections.
    if !a_sec_has_errors {
        let a_ee_walk = |a_bf: &mut BRepOffsetBuildOffsetFaces,
                         the_me_inv: &mut OcctShapeSet,
                         a_mv_inv: &mut OcctIndexedShapeMap,
                         a_mv_check_add: &mut OcctShapeSet,
                         a_dmve: &ShapeIndexedDataMap<Vec<Shape>>| {
            let a_ees = p_ds.interf_ee.clone();
            let _a_nb = a_ees.len();
            for a_ee in a_ees.iter() {
                // OCCT L8708-8709: aE1/aE2 — the DS shapes.
                let a_e1 = p_ds.shape(a_ee.e1).clone();
                let a_e2 = p_ds.shape(a_ee.e2).clone();
                // OCCT L8711: HasIndexNew() — the rcad new_vertex sentinel
                // (usize::MAX when absent); the common part of the
                // point-interference carries the new vertex, the
                // edge-overlap interference has none (the OCCT
                // CommonPart().Type() probe of the !HasIndexNew branch is
                // the EDGE part by construction).
                let b_has_index_new = a_ee.new_vertex != usize::MAX;
                if !b_has_index_new {
                    // OCCT L8713-8717.
                    if set_contains(the_mecheckext, &a_e1) {
                        set_add(the_me_inv, &a_e1);
                    }
                    continue;
                }
                // OCCT L8720-8723.
                if a_bf.my_invalid_edges.contains(&a_e2) {
                    set_add(the_me_inv, &a_e1);
                }
                // OCCT L8725-8736.
                if set_contains(the_me_inv_on_art, &a_e2) {
                    // avoid checking the vertices of the split edge
                    // intersected by the invalid edge from the artificial
                    // face.
                    let (a_v1, a_v2) = super::brep_offset_tool::top_exp_vertices(&a_e2);
                    if a_dmve.contains_key(&shape_key_of(&a_v1))
                        && a_dmve.contains_key(&shape_key_of(&a_v2))
                    {
                        continue;
                    }
                }
                // OCCT L8738-8753: add the vertices of all the images of the
                // edge from the splits for checking.
                let a_le_or = shape_data_map::find(the_dme_or, &a_e1);
                for a_e_or in a_le_or.iter() {
                    let p_le_im = match shape_data_map::seek(the_e_images, a_e_or) {
                        Some(v) => v.clone(),
                        None => continue,
                    };
                    for a_e_im in p_le_im.iter() {
                        for a_v in bat::sub_shapes(a_e_im) {
                            if !set_contains(the_mv_old, &a_v) {
                                a_mv_inv.add(&a_v);
                                set_add(a_mv_check_add, &a_v);
                            }
                        }
                    }
                }
            }
        };
        a_ee_walk(a_bf, the_me_inv, &mut a_mv_inv, &mut a_mv_check_add, &a_dmve);
    }

    // OCCT L8755-8762: the section result (the common blocks are contained
    // in the result of the SECTION operation between the edge sets).
    let a_sec_r = a_sec
        .builder
        .my_shape
        .as_ref()
        .map(super::brep_offset_make_offset_1::brep_root_shape)
        .unwrap_or_else(Shape::null);
    let mut a_ms_sec = OcctIndexedShapeMap::new();
    // OCCT L8762: TopExp::MapShapes(aSecR, aMSSec) — the all-types form.
    for a_s in explorer(&a_sec_r, ShapeType::Shape, ShapeType::Shape) {
        a_ms_sec.add(&a_s);
    }

    // OCCT L8764-8810: the images walk — anIm = aSec.Images().
    for a_e in explorer(the_splits, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L8768-8772: aSec.IsDeleted(aE) — the GAP form takes the
        // OCCT false path (the rcad section keeps no deleted table).
        let b_deleted = false;
        if b_deleted {
            continue;
        }
        let p_le_im = match a_sec_images(&a_sec).get((a_e.ptr_id(), a_e.location)) {
            Some(v) => v.clone(),
            None => continue,
        };
        let mut b_broke = false;
        for a_e_im in p_le_im.iter() {
            if !a_ms_sec.contains(a_e_im) {
                // OCCT L8782-8802: the edge is included in the section only
                // partially — check its vertices; if one of them is new the
                // edge might be removed.
                let (a_v1, a_v2) = super::brep_offset_tool::top_exp_vertices(a_e_im);
                if !set_contains(the_mv_old, &a_v1) || !set_contains(the_mv_old, &a_v2) {
                    // OCCT L8790-8797: make the new vertex in the middle of
                    // the edge and add it for checking.
                    let a_v = middle_vertex(a_e_im);
                    {
                        let entry = a_dmve
                            .entry(shape_key_of(&a_v))
                            .or_insert((a_v.clone(), Vec::new()));
                        entry.1.push(a_e.clone());
                    }
                    a_mv_inv.add(&a_v);
                    b_broke = true;
                    break;
                }
            }
        }
        let _ = b_broke;
    }

    // OCCT L8812-8828: add for check the edges created from the common
    // between the splits of the offset faces edges.
    for a_e in the_mecheckext.values() {
        // OCCT L8818-8821: make the new vertex in the middle of the edge.
        let a_v = middle_vertex(a_e);
        {
            let entry = a_dmve
                .entry(shape_key_of(&a_v))
                .or_insert((a_v.clone(), Vec::new()));
            entry.1.push(a_e.clone());
        }
        a_mv_inv.add(&a_v);
    }

    // OCCT L8830-8868: add for check also the vertices connected only to
    // the new or old edges.
    {
        let a_nb = a_dmve.len();
        for i in 1..=a_nb {
            let (a_v, a_lex) = {
                let (_k, v) = a_dmve.get_index(i - 1).expect("index");
                (v.0.clone(), v.1.clone())
            };
            if set_contains(the_mv_old, &a_v) {
                continue;
            }
            let mut b_new = false;
            let mut b_old = false;
            for a_e in a_lex.iter() {
                if set_contains(the_mecheckext, a_e) {
                    continue;
                }
                if set_contains(the_me_new, a_e) {
                    b_new = true;
                } else {
                    b_old = true;
                }
                if b_new && b_old {
                    break;
                }
            }
            if !b_new || !b_old {
                a_mv_inv.add(&a_v);
                a_mv_check_add.remove(&shape_key_of(&a_v));
            }
        }
    }

    // OCCT L8870-8986: perform the checking of the vertices.
    let a_nb_iv = a_mv_inv.extent();
    for iv in 1..=a_nb_iv {
        let a_v = a_mv_inv.find_key_1(iv).clone();
        if set_contains(the_mv_old, &a_v) {
            continue;
        }
        let p_le_inv = match a_dmve.get(&shape_key_of(&a_v)) {
            Some(v) => v.1.clone(),
            None => continue,
        };
        // OCCT L8878-8890: find the faces by the edges to check the vertex.
        let mut a_mf = OcctIndexedShapeMap::new();
        for a_e in p_le_inv.iter() {
            let a_lf = match shape_data_map::seek(the_melf, a_e) {
                Some(v) => v.clone(),
                None => continue,
            };
            for a_f in a_lf.iter() {
                a_mf.add(a_f);
            }
        }
        // OCCT L8892-8920: check the vertex to belong to some split of the
        // faces.
        let mut b_invalid = true;
        let a_nb_f = a_mf.extent();
        for i in 1..=a_nb_f {
            if !b_invalid {
                break;
            }
            let a_f = a_mf.find_key_1(i).clone();
            let a_lf_im = match a_bf.my_of_images.get(&shape_key_of(&a_f)) {
                Some(v) => v.1.clone(),
                None => continue,
            };
            for a_f_im in a_lf_im.iter() {
                if !b_invalid {
                    break;
                }
                for a_vf in explorer(a_f_im, ShapeType::Vertex, ShapeType::Shape) {
                    if b_invalid {
                        b_invalid = !a_vf.is_same(&a_v);
                    }
                }
            }
            // OCCT L8910-8919: the classification.
            if b_invalid {
                let i_status = compute_vf_gap(&a_v, &a_f);
                if i_status == 0 {
                    // OCCT L8914-8918: classify the point relatively the
                    // faces (the GAP IsPointInOnFace keeps the OCCT
                    // structure).
                    let _a_p2d = glam::DVec2::ZERO;
                    for a_f_im in a_lf_im.iter() {
                        if !b_invalid {
                            break;
                        }
                        b_invalid = !is_point_in_on_face_gap(a_f_im, _a_p2d);
                    }
                }
            }
        }
        // OCCT L8922-8937: check the same vertex for the solids.
        if b_invalid && set_contains(&a_mv_check_add, &a_v) {
            let a_p = bat::brep_tool_pnt(&a_v).unwrap_or(glam::DVec3::ZERO);
            let a_tol_v = bat::brep_tool_tolerance(&a_v);
            for a_sol in explorer(&a_bf.my_solids, ShapeType::Solid, ShapeType::Shape) {
                if !b_invalid {
                    break;
                }
                // OCCT L8929-8933: myContext->SolidClassifier(aSol);
                // aSC.Perform(aP, aTolV); bInvalid = (aSC.State() == OUT).
                let mut a_sc = crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(&a_sol);
                a_sc.perform(a_p, a_tol_v);
                b_invalid = a_sc.state() == 1; // OCCT TopAbs_OUT
            }
        }
        // OCCT L8940-8947.
        if b_invalid {
            set_add(the_verts_to_avoid, &a_v);
            for a_e in p_le_inv.iter() {
                set_add(the_me_inv, a_e);
            }
        }
    }
}

/// OCCT aSec.Images() — the images table of the section builder.
fn a_sec_images<'a>(a_sec: &'a BOPAlgoSection<'a>) -> &'a crate::bop::algo::occt_map::OcctDataMapInt<ShapeKey, Vec<Shape>> {
    &a_sec.builder.my_images
}

/// OCCT `BRep_Builder().MakeVertex(aV, aC->Value((f + l) * 0.5),
/// Precision::Confusion())` (cxx L8791-8796 / L8818-8821) — the vertex in
/// the middle of the edge curve.
fn middle_vertex(the_e: &Shape) -> Shape {
    let a_v = bat::builder_make_vertex();
    let mid = match bat::brep_tool_curve(the_e) {
        Some((c, f, l)) => {
            let p = CurveEval::point_at(&c, 0.5 * (f + l));
            let _ = &c;
            p
        }
        None => {
            let _ = std::marker::PhantomData::<Curve3>;
            glam::DVec3::ZERO
        }
    };
    let mut a_v = a_v;
    bat::builder_update_vertex_point_tol(&mut a_v, mid, rcad_kernel::precision::CONFUSION);
    a_v
}

// ---------------------------------------------------------------------------
// OCCT FilterSplits (cxx L8997-9047).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FilterSplits (cxx L8997-9047) —
/// filter the images of the edges from the invalid edges.
pub(crate) fn filter_splits_impl(
    _a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_le: &[Shape],
    the_me_filter: &OcctShapeSet,
    the_is_inv: bool,
    the_e_images: &mut ShapeDataMap<Vec<Shape>>,
    the_splits: &mut Shape,
) {
    // OCCT L9002-9004.
    let mut a_splits = empty_compound();
    let mut a_m_fence: OcctShapeSet = HashMap::new();

    for a_e in the_le.iter() {
        let p_le_im = match super::brep_offset_make_offset_1::dm_change_seek(the_e_images, a_e) {
            Some(v) => v,
            None => continue,
        };
        let mut kept: Vec<Shape> = Vec::new();
        for a_e_im in p_le_im.iter() {
            // OCCT L9019-9024: the filter keeps/removes by theIsInv flag.
            if set_contains(the_me_filter, a_e_im) == the_is_inv {
                continue;
            }
            kept.push(a_e_im.clone());
            if set_add(&mut a_m_fence, a_e_im) {
                add_to_container_shape(a_e_im, &mut a_splits);
            }
        }
        // OCCT L9036-9040.
        if kept.is_empty() {
            shape_data_map::un_bind(the_e_images, a_e);
        } else {
            if let Some(entry) = the_e_images.get_mut(&shape_key_of(a_e)) {
                entry.1 = kept;
            }
        }
    }
    // OCCT L9042: theSplits = aSplits.
    *the_splits = a_splits;
}

// ---------------------------------------------------------------------------
// OCCT UpdateNewIntersectionEdges (cxx L9048-9296).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::UpdateNewIntersectionEdges (cxx
/// L9048-9296) — updating the maps of images and origins of the offset
/// edges.
pub(crate) fn update_new_intersection_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_le: &[Shape],
    the_melf: &ShapeDataMap<Vec<Shape>>,
    the_e_images: &ShapeDataMap<Vec<Shape>>,
    the_eetrim: &mut ShapeDataMap<Vec<Shape>>,
) {
    let mut my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();
    let mut my_oe_images = a_bf.my_oe_images.clone();
    let mut my_oe_origins = a_bf.my_oe_origins.clone();
    let mut my_modified_edges = a_bf.my_modified_edges.clone();
    let my_as_des = a_bf.my_as_des.clone();
    let mut my_e_trim_e_inf = a_bf.my_e_trim_e_inf.clone().unwrap_or_default();

    let _ = &my_as_des;

    // OCCT L9055: aLEImEmpty.
    let a_le_im_empty: Vec<Shape> = Vec::new();

    // OCCT L9057-9295.
    for a_e in the_le.iter() {
        // OCCT L9060-9080: when the edge has no images — filter the trimmed
        // parts.
        let has_images = shape_data_map::is_bound(the_e_images, a_e);
        if !has_images {
            let p_let = match super::brep_offset_make_offset_1::dm_change_seek(the_eetrim, a_e) {
                Some(v) => v,
                None => continue,
            };
            let mut kept: Vec<Shape> = Vec::new();
            for a_et in p_let.iter() {
                if !a_bf.my_invalid_edges.contains(a_et) && !a_bf.my_inverted_edges.contains(a_et)
                {
                    // OCCT L9069-9071: pLET->Remove(aItLET).
                    continue;
                }
                kept.push(a_et.clone());
            }
            *p_let = kept;
            if p_let.is_empty() {
                continue;
            }
        }

        // OCCT L9082-9083: the new images.
        let a_le_new: Vec<Shape> = if shape_data_map::is_bound(the_e_images, a_e) {
            shape_data_map::find(the_e_images, a_e)
        } else {
            a_le_im_empty.clone()
        };

        // OCCT L9085-9090: save the connection to the untrimmed edge.
        for a_et in a_le_new.iter() {
            shape_data_map::bind(&mut my_e_trim_e_inf, a_et, a_e.clone());
            set_add(&mut my_modified_edges, a_et);
        }

        // OCCT L9092-9138: check if it is an existing edge.
        if !shape_data_map::is_bound(the_eetrim, a_e) {
            let a_lf = shape_data_map::find(the_melf, a_e);
            // OCCT L9097-9102: add this edge to AsDes.
            if let Some(my_as_des) = &my_as_des {
                for a_f in a_lf.iter() {
                    my_as_des.borrow_mut().add(a_f, a_e);
                }
            }
            // OCCT L9104-9106: add aE to the images.
            shape_data_map::bind(&mut my_oe_images, a_e, a_le_new.clone());
            set_add(&mut my_modified_edges, a_e);
            // OCCT L9108-9123: add to the origins.
            for a_e_new in a_le_new.iter() {
                if shape_data_map::is_bound(&my_oe_origins, a_e_new) {
                    let a_e_origins = shape_data_map::change_find(&mut my_oe_origins, a_e_new);
                    append_to_list(a_e_origins, a_e);
                } else {
                    let a_e_origins = vec![a_e.clone()];
                    shape_data_map::bind(&mut my_oe_origins, a_e_new, a_e_origins);
                }
            }
            // OCCT L9125-9138: update the connection to the initial origins.
            if shape_data_map::is_bound(&my_edges_origins, a_e) {
                let a_le_or_init = shape_data_map::find(&my_edges_origins, a_e);
                for a_e_new in a_le_new.iter() {
                    if shape_data_map::is_bound(&my_edges_origins, a_e_new) {
                        let a_le_new_or =
                            shape_data_map::change_find(&mut my_edges_origins, a_e_new);
                        for a_e_or in a_le_or_init.iter() {
                            append_to_list(a_le_new_or, a_e_or);
                        }
                    } else {
                        shape_data_map::bind(
                            &mut my_edges_origins,
                            a_e_new,
                            a_le_or_init.clone(),
                        );
                    }
                }
            }
            continue;
        }

        // OCCT L9140-9240: the old images.
        let a_le_old = shape_data_map::find(the_eetrim, a_e);
        // OCCT L9143: the list of the initial origins.
        let mut an_init_origins: Vec<Shape> = Vec::new();

        for a_e_old in a_le_old.iter() {
            if shape_data_map::is_bound(&my_oe_origins, a_e_old) {
                // OCCT L9148-9151: get its origins.
                let a_e_origins = shape_data_map::find(&my_oe_origins, a_e_old);
                for a_e_or in a_e_origins.iter() {
                    set_add(&mut my_modified_edges, a_e_or);
                    // OCCT L9156-9157: aEImages = myOEImages.ChangeFind(aEOr).
                    let images_bound =
                        shape_data_map::is_bound(&my_oe_images, a_e_or);
                    let mut a_e_images_local: Vec<Shape> = if images_bound {
                        shape_data_map::find(&my_oe_images, a_e_or)
                    } else {
                        Vec::new()
                    };
                    // OCCT L9159-9171: remove the old edge from the images.
                    a_e_images_local.retain(|a_e_im| !a_e_im.is_same(a_e_old));
                    // OCCT L9173-9188: add the new images.
                    for a_e_new in a_le_new.iter() {
                        append_to_list(&mut a_e_images_local, a_e_new);
                        if shape_data_map::is_bound(&my_oe_origins, a_e_new) {
                            let a_e_new_origins =
                                shape_data_map::change_find(&mut my_oe_origins, a_e_new);
                            append_to_list(a_e_new_origins, a_e_or);
                        } else {
                            let a_e_new_origins = vec![a_e_or.clone()];
                            shape_data_map::bind(&mut my_oe_origins, a_e_new, a_e_new_origins);
                        }
                    }
                    if images_bound {
                        if let Some(entry) = my_oe_images.get_mut(&shape_key_of(a_e_or)) {
                            entry.1 = a_e_images_local;
                        }
                    } else {
                        shape_data_map::bind(&mut my_oe_images, a_e_or, a_e_images_local);
                    }
                }
            } else {
                // OCCT L9190-9210: add to the images.
                shape_data_map::bind(&mut my_oe_images, a_e_old, a_le_new.clone());
                set_add(&mut my_modified_edges, a_e_old);
                // OCCT L9196-9208: add to the origins.
                for a_e_new in a_le_new.iter() {
                    if shape_data_map::is_bound(&my_oe_origins, a_e_new) {
                        let a_e_origins =
                            shape_data_map::change_find(&mut my_oe_origins, a_e_new);
                        append_to_list(a_e_origins, a_e_old);
                    } else {
                        let a_e_origins = vec![a_e_old.clone()];
                        shape_data_map::bind(&mut my_oe_origins, a_e_new, a_e_origins);
                    }
                }
            }
            // OCCT L9212-9221: update the connection to the initial shape.
            if shape_data_map::is_bound(&my_edges_origins, a_e_old) {
                let a_le_or_init = shape_data_map::find(&my_edges_origins, a_e_old);
                for a_e_or_init in a_le_or_init.iter() {
                    append_to_list(&mut an_init_origins, a_e_or_init);
                }
            }
        }

        // OCCT L9223-9240.
        if !an_init_origins.is_empty() {
            for a_e_new in a_le_new.iter() {
                if shape_data_map::is_bound(&my_edges_origins, a_e_new) {
                    let a_le_new_or = shape_data_map::change_find(&mut my_edges_origins, a_e_new);
                    for a_e_or in an_init_origins.iter() {
                        append_to_list(a_le_new_or, a_e_or);
                    }
                } else {
                    shape_data_map::bind(&mut my_edges_origins, a_e_new, an_init_origins.clone());
                }
            }
        }
    }

    a_bf.my_oe_images = my_oe_images;
    a_bf.my_oe_origins = my_oe_origins;
    a_bf.my_modified_edges = my_modified_edges;
    a_bf.my_edges_origins = Some(my_edges_origins);
    a_bf.my_e_trim_e_inf = Some(my_e_trim_e_inf);
}

// ---------------------------------------------------------------------------
// OCCT FillGaps (cxx L9297-9412).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FillGaps (cxx L9297-9412) — fill the
/// possible gaps (holes) in the splits of the offset faces to increase the
/// possibility of creating a closed volume from these splits.
#[allow(unused_variables)]
pub(crate) fn fill_gaps_impl(a_bf: &mut BRepOffsetBuildOffsetFaces, _the_range: &ProgressScope) {
    // OCCT L9298-9301.
    let a_nb_f = a_bf.my_of_images.len();
    if a_nb_f == 0 {
        return;
    }

    let _a_ps = ProgressScope::new(&NoopProgress, "Filling gaps", 2 * a_nb_f);

    // OCCT L9306-9320: map the splits of the faces to find the free edges.
    let mut an_efmap: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    for i in 1..=a_nb_f {
        let a_lf_im: Vec<Shape> = a_bf
            .my_of_images
            .get_index(i - 1)
            .expect("myOFImages index")
            .1
             .1
            .clone();
        for a_f_im in a_lf_im.iter() {
            map_shapes_and_ancestors_map(
                a_f_im,
                ShapeType::Edge,
                ShapeType::Face,
                &mut an_efmap,
            );
        }
    }

    // OCCT L9322-9400: analyze the images of each offset face on the
    // presence of the free edges.
    for i in 1..=a_nb_f {
        let (a_f, mut a_lf_images) = {
            let (_k, v) = a_bf.my_of_images.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        if a_lf_images.is_empty() {
            continue;
        }

        // OCCT L9334-9336: collect all the edges from the splits.
        let mut an_edges = empty_compound();
        // OCCT L9338-9340: collect all the free edges into a map with the
        // reverted orientation.
        let mut a_free_edges_map: OcctShapeSet = HashMap::new();
        for a_f_im in a_lf_images.iter() {
            for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                // OCCT L9345-9349: skip the internals.
                if a_e.orientation != Orientation::Forward
                    && a_e.orientation != Orientation::Reversed
                {
                    continue;
                }
                // OCCT L9351-9355: the free edges.
                let a_lf = match an_efmap.get(&shape_key_of(&a_e)) {
                    Some(v) => v.1.clone(),
                    None => continue,
                };
                if a_lf.len() == 1 {
                    let a_e_rev = bat::reversed(&a_e);
                    set_add(&mut a_free_edges_map, &a_e_rev);
                }
                // OCCT L9357: collect the edge.
                add_to_container_shape(&a_e, &mut an_edges);
            }
        }

        // OCCT L9360-9364: no free edges.
        if a_free_edges_map.is_empty() {
            continue;
        }

        // OCCT L9366-9376: free edges are found — build the new splits
        // using all the kept edges.
        let mut a_lf_new: Vec<Shape> = Vec::new();
        let mut a_dummy: ShapeDataMap<Shape> = HashMap::new();
        build_splits_of_face(&a_f, &an_edges, &mut a_dummy, &mut a_lf_new);

        // OCCT L9378-9398: find the faces filling the holes.
        for a_f_new in a_lf_new.iter() {
            let mut b_found = false;
            for a_e in explorer(a_f_new, ShapeType::Edge, ShapeType::Shape) {
                if set_contains(&a_free_edges_map, &a_e) {
                    b_found = true;
                    break;
                }
            }
            if b_found {
                // OCCT L9384-9387: add the face to the splits.
                a_lf_images.push(a_f_new.clone());
            }
        }
        super::brep_offset_make_offset_1::add_of_images(&mut a_bf.my_of_images, &a_f, a_lf_images);
    }
}

// ---------------------------------------------------------------------------
// OCCT FillHistory (cxx L9413-9491).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FillHistory (cxx L9413-9491) — saving
/// the obtained results in the history tools.
pub(crate) fn fill_history_impl(a_bf: &mut BRepOffsetBuildOffsetFaces, the_image: &mut crate::brep_algo::image::BRepAlgoImage) {
    // OCCT L9414-9417.
    let a_nb_f = a_bf.my_of_images.len();
    if a_nb_f == 0 {
        return;
    }

    // OCCT L9424: the map of the kept edges.
    let mut an_edges_map = OcctIndexedShapeMap::new();

    // OCCT L9426-9449: fill the history for the faces.
    for i in 1..=a_nb_f {
        let (a_f, a_lf_images) = {
            let (_k, v) = a_bf.my_of_images.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        if a_lf_images.is_empty() {
            continue;
        }
        // OCCT L9440-9448: Add or Bind the splits to the history map.
        if the_image.has_image(&a_f) {
            the_image.add_list(&a_f, &a_lf_images);
        } else {
            the_image.bind_list(&a_f, &a_lf_images);
        }

        // OCCT L9450-9460: collect the edges from the splits.
        for a_f_im in a_lf_images.iter() {
            map_shapes_indexed(a_f_im, ShapeType::Edge, &mut an_edges_map);
        }
    }

    // OCCT L9463-9490: fill the history for the edges.
    let a_keys: Vec<ShapeKey> = a_bf.my_oe_images.keys().cloned().collect();
    for a_key in a_keys.iter() {
        let (a_e, a_le_im) = match a_bf.my_oe_images.get(a_key) {
            Some(v) => (v.0.clone(), v.1.clone()),
            None => continue,
        };
        let mut b_has_image = the_image.has_image(&a_e);
        for a_e_im in a_le_im.iter() {
            if an_edges_map.contains(a_e_im) {
                if b_has_image {
                    the_image.add(&a_e, a_e_im);
                } else {
                    the_image.bind(&a_e, a_e_im);
                    b_has_image = true;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepOffset_MakeOffset dispatchers (cxx L9493-9533).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_MakeOffset::BuildSplitsOfTrimmedFaces (cxx L9497-9506)
/// — building the splits of the already trimmed faces (the file carries
/// the method; the BRepOffset_MakeOffset class translation calls this
/// form).
#[allow(unused)]
pub(crate) fn brep_offset_make_offset_build_splits_of_trimmed_faces(
    the_lf: &[Shape],
    the_as_des: &std::rc::Rc<std::cell::RefCell<crate::brep_algo::as_des::BRepAlgoAsDes>>,
    the_image: &mut crate::brep_algo::image::BRepAlgoImage,
    the_range: &ProgressScope,
) {
    let mut a_bf_tool = BRepOffsetBuildOffsetFaces::new();
    a_bf_tool.set_faces(the_lf);
    a_bf_tool.set_as_des_info(the_as_des);
    a_bf_tool.build_splits_of_trimmed_faces(the_image, the_range);
}

/// OCCT BRepOffset_MakeOffset::BuildSplitsOfExtendedFaces (cxx
/// L9514-9533) — building the splits of the not-trimmed offset faces; for
/// the cases in which the invalidity will be found, these invalidities
/// will be rebuilt.
#[allow(clippy::too_many_arguments)]
#[allow(unused)]
pub(crate) fn brep_offset_make_offset_build_splits_of_extended_faces(
    the_lf: &[Shape],
    the_analyse: &super::brep_offset_tool_d::BRepOffsetAnalyse,
    the_as_des: &std::rc::Rc<std::cell::RefCell<crate::brep_algo::as_des::BRepAlgoAsDes>>,
    the_edges_origins: &ShapeDataMap<Vec<Shape>>,
    the_faces_origins: &ShapeDataMap<Shape>,
    the_e_trim_e_inf: &ShapeDataMap<Shape>,
    the_image: &mut crate::brep_algo::image::BRepAlgoImage,
    the_range: &ProgressScope,
) {
    let mut a_bf_tool = BRepOffsetBuildOffsetFaces::new();
    a_bf_tool.set_faces(the_lf);
    a_bf_tool.set_as_des_info(the_as_des);
    a_bf_tool.set_analysis(the_analyse);
    a_bf_tool.set_edges_origins(the_edges_origins);
    a_bf_tool.set_faces_origins(the_faces_origins);
    a_bf_tool.set_inf_edges(the_e_trim_e_inf);
    a_bf_tool.build_splits_of_extended_faces(the_image, the_range);
}

// The State import anchor of the classifier comparisons.
#[allow(unused_imports)]
use State as _state_anchor;
