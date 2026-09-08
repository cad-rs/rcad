// OCCT BRepOffset_MakeOffset_1.cxx L3089-3887 — module d of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L3089-3887:
//   - FindFacesInsideHoleWires (cxx L3089-3310)
//   - CheckInverted (cxx L3316-3514)
//   - GetVerticesOnEdges (cxx L3522-3542, file static)
//   - CheckInvertedBlock (cxx L3549-3680)
//   - RemoveInvalidSplitsByInvertedEdges (cxx L3687-3887)
//
// Local dependency carriers of this module:
//   - BOPTools_AlgoTools3D::PointInFace (TKBO) — the GAP carrier with the
//     OCCT error-path structure (iErr != 0 skips the face; the TKBO
//     translation unit has no rcad form yet).
//   - IntTools_Context::IsValidPointForFace — the local reduced re-host
//     over the FClass2dCarrier (the brep_offset_tool_d.rs precedent; the
//     rcad IntToolsContext form is DS-bound and the offset pipeline works
//     on the BRep shapes).
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#48).

use std::collections::HashMap;

use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::pave_filler::PaveFiller;
use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    append_to_list, block_shape, empty_compound, find_shape, map_shapes_and_ancestors_map,
    map_shapes_indexed, shape_key_of, BRepOffsetBuildOffsetFaces,
};
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, shape_indexed_data_map, OcctIndexedShapeMap,
    OcctShapeSet, ShapeDataMap, ShapeIndexedDataMap,
};
use crate::feat::loc_ope_wires_on_shape_b::ShapeKey;
use super::brep_offset_make_offset_1_c::outer_wire_of;

// ---------------------------------------------------------------------------
// Local dependency carriers.
// ---------------------------------------------------------------------------

/// OCCT BOPTools_AlgoTools3D::PointInFace(F, P3D, P2D, context)
/// (BOPTools_AlgoTools3D.cxx) — GAP carrier: the TKBO unit has no rcad
/// translation yet; the carrier takes the OCCT error path (iErr != 0 —
/// the cxx L3252-3256 `if (iErr) continue;` structure is preserved).
pub(crate) fn point_in_face_gap(_the_f: &Shape) -> i32 {
    // GAP: BOPTools_AlgoTools3D::PointInFace (BOPTools_AlgoTools3D not
    // translated) — the OCCT error path is taken.
    1
}

/// OCCT IntTools_Context::IsValidPointForFace(P, F, Tol)
/// (IntTools_Context.cxx) — the reduced re-host: the point must project
/// onto the face surface within Tol and classify IN (the FClass2dCarrier
/// form of brep_offset_tool_d.rs; the rcad IntToolsContext is DS-bound).
pub(crate) fn is_valid_point_for_face(a_p: glam::DVec3, the_f: &Shape, a_tol: f64) -> bool {
    use rcad_kernel::geom::{Surface3, SurfaceEval};
    let surf = match the_f.as_face().and_then(|fd| fd.surface.clone()) {
        Some(s) => s,
        None => return false,
    };
    // The projection of the point onto the surface (the
    // IntTools_Context::ProjPnt reduced form: the plane dot-product
    // parameterization; other surface kinds take the domain center).
    let uv = match &surf {
        Surface3::Plane(pl) => {
            let d = a_p - pl.origin;
            glam::DVec2::new(d.dot(pl.u_dir), d.dot(pl.v_dir))
        }
        _ => {
            let dom = surf.default_domain();
            glam::DVec2::new(0.5 * (dom[0] + dom[1]), 0.5 * (dom[2] + dom[3]))
        }
    };
    let p_surf = surf.point_at(uv.x, uv.y);
    if p_surf.distance(a_p) > a_tol {
        return false;
    }
    // The classification (IntTools_Context::FClass2d Perform).
    let fc = super::brep_offset_tool_d::FClass2dCarrier::new(the_f);
    matches!(
        fc.perform(uv),
        rcad_kernel::topo::topods::State::In | rcad_kernel::topo::topods::State::On
    )
}

// ---------------------------------------------------------------------------
// OCCT FindFacesInsideHoleWires (cxx L3089-3310).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FindFacesInsideHoleWires (cxx
/// L3089-3310) — find faces inside the hole wires of the original face.
pub(crate) fn find_faces_inside_hole_wires_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_f_origin: &Shape,
    the_f_offset: &Shape,
    the_lf_images: &[Shape],
    the_dmeorleim: &ShapeDataMap<Vec<Shape>>,
    the_ef_map: &ShapeIndexedDataMap<Vec<Shape>>,
    the_mf_holes: &mut OcctShapeSet,
) {
    // OCCT L3100-3103.
    if the_lf_images.is_empty() {
        return;
    }
    // OCCT L3105-3116: find all hole wires in the original face.
    let mut a_l_hole_wires: Vec<Shape> = Vec::new();
    let an_outer_wire = outer_wire_of(the_f_origin);
    for a_wire_cur in explorer(the_f_origin, ShapeType::Wire, ShapeType::Shape) {
        let a_hole_wire = a_wire_cur;
        if !a_hole_wire.is_same(&an_outer_wire) && a_hole_wire.orientation != Orientation::Internal
        {
            a_l_hole_wires.push(a_hole_wire);
        }
    }
    // OCCT L3118-3122: no holes in the face.
    if a_l_hole_wires.is_empty() {
        return;
    }

    // OCCT L3124-3129: pLFNewHoles = myFNewHoles.ChangeSeek/Bound(theFOrigin).
    let a_key_origin = shape_key_of(the_f_origin);
    if !a_bf.my_f_new_holes.contains_key(&a_key_origin) {
        a_bf
            .my_f_new_holes
            .insert(a_key_origin, (the_f_origin.clone(), Vec::new()));
    }

    // OCCT L3130-3233: when the new hole faces are not built yet — build
    // them.
    if a_bf
        .my_f_new_holes
        .get(&a_key_origin)
        .map(|v| v.1.is_empty())
        .unwrap_or(true)
    {
        // OCCT L3138-3143: map the edges of the splits.
        let mut a_me_splits = OcctIndexedShapeMap::new();
        for a_f_im in the_lf_images.iter() {
            map_shapes_indexed(a_f_im, ShapeType::Edge, &mut a_me_splits);
        }

        for a_hole_wire in a_l_hole_wires.iter() {
            // OCCT L3149-3169: find the images of all edges of the original
            // wire.
            let mut a_me_im_wire = OcctIndexedShapeMap::new();
            for a_e_or in bat::sub_shapes(a_hole_wire) {
                let p_le_im = match shape_data_map::seek(the_dmeorleim, &a_e_or) {
                    Some(v) if !v.is_empty() => v.clone(),
                    _ => continue,
                };
                for a_e_im in p_le_im.iter() {
                    if a_me_splits.contains(a_e_im) {
                        a_me_im_wire.add(a_e_im);
                    }
                }
            }
            // OCCT L3171-3174.
            if a_me_im_wire.is_empty() {
                continue;
            }
            // OCCT L3176-3183: build a new planar face using these edges —
            // both orientations of each image edge.
            let mut a_le: Vec<Shape> = Vec::new();
            let a_nb_e = a_me_im_wire.extent();
            for i in 1..=a_nb_e {
                let a_e = a_me_im_wire.find_key_1(i).clone();
                a_le.push(bat::oriented(&a_e, Orientation::Forward));
                a_le.push(bat::oriented(&a_e, Orientation::Reversed));
            }
            // OCCT L3185-3194: BOPAlgo_BuilderFace over the FORWARD
            // oriented offset face (architecture difference #42 — the
            // DS-bound rcad BuilderFace driver).
            let mut a_ff = the_f_offset.clone();
            a_ff = bat::oriented(&a_ff, Orientation::Forward);
            let mut a_bf_args: Vec<Shape> = Vec::new();
            a_bf_args.push(a_ff.clone());
            a_bf_args.extend(a_le.iter().cloned());
            let mut a_bf_filler = PaveFiller::new();
            a_bf_filler.set_arguments(a_bf_args);
            let a_prog = rcad_kernel::core::message::NoopProgress;
            let a_ps = rcad_kernel::core::message::ProgressScope::new(
                &a_prog,
                "BOPAlgo_BuilderFace",
                1,
            );
            a_bf_filler.perform(&a_ps);
            let mut a_bface = crate::bop::algo::builder_face::BuilderFace::new(a_bf_filler.ds());
            a_bface.my_face = Some(a_bf_filler.ds().arguments[0].clone());
            a_bface.my_edges = a_bf_filler.ds().arguments[1..].to_vec();
            a_bface.perform();

            let a_lf_new: Vec<Shape> = a_bface.my_areas.clone();
            if a_lf_new.is_empty() {
                continue;
            }

            // OCCT L3196-3207: check that the outer edges in the new faces
            // are not inverted (an inverted edge means the hole has been
            // filled during the offset).
            let mut a_dmef_new: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
            for a_f_new in a_lf_new.iter() {
                map_shapes_and_ancestors_map(
                    a_f_new,
                    ShapeType::Edge,
                    ShapeType::Face,
                    &mut a_dmef_new,
                );
            }
            // OCCT L3209-3220.
            let a_nb_e = a_dmef_new.len();
            let mut b_broke = false;
            for i in 1..=a_nb_e {
                if shape_indexed_data_map::value_1(&a_dmef_new, i).len() == 1 {
                    let a_e = shape_indexed_data_map::find_key_1(&a_dmef_new, i).clone();
                    if a_bf.my_inverted_edges.contains(&a_e) {
                        b_broke = true;
                        break;
                    }
                }
            }
            // OCCT L3222-3225: if (i <= aNbE) continue.
            if b_broke {
                continue;
            }
            // OCCT L3227-3231: pLFNewHoles->Append(aItLFNew.Value()).
            let entry = a_bf
                .my_f_new_holes
                .get_mut(&a_key_origin)
                .expect("myFNewHoles entry");
            for a_f_new in a_lf_new.iter() {
                entry.1.push(a_f_new.clone());
            }
        }
    }

    // OCCT L3235-3240: the Edge-Face maps for the splits and for the holes.
    let mut an_ef_splits_map: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let mut an_ef_holes_map: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();

    // OCCT L3242-3272: among the splits of the offset face find those
    // located inside the hole faces.
    let p_lf_new_holes: Vec<Shape> = a_bf
        .my_f_new_holes
        .get(&a_key_origin)
        .map(|v| v.1.clone())
        .unwrap_or_default();
    for a_f_cur in the_lf_images.iter() {
        let a_f_im = a_f_cur.clone();
        map_shapes_and_ancestors_map(
            &a_f_im,
            ShapeType::Edge,
            ShapeType::Face,
            &mut an_ef_splits_map,
        );
        // OCCT L3250-3256: get a point inside the face (the
        // BOPTools_AlgoTools3D::PointInFace GAP carrier — the OCCT error
        // path skips the face).
        let i_err = point_in_face_gap(&a_f_im);
        if i_err != 0 {
            continue;
        }
        let a_tol = bat::brep_tool_tolerance(&a_f_im);
        for a_f_new in p_lf_new_holes.iter() {
            // OCCT L3264: myContext->IsValidPointForFace(aP3D, aFNew,
            // aTol) — the local reduced re-host.
            if is_valid_point_for_face(glam::DVec3::ZERO, a_f_new, a_tol) {
                // OCCT L3266-3269: the face is classified as IN.
                set_add(the_mf_holes, &a_f_im);
                map_shapes_and_ancestors_map(
                    &a_f_im,
                    ShapeType::Edge,
                    ShapeType::Face,
                    &mut an_ef_holes_map,
                );
                break;
            }
        }
    }

    // OCCT L3274-3309: out of all found holes find those which cannot be
    // removed by checking their connectivity to the splits of other offset
    // faces.
    let a_nb_e = an_ef_holes_map.len();
    for i in 1..=a_nb_e {
        let (an_edge, a_lf_holes) = {
            let (_k, v) = an_ef_holes_map.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        // OCCT L3283-3286: check if the edge is outer for the holes.
        if a_lf_holes.len() != 1 {
            continue;
        }
        let a_f_hole = a_lf_holes[0].clone();
        if !set_contains(the_mf_holes, &a_f_hole) {
            // OCCT L3291: already removed.
            continue;
        }
        // OCCT L3296-3300: check if the edge is not outer for the splits.
        let a_l_splits = shape_indexed_data_map::find(&an_ef_splits_map, &an_edge);
        if a_l_splits.len() == 1 {
            continue;
        }
        // OCCT L3302-3308: check if the edge is connected only to the
        // splits of the current offset face.
        let a_lf_all = shape_indexed_data_map::find(the_ef_map, &an_edge);
        if a_lf_all.len() == 2 {
            // Avoid removal of the hole from the splits.
            the_mf_holes.remove(&shape_key_of(&a_f_hole));
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT CheckInverted (cxx L3316-3514).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::CheckInverted (cxx L3316-3514) —
/// checks if the edge has been inverted: the direction from the first
/// vertex to the last vertex on the original edge is compared with the
/// same direction on the new edge.
#[allow(unused_assignments)]
pub(crate) fn check_inverted_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_e_im: &Shape,
    the_f_or: &Shape,
    the_dmve: &ShapeIndexedDataMap<Vec<Shape>>,
    the_m_edges: &OcctIndexedShapeMap,
) -> bool {
    let my_oe_origins = a_bf.my_oe_origins.clone();
    let my_oe_images = a_bf.my_oe_images.clone();
    let my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();

    // OCCT L3329-3334: vertices on the offset edge.
    let (mut a_vi1, mut a_vi2) = super::brep_offset_tool::top_exp_vertices(the_e_im);
    // OCCT L3330: vertices on the original edge.
    let mut a_vo1 = Shape::null();
    let mut a_vo2 = Shape::null();

    // OCCT L3336-3410: find the images.
    let mut a_le_images: Vec<Shape> = Vec::new();
    if shape_data_map::is_bound(&my_oe_origins, the_e_im) {
        // OCCT L3340-3341: anImages wire.
        let mut an_images = bat::builder_make_wire();
        let mut a_m_im_fence: OcctShapeSet = HashMap::new();
        let a_l_offset_or = shape_data_map::find(&my_oe_origins, the_e_im);
        for a_e_offset_or in a_l_offset_or.iter() {
            let a_l_images = shape_data_map::value(&my_oe_images, a_e_offset_or).clone();
            for an_im in a_l_images.iter() {
                if the_m_edges.contains(an_im) && set_add(&mut a_m_im_fence, an_im) {
                    bat::builder_add_wire_edge(&mut an_images, an_im);
                    a_le_images.push(an_im.clone());
                }
            }
        }
        // OCCT L3363-3380: find the alone vertices.
        let (mut a_vw1, mut a_vw2) = (Shape::null(), Shape::null());
        let mut a_dm_im_ve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
        map_shapes_and_ancestors_map(
            &an_images,
            ShapeType::Vertex,
            ShapeType::Edge,
            &mut a_dm_im_ve,
        );
        let mut a_lv_alone: Vec<Shape> = Vec::new();
        let a_nb = a_dm_im_ve.len();
        for i in 1..=a_nb {
            let a_l_im_e = shape_indexed_data_map::value_1(&a_dm_im_ve, i);
            if a_l_im_e.len() == 1 {
                a_lv_alone.push(shape_indexed_data_map::find_key_1(&a_dm_im_ve, i).clone());
            }
        }
        // OCCT L3382-3405.
        if a_lv_alone.len() > 1 {
            a_vw1 = a_lv_alone[0].clone();
            a_vw2 = a_lv_alone[a_lv_alone.len() - 1].clone();
            // OCCT L3387-3393: check distances.
            let a_pi1 = bat::brep_tool_pnt(&a_vi1).unwrap_or(glam::DVec3::ZERO);
            let a_pw1 = bat::brep_tool_pnt(&a_vw1).unwrap_or(glam::DVec3::ZERO);
            let a_pw2 = bat::brep_tool_pnt(&a_vw2).unwrap_or(glam::DVec3::ZERO);
            let a_dist1 = a_pi1.distance_squared(a_pw1);
            let a_dist2 = a_pi1.distance_squared(a_pw2);
            if a_dist1 < a_dist2 {
                a_vi1 = a_vw1;
                a_vi2 = a_vw2;
            } else {
                a_vi1 = a_vw2;
                a_vi2 = a_vw1;
            }
        }
    } else {
        // OCCT L3409: aLEImages.Append(theEIm).
        a_le_images.push(the_e_im.clone());
    }

    // OCCT L3412-3414: find the edges connected to these vertices.
    let a_lie1 = shape_indexed_data_map::find(the_dmve, &a_vi1);
    let a_lie2 = shape_indexed_data_map::find(the_dmve, &a_vi2);

    // OCCT L3416-3450: find the vertices on the original face corresponding
    // to the vertices on the offset edge — the original edges for both
    // lists.
    let mut a_loe1: Vec<Shape> = Vec::new();
    let mut a_loe2: Vec<Shape> = Vec::new();
    for i in 0..2 {
        let a_lie = if i == 0 { &a_lie1 } else { &a_lie2 };
        let mut a_m_fence: OcctShapeSet = HashMap::new();
        for a_ei in a_lie.iter() {
            if shape_data_map::is_bound(&my_edges_origins, a_ei) {
                let a_le_origins = shape_data_map::find(&my_edges_origins, a_ei);
                for a_eo in a_le_origins.iter() {
                    if a_eo.shape_type() == ShapeType::Edge && set_add(&mut a_m_fence, a_eo) {
                        let mut a_eo_in = Shape::null();
                        if find_shape(a_eo, the_f_or, None, &mut a_eo_in) {
                            if i == 0 {
                                append_to_list(&mut a_loe1, a_eo);
                            } else {
                                append_to_list(&mut a_loe2, a_eo);
                            }
                        }
                    }
                }
            }
        }
    }

    // OCCT L3452-3455.
    if a_loe1.len() < 2 || a_loe2.len() < 2 {
        return false;
    }

    // OCCT L3457-3486: find the vertices common for the max number of edges
    // in the lists.
    for i in 0..2 {
        let a_loe = if i == 0 { &a_loe1 } else { &a_loe2 };
        let mut a_vo = Shape::null();

        let mut a_dmve_loc: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
        for a_e in a_loe.iter() {
            map_shapes_and_ancestors_map(a_e, ShapeType::Vertex, ShapeType::Edge, &mut a_dmve_loc);
        }
        let mut a_nb_e_max = 0usize;
        for j in 1..=a_dmve_loc.len() {
            let a_nb_e = shape_indexed_data_map::value_1(&a_dmve_loc, j).len();
            if a_nb_e > 1 && a_nb_e > a_nb_e_max {
                a_vo = shape_indexed_data_map::find_key_1(&a_dmve_loc, j).clone();
                a_nb_e_max = a_nb_e;
            }
        }
        if a_vo.is_null() {
            return false;
        }
        if i == 0 {
            a_vo1 = a_vo;
        } else {
            a_vo2 = a_vo;
        }
    }

    // OCCT L3488-3491.
    if a_vo1.is_same(&a_vo2) {
        return false;
    }

    // OCCT L3493-3503: check the positions of the offset and original
    // vertices.
    let a_pi1 = bat::brep_tool_pnt(&a_vi1).unwrap_or(glam::DVec3::ZERO);
    let a_pi2 = bat::brep_tool_pnt(&a_vi2).unwrap_or(glam::DVec3::ZERO);
    let a_po1 = bat::brep_tool_pnt(&a_vo1).unwrap_or(glam::DVec3::ZERO);
    let a_po2 = bat::brep_tool_pnt(&a_vo2).unwrap_or(glam::DVec3::ZERO);
    let a_vi = a_pi2 - a_pi1;
    let a_vo = a_po2 - a_po1;
    let an_angle = if a_vi.length() > 0. && a_vo.length() > 0. {
        a_vi.angle_between(a_vo)
    } else {
        0.
    };
    let b_inverted = (an_angle - std::f64::consts::PI).abs() < 1.0e-4;
    if b_inverted {
        // OCCT L3506-3512.
        for a_e_invr in a_le_images.iter() {
            a_bf.my_inverted_edges.add(a_e_invr);
        }
    }
    b_inverted
}

// ---------------------------------------------------------------------------
// OCCT GetVerticesOnEdges (cxx L3522-3542, file static).
// ---------------------------------------------------------------------------

/// OCCT static GetVerticesOnEdges (cxx L3522-3542) — get the vertices from
/// the given shape belonging to the given edges.
pub(crate) fn get_vertices_on_edges(
    the_cb: &Shape,
    the_edges: &OcctIndexedShapeMap,
    the_vertices_on_edges: &mut OcctShapeSet,
    the_all_vertices: &mut OcctShapeSet,
) {
    for a_e in explorer(the_cb, ShapeType::Edge, ShapeType::Shape) {
        let is_on_given_edges = the_edges.contains(&a_e);
        for a_v in bat::sub_shapes(&a_e) {
            set_add(the_all_vertices, &a_v);
            if is_on_given_edges {
                set_add(the_vertices_on_edges, &a_v);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT CheckInvertedBlock (cxx L3549-3680).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::CheckInvertedBlock (cxx L3549-3680) —
/// checks if it is possible to remove the block containing the inverted
/// edges.
pub(crate) fn check_inverted_block_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_cb: &Shape,
    the_lcbf: &[Shape],
    the_dmcbv_inverted: &mut ShapeDataMap<OcctShapeSet>,
    the_dmcbv_all: &mut ShapeDataMap<OcctShapeSet>,
) -> bool {
    let my_oe_origins = a_bf.my_oe_origins.clone();

    // OCCT L3559-3564: 1. there should be more than just one face in the
    // block (NbChildren < 2).
    if bat::sub_shapes(the_cb).len() < 2 {
        return false;
    }

    // OCCT L3566-3588: 2. the block should contain at least two connected
    // inverted edges with different origins.
    let mut a_mecb_inv: OcctShapeSet = HashMap::new();
    let mut a_cecb_inv = empty_compound();
    for a_e in explorer(the_cb, ShapeType::Edge, ShapeType::Shape) {
        if a_bf.my_inverted_edges.contains(&a_e) && set_add(&mut a_mecb_inv, &a_e) {
            bat::builder_add_compound_shape(&mut a_cecb_inv, &a_e);
        }
    }
    if a_mecb_inv.len() < 2 {
        return false;
    }

    // OCCT L3590-3632: check that the edges are connected and different.
    let a_cecb_inv_edges: Vec<Shape> =
        explorer(&a_cecb_inv, ShapeType::Edge, ShapeType::Shape);
    let a_locations = [glam::DAffine3::IDENTITY];
    let a_lcbe = crate::bop::algo::wire_splitter::make_connexity_blocks(
        &a_cecb_inv_edges,
        &a_locations,
    );
    let mut b_found_block = false;
    for a_block in a_lcbe.iter() {
        // OCCT L3597: aCBE — the block (the connexity-block list carrier).
        let a_cbe = block_shape(&a_block.shapes);
        // OCCT L3598-3621: count the unique edges in the block.
        let mut a_nb_unique = 0usize;
        let mut a_me_origins: OcctShapeSet = HashMap::new();
        for a_e in bat::sub_shapes(&a_cbe) {
            let p_le_or = shape_data_map::seek(&my_oe_origins, &a_e);
            match p_le_or {
                None => {
                    set_add(&mut a_me_origins, &a_e);
                    a_nb_unique += 1;
                }
                Some(le_or) => {
                    let le_or = le_or.clone();
                    for a_e_or in le_or.iter() {
                        if set_add(&mut a_me_origins, a_e_or) {
                            a_nb_unique += 1;
                        }
                    }
                }
            }
        }
        // OCCT L3623-3632: if (aNbUnique >= 2) break; the iterator probe.
        if a_nb_unique >= 2 {
            b_found_block = true;
            break;
        }
    }
    // OCCT L3629-3632: if (!aItLCBE.More()) return false.
    if !b_found_block {
        return false;
    }

    // OCCT L3634-3649: 3. the block should not contain inverted edges whose
    // vertices are contained in the other blocks.
    let has_mv_inverted = super::brep_offset_make_offset_1::dm_change_seek(the_dmcbv_inverted, the_cb).is_some();
    if !has_mv_inverted {
        let entry_inv = the_dmcbv_inverted
            .entry(shape_key_of(the_cb))
            .or_insert((the_cb.clone(), HashMap::new()));
        let entry_all = the_dmcbv_all
            .entry(shape_key_of(the_cb))
            .or_insert((the_cb.clone(), HashMap::new()));
        let (inv, all) = (&mut entry_inv.1, &mut entry_all.1);
        get_vertices_on_edges(the_cb, &a_bf.my_inverted_edges, inv, all);
    }

    for a_cb1 in the_lcbf.iter() {
        if a_cb1.is_same(the_cb) {
            continue;
        }
        // OCCT L3660-3671: collect the vertices from the inverted edges.
        let has_mv_inverted1 =
            super::brep_offset_make_offset_1::dm_change_seek(the_dmcbv_inverted, a_cb1).is_some();
        if !has_mv_inverted1 {
            let entry_inv = the_dmcbv_inverted
                .entry(shape_key_of(a_cb1))
                .or_insert((a_cb1.clone(), HashMap::new()));
            let entry_all = the_dmcbv_all
                .entry(shape_key_of(a_cb1))
                .or_insert((a_cb1.clone(), HashMap::new()));
            let (inv, all) = (&mut entry_inv.1, &mut entry_all.1);
            get_vertices_on_edges(a_cb1, &a_bf.my_inverted_edges, inv, all);
        }
        // OCCT L3673-3676: NCollection_MapAlgo::HasIntersection.
        let p_mv_inverted = shape_data_map::value(the_dmcbv_inverted, the_cb);
        let p_mv_all1 = shape_data_map::value(the_dmcbv_all, a_cb1);
        if p_mv_inverted.keys().any(|k| p_mv_all1.contains_key(k)) {
            return false;
        }
    }

    // OCCT L3679.
    true
}

// ---------------------------------------------------------------------------
// OCCT RemoveInvalidSplitsByInvertedEdges (cxx L3687-3887).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::RemoveInvalidSplitsByInvertedEdges
/// (cxx L3687-3887) — looking for the invalid faces containing the
/// inverted edges that can be safely removed.
pub(crate) fn remove_invalid_splits_by_inverted_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_me_removed: &mut OcctIndexedShapeMap,
) {
    // OCCT L3690-3693.
    if a_bf.my_inverted_edges.is_empty() {
        return;
    }
    let my_oe_origins = a_bf.my_oe_origins.clone();

    // OCCT L3699-3701.
    let mut a_me_avoid = OcctIndexedShapeMap::new();
    let mut a_dmvf: ShapeDataMap<Vec<Shape>> = HashMap::new();

    // OCCT L3702-3779: check the faces on regularity — the splits of the
    // same face should not be connected only by a vertex.
    let a_nb = a_bf.my_of_images.len();
    for i in 1..=a_nb {
        let a_lf_im: Vec<Shape> = a_bf
            .my_of_images
            .get_index(i - 1)
            .expect("myOFImages index")
            .1
             .1
            .clone();

        // OCCT L3707-3708: aCFIm compound.
        let mut a_cf_im = empty_compound();

        // OCCT L3710-3711: aDMEF — the map to use only the outer edges.
        let mut a_dmef: ShapeDataMap<Vec<Shape>> = HashMap::new();
        for a_f in a_lf_im.iter() {
            bat::builder_add_compound_shape(&mut a_cf_im, a_f);
            // OCCT L3719-3735.
            for a_e in explorer(a_f, ShapeType::Edge, ShapeType::Shape) {
                match super::brep_offset_make_offset_1::dm_change_seek(&mut a_dmef, &a_e) {
                    Some(p_lf) => {
                        // OCCT L3731-3732: internal edges should not be
                        // used.
                        a_me_avoid.add(&a_e);
                        p_lf.push(a_f.clone());
                    }
                    None => {
                        shape_data_map::bind(&mut a_dmef, &a_e, vec![a_f.clone()]);
                    }
                }
            }
            // OCCT L3737-3749: fill the connection map of the vertices of
            // the inverted edges to the faces.
            for a_v in explorer(a_f, ShapeType::Vertex, ShapeType::Shape) {
                let entry = a_dmvf
                    .entry(shape_key_of(&a_v))
                    .or_insert((a_v.clone(), Vec::new()));
                append_to_list(&mut entry.1, a_f);
            }
        }

        // OCCT L3752-3758: for the splits to be regular they should form
        // only one block (the shell_splitter re-host of the
        // BOPTools_AlgoTools::MakeConnexityBlocks EDGE/FACE form,
        // architecture difference #44).
        let a_cf_faces = explorer(&a_cf_im, ShapeType::Face, ShapeType::Shape);
        let a_lcbf =
            crate::bop::algo::shell_splitter::make_connexity_blocks(&a_cf_faces);
        if a_lcbf.len() == 1 {
            continue;
        }

        // OCCT L3760-3778: check if the inverted edges create the
        // irregularity.
        let mut a_dmcbv_inverted: ShapeDataMap<OcctShapeSet> = HashMap::new();
        let mut a_dmcbv_all: ShapeDataMap<OcctShapeSet> = HashMap::new();
        // OCCT L3766-3778: the block walk over aLCBF (each block carried as
        // the compound of its faces).
        let a_lcbf_shapes: Vec<Shape> =
            a_lcbf.iter().map(|b| block_shape(&b.shapes)).collect();
        for a_cb in a_lcbf_shapes.iter() {
            // OCCT L3771-3777: check if it is possible to remove the block.
            if !a_bf.check_inverted_block(a_cb, &a_lcbf_shapes, &mut a_dmcbv_inverted, &mut a_dmcbv_all)
            {
                // OCCT L3775: none of the edges in this block should be
                // removed.
                map_shapes_indexed(a_cb, ShapeType::Edge, &mut a_me_avoid);
            }
        }
    }

    // OCCT L3781-3800: all edges not included in aMEAvoid can be removed.
    let mut a_me_rem: OcctShapeSet = HashMap::new();
    for i_inverted in 1..=a_bf.my_inverted_edges.extent() {
        let a_e = a_bf.my_inverted_edges.find_key_1(i_inverted).clone();
        if !a_me_avoid.contains(&a_e) {
            for a_v in bat::sub_shapes(&a_e) {
                let p_lf = shape_data_map::seek(&a_dmvf, &a_v);
                if let Some(p_lf) = p_lf {
                    if p_lf.len() > 3 {
                        set_add(&mut a_me_rem, &a_e);
                        break;
                    }
                }
            }
        }
    }

    // OCCT L3802-3805.
    if a_me_rem.is_empty() {
        return;
    }

    // OCCT L3807-3857: all invalid faces containing these edges can be
    // removed.
    let mut a_inv_faces: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let mut a_mf_rem: OcctShapeSet = HashMap::new();
    let mut a_mf_to_update = OcctIndexedShapeMap::new();
    let a_nb = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb {
        let (a_f, mut a_lf_im) = {
            let (_k, v) = a_bf.my_invalid_faces.get_index(i - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        // OCCT L3818-3851: the iterator-removal walk (the index form of
        // aLFIm.Remove(aIt)).
        let mut j = 0usize;
        while j < a_lf_im.len() {
            let a_f_im = a_lf_im[j].clone();
            // OCCT L3825-3835: to be removed the face should have at least
            // two not-connected inverted edges.
            let mut a_ce_inv = empty_compound();
            for a_e in explorer(&a_f_im, ShapeType::Edge, ShapeType::Shape) {
                if set_contains(&a_me_rem, &a_e) {
                    bat::builder_add_compound_shape(&mut a_ce_inv, &a_e);
                }
            }
            // OCCT L3838-3839: check connectivity (the
            // BOPTools_AlgoTools::MakeConnexityBlocks VERTEX/EDGE form).
            let a_ce_inv_edges = explorer(&a_ce_inv, ShapeType::Edge, ShapeType::Shape);
            let a_locations = [glam::DAffine3::IDENTITY];
            let a_lcbe = crate::bop::algo::wire_splitter::make_connexity_blocks(
                &a_ce_inv_edges,
                &a_locations,
            );
            // OCCT L3841-3850.
            if a_lcbe.len() >= 2 {
                a_mf_to_update.add(&a_f);
                set_add(&mut a_mf_rem, &a_f_im);
                a_lf_im.remove(j);
            } else {
                j += 1;
            }
        }
        // OCCT L3853-3856.
        if !a_lf_im.is_empty() {
            indexed_add_vec(&mut a_inv_faces, &a_f, a_lf_im);
        }
    }

    // OCCT L3859-3862.
    if a_mf_rem.is_empty() {
        return;
    }

    // OCCT L3864: myInvalidFaces = aInvFaces.
    a_bf.my_invalid_faces = a_inv_faces;

    // OCCT L3866-3886: remove from the splits.
    let a_nb = a_mf_to_update.extent();
    for i in 1..=a_nb {
        let a_f = a_mf_to_update.find_key_1(i).clone();
        let a_key = shape_key_of(&a_f);
        if !a_bf.my_of_images.contains_key(&a_key) {
            continue;
        }
        let mut a_lf_im: Vec<Shape> = a_bf
            .my_of_images
            .get(&a_key)
            .expect("myOFImages entry")
            .1
            .clone();
        // OCCT L3872-3885: the iterator-removal walk.
        let mut j = 0usize;
        while j < a_lf_im.len() {
            let a_f_im = a_lf_im[j].clone();
            if set_contains(&a_mf_rem, &a_f_im) {
                map_shapes_indexed(&a_f_im, ShapeType::Edge, the_me_removed);
                a_lf_im.remove(j);
            } else {
                j += 1;
            }
        }
        add_of_images_local(&mut a_bf.my_of_images, &a_key, a_lf_im);
    }
    let _ = my_oe_origins;
}

/// The OCCT `aInvFaces.Add(aF, aLFIm)` form over the
/// ShapeIndexedDataMap (the IndexedDataMap::Add replace-on-Add semantics).
fn indexed_add_vec(m: &mut ShapeIndexedDataMap<Vec<Shape>>, k: &Shape, v: Vec<Shape>) {
    let key = shape_key_of(k);
    if let Some(entry) = m.get_mut(&key) {
        entry.1 = v;
    } else {
        m.insert(key, (k.clone(), v));
    }
}

/// The local add_of_images keyed by the precomputed key.
fn add_of_images_local(
    m: &mut ShapeIndexedDataMap<Vec<Shape>>,
    key: &ShapeKey,
    v: Vec<Shape>,
) {
    if let Some(entry) = m.get_mut(key) {
        entry.1 = v;
    }
}
