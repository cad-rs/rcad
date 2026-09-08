// OCCT BRepOffset_MakeOffset_1.cxx L1301-2054 — module b of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L1175-2054 (the cxx order within the module
// keeps the OCCT line anchors):
//   - IntersectTrimmedEdges (cxx L1207-1292)
//   - checkConnectionsOfFace (cxx L1301-1334)
//   - BuildSplitsOfFaces (cxx L1342-1833)
//   - GetEdges (cxx L1839-1954)
//   - CheckIfArtificial (cxx L1960-2053)
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#48).

use std::collections::HashMap;

use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topo::topods::ShapeType;
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::builder::Builder;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, shape_indexed_data_map,
    OcctIndexedShapeMap, OcctShapeSet, ShapeDataMap, ShapeIndexedDataMap,
};
use super::brep_offset_make_offset_1::{
    add_of_images, build_splits_of_face, empty_compound, map_shapes_and_ancestors_map,
    map_shapes_indexed, process_micro_edge, shape_key_of, BRepOffsetBuildOffsetFaces,
    DataMapOfShapeIndexedMapOfShape,
};

// ---------------------------------------------------------------------------
// OCCT UpdateIntersectedEdges (cxx L1175-1203).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::UpdateIntersectedEdges (cxx
/// L1175-1203) — saving connection from trimmed edges to not trimmed ones.
pub(crate) fn update_intersected_edges(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_la: &[Shape],
    the_gf: &Builder,
) {
    // The myETrimEInf member is taken in the owned Option form
    // (architecture difference #40); the writes go through the local
    // snapshot of the map.
    let mut a_e_trim_e_inf = a_bf.my_e_trim_e_inf.clone().unwrap_or_default();
    for a_s in the_la.iter() {
        // OCCT L1182: pEInf = myETrimEInf->Seek(aS).
        let p_e_inf = match shape_data_map::seek(&a_e_trim_e_inf, a_s) {
            Some(v) => v.clone(),
            None => continue,
        };

        // OCCT L1188: aLSIm = theGF.Modified(aS).
        let a_ls_im =
            super::brep_offset_make_offset_1::builder_modified(the_gf, a_s);
        if a_ls_im.is_empty() {
            continue;
        }

        for a_e_im in a_ls_im.iter() {
            // OCCT L1197-1200: if (!myETrimEInf->IsBound(aEIm))
            //   myETrimEInf->Bind(aEIm, *pEInf).
            if !shape_data_map::is_bound(&a_e_trim_e_inf, a_e_im) {
                shape_data_map::bind(&mut a_e_trim_e_inf, a_e_im, p_e_inf.clone());
            }
        }
    }
    a_bf.my_e_trim_e_inf = Some(a_e_trim_e_inf);
}

// ---------------------------------------------------------------------------
// OCCT IntersectTrimmedEdges (cxx L1207-1292).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::IntersectTrimmedEdges (cxx
/// L1207-1292) — intersection of the trimmed edges among themselves.
pub(crate) fn intersect_trimmed_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    _the_range: &ProgressScope,
) {
    // OCCT L1210: get edges to intersect from descendants of the offset
    // faces.
    let mut a_ls: Vec<Shape> = Vec::new();

    // OCCT L1212-1213: Message_ProgressScope aPS(theRange, nullptr, 2);
    // iterator over *myFaces.
    let a_faces: Vec<Shape> = a_bf.my_faces.as_ref().unwrap().clone();
    for a_f in a_faces.iter() {
        let _a_ps = ProgressScope::new(&NoopProgress, "", 2);
        // OCCT L1222: aLE = myAsDes->Descendant(aF).
        let a_le: Vec<Shape> = a_bf
            .my_as_des
            .as_ref()
            .expect("myAsDes")
            .borrow()
            .descendant(a_f)
            .to_vec();
        for a_e_cur in a_le.iter() {
            let a_e = a_e_cur.clone();
            // OCCT L1228-1231: if (ProcessMicroEdge(aE, myContext))
            // continue.
            if process_micro_edge(&a_e) {
                continue;
            }
            // OCCT L1233-1236: if (myModifiedEdges.Add(aE)) aLS.Append(aE).
            if set_add(&mut a_bf.my_modified_edges, &a_e) {
                a_ls.push(a_e);
            }
        }
    }

    // OCCT L1240-1244: if (aLS.Extent() < 2) return — nothing to intersect.
    if a_ls.len() < 2 {
        return;
    }

    // OCCT L1247-1253: BOPAlgo_Builder aGFE; SetArguments(aLS);
    // Perform(aPS.Next()); if (aGFE.HasErrors()) return (architecture
    // difference #41 — the PaveFiller + Builder pipeline).
    let mut a_gf_filler = PaveFiller::new();
    a_gf_filler.set_arguments(a_ls.clone());
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "BOPAlgo_Builder", 1);
    a_gf_filler.perform(&a_ps);
    let mut a_gf = Builder::new(
        a_gf_filler.ds(),
        crate::bop::algo::builder::BooleanOpType::Union,
        a_gf_filler.fuzzy_value(),
    );
    a_gf.my_arguments = a_gf_filler.ds().arguments.clone();
    if a_gf.build().is_err() || a_gf.has_errors() {
        return;
    }

    // OCCT L1255-1288: aLA; fill map with edges images.
    let mut a_la: Vec<Shape> = Vec::new();
    for a_it in a_ls.iter() {
        let a_e = a_it;
        // OCCT L1265: aLEIm = aGFE.Modified(aE).
        let a_le_im = super::brep_offset_make_offset_1::builder_modified(&a_gf, a_e);
        if a_le_im.is_empty() {
            continue;
        }
        // OCCT L1271: aLA.Append(aE).
        a_la.push(a_e.clone());
        // OCCT L1273: myOEImages.Bind(aE, aLEIm).
        shape_data_map::bind(&mut a_bf.my_oe_images, a_e, a_le_im.clone());
        // OCCT L1275-1287: save origins.
        for a_e_im in a_le_im.iter() {
            // OCCT L1279-1286: myOEOrigins.ChangeSeek/Bound + Append.
            super::brep_offset_make_offset_1::add_to_container_map_list(
                &mut a_bf.my_oe_origins,
                a_e_im,
                a_e,
            );
        }
    }

    // OCCT L1290-1291: UpdateOrigins(aLA, *myEdgesOrigins, aGFE);
    // UpdateIntersectedEdges(aLA, aGFE).
    let mut an_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();
    super::brep_offset_make_offset_1::update_origins(&a_la, &mut an_edges_origins, &a_gf);
    a_bf.my_edges_origins = Some(an_edges_origins);
    update_intersected_edges(a_bf, &a_la, &a_gf);
}

// ---------------------------------------------------------------------------
// OCCT checkConnectionsOfFace (cxx L1301-1334).
// ---------------------------------------------------------------------------

/// OCCT static checkConnectionsOfFace (cxx L1301-1334) — checks the number
/// of connections for theFace with theLF; returns true if the number of
/// connections is more than 1.
pub(crate) fn check_connections_of_face(the_face: &Shape, the_lf: &[Shape]) -> bool {
    let mut a_shape_vert = OcctIndexedShapeMap::new();
    for a_shape in the_lf.iter() {
        if a_shape.is_same(the_face) {
            continue;
        }
        // OCCT L1313: TopExp::MapShapes(aShape, TopAbs_VERTEX, aShapeVert).
        map_shapes_indexed(a_shape, ShapeType::Vertex, &mut a_shape_vert);
    }
    let mut a_nb_connections = 0usize;
    let mut a_face_vertices = OcctIndexedShapeMap::new();
    // OCCT L1317: TopExp::MapShapes(theFace, TopAbs_VERTEX, aFaceVertices).
    map_shapes_indexed(the_face, ShapeType::Vertex, &mut a_face_vertices);
    for a_vert in a_face_vertices.iter() {
        if a_shape_vert.contains(a_vert) {
            a_nb_connections += 1;
        }
        if a_nb_connections > 1 {
            return true;
        }
    }
    false
}

// ---------------------------------------------------------------------------
// OCCT BuildSplitsOfFaces (cxx L1342-1833).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::BuildSplitsOfFaces (cxx L1342-1833) —
/// building the splits of offset faces and looking for the invalid splits.
pub(crate) fn build_splits_of_faces_impl(a_bf: &mut BRepOffsetBuildOffsetFaces, _the_range: &ProgressScope) {
    // OCCT L1344: BRep_Builder aBB.
    // OCCT L1345: int i, aNb.
    //
    // OCCT L1348: processed faces.
    let mut a_lf_done: Vec<Shape> = Vec::new();
    // OCCT L1350-1353: extended face - map of neutral edges, i.e. in one
    // split - valid and in other - invalid.
    let mut a_dmfmne: ShapeDataMap<OcctShapeSet> = HashMap::new();
    // OCCT L1355-1358: map of valid edges for each face.
    let mut a_dmfmve: ShapeDataMap<OcctShapeSet> = HashMap::new();
    // OCCT L1360: map of invalid edges for each face.
    let mut a_dmfmie: DataMapOfShapeIndexedMapOfShape = HashMap::new();
    // OCCT L1362-1365: map of valid inverted edges for the face.
    let mut a_dmfmvie: ShapeDataMap<OcctShapeSet> = HashMap::new();
    // OCCT L1367: map of splits to check for internals.
    let mut a_mf_to_check_int = OcctIndexedShapeMap::new();
    // OCCT L1369: map of edges created from vertex and marked as invalid.
    let mut a_m_edge_invalid_by_vertex: OcctShapeSet = HashMap::new();
    // OCCT L1371: map of edges created from vertex and marked as valid.
    let mut a_m_edge_valid_by_vertex: OcctShapeSet = HashMap::new();
    // OCCT L1373-1374: connection map from old edges to new ones.
    let mut a_dmeorleim: ShapeDataMap<Vec<Shape>> = HashMap::new();

    // OCCT L1377-1381: the outer scopes and the face loop.
    let a_faces: Vec<Shape> = a_bf.my_faces.as_ref().unwrap().clone();
    for a_f in a_faces.iter() {
        let _a_psbf = ProgressScope::new(&NoopProgress, "Building faces", 2 * a_faces.len());

        // OCCT L1389-1393: pLFIm = myOFImages.ChangeSeek(aF);
        // if (pLFIm && pLFIm->IsEmpty()) continue.  (The mutable pointer
        // form of the OCCT code is re-looked-up at the use points below —
        // the map is not restructured between them.)
        let a_key = shape_key_of(a_f);
        if let Some(v) = a_bf.my_of_images.get(&a_key) {
            if v.1.is_empty() {
                continue;
            }
        }

        // OCCT L1395-1401: get edges by which the face should be split.
        let mut a_ce = Shape::null();
        let mut a_map_e_inv = OcctIndexedShapeMap::new();
        let b_found = a_bf.get_edges(a_f, &mut a_ce, Some(&mut a_map_e_inv));
        if !b_found {
            continue;
        }

        // OCCT L1412-1413: build splits.
        let mut a_lf_images: Vec<Shape> = Vec::new();
        let mut a_faces_origins = a_bf.my_faces_origins.clone().unwrap_or_default();
        build_splits_of_face(a_f, &a_ce, &mut a_faces_origins, &mut a_lf_images);
        a_bf.my_faces_origins = Some(a_faces_origins);

        // OCCT L1415-1522: the invalid-edge rebuild attempt.
        if !a_map_e_inv.is_empty() {
            // OCCT L1417-1420: check if all possible faces are built.
            let mut a_men_inv: OcctShapeSet = HashMap::new();
            let mut b_artificial_case =
                a_lf_images.is_empty()
                    || a_bf.check_if_artificial(a_f, &a_lf_images, &a_ce, &a_map_e_inv, &mut a_men_inv);

            // OCCT L1422-1429: try to build splits using invalid edges —
            // aCE1 = compound { aCE } + all aMapEInv edges.
            let mut a_ce1 = empty_compound();
            bat::builder_add_compound_shape(&mut a_ce1, &a_ce);
            for i in 1..=a_map_e_inv.extent() {
                let a_e_inv = a_map_e_inv.find_key_1(i).clone();
                bat::builder_add_compound_shape(&mut a_ce1, &a_e_inv);
            }

            let mut a_lf_images1: Vec<Shape> = Vec::new();
            let mut a_faces_origins = a_bf.my_faces_origins.clone().unwrap_or_default();
            build_splits_of_face(a_f, &a_ce1, &mut a_faces_origins, &mut a_lf_images1);
            a_bf.my_faces_origins = Some(a_faces_origins);

            // OCCT L1434-1476: check if the rebuilding has added some new
            // faces to the splits.  (The NCollection_List iterator removal
            // form is the index loop of the OCCT Remove(aItLFIm).)
            let mut i_f = 0usize;
            while i_f < a_lf_images1.len() {
                let mut b_all_inv = true;
                // OCCT L1438-1441: additional check for artificial case —
                // if the current image face consists only of edges from
                // aMapEInv and aMENInv then recheck the current face for
                // the further processing.
                let mut a_to_re_check_face = b_artificial_case;
                let a_f_im = a_lf_images1[i_f].clone();
                // OCCT L1443-1456: the edge walk; b_broke mirrors the
                // aExpE.More() state after the break.
                let mut b_broke = false;
                for a_e in explorer(&a_f_im, ShapeType::Edge, ShapeType::Shape) {
                    if !a_map_e_inv.contains(&a_e) {
                        b_all_inv = false;
                        if !set_contains(&a_men_inv, &a_e) {
                            a_to_re_check_face = false;
                            b_broke = true;
                            break;
                        }
                    }
                }
                // OCCT L1459-1462: if the current image face is to recheck
                // then check the number of connections for this face with
                // the other image faces of the current face.
                if !b_broke && a_to_re_check_face {
                    a_to_re_check_face = check_connections_of_face(&a_f_im, &a_lf_images1);
                }
                // OCCT L1464-1475: do not delete the image face from the
                // further processing if aToReCheckFace is true.
                if !b_broke && !a_to_re_check_face {
                    if b_all_inv {
                        a_mf_to_check_int.add(&a_f_im);
                    }
                    a_lf_images1.remove(i_f);
                } else {
                    i_f += 1;
                }
            }

            // OCCT L1478-1488: if (bArtificialCase) { if (extent equal)
            // bArtificialCase = false; else aLFImages = aLFImages1; }.
            if b_artificial_case {
                if a_lf_images.len() == a_lf_images1.len() {
                    b_artificial_case = false;
                } else {
                    a_lf_images = a_lf_images1.clone();
                }
            }

            // OCCT L1490-1521: if (bArtificialCase) — make the face
            // invalid.
            if b_artificial_case {
                let mut a_me_inv = OcctIndexedShapeMap::new();

                // OCCT L1495: *pLFIm = aLFImages.
                add_of_images(&mut a_bf.my_of_images, a_f, a_lf_images.clone());
                for a_f_im in a_lf_images.iter() {
                    for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                        if a_map_e_inv.contains(&a_e) {
                            // OCCT L1506-1507.
                            a_bf.my_invalid_edges.add(&a_e);
                            a_me_inv.add(&a_e);
                        } else {
                            // OCCT L1511: myValidEdges.Add(aE).
                            a_bf.my_valid_edges.add(&a_e);
                        }
                    }
                }
                // OCCT L1516-1517: myArtInvalidFaces.Bind(aF, aMEInv);
                // aDMFMIE.Bind(aF, aMEInv).
                shape_data_map::bind(&mut a_bf.my_art_invalid_faces, a_f, a_me_inv);
                let a_me_inv_ref = shape_data_map::value(&a_bf.my_art_invalid_faces, a_f);
                let a_me_inv_copy = super::brep_offset_make_offset_1::clone_indexed_map(a_me_inv_ref);
                shape_data_map::bind(&mut a_dmfmie, a_f, a_me_inv_copy);
                // OCCT L1518: aLFDone.Append(aF).
                a_lf_done.push(a_f.clone());
                // OCCT L1520: continue.
                continue;
            }
        }

        // OCCT L1524-1534: find invalid edges.
        a_bf.find_invalid_edges_per_face(
            a_f,
            &a_lf_images,
            &mut a_dmfmve,
            &mut a_dmfmne,
            &mut a_dmfmie,
            &mut a_dmfmvie,
            &mut a_dmeorleim,
            &mut a_m_edge_invalid_by_vertex,
            &mut a_m_edge_valid_by_vertex,
            _the_range,
        );

        // OCCT L1536-1545: save the new splits.
        let a_key = shape_key_of(a_f);
        if !a_bf.my_of_images.contains_key(&a_key) {
            a_bf.my_of_images.insert(a_key, (a_f.clone(), Vec::new()));
        }
        let a_idx = a_bf
            .my_of_images
            .get_full(&a_key)
            .map(|(i, _, _)| i + 1)
            .expect("myOFImages: the face was just added");
        {
            let p_lf_im =
                super::brep_offset_make_offset_1::shape_indexed_data_map_change_find_1(
                    &mut a_bf.my_of_images,
                    a_idx,
                );
            p_lf_im.clear();
            p_lf_im.extend(a_lf_images.iter().cloned());
        }

        // OCCT L1547: aLFDone.Append(aF).
        a_lf_done.push(a_f.clone());
    }

    // OCCT L1550-1553: if (myInvalidEdges.IsEmpty() && myArtInvalidFaces.
    // IsEmpty() && aDMFMIE.IsEmpty()) return.
    if a_bf.my_invalid_edges.is_empty()
        && a_bf.my_art_invalid_faces.is_empty()
        && a_dmfmie.is_empty()
    {
        return;
    }

    // OCCT L1555-1557: additional step to find invalid edges by checking
    // unclassified edges in the splits of SD faces.
    a_bf.find_invalid_edges_list(&a_lf_done, &mut a_dmfmie, &mut a_dmfmve, &mut a_dmfmne);

    // OCCT L1559-1561: additional step to mark inverted edges located
    // inside loops of invalid edges as invalid as well.
    a_bf.make_inverted_edges_invalid(&a_lf_done);

    // OCCT L1563-1595: the OFFSET_DEBUG compound blocks are diagnostics
    // only (the OCCT #ifdef body) — omitted under the rcad no-debug form.

    // OCCT L1597-1613: build Edge-Face connectivity map to find faces whose
    // removal may potentially lead to creation of holes in the faces.
    let mut an_ef_map: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let a_nb_f = a_bf.my_of_images.len();
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
                &mut an_ef_map,
            );
        }
    }

    // OCCT L1615: anEmptyMap.
    let an_empty_map: OcctShapeSet = HashMap::new();
    // OCCT L1617: invalid faces inside the holes.
    let mut a_mf_inv_in_hole = OcctIndexedShapeMap::new();
    // OCCT L1619-1620: all hole faces.
    let mut a_f_holes = empty_compound();
    // OCCT L1622: faces containing only the inverted edges and the invalid
    // ones.
    let mut an_inverted_faces: Vec<Shape> = Vec::new();

    let _a_psif = ProgressScope::new(&NoopProgress, "Checking validity of faces", a_lf_done.len());
    // OCCT L1627-1691: find invalid faces — considering faces containing
    // only invalid edges as invalid.
    for a_f_cur in a_lf_done.iter() {
        let a_f = a_f_cur.clone();
        // OCCT L1635: aLFImages = myOFImages.ChangeFromKey(aF) — the
        // snapshot form (the mutable reference is re-bound after the
        // FindInvalidFaces call that may clear the list).
        let a_key = shape_key_of(&a_f);
        let mut a_lf_images: Vec<Shape> = match a_bf.my_of_images.get(&a_key) {
            Some(v) => v.1.clone(),
            None => continue,
        };

        let mut a_lf_inv: Vec<Shape> = Vec::new();
        // OCCT L1638: bArtificialCase = myArtInvalidFaces.IsBound(aF).
        let b_artificial_case = shape_data_map::is_bound(&a_bf.my_art_invalid_faces, &a_f);
        if b_artificial_case {
            // OCCT L1641: aLFInv = aLFImages.
            a_lf_inv = a_lf_images.clone();
        } else {
            // OCCT L1646-1650: neutral edges.
            let a_mne_in_map = shape_data_map::seek(&a_dmfmne, &a_f).is_some();
            let _a_p_mne: &OcctShapeSet = if a_mne_in_map {
                shape_data_map::value(&a_dmfmne, &a_f)
            } else {
                &an_empty_map
            };

            // OCCT L1652-1654: find faces inside holes wires.
            let mut a_mf_holes: OcctShapeSet = HashMap::new();
            let my_faces_origins = a_bf.my_faces_origins.clone().unwrap_or_default();
            let a_f_or = shape_data_map::find(&my_faces_origins, &a_f);
            a_bf.find_faces_inside_hole_wires(
                &a_f_or,
                &a_f,
                &a_lf_images,
                &a_dmeorleim,
                &an_ef_map,
                &mut a_mf_holes,
            );
            // OCCT L1656-1660: aBB.Add(aFHoles, hole).
            for a_hole in a_mf_holes.values() {
                bat::builder_add_compound_shape(&mut a_f_holes, a_hole);
            }

            // OCCT L1662-1672: find invalid faces.
            let a_mne = if a_mne_in_map {
                let v = shape_data_map::find(&a_dmfmne, &a_f);
                v
            } else {
                HashMap::new()
            };
            a_bf.find_invalid_faces(
                &mut a_lf_images,
                &a_dmfmve,
                &a_dmfmie,
                &a_mne,
                &a_m_edge_invalid_by_vertex,
                &a_m_edge_valid_by_vertex,
                &a_mf_holes,
                &mut a_mf_inv_in_hole,
                &mut a_lf_inv,
                &mut an_inverted_faces,
            );
            // OCCT L1635 (write-back): aLFImages = myOFImages.
            add_of_images(&mut a_bf.my_of_images, &a_f, a_lf_images.clone());
        }

        // OCCT L1675-1690.
        if !a_lf_inv.is_empty() {
            if shape_data_map::is_bound(&a_bf.my_already_inv_faces, &a_f) {
                if *shape_data_map::value(&a_bf.my_already_inv_faces, &a_f) > 2 {
                    if a_lf_inv.len() == a_lf_images.len() && !b_artificial_case {
                        // OCCT L1683: aLFImages.Clear().
                        a_lf_images.clear();
                        add_of_images(&mut a_bf.my_of_images, &a_f, a_lf_images.clone());
                    }
                    // OCCT L1686: aLFInv.Clear().
                    a_lf_inv.clear();
                }
            }
            // OCCT L1689: myInvalidFaces.Add(aF, aLFInv).
            indexed_data_map_add(&mut a_bf.my_invalid_faces, &a_f, a_lf_inv);
        }
    }

    // OCCT L1693-1697: if (myInvalidFaces.IsEmpty()) { myInvalidEdges.
    // Clear(); return; }.
    if a_bf.my_invalid_faces.is_empty() {
        a_bf.my_invalid_edges.clear();
        return;
    }

    // OCCT L1716-1723: remove invalid splits of faces using inverted
    // edges.
    let mut a_me_removed = OcctIndexedShapeMap::new();
    a_bf.remove_invalid_splits_by_inverted_edges(&mut a_me_removed);
    if a_bf.my_invalid_faces.is_empty() {
        a_bf.my_invalid_edges.clear();
        return;
    }

    // OCCT L1726: remove invalid splits from valid splits.
    a_bf.remove_invalid_splits_from_valid(&a_dmfmvie);

    // OCCT L1729-1735: remove inside faces.
    a_bf.remove_inside_faces(
        &an_inverted_faces,
        &a_mf_to_check_int,
        &a_mf_inv_in_hole,
        &a_f_holes,
        &mut a_me_removed,
        _the_range,
    );

    // OCCT L1737-1751: make compound of valid splits.
    let mut a_cf_im = empty_compound();
    let a_nb = a_bf.my_of_images.len();
    for i in 1..=a_nb {
        let a_lf_im: Vec<Shape> = a_bf
            .my_of_images
            .get_index(i - 1)
            .expect("myOFImages index")
            .1
             .1
            .clone();
        for a_f_im in a_lf_im.iter() {
            bat::builder_add_compound_shape(&mut a_cf_im, a_f_im);
        }
    }

    // OCCT L1753-1755: aDMEF; MapShapesAndAncestors(aCFIm, EDGE, FACE,
    // aDMEF).
    let mut a_dmef: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    map_shapes_and_ancestors_map(&a_cf_im, ShapeType::Edge, ShapeType::Face, &mut a_dmef);

    // OCCT L1757-1758: filter maps of images and origins.
    a_bf.filter_edges_images(&a_cf_im);

    // OCCT L1761: aMERemoved.Extent() ? myInsideEdges : aMERemoved.
    if !a_me_removed.is_empty() {
        let my_inside_edges = super::brep_offset_make_offset_1::clone_indexed_map(&a_bf.my_inside_edges);
        a_bf.filter_invalid_faces(&a_dmef, &my_inside_edges);
    } else {
        a_bf.filter_invalid_faces(&a_dmef, &a_me_removed);
    }

    // OCCT L1762-1767.
    let a_nb = a_bf.my_invalid_faces.len();
    if a_nb == 0 {
        a_bf.my_invalid_edges.clear();
        return;
    }

    // OCCT L1790: aMERemoved.Extent() ? myInsideEdges : aMERemoved.
    let mut a_me_use_in_rebuild: OcctShapeSet = HashMap::new();
    if !a_me_removed.is_empty() {
        let my_inside_edges = super::brep_offset_make_offset_1::clone_indexed_map(&a_bf.my_inside_edges);
        a_bf.filter_invalid_edges(
            &a_dmfmie,
            &a_me_removed,
            &my_inside_edges,
            &mut a_me_use_in_rebuild,
        );
    } else {
        a_bf.filter_invalid_edges(
            &a_dmfmie,
            &a_me_removed,
            &a_me_removed,
            &mut a_me_use_in_rebuild,
        );
    }

    // OCCT L1793-1794: check additionally validity of edges originated
    // from vertices.
    a_bf.check_edges_created_by_vertex();

    // OCCT L1808-1818.
    a_bf.my_last_inv_edges.clear();
    let a_nb = a_bf.my_invalid_edges.extent();
    for i in 1..=a_nb {
        let a_e = a_bf.my_invalid_edges.find_key_1(i).clone();
        set_add(&mut a_bf.my_last_inv_edges, &a_e);
        if !set_contains(&a_me_use_in_rebuild, &a_e) {
            a_bf.my_edges_to_avoid.add(&a_e);
        }
    }

    // OCCT L1820-1832.
    let a_nb = a_bf.my_invalid_faces.len();
    for i in 1..=a_nb {
        let a_f = shape_indexed_data_map::find_key_1(&a_bf.my_invalid_faces, i).clone();
        if shape_data_map::is_bound(&a_bf.my_already_inv_faces, &a_f) {
            *shape_data_map::change_find(&mut a_bf.my_already_inv_faces, &a_f) += 1;
        } else {
            shape_data_map::bind(&mut a_bf.my_already_inv_faces, &a_f, 1);
        }
    }
}

/// The OCCT NCollection_IndexedDataMap::Add form — binds when absent and
/// REPLACES the value when the key is already present (the
/// shape_indexed_data_map::add keeps the old value; the cxx sites rely on
/// the OCCT replace-on-Add semantics).
pub(crate) fn indexed_data_map_add<V>(m: &mut ShapeIndexedDataMap<V>, k: &Shape, v: V) {
    let key = shape_key_of(k);
    if let Some(entry) = m.get_mut(&key) {
        entry.1 = v;
    } else {
        m.insert(key, (k.clone(), v));
    }
}

// ---------------------------------------------------------------------------
// OCCT GetEdges (cxx L1839-1954).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::GetEdges (cxx L1839-1954) — getting
/// edges from the AsDes map to build the splits of faces.
pub(crate) fn get_edges_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_face: &Shape,
    the_edges: &mut Shape,
    mut the_inv_map: Option<&mut OcctIndexedShapeMap>,
) -> bool {
    // OCCT L1845-1861: get boundary edges.
    let mut a_mf_bounds: OcctShapeSet = HashMap::new();
    for a_e in explorer(the_face, ShapeType::Edge, ShapeType::Shape) {
        // OCCT L1850: pLEIm = myOEImages.Seek(aE).
        match shape_data_map::seek(&a_bf.my_oe_images, &a_e) {
            Some(p_le_im) => {
                let images = p_le_im.clone();
                for a_it_le in images.iter() {
                    set_add(&mut a_mf_bounds, a_it_le);
                }
            }
            None => {
                // OCCT L1859: aMFBounds.Add(aE).
                set_add(&mut a_mf_bounds, &a_e);
            }
        }
    }

    // OCCT L1863-1864: bool bFound(false), bUpdate(false).
    let mut b_found = false;
    let mut b_update = false;

    // OCCT L1866-1867: the resulting edges — a compound.
    let mut an_edges = empty_compound();
    // OCCT L1869: fence map.
    let mut a_me_fence: OcctShapeSet = HashMap::new();
    // OCCT L1871: the edges by which the offset face should be split.
    let a_le: Vec<Shape> = a_bf
        .my_as_des
        .as_ref()
        .expect("myAsDes")
        .borrow()
        .descendant(the_face)
        .to_vec();
    for a_e_cur in a_le.iter() {
        let a_e = a_e_cur.clone();
        // OCCT L1877-1880.
        if !b_update {
            b_update = set_contains(&a_bf.my_modified_edges, &a_e);
        }

        // OCCT L1882: pLEIm = myOEImages.Seek(aE).
        match shape_data_map::seek(&a_bf.my_oe_images, &a_e).cloned() {
            Some(p_le_im) => {
                for a_e_im_cur in p_le_im.iter() {
                    let a_e_im = a_e_im_cur.clone();
                    // OCCT L1890-1893: if (!aMEFence.Add(aEIm)) continue.
                    if !set_add(&mut a_me_fence, &a_e_im) {
                        continue;
                    }

                    // OCCT L1895-1906: myEdgesToAvoid.
                    if a_bf.my_edges_to_avoid.contains(&a_e_im) {
                        if let Some(the_inv_map) = the_inv_map.as_deref_mut() {
                            the_inv_map.add(&a_e_im);
                        }
                        if !b_update {
                            b_update = set_contains(&a_bf.my_last_inv_edges, &a_e_im);
                        }
                        continue;
                    }
                    // OCCT L1908-1911: check for micro edge.
                    if process_micro_edge(&a_e_im) {
                        continue;
                    }
                    // OCCT L1913: aBB.Add(anEdges, aEIm).
                    bat::builder_add_compound_shape(&mut an_edges, &a_e_im);
                    if !b_found {
                        b_found = !set_contains(&a_mf_bounds, &a_e_im);
                    }
                    // OCCT L1919-1922.
                    if !b_update {
                        b_update = set_contains(&a_bf.my_modified_edges, &a_e_im);
                    }
                }
            }
            None => {
                // OCCT L1927-1938: myEdgesToAvoid.
                if a_bf.my_edges_to_avoid.contains(&a_e) {
                    if let Some(the_inv_map) = the_inv_map.as_deref_mut() {
                        the_inv_map.add(&a_e);
                    }
                    if !b_update {
                        b_update = set_contains(&a_bf.my_last_inv_edges, &a_e);
                    }
                    continue;
                }
                // OCCT L1940-1943: check for micro edge.
                if process_micro_edge(&a_e) {
                    continue;
                }
                // OCCT L1944: aBB.Add(anEdges, aE).
                bat::builder_add_compound_shape(&mut an_edges, &a_e);
                if !b_found {
                    b_found = !set_contains(&a_mf_bounds, &a_e);
                }
            }
        }
    }

    // OCCT L1952-1953: theEdges = anEdges; return bFound && bUpdate.
    *the_edges = an_edges;
    b_found && b_update
}

// ---------------------------------------------------------------------------
// OCCT CheckIfArtificial (cxx L1960-2053).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::CheckIfArtificial (cxx L1960-2053) —
/// checks if the face is artificially invalid.
pub(crate) fn check_if_artificial_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_f: &Shape,
    the_lf_images: &[Shape],
    the_ce: &Shape,
    the_map_e_inv: &OcctIndexedShapeMap,
    the_men_inv: &mut OcctShapeSet,
) -> bool {
    // OCCT L1967-1975: all boundary edges should be used.
    let mut a_me_used = OcctIndexedShapeMap::new();
    for a_f_im in the_lf_images.iter() {
        map_shapes_indexed(a_f_im, ShapeType::Edge, &mut a_me_used);
        map_shapes_indexed(a_f_im, ShapeType::Vertex, &mut a_me_used);
    }

    // OCCT L1977-1979: aMVE — the vertex ancestors in theCE.
    let mut a_mve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    map_shapes_and_ancestors_map(&*the_ce, ShapeType::Vertex, ShapeType::Edge, &mut a_mve);

    // OCCT L1981-2003.
    let a_nb = the_map_e_inv.extent();
    for i in 1..=a_nb {
        let a_e_inv = the_map_e_inv.find_key_1(i).clone();
        for a_ve_inv in explorer(&a_e_inv, ShapeType::Vertex, ShapeType::Shape) {
            let p_le_n_inv = a_mve
                .get(&shape_key_of(&a_ve_inv))
                .map(|e| e.1.clone());
            if let Some(p_le_n_inv) = p_le_n_inv {
                for a_en_inv in p_le_n_inv.iter() {
                    if !a_me_used.contains(a_en_inv) {
                        set_add(the_men_inv, a_en_inv);
                    }
                }
            }
        }
    }

    // OCCT L2005-2008: if (theMENInv.IsEmpty()) return false.
    if the_men_inv.is_empty() {
        return false;
    }

    // OCCT L2010-2011: aMEFound — the edges of theCE.
    let mut a_me_found = OcctIndexedShapeMap::new();
    map_shapes_indexed(the_ce, ShapeType::Edge, &mut a_me_found);

    // OCCT L2013: aLE = myAsDes->Descendant(theF).
    let a_le: Vec<Shape> = a_bf
        .my_as_des
        .as_ref()
        .expect("myAsDes")
        .borrow()
        .descendant(the_f)
        .to_vec();

    // OCCT L2015-2050: the edge walk (the iterator-break form of the OCCT
    // code with the b_broke flags).
    let mut b_broke = false;
    'outer: for a_e_cur in a_le.iter() {
        let a_e = a_e_cur.clone();
        // OCCT L2019: pLEIm = myOEImages.Seek(aE).
        match shape_data_map::seek(&a_bf.my_oe_images, &a_e).cloned() {
            Some(p_le_im) => {
                // OCCT L2021: bChecked = false.
                let mut b_checked = false;
                let mut b_broke_inner = false;
                for a_e_im_cur in p_le_im.iter() {
                    let a_e_im = a_e_im_cur.clone();
                    // OCCT L2026-2029.
                    if !a_me_found.contains(&a_e_im) || set_contains(the_men_inv, &a_e_im) {
                        continue;
                    }
                    // OCCT L2031-2035.
                    b_checked = true;
                    if a_me_used.contains(&a_e_im) {
                        b_broke_inner = true;
                        break;
                    }
                }
                // OCCT L2038-2041: if (bChecked && !aItLEIm.More()) break.
                if b_checked && !b_broke_inner {
                    b_broke = true;
                    break 'outer;
                }
            }
            None => {
                // OCCT L2043-2049.
                if a_me_found.contains(&a_e)
                    && !set_contains(the_men_inv, &a_e)
                    && !a_me_used.contains(&a_e)
                {
                    b_broke = true;
                    break 'outer;
                }
            }
        }
    }

    // OCCT L2052: return aItLE.More().
    b_broke
}
