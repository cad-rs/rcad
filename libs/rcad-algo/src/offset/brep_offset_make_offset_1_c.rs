// OCCT BRepOffset_MakeOffset_1.cxx L2057-3083 — module c of the 1:1
// translation (split from brep_offset_make_offset_1.rs to respect the
// 2000-line file limit).
//
// This module carries cxx L2057-3083:
//   - FindInvalidEdges, the per-face form (cxx L2057-2528)
//   - addAsNeutral (cxx L2536-2561)
//   - FindInvalidEdges, the list form (cxx L2569-2769)
//   - MakeInvertedEdgesInvalid (cxx L2775-2850)
//   - FindInvalidFaces (cxx L2856-3083)
//
// Local dependency carriers of this module:
//   - BOPTools_Set (TKBO, BOPTools_Set.cxx L46-166) — the local
//     BOPToolsSet re-host (Add / IsEqual forms; architecture difference
//     #45 of brep_offset_make_offset_1.rs).
//   - BOPTools_AlgoTools3D::GetNormalToFaceOnEdge (TKBO,
//     BOPTools_AlgoTools3D.cxx L351-376) — the local re-host with the
//     pcurve + surface-derivatives body (the bop Builder re-host form).
//
// The architecture-difference numbering continues in
// brep_offset_make_offset_1.rs (#38-#48).

use std::collections::HashMap;

use rcad_kernel::core::message::ProgressScope;
use rcad_kernel::geom::{Curve2dEval, SurfaceEval};
use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;
use crate::feat::brep_feat_builder::explorer;

use super::brep_offset_make_offset_1::{
    append_to_list, find_shape, get_average_tangent, map_shapes_and_ancestors_map,
    map_shapes_indexed, nb_points, shape_key_of, BRepOffsetBuildOffsetFaces,
    DataMapOfShapeIndexedMapOfShape,
};
use super::brep_offset_tool::{
    set_add, set_contains, shape_data_map, shape_indexed_data_map, OcctIndexedShapeMap,
    OcctShapeSet, ShapeDataMap, ShapeIndexedDataMap,
};

// ---------------------------------------------------------------------------
// Local dependency carriers.
// ---------------------------------------------------------------------------

/// OCCT BOPTools_Set (BOPTools_Set.cxx L46-166) — the set of sub-shapes
/// used for the equality check of the shapes by their content
/// (architecture difference #45).
pub(crate) struct BOPToolsSet {
    my_shapes: Vec<Shape>, // OCCT: myShapes (NCollection_List)
    my_nb_shapes: usize,   // OCCT: myNbShapes
}

impl BOPToolsSet {
    /// OCCT BOPTools_Set::Add(theS, theType) (BOPTools_Set.cxx L124-166).
    pub fn add(&mut self, the_s: &Shape, the_type: ShapeType) {
        self.my_shapes.clear();
        self.my_nb_shapes = 0;
        for a_sx in explorer(the_s, the_type, ShapeType::Shape) {
            // OCCT L137-141: degenerated edges are skipped.
            if the_type == ShapeType::Edge && edge_is_degenerated(&a_sx) {
                continue;
            }
            // OCCT L143-158: INTERNAL shapes are added twice (both
            // orientations).
            let a_or = a_sx.orientation;
            if a_or == Orientation::Internal {
                let mut a_sy = a_sx.clone();
                a_sy = bat::oriented(&a_sy, Orientation::Forward);
                self.my_shapes.push(a_sy.clone());
                a_sy = bat::oriented(&a_sy, Orientation::Reversed);
                self.my_shapes.push(a_sy);
            } else {
                self.my_shapes.push(a_sx);
            }
        }
        self.my_nb_shapes = self.my_shapes.len();
    }

    /// OCCT BOPTools_Set::IsEqual(theOther) (BOPTools_Set.cxx L86-110).
    pub fn is_equal(&self, the_other: &BOPToolsSet) -> bool {
        if the_other.my_nb_shapes != self.my_nb_shapes {
            return false;
        }
        let mut a_m1: OcctShapeSet = HashMap::new();
        for a_s in self.my_shapes.iter() {
            set_add(&mut a_m1, a_s);
        }
        for a_s in the_other.my_shapes.iter() {
            if !set_contains(&a_m1, a_s) {
                return false;
            }
        }
        true
    }
}

/// The BRep_Tool::Degenerated(E) probe of BOPTools_Set.cxx L139.
fn edge_is_degenerated(the_e: &Shape) -> bool {
    the_e
        .as_edge()
        .map(|ed| ed.degenerated)
        .unwrap_or(false)
}

/// OCCT BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(E, F, Dir)
/// (BOPTools_AlgoTools3D.cxx L351-376) — the 3-argument form: the surface
/// normal at the middle parameter of the edge, computed via the edge's
/// pcurve on the face and the surface first derivatives (the local re-host
/// of the bop Builder re-host body).
pub(crate) fn get_normal_to_face_on_edge(the_e: &Shape, the_f: &Shape) -> Option<glam::DVec3> {
    // OCCT L352-355: BRep_Tool::Range(aE, aT1, aT2); aT = (aT1 + aT2)/2.
    let (a_t1, a_t2) = bat::brep_tool_range(the_e);
    let a_t = 0.5 * (a_t1 + a_t2);
    let surf = the_f.as_face().and_then(|fd| fd.surface.clone())?;
    // OCCT L365: aC2D1 = BRep_Tool::CurveOnSurface(aE, aF1, aTolPC).
    let pc = bat::brep_tool_curve_on_surface(the_e, the_f);
    if let Some((pc, _, _)) = pc {
        // OCCT L367-369: aC2D1->D0(aT, aP2D).
        let uv = Curve2dEval::point_at(&pc, a_t);
        // OCCT L371-375: aDNF1 = aDD1U ^ aDD1V.
        let (_p, d1u, d1v) = surf.derivatives(uv.x, uv.y);
        let n = d1u.cross(d1v);
        if n.length_squared() >= 1e-24 {
            return Some(n.normalize());
        }
    }
    None
}

/// OCCT gp_Vec::IsParallel(theOther, theAngularTolerance) — the parallelism
/// probe of cxx L2258.
fn vec_is_parallel(a: glam::DVec3, b: glam::DVec3, the_angular_tolerance: f64) -> bool {
    let la = a.length();
    let lb = b.length();
    if la <= 0. || lb <= 0. {
        return false;
    }
    let cos_ang = a.dot(b) / (la * lb);
    let ang = cos_ang.clamp(-1., 1.).acos();
    ang <= the_angular_tolerance || ang >= std::f64::consts::PI - the_angular_tolerance
}

/// OCCT gp_Dir::IsEqual(theOther, theAngularTolerance) — the direction
/// equality probe of cxx L2713.
fn dir_is_equal(a: glam::DVec3, b: glam::DVec3, the_angular_tolerance: f64) -> bool {
    vec_is_parallel(a, b, the_angular_tolerance) && a.dot(b) > 0.
}

/// OCCT BRepTools::OuterWire(F) — the outer wire of the face (the rcad
/// face data carries the outer wire slot).
pub(crate) fn outer_wire_of(the_f: &Shape) -> Shape {
    the_f
        .as_face()
        .map(|fd| fd.outer_wire.clone())
        .unwrap_or_else(Shape::null)
}

// ---------------------------------------------------------------------------
// OCCT FindInvalidEdges, the per-face form (cxx L2057-2528).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FindInvalidEdges (cxx L2057-2528) —
/// the per-face form.  Edge is considered as invalid in the following
/// cases: 1. its orientation on the face has changed comparing to the
/// original edge and face; 2. the vertices of the edge have changed places
/// comparing to the original edge and face.
#[allow(clippy::too_many_arguments)]
pub(crate) fn find_invalid_edges_per_face_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_f: &Shape,
    the_lf_images: &[Shape],
    the_dmfmve: &mut ShapeDataMap<OcctShapeSet>,
    the_dmfmne: &mut ShapeDataMap<OcctShapeSet>,
    the_dmfmie: &mut DataMapOfShapeIndexedMapOfShape,
    the_dmfmvie: &mut ShapeDataMap<OcctShapeSet>,
    the_dmeorleim: &mut ShapeDataMap<Vec<Shape>>,
    the_edges_invalid_by_vertex: &mut OcctShapeSet,
    the_edges_valid_by_vertex: &mut OcctShapeSet,
    _the_range: &ProgressScope,
) {
    let my_edges_origins = a_bf.my_edges_origins.clone().unwrap_or_default();
    let my_oe_origins = a_bf.my_oe_origins.clone();
    let my_oe_images = a_bf.my_oe_images.clone();
    let my_faces_origins = a_bf.my_faces_origins.clone().unwrap_or_default();

    // OCCT L2085: original face.
    let a_f_or = shape_data_map::find(&my_faces_origins, the_f);
    // OCCT L2087: invalid edges.
    let mut a_me_inv = OcctIndexedShapeMap::new();
    // OCCT L2089: valid edges.
    let mut a_me_val: OcctShapeSet = HashMap::new();
    // OCCT L2091: internal edges.
    let mut a_me_int: OcctShapeSet = HashMap::new();
    //
    // OCCT L2094-2096: maps for checking the inverted edges.
    let mut a_dmve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let mut a_dmef: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    let mut a_m_edges = OcctIndexedShapeMap::new();
    // OCCT L2098-2099: back map from the original shapes to their offset
    // images.
    let mut an_images: ShapeDataMap<Vec<Shape>> = HashMap::new();

    // OCCT L2101-2162: the first walk over the splits.
    for a_f_im in the_lf_images.iter() {
        let _a_ps = ProgressScope::new(&NoopProgressTag, "", 2);
        // OCCT L2111-2114: keep all edges.
        for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
            a_m_edges.add(&a_e);
            // OCCT L2118-2124: keep connection from edges to faces.
            {
                let key = shape_key_of(&a_e);
                let entry = a_dmef.entry(key).or_insert((a_e.clone(), Vec::new()));
                append_to_list(&mut entry.1, a_f_im);
            }
            // OCCT L2126-2138: keep connection from vertices to edges
            // (TopoDS_Iterator over the edge).
            for a_v in bat::sub_shapes(&a_e) {
                let key = shape_key_of(&a_v);
                let entry = a_dmve.entry(key).or_insert((a_v.clone(), Vec::new()));
                append_to_list(&mut entry.1, &a_e);
            }

            // OCCT L2140-2160: back map from original edges to their offset
            // images.
            let p_l_or = match shape_data_map::seek(&my_edges_origins, &a_e) {
                Some(v) => v.clone(),
                None => continue,
            };
            for a_s_or in p_l_or.iter() {
                let mut a_s_in_f = Shape::null();
                if !find_shape(a_s_or, &a_f_or, a_bf.my_analyzer.as_ref(), &mut a_s_in_f) {
                    continue;
                }
                let entry = an_images
                    .entry(shape_key_of(&a_s_in_f))
                    .or_insert((a_s_in_f.clone(), Vec::new()));
                append_to_list(&mut entry.1, &a_e);
            }
        }
    }

    // OCCT L2164-2167: the map used to find the edges on the original face
    // adjacent to the same vertex; filled at first necessity.
    let mut a_dmvef_or: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();

    // OCCT L2169-2466: the second walk — the per-split classification.
    for a_f_im in the_lf_images.iter() {
        let _a_ps = ProgressScope::new(&NoopProgressTag, "", 2);
        // OCCT L2179: valid edges for this split.
        let mut a_mve: OcctShapeSet = HashMap::new();
        // OCCT L2181: invalid edges for this split.
        let mut a_mie = OcctIndexedShapeMap::new();

        for a_e_cur in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
            let a_e_im = a_e_cur.clone();
            // OCCT L2188-2192.
            if a_e_im.orientation == Orientation::Internal {
                set_add(&mut a_me_int, &a_e_im);
                continue;
            }
            // OCCT L2194-2198.
            let p_le_or: Vec<Shape> =
                match shape_data_map::seek(&my_edges_origins, &a_e_im) {
                    Some(v) if !v.is_empty() => v.clone(),
                    _ => continue,
                };
            // OCCT L2200-2212: count the vertex-originated originals.
            let mut a_nb_v_or = 0usize;
            for a_s in p_le_or.iter() {
                if a_s.shape_type() == ShapeType::Vertex {
                    a_nb_v_or += 1;
                }
            }
            if a_nb_v_or > 1 && (p_le_or.len() - a_nb_v_or) > 1 {
                continue;
            }
            // OCCT L2214-2219.
            let mut a_me: OcctShapeSet = HashMap::new();
            let mut a_mv: OcctShapeSet = HashMap::new();
            let mut a_mf: OcctShapeSet = HashMap::new();
            let mut b_invalid = false;
            let mut b_checked = false;
            let a_nb_p = nb_points(&a_e_im);
            let mut a_nb_inv = 0usize;
            let b_use_vertex = if a_nb_v_or == 0 {
                false
            } else {
                a_nb_v_or == 1
                    && shape_indexed_data_map::find(&a_dmef, &a_e_im).len() == 1
                    && !shape_data_map::is_bound(&my_oe_origins, &a_e_im)
            };

            // OCCT L2221-2378: the per-original classification.
            for a_s_or in p_le_or.iter() {
                let b_vertex = a_s_or.shape_type() == ShapeType::Vertex;
                let mut a_e_or_f = Shape::null();
                if b_vertex {
                    // OCCT L2230-2234: for some cases it is impossible to
                    // check the validity of the edge.
                    if !b_use_vertex {
                        continue;
                    }
                    // OCCT L2235-2240: find edges on the original face
                    // adjacent to this vertex.
                    if a_dmvef_or.is_empty() {
                        map_shapes_and_ancestors_map(
                            &a_f_or,
                            ShapeType::Vertex,
                            ShapeType::Edge,
                            &mut a_dmvef_or,
                        );
                    }
                    if let Some(p_lef_or) = a_dmvef_or.get(&shape_key_of(a_s_or)) {
                        let p_lef_or = p_lef_or.1.clone();
                        // OCCT L2245-2246: aCEOr compound.
                        let mut a_ce_or = super::brep_offset_make_offset_1::empty_compound();
                        // OCCT L2247-2262: avoid classification of edges
                        // originated from vertices located between tangent
                        // edges.
                        let mut b_all_tgt = true;
                        let a_v_ref = get_average_tangent(&p_lef_or[0], a_nb_p);
                        for a_e_or in p_lef_or.iter() {
                            super::brep_offset_make_offset_1::add_to_container_shape(
                                a_e_or,
                                &mut a_ce_or,
                            );
                            let a_v_cur = get_average_tangent(a_e_or, a_nb_p);
                            if !vec_is_parallel(a_v_ref, a_v_cur, rcad_kernel::precision::ANGULAR)
                            {
                                b_all_tgt = false;
                            }
                        }
                        if !b_all_tgt {
                            a_e_or_f = a_ce_or;
                        }
                    }
                } else {
                    // OCCT L2271: FindShape(aSOr, aFOr, myAnalyzer, aEOrF).
                    find_shape(a_s_or, &a_f_or, a_bf.my_analyzer.as_ref(), &mut a_e_or_f);
                    // OCCT L2273-2278: theDMEOrLEIm connection.
                    {
                        let entry = the_dmeorleim
                            .entry(shape_key_of(a_s_or))
                            .or_insert((a_s_or.clone(), Vec::new()));
                        append_to_list(&mut entry.1, &a_e_im);
                    }
                }
                // OCCT L2281-2285: the edge has not been found.
                if a_e_or_f.is_null() {
                    continue;
                }

                if b_vertex {
                    // OCCT L2287-2331: check the original edges sharing the
                    // vertex.
                    let mut a_mv_total: OcctShapeSet = HashMap::new();
                    let mut a_nb_checked = 0usize;
                    for a_e_or in bat::sub_shapes(&a_e_or_f) {
                        let a_l_im = match shape_data_map::seek(&an_images, &a_e_or) {
                            Some(v) => v.clone(),
                            None => continue,
                        };
                        a_nb_checked += 1;
                        let mut a_mv_loc: ShapeIndexedDataMap<Vec<Shape>> =
                            indexmap::IndexMap::new();
                        for a_im in a_l_im.iter() {
                            map_shapes_and_ancestors_map(
                                a_im,
                                ShapeType::Vertex,
                                ShapeType::Edge,
                                &mut a_mv_loc,
                            );
                        }
                        for i in 1..=a_mv_loc.len() {
                            let a_key =
                                shape_indexed_data_map::find_key_1(&a_mv_loc, i).clone();
                            let a_list = shape_indexed_data_map::value_1(&a_mv_loc, i);
                            if a_list.len() > 1 && !set_add(&mut a_mv_total, &a_key) {
                                b_invalid = true;
                                set_add(the_edges_invalid_by_vertex, &a_e_im);
                                break;
                            }
                        }
                        if b_invalid {
                            break;
                        }
                    }
                    if !b_invalid && a_nb_checked < 2 {
                        continue;
                    } else {
                        set_add(the_edges_valid_by_vertex, &a_e_im);
                    }
                } else {
                    // OCCT L2333-2375: check orientations of the image edge
                    // and the original edge.
                    let mut a_v_sum1 = get_average_tangent(&a_e_im, a_nb_p);
                    let mut a_v_sum2 = get_average_tangent(&a_e_or_f, a_nb_p);
                    let l1 = a_v_sum1.length();
                    let l2 = a_v_sum2.length();
                    if l1 > 0. {
                        a_v_sum1 /= l1;
                    }
                    if l2 > 0. {
                        a_v_sum2 /= l2;
                    }
                    let a_cos = a_v_sum1.dot(a_v_sum2);
                    if a_cos.abs() < 0.9999 {
                        continue;
                    }
                    // OCCT L2354-2360.
                    set_add(&mut a_me, &a_e_or_f);
                    for a_v in explorer(&a_e_or_f, ShapeType::Vertex, ShapeType::Shape) {
                        set_add(&mut a_mv, &a_v);
                    }
                    if let Some(the_analyzer) = a_bf.my_analyzer.as_ref() {
                        // OCCT L2361-2369: myAnalyzer->Ancestors(aEOrF)
                        // (the BRepOffset_Analyse carrier of
                        // brep_offset_tool_d.rs — the not-yet-translated
                        // Analyse unit).
                        for a_it_fa in the_analyzer.ancestors(&a_e_or_f).iter() {
                            set_add(&mut a_mf, a_it_fa);
                        }
                    }
                    // OCCT L2371-2375.
                    if a_cos < rcad_kernel::precision::CONFUSION {
                        b_invalid = true;
                        a_nb_inv += 1;
                    }
                }
                // OCCT L2377: bChecked = true.
                b_checked = true;
            }

            // OCCT L2380-2383.
            if !b_checked {
                continue;
            }
            // OCCT L2385-2428.
            let mut b_local_only = a_nb_v_or > 1 && (p_le_or.len() - a_nb_v_or) > 1;
            let a_nb_e = a_me.len();
            let a_nb_v = a_mv.len();
            if a_nb_e > 1 && a_nb_v == 2 * a_nb_e {
                let mut b_skip = true;

                // OCCT L2391-2396: allow the edge to be accounted for the
                // local analysis when it is originated from more than two
                // faces / unanimously classified / not a boundary edge.
                if a_mf.len() > 2 && (a_nb_inv == 0 || a_nb_inv == a_nb_e) {
                    if the_lf_images.len() > 2 {
                        // OCCT L2401-2416: the vertex walk over aDMVE /
                        // aDMEF.
                        let mut b_broke_v = false;
                        let mut b_broke_e = false;
                        'vwalk: for a_v in bat::sub_shapes(&a_e_im) {
                            let a_dmve_list =
                                shape_indexed_data_map::find(&a_dmve, &a_v).clone();
                            for a_e2 in a_dmve_list.iter() {
                                if shape_indexed_data_map::find(&a_dmef, a_e2).len() < 2 {
                                    b_broke_e = true;
                                    break;
                                }
                            }
                            if b_broke_e {
                                b_broke_v = true;
                                break 'vwalk;
                            }
                        }
                        // OCCT L2417: bSkip = itV.More().
                        b_skip = b_broke_v;
                    }
                }
                if b_skip {
                    continue;
                } else {
                    b_local_only = true;
                }
            }
            // OCCT L2430-2439.
            if b_invalid {
                if !b_local_only {
                    a_bf.my_invalid_edges.add(&a_e_im);
                }
                a_mie.add(&a_e_im);
                a_me_inv.add(&a_e_im);
                continue;
            }
            // OCCT L2441-2442: check if the edge has been inverted.
            let b_inverted = if a_nb_e == 0 || b_local_only {
                false
            } else {
                a_bf.check_inverted(&a_e_im, &a_f_or, &a_dmve, &a_m_edges)
            };
            // OCCT L2444-2452.
            if !b_inverted || a_nb_v_or == 0 {
                if !b_local_only {
                    a_bf.my_valid_edges.add(&a_e_im);
                }
                set_add(&mut a_mve, &a_e_im);
                set_add(&mut a_me_val, &a_e_im);
            }
        }

        // OCCT L2455-2459: valid edges.
        if !a_mve.is_empty() {
            shape_data_map::bind(the_dmfmve, a_f_im, a_mve);
        }
        // OCCT L2461-2465: invalid edges.
        if !a_mie.is_empty() {
            shape_data_map::bind(the_dmfmie, a_f_im, a_mie);
        }
    }

    // OCCT L2468-2472: process invalid edges — the inverted check and the
    // neutral fill.
    let mut a_mvie: OcctShapeSet = HashMap::new();
    let mut a_mne: OcctShapeSet = HashMap::new();

    // OCCT L2474-2517.
    let a_nb_e_inv = a_me_inv.extent();
    for i in 1..=a_nb_e_inv {
        let a_e_im = a_me_inv.find_key_1(i).clone();
        // OCCT L2479-2485: neutral edges — valid for one split and invalid
        // for the other.
        if set_contains(&a_me_val, &a_e_im) {
            set_add(&mut a_mne, &a_e_im);
            continue;
        }
        // OCCT L2487-2491.
        if !a_bf.my_inverted_edges.contains(&a_e_im) {
            continue;
        }
        // OCCT L2493-2497.
        let my_oe_origins_map = &my_oe_origins;
        let p_loe_or = match shape_data_map::seek(&my_oe_origins_map, &a_e_im) {
            Some(v) => v.clone(),
            None => continue,
        };
        for a_oe_or in p_loe_or.iter() {
            // OCCT L2503: aLEIm1 = myOEImages.Find(aOEOr).
            if !shape_data_map::is_bound(&my_oe_images, a_oe_or) {
                continue;
            }
            let a_le_im1 = shape_data_map::value(&my_oe_images, a_oe_or).clone();
            for a_e_im1 in a_le_im1.iter() {
                if a_m_edges.contains(a_e_im1)
                    && !a_me_inv.contains(a_e_im1)
                    && !set_contains(&a_me_int, a_e_im1)
                    && a_bf.my_inverted_edges.contains(a_e_im1)
                {
                    a_bf.my_invalid_edges.add(a_e_im1);
                    set_add(&mut a_mvie, a_e_im1);
                }
            }
        }
    }

    // OCCT L2519-2527.
    if !a_mne.is_empty() {
        shape_data_map::bind(the_dmfmne, the_f, a_mne);
    }
    if !a_mvie.is_empty() {
        shape_data_map::bind(the_dmfmvie, the_f, a_mvie);
    }
}

// The NoopProgress alias of the progress-flattened scopes (architecture
// difference #39).
use rcad_kernel::core::message::NoopProgress as NoopProgressTag;

// ---------------------------------------------------------------------------
// OCCT addAsNeutral (cxx L2536-2561).
// ---------------------------------------------------------------------------

/// OCCT static addAsNeutral (cxx L2536-2561) — adds the edge into the
/// corresponding maps making it neutral.
pub(crate) fn add_as_neutral(
    the_e: &Shape,
    the_f_inv: &Shape,
    the_f_val: &Shape,
    the_loc_inv_edges: &mut DataMapOfShapeIndexedMapOfShape,
    the_loc_valid_edges: &mut ShapeDataMap<OcctShapeSet>,
) {
    // OCCT L2544-2551.
    {
        let entry = the_loc_inv_edges
            .entry(shape_key_of(the_f_inv))
            .or_insert((the_f_inv.clone(), OcctIndexedShapeMap::new()));
        entry.1.add(the_e);
    }
    // OCCT L2553-2560.
    {
        let entry = the_loc_valid_edges
            .entry(shape_key_of(the_f_val))
            .or_insert((the_f_val.clone(), HashMap::new()));
        set_add(&mut entry.1, the_e);
    }
}

// ---------------------------------------------------------------------------
// OCCT FindInvalidEdges, the list form (cxx L2569-2769).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FindInvalidEdges (cxx L2569-2769) —
/// additional method to look for invalid edges:
/// 1. find edges unclassified in faces;
/// 2. find SD faces in which the same edge is classified;
/// 3. check if the edge is neutral in the face in which it was not
///    classified.
pub(crate) fn find_invalid_edges_list_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_lfoffset: &[Shape],
    the_loc_inv_edges: &mut DataMapOfShapeIndexedMapOfShape,
    the_loc_valid_edges: &mut ShapeDataMap<OcctShapeSet>,
    the_neutral_edges: &mut ShapeDataMap<OcctShapeSet>,
) {
    let my_faces_origins = a_bf.my_faces_origins.clone().unwrap_or_default();

    // OCCT L2583-2586: the unclassified edges per face.
    let mut a_me_unclassified: ShapeIndexedDataMap<OcctShapeSet> = indexmap::IndexMap::new();
    // OCCT L2587: the split face -> offset face map.
    let mut a_f_split_f_offset: ShapeDataMap<Shape> = HashMap::new();

    // OCCT L2589-2608: avoid artificial faces.
    let mut a_new_faces: OcctShapeSet = HashMap::new();
    if let Some(the_analyzer) = a_bf.my_analyzer.as_ref() {
        let mut a_map_new_tmp: OcctShapeSet = HashMap::new();
        // OCCT L2594: myAnalyzer->NewFaces() (the BRepOffset_Analyse
        // carrier of brep_offset_tool_d.rs — the not-yet-translated
        // Analyse unit).
        for it in the_analyzer.new_faces().iter() {
            set_add(&mut a_map_new_tmp, it);
        }
        for a_f_offset in the_lfoffset.iter() {
            let a_f_origin = shape_data_map::find(&my_faces_origins, a_f_offset);
            if set_contains(&a_map_new_tmp, &a_f_origin) {
                set_add(&mut a_new_faces, a_f_offset);
            }
        }
    }

    // OCCT L2610-2656: collect the unclassified edges.
    let mut an_ef_map: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    for a_f in the_lfoffset.iter() {
        if set_contains(&a_new_faces, a_f) {
            continue;
        }
        // OCCT L2620: aLFImages = myOFImages.FindFromKey(aF).
        let a_lf_images = match a_bf.my_of_images.get(&shape_key_of(a_f)) {
            Some(v) => v.1.clone(),
            None => continue,
        };
        for a_f_im in a_lf_images.iter() {
            map_shapes_and_ancestors_map(a_f_im, ShapeType::Edge, ShapeType::Face, &mut an_ef_map);

            // OCCT L2627-2630.
            let p_me_invalid = shape_data_map::seek(the_loc_inv_edges, a_f_im);
            let p_me_valid = shape_data_map::seek(the_loc_valid_edges, a_f_im);

            for a_e in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                // OCCT L2635: the edge is classified in some face.
                if a_bf.my_invalid_edges.contains(&a_e) != a_bf.my_valid_edges.contains(&a_e) {
                    let in_inv = p_me_invalid.map(|m| m.contains(&a_e)).unwrap_or(false);
                    let in_val = p_me_valid.map(|m| set_contains(m, &a_e)).unwrap_or(false);
                    // OCCT L2639-2652: but not in the current one.
                    if !in_inv && !in_val {
                        let entry = a_me_unclassified
                            .entry(shape_key_of(&a_e))
                            .or_insert((a_e.clone(), HashMap::new()));
                        set_add(&mut entry.1, a_f_im);
                        shape_data_map::bind(&mut a_f_split_f_offset, a_f_im, a_f.clone());
                    }
                }
            }
        }
    }

    // OCCT L2658-2661.
    if a_me_unclassified.is_empty() {
        return;
    }

    // OCCT L2663-2768: analyze the unclassified edges.
    let a_nb_e = a_me_unclassified.len();
    for i_e in 1..=a_nb_e {
        let (a_e, a_mf_unclassified) = {
            let (_, v) = a_me_unclassified.get_index(i_e - 1).expect("index");
            (v.0.clone(), v.1.clone())
        };
        // OCCT L2671: aLF = anEFMap.FindFromKey(aE).
        let a_lf = match an_ef_map.get(&shape_key_of(&a_e)) {
            Some(v) => v.1.clone(),
            None => continue,
        };

        for a_f_classified in a_lf.iter() {
            if set_contains(&a_mf_unclassified, a_f_classified) {
                continue;
            }

            // OCCT L2681-2682: BOPTools_Set of the classified face.
            let mut an_edge_set_class = BOPToolsSet {
                my_shapes: Vec::new(),
                my_nb_shapes: 0,
            };
            an_edge_set_class.add(a_f_classified, ShapeType::Edge);

            // OCCT L2684-2686.
            let mut a_e_classified = Shape::null();
            find_shape(&a_e, a_f_classified, None, &mut a_e_classified);
            let an_ori_class = a_e_classified.orientation;

            // OCCT L2688-2691: the normal of the classified face.
            let a_dn_class = get_normal_to_face_on_edge(&a_e_classified, a_f_classified);

            // OCCT L2693-2695.
            let is_invalid = shape_data_map::seek(the_loc_inv_edges, a_f_classified)
                .map(|m| m.contains(&a_e))
                .unwrap_or(false);

            for a_f_unclassified in a_mf_unclassified.values() {
                // OCCT L2703-2704: BOPTools_Set of the unclassified face.
                let mut an_edge_set_unclass = BOPToolsSet {
                    my_shapes: Vec::new(),
                    my_nb_shapes: 0,
                };
                an_edge_set_unclass.add(a_f_unclassified, ShapeType::Edge);

                // OCCT L2706: the content equality.
                if an_edge_set_class.is_equal(&an_edge_set_unclass) {
                    // OCCT L2708-2711.
                    let a_dn_unclass =
                        get_normal_to_face_on_edge(&a_e, a_f_unclassified);
                    let is_same_ori = match (&a_dn_class, &a_dn_unclass) {
                        (Some(c), Some(u)) => dir_is_equal(
                            *c,
                            *u,
                            rcad_kernel::precision::ANGULAR,
                        ),
                        _ => false,
                    };

                    // OCCT L2715-2760: among the other splits of the same
                    // face find those where the edge is contained with a
                    // different orientation.
                    let a_f_offset =
                        shape_data_map::find(&a_f_split_f_offset, a_f_unclassified);
                    let a_lf_splits = match a_bf.my_of_images.get(&shape_key_of(&a_f_offset)) {
                        Some(v) => v.1.clone(),
                        None => continue,
                    };
                    for a_f_sp in a_lf_splits.iter() {
                        if !a_f_sp.is_same(a_f_unclassified)
                            && set_contains(&a_mf_unclassified, a_f_sp)
                        {
                            let mut a_e_unclassified = Shape::null();
                            find_shape(&a_e, a_f_sp, None, &mut a_e_unclassified);

                            let mut an_ori_unclass = a_e_unclassified.orientation;
                            if !is_same_ori {
                                an_ori_unclass = bat::top_abs_reverse(an_ori_unclass);
                            }

                            if an_ori_class != an_ori_unclass {
                                // OCCT L2737-2746: make the edge neutral
                                // for the face.
                                {
                                    let entry = the_neutral_edges
                                        .entry(shape_key_of(&a_f_offset))
                                        .or_insert((a_f_offset.clone(), HashMap::new()));
                                    set_add(&mut entry.1, &a_e);
                                }

                                if is_invalid && is_same_ori {
                                    // OCCT L2748-2752: invalid in
                                    // aFClassified, valid in aFSp.
                                    add_as_neutral(
                                        &a_e,
                                        a_f_classified,
                                        a_f_sp,
                                        the_loc_inv_edges,
                                        the_loc_valid_edges,
                                    );
                                } else {
                                    // OCCT L2753-2757: invalid in aFSp,
                                    // valid in aFClassified.
                                    add_as_neutral(
                                        &a_e,
                                        a_f_sp,
                                        a_f_classified,
                                        the_loc_inv_edges,
                                        the_loc_valid_edges,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT MakeInvertedEdgesInvalid (cxx L2775-2850).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::MakeInvertedEdgesInvalid (cxx
/// L2775-2850) — makes the inverted edges located inside a loop of invalid
/// edges invalid as well.
pub(crate) fn make_inverted_edges_invalid_impl(a_bf: &mut BRepOffsetBuildOffsetFaces, the_lfoffset: &[Shape]) {
    // OCCT L2778-2781.
    if a_bf.my_invalid_edges.is_empty() || a_bf.my_inverted_edges.is_empty() {
        return;
    }

    // OCCT L2783-2789: map all invalid edges into a compound.
    let mut a_cbe_inv = super::brep_offset_make_offset_1::empty_compound();
    for i in 1..=a_bf.my_invalid_edges.extent() {
        let a_e = a_bf.my_invalid_edges.find_key_1(i).clone();
        bat::builder_add_compound_shape(&mut a_cbe_inv, &a_e);
    }

    // OCCT L2791-2793: make loops of invalid edges — BOPTools_AlgoTools::
    // MakeConnexityBlocks(aCBEInv, TopAbs_VERTEX, TopAbs_EDGE, aLCB) (the
    // wire_splitter re-host, architecture difference #44; the shapes carry
    // the identity location of the offset pipeline).
    let a_inv_edges: Vec<Shape> = a_bf.my_invalid_edges.iter().cloned().collect();
    let a_locations = [glam::DAffine3::IDENTITY];
    let a_lcb =
        crate::bop::algo::wire_splitter::make_connexity_blocks(&a_inv_edges, &a_locations);

    // OCCT L2795-2826: analyze each loop on closeness; use only closed
    // ones.
    let mut a_dmvcv: ShapeDataMap<Shape> = HashMap::new();
    for a_block in a_lcb.iter() {
        let a_cb_edges = &a_block.shapes;

        let mut a_dmve: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
        // OCCT L2806: MapShapesAndAncestors(aCB, VERTEX, EDGE, aDMVE) over
        // the block's edges.
        for a_e in a_cb_edges.iter() {
            for a_v in bat::sub_shapes(a_e) {
                let entry = a_dmve
                    .entry(shape_key_of(&a_v))
                    .or_insert((a_v.clone(), Vec::new()));
                append_to_list(&mut entry.1, a_e);
            }
        }
        let mut is_closed = true;
        for i_v in 1..=a_dmve.len() {
            if shape_indexed_data_map::value_1(&a_dmve, i_v).len() != 2 {
                is_closed = false;
                break;
            }
        }
        if !is_closed {
            continue;
        }
        // OCCT L2822-2825: bind the loop to each vertex of the loop.
        for i_v in 1..=a_dmve.len() {
            let a_key_v = shape_key_of(shape_indexed_data_map::find_key_1(&a_dmve, i_v));
            a_dmvcv.insert(
                a_key_v,
                (
                    shape_indexed_data_map::find_key_1(&a_dmve, i_v).clone(),
                    super::brep_offset_make_offset_1::block_shape(a_cb_edges),
                ),
            );
        }
    }

    // OCCT L2828-2849: check if any inverted edges of the offset faces are
    // locked inside the loops of invalid edges.
    for a_f in the_lfoffset.iter() {
        let a_lf_im = match a_bf.my_of_images.get(&shape_key_of(a_f)) {
            Some(v) => v.1.clone(),
            None => continue,
        };
        for a_f_im in a_lf_im.iter() {
            for a_e_cur in explorer(a_f_im, ShapeType::Edge, ShapeType::Shape) {
                let a_e = a_e_cur;
                if !a_bf.my_invalid_edges.contains(&a_e) && a_bf.my_inverted_edges.contains(&a_e)
                {
                    // OCCT L2840-2841: FirstVertex/LastVertex.
                    let (v_first, v_last) = bat::top_exp_vertices_raw(&a_e);
                    let p_cb1 = v_first
                        .as_ref()
                        .and_then(|v| shape_data_map::seek(&a_dmvcv, v));
                    let p_cb2 = v_last
                        .as_ref()
                        .and_then(|v| shape_data_map::seek(&a_dmvcv, v));
                    if let (Some(p_cb1), Some(p_cb2)) = (p_cb1, p_cb2) {
                        if p_cb1.is_same(p_cb2) {
                            a_bf.my_invalid_edges.add(&a_e);
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT FindInvalidFaces (cxx L2856-3083).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_BuildOffsetFaces::FindInvalidFaces (cxx L2856-3083) —
/// looking for the invalid faces by analyzing their invalid edges.
#[allow(clippy::too_many_arguments, unused_assignments)]
pub(crate) fn find_invalid_faces_impl(
    a_bf: &mut BRepOffsetBuildOffsetFaces,
    the_lf_images: &mut Vec<Shape>,
    the_dmfmve: &ShapeDataMap<OcctShapeSet>,
    the_dmfmie: &DataMapOfShapeIndexedMapOfShape,
    the_me_neutral: &OcctShapeSet,
    the_edges_invalid_by_vertex: &OcctShapeSet,
    the_edges_valid_by_vertex: &OcctShapeSet,
    the_mf_holes: &OcctShapeSet,
    the_mf_inv_in_hole: &mut OcctIndexedShapeMap,
    the_inv_faces: &mut Vec<Shape>,
    the_inverted_faces: &mut Vec<Shape>,
) {
    // OCCT L2880-2885: the per-face state flags (the OCCT declaration
    // block).
    let mut b_has_valid: bool;
    let mut b_all_valid: bool;
    let mut b_all_invalid: bool;
    let mut b_has_really_invalid: bool;
    let mut b_all_inv_neutral: bool;
    let mut b_valid: bool;
    let mut b_valid_loc: bool;
    let mut b_invalid: bool;
    let mut b_invalid_loc: bool;
    let mut b_neutral: bool;
    let mut b_inverted: bool;
    let mut b_is_invalid_by_inverted: bool;
    let mut b_has_inverted: bool;
    let mut a_nb_checked: usize;

    // OCCT L2885.
    let b_treat_inverted_as_invalid = the_lf_images.len() == 1;
    // OCCT L2888: neutral edges to remove.
    let mut a_men_rem = OcctIndexedShapeMap::new();
    // OCCT L2891: faces for post treat.
    let mut a_lfpt: Vec<Shape> = Vec::new();

    // OCCT L2893-2900: aDMEF over all splits.
    let mut a_dmef: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
    for a_f_im in the_lf_images.iter() {
        map_shapes_and_ancestors_map(a_f_im, ShapeType::Edge, ShapeType::Face, &mut a_dmef);
    }

    // OCCT L2902-3033: the main classification walk (the iterator-removal
    // form of the OCCT code is the index loop).
    let mut i_f = 0usize;
    while i_f < the_lf_images.len() {
        let a_f_im = the_lf_images[i_f].clone();
        // OCCT L2908-2911: the per-split edge maps.
        let p_mve = shape_data_map::seek(the_dmfmve, &a_f_im);
        let p_mie = shape_data_map::seek(the_dmfmie, &a_f_im);
        //
        b_has_valid = false;
        b_all_valid = true;
        b_all_invalid = true;
        b_has_really_invalid = false;
        b_all_inv_neutral = true;
        b_is_invalid_by_inverted = true;
        b_has_inverted = false;
        a_nb_checked = 0;

        // OCCT L2922: aWIm = BRepTools::OuterWire(aFIm).
        let a_w_im = outer_wire_of(&a_f_im);
        let mut b_broke = false;
        for a_e_im in explorer(&a_w_im, ShapeType::Edge, ShapeType::Shape) {
            b_valid = a_bf.my_valid_edges.contains(&a_e_im);
            b_invalid = a_bf.my_invalid_edges.contains(&a_e_im);
            b_neutral = set_contains(the_me_neutral, &a_e_im);
            // OCCT L2932-2936: the edge has not been checked.
            if !b_valid && !b_invalid && !b_neutral {
                continue;
            }
            // OCCT L2938-2943: skip not-boundary edges originated from a
            // vertex.
            if (set_contains(the_edges_invalid_by_vertex, &a_e_im)
                || set_contains(the_edges_valid_by_vertex, &a_e_im))
                && shape_indexed_data_map::find(&a_dmef, &a_e_im).len() != 1
            {
                continue;
            }
            a_nb_checked += 1;
            // OCCT L2947-2953.
            b_invalid_loc = p_mie.map(|m| m.contains(&a_e_im)).unwrap_or(false);
            b_has_really_invalid = b_invalid
                && b_invalid_loc
                && !b_valid
                && !set_contains(the_edges_invalid_by_vertex, &a_e_im);
            if b_has_really_invalid {
                b_broke = true;
                break;
            }
            // OCCT L2955-2960.
            b_valid_loc = p_mve.map(|m| set_contains(m, &a_e_im)).unwrap_or(false);
            b_inverted = a_bf.my_inverted_edges.contains(&a_e_im);
            if !b_invalid && !b_invalid_loc && b_treat_inverted_as_invalid {
                b_invalid = b_inverted;
            }
            // OCCT L2962-2971.
            if b_valid_loc && b_neutral {
                b_has_valid = true;
            }
            b_all_valid &= b_valid_loc;
            b_all_invalid &= b_invalid || b_invalid_loc;
            b_all_inv_neutral &= b_all_invalid && b_neutral;
            b_is_invalid_by_inverted &= b_invalid_loc || b_inverted;
            b_has_inverted |= b_inverted;
        }
        let _ = b_broke;

        // OCCT L2974-2978.
        if a_nb_checked == 0 {
            i_f += 1;
            continue;
        }
        // OCCT L2980-2994.
        if !b_has_really_invalid && (b_all_inv_neutral && !b_has_valid) && a_nb_checked > 1 {
            if b_has_inverted {
                // The part seems to be filled due to overlapping of parts
                // rather than due to multi-connection of faces.
                i_f += 1;
                continue;
            }
            // OCCT L2990-2992: remove the edges from neutral; remove face.
            map_shapes_indexed(&a_f_im, ShapeType::Edge, &mut a_men_rem);
            the_lf_images.remove(i_f);
            continue;
        }
        // OCCT L2996-3006.
        if b_has_really_invalid
            || (b_all_invalid && !(b_has_valid || b_all_valid) && (!b_all_inv_neutral || a_nb_checked != 1))
        {
            the_inv_faces.push(a_f_im.clone());
            if set_contains(the_mf_holes, &a_f_im) {
                the_mf_inv_in_hole.add(&a_f_im);
            }
            i_f += 1;
            continue;
        }
        // OCCT L3008-3015.
        if set_contains(the_mf_holes, &a_f_im) {
            map_shapes_indexed(&a_f_im, ShapeType::Edge, &mut a_men_rem);
            the_lf_images.remove(i_f);
            continue;
        }
        // OCCT L3017-3021.
        if b_is_invalid_by_inverted && !(b_has_valid || b_all_valid) {
            the_inverted_faces.push(a_f_im.clone());
        }
        // OCCT L3023-3031.
        if !b_all_inv_neutral {
            a_lfpt.push(a_f_im.clone());
        } else {
            map_shapes_indexed(&a_f_im, ShapeType::Edge, &mut a_men_rem);
        }
        i_f += 1;
    }

    // OCCT L3035-3038.
    if a_lfpt.is_empty() || a_men_rem.is_empty() {
        return;
    }

    // OCCT L3040-3082: check the splits once more.
    for a_f_im in a_lfpt.iter() {
        let p_mve = shape_data_map::seek(the_dmfmve, a_f_im);
        b_has_valid = false;
        b_all_valid = true;
        b_all_invalid = true;

        let a_w_im = outer_wire_of(a_f_im);
        for a_e_im in explorer(&a_w_im, ShapeType::Edge, ShapeType::Shape) {
            b_valid = a_bf.my_valid_edges.contains(&a_e_im);
            b_invalid = a_bf.my_invalid_edges.contains(&a_e_im);
            b_neutral = set_contains(the_me_neutral, &a_e_im) && !a_men_rem.contains(&a_e_im);
            b_valid_loc = p_mve.map(|m| set_contains(m, &a_e_im)).unwrap_or(false);

            if !b_invalid && b_treat_inverted_as_invalid {
                b_invalid = a_bf.my_inverted_edges.contains(&a_e_im);
            }
            if b_valid_loc && b_neutral {
                b_has_valid = true;
            }
            b_all_valid = b_all_valid && b_valid_loc;
            b_all_invalid = b_all_invalid && b_invalid;
        }
        if b_all_invalid && !b_has_valid && !b_all_valid {
            the_inv_faces.push(a_f_im.clone());
        }
    }
}
