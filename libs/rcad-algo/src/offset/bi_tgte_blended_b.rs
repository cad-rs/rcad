// OCCT BiTgte_Blend.cxx L1285-2664 — the member-function half of the 1:1
// translation (ComputeCenters / ComputeSurfaces / ComputeShape / Intersect).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BiTgte/BiTgte_Blend.cxx
//
// This file is declared as the child module `blended_b` of
// bi_tgte_blended.rs (the #[path] form), so the private class fields keep
// the OCCT member visibility inside the class module tree. The statics, GAP
// carriers, class fields and the public API live in the parent.

use std::collections::{HashMap, HashSet};

use glam::DVec3;
use indexmap::IndexMap;
use rcad_kernel::geom::{Curve2dEval, Curve3, CurveEval, SurfaceEval};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo::topods::{BRepBuilder, Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::feat::brep_feat_builder::explorer;
use crate::feat::loc_ope_glued_shape::map_shapes_and_ancestors;
use crate::feat::loc_ope_wires_on_shape::{brep_tool_pnt, shape_key};
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_tolerance, ShapeKey,
};

use crate::offset::brep_offset_make_simple_offset::top_exp_vertices_shape;
use crate::offset::brep_offset_offset::{
    update_edge_2d, update_edge_2d_seam, update_face_surface, BRepOffsetStatus,
    GeomAbsShapeKind,
};
use crate::offset::brep_offset_offset_b::BRepOffsetOffset;
use crate::offset::bi_tgte_curve_on_edge::{BiTgteCurveOnEdge, GeomApiProjectPointOnCurve};

use super::{
    add, brep_lib_build_curves3d, brep_lib_make_edge_3d, brep_lib_make_edge_pcurve,
    brep_lib_same_parameter, brep_tools_update, curve_first_parameter, curve_last_parameter,
    find_created_edge, find_vertex, is_in_face, is_on_restriction, k_part_curve_3d, make_curve,
    make_degenerated_edge, orientation_of, oriented_shape, shape_reversed, touched, AncestorsMap,
    BRepOffsetAnalyse, BRepOffsetInter2d, BRepOffsetInter3d, BRepOffsetMakeLoops,
    BRepBuilderAPISewing, BiTgteBlend,
};

use crate::brep_algo::image::BRepAlgoImage;
use crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity;
use crate::topalgo::brep_bnd_lib::BRepBndLib;

/// OCCT myAncestors.FindFromKey(S) — the shared re-host keyed lookup.
fn ancestors_find<'a>(map: &'a AncestorsMap, s: &Shape) -> &'a [Shape] {
    map.get(&shape_key(s))
        .map(|(_, v)| v.as_slice())
        .unwrap_or(&[])
}

/// OCCT Geom_Curve::Transformed(Loc.Transformation()) — the gp_Trsf
/// application is not translated — GAP leaf (arch. diff. #28 family).
fn curve_transformed_gap(_the_c: Curve3, _the_loc: u32) -> Curve3 {
    panic!("GAP: Geom_Curve::Transformed (gp_Trsf application not translated)");
}

/// OCCT gp_Pnt::Transform(L.Transformation()) — GAP leaf.
fn point_transformed_gap(_the_p: DVec3, _the_loc: u32) -> DVec3 {
    panic!("GAP: gp_Pnt::Transform (gp_Trsf application not translated)");
}

impl BiTgteBlend {
    /// OCCT BiTgte_Blend::ComputeCenters() (cxx L1285-1758) — Computes the
    /// center lines.
    pub(crate) fn compute_centers(&mut self) {
        // ------------
        // Preanalyze.
        // ------------
        // OCCT L1290: TolAngle = 2 * asin(myTol / abs(myRadius * 0.5)).
        let tol_angle = 2.0 * (self.my_tol / (self.my_radius * 0.5).abs()).asin();
        self.my_analyse.perform(&self.my_shape, tol_angle);

        // ------------------------------------------
        // calculate faces touched by caps
        // ------------------------------------------
        // OCCT L1296-1297.
        let mut touched_by_cork: HashSet<Shape> = HashSet::new();
        touched(
            &self.my_analyse,
            &self.my_stop_faces,
            &self.my_shape,
            &mut touched_by_cork,
        );

        // -----------------------
        // init of the intersector
        // -----------------------
        // OCCT L1302-1307.
        let mut side = State::In;
        if self.my_radius < 0.0 {
            side = State::Out;
        }
        let mut inter = BRepOffsetInter3d::new(&self.my_as_des, side, self.my_tol);

        // OCCT L1309-1311: MapSBox / Done.
        let mut map_s_box: HashMap<Shape, BndBox> = HashMap::new();
        let done: HashSet<Shape> = HashSet::new();
        let _ = &done;

        // OCCT L1313-1315: BRep_Builder B; TopoDS_Compound Co — to only know
        // on which edges the tubes are made.
        let mut b = BRepBuilder::new();
        let co = b.make_compound(&mut self.my_brep, vec![]);

        // ----------------------------------------
        // Calculate Sections Face/Face + Propagation
        // ----------------------------------------
        let mut jen_rajoute = true;

        while jen_rajoute {
            jen_rajoute = false;

            let mut fini = false;

            // OCCT L1329: DataMap<Shape, Shape> EdgeTgt.
            let mut edge_tgt: HashMap<ShapeKey, Shape> = HashMap::new();

            while !fini {
                // -------------------------------------------------
                // locate in myFaces the Faces connected to myEdges.
                // -------------------------------------------------
                fini = true;
                // OCCT L1339-1368.
                let nb_edges = self.my_edges.len();
                for i in 0..nb_edges {
                    let e = self.my_edges.get_index(i).expect("index").0.clone();
                    if brep_tool_degenerated(&e) {
                        continue;
                    }

                    // OCCT L1347: L = myAncestors.FindFromKey(E).
                    let l = ancestors_find(&self.my_ancestors, &e).to_vec();
                    if l.len() == 1 {
                        // So this is a free border onwhich the ball should
                        // roll.
                        self.my_faces.insert(e.clone(), ());

                        // set in myStopFaces to not propagate the tube on
                        // free border.
                        self.my_stop_faces.insert(e);
                    } else {
                        for sh in &l {
                            if !self.my_stop_faces.contains(sh) {
                                self.my_faces.insert(sh.clone(), ());
                            }
                        }
                    }
                }
                self.my_edges.clear();

                // --------------------------------------------
                // Construction of Offsets of all faces.
                // --------------------------------------------
                // OCCT L1375-1476.
                let nb_faces = self.my_faces.len();
                for i in 0..nb_faces {
                    let as_ = self.my_faces.get_index(i).expect("index").0.clone();
                    if self.my_map_sf.contains_key(&as_) {
                        continue;
                    }

                    let mut of1 = BRepOffsetOffset::new();
                    let big_f = Shape::null();

                    if as_.shape_type() == ShapeType::Face {
                        let f = as_.clone();
                        if touched_by_cork.contains(&f) {
                            // OCCT L1391: BRepOffset_Tool::EnLargeFace(F,
                            // BigF, true) — GAP (arch. diff. #23).
                            let enlarged = brep_offset_tool_en_large_face(&f, true);
                            // OCCT L1392: OF1.Init(BigF, myRadius, EdgeTgt)
                            // — the (F, Offset, Created) Init form (the
                            // defaults OffsetOutside=true, Join=Arc).
                            of1.init_face_created(
                                &mut self.my_brep,
                                &enlarged,
                                self.my_radius,
                                &edge_tgt,
                                true,
                                false,
                            );
                        } else {
                            of1.init_face_created(&mut self.my_brep, &f, self.my_radius, &edge_tgt, true, false);
                        }
                    } else {
                        // So this is a Free Border edge on which the ball
                        // rolls.
                        of1.init_edge(&mut self.my_brep, &as_, self.my_radius);
                    }
                    let _ = big_f;

                    // ------------------------------------
                    // Increment the map of created tangents
                    // ------------------------------------
                    // OCCT L1407-1447.
                    let mut let_: Vec<Shape> = Vec::new();
                    if as_.shape_type() == ShapeType::Face {
                        // OCCT L1410: Analyse.Edges(TopoDS::Face(As), ...) —
                        // the face overload.
                        self.my_analyse.edges_on_face(
                            &as_,
                            ChFiDS_TypeOfConcavity::Tangential,
                            &mut let_,
                        );
                    }

                    for itlet in &let_ {
                        let cur = itlet;
                        if !edge_tgt.contains_key(&shape_key(cur)) {
                            // OCCT L1419-1420: OTE =
                            // TopoDS::Edge(OF1.Generated(Cur)) — the first
                            // Generated call.
                            let ote = of1.generated(cur);
                            // OCCT L1422: EdgeTgt.Bind(Cur,
                            // OF1.Generated(Cur)) — the second Generated
                            // call (kept literally).
                            let a_generated = of1.generated(cur);
                            edge_tgt.insert(shape_key(cur), a_generated);
                            let (v1, v2) = top_exp_vertices_shape(cur);
                            let (ov1, ov2) = top_exp_vertices_shape(&ote);
                            let mut le: Vec<Shape> = Vec::new();
                            if !edge_tgt.contains_key(&shape_key(&v1)) {
                                // OCCT L1426: Analyse.Edges(V1, ...) — the
                                // vertex overload.
                                self.my_analyse.edges_on_vertex(
                                    &v1,
                                    ChFiDS_TypeOfConcavity::Tangential,
                                    &mut le,
                                );
                                let la = self.my_analyse.ancestors(&v1);
                                if le.len() == la.len() {
                                    edge_tgt.insert(shape_key(&v1), ov1.clone());
                                }
                            }
                            if !edge_tgt.contains_key(&shape_key(&v2)) {
                                le.clear();
                                self.my_analyse.edges_on_vertex(
                                    &v2,
                                    ChFiDS_TypeOfConcavity::Tangential,
                                    &mut le,
                                );
                                let la = self.my_analyse.ancestors(&v2);
                                if le.len() == la.len() {
                                    edge_tgt.insert(shape_key(&v2), ov2.clone());
                                }
                            }
                        }
                    }
                    // end of map created tangent

                    // OCCT L1450-1453.
                    if of1.status() == BRepOffsetStatus::Reversed
                        || of1.status() == BRepOffsetStatus::Degenerated
                    {
                        continue;
                    }

                    // OCCT L1455: F1 = OF1.Face().
                    let f1 = of1.face();

                    // increment S D
                    // OCCT L1458-1463.
                    self.my_init_offset_face.set_root(&as_);
                    self.my_init_offset_face.bind(&as_, &f1);

                    let mut box1 = BndBox::new();
                    BRepBndLib::add(&f1, &mut box1, false);
                    map_s_box.insert(f1.clone(), box1);

                    // ---------------------------------------------
                    // intersection with all already created faces.
                    // ---------------------------------------------
                    fini = !self.intersect(&as_, &f1, &map_s_box, &of1, &mut inter);

                    if as_.shape_type() == ShapeType::Face {
                        b.add_to_compound(&mut self.my_brep, co.clone(), as_.clone());
                    }

                    self.my_map_sf.insert(as_, of1);
                }
            } // end of : while ( !Fini)

            //--------------------------------------------------------
            // so the offsets were created and intersected.
            // now the tubes are constructed.
            //--------------------------------------------------------
            // Construction of tubes on edge.
            //--------------------------------------------------------
            // OCCT L1485-1489.
            let mut ot = ChFiDS_TypeOfConcavity::Convex;
            if self.my_radius < 0.0 {
                ot = ChFiDS_TypeOfConcavity::Concave;
            }

            // OCCT L1491-1496.
            let mut map: AncestorsMap = IndexMap::new();
            map_shapes_and_ancestors(&co, ShapeType::Edge, ShapeType::Face, &mut map);
            map_shapes_and_ancestors(&co, ShapeType::Vertex, ShapeType::Edge, &mut map);

            // OCCT L1498-1573.
            for e in explorer(&co, ShapeType::Edge, ShapeType::Shape) {
                if self.my_map_sf.contains_key(&e) {
                    continue;
                }

                // OCCT L1507: Anc = Map.FindFromKey(E).
                let anc = ancestors_find(&map, &e).to_vec();
                if anc.len() == 2 {
                    // OCCT L1510: L = myAnalyse.Type(E).
                    let l = self.my_analyse.type_(&e);
                    if !l.is_empty() && l[0].type_of() == ot {
                        // OCCT L1513-1516.
                        let anc_first = anc[0].clone();
                        let anc_last = anc[1].clone();
                        let e_on1 = self.my_map_sf[&anc_first].generated(&e);
                        let e_on2 = self.my_map_sf[&anc_last].generated(&e);
                        // find if exits tangent edges in the original shape
                        let (v1f, v1l) = top_exp_vertices_shape(&e);
                        let mut tang_e: Vec<Shape> = Vec::new();
                        self.my_analyse.tangent_edges(&e, &v1f, &mut tang_e);
                        // find if the pipe on the tangent edges are soon
                        // created.
                        // OCCT L1526-1537.
                        let mut e1f = Shape::null();
                        let mut find = false;
                        for itl in &tang_e {
                            if find {
                                break;
                            }
                            if let Some(of_itl) = self.my_map_sf.get(itl) {
                                e1f = of_itl.generated(&v1f);
                                find = true;
                            }
                        }
                        // OCCT L1538-1552.
                        tang_e.clear();
                        self.my_analyse.tangent_edges(&e, &v1l, &mut tang_e);
                        let mut e1l = Shape::null();
                        find = false;
                        for itl in &tang_e {
                            if find {
                                break;
                            }
                            if let Some(of_itl) = self.my_map_sf.get(itl) {
                                e1l = of_itl.generated(&v1l);
                                find = true;
                            }
                        }
                        // OCCT L1553: OF1(E, EOn1, EOn2, myRadius, E1f, E1l)
                        // — the defaults Polynomial=false, Tol=1e-4,
                        // Conti=C1.
                        let of1 = BRepOffsetOffset::with_path_first_last(
                            &mut self.my_brep,
                            &e,
                            &e_on1,
                            &e_on2,
                            self.my_radius,
                            &e1f,
                            &e1l,
                            false,
                            1.0e-4,
                            GeomAbsShapeKind::C1,
                        );
                        let f1 = of1.face();

                        // maj S D
                        // OCCT L1557-1562.
                        self.my_init_offset_face.set_root(&e);
                        self.my_init_offset_face.bind(&e, &f1);

                        let mut box1 = BndBox::new();
                        BRepBndLib::add(&f1, &mut box1, false);
                        map_s_box.insert(f1.clone(), box1);

                        // ---------------------------------------------
                        // intersection with all already created faces.
                        // ---------------------------------------------
                        // OCCT L1567-1568.
                        let is_on_rest = self.intersect(&e, &f1, &map_s_box, &of1, &mut inter);
                        jen_rajoute = jen_rajoute || is_on_rest;

                        self.my_map_sf.insert(e, of1);
                    }
                }
            }
        } // end while JenRajoute

        // OCCT L1577-1578.
        self.my_edges.clear();
        let new_edges = inter.new_edges();
        for ne in new_edges {
            self.my_edges.insert(ne, ());
        }

        // -------------------------------------------------------------------
        // now it is necessary to limit edges on the neighbors (otherwise one
        // will go too far and will not be able to construct faces).
        // -------------------------------------------------------------------

        // Proceed with MakeLoops
        // OCCT L1586-1592.
        let mut a_dmvv: AncestorsMap = IndexMap::new();
        let mut ot = ChFiDS_TypeOfConcavity::Concave;
        if self.my_radius < 0.0 {
            ot = ChFiDS_TypeOfConcavity::Convex;
        }

        let mut lof: Vec<Shape> = Vec::new();
        // OCCT L1596-1660.
        let nb_faces = self.my_faces.len();
        for i in 0..nb_faces {
            let cur_s = self.my_faces.get_index(i).expect("index").0.clone();

            // tube on free border, it is undesirable.
            if self.my_stop_faces.contains(&cur_s) {
                continue;
            }

            if !self.my_map_sf.contains_key(&cur_s) {
                continue; // inverted or degenerated
            }

            let cur_of = self.my_map_sf[&cur_s].face();
            lof.push(cur_of.clone());

            if cur_s.shape_type() == ShapeType::Face {
                let cur_f = cur_s.clone();
                // OCCT L1617: expe(CurF.Oriented(TopAbs_FORWARD), EDGE).
                for cur_e in explorer(
                    &oriented_shape(&cur_f, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                ) {
                    // --------------------------------------------------------------
                    // set in myAsDes the edges generated by limitations of the
                    // initial square if the type is correct (The edges that
                    // will disappear are not set)
                    // --------------------------------------------------------------
                    // OCCT L1626: L = myAnalyse.Type(CurE).
                    let l = self.my_analyse.type_(&cur_e);
                    if !l.is_empty() && l[0].type_of() != ot {
                        // a priori doe s not disappear, so it is set
                        let cur_oe = self.my_map_sf[&cur_f].generated(&cur_e);
                        self.my_as_des
                            .add(&cur_of, &oriented_shape(&cur_oe, cur_e.orientation));
                    } else {
                        // OCCT L1638: Lanc = myAnalyse.Ancestors(CurE).
                        let lanc = self.my_analyse.ancestors(&cur_e).clone();
                        let lanc_first = lanc.first().expect("BiTgte_Blend: empty ancestors");
                        let lanc_last = lanc.last().expect("BiTgte_Blend: empty ancestors");
                        if !self.my_faces.contains_key(lanc_first)
                            || !self.my_faces.contains_key(lanc_last)
                            || self.my_stop_faces.contains(lanc_first)
                            || self.my_stop_faces.contains(lanc_last)
                        {
                            let cur_oe = self.my_map_sf[&cur_f].generated(&cur_e);
                            self.my_as_des
                                .add(&cur_of, &oriented_shape(&cur_oe, cur_e.orientation));
                        }
                    }
                }
                // OCCT L1650-1658: the empty EdgeInt map + aDMVV.
                let mut an_empty_map: HashMap<Shape, Vec<Shape>> = HashMap::new();
                BRepOffsetInter2d::compute(
                    &self.my_as_des,
                    &cur_of,
                    &self.my_edges,
                    self.my_tol,
                    &mut an_empty_map,
                    &mut a_dmvv,
                );
            }
        }

        // ----------------------------------------------------------------
        // It is also required to make 2D intersections with generated tubes
        // (Useful for unwinding)
        // ----------------------------------------------------------------
        // OCCT L1666-1705.
        let map_sf_keys: Vec<Shape> = self.my_map_sf.keys().cloned().collect();
        for cur_s in &map_sf_keys {
            if cur_s.shape_type() == ShapeType::Face {
                continue;
            }

            let cur_of = self.my_map_sf[cur_s].face();

            // no unwinding by tubes on free border.
            if self.my_stop_faces.contains(cur_s) {
                continue;
            }

            lof.push(cur_of.clone());

            // --------------------------------------------------------------
            // set in myAsDes the edge restrictions of the square
            // --------------------------------------------------------------
            for cur_oe in explorer(
                &oriented_shape(&cur_of, Orientation::Forward),
                ShapeType::Edge,
                ShapeType::Shape,
            ) {
                self.my_as_des.add(&cur_of, &cur_oe);
            }

            let mut an_empty_map: HashMap<Shape, Vec<Shape>> = HashMap::new();
            BRepOffsetInter2d::compute(
                &self.my_as_des,
                &cur_of,
                &self.my_edges,
                self.my_tol,
                &mut an_empty_map,
                &mut a_dmvv,
            );
        }
        //
        // fuse vertices on edges stored in AsDes
        // OCCT L1708-1709.
        let mut an_empty_image = BRepAlgoImage::new();
        BRepOffsetInter2d::fuse_vertices(&a_dmvv, &self.my_as_des, &mut an_empty_image);
        // ------------
        // unwinding
        // ------------
        // OCCT L1713-1714.
        let mut make_loops = BRepOffsetMakeLoops::new();
        make_loops.build(
            &mut lof,
            &self.my_as_des,
            &self.my_image_offset,
            &mut an_empty_image,
        );

        // ------------------------------------------------------------
        // It is possible to unwind edges at least one ancestor which of
        // is a face of the initial shape, so:
        // the edges generated by intersection tube-tube are missing
        // ------------------------------------------------------------

        // --------------------------------------------------------------
        // Currently set the unwinded surfaces in <myResult>
        // --------------------------------------------------------------
        // OCCT L1725: B.MakeCompound(TopoDS::Compound(myResult)).
        self.my_result = b.make_compound(&mut self.my_brep, vec![]);
        for cur_lof in &lof {
            if !self.my_image_offset.has_image(cur_lof) {
                continue;
            }

            let mut lim: Vec<Shape> = Vec::new();
            self.my_image_offset.last_image(cur_lof, &mut lim);
            for cur_lim in lim {
                // If a face is its own image, it is not set
                if cur_lim.is_same(cur_lof) {
                    break;
                }

                b.add_to_compound(&mut self.my_brep, self.my_result.clone(), cur_lim);
            }
        }
    }

    /// OCCT BiTgte_Blend::ComputeSurfaces() (cxx L1762-2218) — Perform the
    /// generated surfaces.
    pub(crate) fn compute_surfaces(&mut self) {
        // set in myFaces, the faces actually implied in the connection
        self.my_faces.clear();

        // construct
        // 1 - Tubes (True Fillets)
        // 2 - Spheres.

        // OCCT L1771-1773: Empty / EmptyMap.
        let empty: Vec<Shape> = Vec::new();
        let _ = &empty;
        let empty_map: HashMap<Shape, Vec<Shape>> = HashMap::new();
        let _ = &empty_map;

        // OCCT L1775-1776: GS1, GS2 / GC1, GC2 handles.
        let mut gs1: Option<rcad_kernel::geom::Surface3> = None;
        let mut gs2: Option<rcad_kernel::geom::Surface3> = None;
        let mut gc1: Option<rcad_kernel::geom::Curve3> = None;
        let mut gc2: Option<rcad_kernel::geom::Curve3> = None;

        // OCCT L1778-1779.
        let tol_angle = 2.0 * (self.my_tol / (self.my_radius * 0.5).abs()).asin();
        let center_analyse = BRepOffsetAnalyse::with_shape(&self.my_result, tol_angle);

        // -----------------------------------------------------
        // Construction of tubes in myResult
        // -----------------------------------------------------
        // OCCT L1784-1785.
        let mut b = BRepBuilder::new();
        self.my_result = b.make_compound(&mut self.my_brep, vec![]);

        // --------------------------------------------------------------------
        // Dummy: for construction of spheres:
        // Set in Co the center line, then it there are at least 3
        // center lines sharing the same vertex, Sphere on this vertex.
        // --------------------------------------------------------------------
        // OCCT L1792-1793.
        let co = b.make_compound(&mut self.my_brep, vec![]);

        // --------------------------------------------------------------------
        // Iteration on the edges lines of center
        // and their valid valid part is taken after cut and tube construction.
        // --------------------------------------------------------------------

        // OCCT L1802-2218.
        let nb_edges = self.my_edges.len();
        for i in 0..nb_edges {
            let cur_e = self.my_edges.get_index(i).expect("index").0.clone();

            // OCCT L1806: L = myAsDes->Ascendant(CurE).
            let l = self.my_as_des.ascendant(&cur_e).to_vec();
            if l.len() != 2 {
                continue;
            }

            // --------------------------------------------------------------
            // F1 and F2 = 2 parallel faces intersecting in CurE.
            // --------------------------------------------------------------
            let f1 = l[0].clone();
            let f2 = l[1].clone();

            // -----------------------------------------------------
            // find the orientation of edges of intersection
            // in the initial faces.
            // -----------------------------------------------------
            // OCCT L1822-1826.
            let ld1 = self.my_as_des.descendant(&f1).to_vec();
            let ld2 = self.my_as_des.descendant(&f2).to_vec();

            let orien1 = orientation_of(&cur_e, &f1, &ld1);
            let orien2 = orientation_of(&cur_e, &f2, &ld2);

            // ---------------------------------------------------------
            // Or1 and Or2 : the shapes generators of parallel faces
            // ---------------------------------------------------------
            // OCCT L1831-1835.
            let or1 = self.my_init_offset_face.image_from(&f1).clone();
            let or2 = self.my_init_offset_face.image_from(&f2).clone();

            self.my_faces.insert(or1.clone(), ());
            self.my_faces.insert(or2.clone(), ());

            let mut oe1 = Shape::null();
            let mut oe2 = Shape::null();
            let mut of1 = Shape::null();
            let mut of2 = Shape::null();
            // OCCT L1839: TopLoc_Location Loc — the rcad location index.
            let loc: u32 = 0;
            let mut f1_par: f64 = 0.0;
            let mut l1_par: f64 = 0.0;
            let mut f2_par: f64 = 0.0;
            let mut l2_par: f64 = 0.0;

            let mut of1_is_edge = false;

            // OCCT L1844-1855.
            if or1.shape_type() == ShapeType::Edge {
                of1_is_edge = true;
                oe1 = or1.clone();
                // OCCT L1848: GC1 = BRep_Tool::Curve(OE1, Loc, f1, l1) — the
                // rcad location travels with the Shape (the L out form is
                // not carried).
                if let Some((gc, pf, pl)) = crate::offset::brep_offset_make_simple_offset::edge_curve_of(&oe1) {
                    // OCCT L1849: GC1 = Transformed(Loc.Transformation()).
                    gc1 = Some(curve_transformed_gap(gc, loc));
                    f1_par = pf;
                    l1_par = pl;
                }
            } else if or1.shape_type() == ShapeType::Face {
                of1 = or1.clone();
                gs1 = Self::brep_tool_surface(&of1);
            }

            // ----------------------------------------------------------------
            // If a vertex is used in contact, currently nothing is done
            // and the vertexes are not managed (Intersections with sphere);
            // ----------------------------------------------------------------
            // OCCT L1861-1864.
            if of1.is_null() && oe1.is_null() {
                continue;
            }

            let mut of2_is_edge = false;

            // OCCT L1868-1879.
            if or2.shape_type() == ShapeType::Edge {
                of2_is_edge = true;
                oe2 = or2.clone();
                if let Some((gc, pf, pl)) = crate::offset::brep_offset_make_simple_offset::edge_curve_of(&oe2) {
                    gc2 = Some(curve_transformed_gap(gc, loc));
                    f2_par = pf;
                    l2_par = pl;
                }
            } else if or2.shape_type() == ShapeType::Face {
                of2 = or2.clone();
                gs2 = Self::brep_tool_surface(&of2);
            }
            // ----------------------------------------------------------------
            // If a vertex is used in contact, currently nothing is done
            // and the vertexes are not managed (Intersections with sphere);
            // ----------------------------------------------------------------
            // OCCT L1884-1887.
            if of2.is_null() && oe2.is_null() {
                continue;
            }
            let _ = (&f1_par, &l1_par, &f2_par, &l2_par);

            let mut cur_l: Vec<Shape> = Vec::new();

            // OCCT L1891-1905.
            if !self.my_image_offset.has_image(&cur_e) {
                // the tubes are not unwinded
                if of1_is_edge && of2_is_edge {
                    // if I don't have the image, possibly
                    cur_l.push(cur_e.clone()); // I'm on intersection tube-tube
                } // See comment on the call to
                else {
                    // MakeLoops
                    continue;
                }
            } else {
                self.my_image_offset.last_image(&cur_e, &mut cur_l);
            }

            // ---------------------------------------------------------------
            // CurL = List of edges descending from CurE ( = Cuts of CurE)
            // ---------------------------------------------------------------
            // OCCT L1910-2177.
            for itl in cur_l {
                let cur_cut_e = itl;

                // OCCT L1915-1916.
                let pc1 = brep_tool_curve_on_surface(&cur_cut_e, &f1);
                let pc2 = brep_tool_curve_on_surface(&cur_cut_e, &f2);
                if pc1.is_none() || pc2.is_none() {
                    // OCCT L1919-1922 (OCCT_DEBUG print): "No PCurves on
                    // Intersections : No tubes constructed".
                    continue;
                }
                let (pc1, f1_par, l1_par) = match pc1 {
                    Some((c, pf, pl)) => (c, pf, pl),
                    None => continue,
                };
                let (pc2, f2_par, l2_par) = match pc2 {
                    Some((c, pf, pl)) => (c, pf, pl),
                    None => continue,
                };

                // OCCT L1926-1938.
                let mut e1f = Shape::null();
                let mut e1l = Shape::null();
                let mut vf_on_e1 = Shape::null();
                let mut vl_on_e1 = Shape::null();
                let mut vf_on_e2 = Shape::null();
                let mut vl_on_e2 = Shape::null();
                let tang_e: Vec<Shape> = Vec::new();
                let mut map_on_v1f: HashSet<Shape> = HashSet::new();
                let mut map_on_v1l: HashSet<Shape> = HashSet::new();
                let _ = &tang_e;

                let (v1f, v1l) = top_exp_vertices_shape(&cur_cut_e);

                // find if the pipe on the tangent edges are soon created.
                // edges generated by V1f and V1l + Maj MapOnV1f/l
                e1f = find_created_edge(
                    &v1f,
                    &cur_cut_e,
                    &self.my_map_sf,
                    &mut map_on_v1f,
                    &center_analyse,
                    self.my_radius,
                    self.my_tol,
                );

                e1l = find_created_edge(
                    &v1l,
                    &cur_cut_e,
                    &self.my_map_sf,
                    &mut map_on_v1l,
                    &center_analyse,
                    self.my_radius,
                    self.my_tol,
                );

                // OCCT L1940-1989.
                let mut e1 = Shape::null();
                let mut e2 = Shape::null();
                if of1_is_edge {
                    // OCCT L1943-1944.
                    let con_e = BiTgteCurveOnEdge::with_edges(&cur_cut_e, &oe1);
                    let c = make_curve(&con_e);
                    // OCCT L1945-1946: C->Value(C->FirstParameter()) — the
                    // null handle would raise (the MakeCurve result is
                    // assumed non-null).
                    let c_ref = c.as_ref().expect("BiTgte_Blend: null MakeCurve result");
                    let p1 = c_ref.point_at(curve_first_parameter(c_ref));
                    let p2 = c_ref.point_at(curve_last_parameter(c_ref));
                    vf_on_e1 = find_vertex(p1, &map_on_v1f, self.my_tol, &mut self.my_brep);
                    if vf_on_e1.is_null() {
                        vf_on_e1 = find_vertex(p1, &map_on_v1l, self.my_tol, &mut self.my_brep);
                    }
                    vl_on_e1 = find_vertex(p2, &map_on_v1l, self.my_tol, &mut self.my_brep);
                    if vl_on_e1.is_null() {
                        vl_on_e1 = find_vertex(p2, &map_on_v1f, self.my_tol, &mut self.my_brep);
                    }
                    if p1.distance_squared(p2) < self.my_tol * self.my_tol {
                        // BRepOffset_Offset manages degenerated KPart
                        // It is REQUIRED that C should be a circle with ZERO
                        // radius
                        e1 = make_degenerated_edge(c_ref, &vf_on_e1, &mut self.my_brep);
                    } else {
                        // OCCT L1965: E1 = BRepLib_MakeEdge(C, VfOnE1,
                        // VlOnE1) — the NotDone conversion would raise.
                        let (_done, made) =
                            brep_lib_make_edge_3d(&mut self.my_brep, c_ref, &vf_on_e1, &vl_on_e1);
                        e1 = made;
                    }
                } else {
                    // OCCT L1970-1974.
                    let gs1_ref = gs1
                        .as_ref()
                        .expect("BiTgte_Blend: null GS1 in the face-contact branch");
                    let mut p2d = pc1.point_at(f1_par);
                    let p1 = gs1_ref.point_at(p2d.x, p2d.y);
                    p2d = pc1.point_at(l1_par);
                    let p2 = gs1_ref.point_at(p2d.x, p2d.y);
                    vf_on_e1 = find_vertex(p1, &map_on_v1f, self.my_tol, &mut self.my_brep);
                    vl_on_e1 = find_vertex(p2, &map_on_v1l, self.my_tol, &mut self.my_brep);
                    // OCCT L1977-1986: BRepLib_MakeEdge MKE(PC1, GS1, VfOnE1,
                    // VlOnE1, f1, l1) — the face-less pcurve form — GAP leaf
                    // (arch. diff. #29).
                    let (done, made) = brep_lib_make_edge_pcurve(
                        &mut self.my_brep,
                        &pc1,
                        gs1_ref,
                        &vf_on_e1,
                        &vl_on_e1,
                        f1_par,
                        l1_par,
                    );
                    if done {
                        e1 = made;
                    } else {
                        println!("Edge Not Done");
                        e1 = made;
                    }

                    k_part_curve_3d(&e1, &pc1, gs1_ref, &mut self.my_brep);
                }

                // OCCT L1991-2038.
                if of2_is_edge {
                    let con_e = BiTgteCurveOnEdge::with_edges(&cur_cut_e, &oe2);
                    let c = make_curve(&con_e);
                    let c_ref = c.as_ref().expect("BiTgte_Blend: null MakeCurve result");
                    let p1 = c_ref.point_at(curve_first_parameter(c_ref));
                    let p2 = c_ref.point_at(curve_last_parameter(c_ref));
                    vf_on_e2 = find_vertex(p1, &map_on_v1f, self.my_tol, &mut self.my_brep);
                    if vf_on_e2.is_null() {
                        vf_on_e2 = find_vertex(p1, &map_on_v1l, self.my_tol, &mut self.my_brep);
                    }
                    vl_on_e2 = find_vertex(p2, &map_on_v1l, self.my_tol, &mut self.my_brep);
                    if vl_on_e2.is_null() {
                        vl_on_e2 = find_vertex(p2, &map_on_v1f, self.my_tol, &mut self.my_brep);
                    }
                    if p1.distance_squared(p2) < self.my_tol * self.my_tol {
                        // BRepOffset_Offset manages degenerated KParts
                        // It is REQUIRED that C should be a circle with ZERO
                        // radius
                        e2 = make_degenerated_edge(c_ref, &vf_on_e2, &mut self.my_brep);
                    } else {
                        let (_done, made) =
                            brep_lib_make_edge_3d(&mut self.my_brep, c_ref, &vf_on_e2, &vl_on_e2);
                        e2 = made;
                    }
                } else {
                    let gs2_ref = gs2
                        .as_ref()
                        .expect("BiTgte_Blend: null GS2 in the face-contact branch");
                    let mut p2d = pc2.point_at(f2_par);
                    let p1 = gs2_ref.point_at(p2d.x, p2d.y);
                    p2d = pc2.point_at(l2_par);
                    let p2 = gs2_ref.point_at(p2d.x, p2d.y);
                    vf_on_e2 = find_vertex(p1, &map_on_v1f, self.my_tol, &mut self.my_brep);
                    vl_on_e2 = find_vertex(p2, &map_on_v1l, self.my_tol, &mut self.my_brep);
                    let (done, made) = brep_lib_make_edge_pcurve(
                        &mut self.my_brep,
                        &pc2,
                        gs2_ref,
                        &vf_on_e2,
                        &vl_on_e2,
                        f2_par,
                        l2_par,
                    );
                    if done {
                        e2 = made;
                    } else {
                        println!("edge not Done");
                        e2 = made;
                    }
                    k_part_curve_3d(&e2, &pc2, gs2_ref, &mut self.my_brep);
                }
                // OCCT L2039-2049: Increment of the Map of Created if
                // reconstruction of the Shape is required.
                if self.my_build_shape {
                    self.my_created.insert(cur_cut_e.clone(), HashMap::new());
                    if let Some(entry) = self.my_created.get_mut(&cur_cut_e) {
                        entry.insert(or1.clone(), vec![e1.clone()]);
                        entry.insert(or2.clone(), vec![e2.clone()]);
                    }
                }

                // ----------------------------------------------------------
                // try to init E1f, E1l, if not found with Analysis.
                // Should happen only if the THEORETICALLY tangent edges
                // are not actually tangent ( Cf: Approximation of lines
                // of intersection that add noise.)
                // ----------------------------------------------------------
                // OCCT L2057-2093.
                if e1f.is_null() && !vf_on_e1.is_null() && !vf_on_e2.is_null() {
                    for e in map_on_v1f.iter() {
                        if !e.is_null() {
                            let (a_vertex1, a_vertex2) = top_exp_vertices_shape(e);
                            if (a_vertex1.is_same(&vf_on_e1) && a_vertex2.is_same(&vf_on_e2))
                                || (a_vertex2.is_same(&vf_on_e1) && a_vertex1.is_same(&vf_on_e2))
                            {
                                e1f = e.clone();
                                break;
                            }
                        }
                    }
                }
                if e1l.is_null() && !vl_on_e1.is_null() && !vl_on_e2.is_null() {
                    for e in map_on_v1l.iter() {
                        if !e.is_null() {
                            let (a_vertex1, a_vertex2) = top_exp_vertices_shape(e);
                            if (a_vertex1.is_same(&vl_on_e1) && a_vertex2.is_same(&vl_on_e2))
                                || (a_vertex2.is_same(&vl_on_e1) && a_vertex1.is_same(&vl_on_e2))
                            {
                                e1l = e.clone();
                                break;
                            }
                        }
                    }
                }

                // OCCT L2095-2096.
                e1 = oriented_shape(&e1, orien1);
                e2 = oriented_shape(&e2, orien2);

                // OCCT L2098-2101.
                let an_offset = BRepOffsetOffset::with_path(
                    &mut self.my_brep,
                    &cur_cut_e,
                    &e1,
                    &e2,
                    -self.my_radius,
                    self.my_nubs,
                    self.my_tol,
                    GeomAbsShapeKind::C2,
                );
                self.my_map_sf.insert(cur_cut_e.clone(), an_offset);
                self.my_centers.insert(cur_cut_e.clone(), ());
                b.add_to_compound(&mut self.my_brep, co.clone(), cur_cut_e.clone());

                // OCCT L2103-2104.
                let tuyo = self.my_map_sf[&cur_cut_e].face();
                b.add_to_compound(&mut self.my_brep, self.my_result.clone(), tuyo.clone());

                // OCCT L2106-2176.
                if self.my_build_shape {
                    // method based ONLY on the construction of fillet:
                    // the first edge of the tube is exactly on Shape1.
                    // OCCT L2110: GeomAPI_ProjectPointOnCurve Projector.
                    let mut projector = GeomApiProjectPointOnCurve::new();
                    let exp: Vec<Shape> =
                        explorer(&tuyo, ShapeType::Edge, ShapeType::Shape);
                    let mut v1 = Shape::null();
                    let mut v2 = Shape::null();
                    let mut exp_index = 0usize;
                    if of1_is_edge {
                        // Update CutEdges.
                        // OCCT L2115-2116.
                        let e_on_f1 = &exp[exp_index];
                        let (ev1, ev2) = top_exp_vertices_shape(e_on_f1);
                        v1 = ev1;
                        v2 = ev2;

                        // OCCT L2118-2124.
                        let p1 = brep_tool_pnt(&v1).expect("null vertex point");
                        projector.init(p1, gc1.as_ref().expect("null GC1"));
                        let u1 = projector.lower_distance_parameter();

                        let p2 = brep_tool_pnt(&v2).expect("null vertex point");
                        projector.init(p2, gc1.as_ref().expect("null GC1"));
                        let u2 = projector.lower_distance_parameter();

                        // OCCT L2126-2129.
                        let or1_edge = &or1;
                        b.update_vertex_on_edge(
                            &mut self.my_brep,
                            oriented_shape(&v1, Orientation::Internal),
                            u1,
                            or1_edge.clone(),
                            self.my_tol,
                        );
                        b.update_vertex_on_edge(
                            &mut self.my_brep,
                            oriented_shape(&v2, Orientation::Internal),
                            u2,
                            or1_edge.clone(),
                            self.my_tol,
                        );

                        // OCCT L2135-2142.
                        self.my_cut_edges.entry(or1.clone()).or_default();
                        let l1 = self.my_cut_edges.get_mut(&or1).expect("just inserted");
                        l1.push(v1.clone());
                        l1.push(v2.clone());
                    }
                    if of2_is_edge {
                        // Update CutEdges.
                        // OCCT L2146-2147: exp.Next(); EOnF2 = exp.Current().
                        exp_index += 1;
                        let e_on_f2 = &exp[exp_index];
                        let (ev1, ev2) = top_exp_vertices_shape(e_on_f2);
                        v1 = ev1;
                        v2 = ev2;

                        // OCCT L2150-2156.
                        let p1 = brep_tool_pnt(&v1).expect("null vertex point");
                        projector.init(p1, gc2.as_ref().expect("null GC2"));
                        let u1 = projector.lower_distance_parameter();

                        let p2 = brep_tool_pnt(&v2).expect("null vertex point");
                        projector.init(p2, gc2.as_ref().expect("null GC2"));
                        let u2 = projector.lower_distance_parameter();

                        // OCCT L2158-2161.
                        b.update_vertex_on_edge(
                            &mut self.my_brep,
                            oriented_shape(&v1, Orientation::Internal),
                            u1,
                            or2.clone(),
                            self.my_tol,
                        );
                        b.update_vertex_on_edge(
                            &mut self.my_brep,
                            oriented_shape(&v2, Orientation::Internal),
                            u2,
                            or2.clone(),
                            self.my_tol,
                        );

                        // OCCT L2167-2174.
                        self.my_cut_edges.entry(or2.clone()).or_default();
                        let l2 = self.my_cut_edges.get_mut(&or2).expect("just inserted");
                        l2.push(v1.clone());
                        l2.push(v2.clone());
                    }
                }
            }
        }

        // ---------------------------------------------------
        // Construction of spheres,
        // if enough tubes arrive at the vertex
        // ---------------------------------------------------
        // OCCT L2184-2217.
        let mut map: AncestorsMap = IndexMap::new();
        map_shapes_and_ancestors(&co, ShapeType::Vertex, ShapeType::Edge, &mut map);

        for j in 0..map.len() {
            let (v, edges_of_v) = {
                let entry = map.get_index(j).expect("index");
                (entry.1.0.clone(), entry.1.1.clone())
            };
            if edges_of_v.len() != 3 {
                continue;
            }

            let mut loe: Vec<Shape> = Vec::new();

            for it_value in &edges_of_v {
                // OCCT L2201-2209: bool Reverse = true — the reversed branch
                // only.
                let reversed = true;
                if reversed {
                    let a_generated = self.my_map_sf[it_value].generated(&v);
                    loe.push(shape_reversed(&a_generated));
                } else {
                    let a_generated = self.my_map_sf[it_value].generated(&v);
                    loe.push(a_generated);
                }
            }

            // OCCT L2212: OFT(V, LOE, -myRadius, myNubs, myTol, GeomAbs_C2).
            let oft = BRepOffsetOffset::with_vertex(
                &mut self.my_brep,
                &v,
                &loe,
                -self.my_radius,
                self.my_nubs,
                self.my_tol,
                GeomAbsShapeKind::C2,
            );
            self.my_map_sf.insert(v.clone(), oft);
            self.my_centers.insert(v.clone(), ());

            // OCCT L2216: B.Add(myResult, OFT.Face()).
            let oft_face = self.my_map_sf[&v].face();
            b.add_to_compound(&mut self.my_brep, self.my_result.clone(), oft_face);
        }
    }

    /// OCCT BiTgte_Blend::ComputeShape() (cxx L2222-2507) — Build the
    /// resulting shape (all the faces must be computed).
    pub(crate) fn compute_shape(&mut self) {
        // Find in the initial Shapel:
        //  - untouched Faces
        //  - generated tubes
        //  - the faces neighbors of tubes that should be reconstructed
        //    preserving sharing.

        // For Debug : Visualize edges of the initial shape that should be
        // reconstructed.
        // end debug

        // OCCT L2232: DataMap<Shape, Shape> Created — declared and unused in
        // the OCCT body.
        let _created: HashMap<Shape, Shape> = HashMap::new();

        // OCCT L2234-2236.
        let empty: Vec<Shape> = Vec::new();
        let _ = &empty;
        let empty_map: HashMap<Shape, Vec<Shape>> = HashMap::new();
        let _ = &empty_map;

        // OCCT L2238.
        let mut b = BRepBuilder::new();

        // Maj of the Map of created.
        // Update edges that do not change in the resulting shape
        // i.e. invariant edges in the unwinding.
        // OCCT L2243-2291.
        let faces: Vec<Shape> = explorer(&self.my_shape, ShapeType::Face, ShapeType::Shape);
        for cur_f in &faces {
            if !self.my_faces.contains_key(cur_f) {
                continue; // so the face is not touched
            }

            // so the faces are unwinded
            if !self.my_map_sf.contains_key(cur_f) {
                continue; // inverted or degenerated
            }

            let cur_of = self.my_map_sf[cur_f].face();

            if !self.my_image_offset.has_image(&cur_of) {
                // face disappears in unwinding
                continue;
            }

            for cur_e in explorer(cur_f, ShapeType::Edge, ShapeType::Shape) {
                // OCCT L2273: CurOE = Offset.Generated(CurE).
                let cur_oe = self.my_map_sf[cur_f].generated(&cur_e);

                if !self.my_image_offset.has_image(&cur_oe) {
                    continue;
                }
                // CurOE disappears

                // OCCT L2283: ImE = Image(CurOE).First().
                let im_e = self
                    .my_image_offset
                    .image(&cur_oe)
                    .first()
                    .cloned()
                    .expect("BiTgte_Blend: empty image list");
                if im_e.is_same(&cur_oe) {
                    self.my_created.insert(cur_oe.clone(), HashMap::new());
                    if let Some(entry) = self.my_created.get_mut(&cur_oe) {
                        entry.insert(cur_f.clone(), vec![cur_e.clone()]);
                    }
                }
            }
        }

        // The connected faces are already in myResult.
        // So it is necessary to add faces:
        //    - non-touched (so not in myFaces)
        //    - issuing from the unwinding (non degenerated, non inverted,
        //      non disappeared)
        // OCCT L2297-2478.
        for cur_f in &faces {
            if !self.my_faces.contains_key(cur_f) {
                // so the face is not touched
                b.add_to_compound(&mut self.my_brep, self.my_result.clone(), cur_f.clone());
            } else {
                // so the faces are unwindeds

                if !self.my_map_sf.contains_key(cur_f) {
                    continue; // inverted or degenerated
                }

                let cur_of = self.my_map_sf[cur_f].face();

                if !self.my_image_offset.has_image(&cur_of) {
                    // face disappears in unwinding
                    continue;
                }

                // List of faces generated by a face in the unwinding
                let mut lim: Vec<Shape> = Vec::new();
                self.my_image_offset.last_image(&cur_of, &mut lim);
                for debouc_face in lim {
                    // DeboucFace = offset Face unwinded in "Debouc".

                    // OCCT L2332-2333: L + S = BRep_Tool::Surface(CurF, L).
                    let loc = cur_f.location;
                    let s = cur_f
                        .as_face()
                        .and_then(|fd| fd.surface.clone())
                        .expect("BiTgte_Blend::ComputeShape: null face surface");

                    // OCCT L2335-2337: B.MakeFace(NewF);
                    // B.UpdateFace(NewF, S, L, Tolerance(CurF)) — the empty
                    // face + the post-hoc surface update (BRep_Builder.cxx
                    // L564-578; the brep_offset_offset::update_face_surface
                    // carrier).
                    let mut new_f =
                        b.make_face(&mut self.my_brep, None, Shape::null());
                    update_face_surface(&new_f, &s, loc, brep_tool_tolerance(cur_f), &mut self.my_brep);

                    let mut map_ss: HashMap<Shape, Shape> = HashMap::new();

                    // OCCT L2341-2344: Face = DeboucFace.Oriented(FORWARD).
                    let face = oriented_shape(&debouc_face, Orientation::Forward);
                    // OCCT L2345-2377: the created-vertex bookkeeping.
                    for e in explorer(&face, ShapeType::Edge, ShapeType::Shape) {
                        let (v1, v2) = top_exp_vertices_shape(&e);
                        // OCCT L2351: myCreated.IsBound(E).
                        let oe = self
                            .my_created
                            .get(&e)
                            .and_then(|m| m.get(cur_f))
                            .map(|l| l[0].clone());
                        if let Some(oe) = oe {
                            let (ov1, ov2) = top_exp_vertices_shape(&oe);
                            if !self.my_created.contains_key(&v1) {
                                self.my_created.insert(v1.clone(), HashMap::new());
                            }
                            if !self.my_created.contains_key(&v2) {
                                self.my_created.insert(v2.clone(), HashMap::new());
                            }
                            let bound_v1 = self
                                .my_created
                                .get(&v1)
                                .map(|m| m.contains_key(cur_f))
                                .unwrap_or(false);
                            if !bound_v1 {
                                let entry = self.my_created.get_mut(&v1).expect("just inserted");
                                entry.insert(cur_f.clone(), vec![ov1.clone()]);
                            }
                            let bound_v2 = self
                                .my_created
                                .get(&v2)
                                .map(|m| m.contains_key(cur_f))
                                .unwrap_or(false);
                            if !bound_v2 {
                                let entry = self.my_created.get_mut(&v2).expect("just inserted");
                                entry.insert(cur_f.clone(), vec![ov2.clone()]);
                            }
                        }
                    }

                    // OCCT L2379-2470: the wire loop.
                    for w in explorer(&face, ShapeType::Wire, ShapeType::Shape) {
                        // OCCT L2383: expe(W.Oriented(FORWARD), EDGE).
                        let ow = b.make_wire(&mut self.my_brep);

                        for e in explorer(
                            &oriented_shape(&w, Orientation::Forward),
                            ShapeType::Edge,
                            ShapeType::Shape,
                        ) {
                            // OCCT L2391: C2d = BRep_Tool::CurveOnSurface(E,
                            // Face, f, l).
                            let Some((c2d, pf, pl)) = brep_tool_curve_on_surface(&e, &face) else {
                                // OCCT: the null pcurve would raise on the
                                // first use below.
                                panic!("BiTgte_Blend::ComputeShape: null pcurve on Face");
                            };
                            let mut oe = Shape::null();
                            if let Some(found) = map_ss.get(&e) {
                                // this is an edge of cutting
                                // OCCT L2395-2410.
                                oe = found.clone();
                                let reversed_e = shape_reversed(&e);
                                let c2d_1 = brep_tool_curve_on_surface(&reversed_e, &face);
                                if let Some((c2d_1, _, _)) = c2d_1 {
                                    if e.orientation == Orientation::Forward {
                                        update_edge_2d_seam(
                                            &oe,
                                            &c2d,
                                            &c2d_1,
                                            &new_f,
                                            brep_tool_tolerance(&e),
                                            &mut self.my_brep,
                                        );
                                    } else {
                                        update_edge_2d_seam(
                                            &oe,
                                            &c2d_1,
                                            &c2d,
                                            &new_f,
                                            brep_tool_tolerance(&e),
                                            &mut self.my_brep,
                                        );
                                    }
                                }
                                // OCCT L2410: B.Range(OE, f, l).
                                b.set_edge_range(&mut self.my_brep, oe.clone(), pf, pl);
                            } else {
                                // Is there an image in the Map of Created ?
                                // OCCT L2415-2461 — the outer else attaches
                                // to IsBound(E); the UpdateEdge below runs
                                // even for a null OE (the OCCT quirk).
                                if self.my_created.contains_key(&e) {
                                    if let Some(found) = self
                                        .my_created
                                        .get(&e)
                                        .and_then(|m| m.get(cur_f))
                                        .map(|l| l[0].clone())
                                    {
                                        oe = found;
                                    }
                                } else {
                                    // OCCT L2424: B.MakeEdge(OE) — the rcad
                                    // edge constructor binds the oriented
                                    // vertices together with the pcurve range
                                    // (arch. diff. #29).
                                    let (v1, v2) = top_exp_vertices_shape(&e);
                                    // OCCT L2427-2442: OV1.
                                    let ov1 = self.my_created.get(&v1).and_then(|m| m.get(cur_f)).map(|l| l[0].clone());
                                    let ov1 = match ov1 {
                                        Some(found) => found,
                                        None => {
                                            let made = b.add_vertex(
                                                &mut self.my_brep,
                                                DVec3::ZERO,
                                                brep_tool_tolerance(&v1),
                                            );
                                            // OCCT L2434: P2d = C2d->Value(
                                            // BRep_Tool::Parameter(V1, E,
                                            // Face)) — the (V, E) re-host of
                                            // the face-keyed parameter.
                                            let p2d = c2d.point_at(
                                                crate::feat::loc_ope_wires_on_shape_b::brep_tool_parameter(&v1, &e),
                                            );
                                            // OCCT L2435-2437: P = S->D0(X,
                                            // Y); P.Transform(L.Transformation()).
                                            let p = s.point_at(p2d.x, p2d.y);
                                            let p = point_transformed_gap(p, loc);
                                            // OCCT L2438-2441.
                                            b.update_vertex_point(
                                                &mut self.my_brep,
                                                made.clone(),
                                                p,
                                                brep_tool_tolerance(&v1),
                                            );
                                            self.my_created.insert(v1.clone(), HashMap::new());
                                            let entry =
                                                self.my_created.get_mut(&v1).expect("just inserted");
                                            entry.insert(cur_f.clone(), vec![made.clone()]);
                                            made
                                        }
                                    };
                                    // OCCT L2443-2458: OV2.
                                    let ov2 = self.my_created.get(&v2).and_then(|m| m.get(cur_f)).map(|l| l[0].clone());
                                    let ov2 = match ov2 {
                                        Some(found) => found,
                                        None => {
                                            let made = b.add_vertex(
                                                &mut self.my_brep,
                                                DVec3::ZERO,
                                                brep_tool_tolerance(&v2),
                                            );
                                            let p2d = c2d.point_at(
                                                crate::feat::loc_ope_wires_on_shape_b::brep_tool_parameter(&v2, &e),
                                            );
                                            let p = s.point_at(p2d.x, p2d.y);
                                            let p = point_transformed_gap(p, loc);
                                            b.update_vertex_point(
                                                &mut self.my_brep,
                                                made.clone(),
                                                p,
                                                brep_tool_tolerance(&v2),
                                            );
                                            self.my_created.insert(v2.clone(), HashMap::new());
                                            let entry =
                                                self.my_created.get_mut(&v2).expect("just inserted");
                                            entry.insert(cur_f.clone(), vec![made.clone()]);
                                            made
                                        }
                                    };
                                    // OCCT L2459-2460: B.Add(OE, OV1.Oriented(
                                    // V1.Orientation())); B.Add(OE,
                                    // OV2.Oriented(V2.Orientation())) — the
                                    // rcad edge constructor binds the oriented
                                    // vertices.
                                    oe = b.add_edge(
                                        &mut self.my_brep,
                                        None,
                                        oriented_shape(&ov1, v1.orientation),
                                        oriented_shape(&ov2, v2.orientation),
                                        [pf, pl],
                                    );
                                }
                                // OCCT L2462-2463.
                                update_edge_2d(
                                    &oe,
                                    &c2d,
                                    &new_f,
                                    brep_tool_tolerance(&e),
                                    &mut self.my_brep,
                                );
                                b.set_edge_range(&mut self.my_brep, oe.clone(), pf, pl);
                                // OCCT L2465: MapSS.Bind(E, OE).
                                map_ss.insert(e.clone(), oe.clone());
                            }
                            // OCCT L2467: B.Add(OW, OE.Oriented(E.Orientation())).
                            b.add_to_wire(
                                &mut self.my_brep,
                                ow.clone(),
                                oriented_shape(&oe, e.orientation),
                            );
                        }
                        // OCCT L2469: B.Add(NewF, OW.Oriented(W.Orientation())).
                        b.add_to_face(
                            &mut self.my_brep,
                            new_f.clone(),
                            oriented_shape(&ow, w.orientation),
                        );
                    }

                    // OCCT L2472-2475.
                    new_f = oriented_shape(&new_f, debouc_face.orientation);
                    brep_tools_update(&new_f);
                    b.add_to_compound(&mut self.my_brep, self.my_result.clone(), new_f);
                }
            }
        }

        // non-regarding the cause, there always remain greeb borders on this
        // Shape, so it is sewn.
        // OCCT L2481-2506.
        let mut sew = BRepBuilderAPISewing::new(self.my_tol);

        brep_lib_build_curves3d(&self.my_result);

        for a_face in explorer(&self.my_result, ShapeType::Face, ShapeType::Shape) {
            sew.add(&a_face);
        }

        sew.perform();

        // SameParameter is done in case Sew does not do it (Detect that the
        // edges are not sameparameter but does nothing.)

        let sewed_shape = sew.sewed_shape();
        if !sewed_shape.is_null() {
            for sec in explorer(&sewed_shape, ShapeType::Edge, ShapeType::Shape) {
                let tol = brep_tool_tolerance(&sec);
                brep_lib_same_parameter(&sec, tol);
            }
            self.my_result = sewed_shape;
        }
    }

    /// OCCT BiTgte_Blend::Intersect(Init, Face, MapSBox, OF1, Inter)
    /// (cxx L2511-2664) — Computes the intersections with <Face> and all
    /// the OffsetFaces stored in <myMapSF>. Returns <True> if an
    /// intersections ends on a boundary of a Face.
    pub(crate) fn intersect(
        &mut self,
        init: &Shape,
        face: &Shape,
        map_s_box: &HashMap<Shape, BndBox>,
        of1: &BRepOffsetOffset,
        inter: &mut BRepOffsetInter3d,
    ) -> bool {
        let mut jen_rajoute = false;

        // OCCT L2520: Box1 = MapSBox(Face).
        let box1 = map_s_box.get(face).expect("BiTgte_Blend::Intersect: no box for Face");

        // -----------------------------------------------
        // intersection with all already created faces.
        // -----------------------------------------------
        // OCCT L2525-2526.
        let init_shape1 = of1.initial_shape();
        let f1_sur_bord_libre =
            init_shape1.shape_type() == ShapeType::Edge && self.my_stop_faces.contains(&init_shape1);

        // OCCT L2528-2531.
        let mut done: HashSet<Shape> = HashSet::new();
        let map_sf_keys: Vec<Shape> = self.my_map_sf.keys().cloned().collect();
        for it_key in &map_sf_keys {
            let of2 = self
                .my_map_sf
                .get(it_key)
                .expect("BiTgte_Blend::Intersect: key vanished");
            let f2 = of2.face();

            // OCCT L2536-2539.
            let box2 = map_s_box.get(&f2).expect("BiTgte_Blend::Intersect: no box for F2");
            if box1.is_out_box(box2) {
                continue;
            }

            // OCCT L2541-2544.
            if inter.is_done(face, &f2) {
                continue;
            }

            // 2 tubes created on free border are not intersected.
            // OCCT L2547-2561.
            let init_shape2 = of2.initial_shape();
            let f2_sur_bord_libre = init_shape2.shape_type() == ShapeType::Edge
                && self.my_stop_faces.contains(&init_shape2);

            if f1_sur_bord_libre && f2_sur_bord_libre {
                continue;
            }

            // -------------------------------------------------------
            // Tubes are not intersected with neighbor faces.
            // -------------------------------------------------------
            // OCCT L2566-2574.
            if init.shape_type() == ShapeType::Edge {
                if it_key.shape_type() == ShapeType::Face && is_in_face(init, it_key) {
                    continue;
                }
            }

            // OCCT L2576.
            inter.face_inter(face, &f2, &self.my_init_offset_face);

            // ------------------------------------------
            // an edge of F1 or F2 has been touched ?
            // if yes, add faces in myFaces
            //   ==> JenRajoute = True
            // ------------------------------------------
            // OCCT L2583-2661.
            let mut l_int: Vec<Shape> = Vec::new();
            done.clear();
            if self.my_as_des.has_common_descendant(face, &f2, &mut l_int) {
                for cur_e in &l_int {
                    let (v1, v2) = top_exp_vertices_shape(cur_e);

                    // OCCT L2595: Done.Add(V1) — true when newly added.
                    if done.insert(v1.clone()) {
                        let mut e1 = Shape::null();
                        let mut e2 = Shape::null();
                        let is_on_r1 = is_on_restriction(&v1, cur_e, face, &mut e1);
                        let is_on_r2 = is_on_restriction(&v1, cur_e, &f2, &mut e2);
                        if is_on_r1 {
                            if !self.my_stop_faces.contains(init) {
                                add(
                                    &e1,
                                    &mut self.my_edges,
                                    init,
                                    of1,
                                    &self.my_analyse,
                                    is_on_r1 && is_on_r2,
                                );
                                jen_rajoute = true;
                            }
                        }
                        if is_on_r2 {
                            if !self.my_stop_faces.contains(it_key) {
                                let of2_ref = self
                                    .my_map_sf
                                    .get(it_key)
                                    .expect("BiTgte_Blend::Intersect: key vanished");
                                add(
                                    &e2,
                                    &mut self.my_edges,
                                    it_key,
                                    of2_ref,
                                    &self.my_analyse,
                                    is_on_r1 && is_on_r2,
                                );
                                jen_rajoute = true;
                            }
                        }
                    }

                    // OCCT L2625: Done.Add(V2).
                    if done.insert(v2.clone()) {
                        let mut e1 = Shape::null();
                        let mut e2 = Shape::null();
                        let is_on_r1 = is_on_restriction(&v2, cur_e, face, &mut e1);
                        let is_on_r2 = is_on_restriction(&v2, cur_e, &f2, &mut e2);

                        // If IsOnR1 && IsOnR2,
                        // Leave in the same tps on 2 faces, propagate only on
                        // free borders.
                        // A priori, only facet is closed.
                        if is_on_r1 {
                            if !self.my_stop_faces.contains(init) {
                                add(
                                    &e1,
                                    &mut self.my_edges,
                                    init,
                                    of1,
                                    &self.my_analyse,
                                    is_on_r1 && is_on_r2,
                                );
                                jen_rajoute = true;
                            }
                        }
                        if is_on_r2 {
                            if !self.my_stop_faces.contains(it_key) {
                                let of2_ref = self
                                    .my_map_sf
                                    .get(it_key)
                                    .expect("BiTgte_Blend::Intersect: key vanished");
                                add(
                                    &e2,
                                    &mut self.my_edges,
                                    it_key,
                                    of2_ref,
                                    &self.my_analyse,
                                    is_on_r1 && is_on_r2,
                                );
                                jen_rajoute = true;
                            }
                        }
                    }
                }
            }
        }

        jen_rajoute
    }
}

/// OCCT BRep_Builder::UpdateFace(F, S, L, Tol) — the shared real carrier
/// lives in brep_offset_offset::update_face_surface (BRep_Builder.cxx
/// L564-578).

/// OCCT BRepOffset_Tool::EnLargeFace(F, BigF, AddToShape) — GAP static
/// (arch. diff. #23); the out-param BigF becomes the return value.
fn brep_offset_tool_en_large_face(_the_f: &Shape, _add_to_shape: bool) -> Shape {
    panic!("GAP: BRepOffset_Tool::EnLargeFace (TKOffset/BRepOffset not translated)");
}
