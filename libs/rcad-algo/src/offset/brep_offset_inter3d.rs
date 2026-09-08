// OCCT BRepOffset_Inter3d.cxx L1-1554 + BRepOffset_Inter3d.hxx L20-134 —
// 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_Inter3d.cxx / .hxx
//
// OCCT inheritance chain (hxx L41): none — BRepOffset_Inter3d is a
// standalone class.
//
// Architecture differences (numbering continues from brep_offset_tool.rs,
// which carries the shared shape-map forms):
// 32. Handle(BRepAlgo_AsDes) -> the owned BRepAlgoAsDes (the Handle
//     refcount form belongs to the MakeOffset batch; the accessor pair
//     as_des()/as_des_mut() keeps the OCCT myAsDes access surface).
// 33. BRepOffset_Tool::Inter3D / PipeInter / FindCommonShapes /
//     OrientSection / TryProject / CheckBounds / EnLargeFace -> the
//     crate::offset::brep_offset_tool cluster (the OCCT include form).
// 34. BRepOffset_Analyse -> the BRepOffsetAnalyse GAP carrier of
//     brep_offset_tool.rs (staged as its own translation unit; the
//     accessor surface is the one consumed here).
// 35. BRepOffset_Offset (MapSF values) -> the BRepOffsetOffset of
//     brep_offset_offset_b.rs.
// 36. The CompletInt BVH pair selection -> the crate::bop::tools::box_tree
//     BoxTree + PairSelector (the OCCT-aligned BOPTools_BoxTree /
//     BOPTools_BoxPairSelector translation); Bnd_Tools::Bnd2BVH ->
//     the Aabb form built from the face bounding box (the BRepBndLib::Add
//     reduced sampling re-host below).
// 37. Message_ProgressScope / Message_ProgressRange -> the rcad
//     NoopProgress + ProgressScope forms; the OCCT child-scope chaining is
//     flattened to local scopes (the progress machinery is
//     result-inert; annotated at the sites).

use std::collections::HashMap;

use glam::DVec3;
use rcad_kernel::topo::topods::{Orientation, ShapeType, State};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::feat::brep_feat_builder::explorer;

use rcad_kernel::geom::{CurveEval, SurfaceEval};

use super::brep_offset_offset_b::BRepOffsetOffset;
use super::brep_offset_tool::{
    empty_copied, find_common_shapes, oriented, pipe_inter, set_add, set_contains,
    shape_data_map, shape_indexed_data_map, top_exp_vertices, OcctIndexedShapeMap, OcctShapeSet,
    ShapeDataMap, ShapeIndexedDataMap,
};
use super::brep_offset_tool_b::{b_update_vertex_on_edge, inter3d, try_project};
use super::brep_offset_tool_c::{check_bounds, en_large_face};
use super::brep_offset_tool_d::BRepOffsetAnalyse;

use crate::brep_algo::tool as bat;
use crate::feat::loc_ope_wires_on_shape_b::{
    brep_tool_curve, brep_tool_curve_on_surface, brep_tool_degenerated, brep_tool_range,
    brep_tool_tolerance,
};

// ---------------------------------------------------------------------------
// OCCT class (BRepOffset_Inter3d.hxx L41-134).
// ---------------------------------------------------------------------------

/// OCCT BRepOffset_Inter3d (BRepOffset_Inter3d.hxx L41-134) — computes the
/// connection of the offset and not offset faces according to the
/// connection type required; stores the result in the AsDes tool.
pub struct BRepOffsetInter3d {
    my_as_des: BRepAlgoAsDes,        // OCCT: myAsDes (hxx L128)
    my_touched: OcctIndexedShapeMap, // OCCT: myTouched (hxx L129)
    // OCCT: myDone (hxx L130) — DataMap<Face, List<Face>>.
    my_done: ShapeDataMap<Vec<Shape>>,
    my_new_edges: OcctIndexedShapeMap, // OCCT: myNewEdges (hxx L131)
    my_side: State,                    // OCCT: mySide (hxx L132)
    my_tol: f64,                       // OCCT: myTol (hxx L133)
}

impl BRepOffsetInter3d {
    /// OCCT BRepOffset_Inter3d::BRepOffset_Inter3d(AsDes, Side, Tol) (cxx
    /// L52-59).
    pub fn new(as_des: BRepAlgoAsDes, side: State, tol: f64) -> Self {
        BRepOffsetInter3d {
            my_as_des: as_des,
            my_touched: OcctIndexedShapeMap::new(),
            my_done: HashMap::new(),
            my_new_edges: OcctIndexedShapeMap::new(),
            my_side: side,
            my_tol: tol,
        }
    }

    /// OCCT BRepOffset_Inter3d::TouchedFaces() (hxx L109-112).
    pub fn touched_faces(&mut self) -> &mut OcctIndexedShapeMap {
        &mut self.my_touched
    }

    /// OCCT BRepOffset_Inter3d::AsDes() (hxx L115).
    pub fn as_des(&self) -> &BRepAlgoAsDes {
        &self.my_as_des
    }

    /// OCCT BRepOffset_Inter3d::AsDes() — the mutable form (the Handle
    /// dereference; architecture difference #32).
    pub fn as_des_mut(&mut self) -> &mut BRepAlgoAsDes {
        &mut self.my_as_des
    }

    /// OCCT BRepOffset_Inter3d::NewEdges() (hxx L118).
    pub fn new_edges(&mut self) -> &mut OcctIndexedShapeMap {
        &mut self.my_new_edges
    }

    // -----------------------------------------------------------------------
    // OCCT static ExtentEdge (cxx L63-87) — the Inter3d-local (F, E, NE)
    // form.
    // -----------------------------------------------------------------------

    /// OCCT static ExtentEdge(F, E, NE) (cxx L63-87) — the empty-copied
    /// edge extended by two fresh extremity vertices at the +/-100*length
    /// range positions.
    fn extent_edge(_f: &Shape, e: &Shape, ne: &mut Shape) {
        // OCCT L65-66: aLocalShape = E.EmptyCopied(); NE = TopoDS::Edge(..).
        *ne = empty_copied(e);

        // Enough for analytic edges, in general case reconstruct the
        // geometry of the edge recalculating the intersection of surfaces.

        // OCCT L72: NE.Orientation(TopAbs_FORWARD).
        ne.orientation = Orientation::Forward;
        // OCCT L73-77: the range extension (f/l -=/+= 100 * length).
        let (mut f_par, mut l_par) = brep_tool_range(e);
        let length = l_par - f_par;
        f_par -= 100.0 * length;
        l_par += 100.0 * length;

        // OCCT L79-80: BRep_Builder B; B.Range(NE, f, l).
        bat::builder_range_edge(ne, f_par, l_par);
        // OCCT L81: BRepAdaptor_Curve CE(E).
        let ce = brep_tool_curve(e);
        // OCCT L82-83: V1/V2 = BRepLib_MakeVertex(CE.Value(f/l)).
        let v1 = match &ce {
            Some((c, _, _)) => {
                let p = CurveEval::point_at(c, f_par);
                let mut v = bat::builder_make_vertex();
                bat::builder_update_vertex_point_tol(&mut v, p, 0.0);
                v
            }
            None => Shape::null(),
        };
        let v2 = match &ce {
            Some((c, _, _)) => {
                let p = CurveEval::point_at(c, l_par);
                let mut v = bat::builder_make_vertex();
                bat::builder_update_vertex_point_tol(&mut v, p, 0.0);
                v
            }
            None => Shape::null(),
        };
        // OCCT L84-85: B.Add(NE, V1.Oriented(FORWARD));
        // B.Add(NE, V2.Oriented(REVERSED)).
        crate::brep_algo::tool::builder_add_edge_vertex(ne, &oriented(&v1, Orientation::Forward));
        crate::brep_algo::tool::builder_add_edge_vertex(ne, &oriented(&v2, Orientation::Reversed));
        // OCCT L86: NE.Orientation(E.Orientation()).
        ne.orientation = e.orientation;
    }

    // -----------------------------------------------------------------------
    // OCCT BRepOffset_Inter3d::CompletInt (hxx L53-55; cxx L91-148).
    // -----------------------------------------------------------------------

    /// OCCT BRepOffset_Inter3d::CompletInt(SetOfFaces, InitOffsetFace,
    /// theRange) (cxx L91-148) — the complete intersection of the offset
    /// faces (the BVH pair selection).
    pub fn complet_int(
        &mut self,
        set_of_faces: &[Shape],
        init_offset_face: &BRepAlgoImage,
        _the_range: &rcad_kernel::message::ProgressScope,
    ) {
        //---------------------------------------------------------------
        // Calculate the intersections of offset faces
        // Distinction of intersection between faces // tangents.
        //---------------------------------------------------------------

        // Prepare tools for sorting the bounding boxes
        // OCCT L101-102: BOPTools_BoxTree aBBTree; SetSize(Extent).
        let mut a_bb_tree = crate::bop::tools::box_tree::BoxTree::new();
        a_bb_tree.set_size(set_of_faces.len());
        //
        // OCCT L104: aMFaces = the IndexedDataMap<Face, Bnd_Box>.
        let mut a_mfaces: ShapeIndexedDataMap<FaceBndBox> = indexmap::IndexMap::new();
        // Construct bounding boxes for faces and add them to the tree
        // OCCT L106-118.
        for it_value in set_of_faces {
            let a_f = it_value;
            //
            // compute bounding box
            // OCCT L112-113: Bnd_Box aBoxF; BRepBndLib::Add(aF, aBoxF).
            let a_box_f = brep_bnd_lib_add_face(a_f);
            //
            // OCCT L115: int i = aMFaces.Add(aF, aBoxF).
            let i = shape_indexed_data_map::add(&mut a_mfaces, a_f, a_box_f.clone());
            //
            // OCCT L117: aBBTree.Add(i, Bnd_Tools::Bnd2BVH(aBoxF)).
            a_bb_tree.add(i, a_box_f.to_aabb());
        }

        // Build BVH
        // OCCT L121: aBBTree.Build().
        a_bb_tree.build();

        // Perform selection of the pairs
        // OCCT L124-128: BOPTools_BoxPairSelector aSelector; SetBVHSets;
        // SetSame(true); Select(); Sort().
        let mut a_selector = crate::bop::tools::box_tree::PairSelector::new();
        a_selector.set_bvh_sets(&a_bb_tree);
        a_selector.set_same(true);
        a_selector.select();
        a_selector.sort();

        // Treat the selected pairs
        // OCCT L131-147: the pair walk (the std::min/max ordered faces).
        let a_pairs = a_selector.pairs().to_vec();
        let a_nb_pairs = a_pairs.len();
        let a_prog = rcad_kernel::message::NoopProgress;
        let _a_ps = rcad_kernel::message::ProgressScope::new(
            &a_prog,
            "Complete intersection",
            a_nb_pairs,
        );
        for i_pair in 0..a_nb_pairs {
            let a_pair = a_pairs[i_pair];

            let id1 = a_pair.0;
            let id2 = a_pair.1;
            // OCCT L142-143: the min/max key forms.
            let a_f1 = shape_indexed_data_map::find_key_1(&a_mfaces, id1.min(id2) + 1).clone();
            let a_f2 = shape_indexed_data_map::find_key_1(&a_mfaces, id1.max(id2) + 1).clone();

            // intersect faces
            // OCCT L146: FaceInter(aF1, aF2, InitOffsetFace).
            self.face_inter(&a_f1, &a_f2, init_offset_face);
        }
    }

    // -----------------------------------------------------------------------
    // OCCT BRepOffset_Inter3d::FaceInter (cxx L155-254).
    // -----------------------------------------------------------------------

    /// OCCT BRepOffset_Inter3d::FaceInter(F1, F2, InitOffsetFace) (cxx
    /// L155-254) — performs the intersection of the given faces.
    pub fn face_inter(&mut self, f1: &Shape, f2: &Shape, init_offset_face: &BRepAlgoImage) {
        let mut l_int1: Vec<Shape> = Vec::new();
        let mut l_int2: Vec<Shape> = Vec::new();
        let null_edge = Shape::null();
        let null_face = Shape::null();

        if f1.is_same(f2) {
            return;
        }
        if self.is_done(f1, f2) {
            return;
        }

        // OCCT L172-177: InitF1/InitF2 = InitOffsetFace.ImageFrom.
        let init_f1 = init_offset_face.image_from(f1);
        let init_f2 = init_offset_face.image_from(f2);
        if init_f1.is_same(&init_f2) {
            return;
        }

        // OCCT L179-183: the InterPipes / InterFaces probes; LE, LV;
        // LInt1.Clear(); LInt2.Clear().
        let inter_pipes = init_f2.shape_type() == ShapeType::Edge
            && init_f1.shape_type() == ShapeType::Edge;
        let inter_faces = init_f1.shape_type() == ShapeType::Face
            && init_f2.shape_type() == ShapeType::Face;
        let mut le: Vec<Shape> = Vec::new();
        let mut lv: Vec<Shape> = Vec::new();
        l_int1.clear();
        l_int2.clear();
        // OCCT L184-252: the shared-shape decision tree.
        let has_common = find_common_shapes(f1, f2, &mut le, &mut lv)
            || self.my_as_des.has_common_descendant(f1, f2, &mut le);
        if has_common {
            //-------------------------------------------------
            // F1 and F2 share shapes.
            //-------------------------------------------------
            if le.is_empty() && !lv.is_empty() {
                if inter_pipes {
                    //----------------------
                    // tubes share a vertex.
                    //----------------------
                    // OCCT L196-215: the VE1/VE2 vertex scan; the V shared
                    // vertex; the no-sphere PipeInter.
                    let ee1 = &init_f1;
                    let ee2 = &init_f2;
                    let (ve1_0, ve1_1) = top_exp_vertices(ee1);
                    let (ve2_0, ve2_1) = top_exp_vertices(ee2);
                    let ve1 = [ve1_0, ve1_1];
                    let ve2 = [ve2_0, ve2_1];
                    let mut v = Shape::null();
                    for i in 0..2 {
                        for j in 0..2 {
                            if ve1[i].is_same(&ve2[j]) {
                                v = ve1[i].clone();
                            }
                        }
                    }
                    if !init_offset_face.has_image(&v) {
                        // no sphere
                        pipe_inter(f1, f2, &mut l_int1, &mut l_int2, self.my_side);
                    }
                } else {
                    //--------------------------------------------------------
                    // Intersection having only common vertices
                    // and supports having common edges.
                    // UNSUFFICIENT, but a larger criterion shakes too
                    // many sections.
                    //--------------------------------------------------------
                    if inter_faces {
                        let init_f1_face = init_f1.clone();
                        let init_f2_face = init_f2.clone();
                        if find_common_shapes(&init_f1_face, &init_f2_face, &mut le, &mut lv) {
                            if !le.is_empty() {
                                inter3d(
                                    f1,
                                    f2,
                                    &mut l_int1,
                                    &mut l_int2,
                                    self.my_side,
                                    &null_edge,
                                    &null_face,
                                    &null_face,
                                );
                            }
                        } else {
                            inter3d(
                                f1,
                                f2,
                                &mut l_int1,
                                &mut l_int2,
                                self.my_side,
                                &null_edge,
                                &null_face,
                                &null_face,
                            );
                        }
                    }
                }
            }
        } else {
            if inter_pipes {
                pipe_inter(f1, f2, &mut l_int1, &mut l_int2, self.my_side);
            } else {
                inter3d(
                    f1,
                    f2,
                    &mut l_int1,
                    &mut l_int2,
                    self.my_side,
                    &null_edge,
                    &null_face,
                    &null_face,
                );
            }
        }
        // OCCT L253: Store(F1, F2, LInt1, LInt2).
        self.store(f1, f2, &l_int1, &l_int2);
    }

    // -----------------------------------------------------------------------
    // OCCT BRepOffset_Inter3d::ConnexIntByArc (cxx L258-442).
    // -----------------------------------------------------------------------

    /// OCCT BRepOffset_Inter3d::ConnexIntByArc(SetOfFaces, ShapeInit,
    /// Analyse, InitOffsetFace, theRange) (cxx L258-442) — computes the
    /// connections of the offset faces that have to be connected by arcs.
    pub fn connex_int_by_arc(
        &mut self,
        _set_of_faces: &[Shape],
        shape_init: &Shape,
        analyse: &BRepOffsetAnalyse,
        init_offset_face: &BRepAlgoImage,
        _the_range: &rcad_kernel::message::ProgressScope,
    ) {
        // OCCT L264-268: OT = Concave (Convex for TopAbs_OUT).
        let ot = if self.my_side == State::Out {
            crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity::Convex
        } else {
            crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity::Concave
        };
        let mut l_int1: Vec<Shape> = Vec::new();
        let mut l_int2: Vec<Shape> = Vec::new();
        let mut f1 = Shape::null();
        let mut f2 = Shape::null();
        let null_edge = Shape::null();
        let null_face = Shape::null();
        //---------------------------------------------------------------------
        // etape 1 : Intersection of faces // corresponding to the initial faces
        //           separated by a concave edge if offset > 0, otherwise convex.
        //---------------------------------------------------------------------
        // OCCT L269-308: the edge explorer walk.
        for exp_value in explorer(shape_init, ShapeType::Edge, ShapeType::Shape) {
            let e = exp_value;
            let l = analyse.type_(&e);
            if !l.is_empty() && l[0].my_type == ot {
                //-----------------------------------------------------------
                // edge is of the proper type , return adjacent faces.
                //-----------------------------------------------------------
                let anc = analyse.ancestors(&e);
                if anc.len() == 2 {
                    let init_f1 = anc[0].clone();
                    let init_f2 = anc[anc.len() - 1].clone();
                    f1 = init_offset_face.image(&init_f1).first().cloned().unwrap_or_else(Shape::null);
                    f2 = init_offset_face.image(&init_f2).first().cloned().unwrap_or_else(Shape::null);
                    if !self.is_done(&f1, &f2) {
                        inter3d(&f1, &f2, &mut l_int1, &mut l_int2, self.my_side, &e, &init_f1, &init_f2);
                        self.store(&f1, &f2, &l_int1, &l_int2);
                    }
                }
            }
        }
        //---------------------------------------------------------------------
        // etape 2 : Intersections of tubes sharing a vertex without sphere with:
        //           - tubes on each other edge sharing the vertex
        //           - faces containing an edge connected to vertex that has no tubes.
        //---------------------------------------------------------------------
        // OCCT L314-441: the tube/tube walk.
        for exp_value in explorer(shape_init, ShapeType::Edge, ShapeType::Shape) {
            let e1 = exp_value;
            if init_offset_face.has_image(&e1) {
                //---------------------------
                // E1 generated a tube.
                //---------------------------
                f1 = init_offset_face.image(&e1).first().cloned().unwrap_or_else(Shape::null);
                let (v0, v1) = top_exp_vertices(&e1);
                let v = [v0, v1];
                let anc_e1 = analyse.ancestors(&e1);

                for i in 0..2 {
                    if !init_offset_face.has_image(&v[i]) {
                        //-----------------------------
                        // the vertex has no sphere.
                        //-----------------------------
                        let anc = analyse.ancestors(&v[i]);
                        let mut tang_on_v: Vec<Shape> = Vec::new();
                        analyse.tangent_edges(&e1, &v[i], &mut tang_on_v);
                        let mut mtev: OcctShapeSet = HashMap::new();
                        for it_value in &tang_on_v {
                            set_add(&mut mtev, it_value);
                        }
                        for it_value in &anc {
                            let e2 = it_value;
                            //  Modified by skv - Fri Jan 16 16:27:54 2004 OCC4455
                            let mut is_to_skip = false;

                            if !e1.is_same(e2) {
                                let a_l = analyse.type_(e2);

                                is_to_skip = set_contains(&mtev, e2)
                                    && (a_l.is_empty()
                                        || (!a_l.is_empty() && a_l[0].my_type != ot));
                            }

                            if e1.is_same(e2) || is_to_skip {
                                continue;
                            }
                            //  Modified by skv - Fri Jan 16 16:27:54 2004 OCC4455 End
                            if init_offset_face.has_image(e2) {
                                //-----------------------------
                                // E2 generated a tube.
                                //-----------------------------
                                f2 = init_offset_face.image(e2).first().cloned().unwrap_or_else(Shape::null);
                                if !self.is_done(&f1, &f2) {
                                    //---------------------------------------------------------------------
                                    // Intersection tube/tube if the edges are not tangent (AFINIR).
                                    //----------------------------------------------------------------------
                                    pipe_inter(&f1, &f2, &mut l_int1, &mut l_int2, self.my_side);
                                    self.store(&f1, &f2, &l_int1, &l_int2);
                                }
                            } else {
                                //-------------------------------------------------------
                                // Intersection of the tube of E1 with faces //
                                // to face containing E2 if they are not tangent
                                // to the tube or if E2 is not a tangent edge.
                                //-------------------------------------------------------
                                let l = analyse.type_(e2);
                                if !l.is_empty()
                                    && l[0].my_type
                                        == crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity::Tangential
                                {
                                    continue;
                                }
                                let anc_e2 = analyse.ancestors(e2);
                                if anc_e2.len() == 2 {
                                    let mut init_f2 = anc_e2[0].clone();
                                    let mut tangent_faces = init_f2.is_same(&anc_e1[0])
                                        || init_f2.is_same(&anc_e1[anc_e1.len() - 1]);
                                    if !tangent_faces {
                                        f2 = init_offset_face.image(&init_f2).first().cloned().unwrap_or_else(Shape::null);
                                        if !self.is_done(&f1, &f2) {
                                            inter3d(
                                                &f1,
                                                &f2,
                                                &mut l_int1,
                                                &mut l_int2,
                                                self.my_side,
                                                &null_edge,
                                                &null_face,
                                                &null_face,
                                            );
                                            self.store(&f1, &f2, &l_int1, &l_int2);
                                        }
                                    }
                                    init_f2 = anc_e2[anc_e2.len() - 1].clone();
                                    tangent_faces = init_f2.is_same(&anc_e1[0])
                                        || init_f2.is_same(&anc_e1[anc_e1.len() - 1]);
                                    if !tangent_faces {
                                        f2 = init_offset_face.image(&init_f2).first().cloned().unwrap_or_else(Shape::null);
                                        if !self.is_done(&f1, &f2) {
                                            inter3d(
                                                &f1,
                                                &f2,
                                                &mut l_int1,
                                                &mut l_int2,
                                                self.my_side,
                                                &null_edge,
                                                &null_face,
                                                &null_face,
                                            );
                                            self.store(&f1, &f2, &l_int1, &l_int2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // OCCT BRepOffset_Inter3d::ConnexIntByInt (cxx L446-1041).
    // -----------------------------------------------------------------------

    /// OCCT BRepOffset_Inter3d::ConnexIntByInt(SI, MapSF, Analyse, MES,
    /// Build, Failed, theRange, bIsPlanar) (cxx L446-1041) — computes the
    /// intersection of the offset faces that have to be connected by sharp
    /// edges.
    #[allow(clippy::too_many_arguments)]
    pub fn connex_int_by_int(
        &mut self,
        si: &Shape,
        map_sf: &ShapeDataMap<BRepOffsetOffset>,
        analyse: &BRepOffsetAnalyse,
        mes: &mut ShapeDataMap<Shape>,
        build: &mut ShapeDataMap<Shape>,
        failed: &mut Vec<Shape>,
        _the_range: &rcad_kernel::message::ProgressScope,
        b_is_planar: bool,
    ) {
        let mut vemak = OcctIndexedShapeMap::new();
        let mut f1 = Shape::null();
        let mut f2 = Shape::null();
        let mut of1 = Shape::null();
        let mut of2 = Shape::null();
        let mut nf1 = Shape::null();
        let mut nf2 = Shape::null();
        let mut cur_side = self.my_side;
        let mut b_edge;
        //
        // OCCT L464: TopExp::MapShapes(SI, TopAbs_EDGE, VEmap).
        for e in explorer(si, ShapeType::Edge, ShapeType::Shape) {
            vemak.add(&e);
        }
        // Take the vertices for treatment
        // OCCT L466-488: the bIsPlanar pre-pass.
        let mut a_nb_planar_start = 0usize;
        if b_is_planar {
            let a_nb = vemak.extent();
            a_nb_planar_start = a_nb;
            for i in 1..=a_nb {
                let a_e = vemak.at_1(i).clone();
                let a_f_gen = analyse.generated(&a_e);
                if !a_f_gen.is_null() {
                    for e in explorer(&a_f_gen, ShapeType::Edge, ShapeType::Shape) {
                        vemak.add(&e);
                    }
                }
            }

            // Add vertices for treatment
            // OCCT L481: TopExp::MapShapes(SI, TopAbs_VERTEX, VEmap).
            for v in explorer(si, ShapeType::Vertex, ShapeType::Shape) {
                vemak.add(&v);
            }

            // OCCT L483-487: the NewFaces vertex adds.
            for it_nf in analyse.new_faces() {
                for v in explorer(&it_nf, ShapeType::Vertex, ShapeType::Shape) {
                    vemak.add(&v);
                }
            }
        }
        //
        // OCCT L490-493: aDMVLF1 / aDMVLF2 / aDMIntFF (DataMaps) and aDMIntE
        // (IndexedDataMap).
        let mut a_dm_vlf1: ShapeDataMap<Vec<Shape>> = HashMap::new();
        let mut a_dm_vlf2: ShapeDataMap<Vec<Shape>> = HashMap::new();
        let mut a_dm_int_ff: ShapeDataMap<Vec<Shape>> = HashMap::new();
        let mut a_dm_int_e: ShapeIndexedDataMap<Vec<Shape>> = indexmap::IndexMap::new();
        //
        // OCCT L495-649: the bIsPlanar vertex-connexity analysis.
        if b_is_planar {
            // Find internal edges in the faces to skip them while preparing faces
            // for intersection through vertices
            // OCCT L499-522: aDMFEI.
            let mut a_dm_fei: ShapeDataMap<OcctShapeSet> = HashMap::new();
            {
                for a_fx in explorer(si, ShapeType::Face, ShapeType::Shape) {
                    let mut a_mei: OcctShapeSet = HashMap::new();
                    for a_ex in explorer(&a_fx, ShapeType::Edge, ShapeType::Shape) {
                        if a_ex.orientation != Orientation::Forward
                            && a_ex.orientation != Orientation::Reversed
                        {
                            set_add(&mut a_mei, &a_ex);
                        }
                    }
                    if !a_mei.is_empty() {
                        shape_data_map::bind(&mut a_dm_fei, &a_fx, a_mei);
                    }
                }
            }

            // Analyze faces connected through vertices
            // OCCT L525: for (i = aNb + 1, aNb = VEmap.Extent(); i <= aNb; ++i).
            let a_nb = vemak.extent();
            for i in (a_nb_planar_start + 1)..=a_nb {
                let a_s = vemak.at_1(i).clone();
                if a_s.shape_type() != ShapeType::Vertex {
                    continue;
                }

                // Find faces connected to the vertex
                // OCCT L538-552.
                let mut a_lf: Vec<Shape> = Vec::new();
                {
                    let a_le = analyse.ancestors(&a_s);
                    for it_le in &a_le {
                        let a_lea = analyse.ancestors(it_le);
                        for it_lea in &a_lea {
                            if !a_lf.iter().any(|s| s.is_same(it_lea)) {
                                a_lf.push(it_lea.clone());
                            }
                        }
                    }
                }

                if a_lf.len() < 2 {
                    continue;
                }

                // build lists of faces connected to the same vertex by looking for
                // the pairs in which the vertex is alone (not connected to shared edges)
                // OCCT L559-642.
                let mut a_lf1: Vec<Shape> = Vec::new();
                let mut a_lf2: Vec<Shape> = Vec::new();

                for (li, a_fv1) in a_lf.iter().enumerate() {
                    // get edges of first face connected to current vertex
                    let mut a_me: OcctShapeSet = HashMap::new();
                    let p_f1_internal = shape_data_map::seek(&a_dm_fei, a_fv1);
                    let p_le1 = analyse.descendants(a_fv1);
                    let p_le1 = match p_le1 {
                        Some(v) => v,
                        // OCCT L573-575: if (!pLE1) continue.
                        None => continue,
                    };
                    let mut broke = false;
                    for a_e in &p_le1 {
                        if let Some(internal) = p_f1_internal {
                            if set_contains(internal, a_e) {
                                broke = true;
                                break;
                            }
                        }

                        // OCCT L586-593: the TopoDS_Iterator vertex scan.
                        for a_v in bat_sub_shapes(a_e) {
                            if a_s.is_same(&a_v) {
                                set_add(&mut a_me, a_e);
                                break;
                            }
                        }
                    }
                    if broke {
                        // OCCT L595-598: if (itLE1.More()) continue.
                        continue;
                    }
                    let _ = li;

                    // get to the next face in the list
                    // OCCT L600-641.
                    for a_fv2 in &a_lf[(li + 1)..] {
                        let p_f2_internal = shape_data_map::seek(&a_dm_fei, a_fv2);

                        let p_le2 = analyse.descendants(a_fv2);
                        let p_le2 = match p_le2 {
                            Some(v) => v,
                            // OCCT L609-612: if (!pLE2) continue.
                            None => continue,
                        };
                        let mut broke2 = false;
                        for a_ev2 in &p_le2 {
                            if !set_contains(&a_me, a_ev2) {
                                continue;
                            }

                            if let Some(internal) = p_f2_internal {
                                if set_contains(internal, a_ev2) {
                                    // Avoid intersection of faces connected by internal edge
                                    broke2 = true;
                                    break;
                                }
                            }

                            if analyse.has_ancestor(a_ev2) && analyse.ancestors(a_ev2).len() == 2 {
                                // Faces will be intersected through the edge
                                broke2 = true;
                                break;
                            }
                        }

                        if !broke2 {
                            a_lf1.push(a_fv1.clone());
                            a_lf2.push(a_fv2.clone());
                        }
                    }
                }
                //
                // OCCT L644-648.
                if !a_lf1.is_empty() {
                    shape_data_map::bind(&mut a_dm_vlf1, &a_s, a_lf1);
                    shape_data_map::bind(&mut a_dm_vlf2, &a_s, a_lf2);
                }
            }
        }
        //
        // OCCT L652-851: the main intersection walk.
        let a_nb = vemak.extent();
        for i in 1..=a_nb {
            let a_s = vemak.at_1(i).clone();
            //
            let mut e = Shape::null();
            let mut a_lf1: Vec<Shape> = Vec::new();
            let mut a_lf2: Vec<Shape> = Vec::new();
            //
            b_edge = a_s.shape_type() == ShapeType::Edge;
            if b_edge {
                // faces connected by the edge
                // OCCT L669-704.
                e = a_s.clone();
                //
                let l = analyse.type_(&e);
                if l.is_empty() {
                    continue;
                }
                //
                let ot = l[0].my_type;
                if ot != crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity::Convex
                    && ot != crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity::Concave
                {
                    continue;
                }
                //
                if ot == crate::fillet::chfi_ds::ChFiDS_TypeOfConcavity::Concave {
                    cur_side = State::In;
                } else {
                    cur_side = State::Out;
                }
                //-----------------------------------------------------------
                // edge is of the proper type, return adjacent faces.
                //-----------------------------------------------------------
                let anc = analyse.ancestors(&e);
                if anc.len() != 2 {
                    continue;
                }
                //
                f1 = anc[0].clone();
                f2 = anc[anc.len() - 1].clone();
                //
                a_lf1.push(f1.clone());
                a_lf2.push(f2.clone());
            } else {
                // OCCT L706-717: the vertex form.
                if !shape_data_map::is_bound(&a_dm_vlf1, &a_s) {
                    continue;
                }
                //
                a_lf1 = shape_data_map::find(&a_dm_vlf1, &a_s);
                a_lf2 = shape_data_map::find(&a_dm_vlf2, &a_s);
                //
                cur_side = self.my_side;
            }
            //
            // OCCT L719-849: the F1xF2 pair walk.
            for pair in 0..a_lf1.len() {
                f1 = a_lf1[pair].clone();
                f2 = a_lf2[pair].clone();
                //
                // OCCT L726-727: OF1/OF2 = MapSF(F).Face().
                of1 = map_sf_value_face(map_sf, &f1);
                of2 = map_sf_value_face(map_sf, &f2);
                // OCCT L728-739: the MES bind for OF1.
                if !shape_data_map::is_bound(mes, &of1) {
                    let mut enlarge_u = true;
                    let mut enlarge_vfirst = true;
                    let mut enlarge_vlast = true;
                    check_bounds(&f1, analyse, &mut enlarge_u, &mut enlarge_vfirst, &mut enlarge_vlast);
                    en_large_face(&of1, &mut nf1, true, true, enlarge_u, enlarge_vfirst, enlarge_vlast, 1, -1.0, -1.0, -1.0, -1.0);
                    shape_data_map::bind(mes, &of1, nf1.clone());
                } else {
                    nf1 = shape_data_map::find(mes, &of1);
                }
                //
                // OCCT L741-752: the MES bind for OF2.
                if !shape_data_map::is_bound(mes, &of2) {
                    let mut enlarge_u = true;
                    let mut enlarge_vfirst = true;
                    let mut enlarge_vlast = true;
                    check_bounds(&f2, analyse, &mut enlarge_u, &mut enlarge_vfirst, &mut enlarge_vlast);
                    en_large_face(&of2, &mut nf2, true, true, enlarge_u, enlarge_vfirst, enlarge_vlast, 1, -1.0, -1.0, -1.0, -1.0);
                    shape_data_map::bind(mes, &of2, nf2.clone());
                } else {
                    nf2 = shape_data_map::find(mes, &of2);
                }
                //
                // OCCT L754-848.
                if !self.is_done(&nf1, &nf2) {
                    let mut l_int1: Vec<Shape> = Vec::new();
                    let mut l_int2: Vec<Shape> = Vec::new();
                    inter3d(&nf1, &nf2, &mut l_int1, &mut l_int2, cur_side, &e, &f1, &f2);
                    self.set_done(&nf1, &nf2);
                    if !l_int1.is_empty() {
                        self.store(&nf1, &nf2, &l_int1, &l_int2);
                        //
                        // OCCT L763-764: TopoDS_Compound C; B.MakeCompound(C).
                        let mut c = crate::brep_algo::tool::builder_make_compound();
                        //
                        // OCCT L766-775: the prior Build(aS) edge merge.
                        if shape_data_map::is_bound(build, &a_s) {
                            let a_se = shape_data_map::find(build, &a_s);
                            for a_ne in explorer(&a_se, ShapeType::Edge, ShapeType::Shape) {
                                crate::brep_algo::tool::builder_add_compound_shape(&mut c, &a_ne);
                            }
                        }
                        //
                        // OCCT L777-792.
                        for a_ne in &l_int1 {
                            crate::brep_algo::tool::builder_add_compound_shape(&mut c, a_ne);
                            //
                            // keep connection from new edge to shape from which it was created
                            // OCCT L784-786: aDMIntE(aDMIntE.Add(aNE, empty)).
                            let pos = shape_indexed_data_map::add(
                                &mut a_dm_int_e,
                                a_ne,
                                Vec::new(),
                            );
                            shape_indexed_data_map::value_1_mut(&mut a_dm_int_e, pos).push(a_s.clone());
                            // keep connection to faces created the edge as well
                            // OCCT L788-791.
                            shape_data_map::bound(&mut a_dm_int_ff, a_ne).push(f1.clone());
                            shape_data_map::bound(&mut a_dm_int_ff, a_ne).push(f2.clone());
                        }
                        //
                        // OCCT L794: Build.Bind(aS, C).
                        shape_data_map::bind(build, &a_s, c);
                    } else {
                        // OCCT L798: Failed.Append(aS).
                        failed.push(a_s.clone());
                    }
                } else {
                    // IsDone(NF1,NF2)
                    //  Modified by skv - Fri Dec 26 12:20:13 2003 OCC4455
                    // OCCT L804-847.
                    let a_l_int1 = self.my_as_des.descendant(&nf1).to_vec();
                    let a_l_int2 = self.my_as_des.descendant(&nf2).to_vec();

                    if !a_l_int1.is_empty() {
                        let mut c = crate::brep_algo::tool::builder_make_compound();
                        //
                        // OCCT L812-821: the prior Build(aS) edge merge.
                        if shape_data_map::is_bound(build, &a_s) {
                            let a_se = shape_data_map::find(build, &a_s);
                            for a_ne in explorer(&a_se, ShapeType::Edge, ShapeType::Shape) {
                                crate::brep_algo::tool::builder_add_compound_shape(&mut c, &a_ne);
                            }
                        }
                        //
                        // OCCT L823-841.
                        for an_e1 in &a_l_int1 {
                            for an_e2 in &a_l_int2 {
                                if an_e1.is_same(an_e2) {
                                    crate::brep_algo::tool::builder_add_compound_shape(&mut c, an_e1);
                                    //
                                    if let Some(p_ls) =
                                        shape_indexed_data_map::change_seek(&mut a_dm_int_e, an_e1)
                                    {
                                        p_ls.push(a_s.clone());
                                    }
                                }
                            }
                        }
                        // OCCT L842: Build.Bind(aS, C).
                        shape_data_map::bind(build, &a_s, c);
                    } else {
                        // OCCT L846: Failed.Append(aS).
                        failed.push(a_s.clone());
                    }
                }
            }
            //  Modified by skv - Fri Dec 26 12:20:14 2003 OCC4455 End
        }
        //
        // create unique intersection for each localized shared part
        // OCCT L853-1040.
        let a_nb = shape_indexed_data_map::extent(&a_dm_int_e);
        for i in 1..=a_nb {
            let a_ls = shape_indexed_data_map::value_1(&a_dm_int_e, i).clone();
            if a_ls.len() < 2 {
                continue;
            }
            //
            // intersection edge
            // OCCT L869: aE = TopoDS::Edge(aDMIntE.FindKey(i)).
            let a_e = shape_indexed_data_map::find_key_1(&a_dm_int_e, i).clone();
            // faces created the edge
            // OCCT L871-873: aLFF = aDMIntFF.Find(aE); aF1/aF2 = First/Last.
            let a_lff = shape_data_map::find(&a_dm_int_ff, &a_e);
            let a_f1 = a_lff.first().cloned().unwrap_or_else(Shape::null);
            let a_f2 = a_lff.last().cloned().unwrap_or_else(Shape::null);

            // Build really localized blocks from the original shapes in <aLS>:
            // 1. Find edges from original faces connected to two or more shapes in <aLS>;
            // 2. Make connexity blocks from edges in <aLS> and found connection edges;
            // 3. Check if the vertices from <aLS> are not connected by these connection edges:
            //    a. If so - add these vertices to Connexity Block containing the corresponding
            //       connexity edge;
            //    b. If not - add this vertex to list of connexity blocks
            // 4. Create unique intersection edge for each connexity block

            // list of vertices
            // OCCT L885-903.
            let mut a_lv: Vec<Shape> = Vec::new();
            // compound of edges to build connexity blocks
            let mut a_ce = crate::brep_algo::tool::builder_make_compound();
            let mut a_ms: OcctShapeSet = HashMap::new();
            for a_s in &a_ls {
                set_add(&mut a_ms, a_s);
                if a_s.shape_type() == ShapeType::Edge {
                    crate::brep_algo::tool::builder_add_compound_shape(&mut a_ce, a_s);
                } else {
                    a_lv.push(a_s.clone());
                }
            }
            //
            // look for additional edges to connect the shared parts
            // OCCT L906-948.
            let mut a_me_connection: OcctShapeSet = HashMap::new();
            for j in 0..2usize {
                let a_f = if j == 0 { &a_f1 } else { &a_f2 };
                //
                for a_ef in explorer(a_f, ShapeType::Edge, ShapeType::Shape) {
                    if set_contains(&a_ms, &a_ef) || set_contains(&a_me_connection, &a_ef) {
                        continue;
                    }
                    //
                    let (a_v1, a_v2) = top_exp_vertices(&a_ef);
                    //
                    // find parts to which the edge is connected
                    // OCCT L924-940.
                    let mut i_counter = 0i32;
                    for a_s in &a_ls {
                        // iterator is not suitable here, because aS may be a vertex
                        for a_v in explorer(a_s, ShapeType::Vertex, ShapeType::Shape) {
                            if a_v.is_same(&a_v1) || a_v.is_same(&a_v2) {
                                i_counter += 1;
                                break;
                            }
                        }
                    }
                    //
                    // OCCT L942-946.
                    if i_counter >= 2 {
                        crate::brep_algo::tool::builder_add_compound_shape(&mut a_ce, &a_ef);
                        set_add(&mut a_me_connection, &a_ef);
                    }
                }
            }
            //
            // OCCT L950-951: MakeConnexityBlocks(aCE, VERTEX, EDGE, aLCBE).
            let mut a_lcbe = crate::offset::brep_offset_tool::make_connexity_blocks_edges(&bat_sub_shapes(&a_ce));
            //
            // create connexity blocks for alone vertices
            // OCCT L953-986.
            let mut a_lcbv: Vec<Shape> = Vec::new();
            for a_v in &a_lv {
                // check if this vertex is contained in some connexity block of edges
                let mut placed = false;
                for a_cb in &mut a_lcbe {
                    let mut found = false;
                    for a_expv in explorer(a_cb, ShapeType::Vertex, ShapeType::Shape) {
                        if a_v.is_same(&a_expv) {
                            // OCCT L969: B.Add(aCB, aV).
                            crate::brep_algo::tool::builder_add_compound_shape(a_cb, a_v);
                            found = true;
                            break;
                        }
                    }
                    if found {
                        placed = true;
                        break;
                    }
                }
                //
                if !placed {
                    // OCCT L981-985.
                    let mut a_cv = crate::brep_algo::tool::builder_make_compound();
                    crate::brep_algo::tool::builder_add_compound_shape(&mut a_cv, a_v);
                    a_lcbv.push(a_cv);
                }
            }
            //
            // OCCT L988: aLCBE.Append(aLCBV).
            a_lcbe.extend(a_lcbv);
            //
            // OCCT L990-993.
            if a_lcbe.len() == 1 {
                continue;
            }
            //
            // OCCT L995-996: aNF1/aNF2 = MES(MapSF(aF1/aF2).Face()).
            let a_nf1 = shape_data_map::find(mes, &map_sf_value_face(map_sf, &a_f1));
            let a_nf2 = shape_data_map::find(mes, &map_sf_value_face(map_sf, &a_f2));
            //
            // OCCT L998-1039: the new-edge rebinding walk.
            for a_cb in a_lcbe.iter().skip(1) {
                // make new edge with different tedge instance
                // OCCT L1002-1009.
                let (a_v1, a_v2) = top_exp_vertices(&a_e);
                let (a_t1, a_t2) = brep_tool_range(&a_e);
                //
                let mut a_new_edge = Shape::null();
                bop_split_edge_public(&a_e, &a_v1, a_t1, &a_v2, a_t2, &mut a_new_edge);
                //
                // OCCT L1011-1012: myAsDes->Add(aNF1/aNF2, aNewEdge).
                self.my_as_des.add(&a_nf1, &a_new_edge);
                self.my_as_des.add(&a_nf2, &a_new_edge);
                //
                // OCCT L1014-1038.
                for a_s in bat_sub_shapes(a_cb) {
                    if set_contains(&a_me_connection, &a_s) {
                        continue;
                    }
                    // OCCT L1023: aCI = Build.ChangeFind(aS).
                    if !shape_data_map::is_bound(build, &a_s) {
                        continue;
                    }
                    let a_ci = shape_data_map::find(build, &a_s);
                    //
                    // OCCT L1025-1037.
                    let mut a_new_ci = crate::brep_algo::tool::builder_make_compound();
                    for a_sx in explorer(&a_ci, ShapeType::Edge, ShapeType::Shape) {
                        if !a_sx.is_same(&a_e) {
                            crate::brep_algo::tool::builder_add_compound_shape(&mut a_new_ci, &a_sx);
                        }
                    }
                    crate::brep_algo::tool::builder_add_compound_shape(&mut a_new_ci, &a_new_edge);
                    shape_data_map::change_find(build, &a_s).clone_from(&a_new_ci);
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // OCCT BRepOffset_Inter3d::ContextIntByInt (cxx L1045-1272).
    // -----------------------------------------------------------------------

    /// OCCT BRepOffset_Inter3d::ContextIntByInt(ContextFaces,
    /// ExtentContext, MapSF, Analyse, MES, Build, Failed, theRange,
    /// bIsPlanar) (cxx L1045-1272) — computes the intersection with the not
    /// offset faces.
    #[allow(clippy::too_many_arguments)]
    pub fn context_int_by_int(
        &mut self,
        context_faces: &OcctIndexedShapeMap,
        extent_context: bool,
        map_sf: &ShapeDataMap<BRepOffsetOffset>,
        analyse: &BRepOffsetAnalyse,
        mes: &mut ShapeDataMap<Shape>,
        build: &mut ShapeDataMap<Shape>,
        failed: &mut Vec<Shape>,
        _the_range: &rcad_kernel::message::ProgressScope,
        b_is_planar: bool,
    ) {
        let mut mv: OcctShapeSet = HashMap::new();
        let _ = &mut mv;
        let mut of = Shape::null();
        let mut nf = Shape::null();
        let mut wcf = Shape::null();
        let mut b_edge;

        let a_nb = context_faces.extent();
        // OCCT L1067-1076: the Touched / EnLargeFace pre-pass.
        for i in 1..=a_nb {
            let cf = context_faces.at_1(i).clone();
            self.my_touched.add(&cf);
            if extent_context {
                en_large_face(&cf, &mut nf, false, false, true, true, true, 1, -1.0, -1.0, -1.0, -1.0);
                shape_data_map::bind(mes, &cf, nf.clone());
            }
        }
        let side = State::Out;

        // OCCT L1080-1271: the main walk.
        for i in 1..=a_nb {
            let cf = context_faces.at_1(i).clone();
            if extent_context {
                wcf = shape_data_map::find(mes, &cf);
            } else {
                wcf = cf.clone();
            }

            // OCCT L1096-1103: VEmap (edges + the bIsPlanar vertices).
            let mut vemak = OcctIndexedShapeMap::new();
            for e in explorer(&oriented(&cf, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
                vemak.add(&e);
            }
            //
            if b_is_planar {
                for v in explorer(&oriented(&cf, Orientation::Forward), ShapeType::Vertex, ShapeType::Shape) {
                    vemak.add(&v);
                }
            }
            //
            let a_nb_ve = vemak.extent();
            for j in 1..=a_nb_ve {
                let a_s = vemak.at_1(j).clone();
                //
                b_edge = a_s.shape_type() == ShapeType::Edge;
                //
                let mut e = Shape::null();
                let mut anc: Vec<Shape> = Vec::new();
                //
                if b_edge {
                    // faces connected by the edge
                    //
                    // OCCT L1118-1164.
                    e = a_s.clone();
                    if !analyse.has_ancestor(&e) {
                        //----------------------------------------------------------------
                        // the edges of faces of context that are not in the initial shape
                        // can appear in the result.
                        //----------------------------------------------------------------
                        if !extent_context {
                            self.my_as_des.add(&cf, &e);
                            self.my_new_edges.add(&e);
                        } else if !shape_data_map::is_bound(mes, &e) {
                            // OCCT L1134-1153: the ExtentEdge + vertex
                            // parameter form.
                            let mut ne = Shape::null();
                            let (f_par, l_par) =
                                brep_tool_range(&e);
                            let tol = brep_tool_tolerance(&e);
                            Self::extent_edge(&cf, &e, &mut ne);
                            let (v1, v2) = top_exp_vertices(&e);
                            ne.orientation = Orientation::Forward;
                            self.my_as_des.add(&ne, &oriented(&v1, Orientation::Reversed));
                            self.my_as_des.add(&ne, &oriented(&v2, Orientation::Forward));
                            // OCCT L1144-1147: the INTERNAL-vertex
                            // parameter updates.
                            let mut ne_fwd = ne.clone();
                            ne_fwd.orientation = Orientation::Forward;
                            b_update_vertex_on_edge(
                                &mut ne_fwd,
                                &oriented(&v1, Orientation::Internal),
                                f_par,
                                tol,
                            );
                            b_update_vertex_on_edge(
                                &mut ne_fwd,
                                &oriented(&v2, Orientation::Internal),
                                l_par,
                                tol,
                            );
                            // OCCT L1150-1153.
                            ne.orientation = e.orientation;
                            self.my_as_des.add(&cf, &ne);
                            self.my_new_edges.add(&ne);
                            shape_data_map::bind(mes, &e, ne);
                        } else {
                            // OCCT L1155-1161.
                            let ne = shape_data_map::find(mes, &e);
                            let a_local_shape = oriented(&ne, e.orientation);
                            self.my_as_des.add(&cf, &a_local_shape);
                        }
                        continue;
                    }
                    anc = analyse.ancestors(&e);
                } else {
                    // faces connected by the vertex
                    //
                    // OCCT L1167-1209.
                    if !analyse.has_ancestor(&a_s) {
                        continue;
                    }
                    //
                    let a_le = analyse.ancestors(&a_s);
                    for it_value in &a_le {
                        let a_e = it_value;
                        //
                        if brep_tool_degenerated(a_e) {
                            continue;
                        }
                        //
                        if vemak.contains(a_e) {
                            continue;
                        }
                        //
                        let a_lf = analyse.ancestors(a_e);
                        for a_f in &a_lf {
                            let mut b_add = true;
                            for a_ef in explorer(a_f, ShapeType::Edge, ShapeType::Shape) {
                                b_add = !vemak.contains(&a_ef);
                                if !b_add {
                                    break;
                                }
                            }
                            if b_add {
                                anc.push(a_f.clone());
                            }
                        }
                    }
                }
                //
                // OCCT L1212-1269: the F walk.
                for it_value in &anc {
                    let f = it_value;
                    // OCCT L1216-1218: OF = MapSF(F).Face(); OE =
                    // MapSF(F).Generated(E).
                    of = map_sf_value_face(map_sf, f);
                    let oe = map_sf_value_generated(map_sf, f, &e);
                    // OCCT L1220-1228: the MES bind.
                    if !shape_data_map::is_bound(mes, &of) {
                        en_large_face(&of, &mut nf, true, true, true, true, true, 1, -1.0, -1.0, -1.0, -1.0);
                        shape_data_map::bind(mes, &of, nf.clone());
                    } else {
                        nf = shape_data_map::find(mes, &of);
                    }
                    // OCCT L1229-1268.
                    if !self.is_done(&nf, &cf) {
                        let mut l_int1: Vec<Shape> = Vec::new();
                        let mut l_int2: Vec<Shape> = Vec::new();
                        let mut loe: Vec<Shape> = Vec::new();
                        loe.push(oe);
                        inter3d(&wcf, &nf, &mut l_int1, &mut l_int2, side, &e, &cf, f);
                        self.set_done(&nf, &cf);
                        if !l_int1.is_empty() {
                            self.store(&cf, &nf, &l_int1, &l_int2);
                            // OCCT L1239-1262.
                            if l_int1.len() == 1 && !shape_data_map::is_bound(build, &a_s) {
                                shape_data_map::bind(build, &a_s, l_int1[0].clone());
                            } else {
                                let mut c = crate::brep_algo::tool::builder_make_compound();
                                if shape_data_map::is_bound(build, &a_s) {
                                    let a_se = shape_data_map::find(build, &a_s);
                                    for a_ne in explorer(&a_se, ShapeType::Edge, ShapeType::Shape) {
                                        crate::brep_algo::tool::builder_add_compound_shape(&mut c, &a_ne);
                                    }
                                }
                                //
                                for it_l in &l_int1 {
                                    crate::brep_algo::tool::builder_add_compound_shape(&mut c, it_l);
                                }
                                shape_data_map::bind(build, &a_s, c);
                            }
                        } else {
                            // OCCT L1266: Failed.Append(aS).
                            failed.push(a_s.clone());
                        }
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // OCCT BRepOffset_Inter3d::ContextIntByArc (cxx L1276-1496).
    // -----------------------------------------------------------------------

    /// OCCT BRepOffset_Inter3d::ContextIntByArc(ContextFaces, InSide,
    /// Analyse, InitOffsetFace, InitOffsetEdge, theRange) (cxx
    /// L1276-1496) — computes the connections of the not offset faces that
    /// have to be connected by arcs.
    #[allow(clippy::too_many_arguments)]
    pub fn context_int_by_arc(
        &mut self,
        context_faces: &OcctIndexedShapeMap,
        in_side: bool,
        analyse: &BRepOffsetAnalyse,
        init_offset_face: &BRepAlgoImage,
        init_offset_edge: &mut BRepAlgoImage,
        _the_range: &rcad_kernel::message::ProgressScope,
    ) {
        let mut l_int1: Vec<Shape> = Vec::new();
        let mut l_int2: Vec<Shape> = Vec::new();
        let mut mv: OcctShapeSet = HashMap::new();
        let mut of1 = Shape::null();
        let mut oe = Shape::null();
        let null_edge = Shape::null();
        let null_face = Shape::null();

        // OCCT L1294-1298: the Touched pre-pass.
        for j in 1..=context_faces.extent() {
            let cf = context_faces.at_1(j).clone();
            self.my_touched.add(&cf);
        }

        // OCCT L1300-1495: the main walk.
        for j in 1..=context_faces.extent() {
            let cf = context_faces.at_1(j).clone();
            for exp_value in explorer(&oriented(&cf, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
                let e = exp_value;
                if !analyse.has_ancestor(&e) {
                    if in_side {
                        // OCCT L1315: myAsDes->Add(CF, E).
                        self.my_as_des.add(&cf, &e);
                    } else {
                        // OCCT L1319-1346: the InitOffsetEdge forms.
                        let mut ne = Shape::null();
                        if !init_offset_edge.has_image(&e) {
                            let (f_par, l_par) =
                                brep_tool_range(&e);
                            let tol = brep_tool_tolerance(&e);
                            Self::extent_edge(&cf, &e, &mut ne);
                            let (v1, v2) = top_exp_vertices(&e);
                            ne.orientation = Orientation::Forward;
                            self.my_as_des.add(&ne, &oriented(&v1, Orientation::Reversed));
                            self.my_as_des.add(&ne, &oriented(&v2, Orientation::Forward));
                            // OCCT L1331-1334.
                            let mut ne_fwd = ne.clone();
                            ne_fwd.orientation = Orientation::Forward;
                            b_update_vertex_on_edge(
                                &mut ne_fwd,
                                &oriented(&v1, Orientation::Internal),
                                f_par,
                                tol,
                            );
                            b_update_vertex_on_edge(
                                &mut ne_fwd,
                                &oriented(&v2, Orientation::Internal),
                                l_par,
                                tol,
                            );
                            // OCCT L1337-1339.
                            ne.orientation = e.orientation;
                            self.my_as_des.add(&cf, &ne);
                            init_offset_edge.bind(&e, &ne);
                        } else {
                            // OCCT L1341-1345.
                            ne = init_offset_edge.image(&e).first().cloned().unwrap_or_else(Shape::null);
                            self.my_as_des.add(&cf, &oriented(&ne, e.orientation));
                        }
                    }
                    continue;
                }
                // OCCT L1349: OE.Nullify().
                oe = Shape::null();
                //---------------------------------------------------
                // OF1 parallel facee generated by the ancestor of E.
                //---------------------------------------------------
                // OCCT L1353-1355.
                let si = analyse.ancestors(&e).first().cloned().unwrap_or_else(Shape::null);
                of1 = init_offset_face.image(&si).first().cloned().unwrap_or_else(Shape::null);
                oe = init_offset_edge.image(&e).first().cloned().unwrap_or_else(Shape::null);

                {
                    // Check if OE has pcurve in CF
                    // OCCT L1360-1369.
                    let c1 = brep_tool_curve_on_surface(&oe, &cf);
                    let c2 = brep_tool_curve_on_surface(&oe, &of1);

                    if c1.is_none() || c2.is_none() {
                        continue;
                    }
                }

                //--------------------------------------------------
                // MAJ of OE on cap CF.
                //--------------------------------------------------
                // OCCT L1379-1387.
                l_int1.clear();
                l_int1.push(oe.clone());
                l_int2.clear();
                let mut an_ori1 = Orientation::Forward;
                let mut an_ori2 = Orientation::Forward;
                crate::offset::brep_offset_tool::orient_section(&oe, &cf, &of1, &mut an_ori1, &mut an_ori2);
                an_ori1 = crate::brep_algo::tool::top_abs_reverse(an_ori1);
                if let Some(first) = l_int1.first_mut() {
                    first.orientation = an_ori1;
                }
                self.store(&cf, &of1, &l_int1, &l_int2);

                //------------------------------------------------------
                // Processing of offsets on the ancestors of vertices.
                //------------------------------------------------------
                // OCCT L1392-1447.
                let (v0, v1) = top_exp_vertices(&e);
                let v = [v0, v1];
                for i in 0..2 {
                    if !set_add(&mut mv, &v[i]) {
                        continue;
                    }
                    let le = analyse.ancestors(&v[i]);
                    for it_le in &le {
                        let ev = it_le;
                        if init_offset_face.has_image(ev) {
                            //-------------------------------------------------
                            // OF1 parallel face generated by an ancestor edge of V[i].
                            //-------------------------------------------------
                            // OCCT L1411-1412.
                            of1 = init_offset_face.image(ev).first().cloned().unwrap_or_else(Shape::null);
                            oe = init_offset_edge.image(&v[i]).first().cloned().unwrap_or_else(Shape::null);

                            {
                                // Check if OE has pcurve in CF and OF1
                                // OCCT L1417-1426.
                                let c1 = brep_tool_curve_on_surface(&oe, &cf);
                                let c2 = brep_tool_curve_on_surface(&oe, &of1);

                                if c1.is_none() || c2.is_none() {
                                    continue;
                                }
                            }

                            //--------------------------------------------------
                            // MAj of OE on cap CF.
                            //--------------------------------------------------
                            // OCCT L1436-1444.
                            l_int1.clear();
                            l_int1.push(oe.clone());
                            l_int2.clear();
                            let mut o1 = Orientation::Forward;
                            let mut o2 = Orientation::Forward;
                            crate::offset::brep_offset_tool::orient_section(&oe, &cf, &of1, &mut o1, &mut o2);
                            o1 = crate::brep_algo::tool::top_abs_reverse(o1);
                            if let Some(first) = l_int1.first_mut() {
                                first.orientation = o1;
                            }
                            self.store(&cf, &of1, &l_int1, &l_int2);
                            let _ = o2;
                        }
                    }
                }
            }

            // OCCT L1450-1494: the vertex walk.
            for exp_value in explorer(&oriented(&cf, Orientation::Forward), ShapeType::Vertex, ShapeType::Shape) {
                let v = exp_value;
                if !analyse.has_ancestor(&v) {
                    continue;
                }
                let le = analyse.ancestors(&v);
                for it_le in &le {
                    let ev = it_le;
                    let lf = analyse.ancestors(ev);
                    for it_lf in &lf {
                        let fev = it_lf;
                        //-------------------------------------------------
                        // OF1 parallel face generated by uneFace ancestor of V[i].
                        //-------------------------------------------------
                        // OCCT L1470.
                        of1 = init_offset_face.image(fev).first().cloned().unwrap_or_else(Shape::null);
                        if !self.is_done(&of1, &cf) {
                            //-------------------------------------------------------
                            // Find if one of edges of OF1 has no trace in CF.
                            //-------------------------------------------------------
                            // OCCT L1476-1481.
                            let mut loe: Vec<Shape> = Vec::new();
                            for e2 in explorer(&oriented(&of1, Orientation::Forward), ShapeType::Edge, ShapeType::Shape) {
                                loe.push(e2);
                            }
                            //-------------------------------------------------------
                            // If no trace try intersection.
                            //-------------------------------------------------------
                            // OCCT L1485-1489.
                            if !try_project(&cf, &of1, &loe, &mut l_int1, &mut l_int2, self.my_side, self.my_tol)
                                || l_int1.is_empty()
                            {
                                inter3d(&cf, &of1, &mut l_int1, &mut l_int2, self.my_side, &null_edge, &null_face, &null_face);
                            }
                            self.store(&cf, &of1, &l_int1, &l_int2);
                        }
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // OCCT BRepOffset_Inter3d::SetDone / IsDone / Store (cxx L1500-1554).
    // -----------------------------------------------------------------------

    /// OCCT BRepOffset_Inter3d::SetDone(F1, F2) (cxx L1500-1514).
    pub fn set_done(&mut self, f1: &Shape, f2: &Shape) {
        if !shape_data_map::is_bound(&self.my_done, f1) {
            shape_data_map::bind(&mut self.my_done, f1, Vec::new());
        }
        shape_data_map::change_find(&mut self.my_done, f1).push(f2.clone());
        if !shape_data_map::is_bound(&self.my_done, f2) {
            shape_data_map::bind(&mut self.my_done, f2, Vec::new());
        }
        shape_data_map::change_find(&mut self.my_done, f2).push(f1.clone());
    }

    /// OCCT BRepOffset_Inter3d::IsDone(F1, F2) (cxx L1518-1532).
    pub fn is_done(&self, f1: &Shape, f2: &Shape) -> bool {
        if shape_data_map::is_bound(&self.my_done, f1) {
            for it_value in shape_data_map::value(&self.my_done, f1) {
                if it_value.is_same(f2) {
                    return true;
                }
            }
        }
        false
    }

    /// OCCT BRepOffset_Inter3d::Store(F1, F2, LInt1, LInt2) (cxx
    /// L1536-1554).
    fn store(&mut self, f1: &Shape, f2: &Shape, l_int1: &[Shape], l_int2: &[Shape]) {
        if !l_int1.is_empty() {
            self.my_touched.add(f1);
            self.my_touched.add(f2);
            self.my_as_des.add_list(f1, l_int1);
            self.my_as_des.add_list(f2, l_int2);
            for it_value in l_int1 {
                self.my_new_edges.add(it_value);
            }
        }
        self.set_done(f1, f2);
    }
}

// ---------------------------------------------------------------------------
// Local helpers (the cxx statics / reduced re-hosts of this module).
// ---------------------------------------------------------------------------

/// OCCT Bnd_Box (the CompletInt face box; the BRepBndLib::Add reduced
/// sampling re-host — architecture difference #36).
#[derive(Debug, Clone)]
struct FaceBndBox {
    min: DVec3,
    max: DVec3,
}

impl FaceBndBox {
    /// OCCT Bnd_Tools::Bnd2BVH — the Aabb form.
    fn to_aabb(&self) -> crate::bop::tools::box_tree::Aabb {
        crate::bop::tools::box_tree::Aabb {
            min: self.min,
            max: self.max,
            gap: 0.0,
        }
    }
}

/// OCCT BRepBndLib::Add(F, Box) — the reduced sampling re-host (surface
/// domain corners + the edge curve samples; the BndLib gap annotation).
fn brep_bnd_lib_add_face(f: &Shape) -> FaceBndBox {
    let mut min = DVec3::splat(f64::INFINITY);
    let mut max = DVec3::splat(f64::NEG_INFINITY);
    // The face surface domain corners.
    if let Some(s) = crate::offset::brep_offset_tool::face_surface_of(f) {
        let d = s.default_domain();
        for (u, v) in [(d[0], d[2]), (d[1], d[2]), (d[0], d[3]), (d[1], d[3])] {
            if u.is_finite() && v.is_finite() {
                let p = s.point_at(u, v);
                min = min.min(p);
                max = max.max(p);
            }
        }
    }
    // The edge samples.
    for e in explorer(f, ShapeType::Edge, ShapeType::Shape) {
        if let Some((c, f0, l0)) = brep_tool_curve(&e) {
            use rcad_kernel::geom::{CurveEval, SurfaceEval};
            for t in [f0, 0.5 * (f0 + l0), l0] {
                let p = CurveEval::point_at(&c, t);
                min = min.min(p);
                max = max.max(p);
            }
        }
    }
    FaceBndBox { min, max }
}

/// OCCT MapSF(F).Face() — the BRepOffset_Offset::Face() accessor form of
/// the NCollection_DataMap value.
fn map_sf_value_face(map_sf: &ShapeDataMap<BRepOffsetOffset>, f: &Shape) -> Shape {
    map_sf
        .get(&bat::shape_key(f))
        .map(|e| e.1.face())
        .expect("MapSF(F) unbound")
}

/// OCCT MapSF(F).Generated(E) — the BRepOffset_Offset::Generated() accessor
/// form.
fn map_sf_value_generated(map_sf: &ShapeDataMap<BRepOffsetOffset>, f: &Shape, e: &Shape) -> Shape {
    map_sf
        .get(&bat::shape_key(f))
        .map(|o| o.1.generated(e))
        .expect("MapSF(F) unbound")
}

/// OCCT TopoDS_Iterator(aE) — the stored sub-shape walk.
fn bat_sub_shapes(s: &Shape) -> Vec<Shape> {
    crate::brep_algo::tool::sub_shapes(s)
}

/// OCCT BOPTools_AlgoTools::MakeSplitEdge(aE, aV1, aT1, aV2, aT2, aNewE)
/// (BOPTools_AlgoTools.cxx L910-945) — the reduced re-host: the empty-copied
/// edge with the trimmed range and the extremity vertices (the rcad
/// DS-bound re-host lives in bop::tools::algo_tools).
fn bop_split_edge_public(e: &Shape, v1: &Shape, t1: f64, v2: &Shape, t2: f64, new_e: &mut Shape) {
    *new_e = empty_copied(e);
    new_e.orientation = Orientation::Forward;
    if let rcad_kernel::topo::topods::TShape::Edge(ed) =
        std::sync::Arc::make_mut(&mut new_e.data)
    {
        ed.range = [t1, t2];
        ed.first = oriented(v1, Orientation::Forward);
        ed.last = oriented(v2, Orientation::Reversed);
        ed.my_shapes.push(oriented(v1, Orientation::Forward));
        ed.my_shapes.push(oriented(v2, Orientation::Reversed));
    }
}

