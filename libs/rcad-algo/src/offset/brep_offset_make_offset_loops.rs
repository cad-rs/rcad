// OCCT BRepOffset_MakeLoops.cxx L41-660 + BRepOffset_MakeLoops.hxx L33-77 —
// 1:1 translation (the Build / BuildOnContext / BuildFaces member bodies of
// the loop builder, split out of brep_offset_make_offset.rs for the
// 2000-line module rule; the architecture-difference numbering continues
// there: #38 map forms, #40 progress flattening, #41 builder pool, #43 this
// unit, #4/#5 the Arc-payload and no-op-cast bridges, #56 TopoDS_Iterator).

use std::collections::HashMap;

use rcad_kernel::topo::topods::{Orientation, ShapeType};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_tool::{
    oriented, shape_data_map, top_exp_vertices, OcctShapeSet,
};
use super::brep_offset_make_offset::DataMapOfShapeShape;
use super::brep_offset_tool_d::BRepOffsetAnalyse;

use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::r#loop::BRepAlgoLoop;
use crate::brep_algo::tool as bat;
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_tolerance;

/// OCCT BRepOffset_MakeLoops (TKOffset/BRepOffset/BRepOffset_MakeLoops.hxx
/// L33-77) — the loop builder of MakeLoops / MakeFaces (architecture
/// difference #43: the Build leaf stays a GAP; BuildOnContext / BuildFaces
/// carry the cxx L216-660 translation).
pub(crate) struct BRepOffsetMakeLoops {
    /// OCCT hxx L60: myVerVerMap — the vertex substitution map carried
    /// between the Loops passes.
    pub(crate) my_ver_ver_map: DataMapOfShapeShape,
}

impl BRepOffsetMakeLoops {
    /// OCCT BRepOffset_MakeLoops::BRepOffset_MakeLoops().
    pub fn new() -> Self {
        BRepOffsetMakeLoops {
            my_ver_ver_map: HashMap::new(),
        }
    }

    /// OCCT BRepOffset_MakeLoops::Build(LF, AsDes, Image, theImageVV,
    /// theRange) (cxx L41-193).
    pub fn build(
        &mut self,
        lf: &Vec<Shape>,
        as_des: &mut BRepAlgoAsDes,
        image_offset: &mut BRepAlgoImage,
        image_vv: &mut BRepAlgoImage,
    ) {
        // OCCT L46-50: the iterators, Loops + VerticesForSubstitute +
        // SetImageVV (arch. diff. #40: the progress scopes are result-inert
        // and flattened — the `!aPS.More()` user-breaks are the NoopProgress
        // never-break form).
        let mut a_loops = BRepAlgoLoop::new();
        // OCCT L48: Loops.VerticesForSubstitute(myVerVerMap).
        a_loops.vertices_for_substitute(&self.my_ver_ver_map);
        // OCCT L49: Loops.SetImageVV(theImageVV) — the loop only reads the
        // VV image (cxx L547-548 of loop.rs), so the rcad by-value carrier
        // takes a clone of the handle-shared image (arch. diff. #4).
        a_loops.set_image_vv(image_vv.clone());
        // OCCT L55: Message_ProgressScope aPS1 "Init loops".

        // OCCT L56: for (; it.More(); it.Next(), aPS1.Next()).
        for it in lf.iter() {
            // OCCT L63: const TopoDS_Face& F = TopoDS::Face(it.Value()).
            let f = it;
            //---------------------------
            // Initialization of Loops.
            //---------------------------
            // OCCT L66: Loops.Init(F).
            a_loops.init(f);
            //-----------------------------
            // return edges of F.
            //-----------------------------
            // OCCT L70: const NCollection_List& LE = AsDes->Descendant(F).
            let le: Vec<Shape> = as_des.descendant(f).to_vec();
            // OCCT L71: NCollection_List<TopoDS_Shape> AddedEdges.
            let mut added_edges: Vec<Shape> = Vec::new();

            // OCCT L73: for (itl.Initialize(LE); itl.More(); itl.Next()).
            for itl in &le {
                // OCCT L75: TopoDS_Edge E = TopoDS::Edge(itl.Value()).
                let e = itl;
                // OCCT L76: if (Image.HasImage(E)).
                if image_offset.has_image(e) {
                    //-------------------------------------------
                    // E was already cut in another face.
                    // Return the cut edges reorientate them as E.
                    // See pb for the edges that have disappeared?
                    //-------------------------------------------
                    // OCCT L83: const NCollection_List& LCE = Image.Image(E).
                    let lce = image_offset.image(e);
                    // OCCT L84: for (itLCE.Initialize(LCE); ...).
                    for it_lce in &lce {
                        // OCCT L86: CE = itLCE.Value().Oriented(E.Orientation()).
                        let ce = oriented(it_lce, e.orientation);
                        // OCCT L87: Loops.AddConstEdge(TopoDS::Edge(CE)).
                        a_loops.add_const_edge(&ce);
                    }
                } else {
                    // OCCT L92: Loops.AddEdge(E, AsDes->Descendant(E)).
                    let lv = as_des.descendant(e).to_vec();
                    a_loops.add_edge(e, &lv);
                    // OCCT L93: AddedEdges.Append(E).
                    added_edges.push(e.clone());
                }
            }
            //------------------------
            // Unwind.
            //------------------------
            // OCCT L98-99: Loops.Perform(); Loops.WiresToFaces().
            a_loops.perform();
            a_loops.wires_to_faces();
            //------------------------
            // MAJ SD.
            //------------------------
            // OCCT L104: const NCollection_List& NF = Loops.NewFaces().
            let nf = a_loops.new_faces().to_vec();
            //-----------------------
            // F => New faces;
            //-----------------------
            // OCCT L108: Image.Bind(F, NF).
            image_offset.bind_list(f, &nf);

            // OCCT L110: for (itAdded.Initialize(AddedEdges); ...).
            for it_added in &added_edges {
                // OCCT L113: const TopoDS_Edge& E =
                // TopoDS::Edge(itAdded.Value()).
                let e = it_added;
                //-----------------------
                //  E => New edges;
                //-----------------------
                // OCCT L118: const NCollection_List& LoopNE =
                // Loops.NewEdges(E).
                // OCCT L119-125: HasImage ? Image.Add(E, LoopNE) :
                // Image.Bind(E, LoopNE).
                if image_offset.has_image(e) {
                    image_offset.add_list(e, a_loops.new_edges(e));
                } else {
                    image_offset.bind_list(e, a_loops.new_edges(e));
                }
            }
        }
        // OCCT L128: Loops.GetVerticesForSubstitute(myVerVerMap).
        a_loops.get_vertices_for_substitute(&mut self.my_ver_ver_map);
        // OCCT L129-132: if (myVerVerMap.IsEmpty()) return.
        if self.my_ver_ver_map.is_empty() {
            return;
        }
        // OCCT L133: BRep_Builder BB; L135: aPS2 "Building loops".
        // OCCT L136: for (it.Initialize(LF); it.More(); it.Next(), aPS2.Next()).
        for it in lf.iter() {
            // OCCT L139: TopoDS_Shape F = it.Value().
            let f = it;
            // OCCT L140-141: NCollection_List LIF; Image.LastImage(F, LIF).
            let mut l_if: Vec<Shape> = Vec::new();
            image_offset.last_image(f, &mut l_if);
            // OCCT L142: for (itl.Initialize(LIF); itl.More(); itl.Next()).
            for itl in &l_if {
                // OCCT L144: const TopoDS_Shape& IF = itl.Value().
                let i_f = itl;
                // OCCT L145: TopExp_Explorer EdExp(IF, TopAbs_EDGE).
                for ed_exp in bat::explorer(i_f, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT L147: TopoDS_Shape E = EdExp.Current().
                    let mut e = ed_exp;
                    // OCCT L148-153: VList over TopoDS_Iterator(E)
                    // (arch. diff. #56: bat::sub_shapes; the null filter
                    // carries the rcad nullified extremities).
                    let v_list: Vec<Shape> =
                        bat::sub_shapes(&e).into_iter().filter(|v| !v.is_null()).collect();
                    // OCCT L154: for (itlv over VList).
                    for v in &v_list {
                        // OCCT L156: if (myVerVerMap.IsBound(V)).
                        if shape_data_map::is_bound(&self.my_ver_ver_map, v) {
                            // OCCT L158: TopoDS_Shape NewV = myVerVerMap(V).
                            let mut new_v = shape_data_map::find(&self.my_ver_ver_map, v);
                            // OCCT L159: E.Free(true).
                            bat::builder_set_free(&mut e, true);
                            // OCCT L160: NewV.Orientation(V.Orientation()).
                            new_v.orientation = v.orientation;
                            // OCCT L161-166: the TV/NewTV tolerance merge and
                            // the points list append (arch. diff. #4).
                            let (a_tv_tol, a_tv_points) = match v.data.as_ref() {
                                rcad_kernel::topo::topods::TShape::Vertex(a_tv) => {
                                    (a_tv.tolerance, a_tv.points.clone())
                                }
                                _ => (0.0, Vec::new()),
                            };
                            if let rcad_kernel::topo::topods::TShape::Vertex(a_new_tv) =
                                std::sync::Arc::make_mut(&mut new_v.data)
                            {
                                if a_tv_tol > a_new_tv.tolerance {
                                    a_new_tv.tolerance = a_tv_tol;
                                }
                                a_new_tv.points.extend(a_tv_points);
                            }
                            // OCCT L167: AsDes->Replace(V, NewV).
                            as_des.replace(v, &new_v);
                            // OCCT L168-169: BB.Remove(E, V); BB.Add(E, NewV).
                            bat::builder_remove_edge_vertex(&mut e, v);
                            bat::builder_add_edge_vertex(&mut e, &new_v);
                        }
                    }
                }
            }
        }
    }

    /// OCCT BRepOffset_MakeLoops::BuildOnContext(LC, Analyse, AsDes, Image,
    /// InSide, theRange) (cxx L216-451).
    pub fn build_on_context(
        &mut self,
        lc: &Vec<Shape>,
        analyse: &BRepOffsetAnalyse,
        as_des: &mut BRepAlgoAsDes,
        image_offset: &mut BRepAlgoImage,
        in_side: bool,
    ) {
        //-----------------------------------------
        // unwinding of caps.
        //-----------------------------------------
        // OCCT L224-229: the iterators, Loops, MapExtent (arch. diff. #38:
        // NCollection_Map -> OcctShapeSet; arch. diff. #40: the
        // Message_ProgressScope is result-inert and flattened — the
        // `!aPS.More()` user-break is the NoopProgress never-break form).
        let mut a_loops = BRepAlgoLoop::new();
        // OCCT L227: Loops.VerticesForSubstitute(myVerVerMap).
        a_loops.vertices_for_substitute(&self.my_ver_ver_map);
        // OCCT L229: NCollection_Map<TopoDS_Shape> MapExtent.
        let mut map_extent: OcctShapeSet = HashMap::new();

        // OCCT L232: for (; it.More(); it.Next(), aPS.Next()).
        for it in lc.iter() {
            // OCCT L235: const TopoDS_Face& F = TopoDS::Face(it.Value()).
            let f = it;
            // OCCT L236: NCollection_Map<TopoDS_Shape> MBound.
            let mut m_bound: OcctShapeSet = HashMap::new();
            //-----------------------------------------------
            // Initialisation of Loops.
            // F is reversed it will be added in myOffC.
            // and myOffC will be reversed in the final result.
            //-----------------------------------------------
            // OCCT L242: TopoDS_Shape aReversedF = F.Reversed().
            let a_reversed_f = bat::reversed(f);
            // OCCT L243-250: if (InSide) Loops.Init(TopoDS::Face(aReversedF));
            // else Loops.Init(F) (the TopoDS::Face cast is a no-op,
            // arch. diff. #5).
            if in_side {
                a_loops.init(&a_reversed_f);
            } else {
                a_loops.init(f);
            }
            //--------------------------------------------------------
            // return edges of F not modified by definition.
            //--------------------------------------------------------
            // OCCT L255: for (exp.Init(F.Oriented(TopAbs_FORWARD),
            // TopAbs_EDGE); exp.More(); exp.Next()).
            for e_exp in bat::explorer(
                &oriented(f, Orientation::Forward),
                ShapeType::Edge,
                ShapeType::Shape,
            ) {
                // OCCT L257: TopoDS_Edge CE = TopoDS::Edge(exp.Current()).
                let ce = e_exp;
                // OCCT L258: MBound.Add(CE).
                m_bound.insert(bat::shape_key(&ce), ce.clone());
                // OCCT L259: if (Analyse.HasAncestor(CE)).
                if analyse.has_ancestor(&ce) {
                    // the stop of cups except for the connectivity stops
                    // between caps.
                    // OCCT L265-271: AddConstEdge(InSide ? CE : CE.Reversed()).
                    if in_side {
                        a_loops.add_const_edge(&ce);
                    } else {
                        a_loops.add_const_edge(&bat::reversed(&ce));
                    }
                }
            }
            //------------------------------------------------------
            // Trace of offsets + connectivity edge between caps.
            //------------------------------------------------------
            // OCCT L278: const NCollection_List<TopoDS_Shape>& LE =
            // AsDes->Descendant(F).
            let le: Vec<Shape> = as_des.descendant(f).to_vec();
            // OCCT L279: NCollection_List<TopoDS_Shape> AddedEdges.
            let mut added_edges: Vec<Shape> = Vec::new();

            // OCCT L281: for (itl.Initialize(LE); itl.More(); itl.Next()).
            for itl in &le {
                // OCCT L283: TopoDS_Edge E = TopoDS::Edge(itl.Value()).
                let e = itl;
                // OCCT L284: if (Image.HasImage(E)).
                if image_offset.has_image(e) {
                    //-------------------------------------------
                    // E was already cut in another face.
                    // Return cut edges and orientate them as E.
                    // See pb for the edges that have disappeared?
                    //-------------------------------------------
                    // OCCT L290: const NCollection_List& LCE = Image.Image(E).
                    let lce = image_offset.image(e);
                    // OCCT L292: for (itLCE.Initialize(LCE); ...).
                    for it_lce in &lce {
                        // OCCT L294: CE = itLCE.Value().Oriented(E.Orientation()).
                        let mut ce = oriented(it_lce, e.orientation);
                        // OCCT L295: if (MapExtent.Contains(E)).
                        if map_extent.contains_key(&bat::shape_key(e)) {
                            // OCCT L296-298: AddConstEdge(CE); continue.
                            a_loops.add_const_edge(&ce);
                            continue;
                        }
                        // OCCT L300-302: if (!MBound.Contains(E)) CE.Reverse().
                        if !m_bound.contains_key(&bat::shape_key(e)) {
                            ce = bat::reversed(&ce);
                        }
                        // OCCT L303-309: AddConstEdge(InSide ? CE : CE.Reversed()).
                        if in_side {
                            a_loops.add_const_edge(&ce);
                        } else {
                            a_loops.add_const_edge(&bat::reversed(&ce));
                        }
                    }
                } else {
                    // OCCT L314: if (IsBetweenCorks(E, AsDes, LContext) &&
                    // AsDes->HasDescendant(E)).
                    if is_between_corks(e, as_des, lc) && as_des.has_descendant(e) {
                        // connection between 2 caps
                        // OCCT L317: MapExtent.Add(E).
                        map_extent.insert(bat::shape_key(e), e.clone());
                        // OCCT L318: NCollection_List<TopoDS_Shape> LV.
                        if in_side {
                            // OCCT L320-323: LV.Append(descendant Reversed);
                            // Loops.AddEdge(E, LV).
                            let mut lv: Vec<Shape> = Vec::new();
                            for it_lce in as_des.descendant(e) {
                                lv.push(bat::reversed(it_lce));
                            }
                            a_loops.add_edge(e, &lv);
                        } else {
                            // OCCT L326-327: Loops.AddEdge(E,
                            // AsDes->Descendant(E)).
                            let lv = as_des.descendant(e).to_vec();
                            a_loops.add_edge(e, &lv);
                        }
                        // OCCT L329: AddedEdges.Append(E).
                        added_edges.push(e.clone());
                    } else if is_between_corks(e, as_des, lc) {
                        // OCCT L333-343: AddConstEdge(InSide ? E : E.Reversed()).
                        if in_side {
                            a_loops.add_const_edge(e);
                        } else {
                            a_loops.add_const_edge(&bat::reversed(e));
                        }
                    } else {
                        // OCCT L345-355: AddConstEdge(InSide ? E.Reversed() : E).
                        if in_side {
                            a_loops.add_const_edge(&bat::reversed(e));
                        } else {
                            a_loops.add_const_edge(e);
                        }
                    }
                }
            }
            //------------------------
            // Unwind.
            //------------------------
            // OCCT L361-362: Loops.Perform(); Loops.WiresToFaces().
            a_loops.perform();
            a_loops.wires_to_faces();
            //------------------------
            // MAJ SD.
            //------------------------
            // OCCT L367: const NCollection_List& NF = Loops.NewFaces().
            let nf = a_loops.new_faces().to_vec();
            //-----------------------
            // F => New faces;
            //-----------------------
            // OCCT L370: Image.Bind(F, NF).
            image_offset.bind_list(f, &nf);

            // OCCT L373: for (itAdded.Initialize(AddedEdges); ...).
            for it_added in &added_edges {
                // OCCT L376: const TopoDS_Edge& E =
                // TopoDS::Edge(itAdded.Value()).
                let e = it_added;
                //-----------------------
                //  E => New edges;
                //-----------------------
                // OCCT L381-387: HasImage ? Image.Add : Image.Bind.
                if image_offset.has_image(e) {
                    image_offset.add_list(e, a_loops.new_edges(e));
                } else {
                    image_offset.bind_list(e, a_loops.new_edges(e));
                }
            }
        }
        // OCCT L392: Loops.GetVerticesForSubstitute(myVerVerMap).
        a_loops.get_vertices_for_substitute(&mut self.my_ver_ver_map);
        // OCCT L393-396: if (myVerVerMap.IsEmpty()) return.
        if self.my_ver_ver_map.is_empty() {
            return;
        }
        // OCCT L397-451: the vertex substitution tail (arch. diff. #4: the
        // OCCT BRep_Builder edits run through shared TShape handles; the
        // rcad builder helpers mutate the local Arc payload).
        // OCCT L399: for (it.Initialize(LContext); it.More(); it.Next()).
        for it in lc.iter() {
            // OCCT L401: TopoDS_Shape F = it.Value().
            let f = it;
            // OCCT L402: NCollection_List LIF; Image.LastImage(F, LIF).
            let mut l_if: Vec<Shape> = Vec::new();
            image_offset.last_image(f, &mut l_if);
            // OCCT L403: for (itl.Initialize(LIF); itl.More(); itl.Next()).
            for itl in &l_if {
                // OCCT L405: const TopoDS_Shape& IF = itl.Value().
                let i_f = itl;
                // OCCT L406: TopExp_Explorer EdExp(IF, TopAbs_EDGE).
                for ed_exp in bat::explorer(i_f, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT L408: TopoDS_Shape E = EdExp.Current().
                    let mut e = ed_exp;
                    // OCCT L409-414: VList over TopoDS_Iterator(E)
                    // (arch. diff. #56: bat::sub_shapes — the composed-children
                    // walk; the null filter carries the rcad nullified
                    // extremities).
                    let v_list: Vec<Shape> =
                        bat::sub_shapes(&e).into_iter().filter(|v| !v.is_null()).collect();
                    // OCCT L415: for (itlv over VList).
                    for v in &v_list {
                        // OCCT L418: if (myVerVerMap.IsBound(V)).
                        if shape_data_map::is_bound(&self.my_ver_ver_map, v) {
                            // OCCT L420: TopoDS_Shape NewV = myVerVerMap(V).
                            let mut new_v = shape_data_map::find(&self.my_ver_ver_map, v);
                            // OCCT L421: E.Free(true).
                            bat::builder_set_free(&mut e, true);
                            // OCCT L422: NewV.Orientation(V.Orientation()).
                            new_v.orientation = v.orientation;
                            // OCCT L423-428: the TV/NewTV tolerance merge and
                            // the points list append (the handle forms ->
                            // the vertex payloads, arch. diff. #4).
                            let (a_tv_tol, a_tv_points) = match v.data.as_ref() {
                                rcad_kernel::topo::topods::TShape::Vertex(a_tv) => {
                                    (a_tv.tolerance, a_tv.points.clone())
                                }
                                _ => (0.0, Vec::new()),
                            };
                            if let rcad_kernel::topo::topods::TShape::Vertex(a_new_tv) =
                                std::sync::Arc::make_mut(&mut new_v.data)
                            {
                                if a_tv_tol > a_new_tv.tolerance {
                                    a_new_tv.tolerance = a_tv_tol;
                                }
                                a_new_tv.points.extend(a_tv_points);
                            }
                            // OCCT L429: AsDes->Replace(V, NewV).
                            as_des.replace(v, &new_v);
                            // OCCT L430-431: BB.Remove(E, V); BB.Add(E, NewV).
                            bat::builder_remove_edge_vertex(&mut e, v);
                            bat::builder_add_edge_vertex(&mut e, &new_v);
                        }
                    }
                }
            }
        }
    }

    /// OCCT BRepOffset_MakeLoops::BuildFaces(LF, AsDes, Image, theRange)
    /// (cxx L455-660).
    pub fn build_faces(
        &mut self,
        lof: &Vec<Shape>,
        as_des: &mut BRepAlgoAsDes,
        image_offset: &mut BRepAlgoImage,
    ) {
        // OCCT L459-462: the iterators, ToRebuild, Loops +
        // Loops.VerticesForSubstitute(myVerVerMap) (arch. diff. #40: the
        // progress scope is flattened).
        let mut a_loops = BRepAlgoLoop::new();
        a_loops.vertices_for_substitute(&self.my_ver_ver_map);

        //----------------------------------
        // Loop on all faces //.
        //----------------------------------
        // OCCT L469: for (itr.Initialize(LF); itr.More(); itr.Next(), aPS.Next()).
        for itr in lof.iter() {
            // OCCT L474: TopoDS_Face F = TopoDS::Face(itr.Value()).
            let f = itr;
            // OCCT L475: Loops.Init(F).
            a_loops.init(f);
            // OCCT L476: ToRebuild = false.
            let mut to_rebuild = false;
            // OCCT L477: NCollection_List<TopoDS_Shape> AddedEdges.
            let mut added_edges: Vec<Shape> = Vec::new();

            // OCCT L479: if (!Image.HasImage(F)).
            if !image_offset.has_image(f) {
                //----------------------------------
                // Face F not yet reconstructed.
                //----------------------------------
                // OCCT L483: const NCollection_List& LE = AsDes->Descendant(F).
                let le: Vec<Shape> = as_des.descendant(f).to_vec();
                //----------------------------------------------------------------
                // first loop to find if the edges of the face were reconstructed.
                // - maj on map MONV. Some vertices on reconstructed edges
                // coincide geometrically with old but are not IsSame.
                //----------------------------------------------------------------
                // OCCT L488-489: the MONV datamap; OV1, OV2, NV1, NV2
                // (the vertex pairs are loop-local registers).
                let mut monv: DataMapOfShapeShape = HashMap::new();

                // OCCT L491: for (itl.Initialize(LE); itl.More(); itl.Next()).
                for itl in &le {
                    // OCCT L493: TopoDS_Edge E = TopoDS::Edge(itl.Value()).
                    let e = itl;
                    // OCCT L494: if (Image.HasImage(E)).
                    if image_offset.has_image(e) {
                        // OCCT L496: const NCollection_List& LCE = Image.Image(E).
                        let lce = image_offset.image(e);
                        // OCCT L497: if (LCE.Extent() == 1 && LCE.First().IsSame(E)).
                        if lce.len() == 1 && lce[0].is_same(e) {
                            // OCCT L498-502: CE =
                            // TopoDS::Edge(LCE.First().Oriented(E.Orientation())).
                            let ce = oriented(&lce[0], e.orientation);
                            // OCCT L504: Loops.AddConstEdge(CE); continue.
                            a_loops.add_const_edge(&ce);
                            continue;
                        }
                        //----------------------------------
                        // F should be reconstructed.
                        //----------------------------------
                        // OCCT L508: ToRebuild = true.
                        to_rebuild = true;
                        // OCCT L509: for (itLCE.Initialize(LCE); ...).
                        for it_lce in &lce {
                            // OCCT L511: CE =
                            // TopoDS::Edge(itLCE.Value().Oriented(E.Orientation())).
                            let ce = oriented(it_lce, e.orientation);
                            // OCCT L513-514: TopExp::Vertices(E, OV1, OV2);
                            // TopExp::Vertices(CE, NV1, NV2) (CumOri=false).
                            let (o_v1, o_v2) = top_exp_vertices(e);
                            let (n_v1, n_v2) = top_exp_vertices(&ce);
                            // OCCT L515-520: the MONV bindings of the
                            // substituted extremities.
                            if !o_v1.is_same(&n_v1) {
                                shape_data_map::bind(&mut monv, &o_v1, n_v1.clone());
                            }
                            if !o_v2.is_same(&n_v2) {
                                shape_data_map::bind(&mut monv, &o_v2, n_v2.clone());
                            }
                            // OCCT L521: Loops.AddConstEdge(CE).
                            a_loops.add_const_edge(&ce);
                        }
                    }
                }
                // OCCT L524: if (ToRebuild).
                if to_rebuild {
                    //-----------------------------------------------------------
                    // Non-reconstructed edges on other faces are added.
                    // If their vertices were reconstructed they are reconstructed.
                    //-----------------------------------------------------------
                    // OCCT L531: for (itl.Initialize(LE); ...).
                    for itl in &le {
                        // OCCT L533-534: double f, l; BRep_Tool::Range(E, f, l).
                        let (a_f, a_l) = bat::brep_tool_range(itl);
                        // OCCT L536: if (!Image.HasImage(E)).
                        if !image_offset.has_image(itl) {
                            // OCCT L537: TopExp::Vertices(E, OV1, OV2).
                            let (o_v1, o_v2) = top_exp_vertices(itl);
                            // OCCT L538: NCollection_List<TopoDS_Shape> LV.
                            let mut lv: Vec<Shape> = Vec::new();
                            // OCCT L539: if (MONV.IsBound(OV1)).
                            if shape_data_map::is_bound(&monv, &o_v1) {
                                // OCCT L541-543: VV = TopoDS::Vertex(MONV(OV1));
                                // VV.Orientation(TopAbs_FORWARD); LV.Append(VV).
                                let mut vv = shape_data_map::find(&monv, &o_v1);
                                vv.orientation = Orientation::Forward;
                                // OCCT L546-547: B.UpdateVertex(
                                // TopoDS::Vertex(VV.Oriented(TopAbs_INTERNAL)),
                                // f, E, BRep_Tool::Tolerance(VV)).
                                builder_update_vertex_param(
                                    itl,
                                    &oriented(&vv, Orientation::Internal),
                                    a_f,
                                    brep_tool_tolerance(&vv),
                                );
                                lv.push(vv);
                            }
                            // OCCT L549: if (MONV.IsBound(OV2)).
                            if shape_data_map::is_bound(&monv, &o_v2) {
                                // OCCT L551-553: VV = TopoDS::Vertex(MONV(OV2));
                                // VV.Orientation(TopAbs_REVERSED); LV.Append(VV).
                                let mut vv = shape_data_map::find(&monv, &o_v2);
                                vv.orientation = Orientation::Reversed;
                                // OCCT L556-557: B.UpdateVertex(
                                // TopoDS::Vertex(VV.Oriented(TopAbs_INTERNAL)),
                                // l, E, BRep_Tool::Tolerance(VV)).
                                builder_update_vertex_param(
                                    itl,
                                    &oriented(&vv, Orientation::Internal),
                                    a_l,
                                    brep_tool_tolerance(&vv),
                                );
                                lv.push(vv);
                            }
                            // OCCT L562-570: LV.IsEmpty() ?
                            // Loops.AddConstEdge(E) : Loops.AddEdge(E, LV) +
                            // AddedEdges.Append(E).
                            if lv.is_empty() {
                                a_loops.add_const_edge(itl);
                            } else {
                                a_loops.add_edge(itl, &lv);
                                added_edges.push(itl.clone());
                            }
                        }
                    }
                }
            }
            // OCCT L574: if (ToRebuild).
            if to_rebuild {
                //------------------------
                // Reconstruction.
                //------------------------
                // OCCT L577-578: Loops.Perform(); Loops.WiresToFaces().
                a_loops.perform();
                a_loops.wires_to_faces();
                //------------------------
                // MAJ SD.
                //------------------------
                // OCCT L583: const NCollection_List& NF = Loops.NewFaces().
                let nf = a_loops.new_faces().to_vec();
                //-----------------------
                // F => New faces;
                //-----------------------
                // OCCT L586: Image.Bind(F, NF).
                image_offset.bind_list(f, &nf);

                // OCCT L588: for (itAdded.Initialize(AddedEdges); ...).
                for it_added in &added_edges {
                    // OCCT L591: const TopoDS_Edge& E =
                    // TopoDS::Edge(itAdded.Value()).
                    let e = it_added;
                    //-----------------------
                    //  E => New edges;
                    //-----------------------
                    // OCCT L596-602: HasImage ? Image.Add : Image.Bind.
                    if image_offset.has_image(e) {
                        image_offset.add_list(e, a_loops.new_edges(e));
                    } else {
                        image_offset.bind_list(e, a_loops.new_edges(e));
                    }
                }
            }
        }
        // OCCT L607: Loops.GetVerticesForSubstitute(myVerVerMap).
        a_loops.get_vertices_for_substitute(&mut self.my_ver_ver_map);
        // OCCT L608-611: if (myVerVerMap.IsEmpty()) return.
        if self.my_ver_ver_map.is_empty() {
            return;
        }
        // OCCT L612-660: the vertex substitution tail — the OCCT source
        // repeats the BuildOnContext block verbatim; translated in place.
        // OCCT L614: for (itr.Initialize(LF); itr.More(); itr.Next()).
        for itr in lof.iter() {
            // OCCT L616: TopoDS_Shape F = itr.Value().
            let f = itr;
            // OCCT L617-618: NCollection_List LIF; Image.LastImage(F, LIF).
            let mut l_if: Vec<Shape> = Vec::new();
            image_offset.last_image(f, &mut l_if);
            // OCCT L619: for (itl.Initialize(LIF); itl.More(); itl.Next()).
            for itl in &l_if {
                // OCCT L621: const TopoDS_Shape& IF = itl.Value().
                let i_f = itl;
                // OCCT L622: TopExp_Explorer EdExp(IF, TopAbs_EDGE).
                for ed_exp in bat::explorer(i_f, ShapeType::Edge, ShapeType::Shape) {
                    // OCCT L624: TopoDS_Shape E = EdExp.Current().
                    let mut e = ed_exp;
                    // OCCT L625-630: VList over TopoDS_Iterator(E)
                    // (arch. diff. #56: bat::sub_shapes; the null filter
                    // carries the rcad nullified extremities).
                    let v_list: Vec<Shape> =
                        bat::sub_shapes(&e).into_iter().filter(|v| !v.is_null()).collect();
                    // OCCT L631: for (itlv over VList).
                    for v in &v_list {
                        // OCCT L634: if (myVerVerMap.IsBound(V)).
                        if shape_data_map::is_bound(&self.my_ver_ver_map, v) {
                            // OCCT L636: TopoDS_Shape NewV = myVerVerMap(V).
                            let mut new_v = shape_data_map::find(&self.my_ver_ver_map, v);
                            // OCCT L637: E.Free(true).
                            bat::builder_set_free(&mut e, true);
                            // OCCT L638: NewV.Orientation(V.Orientation()).
                            new_v.orientation = v.orientation;
                            // OCCT L639-644: the TV/NewTV tolerance merge and
                            // the points list append (arch. diff. #4).
                            let (a_tv_tol, a_tv_points) = match v.data.as_ref() {
                                rcad_kernel::topo::topods::TShape::Vertex(a_tv) => {
                                    (a_tv.tolerance, a_tv.points.clone())
                                }
                                _ => (0.0, Vec::new()),
                            };
                            if let rcad_kernel::topo::topods::TShape::Vertex(a_new_tv) =
                                std::sync::Arc::make_mut(&mut new_v.data)
                            {
                                if a_tv_tol > a_new_tv.tolerance {
                                    a_new_tv.tolerance = a_tv_tol;
                                }
                                a_new_tv.points.extend(a_tv_points);
                            }
                            // OCCT L645: AsDes->Replace(V, NewV).
                            as_des.replace(v, &new_v);
                            // OCCT L646-647: BB.Remove(E, V); BB.Add(E, NewV).
                            bat::builder_remove_edge_vertex(&mut e, v);
                            bat::builder_add_edge_vertex(&mut e, &new_v);
                        }
                    }
                }
            }
        }
    }
}

/// OCCT IsBetweenCorks(E, AsDes, LContext) (BRepOffset_MakeLoops.cxx
/// L196-213) — the file static: E lies between corks when every ascendant
/// face of E is in LContext (no ascendant counts as between).
fn is_between_corks(e: &Shape, as_des: &BRepAlgoAsDes, l_context: &[Shape]) -> bool {
    // OCCT L197-200: if (!AsDes->HasAscendant(E)) return true.
    if !as_des.has_ascendant(e) {
        return true;
    }
    // OCCT L201: const NCollection_List& LF = AsDes->Ascendant(E).
    let lf = as_des.ascendant(e);
    // OCCT L202-214: every S of LF must be IsSame to a member of LContext.
    for s in lf {
        let mut found = false;
        for it2 in l_context {
            if s.is_same(it2) {
                found = true;
                break;
            }
        }
        if !found {
            return false;
        }
    }
    true
}

/// OCCT BRep_Builder::UpdateVertex(V, P, E, Tol) (BRep_Builder.cxx
/// L1220-1313) — the BuildFaces call passes TopoDS::Vertex(
/// VV.Oriented(TopAbs_INTERNAL)): the copy shares the new-vertex TShape,
/// which is not IsSame to the edge extremities, so the builder search ends
/// INTERNAL and the parameter lands in the edge's vertex parameter map (the
/// kernel BRepBuilder::update_vertex_on_edge INTERNAL branch), followed by
/// the TV->UpdateTolerance(Tol) max-update on the vertex payload.
/// Arch. diff. #4: the OCCT edit is handle-shared; the rcad form mutates
/// the local Arc payload.
fn builder_update_vertex_param(the_e: &Shape, the_v: &Shape, par: f64, the_tol: f64) {
    let mut a_edata = std::sync::Arc::clone(&the_e.data);
    if let rcad_kernel::topo::topods::TShape::Edge(a_ed) = std::sync::Arc::make_mut(&mut a_edata) {
        a_ed.vertex_params.insert(the_v.ptr_id(), par);
    }
    let mut a_vdata = std::sync::Arc::clone(&the_v.data);
    if let rcad_kernel::topo::topods::TShape::Vertex(a_vd) = std::sync::Arc::make_mut(&mut a_vdata)
    {
        a_vd.tolerance = a_vd.tolerance.max(the_tol);
    }
}

