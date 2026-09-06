// OCCT HLRTopoBRep_OutLiner (TKHLR/HLRTopoBRep/HLRTopoBRep_OutLiner.hxx
// L1-77 + .cxx L1-344 + .lxx L1-51) — builds the outLined shape: on every
// face of the original shape the internal outlines (IntL) and isolines
// (IsoL) recorded in the Data structure are added as new INTERNAL wires,
// split edges (SplE) replace their originals, and the result is collected
// into a compound of shells (BuildShape).
//
// Translation notes:
// - Standard_Transient base class (hxx L33): rcad keeps a plain struct; the
//   handle<HLRTopoBRep_OutLiner> semantics land with the consumer stage and
//   DEFINE_STANDARD_RTTIEXT has no rcad counterpart.
// - OCCT mutates TopoDS trees through BRep_Builder; rcad routes every
//   builder operation through the `brep: &mut BRep` arena carried by
//   fill/process_face/build_shape (the ShapeBuild_ReShape precedent) — the
//   OCCT member functions have no such parameter.
// - The MST map (NCollection_DataMap<TopoDS_Shape, BRepTopAdaptor_Tool,
//   TopTools_ShapeMapHasher>) is the contract encoding
//   Vec<(Shape, BRepTopAdaptorTool)> keyed by Shape::ptr_id (the data.rs
//   map precedent; iteration order = insertion order).  Bind copies the
//   Tool value (the two handles are shared), so the tool chain is Clone.
// - TopExp_Explorer / TopExp::MapShapesAndAncestors / TopExp::Vertices are
//   translated locally (TopExp_Explorer.cxx L34-244, TopExp.cxx L80-120 and
//   L214-251): rcad has no shared TopExp module yet and the existing
//   shhealing `topexp_explorer` helper carries no ToAvoid support.

use rcad_kernel::base::extrema::ExtPC;
use rcad_kernel::topods::{tshape_flags, BRep, BRepBuilder, Orientation, Shape, ShapeType};
use glam::DVec3;

use crate::bop::int_tools::bean_face_intersector::{BRepAdaptorCurve, GeomAbsCurveType};
use crate::hlr::algo::projector::Projector;
use crate::hlr::contap::Contour;
use crate::shhealing::shape_build::brep_tool::{builder_add, raw_subshapes, set_flag_inplace};
use crate::topalgo::brep_top_adaptor::tool::BRepTopAdaptorTool;

use super::data::Data;

/// OCCT Extrema_GGExtPC.hxx L88: the two-argument Extrema_ExtPC(P, C)
/// constructor forwards the default `theTolF = 1.0e-10`.
const EXTREMA_2ARG_TOLF: f64 = 1.0e-10;

/// OCCT TopAbs_ShapeEnum ranking (TopAbs_ShapeEnum.hxx L48-58): COMPOUND = 0
/// < ... < VERTEX = 7 < SHAPE = 8 — "more complex" is the LOWER value.  The
/// rcad ShapeType enum is declared in the opposite order, so the OCCT value
/// is mapped explicitly.
fn occt_shape_rank(t: ShapeType) -> i32 {
    match t {
        ShapeType::Compound => 0,
        ShapeType::CompSolid => 1,
        ShapeType::Solid => 2,
        ShapeType::Shell => 3,
        ShapeType::Face => 4,
        ShapeType::Wire => 5,
        ShapeType::Edge => 6,
        ShapeType::Vertex => 7,
        ShapeType::Shape => 8,
    }
}

/// OCCT NCollection_DataMap::IsBound under TopTools_ShapeMapHasher — the
/// TShape-pointer identity (the data.rs map precedent; the HLR shapes carry
/// identity locations, so the ptr_id key equals the TShape + Location key).
fn map_is_bound<V>(map: &[(Shape, V)], s: &Shape) -> bool {
    let id = s.ptr_id();
    map.iter().any(|(k, _)| k.ptr_id() == id)
}

/// OCCT NCollection_DataMap::FindIndex / FindFromKey — the position of the
/// key in the map (FindFromKey raises Standard_NoSuchObject when unbound;
/// the expect mirrors it).
fn map_find_from_key<'a>(map: &'a [(Shape, Vec<Shape>)], s: &Shape) -> &'a Vec<Shape> {
    let id = s.ptr_id();
    let pos = map
        .iter()
        .position(|(k, _)| k.ptr_id() == id)
        .expect("Standard_NoSuchObject: aVEMap key not bound");
    &map[pos].1
}

/// OCCT NCollection_IndexedDataMap::FindIndex + Add (TopExp.cxx L98-102,
/// L113-117): the position of the key, appending the key with an empty list
/// when not yet present.
fn ve_map_pos_or_add(map: &mut Vec<(Shape, Vec<Shape>)>, s: &Shape) -> usize {
    let id = s.ptr_id();
    if let Some(pos) = map.iter().position(|(k, _)| k.ptr_id() == id) {
        pos
    } else {
        map.push((s.clone(), Vec::new()));
        map.len() - 1
    }
}

/// OCCT NCollection_Map::Add — True when the key was not yet present
/// (BuildShape cxx L324 / L334).
fn shape_map_add(shape_map: &mut Vec<Shape>, s: &Shape) -> bool {
    let id = s.ptr_id();
    if shape_map.iter().any(|x| x.ptr_id() == id) {
        false
    } else {
        shape_map.push(s.clone());
        true
    }
}

/// OCCT TopExp_Explorer (TopExp_Explorer.cxx L34-244) — the stack explorer
/// with ToFind/ToAvoid the OutLiner drives (cxx L105, L131, L311).  Each
/// stack level is a TopoDS_Iterator (TopoDS_Iterator.cxx L28-70, the
/// defaults cumOri = cumLoc = true): the composed children snapshot and the
/// current position.
struct TopExpExplorer {
    /// OCCT toFind.
    to_find: ShapeType,
    /// OCCT toAvoid.
    to_avoid: ShapeType,
    /// OCCT hasMore.
    has_more: bool,
    /// OCCT myShape — Current() while the stack is empty.
    my_shape: Shape,
    /// OCCT myStack[0..=myStackTop]: (the composed children, the position).
    my_stack: Vec<(Vec<Shape>, usize)>,
}

impl TopExpExplorer {
    /// OCCT TopExp_Explorer() (cxx L39-45): toFind = toAvoid = TopAbs_SHAPE,
    /// hasMore = false.
    fn new() -> Self {
        TopExpExplorer {
            to_find: ShapeType::Shape,
            to_avoid: ShapeType::Shape,
            has_more: false,
            my_shape: Shape::null(),
            my_stack: Vec::new(),
        }
    }

    /// OCCT shouldAvoid (cxx L25-28).
    fn should_avoid(&self, ty: ShapeType) -> bool {
        self.to_avoid != ShapeType::Shape && self.to_avoid == ty
    }

    /// OCCT pushIterator (cxx L182-193) with the TopoDS_Iterator(S)
    /// construction (TopoDS_Iterator.cxx L28-70): the children carry the
    /// composed orientation (TopAbs::Compose) and location (Move).
    fn push_iterator(&mut self, brep: &mut BRep, s: &Shape) {
        let mut children = raw_subshapes(brep, s);
        for c in children.iter_mut() {
            c.orientation = s.orientation.compose(c.orientation);
            if s.location != 0 {
                // OCCT updateCurrentShape: myShape.Move(myLocation, false).
                let composed = brep.get_location(s.location) * brep.get_location(c.location);
                c.location = brep.add_location(composed);
            }
        }
        self.my_stack.push((children, 0));
    }

    /// OCCT Init (cxx L66-108).
    fn init(&mut self, brep: &mut BRep, s: &Shape, to_find: ShapeType, to_avoid: ShapeType) {
        // OCCT Clear() (cxx L195-204).
        self.my_stack.clear();
        self.has_more = false;
        self.my_shape = s.clone();
        self.to_find = to_find;
        self.to_avoid = to_avoid;

        if s.is_null() {
            self.has_more = false;
            return;
        }

        if to_find == ShapeType::Shape {
            self.has_more = false;
        } else {
            let ty = s.shape_type();
            if occt_shape_rank(ty) > occt_shape_rank(to_find) {
                self.has_more = false;
            } else if ty != to_find {
                self.has_more = true;
                self.next(brep);
            } else {
                self.has_more = true;
            }
        }
    }

    /// OCCT More (TopExp_Explorer.hxx L155-158).
    fn more(&self) -> bool {
        self.has_more
    }

    /// OCCT Current (cxx L156-161): myStackTop < 0 ? myShape :
    /// myStack[myStackTop].Value().
    fn current(&self) -> &Shape {
        if self.my_stack.is_empty() {
            &self.my_shape
        } else {
            let top = self.my_stack.last().expect("TopExp_Explorer: stack");
            &top.0[top.1]
        }
    }

    /// OCCT Next (cxx L110-154).
    fn next(&mut self, brep: &mut BRep) {
        if self.my_stack.is_empty() {
            let ty = self.my_shape.shape_type();
            if ty == self.to_find {
                self.has_more = false;
                return;
            } else if self.should_avoid(ty) {
                self.has_more = false;
                return;
            } else {
                let s = self.my_shape.clone();
                self.push_iterator(brep, &s);
            }
        } else {
            // myStack[myStackTop].Next().
            let top = self.my_stack.len() - 1;
            self.my_stack[top].1 += 1;
        }

        loop {
            let top = self.my_stack.len() - 1;
            if self.my_stack[top].1 < self.my_stack[top].0.len() {
                let shap_top = self.my_stack[top].0[self.my_stack[top].1].clone();
                let ty = shap_top.shape_type();
                if ty == self.to_find {
                    self.has_more = true;
                    return;
                } else if occt_shape_rank(self.to_find) > occt_shape_rank(ty)
                    && !self.should_avoid(ty)
                {
                    // pushIterator(TopoDS_Iterator(aShapTop)) — aTopIter
                    // reference invalid after the push (cxx L152-155).
                    self.push_iterator(brep, &shap_top);
                } else {
                    self.my_stack[top].1 += 1;
                }
            } else {
                // popIterator (cxx L196-201).
                self.my_stack.pop();
                if self.my_stack.is_empty() {
                    break;
                }
                let top = self.my_stack.len() - 1;
                self.my_stack[top].1 += 1;
            }
        }
        self.has_more = false;
    }
}

/// OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L80-120) — the indexed
/// (sub-shape -> ancestors) map; rcad encodes the
/// NCollection_IndexedDataMap as Vec<(Shape, Vec<Shape>)> (insertion order,
/// ptr_id keys).
fn map_shapes_and_ancestors(
    brep: &mut BRep,
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut Vec<(Shape, Vec<Shape>)>,
) {
    // visit ancestors (cxx L89-107).
    let mut exa = TopExpExplorer::new();
    exa.init(brep, s, ta, ShapeType::Shape);
    while exa.more() {
        // const TopoDS_Shape& anc = exa.Current();
        let anc = exa.current().clone();
        let mut exs = TopExpExplorer::new();
        exs.init(brep, &anc, ts, ShapeType::Shape);
        while exs.more() {
            // int index = M.FindIndex(exs.Current()); if (index == 0)
            // index = M.Add(exs.Current(), empty); M(index).Append(anc);
            let cur = exs.current().clone();
            let pos = ve_map_pos_or_add(m, &cur);
            m[pos].1.push(anc.clone());
            exs.next(brep);
        }
        exa.next(brep);
    }

    // visit shapes not under ancestors (cxx L109-119).
    let mut ex = TopExpExplorer::new();
    ex.init(brep, s, ts, ta);
    while ex.more() {
        // int index = M.FindIndex(ex.Current()); if (index == 0)
        // index = M.Add(ex.Current(), empty);
        let cur = ex.current().clone();
        ve_map_pos_or_add(m, &cur);
        ex.next(brep);
    }
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast) (TopExp.cxx L214-251, the default
/// CumOri = false): the FORWARD child is Vfirst, the REVERSED child is Vlast,
/// both stay null when absent (TopoDS_Iterator(E, false) keeps the stored
/// child orientations; the rcad TEdgeData stores them as first / last in
/// storage order).  OCCT Nullify() at the end when undefined — the helper
/// returns the null shape for that case.
fn top_exp_vertices(e: &Shape) -> (Shape, Shape) {
    let mut v_first = Shape::null();
    let mut v_last = Shape::null();
    if let Some(ed) = e.as_edge() {
        for v in [&ed.first, &ed.last] {
            if v.orientation == Orientation::Forward {
                v_first = v.clone();
            } else if v.orientation == Orientation::Reversed {
                v_last = v.clone();
            }
        }
    }
    (v_first, v_last)
}

/// OCCT BRepAdaptor_Curve(E) (BRepAdaptor_Curve.hxx L44-52) — the edge
/// adaptor over the edge's 3D curve and the [First, Last] parameter range.
/// rcad maps it onto BRepAdaptorCurve::with_range over the TEdgeData curve
/// + range (the established bop/int_tools/context.rs encoding); an edge
/// without a 3D curve raises in OCCT on first use, the expect mirrors it
/// (the IntL / IsoL edges always carry 3D curves).
fn brep_adaptor_curve(e: &Shape) -> BRepAdaptorCurve {
    let ed = e
        .as_edge()
        .expect("Standard_NoSuchObject: BRepAdaptor_Curve: not an edge");
    let curve = ed
        .curve
        .clone()
        .expect("StdFail_NotDone: BRepAdaptor_Curve: edge without 3D curve");
    BRepAdaptorCurve::with_range(curve, ed.range[0], ed.range[1])
}

/// OCCT HLRTopoBRep_OutLiner (hxx L33-73).
pub struct OutLiner {
    /// OCCT myOriginalShape (hxx L70).
    my_original_shape: Shape,
    /// OCCT myOutLinedShape (hxx L71).
    my_out_lined_shape: Shape,
    /// OCCT myDS (hxx L72) — the value member.
    my_ds: Data,
}

impl OutLiner {
    /// OCCT HLRTopoBRep_OutLiner() (cxx L44).
    pub fn new() -> Self {
        OutLiner {
            my_original_shape: Shape::null(),
            my_out_lined_shape: Shape::null(),
            my_ds: Data::new(),
        }
    }

    /// OCCT HLRTopoBRep_OutLiner(const TopoDS_Shape& OriS) (cxx L48-51).
    pub fn with_original_shape(ori_s: &Shape) -> Self {
        OutLiner {
            my_original_shape: ori_s.clone(),
            my_out_lined_shape: Shape::null(),
            my_ds: Data::new(),
        }
    }

    /// OCCT HLRTopoBRep_OutLiner(const TopoDS_Shape& OriS, const TopoDS_Shape&
    /// OutS) (cxx L55-59).
    pub fn with_original_and_outlined(ori_s: &Shape, out_s: &Shape) -> Self {
        OutLiner {
            my_original_shape: ori_s.clone(),
            my_out_lined_shape: out_s.clone(),
            my_ds: Data::new(),
        }
    }

    /// OCCT OriginalShape(const TopoDS_Shape& OriS) (lxx L19-22).
    pub fn set_original_shape(&mut self, ori_s: &Shape) {
        self.my_original_shape = ori_s.clone();
    }

    /// OCCT TopoDS_Shape& OriginalShape() (lxx L26-29).
    pub fn original_shape(&mut self) -> &mut Shape {
        &mut self.my_original_shape
    }

    /// OCCT OutLinedShape(const TopoDS_Shape& OutS) (lxx L33-36).
    pub fn set_out_lined_shape(&mut self, out_s: &Shape) {
        self.my_out_lined_shape = out_s.clone();
    }

    /// OCCT TopoDS_Shape& OutLinedShape() (lxx L40-43).
    pub fn out_lined_shape(&mut self) -> &mut Shape {
        &mut self.my_out_lined_shape
    }

    /// OCCT HLRTopoBRep_Data& DataStructure() (lxx L47-50).
    pub fn data_structure(&mut self) -> &mut Data {
        &mut self.my_ds
    }

    /// OCCT Fill (cxx L63-92).  The `brep` arena is the rcad builder surface
    /// (see the module notes).
    pub fn fill(
        &mut self,
        brep: &mut BRep,
        p: &Projector,
        mst: &mut Vec<(Shape, BRepTopAdaptorTool)>,
        nb_iso: usize,
    ) {
        if !self.my_original_shape.is_null() {
            if self.my_out_lined_shape.is_null() {
                // gp_Vec Vecz(0., 0., 1.);
                let mut vecz = DVec3::new(0.0, 0.0, 1.0);
                // gp_Trsf Tr(P.Transformation());
                let mut tr = *p.transformation();
                // Tr.Invert();
                tr.invert();
                // Vecz.Transform(Tr);
                vecz = tr.transform_vec(vecz);
                // Contap_Contour FO;
                let mut fo = Contour::new();
                if p.perspective() {
                    // gp_Pnt Eye; Eye.SetXYZ(P.Focus() * Vecz.XYZ());
                    let eye = p.focus() * vecz;
                    // FO.Init(Eye);
                    fo.init_eye(eye);
                } else {
                    // gp_Dir DirZ(Vecz);
                    // FO.Init(DirZ);
                    fo.init(vecz);
                }
                // HLRTopoBRep_DSFiller::Insert(myOriginalShape, FO, myDS, MST, nbIso);
                // (the landed rcad Insert carries the leading `brep` builder
                // arena — see the module notes.)
                super::ds_filler::insert(brep, &self.my_original_shape, &mut fo, &mut self.my_ds, mst, nb_iso);
                // BuildShape(MST);
                self.build_shape(brep, mst);
            }
        }
    }

    // ====================================================================
    // OCCT ProcessFace (cxx L94-304) — build a Face using myDS and add the
    // new face to a shell.
    // ====================================================================
    fn process_face(
        &mut self,
        brep: &mut BRep,
        f: &Shape,
        s: &Shape,
        mst: &mut Vec<(Shape, BRepTopAdaptorTool)>,
    ) {
        // OCCT L104: BRep_Builder B.
        let mut b = BRepBuilder::new();
        // OCCT L105: TopExp_Explorer exE, exW.
        let mut ex_e = TopExpExplorer::new();
        let mut ex_w = TopExpExplorer::new();
        // OCCT L106: // bool splitted = false;

        // OCCT L108-109: NCollection_IndexedDataMap<TopoDS_Shape,
        // NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher> aVEMap.
        let mut a_ve_map: Vec<(Shape, Vec<Shape>)> = Vec::new();
        // OCCT L110: TopExp::MapShapesAndAncestors(F, TopAbs_VERTEX,
        // TopAbs_EDGE, aVEMap).
        map_shapes_and_ancestors(brep, f, ShapeType::Vertex, ShapeType::Edge, &mut a_ve_map);

        // OCCT L112: TopoDS_Shape NF; (null declaration) —
        // OCCT L124: NF = F.EmptyCopied(); (the null TShape is replaced at
        // once, the copy is unconditional after the commented L115-123 scan).
        let nf = brep.empty_copied(f);

        // OCCT L126: for (exW.Init(F, TopAbs_WIRE); exW.More(); exW.Next())
        ex_w.init(brep, f, ShapeType::Wire, ShapeType::Shape);
        while ex_w.more() {
            // OCCT L128-129: TopoDS_Wire W; B.MakeWire(W);
            let w = b.make_wire(brep);

            // OCCT L131: for (exE.Init(exW.Current(), TopAbs_EDGE); exE.More();
            // exE.Next())
            ex_e.init(brep, ex_w.current(), ShapeType::Edge, ShapeType::Shape);
            while ex_e.more() {
                // OCCT L133: TopoDS_Edge E = TopoDS::Edge(exE.Current());
                let e = ex_e.current().clone();
                // OCCT L134: if (myDS.EdgeHasSplE(E))
                if self.my_ds.edge_has_spl_e(&e) {
                    // OCCT L137-138: for (itS.Initialize(myDS.EdgeSplE(E));
                    // itS.More(); itS.Next())
                    let mut it_s = 0;
                    while it_s < self.my_ds.edge_spl_e(&e).len() {
                        // OCCT L140: TopoDS_Edge newE = TopoDS::Edge(itS.Value());
                        let mut new_e = self.my_ds.edge_spl_e(&e)[it_s].clone();
                        // OCCT L141: newE.Orientation(E.Orientation());
                        new_e.orientation = e.orientation;
                        // OCCT L142: myDS.AddOldS(newE, E);
                        self.my_ds.add_old_s(&new_e, &e);
                        // OCCT L143: B.Add(W, newE);
                        builder_add(brep, &w, &new_e);
                        it_s += 1;
                    }
                } else {
                    // OCCT L148: B.Add(W, E);
                    builder_add(brep, &w, &e);
                }
                // exE.Next().
                ex_e.next(brep);
            }
            // OCCT L151: B.Add(NF, W); — add the new wire in the new face.
            builder_add(brep, &nf, &w);
            // exW.Next().
            ex_w.next(brep);
        }

        // OCCT L154: myDS.AddIntL(F);
        self.my_ds.add_int_l(f);
        // OCCT L155: NCollection_List<TopoDS_Shape>& OutL = myDS.AddOutL(F);
        // (the OCCT map-owned reference is re-derived at each Append site
        // below — Rust cannot hold it across the myDS mutations in between;
        // AddOutL re-binds nothing, F is already bound by AddIntL.)
        self.my_ds.add_out_l(f);

        // OCCT L157: if (myDS.FaceHasIntL(F)) — get the InternalOutLines on
        // face F.
        if self.my_ds.face_has_int_l(f) {
            // OCCT L159: TopoDS_Wire W; (null until MakeWire).
            let mut w = Shape::null();

            // OCCT L161-162: for (itE.Initialize(myDS.FaceIntL(F)); itE.More();
            // itE.Next())
            let mut it_e = 0;
            while it_e < self.my_ds.face_int_l(f).len() {
                // OCCT L164: TopoDS_Edge E = TopoDS::Edge(itE.Value());
                let mut e = self.my_ds.face_int_l(f)[it_e].clone();
                // OCCT L165: E.Orientation(TopAbs_INTERNAL);
                e.orientation = Orientation::Internal;
                // OCCT L166: // Check, if outline edge coincides real edge

                // OCCT L168: BRepAdaptor_Curve C(E);
                let c = brep_adaptor_curve(&e);
                // OCCT L169: double par = 0.34 * C.FirstParameter() + 0.66 *
                // C.LastParameter();
                let par = 0.34 * c.first_parameter() + 0.66 * c.last_parameter();
                // OCCT L170: gp_Pnt P = C.Value(par);
                let p = c.value(par);
                // OCCT L171-172: TopoDS_Vertex V1, V2, aV1, aV2;
                // TopExp::Vertices(E, V1, V2); (the aV1 / aV2 pair is declared
                // before the loop in OCCT and fully rewritten by each
                // TopExp::Vertices call — the initial null is never read, so
                // the rcad bindings live at the call sites).
                let (v1, v2) = top_exp_vertices(&e);

                // OCCT L174: bool SameEdge = false;
                let mut same_edge = false;
                // OCCT L175: if (!V1.IsNull() && aVEMap.Contains(V1))
                if !v1.is_null() && map_is_bound(&a_ve_map, &v1) {
                    // OCCT L177: const NCollection_List<TopoDS_Shape>& aEList
                    // = aVEMap.FindFromKey(V1);
                    let a_e_list = map_find_from_key(&a_ve_map, &v1);
                    // OCCT L178-179: NCollection_List<TopoDS_Shape>::Iterator
                    // it(aEList); for (; it.More(); it.Next())
                    let mut it = 0;
                    while it < a_e_list.len() {
                        // OCCT L181: const TopoDS_Edge& aE = TopoDS::Edge(it.Value());
                        let a_e = a_e_list[it].clone();
                        // OCCT L182: TopExp::Vertices(aE, aV1, aV2);
                        let (a_v1, a_v2) = top_exp_vertices(&a_e);

                        // OCCT L184: if ((V1.IsSame(aV1) && V2.IsSame(aV2)) ||
                        // (V1.IsSame(aV2) && V2.IsSame(aV1)))
                        if (v1.is_same(&a_v1) && v2.is_same(&a_v2))
                            || (v1.is_same(&a_v2) && v2.is_same(&a_v1))
                        {
                            // OCCT L186: BRepAdaptor_Curve aC(aE);
                            let a_c = brep_adaptor_curve(&a_e);
                            // OCCT L187: if ((C.GetType() == GeomAbs_Line) &&
                            // (aC.GetType() == GeomAbs_Line))
                            if c.get_type() == GeomAbsCurveType::Line
                                && a_c.get_type() == GeomAbsCurveType::Line
                            {
                                // OCCT L189-190: SameEdge = true; break;
                                same_edge = true;
                                break;
                            } else {
                                // OCCT L194: // Try to project one point
                                // OCCT L195: Extrema_ExtPC anExt(P, aC);
                                let an_ext = ExtPC::new(
                                    p,
                                    a_c.curve(),
                                    EXTREMA_2ARG_TOLF,
                                    a_c.first_parameter(),
                                    a_c.last_parameter(),
                                );
                                // OCCT L196: if (anExt.IsDone())
                                if an_ext.is_done() {
                                    // OCCT L198: int aNe = anExt.NbExt();
                                    let a_ne = an_ext.nb_ext();
                                    // OCCT L199: if (aNe > 0)
                                    if a_ne > 0 {
                                        // OCCT L201: double dist = RealLast();
                                        let mut dist = f64::MAX;
                                        // OCCT L202-203: for (ec = 1; ec <= aNe; ++ec)
                                        let mut ec = 1;
                                        while ec <= a_ne {
                                            // OCCT L206: dist = std::min(dist,
                                            // anExt.SquareDistance(ec));
                                            dist = dist.min(an_ext.square_distance(ec));
                                            // OCCT L210: if (dist <= 1.e-14)
                                            if dist <= 1.0e-14 {
                                                // OCCT L212-213: SameEdge = true; break;
                                                same_edge = true;
                                                break;
                                            }
                                            ec += 1;
                                        }
                                    }
                                }
                            }
                        }
                        it += 1;
                    }
                }

                // OCCT L222-226: if (SameEdge) { OutL.Append(E); continue; }
                if same_edge {
                    self.my_ds.add_out_l(f).push(e.clone());
                    it_e += 1;
                    continue;
                }

                // OCCT L228: if (myDS.EdgeHasSplE(E))
                if self.my_ds.edge_has_spl_e(&e) {
                    // OCCT L231-232: for (itS.Initialize(myDS.EdgeSplE(E));
                    // itS.More(); itS.Next())
                    let mut it_s = 0;
                    while it_s < self.my_ds.edge_spl_e(&e).len() {
                        // OCCT L234: TopoDS_Shape newE = itS.Value();
                        let mut new_e = self.my_ds.edge_spl_e(&e)[it_s].clone();
                        // OCCT L235: newE.Orientation(TopAbs_INTERNAL);
                        new_e.orientation = Orientation::Internal;
                        // OCCT L236-239: if (W.IsNull()) { B.MakeWire(W); }
                        if w.is_null() {
                            w = b.make_wire(brep);
                        }
                        // OCCT L240: myDS.AddOldS(newE, F);
                        self.my_ds.add_old_s(&new_e, f);
                        // OCCT L241: B.Add(W, newE);
                        builder_add(brep, &w, &new_e);
                        it_s += 1;
                    }
                } else {
                    // OCCT L246-249: if (W.IsNull()) { B.MakeWire(W); }
                    if w.is_null() {
                        w = b.make_wire(brep);
                    }
                    // OCCT L250: myDS.AddOldS(E, F);
                    self.my_ds.add_old_s(&e, f);
                    // OCCT L251: B.Add(W, E);
                    builder_add(brep, &w, &e);
                }
                // itE.Next().
                it_e += 1;
            }
            // OCCT L254-257: if (!W.IsNull()) { B.Add(NF, W); } — add the new
            // wire in the new face.
            if !w.is_null() {
                builder_add(brep, &nf, &w);
            }
        }

        // OCCT L260: if (myDS.FaceHasIsoL(F)) — get the IsoLines on face F.
        if self.my_ds.face_has_iso_l(f) {
            // OCCT L262: TopoDS_Wire W;
            let mut w = Shape::null();

            // OCCT L264-265: for (itE.Initialize(myDS.FaceIsoL(F)); itE.More();
            // itE.Next())
            let mut it_e = 0;
            while it_e < self.my_ds.face_iso_l(f).len() {
                // OCCT L267: TopoDS_Edge E = TopoDS::Edge(itE.Value());
                let mut e = self.my_ds.face_iso_l(f)[it_e].clone();
                // OCCT L268: E.Orientation(TopAbs_INTERNAL);
                e.orientation = Orientation::Internal;
                // OCCT L269: if (myDS.EdgeHasSplE(E)) — normally IsoLines are
                // never split.
                if self.my_ds.edge_has_spl_e(&e) {
                    // OCCT L272-273: for (itS.Initialize(myDS.EdgeSplE(E)); ...)
                    let mut it_s = 0;
                    while it_s < self.my_ds.edge_spl_e(&e).len() {
                        // OCCT L275: TopoDS_Shape newE = itS.Value();
                        let mut new_e = self.my_ds.edge_spl_e(&e)[it_s].clone();
                        // OCCT L276: newE.Orientation(TopAbs_INTERNAL);
                        new_e.orientation = Orientation::Internal;
                        // OCCT L277-280: if (W.IsNull()) { B.MakeWire(W); }
                        if w.is_null() {
                            w = b.make_wire(brep);
                        }
                        // OCCT L281: myDS.AddOldS(newE, F);
                        self.my_ds.add_old_s(&new_e, f);
                        // OCCT L282: B.Add(W, newE);
                        builder_add(brep, &w, &new_e);
                        it_s += 1;
                    }
                } else {
                    // OCCT L287-290: if (W.IsNull()) { B.MakeWire(W); }
                    if w.is_null() {
                        w = b.make_wire(brep);
                    }
                    // OCCT L291: myDS.AddOldS(E, F);
                    self.my_ds.add_old_s(&e, f);
                    // OCCT L292: B.Add(W, E);
                    builder_add(brep, &w, &e);
                }
                // itE.Next().
                it_e += 1;
            }
            // OCCT L295-298: if (!W.IsNull()) { B.Add(NF, W); } — add the new
            // wire in the new face.
            if !w.is_null() {
                builder_add(brep, &nf, &w);
            }
        }

        // OCCT L300: myDS.AddOldS(NF, F);
        self.my_ds.add_old_s(&nf, f);
        // OCCT L301: MST.Bind(NF, MST.ChangeFind(F)); — the Tool value copy
        // (the two handles are shared): the rcad Clone.
        let pos_f = mst
            .iter()
            .position(|(k, _)| k.ptr_id() == f.ptr_id())
            .expect("Standard_NoSuchObject: MST.ChangeFind(F)");
        let tool = mst[pos_f].1.clone();
        mst.push((nf.clone(), tool));

        // OCCT L303: B.Add(S, NF); — add the face in the shell.
        builder_add(brep, s, &nf);
    }

    // ====================================================================
    // OCCT BuildShape (cxx L306-344).
    // ====================================================================
    fn build_shape(
        &mut self,
        brep: &mut BRep,
        mst: &mut Vec<(Shape, BRepTopAdaptorTool)>,
    ) {
        // OCCT L311: TopExp_Explorer exshell, exface, exedge;
        let mut ex_shell = TopExpExplorer::new();
        let mut ex_face = TopExpExplorer::new();
        let mut ex_edge = TopExpExplorer::new();
        // OCCT L312: BRep_Builder B;
        let mut b = BRepBuilder::new();
        // OCCT L313: B.MakeCompound(TopoDS::Compound(myOutLinedShape));
        self.my_out_lined_shape = b.make_compound(brep, Vec::new());
        // OCCT L314: NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>
        // ShapeMap;
        let mut shape_map: Vec<Shape> = Vec::new();

        // OCCT L316: for (exshell.Init(myOriginalShape, TopAbs_SHELL);
        // exshell.More(); exshell.Next()) — faces in a shell (open or close).
        ex_shell.init(
            brep,
            &self.my_original_shape,
            ShapeType::Shell,
            ShapeType::Shape,
        );
        while ex_shell.more() {
            // OCCT L318-319: TopoDS_Shell theShell; B.MakeShell(theShell);
            let the_shell = b.make_shell(brep);
            // OCCT L320: theShell.Closed(exshell.Current().Closed());
            let cur_shell = ex_shell.current().clone();
            let closed = cur_shell
                .as_shell()
                .map(|sd| sd.flags & tshape_flags::CLOSED != 0)
                .unwrap_or(false);
            set_flag_inplace(brep, &the_shell, tshape_flags::CLOSED, closed);

            // OCCT L322: for (exface.Init(exshell.Current(), TopAbs_FACE);
            // exface.More(); exface.Next())
            ex_face.init(brep, &cur_shell, ShapeType::Face, ShapeType::Shape);
            while ex_face.more() {
                let cur_face = ex_face.current().clone();
                // OCCT L324: if (ShapeMap.Add(exface.Current()))
                if shape_map_add(&mut shape_map, &cur_face) {
                    // OCCT L326: ProcessFace(TopoDS::Face(exface.Current()),
                    // theShell, MST);
                    self.process_face(brep, &cur_face, &the_shell, mst);
                }
                // exface.Next().
                ex_face.next(brep);
            }
            // OCCT L329: B.Add(myOutLinedShape, theShell);
            builder_add(brep, &self.my_out_lined_shape, &the_shell);
            // exshell.Next().
            ex_shell.next(brep);
        }

        // OCCT L332: for (exface.Init(myOriginalShape, TopAbs_FACE,
        // TopAbs_SHELL); exface.More(); exface.Next()) — faces not in a shell.
        ex_face.init(
            brep,
            &self.my_original_shape,
            ShapeType::Face,
            ShapeType::Shell,
        );
        while ex_face.more() {
            let cur_face = ex_face.current().clone();
            // OCCT L334: if (ShapeMap.Add(exface.Current()))
            if shape_map_add(&mut shape_map, &cur_face) {
                // OCCT L336: ProcessFace(TopoDS::Face(exface.Current()),
                // myOutLinedShape, MST);
                let s = self.my_out_lined_shape.clone();
                self.process_face(brep, &cur_face, &s, mst);
            }
            // exface.Next().
            ex_face.next(brep);
        }

        // OCCT L340: for (exedge.Init(myOriginalShape, TopAbs_EDGE,
        // TopAbs_FACE); exedge.More(); exedge.Next()) — edges not in a face.
        ex_edge.init(
            brep,
            &self.my_original_shape,
            ShapeType::Edge,
            ShapeType::Face,
        );
        while ex_edge.more() {
            let cur_edge = ex_edge.current().clone();
            // OCCT L342: B.Add(myOutLinedShape, exedge.Current());
            builder_add(brep, &self.my_out_lined_shape, &cur_edge);
            // exedge.Next().
            ex_edge.next(brep);
        }
    }
}

impl Default for OutLiner {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::geom::{Circle3, Curve3, Line3, Plane, Surface3};
    use std::sync::Arc;

    /// A shell of two unit squares (z = 0 and z = 1 planes) plus an extra
    /// circular arc edge on the z = 0 plane from (1,0,0) to (0,0,0) — the
    /// IntL candidate that matches the square edge (v10, v00) by vertices
    /// but not by curve type, so the SameEdge Extrema branch runs and the
    /// wire-append path follows.
    fn two_face_fixture() -> (BRep, Shape, Shape, Shape, Shape) {
        let mut brep = BRep::new();
        let mut b = BRepBuilder::new();

        // The z = 0 square: vertices (0,0), (1,0), (1,1), (0,1).
        let pts = [
            DVec3::new(0.0, 0.0, 0.0),
            DVec3::new(1.0, 0.0, 0.0),
            DVec3::new(1.0, 1.0, 0.0),
            DVec3::new(0.0, 1.0, 0.0),
        ];
        let vs: Vec<Shape> = pts.iter().map(|p| b.add_vertex(&mut brep, *p, 1e-7)).collect();
        let mut edges = Vec::new();
        for i in 0..4 {
            let a = pts[i];
            let bb = pts[(i + 1) % 4];
            let curve = Curve3::Line(Line3 {
                origin: a,
                direction: (bb - a).normalize(),
            });
            // OCCT BRep_Builder::Add(E, V) stores the start vertex FORWARD
            // and the end vertex REVERSED on the edge.
            let v_first = vs[i].clone();
            let mut v_last = vs[(i + 1) % 4].clone();
            v_last.orientation = Orientation::Reversed;
            edges.push(b.add_edge(&mut brep, Some(curve), v_first, v_last, [0.0, 1.0]));
        }
        let wire0 = brep.add_twire(edges.clone());
        let face1 = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::ZERO,
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire0,
            Vec::new(),
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );

        // The z = 1 square.
        let pts1: Vec<DVec3> = pts.iter().map(|p| *p + DVec3::new(0.0, 0.0, 1.0)).collect();
        let vs1: Vec<Shape> = pts1.iter().map(|p| b.add_vertex(&mut brep, *p, 1e-7)).collect();
        let mut edges1 = Vec::new();
        for i in 0..4 {
            let a = pts1[i];
            let bb = pts1[(i + 1) % 4];
            let curve = Curve3::Line(Line3 {
                origin: a,
                direction: (bb - a).normalize(),
            });
            let v_first = vs1[i].clone();
            let mut v_last = vs1[(i + 1) % 4].clone();
            v_last.orientation = Orientation::Reversed;
            edges1.push(b.add_edge(&mut brep, Some(curve), v_first, v_last, [0.0, 1.0]));
        }
        let wire1 = brep.add_twire(edges1);
        let face2 = brep.add_tface(
            Some(Surface3::Plane(Plane {
                origin: DVec3::new(0.0, 0.0, 1.0),
                normal: DVec3::new(0.0, 0.0, 1.0),
                u_dir: DVec3::new(1.0, 0.0, 0.0),
                v_dir: DVec3::new(0.0, 1.0, 0.0),
            })),
            wire1,
            Vec::new(),
            None,
            Some([0.0, 1.0, 0.0, 1.0]),
            Vec::new(),
            true,
        );

        // The circular arc from (1,0,0) (t = 0) to (0,0,0) (t = PI) on z = 0:
        // the semicircle of radius 0.5 around (0.5, 0, 0).
        let v10 = vs[1].clone();
        let mut v00 = vs[0].clone();
        v00.orientation = Orientation::Reversed;
        let arc = b.add_edge(
            &mut brep,
            Some(Curve3::Circle(Circle3::new_with_ref_dir(
                DVec3::new(0.5, 0.0, 0.0),
                DVec3::new(0.0, 0.0, 1.0),
                0.5,
                DVec3::new(1.0, 0.0, 0.0),
            ))),
            v10,
            v00,
            [0.0, std::f64::consts::PI],
        );

        // The original shape: the shell of the two faces.
        let shell = brep.add_tshell(vec![face1.clone(), face2.clone()]);

        (brep, shell, face1, face2, arc)
    }

    /// The MST map pre-bound for both faces (what DSFiller::Insert leaves
    /// behind for already-visited faces).
    fn bound_mst(brep: &BRep, faces: [&Shape; 2]) -> Vec<(Shape, BRepTopAdaptorTool)> {
        let arc_brep = Arc::new(brep.clone());
        faces
            .into_iter()
            .map(|f| {
                (
                    f.clone(),
                    BRepTopAdaptorTool::new_face(arc_brep.clone(), f, 1.0e-7),
                )
            })
            .collect()
    }

    /// OCCT anchor: the three constructors (cxx L44-59), the lxx accessors
    /// (L19-50) round-trip and DataStructure() reachability.
    #[test]
    fn out_liner_ctors_and_accessors() {
        let (brep, shell, face1, face2, arc) = two_face_fixture();

        // HLRTopoBRep_OutLiner() — all members null / default.
        let mut ol = OutLiner::new();
        assert!(ol.original_shape().is_null());
        assert!(ol.out_lined_shape().is_null());
        // DataStructure() is reachable and answers through the Data maps.
        assert!(!ol.data_structure().edge_has_spl_e(&arc));

        // HLRTopoBRep_OutLiner(OriS).
        let mut ol = OutLiner::with_original_shape(&shell);
        assert!(ol.original_shape().is_same(&shell));
        assert!(ol.out_lined_shape().is_null());

        // OutLinedShape(OutS) setter + getter round-trip.
        ol.set_out_lined_shape(&face2);
        assert!(ol.out_lined_shape().is_same(&face2));

        // HLRTopoBRep_OutLiner(OriS, OutS).
        let mut ol = OutLiner::with_original_and_outlined(&shell, &face1);
        assert!(ol.original_shape().is_same(&shell));
        assert!(ol.out_lined_shape().is_same(&face1));
        // OriginalShape(OriS) setter round-trip.
        ol.set_original_shape(&face2);
        assert!(ol.original_shape().is_same(&face2));

        let _ = brep;
        let _ = face1;
    }

    /// OCCT anchor: BuildShape (cxx L306-344) + ProcessFace (cxx L94-304)
    /// over the two-face shell with one IntL arc on face1 — the result is a
    /// compound of one shell of two rebuilt faces, the IntL arc lands as an
    /// INTERNAL edge of an extra wire on the face1 copy, the NF -> F and
    /// E -> F old-shape records are set, and MST gains the NF bindings.
    #[test]
    fn out_liner_process_face_build_shape() {
        let (mut brep, shell, face1, face2, arc) = two_face_fixture();
        let mut mst = bound_mst(&brep, [&face1, &face2]);

        let mut ol = OutLiner::with_original_shape(&shell);
        // The DSFiller-produced state: face1 has one internal outline (the
        // arc); the arc is not split.
        ol.data_structure().add_int_l(&face1).push(arc.clone());

        ol.build_shape(&mut brep, &mut mst);

        // The result is a compound (MakeCompound, cxx L313).
        let result = ol.out_lined_shape().clone();
        assert_eq!(result.shape_type(), ShapeType::Compound);

        // One shell in the compound (B.Add(myOutLinedShape, theShell), L329).
        let children = raw_subshapes(&brep, &result);
        assert_eq!(children.len(), 1);
        let res_shell = &children[0];
        assert_eq!(res_shell.shape_type(), ShapeType::Shell);

        // Two faces in the shell, both fresh copies (EmptyCopied).
        let faces = raw_subshapes(&brep, res_shell);
        assert_eq!(faces.len(), 2);
        assert!(faces.iter().all(|f| !f.is_same(&face1) && !f.is_same(&face2)));

        // Each face copy maps back to its original (AddOldS(NF, F), L300).
        let ds = ol.data_structure();
        let nf1 = faces
            .iter()
            .find(|f| ds.new_s_old_s(f).is_same(&face1))
            .expect("face1 copy");
        let _nf2 = faces
            .iter()
            .find(|f| ds.new_s_old_s(f).is_same(&face2))
            .expect("face2 copy");

        // face1 copy: outer wire (4 edges) + one internal wire carrying the
        // arc as an INTERNAL edge (B.Add(NF, W), L256).
        let nf1_children = raw_subshapes(&brep, nf1);
        assert_eq!(nf1_children.len(), 2);
        let inner_wire = &nf1_children[1];
        assert_eq!(inner_wire.shape_type(), ShapeType::Wire);
        let wire_edges = raw_subshapes(&brep, inner_wire);
        assert_eq!(wire_edges.len(), 1);
        assert!(wire_edges[0].is_same(&arc));
        assert_eq!(wire_edges[0].orientation, Orientation::Internal);
        // AddOldS(E, F) recorded the arc against face1 (L250).
        assert!(ds.new_s_old_s(&arc).is_same(&face1));

        // face2 copy: outer wire only.
        let nf2 = faces.iter().find(|f| ds.new_s_old_s(f).is_same(&face2)).unwrap();
        assert_eq!(raw_subshapes(&brep, nf2).len(), 1);

        // MST.Bind(NF, MST.ChangeFind(F)) (L301): two new bindings whose
        // tools carry the loaded surface.
        assert_eq!(mst.len(), 4);
        for f in &faces {
            let pos = mst
                .iter()
                .position(|(k, _)| k.ptr_id() == f.ptr_id())
                .expect("NF bound in MST");
            assert!(mst[pos].1.get_surface().is_some());
        }
    }

    /// OCCT anchor: the Fill control flow (cxx L63-92) — a null original is
    /// a no-op, an already built outLined shape is a no-op (L70), and the
    /// main path runs Insert + BuildShape producing the compound.
    #[test]
    fn out_liner_fill_control_flow() {
        let (mut brep, shell, face1, face2, _arc) = two_face_fixture();
        let mut mst = bound_mst(&brep, [&face1, &face2]);
        let projector = Projector::new();

        // Null original: the guard at cxx L68 skips everything.
        let mut empty = OutLiner::new();
        empty.fill(&mut brep, &projector, &mut mst, 0);
        assert!(empty.out_lined_shape().is_null());

        // Main path: Insert + BuildShape.
        let mut ol = OutLiner::with_original_shape(&shell);
        ol.fill(&mut brep, &projector, &mut mst, 0);
        let result = ol.out_lined_shape().clone();
        assert!(!result.is_null());
        assert_eq!(result.shape_type(), ShapeType::Compound);
        let children = raw_subshapes(&brep, &result);
        assert_eq!(children.len(), 1);
        assert_eq!(raw_subshapes(&brep, &children[0]).len(), 2);

        // Second Fill: myOutLinedShape is no longer null (cxx L70) — the
        // result is untouched and MST gains no entries.
        let mst_len = mst.len();
        let result_id = result.ptr_id();
        ol.fill(&mut brep, &projector, &mut mst, 0);
        assert_eq!(ol.out_lined_shape().ptr_id(), result_id);
        assert_eq!(mst.len(), mst_len);
    }
}
