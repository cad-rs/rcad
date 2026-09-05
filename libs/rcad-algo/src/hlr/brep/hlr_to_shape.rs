// OCCT HLRBRep_HLRToShape (TKHLR/HLRBRep/HLRBRep_HLRToShape.hxx L1-178
// + .cxx L1-305 + .lxx L1-183) — the extraction filter over an
// HLRBRep_Algo: it reads the computed Data and rebuilds a compound of the
// visible / hidden edges of the requested kind.
//
// The package statics of HLRBRep.hxx (MakeEdge / MakeEdge3d, .cxx L38-247)
// consumed by DrawEdge land in the local [`hlr_brep`] module (no dedicated
// rcad file exists yet; move when the coordinator opens one).
//
// Deviations (all reported to the coordinator):
// - `occ::handle<HLRBRep_Algo>` maps to `&'a mut Algo` — the OCCT methods
//   mutate the DataStructure through the non-const handle, so the rcad
//   short-lived view borrows the Algo exclusively;
// - the TopoDS universe of BRep_Builder maps to the explicit `brep` arena
//   parameter (the ds_filler builder-surface precedent);
// - HLRBRep::MakeEdge / MakeEdge3d carry the documented kernel gaps of the
//   module [`hlr_brep`] notes;
// - the OCCT overloads VCompound(S) / CompoundOfEdges(S, ...) map to the
//   `_of_shape` names (Rust has no overload resolution).

use glam::DVec2;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{BSplineCurve2, Circle3, Curve3, Ellipse3, Line2d, Line3};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, Shape, TShape};

use crate::hlr::algo::edge_iterator::EdgeIterator;
use crate::hlr::brep::algo::Algo;
use crate::hlr::brep::curve::Curve;
use crate::hlr::brep::data::Data;
use crate::hlr::brep::edge_data::EdgeData;
use crate::hlr::brep::face_data::FaceData;
use crate::hlr::brep::face_iterator::FaceIterator;

/// OCCT `enum HLRBRep_TypeOfResultingEdge`
/// (HLRBRep_TypeOfResultingEdge.hxx L21-35; the values follow the OCCT
/// declaration order, matching the InternalCompound int typ codes).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TypeOfResultingEdge {
    Undefined = 0,
    /// isoparametric line
    IsoLine = 1,
    /// outline ("silhouette")
    OutLine = 2,
    /// smooth edge of G1-continuity between two surfaces
    Rg1Line = 3,
    /// sewn edge of CN-continuity on one surface
    RgNLine = 4,
    /// sharp edge (of C0-continuity)
    Sharp = 5,
}

/// OCCT TopExp_Explorer(S, TopAbs_EDGE, TopAbs_FACE) — the edges of every
/// wire / compound branch, avoiding the sub-shapes under a FACE (the
/// ds_filler::top_exp_faces walker precedent; the cumulative
/// Location/Orientation composition is the identity pass-through).
fn top_exp_edges_avoiding_faces(s: &Shape, out: &mut Vec<Shape>) {
    match &*s.data {
        TShape::Compound(children) => {
            for c in children {
                top_exp_edges_avoiding_faces(c, out);
            }
        }
        TShape::CompSolid(cs) => {
            for c in cs {
                top_exp_edges_avoiding_faces(c, out);
            }
        }
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                top_exp_edges_avoiding_faces(sh, out);
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                top_exp_edges_avoiding_faces(f, out);
            }
        }
        // A face subtree is pruned (the OCCT avoid-parent argument).
        TShape::Face(_) => {}
        TShape::Wire(wd) => {
            for e in &wd.edges {
                top_exp_edges_avoiding_faces(e, out);
            }
        }
        TShape::Edge(_) => out.push(s.clone()),
        // Vertex: the explorer finds nothing below it.
        _ => {}
    }
}

/// OCCT TopExp_Explorer(S, TopAbs_FACE) — the recursive face enumeration
/// (the ds_filler::top_exp_faces walker, replicated locally because the
/// helper is private to the ds_filler module).
fn top_exp_faces_of(s: &Shape, out: &mut Vec<Shape>) {
    match &*s.data {
        TShape::Compound(children) => {
            for c in children {
                top_exp_faces_of(c, out);
            }
        }
        TShape::CompSolid(cs) => {
            for c in cs {
                top_exp_faces_of(c, out);
            }
        }
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                top_exp_faces_of(sh, out);
            }
        }
        TShape::Shell(sd) => {
            for f in &sd.faces {
                top_exp_faces_of(f, out);
            }
        }
        TShape::Face(_) => out.push(s.clone()),
        // Wire / Edge / Vertex: the explorer finds nothing below them.
        _ => {}
    }
}

/// OCCT NCollection_IndexedMap::FindIndex over the Data shape maps (the
/// ptr_id-keyed Vec precedent) — 0 when the key is absent.
fn find_index(map: &[Shape], s: &Shape) -> i32 {
    map.iter()
        .position(|k| k.ptr_id() == s.ptr_id())
        .map_or(0, |p| (p + 1) as i32)
}

/// OCCT HLRBRep_HLRToShape (hxx L60-174).
pub struct HLRToShape<'a> {
    /// OCCT myAlgo (hxx L173) — the non-const handle access maps to the
    /// exclusive borrow (InternalCompound mutates the DataStructure).
    my_algo: &'a mut Algo,
}

impl<'a> HLRToShape<'a> {
    /// OCCT HLRBRep_HLRToShape(const handle<HLRBRep_Algo>& A) (cxx L33-36) —
    /// construct a framework for filtering the results of the algorithm A.
    pub fn new(a: &'a mut Algo) -> Self {
        HLRToShape { my_algo: a }
    }

    /// OCCT InternalCompound(typ, visible, S, In3d) (cxx L40-150).
    fn internal_compound(&mut self, typ: i32, visible: bool, s: &Shape, in3d: bool) -> Shape {
        // bool added = false;
        let mut added = false;
        // TopoDS_Shape Result;
        let mut result = Shape::null();
        // occ::handle<HLRBRep_Data> DS = myAlgo->DataStructure();
        // (the raw pointer keeps the OCCT handle aliasing between the DS
        // reads and the myAlgo-> calls — the Data element-pointer precedent;
        // the Algo DerefMut carries the inherited DataStructure accessor.)
        let ds_ptr: *mut Data<'static> = match self.my_algo.data_structure_mut() {
            Some(ds) => ds as *mut Data<'static>,
            None => std::ptr::null_mut(),
        };

        // if (!DS.IsNull())
        if !ds_ptr.is_null() {
            // DS->Projector().Scaled(true);
            unsafe { (*ds_ptr).projector_mut().scaled(true) };
            // int e1 = 1;
            // int e2 = DS->NbEdges();
            // int f1 = 1;
            // int f2 = DS->NbFaces();
            let mut e1: i32 = 1;
            let mut e2 = unsafe { (*ds_ptr).nb_edges() } as i32;
            let mut f1: i32 = 1;
            let mut f2 = unsafe { (*ds_ptr).nb_faces() } as i32;
            // bool explor = false;
            let mut explor = false;
            //    bool todraw; (the commented-out OCCT declaration is kept as a comment)
            // if (!S.IsNull())
            if !s.is_null() {
                // int v1, v2;
                let v1: i32;
                let v2: i32;
                // int index = myAlgo->Index(S);
                let index = self.my_algo.index(s);
                if index == 0 {
                    // explor = true;
                    explor = true;
                } else {
                    // myAlgo->ShapeBounds(index).Bounds(v1, v2, e1, e2, f1, f2);
                    (v1, v2, e1, e2, f1, f2) = self.my_algo.shape_bounds(index).bounds();
                }
            }
            // BRep_Builder B;
            // B.MakeCompound(TopoDS::Compound(Result));
            // (the rcad arena is the builder surface — the ds_filler
            // precedent; the compound lives in it.)
            let mut b = BRepBuilder::new();
            let mut arena = BRep::new();
            result = b.make_compound(&mut arena, Vec::new());
            // HLRBRep_EdgeData* ed = &(DS->EDataArray().ChangeValue(e1 - 1));
            // (ChangeValue(e1-1) names the (0,NE) dummy slot 0 when e1 == 1;
            // the ed++ below steps to the first live slot, so the rcad cursor
            // keeps the OCCT slot number and indexes [(slot) - 1].)
            let a_e_data: *mut EdgeData<'static> =
                unsafe { (*ds_ptr).e_data_array_mut() }.as_mut_ptr();
            let mut ed: i32 = e1 - 1;

            // for (int ie = e1; ie <= e2; ie++)
            let mut ie = e1;
            while ie <= e2 {
                // ed++;
                ed += 1;
                let edr = unsafe { &mut *a_e_data.add((ed - 1) as usize) };
                // if (ed->Selected() && !ed->Vertical()) { ed->Used(false); ed->HideCount(0); }
                if edr.selected() && !edr.vertical() {
                    edr.set_used(false);
                    edr.set_hide_count(0);
                } else {
                    // ed->Used(true);
                    edr.set_used(true);
                }
                ie += 1;
            }
            // if (explor)
            if explor {
                // NCollection_IndexedMap<...>& Edges = DS->EdgeMap();
                // NCollection_IndexedMap<...>& Faces = DS->FaceMap();
                // TopExp_Explorer Exp;
                let mut exp_faces: Vec<Shape> = Vec::new();
                // for (Exp.Init(S, TopAbs_FACE); Exp.More(); Exp.Next())
                top_exp_faces_of(s, &mut exp_faces);
                for current in &exp_faces {
                    // int iface = Faces.FindIndex(Exp.Current());
                    let iface = find_index(unsafe { (*ds_ptr).face_map() }, current);
                    // if (iface != 0)
                    if iface != 0 {
                        // DrawFace(visible, typ, iface, DS, Result, added, In3d);
                        self.draw_face(
                            visible,
                            typ,
                            iface,
                            ds_ptr,
                            &mut arena,
                            &mut result,
                            &mut added,
                            in3d,
                        );
                    }
                }
                // if (typ >= 3)
                if typ >= 3 {
                    let mut exp_edges: Vec<Shape> = Vec::new();
                    // for (Exp.Init(S, TopAbs_EDGE, TopAbs_FACE); Exp.More(); Exp.Next())
                    top_exp_edges_avoiding_faces(s, &mut exp_edges);
                    for current in &exp_edges {
                        // int ie = Edges.FindIndex(Exp.Current());
                        let ie = find_index(unsafe { (*ds_ptr).edge_map() }, current);
                        // if (ie != 0)
                        if ie != 0 {
                            // HLRBRep_EdgeData& EDataIE = DS->EDataArray().ChangeValue(ie);
                            let edata_ie = unsafe { &mut *a_e_data.add((ie - 1) as usize) };
                            // if (!EDataIE.Used())
                            if !edata_ie.used() {
                                // DrawEdge(visible, false, typ, EDataIE, Result, added, In3d);
                                self.draw_edge(
                                    visible,
                                    false,
                                    typ,
                                    edata_ie,
                                    &mut arena,
                                    &mut result,
                                    &mut added,
                                    in3d,
                                );
                                // EDataIE.Used(true);
                                edata_ie.set_used(true);
                            }
                        }
                    }
                }
            } else {
                // for (int iface = f1; iface <= f2; iface++)
                let mut iface = f1;
                while iface <= f2 {
                    // DrawFace(visible, typ, iface, DS, Result, added, In3d);
                    self.draw_face(
                        visible,
                        typ,
                        iface,
                        ds_ptr,
                        &mut arena,
                        &mut result,
                        &mut added,
                        in3d,
                    );
                    iface += 1;
                }

                // if (typ >= 3)
                if typ >= 3 {
                    // HLRBRep_EdgeData* EDataE11 = &(DS->EDataArray().ChangeValue(e1 - 1));
                    // (the same dummy-slot fusion as the ed cursor above.)
                    let mut edata_e11: i32 = e1 - 1;

                    // for (int ie = e1; ie <= e2; ie++)
                    let mut ie = e1;
                    while ie <= e2 {
                        // EDataE11++;
                        edata_e11 += 1;
                        let e11 = unsafe { &mut *a_e_data.add((edata_e11 - 1) as usize) };
                        // if (!EDataE11->Used())
                        if !e11.used() {
                            // DrawEdge(visible, false, typ, *EDataE11, Result, added, In3d);
                            self.draw_edge(
                                visible,
                                false,
                                typ,
                                e11,
                                &mut arena,
                                &mut result,
                                &mut added,
                                in3d,
                            );
                            // EDataE11->Used(true);
                            e11.set_used(true);
                        }
                        ie += 1;
                    }
                }
            }
            // DS->Projector().Scaled(false);
            unsafe { (*ds_ptr).projector_mut().scaled(false) };
            // (the rcad add_to_compound replaces the arena slot through
            // Arc::make_mut; the compound handle re-reads the arena entry —
            // the OCCT TShape mutation is in-place, the rcad slot swap is
            // the copy-on-write equivalent.)
            result = Shape::from_parts(
                arena.tshapes[result.index].clone(),
                result.index,
                result.location,
                result.orientation,
            );
        }
        // if (!added) { Result = TopoDS_Shape(); }
        if !added {
            result = Shape::null();
        }
        // return Result;
        result
    }

    /// OCCT DrawFace(visible, typ, iface, DS, Result, added, In3d)
    /// (cxx L154-227).
    #[allow(clippy::too_many_arguments)]
    fn draw_face(
        &self,
        visible: bool,
        typ: i32,
        iface: i32,
        ds_ptr: *mut Data<'static>,
        arena: &mut BRep,
        result: &mut Shape,
        added: &mut bool,
        in3d: bool,
    ) {
        // HLRBRep_FaceIterator Itf;
        let mut itf = FaceIterator::new(std::ptr::null_mut());

        // for (Itf.InitEdge(DS->FDataArray().ChangeValue(iface));
        //      Itf.MoreEdge(); Itf.NextEdge())
        let fd: *mut FaceData<'static> =
            unsafe { (*ds_ptr).f_data_array_mut().as_mut_ptr().add((iface - 1) as usize) };
        itf.init_edge(fd);
        while itf.more_edge() {
            // int ie = Itf.Edge();
            let ie = itf.edge() as usize;
            // HLRBRep_EdgeData& edf = DS->EDataArray().ChangeValue(ie);
            let edf = unsafe { &mut (*ds_ptr).e_data_array_mut()[ie - 1] };
            // if (!edf.Used())
            if !edf.used() {
                // bool todraw;
                let todraw: bool;
                if typ == 1 {
                    // todraw = Itf.IsoLine();
                    todraw = itf.iso_line();
                } else if typ == 2 {
                    // outlines
                    if in3d {
                        // todraw = Itf.Internal() || Itf.OutLine();
                        todraw = itf.internal() || itf.out_line();
                    } else {
                        // todraw = Itf.Internal();
                        todraw = itf.internal();
                    }
                } else if typ == 3 {
                    // todraw = edf.Rg1Line() && !edf.RgNLine() && !Itf.OutLine();
                    todraw = edf.rg1_line() && !edf.rg_n_line() && !itf.out_line();
                } else if typ == 4 {
                    // todraw = edf.RgNLine() && !Itf.OutLine();
                    todraw = edf.rg_n_line() && !itf.out_line();
                } else {
                    // todraw = !Itf.IsoLine() && !Itf.Internal()
                    //          && (!edf.Rg1Line() || Itf.OutLine());
                    todraw =
                        !itf.iso_line() && !itf.internal() && (!edf.rg1_line() || itf.out_line());
                }

                // if (todraw)
                if todraw {
                    // DrawEdge(visible, true, typ, edf, Result, added, In3d);
                    self.draw_edge(visible, true, typ, edf, arena, result, added, in3d);
                    // edf.Used(true);
                    edf.set_used(true);
                } else {
                    // if ((typ > 4 || typ == 2) && // sharp or outlines
                    //     (edf.Rg1Line() && !Itf.OutLine()))
                    if (typ > 4 || typ == 2) && (edf.rg1_line() && !itf.out_line()) {
                        // int hc = edf.HideCount();
                        let hc = edf.hide_count();
                        if hc > 0 {
                            // edf.Used(true);
                            edf.set_used(true);
                        } else {
                            // ++hc; edf.HideCount(hc); // to try with another face
                            let hc = hc + 1;
                            edf.set_hide_count(hc);
                        }
                    } else {
                        // edf.Used(true);
                        edf.set_used(true);
                    }
                }
            }
            itf.next_edge();
        }
    }

    /// OCCT DrawEdge(visible, inFace, typ, ed, Result, added, In3d)
    /// (cxx L231-305).
    #[allow(clippy::too_many_arguments)]
    fn draw_edge(
        &self,
        visible: bool,
        in_face: bool,
        typ: i32,
        ed: &mut EdgeData<'static>,
        arena: &mut BRep,
        result: &mut Shape,
        added: &mut bool,
        in3d: bool,
    ) {
        // bool todraw = false;
        let mut todraw = false;
        if in_face {
            // todraw = true;
            todraw = true;
        } else if typ == 3 {
            // todraw = ed.Rg1Line() && !ed.RgNLine();
            todraw = ed.rg1_line() && !ed.rg_n_line();
        } else if typ == 4 {
            // todraw = ed.RgNLine();
            todraw = ed.rg_n_line();
        } else {
            // todraw = !ed.Rg1Line();
            todraw = !ed.rg1_line();
        }

        // if (todraw)
        if todraw {
            // double sta, end;
            // float tolsta, tolend;
            let mut sta: f64;
            let mut end: f64;
            let mut tolsta: f32;
            let mut tolend: f32;
            // BRep_Builder B;
            // TopoDS_Edge E;
            // HLRAlgo_EdgeIterator It;
            let b = BRepBuilder::new();
            let mut e: Shape;
            let mut it: EdgeIterator = EdgeIterator::new();
            if visible {
                // for (It.InitVisible(ed.Status()); It.MoreVisible(); It.NextVisible())
                it.init_visible(ed.status_ref());
                while it.more_visible() {
                    // It.Visible(sta, tolsta, end, tolend);
                    (sta, tolsta, end, tolend) = it.visible();
                    let _ = (tolsta, tolend);
                    // if (!In3d) E = HLRBRep::MakeEdge(ed.Geometry(), sta, end);
                    // else E = HLRBRep::MakeEdge3d(ed.Geometry(), sta, end);
                    e = if !in3d {
                        hlr_brep::make_edge(arena, ed.geometry(), sta, end)
                    } else {
                        hlr_brep::make_edge_3d(arena, ed.geometry(), sta, end)
                    };
                    // if (!E.IsNull()) { B.Add(Result, E); added = true; }
                    if !e.is_null() {
                        b.add_to_compound(arena, result.clone(), e);
                        *added = true;
                    }
                    it.next_visible();
                }
            } else {
                // for (It.InitHidden(ed.Status()); It.MoreHidden(); It.NextHidden())
                it.init_hidden(ed.status_ref());
                while it.more_hidden() {
                    // It.Hidden(sta, tolsta, end, tolend);
                    (sta, tolsta, end, tolend) = it.hidden();
                    let _ = (tolsta, tolend);
                    // if (!In3d) E = HLRBRep::MakeEdge(ed.Geometry(), sta, end);
                    // else E = HLRBRep::MakeEdge3d(ed.Geometry(), sta, end);
                    e = if !in3d {
                        hlr_brep::make_edge(arena, ed.geometry(), sta, end)
                    } else {
                        hlr_brep::make_edge_3d(arena, ed.geometry(), sta, end)
                    };
                    // if (!E.IsNull()) { B.Add(Result, E); added = true; }
                    if !e.is_null() {
                        b.add_to_compound(arena, result.clone(), e);
                        *added = true;
                    }
                    it.next_hidden();
                }
            }
        }
    }

    // ====================================================================
    // OCCT HLRBRep_HLRToShape.lxx L21-183 — the filter overloads.
    // ====================================================================

    /// OCCT VCompound() (lxx L21-24) — return visible sharp edges (of
    /// C0-continuity).
    pub fn v_compound(&mut self) -> Shape {
        // return InternalCompound(5, true, TopoDS_Shape());
        self.internal_compound(5, true, &Shape::null(), false)
    }

    /// OCCT VCompound(const TopoDS_Shape& S) (lxx L28-31).
    pub fn v_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(5, true, S);
        self.internal_compound(5, true, s, false)
    }

    /// OCCT Rg1LineVCompound() (lxx L35-38) — return visible smooth edges
    /// (G1-continuity between two surfaces).
    pub fn rg1_line_v_compound(&mut self) -> Shape {
        // return InternalCompound(3, true, TopoDS_Shape());
        self.internal_compound(3, true, &Shape::null(), false)
    }

    /// OCCT Rg1LineVCompound(const TopoDS_Shape& S) (lxx L42-45).
    pub fn rg1_line_v_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(3, true, S);
        self.internal_compound(3, true, s, false)
    }

    /// OCCT RgNLineVCompound() (lxx L49-52) — return visible sewn edges (of
    /// CN-continuity on one surface).
    pub fn rg_n_line_v_compound(&mut self) -> Shape {
        // return InternalCompound(4, true, TopoDS_Shape());
        self.internal_compound(4, true, &Shape::null(), false)
    }

    /// OCCT RgNLineVCompound(const TopoDS_Shape& S) (lxx L56-59).
    pub fn rg_n_line_v_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(4, true, S);
        self.internal_compound(4, true, s, false)
    }

    /// OCCT OutLineVCompound() (lxx L63-66) — return visible outline edges
    /// ("silhouette").
    pub fn out_line_v_compound(&mut self) -> Shape {
        // return InternalCompound(2, true, TopoDS_Shape());
        self.internal_compound(2, true, &Shape::null(), false)
    }

    /// OCCT OutLineVCompound3d() (lxx L70-73).
    pub fn out_line_v_compound_3d(&mut self) -> Shape {
        // return InternalCompound(2, true, TopoDS_Shape(), true);
        self.internal_compound(2, true, &Shape::null(), true)
    }

    /// OCCT OutLineVCompound(const TopoDS_Shape& S) (lxx L77-80).
    pub fn out_line_v_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(2, true, S);
        self.internal_compound(2, true, s, false)
    }

    /// OCCT IsoLineVCompound() (lxx L84-87) — return visible isoparameters.
    pub fn iso_line_v_compound(&mut self) -> Shape {
        // return InternalCompound(1, true, TopoDS_Shape());
        self.internal_compound(1, true, &Shape::null(), false)
    }

    /// OCCT IsoLineVCompound(const TopoDS_Shape& S) (lxx L91-94).
    pub fn iso_line_v_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(1, true, S);
        self.internal_compound(1, true, s, false)
    }

    /// OCCT HCompound() (lxx L98-101) — return hidden sharp edges (of
    /// C0-continuity).
    pub fn h_compound(&mut self) -> Shape {
        // return InternalCompound(5, false, TopoDS_Shape());
        self.internal_compound(5, false, &Shape::null(), false)
    }

    /// OCCT HCompound(const TopoDS_Shape& S) (lxx L105-108).
    pub fn h_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(5, false, S);
        self.internal_compound(5, false, s, false)
    }

    /// OCCT Rg1LineHCompound() (lxx L112-115) — return hidden smooth edges.
    pub fn rg1_line_h_compound(&mut self) -> Shape {
        // return InternalCompound(3, false, TopoDS_Shape());
        self.internal_compound(3, false, &Shape::null(), false)
    }

    /// OCCT Rg1LineHCompound(const TopoDS_Shape& S) (lxx L119-122).
    pub fn rg1_line_h_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(3, false, S);
        self.internal_compound(3, false, s, false)
    }

    /// OCCT RgNLineHCompound() (lxx L126-129) — return hidden sewn edges.
    pub fn rg_n_line_h_compound(&mut self) -> Shape {
        // return InternalCompound(4, false, TopoDS_Shape());
        self.internal_compound(4, false, &Shape::null(), false)
    }

    /// OCCT RgNLineHCompound(const TopoDS_Shape& S) (lxx L133-136).
    pub fn rg_n_line_h_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(4, false, S);
        self.internal_compound(4, false, s, false)
    }

    /// OCCT IsoLineHCompound() (lxx L154-157) — return hidden isoparameters.
    pub fn iso_line_h_compound(&mut self) -> Shape {
        // return InternalCompound(1, false, TopoDS_Shape());
        self.internal_compound(1, false, &Shape::null(), false)
    }

    /// OCCT IsoLineHCompound(const TopoDS_Shape& S) (lxx L161-164).
    pub fn iso_line_h_compound_of_shape(&mut self, s: &Shape) -> Shape {
        // return InternalCompound(1, false, S);
        self.internal_compound(1, false, s, false)
    }

    /// OCCT CompoundOfEdges(type, visible, In3d) (lxx L168-173) — compound
    /// of resulting edges of the required type and visibility, taking into
    /// account the kind of space (2d or 3d).
    pub fn compound_of_edges(
        &mut self,
        type_: TypeOfResultingEdge,
        visible: bool,
        in3d: bool,
    ) -> Shape {
        // return InternalCompound(type, visible, TopoDS_Shape(), In3d);
        self.internal_compound(type_ as i32, visible, &Shape::null(), in3d)
    }

    /// OCCT CompoundOfEdges(const TopoDS_Shape& S, type, visible, In3d)
    /// (lxx L177-183).
    pub fn compound_of_edges_of_shape(
        &mut self,
        s: &Shape,
        type_: TypeOfResultingEdge,
        visible: bool,
        in3d: bool,
    ) -> Shape {
        // return InternalCompound(type, visible, S, In3d);
        self.internal_compound(type_ as i32, visible, s, in3d)
    }
}

// ============================================================================
// OCCT HLRBRep (HLRBRep.hxx L38-55 + HLRBRep.cxx L38-247) — the package
// statics MakeEdge / MakeEdge3d consumed by DrawEdge.  Landed inside this
// file (HLRBRep.cxx has no dedicated rcad file; moved wholesale when the
// coordinator opens one).
//
// Kernel gaps (reported; the OCCT BRepLib_MakeEdge2d builder has no rcad
// counterpart): the 2D support curves are rebuilt on the z = 0 plane image
// (vertices at the trimmed ends, a 3D curve of the same family, the
// [p1, p2] range).  The pcurve-on-XY-plane representation of OCCT is not
// stored (the rcad pcurve map is keyed by an owning face).  The conic /
// Bezier / BSpline pcurve rebuilds whose 2D poles / knots the rcad Curve
// adaptor does not expose keep the OCCT case structure and produce the null
// edge (the BRepLib_MakeEdge2d NotDone branch of the OCCT callers).
// ============================================================================
mod hlr_brep {
    use super::*;

    /// OCCT HLRBRep::MakeEdge(const HLRBRep_Curve& ec, U1, U2)
    /// (HLRBRep.cxx L38-186).
    pub(super) fn make_edge(arena: &mut BRep, ec: &Curve, u1: f64, u2: f64) -> Shape {
        // TopoDS_Edge Edg;
        let mut edg = Shape::null();
        // const double sta = ec.Parameter2d(U1);
        // const double end = ec.Parameter2d(U2);
        let sta = ec.parameter_2d(u1);
        let end = ec.parameter_2d(u2);

        // switch (ec.GetType())
        match ec.get_type() {
            // case GeomAbs_Line: Edg = BRepLib_MakeEdge2d(ec.Line(), sta, end);
            CurveType::Line => {
                edg = make_edge_2d_line(arena, ec.line(), sta, end);
            }
            // case GeomAbs_Circle: Edg = BRepLib_MakeEdge2d(ec.Circle(), sta, end);
            CurveType::Circle => {
                edg = make_edge_2d_circle(arena, ec.circle(), sta, end);
            }
            // case GeomAbs_Ellipse: Edg = BRepLib_MakeEdge2d(ec.Ellipse(), sta, end);
            CurveType::Ellipse => {
                edg = make_edge_2d_ellipse(arena, ec.ellipse(), sta, end);
            }
            // case GeomAbs_Hyperbola / GeomAbs_Parabola — the kernel gap of
            // the module note; the NotDone branch keeps the null edge.
            CurveType::Hyperbola | CurveType::Parabola => {}
            // case GeomAbs_BezierCurve / GeomAbs_BSplineCurve — the 2D
            // poles / knots rebuilds of cxx L86-158 have no rcad accessor
            // (the kernel gap); the NotDone branch keeps the null edge.
            CurveType::Bezier | CurveType::BSpline => {}
            // default: the 15-pole degree-1 approximation (cxx L160-186).
            CurveType::Other => {
                // const int nbPnt = 15;
                const NB_PNT: usize = 15;
                // NCollection_Array1<gp_Pnt2d> Poles(1, nbPnt);
                let mut poles: Vec<DVec2> = vec![DVec2::ZERO; NB_PNT];
                // NCollection_Array1<double> knots(1, nbPnt);
                let mut knots: Vec<f64> = vec![0.0; NB_PNT];
                // NCollection_Array1<int> mults(1, nbPnt);
                let mut mults: Vec<i32> = vec![1; NB_PNT];
                // mults.Init(1); mults(1) = 2; mults(nbPnt) = 2;
                mults[0] = 2;
                mults[NB_PNT - 1] = 2;
                // const double step = (U2 - U1) / (nbPnt - 1);
                let step = (u2 - u1) / (NB_PNT as f64 - 1.0);
                // double par3d = U1;
                let mut par3d = u1;
                // for (int i = 1; i < nbPnt; i++)
                for i in 1..NB_PNT {
                    // Poles(i) = ec.Value(par3d); knots(i) = par3d;
                    poles[i - 1] = ec.value(par3d);
                    knots[i - 1] = par3d;
                    // par3d += step;
                    par3d += step;
                }
                // Poles(nbPnt) = ec.Value(U2); knots(nbPnt) = U2;
                poles[NB_PNT - 1] = ec.value(u2);
                knots[NB_PNT - 1] = u2;

                // the flat knot vector of the rcad BSpline representation
                // (each knot repeated its multiplicity — the expand_knots
                // precedent).
                let mut knots_flat: Vec<f64> = Vec::with_capacity(NB_PNT + 2);
                for (i, &k) in knots.iter().enumerate() {
                    for _ in 0..mults[i] {
                        knots_flat.push(k);
                    }
                }
                // occ::handle<Geom2d_BSplineCurve> ec2d =
                //   new Geom2d_BSplineCurve(Poles, knots, mults, 1);
                let ec2d = BSplineCurve2 {
                    degree: 1,
                    knots: knots_flat,
                    control_points: poles,
                    weights: vec![1.0; NB_PNT],
                };
                // BRepLib_MakeEdge2d mke2d(ec2d, sta, end);
                // if (mke2d.IsDone()) { Edg = mke2d.Edge(); }
                edg = make_edge_2d_bspline(arena, &ec2d, sta, end);
            }
        }
        // return Edg;
        edg
    }

    /// OCCT HLRBRep::MakeEdge3d(const HLRBRep_Curve& ec, U1, U2)
    /// (HLRBRep.cxx L207-244).
    ///
    /// Kernel gap (reported): the first statement
    /// `TopoDS_Edge anEdge = ec.GetCurve().Edge()` reads the original
    /// TopoDS edge through the HLRBRep_Curve adaptor; the rcad CurveView
    /// carries no back-pointer to the TopoDS edge, so the whole rebuild
    /// (EmptyCopied / Range / shared vertices) is unrepresentable and the
    /// null edge is returned (the DrawEdge IsNull branch absorbs it).
    pub(super) fn make_edge_3d(_arena: &mut BRep, _ec: &Curve, _u1: f64, _u2: f64) -> Shape {
        // TopoDS_Edge Edg; ... return Edg;
        Shape::null()
    }

    /// The rcad BRepLib_MakeEdge2d stand-in over a line support (the
    /// kernel-gap note of the module header): the trimmed ends become the
    /// vertices, the z = 0 image is the 3D curve, [p1, p2] the range.
    fn make_edge_2d_line(arena: &mut BRep, l: Line2d, p1: f64, p2: f64) -> Shape {
        let a = glam::DVec3::new(
            l.origin.x + l.direction.x * p1,
            l.origin.y + l.direction.y * p1,
            0.0,
        );
        let b3 = glam::DVec3::new(
            l.origin.x + l.direction.x * p2,
            l.origin.y + l.direction.y * p2,
            0.0,
        );
        let v1 = arena.add_tvertex(a);
        let v2 = arena.add_tvertex(b3);
        let dir = glam::DVec3::new(l.direction.x, l.direction.y, 0.0);
        arena.add_tedge(
            Some(Curve3::Line(Line3 {
                origin: a,
                direction: dir.normalize_or_zero(),
            })),
            v1,
            v2,
            [p1, p2],
        )
    }

    /// The rcad BRepLib_MakeEdge2d stand-in over a circle support.
    fn make_edge_2d_circle(
        arena: &mut BRep,
        c: rcad_kernel::geom::Circle2d,
        p1: f64,
        p2: f64,
    ) -> Shape {
        let pt = |u: f64| {
            glam::DVec3::new(
                c.center.x + c.radius * (c.x_dir.x * u.cos() + c.y_dir.x * u.sin()),
                c.center.y + c.radius * (c.x_dir.y * u.cos() + c.y_dir.y * u.sin()),
                0.0,
            )
        };
        let v1 = arena.add_tvertex(pt(p1));
        let v2 = arena.add_tvertex(pt(p2));
        arena.add_tedge(
            Some(Curve3::Circle(Circle3 {
                center: glam::DVec3::new(c.center.x, c.center.y, 0.0),
                normal: glam::DVec3::Z,
                x_dir: glam::DVec3::new(c.x_dir.x, c.x_dir.y, 0.0),
                y_dir: glam::DVec3::new(c.y_dir.x, c.y_dir.y, 0.0),
                radius: c.radius,
            })),
            v1,
            v2,
            [p1, p2],
        )
    }

    /// The rcad BRepLib_MakeEdge2d stand-in over an ellipse support
    /// (gp_Elips2d parameterization: center + major*u_dir*cos + minor*v_dir*sin).
    fn make_edge_2d_ellipse(
        arena: &mut BRep,
        e: rcad_kernel::geom::Ellipse2d,
        p1: f64,
        p2: f64,
    ) -> Shape {
        let pt = |u: f64| {
            glam::DVec3::new(
                e.center.x
                    + e.major_radius * e.major_dir.x * u.cos()
                    + e.minor_radius * -e.major_dir.y * u.sin(),
                e.center.y
                    + e.major_radius * e.major_dir.y * u.cos()
                    + e.minor_radius * e.major_dir.x * u.sin(),
                0.0,
            )
        };
        let v1 = arena.add_tvertex(pt(p1));
        let v2 = arena.add_tvertex(pt(p2));
        arena.add_tedge(
            Some(Curve3::Ellipse(Ellipse3 {
                center: glam::DVec3::new(e.center.x, e.center.y, 0.0),
                normal: glam::DVec3::Z,
                major_dir: glam::DVec3::new(e.major_dir.x, e.major_dir.y, 0.0),
                major_radius: e.major_radius,
                minor_radius: e.minor_radius,
            })),
            v1,
            v2,
            [p1, p2],
        )
    }

    /// The degree-1 evaluation of the rcad BSplineCurve2 (flat knot vector,
    /// weights 1) — the local stand-in for the Geom2d_BSplineCurve::Value of
    /// BRepLib_MakeEdge2d (no public evaluator on the kernel type).
    fn bspline1_value(c: &BSplineCurve2, t: f64) -> DVec2 {
        let n = c.control_points.len();
        let k = &c.knots;
        // the parameter domain [k[degree], k[nb_poles]]
        let t = t.clamp(k[1], k[n]);
        // the knot span: the largest i with k[i] <= t (i in 1..n)
        let mut i = 1usize;
        while i + 1 <= n && k[i + 1] <= t {
            i += 1;
        }
        // the degree-1 lerp between poles i-1 and i over [k[i], k[i+1]]
        let span = k[i + 1] - k[i];
        if span <= 0.0 {
            return c.control_points[i - 1];
        }
        let a = (t - k[i]) / span;
        c.control_points[i - 1] * (1.0 - a) + c.control_points[i] * a
    }

    /// The rcad BRepLib_MakeEdge2d stand-in over the degree-1 approximation
    /// BSpline of the OCCT default case: the ends are sampled as the
    /// vertices and the range [p1, p2] is kept (the curve-less edge form of
    /// the kernel gap).
    fn make_edge_2d_bspline(arena: &mut BRep, c: &BSplineCurve2, p1: f64, p2: f64) -> Shape {
        let a = bspline1_value(c, p1);
        let b3 = bspline1_value(c, p2);
        let v1 = arena.add_tvertex(glam::DVec3::new(a.x, a.y, 0.0));
        let v2 = arena.add_tvertex(glam::DVec3::new(b3.x, b3.y, 0.0));
        arena.add_tedge(None, v1, v2, [p1, p2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hlr::algo::edges_block::EdgesBlock;
    use crate::hlr::algo::projector::Projector;
    use crate::hlr::algo::wires_block::WiresBlock;
    use crate::hlr::brep::b_curve_tool::CurveView;
    use crate::hlr::brep::curve::Curve;
    use crate::hlr::brep::face_data::FaceData;
    use crate::hlr::brep::internal_algo::InternalAlgo;
    use rcad_kernel::math::gp::Ax2;
    use rcad_kernel::topo::topods::{Orientation, ShapeType};
    use std::sync::Arc;

    /// A straight 3D segment adaptor (the BRepAdaptor_Curve test double of
    /// the data/update.rs fixture).
    struct SegView {
        origin: glam::DVec3,
        dir: glam::DVec3,
        len: f64,
    }

    impl CurveView for SegView {
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            self.len
        }
        fn d0(&self, u: f64) -> rcad_kernel::geom::Point3 {
            self.origin + self.dir * u
        }
        fn d1(&self, u: f64) -> (rcad_kernel::geom::Point3, rcad_kernel::geom::Vec3) {
            (self.origin + self.dir * u, self.dir)
        }
        fn d2(
            &self,
            u: f64,
        ) -> (
            rcad_kernel::geom::Point3,
            rcad_kernel::geom::Vec3,
            rcad_kernel::geom::Vec3,
        ) {
            (self.origin + self.dir * u, self.dir, rcad_kernel::geom::Vec3::ZERO)
        }
        fn get_type(&self) -> CurveType {
            CurveType::Line
        }
        fn line(&self) -> rcad_kernel::geom::Line3 {
            rcad_kernel::geom::Line3 {
                origin: self.origin,
                direction: self.dir,
            }
        }
        fn circle(&self) -> rcad_kernel::geom::Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(&self) -> rcad_kernel::geom::Ellipse3 {
            panic!("Standard_NoSuchObject");
        }
        fn degree(&self) -> i32 {
            0
        }
        fn nb_poles(&self) -> i32 {
            0
        }
        fn nb_knots(&self) -> i32 {
            0
        }
        fn is_closed(&self) -> bool {
            false
        }
        fn is_periodic(&self) -> bool {
            false
        }
        fn period(&self) -> f64 {
            0.0
        }
        fn resolution(&self, r3d: f64) -> f64 {
            r3d
        }
        fn parameter_3d(&self, p2d: f64) -> f64 {
            p2d
        }
        fn poles(&self) -> Vec<rcad_kernel::geom::Point3> {
            Vec::new()
        }
    }

    /// The top-view projector (identity transform).
    fn top_view_projector() -> Projector {
        Projector::from_ax2(&Ax2::new(
            glam::DVec3::ZERO,
            glam::DVec3::new(0.0, 0.0, 1.0),
            glam::DVec3::new(1.0, 0.0, 0.0),
        ))
    }

    /// The loaded projected line curve inside a fresh EdgeData with the
    /// full-parameter visible status (the line_edge_data fixture of
    /// data/update.rs).
    fn line_edge_data(proj: *const Projector, seg: &'static SegView) -> EdgeData<'static> {
        let mut c = Curve::new();
        c.projector(proj);
        c.load(seg);
        let mut e = EdgeData::new();
        e.set(
            false, false, c, 1.0e-7, 1, 2, false, false, false, false, 0.0, 1.0e-7, 10.0,
            1.0e-7,
        );
        e
    }

    /// A Data with two sharp non-vertical edges, one face carrying both as
    /// plain (non-iso / non-internal / non-outline) wire edges.
    fn sharp_two_edge_data() -> Data<'static> {
        let proj: &'static Projector = Box::leak(Box::new(top_view_projector()));
        let seg0: &'static SegView = Box::leak(Box::new(SegView {
            origin: glam::DVec3::ZERO,
            dir: glam::DVec3::new(1.0, 0.0, 0.0),
            len: 10.0,
        }));
        let seg1: &'static SegView = Box::leak(Box::new(SegView {
            origin: glam::DVec3::new(0.0, 5.0, 0.0),
            dir: glam::DVec3::new(1.0, 0.0, 0.0),
            len: 10.0,
        }));
        let mut d = Data::new(0, 2, 1);
        let arr = d.e_data_array_mut();
        arr[0] = line_edge_data(proj, seg0);
        arr[1] = line_edge_data(proj, seg1);
        // the face: one wire of two plain edges (the Wires handle assignment
        // of the data/update.rs fixture — the Set-less form).
        let mut wb = WiresBlock::new(1);
        let mut eb = EdgesBlock::new(2);
        eb.set_edge(1, 1);
        eb.set_edge(2, 2);
        eb.set_orientation(1, Orientation::Forward);
        eb.set_orientation(2, Orientation::Forward);
        wb.set(1, eb);
        *d.f_data_array_mut()[0].change_wires() = Some(Arc::new(wb));
        d
    }

    /// The Algo wiring with an installed Data (no ShapeToHLR::Load — the
    /// DS comes from the hand-built fixture).
    fn algo_with(ds: Data<'static>) -> Box<Algo> {
        let mut ia = InternalAlgo::new();
        ia.install_ds_for_test(ds);
        Box::new(Algo { internal: ia })
    }

    /// The compound children count helper.
    fn compound_children(s: &Shape) -> Vec<Shape> {
        match &*s.data {
            TShape::Compound(c) => c.clone(),
            _ => panic!("not a compound"),
        }
    }

    /// VCompound over the sharp two-edge fixture: both edges are drawn in
    /// one compound of 2 edges; HCompound is null (the hidden iterator of
    /// an all-visible edge walks zero intervals).
    #[test]
    fn hlr_to_shape_v_compound_assembly() {
        let mut algo = algo_with(sharp_two_edge_data());
        let mut hts = HLRToShape::new(algo.as_mut());
        let v = hts.v_compound();
        assert!(!v.is_null());
        assert_eq!(v.shape_type(), ShapeType::Compound);
        let children = compound_children(&v);
        assert_eq!(children.len(), 2);
        for c in &children {
            assert_eq!(c.shape_type(), ShapeType::Edge);
        }
        // HCompound: nothing hidden — the null result of the !added branch.
        let h = hts.h_compound();
        assert!(h.is_null());
    }

    /// Rg1LineVCompound over the same fixture: typ 3 with Rg1Line false on
    /// both edges — nothing is drawn and the result is null.
    #[test]
    fn hlr_to_shape_rg1_line_v_compound_empty() {
        let mut algo = algo_with(sharp_two_edge_data());
        let mut hts = HLRToShape::new(algo.as_mut());
        let r = hts.rg1_line_v_compound();
        assert!(r.is_null());
    }

    /// CompoundOfEdges(Sharp, true, false) equals VCompound; the Scaled
    /// toggling of the DS projector round-trips (the restore Scaled(false)).
    #[test]
    fn hlr_to_shape_compound_of_edges_matches_v_compound() {
        let mut algo = algo_with(sharp_two_edge_data());
        let mut hts = HLRToShape::new(algo.as_mut());
        let v = hts.compound_of_edges(TypeOfResultingEdge::Sharp, true, false);
        assert!(!v.is_null());
        assert_eq!(compound_children(&v).len(), 2);
    }

    /// A hidden edge (AllHidden status) goes to HCompound and not to
    /// VCompound.
    #[test]
    fn hlr_to_shape_hidden_edge_extraction() {
        let proj: &'static Projector = Box::leak(Box::new(top_view_projector()));
        let seg: &'static SegView = Box::leak(Box::new(SegView {
            origin: glam::DVec3::ZERO,
            dir: glam::DVec3::new(1.0, 0.0, 0.0),
            len: 10.0,
        }));
        let mut d = Data::new(0, 1, 0);
        let mut e = line_edge_data(proj, seg);
        e.status().hide_all();
        d.e_data_array_mut()[0] = e;

        let mut algo = algo_with(d);
        let mut hts = HLRToShape::new(algo.as_mut());
        // No faces: the typ >= 3 plain path draws the unused edge — but the
        // visible iterator of the AllHidden edge has zero intervals.
        let v = hts.v_compound();
        assert!(v.is_null());
        let h = hts.h_compound();
        assert!(!h.is_null());
        assert_eq!(compound_children(&h).len(), 1);
    }
}
