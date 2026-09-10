//! OCCT BRepTools_Quilt (TKBRep/BRepTools — `BRepTools_Quilt.hxx`
//! + `BRepTools_Quilt.cxx` L39-590) — a tool to glue faces at common edges
//! and reconstruct shells (the Glue form consumed by
//! BRepOffset_MakeOffset::IsConnectedShell / MakeThickSolid / MakeShells).
//!
//! Architecture bridges (the brep_offset_make_offset.rs numbering style):
//! 1. `NCollection_IndexedDataMap<TopoDS_Shape, TopoDS_Shape,
//!    TopTools_ShapeMapHasher> myBounds` -> `IndexMap<ShapeKey, (Shape,
//!    Shape)>` (key shape + item); `Add` inserts only when the key is absent
//!    (the OCCT IndexedDataMap::Add keeps the existing binding, NCollection_
//!    IndexedDataMap.hxx L498-501 "ignored if Key was already bound").
//!    `NCollection_DataMap` (M / MF of Shells) uses the same carrier with
//!    insert-replaces (DataMap::Bind semantics); `NCollection_Map` -> the
//!    insertion-ordered set (the OCCT set iteration order is bucket-based
//!    and unspecified; the rcad order is the insertion order).
//! 2. `TopTools_ShapeMapHasher` -> `brep_algo::tool::shape_key` (TShape ptr
//!    + Location, orientation ignored).
//! 3. `BRep_Builder` -> the brep_algo::tool in-place `Arc::make_mut`
//!    re-hosts (the no-pool Arc TShape form; the thru_sections_b.rs L412
//!    precedent for the shell forms).
//! 4. `TopoDS_Iterator(S, cumOri = false)` -> `sub_shapes(&oriented(s,
//!    Forward))` (compose with FORWARD is the identity; the rcad sub_shapes
//!    always composes the parent orientation, the OCCT default cumOri = true
//!    maps to the plain `sub_shapes(s)`).
//! 5. Shells() handle aliasing: OCCT `SH` / `M(E)` / `MF(F)` / the result
//!    compound alias the live shell TShape handles, and `TopoDS_Shape::IsSame`
//!    identity follows the handle.  The rcad `Arc::make_mut` fork model would
//!    strand those aliases (stale face lists / broken IsSame) on the first
//!    in-place shell mutation, so the shells live in a local registry and the
//!    `SH` / `M` / `MF` / result slots reference it by index (the
//!    brep_offset_make_offset.rs arch. diff. #39 aliasing-bridge style).

use std::sync::Arc;

use indexmap::IndexMap;
use rcad_kernel::geom::Curve2dEval;
use rcad_kernel::topo::topods::{tshape_flags, Orientation, ShapeType, TShape, TShellData};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;

/// OCCT NCollection_IndexedDataMap / DataMap of (shape -> shape) keyed by
/// TopTools_ShapeMapHasher; the payload is (key shape, item).
type BoundsMap = IndexMap<bat::ShapeKey, (Shape, Shape)>;

/// OCCT NCollection_Map of shapes (insertion-ordered set carrier).
type ShapeSet = IndexMap<bat::ShapeKey, Shape>;

/// OCCT NCollection_Set::Add — inserts when absent (an existing member is
/// kept; the IndexMap replace form is set-equal).
fn set_add(s: &mut ShapeSet, sh: &Shape) {
    s.entry(bat::shape_key(sh)).or_insert_with(|| sh.clone());
}

/// OCCT NCollection_Set::Contains.
fn set_contains(s: &ShapeSet, sh: &Shape) -> bool {
    s.contains_key(&bat::shape_key(sh))
}

/// The Shells() shell wrapper — the OCCT `TopoDS_Shell SH` handle plus its
/// TopoDS_Shape orientation field, carried as a registry index (bridge 5).
#[derive(Clone)]
struct QuiltShellRef {
    shell: usize,
    orientation: Orientation,
}

/// OCCT NCollection_DataMap::Bind of a shell wrapper — replaces the value
/// when the key exists (IndexMap::insert keeps the insertion position).
fn shell_map_bind(
    m: &mut IndexMap<bat::ShapeKey, (Shape, QuiltShellRef)>,
    k: &Shape,
    shell: usize,
    orientation: Orientation,
) {
    m.insert(
        bat::shape_key(k),
        (
            k.clone(),
            QuiltShellRef {
                shell,
                orientation,
            },
        ),
    );
}

/// OCCT BRep_Builder::MakeShell() — the thru_sections_b.rs L412 form without
/// the pool (the no-pool Arc TShape).
fn builder_make_shell() -> Shape {
    Shape {
        data: Arc::new(TShape::Shell(TShellData {
            my_shapes: Vec::new(),
            flags: tshape_flags::DEFAULT,
            faces: Vec::new(),
        })),
        index: usize::MAX,
        location: 0,
        orientation: Orientation::Forward,
    }
}

/// OCCT BRep_Builder::Add(Shell, Face) — the thru_sections_b.rs L417 form.
/// The OCCT call sites pass the FORWARD-oriented parent (arefShape), so the
/// stored child orientation is the child's own (TopoDS_Builder::Add composes
/// Reverse(parent.Ori) into the child; FORWARD is the identity).
fn builder_add_shell_face(the_shell: &mut Shape, the_face: &Shape) {
    if let TShape::Shell(sd) = Arc::make_mut(&mut the_shell.data) {
        sd.faces.push(the_face.clone());
        sd.my_shapes.push(the_face.clone());
    }
}

/// OCCT TopoDS_Shape::Orientable(theIsOrientable) — the ORIENTABLE flag
/// write on the shape's own TShape (the bat::builder_set_closed form).
fn builder_set_orientable(the_s: &mut Shape, flag: bool) {
    let ts = Arc::make_mut(&mut the_s.data);
    let flags = match ts {
        TShape::Vertex(v) => &mut v.flags,
        TShape::Edge(e) => &mut e.flags,
        TShape::Wire(w) => &mut w.flags,
        TShape::Face(f) => &mut f.flags,
        TShape::Shell(sh) => &mut sh.flags,
        TShape::Solid(so) => &mut so.flags,
        TShape::CompSolid(_) | TShape::Compound(_) => return,
    };
    if flag {
        *flags |= tshape_flags::ORIENTABLE;
    } else {
        *flags &= !tshape_flags::ORIENTABLE;
    }
}

/// OCCT BRepTools_Quilt (BRepTools_Quilt.hxx L45-94).
pub struct BRepToolsQuilt {
    /// OCCT: myBounds.
    my_bounds: BoundsMap,
    /// OCCT: hasCopy.
    has_copy: bool,
}

impl BRepToolsQuilt {
    /// OCCT BRepTools_Quilt::BRepTools_Quilt() (BRepTools_Quilt.cxx L39-42).
    pub fn new() -> Self {
        BRepToolsQuilt {
            my_bounds: IndexMap::new(),
            has_copy: false,
        }
    }

    /// OCCT NCollection_IndexedDataMap::Add (hxx L498-501): binds the key
    /// only when absent — an existing binding is kept.
    fn my_bounds_add(&mut self, k: &Shape, item: Shape) {
        let key = bat::shape_key(k);
        self.my_bounds.entry(key).or_insert((k.clone(), item));
    }

    /// OCCT NeedCopied (BRepTools_Quilt.cxx L46-63) — tests if the shape must
    /// be copied, i.e. it contains a bound subshape.
    fn need_copied(the_shape: &Shape, my_bounds: &BoundsMap) -> bool {
        // test if the shape must be copied
        // i.e. it contains a bound subshape
        // OCCT L53: TopoDS_Iterator itv(theShape) — cumOri = true.
        let itv = bat::sub_shapes(the_shape);
        let mut is_copied = false;
        for v in &itv {
            if my_bounds.contains_key(&bat::shape_key(v)) {
                is_copied = true;
                break;
            }
        }
        is_copied
    }

    /// OCCT CopyShape (BRepTools_Quilt.cxx L65-93).
    fn copy_shape(e: &Shape, my_bounds: &mut BoundsMap) {
        // OCCT L69-71: TopoDS_Edge NE = E; NE.EmptyCopy();
        // NE.Orientation(TopAbs_FORWARD).
        let mut ne = bat::empty_copied(e);
        ne.orientation = Orientation::Forward;
        // add the edges
        // OCCT L74-75: TopoDS_Iterator itv; itv.Initialize(E, false) —
        // cumOri = false (bridge 4).
        let itv = bat::sub_shapes(&bat::oriented(e, Orientation::Forward));
        for v in &itv {
            if my_bounds.contains_key(&bat::shape_key(v)) {
                let bound = my_bounds.get(&bat::shape_key(v)).unwrap().1.clone();
                bat::builder_add_edge_vertex(&mut ne, &bat::oriented(&bound, v.orientation));
            } else {
                bat::builder_add_edge_vertex(&mut ne, v);
            }
        }
        // set the 3d range (OCCT L88-91)
        let (f, l) = bat::brep_tool_range(e);
        bat::builder_range_edge(&mut ne, f, l);
        // OCCT L92: myBounds.Add(E, NE.Oriented(TopAbs_FORWARD))
        let key = bat::shape_key(e);
        my_bounds
            .entry(key)
            .or_insert((e.clone(), bat::oriented(&ne, Orientation::Forward)));
    }

    /// OCCT BRepTools_Quilt::Add(S) (BRepTools_Quilt.cxx L114-283).
    pub fn add(&mut self, s: &Shape) {
        // Binds all the faces of S
        //  - to the face itself if it is not copied
        //  - to the copy if it is copied
        if self.my_bounds.contains_key(&bat::shape_key(s)) {
            return;
        }

        // OCCT L126-129: for (TopExp_Explorer wex(S, TopAbs_WIRE, TopAbs_FACE))
        for wex in bat::explorer(s, ShapeType::Wire, ShapeType::Face) {
            self.my_bounds_add(&wex, wex.clone());
        }

        // OCCT L131-134: for (TopExp_Explorer eex(S, TopAbs_EDGE, TopAbs_WIRE))
        for eex in bat::explorer(s, ShapeType::Edge, ShapeType::Wire) {
            self.my_bounds_add(&eex, eex.clone());
        }

        // OCCT L136-139: for (TopExp_Explorer vex(S, TopAbs_VERTEX, TopAbs_EDGE))
        for vex in bat::explorer(s, ShapeType::Vertex, ShapeType::Edge) {
            self.my_bounds_add(&vex, vex.clone());
        }

        // explore the faces (OCCT L142)
        for fex in bat::explorer(s, ShapeType::Face, ShapeType::Shape) {
            // explore the edges of the face and try to copy them
            // if one edge is bound the face must be copied

            let mut copy_face = false;
            // OCCT L149: const TopoDS_Face& F = TopoDS::Face(fex.Current());
            let f = fex;

            if self.has_copy {
                // if their is no binding, do not test for copy
                // OCCT L154: for (TopExp_Explorer fed(F, TopAbs_EDGE))
                for fed in bat::explorer(&f, ShapeType::Edge, ShapeType::Shape) {
                    if self.my_bounds.contains_key(&bat::shape_key(&fed)) {
                        copy_face = true;
                    } else {
                        // test if the edge must be copied
                        // i.e. it contains a bound vertex

                        let copy_edge = Self::need_copied(&fed, &self.my_bounds);
                        // OCCT L168: const TopoDS_Edge& E = TopoDS::Edge(fed.Current());
                        let e = fed;

                        if copy_edge {
                            // copy of an edge
                            copy_face = true;
                            Self::copy_shape(&e, &mut self.my_bounds);
                        }
                    }
                }
            }

            // NF will be the copy of F or F itself (OCCT L212)
            let mut nf = f.clone();

            if copy_face {
                // copy of a face (OCCT L217-220)
                nf = bat::empty_copied(&f);
                nf.orientation = Orientation::Forward;

                // OCCT L222: for (TopoDS_Iterator itw(F, false); ...) — cumOri = false
                let itw = bat::sub_shapes(&bat::oriented(&f, Orientation::Forward));
                for itw_v in &itw {
                    // OCCT L224: const TopoDS_Wire& W = TopoDS::Wire(itw.Value());
                    let w = itw_v;

                    // OCCT L226-227: TopoDS_Wire NW; B.MakeWire(NW);
                    let mut nw = bat::builder_make_wire();
                    // OCCT L228: TopoDS_Iterator ite(W, false) — cumOri = false
                    let ite = bat::sub_shapes(&bat::oriented(w, Orientation::Forward));
                    let mut u_first = 0.0;
                    let mut u_last = 0.0;

                    // Reconstruction of wires. (OCCT L233)
                    for ite_v in &ite {
                        // OCCT L235-236
                        let e = ite_v.clone();
                        let mut oe = e.orientation;
                        if self.my_bounds.contains_key(&bat::shape_key(&e)) {
                            // OCCT L239: const TopoDS_Edge& NE =
                            // TopoDS::Edge(myBounds.FindFromKey(E)); the OCCT
                            // handle aliases the map item and UpdateEdge
                            // mutates it in place; the rcad fork is written
                            // back to the map slot below (bridge 1).
                            let ekey = bat::shape_key(&e);
                            let mut ne = self.my_bounds.get(&ekey).unwrap().1.clone();
                            // pcurve.
                            if ne.orientation == Orientation::Forward {
                                // OCCT L243-246: B.UpdateEdge(NE,
                                // BRep_Tool::CurveOnSurface(E, F, UFirst, ULast),
                                // F, BRep_Tool::Tolerance(E));
                                let (c2d, cf, cl) = bat::brep_tool_curve_on_surface(&e, &f)
                                    .expect("BRep_Tool::CurveOnSurface null");
                                u_first = cf;
                                u_last = cl;
                                bat::builder_update_edge_pcurve(
                                    &mut ne,
                                    &c2d,
                                    &f,
                                    bat::brep_tool_tolerance(&e),
                                );
                            } else {
                                // If NE is REVERSED
                                // => the 3D curves do not have the same orientation.
                                // (This is a convention, cf. BRepTools_Quilt.cdl and the Bind
                                // method.)
                                // => the PCurve of E on F must be reversed.
                                // OCCT L256-262
                                oe = bat::top_abs_reverse(oe);
                                let (ce, cf, cl) = bat::brep_tool_curve_on_surface(&e, &f)
                                    .expect("BRep_Tool::CurveOnSurface null");
                                u_first = cf;
                                u_last = cl;
                                let nce = rcad_kernel::geom::reverse_curve2d(&ce);
                                bat::builder_update_edge_pcurve(
                                    &mut ne,
                                    &nce,
                                    &f,
                                    bat::brep_tool_tolerance(&e),
                                );
                                let tmp = u_first;
                                u_first = ce.reversed_parameter(u_last);
                                u_last = ce.reversed_parameter(tmp);
                            }
                            // pcurve range (OCCT L265)
                            bat::builder_range_edge_on_face(&mut ne, &f, u_first, u_last);
                            // the UpdateEdge / Range write-back (bridge 1)
                            self.my_bounds.get_mut(&ekey).unwrap().1 = ne.clone();
                            // OCCT L267: B.Add(NW, NE.Oriented(OE));
                            bat::builder_add_wire_edge(&mut nw, &bat::oriented(&ne, oe));
                        } else {
                            // OCCT L271: B.Add(NW, E);
                            bat::builder_add_wire_edge(&mut nw, &e);
                        }
                    }
                    // OCCT L274-275
                    nw.orientation = w.orientation;
                    bat::builder_add_face_wire(&mut nf, &nw);
                }
                // OCCT L277
                nf.orientation = f.orientation;
            }

            // binds the face to itself or its copy (OCCT L281)
            self.my_bounds_add(&f, nf);
        }
    }

    /// OCCT BRepTools_Quilt::Bind(Vold, Vnew) (BRepTools_Quilt.cxx L287-293).
    pub fn bind_vertex(&mut self, vold: &Shape, vnew: &Shape) {
        if !self.my_bounds.contains_key(&bat::shape_key(vold)) {
            self.my_bounds_add(vold, vnew.clone());
        }
    }

    /// OCCT BRepTools_Quilt::Bind(Eold, Enew) (BRepTools_Quilt.cxx L297-345).
    pub fn bind_edge(&mut self, eold: &Shape, enew: &Shape) {
        if !self.my_bounds.contains_key(&bat::shape_key(eold)) {
            // OCCT L301: TopoDS_Edge ENew = Enew;
            let mut e_new = enew.clone();
            if self.is_copied(enew) {
                // OCCT L304-305: ENew = TopoDS::Edge(Copy(Enew));
                // ENew.Orientation(Enew.Orientation());
                e_new = self.copy(enew);
                e_new.orientation = enew.orientation;
            }

            // OCCT L308-315
            let key = bat::oriented(eold, Orientation::Forward);
            let item = if eold.orientation != e_new.orientation {
                bat::oriented(&e_new, Orientation::Reversed)
            } else {
                bat::oriented(&e_new, Orientation::Forward)
            };
            self.my_bounds_add(&key, item);
            // if new binding bind also the vertices
            // OCCT L317: TopoDS_Iterator itold(Eold) — cumOri = true
            let itold = bat::sub_shapes(eold);
            for itold_v in &itold {
                if !self.my_bounds.contains_key(&bat::shape_key(itold_v)) {
                    // find the vertex of Enew with same orientation
                    let an_orien = itold_v.orientation;
                    // OCCT L324: TopoDS_Iterator itnew(ENew) — cumOri = true
                    let itnew = bat::sub_shapes(&e_new);
                    for itnew_v in &itnew {
                        if itnew_v.orientation == an_orien {
                            // OCCT L329: TopoDS_Vertex VNew = TopoDS::Vertex(itnew.Value());
                            let mut v_new = itnew_v.clone();
                            if self.is_copied(&v_new) {
                                // if VNew has been copied take the copy
                                v_new = self.copy(&v_new);
                            }
                            self.my_bounds_add(itold_v, v_new);
                            break;
                        }
                    }
                }
            }
            self.has_copy = true;
        }
    }

    /// OCCT BRepTools_Quilt::IsCopied(S) (BRepTools_Quilt.cxx L349-359).
    pub fn is_copied(&self, s: &Shape) -> bool {
        if let Some((_, item)) = self.my_bounds.get(&bat::shape_key(s)) {
            !s.is_same(item)
        } else {
            false
        }
    }

    /// OCCT BRepTools_Quilt::Copy(S) (BRepTools_Quilt.cxx L363-367) —
    /// returns the shape substituted to S in the Quilt.
    pub fn copy(&self, s: &Shape) -> Shape {
        // Standard_NoSuchObject_Raise_if(!IsCopied(S), "BRepTools_Quilt::Copy");
        if !self.is_copied(s) {
            panic!("Standard_NoSuchObject: BRepTools_Quilt::Copy");
        }
        self.my_bounds.get(&bat::shape_key(s)).unwrap().1.clone()
    }

    /// OCCT BRepTools_Quilt::Shells() (BRepTools_Quilt.cxx L371-590) —
    /// returns a Compound of shells made from the current set of faces.
    pub fn shells(&self) -> Shape {
        // Outline of the algorithm
        //
        // In the map M we bind the free edges to their shells
        // We explore all the faces in myBounds
        // For each one we search the edges in the map and either :
        //
        // - Start a new shell if no edge is a free edge.
        // - Add the face to an existing shell
        // - Connect other shells if the face touch more than one shell

        // In the Map M the Shell is bound with the relative orientation of E
        // in the shell
        // In the map MF we binb the face to its shell.
        // In the Map MF the Shell is bound with the relative orientation of F
        // in the shell

        // OCCT L389: NCollection_DataMap M, MF — the registry-index carrier
        // (bridge 5).
        let mut m: IndexMap<bat::ShapeKey, (Shape, QuiltShellRef)> = IndexMap::new();
        let mut mf: IndexMap<bat::ShapeKey, (Shape, QuiltShellRef)> = IndexMap::new();
        // OCCT L390-393: BRep_Builder B; TopoDS_Compound result;
        // B.MakeCompound(result); — the rcad mirror lists (bridge 5): the
        // OCCT compound is mutated live (B.Add at L442 / B.Remove at L524);
        // the rcad assembly reproduces the same final content order
        // (surviving shells in creation order, then the other shapes).
        let mut shell_registry: Vec<Shape> = Vec::new();
        let mut result_shells: Vec<usize> = Vec::new();
        let mut result_other: Vec<Shape> = Vec::new();

        // OCCT L395-396: MapOtherShape / EdgesFaces (gka)
        let mut map_other_shape: ShapeSet = IndexMap::new();
        let mut edges_faces: ShapeSet = IndexMap::new();

        // loop on the face in myBounds (OCCT L403-405: for (int ii = 1;
        // ii <= myBounds.Extent(); ii++) ... FindFromIndex(ii) — the rcad
        // index order is the OCCT index order).
        for ii in 0..self.my_bounds.len() {
            let shape = self.my_bounds[ii].1.clone();
            if shape.shape_type() == ShapeType::Face {
                // OCCT L408-411 (gka)
                for a_exp_edg in bat::explorer(&shape, ShapeType::Edge, ShapeType::Shape) {
                    set_add(&mut edges_faces, &a_exp_edg);
                }

                // OCCT L413-414: TopoDS_Shell SH (null); TopAbs_Orientation NewO;
                let mut sh: Option<QuiltShellRef> = None;
                let mut new_o = Orientation::Forward;

                // OCCT L416-435: TopExp_Explorer itf1(Shape, TopAbs_EDGE)
                let itf1 = bat::explorer(&shape, ShapeType::Edge, ShapeType::Shape);
                for e in &itf1 {
                    if let Some((_, shref)) = m.get(&bat::shape_key(e)) {
                        // OCCT L422: SH = TopoDS::Shell(M(E));
                        sh = Some(shref.clone());
                        let sh = sh.as_ref().unwrap();
                        if sh.orientation == e.orientation {
                            new_o = bat::top_abs_reverse(shape.orientation);
                        } else {
                            new_o = shape.orientation;
                        }
                        shell_map_bind(&mut mf, &shape, sh.shell, new_o);
                        break;
                    }
                }

                // OCCT L437: if (SH.IsNull())
                if sh.is_none() {
                    // Create a new shell, closed. Add it to the result.
                    // OCCT L440-443: B.MakeShell(SH); SH.Closed(true);
                    // B.Add(result, SH);
                    let mut shs = builder_make_shell();
                    bat::builder_set_closed(&mut shs, true);
                    shell_registry.push(shs);
                    let id = shell_registry.len() - 1;
                    sh = Some(QuiltShellRef {
                        shell: id,
                        orientation: Orientation::Forward,
                    });
                    result_shells.push(id);
                    let sh = sh.as_ref().unwrap();
                    // OCCT L443: MF.Bind(Shape, SH.Oriented(Shape.Orientation()));
                    shell_map_bind(&mut mf, &shape, sh.shell, shape.orientation);
                }

                let sh = sh.as_ref().unwrap();

                // Add the face to the shell (OCCT L446-450): SH.Free(true);
                // arefShape = SH.Oriented(TopAbs_FORWARD);
                // B.Add(arefShape, Shape.Oriented(MF(Shape).Orientation()));
                // — the FORWARD parent stores the face with its own
                // orientation (bridge: builder_add_shell_face).
                bat::builder_set_free(&mut shell_registry[sh.shell], true);
                let mf_shape_ori = mf.get(&bat::shape_key(&shape)).unwrap().1.orientation;
                builder_add_shell_face(
                    &mut shell_registry[sh.shell],
                    &bat::oriented(&shape, mf_shape_ori),
                );

                // OCCT L452: TopExp_Explorer itf(Shape.Oriented(TopAbs_FORWARD),
                // TopAbs_EDGE)
                let itf = bat::explorer(
                    &bat::oriented(&shape, Orientation::Forward),
                    ShapeType::Edge,
                    ShapeType::Shape,
                );

                for e in &itf {
                    // OCCT L456: const TopoDS_Edge& E = TopoDS::Edge(itf.Current());
                    if let Some((_, old_shell)) = m.get(&bat::shape_key(e)).cloned() {
                        // OCCT L460-461: const TopoDS_Shape oldShell = M(E);
                        // if (!oldShell.IsSame(SH))
                        if old_shell.shell != sh.shell {
                            // Fuse the old shell with the new one
                            // Compare the orientation of E in SH and in oldshell.
                            // OCCT L465-471
                            let mut an_orien = e.orientation;
                            if mf.get(&bat::shape_key(&shape)).unwrap().1.orientation
                                == Orientation::Reversed
                            {
                                an_orien = bat::top_abs_reverse(an_orien);
                            }

                            // if rev = True oldShell has to be reversed.
                            let rev = an_orien == old_shell.orientation;

                            // Add the faces of oldShell in SH.
                            // OCCT L475: for (TopoDS_Iterator its(oldShell)) —
                            // cumOri = true: the faces carry the oldShell
                            // wrapper orientation composed in.
                            let old_shell_shape = bat::oriented(
                                &shell_registry[old_shell.shell],
                                old_shell.orientation,
                            );
                            for fo in bat::sub_shapes(&old_shell_shape) {
                                // OCCT L477: const TopoDS_Face Fo = TopoDS::Face(its.Value());
                                // update the orientation of Fo in SH.
                                // OCCT L479-487
                                let mf_fo_ori =
                                    mf.get(&bat::shape_key(&fo)).unwrap().1.orientation;
                                let new_o_fo = if rev {
                                    bat::top_abs_reverse(mf_fo_ori)
                                } else {
                                    mf_fo_ori
                                };

                                // OCCT L489: MF.Bind(Fo, SH.Oriented(NewOFo));
                                shell_map_bind(&mut mf, &fo, sh.shell, new_o_fo);
                                // OCCT L491-492: arefShapeFo = SH.Oriented(FORWARD);
                                // B.Add(arefShapeFo, Fo.Oriented(NewOFo));
                                builder_add_shell_face(
                                    &mut shell_registry[sh.shell],
                                    &bat::oriented(&fo, new_o_fo),
                                );
                            }
                            // Rebind the free edges of the old shell to the new shell
                            // gka BUG 6491 (OCCT L494-522): TopExp_Explorer
                            // aexp(SH, TopAbs_EDGE) — the SH wrapper
                            // orientation composes in.
                            let aexp = bat::explorer(
                                &bat::oriented(&shell_registry[sh.shell], sh.orientation),
                                ShapeType::Edge,
                                ShapeType::Shape,
                            );
                            for ae in &aexp {
                                if !m.contains_key(&bat::shape_key(ae)) {
                                    continue;
                                }
                                let a_s = m.get(&bat::shape_key(ae)).unwrap().1.clone();
                                if a_s.shell == old_shell.shell {
                                    // update the orientation of free edges in SH.
                                    // OCCT L510-518
                                    let new_o2 = if rev {
                                        bat::top_abs_reverse(a_s.orientation)
                                    } else {
                                        a_s.orientation
                                    };

                                    // OCCT L520: M.Bind(ae, SH.Oriented(NewO));
                                    shell_map_bind(&mut m, ae, sh.shell, new_o2);
                                }
                            }
                            // remove the old shell from the result (OCCT L524)
                            result_shells.retain(|&id| id != old_shell.shell);
                        }
                        // Test if SH is always orientable. (OCCT L526-536)
                        let mut an_orien = e.orientation;
                        if mf.get(&bat::shape_key(&shape)).unwrap().1.orientation
                            == Orientation::Reversed
                        {
                            an_orien = bat::top_abs_reverse(an_orien);
                        }

                        if m.get(&bat::shape_key(e)).unwrap().1.orientation == an_orien {
                            // OCCT L535: SH.Orientable(false);
                            builder_set_orientable(&mut shell_registry[sh.shell], false);
                        }

                        // remove the edge from M (no more a free edge) (OCCT L539)
                        m.shift_remove(&bat::shape_key(e));
                    } else {
                        // OCCT L541-556
                        let mut new_o2 = e.orientation;
                        if mf.get(&bat::shape_key(&shape)).unwrap().1.orientation
                            == Orientation::Reversed
                        {
                            new_o2 = bat::top_abs_reverse(new_o2);
                        }
                        // OCCT L548: if (!E.IsNull()) — the rcad explorer
                        // never yields null shapes.
                        // OCCT L550: M.Bind(E, SH.Oriented(NewO));
                        shell_map_bind(&mut m, e, sh.shell, new_o2);
                    }
                }

                // freeze the shell (OCCT L560)
                bat::builder_set_free(&mut shell_registry[sh.shell], false);
            } else {
                // OCCT L564: MapOtherShape.Add(Shape);
                set_add(&mut map_other_shape, &shape);
            }
        }

        // Unclose all shells having free edges (OCCT L570-577): the OCCT
        // S.Closed(false) mutates the shell TShape aliased by the result
        // compound; the rcad write goes through the registry slots the
        // compound is assembled from (bridge 5).
        for (_, (_, shref)) in m.iter() {
            bat::builder_set_closed(&mut shell_registry[shref.shell], false);
        }

        // gka version for free edges (OCCT L579-588)
        for (_, key_shape) in map_other_shape.iter() {
            if !set_contains(&edges_faces, key_shape)
                && self.my_bounds.contains_key(&bat::shape_key(key_shape))
            {
                let a_sh = self.my_bounds.get(&bat::shape_key(key_shape)).unwrap().1.clone();
                // OCCT L586: B.Add(result, aSh);
                result_other.push(a_sh);
            }
        }

        // OCCT L589: return result; — the deferred compound assembly
        // (bridge 5).
        let mut result = bat::builder_make_compound();
        for id in &result_shells {
            bat::builder_add_compound_shape(&mut result, &shell_registry[*id]);
        }
        for s in &result_other {
            bat::builder_add_compound_shape(&mut result, s);
        }
        result
    }
}
