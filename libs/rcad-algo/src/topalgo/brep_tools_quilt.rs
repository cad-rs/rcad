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
//! 3. `BRep_Builder` -> the `BRep` pool + its in-place mutators, with the
//!    owning pool as the leading `brep` argument (the
//!    topalgo/brep_tools_substitution.rs convention "the owning pool is the
//!    leading `brep` argument"; the brep_offset_api_thru_sections_b.rs L422
//!    in-pool form of the shell maker).  Every TShape this tool creates (the
//!    copies, wires, faces, the shells and the result compound) is a slot of
//!    the caller's pool, so the products are pool-registered instead of
//!    pool-free (`Shape::index == usize::MAX` reads back as null through
//!    `Shape::is_null()` and panics the pool readers — pit 19).
//!    A source shape that is not a slot of `brep` is materialized first
//!    (`BRep::import_shape_tree`, the cross-arena bridge) so every product
//!    resolves through the pool.
//! 4. `TopoDS_Iterator(S, cumOri = false)` -> `sub_shapes(&oriented(s,
//!    Forward))` (compose with FORWARD is the identity; the rcad sub_shapes
//!    always composes the parent orientation, the OCCT default cumOri = true
//!    maps to the plain `sub_shapes(s)`).
//! 5. Shell handles: OCCT `SH` / `M(E)` / `MF(F)` / the result compound alias
//!    the live shell TShape handles, and every copy observes an in-place
//!    `BRep_Builder` edit.  The rcad pool slot + the Shape's Arc is that shared
//!    TShape, so the shells are carried as ordinary pool-registered Shapes and
//!    the former local `shell_registry` mirror (and its index reference) is
//!    gone — the pool is the registry.

use std::sync::Arc;

use indexmap::IndexMap;
use rcad_kernel::geom::{Curve2d, Curve2dEval};
use rcad_kernel::topo::topods::{tshape_flags, BRep, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::tool as bat;

/// OCCT NCollection_IndexedDataMap / DataMap of (shape -> shape) keyed by
/// TopTools_ShapeMapHasher; the payload is (key shape, item).
type BoundsMap = IndexMap<bat::ShapeKey, (Shape, Shape)>;

/// OCCT NCollection_DataMap<TopoDS_Shape, TopoDS_Shape> — the M / MF maps of
/// Shells(): key shape -> the oriented shell handle.
type ShellMap = IndexMap<bat::ShapeKey, (Shape, Shape)>;

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

/// OCCT NCollection_DataMap::Bind of a shape -> shell wrapper — replaces the
/// value when the key exists (IndexMap::insert keeps the insertion position).
/// The OCCT value is the oriented shell handle (`M(E)` / `MF(F)` return a
/// `TopoDS_Shape`), so the rcad item is the oriented shell Shape.
fn shell_map_bind(m: &mut ShellMap, k: &Shape, shell: &Shape) {
    m.insert(bat::shape_key(k), (k.clone(), shell.clone()));
}

// ---------------------------------------------------------------------------
// The pool bridge and the BRep_Builder carriers over the owning pool.
// ---------------------------------------------------------------------------

/// rcad pool test (architecture difference): true when `r` resolves to its own
/// TShape inside `brep` — the index must be a slot AND carry the same Arc.  A
/// shape built in another arena can carry an in-range index that aliases an
/// unrelated TShape (pit 2), and a pool-free shape carries usize::MAX.
fn pool_owns_shape(brep: &BRep, r: &Shape) -> bool {
    r.index < brep.tshapes.len() && Arc::ptr_eq(&brep.tshapes[r.index], &r.data)
}

/// The OCCT handle semantics (the TopoDS_Shape carries its TShape) over the
/// rcad pool: the shape itself when it is a slot of `brep`, otherwise the
/// same tree materialized into the pool (`BRep::import_shape_tree`).
fn pool_shape(brep: &mut BRep, r: &Shape) -> Shape {
    if pool_owns_shape(brep, r) {
        r.clone()
    } else {
        brep.import_shape_tree(r)
    }
}

/// OCCT TopoDS_Shape::EmptyCopied over the owning pool — the fresh TShape is a
/// slot of `brep` (`BRep::empty_copy`), so the copy is pool-registered.
fn builder_empty_copied(brep: &mut BRep, r: &Shape) -> Shape {
    let r = pool_shape(brep, r);
    brep.empty_copied(&r)
}

/// OCCT BRep_Builder::MakeShell() — the in-pool form (the
/// brep_offset_api_thru_sections_b.rs L422 helper).
fn builder_make_shell(brep: &mut BRep) -> Shape {
    brep.add_tshell(vec![])
}

/// OCCT BRep_Builder::MakeWire(W).
fn builder_make_wire(brep: &mut BRep) -> Shape {
    brep.add_twire(vec![])
}

/// OCCT BRep_Builder::MakeCompound(C).
fn builder_make_compound(brep: &mut BRep) -> Shape {
    brep.add_tcompound(vec![])
}

/// OCCT BRep_Builder::Add(Compound, S) — the in-place compound edit (the
/// BRep::add_to_compound sibling; every handle of the compound observes the
/// added child, as BRep_Builder::Add does).
fn builder_add_compound_shape(brep: &mut BRep, the_c: &Shape, the_s: &Shape) {
    rcad_kernel::topo::topods::BRepBuilder::new().add_to_compound(
        brep,
        the_c.clone(),
        the_s.clone(),
    );
}

/// OCCT BRep_Builder::Remove(Compound, S) (TopoDS_Builder.cxx L106-135).
fn builder_remove_compound_shape(brep: &mut BRep, the_c: &Shape, the_s: &Shape) {
    rcad_kernel::topo::topods::BRepBuilder::new().remove_from_compound(
        brep,
        the_c.clone(),
        the_s.clone(),
    );
}

/// OCCT BRep_Builder::Add(Shell, Face) — the in-place shell edit through the
/// pool slot (the BRep::shell_mut identity contract; Arc::make_mut would
/// strand the face on a fork).  The OCCT call sites pass the FORWARD-oriented
/// parent (arefShape), so the stored child orientation is the child's own
/// (TopoDS_Builder::Add composes Reverse(parent.Ori) into the child; FORWARD
/// is the identity).
fn builder_add_shell_face(brep: &mut BRep, the_shell: &Shape, the_face: &Shape) {
    let sd = brep.shell_mut(the_shell.clone());
    sd.faces.push(the_face.clone());
    sd.my_shapes.push(the_face.clone());
}

/// OCCT BRep_Builder::Add(W, E) — the in-place wire edit (see
/// [`builder_add_shell_face`]).
fn builder_add_wire_edge(brep: &mut BRep, the_w: &Shape, the_e: &Shape) {
    let wd = brep.wire_mut(the_w.clone());
    wd.edges.push(the_e.clone());
}

/// OCCT BRep_Builder::Add(F, W) — the first wire is the outer wire, the
/// following ones are inner wires (BRep_Builder.cxx Add Face branch;
/// tool.rs::builder_add_face_wire in its pool form).
fn builder_add_face_wire(brep: &mut BRep, the_f: &Shape, the_w: &Shape) {
    let fd = brep.face_mut(the_f.clone());
    if fd.outer_wire.is_null() {
        fd.outer_wire = the_w.clone();
    } else {
        fd.inner_wires.push(the_w.clone());
    }
}

/// OCCT BRep_Builder::Add(E, V) — attach the vertex by orientation
/// (FORWARD -> first, REVERSED -> last; tool.rs::builder_add_edge_vertex in
/// its pool form).
fn builder_add_edge_vertex(brep: &mut BRep, the_e: &Shape, the_v: &Shape) {
    let ed = brep.edge_mut_inplace(the_e.clone());
    match the_v.orientation {
        Orientation::Reversed => ed.last = the_v.clone(),
        _ => ed.first = the_v.clone(),
    }
}

/// OCCT BRep_Builder::UpdateEdge(E, C2d, F, Tol) — bind the pcurve on the face
/// (the pcurve range is the edge 3D range).  The OCCT shared-TShape edit goes
/// through the pool slot, so `myBounds.FindFromKey(E)` observes it (no fork
/// write-back); a source edge assembled outside the pool keeps the Arc form
/// (the shape travels in the caller's handle).
fn builder_update_edge_pcurve(
    brep: &mut BRep,
    the_e: &mut Shape,
    the_c2d: &Curve2d,
    the_f: &Shape,
    the_tol: f64,
) {
    if pool_owns_shape(brep, the_e) {
        let ed = brep.edge_mut_inplace(the_e.clone());
        // OCCT static UpdateCurves two-step rule (BRep_Builder.cxx L104-167).
        let [f0, l0] =
            rcad_kernel::topods::update_curves_range(the_c2d.default_domain(), ed);
        ed.pcurves.insert(bat::shape_key(the_f), (the_c2d.clone(), f0, l0));
        ed.tolerance = ed.tolerance.max(the_tol);
    } else {
        bat::builder_update_edge_pcurve(the_e, the_c2d, the_f, the_tol);
    }
}

/// OCCT BRep_Builder::Range(E, First, Last).
fn builder_range_edge(brep: &mut BRep, the_e: &mut Shape, the_first: f64, the_last: f64) {
    if pool_owns_shape(brep, the_e) {
        brep.edge_mut_inplace(the_e.clone()).range = [the_first, the_last];
    } else {
        bat::builder_range_edge(the_e, the_first, the_last);
    }
}

/// OCCT BRep_Builder::Range(E, F, First, Last) — the pcurve range of the edge
/// on the face (BRep_Builder.cxx Range CurveOnSurface branch).
fn builder_range_edge_on_face(
    brep: &mut BRep,
    the_e: &mut Shape,
    the_f: &Shape,
    the_first: f64,
    the_last: f64,
) {
    let key = bat::shape_key(the_f);
    if pool_owns_shape(brep, the_e) {
        let ed = brep.edge_mut_inplace(the_e.clone());
        if let Some(entry) = ed.pcurves.get_mut(&key) {
            entry.1 = the_first;
            entry.2 = the_last;
        }
    } else {
        bat::builder_range_edge_on_face(the_e, the_f, the_first, the_last);
    }
}

/// A per-kind shape-flag write through the pool slot (OCCT
/// TopoDS_Shape::Closed / BRep_Builder::Free / TopoDS_Shape::Orientable —
/// tool.rs::builder_set_closed / builder_set_free / builder_set_orientable in
/// their pool form).  The Shells() call sites pass the shells this tool just
/// created, so the target is always a slot of `brep`; a shape from outside the
/// pool is materialized first (`pool_shape`) rather than aliased by index.
fn builder_set_flag(brep: &mut BRep, the_s: &Shape, flag: u16, on: bool) {
    let the_s = pool_shape(brep, the_s);
    let flags = match &*the_s.data {
        TShape::Vertex(_) => &mut brep.vertex_mut(the_s.clone()).flags,
        TShape::Edge(_) => &mut brep.edge_mut_inplace(the_s.clone()).flags,
        TShape::Wire(_) => &mut brep.wire_mut(the_s.clone()).flags,
        TShape::Face(_) => &mut brep.face_mut(the_s.clone()).flags,
        TShape::Shell(_) => &mut brep.shell_mut(the_s.clone()).flags,
        TShape::Solid(_) => &mut brep.solid_mut(the_s.clone()).flags,
        TShape::CompSolid(_) | TShape::Compound(_) => return,
    };
    if on {
        *flags |= flag;
    } else {
        *flags &= !flag;
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
    fn copy_shape(brep: &mut BRep, e: &Shape, my_bounds: &mut BoundsMap) {
        // OCCT L69-71: TopoDS_Edge NE = E; NE.EmptyCopy();
        // NE.Orientation(TopAbs_FORWARD).
        let mut ne = builder_empty_copied(brep, e);
        ne.orientation = Orientation::Forward;
        // add the edges
        // OCCT L74-75: TopoDS_Iterator itv; itv.Initialize(E, false) —
        // cumOri = false (bridge 4).
        let itv = bat::sub_shapes(&bat::oriented(e, Orientation::Forward));
        for v in &itv {
            if my_bounds.contains_key(&bat::shape_key(v)) {
                let bound = my_bounds.get(&bat::shape_key(v)).unwrap().1.clone();
                builder_add_edge_vertex(brep, &ne, &bat::oriented(&bound, v.orientation));
            } else {
                builder_add_edge_vertex(brep, &ne, v);
            }
        }
        // set the 3d range (OCCT L88-91)
        let (f, l) = bat::brep_tool_range(e);
        builder_range_edge(brep, &mut ne, f, l);
        // OCCT L92: myBounds.Add(E, NE.Oriented(TopAbs_FORWARD))
        let key = bat::shape_key(e);
        my_bounds
            .entry(key)
            .or_insert((e.clone(), bat::oriented(&ne, Orientation::Forward)));
    }

    /// OCCT BRepTools_Quilt::Add(S) (BRepTools_Quilt.cxx L114-283).
    pub fn add(&mut self, brep: &mut BRep, s: &Shape) {
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
                            Self::copy_shape(brep, &e, &mut self.my_bounds);
                        }
                    }
                }
            }

            // NF will be the copy of F or F itself (OCCT L212)
            let mut nf = f.clone();

            if copy_face {
                // copy of a face (OCCT L217-220)
                nf = builder_empty_copied(brep, &f);
                nf.orientation = Orientation::Forward;

                // OCCT L222: for (TopoDS_Iterator itw(F, false); ...) — cumOri = false
                let itw = bat::sub_shapes(&bat::oriented(&f, Orientation::Forward));
                for itw_v in &itw {
                    // OCCT L224: const TopoDS_Wire& W = TopoDS::Wire(itw.Value());
                    let w = itw_v;

                    // OCCT L226-227: TopoDS_Wire NW; B.MakeWire(NW);
                    let mut nw = builder_make_wire(brep);
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
                            // handle aliases the map item, so UpdateEdge / Range
                            // edit the shared TShape the map holds — the rcad
                            // pool slot when the item is a slot of the pool.
                            // The write-back below keeps the map item in step
                            // for a source edge assembled outside the pool (the
                            // Arc form of the two range helpers).
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
                                builder_update_edge_pcurve(
                                    brep,
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
                                builder_update_edge_pcurve(
                                    brep,
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
                            builder_range_edge_on_face(brep, &mut ne, &f, u_first, u_last);
                            // the UpdateEdge / Range write-back (bridge 3)
                            self.my_bounds.get_mut(&ekey).unwrap().1 = ne.clone();
                            // OCCT L267: B.Add(NW, NE.Oriented(OE));
                            builder_add_wire_edge(brep, &nw, &bat::oriented(&ne, oe));
                        } else {
                            // OCCT L271: B.Add(NW, E);
                            builder_add_wire_edge(brep, &nw, &e);
                        }
                    }
                    // OCCT L274-275
                    nw.orientation = w.orientation;
                    builder_add_face_wire(brep, &nf, &nw);
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
    pub fn shells(&self, brep: &mut BRep) -> Shape {
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

        // OCCT L389: NCollection_DataMap M, MF — the value is the oriented
        // shell handle (bridge 5: the pool slot IS the shared TShape).
        let mut m: ShellMap = IndexMap::new();
        let mut mf: ShellMap = IndexMap::new();
        // OCCT L390-393: BRep_Builder B; TopoDS_Compound result;
        // B.MakeCompound(result); — the compound is a pool slot, edited live by
        // B.Add (L442 / L492 / L586) and B.Remove (L524).
        let result = builder_make_compound(brep);

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
                let mut sh: Option<Shape> = None;
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
                        shell_map_bind(&mut mf, &shape, &bat::oriented(sh, new_o));
                        break;
                    }
                }

                // OCCT L437: if (SH.IsNull())
                if sh.is_none() {
                    // Create a new shell, closed. Add it to the result.
                    // OCCT L440-443: B.MakeShell(SH); SH.Closed(true);
                    // B.Add(result, SH);
                    let shs = builder_make_shell(brep);
                    builder_set_flag(brep, &shs, tshape_flags::CLOSED, true);
                    // OCCT L442: B.Add(result, SH).
                    builder_add_compound_shape(brep, &result, &shs);
                    sh = Some(shs);
                    let sh = sh.as_ref().unwrap();
                    // OCCT L443: MF.Bind(Shape, SH.Oriented(Shape.Orientation()));
                    shell_map_bind(&mut mf, &shape, &bat::oriented(sh, shape.orientation));
                }

                let sh = sh.clone().unwrap();

                // Add the face to the shell (OCCT L446-450): SH.Free(true);
                // arefShape = SH.Oriented(TopAbs_FORWARD);
                // B.Add(arefShape, Shape.Oriented(MF(Shape).Orientation()));
                // — the FORWARD parent stores the face with its own
                // orientation.
                builder_set_flag(brep, &sh, tshape_flags::FREE, true);
                let mf_shape_ori = mf.get(&bat::shape_key(&shape)).unwrap().1.orientation;
                builder_add_shell_face(brep, &sh, &bat::oriented(&shape, mf_shape_ori));

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
                        if !old_shell.is_same(&sh) {
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
                            for fo in bat::sub_shapes(&old_shell) {
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
                                shell_map_bind(
                                    &mut mf,
                                    &fo,
                                    &bat::oriented(&sh, new_o_fo),
                                );
                                // OCCT L491-492: arefShapeFo = SH.Oriented(FORWARD);
                                // B.Add(arefShapeFo, Fo.Oriented(NewOFo));
                                builder_add_shell_face(brep, &sh, &bat::oriented(&fo, new_o_fo));
                            }
                            // Rebind the free edges of the old shell to the new shell
                            // gka BUG 6491 (OCCT L494-522): TopExp_Explorer
                            // aexp(SH, TopAbs_EDGE) — the SH wrapper
                            // orientation composes in.
                            let aexp = bat::explorer(
                                &bat::oriented(&sh, sh.orientation),
                                ShapeType::Edge,
                                ShapeType::Shape,
                            );
                            for ae in &aexp {
                                if !m.contains_key(&bat::shape_key(ae)) {
                                    continue;
                                }
                                let a_s = m.get(&bat::shape_key(ae)).unwrap().1.clone();
                                if a_s.is_same(&old_shell) {
                                    // update the orientation of free edges in SH.
                                    // OCCT L510-518
                                    let new_o2 = if rev {
                                        bat::top_abs_reverse(a_s.orientation)
                                    } else {
                                        a_s.orientation
                                    };

                                    // OCCT L520: M.Bind(ae, SH.Oriented(NewO));
                                    shell_map_bind(&mut m, ae, &bat::oriented(&sh, new_o2));
                                }
                            }
                            // remove the old shell from the result (OCCT L524)
                            let a_remove = bat::oriented(&old_shell, Orientation::Forward);
                            builder_remove_compound_shape(brep, &result, &a_remove);
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
                            builder_set_flag(brep, &sh, tshape_flags::ORIENTABLE, false);
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
                        shell_map_bind(&mut m, e, &bat::oriented(&sh, new_o2));
                    }
                }

                // freeze the shell (OCCT L560)
                builder_set_flag(brep, &sh, tshape_flags::FREE, false);
            } else {
                // OCCT L564: MapOtherShape.Add(Shape);
                set_add(&mut map_other_shape, &shape);
            }
        }

        // Unclose all shells having free edges (OCCT L570-577): S.Closed(false)
        // edits the shell TShape the result compound aliases (the pool slot).
        for (_, (_, shref)) in m.iter() {
            builder_set_flag(brep, shref, tshape_flags::CLOSED, false);
        }

        // gka version for free edges (OCCT L579-588)
        for (_, key_shape) in map_other_shape.iter() {
            if !set_contains(&edges_faces, key_shape)
                && self.my_bounds.contains_key(&bat::shape_key(key_shape))
            {
                let a_sh = self.my_bounds.get(&bat::shape_key(key_shape)).unwrap().1.clone();
                // OCCT L586: B.Add(result, aSh);
                builder_add_compound_shape(brep, &result, &a_sh);
            }
        }

        // OCCT L589: return result;
        result
    }
}
