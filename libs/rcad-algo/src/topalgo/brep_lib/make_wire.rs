//! OCCT BRepLib_MakeWire (TKTopAlgo/BRepLib/BRepLib_MakeWire.cxx L1-489,
//! whole file) — 1:1 translation.
//!
//! The local helpers `topexp_vertices_edge` / `topexp_vertices_wire`
//! translate OCCT TopExp::Vertices (TKBRep/TopExp/TopExp.cxx L214-247 and
//! L255-299), which the MakeWire body drives.
//!
//! Architecture adaptations (rcad kernel):
//! - `TopoDS_Iterator` over an edge -> the edge TShape children, encoded in
//!   the kernel `TEdgeData` as the `first` / `last` vertex slots.
//! - `NCollection_Map<TopoDS_Shape, TopTools_ShapeHasher>` (myVertices and
//!   the wire-form vertex map) -> `Vec<Shape>` with a linear `is_equal`
//!   scan (TopTools_ShapeMapHasher::IsEqual is the full
//!   TShape+Location+Orientation equality, i.e. the kernel `is_equal`);
//!   the Vec preserves the map insertion order for the iteration sites.
//! - `BRep_Builder::Add(anEdge, aVertex)` — OCCT appends the vertex to the
//!   edge child list; the kernel edge encodes the two ends in fixed slots,
//!   so the add fills `first` when null and `last` otherwise.
//! - `TopoDS_Shape::Closed(...)` -> the `tshape_flags::CLOSED` bit of the
//!   wire/edge TShape, read and written through the owning BRep.

use rcad_kernel::topo::topods::{tshape_flags, BRep, BRepBuilder, Orientation, Shape, ShapeType};

use glam::DVec3;

// =========================================================================
// OCCT BRepLib_WireError (BRepLib_MakeWire.hxx L27-34).
// =========================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepLibWireError {
    NoError,
    EmptyWire,
    DisconnectedWire,
    NonManifoldWire,
    WireDone,
}

// =========================================================================
// OCCT BRepLib_MakeWire (BRepLib_MakeWire.hxx L40-120).
// =========================================================================
pub struct BRepLibMakeWire {
    /// hxx: TopoDS_Wire myShape.
    my_shape: Shape,
    /// hxx: TopoDS_Edge myEdge.
    my_edge: Shape,
    /// hxx: TopoDS_Vertex myVertex.
    my_vertex: Shape,
    /// hxx: NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> myVertices.
    my_vertices: Vec<Shape>,
    /// BRepLib_MakeShape hxx: TopoDS_Vertex VF.
    vf: Shape,
    /// BRepLib_MakeShape hxx: TopoDS_Vertex VL.
    vl: Shape,
    /// hxx: BRepLib_WireError myError.
    my_error: BRepLibWireError,
    /// OCCT BRepLib_MakeShape::Done()/NotDone() flag.
    done: bool,
}

impl Default for BRepLibMakeWire {
    fn default() -> Self {
        Self::new()
    }
}

impl BRepLibMakeWire {
    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire() (L35-38).
    pub fn new() -> Self {
        BRepLibMakeWire {
            my_shape: Shape::null(),
            my_edge: Shape::null(),
            my_vertex: Shape::null(),
            my_vertices: Vec::new(),
            vf: Shape::null(),
            vl: Shape::null(),
            my_error: BRepLibWireError::EmptyWire,
            done: false,
        }
    }

    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire(const TopoDS_Edge& E)
    /// (L42-45).
    pub fn new_with_edge(brep: &mut BRep, bb: &mut BRepBuilder, e: &Shape) -> Self {
        let mut r = BRepLibMakeWire::new();
        r.add_edge(brep, bb, e);
        r
    }

    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire(E1, E2) (L49-53).
    pub fn new_with_two_edges(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        e1: &Shape,
        e2: &Shape,
    ) -> Self {
        let mut r = BRepLibMakeWire::new();
        r.add_edge(brep, bb, e1);
        r.add_edge(brep, bb, e2);
        r
    }

    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire(E1, E2, E3) (L57-64).
    pub fn new_with_three_edges(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        e1: &Shape,
        e2: &Shape,
        e3: &Shape,
    ) -> Self {
        let mut r = BRepLibMakeWire::new();
        r.add_edge(brep, bb, e1);
        r.add_edge(brep, bb, e2);
        r.add_edge(brep, bb, e3);
        r
    }

    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire(E1, E2, E3, E4) (L68-77).
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_four_edges(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        e1: &Shape,
        e2: &Shape,
        e3: &Shape,
        e4: &Shape,
    ) -> Self {
        let mut r = BRepLibMakeWire::new();
        r.add_edge(brep, bb, e1);
        r.add_edge(brep, bb, e2);
        r.add_edge(brep, bb, e3);
        r.add_edge(brep, bb, e4);
        r
    }

    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire(const TopoDS_Wire& W)
    /// (L81-84).
    pub fn new_with_wire(brep: &mut BRep, bb: &mut BRepBuilder, w: &Shape) -> Self {
        let mut r = BRepLibMakeWire::new();
        r.add_wire(brep, bb, w);
        r
    }

    /// OCCT BRepLib_MakeWire::BRepLib_MakeWire(const TopoDS_Wire& W,
    /// const TopoDS_Edge& E) (L88-92).
    pub fn new_with_wire_and_edge(
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        w: &Shape,
        e: &Shape,
    ) -> Self {
        let mut r = BRepLibMakeWire::new();
        r.add_wire(brep, bb, w);
        r.add_edge(brep, bb, e);
        r
    }

    /// OCCT BRepLib_MakeWire::Add(const TopoDS_Wire& W) (L96-106).
    pub fn add_wire(&mut self, brep: &mut BRep, bb: &mut BRepBuilder, w: &Shape) {
        // for (TopoDS_Iterator it(W); it.More(); it.Next()).
        for v in wire_edges(brep, w) {
            self.add_edge(brep, bb, &v);
            // L101-104.
            if self.my_error != BRepLibWireError::WireDone {
                break;
            }
        }
    }

    /// OCCT BRepLib_MakeWire::Add(const TopoDS_Edge& E) (L110-113).
    pub fn add_edge(&mut self, brep: &mut BRep, bb: &mut BRepBuilder, e: &Shape) {
        self.add_edge_checked(brep, bb, e, true);
    }

    /// OCCT BRepLib_MakeWire::Add(const TopoDS_Edge& E, bool
    /// IsCheckGeometryProximity) (L123-453).
    pub fn add_edge_checked(
        &mut self,
        brep: &mut BRep,
        bb: &mut BRepBuilder,
        e: &Shape,
        is_check_geometry_proximity: bool,
    ) {
        // L126: bool forward = false.
        let mut forward = false;
        // L128: bool reverse = false.
        let mut reverse = false;
        // L130: bool init = false.
        let init;

        if self.my_edge.is_null() {
            // L135-137: first edge, create the wire.
            init = true;
            // L139: B.MakeWire(TopoDS::Wire(myShape)).
            self.my_shape = bb.make_wire(brep);
            // L142: myEdge = E.
            self.my_edge = e.clone();
            // L145-148: add the vertices.
            for v in edge_vertices(brep, &self.my_edge) {
                if !self.my_vertices.iter().any(|x| x.is_equal(&v)) {
                    self.my_vertices.push(v);
                }
            }
        } else {
            // L153: init = myShape.Closed() — if it is closed no control.
            init = shape_closed(brep, &self.my_shape);
            // L154-155: EE = TopoDS::Edge(E.Oriented(TopAbs_FORWARD)).
            let mut ee = e.clone();
            ee.orientation = Orientation::Forward;

            // L160-161.
            let mut connected = false;
            let mut copyedge = false;

            // L163-169.
            if self.my_error != BRepLibWireError::NonManifoldWire {
                if self.vf.is_null() || self.vl.is_null() {
                    self.my_error = BRepLibWireError::NonManifoldWire;
                }
            }

            // L171-286: for (it.Initialize(EE); it.More(); it.Next()).
            for ve in edge_vertices(brep, &ee) {
                // L177: if the vertex is in the wire, ok for the connection.
                if self.my_vertices.iter().any(|x| x.is_equal(&ve)) {
                    connected = true;
                    self.my_vertex = ve.clone();
                    if self.my_error != BRepLibWireError::NonManifoldWire {
                        // L184: is it always so?
                        if self.vf.is_same(&self.vl) {
                            // L187: orientation indetermined (in 3d) —
                            // preserve the initial.
                            if !self.vf.is_same(&ve) {
                                self.my_error = BRepLibWireError::NonManifoldWire;
                            }
                        } else if self.vf.is_same(&ve) {
                            // L194-203.
                            if ve.orientation == Orientation::Forward {
                                reverse = true;
                            } else {
                                forward = true;
                            }
                        } else if self.vl.is_same(&ve) {
                            // L205-215.
                            if ve.orientation == Orientation::Reversed {
                                reverse = true;
                            } else {
                                forward = true;
                            }
                        } else {
                            // L217-218.
                            self.my_error = BRepLibWireError::NonManifoldWire;
                        }
                    }
                } else if is_check_geometry_proximity {
                    // L223-285: search if there is a similar vertex in the
                    // edge.
                    let pe = brep.vertex(ve.clone()).point;
                    for i in 0..self.my_vertices.len() {
                        let vw = self.my_vertices[i].clone();
                        let pw = brep.vertex(vw.clone()).point;
                        let l = pe.distance(pw);

                        // L234.
                        if (l < brep.vertex(ve.clone()).tolerance)
                            || (l < brep.vertex(vw.clone()).tolerance)
                        {
                            copyedge = true;
                            if self.my_error != BRepLibWireError::NonManifoldWire {
                                // L240: is it always so?
                                if self.vf.is_same(&self.vl) {
                                    // L243: orientation indetermined (in 3d)
                                    // — preserve the initial.
                                    if !self.vf.is_same(&vw) {
                                        self.my_error = BRepLibWireError::NonManifoldWire;
                                    }
                                } else if self.vf.is_same(&vw) {
                                    // L250-259.
                                    if ve.orientation == Orientation::Forward {
                                        reverse = true;
                                    } else {
                                        forward = true;
                                    }
                                } else if self.vl.is_same(&vw) {
                                    // L261-271.
                                    if ve.orientation == Orientation::Reversed {
                                        reverse = true;
                                    } else {
                                        forward = true;
                                    }
                                } else {
                                    // L273-274.
                                    self.my_error = BRepLibWireError::NonManifoldWire;
                                }
                            }
                            // L278: break.
                            break;
                        }
                    }
                    // L281-284.
                    if copyedge {
                        connected = true;
                    }
                }
            }

            // L288-293.
            if !connected {
                self.my_error = BRepLibWireError::DisconnectedWire;
                self.not_done();
                return;
            } else {
                if !copyedge {
                    // L298: myEdge = EE.
                    self.my_edge = ee.clone();
                    // L299-302.
                    for v in edge_vertices(brep, &ee) {
                        if !self.my_vertices.iter().any(|x| x.is_equal(&v)) {
                            self.my_vertices.push(v);
                        }
                    }
                } else {
                    // L307: copy the edge — Dummy = EE.EmptyCopied().
                    self.my_edge = brep.empty_copied(&ee);

                    // L310-367.
                    for ve in edge_vertices(brep, &ee) {
                        // L314: PE = BRep_Tool::Pnt(VE).
                        let pe = brep.vertex(ve.clone()).point;

                        // L316: bool newvertex = false.
                        let mut newvertex = false;
                        for i in 0..self.my_vertices.len() {
                            let vw = self.my_vertices[i].clone();
                            // L319-323.
                            let pw = brep.vertex(vw.clone()).point;
                            let l = pe.distance(pw);
                            let tolw = brep.vertex(vw.clone()).tolerance;
                            let tole = brep.vertex(ve.clone()).tolerance;

                            // L325.
                            if (l < tole) || (l < tolw) {
                                // L328-345.
                                let mut maxtol = 0.5 * (tolw + tole + l);
                                let cw;
                                let ce;
                                if maxtol > tolw && maxtol > tole {
                                    cw = (maxtol - tole) / l;
                                    ce = 1.0 - cw;
                                } else if maxtol > tolw {
                                    maxtol = tole;
                                    cw = 0.0;
                                    ce = 1.0;
                                } else {
                                    maxtol = tolw;
                                    cw = 1.0;
                                    ce = 0.0;
                                }

                                // L347-349.
                                let pc = DVec3::new(
                                    cw * pw.x + ce * pe.x,
                                    cw * pw.y + ce * pe.y,
                                    cw * pw.z + ce * pe.z,
                                );

                                // L351: B.UpdateVertex(VW, PC, maxtol).
                                bb.update_vertex_point(brep, vw.clone(), pc, maxtol);

                                // L353-356.
                                newvertex = true;
                                self.my_vertex = vw.clone();
                                self.my_vertex.orientation = ve.orientation;
                                builder_add_edge_vertex(brep, &self.my_edge, &self.my_vertex);
                                // L357: B.Transfert(EE, myEdge, VE, myVertex).
                                bb.transfert_vertex_param(
                                    brep,
                                    ee.clone(),
                                    ve.clone(),
                                    self.my_edge.clone(),
                                    self.my_vertex.clone(),
                                );
                                // L358: break.
                                break;
                            }
                        }
                        // L361-366.
                        if !newvertex {
                            if !self.my_vertices.iter().any(|x| x.is_equal(&ve)) {
                                self.my_vertices.push(ve.clone());
                            }
                            builder_add_edge_vertex(brep, &self.my_edge, &ve);
                            bb.transfert_vertex_param(
                                brep,
                                ee.clone(),
                                ve.clone(),
                                self.my_edge.clone(),
                                ve.clone(),
                            );
                        }
                    }
                }
            }
            // L370-379: make a decision about the orientation of the edge —
            // if there is an ambiguity (in 3d) preserve the orientation
            // given at input.
            if ((forward == reverse) && (e.orientation == Orientation::Reversed))
                || (reverse && !forward)
            {
                // L378: myEdge.Reverse().
                self.my_edge.orientation = match self.my_edge.orientation {
                    Orientation::Forward => Orientation::Reversed,
                    Orientation::Reversed => Orientation::Forward,
                    other => other,
                };
            }
        }

        // L383: add myEdge to myShape.
        bb.add_to_wire(brep, self.my_shape.clone(), self.my_edge.clone());
        // L384: myShape.Closed(false).
        set_shape_closed(brep, &self.my_shape, false);

        // L386-444: initialize VF, VL.
        if init {
            // L389: TopExp::Vertices(TopoDS::Wire(myShape), VF, VL).
            let (vf, vl) = topexp_vertices_wire(brep, &self.my_shape);
            self.vf = vf;
            self.vl = vl;
        } else {
            if self.my_error == BRepLibWireError::WireDone {
                // L394: update only.
                // L395-396: TopExp::Vertices(myEdge, V1, V2).
                let (v1, v2) = topexp_vertices_edge(brep, &self.my_edge, false);
                let mut vref = Shape::null();
                if !v1.is_null() && v1.is_same(&self.my_vertex) {
                    // L398-399.
                    vref = v2;
                } else if !v2.is_null() && v2.is_same(&self.my_vertex) {
                    // L401-403.
                    vref = v1;
                } else {
                    // L410.
                    self.my_error = BRepLibWireError::NonManifoldWire;
                }

                if self.vf.is_same(&self.vl) {
                    // L413-420: particular case — it is required to control
                    // the orientation (the OCCT debug print is compiled out).
                } else {
                    // L422: general case.
                    if self.vf.is_same(&self.my_vertex) {
                        // L424-425.
                        self.vf = vref;
                    } else if self.vl.is_same(&self.my_vertex) {
                        // L427-429.
                        self.vl = vref;
                    } else {
                        // L436.
                        self.my_error = BRepLibWireError::NonManifoldWire;
                    }
                }
            }
            // L440-443.
            if self.my_error == BRepLibWireError::NonManifoldWire {
                // L442: VF = VL = TopoDS_Vertex() — nullify.
                self.vf = Shape::null();
                self.vl = Shape::null();
            }
        }
        // L445-449: test myShape is closed.
        if !self.vf.is_null() && !self.vl.is_null() && self.vf.is_same(&self.vl) {
            set_shape_closed(brep, &self.my_shape, true);
        }

        // L451-452.
        self.my_error = BRepLibWireError::WireDone;
        self.done();
    }

    /// OCCT BRepLib_MakeWire::Wire() (L457-460).
    pub fn wire(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT BRepLib_MakeWire::Edge() (L464-467).
    pub fn edge(&self) -> Shape {
        self.my_edge.clone()
    }

    /// OCCT BRepLib_MakeWire::Vertex() (L471-474).
    pub fn vertex(&self) -> Shape {
        self.my_vertex.clone()
    }

    /// OCCT BRepLib_MakeShape::IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT BRepLib_MakeWire::Error() (L485-488).
    pub fn error(&self) -> BRepLibWireError {
        self.my_error
    }

    /// OCCT BRepLib_MakeShape::Done().
    fn done(&mut self) {
        self.done = true;
    }

    /// OCCT BRepLib_MakeShape::NotDone().
    fn not_done(&mut self) {
        self.done = false;
    }
}

// =========================================================================
// OCCT TopExp::Vertices — the edge form (TopExp.cxx L214-247).
// =========================================================================
pub fn topexp_vertices_edge(brep: &BRep, e: &Shape, _cum_ori: bool) -> (Shape, Shape) {
    let mut vfirst = Shape::null();
    let mut vlast = Shape::null();
    // L222-223: boolean flags (see #27021).
    let mut is_first_defined = false;
    let mut is_last_defined = false;

    // L225: TopoDS_Iterator ite(E, CumOri) over the vertex children.
    for av in edge_vertices(brep, e) {
        if av.orientation == Orientation::Forward {
            vfirst = av.clone();
            is_first_defined = true;
        } else if av.orientation == Orientation::Reversed {
            vlast = av.clone();
            is_last_defined = true;
        }
    }

    if !is_first_defined {
        vfirst = Shape::null();
    }
    if !is_last_defined {
        vlast = Shape::null();
    }
    (vfirst, vlast)
}

// =========================================================================
// OCCT TopExp::Vertices — the wire form (TopExp.cxx L255-299).
// =========================================================================
pub fn topexp_vertices_wire(brep: &BRep, w: &Shape) -> (Shape, Shape) {
    // L256: Vfirst = Vlast = TopoDS_Vertex() — nullify.
    let mut vfirst = Shape::null();
    let mut vlast = Shape::null();

    // L258: NCollection_Map vmap.
    let mut vmap: Vec<Shape> = Vec::new();
    // L261: TopoDS_Vertex V1, V2.
    let mut v2 = Shape::null();

    // L263-280.
    for es in wire_edges(brep, w) {
        let (v1, v1v2) = if es.orientation == Orientation::Reversed {
            // L266-268: TopExp::Vertices(E, V2, V1).
            let (last, first) = topexp_vertices_edge(brep, &es, false);
            (first, last)
        } else {
            // L270-272: TopExp::Vertices(E, V1, V2).
            topexp_vertices_edge(brep, &es, false)
        };
        v2 = v1v2;
        // L282-287: add or remove in the vertex map.
        let mut a = v1;
        a.orientation = Orientation::Forward;
        if let Some(pos) = vmap.iter().position(|x| x.is_equal(&a)) {
            vmap.remove(pos);
        } else {
            vmap.push(a);
        }
        let mut b = v2.clone();
        b.orientation = Orientation::Reversed;
        if let Some(pos) = vmap.iter().position(|x| x.is_equal(&b)) {
            vmap.remove(pos);
        } else {
            vmap.push(b);
        }
    }
    // L281-288: if the map is empty the wire is closed.
    if vmap.is_empty() {
        let mut f = v2.clone();
        f.orientation = Orientation::Forward;
        vfirst = f;
        let mut l = v2.clone();
        l.orientation = Orientation::Reversed;
        vlast = l;
    } else if vmap.len() == 2 {
        // L289-299: open.
        if let Some(k) = vmap
            .iter()
            .find(|k| k.orientation == Orientation::Forward)
            .cloned()
        {
            vfirst = k;
        }
        if let Some(k) = vmap
            .iter()
            .find(|k| k.orientation == Orientation::Reversed)
            .cloned()
        {
            vlast = k;
        }
    }
    (vfirst, vlast)
}

// =========================================================================
// Kernel accessors (the TopoDS_Iterator encodings).
// =========================================================================

/// The wire child list (TWireData::edges).
pub(crate) fn wire_edges(brep: &BRep, w: &Shape) -> Vec<Shape> {
    if w.shape_type() == ShapeType::Wire {
        brep.wire(w.clone()).edges.clone()
    } else {
        Vec::new()
    }
}

/// The edge vertex children (TEdgeData::first / ::last — the TopoDS_Iterator
/// over an edge), skipping the null slots.
fn edge_vertices(brep: &BRep, e: &Shape) -> Vec<Shape> {
    if e.shape_type() != ShapeType::Edge {
        return Vec::new();
    }
    let ed = brep.edge(e.clone());
    let mut out = Vec::new();
    if !ed.first.is_null() {
        out.push(ed.first.clone());
    }
    if !ed.last.is_null() {
        out.push(ed.last.clone());
    }
    out
}

/// The TShape CLOSED flag of a wire/edge.
pub(crate) fn shape_closed(brep: &BRep, s: &Shape) -> bool {
    match s.shape_type() {
        ShapeType::Wire => brep.wire(s.clone()).flags & tshape_flags::CLOSED != 0,
        ShapeType::Edge => brep.edge(s.clone()).flags & tshape_flags::CLOSED != 0,
        _ => false,
    }
}

/// Set the TShape CLOSED flag of a wire/edge.
pub(crate) fn set_shape_closed(brep: &mut BRep, s: &Shape, closed: bool) {
    match s.shape_type() {
        ShapeType::Wire => {
            let wd = brep.wire_mut(s.clone());
            if closed {
                wd.flags |= tshape_flags::CLOSED;
            } else {
                wd.flags &= !tshape_flags::CLOSED;
            }
        }
        ShapeType::Edge => {
            let ed = brep.edge_mut_inplace(s.clone());
            if closed {
                ed.flags |= tshape_flags::CLOSED;
            } else {
                ed.flags &= !tshape_flags::CLOSED;
            }
        }
        _ => {}
    }
}

/// OCCT BRep_Builder::Add(anEdge, aVertex) — the kernel edge encodes the two
/// ends in fixed slots, so the add fills `first` when null and `last`
/// otherwise (architecture adaptation, see the module header).
pub(crate) fn builder_add_edge_vertex(brep: &mut BRep, edge: &Shape, v: &Shape) {
    let ed = brep.edge_mut_inplace(edge.clone());
    if ed.first.is_null() {
        ed.first = v.clone();
    } else {
        ed.last = v.clone();
    }
}
