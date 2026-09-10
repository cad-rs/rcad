//! OCCT ShapeExtend_WireData (TKShHealing): `.hxx` L59-230 and `.cxx`
//! L40-725 — a wire modeled as an ordered list of edges, allowing work with
//! incorrect wires (the data structure every `ShapeFix_Wire` operation
//! operates on).
//!
//! TopoDS correspondence: an edge entry is a `rcad_kernel::topo::topods::Shape`
//! carrying its own Orientation; `Wire()`/`WireAPIMake()` rebuild a wire in
//! the associated `BRep` pool.
//!
//! Architecture mapping notes:
//! - OCCT `IsSame` (same TShape + same Location) maps to
//!   `Shape::is_partner` (ptr_id + location index; the location index is
//!   stable within one BRep pool), following the identity mapping documented
//!   in `shape_build/brep_tool.rs`.  The NCollection maps keyed by
//!   `TopTools_ShapeMapHasher` (IsSame identity) use `(ptr_id, location)`
//!   keys.
//! - OCCT `bool& ManifoldMode()` (a reference accessor) maps to the
//!   `manifold_mode` getter + `set_manifold_mode` setter pair.
//! - `handle(NCollection_HSequence(TopoDS_Shape))` maps to `Vec<Shape>`.

use crate::brep_algo::tool::top_exp_vertices_wire;
use crate::shhealing::shape_build::brep_tool::{
    builder_add, iter_subshapes, set_flag_inplace, shape_is_null,
};
use rcad_kernel::topods::{
    BRep, BRepBuilder, BRepTool, Orientation, Shape, ShapeType, TShape, tshape_flags,
};
use std::collections::HashSet;

/// OCCT ShapeExtend_WireData (ShapeExtend_WireData.hxx L59-230).
pub struct WireData {
    /// OCCT myEdges (hxx L223).
    my_edges: Vec<Shape>,
    /// OCCT myNonmanifoldEdges (hxx L224).
    my_nonmanifold_edges: Vec<Shape>,
    /// OCCT mySeams (hxx L225) — the lazy seam rank list, valid when
    /// `my_seam_f >= 0`.
    my_seams: Vec<i32>,
    /// OCCT mySeamsCache (hxx L226, TColStd_PackedMapOfInteger).
    my_seams_cache: HashSet<i32>,
    /// OCCT mySeamF (hxx L227): -1 = seams not computed, 0 = no seams.
    my_seam_f: i32,
    /// OCCT mySeamR (hxx L228).
    my_seam_r: i32,
    /// OCCT myManifoldMode (hxx L229).
    my_manifold_mode: bool,
}

/// OCCT `TopoDS_Iterator` (defaults cumOri = cumLoc = True) mapped to the
/// `brep_tool::iter_subshapes` primitive: children with the cumulative
/// orientation and location applied.
fn wire_edges(brep: &mut BRep, wire: &Shape) -> Vec<Shape> {
    iter_subshapes(brep, wire, true, true)
}

/// OCCT TopoDS_Iterator over the edge's vertices (cxx L92-103): the
/// FORWARD-oriented one is V1, the REVERSED one is V2 (orientations composed
/// with the edge's own, like the OCCT vertex iterator).
fn edge_vertices(brep: &mut BRep, edge: &Shape) -> (Option<Shape>, Option<Shape>) {
    let mut v1 = None;
    let mut v2 = None;
    for sv in iter_subshapes(brep, edge, true, true) {
        match sv.orientation {
            Orientation::Forward => v1 = Some(sv),
            Orientation::Reversed => v2 = Some(sv),
            _ => {}
        }
    }
    (v1, v2)
}

/// OCCT BRep_Tool::Degenerated(edge).
fn edge_degenerated(edge: &Shape) -> bool {
    matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.degenerated)
}

/// OCCT TopoDS_Shape::Reverse (flips Forward <-> Reversed; INTERNAL and
/// EXTERNAL pass through unchanged).
fn reverse_orientation(s: &mut Shape) {
    s.orientation = match s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    };
}

impl WireData {
    /// OCCT ShapeExtend_WireData() (cxx L40-43): empty constructor, creates
    /// empty wire with no edges (calls Clear()).
    pub fn new() -> Self {
        let mut wd = WireData {
            my_edges: Vec::new(),
            my_nonmanifold_edges: Vec::new(),
            my_seams: Vec::new(),
            my_seams_cache: HashSet::new(),
            my_seam_f: -1,
            my_seam_r: -1,
            my_manifold_mode: true,
        };
        wd.clear();
        wd
    }

    /// OCCT ShapeExtend_WireData(wire, chained, theManifold) (cxx L47-52):
    /// constructor initializing the data from TopoDS_Wire (calls Init; the
    /// bool result is discarded like in OCCT).
    pub fn new_from_wire(brep: &mut BRep, wire: &Shape, chained: bool, the_manifold: bool) -> Self {
        let mut wd = WireData::new();
        wd.init(brep, wire, chained, the_manifold);
        wd
    }

    /// OCCT Init(handle(ShapeExtend_WireData) other) (cxx L56-70): copies
    /// data from another WireData.  The OCCT handle-null check has no rcad
    /// counterpart (a Rust reference cannot be null).
    pub fn init_from_other(&mut self, other: &WireData) {
        self.clear();
        let nb = other.nb_edges();
        for i in 1..=nb {
            self.add_edge(&other.edge(i), 0);
        }
        let nb = other.nb_nonmanifold_edges();
        for i in 1..=nb {
            self.add_edge(&other.nonmanifold_edge(i), 0);
        }
        self.my_manifold_mode = other.manifold_mode();
    }

    /// OCCT Init(wire, chained, theManifold) (cxx L74-148): loads an already
    /// existing wire.  When `chained` is True the edges are added in the
    /// sequence as they are explored by TopoDS_Iterator; else the wire is
    /// explored by BRepTools_WireExplorer and it is guaranteed that edges
    /// will be sequentially connected.
    pub fn init(&mut self, brep: &mut BRep, wire: &Shape, chained: bool, the_manifold: bool) -> bool {
        self.clear();
        self.my_manifold_mode = the_manifold;
        let mut ok = true;
        let mut vlast: Option<Shape> = None;
        for e in wire_edges(brep, wire) {
            // protect against INTERNAL/EXTERNAL edges (cxx L84-89).
            if e.orientation != Orientation::Reversed && e.orientation != Orientation::Forward {
                self.my_nonmanifold_edges.push(e);
                continue;
            }

            let (v1, v2) = edge_vertices(brep, &e);

            // chainage? Si pas bon et chained False on repart sur WireExplorer
            // (cxx L106-113: !Vlast.IsNull() && !Vlast.IsSame(V1) &&
            // theManifold — IsSame against a null V1 is false).
            let vlast_is_null = match &vlast {
                Some(v) => v.is_null(),
                None => true,
            };
            if !vlast_is_null {
                let same = match (&vlast, &v1) {
                    (Some(vl), Some(v1)) => vl.is_partner(v1),
                    _ => false,
                };
                if !same && the_manifold {
                    ok = false;
                    if !chained {
                        break;
                    }
                }
            }
            vlast = v2;
            if wire.orientation == Orientation::Reversed {
                self.my_edges.insert(0, e);
            } else {
                self.my_edges.push(e);
            }
        }

        if !self.my_manifold_mode {
            let nb = self.my_nonmanifold_edges.len() as i32;
            for i in 1..=nb {
                self.my_edges.push(self.my_nonmanifold_edges[(i - 1) as usize].clone());
            }
            self.my_nonmanifold_edges.clear();
        }
        //    refaire chainage ?  Par WireExplorer (cxx L135-145).
        if ok || chained {
            return ok;
        }

        self.clear();
        let mut we = BRepToolsWireExplorer::new(brep, wire);
        while we.more() {
            self.my_edges.push(we.current());
            we.next();
        }

        ok
    }

    /// OCCT Clear() (cxx L152-159): clears data about the wire.
    pub fn clear(&mut self) {
        // OCCT L154-155: myEdges / myNonmanifoldEdges = new sequences.
        self.my_edges.clear();
        self.my_nonmanifold_edges.clear();
        self.my_seam_f = -1;
        self.my_seam_r = -1;
        // OCCT L157: mySeams.Nullify() — the seam rank list is dropped (not
        // cleared to an empty sequence); ComputeSeams re-creates it.
        self.my_seams.clear();
        // OCCT Clear() leaves mySeamsCache untouched (cxx L152-159 has no
        // cache clear); the mySeamF = -1 protocol above forces IsSeam to
        // recompute (and rebuild the cache) before it can be read again.
        self.my_manifold_mode = true;
    }

    /// OCCT ComputeSeams(enforce) (cxx L163-222): computes the list of seam
    /// edges.  By default (direct call) computing is enforced; for indirect
    /// call (from IsSeam) it is redone only if not yet already done or if
    /// the list of edges has changed.  A seam edge is present twice in the
    /// list, once as FORWARD and once as REVERSED.
    pub fn compute_seams(&mut self, enforce: bool) {
        if self.my_seam_f >= 0 && !enforce {
            return;
        }

        // OCCT L170-175: mySeams = new HSequence<int>; mySeamF = mySeamR = 0;
        // ME = IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> (IsSame
        // identity -> (ptr_id, location) key); SE[num] = the rank of the
        // REVERSED edge stored at map entry num.
        self.my_seams.clear();
        self.my_seam_f = 0;
        self.my_seam_r = 0;
        let nb = self.nb_edges();
        let mut me: Vec<(u64, u32)> = Vec::new();
        let mut se = vec![0i32; (nb + 1) as usize];

        //  deux passes : d abord on mappe les Edges REVERSED
        //  Pour chacune, on note aussi son RANG dans la liste
        for i in 1..=nb {
            let s = self.edge(i);
            if s.orientation == Orientation::Reversed {
                // OCCT IndexedMap::Add: appends only when not yet present.
                let key = (s.ptr_id(), s.location);
                let mut num = None;
                for (idx, k) in me.iter().enumerate() {
                    if *k == key {
                        num = Some(idx + 1);
                        break;
                    }
                }
                let num = match num {
                    Some(num) => num,
                    None => {
                        me.push(key);
                        me.len()
                    }
                };
                se[num] = i;
            }
        }

        //  ensuite on voit les Edges FORWARD qui y seraient deja -> on note
        //  leur n0 c-a-d le n0 de la directe ET de la reverse
        for i in 1..=nb {
            let s = self.edge(i);
            if s.orientation == Orientation::Reversed {
                continue;
            }
            let key = (s.ptr_id(), s.location);
            let num = match me.iter().position(|k| *k == key) {
                Some(pos) => pos + 1,
                None => continue, // OCCT: num <= 0 -> continue.
            };
            if self.my_seam_f == 0 {
                self.my_seam_f = i;
                self.my_seam_r = se[num];
            } else {
                self.my_seams.push(i);
                self.my_seams.push(se[num]);
            }
        }

        // OCCT L215-219: mySeamsCache.Clear() then re-filled from mySeams.
        self.my_seams_cache.clear();
        for v in &self.my_seams {
            self.my_seams_cache.insert(*v);
        }
    }

    /// OCCT SetLast(num) (cxx L226-240): does a circular permutation in order
    /// to set <num>th edge last.
    pub fn set_last(&mut self, num: i32) {
        if num == 0 {
            return;
        }
        let nb = self.nb_edges() as i32;
        let mut i = nb;
        while i > num {
            let edge = self.my_edges[(nb - 1) as usize].clone();
            self.my_edges.remove((nb - 1) as usize);
            self.my_edges.insert(0, edge);
            i -= 1;
        }
        self.my_seam_f = -1;
    }

    /// OCCT SetDegeneratedLast() (cxx L244-255): when the wire contains at
    /// least one degenerated edge, sets it as last one.
    pub fn set_degenerated_last(&mut self) {
        let nb = self.nb_edges();
        for i in 1..=nb {
            if edge_degenerated(&self.edge(i)) {
                self.set_last(i);
                return;
            }
        }
    }

    /// OCCT Add(edge, atnum) (cxx L259-281): adds an edge to a wire, being
    /// defined (not yet ended).  <atnum> = 0 (D) appends at end, 1 prepends
    /// at start, else inserts before <atnum>.  A Null Edge is simply ignored.
    pub fn add_edge(&mut self, edge: &Shape, atnum: i32) {
        if edge.orientation != Orientation::Reversed
            && edge.orientation != Orientation::Forward
            && self.my_manifold_mode
        {
            self.my_nonmanifold_edges.push(edge.clone());
            return;
        }

        // OCCT L268-271: a Null Edge is simply ignored.
        if shape_is_null(edge) {
            return;
        }
        if atnum == 0 {
            self.my_edges.push(edge.clone());
        } else {
            self.my_edges.insert((atnum - 1) as usize, edge.clone());
        }
        self.my_seam_f = -1;
    }

    /// OCCT Add(wire, atnum) (cxx L285-324): adds an entire wire, considered
    /// as a list of edges (the wire is assumed to be ordered).
    pub fn add_wire(&mut self, brep: &mut BRep, wire: &Shape, atnum: i32) {
        if wire.is_null() {
            return;
        }
        let mut n = atnum;
        let mut a_nm_edges: Vec<Shape> = Vec::new();
        for edge in wire_edges(brep, wire) {
            if edge.orientation != Orientation::Reversed && edge.orientation != Orientation::Forward
            {
                if self.my_manifold_mode {
                    self.my_nonmanifold_edges.push(edge);
                } else {
                    a_nm_edges.push(edge);
                }
                continue;
            }
            if n == 0 {
                self.my_edges.push(edge);
            } else {
                self.my_edges.insert((n - 1) as usize, edge);
                n += 1;
            }
        }
        let nb = a_nm_edges.len() as i32;
        for i in 1..=nb {
            self.my_edges.push(a_nm_edges[(i - 1) as usize].clone());
        }
        self.my_seam_f = -1;
    }

    /// OCCT Add(handle(ShapeExtend_WireData) wire, atnum) (cxx L328-384):
    /// adds a wire in the form of WireData.
    pub fn add_wire_data(&mut self, wire: &WireData, atnum: i32) {
        // OCCT L330-333: the handle-null check has no rcad counterpart.
        let mut a_nm_edges: Vec<Shape> = Vec::new();
        let mut n = atnum;
        for i in 1..=wire.nb_edges() {
            let a_e = wire.edge(i);
            if a_e.orientation == Orientation::Internal || a_e.orientation == Orientation::External
            {
                a_nm_edges.push(a_e);
                continue;
            }

            if n == 0 {
                self.my_edges.push(wire.edge(i));
            } else {
                self.my_edges.insert((n - 1) as usize, wire.edge(i));
                n += 1;
            }
        }

        // non-manifold edges for non-manifold wire should be added at end.
        for e in a_nm_edges {
            self.my_edges.push(e);
        }

        for i in 1..=wire.nb_nonmanifold_edges() {
            if self.my_manifold_mode {
                self.my_nonmanifold_edges.push(wire.nonmanifold_edge(i));
            } else {
                // OCCT L371-379 appends wire->Edge(i) here (the i-th edge of
                // the source wire data, NOT NonmanifoldEdge(i)); the 1:1
                // translation keeps the OCCT statement as written.
                if n == 0 {
                    self.my_edges.push(wire.edge(i));
                } else {
                    self.my_edges.insert((n - 1) as usize, wire.edge(i));
                    n += 1;
                }
            }
        }

        self.my_seam_f = -1;
    }

    /// OCCT Add(shape, atnum) (cxx L388-398): adds an edge or a wire invoking
    /// the corresponding method Add.
    pub fn add_shape(&mut self, brep: &mut BRep, shape: &Shape, atnum: i32) {
        if shape.shape_type() == ShapeType::Edge {
            self.add_edge(shape, atnum);
        } else if shape.shape_type() == ShapeType::Wire {
            self.add_wire(brep, shape, atnum);
        }
    }

    /// OCCT AddOriented(edge, mode) (cxx L402-414): adds an edge to start or
    /// end, according to <mode>: 0 at end as direct, 1 at end as reversed,
    /// 2 at start as direct, 3 at start as reversed, < 0 no adding.
    pub fn add_oriented_edge(&mut self, edge: &Shape, mode: i32) {
        if shape_is_null(edge) || mode < 0 {
            return;
        }
        let mut e = edge.clone();
        if mode == 1 || mode == 3 {
            reverse_orientation(&mut e);
        }
        self.add_edge(&e, mode / 2); // mode = 0,1 -> 0  mode = 2,3 -> 1
    }

    /// OCCT AddOriented(wire, mode) (cxx L418-430): adds a wire to start or
    /// end, according to <mode>.
    pub fn add_oriented_wire(&mut self, brep: &mut BRep, wire: &Shape, mode: i32) {
        if wire.is_null() || mode < 0 {
            return;
        }
        let mut w = wire.clone();
        if mode == 1 || mode == 3 {
            reverse_orientation(&mut w);
        }
        self.add_wire(brep, &w, mode / 2); // mode = 0,1 -> 0  mode = 2,3 -> 1
    }

    /// OCCT AddOriented(shape, mode) (cxx L432-442): adds an edge or a wire
    /// invoking the corresponding method AddOriented.
    pub fn add_oriented_shape(&mut self, brep: &mut BRep, shape: &Shape, mode: i32) {
        if shape.shape_type() == ShapeType::Edge {
            self.add_oriented_edge(shape, mode);
        } else if shape.shape_type() == ShapeType::Wire {
            self.add_oriented_wire(brep, shape, mode);
        }
    }

    /// OCCT Remove(num) (cxx L446-452): removes an Edge, given its rank (by
    /// default removes the last edge).
    pub fn remove(&mut self, num: i32) {
        let idx = if num > 0 { num } else { self.nb_edges() };
        self.my_edges.remove((idx - 1) as usize);

        self.my_seam_f = -1;
    }

    /// OCCT Set(edge, num) (cxx L456-476): replaces an edge at the given rank
    /// number <num> with a new one (default is last edge, <num> = 0).
    pub fn set_edge(&mut self, edge: &Shape, num: i32) {
        if edge.orientation != Orientation::Reversed
            && edge.orientation != Orientation::Forward
            && self.my_manifold_mode
        {
            if num <= self.my_nonmanifold_edges.len() as i32 {
                // OCCT SetValue(num, edge): num < 1 raises
                // Standard_OutOfRange; the rcad index panics the same way
                // (num == 0 underflows the index).
                self.my_nonmanifold_edges[(num - 1) as usize] = edge.clone();
            } else {
                self.my_nonmanifold_edges.push(edge.clone());
            }
        } else {
            let idx = if num > 0 { num } else { self.nb_edges() };
            self.my_edges[(idx - 1) as usize] = edge.clone();
        }
        self.my_seam_f = -1;
    }

    /// OCCT Reverse() (cxx L483-506): reverses the sense of the list and the
    /// orientation of each Edge.  Should be called when either the wire has
    /// no seam edges or the face is not available.
    pub fn reverse(&mut self) {
        let nb = self.nb_edges();

        // inverser les edges + les permuter pour inverser le wire
        for i in 1..=nb / 2 {
            let mut s1 = self.my_edges[(i - 1) as usize].clone();
            reverse_orientation(&mut s1);
            let mut s2 = self.my_edges[((nb + 1 - i) - 1) as usize].clone();
            reverse_orientation(&mut s2);
            self.my_edges[(i - 1) as usize] = s2;
            self.my_edges[((nb + 1 - i) - 1) as usize] = s1;
        }
        //  nb d edges impair : inverser aussi l edge du milieu (rang inchange)
        if nb % 2 == 1 {
            //  test impair
            let i = (nb + 1) / 2;
            let mut si = self.my_edges[(i - 1) as usize].clone();
            reverse_orientation(&mut si);
            self.my_edges[(i - 1) as usize] = si;
        }
        self.my_seam_f = -1;
    }

    /// OCCT Reverse(face) (cxx L545-572): reverses the sense of the list and
    /// the orientation of each Edge; the face is necessary for swapping
    /// pcurves of seam edges (when the edge is reversed, the pcurves must be
    /// swapped).  A null face performs no swapping.
    pub fn reverse_face(&mut self, brep: &mut BRep, face: &Shape) {
        self.reverse();
        if face.is_null() {
            return;
        }

        //  ATTENTION aux coutures
        //  Une edge de couture est presente deux fois, FWD et REV
        //  Les inverser revient a permuter leur role ... donc ne rien faire
        //  Il faut donc aussi permuter leurs pcurves
        self.compute_seams(true);
        if self.my_seam_f > 0 {
            swap_seam(brep, &self.my_edges[(self.my_seam_f - 1) as usize], face);
        }
        if self.my_seam_r > 0 {
            swap_seam(brep, &self.my_edges[(self.my_seam_r - 1) as usize], face);
        }
        let nb = self.my_seams.len() as i32;
        for i in 1..=nb {
            swap_seam(brep, &self.my_seams_value(i), face);
        }
        self.my_seam_f = -1;
    }

    /// OCCT mySeams->Value(i) (the seam rank list element).
    fn my_seams_value(&self, i: i32) -> Shape {
        self.my_edges[(self.my_seams[(i - 1) as usize] - 1) as usize].clone()
    }

    /// OCCT NbEdges() (cxx L576-579): returns the count of currently
    /// recorded edges.
    pub fn nb_edges(&self) -> i32 {
        self.my_edges.len() as i32
    }

    /// OCCT Edge(num) (cxx L583-592): returns <num>th Edge; a negative rank
    /// returns the reversed edge.
    pub fn edge(&self, num: i32) -> Shape {
        if num < 0 {
            let mut e = self.edge(-num);
            reverse_orientation(&mut e);
            return e;
        }
        self.my_edges[(num - 1) as usize].clone()
    }

    /// OCCT NbNonManifoldEdges() (cxx L596-599).
    pub fn nb_nonmanifold_edges(&self) -> i32 {
        self.my_nonmanifold_edges.len() as i32
    }

    /// OCCT NonmanifoldEdge(num) (cxx L603-612): returns <num>th nonmanifold
    /// Edge (a negative rank returns a Null Edge).
    pub fn nonmanifold_edge(&self, num: i32) -> Shape {
        if num < 0 {
            return Shape::null();
        }

        self.my_nonmanifold_edges[(num - 1) as usize].clone()
    }

    /// OCCT Index(edge) (cxx L616-626): returns the index of the edge; if the
    /// edge is a seam the orientation is also checked; returns 0 if the edge
    /// is not found in the list.
    pub fn index(&mut self, edge: &Shape) -> i32 {
        for i in 1..=self.nb_edges() {
            let e = self.edge(i);
            if e.is_partner(edge) && (e.orientation == edge.orientation || !self.is_seam(i)) {
                return i;
            }
        }
        0
    }

    /// OCCT IsSeam(num) (cxx L630-647): tells if an Edge is seam (see
    /// ComputeSeams).  An edge is considered as seam if it presents twice in
    /// the edge list, once as FORWARD and once as REVERSED.
    pub fn is_seam(&mut self, num: i32) -> bool {
        if self.my_seam_f < 0 {
            self.compute_seams(true);
        }
        if self.my_seam_f == 0 {
            return false;
        }

        if num == self.my_seam_f || num == self.my_seam_r {
            return true;
        }
        // Use hash set for O(1) lookup instead of O(n) linear search
        self.my_seams_cache.contains(&num)
    }

    /// OCCT Wire() (cxx L651-685): makes TopoDS_Wire using BRep_Builder (just
    /// creates the TopoDS_Wire object and adds all edges into it).  Should be
    /// called when the wire is correct and adjacent edges share common
    /// vertices.
    pub fn wire(&self, brep: &mut BRep) -> Shape {
        // OCCT L653-655: TopoDS_Wire W; BRep_Builder B; B.MakeWire(W).
        let w = brep.add_twire(Vec::new());
        let nb = self.nb_edges();
        let mut ismanifold = true;
        for i in 1..=nb {
            let a_e = self.edge(i);
            if a_e.orientation != Orientation::Forward
                && a_e.orientation != Orientation::Reversed
            {
                ismanifold = false;
            }
            // OCCT L665: B.Add(W, aE).
            builder_add(brep, &w, &a_e);
        }
        let mut closed = false;
        if ismanifold {
            // OCCT L669-675: TopExp::Vertices(W, vf, vl); the CLOSED flag is
            // set when the first and last vertices are the same (vf.IsSame(vl)
            // -> rcad is_partner).
            let (vf, vl) = top_exp_vertices_wire(&w);
            closed = match (&vf, &vl) {
                (Some(vf), Some(vl)) => !vf.is_null() && !vl.is_null() && vf.is_partner(vl),
                _ => false,
            };
        }
        // OCCT MakeWire leaves the flag false and only sets it in the
        // ismanifold branch (W.Closed(true), cxx L673); add_twire
        // pre-computes the BRep_Tool::IsClosed rule, so the flag is
        // normalized to the OCCT Wire() value in both branches.
        set_flag_inplace(brep, &w, tshape_flags::CLOSED, closed);
        if self.my_manifold_mode {
            let nb = self.nb_nonmanifold_edges();
            for i in 1..=nb {
                // OCCT L681: B.Add(W, NonmanifoldEdge(i)).
                builder_add(brep, &w, &self.nonmanifold_edge(i));
            }
        }
        w
    }

    /// OCCT WireAPIMake() (cxx L689-711): makes TopoDS_Wire using
    /// BRepBuilderAPI_MakeWire, which merges geometrically coincident
    /// vertices and can disturb the correct order of edges in the wire.  If
    /// the builder fails, a null shape is returned.
    pub fn wire_api_make(&self) -> Shape {
        let mut mw = BRepBuilderAPIMakeWire::new();
        let nb = self.nb_edges();
        for i in 1..=nb {
            mw.add(&self.edge(i));
        }
        if self.my_manifold_mode {
            let nb = self.nb_nonmanifold_edges();
            for i in 1..=nb {
                mw.add(&self.nonmanifold_edge(i));
            }
        }
        let mut w = Shape::null();
        if mw.is_done() {
            w = mw.wire();
        }
        w
    }

    /// OCCT NonmanifoldEdges() (cxx L715-718): returns the sequence of
    /// non-manifold edges (not empty when the wire data was set in manifold
    /// mode but the initial wire has INTERNAL orientation or contains
    /// INTERNAL edges).
    pub fn nonmanifold_edges(&self) -> &[Shape] {
        &self.my_nonmanifold_edges
    }

    /// OCCT ManifoldMode() (cxx L722-725, `bool&` accessor): the getter of
    /// the mode defining manifold wire data or not.
    pub fn manifold_mode(&self) -> bool {
        self.my_manifold_mode
    }

    /// OCCT ManifoldMode() assignment (the `bool&` reference written by the
    /// caller).
    pub fn set_manifold_mode(&mut self, mode: bool) {
        self.my_manifold_mode = mode;
    }
}

// ---------------------------------------------------------------------------
// GAP carriers and static helpers (the OCCT statics of the .cxx file)
// ---------------------------------------------------------------------------

/// OCCT SwapSeam (cxx L511-543, static auxiliary for Reverse): swaps the
/// pcurves of a seam edge on the face.  Only performed once per edge (the
/// REVERSED-oriented occurrence returns early).  When a pcurve is missing
/// the swap is skipped (the `:q0` protection, cxx L534-537).
///
/// `B.UpdateEdge(E, c2dr, c2df, theface, 0.)` + `B.Range(E, theface, uff,
/// ulf)` map to `BRepBuilder::update_edge_pcurve_closed(E, c2dr, c2df, F,
/// uff, ulf, 0.)`: the rcad CurveOnClosedSurface carrier stores one shared
/// range field, so the UpdateEdge + Range pair lands in one call with the
/// same net state (FORWARD side carries c2dr over the range (uff, ulf)).
fn swap_seam(brep: &mut BRep, s: &Shape, f: &Shape) {
    let mut e = s.clone();
    if shape_is_null(&e) || f.is_null() {
        return;
    }
    if e.orientation == Orientation::Reversed {
        return; // ne le faire qu une fois !
    }

    let mut theface = f.clone();
    theface.orientation = Orientation::Forward;

    //: S4136  double Tol = BRep_Tool::Tolerance(theface);
    // d abord FWD puis REV
    let (c2df, uff, ulf) = match brep.curve_on_surface(&e, &theface) {
        Some(v) => v,
        None => return,
    };
    e.orientation = Orientation::Reversed;
    let (c2dr, _ufr, _ulr) = match brep.curve_on_surface(&e, &theface) {
        Some(v) => v,
        None => return, //: q0
    };
    // On permute
    e.orientation = Orientation::Forward;
    let mut b = BRepBuilder::new();
    b.update_edge_pcurve_closed(brep, e, c2dr, c2df, theface, uff, ulf, 0.); //: S4136: Tol
}

/// GAP carrier for OCCT BRepTools_WireExplorer (TKBRep) — the
/// `BRepTools_WireExplorer we(wire)` construction used by Init (cxx L142).
/// The rcad reduction walks the wire's stored edge order (the connected
/// order kept in `TWireData::edges`); OCCT documents that in this mode "not
/// all edges will be found" for disconnected wires and wires with seam
/// edges — the stored-order walk visits exactly the stored order.  GAP:
/// closes with the TKBRep BRepTools batch.
struct BRepToolsWireExplorer {
    my_edges: Vec<Shape>,
    my_index: usize,
}

impl BRepToolsWireExplorer {
    /// OCCT BRepTools_WireExplorer(W) (the wire-only constructor; the face
    /// argument defaults to a Null face).
    fn new(brep: &mut BRep, the_w: &Shape) -> Self {
        BRepToolsWireExplorer {
            my_edges: wire_edges(brep, the_w),
            my_index: 0,
        }
    }

    /// OCCT More().
    fn more(&self) -> bool {
        self.my_index < self.my_edges.len()
    }

    /// OCCT Current().
    fn current(&self) -> Shape {
        self.my_edges[self.my_index].clone()
    }

    /// OCCT Next().
    fn next(&mut self) {
        self.my_index += 1;
    }
}

/// GAP carrier for OCCT BRepBuilderAPI_MakeWire (TKTopAlgo) — Add stores the
/// edge list; IsDone() is false so the callers take the OCCT not-done exit
/// (WireAPIMake returns the null wire).  GAP: closes with the
/// BRepBuilderAPI batch (the same pending-classification as the
/// BRepLibMakeWire carrier in brep_algo/normal_projection.rs).
#[allow(dead_code)]
struct BRepBuilderAPIMakeWire {
    my_edges: Vec<Shape>,
}

#[allow(dead_code)]
impl BRepBuilderAPIMakeWire {
    /// OCCT BRepBuilderAPI_MakeWire().
    fn new() -> Self {
        BRepBuilderAPIMakeWire {
            my_edges: Vec::new(),
        }
    }

    /// OCCT BRepBuilderAPI_MakeWire::Add(edge).
    fn add(&mut self, edge: &Shape) {
        self.my_edges.push(edge.clone());
    }

    /// OCCT BRepBuilderAPI_MakeWire::IsDone() — pending (false).
    fn is_done(&self) -> bool {
        false
    }

    /// OCCT BRepBuilderAPI_MakeWire::Wire() — pending (null).
    fn wire(&self) -> Shape {
        Shape::null()
    }
}
