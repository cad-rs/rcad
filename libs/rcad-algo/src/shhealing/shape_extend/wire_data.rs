//! OCCT ShapeExtend_WireData (TKShHealing ShapeExtend package).
//!
//! 1:1 translation of `ShapeExtend_WireData.hxx` + `.cxx` — a wire modeled as
//! an ordered list of edge `Shape`s, allowing work with incorrect wires
//! (the data structure every `ShapeFix_Wire` operation operates on).
//!
//! TopoDS correspondence: an edge entry is a `rcad_kernel::topo::topods::Shape`
//! carrying its own Orientation; `IsSame` is the TShape-pointer identity;
//! `Wire()` rebuilds a wire in the associated `BRep` pool.

use rcad_kernel::topo::topods::{BRep, Orientation, Shape, TEdgeData, TShape};
use std::collections::HashSet;

/// OCCT ShapeExtend_WireData.
#[derive(Debug, Clone, Default)]
pub struct WireData {
    my_edges: Vec<Shape>,
    my_nonmanifold_edges: Vec<Shape>,
    /// Lazy seam rank list (`mySeams`), valid when `my_seam_f >= 0`.
    my_seams: Vec<i32>,
    my_seams_cache: HashSet<i32>,
    /// -1 = seams not computed, 0 = no seams.
    my_seam_f: i32,
    my_seam_r: i32,
    my_manifold_mode: bool,
}

/// OCCT TopoDS_Iterator over the wire's edges (rcad: the wire TShape edge
/// list, with orientations).
fn wire_edges(wire: &Shape) -> Vec<Shape> {
    match wire.data.as_ref() {
        TShape::Wire(wd) => wd.edges.clone(),
        _ => Vec::new(),
    }
}

/// OCCT TopoDS_Iterator over the edge's vertices: the FORWARD-oriented one is
/// V1, the REVERSED one is V2.
fn edge_vertices(edge: &Shape) -> (Option<Shape>, Option<Shape>) {
    let mut v1 = None;
    let mut v2 = None;
    if let TShape::Edge(ed) = edge.data.as_ref() {
        for sv in [&ed.first, &ed.last] {
            match sv.orientation {
                Orientation::Forward => v1 = Some(sv.clone()),
                Orientation::Reversed => v2 = Some(sv.clone()),
                _ => {}
            }
        }
    }
    (v1, v2)
}

/// OCCT BRep_Tool::Degenerated(edge).
fn edge_degenerated(edge: &Shape) -> bool {
    matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.degenerated)
}

impl WireData {
    /// OCCT ShapeExtend_WireData() — empty wire, no edges.
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

    /// OCCT Init(handle(WireData) other) — copies data from another WireData.
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

    /// OCCT Init(wire, chained, theManifold) — loads an already existing wire
    /// (L74-148).  The `chained == false` fallback re-explores the wire with
    /// BRepTools_WireExplorer; rcad keeps the stored (connected) edge order.
    pub fn init(&mut self, wire: &Shape, chained: bool, the_manifold: bool) -> bool {
        self.clear();
        self.my_manifold_mode = the_manifold;
        let mut ok = true;
        let mut vlast: Option<Shape> = None;
        for e in wire_edges(wire) {
            // protect against INTERNAL/EXTERNAL edges.
            if e.orientation != Orientation::Reversed && e.orientation != Orientation::Forward {
                self.my_nonmanifold_edges.push(e);
                continue;
            }

            let (v1, v2) = edge_vertices(&e);

            // chainage? Si pas bon et chained False on repart sur WireExplorer.
            if let (Some(vl), Some(v1)) = (&vlast, &v1) {
                if !vl.is_same(v1) && the_manifold {
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
            let nb = self.my_nonmanifold_edges.len();
            for i in 0..nb {
                self.my_edges.push(self.my_nonmanifold_edges[i].clone());
            }
            self.my_nonmanifold_edges.clear();
        }
        // refaire chainage? Par WireExplorer.
        if ok || chained {
            return ok;
        }

        self.clear();
        // rcad keeps the wire's stored edge order (BRepTools_WireExplorer
        // equivalent for connected wires).
        for e in wire_edges(wire) {
            self.my_edges.push(e);
        }
        ok
    }

    /// OCCT Clear() — clears data about the wire.
    pub fn clear(&mut self) {
        self.my_edges.clear();
        self.my_nonmanifold_edges.clear();
        self.my_seam_f = -1;
        self.my_seam_r = -1;
        self.my_seams.clear();
        self.my_seams_cache.clear();
        self.my_manifold_mode = true;
    }

    /// OCCT ComputeSeams(enforce) — computes the list of seam edges
    /// (L163-222).  A seam edge is present twice, once FORWARD and once
    /// REVERSED (same TShape).
    pub fn compute_seams(&mut self, enforce: bool) {
        if self.my_seam_f >= 0 && !enforce {
            return;
        }

        self.my_seams.clear();
        self.my_seam_f = 0;
        self.my_seam_r = 0;
        let nb = self.nb_edges();
        // OCCT: NCollection_IndexedMap ME of REVERSED edges (IsSame identity)
        // + SE[num] = rank of that edge. (ptr, rank) pairs.
        let mut me: Vec<(u64, i32)> = Vec::new();

        // deux passes : d'abord on mappe les Edges REVERSED.
        for i in 1..=nb {
            let s = self.edge(i);
            if s.orientation == Orientation::Reversed {
                // IndexedMap::Add: no duplicate entries.
                if !me.iter().any(|(ptr, _)| *ptr == s.ptr_id()) {
                    me.push((s.ptr_id(), i));
                }
            }
        }

        // ensuite on voit les Edges FORWARD qui y seraient deja.
        for i in 1..=nb {
            let s = self.edge(i);
            if s.orientation == Orientation::Reversed {
                continue;
            }
            let Some(pos) = me.iter().position(|(ptr, _)| *ptr == s.ptr_id()) else {
                continue;
            };
            let se_rank = me[pos].1;
            if self.my_seam_f == 0 {
                self.my_seam_f = i;
                self.my_seam_r = se_rank;
            } else {
                self.my_seams.push(i);
                self.my_seams.push(se_rank);
            }
        }

        self.my_seams_cache.clear();
        for v in &self.my_seams {
            self.my_seams_cache.insert(*v);
        }
    }

    /// OCCT SetLast(num) — circular permutation setting the `num`th edge last.
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

    /// OCCT SetDegeneratedLast() — sets the first degenerated edge (if any)
    /// as last one.
    pub fn set_degenerated_last(&mut self) {
        let nb = self.nb_edges();
        for i in 1..=nb {
            if edge_degenerated(&self.edge(i)) {
                self.set_last(i);
                return;
            }
        }
    }

    /// OCCT Add(edge, atnum) — adds an edge; `atnum == 0` appends at end,
    /// `atnum == 1` prepends, else inserts before `atnum`.
    pub fn add_edge(&mut self, edge: &Shape, atnum: i32) {
        if edge.orientation != Orientation::Reversed
            && edge.orientation != Orientation::Forward
            && self.my_manifold_mode
        {
            self.my_nonmanifold_edges.push(edge.clone());
            return;
        }
        if atnum == 0 {
            self.my_edges.push(edge.clone());
        } else {
            self.my_edges.insert((atnum - 1) as usize, edge.clone());
        }
        self.my_seam_f = -1;
    }

    /// OCCT Add(wire, atnum) — adds an entire wire (assumed ordered).
    pub fn add_wire(&mut self, wire: &Shape, atnum: i32) {
        let mut n = atnum;
        let mut nm_edges: Vec<Shape> = Vec::new();
        for edge in wire_edges(wire) {
            if edge.orientation != Orientation::Reversed && edge.orientation != Orientation::Forward
            {
                if self.my_manifold_mode {
                    self.my_nonmanifold_edges.push(edge);
                } else {
                    nm_edges.push(edge);
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
        for e in nm_edges {
            self.my_edges.push(e);
        }
        self.my_seam_f = -1;
    }

    /// OCCT Add(wire, atnum) — adds a wire in the form of WireData (L328-384).
    pub fn add_wire_data(&mut self, wire: &WireData, atnum: i32) {
        let mut nm_edges: Vec<Shape> = Vec::new();
        let mut n = atnum;
        for i in 1..=wire.nb_edges() {
            let a_e = wire.edge(i);
            if a_e.orientation == Orientation::Internal || a_e.orientation == Orientation::External
            {
                nm_edges.push(a_e);
                continue;
            }
            if n == 0 {
                self.my_edges.push(a_e);
            } else {
                self.my_edges.insert((n - 1) as usize, a_e);
                n += 1;
            }
        }

        // non-manifold edges for non-manifold wire should be added at end.
        for e in nm_edges {
            self.my_edges.push(e);
        }

        for i in 1..=wire.nb_nonmanifold_edges() {
            let e = wire.nonmanifold_edge(i);
            if self.my_manifold_mode {
                self.my_nonmanifold_edges.push(e);
            } else if n == 0 {
                self.my_edges.push(e);
            } else {
                self.my_edges.insert((n - 1) as usize, e);
                n += 1;
            }
        }

        self.my_seam_f = -1;
    }

    /// OCCT AddOriented(edge, mode): 0 at end direct, 1 at end reversed,
    /// 2 at start direct, 3 at start reversed, < 0 no adding.
    pub fn add_oriented_edge(&mut self, edge: &Shape, mode: i32) {
        let mut e = edge.clone();
        if mode == 1 || mode == 3 {
            e.orientation = match e.orientation {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
                other => other,
            };
        }
        self.add_edge(&e, mode / 2); // mode = 0,1 -> 0  mode = 2,3 -> 1
    }

    /// OCCT Remove(num) — removes an edge by rank (0 removes the last).
    pub fn remove(&mut self, num: i32) {
        let idx = if num > 0 { num - 1 } else { self.nb_edges() - 1 };
        self.my_edges.remove(idx as usize);
        self.my_seam_f = -1;
    }

    /// OCCT Set(edge, num) — replaces the edge at rank `num` (0 = last).
    pub fn set_edge(&mut self, edge: &Shape, num: i32) {
        if edge.orientation != Orientation::Reversed
            && edge.orientation != Orientation::Forward
            && self.my_manifold_mode
        {
            if num > 0 && (num as usize) <= self.my_nonmanifold_edges.len() {
                self.my_nonmanifold_edges[(num - 1) as usize] = edge.clone();
            } else {
                self.my_nonmanifold_edges.push(edge.clone());
            }
        } else {
            let idx = if num > 0 { num - 1 } else { self.nb_edges() as i32 - 1 };
            self.my_edges[idx as usize] = edge.clone();
        }
        self.my_seam_f = -1;
    }

    /// OCCT Reverse() — reverses the list order and each edge's orientation.
    pub fn reverse(&mut self) {
        let nb = self.nb_edges();
        for i in 1..=nb / 2 {
            let mut s1 = self.my_edges[(i - 1) as usize].clone();
            reverse_orientation(&mut s1);
            let mut s2 = self.my_edges[((nb + 1 - i) - 1) as usize].clone();
            reverse_orientation(&mut s2);
            self.my_edges[(i - 1) as usize] = s2;
            self.my_edges[((nb + 1 - i) - 1) as usize] = s1;
        }
        // nb d'edges impair : inverser aussi l'edge du milieu (rang inchange).
        if nb % 2 == 1 {
            let i = (nb + 1) / 2;
            let mut si = self.my_edges[(i - 1) as usize].clone();
            reverse_orientation(&mut si);
            self.my_edges[(i - 1) as usize] = si;
        }
        self.my_seam_f = -1;
    }

    /// OCCT NbEdges().
    pub fn nb_edges(&self) -> i32 {
        self.my_edges.len() as i32
    }

    /// OCCT Edge(num) — rank access; negative `num` returns the reversed
    /// edge.
    pub fn edge(&self, num: i32) -> Shape {
        if num < 0 {
            let mut e = self.edge(-num);
            e.orientation = match e.orientation {
                Orientation::Forward => Orientation::Reversed,
                Orientation::Reversed => Orientation::Forward,
                other => other,
            };
            return e;
        }
        self.my_edges[(num - 1) as usize].clone()
    }

    /// OCCT NbNonManifoldEdges().
    pub fn nb_nonmanifold_edges(&self) -> i32 {
        self.my_nonmanifold_edges.len() as i32
    }

    /// OCCT NonmanifoldEdge(num).
    pub fn nonmanifold_edge(&self, num: i32) -> Shape {
        self.my_nonmanifold_edges[(num - 1) as usize].clone()
    }

    /// OCCT Index(edge) — the rank of the edge (seam orientation checked);
    /// 0 if not found.
    pub fn index(&mut self, edge: &Shape) -> i32 {
        for i in 1..=self.nb_edges() {
            let e = self.edge(i);
            if e.is_same(edge) && (e.orientation == edge.orientation || !self.is_seam(i)) {
                return i;
            }
        }
        0
    }

    /// OCCT IsSeam(num) — the edge is a seam (present FORWARD and REVERSED).
    pub fn is_seam(&mut self, num: i32) -> bool {
        if self.my_seam_f < 0 {
            self.compute_seams(false);
        }
        if self.my_seam_f == 0 {
            return false;
        }
        if num == self.my_seam_f || num == self.my_seam_r {
            return true;
        }
        self.my_seams_cache.contains(&num)
    }

    /// OCCT Wire() — builds the wire in the given `BRep` pool from the
    /// current edges (OCCT: BRep_Builder::MakeWire/Add + Closed flag).
    pub fn wire(&self, brep: &mut BRep) -> Shape {
        let mut all = Vec::new();
        let nb = self.nb_edges();
        for i in 1..=nb {
            all.push(self.edge(i));
        }
        if self.my_manifold_mode {
            let nb = self.nb_nonmanifold_edges();
            for i in 1..=nb {
                all.push(self.nonmanifold_edge(i));
            }
        }
        brep.add_twire(all)
    }

    /// OCCT ManifoldMode().
    pub fn manifold_mode(&self) -> bool {
        self.my_manifold_mode
    }

    /// OCCT ManifoldMode() setter (returns a mutable reference in OCCT).
    pub fn set_manifold_mode(&mut self, mode: bool) {
        self.my_manifold_mode = mode;
    }

    /// Access to the raw edge list (equivalent of iterating `myEdges`).
    pub fn edges(&self) -> &[Shape] {
        &self.my_edges
    }

    /// Data summary for the ported TEdgeData (helper used by tests).
    pub fn edge_data(edge: &Shape) -> Option<&TEdgeData> {
        match edge.data.as_ref() {
            TShape::Edge(ed) => Some(ed),
            _ => None,
        }
    }
}

/// OCCT TopoDS_Shape::Reverse — flips Forward <-> Reversed.
fn reverse_orientation(s: &mut Shape) {
    s.orientation = match s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    };
}
