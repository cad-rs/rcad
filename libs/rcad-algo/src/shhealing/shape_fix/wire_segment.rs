//! 1:1 translation of OCCT `ShapeFix_WireSegment`
//! (`TKShHealing/ShapeFix/ShapeFix_WireSegment.hxx` L17-166 +
//! `ShapeFix_WireSegment.cxx` L1-310, docket row `ShapeFix_WireSegment`).
//!
//! Function-count equation (OCCT ShapeFix_WireSegment.cxx = rcad
//! wire_segment.rs): `ShapeFix_WireSegment()` / `ShapeFix_WireSegment(wire,
//! ori)` / `Clear` / `Load` / `WireData` / `Orientation` x2 / `FirstVertex` /
//! `LastVertex` / `IsClosed` / `NbEdges` / `Edge` / `SetEdge` / `AddEdge` x2 /
//! `SetPatchIndex` / `DefineIUMin` / `DefineIUMax` / `DefineIVMin` /
//! `DefineIVMax` / `GetPatchIndex` / `CheckPatchIndex` / `SetVertex` /
//! `GetVertex` / `IsVertex` — 25 OCCT functions = 25 rcad functions.
//!
//! Architecture bridges:
//! 1. `BRep` pool argument — `ShapeAnalysis_Edge::FirstVertex/LastVertex`
//!    read the TShape graph through `rcad_kernel::BRep` (bridge #1 of
//!    shape_analysis/edge.rs).
//! 2. `occ::handle<ShapeExtend_WireData>` — the owned
//!    `shape_extend::wire_data::WireData` value.
//! 3. `occ::handle<NCollection_HSequence<int>>` x4 — the owned
//!    `Vec<i32>` sequences (1-based OCCT index -> index - 1).
//! 4. The OCCT_DEBUG prints are compiled out.

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation};

use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_build::brep_tool::occt_is_same;
use crate::shhealing::shape_extend::wire_data::WireData;

/// OCCT `#define MININD -32000` (cxx L116).
pub const MININD: i32 = -32000;
/// OCCT `#define MAXIND 32000` (cxx L117).
pub const MAXIND: i32 = 32000;

/// OCCT ShapeFix_WireSegment (hxx L56-164): auxiliary data-storage class for
/// ComposeShell; represents a segment of a wire (or a whole wire) with the
/// associated patch-index data and the traversal orientation flag.
pub struct ShapeFixWireSegment {
    /// OCCT myWire (hxx L157).
    my_wire: WireData,
    /// OCCT myVertex (hxx L158).
    my_vertex: Shape,
    /// OCCT myOrient (hxx L159).
    my_orient: Orientation,
    /// OCCT myIUMin (hxx L160).
    my_iu_min: Vec<i32>,
    /// OCCT myIUMax (hxx L161).
    my_iu_max: Vec<i32>,
    /// OCCT myIVMin (hxx L162).
    my_iv_min: Vec<i32>,
    /// OCCT myIVMax (hxx L163).
    my_iv_max: Vec<i32>,
}

impl Default for ShapeFixWireSegment {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixWireSegment {
    /// OCCT ShapeFix_WireSegment::ShapeFix_WireSegment() (cxx L25-29):
    /// creates the empty segment.
    pub fn new() -> Self {
        let mut this = ShapeFixWireSegment {
            my_wire: WireData::new(),
            my_vertex: Shape::null(),
            my_orient: Orientation::Forward,
            my_iu_min: Vec::new(),
            my_iu_max: Vec::new(),
            my_iv_min: Vec::new(),
            my_iv_max: Vec::new(),
        };
        this.clear();
        this.my_orient = Orientation::Forward;
        this
    }

    /// OCCT ShapeFix_WireSegment::ShapeFix_WireSegment(wire, ori) (cxx
    /// L33-38): creates the segment and initializes it with the wire and the
    /// orientation (the hxx default ori = TopAbs_EXTERNAL is passed by the
    /// caller — Rust has no default arguments).
    pub fn new_with_orientation(wire: &WireData, ori: Orientation) -> Self {
        let mut this = Self::new();
        this.load(wire);
        this.my_orient = ori;
        this
    }

    /// OCCT ShapeFix_WireSegment::Clear (cxx L42-51): clears all fields.
    pub fn clear(&mut self) {
        self.my_wire = WireData::new();
        self.my_wire.set_manifold_mode(false);
        self.my_iu_min = Vec::new();
        self.my_iu_max = Vec::new();
        self.my_iv_min = Vec::new();
        self.my_iv_max = Vec::new();
        self.my_vertex = Shape::null();
    }

    /// OCCT ShapeFix_WireSegment::Load (cxx L55-64): loads the wire.
    pub fn load(&mut self, wire: &WireData) {
        // L57: (the commented-out `myWire = wire` stays untranslated).
        // L58-59.
        self.clear();
        self.my_wire.set_manifold_mode(wire.manifold_mode());
        // L60-63.
        for i in 1..=wire.nb_edges() {
            let e = wire.edge(i);
            self.add_edge(i, &e);
        }
    }

    /// OCCT ShapeFix_WireSegment::WireData (cxx L68-71): returns the wire.
    pub fn wire_data(&self) -> &WireData {
        &self.my_wire
    }

    /// OCCT ShapeFix_WireSegment::Orientation(ori) (cxx L75-78): sets the
    /// orientation flag.
    pub fn set_orientation(&mut self, ori: Orientation) {
        self.my_orient = ori;
    }

    /// OCCT ShapeFix_WireSegment::Orientation() (cxx L82-85): returns the
    /// orientation flag.
    pub fn orientation(&self) -> Orientation {
        self.my_orient
    }

    /// OCCT ShapeFix_WireSegment::FirstVertex (cxx L89-93): the first vertex
    /// of the first edge in the wire (no dependence on Orientation()).
    pub fn first_vertex(&self, brep: &mut BRep) -> Shape {
        let sae = ShapeAnalysisEdge::new();
        let e = self.my_wire.edge(1);
        sae.first_vertex(brep, &e)
    }

    /// OCCT ShapeFix_WireSegment::LastVertex (cxx L97-101): the last vertex
    /// of the last edge in the wire.
    pub fn last_vertex(&self, brep: &mut BRep) -> Shape {
        let sae = ShapeAnalysisEdge::new();
        let e = self.my_wire.edge(self.my_wire.nb_edges());
        sae.last_vertex(brep, &e)
    }

    /// OCCT ShapeFix_WireSegment::IsClosed (cxx L105-110): FirstVertex() ==
    /// LastVertex().
    pub fn is_closed(&self, brep: &mut BRep) -> bool {
        let v = self.first_vertex(brep);
        let last = self.last_vertex(brep);
        occt_is_same(brep, &v, &last)
    }

    /// OCCT ShapeFix_WireSegment::NbEdges (cxx L121-124): the number of
    /// edges in the wire.
    pub fn nb_edges(&self) -> i32 {
        self.my_wire.nb_edges()
    }

    /// OCCT ShapeFix_WireSegment::Edge (cxx L128-132): the edge by index.
    pub fn edge(&self, i: i32) -> Shape {
        self.my_wire.edge(i)
    }

    /// OCCT ShapeFix_WireSegment::SetEdge (cxx L136-139): replaces the edge
    /// at index i.
    pub fn set_edge(&mut self, i: i32, edge: &Shape) {
        self.my_wire.set_edge(edge, i);
    }

    /// OCCT ShapeFix_WireSegment::AddEdge(i, edge) (cxx L143-146): inserts a
    /// new edge with implicitly defined (indefinite) patch indices.
    pub fn add_edge(&mut self, i: i32, edge: &Shape) {
        // L145: AddEdge(i, edge, MININD, MAXIND, MININD, MAXIND).
        self.add_edge_with_patch_index(i, edge, MININD, MAXIND, MININD, MAXIND);
    }

    /// OCCT ShapeFix_WireSegment::AddEdge(i, edge, iumin, iumax, ivmin,
    /// ivmax) (cxx L150-172): inserts a new edge with explicitly defined
    /// patch indices; i == 0 appends at the end of the wire.
    pub fn add_edge_with_patch_index(
        &mut self,
        i: i32,
        edge: &Shape,
        iumin: i32,
        iumax: i32,
        ivmin: i32,
        ivmax: i32,
    ) {
        // L157.
        self.my_wire.add_edge(edge, i);
        // L158-171.
        if i == 0 {
            // L160-163: myIUMin->Append(iumin) etc.
            self.my_iu_min.push(iumin);
            self.my_iu_max.push(iumax);
            self.my_iv_min.push(ivmin);
            self.my_iv_max.push(ivmax);
        } else {
            // L165-170: myIUMin->InsertBefore(i, iumin) etc.
            self.my_iu_min.insert((i - 1) as usize, iumin);
            self.my_iu_max.insert((i - 1) as usize, iumax);
            self.my_iv_min.insert((i - 1) as usize, ivmin);
            self.my_iv_max.insert((i - 1) as usize, ivmax);
        }
    }

    /// OCCT ShapeFix_WireSegment::SetPatchIndex (cxx L176-186): sets the
    /// patch indices for edge i.
    pub fn set_patch_index(&mut self, i: i32, iumin: i32, iumax: i32, ivmin: i32, ivmax: i32) {
        // L182-185.
        self.my_iu_min[(i - 1) as usize] = iumin;
        self.my_iu_max[(i - 1) as usize] = iumax;
        self.my_iv_min[(i - 1) as usize] = ivmin;
        self.my_iv_max[(i - 1) as usize] = ivmax;
    }

    /// OCCT ShapeFix_WireSegment::DefineIUMin (cxx L190-200): modifies the
    /// minimal U patch index for edge i.
    pub fn define_iu_min(&mut self, i: i32, iumin: i32) {
        // L192-195.
        if self.my_iu_min[(i - 1) as usize] < iumin {
            self.my_iu_min[(i - 1) as usize] = iumin;
        }
        // L196-199: the OCCT_DEBUG print is compiled out.
    }

    /// OCCT ShapeFix_WireSegment::DefineIUMax (cxx L204-215): modifies the
    /// maximal U patch index for edge i.
    pub fn define_iu_max(&mut self, i: i32, iumax: i32) {
        // L206-209.
        if self.my_iu_max[(i - 1) as usize] > iumax {
            self.my_iu_max[(i - 1) as usize] = iumax;
        }
        // L210-214: the OCCT_DEBUG print is compiled out.
    }

    /// OCCT ShapeFix_WireSegment::DefineIVMin (cxx L219-230): modifies the
    /// minimal V patch index for edge i.
    pub fn define_iv_min(&mut self, i: i32, ivmin: i32) {
        // L221-224.
        if self.my_iv_min[(i - 1) as usize] < ivmin {
            self.my_iv_min[(i - 1) as usize] = ivmin;
        }
        // L225-229: the OCCT_DEBUG print is compiled out.
    }

    /// OCCT ShapeFix_WireSegment::DefineIVMax (cxx L234-245): modifies the
    /// maximal V patch index for edge i.
    pub fn define_iv_max(&mut self, i: i32, ivmax: i32) {
        // L236-239.
        if self.my_iv_max[(i - 1) as usize] > ivmax {
            self.my_iv_max[(i - 1) as usize] = ivmax;
        }
        // L240-244: the OCCT_DEBUG print is compiled out.
    }

    /// OCCT ShapeFix_WireSegment::GetPatchIndex (cxx L249-259): returns the
    /// patch indices for edge i.
    pub fn get_patch_index(&self, i: i32, iumin: &mut i32, iumax: &mut i32, ivmin: &mut i32, ivmax: &mut i32) {
        // L255-258.
        *iumin = self.my_iu_min[(i - 1) as usize];
        *iumax = self.my_iu_max[(i - 1) as usize];
        *ivmin = self.my_iv_min[(i - 1) as usize];
        *ivmax = self.my_iv_max[(i - 1) as usize];
    }

    /// OCCT ShapeFix_WireSegment::CheckPatchIndex (cxx L263-274): checks the
    /// patch indices for edge i to satisfy IUMin(i) <= IUMax(i) <=
    /// IUMin(i)+1.
    pub fn check_patch_index(&self, i: i32) -> bool {
        // L265-266.
        let d_u = self.my_iu_max[(i - 1) as usize] - self.my_iu_min[(i - 1) as usize];
        let d_v = self.my_iv_max[(i - 1) as usize] - self.my_iv_min[(i - 1) as usize];
        // L267-273.
        let ok = (d_u == 0 || d_u == 1) && (d_v == 0 || d_v == 1);
        ok
    }

    /// OCCT ShapeFix_WireSegment::SetVertex (cxx L278-282).
    pub fn set_vertex(&mut self, the_vertex: &Shape) {
        // L280-281.
        self.my_vertex = the_vertex.clone();
        // SetVertex(theVertex, MININD, MAXIND, MININD, MAXIND);
    }

    /// OCCT ShapeFix_WireSegment::GetVertex (cxx L300-303).
    pub fn get_vertex(&self) -> Shape {
        self.my_vertex.clone()
    }

    /// OCCT ShapeFix_WireSegment::IsVertex (cxx L307-310).
    pub fn is_vertex(&self) -> bool {
        !self.my_vertex.is_null()
    }
}
