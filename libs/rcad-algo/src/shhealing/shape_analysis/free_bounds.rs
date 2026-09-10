//! OCCT TKShHealing — ShapeAnalysis package class: `ShapeAnalysis_FreeBounds`
//! (`ShapeAnalysis_FreeBounds.hxx` L17-224 + `ShapeAnalysis_FreeBounds.lxx`
//! L22-32 + `ShapeAnalysis_FreeBounds.cxx` L1-690).
//!
//! This class is intended to output free bounds of the shape: the wires
//! consisting of edges referenced by the faces of the shape only once.
//! It works on two distinct types of shapes: a compound of faces (the sewing
//! analyzer forecasts the free bounds) and a compound of shells
//! (`ShapeAnalysis_Shell` outputs the actual free bounds).  It also provides
//! static methods for advanced use: connecting edges/wires to wires,
//! extracting closed sub-wires out of wires, dispatching wires into
//! compounds for closed and open wires.
//!
//! Docket status: `[B-1:1]` (W1-5) — consumer: TKOffset
//! (docs/TKSHHEALING_PATH_B_DOCKET.md §3).
//!
//! Architecture bridges (numbered, referenced by the methods below):
//! 1. `BRep` pool argument - OCCT `BRep_Tool` / `BRep_Builder` / `TopExp`
//!    read and mutate the TShape graph through the shape document; the rcad
//!    equivalents resolve TShape data through `rcad_kernel::BRep`.  Every
//!    method takes `brep` as that stand-in (the shape_build/brep_tool.rs and
//!    shape_extend/wire_data.rs precedents).
//! 2. `occ::handle<NCollection_HSequence<TopoDS_Shape>>` /
//!    `NCollection_HArray1<TopoDS_Shape>` -> `Vec<Shape>`; the OCCT
//!    `IsNull()` handle test maps to `Option<&[Shape]>` where the OCCT
//!    argument may be null (DispatchWires) and to the empty check where the
//!    sequence is only iterated.
//! 3. `NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, ShapeMapHasher>` /
//!    `NCollection_Map<TopoDS_Shape, ...>` -> `HashMap<(u64, u32), Shape>` /
//!    `HashSet<(u64, u32)>` keyed by the `TopTools_ShapeMapHasher` IsSame
//!    identity (TShape pointer + location index; the wire_data.rs mapping
//!    note — the location index is stable within one BRep pool).
//! 4. OCCT `const handle<...>&` arguments mutated through the shared handle
//!    (`edges`, `iwires`, `sewd` aliased into `saw`, `arrwires` aliased into
//!    the BoxBndTree selector) map to `&mut Vec<Shape>` arguments or to
//!    explicit parameters on the carrier methods (Rust cannot express the
//!    handle aliasing; the call sites keep the OCCT statement order).
//! 5. The OCCT 8.0 `Standard_DEPRECATED` out-parameter overloads are kept
//!    with an `_out` suffix (the edge.rs overload-suffix bridge).
//! 6. `gp_Pnt` -> `DVec3`; `Bnd_Box` -> `rcad_kernel::math::bnd::BndBox`.
//!
//! GAP carriers (untranslated classes this file depends on; each keeps the
//! OCCT dependency, call anchors and failure path, closing with its batch):
//! - [`BRepBuilderAPISewing`] (TKTopAlgo BRepBuilderAPI, E3-H W3+ on demand)
//!
//! Retired carriers (Rule 4, replaced by the W2 1:1 real bodies in the same
//! commit): ShapeAnalysis_Shell -> [`shell::ShapeAnalysisShell`],
//! ShapeAnalysis_Wire -> [`wire::ShapeAnalysisWire`],
//! ShapeAnalysis_BoxBndTreeSelector + the NCollection_UBTree re-host ->
//! [`box_bnd_tree`], ShapeAnalysis::FindBounds ->
//! [`analysis::find_bounds`], ShapeBuild_Vertex ->
//! [`shape_build::vertex::ShapeBuildVertex`], ShapeBuild_Edge::
//! CopyReplaceVertices -> [`shape_build::edge::ShapeBuildEdge`].

use std::collections::{HashMap, HashSet};

use glam::DVec3;
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topods::{tshape_flags, BRep, Orientation, Shape, ShapeType, TShape};

use crate::brep_algo::tool::brep_tool_pnt;
use crate::shhealing::shape_analysis::analysis::find_bounds;
use crate::shhealing::shape_analysis::box_bnd_tree::{
    NCollectionUBTree, NCollectionUBTreeFiller, ShapeAnalysisBoxBndTreeSelector,
};
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::shell::ShapeAnalysisShell;
use crate::shhealing::shape_analysis::wire::ShapeAnalysisWire;
use crate::shhealing::shape_build::brep_tool::{
    builder_add, iter_subshapes, set_flag_inplace, topexp_explorer,
};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_build::vertex::ShapeBuildVertex;
use crate::shhealing::shape_extend::explorer::ShapeExtendExplorer;
use crate::shhealing::shape_extend::status::ShapeExtendStatus;
use crate::shhealing::shape_extend::wire_data::WireData;

/// OCCT `TopoDS_Shape::IsSame` (same TShape, same Location; Orientations may
/// differ) — the location-index identity used by the shape-keyed maps in
/// this file (the edge.rs / wire_data.rs precedents).
fn shape_is_same(a: &Shape, b: &Shape) -> bool {
    a.ptr_id() == b.ptr_id() && a.location == b.location
}

/// OCCT TopoDS_Shape::Reverse (flips Forward <-> Reversed; INTERNAL and
/// EXTERNAL pass through unchanged).
fn topods_reverse(s: &mut Shape) {
    s.orientation = match s.orientation {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    };
}

/// OCCT BRep_Tool::Pnt(vertex): the vertex point (the brep_algo tool
/// identity-location precedent).  The OCCT Standard_NoSuchObject on a null
/// vertex maps to a panic with the same contract.
fn brep_tool_pnt_checked(v: &Shape) -> DVec3 {
    brep_tool_pnt(v).expect("BRep_Tool::Pnt: vertex is null")
}

/// OCCT BRep_Tool::Degenerated(edge).
fn brep_tool_degenerated(edge: &Shape) -> bool {
    matches!(edge.data.as_ref(), TShape::Edge(ed) if ed.degenerated)
}

// ---------------------------------------------------------------------------
// GAP carriers (untranslated dependencies; see the module doc).
// ---------------------------------------------------------------------------

/// GAP carrier for OCCT `BRepBuilderAPI_Sewing` (TKTopAlgo/BRepBuilderAPI,
/// untranslated — docket E3-H: W3+ on demand; the `BRepBuilderAPI_Sewing`
/// construction anchor is FreeBounds.cxx L70).  The dependency, the Add /
/// Perform / NbFreeEdges / FreeEdge call form and the OCCT outcome of an
/// empty sewing analysis (no free edges) are kept.  GAP: closes with the
/// BRepBuilderAPI batch.
#[allow(dead_code)]
struct BRepBuilderAPISewing {
    tolerance: f64,
    #[allow(dead_code)]
    option1: bool,
    #[allow(dead_code)]
    option2: bool,
}

#[allow(dead_code)]
impl BRepBuilderAPISewing {
    /// OCCT BRepBuilderAPI_Sewing(tolerance, option1, option2, ...).
    fn new(tolerance: f64, option1: bool, option2: bool) -> Self {
        BRepBuilderAPISewing {
            tolerance,
            option1,
            option2,
        }
    }

    /// OCCT Add(shape).
    fn add(&mut self, _shape: &Shape) {}

    /// OCCT Perform() — pending: the sewing analysis is untranslated, the
    /// carrier outcome equals an analysis that found nothing.
    fn perform(&mut self) {}

    /// OCCT NbFreeEdges() — pending (0 keeps the OCCT no-free-edges outcome).
    fn nb_free_edges(&self) -> i32 {
        0
    }

    /// OCCT FreeEdge(index) — pending (a Null shape).
    fn free_edge(&self, _index: i32) -> Shape {
        Shape::null()
    }
}

// ---------------------------------------------------------------------------
// The class.
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_FreeBounds (ShapeAnalysis_FreeBounds.hxx L60-220).
pub struct ShapeAnalysisFreeBounds {
    /// OCCT myWires (hxx L214): compound of closed wires out of free edges.
    my_wires: Shape,
    /// OCCT myEdges (hxx L215): compound of open wires out of free edges.
    my_edges: Shape,
    /// OCCT myTolerance (hxx L216).
    my_tolerance: f64,
    /// OCCT myShared (hxx L217).
    my_shared: bool,
    /// OCCT mySplitClosed (hxx L218).
    my_split_closed: bool,
    /// OCCT mySplitOpen (hxx L219).
    my_split_open: bool,
}

impl ShapeAnalysisFreeBounds {
    /// OCCT ShapeAnalysis_FreeBounds() (cxx L57): empty constructor.
    pub fn new() -> Self {
        ShapeAnalysisFreeBounds {
            my_wires: Shape::null(),
            my_edges: Shape::null(),
            my_tolerance: 0.0,
            my_shared: false,
            my_split_closed: false,
            my_split_open: false,
        }
    }

    /// OCCT ShapeAnalysis_FreeBounds(shape, toler, splitclosed = false,
    /// splitopen = true) (cxx L61-96): builds FORECASTING free bounds of
    /// the shape (a compound of faces) with the sewing analyzer.
    pub fn new_forecast(
        brep: &mut BRep,
        shape: &Shape,
        toler: f64,
        splitclosed: bool,
        splitopen: bool,
    ) -> Self {
        let mut result = ShapeAnalysisFreeBounds {
            my_wires: Shape::null(),
            my_edges: Shape::null(),
            my_tolerance: toler,
            my_shared: false,
            my_split_closed: splitclosed,
            my_split_open: splitopen,
        };
        // OCCT L70: BRepBuilderAPI_Sewing Sew(toler, false, false).
        let mut sew = BRepBuilderAPISewing::new(toler, false, false);
        // OCCT L71-74: for (TopoDS_Iterator S(shape); S.More(); S.Next()) Sew.Add(S.Value()).
        for s in iter_subshapes(brep, shape, true, true) {
            sew.add(&s);
        }
        // OCCT L75: Sew.Perform().
        sew.perform();
        //
        // Extract free edges.
        //
        // OCCT L79: int nbedge = Sew.NbFreeEdges().
        let nbedge = sew.nb_free_edges();
        // OCCT L80: handle edges = new NCollection_HSequence<TopoDS_Shape>.
        let mut edges: Vec<Shape> = Vec::new();
        // OCCT L82-89: collect the non-degenerated free edges.
        for iedge in 1..=nbedge {
            let an_edge = sew.free_edge(iedge);
            if !brep_tool_degenerated(&an_edge) {
                edges.push(an_edge);
            }
        }
        //
        // Chainage.
        //
        // OCCT L93-95: wires = ConnectEdgesToWires(edges, toler, false);
        // DispatchWires(wires, myWires, myEdges); SplitWires().
        let wires = Self::connect_edges_to_wires(brep, &mut edges, toler, false);
        let mut my_wires = result.my_wires.clone();
        let mut my_edges = result.my_edges.clone();
        Self::dispatch_wires(brep, Some(&wires), &mut my_wires, &mut my_edges);
        result.my_wires = my_wires;
        result.my_edges = my_edges;
        result.split_wires_private(brep);
        result
    }

    /// OCCT ShapeAnalysis_FreeBounds(shape, splitclosed = false, splitopen =
    /// true, checkinternaledges = false) (cxx L100-131): builds ACTUAL free
    /// bounds of the shape (a compound of shells) with ShapeAnalysis_Shell.
    pub fn new_actual(
        brep: &mut BRep,
        shape: &Shape,
        splitclosed: bool,
        splitopen: bool,
        checkinternaledges: bool,
    ) -> Self {
        let mut result = ShapeAnalysisFreeBounds {
            my_wires: Shape::null(),
            my_edges: Shape::null(),
            my_tolerance: 0.0,
            my_shared: true,
            my_split_closed: splitclosed,
            my_split_open: splitopen,
        };
        // OCCT L109-115: aTmpShell = MakeShell; Add every explored face.
        let a_tmp_shell = brep.add_tshell(Vec::new());
        for a_exp_face in topexp_explorer(brep, shape, ShapeType::Face) {
            builder_add(brep, &a_tmp_shell, &a_exp_face);
        }

        // OCCT L117-118: ShapeAnalysis_Shell sas;
        // sas.CheckOrientedShells(aTmpShell, true, checkinternaledges).
        let mut sas = ShapeAnalysisShell::new();
        sas.check_oriented_shells(brep, &a_tmp_shell, true, checkinternaledges);

        // OCCT L120-130.
        if sas.has_free_edges() {
            let see = ShapeExtendExplorer;
            // OCCT L122-124: edges = see.SeqFromCompound(sas.FreeEdges(), false).
            let a_free_edges = sas.free_edges(brep);
            let mut edges = see.seq_from_compound(brep, &a_free_edges, false);

            // OCCT L126-129: wires = ConnectEdgesToWires(edges,
            // Precision::Confusion(), true); DispatchWires; SplitWires.
            let wires =
                Self::connect_edges_to_wires(brep, &mut edges, CONFUSION, true);
            let mut my_wires = result.my_wires.clone();
            let mut my_edges = result.my_edges.clone();
            Self::dispatch_wires(brep, Some(&wires), &mut my_wires, &mut my_edges);
            result.my_wires = my_wires;
            result.my_edges = my_edges;
            result.split_wires_private(brep);
        }
        result
    }

    /// OCCT GetClosedWires() (lxx L22-25): compound of closed wires out of
    /// free edges.
    pub fn get_closed_wires(&self) -> &Shape {
        &self.my_wires
    }

    /// OCCT GetOpenWires() (lxx L29-32): compound of open wires out of free
    /// edges.
    pub fn get_open_wires(&self) -> &Shape {
        &self.my_edges
    }

    /// OCCT ConnectEdgesToWires(edges, toler, shared) (cxx L135-164):
    /// builds a sequence of wires out of a sequence of not sorted edges,
    /// trying to build wires of maximum length.  Bridge #4: `edges` is
    /// mutated in place (the orientations may change when connecting).
    pub fn connect_edges_to_wires(
        brep: &mut BRep,
        edges: &mut Vec<Shape>,
        toler: f64,
        shared: bool,
    ) -> Vec<Shape> {
        // OCCT L140-141: iwires = new HSequence; BRep_Builder B.
        let mut iwires: Vec<Shape> = Vec::new();

        // OCCT L144-150: every edge is wrapped into its own wire.
        for i in 1..=edges.len() as i32 {
            let wire = brep.add_twire(Vec::new());
            builder_add(brep, &wire, &edges[(i - 1) as usize]);
            iwires.push(wire);
        }

        // OCCT L152-153.
        let wires = Self::connect_wires_to_wires(brep, &mut iwires, toler, shared);

        // OCCT L155-161: mirror the resulting wire orientations onto the
        // edges.
        for i in 1..=edges.len() as i32 {
            if iwires[(i - 1) as usize].orientation == Orientation::Reversed {
                topods_reverse(&mut edges[(i - 1) as usize]);
            }
        }

        wires
    }

    /// OCCT ConnectEdgesToWires(edges, toler, shared, wires) — the
    /// Standard_DEPRECATED out-parameter overload (cxx L168-175).
    pub fn connect_edges_to_wires_out(
        brep: &mut BRep,
        edges: &mut Vec<Shape>,
        toler: f64,
        shared: bool,
        wires: &mut Vec<Shape>,
    ) {
        *wires = Self::connect_edges_to_wires(brep, edges, toler, shared);
    }

    /// OCCT ConnectWiresToWires(iwires, toler, shared) (cxx L179-186):
    /// connects wires from the given sequence into longer wires with an
    /// empty vertices map.
    pub fn connect_wires_to_wires(
        brep: &mut BRep,
        iwires: &mut Vec<Shape>,
        toler: f64,
        shared: bool,
    ) -> Vec<Shape> {
        // OCCT L184-185: NCollection_DataMap<...> map; return Connect...(map).
        let mut map: HashMap<(u64, u32), Shape> = HashMap::new();
        Self::connect_wires_to_wires_vertices(brep, iwires, toler, shared, &mut map)
    }

    /// OCCT ConnectWiresToWires(iwires, toler, shared, owires) — the
    /// Standard_DEPRECATED out-parameter overload (cxx L190-197).
    pub fn connect_wires_to_wires_out(
        brep: &mut BRep,
        iwires: &mut Vec<Shape>,
        toler: f64,
        shared: bool,
        owires: &mut Vec<Shape>,
    ) {
        *owires = Self::connect_wires_to_wires(brep, iwires, toler, shared);
    }

    /// OCCT ConnectWiresToWires(iwires, toler, shared, vertices)
    /// (cxx L432-441): also fills the map of original to new connecting
    /// vertices.
    pub fn connect_wires_to_wires_vertices(
        brep: &mut BRep,
        iwires: &mut Vec<Shape>,
        toler: f64,
        shared: bool,
        vertices: &mut HashMap<(u64, u32), Shape>,
    ) -> Vec<Shape> {
        let mut owires: Vec<Shape> = Vec::new();
        connect_wires_to_wires_impl(brep, iwires, toler, shared, &mut owires, vertices);
        owires
    }

    /// OCCT ConnectWiresToWires(iwires, toler, shared, owires, vertices) —
    /// the Standard_DEPRECATED out-parameter overload (cxx L445-453).
    pub fn connect_wires_to_wires_out_vertices(
        brep: &mut BRep,
        iwires: &mut Vec<Shape>,
        toler: f64,
        shared: bool,
        owires: &mut Vec<Shape>,
        vertices: &mut HashMap<(u64, u32), Shape>,
    ) {
        connect_wires_to_wires_impl(brep, iwires, toler, shared, owires, vertices);
    }

    /// OCCT SplitWires(wires, toler, shared, closed, open) (cxx L588-605):
    /// extracts closed sub-wires out of `wires` into `closed`, the remaining
    /// open wires go to `open`.
    pub fn split_wires(
        brep: &mut BRep,
        wires: &[Shape],
        toler: f64,
        shared: bool,
        closed: &mut Vec<Shape>,
        open: &mut Vec<Shape>,
    ) {
        // OCCT L595-596: closed = new HSequence; open = new HSequence.
        closed.clear();
        open.clear();

        // OCCT L598-604.
        for i in 1..=wires.len() as i32 {
            let mut tmpclosed: Vec<Shape> = Vec::new();
            let mut tmpopen: Vec<Shape> = Vec::new();
            split_wire(brep, &wires[(i - 1) as usize], toler, shared, &mut tmpclosed, &mut tmpopen);
            // OCCT L602-603: Append(HSequence) appends every element.
            closed.extend(tmpclosed);
            open.extend(tmpopen);
        }
    }

    /// OCCT DispatchWires(wires, closed, open) (cxx L609-639): dispatches a
    /// sequence of wires into two compounds, `closed` for closed wires and
    /// `open` for open wires.  Bridge #2: the OCCT null handle argument maps
    /// to `Option<&[Shape]>`.
    pub fn dispatch_wires(
        brep: &mut BRep,
        wires: Option<&[Shape]>,
        closed: &mut Shape,
        open: &mut Shape,
    ) {
        // OCCT L615-622: null compounds are created.
        if closed.is_null() {
            *closed = brep.add_tcompound(Vec::new());
        }
        if open.is_null() {
            *open = brep.add_tcompound(Vec::new());
        }
        // OCCT L623-626: a null wires sequence returns after the creation.
        let wires = match wires {
            Some(w) => w,
            None => return,
        };

        // OCCT L628-638.
        for iw in 1..=wires.len() as i32 {
            let wire = &wires[(iw - 1) as usize];
            if brep.has_flag(wire.clone(), tshape_flags::CLOSED) {
                builder_add(brep, closed, wire);
            } else {
                builder_add(brep, open, wire);
            }
        }
    }

    /// OCCT private SplitWires() (cxx L648-690): splits the compounds of
    /// closed (myWires) and open (myEdges) wires according to mySplitClosed
    /// / mySplitOpen and rebuilds the compounds.
    fn split_wires_private(&mut self, brep: &mut BRep) {
        // OCCT L650-653: nothing to do when both splits are off.
        if !self.my_split_closed && !self.my_split_open {
            return;
        }

        let see = ShapeExtendExplorer;
        // OCCT L655-658.
        let closedwires = see.seq_from_compound(brep, &self.my_wires, false);
        let openwires = see.seq_from_compound(brep, &self.my_edges, false);

        let cw1: Vec<Shape>;
        let ow1: Vec<Shape>;
        let cw2: Vec<Shape>;
        let ow2: Vec<Shape>;
        // OCCT L660-668.
        if self.my_split_closed {
            let a_wires = closedwires;
            let mut c = Vec::new();
            let mut o = Vec::new();
            Self::split_wires(brep, &a_wires, self.my_tolerance, self.my_shared, &mut c, &mut o);
            cw1 = c;
            ow1 = o;
        } else {
            cw1 = closedwires;
            ow1 = Vec::new();
        }

        // OCCT L670-678.
        if self.my_split_open {
            let a_wires = openwires;
            let mut c = Vec::new();
            let mut o = Vec::new();
            Self::split_wires(brep, &a_wires, self.my_tolerance, self.my_shared, &mut c, &mut o);
            cw2 = c;
            ow2 = o;
        } else {
            cw2 = Vec::new();
            ow2 = openwires;
        }

        // OCCT L680-683: closedwires = cw1; closedwires->Append(cw2); ...
        let mut closedwires = cw1;
        closedwires.extend(cw2);
        let mut openwires = ow1;
        openwires.extend(ow2);

        // OCCT L686-689: szv#4:S4163:12Mar99 SGI warns.
        let comp_wires = see.compound_from_seq(brep, &closedwires);
        let comp_edges = see.compound_from_seq(brep, &openwires);
        self.my_wires = comp_wires;
        self.my_edges = comp_edges;
    }
}

impl Default for ShapeAnalysisFreeBounds {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT connectWiresToWiresImpl (FreeBounds.cxx L201-428, static): the
/// ConnectWiresToWires engine shared by the public overloads.
fn connect_wires_to_wires_impl(
    brep: &mut BRep,
    iwires: &mut Vec<Shape>,
    toler: f64,
    shared: bool,
    owires: &mut Vec<Shape>,
    vertices: &mut HashMap<(u64, u32), Shape>,
) {
    // OCCT L208-211: a null or empty input returns.
    if iwires.is_empty() {
        return;
    }
    // OCCT L212-218: arrwires = new HArray1(1, iwires->Length()); filled.
    let mut arrwires: Vec<Shape> = Vec::with_capacity(iwires.len());
    for i in 1..=iwires.len() as i32 {
        arrwires.push(iwires[(i - 1) as usize].clone());
    }
    // OCCT L219: owires = new HSequence.
    owires.clear();
    // OCCT L220: double tolerance = std::max(toler, Precision::Confusion()).
    let tolerance = toler.max(CONFUSION);

    // OCCT L222-223: sewd = new ShapeExtend_WireData(arrwires->Value(1))
    // (the single-argument constructor: chained = true, manifold = true).
    let mut sewd = WireData::new_from_wire(brep, &arrwires[0], true, true);

    // OCCT L225: bool isUsedManifoldMode = true.
    let mut is_used_manifold_mode = true;

    // OCCT L227-231: the non-manifold first wire re-loads the data.
    if sewd.nb_edges() < 1 && sewd.nb_nonmanifold_edges() > 0 {
        is_used_manifold_mode = false;
        sewd = WireData::new_from_wire(brep, &arrwires[0], true, is_used_manifold_mode);
    }

    // OCCT L233-235: saw = new ShapeAnalysis_Wire; saw->Load(sewd);
    // SetPrecision(tolerance).  The OCCT Load stores the sewd handle —
    // saw->WireData() IS sewd from here on; the rcad value model transfers
    // the ownership into saw and the sewd reads/writes below go through
    // the accessors (the wire_data_mut aliasing bridge).
    let mut saw = ShapeAnalysisWire::new();
    saw.load_wire_data(sewd);
    saw.set_precision(tolerance);

    // OCCT L237-240: the UBTree + filler + selector; LoadList(1).
    let a_bb_tree = NCollectionUBTree::new();
    let mut a_tree_filler = NCollectionUBTreeFiller::new(a_bb_tree);
    let mut a_sel = ShapeAnalysisBoxBndTreeSelector::new(shared);
    a_sel.load_list(1);

    // OCCT L242-254: fill the tree with the bounds boxes of the wires.
    for inb_w in 2..=arrwires.len() as i32 {
        let tr_w = &arrwires[(inb_w - 1) as usize];
        // OCCT L246-247: TopoDS_Vertex trV1, trV2; FindBounds(trW, trV1, trV2).
        let mut tr_v1: Option<Shape> = None;
        let mut tr_v2: Option<Shape> = None;
        find_bounds(brep, tr_w, &mut tr_v1, &mut tr_v2);
        // OCCT L248-252: trP1/trP2 = BRep_Tool::Pnt(trV1/trV2); aBox.Set(trP1)
        // (the OCCT Set(P) builds the single-point box — BndBox::from_point);
        // aBox.Add(trP2); aBox.SetGap(tolerance).
        let tr_p1 = brep_tool_pnt_checked(&tr_v1.unwrap_or_else(Shape::null));
        let tr_p2 = brep_tool_pnt_checked(&tr_v2.unwrap_or_else(Shape::null));
        let mut a_box = BndBox::from_point(tr_p1);
        a_box.add_point(tr_p2);
        a_box.set_gap(tolerance);
        // OCCT L253: aTreeFiller.Add(inbW, aBox).
        a_tree_filler.add(inb_w, a_box);
    }

    // OCCT L256: aTreeFiller.Fill().
    a_tree_filler.fill();

    // OCCT L259-260: int nsel; ShapeAnalysis_Edge sae.
    let sae = ShapeAnalysisEdge::new();
    let mut done = false;

    while !done {
        // OCCT L264-266: found / tail / direct / lwire.
        let mut found = false;
        let mut tail = false;
        let mut direct = false;
        let mut lwire = 0;
        // OCCT L266: aSel.SetStop().
        a_sel.set_stop();
        // OCCT L267-270: FVBox / LVBox / Vf / Vl (sewd == saw->WireData()).
        let vf = sae.first_vertex(brep, &saw.wire_data().unwrap().edge(1));
        let vl = sae.last_vertex(
            brep,
            &saw.wire_data()
                .unwrap()
                .edge(saw.wire_data().unwrap().nb_edges()),
        );

        // OCCT L272-278: pf / pl points; FVBox.Set(pf); FVBox.SetGap;
        // LVBox.Set(pl); LVBox.SetGap (the OCCT Set(P) single-point box).
        let pf = brep_tool_pnt_checked(&vf);
        let pl = brep_tool_pnt_checked(&vl);
        let mut fv_box = BndBox::from_point(pf);
        fv_box.set_gap(tolerance);
        let mut lv_box = BndBox::from_point(pl);
        lv_box.set_gap(tolerance);

        // OCCT L280: aSel.DefineBoxes(FVBox, LVBox).
        a_sel.define_boxes(&fv_box, &lv_box);

        // OCCT L282-290: the shared-vertex or point-based definition.
        if shared {
            a_sel.define_vertexes(&vf, &vl);
        } else {
            a_sel.define_pnt(pf, pl);
            a_sel.set_tolerance(tolerance);
        }

        // OCCT L292: nsel = aBBTree.Select(aSel).
        let nsel = a_tree_filler.select(brep, &mut a_sel, &arrwires);

        // OCCT L294-301.
        if nsel != 0 && !a_sel.last_check_status(ShapeExtendStatus::Fail) {
            found = true;
            lwire = a_sel.get_nb();
            tail = a_sel.last_check_status(ShapeExtendStatus::Done1)
                || a_sel.last_check_status(ShapeExtendStatus::Done2);
            direct = a_sel.last_check_status(ShapeExtendStatus::Done1)
                || a_sel.last_check_status(ShapeExtendStatus::Done3);
            a_sel.load_list(lwire);
        }

        if found {
            // OCCT L305-308: the non-direct wire is reversed in the array.
            if !direct {
                topods_reverse(&mut arrwires[(lwire - 1) as usize]);
            }

            // OCCT L310-311: acurwd = new ShapeExtend_WireData(..., true,
            // isUsedManifoldMode).
            let acurwd = WireData::new_from_wire(
                brep,
                &arrwires[(lwire - 1) as usize],
                true,
                is_used_manifold_mode,
            );
            // OCCT L312-315: an empty continuation restarts the search.
            if acurwd.nb_edges() == 0 {
                continue;
            }
            // OCCT L316: sewd->Add(acurwd, (tail ? 0 : 1)).
            saw.wire_data_mut()
                .unwrap()
                .add_wire_data(&acurwd, if tail { 0 } else { 1 });
        } else {
            // OCCT L320-359: the CheckConnected fallback walk (the OCCT
            // for-condition re-reads saw->NbEdges() every iteration; the
            // CheckConnected prec argument keeps the OCCT default 0.0).
            let mut i = 1;
            while i <= saw.nb_edges() {
                if saw.check_connected(brep, i, 0.0) {
                    // OCCT L324-327: n2 / n1 / E1 / E2.
                    let n2 = i;
                    let n1 = if n2 > 1 { n2 - 1 } else { saw.nb_edges() };
                    let e1 = saw.wire_data().unwrap().edge(n1);
                    let e2 = saw.wire_data().unwrap().edge(n2);

                    // OCCT L329-331: Vprev / Vfol.
                    let vprev = sae.last_vertex(brep, &e1);
                    let vfol = sae.first_vertex(brep, &e2);

                    // OCCT L333-341: V = Vprev or the combined vertex
                    // (OCCT L339-340: ShapeBuild_Vertex sbv; V =
                    // sbv.CombineVertex(Vprev, Vfol) — the default
                    // tolFactor = 1.0001).
                    let v = if saw.last_check_status(ShapeExtendStatus::Done1) {
                        vprev.clone()
                    } else {
                        let sbv = ShapeBuildVertex;
                        sbv.combine_vertex(brep, &vprev, &vfol, 1.0001)
                    };
                    // OCCT L342-343: vertices.Bind(Vprev, V); Bind(Vfol, V).
                    vertices.insert((vprev.ptr_id(), vprev.location), v.clone());
                    vertices.insert((vfol.ptr_id(), vfol.location), v.clone());

                    // OCCT L345-357: ShapeBuild_Edge sbe; the Set calls
                    // (the OCCT null TopoDS_Vertex() argument maps to the
                    // null shape).
                    let sbe = ShapeBuildEdge;
                    if saw.nb_edges() < 2 {
                        let e2r = sbe.copy_replace_vertices(brep, &e2, &v, &v);
                        saw.wire_data_mut().unwrap().set_edge(&e2r, n2);
                    } else {
                        let e2r = sbe.copy_replace_vertices(brep, &e2, &v, &Shape::null());
                        saw.wire_data_mut().unwrap().set_edge(&e2r, n2);
                        if !saw.last_check_status(ShapeExtendStatus::Done1) {
                            let e1r = sbe.copy_replace_vertices(brep, &e1, &Shape::null(), &v);
                            saw.wire_data_mut().unwrap().set_edge(&e1r, n1);
                        }
                    }
                }
                i += 1;
            }

            // OCCT L361: TopoDS_Wire wire = sewd->Wire() (the Closed flag
            // write below mutates the TShape in place, not the wrapper).
            let wire = saw.wire_data().unwrap().wire(brep);
            // OCCT L362-368: the closed-flag checks.
            if is_used_manifold_mode {
                if !saw.check_connected(brep, 1, 0.0)
                    && saw.last_check_status(ShapeExtendStatus::Ok)
                {
                    set_flag_inplace(brep, &wire, tshape_flags::CLOSED, true);
                }
            } else {
                // OCCT L371-394: the vmap walk (the doubled vertices test).
                let mut vmap: HashSet<(u64, u32)> = HashSet::new();
                for e in iter_subshapes(brep, &wire, true, true) {
                    // OCCT L377: TopoDS_Iterator ite(E, false, true).
                    for v in iter_subshapes(brep, &e, false, true) {
                        if v.orientation == Orientation::Forward
                            || v.orientation == Orientation::Reversed
                        {
                            // OCCT L383-386: !vmap.Add(V) -> vmap.Remove(V).
                            let key = (v.ptr_id(), v.location);
                            if !vmap.insert(key) {
                                vmap.remove(&key);
                            }
                        }
                    }
                }
                if vmap.is_empty() {
                    set_flag_inplace(brep, &wire, tshape_flags::CLOSED, true);
                }
            }

            // OCCT L396-398: owires->Append(wire); sewd->Clear();
            // ManifoldMode() = isUsedManifoldMode.
            owires.push(wire);
            saw.wire_data_mut().unwrap().clear();
            saw.wire_data_mut()
                .unwrap()
                .set_manifold_mode(is_used_manifold_mode);

            // OCCT L400-415: pick the next not-yet-connected wire.
            lwire = -1;
            for i in 1..=arrwires.len() as i32 {
                if !a_sel.cont_wire(i) {
                    lwire = i;
                    saw.wire_data_mut()
                        .unwrap()
                        .add_wire(brep, &arrwires[(lwire - 1) as usize], 0);
                    a_sel.load_list(lwire);

                    if saw.wire_data().unwrap().nb_edges() > 0 {
                        break;
                    }
                    saw.wire_data_mut().unwrap().clear();
                }
            }

            // OCCT L417-420.
            if lwire == -1 {
                done = true;
            }
        }
    }

    // OCCT L424-427: iwires->SetValue(i, arrwires->Value(i)).
    for i in 1..=iwires.len() as i32 {
        iwires[(i - 1) as usize] = arrwires[(i - 1) as usize].clone();
    }
}

/// OCCT SplitWire (FreeBounds.cxx L455-586, static): splits a wire into its
/// closed sub-wires and the remaining open wire chain.
fn split_wire(
    brep: &mut BRep,
    wire: &Shape,
    toler: f64,
    shared: bool,
    closed: &mut Vec<Shape>,
    open: &mut Vec<Shape>,
) {
    // OCCT L461-463: closed = new HSequence; open = new HSequence;
    // tolerance = std::max(toler, Precision::Confusion()).
    closed.clear();
    open.clear();
    let tolerance = toler.max(CONFUSION);

    // OCCT L465-466: BRep_Builder B; ShapeAnalysis_Edge sae.
    let sae = ShapeAnalysisEdge::new();

    // OCCT L468-469: sewd = new ShapeExtend_WireData(wire); nbedges.
    let sewd = WireData::new_from_wire(brep, wire, true, true);
    let nbedges = sewd.nb_edges();

    // OCCT L472: ces — list of indices of connected edges to build a wire.
    let mut ces: Vec<i32> = Vec::new();
    // OCCT L474-477: statuses — 0-free, 1-in CES, 2-already in wire,
    // 3-no closed wire can be produced starting at this edge.
    let mut statuses: Vec<i32> = vec![0; nbedges as usize];

    // building closed wires
    // OCCT L481-573.
    for i in 1..=nbedges {
        if statuses[(i - 1) as usize] == 0 {
            ces.push(i);
            statuses[(i - 1) as usize] = 1; // putting into CES
            let mut search_backward = true;

            loop {
                let mut found: bool;
                let mut lvertex: Shape;
                let mut lpoint: DVec3;

                // searching for connection in ces
                if search_backward {
                    search_backward = false;
                    found = false;
                    let edge = sewd.edge(*ces.last().unwrap());
                    lvertex = sae.last_vertex(brep, &edge);
                    lpoint = brep_tool_pnt_checked(&lvertex);
                    // OCCT L505-513: for (j = ces.Length(); j >= 1 && !found; j--).
                    let mut j = ces.len() as i32;
                    while j >= 1 && !found {
                        let fv = sae.first_vertex(brep, &sewd.edge(ces[(j - 1) as usize]));
                        let fp = brep_tool_pnt_checked(&fv);
                        if (shared && shape_is_same(&lvertex, &fv))
                            || (!shared && lpoint.distance(fp) <= tolerance)
                        {
                            found = true;
                        }
                        j -= 1;
                    }

                    if found {
                        // OCCT L517: j++ — because of decreasing last iteration.
                        j += 1;
                        // making closed wire
                        // OCCT L519-528.
                        let wire1 = brep.add_twire(Vec::new());
                        for cesindex in j..=*ces.last().unwrap() {
                            builder_add(brep, &wire1, &sewd.edge(ces[(cesindex - 1) as usize]));
                            statuses[(ces[(cesindex - 1) as usize] - 1) as usize] = 2;
                        }
                        set_flag_inplace(brep, &wire1, tshape_flags::CLOSED, true);
                        closed.push(wire1);
                        // OCCT L528: ces.Remove(j, ces.Length()).
                        ces.truncate((j - 1) as usize);
                        if ces.is_empty() {
                            break;
                        }
                    }
                }

                // searching for connection among free edges
                // OCCT L537-540.
                found = false;
                let edge = sewd.edge(*ces.last().unwrap());
                lvertex = sae.last_vertex(brep, &edge);
                lpoint = brep_tool_pnt_checked(&lvertex);
                // OCCT L542-553.
                let mut j = 1i32;
                while j <= nbedges && !found {
                    if statuses[(j - 1) as usize] == 0 {
                        let fv = sae.first_vertex(brep, &sewd.edge(j));
                        let fp = brep_tool_pnt_checked(&fv);
                        if (shared && shape_is_same(&lvertex, &fv))
                            || (!shared && lpoint.distance(fp) <= tolerance)
                        {
                            found = true;
                        }
                    }
                    j += 1;
                }

                if found {
                    // OCCT L557: j-- — because of last iteration (the OCCT
                    // for-increment ran before the exit test).
                    j -= 1;
                    ces.push(j);
                    statuses[(j - 1) as usize] = 1; // putting into CES
                    search_backward = true;
                    continue;
                }

                // no edges found - mark the branch as open (use status 3)
                // OCCT L565-570.
                let last = *ces.last().unwrap();
                statuses[(last - 1) as usize] = 3;
                ces.remove((ces.len() - 1) as usize);
                if ces.is_empty() {
                    break;
                }
            }
        }
    }

    // building open wires
    // OCCT L576-583: collect the edges not consumed by closed wires.
    let mut edges: Vec<Shape> = Vec::new();
    for i in 1..=nbedges {
        if statuses[(i - 1) as usize] != 2 {
            edges.push(sewd.edge(i));
        }
    }

    // OCCT L585.
    let w_open = ShapeAnalysisFreeBounds::connect_edges_to_wires(brep, &mut edges, toler, shared);
    *open = w_open;
}
