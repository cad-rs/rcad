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
//! - [`ShapeAnalysisShell`] (ShapeAnalysis, W2)
//! - [`ShapeAnalysisWire`] (ShapeAnalysis, W2)
//! - [`ShapeAnalysisBoxBndTreeSelector`] + [`NCollectionUBTree`] (W2 +
//!   the TKFoundation NCollection_UBTree traversal)
//! - [`shape_analysis_find_bounds`] (ShapeAnalysis statics, W2) — bridged
//!   faithfully over the available TopExp/ShapeAnalysis_Edge primitives
//! - [`ShapeBuildVertex`] / [`shape_build_edge_copy_replace_vertices`]
//!   (ShapeBuild, W1-2) — bridged faithfully over the brep_algo tool
//!   primitives; W1-2 replaces the carriers with the 1:1 classes.

use std::collections::{HashMap, HashSet};

use glam::DVec3;
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topods::{tshape_flags, BRep, Orientation, Shape, ShapeType, TShape};

use crate::brep_algo::tool::{
    brep_tool_pnt, brep_tool_tolerance, builder_add_edge_vertex, builder_update_vertex_point_tol,
    empty_copied, top_exp_vertices_wire,
};
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_build::brep_tool::{
    builder_add, iter_subshapes, set_flag_inplace, shape_is_null, topexp_explorer,
};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_extend::explorer::ShapeExtendExplorer;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_extend::wire_data::WireData;

// OCCT Standard_Real.hxx RealSmall(): the smallest positive representable real.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;

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

/// GAP carrier for OCCT `ShapeAnalysis_Shell` (ShapeAnalysis, untranslated
/// W2 batch).  The dependency and the CheckOrientedShells / HasFreeEdges /
/// FreeEdges call form are kept; the carrier keeps the OCCT failure path
/// (`HasFreeEdges()` false — the FreeBounds constructor then leaves the
/// result compounds null).  GAP: closes with the W2 ShapeAnalysis batch.
#[allow(dead_code)]
struct ShapeAnalysisShell {
    #[allow(dead_code)]
    my_has_free_edges: bool,
}

#[allow(dead_code)]
impl ShapeAnalysisShell {
    /// OCCT ShapeAnalysis_Shell().
    fn new() -> Self {
        ShapeAnalysisShell {
            my_has_free_edges: false,
        }
    }

    /// OCCT CheckOrientedShells(shell, check_all, check_internal_edges)
    /// — pending: the shell check is untranslated; the carrier performs no
    /// analysis, so HasFreeEdges stays false (the OCCT not-checked outcome).
    fn check_oriented_shells(&mut self, _shell: &Shape, _check_all: bool, _check_internal_edges: bool) {
    }

    /// OCCT HasFreeEdges() — pending (false keeps the OCCT failure path).
    fn has_free_edges(&self) -> bool {
        false
    }

    /// OCCT FreeEdges() — pending (a Null compound).
    fn free_edges(&self) -> Shape {
        Shape::null()
    }
}

/// GAP carrier for OCCT `ShapeAnalysis_Wire` (ShapeAnalysis, untranslated
/// W2 batch) — the Load / SetPrecision / NbEdges / CheckConnected /
/// LastCheckStatus call form used by connectWiresToWiresImpl.  Architecture
/// bridge #4: OCCT `Load` stores a handle sharing the caller's WireData;
/// the rcad value model expresses the aliasing by passing the loaded wire
/// data to the query methods (NbEdges delegates to it exactly like the
/// OCCT inline `WireData()->NbEdges()`, ShapeAnalysis_Wire.hxx L162).
/// CheckConnected is the GAP: it reports "not connected" and keeps the
/// OK status — the OCCT outcome when no connection is found.  GAP: closes
/// with the W2 ShapeAnalysis batch.
#[allow(dead_code)]
struct ShapeAnalysisWire {
    #[allow(dead_code)]
    my_precision: f64,
    my_status: i32,
}

#[allow(dead_code)]
impl ShapeAnalysisWire {
    /// OCCT ShapeAnalysis_Wire().
    fn new() -> Self {
        ShapeAnalysisWire {
            my_precision: 0.0,
            my_status: encode_status(ShapeExtendStatus::Ok),
        }
    }

    /// OCCT Load(sbwd) (ShapeAnalysis_Wire.cxx L142-148): ClearStatuses()
    /// then myWire = sbwd (the aliasing bridge, see the struct doc).
    fn load(&mut self, _sewd: &WireData) {
        self.my_status = encode_status(ShapeExtendStatus::Ok);
    }

    /// OCCT SetPrecision(precision) (ShapeAnalysis_Wire.cxx L201-204).
    fn set_precision(&mut self, precision: f64) {
        self.my_precision = precision;
    }

    /// OCCT NbEdges() (ShapeAnalysis_Wire.hxx L162): WireData()->NbEdges().
    fn nb_edges(&self, sewd: &WireData) -> i32 {
        sewd.nb_edges()
    }

    /// OCCT CheckConnected(num, prec) (ShapeAnalysis_Wire.cxx L693) —
    /// GAP: pending; false keeps the OCCT "no connection found" outcome.
    fn check_connected(&mut self, _sewd: &WireData, _num: i32) -> bool {
        false
    }

    /// OCCT LastCheckStatus(Status) (ShapeAnalysis_Wire.hxx L554).
    fn last_check_status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }
}

/// GAP carrier for OCCT `ShapeAnalysis_BoxBndTreeSelector`
/// (ShapeAnalysis_BoxBndTree.hxx L33-110 + .cxx L29-195, untranslated W2
/// batch).  The hxx inline methods are kept 1:1; Reject/Accept keep the
/// .cxx walk over `ShapeAnalysis::FindBounds` + `BRep_Tool::Pnt` (bridged
/// below).  Architecture bridge #4: the OCCT selector holds the arrwires
/// HArray1 handle; the rcad carrier receives the array on Accept (the
/// `arrwires->ChangeValue(lwire).Reverse()` aliasing must stay visible).
/// GAP: closes with the W2 ShapeAnalysis batch (the BoxBndTree docket row
/// lists FreeBounds as its consumer).
struct ShapeAnalysisBoxBndTreeSelector {
    my_fbox: BndBox,
    my_lbox: BndBox,
    my_shared: bool,
    my_nb: i32,
    my_fvertex: Shape,
    my_lvertex: Shape,
    my_f_pnt: DVec3,
    my_l_pnt: DVec3,
    /// OCCT NCollection_Map<int> myList.
    my_list: HashSet<i32>,
    my_tol: f64,
    my_min3d: f64,
    /// OCCT NCollection_Array1<int> myArrIndices(1, 2) — [First, Last].
    my_arr_indices: [i32; 2],
    my_status: i32,
    my_stop: bool,
}

/// OCCT Accept's local enum { First = 1, Last = 2 } (ShapeAnalysis_BoxBndTree.cxx L51-55).
const SELECTOR_FIRST: usize = 1;
const SELECTOR_LAST: usize = 2;

impl ShapeAnalysisBoxBndTreeSelector {
    /// OCCT ShapeAnalysis_BoxBndTreeSelector(theSeq, theShared)
    /// (ShapeAnalysis_BoxBndTree.hxx L36-47).  Bridge #4: theSeq is passed
    /// on Accept instead of being stored.
    fn new(the_shared: bool) -> Self {
        ShapeAnalysisBoxBndTreeSelector {
            my_fbox: BndBox::new(),
            my_lbox: BndBox::new(),
            my_shared: the_shared,
            my_nb: 0,
            my_fvertex: Shape::null(),
            my_lvertex: Shape::null(),
            my_f_pnt: DVec3::ZERO,
            my_l_pnt: DVec3::ZERO,
            my_list: HashSet::new(),
            my_tol: 1e-7,
            my_min3d: 1e-7,
            my_arr_indices: [0; 2],
            my_status: encode_status(ShapeExtendStatus::Ok),
            my_stop: false,
        }
    }

    /// OCCT DefineBoxes (hxx L49-54).
    fn define_boxes(&mut self, the_fbox: &BndBox, the_lbox: &BndBox) {
        self.my_fbox = the_fbox.clone();
        self.my_lbox = the_lbox.clone();
        self.my_arr_indices = [0; 2];
    }

    /// OCCT DefineVertexes (hxx L56-61).
    fn define_vertexes(&mut self, the_vf: &Shape, the_vl: &Shape) {
        self.my_fvertex = the_vf.clone();
        self.my_lvertex = the_vl.clone();
        self.my_status = encode_status(ShapeExtendStatus::Ok);
    }

    /// OCCT DefinePnt (hxx L63-68).
    fn define_pnt(&mut self, the_f_pnt: DVec3, the_l_pnt: DVec3) {
        self.my_f_pnt = the_f_pnt;
        self.my_l_pnt = the_l_pnt;
        self.my_status = encode_status(ShapeExtendStatus::Ok);
    }

    /// OCCT GetNb (hxx L70).
    fn get_nb(&self) -> i32 {
        self.my_nb
    }

    /// OCCT SetNb (hxx L72).
    fn set_nb(&mut self, the_nb: i32) {
        self.my_nb = the_nb;
    }

    /// OCCT LoadList (hxx L74).
    fn load_list(&mut self, elem: i32) {
        self.my_list.insert(elem);
    }

    /// OCCT SetStop (hxx L76).
    fn set_stop(&mut self) {
        self.my_stop = false;
    }

    /// OCCT SetTolerance (hxx L78-83).
    fn set_tolerance(&mut self, the_tol: f64) {
        self.my_tol = the_tol;
        self.my_min3d = the_tol;
        self.my_status = encode_status(ShapeExtendStatus::Ok);
    }

    /// OCCT ContWire (hxx L85).
    fn cont_wire(&self, nb_wire: i32) -> bool {
        self.my_list.contains(&nb_wire)
    }

    /// OCCT LastCheckStatus (hxx L87-90).
    fn last_check_status(&self, the_status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, the_status)
    }

    /// OCCT Reject (ShapeAnalysis_BoxBndTree.cxx L29-34).
    fn reject(&self, the_bnd: &BndBox) -> bool {
        let fch = self.my_fbox.is_out_box(the_bnd);
        let lch = self.my_lbox.is_out_box(the_bnd);
        fch && lch
    }

    /// OCCT Accept (ShapeAnalysis_BoxBndTree.cxx L38-195).  Bridge #1/#4:
    /// the pool and the aliased arrwires array are passed (the OCCT
    /// selector reads them through stored handles).  The unused_assignments
    /// allowance keeps the OCCT L168-169 `dm1 = dm2` statement, which is
    /// never re-read in OCCT either.
    #[allow(unused_assignments)]
    fn accept(&mut self, brep: &BRep, arrwires: &[Shape], the_obj: i32) -> bool {
        // OCCT L40-44: the range check (Standard_NoSuchObject; rcad panics).
        if the_obj < 1 || the_obj > arrwires.len() as i32 {
            panic!("ShapeAnalysis_BoxBndTreeSelector::Accept : no such object for current index");
        }
        let mut is_accept = false;
        // OCCT L46-49.
        if self.my_list.contains(&the_obj) {
            return false;
        }

        // OCCT L57-59: W = mySeq->Value(theObj); FindBounds(W, V1, V2).
        let w = &arrwires[(the_obj - 1) as usize];
        let (v1, v2) = shape_analysis_find_bounds(brep, w);

        if self.my_shared {
            // OCCT L60-98: the IsSame chain.
            if shape_is_same(&self.my_lvertex, &v1) {
                self.my_status = encode_status(ShapeExtendStatus::Done1);
                is_accept = true;
                self.my_arr_indices[SELECTOR_LAST - 1] = the_obj;
            } else if shape_is_same(&self.my_lvertex, &v2) {
                self.my_status = encode_status(ShapeExtendStatus::Done2);
                is_accept = true;
                self.my_arr_indices[SELECTOR_LAST - 1] = the_obj;
            } else if shape_is_same(&self.my_fvertex, &v2) {
                self.my_status = encode_status(ShapeExtendStatus::Done3);
                is_accept = true;
                self.my_arr_indices[SELECTOR_FIRST - 1] = the_obj;
            } else if shape_is_same(&self.my_fvertex, &v1) {
                self.my_status = encode_status(ShapeExtendStatus::Done4);
                is_accept = true;
                self.my_arr_indices[SELECTOR_FIRST - 1] = the_obj;
            } else {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
            }

            // OCCT L100-112.
            if is_accept {
                self.set_nb(the_obj);
                if self.my_arr_indices[SELECTOR_LAST - 1] != 0 {
                    self.my_stop = true;
                }
                return true;
            } else {
                self.my_stop = false;
            }
        } else {
            // OCCT L117-118: p1/p2 = BRep_Tool::Pnt(V1/V2).
            let p1 = brep_tool_pnt_checked(&v1);
            let p2 = brep_tool_pnt_checked(&v2);

            // OCCT L120-136: the four distances and the res1/res2 picks.
            let tailhead = p1.distance(self.my_l_pnt);
            let tailtail = p2.distance(self.my_l_pnt);
            let headhead = p1.distance(self.my_f_pnt);
            let headtail = p2.distance(self.my_f_pnt);
            let mut dm1 = tailhead;
            let mut dm2 = headtail;
            let mut res1 = 0;
            let mut res2 = 0;
            if tailhead > tailtail {
                res1 = 1;
                dm1 = tailtail;
            }
            if headtail > headhead {
                res2 = 1;
                dm2 = headhead;
            }
            let mut result = res1;
            let min3d = dm1.min(dm2);
            // OCCT L140-143.
            if min3d > self.my_min3d {
                return false;
            }

            // OCCT L145-151.
            let min_ind = if dm1 > dm2 { SELECTOR_FIRST } else { SELECTOR_LAST };
            let max_ind = if dm1 > dm2 { SELECTOR_LAST } else { SELECTOR_FIRST };
            self.my_arr_indices[min_ind - 1] = the_obj;
            if (min3d - self.my_min3d) > REAL_SMALL {
                self.my_arr_indices[max_ind - 1] = 0;
            }

            // OCCT L153-158.
            self.my_min3d = min3d;
            if min3d > self.my_tol {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
                return false;
            }

            // OCCT L160-161.
            let an_obj = if self.my_arr_indices[SELECTOR_LAST - 1] != 0 {
                self.my_arr_indices[SELECTOR_LAST - 1]
            } else {
                self.my_arr_indices[SELECTOR_FIRST - 1]
            };
            self.set_nb(an_obj);

            // OCCT L163-166.
            if min3d == 0.0 && min_ind == SELECTOR_LAST {
                self.my_stop = true;
            }

            // OCCT L168-190: the status/result switch.
            if dm1 > dm2 {
                dm1 = dm2;
                result = res2 + 2;
            }
            if an_obj == the_obj {
                match result {
                    0 => self.my_status = encode_status(ShapeExtendStatus::Done1),
                    1 => self.my_status = encode_status(ShapeExtendStatus::Done2),
                    2 => self.my_status = encode_status(ShapeExtendStatus::Done3),
                    3 => self.my_status = encode_status(ShapeExtendStatus::Done4),
                    _ => {}
                }
            }
            return true;
        }

        // OCCT L194.
        false
    }
}

/// GAP carrier for OCCT `NCollection_UBTree<int, Bnd_Box>` +
/// `NCollection_UBTreeFiller<int, Bnd_Box>` (TKFoundation, untranslated):
/// the int/Bnd_Box tree instantiation used by connectWiresToWiresImpl
/// (FreeBounds.cxx L237-238).  Add/Fill keep the OCCT filler form; Select
/// walks the filled entries in insertion order invoking Reject/Accept and
/// counting the accepted ones — the documented reduction of the OCCT tree
/// traversal (the traversal order is the observable difference; the
/// selector semantics are unchanged).  GAP: closes with the TKFoundation
/// NCollection_UBTree batch.
#[allow(dead_code)]
struct NCollectionUBTree {
    entries: Vec<(i32, BndBox)>,
}

#[allow(dead_code)]
impl NCollectionUBTree {
    /// OCCT NCollection_UBTree<int, Bnd_Box>().
    fn new() -> Self {
        NCollectionUBTree { entries: Vec::new() }
    }
}

/// OCCT NCollection_UBTreeFiller<int, Bnd_Box>(aBBTree) — Add queues the
/// (object, box) pair; Fill() commits them into the tree.
#[allow(dead_code)]
struct NCollectionUBTreeFiller {
    tree: NCollectionUBTree,
}

#[allow(dead_code)]
impl NCollectionUBTreeFiller {
    /// OCCT NCollection_UBTreeFiller(aBBTree).
    fn new(a_bb_tree: NCollectionUBTree) -> Self {
        NCollectionUBTreeFiller { tree: a_bb_tree }
    }

    /// OCCT Add(theObj, theBnd).
    fn add(&mut self, the_obj: i32, the_bnd: BndBox) {
        self.tree.entries.push((the_obj, the_bnd));
    }

    /// OCCT Fill() — commits the queued pairs (the carrier tree stores them
    /// in insertion order).
    fn fill(&mut self) {}

    /// OCCT aBBTree.Select(aSel): the tree traversal invoking Reject/Accept;
    /// returns the number of accepted objects.
    fn select(&mut self, brep: &BRep, a_sel: &mut ShapeAnalysisBoxBndTreeSelector, arrwires: &[Shape]) -> i32 {
        let mut nsel = 0;
        for (the_obj, the_bnd) in &self.tree.entries {
            if a_sel.reject(the_bnd) {
                continue;
            }
            if a_sel.accept(brep, arrwires, *the_obj) {
                nsel += 1;
            }
        }
        nsel
    }
}

/// GAP carrier for OCCT `ShapeAnalysis::FindBounds(shape, V1, V2)`
/// (ShapeAnalysis.cxx L73-102, ShapeAnalysis statics — untranslated W2
/// batch).  The .cxx walk is kept 1:1 over the available primitives
/// (TopExp::Vertices via the brep_algo tool, ShapeAnalysis_Edge via the
/// W1-1 class).  GAP: closes with the W2 ShapeAnalysis batch (the statics
/// row), which replaces this carrier with the 1:1 static.
fn shape_analysis_find_bounds(brep: &BRep, shape: &Shape) -> (Shape, Shape) {
    // OCCT L75-76: V1.Nullify(); V2.Nullify().
    let mut v1 = Shape::null();
    let mut v2 = Shape::null();
    // OCCT L77: ShapeAnalysis_Edge EA.
    let ea = ShapeAnalysisEdge::new();
    if shape.shape_type() == ShapeType::Wire {
        // OCCT L80-82: TopExp::Vertices(W, V1, V2) — invalid work with
        // reversed wires replaced on TopExp.
        let (w_v1, w_v2) = top_exp_vertices_wire(shape);
        v1 = w_v1.unwrap_or_else(Shape::null);
        v2 = w_v2.unwrap_or_else(Shape::null);
    } else if shape.shape_type() == ShapeType::Edge {
        // OCCT L95-96.
        v1 = ea.first_vertex(brep, shape);
        v2 = ea.last_vertex(brep, shape);
    } else if shape.shape_type() == ShapeType::Vertex {
        // OCCT L100: V1 = V2 = shape.
        v1 = shape.clone();
        v2 = shape.clone();
    }
    (v1, v2)
}

/// GAP carrier for OCCT `ShapeBuild_Vertex::CombineVertex` (W1-2 ShapeBuild
/// batch).  Both .cxx overloads are kept 1:1 over the available primitives
/// (BRep_Tool::Pnt/Tolerance via the brep_algo tool, BRep_Builder
/// MakeVertex/UpdateVertex via the tool builders).  GAP: closes with the
/// W1-2 delivery, which replaces the carrier with the 1:1 class.
#[allow(dead_code)]
struct ShapeBuildVertex;

#[allow(dead_code)]
impl ShapeBuildVertex {
    /// OCCT ShapeBuild_Vertex::CombineVertex(V1, V2, tolFactor = 1.0001)
    /// (ShapeBuild_Vertex.cxx L25-34): delegates to the points overload.
    fn combine_vertex(brep: &mut BRep, v1: &Shape, v2: &Shape) -> Shape {
        Self::combine_vertex_tol(brep, v1, v2, 1.0001)
    }

    /// OCCT CombineVertex(V1, V2, tolFactor) — the vertex overload
    /// (ShapeBuild_Vertex.cxx L25-34).
    fn combine_vertex_tol(brep: &mut BRep, v1: &Shape, v2: &Shape, tol_factor: f64) -> Shape {
        Self::combine_vertex_pnts(
            brep,
            &brep_tool_pnt_checked(v1),
            &brep_tool_pnt_checked(v2),
            brep_tool_tolerance(v1),
            brep_tool_tolerance(v2),
            tol_factor,
        )
    }

    /// OCCT CombineVertex(pnt1, pnt2, tol1, tol2, tolFactor)
    /// (ShapeBuild_Vertex.cxx L38-70).
    #[allow(clippy::too_many_arguments)]
    fn combine_vertex_pnts(
        brep: &mut BRep,
        pnt1: &DVec3,
        pnt2: &DVec3,
        tol1: f64,
        tol2: f64,
        tol_factor: f64,
    ) -> Shape {
        // OCCT L41-43: gp_Vec v = pnt2.XYZ() - pnt1.XYZ(); dist = v.Magnitude().
        let v = *pnt2 - *pnt1;
        let dist = v.length();
        // OCCT L45-58: the containment / barycenter cases.
        let pos;
        let tol;
        if dist + tol2 <= tol1 {
            pos = *pnt1;
            tol = tol1;
        } else if dist + tol1 <= tol2 {
            pos = *pnt2;
            tol = tol2;
        } else {
            tol = 0.5 * (dist + tol1 + tol2);
            // szv#4:S4163:12Mar99 anti-exception
            let s = if dist > 0.0 { (tol2 - tol1) / dist } else { 0.0 };
            pos = 0.5 * ((1.0 - s) * *pnt1 + (1.0 + s) * *pnt2);
        }

        // OCCT L61-64: TopoDS_Vertex V; BRep_Builder B; B.MakeVertex(V, pos,
        // tolFactor * tol) — the builder MakeVertex(P, Tol) is the plain
        // MakeVertex + UpdateVertex (tolerance keeps the max).
        let mut nv = brep.add_tvertex_unique(pos);
        builder_update_vertex_point_tol(&mut nv, pos, tol_factor * tol);
        nv
    }
}

/// GAP carrier for OCCT `ShapeBuild_Edge::CopyReplaceVertices(edge, V1, V2)`
/// (ShapeBuild_Edge.cxx L59-152, W1-2 ShapeBuild batch).  The .cxx walk is
/// kept 1:1 over the available primitives (TopoDS_Iterator via brep_tool,
/// EmptyCopied via the brep_algo tool, BRep_Builder Add via the tool
/// builders, CopyRanges via the existing shape_build/edge.rs compat API —
/// the OCCT defaults alpha = 0, beta = 1).  GAP: closes with the W1-2
/// delivery, which replaces the carrier with the 1:1 method.
fn shape_build_edge_copy_replace_vertices(
    brep: &mut BRep,
    edge: &Shape,
    v1: Option<&Shape>,
    v2: Option<&Shape>,
) -> Shape {
    // OCCT L63: NCollection_Sequence<TopoDS_Shape> aNMVertices.
    let mut a_nm_vertices: Vec<Shape> = Vec::new();
    // OCCT L64: newV1 = V1, newV2 = V2.
    let mut new_v1: Option<Shape> = v1.cloned();
    let mut new_v2: Option<Shape> = v2.cloned();
    // OCCT L65-98: when an extremity is missing, take it from the edge.
    if new_v1.is_none() || new_v2.is_none() {
        // OCCT L67-75: TopoDS_Iterator it; cumOri = the edge is FWD or REV.
        let cum_ori = edge.orientation == Orientation::Forward
            || edge.orientation == Orientation::Reversed;
        for sv in iter_subshapes(brep, edge, cum_ori, true) {
            if sv.shape_type() != ShapeType::Vertex {
                continue;
            }
            // OCCT L79-96.
            if sv.orientation == Orientation::Forward {
                if new_v1.is_none() {
                    new_v1 = Some(sv.clone());
                }
            } else if sv.orientation == Orientation::Reversed {
                if new_v2.is_none() {
                    new_v2 = Some(sv.clone());
                }
            } else if v1.is_none() && v2.is_none() {
                a_nm_vertices.push(sv.clone());
            }
        }
    }
    // OCCT L99-100: newV1.Orientation(TopAbs_FORWARD);
    // newV2.Orientation(TopAbs_REVERSED) (a no-op on the OCCT null handle —
    // the rcad None stays None).
    let new_v1 = new_v1.map(|mut v| {
        v.orientation = Orientation::Forward;
        v
    });
    let new_v2 = new_v2.map(|mut v| {
        v.orientation = Orientation::Reversed;
        v
    });

    // OCCT L103-105: sh = edge.EmptyCopied(); E = TopoDS::Edge(sh).
    let mut e = empty_copied(edge);

    // OCCT L107-117: B.Add(E, newV1 / newV2).
    if let Some(v) = &new_v1 {
        if !shape_is_null(v) {
            builder_add_edge_vertex(&mut e, v);
        }
    }
    if let Some(v) = &new_v2 {
        if !shape_is_null(v) {
            builder_add_edge_vertex(&mut e, v);
        }
    }

    // OCCT L120-124: the internal/external vertices.
    for i in 1..=a_nm_vertices.len() {
        builder_add_edge_vertex(&mut e, &a_nm_vertices[i - 1]);
    }

    // OCCT L130: CopyRanges(E, edge) — the defaults alpha = 0, beta = 1.
    ShapeBuildEdge.copy_ranges(brep, &e, edge, 0.0, 1.0);

    e
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
        sas.check_oriented_shells(&a_tmp_shell, true, checkinternaledges);

        // OCCT L120-130.
        if sas.has_free_edges() {
            let see = ShapeExtendExplorer;
            // OCCT L122-124: edges = see.SeqFromCompound(sas.FreeEdges(), false).
            let mut edges = see.seq_from_compound(brep, &sas.free_edges(), false);

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

    // OCCT L233-235: saw = new ShapeAnalysis_Wire; Load; SetPrecision.
    let mut saw = ShapeAnalysisWire::new();
    saw.load(&sewd);
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
        let (tr_v1, tr_v2) = shape_analysis_find_bounds(brep, tr_w);
        // OCCT L248-252: trP1/trP2 = BRep_Tool::Pnt(trV1/trV2); aBox.Set(trP1)
        // (the OCCT Set(P) builds the single-point box — BndBox::from_point);
        // aBox.Add(trP2); aBox.SetGap(tolerance).
        let tr_p1 = brep_tool_pnt_checked(&tr_v1);
        let tr_p2 = brep_tool_pnt_checked(&tr_v2);
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
        // OCCT L267-270: FVBox / LVBox / Vf / Vl.
        let vf = sae.first_vertex(brep, &sewd.edge(1));
        let vl = sae.last_vertex(brep, &sewd.edge(sewd.nb_edges()));

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
            sewd.add_wire_data(&acurwd, if tail { 0 } else { 1 });
        } else {
            // OCCT L320-359: the CheckConnected fallback walk.
            for i in 1..=saw.nb_edges(&sewd) {
                if saw.check_connected(&sewd, i) {
                    // OCCT L324-327: n2 / n1 / E1 / E2.
                    let n2 = i;
                    let n1 = if n2 > 1 { n2 - 1 } else { saw.nb_edges(&sewd) };
                    let e1 = sewd.edge(n1);
                    let e2 = sewd.edge(n2);

                    // OCCT L329-331: Vprev / Vfol.
                    let vprev = sae.last_vertex(brep, &e1);
                    let vfol = sae.first_vertex(brep, &e2);

                    // OCCT L333-341: V = Vprev or the combined vertex
                    // (OCCT L339-340: ShapeBuild_Vertex sbv; V =
                    // sbv.CombineVertex(Vprev, Vfol) — the rcad unit struct
                    // has no instance, the call is through the type).
                    let v = if saw.last_check_status(ShapeExtendStatus::Done1) {
                        vprev.clone()
                    } else {
                        ShapeBuildVertex::combine_vertex(brep, &vprev, &vfol)
                    };
                    // OCCT L342-343: vertices.Bind(Vprev, V); Bind(Vfol, V).
                    vertices.insert((vprev.ptr_id(), vprev.location), v.clone());
                    vertices.insert((vfol.ptr_id(), vfol.location), v.clone());

                    // OCCT L345-357: ShapeBuild_Edge sbe; the Set calls
                    // (the OCCT null TopoDS_Vertex() argument maps to None).
                    if saw.nb_edges(&sewd) < 2 {
                        let e2r =
                            shape_build_edge_copy_replace_vertices(brep, &e2, Some(&v), Some(&v));
                        sewd.set_edge(&e2r, n2);
                    } else {
                        let e2r = shape_build_edge_copy_replace_vertices(
                            brep,
                            &e2,
                            Some(&v),
                            None,
                        );
                        sewd.set_edge(&e2r, n2);
                        if !saw.last_check_status(ShapeExtendStatus::Done1) {
                            let e1r = shape_build_edge_copy_replace_vertices(
                                brep,
                                &e1,
                                None,
                                Some(&v),
                            );
                            sewd.set_edge(&e1r, n1);
                        }
                    }
                }
            }

            // OCCT L361: TopoDS_Wire wire = sewd->Wire() (the Closed flag
            // write below mutates the TShape in place, not the wrapper).
            let wire = sewd.wire(brep);
            // OCCT L362-368: the closed-flag checks.
            if is_used_manifold_mode {
                if !saw.check_connected(&sewd, 1)
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
            sewd.clear();
            sewd.set_manifold_mode(is_used_manifold_mode);

            // OCCT L400-415: pick the next not-yet-connected wire.
            lwire = -1;
            for i in 1..=arrwires.len() as i32 {
                if !a_sel.cont_wire(i) {
                    lwire = i;
                    sewd.add_wire(brep, &arrwires[(lwire - 1) as usize], 0);
                    a_sel.load_list(lwire);

                    if sewd.nb_edges() > 0 {
                        break;
                    }
                    sewd.clear();
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
