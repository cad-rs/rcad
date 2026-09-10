//! 1:1 translation of OCCT `ShapeFix_Wire`
//! (`TKShHealing/ShapeFix/ShapeFix_Wire.hxx` L17-548 + `ShapeFix_Wire.cxx`
//! L1-4479 + `ShapeFix_Wire_1.cxx` L1-2006 + `ShapeFix_Wire.lxx` L1-427,
//! docket row `ShapeFix_Wire`).
//!
//! Function-count equation (OCCT Wire.cxx + Wire_1.cxx + lxx = rcad wire/):
//! constructors + section 1 (cxx L117-281): `ShapeFix_Wire()` /
//! `ShapeFix_Wire(wire,face,prec)` / `SetPrecision` / `SetMaxTailAngle` /
//! `SetMaxTailWidth` / `ClearModes` / `ClearStatuses` / `Init(w,f,p)` /
//! `Init(saw)` / `Load(wire)` / `Load(sbwd)` / `NbEdges` — 12; API level
//! (cxx L297-1343): `Perform` / `FixReorder(bool)` / `FixSmall(bool,prec)` /
//! `FixConnected(prec)` / `FixEdgeCurves` / `FixDegenerated()` /
//! `FixSelfIntersection` / `FixLacking(bool)` / `FixClosed` — 9; advanced
//! level (cxx L1351-4289 + Wire_1.cxx): `FixReorder(wi)` / `FixSmall(num)` /
//! `FixConnected(num)` / `FixSeam` / `FixShifted` / `FixDegenerated(num)` /
//! `FixSelfIntersectingEdge` / `FixIntersectingEdges(num)` /
//! `FixIntersectingEdges(num1,num2)` / `FixLacking(num,force)` /
//! `FixNotchedEdges` / `FixTails` / `UpdateWire` / `FixDummySeam` /
//! `FixGaps3d` / `FixGaps2d` / `FixGap3d` / `FixGap2d` — 18; file statics
//! (cxx + Wire_1.cxx): `UpdateEdgeUVPoints` / `TryNewPCurve` /
//! `howMuchPCurves` / `RemoveLoop(tolfact,prec,RemoveLoop3d)` /
//! `RemoveLoop(E1,E2)` / `ComputeLocalDeviation` / `TryBendingPCurve` /
//! `CopyReversePcurves` / `AdjustOnPeriodic3d` / `AdjustOnPeriodic2d` — 10;
//! lxx inlines (lxx L24-427): `SetFace` x2 / `SetSurface` x3 / `IsLoaded` /
//! `IsReady` / `Wire` / `WireAPIMake` / `Analyzer` / `WireData` / `Face` /
//! `ModifyTopologyMode` / `ModifyGeometryMode` / `ModifyRemoveLoopMode` /
//! `ClosedWireMode` / `PreferencePCurveMode` / `FixGapsByRangesMode` /
//! `FixReorderMode` / `FixSmallMode` / `FixConnectedMode` /
//! `FixEdgeCurvesMode` / `FixDegeneratedMode` / `FixReversed2dMode` /
//! `FixRemovePCurveMode` / `FixRemoveCurve3dMode` / `FixAddPCurveMode` /
//! `FixAddCurve3dMode` / `FixSeamMode` / `FixShiftedMode` /
//! `FixSameParameterMode` / `FixVertexToleranceMode` / `FixLackingMode` /
//! `FixSelfIntersectionMode` / `FixGaps3dMode` / `FixGaps2dMode` /
//! `FixNotchedEdgesMode` / `FixSelfIntersectingEdgeMode` /
//! `FixIntersectingEdgesMode` / `FixNonAdjacentIntersectingEdgesMode` /
//! `FixTailMode` / `StatusReorder` / `StatusSmall` / `StatusConnected` /
//! `StatusEdgeCurves` / `StatusDegenerated` / `StatusLacking` /
//! `StatusSelfIntersection` / `StatusGaps3d` / `StatusGaps2d` /
//! `StatusClosed` / `StatusNotches` / `StatusFixTails` / `LastFixStatus` /
//! `FixEdgeTool` / `StatusRemovedSegment` — 56.
//! 105 OCCT functions = 105 rcad functions (the one extra rcad method,
//! `fix_intersecting_edges_impl`, is the Rust-internal split of the
//! `FixIntersectingEdges(num)` body documented in `fix_intersect.rs`).
//!
//! Architecture bridges:
//! 1. `ShapeFix_Root` inheritance — the rcad composition field `base` (Rust
//!    has no inheritance); the Context()/Precision()/MinTolerance()/
//!    MaxTolerance()/SendWarning root members are reached through it.
//! 2. `occ::handle<ShapeAnalysis_Wire>` — the rcad owned value
//!    `my_analyzer` (the handle copy in `Init(saw)` becomes a move).
//! 3. `BRep` pool argument — OCCT `BRep_Tool`/`BRep_Builder` read and write
//!    TShapes through global accessors; the rcad equivalents take
//!    `brep: &mut BRep` (the ShapeFix_Edge precedent).
//! 4. OCCT overloads are disambiguated with descriptive suffixes (the
//!    ShapeAnalysis_Wire precedent): `Init(saw)` -> `init_analyzer`,
//!    `Load(sbwd)` -> `load_wire_data`, `SetFace(face,sa)` ->
//!    `set_face_with_surface`, `SetSurface(surf)` -> `set_surface_geom`,
//!    `SetSurface(surf,loc)` -> `set_surface_geom_loc`,
//!    `FixReorder(wi)` -> `fix_reorder_ordered`, `FixSmall(num,...)` ->
//!    `fix_small_edge`, `FixConnected(num,...)` -> `fix_connected_edge`,
//!    `FixDegenerated(num)` -> `fix_degenerated_edge`, `FixLacking(num,..)`
//!    -> `fix_lacking_edge`, `FixIntersectingEdges(num1,num2)` ->
//!    `fix_intersecting_edges_pair`, `RemoveLoop(E1,E2)` ->
//!    `remove_loop_split`.
//! 5. OCCT `try { OCC_CATCH_SIGNALS } catch (Standard_Failure)` blocks —
//!    the failure arms are annotated at site (the tranche-1 precedent); rcad
//!    raises no exceptions in these bodies, so the OK-arm flow is preserved.
//! 6. `ShapeFix_IntersectionTool` / `ShapeFix_SplitTool` — W3 docket rows
//!    not yet landed; the GAP re-hosts at the bottom of `fix_adv.rs` keep
//!    the call shapes and OCCT's failure/no-op paths (module doc of
//!    `fix_adv.rs`).
//! 7. `GeomAPI_ExtremaCurveCurve` / `GeomAPI_ProjectPointOnCurve` /
//!    `Geom2dAPI_*` / `Geom2dInt_GInter` (the general curve-curve
//!    intersection) / `Approx_Curve3d` / `Approx_Curve2d` /
//!    `GeomConvert_CompCurveToBSplineCurve` / `GeomAPI::To2d`/`To3d` —
//!    TKGeomBase/TKG3d leaves not translated; the GAP re-hosts at the
//!    bottom of `fix_gaps.rs` / `fix_adv.rs` keep the call shapes and the
//!    OCCT `!IsDone()` failure branches (the `wire_checks.rs` bridge #5
//!    precedent).

mod fix_adv;
mod fix_api;
mod fix_gaps;
mod fix_intersect;
mod wire_statics;

use glam::{DVec2, DVec3};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, Orientation, TShape, tshape_flags};

use crate::shhealing::shape_analysis::wire::ShapeAnalysisWire;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::edge::ShapeFixEdge;
use crate::shhealing::shape_fix::root::ShapeFixRoot;

/// OCCT `RealLast()` — the greatest finite double.
pub(crate) const REAL_LAST: f64 = f64::MAX;

// ---------------------------------------------------------------------------
// BRep_Tool re-hosts shared by the wire submodules (bridge #3; the
// shape_analysis/wire.rs precedents).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Pnt(V).
pub(crate) fn brep_tool_pnt(v: &Shape) -> DVec3 {
    match v.data.as_ref() {
        TShape::Vertex(vd) => vd.point,
        _ => DVec3::ZERO,
    }
}

/// OCCT BRep_Tool::Tolerance(shape).
pub(crate) fn brep_tool_tolerance(s: &Shape) -> f64 {
    match s.data.as_ref() {
        TShape::Vertex(vd) => vd.tolerance,
        TShape::Edge(ed) => ed.tolerance,
        _ => 0.0,
    }
}

/// OCCT BRep_Tool::Degenerated(E).
pub(crate) fn brep_tool_degenerated(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.degenerated,
        _ => false,
    }
}

/// OCCT BRep_Tool::SameParameter(E).
pub(crate) fn brep_tool_same_parameter(e: &Shape) -> bool {
    match e.data.as_ref() {
        TShape::Edge(ed) => ed.same_parameter,
        _ => false,
    }
}

/// OCCT `TopoDS_Shape::Oriented(o)` — a copy with the given orientation.
pub(crate) fn shape_oriented(s: &Shape, o: Orientation) -> Shape {
    let mut r = s.clone();
    r.orientation = o;
    r
}

/// OCCT `TopoDS_Shape::Free()` — the TShape FREE flag (carried per-variant
/// in the rcad TShape data records).
pub(crate) fn shape_free(brep: &BRep, s: &Shape) -> bool {
    match brep.tshapes.get(s.index).map(|ts| ts.as_ref()) {
        Some(TShape::Vertex(vd)) => vd.flags & tshape_flags::FREE != 0,
        Some(TShape::Edge(ed)) => ed.flags & tshape_flags::FREE != 0,
        Some(TShape::Wire(wd)) => wd.flags & tshape_flags::FREE != 0,
        Some(TShape::Face(fd)) => fd.flags & tshape_flags::FREE != 0,
        _ => false,
    }
}

/// OCCT ShapeFix_Wire (hxx L87-543): a set of tools for repairing a wire.
pub struct ShapeFixWire {
    /// OCCT the ShapeFix_Root base subobject (bridge #1).
    pub base: ShapeFixRoot,
    /// OCCT myFixEdge (hxx L469).
    pub(crate) my_fix_edge: ShapeFixEdge,
    /// OCCT myAnalyzer (hxx L470).
    pub(crate) my_analyzer: ShapeAnalysisWire,
    /// OCCT myGeomMode (hxx L471).
    pub(crate) my_geom_mode: bool,
    /// OCCT myTopoMode (hxx L472).
    pub(crate) my_topo_mode: bool,
    /// OCCT myClosedMode (hxx L473).
    pub(crate) my_closed_mode: bool,
    /// OCCT myPreference2d (hxx L474).
    pub(crate) my_preference2d: bool,
    /// OCCT myFixGapsByRanges (hxx L475).
    pub(crate) my_fix_gaps_by_ranges: bool,
    /// OCCT myFixReversed2dMode (hxx L476).
    pub(crate) my_fix_reversed2d_mode: i32,
    /// OCCT myFixRemovePCurveMode (hxx L477).
    pub(crate) my_fix_remove_pcurve_mode: i32,
    /// OCCT myFixAddPCurveMode (hxx L478).
    pub(crate) my_fix_add_pcurve_mode: i32,
    /// OCCT myFixRemoveCurve3dMode (hxx L479).
    pub(crate) my_fix_remove_curve3d_mode: i32,
    /// OCCT myFixAddCurve3dMode (hxx L480).
    pub(crate) my_fix_add_curve3d_mode: i32,
    /// OCCT myFixSeamMode (hxx L481).
    pub(crate) my_fix_seam_mode: i32,
    /// OCCT myFixShiftedMode (hxx L482).
    pub(crate) my_fix_shifted_mode: i32,
    /// OCCT myFixSameParameterMode (hxx L483).
    pub(crate) my_fix_same_parameter_mode: i32,
    /// OCCT myFixVertexToleranceMode (hxx L484).
    pub(crate) my_fix_vertex_tolerance_mode: i32,
    /// OCCT myFixNotchedEdgesMode (hxx L485).
    pub(crate) my_fix_notched_edges_mode: i32,
    /// OCCT myFixSelfIntersectingEdgeMode (hxx L486).
    pub(crate) my_fix_self_intersecting_edge_mode: i32,
    /// OCCT myFixIntersectingEdgesMode (hxx L487).
    pub(crate) my_fix_intersecting_edges_mode: i32,
    /// OCCT myFixNonAdjacentIntersectingEdgesMode (hxx L488).
    pub(crate) my_fix_non_adjacent_intersecting_edges_mode: i32,
    /// OCCT myFixTailMode (hxx L489).
    pub(crate) my_fix_tail_mode: i32,
    /// OCCT myRemoveLoopMode (hxx L490).
    pub(crate) my_remove_loop_mode: i32,
    /// OCCT myFixReorderMode (hxx L491).
    pub(crate) my_fix_reorder_mode: i32,
    /// OCCT myFixSmallMode (hxx L492).
    pub(crate) my_fix_small_mode: i32,
    /// OCCT myFixConnectedMode (hxx L493).
    pub(crate) my_fix_connected_mode: i32,
    /// OCCT myFixEdgeCurvesMode (hxx L494).
    pub(crate) my_fix_edge_curves_mode: i32,
    /// OCCT myFixDegeneratedMode (hxx L495).
    pub(crate) my_fix_degenerated_mode: i32,
    /// OCCT myFixSelfIntersectionMode (hxx L496).
    pub(crate) my_fix_self_intersection_mode: i32,
    /// OCCT myFixLackingMode (hxx L497).
    pub(crate) my_fix_lacking_mode: i32,
    /// OCCT myFixGaps3dMode (hxx L498).
    pub(crate) my_fix_gaps3d_mode: i32,
    /// OCCT myFixGaps2dMode (hxx L499).
    pub(crate) my_fix_gaps2d_mode: i32,
    /// OCCT myLastFixStatus (hxx L500).
    pub(crate) my_last_fix_status: i32,
    /// OCCT myStatusReorder (hxx L501).
    pub(crate) my_status_reorder: i32,
    /// OCCT myStatusSmall (hxx L502).
    pub(crate) my_status_small: i32,
    /// OCCT myStatusConnected (hxx L503).
    pub(crate) my_status_connected: i32,
    /// OCCT myStatusEdgeCurves (hxx L504).
    pub(crate) my_status_edge_curves: i32,
    /// OCCT myStatusDegenerated (hxx L505).
    pub(crate) my_status_degenerated: i32,
    /// OCCT myStatusClosed (hxx L506).
    pub(crate) my_status_closed: i32,
    /// OCCT myStatusSelfIntersection (hxx L507).
    pub(crate) my_status_self_intersection: i32,
    /// OCCT myStatusLacking (hxx L508).
    pub(crate) my_status_lacking: i32,
    /// OCCT myStatusGaps3d (hxx L509).
    pub(crate) my_status_gaps3d: i32,
    /// OCCT myStatusGaps2d (hxx L510).
    pub(crate) my_status_gaps2d: i32,
    /// OCCT myStatusRemovedSegment (hxx L511).
    pub(crate) my_status_removed_segment: bool,
    /// OCCT myStatusNotches (hxx L512).
    pub(crate) my_status_notches: i32,
    /// OCCT myStatusFixTails (hxx L513).
    pub(crate) my_status_fix_tails: i32,
    /// OCCT myMaxTailAngleSine (hxx L514).
    pub(crate) my_max_tail_angle_sine: f64,
    /// OCCT myMaxTailWidth (hxx L515).
    pub(crate) my_max_tail_width: f64,
}

impl Default for ShapeFixWire {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixWire {
    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Wire.cxx L117-126 — constructor.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Wire::ShapeFix_Wire() (cxx L117-126): empty constructor,
    /// creates a clear object with default flags.
    pub fn new() -> Self {
        let mut this = ShapeFixWire {
            base: ShapeFixRoot::new(),
            my_fix_edge: ShapeFixEdge::new(),
            my_analyzer: ShapeAnalysisWire::new(),
            my_geom_mode: false,
            my_topo_mode: false,
            my_closed_mode: false,
            my_preference2d: false,
            my_fix_gaps_by_ranges: false,
            my_fix_reversed2d_mode: -1,
            my_fix_remove_pcurve_mode: -1,
            my_fix_add_pcurve_mode: -1,
            my_fix_remove_curve3d_mode: -1,
            my_fix_add_curve3d_mode: -1,
            my_fix_seam_mode: -1,
            my_fix_shifted_mode: -1,
            my_fix_same_parameter_mode: -1,
            my_fix_vertex_tolerance_mode: -1,
            my_fix_notched_edges_mode: -1,
            my_fix_self_intersecting_edge_mode: -1,
            my_fix_intersecting_edges_mode: -1,
            my_fix_non_adjacent_intersecting_edges_mode: -1,
            my_fix_tail_mode: 0,
            my_remove_loop_mode: -1,
            my_fix_reorder_mode: -1,
            my_fix_small_mode: -1,
            my_fix_connected_mode: -1,
            my_fix_edge_curves_mode: -1,
            my_fix_degenerated_mode: -1,
            my_fix_self_intersection_mode: -1,
            my_fix_lacking_mode: -1,
            my_fix_gaps3d_mode: -1,
            my_fix_gaps2d_mode: -1,
            my_last_fix_status: 0,
            my_status_reorder: 0,
            my_status_small: 0,
            my_status_connected: 0,
            my_status_edge_curves: 0,
            my_status_degenerated: 0,
            my_status_closed: 0,
            my_status_self_intersection: 0,
            my_status_lacking: 0,
            my_status_gaps3d: 0,
            my_status_gaps2d: 0,
            my_status_removed_segment: false,
            my_status_notches: 0,
            my_status_fix_tails: 0,
            my_max_tail_angle_sine: 0.0,
            my_max_tail_width: -1.0,
        };
        this.clear_modes();
        this.clear_statuses();
        this.my_status_removed_segment = false;
        this
    }

    /// OCCT ShapeFix_Wire::ShapeFix_Wire(wire, face, prec) (cxx L130-140):
    /// creates a new object with default flags and prepares it for use.
    pub fn with_wire_face(brep: &mut BRep, wire: &Shape, face: &Shape, prec: f64) -> Self {
        let mut this = ShapeFixWire {
            base: ShapeFixRoot::new(),
            my_fix_edge: ShapeFixEdge::new(),
            my_analyzer: ShapeAnalysisWire::new(),
            my_geom_mode: false,
            my_topo_mode: false,
            my_closed_mode: false,
            my_preference2d: false,
            my_fix_gaps_by_ranges: false,
            my_fix_reversed2d_mode: -1,
            my_fix_remove_pcurve_mode: -1,
            my_fix_add_pcurve_mode: -1,
            my_fix_remove_curve3d_mode: -1,
            my_fix_add_curve3d_mode: -1,
            my_fix_seam_mode: -1,
            my_fix_shifted_mode: -1,
            my_fix_same_parameter_mode: -1,
            my_fix_vertex_tolerance_mode: -1,
            my_fix_notched_edges_mode: -1,
            my_fix_self_intersecting_edge_mode: -1,
            my_fix_intersecting_edges_mode: -1,
            my_fix_non_adjacent_intersecting_edges_mode: -1,
            my_fix_tail_mode: 0,
            my_remove_loop_mode: -1,
            my_fix_reorder_mode: -1,
            my_fix_small_mode: -1,
            my_fix_connected_mode: -1,
            my_fix_edge_curves_mode: -1,
            my_fix_degenerated_mode: -1,
            my_fix_self_intersection_mode: -1,
            my_fix_lacking_mode: -1,
            my_fix_gaps3d_mode: -1,
            my_fix_gaps2d_mode: -1,
            my_last_fix_status: 0,
            my_status_reorder: 0,
            my_status_small: 0,
            my_status_connected: 0,
            my_status_edge_curves: 0,
            my_status_degenerated: 0,
            my_status_closed: 0,
            my_status_self_intersection: 0,
            my_status_lacking: 0,
            my_status_gaps3d: 0,
            my_status_gaps2d: 0,
            my_status_removed_segment: false,
            my_status_notches: 0,
            my_status_fix_tails: 0,
            my_max_tail_angle_sine: 0.0,
            my_max_tail_width: -1.0,
        };
        this.clear_modes();
        this.base.set_max_tolerance(prec);
        this.my_status_removed_segment = false;
        this.init(brep, wire, face, prec);
        this
    }

    /// OCCT ShapeFix_Wire::SetPrecision (cxx L144-148): sets the working
    /// precision (to root and to analyzer).
    pub fn set_precision(&mut self, prec: f64) {
        self.base.set_precision(prec);
        self.my_analyzer.set_precision(prec);
    }

    /// OCCT ShapeFix_Wire::SetMaxTailAngle (cxx L152-156): sets the maximal
    /// allowed angle of the tails in radians.
    pub fn set_max_tail_angle(&mut self, the_max_tail_angle: f64) {
        self.my_max_tail_angle_sine = the_max_tail_angle.sin();
        self.my_max_tail_angle_sine = if self.my_max_tail_angle_sine >= 0.0 {
            self.my_max_tail_angle_sine
        } else {
            0.0
        };
    }

    /// OCCT ShapeFix_Wire::SetMaxTailWidth (cxx L160-163): sets the maximal
    /// allowed width of the tails.
    pub fn set_max_tail_width(&mut self, the_max_tail_width: f64) {
        self.my_max_tail_width = the_max_tail_width;
    }

    /// OCCT ShapeFix_Wire::ClearModes (cxx L167-202): sets all modes to
    /// default.
    pub fn clear_modes(&mut self) {
        self.my_topo_mode = false;
        self.my_geom_mode = true;
        self.my_closed_mode = true;
        self.my_preference2d = true;
        self.my_fix_gaps_by_ranges = false;

        self.my_remove_loop_mode = -1;

        self.my_fix_reversed2d_mode = -1;
        self.my_fix_remove_pcurve_mode = -1;
        self.my_fix_remove_curve3d_mode = -1;
        self.my_fix_add_pcurve_mode = -1;
        self.my_fix_add_curve3d_mode = -1;
        self.my_fix_seam_mode = -1;
        self.my_fix_shifted_mode = -1;
        self.my_fix_same_parameter_mode = -1;
        self.my_fix_vertex_tolerance_mode = -1;

        self.my_fix_notched_edges_mode = -1;
        self.my_fix_self_intersecting_edge_mode = -1;
        self.my_fix_intersecting_edges_mode = -1;
        self.my_fix_non_adjacent_intersecting_edges_mode = -1;
        self.my_fix_tail_mode = 0;

        self.my_fix_reorder_mode = -1;
        self.my_fix_small_mode = -1;
        self.my_fix_connected_mode = -1;
        self.my_fix_edge_curves_mode = -1;
        self.my_fix_degenerated_mode = -1;
        self.my_fix_self_intersection_mode = -1;
        self.my_fix_lacking_mode = -1;
        self.my_fix_gaps3d_mode = -1;
        self.my_fix_gaps2d_mode = -1;
    }

    /// OCCT ShapeFix_Wire::ClearStatuses (cxx L206-224): clears all statuses.
    pub fn clear_statuses(&mut self) {
        let empty_status = encode_status(ShapeExtendStatus::Ok);

        self.my_last_fix_status = empty_status;

        self.my_status_reorder = empty_status;
        self.my_status_small = empty_status;
        self.my_status_connected = empty_status;
        self.my_status_edge_curves = empty_status;
        self.my_status_degenerated = empty_status;
        self.my_status_self_intersection = empty_status;
        self.my_status_lacking = empty_status;
        self.my_status_gaps3d = empty_status;
        self.my_status_gaps2d = empty_status;
        self.my_status_closed = empty_status;
        self.my_status_notches = empty_status;
        self.my_status_fix_tails = empty_status;
    }

    /// OCCT ShapeFix_Wire::Init(wire, face, prec) (cxx L228-233): loads the
    /// analyzer with all the data for the wire and face and drops all fixing
    /// statuses.
    pub fn init(&mut self, brep: &mut BRep, wire: &Shape, face: &Shape, prec: f64) {
        self.load(brep, wire);
        self.set_face(brep, face);
        self.set_precision(prec);
    }

    /// OCCT ShapeFix_Wire::Init(saw) (cxx L237-243): loads the analyzer with
    /// all the data already prepared and drops all fixing statuses.
    pub fn init_analyzer(&mut self, saw: ShapeAnalysisWire) {
        self.clear_statuses();
        self.my_analyzer = saw;
        self.base.my_shape = Shape::null();
    }

    /// OCCT ShapeFix_Wire::Load(wire) (cxx L247-260): loads data for the
    /// wire and drops all fixing statuses.
    pub fn load(&mut self, brep: &mut BRep, wire: &Shape) {
        self.clear_statuses();

        let mut w = wire.clone();
        if let Some(ctx) = self.base.my_context.as_mut() {
            let s = ctx.apply(brep, wire, rcad_kernel::topods::ShapeType::Shape);
            w = s;
        }

        self.my_analyzer.load(brep, &w);
        self.base.my_shape = wire.clone();
    }

    /// OCCT ShapeFix_Wire::Load(sbwd) (cxx L264-273): loads data for the
    /// wire and drops all fixing statuses.
    pub fn load_wire_data(&mut self, brep: &mut BRep, sbwd: crate::shhealing::shape_extend::wire_data::WireData) {
        self.clear_statuses();
        self.my_analyzer.load_wire_data(sbwd);
        if self.base.my_context.is_some() {
            self.update_wire(brep);
        }
        self.base.my_shape = Shape::null();
    }

    /// OCCT ShapeFix_Wire::NbEdges (cxx L277-281): returns the number of
    /// edges in the working wire.
    pub fn nb_edges(&self) -> i32 {
        match self.my_analyzer.wire_data() {
            Some(sbwd) => sbwd.nb_edges(),
            None => 0,
        }
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Wire.lxx — the inline accessors (lxx L24-427).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Wire::SetFace(face) (lxx L24-27): sets the working face
    /// for the wire.
    pub fn set_face(&mut self, brep: &mut BRep, face: &Shape) {
        self.my_analyzer.set_face(brep, face);
    }

    /// OCCT ShapeFix_Wire::SetSurface(surf) (lxx L31-34): sets the surface
    /// for the wire.
    pub fn set_surface_geom(
        &mut self,
        brep: &mut BRep,
        surf: &rcad_kernel::geom::Surface3,
    ) {
        self.my_analyzer.set_surface_geom(brep, surf);
    }

    /// OCCT ShapeFix_Wire::SetFace(theFace, theSurfaceAnalysis) (lxx L38-42):
    /// sets the working face for the wire and the surface analysis object.
    pub fn set_face_with_surface(
        &mut self,
        brep: &mut BRep,
        the_face: &Shape,
        the_surface_analysis: crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface,
    ) {
        self.my_analyzer
            .set_face_with_surface(the_face, the_surface_analysis);
    }

    /// OCCT ShapeFix_Wire::SetSurface(theSurfaceAnalysis) (lxx L46-49): sets
    /// the surface analysis for the wire.
    pub fn set_surface(
        &mut self,
        the_surface_analysis: crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface,
    ) {
        self.my_analyzer.set_surface(the_surface_analysis);
    }

    /// OCCT ShapeFix_Wire::SetSurface(surf, loc) (lxx L53-57): sets the
    /// surface for the wire.
    pub fn set_surface_geom_loc(
        &mut self,
        brep: &mut BRep,
        surf: &rcad_kernel::geom::Surface3,
        loc: u32,
    ) {
        self.my_analyzer.set_surface_geom_loc(brep, surf, loc);
    }

    /// OCCT ShapeFix_Wire::IsLoaded (lxx L61-64): tells if the wire is
    /// loaded.
    pub fn is_loaded(&self) -> bool {
        self.my_analyzer.is_loaded()
    }

    /// OCCT ShapeFix_Wire::IsReady (lxx L68-71): tells if the wire and face
    /// are loaded.
    pub fn is_ready(&self) -> bool {
        self.my_analyzer.is_ready()
    }

    /// OCCT ShapeFix_Wire::Wire (lxx L75-78): makes the resulting Wire (by
    /// the basic BRep_Builder).
    pub fn wire(&self, brep: &mut BRep) -> Shape {
        self.my_analyzer
            .wire_data()
            .map(|wd| wd.wire(brep))
            .unwrap_or_else(Shape::null)
    }

    /// OCCT ShapeFix_Wire::WireAPIMake (lxx L82-85): makes the resulting
    /// Wire (by BRepAPI_MakeWire).
    pub fn wire_api_make(&self) -> Shape {
        self.my_analyzer
            .wire_data()
            .map(|wd| wd.wire_api_make())
            .unwrap_or_else(Shape::null)
    }

    /// OCCT ShapeFix_Wire::Analyzer (lxx L89-92): returns the Analyzer
    /// (working tool).
    pub fn analyzer(&self) -> &ShapeAnalysisWire {
        &self.my_analyzer
    }

    /// OCCT ShapeFix_Wire::WireData (lxx L96-99): returns the working wire.
    pub fn wire_data(&self) -> Option<&crate::shhealing::shape_extend::wire_data::WireData> {
        self.my_analyzer.wire_data()
    }

    /// OCCT ShapeFix_Wire::Face (lxx L103-106): returns the working face
    /// (Analyzer.Face()).
    pub fn face(&self) -> Shape {
        self.my_analyzer.face().clone()
    }

    /// OCCT ShapeFix_Wire::ModifyTopologyMode (lxx L110-113).
    pub fn modify_topology_mode(&mut self) -> &mut bool {
        &mut self.my_topo_mode
    }

    /// OCCT ShapeFix_Wire::ModifyGeometryMode (lxx L117-120).
    pub fn modify_geometry_mode(&mut self) -> &mut bool {
        &mut self.my_geom_mode
    }

    /// OCCT ShapeFix_Wire::ModifyRemoveLoopMode (lxx L124-127).
    pub fn modify_remove_loop_mode(&mut self) -> &mut i32 {
        &mut self.my_remove_loop_mode
    }

    /// OCCT ShapeFix_Wire::ClosedWireMode (lxx L131-134).
    pub fn closed_wire_mode(&mut self) -> &mut bool {
        &mut self.my_closed_mode
    }

    /// OCCT ShapeFix_Wire::PreferencePCurveMode (lxx L138-141).
    pub fn preference_pcurve_mode(&mut self) -> &mut bool {
        &mut self.my_preference2d
    }

    /// OCCT ShapeFix_Wire::FixGapsByRangesMode (lxx L145-148).
    pub fn fix_gaps_by_ranges_mode(&mut self) -> &mut bool {
        &mut self.my_fix_gaps_by_ranges
    }

    /// OCCT ShapeFix_Wire::FixReorderMode (lxx L155-158).
    pub fn fix_reorder_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_reorder_mode
    }

    /// OCCT ShapeFix_Wire::FixSmallMode (lxx L162-165).
    pub fn fix_small_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_small_mode
    }

    /// OCCT ShapeFix_Wire::FixConnectedMode (lxx L169-172).
    pub fn fix_connected_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_connected_mode
    }

    /// OCCT ShapeFix_Wire::FixEdgeCurvesMode (lxx L176-179).
    pub fn fix_edge_curves_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_edge_curves_mode
    }

    /// OCCT ShapeFix_Wire::FixDegeneratedMode (lxx L183-186).
    pub fn fix_degenerated_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_degenerated_mode
    }

    /// OCCT ShapeFix_Wire::FixReversed2dMode (lxx L193-196).
    pub fn fix_reversed2d_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_reversed2d_mode
    }

    /// OCCT ShapeFix_Wire::FixRemovePCurveMode (lxx L200-203).
    pub fn fix_remove_pcurve_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_remove_pcurve_mode
    }

    /// OCCT ShapeFix_Wire::FixRemoveCurve3dMode (lxx L207-210).
    pub fn fix_remove_curve3d_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_remove_curve3d_mode
    }

    /// OCCT ShapeFix_Wire::FixAddPCurveMode (lxx L214-217).
    pub fn fix_add_pcurve_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_add_pcurve_mode
    }

    /// OCCT ShapeFix_Wire::FixAddCurve3dMode (lxx L221-224).
    pub fn fix_add_curve3d_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_add_curve3d_mode
    }

    /// OCCT ShapeFix_Wire::FixSeamMode (lxx L228-231).
    pub fn fix_seam_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_seam_mode
    }

    /// OCCT ShapeFix_Wire::FixShiftedMode (lxx L235-238).
    pub fn fix_shifted_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_shifted_mode
    }

    /// OCCT ShapeFix_Wire::FixSameParameterMode (lxx L242-245).
    pub fn fix_same_parameter_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_same_parameter_mode
    }

    /// OCCT ShapeFix_Wire::FixVertexToleranceMode (lxx L249-252).
    pub fn fix_vertex_tolerance_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_vertex_tolerance_mode
    }

    /// OCCT ShapeFix_Wire::FixLackingMode (lxx L256-259).
    pub fn fix_lacking_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_lacking_mode
    }

    /// OCCT ShapeFix_Wire::FixSelfIntersectionMode (lxx L263-266).
    pub fn fix_self_intersection_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_self_intersection_mode
    }

    /// OCCT ShapeFix_Wire::FixGaps3dMode (lxx L270-273).
    pub fn fix_gaps3d_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_gaps3d_mode
    }

    /// OCCT ShapeFix_Wire::FixGaps2dMode (lxx L277-280).
    pub fn fix_gaps2d_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_gaps2d_mode
    }

    /// OCCT ShapeFix_Wire::FixNotchedEdgesMode (lxx L284-287).
    pub fn fix_notched_edges_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_notched_edges_mode
    }

    /// OCCT ShapeFix_Wire::FixSelfIntersectingEdgeMode (lxx L291-294).
    pub fn fix_self_intersecting_edge_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_self_intersecting_edge_mode
    }

    /// OCCT ShapeFix_Wire::FixIntersectingEdgesMode (lxx L298-301).
    pub fn fix_intersecting_edges_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_intersecting_edges_mode
    }

    /// OCCT ShapeFix_Wire::FixNonAdjacentIntersectingEdgesMode
    /// (lxx L305-308).
    pub fn fix_non_adjacent_intersecting_edges_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_non_adjacent_intersecting_edges_mode
    }

    /// OCCT ShapeFix_Wire::FixTailMode (lxx L312-315).
    pub fn fix_tail_mode(&mut self) -> &mut i32 {
        &mut self.my_fix_tail_mode
    }

    /// OCCT ShapeFix_Wire::StatusReorder (lxx L322-325).
    pub fn status_reorder(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_reorder, status)
    }

    /// OCCT ShapeFix_Wire::StatusSmall (lxx L329-332).
    pub fn status_small(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_small, status)
    }

    /// OCCT ShapeFix_Wire::StatusConnected (lxx L336-339).
    pub fn status_connected(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_connected, status)
    }

    /// OCCT ShapeFix_Wire::StatusEdgeCurves (lxx L343-346).
    pub fn status_edge_curves(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_edge_curves, status)
    }

    /// OCCT ShapeFix_Wire::StatusDegenerated (lxx L350-353).
    pub fn status_degenerated(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_degenerated, status)
    }

    /// OCCT ShapeFix_Wire::StatusLacking (lxx L357-360).
    pub fn status_lacking(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_lacking, status)
    }

    /// OCCT ShapeFix_Wire::StatusSelfIntersection (lxx L364-367).
    pub fn status_self_intersection(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_self_intersection, status)
    }

    /// OCCT ShapeFix_Wire::StatusGaps3d (lxx L371-374).
    pub fn status_gaps3d(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_gaps3d, status)
    }

    /// OCCT ShapeFix_Wire::StatusGaps2d (lxx L378-381).
    pub fn status_gaps2d(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_gaps2d, status)
    }

    /// OCCT ShapeFix_Wire::StatusClosed (lxx L385-388).
    pub fn status_closed(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_closed, status)
    }

    /// OCCT ShapeFix_Wire::StatusNotches (lxx L392-395).
    pub fn status_notches(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_notches, status)
    }

    /// OCCT ShapeFix_Wire::StatusFixTails (lxx L399-402).
    pub fn status_fix_tails(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status_fix_tails, status)
    }

    /// OCCT ShapeFix_Wire::LastFixStatus (lxx L409-412): the status for
    /// low-level methods (common).
    pub fn last_fix_status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_last_fix_status, status)
    }

    /// OCCT ShapeFix_Wire::FixEdgeTool (lxx L416-419): returns the tool for
    /// fixing wires.
    pub fn fix_edge_tool(&mut self) -> &mut ShapeFixEdge {
        &mut self.my_fix_edge
    }

    /// OCCT ShapeFix_Wire::StatusRemovedSegment (lxx L423-426).
    pub fn status_removed_segment(&self) -> bool {
        self.my_status_removed_segment
    }

    // -----------------------------------------------------------------------
    // Internal helpers shared by the submodules.
    // -----------------------------------------------------------------------

    /// OCCT `Context()->Apply(shape)` — the root-context application; a None
    /// context returns the shape unchanged (the callers guard with
    /// `Context().IsNull()` exactly as OCCT does).
    pub(crate) fn context_apply(
        &mut self,
        brep: &mut BRep,
        shape: &Shape,
    ) -> Shape {
        match self.base.my_context.as_mut() {
            Some(ctx) => ctx.apply(brep, shape, rcad_kernel::topods::ShapeType::Shape),
            None => shape.clone(),
        }
    }
}
