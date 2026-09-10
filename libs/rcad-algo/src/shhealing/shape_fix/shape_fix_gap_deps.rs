//! GAP carriers for the W1-6 batch (ShapeFix package statics +
//! ShapeUpgrade_UnifySameDomain).
//!
//! These stand-ins cover dependencies owned by *other, not-yet-landed*
//! TKShHealing batches.  Each carrier names its OCCT anchor and the owning
//! batch, keeps OCCT's failure/no-op path where the real algorithm is missing,
//! and is replaced wholesale when the owning batch lands:
//!
//! - [`ShapeFixFaceGap`] — OCCT `ShapeFix_Face` (ShapeFix_Face.cxx, 3,259
//!   LOC) reduced to the accessor/Perform surface consumed by
//!   `ShapeUpgrade_UnifySameDomain::UnifyEdges` (cxx L4404-4421).  W3 docket
//!   row; `Perform` keeps OCCT's "nothing fixed" path (the face is returned
//!   unchanged).
//! - [`ShapeFixShellGap`] — OCCT `ShapeFix_Shell` (ShapeFix_Shell.cxx, 1,727
//!   LOC) reduced to `FixFaceOrientation`/`Shell` consumed by UnifyEdges
//!   (cxx L4433-4435).  W3 docket row; keeps OCCT's "shell unchanged" path.
//! - [`ShapeFixWireGap`] — OCCT `ShapeFix_Wire` (ShapeFix_Wire.cxx + _1,
//!   6,483 LOC) reduced to the mode-flag accessors consumed by the
//!   `SetFixWireModes` static (UnifySameDomain.cxx L3133-3144).  W3 docket
//!   row; the flags are stored but drive no behavior.
//! - [`ShapeFixShapeGap`] — OCCT `ShapeFix_Shape` (ShapeFix_Shape.cxx, 358
//!   LOC) reduced to the surface consumed by `ShapeFix::RemoveSmallEdges`
//!   (ShapeFix.cxx L291-308).  W3 docket row; `Perform` keeps OCCT's
//!   "nothing fixed" path.
//! - [`brep_lib_same_parameter_edge`] — OCCT
//!   `BRepLib::SameParameter(edge, tol)` (BRepLib.cxx, BRepLib.hxx L161)
//!   reduced to the sampled-deviation tolerance update (the docket section 4
//!   gap 3, kernel completion item); consumed by `GlueEdgesWith3DCurves`
//!   (UnifySameDomain.cxx L1748) and, since the W3 tranche 1, by the 1:1
//!   `ShapeFix_Edge::FixSameParameter` (ShapeFix_Edge.cxx L850, its true
//!   OCCT consumer).
//!
//! Retired by the W3 tranche 1 (ShapeFix_Edge landed 1:1 in
//! `shape_fix/edge.rs`): the former `ShapeFixEdgeGap` carrier (the
//! FixSameParameter + Status reduction consumed by `ShapeFix::SameParameter`)
//! — deleted, Rule 4.

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, TShape};
use std::sync::Arc;

use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face — W3 docket row (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `BRep_Builder::UpdateEdge(E, Tol)` (BRep_Builder.cxx L578-600) —
/// the max-with-existing edge tolerance update.
fn set_edge_tolerance_value(brep: &mut BRep, edge: &Shape, tol: f64) {
    // SAFETY: same in-place pattern as BRep::edge_mut_inplace (the healing
    // chain mutates the TShape through the shared handle).
    let ptr = Arc::as_ptr(&brep.tshapes[edge.index]) as *mut TShape;
    let ts = unsafe { &mut *ptr };
    if let TShape::Edge(ed) = ts {
        ed.tolerance = ed.tolerance.max(tol);
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face — W3 docket row (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `ShapeFix_Face` reduced to the accessor/Perform surface consumed by
/// `ShapeUpgrade_UnifySameDomain::UnifyEdges` (UnifySameDomain.cxx
/// L4404-4421).  `Perform` keeps OCCT's "nothing fixed" path: the face is
/// returned unchanged.  Replaced wholesale by the W3 1:1 translation.
pub struct ShapeFixFaceGap {
    my_face: Shape,
    my_precision: f64,
    my_min_tol: f64,
    my_max_tol: f64,
    my_fix_orientation_mode: bool,
    my_fix_missing_seam_mode: bool,
    my_fix_small_area_wire_mode: bool,
    my_fix_add_natural_bound_mode: bool,
    my_fix_intersecting_wires_mode: bool,
    my_fix_loop_wires_mode: bool,
    my_fix_split_face_mode: bool,
    my_fix_periodic_degenerated_mode: bool,
    my_fix_wire: ShapeFixWireGap,
}

impl Default for ShapeFixFaceGap {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixFaceGap {
    /// OCCT ShapeFix_Face() (ShapeFix_Face.cxx L50-62).
    pub fn new() -> Self {
        ShapeFixFaceGap {
            my_face: Shape::null(),
            my_precision: 0.0,
            my_min_tol: 0.0,
            my_max_tol: 0.0,
            my_fix_orientation_mode: true,
            my_fix_missing_seam_mode: true,
            my_fix_small_area_wire_mode: false,
            my_fix_add_natural_bound_mode: true,
            my_fix_intersecting_wires_mode: false,
            my_fix_loop_wires_mode: false,
            my_fix_split_face_mode: true,
            my_fix_periodic_degenerated_mode: true,
            my_fix_wire: ShapeFixWireGap::new(),
        }
    }

    /// OCCT ShapeFix_Face::Init(f) (ShapeFix_Face.cxx L86-107) — the
    /// (face, precision) form used at UnifySameDomain.cxx L4404.
    pub fn with_face(the_face: &Shape) -> Self {
        let mut s = Self::new();
        s.my_face = the_face.clone();
        s
    }

    /// OCCT ShapeFix_ShapeShell getters: Face().
    pub fn face(&self) -> Shape {
        self.my_face.clone()
    }

    /// OCCT ShapeFix_Root::SetPrecision.
    pub fn set_precision(&mut self, prec: f64) {
        self.my_precision = prec;
    }

    /// OCCT ShapeFix_Root::SetMinTolerance.
    pub fn set_min_tolerance(&mut self, val: f64) {
        self.my_min_tol = val;
    }

    /// OCCT ShapeFix_Root::SetMaxTolerance.
    pub fn set_max_tolerance(&mut self, val: f64) {
        self.my_max_tol = val;
    }

    /// OCCT FixOrientationMode flag.
    pub fn fix_orientation_mode(&mut self) -> &mut bool {
        &mut self.my_fix_orientation_mode
    }

    /// OCCT FixMissingSeamMode flag.
    pub fn fix_missing_seam_mode(&mut self) -> &mut bool {
        &mut self.my_fix_missing_seam_mode
    }

    /// OCCT FixSmallAreaWireMode flag.
    pub fn fix_small_area_wire_mode(&mut self) -> &mut bool {
        &mut self.my_fix_small_area_wire_mode
    }

    /// OCCT FixAddNaturalBoundMode flag.
    pub fn fix_add_natural_bound_mode(&mut self) -> &mut bool {
        &mut self.my_fix_add_natural_bound_mode
    }

    /// OCCT FixIntersectingWiresMode flag.
    pub fn fix_intersecting_wires_mode(&mut self) -> &mut bool {
        &mut self.my_fix_intersecting_wires_mode
    }

    /// OCCT FixLoopWiresMode flag.
    pub fn fix_loop_wires_mode(&mut self) -> &mut bool {
        &mut self.my_fix_loop_wires_mode
    }

    /// OCCT FixSplitFaceMode flag.
    pub fn fix_split_face_mode(&mut self) -> &mut bool {
        &mut self.my_fix_split_face_mode
    }

    /// OCCT FixPeriodicDegeneratedMode flag.
    pub fn fix_periodic_degenerated_mode(&mut self) -> &mut bool {
        &mut self.my_fix_periodic_degenerated_mode
    }

    /// OCCT FixWireTool() — the ShapeFix_Wire tool handle.
    pub fn fix_wire_tool(&mut self) -> &mut ShapeFixWireGap {
        &mut self.my_fix_wire
    }

    /// OCCT ShapeFix_Face::Perform (ShapeFix_Face.cxx L446-470) — GAP:
    /// returns the "nothing fixed" path (the face passes through
    /// unchanged); the orientation/wire-fix machinery is W3 scope.
    pub fn perform(&mut self, _brep: &mut BRep) -> bool {
        false
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shell — W3 docket row (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `ShapeFix_Shell` reduced to `FixFaceOrientation`/`Shell` consumed by
/// UnifySameDomain.cxx L4433-4435.  Keeps OCCT's "shell unchanged" path.
pub struct ShapeFixShellGap {
    my_shell: Shape,
}

impl Default for ShapeFixShellGap {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixShellGap {
    /// OCCT ShapeFix_Shell() (ShapeFix_Shell.cxx L37-40).
    pub fn new() -> Self {
        ShapeFixShellGap {
            my_shell: Shape::null(),
        }
    }

    /// OCCT ShapeFix_Shell::Shell() — the resulting shell.
    pub fn shell(&self) -> Shape {
        self.my_shell.clone()
    }

    /// OCCT ShapeFix_Shell::FixFaceOrientation(shell)
    /// (ShapeFix_Shell.cxx L120-256) — GAP: keeps the "shell unchanged"
    /// path (the face-orientation machinery is W3 scope).
    pub fn fix_face_orientation(&mut self, the_shell: &Shape) -> bool {
        self.my_shell = the_shell.clone();
        false
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Wire — W3 docket row (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `ShapeFix_Wire` reduced to the mode-flag accessors consumed by the
/// `SetFixWireModes` static (UnifySameDomain.cxx L3133-3144).  The flags are
/// stored but drive no behavior; replaced wholesale by the W3 translation.
pub struct ShapeFixWireGap {
    my_fix_self_intersection_mode: bool,
    my_fix_non_adjacent_intersecting_edges_mode: bool,
    my_fix_lacking_mode: bool,
    my_fix_notched_edges_mode: bool,
    my_modify_topology_mode: bool,
    my_modify_remove_loop_mode: bool,
    my_fix_gaps_by_ranges_mode: bool,
    my_fix_small_mode: bool,
    my_fix_connected_mode: bool,
    my_fix_edge_curves_mode: bool,
    my_fix_degenerated_mode: bool,
}

impl Default for ShapeFixWireGap {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixWireGap {
    /// OCCT ShapeFix_Wire::ShapeFix_Wire() (ShapeFix_Wire.cxx L60-67).
    pub fn new() -> Self {
        ShapeFixWireGap {
            my_fix_self_intersection_mode: true,
            my_fix_non_adjacent_intersecting_edges_mode: false,
            my_fix_lacking_mode: true,
            my_fix_notched_edges_mode: true,
            my_modify_topology_mode: false,
            my_modify_remove_loop_mode: false,
            my_fix_gaps_by_ranges_mode: false,
            my_fix_small_mode: true,
            my_fix_connected_mode: true,
            my_fix_edge_curves_mode: true,
            my_fix_degenerated_mode: true,
        }
    }

    /// OCCT FixSelfIntersectionMode flag.
    pub fn fix_self_intersection_mode(&mut self) -> &mut bool {
        &mut self.my_fix_self_intersection_mode
    }

    /// OCCT FixNonAdjacentIntersectingEdgesMode flag.
    pub fn fix_non_adjacent_intersecting_edges_mode(&mut self) -> &mut bool {
        &mut self.my_fix_non_adjacent_intersecting_edges_mode
    }

    /// OCCT FixLackingMode flag.
    pub fn fix_lacking_mode(&mut self) -> &mut bool {
        &mut self.my_fix_lacking_mode
    }

    /// OCCT FixNotchedEdgesMode flag.
    pub fn fix_notched_edges_mode(&mut self) -> &mut bool {
        &mut self.my_fix_notched_edges_mode
    }

    /// OCCT ModifyTopologyMode flag.
    pub fn modify_topology_mode(&mut self) -> &mut bool {
        &mut self.my_modify_topology_mode
    }

    /// OCCT ModifyRemoveLoopMode flag.
    pub fn modify_remove_loop_mode(&mut self) -> &mut bool {
        &mut self.my_modify_remove_loop_mode
    }

    /// OCCT FixGapsByRangesMode flag.
    pub fn fix_gaps_by_ranges_mode(&mut self) -> &mut bool {
        &mut self.my_fix_gaps_by_ranges_mode
    }

    /// OCCT FixSmallMode flag.
    pub fn fix_small_mode(&mut self) -> &mut bool {
        &mut self.my_fix_small_mode
    }

    /// OCCT FixConnectedMode flag.
    pub fn fix_connected_mode(&mut self) -> &mut bool {
        &mut self.my_fix_connected_mode
    }

    /// OCCT FixEdgeCurvesMode flag.
    pub fn fix_edge_curves_mode(&mut self) -> &mut bool {
        &mut self.my_fix_edge_curves_mode
    }

    /// OCCT FixDegeneratedMode flag.
    pub fn fix_degenerated_mode(&mut self) -> &mut bool {
        &mut self.my_fix_degenerated_mode
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Shape — W3 docket row (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `ShapeFix_Shape` reduced to the surface consumed by
/// `ShapeFix::RemoveSmallEdges` (ShapeFix.cxx L291-308).  `Perform` keeps
/// OCCT's "nothing fixed" path; replaced wholesale by the W3 translation.
pub struct ShapeFixShapeGap {
    my_shape: Shape,
    my_precision: f64,
    my_context: Option<()>,
    my_face_tool: ShapeFixFaceGap,
    my_wire_tool: ShapeFixWireGap,
}

impl Default for ShapeFixShapeGap {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixShapeGap {
    /// OCCT ShapeFix_Shape() (ShapeFix_Shape.cxx L40-45).
    pub fn new() -> Self {
        ShapeFixShapeGap {
            my_shape: Shape::null(),
            my_precision: 0.0,
            my_context: None,
            my_face_tool: ShapeFixFaceGap::new(),
            my_wire_tool: ShapeFixWireGap::new(),
        }
    }

    /// OCCT ShapeFix_Shape::Init(S) (ShapeFix_Shape.cxx L52-112) — GAP: the
    /// shape is stored and returned as is (the sub-shape fix chains are W3
    /// scope).
    pub fn init(&mut self, s: &Shape) {
        self.my_shape = s.clone();
    }

    /// OCCT ShapeFix_Root::SetPrecision.
    pub fn set_precision(&mut self, prec: f64) {
        self.my_precision = prec;
    }

    /// OCCT FixFaceTool() — the ShapeFix_Face tool handle.
    pub fn fix_face_tool(&mut self) -> &mut ShapeFixFaceGap {
        &mut self.my_face_tool
    }

    /// OCCT FixWireTool() — the ShapeFix_Wire tool handle.
    pub fn fix_wire_tool(&mut self) -> &mut ShapeFixWireGap {
        &mut self.my_wire_tool
    }

    /// OCCT ShapeFix_Shape::Perform — GAP: the "nothing fixed" path.
    pub fn perform(&mut self) -> bool {
        false
    }

    /// OCCT ShapeFix_Shape::Shape() — the resulting shape.
    pub fn shape(&self) -> Shape {
        self.my_shape.clone()
    }

    /// OCCT ShapeFix_Shape::Context() — the reshape context (None: the
    /// carrier performs no replacements).
    pub fn context(&self) -> Option<()> {
        self.my_context
    }
}

// ---------------------------------------------------------------------------
// OCCT BRepLib::SameParameter(edge) — docket section 4 gap 3 (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `BRepLib::SameParameter(const TopoDS_Edge&, const Standard_Real tol,
/// const Standard_Boolean exact)` (BRepLib.cxx L375-457) — GAP-carrier
/// reduction (the docket section 4 gap 3: the kernel completion item).
///
/// The OCCT body computes the maximal deviation between the 3D curve and
/// every pcurve of the edge and raises the edge (and vertex) tolerances to
/// it, setting the SameParameter flag.  The carrier re-hosts exactly that
/// contract on the W1-1 `ShapeAnalysis_Edge::CheckSameParameter` sampler
/// (NbControl = 23, the OCCT BRepLib default); the `exact` branch that
/// enforces same parametrization by rebuilding the pcurves is kernel scope.
/// Consumed by `GlueEdgesWith3DCurves` (UnifySameDomain.cxx L1748).
pub fn brep_lib_same_parameter_edge(brep: &mut BRep, edge: &Shape, tol: f64, _exact: bool) {
    // OCCT L380-382: if (BRep_Tool::SameParameter(E) && BRep_Tool::SameRange(E)
    // && tol <= BRep_Tool::Tolerance(E)) return — nothing to do.
    let (sp, sr, etol) = match edge.data.as_ref() {
        rcad_kernel::topods::TShape::Edge(ed) => (ed.same_parameter, ed.same_range, ed.tolerance),
        _ => return,
    };
    if sp && sr && tol <= etol {
        return;
    }

    // OCCT L389-424: compute the maximal deviation over the pcurves.
    let mut sae = ShapeAnalysisEdge::new();
    let mut maxdev = 0.0;
    sae.check_same_parameter(brep, edge, &mut maxdev, 23);
    let new_tol = if maxdev > tol { maxdev } else { tol };

    // OCCT L426-457: update the edge tolerance and set the flags.
    set_edge_tolerance_value(brep, edge, new_tol);
    let mut builder = BRepBuilder::new();
    builder.set_edge_same_parameter(brep, edge.clone(), true);
    builder.set_edge_same_range(brep, edge.clone(), true);
}
