//! GAP carriers for the W1-6 batch (ShapeFix package statics +
//! ShapeUpgrade_UnifySameDomain).
//!
//! These stand-ins cover dependencies owned by *other, not-yet-landed*
//! TKShHealing batches.  Each carrier names its OCCT anchor and the owning
//! batch, keeps OCCT's failure/no-op path where the real algorithm is missing,
//! and is replaced wholesale when the owning batch lands:
//!
//! - [`ShapeFixEdgeGap`] — OCCT `ShapeFix_Edge` (ShapeFix_Edge.cxx, 957 LOC)
//!   reduced to `FixSameParameter` + `Status`, the only members consumed by
//!   `ShapeFix::SameParameter` (ShapeFix.cxx L101/L143/L148/L165).  W3 docket
//!   row; the reduced body re-hosts on the landed W1-1
//!   `ShapeAnalysis_Edge::CheckSameParameter`.
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
//!   `BRepLib::SameParameter(edge, tol, exact)` (BRepLib.cxx) reduced to the
//!   sampled-deviation tolerance update (the docket section 4 gap 3, kernel
//!   completion item); consumed by `GlueEdgesWith3DCurves`
//!   (UnifySameDomain.cxx L1748).

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, TShape};
use std::sync::Arc;

use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};
use rcad_kernel::topods::BRepTool;

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Edge — W3 docket row (GAP carrier).
// ---------------------------------------------------------------------------

/// OCCT `ShapeFix_Edge` reduced to the `FixSameParameter` + `Status` surface
/// consumed by `ShapeFix::SameParameter` (ShapeFix.cxx L101, L143, L148,
/// L165).  Replaced wholesale by the W3 1:1 translation.
pub struct ShapeFixEdgeGap {
    /// OCCT myStatus (ShapeFix_Root).
    my_status: i32,
}

impl Default for ShapeFixEdgeGap {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixEdgeGap {
    /// OCCT ShapeFix_Edge() / ShapeFix_Root constructor (ShapeFix_Root.cxx
    /// L31-34): myStatus = ShapeExtend_OK.
    pub fn new() -> Self {
        ShapeFixEdgeGap {
            my_status: encode_status(ShapeExtendStatus::Ok),
        }
    }

    /// OCCT ShapeFix_Edge::Status (ShapeFix_Edge.cxx L940-943) /
    /// ShapeFix_Root::Status — the status bits of the last Fix.
    pub fn status(&self, status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, status)
    }

    /// OCCT ShapeFix_Edge::FixSameParameter(edge, face)
    /// (ShapeFix_Edge.cxx L798-936) — GAP-carrier reduction.
    ///
    /// The OCCT body runs `BRepLib::SameParameter` on a copy of the edge and
    /// compares it against the direct pcurve deviation (choosing the best).
    /// The carrier keeps the observable contract of the walk in
    /// `ShapeFix::SameParameter`: the edge SameParameter flag is set, the
    /// maximal deviation is computed over every pcurve (W1-1
    /// `ShapeAnalysis_Edge::CheckSameParameter`, NbControl = 23), and the
    /// edge/vertex tolerances are raised to it when it exceeds the current
    /// edge tolerance (OCCT L914-929).  The copy-and-compare refinement and
    /// the `TempSameRange` repair are W3 scope.
    pub fn fix_same_parameter(&mut self, brep: &mut BRep, edge: &Shape, face: &Shape) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        // OCCT L804-813: a degenerated edge gets SameRange/SameParameter
        // flags and reports no fix.
        if let Some(ed) = edge.as_edge() {
            if ed.degenerated {
                let mut builder = BRepBuilder::new();
                if !ed.same_range {
                    // OCCT L809: TempSameRange(edge, Precision::PConfusion())
                    // — the flag update only (the range repair is W3 scope).
                    let mut ed_flags = ed.clone();
                    ed_flags.same_range = true;
                    set_edge_data(brep, edge, ed_flags);
                }
                // OCCT L811: B.SameParameter(edge, true).
                builder.set_edge_same_parameter(brep, edge.clone(), true);
                return false;
            }
        }

        // OCCT L820-824: the extremity vertices and the current tolerances.
        let mut sae = ShapeAnalysisEdge::new();
        let v1 = sae.first_vertex(brep, edge);
        let v2 = sae.last_vertex(brep, edge);
        let tol_fv = if v1.is_null() {
            0.0
        } else {
            brep.tolerance(&v1)
        };
        let tol_lv = if v2.is_null() {
            0.0
        } else {
            brep.tolerance(&v2)
        };
        let tol = brep.tolerance(edge);

        // OCCT L826: wasSP = BRep_Tool::SameParameter(edge).
        let was_sp = match edge.data.as_ref() {
            rcad_kernel::topods::TShape::Edge(ed) => ed.same_parameter,
            _ => true,
        };
        // OCCT L872: B.SameParameter(edge, true) before the deviation walk.
        // (The carrier performs no BRepLib copy run, so the flag is set on
        // the edge itself at the same point of the flow.)
        {
            let mut builder = BRepBuilder::new();
            builder.set_edge_same_parameter(brep, edge.clone(), true);
        }

        // OCCT L875-882: compute the deviation on the pcurves (all faces
        // when the input edge was not SameParameter, the given face only
        // otherwise).  NbControl = 23 (the OCCT default, ShapeFix.cxx L231).
        let a_face = if !was_sp { Shape::null() } else { face.clone() };
        let mut maxdev = 0.0;
        sae.check_same_parameter_face(brep, edge, &a_face, &mut maxdev, 23);
        if sae.status(ShapeExtendStatus::Fail2) {
            // OCCT L883-886: FAIL2 of the check maps to FAIL1 here.
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
        }

        // OCCT L914-922: restore the vertex tolerances to at least maxdev.
        let mut builder = BRepBuilder::new();
        if !v1.is_null() {
            update_vertex_tolerance_max(brep, &mut builder, &v1, maxdev.max(tol_fv));
        }
        if !v2.is_null() {
            update_vertex_tolerance_max(brep, &mut builder, &v2, maxdev.max(tol_lv));
        }

        // OCCT L924-929: raise the edge tolerance when the deviation
        // exceeds it (B.UpdateEdge(edge, maxdev) + FixVertexTolerance —
        // the vertex floors above cover the FixVertexTolerance step).
        if maxdev > tol {
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
            set_edge_tolerance_value(brep, edge, maxdev);
        }

        // OCCT L931-934: !wasSP && !SP -> DONE2 (the BRepLib branch did not
        // reach a SameParameter result; the carrier has no copy run).
        if !was_sp {
            self.my_status |= encode_status(ShapeExtendStatus::Done2);
        }
        self.status(ShapeExtendStatus::Done)
    }
}

/// OCCT `BRep_Builder::UpdateVertex(V, Tol)` — the max-with-existing
/// tolerance update (BRep_Builder.cxx L1194-1216).
fn update_vertex_tolerance_max(brep: &mut BRep, builder: &mut BRepBuilder, v: &Shape, tol: f64) {
    builder.update_vertex_tolerance(brep, v.clone(), tol);
}

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

/// In-place TEdgeData replace helper for the carrier (the
/// `BRep_Builder::SameRange` flag writes through the shared handle).
fn set_edge_data(brep: &mut BRep, edge: &Shape, ed: rcad_kernel::topods::TEdgeData) {
    let ptr = Arc::as_ptr(&brep.tshapes[edge.index]) as *mut TShape;
    let ts = unsafe { &mut *ptr };
    *ts = TShape::Edge(ed);
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
