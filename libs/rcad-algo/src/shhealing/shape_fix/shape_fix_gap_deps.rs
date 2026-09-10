//! GAP carriers for the ShapeFix batch (ShapeFix package statics +
//! ShapeUpgrade_UnifySameDomain).
//!
//! These stand-ins cover dependencies owned by *other, not-yet-landed*
//! TKShHealing batches.  Each carrier names its OCCT anchor and the owning
//! batch, keeps OCCT's failure/no-op path where the real algorithm is missing,
//! and is replaced wholesale when the owning batch lands:
//!
//! - [`ShapeFixSolidGap`] — OCCT `ShapeFix_Solid` (ShapeFix_Solid.cxx, 749
//!   LOC) reduced to the surface consumed by `ShapeFix_Shape`
//!   (ShapeFix_Shape.cxx L49/L62/L165-166/L173/L323-351 and the lxx tool
//!   chain L24-55).  W3 docket row; `Perform` keeps OCCT's "nothing fixed"
//!   path; the real `ShapeFixShell` is embedded so the
//!   FixShellTool/FixFaceTool/FixWireTool/FixEdgeTool chain stays real.
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
//!
//! Retired by the W3 tranche 3 (ShapeFix_Face / ShapeFix_Shell landed 1:1 in
//! `shape_fix/face_a.rs` + `shape_fix/shell.rs`): the former
//! `ShapeFixFaceGap` (the accessor/Perform reduction consumed by
//! UnifySameDomain.cxx L4404-4421) and `ShapeFixShellGap` (the
//! FixFaceOrientation/Shell reduction consumed by UnifySameDomain.cxx
//! L4433-4445) carriers — deleted, Rule 4; every consumer now calls the real
//! classes.
//!
//! Retired by the W3 tranche 4 (ShapeFix_Shape landed 1:1 in
//! `shape_fix/shape_fix_shape.rs`): the former `ShapeFixShapeGap` carrier
//! (the Init/Perform/Shape reduction consumed by `ShapeFix::RemoveSmallEdges`,
//! ShapeFix.cxx L291-308) — deleted, Rule 4; the static now drives the real
//! class.

use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, TShape};
use std::sync::Arc;

use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::root::{MsgRegistratorHandle, ShapeFixRoot};
use crate::shhealing::shape_fix::shape_fix::MessageProgressRange;
use crate::shhealing::shape_fix::shell::ShapeFixShell;

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
// OCCT ShapeFix_Solid — W3 docket row (GAP carrier).
// ---------------------------------------------------------------------------

// Retired by the W3 tranche 4 (ShapeFix_Shape landed 1:1 in
// `shape_fix/shape_fix_shape.rs`): the former `ShapeFixShapeGap` carrier
// (the Init/Perform/Shape reduction consumed by `ShapeFix::RemoveSmallEdges`,
// ShapeFix.cxx L291-308) — deleted, Rule 4; the static drives the real class.

/// OCCT `ShapeFix_Solid` — GAP carrier for the not-yet-landed W3 docket row
/// (ShapeFix_Solid.cxx, 749 LOC).  The carrier hosts exactly the surface
/// consumed by `ShapeFix_Shape` (the constructor `new ShapeFix_Solid`,
/// ShapeFix_Shape.cxx L49/L62; the SOLID case Init/SetContext/Perform calls,
/// ShapeFix_Shape.cxx L165-171; the Set* forwarders, ShapeFix_Shape.cxx
/// L323-351): the empty constructor (cxx L52-60), Init (cxx L75-84), the
/// inline `FixShellTool` (hxx L69), and the Root setters that forward to the
/// shell tool (cxx L721-748).  The real `ShapeFixShell` (the W3 tranche 3
/// class) is embedded, so the FixShellTool/FixFaceTool/FixWireTool/
/// FixEdgeTool chain of ShapeFix_Shape stays real.  `Perform` (cxx L460-645)
/// keeps OCCT's "nothing fixed" path — returns false and leaves the shape
/// unchanged; replaced wholesale when the ShapeFix_Solid batch lands.
// The myStatus/myFixShellMode/myFixShellOrientationMode/myCreateOpenSolidMode
// fields are stored per the OCCT hxx L133-136 but not consumed by the hosted
// surface (the carrier does not translate Status()/Perform()).
#[allow(dead_code)]
pub struct ShapeFixSolidGap {
    /// The OCCT ShapeFix_Root base subobject.
    pub base: ShapeFixRoot,
    /// OCCT mySolid (hxx L131).
    pub(crate) my_solid: Shape,
    /// OCCT myFixShell (hxx L132) — the real W3 tranche 3 class.
    pub(crate) my_fix_shell: ShapeFixShell,
    /// OCCT myStatus (hxx L133).
    pub(crate) my_status: i32,
    /// OCCT myFixShellMode (hxx L134).
    pub(crate) my_fix_shell_mode: i32,
    /// OCCT myFixShellOrientationMode (hxx L135).
    pub(crate) my_fix_shell_orientation_mode: i32,
    /// OCCT myCreateOpenSolidMode (hxx L136).
    pub(crate) my_create_open_solid_mode: bool,
}

impl Default for ShapeFixSolidGap {
    fn default() -> Self {
        Self::new()
    }
}

impl ShapeFixSolidGap {
    /// OCCT ShapeFix_Solid::ShapeFix_Solid() (cxx L52-60): empty constructor.
    pub fn new() -> Self {
        ShapeFixSolidGap {
            base: ShapeFixRoot::new(),
            // OCCT L131: mySolid — the default-constructed null handle.
            my_solid: Shape::null(),
            // OCCT L57: myFixShell = new ShapeFix_Shell.
            my_fix_shell: ShapeFixShell::new(),
            // OCCT L53: myStatus = EncodeStatus(ShapeExtend_OK).
            my_status: encode_status(ShapeExtendStatus::Ok),
            // OCCT L54-56.
            my_fix_shell_mode: -1,
            my_fix_shell_orientation_mode: -1,
            my_create_open_solid_mode: false,
        }
    }

    /// OCCT ShapeFix_Solid::Init (cxx L75-84): initializes by a solid
    /// (mySolid = solid; myShape = solid).
    pub fn init(&mut self, _brep: &mut BRep, solid: &Shape) {
        // OCCT L76: mySolid = solid.
        self.my_solid = solid.clone();
        // OCCT L83: myShape = solid.
        self.base.my_shape = solid.clone();
    }

    /// OCCT ShapeFix_Solid::FixShellTool (hxx L69, inline): returns the tool
    /// for fixing shells (the real W3 tranche 3 class).
    pub fn fix_shell_tool(&mut self) -> &mut ShapeFixShell {
        &mut self.my_fix_shell
    }

    /// OCCT ShapeFix_Solid::Perform (cxx L460-645) — GAP: the "nothing
    /// fixed" path (the shell-per-shell `ShapeFix_Shell::Perform` loop and
    /// the shell-orientation stage are the ShapeFix_Solid W3 row scope); the
    /// shape passes through unchanged and no DONE status is recorded.
    pub fn perform(&mut self, _brep: &mut BRep, _the_progress: MessageProgressRange) -> bool {
        false
    }

    /// OCCT ShapeFix_Solid::SetMsgRegistrator (cxx L721-725).
    pub fn set_msg_registrator(&mut self, msgreg: MsgRegistratorHandle) {
        // OCCT L723: ShapeFix_Root::SetMsgRegistrator(msgreg).
        self.base.set_msg_registrator(msgreg.clone());
        // OCCT L724: myFixShell->SetMsgRegistrator(msgreg).
        self.my_fix_shell.set_msg_registrator(msgreg);
    }

    /// OCCT ShapeFix_Solid::SetPrecision (cxx L729-733).
    pub fn set_precision(&mut self, preci: f64) {
        // OCCT L731.
        self.base.set_precision(preci);
        // OCCT L732: myFixShell->SetPrecision(preci).
        self.my_fix_shell.set_precision(preci);
    }

    /// OCCT ShapeFix_Solid::SetMinTolerance (cxx L737-741).
    pub fn set_min_tolerance(&mut self, mintol: f64) {
        // OCCT L739.
        self.base.set_min_tolerance(mintol);
        // OCCT L740: myFixShell->SetMinTolerance(mintol).
        self.my_fix_shell.set_min_tolerance(mintol);
    }

    /// OCCT ShapeFix_Solid::SetMaxTolerance (cxx L745-748).
    pub fn set_max_tolerance(&mut self, maxtol: f64) {
        // OCCT L747.
        self.base.set_max_tolerance(maxtol);
        // OCCT L748: myFixShell->SetMaxTolerance(maxtol).
        self.my_fix_shell.set_max_tolerance(maxtol);
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
