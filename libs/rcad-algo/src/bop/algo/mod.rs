//! OCCT BOPAlgo ?Boolean Operation Algorithms.
//!
//! | Module               | OCCT class                  | Description                     |
//! |----------------------|-----------------------------|---------------------------------|
//! | alerts               | BOPAlgo_Alerts              | Alert types                     |
//! | check_result         | BOPAlgo_CheckResult         | Argument analysis result        |
//! | checker_si           | BOPAlgo_CheckerSI           | Self-interference check         |
//! | glue_enum            | BOPAlgo_GlueEnum            | Glue mode enumeration           |
//! | operation            | BOPAlgo_Operation           | Operation type enumeration      |
//! | pave_filler          | BOPAlgo_PaveFiller          | Intersection computation        |
//! | builder              | BOPAlgo_Builder             | Result construction             |
//! | builder_area         | BOPAlgo_BuilderArea         | Area building                   |
//! | builder_face         | BOPAlgo_BuilderFace         | Face splitting                  |
//! | builder_solid        | BOPAlgo_BuilderSolid        | Solid building                  |
//! | shell_splitter       | BOPAlgo_ShellSplitter       | Shell partitioning              |
//! | section_attribute    | BOPAlgo_SectionAttribute    | Section parameters              |

pub mod pave_filler;
pub mod pave_filler_make_blocks;
pub mod occt_map;
pub mod builder;
pub mod builder_area;
pub mod builder_face;
pub mod builder_solid;
pub mod shell_splitter;
pub mod wire_splitter;
pub mod argument_analyzer;
pub mod check_result;
pub mod checker_si;
pub mod section_attribute;

pub use check_result::{CheckResult, CheckStatus};

// Re-export shared types
pub use builder::{Builder, BooleanError};

/// OCCT TopoDS_Iterator cumLoc semantics: the effective Location of an edge's
/// vertex is `edge.Location * vertex.Location` (TopoDS_Iterator.cxx L76-78,
/// TopoDS_Shape::Move L193-201 — the parent Location is pre-multiplied).
/// rcad stores Locations as indices into the DS table; the composed transform
/// is looked up in the same table. Identity fast paths (the common case —
/// every non-folded shape carries location 0) return the original index
/// without touching the table, so no registration is ever needed for them.
pub fn compose_edge_vertex_location(
    edge_loc: u32,
    vertex_loc: u32,
    locations: &[glam::DAffine3],
) -> u32 {
    if edge_loc == 0 {
        return vertex_loc;
    }
    if vertex_loc == 0 {
        return edge_loc;
    }
    // Both non-identity: look the composed transform up in the table. The
    // fused inputs never produce a composed transform absent from the table
    // (single-fold shapes compose to one of the two operands); fall back to
    // the edge location if it is ever missing.
    let edge_tr = locations.get(edge_loc as usize).copied().unwrap_or(glam::DAffine3::IDENTITY);
    let vertex_tr = locations
        .get(vertex_loc as usize)
        .copied()
        .unwrap_or(glam::DAffine3::IDENTITY);
    let composed = edge_tr * vertex_tr;
    locations
        .iter()
        .position(|l| *l == composed)
        .map(|i| i as u32)
        .unwrap_or(edge_loc)
}

/// OCCT BRep_Builder::UpdateEdge / BRep_Tool::CurveOnSurface
/// (BRep_Builder.cxx L692, BRep_Tool.cxx L345): the pcurve of an edge on a
/// face is keyed by `(face TShape, aLoc)` with
/// `aLoc = L.Predivided(E.Location())` — the face location divided by the
/// edge's location. A located edge (e.g. the top edge of a prism, sharing its
/// TShape with the base edge) therefore has its OWN pcurve key, distinct from
/// the base edge's.
///
/// rcad stores Locations as indices into the DS table. The composed transform
/// face_loc * edge_loc^-1 may not be present in the table (only the forward
/// folds are registered); fall back to the edge location index when the face
/// has no location, and to the face location otherwise — the key still
/// separates located edges from their base copies, which is the OCCT
/// semantics that matters for the shared-TShape case.
pub fn compose_face_edge_pcurve_location(
    face_loc: u32,
    edge_loc: u32,
    locations: &[glam::DAffine3],
) -> u32 {
    let tr = |idx: u32| -> glam::DAffine3 {
        // DS location table: slot 0 stores identity, real transforms start
        // at index 1 (DS::new / brep_top_shapes_with_locations).
        locations
            .get(idx as usize)
            .copied()
            .unwrap_or(glam::DAffine3::IDENTITY)
    };
    let composed = tr(face_loc) * tr(edge_loc).inverse();
    rcad_kernel::topo::topods::pcurve_location_id(&composed)
}

// ===
// BOPAlgo_Report ?collects alerts during algorithm execution
// ===
#[derive(Debug, Clone, Default)]
pub struct Report {
    alerts: Vec<Alert>,
    has_errors: bool,
}
impl Report {
    pub fn new() -> Self { Report { alerts: Vec::new(), has_errors: false } }
    pub fn add_alert(&mut self, a: Alert) { self.alerts.push(a); }
    pub fn add_error(&mut self, a: Alert) { self.alerts.push(a); self.has_errors = true; }
    pub fn add_warning(&mut self, a: Alert) { self.alerts.push(a); }
    pub fn has_errors(&self) -> bool { self.has_errors }
    pub fn errors(&self) -> &[Alert] { &self.alerts }
    pub fn clear(&mut self) { self.alerts.clear(); self.has_errors = false; }
    /// OCCT Message_Report::Merge — append the alerts of another report.
    pub fn merge(&mut self, other: Report) {
        self.has_errors |= other.has_errors;
        self.alerts.extend(other.alerts);
    }
}

/// OCCT BOPAlgo_GlueEnum — glue mode for coincident geometry.
///
/// Three modes:
/// - GlueOff: default, full intersection computation
/// - GlueShift: for partially coincident shapes (skips FACE/FACE intersections)
/// - GlueFull: for fully coincident shapes (skips VERTEX/FACE, EDGE/FACE, FACE/FACE)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlueEnum {
    GlueOff = 0,
    GlueShift = 1,
    GlueFull = 2,
}

impl Default for GlueEnum {
    fn default() -> Self { GlueEnum::GlueOff }
}

// OCCT BOPAlgo_Operation — boolean operation type enumeration.
//
// OCCT BOPAlgo_Operation.hxx
// Includes: COMMON, FUSE, CUT, CUT21, SECTION, UNKNOWN.

/// OCCT BOPAlgo_Operation — type of boolean operation.
/// OCCT BOPAlgo_Operation.hxx: COMMON, FUSE, CUT, CUT21, SECTION, UNKNOWN.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    /// Boolean union (OCCT BOPAlgo_FUSE).
    Union,
    /// Boolean intersection (OCCT BOPAlgo_COMMON).
    Intersection,
    /// Boolean subtraction a - b (OCCT BOPAlgo_CUT).
    Cut,
    /// Boolean subtraction b - a (OCCT BOPAlgo_CUT21).
    Cut21,
    /// Boolean section / intersection curves (OCCT BOPAlgo_SECTION).
    Section,
    /// Undefined operation (OCCT BOPAlgo_UNKNOWN).
    Unknown,
}

// OCCT BOPAlgo_Alerts — diagnostic alert types for Boolean Operations.
//
// OCCT BOPAlgo_Alerts.hxx

/// OCCT BOPAlgo_Alert — diagnostic alerts for Boolean Operation algorithms.
#[derive(Debug, Clone)]
pub enum Alert {
    TooSmallRange(usize, f64),
    BuildingPCurveFailed(usize, usize),
    SolidBuilderUnusedFaces(Vec<usize>),
    EdgeWithoutCurve(usize),
    TooFewArguments,
    NoFiller,
    BOPNotAllowed,
    BOPNotSet,
    EmptyShape,
    IntersectionFailed(usize, usize),
    /// BOPAlgo_AlertSelfInterferingShape — vertex/edge from same argument
    /// with unexpected interference.
    SelfInterferingShape(Vec<usize>),
    /// BOPAlgo_AlertAcquiredSelfIntersection — acquired self-intersection
    /// detected during CheckSelfInterference.
    AcquiredSelfIntersection(Vec<usize>),
    /// BOPAlgo_AlertTooSmallEdge — edge range too small for splitting.
    TooSmallEdge(usize),
    /// BOPAlgo_AlertNotSplittableEdge (stub).
    NotSplittableEdge(usize),
    /// BOPAlgo_AlertBadPositioning (stub).
    BadPositioning(Vec<usize>),
    /// BOPAlgo_AlertPostTreatFF — error in the MakeBlocks post-treatment.
    PostTreatFF,
    /// BOPAlgo_AlertSolidBuilderFailed — solid builder failed to build solids.
    SolidBuilderFailed,
    /// BOPAlgo_AlertUnableToMakeClosedEdgeOnFace — a seam (closed-surface) edge
    /// could not be split on the face. Carries the (face, split edge) shapes.
    UnableToMakeClosedEdgeOnFace(Vec<rcad_kernel::topo_shape::Shape>),
}
