//! OCCT BOPAlgo ?Boolean Operation Algorithms.
//!
//! | Module          | OCCT class              | Description                     |
//! |-----------------|-------------------------|---------------------------------|
//! | pave_filler     | BOPAlgo_PaveFiller      | Intersection computation        |
//! | builder         | BOPAlgo_Builder         | Result construction             |
//! | builder_face    | BOPAlgo_BuilderFace     | Face splitting                  |
//! | builder_solid   | BOPAlgo_BuilderSolid    | Solid building                  |
//! | shell_splitter  | BOPAlgo_ShellSplitter   | Shell partitioning              |
//! | checker_si      | BOPAlgo_CheckerSI       | Self-interference check         |

pub mod pave_filler;
pub mod builder;
pub mod builder_face;
pub mod builder_solid;
pub mod shell_splitter;
pub mod checker_si;

// Re-export shared types
pub use builder::{BooleanBuilder, BooleanError};

// ===
// BOPAlgo_GlueEnum ?glue mode for coincident geometry
// ===
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GlueEnum {
    GlueOff = 0,
    GlueShift = 1,
    GlueFull = 2,
}
impl Default for GlueEnum { fn default() -> Self { GlueEnum::GlueOff } }

// ===
// BOPAlgo_Operation ?boolean operation type
// ===
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOpType {
    Union,
    Intersection,
    Difference,
}

// ===
// BOPAlgo_Alert ?diagnostic alerts
// ===
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
    /// BOPAlgo_AlertSelfInterferingShape ?vertex/edge from same argument
    /// with unexpected interference.
    SelfInterferingShape(Vec<usize>),
    /// BOPAlgo_AlertAcquiredSelfIntersection ?acquired self-intersection
    /// detected during CheckSelfInterference.
    AcquiredSelfIntersection(Vec<usize>),
    /// BOPAlgo_AlertTooSmallEdge ?edge range too small for splitting.
    TooSmallEdge(usize),
    /// BOPAlgo_AlertNotSplittableEdge (stub).
    NotSplittableEdge(usize),
    /// BOPAlgo_AlertBadPositioning (stub).
    BadPositioning(Vec<usize>),
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
}
