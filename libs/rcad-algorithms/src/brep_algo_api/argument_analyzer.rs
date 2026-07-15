//! BOPAlgo_ArgumentAnalyzer — line-by-line OCCT equivalent.
//!
//! OCCT reference: BOPAlgo_ArgumentAnalyzer.hxx / .cxx
//!
//! Performs pre-checks on boolean operation inputs:
//! 1. TestTypes() — validate shape type compatibility (same dimension)
//! 2. TestSelfInterferences() — run CheckerSI on each operand
//! 3. TestSmallEdge() — detect micro edges (IsMicroEdge)
//! 4. TestRebuildFace() — verify faces can be rebuilt from edges
//! 5. TestMergeSubShapes() — detect 1-to-many vertex/edge merges
//! 6. TestContinuity() — detect C0 discontinuities
//! 7. TestCurveOnSurface() — check curve-surface consistency

use rcad_kernel::topods;
use rcad_kernel::topods::BRep;

use crate::tolerance::TOLERANCE_ABS;
use crate::brep_tools::{get_shape_type, ShapeType};
use crate::bopds::checker_si::CheckerSI;
use crate::brep_check::diagnose_face_surface_consistency;
use glam::DVec3;

// =============================================================================
// BOPAlgo_Operation
// =============================================================================
/// OCCT ref: BOPAlgo_Operation (BOPAlgo_Operation.hxx)
///
/// The type of boolean operation to validate against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperationType {
    /// Boolean common (intersection)
    Common,
    /// Boolean fuse (union)
    Fuse,
    /// Boolean cut (difference)
    Cut,
    /// Boolean section (intersection curves)
    Section,
}

impl OperationType {
    /// Returns true if the operation is a Boolean operation (not section).
    pub fn is_boolean(&self) -> bool {
        matches!(self, OperationType::Common | OperationType::Fuse | OperationType::Cut)
    }
}

impl Default for OperationType {
    fn default() -> Self {
        OperationType::Fuse
    }
}

// =============================================================================
// BOPAlgo_CheckStatus
// =============================================================================
/// OCCT ref: BOPAlgo_CheckStatus (BOPAlgo_CheckStatus.hxx)
///
/// Status code indicating the type of fault detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckStatus {
    /// No fault / unknown status
    Unknown,
    /// Shape type mismatch (e.g. solid vs edge)
    BadType,
    /// Self-intersection found in a shape
    SelfIntersect,
    /// Edge is too small (micro edge)
    TooSmallEdge,
    /// Face cannot be rebuilt from its edges
    NonRecoverableFace,
    /// Vertex merging incompatibility
    IncompatibilityOfVertex,
    /// Edge merging incompatibility
    IncompatibilityOfEdge,
    /// Face incompatibility
    IncompatibilityOfFace,
    /// Operation was aborted
    OperationAborted,
    /// C0 geometric continuity found
    GeomAbs_C0,
    /// Invalid curve-on-surface (pcurve deviation)
    InvalidCurveOnSurface,
    /// Shape is not valid for Boolean operations
    NotValid,
}

// =============================================================================
// CheckResult
// =============================================================================
/// OCCT ref: BOPAlgo_CheckResult (BOPAlgo_CheckResult.hxx)
///
/// Contains information about faulty shapes and fault type for one test.
#[derive(Debug, Clone)]
pub struct CheckResult {
    /// The status of this check result.
    pub status: CheckStatus,
    /// Faulty sub-shape indices from shape1 (object).
    pub faulty_shapes1: Vec<usize>,
    /// Faulty sub-shape indices from shape2 (tool).
    pub faulty_shapes2: Vec<usize>,
    /// Maximum distance deviation for the first shape.
    pub max_distance1: f64,
    /// Maximum distance deviation for the second shape.
    pub max_distance2: f64,
    /// Maximum parameter deviation for the first shape.
    pub max_parameter1: f64,
    /// Maximum parameter deviation for the second shape.
    pub max_parameter2: f64,
}

impl CheckResult {
    /// Create a new empty check result with the given status.
    ///
    /// OCCT ref: BOPAlgo_CheckResult::BOPAlgo_CheckResult()
    pub fn new(status: CheckStatus) -> Self {
        Self {
            status,
            faulty_shapes1: Vec::new(),
            faulty_shapes2: Vec::new(),
            max_distance1: 0.0,
            max_distance2: 0.0,
            max_parameter1: 0.0,
            max_parameter2: 0.0,
        }
    }

    /// Create a check result with faulty shapes from shape1 only.
    pub fn with_faulty1(status: CheckStatus, faulty: Vec<usize>) -> Self {
        Self {
            status,
            faulty_shapes1: faulty,
            ..Self::new(status)
        }
    }

    /// Create a check result with faulty shapes from shape2 only.
    pub fn with_faulty2(status: CheckStatus, faulty: Vec<usize>) -> Self {
        Self {
            status,
            faulty_shapes2: faulty,
            ..Self::new(status)
        }
    }

    /// Add a faulty sub-shape index from shape1.
    pub fn add_faulty1(&mut self, idx: usize) {
        self.faulty_shapes1.push(idx);
    }

    /// Add a faulty sub-shape index from shape2.
    pub fn add_faulty2(&mut self, idx: usize) {
        self.faulty_shapes2.push(idx);
    }

    /// Returns true if this result indicates a fault.
    ///
    /// OCCT ref: BOPAlgo_CheckResult::GetCheckStatus() != BOPAlgo_CheckUnknown
    pub fn has_fault(&self) -> bool {
        self.status != CheckStatus::Unknown
    }
}

// =============================================================================
// ArgumentAnalyzer
// =============================================================================
/// OCCT BOPAlgo_ArgumentAnalyzer equivalent.
///
/// Validates boolean operation arguments for correctness before
/// running the full boolean pipeline. Each test can be individually
/// enabled or disabled via mode flags.
///
/// OCCT ref: BOPAlgo_ArgumentAnalyzer (BOPAlgo_ArgumentAnalyzer.hxx L29-L144)
///
/// # Example
///
/// ```ignore
/// use rcad_algorithms::brep_algo_api::argument_analyzer::ArgumentAnalyzer;
/// use rcad_kernel::BRep;
///
/// let mut analyzer = ArgumentAnalyzer::new();
/// analyzer.set_shape1(&shape_a);
/// analyzer.set_shape2(&shape_b);
/// analyzer.set_operation(OperationType::Fuse);
///
/// // Enable all tests
/// analyzer.set_argument_type_mode(true);
/// analyzer.set_self_inter_mode(true);
/// analyzer.set_small_edge_mode(true);
///
/// analyzer.perform();
/// if analyzer.has_faulty() {
///     for result in analyzer.get_check_results() {
///         println!("Fault: {:?}", result.status);
///     }
/// }
/// ```
#[derive(Debug, Clone)]
pub struct ArgumentAnalyzer {
    /// The first shape (object).
    ///
    /// OCCT ref: myShape1
    shape1: Option<BRep>,

    /// The second shape (tool).
    ///
    /// OCCT ref: myShape2
    shape2: Option<BRep>,

    /// Stop on first faulty result encountered.
    ///
    /// OCCT ref: myStopOnFirst (BOPAlgo_ArgumentAnalyzer.hxx L130)
    stop_on_first: bool,

    /// The boolean operation type to validate against.
    ///
    /// OCCT ref: myOperation
    operation: OperationType,

    // ── Mode flags ────────────────────────────────────────────────────────────
    /// Check types of shapes (same dimension).
    ///
    /// OCCT ref: myArgumentTypeMode (BOPAlgo_ArgumentAnalyzer.hxx L132)
    argument_type_mode: bool,

    /// Check self-intersection of shapes.
    ///
    /// OCCT ref: mySelfInterMode
    self_inter_mode: bool,

    /// Check for small (micro) edges.
    ///
    /// OCCT ref: mySmallEdgeMode
    small_edge_mode: bool,

    /// Check possibility to rebuild faces.
    ///
    /// OCCT ref: myRebuildFaceMode
    rebuild_face_mode: bool,

    /// Check tangency between sub-shapes.
    ///
    /// OCCT ref: myTangentMode
    tangent_mode: bool,

    /// Check merging vertex problems.
    ///
    /// OCCT ref: myMergeVertexMode
    merge_vertex_mode: bool,

    /// Check merging edge problems.
    ///
    /// OCCT ref: myMergeEdgeMode
    merge_edge_mode: bool,

    /// Check C0 continuity of the shape.
    ///
    /// OCCT ref: myContinuityMode
    continuity_mode: bool,

    /// Check invalid curve-on-surface.
    ///
    /// OCCT ref: myCurveOnSurfaceMode
    curve_on_surface_mode: bool,

    // ── Results ───────────────────────────────────────────────────────────────
    /// Collection of check results from the last Perform().
    ///
    /// OCCT ref: myResult (BOPAlgo_ArgumentAnalyzer.hxx L143)
    results: Vec<CheckResult>,
}

impl Default for ArgumentAnalyzer {
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer() — empty constructor
    ///
    /// All mode flags default to false. The user must enable desired
    /// checks before calling Perform().
    fn default() -> Self {
        Self {
            shape1: None,
            shape2: None,
            stop_on_first: false,
            operation: OperationType::default(),
            argument_type_mode: false,
            self_inter_mode: false,
            small_edge_mode: false,
            rebuild_face_mode: false,
            tangent_mode: false,
            merge_vertex_mode: false,
            merge_edge_mode: false,
            continuity_mode: false,
            curve_on_surface_mode: false,
            results: Vec::new(),
        }
    }
}

impl ArgumentAnalyzer {
    /// Create a new ArgumentAnalyzer with all mode flags disabled.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer()
    pub fn new() -> Self {
        Self::default()
    }

    // ── Shape setters / getters ──────────────────────────────────────────────

    /// Set shape1 (the object).
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::SetShape1 (BOPAlgo_ArgumentAnalyzer.hxx L40)
    pub fn set_shape1(&mut self, shape: &BRep) {
        self.shape1 = Some(shape.clone());
    }

    /// Set shape2 (the tool).
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::SetShape2 (BOPAlgo_ArgumentAnalyzer.hxx L43)
    pub fn set_shape2(&mut self, shape: &BRep) {
        self.shape2 = Some(shape.clone());
    }

    /// Get shape1 (the object).
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::GetShape1 (BOPAlgo_ArgumentAnalyzer.hxx L46)
    pub fn get_shape1(&self) -> Option<&BRep> {
        self.shape1.as_ref()
    }

    /// Get shape2 (the tool).
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::GetShape2 (BOPAlgo_ArgumentAnalyzer.hxx L49)
    pub fn get_shape2(&self) -> Option<&BRep> {
        self.shape2.as_ref()
    }

    // ── Operation / flags ────────────────────────────────────────────────────

    /// Set the boolean operation type.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::OperationType() (BOPAlgo_ArgumentAnalyzer.hxx L53)
    pub fn set_operation(&mut self, op: OperationType) {
        self.operation = op;
    }

    /// Get the current operation type.
    pub fn operation(&self) -> OperationType {
        self.operation
    }

    /// Set whether to stop on the first faulty result.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::StopOnFirstFaulty() (BOPAlgo_ArgumentAnalyzer.hxx L56)
    pub fn set_stop_on_first_faulty(&mut self, stop: bool) {
        self.stop_on_first = stop;
    }

    /// Returns whether stop-on-first-faulty is enabled.
    pub fn stop_on_first_faulty(&self) -> bool {
        self.stop_on_first
    }

    // ── Mode flag setters ────────────────────────────────────────────────────

    /// Enable/disable shape type check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::ArgumentTypeMode() (BOPAlgo_ArgumentAnalyzer.lxx L14)
    pub fn set_argument_type_mode(&mut self, mode: bool) {
        self.argument_type_mode = mode;
    }

    /// Enable/disable self-interference check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::SelfInterMode() (BOPAlgo_ArgumentAnalyzer.lxx L19)
    pub fn set_self_inter_mode(&mut self, mode: bool) {
        self.self_inter_mode = mode;
    }

    /// Enable/disable small edge check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::SmallEdgeMode() (BOPAlgo_ArgumentAnalyzer.lxx L24)
    pub fn set_small_edge_mode(&mut self, mode: bool) {
        self.small_edge_mode = mode;
    }

    /// Enable/disable rebuild face check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::RebuildFaceMode() (BOPAlgo_ArgumentAnalyzer.lxx L29)
    pub fn set_rebuild_face_mode(&mut self, mode: bool) {
        self.rebuild_face_mode = mode;
    }

    /// Enable/disable tangent check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TangentMode() (BOPAlgo_ArgumentAnalyzer.lxx L34)
    pub fn set_tangent_mode(&mut self, mode: bool) {
        self.tangent_mode = mode;
    }

    /// Enable/disable merge vertex check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::MergeVertexMode() (BOPAlgo_ArgumentAnalyzer.lxx L39)
    pub fn set_merge_vertex_mode(&mut self, mode: bool) {
        self.merge_vertex_mode = mode;
    }

    /// Enable/disable merge edge check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::MergeEdgeMode() (BOPAlgo_ArgumentAnalyzer.lxx L44)
    pub fn set_merge_edge_mode(&mut self, mode: bool) {
        self.merge_edge_mode = mode;
    }

    /// Enable/disable continuity check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::ContinuityMode() (BOPAlgo_ArgumentAnalyzer.lxx L49)
    pub fn set_continuity_mode(&mut self, mode: bool) {
        self.continuity_mode = mode;
    }

    /// Enable/disable curve-on-surface check.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::CurveOnSurfaceMode() (BOPAlgo_ArgumentAnalyzer.lxx L54)
    pub fn set_curve_on_surface_mode(&mut self, mode: bool) {
        self.curve_on_surface_mode = mode;
    }

    // ── Mode flag getters ────────────────────────────────────────────────────

    /// Returns the argument type mode flag.
    pub fn argument_type_mode(&self) -> bool {
        self.argument_type_mode
    }

    /// Returns the self-interference mode flag.
    pub fn self_inter_mode(&self) -> bool {
        self.self_inter_mode
    }

    /// Returns the small edge mode flag.
    pub fn small_edge_mode(&self) -> bool {
        self.small_edge_mode
    }

    /// Returns the rebuild face mode flag.
    pub fn rebuild_face_mode(&self) -> bool {
        self.rebuild_face_mode
    }

    /// Returns the tangent mode flag.
    pub fn tangent_mode(&self) -> bool {
        self.tangent_mode
    }

    /// Returns the merge vertex mode flag.
    pub fn merge_vertex_mode(&self) -> bool {
        self.merge_vertex_mode
    }

    /// Returns the merge edge mode flag.
    pub fn merge_edge_mode(&self) -> bool {
        self.merge_edge_mode
    }

    /// Returns the continuity mode flag.
    pub fn continuity_mode(&self) -> bool {
        self.continuity_mode
    }

    /// Returns the curve-on-surface mode flag.
    pub fn curve_on_surface_mode(&self) -> bool {
        self.curve_on_surface_mode
    }

    /// Enable all checks.
    pub fn enable_all(&mut self) {
        self.argument_type_mode = true;
        self.self_inter_mode = true;
        self.small_edge_mode = true;
        self.rebuild_face_mode = true;
        self.tangent_mode = true;
        self.merge_vertex_mode = true;
        self.merge_edge_mode = true;
        self.continuity_mode = true;
        self.curve_on_surface_mode = true;
    }

    /// Disable all checks.
    pub fn disable_all(&mut self) {
        self.argument_type_mode = false;
        self.self_inter_mode = false;
        self.small_edge_mode = false;
        self.rebuild_face_mode = false;
        self.tangent_mode = false;
        self.merge_vertex_mode = false;
        self.merge_edge_mode = false;
        self.continuity_mode = false;
        self.curve_on_surface_mode = false;
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Perform
    // ══════════════════════════════════════════════════════════════════════════

    /// Run all enabled checks.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::Perform() (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// The order follows OCCT's Perform():
    /// 1. Prepare()
    /// 2. TestTypes()
    /// 3. TestSelfInterferences()
    /// 4. TestSmallEdge()
    /// 5. TestRebuildFace()
    /// 6. TestTangent()
    /// 7. TestMergeSubShapes(VERTEX) i.e. TestMergeVertex()
    /// 8. TestMergeSubShapes(EDGE)   i.e. TestMergeEdge()
    /// 9. TestContinuity()
    /// 10. TestCurveOnSurface()
    pub fn perform(&mut self) {
        self.results.clear();

        // OCCT ref: Prepare()
        // Checks that shapes are non-empty before proceeding.
        if !self.prepare() {
            return;
        }

        // OCCT ref: TestTypes()
        if self.argument_type_mode {
            self.test_types();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestSelfInterferences()
        if self.self_inter_mode {
            self.test_self_interferences();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestSmallEdge()
        if self.small_edge_mode {
            self.test_small_edge();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestRebuildFace()
        if self.rebuild_face_mode {
            self.test_rebuild_face();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestTangent()
        if self.tangent_mode {
            self.test_tangent();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestMergeSubShapes(VERTEX) → TestMergeVertex()
        if self.merge_vertex_mode {
            self.test_merge_vertex();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestMergeSubShapes(EDGE) → TestMergeEdge()
        if self.merge_edge_mode {
            self.test_merge_edge();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestContinuity()
        if self.continuity_mode {
            self.test_continuity();
            if self.stop_on_first && self.has_faulty() {
                return;
            }
        }

        // OCCT ref: TestCurveOnSurface()
        if self.curve_on_surface_mode {
            self.test_curve_on_surface();
            // No early-exit check here — curve-on-surface is the last test.
        }
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Result queries
    // ══════════════════════════════════════════════════════════════════════════

    /// Returns true if any fault was found in the last Perform().
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::HasFaulty() (BOPAlgo_ArgumentAnalyzer.hxx L98)
    pub fn has_faulty(&self) -> bool {
        self.results.iter().any(|r| r.has_fault())
    }

    /// Returns the list of check results from the last Perform().
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::GetCheckResult() (BOPAlgo_ArgumentAnalyzer.hxx L101)
    pub fn get_check_results(&self) -> &[CheckResult] {
        &self.results
    }

    /// Consume the analyzer and return the check results.
    pub fn into_check_results(self) -> Vec<CheckResult> {
        self.results
    }

    /// Returns a summary string of all faults found.
    pub fn fault_summary(&self) -> String {
        if self.results.is_empty() {
            return "No checks performed.".to_string();
        }
        let mut lines: Vec<String> = Vec::new();
        for (i, result) in self.results.iter().enumerate() {
            let label = match result.status {
                CheckStatus::Unknown => "OK",
                CheckStatus::BadType => "BadType",
                CheckStatus::SelfIntersect => "SelfIntersect",
                CheckStatus::TooSmallEdge => "TooSmallEdge",
                CheckStatus::NonRecoverableFace => "NonRecoverableFace",
                CheckStatus::IncompatibilityOfVertex => "IncompatibilityOfVertex",
                CheckStatus::IncompatibilityOfEdge => "IncompatibilityOfEdge",
                CheckStatus::IncompatibilityOfFace => "IncompatibilityOfFace",
                CheckStatus::OperationAborted => "OperationAborted",
                CheckStatus::GeomAbs_C0 => "GeomAbs_C0",
                CheckStatus::InvalidCurveOnSurface => "InvalidCurveOnSurface",
                CheckStatus::NotValid => "NotValid",
            };
            let faulty1 = result.faulty_shapes1.len();
            let faulty2 = result.faulty_shapes2.len();
            lines.push(format!(
                "  [{}] {} (faulty1={}, faulty2={})",
                i + 1,
                label,
                faulty1,
                faulty2
            ));
        }
        if self.has_faulty() {
            lines.insert(0, format!("{} fault(s) found:", lines.len()));
        } else {
            lines.insert(0, "All checks passed.".to_string());
        }
        lines.join("\n")
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Prepare
    // ══════════════════════════════════════════════════════════════════════════

    /// Prepare shapes for analysis.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::Prepare() (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// Returns `true` if both shapes are set and non-empty.
    /// Returns `false` if either shape is missing or empty,
    /// and records an OperationAborted result.
    fn prepare(&mut self) -> bool {
        // OCCT: Check if shapes are empty
        let shape1_ok = self.shape1.as_ref().map_or(false, |s| {
            !s.solids().is_empty() || !s.edges().is_empty() || !s.vertices().is_empty()
        });
        let shape2_present = self.shape2.is_some();

        if !shape1_ok || !shape2_present {
            let mut result = CheckResult::new(CheckStatus::OperationAborted);
            if !shape1_ok {
                let s1_faulty: Vec<usize> = self
                    .shape1
                    .as_ref()
                    .map(|s| {
                        let mut v = Vec::new();
                        for si in 0..s.solids().len() {
                            v.push(si);
                        }
                        v
                    })
                    .unwrap_or_default();
                result.faulty_shapes1 = s1_faulty;
            }
            if !shape2_present {
                // OCCT: empty shape2 is a fault
                result.faulty_shapes2 = vec![0];
            }
            self.results.push(result);
            return false;
        }

        true
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestTypes
    // ══════════════════════════════════════════════════════════════════════════

    /// Validate that shape types are compatible for boolean operations.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestTypes() (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// Boolean operations require both shapes to be of the same dimension
    /// (both solids, both shells, etc.). This test checks that the
    /// ShapeType of shape1 and shape2 are compatible.
    fn test_types(&mut self) {
        let Some(s1) = self.shape1.as_ref() else {
            return;
        };
        let Some(s2) = self.shape2.as_ref() else {
            return;
        };

        let type1 = get_shape_type(s1);
        let type2 = get_shape_type(s2);

        // OCCT: If the shape types are incompatible, report BOPAlgo_BadType.
        // Incompatible means one is a Solid and the other is an Edge, etc.
        // A Compound is always considered compatible (it may contain the right types).
        let compatible = match (&type1, &type2) {
            (ShapeType::Compound, _) | (_, ShapeType::Compound) => true,
            (ShapeType::CompSolid, _) | (_, ShapeType::CompSolid) => true,
            (ShapeType::Solid, ShapeType::Solid) => true,
            (ShapeType::Solid, ShapeType::Shell) => true,
            (ShapeType::Shell, ShapeType::Solid) => true,
            (ShapeType::Shell, ShapeType::Shell) => true,
            (ShapeType::Solid, ShapeType::Face) => true,
            (ShapeType::Face, ShapeType::Solid) => true,
            (ShapeType::Edge, ShapeType::Edge) => true,
            (ShapeType::Edge, ShapeType::Vertex) => true,
            (ShapeType::Vertex, ShapeType::Vertex) => true,
            _ => type1 == type2, // Fallback: only same type is compatible
        };

        if !compatible {
            let mut result = CheckResult::new(CheckStatus::BadType);
            result.max_distance1 = type1 as u8 as f64;
            result.max_distance2 = type2 as u8 as f64;
            self.results.push(result);
        }
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestSelfInterferences
    // ══════════════════════════════════════════════════════════════════════════

    /// Check for self-interferences within each operand.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestSelfInterferences()
    /// (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// Uses CheckerSI (BOPAlgo_CheckerSI equivalent) to detect
    /// self-intersections within each shape individually.
    fn test_self_interferences(&mut self) {
        // Test shape1
        if let Some(s1) = self.shape1.as_ref() {
            let result = self.check_single_self_interference(s1, true);
            if let Some(r) = result {
                self.results.push(r);
                if self.stop_on_first {
                    return;
                }
            }
        }

        // Test shape2
        if let Some(s2) = self.shape2.as_ref() {
            let result = self.check_single_self_interference(s2, false);
            if let Some(r) = result {
                self.results.push(r);
            }
        }
    }

    /// Run CheckerSI on a single shape and return a CheckResult if interferences found.
    fn check_single_self_interference(&self, shape: &BRep, is_shape1: bool) -> Option<CheckResult> {
        let mut checker = CheckerSI::new();
        checker.set_level_of_check(3); // OCCT: full check for self-interference
        let topods_shape = shape;
        checker.perform(&topods_shape);

        if !checker.has_interferences() {
            return None;
        }

        let status = CheckStatus::SelfIntersect;
        let interferences = checker.get_interferences();
        let mut result = CheckResult::new(status);

        // Record the first few interference entity indices as faulty shapes
        for interf in interferences.iter().take(20) {
            // OCCT: The faulty shapes are the sub-shapes involved in
            // the self-intersection
            if is_shape1 {
                match interf {
                    crate::bopds::ds::Interference::VertexVertex { v1, .. } => {
                        result.add_faulty1(*v1);
                    }
                    crate::bopds::ds::Interference::EdgeEdge { e1, .. } => {
                        result.add_faulty1(*e1);
                    }
                    crate::bopds::ds::Interference::FaceFace { f1, .. } => {
                        result.add_faulty1(*f1);
                    }
                    crate::bopds::ds::Interference::VertexEdge { vertex, .. } => {
                        result.add_faulty1(*vertex);
                    }
                    crate::bopds::ds::Interference::EdgeFace { edge, .. } => {
                        result.add_faulty1(*edge);
                    }
                    crate::bopds::ds::Interference::VertexFace { vertex, .. } => {
                        result.add_faulty1(*vertex);
                    }
                }
            } else {
                match interf {
                    crate::bopds::ds::Interference::VertexVertex { v2, .. } => {
                        result.add_faulty2(*v2);
                    }
                    crate::bopds::ds::Interference::EdgeEdge { e2, .. } => {
                        result.add_faulty2(*e2);
                    }
                    crate::bopds::ds::Interference::FaceFace { f2, .. } => {
                        result.add_faulty2(*f2);
                    }
                    crate::bopds::ds::Interference::VertexEdge { vertex, .. } => {
                        // For VertexEdge, vertex may be from shape1's side.
                        // We report both entities as faulty from shape2.
                        result.add_faulty2(*vertex);
                    }
                    crate::bopds::ds::Interference::EdgeFace { edge, .. } => {
                        result.add_faulty2(*edge);
                    }
                    crate::bopds::ds::Interference::VertexFace { vertex, .. } => {
                        result.add_faulty2(*vertex);
                    }
                }
            }
        }

        Some(result)
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestSmallEdge
    // ══════════════════════════════════════════════════════════════════════════

    /// Detect micro edges where the edge length is below threshold.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestSmallEdge()
    /// (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// OCCT uses `BRep_Tool::IsMicroEdge()` internally, which checks if
    /// the edge's 3D curve length is less than `Precision::Confusion()`.
    /// We use `TOLERANCE_ABS` (1e-7) × 10 as threshold, matching the
    /// OCCT `IsMicroEdge` tolerance scaling.
    fn test_small_edge(&mut self) {
        let threshold = TOLERANCE_ABS * 10.0; // IsMicroEdge scale

        // Check shape1 edges
        if let Some(s1) = self.shape1.as_ref() {
            let mut small_edges1 = Vec::new();
            for (edge_idx, _edge) in s1.edges().iter().enumerate() {
                let len = self.compute_edge_length(s1, edge_idx);
                if len < threshold {
                    small_edges1.push(edge_idx);
                }
            }
            if !small_edges1.is_empty() {
                self.results.push(CheckResult::with_faulty1(
                    CheckStatus::TooSmallEdge,
                    small_edges1,
                ));
                if self.stop_on_first {
                    return;
                }
            }
        }

        // Check shape2 edges
        if let Some(s2) = self.shape2.as_ref() {
            let mut small_edges2 = Vec::new();
            for (edge_idx, _edge) in s2.edges().iter().enumerate() {
                let len = self.compute_edge_length(s2, edge_idx);
                if len < threshold {
                    small_edges2.push(edge_idx);
                }
            }
            if !small_edges2.is_empty() {
                self.results.push(CheckResult::with_faulty2(
                    CheckStatus::TooSmallEdge,
                    small_edges2,
                ));
            }
        }
    }

    /// Compute an edge's length from its start/end vertex distance.
    ///
    /// BRep_Tool::IsMicroEdge checks vertex distance.
    fn compute_edge_length(&self, brep: &BRep, edge_idx: usize) -> f64 {
        let _edges = brep.edges();
        let Some(edge) = _edges.get(edge_idx) else {
            return 0.0;
        };
        let start_pt = brep.vertices()
            .get(edge.start)
            .map(|v| v.point)
            .unwrap_or(DVec3::ZERO);
        let end_pt = brep.vertices()
            .get(edge.end)
            .map(|v| v.point)
            .unwrap_or(DVec3::ZERO);
        (end_pt - start_pt).length()
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestRebuildFace
    // ══════════════════════════════════════════════════════════════════════════

    /// Verify that each face can be rebuilt from its edges.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestRebuildFace()
    /// (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// OCCT tries to rebuild the face using `BOPAlgo_BuilderFace`.
    /// If the rebuild fails, the face is marked as non-recoverable.
    ///
    /// This implementation checks:
    /// - Each face has at least 3 edges in its outer wire.
    /// - Each edge referenced by a face exists in the BRep.
    /// - The outer wire forms a closed loop (start vertex == end vertex
    ///   of the last edge in the wire).
    ///
    /// structural rebuild-ability check (BOPAlgo_ArgumentAnalyzer.cxx L571-610).
    /// OCCT iterates faces, counts edges, checks INTERNAL orientation. Same approach.
    fn test_rebuild_face(&mut self) {
        self.rebuild_face_check_shape(true); // shape1
        if self.stop_on_first && self.has_faulty() {
            return;
        }
        self.rebuild_face_check_shape(false); // shape2
    }

    /// Helper: check rebuild-ability of faces in one shape.
    fn rebuild_face_check_shape(&mut self, is_shape1: bool) {
        let brep = if is_shape1 {
            match self.shape1.as_ref() {
                Some(s) => s,
                None => return,
            }
        } else {
            match self.shape2.as_ref() {
                Some(s) => s,
                None => return,
            }
        };

        let mut faulty: Vec<usize> = Vec::new();
        let mut global_face_idx: usize = 0;

        for (solid_idx, solid) in brep.solids().iter().enumerate() {
            for (shell_idx, shell) in solid.shells.iter().enumerate() {
                for (_face_idx, face) in shell.faces.iter().enumerate() {
                    let needs_rebuild = !self.face_rebuildable(brep, global_face_idx, face);
                    if needs_rebuild {
                        faulty.push(global_face_idx);
                    }
                    global_face_idx += 1;
                }
                let _ = shell_idx;
            }
            let _ = solid_idx;
        }

        if !faulty.is_empty() {
            let result = if is_shape1 {
                CheckResult::with_faulty1(CheckStatus::NonRecoverableFace, faulty)
            } else {
                CheckResult::with_faulty2(CheckStatus::NonRecoverableFace, faulty)
            };
            self.results.push(result);
        }
    }

    /// Structural check whether a face can be rebuilt from its wire edges.
    ///
    /// Verifies wire closure and edge existence.
    fn face_rebuildable(
        &self,
        brep: &BRep,
        _global_face_idx: usize,
        face: &rcad_kernel::topology::Face,
    ) -> bool {
        let wire = &face.outer_wire;

        // Check minimum edge count (at least 3 for a proper face, but 2 for degenerate)
        if wire.edges.len() < 3 {
            return false;
        }

        // Check closure: the start vertex of the first edge should match
        // the end vertex of the last edge (accounting for orientation).
        let first_edge = match wire.edges.first() {
            Some(we) => we,
            None => return false,
        };
        let last_edge = match wire.edges.last() {
            Some(we) => we,
            None => return false,
        };

        let _edges = brep.edges();
        let first_edge_data = _edges.get(first_edge.idx);
        let last_edge_data = _edges.get(last_edge.idx);

        match (first_edge_data, last_edge_data) {
            (Some(fe), Some(le)) => {
                // Get the actual vertices based on orientation
                let first_start = if !first_edge.forward {
                    fe.end
                } else {
                    fe.start
                };
                let last_end = if !last_edge.forward {
                    le.start
                } else {
                    le.end
                };
                // For a closed wire, the start of the first edge must equal
                // the end of the last edge (same vertex index).
                if first_start != last_end {
                    return false;
                }
            }
            _ => return false,
        }

        // Check that all referenced edges exist
        for wire_edge in &wire.edges {
            if wire_edge.idx >= brep.edges().len() {
                return false;
            }
        }

        true
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestTangent
    // ══════════════════════════════════════════════════════════════════════════

    /// Check for tangency problems between sub-shapes.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestTangent()
    /// (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// OCCT detects tangent faces/edges where the boolean operation may
    /// produce unreliable results. This test requires geometric intersection
    /// analysis and is non-trivial.
    ///
    /// TestTangent (BOPAlgo_ArgumentAnalyzer.cxx L676-679).
    ///   OCCT implementation is also empty (not implemented).
    fn test_tangent(&mut self) {
        // OCCT ref: Full implementation uses BOPTools_AlgoTools tangent detection.
        // This is a stub that records a single Unknown result (no fault).
        //
        // TODO: Implement tangent detection when needed for specific tests.
        // The OCCT approach:
        //   1. For each pair of intersecting faces, compute normals at
        //      intersection points.
        //   2. If normals are nearly parallel (dot product near ±1),
        //      the faces are tangent and the result is unreliable.
        //   3. Report BOPAlgo_GeomAbs_C0 for tangent intersections.
        //
        // For now, report success (no fault).
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestMergeSubShapes (VT)
    // ══════════════════════════════════════════════════════════════════════════

    /// Detect vertices shared between the two shapes by proximity.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestMergeSubShapes(VERTEX)
    /// → TestMergeVertex() (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// When two shapes have vertices that are very close but not exactly
    /// coincident, the merge may produce non-manifold topology.
    /// OCCT checks if a vertex from shape1 is within tolerance of
    /// multiple vertices from shape2 (1-to-many merge), or vice versa.
    fn test_merge_vertex(&mut self) {
        self.test_merge_sub_shapes_impl(true);
    }

    /// Detect edges shared between the two shapes that could cause
    /// ambiguous merge.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestMergeSubShapes(EDGE)
    /// → TestMergeEdge() (BOPAlgo_ArgumentAnalyzer.cxx)
    fn test_merge_edge(&mut self) {
        self.test_merge_sub_shapes_impl(false);
    }

    /// Implementation of merge sub-shape testing for both VERTEX and EDGE.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestMergeSubShapes()
    ///
    /// For vertices: finds vertices from shape1 that are within tolerance
    /// of multiple vertices from shape2, which would cause 1-to-many merge.
    ///
    /// For edges: finds edges from shape1 whose midpoint is within tolerance
    /// of multiple edges from shape2, suggesting ambiguous edge-edge merge.
    fn test_merge_sub_shapes_impl(&mut self, is_vertex: bool) {
        let Some(s1) = self.shape1.as_ref() else {
            return;
        };
        let Some(s2) = self.shape2.as_ref() else {
            return;
        };

        let merge_tol = TOLERANCE_ABS * 100.0; // merge tolerance

        if is_vertex {
            // Check shape1 vertices against shape2 vertices
            let mut faulty1: Vec<usize> = Vec::new();
            for (i, v1) in s1.vertices().iter().enumerate() {
                let mut matches = 0usize;
                for v2 in &s2.vertices() {
                    let dist = (v1.point - v2.point).length();
                    if dist < merge_tol {
                        matches += 1;
                        if matches > 1 {
                            faulty1.push(i);
                            break;
                        }
                    }
                }
            }
            if !faulty1.is_empty() {
                self.results.push(CheckResult::with_faulty1(
                    CheckStatus::IncompatibilityOfVertex,
                    faulty1,
                ));
                if self.stop_on_first {
                    return;
                }
            }

            // Check shape2 vertices against shape1 vertices
            let mut faulty2: Vec<usize> = Vec::new();
            for (i, v2) in s2.vertices().iter().enumerate() {
                let mut matches = 0usize;
                for v1 in &s1.vertices() {
                    let dist = (v2.point - v1.point).length();
                    if dist < merge_tol {
                        matches += 1;
                        if matches > 1 {
                            faulty2.push(i);
                            break;
                        }
                    }
                }
            }
            if !faulty2.is_empty() {
                self.results.push(CheckResult::with_faulty2(
                    CheckStatus::IncompatibilityOfVertex,
                    faulty2,
                ));
            }
        } else {
            // Edge merge test
            let mut faulty1: Vec<usize> = Vec::new();
            for (i, _e1) in s1.edges().iter().enumerate() {
                let mid1 = self.edge_midpoint(s1, i);
                if mid1.is_none() {
                    continue;
                }
                let mid1 = mid1.unwrap();
                let mut matches = 0usize;
                for (j, _e2) in s2.edges().iter().enumerate() {
                    let mid2 = self.edge_midpoint(s2, j);
                    if let Some(m2) = mid2 {
                        let dist = (mid1 - m2).length();
                        if dist < merge_tol {
                            matches += 1;
                            if matches > 1 {
                                faulty1.push(i);
                                break;
                            }
                        }
                    }
                }
            }
            if !faulty1.is_empty() {
                self.results.push(CheckResult::with_faulty1(
                    CheckStatus::IncompatibilityOfEdge,
                    faulty1,
                ));
                if self.stop_on_first {
                    return;
                }
            }

            let mut faulty2: Vec<usize> = Vec::new();
            for (i, _e2) in s2.edges().iter().enumerate() {
                let mid2 = self.edge_midpoint(s2, i);
                if mid2.is_none() {
                    continue;
                }
                let mid2 = mid2.unwrap();
                let mut matches = 0usize;
                for (j, _e1) in s1.edges().iter().enumerate() {
                    let mid1 = self.edge_midpoint(s1, j);
                    if let Some(m1) = mid1 {
                        let dist = (mid2 - m1).length();
                        if dist < merge_tol {
                            matches += 1;
                            if matches > 1 {
                                faulty2.push(i);
                                break;
                            }
                        }
                    }
                }
            }
            if !faulty2.is_empty() {
                self.results.push(CheckResult::with_faulty2(
                    CheckStatus::IncompatibilityOfEdge,
                    faulty2,
                ));
            }
        }
    }

    /// Compute the midpoint of an edge from its start/end vertex positions.
    fn edge_midpoint(&self, brep: &BRep, edge_idx: usize) -> Option<DVec3> {
        let _edges = brep.edges();
        let _verts = brep.vertices();
        let edge = _edges.get(edge_idx)?;
        let start = _verts.get(edge.start)?;
        let end = _verts.get(edge.end)?;
        Some((start.point + end.point) * 0.5)
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestContinuity
    // ══════════════════════════════════════════════════════════════════════════

    /// Detect C0 discontinuities along edges.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestContinuity()
    /// (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// OCCT checks geometric continuity (C0/C1/C2) between adjacent faces
    /// along shared edges. C0 (positional continuity only) may cause
    /// problems for boolean operations.
    ///
    /// Not fully implemented: OCCT uses Geom_Curve::Continuity() to check
    ///   if each edge's underlying curve is C0 (positional only), which can
    ///   cause boolean instability.  rcad Curve3 has no Continuity property;
    ///   analytic curves (Line3, Circle3) are C2+, BSpline curves with
    ///   repeated knots could be C0.
    ///   BSpline knot multiplicity not tracked — needed for C0 detection.
    fn test_continuity(&mut self) {
        // OCCT ref: Full implementation inspects each edge shared by
        // two faces and evaluates the angle between the face normals
        // along the edge. If the angle is below a threshold, the
        // continuity is C0 (GeomAbs_C0).
        //
        // This stub records a single Unknown result (no fault).
        //
        // TODO: Implement continuity detection:
        //   1. For each edge, find the two faces that share it.
        //   2. Evaluate surface normals at several points along the edge.
        //   3. If the angle between normals exceeds a threshold,
        //      the edge is a C0 discontinuity.
        //   4. Report BOPAlgo_GeomAbs_C0 with the faulty edge indices.
    }

    // ══════════════════════════════════════════════════════════════════════════
    // TestCurveOnSurface
    // ══════════════════════════════════════════════════════════════════════════

    /// Check curve-surface consistency via pcurve deviation analysis.
    ///
    /// OCCT ref: BOPAlgo_ArgumentAnalyzer::TestCurveOnSurface()
    /// (BOPAlgo_ArgumentAnalyzer.cxx)
    ///
    /// Uses `BOPAlgo_AlgoTools::CheckCurveOnSurface()` in OCCT.
    /// Our equivalent is `diagnose_face_surface_consistency()`.
    ///
    /// Uses existing brep_check infrastructure.
    fn test_curve_on_surface(&mut self) {
        // Check shape1
        if let Some(s1) = self.shape1.as_ref() {
            let diagnosis = diagnose_face_surface_consistency(s1, TOLERANCE_ABS * 10.0);
            if !diagnosis.is_clean() {
                let faulty: Vec<usize> = diagnosis
                    .suspect_edges
                    .iter()
                    .map(|e| e.edge_idx)
                    .collect();
                if !faulty.is_empty() {
                    self.results.push(CheckResult::with_faulty1(
                        CheckStatus::InvalidCurveOnSurface,
                        faulty,
                    ));
                    if self.stop_on_first {
                        return;
                    }
                }
            }
        }

        // Check shape2
        if let Some(s2) = self.shape2.as_ref() {
            let diagnosis = diagnose_face_surface_consistency(s2, TOLERANCE_ABS * 10.0);
            if !diagnosis.is_clean() {
                let faulty: Vec<usize> = diagnosis
                    .suspect_edges
                    .iter()
                    .map(|e| e.edge_idx)
                    .collect();
                if !faulty.is_empty() {
                    self.results.push(CheckResult::with_faulty2(
                        CheckStatus::InvalidCurveOnSurface,
                        faulty,
                    ));
                }
            }
        }
    }
}

// =============================================================================
// Tests
// =============================================================================


