// OCCT BOPAlgo_ArgumentAnalyzer — input validation for Boolean Operations.
//
// OCCT BOPAlgo_ArgumentAnalyzer.cxx / .hxx (~1000 lines).
// Checks shapes for validity before boolean operations.

/// OCCT BOPAlgo_ArgumentAnalyzer — checks shape validity for Boolean Ops.
///
/// Flags control which checks are performed:
/// - ArgumentTypeMode: shape type compatibility with operation
/// - SelfInterMode: self-intersection of each argument
/// - SmallEdgeMode: micro/small edges
/// - RebuildFaceMode: faces that cannot be rebuilt from their edges
/// - TangentMode: tangent sub-shapes
/// - MergeVertexMode: vertices that should be merged
/// - MergeEdgeMode: edges that should be merged
/// - ContinuityMode: C0 continuity issues
/// - CurveOnSurfaceMode: pcurve deviation from 3D curve
#[derive(Debug, Clone)]
pub struct ArgumentAnalyzer {
    // Input shapes
    my_shape1: Option<TopoShape>,
    my_shape2: Option<TopoShape>,
    // Options
    my_stop_on_first: bool,
    my_operation: i32, // BOPAlgo_Operation enum
    // Mode flags
    my_argument_type_mode: bool,
    my_self_inter_mode: bool,
    my_small_edge_mode: bool,
    my_rebuild_face_mode: bool,
    my_tangent_mode: bool,
    my_merge_vertex_mode: bool,
    my_merge_edge_mode: bool,
    my_continuity_mode: bool,
    my_curve_on_surface_mode: bool,
    // Internal state
    my_empty1: bool,
    my_empty2: bool,
    // Results
    my_result: Vec<CheckResult>,
}

/// OCCT BOPAlgo_CheckResult — result of a single check.
#[derive(Debug, Clone)]
pub struct CheckResult {
    pub check_status: CheckStatus,
}

/// OCCT BOPAlgo_CheckStatus — type of check failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckStatus {
    BadType,
    SelfIntersect,
    TooSmallEdge,
    NonRecoverableFace,
    IncompatibilityOfVertex,
    IncompatibilityOfEdge,
    IncompatibilityOfFace,
    GeomAbsC0,
    InvalidCurveOnSurface,
    OperationAborted,
    CheckUnknown,
}

/// rcad: placeholder shape type for the analyzer.
#[derive(Debug, Clone)]
pub enum TopoShape {
    Null,
    Shape(usize),
}

impl ArgumentAnalyzer {
    /// OCCT BOPAlgo_ArgumentAnalyzer() — empty constructor.
    pub fn new() -> Self {
        ArgumentAnalyzer {
            my_shape1: None,
            my_shape2: None,
            my_stop_on_first: false,
            my_operation: 0, // BOPAlgo_UNKNOWN
            my_argument_type_mode: false,
            my_self_inter_mode: false,
            my_small_edge_mode: false,
            my_rebuild_face_mode: false,
            my_tangent_mode: false,
            my_merge_vertex_mode: false,
            my_merge_edge_mode: false,
            my_continuity_mode: false,
            my_curve_on_surface_mode: false,
            my_empty1: false,
            my_empty2: false,
            my_result: Vec::new(),
        }
    }

    // ── Setters / Getters (OCCT 1:1) ──────────────────────────────────

    pub fn set_shape1(&mut self, shape: TopoShape) { self.my_shape1 = Some(shape); }
    pub fn set_shape2(&mut self, shape: TopoShape) { self.my_shape2 = Some(shape); }
    pub fn get_shape1(&self) -> &Option<TopoShape> { &self.my_shape1 }
    pub fn get_shape2(&self) -> &Option<TopoShape> { &self.my_shape2 }

    pub fn operation_type(&mut self) -> &mut i32 { &mut self.my_operation }
    pub fn stop_on_first_faulty(&mut self) -> &mut bool { &mut self.my_stop_on_first }

    pub fn argument_type_mode(&mut self) -> &mut bool { &mut self.my_argument_type_mode }
    pub fn self_inter_mode(&mut self) -> &mut bool { &mut self.my_self_inter_mode }
    pub fn small_edge_mode(&mut self) -> &mut bool { &mut self.my_small_edge_mode }
    pub fn rebuild_face_mode(&mut self) -> &mut bool { &mut self.my_rebuild_face_mode }
    pub fn tangent_mode(&mut self) -> &mut bool { &mut self.my_tangent_mode }
    pub fn merge_vertex_mode(&mut self) -> &mut bool { &mut self.my_merge_vertex_mode }
    pub fn merge_edge_mode(&mut self) -> &mut bool { &mut self.my_merge_edge_mode }
    pub fn continuity_mode(&mut self) -> &mut bool { &mut self.my_continuity_mode }
    pub fn curve_on_surface_mode(&mut self) -> &mut bool { &mut self.my_curve_on_surface_mode }

    // ── Perform ───────────────────────────────────────────────────────

    /// OCCT BOPAlgo_ArgumentAnalyzer::Perform().
    /// Runs all enabled checks.
    pub fn perform(&mut self) {
        self.my_result.clear();

        // 1. Prepare
        self.prepare();

        // 2-10. Run tests based on mode flags
        if self.my_argument_type_mode {
            self.test_types();
        }
        if self.my_self_inter_mode {
            self.test_self_interferences();
            if self.has_stopper() { return; }
        }
        if self.my_small_edge_mode && self.should_run() {
            self.test_small_edge();
            if self.has_stopper() { return; }
        }
        if self.my_rebuild_face_mode && self.should_run() {
            self.test_rebuild_face();
            if self.has_stopper() { return; }
        }
        if self.my_tangent_mode && self.should_run() {
            self.test_tangent();
            if self.has_stopper() { return; }
        }
        if self.my_merge_vertex_mode && self.should_run() {
            self.test_merge_vertex();
            if self.has_stopper() { return; }
        }
        if self.my_merge_edge_mode && self.should_run() {
            self.test_merge_edge();
            if self.has_stopper() { return; }
        }
        if self.my_continuity_mode && self.should_run() {
            self.test_continuity();
            if self.has_stopper() { return; }
        }
        if self.my_curve_on_surface_mode && self.should_run() {
            self.test_curve_on_surface();
        }
    }

    /// Returns true if there is a faulty result.
    pub fn has_faulty(&self) -> bool {
        !self.my_result.is_empty()
    }

    /// Returns the check results.
    pub fn get_check_result(&self) -> &[CheckResult] {
        &self.my_result
    }

    // ── Protected methods (OCCT naming) ──────────────────────────────

    fn prepare(&mut self) {
        // rcad: stub — OCCT checks IsEmptyShape via BOPTools_AlgoTools3D
    }

    fn test_types(&mut self) {
        // rcad: stub — OCCT checks dimension compatibility with operation
    }

    fn test_self_interferences(&mut self) {
        // rcad: stub — requires BOPAlgo_CheckerSI
    }

    fn test_small_edge(&mut self) {
        // rcad: stub — requires BRepExtrema_DistShapeShape
    }

    fn test_rebuild_face(&mut self) {
        // rcad: stub — requires BOPAlgo_BuilderFace
    }

    fn test_tangent(&mut self) {
        // OCCT: not implemented
    }

    fn test_merge_sub_shapes(&mut self, _the_type: u32) {
        // rcad: stub — iterates sub-shapes between shape1/shape2
    }

    fn test_merge_vertex(&mut self) {
        self.test_merge_sub_shapes(0); // TopAbs_VERTEX
    }

    fn test_merge_edge(&mut self) {
        self.test_merge_sub_shapes(1); // TopAbs_EDGE
    }

    fn test_continuity(&mut self) {
        // rcad: stub — checks GeomAbs_C0 edges/faces
    }

    fn test_curve_on_surface(&mut self) {
        // rcad: stub — requires ComputeTolerance on face/edge pairs
    }

    // ── Helpers ──────────────────────────────────────────────────────

    fn has_stopper(&self) -> bool {
        self.my_stop_on_first && self.has_faulty()
    }

    fn should_run(&self) -> bool {
        self.my_result.is_empty() || !self.my_stop_on_first
    }
}

impl Default for ArgumentAnalyzer {
    fn default() -> Self {
        Self::new()
    }
}
