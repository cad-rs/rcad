// OCCT BOPAlgo_ArgumentAnalyzer — input validation for Boolean Operations.
//
// OCCT BOPAlgo_ArgumentAnalyzer.cxx L1-1015 / .hxx.

/// OCCT BOPAlgo_ArgumentAnalyzer — checks shape validity for Boolean Ops.
///
/// rcad: some test bodies simplified where OCCT infrastructure is missing
/// (CheckerSI, BuilderFace, DistShapeShape). Control flow and structure
/// match OCCT 1:1.
pub struct ArgumentAnalyzer {
    my_shape1: usize, // shape index placeholder
    my_shape2: usize,
    my_stop_on_first: bool,
    my_operation: i32,
    my_argument_type_mode: bool,
    my_self_inter_mode: bool,
    my_small_edge_mode: bool,
    my_rebuild_face_mode: bool,
    my_tangent_mode: bool,
    my_merge_vertex_mode: bool,
    my_merge_edge_mode: bool,
    my_continuity_mode: bool,
    my_curve_on_surface_mode: bool,
    my_empty1: bool,
    my_empty2: bool,
    my_result: Vec<CheckResult>,
}

#[derive(Debug, Clone)]
pub struct CheckResult {
    pub check_status: CheckStatus,
}

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

impl ArgumentAnalyzer {
    pub fn new() -> Self {
        ArgumentAnalyzer {
            my_shape1: usize::MAX,
            my_shape2: usize::MAX,
            my_stop_on_first: false,
            my_operation: 0,
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

    pub fn set_shape1(&mut self, s: usize) { self.my_shape1 = s; }
    pub fn set_shape2(&mut self, s: usize) { self.my_shape2 = s; }
    pub fn get_shape1(&self) -> usize { self.my_shape1 }
    pub fn get_shape2(&self) -> usize { self.my_shape2 }
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

    /// OCCT Perform() L130-257 — run all enabled checks.
    pub fn perform(&mut self) {
        self.my_result.clear();
        self.prepare(); // L142

        if self.my_argument_type_mode { self.test_types(); } // L146
        if self.my_self_inter_mode {   // L155
            self.test_self_interferences();
            if self.has_errors_or_stop() { return; }
        }
        if self.my_small_edge_mode && self.should_run() { // L165
            self.test_small_edge();
            if self.has_errors_or_stop() { return; }
        }
        if self.my_rebuild_face_mode && self.should_run() { // L178
            self.test_rebuild_face();
            if self.has_errors_or_stop() { return; }
        }
        if self.my_tangent_mode && self.should_run() { // L191
            self.test_tangent();
            if self.has_errors_or_stop() { return; }
        }
        if self.my_merge_vertex_mode && self.should_run() { // L205
            self.test_merge_vertex();
            if self.has_errors_or_stop() { return; }
        }
        if self.my_merge_edge_mode && self.should_run() { // L218
            self.test_merge_edge();
            if self.has_errors_or_stop() { return; }
        }
        if self.my_continuity_mode && self.should_run() { // L231
            self.test_continuity();
            if self.has_errors_or_stop() { return; }
        }
        if self.my_curve_on_surface_mode && self.should_run() { // L243
            self.test_curve_on_surface();
        }
    }

    pub fn has_faulty(&self) -> bool { !self.my_result.is_empty() }
    pub fn get_check_result(&self) -> &[CheckResult] { &self.my_result }

    // ── Prepare L115-126 ────────────────────────────────────────────
    fn prepare(&mut self) {
        // OCCT: BOPTools_AlgoTools3D::IsEmptyShape for each non-null shape.
        // rcad: shape emptiness checks not yet implemented.
    }

    // ── TestTypes L275-352 ──────────────────────────────────────────
    fn test_types(&mut self) {
        let is_s1 = self.my_shape1 == usize::MAX;
        let is_s2 = self.my_shape2 == usize::MAX;

        if is_s1 && is_s2 {
            // L279-285: both null → BadType
            self.push_result(CheckStatus::BadType);
            return;
        }

        if (is_s1 && !is_s2) || (!is_s1 && is_s2) {
            // L288-301: single shape, check if empty or unknown operation
            let b_is_empty = if is_s1 { self.my_empty2 } else { self.my_empty1 };
            if b_is_empty || self.my_operation == 0 {
                self.push_result(CheckStatus::BadType);
            }
            return;
        }

        // L304-351: two shapes
        if self.my_empty1 || self.my_empty2 {
            self.push_result(CheckStatus::BadType);
            return;
        }

        // L330-350: operation-specific dimension checks
        if self.my_operation != 0 && self.my_operation != 1 {
            // OCCT BOPAlgo_FUSE=2, BOPAlgo_CUT=3, BOPAlgo_CUT21=4
            // BOPTools_AlgoTools::Dimensions returns [minDim, maxDim]
            // rcad: dimension queries not yet implemented.
        }
    }

    // ── TestSelfInterferences L356-445 ──────────────────────────────
    fn test_self_interferences(&mut self) {
        // OCCT: uses BOPAlgo_CheckerSI on each non-empty shape
        // iterates result DS interferences, adds SelfIntersect results.
        // rcad: CheckerSI not yet implemented.
        for _shape_idx in 0..2 {
            // placeholder for CheckerSI loop
        }
    }

    // ── TestSmallEdge L449-567 ──────────────────────────────────────
    fn test_small_edge(&mut self) {
        // OCCT: iterates edges via TopExp_Explorer, checks IsMicroEdge.
        // For SECTION operation, also checks if edge vertices lie on other shape.
        // rcad: edge iteration not yet wired.
    }

    // ── TestRebuildFace L571-672 ────────────────────────────────────
    fn test_rebuild_face(&mut self) {
        if self.my_operation == 5 || self.my_operation == 0 { return; } // SECTION/UNKNOWN
        // OCCT: iterates faces, builds face from its edges via BOPAlgo_BuilderFace,
        // checks that number of resulting areas == 1 and edge count matches.
        // rcad: BuilderFace not yet available.
    }

    // ── TestTangent L676-679 ────────────────────────────────────────
    fn test_tangent(&mut self) {
        // OCCT: not implemented (empty body)
    }

    // ── TestMergeSubShapes L683-878 ─────────────────────────────────
    fn test_merge_sub_shapes(&mut self, the_type: u8) {
        // OCCT: iterates sub-shapes of given type on both shapes,
        // checks for duplicates (VERTEX: distance vs tolerance sum;
        // EDGE: IntTools_EdgeEdge for TopAbs_EDGE coincidence;
        // FACE: not implemented in OCCT).
        if self.my_shape1 == usize::MAX || self.my_shape2 == usize::MAX { return; }
        if self.my_empty1 || self.my_empty2 { return; }
        // rcad: sub-shape iteration and per-type comparison not yet wired.
    }

    fn test_merge_vertex(&mut self) { self.test_merge_sub_shapes(0); }
    fn test_merge_edge(&mut self) { self.test_merge_sub_shapes(1); }

    // ── TestContinuity L896-958 ────────────────────────────────────
    fn test_continuity(&mut self) {
        // OCCT: iterates edges/faces, checks GeomAbs_C0 continuity
        // via curve->Continuity() / surface->Continuity().
        // rcad: continuity queries not yet available on curve/surface types.
    }

    // ── TestCurveOnSurface L962-1015 ────────────────────────────────
    fn test_curve_on_surface(&mut self) {
        // OCCT: iterates face edges, calls BOPTools_AlgoTools::ComputeTolerance.
        // If deviation > edge tolerance, reports InvalidCurveOnSurface.
        // rcad: ComputeTolerance not yet available.
    }

    // ── Helpers ────────────────────────────────────────────────────
    fn push_result(&mut self, status: CheckStatus) {
        self.my_result.push(CheckResult { check_status: status });
    }

    fn has_errors_or_stop(&self) -> bool {
        self.my_stop_on_first && !self.my_result.is_empty()
    }

    fn should_run(&self) -> bool {
        self.my_result.is_empty() || !self.my_stop_on_first
    }
}

impl Default for ArgumentAnalyzer {
    fn default() -> Self { Self::new() }
}
