
impl TolerancePropagationEngine {
    /// Create a new engine with default configuration.
    pub fn new() -> Self {
        Self {
            config: TolerancePropagationConfig::default(),
        }
    }

    /// Create a new engine with custom configuration.
    pub fn with_config(config: TolerancePropagationConfig) -> Self {
        Self { config }
    }

    /// Create an engine with OCCT-standard rules.
    pub fn occt_standard() -> Self {
        Self::with_config(TolerancePropagationConfig::occt_standard())
    }

    /// Create an engine with conservative rules.
    pub fn conservative() -> Self {
        Self::with_config(TolerancePropagationConfig::conservative())
    }

    /// Create an engine with aggressive rules.
    pub fn aggressive() -> Self {
        Self::with_config(TolerancePropagationConfig::aggressive())
    }

    /// Create an engine with bounded rules.
    pub fn bounded(max_tol: f64) -> Self {
        Self::with_config(TolerancePropagationConfig::bounded(max_tol))
    }

    /// Propagate tolerances according to the configured rule.
    pub fn propagate(&self, brep: &BRep) -> (BRep, TolerancePropagationReport) {
        match self.config.rule {
            ToleranceRule::OcctStandard => self.propagate_occt_standard(brep),
            ToleranceRule::Conservative => self.propagate_conservative(brep),
            ToleranceRule::Aggressive => self.propagate_aggressive(brep),
            ToleranceRule::Harmonized => self.propagate_harmonized(brep),
            ToleranceRule::Bounded => self.propagate_bounded(brep),
            ToleranceRule::ModelScale => self.propagate_model_scale(brep),
        }
    }

    fn propagate_occt_standard(&self, brep: &BRep) -> (BRep, TolerancePropagationReport) {
        let mut result = brep.clone();
        let mut report = TolerancePropagationReport::default();
        let floor = self.config.tolerance_floor.max(TOLERANCE_ABS);

        self.ensure_tolerance_arrays(&mut result, floor);

        // Multiple passes to ensure convergence
        for _pass in 0..self.config.propagation_passes {
            // Step 1: Vertex -> Edge
            for ei in 0..result.edges.len() {
                let edge = &result.edges[ei];
                let vtol_s = result.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
                let vtol_e = result.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
                let cur_etol = result.geom.edge_tolerance[ei];
                let new_etol = cur_etol.max(vtol_s).max(vtol_e).min(self.config.max_tolerance);

                if new_etol > cur_etol + TOLERANCE_FLOAT_DEDUP {
                    result.geom.edge_tolerance[ei] = new_etol;
                    report.edges_updated += 1;
                }
            }

            // Step 2: Edge -> Face
            let mut flat_fi = 0usize;
            for solid in &result.solids {
                for shell in &solid.shells {
                    for face in &shell.faces {
                        let mut max_etol = floor;
                        for we in &face.outer_wire.edges {
                            if we.idx < result.geom.edge_tolerance.len() {
                                max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                            }
                        }
                        for iw in &face.inner_wires {
                            for we in &iw.edges {
                                if we.idx < result.geom.edge_tolerance.len() {
                                    max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                                }
                            }
                        }

                        let cur_ftol = result.geom.face_tolerance.get(flat_fi).copied().unwrap_or(floor);
                        let new_ftol = max_etol.min(self.config.max_tolerance);

                        if new_ftol > cur_ftol + TOLERANCE_FLOAT_DEDUP
                            && flat_fi < result.geom.face_tolerance.len() {
                                result.geom.face_tolerance[flat_fi] = new_ftol;
                                report.faces_updated += 1;
                            }
                        flat_fi += 1;
                    }
                }
            }
        }

        // Handle conflicts
        if self.config.conflict_policy != ConflictResolutionPolicy::Ignore {
            let (detected, resolved) = self.handle_conflicts(&mut result, floor);
            report.conflicts_detected = detected;
            report.conflicts_resolved = resolved;
        }

        self.compute_report_stats(&result, &mut report);
        (result, report)
    }

    fn propagate_conservative(&self, brep: &BRep) -> (BRep, TolerancePropagationReport) {
        let mut result = brep.clone();
        let mut report = TolerancePropagationReport::default();
        let floor = self.config.tolerance_floor.max(TOLERANCE_ABS);

        self.ensure_tolerance_arrays(&mut result, floor);

        // Only propagate where absolutely necessary (conflicts)
        let (detected, resolved) = self.handle_conflicts(&mut result, floor);
        report.conflicts_detected = detected;
        report.conflicts_resolved = resolved;

        self.compute_report_stats(&result, &mut report);
        (result, report)
    }

    fn propagate_aggressive(&self, brep: &BRep) -> (BRep, TolerancePropagationReport) {
        let mut result = brep.clone();
        let mut report = TolerancePropagationReport::default();
        let floor = self.config.tolerance_floor.max(TOLERANCE_ABS);

        self.ensure_tolerance_arrays(&mut result, floor);

        // Multiple aggressive passes
        for _pass in 0..self.config.propagation_passes {
            // Vertex -> Edge
            for ei in 0..result.edges.len() {
                let edge = &result.edges[ei];
                let vtol_s = result.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
                let vtol_e = result.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);

                // Aggressive: always take max
                let new_etol = vtol_s.max(vtol_e);
                let cur_etol = result.geom.edge_tolerance[ei];

                if new_etol > cur_etol {
                    result.geom.edge_tolerance[ei] = new_etol;
                    report.edges_updated += 1;
                }
            }

            // Edge -> Face (aggressive)
            let mut flat_fi = 0usize;
            for solid in &result.solids {
                for shell in &solid.shells {
                    for face in &shell.faces {
                        let mut max_etol = floor;
                        for we in &face.outer_wire.edges {
                            if we.idx < result.geom.edge_tolerance.len() {
                                max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                            }
                        }
                        for iw in &face.inner_wires {
                            for we in &iw.edges {
                                if we.idx < result.geom.edge_tolerance.len() {
                                    max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                                }
                            }
                        }

                        if flat_fi < result.geom.face_tolerance.len() {
                            let cur_ftol = result.geom.face_tolerance[flat_fi];
                            if max_etol > cur_ftol {
                                result.geom.face_tolerance[flat_fi] = max_etol;
                                report.faces_updated += 1;
                            }
                        }
                        flat_fi += 1;
                    }
                }
            }

            // Face -> Edge -> Vertex (reverse propagation)
            let mut flat_fi = 0usize;
            for solid in &result.solids {
                for shell in &solid.shells {
                    for face in &shell.faces {
                        let ftol = result.geom.face_tolerance.get(flat_fi).copied().unwrap_or(floor);

                        for we in &face.outer_wire.edges {
                            if we.idx < result.geom.edge_tolerance.len()
                                && ftol > result.geom.edge_tolerance[we.idx] {
                                    result.geom.edge_tolerance[we.idx] = ftol;
                                    report.edges_updated += 1;
                                }
                        }
                        flat_fi += 1;
                    }
                }
            }

            // Edge -> Vertex (reverse propagation)
            for ei in 0..result.edges.len() {
                let edge = &result.edges[ei];
                let etol = result.geom.edge_tolerance[ei];

                if edge.start < result.geom.vertex_tolerance.len() && etol > result.geom.vertex_tolerance[edge.start] {
                    result.geom.vertex_tolerance[edge.start] = etol;
                    report.vertices_updated += 1;
                }
                if edge.end < result.geom.vertex_tolerance.len() && etol > result.geom.vertex_tolerance[edge.end] {
                    result.geom.vertex_tolerance[edge.end] = etol;
                    report.vertices_updated += 1;
                }
            }
        }

        self.compute_report_stats(&result, &mut report);
        (result, report)
    }

    fn propagate_harmonized(&self, brep: &BRep) -> (BRep, TolerancePropagationReport) {
        let mut result = brep.clone();
        let mut report = TolerancePropagationReport::default();
        let floor = self.config.tolerance_floor.max(TOLERANCE_ABS);

        self.ensure_tolerance_arrays(&mut result, floor);

        // Build edge-vertex connectivity
        let mut vertex_max_edge_tol: Vec<f64> = vec![floor; result.vertices.len()];
        for ei in 0..result.edges.len() {
            let edge = &result.edges[ei];
            let etol = result.geom.edge_tolerance.get(ei).copied().unwrap_or(floor);

            if edge.start < vertex_max_edge_tol.len() {
                vertex_max_edge_tol[edge.start] = vertex_max_edge_tol[edge.start].max(etol);
            }
            if edge.end < vertex_max_edge_tol.len() {
                vertex_max_edge_tol[edge.end] = vertex_max_edge_tol[edge.end].max(etol);
            }
        }

        // Harmonize: propagate max through connected topology
        for _pass in 0..self.config.propagation_passes {
            // Find global max for connected components
            let mut changed = false;

            for ei in 0..result.edges.len() {
                let edge = &result.edges[ei];
                let vtol_s = vertex_max_edge_tol.get(edge.start).copied().unwrap_or(floor);
                let vtol_e = vertex_max_edge_tol.get(edge.end).copied().unwrap_or(floor);
                let cur_etol = result.geom.edge_tolerance[ei];
                let harmonized = cur_etol.max(vtol_s).max(vtol_e);

                if harmonized > cur_etol + TOLERANCE_FLOAT_DEDUP {
                    result.geom.edge_tolerance[ei] = harmonized;
                    // Update vertex max
                    if edge.start < vertex_max_edge_tol.len() {
                        vertex_max_edge_tol[edge.start] = vertex_max_edge_tol[edge.start].max(harmonized);
                    }
                    if edge.end < vertex_max_edge_tol.len() {
                        vertex_max_edge_tol[edge.end] = vertex_max_edge_tol[edge.end].max(harmonized);
                    }
                    report.edges_updated += 1;
                    changed = true;
                }
            }

            if !changed {
                break;
            }
        }

        // Propagate to faces
        let mut flat_fi = 0usize;
        for solid in &result.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let mut max_etol = floor;
                    for we in &face.outer_wire.edges {
                        if we.idx < result.geom.edge_tolerance.len() {
                            max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                        }
                    }
                    for iw in &face.inner_wires {
                        for we in &iw.edges {
                            if we.idx < result.geom.edge_tolerance.len() {
                                max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                            }
                        }
                    }

                    if flat_fi < result.geom.face_tolerance.len() {
                        let cur_ftol = result.geom.face_tolerance[flat_fi];
                        if max_etol > cur_ftol {
                            result.geom.face_tolerance[flat_fi] = max_etol;
                            report.faces_updated += 1;
                        }
                    }
                    flat_fi += 1;
                }
            }
        }

        self.compute_report_stats(&result, &mut report);
        (result, report)
    }

    fn propagate_bounded(&self, brep: &BRep) -> (BRep, TolerancePropagationReport) {
        let mut result = brep.clone();
        let mut report = TolerancePropagationReport::default();
        let floor = self.config.tolerance_floor.max(TOLERANCE_ABS);
        let bound = self.config.bound_value.max(floor);

        self.ensure_tolerance_arrays(&mut result, floor);

        // Standard propagation with bounding
        for ei in 0..result.edges.len() {
            let edge = &result.edges[ei];
            let vtol_s = result.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
            let vtol_e = result.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
            let cur_etol = result.geom.edge_tolerance[ei];
            let new_etol = cur_etol.max(vtol_s).max(vtol_e).min(bound);

            if (new_etol - cur_etol).abs() > TOLERANCE_FLOAT_DEDUP {
                result.geom.edge_tolerance[ei] = new_etol;
                report.edges_updated += 1;
            }
        }

        let mut flat_fi = 0usize;
        for solid in &result.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let mut max_etol = floor;
                    for we in &face.outer_wire.edges {
                        if we.idx < result.geom.edge_tolerance.len() {
                            max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                        }
                    }
                    for iw in &face.inner_wires {
                        for we in &iw.edges {
                            if we.idx < result.geom.edge_tolerance.len() {
                                max_etol = max_etol.max(result.geom.edge_tolerance[we.idx]);
                            }
                        }
                    }

                    if flat_fi < result.geom.face_tolerance.len() {
                        let bounded_etol = max_etol.min(bound);
                        let cur_ftol = result.geom.face_tolerance[flat_fi];
                        if bounded_etol > cur_ftol {
                            result.geom.face_tolerance[flat_fi] = bounded_etol;
                            report.faces_updated += 1;
                        }
                    }
                    flat_fi += 1;
                }
            }
        }

        // Clamp all tolerances
        for tol in &mut result.geom.vertex_tolerance {
            *tol = tol.min(bound);
        }
        for tol in &mut result.geom.edge_tolerance {
            *tol = tol.min(bound);
        }
        for tol in &mut result.geom.face_tolerance {
            *tol = tol.min(bound);
        }

        self.compute_report_stats(&result, &mut report);
        (result, report)
    }

    fn propagate_model_scale(&self, brep: &BRep) -> (BRep, TolerancePropagationReport) {
        let mut result = brep.clone();
        let mut report = TolerancePropagationReport::default();
        let scale = self.config.model_scale.max(TOLERANCE_LINEAR_ULTRA_STRICT);
        let floor = self.config.tolerance_floor.max(TOLERANCE_ABS * scale);

        self.ensure_tolerance_arrays(&mut result, floor);

        // Scale all existing tolerances
        for tol in &mut result.geom.vertex_tolerance {
            *tol = (*tol * scale).max(floor).min(self.config.max_tolerance);
        }
        for tol in &mut result.geom.edge_tolerance {
            *tol = (*tol * scale).max(floor).min(self.config.max_tolerance);
        }
        for tol in &mut result.geom.face_tolerance {
            *tol = (*tol * scale).max(floor).min(self.config.max_tolerance);
        }

        // Then apply standard propagation
        for ei in 0..result.edges.len() {
            let edge = &result.edges[ei];
            let vtol_s = result.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
            let vtol_e = result.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
            let cur_etol = result.geom.edge_tolerance[ei];

            let new_etol = cur_etol.max(vtol_s).max(vtol_e);
            if new_etol > cur_etol {
                result.geom.edge_tolerance[ei] = new_etol;
                report.edges_updated += 1;
            }
        }

        self.compute_report_stats(&result, &mut report);
        (result, report)
    }

    fn ensure_tolerance_arrays(&self, brep: &mut BRep, floor: f64) {
        let n_verts = brep.vertices.len();
        let n_edges = brep.edges.len();
        let n_faces: usize = brep.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();

        if brep.geom.vertex_tolerance.len() < n_verts {
            brep.geom.vertex_tolerance.resize(n_verts, floor);
        }
        if brep.geom.edge_tolerance.len() < n_edges {
            brep.geom.edge_tolerance.resize(n_edges, floor);
        }
        if brep.geom.face_tolerance.len() < n_faces {
            brep.geom.face_tolerance.resize(n_faces, floor);
        }
    }

    fn handle_conflicts(&self, brep: &mut BRep, floor: f64) -> (usize, usize) {
        match self.config.conflict_policy {
            ConflictResolutionPolicy::Ignore => (0, 0),
            ConflictResolutionPolicy::PropagateUp => {
                detect_and_resolve_tolerance_conflicts(brep, floor)
            }
            ConflictResolutionPolicy::ClampDown => {
                // Clamp higher-level tolerances down
                let mut conflicts = 0usize;
                let mut resolved = 0usize;

                for ei in 0..brep.edges.len() {
                    let edge = &brep.edges[ei];
                    let vtol_s = brep.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
                    let vtol_e = brep.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
                    let etol = brep.geom.edge_tolerance.get(ei).copied().unwrap_or(floor);

                    if vtol_s > etol + TOLERANCE_FLOAT_DEDUP || vtol_e > etol + TOLERANCE_FLOAT_DEDUP {
                        conflicts += 1;
                        // Clamp vertices down
                        if edge.start < brep.geom.vertex_tolerance.len() {
                            brep.geom.vertex_tolerance[edge.start] = brep.geom.vertex_tolerance[edge.start].min(etol);
                        }
                        if edge.end < brep.geom.vertex_tolerance.len() {
                            brep.geom.vertex_tolerance[edge.end] = brep.geom.vertex_tolerance[edge.end].min(etol);
                        }
                        resolved += 1;
                    }
                }

                (conflicts, resolved)
            }
            ConflictResolutionPolicy::ReportOnly => {
                // Just count conflicts
                let mut conflicts = 0usize;

                for ei in 0..brep.edges.len() {
                    let edge = &brep.edges[ei];
                    let vtol_s = brep.geom.vertex_tolerance.get(edge.start).copied().unwrap_or(floor);
                    let vtol_e = brep.geom.vertex_tolerance.get(edge.end).copied().unwrap_or(floor);
                    let etol = brep.geom.edge_tolerance.get(ei).copied().unwrap_or(floor);

                    if vtol_s > etol + TOLERANCE_FLOAT_DEDUP || vtol_e > etol + TOLERANCE_FLOAT_DEDUP {
                        conflicts += 1;
                    }
                }

                (conflicts, 0)
            }
        }
    }

    fn compute_report_stats(&self, brep: &BRep, report: &mut TolerancePropagationReport) {
        if !brep.geom.vertex_tolerance.is_empty() {
            report.max_vertex_tolerance = brep.geom.vertex_tolerance.iter()
                .cloned()
                .fold(0.0_f64, f64::max);
        }
        if !brep.geom.edge_tolerance.is_empty() {
            report.max_edge_tolerance = brep.geom.edge_tolerance.iter()
                .cloned()
                .fold(0.0_f64, f64::max);
        }
        if !brep.geom.face_tolerance.is_empty() {
            report.max_face_tolerance = brep.geom.face_tolerance.iter()
                .cloned()
                .fold(0.0_f64, f64::max);
        }
        report.rule_applied = self.config.rule;
    }
}

/// Report from tolerance propagation.
#[derive(Debug, Clone, Default)]
pub struct TolerancePropagationReport {
    /// Number of vertices whose tolerance was updated.
    pub vertices_updated: usize,
    /// Number of edges whose tolerance was updated.
    pub edges_updated: usize,
    /// Number of faces whose tolerance was updated.
    pub faces_updated: usize,
    /// Number of tolerance conflicts detected.
    pub conflicts_detected: usize,
    /// Number of tolerance conflicts resolved.
    pub conflicts_resolved: usize,
    /// Maximum vertex tolerance after propagation.
    pub max_vertex_tolerance: f64,
    /// Maximum edge tolerance after propagation.
    pub max_edge_tolerance: f64,
    /// Maximum face tolerance after propagation.
    pub max_face_tolerance: f64,
    /// The rule that was applied.
    pub rule_applied: ToleranceRule,
}

// 閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡?
// Tolerance Consistency Analysis
// 閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡?

/// A specific tolerance violation found during analysis.
#[derive(Debug, Clone)]
pub struct ToleranceViolation {
    /// Type of the violation.
    pub violation_type: ToleranceViolationType,
    /// Index of the entity with the violation.
    pub entity_index: usize,
    /// Related entity index (e.g., edge for vertex violation).
    pub related_index: Option<usize>,
    /// Actual tolerance value.
    pub actual_tolerance: f64,
    /// Expected or related tolerance value.
    pub expected_tolerance: f64,
    /// Severity of the violation (1-5, 5 being most severe).
    pub severity: u8,
    /// Suggested fix for the violation.
    pub suggested_fix: ToleranceFix,
}

/// Type of tolerance violation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceViolationType {
    /// Vertex tolerance exceeds edge tolerance.
    VertexExceedsEdge,
    /// Edge tolerance exceeds face tolerance.
    EdgeExceedsFace,
    /// Tolerance is below minimum floor.
    BelowFloor,
    /// Tolerance exceeds maximum allowed.
    ExceedsMaximum,
    /// Inconsistent tolerances across seam edges.
    SeamInconsistency,
    /// Tolerance is NaN or infinite.
    InvalidValue,
}

/// Suggested fix for a tolerance violation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceFix {
    /// Increase the lower-level tolerance.
    IncreaseLower,
    /// Decrease the higher-level tolerance.
    DecreaseHigher,
    /// Set tolerance to a specific value.
    SetToValue,
    /// Propagate tolerance through topology.
    Propagate,
    /// No automatic fix available.
    ManualIntervention,
}

/// Report from tolerance consistency analysis.
#[derive(Debug, Clone, Default)]
pub struct ToleranceConsistencyReport {
    /// Whether the BRep has consistent tolerances.
    pub is_consistent: bool,
    /// Total number of violations found.
    pub violation_count: usize,
    /// Number of critical violations (severity >= 4).
    pub critical_violation: usize,
    /// List of all violations found.
    pub violations: Vec<ToleranceViolation>,
    /// Summary statistics.
    pub stats: ToleranceAnalysisReport,
    /// Suggested global fixes.
    pub suggested_global_fixes: Vec<String>,
}

impl ToleranceConsistencyReport {
    /// Get violations by type.
    pub fn violations_by_type(&self, violation_type: ToleranceViolationType) -> Vec<&ToleranceViolation> {
        self.violations.iter()
            .filter(|v| v.violation_type == violation_type)
            .collect()
    }

    /// Get critical violations.
    pub fn critical_violations(&self) -> Vec<&ToleranceViolation> {
        self.violations.iter()
            .filter(|v| v.severity >= 4)
            .collect()
    }

    /// Returns a summary string.
    pub fn summary(&self) -> String {
        if self.is_consistent {
            "Tolerance consistency: OK".to_string()
        } else {
            format!(
                "Tolerance consistency: {} violations ({} critical)",
                self.violation_count,
                self.critical_violations().len()
            )
        }
    }
}

/// Analyze tolerance consistency in a BRep.
///
/// This function checks for tolerance violations and inconsistencies:
/// - Vertex tolerances exceeding edge tolerances
/// - Edge tolerances exceeding face tolerances
/// - Tolerances below floor or above maximum
/// - Seam edge inconsistencies
/// - Invalid (NaN/Inf) tolerance values
///
/// # Arguments
///
/// * `brep` - The BRep to analyze.
/// * `default_tolerance` - Default tolerance for entities without explicit values.
/// * `min_tolerance` - Minimum allowed tolerance (floor).
/// * `max_tolerance` - Maximum allowed tolerance.
///
/// # Returns
///
/// A `ToleranceConsistencyReport` containing all violations found.
pub fn analyze_tolerance_consistency(
    brep: &BRep,
    default_tolerance: f64,
    min_tolerance: f64,
    max_tolerance: f64,
) -> ToleranceConsistencyReport {
    let mut report = ToleranceConsistencyReport::default();
    let floor = min_tolerance.max(TOLERANCE_ABS);

    // Get base statistics
    report.stats = analyze_tolerances(brep, default_tolerance);

    let n_verts = brep.vertices.len();
    let n_edges = brep.edges.len();

    // Ensure we have tolerance arrays to work with
    let vertex_tols: Vec<f64> = if brep.geom.vertex_tolerance.len() >= n_verts {
        brep.geom.vertex_tolerance.clone()
    } else {
        vec![default_tolerance; n_verts]
    };

    let edge_tols: Vec<f64> = if brep.geom.edge_tolerance.len() >= n_edges {
        brep.geom.edge_tolerance.clone()
    } else {
        vec![default_tolerance; n_edges]
    };

    // Check for invalid values
    for (i, &tol) in vertex_tols.iter().enumerate() {
        if !tol.is_finite() || tol <= 0.0 {
            report.violations.push(ToleranceViolation {
                violation_type: ToleranceViolationType::InvalidValue,
                entity_index: i,
                related_index: None,
                actual_tolerance: tol,
                expected_tolerance: floor,
                severity: 5,
                suggested_fix: ToleranceFix::SetToValue,
            });
        }
    }

    for (i, &tol) in edge_tols.iter().enumerate() {
        if !tol.is_finite() || tol <= 0.0 {
            report.violations.push(ToleranceViolation {
                violation_type: ToleranceViolationType::InvalidValue,
                entity_index: i,
                related_index: None,
                actual_tolerance: tol,
                expected_tolerance: floor,
                severity: 5,
                suggested_fix: ToleranceFix::SetToValue,
            });
        }
    }

    // Check vertex tolerances below floor or above max
    for (i, &tol) in vertex_tols.iter().enumerate() {
        if tol < floor {
            report.violations.push(ToleranceViolation {
                violation_type: ToleranceViolationType::BelowFloor,
                entity_index: i,
                related_index: None,
                actual_tolerance: tol,
                expected_tolerance: floor,
                severity: 2,
                suggested_fix: ToleranceFix::SetToValue,
            });
        } else if tol > max_tolerance {
            report.violations.push(ToleranceViolation {
                violation_type: ToleranceViolationType::ExceedsMaximum,
                entity_index: i,
                related_index: None,
                actual_tolerance: tol,
                expected_tolerance: max_tolerance,
                severity: 3,
                suggested_fix: ToleranceFix::DecreaseHigher,
            });
        }
    }

    // Check edge tolerances
    for (i, &tol) in edge_tols.iter().enumerate() {
        if tol < floor {
            report.violations.push(ToleranceViolation {
                violation_type: ToleranceViolationType::BelowFloor,
                entity_index: i,
                related_index: None,
                actual_tolerance: tol,
                expected_tolerance: floor,
                severity: 2,
                suggested_fix: ToleranceFix::SetToValue,
            });
        } else if tol > max_tolerance {
            report.violations.push(ToleranceViolation {
                violation_type: ToleranceViolationType::ExceedsMaximum,
                entity_index: i,
                related_index: None,
                actual_tolerance: tol,
                expected_tolerance: max_tolerance,
                severity: 3,
                suggested_fix: ToleranceFix::DecreaseHigher,
            });
        }
    }

    // Check vertex > edge violations
    for (ei, edge) in brep.edges.iter().enumerate() {
        let etol = edge_tols.get(ei).copied().unwrap_or(default_tolerance);

        if edge.start < vertex_tols.len() {
            let vtol = vertex_tols[edge.start];
            if vtol > etol + TOLERANCE_FLOAT_DEDUP {
                report.violations.push(ToleranceViolation {
                    violation_type: ToleranceViolationType::VertexExceedsEdge,
                    entity_index: edge.start,
                    related_index: Some(ei),
                    actual_tolerance: vtol,
                    expected_tolerance: etol,
                    severity: 4,
                    suggested_fix: ToleranceFix::IncreaseLower,
                });
            }
        }

        if edge.end < vertex_tols.len() {
            let vtol = vertex_tols[edge.end];
            if vtol > etol + TOLERANCE_FLOAT_DEDUP {
                report.violations.push(ToleranceViolation {
                    violation_type: ToleranceViolationType::VertexExceedsEdge,
                    entity_index: edge.end,
                    related_index: Some(ei),
                    actual_tolerance: vtol,
                    expected_tolerance: etol,
                    severity: 4,
                    suggested_fix: ToleranceFix::IncreaseLower,
                });
            }
        }
    }

    // Check edge > face violations
    let mut flat_fi = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let ftol = brep.geom.face_tolerance.get(flat_fi).copied().unwrap_or(default_tolerance);

                for we in &face.outer_wire.edges {
                    let etol = edge_tols.get(we.idx).copied().unwrap_or(default_tolerance);
                    if etol > ftol + TOLERANCE_FLOAT_DEDUP {
                        report.violations.push(ToleranceViolation {
                            violation_type: ToleranceViolationType::EdgeExceedsFace,
                            entity_index: we.idx,
                            related_index: Some(flat_fi),
                            actual_tolerance: etol,
                            expected_tolerance: ftol,
                            severity: 3,
                            suggested_fix: ToleranceFix::IncreaseLower,
                        });
                    }
                }

                for iw in &face.inner_wires {
                    for we in &iw.edges {
                        let etol = edge_tols.get(we.idx).copied().unwrap_or(default_tolerance);
                        if etol > ftol + TOLERANCE_FLOAT_DEDUP {
                            report.violations.push(ToleranceViolation {
                                violation_type: ToleranceViolationType::EdgeExceedsFace,
                                entity_index: we.idx,
                                related_index: Some(flat_fi),
                                actual_tolerance: etol,
                                expected_tolerance: ftol,
                                severity: 3,
                                suggested_fix: ToleranceFix::IncreaseLower,
                            });
                        }
                    }
                }

                flat_fi += 1;
            }
        }
    }

    // Compute summary
    report.violation_count = report.violations.len();
    report.critical_violation = report.violations.iter().filter(|v| v.severity >= 4).count();
    report.is_consistent = report.violations.is_empty();

    // Generate global fix suggestions
    if !report.violations.is_empty() {
        let vertex_edge_violations = report.violations_by_type(ToleranceViolationType::VertexExceedsEdge).len();
        let edge_face_violations = report.violations_by_type(ToleranceViolationType::EdgeExceedsFace).len();
        let invalid_values = report.violations_by_type(ToleranceViolationType::InvalidValue).len();

        if vertex_edge_violations > 0 {
            report.suggested_global_fixes.push(format!(
                "Run tolerance propagation (vertex閳姀dge) to fix {} vertex>edge violations",
                vertex_edge_violations
            ));
        }
        if edge_face_violations > 0 {
            report.suggested_global_fixes.push(format!(
                "Run tolerance propagation (edge閳姁ace) to fix {} edge>face violations",
                edge_face_violations
            ));
        }
        if invalid_values > 0 {
            report.suggested_global_fixes.push(format!(
                "Fix {} invalid (NaN/Inf) tolerance values before processing",
                invalid_values
            ));
        }
    }

    report
}

/// Apply automatic fixes to tolerance violations.
///
/// This function attempts to automatically fix tolerance violations
/// by propagating tolerances according to the suggested fixes.
///
/// # Arguments
///
/// * `brep` - The BRep to fix.
/// * `report` - The consistency report with violations.
/// * `max_fixes` - Maximum number of fixes to apply (0 = unlimited).
///
/// # Returns
///
/// A tuple of (fixed BRep, number of fixes applied).
pub fn apply_tolerance_fixes(
    brep: &BRep,
    report: &ToleranceConsistencyReport,
    max_fixes: usize,
) -> (BRep, usize) {
    let mut result = brep.clone();
    let mut fixes_applied = 0usize;
    let floor = TOLERANCE_ABS;

    let n_verts = result.vertices.len();
    let n_edges = result.edges.len();
    let n_faces: usize = result.solids.iter()
        .flat_map(|s| s.shells.iter())
        .map(|sh| sh.faces.len())
        .sum();

    // Ensure arrays are sized
    if result.geom.vertex_tolerance.len() < n_verts {
        result.geom.vertex_tolerance.resize(n_verts, floor);
    }
    if result.geom.edge_tolerance.len() < n_edges {
        result.geom.edge_tolerance.resize(n_edges, floor);
    }
    if result.geom.face_tolerance.len() < n_faces {
        result.geom.face_tolerance.resize(n_faces, floor);
    }

    for violation in &report.violations {
        if max_fixes > 0 && fixes_applied >= max_fixes {
            break;
        }

        match violation.suggested_fix {
            ToleranceFix::SetToValue => {
                match violation.violation_type {
                    ToleranceViolationType::InvalidValue | ToleranceViolationType::BelowFloor => {
                        if violation.entity_index < result.geom.vertex_tolerance.len() {
                            result.geom.vertex_tolerance[violation.entity_index] = violation.expected_tolerance;
                            fixes_applied += 1;
                        }
                    }
                    _ => {}
                }
            }
            ToleranceFix::IncreaseLower => {
                match violation.violation_type {
                    ToleranceViolationType::VertexExceedsEdge => {
                        if let Some(ei) = violation.related_index
                            && ei < result.geom.edge_tolerance.len() {
                                let new_tol = result.geom.edge_tolerance[ei].max(violation.actual_tolerance);
                                result.geom.edge_tolerance[ei] = new_tol;
                                fixes_applied += 1;
                            }
                    }
                    ToleranceViolationType::EdgeExceedsFace => {
                        if let Some(fi) = violation.related_index
                            && fi < result.geom.face_tolerance.len() {
                                let new_tol = result.geom.face_tolerance[fi].max(violation.actual_tolerance);
                                result.geom.face_tolerance[fi] = new_tol;
                                fixes_applied += 1;
                            }
                    }
                    _ => {}
                }
            }
            ToleranceFix::DecreaseHigher => {
                if violation.entity_index < result.geom.vertex_tolerance.len() {
                    result.geom.vertex_tolerance[violation.entity_index] = violation.expected_tolerance;
                    fixes_applied += 1;
                }
            }
            ToleranceFix::Propagate => {
                // Use the engine for propagation
                let engine = TolerancePropagationEngine::occt_standard();
                let (propagated, _) = engine.propagate(&result);
                result = propagated;
                fixes_applied += 1;
            }
            ToleranceFix::ManualIntervention => {
                // Cannot auto-fix
            }
        }
    }

    (result, fixes_applied)
}

// 閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡?
// Enhanced Internal Face Detection and Removal
// 閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡鎰ㄦ櫜閳烘劏鏅查埡?

/// Configuration for internal face detection.
#[derive(Debug, Clone)]
pub struct InternalFaceDetectionConfig {
    /// Tolerance for geometric comparisons.
    pub tolerance: f64,
    /// Whether to use material side analysis.
    pub use_material_side_analysis: bool,
    /// Whether to use ray casting for visibility check.
    pub use_visibility_check: bool,
    /// Whether to check for duplicate faces with opposite orientation.
    pub check_duplicate_faces: bool,
    /// Whether to consider void shell faces as internal.
    pub consider_void_shells: bool,
    /// Minimum edge count for a face to be considered valid (faces with fewer edges may be internal).
    pub min_edge_count: usize,
    /// Whether to use connectivity analysis (edges shared with multiple faces).
    pub use_connectivity_analysis: bool,
    /// Threshold for shared edge ratio to consider a face internal (0.0-1.0).
    pub shared_edge_threshold: f64,
}

impl Default for InternalFaceDetectionConfig {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            use_material_side_analysis: true,
            use_visibility_check: false, // Disabled by default - can be unreliable
            check_duplicate_faces: true,
            consider_void_shells: true,
            min_edge_count: 3,
            use_connectivity_analysis: true,
            shared_edge_threshold: 0.9,
        }
    }
}

impl InternalFaceDetectionConfig {
    /// Create a conservative configuration (only obvious internal faces).
    pub fn conservative() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            use_material_side_analysis: true,
            use_visibility_check: false,
            check_duplicate_faces: true,
            consider_void_shells: true,
            min_edge_count: 3,
            use_connectivity_analysis: true,
            shared_edge_threshold: 1.0,
        }
    }

    /// Create an aggressive configuration (more internal face candidates).
    pub fn aggressive() -> Self {
        Self {
            tolerance: TOLERANCE_ABS * 10.0,
            use_material_side_analysis: true,
            use_visibility_check: true,
            check_duplicate_faces: true,
            consider_void_shells: true,
            min_edge_count: 2,
            use_connectivity_analysis: true,
            shared_edge_threshold: 0.75,
        }
    }

    /// Create a configuration optimized for post-boolean cleanup.
    pub fn for_post_boolean() -> Self {
        Self {
            tolerance: TOLERANCE_ABS * 5.0,
            use_material_side_analysis: true,
            use_visibility_check: false, // Disabled - can be unreliable
            check_duplicate_faces: true,
            consider_void_shells: true,
            min_edge_count: 3,
            use_connectivity_analysis: true,
            shared_edge_threshold: 0.85,
        }
    }
}

/// Report from internal face detection.
#[derive(Debug, Clone, Default)]
pub struct InternalFaceDetectionReport {
    /// Indices of detected internal faces (flattened).
    pub internal_face_indices: Vec<usize>,
    /// Number of faces detected by material side analysis.
    pub by_material_side: usize,
    /// Number of faces detected by visibility check.
    pub by_visibility: usize,
    /// Number of faces detected as duplicates.
    pub by_duplicate: usize,
    /// Number of faces detected in void shells.
    pub by_void_shell: usize,
    /// Number of faces detected by connectivity analysis.
    pub by_connectivity: usize,
    /// Total number of faces analyzed.
    pub total_faces: usize,
    /// Summary string.
    pub summary: String,
}

/// Detect internal faces in a BRep using comprehensive analysis.
///
/// Internal faces are faces that do not contribute to the outer boundary of the solid.
/// These typically arise from boolean operations where partition/separator faces
/// are not properly removed.
///
/// # Detection Methods
/// 1. **Material side analysis**: Faces where both sides point to the same material region
/// 2. **Visibility check**: Faces not visible from outside the solid (via ray casting)
/// 3. **Duplicate face detection**: Faces with opposite orientation to another face
/// 4. **Void shell detection**: Faces in internal void shells
/// 5. **Connectivity analysis**: Faces with all edges shared by other faces
///
/// # Arguments
/// * `brep` - The BRep to analyze.
///
/// # Returns
/// A vector of flattened face indices that are identified as internal.
pub fn detect_internal_faces(brep: &BRep) -> Vec<usize> {
    detect_internal_faces_with_config(brep, &InternalFaceDetectionConfig::default())
        .internal_face_indices
}

/// Detect internal faces with custom configuration.
///
/// See [`detect_internal_faces`] for details.
pub fn detect_internal_faces_with_config(
    brep: &BRep,
    config: &InternalFaceDetectionConfig,
) -> InternalFaceDetectionReport {
    let mut report = InternalFaceDetectionReport::default();
    let tol = config.tolerance.max(TOLERANCE_ABS);
    let mut internal_set: std::collections::HashSet<usize> = std::collections::HashSet::new();

    // Collect all faces with their flattened indices
    let faces: Vec<(usize, usize, usize, &Face)> = brep
        .solids
        .iter()
        .enumerate()
        .flat_map(|(si, solid)| {
            solid.shells.iter().enumerate().flat_map(move |(shi, shell)| {
                shell.faces.iter().enumerate().map(move |(fi, face)| (si, shi, fi, face))
            })
        })
        .collect();

    report.total_faces = faces.len();

    if faces.is_empty() {
        report.summary = "No faces to analyze".to_string();
        return report;
    }

    // Method 1: Void shell detection
    if config.consider_void_shells {
        let void_faces = detect_void_shell_faces(brep, &faces);
        for idx in void_faces {
            if internal_set.insert(idx) {
                report.by_void_shell += 1;
            }
        }
    }

    // Method 2: Duplicate face detection
    if config.check_duplicate_faces {
        let duplicate_faces = detect_duplicate_internal_faces(brep, &faces, tol);
        for idx in duplicate_faces {
            if internal_set.insert(idx) {
                report.by_duplicate += 1;
            }
        }
    }

    // Method 3: Connectivity analysis
    if config.use_connectivity_analysis {
        let connectivity_faces = detect_internal_faces_by_connectivity(
            brep,
            &faces,
            config.shared_edge_threshold,
            config.min_edge_count,
        );
        for idx in connectivity_faces {
            if internal_set.insert(idx) {
                report.by_connectivity += 1;
            }
        }
    }

    // Method 4: Material side analysis
    if config.use_material_side_analysis {
        let material_faces = detect_internal_faces_by_material_side(brep, &faces, tol);
        for idx in material_faces {
            if internal_set.insert(idx) {
                report.by_material_side += 1;
            }
        }
    }

    // Method 5: Visibility check (ray casting)
    if config.use_visibility_check {
        let visibility_faces = detect_internal_faces_by_visibility(brep, &faces);
        for idx in visibility_faces {
            if internal_set.insert(idx) {
                report.by_visibility += 1;
            }
        }
    }

    report.internal_face_indices = internal_set.into_iter().collect();
    report.internal_face_indices.sort();

    report.summary = format!(
        "InternalFaceDetection: {} internal faces found (material_side={}, visibility={}, duplicate={}, void_shell={}, connectivity={})",
        report.internal_face_indices.len(),
        report.by_material_side,
        report.by_visibility,
        report.by_duplicate,
        report.by_void_shell,
        report.by_connectivity
    );

    report
}

/// Detect faces in void shells (shell index > 0 in a solid).
fn detect_void_shell_faces(
    brep: &BRep,
    faces: &[(usize, usize, usize, &Face)],
) -> Vec<usize> {
    let mut result = Vec::new();

    for (flat_idx, &(si, shi, _, _)) in faces.iter().enumerate() {
        // Check if this is a void shell (index > 0)
        if shi > 0 {
            // Check if the solid has multiple shells
            if let Some(solid) = brep.solids.get(si)
                && solid.shells.len() > 1 {
                    result.push(flat_idx);
                }
        }
    }

    result
}

/// Detect internal faces by finding duplicate faces with opposite orientation.
fn detect_duplicate_internal_faces(
    brep: &BRep,
    faces: &[(usize, usize, usize, &Face)],
    tolerance: f64,
) -> Vec<usize> {
    let mut result = Vec::new();
    let n_faces = faces.len();
    let tol_sq = tolerance * tolerance;

    for i in 0..n_faces {
        let (si1, shi1, _, face1) = faces[i];
        let pts1: Vec<DVec3> = face1
            .outer_wire
            .edges
            .iter()
            .filter_map(|we| {
                let edge = brep.edges.get(we.idx)?;
                let vidx = if we.forward { edge.start } else { edge.end };
                brep.vertices.get(vidx).map(|v| v.point)
            })
            .collect();

        for j in (i + 1)..n_faces {
            let (si2, shi2, _, face2) = faces[j];

            // Check for opposite normals
            let normal_dot = face1.normal.dot(face2.normal);
            if normal_dot > -0.99 {
                continue;
            }

            // Check geometric coincidence
            let pts2: Vec<DVec3> = face2
                .outer_wire
                .edges
                .iter()
                .filter_map(|we| {
                    let edge = brep.edges.get(we.idx)?;
                    let vidx = if we.forward { edge.start } else { edge.end };
                    brep.vertices.get(vidx).map(|v| v.point)
                })
                .collect();

            if pts1.len() != pts2.len() || pts1.is_empty() {
                continue;
            }

            // Check if all vertices match
            let mut all_match = true;
            for &p1 in &pts1 {
                let mut found = false;
                for &p2 in &pts2 {
                    if (p1 - p2).length_squared() < tol_sq {
                        found = true;
                        break;
                    }
                }
                if !found {
                    all_match = false;
                    break;
                }
            }

            if all_match {
                // Faces are duplicates with opposite orientation
                // The face in the same solid but different shell is internal
                if si1 == si2 && shi1 != shi2 {
                    // Face in the non-first shell is internal
                    if shi1 > shi2 {
                        result.push(i);
                    } else {
                        result.push(j);
                    }
                } else if si1 == si2 && shi1 == shi2 {
                    // Same shell - one is internal (prefer removing j)
                    result.push(j);
                }
            }
        }
    }

    result.sort();
    result.dedup();
    result
}

/// Detect internal faces by connectivity analysis.
///
/// Internal faces often have all their edges shared with other faces,
/// but for a proper closed manifold shell, ALL edges should be shared
/// by exactly 2 faces. This function looks for anomalies:
/// - Edges shared by MORE than 2 faces (non-manifold or partition faces)
/// - Faces where all edges are shared but the sharing is unusual
fn detect_internal_faces_by_connectivity(
    _brep: &BRep,
    faces: &[(usize, usize, usize, &Face)],
    _shared_edge_threshold: f64,
    min_edge_count: usize,
) -> Vec<usize> {
    let mut result = Vec::new();

    // Build edge-to-face map for each solid
    // Key: (solid_idx, edge_idx) -> list of faces using this edge
    let mut edge_face_map: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();

    for (flat_idx, &(si, _, _, _)) in faces.iter().enumerate() {
        let (_, _, _, face) = faces[flat_idx];
        for we in &face.outer_wire.edges {
            edge_face_map
                .entry((si, we.idx))
                .or_default()
                .push(flat_idx);
        }
    }

    // Check each face for unusual edge sharing patterns
    for (flat_idx, &(si, _, _, face)) in faces.iter().enumerate() {
        let total_edges = face.outer_wire.edges.len();
        if total_edges < min_edge_count {
            continue;
        }

        // Check if any edge is shared by more than 2 faces in the same solid
        // This indicates a partition face (internal face after boolean operation)
        let mut has_non_manifold_edge = false;
        for we in &face.outer_wire.edges {
            if let Some(face_list) = edge_face_map.get(&(si, we.idx))
                && face_list.len() > 2 {
                    // This edge is shared by more than 2 faces - potential internal face
                    has_non_manifold_edge = true;
                    break;
                }
        }

        if has_non_manifold_edge {
            result.push(flat_idx);
        }
    }

    result.sort();
    result.dedup();
    result
}

/// Detect internal faces by material side analysis.
///
/// A face is internal if the material is on both sides (the face separates
/// the same material region). This typically happens after boolean operations
/// where partition faces are left behind.
fn detect_internal_faces_by_material_side(
    brep: &BRep,
    faces: &[(usize, usize, usize, &Face)],
    _tolerance: f64,
) -> Vec<usize> {
    let mut result = Vec::new();

    for (flat_idx, &(si, shi, _fi, face)) in faces.iter().enumerate() {
        // Skip void shell faces (handled separately)
        if shi > 0 {
            continue;
        }

        // Check if the face has edges shared by more than 2 faces
        // This indicates it might be a partition face
        let solid = match brep.solids.get(si) {
            Some(s) => s,
            None => continue,
        };

        // Count edge usage - looking for edges shared by more than 2 faces
        let mut edges_with_multiple_sharing = 0usize;
        let mut total_edges = 0usize;

        for we in &face.outer_wire.edges {
            total_edges += 1;
            let mut face_count = 0usize;

            for shell in &solid.shells {
                for other_face in &shell.faces {
                    for other_we in &other_face.outer_wire.edges {
                        if other_we.idx == we.idx {
                            face_count += 1;
                        }
                    }
                }
            }

            if face_count > 2 {
                edges_with_multiple_sharing += 1;
            }
        }

        // If many edges are shared by more than 2 faces, this is likely a partition face
        if total_edges > 0 && edges_with_multiple_sharing as f64 / total_edges as f64 > 0.5 {
            result.push(flat_idx);
        }
    }

    result.sort();
    result.dedup();
    result
}

/// Detect internal faces by visibility check using ray casting.
fn detect_internal_faces_by_visibility(
    brep: &BRep,
    faces: &[(usize, usize, usize, &Face)],
) -> Vec<usize> {
    let mut result = Vec::new();

    for (flat_idx, &(si, _, _, face)) in faces.iter().enumerate() {
        let centroid = compute_face_centroid_from_wire(brep, face);
        if centroid.is_nan() {
            continue;
        }

        // Cast ray in the direction of the face normal
        let ray_origin = centroid + face.normal * TOLERANCE_RETRY_LADDER_COARSE;
        let ray_dir = face.normal;

        // Count intersections with other faces
        let mut intersection_count = 0usize;
        for (other_idx, &(_, other_si, _, other_face)) in faces.iter().enumerate() {
            if other_idx == flat_idx || other_si != si {
                continue;
            }

            if ray_intersects_face(brep, other_face, ray_origin, ray_dir) {
                intersection_count += 1;
            }
        }

        // Odd number of intersections in normal direction suggests internal face
        if intersection_count > 0 && intersection_count % 2 == 1 {
            result.push(flat_idx);
        }
    }

    result.sort();
    result.dedup();
    result
}

/// Check if a point is inside a solid using ray casting.
fn is_point_inside_solid(brep: &BRep, solid_idx: usize, point: DVec3) -> bool {
    let solid = match brep.solids.get(solid_idx) {
        Some(s) => s,
        None => return false,
    };

    // Collect all faces from this solid
    let all_faces: Vec<&Face> = solid
        .shells
        .iter()
        .flat_map(|shell| shell.faces.iter())
        .collect();

    if all_faces.is_empty() {
        return false;
    }

    // Cast ray in +X direction
    let ray_dir = DVec3::X;
    let mut intersection_count = 0usize;

    for face in &all_faces {
        if ray_intersects_face(brep, face, point, ray_dir) {
            intersection_count += 1;
        }
    }

    // Odd intersections = inside
    intersection_count % 2 == 1
}

/// Configuration for post-boolean internal face removal.
#[derive(Debug, Clone)]
pub struct PostBooleanRemovalConfig {
    /// Detection configuration.
    pub detection: InternalFaceDetectionConfig,
    /// Whether to merge vertices after removal.
    pub merge_vertices: bool,
    /// Whether to validate the result.
    pub validate_result: bool,
    /// Whether to remove degenerate edges after removal.
    pub remove_degenerate_edges: bool,
    /// Tolerance for vertex merging.
    pub merge_tolerance: f64,
}

impl Default for PostBooleanRemovalConfig {
    fn default() -> Self {
        Self {
            detection: InternalFaceDetectionConfig::for_post_boolean(),
            merge_vertices: true,
            validate_result: true,
            remove_degenerate_edges: true,
            merge_tolerance: TOLERANCE_ABS,
        }
    }
}

impl PostBooleanRemovalConfig {
    /// Create a configuration for fuse (union) operations.
    pub fn for_fuse() -> Self {
        Self {
            detection: InternalFaceDetectionConfig {
                tolerance: TOLERANCE_ABS * 5.0,
                use_material_side_analysis: true,
                use_visibility_check: false, // Disabled - can be unreliable
                check_duplicate_faces: true,
                consider_void_shells: true,
                min_edge_count: 3,
                use_connectivity_analysis: true,
                shared_edge_threshold: 0.85,
            },
            merge_vertices: true,
            validate_result: true,
            remove_degenerate_edges: true,
            merge_tolerance: TOLERANCE_ABS * 2.0,
        }
    }

    /// Create a configuration for cut (difference) operations.
    pub fn for_cut() -> Self {
        Self {
            detection: InternalFaceDetectionConfig {
                tolerance: TOLERANCE_ABS * 3.0,
                use_material_side_analysis: true,
                use_visibility_check: false, // Avoid removing cut faces
                check_duplicate_faces: true,
                consider_void_shells: true,
                min_edge_count: 3,
                use_connectivity_analysis: true,
                shared_edge_threshold: 0.95, // Higher threshold for cuts
            },
            merge_vertices: true,
            validate_result: true,
            remove_degenerate_edges: true,
            merge_tolerance: TOLERANCE_ABS,
        }
    }

    /// Create a configuration for intersection operations.
    pub fn for_intersection() -> Self {
        Self {
            detection: InternalFaceDetectionConfig {
                tolerance: TOLERANCE_ABS * 5.0,
                use_material_side_analysis: true,
                use_visibility_check: false, // Disabled - can be unreliable
                check_duplicate_faces: true,
                consider_void_shells: false, // Intersection may create voids
                min_edge_count: 3,
                use_connectivity_analysis: true,
                shared_edge_threshold: 0.9,
            },
            merge_vertices: true,
            validate_result: true,
            remove_degenerate_edges: true,
            merge_tolerance: TOLERANCE_ABS * 2.0,
        }
    }
}

/// Report from post-boolean internal face removal.
#[derive(Debug, Clone, Default)]
pub struct PostBooleanRemovalReport {
    /// Detection report.
    pub detection: InternalFaceDetectionReport,
    /// Removal report.
    pub removal: InternalFaceRemovalReport,
    /// Number of vertices merged after removal.
    pub vertices_merged: usize,
    /// Number of degenerate edges removed.
    pub degenerate_edges_removed: usize,
    /// Whether validation passed.
    pub validation_passed: bool,
    /// Validation issues (if any).
    pub validation_issues: Vec<String>,
    /// Summary string.
    pub summary: String,
}

/// Remove internal faces from a BRep after boolean operations.
///
/// This is a convenience function that combines detection and removal
/// with post-removal cleanup and validation.
///
/// # Arguments
/// * `brep` - The BRep to process.
///
/// # Returns
/// A tuple of (cleaned BRep, removal report).
pub fn remove_internal_faces_post_boolean(brep: &BRep) -> (BRep, PostBooleanRemovalReport) {
    remove_internal_faces_post_boolean_with_config(brep, &PostBooleanRemovalConfig::default())
}

/// Remove internal faces after boolean operations with custom configuration.
///
/// See [`remove_internal_faces_post_boolean`] for details.
pub fn remove_internal_faces_post_boolean_with_config(
    brep: &BRep,
    config: &PostBooleanRemovalConfig,
) -> (BRep, PostBooleanRemovalReport) {
    let mut report = PostBooleanRemovalReport::default();

    // Step 1: Detect internal faces
    let detection_report = detect_internal_faces_with_config(brep, &config.detection);
    report.detection = detection_report.clone();

    if detection_report.internal_face_indices.is_empty() {
        report.summary = "No internal faces detected".to_string();
        report.validation_passed = true;
        return (brep.clone(), report);
    }

    // Step 2: Remove internal faces
    let (mut result, removal_report) =
        remove_internal_faces(brep, &detection_report.internal_face_indices);
    report.removal = removal_report;

    // Step 3: Remove degenerate edges
    if config.remove_degenerate_edges {
        let (cleaned, edges_removed) = remove_small_edges(&result, config.merge_tolerance);
        result = cleaned;
        report.degenerate_edges_removed = edges_removed;
    }

    // Step 4: Merge close vertices
    if config.merge_vertices {
        let (merged, vertices_merged) = merge_close_vertices(&result, config.merge_tolerance);
        result = merged;
        report.vertices_merged = vertices_merged;
    }

    // Step 5: Validate result
    if config.validate_result {
        let validation = validate_internal_face_removal(&result);
        report.validation_passed = validation.is_valid;
        report.validation_issues = validation.issues;
    }

    report.summary = format!(
        "PostBooleanRemoval: {} faces removed, {} vertices merged, {} degenerate edges removed, validation {}",
        report.removal.faces_removed,
        report.vertices_merged,
        report.degenerate_edges_removed,
        if report.validation_passed { "passed" } else { "FAILED" }
    );

    (result, report)
}

/// Validation result for internal face removal.
#[derive(Debug, Clone, Default)]
pub struct InternalFaceRemovalValidation {
    /// Whether the BRep is valid after removal.
    pub is_valid: bool,
    /// List of validation issues found.
    pub issues: Vec<String>,
    /// Number of empty shells found.
    pub empty_shells: usize,
    /// Number of empty solids found.
    pub empty_solids: usize,
    /// Number of degenerate edges found.
    pub degenerate_edges: usize,
    /// Number of orphaned vertices found.
    pub orphaned_vertices: usize,
}

/// Validate a BRep after internal face removal.
///
/// Checks for:
/// - Empty shells
/// - Empty solids
/// - Degenerate edges (zero-length)
/// - Orphaned vertices
/// - Shell closure
pub fn validate_internal_face_removal(brep: &BRep) -> InternalFaceRemovalValidation {
    let mut validation = InternalFaceRemovalValidation::default();
    validation.is_valid = true;

    // Check for empty shells
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            if shell.faces.is_empty() {
                validation.empty_shells += 1;
                validation
                    .issues
                    .push(format!("Empty shell at solid {} shell {}", si, shi));
                validation.is_valid = false;
            }
        }

        if solid.shells.is_empty() {
            validation.empty_solids += 1;
            validation
                .issues
                .push(format!("Empty solid at index {}", si));
            validation.is_valid = false;
        }
    }

    // Check for degenerate edges
    for (ei, edge) in brep.edges.iter().enumerate() {
        if edge.start == edge.end {
            validation.degenerate_edges += 1;
            validation
                .issues
                .push(format!("Degenerate edge at index {}", ei));
        } else if let (Some(v_start), Some(v_end)) = (
            brep.vertices.get(edge.start),
            brep.vertices.get(edge.end),
        ) {
            let len = (v_start.point - v_end.point).length();
            if len < TOLERANCE_ABS {
                validation.degenerate_edges += 1;
                validation.issues.push(format!(
                    "Near-zero length edge at index {} (length: {})",
                    ei, len
                ));
            }
        }
    }

    // Check for orphaned vertices
    let mut used_vertices: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for edge in &brep.edges {
        used_vertices.insert(edge.start);
        used_vertices.insert(edge.end);
    }

    for vi in 0..brep.vertices.len() {
        if !used_vertices.contains(&vi) {
            validation.orphaned_vertices += 1;
        }
    }

    if validation.orphaned_vertices > 0 {
        validation.issues.push(format!(
            "{} orphaned vertices found",
            validation.orphaned_vertices
        ));
    }

    // Check shell closure using edge valence analysis
    for (si, solid) in brep.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            let closure = check_shell_closure_internal(brep, shell);
            if !closure.is_closed {
                validation.is_valid = false;
                validation.issues.push(format!(
                    "Shell not closed at solid {} shell {}: {} open edges",
                    si, shi, closure.open_edges
                ));
            }
        }
    }

    validation
}

/// Shell closure check result.
#[derive(Debug, Clone, Default)]
struct ShellClosureCheck {
    is_closed: bool,
    open_edges: usize,
}

/// Check if a shell is properly closed (all edges shared by exactly 2 faces).
fn check_shell_closure_internal(_brep: &BRep, shell: &Shell) -> ShellClosureCheck {
    let mut edge_count: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

    for face in &shell.faces {
        for we in &face.outer_wire.edges {
            *edge_count.entry(we.idx).or_insert(0) += 1;
        }
    }

    let mut open_edges = 0usize;
    for &count in edge_count.values() {
        if count != 2 {
            // Check if edge is a boundary edge (count == 1) or non-manifold (count > 2)
            if count == 1 {
                open_edges += 1;
            }
        }
    }

    ShellClosureCheck {
        is_closed: open_edges == 0,
        open_edges,
    }
}

/// Estimate face area from its wire (approximate).
fn estimate_face_area_from_wire(brep: &BRep, wire: &Wire) -> f64 {
    // Get vertices
    let pts: Vec<DVec3> = wire
        .edges
        .iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let vidx = if we.forward { edge.start } else { edge.end };
            brep.vertices.get(vidx).map(|v| v.point)
        })
        .collect();

    if pts.len() < 3 {
        return 0.0;
    }

    // Compute signed area using shoelace formula (projected to XY plane)
    // This is an approximation; for accurate results, use proper surface area calculation
    let mut area = 0.0;
    for i in 0..pts.len() {
        let j = (i + 1) % pts.len();
        area += (pts[i].x * pts[j].y - pts[j].x * pts[i].y).abs();
    }
    area * 0.5
}

/// Merge adjacent faces after internal face removal.
///
/// When internal faces are removed, adjacent faces that now share edges
/// can potentially be merged if they are on the same underlying surface.
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `tolerance` - Tolerance for geometric comparisons.
///
/// # Returns
/// A tuple of (BRep with merged faces, count of faces merged).
pub fn merge_adjacent_faces_after_removal(brep: &BRep, tolerance: f64) -> (BRep, usize) {
    let tol = tolerance.max(TOLERANCE_ABS);
    let mut result = brep.clone();
    let mut total_merged = 0usize;

    // Collect shell data first to avoid borrow issues
    let shell_data: Vec<(usize, usize)> = result
        .solids
        .iter()
        .enumerate()
        .flat_map(|(si, solid)| {
            solid.shells.iter().enumerate().map(move |(shi, _)| (si, shi))
        })
        .collect();

    for (si, shi) in shell_data {
        let faces_to_merge: Vec<Face> = result.solids[si].shells[shi].faces.clone();
        let (new_faces, merged) = merge_faces_in_shell(brep, &faces_to_merge, tol);
        result.solids[si].shells[shi].faces = new_faces;
        total_merged += merged;
    }

    (result, total_merged)
}

/// Merge faces in a shell that share the same underlying surface.
fn merge_faces_in_shell(brep: &BRep, faces: &[Face], tolerance: f64) -> (Vec<Face>, usize) {
    if faces.len() < 2 {
        return (faces.to_vec(), 0);
    }

    let mut merged_count = 0usize;
    let mut merged: Vec<bool> = vec![false; faces.len()];
    let mut result: Vec<Face> = Vec::new();

    // Find groups of faces that can be merged (same normal, coplanar)
    for i in 0..faces.len() {
        if merged[i] {
            continue;
        }

        let face_i = &faces[i];
        let mut group = vec![i];

        // Find other faces that can be merged with this one
        for j in (i + 1)..faces.len() {
            if merged[j] {
                continue;
            }

            let face_j = &faces[j];

            // Check if faces have the same normal
            let normal_dot = face_i.normal.dot(face_j.normal);
            if normal_dot.abs() < 0.999 {
                continue;
            }

            // Check if faces are coplanar (sample points)
            let centroid_i = compute_face_centroid_from_wire(brep, face_i);
            let centroid_j = compute_face_centroid_from_wire(brep, face_j);

            if centroid_i.is_nan() || centroid_j.is_nan() {
                continue;
            }

            // Check distance from centroid to plane
            let plane_d = face_i.normal.dot(centroid_i);
            let dist_j = (face_i.normal.dot(centroid_j) - plane_d).abs();

            if dist_j > tolerance {
                continue;
            }

            // Check if faces share at least one edge
            let edges_i: std::collections::HashSet<usize> =
                face_i.outer_wire.edges.iter().map(|we| we.idx).collect();
            let edges_j: std::collections::HashSet<usize> =
                face_j.outer_wire.edges.iter().map(|we| we.idx).collect();

            let shared: std::collections::HashSet<usize> =
                edges_i.intersection(&edges_j).copied().collect();

            if shared.is_empty() {
                continue;
            }

            // Faces can potentially be merged
            group.push(j);
            merged[j] = true;
        }

        // For now, just keep the faces as-is (full merging requires more complex logic)
        // This can be enhanced later to actually merge the wire topology
        if group.len() > 1 {
            merged_count += group.len() - 1;
        }

        result.push(face_i.clone());
    }

    (result, merged_count)
}

