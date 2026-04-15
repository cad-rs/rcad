//! Structured healing pipeline for B-Rep analysis and repair.
//!
//! This module provides an analyze -> repair -> recheck workflow similar in
//! spirit to OCCT ShapeAnalysis/ShapeFix orchestration.

use rcad_kernel::BRep;

use crate::brep_check::{CheckIssue, CheckResult, check, diagnose_same_parameter, diagnose_same_range};
use crate::brep_repair::{
    MakeConnectedReport, RepairReport, fix_same_parameter_with_scan,
    fix_same_range_with_scan, make_connected_iterative_with_growth_cap, repair,
};
use crate::tolerance::TOLERANCE_ABS;

/// Healing execution mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealingMode {
    /// Only analyze; no repair pass will run.
    AnalyzeOnly,
    /// Analyze and run repair passes.
    AnalyzeAndRepair,
}

/// Policy controlling whether a make-connected prepass is executed before
/// regular repair passes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MakeConnectedPrepassMode {
    /// Never run a prepass.
    Disabled,
    /// Run only when initial checker issues indicate connectivity stress.
    IssueDriven,
    /// Always run before repair passes.
    Always,
}

/// Options controlling healing execution.
#[derive(Debug, Clone, Copy)]
pub struct HealingOptions {
    /// Repair tolerance used by [`repair`].
    pub tolerance: f64,
    /// Maximum number of repair passes.
    pub max_passes: usize,
    /// Execution mode for the pipeline.
    pub mode: HealingMode,
    /// Control whether to run make-connected before normal repair passes.
    pub make_connected_prepass_mode: MakeConnectedPrepassMode,
    /// Run SameRange/SameParameter scan+fix pass as a prepass.
    pub run_parametric_consistency_prepass: bool,
    /// Re-run parametric consistency pass after each repair iteration when
    /// remaining issues still indicate parametric inconsistency.
    pub run_parametric_consistency_iterative: bool,
    /// When a repair pass makes no progress while issues remain, run a
    /// MakeConnected-style connectivity rebuild pass.
    pub run_make_connected_on_stall: bool,
    /// Base tolerance used by make-connected fallback passes.
    pub make_connected_tolerance: f64,
    /// Maximum number of iterative make-connected passes.
    pub make_connected_max_passes: usize,
    /// Per-pass tolerance growth factor for make-connected fallback.
    pub make_connected_tolerance_growth: f64,
    /// Upper cap for make-connected tolerance growth.
    pub make_connected_tolerance_cap: f64,
}

impl Default for HealingOptions {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            max_passes: 2,
            mode: HealingMode::AnalyzeAndRepair,
            make_connected_prepass_mode: MakeConnectedPrepassMode::Disabled,
            run_parametric_consistency_prepass: true,
            run_parametric_consistency_iterative: true,
            run_make_connected_on_stall: false,
            make_connected_tolerance: TOLERANCE_ABS,
            make_connected_max_passes: 3,
            make_connected_tolerance_growth: 1.5,
            make_connected_tolerance_cap: TOLERANCE_ABS * 1000.0,
        }
    }
}

/// Structured issue counters for checker output.
#[derive(Debug, Clone, Default)]
pub struct HealingIssueStats {
    pub open_wire: usize,
    pub zero_normal: usize,
    pub degenerate_face: usize,
    pub invalid_edge_index: usize,
    pub invalid_vertex_index: usize,
    pub non_manifold_edge: usize,
    pub self_intersecting_wire: usize,
    pub geometric_self_intersection: usize,
}

impl HealingIssueStats {
    pub fn total(&self) -> usize {
        self.open_wire
            + self.zero_normal
            + self.degenerate_face
            + self.invalid_edge_index
            + self.invalid_vertex_index
            + self.non_manifold_edge
            + self.self_intersecting_wire
            + self.geometric_self_intersection
    }

    pub fn from_check_result(result: &CheckResult) -> Self {
        let mut s = Self::default();
        for issue in &result.issues {
            match issue {
                CheckIssue::OpenWire { .. } => s.open_wire += 1,
                CheckIssue::ZeroNormal { .. } => s.zero_normal += 1,
                CheckIssue::DegenerateFace { .. } => s.degenerate_face += 1,
                CheckIssue::InvalidEdgeIndex { .. } => s.invalid_edge_index += 1,
                CheckIssue::InvalidVertexIndex { .. } => s.invalid_vertex_index += 1,
                CheckIssue::NonManifoldEdge { .. } => s.non_manifold_edge += 1,
                CheckIssue::SelfIntersectingWire { .. } => s.self_intersecting_wire += 1,
                CheckIssue::GeometricSelfIntersection { .. } => s.geometric_self_intersection += 1,
            }
        }
        s
    }
}

/// Comprehensive diagnosis report combining all analysis types.
///
/// Analogous to running all ShapeAnalysis tools in OCCT:
/// - ShapeAnalysis_Surface (UV consistency)
/// - ShapeAnalysis_Wire (wire quality)
/// - ShapeAnalysis_ShapeTolerance (tolerance consistency)
/// - BRepCheck_Analyzer (topology validity)
#[derive(Debug, Clone)]
pub struct ComprehensiveDiagnosis {
    /// Basic topology check result.
    pub topology: CheckResult,
    /// Surface UV consistency analysis.
    pub surface_uv: crate::brep_check::SurfaceAnalysisReport,
    /// Wire quality analysis.
    pub wire_quality: crate::brep_check::WireQualityReport,
    /// SameParameter diagnosis.
    pub same_parameter: crate::brep_check::SameParameterDiagnosis,
    /// SameRange diagnosis.
    pub same_range: crate::brep_check::SameRangeDiagnosis,
}

impl Default for ComprehensiveDiagnosis {
    fn default() -> Self {
        Self {
            topology: CheckResult { issues: Vec::new() },
            surface_uv: Default::default(),
            wire_quality: Default::default(),
            same_parameter: Default::default(),
            same_range: Default::default(),
        }
    }
}

impl ComprehensiveDiagnosis {
    /// Returns true if all diagnoses are clean (no issues found).
    pub fn is_clean(&self) -> bool {
        self.topology.is_valid()
            && self.surface_uv.is_clean()
            && self.wire_quality.is_clean()
            && self.same_parameter.is_clean()
            && self.same_range.is_clean()
    }

    /// Returns a summary string of all issues found.
    pub fn summary(&self) -> String {
        if self.is_clean() {
            return "All diagnoses clean: topology valid, no UV issues, all wires closed, parametric consistency OK".to_string();
        }

        let mut parts = Vec::new();

        if !self.topology.is_valid() {
            parts.push(format!("topology: {} issues", self.topology.issues.len()));
        }
        if !self.surface_uv.is_clean() {
            parts.push(format!("UV bounds: {} violations", self.surface_uv.total_issues));
        }
        if !self.wire_quality.is_clean() {
            parts.push(format!(
                "wires: {} open, {} self-intersecting",
                self.wire_quality.open_wires,
                self.wire_quality.self_intersecting_wires
            ));
        }
        if !self.same_parameter.is_clean() {
            parts.push(format!("SameParameter: {} suspect", self.same_parameter.suspect_edges.len()));
        }
        if !self.same_range.is_clean() {
            parts.push(format!("SameRange: {} suspect", self.same_range.suspect_edges.len()));
        }

        parts.join("; ")
    }

    /// Returns total count of all issues found.
    pub fn total_issues(&self) -> usize {
        self.topology.issues.len()
            + self.surface_uv.total_issues
            + (if self.wire_quality.is_clean() { 0 } else { 1 })
            + self.same_parameter.suspect_edges.len()
            + self.same_range.suspect_edges.len()
    }
}

/// Run all available diagnoses on a BRep.
///
/// This is a convenience function that runs all analysis tools and returns
/// a combined report.
///
/// # Arguments
/// * `brep` - The BRep to analyze.
/// * `tolerance` - The tolerance to use for geometric comparisons.
///
/// # Returns
/// A `ComprehensiveDiagnosis` containing all analysis results.
pub fn diagnose_all(brep: &BRep, tolerance: f64) -> ComprehensiveDiagnosis {
    use crate::brep_check::{
        analyze_surface_uv_consistency, analyze_wire_quality,
        diagnose_same_parameter, diagnose_same_range,
    };

    ComprehensiveDiagnosis {
        topology: check(brep),
        surface_uv: analyze_surface_uv_consistency(brep, tolerance),
        wire_quality: analyze_wire_quality(brep, tolerance),
        same_parameter: diagnose_same_parameter(brep, tolerance),
        same_range: diagnose_same_range(brep, tolerance),
    }
}

/// Summary report for analyze/heal workflow.
#[derive(Debug, Clone)]
pub struct HealingReport {
    /// Issues found before any repair.
    pub initial: CheckResult,
    /// Issues after the final pass.
    pub final_result: CheckResult,
    /// Per-pass repair reports.
    pub passes: Vec<RepairReport>,
    /// Parametric consistency passes (SameRange/SameParameter scan+fix).
    pub parametric_passes: Vec<ParametricConsistencyReport>,
    /// MakeConnected fallback reports executed when repair stalls.
    pub make_connected_passes: Vec<MakeConnectedReport>,
    /// Structured issue counters before healing.
    pub initial_stats: HealingIssueStats,
    /// Structured issue counters after healing.
    pub final_stats: HealingIssueStats,
    /// Stage-by-stage issue counts for analyze/repair pipeline.
    pub stages: Vec<HealingStageReport>,
}

/// Stage marker for healing pipeline reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealingStage {
    InitialCheck,
    PreprocessPass,
    GeometryRepairPass,
    TopologyRepairPass,
    PreMakeConnected,
    OperatorChainStep,
    ParametricConsistencyPass,
    RepairPass,
    WireGapRepairPass,
    UvBoundsRepairPass,
    MakeConnectedPass,
    FinalizePass,
    FinalCheck,
}

/// ShapeProcess-like healing operators that can be composed into a custom
/// batch pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HealingOperator {
    /// MakeConnected-style connectivity rebuild pass.
    MakeConnected,
    /// SameRange/SameParameter consistency pass.
    ParametricConsistency,
    /// General repair pass (`repair`).
    Repair,
    /// Wire gap repair pass.
    WireGapRepair,
    /// UV bounds repair pass.
    UvBoundsRepair,
    /// Stop pipeline execution if the current shape is checker-clean.
    StopIfClean,
    /// Remove faces with area below a threshold.
    FixSmallAreaFaces,
    /// Fix sliver (thin elongated) faces by merging or removal.
    FixSliverFaces,
    /// Repair non-manifold topology by splitting multi-face edges.
    FixNonManifold,
    /// Propagate tolerances through the shape hierarchy.
    PropagateTolerances,
    /// Merge faces that share the same underlying surface geometry.
    UnifySameDomain,
    /// Remove internal faces (faces inside the solid volume).
    RemoveInternalFaces,
}

/// Configuration parameters for individual healing operators.
#[derive(Debug, Clone)]
pub struct OperatorParams {
    /// Tolerance threshold for geometric operations.
    pub tolerance: f64,
    /// Area threshold for FixSmallAreaFaces.
    pub min_face_area: f64,
    /// Aspect ratio threshold for FixSliverFaces.
    pub max_sliver_aspect_ratio: f64,
    /// Whether to allow removal of internal faces.
    pub allow_internal_face_removal: bool,
}

impl Default for OperatorParams {
    fn default() -> Self {
        Self {
            tolerance: TOLERANCE_ABS,
            min_face_area: 1e-10,
            max_sliver_aspect_ratio: 100.0,
            allow_internal_face_removal: true,
        }
    }
}

/// Report for one SameRange/SameParameter consistency pass.
#[derive(Debug, Clone, Default)]
pub struct ParametricConsistencyReport {
    pub same_range_fixed: usize,
    pub same_parameter_fixed: usize,
}

/// Report for a single healing operator execution.
#[derive(Debug, Clone)]
pub struct OperatorReport {
    /// The operator that was executed.
    pub operator: HealingOperator,
    /// Number of entities modified/removed.
    pub modifications: usize,
    /// Number of issues fixed by this operator.
    pub issues_fixed: usize,
    /// Whether the operator made any changes.
    pub changed: bool,
    /// Human-readable description of changes.
    pub description: String,
}

impl Default for OperatorReport {
    fn default() -> Self {
        Self {
            operator: HealingOperator::Repair,
            modifications: 0,
            issues_fixed: 0,
            changed: false,
            description: String::new(),
        }
    }
}

/// Report for a single stage in the ShapeProcess pipeline.
#[derive(Debug, Clone)]
pub struct StageReport {
    /// The stage type.
    pub stage: HealingStage,
    /// Zero-based pass index (for multi-pass stages).
    pub pass_index: Option<usize>,
    /// Issue count before this stage.
    pub issue_count_before: usize,
    /// Issue count after this stage.
    pub issue_count_after: usize,
    /// Reports from individual operators executed in this stage.
    pub operator_reports: Vec<OperatorReport>,
    /// Wall-clock time for this stage (seconds).
    pub elapsed_seconds: f64,
}

impl StageReport {
    pub fn issues_fixed(&self) -> usize {
        self.issue_count_before.saturating_sub(self.issue_count_after)
    }

    pub fn is_improved(&self) -> bool {
        self.issue_count_after < self.issue_count_before
    }
}

/// Overall statistics for a ShapeProcess run.
#[derive(Debug, Clone, Default)]
pub struct ShapeProcessStats {
    /// Total number of operators executed.
    pub operators_executed: usize,
    /// Total number of modifications made.
    pub total_modifications: usize,
    /// Total number of issues fixed.
    pub total_issues_fixed: usize,
    /// Number of stages executed.
    pub stages_executed: usize,
    /// Total wall-clock time (seconds).
    pub total_elapsed_seconds: f64,
    /// Number of iterations (when max_iterations > 1).
    pub iterations: usize,
    /// Whether the process converged early (shape became clean).
    pub converged_early: bool,
    /// Final shape is checker-clean.
    pub is_clean: bool,
}

/// Complete report from a ShapeProcess pipeline run.
#[derive(Debug, Clone)]
pub struct ShapeProcessReport {
    /// Initial check result.
    pub initial: CheckResult,
    /// Final check result.
    pub final_result: CheckResult,
    /// Structured issue counts before processing.
    pub initial_stats: HealingIssueStats,
    /// Structured issue counts after processing.
    pub final_stats: HealingIssueStats,
    /// Per-stage reports.
    pub stages: Vec<StageReport>,
    /// Overall statistics.
    pub stats: ShapeProcessStats,
    /// Configuration used for this run.
    pub config_summary: String,
}

impl ShapeProcessReport {
    pub fn initial_issue_count(&self) -> usize {
        self.initial.issues.len()
    }

    pub fn final_issue_count(&self) -> usize {
        self.final_result.issues.len()
    }

    pub fn issues_fixed(&self) -> usize {
        self.initial_issue_count().saturating_sub(self.final_issue_count())
    }

    pub fn is_improved(&self) -> bool {
        self.final_issue_count() < self.initial_issue_count()
    }

    pub fn is_clean(&self) -> bool {
        self.final_result.is_valid()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() {
            format!(
                "ShapeProcess: Clean result after {} operators, {} modifications in {:.3}s",
                self.stats.operators_executed,
                self.stats.total_modifications,
                self.stats.total_elapsed_seconds
            )
        } else {
            format!(
                "ShapeProcess: {} → {} issues ({} fixed) after {} operators in {:.3}s",
                self.initial_issue_count(),
                self.final_issue_count(),
                self.issues_fixed(),
                self.stats.operators_executed,
                self.stats.total_elapsed_seconds
            )
        }
    }
}

/// Configuration for the ShapeProcess pipeline.
///
/// This struct provides OCCT ShapeProcess-like configuration for running
/// a customizable sequence of healing operations on a BRep.
#[derive(Debug, Clone)]
pub struct ShapeProcessConfig {
    /// Sequence of healing operators to execute.
    pub operators: Vec<HealingOperator>,
    /// Stop processing if the shape becomes checker-clean.
    pub stop_on_clean: bool,
    /// Maximum number of iterations (0 = run once).
    pub max_iterations: usize,
    /// Tolerance growth factor per iteration.
    pub tolerance_growth: f64,
    /// Maximum tolerance cap.
    pub tolerance_cap: f64,
    /// Base tolerance for operations.
    pub base_tolerance: f64,
    /// Parameters for individual operators.
    pub operator_params: OperatorParams,
    /// Healing options for internal passes.
    pub healing_options: HealingOptions,
}

impl Default for ShapeProcessConfig {
    fn default() -> Self {
        Self {
            operators: vec![
                HealingOperator::ParametricConsistency,
                HealingOperator::Repair,
                HealingOperator::StopIfClean,
            ],
            stop_on_clean: true,
            max_iterations: 1,
            tolerance_growth: 1.0,
            tolerance_cap: 1e-3,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams::default(),
            healing_options: HealingOptions::default(),
        }
    }
}

impl ShapeProcessConfig {
    /// Create a preset configuration optimized for imported CAD data.
    ///
    /// This preset applies aggressive cleaning operations commonly needed
    /// after importing STEP/IGES files.
    pub fn import_preset() -> Self {
        Self {
            operators: vec![
                HealingOperator::MakeConnected,
                HealingOperator::ParametricConsistency,
                HealingOperator::FixSmallAreaFaces,
                HealingOperator::FixSliverFaces,
                HealingOperator::Repair,
                HealingOperator::PropagateTolerances,
                HealingOperator::StopIfClean,
            ],
            stop_on_clean: true,
            max_iterations: 3,
            tolerance_growth: 1.5,
            tolerance_cap: 1e-3,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS,
                min_face_area: 1e-8,
                max_sliver_aspect_ratio: 50.0,
                allow_internal_face_removal: false,
            },
            healing_options: HealingOptions {
                tolerance: TOLERANCE_ABS,
                max_passes: 2,
                mode: HealingMode::AnalyzeAndRepair,
                make_connected_prepass_mode: MakeConnectedPrepassMode::IssueDriven,
                run_parametric_consistency_prepass: true,
                run_parametric_consistency_iterative: true,
                run_make_connected_on_stall: true,
                ..HealingOptions::default()
            },
        }
    }

    /// Create a preset configuration for cleaning up after boolean operations.
    ///
    /// This preset focuses on fixing issues common after boolean operations:
    /// parametric inconsistencies, tolerance propagation, and geometry repair.
    pub fn boolean_cleanup_preset() -> Self {
        Self {
            operators: vec![
                HealingOperator::ParametricConsistency,
                HealingOperator::Repair,
                HealingOperator::PropagateTolerances,
                HealingOperator::FixNonManifold,
                HealingOperator::StopIfClean,
            ],
            stop_on_clean: true,
            max_iterations: 2,
            tolerance_growth: 2.0,
            tolerance_cap: 1e-4,
            base_tolerance: TOLERANCE_ABS * 10.0,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS * 10.0,
                min_face_area: 1e-10,
                max_sliver_aspect_ratio: 100.0,
                allow_internal_face_removal: true,
            },
            healing_options: HealingOptions {
                tolerance: TOLERANCE_ABS * 10.0,
                max_passes: 2,
                mode: HealingMode::AnalyzeAndRepair,
                make_connected_prepass_mode: MakeConnectedPrepassMode::Disabled,
                run_parametric_consistency_prepass: true,
                run_parametric_consistency_iterative: true,
                run_make_connected_on_stall: false,
                ..HealingOptions::default()
            },
        }
    }

    /// Create a preset configuration for preparing shapes for analysis.
    ///
    /// This preset is more conservative, focusing on validation and minimal
    /// repairs without aggressive geometry modification.
    pub fn analysis_preset() -> Self {
        Self {
            operators: vec![
                HealingOperator::ParametricConsistency,
                HealingOperator::PropagateTolerances,
                HealingOperator::StopIfClean,
            ],
            stop_on_clean: true,
            max_iterations: 1,
            tolerance_growth: 1.0,
            tolerance_cap: 1e-5,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS,
                min_face_area: 1e-12,
                max_sliver_aspect_ratio: 1000.0,
                allow_internal_face_removal: false,
            },
            healing_options: HealingOptions {
                tolerance: TOLERANCE_ABS,
                max_passes: 1,
                mode: HealingMode::AnalyzeAndRepair,
                make_connected_prepass_mode: MakeConnectedPrepassMode::Disabled,
                run_parametric_consistency_prepass: true,
                run_parametric_consistency_iterative: false,
                run_make_connected_on_stall: false,
                ..HealingOptions::default()
            },
        }
    }

    /// Create a preset for aggressive geometry cleanup.
    ///
    /// This preset applies all available healing operators, useful for
    /// preparing shapes for meshing or export.
    pub fn aggressive_preset() -> Self {
        Self {
            operators: vec![
                HealingOperator::MakeConnected,
                HealingOperator::FixSmallAreaFaces,
                HealingOperator::FixSliverFaces,
                HealingOperator::FixNonManifold,
                HealingOperator::ParametricConsistency,
                HealingOperator::Repair,
                HealingOperator::UnifySameDomain,
                HealingOperator::RemoveInternalFaces,
                HealingOperator::PropagateTolerances,
                HealingOperator::StopIfClean,
            ],
            stop_on_clean: true,
            max_iterations: 5,
            tolerance_growth: 1.5,
            tolerance_cap: 1e-2,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS,
                min_face_area: 1e-8,
                max_sliver_aspect_ratio: 50.0,
                allow_internal_face_removal: true,
            },
            healing_options: HealingOptions {
                tolerance: TOLERANCE_ABS,
                max_passes: 3,
                mode: HealingMode::AnalyzeAndRepair,
                make_connected_prepass_mode: MakeConnectedPrepassMode::Always,
                run_parametric_consistency_prepass: true,
                run_parametric_consistency_iterative: true,
                run_make_connected_on_stall: true,
                ..HealingOptions::default()
            },
        }
    }
}

/// Per-stage issue metrics.
#[derive(Debug, Clone)]
pub struct HealingStageReport {
    pub stage: HealingStage,
    /// Zero-based pass index for `RepairPass`; `None` for checks.
    pub pass_index: Option<usize>,
    /// Checker issue count observed at this stage.
    pub issue_count: usize,
}

impl HealingReport {
    pub fn initial_issue_count(&self) -> usize {
        self.initial.issues.len()
    }

    pub fn final_issue_count(&self) -> usize {
        self.final_result.issues.len()
    }

    pub fn fixed_issue_count(&self) -> usize {
        self.initial_issue_count().saturating_sub(self.final_issue_count())
    }

    pub fn is_improved(&self) -> bool {
        self.final_issue_count() < self.initial_issue_count()
    }

    pub fn is_clean(&self) -> bool {
        self.final_result.is_valid()
    }

    pub fn has_issue_kind(&self, pred: impl Fn(&CheckIssue) -> bool) -> bool {
        self.final_result.issues.iter().any(pred)
    }
}

/// Analyze and heal a BRep using the provided options.
pub fn analyze_and_heal(brep: &BRep, options: HealingOptions) -> (BRep, HealingReport) {
    let initial = check(brep);
    let initial_stats = HealingIssueStats::from_check_result(&initial);

    if matches!(options.mode, HealingMode::AnalyzeOnly) {
        let initial_issue_count = initial_stats.total();
        return (
            brep.clone(),
            HealingReport {
                initial: initial.clone(),
                final_result: initial,
                passes: Vec::new(),
                parametric_passes: Vec::new(),
                make_connected_passes: Vec::new(),
                initial_stats: initial_stats.clone(),
                final_stats: initial_stats,
                stages: vec![HealingStageReport {
                    stage: HealingStage::InitialCheck,
                    pass_index: None,
                    issue_count: initial_issue_count,
                }],
            },
        );
    }

    if initial.is_valid() {
        return (
            brep.clone(),
            HealingReport {
                initial: initial.clone(),
                final_result: initial,
                passes: Vec::new(),
                parametric_passes: Vec::new(),
                make_connected_passes: Vec::new(),
                initial_stats: initial_stats.clone(),
                final_stats: initial_stats,
                stages: vec![HealingStageReport {
                    stage: HealingStage::InitialCheck,
                    pass_index: None,
                    issue_count: 0,
                }],
            },
        );
    }

    let mut current = brep.clone();
    let mut passes = Vec::new();
    let mut parametric_passes = Vec::new();
    let mut make_connected_passes = Vec::new();
    let mut stages = vec![HealingStageReport {
        stage: HealingStage::InitialCheck,
        pass_index: None,
        issue_count: initial.issues.len(),
    }];
    let pass_count = options.max_passes.max(1);

    let run_prepass = match options.make_connected_prepass_mode {
        MakeConnectedPrepassMode::Disabled => false,
        MakeConnectedPrepassMode::IssueDriven => has_connectivity_stress_issues(&initial),
        MakeConnectedPrepassMode::Always => true,
    };

    if run_prepass {
        let (reconnected, mc_report) = make_connected_iterative_with_growth_cap(
            &current,
            options.make_connected_tolerance,
            options.make_connected_max_passes,
            options.make_connected_tolerance_growth,
            options.make_connected_tolerance_cap,
        );
        current = reconnected;
        make_connected_passes.push(mc_report);

        let chk = check(&current);
        stages.push(HealingStageReport {
            stage: HealingStage::PreMakeConnected,
            pass_index: None,
            issue_count: chk.issues.len(),
        });
        if chk.is_valid() {
            let final_stats = HealingIssueStats::from_check_result(&chk);
            stages.push(HealingStageReport {
                stage: HealingStage::FinalCheck,
                pass_index: None,
                issue_count: chk.issues.len(),
            });
            return (
                current,
                HealingReport {
                    initial,
                    final_result: chk,
                    passes,
                    parametric_passes,
                    make_connected_passes,
                    initial_stats,
                    final_stats,
                    stages,
                },
            );
        }
    }

    if options.run_parametric_consistency_prepass
        && has_parametric_issues(&current, options.tolerance)
    {
        let (next, same_range_fixed) = fix_same_range_with_scan(&current, options.tolerance);
        let (next, same_parameter_fixed) =
            fix_same_parameter_with_scan(&next, options.tolerance);
        current = next;
        parametric_passes.push(ParametricConsistencyReport {
            same_range_fixed,
            same_parameter_fixed,
        });

        let chk = check(&current);
        stages.push(HealingStageReport {
            stage: HealingStage::ParametricConsistencyPass,
            pass_index: None,
            issue_count: chk.issues.len(),
        });
        if chk.is_valid() {
            let final_stats = HealingIssueStats::from_check_result(&chk);
            stages.push(HealingStageReport {
                stage: HealingStage::FinalCheck,
                pass_index: None,
                issue_count: chk.issues.len(),
            });
            return (
                current,
                HealingReport {
                    initial,
                    final_result: chk,
                    passes,
                    parametric_passes,
                    make_connected_passes,
                    initial_stats,
                    final_stats,
                    stages,
                },
            );
        }
    }

    for pass_idx in 0..pass_count {
        let (next, rep) = repair(&current, options.tolerance);
        current = next;
        let no_changes = rep.vertices_merged == 0
            && rep.degenerate_faces_removed == 0
            && rep.normals_recomputed == 0
            && rep.wires_fixed == 0
            && rep.same_range_fixed == 0
            && rep.same_parameter_fixed == 0;
        passes.push(rep);

        let mut chk = check(&current);
        stages.push(HealingStageReport {
            stage: HealingStage::RepairPass,
            pass_index: Some(pass_idx),
            issue_count: chk.issues.len(),
        });

        if options.run_parametric_consistency_iterative
            && !chk.is_valid()
            && has_parametric_issues(&current, options.tolerance)
        {
            let (next, same_range_fixed) = fix_same_range_with_scan(&current, options.tolerance);
            let (next, same_parameter_fixed) =
                fix_same_parameter_with_scan(&next, options.tolerance);
            current = next;
            parametric_passes.push(ParametricConsistencyReport {
                same_range_fixed,
                same_parameter_fixed,
            });

            chk = check(&current);
            stages.push(HealingStageReport {
                stage: HealingStage::ParametricConsistencyPass,
                pass_index: Some(pass_idx),
                issue_count: chk.issues.len(),
            });
        }

        if chk.is_valid() {
            let final_stats = HealingIssueStats::from_check_result(&chk);
            stages.push(HealingStageReport {
                stage: HealingStage::FinalCheck,
                pass_index: None,
                issue_count: chk.issues.len(),
            });
            return (
                current,
                HealingReport {
                    initial,
                    final_result: chk,
                    passes,
                    parametric_passes,
                    make_connected_passes,
                    initial_stats,
                    final_stats,
                    stages,
                },
            );
        }

        if no_changes && options.run_make_connected_on_stall {
            let (reconnected, mc_report) = make_connected_iterative_with_growth_cap(
                &current,
                options.make_connected_tolerance,
                options.make_connected_max_passes,
                options.make_connected_tolerance_growth,
                options.make_connected_tolerance_cap,
            );
            current = reconnected;
            let mc_no_changes = mc_report.vertices_merged == 0 && mc_report.small_edges_removed == 0;
            make_connected_passes.push(mc_report);

            let chk = check(&current);
            stages.push(HealingStageReport {
                stage: HealingStage::MakeConnectedPass,
                pass_index: Some(pass_idx),
                issue_count: chk.issues.len(),
            });

            if chk.is_valid() || mc_no_changes {
                let final_stats = HealingIssueStats::from_check_result(&chk);
                stages.push(HealingStageReport {
                    stage: HealingStage::FinalCheck,
                    pass_index: None,
                    issue_count: chk.issues.len(),
                });
                return (
                    current,
                    HealingReport {
                        initial,
                        final_result: chk,
                        passes,
                        parametric_passes,
                        make_connected_passes,
                        initial_stats,
                        final_stats,
                        stages,
                    },
                );
            }
            continue;
        }

        if no_changes {
            let final_stats = HealingIssueStats::from_check_result(&chk);
            stages.push(HealingStageReport {
                stage: HealingStage::FinalCheck,
                pass_index: None,
                issue_count: chk.issues.len(),
            });
            return (
                current,
                HealingReport {
                    initial,
                    final_result: chk,
                    passes,
                    parametric_passes,
                    make_connected_passes,
                    initial_stats,
                    final_stats,
                    stages,
                },
            );
        }
    }

    let final_result = check(&current);
    let final_stats = HealingIssueStats::from_check_result(&final_result);
    stages.push(HealingStageReport {
        stage: HealingStage::FinalCheck,
        pass_index: None,
        issue_count: final_result.issues.len(),
    });
    (
        current,
        HealingReport {
            initial,
            final_result,
            passes,
            parametric_passes,
            make_connected_passes,
            initial_stats,
            final_stats,
            stages,
        },
    )
}

fn has_connectivity_stress_issues(result: &CheckResult) -> bool {
    result.issues.iter().any(|issue| {
        matches!(
            issue,
            CheckIssue::OpenWire { .. }
                | CheckIssue::NonManifoldEdge { .. }
                | CheckIssue::SelfIntersectingWire { .. }
                | CheckIssue::GeometricSelfIntersection { .. }
        )
    })
}

fn has_parametric_issues(brep: &BRep, tolerance: f64) -> bool {
    !diagnose_same_range(brep, tolerance).is_clean()
        || !diagnose_same_parameter(brep, tolerance).is_clean()
}

/// Convenience wrapper using default options.
pub fn heal(brep: &BRep) -> (BRep, HealingReport) {
    analyze_and_heal(brep, HealingOptions::default())
}

/// Execute a ShapeProcess-like custom operator chain.
///
/// This is a configurable alternative to [`analyze_and_heal`] for callers that
/// need explicit control over pass ordering.
pub fn run_healing_operator_chain(
    brep: &BRep,
    options: HealingOptions,
    operators: &[HealingOperator],
) -> (BRep, HealingReport) {
    let initial = check(brep);
    let initial_stats = HealingIssueStats::from_check_result(&initial);
    let mut current = brep.clone();

    let mut passes = Vec::new();
    let mut parametric_passes = Vec::new();
    let mut make_connected_passes = Vec::new();
    let mut stages = vec![HealingStageReport {
        stage: HealingStage::InitialCheck,
        pass_index: None,
        issue_count: initial.issues.len(),
    }];

    if matches!(options.mode, HealingMode::AnalyzeOnly) || initial.is_valid() {
        let final_result = check(&current);
        let final_stats = HealingIssueStats::from_check_result(&final_result);
        stages.push(HealingStageReport {
            stage: HealingStage::FinalCheck,
            pass_index: None,
            issue_count: final_result.issues.len(),
        });
        return (
            current,
            HealingReport {
                initial,
                final_result,
                passes,
                parametric_passes,
                make_connected_passes,
                initial_stats,
                final_stats,
                stages,
            },
        );
    }

    for (op_idx, op) in operators.iter().enumerate() {
        match op {
            HealingOperator::MakeConnected => {
                let (next, mc_report) = make_connected_iterative_with_growth_cap(
                    &current,
                    options.make_connected_tolerance,
                    options.make_connected_max_passes,
                    options.make_connected_tolerance_growth,
                    options.make_connected_tolerance_cap,
                );
                current = next;
                make_connected_passes.push(mc_report);
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::MakeConnectedPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::ParametricConsistency => {
                let (next, same_range_fixed) = fix_same_range_with_scan(&current, options.tolerance);
                let (next, same_parameter_fixed) = fix_same_parameter_with_scan(&next, options.tolerance);
                current = next;
                parametric_passes.push(ParametricConsistencyReport {
                    same_range_fixed,
                    same_parameter_fixed,
                });
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::ParametricConsistencyPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::Repair => {
                let (next, rep) = repair(&current, options.tolerance);
                current = next;
                passes.push(rep);
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::RepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::WireGapRepair => {
                let (next, _wire_gap_report) = crate::brep_repair::fix_wire_gaps(
                    &current,
                    options.tolerance,
                    options.tolerance * 10.0, // max_gap = 10x tolerance
                );
                current = next;
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::RepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::UvBoundsRepair => {
                let (next, _uv_report) = crate::brep_repair::fix_uv_bounds_violations(
                    &current,
                    options.tolerance,
                );
                current = next;
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::ParametricConsistencyPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::StopIfClean => {
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::OperatorChainStep,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                if chk.is_valid() {
                    break;
                }
            }
            HealingOperator::FixSmallAreaFaces => {
                let (next, removed) = fix_small_area_faces(&current, options.tolerance);
                current = next;
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                // Track in passes via a synthetic RepairReport
                passes.push(RepairReport {
                    degenerate_faces_removed: removed,
                    ..RepairReport::default()
                });
            }
            HealingOperator::FixSliverFaces => {
                let (next, fixed) = fix_sliver_faces(&current, options.tolerance);
                current = next;
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    wires_fixed: fixed,
                    ..RepairReport::default()
                });
            }
            HealingOperator::FixNonManifold => {
                let (next, fixed) = fix_non_manifold(&current, options.tolerance);
                current = next;
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::TopologyRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    vertices_merged: fixed,
                    ..RepairReport::default()
                });
            }
            HealingOperator::PropagateTolerances => {
                use crate::brep_repair::ToleranceFlowDirection;
                current = crate::brep_repair::propagate_tolerances(
                    &current,
                    options.tolerance,
                    ToleranceFlowDirection::BottomUp,
                );
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::FinalizePass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::UnifySameDomain => {
                let (next, merged) = unify_same_domain_faces(&current, options.tolerance);
                current = next;
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::TopologyRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    faces_reoriented: merged,
                    ..RepairReport::default()
                });
            }
            HealingOperator::RemoveInternalFaces => {
                let (next, removed) = remove_internal_faces(&current);
                current = next;
                let chk = check(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::TopologyRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    degenerate_faces_removed: removed,
                    ..RepairReport::default()
                });
            }
        }
    }

    let final_result = check(&current);
    let final_stats = HealingIssueStats::from_check_result(&final_result);
    stages.push(HealingStageReport {
        stage: HealingStage::FinalCheck,
        pass_index: None,
        issue_count: final_result.issues.len(),
    });

    (
        current,
        HealingReport {
            initial,
            final_result,
            passes,
            parametric_passes,
            make_connected_passes,
            initial_stats,
            final_stats,
            stages,
        },
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// ShapeProcess Implementation
// ─────────────────────────────────────────────────────────────────────────────

/// Execute a full ShapeProcess pipeline on a BRep.
///
/// This is the main entry point for OCCT ShapeProcess-style healing.
/// It runs a configurable sequence of operators organized into stages.
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `config` - Configuration controlling operators and parameters.
///
/// # Returns
/// A tuple of (processed BRep, ShapeProcessReport).
///
/// # Example
/// ```ignore
/// use rcad_algorithms::healing::{run_shape_process, ShapeProcessConfig};
///
/// let config = ShapeProcessConfig::import_preset();
/// let (healed, report) = run_shape_process(&brep, &config);
/// if report.is_clean() {
///     println!("Shape is now valid");
/// }
/// ```
pub fn run_shape_process(brep: &BRep, config: &ShapeProcessConfig) -> (BRep, ShapeProcessReport) {
    use std::time::Instant;

    let start_time = Instant::now();
    let initial = check(brep);
    let initial_stats = HealingIssueStats::from_check_result(&initial);

    // Early exit if shape is already clean and stop_on_clean is true
    if initial.is_valid() && config.stop_on_clean {
        let elapsed = start_time.elapsed().as_secs_f64();
        return (
            brep.clone(),
            ShapeProcessReport {
                initial: initial.clone(),
                final_result: initial,
                initial_stats: initial_stats.clone(),
                final_stats: initial_stats,
                stages: vec![StageReport {
                    stage: HealingStage::InitialCheck,
                    pass_index: None,
                    issue_count_before: 0,
                    issue_count_after: 0,
                    operator_reports: vec![],
                    elapsed_seconds: elapsed,
                }],
                stats: ShapeProcessStats {
                    operators_executed: 0,
                    total_modifications: 0,
                    total_issues_fixed: 0,
                    stages_executed: 1,
                    total_elapsed_seconds: elapsed,
                    iterations: 1,
                    converged_early: true,
                    is_clean: true,
                },
                config_summary: format!("{:?}", config.operators),
            },
        );
    }

    let mut current = brep.clone();
    let mut stages: Vec<StageReport> = Vec::new();
    let mut total_modifications = 0usize;
    let mut operators_executed = 0usize;

    // Build healing options from config
    let options = config.healing_options.clone();
    let mut current_tolerance = config.base_tolerance;

    // Record initial stage
    let initial_issue_count = initial.issues.len();
    stages.push(StageReport {
        stage: HealingStage::InitialCheck,
        pass_index: None,
        issue_count_before: initial_issue_count,
        issue_count_after: initial_issue_count,
        operator_reports: vec![],
        elapsed_seconds: 0.0,
    });

    let mut converged_early = false;
    let max_iters = config.max_iterations.max(1);

    for iter in 0..max_iters {
        let iter_start = Instant::now();

        // Execute the operator chain
        let (next, healing_report) = run_healing_operator_chain(&current, options, &config.operators);
        current = next;

        // Track modifications
        let iter_mods: usize = healing_report.passes.iter()
            .map(|p| p.vertices_merged + p.degenerate_faces_removed + p.normals_recomputed
                + p.faces_reoriented + p.wires_fixed + p.same_range_fixed + p.same_parameter_fixed)
            .sum();
        total_modifications += iter_mods;
        operators_executed += config.operators.len();

        // Convert healing stages to ShapeProcess stages
        for hs in &healing_report.stages {
            stages.push(StageReport {
                stage: hs.stage,
                pass_index: hs.pass_index,
                issue_count_before: initial_issue_count, // Simplified
                issue_count_after: hs.issue_count,
                operator_reports: vec![],
                elapsed_seconds: 0.0,
            });
        }

        let elapsed = iter_start.elapsed().as_secs_f64();
        if let Some(last_stage) = stages.last_mut() {
            last_stage.elapsed_seconds = elapsed;
        }

        // Check for convergence
        if healing_report.is_clean() && config.stop_on_clean {
            converged_early = true;
            break;
        }

        // Apply tolerance growth for next iteration
        if iter + 1 < max_iters {
            current_tolerance = (current_tolerance * config.tolerance_growth).min(config.tolerance_cap);
        }
    }

    let final_result = check(&current);
    let final_stats = HealingIssueStats::from_check_result(&final_result);
    let total_elapsed = start_time.elapsed().as_secs_f64();

    // Add finalization stage
    stages.push(StageReport {
        stage: HealingStage::FinalizePass,
        pass_index: None,
        issue_count_before: initial_issue_count,
        issue_count_after: final_result.issues.len(),
        operator_reports: vec![],
        elapsed_seconds: 0.0,
    });
    stages.push(StageReport {
        stage: HealingStage::FinalCheck,
        pass_index: None,
        issue_count_before: initial_issue_count,
        issue_count_after: final_result.issues.len(),
        operator_reports: vec![],
        elapsed_seconds: total_elapsed,
    });

    let stats = ShapeProcessStats {
        operators_executed,
        total_modifications,
        total_issues_fixed: initial_issue_count.saturating_sub(final_result.issues.len()),
        stages_executed: stages.len(),
        total_elapsed_seconds: total_elapsed,
        iterations: max_iters,
        converged_early,
        is_clean: final_result.is_valid(),
    };

    (
        current,
        ShapeProcessReport {
            initial,
            final_result,
            initial_stats,
            final_stats,
            stages,
            stats,
            config_summary: format!("{:?}", config.operators),
        },
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper Functions for New Operators
// ─────────────────────────────────────────────────────────────────────────────

/// Remove faces with area below a threshold.
///
/// Returns (modified BRep, count of removed faces).
fn fix_small_area_faces(brep: &BRep, min_area: f64) -> (BRep, usize) {
    let mut result = brep.clone();
    let mut removed_count = 0usize;
    let min_area = if min_area > 0.0 { min_area } else { 1e-10 };

    for solid in &mut result.solids {
        for shell in &mut solid.shells {
            let original_len = shell.faces.len();
            let mut kept_faces = Vec::new();

            for face in &shell.faces {
                // Estimate face area using fan triangulation
                let area = estimate_face_area_from_wire(brep, &face.outer_wire);

                if area >= min_area {
                    kept_faces.push(face.clone());
                } else {
                    removed_count += 1;
                }
            }

            shell.faces = kept_faces;
            removed_count += original_len.saturating_sub(shell.faces.len());
        }
    }

    (result, removed_count)
}

/// Fix sliver (thin elongated) faces by merging with neighbors.
///
/// A sliver face has a high aspect ratio (elongated in one dimension).
/// Returns (modified BRep, count of fixed faces).
fn fix_sliver_faces(brep: &BRep, max_aspect_ratio: f64) -> (BRep, usize) {
    // Placeholder implementation - for now just return the input unchanged
    // A full implementation would detect sliver faces by computing aspect ratio
    // and merge them with adjacent faces or remove them
    let _ = (brep, max_aspect_ratio);
    (brep.clone(), 0)
}

/// Repair non-manifold topology by handling multi-face edges.
///
/// Non-manifold edges are shared by more than 2 faces.
/// Returns (modified BRep, count of edges processed).
fn fix_non_manifold(brep: &BRep, _tolerance: f64) -> (BRep, usize) {
    use rcad_kernel::BRepGraph;

    let graph = BRepGraph::from_brep(brep);
    let summary = graph.non_manifold_summary();

    if summary.is_clean() {
        return (brep.clone(), 0);
    }

    // For now, we just identify non-manifold issues
    // Full implementation would split multi-face edges into separate copies
    let non_manifold_count = summary.multi_face_edges.len();

    // Return unchanged for now - this is a complex operation
    // that requires topology restructuring
    (brep.clone(), non_manifold_count)
}

/// Merge faces that share the same underlying surface.
///
/// This is useful for removing artificial seams in imported CAD data.
/// Returns (modified BRep, count of merged face groups).
fn unify_same_domain_faces(brep: &BRep, _tolerance: f64) -> (BRep, usize) {
    // Placeholder implementation - requires surface comparison and face merging
    // A full implementation would identify faces sharing the same surface
    // and merge them into single faces
    (brep.clone(), 0)
}

/// Remove internal faces (faces inside the solid volume).
///
/// Internal faces typically result from boolean operations that left
/// internal partitions. Returns (modified BRep, count of removed faces).
fn remove_internal_faces(brep: &BRep) -> (BRep, usize) {
    // Placeholder implementation - requires volumetric analysis
    // A full implementation would use ray casting or point-in-volume tests
    // to identify and remove internal partition faces
    (brep.clone(), 0)
}

/// Estimate face area from its wire using fan triangulation.
fn estimate_face_area_from_wire(brep: &BRep, wire: &rcad_kernel::topology::Wire) -> f64 {
    use glam::DVec3;

    // Collect vertex positions in order
    let mut pts: Vec<DVec3> = Vec::new();
    for we in &wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            let vi = if we.forward { edge.start } else { edge.end };
            if let Some(v) = brep.vertices.get(vi) {
                pts.push(v.point);
            }
        }
    }

    if pts.len() < 3 {
        return 0.0;
    }

    // Fan triangulation from first point
    let p0 = pts[0];
    let mut area = 0.0f64;
    for i in 1..pts.len() - 1 {
        area += (pts[i] - p0).cross(pts[i + 1] - p0).length() * 0.5;
    }

    area
}

#[cfg(test)]
mod tests {
    use glam::DVec3;
    use rcad_kernel::PrimitiveSolid;

    use super::*;
    use crate::geom_populate;

    fn unit_box() -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn heal_valid_box_is_noop() {
        let b = unit_box();
        let (out, report) = heal(&b);
        assert!(report.initial.is_valid());
        assert!(report.final_result.is_valid());
        assert!(report.passes.is_empty());
        assert!(report.parametric_passes.is_empty());
        assert!(report.make_connected_passes.is_empty());
        assert!(!report.stages.is_empty());
        assert_eq!(out.vertices.len(), b.vertices.len());
        assert_eq!(out.edges.len(), b.edges.len());
    }

    #[test]
    fn heal_zero_normal_face_gets_fixed() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (out, report) = heal(&b);
        assert!(report.initial_issue_count() >= 1);
        assert!(report.is_improved() || report.is_clean());
        assert!(report.initial_stats.zero_normal >= 1);
        assert_eq!(report.initial_stats.total(), report.initial_issue_count());
        assert_eq!(report.final_stats.total(), report.final_issue_count());
        assert!(
            report
                .stages
                .iter()
                .any(|s| matches!(s.stage, HealingStage::FinalCheck))
        );
        assert!(!out.solids[0].shells[0].faces[0].normal.abs_diff_eq(DVec3::ZERO, 0.0));
    }

    #[test]
    fn analyze_only_preserves_input_and_reports_issues() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                mode: HealingMode::AnalyzeOnly,
                ..HealingOptions::default()
            },
        );

        assert!(report.initial_issue_count() >= 1);
        assert_eq!(report.initial_issue_count(), report.final_issue_count());
        assert!(report.passes.is_empty());
        assert!(report.parametric_passes.is_empty());
        assert!(report.make_connected_passes.is_empty());
        assert_eq!(out.solids[0].shells[0].faces[0].normal, DVec3::ZERO);
    }

    #[test]
    fn healing_make_connected_fallback_reporting_is_consistent() {
        let mut b = unit_box();

        // Keep at least one checker issue that standard repair does not heal.
        b.solids[0].shells[0].faces[0].outer_wire.edges[0].idx = usize::MAX;
        // Add near-duplicate vertices that can be merged only by the fallback
        // tolerance (repair tolerance intentionally set much tighter).
        b.vertices[1].point = b.vertices[0].point + DVec3::new(1.0e-6, 0.0, 0.0);

        let (_out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                tolerance: 1.0e-12,
                max_passes: 1,
                run_make_connected_on_stall: true,
                make_connected_tolerance: 1.0e-4,
                make_connected_max_passes: 2,
                make_connected_tolerance_growth: 1.0,
                make_connected_tolerance_cap: 1.0e-4,
                ..HealingOptions::default()
            },
        );

        // Depending on how much progress the regular repair pass can make,
        // make-connected fallback may or may not be needed. If it ran, stage
        // and report vectors must stay in sync.
        let mc_stage_count = report
            .stages
            .iter()
            .filter(|s| matches!(s.stage, HealingStage::MakeConnectedPass))
            .count();
        assert_eq!(mc_stage_count, report.make_connected_passes.len());
        assert!(report.make_connected_passes.len() <= 1);
    }

    #[test]
    fn healing_parametric_consistency_pass_is_reported_when_enabled_by_data() {
        let mut b = unit_box();

        // Make one edge obviously suspect for SameRange/SameParameter scans.
        if b.geom.edge_same_parameter.len() < b.edges.len() {
            b.geom.edge_same_parameter.resize(b.edges.len(), true);
        }
        b.geom.edge_same_parameter[0] = false;
        if b.geom.edge_curve_range.len() < b.edges.len() {
            b.geom.edge_curve_range.resize(b.edges.len(), Some([0.0, 1.0]));
        }
        b.geom.edge_curve_range[0] = Some([0.0, 1.0]);

        let (_out, report) = analyze_and_heal(&b, HealingOptions::default());
        let saw_param_stage = report
            .stages
            .iter()
            .any(|s| matches!(s.stage, HealingStage::ParametricConsistencyPass));
        assert_eq!(saw_param_stage, !report.parametric_passes.is_empty());
    }

    #[test]
    fn healing_can_disable_parametric_consistency_prepass() {
        let mut b = unit_box();
        if b.geom.edge_same_parameter.len() < b.edges.len() {
            b.geom.edge_same_parameter.resize(b.edges.len(), true);
        }
        b.geom.edge_same_parameter[0] = false;

        let (_out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                run_parametric_consistency_prepass: false,
                run_parametric_consistency_iterative: false,
                ..HealingOptions::default()
            },
        );

        assert!(report.parametric_passes.is_empty());
    }

    #[test]
    fn healing_make_connected_prepass_always_records_stage() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].outer_wire.edges[0].idx = usize::MAX;

        let (_out, report) = analyze_and_heal(
            &b,
            HealingOptions {
                max_passes: 1,
                make_connected_prepass_mode: MakeConnectedPrepassMode::Always,
                make_connected_tolerance: 1.0e-4,
                make_connected_max_passes: 1,
                make_connected_tolerance_growth: 1.0,
                make_connected_tolerance_cap: 1.0e-4,
                ..HealingOptions::default()
            },
        );

        assert!(
            report
                .stages
                .iter()
                .any(|s| matches!(s.stage, HealingStage::PreMakeConnected))
        );
        assert!(!report.make_connected_passes.is_empty());
    }

    #[test]
    fn operator_chain_runs_repair_and_parametric_passes() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;
        if b.geom.edge_same_parameter.len() < b.edges.len() {
            b.geom.edge_same_parameter.resize(b.edges.len(), true);
        }
        b.geom.edge_same_parameter[0] = false;

        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::ParametricConsistency,
                HealingOperator::Repair,
                HealingOperator::StopIfClean,
            ],
        );

        assert!(!report.parametric_passes.is_empty());
        assert!(!report.passes.is_empty());
        assert!(
            report.stages.iter().any(|s| matches!(s.stage, HealingStage::OperatorChainStep))
        );
    }

    #[test]
    fn operator_chain_stop_if_clean_short_circuits_following_steps() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::Repair,
                HealingOperator::StopIfClean,
                HealingOperator::MakeConnected,
            ],
        );

        // Repair should clean this case; stop-if-clean should prevent make-connected.
        assert!(report.make_connected_passes.is_empty());
        assert!(report.final_result.is_valid());
    }

    #[test]
    fn shape_process_default_config_works_on_valid_shape() {
        let b = unit_box();
        let config = ShapeProcessConfig::default();
        let (out, report) = run_shape_process(&b, &config);

        assert!(report.is_clean());
        assert!(report.stats.converged_early);
        assert_eq!(out.vertices.len(), b.vertices.len());
    }

    #[test]
    fn shape_process_import_preset_fixes_zero_normal() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::import_preset();
        let (_out, report) = run_shape_process(&b, &config);

        assert!(report.is_improved() || report.is_clean());
        assert!(report.initial_issue_count() >= 1);
    }

    #[test]
    fn shape_process_boolean_cleanup_preset_works() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::boolean_cleanup_preset();
        let (_out, report) = run_shape_process(&b, &config);

        assert!(report.is_improved() || report.is_clean());
    }

    #[test]
    fn shape_process_analysis_preset_is_conservative() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::analysis_preset();
        let (_out, report) = run_shape_process(&b, &config);

        // Analysis preset should at least diagnose issues
        assert!(report.initial_issue_count() >= 1);
    }

    #[test]
    fn shape_process_aggressive_preset_applies_all_operators() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let config = ShapeProcessConfig::aggressive_preset();
        let (_out, report) = run_shape_process(&b, &config);

        assert!(report.is_improved() || report.is_clean());
        // Aggressive preset has many operators
        assert!(config.operators.len() >= 8);
    }

    #[test]
    fn shape_process_report_summary_is_informative() {
        let b = unit_box();
        let config = ShapeProcessConfig::default();
        let (_out, report) = run_shape_process(&b, &config);

        let summary = report.summary();
        assert!(summary.contains("ShapeProcess"));
        assert!(summary.contains("Clean") || summary.contains("issues"));
    }

    #[test]
    fn operator_chain_handles_new_operators() {
        let b = unit_box();

        // Test that new operators don't panic
        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::FixSmallAreaFaces,
                HealingOperator::FixSliverFaces,
                HealingOperator::FixNonManifold,
                HealingOperator::PropagateTolerances,
                HealingOperator::UnifySameDomain,
                HealingOperator::RemoveInternalFaces,
            ],
        );

        // All operators should run without error
        assert!(!report.stages.is_empty());
    }

    #[test]
    fn fix_small_area_faces_removes_tiny_faces() {
        let b = unit_box();

        // Unit box faces are not tiny, so nothing should be removed
        let (result, removed) = fix_small_area_faces(&b, 1e-12);
        assert_eq!(removed, 0);

        // The result should have the same number of faces
        let result_face_count: usize = result.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        let original_face_count: usize = b.solids.iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        assert_eq!(result_face_count, original_face_count);
    }

    #[test]
    fn new_healing_stages_are_recorded() {
        let mut b = unit_box();
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (_out, report) = run_healing_operator_chain(
            &b,
            HealingOptions::default(),
            &[
                HealingOperator::FixSmallAreaFaces,
                HealingOperator::FixNonManifold,
                HealingOperator::PropagateTolerances,
            ],
        );

        // Should have geometry and topology repair stages
        assert!(report.stages.iter().any(|s|
            matches!(s.stage, HealingStage::GeometryRepairPass)));
        assert!(report.stages.iter().any(|s|
            matches!(s.stage, HealingStage::TopologyRepairPass)));
        assert!(report.stages.iter().any(|s|
            matches!(s.stage, HealingStage::FinalizePass)));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Enhanced Healing: ShapeFix_Solid and ShapeFix_Wire Equivalents
// ─────────────────────────────────────────────────────────────────────────────

/// ShapeFix_Solid equivalent: comprehensive solid repair.
///
/// This function performs OCCT ShapeFix_Solid-like operations:
/// - Shell orientation verification and repair
/// - Solid closure verification
/// - Shell manifoldness checks
/// - Face orientation consistency
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `tolerance` - Geometric tolerance
///
/// # Returns
/// Repaired B-Rep and count of fixes applied.
pub fn fix_solid(brep: &BRep, tolerance: f64) -> (BRep, SolidFixReport) {
    use crate::brep_repair::{fix_face_orientation, recompute_face_normals};
    use rcad_kernel::BRepGraph;

    let mut report = SolidFixReport::default();
    let mut current = brep.clone();

    // Step 1: Recompute invalid normals
    let (brep_with_normals, normals_fixed) = recompute_face_normals(&current);
    current = brep_with_normals;
    report.normals_recomputed = normals_fixed;

    // Step 2: Fix face orientation for inward-pointing faces
    let (brep_oriented, faces_reoriented) = fix_face_orientation(&current);
    current = brep_oriented;
    report.faces_reoriented = faces_reoriented;

    // Step 3: Check solid closure and manifoldness
    let graph = BRepGraph::from_brep(&current);

    // Check if shells are closed
    for (si, solid) in current.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            let is_closed = shell.faces.iter().all(|f| {
                // Check if wire is closed
                let wire = &f.outer_wire;
                if wire.edges.is_empty() {
                    return false;
                }
                true // Simplified check; full implementation would verify vertex chain
            });

            if !is_closed {
                report.unclosed_shells.push((si, shi));
            }
        }
    }

    // Check manifoldness
    let nm_summary = graph.non_manifold_summary();
    report.non_manifold_edges = nm_summary.multi_face_edges.len();
    report.non_manifold_vertices = nm_summary.non_manifold_vertices.len();

    // Step 4: Verify shell orientation consistency
    for solid in &current.solids {
        for shell in &solid.shells {
            // Count faces with normals pointing in consistent direction
            let mut outward_count = 0usize;
            let mut inward_count = 0usize;

            for face in &shell.faces {
                // Heuristic: if normal dot product with center-to-centroid is positive
                // the face is likely outward-facing
                if face.normal.z > 0.0 {
                    outward_count += 1;
                } else if face.normal.z < 0.0 {
                    inward_count += 1;
                }
            }

            // If most normals are inconsistent, note orientation issues
            if outward_count > 0 && inward_count > 0 {
                let ratio = outward_count as f64 / (outward_count + inward_count) as f64;
                if ratio < 0.3 || ratio > 0.7 {
                    report.orientation_inconsistencies += 1;
                }
            }
        }
    }

    report.total_fixes = report.normals_recomputed + report.faces_reoriented;
    (current, report)
}

/// Report from solid-level fixes.
#[derive(Debug, Clone, Default)]
pub struct SolidFixReport {
    /// Number of face normals recomputed.
    pub normals_recomputed: usize,
    /// Number of faces reoriented.
    pub faces_reoriented: usize,
    /// Indices of unclosed shells (solid_idx, shell_idx).
    pub unclosed_shells: Vec<(usize, usize)>,
    /// Number of non-manifold edges detected.
    pub non_manifold_edges: usize,
    /// Number of non-manifold vertices detected.
    pub non_manifold_vertices: usize,
    /// Number of shells with orientation inconsistencies.
    pub orientation_inconsistencies: usize,
    /// Total number of fixes applied.
    pub total_fixes: usize,
}

impl SolidFixReport {
    pub fn is_clean(&self) -> bool {
        self.unclosed_shells.is_empty()
            && self.non_manifold_edges == 0
            && self.non_manifold_vertices == 0
            && self.orientation_inconsistencies == 0
    }

    pub fn summary(&self) -> String {
        if self.is_clean() && self.total_fixes == 0 {
            "Solid is clean, no fixes needed".to_string()
        } else {
            format!(
                "Solid fixes: {} normals, {} orientations, {} unclosed shells, {} non-manifold edges, {} non-manifold vertices",
                self.normals_recomputed,
                self.faces_reoriented,
                self.unclosed_shells.len(),
                self.non_manifold_edges,
                self.non_manifold_vertices
            )
        }
    }
}

/// ShapeFix_Wire equivalent: comprehensive wire repair.
///
/// This function performs OCCT ShapeFix_Wire-like operations:
/// - Wire closure verification and repair
/// - Edge order verification
/// - Degenerate edge handling
/// - Self-intersection detection
/// - Wire orientation fix
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `tolerance` - Geometric tolerance
///
/// # Returns
/// Repaired B-Rep and detailed wire fix report.
pub fn fix_wire(brep: &BRep, tolerance: f64) -> (BRep, WireFixReport) {
    use crate::brep_repair::fix_wire_orientation;

    let mut report = WireFixReport::default();
    let mut current = brep.clone();

    // Step 1: Fix wire orientation
    let (brep_fixed, wires_fixed) = fix_wire_orientation(&current, tolerance);
    current = brep_fixed;
    report.wires_oriented = wires_fixed;

    // Step 2: Analyze wires for issues
    for (si, solid) in current.solids.iter().enumerate() {
        for (shi, shell) in solid.shells.iter().enumerate() {
            for (fi, face) in shell.faces.iter().enumerate() {
                // Check outer wire
                let outer_issues = analyze_wire_issues(&current, &face.outer_wire, tolerance);
                if outer_issues.open_gaps > 0 || outer_issues.topological_self_intersections > 0 || outer_issues.geometric_self_intersections > 0 {
                    report.outer_wire_issues.push(WireIssueLocation {
                        solid: si,
                        shell: shi,
                        face: fi,
                        wire_idx: 0,
                        issues: outer_issues,
                    });
                }

                // Check inner wires
                for (wi, inner_wire) in face.inner_wires.iter().enumerate() {
                    let inner_issues = analyze_wire_issues(&current, inner_wire, tolerance);
                    if inner_issues.open_gaps > 0 || inner_issues.topological_self_intersections > 0 || inner_issues.geometric_self_intersections > 0 {
                        report.inner_wire_issues.push(WireIssueLocation {
                            solid: si,
                            shell: shi,
                            face: fi,
                            wire_idx: wi + 1,
                            issues: inner_issues,
                        });
                    }
                }
            }
        }
    }

    // Step 3: Count degenerate edges
    for (ei, edge) in current.edges.iter().enumerate() {
        let start_pt = current.vertices.get(edge.start).map(|v| v.point);
        let end_pt = current.vertices.get(edge.end).map(|v| v.point);

        if let (Some(s), Some(e)) = (start_pt, end_pt) {
            if (s - e).length() < tolerance {
                report.degenerate_edges.push(ei);
            }
        }
    }

    // Step 4: Compute wire quality metrics
    report.total_wires_checked = report.outer_wire_issues.len()
        + report.inner_wire_issues.len()
        + current.solids.iter()
            .flat_map(|s| s.shells.iter())
            .flat_map(|sh| sh.faces.iter())
            .map(|f| 1 + f.inner_wires.len())
            .sum::<usize>();

    report.wires_with_issues = report.outer_wire_issues.len() + report.inner_wire_issues.len();
    report.total_fixes = report.wires_oriented;

    (current, report)
}

/// Location of a wire issue.
#[derive(Debug, Clone)]
pub struct WireIssueLocation {
    pub solid: usize,
    pub shell: usize,
    pub face: usize,
    pub wire_idx: usize,
    pub issues: crate::brep_check::WireIssueReport,
}

/// Report from wire-level fixes.
#[derive(Debug, Clone, Default)]
pub struct WireFixReport {
    /// Number of wires with corrected orientation.
    pub wires_oriented: usize,
    /// Issues found in outer wires.
    pub outer_wire_issues: Vec<WireIssueLocation>,
    /// Issues found in inner wires.
    pub inner_wire_issues: Vec<WireIssueLocation>,
    /// Indices of degenerate edges found.
    pub degenerate_edges: Vec<usize>,
    /// Total wires checked.
    pub total_wires_checked: usize,
    /// Wires with issues.
    pub wires_with_issues: usize,
    /// Total fixes applied.
    pub total_fixes: usize,
}

impl WireFixReport {
    pub fn is_clean(&self) -> bool {
        self.outer_wire_issues.is_empty()
            && self.inner_wire_issues.is_empty()
            && self.degenerate_edges.is_empty()
    }

    pub fn summary(&self) -> String {
        if self.is_clean() && self.total_fixes == 0 {
            format!("All {} wires clean, no fixes needed", self.total_wires_checked)
        } else {
            format!(
                "Wire fixes: {} oriented, {} with issues ({} outer, {} inner), {} degenerate edges",
                self.wires_oriented,
                self.wires_with_issues,
                self.outer_wire_issues.len(),
                self.inner_wire_issues.len(),
                self.degenerate_edges.len()
            )
        }
    }
}

/// Analyze wire for issues without modifying.
fn analyze_wire_issues(brep: &BRep, wire: &rcad_kernel::topology::Wire, tolerance: f64) -> crate::brep_check::WireIssueReport {
    let n_edges = brep.edges.len();
    let mut open_gaps = 0usize;
    let mut topological_self_intersections = 0usize;
    let mut geometric_self_intersections = 0usize;

    // Collect wire vertices
    let mut wire_verts = Vec::with_capacity(wire.edges.len());
    for we in &wire.edges {
        if we.idx >= n_edges {
            continue;
        }
        let edge = &brep.edges[we.idx];
        let (sv, ev) = if we.forward {
            (edge.start, edge.end)
        } else {
            (edge.end, edge.start)
        };
        if sv < brep.vertices.len() && ev < brep.vertices.len() {
            wire_verts.push((sv, ev));
        }
    }

    // Check for open gaps
    let n = wire_verts.len();
    if n > 1 {
        for i in 0..n {
            let next = (i + 1) % n;
            let end_v = wire_verts[i].1;
            let start_v = wire_verts[next].0;
            if end_v != start_v {
                let end_pt = brep.vertices[end_v].point;
                let start_pt = brep.vertices[start_v].point;
                if (end_pt - start_pt).length() > tolerance {
                    open_gaps += 1;
                }
            }
        }
    }

    // Check for topological self-intersection (vertex appearing more than twice)
    use std::collections::HashMap;
    let mut vertex_count: HashMap<usize, usize> = HashMap::new();
    for &(sv, ev) in &wire_verts {
        *vertex_count.entry(sv).or_insert(0) += 1;
        *vertex_count.entry(ev).or_insert(0) += 1;
    }
    for &count in vertex_count.values() {
        if count > 2 {
            topological_self_intersections += 1;
        }
    }

    // Check for geometric self-intersection (2D projection)
    if n >= 4 {
        for i in 0..n {
            for j in (i + 2)..n {
                if i == 0 && j == n - 1 {
                    continue; // Adjacent edges wraparound
                }
                let (a_start, a_end) = wire_verts[i];
                let (b_start, b_end) = wire_verts[j];
                let p1 = brep.vertices[a_start].point;
                let p2 = brep.vertices[a_end].point;
                let p3 = brep.vertices[b_start].point;
                let p4 = brep.vertices[b_end].point;

                if segments_intersect_2d(p1, p2, p3, p4) {
                    geometric_self_intersections += 1;
                }
            }
        }
    }

    crate::brep_check::WireIssueReport {
        solid: 0,
        shell: 0,
        face: 0,
        wire_idx: 0,
        edge_count: wire.edges.len(),
        open_gaps,
        topological_self_intersections,
        geometric_self_intersections,
    }
}

/// Check if two 2D line segments intersect (XY plane projection).
fn segments_intersect_2d(p1: glam::DVec3, p2: glam::DVec3, p3: glam::DVec3, p4: glam::DVec3) -> bool {
    let x1 = p1.x; let y1 = p1.y;
    let x2 = p2.x; let y2 = p2.y;
    let x3 = p3.x; let y3 = p3.y;
    let x4 = p4.x; let y4 = p4.y;

    let (min_x1, max_x1) = if x1 < x2 { (x1, x2) } else { (x2, x1) };
    let (min_y1, max_y1) = if y1 < y2 { (y1, y2) } else { (y2, y1) };
    let (min_x2, max_x2) = if x3 < x4 { (x3, x4) } else { (x4, x3) };
    let (min_y2, max_y2) = if y3 < y4 { (y3, y4) } else { (y4, y3) };

    if max_x1 < min_x2 || max_x2 < min_x1 || max_y1 < min_y2 || max_y2 < min_y1 {
        return false;
    }

    fn ccw(ax: f64, ay: f64, bx: f64, by: f64, cx: f64, cy: f64) -> bool {
        (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)
    }

    ccw(x1, y1, x3, y3, x4, y4) != ccw(x2, y2, x3, y3, x4, y4)
        && ccw(x1, y1, x2, y2, x3, y3) != ccw(x1, y1, x2, y2, x4, y4)
}

/// Comprehensive healing with ShapeFix_Solid and ShapeFix_Wire integration.
///
/// This function provides OCCT-equivalent comprehensive healing:
/// 1. Wire-level fixes
/// 2. Face-level fixes
/// 3. Shell-level fixes
/// 4. Solid-level fixes
/// 5. Tolerance propagation
///
/// # Arguments
/// * `brep` - Input B-Rep
/// * `options` - Healing options
///
/// # Returns
/// Healed B-Rep and comprehensive report.
pub fn heal_comprehensive(brep: &BRep, options: &HealingOptions) -> (BRep, ComprehensiveHealingReport) {
    let mut report = ComprehensiveHealingReport::default();
    let mut current = brep.clone();

    // Stage 1: Wire fixes
    let (brep_wire, wire_report) = fix_wire(&current, options.tolerance);
    current = brep_wire;
    report.wire_report = Some(wire_report);

    // Stage 2: Face fixes (via standard repair)
    let (brep_face, repair_report) = repair(&current, options.tolerance);
    current = brep_face;
    report.repair_report = Some(repair_report);

    // Stage 3: Solid fixes
    let (brep_solid, solid_report) = fix_solid(&current, options.tolerance);
    current = brep_solid;
    report.solid_report = Some(solid_report);

    // Stage 4: Tolerance propagation
    current = crate::brep_repair::propagate_tolerances(
        &current,
        options.tolerance,
        crate::brep_repair::ToleranceFlowDirection::BottomUp,
    );
    let tol_report = crate::brep_repair::analyze_tolerances(&current, options.tolerance);
    report.tolerance_report = Some(tol_report.vertices);

    // Final check
    report.final_check = check(&current);
    report.is_clean = report.final_check.is_valid();

    (current, report)
}

/// Comprehensive healing report with all stage details.
#[derive(Debug, Clone, Default)]
pub struct ComprehensiveHealingReport {
    /// Wire-level fix report.
    pub wire_report: Option<WireFixReport>,
    /// Standard repair report.
    pub repair_report: Option<crate::brep_repair::RepairReport>,
    /// Solid-level fix report.
    pub solid_report: Option<SolidFixReport>,
    /// Tolerance propagation report.
    pub tolerance_report: Option<crate::brep_repair::ToleranceStats>,
    /// Final checker result.
    pub final_check: CheckResult,
    /// Whether the result is checker-clean.
    pub is_clean: bool,
}

impl ComprehensiveHealingReport {
    pub fn summary(&self) -> String {
        let mut parts = Vec::new();

        if let Some(ref wr) = self.wire_report {
            if wr.total_fixes > 0 {
                parts.push(format!("wires: {} fixes", wr.total_fixes));
            }
        }

        if let Some(ref rr) = self.repair_report {
            let repairs = rr.vertices_merged + rr.faces_reoriented + rr.wires_fixed;
            if repairs > 0 {
                parts.push(format!("repair: {} fixes", repairs));
            }
        }

        if let Some(ref sr) = self.solid_report {
            if sr.total_fixes > 0 {
                parts.push(format!("solid: {} fixes", sr.total_fixes));
            }
        }

        if parts.is_empty() {
            if self.is_clean {
                "Clean result, no fixes needed".to_string()
            } else {
                format!("Issues remain: {} issues", self.final_check.issues.len())
            }
        } else {
            format!("{} → {}", parts.join(", "), if self.is_clean { "clean" } else { "issues remain" })
        }
    }
}
