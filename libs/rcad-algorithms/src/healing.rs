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
    PreMakeConnected,
    OperatorChainStep,
    ParametricConsistencyPass,
    RepairPass,
    MakeConnectedPass,
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
    /// Stop pipeline execution if the current shape is checker-clean.
    StopIfClean,
}

/// Report for one SameRange/SameParameter consistency pass.
#[derive(Debug, Clone, Default)]
pub struct ParametricConsistencyReport {
    pub same_range_fixed: usize,
    pub same_parameter_fixed: usize,
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
}
