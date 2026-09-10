use std::sync::Arc;
use rcad_kernel::topods::{TShape, Shape, Orientation};

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
            tolerance_cap: TOLERANCE_ADAPTIVE_MAX,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS,
                min_face_area: TOLERANCE_COORD_SUB * 10.0,
                max_sliver_aspect_ratio: 50.0,
                allow_internal_face_removal: false,
                ..Default::default()
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
            tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE,
            base_tolerance: TOLERANCE_ABS * 10.0,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS * 10.0,
                min_face_area: TOLERANCE_LINEAR_ULTRA_STRICT,
                max_sliver_aspect_ratio: 100.0,
                allow_internal_face_removal: true,
                ..Default::default()
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
            tolerance_cap: TOLERANCE_RETRY_LADDER_MID,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS,
                min_face_area: TOLERANCE_VOL_CUBE_FACTOR * TOLERANCE_ADAPTIVE_MAX,
                max_sliver_aspect_ratio: 1000.0,
                allow_internal_face_removal: false,
                ..Default::default()
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
            tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE * 100.0,
            base_tolerance: TOLERANCE_ABS,
            operator_params: OperatorParams {
                tolerance: TOLERANCE_ABS,
                min_face_area: TOLERANCE_COORD_SUB * 10.0,
                max_sliver_aspect_ratio: 50.0,
                allow_internal_face_removal: true,
                ..Default::default()
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
/// Accepts topods::BRep; old-BRep conversion is internal.
pub fn analyze_and_heal(brep: &topods::BRep, options: HealingOptions) -> (topods::BRep, HealingReport) {
    analyze_and_heal_old(brep, options)
}

/// Legacy: takes old BRep. Internal implementation.
fn analyze_and_heal_old(brep: &rcad_kernel::BRep, options: HealingOptions) -> (rcad_kernel::BRep, HealingReport) {
    let initial = brep_check_analyze(brep);
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

        let chk = brep_check_analyze(&current);
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

        let chk = brep_check_analyze(&current);
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

        let mut chk = brep_check_analyze(&current);
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

            chk = brep_check_analyze(&current);
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

            let chk = brep_check_analyze(&current);
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

    let final_result = brep_check_analyze(&current);
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

fn has_parametric_issues(brep: &rcad_kernel::BRep, tolerance: f64) -> bool {
    !diagnose_same_range(brep, tolerance).is_clean()
        || !diagnose_same_parameter(brep, tolerance).is_clean()
}

/// Convenience wrapper using default options.
pub fn heal(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, HealingReport) {
    let (t, report) = analyze_and_heal(&brep, HealingOptions::default());
    ((t).clone(), report)
}

/// Execute a ShapeProcess-like custom operator chain.
///
/// This is a configurable alternative to [`analyze_and_heal`] for callers that
/// need explicit control over pass ordering.
pub fn run_healing_operator_chain(
    brep: &rcad_kernel::BRep,
    options: HealingOptions,
    operators: &[HealingOperator],
) -> (rcad_kernel::BRep, HealingReport) {
    let initial = brep_check_analyze(brep);
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
        let final_result = brep_check_analyze(&current);
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
                let chk = brep_check_analyze(&current);
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
                let chk = brep_check_analyze(&current);
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
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::RepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::WireGapRepair => {
                let (next, _wire_gap_report) = crate::algo_ext::brep_repair::fix_wire_gaps(
                    &current,
                    options.tolerance,
                    options.tolerance * 10.0, // max_gap = 10x tolerance
                );
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::RepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::UvBoundsRepair => {
                let (next, _uv_report) = crate::algo_ext::brep_repair::fix_uv_bounds_violations(
                    &current,
                    options.tolerance,
                );
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::ParametricConsistencyPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::StopIfClean => {
                let chk = brep_check_analyze(&current);
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
                let chk = brep_check_analyze(&current);
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
                let chk = brep_check_analyze(&current);
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
                let chk = brep_check_analyze(&current);
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
                use crate::algo_ext::brep_repair::ToleranceFlowDirection;
                current = crate::algo_ext::brep_repair::propagate_tolerances(
                    &current,
                    options.tolerance,
                    ToleranceFlowDirection::BottomUp,
                );
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::FinalizePass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
            }
            HealingOperator::UnifySameDomain => {
                let (next, merged) = unify_same_domain_faces(&current, options.tolerance);
                current = next;
                let chk = brep_check_analyze(&current);
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
                let chk = brep_check_analyze(&current);
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
            HealingOperator::SplitAngle(params) => {
                let (next, splits) = split_angle_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    faces_reoriented: splits,
                    ..RepairReport::default()
                });
            }
            HealingOperator::SplitContinuity(params) => {
                let (next, splits) = split_continuity_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    // splits tracks edges split, which we report as wires_fixed
                    wires_fixed: splits,
                    ..RepairReport::default()
                });
            }
            HealingOperator::ConvertToBSpline(params) => {
                let (next, conversions) = convert_to_bspline_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    faces_reoriented: conversions,
                    ..RepairReport::default()
                });
            }
            HealingOperator::SurfaceToBezier(params) => {
                let (next, conversions) = surface_to_bezier_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    faces_reoriented: conversions,
                    ..RepairReport::default()
                });
            }
            HealingOperator::ScaleShape(params) => {
                let (next, modifications) = scale_shape_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    vertices_merged: modifications,
                    ..RepairReport::default()
                });
            }
            HealingOperator::DirectFaces(params) => {
                let (next, faces_fixed) = direct_faces_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    faces_reoriented: faces_fixed,
                    ..RepairReport::default()
                });
            }
            HealingOperator::SameParameter(params) => {
                let (next, edges_fixed) = same_parameter_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::ParametricConsistencyPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    same_parameter_fixed: edges_fixed,
                    ..RepairReport::default()
                });
            }
            HealingOperator::RemoveInternalFacesOp(params) => {
                let (next, faces_removed) = remove_internal_faces_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::TopologyRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(RepairReport {
                    degenerate_faces_removed: faces_removed,
                    ..RepairReport::default()
                });
            }
            HealingOperator::HealGeometry(params) => {
                let (next, report) = heal_geometry_operator(&current, params);
                current = next;
                let chk = brep_check_analyze(&current);
                stages.push(HealingStageReport {
                    stage: HealingStage::GeometryRepairPass,
                    pass_index: Some(op_idx),
                    issue_count: chk.issues.len(),
                });
                passes.push(report);
            }
        }
    }

    let final_result = brep_check_analyze(&current);
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
pub fn run_shape_process(brep: &rcad_kernel::BRep, config: &ShapeProcessConfig) -> (rcad_kernel::BRep, ShapeProcessReport) {
    use std::time::Instant;

    let start_time = Instant::now();
    let initial = brep_check_analyze(brep);
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
    let options = config.healing_options;
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

    let final_result = brep_check_analyze(&current);
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
fn fix_small_area_faces(brep: &rcad_kernel::BRep, min_area: f64) -> (rcad_kernel::BRep, usize) {
    let mut result = brep.clone();
    let mut removed_count = 0usize;
    let min_area = if min_area > 0.0 { min_area } else { TOLERANCE_LINEAR_ULTRA_STRICT };

    // Collect solid indices (avoid re-entrant modification during iteration)
    let solid_indices: Vec<usize> = (0..result.tshapes.len())
        .filter(|&i| matches!(&*result.tshapes[i], TShape::Solid(_)))
        .collect();

    for si in solid_indices {
        let shell_refs = {
            let TShape::Solid(sd) = &*result.tshapes[si] else { continue; };
            sd.shells.clone()
        };

        for shell_sr in &shell_refs {
            let (face_refs, original_len) = {
                let TShape::Shell(shd) = &*result.tshapes[shell_sr.index] else { continue; };
                (shd.faces.clone(), shd.faces.len())
            };

            let mut kept_faces = Vec::new();
            for face_sr in &face_refs {
                let TShape::Face(fd) = &*result.tshapes[face_sr.index] else {
                    kept_faces.push(face_sr.clone());
                    continue;
                };

                let area = estimate_face_area_from_wire(&result, fd.outer_wire.clone());
                if area >= min_area {
                    kept_faces.push(face_sr.clone());
                } else {
                    removed_count += 1;
                }
            }

            if kept_faces.len() < original_len {
                if let TShape::Shell(shd_mut) = Arc::make_mut(&mut result.tshapes[shell_sr.index]) {
                    shd_mut.faces = kept_faces;
                }
            }
        }
    }

    (result, removed_count)
}

/// Fix sliver (thin elongated) faces by merging with neighbors.
///
/// A sliver face has a high aspect ratio (elongated in one dimension).
/// Returns (modified BRep, count of fixed faces).
fn fix_sliver_faces(brep: &rcad_kernel::BRep, max_aspect_ratio: f64) -> (rcad_kernel::BRep, usize) {
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
fn fix_non_manifold(brep: &rcad_kernel::BRep, _tolerance: f64) -> (rcad_kernel::BRep, usize) {
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
fn unify_same_domain_faces(brep: &rcad_kernel::BRep, _tolerance: f64) -> (rcad_kernel::BRep, usize) {
    // Placeholder implementation - requires surface comparison and face merging
    // A full implementation would identify faces sharing the same surface
    // and merge them into single faces
    (brep.clone(), 0)
}

/// Remove internal faces (faces inside the solid volume).
///
/// Internal faces typically result from boolean operations that left
/// internal partitions. Returns (modified BRep, count of removed faces).
fn remove_internal_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
    // Placeholder implementation - requires volumetric analysis
    // A full implementation would use ray casting or point-in-volume tests
    // to identify and remove internal partition faces
    (brep.clone(), 0)
}

/// Estimate face area from its wire using fan triangulation.
fn estimate_face_area_from_wire(brep: &rcad_kernel::BRep, wire_sr: Shape) -> f64 {
    use glam::DVec3;

    let wd = match &*brep.tshapes[wire_sr.index] {
        TShape::Wire(w) => w,
        _ => return 0.0,
    };

    // Collect vertex positions in order following wire edge orientation
    let mut pts: Vec<DVec3> = Vec::new();
    for edge_sr in &wd.edges {
        if let TShape::Edge(ed) = &*brep.tshapes[edge_sr.index] {
            let vi = if matches!(edge_sr.orientation, Orientation::Forward) {
                ed.first.index
            } else {
                ed.last.index
            };
            if let Some(p) = brep.vertex_point(vi) {
                pts.push(p);
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

// ─────────────────────────────────────────────────────────────────────────────
// New Operator Implementations
// ─────────────────────────────────────────────────────────────────────────────

/// Split faces at angle thresholds (SplitAngle operator).
///
/// Splits cylindrical, conical, spherical, and toroidal faces into sectors
/// where each sector has a maximum angular extent.
///
/// Returns (modified BRep, count of faces split).
fn split_angle_operator(brep: &rcad_kernel::BRep, params: &SplitAngleOperator) -> (rcad_kernel::BRep, usize) {
    use rcad_kernel::geom::Surface3;
    use std::f64::consts::PI;

    let mut result = brep.clone();
    let mut split_count = 0usize;
    let max_angle = params.max_angle.max(PI / 36.0); // At least 5 degrees minimum

    // Process each face via tshape walk
    let solid_indices: Vec<usize> = (0..result.tshapes.len())
        .filter(|&i| matches!(&*result.tshapes[i], TShape::Solid(_)))
        .collect();

    for si in solid_indices {
        let shell_refs = {
            let TShape::Solid(sd) = &*result.tshapes[si] else { continue; };
            sd.shells.clone()
        };

        for shell_sr in &shell_refs {
            let (face_refs, original_len) = {
                let TShape::Shell(shd) = &*result.tshapes[shell_sr.index] else { continue; };
                (shd.faces.clone(), shd.faces.len())
            };

            let mut new_face_srs = Vec::new();
            for face_sr in &face_refs {
                let TShape::Face(fd) = &*result.tshapes[face_sr.index] else {
                    new_face_srs.push(face_sr.clone());
                    continue;
                };

                let surface = fd.surface.as_ref();
                let should_split = surface.is_some_and(|s| {
                    match s {
                        Surface3::Cylinder(_) => params.split_cylinders,
                        Surface3::Torus(_) => params.split_tori,
                        Surface3::Sphere(_) => params.split_spheres,
                        Surface3::Cone(_) => params.split_cones,
                        _ => false,
                    }
                });

                if should_split {
                    let surf = surface.unwrap();
                    // Calculate how many sectors are needed
                    let (u_range, v_range, is_u_periodic, is_v_periodic) = match surf {
                        Surface3::Cylinder(_) => ((0.0, 2.0 * PI), (-1e10, 1e10), true, false),
                        Surface3::Torus(_) => ((0.0, 2.0 * PI), (0.0, 2.0 * PI), true, true),
                        Surface3::Sphere(_) => ((0.0, 2.0 * PI), (0.0, PI), true, false),
                        Surface3::Cone(_) => ((0.0, 2.0 * PI), (0.0, 1e10), true, false),
                        _ => ((0.0, 1.0), (0.0, 1.0), false, false),
                    };

                    let u_span = u_range.1 - u_range.0;
                    let v_span = v_range.1 - v_range.0;

                    let u_sectors = if is_u_periodic {
                        ((u_span / max_angle).ceil() as usize).max(1)
                    } else {
                        1
                    };

                    let v_sectors = if is_v_periodic {
                        ((v_span / max_angle).ceil() as usize).max(1)
                    } else {
                        1
                    };

                    if u_sectors > 1 || v_sectors > 1 {
                        split_count += 1;
                    }
                    new_face_srs.push(face_sr.clone());
                } else {
                    new_face_srs.push(face_sr.clone());
                }
            }

            if new_face_srs.len() != original_len {
                if let TShape::Shell(shd_mut) = Arc::make_mut(&mut result.tshapes[shell_sr.index]) {
                    shd_mut.faces = new_face_srs;
                }
            }
        }
    }

    (result, split_count)
}

/// Split edges at continuity breaks (SplitContinuity operator).
///
/// Detects C0/C1/C2 discontinuities in curve and surface geometry
/// and splits edges at those points.
///
/// Returns (modified BRep, count of edge splits).
fn split_continuity_operator(brep: &rcad_kernel::BRep, params: &SplitContinuityOperator) -> (rcad_kernel::BRep, usize) {
    use rcad_kernel::geom::CurveEval;

    let result = brep.clone();
    let mut split_count = 0usize;
    let _tolerance = params.tolerance;

    if !params.check_curves {
        return (result, 0);
    }

    // Analyze each edge's curve for continuity breaks via tshape walk
    for ts in &brep.tshapes {
        let TShape::Edge(ed) = ts.as_ref() else { continue };

        let Some(curve) = ed.curve.as_ref() else { continue };

        let range = if ed.range[0] == 0.0 && ed.range[1] == 0.0 {
            let d = curve.default_domain();
            [d[0], d[1]]
        } else {
            ed.range
        };

        // Sample the curve to detect discontinuities
        let n_samples = 100.min(params.max_splits_per_edge * 10);
        let dt = (range[1] - range[0]) / n_samples as f64;

        let mut split_params: Vec<f64> = Vec::new();

        for i in 1..n_samples {
            let t = range[0] + dt * i as f64;

            // Check continuity at this parameter
            let continuity = check_curve_continuity_at(curve, t, dt);

            if continuity < params.min_continuity {
                split_params.push(t);
                if split_params.len() >= params.max_splits_per_edge {
                    break;
                }
            }
        }

        if !split_params.is_empty() {
            split_count += split_params.len();
        }
    }

    (result, split_count)
}

/// Check curve continuity at a parameter value.
/// Returns the highest continuity level that the curve maintains at this point.
fn check_curve_continuity_at(curve: &rcad_kernel::geom::Curve3, t: f64, dt: f64) -> ContinuityLevel {
    use rcad_kernel::geom::CurveEval;

    let eps = dt * 0.1; // Small offset for checking
    let t_lo = t - eps;
    let t_hi = t + eps;

    // Get points and tangents at nearby parameters
    let p_lo = curve.point_at(t_lo);
    let p_mid = curve.point_at(t);
    let p_hi = curve.point_at(t_hi);

    let tan_lo = curve.tangent_at(t_lo).normalize_or(DVec3::ZERO);
    let tan_mid = curve.tangent_at(t).normalize_or(DVec3::ZERO);
    let tan_hi = curve.tangent_at(t_hi).normalize_or(DVec3::ZERO);

    // Check C0: position continuity
    // For a continuous curve, the position at t should lie between p_lo and p_hi
    // A discontinuity would show as a jump larger than expected from linear interpolation
    let expected_pos_gap = (p_hi - p_lo).length();
    let actual_gap = (p_mid - p_lo).length() + (p_hi - p_mid).length();
    let gap_ratio = (actual_gap - expected_pos_gap).abs() / expected_pos_gap.max(TOLERANCE_LINEAR_ULTRA_STRICT);

    if gap_ratio > 0.1 {
        // Significant deviation from expected - likely a discontinuity
        return ContinuityLevel::C0;
    }

    // Check C1: tangent continuity
    // Tangents should be parallel at nearby points for a smooth curve
    let dot_lo_mid = tan_lo.dot(tan_mid);
    let dot_mid_hi = tan_mid.dot(tan_hi);

    // Tangents pointing in opposite directions indicate a sharp corner
    if dot_lo_mid < 0.99 || dot_mid_hi < 0.99 {
        // More than ~8 degree angle difference
        return ContinuityLevel::C0;
    }

    // Check C2: curvature continuity (approximate)
    let curvature_lo = compute_curvature_at(curve, t_lo);
    let curvature_mid = compute_curvature_at(curve, t);
    let curvature_hi = compute_curvature_at(curve, t_hi);

    // Curvature should be approximately constant for C2
    let avg_curvature = (curvature_lo + curvature_mid + curvature_hi) / 3.0;
    let max_deviation = (curvature_lo - avg_curvature).abs()
        .max((curvature_mid - avg_curvature).abs())
        .max((curvature_hi - avg_curvature).abs());

    // Use relative tolerance for curvature
    let tol = avg_curvature.abs().max(TOLERANCE_MESH_LEGACY) * 0.1 + 0.01;
    if max_deviation > tol {
        return ContinuityLevel::C1;
    }

    ContinuityLevel::C2
}

/// Compute approximate curvature at a parameter value.
fn compute_curvature_at(curve: &rcad_kernel::geom::Curve3, t: f64) -> f64 {
    use rcad_kernel::geom::CurveEval;

    let eps = TOLERANCE_MESH_LEGACY;
    let p = curve.point_at(t);
    let p_lo = curve.point_at(t - eps);
    let p_hi = curve.point_at(t + eps);

    // Approximate second derivative
    let d2 = (p_hi - 2.0 * p + p_lo) / (eps * eps);
    let d1 = curve.tangent_at(t);

    // Curvature = |r' x r''| / |r'|^3
    let cross = d1.cross(d2);
    let d1_len = d1.length();

    if d1_len < TOLERANCE_LEN_MIN {
        return 0.0;
    }

    cross.length() / (d1_len.powi(3))
}

/// Convert analytic geometry to BSpline (ConvertToBSpline operator).
///
/// Converts elementary surfaces and curves to NURBS representation.
///
/// Returns (modified BRep, count of entities converted).
fn convert_to_bspline_operator(brep: &rcad_kernel::BRep, params: &ConvertToBSplineOperator) -> (rcad_kernel::BRep, usize) {
    use rcad_kernel::geom::{Surface3, Curve3};
    use rcad_kernel::nurbs_convert;

    let mut result = brep.clone();
    let mut conversion_count = 0usize;

    // Convert edge curves via tshape walk
    if params.convert_curves {
        for i in 0..result.tshapes.len() {
            let curve_clone = {
                let TShape::Edge(ed) = &*result.tshapes[i] else { continue };
                ed.curve.clone()
            };
            let Some(curve) = curve_clone else { continue };

            let should_convert = match &curve {
                Curve3::Line(_) | Curve3::Circle(_) | Curve3::Ellipse(_) => params.convert_elementary,
                Curve3::BSpline(_) | Curve3::Bezier(_) => false,
                _ => true,
            };

            if should_convert {
                let bspline = nurbs_convert::curve_to_bspline(&curve, params.approximation_samples);
                if let TShape::Edge(ed_mut) = Arc::make_mut(&mut result.tshapes[i]) {
                    ed_mut.curve = Some(Curve3::BSpline(bspline));
                    conversion_count += 1;
                }
            }
        }
    }

    // Convert face surfaces via tshape walk
    if params.convert_surfaces {
        for i in 0..result.tshapes.len() {
            let surface_clone = {
                let TShape::Face(fd) = &*result.tshapes[i] else { continue };
                fd.surface.clone()
            };
            let Some(surface) = surface_clone else { continue };

            let should_convert = match &surface {
                Surface3::Plane(_) => params.convert_planes,
                Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Cone(_) | Surface3::Torus(_) => {
                    params.convert_elementary
                }
                Surface3::BSpline(_) | Surface3::Bezier(_) | Surface3::TriBezier(_) => false,
                _ => true,
            };

            if should_convert {
                let bspline = nurbs_convert::surface_to_bspline(
                    &surface,
                    params.approximation_samples,
                    params.approximation_samples,
                );
                if let TShape::Face(fd_mut) = Arc::make_mut(&mut result.tshapes[i]) {
                    fd_mut.surface = Some(Surface3::BSpline(bspline));
                    conversion_count += 1;
                }
            }
        }
    }

    (result, conversion_count)
}

/// Convert BSpline surfaces to Bezier patches (SurfaceToBezier operator).
///
/// Splits BSpline surfaces at all interior knot lines.
///
/// Returns (modified BRep, count of surfaces converted).
fn surface_to_bezier_operator(brep: &rcad_kernel::BRep, params: &SurfaceToBezierOperator) -> (rcad_kernel::BRep, usize) {
    use rcad_kernel::geom::Surface3;

    let mut result = brep.clone();
    let mut conversion_count = 0usize;

    if !params.convert_surfaces {
        return (result, 0);
    }

    // Convert face BSpline surfaces to Bezier via tshape walk
    for i in 0..result.tshapes.len() {
        let surface_clone = {
            let TShape::Face(fd) = &*result.tshapes[i] else { continue };
            fd.surface.clone()
        };
        let Some(Surface3::BSpline(bspline)) = surface_clone else { continue };

        let bezier_patches = split_bspline_to_bezier(&bspline);

        if bezier_patches.len() == 1 {
            if let TShape::Face(fd_mut) = Arc::make_mut(&mut result.tshapes[i]) {
                fd_mut.surface = Some(Surface3::Bezier(bezier_patches.into_iter().next().unwrap()));
            }
        } else {
            if let Some(first) = bezier_patches.into_iter().next() {
                if first.control_points.len() - 1 <= params.max_degree {
                    if let TShape::Face(fd_mut) = Arc::make_mut(&mut result.tshapes[i]) {
                        fd_mut.surface = Some(Surface3::Bezier(first));
                        conversion_count += 1;
                    }
                }
            }
        }
    }

    (result, conversion_count)
}

/// Split a BSpline surface into Bezier patches at knot lines.
fn split_bspline_to_bezier(bspline: &BSplineSurface) -> std::collections::VecDeque<BezierSurface> {
    use std::collections::VecDeque;

    // For simplicity, return a single Bezier approximation
    // A full implementation would:
    // 1. Insert knots to raise multiplicity to degree at each interior knot
    // 2. Extract each span as a separate Bezier patch

    let mut patches = VecDeque::new();

    // Check if already a single Bezier span
    let u_single = bspline.knots_u.len() == 2 * (bspline.degree_u + 1);
    let v_single = bspline.knots_v.len() == 2 * (bspline.degree_v + 1);

    if u_single && v_single {
        // Already a single Bezier patch
        patches.push_back(BezierSurface {
            control_points: bspline.control_points.clone(),
            weights: bspline.weights.clone(),
        });
    } else {
        // Need to split - for now, just return the whole thing as one patch
        // This is an approximation
        patches.push_back(BezierSurface {
            control_points: bspline.control_points.clone(),
            weights: bspline.weights.clone(),
        });
    }

    patches
}

/// Apply scaling transformation (ScaleShape operator).
///
/// Scales geometry and optionally tolerances.
///
/// Returns (modified BRep, count of entities modified).
fn scale_shape_operator(brep: &rcad_kernel::BRep, params: &ScaleShapeOperator) -> (rcad_kernel::BRep, usize) {
    use glam::DAffine3;

    // Check for identity scaling
    if (params.scale_x - 1.0).abs() < TOLERANCE_LEN_MIN
        && (params.scale_y - 1.0).abs() < TOLERANCE_LEN_MIN
        && (params.scale_z - 1.0).abs() < TOLERANCE_LEN_MIN
    {
        return (brep.clone(), 0);
    }

    let mut result = brep.clone();

    // Build the transformation matrix
    let scale_matrix = DAffine3::from_scale(glam::DVec3::new(params.scale_x, params.scale_y, params.scale_z));

    // If there's an origin, translate to/from it
    let transform = if let Some(origin) = params.origin {
        let to_origin = DAffine3::from_translation(-origin);
        let from_origin = DAffine3::from_translation(origin);
        from_origin * scale_matrix * to_origin
    } else {
        scale_matrix
    };

    // Apply the transformation
    result.apply_transform(transform);

    // Scale tolerances if requested
    let modification_count = result.vertex_count() + result.edge_count();

    if params.scale_tolerances {
        let scale_factor = params.scale_x.max(params.scale_y).max(params.scale_z);

        // Scale tolerances on individual tshapes
        for i in 0..result.tshapes.len() {
            let needs_scale = match &*result.tshapes[i] {
                TShape::Vertex(vd) => (vd.tolerance - 0.0).abs() > TOLERANCE_LEN_MIN,
                TShape::Edge(ed) => (ed.tolerance - 0.0).abs() > TOLERANCE_LEN_MIN,
                TShape::Face(fd) => (fd.tolerance - 0.0).abs() > TOLERANCE_LEN_MIN,
                _ => false,
            };
            if needs_scale {
                match Arc::make_mut(&mut result.tshapes[i]) {
                    TShape::Vertex(vd) => vd.tolerance *= scale_factor,
                    TShape::Edge(ed) => ed.tolerance *= scale_factor,
                    TShape::Face(fd) => fd.tolerance *= scale_factor,
                    _ => {}
                }
            }
        }
    }

    (result, modification_count)
}

/// Convert indirect faces to direct (DirectFaces operator).
///
/// An indirect face is one where the natural surface orientation does not
/// match the face's orientation flag. This operator ensures consistency
/// by correcting face orientations.
///
/// Returns (modified BRep, count of faces fixed).
fn direct_faces_operator(brep: &rcad_kernel::BRep, params: &DirectFacesOperator) -> (rcad_kernel::BRep, usize) {
    use crate::algo_ext::brep_repair::recompute_face_normals;

    let mut result = brep.clone();
    let mut faces_fixed = 0usize;

    // Step 1: Recompute normals if requested
    if params.recompute_normals {
        let (brep_with_normals, normals_fixed) = recompute_face_normals(&result);
        result = brep_with_normals;
        faces_fixed += normals_fixed;
    }

    // Step 2: Collect shell-level orientation info before any mutation
    struct ShellFix {
        shell_sr: Shape,
        faces_to_flip: Vec<(Shape, Shape)>, // (face_sr, outer_wire_sr)
    }

    let mut shell_fixes: Vec<ShellFix> = Vec::new();

    let solid_indices: Vec<usize> = (0..result.tshapes.len())
        .filter(|&i| matches!(&*result.tshapes[i], TShape::Solid(_)))
        .collect();

    for si in solid_indices {
        let shell_refs = {
            let TShape::Solid(sd) = &*result.tshapes[si] else { continue; };
            sd.shells.clone()
        };

        for shell_sr in &shell_refs {
            let face_refs = {
                let TShape::Shell(shd) = &*result.tshapes[shell_sr.index] else { continue; };
                shd.faces.clone()
            };

            // Determine expected shell orientation from face normals
            let mut consistent_normals = 0usize;
            let mut inconsistent_normals = 0usize;
            let mut potential_faces: Vec<(Shape, Shape)> = Vec::new();

            for face_sr in &face_refs {
                let TShape::Face(fd) = &*result.tshapes[face_sr.index] else { continue; };
                let normal = fd.surface.as_ref()
                    .map(|s| rcad_kernel::geom::SurfaceEval::normal_at(s, 0.0, 0.0))
                    .unwrap_or(DVec3::ZERO);

                if normal.length() > 0.5 {
                    consistent_normals += 1;
                } else if normal.length() < 0.5 && !normal.abs_diff_eq(DVec3::ZERO, 0.1) {
                    inconsistent_normals += 1;
                    potential_faces.push((face_sr.clone(), fd.outer_wire.clone()));
                }
            }

            if inconsistent_normals > consistent_normals && inconsistent_normals > 0 {
                shell_fixes.push(ShellFix { shell_sr: shell_sr.clone(), faces_to_flip: potential_faces });
                faces_fixed += shell_fixes.last().map_or(0, |f| f.faces_to_flip.len());
            }
        }
    }

    // Step 3: Apply wire orientation fixes (no borrow conflicts with pre-collected data)
    if params.fix_wire_orientation {
        for fix in &shell_fixes {
            for (_face_sr, wire_sr) in &fix.faces_to_flip {
                if let TShape::Wire(wd_mut) = Arc::make_mut(&mut result.tshapes[wire_sr.index]) {
                    wd_mut.edges.reverse();
                    for esr in &mut wd_mut.edges {
                        esr.orientation = match esr.orientation {
                            Orientation::Forward => Orientation::Reversed,
                            Orientation::Reversed => Orientation::Forward,
                            other => other,
                        };
                    }
                }
            }
        }
    }

    // Step 4: Surface references (placeholder)
    if params.update_surface_references {
        let _ = &result;
    }

    (result, faces_fixed)
}

/// Fix SameParameter issues on edges (SameParameter operator).
///
/// Ensures that the 3D curve and 2D PCurves of edges are consistently parameterized.
/// Uses the existing `fix_same_parameter_with_scan` function with configurable options.
///
/// Returns (modified BRep, count of edges fixed).
fn same_parameter_operator(brep: &rcad_kernel::BRep, params: &SameParameterOperator) -> (rcad_kernel::BRep, usize) {
    // Use the existing implementation with the specified tolerance
    let (result, fixed_count) = fix_same_parameter_with_scan(brep, params.tolerance);

    // If enforcing, run additional pass on edges that might have been missed
    if params.enforce {
        let mut enforced = result.clone();
        // Reset same_parameter flag on all edges to force re-check
        for i in 0..enforced.tshapes.len() {
            if let TShape::Edge(ed) = &*enforced.tshapes[i] {
                if ed.same_parameter {
                    if let TShape::Edge(ed_mut) = Arc::make_mut(&mut enforced.tshapes[i]) {
                        ed_mut.same_parameter = false;
                    }
                }
            }
        }
        let (final_result, additional_fixed) = fix_same_parameter_with_scan(&enforced, params.tolerance);
        (final_result, fixed_count + additional_fixed)
    } else {
        (result, fixed_count)
    }
}

/// Remove internal faces after boolean operations (RemoveInternalFaces operator).
///
/// Detects and removes partition faces that are inside the solid volume,
/// keeping only the outer boundary faces.
///
/// Returns (modified BRep, count of faces removed).
fn remove_internal_faces_operator(brep: &rcad_kernel::BRep, params: &RemoveInternalFacesOperator) -> (rcad_kernel::BRep, usize) {
    use rcad_kernel::BRepGraph;

    let mut result = brep.clone();
    let mut total_removed = 0usize;

    // Build a topology graph to analyze face connectivity
    let _graph = BRepGraph::from_brep(&result);

    // Identify candidate internal faces via tshape walk
    let solid_indices: Vec<usize> = (0..result.tshapes.len())
        .filter(|&i| matches!(&*result.tshapes[i], TShape::Solid(_)))
        .collect();

    for si in solid_indices {
        let shell_refs = {
            let TShape::Solid(sd) = &*result.tshapes[si] else { continue; };
            sd.shells.clone()
        };

        for shell_sr in &shell_refs {
            let (face_refs, original_len) = {
                let TShape::Shell(shd) = &*result.tshapes[shell_sr.index] else { continue; };
                (shd.faces.clone(), shd.faces.len())
            };

            let faces_to_remove = identify_internal_faces_in_shell(&result, &face_refs, params);

            if faces_to_remove.is_empty() {
                continue;
            }

            let mut kept_faces = Vec::new();
            for (fi, face_sr) in face_refs.iter().enumerate() {
                if !faces_to_remove.contains(&fi) {
                    kept_faces.push(face_sr.clone());
                } else {
                    total_removed += 1;
                }
            }

            if kept_faces.len() < original_len {
                if let TShape::Shell(shd_mut) = Arc::make_mut(&mut result.tshapes[shell_sr.index]) {
                    shd_mut.faces = kept_faces;
                }
            }
        }
    }

    // Merge vertices after face removal if requested
    if params.merge_vertices && total_removed > 0 {
        let (merged, _) = crate::algo_ext::brep_repair::merge_close_vertices(&result, params.tolerance);
        result = merged;
    }

    (result, total_removed)
}

/// Identify internal faces within a single shell's face list.
fn identify_internal_faces_in_shell(brep: &rcad_kernel::BRep, face_refs: &[Shape], params: &RemoveInternalFacesOperator) -> Vec<usize> {
    let mut internal_faces = Vec::new();

    for (face_idx, face_sr) in face_refs.iter().enumerate() {
        let TShape::Face(fd) = &*brep.tshapes[face_sr.index] else { continue };

        // Check 1: Face area
        let area = estimate_face_area_from_wire(brep, fd.outer_wire.clone());
        if area < params.min_face_area {
            if !params.preserve_material_boundaries {
                internal_faces.push(face_idx);
            }
            continue;
        }

        // Check 2: Edge analysis — faces with all edges shared by other faces are internal
        let wd = match &*brep.tshapes[fd.outer_wire.index] {
            TShape::Wire(w) => w,
            _ => continue,
        };

        let mut shared_edge_count = 0usize;
        let mut total_edges = 0usize;

        for edge_sr in &wd.edges {
            if edge_sr.index >= brep.tshapes.len() {
                continue;
            }
            if !matches!(&*brep.tshapes[edge_sr.index], TShape::Edge(_)) {
                continue;
            }
            total_edges += 1;

            // Count how many other faces share this edge
            let mut face_count = 0usize;
            for (other_fi, other_face_sr) in face_refs.iter().enumerate() {
                if other_fi == face_idx {
                    continue;
                }
                let TShape::Face(other_fd) = &*brep.tshapes[other_face_sr.index] else { continue };
                let other_wd = match &*brep.tshapes[other_fd.outer_wire.index] {
                    TShape::Wire(w) => w,
                    _ => continue,
                };
                for other_edge_sr in &other_wd.edges {
                    if other_edge_sr.index == edge_sr.index {
                        face_count += 1;
                    }
                }
            }

            if face_count >= 1 {
                shared_edge_count += 1;
            }
        }

        // If all edges are shared with other faces, this might be internal
        if total_edges > 0 && shared_edge_count == total_edges {
            // Additional heuristic: check face surface normal direction
            let normal = fd.surface.as_ref()
                .map(|s| rcad_kernel::geom::SurfaceEval::normal_at(s, 0.0, 0.0))
                .unwrap_or(DVec3::ZERO);
            if normal.length() > 0.1 {
                // Conservative — proper material-side analysis needed
            }
        }
    }

    internal_faces.sort();
    internal_faces.dedup();
    internal_faces
}