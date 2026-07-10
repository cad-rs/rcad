
/// Comprehensive geometry healing (HealGeometry operator).
///
/// Combines multiple repair operations into a single configurable pass.
///
/// Returns (modified BRep, repair report).
fn heal_geometry_operator(brep: &rcad_kernel::BRep, params: &HealGeometryOperator) -> (rcad_kernel::BRep, RepairReport) {
    use crate::brep_repair::{
        fix_face_orientation, fix_wire_gaps, fix_uv_bounds_violations,
        recompute_face_normals, remove_degenerate_faces, propagate_tolerances,
        ToleranceFlowDirection,
    };

    let mut current = brep.clone();
    let mut total_report = RepairReport::default();

    let sequence = params.get_sequence();

    for _pass in 0..params.max_passes {
        let pass_start_totals = (
            total_report.vertices_merged,
            total_report.faces_reoriented,
            total_report.wires_fixed,
            total_report.same_parameter_fixed,
            total_report.same_range_fixed,
        );

        for step in &sequence {
            match step {
                HealGeometryStep::RecomputeNormals => {
                    let (next, fixed) = recompute_face_normals(&current);
                    current = next;
                    total_report.normals_recomputed += fixed;
                }
                HealGeometryStep::FixSameRange => {
                    let (next, fixed) = fix_same_range_with_scan(&current, params.tolerance);
                    current = next;
                    total_report.same_range_fixed += fixed;
                }
                HealGeometryStep::FixSameParameter => {
                    let (next, fixed) = fix_same_parameter_with_scan(&current, params.tolerance);
                    current = next;
                    total_report.same_parameter_fixed += fixed;
                }
                HealGeometryStep::FixFaceOrientation => {
                    let (next, fixed) = fix_face_orientation(&current);
                    current = next;
                    total_report.faces_reoriented += fixed;
                }
                HealGeometryStep::FixWireGaps => {
                    let (next, report) = fix_wire_gaps(&current, params.tolerance, params.tolerance * 10.0);
                    current = next;
                    total_report.wires_fixed += report.wires_fixed;
                }
                HealGeometryStep::FixUvBounds => {
                    let (next, report) = fix_uv_bounds_violations(&current, params.tolerance);
                    current = next;
                    total_report.faces_reoriented += report.faces_adjusted;
                }
                HealGeometryStep::RemoveDegenerateFaces => {
                    let (next, fixed) = remove_degenerate_faces(&current);
                    current = next;
                    total_report.degenerate_faces_removed += fixed;
                }
                HealGeometryStep::RemoveSmallEdges => {
                    let (next, fixed) = crate::brep_repair::remove_small_edges(&current, params.min_edge_length);
                    current = next;
                    total_report.vertices_merged += fixed;
                }
                HealGeometryStep::PropagateTolerances => {
                    current = propagate_tolerances(&current, params.tolerance, ToleranceFlowDirection::BottomUp);
                }
            }
        }

        // Check if this pass made any changes
        let pass_end_totals = (
            total_report.vertices_merged,
            total_report.faces_reoriented,
            total_report.wires_fixed,
            total_report.same_parameter_fixed,
            total_report.same_range_fixed,
        );

        if pass_end_totals == pass_start_totals {
            // No changes this pass - stop iterating
            break;
        }
    }

    (current, total_report)
}

/// Run a healing pipeline with rollback support and progress callbacks.
///
/// This is an enhanced version of `run_healing_operator_chain` that supports:
/// - Automatic rollback on failure
/// - Progress callbacks for monitoring
/// - Result aggregation
///
/// # Arguments
/// * `brep` - The BRep to process.
/// * `operators` - The sequence of operators to execute.
/// * `options` - Healing options.
/// * `rollback_config` - Configuration for rollback behavior.
/// * `progress_callback` - Optional callback for progress monitoring.
///
/// # Returns
/// A tuple of (processed BRep, PipelineExecutionReport).
pub fn run_healing_pipeline_with_rollback(
    brep: &rcad_kernel::BRep,
    operators: &[HealingOperator],
    options: HealingOptions,
    rollback_config: RollbackConfig,
    progress_callback: Option<&dyn ProgressCallback>,
) -> (rcad_kernel::BRep, PipelineExecutionReport) {
    use std::time::Instant;

    let start_time = Instant::now();
    let mut current = brep.clone();
    let mut aggregation = OperatorResultAggregation::new();
    let mut snapshots: Vec<BRepSnapshot> = Vec::new();

    // Create initial snapshot
    if rollback_config.enabled {
        snapshots.push(BRepSnapshot::new(
            brep,
            0,
            "initial",
            start_time.elapsed().as_secs_f64(),
        ));
    }

    let initial_issues = brep_check_analyze(brep).issues.len();
    let mut best_state = (brep.clone(), initial_issues, 0); // (brep, issues, operator_index)

    for (op_idx, op) in operators.iter().enumerate() {
        // Check for cancellation
        if let Some(cb) = progress_callback
            && cb.is_cancelled() {
                let final_brep = current.clone();
                let report = PipelineExecutionReport {
                    aggregation,
                    snapshots,
                    final_brep: final_brep.clone(),
                    completed: false,
                    failure_reason: Some("Cancelled by user".to_string()),
                    rollback_index: None,
                };
                return (final_brep, report);
            }

        // Notify progress callback
        if let Some(cb) = progress_callback {
            cb.on_operator_start(op_idx, op);
            let progress = (op_idx as f64) / (operators.len() as f64);
            cb.on_progress(progress, &format!("Executing operator {}/{}", op_idx + 1, operators.len()));
        }

        // Create snapshot if configured
        if rollback_config.enabled
            && (rollback_config.snapshot_before_each
                || rollback_config.snapshot_indices.contains(&op_idx))
        {
            // Limit number of snapshots
            if snapshots.len() >= rollback_config.max_snapshots {
                snapshots.remove(0);
            }
            snapshots.push(BRepSnapshot::new(
                &current,
                op_idx,
                format!("before_operator_{}", op_idx),
                start_time.elapsed().as_secs_f64(),
            ));
        }

        // Execute the operator
        let op_start = Instant::now();
        let issues_before = brep_check_analyze(&current).issues.len();

        let (next, healing_report) = run_healing_operator_chain(&current, options, std::slice::from_ref(op));
        current = next;

        let issues_after = brep_check_analyze(&current).issues.len();
        let op_elapsed = op_start.elapsed().as_secs_f64();

        // Build operator result
        let changed = issues_before != issues_after;
        let issues_fixed = issues_before.saturating_sub(issues_after);
        let modifications = healing_report.passes.iter()
            .map(|p| p.vertices_merged + p.degenerate_faces_removed + p.normals_recomputed
                + p.faces_reoriented + p.wires_fixed + p.same_range_fixed + p.same_parameter_fixed)
            .sum();

        let result = OperatorResult {
            operator: op.clone(),
            changed,
            modifications,
            issues_fixed,
            description: if changed {
                format!("Fixed {} issues", issues_fixed)
            } else {
                "No changes".to_string()
            },
            elapsed_seconds: op_elapsed,
            skipped: false,
            skip_reason: None,
        };

        // Check for rollback conditions
        let mut should_rollback = false;
        let mut rollback_reason = None;

        if rollback_config.enabled {
            // Check for regression
            if rollback_config.rollback_on_regression && issues_after > issues_before {
                should_rollback = true;
                rollback_reason = Some(format!(
                    "Issue regression: {} -> {} issues",
                    issues_before, issues_after
                ));
            }

            // Check threshold
            if rollback_config.max_issues_threshold > 0 && issues_after > rollback_config.max_issues_threshold {
                should_rollback = true;
                rollback_reason = Some(format!(
                    "Issues exceed threshold: {} > {}",
                    issues_after, rollback_config.max_issues_threshold
                ));
            }
        }

        // Track best state for potential rollback
        if issues_after < best_state.1 {
            best_state = (current.clone(), issues_after, op_idx);
        }

        // Notify progress callback
        if let Some(cb) = progress_callback {
            cb.on_operator_complete(op_idx, &result);
        }

        aggregation.add_result(result);

        // Handle rollback
        if should_rollback {
            if let Some(ref reason) = rollback_reason
                && let Some(cb) = progress_callback {
                    cb.on_error(op_idx, reason);
                }

            // Find the best snapshot to rollback to
            let rollback_idx = if issues_before <= issues_after {
                // Rollback to before this operator
                op_idx.saturating_sub(1)
            } else {
                // Keep current state but note the issue
                best_state.2
            };

            // Find snapshot for rollback
            let rollback_snapshot = snapshots.iter().rev().find(|s| s.operator_index <= rollback_idx).cloned();
            if let Some(snapshot) = rollback_snapshot {
                current = snapshot.brep.clone();
                aggregation.rollback_triggered = true;
                aggregation.rollback_reason = rollback_reason.clone();

                let final_brep = current.clone();
                let report = PipelineExecutionReport {
                    aggregation,
                    snapshots,
                    final_brep: final_brep.clone(),
                    completed: false,
                    failure_reason: rollback_reason,
                    rollback_index: Some(snapshot.operator_index),
                };
                return (final_brep, report);
            }
        }
    }

    let report = PipelineExecutionReport {
        aggregation,
        snapshots,
        final_brep: current.clone(),
        completed: true,
        failure_reason: None,
        rollback_index: None,
    };

    (current, report)
}

// ─────────────────────────────────────────────────────────────────────────────
// Advanced Operator Chain Execution
// ─────────────────────────────────────────────────────────────────────────────

/// Run an advanced operator chain with conditions and dependencies.
///
/// This provides the enhanced chaining capabilities including:
/// - Conditional execution
/// - Operator dependencies
/// - Result propagation
///
/// # Example
/// ```ignore
/// use rcad_algorithms::healing::{
///     run_advanced_operator_chain, OperatorChainConfig,
///     HealingOperatorWithCondition, OperatorCondition, HealingOperator,
/// };
///
/// let config = OperatorChainConfig {
///     operators: vec![
///         HealingOperatorWithCondition::new(HealingOperator::ParametricConsistency),
///         HealingOperatorWithCondition::with_condition(
///             HealingOperator::Repair,
///             OperatorCondition::OnlyIfIssues,
///         ),
///     ],
///     ..Default::default()
/// };
///
/// let (result, report) = run_advanced_operator_chain(&brep, &config);
/// ```
pub fn run_advanced_operator_chain(brep: &rcad_kernel::BRep, config: &OperatorChainConfig) -> (rcad_kernel::BRep, OperatorChainReport) {
    use std::time::Instant;

    let start_time = Instant::now();
    let initial = brep_check_analyze(brep);
    let initial_stats = HealingIssueStats::from_check_result(&initial);

    let mut current = brep.clone();
    let mut operator_results: Vec<OperatorResult> = Vec::new();
    let mut operators_executed = 0usize;
    let mut operators_skipped = 0usize;

    // Build options from config
    let options = HealingOptions {
        tolerance: config.base_tolerance,
        ..config.healing_options
    };

    for op_with_cond in config.operators.iter() {
        // Check dependencies
        let mut skip = false;
        let mut skip_reason = None;

        for &dep_idx in &op_with_cond.dependencies {
            if let Some(dep_result) = operator_results.get(dep_idx)
                && !dep_result.changed && op_with_cond.skip_on_dependency_failure {
                    skip = true;
                    skip_reason = Some(format!("Dependency {} made no changes", dep_idx));
                    break;
                }
        }

        // Check condition if dependencies passed
        if !skip
            && let Some(ref condition) = op_with_cond.condition {
                let (_, temp_report) = analyze_and_heal(&current, HealingOptions {
                    mode: HealingMode::AnalyzeOnly,
                    ..options
                });
                if !condition.evaluate(&current, &temp_report, &operator_results) {
                    skip = true;
                    skip_reason = Some("Condition not met".to_string());
                }
            }

        if skip {
            operator_results.push(OperatorResult {
                operator: op_with_cond.operator.clone(),
                changed: false,
                modifications: 0,
                issues_fixed: 0,
                description: String::new(),
                elapsed_seconds: 0.0,
                skipped: true,
                skip_reason,
            });
            operators_skipped += 1;
            continue;
        }

        // Execute the operator
        let op_start = Instant::now();
        let issues_before = brep_check_analyze(&current).issues.len();

        // Convert HealingOperatorWithCondition's operator to simple operator
        let simple_op = op_with_cond.operator.clone();

        // Run the operator
        let (next, _) = run_healing_operator_chain(
            &current,
            options,
            &[simple_op.clone()],
        );
        current = next;

        let issues_after = brep_check_analyze(&current).issues.len();
        let op_elapsed = op_start.elapsed().as_secs_f64();

        let changed = issues_before != issues_after;
        let issues_fixed = issues_before.saturating_sub(issues_after);

        operator_results.push(OperatorResult {
            operator: simple_op,
            changed,
            modifications: issues_fixed,
            issues_fixed,
            description: if changed {
                format!("Fixed {} issues", issues_fixed)
            } else {
                "No changes".to_string()
            },
            elapsed_seconds: op_elapsed,
            skipped: false,
            skip_reason: None,
        });
        operators_executed += 1;

        // Check stop condition
        if config.stop_on_clean && brep_check_analyze(&current).is_valid() {
            break;
        }
    }

    let final_result = brep_check_analyze(&current);
    let final_stats = HealingIssueStats::from_check_result(&final_result);
    let total_elapsed = start_time.elapsed().as_secs_f64();

    let is_clean = final_result.is_valid();
    let summary = if is_clean {
        format!("Shape is clean after {} operators ({:.3}s)", operators_executed, total_elapsed)
    } else {
        format!(
            "{} issues remain after {} operators ({:.3}s)",
            final_result.issues.len(),
            operators_executed,
            total_elapsed
        )
    };

    (
        current,
        OperatorChainReport {
            operator_results,
            initial,
            final_result,
            initial_stats,
            final_stats,
            total_elapsed_seconds: total_elapsed,
            operators_executed,
            operators_skipped,
            is_clean,
            summary,
        },
    )
}
