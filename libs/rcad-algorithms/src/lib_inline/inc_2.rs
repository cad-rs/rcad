/// Like [`boolean_op`] — standard path with tolerance correction.
///
/// Calls [`boolean_op`] internally, then runs tolerance correction
/// on the result before wrapping it as old rcad_kernel::BRep for backward compatibility.
pub fn boolean_op_with_retry(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
) -> Result<rcad_kernel::BRep, BooleanError> {
 let t = boolean_op(op, a, b)?;
 Ok(t)
}

/// Perform a boolean operation with advanced execution options and report.
pub fn boolean_op_with_options(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
 mut options: BooleanOptions,
) -> Result<(rcad_kernel::BRep, BooleanExecutionReport), BooleanError> {
 merge_pairwise_model_tol_into_boolean_options(&mut options, a, b);

 let input_faces_a = face_count_of(a);
 let input_faces_b = face_count_of(b);
 let used_bvh = options.use_bvh && has_faces(a) && has_faces(b);

 let (mut out, mut report, history_opt) = if options.include_history {
 let (result, history) = if options.use_bvh {
 if options.fuzzy_tol <= 0.0 && !options.use_glue {
 let (r, h) = boolean_op_with_history(op, a, b)?; (r, h)
 } else {
 let ds_tol = if options.fuzzy_tol > 0.0 { options.fuzzy_tol } else { TOLERANCE_ABS };
 let mut ds = bopds::ds::DS::new_from_topods(a, b, ds_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.configure_fuzzy(options.fuzzy_tol);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let (t, h) = builder.build_with_history()?;
 (t, h)
 }
 } else {
 let ds_tol = if options.fuzzy_tol > 0.0 { options.fuzzy_tol } else { TOLERANCE_ABS };
 let mut ds = bopds::ds::DS::new_from_topods(a, b, ds_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let mut filler = pave_filler::PaveFiller::new(&mut ds);
 filler.brep = Some(&mut brep);
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.configure_fuzzy(options.fuzzy_tol);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let (t, h) = builder.build_with_history()?;
 (t, h)
 };
 (
 result,
 BooleanExecutionReport {
 input_faces_a,
 input_faces_b,
 used_bvh,
 ..BooleanExecutionReport::default()
 },
 Some(history),
 )
 } else {
 let result = if options.use_bvh {
 if options.fuzzy_tol > 0.0 || options.use_glue {
 let mut ds = bopds::ds::DS::new_from_topods(a, b, options.fuzzy_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let r = builder.build()?;
 boolean_postprocess_pave_result_topods(op, a, b, r)?
 } else {
 boolean_op(op, a, b)?
 }
 } else {
 let ds_tol = if options.fuzzy_tol > 0.0 { options.fuzzy_tol } else { TOLERANCE_ABS };
 let mut ds = bopds::ds::DS::new_from_topods(a, b, ds_tol);
 let mut brep = rcad_kernel::topods::BRep::new();
 let (face_refs, ic_edge_map) = {
 let mut filler = pave_filler::PaveFiller::new(&mut ds);
 filler.brep = Some(&mut brep);
 filler.configure_glue(options.use_glue, options.glue_tolerance);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map)
 .with_glue(options.use_glue, options.glue_tolerance);
 let r = builder.build()?;
 boolean_postprocess_pave_result_topods(op, a, b, r)?
 };
 (
 result,
 BooleanExecutionReport {
 input_faces_a,
 input_faces_b,
 used_bvh,
 ..BooleanExecutionReport::default()
 },
 None,
 )
 };

 if options.run_healing {
 let mut healing_options = options.healing;
 // If boolean make-connected is enabled, allow healing to use the same
 // connectivity rebuild policy when repair passes stall.
 if options.run_make_connected {
 healing_options.make_connected_prepass_mode = MakeConnectedPrepassMode::IssueDriven;
 healing_options.run_make_connected_on_stall = true;
 healing_options.make_connected_tolerance = options.make_connected_tolerance;
 healing_options.make_connected_max_passes = options.make_connected_max_passes;
 healing_options.make_connected_tolerance_growth =
 options.make_connected_tolerance_growth;
 healing_options.make_connected_tolerance_cap = options.make_connected_tolerance_cap;
 }
 let (healed, heal_report) = analyze_and_heal(&out, healing_options);
 out = healed;
 report.healed = true;
 report.healing_report = Some(heal_report);
 }

 if options.run_make_connected {
 let old_for_mc = out.clone();
 let (connected, connected_report) = run_make_connected_for_boolean_output(
 &old_for_mc,
 history_opt.as_ref(),
 &options,
 &mut report,
 );
 out = connected;
 report.made_connected = true;
 report.make_connected_report = Some(connected_report);
 }

 if options.run_simplify {
 let (simplified, simp_report) = simplify_brep_post_ops(&out, options.simplify);
 out = simplified;
 report.simplified = true;
 report.simplify_report = Some(simp_report);
 }

 if options.run_propagate_geom_tolerances {
 let floor = resolved_boolean_fuzzy_tol_for_ds(options.fuzzy_tol);
 let old_out = out.clone();
 out = propagate_tolerances(&old_out, floor, ToleranceFlowDirection::BottomUp);
 report.propagated_geom_tolerances = true;
 }

 report.output_faces = face_count_of(&out);
 report.configured_fuzzy_tol = options.fuzzy_tol;
 report.effective_fuzzy_tol = resolved_boolean_fuzzy_tol_for_ds(options.fuzzy_tol);
 report.boolean_history = history_opt.as_ref().cloned();
 if let Some(history) = history_opt {
 report.history_faces = history.len();
 report.history_edges = history.edge_origins.len();
 report.history_vertices = history.vertex_origins.len();
 report.history_shells = history.shell_origins.len();
 report.history_solids = history.solid_origins.len();
 report.persistent_face_labels = persistent_face_labels_from_history(&history);
 report.persistent_edge_labels = persistent_edge_labels_from_history(&history);
 report.persistent_shell_labels = persistent_shell_labels_from_history(&history);
 report.persistent_solid_labels = persistent_solid_labels_from_history(&history);
 }

 Ok((out, report))
}

/// Robust boolean operation with automatic fuzzy-tolerance retries.
///
/// Attempts run in this order:
/// 1. `options.base.fuzzy_tol`
/// 2. each value in `options.fuzzy_retry_ladder`
///
/// The first successful attempt is returned, with retry metadata in
/// [`BooleanExecutionReport`].
pub fn boolean_op_robust(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
 options: BooleanRobustOptions,
) -> Result<(rcad_kernel::BRep, BooleanExecutionReport), BooleanError> {
 const MAX_RETRY_ESCALATION_ROUNDS: usize = 2;

 let mut pending = std::collections::VecDeque::new();
 pending.push_back((options.base.fuzzy_tol.max(0.0), None, 0usize));
 let mut tried: Vec<(f64, Option<BooleanRetryClass>, usize)> = Vec::new();
 let mut attempt_reports: Vec<BooleanRobustAttemptReport> = Vec::new();
 let mut last_err: Option<BooleanError> = None;

 while let Some((fuzzy, origin_retry_class, retry_round)) = pending.pop_front() {
 if tried.iter().any(|(v, cls, round)| {
 (*v - fuzzy).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP && *cls == origin_retry_class && *round == retry_round
 }) {
 continue;
 }
 tried.push((fuzzy, origin_retry_class, retry_round));

 let mut attempt_options = options.base;
 attempt_options.fuzzy_tol = fuzzy;
 tune_boolean_options_for_retry_class(&mut attempt_options, origin_retry_class, retry_round);
 let attempt_make_connected_scoped_enabled =
 attempt_options.run_make_connected && attempt_options.make_connected_scoped;
 let attempt_scope_seed_mode =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_seed_mode)
 } else {
 None
 };
 let attempt_scope_history_ring_depth =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_history_ring_depth)
 } else {
 None
 };
 let attempt_scope_seed_length =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_seed_length)
 } else {
 None
 };
 let attempt_scope_min_history_edges =
 if attempt_options.run_make_connected && attempt_options.make_connected_scoped {
 Some(attempt_options.make_connected_scope_min_history_edges)
 } else {
 None
 };
 match boolean_op_with_options(op, a, b, attempt_options) {
 Ok((brep, mut report)) => {
 attempt_reports.push(BooleanRobustAttemptReport {
 fuzzy_tol: fuzzy,
 success: true,
 retry_round,
 origin_retry_class,
 make_connected_scoped_enabled: attempt_make_connected_scoped_enabled,
 make_connected_scope_seed_mode: report.make_connected_scope_seed_mode,
 make_connected_scope_history_ring_depth: report
 .make_connected_scope_history_ring_depth,
 make_connected_scope_seed_length: attempt_scope_seed_length,
 make_connected_scope_min_history_edges: attempt_scope_min_history_edges,
 make_connected_scope_seed_source: report.make_connected_scope_seed_source,
 make_connected_scope_history_seed_edge_count: Some(
 report.make_connected_scope_history_seed_edge_count,
 ),
 make_connected_scope_heuristic_seed_edge_count: Some(
 report.make_connected_scope_heuristic_seed_edge_count,
 ),
 make_connected_scope_seed_vertex_count: Some(
 report.make_connected_scope_seed_vertices.len(),
 ),
 make_connected_scope_seed_edge_count: Some(
 report.make_connected_scope_seed_edges.len(),
 ),
 used_glue: attempt_options.use_glue,
 glue_tolerance: attempt_options.glue_tolerance,
 retry_class: None,
 error_message: None,
 output_faces: Some(report.output_faces),
 made_connected: report.made_connected,
 make_connected_scope_fallback_applied: report
 .make_connected_scope_fallback_applied,
 make_connected_scope_fallback_reason: report
 .make_connected_scope_fallback_reason,
 make_connected_scope_seed_edge_coverage: report
 .make_connected_scope_seed_edge_coverage,
 make_connected_scope_seed_face_coverage: report
 .make_connected_scope_seed_face_coverage,
 make_connected_scope_global_fallback_initial_tolerance: report
 .make_connected_scope_global_fallback_initial_tolerance,
 make_connected_scope_global_fallback_max_passes: report
 .make_connected_scope_global_fallback_max_passes,
 });
 report.robust_attempts = attempt_reports;
 report.retry_count = tried.len().saturating_sub(1);
 report.configured_fuzzy_tol = fuzzy;
 report.effective_fuzzy_tol = resolved_boolean_fuzzy_tol_for_ds(fuzzy);
 return Ok((brep, report));
 }
 Err(err) => {
 let retry_class = classify_boolean_retry(&err);
 attempt_reports.push(BooleanRobustAttemptReport {
 fuzzy_tol: fuzzy,
 success: false,
 retry_round,
 origin_retry_class,
 make_connected_scoped_enabled: attempt_make_connected_scoped_enabled,
 make_connected_scope_seed_mode: attempt_scope_seed_mode,
 make_connected_scope_history_ring_depth: attempt_scope_history_ring_depth,
 make_connected_scope_seed_length: attempt_scope_seed_length,
 make_connected_scope_min_history_edges: attempt_scope_min_history_edges,
 make_connected_scope_seed_source: None,
 make_connected_scope_history_seed_edge_count: None,
 make_connected_scope_heuristic_seed_edge_count: None,
 make_connected_scope_seed_vertex_count: None,
 make_connected_scope_seed_edge_count: None,
 used_glue: attempt_options.use_glue,
 glue_tolerance: attempt_options.glue_tolerance,
 retry_class: Some(retry_class),
 error_message: Some(format!("{err:?}")),
 output_faces: None,
 made_connected: false,
 make_connected_scope_fallback_applied: false,
 make_connected_scope_fallback_reason: None,
 make_connected_scope_seed_edge_coverage: None,
 make_connected_scope_seed_face_coverage: None,
 make_connected_scope_global_fallback_initial_tolerance: None,
 make_connected_scope_global_fallback_max_passes: None,
 });
 for candidate in boolean_retry_followup_attempts(
 fuzzy,
 &options.fuzzy_retry_ladder,
 &err,
 options.retry_policy,
 origin_retry_class,
 retry_round,
 MAX_RETRY_ESCALATION_ROUNDS,
 attempt_make_connected_scoped_enabled,
 ) {
 let seen = tried.iter().any(|(v, cls, round)| {
 (*v - candidate.0).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP
 && *cls == candidate.1
 && *round == candidate.2
 }) || pending.iter().any(|(v, cls, round)| {
 (*v - candidate.0).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP
 && *cls == candidate.1
 && *round == candidate.2
 });
 if !seen {
 pending.push_back(candidate);
 }
 }
 last_err = Some(err);
 }
 }
 }

 Err(last_err.unwrap_or(BooleanError::DegenerateResult))
}

/// Run post-operation simplification passes on a rcad_kernel::BRep.
pub fn simplify_brep_post_ops(brep: &topods::BRep, options: SimplifyOptions) -> (topods::BRep, SimplifyReport) {
 let old = brep.clone();
 let (result, report) = simplify_brep_post_ops_old(&old, options);
 (result, report)
}

/// Legacy: takes old rcad_kernel::BRep.
fn simplify_brep_post_ops_old(brep: &rcad_kernel::BRep, options: SimplifyOptions) -> (rcad_kernel::BRep, SimplifyReport) {
 fn closure_score(brep: &rcad_kernel::BRep) -> usize {
 let report = crate::brep_check::validate_solid_closure(brep);
 report
 .issues
 .iter()
 .map(|iss| match iss {
 crate::CheckIssue::SolidNotClosed {
 boundary_edge_count,
 ..
 } => *boundary_edge_count,
 _ => 1,
 })
 .sum()
 }

 let before = brep_check_analyze(brep);
 let mut out = brep.clone();
 let mut report = SimplifyReport {
 issues_before: before.issues.len(),
 ..SimplifyReport::default()
 };

 if options.merge_vertices {
 let (next, merged) = merge_close_vertices(&out, options.merge_tolerance);
 out = next;
 report.vertices_merged = merged;
 }
 if options.recompute_normals {
 let (next, n) = recompute_face_normals(&out);
 out = next;
 report.normals_recomputed = n;
 }
 if options.remove_degenerate_faces {
 let (next, n) = remove_degenerate_faces(&out);
 out = next;
 report.degenerate_faces_removed = n;
 }
 if options.remove_internal_faces {
 let (next, n) = remove_internal_faces(&out);
 out = next;
 report.internal_faces_removed = n;
 }
 if options.fix_wire_orientation {
 let (next, n) = fix_wire_orientation(&out, options.merge_tolerance);
 out = next;
 report.wires_fixed = n;
 }
 if options.unify_same_domain_faces {
 let cur_score = closure_score(&out);
 let (next, n) = unify_same_domain_faces(&out);
 let next_score = closure_score(&next);
 if next_score <= cur_score {
 out = next;
 report.same_domain_face_merges = n;
 }
 }
  // After same-domain unification, run same-domain unification once more to
 // absorb newly adjacent coplanar patches produced by the fuse pass.
 if options.unify_same_domain_faces {
 let cur_score = closure_score(&out);
 let (next, n) = unify_same_domain_faces(&out);
 let next_score = closure_score(&next);
 if next_score <= cur_score {
 out = next;
 report.same_domain_face_merges += n;
 }
 }

 // Kernel-level wire cleanup: collapse consecutive collinear segments so
 // post-boolean faces do not keep fragmented edge chains.
 let collinear_edge_merges = rcad_kernel::merge_collinear_edges_in_wires(
 &mut out,
 options.merge_tolerance.max(tolerance::TOLERANCE_ABS),
 );
 report.wires_fixed += collinear_edge_merges;

 if options.remove_small_edges {
 let cur_score = closure_score(&out);
 let (next, n) = remove_small_edges(&out, options.small_edge_min_length);
 let next_score = closure_score(&next);
 if next_score <= cur_score {
 out = next;
 report.small_edges_removed = n;
 }
 }

 // Final safety net: never return an open solid from simplification if it
 // can be repaired into a closed one with the standard solid fixer.
 if !crate::brep_check::validate_solid_closure(&out).is_clean() {
 let (fixed, _fix_report) =
 fix_solid(&out, options.merge_tolerance.max(tolerance::TOLERANCE_ABS));
 if crate::brep_check::validate_solid_closure(&fixed).is_clean() {
 out = fixed;
 } else {
 let (healed, _heal_report) = heal_comprehensive(&out, &HealingOptions::default());
 if crate::brep_check::validate_solid_closure(&healed).is_clean() {
 out = healed;
 }
 }
 }

 // Mesh check removed — topods::BRep handles triangulation via separate tessellation pipeline.

 report.issues_after = brep_check_analyze(&out).issues.len();
 (out, report)
}

/// Boolean + simplification convenience pipeline.
///
/// Mirrors OCCT's `BRepAlgoAPI::SimplifyResult()` which runs
/// `BRepLib_MakeConnected` before simplification to merge coincident
/// vertices and edges, ensuring clean topological connectivity.
pub fn boolean_op_simplified(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
 options: SimplifyOptions,
) -> Result<(rcad_kernel::BRep, SimplifyReport), BooleanError> {
 let t = boolean_op_topods_simplified(op, a, b, options)?;
 Ok((t.0, t.1))
}

/// Same as boolean_op_simplified but returns topods::BRep.
pub fn boolean_op_topods_simplified(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
 options: SimplifyOptions,
) -> Result<(topods::BRep, SimplifyReport), BooleanError> {
 let raw = boolean_op(op, a, b)?;
 let (connected, _mc_report) = make_connected_enhanced(
 &raw,
 tolerance::TOLERANCE_ABS,
 3,
 );
 let (simplified, report) = simplify_brep_post_ops(&connected, options);
 Ok((simplified, report))
}

/// Split `target` by one or more `tools` without boolean classification.
///
/// This is a first-stage splitter built on top of [`imprint_shape`]. It keeps
/// target material and iteratively imprints tool boundaries onto the evolving
/// target shape.
pub fn split_shape(target: &rcad_kernel::BRep, tools: &[rcad_kernel::BRep]) -> (rcad_kernel::BRep, SplitterReport) {
 split_shape_with_options(target, tools, SplitterOptions::default())
}

/// Like [`split_shape`] with advanced options.
pub fn split_shape_with_options(
 target: &rcad_kernel::BRep,
 tools: &[rcad_kernel::BRep],
 options: SplitterOptions,
) -> (rcad_kernel::BRep, SplitterReport) {
 let (result, report) = split_brep_internal_with_partial_report(target, tools, options, false);
 match result {
 Ok(brep) => (brep, report),
 Err(_) => unreachable!("unchecked splitter path should not fail"),
 }
}

/// Split `target` by tools and validate each executed step.
///
/// Returns a step-indexed error if an intermediate split result has structural
/// validity issues, excluding `NonManifoldEdge` (which can be expected for
/// split-first intermediate topology).
pub fn split_shape_checked_with_options(
 target: &rcad_kernel::BRep,
 tools: &[rcad_kernel::BRep],
 options: SplitterOptions,
) -> Result<(rcad_kernel::BRep, SplitterReport), SplitterError> {
 let (result, report) = split_brep_internal_with_partial_report(target, tools, options, true);
 result.map(|brep| (brep, report))
}

fn split_brep_internal_with_partial_report(
 target: &rcad_kernel::BRep,
 tools: &[rcad_kernel::BRep],
 options: SplitterOptions,
 validate_each_step: bool,
) -> (Result<rcad_kernel::BRep, SplitterError>, SplitterReport) {
 let mut acc = target.clone();
 let mut report = SplitterReport::default();

 for (step_index, tool) in tools.iter().enumerate() {
 let input_faces = face_count_of(&acc);
 let fuzzy = options.fuzzy_tolerance.max(0.0);
 let skipped_by_broad_phase =
 options.broad_phase_pruning && breps_farther_than_tolerance(&acc, tool, fuzzy);

 if skipped_by_broad_phase {
 report.steps.push(SplitterStepReport {
 step_index,
 input_faces,
 seam_edges: 0,
 output_faces: input_faces,
 healed: false,
 skipped_by_broad_phase: true,
 validation_issue_count: if validate_each_step { Some(0) } else { None },
 validation_first_issue: None,
 });
 continue;
 }

 let mut step = imprint_shape(&acc, &tool);
 let seam_edges = step.seam_edges.len();

 if options.heal_after_each_step {
 let mut healing = options.healing;
 align_healing_options_with_boolean_operands(
 &mut healing,
 &acc,
 tool,
 options.fuzzy_tolerance,
 );
 let (healed, _) = analyze_and_heal(&step.brep, healing);
 step.brep = healed;
 }

 let step_old = step.brep.clone();

 let mut validation_issue_count = None;
 let mut validation_first_issue = None;
 let output_faces = face_count_of(&step_old);
 if validate_each_step {
 let validity = brep_check_analyze(&step_old);
 let (issue_count, first_issue) =
 splitter_issues_by_level(&validity, options.validation_level);
 validation_issue_count = Some(issue_count);
 validation_first_issue = first_issue.clone();
 if issue_count > 0 {
 report.steps.push(SplitterStepReport {
 step_index,
 input_faces,
 seam_edges,
 output_faces,
 healed: options.heal_after_each_step,
 skipped_by_broad_phase: false,
 validation_issue_count,
 validation_first_issue,
 });
 return (
 Err(SplitterError::StepInvalid {
 step_index,
 issue_count,
 first_issue,
 }),
 report,
 );
 }
 }

 report.total_seam_edges += seam_edges;
 report.steps.push(SplitterStepReport {
 step_index,
 input_faces,
 seam_edges,
 output_faces,
 healed: options.heal_after_each_step,
 skipped_by_broad_phase: false,
 validation_issue_count,
 validation_first_issue,
 });

 acc = step_old;
 }

 (Ok(acc), report)
}

fn brep_bounds(brep: &rcad_kernel::BRep) -> Option<(glam::DVec3, glam::DVec3)> {
 let vi = rcad_kernel::vertex_indices(brep);
 let first = brep.vertex_point(*vi.first()?)?;
 let mut min = first;
 let mut max = first;
 for &v in &vi {
 if let Some(pt) = brep.vertex_point(v) {
 min = min.min(pt);
 max = max.max(pt);
 }
 }
 Some((min, max))
}

fn aabb_distance(
 min_a: glam::DVec3,
 max_a: glam::DVec3,
 min_b: glam::DVec3,
 max_b: glam::DVec3,
) -> f64 {
 let dx = if max_a.x < min_b.x {
 min_b.x - max_a.x
 } else if max_b.x < min_a.x {
 min_a.x - max_b.x
 } else {
 0.0
 };
 let dy = if max_a.y < min_b.y {
 min_b.y - max_a.y
 } else if max_b.y < min_a.y {
 min_a.y - max_b.y
 } else {
 0.0
 };
 let dz = if max_a.z < min_b.z {
 min_b.z - max_a.z
 } else if max_b.z < min_a.z {
 min_a.z - max_b.z
 } else {
 0.0
 };
 (dx * dx + dy * dy + dz * dz).sqrt()
}

fn breps_farther_than_tolerance(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep, tol: f64) -> bool {
 let Some((min_a, max_a)) = brep_bounds(a) else {
 return false;
 };
 let Some((min_b, max_b)) = brep_bounds(b) else {
 return false;
 };
 aabb_distance(min_a, max_a, min_b, max_b) > tol
}

fn splitter_issues_by_level(
 validity: &CheckResult,
 level: SplitterValidationLevel,
) -> (usize, Option<String>) {
 let filtered: Vec<&CheckIssue> = match level {
 SplitterValidationLevel::Relaxed => validity
 .issues
 .iter()
 .filter(|issue| !matches!(issue, CheckIssue::NonManifoldEdge { .. }))
 .collect(),
 SplitterValidationLevel::Strict => validity.issues.iter().collect(),
 };
 (filtered.len(), filtered.first().map(|it| it.to_string()))
}

/// Split each object by a shared set of tools.
///
/// This is a grouped splitter API similar to object/tool workflows in mature
/// boolean kernels: every input object is split against all tools, and results
/// are returned in object order.
pub fn split_objects_with_tools(
 objects: &[rcad_kernel::BRep],
 tools: &[rcad_kernel::BRep],
) -> (Vec<rcad_kernel::BRep>, SplitterObjectsReport) {
 split_objects_with_tools_options(objects, tools, SplitterOptions::default())
}

/// Like [`split_objects_with_tools`] but with advanced options.
pub fn split_objects_with_tools_options(
 objects: &[rcad_kernel::BRep],
 tools: &[rcad_kernel::BRep],
 options: SplitterOptions,
) -> (Vec<rcad_kernel::BRep>, SplitterObjectsReport) {
 let mut outputs = Vec::with_capacity(objects.len());
 let mut objects_report = Vec::with_capacity(objects.len());

 for (object_index, object) in objects.iter().enumerate() {
 let (split, report) = split_shape_with_options(object, tools, options);
 outputs.push(split);
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: true,
 error: None,
 });
 }

 (
 outputs,
 SplitterObjectsReport {
 objects: objects_report,
 },
 )
}

/// Checked grouped splitter variant.
///
/// Validates each split step for each object and returns the first error.
pub fn split_objects_with_tools_checked_options(
 objects: &[rcad_kernel::BRep],
 tools: &[rcad_kernel::BRep],
 options: SplitterOptions,
) -> Result<(Vec<rcad_kernel::BRep>, SplitterObjectsReport), SplitterError> {
 let mut outputs = Vec::with_capacity(objects.len());
 let mut objects_report = Vec::with_capacity(objects.len());

 for (object_index, object) in objects.iter().enumerate() {
 let (split, report) = split_shape_checked_with_options(object, tools, options)?;
 outputs.push(split);
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: true,
 error: None,
 });
 }

 Ok((
 outputs,
 SplitterObjectsReport {
 objects: objects_report,
 },
 ))
}

/// Checked grouped splitter with per-object failure collection.
///
/// Unlike [`split_objects_with_tools_checked_options`], this function does not
/// fail fast. It records per-object errors in the returned report and keeps
/// processing remaining objects.
pub fn split_objects_with_tools_checked_collect_options(
 objects: &[rcad_kernel::BRep],
 tools: &[rcad_kernel::BRep],
 options: SplitterOptions,
) -> (Vec<Option<rcad_kernel::BRep>>, SplitterObjectsReport) {
 let mut outputs = Vec::with_capacity(objects.len());
 let mut objects_report = Vec::with_capacity(objects.len());

 for (object_index, object) in objects.iter().enumerate() {
 let (result, report) =
 split_brep_internal_with_partial_report(object, tools, options, true);
 match result {
 Ok(split) => {
 outputs.push(Some(split));
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: true,
 error: None,
 });
 }
 Err(err) => {
 outputs.push(None);
 objects_report.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: false,
 error: Some(err),
 });
 }
 }
 }

 (
 outputs,
 SplitterObjectsReport {
 objects: objects_report,
 },
 )
}

/// Like [`boolean_op`] but also returns a [`BooleanHistory`] mapping each result
/// face back to its source in solid A or B.
pub fn boolean_op_with_history(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
) -> Result<(rcad_kernel::BRep, BooleanHistory), BooleanError> {
 if matches!(op, BooleanOpType::Union) {
 let (t, hist) = bop_occt_ops::fuse_with_history(a, b)?;
 return Ok((t, hist));
 }

 let mut ds = bopds::ds::DS::new_from_topods(a, b, TOLERANCE_ABS);
 let fuzzy_tol = ds.fuzzy_tol;
 let mut brep = rcad_kernel::topods::BRep::new();
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let (face_refs, ic_edge_map) = {
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.set_run_parallel(false);
 filler.configure_fuzzy(fuzzy_tol);
 filler.set_non_destructive(false);
 filler.configure_glue(false, TOLERANCE_ABS);
 filler.set_use_obb(false);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 ds.build_container_images();
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
 let (t, history) = builder.build_with_history()?;
 Ok((t, history))
}

pub fn boolean_op_par(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
) -> Result<(rcad_kernel::BRep, BooleanHistory), BooleanError> {
 if matches!(op, BooleanOpType::Union) {
 let (t, h) = bop_occt_ops::fuse_with_history_par(a, b)?;
 return Ok((t, h));
 }

 let mut ds = bopds::ds::DS::new_from_topods(a, b, TOLERANCE_ABS);
 let fuzzy_tol = ds.fuzzy_tol;
 let mut brep = rcad_kernel::topods::BRep::new();
 let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
 let (face_refs, ic_edge_map) = {
 let mut filler = match (&bvh_a, &bvh_b) {
 (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh_and_brep(&mut ds, ba, bb, &mut brep),
 _ => {
 let mut f = pave_filler::PaveFiller::new(&mut ds);
 f.brep = Some(&mut brep);
 f
 }
 };
 filler.set_run_parallel(true);
 filler.configure_fuzzy(fuzzy_tol);
 filler.set_non_destructive(false);
 filler.configure_glue(false, TOLERANCE_ABS);
 filler.set_use_obb(false);
 filler.perform();
 (std::mem::take(&mut filler.face_refs), std::mem::take(&mut filler.ic_edge_map))
 };
 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
 let (t, history) = builder.build_with_history()?;
 Ok((t, history))
}

/// Check if any solid in the rcad_kernel::BRep has at least one face (deep check across all solids).
fn has_any_face(brep: &rcad_kernel::BRep) -> bool {
 rcad_kernel::face_count(brep) > 0
}

/// Build BVHs for both BReps if they have faces; returns None for empty BReps.
fn build_optional_bvhs(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> (Option<bvh::Bvh>, Option<bvh::Bvh>) {
 let has_faces_a = rcad_kernel::face_count(a) > 0;
 let has_faces_b = rcad_kernel::face_count(b) > 0;
 (
 if has_faces_a {
 Some(bvh::Bvh::build(a))
 } else {
 None
 },
 if has_faces_b {
 Some(bvh::Bvh::build(b))
 } else {
 None
 },
 )
}

fn has_faces(brep: &rcad_kernel::BRep) -> bool {
 rcad_kernel::face_count(brep) > 0
}

fn make_connected_seed_vertices_from_short_edges(brep: &rcad_kernel::BRep, seed_length: f64) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 let threshold = seed_length.max(tolerance::TOLERANCE_ABS);
 for (ei, ts) in brep.tshapes.iter().enumerate() {
 let topods::TShape::Edge(ed) = &**ts else { continue; };
 let Some(ps) = brep.vertex_point(ed.first.index) else { continue; };
 let Some(pe) = brep.vertex_point(ed.last.index) else { continue; };
 if (pe - ps).length() <= threshold {
 out.insert(ed.first.index);
 out.insert(ed.last.index);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_near_duplicates(brep: &rcad_kernel::BRep, seed_length: f64) -> Vec<usize> {
 let vi = rcad_kernel::vertex_indices(brep);
 let mut out = std::collections::BTreeSet::new();
 let threshold = seed_length.max(tolerance::TOLERANCE_ABS);
 let threshold2 = threshold * threshold;
 for i in 0..vi.len() {
 let pi = brep.vertex_point(vi[i]).unwrap_or(glam::DVec3::ZERO);
 for j in (i + 1)..vi.len() {
 let pj = brep.vertex_point(vi[j]).unwrap_or(glam::DVec3::ZERO);
 let d2 = (pi - pj).length_squared();
 if d2 <= threshold2 {
 out.insert(vi[i]);
 out.insert(vi[j]);
 }
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_tolerance_tagged_edges(
 brep: &rcad_kernel::BRep,
 tolerance_threshold: f64,
) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 let threshold = tolerance_threshold.max(tolerance::TOLERANCE_ABS);
 for (ei, ts) in brep.tshapes.iter().enumerate() {
 let topods::TShape::Edge(ed) = &**ts else { continue; };
 let edge_tol = ed.tolerance;
 if edge_tol >= threshold {
 out.insert(ed.first.index);
 out.insert(ed.last.index);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_multi_pcurve_edges(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 for (ei, ts) in brep.tshapes.iter().enumerate() {
 let topods::TShape::Edge(ed) = &**ts else { continue; };
 if ed.pcurves.len() >= 2 {
 out.insert(ed.first.index);
 out.insert(ed.last.index);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_vertices_from_topology_seam_candidates(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let mut out = std::collections::BTreeSet::new();
 for ei in rcad_kernel::periodic_seam_edge_indices(brep) {
 if let topods::TShape::Edge(ed) = &*brep.tshapes[ei] {
 out.insert(ed.first.index);
 out.insert(ed.last.index);
 }
 }
 out.into_iter().collect()
}

fn make_connected_seed_edges_from_short_edges(brep: &rcad_kernel::BRep, seed_length: f64) -> Vec<usize> {
 let mut out = Vec::new();
 let threshold = seed_length.max(tolerance::TOLERANCE_ABS);
 for (ei, ts) in brep.tshapes.iter().enumerate() {
 let topods::TShape::Edge(ed) = &**ts else { continue; };
 let Some(ps) = brep.vertex_point(ed.first.index) else { continue; };
 let Some(pe) = brep.vertex_point(ed.last.index) else { continue; };
 if (pe - ps).length() <= threshold {
 out.push(ei);
 }
 }
 out
}

fn make_connected_seed_edges_from_near_duplicates(brep: &rcad_kernel::BRep, seed_length: f64) -> Vec<usize> {
 let dup_vertices: std::collections::HashSet<usize> =
 make_connected_seed_vertices_from_near_duplicates(brep, seed_length)
 .into_iter()
 .collect();
 let mut out = Vec::new();
 for (ei, ts) in brep.tshapes.iter().enumerate() {
 let topods::TShape::Edge(ed) = &**ts else { continue; };
 if dup_vertices.contains(&ed.first.index) || dup_vertices.contains(&ed.last.index) {
 out.push(ei);
 }
 }
 out
}

fn make_connected_seed_edges_from_tolerance_tagged_edges(
 brep: &rcad_kernel::BRep,
 tolerance_threshold: f64,
) -> Vec<usize> {
 let threshold = tolerance_threshold.max(tolerance::TOLERANCE_ABS);
 let mut out = Vec::new();
 for (ei, ts) in brep.tshapes.iter().enumerate() {
 let topods::TShape::Edge(ed) = &**ts else { continue; };
 if ed.tolerance >= threshold {
 out.push(ei);
 }
 }
 out
}

fn make_connected_seed_edges_from_multi_pcurve_edges(brep: &rcad_kernel::BRep) -> Vec<usize> {
 let mut out = Vec::new();
 for (ei, ts) in brep.tshapes.iter().enumerate() {
 let topods::TShape::Edge(ed) = &**ts else { continue; };
 if ed.pcurves.len() >= 2 {
 out.push(ei);
 }
 }
 out
}

fn make_connected_seed_edges_from_topology_seam_candidates(brep: &rcad_kernel::BRep) -> Vec<usize> {
 rcad_kernel::periodic_seam_edge_indices(brep)
}

fn make_connected_seed_edges(
 brep: &rcad_kernel::BRep,
 seed_length: f64,
 mode: MakeConnectedScopeSeedMode,
) -> Vec<usize> {
 match mode {
 MakeConnectedScopeSeedMode::ShortEdges => {
 make_connected_seed_edges_from_short_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::NearDuplicateVertices => {
 make_connected_seed_edges_from_near_duplicates(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::ToleranceTaggedEdges => {
 make_connected_seed_edges_from_tolerance_tagged_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::MultiPcurveEdges => {
 make_connected_seed_edges_from_multi_pcurve_edges(brep)
 }
 MakeConnectedScopeSeedMode::TopologySeamCandidates => {
 make_connected_seed_edges_from_topology_seam_candidates(brep)
 }
 MakeConnectedScopeSeedMode::Hybrid => {
 let mut set = std::collections::BTreeSet::new();
 for ei in make_connected_seed_edges_from_short_edges(brep, seed_length) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_near_duplicates(brep, seed_length) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_tolerance_tagged_edges(brep, seed_length) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_multi_pcurve_edges(brep) {
 set.insert(ei);
 }
 for ei in make_connected_seed_edges_from_topology_seam_candidates(brep) {
 set.insert(ei);
 }
 set.into_iter().collect()
 }
 }
}

fn make_connected_seed_vertices_from_edge_ids(brep: &rcad_kernel::BRep, edge_ids: &[usize]) -> Vec<usize> {
 let mut set = std::collections::BTreeSet::new();
 for &ei in edge_ids {
 if let topods::TShape::Edge(ed) = &*brep.tshapes[ei] {
 set.insert(ed.first.index);
 set.insert(ed.last.index);
 }
 }
 set.into_iter().collect()
}

fn select_scoped_seed_edges(
 brep: &rcad_kernel::BRep,
 history: Option<&BooleanHistory>,
 seed_length: f64,
 mode: MakeConnectedScopeSeedMode,
 history_ring_depth: usize,
 min_history_edges: usize,
) -> (Vec<usize>, usize, usize, MakeConnectedScopeSeedSource) {
 let history_seed_edges_raw = history
 .map(|h| make_connected_seed_edges_from_boolean_history(brep, h))
 .unwrap_or_default();
 // Expand history-derived seeds to configurable ring depth around boolean
 // interface topology while preserving raw-history count semantics for reports.
 let history_seed_edges =
 expand_seed_edges_with_ring_depth(brep, &history_seed_edges_raw, history_ring_depth);
 let heuristic_seed_edges = make_connected_seed_edges(brep, seed_length, mode);

 if history_seed_edges_raw.is_empty() {
 return (
 heuristic_seed_edges.clone(),
 0,
 heuristic_seed_edges.len(),
 MakeConnectedScopeSeedSource::Heuristic,
 );
 }

 if history_seed_edges_raw.len() < min_history_edges {
 let mut set = std::collections::BTreeSet::new();
 for ei in &history_seed_edges {
 set.insert(*ei);
 }
 for ei in &heuristic_seed_edges {
 set.insert(*ei);
 }
 return (
 set.into_iter().collect(),
 history_seed_edges_raw.len(),
 heuristic_seed_edges.len(),
 MakeConnectedScopeSeedSource::HistoryAugmentedHeuristic,
 );
 }

 (
 history_seed_edges.clone(),
 history_seed_edges_raw.len(),
 heuristic_seed_edges.len(),
 MakeConnectedScopeSeedSource::History,
 )
}

fn expand_seed_edges_with_ring_depth(
 brep: &rcad_kernel::BRep,
 seed_edges: &[usize],
 ring_depth: usize,
) -> Vec<usize> {
 let mut out: std::collections::BTreeSet<usize> = seed_edges.iter().copied().collect();
 if ring_depth == 0 || seed_edges.is_empty() {
 return out.into_iter().collect();
 }

 let mut visited_faces = std::collections::BTreeSet::new();
 let mut frontier = std::collections::BTreeSet::new();
 for &ei in seed_edges {
 for fi in rcad_kernel::edge_adjacent_faces(brep, ei) {
 if visited_faces.insert(fi) {
 frontier.insert(fi);
 }
 }
 }

 for _ in 0..ring_depth {
 if frontier.is_empty() {
 break;
 }
 let current: Vec<usize> = frontier.iter().copied().collect();
 frontier.clear();

 for fi in current {
 for fei in rcad_kernel::face_edges(brep, fi) {
 out.insert(fei);
 for nfi in rcad_kernel::edge_adjacent_faces(brep, fei) {
 if visited_faces.insert(nfi) {
 frontier.insert(nfi);
 }
 }
 }
 }
 }

 out.into_iter().collect()
}

fn make_connected_seed_edges_from_boolean_history(
 brep: &rcad_kernel::BRep,
 history: &BooleanHistory,
) -> Vec<usize> {
 let mut seed_edges = std::collections::BTreeSet::new();
 let edge_count = brep.edge_count();

 // If edge history is available, prefer boundary-like generated/split edges.
 for (ei, origin) in history.edge_origins.iter().enumerate() {
 if ei >= edge_count {
 break;
 }
 if matches!(
 origin,
 EdgeOrigin::Generated | EdgeOrigin::SplitFromA(_) | EdgeOrigin::SplitFromB(_)
 ) {
 seed_edges.insert(ei);
 }
 }

 // Fallback semantic extraction from face history: edges adjacent to both A and B faces
 // are strong candidates for boolean interface cleanup.
 for ei in 0..edge_count {
 let adjacent = rcad_kernel::edge_adjacent_faces(brep, ei);
 if adjacent.is_empty() {
 continue;
 }
 let mut has_a = false;
 let mut has_b = false;
 let mut has_generated = false;
 for fi in adjacent {
 if fi >= history.face_origins.len() {
 continue;
 }
 match history.face_origins[fi] {
 FaceOrigin::FromA(_) => has_a = true,
 FaceOrigin::FromB(_) => has_b = true,
 FaceOrigin::Generated => has_generated = true,
 }
 }
 if has_generated || (has_a && has_b) {
 seed_edges.insert(ei);
 }
 }

 seed_edges.into_iter().collect()
}

fn make_connected_seed_edge_labels(brep: &rcad_kernel::BRep, edge_ids: &[usize]) -> Vec<String> {
 edge_ids
 .iter()
 .map(|&ei| {
 let topods::TShape::Edge(ed) = &*brep.tshapes[ei] else { return format!("edge.{ei}.invalid-edge") };
 let pa = brep.vertex_point(ed.first.index);
 let pb = brep.vertex_point(ed.last.index);
 match (pa, pb) {
 (Some(a), Some(b)) => {
 let a_label = format!("{:.9},{:.9},{:.9}", a.x, a.y, a.z);
 let b_label = format!("{:.9},{:.9},{:.9}", b.x, b.y, b.z);
 if a_label <= b_label {
 format!("edge.{ei}.{a_label}->{b_label}")
 } else {
 format!("edge.{ei}.{b_label}->{a_label}")
 }
 }
 _ => format!("edge.{ei}.invalid-vertex"),
 }
 })
 .collect()
}

fn make_connected_seed_vertices(
 brep: &rcad_kernel::BRep,
 seed_length: f64,
 mode: MakeConnectedScopeSeedMode,
) -> Vec<usize> {
 match mode {
 MakeConnectedScopeSeedMode::ShortEdges => {
 make_connected_seed_vertices_from_short_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::NearDuplicateVertices => {
 make_connected_seed_vertices_from_near_duplicates(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::ToleranceTaggedEdges => {
 make_connected_seed_vertices_from_tolerance_tagged_edges(brep, seed_length)
 }
 MakeConnectedScopeSeedMode::MultiPcurveEdges => {
 make_connected_seed_vertices_from_multi_pcurve_edges(brep)
 }
 MakeConnectedScopeSeedMode::TopologySeamCandidates => {
 make_connected_seed_vertices_from_topology_seam_candidates(brep)
 }
 MakeConnectedScopeSeedMode::Hybrid => {
 let mut set = std::collections::BTreeSet::new();
 for v in make_connected_seed_vertices_from_short_edges(brep, seed_length) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_near_duplicates(brep, seed_length) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_tolerance_tagged_edges(brep, seed_length) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_multi_pcurve_edges(brep) {
 set.insert(v);
 }
 for v in make_connected_seed_vertices_from_topology_seam_candidates(brep) {
 set.insert(v);
 }
 set.into_iter().collect()
 }
 }
}

/// Create stable per-face labels from boolean history.
pub fn persistent_face_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .face_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 FaceOrigin::FromA(src) => format!("face.{idx}.A.{src}"),
 FaceOrigin::FromB(src) => format!("face.{idx}.B.{src}"),
 FaceOrigin::Generated => format!("face.{idx}.G"),
 })
 .collect()
}

/// Create stable per-edge labels from boolean history.
pub fn persistent_edge_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .edge_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 EdgeOrigin::FromA(src) => format!("edge.{idx}.A.{src}"),
 EdgeOrigin::FromB(src) => format!("edge.{idx}.B.{src}"),
 EdgeOrigin::Generated => format!("edge.{idx}.G"),
 EdgeOrigin::SplitFromA(src) => format!("edge.{idx}.A.split.{src}"),
 EdgeOrigin::SplitFromB(src) => format!("edge.{idx}.B.split.{src}"),
 })
 .collect()
}

/// Create stable per-shell labels from boolean history.
pub fn persistent_shell_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .shell_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 ShellOrigin::FromA => format!("shell.{idx}.A"),
 ShellOrigin::FromB => format!("shell.{idx}.B"),
 ShellOrigin::Generated => format!("shell.{idx}.G"),
 ShellOrigin::Mixed => format!("shell.{idx}.M"),
 })
 .collect()
}

/// Create stable per-solid labels from boolean history.
pub fn persistent_solid_labels_from_history(history: &BooleanHistory) -> Vec<String> {
 history
 .solid_origins
 .iter()
 .enumerate()
 .map(|(idx, origin)| match origin {
 SolidOrigin::FromA => format!("solid.{idx}.A"),
 SolidOrigin::FromB => format!("solid.{idx}.B"),
 SolidOrigin::Generated => format!("solid.{idx}.G"),
 SolidOrigin::Mixed => format!("solid.{idx}.M"),
 })
 .collect()
}

/// Union two BReps and return both the result and face origin history.
pub fn union_with_history(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<(rcad_kernel::BRep, BooleanHistory), BooleanError> {
 boolean_op_with_history(BooleanOpType::Union, a, b)
}

/// Intersect two BReps and return both the result and face origin history.
pub fn intersection_with_history(
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
) -> Result<(rcad_kernel::BRep, BooleanHistory), BooleanError> {
 boolean_op_with_history(BooleanOpType::Intersection, a, b)
}

/// Subtract solid B from solid A and return both the result and face origin history.
pub fn difference_with_history(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<(rcad_kernel::BRep, BooleanHistory), BooleanError> {
 boolean_op_with_history(BooleanOpType::Difference, a, b)
}

/// Run boolean operation followed by structured healing using default options.
pub fn boolean_op_healed(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
) -> Result<(rcad_kernel::BRep, HealingReport), BooleanError> {
 let raw = boolean_op(op, a, b)?;
 let mut healing = HealingOptions::default();
 align_healing_options_with_boolean_operands(&mut healing, a, b, 0.0);
 let (healed, report) = analyze_and_heal(&raw, healing);
 Ok((healed, report))
}

/// Run boolean operation followed by structured healing using custom options.
pub fn boolean_op_healed_with_options(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
 mut options: HealingOptions,
) -> Result<(rcad_kernel::BRep, HealingReport), BooleanError> {
 let raw = boolean_op(op, a, b)?;
 align_healing_options_with_boolean_operands(&mut options, a, b, 0.0);
 let (healed, report) = analyze_and_heal(&raw, options);
 Ok((healed, report))
}

/// Multi-body boolean fuse (union) over a list of solids.
///
/// Delegates to [`general_fuse_with_options`] with [`BooleanOptions::default`]. Each fold step
/// uses [`boolean_op_with_options`], so pairwise [`merge_pairwise_model_tol_into_boolean_options`]
/// runs on every `(accumulator, part)` pair.
pub fn general_fuse(parts: &[rcad_kernel::BRep]) -> Result<rcad_kernel::BRep, BooleanError> {
 general_fuse_with_options(parts, BooleanOptions::default())
}
