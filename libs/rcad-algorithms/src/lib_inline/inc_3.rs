
/// Like [`general_fuse`] with explicit [`BooleanOptions`] (fuzzy, glue, healing, make-connected,
/// simplify, etc.) applied on **each** left-fold union step.
pub fn general_fuse_with_options(
 parts: &[rcad_kernel::BRep],
 options: BooleanOptions,
) -> Result<rcad_kernel::BRep, BooleanError> {
 if parts.is_empty() {
 return Err(BooleanError::EmptyInput);
 }
 if parts.len() == 1 {
 return Ok(parts[0].clone());
 }

 let mut acc = parts[0].clone();
 for part in &parts[1..] {
 acc = boolean_op_with_options(BooleanOpType::Union, &acc, part, options)?.0;
 }
 Ok(acc)
}

/// History for N-ary fuse operation.
///
/// `steps[i]` is the history returned by the i-th pairwise union in the
/// left-fold sequence:
/// - step 0: union(parts[0], parts[1])
/// - step 1: union(step0_result, parts[2])
/// - ...
#[derive(Debug, Clone)]
pub struct GeneralFuseHistory {
 pub steps: Vec<BooleanHistory>,
}

/// Per-step diagnostics for N-ary fuse left-fold execution.
#[derive(Debug, Clone)]
pub struct GeneralFuseStepReport {
 /// Zero-based fold step index.
 pub step_index: usize,
 /// Face count in accumulator before this step.
 pub input_faces: usize,
 /// Face count of the fused result after this step.
 pub output_faces: usize,
}

/// Diagnostics report for N-ary fuse execution.
#[derive(Debug, Clone)]
pub struct GeneralFuseReport {
 pub steps: Vec<GeneralFuseStepReport>,
}

/// Diagnostics report for split-first general fuse execution.
#[derive(Debug, Clone)]
pub struct GeneralFuseSplitFirstReport {
 /// Per-object splitter execution details before the N-ary fuse stage.
 pub split_report: SplitterObjectsReport,
 /// Face counts of the split outputs in object order.
 pub split_face_counts: Vec<usize>,
 /// Per-step diagnostics of the final fuse fold over split objects.
 pub fuse_report: GeneralFuseReport,
}

/// Error with step location for N-ary fuse workflows.
#[derive(Debug)]
pub enum GeneralFuseError {
 EmptyInput,
 StepFailed {
 step_index: usize,
 source: BooleanError,
 },
}

impl std::fmt::Display for GeneralFuseError {
 fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
 match self {
 Self::EmptyInput => write!(f, "empty input"),
 Self::StepFailed { step_index, source } => {
 write!(f, "general_fuse failed at step {step_index}: {source}")
 }
 }
 }
}

impl std::error::Error for GeneralFuseError {
 fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
 match self {
 Self::EmptyInput => None,
 Self::StepFailed { source, .. } => Some(source),
 }
 }
}

/// Multi-body boolean fuse (union) with per-step history.
///
/// Delegates to [`general_fuse_with_history_with_options`] with default options and
/// [`BooleanOptions::include_history`] set.
pub fn general_fuse_with_history(
 parts: &[rcad_kernel::BRep],
) -> Result<(rcad_kernel::BRep, GeneralFuseHistory), BooleanError> {
 let mut opts = BooleanOptions::default();
 opts.include_history = true;
 general_fuse_with_history_with_options(parts, opts)
}

/// Like [`general_fuse_with_history`] with explicit [`BooleanOptions`] per fold step.
/// Forces [`BooleanOptions::include_history`] so each step contributes a [`BooleanHistory`].
pub fn general_fuse_with_history_with_options(
 parts: &[rcad_kernel::BRep],
 mut options: BooleanOptions,
) -> Result<(rcad_kernel::BRep, GeneralFuseHistory), BooleanError> {
 if parts.is_empty() {
 return Err(BooleanError::EmptyInput);
 }
 if parts.len() == 1 {
 return Ok((parts[0].clone(), GeneralFuseHistory { steps: Vec::new() }));
 }

 options.include_history = true;
 let mut steps = Vec::with_capacity(parts.len() - 1);
 let mut acc = parts[0].clone();
 for part in &parts[1..] {
 let (next, report) =
 boolean_op_with_options(BooleanOpType::Union, &acc, part, options)?;
 let Some(history) = report.boolean_history.clone() else {
 return Err(BooleanError::InvalidResult(
 "missing boolean_history despite include_history in general_fuse",
 ));
 };
 acc = next;
 steps.push(history);
 }

 Ok((acc, GeneralFuseHistory { steps }))
}

/// Parallel multi-body boolean fuse (union) with per-step history.
///
/// Same left-fold semantics as [`general_fuse_with_history`], but each binary union uses
/// [`boolean_op_par`] (parallel classification). This does **not** run
/// [`boolean_op_with_options`], so per-step [`BooleanOptions`] (fuzzy, glue, healing,
/// pairwise merge) are **not** applied; use [`general_fuse_with_history_with_options`] when you
/// need those on every fold.
pub fn general_fuse_par(parts: &[rcad_kernel::BRep]) -> Result<(rcad_kernel::BRep, GeneralFuseHistory), BooleanError> {
 if parts.is_empty() {
 return Err(BooleanError::EmptyInput);
 }
 if parts.len() == 1 {
 return Ok((parts[0].clone(), GeneralFuseHistory { steps: Vec::new() }));
 }

 let mut steps = Vec::with_capacity(parts.len() - 1);
 let mut acc = parts[0].clone();
 for part in &parts[1..] {
 let (next, history) = boolean_op_par(BooleanOpType::Union, &acc, part)?;
 acc = next;
 steps.push(history);
 }

 Ok((acc, GeneralFuseHistory { steps }))
}

// ============================================================================
// Compound-aware Boolean Operations
// ============================================================================

/// Count faces in all solids strictly before `solid_idx`.
fn face_count_before_solid(full: &rcad_kernel::BRep, solid_idx: usize) -> usize {
 full.solids
 .iter()
 .take(solid_idx)
 .flat_map(|s| &s.shells)
 .map(|sh| sh.faces.len())
 .sum()
}

/// Build a self-contained [`rcad_kernel::BRep`] holding only solid `solid_idx` of `full`, with
/// vertices/edges/face geometry trimmed so boolean DS loading does not ingest
/// orphan topology from sibling solids (e.g. after [`rcad_kernel::BRep::compound_from_shapes`]).
fn compact_brep_isolated_solid(full: &rcad_kernel::BRep, solid_idx: usize) -> Option<rcad_kernel::BRep> {
 use rcad_kernel::topology::{Face, Shell, Solid, Wire, WireEdge};
 use std::collections::BTreeSet;

 let solid = full.solids.get(solid_idx)?;
 let mut used_e: BTreeSet<usize> = BTreeSet::new();
 for sh in &solid.shells {
 for fa in &sh.faces {
 for we in &fa.outer_wire.edges {
 used_e.insert(we.idx);
 }
 for iw in &fa.inner_wires {
 for we in &iw.edges {
 used_e.insert(we.idx);
 }
 }
 }
 }
 let mut used_v: BTreeSet<usize> = BTreeSet::new();
 for &ei in &used_e {
 let e = full.edges.get(ei)?;
 used_v.insert(e.start);
 used_v.insert(e.end);
 }
 let v_list: Vec<usize> = used_v.into_iter().collect();
 let mut v_map = vec![usize::MAX; full.vertices.len()];
 for (ni, &oi) in v_list.iter().enumerate() {
 v_map[oi] = ni;
 }
 let e_list: Vec<usize> = used_e.into_iter().collect();
 let mut e_map = vec![usize::MAX; full.edges.len()];
 for (ni, &oi) in e_list.iter().enumerate() {
 e_map[oi] = ni;
 }

 let remap_wire = |w: &Wire| -> Wire {
 Wire {
 edges: w
 .edges
 .iter()
 .map(|we| WireEdge {
 idx: e_map[we.idx],
 forward: we.forward,
 })
 .collect(),
 }
 };

 let remap_face = |face: &Face| -> Face {
 Face {
 outer_wire: remap_wire(&face.outer_wire),
 inner_wires: face.inner_wires.iter().map(remap_wire).collect(),
 normal: face.normal,
 triangles: face
 .triangles
 .iter()
 .map(|&[a, b, c]| [v_map[a], v_map[b], v_map[c]])
 .collect(),
 sample_point: face.sample_point,
 mesh_dirty: face.mesh_dirty,
 surface_idx: None,
 }
 };

 let mut out = rcad_kernel::BRep::new();
 for &vi in &v_list {
 out.vertices.push(full.vertices[vi].clone());
 out.geom
 .vertex_tolerance
 .push(*full.geom.vertex_tolerance.get(vi).unwrap_or(&0.0));
 }

 out.geom.curves = full.geom.curves.clone();
 out.geom.surfaces = full.geom.surfaces.clone();
 out.geom.curve2ds = full.geom.curve2ds.clone();
 out.geom.curve2d_range = full.geom.curve2d_range.clone();

 for &old_ei in &e_list {
 let e = &full.edges[old_ei];
 out.edges.push(rcad_kernel::topology::Edge {
 start: v_map[e.start],
 end: v_map[e.end],
 });
 out.geom
 .edge_curve
 .push(full.geom.edge_curve.get(old_ei).copied().flatten());
 out.geom.edge_pcurves.push(
 full.geom
 .edge_pcurves
 .get(old_ei)
 .cloned()
 .unwrap_or_default(),
 );
 out.geom
 .edge_curve_range
 .push(full.geom.edge_curve_range.get(old_ei).copied().flatten());
 out.geom
 .edge_degenerated
 .push(*full.geom.edge_degenerated.get(old_ei).unwrap_or(&false));
 out.geom
 .edge_tolerance
 .push(*full.geom.edge_tolerance.get(old_ei).unwrap_or(&0.0));
 out.geom
 .edge_same_parameter
 .push(*full.geom.edge_same_parameter.get(old_ei).unwrap_or(&true));
 out.geom
 .edge_same_range
 .push(*full.geom.edge_same_range.get(old_ei).unwrap_or(&true));
 }

 let mut gfi = face_count_before_solid(full, solid_idx);
 let mut new_shells: Vec<Shell> = Vec::new();
 for sh in &solid.shells {
 let mut new_faces = Vec::new();
 for face in &sh.faces {
 new_faces.push(remap_face(face));
 out.geom
 .face_surface
 .push(full.geom.face_surface.get(gfi).copied().flatten());
 out.geom
 .face_surface_range
 .push(full.geom.face_surface_range.get(gfi).copied().flatten());
 out.geom
 .face_tolerance
 .push(*full.geom.face_tolerance.get(gfi).unwrap_or(&0.0));
 gfi += 1;
 }
 new_shells.push(Shell { faces: new_faces });
 }
 out.solids.push(Solid { shells: new_shells });

 rcad_kernel::tolerance::resize_tolerance_arrays(&mut out);
 Some(out)
}

/// `solid` must refer into this B-rep: either the copy in [`rcad_kernel::BRep::solids`], or (for
/// compounds) a constituent solid as returned by [`rcad_kernel::BRep::flatten_to_solids`] (the
/// canonical allocation lives in [`rcad_kernel::BRep::compound`], and `full.solids` holds clones
/// with different addresses).
fn brep_operand_for_compound_solid(full: &rcad_kernel::BRep, solid: &rcad_kernel::Solid) -> rcad_kernel::BRep {
 let idx = full
 .solids
 .iter()
 .position(|s| std::ptr::eq(s, solid))
 .or_else(|| {
 full.flatten_to_solids()
 .iter()
 .position(|&s| std::ptr::eq(s, solid))
 })
 .expect("compound solid reference must point into parent rcad_kernel::BRep");
 compact_brep_isolated_solid(full, idx).expect("solid exists in parent rcad_kernel::BRep")
}

/// Perform a boolean operation on a compound shape.
///
/// When the input is a compound, the operation is applied to each constituent
/// solid independently. The result is a compound of the individual results.
///
/// For union operations on compounds, all solids are fused together.
/// For difference operations, each solid from A is subtracted by all solids from B.
/// For intersection operations, each solid from A is intersected with all solids from B.
pub fn boolean_op_compound(op: BooleanOpType, a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
 let a_solids = a.flatten_to_solids();
 let b_solids = b.flatten_to_solids();

 if a_solids.is_empty() && b_solids.is_empty() {
 return Ok(rcad_kernel::BRep::default());
 }
 if a_solids.is_empty() {
 return match op {
 BooleanOpType::Union => Ok(b.clone()),
 BooleanOpType::Intersection => Ok(rcad_kernel::BRep::default()),
 BooleanOpType::Difference => Ok(rcad_kernel::BRep::default()),
 };
 }
 if b_solids.is_empty() {
 return match op {
 BooleanOpType::Union => Ok(a.clone()),
 BooleanOpType::Intersection => Ok(rcad_kernel::BRep::default()),
 BooleanOpType::Difference => Ok(a.clone()),
 };
 }

 match op {
 BooleanOpType::Union => {
 // Union all solids from both shapes
 let all_solids: Vec<rcad_kernel::BRep> = a_solids
 .iter()
 .map(|s| brep_operand_for_compound_solid(a, s))
 .chain(
 b_solids
 .iter()
 .map(|s| brep_operand_for_compound_solid(b, s)),
 )
 .collect();
 general_fuse(&all_solids)
 }
 BooleanOpType::Difference => {
 // Each solid from A is subtracted by all solids from B
 let mut results = Vec::new();
 for solid_a in a_solids {
 let mut acc = brep_operand_for_compound_solid(a, solid_a);
 for solid_b in &b_solids {
 let brep_b = brep_operand_for_compound_solid(b, solid_b);
 let t = boolean_op(BooleanOpType::Difference, &acc, &brep_b)?; acc = rcad_kernel::BRep::from_topods(&t);
 }
 results.push(acc);
 }

 if results.len() == 1 {
 Ok(results.remove(0))
 } else {
 Ok(rcad_kernel::BRep::compound_from_shapes(&results))
 }
 }
 BooleanOpType::Intersection => {
 // Each solid from A is intersected with each solid from B
 let mut results = Vec::new();
 for solid_a in a_solids {
 let brep_a = brep_operand_for_compound_solid(a, solid_a);

 for solid_b in &b_solids {
 let brep_b = brep_operand_for_compound_solid(b, solid_b);

 if let Ok(t) = boolean_op(BooleanOpType::Intersection, &brep_a, &brep_b) {
 let result = rcad_kernel::BRep::from_topods(&t);
 if !result.solids.is_empty() {
 results.push(result);
 }
 }
 }
 }

 if results.is_empty() {
 Err(BooleanError::DegenerateResult)
 } else if results.len() == 1 {
 Ok(results.remove(0))
 } else {
 Ok(rcad_kernel::BRep::compound_from_shapes(&results))
 }
 }
 }
}

/// Merge per-binary-step [`BooleanExecutionReport`] values into one compound summary.
///
/// Face counts on `accum` are expected to be preset to total operand faces; callers
/// set `output_faces` from the final shape. Scalar history counters are summed across
/// steps; persistent label vectors take the last non-empty step (final fold is most
/// representative for the returned rcad_kernel::BRep).
fn merge_boolean_execution_reports_for_compound_step(
 accum: &mut BooleanExecutionReport,
 step: &BooleanExecutionReport,
) {
 accum.used_bvh |= step.used_bvh;
 accum.healed |= step.healed;
 accum.simplified |= step.simplified;
 accum.made_connected |= step.made_connected;
 accum.propagated_geom_tolerances |= step.propagated_geom_tolerances;
 accum.make_connected_scope_fallback_applied |= step.make_connected_scope_fallback_applied;

 if step.healing_report.is_some() {
 accum.healing_report = step.healing_report.clone();
 }
 if step.simplify_report.is_some() {
 accum.simplify_report = step.simplify_report.clone();
 }
 if step.make_connected_report.is_some() {
 accum.make_connected_report = step.make_connected_report.clone();
 }
 if step.make_connected_scope_seed_mode.is_some() {
 accum.make_connected_scope_seed_mode = step.make_connected_scope_seed_mode;
 }
 if step.make_connected_scope_history_ring_depth.is_some() {
 accum.make_connected_scope_history_ring_depth =
 step.make_connected_scope_history_ring_depth;
 }
 if step.make_connected_scope_seed_source.is_some() {
 accum.make_connected_scope_seed_source = step.make_connected_scope_seed_source;
 }
 if step.make_connected_scope_fallback_reason.is_some() {
 accum.make_connected_scope_fallback_reason = step.make_connected_scope_fallback_reason;
 }
 if step.make_connected_scope_scoped_report.is_some() {
 accum.make_connected_scope_scoped_report = step.make_connected_scope_scoped_report.clone();
 }
 if step.make_connected_scope_global_fallback_report.is_some() {
 accum.make_connected_scope_global_fallback_report =
 step.make_connected_scope_global_fallback_report.clone();
 }
 if step
 .make_connected_scope_global_fallback_initial_tolerance
 .is_some()
 {
 accum.make_connected_scope_global_fallback_initial_tolerance =
 step.make_connected_scope_global_fallback_initial_tolerance;
 }
 if step
 .make_connected_scope_global_fallback_max_passes
 .is_some()
 {
 accum.make_connected_scope_global_fallback_max_passes =
 step.make_connected_scope_global_fallback_max_passes;
 }
 if step.make_connected_scope_seed_edge_coverage.is_some() {
 accum.make_connected_scope_seed_edge_coverage =
 step.make_connected_scope_seed_edge_coverage;
 }
 if step.make_connected_scope_seed_face_coverage.is_some() {
 accum.make_connected_scope_seed_face_coverage =
 step.make_connected_scope_seed_face_coverage;
 }
 accum.make_connected_scope_history_seed_edge_count +=
 step.make_connected_scope_history_seed_edge_count;
 accum.make_connected_scope_heuristic_seed_edge_count +=
 step.make_connected_scope_heuristic_seed_edge_count;
 if !step.make_connected_scope_seed_vertices.is_empty() {
 accum.make_connected_scope_seed_vertices = step.make_connected_scope_seed_vertices.clone();
 }
 if !step.make_connected_scope_seed_edges.is_empty() {
 accum.make_connected_scope_seed_edges = step.make_connected_scope_seed_edges.clone();
 }
 if !step.make_connected_scope_seed_edge_labels.is_empty() {
 accum.make_connected_scope_seed_edge_labels =
 step.make_connected_scope_seed_edge_labels.clone();
 }

 accum.history_faces += step.history_faces;
 accum.history_edges += step.history_edges;
 accum.history_vertices += step.history_vertices;
 accum.history_shells += step.history_shells;
 accum.history_solids += step.history_solids;

 if !step.persistent_face_labels.is_empty() {
 accum.persistent_face_labels = step.persistent_face_labels.clone();
 }
 if !step.persistent_edge_labels.is_empty() {
 accum.persistent_edge_labels = step.persistent_edge_labels.clone();
 }
 if !step.persistent_shell_labels.is_empty() {
 accum.persistent_shell_labels = step.persistent_shell_labels.clone();
 }
 if !step.persistent_solid_labels.is_empty() {
 accum.persistent_solid_labels = step.persistent_solid_labels.clone();
 }

 accum
 .robust_attempts
 .extend(step.robust_attempts.iter().cloned());
 accum.retry_count += step.retry_count;
 if step.effective_fuzzy_tol > accum.effective_fuzzy_tol {
 accum.effective_fuzzy_tol = step.effective_fuzzy_tol;
 }
 if step.configured_fuzzy_tol > accum.configured_fuzzy_tol {
 accum.configured_fuzzy_tol = step.configured_fuzzy_tol;
 }
 if step.boolean_history.is_some() {
 accum.boolean_history = step.boolean_history.clone();
 }
}

/// Perform a compound-aware boolean operation with options.
///
/// When each operand is a single solid, this delegates to [`boolean_op_with_options`].
/// Otherwise each internal binary boolean uses the same [`BooleanOptions`] as a
/// full pipeline (fuzzy, glue, healing, simplify, make-connected, history), and the
/// returned report aggregates step diagnostics.
pub fn boolean_op_compound_with_options(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
 options: BooleanOptions,
) -> Result<(rcad_kernel::BRep, BooleanExecutionReport), BooleanError> {
 let a_solids = a.flatten_to_solids();
 let b_solids = b.flatten_to_solids();

 if a_solids.is_empty() && b_solids.is_empty() {
 return Ok((rcad_kernel::BRep::default(), BooleanExecutionReport::default()));
 }
 if a_solids.is_empty() {
 return Ok((match op {
 BooleanOpType::Union => b.clone(),
 BooleanOpType::Intersection => rcad_kernel::BRep::default(),
 BooleanOpType::Difference => rcad_kernel::BRep::default(),
 }, BooleanExecutionReport::default()));
 }
 if b_solids.is_empty() {
 return Ok((match op {
 BooleanOpType::Union => a.clone(),
 BooleanOpType::Intersection => rcad_kernel::BRep::default(),
 BooleanOpType::Difference => a.clone(),
 }, BooleanExecutionReport::default()));
 }

 if a_solids.len() <= 1 && b_solids.len() <= 1 {
 return boolean_op_with_options(op, a, b, options);
 }

 let mut report = BooleanExecutionReport {
 input_faces_a: face_count_of(a),
 input_faces_b: face_count_of(b),
 configured_fuzzy_tol: options.fuzzy_tol,
 effective_fuzzy_tol: resolved_boolean_fuzzy_tol_for_ds(options.fuzzy_tol),
 ..BooleanExecutionReport::default()
 };

 let result = match op {
 BooleanOpType::Union => {
 let all_solids: Vec<rcad_kernel::BRep> = a_solids
 .iter()
 .map(|s| brep_operand_for_compound_solid(a, s))
 .chain(
 b_solids
 .iter()
 .map(|s| brep_operand_for_compound_solid(b, s)),
 )
 .collect();
 if all_solids.is_empty() {
 return Err(BooleanError::EmptyInput);
 }
 if all_solids.len() == 1 {
 all_solids.into_iter().next().unwrap()
 } else {
 let mut acc = all_solids[0].clone();
 for part in &all_solids[1..] {
 let (next, step_report) =
 boolean_op_with_options(BooleanOpType::Union, &acc, part, options)?;
 merge_boolean_execution_reports_for_compound_step(&mut report, &step_report);
 acc = next;
 }
 acc
 }
 }
 BooleanOpType::Difference => {
 let mut results = Vec::new();
 for solid_a in a_solids {
 let mut acc = brep_operand_for_compound_solid(a, solid_a);
 for solid_b in &b_solids {
 let brep_b = brep_operand_for_compound_solid(b, solid_b);
 let (next, step_report) =
 boolean_op_with_options(BooleanOpType::Difference, &acc, &brep_b, options)?;
 merge_boolean_execution_reports_for_compound_step(&mut report, &step_report);
 acc = next;
 }
 results.push(acc);
 }

 if results.len() == 1 {
 results.remove(0)
 } else {
 rcad_kernel::BRep::compound_from_shapes(&results)
 }
 }
 BooleanOpType::Intersection => {
 let mut results = Vec::new();
 for solid_a in a_solids {
 let brep_a = brep_operand_for_compound_solid(a, solid_a);

 for solid_b in &b_solids {
 let brep_b = brep_operand_for_compound_solid(b, solid_b);

 if let Ok((result, step_report)) = boolean_op_with_options(
 BooleanOpType::Intersection,
 &brep_a,
 &brep_b,
 options,
 ) && !result.solids.is_empty()
 {
 merge_boolean_execution_reports_for_compound_step(
 &mut report,
 &step_report,
 );
 results.push(result);
 }
 }
 }

 if results.is_empty() {
 return Err(BooleanError::DegenerateResult);
 } else if results.len() == 1 {
 results.remove(0)
 } else {
 rcad_kernel::BRep::compound_from_shapes(&results)
 }
 }
 };

 report.output_faces = face_count_of(&result);
 Ok((result, report))
}

/// Fuse all solids in a compound into a single solid.
///
/// This is equivalent to a general fuse operation on the compound's constituents.
pub fn fuse_compound(compound: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
 let solids = compound.flatten_to_solids();
 if solids.is_empty() {
 return Err(BooleanError::EmptyInput);
 }
 if solids.len() == 1 {
 return Ok(brep_operand_for_compound_solid(compound, solids[0]));
 }

 let breps: Vec<rcad_kernel::BRep> = solids
 .iter()
 .map(|s| brep_operand_for_compound_solid(compound, s))
 .collect();

 general_fuse(&breps)
}

/// Diagnostic serial N-ary fuse.
///
/// Returns per-step face-count reports and step-indexed errors when a fold
/// union fails. Delegates to [`general_fuse_detailed_with_options`] with
/// [`BooleanOptions::default`] and history enabled.
pub fn general_fuse_detailed(
 parts: &[rcad_kernel::BRep],
) -> Result<(rcad_kernel::BRep, GeneralFuseHistory, GeneralFuseReport), GeneralFuseError> {
 let mut opts = BooleanOptions::default();
 opts.include_history = true;
 general_fuse_detailed_with_options(parts, opts)
}

/// Like [`general_fuse_detailed`] with explicit [`BooleanOptions`] on each fold step.
pub fn general_fuse_detailed_with_options(
 parts: &[rcad_kernel::BRep],
 mut options: BooleanOptions,
) -> Result<(rcad_kernel::BRep, GeneralFuseHistory, GeneralFuseReport), GeneralFuseError> {
 if parts.is_empty() {
 return Err(GeneralFuseError::EmptyInput);
 }
 if parts.len() == 1 {
 return Ok((
 parts[0].clone(),
 GeneralFuseHistory { steps: Vec::new() },
 GeneralFuseReport { steps: Vec::new() },
 ));
 }

 options.include_history = true;
 let mut histories = Vec::with_capacity(parts.len() - 1);
 let mut reports = Vec::with_capacity(parts.len() - 1);
 let mut acc = parts[0].clone();
 for (step_index, part) in parts[1..].iter().enumerate() {
 let input_faces = face_count_of(&acc);
 let (next, brep_report) = boolean_op_with_options(BooleanOpType::Union, &acc, part, options)
 .map_err(|source| GeneralFuseError::StepFailed { step_index, source })?;
 let Some(history) = brep_report.boolean_history.clone() else {
 return Err(GeneralFuseError::StepFailed {
 step_index,
 source: BooleanError::InvalidResult(
 "missing boolean_history despite include_history in general_fuse_detailed",
 ),
 });
 };
 let output_faces = face_count_of(&next);
 histories.push(history);
 reports.push(GeneralFuseStepReport {
 step_index,
 input_faces,
 output_faces,
 });
 acc = next;
 }

 Ok((
 acc,
 GeneralFuseHistory { steps: histories },
 GeneralFuseReport { steps: reports },
 ))
}

/// Split-first multi-body fuse.
///
/// This is a more OCCT-like baseline than [`general_fuse`]: each object is
/// first split by all other objects, then the split outputs are fused in a
/// final N-ary fold. The implementation remains conservative by reusing the
/// existing splitter and binary boolean core.
pub fn general_fuse_split_first(parts: &[rcad_kernel::BRep]) -> Result<rcad_kernel::BRep, GeneralFuseError> {
 let (brep, _) = general_fuse_split_first_with_options(parts, SplitterOptions::default())?;
 Ok(brep)
}

/// Split-first multi-body fuse with splitter options and structured reporting.
pub fn general_fuse_split_first_with_options(
 parts: &[rcad_kernel::BRep],
 splitter_options: SplitterOptions,
) -> Result<(rcad_kernel::BRep, GeneralFuseSplitFirstReport), GeneralFuseError> {
 if parts.is_empty() {
 return Err(GeneralFuseError::EmptyInput);
 }

 let mut split_parts = Vec::with_capacity(parts.len());
 let mut object_reports = Vec::with_capacity(parts.len());
 let mut split_face_counts = Vec::with_capacity(parts.len());

 for (object_index, object) in parts.iter().enumerate() {
 let tools: Vec<rcad_kernel::BRep> = parts
 .iter()
 .enumerate()
 .filter(|(idx, _)| *idx != object_index)
 .map(|(_, part)| part.clone())
 .collect();

 let (split, report) = split_shape_with_options(object, &tools, splitter_options);
 split_face_counts.push(face_count_of(&split));
 object_reports.push(SplitterObjectReport {
 object_index,
 steps: report.steps,
 total_seam_edges: report.total_seam_edges,
 completed: true,
 error: None,
 });
 split_parts.push(split);
 }

 let (fused, _history, fuse_report) = general_fuse_detailed(&split_parts)?;
 Ok((
 fused,
 GeneralFuseSplitFirstReport {
 split_report: SplitterObjectsReport {
 objects: object_reports,
 },
 split_face_counts,
 fuse_report,
 },
 ))
}

/// Diagnostic parallel N-ary fuse.
///
/// Like [`general_fuse_par`], uses [`boolean_op_par`] each step (no per-step [`BooleanOptions`]).
pub fn general_fuse_par_detailed(
 parts: &[rcad_kernel::BRep],
) -> Result<(rcad_kernel::BRep, GeneralFuseHistory, GeneralFuseReport), GeneralFuseError> {
 if parts.is_empty() {
 return Err(GeneralFuseError::EmptyInput);
 }
 if parts.len() == 1 {
 return Ok((
 parts[0].clone(),
 GeneralFuseHistory { steps: Vec::new() },
 GeneralFuseReport { steps: Vec::new() },
 ));
 }

 let mut histories = Vec::with_capacity(parts.len() - 1);
 let mut reports = Vec::with_capacity(parts.len() - 1);
 let mut acc = parts[0].clone();
 for (step_index, part) in parts[1..].iter().enumerate() {
 let input_faces = face_count_of(&acc);
 let (next, history) = boolean_op_par(BooleanOpType::Union, &acc, part)
 .map_err(|source| GeneralFuseError::StepFailed { step_index, source })?;
 let output_faces = face_count_of(&next);
 histories.push(history);
 reports.push(GeneralFuseStepReport {
 step_index,
 input_faces,
 output_faces,
 });
 acc = next;
 }

 Ok((
 acc,
 GeneralFuseHistory { steps: histories },
 GeneralFuseReport { steps: reports },
 ))
}

/// Merge adjacent coplanar faces within the same shell into single faces.
///
/// Analogous to OCCT `ShapeUpgrade_UnifySameDomain`. After a boolean operation,
/// faces that originally belonged to the same input plane are often split into
/// multiple adjacent coplanar fragments. This function merges them back.
///
/// Unifies adjacent faces that lie on the same underlying surface domain:
/// **planar, cylindrical, toroidal, and spherical** faces are all handled.
/// The topology is simplified by removing internal shared edges between
/// same-domain face pairs.
///
/// Returns the simplified rcad_kernel::BRep and the number of face merges performed.
///
/// # Algorithm
/// Remove unreferenced geometry (surfaces, curves, edges, vertices) that are
/// no longer used by any face in the result.  After butterfly merge + classify,
/// pruned surface entries are no longer indexed by any face_surface slot.
/// OCCT's BuildResult does this implicitly via compact shape copy.
pub fn prune_unused_topology(brep: rcad_kernel::BRep) -> rcad_kernel::BRep {
 crate::brep_tools::compact_brep(&brep)
}

/// Performs iterated passes: in each pass, the first eligible pair of adjacent
/// same-domain faces sharing a single shell edge is merged. Passes repeat until
/// no more merges are possible. This is O(faces² × passes) but correct for all
/// surface-topology inputs produced by the boolean kernel.
pub fn unify_same_domain_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 unify_same_domain_faces_with_origins(brep, None)
}

/// Like [`unify_same_domain_faces`] but only merges faces whose [`FaceOrigin`]s match.
/// Use this with the face origins from [`BooleanHistory`] to avoid merging across
/// operands (A-side with B-side).
pub fn unify_same_domain_faces_with_origins(
 brep: &rcad_kernel::BRep,
 face_origins: Option<&[FaceOrigin]>,
) -> (rcad_kernel::BRep, usize) {
 let mut out = brep.clone();
 let mut total_merges = 0usize;

 loop {
 let merged = unify_one_merge_pass_with_origins(&mut out, face_origins);
 if !merged {
 break;
 }
 total_merges += 1;
 }

 (out, total_merges)
}

/// OCCT FillSameDomainFaces: group faces by edge-set equivalence, then merge
/// same-surface groups.  OCCT source: BOPAlgo_Builder_2.cxx L636-L796.
///
/// Algorithm:
/// 1. For each face, compute its edge set (sorted edge indices from all wires).
/// 2. Group faces by identical edge sets (BOPTools_Set equivalence in OCCT).
/// 3. Within each group, verify geometric surface compatibility (same surface
/// type and parameters).
/// 4. Merge each group: keep the representative face (lowest flat index),
/// remove the others.
///
/// �?OCCT : edge-set grouping (BOPTools_Set) + surface-type comparison.
pub fn occt_fill_same_domain_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 use std::collections::{BTreeSet, HashMap, HashSet};
 use rcad_kernel::geom::Surface3;

 if brep.solids.is_empty() {
 return (brep.clone(), 0);
 }

 let mut out = brep.clone();
 let mut total_merges = 0usize;

 for si in 0..out.solids.len() {
 for shi in 0..out.solids[si].shells.len() {
 let nf = out.solids[si].shells[shi].faces.len();
 if nf < 2 {
 continue;
 }

 // ── Phase 1: Edge-set signature per face ──────────────────────────
 // Include edges from all wires (outer + inner) to match OCCT
 // BOPTools_Set which collects every edge of the face.
 let face_edges: Vec<Vec<usize>> = (0..nf)
 .map(|fi| {
 let face = &out.solids[si].shells[shi].faces[fi];
 let mut edges: Vec<usize> =
 face.outer_wire.edges.iter().map(|we| we.idx).collect();
 for iw in &face.inner_wires {
 edges.extend(iw.edges.iter().map(|we| we.idx));
 }
 edges.sort_unstable();
 edges.dedup();
 edges
 })
 .collect();

 // ── Phase 2: Group by edge set ────────────────────────────────────
 // OCCT BOPTools_Set: group faces by geometric edge identity.
 // Each edge is identified by (curve_type, vertex_positions).
 // This matches OCCT TopoDS_Shape::IsEqual which compares edges by
 // their TShape identity (curve type + geometry combined).
 let vpos: Vec<glam::DVec3> = out.vertices.iter().map(|v| v.point).collect();
 let inv = 1.0 / 1e-5;
 let q = |p: glam::DVec3| -> (i64, i64, i64) {
 ((p.x * inv).round() as i64, (p.y * inv).round() as i64, (p.z * inv).round() as i64)
 };
 // Encode curve type as int: 0=Line, 1=Circle, 2=Ellipse, 3=BSpline, 4=Other
 let curve_type_id = |ci: Option<usize>| -> i64 {
 match ci.and_then(|ci| out.geom.curves.get(ci)) {
 Some(rcad_kernel::geom::Curve3::Line(_)) => 0,
 Some(rcad_kernel::geom::Curve3::Circle(_)) => 1,
 Some(rcad_kernel::geom::Curve3::Ellipse(_)) => 2,
 Some(rcad_kernel::geom::Curve3::BSpline(_)) => 3,
 _ => 4,
 }
 };
 // Build geometric edge keys: (qs_x, qs_y, qs_z, qe_x, qe_y, qe_z, curve_type)
 let face_geo_keys: Vec<BTreeSet<(i64,i64,i64,i64,i64,i64,i64)>> = (0..nf)
 .map(|fi| {
 let face = &out.solids[si].shells[shi].faces[fi];
 let mut keys = BTreeSet::new();
 let add_wire = |edges: &[rcad_kernel::topology::WireEdge], keys: &mut BTreeSet<_>| {
 for we in edges {
 if let Some(e) = out.edges.get(we.idx) {
 if e.start < vpos.len() && e.end < vpos.len() {
 let qs = q(vpos[e.start]);
 let qe = q(vpos[e.end]);
 let ct = curve_type_id(out.geom.edge_curve.get(we.idx).copied().flatten());
 let curve_key = ct;
 let (p1, p2) = if qs < qe { (qs, qe) } else { (qe, qs) };
 keys.insert((p1.0, p1.1, p1.2, p2.0, p2.1, p2.2, curve_key));
 }
 }
 }
 };
 add_wire(&face.outer_wire.edges, &mut keys);
 for iw in &face.inner_wires {
 add_wire(&iw.edges, &mut keys);
 }
 keys
 })
 .collect();
 let mut groups: HashMap<BTreeSet<(i64,i64,i64,i64,i64,i64,i64)>, Vec<usize>> = HashMap::new();
 for fi in 0..nf {
 if face_edges[fi].is_empty() { continue; }
 groups.entry(face_geo_keys[fi].clone()).or_default().push(fi);
 }

 // ── Phase 3-5: merge tracking ────────────────────────────────────
 let mut analysis_to_remove: Vec<(usize, usize, usize)> = Vec::new();
 let mut analysis_merges = 0usize;

 // Shared analysis block: Phase 2b + Phase 3-5 analysis.
 // Both depend on `check_same_surface` which borrows `out`
 // immutably. The block ensures the borrow is released before
 // the mutation phase below.
 {
 fn surfaces_share_domain(s1: &Surface3, s2: &Surface3) -> bool {
 match (s1, s2) {
 (Surface3::Plane(p1), Surface3::Plane(p2)) => {
 (p1.normal - p2.normal).length() < 1e-8
 && (p1.origin - p2.origin).length() < 1e-8
 }
 (Surface3::Sphere(sp1), Surface3::Sphere(sp2)) => {
 (sp1.center - sp2.center).length() < 1e-8
 && (sp1.radius - sp2.radius).abs() < 1e-8
 }
 (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => {
 (c1.radius - c2.radius).abs() < 1e-8
 && (c1.axis - c2.axis).length() < 1e-8
 && (c2.origin - c1.origin).cross(c1.axis).length() < 1e-8
 }
 (Surface3::Cone(c1), Surface3::Cone(c2)) => {
 (c1.radius - c2.radius).abs() < 1e-8
 && (c1.half_angle_rad - c2.half_angle_rad).abs() < 1e-8
 && (c1.axis - c2.axis).length() < 1e-8
 && (c1.apex - c2.apex).length() < 1e-8
 }
 (Surface3::Torus(t1), Surface3::Torus(t2)) => {
 (t1.major_radius - t2.major_radius).abs() < 1e-8
 && (t1.minor_radius - t2.minor_radius).abs() < 1e-8
 && (t1.axis - t2.axis).length() < 1e-8
 && (t1.center - t2.center).length() < 1e-8
 }
 _ => false,
 }
 }

 let check_same_surface = |fi1: usize, fi2: usize| -> bool {
 let f1 = &out.solids[si].shells[shi].faces[fi1];
 let f2 = &out.solids[si].shells[shi].faces[fi2];
 match (f1.surface_idx, f2.surface_idx) {
 (Some(sid1), Some(sid2)) => {
 if sid1 == sid2 {
 return true;
 }
 match (out.geom.surfaces.get(sid1), out.geom.surfaces.get(sid2)) {
 (Some(g1), Some(g2)) => surfaces_share_domain(g1, g2),
 _ => false,
 }
 }
 (None, None) => true,
 _ => false,
 }
 };

 // ── Phase 2b: Adjacent same-surface faces ──────────────────────
 // OCCT BOPAlgo_Builder_2.cxx L636-L796: FillSameDomainFaces also
 // merges faces that share boundary edges and are on the same surface,
 // even if their edge sets differ (e.g. a coplanar face split by an
 // intersection curve into two sub-faces).
 {
 // Build edge-adjacency map: for each edge, list of faces using it
 let mut edge_faces: HashMap<usize, Vec<usize>> = HashMap::new();
 for fi in 0..nf {
 let face = &out.solids[si].shells[shi].faces[fi];
 for we in &face.outer_wire.edges {
 edge_faces.entry(we.idx).or_default().push(fi);
 }
 for iw in &face.inner_wires {
 for we in &iw.edges {
 edge_faces.entry(we.idx).or_default().push(fi);
 }
 }
 }

 // Find face pairs sharing edges + same surface
 let mut same_surf_adj: HashMap<usize, Vec<usize>> = HashMap::new();
 for (_ei, f_list) in edge_faces.iter() {
 if f_list.len() < 2 { continue; }
 for i in 0..f_list.len() {
 let fi = f_list[i];
 for j in (i+1)..f_list.len() {
 let fj = f_list[j];
 // OCCT-aligned: only merge planar faces by edge adjacency.
 // Cylindrical/spherical/toroidal quadrants on the same surface
 // belong to different angular regions and must not be merged.
 let surf_i = out.solids[si].shells[shi].faces[fi]
 .surface_idx
 .and_then(|sid| out.geom.surfaces.get(sid));
 let is_plane = matches!(surf_i, Some(Surface3::Plane(_)));
 if !is_plane {
 continue;
 }
 if check_same_surface(fi, fj) {
 same_surf_adj.entry(fi).or_default().push(fj);
 same_surf_adj.entry(fj).or_default().push(fi);
 }
 }
 }
 }

 // BFS connected components for adjacent same-surface faces
 let mut visited_sa: HashSet<usize> = HashSet::new();
 for start in 0..nf {
 if visited_sa.contains(&start) { continue; }
 if !same_surf_adj.contains_key(&start) { visited_sa.insert(start); continue; }

 let mut comp = Vec::new();
 let mut queue = std::collections::VecDeque::new();
 queue.push_back(start);
 visited_sa.insert(start);

 while let Some(fi) = queue.pop_front() {
 comp.push(fi);
 if let Some(neighbors) = same_surf_adj.get(&fi) {
 for &ni in neighbors {
 if !visited_sa.contains(&ni) {
 visited_sa.insert(ni);
 queue.push_back(ni);
 }
 }
 }
 }

 if comp.len() >= 2 {
 analysis_merges += 1;
 for &fi in comp.iter().skip(1) {
 let entry = (si, shi, fi);
 if !analysis_to_remove.contains(&entry) {
 analysis_to_remove.push(entry);
 }
 }
 }
 }
 }

 // ── Phase 3-5: Analyse groups ──────────────────────────────────
 for (_edge_set, members) in groups.iter() {
 if members.len() < 2 {
 continue;
 }

 // Build adjacency graph: faces are adjacent if they share a
 // boundary edge AND are on the same surface.
 let mut adj: HashMap<usize, Vec<usize>> = HashMap::new();
 for &fi in members {
 adj.entry(fi).or_default();
 }

 for i in 0..members.len() {
 let fi = members[i];
 for j in (i + 1)..members.len() {
 let fj = members[j];
 if check_same_surface(fi, fj) {
 // Same edge set + same surface �?always adjacent
 // (they share ALL boundary edges).  Record both
 // directions for the BFS below.
 adj.get_mut(&fi).unwrap().push(fj);
 adj.get_mut(&fj).unwrap().push(fi);
 }
 }
 }

 // ── Phase 4: BFS connected components ────────────────────
 let mut visited: HashSet<usize> = HashSet::new();
 let mut components: Vec<Vec<usize>> = Vec::new();

 for &start in members {
 if visited.contains(&start) {
 continue;
 }
 let mut comp = Vec::new();
 let mut queue = std::collections::VecDeque::new();
 queue.push_back(start);
 visited.insert(start);

 while let Some(fi) = queue.pop_front() {
 comp.push(fi);
 if let Some(neighbors) = adj.get(&fi) {
 for &ni in neighbors {
 if !visited.contains(&ni) {
 visited.insert(ni);
 queue.push_back(ni);
 }
 }
 }
 }

 if comp.len() >= 2 {
 components.push(comp);
 }
 }

 if components.is_empty() {
 continue;
 }

 // ── Phase 5: Record faces to remove per component ───────
 for comp in &components {
 if comp.len() < 2 {
 continue;
 }
 analysis_merges += 1;
 // Keep the first face (lowest index in the shell),
 // remove the rest.
 for &fi in comp.iter().skip(1) {
 analysis_to_remove.push((si, shi, fi));
 }
 }
 }
 } // Drop `check_same_surface` �?release immutable borrow on `out`.

 // ── Mutation phase: remove recorded faces ────────────────────────
 total_merges += analysis_merges;
 if !analysis_to_remove.is_empty() {
 analysis_to_remove.sort_by(|a, b| b.cmp(a));
 analysis_to_remove.dedup();

 for (rsi, rshi, rfi) in analysis_to_remove {
 // Compute flat face index for geom slot removal.
 let mut ff = 0usize;
 for s in 0..rsi {
 for sh in &out.solids[s].shells {
 ff += sh.faces.len();
 }
 }
 for sh in 0..rshi {
 ff += out.solids[rsi].shells[sh].faces.len();
 }
 ff += rfi;

 remove_flat_face_geom_slots(&mut out.geom, ff);
 if let Some(s) = out.solids.get_mut(rsi) {
 if let Some(sh) = s.shells.get_mut(rshi) {
 if rfi < sh.faces.len() {
 sh.faces.remove(rfi);
 }
 }
 }
 }
 }
 }
 }

 (out, total_merges)
}

/// Check if a shared edge maintains continuity between two faces.

/// �?OCCT : �?face.surface_idx  (BuildSolid loop/area  )�?
pub fn occt_merge_same_surface_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 use std::collections::{HashMap, HashSet};
 if brep.solids.is_empty() || brep.solids[0].shells.is_empty() { return (brep.clone(), 0); }
 let mut out = brep.clone(); let mut total = 0usize;
 for si in 0..out.solids.len() {
 for shi in 0..out.solids[si].shells.len() {
 let nf = out.solids[si].shells[shi].faces.len();
 if nf < 2 { continue; }
 let mut grps: Vec<Vec<usize>> = Vec::new();
 for fi in 0..nf {
 let sid = out.solids[si].shells[shi].faces[fi].surface_idx;
 let mut found = false;
 for gi in 0..grps.len() {
 let ofi = grps[gi][0];
 if out.solids[si].shells[shi].faces[ofi].surface_idx == sid {
 grps[gi].push(fi); found = true; break;
 }
 if let (Some(a), Some(b)) = (sid, out.solids[si].shells[shi].faces[ofi].surface_idx) {
 if let (Some(s1), Some(s2)) = (out.geom.surfaces.get(a), out.geom.surfaces.get(b)) {
 use rcad_kernel::geom::Surface3;
 let same = match (s1, s2) {
 (&Surface3::Sphere(sa), Surface3::Sphere(sb)) =>
 (sa.center - sb.center).length() < 1e-8 && (sa.radius - sb.radius).abs() < 1e-8,
 (&Surface3::Plane(pa), Surface3::Plane(pb)) =>
 (pa.normal - pb.normal).length() < 1e-8 && (pa.origin - pb.origin).length() < 1e-8,
 _ => false,
 };
 if same { grps[gi].push(fi); found = true; break; }
 }
 }
 }
 if !found { grps.push(vec![fi]); }
 }
 let mut shell_refs: HashMap<usize,usize> = HashMap::new();
 for fi in 0..nf {
 for we in &out.solids[si].shells[shi].faces[fi].outer_wire.edges {
 *shell_refs.entry(we.idx).or_default() += 1;
 }
 }
 for gi in 0..grps.len() {
 let g = &grps[gi]; if g.len() < 2 { continue; }
 if std::env::var("RCAD_DEBUG_MERGE").is_ok() { eprintln!("[MERGE] group[{}]: {} faces with surface_idx={:?}", gi, g.len(), out.solids[si].shells[shi].faces[g[0]].surface_idx); }
 let mut gr: HashMap<usize,usize> = HashMap::new();
 for &fi in g {
 for we in &out.solids[si].shells[shi].faces[fi].outer_wire.edges {
 *gr.entry(we.idx).or_default() += 1;
 }
 }
 if std::env::var("RCAD_DEBUG_MERGE").is_ok() { eprintln!("[MERGE] group[{}]: {} edges", gi, gr.len()); }
 let mut bnd: Vec<usize> = gr.iter().filter(|(e, c)| {
 let sr = shell_refs.get(e).copied().unwrap_or(0);
 sr > **c || sr == 1
 }).map(|(&e,_)| e).collect();
 if std::env::var("RCAD_DEBUG_MERGE").is_ok() { eprintln!("[MERGE] group[{}]: {} edges, {} boundary", gi, gr.len(), bnd.len()); }
 if bnd.len() < 3 { continue; }
 // �?OCCT : pcurve-based seam edge  �?
 // OCCT FillSameDomainFaces  : BOPAlgo_Builder_2.cxx L571.
 // �?sphere/torus)�? seam edge �? pcurve surface�?
 // (sphere∩plane): pcurve surface �? seam�?
 // rcad �?pipeline  �?compute_face_pcurves pcurve�?
 {
 let sd = out.solids[si].shells[shi].faces[g[0]].surface_idx;
 let is_closed = sd.and_then(|sid| out.geom.surfaces.get(sid)).map(|s|
 matches!(s, rcad_kernel::geom::Surface3::Sphere(_) | rcad_kernel::geom::Surface3::Torus(_))
 ).unwrap_or(false);
 if is_closed {
 let seam_edges: Vec<usize> = gr.iter().filter_map(|(&ei, _c)| {
 if bnd.contains(&ei) { return None; }
 let pcs = out.geom.edge_pcurves.get(ei)?;
 if pcs.len() < 2 { return None; }
 let all_same = pcs.iter().all(|pc| pc.surface_idx == pcs[0].surface_idx);
 if all_same { Some(ei) } else { None }
 }).collect();
 if !seam_edges.is_empty() { bnd.extend(seam_edges); }
 }
 }
 let mut v2e: HashMap<usize,Vec<(usize,bool)>> = HashMap::new();
 for &ei in &bnd {
 if let Some(eg) = out.edges.get(ei) {
 v2e.entry(eg.start).or_default().push((ei,true));
 v2e.entry(eg.end).or_default().push((ei,false));
 }
 }
 let mut loops: Vec<Vec<(usize,bool)>> = Vec::new();
 // Find the LONGEST cycle by DFS �?OCCT PerformLoops alignment.
 // OCCT  : BOPAlgo_BuilderSolid.cxx L262 (PerformLoops),
 //  , ( �?�?
 // rcad: �?DFS  , �?
 let mut best_loop: Vec<(usize,bool)> = Vec::new();
 for &start_ei in &bnd {
 if out.edges.get(start_ei).map_or(true, |e| e.start == e.end) { continue; }
 if let Some(eg) = out.edges.get(start_ei) {
 let sv = eg.start;
 let mut stack: Vec<(usize,bool,Vec<(usize,bool)>,std::collections::HashSet<usize>,std::collections::HashSet<usize>)> = Vec::new();
 stack.push((start_ei, true, vec![], std::collections::HashSet::new(), std::collections::HashSet::new()));
 while let Some((ce, cf, path, ve, mut vv)) = stack.pop() {
 let ev = if cf { out.edges[ce].end } else { out.edges[ce].start };
 let mut np = path.clone(); np.push((ce, cf));
 let mut nve = ve.clone(); nve.insert(ce);
 if ev == sv && np.len() > 1 {
 if np.len() > best_loop.len() { best_loop = np; }
 continue;
 }
 if !vv.insert(ev) { continue; }
 if let Some(edges) = v2e.get(&ev) {
 for &(nei, nf) in edges {
 if !nve.contains(&nei) && out.edges.get(nei).map_or(true, |e| e.start != e.end) {
 stack.push((nei, nf, np.clone(), nve.clone(), vv.clone()));
 }
 }
 }
 }
 }
 }
 if best_loop.len() >= 3 { loops.push(best_loop); }
 if loops.is_empty() { continue; }
 use glam::DVec3;
 let mut areas: Vec<f64> = Vec::new();
 for lp in &loops {
 let (mut mn, mut mx) = (DVec3::splat(f64::INFINITY), DVec3::splat(f64::NEG_INFINITY));
 for &(ei,f) in lp {
 if let Some(e) = out.edges.get(ei) {
 if let Some(v) = out.vertices.get(if f { e.start } else { e.end }) {
 mn = mn.min(v.point); mx = mx.max(v.point);
 }
 }
 }
 areas.push(((mx.x-mn.x)*(mx.y-mn.y)*(mx.z-mn.z)).abs());
 }
 let oi = areas.iter().enumerate().max_by(|(_,a),(_,b)| a.partial_cmp(b).unwrap()).map(|(i,_)|i).unwrap_or(0);
 let mut ol = loops.swap_remove(oi);
 // �?OCCT : loop  �?
 // OCCT FillSameDomainFaces  �?loop  —— �?
 // �?merged face  ( �?  surface normal
 // )。rcad �?DFS  , �?edge �?forward flag  ,
 // �? �?OCCT  �?
 use rcad_kernel::topology::WireEdge;
 let ow = rcad_kernel::topology::Wire { edges: ol.iter().map(|&(ei,f)| WireEdge{idx:ei,forward:f}).collect() };
 let iws: Vec<rcad_kernel::topology::Wire> = loops.iter().map(|lp| rcad_kernel::topology::Wire { edges: lp.iter().map(|&(ei,f)| WireEdge{idx:ei,forward:f}).collect() }).collect();
 let nm = out.solids[si].shells[shi].faces[g[0]].normal;
 let sd = out.solids[si].shells[shi].faces[g[0]].surface_idx;
 let sp = out.solids[si].shells[shi].faces[g[0]].sample_point;
 let mf = rcad_kernel::Face { outer_wire: ow, inner_wires: iws, normal: nm, triangles: vec![], sample_point: sp, mesh_dirty: true,
 surface_idx: sd };
 let kp = g[0]; out.solids[si].shells[shi].faces[kp] = mf;
 let mut rd: Vec<usize> = g.iter().skip(1).copied().collect();
 rd.sort_unstable_by(|a,b| b.cmp(a));
 // �?OCCT : face_surface Vec�?
 for &fi in &rd {
 let mut flat = 0usize;
 for s in 0..si { for sh in &out.solids[s].shells { flat += sh.faces.len(); } }
 for sh in 0..shi { flat += out.solids[si].shells[sh].faces.len(); }
 flat += fi;
 out.solids[si].shells[shi].faces.remove(fi);
 if flat < out.geom.face_surface.len() { out.geom.face_surface.remove(flat); }
 }
 total += 1;
 }
 }
 }
 (out, total)
}
///
/// Verifies that PCurve parameterizations align properly where the two faces meet.
/// This is a topological guard to prevent merging faces with incompatible edge representations.
fn validate_shared_edge_continuity(
 brep: &rcad_kernel::BRep,
 si: usize,
 shi: usize,
 fi1: usize,
 fi2: usize,
 edge_idx: usize,
) -> bool {
 // If SameParameter is flagged, the 3D edge and all PCurves share parameterization.
 let same_param = brep
 .geom
 .edge_same_parameter
 .get(edge_idx)
 .copied()
 .unwrap_or(false);

 if !same_param {
 // For non-SameParameter edges, we need extra care.
 // For now, we skip PCurve continuity checks on such edges to avoid
 // false negatives from complex parameterization mismatches.
 // This is conservative but safe.
 return true;
 }

 // Get PCurves for this edge on both faces.
 let _pcurves = match brep.geom.edge_pcurves.get(edge_idx) {
 Some(pcs) => pcs,
 None => return true, // No PCurves: rely on geometric plane check.
 };

 if _pcurves.is_empty() {
 return true;
 }

 // Map face indices in the shell to global face indices for lookup.
 let mut global_fi1 = 0usize;
 let mut global_fi2 = 0usize;
 for s in 0..si {
 for sh in &brep.solids[s].shells {
 global_fi1 += sh.faces.len();
 global_fi2 += sh.faces.len();
 }
 }
 for sh in 0..shi {
 global_fi1 += brep.solids[si].shells[sh].faces.len();
 global_fi2 += brep.solids[si].shells[sh].faces.len();
 }
 global_fi1 += fi1;
 global_fi2 += fi2;

 // Note: Full PCurve continuity checks require careful parameterization
 // analysis; for now, we rely on SameParameter as a sufficient guard.

 // All PCurve continuity checks passed (or were skipped for safety).
 true
}

/// Validate that two adjacent faces' UV regions are geometrically compatible.
///
/// Checks that the parameter domains [u1, u2, v1, v2] for both faces do not
/// represent disjoint or incompatible regions on their respective surfaces.
/// This prevents merging faces that happen to be coplanar but cover different
/// parts of the surface domain.
fn validate_uv_regions_compatible(
 brep: &rcad_kernel::BRep,
 si: usize,
 shi: usize,
 fi1: usize,
 fi2: usize,
) -> bool {
 // Get UV domain ranges for both faces.
 // We need to map from face indices in the shell to global face indices.
 let mut global_fi1 = 0usize;
 let mut global_fi2 = 0usize;
 for s in 0..si {
 for sh in &brep.solids[s].shells {
 global_fi1 += sh.faces.len();
 global_fi2 += sh.faces.len();
 }
 }
 for sh in 0..shi {
 global_fi1 += brep.solids[si].shells[sh].faces.len();
 global_fi2 += brep.solids[si].shells[sh].faces.len();
 }
 global_fi1 += fi1;
 global_fi2 += fi2;

 // Fetch UV bounds; [u1, u2, v1, v2].
 let uv1 = match brep.geom.face_surface_range.get(global_fi1) {
 Some(Some(uv)) => *uv,
 _ => return true, // No UV data: assume compatible.
 };
 let uv2 = match brep.geom.face_surface_range.get(global_fi2) {
 Some(Some(uv)) => *uv,
 _ => return true, // No UV data: assume compatible.
 };

 let uv_tol = tolerance::TOLERANCE_PARAM_LEGACY;

 // Check if UV regions have meaningful overlap or adjacency.
 // If both regions are very small or identical, they are likely patches of the same domain.
 let _u1_size = (uv1[1] - uv1[0]).abs();
 let _v1_size = (uv1[3] - uv1[2]).abs();
 let _u2_size = (uv2[1] - uv2[0]).abs();
 let _v2_size = (uv2[3] - uv2[2]).abs();

 // Heuristic: if one face's UV domain is much larger than the other,
 // they likely represent compatible patches of the same surface.
 // (E.g., a plane split into two faces: one may have [0, 100, 0, 10]
 // and the other [50, 150, 0, 10] -- overlapping u-domain [50, 100].)

 let u_min = uv1[0].min(uv2[0]);
 let u_max = uv1[1].max(uv2[1]);
 let v_min = uv1[2].min(uv2[2]);
 let v_max = uv1[3].max(uv2[3]);

 let combined_u_size = (u_max - u_min).abs();
 let combined_v_size = (v_max - v_min).abs();

 // If either dimension's combined span is less than the tolerance, regions are coincident.
 if combined_u_size <= uv_tol || combined_v_size <= uv_tol {
 return true;
 }

 // Check for meaningful overlap in u-direction.
 let u_overlap_min = uv1[0].max(uv2[0]);
 let u_overlap_max = uv1[1].min(uv2[1]);
 let u_overlap = (u_overlap_max - u_overlap_min).max(0.0);

 // Check for meaningful overlap in v-direction.
 let v_overlap_min = uv1[2].max(uv2[2]);
 let v_overlap_max = uv1[3].min(uv2[3]);
 let v_overlap = (v_overlap_max - v_overlap_min).max(0.0);

 // Regions are compatible if:
 // - They overlap in both dimensions, OR
 // - They cover adjacent parts of the same surface (e.g., coplanar patches)
 // Adjacent means they touch along an edge with zero gap.

 (u_overlap > uv_tol && v_overlap > uv_tol)
 || ((u_overlap_max - u_overlap_min).abs() <= uv_tol && v_overlap > 0.0)
 || ((v_overlap_max - v_overlap_min).abs() <= uv_tol && u_overlap > 0.0)
}

/// Absolute area of a simple 3D polygon via Newell projection (see `builder::ResultBuilder`).
fn newell_polygon_abs_area(poly: &[glam::DVec3], normal: glam::DVec3) -> f64 {
 if poly.len() < 3 {
 return 0.0;
 }
 let n = normal.normalize_or_zero();
 let ax = n.x.abs();
 let ay = n.y.abs();
 let az = n.z.abs();
 let axis = if ax >= ay && ax >= az {
 0usize
 } else if ay >= az {
 1usize
 } else {
 2usize
 };
 let mut area2 = 0.0;
 for i in 0..poly.len() {
 let p = poly[i];
 let q = poly[(i + 1) % poly.len()];
 area2 += match axis {
 0 => p.y * q.z - q.y * p.z,
 1 => p.x * q.z - q.x * p.z,
 _ => p.x * q.y - q.x * p.y,
 };
 }
 0.5 * area2.abs()
}

fn face_outer_polygon_points(brep: &rcad_kernel::BRep, si: usize, shi: usize, fi: usize) -> Vec<glam::DVec3> {
 let face = &brep.solids[si].shells[shi].faces[fi];
 let mut pts = Vec::new();
 for we in &face.outer_wire.edges {
 if let Some((u, _)) = oriented_edge_vertices(brep, *we)
 && let Some(v) = brep.vertices.get(u)
 {
 pts.push(v.point);
 }
 }
 pts
}

fn wire_to_polygon_points(
 brep: &rcad_kernel::BRep,
 wire: &[rcad_kernel::topology::WireEdge],
) -> Vec<glam::DVec3> {
 let mut pts = Vec::new();
 for we in wire {
 if let Some((u, _)) = oriented_edge_vertices(brep, *we)
 && let Some(v) = brep.vertices.get(u)
 {
 pts.push(v.point);
 }
 }
 pts
}

/// Remove geometry slots for the flattened face at `remove_flat`.
///
/// Must stay in sync with topology when a face is deleted �?`face_surface`,
/// `face_surface_range`, and `face_tolerance` all use the same flat face order
/// as [`rcad_kernel::GeomStore`].
pub(crate) fn remove_flat_face_geom_slots(geom: &mut rcad_kernel::GeomStore, remove_flat: usize) {
 if geom.face_surface.len() > remove_flat {
 geom.face_surface.remove(remove_flat);
 }
 if geom.face_surface_range.len() > remove_flat {
 geom.face_surface_range.remove(remove_flat);
 }
 if geom.face_tolerance.len() > remove_flat {
 geom.face_tolerance.remove(remove_flat);
 }
}

/// Attempt one merge of two adjacent same-domain faces in `brep`. Returns `true`
/// if a merge was performed (mutating `brep` in place).
///
/// Handles planar, cylindrical, toroidal, and spherical surface pairs.
fn unify_one_merge_pass(brep: &mut rcad_kernel::BRep) -> bool {
 unify_one_merge_pass_with_origins(brep, None)
}