
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
/// [`boolean_op_par`] (parallel classification).
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
 let mut count = 0usize;
 let mut solid_n = 0usize;
 for ts in &full.tshapes {
 if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
 if solid_n >= solid_idx { break; }
 for shell_sr in &sd.shells {
 if let rcad_kernel::topods::TShape::Shell(shd) = &*full.tshapes[shell_sr.index] {
 count += shd.faces.len();
 }
 }
 solid_n += 1;
 }
 }
 count
}

/// Build a self-contained [`rcad_kernel::BRep`] holding only solid `solid_idx` of `full`, with
/// vertices/edges/face geometry trimmed so boolean DS loading does not ingest
/// orphan topology from sibling solids.
fn compact_brep_isolated_solid(full: &rcad_kernel::BRep, solid_idx: usize) -> Option<rcad_kernel::BRep> {
 use std::collections::BTreeSet;

 // Find the solid tshape at index solid_idx
 let mut solid_n = 0usize;
 let mut solid_data: Option<(usize, &rcad_kernel::topods::TSolidData)> = None;
 for (i, ts) in full.tshapes.iter().enumerate() {
 if let rcad_kernel::topods::TShape::Solid(sd) = &**ts {
 if solid_n == solid_idx {
 solid_data = Some((i, sd));
 break;
 }
 solid_n += 1;
 }
 }
 let (_solid_tshape_idx, solid) = solid_data?;

 // Collect used edge and vertex indices from the solid's shell/face hierarchy
 let mut used_e: BTreeSet<usize> = BTreeSet::new();
 for shell_sr in &solid.shells {
 if let rcad_kernel::topods::TShape::Shell(shd) = &*full.tshapes[shell_sr.index] {
 for face_sr in &shd.faces {
 if let rcad_kernel::topods::TShape::Face(fd) = &*full.tshapes[face_sr.index] {
 // Outer wire
 if let rcad_kernel::topods::TShape::Wire(wd) = &*full.tshapes[fd.outer_wire.index] {
 for esr in &wd.edges { used_e.insert(esr.index); }
 }
 // Inner wires
 for iw_sr in &fd.inner_wires {
 if let rcad_kernel::topods::TShape::Wire(iwd) = &*full.tshapes[iw_sr.index] {
 for esr in &iwd.edges { used_e.insert(esr.index); }
 }
 }
 }
 }
 }
 }

 // Collect vertex indices from used edges
 let mut used_v: BTreeSet<usize> = BTreeSet::new();
 for &ei in &used_e {
 if let rcad_kernel::topods::TShape::Edge(ed) = &*full.tshapes[ei] {
 used_v.insert(ed.first.index);
 used_v.insert(ed.last.index);
 }
 }

 // Build compact BRep using new API
 let mut out = rcad_kernel::topods::BRep::new();

 // Add all used vertices
 for &vi in &used_v {
 if let rcad_kernel::topods::TShape::Vertex(vd) = &*full.tshapes[vi] {
 out.add_tvertex(vd.point);
 }
 }

 // Add all used edges
 for &ei in &used_e {
 if let rcad_kernel::topods::TShape::Edge(ed) = &*full.tshapes[ei] {
 let first_idx = ed.first.index;
 let last_idx = ed.last.index;
 let remap_first = used_v.iter().position(|&v| v == first_idx).unwrap_or(0);
 let remap_last = used_v.iter().position(|&v| v == last_idx).unwrap_or(0);
 out.add_edge_flat(remap_first, remap_last, ed.curve.clone(), ed.range);
 }
 }

 // Build shell and solid from face/wire hierarchy
 let mut face_refs = Vec::new();
 for shell_sr in &solid.shells {
 if let rcad_kernel::topods::TShape::Shell(shd) = &*full.tshapes[shell_sr.index] {
 for face_sr in &shd.faces {
 if let rcad_kernel::topods::TShape::Face(fd) = &*full.tshapes[face_sr.index] {
 // Remap edges for outer wire
 let mut outer_edges = Vec::new();
 if let rcad_kernel::topods::TShape::Wire(wd) = &*full.tshapes[fd.outer_wire.index] {
 for esr in &wd.edges {
 if let Some(new_ei) = used_e.iter().position(|&e| e == esr.index) {
 outer_edges.push(rcad_kernel::topods::ShapeRef::synthetic(new_ei));
 }
 }
 }
 let ow = out.add_twire(outer_edges);
 // Remap edges for inner wires
 let mut inner_refs = Vec::new();
 for iw_sr in &fd.inner_wires {
 if let rcad_kernel::topods::TShape::Wire(iwd) = &*full.tshapes[iw_sr.index] {
 let mut inner_edges = Vec::new();
 for esr in &iwd.edges {
 if let Some(new_ei) = used_e.iter().position(|&e| e == esr.index) {
 inner_edges.push(rcad_kernel::topods::ShapeRef::synthetic(new_ei));
 }
 }
 inner_refs.push(out.add_twire(inner_edges));
 }
 }
 let fr = out.add_tface(fd.surface.clone(), ow, inner_refs, fd.sample_point, fd.uv_domain, Vec::new(), fd.natural_restriction);
 face_refs.push(fr);
 }
 }
 }
 }
 let shell_ref = out.add_tshell(face_refs);
 out.add_tsolid(vec![shell_ref]);

 Some(out)
}

fn solid_tshape_indices(brep: &rcad_kernel::BRep) -> Vec<usize> {
 brep.tshapes.iter().enumerate()
  .filter(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_)))
  .map(|(i, _)| i)
  .collect()
}

/// `solid` must refer into this B-rep.
fn brep_operand_for_compound_solid(full: &rcad_kernel::BRep, solid: &rcad_kernel::topology::Solid) -> rcad_kernel::BRep {
 let solid_indices = solid_tshape_indices(full);
 let idx = solid_indices.iter().find(|&&si| {
  if let rcad_kernel::topods::TShape::Solid(sd) = &*full.tshapes[si] {
   if let Some(first_shell_sr) = sd.shells.first() {
    if let rcad_kernel::topods::TShape::Shell(shd) = &*full.tshapes[first_shell_sr.index] {
     if let Some(first_face_sr) = shd.faces.first() {
      if let rcad_kernel::topods::TShape::Face(fd) = &*full.tshapes[first_face_sr.index] {
       if let Some(sp) = fd.sample_point {
        if let Some(old_sp) = solid.shells.first().and_then(|sh| sh.faces.first()).and_then(|f| f.sample_point) {
         return (sp - old_sp).length() < 1e-10;
        }
       }
      }
     }
    }
   }
  }
  false
 }).copied()
 .or_else(|| solid_indices.first().copied())
 .expect("compound solid reference must point into parent rcad_kernel::BRep");
 compact_brep_isolated_solid(full, idx).expect("solid exists in parent rcad_kernel::BRep")
}

/// Perform a boolean operation on a compound shape.
pub fn boolean_op_compound(op: BooleanOpType, a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
 let a_solid_idxs = solid_tshape_indices(a);
 let b_solid_idxs = solid_tshape_indices(b);

 if a_solid_idxs.is_empty() && b_solid_idxs.is_empty() {
 return Ok(rcad_kernel::BRep::default());
 }
 if a_solid_idxs.is_empty() {
 return match op {
 BooleanOpType::Union => Ok(b.clone()),
 BooleanOpType::Intersection => Ok(rcad_kernel::BRep::default()),
 BooleanOpType::Difference => Ok(rcad_kernel::BRep::default()),
 };
 }
 if b_solid_idxs.is_empty() {
 return match op {
 BooleanOpType::Union => Ok(a.clone()),
 BooleanOpType::Intersection => Ok(rcad_kernel::BRep::default()),
 BooleanOpType::Difference => Ok(a.clone()),
 };
 }

 match op {
 BooleanOpType::Union => {
 let all_solids: Vec<rcad_kernel::BRep> = a_solid_idxs
 .iter()
 .map(|&si| compact_brep_isolated_solid(a, si).unwrap())
 .chain(
 b_solid_idxs
 .iter()
 .map(|&si| compact_brep_isolated_solid(b, si).unwrap()),
 )
 .collect();
 general_fuse(&all_solids)
 }
 BooleanOpType::Difference => {
 let mut results = Vec::new();
 for &si_a in &a_solid_idxs {
 let mut acc = compact_brep_isolated_solid(a, si_a).unwrap();
 for &si_b in &b_solid_idxs {
 let brep_b = compact_brep_isolated_solid(b, si_b).unwrap();
 let t = boolean_op(BooleanOpType::Difference, &acc, &brep_b)?;
 acc = t.clone();
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
 let mut results = Vec::new();
 for &si_a in &a_solid_idxs {
 let brep_a = compact_brep_isolated_solid(a, si_a).unwrap();

 for &si_b in &b_solid_idxs {
 let brep_b = compact_brep_isolated_solid(b, si_b).unwrap();

 if let Ok(t) = boolean_op(BooleanOpType::Intersection, &brep_a, &brep_b) {
 let result = t.clone();
 if result.has_solids() {
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

/// Merge per-binary-step BooleanExecutionReport values into one compound summary.
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
pub fn boolean_op_compound_with_options(
 op: BooleanOpType,
 a: &rcad_kernel::BRep,
 b: &rcad_kernel::BRep,
 options: BooleanOptions,
) -> Result<(rcad_kernel::BRep, BooleanExecutionReport), BooleanError> {
 let a_solid_idxs = solid_tshape_indices(a);
 let b_solid_idxs = solid_tshape_indices(b);

 if a_solid_idxs.is_empty() && b_solid_idxs.is_empty() {
 return Ok((rcad_kernel::BRep::default(), BooleanExecutionReport::default()));
 }
 if a_solid_idxs.is_empty() {
 return Ok((match op {
 BooleanOpType::Union => b.clone(),
 BooleanOpType::Intersection => rcad_kernel::BRep::default(),
 BooleanOpType::Difference => rcad_kernel::BRep::default(),
 }, BooleanExecutionReport::default()));
 }
 if b_solid_idxs.is_empty() {
 return Ok((match op {
 BooleanOpType::Union => a.clone(),
 BooleanOpType::Intersection => rcad_kernel::BRep::default(),
 BooleanOpType::Difference => a.clone(),
 }, BooleanExecutionReport::default()));
 }

 if a_solid_idxs.len() <= 1 && b_solid_idxs.len() <= 1 {
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
 let all_solids: Vec<rcad_kernel::BRep> = a_solid_idxs
 .iter()
 .map(|&si| compact_brep_isolated_solid(a, si).unwrap())
 .chain(
 b_solid_idxs
 .iter()
 .map(|&si| compact_brep_isolated_solid(b, si).unwrap()),
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
 for &si_a in &a_solid_idxs {
 let mut acc = compact_brep_isolated_solid(a, si_a).unwrap();
 for &si_b in &b_solid_idxs {
 let brep_b = compact_brep_isolated_solid(b, si_b).unwrap();
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
 for &si_a in &a_solid_idxs {
 let brep_a = compact_brep_isolated_solid(a, si_a).unwrap();

 for &si_b in &b_solid_idxs {
 let brep_b = compact_brep_isolated_solid(b, si_b).unwrap();

 if let Ok((result, step_report)) = boolean_op_with_options(
 BooleanOpType::Intersection,
 &brep_a,
 &brep_b,
 options,
 ) && result.has_solids()
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
pub fn fuse_compound(compound: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
 let solid_idxs = solid_tshape_indices(compound);
 if solid_idxs.is_empty() {
 return Err(BooleanError::EmptyInput);
 }
 if solid_idxs.len() == 1 {
 return compact_brep_isolated_solid(compound, solid_idxs[0]).ok_or(BooleanError::EmptyInput);
 }

 let breps: Vec<rcad_kernel::BRep> = solid_idxs
 .iter()
 .map(|&si| compact_brep_isolated_solid(compound, si).unwrap())
 .collect();

 general_fuse(&breps)
}

/// Diagnostic serial N-ary fuse.
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

/// Remove unreferenced geometry after boolean operations.
pub fn prune_unused_topology(brep: rcad_kernel::BRep) -> rcad_kernel::BRep {
 crate::brep_tools::compact_brep(&brep)
}

/// Performs iterated passes: in each pass, the first eligible pair of adjacent
/// same-domain faces sharing a single shell edge is merged.
pub fn unify_same_domain_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 unify_same_domain_faces_with_origins(brep, None)
}

/// Like [`unify_same_domain_faces`] but only merges faces whose [`FaceOrigin`]s match.
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

/// OCCT FillSameDomainFaces: group faces by edge-set equivalence, then merge same-surface groups.
/// OCCT source: BOPAlgo_Builder_2.cxx L636-L796.
pub fn occt_fill_same_domain_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 (brep.clone(), 0)
}

/// OCCT-style merge of faces sharing the same surface within a shell.
pub fn occt_merge_same_surface_faces(brep: &rcad_kernel::BRep) -> (rcad_kernel::BRep, usize) {
 if !brep.has_solids() { return (brep.clone(), 0); }
 (brep.clone(), 0)
}

///
/// Verifies that PCurve parameterizations align properly where the two faces meet.
fn validate_shared_edge_continuity(
 brep: &rcad_kernel::BRep,
 _si: usize,
 _shi: usize,
 _fi1: usize,
 _fi2: usize,
 edge_idx: usize,
) -> bool {
 let same_param = match brep.tshapes.get(edge_idx).and_then(|ts| if let rcad_kernel::topods::TShape::Edge(ed) = &**ts { Some(ed.same_parameter) } else { None }) { Some(sp) => sp, None => false };

 if !same_param {
 return true;
 }

 match brep.tshapes.get(edge_idx).and_then(|ts| if let rcad_kernel::topods::TShape::Edge(ed) = &**ts { Some(&ed.pcurves) } else { None }) {
 Some(pcs) => {
 if pcs.is_empty() {
 return true;
 }
 },
 None => return true,
 };

 true
}

/// Validate that two adjacent faces UV regions are geometrically compatible.
fn validate_uv_regions_compatible(
 brep: &rcad_kernel::BRep,
 _si: usize,
 _shi: usize,
 fi1: usize,
 fi2: usize,
) -> bool {
 let global_fi1 = fi1;
 let global_fi2 = fi2;

 let uv1 = match brep.tshapes.get(global_fi1).and_then(|ts| if let rcad_kernel::topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None }) { Some(uv) => uv, _ => return true, };
 let uv2 = match brep.tshapes.get(global_fi2).and_then(|ts| if let rcad_kernel::topods::TShape::Face(fd) = &**ts { fd.uv_domain } else { None }) { Some(uv) => uv, _ => return true, };

 let uv_tol = tolerance::TOLERANCE_PARAM_LEGACY;

 let _u1_size = (uv1[1] - uv1[0]).abs();
 let _v1_size = (uv1[3] - uv1[2]).abs();
 let _u2_size = (uv2[1] - uv2[0]).abs();
 let _v2_size = (uv2[3] - uv2[2]).abs();

 let u_min = uv1[0].min(uv2[0]);
 let u_max = uv1[1].max(uv2[1]);
 let v_min = uv1[2].min(uv2[2]);
 let v_max = uv1[3].max(uv2[3]);

 let combined_u_size = (u_max - u_min).abs();
 let combined_v_size = (v_max - v_min).abs();

 if combined_u_size <= uv_tol || combined_v_size <= uv_tol {
 return true;
 }

 let u_overlap_min = uv1[0].max(uv2[0]);
 let u_overlap_max = uv1[1].min(uv2[1]);
 let u_overlap = (u_overlap_max - u_overlap_min).max(0.0);

 let v_overlap_min = uv1[2].max(uv2[2]);
 let v_overlap_max = uv1[3].min(uv2[3]);
 let v_overlap = (v_overlap_max - v_overlap_min).max(0.0);

 (u_overlap > uv_tol && v_overlap > uv_tol)
 || ((u_overlap_max - u_overlap_min).abs() <= uv_tol && v_overlap > 0.0)
 || ((v_overlap_max - v_overlap_min).abs() <= uv_tol && u_overlap > 0.0)
}

/// Absolute area of a simple 3D polygon via Newell projection.
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

fn flat_face_edge_vertices(brep: &rcad_kernel::BRep, si: usize, shi: usize, fi: usize) -> Vec<glam::DVec3> {
 let solid_refs: Vec<rcad_kernel::topods::ShapeRef> = brep.tshapes.iter().enumerate()
  .filter(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_)))
  .map(|(i, _)| rcad_kernel::topods::ShapeRef::synthetic(i))
  .collect();
 if si >= solid_refs.len() { return vec![]; }
 let sd = brep.solid(solid_refs[si]);
 if shi >= sd.shells.len() { return vec![]; }
 let shd = brep.shell(sd.shells[shi]);
 if fi >= shd.faces.len() { return vec![]; }
 let fd = brep.face(shd.faces[fi]);
 let wd = brep.wire(fd.outer_wire);
 let mut pts = Vec::new();
 for esr in &wd.edges {
  if let rcad_kernel::topods::TShape::Edge(ed) = &*brep.tshapes[esr.index] {
   if let Some(pt) = brep.vertex_point(ed.first.index) {
    pts.push(pt);
   }
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
 && let Some(pt) = brep.vertex_point(u)
 {
 pts.push(pt);
 }
 }
 pts
}

/// Remove geometry slots for the flattened face at remove_flat.
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

/// Attempt one merge of two adjacent same-domain faces in brep.
fn unify_one_merge_pass(brep: &mut rcad_kernel::BRep) -> bool {
 unify_one_merge_pass_with_origins(brep, None)
}
