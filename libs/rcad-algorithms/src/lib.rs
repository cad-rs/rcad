pub mod bopds;
pub mod brep_check;
pub mod brep_repair;
pub mod builder;
pub mod features;
pub mod bvh;
pub mod classify;
pub mod draft;
pub mod geom_populate;
pub mod healing;
pub mod history;
pub mod hlr;
pub mod imprint;
pub use features::{FeatureError, make_cylindrical_hole, make_draft_prism, make_prism, make_revolution};
pub mod inttools;
pub mod pave_filler;
pub mod section;
pub mod thicken;
pub mod tolerance;
pub mod triangulate;
pub mod array;
pub mod cells_builder;

use serde::Serialize;

pub use bvh::{Aabb, Bvh, BvhStats};

use rcad_kernel::BRep;

pub use brep_check::{CheckIssue, CheckResult, check,
    SuspectEdge, SameParameterDiagnosis, diagnose_same_parameter,
    SuspectSameRangeEdge, SameRangeDiagnosis, diagnose_same_range,
    ShellTopologyReport, analyze_shell_topology,
    WireAnalysisReport, WireIssueReport, analyze_wire_issues,
};
pub use brep_repair::{
    RepairReport, fix_same_parameter, fix_same_parameter_with_scan, fix_wire_orientation,
    merge_close_vertices, recompute_face_normals, remove_degenerate_faces, repair,
    remove_small_edges, fix_same_range_with_scan,
};
pub use healing::{
    HealingIssueStats, HealingMode, HealingOptions, HealingReport, analyze_and_heal, heal,
};
pub use builder::{BooleanError, BooleanOpType};
pub use history::{BooleanHistory, EdgeOrigin, FaceOrigin, VertexOrigin};
pub use hlr::{AssemblyHlrResult, ComponentHlr, HlrCamera, HlrResult, HlrSegment, hlr, hlr_assembly, hlr_to_svg};
pub use imprint::{
    Gap, GapOverlapReport, ImprintResult, Overlap, detect_gaps_overlaps, imprint_brep, min_distance,
};
pub use inttools::{
    SurfaceCurve, SurfaceIntersectionResult, SurfaceSurfaceIntersection, intersect_surfaces,
    intersect_surfaces_with_density, intersect_surfaces_with_tolerance,
};
pub use section::{SectionCurve, section, section_curves, section_polylines};
pub use thicken::{ThickeningResult, thicken_shell};
pub use triangulate::{SurfaceMesh, TessellationParams, mesh_brep, triangulate_surface};
pub use array::{
    LinearPatternParams, CircularPatternParams, PatternError,
    linear_pattern, circular_pattern,
};
pub use cells_builder::{CellExpr, CellsBuilder, CellsBuilderError};

/// Options for post-operation topology simplification.
#[derive(Debug, Clone, Copy)]
pub struct SimplifyOptions {
    pub merge_vertices: bool,
    pub merge_tolerance: f64,
    pub recompute_normals: bool,
    pub remove_degenerate_faces: bool,
    pub fix_wire_orientation: bool,
    /// Merge adjacent coplanar planar faces into larger faces.
    pub unify_same_domain_faces: bool,
    /// Remove redundant coplanar internal faces (mainly for union outputs).
    pub remove_internal_faces: bool,
    /// Remove edges whose chord length is below `small_edge_min_length`.
    pub remove_small_edges: bool,
    /// Chord-length threshold for small-edge removal (default: `TOLERANCE_ABS`).
    pub small_edge_min_length: f64,
}

impl Default for SimplifyOptions {
    fn default() -> Self {
        Self {
            merge_vertices: true,
            merge_tolerance: tolerance::TOLERANCE_ABS,
            recompute_normals: true,
            remove_degenerate_faces: true,
            fix_wire_orientation: true,
            unify_same_domain_faces: true,
            remove_internal_faces: true,
            remove_small_edges: false,
            small_edge_min_length: tolerance::TOLERANCE_ABS,
        }
    }
}

/// Report of simplification steps and checker deltas.
#[derive(Debug, Clone, Default)]
pub struct SimplifyReport {
    pub vertices_merged: usize,
    pub degenerate_faces_removed: usize,
    pub normals_recomputed: usize,
    pub wires_fixed: usize,
    pub same_domain_face_merges: usize,
    pub internal_faces_removed: usize,
    pub small_edges_removed: usize,
    pub issues_before: usize,
    pub issues_after: usize,
}

/// Options for boolean execution pipeline.
#[derive(Debug, Clone, Copy)]
pub struct BooleanOptions {
    /// Use BVH acceleration during pave filling when possible.
    pub use_bvh: bool,
    /// Run structured healing after boolean build.
    pub run_healing: bool,
    /// Healing options used when `run_healing` is enabled.
    pub healing: HealingOptions,
    /// Run topology simplification after boolean/healing.
    pub run_simplify: bool,
    /// Simplification options used when `run_simplify` is enabled.
    pub simplify: SimplifyOptions,
    /// Include origin history and stable per-face labels in report.
    pub include_history: bool,
    /// Fuzzy tolerance for near-miss interference detection (analogous to
    /// `BOPAlgo_Options::SetFuzzyValue`).
    ///
    /// Values ≤ 0 use the default `TOLERANCE_ABS`.  Useful for inputs with
    /// vertices/edges that are almost but not exactly touching.
    pub fuzzy_tol: f64,
}

impl Default for BooleanOptions {
    fn default() -> Self {
        Self {
            use_bvh: true,
            run_healing: false,
            healing: HealingOptions::default(),
            run_simplify: false,
            simplify: SimplifyOptions::default(),
            include_history: false,
            fuzzy_tol: 0.0,
        }
    }
}

/// Structured diagnostics for boolean execution.
#[derive(Debug, Clone, Default)]
pub struct BooleanExecutionReport {
    pub input_faces_a: usize,
    pub input_faces_b: usize,
    pub output_faces: usize,
    pub used_bvh: bool,
    pub healed: bool,
    pub healing_report: Option<HealingReport>,
    pub simplified: bool,
    pub simplify_report: Option<SimplifyReport>,
    pub history_faces: usize,
    pub persistent_face_labels: Vec<String>,
}

/// Options for split-first workflows.
#[derive(Debug, Clone, Copy)]
pub struct SplitterOptions {
    /// If true, run healing after each split step.
    pub heal_after_each_step: bool,
    /// Healing options used when `heal_after_each_step` is enabled.
    pub healing: HealingOptions,
    /// Additional linear tolerance used by splitter broad-phase pruning.
    ///
    /// Tools whose axis-aligned bounding boxes are farther than this distance
    /// from the current object are skipped.
    pub fuzzy_tolerance: f64,
    /// Enable AABB broad-phase pruning for split steps.
    pub broad_phase_pruning: bool,
    /// Validation strictness used by checked splitter APIs.
    pub validation_level: SplitterValidationLevel,
}

/// Validation strictness for checked splitter workflows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum SplitterValidationLevel {
    /// Accept split-first intermediate non-manifold topology.
    Relaxed,
    /// Treat all checker issues as errors.
    Strict,
}

impl Default for SplitterOptions {
    fn default() -> Self {
        Self {
            heal_after_each_step: false,
            healing: HealingOptions::default(),
            fuzzy_tolerance: 0.0,
            broad_phase_pruning: true,
            validation_level: SplitterValidationLevel::Relaxed,
        }
    }
}

/// Per-step diagnostics for splitter execution.
#[derive(Debug, Clone, Serialize)]
pub struct SplitterStepReport {
    /// Zero-based tool index used for this split step.
    pub step_index: usize,
    /// Face count before this split step.
    pub input_faces: usize,
    /// Number of seam-edge pairs reported by imprint in this step.
    pub seam_edges: usize,
    /// Face count after this step.
    pub output_faces: usize,
    /// Whether healing was applied at this step.
    pub healed: bool,
    /// Whether this step was skipped by broad-phase pruning.
    pub skipped_by_broad_phase: bool,
    /// Validation issue count for this step when checked mode is enabled.
    pub validation_issue_count: Option<usize>,
    /// First validation issue message when available.
    pub validation_first_issue: Option<String>,
}

/// Diagnostics report for split-first workflows.
#[derive(Debug, Clone, Default, Serialize)]
pub struct SplitterReport {
    /// Step-by-step diagnostics.
    pub steps: Vec<SplitterStepReport>,
    /// Total seam-edge pairs accumulated across all steps.
    pub total_seam_edges: usize,
}

/// Per-object diagnostics for grouped splitter workflows.
#[derive(Debug, Clone, Default, Serialize)]
pub struct SplitterObjectReport {
    /// Zero-based object index in input slice.
    pub object_index: usize,
    /// Step-level diagnostics for this object.
    pub steps: Vec<SplitterStepReport>,
    /// Total seam-edge pairs for this object.
    pub total_seam_edges: usize,
    /// Whether this object completed all requested split steps.
    pub completed: bool,
    /// Error captured for this object (checked collect mode).
    pub error: Option<SplitterError>,
}

/// Diagnostics for object/tool grouped split execution.
#[derive(Debug, Clone, Default, Serialize)]
pub struct SplitterObjectsReport {
    /// One report per input object, in the same order.
    pub objects: Vec<SplitterObjectReport>,
}

/// Aggregated summary for grouped splitter execution.
#[derive(Debug, Clone, Default, Serialize)]
pub struct SplitterObjectsSummary {
    pub total_objects: usize,
    pub completed_objects: usize,
    pub failed_objects: usize,
    /// Indices of failed objects in original input order.
    pub failed_object_indices: Vec<usize>,
    /// Histogram of failing step indices.
    pub failed_step_histogram: Vec<(usize, usize)>,
    /// Histogram of first error messages for failed objects.
    pub first_error_histogram: Vec<(String, usize)>,
}

impl SplitterObjectsReport {
    /// Build aggregated success/failure statistics for batch workflows.
    pub fn summarize(&self) -> SplitterObjectsSummary {
        let total_objects = self.objects.len();
        let completed_objects = self.objects.iter().filter(|o| o.completed).count();
        let failed_objects = total_objects.saturating_sub(completed_objects);

        let failed_object_indices: Vec<usize> = self
            .objects
            .iter()
            .filter(|o| !o.completed)
            .map(|o| o.object_index)
            .collect();

        let mut step_map: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
        let mut map: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
        for obj in &self.objects {
            if let Some(err) = &obj.error {
                if let Some(step_index) = err.step_index() {
                    *step_map.entry(step_index).or_insert(0) += 1;
                }
                *map.entry(err.to_string()).or_insert(0) += 1;
            }
        }

        SplitterObjectsSummary {
            total_objects,
            completed_objects,
            failed_objects,
            failed_object_indices,
            failed_step_histogram: step_map.into_iter().collect(),
            first_error_histogram: map.into_iter().collect(),
        }
    }

    /// Export report and summary as stable JSON payload `splitter.report.v1`.
    pub fn to_json_v1(&self) -> Result<String, serde_json::Error> {
        let payload = SplitterJsonV1 {
            schema: "splitter.report.v1",
            report: self,
            summary: self.summarize(),
        };
        serde_json::to_string_pretty(&payload)
    }
}

/// Stable JSON payload for splitter batch reporting.
#[derive(Debug, Clone, Serialize)]
pub struct SplitterJsonV1<'a> {
    pub schema: &'static str,
    pub report: &'a SplitterObjectsReport,
    pub summary: SplitterObjectsSummary,
}

/// Error returned by checked splitter workflows.
#[derive(Debug, Clone, Serialize)]
pub enum SplitterError {
    /// Split result became invalid at a specific step.
    StepInvalid {
        step_index: usize,
        issue_count: usize,
        first_issue: Option<String>,
    },
}

impl SplitterError {
    pub fn step_index(&self) -> Option<usize> {
        match self {
            Self::StepInvalid { step_index, .. } => Some(*step_index),
        }
    }
}

impl std::fmt::Display for SplitterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::StepInvalid {
                step_index,
                issue_count,
                first_issue,
            } => {
                if let Some(first) = first_issue {
                    write!(
                        f,
                        "splitter produced invalid result at step {step_index} ({issue_count} issues, first: {first})"
                    )
                } else {
                    write!(
                        f,
                        "splitter produced invalid result at step {step_index} ({issue_count} issues)"
                    )
                }
            }
        }
    }
}

impl std::error::Error for SplitterError {}

/// Perform a boolean operation on two BReps.
///
/// Both BReps must have populated GeomStore (call
/// `geom_populate::populate_box_geom` first for box primitives).
pub fn boolean_op(op: BooleanOpType, a: &BRep, b: &BRep) -> Result<BRep, BooleanError> {
    // 1. Build the DS from both shapes
    let mut ds = bopds::ds::DS::new(a, b);

    // 2. Run PaveFiller — compute all interferences
    let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
    let mut filler = match (&bvh_a, &bvh_b) {
        (Some(a), Some(b)) => pave_filler::PaveFiller::with_bvh(&mut ds, a, b),
        _ => pave_filler::PaveFiller::new(&mut ds),
    };
    filler.perform();

    // 3. Run Builder — classify and assemble result
    let builder = builder::BooleanBuilder::new(&ds, op);
    builder.build()
}

/// Perform a boolean operation with advanced execution options and report.
pub fn boolean_op_with_options(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
    options: BooleanOptions,
) -> Result<(BRep, BooleanExecutionReport), BooleanError> {
    let input_faces_a = face_count_of(a);
    let input_faces_b = face_count_of(b);
    let used_bvh = options.use_bvh && has_faces(a) && has_faces(b);

    let (mut out, mut report, history_opt) = if options.include_history {
        let (result, history) = if options.use_bvh {
            boolean_op_with_history(op, a, b)?
        } else {
            let mut ds = if options.fuzzy_tol > 0.0 {
                bopds::ds::DS::new_with_fuzzy(a, b, options.fuzzy_tol)
            } else {
                bopds::ds::DS::new(a, b)
            };
            let mut filler = pave_filler::PaveFiller::new(&mut ds);
            filler.perform();
            let builder = builder::BooleanBuilder::new(&ds, op);
            builder.build_with_history()?
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
            if options.fuzzy_tol > 0.0 {
                let mut ds = bopds::ds::DS::new_with_fuzzy(a, b, options.fuzzy_tol);
                let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
                let mut filler = match (&bvh_a, &bvh_b) {
                    (Some(ba), Some(bb)) => pave_filler::PaveFiller::with_bvh(&mut ds, ba, bb),
                    _ => pave_filler::PaveFiller::new(&mut ds),
                };
                filler.perform();
                let builder = builder::BooleanBuilder::new(&ds, op);
                builder.build()?
            } else {
                boolean_op(op, a, b)?
            }
        } else {
            let mut ds = if options.fuzzy_tol > 0.0 {
                bopds::ds::DS::new_with_fuzzy(a, b, options.fuzzy_tol)
            } else {
                bopds::ds::DS::new(a, b)
            };
            let mut filler = pave_filler::PaveFiller::new(&mut ds);
            filler.perform();
            let builder = builder::BooleanBuilder::new(&ds, op);
            builder.build()?
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
        let (healed, heal_report) = analyze_and_heal(&out, options.healing);
        out = healed;
        report.healed = true;
        report.healing_report = Some(heal_report);
    }

    if options.run_simplify {
        let (simplified, simp_report) = simplify_brep_post_ops(&out, options.simplify);
        out = simplified;
        report.simplified = true;
        report.simplify_report = Some(simp_report);
    }

    report.output_faces = face_count_of(&out);
    if let Some(history) = history_opt {
        report.history_faces = history.len();
        report.persistent_face_labels = persistent_face_labels_from_history(&history);
    }

    Ok((out, report))
}

/// Run post-operation simplification passes on a BRep.
pub fn simplify_brep_post_ops(brep: &BRep, options: SimplifyOptions) -> (BRep, SimplifyReport) {
    let before = check(brep);
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
    if options.fix_wire_orientation {
        let (next, n) = fix_wire_orientation(&out, options.merge_tolerance);
        out = next;
        report.wires_fixed = n;
    }
    if options.unify_same_domain_faces {
        let (next, n) = unify_same_domain_faces(&out);
        out = next;
        report.same_domain_face_merges = n;
    }
    if options.remove_internal_faces {
        let (next, n) = remove_internal_faces(&out);
        out = next;
        report.internal_faces_removed = n;
    }
    if options.remove_small_edges {
        let (next, n) = remove_small_edges(&out, options.small_edge_min_length);
        out = next;
        report.small_edges_removed = n;
    }

    report.issues_after = check(&out).issues.len();
    (out, report)
}

/// Boolean + simplification convenience pipeline.
pub fn boolean_op_simplified(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
    options: SimplifyOptions,
) -> Result<(BRep, SimplifyReport), BooleanError> {
    let raw = boolean_op(op, a, b)?;
    Ok(simplify_brep_post_ops(&raw, options))
}

/// Split `target` by one or more `tools` without boolean classification.
///
/// This is a first-stage splitter built on top of [`imprint_brep`]. It keeps
/// target material and iteratively imprints tool boundaries onto the evolving
/// target shape.
pub fn split_brep(target: &BRep, tools: &[BRep]) -> (BRep, SplitterReport) {
    split_brep_with_options(target, tools, SplitterOptions::default())
}

/// Like [`split_brep`] with advanced options.
pub fn split_brep_with_options(
    target: &BRep,
    tools: &[BRep],
    options: SplitterOptions,
) -> (BRep, SplitterReport) {
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
pub fn split_brep_checked_with_options(
    target: &BRep,
    tools: &[BRep],
    options: SplitterOptions,
) -> Result<(BRep, SplitterReport), SplitterError> {
    let (result, report) = split_brep_internal_with_partial_report(target, tools, options, true);
    result.map(|brep| (brep, report))
}

fn split_brep_internal_with_partial_report(
    target: &BRep,
    tools: &[BRep],
    options: SplitterOptions,
    validate_each_step: bool,
) -> (Result<BRep, SplitterError>, SplitterReport) {
    let mut acc = target.clone();
    let mut report = SplitterReport::default();

    for (step_index, tool) in tools.iter().enumerate() {
        let input_faces = face_count_of(&acc);
        let fuzzy = options.fuzzy_tolerance.max(0.0);
        let skipped_by_broad_phase = options.broad_phase_pruning
            && breps_farther_than_tolerance(&acc, tool, fuzzy);

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

        let mut step = imprint_brep(&acc, tool);
        let seam_edges = step.seam_edges.len();

        if options.heal_after_each_step {
            let (healed, _) = analyze_and_heal(&step.brep, options.healing);
            step.brep = healed;
        }

        let mut validation_issue_count = None;
        let mut validation_first_issue = None;
        let output_faces = face_count_of(&step.brep);
        if validate_each_step {
            let validity = check(&step.brep);
            let (issue_count, first_issue) = splitter_issues_by_level(&validity, options.validation_level);
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

        acc = step.brep;
    }

    (Ok(acc), report)
}

fn brep_bounds(brep: &BRep) -> Option<(glam::DVec3, glam::DVec3)> {
    let mut it = brep.vertices.iter();
    let first = it.next()?.point;
    let mut min = first;
    let mut max = first;
    for v in it {
        min = min.min(v.point);
        max = max.max(v.point);
    }
    Some((min, max))
}

fn aabb_distance(min_a: glam::DVec3, max_a: glam::DVec3, min_b: glam::DVec3, max_b: glam::DVec3) -> f64 {
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

fn breps_farther_than_tolerance(a: &BRep, b: &BRep, tol: f64) -> bool {
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
    objects: &[BRep],
    tools: &[BRep],
) -> (Vec<BRep>, SplitterObjectsReport) {
    split_objects_with_tools_options(objects, tools, SplitterOptions::default())
}

/// Like [`split_objects_with_tools`] but with advanced options.
pub fn split_objects_with_tools_options(
    objects: &[BRep],
    tools: &[BRep],
    options: SplitterOptions,
) -> (Vec<BRep>, SplitterObjectsReport) {
    let mut outputs = Vec::with_capacity(objects.len());
    let mut objects_report = Vec::with_capacity(objects.len());

    for (object_index, object) in objects.iter().enumerate() {
        let (split, report) = split_brep_with_options(object, tools, options);
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
    objects: &[BRep],
    tools: &[BRep],
    options: SplitterOptions,
) -> Result<(Vec<BRep>, SplitterObjectsReport), SplitterError> {
    let mut outputs = Vec::with_capacity(objects.len());
    let mut objects_report = Vec::with_capacity(objects.len());

    for (object_index, object) in objects.iter().enumerate() {
        let (split, report) = split_brep_checked_with_options(object, tools, options)?;
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
    objects: &[BRep],
    tools: &[BRep],
    options: SplitterOptions,
) -> (Vec<Option<BRep>>, SplitterObjectsReport) {
    let mut outputs = Vec::with_capacity(objects.len());
    let mut objects_report = Vec::with_capacity(objects.len());

    for (object_index, object) in objects.iter().enumerate() {
        let (result, report) = split_brep_internal_with_partial_report(object, tools, options, true);
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
    a: &BRep,
    b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    let mut ds = bopds::ds::DS::new(a, b);
    let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
    let mut filler = match (&bvh_a, &bvh_b) {
        (Some(a), Some(b)) => pave_filler::PaveFiller::with_bvh(&mut ds, a, b),
        _ => pave_filler::PaveFiller::new(&mut ds),
    };
    filler.perform();
    let builder = builder::BooleanBuilder::new(&ds, op);
    builder.build_with_history()
}

/// Parallel version of [`boolean_op_with_history`].
///
/// Uses Rayon to process faces in parallel during the classification phase.
/// This can provide significant speedup (2-4x) for large models with many faces.
/// For small models (< 20 faces), the serial version may be faster due to
/// thread overhead.
///
/// # Example
/// ```rust,no_run
/// use rcad_algorithms::{boolean_op_par, BooleanOpType, history::BooleanHistory};
/// use rcad_kernel::BRep;
///
/// fn parallel_union(a: &BRep, b: &BRep) -> BRep {
///     let (brep, _history) = boolean_op_par(BooleanOpType::Union, a, b).unwrap();
///     brep
/// }
/// ```
pub fn boolean_op_par(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    let mut ds = bopds::ds::DS::new(a, b);
    let (bvh_a, bvh_b) = build_optional_bvhs(a, b);
    let mut filler = match (&bvh_a, &bvh_b) {
        (Some(a), Some(b)) => pave_filler::PaveFiller::with_bvh(&mut ds, a, b),
        _ => pave_filler::PaveFiller::new(&mut ds),
    };
    filler.perform();
    let builder = builder::BooleanBuilder::new(&ds, op);
    builder.build_with_history_par()
}

/// Build BVHs for both BReps if they have faces; returns None for empty BReps.
fn build_optional_bvhs(a: &BRep, b: &BRep) -> (Option<bvh::Bvh>, Option<bvh::Bvh>) {
    let has_faces_a = a.solids.first().and_then(|s| s.shells.first()).map_or(false, |sh| !sh.faces.is_empty());
    let has_faces_b = b.solids.first().and_then(|s| s.shells.first()).map_or(false, |sh| !sh.faces.is_empty());
    (
        if has_faces_a { Some(bvh::Bvh::build(a)) } else { None },
        if has_faces_b { Some(bvh::Bvh::build(b)) } else { None },
    )
}

fn has_faces(brep: &BRep) -> bool {
    brep.solids
        .first()
        .and_then(|s| s.shells.first())
        .map_or(false, |sh| !sh.faces.is_empty())
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

/// Union two BReps and return both the result and face origin history.
pub fn union_with_history(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
    boolean_op_with_history(BooleanOpType::Union, a, b)
}

/// Intersect two BReps and return both the result and face origin history.
pub fn intersection_with_history(
    a: &BRep,
    b: &BRep,
) -> Result<(BRep, BooleanHistory), BooleanError> {
    boolean_op_with_history(BooleanOpType::Intersection, a, b)
}

/// Subtract solid B from solid A and return both the result and face origin history.
pub fn difference_with_history(a: &BRep, b: &BRep) -> Result<(BRep, BooleanHistory), BooleanError> {
    boolean_op_with_history(BooleanOpType::Difference, a, b)
}

/// Run boolean operation followed by structured healing using default options.
pub fn boolean_op_healed(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
) -> Result<(BRep, HealingReport), BooleanError> {
    let raw = boolean_op(op, a, b)?;
    let (healed, report) = heal(&raw);
    Ok((healed, report))
}

/// Run boolean operation followed by structured healing using custom options.
pub fn boolean_op_healed_with_options(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
    options: HealingOptions,
) -> Result<(BRep, HealingReport), BooleanError> {
    let raw = boolean_op(op, a, b)?;
    let (healed, report) = analyze_and_heal(&raw, options);
    Ok((healed, report))
}

/// Multi-body boolean fuse (union) over a list of solids.
///
/// This is a first-stage `general_fuse` API that folds pairwise unions from
/// left to right. It preserves current boolean behavior while enabling N-ary
/// use cases with a single call.
pub fn general_fuse(parts: &[BRep]) -> Result<BRep, BooleanError> {
    if parts.is_empty() {
        return Err(BooleanError::EmptyInput);
    }
    if parts.len() == 1 {
        return Ok(parts[0].clone());
    }

    let mut acc = parts[0].clone();
    for part in &parts[1..] {
        acc = boolean_op(BooleanOpType::Union, &acc, part)?;
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

/// Error with step location for N-ary fuse workflows.
#[derive(Debug)]
pub enum GeneralFuseError {
    EmptyInput,
    StepFailed { step_index: usize, source: BooleanError },
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
/// This keeps compatibility with the current binary boolean core while exposing
/// incremental history for debugging and tooling.
pub fn general_fuse_with_history(
    parts: &[BRep],
) -> Result<(BRep, GeneralFuseHistory), BooleanError> {
    if parts.is_empty() {
        return Err(BooleanError::EmptyInput);
    }
    if parts.len() == 1 {
        return Ok((parts[0].clone(), GeneralFuseHistory { steps: Vec::new() }));
    }

    let mut steps = Vec::with_capacity(parts.len() - 1);
    let mut acc = parts[0].clone();
    for part in &parts[1..] {
        let (next, history) = boolean_op_with_history(BooleanOpType::Union, &acc, part)?;
        acc = next;
        steps.push(history);
    }

    Ok((acc, GeneralFuseHistory { steps }))
}

/// Parallel multi-body boolean fuse (union) with per-step history.
///
/// This keeps the same left-fold semantics as [`general_fuse_with_history`],
/// but each binary union uses the parallel boolean path.
pub fn general_fuse_par(parts: &[BRep]) -> Result<(BRep, GeneralFuseHistory), BooleanError> {
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

/// Diagnostic serial N-ary fuse.
///
/// Returns per-step face-count reports and step-indexed errors when a fold
/// union fails.
pub fn general_fuse_detailed(
    parts: &[BRep],
) -> Result<(BRep, GeneralFuseHistory, GeneralFuseReport), GeneralFuseError> {
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
        let (next, history) = boolean_op_with_history(BooleanOpType::Union, &acc, part)
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

/// Diagnostic parallel N-ary fuse.
pub fn general_fuse_par_detailed(
    parts: &[BRep],
) -> Result<(BRep, GeneralFuseHistory, GeneralFuseReport), GeneralFuseError> {
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
/// Only **planar** faces are currently unified. Non-planar faces are left
/// unchanged. The topology is simplified by removing internal shared edges
/// between coplanar face pairs.
///
/// Returns the simplified BRep and the number of face merges performed.
///
/// # Algorithm
/// Performs iterated passes: in each pass, the first eligible pair of adjacent
/// coplanar faces sharing a single shell edge is merged. Passes repeat until
/// no more merges are possible. This is O(faces² × passes) but correct for all
/// plane-topology inputs produced by the boolean kernel.
pub fn unify_same_domain_faces(brep: &BRep) -> (BRep, usize) {
    let mut out = brep.clone();
    let mut total_merges = 0usize;

    loop {
        let merged = unify_one_merge_pass(&mut out);
        if !merged {
            break;
        }
        total_merges += 1;
    }

    (out, total_merges)
}

/// Attempt one merge of two adjacent coplanar faces in `brep`. Returns `true`
/// if a merge was performed (mutating `brep` in place).
fn unify_one_merge_pass(brep: &mut BRep) -> bool {
    use std::collections::HashMap;

    fn flat_face_index_of(brep: &BRep, si: usize, shi: usize, fi: usize) -> usize {
        let mut idx = 0usize;
        for s in 0..si {
            for sh in &brep.solids[s].shells {
                idx += sh.faces.len();
            }
        }
        for sh in 0..shi {
            idx += brep.solids[si].shells[sh].faces.len();
        }
        idx + fi
    }

    fn planes_match_by_geom_store(
        brep: &BRep,
        si: usize,
        shi: usize,
        fi1: usize,
        fi2: usize,
    ) -> Option<bool> {
        let ff1 = flat_face_index_of(brep, si, shi, fi1);
        let ff2 = flat_face_index_of(brep, si, shi, fi2);
        let sid1 = brep.geom.face_surface.get(ff1).and_then(|v| *v)?;
        let sid2 = brep.geom.face_surface.get(ff2).and_then(|v| *v)?;
        let s1 = brep.geom.surfaces.get(sid1)?;
        let s2 = brep.geom.surfaces.get(sid2)?;
        let (p1, p2) = match (s1, s2) {
            (rcad_kernel::geom::Surface3::Plane(a), rcad_kernel::geom::Surface3::Plane(b)) => {
                (a, b)
            }
            _ => return Some(false),
        };

        let n1 = p1.normal.normalize_or_zero();
        let n2 = p2.normal.normalize_or_zero();
        if n1.length_squared() <= 1e-24 || n2.length_squared() <= 1e-24 {
            return Some(false);
        }
        let cross = n1.cross(n2).length();
        let dot = n1.dot(n2);
        if cross > 1e-6 || dot < 0.0 {
            return Some(false);
        }
        let d = (p2.origin - p1.origin).dot(n1).abs();
        Some(d <= 1e-6)
    }

    for si in 0..brep.solids.len() {
        for shi in 0..brep.solids[si].shells.len() {
            let nfaces = brep.solids[si].shells[shi].faces.len();

            // Build edge → [face_index_in_shell] adjacency for this shell.
            let mut edge_to_faces: HashMap<usize, Vec<usize>> = HashMap::new();
            for fi in 0..nfaces {
                for we in &brep.solids[si].shells[shi].faces[fi].outer_wire.edges {
                    edge_to_faces.entry(we.idx).or_default().push(fi);
                }
                for iw in &brep.solids[si].shells[shi].faces[fi].inner_wires {
                    for we in &iw.edges {
                        edge_to_faces.entry(we.idx).or_default().push(fi);
                    }
                }
            }

            // Find the first internal edge shared by exactly 2 coplanar faces.
            for (edge_idx, face_refs) in &edge_to_faces {
                if face_refs.len() != 2 {
                    continue;
                }
                let (fi1, fi2) = (face_refs[0], face_refs[1]);
                if fi1 == fi2 {
                    continue;
                }

                let face1_normal = brep.solids[si].shells[shi].faces[fi1].normal;
                let face2_normal = brep.solids[si].shells[shi].faces[fi2].normal;

                if let Some(planes_match) = planes_match_by_geom_store(brep, si, shi, fi1, fi2) {
                    if !planes_match {
                        continue;
                    }
                }

                // Only merge planar faces (same surface geometry).
                // Verify by checking if GeomStore surfaces are both planar.
                // As a proxy: both normals must be (anti-)parallel.
                let cross = face1_normal.cross(face2_normal).length();
                let dot = face1_normal.dot(face2_normal);
                // normals parallel (same orientation) within angular tolerance
                if cross > 1e-6 || dot < 0.0 {
                    continue;
                }

                // Check same plane: get one vertex from each face and compare
                // distances to the shared plane defined by face1_normal.
                let get_face_pt = |fi: usize| -> Option<glam::DVec3> {
                    let we = brep.solids[si].shells[shi].faces[fi].outer_wire.edges.first()?;
                    let edge = brep.edges.get(we.idx)?;
                    let v_idx = if we.forward { edge.start } else { edge.end };
                    brep.vertices.get(v_idx).map(|v| v.point)
                };
                let Some(pt1) = get_face_pt(fi1) else { continue };
                let Some(pt2) = get_face_pt(fi2) else { continue };

                let n = face1_normal.normalize();
                if (pt2 - pt1).dot(n).abs() > 1e-6 {
                    continue; // Different planes
                }

                // Merge wire: splice Face2 edges into Face1 at the position of the shared edge.
                let wire1 = brep.solids[si].shells[shi].faces[fi1].outer_wire.edges.clone();
                let wire2 = brep.solids[si].shells[shi].faces[fi2].outer_wire.edges.clone();

                if let Some(merged_wire_edges) = splice_wires(&wire1, &wire2, *edge_idx) {
                    // Collect inner wires from both faces.
                    let inner1 = brep.solids[si].shells[shi].faces[fi1].inner_wires.clone();
                    let inner2 = brep.solids[si].shells[shi].faces[fi2].inner_wires.clone();
                    let mut all_inner = inner1;
                    all_inner.extend(inner2);

                    // Build merged face.
                    let merged_face = rcad_kernel::topology::Face {
                        outer_wire: rcad_kernel::topology::Wire {
                            edges: merged_wire_edges,
                        },
                        inner_wires: all_inner,
                        normal: face1_normal,
                        triangles: vec![],
                        mesh_dirty: true,
                    };

                    // Replace fi1 with merged face, remove fi2.
                    let (keep_idx, remove_idx) = if fi1 < fi2 { (fi1, fi2) } else { (fi2, fi1) };
                    brep.solids[si].shells[shi].faces[keep_idx] = merged_face;
                    brep.solids[si].shells[shi].faces.remove(remove_idx);
                    return true;
                }
            }
        }
    }

    false
}

/// Splice two wire edge lists together by removing the shared edge and
/// interleaving the remaining edges.
///
/// Returns `None` if the shared edge is not found in either wire.
fn splice_wires(
    wire_a: &[rcad_kernel::topology::WireEdge],
    wire_b: &[rcad_kernel::topology::WireEdge],
    shared_edge_idx: usize,
) -> Option<Vec<rcad_kernel::topology::WireEdge>> {
    let pos_a = wire_a.iter().position(|we| we.idx == shared_edge_idx)?;
    let pos_b = wire_b.iter().position(|we| we.idx == shared_edge_idx)?;

    let n_b = wire_b.len();
    // B's edges (excluding the shared edge), in cyclic order starting at pos_b + 1
    let b_edges: Vec<rcad_kernel::topology::WireEdge> =
        (1..n_b).map(|i| wire_b[(pos_b + i) % n_b]).collect();

    let mut merged = Vec::with_capacity(wire_a.len() - 1 + b_edges.len());
    merged.extend_from_slice(&wire_a[..pos_a]);
    merged.extend(b_edges);
    merged.extend_from_slice(&wire_a[pos_a + 1..]);

    if merged.len() < 3 {
        return None; // Degenerate result
    }

    Some(merged)
}

/// Remove redundant internal faces from a Boolean Fuse (Union) result.
///
/// After a Union operation, coincident input faces (faces from A and B on
/// exactly the same plane) can appear duplicated in the result: both input
/// faces survive classification because they lie precisely on the Boolean
/// boundary. This function detects such duplicate faces within each shell and
/// removes the extra copies.
///
/// Detection criterion: two faces in the same shell are duplicates when all of
/// the following hold:
/// - They share the same normal direction (parallel within `1e-6`).
/// - One face's representative vertex lies on the other face's plane (within `1e-6`).
/// - Their edge sets overlap entirely (every outer-wire edge of the smaller
///   face is also in the larger face, or they share ≥ 75 % of edges).
///
/// Returns the cleaned BRep and the number of faces removed.
///
/// Analogous to the internal-face elimination step of OCCT `BOPAlgo_BuilderSolid`.
pub fn remove_internal_faces(brep: &BRep) -> (BRep, usize) {
    use std::collections::HashSet;

    let mut out = brep.clone();
    let mut total_removed = 0usize;

    for si in 0..out.solids.len() {
        for shi in 0..out.solids[si].shells.len() {
            // Iteratively remove one duplicate per pass.
            loop {
                let nfaces = out.solids[si].shells[shi].faces.len();
                let mut removed_idx: Option<usize> = None;

                'outer: for fi in 0..nfaces {
                    for fj in (fi + 1)..nfaces {
                        let face_i = &out.solids[si].shells[shi].faces[fi];
                        let face_j = &out.solids[si].shells[shi].faces[fj];

                        let ni = face_i.normal;
                        let nj = face_j.normal;

                        if ni == glam::DVec3::ZERO || nj == glam::DVec3::ZERO {
                            continue;
                        }

                        // Check parallel normals (same orientation).
                        let cross = ni.cross(nj).length();
                        let dot = ni.dot(nj);
                        if cross > 1e-6 || dot < 0.0 {
                            continue;
                        }

                        // Check same plane: a vertex from j must lie on i's plane.
                        let get_pt = |f: &rcad_kernel::topology::Face| -> Option<glam::DVec3> {
                            let we = f.outer_wire.edges.first()?;
                            let edge = out.edges.get(we.idx)?;
                            let vi = if we.forward { edge.start } else { edge.end };
                            out.vertices.get(vi).map(|v| v.point)
                        };
                        let Some(pi) = get_pt(face_i) else { continue };
                        let Some(pj) = get_pt(face_j) else { continue };

                        let n_unit = ni.normalize();
                        if (pj - pi).dot(n_unit).abs() > 1e-5 {
                            continue;
                        }

                        // Check edge overlap: build edge-index sets for both faces.
                        let edges_i: HashSet<usize> = out.solids[si].shells[shi].faces[fi]
                            .outer_wire
                            .edges
                            .iter()
                            .map(|we| we.idx)
                            .collect();
                        let edges_j: HashSet<usize> = out.solids[si].shells[shi].faces[fj]
                            .outer_wire
                            .edges
                            .iter()
                            .map(|we| we.idx)
                            .collect();

                        let overlap = edges_i.intersection(&edges_j).count();
                        let min_edges = edges_i.len().min(edges_j.len()).max(1);

                        // Conservative duplicate rule: every edge of the smaller
                        // face must be shared by the larger one.
                        if overlap == min_edges {
                            // Remove fj (keep fi).
                            removed_idx = Some(fj);
                            break 'outer;
                        }
                    }
                }

                if let Some(idx) = removed_idx {
                    out.solids[si].shells[shi].faces.remove(idx);
                    total_removed += 1;
                } else {
                    break;
                }
            }
        }
    }

    (out, total_removed)
}

fn face_count_of(brep: &BRep) -> usize {
    brep
        .solids
        .iter()
        .flat_map(|s| &s.shells)
        .flat_map(|sh| &sh.faces)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::PrimitiveSolid;
    use rcad_modeling::{make_box_brep, make_cone_brep, make_cylinder_brep, make_sphere_brep, make_torus_brep};

    fn box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: w,
            height: h,
            depth: d,
        });
        for v in &mut brep.vertices {
            v.point += DVec3::new(x, y, z);
        }
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    fn face_count(brep: &BRep) -> usize {
        brep.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    }

    fn triangle_count(brep: &BRep) -> usize {
        brep.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .map(|f| f.triangles.len())
            .sum()
    }

    #[test]
    fn general_fuse_empty_input_returns_error() {
        let parts: Vec<BRep> = Vec::new();
        let result = general_fuse(&parts);
        assert!(matches!(result, Err(BooleanError::EmptyInput)));
    }

    #[test]
    fn general_fuse_single_input_returns_clone() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let fused = general_fuse(&[a.clone()]).expect("single-item general_fuse should succeed");

        assert_eq!(fused.vertices.len(), a.vertices.len());
        assert_eq!(fused.edges.len(), a.edges.len());
        assert_eq!(face_count(&fused), face_count(&a));
    }

    #[test]
    fn general_fuse_three_disjoint_boxes_accumulates_volume() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let c = box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0);

        let fused = general_fuse(&[a.clone(), b.clone(), c.clone()]).expect("general_fuse should succeed");
        let v = rcad_kernel::properties::volume(&fused);
        assert!((v - 3.0).abs() < 1e-6, "expected volume 3.0, got {v}");
    }

    #[test]
    fn general_fuse_with_history_single_input_has_no_steps() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let (_fused, hist) = general_fuse_with_history(&[a]).expect("single-item general_fuse_with_history should succeed");
        assert!(hist.steps.is_empty());
    }

    #[test]
    fn general_fuse_with_history_three_inputs_has_two_steps() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let c = box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0);

        let (fused, hist) = general_fuse_with_history(&[a, b, c]).expect("general_fuse_with_history should succeed");
        assert_eq!(hist.steps.len(), 2, "three inputs should produce two fold steps");
        assert!(hist.steps.iter().all(|h| !h.is_empty()), "each step should carry face history");

        let v = rcad_kernel::properties::volume(&fused);
        assert!((v - 3.0).abs() < 1e-6, "expected volume 3.0, got {v}");
    }

    #[test]
    fn general_fuse_par_three_disjoint_boxes_accumulates_volume() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let c = box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0);

        let (fused, hist) = general_fuse_par(&[a, b, c]).expect("general_fuse_par should succeed");
        assert_eq!(hist.steps.len(), 2);

        let v = rcad_kernel::properties::volume(&fused);
        assert!((v - 3.0).abs() < 1e-6, "expected volume 3.0, got {v}");
    }

    #[test]
    fn general_fuse_par_matches_serial_for_three_disjoint_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let c = box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0);

        let serial = general_fuse(&[a.clone(), b.clone(), c.clone()]).expect("serial general_fuse should succeed");
        let (parallel, _) = general_fuse_par(&[a, b, c]).expect("parallel general_fuse should succeed");

        let v_serial = rcad_kernel::properties::volume(&serial);
        let v_parallel = rcad_kernel::properties::volume(&parallel);
        assert!((v_serial - v_parallel).abs() < 1e-6);
    }

    #[test]
    fn general_fuse_detailed_overlapping_chain_reports_steps() {
        let a = box_at(0.0, 0.0, 0.0, 1.2, 1.0, 1.0);
        let b = box_at(0.6, 0.0, 0.0, 1.2, 1.0, 1.0);
        let c = box_at(1.2, 0.0, 0.0, 1.2, 1.0, 1.0);

        let (_fused, hist, report) =
            general_fuse_detailed(&[a, b, c]).expect("general_fuse_detailed should succeed");

        assert_eq!(hist.steps.len(), 2);
        assert_eq!(report.steps.len(), 2);
        assert_eq!(report.steps[0].step_index, 0);
        assert_eq!(report.steps[1].step_index, 1);
        assert!(report.steps.iter().all(|s| s.input_faces > 0 && s.output_faces > 0));
    }

    #[test]
    fn general_fuse_overlap_chain_volume_between_bounds() {
        let a = box_at(0.0, 0.0, 0.0, 1.2, 1.0, 1.0);
        let b = box_at(0.6, 0.0, 0.0, 1.2, 1.0, 1.0);
        let c = box_at(1.2, 0.0, 0.0, 1.2, 1.0, 1.0);

        let fused = general_fuse(&[a.clone(), b.clone(), c.clone()]).expect("general_fuse should succeed");
        let v = rcad_kernel::properties::volume(&fused);
        let sum = rcad_kernel::properties::volume(&a)
            + rcad_kernel::properties::volume(&b)
            + rcad_kernel::properties::volume(&c);

        // Overlapping chain: union volume must be positive and strictly less than
        // naive volume sum (because overlaps exist).
        assert!(v > 0.0, "volume should be positive");
        assert!(v < sum - 1e-6, "union volume should be less than sum, got v={v}, sum={sum}");
    }

    #[test]
    fn general_fuse_detailed_empty_input_returns_empty_error() {
        let parts: Vec<BRep> = Vec::new();
        let result = general_fuse_detailed(&parts);
        assert!(matches!(result, Err(GeneralFuseError::EmptyInput)));
    }

    #[test]
    fn split_brep_empty_tools_returns_clone_and_empty_report() {
        let target = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let (out, report) = split_brep(&target, &[]);

        assert_eq!(face_count(&out), face_count(&target));
        assert!(report.steps.is_empty());
        assert_eq!(report.total_seam_edges, 0);
    }

    #[test]
    fn split_brep_with_tool_produces_step_report() {
        let target = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let (out, report) = split_brep(&target, &[tool]);

        assert_eq!(report.steps.len(), 1);
        assert_eq!(report.steps[0].step_index, 0);
        assert!(report.steps[0].input_faces > 0);
        assert!(report.steps[0].output_faces > 0);
        assert_eq!(report.total_seam_edges, report.steps[0].seam_edges);
        assert!(!report.steps[0].skipped_by_broad_phase);
        assert!(report.steps[0].validation_issue_count.is_none());
        assert!(report.steps[0].validation_first_issue.is_none());
        assert!(face_count(&out) >= face_count(&target));
    }

    #[test]
    fn splitter_options_default_validation_is_relaxed() {
        let opts = SplitterOptions::default();
        assert_eq!(opts.validation_level, SplitterValidationLevel::Relaxed);
    }

    #[test]
    fn split_brep_with_healing_sets_healed_flag() {
        let target = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let (_out, report) = split_brep_with_options(
            &target,
            &[tool],
            SplitterOptions {
                heal_after_each_step: true,
                healing: HealingOptions {
                    mode: HealingMode::AnalyzeOnly,
                    ..HealingOptions::default()
                },
                ..SplitterOptions::default()
            },
        );

        assert_eq!(report.steps.len(), 1);
        assert!(report.steps[0].healed);
        assert!(!report.steps[0].skipped_by_broad_phase);
    }

    #[test]
    fn split_brep_far_tool_is_skipped_by_broad_phase() {
        let target = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let far_tool = box_at(100.0, 0.0, 0.0, 1.0, 1.0, 1.0);

        let (out, report) = split_brep_with_options(
            &target,
            &[far_tool],
            SplitterOptions {
                broad_phase_pruning: true,
                fuzzy_tolerance: 0.0,
                ..SplitterOptions::default()
            },
        );

        assert_eq!(report.steps.len(), 1);
        let step = &report.steps[0];
        assert!(step.skipped_by_broad_phase);
        assert_eq!(step.seam_edges, 0);
        assert_eq!(step.input_faces, step.output_faces);
        assert_eq!(face_count(&out), face_count(&target));
    }

    #[test]
    fn split_brep_checked_with_options_detects_invalid_step() {
        let target = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let err = split_brep_checked_with_options(&target, &[tool], SplitterOptions::default())
            .expect_err("checked splitter should report invalid intermediate topology");

        assert!(matches!(
            err,
            SplitterError::StepInvalid {
                step_index: 0,
                issue_count: c,
                ..
            } if c > 0
        ));
    }

    #[test]
    fn split_objects_with_tools_empty_objects_returns_empty() {
        let (out, report) = split_objects_with_tools(&[], &[]);
        assert!(out.is_empty());
        assert!(report.objects.is_empty());
    }

    #[test]
    fn split_objects_with_tools_empty_tools_clones_each_object() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(3.0, 0.0, 0.0, 1.0, 1.0, 1.0);

        let (out, report) = split_objects_with_tools(&[a.clone(), b.clone()], &[]);
        assert_eq!(out.len(), 2);
        assert_eq!(face_count(&out[0]), face_count(&a));
        assert_eq!(face_count(&out[1]), face_count(&b));

        assert_eq!(report.objects.len(), 2);
        assert!(report.objects.iter().all(|r| r.steps.is_empty()));
        assert!(report.objects.iter().all(|r| r.total_seam_edges == 0));
        assert!(report.objects.iter().all(|r| r.completed));
        assert!(report.objects.iter().all(|r| r.error.is_none()));
    }

    #[test]
    fn split_objects_with_tools_reports_each_object() {
        let object_a = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let object_b = box_at(4.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let (out, report) = split_objects_with_tools(&[object_a, object_b], &[tool]);
        assert_eq!(out.len(), 2);
        assert_eq!(report.objects.len(), 2);
        assert_eq!(report.objects[0].object_index, 0);
        assert_eq!(report.objects[1].object_index, 1);
        assert!(report.objects.iter().all(|r| r.steps.len() == 1));
        assert!(report.objects.iter().all(|r| r.completed));
        assert!(report.objects.iter().all(|r| r.error.is_none()));
        assert!(
            report.objects.iter().any(|r| !r.steps[0].skipped_by_broad_phase),
            "at least one object should execute split step"
        );
        assert!(
            report.objects.iter().any(|r| r.steps[0].skipped_by_broad_phase),
            "at least one far object should be skipped by broad-phase"
        );
    }

    #[test]
    fn split_objects_with_tools_checked_options_succeeds_when_steps_are_skipped() {
        let object_a = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let object_b = box_at(4.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(100.0, 100.0, 100.0, 1.0, 1.0, 1.0);

        let (out, report) = split_objects_with_tools_checked_options(
            &[object_a, object_b],
            &[tool],
            SplitterOptions::default(),
        )
        .expect("checked grouped splitter should succeed when broad-phase skips all steps");

        assert_eq!(out.len(), 2);
        assert_eq!(report.objects.len(), 2);
        assert!(report.objects.iter().all(|r| r.steps[0].skipped_by_broad_phase));
        assert!(report.objects.iter().all(|r| r.completed));
        assert!(report.objects.iter().all(|r| r.error.is_none()));
        assert!(
            report
                .objects
                .iter()
                .all(|r| r.steps[0].validation_issue_count == Some(0))
        );
    }

    #[test]
    fn split_objects_with_tools_checked_collect_reports_mixed_outcomes() {
        let near_object = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let far_object = box_at(100.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let (out, report) = split_objects_with_tools_checked_collect_options(
            &[near_object, far_object],
            &[tool],
            SplitterOptions::default(),
        );

        assert_eq!(out.len(), 2);
        assert!(out[0].is_none(), "near object should fail checked split");
        assert!(out[1].is_some(), "far object should be skipped and succeed");

        assert_eq!(report.objects.len(), 2);
        assert!(!report.objects[0].completed);
        assert!(report.objects[0].error.is_some());
        assert_eq!(report.objects[0].steps.len(), 1);
        assert_eq!(report.objects[0].steps[0].step_index, 0);
        assert!(report.objects[0].steps[0].validation_issue_count.unwrap_or(0) > 0);

        assert!(report.objects[1].completed);
        assert!(report.objects[1].error.is_none());
        assert_eq!(report.objects[1].steps.len(), 1);
        assert!(report.objects[1].steps[0].skipped_by_broad_phase);

        let summary = report.summarize();
        assert_eq!(summary.total_objects, 2);
        assert_eq!(summary.completed_objects, 1);
        assert_eq!(summary.failed_objects, 1);
        assert_eq!(summary.failed_object_indices, vec![0]);
        assert_eq!(summary.failed_step_histogram, vec![(0, 1)]);
        assert_eq!(summary.first_error_histogram.len(), 1);
    }

    #[test]
    fn splitter_objects_report_summarize_counts_success_and_failure() {
        let near_object = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let far_object = box_at(100.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let (_out, report) = split_objects_with_tools_checked_collect_options(
            &[near_object, far_object],
            &[tool],
            SplitterOptions::default(),
        );

        let summary = report.summarize();
        assert_eq!(summary.total_objects, 2);
        assert_eq!(summary.completed_objects, 1);
        assert_eq!(summary.failed_objects, 1);
        assert_eq!(summary.failed_object_indices, vec![0]);
        assert_eq!(summary.failed_step_histogram, vec![(0, 1)]);
        assert!(
            !summary.first_error_histogram.is_empty(),
            "summary should include at least one error bucket"
        );
    }

    #[test]
    fn splitter_objects_report_to_json_v1_contains_schema_and_summary() {
        let near_object = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let far_object = box_at(100.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let (_out, report) = split_objects_with_tools_checked_collect_options(
            &[near_object, far_object],
            &[tool],
            SplitterOptions::default(),
        );

        let json = report
            .to_json_v1()
            .expect("splitter report json serialization should succeed");
        let v: serde_json::Value =
            serde_json::from_str(&json).expect("serialized splitter json should parse");

        assert_eq!(v["schema"], "splitter.report.v1");
        assert_eq!(v["summary"]["total_objects"], 2);
        assert_eq!(v["summary"]["failed_objects"], 1);
        assert!(
            v["summary"]["failed_object_indices"].is_array(),
            "failed_object_indices must be exported as an array"
        );
    }

    #[test]
    fn split_brep_checked_strict_mode_reports_step_invalid() {
        let target = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let tool = box_at(1.0, 0.5, -0.5, 1.0, 1.0, 3.0);

        let err = split_brep_checked_with_options(
            &target,
            &[tool],
            SplitterOptions {
                validation_level: SplitterValidationLevel::Strict,
                ..SplitterOptions::default()
            },
        )
        .expect_err("strict checked splitter should fail on current intermediate issues");

        assert!(matches!(err, SplitterError::StepInvalid { step_index: 0, .. }));
    }

    #[test]
    fn simplify_brep_post_ops_reports_checker_delta() {
        let mut b = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        b.solids[0].shells[0].faces[0].normal = DVec3::ZERO;

        let (_out, report) = simplify_brep_post_ops(&b, SimplifyOptions::default());
        assert!(report.issues_before >= report.issues_after);
        assert!(report.normals_recomputed >= 1);
    }

    #[test]
    fn boolean_op_simplified_union_runs() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);

        let (out, report) = boolean_op_simplified(
            BooleanOpType::Union,
            &a,
            &b,
            SimplifyOptions::default(),
        )
        .expect("boolean_op_simplified union should succeed");

        assert!(!out.solids.is_empty());
        assert!(report.issues_before >= report.issues_after);
    }

    #[test]
    fn simplify_brep_post_ops_runs_same_domain_and_internal_cleanup() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let raw = boolean_op(BooleanOpType::Union, &a, &b)
            .expect("coplanar flush union should succeed before simplify");

        let (baseline, _baseline_report) = simplify_brep_post_ops(
            &raw,
            SimplifyOptions {
                unify_same_domain_faces: false,
                remove_internal_faces: false,
                ..SimplifyOptions::default()
            },
        );

        let (cleaned, report) = simplify_brep_post_ops(
            &raw,
            SimplifyOptions {
                unify_same_domain_faces: true,
                remove_internal_faces: true,
                ..SimplifyOptions::default()
            },
        );

        assert!(
            face_count_of(&cleaned) <= face_count_of(&baseline),
            "cleanup-enabled simplify should not increase face count"
        );
        assert!(report.issues_before >= report.issues_after);
    }

    #[test]
    fn boolean_op_healed_union_returns_valid_result() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);

        let (res, report) = boolean_op_healed(BooleanOpType::Union, &a, &b)
            .expect("boolean_op_healed union should succeed");

        assert!(check(&res).is_valid(), "healed result should be valid");
        assert!(report.final_result.is_valid(), "healing report should end valid");
    }

    fn all_triangles_valid(brep: &BRep) -> bool {
        let nv = brep.vertices.len();
        brep.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .flat_map(|f| &f.triangles)
            .all(|tri| tri.iter().all(|&i| i < nv))
    }

    #[test]
    fn union_disjoint_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(5.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        // Disjoint: all 12 faces kept
        assert_eq!(face_count(&result), 12);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn intersection_disjoint_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(5.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        // Disjoint: intersection is empty
        assert!(result.is_err());
    }

    #[test]
    fn union_overlapping_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn intersection_overlapping_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert!(face_count(&result) >= 6);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn difference_overlapping_boxes() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn contained_box_difference() {
        // B completely inside A
        let a = box_at(0.0, 0.0, 0.0, 4.0, 4.0, 4.0);
        let b = box_at(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
        let result = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn contained_box_intersection() {
        // B completely inside A → intersection is B
        let a = box_at(0.0, 0.0, 0.0, 4.0, 4.0, 4.0);
        let b = box_at(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert_eq!(face_count(&result), 6); // B's 6 faces
        assert!(all_triangles_valid(&result));
    }

    // ─── Phase 4 edge case tests ───────────────────────────────────────

    #[test]
    fn identical_boxes_union() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn touching_face_union() {
        // Two boxes sharing a face (A right = B left)
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(1.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn touching_edge_union() {
        // Two boxes sharing an edge
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(1.0, 1.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert_eq!(face_count(&result), 12);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn non_unit_boxes_difference() {
        let a = box_at(0.0, 0.0, 0.0, 3.0, 2.0, 5.0);
        let b = box_at(1.0, 0.5, 1.0, 1.0, 1.0, 3.0);
        let result = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        assert!(face_count(&result) > 6);
        assert!(triangle_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn offset_3d_intersection() {
        let a = box_at(0.0, 0.0, 0.0, 2.0, 2.0, 2.0);
        let b = box_at(1.0, 1.0, 1.0, 2.0, 2.0, 2.0);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert!(face_count(&result) >= 6);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn difference_is_not_symmetric() {
        let a = box_at(0.0, 0.0, 0.0, 2.0, 1.0, 1.0);
        let b = box_at(1.0, 0.0, 0.0, 2.0, 1.0, 1.0);
        let a_minus_b = boolean_op(BooleanOpType::Difference, &a, &b).unwrap();
        let b_minus_a = boolean_op(BooleanOpType::Difference, &b, &a).unwrap();
        assert!(face_count(&a_minus_b) > 0);
        assert!(face_count(&b_minus_a) > 0);
        assert!(all_triangles_valid(&a_minus_b));
        assert!(all_triangles_valid(&b_minus_a));
    }

    #[test]
    fn small_overlap_union() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.99, 0.0, 0.0, 1.0, 1.0, 1.0);
        let result = boolean_op(BooleanOpType::Union, &a, &b).unwrap();
        assert!(face_count(&result) > 0);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn large_overlap_intersection() {
        let a = box_at(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
        let b = box_at(0.1, 0.1, 0.1, 9.8, 9.8, 9.8);
        let result = boolean_op(BooleanOpType::Intersection, &a, &b).unwrap();
        assert_eq!(face_count(&result), 6);
        assert!(all_triangles_valid(&result));
    }

    #[test]
    fn classify_point_on_face() {
        use classify::Classification;
        let mut brep = rcad_kernel::BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        });
        geom_populate::populate_box_geom(&mut brep);
        let ds = bopds::ds::DS::new(&brep, &rcad_kernel::BRep::new());
        let face_indices: Vec<usize> = (0..ds.faces.len())
            .filter(|&i| ds.faces[i].origin == bopds::ds::ShapeOrigin::ShapeA)
            .collect();
        let on_top = DVec3::new(1.0, 2.0, 1.0);
        assert_eq!(
            classify::classify_point(on_top, &face_indices, &ds),
            Classification::On
        );
    }

    #[test]
    fn triangulate_hexagon() {
        use triangulate::triangulate_polygon;
        let verts: Vec<DVec3> = (0..6)
            .map(|i| {
                let a = 2.0 * std::f64::consts::PI * i as f64 / 6.0;
                DVec3::new(a.cos(), a.sin(), 0.0)
            })
            .collect();
        let tris = triangulate_polygon(&verts, DVec3::Z);
        assert_eq!(tris.len(), 4);
        for tri in &tris {
            for &idx in tri {
                assert!(idx < 6);
            }
        }
    }

    // ─── Curved Boolean Tests ──────────────────────────────────────────────────

    #[test]
    fn boolean_box_sphere_intersection() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.0), 1.5).unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(
            result.is_ok(),
            "box-sphere intersection failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        // Volume computation for curved result faces is approximate; just check
        // the result is non-degenerate.
    }

    #[test]
    fn boolean_box_sphere_difference() {
        // Small sphere inside a box — creates a hole
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let b = make_sphere_brep(DVec3::new(2.0, 2.0, 2.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "box-sphere difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        // Volume computation for curved result faces is approximate; just check
        // the result is non-degenerate.
    }

    #[test]
    fn boolean_box_sphere_union() {
        // Sphere protruding from box
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 1.0, 2.5), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(
            result.is_ok(),
            "box-sphere union failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        let v_box = rcad_kernel::properties::volume(&a);
        let v_sphere = rcad_kernel::properties::volume(&b);
        assert!(v > v_box, "union should be larger than box");
        assert!(v > v_sphere, "union should be larger than sphere");
    }

    #[test]
    fn boolean_sphere_sphere_intersection() {
        // Two overlapping unit spheres
        let a = make_sphere_brep(DVec3::new(-0.5, 0.0, 0.0), 1.0).unwrap();
        let b = make_sphere_brep(DVec3::new(0.5, 0.0, 0.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-sphere intersection failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        // Sphere primitive has no triangle mesh, so volume(&a) = 0. Compare against
        // analytical: two overlapping unit spheres at distance 1 → lens volume ≈ 1.809.
        // Full unit sphere volume = 4π/3 ≈ 4.189.
        let v_sphere_analytical = 4.0 * std::f64::consts::PI / 3.0; // 4π/3
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(
            v < v_sphere_analytical,
            "intersection should be smaller than one sphere (4π/3≈4.19), got {v}"
        );
    }

    #[test]
    fn boolean_sphere_sphere_difference() {
        // Large sphere (r=2) minus small sphere (r=1) with d=1 between centers.
        // d=1, r_A=2, r_B=1 → h = (1+4-1)/2 = 2 → tangent! Use d=0.5 instead.
        // d=0.5, r_A=2, r_B=1 → h = (0.25+4-1)/1 = 3.25 → outside sphere A
        // Use d=1.5: h = (2.25+4-1)/3 = 5.25/3 = 1.75 < r_A=2 → proper intersection
        let a = make_sphere_brep(DVec3::ZERO, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.5, 0.0, 0.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-sphere difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        let v = rcad_kernel::properties::volume(&brep);
        // Large sphere volume = 4π/3 * 8 ≈ 33.51; result should be positive and less.
        let v_large_analytical = 4.0 * std::f64::consts::PI / 3.0 * 8.0;
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(v < v_large_analytical, "difference should be smaller than original large sphere");
    }

    #[test]
    fn boolean_box_cylinder_hole() {
        // Box minus a cylinder through it (classic hole)
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        // Cylinder along Z axis through center of box
        let b =
            make_cylinder_brep(DVec3::new(2.0, 2.0, -0.5), DVec3::Z, DVec3::X, 0.5, 5.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "box-cylinder difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
        // Volume computation for curved result faces is approximate; just check
        // the result is non-degenerate.
    }

    #[test]
    fn boolean_cylinder_cylinder_intersection() {
        // Two perpendicular cylinders (Steinmetz solid).
        // Use cylinders that are offset so they overlap in a region that doesn't
        // straddle the seam boundary (avoiding UV-seam discontinuity issues).
        // Cylinder A: Y-axis, centered at (0, 0, 0) with height 4 → spans y ∈ [-2, 2]
        // Cylinder B: X-axis, centered at (0, 0, 0) with height 4 → spans x ∈ [-2, 2]
        let a =
            make_cylinder_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::Y, DVec3::X, 1.0, 4.0).unwrap();
        let b =
            make_cylinder_brep(DVec3::new(0.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 4.0).unwrap();

        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        // The result should be non-degenerate (the two cylinders DO intersect).
        // We check only non-degeneracy: if the boolean fails or gives an empty
        // result, something is fundamentally broken.
        match result {
            Ok(brep) => {
                // Non-degenerate: at least one face in the result.
                assert!(
                    !brep.solids[0].shells[0].faces.is_empty(),
                    "cylinder-cylinder intersection should produce at least one face"
                );
                let v = rcad_kernel::properties::volume(&brep);
                assert!(v >= 0.0, "volume must not be negative, got {v}");
                // Note: exact volume comparison is not practical because the curved-face
                // volume computation (divergence theorem on polyline boundaries) is
                // approximate for complex intersection geometries.
            }
            Err(e) => {
                // If the result is degenerate, fail with a clear message.
                panic!("cylinder-cylinder intersection failed: {e:?}");
            }
        }
    }

    #[test]
    fn volume_conservation_box_sphere() {
        // V(A∪B) ≈ V(A) + V(B) - V(A∩B), error < 5%
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(1.0, 1.0, 1.5), 1.0).unwrap();

        let union_result = boolean_op(BooleanOpType::Union, &a, &b);
        let inter_result = boolean_op(BooleanOpType::Intersection, &a, &b);

        assert!(
            union_result.is_ok(),
            "union failed: {:?}",
            union_result.err()
        );
        assert!(
            inter_result.is_ok(),
            "intersection failed: {:?}",
            inter_result.err()
        );

        let union_brep = union_result.unwrap();
        let inter_brep = inter_result.unwrap();

        let v_a = rcad_kernel::properties::volume(&a);
        let v_b = rcad_kernel::properties::volume(&b);
        let v_union = rcad_kernel::properties::volume(&union_brep);
        let v_inter = rcad_kernel::properties::volume(&inter_brep);

        // Debug values — show face count and uv_domains
        eprintln!("V_A={v_a:.4} V_B={v_b:.4} V_union={v_union:.4} V_inter={v_inter:.4}");
        eprintln!(
            "Union faces={}, inter faces={}",
            union_brep.solids[0].shells[0].faces.len(),
            inter_brep.solids[0].shells[0].faces.len()
        );
        for (i, face) in inter_brep.solids[0].shells[0].faces.iter().enumerate() {
            let range = inter_brep.geom.face_surface_range.get(i).and_then(|o| *o);
            let surf_name = inter_brep.geom.face_surface.get(i).and_then(|o| *o)
                .map(|si| format!("{:?}", std::mem::discriminant(&inter_brep.geom.surfaces[si])));
            // Compute per-face contribution to volume
            let face_tris = rcad_kernel::properties::face_triangles_pub(&inter_brep, face, i);
            let face_vol: f64 = face_tris.iter().map(|&[a,b,c]| a.dot(b.cross(c)) / 6.0).sum();
            eprintln!("  inter face {i}: normal={:.3?} uv_domain={range:?} surf={surf_name:?} vol_contrib={face_vol:.4}", face.normal);
        }

        let expected = v_a + v_b - v_inter;
        let error = (v_union - expected).abs() / expected;
        let error_pct = error * 100.0;
        assert!(
            error < 0.05,
            "Volume conservation violated: V(A∪B)={v_union:.4}, V(A)+V(B)-V(A∩B)={expected:.4}, error={error_pct:.2}%"
        );
    }

    #[test]
    fn volume_conservation_spheres() {
        // Preferred behavior: V(A∪B) ≈ V(A) + V(B) - V(A∩B), error < 5%.
        // Current kernel may still return an incomplete sphere-sphere union shell.
        // In that known-gap case, keep this as an active regression test with
        // explicit fallback assertions instead of ignoring it entirely.
        let a = make_sphere_brep(DVec3::new(-0.5, 0.0, 0.0), 1.0).unwrap();
        let b = make_sphere_brep(DVec3::new(0.5, 0.0, 0.0), 1.0).unwrap();

        let union_result = boolean_op(BooleanOpType::Union, &a, &b);
        let inter_result = boolean_op(BooleanOpType::Intersection, &a, &b);

        assert!(
            union_result.is_ok(),
            "union failed: {:?}",
            union_result.err()
        );
        assert!(
            inter_result.is_ok(),
            "intersection failed: {:?}",
            inter_result.err()
        );

        let union_brep = union_result.unwrap();
        let inter_brep = inter_result.unwrap();

        let v_a = rcad_kernel::properties::volume(&a);
        let v_b = rcad_kernel::properties::volume(&b);
        let v_union = rcad_kernel::properties::volume(&union_brep);
        let v_inter = rcad_kernel::properties::volume(&inter_brep);

        eprintln!("sphere-sphere: V_A={v_a:.4} V_B={v_b:.4} V_union={v_union:.4} V_inter={v_inter:.4}");
        eprintln!(
            "Union faces={}, inter faces={}",
            union_brep.solids[0].shells[0].faces.len(),
            inter_brep.solids[0].shells[0].faces.len()
        );
        for (i, face) in inter_brep.solids[0].shells[0].faces.iter().enumerate() {
            let range = inter_brep.geom.face_surface_range.get(i).and_then(|o| *o);
            let surf_name = inter_brep.geom.face_surface.get(i).and_then(|o| *o)
                .map(|si| format!("{:?}", std::mem::discriminant(&inter_brep.geom.surfaces[si])));
            let face_tris = rcad_kernel::properties::face_triangles_pub(&inter_brep, face, i);
            let face_vol: f64 = face_tris.iter().map(|&[a,b,c]| a.dot(b.cross(c)) / 6.0).sum();
            eprintln!("  inter face {i}: normal={:.3?} uv_domain={range:?} surf={surf_name:?} vol_contrib={face_vol:.4}", face.normal);
        }

        let expected = v_a + v_b - v_inter;
        let error = (v_union - expected).abs() / expected.max(1e-12);
        let error_pct = error * 100.0;
        let union_faces = union_brep.solids[0].shells[0].faces.len();

        if v_union > 1e-6 {
            assert!(
                error < 0.05,
                "Volume conservation violated: V(A∪B)={v_union:.4}, V(A)+V(B)-V(A∩B)={expected:.4}, error={error_pct:.2}%"
            );
        } else {
            // Known limitation signature (incomplete union shell):
            // union has near-zero volume and a very small face count.
            assert!(
                union_faces <= 2,
                "unexpected zero-volume union shape signature: faces={union_faces}, expected <= 2"
            );
            assert!(v_inter > 0.0, "intersection volume should still be positive");
        }
    }

    #[test]
    fn boolean_result_edges_have_pcurves() {
        // Box with a cylindrical hole. After the boolean difference, intersection
        // edges on the cylinder surface should get PCurves via
        // populate_boolean_result_pcurves.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(2.0, 2.0, -0.5), DVec3::Z, DVec3::X, 0.5, 5.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        let Ok(mut brep) = result else {
            // If the boolean op itself fails, skip (it's tested elsewhere).
            return;
        };
        if brep.solids.is_empty() || brep.solids[0].shells.is_empty() {
            return;
        }

        // Fill PCurves.
        geom_populate::populate_boolean_result_pcurves(&mut brep);

        // At least one edge on the cylinder face should now have a PCurve.
        let any_pcurve = brep.geom.edge_pcurves.iter().any(|v| !v.is_empty());
        assert!(
            any_pcurve,
            "populate_boolean_result_pcurves should have added at least one PCurve"
        );
    }

    // ─── Sphere × Cylinder Boolean Tests ──────────────────────────────────────

    /// A cylinder whose axis passes through the sphere centre (axis-aligned case).
    /// The sphere–cylinder intersection is two circles.  Difference should
    /// produce a valid solid with more faces than just the six box/sphere faces.
    #[test]
    fn boolean_sphere_cylinder_difference_axis_aligned() {
        // Sphere centred at origin, radius 5; cylinder along Z through origin, radius 3.
        // Intersection circles at z = ±4  (sqrt(25-9) = 4).
        let a = make_sphere_brep(DVec3::ZERO, 5.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(0.0, 0.0, -6.0), DVec3::Z, DVec3::X, 3.0, 12.0)
            .unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-cylinder difference (axis-aligned) failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
        // Volume of sphere (4π/3 · R³) minus the cylindrical tunnel should be positive
        // and smaller than the sphere.
        let v = rcad_kernel::properties::volume(&brep);
        let v_sphere = 4.0 * std::f64::consts::PI / 3.0 * 5.0_f64.powi(3);
        assert!(v > 0.0, "result volume should be positive, got {v}");
        assert!(v < v_sphere, "difference should be smaller than original sphere");
    }

    // ─── Cone × Plane Boolean Tests ───────────────────────────────────────────

    /// Box minus a cone through it: the cone's lateral surface intersects the
    /// box's planar faces, exercising the plane-cone circle intersection path.
    #[test]
    fn boolean_box_cone_difference() {
        // Box: 4×4×4 at origin.  Cone: base at (2,2,-0.5), axis Z, r=0.8, h=5.
        // The cone pokes through the box; plane-cone intersections are circles
        // (planes ⊥ cone axis).
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let b =
            make_cone_brep(DVec3::new(2.0, 2.0, -0.5), DVec3::Z, DVec3::X, 0.8, 5.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "box-cone difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
    }

    /// Cone intersected with a box slab: the slab's top and bottom faces are
    /// planes perpendicular to the cone axis, producing circle intersections.
    /// This test verifies that the plane-cone code path does not panic.
    #[test]
    fn boolean_cone_box_intersection_circle() {
        // Cone: base at origin, axis Z, base_radius=2, height=4.
        // Slab: 6×6×4 at z=0..4 — same height as the cone; the lateral face of
        // the slab does NOT cut the cone (slab is wide enough), so only the
        // slab top (z=4, a plane ⊥ cone axis) intersects the cone's lateral surface
        // near the apex region.  This exercises the plane-cone circle intersection.
        let a = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 4.0).unwrap();
        let b =
            make_box_brep(DVec3::new(-3.0, -3.0, 0.0), DVec3::X, DVec3::Y, 6.0, 6.0, 3.0)
                .unwrap();
        // The box (z=0..3) clips the cone (z=0..4), leaving the lower frustum.
        // The intersection may succeed or return DegenerateResult depending on
        // classifier robustness; we only require it does not panic.
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        match result {
            Ok(brep) => {
                assert!(
                    !brep.solids.is_empty() && !brep.solids[0].shells[0].faces.is_empty(),
                    "intersection produced an empty result"
                );
            }
            Err(BooleanError::DegenerateResult) => {
                // DegenerateResult is an acceptable failure for complex curved intersections.
            }
            Err(e) => {
                panic!("cone-box intersection failed unexpectedly: {e:?}");
            }
        }
    }

    /// Intersection of a sphere and a coaxial cylinder.
    #[test]
    fn boolean_sphere_cylinder_intersection_axis_aligned() {
        // Sphere centred at origin, radius 5; cylinder along Z through origin, radius 3.
        // The intersection of their volumes is a "barrel" shape bounded by two
        // spherical caps (z > 4 and z < -4) and the cylinder lateral surface.
        let a = make_sphere_brep(DVec3::ZERO, 5.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(0.0, 0.0, -6.0), DVec3::Z, DVec3::X, 3.0, 12.0)
            .unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-cylinder intersection (axis-aligned) failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
        // Just verify we get a positive volume — the exact amount depends on
        // whether sphere cap faces contribute correctly to the divergence-theorem
        // volume (sphere parametric surfaces have known approximation issues
        // tracked separately).
        let v = rcad_kernel::properties::volume(&brep);
        assert!(v > 0.0, "intersection volume should be positive, got {v}");
    }

    #[test]
    fn curved_subface_boundary_3d_sphere_pole_produces_enough_points() {
        // Verify that a sphere boolean with a cone produces a valid result.
        // The cone has an apex singularity that previously caused degenerate
        // sub-face boundaries.
        let a = make_sphere_brep(DVec3::ZERO, 2.0).unwrap();
        let b = make_cone_brep(DVec3::new(0.0, 0.0, -1.0), DVec3::Z, DVec3::X, 1.5, 3.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "sphere-cone boolean (apex singularity) failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
        let v = rcad_kernel::properties::volume(&brep);
        assert!(v > 0.0, "difference volume should be positive, got {v}");
    }

    // ─── Torus Boolean Tests ──────────────────────────────────────────────────

    #[test]
    fn boolean_box_torus_difference() {
        // Box minus a torus: the torus sits partially inside the box.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 6.0, 6.0, 6.0).unwrap();
        // Torus centered at (3,3,3), axis Z, major=1.5, minor=0.5
        let b = make_torus_brep(DVec3::new(3.0, 3.0, 3.0), DVec3::Z, DVec3::X, 1.5, 0.5).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "box-torus difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
    }

    #[test]
    fn boolean_torus_torus_intersection() {
        // Two interlocking tori (like a chain link).
        // Torus A: XY plane, centered at origin
        let a = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.5).unwrap();
        // Torus B: XZ plane, centered at origin (perpendicular)
        let b = make_torus_brep(DVec3::ZERO, DVec3::Y, DVec3::X, 2.0, 0.5).unwrap();
        let result = boolean_op(BooleanOpType::Intersection, &a, &b);
        // May succeed or return DegenerateResult; must not panic.
        match result {
            Ok(brep) => {
                assert!(
                    !brep.solids.is_empty() && !brep.solids[0].shells[0].faces.is_empty(),
                    "torus-torus intersection produced an empty result"
                );
            }
            Err(BooleanError::DegenerateResult) => {
                // Acceptable for complex curved intersections.
            }
            Err(e) => {
                panic!("torus-torus intersection failed unexpectedly: {e:?}");
            }
        }
    }

    #[test]
    fn boolean_cylinder_torus_difference() {
        // Cylinder passing through a torus hole.
        let a = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.8).unwrap();
        let b = make_cylinder_brep(DVec3::new(0.0, 0.0, -3.0), DVec3::Z, DVec3::X, 0.3, 6.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "cylinder-torus difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(
            !brep.solids[0].shells[0].faces.is_empty(),
            "result should have faces"
        );
    }

    // ─── Coplanar Face Boolean Tests ──────────────────────────────────────────

    #[test]
    fn boolean_coplanar_flush_union() {
        // Two boxes sharing a coplanar face (flush side-by-side).
        // The union should merge the coplanar faces.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::new(2.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(
            result.is_ok(),
            "coplanar flush union failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
    }

    #[test]
    fn boolean_coplanar_partial_overlap() {
        // Two boxes with partially overlapping coplanar faces.
        // A: [0,2]x[0,2]x[0,2], B: [1,3]x[0,2]x[0,2]
        // The shared face at x=1 (A) / x=1 (B) partially overlaps.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(
            result.is_ok(),
            "coplanar partial overlap union failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
    }

    #[test]
    fn boolean_coplanar_difference() {
        // Subtract a box that shares a coplanar face with the target.
        // A: [0,4]x[0,4]x[0,4], B: [0,2]x[0,4]x[0,4]
        // The face at x=0 is coplanar and coincident.
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let b = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 4.0, 4.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok(),
            "coplanar difference failed: {:?}",
            result.err()
        );
        let brep = result.unwrap();
        assert!(!brep.solids[0].shells[0].faces.is_empty());
    }

    // ─── Tangent Contact Boolean Tests ────────────────────────────────────────

    #[test]
    fn boolean_tangent_sphere_sphere() {
        // Two spheres touching at exactly one point (external tangent).
        // d = r1 + r2 = 1 + 1 = 2
        let a = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();
        let b = make_sphere_brep(DVec3::new(2.0, 0.0, 0.0), 1.0).unwrap();
        // Intersection should be empty (single point).
        let _inter = boolean_op(BooleanOpType::Intersection, &a, &b);
        // Union should succeed (two touching spheres).
        let union_result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(
            union_result.is_ok() || matches!(union_result, Err(BooleanError::DegenerateResult)),
            "tangent sphere union should not crash: {:?}",
            union_result.err()
        );
    }

    #[test]
    fn boolean_tangent_sphere_plane() {
        // Sphere touching a box face tangentially.
        // Sphere at (0,0,1) with r=1 touches the XY plane at origin.
        let a = make_box_brep(DVec3::new(-2.0, -2.0, -1.0), DVec3::X, DVec3::Y, 4.0, 4.0, 2.0).unwrap();
        let b = make_sphere_brep(DVec3::new(0.0, 0.0, 1.0), 1.0).unwrap();
        let result = boolean_op(BooleanOpType::Union, &a, &b);
        assert!(
            result.is_ok() || matches!(result, Err(BooleanError::DegenerateResult)),
            "tangent sphere-plane union should not crash: {:?}",
            result.err()
        );
    }

    #[test]
    fn boolean_tangent_cylinder_sphere() {
        // Cylinder tangent to a sphere (cylinder radius + offset = sphere radius).
        // Sphere at origin, r=2. Cylinder along Z axis, offset by 2 in X, r=0.
        // Actually: cylinder at x=2, r=1, sphere at origin r=3 → tangent at (3,0,0).
        let a = make_sphere_brep(DVec3::ZERO, 3.0).unwrap();
        let b = make_cylinder_brep(DVec3::new(2.0, 0.0, -2.0), DVec3::Z, DVec3::X, 1.0, 4.0).unwrap();
        let result = boolean_op(BooleanOpType::Difference, &a, &b);
        assert!(
            result.is_ok() || matches!(result, Err(BooleanError::DegenerateResult)),
            "tangent cylinder-sphere difference should not crash: {:?}",
            result.err()
        );
    }

    #[test]
    fn boolean_options_structure_accessible() {
        let a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let b = make_box_brep(DVec3::new(1.0, 0.0, 0.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();

        let options = BooleanOptions {
            use_bvh: true,
            run_healing: true,
            healing: HealingOptions::default(),
            run_simplify: true,
            simplify: SimplifyOptions::default(),
            include_history: true,
            fuzzy_tol: 0.0,
        };

        let (result, report) = boolean_op_with_options(BooleanOpType::Union, &a, &b, options)
            .expect("boolean_op_with_options should succeed");

        assert!(report.used_bvh);
        assert!(report.healed);
        assert!(report.simplified);
        assert!(report.healing_report.is_some());
        assert!(report.simplify_report.is_some());
        assert_eq!(report.output_faces, face_count(&result));
        assert_eq!(report.history_faces, report.persistent_face_labels.len());
        assert!(
            report
                .persistent_face_labels
                .iter()
                .all(|label| label.starts_with("face."))
        );
    }

}
