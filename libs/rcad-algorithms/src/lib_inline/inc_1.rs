


use rcad_kernel::BRep;

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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MakeConnectedScopeSeedMode {
 ShortEdges,
 NearDuplicateVertices,
 ToleranceTaggedEdges,
 MultiPcurveEdges,
 TopologySeamCandidates,
 Hybrid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MakeConnectedScopeSeedSource {
 Heuristic,
 History,
 HistoryAugmentedHeuristic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MakeConnectedScopeFallbackReason {
 InsufficientSeedCoverage,
 NoScopedChanges,
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
 /// Run baseline connectivity rebuilding (MakeConnected-style) after boolean.
 pub run_make_connected: bool,
 /// Tolerance used by connectivity rebuilding.
 pub make_connected_tolerance: f64,
 /// Maximum number of iterative make-connected passes.
 pub make_connected_max_passes: usize,
 /// Per-pass tolerance growth factor for iterative make-connected.
 pub make_connected_tolerance_growth: f64,
 /// Upper bound for make-connected tolerance growth.
 pub make_connected_tolerance_cap: f64,
 /// Enable scoped make-connected mode (local region only).
 pub make_connected_scoped: bool,
 /// Seed edge length threshold used to derive local scope vertices.
 pub make_connected_scope_seed_length: f64,
 /// Ring depth used when expanding history-derived seed edges in scoped mode.
 ///
 /// `0` keeps raw history edges only.
 /// `1` includes edges on faces adjacent to history edges (previous behavior).
 pub make_connected_scope_history_ring_depth: usize,
 /// When scoped make-connected makes no changes, retry with global scope.
 ///
 /// This keeps localized cleanup as the first attempt while preserving a
 /// broader recovery path for cases where scoped seeds miss the stressed
 /// region.
 pub make_connected_scope_fallback_to_global: bool,
 /// Minimum number of scoped seed vertices required before running the
 /// scoped pass.
 ///
 /// Values of `0` disable coverage-based fallback. Values `> 0` escalate
 /// directly to global make-connected when scoped seed coverage is smaller
 /// than this threshold.
 pub make_connected_scope_fallback_min_seed_vertices: usize,
 /// Minimum fraction of edges that must be covered by scoped seed edges
 /// before running the scoped pass.
 ///
 /// Values `<= 0` disable edge-ratio-based fallback. Values are clamped to
 /// the range `[0, 1]` when evaluated.
 pub make_connected_scope_fallback_min_seed_edge_coverage: f64,
 /// Minimum fraction of faces that must be touched by scoped seed edges
 /// before running the scoped pass.
 ///
 /// Values `<= 0` disable face-ratio-based fallback. Values are clamped to
 /// the range `[0, 1]` when evaluated.
 pub make_connected_scope_fallback_min_seed_face_coverage: f64,
 /// Multiplier applied to the base make-connected tolerance when scoped
 /// execution escalates to a global fallback pass.
 ///
 /// Values below `1.0` are clamped to `1.0`.
 pub make_connected_scope_global_fallback_tolerance_multiplier: f64,
 /// Maximum number of iterative passes used by global fallback.
 ///
 /// Values of `0` inherit `make_connected_max_passes`.
 pub make_connected_scope_global_fallback_max_passes: usize,
 /// Per-pass tolerance growth factor used by global fallback.
 ///
 /// Values `<= 0` inherit `make_connected_tolerance_growth`.
 pub make_connected_scope_global_fallback_tolerance_growth: f64,
 /// Upper cap for tolerance growth used by global fallback.
 ///
 /// Values `<= 0` inherit `make_connected_tolerance_cap`.
 pub make_connected_scope_global_fallback_tolerance_cap: f64,
 /// Seed derivation strategy for scoped mode.
 pub make_connected_scope_seed_mode: MakeConnectedScopeSeedMode,
 /// Minimum history-seed edge count before skipping heuristic augmentation.
 ///
 /// In scoped mode, if history-derived seed edges are fewer than this value,
 /// heuristic seed edges are unioned in to improve local coverage.
 pub make_connected_scope_min_history_edges: usize,
 /// Fuzzy tolerance for near-miss interference detection (analogous to
 /// `BOPAlgo_Options::SetFuzzyValue`).
 ///
 /// Values at or below zero select the default floor [`tolerance::TOLERANCE_ABS`] inside
 /// [`bopds::ds::DS::new_with_fuzzy`]. [`resolved_boolean_fuzzy_tol_for_ds`] matches that
 /// clamp for [`BooleanExecutionReport::effective_fuzzy_tol`]. For FEA / large-scale mechanical
 /// workflows, prefer [`BooleanRobustOptions::for_fea`] or [`BooleanRobustOptions::for_mechanical_multiscale`].
 pub fuzzy_tol: f64,
 /// Enable glue detection and fast-path merging for shared faces.
 ///
 /// Glue mode detects face pairs with identical geometry and opposite normals,
 /// then merges them directly without pave-filling. This is faster for
 /// contact/assembly scenarios.
 pub use_glue: bool,
 /// Tolerance for shared-face detection in glue mode.
 ///
 /// Controls how close edges must be to be considered "shared" (coplanar,
 /// coincident vertices, etc.). Defaults to `TOLERANCE_ABS`.
 ///
 /// [`boolean_op_with_options`] also raises this toward
 /// [`tolerance::combined_linear_tol_models`] when both operands are known (paired model bound;
 /// includes [`Self::fuzzy_tol`] when it is strictly positive).
 pub glue_tolerance: f64,
 /// After healing, make-connected, and simplify, run [`propagate_tolerances`] bottom-up
 /// with floor [`resolved_boolean_fuzzy_tol_for_ds`] so `GeomStore` tolerance arrays
 /// are sized and consistent with the effective pave fuzzy (FEA / multiscale preset: on).
 pub run_propagate_geom_tolerances: bool,
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
 run_make_connected: false,
 make_connected_tolerance: tolerance::TOLERANCE_ABS,
 make_connected_max_passes: 3,
 make_connected_tolerance_growth: 1.0,
 make_connected_tolerance_cap: tolerance::TOLERANCE_ABS * 1000.0,
 make_connected_scoped: false,
 make_connected_scope_seed_length: tolerance::TOLERANCE_ABS * 10.0,
 make_connected_scope_history_ring_depth: 1,
 make_connected_scope_fallback_to_global: true,
 make_connected_scope_fallback_min_seed_vertices: 1,
 make_connected_scope_fallback_min_seed_edge_coverage: 0.0,
 make_connected_scope_fallback_min_seed_face_coverage: 0.0,
 make_connected_scope_global_fallback_tolerance_multiplier: 1.0,
 make_connected_scope_global_fallback_max_passes: 0,
 make_connected_scope_global_fallback_tolerance_growth: 0.0,
 make_connected_scope_global_fallback_tolerance_cap: 0.0,
 make_connected_scope_seed_mode: MakeConnectedScopeSeedMode::Hybrid,
 make_connected_scope_min_history_edges: 2,
 fuzzy_tol: 0.0,
 use_glue: false,
 glue_tolerance: tolerance::TOLERANCE_ABS,
 run_propagate_geom_tolerances: false,
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
 pub made_connected: bool,
 pub make_connected_report: Option<MakeConnectedReport>,
 /// Seed mode used for scoped make-connected, if scoped mode was enabled.
 pub make_connected_scope_seed_mode: Option<MakeConnectedScopeSeedMode>,
 /// Configured history-ring depth used in scoped mode.
 pub make_connected_scope_history_ring_depth: Option<usize>,
 /// Seed source used in scoped mode.
 pub make_connected_scope_seed_source: Option<MakeConnectedScopeSeedSource>,
 /// Whether scoped make-connected escalated to a global fallback pass.
 pub make_connected_scope_fallback_applied: bool,
 /// Why scoped make-connected escalated to a global fallback pass.
 pub make_connected_scope_fallback_reason: Option<MakeConnectedScopeFallbackReason>,
 /// Report for the scoped make-connected phase when it was executed.
 pub make_connected_scope_scoped_report: Option<MakeConnectedReport>,
 /// Report for the global fallback make-connected phase when it was executed.
 pub make_connected_scope_global_fallback_report: Option<MakeConnectedReport>,
 /// Initial tolerance used for the global fallback phase, when executed.
 pub make_connected_scope_global_fallback_initial_tolerance: Option<f64>,
 /// Maximum passes configured for the global fallback phase, when executed.
 pub make_connected_scope_global_fallback_max_passes: Option<usize>,
 /// Ratio of scoped seed edges to total edges in the candidate shape.
 pub make_connected_scope_seed_edge_coverage: Option<f64>,
 /// Ratio of faces touched by scoped seed edges to total faces.
 pub make_connected_scope_seed_face_coverage: Option<f64>,
 /// Number of history-derived seed edges before union.
 pub make_connected_scope_history_seed_edge_count: usize,
 /// Number of heuristic-derived seed edges before union.
 pub make_connected_scope_heuristic_seed_edge_count: usize,
 /// Seed vertices used for scoped make-connected.
 pub make_connected_scope_seed_vertices: Vec<usize>,
 /// Seed edges used for scoped make-connected.
 pub make_connected_scope_seed_edges: Vec<usize>,
 /// Stable labels for scoped seed edges (orientation-insensitive).
 pub make_connected_scope_seed_edge_labels: Vec<String>,
 pub history_faces: usize,
 pub history_edges: usize,
 pub history_vertices: usize,
 pub history_shells: usize,
 pub history_solids: usize,
 pub persistent_face_labels: Vec<String>,
 pub persistent_edge_labels: Vec<String>,
 pub persistent_shell_labels: Vec<String>,
 pub persistent_solid_labels: Vec<String>,
 /// Full face/edge/vertex history when [`BooleanOptions::include_history`] was enabled.
 ///
 /// Populated from the boolean builder **before** optional healing / simplify; if those change
 /// topology, indices may not match the final [`BRep`] (same caveat as derived label fields).
 pub boolean_history: Option<BooleanHistory>,
 /// Per-attempt diagnostics recorded by `boolean_op_robust`.
 pub robust_attempts: Vec<BooleanRobustAttemptReport>,
 /// Number of retry attempts performed before success.
 pub retry_count: usize,
 /// Configured pave fuzzy ([`BooleanOptions::fuzzy_tol`]) for this run, **before**
 /// [`resolved_boolean_fuzzy_tol_for_ds`] clamp used inside [`bopds::ds::DS`].
 ///
 /// Use this (not [`Self::effective_fuzzy_tol`]) when re-merging
 /// [`HealingOptions`] so `combined_linear_tol_models` workspace pairing matches
 /// the boolean attempt (`fuzzy_tol > 0` vs `0`).
 pub configured_fuzzy_tol: f64,
 /// Fuzzy tolerance value that produced the final result.
 pub effective_fuzzy_tol: f64,
 /// Whether [`propagate_tolerances`] (bottom-up) ran after the boolean pipeline.
 pub propagated_geom_tolerances: bool,
}

/// Robust boolean retry controls.
#[derive(Debug, Clone)]
pub struct BooleanRobustOptions {
 /// Base execution options for each attempt.
 pub base: BooleanOptions,
 /// Additional fuzzy tolerance values to try when an attempt fails.
 pub fuzzy_retry_ladder: Vec<f64>,
 /// Retry policy controlling candidate generation after each failure.
 pub retry_policy: BooleanRetryPolicy,
 /// Configuration for extreme geometry handling.
 pub extreme_geometry: ExtremeGeometryRetryConfig,
}

/// Retry classes used by adaptive robust-boolean retry policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanRetryClass {
 /// Input is structurally invalid for retry (e.g. empty input).
 FatalInput,
 /// Missing geometry payload cannot be fixed by fuzzy escalation.
 IncompleteData,
 /// Topology degeneracy may be resolved by increased fuzzy tolerance.
 DegenerateTopology,
 /// Numeric instability often needs stronger fuzzy escalation first.
 NumericalInstability,
}

/// Retry-policy presets for robust boolean fuzzy escalation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanRetryPolicy {
 /// Conservative: only retry with ladder values larger than attempted fuzzy.
 Conservative,
 /// Adaptive: classify failures and choose escalation candidates by class.
 AdaptiveByFailureClass,
 /// Aggressive: retry ladder values plus multiplicative fuzzy boosts.
 Aggressive,
}

/// Retry strategy for extreme geometry conditions.
///
/// This policy extends the base retry mechanism to account for geometric
/// conditions that require specialized tolerance adjustments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtremeGeometryRetryPolicy {
 /// No extreme geometry handling (use base retry policy only).
 None,
 /// Detect extreme geometry and adjust tolerances before first attempt.
 PreAnalyze,
 /// Detect extreme geometry and use specialized retry ladder.
 AdaptiveTolerance,
 /// Full extreme geometry analysis with geometry-aware retry strategy.
 GeometryAware,
}

/// Configuration for extreme geometry retry handling.
#[derive(Debug, Clone)]
pub struct ExtremeGeometryRetryConfig {
 /// Policy to use for extreme geometry.
 pub policy: ExtremeGeometryRetryPolicy,
 /// Whether to check for near-tangent configurations.
 pub check_near_tangent: bool,
 /// Whether to check for high aspect ratio geometry.
 pub check_aspect_ratio: bool,
 /// Whether to check for degenerate geometry.
 pub check_degenerate: bool,
 /// Whether to check for size differences between inputs.
 pub check_size_difference: bool,
 /// Maximum fuzzy tolerance multiplier for extreme geometry.
 pub max_fuzzy_multiplier: f64,
 /// Number of additional retry steps to add for extreme geometry.
 pub extra_retry_steps: usize,
}

impl Default for ExtremeGeometryRetryConfig {
 fn default() -> Self {
 Self {
 policy: ExtremeGeometryRetryPolicy::AdaptiveTolerance,
 check_near_tangent: true,
 check_aspect_ratio: true,
 check_degenerate: true,
 check_size_difference: true,
 max_fuzzy_multiplier: 1000.0,
 extra_retry_steps: 2,
 }
 }
}

impl ExtremeGeometryRetryConfig {
 /// Create a configuration that skips all extreme geometry checks.
 pub fn none() -> Self {
 Self {
 policy: ExtremeGeometryRetryPolicy::None,
 check_near_tangent: false,
 check_aspect_ratio: false,
 check_degenerate: false,
 check_size_difference: false,
 max_fuzzy_multiplier: 1.0,
 extra_retry_steps: 0,
 }
 }

 /// Create a configuration for geometry-aware retry.
 pub fn geometry_aware() -> Self {
 Self {
 policy: ExtremeGeometryRetryPolicy::GeometryAware,
 ..Default::default()
 }
 }

 /// Build a specialized retry ladder based on extreme geometry analysis.
 pub fn build_retry_ladder(
 &self,
 base_ladder: &[f64],
 analysis: &ExtremeGeometryAnalysis,
 ) -> Vec<f64> {
 if self.policy == ExtremeGeometryRetryPolicy::None {
 return base_ladder.to_vec();
 }

 let mut ladder = base_ladder.to_vec();

 // Add tolerance adjustments for near-tangent configurations
 if self.check_near_tangent && !analysis.near_tangent_configs.is_empty() {
 for config in &analysis.near_tangent_configs {
 let tol = config.suggested_fuzzy_adjustment;
 if !ladder
 .iter()
 .any(|&t| (t - tol).abs() < tolerance::TOLERANCE_ABS)
 {
 ladder.push(tol);
 }
 }
 }

 // Add tolerance adjustments for high aspect ratio edges
 if self.check_aspect_ratio {
 for edge in &analysis.high_aspect_ratio_edges {
 if edge.is_problematic {
 let tol = tolerance::TOLERANCE_ABS * edge.suggested_tolerance_multiplier;
 if !ladder
 .iter()
 .any(|&t| (t - tol).abs() < tolerance::TOLERANCE_ABS)
 {
 ladder.push(tol);
 }
 }
 }
 }

 // Add tolerance adjustments for size difference
 if self.check_size_difference
 && let Some(ref sd) = analysis.size_difference
 && sd.is_extreme
 {
 let tol = tolerance::TOLERANCE_ABS * sd.suggested_tolerance_multiplier;
 if !ladder
 .iter()
 .any(|&t| (t - tol).abs() < tolerance::TOLERANCE_ABS)
 {
 ladder.push(tol);
 }
 }

 // Add the recommended fuzzy tolerance from the analysis
 if analysis.recommended_fuzzy_tolerance > tolerance::TOLERANCE_ABS {
 let tol = analysis
 .recommended_fuzzy_tolerance
 .min(tolerance::TOLERANCE_ABS * self.max_fuzzy_multiplier);
 if !ladder
 .iter()
 .any(|&t| (t - tol).abs() < tolerance::TOLERANCE_ABS)
 {
 ladder.push(tol);
 }
 }

 // Sort and deduplicate
 ladder.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
 ladder.dedup_by(|a, b| (*a - *b).abs() < tolerance::TOLERANCE_ABS);

 // Cap the ladder
 ladder.truncate(base_ladder.len() + self.extra_retry_steps + 1);

 ladder
 }
}

/// Per-attempt diagnostics for robust boolean retry execution.
#[derive(Debug, Clone)]
pub struct BooleanRobustAttemptReport {
 /// Fuzzy tolerance used for this attempt.
 pub fuzzy_tol: f64,
 /// Whether this attempt succeeded.
 pub success: bool,
 /// Escalation round used for this attempt.
 pub retry_round: usize,
 /// Failure class that scheduled this retry attempt.
 pub origin_retry_class: Option<BooleanRetryClass>,
 /// Whether scoped make-connected was enabled for this attempt.
 pub make_connected_scoped_enabled: bool,
 /// Effective scoped seed mode configured for this attempt.
 pub make_connected_scope_seed_mode: Option<MakeConnectedScopeSeedMode>,
 /// Effective history ring depth configured for this attempt.
 pub make_connected_scope_history_ring_depth: Option<usize>,
 /// Effective scoped seed length configured for this attempt.
 pub make_connected_scope_seed_length: Option<f64>,
 /// Effective minimum history-edge threshold before heuristic augmentation.
 pub make_connected_scope_min_history_edges: Option<usize>,
 /// Effective scoped seed source observed during this attempt.
 pub make_connected_scope_seed_source: Option<MakeConnectedScopeSeedSource>,
 /// Number of history-derived scoped seed edges observed during this attempt.
 pub make_connected_scope_history_seed_edge_count: Option<usize>,
 /// Number of heuristic-derived scoped seed edges observed during this attempt.
 pub make_connected_scope_heuristic_seed_edge_count: Option<usize>,
 /// Number of scoped seed vertices observed during this attempt.
 pub make_connected_scope_seed_vertex_count: Option<usize>,
 /// Number of scoped seed edges observed during this attempt.
 pub make_connected_scope_seed_edge_count: Option<usize>,
 /// Whether glue mode was enabled for this attempt.
 pub used_glue: bool,
 /// Effective glue tolerance configured for this attempt.
 pub glue_tolerance: f64,
 /// Retry classification for a failed attempt.
 pub retry_class: Option<BooleanRetryClass>,
 /// Debug message for a failed attempt.
 pub error_message: Option<String>,
 /// Face count of the successful result.
 pub output_faces: Option<usize>,
 /// Whether make-connected ran during this attempt.
 pub made_connected: bool,
 /// Whether scoped make-connected escalated to global fallback.
 pub make_connected_scope_fallback_applied: bool,
 /// Scoped fallback reason, when present.
 pub make_connected_scope_fallback_reason: Option<MakeConnectedScopeFallbackReason>,
 /// Scoped seed edge coverage ratio for this attempt.
 pub make_connected_scope_seed_edge_coverage: Option<f64>,
 /// Scoped seed face coverage ratio for this attempt.
 pub make_connected_scope_seed_face_coverage: Option<f64>,
 /// Global fallback initial tolerance used in this attempt, when present.
 pub make_connected_scope_global_fallback_initial_tolerance: Option<f64>,
 /// Global fallback max-passes used in this attempt, when present.
 pub make_connected_scope_global_fallback_max_passes: Option<usize>,
}

impl Default for BooleanRobustOptions {
 fn default() -> Self {
 Self {
 base: BooleanOptions::default(),
 fuzzy_retry_ladder: boolean_fuzzy_ladder_scaled(tolerance::TOLERANCE_ABS, None),
 retry_policy: BooleanRetryPolicy::AdaptiveByFailureClass,
 extreme_geometry: ExtremeGeometryRetryConfig::default(),
 }
 }
}

impl BooleanRobustOptions {
 /// Preset for **FEA-oriented** booleans: scale-aware fuzzy/glue, glue, healing,
 /// make-connected, and **bottom-up `propagate_tolerances`** enabled. Use with
 /// [`boolean_op_robust`] for mesh-friendly watertight recovery.
 pub fn for_fea(a: &BRep, b: &BRep) -> Self {
 let ctx = tolerance::ToleranceContext::from_two_breps(a, b);
 let fuzzy = ctx.adaptive_linear(tolerance::ToleranceLevel::Normal);
 let glue = ctx.adaptive_linear(tolerance::ToleranceLevel::Normal);
 let mut base = BooleanOptions::default();
 base.use_glue = true;
 base.glue_tolerance = glue;
 base.fuzzy_tol = fuzzy;
 base.run_make_connected = true;
 base.run_healing = true;
 base.run_propagate_geom_tolerances = true;
 base.make_connected_tolerance = ctx.adaptive_linear(tolerance::ToleranceLevel::Normal);
 base.make_connected_tolerance_cap = ctx.adaptive_linear(tolerance::ToleranceLevel::Coarse);
 Self {
 base,
 fuzzy_retry_ladder: tolerance::boolean_fuzzy_ladder_scaled(fuzzy, None),
 retry_policy: BooleanRetryPolicy::AdaptiveByFailureClass,
 extreme_geometry: ExtremeGeometryRetryConfig::default(),
 }
 }

 /// Preset for **mechanical multi-scale** assemblies: relaxed starting fuzzy, wider retry
 /// ladder, and geometry-aware extreme-geometry escalation.
 pub fn for_mechanical_multiscale(a: &BRep, b: &BRep) -> Self {
 let ctx = tolerance::ToleranceContext::from_two_breps(a, b);
 let fuzzy = ctx.adaptive_linear(tolerance::ToleranceLevel::Relaxed);
 let glue = ctx.adaptive_linear(tolerance::ToleranceLevel::Normal);
 let coarse = ctx.adaptive_linear(tolerance::ToleranceLevel::Coarse);
 let mut base = BooleanOptions::default();
 base.use_glue = true;
 base.glue_tolerance = glue;
 base.fuzzy_tol = fuzzy;
 base.run_make_connected = true;
 base.run_healing = true;
 base.run_propagate_geom_tolerances = true;
 base.make_connected_tolerance = ctx.adaptive_linear(tolerance::ToleranceLevel::Relaxed);
 base.make_connected_tolerance_cap = ctx.adaptive_linear(tolerance::ToleranceLevel::Coarse);
 Self {
 base,
 fuzzy_retry_ladder: tolerance::boolean_fuzzy_ladder_scaled(fuzzy, Some(coarse)),
 retry_policy: BooleanRetryPolicy::AdaptiveByFailureClass,
 extreme_geometry: ExtremeGeometryRetryConfig::geometry_aware(),
 }
 }
}

/// Effective fuzzy tolerance used inside [`bopds::ds::DS`] (see [`bopds::ds::DS::new_with_fuzzy`]).
///
/// Values below [`tolerance::TOLERANCE_ABS`] clamp up to that floor. Use this for diagnostics
/// ([`BooleanExecutionReport::effective_fuzzy_tol`]) so reports match runtime behavior when
/// [`BooleanOptions::fuzzy_tol`] is `0.0` (閳ユ竸efault fuzzy閳?.
#[inline]
pub fn resolved_boolean_fuzzy_tol_for_ds(configured_fuzzy: f64) -> f64 {
 configured_fuzzy.max(tolerance::TOLERANCE_ABS)
}

/// Raises [`BooleanOptions`] glue / make-connected bands by OCCT-style
/// [`tolerance::combined_linear_tol_models`] over the two operands (and optional fuzzy workspace).
///
/// Also lifts [`BooleanOptions::healing`] so post-boolean [`analyze_and_heal`] uses at least the
/// same linear floor as glue / make-connected and [`resolved_boolean_fuzzy_tol_for_ds`], avoiding
/// repair passes that stay tighter than the pave / fuzzy context.
///
/// Idempotent (`max`). Call at binary boolean entry whenever both [`BRep`] operands are known.
fn merge_pairwise_model_tol_into_boolean_options(options: &mut BooleanOptions, a: &BRep, b: &BRep) {
 let base_ctx = tolerance::ToleranceContext::from_two_breps(a, b);
 let fuzzy_user = options.fuzzy_tol.max(0.0);
 let ctx = tolerance::ToleranceContext::new(base_ctx.adaptive, fuzzy_user);
 let use_workspace = options.fuzzy_tol > 0.0;

 let glue_floor = tolerance::combined_linear_tol_models(
 &ctx,
 tolerance::ToleranceLevel::Strict,
 use_workspace,
 a,
 b,
 );
 let mc_floor = tolerance::combined_linear_tol_models(
 &ctx,
 tolerance::ToleranceLevel::Normal,
 use_workspace,
 a,
 b,
 );
 let mc_cap_floor = tolerance::combined_linear_tol_models(
 &ctx,
 tolerance::ToleranceLevel::Coarse,
 use_workspace,
 a,
 b,
 );

 options.glue_tolerance = options.glue_tolerance.max(glue_floor);
 options.make_connected_tolerance = options.make_connected_tolerance.max(mc_floor);
 options.make_connected_tolerance_cap = options.make_connected_tolerance_cap.max(mc_cap_floor);

 let heal_floor = mc_floor
 .max(glue_floor)
 .max(resolved_boolean_fuzzy_tol_for_ds(options.fuzzy_tol));
 let mut h = options.healing;
 h.tolerance = h.tolerance.max(heal_floor);
 h.make_connected_tolerance = h.make_connected_tolerance.max(options.make_connected_tolerance);
 h.make_connected_tolerance_cap = h
 .make_connected_tolerance_cap
 .max(options.make_connected_tolerance_cap);
 options.healing = h;
}

/// Lift [`HealingOptions`]'s repair / make-connected tolerances using pairwise
/// [`combined_linear_tol_models`] over `a` and `b` and an optional pave fuzzy (`fuzzy_tol`),
/// matching the healing branch inside [`merge_pairwise_model_tol_into_boolean_options`].
///
/// Caller fields are preserved via `max` against computed floors. Use when running
/// [`analyze_and_heal`] after an operation whose operands are known but options were not merged
/// through [`BooleanOptions`] (for example [`boolean_op_healed_with_options`] or split imprint steps).
pub fn align_healing_options_with_boolean_operands(
 healing: &mut HealingOptions,
 a: &BRep,
 b: &BRep,
 fuzzy_tol: f64,
) {
 let mut bridge = BooleanOptions::default();
 bridge.fuzzy_tol = fuzzy_tol.max(0.0);
 bridge.healing = *healing;
 merge_pairwise_model_tol_into_boolean_options(&mut bridge, a, b);
 *healing = bridge.healing;
}

/// Like [`align_healing_options_with_boolean_operands`], but uses
/// [`BooleanExecutionReport::configured_fuzzy_tol`] so post-boolean healing stays consistent
/// with the attempt閳ユ獨 workspace flag (e.g. `fuzzy_tol == 0` vs strictly positive user fuzzy).
pub fn align_healing_options_after_boolean_execution(
 healing: &mut HealingOptions,
 a: &BRep,
 b: &BRep,
 execution: &BooleanExecutionReport,
) {
 align_healing_options_with_boolean_operands(
 healing,
 a,
 b,
 execution.configured_fuzzy_tol,
 );
}

/// Build ordered fuzzy values for robust retry.
///
/// First element is always the initial fuzzy value (clamped to >= 0).
/// Ladder values <= 0 are skipped; duplicates (within epsilon) are removed.
pub fn boolean_retry_fuzzy_values(initial: f64, ladder: &[f64]) -> Vec<f64> {
 let mut values = vec![initial.max(0.0)];
 for &v in ladder {
 if v <= 0.0 {
 continue;
 }
 if !values.iter().any(|e| (*e - v).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP) {
 values.push(v);
 }
 }
 values
}

/// Classify boolean execution failures for adaptive retry policies.
pub fn classify_boolean_retry(err: &BooleanError) -> BooleanRetryClass {
 match err {
 BooleanError::InvalidOperation => BooleanRetryClass::FatalInput,
 BooleanError::TooFewArguments => BooleanRetryClass::FatalInput,
 BooleanError::NoFiller => BooleanRetryClass::FatalInput,
 BooleanError::BOPNotAllowed => BooleanRetryClass::FatalInput,
 BooleanError::BOPNotSet => BooleanRetryClass::FatalInput,
 BooleanError::EmptyShape => BooleanRetryClass::FatalInput,
 BooleanError::EmptyInput => BooleanRetryClass::FatalInput,
 BooleanError::MissingGeometry(_) => BooleanRetryClass::IncompleteData,
 BooleanError::DegenerateResult => BooleanRetryClass::DegenerateTopology,
 BooleanError::NumericalFailure(_) => BooleanRetryClass::NumericalInstability,
 BooleanError::EmptyCollection(_) => BooleanRetryClass::DegenerateTopology,
 BooleanError::InvalidResult(_) => BooleanRetryClass::DegenerateTopology,
 BooleanError::IncompleteIntersection(_) => BooleanRetryClass::DegenerateTopology,
 BooleanError::SelfIntersection(_) => BooleanRetryClass::DegenerateTopology,
 BooleanError::OpenShell { .. } => BooleanRetryClass::DegenerateTopology,
 }
}

/// Build next fuzzy values based on the last failure type.
///
/// Returned values are positive, deduplicated, and ordered from smaller to
/// larger escalation.
pub fn boolean_retry_ladder_for_error(
 attempted_fuzzy: f64,
 ladder: &[f64],
 err: &BooleanError,
) -> Vec<f64> {
 let class = classify_boolean_retry(err);
 let mut out: Vec<f64> = Vec::new();
 let mut push_unique = |v: f64| {
 if v <= 0.0 {
 return;
 }
 if !out.iter().any(|e| (*e - v).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP) {
 out.push(v);
 }
 };

 match class {
 BooleanRetryClass::FatalInput | BooleanRetryClass::IncompleteData => {}
 BooleanRetryClass::DegenerateTopology => {
 for &v in ladder {
 if v > attempted_fuzzy {
 push_unique(v);
 }
 }
 }
 BooleanRetryClass::NumericalInstability => {
 let baseline = if attempted_fuzzy > 0.0 {
 attempted_fuzzy
 } else {
 tolerance::TOLERANCE_ABS
 };
 push_unique(baseline * 10.0);
 push_unique(baseline * 100.0);
 for &v in ladder {
 if v > attempted_fuzzy {
 push_unique(v);
 }
 }
 }
 }

 out
}

/// Build next fuzzy values using the configured retry policy.
pub fn boolean_retry_ladder_for_error_with_policy(
 attempted_fuzzy: f64,
 ladder: &[f64],
 err: &BooleanError,
 policy: BooleanRetryPolicy,
) -> Vec<f64> {
 let mut out: Vec<f64> = Vec::new();
 let mut push_unique = |v: f64| {
 if v <= 0.0 {
 return;
 }
 if !out.iter().any(|e| (*e - v).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP) {
 out.push(v);
 }
 };

 match policy {
 BooleanRetryPolicy::AdaptiveByFailureClass => {
 return boolean_retry_ladder_for_error(attempted_fuzzy, ladder, err);
 }
 BooleanRetryPolicy::Conservative => {
 match classify_boolean_retry(err) {
 BooleanRetryClass::FatalInput | BooleanRetryClass::IncompleteData => return out,
 _ => {}
 }
 for &v in ladder {
 if v > attempted_fuzzy {
 push_unique(v);
 }
 }
 }
 BooleanRetryPolicy::Aggressive => {
 match classify_boolean_retry(err) {
 BooleanRetryClass::FatalInput | BooleanRetryClass::IncompleteData => return out,
 _ => {}
 }
 let baseline = if attempted_fuzzy > 0.0 {
 attempted_fuzzy
 } else {
 tolerance::TOLERANCE_ABS
 };
 for &v in ladder {
 if v > attempted_fuzzy {
 push_unique(v);
 }
 }
 push_unique(baseline * 10.0);
 push_unique(baseline * 100.0);
 }
 }

 out
}

fn boolean_retry_followup_attempts(
 attempted_fuzzy: f64,
 ladder: &[f64],
 err: &BooleanError,
 policy: BooleanRetryPolicy,
 origin_retry_class: Option<BooleanRetryClass>,
 retry_round: usize,
 max_retry_escalation_rounds: usize,
 attempted_scoped_cleanup_enabled: bool,
) -> Vec<(f64, Option<BooleanRetryClass>, usize)> {
 let retry_class = classify_boolean_retry(err);
 if matches!(
 retry_class,
 BooleanRetryClass::FatalInput | BooleanRetryClass::IncompleteData
 ) {
 return Vec::new();
 }

 let fuzzy_candidate_round = if origin_retry_class == Some(retry_class) {
 (retry_round + 1).min(max_retry_escalation_rounds)
 } else {
 0
 };
 let strategy_candidate_round = if origin_retry_class == Some(retry_class) {
 retry_round + 1
 } else {
 1
 };
 let can_escalate_strategy = retry_round < max_retry_escalation_rounds;
 let strategy_already_global_biased =
 origin_retry_class.is_some() && !attempted_scoped_cleanup_enabled;
 let fuzzy_candidates =
 boolean_retry_ladder_for_error_with_policy(attempted_fuzzy, ladder, err, policy);

 let mut out: Vec<(f64, Option<BooleanRetryClass>, usize)> = Vec::new();
 let mut push_unique = |candidate: (f64, Option<BooleanRetryClass>, usize)| {
 if candidate.0 <= 0.0 {
 return;
 }
 if !out.iter().any(|existing| {
 (existing.0 - candidate.0).abs() <= tolerance::TOLERANCE_FLOAT_DEDUP
 && existing.1 == candidate.1
 && existing.2 == candidate.2
 }) {
 out.push(candidate);
 }
 };

 if matches!(retry_class, BooleanRetryClass::DegenerateTopology)
 && can_escalate_strategy
 && !strategy_already_global_biased
 {
 push_unique((attempted_fuzzy, Some(retry_class), strategy_candidate_round));
 }

 for candidate in fuzzy_candidates {
 push_unique((candidate, Some(retry_class), fuzzy_candidate_round));
 }

 if matches!(retry_class, BooleanRetryClass::NumericalInstability)
 && can_escalate_strategy
 && !strategy_already_global_biased
 {
 push_unique((attempted_fuzzy, Some(retry_class), strategy_candidate_round));
 }

 out
}

fn tune_boolean_options_for_retry_class(
 options: &mut BooleanOptions,
 retry_class: Option<BooleanRetryClass>,
 retry_round: usize,
) {
 let Some(retry_class) = retry_class else {
 return;
 };

 let base_tol = options
 .make_connected_tolerance
 .max(options.glue_tolerance)
 .max(tolerance::TOLERANCE_ABS);

 match retry_class {
 BooleanRetryClass::FatalInput | BooleanRetryClass::IncompleteData => {}
 BooleanRetryClass::DegenerateTopology => {
 options.use_glue = true;
 options.glue_tolerance = options
 .glue_tolerance
 .max(base_tol * 10.0 * (retry_round as f64 + 1.0));

 if !options.run_make_connected {
 return;
 }

 options.make_connected_max_passes =
 options.make_connected_max_passes.max(4 + retry_round);
 options.make_connected_tolerance_growth = options
 .make_connected_tolerance_growth
 .max(2.0 + retry_round as f64);
 options.make_connected_tolerance_cap = options
 .make_connected_tolerance_cap
 .max(base_tol * 1000.0 * (retry_round as f64 + 1.0));

 if options.make_connected_scoped && retry_round >= 2 {
 options.make_connected_scoped = false;
 }

 if options.make_connected_scoped {
 options.make_connected_scope_seed_length = options
 .make_connected_scope_seed_length
 .max(base_tol * 10.0 * (retry_round as f64 + 1.0));
 options.make_connected_scope_history_ring_depth = options
 .make_connected_scope_history_ring_depth
 .max(2 + retry_round);
 options.make_connected_scope_min_history_edges = options
 .make_connected_scope_min_history_edges
 .max(2 + retry_round);
 options.make_connected_scope_seed_mode =
 match options.make_connected_scope_seed_mode {
 MakeConnectedScopeSeedMode::ShortEdges
 | MakeConnectedScopeSeedMode::NearDuplicateVertices
 | MakeConnectedScopeSeedMode::ToleranceTaggedEdges => {
 MakeConnectedScopeSeedMode::TopologySeamCandidates
 }
 MakeConnectedScopeSeedMode::MultiPcurveEdges => {
 MakeConnectedScopeSeedMode::Hybrid
 }
 mode => mode,
 };
 options.make_connected_scope_fallback_to_global = true;
 options.make_connected_scope_fallback_min_seed_vertices = options
 .make_connected_scope_fallback_min_seed_vertices
 .max(2 + retry_round);
 options.make_connected_scope_fallback_min_seed_edge_coverage = options
 .make_connected_scope_fallback_min_seed_edge_coverage
 .max((0.25 + 0.1 * retry_round as f64).min(1.0));
 options.make_connected_scope_fallback_min_seed_face_coverage = options
 .make_connected_scope_fallback_min_seed_face_coverage
 .max((0.25 + 0.1 * retry_round as f64).min(1.0));
 options.make_connected_scope_global_fallback_tolerance_multiplier = options
 .make_connected_scope_global_fallback_tolerance_multiplier
 .max(10.0 * (retry_round as f64 + 1.0));
 options.make_connected_scope_global_fallback_max_passes = options
 .make_connected_scope_global_fallback_max_passes
 .max(4 + retry_round);
 options.make_connected_scope_global_fallback_tolerance_growth = options
 .make_connected_scope_global_fallback_tolerance_growth
 .max(2.0 + retry_round as f64);
 options.make_connected_scope_global_fallback_tolerance_cap = options
 .make_connected_scope_global_fallback_tolerance_cap
 .max(base_tol * 1000.0 * (retry_round as f64 + 1.0));
 }
 }
 BooleanRetryClass::NumericalInstability => {
 options.use_glue = true;
 options.glue_tolerance = options
 .glue_tolerance
 .max(base_tol * 100.0 * (retry_round as f64 + 1.0));

 if !options.run_make_connected {
 return;
 }

 options.make_connected_max_passes =
 options.make_connected_max_passes.max(5 + retry_round);
 options.make_connected_tolerance_growth = options
 .make_connected_tolerance_growth
 .max(10.0 + 5.0 * retry_round as f64);
 options.make_connected_tolerance_cap = options
 .make_connected_tolerance_cap
 .max(base_tol * 10_000.0 * (retry_round as f64 + 1.0));

 if options.make_connected_scoped && retry_round >= 2 {
 options.make_connected_scoped = false;
 }

 if options.make_connected_scoped {
 options.make_connected_scope_seed_length = options
 .make_connected_scope_seed_length
 .max(base_tol * 100.0 * (retry_round as f64 + 1.0));
 options.make_connected_scope_history_ring_depth = options
 .make_connected_scope_history_ring_depth
 .max(3 + retry_round);
 options.make_connected_scope_min_history_edges = options
 .make_connected_scope_min_history_edges
 .max(3 + retry_round);
 options.make_connected_scope_seed_mode = MakeConnectedScopeSeedMode::Hybrid;
 options.make_connected_scope_fallback_to_global = true;
 options.make_connected_scope_fallback_min_seed_vertices = options
 .make_connected_scope_fallback_min_seed_vertices
 .max(2 + retry_round);
 options.make_connected_scope_fallback_min_seed_edge_coverage = options
 .make_connected_scope_fallback_min_seed_edge_coverage
 .max((0.5 + 0.1 * retry_round as f64).min(1.0));
 options.make_connected_scope_fallback_min_seed_face_coverage = options
 .make_connected_scope_fallback_min_seed_face_coverage
 .max((0.5 + 0.1 * retry_round as f64).min(1.0));
 options.make_connected_scope_global_fallback_tolerance_multiplier = options
 .make_connected_scope_global_fallback_tolerance_multiplier
 .max(100.0 * (retry_round as f64 + 1.0));
 options.make_connected_scope_global_fallback_max_passes = options
 .make_connected_scope_global_fallback_max_passes
 .max(5 + retry_round);
 options.make_connected_scope_global_fallback_tolerance_growth = options
 .make_connected_scope_global_fallback_tolerance_growth
 .max(10.0 + 5.0 * retry_round as f64);
 options.make_connected_scope_global_fallback_tolerance_cap = options
 .make_connected_scope_global_fallback_tolerance_cap
 .max(base_tol * 10_000.0 * (retry_round as f64 + 1.0));
 }
 }
 }
}

/// Tune boolean options for a specific detailed failure class.
///
fn merge_make_connected_reports(
 mut initial: MakeConnectedReport,
 fallback: MakeConnectedReport,
) -> MakeConnectedReport {
 initial.vertices_merged += fallback.vertices_merged;
 initial.small_edges_removed += fallback.small_edges_removed;
 initial.passes_run += fallback.passes_run;
 initial.converged = fallback.converged;
 initial.final_tolerance = fallback.final_tolerance;
 initial.tolerance_cap_applied |= fallback.tolerance_cap_applied;
 initial
}

fn run_make_connected_for_boolean_output(
 brep: &BRep,
 history: Option<&BooleanHistory>,
 options: &BooleanOptions,
 report: &mut BooleanExecutionReport,
) -> (BRep, MakeConnectedReport) {
 let global_fallback_tolerance = options
 .make_connected_tolerance
 .max(tolerance::TOLERANCE_ABS)
 * options
 .make_connected_scope_global_fallback_tolerance_multiplier
 .max(1.0);
 let global_fallback_max_passes = if options.make_connected_scope_global_fallback_max_passes > 0
 {
 options.make_connected_scope_global_fallback_max_passes
 } else {
 options.make_connected_max_passes
 };
 let global_fallback_tolerance_growth =
 if options.make_connected_scope_global_fallback_tolerance_growth > 0.0 {
 options.make_connected_scope_global_fallback_tolerance_growth
 } else {
 options.make_connected_tolerance_growth
 };
 let global_fallback_tolerance_cap =
 if options.make_connected_scope_global_fallback_tolerance_cap > 0.0 {
 options.make_connected_scope_global_fallback_tolerance_cap
 } else {
 options.make_connected_tolerance_cap
 };

 if !options.make_connected_scoped {
 return make_connected_iterative_with_growth_cap(
 brep,
 options.make_connected_tolerance,
 options.make_connected_max_passes,
 options.make_connected_tolerance_growth,
 options.make_connected_tolerance_cap,
 );
 }

 let seed = options
 .make_connected_scope_seed_length
 .max(options.make_connected_tolerance);
 let (mut scope_seed_edges, history_seed_edges, heuristic_seed_edges, seed_source) =
 select_scoped_seed_edges(
 brep,
 history,
 seed,
 options.make_connected_scope_seed_mode,
 options.make_connected_scope_history_ring_depth,
 options.make_connected_scope_min_history_edges,
 );
 let mut scope_vertices =
 make_connected_seed_vertices(brep, seed, options.make_connected_scope_seed_mode);
 scope_vertices.extend(make_connected_seed_vertices_from_edge_ids(
 brep,
 &scope_seed_edges,
 ));
 scope_vertices.sort_unstable();
 scope_vertices.dedup();
 scope_seed_edges.sort_unstable();
 scope_seed_edges.dedup();

 report.make_connected_scope_seed_mode = Some(options.make_connected_scope_seed_mode);
 report.make_connected_scope_history_ring_depth =
 Some(options.make_connected_scope_history_ring_depth);
 report.make_connected_scope_seed_source = Some(seed_source);
 report.make_connected_scope_history_seed_edge_count = history_seed_edges;
 report.make_connected_scope_heuristic_seed_edge_count = heuristic_seed_edges;
 report.make_connected_scope_seed_vertices = scope_vertices.clone();
 report.make_connected_scope_seed_edge_labels =
 make_connected_seed_edge_labels(brep, &scope_seed_edges);
 report.make_connected_scope_seed_edges = scope_seed_edges;
 let seed_edge_coverage = if brep.edges.is_empty() {
 0.0
 } else {
 report.make_connected_scope_seed_edges.len() as f64 / brep.edges.len() as f64
 };
 report.make_connected_scope_seed_edge_coverage = Some(seed_edge_coverage);
 let mut seed_face_set = std::collections::BTreeSet::new();
 for &ei in &report.make_connected_scope_seed_edges {
 for fi in rcad_kernel::edge_adjacent_faces(brep, ei) {
 seed_face_set.insert(fi);
 }
 }
 let total_faces = face_count_of(brep);
 let seed_face_coverage = if total_faces == 0 {
 0.0
 } else {
 seed_face_set.len() as f64 / total_faces as f64
 };
 report.make_connected_scope_seed_face_coverage = Some(seed_face_coverage);

 let min_seed_vertices = options.make_connected_scope_fallback_min_seed_vertices;
 let min_seed_edge_coverage = options
 .make_connected_scope_fallback_min_seed_edge_coverage
 .clamp(0.0, 1.0);
 let min_seed_face_coverage = options
 .make_connected_scope_fallback_min_seed_face_coverage
 .clamp(0.0, 1.0);
 if options.make_connected_scope_fallback_to_global
 && ((min_seed_vertices > 0 && scope_vertices.len() < min_seed_vertices)
 || (min_seed_edge_coverage > 0.0 && seed_edge_coverage < min_seed_edge_coverage)
 || (min_seed_face_coverage > 0.0 && seed_face_coverage < min_seed_face_coverage))
 {
 let (global_connected, global_report) = make_connected_iterative_with_growth_cap(
 brep,
 global_fallback_tolerance,
 global_fallback_max_passes,
 global_fallback_tolerance_growth,
 global_fallback_tolerance_cap,
 );
 report.make_connected_scope_fallback_applied = true;
 report.make_connected_scope_fallback_reason =
 Some(MakeConnectedScopeFallbackReason::InsufficientSeedCoverage);
 report.make_connected_scope_global_fallback_initial_tolerance =
 Some(global_fallback_tolerance);
 report.make_connected_scope_global_fallback_max_passes = Some(global_fallback_max_passes);
 report.make_connected_scope_global_fallback_report = Some(global_report.clone());
 return (global_connected, global_report);
 }

 let (scoped_connected, scoped_report) = make_connected_iterative_scoped_with_growth_cap(
 brep,
 &scope_vertices,
 options.make_connected_tolerance,
 options.make_connected_max_passes,
 options.make_connected_tolerance_growth,
 options.make_connected_tolerance_cap,
 );
 report.make_connected_scope_scoped_report = Some(scoped_report.clone());
 let scoped_no_changes =
 scoped_report.vertices_merged == 0 && scoped_report.small_edges_removed == 0;

 if options.make_connected_scope_fallback_to_global && scoped_no_changes {
 let (global_connected, global_report) = make_connected_iterative_with_growth_cap(
 &scoped_connected,
 global_fallback_tolerance,
 global_fallback_max_passes,
 global_fallback_tolerance_growth,
 global_fallback_tolerance_cap,
 );
 report.make_connected_scope_fallback_applied = true;
 report.make_connected_scope_fallback_reason =
 Some(MakeConnectedScopeFallbackReason::NoScopedChanges);
 report.make_connected_scope_global_fallback_initial_tolerance =
 Some(global_fallback_tolerance);
 report.make_connected_scope_global_fallback_max_passes = Some(global_fallback_max_passes);
 report.make_connected_scope_global_fallback_report = Some(global_report.clone());
 return (
 global_connected,
 merge_make_connected_reports(scoped_report, global_report),
 );
 }

 (scoped_connected, scoped_report)
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

 let mut step_map: std::collections::BTreeMap<usize, usize> =
 std::collections::BTreeMap::new();
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

fn brep_shell_face_count(brep: &BRep) -> usize {
 brep.solids
 .iter()
 .flat_map(|s| &s.shells)
 .flat_map(|sh| &sh.faces)
 .count()
}

/// OCCT `BRepAlgoAPI_Cut` yields an empty shape when operands coincide (e.g. two identical
/// `box` definitions in `bopcut`). Match that without forcing every `DegenerateResult` from the
/// builder to mean 閳ユ竼mpty閳?
/// True when every face has geometry and every face surface is a plane (e.g. `make_box_brep` solids).
///
/// Used to gate 閳ユ笡lanar zero-volume sliver 閳?empty intersection閳?heuristics: operands that include
/// spheres/cylinders etc. can still yield all-plane *wrong* shells with `volume 閳?0`; we must not
/// collapse those to empty or OCCT sphere閳ユ彽ox cases regress to `total_surface_area == 0`.
fn brep_is_pure_plane_solid(brep: &BRep) -> bool {
 let nf = face_count_of(brep);
 if nf == 0 {
 return false;
 }
 if brep.geom.face_surface.len() != nf {
 return false;
 }
 for slot in &brep.geom.face_surface {
 let Some(si) = *slot else {
 return false;
 };
 match brep.geom.surfaces.get(si) {
 Some(rcad_kernel::geom::Surface3::Plane(_)) => {}
 _ => return false,
 }
 }
 true
}

/// True when every face normal is (approximately) 鍗, 鍗, or 鍗 in world space.
///
/// Gates post-intersection [`orthogonal_face_fuse::remove_axis_coplanar_redundant_child_faces`]:
/// True when every face normal is (approximately) +/-X, +/-Y, or +/-Z in world space.
/// Used to gate axis-aligned optimization: for two world-axis-aligned planar solids,
/// **smaller** 2D bbox on a shared plane, but in nested **box閳倻ox** the smaller patch is often the
/// true external face and the larger one is the untrimmed remainder 閳?yielding too-low
/// [`rcad_kernel::surface_area`] (OCCT `bcommon_simple/B1`). Rotated operands (`bcommon_simple/C8`)
/// still need the cleanup, so we only skip when **both** sides satisfy this predicate.
fn brep_is_world_axis_aligned_plane_solid(brep: &BRep) -> bool {
 let is_axis_unit = |n: glam::DVec3| -> bool {
 let n = n.normalize_or_zero();
 if n.length_squared() < tolerance::TOLERANCE_VEC_SQ_MIN {
 return false;
 }
 let ae = tolerance::TOLERANCE_AXIS_ALIGN;
 (n.x.abs() - 1.0).abs() < ae && n.y.abs() < ae && n.z.abs() < ae
 || (n.y.abs() - 1.0).abs() < ae && n.x.abs() < ae && n.z.abs() < ae
 || (n.z.abs() - 1.0).abs() < ae && n.x.abs() < ae && n.y.abs() < ae
 };
 if !brep_is_pure_plane_solid(brep) {
 return false;
 }
 for solid in &brep.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 if !is_axis_unit(face.normal) {
 return false;
 }
 }
 }
 }
 true
}

fn boolean_difference_empty_coincident(a: &BRep, b: &BRep) -> bool {
 if brep_shell_face_count(a) != brep_shell_face_count(b) {
 return false;
 }
 let Some([amin, amax]) = a.bounding_box() else {
 return false;
 };
 let Some([bmin, bmax]) = b.bounding_box() else {
 return false;
 };
 let scale = (amax - amin).length().max((bmax - bmin).length()).max(1.0);
 let tol = tolerance::TOLERANCE_ABS.max(tolerance::TOLERANCE_LEN_MIN * scale);
 if (amin - bmin).length() > tol || (amax - bmax).length() > tol {
 return false;
 }
 // Bbox + face count is not sufficient 閳?an inscribed rotated box shares
 // the same bbox as its container (e.g. bopcut_simple/F5).
 // Also check that vertex sets match (identical shapes have identical vertices).
 if a.vertices.len() != b.vertices.len() {
 return false;
 }
 let a_pts: Vec<glam::DVec3> = a.vertices.iter().map(|v| v.point).collect();
 let b_pts: Vec<glam::DVec3> = b.vertices.iter().map(|v| v.point).collect();
 a_pts.iter().all(|pa| b_pts.iter().any(|pb| (pa - pb).length() <= tol))
 && b_pts.iter().all(|pb| a_pts.iter().any(|pa| (pa - pb).length() <= tol))
}

fn intersection_planar_sliver_should_be_empty(result: &BRep, a: &BRep, b: &BRep) -> bool {
 let nf = face_count_of(result);
 if nf == 0 {
 return true;
 }
 let Some([amin, amax]) = a.bounding_box() else {
 return false;
 };
 let Some([bmin, bmax]) = b.bounding_box() else {
 return false;
 };
 let scale = (amax - amin).length().max((bmax - bmin).length()).max(1.0);
 let vol_tol = tolerance::TOLERANCE_VOL_CUBE_FACTOR * scale * scale * scale;
 let vol = rcad_kernel::properties::volume(result);
 if !vol.is_finite() || vol > vol_tol {
 return false;
 }

 // Require one surface slot per face and all planes 閳?`Iterator::all` on an empty iterator is
 // `true`, and skipping `None` slots could wrongly classify incomplete geom as 閳ユ竵ll planes閳?
 if result.geom.face_surface.len() != nf {
 return false;
 }
 for slot in &result.geom.face_surface {
 let Some(si) = *slot else {
 return false;
 };
 match result.geom.surfaces.get(si) {
 Some(rcad_kernel::geom::Surface3::Plane(_)) => {}
 _ => return false,
 }
 }
 true
}

/// Check if an intersection result is degenerate: all faces planar and
/// all vertices co-planar (zero thickness), meaning the solids only touch
/// at a face without volumetric overlap.
fn intersection_result_is_degenerate_sliver(result: &BRep) -> bool {
 let nf = result.solids.iter().flat_map(|s| s.shells.iter()).flat_map(|sh| sh.faces.iter()).count();
 if nf == 0 { return false; }
 // All faces must be planar
 if result.geom.face_surface.len() < nf { return false; }
 for slot in result.geom.face_surface.iter().take(nf) {
 let Some(si) = *slot else { return false };
 match result.geom.surfaces.get(si) {
 Some(rcad_kernel::geom::Surface3::Plane(_)) => {}
 _ => return false,
 }
 }
 // Check all vertices are co-planar
 let verts: Vec<glam::DVec3> = result.vertices.iter().map(|v| v.point).collect();
 if verts.len() < 3 { return false; }
 // Find first 3 non-collinear vertices to define a reference plane
 let mut ref_normal = glam::DVec3::ZERO;
 'outer: for i in 1..verts.len() {
 let d1 = verts[i] - verts[0];
 if d1.length_squared() < 1e-20 { continue; }
 for j in (i + 1)..verts.len() {
 let d2 = verts[j] - verts[0];
 let n = d1.cross(d2);
 if n.length_squared() > 1e-20 {
 ref_normal = n.normalize();
 break 'outer;
 }
 }
 }
 if ref_normal.length_squared() < 0.5 { return false; }
 let tol = tolerance::TOLERANCE_ABS;
 for v in &verts {
 let d = (*v - verts[0]).dot(ref_normal);
 if d.abs() > tol { return false; }
 }
 true
}

/// Plane recompute and planar-intersection cleanup after [`builder::BooleanBuilder::build`],
/// matching [`boolean_op_pave_fill_build`].
pub(crate) fn boolean_postprocess_pave_result(
 _op: BooleanOpType,
 _a: &BRep,
 _b: &BRep,
 result: BRep,
) -> Result<BRep, BooleanError> {
 // rcad-specific post-processing removed: recompute_plane_surfaces, coplanar clipping,
 // redundant face removal, spurious face removal, and degenerate sliver detection.
 // OCCT does not perform any post-processing after the Builder.
 // If these were needed, the root cause is in the Builder/PaveFiller pipeline.
 if !result.solids.is_empty() && !result.solids[0].shells.is_empty() {
 eprintln!("Post-process result: {} faces", result.solids[0].shells[0].faces.len());
 if std::env::var("RCAD_DEBUG_RESULT_FACES").is_ok() {
 let mut flat_idx = 0usize;
 for solid in &result.solids {
 for shell in &solid.shells {
 for face in &shell.faces {
 let surf_name = result
 .geom
 .face_surface
 .get(flat_idx)
 .and_then(|entry| *entry)
 .and_then(|surface_idx| result.geom.surfaces.get(surface_idx))
 .map(|surface| match surface {
 rcad_kernel::geom::Surface3::Plane(_) => "Plane",
 rcad_kernel::geom::Surface3::Cylinder(_) => "Cylinder",
 rcad_kernel::geom::Surface3::Cone(_) => "Cone",
 rcad_kernel::geom::Surface3::Sphere(_) => "Sphere",
 rcad_kernel::geom::Surface3::Torus(_) => "Torus",
 rcad_kernel::geom::Surface3::BSpline(_) => "BSpline",
 _ => "Other",
 })
 .unwrap_or("None");
 let area = rcad_kernel::properties::face_surface_area(&result, face, flat_idx);
 let uv_range = result
 .geom
 .face_surface_range
 .get(flat_idx)
 .and_then(|entry| *entry)
 .map(|[u0, u1, v0, v1]| {
 format!(" uv=[{u0:.4},{u1:.4}]x[{v0:.4},{v1:.4}]")
 })
 .unwrap_or_default();
 let sample = face
 .sample_point
 .map(|p| format!(" sample=({:.4},{:.4},{:.4})", p.x, p.y, p.z))
 .unwrap_or_default();
 eprintln!(
 "[RESULT_FACE] face[{flat_idx}] surf={surf_name} area={area:.6} outer_edges={} inner_wires={} tris={}{}{}",
 face.outer_wire.edges.len(),
 face.inner_wires.len(),
 face.triangles.len(),
 uv_range,
 sample,
 );
 flat_idx += 1;
 }
 }
 }
 }
 }
 Ok(result)
}

/// Topods variant: no-op post-process, returns result as-is (debug logging removed).
/// The old function only did debug logging 閳?this variant skips it.
pub(crate) fn boolean_postprocess_pave_result_topods(
 _op: BooleanOpType, _a: &BRep, _b: &BRep,
 result: topods::BRep,
) -> Result<topods::BRep, BooleanError> {
 Ok(result)
}

/// DS 閳?[`pave_filler::PaveFiller`] 閳?[`builder::BooleanBuilder`] 閳?plane surface recompute.
///
/// Used internally when a coaxial shortcut must call difference without re-entering other coaxial
/// difference branches (e.g. cylinder 閳?loft frustum after `cone 閳?cylinder`).
/// Direct PaveFiller + BooleanBuilder pipeline, no post-processing.
/// OCCT-aligned: BOPAlgo_BOP::Perform.
pub fn boolean_op(op: BooleanOpType, a: &BRep, b: &BRep) -> Result<topods::BRep, BooleanError> {
 boolean_op_pave_fill_build(op, a, b)
}

pub(crate) fn boolean_op_pave_fill_build(op: BooleanOpType, a: &BRep, b: &BRep) -> Result<topods::BRep, BooleanError> {
 let a_t = a.to_topods();
 let b_t = b.to_topods();
 let mut ds = bopds::ds::DS::new_from_topods(&a_t, &b_t, crate::tolerance::TOLERANCE_ABS);
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

 let builder = builder::BooleanBuilder::with_brep(&ds, op, brep, face_refs, ic_edge_map);
 let (result, _history) = builder.build_with_history_topods()?;
 Ok(result)
}
