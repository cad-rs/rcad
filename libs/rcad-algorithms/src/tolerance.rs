//! Geometric tolerances, adaptive context, and small predicate helpers.
//!
//! # Roadmap
//!
//! Phases B–D: repository file `docs/tolerance-system-improvement-plan.md` (workspace root).
//!
//! # Phase C — numerical IntSS floor & mesh merge
//!
//! Use [`intss_geom_tol_floor`] / [`intss_geom_tol_floor_for_brep_bounds`] when passing
//! [`geom_tol_floor`](crate::inttools::intss::intersect_surfaces_with_density_tol) without a full
//! [`ToleranceContext`]. When fuzz / workspace is relevant, reuse [`combined_linear_tol_models`] or
//! [`combined_linear_tol_for_faces`] — same OCCT-style `Max` semantics.
//!
//! Triangle-soup chaining: [`tessellation_merge_linear_from_brep`] / [`tessellation_merge_linear_from_two_breps`].
//! UV trim closure (phase C): scale with face UV bbox via [`uv_polyline_trim_closed_len_sq_from_uv_poly`];
//! ceiling matches historic [`UV_POLYLINE_TRIM_LEGACY_SQ`] linear-equivalent slack.
//!
//! # Phase B — [`ToleranceContext`]
//!
//! Use [`ToleranceContext`] to bundle scale-aware [`AdaptiveTolerance`] with an optional
//! workspace linear band (e.g. boolean fuzzy). Pairwise BRep work tolerances use
//! [`combined_linear_tol_for_faces`], [`combined_linear_tol_for_edges`],
//! [`combined_linear_tol_for_vertices`], or [`combined_linear_tol_models`] (OCCT-style `max`).
//! Point-in-solid **`classify`** folds [`DSFace::geom_tol`](crate::bopds::ds::DSFace) into relaxed
//! thresholds using [`effective_linear_with_geom_tol`]. Binary **`boolean_op_with_options`** raises
//! glue / make-connected toward [`combined_linear_tol_models`] for the paired operands (`lib.rs`).
//!
//! # Constant taxonomy (phase A)
//!
//! Pick **one row** for the kind of check; do not mix unrelated symbols.
//!
//! | Family | Primary symbols | Typical use |
//! |--------|-----------------|-------------|
//! | Point / length (strict) | [`TOLERANCE_ABS`], [`TOLERANCE_ABS_SQ`] | Coincidence, param equality in analytic code, kernel-scale residuals |
//! | Point / length (mesh / merge) | [`TOLERANCE_MESH_LEGACY`], [`TOLERANCE_PARAM_LEGACY`], [`UV_POLYLINE_TRIM_LEGACY_SQ`], [`UV_TRIM_CLOSED_REL_DOM`], [`uv_polyline_trim_closed_len_sq_from_uv_poly`] | Triangle merge, UV legs, legacy trim bounds; scaled UV closure (phase C) |
//! | Point / length (slacks) | [`TOLERANCE_PLANE_DIST_RELAX`], [`TOLERANCE_COORD_SUB`], [`TOLERANCE_LINEAR_RELAX_8`] | Coplanar slack, test nudges, one-decade-relaxed linear |
//! | Angle / direction (numeric) | [`TOLERANCE_ANG`], [`vectors_parallel`] | Parallel axes, `sin(angle)` near , cross-product magnitude in **intersection / construction** (FP noise) |
//! | Angle (heuristic radians) | [`TOLERANCE_ANG_HEURISTIC_RAD`] | Coarse “same cone / same domain” angles in boolean-style code **not** tied to cross magnitude |
//! | Direction dot heuristics | [`TOLERANCE_DOT_NEARLY_PARALLEL`] | Almost-parallel **unit** normals via `n1·n2` |
//! | Boolean fuzzy ladder | [`TOLERANCE_RETRY_LADDER_MID`], [`TOLERANCE_RETRY_LADDER_COARSE`], [`TOLERANCE_AREA_REL`] | Robust BOP retry enlargement steps |
//! | Underflow / numerical guards | [`TOLERANCE_VEC_SQ_MIN`], [`TOLERANCE_LEN_SQ_DIV_SAFE`], [`TOLERANCE_FLOAT_DEDUP`] | Degenerate normalization, dedup of ladder values |
//! | Adaptive scaling | [`AdaptiveTolerance`], [`ToleranceContext`], [`ToleranceLevel`], **`combined_linear_tol_*`**, **[`effective_linear_with_geom_tol`]**, **`intss_geom_tol_floor*`**, **`tessellation_merge_linear_*`** | Scale + workspace; pairwise `max`; classify / IntSS floors; mesh chaining |
//!
//! # `rcad_kernel::ANGULAR` vs [`TOLERANCE_ANG`] (phase A)
//!
//! | | `rcad_kernel::tolerance::ANGULAR` (~1e-12 rad) | [`TOLERANCE_ANG`] (~1e-9, cross-based) |
//! |-|------------------|----------------|
//! | Role | OCCT **nominal** angular confusion | **Algorithms** slack for accumulated FP error |
//! | Prefer when | New kernel-level topology invariants (rare in `rcad-algorithms`) | `inttools`, marching, analytic intersection branch tests, [`vectors_parallel`] |
//!
//! **Rule of thumb:** inside `rcad-algorithms`, keep using [`TOLERANCE_ANG`] / [`vectors_parallel`] unless you are intentionally matching kernel invariants and import `ANGULAR` explicitly.
//! Third kind: [`TOLERANCE_ANG_HEURISTIC_RAD`] (~1e-6 rad) for coarse angle **differences** in param space—do not substitute for [`TOLERANCE_ANG`] on cross/sin tests.
//!
//! Files that heavily use [`TOLERANCE_ANG`] (audit when tightening angles): `inttools/*.rs` (cone/cylinder/torus/plane intersections), `int_ana.rs`, `shape_construct.rs`, `classify.rs`.
//!
//! # [`AdaptiveTolerance`] vs bare constants (phase A)
//!
//! - Use **[`AdaptiveTolerance::from_brep`] / [`from_two_breps`] / [`from_scale`]** when the caller has a **BRep or known model extent** and compares world-space distances (classification, booleans, mesh eps derived from model).
//! - Use **bare [`TOLERANCE_ABS`] / [`TOLERANCE_MESH_LEGACY`]** only for **local analytic** code paths with no scale context (pure `Surface3` intersection in unit-agnostic math, unit tests, leaf predicates).
//! - Prefer **[`max_face_tolerance_or_abs`] / [`max_face_tolerance_or_abs_pair`]** when epsilon must reflect **stored face tolerances** on a BRep (e.g. soup intersection) without scale context.
//! - Prefer **[`tessellation_merge_linear_from_two_breps`]** (or [`tessellation_merge_linear_from_brep`]) when chaining mesh segments and you want **Relaxed adaptive + `TOLERANCE_MESH_LEGACY` minimum + `model_tolerance`**.
//! - Use **[`ToleranceContext`]** when you need both [`AdaptiveTolerance`] and an optional **workspace linear floor** (e.g. user boolean fuzzy); pair with [`combined_linear_tol_for_faces`] for OCCT-style `max` chains.
//! - [`AdaptiveTolerance::angular_tolerance`] intentionally scales [`TOLERANCE_ANG`] by [`ToleranceLevel`], not kernel `ANGULAR`.

use glam::{DVec2, DVec3};
use tracing::debug;

// re-export kernel Precision constants and tolerance helpers.
pub use rcad_kernel::{
 ANGULAR, APPROXIMATION, COMPUTATIONAL, CONFUSION, edge_tolerance as kernel_edge_tolerance,
 face_tolerance as kernel_face_tolerance, INFINITE_VALUE, INTERSECTION, model_tolerance,
 SQUARE_COMPUTATIONAL, SQUARE_CONFUSION, topods, vertex_tolerance as kernel_vertex_tolerance,
};

/// Absolute tolerance for point coincidence.
///
/// Matches `rcad_kernel::tolerance::CONFUSION` = `Precision::Confusion()` in OCCT.
/// Two points are considered coincident when their distance is below this value.
pub const TOLERANCE_ABS: f64 = 1e-7;

/// Angular tolerance for parallel/perpendicular checks (radians, as cross-product magnitude).
///
/// This is intentionally **looser** than `rcad_kernel::tolerance::ANGULAR` (1e-12):
/// the algorithms layer needs to tolerate slightly imperfect parallelism that
/// arises from floating-point accumulation during intersection computation.
/// Used in [`vectors_parallel`] as `cross(a,b).length_squared() < TOLERANCE_ANG²`.
pub const TOLERANCE_ANG: f64 = 1e-9;

/// Tolerance squared — avoids `sqrt` in distance checks.
pub const TOLERANCE_ABS_SQ: f64 = TOLERANCE_ABS * TOLERANCE_ABS;

/// Default pair/merge epsilon for mesh / triangle-soup paths that historically used `1e-6`.
///
/// Kept as `10 × TOLERANCE_ABS` so all legacy `1e-6` defaults track the kernel confusion value.
pub const TOLERANCE_MESH_LEGACY: f64 = TOLERANCE_ABS * 10.0;

/// Floor for chained-segment merge and similar clamps (avoid exact zero).
pub const TOLERANCE_CLAMP_MIN: f64 = 1e-15;

/// Float dedup / “same ladder value” checks in robust boolean retry (`1e-15`).
pub const TOLERANCE_FLOAT_DEDUP: f64 = 1e-15;

/// Looser float equality for coarse regression asserts (`1e-14`).
pub const TOLERANCE_FLOAT_LOOSE: f64 = 1e-14;

/// Ultra-tight float compare for high-dynamic-range test assertions (`1e-18`).
pub const TOLERANCE_FLOAT_ULTRA: f64 = 1e-18;

/// Squared metric floor treated as vanishing in area-like heuristics (`1e-20`).
pub const TOLERANCE_METRIC_SQ_NEAR_ZERO: f64 = TOLERANCE_FLOAT_ULTRA * 0.01;

// ── Legacy numeric tiers (all derivable from `TOLERANCE_ABS` at unit scale) ─────────────────────

/// Parametric / UV span and linear merge heuristics that historically used `1e-6`.
pub const TOLERANCE_PARAM_LEGACY: f64 = TOLERANCE_MESH_LEGACY;

/// Legacy UV trim loop test upper bound (**squared**) used as the **ceiling** for phase-C scaling.
///
/// Historically **`(Δ).length_squared() < Self`** with [`TOLERANCE_MESH_LEGACY`] numerically (**`~√Self`**
/// effective linear slack). Prefer [`uv_polyline_trim_closed_len_sq_from_uv_poly`] for imprint/builder trims.
pub const UV_POLYLINE_TRIM_LEGACY_SQ: f64 = TOLERANCE_MESH_LEGACY;

/// World-axis normal alignment (`|n·e|-1|`, off-diagonal) historically `1e-6`.
pub const TOLERANCE_AXIS_ALIGN: f64 = TOLERANCE_MESH_LEGACY;

/// Merge / same-domain angle heuristic in **radians** (not [`TOLERANCE_ANG`]).
pub const TOLERANCE_ANG_HEURISTIC_RAD: f64 = 1e-6;

/// Coplanar point–plane distance slack (`1e-5`); `50 × TOLERANCE_ABS`.
pub const TOLERANCE_PLANE_DIST_RELAX: f64 = TOLERANCE_ABS * 50.0;

/// Relative area / fraction comparisons (`1e-4`); `1000 × TOLERANCE_ABS`.
pub const TOLERANCE_AREA_REL: f64 = TOLERANCE_ABS * 1000.0;

/// Boolean fuzzy retry ladder mid-step (`1e-5`); `100 × TOLERANCE_ABS`.
pub const TOLERANCE_RETRY_LADDER_MID: f64 = TOLERANCE_ABS * 100.0;

/// Boolean fuzzy retry ladder coarse step (`1e-4`); [`TOLERANCE_AREA_REL`].
pub const TOLERANCE_RETRY_LADDER_COARSE: f64 = TOLERANCE_AREA_REL;

/// Planar sliver volume threshold = factor × scale³ (historically `1e-9`).
pub const TOLERANCE_VOL_CUBE_FACTOR: f64 = 1e-9;

/// Extra-strict linear residual (e.g. spurious face removal), historically `1e-10`.
pub const TOLERANCE_LINEAR_ULTRA_STRICT: f64 = 1e-10;

/// Near-zero **squared** length before treating a direction as degenerate (`1e-24`).
pub const TOLERANCE_VEC_SQ_MIN: f64 = 1e-24;

/// Squared length floor for normalization / AABB underflow guards (`1e-30`).
pub const TOLERANCE_LEN_SQ_DIV_SAFE: f64 = 1e-30;

/// Ultra-tight squared-length guard for UV-space direction checks (`1e-40`).
///
/// Stricter than [`TOLERANCE_LEN_SQ_DIV_SAFE`] — used when even trace
/// squared-length values near `1e-30` are still degenerate in 2D param space
/// (e.g. UV direction calculations in wire splitting).
pub const TOLERANCE_UV_DIR_SQ_MIN: f64 = 1e-40;

/// Absolute floor for squared tangent length before treating it as zero (`1e-60`).
///
/// Equivalently `TOLERANCE_LEN_SQ_DIV_SAFE²`. Used in tangent-vector underflow
/// guards in edge-edge interference detection, where any floating noise below
/// this threshold should be treated as zero.
pub const TOLERANCE_TANGENT_SQ_ABSOLUTE_MIN: f64 = 1e-60;

/// Cosine threshold for "near parallel" line-line direction check (`0.9063 ≈ cos 25°`).
///
/// When two edges (or their tangents) have a direction dot product at or above
/// this value, they are considered near-parallel enough to warrant additional
/// fuzzy tolerance in edge-edge interference detection.
pub const TOLERANCE_COS_LINE_ANGLE: f64 = 0.9063;

/// Minimum edge / segment length floor (`1e-12`).
pub const TOLERANCE_LEN_MIN: f64 = 1e-12;

/// Tiny coordinate offset (`1e-9` = `0.01 × TOLERANCE_ABS`) for near-degenerate geometry tests.
pub const TOLERANCE_COORD_SUB: f64 = TOLERANCE_ABS * 0.01;

/// Linear epsilon one decade looser than [`TOLERANCE_ABS`] (`1e-8` at default kernel scale).
pub const TOLERANCE_LINEAR_RELAX_8: f64 = TOLERANCE_ABS * 0.1;

/// Lower clamp for adaptive `model_scale` (`1e-10`).
pub const TOLERANCE_MODEL_SCALE_MIN: f64 = 1e-10;

/// Normals “almost parallel” dot-product bound (historically `0.999`).
pub const TOLERANCE_DOT_NEARLY_PARALLEL: f64 = 0.999;

/// Dimensionless factor: `tol * k` sub-epsilons (historically `1e-6` of a base tolerance).
pub const TOLERANCE_TOL_SCALE_MICRO: f64 = 1e-6;

/// Phase C UV trim closure: fractional slack **`max(|Δu|,|Δv|)` · Self** vs [`UV_POLYLINE_TRIM_LEGACY_SQ`] ceiling (`imprint`, `builder`).
pub const UV_TRIM_CLOSED_REL_DOM: f64 = TOLERANCE_TOL_SCALE_MICRO;

/// Upper clamp in [`AdaptiveTolerance::default`] (`1e-3`).
pub const TOLERANCE_ADAPTIVE_MAX: f64 = 1e-3;

/// Coarse axis / corner slack for classification (`2e-2`).
pub const TOLERANCE_AXIS_CORNER_SLACK: f64 = 2e-2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToleranceLevel {
 /// Strict tolerance for high-precision operations (e.g., intersection points).
 /// Scale factor: 1.0
 Strict,
 /// Normal tolerance for general operations (e.g., point classification).
 /// Scale factor: 10.0
 Normal,
 /// Relaxed tolerance for approximate operations (e.g., bounding box checks).
 /// Scale factor: 100.0
 Relaxed,
 /// Very relaxed tolerance for coarse operations (e.g., AABB pre-filter).
 /// Scale factor: 1000.0
 Coarse,
}

impl ToleranceLevel {
 /// Get the scale factor for this tolerance level.
 pub fn scale_factor(self) -> f64 {
 match self {
 ToleranceLevel::Strict => 1.0,
 ToleranceLevel::Normal => 10.0,
 ToleranceLevel::Relaxed => 100.0,
 ToleranceLevel::Coarse => 1000.0,
 }
 }
}

/// Adaptive tolerance context based on model scale.
///
/// Instead of using hard-coded absolute tolerances, this context computes
/// tolerances relative to the model's bounding box size. This ensures that
/// models at different scales (e.g., nanometer vs kilometer) are handled
/// appropriately.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveTolerance {
 /// Base tolerance (typically TOLERANCE_ABS for unit-scale models).
 pub base_tolerance: f64,
 /// Model scale factor (e.g., bounding box diagonal).
 pub model_scale: f64,
 /// Minimum tolerance to prevent excessive precision requirements.
 pub min_tolerance: f64,
 /// Maximum tolerance to prevent excessive looseness.
 pub max_tolerance: f64,
}

impl Default for AdaptiveTolerance {
 fn default() -> Self {
 Self {
 base_tolerance: TOLERANCE_ABS,
 model_scale: 1.0,
 min_tolerance: TOLERANCE_LEN_MIN,
 max_tolerance: TOLERANCE_ADAPTIVE_MAX,
 }
 }
}

impl AdaptiveTolerance {
 /// Create a new adaptive tolerance with default base tolerance.
 pub fn new() -> Self {
 Self::default()
 }

 /// Create a new adaptive tolerance from a BRep's bounding box.
 pub fn from_brep(brep: &rcad_kernel::BRep) -> Self {
 let scale = compute_model_scale(brep);
 Self::from_scale(scale)
 }

 /// Create a new adaptive tolerance from two BReps' combined bounding box.
 pub fn from_two_breps(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Self {
 let scale_a = compute_model_scale(a);
 let scale_b = compute_model_scale(b);
 Self::from_scale(scale_a.max(scale_b))
 }

 /// Create a new adaptive tolerance from a known scale.
 pub fn from_scale(model_scale: f64) -> Self {
 let mut ctx = Self::default();
 ctx.model_scale = model_scale.max(TOLERANCE_MODEL_SCALE_MIN);
 ctx
 }

 /// Get the effective tolerance for a specific level.
 pub fn tolerance(self, level: ToleranceLevel) -> f64 {
 let raw = self.base_tolerance * level.scale_factor() * self.model_scale;
 raw.clamp(self.min_tolerance, self.max_tolerance)
 }

 /// Get the squared tolerance for a specific level.
 pub fn tolerance_sq(self, level: ToleranceLevel) -> f64 {
 let t = self.tolerance(level);
 t * t
 }

 /// Get the angular tolerance (not affected by model scale).
 pub fn angular_tolerance(self, level: ToleranceLevel) -> f64 {
 TOLERANCE_ANG * level.scale_factor()
 }

 /// Check if two points coincide at the given tolerance level.
 pub fn points_coincide_at(self, a: DVec3, b: DVec3, level: ToleranceLevel) -> bool {
 (a - b).length_squared() < self.tolerance_sq(level)
 }

 /// Check if a vector is zero at the given tolerance level.
 pub fn is_zero_vec_at(self, v: DVec3, level: ToleranceLevel) -> bool {
 v.length_squared() < self.tolerance_sq(level)
 }

 /// Check if two parameters are equal at the given tolerance level.
 pub fn params_equal_at(self, a: f64, b: f64, level: ToleranceLevel) -> bool {
 (a - b).abs() < self.tolerance(level)
 }

 /// Get tolerance for point coincidence (strict level).
 pub fn coincidence(self) -> f64 {
 self.tolerance(ToleranceLevel::Strict)
 }

 /// Get tolerance for classification operations (normal level).
 pub fn classification(self) -> f64 {
 self.tolerance(ToleranceLevel::Normal)
 }

 /// Get tolerance for boundary checks (relaxed level).
 pub fn boundary(self) -> f64 {
 self.tolerance(ToleranceLevel::Relaxed)
 }

 /// Get tolerance for AABB pre-filtering (coarse level).
 pub fn coarse(self) -> f64 {
 self.tolerance(ToleranceLevel::Coarse)
 }
}

/// Scale-aware tolerances plus an optional workspace linear band (phase B).
///
/// OCCT-style pairing with stored shape tolerances: use [`combined_linear_tol_for_faces`] or
/// [`combined_linear_tol_models`] to form `max(workspace, adaptive, Tol(…))`.
#[derive(Debug, Clone, Copy)]
pub struct ToleranceContext {
 /// Characteristic-size–scaled tolerances.
 pub adaptive: AdaptiveTolerance,
 /// Extra linear band (e.g. user boolean fuzzy). Negative values are clamped in accessors.
 ///
 /// [`Self::workspace_linear`] returns `max(adaptive.tolerance(level), workspace_fuzzy)`.
 pub workspace_fuzzy: f64,
}

impl Default for ToleranceContext {
 fn default() -> Self {
 Self {
 adaptive: AdaptiveTolerance::default(),
 workspace_fuzzy: 0.0,
 }
 }
}

impl ToleranceContext {
 pub fn new(adaptive: AdaptiveTolerance, workspace_fuzzy: f64) -> Self {
 Self {
 adaptive,
 workspace_fuzzy,
 }
 }

 pub fn from_adaptive(adaptive: AdaptiveTolerance) -> Self {
 Self::new(adaptive, 0.0)
 }

 pub fn from_brep(brep: &rcad_kernel::BRep) -> Self {
 Self::from_adaptive(AdaptiveTolerance::from_brep(brep))
 }

 pub fn from_two_breps(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Self {
 Self::from_adaptive(AdaptiveTolerance::from_two_breps(a, b))
 }

 pub fn from_scale(model_scale: f64) -> Self {
 Self::from_adaptive(AdaptiveTolerance::from_scale(model_scale))
 }

 #[inline]
 fn wf_nonnegative(&self) -> f64 {
 self.workspace_fuzzy.max(0.0)
 }

 /// Linear tolerance from [`AdaptiveTolerance`] only (ignores [`Self::workspace_fuzzy`]).
 #[inline]
 pub fn adaptive_linear(&self, level: ToleranceLevel) -> f64 {
 self.adaptive.tolerance(level)
 }

 /// `max(adaptive_linear(level), workspace_fuzzy)` — for fuzzy-capable pipelines (booleans).
 #[inline]
 pub fn workspace_linear(&self, level: ToleranceLevel) -> f64 {
 let adaptive = self.adaptive.tolerance(level);
 let wf = self.wf_nonnegative();
 let result = adaptive.max(wf);
 if wf > adaptive {
 debug!(
 "ToleranceContext::workspace_linear: workspace fuzzy ({:.3e}) dominates adaptive ({:.3e}) at level {:?} → {:.3e}",
 wf, adaptive, level, result,
 );
 }
 result
 }

 /// Angular tolerance at `level` (algorithms-layer [`TOLERANCE_ANG`] scaling).
 #[inline]
 pub fn angular(&self, level: ToleranceLevel) -> f64 {
 self.adaptive.angular_tolerance(level)
 }
}

/// OCCT-style `max(base, Tol(face_a), Tol(face_b))` with optional operands.
///
/// `base` is [`ToleranceContext::workspace_linear`] when `use_workspace` is true, else
/// [`ToleranceContext::adaptive_linear`]. A missing face uses [`CONFUSION`] for that side only.
#[inline]
pub fn combined_linear_tol_for_faces(
 ctx: &ToleranceContext,
 level: ToleranceLevel,
 use_workspace: bool,
 brep_a: &rcad_kernel::BRep,
 face_a: Option<usize>,
 brep_b: &rcad_kernel::BRep,
 face_b: Option<usize>,
) -> f64 {
 let base = if use_workspace {
 ctx.workspace_linear(level)
 } else {
 ctx.adaptive_linear(level)
 };
 let ta = face_a
 .map(|i| kernel_face_tolerance(brep_a, i))
 .unwrap_or(CONFUSION);
 let tb = face_b
 .map(|i| kernel_face_tolerance(brep_b, i))
 .unwrap_or(CONFUSION);
 base.max(ta).max(tb)
}

/// `max(base, model_tolerance(a), model_tolerance(b))` where `base` is chosen like
/// [`combined_linear_tol_for_faces`].
#[inline]
pub fn combined_linear_tol_models(
 ctx: &ToleranceContext,
 level: ToleranceLevel,
 use_workspace: bool,
 brep_a: &rcad_kernel::BRep,
 brep_b: &rcad_kernel::BRep,
) -> f64 {
 let base = if use_workspace {
 ctx.workspace_linear(level)
 } else {
 ctx.adaptive_linear(level)
 };
 base.max(model_tolerance(brep_a)).max(model_tolerance(brep_b))
}

/// `max(base_linear, geom_tol)` for checks that consume **stored** geometric tolerance
/// (e.g. [`crate::bopds::ds::DSFace::geom_tol`]). Non-finite `entity_geom_tol` is ignored.
#[inline]
pub fn effective_linear_with_geom_tol(base_linear: f64, entity_geom_tol: f64) -> f64 {
 if entity_geom_tol.is_finite() && entity_geom_tol > 0.0 {
 base_linear.max(entity_geom_tol)
 } else {
 base_linear
 }
}

/// Geometric tolerance **floor** for numerical surface–surface intersection (IntSS).
///
/// Combines a caller **baseline** (adaptive / workspace linear band) with the maximum stored
/// tolerance on participating entities — OCCT-style `max`, with at least [`TOLERANCE_ABS`] on the
/// baseline before folding `participant_tolerance_max`.
///
/// See [`intersect_surfaces_with_density_tol`](crate::inttools::intss::intersect_surfaces_with_density_tol).
#[inline]
pub fn intss_geom_tol_floor(baseline_linear: f64, participant_tolerance_max: f64) -> f64 {
 let base = baseline_linear.max(TOLERANCE_ABS);
 effective_linear_with_geom_tol(base, participant_tolerance_max)
}

/// [`intss_geom_tol_floor`] with `participant_tolerance_max = max(model_tol(A), model_tol(B))`.
#[inline]
pub fn intss_geom_tol_floor_for_brep_bounds(
 baseline_linear: f64,
 brep_a: &rcad_kernel::BRep,
 brep_b: &rcad_kernel::BRep,
) -> f64 {
 intss_geom_tol_floor(
 baseline_linear,
 model_tolerance(brep_a).max(model_tolerance(brep_b)),
 )
}

/// Triangle-soup pair/merge epsilon derived from **[`AdaptiveTolerance::from_brep`]** at
/// [`ToleranceLevel::Relaxed`], at least [`TOLERANCE_MESH_LEGACY`], folded with [`model_tolerance`].
#[inline]
pub fn tessellation_merge_linear_from_brep(brep: &rcad_kernel::BRep) -> f64 {
 let adaptive = AdaptiveTolerance::from_brep(brep);
 let base = adaptive.tolerance(ToleranceLevel::Relaxed).max(TOLERANCE_MESH_LEGACY);
 effective_linear_with_geom_tol(base, model_tolerance(brep))
}

/// Like [`tessellation_merge_linear_from_brep`] using [`AdaptiveTolerance::from_two_breps`] and
/// `max(model_tolerance(a), model_tolerance(b))` for the topological fold.
#[inline]
pub fn tessellation_merge_linear_from_two_breps(brep_a: &rcad_kernel::BRep, brep_b: &rcad_kernel::BRep) -> f64 {
 let adaptive = AdaptiveTolerance::from_two_breps(brep_a, brep_b);
 let base = adaptive.tolerance(ToleranceLevel::Relaxed).max(TOLERANCE_MESH_LEGACY);
 let geom = model_tolerance(brep_a).max(model_tolerance(brep_b));
 effective_linear_with_geom_tol(base, geom)
}

/// Squared UV distance threshold for treating trim endpoints as coincident given face UV extents
/// **`u_extent`**, **`v_extent`** (typically `Δu`, `Δv` of the outer UV polygon bbox).
///
/// Let **`ext = max(|u_extent|, |v_extent|)`**. Linear slack is **`min(√`**[`UV_POLYLINE_TRIM_LEGACY_SQ`]**, **`max([`TOLERANCE_ABS`]**, **`ext · [`UV_TRIM_CLOSED_REL_DOM`]`**))`** (then squared). That caps looseness at the historic mesh-tier analogue while tightening on modest UV domains (phase C).
#[inline]
pub fn uv_polyline_trim_closed_len_sq(u_extent: f64, v_extent: f64) -> f64 {
 let ext = u_extent.abs().max(v_extent.abs()).max(TOLERANCE_MODEL_SCALE_MIN);
 let rel_lin = ext * UV_TRIM_CLOSED_REL_DOM;
 let cap_lin = UV_POLYLINE_TRIM_LEGACY_SQ.sqrt();
 let lin = rel_lin.max(TOLERANCE_ABS).min(cap_lin);
 lin * lin
}

/// [`uv_polyline_trim_closed_len_sq`] from the axis-aligned bbox of **`poly`** in UV; falls back to [`UV_POLYLINE_TRIM_LEGACY_SQ`]
/// when **`poly`** is empty.
#[inline]
pub fn uv_polyline_trim_closed_len_sq_from_uv_poly(poly: &[DVec2]) -> f64 {
 if poly.is_empty() {
 return UV_POLYLINE_TRIM_LEGACY_SQ;
 }
 let mut u_min = f64::INFINITY;
 let mut u_max = f64::NEG_INFINITY;
 let mut v_min = f64::INFINITY;
 let mut v_max = f64::NEG_INFINITY;
 for p in poly {
 u_min = u_min.min(p.x);
 u_max = u_max.max(p.x);
 v_min = v_min.min(p.y);
 v_max = v_max.max(p.y);
 }
 uv_polyline_trim_closed_len_sq(u_max - u_min, v_max - v_min)
}

/// OCCT-style `max(base, Tol(edge_a), Tol(edge_b))` with flattened edge indices.
#[inline]
pub fn combined_linear_tol_for_edges(
 ctx: &ToleranceContext,
 level: ToleranceLevel,
 use_workspace: bool,
 brep_a: &rcad_kernel::BRep,
 edge_a: Option<usize>,
 brep_b: &rcad_kernel::BRep,
 edge_b: Option<usize>,
) -> f64 {
 let base = if use_workspace {
 ctx.workspace_linear(level)
 } else {
 ctx.adaptive_linear(level)
 };
 let ta = edge_a
 .map(|i| kernel_edge_tolerance(brep_a, i))
 .unwrap_or(CONFUSION);
 let tb = edge_b
 .map(|i| kernel_edge_tolerance(brep_b, i))
 .unwrap_or(CONFUSION);
 base.max(ta).max(tb)
}

/// OCCT-style `max(base, Tol(vertex_a), Tol(vertex_b))`.
#[inline]
pub fn combined_linear_tol_for_vertices(
 ctx: &ToleranceContext,
 level: ToleranceLevel,
 use_workspace: bool,
 brep_a: &rcad_kernel::BRep,
 vertex_a: Option<usize>,
 brep_b: &rcad_kernel::BRep,
 vertex_b: Option<usize>,
) -> f64 {
 let base = if use_workspace {
 ctx.workspace_linear(level)
 } else {
 ctx.adaptive_linear(level)
 };
 let ta = vertex_a
 .map(|i| kernel_vertex_tolerance(brep_a, i))
 .unwrap_or(CONFUSION);
 let tb = vertex_b
 .map(|i| kernel_vertex_tolerance(brep_b, i))
 .unwrap_or(CONFUSION);
 base.max(ta).max(tb)
}

/// Scale the default boolean fuzzy retry ladder (`10×`, `100×`, `1000×` [`TOLERANCE_ABS`]) by
/// `initial_fuzzy / TOLERANCE_ABS` (minimum scale 1).
///
/// With `extra_coarse_base = Some(c)`, appends `10·c` and `100·c` (multiscale mechanical preset).
pub fn boolean_fuzzy_ladder_scaled(initial_fuzzy: f64, extra_coarse_base: Option<f64>) -> Vec<f64> {
 let scale = (initial_fuzzy / TOLERANCE_ABS).max(1.0);
 let mut out = vec![
 TOLERANCE_ABS * 10.0 * scale,
 TOLERANCE_ABS * 100.0 * scale,
 TOLERANCE_ABS * 1000.0 * scale,
 ];
 if let Some(c) = extra_coarse_base {
 let c = c.max(0.0);
 if c > 0.0 {
 out.push(c * 10.0);
 out.push(c * 100.0);
 }
 }
 out
}

/// Compute the characteristic scale of a model from its bounding box.
/// Returns the diagonal of the bounding box, or 1.0 if the model is empty.
pub fn compute_model_scale(brep: &rcad_kernel::BRep) -> f64 {
 let mut min_pt = DVec3::splat(f64::INFINITY);
 let mut max_pt = DVec3::splat(f64::NEG_INFINITY);
 let mut has_vertices = false;

 for ts in &brep.tshapes {
 if let topods::TShape::Vertex(vd) = ts.as_ref() {
 min_pt = min_pt.min(vd.point);
 max_pt = max_pt.max(vd.point);
 has_vertices = true;
 }
 }

 if !has_vertices {
 return 1.0;
 }

 let diagonal = (max_pt - min_pt).length();
 diagonal.max(TOLERANCE_MODEL_SCALE_MIN)
}

/// Compute the characteristic scale from a collection of points.
pub fn compute_scale_from_points(points: &[DVec3]) -> f64 {
 if points.is_empty() {
 return 1.0;
 }

 let mut min_pt = DVec3::splat(f64::INFINITY);
 let mut max_pt = DVec3::splat(f64::NEG_INFINITY);

 for &p in points {
 min_pt = min_pt.min(p);
 max_pt = max_pt.max(p);
 }

 let diagonal = (max_pt - min_pt).length();
 diagonal.max(TOLERANCE_MODEL_SCALE_MIN)
}

/// Face count in BRep traversal order (matches `geom.face_tolerance` / `face_surface` flat indices).
pub fn flat_face_count(brep: &rcad_kernel::BRep) -> usize {
 let mut count = 0;
 for ts in &brep.tshapes {
 if let topods::TShape::Solid(sd) = ts.as_ref() {
 for sr in &sd.shells {
 if let topods::TShape::Shell(shd) = &*brep.tshapes[sr.index] {
 count += shd.faces.len();
 }
 }
 }
 }
 count
}

/// Maximum positive finite `geom.face_tolerance` entry for the first [`flat_face_count`] slots.
///
/// If the array is missing entries or has no valid values, returns [`TOLERANCE_ABS`].
/// Use to derive mesh / triangle-soup intersection epsilons when geometry came from this BRep.
pub fn max_face_tolerance_or_abs(brep: &rcad_kernel::BRep) -> f64 {
 let mut m = TOLERANCE_ABS;
 for ts in &brep.tshapes {
 if let topods::TShape::Face(fd) = ts.as_ref() {
 let t = fd.tolerance;
 if t.is_finite() && t > 0.0 {
 m = m.max(t);
 }
 }
 }
 m
}

/// [`max_face_tolerance_or_abs`] applied to both inputs; suitable for mesh intersection of two sources.
#[inline]
pub fn max_face_tolerance_or_abs_pair(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> f64 {
 max_face_tolerance_or_abs(a).max(max_face_tolerance_or_abs(b))
}

#[inline]
pub fn points_coincide(a: DVec3, b: DVec3) -> bool {
 (a - b).length_squared() < TOLERANCE_ABS_SQ
}

#[inline]
pub fn is_zero_vec(v: DVec3) -> bool {
 v.length_squared() < TOLERANCE_ABS_SQ
}

/// Returns true if two unit vectors are parallel (or anti-parallel).
#[inline]
pub fn vectors_parallel(a: DVec3, b: DVec3) -> bool {
 a.cross(b).length_squared() < TOLERANCE_ANG * TOLERANCE_ANG
}

#[inline]
pub fn params_equal(a: f64, b: f64) -> bool {
 (a - b).abs() < TOLERANCE_ABS
}

/// Check if two vectors are parallel using adaptive tolerance.
pub fn vectors_parallel_adaptive(a: DVec3, b: DVec3, tol: AdaptiveTolerance) -> bool {
 let ang_tol = tol.angular_tolerance(ToleranceLevel::Normal);
 a.cross(b).length_squared() < ang_tol * ang_tol
}

/// Precision::IsInfinite (Precision.hxx L350-353).
/// Returns true if `R` may be considered as an infinite number.
/// OCCT: std::abs(R) >= 0.5 * Precision::Infinite() where Precision::Infinite() = 2e100.
/// Delegates to rcad_kernel::is_infinite_value().
#[inline]
pub fn precision_is_infinite(r: f64) -> bool {
 rcad_kernel::tolerance::is_infinite_value(r)
}

/// Precision::IsPositiveInfinite (Precision.hxx L357-360).
/// Returns true if R may be considered as a positive infinite number.
#[inline]
pub fn precision_is_positive_infinite(r: f64) -> bool {
 rcad_kernel::tolerance::is_positive_infinite_value(r)
}

/// Precision::IsNegativeInfinite (Precision.hxx L364-367).
/// Returns true if R may be considered as a negative infinite number.
#[inline]
pub fn precision_is_negative_infinite(r: f64) -> bool {
 rcad_kernel::tolerance::is_negative_infinite_value(r)
}

// Re-export for test module convenience (any_perpendicular is in rcad-kernel).
pub use rcad_kernel::geom::any_perpendicular;


