//! Advanced Boolean Retry Strategy Enhancement
//!
//! This module provides detailed failure classification and targeted recovery
//! strategies for boolean operations. It builds on the basic retry mechanism
//! with more sophisticated failure detection and recovery.

use std::fmt;

use crate::tolerance::{
    TOLERANCE_ABS, TOLERANCE_ADAPTIVE_MAX, TOLERANCE_CLAMP_MIN, TOLERANCE_FLOAT_DEDUP,
    TOLERANCE_LEN_MIN, TOLERANCE_MESH_LEGACY, TOLERANCE_RETRY_LADDER_COARSE,
    TOLERANCE_RETRY_LADDER_MID,
};
use crate::{
    BRep, BooleanError, BooleanOpType, BooleanOptions, BooleanRetryClass,
    BooleanRobustOptions, BooleanRetryPolicy, ExtremeGeometryRetryConfig,
    boolean_op_robust, classify_boolean_failure,
};

/// Detailed failure classification for boolean operations.
///
/// This enum provides more specific failure types than `BooleanRetryClass`,
/// enabling targeted recovery strategies for each failure mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[derive(Default)]
pub enum BooleanFailureClass {
    /// Result has degenerate edges or faces (zero-length edges, degenerate triangles).
    DegenerateTopology,
    /// Numerical errors during computation (NaN, infinity, precision loss).
    NumericalInstability,
    /// Result fails validity checks (non-manifold, open shells, invalid orientation).
    InvalidResult,
    /// Missing intersection curves between surfaces that should intersect.
    IncompleteIntersection,
    /// Result contains self-intersecting geometry.
    SelfIntersection,
    /// Input geometry is structurally invalid (empty, missing data).
    InvalidInput,
    /// Unknown or unclassified failure.
    #[default]
    Unknown,
}

impl BooleanFailureClass {
    /// Returns a human-readable description of the failure class.
    pub fn description(&self) -> &'static str {
        match self {
            Self::DegenerateTopology => "Result contains degenerate topology (zero-length edges or degenerate faces)",
            Self::NumericalInstability => "Numerical errors during computation (NaN, infinity, or precision loss)",
            Self::InvalidResult => "Result fails validity checks (non-manifold, open shells, or invalid orientation)",
            Self::IncompleteIntersection => "Missing intersection curves between surfaces",
            Self::SelfIntersection => "Result contains self-intersecting geometry",
            Self::InvalidInput => "Input geometry is structurally invalid",
            Self::Unknown => "Unknown or unclassified failure",
        }
    }

    /// Returns whether this failure class can potentially be recovered by retry.
    pub fn is_recoverable(&self) -> bool {
        matches!(
            self,
            Self::DegenerateTopology
                | Self::NumericalInstability
                | Self::InvalidResult
                | Self::IncompleteIntersection
                | Self::SelfIntersection
        )
    }

    /// Returns the suggested recovery strategy for this failure class.
    pub fn suggested_recovery(&self) -> RecoveryStrategy {
        match self {
            Self::DegenerateTopology => RecoveryStrategy::MakeConnectedCleanup,
            Self::NumericalInstability => RecoveryStrategy::IncreaseFuzzyTolerance,
            Self::InvalidResult => RecoveryStrategy::AlgorithmVariant,
            Self::IncompleteIntersection => RecoveryStrategy::EnableGlueMode,
            Self::SelfIntersection => RecoveryStrategy::MakeConnectedCleanup,
            Self::InvalidInput => RecoveryStrategy::None,
            Self::Unknown => RecoveryStrategy::None,
        }
    }
}

impl fmt::Display for BooleanFailureClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DegenerateTopology => write!(f, "degenerate topology"),
            Self::NumericalInstability => write!(f, "numerical instability"),
            Self::InvalidResult => write!(f, "invalid result"),
            Self::IncompleteIntersection => write!(f, "incomplete intersection"),
            Self::SelfIntersection => write!(f, "self-intersection"),
            Self::InvalidInput => write!(f, "invalid input"),
            Self::Unknown => write!(f, "unknown failure"),
        }
    }
}


/// Recovery strategy to apply for a specific failure class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecoveryStrategy {
    /// No recovery possible.
    None,
    /// Run MakeConnected cleanup to fix topology issues.
    MakeConnectedCleanup,
    /// Increase fuzzy tolerance and retry.
    IncreaseFuzzyTolerance,
    /// Try a different algorithm variant.
    AlgorithmVariant,
    /// Enable Glue mode for better intersection handling.
    EnableGlueMode,
    /// Combine multiple strategies.
    Combined {
        use_glue: bool,
        run_make_connected: bool,
        increase_fuzzy: bool,
    },
}

impl RecoveryStrategy {
    /// Returns a description of this recovery strategy.
    pub fn description(&self) -> &'static str {
        match self {
            Self::None => "No recovery available",
            Self::MakeConnectedCleanup => "Run MakeConnected cleanup to fix topology",
            Self::IncreaseFuzzyTolerance => "Increase fuzzy tolerance and retry",
            Self::AlgorithmVariant => "Try different algorithm variant",
            Self::EnableGlueMode => "Enable Glue mode for intersection handling",
            Self::Combined { .. } => "Combine multiple recovery strategies",
        }
    }
}

/// Configurable retry policy for boolean operations.
///
/// This struct provides fine-grained control over retry behavior,
/// including tolerance growth, glue mode activation, and cleanup aggressiveness.
#[derive(Debug, Clone)]
pub struct RetryPolicy {
    /// Maximum number of retry attempts.
    pub max_attempts: usize,
    /// Factor by which to multiply fuzzy tolerance on each retry.
    pub fuzzy_growth_factor: f64,
    /// Maximum fuzzy tolerance cap.
    pub fuzzy_tolerance_cap: f64,
    /// Number of failures after which to enable Glue mode.
    pub enable_glue_after_n_failures: usize,
    /// Glue tolerance for shared-face detection.
    pub glue_tolerance: f64,
    /// Cleanup aggressiveness level (1-10, higher = more aggressive).
    pub make_connected_aggressiveness: u32,
    /// Maximum passes for MakeConnected cleanup.
    pub make_connected_max_passes: usize,
    /// Initial tolerance for MakeConnected cleanup.
    pub make_connected_initial_tolerance: f64,
    /// Tolerance growth factor for MakeConnected passes.
    pub make_connected_tolerance_growth: f64,
    /// Whether to use scoped MakeConnected when possible.
    pub use_scoped_make_connected: bool,
    /// Whether to fall back to global cleanup when scoped fails.
    pub fallback_to_global_cleanup: bool,
    /// Whether to try algorithm variants on InvalidResult failures.
    pub try_algorithm_variants: bool,
    /// Whether to enable verbose diagnostic output.
    pub verbose_diagnostics: bool,
}

impl Default for RetryPolicy {
    fn default() -> Self {
        Self {
            max_attempts: 5,
            fuzzy_growth_factor: 10.0,
            fuzzy_tolerance_cap: TOLERANCE_ADAPTIVE_MAX,
            enable_glue_after_n_failures: 2,
            glue_tolerance: TOLERANCE_MESH_LEGACY,
            make_connected_aggressiveness: 5,
            make_connected_max_passes: 5,
            make_connected_initial_tolerance: TOLERANCE_MESH_LEGACY,
            make_connected_tolerance_growth: 2.0,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: true,
            verbose_diagnostics: false,
        }
    }
}

impl RetryPolicy {
    /// Creates a conservative retry policy with minimal intervention.
    pub fn conservative() -> Self {
        Self {
            max_attempts: 3,
            fuzzy_growth_factor: 5.0,
            fuzzy_tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE,
            enable_glue_after_n_failures: 3,
            glue_tolerance: TOLERANCE_ABS,
            make_connected_aggressiveness: 3,
            make_connected_max_passes: 3,
            make_connected_initial_tolerance: TOLERANCE_ABS,
            make_connected_tolerance_growth: 1.5,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Creates an aggressive retry policy for difficult geometry.
    pub fn aggressive() -> Self {
        Self {
            max_attempts: 10,
            fuzzy_growth_factor: 20.0,
            fuzzy_tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE * 100.0,
            enable_glue_after_n_failures: 1,
            glue_tolerance: TOLERANCE_RETRY_LADDER_MID,
            make_connected_aggressiveness: 8,
            make_connected_max_passes: 10,
            make_connected_initial_tolerance: TOLERANCE_RETRY_LADDER_MID,
            make_connected_tolerance_growth: 3.0,
            use_scoped_make_connected: false,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: true,
            verbose_diagnostics: true,
        }
    }

    /// Creates a retry policy tuned for numerical instability cases.
    pub fn for_numerical_instability() -> Self {
        Self {
            max_attempts: 8,
            fuzzy_growth_factor: 15.0,
            fuzzy_tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE * 50.0,
            enable_glue_after_n_failures: 1,
            glue_tolerance: TOLERANCE_RETRY_LADDER_MID,
            make_connected_aggressiveness: 6,
            make_connected_max_passes: 7,
            make_connected_initial_tolerance: TOLERANCE_MESH_LEGACY,
            make_connected_tolerance_growth: 2.5,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Creates a retry policy tuned for degenerate topology cases.
    pub fn for_degenerate_topology() -> Self {
        Self {
            max_attempts: 6,
            fuzzy_growth_factor: 5.0,
            fuzzy_tolerance_cap: TOLERANCE_ADAPTIVE_MAX,
            enable_glue_after_n_failures: 2,
            glue_tolerance: TOLERANCE_MESH_LEGACY,
            make_connected_aggressiveness: 9,
            make_connected_max_passes: 10,
            make_connected_initial_tolerance: TOLERANCE_RETRY_LADDER_MID,
            make_connected_tolerance_growth: 2.0,
            use_scoped_make_connected: false,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Creates a retry policy tuned for incomplete intersection cases.
    pub fn for_incomplete_intersection() -> Self {
        Self {
            max_attempts: 6,
            fuzzy_growth_factor: 10.0,
            fuzzy_tolerance_cap: TOLERANCE_RETRY_LADDER_COARSE * 50.0,
            enable_glue_after_n_failures: 0, // Enable immediately
            glue_tolerance: TOLERANCE_RETRY_LADDER_MID,
            make_connected_aggressiveness: 5,
            make_connected_max_passes: 5,
            make_connected_initial_tolerance: TOLERANCE_MESH_LEGACY,
            make_connected_tolerance_growth: 2.0,
            use_scoped_make_connected: true,
            fallback_to_global_cleanup: true,
            try_algorithm_variants: false,
            verbose_diagnostics: false,
        }
    }

    /// Computes the next fuzzy tolerance based on the policy.
    pub fn next_fuzzy_tolerance(&self, current: f64) -> f64 {
        let next = current * self.fuzzy_growth_factor;
        next.min(self.fuzzy_tolerance_cap)
    }

    /// Determines whether glue mode should be enabled after N failures.
    pub fn should_enable_glue(&self, failure_count: usize) -> bool {
        failure_count >= self.enable_glue_after_n_failures
    }

    /// Computes the MakeConnected tolerance based on aggressiveness.
    pub fn make_connected_tolerance(&self, pass: usize) -> f64 {
        let base = self.make_connected_initial_tolerance;
        let growth = self.make_connected_tolerance_growth.powi(pass as i32);
        base * growth
    }
}

/// Diagnostic information for a single retry attempt.
#[derive(Debug, Clone, Default)]
pub struct BooleanAttemptDiagnostic {
    /// Attempt number (1-indexed).
    pub attempt_number: usize,
    /// Fuzzy tolerance used for this attempt.
    pub fuzzy_tolerance: f64,
    /// Whether glue mode was enabled.
    pub glue_enabled: bool,
    /// Glue tolerance used.
    pub glue_tolerance: f64,
    /// Whether MakeConnected cleanup was run.
    pub make_connected_run: bool,
    /// Number of MakeConnected passes.
    pub make_connected_passes: usize,
    /// Whether this attempt succeeded.
    pub success: bool,
    /// Failure class if the attempt failed.
    pub failure_class: Option<BooleanFailureClass>,
    /// Recovery strategy applied before this attempt.
    pub recovery_strategy: Option<RecoveryStrategy>,
    /// Error message if the attempt failed.
    pub error_message: Option<String>,
    /// Time taken for this attempt (in microseconds).
    pub duration_us: Option<u64>,
    /// Number of faces in the result (if successful).
    pub result_faces: Option<usize>,
    /// Whether scoped make-connected was used.
    pub scoped_make_connected: bool,
    /// Whether fallback to global cleanup occurred.
    pub global_fallback: bool,
}

impl BooleanAttemptDiagnostic {
    /// Creates a new diagnostic for an attempt.
    pub fn new(attempt_number: usize, fuzzy_tolerance: f64) -> Self {
        Self {
            attempt_number,
            fuzzy_tolerance,
            ..Self::default()
        }
    }
}

/// Comprehensive diagnostic report for a boolean operation with retries.
#[derive(Debug, Clone, Default)]
pub struct BooleanDiagnosticReport {
    /// All attempt diagnostics.
    pub attempts: Vec<BooleanAttemptDiagnostic>,
    /// Total number of attempts.
    pub total_attempts: usize,
    /// Number of successful attempts (should be 0 or 1).
    pub successful_attempts: usize,
    /// The failure class that was most common (if any failures).
    pub dominant_failure_class: Option<BooleanFailureClass>,
    /// The final successful configuration (if successful).
    pub final_config: Option<FinalSuccessfulConfig>,
    /// Total time taken across all attempts (in microseconds).
    pub total_duration_us: u64,
    /// The retry policy used.
    pub retry_policy: Option<RetryPolicy>,
    /// Whether glue mode was ultimately needed.
    pub glue_mode_needed: bool,
    /// Whether MakeConnected cleanup was ultimately needed.
    pub make_connected_needed: bool,
    /// Summary of recovery strategies applied.
    pub recovery_strategies_applied: Vec<RecoveryStrategy>,
}

impl BooleanDiagnosticReport {
    /// Creates a new empty diagnostic report.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds an attempt diagnostic to the report.
    pub fn add_attempt(&mut self, attempt: BooleanAttemptDiagnostic) {
        if attempt.success {
            self.successful_attempts += 1;
        }
        self.total_attempts += 1;
        self.attempts.push(attempt);
    }

    /// Computes the dominant failure class from failed attempts.
    pub fn compute_dominant_failure_class(&mut self) {
        use std::collections::HashMap;
        let mut counts: HashMap<BooleanFailureClass, usize> = HashMap::new();

        for attempt in &self.attempts {
            if let Some(class) = attempt.failure_class {
                *counts.entry(class).or_insert(0) += 1;
            }
        }

        self.dominant_failure_class = counts
            .into_iter()
            .max_by_key(|(_, count)| *count)
            .map(|(class, _)| class);
    }

    /// Finalizes the report after all attempts are complete.
    pub fn finalize(&mut self) {
        self.compute_dominant_failure_class();

        // Determine if glue mode or MakeConnected were needed
        for attempt in &self.attempts {
            if attempt.success {
                self.glue_mode_needed = attempt.glue_enabled;
                self.make_connected_needed = attempt.make_connected_run;

                // Record the final successful configuration
                self.final_config = Some(FinalSuccessfulConfig {
                    fuzzy_tolerance: attempt.fuzzy_tolerance,
                    glue_enabled: attempt.glue_enabled,
                    glue_tolerance: attempt.glue_tolerance,
                    make_connected_run: attempt.make_connected_run,
                    make_connected_passes: attempt.make_connected_passes,
                    scoped_make_connected: attempt.scoped_make_connected,
                });
            }

            if let Some(strategy) = attempt.recovery_strategy
                && !self.recovery_strategies_applied.contains(&strategy) {
                    self.recovery_strategies_applied.push(strategy);
                }
        }
    }

    /// Returns whether the operation ultimately succeeded.
    pub fn is_success(&self) -> bool {
        self.successful_attempts > 0
    }

    /// Returns a summary string for logging.
    pub fn summary(&self) -> String {
        if self.is_success() {
            format!(
                "Boolean operation succeeded after {} attempts (dominant failure: {})",
                self.total_attempts,
                self.dominant_failure_class
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| "none".to_string())
            )
        } else {
            format!(
                "Boolean operation failed after {} attempts (dominant failure: {})",
                self.total_attempts,
                self.dominant_failure_class
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| "unknown".to_string())
            )
        }
    }
}

/// Configuration that produced a successful result.
#[derive(Debug, Clone)]
pub struct FinalSuccessfulConfig {
    /// Fuzzy tolerance that succeeded.
    pub fuzzy_tolerance: f64,
    /// Whether glue mode was enabled.
    pub glue_enabled: bool,
    /// Glue tolerance used.
    pub glue_tolerance: f64,
    /// Whether MakeConnected was run.
    pub make_connected_run: bool,
    /// Number of MakeConnected passes.
    pub make_connected_passes: usize,
    /// Whether scoped make-connected was used.
    pub scoped_make_connected: bool,
}

/// Failure analyzer for detailed classification of boolean operation failures.
#[derive(Debug, Clone, Default)]
pub struct FailureAnalyzer {
    /// Threshold for detecting degenerate edges (squared length).
    pub degenerate_edge_threshold: f64,
    /// Threshold for detecting degenerate triangles (minimum area).
    pub degenerate_triangle_threshold: f64,
    /// Maximum allowed self-intersection distance.
    pub self_intersection_threshold: f64,
}

impl FailureAnalyzer {
    /// Creates a new failure analyzer with default thresholds.
    pub fn new() -> Self {
        Self {
            degenerate_edge_threshold: TOLERANCE_LEN_MIN,
            degenerate_triangle_threshold: TOLERANCE_CLAMP_MIN,
            self_intersection_threshold: TOLERANCE_MESH_LEGACY,
        }
    }

    /// Classifies a failure based on error message and context.
    pub fn classify_from_error(error_message: &str) -> BooleanFailureClass {
        let msg = error_message.to_lowercase();

        // Check self-intersection first since it contains "intersection"
        if msg.contains("self-intersect") {
            BooleanFailureClass::SelfIntersection
        } else if msg.contains("empty") || msg.contains("missing geometry") {
            BooleanFailureClass::InvalidInput
        } else if msg.contains("degenerate") {
            BooleanFailureClass::DegenerateTopology
        } else if msg.contains("nan") || msg.contains("infinity") || msg.contains("numerical") {
            BooleanFailureClass::NumericalInstability
        } else if msg.contains("invalid") || msg.contains("non-manifold") || msg.contains("open shell") {
            BooleanFailureClass::InvalidResult
        } else if msg.contains("intersection") || msg.contains("missing curve") {
            BooleanFailureClass::IncompleteIntersection
        } else {
            BooleanFailureClass::Unknown
        }
    }

    /// Determines the appropriate recovery strategy for a failure class.
    pub fn determine_recovery_strategy(
        failure_class: BooleanFailureClass,
        attempt_number: usize,
        policy: &RetryPolicy,
    ) -> RecoveryStrategy {
        if !failure_class.is_recoverable() {
            return RecoveryStrategy::None;
        }

        match failure_class {
            BooleanFailureClass::DegenerateTopology => {
                RecoveryStrategy::Combined {
                    use_glue: policy.should_enable_glue(attempt_number),
                    run_make_connected: true,
                    increase_fuzzy: attempt_number > 1,
                }
            }
            BooleanFailureClass::NumericalInstability => {
                RecoveryStrategy::Combined {
                    use_glue: policy.should_enable_glue(attempt_number),
                    run_make_connected: attempt_number > 2,
                    increase_fuzzy: true,
                }
            }
            BooleanFailureClass::InvalidResult => {
                if policy.try_algorithm_variants && attempt_number > 2 {
                    RecoveryStrategy::AlgorithmVariant
                } else {
                    RecoveryStrategy::MakeConnectedCleanup
                }
            }
            BooleanFailureClass::IncompleteIntersection => {
                RecoveryStrategy::EnableGlueMode
            }
            BooleanFailureClass::SelfIntersection => {
                RecoveryStrategy::MakeConnectedCleanup
            }
            _ => RecoveryStrategy::None,
        }
    }
}

/// Builder for creating customized retry policies.
#[derive(Debug, Clone, Default)]
pub struct RetryPolicyBuilder {
    policy: RetryPolicy,
}

impl RetryPolicyBuilder {
    /// Creates a new builder with default settings.
    pub fn new() -> Self {
        Self {
            policy: RetryPolicy::default(),
        }
    }

    /// Starts from a conservative policy.
    pub fn conservative() -> Self {
        Self {
            policy: RetryPolicy::conservative(),
        }
    }

    /// Starts from an aggressive policy.
    pub fn aggressive() -> Self {
        Self {
            policy: RetryPolicy::aggressive(),
        }
    }

    /// Sets the maximum number of attempts.
    pub fn max_attempts(mut self, max: usize) -> Self {
        self.policy.max_attempts = max;
        self
    }

    /// Sets the fuzzy tolerance growth factor.
    pub fn fuzzy_growth_factor(mut self, factor: f64) -> Self {
        self.policy.fuzzy_growth_factor = factor;
        self
    }

    /// Sets the fuzzy tolerance cap.
    pub fn fuzzy_tolerance_cap(mut self, cap: f64) -> Self {
        self.policy.fuzzy_tolerance_cap = cap;
        self
    }

    /// Sets when to enable glue mode.
    pub fn enable_glue_after(mut self, n_failures: usize) -> Self {
        self.policy.enable_glue_after_n_failures = n_failures;
        self
    }

    /// Sets the glue tolerance.
    pub fn glue_tolerance(mut self, tol: f64) -> Self {
        self.policy.glue_tolerance = tol;
        self
    }

    /// Sets the MakeConnected aggressiveness (1-10).
    pub fn make_connected_aggressiveness(mut self, level: u32) -> Self {
        self.policy.make_connected_aggressiveness = level.clamp(1, 10);
        self
    }

    /// Sets the maximum MakeConnected passes.
    pub fn make_connected_max_passes(mut self, passes: usize) -> Self {
        self.policy.make_connected_max_passes = passes;
        self
    }

    /// Sets the initial MakeConnected tolerance.
    pub fn make_connected_initial_tolerance(mut self, tol: f64) -> Self {
        self.policy.make_connected_initial_tolerance = tol;
        self
    }

    /// Sets whether to use scoped MakeConnected.
    pub fn use_scoped_make_connected(mut self, use_scoped: bool) -> Self {
        self.policy.use_scoped_make_connected = use_scoped;
        self
    }

    /// Sets whether to fall back to global cleanup.
    pub fn fallback_to_global_cleanup(mut self, fallback: bool) -> Self {
        self.policy.fallback_to_global_cleanup = fallback;
        self
    }

    /// Sets whether to try algorithm variants.
    pub fn try_algorithm_variants(mut self, try_variants: bool) -> Self {
        self.policy.try_algorithm_variants = try_variants;
        self
    }

    /// Enables verbose diagnostics.
    pub fn verbose_diagnostics(mut self, verbose: bool) -> Self {
        self.policy.verbose_diagnostics = verbose;
        self
    }

    /// Builds the final retry policy.
    pub fn build(self) -> RetryPolicy {
        self.policy
    }
}

/// Map a `BooleanRetryClass` (used by the robust pipeline) to a `BooleanFailureClass`
/// (used by this module) for diagnostic reporting.
fn boolean_retry_class_to_failure_class(rc: BooleanRetryClass) -> BooleanFailureClass {
    match rc {
        BooleanRetryClass::FatalInput | BooleanRetryClass::IncompleteData => {
            BooleanFailureClass::InvalidInput
        }
        BooleanRetryClass::DegenerateTopology => BooleanFailureClass::DegenerateTopology,
        BooleanRetryClass::NumericalInstability => BooleanFailureClass::NumericalInstability,
    }
}

/// Convert a [`RetryPolicy`] into [`BooleanRobustOptions`] for use with
/// [`crate::boolean_op_robust`].
///
/// This bridges the high-level retry policy in this module with the robust boolean
/// execution pipeline. The generated options use the policy's fuzzy growth factor,
/// glue settings, and MakeConnected configuration.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::boolean::{RetryPolicy, retry_policy_to_robust_options};
/// use rcad_algorithms::BooleanOptions;
///
/// let policy = RetryPolicy::default();
/// let base = BooleanOptions::default();
/// let robust_opts = retry_policy_to_robust_options(&policy, base);
/// assert!(robust_opts.fuzzy_retry_ladder.len() <= policy.max_attempts);
/// ```
pub fn retry_policy_to_robust_options(
    policy: &RetryPolicy,
    base_options: BooleanOptions,
) -> BooleanRobustOptions {
    let start_tol = if base_options.fuzzy_tol > 0.0 {
        base_options.fuzzy_tol
    } else {
        TOLERANCE_ABS
    };
    let mut ladder = Vec::new();
    let mut tol = start_tol;
    for _ in 0..policy.max_attempts.saturating_sub(1) {
        tol = policy.next_fuzzy_tolerance(tol);
        tol = tol.min(policy.fuzzy_tolerance_cap);
        if !ladder.iter().any(|&v: &f64| (v - tol).abs() <= TOLERANCE_FLOAT_DEDUP) {
            ladder.push(tol);
        }
        if tol >= policy.fuzzy_tolerance_cap {
            break;
        }
    }

    let mut opts = base_options;
    opts.run_make_connected = true;
    opts.make_connected_max_passes = policy.make_connected_max_passes;
    opts.make_connected_tolerance = policy.make_connected_initial_tolerance;
    opts.make_connected_tolerance_growth = policy.make_connected_tolerance_growth;
    opts.make_connected_scoped = policy.use_scoped_make_connected;
    opts.use_glue = true;
    opts.glue_tolerance = policy.glue_tolerance;

    BooleanRobustOptions {
        base: opts,
        fuzzy_retry_ladder: ladder,
        retry_policy: if policy.try_algorithm_variants {
            BooleanRetryPolicy::AdaptiveByFailureClass
        } else {
            BooleanRetryPolicy::Conservative
        },
        extreme_geometry: ExtremeGeometryRetryConfig::default(),
    }
}

/// Perform a boolean operation with a retry policy and receive a diagnostic report.
///
/// ## OCCT alignment
///
/// The core boolean logic (DS + PaveFiller + BooleanBuilder) invoked by this
/// function's inner pipeline (`boolean_op_robust` 鈫?`boolean_op_with_options`)
/// follows the OCCT BRepAlgoAPI pipeline:
///
/// - **BRepAlgoAPI_Common::Build()** (BRepAlgoAPI.cxx) for `BooleanOpType::Intersection`
/// - **BRepAlgoAPI_Fuse::Build()** (BRepAlgoAPI.cxx) for `BooleanOpType::Union`
/// - **BRepAlgoAPI_Cut::Build()** (BRepAlgoAPI.cxx) for `BooleanOpType::Difference`
///
/// Each OCCT operation calls `BOPAlgo_BOP::Perform()` which runs:
/// 1. 鉁?`BOPAlgo_PaveFiller::Perform()` 鈥?EE/EF/FF intersection detection
/// 2. 鉁?`BOPAlgo_Builder::Build()` 鈥?face splitting + result assembly
///
/// ## rcad-specific enhancement
///
/// **`RetryPolicy`** and the retry ladder (`BooleanRobustOptions::fuzzy_retry_ladder`)
/// are **rcad-specific** 鈥?OCCT's `BRepAlgoAPI` does not implement retry logic.
/// OCCT relies on a single pass with deterministic tolerances (`Precision::Confusion()`).
/// The retry mechanism here handles edge cases where the initial tolerance fails
/// due to near-degenerate geometry, by gradually increasing the fuzzy tolerance.
///
/// This wraps [`crate::boolean_op_robust`] with the high-level [`RetryPolicy`],
/// returning a [`BooleanDiagnosticReport`] that provides detailed per-attempt
/// diagnostics. Use this when you need observable retry data or want to customise
/// the retry behaviour beyond the built-in presets.
///
/// # Example
///
/// ```rust
/// use rcad_algorithms::boolean::{RetryPolicy, boolean_op_with_retry_policy};
/// use rcad_algorithms::{BooleanOpType, BooleanOptions};
/// use rcad_kernel::{BRep, PrimitiveSolid};
///
/// let a = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 1.0 });
/// let b = BRep::from_primitive(PrimitiveSolid::Sphere { radius: 0.8 });
/// let policy = RetryPolicy::conservative();
/// let result = boolean_op_with_retry_policy(
///     BooleanOpType::Intersection, &a, &b, &policy, BooleanOptions::default()
/// );
/// match result {
///     Ok((brep, report)) => {
///         println!("Succeeded after {} attempts", report.total_attempts);
///     }
///     Err(err) => {
///         println!("Operation failed: {}", err);
///     }
/// }
/// ```
pub fn boolean_op_with_retry_policy(
    op: BooleanOpType,
    a: &BRep,
    b: &BRep,
    policy: &RetryPolicy,
    base_options: BooleanOptions,
) -> Result<(BRep, BooleanDiagnosticReport), BooleanError> {
    let start = std::time::Instant::now();
    let robust_options = retry_policy_to_robust_options(policy, base_options);
    let mut diagnostic = BooleanDiagnosticReport::new();
    diagnostic.retry_policy = Some(policy.clone());

    match boolean_op_robust(op, a, b, robust_options) {
        Ok((brep, report)) => {
            for (i, attempt) in report.robust_attempts.iter().enumerate() {
                let mut ad = BooleanAttemptDiagnostic::new(i + 1, attempt.fuzzy_tol);
                ad.glue_enabled = attempt.used_glue;
                ad.glue_tolerance = attempt.glue_tolerance;
                ad.success = attempt.success;

                if let Some(rc) = attempt.retry_class {
                    // Use FailureAnalyzer for finer-grained classification when possible
                    let failure_class = attempt
                        .error_message
                        .as_deref()
                        .map(FailureAnalyzer::classify_from_error)
                        .unwrap_or_else(|| boolean_retry_class_to_failure_class(rc));
                    ad.failure_class = Some(failure_class);
                    ad.error_message = attempt.error_message.clone();
                }

                diagnostic.add_attempt(ad);
            }

            diagnostic.total_duration_us = start.elapsed().as_micros() as u64;
            diagnostic.finalize();
            Ok((brep, diagnostic))
        }
        Err(err) => {
            // Record a single diagnostic entry for the final failure.
            let mut ad = BooleanAttemptDiagnostic::new(1, 0.0);
            ad.success = false;
            ad.failure_class = Some(classify_boolean_failure(&err));
            ad.error_message = Some(format!("{err:?}"));
            diagnostic.add_attempt(ad);

            diagnostic.total_duration_us = start.elapsed().as_micros() as u64;
            diagnostic.finalize();
            Err(err)
        }
    }
}

