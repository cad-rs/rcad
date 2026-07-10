//! Minimal stubs for types still referenced by legacy retry/tuning code.
//! The retry system was deleted -- these remain only to satisfy existing
//! references until the callers are cleaned up.

use crate::BooleanError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BooleanFailureClass {
    DegenerateTopology,
    NumericalInstability,
    InvalidResult,
    IncompleteIntersection,
    SelfIntersection,
    InvalidInput,
    #[default]
    Unknown,
}

impl BooleanFailureClass {
    pub fn is_recoverable(&self) -> bool { false }
    pub fn suggested_recovery(&self) -> RecoveryStrategy { RecoveryStrategy::None }
    pub fn description(&self) -> &'static str { "retry system removed" }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecoveryStrategy {
    None,
    MakeConnectedCleanup,
    IncreaseFuzzyTolerance,
    AlgorithmVariant,
    EnableGlueMode,
    Combined { use_glue: bool, run_make_connected: bool, increase_fuzzy: bool },
}

impl RecoveryStrategy {
    pub fn description(&self) -> &'static str { "retry system removed" }
}

#[derive(Debug, Clone)]
pub struct RetryPolicy;
impl RetryPolicy {
    pub fn new() -> Self { Self }
    pub fn default() -> Self { Self }
    pub fn conservative() -> Self { Self }
    pub fn aggressive() -> Self { Self }
    pub fn for_numerical_instability() -> Self { Self }
    pub fn for_degenerate_topology() -> Self { Self }
    pub fn for_incomplete_intersection() -> Self { Self }
    pub fn next_fuzzy_tolerance(&self, _c: f64) -> f64 { 0.0 }
    pub fn should_enable_glue(&self, _n: usize) -> bool { false }
    pub fn make_connected_tolerance(&self, _p: usize) -> f64 { 0.0 }
}

#[derive(Debug, Clone, Default)]
pub struct BooleanAttemptDiagnostic { pub attempt_number: usize }
#[derive(Debug, Clone, Default)]
pub struct BooleanDiagnosticReport { pub total_attempts: usize }
#[derive(Debug, Clone)]
pub struct FinalSuccessfulConfig { pub fuzzy_tolerance: f64, pub glue_enabled: bool, pub glue_tolerance: f64, pub make_connected_run: bool, pub make_connected_passes: usize, pub scoped_make_connected: bool, }
#[derive(Debug, Clone, Default)]
pub struct FailureAnalyzer;
#[derive(Debug, Clone, Default)]
pub struct RetryPolicyBuilder;
impl RetryPolicyBuilder {
    pub fn new() -> Self { Self }
    pub fn conservative() -> Self { Self }
    pub fn aggressive() -> Self { Self }
    pub fn build(self) -> RetryPolicy { RetryPolicy }
}

// Legacy types from removed lib_inline -- kept as stubs for existing references.
// TODO: remove when all callers are migrated.

#[derive(Debug, Clone)]
pub struct HealingOptions;
impl Default for HealingOptions {
    fn default() -> Self { Self }
}

#[derive(Debug, Clone)]
pub struct BooleanOptions {
    pub use_bvh: bool,
    pub run_healing: bool,
    pub healing: HealingOptions,
    pub fuzzy_tol: f64,
    pub glue_tol: f64,
    pub glue_tolerance: f64,
    pub glue_mode: bool,
    pub use_glue: bool,
    pub run_make_connected: bool,
    pub make_connected_passes: usize,
    pub include_history: bool,
    pub run_fix_self_intersection: bool,
}
impl Default for BooleanOptions {
    fn default() -> Self {
        Self {
            use_bvh: true, run_healing: false, healing: HealingOptions,
            fuzzy_tol: 0.0, glue_tol: 0.0, glue_tolerance: 0.0,
            glue_mode: false, use_glue: false,
            run_make_connected: false, make_connected_passes: 0,
            include_history: false, run_fix_self_intersection: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanRetryPolicy {
    Conservative,
    AdaptiveByFailureClass,
    Aggressive,
}

#[derive(Debug, Clone)]
pub struct ExtremeGeometryRetryConfig {
    pub policy: ExtremeGeometryRetryPolicy,
    pub check_near_tangent: bool,
    pub check_high_aspect_ratio: bool,
    pub max_iterations: u32,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExtremeGeometryRetryPolicy { #[default] Default, Aggressive, Disabled }
impl Default for ExtremeGeometryRetryConfig {
    fn default() -> Self {
        Self { policy: ExtremeGeometryRetryPolicy::Default,
               check_near_tangent: true, check_high_aspect_ratio: true,
               max_iterations: 5 }
    }
}
impl ExtremeGeometryRetryConfig {
    pub fn geometry_aware() -> Self { Self::default() }
}

#[derive(Debug, Clone)]
pub struct BooleanRobustOptions {
    pub base: BooleanOptions,
    pub fuzzy_retry_ladder: Vec<f64>,
    pub retry_policy: BooleanRetryPolicy,
    pub extreme_geometry: ExtremeGeometryRetryConfig,
}
impl Default for BooleanRobustOptions {
    fn default() -> Self {
        Self {
            base: BooleanOptions::default(),
            fuzzy_retry_ladder: vec![],
            retry_policy: BooleanRetryPolicy::AdaptiveByFailureClass,
            extreme_geometry: ExtremeGeometryRetryConfig::default(),
        }
    }
}
impl BooleanRobustOptions {
    pub fn for_fea() -> Self { Self::default() }
    pub fn for_mechanical_multiscale() -> Self { Self::default() }
}

pub fn retry_policy_to_robust_options(
    _policy: &RetryPolicy, _base: BooleanOptions,
) -> BooleanRobustOptions {
    BooleanRobustOptions::default()
}

// Additional stubs for removed old API
use crate::BooleanOpType;
use rcad_kernel::topods;

pub struct GeneralFuseHistory { pub steps: Vec<crate::history::BooleanHistory> }
impl Default for GeneralFuseHistory { fn default() -> Self { Self { steps: vec![] } } }

pub fn general_fuse(_parts: &[topods::BRep]) -> Result<topods::BRep, crate::BooleanError> {
    Err(crate::BooleanError::DegenerateResult("stub"))
}
pub fn general_fuse_with_history(_parts: &[topods::BRep]) -> Result<(topods::BRep, GeneralFuseHistory), crate::BooleanError> {
    Err(crate::BooleanError::DegenerateResult("stub"))
}
pub fn boolean_op_pave_fill_build(_op: BooleanOpType, _a: &topods::BRep, _b: &topods::BRep) -> Result<topods::BRep, crate::BooleanError> {
    crate::bop_occt_union::boolean_op_generic(_op, _a, _b)
}

pub mod shape_analysis {
    pub fn detect_surface_self_intersection(_brep: &super::topods::BRep, _face_idx: usize) -> bool { false }
}
