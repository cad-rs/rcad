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
    pub fn is_recoverable(&self) -> bool {
        false
    }
    pub fn suggested_recovery(&self) -> RecoveryStrategy {
        RecoveryStrategy::None
    }
    pub fn description(&self) -> &'static str {
        "retry system removed"
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecoveryStrategy {
    None,
    MakeConnectedCleanup,
    IncreaseFuzzyTolerance,
    AlgorithmVariant,
    EnableGlueMode,
    Combined {
        use_glue: bool,
        run_make_connected: bool,
        increase_fuzzy: bool,
    },
}

impl RecoveryStrategy {
    pub fn description(&self) -> &'static str {
        "retry system removed"
    }
}

#[derive(Debug, Clone)]
pub struct RetryPolicy;
impl RetryPolicy {
    pub fn new() -> Self {
        Self
    }
    pub fn default() -> Self {
        Self
    }
    pub fn conservative() -> Self {
        Self
    }
    pub fn aggressive() -> Self {
        Self
    }
    pub fn for_numerical_instability() -> Self {
        Self
    }
    pub fn for_degenerate_topology() -> Self {
        Self
    }
    pub fn for_incomplete_intersection() -> Self {
        Self
    }
    pub fn next_fuzzy_tolerance(&self, _c: f64) -> f64 {
        0.0
    }
    pub fn should_enable_glue(&self, _n: usize) -> bool {
        false
    }
    pub fn make_connected_tolerance(&self, _p: usize) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone, Default)]
pub struct BooleanAttemptDiagnostic {
    pub attempt_number: usize,
}
#[derive(Debug, Clone, Default)]
pub struct BooleanDiagnosticReport {
    pub total_attempts: usize,
}
#[derive(Debug, Clone)]
pub struct FinalSuccessfulConfig {
    pub fuzzy_tolerance: f64,
    pub glue_enabled: bool,
    pub glue_tolerance: f64,
    pub make_connected_run: bool,
    pub make_connected_passes: usize,
    pub scoped_make_connected: bool,
}
#[derive(Debug, Clone, Default)]
pub struct FailureAnalyzer;
#[derive(Debug, Clone, Default)]
pub struct RetryPolicyBuilder;
impl RetryPolicyBuilder {
    pub fn new() -> Self {
        Self
    }
    pub fn conservative() -> Self {
        Self
    }
    pub fn aggressive() -> Self {
        Self
    }
    pub fn build(self) -> RetryPolicy {
        RetryPolicy
    }
}
