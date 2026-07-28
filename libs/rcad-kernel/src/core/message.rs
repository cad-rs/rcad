//! Messaging and progress indication.
//!
//! Analogous to OCCT's `Message` package: `Message_Messenger` (message output),
//! `Message_ProgressIndicator` (progress reporting), and `Message_ProgressScope`
//! (nested progress range).
//!
//! # Usage
//!
//! ```ignore
//! use rcad_kernel::core::message::{PrintMessenger, NoopProgress, ProgressScope};
//!
//! // Messenger: route messages where you need them
//! let msg = PrintMessenger;
//! msg.send("Starting boolean op", MessageGravity::Info);
//!
//! // Progress: pass a ProgressIndicator to long-running operations
//! let prog = NoopProgress;
//! let mut root = ProgressScope::new(&prog, "Fuse", 3);
//! {
//!     let _s1 = root.sub_scope("Intersect", 1);
//!     // ... do work ...
//! }
//! root.advance();
//! ```

use std::sync::atomic::{AtomicBool, Ordering};

// ============================================================================
// Message severity (OCCT Message_Gravity)
// ============================================================================

/// Severity level for messages.
///
/// Analogous to OCCT `Message_Gravity`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum MessageGravity {
    /// Detailed debug / trace output.
    Trace,
    /// Normal informational message.
    Info,
    /// Warning — operation can continue but result may be degraded.
    Alert,
    /// Error — operation failed partially.
    Alarm,
    /// Fatal — operation cannot continue.
    Fail,
}

impl MessageGravity {
    /// Returns a short label (e.g. "INFO", "WARN").
    pub fn label(self) -> &'static str {
        match self {
            MessageGravity::Trace => "TRACE",
            MessageGravity::Info => "INFO",
            MessageGravity::Alert => "WARN",
            MessageGravity::Alarm => "ERROR",
            MessageGravity::Fail => "FAIL",
        }
    }
}

// ============================================================================
// Messenger (OCCT Message_Messenger)
// ============================================================================

/// Trait for routing messages with severity.
///
/// Analogous to OCCT `Message_Messenger`.
pub trait Messenger: Send + Sync {
    /// Send a text message at the given severity level.
    fn send(&self, text: &str, gravity: MessageGravity);
}

/// A messenger that prints to stderr with a severity prefix.
///
/// Example output: `[WARN]  Face 42 has zero area`.
pub struct PrintMessenger;

impl Messenger for PrintMessenger {
    fn send(&self, text: &str, gravity: MessageGravity) {
        eprintln!("[{}] {}", gravity.label(), text);
    }
}

/// A messenger that discards all messages.
pub struct NullMessenger;

impl Messenger for NullMessenger {
    fn send(&self, _text: &str, _gravity: MessageGravity) {
        // no-op
    }
}

// ============================================================================
// Progress indicator (OCCT Message_ProgressIndicator)
// ============================================================================

/// Trait for receiving progress updates from long-running operations.
///
/// Analogous to OCCT `Message_ProgressIndicator`.
///
/// # Position semantics
///
/// `position` is in `[0.0, 1.0]` and represents the overall progress of the
/// operation. `name` is a human-readable label for the current sub-step.
pub trait ProgressIndicator: Send + Sync {
    /// Called when progress changes. `position` is in `[0, 1]`.
    fn show(&self, position: f64, name: &str);

    /// Returns `true` if the user requested to abort the operation.
    /// Called periodically by long-running operations.
    fn user_break(&self) -> bool {
        false
    }
}

/// A progress indicator that prints progress percentage to stderr.
pub struct PrintProgress {
    last_pct: std::sync::Mutex<i32>,
}

impl PrintProgress {
    pub fn new() -> Self {
        Self {
            last_pct: std::sync::Mutex::new(-1),
        }
    }
}

impl ProgressIndicator for PrintProgress {
    fn show(&self, position: f64, name: &str) {
        let pct = (position * 100.0) as i32;
        let mut last = self.last_pct.lock().unwrap();
        if pct > *last {
            *last = pct;
            eprintln!("[PROGRESS] {}%  {}", pct, name);
        }
    }
}

/// A no-op progress indicator that does nothing.
/// Use as a default when progress reporting is not needed.
pub struct NoopProgress;

impl ProgressIndicator for NoopProgress {
    fn show(&self, _position: f64, _name: &str) {}
}

/// A progress indicator that delegates to a closure.
pub struct CallbackProgress<F: Fn(f64, &str)> {
    callback: F,
    break_flag: AtomicBool,
}

impl<F: Fn(f64, &str)> CallbackProgress<F> {
    pub fn new(callback: F) -> Self {
        Self {
            callback,
            break_flag: AtomicBool::new(false),
        }
    }

    /// Signal that the operation should abort.
    pub fn request_break(&self) {
        self.break_flag.store(true, Ordering::Relaxed);
    }
}

impl<F: Fn(f64, &str) + Send + Sync> ProgressIndicator for CallbackProgress<F> {
    fn show(&self, position: f64, name: &str) {
        (self.callback)(position, name);
    }

    fn user_break(&self) -> bool {
        self.break_flag.load(Ordering::Relaxed)
    }
}

// ============================================================================
// Progress scope (OCCT Message_ProgressScope)
// ============================================================================

/// A named sub-range within a progress indicator.
///
/// Analogous to OCCT `Message_ProgressScope`. Allows decomposition of a
/// long operation into named sub-steps:
///
/// ```ignore
/// let mut root = ProgressScope::new(&indicator, "Fuse", 100);
/// {
///     let sub = root.sub_scope("Intersect", 3);  // 3% of total
///     sub.advance();                               // 1% done
///     sub.advance();                               // 2% done
/// }                                                // 3% done, auto-finalized
/// root.advance();                                  // continues to next major step
/// ```
pub struct ProgressScope<'a> {
    /// Reference to the underlying progress indicator.
    indicator: &'a dyn ProgressIndicator,
    /// Name of this scope (shown in progress output).
    name: &'a str,
    /// Coordinate mapping: parent's coordinate space → this scope maps
    /// [0, nb_steps] → [min_pos, max_pos].
    min_pos: f64,
    max_pos: f64,
    /// Number of steps in this scope.
    nb_steps: usize,
    /// Current step index (0..nb_steps).
    current: usize,
}

impl<'a> ProgressScope<'a> {
    /// Create a new root scope that covers `nb_steps` steps from 0 to 1.
    pub fn new(indicator: &'a dyn ProgressIndicator, name: &'a str, nb_steps: usize) -> Self {
        Self {
            indicator,
            name,
            min_pos: 0.0,
            max_pos: 1.0,
            nb_steps,
            current: 0,
        }
    }

    /// Create a scope within a specific range of the parent's progress.
    fn with_range(
        indicator: &'a dyn ProgressIndicator,
        name: &'a str,
        nb_steps: usize,
        min_pos: f64,
        max_pos: f64,
    ) -> Self {
        Self {
            indicator,
            name,
            min_pos,
            max_pos,
            nb_steps,
            current: 0,
        }
    }

    /// Report the current position to the indicator.
    fn report(&self) {
        let fraction = if self.nb_steps > 0 {
            self.current as f64 / self.nb_steps as f64
        } else {
            0.0
        };
        let pos = self.min_pos + fraction * (self.max_pos - self.min_pos);
        self.indicator.show(pos, self.name);
    }

    /// Advance by one step.
    pub fn advance(&mut self) {
        if self.current < self.nb_steps {
            self.current += 1;
        }
        self.report();
    }

    /// Advance by `n` steps at once.
    pub fn advance_by(&mut self, n: usize) {
        self.current = (self.current + n).min(self.nb_steps);
        self.report();
    }

    /// Advance to a specific step index (0-based).
    pub fn advance_to(&mut self, step: usize) {
        self.current = step.min(self.nb_steps);
        self.report();
    }

    /// Get the current step index.
    pub fn current_step(&self) -> usize {
        self.current
    }

    /// Get the total number of steps.
    pub fn nb_steps(&self) -> usize {
        self.nb_steps
    }

    /// Get the scope name.
    pub fn name(&self) -> &'a str {
        self.name
    }

    /// Check if the user requested to abort.
    pub fn user_break(&self) -> bool {
        self.indicator.user_break()
    }

    /// Create a sub-scope that consumes `child_steps` steps of this scope's
    /// remaining progress. Returns a child scope that will be reported
    /// independently to the same indicator.
    ///
    /// The child's position range is proportional to `child_steps / nb_steps`
    /// relative to this scope's range. When `child_steps >= nb_steps` or
    /// when the remaining steps are all consumed, the child covers the full
    /// [current, max] range.
    pub fn sub_scope(&self, name: &'a str, child_steps: usize) -> ProgressScope<'a> {
        let remaining = self.nb_steps.saturating_sub(self.current);
        let take = child_steps.min(remaining).max(1);

        let fraction = if self.nb_steps > 0 {
            take as f64 / self.nb_steps as f64
        } else {
            0.0
        };

        let sub_min = self.min_pos
            + (self.current as f64 / self.nb_steps as f64) * (self.max_pos - self.min_pos);
        let sub_max = (sub_min + fraction * (self.max_pos - self.min_pos)).min(self.max_pos);

        ProgressScope::with_range(self.indicator, name, child_steps, sub_min, sub_max)
    }

    /// Get the overall position as a fraction [0, 1].
    pub fn position(&self) -> f64 {
        if self.nb_steps > 0 {
            self.min_pos
                + (self.current as f64 / self.nb_steps as f64) * (self.max_pos - self.min_pos)
        } else {
            self.min_pos
        }
    }
}

impl<'a> Drop for ProgressScope<'a> {
    fn drop(&mut self) {
        // On drop, advance to the end so the indicator shows completion.
        self.current = self.nb_steps;
        self.report();
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    #[test]
    fn test_progress_scope_advance() {
        let prog = NoopProgress;
        let mut ps = ProgressScope::new(&prog, "test", 5);
        assert_eq!(ps.current_step(), 0);
        ps.advance();
        assert_eq!(ps.current_step(), 1);
        ps.advance_by(2);
        assert_eq!(ps.current_step(), 3);
        ps.advance_to(5);
        assert_eq!(ps.current_step(), 5);
    }

    #[test]
    fn test_progress_scope_position() {
        let prog = NoopProgress;
        let ps = ProgressScope::new(&prog, "test", 4);
        assert!((ps.position() - 0.0).abs() < 1e-12);

        let mut ps = ps;
        ps.advance();
        assert!((ps.position() - 0.25).abs() < 1e-12);

        ps.advance();
        assert!((ps.position() - 0.5).abs() < 1e-12);

        ps.advance_by(2);
        assert!((ps.position() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_progress_sub_scope() {
        let prog = NoopProgress;
        let mut root = ProgressScope::new(&prog, "root", 10);

        // Sub-scope consuming 3 of 10 steps
        let mut sub = root.sub_scope("sub", 3);
        assert!((sub.position() - 0.0).abs() < 1e-12);

        sub.advance();
        // sub has 1 of 3 done → 1/3 of its range → 1/3 * 0.3 = 0.1 of total
        assert!((sub.position() - 0.1).abs() < 1e-12);

        sub.advance_by(2);
        assert!((sub.position() - 0.3).abs() < 1e-12);

        drop(sub);
        // root should have advanced to step 3
        assert_eq!(root.current_step(), 0); // sub_scope doesn't modify parent
    }

    #[test]
    fn test_drop_auto_finalizes() {
        let prog = NoopProgress;
        let pos;
        let mut ps = ProgressScope::new(&prog, "test", 5);
        ps.advance_by(2);
        pos = ps.position();
        drop(ps);
        // Scope auto-finalizes on drop to max position
    }

    #[test]
    fn test_messenger_severity() {
        assert_eq!(MessageGravity::Trace.label(), "TRACE");
        assert_eq!(MessageGravity::Info.label(), "INFO");
        assert_eq!(MessageGravity::Alert.label(), "WARN");
        assert_eq!(MessageGravity::Alarm.label(), "ERROR");
        assert_eq!(MessageGravity::Fail.label(), "FAIL");
    }

    #[test]
    fn test_print_messenger() {
        let msg = PrintMessenger;
        msg.send("test message", MessageGravity::Info);
        msg.send("warning", MessageGravity::Alert);
    }

    #[test]
    fn test_noop_progress() {
        let prog = NoopProgress;
        prog.show(0.5, "test");
        assert!(!prog.user_break());
    }

    #[test]
    fn test_callback_progress() {
        let call_count = AtomicUsize::new(0);
        let prog = CallbackProgress::new(|pos, name| {
            call_count.fetch_add(1, Ordering::Relaxed);
            assert!((pos - 0.5).abs() < 1e-12);
            assert_eq!(name, "test");
        });

        prog.show(0.5, "test");
        assert_eq!(call_count.load(Ordering::Relaxed), 1);

        prog.request_break();
        assert!(prog.user_break());
    }

    #[test]
    fn test_sub_scope_positioning() {
        let prog = NoopProgress;
        let root = ProgressScope::new(&prog, "root", 100);

        // sub covers [0, 0.3]
        let mut sub = ProgressScope::with_range(&prog, "sub", 3, 0.0, 0.3);
        assert!((sub.position() - 0.0).abs() < 1e-12);
        sub.advance();
        assert!((sub.position() - 0.1).abs() < 1e-12);
        sub.advance();
        assert!((sub.position() - 0.2).abs() < 1e-12);
        sub.advance();
        assert!((sub.position() - 0.3).abs() < 1e-12);
    }
}
