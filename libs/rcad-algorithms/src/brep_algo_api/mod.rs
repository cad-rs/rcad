//! High-level boolean algorithms.
//!
//! Fuse (Union), Common (Intersection), Cut (Difference), Section
//! with `build`/`shape`/`is_done`/`history` pattern.
//! Uses `topods::BRep` for input/output.
//!
//! ## OCCT mapping
//!
//! | Rust (this module)          | OCCT class                        | Align |
//! |------------------------------|-----------------------------------|-------|
//! | `BooleanOptions`             | `BOPAlgo_Options`                 |   |
//! | `BooleanOp`                  | `BRepAlgoAPI_Fuse` / Common / Cut |   |
//! | `BooleanOp::build()`         | `BRepAlgoAPI_BuilderShape::Build` |   |
//! | `BooleanOp::shape()`         | `BRepAlgoAPI_BuilderShape::Shape` |   |
//! | `BooleanOp::is_done()`       | `BRepAlgoAPI_Algo::IsDone`        |   |
//! | `SectionOp`                  | `BRepAlgoAPI_Section`             |   |
//! | `fuse()` / `common()` / `cut()` | convenience wrapper (no OCCT equivalent) | |

pub mod argument_analyzer;
pub mod builder_operation;
pub mod section;

use rcad_kernel::topods;

use crate::builder::{BooleanError, BooleanOpType};
use crate::history::BooleanHistory;

/// Options for boolean operations.
///
/// `BOPAlgo_Options` — fuzzy tolerance, parallel mode, BVH, etc.
#[derive(Debug, Clone)]
pub struct BooleanOptions {
    /// `BOPAlgo_Options::SetFuzzyValue`
    pub fuzzy_value: f64,
    /// `BOPAlgo_Options::SetParallel`
    pub parallel: bool,
    /// rcad-specific: BVH acceleration toggle (no direct OCCT option)
    pub use_bvh: bool,
    /// `BOPAlgo_Options::SetGlue` (GlueOff/GlueShift/GlueFull)
    pub glue_enabled: bool,
    /// glue tolerance for Shift/Full modes
    pub glue_tolerance: f64,
    /// `BRepAlgoAPI_BuilderOperation` history tracking
    pub track_history: bool,
    /// rcad-specific: auto-heal result after boolean
    pub heal_result: bool,
}

impl Default for BooleanOptions {
    fn default() -> Self {
        Self {
            fuzzy_value: 0.0, parallel: false, use_bvh: false,
            glue_enabled: false, glue_tolerance: 1e-7,
            track_history: true, heal_result: false,
        }
    }
}

impl BooleanOptions {
    pub fn with_fuzzy_value(mut self, value: f64) -> Self { self.fuzzy_value = value; self }
    pub fn with_parallel(mut self, parallel: bool) -> Self { self.parallel = parallel; self }
    pub fn with_bvh(mut self, use_bvh: bool) -> Self { self.use_bvh = use_bvh; self }
    pub fn with_history(mut self, track: bool) -> Self { self.track_history = track; self }
}

/// A boolean operation between two shapes.
///
/// replaces `BRepAlgoAPI_Fuse`, `BRepAlgoAPI_Common`,
/// `BRepAlgoAPI_Cut` with a single struct + `BooleanOpType` discriminator.
/// API pattern (`new` -> `build` -> `shape` -> `is_done`) matches
/// `BRepAlgoAPI_BuilderShape`.
///
/// After `build()`, call `shape()` for the result or `error()` on failure.
pub struct BooleanOp<'a> {
    shape1: &'a topods::BRep,
    shape2: &'a topods::BRep,
    op_type: BooleanOpType,
    options: BooleanOptions,
    result: Option<topods::BRep>,
    history: Option<BooleanHistory>,
    err: Option<BooleanError>,
}

impl<'a> BooleanOp<'a> {
    /// constructor matching `BRepAlgoAPI_Fuse(a, b)` / `Common` / `Cut`,
    /// with `BooleanOpType` selecting the operation type.
    pub fn new(op: BooleanOpType, shape1: &'a topods::BRep, shape2: &'a topods::BRep) -> Self {
        Self { shape1, shape2, op_type: op, options: BooleanOptions::default(), result: None, history: None, err: None }
    }

    /// `BOPAlgo_Options::SetFuzzyValue` etc.
    pub fn set_options(&mut self, options: BooleanOptions) { self.options = options; }

    /// `BRepAlgoAPI_BuilderShape::Build()`.
    /// Pipeline: DS → PaveFiller → BooleanBuilder → result.
    /// Delegates to [`crate::bop_occt_ops::boolean_op_with_history_generic`] which mirrors
    /// OCCT `BOPAlgo_BOP::Perform` → `BOPAlgo_Builder::PerformInternal1` flow.
    pub fn build(&mut self) -> bool {
        self.result = None; self.err = None; self.history = None;
        match crate::bop_occt_ops::boolean_op_with_history_generic(
            self.op_type, self.shape1, self.shape2,
        ) {
            Ok((t, h)) => {
                self.result = Some(t);
                self.history = Some(h);
                true
            }
            Err(e) => { self.err = Some(e); false }
        }
    }

    /// `BRepAlgoAPI_BuilderShape::Shape()`.
    /// Returns the result shape. Panics if `build()` not called or failed.
    pub fn shape(&self) -> &topods::BRep {
        self.result.as_ref().expect("build() not called or failed")
    }

    /// `BRepAlgoAPI_BuilderShape::History()`.
    pub fn history(&self) -> Option<&BooleanHistory> {
        self.history.as_ref()
    }

    /// `BRepAlgoAPI_Algo::IsDone()`.
    pub fn is_done(&self) -> bool { self.result.is_some() }

    /// `BRepAlgoAPI_Algo::Error()`.
    pub fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
}

/// Section operation — computes intersection curves between two shapes.
///
/// corresponds to `BRepAlgoAPI_Section`.
/// Currently a stub (returns clone of first shape).
pub struct SectionOp<'a> {
    shape1: &'a topods::BRep,
    shape2: &'a topods::BRep,
    result: Option<topods::BRep>,
    err: Option<BooleanError>,
}

impl<'a> SectionOp<'a> {
    /// `BRepAlgoAPI_Section(a, b)` constructor.
    pub fn new(shape1: &'a topods::BRep, shape2: &'a topods::BRep) -> Self {
        Self { shape1, shape2, result: None, err: None }
    }

    /// `BRepAlgoAPI_Section::Build()`.
    /// Stub: does not compute actual section curves.
    pub fn build(&mut self) -> bool {
        self.result = None; self.err = None;
        self.result = Some(self.shape1.clone());
        true
    }

    /// `BRepAlgoAPI_Section::Shape()`.
    pub fn shape(&self) -> &topods::BRep {
        self.result.as_ref().expect("build() not called or failed")
    }

    /// `BRepAlgoAPI_Algo::IsDone()`.
    pub fn is_done(&self) -> bool { self.result.is_some() }
}

// ── Convenience functions ──────────────────────────────────────────────
// No direct OCCT equivalents; thin wrappers for the common three operations.

/// Fuse (Union): combine two shapes into one.
///
/// OCCT-equivalent: `BRepAlgoAPI_Fuse(a, b).Shape()`
pub fn fuse(a: &topods::BRep, b: &topods::BRep) -> Result<topods::BRep, BooleanError> {
    let mut op = BooleanOp::new(BooleanOpType::Union, a, b);
    if op.build() { Ok(op.shape().clone()) } else { Err(op.err.unwrap_or(BooleanError::InvalidOperation)) }
}

/// Common (Intersection): the shared volume of two shapes.
///
/// OCCT-equivalent: `BRepAlgoAPI_Common(a, b).Shape()`
pub fn common(a: &topods::BRep, b: &topods::BRep) -> Result<topods::BRep, BooleanError> {
    let mut op = BooleanOp::new(BooleanOpType::Intersection, a, b);
    if op.build() { Ok(op.shape().clone()) } else { Err(op.err.unwrap_or(BooleanError::InvalidOperation)) }
}

/// Cut (Difference): subtract shape b from shape a.
///
/// OCCT-equivalent: `BRepAlgoAPI_Cut(a, b).Shape()`
pub fn cut(a: &topods::BRep, b: &topods::BRep) -> Result<topods::BRep, BooleanError> {
    let mut op = BooleanOp::new(BooleanOpType::Difference, a, b);
    if op.build() { Ok(op.shape().clone()) } else { Err(op.err.unwrap_or(BooleanError::InvalidOperation)) }
}
