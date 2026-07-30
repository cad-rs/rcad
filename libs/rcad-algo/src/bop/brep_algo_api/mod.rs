//! High-level boolean algorithms.
//!
//! Fuse (Union), Common (Intersection), Cut (Difference), Section
//! with `build`/`shape`/`is_done`/`error` pattern.
//! Uses `Shape` for input/output (OCCT: TopoDS_Shape).
//!
//! ## OCCT mapping
//!
//! | Rust (this module)          | OCCT class                        |
//! |------------------------------|-----------------------------------|
//! | `BooleanOptions`             | `BOPAlgo_Options`                 |
//! | `BooleanOp`                  | `BRepAlgoAPI_Fuse` / Common / Cut |
//! | `BooleanOp::build()`         | `BRepAlgoAPI_BuilderShape::Build` |
//! | `BooleanOp::shape()`         | `BRepAlgoAPI_BuilderShape::Shape` |
//! | `BooleanOp::is_done()`       | `BRepAlgoAPI_Algo::IsDone`        |
//! | `SectionOp`                  | `BRepAlgoAPI_Section`             |
//! | `fuse()` / `common()` / `cut()` | convenience wrapper (no OCCT equivalent) |

pub mod argument_analyzer;
pub mod section;

use rcad_kernel::topods;
use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::ds::DS;

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
            fuzzy_value: 0.0,
            parallel: false,
            use_bvh: false,
            glue_enabled: false,
            glue_tolerance: 1e-7,
            track_history: true,
            heal_result: false,
        }
    }
}

impl BooleanOptions {
    pub fn with_fuzzy_value(mut self, value: f64) -> Self {
        self.fuzzy_value = value;
        self
    }
    pub fn with_parallel(mut self, parallel: bool) -> Self {
        self.parallel = parallel;
        self
    }
    pub fn with_bvh(mut self, use_bvh: bool) -> Self {
        self.use_bvh = use_bvh;
        self
    }
    pub fn with_history(mut self, track: bool) -> Self {
        self.track_history = track;
        self
    }
}

/// A boolean operation between two shapes.
///
/// Replaces `BRepAlgoAPI_Fuse`, `BRepAlgoAPI_Common`,
/// `BRepAlgoAPI_Cut` with a single struct + `BooleanOpType` discriminator.
/// API pattern (`new` -> `build` -> `shape` -> `is_done`) matches
/// `BRepAlgoAPI_BuilderShape`.
///
/// After `build()`, call `shape()` for the result or `error()` on failure.
pub struct BooleanOp {
    shape1: Shape,
    shape2: Shape,
    op_type: BooleanOpType,
    options: BooleanOptions,
    result: Option<Shape>,
    err: Option<BooleanError>,
}

impl<'a> BooleanOp {
    /// OCCT: BRepAlgoAPI_Fuse(const TopoDS_Shape& theS1, const TopoDS_Shape& theS2)
    pub fn new(op: BooleanOpType, shape1: Shape, shape2: Shape) -> Self {
        Self {
            shape1,
            shape2,
            op_type: op,
            options: BooleanOptions::default(),
            result: None,
            err: None,
        }
    }

    pub fn set_options(&mut self, options: BooleanOptions) {
        self.options = options;
    }

    /// OCCT: BRepAlgoAPI_BuilderShape::Build()
    pub fn build(&mut self) -> bool {
        self.result = None;
        self.err = None;

        let arg1 = self.shape1.clone();
        let arg2 = self.shape2.clone();
        let mut ds = DS::new();
        ds.set_arguments(vec![arg1, arg2]);
        ds.init(self.options.fuzzy_value.max(1e-7));

        let mut filler = PaveFiller::new(&mut ds);
        filler.perform();
        let fuzz = filler.fuzzy_value();
        drop(filler);
        let mut builder = BooleanBuilder::new(&ds, self.op_type, fuzz);
        builder.my_arguments = ds.arguments.clone();
        // OCCT L330-332: BRepAlgoAPI_BuilderShape::Shape() returns TopoDS_Shape
        match builder.build() {
            Ok(brep) => {
                // Extract root shape from result BRep (OCCT: builder generates TopoDS_Shape directly)
                let root = brep.tshapes.iter().enumerate().rev()
                    .find(|(_, ts)| matches!(ts.as_ref(), topods::TShape::Solid(_) | topods::TShape::Shell(_)))
                    .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, topods::Orientation::Forward));
                self.result = root;
                self.result.is_some()
            }
            Err(_) => { self.err = Some(BooleanError::InvalidOperation); false }
        }
    }

    /// OCCT: BRepAlgoAPI_BuilderShape::Shape()
    pub fn shape(&self) -> &Shape {
        self.result.as_ref().expect("build() not called or failed")
    }

    /// `BRepAlgoAPI_Algo::IsDone()`.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }

    /// `BRepAlgoAPI_Algo::Error()`.
    pub fn error(&self) -> Option<&BooleanError> {
        self.err.as_ref()
    }
}

/// Section operation — computes intersection curves between two shapes.
///
/// Corresponds to `BRepAlgoAPI_Section`.
/// OCCT: BRepAlgoAPI_Section(const TopoDS_Shape& theS1, const TopoDS_Shape& theS2)
pub struct SectionOp {
    shape1: Shape,
    shape2: Shape,
    result: Option<Shape>,
}

impl SectionOp {
    pub fn new(shape1: Shape, shape2: Shape) -> Self {
        Self { shape1, shape2, result: None }
    }
    /// OCCT: BRepAlgoAPI_Section::Build() — stub.
    pub fn build(&mut self) -> bool {
        self.result = self.shape1.clone().into(); true
    }
    pub fn shape(&self) -> &Shape {
        self.result.as_ref().expect("build() not called or failed")
    }
    pub fn is_done(&self) -> bool { self.result.is_some() }
}

// --- Convenience functions ---
// No direct OCCT equivalents; thin wrappers for the common three operations.

/// Fuse (Union): combine two shapes into one.
///
/// OCCT-equivalent: `BRepAlgoAPI_Fuse(a, b).Shape()`
/// OCCT: BRepAlgoAPI_Fuse(a, b).Shape()
pub fn fuse(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = BooleanOp::new(BooleanOpType::Union, a, b);
    if op.build() { Ok(op.shape().clone()) }
    else { Err(op.err.unwrap_or(BooleanError::InvalidOperation)) }
}

/// OCCT: BRepAlgoAPI_Common(a, b).Shape()
pub fn common(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = BooleanOp::new(BooleanOpType::Intersection, a, b);
    if op.build() { Ok(op.shape().clone()) }
    else { Err(op.err.unwrap_or(BooleanError::InvalidOperation)) }
}

/// OCCT: BRepAlgoAPI_Cut(a, b).Shape()
pub fn cut(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = BooleanOp::new(BooleanOpType::Difference, a, b);
    if op.build() { Ok(op.shape().clone()) }
    else { Err(op.err.unwrap_or(BooleanError::InvalidOperation)) }
}
