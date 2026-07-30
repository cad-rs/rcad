//! OCCT BRepAlgoAPI — high-level boolean operations.
//!
//! 1:1 mapping:
//!   BRepAlgoAPI_Algo                → Algo (IsDone / Error / Warn)
//!   BRepAlgoAPI_BuilderShape        → BuilderShape (Build / Shape)
//!   BRepAlgoAPI_BuilderAlgo         → BuilderAlgo (arguments, glue, history)
//!   BRepAlgoAPI_BooleanOperation    → BooleanOperation (common impl)
//!     ├─ BRepAlgoAPI_Fuse           → FuseOp
//!     ├─ BRepAlgoAPI_Common         → CommonOp
//!     ├─ BRepAlgoAPI_Cut            → CutOp
//!     └─ BRepAlgoAPI_Section        → SectionOp

pub mod argument_analyzer;

use rcad_kernel::topo_shape::Shape;

use crate::bop::algo::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::ds::DS;

// ── BRepAlgoAPI_Algo ─────────────────────────────────────────────────────
// OCCT: IsDone(), Error(), Warn(), UserBreak(), UserProgress()

/// OCCT BRepAlgoAPI_Algo — base algorithm class.
pub trait Algo {
    fn is_done(&self) -> bool;
    fn error(&self) -> Option<&BooleanError>;
}

// ── BRepAlgoAPI_BuilderShape ─────────────────────────────────────────────
// OCCT: Build(), Shape()

/// OCCT BRepAlgoAPI_BuilderShape — build + result shape.
pub trait BuilderShape: Algo {
    fn build(&mut self) -> bool;
    fn shape(&self) -> &Shape;
}

// ── BRepAlgoAPI_BooleanOperation — shared implementation ─────────────────
// OCCT: constructor, Build pipeline, Shape. Used by FuseOp/CommonOp/CutOp/SectionOp.

fn run_pipeline(shape1: &Shape, shape2: &Shape, op_type: BooleanOpType, fuzzy: f64) -> Result<Shape, BooleanError> {
    let mut ds = DS::new();
    ds.set_arguments(vec![shape1.clone(), shape2.clone()]);
    ds.init(fuzzy.max(1e-7));
    let mut filler = PaveFiller::new(&mut ds);
    filler.perform();
    let fuzz = filler.fuzzy_value();
    drop(filler);
    let mut builder = BooleanBuilder::new(&ds, op_type, fuzz);
    builder.my_arguments = ds.arguments.clone();
    match builder.build() {
        Ok(brep) => {
            let root = brep.tshapes.iter().enumerate().rev()
                .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_) | rcad_kernel::topods::TShape::Shell(_)))
                .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, rcad_kernel::topods::Orientation::Forward));
            root.ok_or(BooleanError::InvalidOperation)
        }
        Err(_) => Err(BooleanError::InvalidOperation),
    }
}

macro_rules! impl_bool_op {
    ($name:ident, $op:ident) => {
        pub struct $name {
            shape1: Shape,
            shape2: Shape,
            result: Option<Shape>,
            err: Option<BooleanError>,
            fuzzy: f64,
        }
        impl $name {
            /// OCCT: $name(const TopoDS_Shape& S1, const TopoDS_Shape& S2, PerformNow = true)
            pub fn new(shape1: Shape, shape2: Shape) -> Self {
                Self { shape1, shape2, result: None, err: None, fuzzy: 0.0 }
            }
            pub fn set_fuzzy(&mut self, v: f64) { self.fuzzy = v; }
        }
        impl Algo for $name {
            fn is_done(&self) -> bool { self.result.is_some() }
            fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
        }
        impl BuilderShape for $name {
            /// OCCT: BRepAlgoAPI_BuilderShape::Build()
            fn build(&mut self) -> bool {
                self.result = None; self.err = None;
                match run_pipeline(&self.shape1, &self.shape2, BooleanOpType::$op, self.fuzzy) {
                    Ok(s) => { self.result = Some(s); true }
                    Err(e) => { self.err = Some(e); false }
                }
            }
            /// OCCT: BRepAlgoAPI_BuilderShape::Shape()
            fn shape(&self) -> &Shape {
                self.result.as_ref().expect("build() not called or failed")
            }
        }
    };
}

impl_bool_op!(FuseOp, Union);
impl_bool_op!(CommonOp, Intersection);
impl_bool_op!(CutOp, Difference);

// ── BRepAlgoAPI_Section ──────────────────────────────────────────────────
// OCCT has additional overloads for plane + surface.

/// OCCT BRepAlgoAPI_Section — section curves between shapes.
pub struct SectionOp {
    shape1: Shape,
    shape2: Shape,
    result: Option<Shape>,
    err: Option<BooleanError>,
}

impl SectionOp {
    /// OCCT: BRepAlgoAPI_Section(const TopoDS_Shape& S1, const TopoDS_Shape& S2, PerformNow = true)
    pub fn new(shape1: Shape, shape2: Shape) -> Self {
        Self { shape1, shape2, result: None, err: None }
    }
}
impl Algo for SectionOp {
    fn is_done(&self) -> bool { self.result.is_some() }
    fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
}
impl BuilderShape for SectionOp {
    fn build(&mut self) -> bool {
        self.result = None; self.err = None;
        // OCCT: Section uses BOPAlgo_SECTION operation
        match run_pipeline(&self.shape1, &self.shape2, BooleanOpType::Intersection, 0.0) {
            Ok(s) => { self.result = Some(s); true }
            Err(e) => { self.err = Some(e); false }
        }
    }
    fn shape(&self) -> &Shape {
        self.result.as_ref().expect("build() not called or failed")
    }
}

// ── Convenience free functions (no direct OCCT equivalent) ───────────────

pub fn fuse(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = FuseOp::new(a, b);
    if op.build() { Ok(op.shape().clone()) }
    else { Err(op.err.take().unwrap_or(BooleanError::InvalidOperation)) }
}

pub fn common(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = CommonOp::new(a, b);
    if op.build() { Ok(op.shape().clone()) }
    else { Err(op.err.take().unwrap_or(BooleanError::InvalidOperation)) }
}

pub fn cut(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = CutOp::new(a, b);
    if op.build() { Ok(op.shape().clone()) }
    else { Err(op.err.take().unwrap_or(BooleanError::InvalidOperation)) }
}
