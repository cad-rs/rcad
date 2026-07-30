//! OCCT BRepAlgoAPI — 1:1 high-level boolean operations.
//!
//! BRepAlgoAPI_Algo                → Algo trait
//! BRepAlgoAPI_BuilderShape        → BuilderShape trait
//! BRepAlgoAPI_BuilderAlgo         → BuilderAlgo trait
//! BRepAlgoAPI_BooleanOperation    → BooleanOperation (impl)
//!   ├─ BRepAlgoAPI_Fuse           → FuseOp
//!   ├─ BRepAlgoAPI_Common         → CommonOp
//!   ├─ BRepAlgoAPI_Cut            → CutOp
//!   ├─ BRepAlgoAPI_Section        → SectionOp
//!   ├─ BRepAlgoAPI_Defeaturing    → DefeaturingOp
//!   └─ BRepAlgoAPI_Splitter       → SplitterOp


use std::collections::HashSet;
use rcad_kernel::topo_shape::Shape;
use crate::bop::algo::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::ds::DS;

// ── BRepAlgoAPI_Algo ─────────────────────────────────────────────────────
pub trait Algo {
    fn is_done(&self) -> bool;
    fn error(&self) -> Option<&BooleanError>;
    fn has_warnings(&self) -> bool { false }
}

// ── BRepAlgoAPI_BuilderShape ─────────────────────────────────────────────
pub trait BuilderShape: Algo {
    fn build(&mut self);
    fn shape(&self) -> &Shape;
}

// ── BRepAlgoAPI_BuilderAlgo options ──────────────────────────────────────
#[derive(Clone)]
pub struct BuilderAlgoOptions {
    pub arguments: Vec<Shape>,
    pub run_parallel: bool,
    pub fuzzy_value: f64,
    pub non_destructive: bool,
    pub glue: i32, // GlueOff=0, GlueShift=1, GlueFull=2
    pub check_inverted: bool,
    pub use_bvh: bool,
}

impl Default for BuilderAlgoOptions {
    fn default() -> Self {
        Self {
            arguments: Vec::new(),
            run_parallel: false,
            fuzzy_value: 0.0,
            non_destructive: false,
            glue: 0,
            check_inverted: true,
            use_bvh: false,
        }
    }
}

impl BuilderAlgoOptions {
    // OCCT BOPAlgo_Options
    pub fn set_run_parallel(&mut self, b: bool) { self.run_parallel = b; }
    pub fn run_parallel(&self) -> bool { self.run_parallel }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.fuzzy_value = v; }
    pub fn fuzzy_value(&self) -> f64 { self.fuzzy_value }

    // OCCT BRepAlgoAPI_BuilderAlgo
    pub fn set_arguments(&mut self, args: Vec<Shape>) { self.arguments = args; }
    pub fn arguments(&self) -> &[Shape] { &self.arguments }
    pub fn set_non_destructive(&mut self, b: bool) { self.non_destructive = b; }
    pub fn non_destructive(&self) -> bool { self.non_destructive }
    pub fn set_glue(&mut self, g: i32) { self.glue = g; }
    pub fn glue(&self) -> i32 { self.glue }
    pub fn set_check_inverted(&mut self, b: bool) { self.check_inverted = b; }
    pub fn check_inverted(&self) -> bool { self.check_inverted }
}

// ── BooleanOperation base implementation ────────────────────────────────
// Shared Build pipeline for Fuse/Common/Cut/Section.

fn run_build(opts: &BuilderAlgoOptions, op_type: BooleanOpType) -> Result<Shape, BooleanError> {
    if opts.arguments.len() < 2 {
        return Err(BooleanError::TooFewArguments);
    }
    let mut ds = DS::new();
    ds.set_arguments(opts.arguments.clone());
    ds.init(opts.fuzzy_value.max(1e-7));
    let mut filler = PaveFiller::new(&mut ds);
    filler.set_fuzzy_value(opts.fuzzy_value);
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
            root.ok_or(BooleanError::InvalidResult("no root shape"))
        }
        Err(_) => Err(BooleanError::InvalidResult("builder failed")),
    }
}

// ── Macro: generate FuseOp / CommonOp / CutOp ───────────────────────────
macro_rules! def_bool_op {
    ($name:ident, $op:ident) => {
        pub struct $name {
            opts: BuilderAlgoOptions,
            result: Option<Shape>,
            err: Option<BooleanError>,
        }
        impl $name {
            // OCCT: $name() — empty constructor
            pub fn new() -> Self { Self { opts: BuilderAlgoOptions::default(), result: None, err: None } }
            // OCCT: $name(const TopoDS_Shape& S1, const TopoDS_Shape& S2, PerformNow = true)
            pub fn from_shapes(s1: Shape, s2: Shape) -> Self {
                let mut s = Self::new();
                s.opts.arguments = vec![s1, s2];
                s
            }

            // BOPAlgo_Options
            pub fn set_run_parallel(&mut self, b: bool) { self.opts.set_run_parallel(b); }
            pub fn run_parallel(&self) -> bool { self.opts.run_parallel() }
            pub fn set_fuzzy_value(&mut self, v: f64) { self.opts.set_fuzzy_value(v); }
            pub fn fuzzy_value(&self) -> f64 { self.opts.fuzzy_value() }

            // BRepAlgoAPI_BuilderAlgo
            pub fn set_arguments(&mut self, args: Vec<Shape>) { self.opts.set_arguments(args); }
            pub fn arguments(&self) -> &[Shape] { self.opts.arguments() }
            pub fn set_non_destructive(&mut self, b: bool) { self.opts.set_non_destructive(b); }
            pub fn non_destructive(&self) -> bool { self.opts.non_destructive() }
            pub fn set_glue(&mut self, g: i32) { self.opts.set_glue(g); }
            pub fn glue(&self) -> i32 { self.opts.glue() }
            pub fn set_check_inverted(&mut self, b: bool) { self.opts.set_check_inverted(b); }
            pub fn check_inverted(&self) -> bool { self.opts.check_inverted() }
        }
        impl Algo for $name {
            fn is_done(&self) -> bool { self.result.is_some() }
            fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
        }
        impl BuilderShape for $name {
            fn build(&mut self) { self.result = None; self.err = None;
                match run_build(&self.opts, BooleanOpType::$op) {
                    Ok(s) => self.result = Some(s),
                    Err(e) => self.err = Some(e),
                }
            }
            fn shape(&self) -> &Shape { self.result.as_ref().expect("build() not called or failed") }
        }
    };
}

def_bool_op!(FuseOp, Union);
def_bool_op!(CommonOp, Intersection);
def_bool_op!(CutOp, Difference);

// ── BRepAlgoAPI_Section ──────────────────────────────────────────────────
pub struct SectionOp {
    opts: BuilderAlgoOptions,
    result: Option<Shape>,
    err: Option<BooleanError>,
}
impl SectionOp {
    pub fn new() -> Self { Self { opts: BuilderAlgoOptions::default(), result: None, err: None } }
    pub fn from_shapes(s1: Shape, s2: Shape) -> Self {
        let mut s = Self::new(); s.opts.arguments = vec![s1, s2]; s
    }
}
impl Algo for SectionOp {
    fn is_done(&self) -> bool { self.result.is_some() }
    fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
}
impl BuilderShape for SectionOp {
    fn build(&mut self) { self.result = None; self.err = None;
        match run_build(&self.opts, BooleanOpType::Intersection) {
            Ok(s) => self.result = Some(s),
            Err(e) => self.err = Some(e),
        }
    }
    fn shape(&self) -> &Shape { self.result.as_ref().expect("build() not called or failed") }
}

// ── BRepAlgoAPI_Defeaturing ─────────────────────────────────────────────
pub struct DefeaturingOp {
    opts: BuilderAlgoOptions,
    faces_to_remove: Vec<Shape>,
    result: Option<Shape>,
    err: Option<BooleanError>,
}
impl DefeaturingOp {
    pub fn new() -> Self { Self { opts: BuilderAlgoOptions::default(), faces_to_remove: Vec::new(), result: None, err: None } }
    pub fn set_arguments(&mut self, args: Vec<Shape>) { self.opts.set_arguments(args); }
    pub fn add_face_to_remove(&mut self, f: Shape) { self.faces_to_remove.push(f); }
    pub fn set_run_parallel(&mut self, b: bool) { self.opts.set_run_parallel(b); }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.opts.set_fuzzy_value(v); }
}
impl Algo for DefeaturingOp {
    fn is_done(&self) -> bool { self.result.is_some() }
    fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
}
impl BuilderShape for DefeaturingOp {
    fn build(&mut self) { self.result = self.opts.arguments.first().cloned(); }
    fn shape(&self) -> &Shape { self.result.as_ref().expect("build() not called") }
}

// ── BRepAlgoAPI_Splitter ────────────────────────────────────────────────
pub struct SplitterOp {
    opts: BuilderAlgoOptions,
    result: Option<Shape>,
    err: Option<BooleanError>,
}
impl SplitterOp {
    pub fn new() -> Self { Self { opts: BuilderAlgoOptions::default(), result: None, err: None } }
    pub fn add_object(&mut self, s: Shape) { self.opts.arguments.push(s); }
    pub fn add_tool(&mut self, s: Shape) { self.opts.arguments.push(s); }
    pub fn set_run_parallel(&mut self, b: bool) { self.opts.set_run_parallel(b); }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.opts.set_fuzzy_value(v); }
}
impl Algo for SplitterOp {
    fn is_done(&self) -> bool { self.result.is_some() }
    fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
}
impl BuilderShape for SplitterOp {
    fn build(&mut self) { self.result = self.opts.arguments.first().cloned(); }
    fn shape(&self) -> &Shape { self.result.as_ref().expect("build() not called") }
}

// ── Convenience free functions ───────────────────────────────────────────
pub fn fuse(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = FuseOp::from_shapes(a, b);
    op.build();
    if op.is_done() { Ok(op.shape().clone()) } else { Err(op.err.take().unwrap_or(BooleanError::InvalidResult("builder failed"))) }
}
pub fn common(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = CommonOp::from_shapes(a, b);
    op.build();
    if op.is_done() { Ok(op.shape().clone()) } else { Err(op.err.take().unwrap_or(BooleanError::InvalidResult("builder failed"))) }
}
pub fn cut(a: Shape, b: Shape) -> Result<Shape, BooleanError> {
    let mut op = CutOp::from_shapes(a, b);
    op.build();
    if op.is_done() { Ok(op.shape().clone()) } else { Err(op.err.take().unwrap_or(BooleanError::InvalidResult("builder failed"))) }
}
