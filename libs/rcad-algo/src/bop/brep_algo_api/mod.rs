use rcad_kernel::topo_shape::Shape;
use crate::bop::algo::builder::{Builder, BooleanError, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::ds::DS;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};

// ── BRepAlgoAPI_Algo ─────────────────────────────────────────────────────
// OCCT: IsDone(), Error(), Warn() — pure interface
pub trait Algo {
    fn is_done(&self) -> bool;
    fn error(&self) -> Option<&BooleanError>;
}

// ── BRepAlgoAPI_BuilderShape ─────────────────────────────────────────────
// OCCT: concrete class with Shape(), result storage
pub struct BuilderShape {
    pub result: Option<Shape>,
    pub err: Option<BooleanError>,
}
impl BuilderShape {
    pub fn shape(&self) -> &Shape { self.result.as_ref().expect("build() not called or failed") }
}
impl Algo for BuilderShape {
    fn is_done(&self) -> bool { self.result.is_some() }
    fn error(&self) -> Option<&BooleanError> { self.err.as_ref() }
}

// ── BRepAlgoAPI_BuilderAlgo ─────────────────────────────────────────────
// OCCT: SetArguments, SetGlue, SetNonDestructive, SetFuzzyValue, Build, Shape
pub struct BuilderAlgo {
    pub bs: BuilderShape,
    pub arguments: Vec<Shape>,
    pub run_parallel: bool,
    pub fuzzy_value: f64,
    pub non_destructive: bool,
    pub glue: i32,
    pub check_inverted: bool,
    pub use_bvh: bool,
}
impl BuilderAlgo {
    pub fn new() -> Self {
        Self {
            bs: BuilderShape { result: None, err: None },
            arguments: Vec::new(),
            run_parallel: false, fuzzy_value: 0.0,
            non_destructive: false, glue: 0, check_inverted: true, use_bvh: false,
        }
    }
    pub fn set_run_parallel(&mut self, b: bool) { self.run_parallel = b; }
    pub fn get_run_parallel(&self) -> bool { self.run_parallel }
    pub fn set_fuzzy_value(&mut self, v: f64) { self.fuzzy_value = v; }
    pub fn get_fuzzy_value(&self) -> f64 { self.fuzzy_value }
    pub fn set_arguments(&mut self, args: Vec<Shape>) { self.arguments = args; }
    pub fn get_arguments(&self) -> &[Shape] { &self.arguments }
    pub fn set_non_destructive(&mut self, b: bool) { self.non_destructive = b; }
    pub fn get_non_destructive(&self) -> bool { self.non_destructive }
    pub fn set_glue(&mut self, g: i32) { self.glue = g; }
    pub fn get_glue(&self) -> i32 { self.glue }
    pub fn set_check_inverted(&mut self, b: bool) { self.check_inverted = b; }
    pub fn get_check_inverted(&self) -> bool { self.check_inverted }
}
impl Algo for BuilderAlgo {
    fn is_done(&self) -> bool { self.bs.is_done() }
    fn error(&self) -> Option<&BooleanError> { self.bs.error() }
}

// ── BooleanOperation — base for Fuse/Common/Cut/Section ────────────────
fn run_build(algo: &BuilderAlgo, op_type: BooleanOpType) -> Result<Shape, BooleanError> {
    if algo.arguments.len() < 2 { return Err(BooleanError::TooFewArguments); }
    let mut filler = PaveFiller::new();
    filler.set_arguments(algo.arguments.clone());
    filler.set_fuzzy_value(algo.fuzzy_value);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
    filler.perform(&a_ps);
    let fuzz = filler.fuzzy_value();
    // builder borrows the DS from filler; both live in the same scope
    let mut builder = Builder::new(filler.ds(), op_type, fuzz);
    builder.my_arguments = filler.ds().arguments.clone();
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

macro_rules! def_bool_op {
    ($name:ident, $op:ident) => {
        pub struct $name { pub algo: BuilderAlgo }
        impl $name {
            pub fn new() -> Self { Self { algo: BuilderAlgo::new() } }
            pub fn from_shapes(s1: Shape, s2: Shape) -> Self {
                let mut s = Self::new(); s.algo.arguments = vec![s1, s2]; s
            }
            pub fn set_arguments(&mut self, args: Vec<Shape>) { self.algo.set_arguments(args); }
            pub fn get_arguments(&self) -> &[Shape] { self.algo.get_arguments() }
            pub fn set_run_parallel(&mut self, b: bool) { self.algo.set_run_parallel(b); }
            pub fn get_run_parallel(&self) -> bool { self.algo.get_run_parallel() }
            pub fn set_fuzzy_value(&mut self, v: f64) { self.algo.set_fuzzy_value(v); }
            pub fn get_fuzzy_value(&self) -> f64 { self.algo.get_fuzzy_value() }
            pub fn set_non_destructive(&mut self, b: bool) { self.algo.set_non_destructive(b); }
            pub fn get_non_destructive(&self) -> bool { self.algo.get_non_destructive() }
            pub fn set_glue(&mut self, g: i32) { self.algo.set_glue(g); }
            pub fn get_glue(&self) -> i32 { self.algo.get_glue() }
            pub fn set_check_inverted(&mut self, b: bool) { self.algo.set_check_inverted(b); }
            pub fn get_check_inverted(&self) -> bool { self.algo.get_check_inverted() }
            // OCCT BRepAlgoAPI_BuilderShape
            pub fn build(&mut self) {
                self.algo.bs.result = None; self.algo.bs.err = None;
                match run_build(&self.algo, BooleanOpType::$op) {
                    Ok(s) => self.algo.bs.result = Some(s),
                    Err(e) => self.algo.bs.err = Some(e),
                }
            }
            pub fn shape(&self) -> &Shape { self.algo.bs.shape() }
        }
        impl Algo for $name {
            fn is_done(&self) -> bool { self.algo.is_done() }
            fn error(&self) -> Option<&BooleanError> { self.algo.error() }
        }
    };
}

def_bool_op!(FuseOp, Union);
def_bool_op!(CommonOp, Intersection);
def_bool_op!(CutOp, Cut);
def_bool_op!(SectionOp, Section);

// ── BRepAlgoAPI_Defeaturing ─────────────────────────────────────────────
pub struct DefeaturingOp {
    pub algo: BuilderAlgo,
    pub faces_to_remove: Vec<Shape>,
}
impl DefeaturingOp {
    pub fn new() -> Self { Self { algo: BuilderAlgo::new(), faces_to_remove: Vec::new() } }
    pub fn add_face_to_remove(&mut self, f: Shape) { self.faces_to_remove.push(f); }
    pub fn build(&mut self) { self.algo.bs.result = self.algo.arguments.first().cloned(); }
    pub fn shape(&self) -> &Shape { self.algo.bs.shape() }
}
impl Algo for DefeaturingOp {
    fn is_done(&self) -> bool { self.algo.bs.is_done() }
    fn error(&self) -> Option<&BooleanError> { self.algo.bs.error() }
}

// ── BRepAlgoAPI_Splitter ────────────────────────────────────────────────
pub struct SplitterOp { pub algo: BuilderAlgo }
impl SplitterOp {
    pub fn new() -> Self { Self { algo: BuilderAlgo::new() } }
    pub fn add_object(&mut self, s: Shape) { self.algo.arguments.push(s); }
    pub fn add_tool(&mut self, s: Shape) { self.algo.arguments.push(s); }
    pub fn build(&mut self) { self.algo.bs.result = self.algo.arguments.first().cloned(); }
    pub fn shape(&self) -> &Shape { self.algo.bs.shape() }
}
impl Algo for SplitterOp {
    fn is_done(&self) -> bool { self.algo.bs.is_done() }
    fn error(&self) -> Option<&BooleanError> { self.algo.bs.error() }
}

// ── Convenience free functions (BRep form) ──────────────────────────────

/// Collect the top-level shapes of a BRep pool (shapes not referenced by any
/// other TShape).  Analogous to OCCT `TopExp::MapShapes(theShape, TopAbs_SHAPE)`
/// over the compound root.
fn brep_top_shapes(brep: &rcad_kernel::BRep) -> Vec<Shape> {
    use rcad_kernel::topods::TShape;
    let mut referenced = vec![false; brep.tshapes.len()];
    fn mark(sr: &Shape, referenced: &mut Vec<bool>) {
        let i = sr.index;
        if i >= referenced.len() || referenced[i] { return; }
        referenced[i] = true;
        match &*sr.data {
            TShape::Solid(sd) => {
                for sh in &sd.shells { mark(sh, referenced); }
                for v in &sd.internal_vertices { mark(v, referenced); }
                for e in &sd.internal_edges { mark(e, referenced); }
            }
            TShape::Shell(sd) => { for f in &sd.faces { mark(f, referenced); } }
            TShape::Face(fd) => {
                mark(&fd.outer_wire, referenced);
                for w in &fd.inner_wires { mark(w, referenced); }
                for v in &fd.internal_vertices { mark(v, referenced); }
            }
            TShape::Wire(wd) => { for e in &wd.edges { mark(e, referenced); } }
            TShape::Edge(ed) => {
                mark(&ed.first, referenced);
                mark(&ed.last, referenced);
            }
            TShape::CompSolid(cs) => { for s in cs { mark(s, referenced); } }
            TShape::Compound(cd) => { for s in cd { mark(s, referenced); } }
            _ => {}
        }
    }
    for ts in &brep.tshapes {
        match ts.as_ref() {
            TShape::Solid(sd) => {
                for sr in &sd.shells { mark(sr, &mut referenced); }
                for sr in &sd.internal_vertices { mark(sr, &mut referenced); }
                for sr in &sd.internal_edges { mark(sr, &mut referenced); }
            }
            TShape::Shell(sd) => { for sr in &sd.faces { mark(sr, &mut referenced); } }
            TShape::Face(fd) => {
                mark(&fd.outer_wire, &mut referenced);
                for w in &fd.inner_wires { mark(w, &mut referenced); }
                for v in &fd.internal_vertices { mark(v, &mut referenced); }
            }
            TShape::Wire(wd) => { for e in &wd.edges { mark(e, &mut referenced); } }
            TShape::Edge(ed) => {
                mark(&ed.first, &mut referenced);
                mark(&ed.last, &mut referenced);
            }
            TShape::CompSolid(shapes) => { for sr in shapes { mark(sr, &mut referenced); } }
            TShape::Compound(shapes) => { for sr in shapes { mark(sr, &mut referenced); } }
            _ => {}
        }
    }
    brep
        .tshapes
        .iter()
        .enumerate()
        .filter(|(i, _)| !referenced[*i])
        .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, rcad_kernel::topods::Orientation::Forward))
        .collect()
}

/// BRep-form build: run the full PaveFiller + Builder pipeline and return the
/// whole result `BRep` pool (not just the root shape).
fn run_build_brep(algo: &BuilderAlgo, op_type: BooleanOpType) -> Result<rcad_kernel::BRep, BooleanError> {
    if algo.arguments.len() < 2 {
        return Err(BooleanError::TooFewArguments);
    }
    let mut filler = PaveFiller::new();
    filler.set_arguments(algo.arguments.clone());
    filler.set_fuzzy_value(algo.fuzzy_value);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
    filler.perform(&a_ps);
    let fuzz = filler.fuzzy_value();
    // builder borrows the DS from filler; both live in the same scope
    let mut builder = Builder::new(filler.ds(), op_type, fuzz);
    builder.my_arguments = filler.ds().arguments.clone();
    builder.build().map_err(|_| BooleanError::InvalidResult("builder failed"))
}

/// OCCT shortcut: `BRepAlgoAPI_Fuse(a, b).Shape()`.
pub fn fuse(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    op.set_arguments([brep_top_shapes(a), brep_top_shapes(b)].concat());
    run_build_brep(&op, BooleanOpType::Union)
}

/// OCCT shortcut: `BRepAlgoAPI_Common(a, b).Shape()`.
pub fn common(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    op.set_arguments([brep_top_shapes(a), brep_top_shapes(b)].concat());
    run_build_brep(&op, BooleanOpType::Intersection)
}

/// OCCT shortcut: `BRepAlgoAPI_Cut(a, b).Shape()`.
pub fn cut(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    op.set_arguments([brep_top_shapes(a), brep_top_shapes(b)].concat());
    run_build_brep(&op, BooleanOpType::Cut)
}

/// OCCT shortcut: `BRepAlgoAPI_Cut21(a, b).Shape()` — `b` minus `a`.
pub fn cut21(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    cut(b, a) // swap args → b - a
}

/// Dispatch a boolean operation by [`BooleanOpType`] (legacy convenience API).
/// OCCT: BRepAlgoAPI_BOP::SetOperation(BOPAlgo_Operation) — COMMON/FUSE/CUT/
/// CUT21/SECTION.
pub fn boolean_op(op: BooleanOpType, a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    match op {
        BooleanOpType::Union => fuse(a, b),
        BooleanOpType::Intersection => common(a, b),
        BooleanOpType::Cut => cut(a, b),
        BooleanOpType::Cut21 => cut21(a, b),
        BooleanOpType::Section => cut(a, b),
        BooleanOpType::Unknown => Err(BooleanError::TooFewArguments),
    }
}

/// Legacy `boolean_op_with_retry` — the current pipeline already includes the
/// OCCT-style retry ladder, so this is a plain dispatch.
pub fn boolean_op_with_retry(op: BooleanOpType, a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    boolean_op(op, a, b)
}
