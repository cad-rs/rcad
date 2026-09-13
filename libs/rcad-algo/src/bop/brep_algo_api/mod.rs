use rcad_kernel::topo_shape::Shape;
use crate::bop::algo::builder::{Builder, BooleanError, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::algo::section::BOPAlgoSection;
use crate::bop::algo::section_attribute::SectionAttribute;
use crate::bop::ds::DS;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::topods::{TEdgeData, TFaceData, TShape, TShellData, TSolidData, TWireData};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

// 閳光偓閳光偓 BRepAlgoAPI_Algo 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// OCCT: IsDone(), Error(), Warn() 閳?pure interface
pub trait Algo {
    fn is_done(&self) -> bool;
    fn error(&self) -> Option<&BooleanError>;
}

// 閳光偓閳光偓 BRepAlgoAPI_BuilderShape 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
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

// 閳光偓閳光偓 BRepAlgoAPI_BuilderAlgo 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
// OCCT: SetArguments, SetTools, SetGlue, SetNonDestructive, SetFuzzyValue, Build, Shape
pub struct BuilderAlgo {
    pub bs: BuilderShape,
    pub arguments: Vec<Shape>,
    /// Tools of the operation (BOPAlgo_BOP::myTools). OCCT
    /// BRepAlgoAPI_BooleanOperation stores theS1 in myArguments and theS2 in
    /// myTools; BOPAlgo_BOP::CheckData requires both lists non-empty.
    pub tools: Vec<Shape>,
    /// Merged TopLoc_Location table (index 0 = identity) for `arguments`.
    pub locations: Vec<glam::DAffine3>,
    pub run_parallel: bool,
    pub fuzzy_value: f64,
    pub non_destructive: bool,
    pub glue: i32,
    pub check_inverted: bool,
    pub use_bvh: bool,
    /// OCCT BRepAlgoAPI_BuilderAlgo::myIsIntersectionNeeded
    /// (BRepAlgoAPI_BuilderAlgo.hxx L233): TRUE for the empty constructor,
    /// FALSE for the BRepAlgoAPI_BuilderAlgo(const BOPAlgo_PaveFiller&) form
    /// (BRepAlgoAPI_BuilderAlgo.cxx L30, L43). When FALSE, IntersectShapes
    /// returns at once (L110-113) and the Build steps bind the supplied filler.
    pub my_is_intersection_needed: bool,
    /// OCCT BRepAlgoAPI_BuilderAlgo::myFillHistory
    /// (BRepAlgoAPI_BuilderAlgo.hxx L230): TRUE by default (L29).
    pub my_fill_history: bool,
    /// OCCT BRepAlgoAPI_BuilderAlgo::myHistory (hxx L239) — the general
    /// history tool of the operation (BRepAlgoAPI_Algo::myHistory).
    pub my_history: Option<crate::bop::history::BRepToolsHistory>,
}
impl BuilderAlgo {
    pub fn new() -> Self {
        Self {
            bs: BuilderShape { result: None, err: None },
            arguments: Vec::new(),
            tools: Vec::new(),
            locations: vec![glam::DAffine3::IDENTITY],
            run_parallel: false,
            // OCCT BRepAlgoAPI_BuilderAlgorithm (BOPAlgo_Algo base):
            // myFuzzyValue(Precision::Confusion()) — 1e-7 default.
            fuzzy_value: rcad_kernel::precision::CONFUSION,
            non_destructive: false, glue: 0, check_inverted: true, use_bvh: false,
            my_is_intersection_needed: true,
            my_fill_history: true,
            my_history: None,
        }
    }
    pub fn set_run_parallel(&mut self, b: bool) { self.run_parallel = b; }
    pub fn get_run_parallel(&self) -> bool { self.run_parallel }
    /// OCCT BOPAlgo_Options::SetFuzzyValue (BOPAlgo_Options.cxx L107):
    /// myFuzzyValue = max(theFuzz, Precision::Confusion()).
    pub fn set_fuzzy_value(&mut self, v: f64) {
        self.fuzzy_value = v.max(rcad_kernel::precision::CONFUSION);
    }
    pub fn get_fuzzy_value(&self) -> f64 { self.fuzzy_value }
    pub fn set_arguments(&mut self, args: Vec<Shape>) { self.arguments = args; }
    pub fn get_arguments(&self) -> &[Shape] { &self.arguments }
    pub fn set_tools(&mut self, tools: Vec<Shape>) { self.tools = tools; }
    pub fn get_tools(&self) -> &[Shape] { &self.tools }
    pub fn set_non_destructive(&mut self, b: bool) { self.non_destructive = b; }
    pub fn get_non_destructive(&self) -> bool { self.non_destructive }
    pub fn set_glue(&mut self, g: i32) { self.glue = g; }
    pub fn get_glue(&self) -> i32 { self.glue }
    pub fn set_check_inverted(&mut self, b: bool) { self.check_inverted = b; }
    pub fn get_check_inverted(&self) -> bool { self.check_inverted }

    /// OCCT BRepAlgoAPI_BuilderAlgo::SetToFillHistory (BRepAlgoAPI_BuilderAlgo
    /// .hxx L185).
    pub fn set_to_fill_history(&mut self, b: bool) { self.my_fill_history = b; }

    /// OCCT BRepAlgoAPI_BuilderAlgo::HasHistory (BRepAlgoAPI_BuilderAlgo.hxx
    /// L188).
    pub fn has_history(&self) -> bool { self.my_fill_history }

    // ------------------------------------------------------------------
    // OCCT BRepAlgoAPI_BuilderAlgo(const BOPAlgo_PaveFiller&) + Build
    // (BRepAlgoAPI_BuilderAlgo.cxx L38-47, L81-101, L107-134, L138-164)
    // ------------------------------------------------------------------

    /// OCCT BRepAlgoAPI_BuilderAlgo(const BOPAlgo_PaveFiller& thePF)
    /// (BRepAlgoAPI_BuilderAlgo.cxx L38-47):
    ///   myNonDestructive(false), myGlue(BOPAlgo_GlueOff),
    ///   myCheckInverted(true), myFillHistory(true),
    ///   myIsIntersectionNeeded(false), myBuilder(nullptr),
    ///   myDSFiller = (BOPAlgo_PaveFiller*)&aPF.
    ///
    /// Architecture difference: a Rust value type cannot hold the borrowed
    /// filler, so the filler is handed to [`BuilderAlgo::build_with_filler`]
    /// (the Build step that consumes it) and this constructor carries the
    /// constructor's flag assignment (myIsIntersectionNeeded = false).
    pub fn new_with_filler() -> Self {
        let mut a_builder = Self::new();
        a_builder.my_is_intersection_needed = false;
        a_builder
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::Build (BRepAlgoAPI_BuilderAlgo.cxx
    /// L81-101) over the supplied filler: IntersectShapes returns immediately
    /// (L110-113 — myIsIntersectionNeeded is false), the BOPAlgo_Builder is
    /// created over myArguments (L96-98) and BuildResult (L138-164) binds the
    /// filler through myBuilder->PerformWithFiller(*myDSFiller) — the Builder
    /// pass of BOPAlgo_Builder::PerformInternal1, with NO BuildShape step.
    pub fn build_with_filler(&mut self, the_pf: &PaveFiller) {
        // OCCT L84: NotDone();
        self.bs.result = None;
        self.bs.err = None;
        // OCCT L86: Clear(); — BRepAlgoAPI_BuilderAlgo.cxx L58-77:
        // myHistory.Nullify().
        self.my_history = None;
        // OCCT L95-96: myBuilder = new BOPAlgo_Builder(myAllocator);
        // The rcad Builder borrows the DS of the supplied filler — the same
        // binding BOPAlgo_Builder::PerformInternal1 does (L313-315).
        let mut a_builder = Builder::new(
            the_pf.ds(),
            BooleanOpType::Unknown,
            the_pf.fuzzy_value(),
        );
        // OCCT L98: myBuilder->SetArguments(myArguments);
        // OCCT BOPAlgo_Builder::SetArguments (BOPAlgo_Builder.cxx L113-125)
        // appends under the myMapFence; myArguments is the DS argument list
        // (its constructor form has no myTools).
        a_builder.my_arguments = the_pf.ds().arguments.clone();
        // OCCT L141-144: SetRunParallel / SetCheckInverted / SetToFillHistory.
        a_builder.my_run_parallel = self.run_parallel;
        a_builder.my_check_inverted = self.check_inverted;
        a_builder.my_fill_history = self.my_fill_history;
        // OCCT L146: myBuilder->PerformWithFiller(*myDSFiller, theRange);
        a_builder.perform_with_filler(the_pf);
        // OCCT L148: GetReport()->Merge(myBuilder->GetReport());
        // rcad carries the builder's report in the Builder; the rcad
        // BuilderAlgo has no report of its own (interface gap of this facade).
        // OCCT L150-153: if (myBuilder->HasErrors()) return;
        if a_builder.has_errors() {
            self.bs.err = Some(BooleanError::InvalidResult("builder failed"));
            return;
        }
        // OCCT L155: Done();
        // OCCT L157: myShape = myBuilder->Shape();
        self.bs.result = Some(builder_pass_root(&a_builder));
        // OCCT L159-163: if (myFillHistory) { myHistory = new
        // BRepTools_History; myHistory->Merge(myBuilder->History()); }
        if self.my_fill_history {
            let mut a_history = crate::bop::history::BRepToolsHistory::new();
            if let Some(a_builder_history) = a_builder.history() {
                a_history.merge(a_builder_history);
            }
            self.my_history = Some(a_history);
        }
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::Modified (BRepAlgoAPI_BuilderAlgo.cxx
    /// L204-212).
    pub fn modified(&self, the_s: &Shape) -> Vec<Shape> {
        if self.my_fill_history {
            if let Some(a_history) = &self.my_history {
                return a_history.modified(the_s);
            }
        }
        Vec::new()
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::Generated (BRepAlgoAPI_BuilderAlgo.cxx
    /// L216-224).
    pub fn generated(&self, the_s: &Shape) -> Vec<Shape> {
        if self.my_fill_history {
            if let Some(a_history) = &self.my_history {
                return a_history.generated(the_s);
            }
        }
        Vec::new()
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::IsDeleted (BRepAlgoAPI_BuilderAlgo.cxx
    /// L228-231).
    pub fn is_deleted(&self, the_s: &Shape) -> bool {
        if self.my_fill_history {
            if let Some(a_history) = &self.my_history {
                return a_history.is_removed(the_s);
            }
        }
        false
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::HasModified (BRepAlgoAPI_BuilderAlgo.cxx
    /// L235-238).
    pub fn has_modified(&self) -> bool {
        if self.my_fill_history {
            if let Some(a_history) = &self.my_history {
                return a_history.has_modified();
            }
        }
        false
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::HasGenerated (BRepAlgoAPI_BuilderAlgo.cxx
    /// L242-245).
    pub fn has_generated(&self) -> bool {
        if self.my_fill_history {
            if let Some(a_history) = &self.my_history {
                return a_history.has_generated();
            }
        }
        false
    }

    /// OCCT BRepAlgoAPI_BuilderAlgo::HasDeleted (BRepAlgoAPI_BuilderAlgo.cxx
    /// L249-252).
    pub fn has_deleted(&self) -> bool {
        if self.my_fill_history {
            if let Some(a_history) = &self.my_history {
                return a_history.has_removed();
            }
        }
        false
    }
}
impl Algo for BuilderAlgo {
    fn is_done(&self) -> bool { self.bs.is_done() }
    fn error(&self) -> Option<&BooleanError> { self.bs.error() }
}

/// OCCT history read surface of the rcad BRepAlgoAPI_BuilderAlgo facade (the
/// `TheAlgo&` type parameter of the BRepTools_History template constructor /
/// Merge template, BRepTools_History.hxx L99-132, L212-217).
impl crate::bop::history::HistoryAlgo for BuilderAlgo {
    fn is_deleted(&self, the_s: &Shape) -> bool {
        BuilderAlgo::is_deleted(self, the_s)
    }

    fn modified(&self, the_s: &Shape) -> Vec<Shape> {
        BuilderAlgo::modified(self, the_s)
    }

    fn generated(&self, the_s: &Shape) -> Vec<Shape> {
        BuilderAlgo::generated(self, the_s)
    }
}

/// OCCT BOPAlgo_BuilderShape::Shape() for the Builder pass — the compound made
/// by Prepare (BOPAlgo_Builder.cxx L156-164) and filled by the BuildResult
/// steps (BOPAlgo_Builder_1.cxx L130-168), i.e. the splits of all arguments.
///
/// Architecture difference: the rcad Builder keeps a flat BRep pool instead of
/// the compound, so the root is rebuilt from the pool's top-level shapes — the
/// shapes no other pool shape references, in pool (= OCCT Add) order.
fn builder_pass_root(the_builder: &Builder) -> Shape {
    let Some(brep) = the_builder.my_shape.as_ref() else {
        return Shape::null();
    };
    Shape::new(
        Arc::new(TShape::Compound(brep_top_shapes(brep))),
        0,
        rcad_kernel::topods::Orientation::Forward,
    )
}

// 閳光偓閳光偓 BooleanOperation 閳?base for Fuse/Common/Cut/Section 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
/// OCCT BRepAlgoAPI_BooleanOperation::Build: aLArgs = myArguments + myTools
/// combined for the intersection; the builder receives myArguments (objects)
/// and myTools separately (BOPAlgo_BOP::SetArguments/SetTools).
fn run_build(algo: &BuilderAlgo, op_type: BooleanOpType) -> Result<Shape, BooleanError> {
    if algo.arguments.len() < 1 || algo.tools.len() < 1 {
        return Err(BooleanError::TooFewArguments);
    }
    let mut all_args = algo.arguments.clone();
    all_args.extend(algo.tools.iter().cloned());
    let mut filler = PaveFiller::new();
    filler.set_arguments(all_args);
    filler.set_fuzzy_value(algo.fuzzy_value);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
    filler.perform(&a_ps);
    let fuzz = filler.fuzzy_value();
    // builder borrows the DS from filler; both live in the same scope
    let mut builder = Builder::new(filler.ds(), op_type, fuzz);
    // OCCT BOPAlgo_BOP holds the tools by the same TShape identity as the DS
    // (SetTools appends them to myArguments, which the DS references). rcad's
    // DS deep-clones the inputs, so the builder's tool list must carry the
    // DS-cloned shapes (the tail of ds.arguments) — the original BRep shapes
    // have different TShape identities and would never match the split
    // results in BuildRC (bcommon_simple G9: common of box and contained
    // prism).
    builder.my_arguments = filler.ds().arguments.clone();
    let n_tools = algo.tools.len();
    let n_objs = builder.my_arguments.len().saturating_sub(n_tools);
    builder.my_tools = builder.my_arguments[n_objs..].to_vec();
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
                // OCCT BRepAlgoAPI_BooleanOperation(S1, S2, op):
                // myArguments.Append(theS1); myTools.Append(theS2);
                let mut s = Self::new();
                s.algo.arguments = vec![s1];
                s.algo.tools = vec![s2];
                s
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

// 閳光偓閳光偓 BRepAlgoAPI_Section 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
/// OCCT BRepAlgoAPI_Section — the true SECTION operation, driven by
/// BOPAlgo_Section (BOPAlgo_Section.cxx). Replaces the former degraded path
/// that mapped Section onto Cut.
pub struct SectionOp {
    pub algo: BuilderAlgo,
    /// OCCT BRepAlgoAPI_Section::myApprox (Init, BRepAlgoAPI_Section.cxx
    /// L146): default false.
    pub my_approx: bool,
    /// OCCT myComputePCurveOn1 (Init L147): default false.
    pub my_compute_pcurve1: bool,
    /// OCCT myComputePCurveOn2 (Init L148): default false.
    pub my_compute_pcurve2: bool,
}
impl SectionOp {
    pub fn new() -> Self {
        // OCCT BRepAlgoAPI_Section::Init (L142-153): myApprox =
        // myComputePCurve1 = myComputePCurve2 = false.
        Self {
            algo: BuilderAlgo::new(),
            my_approx: false,
            my_compute_pcurve1: false,
            my_compute_pcurve2: false,
        }
    }
    pub fn from_shapes(s1: Shape, s2: Shape) -> Self {
        // OCCT BRepAlgoAPI_Section(Sh1, Sh2, PerformNow) :
        // BRepAlgoAPI_BooleanOperation(Sh1, Sh2, BOPAlgo_SECTION):
        // myArguments.Append(theS1); myTools.Append(theS2);
        let mut s = Self::new();
        s.algo.arguments = vec![s1];
        s.algo.tools = vec![s2];
        s
    }

    /// OCCT BRepAlgoAPI_Section(S1, S2, thePF) — the filler-carrying section
    /// (BRepAlgoAPI_Section.hxx L52 + BRepAlgoAPI_BooleanOperation.cxx
    /// L111-118) which delegates to BRepAlgoAPI_BuilderAlgo(thePF)
    /// (BRepAlgoAPI_BuilderAlgo.cxx L38-47): myIsIntersectionNeeded = false,
    /// so Build() skips IntersectShapes (L110-113) and builds the section over
    /// the DS the supplied filler already computed
    /// (BRepAlgoAPI_BooleanOperation.cxx L199-200: myBuilder = new
    /// BOPAlgo_Section(myAllocator); myBuilder->SetArguments(
    /// myDSFiller->Arguments())).
    ///
    /// Architecture difference: the borrowed filler cannot be stored in the
    /// value type, so it is handed to [`SectionOp::build_with_filler`].
    pub fn from_shapes_with_filler(s1: Shape, s2: Shape) -> Self {
        let mut s = Self::from_shapes(s1, s2);
        s.algo.my_is_intersection_needed = false;
        s
    }

    /// OCCT BRepAlgoAPI_Section::Build over a supplied filler — the common
    /// BRepAlgoAPI_BuilderAlgo::Build (BRepAlgoAPI_BuilderAlgo.cxx L81-101)
    /// with myIsIntersectionNeeded = false and, for BOPAlgo_SECTION, the
    /// builder created over myDSFiller->Arguments()
    /// (BRepAlgoAPI_BooleanOperation.cxx L199-200).
    pub fn build_with_filler(&mut self, the_pf: &PaveFiller) {
        // OCCT BRepAlgoAPI_BuilderAlgo::Build L84-86: NotDone(); Clear();
        self.algo.bs.result = None;
        self.algo.bs.err = None;
        self.algo.my_history = None;
        // OCCT BRepAlgoAPI_BooleanOperation::Build L174-193: with
        // myIsIntersectionNeeded = false the intersection step is skipped —
        // IntersectShapes returns at once (BRepAlgoAPI_BuilderAlgo.cxx
        // L110-113) and the filler's DS is the one already computed.
        match run_build_section_brep_with_filler(&self.algo, the_pf) {
            Ok(brep) => {
                // OCCT BOPAlgo_Section::myShape is the result compound.
                let root = brep.tshapes.iter().enumerate().rev()
                    .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Compound(_)))
                    .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, rcad_kernel::topods::Orientation::Forward));
                match root {
                    Some(s) => self.algo.bs.result = Some(s),
                    None => self.algo.bs.err = Some(BooleanError::InvalidResult("no root shape")),
                }
            }
            Err(e) => self.algo.bs.err = Some(e),
        }
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
    // OCCT BRepAlgoAPI_Section::Approximation (L171-174).
    pub fn approximation(&mut self, b: bool) { self.my_approx = b; }
    // OCCT BRepAlgoAPI_Section::ComputePCurveOn1 (L176-179).
    pub fn compute_pcurve_on1(&mut self, b: bool) { self.my_compute_pcurve1 = b; }
    // OCCT BRepAlgoAPI_Section::ComputePCurveOn2 (L181-184).
    pub fn compute_pcurve_on2(&mut self, b: bool) { self.my_compute_pcurve2 = b; }
    // OCCT BRepAlgoAPI_BuilderShape
    pub fn build(&mut self) {
        self.algo.bs.result = None; self.algo.bs.err = None;
        match run_build_section_brep(&self.algo, self.my_approx, self.my_compute_pcurve1, self.my_compute_pcurve2)
        {
            Ok(brep) => {
                // OCCT BOPAlgo_Section::myShape is the result compound.
                let root = brep.tshapes.iter().enumerate().rev()
                    .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Compound(_)))
                    .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, rcad_kernel::topods::Orientation::Forward));
                match root {
                    Some(s) => self.algo.bs.result = Some(s),
                    None => self.algo.bs.err = Some(BooleanError::InvalidResult("no root shape")),
                }
            }
            Err(e) => self.algo.bs.err = Some(e),
        }
    }
    pub fn shape(&self) -> &Shape { self.algo.bs.shape() }
}
impl Algo for SectionOp {
    fn is_done(&self) -> bool { self.algo.is_done() }
    fn error(&self) -> Option<&BooleanError> { self.algo.error() }
}

// 閳光偓閳光偓 BRepAlgoAPI_Defeaturing 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
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

// 閳光偓閳光偓 BRepAlgoAPI_Splitter 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓
pub struct SplitterOp { pub algo: BuilderAlgo }
impl SplitterOp {
    pub fn new() -> Self { Self { algo: BuilderAlgo::new() } }
    // OCCT BRepAlgoAPI_Splitter::SetArguments (objects) — myArguments.
    pub fn add_object(&mut self, s: Shape) { self.algo.arguments.push(s); }
    // OCCT BRepAlgoAPI_Splitter::AddTool (tools) — myTools.
    pub fn add_tool(&mut self, s: Shape) { self.algo.tools.push(s); }
    // OCCT BRepAlgoAPI_Splitter::Build (BRepAlgoAPI_Splitter.cxx L35-76):
    // aLArgs = myArguments + myTools for the intersection; the builder
    // receives myArguments (objects) and myTools separately
    // (BOPAlgo_Splitter::SetArguments/SetTools).
    pub fn build(&mut self) {
        self.algo.bs.result = None; self.algo.bs.err = None;
        match run_build_splitter_brep(&self.algo) {
            Ok(brep) => {
                let root = brep.tshapes.iter().enumerate().rev()
                    .find(|(_, ts)| matches!(ts.as_ref(), rcad_kernel::topods::TShape::Solid(_) | rcad_kernel::topods::TShape::Shell(_)))
                    .map(|(i, ts)| Shape::from_parts(ts.clone(), i, 0, rcad_kernel::topods::Orientation::Forward));
                self.algo.bs.result = root;
                if self.algo.bs.result.is_none() {
                    self.algo.bs.err = Some(BooleanError::InvalidResult("no root shape"));
                }
            }
            Err(e) => self.algo.bs.err = Some(e),
        }
    }
    pub fn shape(&self) -> &Shape { self.algo.bs.shape() }
}
impl Algo for SplitterOp {
    fn is_done(&self) -> bool { self.algo.bs.is_done() }
    fn error(&self) -> Option<&BooleanError> { self.algo.bs.error() }
}

// 閳光偓閳光偓 Convenience free functions (BRep form) 閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓閳光偓

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

/// Recursively remap every `Shape.location` index through `map` (index 0,
/// identity, is preserved). Shared TShapes are rebuilt once per TShape pointer
/// via the cache, mirroring `clone_arguments_private`.
fn remap_location_tree(
    sr: &Shape,
    map: &HashMap<u32, u32>,
    cache: &mut HashMap<u64, Arc<TShape>>,
) -> Shape {
    let new_loc = if sr.location != 0 {
        map.get(&sr.location).copied().unwrap_or(sr.location)
    } else {
        0
    };
    let ptr = sr.ptr_id();
    if let Some(ts) = cache.get(&ptr) {
        return Shape {
            data: ts.clone(),
            index: sr.index,
            location: new_loc,
            orientation: sr.orientation,
        };
    }
    let new_ts = match &*sr.data {
        TShape::Vertex(vd) => TShape::Vertex(vd.clone()),
        TShape::Edge(ed) => {
            let first = remap_location_tree(&ed.first, map, cache);
            let last = remap_location_tree(&ed.last, map, cache);
            // Phase 1 only: clone identity-keyed maps verbatim. An edge's
            // owning face is its ANCESTOR in this walk, so it is never in the
            // cache yet at this point; keys are rewritten in a second pass
            // (rewrite_identity_keys) once every pointer is final.
            TShape::Edge(TEdgeData {
                first,
                last,
                ..ed.clone()
            })
        }
        TShape::Wire(wd) => TShape::Wire(TWireData {
            edges: wd
                .edges
                .iter()
                .map(|e| remap_location_tree(e, map, cache))
                .collect(),
            ..wd.clone()
        }),
        TShape::Face(fd) => TShape::Face(TFaceData {
            outer_wire: remap_location_tree(&fd.outer_wire, map, cache),
            inner_wires: fd
                .inner_wires
                .iter()
                .map(|w| remap_location_tree(w, map, cache))
                .collect(),
            internal_vertices: fd
                .internal_vertices
                .iter()
                .map(|v| remap_location_tree(v, map, cache))
                .collect(),
            ..fd.clone()
        }),
        TShape::Shell(sd) => TShape::Shell(TShellData {
            faces: sd
                .faces
                .iter()
                .map(|f| remap_location_tree(f, map, cache))
                .collect(),
            ..sd.clone()
        }),
        TShape::Solid(sd) => TShape::Solid(TSolidData {
            shells: sd
                .shells
                .iter()
                .map(|s| remap_location_tree(s, map, cache))
                .collect(),
            internal_vertices: sd
                .internal_vertices
                .iter()
                .map(|v| remap_location_tree(v, map, cache))
                .collect(),
            internal_edges: sd
                .internal_edges
                .iter()
                .map(|e| remap_location_tree(e, map, cache))
                .collect(),
            ..sd.clone()
        }),
        TShape::CompSolid(cd) => TShape::CompSolid(
            cd.iter()
                .map(|s| remap_location_tree(s, map, cache))
                .collect(),
        ),
        TShape::Compound(cd) => TShape::Compound(
            cd.iter()
                .map(|s| remap_location_tree(s, map, cache))
                .collect(),
        ),
    };
    let new_ts = Arc::new(new_ts);
    cache.insert(ptr, new_ts.clone());
    Shape {
        data: new_ts,
        index: sr.index,
        location: new_loc,
        orientation: sr.orientation,
    }
}

/// `brep_top_shapes` + merge this BRep's location table into `global_locs`
/// (appending each entry and recording old-index 閳?new-index) and remap every
/// returned shape's `location` to the merged table. Index 0 (identity) is
/// shared; BRep location tables start at index 1.
///
/// `cache` is SHARED across every argument/tools call of one boolean
/// operation: OCCT's BRepAlgoAPI_BooleanOperation::SetArguments/SetTools hold
/// TopoDS_Shape handles, so a TShape shared between an argument and a tool
/// (the BRepSweep prism of a face taken from the other argument) keeps ONE
/// identity in the DS.  Remapping each BRep with a per-call cache would
/// duplicate the shared TShape and break the vertex-count semantics
/// (boptuc_simple ZP3: the prism's located vertices must stay identical to
/// the cone's).
fn brep_top_shapes_with_locations(
    brep: &rcad_kernel::BRep,
    global_locs: &mut Vec<glam::DAffine3>,
    cache: &mut HashMap<u64, Arc<TShape>>,
) -> Vec<Shape> {
    let mut map: HashMap<u32, u32> = HashMap::new();
    for (i, loc) in brep.locations.iter().enumerate() {
        let old = (i + 1) as u32; // BRep table index (0 = identity)
        // OCCT TopLoc_Location items are deduplicated by transformation:
        // identical Trsf matrices share one TopLoc_Datum3D item
        // (TopLoc_Location::Location() hash/IsEqual), so two shapes built
        // with the same translation resolve to the same table index. Merge
        // only when the transform is absent from the global table.
        let new = match global_locs.iter().position(|l| *l == *loc) {
            Some(existing) => existing as u32,
            None => {
                global_locs.push(*loc);
                (global_locs.len() - 1) as u32
            }
        };
        map.insert(old, new);
    }
    // Both BReps go through the same cache even when one has no located
    // sub-shapes (an empty map remaps locations identically): the cache hit
    // on a TShape shared with the other BRep returns the SAME remapped Arc,
    // preserving the cross-argument TShape identity.
    let tops: Vec<Shape> = brep_top_shapes(brep)
        .into_iter()
        .map(|s| remap_location_tree(&s, &map, cache))
        .collect();
    // Second pass: all pointers are final now; rewrite every edge's identity
    // keys against the complete cache (old ptr -> new Arc ptr).
    rewrite_identity_keys(&tops, cache);
    tops
}

/// Rewrite the face-pointer identity keys of every edge reachable from the
/// rebuilt top shapes using the completed clone cache.  In-place on the shared
/// Arcs; unknown owners keep their pointer.
fn rewrite_identity_keys(tops: &[Shape], cache: &HashMap<u64, Arc<TShape>>) {
    let mut visited: HashSet<u64> = HashSet::new();
    let mut stack: Vec<Shape> = tops.to_vec();
    while let Some(sh) = stack.pop() {
        if !visited.insert(sh.ptr_id()) {
            continue;
        }
        match &*sh.data {
            TShape::Edge(ed) => {
                let raw = Arc::as_ptr(&sh.data) as *mut TShape;
                // SAFETY: single-threaded build; no other &TShape borrow is
                // alive at this point.
                unsafe {
                    if let TShape::Edge(edm) = &mut *raw {
                        edm.pcurves = ed
                            .pcurves
                            .iter()
                            .map(|(&(p, l), v)| {
                                let np =
                                    cache.get(&p).map(|a| Arc::as_ptr(a) as u64).unwrap_or(p);
                                ((np, l), v.clone())
                            })
                            .collect();
                        edm.representations = ed
                            .representations
                            .iter()
                            .map(|r| match r {
                                rcad_kernel::topods::CurveRepresentation::CurveOnSurface { face, pcurve, range } => {
                                    rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                                        face: (
                                            cache.get(&face.0).map(|a| Arc::as_ptr(a) as u64).unwrap_or(face.0),
                                            face.1,
                                        ),
                                        pcurve: pcurve.clone(),
                                        range: *range,
                                    }
                                }
                                rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface { face, pcurve1, pcurve2, range } => {
                                    rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                                        face: (
                                            cache.get(&face.0).map(|a| Arc::as_ptr(a) as u64).unwrap_or(face.0),
                                            face.1,
                                        ),
                                        pcurve1: pcurve1.clone(),
                                        pcurve2: pcurve2.clone(),
                                        range: *range,
                                    }
                                }
                                other => other.clone(),
                            })
                            .collect();
                        edm.vertex_params = ed
                            .vertex_params
                            .iter()
                            .map(|(&k, &v)| {
                                let nk =
                                    cache.get(&k).map(|a| Arc::as_ptr(a) as u64).unwrap_or(k);
                                (nk, v)
                            })
                            .collect();
                    }
                }
            }
            TShape::Wire(wd) => stack.extend(wd.edges.iter().cloned()),
            TShape::Face(fd) => {
                stack.push(fd.outer_wire.clone());
                stack.extend(fd.inner_wires.iter().cloned());
                stack.extend(fd.internal_vertices.iter().cloned());
            }
            TShape::Shell(sd) => stack.extend(sd.faces.iter().cloned()),
            TShape::Solid(sd) => {
                stack.extend(sd.shells.iter().cloned());
                stack.extend(sd.internal_vertices.iter().cloned());
                stack.extend(sd.internal_edges.iter().cloned());
            }
            TShape::CompSolid(cd) => stack.extend(cd.iter().cloned()),
            TShape::Compound(cd) => stack.extend(cd.iter().cloned()),
            TShape::Vertex(_) => {}
        }
    }
}

/// BRep-form build: run the full PaveFiller + Builder pipeline and return the
/// whole result `BRep` pool (not just the root shape).
/// OCCT BRepAlgoAPI_BooleanOperation::Build: aLArgs = myArguments + myTools
/// combined for the intersection; the builder receives myArguments (objects)
/// and myTools separately (BOPAlgo_BOP::SetArguments/SetTools).
fn run_build_brep(algo: &BuilderAlgo, op_type: BooleanOpType) -> Result<rcad_kernel::BRep, BooleanError> {
    if algo.arguments.len() < 1 || algo.tools.len() < 1 {
        return Err(BooleanError::TooFewArguments);
    }
    let mut all_args = algo.arguments.clone();
    all_args.extend(algo.tools.iter().cloned());
    let mut filler = PaveFiller::new();
    filler.set_arguments(all_args);
    filler.ds_mut().set_locations(algo.locations.clone());
    filler.set_fuzzy_value(algo.fuzzy_value);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
    filler.perform(&a_ps);
    let fuzz = filler.fuzzy_value();
    // builder borrows the DS from filler; both live in the same scope
    let mut builder = Builder::new(filler.ds(), op_type, fuzz);
    // OCCT BOPAlgo_BOP holds the tools by the same TShape identity as the DS
    // (SetTools appends them to myArguments, which the DS references). rcad's
    // DS deep-clones the inputs, so the builder's tool list must carry the
    // DS-cloned shapes (the tail of ds.arguments) — the original BRep shapes
    // have different TShape identities and would never match the split
    // results in BuildRC (bcommon_simple G9: common of box and contained
    // prism).
    builder.my_arguments = filler.ds().arguments.clone();
    let n_tools = algo.tools.len();
    let n_objs = builder.my_arguments.len().saturating_sub(n_tools);
    builder.my_tools = builder.my_arguments[n_objs..].to_vec();
    builder.build().map_err(|_| BooleanError::InvalidResult("builder failed"))
}

/// BRep-form splitter build.
/// OCCT BOPAlgo_Splitter::Perform (BOPAlgo_Splitter.cxx L54-93): aLS =
/// myArguments (objects) + myTools (tools) combined into ONE PaveFiller, then
/// PerformInternal -> BOPAlgo_Builder::PerformInternal1 (GF pipeline, no
/// BuildShape).  BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx L130-168)
/// iterates myArguments (objects only), so only the split parts of the
/// OBJECTS enter the result; tool split parts are excluded.
fn run_build_splitter_brep(algo: &BuilderAlgo) -> Result<rcad_kernel::BRep, BooleanError> {
    // OCCT BRepAlgoAPI_Splitter::Build (BRepAlgoAPI_Splitter.cxx L42-46).
    if algo.arguments.is_empty() || (algo.arguments.len() + algo.tools.len()) < 2 {
        return Err(BooleanError::TooFewArguments);
    }
    // OCCT BOPAlgo_Splitter::Perform L64-77: aLS = myArguments + myTools.
    let mut all_args = algo.arguments.clone();
    all_args.extend(algo.tools.iter().cloned());
    let mut filler = PaveFiller::new();
    filler.set_arguments(all_args);
    filler.ds_mut().set_locations(algo.locations.clone());
    filler.set_fuzzy_value(algo.fuzzy_value);
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
    filler.perform(&a_ps);
    let fuzz = filler.fuzzy_value();
    // builder borrows the DS from filler; both live in the same scope
    let mut builder = Builder::new(filler.ds(), BooleanOpType::Union, fuzz);
    // OCCT BRepAlgoAPI_Splitter::Build L71-72: myBuilder->SetArguments
    // (objects); SetTools(tools).  rcad's DS deep-clones the inputs, so the
    // builder's argument list must carry the DS-cloned shapes: the first
    // n_objs entries of ds.arguments (objects), the rest are tools.
    let n_objs = algo.arguments.len();
    builder.my_arguments = filler.ds().arguments[..n_objs].to_vec();
    builder.my_tools = filler.ds().arguments[n_objs..].to_vec();
    builder.my_is_splitter = true;
    builder.build().map_err(|_| BooleanError::InvalidResult("splitter failed"))
}

/// BRep-form SECTION build — the true BOPAlgo_Section pipeline (replaces the
/// former degraded path that mapped Section onto Cut).
///
/// OCCT BRepAlgoAPI_BooleanOperation::Build (BRepAlgoAPI_BooleanOperation.cxx
/// L177-184): aLArgs = myArguments + myTools -> IntersectShapes (the
/// PaveFiller); then L199-200: myBuilder = new BOPAlgo_Section;
/// myBuilder->SetArguments(myDSFiller->Arguments()) — BOPAlgo_Section
/// inherits BOPAlgo_Builder and has no myTools: objects and tools form ONE
/// argument list.
fn run_build_section_brep(
    algo: &BuilderAlgo,
    my_approx: bool,
    my_compute_pcurve1: bool,
    my_compute_pcurve2: bool,
) -> Result<rcad_kernel::BRep, BooleanError> {
    if algo.arguments.is_empty() || algo.tools.is_empty() {
        return Err(BooleanError::TooFewArguments);
    }
    let mut all_args = algo.arguments.clone();
    all_args.extend(algo.tools.iter().cloned());
    let mut filler = PaveFiller::new();
    filler.set_arguments(all_args);
    filler.ds_mut().set_locations(algo.locations.clone());
    filler.set_fuzzy_value(algo.fuzzy_value);
    // OCCT BRepAlgoAPI_Section::Init (BRepAlgoAPI_Section.cxx L143-152):
    // myApprox = myComputePCurve1 = myComputePCurve2 = false (the caller
    // flags override them); SetAttributes (L196-200):
    // myDSFiller->SetSectionAttribute(BOPAlgo_SectionAttribute(myApprox,
    // myComputePCurve1, myComputePCurve2)).
    filler.my_section_attribute = SectionAttribute {
        approximation: my_approx,
        pcurve_on_s1: my_compute_pcurve1,
        pcurve_on_s2: my_compute_pcurve2,
        ..Default::default()
    };
    let a_prog = NoopProgress;
    let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
    filler.perform(&a_ps);
    let fuzz = filler.fuzzy_value();
    // OCCT BRepAlgoAPI_BooleanOperation::Build L199-200: BOPAlgo_Section over
    // the DS arguments (objects + tools as one list).
    let mut a_section = BOPAlgoSection::new(filler.ds(), fuzz);
    a_section.set_arguments(filler.ds().arguments.clone());
    a_section.perform();
    if a_section.has_errors() {
        return Err(BooleanError::InvalidResult("section failed"));
    }
    a_section
        .result_brep()
        .ok_or(BooleanError::InvalidResult("no section result"))
}

/// BRep-form SECTION build over a supplied filler — the
/// BRepAlgoAPI_BuilderAlgo(const BOPAlgo_PaveFiller&) form of the section
/// pipeline (BRepAlgoAPI_BuilderAlgo.cxx L38-47: myIsIntersectionNeeded =
/// false; L107-113: IntersectShapes returns immediately).
///
/// OCCT BRepAlgoAPI_BooleanOperation::Build L199-200:
///   myBuilder = new BOPAlgo_Section(myAllocator);
///   myBuilder->SetArguments(myDSFiller->Arguments());
/// BOPAlgo_Section inherits BOPAlgo_Builder and has no myTools: the objects and
/// the tools form ONE argument list — the filler's argument list.
///
/// The section attribute is NOT re-applied: the filler was performed by the
/// caller with its own attribute (the OCCT SetAttributes step belongs to
/// IntersectShapes, which is skipped here).
fn run_build_section_brep_with_filler(
    algo: &BuilderAlgo,
    the_pf: &PaveFiller,
) -> Result<rcad_kernel::BRep, BooleanError> {
    let _ = algo;
    // OCCT BOPAlgo_Section::CheckData requires at least one argument
    // (BOPAlgo_Builder::CheckData, BOPAlgo_Builder.cxx L130-140 with the
    // Section form of BOPAlgo_Section::CheckData).
    if the_pf.ds().arguments.is_empty() {
        return Err(BooleanError::TooFewArguments);
    }
    let fuzz = the_pf.fuzzy_value();
    let mut a_section = BOPAlgoSection::new(the_pf.ds(), fuzz);
    a_section.set_arguments(the_pf.ds().arguments.clone());
    a_section.perform();
    if a_section.has_errors() {
        return Err(BooleanError::InvalidResult("section failed"));
    }
    a_section
        .result_brep()
        .ok_or(BooleanError::InvalidResult("no section result"))
}

/// OCCT shortcut: `BRepAlgoAPI_Splitter(objects, tools).Shape()`.
/// BRep form: returns the compound of the split parts of the OBJECTS.
pub fn splitter(
    objects: &[rcad_kernel::BRep],
    tools: &[rcad_kernel::BRep],
) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    let mut cache = std::collections::HashMap::new();
    for o in objects {
        op.arguments.extend(brep_top_shapes_with_locations(o, &mut global_locs, &mut cache));
    }
    for t in tools {
        op.tools.extend(brep_top_shapes_with_locations(t, &mut global_locs, &mut cache));
    }
    op.locations = global_locs;
    run_build_splitter_brep(&op)
}

/// OCCT shortcut: `BRepAlgoAPI_Fuse(a, b).Shape()`.
pub fn fuse(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    // One shared clone cache across arguments+tools: OCCT keeps cross-argument
    // TShape identity (boptuc_simple ZP3).
    let mut cache = std::collections::HashMap::new();
    op.arguments = brep_top_shapes_with_locations(a, &mut global_locs, &mut cache);
    op.tools = brep_top_shapes_with_locations(b, &mut global_locs, &mut cache);
    op.locations = global_locs;
    run_build_brep(&op, BooleanOpType::Union)
}

/// OCCT shortcut: `BRepAlgoAPI_Common(a, b).Shape()`.
pub fn common(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    let mut cache = std::collections::HashMap::new();
    op.arguments = brep_top_shapes_with_locations(a, &mut global_locs, &mut cache);
    op.tools = brep_top_shapes_with_locations(b, &mut global_locs, &mut cache);
    op.locations = global_locs;
    run_build_brep(&op, BooleanOpType::Intersection)
}

/// OCCT shortcut: `BRepAlgoAPI_Cut(a, b).Shape()`.
pub fn cut(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    let mut cache = std::collections::HashMap::new();
    op.arguments = brep_top_shapes_with_locations(a, &mut global_locs, &mut cache);
    op.tools = brep_top_shapes_with_locations(b, &mut global_locs, &mut cache);
    op.locations = global_locs;
    run_build_brep(&op, BooleanOpType::Cut)
}

/// OCCT shortcut: `BRepAlgoAPI_Section(a, b).Shape()` — the true SECTION
/// operation (BOPAlgo_Section::BuildSection); PCurve options take the
/// BRepAlgoAPI_Section defaults (off).
pub fn section(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    let mut cache = std::collections::HashMap::new();
    op.arguments = brep_top_shapes_with_locations(a, &mut global_locs, &mut cache);
    op.tools = brep_top_shapes_with_locations(b, &mut global_locs, &mut cache);
    op.locations = global_locs;
    run_build_section_brep(&op, false, false, false)
}

/// OCCT shortcut: `BRepAlgoAPI_Cut21(a, b).Shape()` 閳?`b` minus `a`.
pub fn cut21(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    cut(b, a) // swap args 閳?b - a
}

/// Dispatch a boolean operation by [`BooleanOpType`] (legacy convenience API).
/// OCCT: BRepAlgoAPI_BOP::SetOperation(BOPAlgo_Operation) 閳?COMMON/FUSE/CUT/
/// CUT21/SECTION.
pub fn boolean_op(op: BooleanOpType, a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    match op {
        BooleanOpType::Union => fuse(a, b),
        BooleanOpType::Intersection => common(a, b),
        BooleanOpType::Cut => cut(a, b),
        BooleanOpType::Cut21 => cut21(a, b),
        BooleanOpType::Section => section(a, b),
        BooleanOpType::Unknown => Err(BooleanError::TooFewArguments),
    }
}

/// Legacy `boolean_op_with_retry` 閳?the current pipeline already includes the
/// OCCT-style retry ladder, so this is a plain dispatch.
pub fn boolean_op_with_retry(op: BooleanOpType, a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    boolean_op(op, a, b)
}
