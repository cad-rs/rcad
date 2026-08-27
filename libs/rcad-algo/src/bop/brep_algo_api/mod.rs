use rcad_kernel::topo_shape::Shape;
use crate::bop::algo::builder::{Builder, BooleanError, BooleanOpType};
use crate::bop::algo::pave_filler::PaveFiller;
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
}
impl Algo for BuilderAlgo {
    fn is_done(&self) -> bool { self.bs.is_done() }
    fn error(&self) -> Option<&BooleanError> { self.bs.error() }
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
def_bool_op!(SectionOp, Section);

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
    pub fn add_object(&mut self, s: Shape) { self.algo.arguments.push(s); }
    pub fn add_tool(&mut self, s: Shape) { self.algo.arguments.push(s); }
    pub fn build(&mut self) { self.algo.bs.result = self.algo.arguments.first().cloned(); }
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
fn brep_top_shapes_with_locations(
    brep: &rcad_kernel::BRep,
    global_locs: &mut Vec<glam::DAffine3>,
) -> Vec<Shape> {
    let mut map: HashMap<u32, u32> = HashMap::new();
    for (i, loc) in brep.locations.iter().enumerate() {
        let old = (i + 1) as u32; // BRep table index (0 = identity)
        let new = global_locs.len() as u32;
        global_locs.push(*loc);
        map.insert(old, new);
    }
    if map.is_empty() {
        // No located sub-shapes: keep the original TShape graph untouched
        // (deep-copying would remap every TShape pointer and break the
        // vertex_params/pcurve identity keys for no benefit).
        return brep_top_shapes(brep);
    }
    let mut cache: HashMap<u64, Arc<TShape>> = HashMap::new();
    let tops: Vec<Shape> = brep_top_shapes(brep)
        .into_iter()
        .map(|s| remap_location_tree(&s, &map, &mut cache))
        .collect();
    // Second pass: all pointers are final now; rewrite every edge's identity
    // keys against the complete cache (old ptr -> new Arc ptr).
    rewrite_identity_keys(&tops, &cache);
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

/// OCCT shortcut: `BRepAlgoAPI_Fuse(a, b).Shape()`.
pub fn fuse(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    op.arguments = brep_top_shapes_with_locations(a, &mut global_locs);
    op.tools = brep_top_shapes_with_locations(b, &mut global_locs);
    op.locations = global_locs;
    run_build_brep(&op, BooleanOpType::Union)
}

/// OCCT shortcut: `BRepAlgoAPI_Common(a, b).Shape()`.
pub fn common(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    op.arguments = brep_top_shapes_with_locations(a, &mut global_locs);
    op.tools = brep_top_shapes_with_locations(b, &mut global_locs);
    op.locations = global_locs;
    run_build_brep(&op, BooleanOpType::Intersection)
}

/// OCCT shortcut: `BRepAlgoAPI_Cut(a, b).Shape()`.
pub fn cut(a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    let mut op = BuilderAlgo::new();
    let mut global_locs = vec![glam::DAffine3::IDENTITY];
    op.arguments = brep_top_shapes_with_locations(a, &mut global_locs);
    op.tools = brep_top_shapes_with_locations(b, &mut global_locs);
    op.locations = global_locs;
    run_build_brep(&op, BooleanOpType::Cut)
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
        BooleanOpType::Section => cut(a, b),
        BooleanOpType::Unknown => Err(BooleanError::TooFewArguments),
    }
}

/// Legacy `boolean_op_with_retry` 閳?the current pipeline already includes the
/// OCCT-style retry ladder, so this is a plain dispatch.
pub fn boolean_op_with_retry(op: BooleanOpType, a: &rcad_kernel::BRep, b: &rcad_kernel::BRep) -> Result<rcad_kernel::BRep, BooleanError> {
    boolean_op(op, a, b)
}
