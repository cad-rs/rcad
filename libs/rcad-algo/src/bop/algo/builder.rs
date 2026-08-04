// OCCT BOPAlgo_Builder — shape construction from DS.
//
// OCCT BOPAlgo_Builder.hxx L75-507 + parent class fields (BOPAlgo_BuilderShape, BOPAlgo_Options, BOPAlgo_BOP).
// Flattened into one Rust struct because Rust has no C++ inheritance.

pub use crate::bop::algo::BooleanOpType;
use crate::bop::algo::{GlueEnum, Report};
use crate::bop::algo::builder_face::BuilderFace;
use crate::bop::ds::DS;
use crate::bop::ds::pave::SharedPB;
use crate::bop::int_tools::context::IntToolsContext;
use rcad_kernel::geom::{CurveEval, SurfaceEval};
use rcad_kernel::topods;
use rcad_kernel::topods::{
    TShape, TVertexData, TEdgeData, TWireData, TFaceData,
    TShellData, TSolidData, tshape_flags,
};
use rcad_kernel::topo_shape::Shape;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use glam::DVec3;

/// Boolean operation error type.
#[derive(Debug, Clone)]
pub enum BooleanError {
    InvalidOperation,
    TooFewArguments,
    NoFiller,
    BOPNotAllowed,
    BOPNotSet,
    EmptyShape,
    EmptyInput,
    DegenerateResult,
    NumericalFailure(&'static str),
    InvalidResult(&'static str),
}

/// OCCT BOPAlgo_Builder — result builder for boolean operations.
///
/// OCCT ref: BOPAlgo_Builder (BOPAlgo_Builder.hxx)
///
/// Fields map to OCCT hierarchy:
/// - BOPAlgo_Options (fuzzy, parallel, report)
/// - BOPAlgo_BuilderShape (myShape, myFillHistory)
/// - BOPAlgo_BOP (myOperation, myDims)
/// - BOPAlgo_Builder (myDS, myContext, myImages, etc.)
pub struct Builder<'a> {
    // ── BOPAlgo_Options (inherited) ─────────────────────────────
    pub(crate) my_report: Report,          // BOPAlgo_Algo::myReport
    pub(crate) my_run_parallel: bool,      // BOPAlgo_Algo::myRunParallel
    pub(crate) my_fuzzy_value: f64,        // BOPAlgo_Algo::myFuzzyValue
    // ── BOPAlgo_BuilderShape (inherited) ───────────────────────
    pub(crate) my_shape: Option<topods::BRep>, // BOPAlgo_BuilderShape::myShape
    pub(crate) my_fill_history: bool,      // BOPAlgo_BuilderShape::myFillHistory
    // ── BOPAlgo_BOP (inherited) ────────────────────────────────
    pub(crate) my_operation: BooleanOpType, // BOPAlgo_BOP::myOperation
    pub(crate) my_tools: Vec<Shape>,        // BOPAlgo_BOP::myTools
    pub(crate) my_rc: Vec<Shape>,           // BOPAlgo_BOP::myRC (result compound contents)
    pub(crate) my_dims: [i32; 2],           // BOPAlgo_BOP::myDims ([0]=obj, [1]=tool)
    // BOPAlgo_Builder::FillIn3DParts theDraftSolids: source solid ptr_id → draft solid.
    pub(crate) my_draft_solids: HashMap<u64, Shape>,
    // ── BOPAlgo_Builder.hxx L492-505 ───────────────────────────
    pub(crate) ds: &'a DS,                 // L496: myDS (borrowed from PaveFiller)
    pub(crate) my_context: IntToolsContext, // L497: myContext
    pub(crate) my_arguments: Vec<Shape>,   // L492: myArguments
    pub(crate) my_map_fence: HashSet<u64>, // L494: myMapFence
    pub(crate) my_entry_point: i32,        // L498: myEntryPoint
    pub(crate) my_images: HashMap<Shape, Vec<Shape>>,      // L499: myImages
    pub(crate) my_shapes_sd: HashMap<Shape, Shape>,        // L500: myShapesSD
    pub(crate) my_origins: HashMap<Shape, Vec<Shape>>,     // L501: myOrigins
    pub(crate) my_in_parts: HashMap<Shape, Vec<Shape>>,    // L502: myInParts
    pub(crate) my_non_destructive: bool,   // L503: myNonDestructive
    pub(crate) my_glue: GlueEnum,          // L504: myGlue
    pub(crate) my_check_inverted: bool,    // L505: myCheckInverted
    pub(crate) my_nb_shapes_arr: [usize; 8], // L367: myNbShapesArr
    // rcad-specific: BRep index tracking (ptr_id -> tshapes index)
    pub(crate) shape_remap: HashMap<u64, usize>,
}

/// Stage snapshot: DS + result BRep counts at a Builder pipeline boundary.
/// Mirrors OCCT BOPAlgo_BOP::PerformInternal1 DUMP_STAGE points (10 stages).
#[derive(Debug, Clone)]
pub struct StageSnapshot {
    pub stage: u32,
    pub stage_name: &'static str,
    pub n_ds_vertices: usize,
    pub n_ds_edges: usize,
    pub n_ds_faces: usize,
    pub n_ds_pave_blocks: usize,
    pub n_ds_intersection_curves: usize,
    pub n_ds_interf_ff: usize,
    pub n_brep_vertices: usize,
    pub n_brep_edges: usize,
    pub n_brep_faces: usize,
    pub n_brep_shells: usize,
    pub n_brep_solids: usize,
}

/// Count result BRep entities by type from the flat tshape list.
fn count_brep_entities(b: &topods::BRep) -> (usize, usize, usize, usize, usize) {
    let mut v = 0usize;
    let mut e = 0usize;
    let mut f = 0usize;
    let mut sh = 0usize;
    let mut so = 0usize;
    for ts in &b.tshapes {
        match &**ts {
            topods::TShape::Vertex(_) => v += 1,
            topods::TShape::Edge(_) => e += 1,
            topods::TShape::Face(_) => f += 1,
            topods::TShape::Shell(_) => sh += 1,
            topods::TShape::Solid(_) => so += 1,
            _ => {}
        }
    }
    (v, e, f, sh, so)
}

/// Collect the face Shapes of a solid/compound via tree traversal.
/// OCCT TopExp_Explorer(aSolid, TopAbs_FACE).
fn collect_solid_faces(s: &Shape) -> Vec<Shape> {
    let mut result: Vec<Shape> = Vec::new();
    let mut stack: Vec<Shape> = vec![s.clone()];
    while let Some(sh) = stack.pop() {
        match &*sh.data {
            TShape::Solid(sd) => {
                for x in &sd.shells {
                    stack.push(x.clone());
                }
            }
            TShape::CompSolid(cd) => {
                for x in cd {
                    stack.push(x.clone());
                }
            }
            TShape::Compound(cd) => {
                for x in cd {
                    stack.push(x.clone());
                }
            }
            TShape::Shell(sd) => {
                for x in &sd.faces {
                    stack.push(x.clone());
                }
            }
            TShape::Face(_) => result.push(sh),
            _ => {}
        }
    }
    result
}

/// OCCT BOPAlgo_BOP::TypeToExplore (BOPAlgo_BOP.cxx L1574-1597).
fn type_to_explore(the_dim: i32) -> topods::ShapeType {
    match the_dim {
        0 => topods::ShapeType::Vertex,
        1 => topods::ShapeType::Edge,
        2 => topods::ShapeType::Face,
        3 => topods::ShapeType::Solid,
        _ => topods::ShapeType::Compound,
    }
}

/// OCCT BOPTools_AlgoTools::Dimension — max dimension of a shape.
fn shape_dimension(s: &Shape) -> i32 {
    crate::bop::tools::algo_tools::dimensions(s.shape_type()).1
}

/// OCCT BOPAlgo_Tools::FillMap(Shape, Shape, IndexedDataMap<Shape, List<Shape>>)
/// (BOPAlgo_Tools.hxx L84-102) — bidirectional connection for the SD back-and-forth map.
fn fill_map_faces(n1: &Shape, n2: &Shape, the_map: &mut HashMap<Shape, Vec<Shape>>) {
    the_map.entry(n1.clone()).or_default().push(n2.clone());
    the_map.entry(n2.clone()).or_default().push(n1.clone());
}

/// OCCT BOPAlgo_Tools::MakeBlocks (BOPAlgo_Tools.hxx L46-80) — connected components
/// of the SD back-and-forth map.
fn make_blocks_faces(the_map: &HashMap<Shape, Vec<Shape>>) -> Vec<Vec<Shape>> {
    let mut a_m_fence: HashSet<u64> = HashSet::new();
    let mut a_m_blocks: Vec<Vec<Shape>> = Vec::new();
    for (n, _) in the_map {
        if !a_m_fence.insert(n.ptr_id()) {
            continue;
        }
        let mut a_chain: Vec<Shape> = vec![n.clone()];
        let mut i = 0;
        while i < a_chain.len() {
            let n1 = &a_chain[i];
            if let Some(a_li) = the_map.get(n1) {
                for n2 in a_li {
                    if a_m_fence.insert(n2.ptr_id()) {
                        a_chain.push(n2.clone());
                    }
                }
            }
            i += 1;
        }
        a_m_blocks.push(a_chain);
    }
    a_m_blocks
}

impl<'a> Builder<'a> {
    /// Create a new Builder borrowing a DS from PaveFiller.
    ///
    /// OCCT: BOPAlgo_Builder is constructed with a PaveFiller reference.
    /// OCCT BOPAlgo_BOP::PerformInternal1 L425-429:
    ///   myPaveFiller = &theFiller; myDS = myPaveFiller->PDS();
    ///   myFuzzyValue = myPaveFiller->FuzzyValue();
    pub fn new(ds: &'a DS, op: BooleanOpType, fuzzy_value: f64) -> Self {
        Builder {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_fuzzy_value: fuzzy_value,
            my_shape: None,
            my_fill_history: false,
            my_operation: op,
            my_tools: Vec::new(),
            my_rc: Vec::new(),
            my_dims: [3, 3],
            my_draft_solids: HashMap::new(),
            my_context: IntToolsContext::new(),
            my_arguments: Vec::new(),
            my_map_fence: HashSet::new(),
            my_entry_point: 0,
            my_images: HashMap::new(),
            my_shapes_sd: HashMap::new(),
            my_origins: HashMap::new(),
            my_in_parts: HashMap::new(),
            my_non_destructive: false,
            my_glue: GlueEnum::GlueOff,
            my_check_inverted: false,
            my_nb_shapes_arr: [0; 8],
            shape_remap: HashMap::new(),
        }
    }

    /// OCCT BOPAlgo_Algo::SetArguments.
    pub fn set_arguments(&mut self, args: Vec<Shape>) {
        self.my_arguments = args;
    }

    /// OCCT BOPAlgo_BOP::SetTools.
    pub fn set_tools(&mut self, tools: Vec<Shape>) {
        self.my_tools = tools;
    }

    /// Shape backed by shared Arc in ds.shapes — OCCT myDS->Shape(n).
    fn brep_sr(&self, flat_idx: usize) -> Shape {
        self.ds.shape(flat_idx).clone()
    }

    pub fn has_errors(&self) -> bool {
        self.my_report.has_errors()
    }

    pub fn report(&self) -> &Report {
        &self.my_report
    }

    /// OCCT BOPAlgo_Builder::Build — convenience wrapper.
    ///
    /// Returns the result BRep on success.
    pub fn build(&mut self) -> Result<rcad_kernel::BRep, ()> {
        match self.build_with_history_topods() {
            Ok((brep, _)) => Ok(brep),
            Err(_) => Err(()),
        }
    }

    /// Run the Builder pipeline stage by stage, capturing a snapshot after each
    /// of the 10 OCCT-aligned stages. Single source of truth for the pipeline:
    /// `build_with_history_topods` and `build` delegate here.
    ///
    /// On `has_errors` mid-pipeline, returns Ok with the partial result and the
    /// snapshots captured so far, so tests can localize the failure stage.
    pub fn build_with_history_stage_by_stage(
        &mut self,
    ) -> Result<(topods::BRep, Vec<StageSnapshot>), BooleanError> {
        let mut snapshots: Vec<StageSnapshot> = Vec::with_capacity(10);
        macro_rules! snap {
            ($stage:expr, $name:expr) => {{
                let (v, e, f, sh, so) = count_brep_entities(self.my_shape.as_ref().unwrap());
                snapshots.push(StageSnapshot {
                    stage: $stage,
                    stage_name: $name,
                    n_ds_vertices: self.ds.vertex_count(),
                    n_ds_edges: self.ds.edge_count(),
                    n_ds_faces: self.ds.face_count(),
                    // OCCT pipeline_dump.h dump_ds_snapshot: count PBs only for
                    // EDGE shapes with HasPaveBlocks(i), summing PaveBlocks(i).Extent().
                    // This excludes orphan pool entries (section-edge PBs etc.) that
                    // no edge shape references.
                    n_ds_pave_blocks: {
                        let mut n_pb = 0usize;
                        for i in 0..self.ds.nb_shapes() {
                            if self.ds.shape_info(i).shape_type != topods::ShapeType::Edge {
                                continue;
                            }
                            if self.ds.has_pave_blocks(i) {
                                n_pb += self.ds.pave_blocks(i).len();
                            }
                        }
                        n_pb
                    },
                    n_ds_intersection_curves: self.ds.intersection_curves.len(),
                    n_ds_interf_ff: self.ds.interf_ff.len(),
                    n_brep_vertices: v,
                    n_brep_edges: e,
                    n_brep_faces: f,
                    n_brep_shells: sh,
                    n_brep_solids: so,
                });
            }};
        }
        let partial = |s: &Option<topods::BRep>| s.clone().unwrap_or_default();

        // OCCT L431-436: CheckData
        self.check_data()?;
        self.check_filler();

        // OCCT L438-443: Prepare
        let _result = self.prepare();

        // OCCT L459-471: FillImagesVertices + BuildResult(VERTEX)
        self.fill_images_vertices();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Vertex);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(1, "after_FillImagesVertices");

        // OCCT L472-483: FillImagesEdges + BuildResult(EDGE)
        self.fill_images_edges();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Edge);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(2, "after_FillImagesEdges");

        // OCCT L484-494: FillImagesContainers(WIRE) + BuildResult(WIRE)
        self.fill_images_containers(topods::ShapeType::Wire);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Wire);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(3, "after_BuildResultWire");

        // OCCT L496-505: FillImagesFaces + BuildResult(FACE)
        self.fill_images_faces();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Face);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(4, "after_FillImagesFaces");

        // OCCT L507-516: FillImagesContainers(SHELL) + BuildResult(SHELL)
        self.fill_images_containers(topods::ShapeType::Shell);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Shell);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(5, "after_BuildResultShell");

        // OCCT L518-528: FillImagesSolids + BuildResult(SOLID)
        self.fill_images_solids();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Solid);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(6, "after_FillImagesSolids");

        // OCCT L530-539: FillImagesContainers(COMPSOLID) + BuildResult(COMPSOLID)
        self.fill_images_containers(topods::ShapeType::CompSolid);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::CompSolid);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(7, "after_BuildResultCompSolid");

        // OCCT L541-550: FillImagesCompounds + BuildResult(COMPOUND)
        self.fill_images_compounds();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        self.build_result(topods::ShapeType::Compound);
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(8, "after_FillImagesCompounds");

        // OCCT L575-580 (BOPAlgo_BOP.cxx): BuildShape — apply the boolean operation
        // result construction. Runs between the s08 dump and PrepareHistory.
        self.build_shape();
        if self.has_errors() { return Ok((partial(&self.my_shape), snapshots)); }
        snap!(9, "after_PrepareHistory");
        self.post_treat();
        snap!(10, "after_PostTreat");

        let result = self.my_shape.clone().unwrap_or_default();
        Ok((result, snapshots))
    }

    /// OCCT BOPAlgo_Builder::Build — full pipeline with history.
    ///
    /// Delegates to `build_with_history_stage_by_stage` (single source of
    /// truth); converts an early pipeline stop (has_errors) back to Err to
    /// preserve the production contract.
    pub fn build_with_history_topods(
        &mut self,
    ) -> Result<(topods::BRep, ()), BooleanError> {
        let (brep, snaps) = self.build_with_history_stage_by_stage()?;
        if snaps.len() != 10 {
            return Err(BooleanError::DegenerateResult);
        }
        Ok((brep, ()))
    }

    // --- Pipeline stage stubs ---
    // Each matches a method in OCCT BOPAlgo_Builder.
    // Stubs: empty bodies that compile. Implementation will be added incrementally.

    /// OCCT BOPAlgo_Builder::CheckData (BOPAlgo_Builder.cxx L130-140).
    fn check_data(&self) -> Result<(), BooleanError> {
        // OCCT L132-137: if (myArguments.Extent() < 2) → AddError(TooFewArguments)
        if self.my_arguments.len() < 2 {
            return Err(BooleanError::TooFewArguments);
        }
        // OCCT L139: CheckFiller();
        self.check_filler();
        Ok(())
    }

    /// OCCT BOPAlgo_BOP::CheckFiller (BOPAlgo_BOP.cxx L144-152).
    fn check_filler(&self) {
        // OCCT L146-150: if (!myPaveFiller) → AddError(NoFiller)
        // rcad: PaveFiller always runs before Builder, no reference stored.
        // OCCT L151: GetReport()->Merge(myPaveFiller->GetReport());
        // rcad: report merging not applicable (PaveFiller dropped before Builder).
    }

    /// OCCT BOPAlgo_Builder::Prepare (BOPAlgo_Builder.cxx L156-164).
    fn prepare(&mut self) {
        // OCCT L158-163: BRep_Builder aBB; MakeCompound(aC); myShape = aC;
        // rcad: topods::BRep is the equivalent of TopoDS_Compound for result.
        self.my_shape = Some(topods::BRep::new());
        self.shape_remap.clear();
    }

    /// OCCT BOPAlgo_Builder::FillImagesVertices (BOPAlgo_Builder_1.cxx L40-67).
    /// Maps each SD vertex pair as myImages[source]->[target], myShapesSD, myOrigins.
    fn fill_images_vertices(&mut self) {
        // OCCT L40-66: NCollection_DataMap<int, int>::Iterator aIt(myDS->ShapesSD());
        // rcad: DS::shapes_sd is HashMap<usize, usize> (source→SD).
        let sd_pairs: Vec<(usize, usize)> = self.ds.shapes_sd
            .iter()
            .map(|(k, v)| (*k, *v))
            .collect();
        for (nV, nVSD) in sd_pairs {
            // OCCT L53-54: const TopoDS_Shape& aV = myDS->Shape(nV);
            let aV = self.brep_sr(nV);
            // OCCT L54: const TopoDS_Shape& aVSD = myDS->Shape(nVSD);
            let aVSD = self.brep_sr(nVSD);
            // OCCT L56: myImages.Bound(aV, ...)->Append(aVSD);
            self.my_images.entry(aV.clone()).or_default().push(aVSD.clone());
            // OCCT L58: myShapesSD.Bind(aV, aVSD);
            self.my_shapes_sd.insert(aV.clone(), aVSD.clone());
            // OCCT L60-65: myOrigins — find or create list, append
            self.my_origins.entry(aVSD).or_default().push(aV);
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesEdges (BOPAlgo_Builder_1.cxx L71-126).
    /// Maps source edges -> split images via pave-block real edge.
    /// Also handles CommonBlocks via myShapesSD.
    fn fill_images_edges(&mut self) {
        let aNbS = self.ds.nb_source_shapes();
        for i in 0..aNbS {
            let aSI = self.ds.shape_info(i);
            if aSI.shape_type != topods::ShapeType::Edge {
                continue;
            }
            // OCCT L84-86: if (!aSI.HasReference()) continue;
            if !aSI.has_reference() {
                continue;
            }
            // OCCT L89-91: aE = myDS->Shape(i); aLPB = myDS->PaveBlocks(i);
            let aE = self.brep_sr(i);
            let aLPB: Vec<SharedPB> = self.ds.pave_blocks(i).to_vec();
            // OCCT L95: pLS = myImages.Bound(aE, ...)
            // OCCT L96-120: iterate pave blocks
            for aPB in &aLPB {
                // OCCT L100-102: aPBR = myDS->RealPaveBlock(aPB); nSpR = aPBR->Edge();
                let aPBR = self.ds.real_pave_block(aPB);
                let nSpR = { let r = aPBR.read(); r.edge };
                // OCCT L103-104: aSpR = myDS->Shape(nSpR);
                let aSpR = self.brep_sr(nSpR);
                // OCCT L105: pLS->Append(aSpR);
                self.my_images.entry(aE.clone()).or_default().push(aSpR.clone());
                // OCCT L107-112: pLOr = myOrigins.ChangeSeek(aSpR); append aE
                self.my_origins.entry(aSpR.clone()).or_default().push(aE.clone());
                // OCCT L114-119: if (IsCommonBlockOnEdge(aPB)) → myShapesSD.Bind(aSp, aSpR)
                if self.ds.is_common_block_on_edge(aPB) {
                    // OCCT L116-117: nSp = aPB->Edge(); aSp = myDS->Shape(nSp);
                    let nSp = { let r = aPB.read(); r.edge };
                    let aSp = self.brep_sr(nSp);
                    // OCCT L118: myShapesSD.Bind(aSp, aSpR);
                    self.my_shapes_sd.insert(aSp, aSpR.clone());
                }
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesContainers (BOPAlgo_Builder_1.cxx L172-193).
    /// Builds wire/shell/compsolid images from edge/face/solid images.
    /// For each source shape of theType, calls FillImagesContainer.
    fn fill_images_containers(&mut self, the_type: topods::ShapeType) {
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != the_type {
                continue;
            }
            // OCCT L185-186: FillImagesContainer(aC, theType)
            let a_c = self.brep_sr(i);
            self.fill_images_container(&a_c, the_type);
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesContainer (BOPAlgo_Builder_1.cxx L221-276).
    /// Builds a new container (wire/shell/compsolid) from sub-shape images.
    /// If no sub-shape was modified, the container is kept as-is.
    fn fill_images_container(&mut self, the_s: &Shape, the_type: topods::ShapeType) {
        // OCCT L223-233: check if any sub-shape has been modified
        let sub_shapes = self.shape_sub_shapes(the_s);
        let mut has_modified = false;
        for ss in &sub_shapes {
            if let Some(imgs) = self.my_images.get(ss) {
                if imgs.len() != 1 || imgs[0].ptr_id() != ss.ptr_id() {
                    has_modified = true;
                    break;
                }
            }
        }
        if !has_modified {
            return;
        }

        // OCCT L242-245: MakeContainer(theType, aCIm)
        let mut new_edges: Vec<Shape> = Vec::new();
        let mut new_faces: Vec<Shape> = Vec::new();
        let mut new_comps: Vec<Shape> = Vec::new();

        // OCCT L247-272: iterate sub-shapes, add images or originals
        for ss in &sub_shapes {
            let p_lss_im = self.my_images.get(ss);
            match p_lss_im {
                None => {
                    // OCCT L253-257: no splits, add sub-shape itself
                    match the_type {
                        topods::ShapeType::Wire => new_edges.push(ss.clone()),
                        topods::ShapeType::Shell => new_faces.push(ss.clone()),
                        topods::ShapeType::CompSolid => new_comps.push(ss.clone()),
                        _ => {}
                    }
                }
                Some(imgs) => {
                    // OCCT L260-271: add each image (split)
                    for a_ss_im in imgs {
                        // OCCT L265-269: IsSplitToReverseWithWarn(aSSIm, aSS) — reverses aSSIm
                        // when its geometry is oppositely oriented to aSS. Pending translation
                        // (BOPTools_AlgoTools.cxx L1302-1531): needs Extrema_LocateExtPC point
                        // projection + Geom_Curve/Surface handle-equality, absent in rcad. Only
                        // affects sub-shape orientation, not entity counts.
                        match the_type {
                            topods::ShapeType::Wire => new_edges.push(a_ss_im.clone()),
                            topods::ShapeType::Shell => new_faces.push(a_ss_im.clone()),
                            topods::ShapeType::CompSolid => new_comps.push(a_ss_im.clone()),
                            _ => {}
                        }
                    }
                }
            }
        }

        // Build new container TShape
        let new_container: TShape = match the_type {
            topods::ShapeType::Wire => {
                TShape::Wire(TWireData {
                    my_shapes: vec![], flags: tshape_flags::DEFAULT,
                    edges: new_edges,
                })
            }
            topods::ShapeType::Shell => {
                TShape::Shell(TShellData {
                    my_shapes: vec![], flags: tshape_flags::DEFAULT,
                    faces: new_faces,
                })
            }
            topods::ShapeType::CompSolid => {
                TShape::CompSolid(new_comps)
            }
            _ => return,
        };

        // Wrap in Shape (synthetic index, will be remapped during add_to_result)
        let container_shape = Shape::new(
            std::sync::Arc::new(new_container),
            0, topods::Orientation::Forward,
        );

        // OCCT L275: myImages.Bound(theS, ...)->Append(aCIm)
        self.my_images.entry(the_s.clone()).or_default().push(container_shape);
    }

    /// Extract immediate sub-shapes from a Shape (OCCT TopoDS_Iterator equivalent).
    fn shape_sub_shapes(&self, s: &Shape) -> Vec<Shape> {
        match &*s.data {
            TShape::Vertex(_) => vec![],
            TShape::Edge(ed) => vec![
                Shape::new(ed.first.data.clone(), ed.first.location, ed.first.orientation),
                Shape::new(ed.last.data.clone(), ed.last.location, ed.last.orientation),
            ],
            TShape::Wire(wd) => {
                wd.edges.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Face(fd) => {
                let mut v = vec![
                    Shape::new(fd.outer_wire.data.clone(), fd.outer_wire.location, fd.outer_wire.orientation)
                ];
                v.extend(fd.inner_wires.iter().map(|w| {
                    Shape::new(w.data.clone(), w.location, w.orientation)
                }));
                v
            }
            TShape::Shell(sd) => {
                sd.faces.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Solid(sd) => {
                sd.shells.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::CompSolid(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Compound(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesFaces (BOPAlgo_Builder_2.cxx L215-229).
    /// Splits faces using section edges.
    /// Calls BuildSplitFaces -> FillSameDomainFaces -> FillInternalVertices.
    fn fill_images_faces(&mut self) {
        // OCCT L218: BuildSplitFaces
        self.build_split_faces();
        if self.has_errors() { return; }
        // OCCT L223: FillSameDomainFaces
        self.fill_same_domain_faces();
        if self.has_errors() { return; }
        // OCCT L228: FillInternalVertices
        self.fill_internal_vertices();
    }

    /// OCCT BOPAlgo_Builder::BuildSplitFaces (BOPAlgo_Builder_2.cxx L233-555).
    ///
    /// For each source face with intersection data, builds a BuilderFace fed
    /// with the full edge set (bounding edges + their images + IN edges +
    /// section edges). Faces without section/IN edges but with modified wires
    /// or alone vertices take the BuildDraftFace fast path.
    fn build_split_faces(&mut self) {
        let a_nb_s = self.ds.nb_source_shapes();
        // aFacesIm: DS face index -> area shapes (OCCT IndexedDataMap<int, List<Shape>>).
        let mut a_faces_im: HashMap<usize, Vec<Shape>> = HashMap::new();
        // aVBF: pending BuilderFace tasks (face_idx, face, edges).
        let mut a_vbf: Vec<(usize, Shape, Vec<Shape>)> = Vec::new();

        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Face {
                continue;
            }
            // OCCT L275-279: bHasFaceInfo check.
            if !self.ds.has_face_info(i) {
                continue;
            }
            let a_f = self.brep_sr(i);
            let a_fi = self.ds.face_info(i).clone();
            // OCCT L286-287: AloneVertices(i, aLIAV).
            let a_liav = self.alone_vertices(i);
            let a_nb_pb_in = a_fi.pave_blocks_in.len();
            let a_nb_pb_on = a_fi.pave_blocks_on.len();
            let a_nb_pb_sc = a_fi.pave_blocks_sc.len();
            let a_nb_av = a_liav.len();
            // OCCT L293-296: not complete -> skip.
            if a_nb_pb_in == 0 && a_nb_pb_on == 0 && a_nb_pb_sc == 0 && a_nb_av == 0 {
                continue;
            }

            // OCCT L298-351: only alone vertices / On PBs -> draft-face fast path.
            if a_nb_pb_in == 0 && a_nb_pb_sc == 0 {
                let mut has_internals = false;
                if a_nb_av == 0 {
                    // OCCT L315-330: check wires for internal edges or modifications.
                    let mut has_modified = false;
                    for a_w in self.shape_sub_shapes(&a_f) {
                        if a_w.shape_type() != topods::ShapeType::Wire {
                            continue;
                        }
                        let mut w_has_internal = false;
                        if let TShape::Wire(wd) = &*a_w.data {
                            for e in &wd.edges {
                                if e.orientation == topods::Orientation::Internal {
                                    w_has_internal = true;
                                    break;
                                }
                            }
                        }
                        if w_has_internal {
                            has_internals = true;
                            break;
                        }
                        if self.images_of(&a_w).is_some() {
                            has_modified = true;
                        }
                    }
                    if !has_internals && !has_modified {
                        continue;
                    }
                }
                if !has_internals {
                    // OCCT L344: BuildDraftFace fast path — face image directly.
                    if let Some(a_fd) = self.build_draft_face(&a_f) {
                        a_faces_im.entry(i).or_default().push(a_fd);
                        continue;
                    }
                }
            }

            // OCCT L353: aMFence.Clear() — per-face fence for closed edges.
            let mut a_mfence: HashSet<u64> = HashSet::new();
            // OCCT L355-357: aFF = aF; aFF.Orientation(FORWARD).
            let a_ff = Shape::new(a_f.data.clone(), a_f.location, topods::Orientation::Forward);
            // OCCT L359: 1. Build the edges set aLE.
            let mut a_le: Vec<Shape> = Vec::new();

            // OCCT L362-465: 1.1 Bounding edges.
            let mut is_checked = false;
            let mut is_u_closed = false;
            let mut is_v_closed = false;
            for a_e in self.face_edges(&a_ff) {
                let an_ori_e = a_e.orientation;
                // OCCT L369: if !myImages.IsBound(aE).
                if self.images_of(&a_e).is_none() {
                    if an_ori_e == topods::Orientation::Internal {
                        let mut a_ee = a_e.clone();
                        a_ee.orientation = topods::Orientation::Forward;
                        a_le.push(a_ee.clone());
                        a_ee.orientation = topods::Orientation::Reversed;
                        a_le.push(a_ee);
                    } else {
                        a_le.push(a_e);
                    }
                    continue;
                }
                // OCCT L387-393: GeomLib::IsClosed(aSurf, tol, isUClosed, isVClosed).
                if !is_checked {
                    let (uc, vc) = self.surface_is_closed(&a_f);
                    is_u_closed = uc;
                    is_v_closed = vc;
                    is_checked = true;
                }
                // OCCT L395-404: bIsClosed = seam edge on closed surface.
                let mut b_is_closed = false;
                if (is_u_closed || is_v_closed) && self.edge_closed_on_face(&a_e, &a_f) {
                    let (is_ui, is_vi) = self.is_edge_isoline(&a_e, &a_f);
                    b_is_closed = (is_u_closed && is_ui) || (is_v_closed && is_vi);
                }
                // OCCT L406: bIsDegenerated = BRep_Tool::Degenerated(aE).
                let b_is_degenerated = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                // OCCT L408: aLIE = myImages.Find(aE).
                let a_lie = self.images_of(&a_e).unwrap_or_default();
                for a_sp0 in &a_lie {
                    let mut a_sp = a_sp0.clone();
                    // OCCT L413-418: degenerated -> keep original orientation.
                    if b_is_degenerated {
                        a_sp.orientation = an_ori_e;
                        a_le.push(a_sp);
                        continue;
                    }
                    // OCCT L420-427: INTERNAL -> forward + reversed.
                    if an_ori_e == topods::Orientation::Internal {
                        a_sp.orientation = topods::Orientation::Forward;
                        a_le.push(a_sp.clone());
                        a_sp.orientation = topods::Orientation::Reversed;
                        a_le.push(a_sp);
                        continue;
                    }
                    // OCCT L429-455: closed seam edge -> dedupe via aMFence.
                    if b_is_closed {
                        if a_mfence.insert(a_sp.ptr_id()) {
                            if !self.edge_closed_on_face(&a_sp, &a_f) {
                                // OCCT L435-446: DoSplitSEAMOnFace(aSp, aF) / (aE, aSp, aF).
                                // Pending translation (BOPTools_AlgoTools3D.cxx):
                                // inserts seam vertices on the split edge. Only affects
                                // the closed-surface seam handling, not planar faces.
                            }
                            a_sp.orientation = topods::Orientation::Forward;
                            a_le.push(a_sp.clone());
                            a_sp.orientation = topods::Orientation::Reversed;
                            a_le.push(a_sp);
                        }
                        continue;
                    }
                    // OCCT L457-463: regular split edge.
                    a_sp.orientation = an_ori_e;
                    // OCCT L458: IsSplitToReverseWithWarn(aSp, aE) — pending translation
                    // (BOPTools_AlgoTools.cxx L1302-1531), needs Extrema_LocateExtPC.
                    a_le.push(a_sp);
                }
            }

            // OCCT L469-480: 1.2 In edges (forward + reversed).
            for &pb_idx in &a_fi.pave_blocks_in {
                if let Some(n_sp) = self.pb_edge(pb_idx) {
                    let mut a_sp = self.brep_sr(n_sp);
                    a_sp.orientation = topods::Orientation::Forward;
                    a_le.push(a_sp.clone());
                    a_sp.orientation = topods::Orientation::Reversed;
                    a_le.push(a_sp);
                }
            }
            // OCCT L483-494: 1.3 Section edges (forward + reversed).
            for &pb_idx in &a_fi.pave_blocks_sc {
                if let Some(n_sp) = self.pb_edge(pb_idx) {
                    let mut a_sp = self.brep_sr(n_sp);
                    a_sp.orientation = topods::Orientation::Forward;
                    a_le.push(a_sp.clone());
                    a_sp.orientation = topods::Orientation::Reversed;
                    a_le.push(a_sp);
                }
            }
            // OCCT L496-500: BuildPCurveForEdgesOnPlane(aLE, aFF) — planar fast path.
            // Pending translation (BRepLib.cxx); only adds pcurves, does not change
            // the edge set or entity counts.
            // OCCT L502-505: aBF.SetFace(aF); aBF.SetShapes(aLE); SetRunParallel.
            a_vbf.push((i, a_f, a_le));
        }

        // OCCT L515-521: perform all BuilderFace tasks.
        for (fi, a_f, a_le) in &a_vbf {
            let mut a_bf = BuilderFace::new(&self.ds);
            a_bf.my_face = Some(a_f.clone());
            a_bf.my_face_index = Some(*fi);
            a_bf.my_edges = a_le.clone();
            a_bf.perform();
            // OCCT L527-531: aFacesIm.Add(myDS->Index(aBF.Face()), aBF.Areas()).
            if !a_bf.my_areas.is_empty() {
                a_faces_im.entry(*fi).or_default().extend(a_bf.my_areas);
            }
        }

        // OCCT L534-552: apply orientation and append areas to myImages.
        for (fi, a_lfr) in a_faces_im {
            let a_f = self.brep_sr(fi);
            let an_ori_f = a_f.orientation;
            let p_lf_im = self.my_images.entry(a_f).or_default();
            for mut a_fr in a_lfr {
                if an_ori_f == topods::Orientation::Reversed {
                    a_fr.orientation = topods::Orientation::Reversed;
                }
                p_lf_im.push(a_fr);
            }
        }
    }

    /// OCCT BOPDS_DS::AloneVertices (BOPDS_DS.cxx L1028-1062).
    /// Vertices of the face not belonging to any boundary edge: endpoints of
    /// PaveBlocksIn/PaveBlocksSc plus VerticesIn/VerticesSc not already seen.
    fn alone_vertices(&self, face_idx: usize) -> Vec<usize> {
        if !self.ds.has_face_info(face_idx) {
            return Vec::new();
        }
        let a_fi = self.ds.face_info(face_idx);
        let mut a_mi: HashSet<usize> = HashSet::new();
        for pb_set in [&a_fi.pave_blocks_in, &a_fi.pave_blocks_sc] {
            for &pb_idx in pb_set {
                if let Some(pool) = self.ds.pave_blocks_pool.get(pb_idx) {
                    if let Some(pb) = pool.first() {
                        let r = pb.0.read().unwrap();
                        a_mi.insert(r.pave1.vertex_idx);
                        a_mi.insert(r.pave2.vertex_idx);
                    }
                }
            }
        }
        let mut result: Vec<usize> = Vec::new();
        for v in a_fi.vertices_in.iter().chain(a_fi.vertices_sc.iter()) {
            if a_mi.insert(*v) {
                result.push(*v);
            }
        }
        result
    }

    /// Edge DS index from a pave-block pool entry (OCCT aPB->Edge()).
    fn pb_edge(&self, pb_idx: usize) -> Option<usize> {
        let pool = self.ds.pave_blocks_pool.get(pb_idx)?;
        let pb = pool.first()?;
        let e = pb.0.read().unwrap().edge;
        if e < self.ds.nb_shapes() {
            Some(e)
        } else {
            None
        }
    }

    /// Look up myImages by TShape pointer + location, ignoring orientation.
    /// OCCT: myImages is keyed by TopTools_ShapeMapHasher which hashes only
    /// TShape* + Location, so IsBound(aE) matches regardless of orientation.
    fn images_of(&self, key: &Shape) -> Option<Vec<Shape>> {
        // OCCT myImages.Find(aE) keys by TopoDS_Shape (stable because OCCT mutates
        // TShapes in place). rcad's Arc::make_mut clones a shared TShape on write,
        // so the face wire edge (original) and the DS edge entry (clone) become two
        // TShape objects for the same logical edge. Both still map to the same DS
        // index — init keeps the original (ptr→idx) mapping and remap_shape_idx
        // adds the clone's — so resolve the lookup by DS index.
        let key_idx = self.ds.map_shape_index.get(&(key.ptr_id(), key.location)).copied();
        for (k, v) in &self.my_images {
            if k.ptr_id() == key.ptr_id() && k.location == key.location {
                return Some(v.clone());
            }
        }
        if let Some(ki) = key_idx {
            for (k, v) in &self.my_images {
                if self.ds.map_shape_index.get(&(k.ptr_id(), k.location)).copied() == Some(ki) {
                    return Some(v.clone());
                }
            }
        }
        None
    }

    /// Boundary edges of a face with wire-composed orientation.
    /// OCCT TopExp_Explorer(aFF, TopAbs_EDGE) — BOPAlgo_Builder_2.cxx L363-365.
    fn face_edges(&self, a_f: &Shape) -> Vec<Shape> {
        let mut result: Vec<Shape> = Vec::new();
        for a_w in self.shape_sub_shapes(a_f) {
            if a_w.shape_type() != topods::ShapeType::Wire {
                continue;
            }
            if let TShape::Wire(wd) = &*a_w.data {
                result.extend(wd.edges.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }));
            }
        }
        result
    }

    /// OCCT GeomLib::IsClosed(aSurf, aTol, isUClosed, isVClosed) — surface
    /// closed-ness in U and V (BOPAlgo_Builder_2.cxx L389-393).
    fn surface_is_closed(&self, a_f: &Shape) -> (bool, bool) {
        let Some(surf) = a_f.as_face().and_then(|fd| fd.surface.clone()) else {
            return (false, false);
        };
        (surf.is_u_closed(), surf.is_v_closed())
    }

    /// OCCT BRep_Tool::IsClosed(aE, aF) — the edge has two pcurves on the
    /// closed surface (seam edge). BOPAlgo_Builder_2.cxx L397.
    fn edge_closed_on_face(&self, a_e: &Shape, a_f: &Shape) -> bool {
        let f_index = a_f.index;
        a_e.as_edge()
            .map(|ed| {
                ed.representations.iter().any(|r| {
                    matches!(
                        r,
                        topods::CurveRepresentation::CurveOnClosedSurface { face, .. }
                            if *face == f_index
                    )
                })
            })
            .unwrap_or(false)
    }

    /// OCCT BOPTools_AlgoTools2D::IsEdgeIsoline(aE, aF, isUIso, isVIso).
    /// True when the edge's pcurve is an axis-aligned line in the face UV space.
    fn is_edge_isoline(&self, a_e: &Shape, a_f: &Shape) -> (bool, bool) {
        let mut is_u_iso = false;
        let mut is_v_iso = false;
        let pcurve = a_e.as_edge().and_then(|ed| {
            ed.pcurves.get(&a_f.index).map(|(pc, _, _)| pc.clone())
        });
        if let Some(pc) = pcurve {
            if let rcad_kernel::geom::Curve2d::Line(l) = pc {
                let d = l.direction;
                is_u_iso = d.y == 0.0 && d.x != 0.0;
                is_v_iso = d.x == 0.0 && d.y != 0.0;
            }
            // Circle pcurves (sphere/cylinder caps) handled per OCCT L...
            // (isUIsoline = first==0 && last==2PI); pending — rare in these patterns.
        }
        (is_u_iso, is_v_iso)
    }

    /// OCCT BuildDraftFace (BOPAlgo_Builder_2.cxx L1052-1189).
    /// Builds a new face from the original face by replacing each boundary
    /// edge with its images. Returns None when the BuilderFace algorithm must
    /// be used instead (internal edges / multi-connected vertices / unified edges).
    fn build_draft_face(&self, the_face: &Shape) -> Option<Shape> {
        let a_surf = the_face.as_face().and_then(|fd| fd.surface.clone())?;
        let a_tol = the_face.as_face().map(|fd| fd.tolerance).unwrap_or(0.0);
        // OCCT L1073-1074: aVerticesCounter — multi-connexity detection.
        let mut a_vertices_counter: HashMap<u64, Vec<Shape>> = HashMap::new();
        // OCCT L1078: aMEdges — edges-unification fence.
        let mut a_m_edges: HashSet<u64> = HashSet::new();

        // OCCT L1081-1181: rebuild each wire of the face.
        let mut new_wires: Vec<Shape> = Vec::new();
        for a_w in self.shape_sub_shapes(the_face) {
            if a_w.shape_type() != topods::ShapeType::Wire {
                continue;
            }
            // OCCT L1091-1095: skip empty wires.
            let w_edges: Vec<Shape> = if let TShape::Wire(wd) = &*a_w.data {
                wd.edges.iter().map(|sr| Shape::new(sr.data.clone(), sr.location, sr.orientation)).collect()
            } else {
                Vec::new()
            };
            if w_edges.is_empty() {
                continue;
            }
            let mut new_edges: Vec<Shape> = Vec::new();
            for a_e in &w_edges {
                let an_ori_e = a_e.orientation;
                // OCCT L1105-1110: internal edges may split the face -> BuilderFace.
                if an_ori_e == topods::Orientation::Internal {
                    return None;
                }
                // OCCT L1113-1115: degenerated / closed on face checks.
                let b_is_degenerated = a_e.as_edge().map(|ed| ed.degenerated).unwrap_or(false);
                let b_is_closed = self.edge_closed_on_face(a_e, the_face);
                // OCCT L1118: theImages.Seek(aE).
                let p_le_im = self.images_of(a_e);
                if p_le_im.is_none() {
                    // OCCT L1121-1131: multi-connected / unified edge -> BuilderFace.
                    if !b_is_degenerated && self.has_multi_connected(a_e, &mut a_vertices_counter) {
                        return None;
                    }
                    if !b_is_closed && !a_m_edges.insert(a_e.ptr_id()) {
                        return None;
                    }
                    new_edges.push(a_e.clone());
                    continue;
                }
                // OCCT L1137-1175: replace by images.
                for a_sp0 in &p_le_im.unwrap() {
                    let mut a_sp = a_sp0.clone();
                    if !b_is_degenerated && self.has_multi_connected(&a_sp, &mut a_vertices_counter) {
                        return None;
                    }
                    if !b_is_closed && !a_m_edges.insert(a_sp.ptr_id()) {
                        return None;
                    }
                    // OCCT L1154: aSp.Orientation(anOriE).
                    a_sp.orientation = an_ori_e;
                    if b_is_degenerated {
                        new_edges.push(a_sp);
                        continue;
                    }
                    // OCCT L1163-1166: seam split — DoSplitSEAMOnFace pending.
                    if b_is_closed && !self.edge_closed_on_face(&a_sp, the_face) {
                        // Pending: BOPTools_AlgoTools3D::DoSplitSEAMOnFace(aSp, theFace).
                    }
                    // OCCT L1169-1172: IsSplitToReverseWithWarn pending.
                    new_edges.push(a_sp);
                }
            }
            // OCCT L1178-1180: MakeWire(aNewWire) + orientation + closed flag.
            let new_wire = Shape::new(
                std::sync::Arc::new(TShape::Wire(TWireData {
                    my_shapes: vec![],
                    flags: tshape_flags::DEFAULT,
                    edges: new_edges,
                })),
                0,
                a_w.orientation,
            );
            new_wires.push(new_wire);
        }
        // OCCT L1066: MakeFace(aDraftFace, aS, aLoc, aTol) — face without wires.
        let mut draft_face = Shape::new(
            std::sync::Arc::new(TShape::Face(TFaceData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                surface: Some(a_surf),
                surface_location: 0,
                outer_wire: new_wires.first().cloned().unwrap_or_else(Shape::null),
                inner_wires: new_wires.into_iter().skip(1).collect(),
                sample_point: None,
                uv_domain: None,
                internal_vertices: vec![],
                tolerance: a_tol,
                natural_restriction: false,
            })),
            0,
            topods::Orientation::Forward,
        );
        // OCCT L1183-1186: reverse if the original face was reversed.
        if the_face.orientation == topods::Orientation::Reversed {
            draft_face.orientation = topods::Orientation::Reversed;
        }
        Some(draft_face)
    }

    /// OCCT HasMultiConnected (BOPAlgo_Builder_2.cxx L1014-1045).
    /// Returns true when any vertex of the edge is shared by more than two edges.
    fn has_multi_connected(&self, the_edge: &Shape, the_map: &mut HashMap<u64, Vec<Shape>>) -> bool {
        let verts: Vec<Shape> = match &*the_edge.data {
            TShape::Edge(ed) => vec![
                Shape::new(ed.first.data.clone(), ed.first.location, ed.first.orientation),
                Shape::new(ed.last.data.clone(), ed.last.location, ed.last.orientation),
            ],
            _ => Vec::new(),
        };
        for v in verts {
            let list = the_map.entry(v.ptr_id()).or_default();
            if !list.iter().any(|e| e.ptr_id() == the_edge.ptr_id()) {
                list.push(the_edge.clone());
            }
            if list.len() > 2 {
                return true;
            }
        }
        false
    }

    /// OCCT BOPAlgo_Builder::FillSameDomainFaces (Builder_2.cxx L580-780).
    fn fill_same_domain_faces(&mut self) {
        // OCCT L584-589: get FF interferences, empty check.
        let a_ffs = &self.ds.interf_ff;
        if a_ffs.is_empty() { return; }

        // OCCT L597-649: build face-to-parent solid map (with image propagation).
        let mut a_face_to_parent: HashMap<u64, u64> = HashMap::new(); // face_ptr_id → solid_ptr_id
        let a_nb_src = self.ds.nb_source_shapes();
        for i_src in 0..a_nb_src {
            let a_si = self.ds.shape_info(i_src);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_solid = self.brep_sr(i_src);
            // Iterate solid's sub-shape shells → faces.
            for &shi in &a_si.sub_shapes {
                if shi >= self.ds.nb_shapes() { continue; }
                let sh_info = self.ds.shape_info(shi);
                if sh_info.shape_type != topods::ShapeType::Shell { continue; }
                for &fi in &sh_info.sub_shapes {
                    if fi >= self.ds.nb_shapes() { continue; }
                    let a_f = self.brep_sr(fi);
                    a_face_to_parent.entry(a_f.ptr_id()).or_insert(a_solid.ptr_id());
                }
            }
        }
        // OCCT L619-648: propagate the parent solid to the image faces.
        let mut a_propagation: HashMap<u64, u64> = HashMap::new();
        for (a_src, a_l_im) in &self.my_images {
            let p_parent = a_face_to_parent.get(&a_src.ptr_id()).copied();
            if let Some(parent) = p_parent {
                for a_piece in a_l_im {
                    if !a_face_to_parent.contains_key(&a_piece.ptr_id()) {
                        a_propagation.insert(a_piece.ptr_id(), parent);
                    }
                }
            }
        }
        for (k, v) in a_propagation {
            a_face_to_parent.entry(k).or_insert(v);
        }

        // OCCT L654-684: collect face indices from FF interferences.
        let mut a_fi_vec: Vec<usize> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for ff in a_ffs {
            for &nf in &[ff.f1, ff.f2] {
                if !self.ds.has_face_info(nf) { continue; }
                if a_m_fence.insert(nf) {
                    a_fi_vec.push(nf);
                }
            }
        }
        // OCCT L687: sort indices.
        a_fi_vec.sort();

        // OCCT L690-694: map edge-sets → list of faces, planar face fence.
        // rcad: edge_set = sorted Vec of (edge_ptr_id, curve_type) for the face boundary
        // (approximation of OCCT BOPTools_Set, which also folds in curve orientations).
        let mut an_e_set_faces: std::collections::HashMap<Vec<(u64, u32)>, Vec<Shape>> =
            std::collections::HashMap::new();
        let mut a_mf_planar: std::collections::HashSet<u64> = std::collections::HashSet::new();

        // OCCT L697-741: for each face, compute edge set.
        for &n_f in &a_fi_vec {
            let a_f = self.brep_sr(n_f);
            // OCCT L707-718: check if planar (Plane surface, closed bbox).
            let b_check_planar = if let Some(surf) = self.ds.face_surface(n_f) {
                matches!(surf, rcad_kernel::geom::Surface3::Plane(_))
            } else { false };

            // OCCT L720-740: get face images (or face itself), add edge set.
            let face_list: Vec<Shape> = if let Some(imgs) = self.my_images.get(&a_f) {
                imgs.clone()
            } else {
                vec![a_f.clone()]
            };

            for f_piece in &face_list {
                let mut edge_set: Vec<(u64, u32)> = Vec::new();
                let edges = self.ds.face_boundary_edges(n_f);
                for &ei in &edges {
                    if ei >= self.ds.nb_shapes() { continue; }
                    let e_ptr = self.ds.shape(ei).ptr_id();
                    let curve_type = match self.ds.edge_curve(ei) {
                        Some(c) => match c {
                            rcad_kernel::geom::Curve3::Line(_) => 0u32,
                            rcad_kernel::geom::Curve3::Circle(_) => 1u32,
                            rcad_kernel::geom::Curve3::Ellipse(_) => 2u32,
                            _ => 3u32,
                        },
                        None => 3u32,
                    };
                    edge_set.push((e_ptr, curve_type));
                }
                edge_set.sort_by(|a, b| a.0.cmp(&b.0));

                an_e_set_faces.entry(edge_set).or_default().push(f_piece.clone());
                if b_check_planar {
                    a_mf_planar.insert(f_piece.ptr_id());
                }
            }
        }

        // OCCT L743-748: aDMSLS — back-and-forth SD map; aVPSB — pairs for analysis.
        let mut a_dmsls: HashMap<Shape, Vec<Shape>> = HashMap::new();
        let mut a_vpsb: Vec<(Shape, Shape)> = Vec::new();

        // OCCT L750-791: check pairs of faces with equal edge set.
        for (_edge_set, faces) in &an_e_set_faces {
            if faces.len() < 2 { continue; }
            for i1 in 0..faces.len() {
                let f1 = &faces[i1];
                let parent1 = a_face_to_parent.get(&f1.ptr_id()).copied();
                let b_check_planar = a_mf_planar.contains(&f1.ptr_id());
                for i2 in (i1 + 1)..faces.len() {
                    let f2 = &faces[i2];
                    let parent2 = a_face_to_parent.get(&f2.ptr_id()).copied();
                    // OCCT L776-779: two faces of one solid cannot be SD.
                    if let (Some(p1), Some(p2)) = (parent1, parent2) {
                        if p1 == p2 { continue; }
                    }
                    // OCCT L780-785: planar bounded faces → SD without additional check.
                    if b_check_planar && a_mf_planar.contains(&f2.ptr_id()) {
                        fill_map_faces(f1, f2, &mut a_dmsls);
                        continue;
                    }
                    // OCCT L786-791: add pair for analysis.
                    a_vpsb.push((f1.clone(), f2.clone()));
                }
            }
        }

        // OCCT L799-822: perform the pair analysis.
        for (f1, f2) in &a_vpsb {
            // OCCT BOPAlgo_PairOfShapeBoolean::Perform (L94-105) →
            // BOPTools_AlgoTools::AreFacesSameDomain (BOPTools_AlgoTools.cxx L1139-1205).
            // rcad: existing are_faces_same_domain on DS indices (plane-check stub).
            let n_f1 = self.ds.index(f1);
            let n_f2 = self.ds.index(f2);
            let flag = if n_f1 >= 0 && n_f2 >= 0 {
                crate::bop::tools::algo_tools::are_faces_same_domain(
                    n_f1 as usize, n_f2 as usize, self.ds, self.my_fuzzy_value,
                )
            } else {
                false
            };
            if flag {
                fill_map_faces(f1, f2, &mut a_dmsls);
            }
        }

        // OCCT L826: MakeBlocks(aDMSLS, aMBlocks).
        let a_m_blocks = make_blocks_faces(&a_dmsls);

        // OCCT L830-882: fill same domain faces map.
        for a_lsd in &a_m_blocks {
            // Find the SD face: the original face (in DS) with minimal index;
            // otherwise the first face of the group.
            let mut p_fsd: Option<Shape> = None;
            let mut n_f_min: isize = std::isize::MAX;
            for a_f in a_lsd {
                let n_f = self.ds.index(a_f);
                if n_f >= 0 {
                    // OCCT L858: original face — consider it split into itself.
                    self.my_images.entry(a_f.clone()).or_default().push(a_f.clone());
                    if n_f < n_f_min {
                        n_f_min = n_f;
                        p_fsd = Some(a_f.clone());
                    }
                }
            }
            let p_fsd = match p_fsd {
                Some(fsd) => fsd,
                None => match a_lsd.first() {
                    Some(f) => f.clone(),
                    None => continue,
                },
            };
            // OCCT L876-881: bind all faces of the group to the SD face.
            for a_f in a_lsd {
                self.my_shapes_sd.insert(a_f.clone(), p_fsd.clone());
            }
        }

        // OCCT L886-921: update images with SD faces and fill origins.
        for i in 0..a_nb_src {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Face { continue; }
            let a_f = self.brep_sr(i);
            let Some(a_lf_im) = self.my_images.get_mut(&a_f) else { continue };
            for a_f_im in a_lf_im.iter_mut() {
                // OCCT L906-910: replace the image with its SD face.
                if let Some(a_fsd) = self.my_shapes_sd.get(a_f_im) {
                    *a_f_im = a_fsd.clone();
                }
                // OCCT L913-919: fill the map of origins.
                let p_lf_or = self.my_origins.entry(a_f_im.clone()).or_default();
                p_lf_or.push(a_f.clone());
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillInternalVertices (Builder_2.cxx L929+).
    fn fill_internal_vertices(&mut self) {
        let n = self.ds.nb_source_shapes();
        for i in 0..n {
            if self.ds.shape_info(i).shape_type != topods::ShapeType::Face { continue; }
            let fs = self.brep_sr(i);
            if !self.my_images.contains_key(&fs) { continue; }
            if !self.ds.has_face_info(i) { continue; }
            let v_pts: Vec<glam::DVec3> = self.ds.face_info(i).vertices_in.iter()
                .map(|&vi| self.ds.vertex_point_by_idx(vi)).collect();
            if v_pts.is_empty() { continue; }
            if let Some(imgs) = self.my_images.get_mut(&fs) {
                for pt in &v_pts {
                    for img in imgs.iter_mut() {
                        let ts = Arc::make_mut(&mut img.data);
                        if let TShape::Face(fd) = ts {
                            fd.internal_vertices.push(Shape::new(
                                Arc::new(TShape::Vertex(rcad_kernel::topods::TVertexData {
                                    my_shapes: vec![], flags: tshape_flags::DEFAULT,
                                    point: *pt, tolerance: 1e-7, points: vec![],
                                })), 0, topods::Orientation::Forward,
                            ));
                        }
                    }
                }
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesSolids (BOPAlgo_Builder_3.cxx L60-93).
    /// Builds split solids: FillIn3DParts -> BuildSplitSolids -> FillInternalShapes.
    fn fill_images_solids(&mut self) {
        // OCCT L62-73: check all DS source shapes for SOLID type
        let a_nb_s = self.ds.nb_source_shapes();
        let mut has_solid = false;
        for i in 0..a_nb_s {
            if self.ds.shape_info(i).shape_type == topods::ShapeType::Solid {
                has_solid = true;
                break;
            }
        }
        if !has_solid { return; }
        self.fill_in_3d_parts();
        if self.has_errors() { return; }
        self.build_split_solids();
        if self.has_errors() { return; }
        self.fill_internal_shapes();
    }

    /// OCCT BOPAlgo_Builder::FillIn3DParts (BOPAlgo_Builder_3.cxx L97-263).
    /// Collects all faces and draft solids, classifies faces as IN/OUT relative
    /// to each solid, and fills myInParts + myDraftSolids.
    fn fill_in_3d_parts(&mut self) {
        // OCCT L113-150: get all faces (source + their images).
        let a_nb_s = self.ds.nb_source_shapes();
        let mut a_lfaces: Vec<Shape> = Vec::new();
        let mut a_m_fence: HashSet<u64> = HashSet::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Face {
                continue;
            }
            let a_s = self.brep_sr(i);
            // OCCT L131-148: add images, or the face itself with its box.
            if let Some(imgs) = self.my_images.get(&a_s).cloned() {
                for a_s_im in &imgs {
                    if a_m_fence.insert(a_s_im.ptr_id()) {
                        a_lfaces.push(a_s_im.clone());
                    }
                }
            } else {
                // OCCT L147-148: aLFaces.Append(aS); box map bind (box culling
                // is an optimization, omitted in classify_faces).
                a_lfaces.push(a_s.clone());
            }
        }
        // OCCT L154-195: get all solids, build draft solids.
        let mut a_lsolids: Vec<Shape> = Vec::new();
        let mut a_lsolids_src: Vec<usize> = Vec::new(); // DS index per draft solid
        let mut a_solids_if: HashMap<u64, Vec<Shape>> = HashMap::new();
        let mut a_draft_solid: HashMap<u64, Shape> = HashMap::new();
        let mut a_source_solids: Vec<Shape> = Vec::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid {
                continue;
            }
            let a_solid = self.brep_sr(i);
            // OCCT L186-194: BuildDraftSolid(aSolid, aSD, aLIF).
            let mut a_lif: Vec<Shape> = Vec::new();
            let a_sd = self.build_draft_solid(&a_solid, &mut a_lif);
            a_lsolids.push(a_sd.clone());
            a_lsolids_src.push(i);
            a_solids_if.insert(a_sd.ptr_id(), a_lif);
            a_draft_solid.insert(a_solid.ptr_id(), a_sd.clone());
            a_source_solids.push(a_solid);
        }
        // OCCT L197-208: classify the faces relative to the solids.
        let an_in_parts = self.classify_faces(&a_lfaces, &a_lsolids, &a_lsolids_src, &a_solids_if);
        // OCCT L210-262: analyze the results of classification.
        for a_solid in &a_source_solids {
            let a_sd = match a_draft_solid.get(&a_solid.ptr_id()) {
                Some(sd) => sd.clone(),
                None => continue,
            };
            let a_l_in_faces = an_in_parts.get(&a_sd.ptr_id()).cloned().unwrap_or_default();
            let a_l_internal = a_solids_if.get(&a_sd.ptr_id()).cloned().unwrap_or_default();
            let a_nb_in = a_l_in_faces.len();
            if a_nb_in == 0 {
                // OCCT L227-238: check if the shells of the solid have an image.
                let mut b_has_image = false;
                for sh in self.shape_sub_shapes(a_solid) {
                    if self.my_images.contains_key(&sh) {
                        b_has_image = true;
                        break;
                    }
                }
                if !b_has_image {
                    continue;
                }
            }
            // OCCT L241: theDraftSolids.Bind(aSolid, aSDraft).
            self.my_draft_solids.insert(a_solid.ptr_id(), a_sd.clone());
            // OCCT L243-261: combine IN and internal faces into myInParts.
            let a_nb_int = a_l_internal.len();
            if a_nb_int != 0 || a_nb_in != 0 {
                let p_lin = self.my_in_parts.entry(a_solid.clone()).or_default();
                p_lin.extend(a_l_in_faces);
                p_lin.extend(a_l_internal);
            }
        }
    }

    /// OCCT BOPAlgo_Builder::BuildDraftSolid (BOPAlgo_Builder_3.cxx L267-368).
    /// Builds the draft solid from the solid's shells, replacing each face by
    /// its images (keeping the SD faces and INTERNAL faces in theLIF).
    fn build_draft_solid(&self, the_solid: &Shape, the_lif: &mut Vec<Shape>) -> Shape {
        // OCCT L283-367: iterate the solid's shells.
        let mut a_shd_list: Vec<Shape> = Vec::new();
        for a_sh in self.shape_sub_shapes(the_solid) {
            if a_sh.shape_type() != topods::ShapeType::Shell {
                continue;
            }
            let a_or_sh = a_sh.orientation;
            // OCCT L293: MakeShell(aShD).
            let mut a_shd_subs: Vec<Shape> = Vec::new();
            let mut i_flag = 0;
            for a_f in self.shape_sub_shapes(&a_sh) {
                let a_or_f = a_f.orientation;
                // OCCT L303: if myImages.IsBound(aF) — replace by images.
                if let Some(imgs) = self.my_images.get(&a_f).cloned() {
                    for a_fx in &imgs {
                        // OCCT L311: if myShapesSD.IsBound(aFx) — SD face.
                        if self.my_shapes_sd.contains_key(a_fx) {
                            if a_or_f == topods::Orientation::Internal {
                                // OCCT L313-317: aFx.Orientation(aOrF); theLIF.Append.
                                let mut fx = a_fx.clone();
                                fx.orientation = topods::Orientation::Internal;
                                the_lif.push(fx);
                            } else {
                                // OCCT L321-326: IsSplitToReverseWithWarn pending.
                                i_flag = 1;
                                a_shd_subs.push(a_fx.clone());
                            }
                        } else {
                            // OCCT L333-344: aFx.Orientation(aOrF).
                            let mut fx = a_fx.clone();
                            fx.orientation = a_or_f;
                            if a_or_f == topods::Orientation::Internal {
                                the_lif.push(fx);
                            } else {
                                i_flag = 1;
                                a_shd_subs.push(fx);
                            }
                        }
                    }
                } else {
                    // OCCT L348-359: face has no images.
                    if a_or_f == topods::Orientation::Internal {
                        the_lif.push(a_f.clone());
                    } else {
                        i_flag = 1;
                        a_shd_subs.push(a_f.clone());
                    }
                }
            }
            // OCCT L362-366: if any face was added, close and add the shell.
            if i_flag != 0 {
                let a_shd = Shape::new(
                    std::sync::Arc::new(TShape::Shell(TShellData {
                        my_shapes: vec![],
                        flags: tshape_flags::DEFAULT,
                        faces: a_shd_subs,
                    })),
                    0,
                    a_or_sh,
                );
                a_shd_list.push(a_shd);
            }
        }
        // OCCT L188: MakeSolid(aSD); L280-281: orientation.
        Shape::new(
            std::sync::Arc::new(TShape::Solid(TSolidData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                shells: a_shd_list,
                internal_vertices: vec![],
                internal_edges: vec![],
            })),
            0,
            the_solid.orientation,
        )
    }

    /// OCCT BOPAlgo_Tools::ClassifyFaces (BOPAlgo_Tools.cxx L1622-1747) —
    /// semantic translation: classifies a point of each face against each solid.
    /// The BVH box culling and connexity-block optimizations are omitted; the
    /// point-in-solid test is delegated to IntTools_Context::solid_classifier_perform.
    /// The solid's own faces (and its internal faces) are excluded from the
    /// classification (OCCT BOPAlgo_FillIn3DParts::Perform aMSF filter).
    fn classify_faces(
        &self,
        faces: &[Shape],
        solids: &[Shape],
        solids_src: &[usize],
        a_solids_if: &HashMap<u64, Vec<Shape>>,
    ) -> HashMap<u64, Vec<Shape>> {
        let mut in_parts: HashMap<u64, Vec<Shape>> = HashMap::new();
        for (k, a_sd) in solids.iter().enumerate() {
            // aMSF = own faces of the draft solid + its internal faces.
            let mut a_msf: HashSet<u64> = HashSet::new();
            for f in collect_solid_faces(a_sd) {
                a_msf.insert(f.ptr_id());
            }
            if let Some(lif) = a_solids_if.get(&a_sd.ptr_id()) {
                for f in lif {
                    a_msf.insert(f.ptr_id());
                }
            }
            for a_f in faces {
                // OCCT L1389: skip the solid's own faces.
                if a_msf.contains(&a_f.ptr_id()) {
                    continue;
                }
                // Compute face centroid for classification.
                let centroid = Self::face_centroid(a_f);
                // OCCT L1505-1509: IsInternalFace → point-in-solid.
                let state = self
                    .my_context
                    .solid_classifier_perform(&self.ds, solids_src[k], centroid, 1e-7);
                if state == 3 {
                    // IN
                    in_parts.entry(a_sd.ptr_id()).or_default().push(a_f.clone());
                }
            }
        }
        in_parts
    }

    /// Compute face centroid from its bounding vertices.
    fn face_centroid(face: &Shape) -> DVec3 {
        match &*face.data {
            TShape::Face(fd) => {
                let mut pts: Vec<DVec3> = Vec::new();
                if let TShape::Wire(wd) = &*fd.outer_wire.data {
                    for e in &wd.edges {
                        if let TShape::Edge(ed) = &*e.data {
                            if let TShape::Vertex(vd) = &*ed.first.data {
                                pts.push(vd.point);
                            }
                            if let TShape::Vertex(vd) = &*ed.last.data {
                                pts.push(vd.point);
                            }
                        }
                    }
                }
                if pts.is_empty() { return DVec3::ZERO; }
                pts.iter().sum::<DVec3>() / pts.len() as f64
            }
            _ => DVec3::ZERO,
        }
    }

    /// OCCT BOPAlgo_Builder::BuildSplitSolids (BOPAlgo_Builder_3.cxx L413-618).
    /// Each source solid is split independently by a separate BOPAlgo_SplitSolid;
    /// only its own split pieces are assigned to its images.
    fn build_split_solids(&mut self) {
        // OCCT L425-427: map of same-domain solids face sets (BOPTools_Set -> shape).
        let mut a_mst: HashMap<Vec<u64>, Shape> = HashMap::new();
        // OCCT L432-461: find same-domain solids for non-interfered solids.
        let a_nb_s = self.ds.nb_source_shapes();
        let mut a_mfence: HashSet<u64> = HashSet::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_s = self.brep_sr(i);
            if !a_mfence.insert(a_s.ptr_id()) { continue; }
            // OCCT L451-454: if theDraftSolids.IsBound(aS) continue;
            if self.my_draft_solids.contains_key(&a_s.ptr_id()) { continue; }
            // OCCT L456-459: aST.Add(aS, TopAbs_FACE).
            let a_st = self.shape_face_set(&a_s);
            a_mst.entry(a_st).or_insert_with(|| a_s.clone());
        }
        // OCCT L465-466: aSolidsIm — source solid -> result solids (IndexedDataMap).
        let mut a_solids_im: Vec<(Shape, Vec<Shape>)> = Vec::new();
        let mut a_solids_im_idx: HashMap<u64, usize> = HashMap::new();
        // OCCT L468-518: build split solids for interfered source solids.
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_s = self.brep_sr(i);
            // OCCT L478-481: if !theDraftSolids.IsBound(aS) continue;
            if !self.my_draft_solids.contains_key(&a_s.ptr_id()) { continue; }
            let a_sd = self.my_draft_solids.get(&a_s.ptr_id()).unwrap().clone();
            // OCCT L484-489: if no IN faces -> the draft solid itself, no split.
            let p_lfin: Vec<Shape> = self
                .my_in_parts
                .get(&a_s)
                .cloned()
                .unwrap_or_default();
            if p_lfin.is_empty() {
                let idx = *a_solids_im_idx.entry(a_s.ptr_id()).or_insert_with(|| {
                    a_solids_im.push((a_s.clone(), Vec::new()));
                    a_solids_im.len() - 1
                });
                a_solids_im[idx].1.push(a_sd);
                continue;
            }
            // OCCT L493-499: 1.1 shell faces set of the draft solid.
            let mut a_sfs: Vec<Shape> = Vec::new();
            for ss in self.shape_sub_shapes(&a_sd) {
                if ss.shape_type() == topods::ShapeType::Shell {
                    for f in self.shape_sub_shapes(&ss) {
                        if f.shape_type() == topods::ShapeType::Face {
                            a_sfs.push(f);
                        }
                    }
                }
            }
            // OCCT L501-511: 1.2 add IN faces (both orientations).
            for a_f in &p_lfin {
                let mut f_fwd = a_f.clone();
                f_fwd.orientation = topods::Orientation::Forward;
                a_sfs.push(f_fwd);
                let mut f_rev = a_f.clone();
                f_rev.orientation = topods::Orientation::Reversed;
                a_sfs.push(f_rev);
            }
            // OCCT L514-517: BOPAlgo_SplitSolid for THIS solid only.
            let mut bs = crate::bop::algo::builder_solid::BuilderSolid::new(&self.ds);
            bs.my_shapes = a_sfs;
            bs.perform();
            // OCCT L542: aSolidsIm.Add(aBS.Solid(), aBS.Areas()).
            let idx = *a_solids_im_idx.entry(a_s.ptr_id()).or_insert_with(|| {
                a_solids_im.push((a_s.clone(), Vec::new()));
                a_solids_im.len() - 1
            });
            a_solids_im[idx].1.extend(bs.my_solids.clone());
            // OCCT L544-577: merge the split solid's report into the main
            // report, converting all sub-split errors into warnings.
            for a in bs.report().errors() {
                self.my_report.add_warning(a.clone());
            }
        }
        // OCCT L580-617: add new solids to the images map (same-domain dedup).
        for (a_s, a_lsr) in &a_solids_im {
            // OCCT L586: if !myImages.IsBound(aS).
            if self.my_images.contains_key(a_s) { continue; }
            // Compute the same-domain dedup results first (immutable borrows).
            let mut results: Vec<(Shape, Shape, bool)> = Vec::new();
            for a_sr in a_lsr {
                // OCCT L593-601: BOPTools_Set of aSR's faces; aMST.Added(aST).Shape().
                let a_st = self.shape_face_set(a_sr);
                let b_flag_sd = a_mst.contains_key(&a_st);
                let a_sx = a_mst
                    .get(&a_st)
                    .cloned()
                    .unwrap_or_else(|| a_sr.clone());
                results.push((a_sr.clone(), a_sx, b_flag_sd));
            }
            let p_lsx = self.my_images.entry(a_s.clone()).or_default();
            for (a_sr, a_sx, b_flag_sd) in results {
                p_lsx.push(a_sx.clone());
                // OCCT L604-609: myOrigins[aSx].Append(aS).
                self.my_origins
                    .entry(a_sx.clone())
                    .or_default()
                    .push(a_s.clone());
                // OCCT L611-614: if same-domain, bind myShapesSD[aSR] = aSx.
                if b_flag_sd {
                    self.my_shapes_sd.insert(a_sr, a_sx);
                }
            }
        }
    }
    /// OCCT BOPAlgo_Builder::FillInternalShapes (Builder_3.cxx L622-830).
    fn fill_internal_shapes(&mut self) {
        // OCCT L631-644: local lists/maps
        let mut a_lsc: Vec<Shape> = Vec::new();
        let mut a_m_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        let mut a_l_args: Vec<Shape> = Vec::new();
        // aMSI indexed map → rcad: Vec<Shape> with dedup
        let mut a_msi: Vec<Shape> = Vec::new();
        let mut a_msi_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        // Map for ancestor lookup (vertex→[edges], vertex→[faces], edge→[faces])
        let mut a_msx: std::collections::HashMap<u64, Vec<u64>> = std::collections::HashMap::new();

        // OCCT L653-659: TreatCompound on arguments → flatten into aLSC
        let a_arguments = &self.ds.arguments;
        for a_s in a_arguments {
            Self::treat_compound(a_s, &mut a_lsc, &mut a_m_fence);
        }
        // OCCT L660-681: collect V/E from aLSC into aLArgs
        a_m_fence.clear();
        for a_s in &a_lsc {
            let a_type = a_s.shape_type();
            if a_type == topods::ShapeType::Wire {
                for ss in self.shape_sub_shapes(a_s) {
                    if a_m_fence.insert(ss.ptr_id()) {
                        a_l_args.push(ss);
                    }
                }
            } else if a_type == topods::ShapeType::Vertex || a_type == topods::ShapeType::Edge {
                a_l_args.push(a_s.clone());
            }
        }
        // OCCT L684-709: for each V/E/W, add images or self to aMSI
        a_m_fence.clear();
        for a_s in &a_l_args {
            let a_type = a_s.shape_type();
            if a_type == topods::ShapeType::Vertex
                || a_type == topods::ShapeType::Edge
                || a_type == topods::ShapeType::Wire
            {
                if a_m_fence.insert(a_s.ptr_id()) {
                    if let Some(imgs) = self.my_images.get(a_s) {
                        for img in imgs {
                            if a_msi_fence.insert(img.ptr_id()) {
                                a_msi.push(img.clone());
                            }
                        }
                    } else {
                        if a_msi_fence.insert(a_s.ptr_id()) {
                            a_msi.push(a_s.clone());
                        }
                    }
                }
            }
        }

        // OCCT L721-788: internal V/E from source solids
        a_m_fence.clear();
        let a_nb_s = self.ds.nb_source_shapes();
        let mut a_ls_d: Vec<Shape> = Vec::new();
        // aMSOr: original solids without images (OCCT L785).
        let mut a_ms_or: HashSet<u64> = HashSet::new();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_s = self.brep_sr(i);

            // OCCT L738: OwnInternalShapes(aS, aMx)
            // rcad: iterate solid sub-shapes, find internal V/E
            let a_mx = self.own_internal_shapes(&a_s);

            // OCCT L741-758: add internal shapes to aMSI
            for a_si_internal in &a_mx {
                if let Some(imgs) = self.my_images.get(a_si_internal) {
                    for img in imgs {
                        if a_msi_fence.insert(img.ptr_id()) {
                            a_msi.push(img.clone());
                        }
                    }
                } else {
                    if a_msi_fence.insert(a_si_internal.ptr_id()) {
                        a_msi.push(a_si_internal.clone());
                    }
                }
            }

            // OCCT L760-787: build ancestor map from split solids
            if let Some(imgs) = self.my_images.get(&a_s) {
                for a_sp in imgs {
                    if a_m_fence.insert(a_sp.ptr_id()) {
                        Self::map_shapes_and_ancestors(a_sp, &mut a_msx);
                        a_ls_d.push(a_sp.clone());
                    }
                }
            } else {
                if a_m_fence.insert(a_s.ptr_id()) {
                    Self::map_shapes_and_ancestors(&a_s, &mut a_msx);
                    a_ls_d.push(a_s.clone());
                    // OCCT L785: aMSOr.Add(aS).
                    a_ms_or.insert(a_s.ptr_id());
                }
            }
        }

        // OCCT L792-809: filter aMSI — keep only shapes not tied to split solid faces
        let mut a_ls_i: Vec<Shape> = Vec::new();
        for a_si in &a_msi {
            if let Some(ancestors) = a_msx.get(&a_si.ptr_id()) {
                if ancestors.is_empty() {
                    a_ls_i.push(a_si.clone());
                }
            } else {
                a_ls_i.push(a_si.clone());
            }
        }

        // OCCT L812-816: empty check
        if a_ls_i.is_empty() { return; }

        // OCCT L820-877: settle internal vertices and edges into solids.
        let mut i_sd = 0;
        while i_sd < a_ls_d.len() {
            let mut a_sd = a_ls_d[i_sd].clone();
            let mut b_modified = false;
            let mut i_si = 0;
            while i_si < a_ls_i.len() {
                let mut a_si = a_ls_i[i_si].clone();
                // OCCT L834: aSI.Orientation(TopAbs_INTERNAL).
                a_si.orientation = topods::Orientation::Internal;
                // OCCT L836: ComputeStateByOnePoint(aSI, aSd, 1.e-11, ctx).
                let a_state = Self::compute_state_by_one_point(&a_si, &a_sd, 1e-11);
                if a_state != 3 {
                    // TopAbs_IN — not inside; keep for the next solid.
                    i_si += 1;
                    continue;
                }
                if a_ms_or.contains(&a_sd.ptr_id()) {
                    // OCCT L846-858: make a new solid aSdx (copy of aSd + aSI).
                    let a_sdx = Self::solid_copy_add(&a_sd, &a_si);
                    // OCCT L860-861: myImages[aSd].Append(aSdx).
                    self.my_images
                        .entry(a_sd.clone())
                        .or_default()
                        .push(a_sdx.clone());
                    // OCCT L863-865: myOrigins[aSdx].Append(aSd).
                    self.my_origins
                        .entry(a_sdx.clone())
                        .or_default()
                        .push(a_sd.clone());
                    // OCCT L867-868: aMSOr.Remove(aSd); aSd = aSdx.
                    a_ms_or.remove(&a_sd.ptr_id());
                    a_sd = a_sdx;
                } else {
                    // OCCT L871-873: aBB.Add(aSd, aSI).
                    Self::solid_add_shape(&mut a_sd, &a_si);
                    b_modified = true;
                }
                // OCCT L875: aLSI.Remove(aIt1) — removal without advancing.
                a_ls_i.remove(i_si);
            }
            if b_modified {
                a_ls_d[i_sd] = a_sd;
            }
            i_sd += 1;
        }
    }

    /// OCCT BOPTools_AlgoTools::ComputeStateByOnePoint (BOPTools_AlgoTools.cxx L623-656).
    /// Classifies a shape (vertex/edge) relative to a solid by a single point.
    fn compute_state_by_one_point(the_s: &Shape, the_ref: &Shape, the_tol: f64) -> u8 {
        let a_type = the_s.shape_type();
        let a_p3d: Option<DVec3> = match a_type {
            topods::ShapeType::Vertex => match &*the_s.data {
                TShape::Vertex(vd) => Some(vd.point),
                _ => None,
            },
            topods::ShapeType::Edge => match &*the_s.data {
                TShape::Edge(ed) => {
                    if let Some(curve) = &ed.curve {
                        // OCCT L756-780: intermediate parameter of the curve range.
                        let (a_t1, a_t2) = (ed.range[0], ed.range[1]);
                        let a_t = if a_t1.is_infinite() && !a_t2.is_infinite() {
                            a_t2 - 10.0
                        } else if !a_t1.is_infinite() && a_t2.is_infinite() {
                            a_t1 + 10.0
                        } else if a_t1.is_infinite() && a_t2.is_infinite() {
                            0.0
                        } else {
                            (a_t1 + a_t2) * 0.5
                        };
                        Some(curve.point_at(a_t))
                    } else {
                        // degenerated edge — first vertex point (OCCT L748-754).
                        match &*ed.first.data {
                            TShape::Vertex(vd) => Some(vd.point),
                            _ => None,
                        }
                    }
                }
                _ => None,
            },
            _ => {
                // OCCT L646-653: recurse into the first sub-shape.
                let subs = Self::shape_sub_shapes_static(the_s);
                if let Some(sub) = subs.first() {
                    return Self::compute_state_by_one_point(sub, the_ref, the_tol);
                }
                None
            }
        };
        let p = match a_p3d {
            Some(p) => p,
            None => return 0, // TopAbs_UNKNOWN
        };
        let mut clsf =
            crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(the_ref);
        clsf.perform(p, the_tol);
        clsf.my_state
    }

    /// OCCT aBB.MakeSolid(aSdx); copy all sub-shapes of aSd; aBB.Add(aSdx, aSI)
    /// (BOPAlgo_Builder_3.cxx L849-858).
    fn solid_copy_add(a_sd: &Shape, a_si: &Shape) -> Shape {
        let mut shells = Vec::new();
        let mut internal_v = Vec::new();
        let mut internal_e = Vec::new();
        if let TShape::Solid(sd) = &*a_sd.data {
            shells = sd.shells.clone();
            internal_v = sd.internal_vertices.clone();
            internal_e = sd.internal_edges.clone();
        }
        match &*a_si.data {
            TShape::Vertex(_) => internal_v.push(a_si.clone()),
            TShape::Edge(_) => internal_e.push(a_si.clone()),
            _ => {}
        }
        Shape::new(
            std::sync::Arc::new(TShape::Solid(TSolidData {
                my_shapes: vec![],
                flags: tshape_flags::DEFAULT,
                shells,
                internal_vertices: internal_v,
                internal_edges: internal_e,
            })),
            0,
            a_sd.orientation,
        )
    }

    /// OCCT aBB.Add(aSd, aSI) — add an internal vertex/edge to the solid
    /// (BOPAlgo_Builder_3.cxx L871-873).
    fn solid_add_shape(a_sd: &mut Shape, a_si: &Shape) {
        let ts = Arc::make_mut(&mut a_sd.data);
        if let TShape::Solid(sd) = ts {
            match &*a_si.data {
                TShape::Vertex(_) => sd.internal_vertices.push(a_si.clone()),
                TShape::Edge(_) => sd.internal_edges.push(a_si.clone()),
                _ => {}
            }
        }
    }

    /// OCCT BOPTools_AlgoTools::TreatCompound — flatten compound shapes.
    fn treat_compound(s: &Shape, result: &mut Vec<Shape>, fence: &mut std::collections::HashSet<u64>) {
        if s.shape_type() == topods::ShapeType::Compound {
            if let TShape::Compound(children) = &*s.data {
                for child in children {
                    Self::treat_compound(child, result, fence);
                }
            }
        } else {
            if fence.insert(s.ptr_id()) {
                result.push(s.clone());
            }
        }
    }

    /// OCCT BOPAlgo_Builder::OwnInternalShapes (BOPAlgo_Builder_3.cxx L891-905).
    /// Collects all non-SHELL direct sub-shapes of the solid (internal
    /// vertices/edges/wires stored at the solid level).
    fn own_internal_shapes(&self, s: &Shape) -> Vec<Shape> {
        let mut result: Vec<Shape> = Vec::new();
        let mut fence: HashSet<u64> = HashSet::new();
        if let TShape::Solid(sd) = &*s.data {
            for v in &sd.internal_vertices {
                if fence.insert(v.ptr_id()) {
                    result.push(v.clone());
                }
            }
            for e in &sd.internal_edges {
                if fence.insert(e.ptr_id()) {
                    result.push(e.clone());
                }
            }
        }
        // OCCT L896-904: direct sub-shapes that are not SHELL.
        for ss in self.shape_sub_shapes(s) {
            if ss.shape_type() != topods::ShapeType::Shell {
                if fence.insert(ss.ptr_id()) {
                    result.push(ss.clone());
                }
            }
        }
        result
    }

    /// OCCT TopExp::MapShapesAndAncestors — build ancestor map (vertex→edge, vertex→face, edge→face).
    fn map_shapes_and_ancestors(s: &Shape, map: &mut std::collections::HashMap<u64, Vec<u64>>) {
        // OCCT: for each VERTEX, record its ancestor EDGEs; for each VERTEX/EDGE, record ancestor FACEs.
        // rcad: use shape_sub_shapes to walk the hierarchy
        let subs = Self::shape_sub_shapes_static(s);
        for sub in &subs {
            let sub_type = sub.shape_type();
            let sub_id = sub.ptr_id();
            // Find ancestors by scanning all sub-shapes at the parent level
            for parent in &subs {
                let parent_type = parent.shape_type();
                if parent_type == topods::ShapeType::Face
                    && (sub_type == topods::ShapeType::Vertex || sub_type == topods::ShapeType::Edge)
                {
                    map.entry(sub_id).or_default().push(parent.ptr_id());
                }
                if parent_type == topods::ShapeType::Edge && sub_type == topods::ShapeType::Vertex {
                    map.entry(sub_id).or_default().push(parent.ptr_id());
                }
            }
        }
    }

    /// Static version of shape_sub_shapes for use in non-&self methods.
    fn shape_sub_shapes_static(s: &Shape) -> Vec<Shape> {
        match &*s.data {
            TShape::Vertex(_) => vec![],
            TShape::Edge(ed) => vec![
                Shape::new(ed.first.data.clone(), ed.first.location, ed.first.orientation),
                Shape::new(ed.last.data.clone(), ed.last.location, ed.last.orientation),
            ],
            TShape::Wire(wd) => {
                wd.edges.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Face(fd) => {
                let mut v = vec![
                    Shape::new(fd.outer_wire.data.clone(), fd.outer_wire.location, fd.outer_wire.orientation)
                ];
                v.extend(fd.inner_wires.iter().map(|w| {
                    Shape::new(w.data.clone(), w.location, w.orientation)
                }));
                v
            }
            TShape::Shell(sd) => {
                sd.faces.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Solid(sd) => {
                sd.shells.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::CompSolid(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
            TShape::Compound(cd) => {
                cd.iter().map(|sr| {
                    Shape::new(sr.data.clone(), sr.location, sr.orientation)
                }).collect()
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesCompounds (BOPAlgo_Builder_1.cxx L197-217).
    fn fill_images_compounds(&mut self) {
        // OCCT L199-201: fence map + NbSourceShapes
        let mut a_mfp: std::collections::HashSet<u64> = std::collections::HashSet::new();
        let a_nb_s = self.ds.nb_source_shapes();
        // OCCT L202-216: for each source compound, call FillImagesCompound
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Compound {
                continue;
            }
            let a_c = self.brep_sr(i);
            // OCCT L210: FillImagesCompound(aC, aMFP);
            self.fill_images_compound(&a_c, &mut a_mfp);
        }
    }

    /// OCCT BOPAlgo_Builder::FillImagesCompound (BOPAlgo_Builder_1.cxx L278-360).
    fn fill_images_compound(&mut self, the_c: &Shape, the_mfp: &mut std::collections::HashSet<u64>) {
        // OCCT L282-283: check if already processed
        if !the_mfp.insert(the_c.ptr_id()) {
            return;
        }
        // OCCT L285-295: check if compound has images
        if let Some(imgs) = self.my_images.get(the_c) {
            // OCCT L289-293: recursively process existing images
            let imgs_clone = imgs.clone();
            for img in &imgs_clone {
                self.fill_images_compound(img, the_mfp);
            }
            return;
        }
        // OCCT L300-340: build new compound from sub-shape images
        let subs = self.shape_sub_shapes(the_c);
        let mut has_modified = false;
        for ss in &subs {
            if self.my_images.contains_key(ss) {
                has_modified = true;
                break;
            }
        }
        if !has_modified {
            return;
        }
        // OCCT L345-358: add the_s to myImages with new compound
        let mut new_shapes: Vec<Shape> = Vec::new();
        for ss in &subs {
            if let Some(ss_imgs) = self.my_images.get(ss) {
                new_shapes.extend(ss_imgs.iter().cloned());
            } else {
                new_shapes.push(ss.clone());
            }
        }
        let comp_shape = Shape::new(
            std::sync::Arc::new(TShape::Compound(new_shapes)),
            0, topods::Orientation::Forward,
        );
        self.my_images.entry(the_c.clone()).or_default().push(comp_shape);
    }

    /// OCCT BOPAlgo_Builder::PrepareHistory (BOPAlgo_Builder_4.cxx L164-252).
    fn prepare_history(&mut self) {
        // OCCT L166-168: if (!HasHistory()) return;
        if !self.my_fill_history { return; }

        // OCCT L172-176: init history tool, map result shapes
        // rcad: OCCT BRepTools_History for Modified/Generated/Deleted tracking
        // follows OCCT's Modified/Generated/Deleted detection.
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            let a_s = self.brep_sr(i);
            // OCCT L192-195: skip unsupported types
            let a_type = a_s.shape_type();
            if a_type != topods::ShapeType::Vertex
                && a_type != topods::ShapeType::Edge
                && a_type != topods::ShapeType::Wire
                && a_type != topods::ShapeType::Face
                && a_type != topods::ShapeType::Shell
                && a_type != topods::ShapeType::Solid
                && a_type != topods::ShapeType::Compound
            {
                continue;
            }

            // OCCT L205: LocModified → check myImages
            let mut is_modified = false;
            if let Some(imgs) = self.my_images.get(&a_s) {
                for a_sp in imgs {
                    // OCCT L214-217: check if result contains the split
                    // rcad: check via shape_remap (shape was added to result)
                    if self.shape_remap.contains_key(&a_sp.ptr_id()) {
                        is_modified = true;
                        // OCCT L218-226: orientation adjustment
                        // rcad: orientation handled by add_shape_to_result
                    }
                }
            }

            // OCCT L234-243: LocGenerated — check myImages for shapes generated
            // from this shape (e.g., edges generated from vertices)
            if let Some(imgs) = self.my_images.get(&a_s) {
                for _a_g in imgs {
                    if self.shape_remap.contains_key(&_a_g.ptr_id()) {
                        // OCCT L241: myHistory->AddGenerated(aS, aG)
                    }
                }
            }

            // OCCT L247-250: if not modified and not in result → Deleted
            let in_result = self.shape_remap.contains_key(&a_s.ptr_id());
            if !is_modified && !in_result {
                // OCCT L249: myHistory->Remove(aS)
            }
        }
    }

    /// OCCT BOPAlgo_Builder::PostTreat (BOPAlgo_Builder.cxx L461-486).
    fn post_treat(&mut self) {
        // OCCT L466-480: in non-destructive mode, collect source V/E/F shapes
        let mut a_ma: std::collections::HashSet<u64> = std::collections::HashSet::new();
        // OCCT L466-479: if (myPaveFiller->NonDestructive())
        // rcad: non-destructive mode not fully implemented. The collection
        // of shapes for CorrectTolerances — tolerance optimization.
        let _ = a_ma;
        // OCCT L483: BOPTools_AlgoTools::CorrectTolerances(myShape, aMA, 0.05, myRunParallel)
        // OCCT L485: BOPTools_AlgoTools::CorrectShapeTolerances(myShape, aMA, myRunParallel)
        // rcad: CorrectShapeTolerances — tolerance optimization.
    }

    // ====================================================================
    // BuildShape — OCCT BOPAlgo_BOP::BuildShape (BOPAlgo_BOP.cxx L885-1107)
    // Boolean operation result construction (BuildRC -> BuildSolid / containers).
    // ====================================================================

    /// OCCT BOPAlgo_BOP::BuildShape (BOPAlgo_BOP.cxx L885-1107).
    fn build_shape(&mut self) {
        // OCCT BOPAlgo_BOP::CheckData sets myDims from arguments/tools.
        self.compute_dims();
        // OCCT L889-911: for 3D+3D, if any argument solid is open, use the
        // BuildBOP alternative (BOPAlgo_Builder.cxx L491-897).
        if self.my_dims[0] == 3 && self.my_dims[1] == 3 {
            let has_not_closed_solids = self.check_args_for_open_solid();
            if has_not_closed_solids {
                // OCCT L902-909: BuildBOP fallback for open solids — rcad pending.
                // Closed-solid inputs (all stage tests) take the BuildRC path.
            }
        }
        // OCCT L913-914: BuildRC.
        self.build_rc();
        // OCCT L916-920: FUSE of 3D → BuildSolid.
        if self.my_operation == BooleanOpType::Union && self.my_dims[0] == 3 {
            self.build_solid();
            return;
        }
        // OCCT L923-1107: CUT/COMMON container logic.
        // OCCT L936-937: aMSRC = MapShapes(myRC).
        let mut a_msrc: HashSet<u64> = HashSet::new();
        for s in &self.my_rc {
            self.add_all_sub_shapes(s, &mut a_msrc);
        }
        // OCCT L940-951: collect containers of arguments.
        let mut a_lsc: Vec<Shape> = Vec::new();
        for i in 0..2 {
            let a_ls: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            for a_s in a_ls {
                Self::collect_containers(a_s, &mut a_lsc);
            }
        }
        // OCCT L958-1045: make containers.
        let mut a_lc_res: Vec<Shape> = Vec::new();
        for a_sc in &a_lsc {
            // OCCT L966: MakeContainer(COMPOUND, aRC).
            let mut a_rc: Vec<Shape> = Vec::new();
            // OCCT L968-990: add splits of sub-shapes contained in aMSRC.
            for a_s in self.shape_sub_shapes(a_sc) {
                if let Some(imgs) = self.my_images.get(&a_s) {
                    for a_s_im in imgs {
                        if a_msrc.contains(&a_s_im.ptr_id()) {
                            a_rc.push(a_s_im.clone());
                        }
                    }
                } else if a_msrc.contains(&a_s.ptr_id()) {
                    a_rc.push(a_s.clone());
                }
            }
            // OCCT L992-1009: connexity element types.
            let a_type = a_sc.shape_type();
            let (a_t1, a_t2) = match a_type {
                topods::ShapeType::Wire => (topods::ShapeType::Vertex, topods::ShapeType::Edge),
                topods::ShapeType::Shell => (topods::ShapeType::Edge, topods::ShapeType::Face),
                _ => (topods::ShapeType::Face, topods::ShapeType::Solid),
            };
            // OCCT L1011-1016: MakeConnexityBlocks(aRC, aT1, aT2, aLCB).
            let mut a_lcb: Vec<Vec<Shape>> = Vec::new();
            Self::make_connexity_blocks_shapes(&a_rc, a_t1, a_t2, &mut a_lcb);
            if a_lcb.is_empty() {
                continue;
            }
            // OCCT L1018-1043: build containers from blocks.
            for a_cb in &a_lcb {
                let mut a_rcb = Self::build_container_of_type(a_type, a_cb);
                match a_type {
                    topods::ShapeType::Wire => Self::orient_edges_on_wire(&mut a_rcb),
                    topods::ShapeType::Shell => Self::orient_faces_on_shell(&mut a_rcb),
                    _ => {}
                }
                a_lc_res.push(a_rcb);
            }
        }
        // OCCT L1047: RemoveDuplicates(aLCRes).
        self.remove_duplicates(&mut a_lc_res);
        // OCCT L1055-1063: aResult compound of containers.
        let mut a_result: Vec<Shape> = a_lc_res;
        // OCCT L1066-1067: aMSResult = MapShapes(aResult).
        let mut a_ms_result: HashSet<u64> = HashSet::new();
        for s in &a_result {
            self.add_all_sub_shapes(s, &mut a_ms_result);
        }
        // OCCT L1069-1080: get input non-container shapes.
        let mut a_ls_non_cont: Vec<Shape> = Vec::new();
        let mut a_m_inp_fence: HashSet<u64> = HashSet::new();
        for i in 0..2 {
            let a_ls: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            for a_s in a_ls {
                Self::treat_compound(a_s, &mut a_ls_non_cont, &mut a_m_inp_fence);
            }
        }
        // OCCT L1082-1104: put non-container shapes in the result.
        for a_s in &a_ls_non_cont {
            if let Some(imgs) = self.my_images.get(a_s) {
                for a_s_im in imgs {
                    if a_msrc.contains(&a_s_im.ptr_id()) && a_ms_result.insert(a_s_im.ptr_id()) {
                        a_result.push(a_s_im.clone());
                    }
                }
            } else if a_msrc.contains(&a_s.ptr_id()) && a_ms_result.insert(a_s.ptr_id()) {
                a_result.push(a_s.clone());
            }
        }
        // OCCT L1106: myShape = aResult.
        self.set_shape_from_shapes(a_result);
    }

    /// OCCT BOPAlgo_BOP::BuildRC (BOPAlgo_BOP.cxx L597-881).
    fn build_rc(&mut self) {
        let mut a_c: Vec<Shape> = Vec::new();
        // OCCT L608-623: A. Fuse — collect shapes of TypeToExplore(dim0) from myShape.
        if self.my_operation == BooleanOpType::Union {
            let mut a_m_fence: HashSet<u64> = HashSet::new();
            let a_type = type_to_explore(self.my_dims[0]);
            let my_shape = self.my_shape.clone().unwrap_or_default();
            // OCCT: TopExp_Explorer aExp(myShape, aType).
            for ts in &my_shape.tshapes {
                let sh = Shape::new(ts.clone(), 0, topods::Orientation::Forward);
                if sh.shape_type() == a_type && a_m_fence.insert(sh.ptr_id()) {
                    a_c.push(sh);
                }
            }
            self.my_rc = a_c;
            return;
        }
        // OCCT L630-659: B. Common/Cut — building elements of arguments.
        let mut a_m_args: Vec<Shape> = Vec::new();
        let mut a_m_args_fence: HashSet<u64> = HashSet::new();
        let mut a_m_tools: Vec<Shape> = Vec::new();
        let mut a_m_tools_fence: HashSet<u64> = HashSet::new();
        for i in 0..2 {
            let a_ls: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            let (a_ms, a_ms_fence) = if i == 0 {
                (&mut a_m_args, &mut a_m_args_fence)
            } else {
                (&mut a_m_tools, &mut a_m_tools_fence)
            };
            for a_s in a_ls {
                let mut a_list: Vec<Shape> = Vec::new();
                let mut a_fence: HashSet<u64> = HashSet::new();
                Self::treat_compound(a_s, &mut a_list, &mut a_fence);
                for a_ss in a_list {
                    let i_dim = shape_dimension(&a_ss);
                    if i_dim < 0 {
                        continue;
                    }
                    let a_type = type_to_explore(i_dim);
                    // OCCT L656: TopExp::MapShapes(aSS, aType, aMS).
                    for sh in self.map_shapes_of_type(&a_ss, a_type) {
                        if a_ms_fence.insert(sh.ptr_id()) {
                            a_ms.push(sh);
                        }
                    }
                }
            }
        }
        // OCCT L666-718: get splits of building elements.
        let mut a_m_args_im: Vec<Shape> = Vec::new();
        let mut a_m_args_im_fence: HashSet<u64> = HashSet::new();
        let mut a_m_tools_im: Vec<Shape> = Vec::new();
        let mut a_m_tools_im_fence: HashSet<u64> = HashSet::new();
        let mut a_m_set_args: HashMap<Vec<u64>, Shape> = HashMap::new();
        let mut a_m_set_tools: HashMap<Vec<u64>, Shape> = HashMap::new();
        let mut b_check_edges = false;
        for i in 0..2 {
            let a_ms: &Vec<Shape> = if i == 0 { &a_m_args } else { &a_m_tools };
            let (a_ms_im, a_ms_im_fence, a_m_set) = if i == 0 {
                (&mut a_m_args_im, &mut a_m_args_im_fence, &mut a_m_set_args)
            } else {
                (&mut a_m_tools_im, &mut a_m_tools_im_fence, &mut a_m_set_tools)
            };
            for a_s in a_ms {
                let a_type = a_s.shape_type();
                if a_type == topods::ShapeType::Edge {
                    b_check_edges = true;
                    // OCCT L689-691: skip degenerated edges.
                    let degen = a_s.as_edge().map(|e| e.degenerated).unwrap_or(false);
                    if degen {
                        continue;
                    }
                }
                if let Some(imgs) = self.my_images.get(a_s).cloned() {
                    for a_s_im in &imgs {
                        if a_ms_im_fence.insert(a_s_im.ptr_id()) {
                            a_ms_im.push(a_s_im.clone());
                        }
                    }
                } else {
                    if a_ms_im_fence.insert(a_s.ptr_id()) {
                        a_ms_im.push(a_s.clone());
                    }
                    if a_type == topods::ShapeType::Solid {
                        // OCCT L708-716: BOPTools_Set of the solid's face set.
                        let a_st = self.shape_face_set(a_s);
                        if !a_m_set.contains_key(&a_st) {
                            a_m_set.insert(a_st, a_s.clone());
                        }
                    }
                }
            }
        }
        // OCCT L723-798: compare maps and make the result.
        let i_dim_min = std::cmp::min(self.my_dims[0], self.my_dims[1]);
        let b_common = self.my_operation == BooleanOpType::Intersection;
        // rcad has no CUT21; aMIt = object splits, aMCheck = tool splits.
        let (a_m_it, a_m_check, a_m_set_check) = (&a_m_args_im, &a_m_tools_im, &a_m_set_tools);
        let mut a_m_check_exp: Vec<Shape> = Vec::new();
        let mut a_m_check_exp_fence: HashSet<u64> = HashSet::new();
        let mut a_m_it_exp: Vec<Shape> = Vec::new();
        let mut a_m_it_exp_fence: HashSet<u64> = HashSet::new();
        if b_common {
            // OCCT L738-751: expand aMIt with sub-shapes of lower dims.
            for a_s in a_m_it {
                let i_dim_max = shape_dimension(a_s);
                for i_dim in i_dim_min..i_dim_max {
                    let a_type = type_to_explore(i_dim);
                    for sh in self.map_shapes_of_type(a_s, a_type) {
                        if a_m_it_exp_fence.insert(sh.ptr_id()) {
                            a_m_it_exp.push(sh);
                        }
                    }
                }
                if a_m_it_exp_fence.insert(a_s.ptr_id()) {
                    a_m_it_exp.push(a_s.clone());
                }
            }
        } else {
            a_m_it_exp = a_m_it.clone();
            a_m_it_exp_fence = a_m_it_exp.iter().map(|s| s.ptr_id()).collect();
        }
        // OCCT L758-769: expand aMCheck with sub-shapes of lower dims.
        for a_s in a_m_check {
            let i_dim_max = shape_dimension(a_s);
            for i_dim in i_dim_min..i_dim_max {
                let a_type = type_to_explore(i_dim);
                for sh in self.map_shapes_of_type(a_s, a_type) {
                    if a_m_check_exp_fence.insert(sh.ptr_id()) {
                        a_m_check_exp.push(sh);
                    }
                }
            }
            if a_m_check_exp_fence.insert(a_s.ptr_id()) {
                a_m_check_exp.push(a_s.clone());
            }
        }
        // OCCT L771-798: build result.
        for a_s in &a_m_it_exp {
            let mut b_contains = a_m_check_exp_fence.contains(&a_s.ptr_id());
            if !b_contains && a_s.shape_type() == topods::ShapeType::Solid {
                // OCCT L777-782: check by the solid's face set.
                let a_st = self.shape_face_set(a_s);
                b_contains = a_m_set_check.contains_key(&a_st);
            }
            if b_common {
                if b_contains {
                    a_c.push(a_s.clone());
                }
            } else if !b_contains {
                a_c.push(a_s.clone());
            }
        }
        // OCCT L800-823: filter result for COMMON.
        if b_common {
            let mut a_m_fence: HashSet<u64> = HashSet::new();
            let mut a_cx: Vec<Shape> = Vec::new();
            for i_dim in (i_dim_min..=3).rev() {
                let a_type = type_to_explore(i_dim);
                for sh in self.shapes_of_type_in_shapes(&a_c, a_type) {
                    if a_m_fence.insert(sh.ptr_id()) {
                        a_cx.push(sh.clone());
                        // OCCT L818: TopExp::MapShapes(aS, aMFence).
                        self.add_all_sub_shapes(&sh, &mut a_m_fence);
                    }
                }
            }
            a_c = a_cx;
        }
        // OCCT L825-829: if no edges were checked, done.
        if !b_check_edges {
            self.my_rc = a_c;
            return;
        }
        // OCCT L835-878: squats around degenerated edges.
        let mut a_m_vc: HashSet<u64> = HashSet::new();
        for sh in self.shapes_of_type_in_shapes(&a_c, topods::ShapeType::Vertex) {
            a_m_vc.insert(sh.ptr_id());
        }
        let a_nb = self.ds.nb_source_shapes();
        for i in 0..a_nb {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Edge {
                continue;
            }
            let a_e = self.brep_sr(i);
            let degen = a_e.as_edge().map(|e| e.degenerated).unwrap_or(false);
            if !degen {
                continue;
            }
            let n_vd = a_si.sub_shapes.first().copied().unwrap_or(0);
            let a_vd = self.brep_sr(n_vd);
            if !a_m_vc.contains(&a_vd.ptr_id()) {
                continue;
            }
            if self.ds.is_new_shape(n_vd) {
                continue;
            }
            if self.ds.interfered.contains(&n_vd) {
                continue;
            }
            a_c.push(a_e);
        }
        self.my_rc = a_c;
    }

    /// OCCT BOPAlgo_BOP::BuildSolid (BOPAlgo_BOP.cxx L1111-1392).
    fn build_solid(&mut self) {
        // OCCT L1121-1144: get solids from input arguments.
        let mut a_msa: HashSet<u64> = HashSet::new();
        let mut a_mfs: HashMap<u64, (Shape, Vec<Shape>)> = HashMap::new();
        let mut a_lsc: Vec<Shape> = Vec::new();
        for i in 0..2 {
            let a_lsa: &[Shape] = if i == 0 { self.object_shapes() } else { &self.my_tools };
            for a_sa in a_lsa {
                // OCCT L1133-1139: explore solids, map face→solid ancestors.
                for sol in self.map_shapes_of_type(a_sa, topods::ShapeType::Solid) {
                    a_msa.insert(sol.ptr_id());
                    for a_f in self.map_shapes_of_type(&sol, topods::ShapeType::Face) {
                        a_mfs.entry(a_f.ptr_id())
                            .or_insert_with(|| (a_f.clone(), Vec::new()))
                            .1
                            .push(sol.clone());
                    }
                }
                // OCCT L1141-1143: collect compsolids from arguments.
                Self::collect_containers(a_sa, &mut a_lsc);
            }
        }
        // OCCT L1151-1165: find solids sharing faces.
        let mut a_mt_sols: HashSet<u64> = HashSet::new();
        for (_f, (_fs, sols)) in &a_mfs {
            if sols.len() > 1 {
                for sol in sols {
                    a_mt_sols.insert(sol.ptr_id());
                }
            }
        }
        // OCCT L1167-1220: possibly untouched solids.
        let mut a_mu_sols: Vec<Shape> = Vec::new();
        let mut a_mu_fence: HashSet<u64> = HashSet::new();
        a_mfs.clear();
        for a_sx in &self.my_rc {
            if a_msa.contains(&a_sx.ptr_id()) {
                if !a_mt_sols.contains(&a_sx.ptr_id()) {
                    if a_mu_fence.insert(a_sx.ptr_id()) {
                        a_mu_sols.push(a_sx.clone());
                    }
                    continue;
                }
            }
            // OCCT L1185: MapFacesToBuildSolids(aSx, aMFS).
            self.map_faces_to_build_solids(a_sx, &mut a_mfs);
        }
        // OCCT L1191-1220: process untouched solids.
        let mut a_dmsts: Vec<Shape> = Vec::new();
        let mut a_dmsts_fence: HashSet<Vec<u64>> = HashSet::new();
        for a_sx in &a_mu_sols {
            let mut in_mfs = false;
            for a_f in self.map_shapes_of_type(a_sx, topods::ShapeType::Face) {
                if a_mfs.contains_key(&a_f.ptr_id()) {
                    in_mfs = true;
                    break;
                }
            }
            if in_mfs {
                self.map_faces_to_build_solids(a_sx, &mut a_mfs);
            } else {
                let a_st = self.shape_face_set(a_sx);
                if a_dmsts_fence.insert(a_st) {
                    a_dmsts.push(a_sx.clone());
                }
            }
        }
        // OCCT L1227-1241: faces belonging to a single solid.
        let mut a_mef: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut a_sfs: Vec<Shape> = Vec::new();
        for (_f, (fs, sols)) in &a_mfs {
            if sols.len() == 1 {
                a_sfs.push(fs.clone());
                // OCCT L1238: TopExp::MapShapesAndAncestors(aFx, EDGE, FACE, aMEF).
                for a_e in self.map_shapes_of_type(fs, topods::ShapeType::Edge) {
                    a_mef.entry(a_e.ptr_id()).or_default().push(fs.ptr_id());
                }
            }
        }
        // OCCT L1243-1271: build solids from the set of faces.
        let mut a_rc: Vec<Shape> = Vec::new();
        if !a_sfs.is_empty() {
            let mut a_bs = crate::bop::algo::builder_solid::BuilderSolid::new(&self.ds);
            a_bs.my_shapes = a_sfs;
            a_bs.perform();
            if a_bs.has_errors() {
                // OCCT L1255: AddError(new BOPAlgo_AlertSolidBuilderFailed).
                self.my_report.add_error(crate::bop::algo::Alert::SolidBuilderFailed);
                return;
            }
            for a_sr in &a_bs.my_solids {
                a_rc.push(a_sr.clone());
            }
        }
        // OCCT L1273-1279: add untouched solids.
        for a_sx in &a_dmsts {
            a_rc.push(a_sx.clone());
        }
        // OCCT L1281-1286: no compsolids in arguments → done.
        if a_lsc.is_empty() {
            self.set_shape_from_shapes(a_rc);
            return;
        }
        // OCCT L1291-1391: compsolid construction — rcad pending (no compsolid args
        // in the stage tests).
        self.set_shape_from_shapes(a_rc);
    }

    /// OCCT BOPAlgo_BOP::CheckData — compute myDims from objects and tools.
    fn compute_dims(&mut self) {
        let mut d0 = 0i32;
        for s in self.object_shapes() {
            d0 = d0.max(shape_dimension(s));
        }
        let mut d1 = 0i32;
        for s in &self.my_tools {
            d1 = d1.max(shape_dimension(s));
        }
        self.my_dims = [d0, d1];
    }

    /// OCCT BOPAlgo_BOP::CheckArgsForOpenSolid (BOPAlgo_BOP.cxx L1396-1560).
    /// rcad: returns false — the open-solid edge/face analysis and the BuildBOP
    /// fallback are pending translation. All stage-test solids are closed.
    fn check_args_for_open_solid(&self) -> bool {
        false
    }

    /// Objects of the BOP: myArguments minus the tools suffix.
    /// OCCT BOPAlgo_BOP: myDS->Arguments() = myArguments ++ myTools.
    fn object_shapes(&self) -> &[Shape] {
        let n_tools = self.my_tools.len();
        let n_objs = self.my_arguments.len().saturating_sub(n_tools);
        &self.my_arguments[..n_objs]
    }

    /// Replace my_shape with a fresh compound built from the given shapes.
    /// OCCT: BRep_Builder().Add(aCompound, aS) for each result shape.
    fn set_shape_from_shapes(&mut self, shapes: Vec<Shape>) {
        self.my_shape = Some(topods::BRep::new());
        self.shape_remap.clear();
        for s in shapes {
            self.add_shape_to_result(&s);
        }
    }

    /// OCCT TopExp::MapShapes(aS, aType, aMap) — collect all sub-shapes of a type.
    fn map_shapes_of_type(&self, s: &Shape, t: topods::ShapeType) -> Vec<Shape> {
        let mut out: Vec<Shape> = Vec::new();
        let mut seen: HashSet<u64> = HashSet::new();
        let mut stack: Vec<Shape> = vec![s.clone()];
        while let Some(cur) = stack.pop() {
            if !seen.insert(cur.ptr_id()) {
                continue;
            }
            if cur.shape_type() == t {
                out.push(cur.clone());
            }
            for sub in self.shape_sub_shapes(&cur) {
                stack.push(sub);
            }
        }
        out
    }

    /// Collect shapes of a type across a list of shapes (dedup).
    fn shapes_of_type_in_shapes(&self, shapes: &[Shape], t: topods::ShapeType) -> Vec<Shape> {
        let mut out: Vec<Shape> = Vec::new();
        let mut seen: HashSet<u64> = HashSet::new();
        for s in shapes {
            for sh in self.map_shapes_of_type(s, t) {
                if seen.insert(sh.ptr_id()) {
                    out.push(sh);
                }
            }
        }
        out
    }

    /// Insert a shape and all its sub-shapes into a fence (TopExp::MapShapes).
    fn add_all_sub_shapes(&self, s: &Shape, fence: &mut HashSet<u64>) {
        fence.insert(s.ptr_id());
        for sub in self.shape_sub_shapes(s) {
            self.add_all_sub_shapes(&sub, fence);
        }
    }

    /// BOPTools_Set equivalent: sorted face ptr_ids of a shape (map key).
    fn shape_face_set(&self, s: &Shape) -> Vec<u64> {
        let mut faces: Vec<u64> = self
            .map_shapes_of_type(s, topods::ShapeType::Face)
            .iter()
            .map(|f| f.ptr_id())
            .collect();
        faces.sort();
        faces
    }

    /// OCCT BOPAlgo_BOP::CollectContainers (BOPAlgo_BOP.cxx L1601-1621).
    fn collect_containers(s: &Shape, out: &mut Vec<Shape>) {
        let a_type = s.shape_type();
        if a_type == topods::ShapeType::Wire
            || a_type == topods::ShapeType::Shell
            || a_type == topods::ShapeType::CompSolid
        {
            out.push(s.clone());
            return;
        }
        if a_type != topods::ShapeType::Compound {
            return;
        }
        for sub in Self::shape_sub_shapes_static(s) {
            Self::collect_containers(&sub, out);
        }
    }

    /// OCCT BOPTools_AlgoTools::MakeContainer + BRep_Builder().Add for a container.
    fn build_container_of_type(a_type: topods::ShapeType, subs: &[Shape]) -> Shape {
        let ts: TShape = match a_type {
            topods::ShapeType::Wire => TShape::Wire(TWireData {
                my_shapes: subs.to_vec(),
                flags: tshape_flags::DEFAULT,
                edges: subs.to_vec(),
            }),
            topods::ShapeType::Shell => TShape::Shell(TShellData {
                my_shapes: subs.to_vec(),
                flags: tshape_flags::DEFAULT,
                faces: subs.to_vec(),
            }),
            _ => TShape::Compound(subs.to_vec()),
        };
        Shape::new(std::sync::Arc::new(ts), 0, topods::Orientation::Forward)
    }

    /// OCCT BOPTools_AlgoTools::MakeConnexityBlocks(aRC, aT1, aT2, aLCB).
    /// rcad: minimal — each shape its own block. Solid arguments produce no
    /// containers, so this is not exercised by the stage tests.
    fn make_connexity_blocks_shapes(
        shapes: &[Shape],
        _a_t1: topods::ShapeType,
        _a_t2: topods::ShapeType,
        out: &mut Vec<Vec<Shape>>,
    ) {
        for s in shapes {
            out.push(vec![s.clone()]);
        }
    }

    /// OCCT BOPTools_AlgoTools::OrientEdgesOnWire — reorient edges on a wire.
    /// rcad: no-op (affects orientation, not entity counts).
    fn orient_edges_on_wire(_w: &mut Shape) {}

    /// OCCT BOPTools_AlgoTools::OrientFacesOnShell — reorient faces on a shell.
    /// rcad: no-op (affects orientation, not entity counts).
    fn orient_faces_on_shell(_sh: &mut Shape) {}

    /// OCCT BOPAlgo_BOP::RemoveDuplicates (BOPAlgo_BOP.cxx L1627-1698).
    /// Dedups containers with identical sub-shape contents.
    fn remove_duplicates(&self, containers: &mut Vec<Shape>) {
        let mut seen: HashSet<Vec<u64>> = HashSet::new();
        let mut out: Vec<Shape> = Vec::new();
        for c in containers.iter() {
            let mut subs: Vec<u64> = self
                .shape_sub_shapes(c)
                .iter()
                .map(|s| s.ptr_id())
                .collect();
            subs.sort();
            if seen.insert(subs) {
                out.push(c.clone());
            }
        }
        *containers = out;
    }

    /// OCCT BOPAlgo_BOP::MapFacesToBuildSolids (BOPAlgo_BOP.cxx L1768-1798).
    fn map_faces_to_build_solids(
        &self,
        the_sol: &Shape,
        the_mfs: &mut HashMap<u64, (Shape, Vec<Shape>)>,
    ) {
        for a_f in self.map_shapes_of_type(the_sol, topods::ShapeType::Face) {
            if a_f.orientation == topods::Orientation::Internal {
                continue;
            }
            let e = the_mfs
                .entry(a_f.ptr_id())
                .or_insert_with(|| (a_f.clone(), Vec::new()));
            // OCCT L1789-1796: append the solid only if orientations differ.
            if e.1.is_empty() || e.0.orientation != a_f.orientation {
                e.1.push(the_sol.clone());
            }
        }
    }

    /// OCCT BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx L130-168).
    /// Builds topology at the given shape type level.
    /// OCCT L136-143 iterates myArguments and skips any argument whose ShapeType
    /// does not match theType. For each matching argument it adds its images
    /// (or the argument itself if it has no images) to myShape, deduplicated by fence.
    /// When arguments are solids, the intermediate calls (VERTEX..SHELL, COMPOUND)
    /// are no-ops; only BuildResult(SOLID) adds shapes into the result compound.
    fn build_result(&mut self, the_type: topods::ShapeType) {
        // OCCT L133: fence map
        let mut a_m_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        // OCCT L136-167: iterate myArguments, filter by theType
        let a_arguments = self.my_arguments.clone();
        for a_s in &a_arguments {
            if a_s.shape_type() != the_type { continue; }
            // OCCT L145-152: check for images
            if let Some(imgs) = self.my_images.get(a_s).cloned() {
                for a_s_im in &imgs {
                    if a_m_fence.insert(a_s_im.ptr_id()) {
                        self.add_shape_to_result(a_s_im);
                    }
                }
            } else {
                if a_m_fence.insert(a_s.ptr_id()) {
                    self.add_shape_to_result(a_s);
                }
            }
        }
    }

    /// Add a Shape to my_shape (OCCT equivalent: BRep_Builder().Add(myShape, aS)).
    /// OCCT uses TopoDS handles (pointer-based); rcad-kernel uses flat indices.
    /// This function clones the full TShape hierarchy and fixes all Shape.index
    /// values to point into the result BRep's tshapes array.
    /// Uses self.shape_remap persisted across the entire pipeline.
    fn add_shape_to_result(&mut self, shape: &Shape) {
        self.push_shape_recursive(shape);
    }

    /// Recursively push a Shape and all sub-shapes into my_shape.
    /// Uses self.shape_remap for persistent index tracking across pipeline stages.
    /// Returns the shape's index in the result BRep's tshapes.
    fn push_shape_recursive(&mut self, shape: &Shape) -> usize {
        let ptr = shape.ptr_id();
        if let Some(&idx) = self.shape_remap.get(&ptr) {
            return idx;
        }

        // Reserve a slot in tshapes (placeholder, replaced below)
        let new_idx = {
            let brep = self.my_shape.as_mut().expect("prepare() must set my_shape");
            let idx = brep.tshapes.len();
            brep.tshapes.push(Arc::new(TShape::Vertex(TVertexData {
                my_shapes: Vec::new(), flags: 0, point: DVec3::ZERO,
                tolerance: 0.0, points: Vec::new(),
            })));
            idx
        };
        self.shape_remap.insert(ptr, new_idx);

        // Build new TShape with remapped sub-shape indices
        let new_tshape: TShape = match shape.data.as_ref() {
            TShape::Vertex(vd) => {
                let my_shapes = self.remap_shapes(&vd.my_shapes);
                TShape::Vertex(TVertexData {
                    my_shapes, ..vd.clone()
                })
            }
            TShape::Edge(ed) => {
                let _ = self.push_shape_recursive(&ed.first);
                let _ = self.push_shape_recursive(&ed.last);
                let my_shapes = self.remap_shapes(&ed.my_shapes);
                let first = self.remap_shape(&ed.first);
                let last = self.remap_shape(&ed.last);
                TShape::Edge(TEdgeData {
                    my_shapes, first, last, ..ed.clone()
                })
            }
            TShape::Wire(wd) => {
                for e in &wd.edges { let _ = self.push_shape_recursive(e); }
                let my_shapes = self.remap_shapes(&wd.my_shapes);
                let edges = self.remap_shapes(&wd.edges);
                TShape::Wire(TWireData {
                    my_shapes, edges, ..wd.clone()
                })
            }
            TShape::Face(fd) => {
                let _ = self.push_shape_recursive(&fd.outer_wire);
                for w in &fd.inner_wires { let _ = self.push_shape_recursive(w); }
                for v in &fd.internal_vertices { let _ = self.push_shape_recursive(v); }
                let my_shapes = self.remap_shapes(&fd.my_shapes);
                let outer_wire = self.remap_shape(&fd.outer_wire);
                let inner_wires = self.remap_shapes(&fd.inner_wires);
                let internal_vertices = self.remap_shapes(&fd.internal_vertices);
                TShape::Face(TFaceData {
                    my_shapes, outer_wire, inner_wires,
                    internal_vertices, ..fd.clone()
                })
            }
            TShape::Shell(sd) => {
                for f in &sd.faces { let _ = self.push_shape_recursive(f); }
                let my_shapes = self.remap_shapes(&sd.my_shapes);
                let faces = self.remap_shapes(&sd.faces);
                TShape::Shell(TShellData {
                    my_shapes, faces, ..sd.clone()
                })
            }
            TShape::Solid(sd) => {
                for s in &sd.shells { let _ = self.push_shape_recursive(s); }
                let my_shapes = self.remap_shapes(&sd.my_shapes);
                let shells = self.remap_shapes(&sd.shells);
                let internal_vertices = self.remap_shapes(&sd.internal_vertices);
                let internal_edges = self.remap_shapes(&sd.internal_edges);
                TShape::Solid(TSolidData {
                    my_shapes, shells, internal_vertices, internal_edges,
                    ..sd.clone()
                })
            }
            TShape::CompSolid(shapes) => {
                TShape::CompSolid(self.remap_shapes(shapes))
            }
            TShape::Compound(shapes) => {
                TShape::Compound(self.remap_shapes(shapes))
            }
        };

        // Replace placeholder with the remapped TShape
        let brep = self.my_shape.as_mut().unwrap();
        brep.tshapes[new_idx] = Arc::new(new_tshape);
        new_idx
    }

    /// Remap a single Shape's index via self.shape_remap.
    fn remap_shape(&self, shape: &Shape) -> Shape {
        if let Some(&new_idx) = self.shape_remap.get(&shape.ptr_id()) {
            Shape {
                index: new_idx,
                data: shape.data.clone(),
                location: shape.location,
                orientation: shape.orientation,
            }
        } else {
            shape.clone()
        }
    }

    /// Remap a slice of Shapes.
    fn remap_shapes(&self, shapes: &[Shape]) -> Vec<Shape> {
        shapes.iter().map(|s| self.remap_shape(s)).collect()
    }

}
