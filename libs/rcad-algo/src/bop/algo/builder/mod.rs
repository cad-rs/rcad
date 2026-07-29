// OCCT BOPAlgo_Builder — shape construction from DS.
//
// OCCT BOPAlgo_Builder.hxx L75-507 + parent class fields (BOPAlgo_BuilderShape, BOPAlgo_Options, BOPAlgo_BOP).
// Flattened into one Rust struct because Rust has no C++ inheritance.

pub use crate::bop::algo::BooleanOpType;
use crate::bop::algo::{GlueEnum, Report};
use crate::bop::algo::builder_face::BuilderFace;
use crate::bop::ds::DS;
use crate::bop::int_tools::context::IntToolsContext;
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
pub struct BooleanBuilder<'a> {
    // ── BOPAlgo_Options (inherited) ─────────────────────────────
    pub(crate) my_report: Report,          // BOPAlgo_Algo::myReport
    pub(crate) my_run_parallel: bool,      // BOPAlgo_Algo::myRunParallel
    pub(crate) my_fuzzy_value: f64,        // BOPAlgo_Algo::myFuzzyValue
    // ── BOPAlgo_BuilderShape (inherited) ───────────────────────
    pub(crate) my_shape: Option<topods::BRep>, // BOPAlgo_BuilderShape::myShape
    pub(crate) my_fill_history: bool,      // BOPAlgo_BuilderShape::myFillHistory
    // ── BOPAlgo_BOP (inherited) ────────────────────────────────
    pub(crate) my_operation: BooleanOpType, // BOPAlgo_BOP::myOperation
    // ── BOPAlgo_Builder.hxx L492-505 ───────────────────────────
    pub(crate) ds: &'a DS,                 // L496: myDS
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

impl<'a> BooleanBuilder<'a> {
    /// Create a new BooleanBuilder referencing a DS.
    ///
    /// OCCT: BOPAlgo_Builder is constructed with a PaveFiller reference.
    /// rcad: PaveFiller runs before Builder; only the DS is passed.
    /// OCCT BOPAlgo_BOP::PerformInternal1 L425-429:
    ///   myPaveFiller = &theFiller; myDS = myPaveFiller->PDS();
    ///   myFuzzyValue = myPaveFiller->FuzzyValue();
    /// rcad: PaveFiller dropped before Builder, DS+fuzzy passed explicitly.
    pub fn new(ds: &'a DS, op: BooleanOpType, fuzzy_value: f64) -> Self {
        BooleanBuilder {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_fuzzy_value: fuzzy_value,
            my_shape: None,
            my_fill_history: false,
            my_operation: op,
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

    /// OCCT BOPAlgo_Builder::Build — full pipeline with history.
    ///
    /// Mirrors OCCT BOPAlgo_Builder::PerformInternal1 (BOPAlgo_Builder_3.cxx).
    pub fn build_with_history_topods(
        &mut self,
    ) -> Result<(topods::BRep, ()), BooleanError> {
        // OCCT L425-429: setup from PaveFiller (rcad: done in constructor)

        // OCCT L431-436: CheckData
        self.check_data()?;
        self.check_filler();

        // OCCT L438-443: Prepare
        let _result = self.prepare();

        eprintln!("[BUILDER] DS: {} shapes ({} src, {}V/{}E/{}F)",
            self.ds.nb_shapes(), self.ds.nb_source_shapes(),
            self.ds.vertex_count(), self.ds.edge_count(), self.ds.face_count());
        let (nsh, nso) = {
            let mut a = 0usize; let mut b = 0usize;
            for i in 0..self.ds.nb_shapes() {
                match self.ds.shape_info(i).shape_type {
                    topods::ShapeType::Shell => a += 1,
                    topods::ShapeType::Solid => b += 1,
                    _ => {}
                }
            } (a, b)
        };
        eprintln!("[BUILDER] DS: {} Shell, {} Solid", nsh, nso);

        // OCCT L445-453: TreatEmptyShape
        // (skipped — rcad handles in boolean_op_with_retry)

        // OCCT L459-471: FillImagesVertices + BuildResult(VERTEX)
        self.fill_images_vertices();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Vertex);

        // OCCT L472-483: FillImagesEdges + BuildResult(EDGE)
        self.fill_images_edges();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Edge);

        // OCCT L484-494: FillImagesContainers(WIRE) + BuildResult(WIRE)
        self.fill_images_containers(topods::ShapeType::Wire);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Wire);

        // OCCT L496-505: FillImagesFaces + BuildResult(FACE)
        self.fill_images_faces();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Face);

        // OCCT L507-516: FillImagesContainers(SHELL) + BuildResult(SHELL)
        self.fill_images_containers(topods::ShapeType::Shell);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Shell);

        // OCCT L518-528: FillImagesSolids + BuildResult(SOLID)
        self.fill_images_solids();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Solid);

        // OCCT L530-539: FillImagesContainers(COMPSOLID) + BuildResult(COMPSOLID)
        self.fill_images_containers(topods::ShapeType::CompSolid);
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::CompSolid);

        // OCCT L541-550: FillImagesCompounds + BuildResult(COMPOUND)
        self.fill_images_compounds();
        if self.has_errors() {
            return Err(BooleanError::DegenerateResult);
        }
        self.build_result(topods::ShapeType::Compound);

        // OCCT L552-570: PrepareHistory + PostTreat
        self.prepare_history();
        self.post_treat();

        let result = self.my_shape.clone().unwrap_or_default();
        Ok((result, ()))
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
    /// Maps source edges -> split images via pave-block new_edge.
    /// Also handles CommonBlocks via myShapesSD.
    fn fill_images_edges(&mut self) {
        let aNbS = self.ds.nb_source_shapes();
        // Collect (edge_index, pave_block_indices) for all edges first
        let mut edge_splits: Vec<(usize, Vec<usize>)> = Vec::new();
        for i in 0..aNbS {
            let aSI = self.ds.shape_info(i);
            if aSI.shape_type != topods::ShapeType::Edge {
                continue;
            }
            if !aSI.has_reference() {
                continue;
            }
            let aLPB = self.ds.pave_blocks(i);
            let splits: Vec<usize> = aLPB.iter().map(|pb| pb.0.read().unwrap().edge).collect();
            edge_splits.push((i, splits));
        }
        // Process collected data
        for (i, splits) in edge_splits {
            let aE = self.brep_sr(i);
            for &nSpR in &splits {
                let aSpR = self.brep_sr(nSpR);
                // myImages[aE] -> append aSpR  (OCCT L95, L105)
                self.my_images.entry(aE.clone()).or_default().push(aSpR.clone());
                // myOrigins[aSpR] -> append aE  (OCCT L107-112)
                self.my_origins.entry(aSpR.clone()).or_default().push(aE.clone());
                // CommonBlock handling  (OCCT L114-119)
                let pb = &self.ds.pave_blocks(i);
                for apb in pb.iter() {
                    if self.ds.is_common_block_on_edge(apb) {
                        let pbb = apb.0.read().unwrap();
                        let aSp = self.brep_sr(pbb.edge);
                        self.my_shapes_sd.insert(aSp, aSpR.clone());
                    }
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
        eprintln!("[FIC] {:?} subs={} mod={} imgs={}", the_type, sub_shapes.len(), has_modified, self.my_images.len());
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
                        // OCCT L265-269: IsSplitToReverseWithWarn check (skipped)
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

    /// OCCT BOPAlgo_Builder::BuildSplitFaces (BOPAlgo_Builder_2.cxx L233-380).
    /// For each source face, builds a BuilderFace with section edges.
    fn build_split_faces(&mut self) {
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != topods::ShapeType::Face {
                continue;
            }
            // OCCT L275-279: Check HasFaceInfo
            if !self.ds.has_face_info(i) {
                continue;
            }
            // OCCT L275-279: Check HasFaceInfo
            if !self.ds.has_face_info(i) { continue; }
            let a_fi = self.ds.face_info(i).clone();
            let has_sc = !a_fi.pave_blocks_sc.is_empty();
            let has_in = !a_fi.pave_blocks_in.is_empty();
            let has_av = !a_fi.vertices_in.is_empty();
            if !has_sc && !has_in && !has_av { continue; }
            // OCCT: section edges from PaveBlocksSc
            let section_edges: Vec<Shape> = a_fi.pave_blocks_sc.iter()
                .filter_map(|&pb_idx| {
                    if pb_idx >= self.ds.pave_blocks_pool.len() { return None; }
                    let ei = self.ds.pave_blocks_pool[pb_idx].first()
                        .and_then(|pb| {
                            let e = pb.0.read().unwrap().edge;
                            if e < self.ds.nb_shapes() { Some(e) } else { None }
                        })?;
                    Some(self.brep_sr(ei))
                })
                .collect();
            if section_edges.is_empty() { continue; }
            let face_s = self.brep_sr(i);
            let mut bf = BuilderFace::new(self.ds);
            bf.my_face = Some(face_s.clone());
            bf.my_face_index = Some(i);
            bf.my_edges = section_edges;
            bf.perform();
            if !bf.my_areas.is_empty() {
                for img in &bf.my_areas {
                    self.my_images.entry(face_s.clone())
                        .or_default()
                        .push(img.clone());
                }
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillSameDomainFaces (Builder_2.cxx L580-780).
    fn fill_same_domain_faces(&mut self) {
        // OCCT L584-589: get FF interferences, empty check
        let a_ffs = &self.ds.interf_ff;
        if a_ffs.is_empty() { return; }

        // OCCT L597-649: Build face-to-parent solid map
        // Maps each source face to its parent solid, including propagated image faces
        let mut a_face_to_parent: HashMap<u64, u64> = HashMap::new(); // face_ptr_id → solid_ptr_id
        let a_nb_src = self.ds.nb_source_shapes();
        for i_src in 0..a_nb_src {
            let a_si = self.ds.shape_info(i_src);
            if a_si.shape_type != topods::ShapeType::Solid { continue; }
            let a_solid = self.brep_sr(i_src);
            // Iterate solid's sub-shape shells → faces
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

        // OCCT L654-684: collect face indices from FF interferences
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

        // OCCT L687: sort indices
        a_fi_vec.sort();

        // OCCT L690-694: map edge-sets → list of faces, planar face fence
        // rcad: edge_set = sorted Vec of (edge_ptr_id, curve_type) for face boundary edges
        let mut an_e_set_faces: std::collections::HashMap<Vec<(u64, u32)>, Vec<Shape>> =
            std::collections::HashMap::new();
        let mut a_mf_planar: std::collections::HashSet<u64> = std::collections::HashSet::new();

        // OCCT L697-741: for each face, compute edge set
        for &n_f in &a_fi_vec {
            let a_f = self.brep_sr(n_f);

            // OCCT L707-718: check if planar
            let b_check_planar = {
                // Get surface type from face_info
                if let Some(surf) = self.ds.face_surface(n_f) {
                    matches!(surf, rcad_kernel::geom::Surface3::Plane(_))
                } else { false }
            };

            // OCCT L720-740: get face images (or face itself), add edge set
            let face_list: Vec<Shape> = if let Some(imgs) = self.my_images.get(&a_f) {
                imgs.clone()
            } else {
                vec![a_f.clone()]
            };

            for f_piece in &face_list {
                // Build edge set: (edge_ptr_id, edge_type_code) for each boundary edge
                let mut edge_set: Vec<(u64, u32)> = Vec::new();
                let edges = self.ds.face_boundary_edges(n_f);
                for &ei in &edges {
                    if ei >= self.ds.nb_shapes() { continue; }
                    let e_ptr = self.ds.shape(ei).ptr_id();
                    // Get curve type code (0=Line, 1=Circle, 2=Ellipse, 3=Other)
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

        // OCCT L750-780: find faces with same edge set → check SD pairs
        let mut to_remove: Vec<Shape> = Vec::new();
        for (_edge_set, faces) in &an_e_set_faces {
            if faces.len() < 2 { continue; }

            for i1 in 0..faces.len() {
                let f1 = &faces[i1];
                let parent1 = a_face_to_parent.get(&f1.ptr_id()).copied();
                let b_planar1 = a_mf_planar.contains(&f1.ptr_id());

                for i2 in (i1 + 1)..faces.len() {
                    let f2 = &faces[i2];
                    let parent2 = a_face_to_parent.get(&f2.ptr_id()).copied();

                    // OCCT L776-779: skip if both from same parent solid
                    if let (Some(p1), Some(p2)) = (parent1, parent2) {
                        if p1 == p2 { continue; }
                    }

                    // OCCT L782: if planar face → accept as SD
                    if b_planar1 {
                        to_remove.push(f2.clone());
                    }
                }
            }
        }

        // Remove SD faces from images (keep one per group)
        for r in &to_remove {
            self.my_images.remove(r);
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
        // OCCT L62-73: check all DS shapes for SOLID type
        let a_nb_s = self.ds.nb_shapes();
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

    /// OCCT BOPAlgo_Builder::FillIn3DParts (Builder_3.cxx L97+).
    fn fill_in_3d_parts(&mut self) {
        let n = self.ds.nb_shapes();
        for si in 0..n {
            if self.ds.shape_info(si).shape_type != topods::ShapeType::Solid { continue; }
            for &shi in &self.ds.shape_info(si).sub_shapes {
                if shi >= n { continue; }
                if self.ds.shape_info(shi).shape_type != topods::ShapeType::Shell { continue; }
                for &fi in &self.ds.shape_info(shi).sub_shapes {
                    if fi >= n { continue; }
                    if self.ds.shape_info(fi).shape_type != topods::ShapeType::Face { continue; }
                    let f_shape = self.brep_sr(fi);
                    if self.my_images.contains_key(&f_shape) {
                        let s_shape = self.brep_sr(si);
                        self.my_in_parts.entry(s_shape).or_default().push(f_shape);
                    }
                }
            }
        }
        eprintln!("[F3P] in_parts={}", self.my_in_parts.len());
    }

    /// OCCT BOPAlgo_Builder::BuildSplitSolids (Builder_3.cxx L400-550).
    fn build_split_solids(&mut self) {
        let in_parts: Vec<(Shape, Vec<Shape>)> = self.my_in_parts.iter()
            .map(|(k, v)| (k.clone(), v.clone())).collect();
        for (solid_src, face_shapes) in &in_parts {
            if face_shapes.is_empty() { continue; }
            let mut bs = crate::bop::algo::builder_solid::BuilderSolid::new(self.ds);
            bs.my_shapes = face_shapes.clone();
            bs.perform();
            for solid_img in &bs.my_solids {
                self.my_images.entry(solid_src.clone())
                    .or_default().push(solid_img.clone());
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

        // OCCT L820-830: classify and add internal shapes
        // rcad: BRepClass3d_SolidClassifier for V/E classification
        // pending in rcad-algo.
        let _ = a_ls_d;
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

    /// OCCT BOPAlgo_Builder::OwnInternalShapes — get internal V/E of a solid.
    fn own_internal_shapes(&self, s: &Shape) -> Vec<Shape> {
        let mut result = Vec::new();
        for ss in self.shape_sub_shapes(s) {
            if ss.orientation == topods::Orientation::Internal {
                result.push(Shape::new(ss.data.clone(), ss.location, ss.orientation));
            }
            for sub in self.shape_sub_shapes(&ss) {
                if sub.orientation == topods::Orientation::Internal {
                    result.push(Shape::new(sub.data.clone(), sub.location, sub.orientation));
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

    /// OCCT BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx L130-168).
    /// Builds topology at the given shape type level.
    /// For each argument shape of `the_type` (found via TopExp_Explorer equivalent:
    /// iterate all DS shapes), adds its images (or itself if no images)
    /// to my_shape, deduplicated by fence.
    /// OCCT BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx L130-168).
    fn build_result(&mut self, the_type: topods::ShapeType) {
        // OCCT L133: fence map
        let mut a_m_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        // OCCT L136-167: iterate all source shapes of given type
        // rcad: iterate DS source shapes (equivalent to TopExp_Explorer over arguments)
        let a_nb_s = self.ds.nb_source_shapes();
        for i in 0..a_nb_s {
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != the_type { continue; }
            let a_s = self.brep_sr(i);
            // OCCT L145-152: check for images
            if let Some(imgs) = self.my_images.get(&a_s).cloned() {
                for a_s_im in &imgs {
                    if a_m_fence.insert(a_s_im.ptr_id()) {
                        self.add_shape_to_result(a_s_im);
                    }
                }
            } else {
                if a_m_fence.insert(a_s.ptr_id()) {
                    self.add_shape_to_result(&a_s);
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
