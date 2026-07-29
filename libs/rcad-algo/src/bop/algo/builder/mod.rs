// OCCT BOPAlgo_Builder — shape construction from DS.
//
// OCCT BOPAlgo_Builder.hxx L75-507 + parent class fields (BOPAlgo_BuilderShape, BOPAlgo_Options, BOPAlgo_BOP).
// Flattened into one Rust struct because Rust has no C++ inheritance.

pub use crate::bop::algo::BooleanOpType;
use crate::bop::algo::{GlueEnum, Report};
use crate::bop::algo::builder_face::BuilderFace;
use crate::bop::ds::DS;
use rcad_kernel::topods;
use rcad_kernel::topods::{
    TShape, TVertexData, TEdgeData, TWireData, TFaceData,
    TShellData, TSolidData, tshape_flags,
};
use rcad_kernel::topo_shape::Shape;
use std::collections::HashMap;
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
    // BOPAlgo_Builder — primary data source
    pub(crate) ds: &'a DS,

    // BOPAlgo_Options
    pub(crate) my_report: Report,
    pub(crate) my_run_parallel: bool,
    pub(crate) my_fuzzy_value: f64,

    // BOPAlgo_BuilderShape
    pub(crate) my_shape: Option<topods::BRep>,
    pub(crate) my_fill_history: bool,

    // BOPAlgo_BOP
    pub(crate) my_operation: BooleanOpType,

    // BOPAlgo_Builder — shape images and SD
    // OCCT: myImages (DataMap<Shape, List<Shape>>)
    pub(crate) my_images: HashMap<Shape, Vec<Shape>>,
    // OCCT: myShapesSD (DataMap<Shape, Shape>) — same-defined pairs
    pub(crate) my_shapes_sd: HashMap<Shape, Shape>,
    // OCCT: myOrigins (DataMap<Shape, List<Shape>>) — reverse map of images
    pub(crate) my_origins: HashMap<Shape, Vec<Shape>>,
    // OCCT: myInParts (DataMap<Shape, List<Shape>>) — IN faces of argument solids
    pub(crate) my_in_parts: HashMap<usize, Vec<usize>>,

    // BOPAlgo_Builder — glue mode
    pub(crate) my_glue: GlueEnum,
    // Persistent remap for BRep index tracking (ptr_id -> tshapes index)
    pub(crate) shape_remap: HashMap<u64, usize>,
}

impl<'a> BooleanBuilder<'a> {
    /// Create a new BooleanBuilder referencing a DS.
    ///
    /// OCCT: BOPAlgo_Builder is constructed with a PaveFiller reference.
    /// rcad: PaveFiller runs before Builder; only the DS is passed.
    pub fn new(ds: &'a DS, op: BooleanOpType) -> Self {
        BooleanBuilder {
            ds,
            my_report: Report::new(),
            my_run_parallel: false,
            my_fuzzy_value: 1e-7,
            my_shape: None,
            my_fill_history: false,
            my_operation: op,
            my_images: HashMap::new(),
            my_shapes_sd: HashMap::new(),
            my_origins: HashMap::new(),
            my_in_parts: HashMap::new(),
            my_glue: GlueEnum::GlueOff,
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

    /// OCCT BOPAlgo_Builder::CheckData (BOPAlgo_Builder_3.cxx L61-93).
    fn check_data(&self) -> Result<(), BooleanError> {
        Ok(())
    }

    /// OCCT BOPAlgo_BOP::CheckFiller (BOPAlgo_BOP.cxx L138).
    fn check_filler(&self) {
        // OCCT L138: verifies PaveFiller was executed
        // rcad: PaveFiller always runs before Builder
    }

    /// OCCT BOPAlgo_Builder::Prepare (BOPAlgo_Builder_3.cxx L95-165).
    /// Returns a ResultBuilder (rcad-specific, no OCCT equivalent).
    fn prepare(&mut self) {
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
            // OCCT L281-287: get FaceInfo — PaveBlocksIn/On/Sc, AloneVertices
            let a_fi = self.ds.face_info(i).clone();
            let a_nb_pb_in = a_fi.pave_blocks_in.len();
            let a_nb_pb_on = a_fi.pave_blocks_on.len();
            let a_nb_pb_sc = a_fi.pave_blocks_sc.len();
            let a_nb_av = a_fi.vertices_in.len();
            // OCCT L293-296: skip if no IN/ON/SC edges and no alone vertices
            if a_nb_pb_in == 0 && a_nb_pb_on == 0 && a_nb_pb_sc == 0 && a_nb_av == 0 {
                continue;
            }
            // OCCT L298-310: if no IN and no SC edges → check wire modifications
            if a_nb_pb_in == 0 && a_nb_pb_sc == 0 {
                let face_s = self.brep_sr(i);
                let sub_shapes = self.shape_sub_shapes(&face_s);
                let has_modified_wires = sub_shapes.iter().any(|ss| {
                    self.my_images.contains_key(ss)
                });
                if !has_modified_wires && a_nb_av == 0 {
                    continue;
                }
            }
            // OCCT: Create BuilderFace, process face
            let face_s = self.brep_sr(i);
            let sub_edges = self.shape_sub_shapes(&face_s);
            let mut bf = BuilderFace::new(self.ds);
            bf.my_face = Some(face_s.clone());
            bf.my_edges = sub_edges;
            bf.perform();
            if !bf.my_images.is_empty() {
                for img in &bf.my_images {
                    self.my_images.entry(face_s.clone())
                        .or_default()
                        .push(img.clone());
                }
            }
        }
    }

    /// OCCT BOPAlgo_Builder::FillSameDomainFaces (Builder_2.cxx L580-780).
    fn fill_same_domain_faces(&mut self) {
        let a_ffs = &self.ds.interf_ff;
        if a_ffs.is_empty() { return; }
        let mut to_remove: Vec<Shape> = Vec::new();
        for ff in a_ffs {
            let f1s = self.brep_sr(ff.f1);
            let f2s = self.brep_sr(ff.f2);
            if !self.my_images.contains_key(&f1s) { continue; }
            if !self.my_images.contains_key(&f2s) { continue; }
            if ff.tangent_faces {
                to_remove.push(f2s);
            }
        }
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
                    if self.my_images.contains_key(&self.brep_sr(fi)) {
                        self.my_in_parts.entry(si).or_default().push(fi);
                    }
                }
            }
        }
        eprintln!("[F3P] in_parts={}", self.my_in_parts.len());
    }

    /// OCCT BOPAlgo_Builder::BuildSplitSolids (Builder_3.cxx L400-550).
    fn build_split_solids(&mut self) {
        for (&solid_idx, face_indices) in &self.my_in_parts {
            let face_shapes: Vec<Shape> = face_indices.iter()
                .filter_map(|&fi| {
                    if fi < self.ds.nb_shapes()
                        && self.ds.shape_info(fi).shape_type == topods::ShapeType::Face {
                        Some(self.brep_sr(fi))
                    } else { None }
                })
                .collect();
            if face_shapes.is_empty() { continue; }
            let mut bs = crate::bop::algo::builder_solid::BuilderSolid::new(self.ds);
            bs.my_shapes = face_shapes;
            bs.perform();
            let solid_src = self.brep_sr(solid_idx);
            for solid_img in &bs.my_solids {
                self.my_images.entry(solid_src.clone())
                    .or_default().push(solid_img.clone());
            }
        }
    }
    fn fill_internal_shapes(&mut self) {
        // OCCT: adds internal vertices/edges. Stub.
    }

    /// OCCT BOPAlgo_Builder::FillImagesCompounds (BOPAlgo_Builder_1.cxx L197-217).
    fn fill_images_compounds(&mut self) {
        // OCCT: for each source compound, FillImagesCompound. Stub.
    }

    /// OCCT BOPAlgo_Builder::PrepareHistory (BOPAlgo_Builder_4.cxx L164-252).
    fn prepare_history(&mut self) {
        if !self.my_fill_history { return; }
        // OCCT: for each source shape, check images → Modified/Generated/Deleted
    }

    /// OCCT BOPAlgo_Builder::PostTreat (BOPAlgo_Builder.cxx L461+).
    fn post_treat(&mut self) {
        // OCCT: non-destructive mode. Stub.
    }

    /// OCCT BOPAlgo_Builder::BuildResult (BOPAlgo_Builder_1.cxx L130-168).
    /// Builds topology at the given shape type level.
    /// For each argument shape of `the_type` (found via TopExp_Explorer equivalent:
    /// iterate all DS shapes), adds its images (or itself if no images)
    /// to my_shape, deduplicated by fence.
    fn build_result(&mut self, the_type: topods::ShapeType) {
        // OCCT L132-133: fence map for dedup
        let mut a_m_fence: std::collections::HashSet<u64> = std::collections::HashSet::new();
        // OCCT L136-167: TopExp_Explorer over arguments — iterate all DS shapes
        let a_nb_s = self.ds.nb_shapes();
        let mut n_added = 0usize;
        for i in 0..a_nb_s {
            // OCCT L140-142: filter by shape type (TopExp_Explorer does this)
            let a_si = self.ds.shape_info(i);
            if a_si.shape_type != the_type {
                continue;
            }
            // OCCT L145: aS = argument Shape
            let a_s = self.brep_sr(i);
            // OCCT L145-152: check for images
            let p_ls_im = self.my_images.get(&a_s).cloned();
            if let Some(imgs) = p_ls_im {
                // OCCT L156-164: add images
                for a_s_im in &imgs {
                    if a_m_fence.insert(a_s_im.ptr_id()) {
                        self.add_shape_to_result(&a_s_im);
                        n_added += 1;
                    }
                }
            } else {
                // OCCT L148-151: no images — add self
                if a_m_fence.insert(a_s.ptr_id()) {
                    self.add_shape_to_result(&a_s);
                    n_added += 1;
                }
            }
        }
        eprintln!("[BUILDER] build_result({:?}): {} added (tshapes={})",
            the_type, n_added,
            self.my_shape.as_ref().map_or(0, |b| b.tshapes.len()));
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
