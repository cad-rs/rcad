// OCCT BOPAlgo_MakerVolume — builds solids from a set of shapes.
//
// OCCT source: BOPAlgo_MakerVolume.hxx L108-205 + BOPAlgo_MakerVolume.cxx
// L1-414. The algorithm intersects the given shapes (BOPAlgo_Builder image
// pipeline, when myIntersect is TRUE), collects the resulting faces together
// with the faces of a covering bounding-box solid, builds solids from them
// with BOPAlgo_BuilderSolid, removes the covering solid and fills the created
// solids with the internal vertices/edges of the arguments.
//
// OCCT inheritance chain (Rust has no inheritance -> composition + delegation):
//   BOPAlgo_MakerVolume : BOPAlgo_Builder
//     : BOPAlgo_BuilderShape : BOPAlgo_Algo : BOPAlgo_Options
// rcad mapping:
//   - The flattened BOPAlgo_Options / BOPAlgo_BuilderShape members
//     (myReport, myRunParallel, myFuzzyValue, myShape, myFillHistory) live
//     directly in this struct (same pattern as Builder in builder.rs).
//   - BOPAlgo_Builder's shared state used by the inherited pipeline steps
//     (myDS, myContext, myImages, myShapesSD, myOrigins) lives in the
//     `Builder` instance composed inside `perform_internal` (the rcad
//     BOPAlgo_Builder translation); its fields and the inherited pipeline
//     methods (prepare, fill_images_vertices, fill_images_edges,
//     fill_images_containers, fill_images_faces, ...) are pub(crate) and are
//     called / read here exactly where OCCT calls / reads the inherited
//     members.
//   - OCCT BOPAlgo_MakerVolume::myPaveFiller is a raw pointer created in
//     Perform(); rcad owns it (Option<PaveFiller>). Because the owned
//     PaveFiller lends its DS to the composed Builder, the pipeline runs
//     while the PaveFiller is a local of Perform() and is stored back into
//     my_pave_filler afterwards (a self-referential struct is impossible in
//     safe Rust).
// Architecture differences vs OCCT (no equivalent in Rust):
//   - Message_ProgressScope / BOPAlgo_PISteps weighting
//     (BOPAlgo_MakerVolume.cxx L48-50, L105, L124-125, L190-209) has no rcad
//     translation; rcad scopes are created with the OCCT step counts but the
//     per-stage weights (analyzeProgress / fillPISteps) are not translated.
//   - NCollection_BaseAllocator (L60-62) has no rcad equivalent; PaveFiller
//     owns its containers.

use crate::bop::algo::builder::Builder;
use crate::bop::algo::builder_solid::BuilderSolid;
use crate::bop::algo::pave_filler::PaveFiller;
use crate::bop::algo::{Alert, BooleanOpType, GlueEnum, Report};
use crate::bop::ds::DS;
use rcad_kernel::core::message::{NoopProgress, ProgressScope};
use rcad_kernel::math::bnd::BndBox;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{Orientation, ShapeType, TShape};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// OCCT BOPAlgo_MakerVolume — solid building from solids/faces.
pub struct MakerVolume {
    // --- BOPAlgo_Options / BOPAlgo_Algo (inherited) ---
    my_report: Report,           // BOPAlgo_Algo::myReport
    my_run_parallel: bool,       // BOPAlgo_Algo::myRunParallel
    my_fuzzy_value: f64,         // BOPAlgo_Algo::myFuzzyValue
    // --- BOPAlgo_BuilderShape (inherited) ---
    my_shape: Option<Shape>,     // BOPAlgo_BuilderShape::myShape
    my_fill_history: bool,       // BOPAlgo_BuilderShape::myFillHistory
    // --- BOPAlgo_Builder (inherited; myArguments only — the rest of the
    //     Builder state is carried by the composed Builder instance) ---
    my_arguments: Vec<Shape>,    // BOPAlgo_Builder.hxx L492: myArguments
    my_entry_point: i32,         // BOPAlgo_Builder.hxx L498: myEntryPoint
    // --- BOPAlgo_MakerVolume.hxx L200-204 ---
    my_intersect: bool,          // L200: myIntersect
    my_bbox: BndBox,             // L201: myBBox
    my_sbox: Option<Shape>,      // L202: mySBox (TopoDS_Solid)
    my_faces: Vec<Shape>,        // L203: myFaces
    my_avoid_internal_shapes: bool, // L204: myAvoidInternalShapes
    // --- BOPAlgo_Options (inherited options read by Perform L89-93) ---
    my_non_destructive: bool,    // BOPAlgo_Options::myNonDestructive
    my_glue: GlueEnum,           // BOPAlgo_Options::myGlue
    my_use_obb: bool,            // BOPAlgo_Options::myUseOBB
    // OCCT: BOPAlgo_MakerVolume::myPaveFiller (raw pointer, L56/L62/L96).
    // rcad owns it; it lends its DS to the composed Builder during perform().
    my_pave_filler: Option<PaveFiller>,
}

impl MakerVolume {
    /// OCCT BOPAlgo_MakerVolume::BOPAlgo_MakerVolume (empty constructor).
    /// OCCT BOPAlgo_MakerVolume.cxx L31-33:
    ///   myIntersect(true), myAvoidInternalShapes(false).
    pub fn new() -> Self {
        MakerVolume {
            my_report: Report::new(),
            my_run_parallel: false,
            // BOPAlgo_Options: myFuzzyValue(Precision::Confusion()).
            my_fuzzy_value: rcad_kernel::precision::CONFUSION,
            my_shape: None,
            my_fill_history: false,
            my_arguments: Vec::new(),
            my_entry_point: 0,
            my_intersect: true,
            my_bbox: BndBox::new(),
            my_sbox: None,
            my_faces: Vec::new(),
            my_avoid_internal_shapes: false,
            my_non_destructive: false,
            my_glue: GlueEnum::GlueOff,
            my_use_obb: false,
            my_pave_filler: None,
        }
    }

    /// OCCT BOPAlgo_Algo::SetRunParallel.
    pub fn set_run_parallel(&mut self, the_flag: bool) {
        self.my_run_parallel = the_flag;
    }

    /// OCCT BOPAlgo_Options::SetFuzzyValue
    /// (myFuzzyValue = max(theFuzz, Precision::Confusion())).
    pub fn set_fuzzy_value(&mut self, the_fuzzy_value: f64) {
        self.my_fuzzy_value = the_fuzzy_value.max(rcad_kernel::precision::CONFUSION);
    }

    /// OCCT BOPAlgo_Algo::SetArguments.
    pub fn set_arguments(&mut self, the_ls: Vec<Shape>) {
        self.my_arguments = the_ls;
    }

    /// OCCT BOPAlgo_MakerVolume::SetIntersect (hxx L126).
    pub fn set_intersect(&mut self, b_intersect: bool) {
        self.my_intersect = b_intersect;
    }

    /// OCCT BOPAlgo_MakerVolume::IsIntersect (hxx L129).
    pub fn is_intersect(&self) -> bool {
        self.my_intersect
    }

    /// OCCT BOPAlgo_MakerVolume::Box (hxx L132) — returns the solid box.
    pub fn box_shape(&self) -> Option<&Shape> {
        self.my_sbox.as_ref()
    }

    /// OCCT BOPAlgo_MakerVolume::Faces (hxx L135) — the processed faces.
    pub fn faces(&self) -> &[Shape] {
        &self.my_faces
    }

    /// OCCT BOPAlgo_MakerVolume::SetAvoidInternalShapes (hxx L139-142).
    pub fn set_avoid_internal_shapes(&mut self, the_avoid_internal: bool) {
        self.my_avoid_internal_shapes = the_avoid_internal;
    }

    /// OCCT BOPAlgo_MakerVolume::IsAvoidInternalShapes (hxx L145).
    pub fn is_avoid_internal_shapes(&self) -> bool {
        self.my_avoid_internal_shapes
    }

    /// OCCT BOPAlgo_Options::SetNonDestructive.
    pub fn set_non_destructive(&mut self, the_flag: bool) {
        self.my_non_destructive = the_flag;
    }

    /// OCCT BOPAlgo_Options::SetGlue.
    pub fn set_glue(&mut self, the_glue: GlueEnum) {
        self.my_glue = the_glue;
    }

    /// OCCT BOPAlgo_Options::SetUseOBB.
    pub fn set_use_obb(&mut self, the_flag: bool) {
        self.my_use_obb = the_flag;
    }

    /// OCCT BOPAlgo_Algo::HasErrors.
    pub fn has_errors(&self) -> bool {
        self.my_report.has_errors()
    }

    /// The result shape (OCCT BOPAlgo_BuilderShape::Shape):
    /// - empty compound - if no solids were created;
    /// - solid - if created only one solid;
    /// - compound of solids - if created more than one solid.
    pub fn shape(&self) -> Option<&Shape> {
        self.my_shape.as_ref()
    }

    /// OCCT BOPAlgo_MakerVolume::CheckData (BOPAlgo_MakerVolume.cxx L33-42).
    fn check_data(&mut self) {
        // OCCT L35-39: if (myArguments.IsEmpty()) -> TooFewArguments.
        if self.my_arguments.is_empty() {
            self.my_report.add_error(Alert::TooFewArguments); // no arguments to process
            return;
        }
        //
        // OCCT L41: CheckFiller();
        // rcad: the PaveFiller is always constructed by perform() before the
        // pipeline runs and its report is merged by the caller; there is no
        // nullable myPaveFiller to check (same as Builder::check_filler in
        // builder.rs).
    }

    /// OCCT BOPAlgo_MakerVolume::Perform (BOPAlgo_MakerVolume.cxx L46-98).
    pub fn perform(&mut self) {
        // OCCT L48-50: progress scope "Performing MakeVolume operation" (10),
        // anInterPart = myIntersect ? 9 : 0.5, aBuildPart = 10 - anInterPart.
        // rcad: progress weighting is not translated (see file header).
        let a_prog = NoopProgress;
        let _a_ps = ProgressScope::new(&a_prog, "Performing MakeVolume operation", 10);

        // OCCT L52: GetReport()->Clear();
        self.my_report.clear();
        //
        // OCCT L54-58: if (myEntryPoint == 1) { delete myPaveFiller; ... }
        if self.my_entry_point == 1 {
            self.my_pave_filler = None;
        }
        //
        // OCCT L60-62: allocator + new BOPAlgo_PaveFiller.
        // rcad: no allocators; the PaveFiller is a plain value.
        let mut p_pf = PaveFiller::new();
        //
        // OCCT L64-83: if (!myIntersect) wrap all arguments into one compound
        // and use it as the single argument of the PaveFiller.
        if !self.my_intersect {
            let mut an_args_children: Vec<Shape> = Vec::new();
            for a_s in &self.my_arguments {
                // OCCT L78: aBB.Add(anArgs, aS);
                an_args_children.push(a_s.clone());
            }
            let an_args =
                Shape::new(Arc::new(TShape::Compound(an_args_children)), 0, Orientation::Forward);
            let mut a_ls: Vec<Shape> = Vec::new();
            a_ls.push(an_args); // OCCT L80: aLS.Append(anArgs);
            //
            p_pf.set_arguments(a_ls); // OCCT L82: pPF->SetArguments(aLS);
        } else {
            // OCCT L86: pPF->SetArguments(myArguments);
            p_pf.set_arguments(self.my_arguments.clone());
        }
        //
        // OCCT L89: pPF->SetRunParallel(myRunParallel);
        // Interface gap: rcad PaveFiller has no public SetRunParallel.
        // OCCT L90: pPF->SetFuzzyValue(myFuzzyValue);
        p_pf.set_fuzzy_value(self.my_fuzzy_value);
        // OCCT L91: pPF->SetNonDestructive(myNonDestructive);
        // Interface gap: rcad PaveFiller's non-destructive flag is internal.
        // OCCT L92: pPF->SetGlue(myGlue);
        // Interface gap: rcad PaveFiller::set_glue has a different signature
        // (bool + tolerance) and is not called here to avoid inventing
        // semantics. OCCT L93: pPF->SetUseOBB(myUseOBB); — rcad has no OBB.
        //
        // OCCT L94: pPF->Perform(aPS.Next(anInterPart));
        let a_ps = ProgressScope::new(&a_prog, "intersect", 100);
        p_pf.perform(&a_ps);
        //
        // OCCT L96-97: myEntryPoint = 1; PerformInternal(*pPF, ...);
        self.my_entry_point = 1;
        {
            // The composed Builder borrows the PaveFiller's DS; both are
            // locals of this scope (see the struct-level architecture note).
            self.perform_internal(&p_pf);
        }
        self.my_pave_filler = Some(p_pf);
    }

    /// OCCT BOPAlgo_Algo::PerformInternal -> BOPAlgo_MakerVolume::PerformInternal1
    /// (BOPAlgo_MakerVolume.cxx L102-194).
    ///
    /// `the_filler` plays the role of OCCT's `theFiller` argument; the
    /// inherited Builder state (myDS, myContext, myImages) is carried by the
    /// composed `Builder` constructed here from the filler's DS.
    fn perform_internal(&mut self, the_filler: &PaveFiller) {
        // OCCT L105: Message_ProgressScope aPS(theRange, "Building volumes", 100);
        // rcad: progress weighting is not translated.
        // OCCT L106-108: myPaveFiller = &theFiller; myDS = ...PDS();
        // myContext = ...Context(); — carried by the composed Builder over
        // the filler's DS.
        let ds: &DS = the_filler.ds();
        let mut builder = Builder::new(ds, BooleanOpType::Unknown, self.my_fuzzy_value);
        //
        // OCCT L111: CheckData
        self.check_data();
        if self.has_errors() {
            return;
        }
        //
        // OCCT L118: Prepare();
        // OCCT BOPAlgo_Builder::Prepare (BOPAlgo_Builder.cxx L156-164):
        // makes an empty COMPOUND as myShape. Builder::prepare is the rcad
        // translation for the composed Builder's own BRep result pool; the
        // inherited myShape of THIS object is the OCCT-style Shape below
        // (written by build_shape, read by Shape()/post_treat), so the
        // three-line effect of Prepare on myShape is inlined here.
        self.my_shape = Some(Shape::new(Arc::new(TShape::Compound(Vec::new())), 0, Orientation::Forward));
        if self.has_errors() {
            return;
        }
        //
        // OCCT L124-125: BOPAlgo_PISteps aSteps(PIOperation_Last);
        // analyzeProgress(100., aSteps); — progress weighting is not
        // translated (see file header).
        //
        // OCCT L127-154: Fill Images (only when myIntersect)
        if self.my_intersect {
            // OCCT L131-153: FillImagesVertices, FillImagesEdges,
            // FillImagesContainers(TopAbs_WIRE), FillImagesFaces.
            self.fill_images(&mut builder);
            if self.has_errors() {
                return;
            }
        }
        //
        // OCCT L157: CollectFaces();
        self.collect_faces(&builder, ds);
        if self.has_errors() {
            return;
        }
        //
        // OCCT L163-164: aBoxFaces map + aLSR list.
        let mut a_box_faces: HashSet<(u64, u32)> = HashSet::new();
        let mut a_lsr: Vec<Shape> = Vec::new();
        //
        // OCCT L167: MakeBox(aBoxFaces);
        self.make_box(&mut a_box_faces);
        //
        // OCCT L170: BuildSolids(aLSR, ...);
        self.build_solids(&mut a_lsr, ds);
        if self.has_errors() {
            return;
        }
        //
        // OCCT L177: RemoveBox(aLSR, aBoxFaces);
        self.remove_box(&mut a_lsr, &a_box_faces);
        //
        // OCCT L180: FillInternalShapes(aLSR);
        self.fill_internal_shapes(&a_lsr, ds, &builder);
        //
        // OCCT L183: BuildShape(aLSR);
        self.build_shape(&a_lsr);
        //
        // OCCT L186: PrepareHistory(...);
        self.prepare_history();
        if self.has_errors() {
            return;
        }
        //
        // OCCT L193: PostTreat(...);
        self.post_treat();
    }

    /// OCCT BOPAlgo_MakerVolume::PerformInternal1 L127-154: the inherited
    /// BOPAlgo_Builder image stages.
    ///
    /// The four stages are direct calls into the composed Builder's
    /// pub(crate) translations. MakerVolume's override runs NO BuildResult
    /// steps — BuildResult belongs to BOPAlgo_Builder::PerformInternal1 only
    /// (BOPAlgo_Builder.cxx L459-550); MakerVolume consumes myImages in
    /// CollectFaces instead. OCCT's inherited stages write into this same
    /// object's myReport; rcad keeps the composed Builder's report, merged
    /// after every stage (each merge carries exactly the alerts of the stage
    /// just run, so the per-stage HasErrors checks below read the same state
    /// as OCCT's single shared report).
    fn fill_images(&mut self, builder: &mut Builder) {
        // OCCT L131: FillImagesVertices(...);
        builder.fill_images_vertices();
        self.my_report.merge(builder.report().clone());
        // OCCT L132-135: if (HasErrors()) return;
        if self.has_errors() {
            return;
        }
        // OCCT L137: FillImagesEdges(...);
        builder.fill_images_edges();
        self.my_report.merge(builder.report().clone());
        // OCCT L138-141: if (HasErrors()) return;
        if self.has_errors() {
            return;
        }
        // OCCT L143: FillImagesContainers(TopAbs_WIRE, ...);
        builder.fill_images_containers(ShapeType::Wire);
        self.my_report.merge(builder.report().clone());
        // OCCT L144-147: if (HasErrors()) return;
        if self.has_errors() {
            return;
        }
        // OCCT L149: FillImagesFaces(...);
        builder.fill_images_faces();
        self.my_report.merge(builder.report().clone());
        // OCCT L150-153: if (HasErrors()) return;
        if self.has_errors() {
            return;
        }
    }

    /// OCCT BOPAlgo_MakerVolume::CollectFaces (BOPAlgo_MakerVolume.cxx
    /// L213-251).
    fn collect_faces(&mut self, builder: &Builder, ds: &DS) {
        //
        // OCCT L218: aMFence — TopTools_ShapeMapHasher (TShape + Location).
        let mut a_m_fence: HashSet<(u64, u32)> = HashSet::new();
        //
        // OCCT L220: aNbShapes = myDS->NbSourceShapes();
        let a_nb_shapes = ds.nb_source_shapes();
        for i in 0..a_nb_shapes {
            // OCCT L223-227: skip non-FACE source shapes.
            let a_si = ds.shape_info(i);
            if a_si.shape_type != ShapeType::Face {
                continue;
            }
            //
            // OCCT L229-230: myBBox.Add(aSI.Box());
            self.my_bbox.add_box(&a_si.bbox);
            //
            // OCCT L232-249: add the images when the face is bound in
            // myImages, otherwise add the face itself.
            let a_f = a_si.shape();
            if let Some(a_lf_im) = builder.my_images.get((a_f.ptr_id(), a_f.location)) {
                for a_f_im in a_lf_im.clone() {
                    // OCCT L240: if (aMFence.Add(aFIm))
                    if a_m_fence.insert((a_f_im.ptr_id(), a_f_im.location)) {
                        add_face(&a_f_im, &mut self.my_faces);
                    }
                }
            } else {
                add_face(a_f, &mut self.my_faces);
            }
        }
    }

    /// OCCT BOPAlgo_MakerVolume::MakeBox (BOPAlgo_MakerVolume.cxx L255-276).
    fn make_box(&mut self, the_box_faces: &mut HashSet<(u64, u32)>) {
        //
        // OCCT L261-263: anExt = sqrt(myBBox.SquareExtent()) * 0.5;
        // myBBox.Enlarge(anExt); myBBox.Get(...);
        let an_ext = square_extent(&self.my_bbox).sqrt() * 0.5;
        self.my_bbox.enlarge(an_ext);
        // OCCT L263: myBBox.Get(aXmin, aYmin, aZmin, aXmax, aYmax, aZmax);
        // A void box has no corners: OCCT's MakeBox would raise
        // Standard_ConstructionError; rcad panics (same failure class).
        let (a_xmin, a_ymin, a_zmin, a_xmax, a_ymax, a_zmax) =
            self.my_bbox.get().expect("MakeBox: the covering box is void");
        //
        // OCCT L265: gp_Pnt aPMin(...), aPMax(...);
        let a_p_min = glam::DVec3::new(a_xmin, a_ymin, a_zmin);
        let a_p_max = glam::DVec3::new(a_xmax, a_ymax, a_zmax);
        //
        // OCCT L267: mySBox = BRepPrimAPI_MakeBox(aPMin, aPMax).Solid();
        // rcad: the primitive returns a flat result BRep; the solid TShape is
        // extracted as the result Shape (same extraction as the boolean
        // pipeline's root-shape extraction). OCCT raises
        // Standard_ConstructionError for a degenerate box; rcad panics.
        let a_box_brep = rcad_modeling::make_box_brep(
            a_p_min,
            glam::DVec3::X,
            glam::DVec3::Y,
            a_p_max.x - a_p_min.x,
            a_p_max.y - a_p_min.y,
            a_p_max.z - a_p_min.z,
        )
        .expect("BRepPrimAPI_MakeBox failed");
        let a_solid_pos = a_box_brep
            .tshapes
            .iter()
            .rposition(|ts| matches!(ts.as_ref(), TShape::Solid(_)))
            .expect("MakeBox: no solid in the box BRep");
        self.my_sbox =
            Some(Shape::new(a_box_brep.tshapes[a_solid_pos].clone(), 0, Orientation::Forward));
        //
        // OCCT L269-275: explore the faces of mySBox, append them to myFaces
        // and add them to theBoxFaces.
        let a_sbox = self.my_sbox.as_ref().unwrap().clone();
        for a_f in explore_faces(&a_sbox) {
            // OCCT L273: myFaces.Append(aF);
            self.my_faces.push(a_f.clone());
            // OCCT L274: theBoxFaces.Add(aF);
            the_box_faces.insert((a_f.ptr_id(), a_f.location));
        }
    }

    /// OCCT BOPAlgo_MakerVolume::BuildSolids (BOPAlgo_MakerVolume.cxx
    /// L280-298).
    fn build_solids(&mut self, the_lsr: &mut Vec<Shape>, ds: &DS) {
        // OCCT L283: BOPAlgo_BuilderSolid aBS;
        let mut a_bs = BuilderSolid::new(ds);
        //
        // OCCT L285-287: aBS.SetShapes(myFaces); aBS.SetRunParallel(...);
        // aBS.SetAvoidInternalShapes(...);
        a_bs.my_shapes = self.my_faces.clone();
        // Interface gap: rcad BuilderSolid has no SetRunParallel.
        a_bs.set_avoid_internal_shapes(self.my_avoid_internal_shapes);
        // OCCT L288: aBS.Perform(theRange);
        a_bs.perform();
        // OCCT L289-293: if (aBS.HasErrors()) -> AlertSolidBuilderFailed.
        if a_bs.has_errors() {
            self.my_report.add_error(Alert::SolidBuilderFailed); // SolidBuilder failed
            return;
        }
        //
        // OCCT L295: myReport->Merge(aBS.GetReport());
        let a_bs_report = a_bs.report().clone();
        self.my_report.merge(a_bs_report);
        //
        // OCCT L297: theLSR = aBS.Areas();
        *the_lsr = a_bs.my_solids.clone();
    }

    /// OCCT BOPAlgo_MakerVolume::RemoveBox (BOPAlgo_MakerVolume.cxx
    /// L302-333).
    fn remove_box(&mut self, the_lsr: &mut Vec<Shape>, the_box_faces: &HashSet<(u64, u32)>) {
        //
        // OCCT L311: bFound = false;
        let mut b_found = false;
        // OCCT L312-332: iterate theLSR; remove the FIRST solid that contains
        // a face of theBoxFaces, then stop.
        let mut i = 0;
        while i < the_lsr.len() {
            let a_sr = the_lsr[i].clone();
            //
            // OCCT L317-327: TopExp_Explorer(aSR, TopAbs_FACE).
            for a_f in explore_faces(&a_sr) {
                if the_box_faces.contains(&(a_f.ptr_id(), a_f.location)) {
                    b_found = true;
                    // OCCT L324: theLSR.Remove(aIt);
                    the_lsr.remove(i);
                    break;
                }
            }
            if b_found {
                // OCCT L328-331: if (bFound) break;
                break;
            }
            // OCCT L313: aIt.Next();
            i += 1;
        }
    }

    /// OCCT BOPAlgo_MakerVolume::FillInternalShapes (BOPAlgo_MakerVolume.cxx
    /// L359-403).
    fn fill_internal_shapes(&mut self, the_lsr: &[Shape], ds: &DS, builder: &Builder) {
        // OCCT L361-364: if (myAvoidInternalShapes) return;
        if self.my_avoid_internal_shapes {
            return;
        }
        //
        // OCCT L367: aLSC — all non-compound shapes of the arguments.
        let mut a_lsc: Vec<Shape> = Vec::new();
        // OCCT L369: aMFence — shared by TreatCompound and the WIRE walk.
        let mut a_m_fence: HashSet<(u64, u32)> = HashSet::new();
        //
        // OCCT L371-375: TreatCompound of every argument into aLSC.
        for a_s in ds.arguments() {
            // OCCT L374: BOPTools_AlgoTools::TreatCompound(aS, aLSC, &aMFence);
            // Interface gap: rcad treat_compound has no fence parameter; the
            // fence check is applied at the call site with the same map.
            for a_non_c in crate::bop::tools::algo_tools::treat_compound(a_s) {
                if a_m_fence.insert((a_non_c.ptr_id(), a_non_c.location)) {
                    a_lsc.push(a_non_c);
                }
            }
        }
        //
        // OCCT L378: aLVE — only the edges and vertices from the arguments.
        let mut a_lve: Vec<Shape> = Vec::new();
        //
        // OCCT L380-400.
        for a_s in &a_lsc {
            let a_type = a_s.shape_type();
            if a_type == ShapeType::Wire {
                // OCCT L387-395: TopoDS_Iterator(aS) — the direct sub-shapes
                // of the wire (its edges), deduplicated by the same fence.
                if let TShape::Wire(wd) = &*a_s.data {
                    for a_ss in &wd.edges {
                        // OCCT L390: if (aMFence.Add(aSS))
                        if a_m_fence.insert((a_ss.ptr_id(), a_ss.location)) {
                            a_lve.push(a_ss.clone());
                        }
                    }
                }
            } else if a_type == ShapeType::Vertex || a_type == ShapeType::Edge {
                // OCCT L396-399: aLVE.Append(aS);
                a_lve.push(a_s.clone());
            }
        }
        //
        // OCCT L402: BOPAlgo_Tools::FillInternals(theLSR, aLVE, myImages,
        // myContext);
        // OCCT BOPAlgo_Tools.cxx L1755-1758: FillInternals returns immediately
        // when theSolids or theParts is empty; that early return is mirrored
        // here. INTERFACE GAP: the rest of BOPAlgo_Tools::FillInternals
        // (classification of the internal V/E/F parts and their addition into
        // the solids, BOPAlgo_Tools.cxx L1760-1900) is not translated in rcad;
        // the non-empty path is a no-op pending that translation.
        if the_lsr.is_empty() || a_lve.is_empty() {
            return;
        }
        let _ = builder;
    }

    /// OCCT BOPAlgo_MakerVolume::BuildShape (BOPAlgo_MakerVolume.cxx
    /// L337-355).
    fn build_shape(&mut self, the_lsr: &[Shape]) {
        if the_lsr.len() == 1 {
            // OCCT L339-341: myShape = theLSR.First();
            self.my_shape = Some(the_lsr[0].clone());
        } else {
            // OCCT L344-354: BRep_Builder aBB; for every solid in theLSR:
            // aBB.Add(myShape, aSol); — myShape is the empty compound made by
            // Prepare(), so adding the solids turns it into a compound.
            let mut a_children: Vec<Shape> = Vec::new();
            for a_sol in the_lsr {
                a_children.push(a_sol.clone());
            }
            self.my_shape =
                Some(Shape::new(Arc::new(TShape::Compound(a_children)), 0, Orientation::Forward));
        }
    }

    /// OCCT BOPAlgo_BuilderShape::PrepareHistory, called at
    /// BOPAlgo_MakerVolume.cxx L186.
    fn prepare_history(&mut self) {
        // OCCT BOPAlgo_BuilderShape.cxx L166-168:
        // if (!HasHistory()) return; — the default myFillHistory is FALSE.
        // Interface gap: the history part of PrepareHistory (myOrigins /
        // myImages -> BRepTools_History, needed by BRepOffset's UpdateHistory)
        // is not translated in rcad; Builder::prepare_history is private and
        // works on the topods::BRep result pool of builder.rs.
        if !self.my_fill_history {
            return;
        }
    }

    /// OCCT BOPAlgo_Algo::PostTreat, called at BOPAlgo_MakerVolume.cxx L193
    /// (BOPAlgo_Algo.cxx L466-486).
    fn post_treat(&mut self) {
        // OCCT L466-480: aMapToAvoid is empty (non-destructive mode is off).
        let a_ma: HashSet<(u64, u32)> = HashSet::new();
        if let Some(a_shape) = &self.my_shape {
            // OCCT L483: CorrectTolerances(myShape, aMA, 0.05, myRunParallel);
            // OCCT L485: CorrectShapeTolerances(myShape, aMA, myRunParallel);
            // rcad: the corrections take a result BRep pool; the result Shape
            // is flattened into a view BRep, corrected in place (the TShape
            // Arcs are shared, so the corrections persist) and dropped.
            let mut a_brep = flatten_to_brep(a_shape);
            crate::bop::tools::algo_tools::correct_tolerances(&mut a_brep, &a_ma, 0.05);
            crate::bop::tools::algo_tools::correct_shape_tolerances(&mut a_brep, &a_ma);
        }
    }
}

/// OCCT AddFace (BOPAlgo_MakerVolume.cxx L407-414) — appends the face twice,
/// with FORWARD and REVERSED orientation.
fn add_face(the_f: &Shape, the_lf: &mut Vec<Shape>) {
    let mut a_ff = the_f.clone();
    a_ff.orientation = Orientation::Forward; // OCCT L410
    the_lf.push(a_ff); // OCCT L411
    let mut a_ff = the_f.clone();
    a_ff.orientation = Orientation::Reversed; // OCCT L412
    the_lf.push(a_ff); // OCCT L413
}

/// OCCT Bnd_Box::SquareExtent (Bnd_Box.hxx) — the squared diagonal length of
/// the gap-included box; 0. for a void box. rcad BndBox has no
/// square_extent, so it is computed from get() (which returns the
/// gap-included corners like OCCT's Get()).
fn square_extent(the_box: &BndBox) -> f64 {
    if the_box.is_void() {
        return 0.0;
    }
    match the_box.get() {
        Some((a_xmin, a_ymin, a_zmin, a_xmax, a_ymax, a_zmax)) => {
            let a_dx = a_xmax - a_xmin;
            let a_dy = a_ymax - a_ymin;
            let a_dz = a_zmax - a_zmin;
            a_dx * a_dx + a_dy * a_dy + a_dz * a_dz
        }
        None => 0.0,
    }
}

/// OCCT TopExp_Explorer(theShape, TopAbs_FACE) — all faces of a
/// solid/compound/shell, each with the cumulative orientation composed
/// (TopoDS_Iterator.cxx L72-80). The stored order is preserved (FIFO), the
/// same traversal model as collect_solid_faces in builder.rs.
fn explore_faces(the_shape: &Shape) -> Vec<Shape> {
    let mut a_result: Vec<Shape> = Vec::new();
    let mut a_queue: std::collections::VecDeque<(Shape, Orientation)> =
        std::collections::VecDeque::new();
    a_queue.push_back((the_shape.clone(), Orientation::Forward));
    while let Some((a_sh, a_cum_or)) = a_queue.pop_front() {
        let a_or = a_cum_or.compose(a_sh.orientation);
        match &*a_sh.data {
            TShape::Solid(a_sd) => {
                for a_x in &a_sd.shells {
                    a_queue.push_back((a_x.clone(), a_or));
                }
            }
            TShape::CompSolid(a_cd) => {
                for a_x in a_cd {
                    a_queue.push_back((a_x.clone(), a_or));
                }
            }
            TShape::Compound(a_cd) => {
                for a_x in a_cd {
                    a_queue.push_back((a_x.clone(), a_or));
                }
            }
            TShape::Shell(a_sd) => {
                for a_x in &a_sd.faces {
                    a_queue.push_back((a_x.clone(), a_or));
                }
            }
            TShape::Face(_) => {
                let mut a_f = a_sh.clone();
                a_f.orientation = a_or;
                a_result.push(a_f);
            }
            _ => {}
        }
    }
    a_result
}

/// Flatten a Shape tree into a view BRep pool (each reachable TShape once,
/// deduplicated by TShape identity, TShape-internal reference indices
/// re-pointed to the view pool in place).
///
/// Mirror of builder.rs Builder::push_shape_recursive (private there);
/// representation bridge for the tolerance corrections of post_treat and for
/// the volume checks in the tests. The view shares the TShape Arcs, so
/// in-place edits persist (single-threaded pipeline; OCCT BRep_Builder::Add
/// references the source TShape, it never clones it).
fn flatten_to_brep(the_shape: &Shape) -> rcad_kernel::topods::BRep {
    fn push_recursive(
        a_shape: &Shape,
        a_brep: &mut rcad_kernel::topods::BRep,
        a_remap: &mut HashMap<u64, usize>,
    ) -> usize {
        let a_ptr = a_shape.ptr_id();
        if let Some(&a_idx) = a_remap.get(&a_ptr) {
            return a_idx;
        }
        let a_new_idx = a_brep.tshapes.len();
        a_brep.tshapes.push(a_shape.data.clone());
        a_remap.insert(a_ptr, a_new_idx);
        // Recursively push the sub-shapes first so their view indices exist.
        match a_shape.data.as_ref() {
            TShape::Edge(a_ed) => {
                let _ = push_recursive(&a_ed.first, a_brep, a_remap);
                let _ = push_recursive(&a_ed.last, a_brep, a_remap);
            }
            TShape::Wire(a_wd) => {
                for a_e in &a_wd.edges {
                    let _ = push_recursive(a_e, a_brep, a_remap);
                }
            }
            TShape::Face(a_fd) => {
                let _ = push_recursive(&a_fd.outer_wire, a_brep, a_remap);
                for a_w in &a_fd.inner_wires {
                    let _ = push_recursive(a_w, a_brep, a_remap);
                }
                for a_v in &a_fd.internal_vertices {
                    let _ = push_recursive(a_v, a_brep, a_remap);
                }
            }
            TShape::Shell(a_sd) => {
                for a_f in &a_sd.faces {
                    let _ = push_recursive(a_f, a_brep, a_remap);
                }
            }
            TShape::Solid(a_sd) => {
                for a_s in &a_sd.shells {
                    let _ = push_recursive(a_s, a_brep, a_remap);
                }
                for a_v in &a_sd.internal_vertices {
                    let _ = push_recursive(a_v, a_brep, a_remap);
                }
                for a_e in &a_sd.internal_edges {
                    let _ = push_recursive(a_e, a_brep, a_remap);
                }
            }
            TShape::CompSolid(a_shapes) => {
                for a_s in a_shapes {
                    let _ = push_recursive(a_s, a_brep, a_remap);
                }
            }
            TShape::Compound(a_shapes) => {
                for a_s in a_shapes {
                    let _ = push_recursive(a_s, a_brep, a_remap);
                }
            }
            TShape::Vertex(_) => {}
        }
        // Re-point the TShape-internal reference indices from their original
        // pool positions to the view pool positions, in place on the shared
        // TShape.
        let a_raw = Arc::as_ptr(&a_shape.data) as *mut TShape;
        let a_lookup = |a_remap: &HashMap<u64, usize>, a_s: &Shape| {
            a_remap.get(&a_s.ptr_id()).copied().unwrap_or(a_s.index)
        };
        unsafe {
            match &mut *a_raw {
                TShape::Vertex(a_vd) => {
                    for a_s in a_vd.my_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                }
                TShape::Edge(a_ed) => {
                    for a_s in a_ed.my_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                    a_ed.first.index = a_lookup(a_remap, &a_ed.first);
                    a_ed.last.index = a_lookup(a_remap, &a_ed.last);
                }
                TShape::Wire(a_wd) => {
                    for a_s in a_wd.my_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                    for a_e in a_wd.edges.iter_mut() {
                        a_e.index = a_lookup(a_remap, a_e);
                    }
                }
                TShape::Face(a_fd) => {
                    for a_s in a_fd.my_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                    a_fd.outer_wire.index = a_lookup(a_remap, &a_fd.outer_wire);
                    for a_w in a_fd.inner_wires.iter_mut() {
                        a_w.index = a_lookup(a_remap, a_w);
                    }
                    for a_v in a_fd.internal_vertices.iter_mut() {
                        a_v.index = a_lookup(a_remap, a_v);
                    }
                }
                TShape::Shell(a_sd) => {
                    for a_s in a_sd.my_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                    for a_f in a_sd.faces.iter_mut() {
                        a_f.index = a_lookup(a_remap, a_f);
                    }
                }
                TShape::Solid(a_sd) => {
                    for a_s in a_sd.my_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                    for a_sh in a_sd.shells.iter_mut() {
                        a_sh.index = a_lookup(a_remap, a_sh);
                    }
                    for a_v in a_sd.internal_vertices.iter_mut() {
                        a_v.index = a_lookup(a_remap, a_v);
                    }
                    for a_e in a_sd.internal_edges.iter_mut() {
                        a_e.index = a_lookup(a_remap, a_e);
                    }
                }
                TShape::CompSolid(a_shapes) => {
                    for a_s in a_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                }
                TShape::Compound(a_shapes) => {
                    for a_s in a_shapes.iter_mut() {
                        a_s.index = a_lookup(a_remap, a_s);
                    }
                }
            }
        }
        a_new_idx
    }
    let mut a_brep = rcad_kernel::topods::BRep::new();
    let mut a_remap: HashMap<u64, usize> = HashMap::new();
    let _ = push_recursive(the_shape, &mut a_brep, &mut a_remap);
    a_brep
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::Plane;

    /// Wrap a result BRep's tshapes into one compound argument Shape (the
    /// same argument wrapping as bop::ds::topods_builder::new_from_topods).
    fn brep_arg_shape(brep: &rcad_kernel::topods::BRep) -> Shape {
        Shape::new(
            Arc::new(TShape::Compound(
                brep.tshapes
                    .iter()
                    .map(|ts| Shape::new(ts.clone(), 0, Orientation::Forward))
                    .collect(),
            )),
            0,
            Orientation::Forward,
        )
    }

    fn box_arg(origin: DVec3, dx: f64, dy: f64, dz: f64) -> Shape {
        let brep = rcad_modeling::make_box_brep(origin, DVec3::X, DVec3::Y, dx, dy, dz)
            .expect("box creation failed");
        brep_arg_shape(&brep)
    }

    fn face_arg(plane: &Plane, u1: f64, u2: f64, v1: f64, v2: f64) -> Shape {
        let brep = rcad_modeling::make_face_plane_bounds_brep(plane, u1, u2, v1, v2)
            .expect("face creation failed");
        brep_arg_shape(&brep)
    }

    /// Collect the solids of a result shape (compound of solids or a solid).
    fn collect_solids(s: &Shape) -> Vec<Shape> {
        let mut out = Vec::new();
        let mut stack = vec![s.clone()];
        while let Some(sh) = stack.pop() {
            match &*sh.data {
                TShape::Compound(children) => stack.extend(children.iter().cloned()),
                TShape::Solid(_) => out.push(sh),
                _ => {}
            }
        }
        out
    }

    /// Flatten a shape into a view BRep and compute its volume
    /// (BRepGProp_Vinert equivalent, rcad_kernel::volume).
    fn shape_volume(s: &Shape) -> f64 {
        rcad_kernel::volume(&flatten_to_brep(s))
    }

    /// OCCT-faithful volume of a solid: rebuild the solid keeping only
    /// FORWARD/REVERSED faces before the integration.
    /// OCCT BRepGProp::volumePropertiesFaces skips INTERNAL/EXTERNAL faces
    /// (BRepGProp.cxx L352-355 + L369 `if (isFwd || isRvs)`); the rcad
    /// kernel volume integrates every face occurrence, so the INTERNAL
    /// faces of the internal shells must be dropped from the view first.
    fn solid_volume_skip_internal(solid: &Shape) -> f64 {
        use rcad_kernel::topods::{TSolidData, TShellData, tshape_flags};
        let TShape::Solid(a_sd) = &*solid.data else { return 0.0 };
        let mut a_shells: Vec<Shape> = Vec::new();
        for a_sh in &a_sd.shells {
            let TShape::Shell(a_shd) = &*a_sh.data else { continue };
            let a_faces: Vec<Shape> = a_shd
                .faces
                .iter()
                .filter(|a_f| {
                    a_f.orientation == Orientation::Forward
                        || a_f.orientation == Orientation::Reversed
                })
                .cloned()
                .collect();
            if a_faces.is_empty() {
                continue;
            }
            a_shells.push(Shape::new(
                Arc::new(TShape::Shell(TShellData {
                    my_shapes: vec![],
                    flags: a_shd.flags | tshape_flags::CLOSED,
                    faces: a_faces,
                })),
                0,
                Orientation::Forward,
            ));
        }
        let a_view = Shape::new(
            Arc::new(TShape::Solid(TSolidData {
                my_shapes: vec![],
                flags: a_sd.flags,
                shells: a_shells,
                internal_vertices: vec![],
                internal_edges: vec![],
            })),
            0,
            Orientation::Forward,
        );
        rcad_kernel::volume(&flatten_to_brep(&a_view))
    }

    fn count_faces(s: &Shape) -> usize {
        flatten_to_brep(s)
            .tshapes
            .iter()
            .filter(|ts| matches!(ts.as_ref(), TShape::Face(_)))
            .count()
    }

    /// Two intersecting boxes, default options (myIntersect = TRUE).
    ///
    /// OCCT ground truth (DRAWEXE, OCCT 8.0):
    ///   box b1 0 0 0 2 2 2; box b2 1 1 1 2 2 2; mkvolume r b1 b2
    ///   -> compound of 3 SOLIDs, volumes {7, 7, 1}: the two boxes split by
    ///   their overlap cell (1x1x1). 18 faces, 3 shells. MakerVolume produces
    ///   the non-manifold cell decomposition; it does NOT glue the cells into
    ///   a single manifold union solid.
    #[test]
    fn make_volume_two_intersecting_boxes() {
        let mut a_mv = MakerVolume::new();
        a_mv.set_arguments(vec![
            box_arg(DVec3::new(0.0, 0.0, 0.0), 2.0, 2.0, 2.0),
            box_arg(DVec3::new(1.0, 1.0, 1.0), 2.0, 2.0, 2.0),
        ]);
        assert!(a_mv.is_intersect());
        a_mv.perform();
        assert!(!a_mv.has_errors(), "MakeVolume failed");
        let a_r = a_mv.shape().expect("no result shape");
        let a_solids = collect_solids(a_r);
        assert_eq!(a_solids.len(), 3, "expected 3 solids (7, 7, 1)");
        let mut a_vols: Vec<f64> = a_solids.iter().map(|s| shape_volume(s)).collect();
        a_vols.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((a_vols[0] - 1.0).abs() < 1.0e-6, "overlap cell volume {}, expected 1", a_vols[0]);
        assert!((a_vols[1] - 7.0).abs() < 1.0e-6, "L-cell volume {}, expected 7", a_vols[1]);
        assert!((a_vols[2] - 7.0).abs() < 1.0e-6, "L-cell volume {}, expected 7", a_vols[2]);
        assert_eq!(count_faces(a_r), 18, "expected 18 faces");
        // Box() must not be part of the result: every remaining solid is
        // free of the covering box faces (OCCT RemoveBox).
        let a_box_faces: HashSet<(u64, u32)> = a_mv
            .box_shape()
            .map(|sbox| explore_faces(sbox).iter().map(|f| (f.ptr_id(), f.location)).collect())
            .unwrap_or_default();
        assert!(!a_box_faces.is_empty(), "the covering box must exist");
        for a_s in &a_solids {
            for a_f in explore_faces(a_s) {
                assert!(
                    !a_box_faces.contains(&(a_f.ptr_id(), a_f.location)),
                    "a covering-box face leaked into the result"
                );
            }
        }
    }

    /// A box with a face patch strictly inside it (no interference between
    /// the arguments).
    ///
    /// OCCT ground truth (DRAWEXE, OCCT 8.0):
    ///   box b 0 0 0 1 1 1 + mkplane f (0.25..0.75 square at z=0.5);
    ///   mkvolume r b f
    ///   -> 1 SOLID (not a compound), volume 1, 7 faces, 2 shells: the inner
    ///   face becomes an INTERNAL shell of the box solid.
    #[test]
    fn make_volume_box_with_internal_face() {
        let a_plane = Plane::new(DVec3::new(0.0, 0.0, 0.5), DVec3::Z);
        let mut a_mv = MakerVolume::new();
        a_mv.set_arguments(vec![
            box_arg(DVec3::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0),
            face_arg(&a_plane, 0.25, 0.75, 0.25, 0.75),
        ]);
        a_mv.perform();
        assert!(!a_mv.has_errors(), "MakeVolume failed");
        let a_r = a_mv.shape().expect("no result shape");
        // Result is the alone solid itself (OCCT BuildShape extent == 1).
        assert!(matches!(&*a_r.data, TShape::Solid(_)), "result must be a single solid");
        assert_eq!(count_faces(a_r), 7, "expected 6 box faces + 1 internal face");
        if let TShape::Solid(a_sd) = &*a_r.data {
            assert_eq!(a_sd.shells.len(), 2, "expected outer shell + internal shell");
            // The second shell holds the internal face with INTERNAL
            // orientation (BOPAlgo_BuilderSolid::MakeInternalShells).
            if let TShape::Shell(a_inner) = &*a_sd.shells[1].data {
                assert_eq!(a_inner.faces.len(), 1, "the internal shell holds the face patch");
                assert_eq!(a_inner.faces[0].orientation, Orientation::Internal);
            } else {
                panic!("second shell expected");
            }
        }
        // Volume 1 (OCCT checkprops -v ground truth): the INTERNAL face is
        // skipped by OCCT's volume integration.
        let a_vol = solid_volume_skip_internal(a_r);
        assert!((a_vol - 1.0).abs() < 1.0e-6, "volume {}, expected 1", a_vol);
        // With SetAvoidInternalShapes(true) the internal face is dropped.
        let a_plane2 = Plane::new(DVec3::new(0.0, 0.0, 0.5), DVec3::Z);
        let mut a_mv2 = MakerVolume::new();
        a_mv2.set_avoid_internal_shapes(true);
        a_mv2.set_arguments(vec![
            box_arg(DVec3::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0),
            face_arg(&a_plane2, 0.25, 0.75, 0.25, 0.75),
        ]);
        a_mv2.perform();
        assert!(!a_mv2.has_errors(), "MakeVolume (avoid internal) failed");
        let a_r2 = a_mv2.shape().expect("no result shape");
        assert!(matches!(&*a_r2.data, TShape::Solid(_)));
        assert_eq!(count_faces(a_r2), 6, "internal face must be avoided");
        let a_vol2 = solid_volume_skip_internal(a_r2);
        assert!((a_vol2 - 1.0).abs() < 1.0e-6);
    }

    /// myIntersect = FALSE: the arguments are wrapped into one compound and
    /// are not intersected (OCCT Perform L64-83); used by BRepOffset's
    /// aMV2/aMV3 calls.
    #[test]
    fn make_volume_no_intersect_single_cell() {
        let mut a_mv = MakerVolume::new();
        a_mv.set_intersect(false);
        a_mv.set_arguments(vec![
            box_arg(DVec3::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0),
            box_arg(DVec3::new(2.0, 0.0, 0.0), 1.0, 1.0, 1.0),
        ]);
        a_mv.perform();
        assert!(!a_mv.has_errors(), "MakeVolume failed");
        let a_r = a_mv.shape().expect("no result shape");
        let a_solids = collect_solids(a_r);
        assert_eq!(a_solids.len(), 2, "two disjoint boxes stay two solids");
        for a_s in &a_solids {
            assert!((shape_volume(a_s) - 1.0).abs() < 1.0e-6);
        }
    }
}
