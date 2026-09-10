// OCCT BOPAlgo_Section.cxx + BOPAlgo_Section.hxx — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKBO/BOPAlgo/BOPAlgo_Section.cxx
//         $OCCT_SRC/src/ModelingAlgorithms/TKBO/BOPAlgo/BOPAlgo_Section.hxx
//
// OCCT inheritance chain (BOPAlgo_Section.hxx L33):
//   BOPAlgo_Section : BOPAlgo_Builder
//     BOPAlgo_Builder : BOPAlgo_BuilderShape : BOPAlgo_Algo : BOPAlgo_Options
// Rust has no inheritance -> composition + delegation: the BOPAlgo_Builder
// base sub-object is rcad's `Builder` (crate::bop::algo::builder::Builder);
// the inherited data members (myDS, myContext, myArguments, myImages,
// myShape, myReport, ...) live in it and are reached through `self.builder`.
//
// Architecture differences (referenced from the affected functions):
// 1. OCCT binds myPaveFiller/myDS/myContext inside PerformInternal1
//    (Section.cxx L104-106); rcad borrows the DS at construction
//    (Builder::new) because the PaveFiller must outlive the Builder borrow.
// 2. TopExp_Explorer / TopoDS_Iterator are modelled by explorer() and
//    sub_shapes() below (TopExp_Explorer.cxx L68-171). Child order is the
//    rcad container field order used by builder.rs push_shape_recursive;
//    orientation is composed like TopAbs::Compose(parent, child) of
//    TopoDS_Iterator (cumOri=true). TopoDS_Iterator cumLoc=true is implicit:
//    rcad Shape nodes already carry absolute location indices.
// 3. The result pool construction mirrors builder.rs
//    Builder::push_shape_recursive (private there): BRep_Builder::Add over a
//    flat TShape pool, dedup by TShape identity, TShape-internal reference
//    indices re-pointed in place (see ResultPool below).
// 4. Progress indicator plumbing has no rcad Builder counterpart:
//    fillPIConstants / fillPISteps (Section.cxx L74-96), analyzeProgress
//    (L122-123) and the Message_ProgressScope / UserBreak checks (L196-199,
//    L289-292, L383-386) are therefore reduced to comments at the translated
//    spots.

use crate::bop::algo::builder::Builder;
use crate::bop::algo::{Alert, BooleanOpType};
use crate::bop::ds::DS;
use glam::DAffine3;
use indexmap::IndexMap;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{self, Orientation, ShapeType, TShape};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// OCCT BOPAlgo_Section — the algorithm to build a Section between the
/// arguments (BOPAlgo_Section.hxx L26-33).
///
/// The Section consists of vertices and edges:
/// 1. new vertices that are subjects of V/V, E/E, E/F, F/F interferences
/// 2. vertices that are subjects of V/E, V/F interferences
/// 3. new edges that are subjects of F/F interferences
/// 4. edges that are Common Blocks
pub struct BOPAlgoSection<'a> {
    /// OCCT BOPAlgo_Section : BOPAlgo_Builder base sub-object.
    pub(crate) builder: Builder<'a>,
}

impl<'a> BOPAlgoSection<'a> {
    /// OCCT BOPAlgo_Section::BOPAlgo_Section() (Section.cxx L38-42) and
    /// BOPAlgo_Section(theAllocator) (L46-50): both call Clear().
    ///
    /// Architecture difference #1: OCCT binds myDS from the PaveFiller inside
    /// PerformInternal1; rcad borrows the DS here (the PaveFiller must
    /// outlive the borrow), so the constructor takes it.
    pub fn new(ds: &'a DS, fuzzy_value: f64) -> Self {
        BOPAlgoSection {
            builder: Builder::new(ds, BooleanOpType::Section, fuzzy_value),
        }
    }

    /// OCCT BOPAlgo_Builder::SetArguments — myArguments.
    ///
    /// BOPAlgo_Section inherits BOPAlgo_Builder and has NO myTools (OCCT
    /// BRepAlgoAPI_BooleanOperation::Build L199-200 passes
    /// myDSFiller->Arguments() — objects and tools combined — as ONE list).
    pub fn set_arguments(&mut self, args: Vec<Shape>) {
        self.builder.set_arguments(args);
    }

    /// OCCT BOPAlgo_Algo::HasErrors.
    pub fn has_errors(&self) -> bool {
        self.builder.my_report.has_errors()
    }

    /// OCCT BOPAlgo_Algo::Perform — entry point. The PaveFiller has already
    /// run (rcad binds the DS at construction, architecture difference #1),
    /// so Perform goes straight to PerformInternal1.
    pub fn perform(&mut self) {
        self.perform_internal1();
    }

    /// OCCT BOPAlgo_Section::CheckData (Section.cxx L58-70).
    fn check_data(&mut self) {
        // OCCT L60-62: aNbArgs = myArguments.Extent();
        let a_nb_args = self.builder.my_arguments.len();
        // OCCT L63-67: if (!aNbArgs) { AddError(AlertTooFewArguments); return; }
        if a_nb_args == 0 {
            self.builder.my_report.add_error(Alert::TooFewArguments);
            return;
        }
        // OCCT L69: CheckFiller();
        self.check_filler();
    }

    /// OCCT BOPAlgo_Builder::CheckFiller (BOPAlgo_Builder.cxx L143-152),
    /// called from Section::CheckData.
    fn check_filler(&mut self) {
        // OCCT L145-149: if (!myPaveFiller) { AddError(AlertNoFiller); return; }
        // rcad: the PaveFiller always runs before the Builder (architecture
        // difference #1); there is no optional filler reference.
        // OCCT L151: GetReport()->Merge(myPaveFiller->GetReport());
        // rcad: the PaveFiller report is not carried over (the PaveFiller is
        // dropped before the Builder runs).
    }

    /// OCCT BOPAlgo_Section::PerformInternal1 (Section.cxx L100-163) — the
    /// literal step sequence.
    fn perform_internal1(&mut self) {
        // OCCT L104-106: myPaveFiller = &theFiller; myDS = PDS(); myContext =
        // Context() — bound at construction (architecture difference #1).
        //
        // 1. CheckData (OCCT L108-113)
        self.check_data();
        if self.has_errors() {
            return;
        }
        //
        // 2. Prepare (OCCT L115-120)
        self.builder.prepare();
        if self.has_errors() {
            return;
        }
        //
        // OCCT L122-123: BOPAlgo_PISteps aSteps(PIOperation_Last);
        // analyzeProgress(100., aSteps); — progress plumbing has no rcad
        // Builder counterpart (architecture difference #4); the aPS.Next(...)
        // progress-range arguments of the calls below are dropped for the
        // same reason.
        //
        // 3. Fill Images
        // 3.1 Vertices (OCCT L126-130)
        self.builder.fill_images_vertices();
        if self.has_errors() {
            return;
        }
        //
        // OCCT L132: BuildResult(TopAbs_VERTEX);
        self.builder.build_result(topods::ShapeType::Vertex);
        if self.has_errors() {
            return;
        }
        // 3.2 Edges (OCCT L138-142)
        self.builder.fill_images_edges();
        if self.has_errors() {
            return;
        }
        //
        // OCCT L144: BuildResult(TopAbs_EDGE);
        self.builder.build_result(topods::ShapeType::Edge);
        if self.has_errors() {
            return;
        }
        // 4. Section (OCCT L149-154)
        self.build_section();
        if self.has_errors() {
            return;
        }
        // 5. History (OCCT L155-160)
        self.builder.prepare_history();
        if self.has_errors() {
            return;
        }
        // 6. Post-treatment (OCCT L161-163)
        self.builder.post_treat();
    }

    /// OCCT BOPAlgo_Section::BuildSection (Section.cxx L167-414).
    fn build_section(&mut self) {
        // OCCT L169: Message_ProgressScope aPS(...) — no rcad plumbing (file
        // header difference #4).
        // OCCT L170-181: local declarations (aRC, aRC1, aBB, aExp, aLSA, aLS,
        // iterators, aMSI, aMS, aMFence, aMVE). rcad: the BRep_Builder +
        // TopoDS compound pair becomes the local ResultPool; the OCCT maps
        // keyed by TopTools_ShapeMapHasher (TShape + Location, orientation
        // ignored) become IndexMap<(u64, u32), _> preserving the OCCT
        // IndexedMap/DataMap insertion order.
        let mut a_ls_a: Vec<Shape> = Vec::new();
        let mut a_ls: Vec<Shape> = Vec::new();
        let mut a_ms_i: IndexMap<(u64, u32), (Shape, i32)> = IndexMap::new(); // OCCT aMSI
        let mut a_ms: IndexMap<(u64, u32), Shape> = IndexMap::new(); // OCCT aMS
        let mut a_m_fence: HashSet<(u64, u32)> = HashSet::new(); // OCCT aMFence
        let mut a_rc = ResultPool::new(self.builder.ds.locations.clone());
        let ds: &'a DS = self.builder.ds;
        //
        // OCCT L183: GetReport()->Clear();
        self.builder.my_report.clear();
        //
        // OCCT L185: BOPTools_AlgoTools::MakeContainer(TopAbs_COMPOUND, aRC1);
        // rcad: the compound children are collected first (aBB.Add below),
        // the container TShape is wrapped where aRC1 is first used (step 4).
        let mut a_rc1_children: Vec<Shape> = Vec::new();
        //
        // 1. aRC1 (OCCT L187-241)
        let a_nb = ds.nb_source_shapes(); // OCCT L188: aNb = myDS->NbSourceShapes();
        for i in 0..a_nb {
            // OCCT L191: const BOPDS_ShapeInfo& aSI = myDS->ShapeInfo(i);
            let a_si = ds.shape_info(i);
            // OCCT L192-195: if (aSI.ShapeType() != TopAbs_FACE) { continue; }
            if a_si.shape_type() != ShapeType::Face {
                continue;
            }
            // OCCT L196-199: if (UserBreak(aPS)) { return; } — no rcad
            // progress plumbing (file header difference #4).
            //
            // OCCT L201: const BOPDS_FaceInfo& aFI = myDS->FaceInfo(i);
            let a_fi = ds.face_info(i);
            //
            // 1.1 Vertices that are section vertices (OCCT L203-211)
            for n_v in a_fi.vertices_sc() {
                let n_v = *n_v;
                // OCCT L208-209: nV = aItMI.Key(); aV = myDS->Shape(nV);
                let a_v = ds.shape(n_v).clone();
                // OCCT L210: aBB.Add(aRC1, aV);
                a_rc1_children.push(a_rc.add(&a_v));
            }
            //
            // 1.2 Vertices that are in a face (OCCT L213-228)
            for n_v in a_fi.vertices_in() {
                let n_v = *n_v;
                // OCCT L218-222: if (nV < 0) { continue; } — rcad DS indices
                // are usize; the negative OCCT index has no representation
                // (architecture difference).
                // OCCT L223: if (myDS->IsNewShape(nV) || myDS->HasInterf(nV))
                if ds.is_new_shape(n_v) || ds.has_interf_single(n_v) {
                    // OCCT L225-226: aV = myDS->Shape(nV); aBB.Add(aRC1, aV);
                    let a_v = ds.shape(n_v).clone();
                    a_rc1_children.push(a_rc.add(&a_v));
                }
            }
            //
            // 1.3 Section edges (OCCT L230-240)
            // OCCT L231: const NCollection_IndexedMap<handle<PB>>& aMPBSc =
            //   aFI.PaveBlocksSc(); — rcad stores PB pointer ids; resolve
            //   them through the DS pools (pb_from_ptr).
            for a_pb_key in a_fi.pave_blocks_sc() {
                let Some(a_pb) = ds.pb_from_ptr(*a_pb_key) else {
                    continue;
                };
                // OCCT L237: nE = aPB->Edge();
                let n_e = a_pb.read().edge();
                // OCCT has no guard (myDS->Shape(nE) assumes the edge was
                // built); rcad marks edge-less PBs with NO_EDGE = usize::MAX,
                // so the DS lookup is bounds-checked like builder.rs
                // loc_generated.
                if n_e >= ds.nb_shapes() {
                    continue;
                }
                // OCCT L238-239: aE = myDS->Shape(nE); aBB.Add(aRC1, aE);
                let a_e = ds.shape(n_e).clone();
                a_rc1_children.push(a_rc.add(&a_e));
            }
        }
        //
        // 2. Common blocks between an edge and a face (OCCT L243-269)
        // OCCT L244-245: aPBP = myDS->PaveBlocksPool(); rcad pool is a
        // HashMap keyed by the pool index — sorted keys reproduce the OCCT
        // NCollection_DynamicArray index order 0..Length.
        let mut a_pool_keys: Vec<usize> =
            ds.pave_blocks_pool().keys().copied().collect();
        a_pool_keys.sort_unstable();
        for a_i in &a_pool_keys {
            // OCCT L250: const NCollection_List<handle<PB>>& aLPB = aPBP(i);
            for a_pb in &ds.pave_blocks_pool()[a_i] {
                // OCCT L254-255: aPB = aItPB.Value(); aCB = myDS->CommonBlock(aPB);
                let Some(a_cb_idx) = ds.common_block(a_pb) else {
                    continue; // OCCT L256: if (!aCB.IsNull())
                };
                let a_cb = &ds.common_blocks[a_cb_idx];
                // OCCT L258-259: aLF = aCB->Faces(); aNbF = aLF.Extent();
                let a_nb_f = a_cb.faces().len();
                // OCCT L260: if (aNbF)
                if a_nb_f > 0 {
                    // OCCT L262: aPBR = aCB->PaveBlock1();
                    let Some(a_pbr) = a_cb.pave_block1() else {
                        continue;
                    };
                    // OCCT L263: nE = aPBR->Edge();
                    let n_e = a_pbr.read().edge();
                    // rcad NO_EDGE bounds check (same as 1.3).
                    if n_e >= ds.nb_shapes() {
                        continue;
                    }
                    // OCCT L264-265: aE = myDS->Shape(nE); aBB.Add(aRC1, aE);
                    let a_e = ds.shape(n_e).clone();
                    a_rc1_children.push(a_rc.add(&a_e));
                }
            }
        }
        //
        // OCCT L271-279: dedup myArguments into aLSA (aMFence).
        for a_sa in &self.builder.my_arguments {
            if a_m_fence.insert((a_sa.ptr_id(), a_sa.location)) {
                a_ls_a.push(a_sa.clone());
            }
        }
        //
        // OCCT L281: aMFence.Clear();
        a_m_fence.clear();
        //
        // 3. Treatment boundaries of arguments (OCCT L283-356)
        // 3.1 Set to treat => aLS (OCCT L285-316)
        for a_sa in &a_ls_a {
            a_ls.clear();
            a_ms.clear();
            a_m_fence.clear();
            // OCCT L299-307: aExp.Init(aSA, TopAbs_EDGE); ...
            for a_e in explorer(a_sa, ShapeType::Edge, ShapeType::Shape) {
                if a_m_fence.insert((a_e.ptr_id(), a_e.location)) {
                    a_ls.push(a_e);
                }
            }
            // OCCT L308-316: aExp.Init(aSA, TopAbs_VERTEX); ...
            for a_e in explorer(a_sa, ShapeType::Vertex, ShapeType::Shape) {
                if a_m_fence.insert((a_e.ptr_id(), a_e.location)) {
                    a_ls.push(a_e);
                }
            }
            //
            // 3.2 aMSI (OCCT L318-355)
            for a_s in &a_ls {
                // OCCT L324: if (myImages.IsBound(aS))
                if let Some(a_ls_im) = self
                    .builder
                    .my_images
                    .get((a_s.ptr_id(), a_s.location))
                {
                    // OCCT L326-333: aLSIm = myImages.Find(aS); for each
                    // aSIm: MapShapes(aSIm, TopAbs_VERTEX, aMS);
                    // MapShapes(aSIm, TopAbs_EDGE, aMS);
                    for a_s_im in a_ls_im {
                        map_shapes(a_s_im, ShapeType::Vertex, &mut a_ms);
                        map_shapes(a_s_im, ShapeType::Edge, &mut a_ms);
                    }
                } else {
                    // OCCT L335-339: MapShapes(aS, TopAbs_VERTEX, aMS);
                    // MapShapes(aS, TopAbs_EDGE, aMS);
                    map_shapes(a_s, ShapeType::Vertex, &mut a_ms);
                    map_shapes(a_s, ShapeType::Edge, &mut a_ms);
                }
            }
            // OCCT L342-355: count per-argument occurrences in aMSI.
            let a_nb_ms = a_ms.len();
            for a_idx in 0..a_nb_ms {
                let (&a_key, a_s) = a_ms.get_index(a_idx).unwrap();
                match a_ms_i.get_mut(&a_key) {
                    // OCCT L348-349: int& iCnt = aMSI.ChangeFromKey(aS); ++iCnt;
                    Some(entry) => entry.1 += 1,
                    // OCCT L353: aMSI.Add(aS, 1);
                    None => {
                        a_ms_i.insert(a_key, (a_s.clone(), 1));
                    }
                }
            }
        }
        //
        // OCCT L358-359: aMS.Clear(); aMFence.Clear();
        a_ms.clear();
        a_m_fence.clear();
        //
        // 4. Build the result (OCCT L361-411)
        let mut a_mve: IndexMap<(u64, u32), (Shape, Vec<Shape>)> =
            IndexMap::new(); // OCCT L362-363: aMVE
                             // OCCT L365: MapShapesAndAncestors(aRC1, TopAbs_VERTEX, TopAbs_EDGE, aMVE);
        // aRC1 is a standalone working compound (OCCT L185): it is never added
        // to myShape, so it is NOT pushed into the result pool — a flat-pool
        // wrapper would leave an unreachable top compound behind. The child
        // Shapes keep their result-pool indices (exploration walks the Arc
        // topology only).
        let a_rc1 = Shape {
            data: Arc::new(TShape::Compound(std::mem::take(&mut a_rc1_children))),
            index: usize::MAX,
            location: 0,
            orientation: Orientation::Forward,
        };
        map_shapes_and_ancestors(
            &a_rc1,
            ShapeType::Vertex,
            ShapeType::Edge,
            &mut a_mve,
        );
        //
        // OCCT L367-376: vertices/edges common to more than one argument.
        let a_nb_ms = a_ms_i.len();
        for a_idx in 0..a_nb_ms {
            let (&_a_key, a_entry) = a_ms_i.get_index(a_idx).unwrap();
            // OCCT L372: if (iCnt > 1)
            if a_entry.1 > 1 {
                // OCCT L374: MapShapesAndAncestors(aV, TopAbs_VERTEX,
                // TopAbs_EDGE, aMVE);
                map_shapes_and_ancestors(
                    &a_entry.0,
                    ShapeType::Vertex,
                    ShapeType::Edge,
                    &mut a_mve,
                );
            }
        }
        //
        // OCCT L378: BOPTools_AlgoTools::MakeContainer(TopAbs_COMPOUND, aRC);
        let mut a_rc_children: Vec<Shape> = Vec::new();
        //
        // OCCT L380-411
        let a_nb_ms = a_mve.len();
        for a_idx in 0..a_nb_ms {
            // OCCT L196/L289/L383: if (UserBreak(aPS)) { return; } — no rcad
            // progress plumbing (file header difference #4).
            let (&_a_key, a_entry) = a_mve.get_index(a_idx).unwrap();
            let a_nb_le = a_entry.1.len();
            // OCCT L390-397
            if a_nb_le == 0 {
                // alone vertices
                if a_m_fence.insert((a_entry.0.ptr_id(), a_entry.0.location)) {
                    // OCCT L395: aBB.Add(aRC, aV);
                    a_rc_children.push(a_rc.add(&a_entry.0));
                }
            } else {
                // edges (OCCT L398-410)
                for a_e in &a_entry.1 {
                    if a_m_fence.insert((a_e.ptr_id(), a_e.location)) {
                        // OCCT L407: aBB.Add(aRC, aE);
                        a_rc_children.push(a_rc.add(a_e));
                    }
                }
            }
        }
        //
        // OCCT L413: myShape = aRC;
        let _ = a_rc.make_container(std::mem::take(&mut a_rc_children));
        self.builder.my_shape = Some(a_rc.brep);
    }

    /// OCCT BOPAlgo_BuilderShape::Shape — the result accessor.
    ///
    /// Normalizes the locations pool to the BRep convention (drop the leading
    /// identity entry) exactly like Builder::build_with_history_topods
    /// (builder.rs L1786-1788). The BRep-convention pcurve materialization of
    /// build_with_history_topods is not needed here: the section result
    /// contains no faces, so no pcurve rows can be re-keyed.
    pub fn result_brep(&self) -> Option<topods::BRep> {
        let mut brep = self.builder.my_shape.clone()?;
        if brep.locations.first() == Some(&DAffine3::IDENTITY) {
            brep.locations.remove(0);
        }
        Some(brep)
    }
}

// OCCT TopAbs_ShapeEnum ordinal (TopAbs_ShapeEnum.hxx): COMPOUND=0,
// COMPSOLID=1, SOLID=2, SHELL=3, FACE=4, WIRE=5, EDGE=6, VERTEX=7, SHAPE=8.
// The rcad ShapeType declaration order is reversed, so the explorer's
// "more complex than" comparison (TopExp_Explorer.cxx L31-35) uses this
// explicit mapping.
fn topabs_value(t: ShapeType) -> i32 {
    match t {
        ShapeType::Compound => 0,
        ShapeType::CompSolid => 1,
        ShapeType::Solid => 2,
        ShapeType::Shell => 3,
        ShapeType::Face => 4,
        ShapeType::Wire => 5,
        ShapeType::Edge => 6,
        ShapeType::Vertex => 7,
        ShapeType::Shape => 8,
    }
}

/// OCCT TopoDS_Iterator (cumOri=true): the stored sub-shapes with
/// TopAbs::Compose(parent.Orientation(), child.Orientation()) applied
/// (TopoDS_Iterator.cxx L72-80). The child order is the rcad container field
/// order used by builder.rs push_shape_recursive; location composition is
/// implicit (rcad Shape nodes carry absolute location indices).
fn sub_shapes(sh: &Shape) -> Vec<Shape> {
    let composed = |c: &Shape| {
        let mut c = c.clone();
        c.orientation = sh.orientation.compose(c.orientation);
        c
    };
    match sh.data.as_ref() {
        TShape::Vertex(_) => Vec::new(),
        TShape::Edge(ed) => vec![composed(&ed.first), composed(&ed.last)],
        TShape::Wire(wd) => wd.edges.iter().map(composed).collect(),
        TShape::Face(fd) => {
            let mut out = Vec::with_capacity(
                1 + fd.inner_wires.len() + fd.internal_vertices.len(),
            );
            out.push(composed(&fd.outer_wire));
            out.extend(fd.inner_wires.iter().map(composed));
            out.extend(fd.internal_vertices.iter().map(composed));
            out
        }
        TShape::Shell(sd) => sd.faces.iter().map(composed).collect(),
        TShape::Solid(sd) => {
            let mut out = Vec::new();
            out.extend(sd.shells.iter().map(composed));
            out.extend(sd.internal_vertices.iter().map(composed));
            out.extend(sd.internal_edges.iter().map(composed));
            out
        }
        TShape::CompSolid(cs) => cs.iter().map(composed).collect(),
        TShape::Compound(cd) => cd.iter().map(composed).collect(),
    }
}

/// OCCT TopExp_Explorer (TopExp_Explorer.cxx L68-171): yields the sub-shapes
/// of `s` of type `to_find`, descending into coarser container types, never
/// entering types of `to_avoid` (shouldAvoid: toAvoid != TopAbs_SHAPE, hence
/// ShapeType::Shape as the "avoid nothing" sentinel). The root itself is
/// yielded when its type matches (Init L96-104), even when toAvoid equals
/// toFind; a non-matching root of the avoided type yields nothing (Next
/// L121-125).
fn explorer(s: &Shape, to_find: ShapeType, to_avoid: ShapeType) -> Vec<Shape> {
    let mut out: Vec<Shape> = Vec::new();
    // OCCT Init L84-87: if (toFind == TopAbs_SHAPE) { hasMore = false; }
    if to_find == ShapeType::Shape {
        return out;
    }
    // OCCT Init L90-95: if (ty > toFind) { hasMore = false; }
    if topabs_value(s.shape_type()) > topabs_value(to_find) {
        return out;
    }
    explore_rec(s, to_find, to_avoid, &mut out);
    out
}

/// The iterative walk of TopExp_Explorer::Next (TopExp_Explorer.cxx
/// L110-171) in recursive form: yield on isSameType (L145-149), descend on
/// isMoreComplex && !shouldAvoid (L150-152), skip otherwise (L153-158).
fn explore_rec(
    sh: &Shape,
    to_find: ShapeType,
    to_avoid: ShapeType,
    out: &mut Vec<Shape>,
) {
    let ty = sh.shape_type();
    if ty == to_find {
        out.push(sh.clone());
        return;
    }
    if topabs_value(ty) > topabs_value(to_find) {
        return;
    }
    // shouldAvoid (TopExp_Explorer.cxx L26-29): toAvoid != TopAbs_SHAPE.
    if to_avoid != ShapeType::Shape && ty == to_avoid {
        return;
    }
    for c in sub_shapes(sh) {
        explore_rec(&c, to_find, to_avoid, out);
    }
}

/// OCCT TopExp::MapShapes (TopExp.cxx L35-45): explore S for T and add every
/// yielded shape to the indexed map (TopTools_ShapeMapHasher identity:
/// TShape + Location, orientation ignored — existing entries are kept).
fn map_shapes(
    s: &Shape,
    t: ShapeType,
    m: &mut IndexMap<(u64, u32), Shape>,
) {
    for a_c in explorer(s, t, ShapeType::Shape) {
        m.entry((a_c.ptr_id(), a_c.location)).or_insert(a_c);
    }
}

/// OCCT TopExp::MapShapesAndAncestors (TopExp.cxx L80-120).
fn map_shapes_and_ancestors(
    s: &Shape,
    ts: ShapeType,
    ta: ShapeType,
    m: &mut IndexMap<(u64, u32), (Shape, Vec<Shape>)>,
) {
    // OCCT L87: NCollection_List<TopoDS_Shape> empty;
    // OCCT L89-107: visit ancestors — exa(S, TA); for each ancestor, exs(anc,
    // TS): find or add the sub-shape, append the ancestor.
    for a_anc in explorer(s, ta, ShapeType::Shape) {
        for a_exs in explorer(&a_anc, ts, ShapeType::Shape) {
            let key = (a_exs.ptr_id(), a_exs.location);
            let entry = m.entry(key).or_insert((a_exs.clone(), Vec::new()));
            entry.1.push(a_anc.clone());
        }
    }
    // OCCT L109-119: visit shapes not under ancestors — ex(S, TS, TA).
    for a_ex in explorer(s, ts, ta) {
        let key = (a_ex.ptr_id(), a_ex.location);
        m.entry(key).or_insert((a_ex, Vec::new()));
    }
}

/// rcad result pool for BuildSection — the local equivalent of the
/// BRep_Builder + TopoDS compound pair of OCCT (architecture glue, no direct
/// OCCT function). `add` mirrors builder.rs `Builder::push_shape_recursive`
/// (private there): sub-shapes pushed first, TShape Arc shared (OCCT
/// BRep_Builder::Add references the source TShape, TopoDS_Builder.cxx
/// L57-59), dedup by TShape identity, and the TShape-internal reference
/// indices re-pointed in place to this pool (single-threaded pipeline; the
/// DS is not read after the section result has left the builder).
struct ResultPool {
    brep: topods::BRep,
    remap: HashMap<u64, usize>,
}

impl ResultPool {
    /// The locations pool mirrors the DS pool 1:1 (same convention as
    /// Builder::prepare, builder.rs L1977-1983).
    fn new(ds_locations: Vec<DAffine3>) -> Self {
        let mut brep = topods::BRep::default();
        brep.locations = ds_locations;
        ResultPool {
            brep,
            remap: HashMap::new(),
        }
    }

    /// Mirror of builder.rs Builder::push_shape_recursive.
    fn push_recursive(&mut self, shape: &Shape) -> usize {
        let ptr = shape.ptr_id();
        if let Some(&idx) = self.remap.get(&ptr) {
            return idx;
        }
        // Reserve a slot in tshapes, sharing the source TShape Arc.
        let new_idx = self.brep.tshapes.len();
        self.brep.tshapes.push(shape.data.clone());
        self.remap.insert(ptr, new_idx);
        // Recursively push the sub-shapes first so their result indices exist.
        match shape.data.as_ref() {
            TShape::Edge(ed) => {
                let _ = self.push_recursive(&ed.first);
                let _ = self.push_recursive(&ed.last);
            }
            TShape::Wire(wd) => {
                for e in &wd.edges {
                    let _ = self.push_recursive(e);
                }
            }
            TShape::Face(fd) => {
                let _ = self.push_recursive(&fd.outer_wire);
                for w in &fd.inner_wires {
                    let _ = self.push_recursive(w);
                }
                for v in &fd.internal_vertices {
                    let _ = self.push_recursive(v);
                }
            }
            TShape::Shell(sd) => {
                for f in &sd.faces {
                    let _ = self.push_recursive(f);
                }
            }
            TShape::Solid(sd) => {
                for s in &sd.shells {
                    let _ = self.push_recursive(s);
                }
                for v in &sd.internal_vertices {
                    let _ = self.push_recursive(v);
                }
                for e in &sd.internal_edges {
                    let _ = self.push_recursive(e);
                }
            }
            TShape::CompSolid(shapes) => {
                for s in shapes {
                    let _ = self.push_recursive(s);
                }
            }
            TShape::Compound(shapes) => {
                for s in shapes {
                    let _ = self.push_recursive(s);
                }
            }
            TShape::Vertex(_) => {}
        }
        // Re-point the TShape-internal reference indices from the DS pool
        // positions to this pool's positions, in place on the shared TShape.
        let raw = Arc::as_ptr(&shape.data) as *mut TShape;
        // SAFETY: single-threaded build; no other &TShape borrow of this Arc
        // is alive at this point (same shared-mutation model as builder.rs).
        unsafe {
            match &mut *raw {
                TShape::Vertex(vd) => {
                    for s in vd.my_shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                }
                TShape::Edge(ed) => {
                    for s in ed.my_shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    if let Some(&i) = self.remap.get(&ed.first.ptr_id()) {
                        ed.first.index = i;
                    }
                    if let Some(&i) = self.remap.get(&ed.last.ptr_id()) {
                        ed.last.index = i;
                    }
                }
                TShape::Wire(wd) => {
                    for s in wd.my_shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    for e in wd.edges.iter_mut() {
                        if let Some(&i) = self.remap.get(&e.ptr_id()) {
                            e.index = i;
                        }
                    }
                }
                TShape::Face(fd) => {
                    for s in fd.my_shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    if let Some(&i) = self.remap.get(&fd.outer_wire.ptr_id()) {
                        fd.outer_wire.index = i;
                    }
                    for w in fd.inner_wires.iter_mut() {
                        if let Some(&i) = self.remap.get(&w.ptr_id()) {
                            w.index = i;
                        }
                    }
                    for v in fd.internal_vertices.iter_mut() {
                        if let Some(&i) = self.remap.get(&v.ptr_id()) {
                            v.index = i;
                        }
                    }
                }
                TShape::Shell(sd) => {
                    for s in sd.my_shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    for f in sd.faces.iter_mut() {
                        if let Some(&i) = self.remap.get(&f.ptr_id()) {
                            f.index = i;
                        }
                    }
                }
                TShape::Solid(sd) => {
                    for s in sd.my_shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                    for sh in sd.shells.iter_mut() {
                        if let Some(&i) = self.remap.get(&sh.ptr_id()) {
                            sh.index = i;
                        }
                    }
                    for v in sd.internal_vertices.iter_mut() {
                        if let Some(&i) = self.remap.get(&v.ptr_id()) {
                            v.index = i;
                        }
                    }
                    for e in sd.internal_edges.iter_mut() {
                        if let Some(&i) = self.remap.get(&e.ptr_id()) {
                            e.index = i;
                        }
                    }
                }
                TShape::CompSolid(shapes) => {
                    for s in shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                }
                TShape::Compound(shapes) => {
                    for s in shapes.iter_mut() {
                        if let Some(&i) = self.remap.get(&s.ptr_id()) {
                            s.index = i;
                        }
                    }
                }
            }
        }
        new_idx
    }

    /// OCCT aBB.Add(container, theShape) — returns the pooled child Shape
    /// (the index points into this pool).
    fn add(&mut self, shape: &Shape) -> Shape {
        let idx = self.push_recursive(shape);
        let mut s = shape.clone();
        s.index = idx;
        s
    }

    /// BOPTools_AlgoTools::MakeContainer(TopAbs_COMPOUND, theRC) — wraps the
    /// accumulated children into a compound of this pool.
    fn make_container(&mut self, children: Vec<Shape>) -> Shape {
        let ts = Arc::new(TShape::Compound(children));
        let idx = self.brep.tshapes.len();
        self.brep.tshapes.push(ts.clone());
        Shape {
            data: ts,
            index: idx,
            location: 0,
            orientation: Orientation::Forward,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bop::brep_algo_api::section;
    use glam::DVec3;
    use rcad_kernel::geom::Curve3;

    /// Count (vertices, edges, compounds) reachable from the result root —
    /// the `nbshapes` walk (shapes unreachable from myShape, like the aRC1
    /// working container that shares the flat pool, are not counted; dedup by
    /// TShape identity, mirroring nbshapes without -t / builder.rs).
    fn count_result(brep: &topods::BRep) -> (usize, usize, usize) {
        let (mut v, mut e, mut c) = (0usize, 0usize, 0usize);
        // Root = the last compound in the pool (BuildSection's aRC).
        let Some(root) = brep.tshapes.iter().enumerate().rev().find_map(
            |(i, ts)| match ts.as_ref() {
                TShape::Compound(_) => Some(i),
                _ => None,
            },
        ) else {
            return (0, 0, 0);
        };
        let mut seen: HashSet<u64> = HashSet::new();
        let mut stack: Vec<Shape> = vec![Shape {
            data: brep.tshapes[root].clone(),
            index: root,
            location: 0,
            orientation: Orientation::Forward,
        }];
        while let Some(sh) = stack.pop() {
            if !seen.insert(sh.ptr_id()) {
                continue;
            }
            match sh.data.as_ref() {
                TShape::Vertex(_) => v += 1,
                TShape::Edge(ed) => {
                    e += 1;
                    stack.push(ed.first.clone());
                    stack.push(ed.last.clone());
                }
                TShape::Compound(cd) => {
                    c += 1;
                    for x in cd {
                        stack.push(x.clone());
                    }
                }
                TShape::Wire(wd) => {
                    for x in &wd.edges {
                        stack.push(x.clone());
                    }
                }
                _ => {}
            }
        }
        (v, e, c)
    }

    /// 3D curve types of the section edges, in pool order.
    fn edge_curve_kinds(brep: &topods::BRep) -> Vec<Option<&'static str>> {
        brep.tshapes
            .iter()
            .filter_map(|ts| match ts.as_ref() {
                TShape::Edge(ed) => Some(ed.curve.as_ref().map(|c| match c {
                    Curve3::Line(_) => "line",
                    Curve3::Circle(_) => "circle",
                    Curve3::Ellipse(_) => "ellipse",
                    Curve3::BSpline(_) => "bspline",
                    Curve3::Bezier(_) => "bezier",
                    _ => "other",
                })),
                _ => None,
            })
            .collect()
    }

    /// Anchor 1 — OCCT reference (DRAWEXE):
    ///   box b1 2 2 2
    ///   box b2 1 1 1 2 2 2
    ///   bsection ss b1 b2 -n2d -na
    ///   nbshapes ss -> VERTEX 6, EDGE 6, WIRE 0, FACE 0, COMPOUND 1
    /// The section of the two partially overlapping boxes is the closed
    /// hexagonal loop where the surfaces of the two boxes cross.
    #[test]
    fn section_of_two_boxes_is_closed_hexagon() {
        let b1 = rcad_modeling::make_box_brep(
            DVec3::ZERO,
            DVec3::X,
            DVec3::Y,
            2.0,
            2.0,
            2.0,
        )
        .expect("box1");
        let b2 = rcad_modeling::make_box_brep(
            DVec3::new(1.0, 1.0, 1.0),
            DVec3::X,
            DVec3::Y,
            2.0,
            2.0,
            2.0,
        )
        .expect("box2");
        let result = section(&b1, &b2).expect("section of two boxes");
        let (v, e, c) = count_result(&result);
        assert_eq!(
            (v, e, c),
            (6, 6, 1),
            "OCCT reference: VERTEX 6, EDGE 6, COMPOUND 1"
        );
        // plane x plane section curves are 3D lines.
        let kinds = edge_curve_kinds(&result);
        assert_eq!(kinds.len(), 6, "6 section edges");
        assert!(
            kinds.iter().all(|k| *k == Some("line")),
            "all section edges must be lines, got {kinds:?}"
        );
    }

    /// Anchor 2 — OCCT reference (DRAWEXE):
    ///   box b1 -2 -2 0 4 4 2
    ///   pcylinder c 1 4
    ///   bsection ss b1 c -n2d -na
    ///   nbshapes ss -> VERTEX 2, EDGE 2, WIRE 0, FACE 0, COMPOUND 1
    /// The cylinder (r=1, z in [0,4]) pierces the box bottom (z=0) and top
    /// (z=2) faces; each piercing curve is a closed circle edge carrying its
    /// seam vertex.
    #[test]
    fn section_of_box_and_cylinder_is_two_circles() {
        let b1 = rcad_modeling::make_box_brep(
            DVec3::new(-2.0, -2.0, 0.0),
            DVec3::X,
            DVec3::Y,
            4.0,
            4.0,
            2.0,
        )
        .expect("box");
        let cyl =
            rcad_modeling::make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 4.0)
                .expect("cylinder");
        let result = section(&b1, &cyl).expect("section of box and cylinder");
        let (v, e, c) = count_result(&result);
        assert_eq!(
            (v, e, c),
            (2, 2, 1),
            "OCCT reference: VERTEX 2, EDGE 2, COMPOUND 1"
        );
        // plane x cylinder section curves are circles of radius 1.
        let kinds = edge_curve_kinds(&result);
        assert_eq!(kinds.len(), 2, "2 section edges");
        let mut circles = 0usize;
        for ts in &result.tshapes {
            if let TShape::Edge(ed) = ts.as_ref() {
                if let Some(Curve3::Circle(circ)) = ed.curve.as_ref() {
                    circles += 1;
                    assert!(
                        (circ.radius - 1.0).abs() < 1.0e-6,
                        "circle radius must be 1, got {}",
                        circ.radius
                    );
                }
            }
        }
        assert_eq!(circles, 2, "both section edges must be circles");
    }
}
