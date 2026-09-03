//! OCCT TopOpeBRepBuild package (TKBool/TopOpeBRepBuild) — 1:1 translation
//! of the reconstruction builder used by the ChFi3d Compute flow.
//!
//! Sources (recursion discipline: top-down, one layer per unit, pending
//! recursion points carry their OCCT file:line):
//!   - TopOpeBRepBuild_Builder.hxx (fields L855-910)
//!   - TopOpeBRepBuild_Builder.cxx (Clear L174-215, Perform L116-133)
//!   - TopOpeBRepBuild_Merge.cxx (MergeShapes L174-338, MergeSolids
//!     L356-364, MergeSolid L366-371)
//!
//! Pending recursion (next units, in call order from merge_shapes):
//!   - MapShapes / SplitSectionEdges / IsKPart / MergeKPart / Reverse
//!     (TopOpeBRepBuild_Merge.cxx / Builder1.cxx)
//!   - SplitShapes -> TopOpeBRepBuild_ShellFaceSet / ShapeSet /
//!     AreaBuilder / BlockBuilder / FaceBuilder / WireBuilder chain
//!     (TopOpeBRepBuild_Builder1*.cxx, ~7k lines)
//!   - TopOpeBRepDS_BuildTool / TopOpeBRepDS_Filter / Reducer
//!     (the Perform() pre-processing steps)

use std::collections::HashMap;

use rcad_kernel::geom::{CurveEval as _, SurfaceEval as _};
use rcad_kernel::topo::topods::{BRepTool as _, Orientation, Shape};

use super::topopebrepds::TopOpeBRepDSHDataStructure;

// =========================================================================
// OCCT TopAbs_State (TopAbs_State.hxx).
// =========================================================================
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopAbsState {
    In,
    Out,
    On,
    Unknown,
}

// =========================================================================
// OCCT TopOpeBRepDS_ListOfShapeOn1State (TopOpeBRepDS_ListOfShapeOn1State.hxx)
// — a shape list tagged with the state it was merged/split for.
// =========================================================================
#[derive(Debug, Clone, Default)]
pub struct ListOfShapeOn1State {
    pub has_shape: bool,
    pub state: Option<TopAbsState>,
    pub list: Vec<Shape>,
}

impl ListOfShapeOn1State {
    pub fn set_shape(&mut self, s: TopAbsState) {
        self.has_shape = true;
        self.state = Some(s);
    }

    pub fn clear(&mut self) {
        self.has_shape = false;
        self.state = None;
        self.list.clear();
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Builder (Builder.hxx L855-910 private fields).
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepBuildBuilder {
    /// OCCT: occ::handle<TopOpeBRepDS_HDataStructure> myDataStructure.
    pub my_data_structure: Option<TopOpeBRepDSHDataStructure>,
    /// OCCT: TopOpeBRepDS_BuildTool myBuildTool (the topo-building facade;
    /// its operations map onto the rcad BRep builder — pending).
    pub my_build_tool: TopOpeBRepDSBuildTool,
    /// OCCT: TopOpeBRepDS_DataMapOfShapeListOfShapeOn1State mySplitIN/ON/OUT.
    pub my_split_in: HashMap<u64, ListOfShapeOn1State>,
    pub my_split_on: HashMap<u64, ListOfShapeOn1State>,
    pub my_split_out: HashMap<u64, ListOfShapeOn1State>,
    /// OCCT: myMergedIN / myMergedON / myMergedOUT.
    pub my_merged_in: HashMap<u64, ListOfShapeOn1State>,
    pub my_merged_on: HashMap<u64, ListOfShapeOn1State>,
    pub my_merged_out: HashMap<u64, ListOfShapeOn1State>,
    /// OCCT: TopAbs_State myState1 / myState2; TopoDS_Shape myShape1/2.
    pub my_state1: Option<TopAbsState>,
    pub my_state2: Option<TopAbsState>,
    pub my_shape1: Shape,
    pub my_shape2: Shape,
    /// OCCT: int myIsKPart; bool myClassifyDef / myClassifyVal.
    pub my_is_kpart: i32,
    pub my_classify_def: bool,
    pub my_classify_val: bool,
    /// OCCT: TopTools_IndexedMapOfShape myMAP1 / myMAP2 (Builder.hxx).
    pub my_map1: IndexedShapeMap,
    pub my_map2: IndexedShapeMap,
    /// OCCT: the 2d-pur face list (empty in the fillet 3d flow).
    pub my_list_of_face: Vec<Shape>,
    /// OCCT: static int STATIC_SOLIDINDEX (Builder.cxx file static; set by
    /// SplitSolid — pending unit).
    pub static_solid_index: i32,
    /// rcad architecture: the BRep the shapes live in (OCCT reads geometry
    /// from shape handles).
    pub build_brep: rcad_kernel::topods::BRep,
}

impl Default for TopOpeBRepBuildBuilder {
    fn default() -> Self {
        TopOpeBRepBuildBuilder {
            my_data_structure: None,
            my_build_tool: TopOpeBRepDSBuildTool,
            my_split_in: HashMap::new(),
            my_split_on: HashMap::new(),
            my_split_out: HashMap::new(),
            my_merged_in: HashMap::new(),
            my_merged_on: HashMap::new(),
            my_merged_out: HashMap::new(),
            my_state1: None,
            my_state2: None,
            my_shape1: Shape::null(),
            my_shape2: Shape::null(),
            my_is_kpart: 0,
            my_classify_def: false,
            my_classify_val: false,
            my_map1: IndexedShapeMap::default(),
            my_map2: IndexedShapeMap::default(),
            my_list_of_face: Vec::new(),
            static_solid_index: 1,
            build_brep: rcad_kernel::topods::BRep::default(),
        }
    }
}

/// OCCT TopOpeBRepDS_BuildTool — the facade over BRep_Builder topology
/// construction (MakeSolid/MakeShell/CopyFace/UpdateSurface/...).  The
/// individual operations are translated with the BRep builder units.
#[derive(Debug, Clone, Default)]
pub struct TopOpeBRepDSBuildTool;

impl TopOpeBRepBuildBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    // OCCT accessors used by the ChFi3d call surface (Builder.hxx):
    //   Merged(S, ToBuild) / IsSplit(S, ToBuild) / Splits(S, ToBuild) /
    //   NewFaces(I) — translated alongside the maps above.

    /// OCCT Builder.hxx Merged(S, ToBuild) — read accessor.
    pub fn merged(&self, s: &Shape, tobuild: TopAbsState) -> &[Shape] {
        let m = match tobuild {
            TopAbsState::In => &self.my_merged_in,
            TopAbsState::On => &self.my_merged_on,
            _ => &self.my_merged_out,
        };
        m.get(&s.ptr_id())
            .map(|l| l.list.as_slice())
            .unwrap_or(&[])
    }

    /// OCCT Builder.hxx IsSplit(S, ToBuild).
    pub fn is_split(&self, s: &Shape, tobuild: TopAbsState) -> bool {
        let m = match tobuild {
            TopAbsState::In => &self.my_split_in,
            TopAbsState::On => &self.my_split_on,
            _ => &self.my_split_out,
        };
        m.get(&s.ptr_id()).is_some_and(|l| !l.list.is_empty())
    }

    /// OCCT Builder.hxx Splits(S, ToBuild) — read accessor.
    pub fn splits(&self, s: &Shape, tobuild: TopAbsState) -> &[Shape] {
        let m = match tobuild {
            TopAbsState::In => &self.my_split_in,
            TopAbsState::On => &self.my_split_on,
            _ => &self.my_split_out,
        };
        m.get(&s.ptr_id())
            .map(|l| l.list.as_slice())
            .unwrap_or(&[])
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L174-215 — Clear().
    // =========================================================================
    pub fn clear(&mut self) {
        if self.my_data_structure.is_none() {
            self.my_merged_out.clear();
            self.my_merged_in.clear();
            self.my_merged_on.clear();
            return;
        }
        // OCCT L194-215: the split maps of EDGE shapes with no geometry in
        // the DS are cleared; the whole DataMap iteration relies on the DS
        // shape info (translated with the Filter/Reducer unit).
        self.my_split_out.clear();
        self.my_split_on.clear();
        self.my_split_in.clear();
        self.my_merged_out.clear();
        self.my_merged_in.clear();
        self.my_merged_on.clear();
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L116-133 — Perform(HDS).
    // =========================================================================
    pub fn perform(&mut self, hds: TopOpeBRepDSHDataStructure) {
        self.clear();
        self.my_data_structure = Some(hds);
        // BuildVertices(HDS);            (BuildVertices.cxx — pending unit)
        // SplitEvisoONperiodicF();       (Builder1_2.cxx — pending unit)
        // BuildEdges(HDS);               (BuildEdges.cxx — pending unit)
        // BuildFaces(HDS);               (BuildFaces.cxx — pending unit)
        self.my_is_kpart = 0;
        // InitSection(); SplitSectionEdges();  (Builder1.cxx — pending unit)
        // TopOpeBRepDS_Filter F(HDS, &myShapeClassifier);
        //   F.ProcessFaceInterferences(mySplitON);   (pending unit)
        // TopOpeBRepDS_Reducer R(HDS);
        //   R.ProcessFaceInterferences(mySplitON);   (pending unit)
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Merge.cxx L356-371 — MergeSolids / MergeSolid.
    // =========================================================================
    pub fn merge_solids(&mut self, s1: &Shape, tobuild1: TopAbsState, s2: &Shape, tobuild2: TopAbsState) {
        self.merge_shapes(s1, tobuild1, s2, tobuild2);
    }

    pub fn merge_solid(&mut self, s: &Shape, tobuild: TopAbsState) {
        let snull = Shape::null();
        self.merge_shapes(s, tobuild, &snull, tobuild);
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Merge.cxx L174-338 — MergeShapes.
    // =========================================================================
    pub fn merge_shapes(&mut self, s1: &Shape, tobuild1: TopAbsState, s2: &Shape, tobuild2: TopAbsState) {
        // OCCT L177: lesmemes = S1.IsEqual(S2) — same TShape + location +
        // orientation.
        let lesmemes = s1.is_same(s2) && s1.orientation == s2.orientation;
        if lesmemes {
            return;
        }

        self.my_state1 = Some(tobuild1);
        self.my_state2 = Some(tobuild2);
        self.my_shape1 = s1.clone();
        self.my_shape2 = s2.clone();
        let s1null = s1.is_null();
        let s2null = s2.is_null();

        // MapShapes(S1, S2);       (pending unit — Builder1.cxx)
        // SplitSectionEdges();     (pending unit)
        //======================== debut KPart
        // if (IsKPart()) { MergeKPart(); ClearMaps(); return; }
        // (the fillet flow calls Perform(HDS) so myIsKPart stays 0 and this
        // branch is not taken; translated with the MergeKPart unit.)
        //======================== fin KPart

        // bool RevOri1 = Reverse(ToBuild1, ToBuild2);
        // bool RevOri2 = Reverse(ToBuild2, ToBuild1);
        // TopOpeBRepBuild_ShellFaceSet SFS;
        // SplitShapes(ex1, ToBuild1, ToBuild2, SFS, RevOri1);
        // SplitShapes(ex2, ToBuild2, ToBuild1, SFS, RevOri2);
        // ... (SplitShapes + the ShellFaceSet/ShapeSet/AreaBuilder chain is
        // the next recursion unit, Builder1.cxx L~700+ / Merge.cxx L~400)
        let _ = (s1null, s2null);

        // ClearMaps();
    }
}


// =========================================================================
// OCCT TopAbs::Complement (TopAbs_Orientation.hxx).
// =========================================================================
fn topabs_complement(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Reversed,
        Orientation::Reversed => Orientation::Forward,
        other => other,
    }
}

/// OCCT TopTools_IndexedMapOfShape stand-in — insertion-ordered unique
/// shapes keyed by TShape pointer.
#[derive(Debug, Clone, Default)]
pub struct IndexedShapeMap {
    pub keys: Vec<Shape>,
    pub index: HashMap<u64, usize>,
}

impl IndexedShapeMap {
    pub fn add(&mut self, s: &Shape) {
        if !self.index.contains_key(&s.ptr_id()) {
            self.keys.push(s.clone());
            self.index.insert(s.ptr_id(), self.keys.len());
        }
    }

    pub fn clear(&mut self) {
        self.keys.clear();
        self.index.clear();
    }
}

impl TopOpeBRepBuildBuilder {
    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L621-632 — Reverse.
    // =========================================================================
    pub fn reverse_states(tobuild1: TopAbsState, tobuild2: TopAbsState) -> bool {
        if tobuild1 == TopAbsState::In && tobuild2 == TopAbsState::In {
            false
        } else {
            tobuild1 == TopAbsState::In
        }
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L637-640 — Orient.
    // =========================================================================
    pub fn orient(ori: Orientation, rev: bool) -> Orientation {
        if !rev {
            ori
        } else {
            topabs_complement(ori)
        }
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L766-780 — MapShapes(S1, S2):
    // ClearMaps then TopExp::MapShapes of each non-null shape into
    // myMAP1/myMAP2.
    // =========================================================================
    pub fn map_shapes(&mut self, s1: &Shape, s2: &Shape) {
        let s1null = s1.is_null();
        let s2null = s2.is_null();
        self.clear_maps();
        if !s1null {
            self.my_map1.add(s1);
        }
        if !s2null {
            self.my_map2.add(s2);
        }
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L783-787 — ClearMaps.
    // =========================================================================
    pub fn clear_maps(&mut self) {
        self.my_map1.clear();
        self.my_map2.clear();
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L1713-1908 — SplitShapes: walk the
    // explorer of one input shape, dispatch the Split* machinery, and feed
    // the resulting shapes into the ShapeSet as start elements (splits) or
    // elements (kept originals).
    // =========================================================================
    pub fn split_shapes(
        &mut self,
        shapes: &[Shape],
        tobuild1: TopAbsState,
        tobuild2: TopAbsState,
        start_elements: &mut Vec<Shape>,
        rev_ori: bool,
    ) {
        for a_shape in shapes {
            // compute new orientation <newori> to give to the new shapes.
            let newori = Self::orient(a_shape.orientation, rev_ori);
            let t = a_shape.shape_type();

            match t {
                rcad_kernel::topods::ShapeType::Solid | rcad_kernel::topods::ShapeType::Shell => {
                    self.split_solid(a_shape, tobuild1, tobuild2);
                }
                rcad_kernel::topods::ShapeType::Face => {
                    self.split_face(a_shape, tobuild1, tobuild2);
                }
                rcad_kernel::topods::ShapeType::Edge => {
                    self.split_edge(a_shape, tobuild1, tobuild2);
                }
                _ => {
                    continue;
                }
            }

            if self.is_split(a_shape, tobuild1) {
                let mut is_lson = false;
                let mut ls: Vec<Shape> = self.splits(a_shape, tobuild1).to_vec();
                if matches!(t, rcad_kernel::topods::ShapeType::Edge)
                    && tobuild1 == TopAbsState::In
                    && ls.is_empty()
                {
                    ls = self.splits(a_shape, TopAbsState::On).to_vec();
                    is_lson = true;
                }
                for new_shape in &ls {
                    let mut new_shape = new_shape.clone();
                    new_shape.orientation = newori;
                    if is_lson {
                        let mut add = true;
                        if !self.my_list_of_face.is_empty() {
                            // 2d pur: KeepShape (Builder1.cxx — pending unit).
                            add = true;
                        }
                        if add {
                            start_elements.push(new_shape);
                        }
                    } else {
                        start_elements.push(new_shape);
                    }
                }
            } else {
                // aShape n'a pas de devenir de split par ToBuild1: on
                // construit les parties ToBuild1 de aShape (de S1).
                let mut add = true;
                let isedge = matches!(t, rcad_kernel::topods::ShapeType::Edge);
                let ds = self.my_data_structure.as_ref().expect("DS");
                let hs = ds.shape_index.contains_key(&a_shape.ptr_id());
                let hg = ds.has_geometry(a_shape);

                let mut testkeep = isedge && hs && (!hg);

                // xpu010399 (USA60299): touched-vertex edges without DS
                // geometry also qualify (FUN_touched — Builder1_2.cxx,
                // pending unit; the stand-in reports untouched).
                let mut istouched = isedge && (!hs) && (!hg);
                if istouched {
                    istouched = false; // FUN_touched pending
                }
                testkeep = testkeep || istouched;

                if testkeep {
                    if !self.my_list_of_face.is_empty() {
                        // 2d pur: KeepShape — pending unit.
                        add = true;
                    } else {
                        // on classifie en solide uniquement si E est dans la
                        // DS et E a ete purgee de ses interfs car en bout.
                        let sol = if self.static_solid_index == 1 {
                            self.my_shape2.clone()
                        } else {
                            self.my_shape1.clone()
                        };
                        if !sol.is_null() {
                            let (c3d, crange) = {
                                let b = &self.build_brep;
                                let r: Option<(rcad_kernel::geom::Curve3, [f64; 2])> =
                                    b.edge_curve_world(a_shape);
                                match r {
                                    Some(r) => r,
                                    None => panic!(
                                        "Standard_ProgramError: SplitShapes no 3D curve on edge"
                                    ),
                                }
                            };
                            let (first, last) = (crange[0], crange[1]);
                            let tt: f64 = 0.127956477;
                            let par = (1.0 - tt) * first + tt * last;
                            let p3d = c3d.point_at(par);
                            let tol3d = rcad_kernel::core::precision::CONFUSION;
                            let mut scl =
                                crate::topalgo::brep_class3d::solid_classifier::SolidClassifier::from_shape(&sol);
                            scl.perform(p3d, tol3d);
                            let state = scl.state(); // 0=IN 1=OUT 2=ON
                            add = match (state, tobuild1) {
                                (0, TopAbsState::In) => true,
                                (1, TopAbsState::Out) => true,
                                (2, TopAbsState::On) => true,
                                _ => false,
                            };
                        } else {
                            // sol.IsNull
                            add = true;
                        }
                    }
                }
                if add {
                    let mut a_shape = a_shape.clone();
                    a_shape.orientation = newori;
                    start_elements.push(a_shape);
                }
            }
        } // Ex.More
    }

    // =========================================================================
    // OCCT SplitSolid / SplitFace / SplitEdge (TopOpeBRepBuild_Builder1.cxx)
    // — the split machinery; translated with the ShapeSet/FaceBuilder/
    // AreaBuilder/BlockBuilder units.  Until then no splits are registered,
    // matching a DS without split results.
    // =========================================================================
    fn split_solid(&mut self, _s: &Shape, _tobuild1: TopAbsState, _tobuild2: TopAbsState) {
        // Builder1.cxx SplitSolid — pending unit.
    }

    fn split_face(&mut self, _s: &Shape, _tobuild1: TopAbsState, _tobuild2: TopAbsState) {
        // Builder1.cxx SplitFace — pending unit.
    }

    fn split_edge(&mut self, _s: &Shape, _tobuild1: TopAbsState, _tobuild2: TopAbsState) {
        // Builder1.cxx SplitEdge — pending unit.
    }
}

impl TopOpeBRepBuildBuilder {
    // =========================================================================
    // OCCT Builder.hxx map accessors used by the split/merge flow:
    //   MarkSplit(S, ToBuild) / ChangeSplit(S, ToBuild) /
    //   ChangeMerged(S, ToBuild).
    // =========================================================================
    fn split_map_mut(&mut self, tobuild: TopAbsState) -> &mut HashMap<u64, ListOfShapeOn1State> {
        match tobuild {
            TopAbsState::In => &mut self.my_split_in,
            TopAbsState::On => &mut self.my_split_on,
            _ => &mut self.my_split_out,
        }
    }

    fn merged_map_mut(&mut self, tobuild: TopAbsState) -> &mut HashMap<u64, ListOfShapeOn1State> {
        match tobuild {
            TopAbsState::In => &mut self.my_merged_in,
            TopAbsState::On => &mut self.my_merged_on,
            _ => &mut self.my_merged_out,
        }
    }

    /// OCCT Builder.hxx MarkSplit(S, ToBuild).
    fn mark_split(&mut self, s: &Shape, tobuild: TopAbsState) {
        self.split_map_mut(tobuild)
            .entry(s.ptr_id())
            .or_default()
            .set_shape(tobuild);
    }

    /// OCCT Builder.hxx ChangeSplit(S, ToBuild).
    fn change_split(&mut self, s: &Shape, tobuild: TopAbsState) -> &mut Vec<Shape> {
        &mut self
            .split_map_mut(tobuild)
            .entry(s.ptr_id())
            .or_default()
            .list
    }

    /// OCCT Builder.hxx ChangeMerged(S, ToBuild).
    fn change_merged(&mut self, s: &Shape, tobuild: TopAbsState) -> &mut Vec<Shape> {
        &mut self
            .merged_map_mut(tobuild)
            .entry(s.ptr_id())
            .or_default()
            .list
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L1155-1167 — SplitFace (the
    // SplitFace2 context flag is a debug-only variant: SplitFace1 runs).
    // =========================================================================
    pub fn split_face_public(&mut self, foriented: &Shape, tobuild1: TopAbsState, tobuild2: TopAbsState) {
        self.split_face1(foriented, tobuild1, tobuild2);
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L1171-1300 — SplitFace1.
    // =========================================================================
    pub fn split_face1(&mut self, foriented: &Shape, tobuild1: TopAbsState, tobuild2: TopAbsState) {
        // operation tobuild1 tobuild2  process face F  connect to 1  connect to 2
        // common    IN       IN        yes             yes           yes
        // fuse      OUT      OUT       yes             yes           yes
        // cut 1-2   OUT      IN        yes             yes           no
        // cut 2-1   IN       OUT       yes             yes           no
        let tosplit = self.to_split(foriented, tobuild1);
        if !tosplit {
            return;
        }

        let mut rev_ori1 = Self::reverse_states(tobuild1, tobuild2);
        let rev_ori2 = Self::reverse_states(tobuild2, tobuild1);
        let connect_to1 = true;
        let connect_to2 = false;

        // work on a FORWARD face <Fforward>.
        let mut fforward = foriented.clone();
        fforward.orientation = Orientation::Forward;

        // build the list of faces to split : LF1, LF2.
        let mut lf1: Vec<Shape> = vec![fforward.clone()];
        let mut lf2: Vec<Shape> = Vec::new();
        self.find_same_domain(&mut lf1, &mut lf2);
        let n1 = lf1.len();
        let n2 = lf2.len();

        if n2 == 0 {
            rev_ori1 = false;
        }
        if n1 == 0 {
            // OCCT sets RevOri2 = false here.
        }

        // Create an edge set <WES> connected by vertices
        // (TopOpeBRepBuild_WireEdgeSet — pending unit).
        let mut wes_start_elements: Vec<Shape> = Vec::new();

        for fcur in &lf1 {
            // FillFace(Fcur, ToBuild1, LF2, ToBuild2, WES, RevOri1)
            // (Builder1.cxx — pending unit).
            let _ = fcur;
        }
        for fcur in &lf2 {
            // FillFace(Fcur, ToBuild2, LF1, ToBuild1, WES, RevOri2)
            // (Builder1.cxx — pending unit).
            let _ = fcur;
        }

        // Add the intersection edges to edge set WES
        // (AddIntersectionEdges, Builder.cxx L143-171 — pending unit).
        let _ = (rev_ori1, &mut wes_start_elements);

        // Create a Face Builder FBU; FBU.InitFaceBuilder(WES, Fforward,
        // false) (TopOpeBRepBuild_FaceBuilder.cxx — pending unit).
        //
        // Build the new faces: MakeFaces(Fforward, FBU, FaceList)
        // (Builder.cxx L431+ — pending unit).
        let face_list: Vec<Shape> = Vec::new();
        self.change_merged(&fforward, tobuild1).extend(face_list.iter().cloned());

        // connect new faces as faces built <ToBuild1> on LF1 faces.
        for fcur in &lf1 {
            self.mark_split(fcur, tobuild1);
            if connect_to1 {
                let fl = self.change_split(fcur, tobuild1);
                *fl = face_list.clone();
            }
        }

        // connect new faces as faces built <ToBuild2> on LF2 faces.
        for fcur in &lf2 {
            self.mark_split(fcur, tobuild2);
            if connect_to2 {
                let fl = self.change_split(fcur, tobuild2);
                *fl = face_list.clone();
            }
        }
    }

    // =========================================================================
    // OCCT Builder.cxx ToSplit(S, ToBuild): a shape is to-split when the DS
    // holds split shapes for it (Builder.hxx inline: mySplitS has the shape
    // in the ToBuild map).  The DS-query core is translated with the
    // FillFace unit.
    // =========================================================================
    fn to_split(&self, _s: &Shape, _tobuild: TopAbsState) -> bool {
        // OCCT: return ToSplit(S, ToBuild) checks the DS split info;
        // pending the FillFace/DS-split unit.
        false
    }

    // =========================================================================
    // OCCT Builder.cxx FindSameDomain(L1, L2) (L650-764): complete L1/L2
    // with the DS same-domain shapes.  The rcad DS carries no same-domain
    // table yet (TopOpeBRepDS_AddShapeSameDomain — pending unit), so the
    // completion finds nothing, which matches a DS without same-domain
    // records.
    // =========================================================================
    fn find_same_domain(&self, _l1: &mut Vec<Shape>, _l2: &mut Vec<Shape>) {
        // OCCT L650-764 — pending the DS same-domain table unit.
    }
}
