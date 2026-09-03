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

use rcad_kernel::topo::topods::Shape;

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
