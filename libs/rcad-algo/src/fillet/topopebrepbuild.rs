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

use super::chfi3d_builder_0::brep_tool_parameter;
use super::topopebrepds::{TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference};

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
                    self.split_edge_public(a_shape, tobuild1, tobuild2);
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

    // SplitEdge / SplitEdge1 are translated in the impl block at the end of
    // this file (Builder.cxx L925-1079).
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

impl TopOpeBRepBuildBuilder {
    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L925-937 — SplitEdge (the
    // SplitEdge2 context flag is a debug-only variant: SplitEdge1 runs).
    // =========================================================================
    pub fn split_edge_public(&mut self, e: &Shape, tobuild1: TopAbsState, tobuild2: TopAbsState) {
        self.split_edge1(e, tobuild1, tobuild2);
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L941-1079 — SplitEdge1.
    // =========================================================================
    pub fn split_edge1(&mut self, eoriented: &Shape, tobuild1: TopAbsState, tobuild2: TopAbsState) {
        // work on a FORWARD edge <Eforward>.
        let mut eforward = eoriented.clone();
        eforward.orientation = Orientation::Forward;

        let tosplit = self.to_split(eoriented, tobuild1);
        if !tosplit {
            return;
        }

        Self::reverse_states(tobuild1, tobuild2);
        Self::reverse_states(tobuild2, tobuild1);
        let connect_to1 = true;
        let connect_to2 = false;

        // build the list of edges to split : LE1, LE2.
        let mut le1: Vec<Shape> = vec![eforward.clone()];
        let mut le2: Vec<Shape> = Vec::new();
        self.find_same_domain(&mut le1, &mut le2);

        // Make a PaveSet <PVS> on edge <Eforward>
        // (TopOpeBRepBuild_PaveSet.cxx — pending unit).
        //
        // Add the points/vertices found on edge <Eforward> in <PVS>:
        // TopOpeBRepDS_PointIterator EPIT(myDataStructure->EdgePoints(...));
        // FillVertexSet(EPIT, ToBuild1, PVS) (Builder.cxx L~850 — pending
        // unit).  The DS edge-points list (point interferences on the DS
        // shape of the edge) is scanned to know whether any loop entry
        // exists.
        let ds = self.my_data_structure.as_ref().expect("DS");
        let edge_index = ds.shape_index.get(&eforward.ptr_id()).copied();
        let has_edge_points = match edge_index {
            Some(i) => ds.shape_interferences(i).iter().any(
                |it| matches!(it, TopOpeBRepDSInterference::CurvePoint(_)),
            ),
            None => false,
        };

        // TopOpeBRepBuild_PaveClassifier VCL(Eforward);
        // bool equalpar = PVS.HasEqualParameters();
        // if (equalpar) VCL.SetFirstParameter(PVS.EqualParameters());
        // (PaveSet/PaveClassifier — pending units.)

        // before return if PVS has no vertices, mark <Eforward> as split
        // <ToBuild1>.
        self.mark_split(&eforward, tobuild1);

        // PVS.InitLoop(); if (!PVS.MoreLoop()) return;  — the no-vertex
        // outcome leaves the edge marked split with an empty split list
        // (filled by the LE1 connection below).
        let _pvs_more_loop = has_edge_points;

        // build the new edges:
        //   TopOpeBRepBuild_EdgeBuilder EBU(PVS, VCL);
        //   MakeEdges(Eforward, EBU, EdgeList);
        // (EdgeBuilder.cxx / BuildEdges.cxx — pending units.)
        let edge_list: Vec<Shape> = self.merged(&eforward, tobuild1).to_vec();

        // connect new edges as edges built <ToBuild1> on LE1 edge.
        for ecur in &le1 {
            self.mark_split(ecur, tobuild1);
            if connect_to1 {
                let el = self.change_split(ecur, tobuild1);
                *el = edge_list.clone();
            }
        }

        // connect new edges as edges built <ToBuild2> on LE2 edges.
        for ecur in &le2 {
            self.mark_split(ecur, tobuild2);
            if connect_to2 {
                let el = self.change_split(ecur, tobuild2);
                *el = edge_list.clone();
            }
        }
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Pave (Pave.hxx L30-60 + Pave.cxx) — a vertex and
// its parameter on the edge, with the bound flag (old/new vertex).
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepBuildPave {
    vertex: Shape,
    parameter: f64,
    bound: bool,
    has_same_domain: bool,
    same_domain: Option<Shape>,
}

impl TopOpeBRepBuildPave {
    /// OCCT Pave.cxx L20-30 — Pave(V, P, bound).
    pub fn new(v: Shape, p: f64, bound: bool) -> Self {
        TopOpeBRepBuildPave {
            vertex: v,
            parameter: p,
            bound,
            has_same_domain: false,
            same_domain: None,
        }
    }

    pub fn vertex(&self) -> &Shape {
        &self.vertex
    }

    pub fn change_vertex(&mut self) -> &mut Shape {
        &mut self.vertex
    }

    pub fn parameter(&self) -> f64 {
        self.parameter
    }

    pub fn set_parameter(&mut self, par: f64) {
        self.parameter = par;
    }

    pub fn is_bound(&self) -> bool {
        self.bound
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_PaveSet (PaveSet.hxx L40-80 + PaveSet.cxx L35-490)
// — the ordered set of paves on one edge; the Loop iteration runs over the
// paves after Prepare() (the sort/clean pass — PaveSet.cxx L~200-490,
// pending unit; until then the insertion order is kept).
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepBuildPaveSet {
    pub edge: Shape,
    pub vertices: Vec<TopOpeBRepBuildPave>,
    pub has_equal_parameters: bool,
    pub equal_parameters: f64,
    loop_index: usize,
    /// OCCT: bool myPrepareDone / bool myRemovePV (true by default).
    prepare_done: bool,
    remove_pv: bool,
}

impl TopOpeBRepBuildPaveSet {
    /// OCCT PaveSet.cxx L35-45 — PaveSet(E).
    pub fn new(e: &Shape) -> Self {
        TopOpeBRepBuildPaveSet {
            edge: e.clone(),
            vertices: Vec::new(),
            has_equal_parameters: false,
            equal_parameters: 0.0,
            loop_index: 0,
            prepare_done: false,
            remove_pv: true,
        }
    }

    /// OCCT PaveSet.cxx L47-52 — Append(PV).
    pub fn append(&mut self, pv: TopOpeBRepBuildPave) {
        self.vertices.push(pv);
    }

    /// OCCT PaveSet.cxx InitLoop: Prepare() then iterate.
    pub fn init_loop(&mut self, brep: &rcad_kernel::topods::BRep) {
        if !self.prepare_done {
            self.prepare(brep);
        }
        self.loop_index = 0;
    }

    pub fn more_loop(&self) -> bool {
        self.loop_index < self.vertices.len()
    }

    pub fn loop_pave(&self) -> &TopOpeBRepBuildPave {
        &self.vertices[self.loop_index]
    }

    pub fn next_loop(&mut self) {
        self.loop_index += 1;
    }

    /// OCCT PaveSet.cxx HasEqualParameters: two distinct vertices with
    /// parameters closer than Precision::PConfusion().
    pub fn has_equal_parameters(&mut self) -> bool {
        self.has_equal_parameters = false;
        for i in 0..self.vertices.len() {
            let p1 = self.vertices[i].parameter();
            for j in 0..self.vertices.len() {
                if j == i {
                    continue;
                }
                if self.vertices[j].vertex.is_same(&self.vertices[i].vertex) {
                    continue;
                }
                let p2 = self.vertices[j].parameter();
                if (p1 - p2).abs() < 1.0e-9 /* OCCT Precision::PConfusion() = Confusion()*0.01 */ {
                    self.has_equal_parameters = true;
                    self.equal_parameters = p1;
                }
            }
        }
        self.has_equal_parameters
    }

    pub fn equal_parameters(&self) -> f64 {
        self.equal_parameters
    }
}

impl TopOpeBRepBuildBuilder {
    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L1991-1999 — FillVertexSet.
    // =========================================================================
    pub fn fill_vertex_set(
        &self,
        points: &[(&TopOpeBRepDSInterference, i32)],
        is_point_flags: &[bool],
        tobuild: TopAbsState,
        pvs: &mut TopOpeBRepBuildPaveSet,
    ) {
        for (i, (it, _ind)) in points.iter().enumerate() {
            self.fill_vertex_set_on_value(it, is_point_flags[i], tobuild, pvs);
        }
    }

    // =========================================================================
    // OCCT TopOpeBRepBuild_Builder.cxx L2003-2055 — FillVertexSetOnValue.
    // =========================================================================
    fn fill_vertex_set_on_value(
        &self,
        it: &TopOpeBRepDSInterference,
        ispoint: bool,
        tobuild: TopAbsState,
        pvs: &mut TopOpeBRepBuildPaveSet,
    ) {
        let TopOpeBRepDSInterference::CurvePoint(cpi) = it else {
            return;
        };
        // ind = index of new point or existing vertex.
        let ind = cpi.index_g;
        let v = if ispoint && ind <= self.my_data_structure.as_ref().map(|d| d.points.len() as i32).unwrap_or(0) {
            // OCCT: V = NewVertex(ind) — the vertex built by BuildVertices
            // (pending unit); the DS shape table carries the same entry.
            self.my_data_structure
                .as_ref()
                .expect("DS")
                .shape(ind)
                .clone()
        } else {
            self.my_data_structure
                .as_ref()
                .expect("DS")
                .shape(ind)
                .clone()
        };
        let par = cpi.parameter;
        // OCCT: IT.Orientation(ToBuild) — the TopOpeBRepDS_PointIterator
        // orientation for the state (pending PointIterator unit; the
        // interference transition orientation is the payload).
        let ori = it.transition().orientation_in();
        let _ = tobuild;

        let keep = true;
        if keep {
            let mut v = v;
            v.orientation = ori;
            let pv = TopOpeBRepBuildPave::new(v, par, false);
            pvs.append(pv);
        }
    }
}

// OCCT PaveSet.cxx L131-141 — FUN_islook.
fn fun_islook(brep: &rcad_kernel::topods::BRep, e: &Shape) -> bool {
    let ed = e.as_edge().expect("not an edge");
    let p1 = brep.vertex_position(&ed.first);
    let p2 = brep.vertex_position(&ed.last);
    let dp1p2 = p1.distance(p2);
    dp1p2.abs() > 1.0e-8
}

impl TopOpeBRepBuildPaveSet {
    /// OCCT PaveSet.cxx L64-129 — SortPave(List, SortedList): the n^2
    /// selection sort on Parameter(), then the FORWARD head move
    /// (tete = FORWARD).
    pub fn sort_pave(list: &[TopOpeBRepBuildPave]) -> Vec<TopOpeBRepBuildPave> {
        let n_pv = list.len();
        let mut taken = vec![false; n_pv];
        let mut sorted: Vec<TopOpeBRepBuildPave> = Vec::new();

        for _ in 0..n_pv {
            let mut parmin = f64::MAX;
            let mut chosen: Option<usize> = None;
            for (itest, pv) in list.iter().enumerate() {
                if !taken[itest] {
                    let par = pv.parameter();
                    if par < parmin {
                        parmin = par;
                        chosen = Some(itest);
                    }
                }
            }
            if let Some(i) = chosen {
                sorted.push(list[i].clone());
                taken[i] = true;
            }
        }

        // tete = FORWARD.
        let mut found = false;
        let mut l1: Vec<TopOpeBRepBuildPave> = Vec::new();
        let mut l2: Vec<TopOpeBRepBuildPave> = Vec::new();
        for pv in &sorted {
            if !found {
                if pv.vertex.orientation == Orientation::Forward {
                    found = true;
                    l1.push(pv.clone());
                } else {
                    l2.push(pv.clone());
                }
            } else {
                l1.push(pv.clone());
            }
        }
        l1.extend(l2);
        l1
    }

    /// OCCT PaveSet.cxx L142-297 — Prepare(): add the edge vertices to the
    /// list of paves; an edge vertex already present as an interference
    /// vertex is merged (INTERNAL keeps/adopts the edge orientation,
    /// EXTERNAL or opposite-orientation entries are removed); then the
    /// list is sorted on Parameter().
    pub fn prepare(&mut self, brep: &rcad_kernel::topods::BRep) {
        if self.prepare_done {
            return;
        }

        let is_ed = brep.is_edge_degenerated(&self.edge);
        let mut edge_vertex_count = 0usize;

        // myRemovePV is true by default (jyl + 980217).
        {
            // OCCT TopExp_Explorer(myEdge, VERTEX): first vertex FORWARD,
            // last REVERSED.
            let ed = self.edge.as_edge().expect("not an edge");
            let explorer = [
                (ed.first.clone(), Orientation::Forward),
                (ed.last.clone(), Orientation::Reversed),
            ];
            for (ve, veori) in explorer {
                let vebound = veori == Orientation::Forward || veori == Orientation::Reversed;

                let mut edge_vertex_index = 0usize;
                let mut add_ve = true;
                let mut add = false; // ofv
                let mut remove_idx: Option<usize> = None;
                let mut set_ori_idx: Option<(usize, Orientation)> = None;

                for (idx, pv) in self.vertices.iter().enumerate() {
                    edge_vertex_index += 1; // skip edge vertices inserted at the head
                    if edge_vertex_index <= edge_vertex_count {
                        continue;
                    }

                    // PV = Parametrized vertex, VI = interference vertex.
                    let vi = pv.vertex();
                    let has_vsd = pv.has_same_domain;
                    let vi_ori = vi.orientation;
                    let visameve = vi.is_same(&ve);
                    let mut vsdsameve = false;
                    if has_vsd {
                        if let Some(vsd) = &pv.same_domain {
                            vsdsameve = vsd.is_same(&ve);
                        }
                    }
                    let samevertexprocessing = (visameve || vsdsameve) && !is_ed;

                    if samevertexprocessing && (vebound || vsdsameve) {
                        match vi_ori {
                            Orientation::External => {
                                remove_idx = Some(idx);
                            }
                            Orientation::Internal => {
                                // OCCT: VI orientation adopts the edge
                                // vertex orientation.
                                set_ori_idx = Some((idx, veori));
                            }
                            _ => {
                                if vi_ori != veori {
                                    remove_idx = Some(idx);
                                    let islook = fun_islook(brep, &self.edge);
                                    if (vebound && (vsdsameve || visameve)) && islook {
                                        add = true; // ofv
                                    }
                                }
                            }
                        }
                    }
                    // ofv: addVE = add; break.
                    add_ve = add;
                    break;
                }

                if let Some(idx) = remove_idx {
                    self.vertices.remove(idx);
                }
                if let Some((idx, o)) = set_ori_idx {
                    self.vertices[idx].change_vertex().orientation = o;
                }
                // if VE not found in the list, add it.
                if add_ve {
                    let par_ve = {
                        let b = brep;
                        super::chfi3d_builder_0::brep_tool_parameter(b, &ve, &self.edge)
                    };
                    let new_pv = TopOpeBRepBuildPave::new(ve, par_ve, true);
                    self.vertices.insert(0, new_pv);
                    edge_vertex_count += 1;
                }
            }
        } // myRemovePV

        let ll = self.vertices.len();

        // if no more interference vertices, clear the list.
        if ll == edge_vertex_count {
            self.vertices.clear();
        } else if ll >= 2 {
            // sort the parametrized vertices on Parameter() value.
            let list = self.vertices.clone();
            self.vertices = Self::sort_pave(&list);
        }

        self.prepare_done = true;
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_PaveClassifier (PaveClassifier.hxx + .cxx L35-392).
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepBuildPaveClassifier {
    my_edge: Shape,
    my_edge_periodic: bool,
    my_same_parameters: bool,
    my_closed_vertices: bool,
    my_first: f64,
    my_period: f64,
    my_o1: Orientation,
    my_o2: Orientation,
    my_p1: f64,
    my_p2: f64,
    my_cas1: i32,
    my_cas2: i32,
}

impl TopOpeBRepBuildPaveClassifier {
    /// OCCT PaveClassifier.cxx L35-98 — PaveClassifier(E).
    pub fn new(brep: &rcad_kernel::topods::BRep, e: &Shape) -> Self {
        let mut pc = TopOpeBRepBuildPaveClassifier {
            my_edge: e.clone(),
            my_edge_periodic: false,
            my_same_parameters: false,
            my_closed_vertices: false,
            my_first: 0.0,
            my_period: 0.0,
            my_o1: Orientation::Forward,
            my_o2: Orientation::Forward,
            my_p1: 0.0,
            my_p2: 0.0,
            my_cas1: 0,
            my_cas2: 0,
        };

        if !brep.is_edge_degenerated(&pc.my_edge) {
            let r: Option<(rcad_kernel::geom::Curve3, [f64; 2])> = brep.edge_curve_world(&pc.my_edge);
            let Some((c, fl)) = r else {
                return pc;
            };
            let (f, l) = (fl[0], fl[1]);
            if c.is_periodic() {
                let ed = pc.my_edge.as_edge().expect("not an edge");
                let v1 = ed.first.clone();
                let v2 = ed.last.clone(); // v1 FORWARD, v2 REVERSED
                if !v1.is_null() && !v2.is_null() {
                    // --- the edge has vertices.
                    pc.my_first = f;
                    let domain = c.default_domain();
                    let f_c = domain[0];
                    let l_c = domain[1];
                    pc.my_period = l_c - f_c;
                    pc.my_edge_periodic = v1.is_same(&v2);
                    pc.my_same_parameters = v1.is_same(&v2);
                    if pc.my_same_parameters {
                        pc.my_first = brep_tool_parameter(brep, &v1, &pc.my_edge);
                    }
                } else {
                    // --- the edge has no vertices.
                    pc.my_first = f;
                    pc.my_period = l - f;
                    pc.my_edge_periodic = true;
                    pc.my_same_parameters = false;
                }
            }
        } // ! degenerated
        pc
    }

    /// OCCT PaveClassifier.cxx L102-177 — CompareOnNonPeriodic.
    pub fn compare_on_non_periodic(&self) -> TopAbsState {
        let mut state = TopAbsState::Unknown;
        let lower;
        match self.my_o2 {
            Orientation::Forward => lower = false,
            Orientation::Reversed => lower = true,
            Orientation::Internal => {
                state = TopAbsState::In;
                lower = false;
            }
            Orientation::External => {
                state = TopAbsState::Out;
                lower = false;
            }
        }

        if state == TopAbsState::Unknown {
            if self.my_p1 == self.my_p2 {
                if self.my_o1 == self.my_o2 {
                    state = TopAbsState::In;
                } else {
                    state = TopAbsState::Out;
                }
            } else if self.my_p1 < self.my_p2 {
                state = if lower {
                    TopAbsState::In
                } else {
                    TopAbsState::Out
                };
            } else {
                state = if lower {
                    TopAbsState::Out
                } else {
                    TopAbsState::In
                };
            }
        }

        state
    }

    /// OCCT PaveClassifier.cxx L181-217 — AdjustCase.
    #[allow(clippy::too_many_arguments)]
    pub fn adjust_case(p1: f64, o: Orientation, first: f64, period: f64, tol: f64, cas: &mut i32) -> f64 {
        let p2;
        if (p1 - first).abs() < tol {
            // p1 is first.
            if o == Orientation::Reversed {
                p2 = p1 + period;
                *cas = 1;
            } else {
                p2 = p1;
                *cas = 2;
            }
        } else {
            // p1 is not on first.
            let last = first + period;
            if (p1 - last).abs() < tol {
                // p1 is on last.
                p2 = p1;
                *cas = 3;
            } else {
                // p1 is not on last.
                p2 = super::chfi_ds::elclib_in_period(p1, first, last);
                *cas = 4;
            }
        }
        p2
    }

    /// OCCT PaveClassifier.cxx L266-270 — ToAdjustOnPeriodic.
    pub fn to_adjust_on_periodic(&self) -> bool {
        self.my_same_parameters || (self.my_o1 != self.my_o2)
    }

    /// OCCT PaveClassifier.cxx L221-262 — AdjustOnPeriodic.
    pub fn adjust_on_periodic(&mut self) {
        if !self.to_adjust_on_periodic() {
            return;
        }

        let tol = 1.0e-9; // OCCT Precision::PConfusion()

        if self.my_same_parameters {
            let mut p1 = self.my_p1;
            let mut p2 = self.my_p2;
            p1 = Self::adjust_case(p1, self.my_o1, self.my_first, self.my_period, tol, &mut self.my_cas1);
            p2 = Self::adjust_case(p2, self.my_o2, self.my_first, self.my_period, tol, &mut self.my_cas2);
            self.my_p1 = p1;
            self.my_p2 = p2;
        } else if self.my_o1 != self.my_o2 {
            if self.my_o1 == Orientation::Forward {
                self.my_p2 = Self::adjust_case(
                    self.my_p2,
                    self.my_o2,
                    self.my_p1,
                    self.my_period,
                    tol,
                    &mut self.my_cas2,
                );
            }
            if self.my_o2 == Orientation::Forward {
                self.my_p1 = Self::adjust_case(
                    self.my_p1,
                    self.my_o1,
                    self.my_p2,
                    self.my_period,
                    tol,
                    &mut self.my_cas1,
                );
            }
        }
    }

    /// OCCT PaveClassifier.cxx L274-309 — CompareOnPeriodic.
    pub fn compare_on_periodic(&mut self) -> TopAbsState {
        let state;
        if self.to_adjust_on_periodic() {
            state = self.compare_on_non_periodic();
        } else if self.my_o1 == Orientation::Forward {
            state = TopAbsState::Out;
            self.my_cas1 = 5;
            self.my_cas2 = 5;
        } else if self.my_o1 == Orientation::Reversed {
            state = TopAbsState::Out;
            self.my_cas1 = 6;
            self.my_cas2 = 6;
        } else {
            state = TopAbsState::Out;
            self.my_cas1 = 7;
            self.my_cas2 = 7;
        }
        state
    }

    /// OCCT PaveClassifier.cxx L313-362 — Compare(L1, L2).
    pub fn compare(&mut self, pv1: &TopOpeBRepBuildPave, pv2: &TopOpeBRepBuildPave) -> TopAbsState {
        self.my_cas1 = 0; // debug
        self.my_cas2 = 0; // debug
        self.my_o1 = pv1.vertex().orientation;
        self.my_o2 = pv2.vertex().orientation;
        self.my_p1 = pv1.parameter();
        self.my_p2 = pv2.parameter();

        if self.my_edge_periodic && self.to_adjust_on_periodic() {
            self.adjust_on_periodic();
        }

        if self.my_edge_periodic {
            self.compare_on_periodic()
        } else {
            self.compare_on_non_periodic()
        }
    }

    /// OCCT PaveClassifier.cxx L366-375 — SetFirstParameter(P).
    pub fn set_first_parameter(&mut self, p: f64) {
        self.my_first = p;
        self.my_same_parameters = true;
    }

    /// OCCT PaveClassifier.cxx L379-392 — ClosedVertices(Closed).
    pub fn set_closed_vertices(&mut self, closed: bool) {
        self.my_closed_vertices = closed;
        if closed {
            self.my_edge_periodic = true;
        }
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Loop / LoopSet / LoopClassifier abstract interfaces
// (TopOpeBRepBuild_Loop.hxx / _LoopSet.hxx / _LoopClassifier.hxx) —
// rcad: traits; TopOpeBRepBuild_Pave implements Loop (IsShape() = the
// bound flag: an old edge vertex is a boundary loop, a new interference
// vertex is a block loop), PaveSet implements LoopSet, PaveClassifier
// implements LoopClassifier.
// =========================================================================
pub trait AreaLoop {
    /// OCCT TopOpeBRepBuild_Loop::IsShape().
    fn is_shape(&self) -> bool;
}

pub trait AreaLoopSet {
    /// OCCT TopOpeBRepBuild_LoopSet::InitLoop/MoreLoop/NextLoop/Loop.
    fn loops(&self) -> &[TopOpeBRepBuildPave];
}

pub trait AreaLoopClassifier {
    /// OCCT TopOpeBRepBuild_LoopClassifier::Compare(L1, L2) -> TopAbs_State.
    fn compare(&mut self, l1: &TopOpeBRepBuildPave, l2: &TopOpeBRepBuildPave) -> TopAbsState;
}

impl AreaLoop for TopOpeBRepBuildPave {
    /// OCCT Pave.hxx IsShape() — the bound flag (old vertex = shape).
    fn is_shape(&self) -> bool {
        self.is_bound()
    }
}

impl AreaLoopSet for TopOpeBRepBuildPaveSet {
    fn loops(&self) -> &[TopOpeBRepBuildPave] {
        &self.vertices
    }
}

impl AreaLoopClassifier for TopOpeBRepBuildPaveClassifier {
    fn compare(&mut self, l1: &TopOpeBRepBuildPave, l2: &TopOpeBRepBuildPave) -> TopAbsState {
        self.compare(l1, l2)
    }
}

/// OCCT TopOpeBRepBuild_LoopEnum (TopOpeBRepBuild_AreaBuilder.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoopEnum {
    AnyLoop,
    Boundary,
    Block,
}

// =========================================================================
// OCCT TopOpeBRepBuild_AreaBuilder (AreaBuilder.hxx + .cxx L30-429).
// =========================================================================
#[derive(Debug, Clone)]
pub struct AreaBuilder {
    /// OCCT: bool myUNKNOWNRaise — false by default.
    my_unknown_raise: bool,
    /// OCCT: TopOpeBRepBuild_ListOfListOfLoop myArea.
    my_area: Vec<Vec<usize>>,
    area_index: usize,
    loop_index: usize,
}

impl Default for AreaBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl AreaBuilder {
    /// OCCT AreaBuilder.cxx L30-42.
    pub fn new() -> Self {
        AreaBuilder {
            my_unknown_raise: false, // no raise if UNKNOWN state found
            my_area: Vec::new(),
            area_index: 0,
            loop_index: 0,
        }
    }

    /// OCCT AreaBuilder.cxx L56-103 — CompareLoopWithListOfLoop: compare
    /// the position of Loop <L> with the Area <LOL> using the classifier;
    /// according to <what>, loops of <LOL> are selected or not.
    ///   TopAbs_OUT if <LOL> is empty; UNKNOWN if undefined;
    ///   IN if <L> is inside all the selected Loops;
    ///   OUT if <L> is outside one of the selected Loops.
    fn compare_loop_with_list_of_loop<C: AreaLoopClassifier>(
        &self,
        lc: &mut C,
        l: usize,
        lol: &[usize],
        loops: &[TopOpeBRepBuildPave],
        what: LoopEnum,
    ) -> TopAbsState {
        let mut state = TopAbsState::Unknown;
        if lol.is_empty() {
            return TopAbsState::Out;
        }

        for &cur_l in lol {
            let totest = match what {
                LoopEnum::AnyLoop => true,
                LoopEnum::Boundary => loops[cur_l].is_shape(),
                LoopEnum::Block => !loops[cur_l].is_shape(),
            };
            if totest {
                state = lc.compare(&loops[l], &loops[cur_l]);
                if state == TopAbsState::Out {
                    // <L> is out of at least one Loop of <LOL>: stop.
                    break;
                }
            }
        }
        state
    }

    /// OCCT AreaBuilder.cxx L111-122 — Atomize.
    fn atomize(&self, state: &mut TopAbsState, newstate: TopAbsState) {
        if self.my_unknown_raise {
            if *state == TopAbsState::Unknown {
                panic!("Standard_DomainError: AreaBuilder : Position Unknown");
            }
        } else {
            *state = newstate;
        }
    }

    /// OCCT AreaBuilder.cxx L125-330 — InitAreaBuilder(LS, LC, ForceClass):
    /// the area construction over the loop set (block loops open areas,
    /// boundary loops are absorbed into the areas containing them).
    /// `loops` is the flat loop list; areas reference loops by index.
    pub fn init_area_builder<C: AreaLoopClassifier>(
        &mut self,
        loops: &[TopOpeBRepBuildPave],
        lc: &mut C,
        force_class: bool,
    ) {
        // boundaryloops : list of boundary loops out of the areas.
        let mut my_area: Vec<Vec<usize>> = Vec::new();
        let mut boundaryloops: Vec<usize> = Vec::new();

        for l in 0..loops.len() {
            // process a new loop : L is the new current Loop.
            let boundary_l = loops[l].is_shape();

            // L = Shape et ForceClass : on traite L comme un block
            // L = Shape et !ForceClass : on traite L comme un pur Shape
            // L = !Shape : on traite L comme un block
            let traitercommeblock = (!boundary_l) || force_class;
            if !traitercommeblock {
                // the loop L is a boundary loop: try to insert it in an
                // existing area such as L is inside all the block loops.
                let mut loopinside = false;
                let mut area_hit = 0usize;
                for (ai, a_area) in my_area.iter().enumerate() {
                    if a_area.is_empty() {
                        continue;
                    }
                    let mut state =
                        self.compare_loop_with_list_of_loop(lc, l, a_area, loops, LoopEnum::Block);
                    if state == TopAbsState::Unknown {
                        self.atomize(&mut state, TopAbsState::In);
                    }
                    loopinside = state == TopAbsState::In;
                    if loopinside {
                        area_hit = ai;
                        break;
                    }
                } // end of Area scan

                if loopinside {
                    // OCCT: ADD_Loop_TO_LISTOFLoop(L, aArea) — aArea is the
                    // area the scan stopped on.
                    my_area[area_hit].push(l);
                } else {
                    boundaryloops.push(l);
                }
            } else {
                // the loop L is a block loop.
                let mut loopinside = false;
                let mut area_hit = 0usize;
                for (ai, a_area) in my_area.iter().enumerate() {
                    if a_area.is_empty() {
                        continue;
                    }
                    let mut state =
                        self.compare_loop_with_list_of_loop(lc, l, a_area, loops, LoopEnum::AnyLoop);
                    if state == TopAbsState::Unknown {
                        self.atomize(&mut state, TopAbsState::In);
                    }
                    loopinside = state == TopAbsState::In;
                    if loopinside {
                        area_hit = ai;
                        break;
                    }
                } // end of Area scan

                if loopinside {
                    let a_area = &mut my_area[area_hit];
                    let mut all_shape = true;
                    let mut removed_loops: Vec<usize> = Vec::new();
                    let mut li = 0usize;
                    while li < a_area.len() {
                        let mut state = lc.compare(&loops[a_area[li]], &loops[l]);
                        if state == TopAbsState::Unknown {
                            self.atomize(&mut state, TopAbsState::In); // not OUT
                        }
                        let loopoutside = state == TopAbsState::Out;
                        if loopoutside {
                            let cur_l = a_area[li];
                            removed_loops.push(cur_l);
                            all_shape = all_shape && loops[cur_l].is_shape();
                            a_area.remove(li);
                        } else {
                            li += 1;
                        }
                    }
                    // insert the loop in the area.
                    a_area.push(l);
                    if !removed_loops.is_empty() {
                        if all_shape {
                            boundaryloops.extend(removed_loops);
                        } else {
                            // make a new area with the removed loops.
                            my_area.push(removed_loops);
                        }
                    }
                } else {
                    // create a new area with L; insert boundary loops that
                    // are IN the new area (and remove them from
                    // 'boundaryloops').
                    let mut new_area0: Vec<usize> = vec![l];
                    let mut bi = 0usize;
                    while bi < boundaryloops.len() {
                        let cur_l = boundaryloops[bi];
                        let mut state = lc.compare(&loops[cur_l], &loops[l]);
                        if state == TopAbsState::Unknown {
                            self.atomize(&mut state, TopAbsState::In);
                        }
                        let ashapeinside = state == TopAbsState::In && loops[cur_l].is_shape();
                        let mut ablockinside = false;
                        if ashapeinside {
                            let mut state2 = lc.compare(&loops[l], &loops[cur_l]);
                            if state2 == TopAbsState::Unknown {
                                self.atomize(&mut state2, TopAbsState::In);
                            }
                            ablockinside = state2 == TopAbsState::In;
                        }
                        if ashapeinside && ablockinside {
                            new_area0.push(cur_l);
                            boundaryloops.remove(bi);
                        } else {
                            bi += 1;
                        }
                    } // end of boundaryloops scan
                    my_area.push(new_area0);
                } // Loopinside == False
            } // end of block loop
        } // end of LoopSet LS scan

        self.my_area = my_area;
        self.init_area();
    }

    /// OCCT AreaBuilder.cxx L334-341 — InitArea.
    pub fn init_area(&mut self) -> usize {
        self.area_index = 0;
        self.init_loop();
        self.my_area.len()
    }

    /// OCCT AreaBuilder.cxx L344-349 — MoreArea.
    pub fn more_area(&self) -> bool {
        self.area_index < self.my_area.len()
    }

    /// OCCT AreaBuilder.cxx L352-357 — NextArea.
    pub fn next_area(&mut self) {
        self.area_index += 1;
        self.init_loop();
    }

    /// OCCT AreaBuilder.cxx L360-374 — InitLoop.
    pub fn init_loop(&mut self) -> usize {
        if self.area_index < self.my_area.len() {
            self.my_area[self.area_index].len()
        } else {
            self.loop_index = 0;
            0
        }
    }

    /// OCCT AreaBuilder.cxx L378-383 — MoreLoop.
    pub fn more_loop(&self) -> bool {
        self.loop_index < self.my_area[self.area_index.min(self.my_area.len() - 1)].len()
    }

    /// OCCT AreaBuilder.cxx L386-390 — NextLoop.
    pub fn next_loop(&mut self) {
        self.loop_index += 1;
    }

    /// OCCT AreaBuilder.cxx L393-398 — Loop().
    pub fn loop_pave_index(&self) -> usize {
        self.my_area[self.area_index][self.loop_index]
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_EdgeBuilder (EdgeBuilder.hxx + EdgeBuilder.cxx
// L25-105) — a thin delegate over AreaBuilder for the edge case: the areas
// are built over the PaveSet (loops) with the PaveClassifier.
// =========================================================================
#[derive(Debug, Clone)]
pub struct TopOpeBRepBuildEdgeBuilder {
    pub area: AreaBuilder,
    /// OCCT: the PaveSet (LoopSet) and Pave loops.
    pub loops: Vec<TopOpeBRepBuildPave>,
    pub pave_set_edge: Shape,
}

impl Default for TopOpeBRepBuildEdgeBuilder {
    fn default() -> Self {
        TopOpeBRepBuildEdgeBuilder {
            area: AreaBuilder::new(),
            loops: Vec::new(),
            pave_set_edge: Shape::null(),
        }
    }
}

impl TopOpeBRepBuildEdgeBuilder {
    /// OCCT EdgeBuilder.cxx L34-45 — ctor + InitEdgeBuilder(LS, LC,
    /// ForceClass) -> InitAreaBuilder.
    pub fn init_edge_builder(
        &mut self,
        pvs: &TopOpeBRepBuildPaveSet,
        lc: &mut TopOpeBRepBuildPaveClassifier,
        force_class: bool,
        brep: &rcad_kernel::topods::BRep,
    ) {
        self.loops = pvs.vertices.clone();
        self.pave_set_edge = pvs.edge.clone();
        let _ = brep;
        self.area.init_area_builder(&self.loops, lc, force_class);
    }

    /// OCCT EdgeBuilder.cxx L47-50 — InitEdge() -> InitArea().
    pub fn init_edge(&mut self) {
        self.area.init_area();
    }

    /// OCCT EdgeBuilder.cxx L52-55 — MoreEdge() -> MoreArea().
    pub fn more_edge(&self) -> bool {
        self.area.more_area()
    }

    /// OCCT EdgeBuilder.cxx L57-60 — NextEdge() -> NextArea().
    pub fn next_edge(&mut self) {
        self.area.next_area();
    }

    /// OCCT EdgeBuilder.cxx L62-65 — InitVertex() -> InitLoop().
    pub fn init_vertex(&mut self) {
        self.area.init_loop();
    }

    /// OCCT EdgeBuilder.cxx L67-70 — MoreVertex() -> MoreLoop().
    pub fn more_vertex(&self) -> bool {
        self.area.more_loop()
    }

    /// OCCT EdgeBuilder.cxx L72-75 — NextVertex() -> NextLoop().
    pub fn next_vertex(&mut self) {
        self.area.next_loop();
    }

    /// OCCT EdgeBuilder.cxx L77-84 — Vertex(): the pave vertex of the
    /// current loop.
    pub fn vertex(&self) -> &Shape {
        let idx = self.area.loop_pave_index();
        self.loops[idx].vertex()
    }

    /// OCCT EdgeBuilder.cxx L86-93 — Parameter().
    pub fn parameter(&self) -> f64 {
        let idx = self.area.loop_pave_index();
        self.loops[idx].parameter()
    }
}
