//! OCCT TopOpeBRepBuild WireEdgeSet / FaceBuilder / SplitFace1 machinery —
//! the SplitFace1 support-face reconstruction batch of the
//! TopOpeBRepBuild_HBuilder D6 carrier (fillet::hbuilder).
//!
//! D6 ruling (2026-09-07): TKBool code (TopOpeBRepBuild / TopOpeBRepDS /
//! TopOpeBRepTool) is NOT translated as a standalone package; wherever
//! TKFillet / TKOffset / TKFeat depend on it, the translation lands as a
//! component of the consuming carrier.  These classes are physically in
//! TKBool/TopOpeBRepBuild but per the D6 ruling they translate here, as
//! the SplitFace1 arm of the HBuilder carrier (the build_faces precedent).
//!
//! Translated here (all paths under
//! OCCT src/ModelingAlgorithms/TKBool/TopOpeBRepBuild/):
//! - `TopOpeBRepBuild_WireEdgeSet` (WireEdgeSet.cxx L51-575 + .hxx)
//! - `TopOpeBRepBuild_FaceBuilder` (FaceBuilder.cxx L59-661 + .hxx)
//! - `TopOpeBRepBuild_Builder::SplitFace / SplitFace1` (Builder.cxx
//!   L1155-1300), `FillFace` (L1951-1973), `FillShape` (L1876-1947),
//!   `SplitShapes` (L1713-1872) + `FUN_touched` (L1702-1709),
//!   `AddIntersectionEdges` (L150-178), `ToSplit` (L308-326),
//!   `Reverse` (L621-633), `Orient` (L637-640), `KeepShape` (L540-554),
//!   `ShapePosition` (L499-531), `FindSameDomain` (L650-709),
//!   `Contains` (L836-847), `MarkSplit` (L330-374), `ChangeSplit`
//!   (L449-491)
//! - `TopOpeBRepBuild_Builder::MakeFaces` (Merge.cxx L435-512) +
//!   `CorrectEdgeOrientation` (Merge.cxx L48-136) +
//!   `CorrectUnclosedWire` (Merge.cxx L138-170), `ChangeMerged`
//!   (Merge.cxx L700-728)
//!
//! # Wiring point (for the main agent — NOT wired in this batch)
//!
//! The OCCT call site of the whole arm is the solid reconstruction walk:
//! MergeSolid -> MergeShapes (Merge.cxx L174+) -> FillSolid
//! (L1979-1987) -> FillShape (SOLID) explores shells, and for every shell
//! present in the DS, SplitShapes walks its faces and calls
//! SplitFace(face, ToBuild1, ToBuild2) (Builder.cxx L1737), which lands in
//! SplitFace1.  The expected rcad wiring (the main agent owns it) is
//! inside `TopOpeBRepBuildHBuilder::merge_solid` (fillet/hbuilder.rs),
//! replacing the pending-boundary block that currently feeds the as-found
//! faces to the BOPAlgo_BuilderSolid: for every face of the object solid
//! that passes the ToSplit gate, call
//!
//! ```text
//! coup.split_face1(&mut my_brep, ds, &face, TopAbsState::In, TopAbsState::In)
//! ```
//!
//! and then feed `splits(face, ToBuild1)` (the reconstructed faces
//! recorded by split_face1 into ChangeSplit) to the solid builder where
//! merge_solid today pushes the as-found face.  The FaceBuilder /
//! FillFace orientation-compound semantics (piece orientation x
//! wire-edge orientation, FillFace RevOri) are the core of the mechanism —
//! see `split_face1` below; every classifier judgement follows the OCCT
//! statement, no edge-count-style approximations anywhere.

// The carrier is unwired by design (the merge_solid wiring point belongs
// to the main agent — see the section above), so most of its surface is
// dead code until then.
#![allow(dead_code)]

mod classify;

use std::sync::Arc;

use rcad_kernel::core::precision::{is_infinite_value, CONFUSION};
use rcad_kernel::geom::{Curve2dEval as _, CurveEval as _};
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool as _, Orientation, Shape, ShapeType, TShape};

use super::chfi3d_builder_2::TopAbsState;
use super::chfi3d_ds::{TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference, TopOpeBRepDSKind};
use super::hbuilder::{ShapeStateMap, TopOpeBRepBuildHBuilder};

use classify::{
    bb_add_face_wire, shape_explore_children, shape_key, shape_oriented, shapes_same,
    top_abs_complement, top_exp_vertices, BlockBuilder, BlockIterator, FaceAreaBuilder, Loop,
    LoopSet, ShapeSet, WireEdgeClassifier,
};

// =========================================================================
// OCCT TopOpeBRepBuild_Builder.cxx L73 — static int STATIC_SOLIDINDEX = 0;
// (set to 1 / 2 by SplitSolid, Builder.cxx L1609 / L1617; the SplitFace1
// arm only reads it in the SplitShapes testkeep branch L1828).
// =========================================================================
static STATIC_SOLIDINDEX: std::sync::atomic::AtomicI32 = std::sync::atomic::AtomicI32::new(0);

// The WireEdgeSet (WireEdgeSet.cxx L51-575) continuation file.
mod wes;

pub(crate) use wes::WireEdgeSet;

// The SplitEdge1 chain (PaveSet / PaveClassifier / Area1dBuilder /
// EdgeBuilder / SplitEdge / MakeEdges) continuation file.
mod split_edge;

// =========================================================================
// OCCT TopOpeBRepBuild_FaceBuilder (FaceBuilder.cxx L59-661 + .hxx).
// =========================================================================

// FaceBuilder.cxx L118-122.
const ISUNKNOWN: i32 = -1;
const ISVERTEX: i32 = 0;
const GCLOSEDW: i32 = 1;
const UNCLOSEDW: i32 = 2;
const CLOSEDW: i32 = 10;

pub(crate) struct FaceBuilder {
    pub(crate) my_face: Shape,
    pub(crate) my_loop_set: LoopSet,
    pub(crate) my_block_iterator: BlockIterator,
    pub(crate) my_block_builder: BlockBuilder,
    pub(crate) my_face_area_builder: FaceAreaBuilder,
}

impl Default for FaceBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl FaceBuilder {
    /// FaceBuilder.cxx L59.
    pub(crate) fn new() -> Self {
        FaceBuilder {
            my_face: Shape::null(),
            my_loop_set: LoopSet::new(),
            my_block_iterator: BlockIterator::new(),
            my_block_builder: BlockBuilder::new(),
            my_face_area_builder: FaceAreaBuilder::new(),
        }
    }

    /// FaceBuilder.cxx L63-68.
    #[allow(dead_code)]
    pub(crate) fn new_init(
        brep: &mut BRep,
        wes: &mut WireEdgeSet,
        f: &Shape,
        force_class: bool,
    ) -> Self {
        let mut fb = FaceBuilder::new();
        fb.init_face_builder(brep, wes, f, force_class);
        fb
    }

    /// FaceBuilder.cxx L72-82.
    pub(crate) fn init_face_builder(
        &mut self,
        brep: &mut BRep,
        wes: &mut WireEdgeSet,
        f: &Shape,
        force_class: bool,
    ) {
        self.my_face = f.clone(); // TopoDS::Face(F)
        self.make_loops(brep, wes);
        // TopOpeBRepBuild_BlockBuilder& BB = myBlockBuilder;
        // TopOpeBRepBuild_WireEdgeClassifier WEC(F, BB);
        // TopOpeBRepBuild_LoopSet& LS = myLoopSet;
        // myFaceAreaBuilder.InitFaceAreaBuilder(LS, WEC, ForceClass).
        let wes_shapes: Vec<Shape> = wes.base.my_shapes.clone();
        let bb = &self.my_block_builder;
        let mut wec = WireEdgeClassifier::new(brep, f, bb, &wes_shapes);
        let ls = &mut self.my_loop_set;
        self.my_face_area_builder
            .init_face_area_builder(brep, ls, &mut wec, force_class);
    }

    /// FaceBuilder.cxx L451-454.
    pub(crate) fn face(&self) -> &Shape {
        &self.my_face
    }

    /// FaceBuilder.cxx L458-462.
    pub(crate) fn init_face(&mut self) -> i32 {
        let n = self.my_face_area_builder.init_area();
        n
    }

    /// FaceBuilder.cxx L466-470.
    pub(crate) fn more_face(&self) -> bool {
        let b = self.my_face_area_builder.more_area();
        b
    }

    /// FaceBuilder.cxx L474-477.
    pub(crate) fn next_face(&mut self) {
        self.my_face_area_builder.next_area();
    }

    /// FaceBuilder.cxx L481-485.
    pub(crate) fn init_wire(&mut self) -> i32 {
        let n = self.my_face_area_builder.init_loop();
        n
    }

    /// FaceBuilder.cxx L489-493.
    pub(crate) fn more_wire(&self) -> bool {
        let b = self.my_face_area_builder.more_loop();
        b
    }

    /// FaceBuilder.cxx L497-500.
    pub(crate) fn next_wire(&mut self) {
        self.my_face_area_builder.next_loop();
    }

    /// FaceBuilder.cxx L504-509.
    pub(crate) fn is_old_wire(&self) -> bool {
        let l = self.my_face_area_builder.loop_();
        let b = l.is_shape();
        b
    }

    /// FaceBuilder.cxx L513-518.
    pub(crate) fn old_wire(&self) -> Shape {
        let l = self.my_face_area_builder.loop_();
        let b = l.shape().clone();
        b
    }

    /// FaceBuilder.cxx L522-541.
    pub(crate) fn find_next_valid_element(&mut self) {
        // prerequisites : myBlockIterator.Initialize
        let _ = self.my_face_area_builder.loop_();
        let mut found = false;

        while self.my_block_iterator.more() {
            let i = self.my_block_iterator.value();
            found = self.my_block_builder.element_is_valid_index(i);
            if found {
                break;
            } else {
                self.my_block_iterator.next();
            }
        }
    }

    /// FaceBuilder.cxx L545-560.
    pub(crate) fn init_edge(&mut self) -> i32 {
        let l = self.my_face_area_builder.loop_();
        if l.is_shape() {
            panic!("TopOpeBRepBuild_FaceBuilder:InitEdge");
        } else {
            self.my_block_iterator = *l.block_iterator();
            self.my_block_iterator.initialize();
            self.find_next_valid_element();
        }
        let n = self.my_block_iterator.extent();
        n
    }

    /// FaceBuilder.cxx L564-568.
    pub(crate) fn more_edge(&self) -> bool {
        let b = self.my_block_iterator.more();
        b
    }

    /// FaceBuilder.cxx L572-576.
    pub(crate) fn next_edge(&mut self) {
        self.my_block_iterator.next();
        self.find_next_valid_element();
    }

    /// FaceBuilder.cxx L580-596.
    pub(crate) fn edge(&self) -> Shape {
        if !self.my_block_iterator.more() {
            panic!("OutOfRange");
        }

        let i = self.my_block_iterator.value();
        let isvalid = self.my_block_builder.element_is_valid_index(i);
        if !isvalid {
            panic!("Edge not Valid");
        }

        let e = self.my_block_builder.element_index(i);
        e
    }

    /// FaceBuilder.cxx L600-610 — EdgeConnexity management disactivated
    /// (returns 0; the myMOSI debug map has no live consumer in OCCT 8).
    pub(crate) fn edge_connexity(&self, _e: &Shape) -> i32 {
        0
    }

    /// FaceBuilder.cxx L614-629.
    pub(crate) fn add_edge_wire(&self, brep: &mut BRep, e: &Shape, w: &mut Shape) -> i32 {
        let mut nadd = 0;
        let mut bb = BRepBuilder::new();
        bb.add_to_wire(brep, w.clone(), e.clone());
        nadd += 1;
        nadd
    }

    /// FaceBuilder.cxx L633-661.
    pub(crate) fn make_loops(&mut self, brep: &mut BRep, ss: &mut WireEdgeSet) {
        // TopOpeBRepBuild_BlockBuilder& BB = myBlockBuilder;
        // NCollection_List<handle(TopOpeBRepBuild_Loop)>& LL =
        //   myLoopSet.ChangeListOfLoop();

        // Build blocks on elements of SS
        let mut bb = std::mem::take(&mut self.my_block_builder);
        bb.make_block(brep, ss);
        self.my_block_builder = bb;

        // make list of loop (LL) of the LoopSet
        // - on shapes of the ShapeSet (SS)
        // - on blocks of the BlockBuilder (BB)

        // Add shapes of SS as shape loops
        let ll = self.my_loop_set.change_list_of_loop();
        ll.clear();
        ss.init_shapes();
        while ss.more_shapes() {
            let s = ss.shape();
            let shape_loop = Arc::new(Loop::new_shape(&s));
            ll.push(shape_loop);
            ss.next_shape();
        }

        // Add blocks of BB as block loops
        let bb = &mut self.my_block_builder;
        bb.init_block();
        while bb.more_block() {
            let bi = bb.block_iterator();
            let block_loop = Arc::new(Loop::new_block(bi));
            ll.push(block_loop);
            bb.next_block();
        }
    }
}

/// FaceBuilder.cxx L85-116 — FUN_DetectVerticesOn1Edge: fills the map
/// <mapVon1edge> with vertices of <W> of 3d connexity 1.
fn fun_detect_vertices_on1edge(brep: &BRep, w: &Shape, map_von1e: &mut Vec<(Shape, Shape)>) {
    // OCCT: TopExp::MapShapesAndAncestors(W, TopAbs_VERTEX, TopAbs_EDGE,
    // mapVedges) — insertion order preserved (the IndexedDataMap walk).
    let mut map_vedges: Vec<(Shape, Vec<Shape>)> = Vec::new();
    let mut map_idx: std::collections::HashMap<(u64, u32), usize> =
        std::collections::HashMap::new();
    for w_child in shape_explore_children(brep, w, ShapeType::Edge) {
        for v in shape_explore_children(brep, &w_child, ShapeType::Vertex) {
            if !map_idx.contains_key(&shape_key(&v)) {
                map_idx.insert(shape_key(&v), map_vedges.len());
                map_vedges.push((v.clone(), Vec::new()));
            }
            let i = map_idx[&shape_key(&v)];
            map_vedges[i].1.push(w_child.clone());
        }
    }
    let nv = map_vedges.len();

    for i in 0..nv {
        let (v, loe) = (&map_vedges[i].0, &map_vedges[i].1);
        if v.orientation == Orientation::Internal {
            continue;
        }

        if loe.len() < 2 {
            // Keeping INTERNAL or EXTERNAL edges
            let e = &loe[0];
            let ori_e = e.orientation;
            if (ori_e == Orientation::Internal) || (ori_e == Orientation::External) {
                continue;
            }
            map_von1e.push((v.clone(), e.clone()));
        }
    }
}

/// OCCT BRep_Tool::IsClosed(S) (BRep_Tool.cxx L1707-1762) — the
/// single-argument closure test.  The FACE... wait, the OCCT branches are
/// SHELL / WIRE / EDGE / the Closed() flag fall-through.  GAP carrier of
/// the SHELL and WIRE branches (the batch consumes the EDGE branch):
/// translated per the OCCT body.
fn brep_tool_is_closed(brep: &BRep, the_shape: &Shape) -> bool {
    if the_shape.shape_type() == ShapeType::Shell {
        // BRep_Tool.cxx L1711-1733: every non-degenerated, non-internal
        // edge appears exactly twice.
        let mut seen: std::collections::HashMap<(u64, u32), i32> =
            std::collections::HashMap::new();
        let mut has_bound = false;
        for shell in shape_explore_children(brep, the_shape, ShapeType::Face) {
            let _ = shell;
            break;
        }
        let mut faces: Vec<Shape> = Vec::new();
        collect_faces(brep, the_shape, &mut faces);
        for f in &faces {
            for w in shape_explore_children(brep, f, ShapeType::Wire) {
                for e in shape_explore_children(brep, &w, ShapeType::Edge) {
                    let degenerated = brep.is_edge_degenerated(&e);
                    if degenerated
                        || e.orientation == Orientation::Internal
                        || e.orientation == Orientation::External
                    {
                        continue;
                    }
                    has_bound = true;
                    let c = seen.entry(shape_key(&e)).or_insert(0);
                    *c += 1;
                }
            }
        }
        return has_bound && seen.values().all(|&c| c == 2);
    } else if the_shape.shape_type() == ShapeType::Wire {
        // BRep_Tool.cxx L1734-1751: every vertex appears exactly twice.
        let mut seen: std::collections::HashMap<(u64, u32), i32> =
            std::collections::HashMap::new();
        let mut has_bound = false;
        for w in shape_explore_children(brep, the_shape, ShapeType::Wire) {
            let _ = w;
            break;
        }
        for w_child in shape_explore_children(brep, the_shape, ShapeType::Edge) {
            for v in shape_explore_children(brep, &w_child, ShapeType::Vertex) {
                if v.orientation == Orientation::Internal || v.orientation == Orientation::External
                {
                    continue;
                }
                has_bound = true;
                let c = seen.entry(shape_key(&v)).or_insert(0);
                *c += 1;
            }
        }
        return has_bound && seen.values().all(|&c| c == 2);
    } else if the_shape.shape_type() == ShapeType::Edge {
        // BRep_Tool.cxx L1752-1759.
        let (v_first, v_last) = top_exp_vertices(brep, the_shape);
        return !v_first.is_null() && shapes_same(&v_first, &v_last);
    }
    // TopoDS_Shape::Closed flag fall-through.
    brep.is_closed(the_shape)
}

/// Collect the faces under a shape (the OCCT TopExp_Explorer over
/// TopAbs_FACE — data-side walk).
fn collect_faces(brep: &BRep, s: &Shape, faces: &mut Vec<Shape>) {
    let _ = brep;
    match &*s.data {
        TShape::Shell(sd) => {
            for f in &sd.faces {
                faces.push(f.clone());
            }
        }
        TShape::Solid(sd) => {
            for sh in &sd.shells {
                collect_faces(brep, sh, faces);
            }
        }
        TShape::Compound(cs) | TShape::CompSolid(cs) => {
            for c in cs {
                collect_faces(brep, c, faces);
            }
        }
        TShape::Face(_) => faces.push(s.clone()),
        _ => {}
    }
}

/// FaceBuilder.cxx L124-200 — FUN_AnalyzemapVon1E.
fn fun_analyzemap_von1e(brep: &BRep, map_von1e: &[(Shape, Shape)], map_vv: &mut Vec<(Shape, Shape)>) -> i32 {
    let mut res = ISUNKNOWN;

    let nv = map_von1e.len();
    if nv == 0 {
        res = CLOSEDW;
    } else if nv == 1 {
        let e = &map_von1e[0].1;
        let eclosed = brep_tool_is_closed(brep, e);
        let dge = brep.is_edge_degenerated(e);
        if dge {
            res = ISVERTEX;
        } else if eclosed {
            res = CLOSEDW;
        } else {
            res = UNCLOSEDW;
        }
    } else {
        // Finding among all vertices, couple of vertices falling on same
        // geometry.  Filling up map <mapVV>, with (vi,vj), vi and vj are
        // on same point.
        let tol = CONFUSION;
        for i in 0..nv {
            let vi = &map_von1e[i].0;
            let pi = brep.vertex_position(vi);
            for j in (i + 1)..nv {
                let vj = &map_von1e[j].0;
                let pj = brep.vertex_position(vj);
                let same = pi.distance(pj) <= tol;
                if same {
                    map_vv.push((vi.clone(), vj.clone()));
                    map_vv.push((vj.clone(), vi.clone()));
                    break;
                }
            } // j
        } // i
        let nvv = map_vv.len();
        // (RM_HANGING is #undef'd in OCCT 8 — the #else form L188-196.)
        if nvv == nv {
            res = GCLOSEDW;
        } else {
            res = UNCLOSEDW;
        }
    }

    res
}

impl FaceBuilder {
    /// FaceBuilder.cxx L204-344 (DetectUnclosedWire; RM_HANGING is
    /// #undef'd in OCCT 8 — the #else purge branch L332-340 is the live
    /// one).
    pub(crate) fn detect_unclosed_wire(
        &mut self,
        brep: &mut BRep,
        map_vvsameg: &mut Vec<(Shape, Shape)>,
        map_von1edge: &mut Vec<(Shape, Shape)>,
    ) {
        map_vvsameg.clear();
        map_von1edge.clear();

        self.init_face();
        while self.more_face() {
            self.init_wire();
            while self.more_wire() {
                let isold = self.is_old_wire();
                if isold {
                    self.next_wire();
                    continue;
                }

                // TopoDS_Compound cmp; BB.MakeCompound(cmp) — the wire
                // content compound (the AddEdgeWire target).
                let cmp: Shape;
                {
                    let mut bb = BRepBuilder::new();
                    cmp = bb.make_compound(brep, Vec::new());
                }
                self.init_edge();
                while self.more_edge() {
                    let e = self.edge();
                    // AddEdgeWire(Edge(), cmp): BB.Add(W, E) — the
                    // compound append.
                    {
                        let b1 = BRepBuilder::new();
                        b1.add_to_compound(brep, cmp.clone(), e.clone());
                    }
                    self.next_edge();
                }
                let w = cmp;

                // <mapVon1E> binds vertices of connexity 1 attached to one
                // non-closed, non-degenerated edge.
                let mut map_von1e: Vec<(Shape, Shape)> = Vec::new();
                fun_detect_vertices_on1edge(brep, &w, &mut map_von1e);

                let mut map_vv: Vec<(Shape, Shape)> = Vec::new();
                let res = fun_analyzemap_von1e(brep, &map_von1e, &mut map_vv);

                if res == ISVERTEX {
                    self.next_wire();
                    continue;
                } else if res == CLOSEDW {
                    self.next_wire();
                    continue;
                } else if res == GCLOSEDW {
                    for (k, val) in &map_vv {
                        map_vvsameg.push((k.clone(), val.clone()));
                    }
                    for (k, val) in &map_von1e {
                        map_von1edge.push((k.clone(), val.clone()));
                    }
                } else if res == UNCLOSEDW {
                    // (RM_HANGING #undef'd — FaceBuilder.cxx L332-340.)
                    for ex_cur in shape_explore_children(brep, &w, ShapeType::Edge) {
                        let i = self.my_block_builder.element_of(&ex_cur);
                        self.my_block_builder.set_valid_index(i, false);
                    }
                }
                self.next_wire();
            } // MoreWire
            self.next_face();
        } // MoreFace
    }

    /// FaceBuilder.cxx L348-385.
    pub(crate) fn correct_gclosed_wire(
        &mut self,
        brep: &mut BRep,
        map_vvref: &[(Shape, Shape)],
        map_von1edge: &[(Shape, Shape)],
    ) {
        // prequesitory : edges described by <mapVon1Edge> are not closed,
        // not degenerated
        let nvv = map_vvref.len();
        for i in 0..nvv {
            let (vk, vrefk) = (&map_vvref[i].0, &map_vvref[i].1);
            let v = vk.clone();
            let vref = vrefk.clone();

            if shapes_same(&v, &vref) {
                continue;
            }

            // E = TopoDS::Edge(mapVon1Edge.FindFromKey(V)).
            let Some((_k, e)) = map_von1edge.iter().find(|(k, _)| shapes_same(k, &v)) else {
                // OCCT FindFromKey raises when unbound; the batch only
                // feeds maps built by DetectUnclosedWire (total by
                // construction).
                continue;
            };
            let e = e.clone();
            let parone = classify::brep_tool_parameter(brep, &v, &e);

            // E.Free(true).
            {
                let ed = brep.edge_mut_inplace(e.clone());
                ed.flags |= rcad_kernel::topods::tshape_flags::FREE;
            }
            // BB.Remove(E, V) — remove the vertex from the edge.
            {
                let ed = brep.edge_mut_inplace(e.clone());
                ed.my_shapes.retain(|c| !shapes_same(c, &v));
            }
            // newVref = TopoDS::Vertex(Vref.Oriented(V.Orientation())).
            let new_vref = shape_oriented(&vref, v.orientation);
            // BB.Add(E, newVref).
            {
                let ed = brep.edge_mut_inplace(e.clone());
                ed.my_shapes.push(new_vref.clone());
            }
            // TopOpeBRepDS_BuildTool BT; BT.Parameter(E, newVref, paronE).
            build_tool_parameter(brep, &e, &new_vref, parone);
        }
    }

    /// FaceBuilder.cxx L389-447.
    pub(crate) fn detect_pseudo_internal_edge(
        &mut self,
        brep: &mut BRep,
        map_e: &mut Vec<Shape>,
    ) {
        // TopoDS_Compound cmp; BB.MakeCompound(cmp).
        let cmp: Shape;
        {
            let mut bb = BRepBuilder::new();
            cmp = bb.make_compound(brep, Vec::new());
        }
        self.init_face();
        while self.more_face() {
            self.init_wire();
            while self.more_wire() {
                let isold = self.is_old_wire();
                if isold {
                    self.next_wire();
                    continue;
                }
                self.init_edge();
                while self.more_edge() {
                    let e = self.edge();
                    let b1 = BRepBuilder::new();
                    b1.add_to_compound(brep, cmp.clone(), e.clone());
                    self.next_edge();
                }
                self.next_wire();
            } // MoreWire
            self.next_face();
        } // MoreFace

        // TopExp::MapShapesAndAncestors(cmp, TopAbs_VERTEX, TopAbs_EDGE,
        // mapVOE).
        let mut map_voe: Vec<(Shape, Vec<Shape>)> = Vec::new();
        let mut map_idx: std::collections::HashMap<(u64, u32), usize> =
            std::collections::HashMap::new();
        for c in shape_explore_children(brep, &cmp, ShapeType::Shape) {
            if c.shape_type() != ShapeType::Edge {
                continue;
            }
            for v in shape_explore_children(brep, &c, ShapeType::Vertex) {
                if !map_idx.contains_key(&shape_key(&v)) {
                    map_idx.insert(shape_key(&v), map_voe.len());
                    map_voe.push((v.clone(), Vec::new()));
                }
                let i = map_idx[&shape_key(&v)];
                map_voe[i].1.push(c.clone());
            }
        }
        let nv = map_voe.len();

        map_e.clear();
        for i in 0..nv {
            let le = &map_voe[i].1;
            let ne = le.len();
            if ne == 2 {
                let e1 = &le[0];
                let e2 = &le[1];
                let same = shapes_same(e1, e2);
                let o1 = e1.orientation;
                let o2 = e2.orientation;
                let o1co2 = o1 == top_abs_complement(o2);

                if same && o1co2 {
                    map_e.push(e1.clone());

                    let ie1 = self.my_block_builder.element_of(e1);
                    self.my_block_builder.set_valid_index(ie1, false);

                    let ie2 = self.my_block_builder.element_of(e2);
                    self.my_block_builder.set_valid_index(ie2, false);
                }
            }
        }
    }
}

/// OCCT TopOpeBRepDS_BuildTool::Parameter(E, V, P) (BuildTool.cxx
/// L1239-1287) — the periodic-adjust tail + UpdateVertex(v, p, e, 0).
fn build_tool_parameter(brep: &mut BRep, e: &Shape, v: &Shape, p: f64) {
    let mut p = p;

    // 13/07/95: the periodic adjustment.
    if let Some((c, range)) = brep.edge_curve_world(e) {
        if c.is_periodic() {
            // OCCT C->Period() — the natural domain length (the circle /
            // ellipse period 2*PI is the default domain span).
            let domain = c.default_domain();
            let per = domain[1] - domain[0];
            let mut ov = Orientation::Forward;
            for vofe in shape_explore_children(brep, e, ShapeType::Vertex) {
                if shapes_same(&vofe, v) {
                    ov = vofe.orientation;
                    break;
                }
            }
            let f = range[0];
            if ov == Orientation::Reversed {
                if p < f {
                    let pp = {
                        // ElCLib::InPeriod(p, f, f + per).
                        let mut u = p;
                        let eps = 1.0e-12;
                        while u < f - eps {
                            u += per;
                        }
                        while u > f + per + eps {
                            u -= per;
                        }
                        u
                    };
                    p = pp;
                }
            }
        }
    }

    // myBuilder.UpdateVertex(v, p, e, 0).
    let mut b1 = BRepBuilder::new();
    b1.set_vertex_param(brep, e.clone(), v.clone(), p);
}

// =========================================================================
// OCCT TopOpeBRepBuild_Builder SplitFace1 arm — the D6 carrier methods on
// the TopOpeBRepBuildHBuilder (the myBuildTool calls map to the kernel
// BRepBuilder primitives at each statement, the hbuilder.rs precedent).
// =========================================================================

/// OCCT TopOpeBRepDS_DataStructure::HasShape(S) (DataStructure.cxx) — the
/// shape is present in the DS; the facade routes the lookup through the
/// BOPDS index (the add_shape registration identity).
fn ds_has_shape(ds: &TopOpeBRepDSHDataStructure, s: &Shape) -> bool {
    ds.bopds.index(s) >= 0
}

/// GAP carrier of TopOpeBRepDS_DataStructure::HasSameDomain /
/// SameDomain(S) (DataStructure.cxx): the ChFi3d D6 facade carries no
/// SameDomain table (ChFi3d fills none — the hbuilder.rs ToSplit gate
/// precedent), so the OCCT empty-table path is what returns.
fn ds_same_domain(_ds: &TopOpeBRepDSHDataStructure, _s: &Shape) -> Vec<Shape> {
    Vec::new()
}

/// GAP carrier of TopOpeBRepDS_DataStructure::HasNewSurface /
/// NewSurface(F) (DataStructure.cxx): the face -> new surface map
/// (myNewSurface) is not carried by the ChFi3d D6 facade, so the OCCT
/// unbound path (hns = false, the MakeFaces surface-update skip) is what
/// returns.
fn ds_has_new_surface(_ds: &TopOpeBRepDSHDataStructure, _f: &Shape) -> bool {
    false
}

impl TopOpeBRepBuildHBuilder {
    /// OCCT TopOpeBRepBuild_Builder::ToSplit (Builder.cxx L308-326).
    pub fn to_split(&self, ds: &TopOpeBRepDSHDataStructure, s: &Shape, to_build: TopAbsState) -> bool {
        let issplit = self.is_split(s, to_build);
        let hasgeom = ds.has_geometry(s);
        let hassame = !ds_same_domain(ds, s).is_empty();
        let tosplit = (!issplit) && (hasgeom || hassame);
        tosplit
    }

    /// OCCT TopOpeBRepBuild_Builder::Reverse (Builder.cxx L621-633).
    fn builder_reverse(to_build1: TopAbsState, to_build2: TopAbsState) -> bool {
        let rev;
        if to_build1 == TopAbsState::In && to_build2 == TopAbsState::In {
            rev = false;
        } else {
            rev = to_build1 == TopAbsState::In;
        }
        rev
    }

    /// OCCT TopOpeBRepBuild_Builder::Orient (Builder.cxx L637-640).
    fn builder_orient(ori: Orientation, reverse: bool) -> Orientation {
        if !reverse {
            ori
        } else {
            top_abs_complement(ori)
        }
    }

    /// OCCT TopOpeBRepBuild_Builder::Contains (Builder.cxx L836-847).
    fn builder_contains(s: &Shape, l: &[Shape]) -> bool {
        for it in l {
            if shapes_same(s, it) && s.orientation == it.orientation {
                return true;
            }
        }
        false
    }

    /// OCCT TopOpeBRepBuild_Builder::FindSameDomain (Builder.cxx L650-709).
    fn find_same_domain(
        &self,
        ds: &TopOpeBRepDSHDataStructure,
        l1: &mut Vec<Shape>,
        l2: &mut Vec<Shape>,
    ) {
        let mut nl1 = l1.len() as i32;
        let mut nl2 = l2.len() as i32;

        while nl1 > 0 || nl2 > 0 {
            let mut i = 0usize;
            let count1 = nl1 as usize;
            while i < l1.len() && i < count1 {
                let s1 = l1[i].clone();
                for s2 in ds_same_domain(ds, &s1) {
                    let found = Self::builder_contains(&s2, l2);
                    if !found {
                        // L2.Prepend(S2) — the OCCT prepend puts the shape
                        // at the front of the list.
                        l2.insert(0, s2);
                        nl2 += 1;
                    }
                }
                i += 1;
            }
            nl1 = 0;

            let mut i = 0usize;
            let count2 = nl2 as usize;
            while i < l2.len() && i < count2 {
                let s2 = l2[i].clone();
                for s1 in ds_same_domain(ds, &s2) {
                    let found = Self::builder_contains(&s1, l1);
                    if !found {
                        l1.insert(0, s1);
                        nl1 += 1;
                    }
                }
                i += 1;
            }
            nl2 = 0;
        }
    }

    /// OCCT TopOpeBRepBuild_Builder::ShapePosition (Builder.cxx L499-531).
    ///
    /// GAP carrier of myShapeClassifier.StateShapeShape (the
    /// TopOpeBRepDS_ShapeClassifier — TKBool/TopOpeBRepDS external
    /// machinery): the classifier is not part of the D6 carrier, so the
    /// OCCT UNKNOWN initial state is what returns (ShapePosition only
    /// runs when LS2 is non-empty, which cannot happen on this carrier —
    /// FindSameDomain is a no-op with no SameDomain table).
    fn shape_position(&self, _s: &Shape, _ls: &[Shape]) -> TopAbsState {
        let state = TopAbsState::Unknown;
        state
    }

    /// OCCT TopOpeBRepBuild_Builder::KeepShape (Builder.cxx L540-554).
    fn keep_shape(&self, s1: &Shape, ls2: &[Shape], to_build1: TopAbsState) -> bool {
        let mut keep = true;
        if !ls2.is_empty() {
            let pos2 = self.shape_position(s1, ls2);
            if pos2 != to_build1 {
                keep = false;
            }
        }
        keep
    }

    /// OCCT TopOpeBRepBuild_Builder::MarkSplit (Builder.cxx L330-374).
    fn mark_split(&mut self, s: &Shape, to_build: TopAbsState, bval: bool) {
        let key = shape_key(s);
        let p: Option<&mut ShapeStateMap> = match to_build {
            TopAbsState::Out => Some(&mut self.my_split_out),
            TopAbsState::In => Some(&mut self.my_split_in),
            TopAbsState::On => Some(&mut self.my_split_on),
            TopAbsState::Unknown => None,
        };
        let Some(p) = p else { return };

        let losos = p.entry(key).or_default();
        // TopOpeBRepDS_ListOfShapeOn1State::Split(Bval).
        losos.is_split = bval;
    }

    /// OCCT TopOpeBRepBuild_Builder::ChangeSplit (Builder.cxx L449-491) —
    /// the ListOnState reference of the mySplit(ToBuild) entry.
    fn change_split(&mut self, s: &Shape, to_build: TopAbsState) -> &mut Vec<Shape> {
        let key = shape_key(s);
        let p: Option<&mut ShapeStateMap> = match to_build {
            TopAbsState::Out => Some(&mut self.my_split_out),
            TopAbsState::In => Some(&mut self.my_split_in),
            TopAbsState::On => Some(&mut self.my_split_on),
            TopAbsState::Unknown => None,
        };
        let Some(p) = p else {
            // OCCT: return myEmptyShapeList — the static empty list (the
            // arm is unreachable on the carrier: to_build is IN/OUT/ON at
            // every call site); a leaked empty Vec carries the same
            // never-mutated semantics behind &mut.
            return Box::leak(Box::new(Vec::new()));
        };
        &mut p.entry(key).or_default().list_on_state
    }

    /// OCCT TopOpeBRepBuild_Builder::ChangeMerged (Merge.cxx L700-728) —
    /// the ListOnState reference of the myMerged(ToBuild) entry.
    fn change_merged(&mut self, s: &Shape, to_build: TopAbsState) -> &mut Vec<Shape> {
        let key = shape_key(s);
        let p: Option<&mut ShapeStateMap> = match to_build {
            TopAbsState::Out => Some(&mut self.my_merged_out),
            TopAbsState::In => Some(&mut self.my_merged_in),
            TopAbsState::On => Some(&mut self.my_merged_on),
            TopAbsState::Unknown => None,
        };
        let Some(p) = p else {
            // OCCT: return myEmptyShapeList — see change_split.
            return Box::leak(Box::new(Vec::new()));
        };
        &mut p.entry(key).or_default().list_on_state
    }

    /// OCCT TopOpeBRepBuild_Builder::AddIntersectionEdges (Builder.cxx
    /// L150-178).
    fn add_intersection_edges(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        a_face: &Shape,
        to_build1: TopAbsState,
        rev_ori1: bool,
        wes: &mut WireEdgeSet,
    ) {
        // The OCCT body consumes ToBuild1 at L164
        // (ori = FCurves.Orientation(ToBuild1)); the D6 facade transition
        // carries the IN-state orientation only, so the read re-routes to
        // transition.orientation_in() below.
        let _ = to_build1;
        // TopOpeBRepDS_CurveIterator FCurves =
        //   myDataStructure->FaceCurves(aFace) — the shape interferences
        //   of the face whose GeometryType is TopOpeBRepDS_CURVE
        //   (CurveIterator.cxx MatchInterference L20-27).
        let fcurves: Vec<&TopOpeBRepDSSurfaceCurveRef> = ds
            .shape_interferences_of(a_face)
            .iter()
            .filter_map(|i| match i {
                TopOpeBRepDSInterference::SurfaceCurve(sc)
                    if sc.kind_g == TopOpeBRepDSKind::Curve =>
                {
                    Some(sc)
                }
                _ => None,
            })
            .collect();
        for fcurve in fcurves {
            let ic = fcurve.index_g;
            let lnewe = self.new_edges(ic);
            for an_edge in lnewe {
                // OCCT: ori = FCurves.Orientation(ToBuild1) — the
                // interference transition orientation.  The D6 facade
                // transition carries the IN-state orientation only
                // (chfi3d_ds.rs TopOpeBRepDSTransition).
                let ori = fcurve.transition.orientation_in();
                let newori = Self::builder_orient(ori, rev_ori1);

                if newori == Orientation::External {
                    continue;
                }

                // myBuildTool.Orientation(anEdge, newori).
                let mut an_edge = an_edge;
                an_edge.orientation = newori;

                // myBuildTool.PCurve(aFace, anEdge, PC) — the pcurve store
                // under the face key (the BuildTool::PCurve
                // TopOpeBRepDS_SetThePCurve tail; the periodic-translate
                // adjustments of BuildTool.cxx L1200-1290 are D6-simplified
                // to the store, the build_faces precedent).
                if let Some(pc) = fcurve.pcurve.clone() {
                    let key_loc = brep.compose_pcurve_location(a_face.location, an_edge.location);
                    let key = (a_face.ptr_id(), key_loc);
                    // OCCT BRep_Builder.cxx UpdateCurves (L104-167, reached
                    // through TopOpeBRepDS_BuildTool::PCurve L1188-1239 ->
                    // BRep_Builder::UpdateEdge E,C,S,L,Tol L655-671): the new
                    // pcurve representation starts from the 2D curve's own
                    // range (COS->Range) and then INHERITS the range of the
                    // edge's 3D-curve representation (GC->Range for the
                    // IsCurve3D entry) whenever that range is finite.  Only
                    // an edge with no 3D curve keeps the 2D natural range.
                    let [mut a_f, mut a_l] = pc.default_domain();
                    let ed = brep.edge_mut_inplace(an_edge.clone());
                    if ed.curve.is_some() {
                        if !is_infinite_value(ed.range[0]) {
                            a_f = ed.range[0];
                        }
                        if !is_infinite_value(ed.range[1]) {
                            a_l = ed.range[1];
                        }
                    }
                    ed.pcurves.insert(key, (pc, a_f, a_l));
                }
                // WES.AddStartElement(anEdge).
                wes.add_start_element(brep, &an_edge);
            }
        }
    }

    /// OCCT TopOpeBRepBuild_Builder::FillFace (Builder.cxx L1951-1973).
    fn fill_face(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        f1: &Shape,
        to_build1: TopAbsState,
        lf2: &[Shape],
        to_build2: TopAbsState,
        wes: &mut WireEdgeSet,
        rev_ori: bool,
        my_list_of_face: &mut Vec<Shape>,
    ) {
        // myListOfFace = LF2.
        *my_list_of_face = lf2.to_vec();
        self.fill_shape(
            brep,
            ds,
            f1,
            to_build1,
            lf2,
            to_build2,
            wes,
            rev_ori,
            my_list_of_face,
        );
        // myListOfFace.Clear().
        my_list_of_face.clear();
    }

    /// OCCT TopOpeBRepBuild_Builder::FillShape (Builder.cxx L1876-1947).
    fn fill_shape(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        s1: &Shape,
        to_build1: TopAbsState,
        ls2: &[Shape],
        to_build2: TopAbsState,
        a_set: &mut WireEdgeSet,
        in_rev_ori: bool,
        my_list_of_face: &mut Vec<Shape>,
    ) {
        let mut rev_ori = in_rev_ori;
        let t = s1.shape_type();
        let (t1, t11) = match t {
            ShapeType::Face => (ShapeType::Wire, ShapeType::Edge),
            ShapeType::Solid | ShapeType::Shell => (ShapeType::Shell, ShapeType::Face),
            _ => (ShapeType::Compound, ShapeType::Compound),
        };

        // if the shape S1 is a SameDomain one, get its orientation
        // compared with the shape taken as reference.
        let hsd = !ds_same_domain(ds, s1).is_empty();
        if hsd {
            // TopOpeBRepDS_Config ssc = SameDomainOrientation(S1);
            // if (ssc == TopOpeBRepDS_DIFFORIENTED) RevOri = !RevOri.
            // Unreachable on the D6 carrier (no SameDomain table) — the
            // OCCT statement is preserved for form.
            let ssc_difforiented = false;
            if ssc_difforiented {
                rev_ori = !rev_ori;
            }
        }

        // work on a FORWARD shape <aShape>
        let mut a_shape = s1.clone();
        a_shape.orientation = Orientation::Forward;

        // Explore the SubShapes of type <t1>
        let ex1 = shape_explore_children(brep, &a_shape, t1);
        for a_sub_shape in &ex1 {
            if !ds_has_shape(ds, a_sub_shape) {
                // SubShape is not in DS : classify it with shapes of LS2
                let keep = self.keep_shape(a_sub_shape, ls2, to_build1);
                if keep {
                    let sub_ori = a_sub_shape.orientation;
                    let newori = Self::builder_orient(sub_ori, rev_ori);
                    let a_sub_shape = shape_oriented(a_sub_shape, newori);
                    a_set.add_shape(brep, &a_sub_shape);
                }
            } else {
                // SubShape has geometry : split the <t11> SubShapes of the
                // SubShape
                let ex11 = shape_explore_children(brep, a_sub_shape, t11);
                self.split_shapes(
                    brep,
                    ds,
                    &ex11,
                    to_build1,
                    to_build2,
                    a_set,
                    rev_ori,
                    my_list_of_face,
                );
            }
        } // exploration of SubShapes of type <t1> of shape <S1>
    }

    /// OCCT TopOpeBRepBuild_Builder::SplitShapes (Builder.cxx L1713-1872).
    fn split_shapes(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        ex: &[Shape],
        to_build1: TopAbsState,
        to_build2: TopAbsState,
        a_set: &mut WireEdgeSet,
        rev_ori: bool,
        my_list_of_face: &mut Vec<Shape>,
    ) {
        for a_shape in ex {
            let a_shape = a_shape.clone();

            // compute new orientation <newori> to give to the new shapes
            let newori = Self::builder_orient(a_shape.orientation, rev_ori);

            let t = a_shape.shape_type();

            if t == ShapeType::Solid || t == ShapeType::Shell {
                // OCCT: SplitSolid(aShape, ToBuild1, ToBuild2) — the solid
                // arm is not part of this batch (merge_solid carries the
                // SplitSolid semantics); unreachable from SplitFace1's
                // face fill.
                continue;
            } else if t == ShapeType::Face {
                self.split_face1(brep, ds, &a_shape, to_build1, to_build2);
            } else if t == ShapeType::Edge {
                // OCCT L1741: SplitEdge(aShape, ToBuild1, ToBuild2) — the
                // SplitEdge1 body (Builder.cxx L941-1079): the pave walk +
                // the 1d area walk select the per-state split pieces
                // (hbuilder_face/split_edge.rs).
                self.split_edge(brep, ds, &a_shape, to_build1, to_build2);
            } else {
                continue;
            }

            if self.is_split(&a_shape, to_build1) {
                //----------------------- IFV
                let mut is_lson = false;
                //----------------------- IFV
                let mut ls = self.splits(&a_shape, to_build1);
                //----------------------- IFV
                if t == ShapeType::Edge && to_build1 == TopAbsState::In && ls.is_empty() {
                    let lson = self.splits(&a_shape, TopAbsState::On);
                    ls = lson;
                    is_lson = true;
                }
                //----------------------- IFV
                for new_shape in &ls {
                    let mut new_shape = new_shape.clone();
                    // myBuildTool.Orientation(newShape, newori).
                    new_shape.orientation = newori;
                    //----------------------- IFV
                    if is_lson {
                        let mut add = true;
                        if !my_list_of_face.is_empty() {
                            // 2d pur
                            add = self.keep_shape(&new_shape, my_list_of_face, to_build1);
                        }
                        if add {
                            a_set.add_start_element(brep, &new_shape);
                        }
                    } else {
                        //----------------------- IFV
                        a_set.add_start_element(brep, &new_shape);
                    }
                }
            } else {
                // aShape n'a pas de devenir de split par ToBuild1
                // on construit les parties ToBuild1 de aShape (de S1)
                let mut add = true;
                let mut testkeep = false;
                let isedge = t == ShapeType::Edge;
                let hs = ds_has_shape(ds, &a_shape);
                let hg = ds.has_geometry(&a_shape);

                testkeep = isedge && hs && (!hg);

                // xpu010399 : USA60299 (!hs)&&(!hg), but vertex on bound is
                // touched (v7) -> testkeep
                let mut istouched = isedge && (!hs) && (!hg);
                if istouched {
                    istouched = fun_touched(brep, ds, &a_shape);
                }
                testkeep = testkeep || istouched;


                if testkeep {
                    if !my_list_of_face.is_empty() {
                        // 2d pur
                        let keep = self.keep_shape(&a_shape, my_list_of_face, to_build1);
                        add = keep;
                    } else {
                        // 3d : on classifie en solide uniqt si E dans la DS
                        // et E a ete purgee de ses interfs car en bout
                        //
                        // OCCT L1827-1835: sol = (STATIC_SOLIDINDEX == 1) ?
                        //   myShape2 : myShape1 — the Builder members
                        //   myShape1 / myShape2 are set by MergeShapes /
                        //   Perform(HDS,S1,S2), neither runs on the D6
                        //   facade; the handles are null, so the OCCT
                        //   sol.IsNull() branch (L1857-1860, add = true)
                        //   is the live path.
                        let static_solidindex = STATIC_SOLIDINDEX.load(
                            std::sync::atomic::Ordering::Relaxed,
                        );
                        let _my_shape1 = Shape::null();
                        let _my_shape2 = Shape::null();
                        let sol = if static_solidindex == 1 {
                            _my_shape2.clone()
                        } else {
                            _my_shape1.clone()
                        };
                        if !sol.is_null() {
                            // OCCT L1838-1853: BRepClass3d_SolidClassifier
                            // SCL(sol, P3D, tol3d) at
                            // par = (1 - 0.127956477) * first + 0.127956477
                            // * last; add = (state == ToBuild1).  Unreachable
                            // on this carrier (sol is null); the statement is
                            // preserved for form.
                            let tt = 0.127956477;
                            let _ = tt;
                            let _ = to_build1;
                            add = false;
                        } else {
                            // sol.IsNull
                            add = true;
                        }
                    }
                }
                if add {
                    let a_shape = shape_oriented(&a_shape, newori);
                    a_set.add_element(brep, &a_shape);
                }
            }
        } // Ex.More
    }

    /// OCCT TopOpeBRepBuild_Builder::SplitFace (Builder.cxx L1155-1167)
    /// — the GetcontextSF2 debug switch is compiled out, so SplitFace1 is
    /// the body.
    pub fn split_face(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        foriented: &Shape,
        to_build1: TopAbsState,
        to_build2: TopAbsState,
    ) {
        self.split_face1(brep, ds, foriented, to_build1, to_build2);
    }

    /// OCCT TopOpeBRepBuild_Builder::SplitFace1 (Builder.cxx L1171-1300).
    pub fn split_face1(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        foriented: &Shape,
        to_build1: TopAbsState,
        to_build2: TopAbsState,
    ) {
        //                              process  connect  connect
        // operation tobuild1 tobuild2  face F   to 1     to 2
        // --------- -------- --------  -------  -------  -------
        // common    IN       IN        yes      yes      yes
        // fuse      OUT      OUT       yes      yes      yes
        // cut 1-2   OUT      IN        yes      yes      no
        // cut 2-1   IN       OUT       yes      yes      no
        //
        let tosplit = self.to_split(ds, foriented, to_build1);
        if !tosplit {
            return;
        }

        let mut rev_ori1 = Self::builder_reverse(to_build1, to_build2);
        let mut rev_ori2 = Self::builder_reverse(to_build2, to_build1);
        let connect_to1 = true;
        let connect_to2 = false;

        // work on a FORWARD face <Fforward>
        let mut fforward = foriented.clone();
        fforward.orientation = Orientation::Forward;

        // build the list of faces to split : LF1, LF2
        let mut lf1: Vec<Shape> = Vec::new();
        let mut lf2: Vec<Shape> = Vec::new();
        lf1.push(fforward.clone());
        self.find_same_domain(ds, &mut lf1, &mut lf2);
        let n1 = lf1.len();
        let n2 = lf2.len();

        // SplitFace on a face having other same domained faces on the
        // other shape : do not reverse orientation of faces in FillFace
        if n2 == 0 {
            rev_ori1 = false;
        }
        if n1 == 0 {
            rev_ori2 = false;
        }

        // Create an edge set <WES> connected by vertices
        // ----------------------------------------------
        let mut wes = WireEdgeSet::new(&fforward);

        let mut my_list_of_face: Vec<Shape> = Vec::new();

        for fcur in lf1.clone() {
            self.fill_face(
                brep,
                ds,
                &fcur,
                to_build1,
                &lf2,
                to_build2,
                &mut wes,
                rev_ori1,
                &mut my_list_of_face,
            );
        }

        for fcur in lf2.clone() {
            self.fill_face(
                brep,
                ds,
                &fcur,
                to_build2,
                &lf1,
                to_build1,
                &mut wes,
                rev_ori2,
                &mut my_list_of_face,
            );
        }

        // Add the intersection edges to edge set WES
        // -----------------------------------------
        self.add_intersection_edges(brep, ds, &fforward, to_build1, rev_ori1, &mut wes);

        // Create a Face Builder FBU
        // ------------------------
        let mut fbu = FaceBuilder::new();
        fbu.init_face_builder(brep, &mut wes, &fforward, false); // forceclass = False

        // Build the new faces
        // -------------------
        // OCCT: NCollection_List<TopoDS_Shape>& FaceList =
        //   ChangeMerged(Fforward, ToBuild1); MakeFaces(Fforward, FBU,
        //   FaceList).  The OCCT list reference becomes: MakeFaces fills a
        //   local list, the ChangeMerged entry receives it below (the same
        //   data flow).
        let mut face_list: Vec<Shape> = Vec::new();
        self.make_faces(brep, ds, &fforward, &mut fbu, &mut face_list);
        self.change_merged(&fforward, to_build1).clone_from(&face_list);

        // connect new faces as faces built <ToBuild1> on LF1 faces
        // --------------------------------------------------------
        for fcur in &lf1 {
            self.mark_split(fcur, to_build1, true);
            let fl = self.change_split(fcur, to_build1);
            if connect_to1 {
                fl.clone_from(&face_list);
            }
        }

        // connect new faces as faces built <ToBuild2> on LF2 faces
        // --------------------------------------------------------
        for fcur in &lf2 {
            self.mark_split(fcur, to_build2, true);
            let fl = self.change_split(fcur, to_build2);
            if connect_to2 {
                fl.clone_from(&face_list);
            }
        }
    } // SplitFace1

    /// OCCT TopOpeBRepBuild_Builder::MakeFaces (Merge.cxx L435-512).
    fn make_faces(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        a_face: &Shape,
        fabu: &mut FaceBuilder,
        l: &mut Vec<Shape>,
    ) {
        let hashds = true; // myDataStructure is stored by Perform
        let mut b1 = BRepBuilder::new();

        fabu.init_face();
        while fabu.more_face() {
            // myBuildTool.CopyFace(aFace, newFace) — Fou =
            // Fin.EmptyCopied() (BuildTool.cxx CopyFace).
            let new_face = brep.empty_copied(a_face);
            let mut hns = false;
            if hashds {
                // BDS.HasNewSurface(aFace) — the D6 facade carries no
                // myNewSurface map (GAP carrier ds_has_new_surface), the
                // OCCT unbound path.
                hns = ds_has_new_surface(ds, a_face);
                if hns {
                    // myBuildTool.UpdateSurface(newFace, SU) — unreachable
                    // on the carrier (hns false); the statement is
                    // preserved for form.
                }
            }

            fabu.init_wire();
            while fabu.more_wire() {
                let isold = fabu.is_old_wire();
                let mut new_wire: Shape;
                if isold {
                    new_wire = fabu.old_wire();
                } else {
                    // myBuildTool.MakeWire(newWire).
                    new_wire = b1.make_wire(brep);
                    fabu.init_edge();
                    while fabu.more_edge() {
                        let e = fabu.edge();
                        if hns {
                            // myBuildTool.UpdateSurface(E, aFace, newFace)
                            // — the pcurve transfer of the old face key to
                            // the new face key (BuildTool.cxx
                            // UpdateSurface(E, oldF, newF)).  Unreachable
                            // on the carrier (hns is false — see the GAP
                            // carrier ds_has_new_surface); the statement
                            // is preserved for form.
                        }
                        // myBuildTool.AddWireEdge(newWire, E).
                        b1.add_to_wire(brep, new_wire.clone(), e);
                        fabu.next_edge();
                    }
                }
                //----------- IFV
                if !isold {
                    // OCCT L489-505: BRepCheck_Analyzer bca(newWire, false);
                    // if (!bca.IsValid()) { newWire.Free(true);
                    // CorrectUnclosedWire(newWire); ... BadOrientationOfSubshape
                    // -> CorrectEdgeOrientation(newWire); }.
                    //
                    // GAP carrier of BRepCheck_Analyzer (TKBRep/BRepCheck —
                    // external untranslated dependency): the analyzer is not
                    // translated; the OCCT invalid-wire path is the one taken
                    // (FaceBuilder block wires are the population this
                    // correction exists for) — the BadOrientationOfSubshape
                    // status scan has no carrier (no per-shape status list)
                    // and stays unreached.
                    let bca_is_valid = false;
                    if !bca_is_valid {
                        // newWire.Free(true).
                        {
                            let wd = brep.wire_mut(new_wire.clone());
                            wd.flags |= rcad_kernel::topods::tshape_flags::FREE;
                        }
                        correct_unclosed_wire(brep, &mut new_wire);
                        // (the BRepCheck_BadOrientationOfSubshape scan:
                        // no carrier — see the GAP note)
                    }
                }
                // myBuildTool.Closed(newWire, true).  // NYI : check exact du
                // caractere closed du wire
                {
                    let wd = brep.wire_mut(new_wire.clone());
                    wd.flags |= rcad_kernel::topods::tshape_flags::CLOSED;
                }
                // myBuildTool.AddFaceWire(newFace, newWire).
                bb_add_face_wire(brep, &new_face, &new_wire);
                fabu.next_wire();
            }

            fabu.next_face();

            // L.Append(newFace).
            l.push(new_face);
        }
    }
}

/// OCCT Builder.cxx L1702-1709 — FUN_touched.
fn fun_touched(brep: &BRep, ds: &TopOpeBRepDSHDataStructure, eor: &Shape) -> bool {
    let (vf, vl) = top_exp_vertices(brep, eor);
    let hvf = ds.bopds.index(&vf) >= 0;
    let hvl = ds.bopds.index(&vl) >= 0;
    hvf || hvl
}

/// OCCT Merge.cxx L48-136 — CorrectEdgeOrientation (the IFV wire repair).
fn correct_edge_orientation(brep: &mut BRep, a_wire: &mut Shape) {
    let mut an_edge_list: Vec<Shape> = Vec::new();
    let mut an_aux_list: Vec<Shape> = Vec::new();
    let mut a_true_edge_list: Vec<Shape> = Vec::new();
    let mut append = true;

    // TopoDS_Iterator tdi(aWire, false, false) — the direct children.
    an_edge_list.extend(shape_explore_children(brep, a_wire, ShapeType::Edge));

    let mut n = an_edge_list.len() as i32;
    if n <= 1 {
        return;
    }

    let mut idx = 0usize;
    let mut an_cur_edge = an_edge_list[idx].clone();
    // TopExp::Vertices(E, vf, vl, true) — CumOri = true.
    let (mut vf, mut vl) = top_exp_vertices_cumori(brep, &an_cur_edge);
    a_true_edge_list.push(an_cur_edge.clone());
    idx += 1;

    while n > 0 && append {
        append = false;
        while idx < an_edge_list.len() {
            an_cur_edge = an_edge_list[idx].clone();
            let (v1f, v1l) = top_exp_vertices_cumori(brep, &an_cur_edge);
            if shapes_same(&v1f, &vl) {
                a_true_edge_list.push(an_cur_edge.clone());
                vl = v1l;
                append = true;
                idx += 1;
                continue;
            }
            if shapes_same(&v1l, &vf) {
                a_true_edge_list.push(an_cur_edge.clone());
                vf = v1f;
                append = true;
                idx += 1;
                continue;
            }

            if shapes_same(&v1l, &vl) {
                let an_rev_edge = shape_oriented(&an_cur_edge, Orientation::Reversed);
                a_true_edge_list.push(an_rev_edge);
                vl = v1f;
                append = true;
                idx += 1;
                continue;
            }
            if shapes_same(&v1f, &vf) {
                let an_rev_edge = shape_oriented(&an_cur_edge, Orientation::Reversed);
                a_true_edge_list.push(an_rev_edge);
                vf = v1l;
                append = true;
                idx += 1;
                continue;
            }

            an_aux_list.push(an_cur_edge.clone());
            idx += 1;
        }

        an_edge_list = std::mem::take(&mut an_aux_list);
        // anAuxList.Clear() — something wrong in Assign when list contains
        // 1 element.
        n = an_edge_list.len() as i32;
        idx = 0;
    }

    if n > 0 {
        a_true_edge_list.extend(an_edge_list.iter().cloned());
    }

    // aWire.Nullify(); BB.MakeWire(TopoDS::Wire(aWire));
    // for each edge: BB.Add(aWire, edge).
    let mut bb = BRepBuilder::new();
    let new_wire = bb.make_wire(brep);
    for e in &a_true_edge_list {
        bb.add_to_wire(brep, new_wire.clone(), e.clone());
    }
    *a_wire = new_wire;
}

/// OCCT TopExp::Vertices(E, Vfirst, Vlast, CumOri = true)
/// (TopExp.cxx L182-210): the iterator composes the edge orientation into
/// the stored vertices, then Vfirst is the FORWARD one and Vlast the
/// REVERSED one — for a REVERSED edge the canonical last vertex comes out
/// FORWARD and the canonical first REVERSED.
fn top_exp_vertices_cumori(brep: &BRep, e: &Shape) -> (Shape, Shape) {
    let (cf, cl) = (brep.first_vertex(e), brep.last_vertex(e));
    if e.orientation == Orientation::Reversed {
        (
            shape_oriented(&cl, Orientation::Forward),
            shape_oriented(&cf, Orientation::Reversed),
        )
    } else {
        (
            shape_oriented(&cf, Orientation::Forward),
            shape_oriented(&cl, Orientation::Reversed),
        )
    }
}

/// OCCT Merge.cxx L138-170 — CorrectUnclosedWire (the IFV wire repair).
fn correct_unclosed_wire(brep: &mut BRep, a_wire: &mut Shape) {
    // TopoDS_Iterator tdi(aWire, false, false): ed.NbChildren() <= 1 ->
    // BB.Remove(aWire, ed).
    let mut to_remove: Vec<Shape> = Vec::new();
    let edges_now = shape_explore_children(brep, a_wire, ShapeType::Edge);
    for ed in &edges_now {
        // ed.NbChildren(): the edge's subshape (vertex) count.
        let nbv = brep.edge(ed.clone()).my_shapes.len();
        if nbv <= 1 {
            to_remove.push(ed.clone());
        }
    }

    // TopExp::MapShapesAndAncestors(aWire, TopAbs_VERTEX, TopAbs_EDGE,
    // VElists); the single-edge vertices lose their edge.
    let mut ve_lists: Vec<(Shape, Vec<Shape>)> = Vec::new();
    let mut map_idx: std::collections::HashMap<(u64, u32), usize> = std::collections::HashMap::new();
    let edges_kept = shape_explore_children(brep, a_wire, ShapeType::Edge);
    for e in &edges_kept {
        if to_remove.iter().any(|r| shapes_same(r, e)) {
            continue;
        }
        for v in shape_explore_children(brep, e, ShapeType::Vertex) {
            if !map_idx.contains_key(&shape_key(&v)) {
                map_idx.insert(shape_key(&v), ve_lists.len());
                ve_lists.push((v.clone(), Vec::new()));
            }
            let i = map_idx[&shape_key(&v)];
            ve_lists[i].1.push(e.clone());
        }
    }

    for (_v, elist) in &ve_lists {
        if elist.len() == 1 {
            let an_edge = &elist[0];
            // BB.Remove(aWire, anEdge).
            to_remove.push(an_edge.clone());
        }
    }

    if !to_remove.is_empty() {
        let wd = brep.wire_mut(a_wire.clone());
        wd.edges
            .retain(|e| !to_remove.iter().any(|r| shapes_same(r, e)));
    }
}

/// The interference reference type AddIntersectionEdges walks (the
/// SurfaceCurve quadruplet).
type TopOpeBRepDSSurfaceCurveRef = super::chfi3d_ds::TopOpeBRepDSSurfaceCurveInterference;
