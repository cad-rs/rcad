//! OCCT TopOpeBRepBuild SplitEdge1 chain (the SplitEdge arm of
//! TopOpeBRepBuild_Builder, Builder.cxx L925-1079 + Merge.cxx L516-622 +
//! Builder.cxx L1991-2055) — continuation file of `fillet::hbuilder_face`
//! (the D6 carrier of the TKBool TopOpeBRepBuild classes as components of
//! the TopOpeBRepBuild_HBuilder translation):
//!
//! - `TopOpeBRepBuild_Pave` (Pave.hxx L30-58)
//! - `TopOpeBRepBuild_PaveSet` (PaveSet.cxx L35-491)
//! - `TopOpeBRepBuild_PaveClassifier` (PaveClassifier.cxx L35-393)
//! - `TopOpeBRepBuild_Area1dBuilder::InitAreaBuilder` (Area1dBuilder.cxx
//!   L67-298) + the AreaBuilder base helpers (AreaBuilder.cxx L59-121,
//!   L334-429)
//! - `TopOpeBRepBuild_EdgeBuilder` (EdgeBuilder.cxx L23-105)
//! - `Builder::FillVertexSet / FillVertexSetOnValue` (Builder.cxx
//!   L1991-2055)
//! - `Builder::SplitEdge / SplitEdge1` (Builder.cxx L925-937 / L941-1079)
//! - `Builder::MakeEdges` (Merge.cxx L516-622)
//!
//! The chain is the per-state selection of the edge split pieces: the
//! pave walk (FillVertexSetOnValue) orients each interference pave with
//! the transition orientation of the state being built
//! (TopOpeBRepDS_Transition::Orientation(S, T)), and the 1d area walk
//! (Area1dBuilder over the PaveSet with the PaveClassifier) groups each
//! interference pave with the bound pave that bounds the piece valid for
//! that state — the paves left in boundaryloops are never built, so
//! Splits(E, ToBuild) carries only the pieces on the state (the OCCT
//! guarantee SplitShapes L1755 relies on).
//!
//! Architecture differences (Rust <-> C++), per the D6 carrier notes:
//! - OCCT pave handles (`occ::handle<TopOpeBRepBuild_Pave>`) are shared
//!   across the area / boundaryloops lists; rcad clones the pave value
//!   (the Shape handle inside is the shared payload; no pave is mutated
//!   after insertion, so the aliasing semantics are preserved).
//! - OCCT reads the BRep through `BRep_Tool` globals; rcad threads the
//!   owning pool (`brep: &BRep`) into PaveSet::Prepare and the
//!   PaveClassifier ctor (the same read-only pool the caller holds).
//! - `BRep_Builder` empty-copy + AddEdgeVertex (BuildTool.cxx L293-316 /
//!   L915-918) maps to the `add_tedge` primitive + the in-place TEdge
//!   edits (the hbuilder.rs split_ds_edges precedent).
//! - The OCCT one-surviving-vertex piece append (Merge.cxx L617 appends a
//!   one-vertex edge) cannot be represented by the rcad TEdge (two bound
//!   vertices required); on this carrier every area pairs one interference
//!   pave with one bound pave, both non-EXTERNAL after PaveSet::Prepare,
//!   so the arm is unreachable — the piece is skipped with the note.

#![allow(dead_code)]

use rcad_kernel::core::precision::PCONFUSION;
use rcad_kernel::geom::CurveEval as _;
use rcad_kernel::math::el::in_period;
use rcad_kernel::topods::{BRep, Orientation, Shape, TShape};

use super::super::chfi3d_builder_0::brep_tool_parameter;
use super::super::chfi3d_builder_2::TopAbsState;
use super::super::chfi3d_ds::{TopOpeBRepDSHDataStructure, TopOpeBRepDSInterference, TopOpeBRepDSKind};
use super::super::hbuilder::TopOpeBRepBuildHBuilder;
use super::classify::{
    shape_oriented, shapes_equal, shapes_same, top_exp_vertices, LoopEnum,
};

/// BRep_Tool::Pnt(V) read off the TShape data.
fn vertex_point_of(s: &Shape) -> glam::DVec3 {
    match &*s.data {
        TShape::Vertex(vd) => vd.point,
        _ => glam::DVec3::ZERO,
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Pave (Pave.hxx L30-58): the vertex + parameter
// pair, Loop subclass (IsShape() = the bound flag).
// =========================================================================
#[derive(Debug, Clone)]
pub(crate) struct Pave {
    /// OCCT: TopoDS_Shape myVertex.
    pub(crate) my_vertex: Shape,
    /// OCCT: double myParam.
    pub(crate) my_param: f64,
    /// OCCT: bool myIsShape (the Loop::IsShape override — true for the
    /// edge bound vertices).
    pub(crate) my_is_shape: bool,
    /// OCCT: bool myHasSameDomain (ChFi3d fills no SameDomain table —
    /// stays false on this carrier).
    #[allow(dead_code)]
    pub(crate) my_has_same_domain: bool,
    /// OCCT: TopoDS_Shape mySameDomain.
    #[allow(dead_code)]
    pub(crate) my_same_domain: Shape,
}

impl Pave {
    /// Pave.hxx L34-37: V = vertex, P = parameter, bound = true if V is an
    /// old vertex / false if V is a new vertex.
    pub(crate) fn new(v: &Shape, p: f64, bound: bool) -> Self {
        Pave {
            my_vertex: v.clone(),
            my_param: p,
            my_is_shape: bound,
            my_has_same_domain: false,
            my_same_domain: Shape::null(),
        }
    }

    /// Pave.hxx L62-66.
    pub(crate) fn vertex(&self) -> &Shape {
        &self.my_vertex
    }

    /// Pave.hxx L74-76: Parameter().
    pub(crate) fn parameter(&self) -> f64 {
        self.my_param
    }

    /// Pave.hxx L55 (the Loop::IsShape override).
    pub(crate) fn is_shape(&self) -> bool {
        self.my_is_shape
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_PaveSet (PaveSet.cxx L35-491).
// =========================================================================
pub(crate) struct PaveSet {
    /// OCCT: TopoDS_Edge myEdge.
    my_edge: Shape,
    /// OCCT: TopOpeBRepBuild_ListOfPave myVertices (the myVerticesIt
    /// iterator becomes the index cursor).
    my_vertices: Vec<Pave>,
    my_vertices_it: usize,
    /// OCCT: bool myHasEqualParameters.
    my_has_equal_parameters: bool,
    /// OCCT: double myEqualParameters.
    my_equal_parameters: f64,
    /// OCCT: bool myClosed.
    my_closed: bool,
    /// OCCT: bool myPrepareDone.
    my_prepare_done: bool,
    /// OCCT: bool myRemovePV.
    my_remove_pv: bool,
}

impl PaveSet {
    /// PaveSet.cxx L35-42.
    pub(crate) fn new(e: &Shape) -> Self {
        PaveSet {
            my_edge: e.clone(),
            my_vertices: Vec::new(),
            my_vertices_it: 0,
            my_has_equal_parameters: false,
            my_equal_parameters: 0.0,
            my_closed: false,
            my_prepare_done: false,
            my_remove_pv: true,
        }
    }

    /// PaveSet.cxx L46-48: RemovePV(B).
    #[allow(dead_code)]
    pub(crate) fn remove_pv(&mut self, b: bool) {
        self.my_remove_pv = b;
    }

    /// PaveSet.cxx L53-58: Append(PV).
    pub(crate) fn append(&mut self, pv: Pave) {
        self.my_vertices.push(pv);
        self.my_prepare_done = false;
    }

    /// PaveSet.cxx L64-129: SortPave — sort the paves on Parameter() value,
    /// then rotate so the list starts on the first FORWARD pave ("tete =
    /// FORWARD"); the leading REVERSED paves wrap to the end.
    fn sort_pave(list: Vec<Pave>, sorted_list: &mut Vec<Pave>) {
        let n_pv = list.len();
        let mut taken = vec![false; n_pv];

        let mut sorted: Vec<Pave> = Vec::with_capacity(n_pv);
        for _ in 1..=n_pv {
            let mut parmin = f64::INFINITY;
            let mut i_pv = 0usize;
            let mut pv1: Option<&Pave> = None;
            for (itest, pv2) in list.iter().enumerate() {
                if !taken[itest] {
                    let par = pv2.parameter();
                    if par < parmin {
                        parmin = par;
                        pv1 = Some(pv2);
                        i_pv = itest;
                    }
                }
            }
            if let Some(pv1) = pv1 {
                sorted.push(pv1.clone());
            }
            taken[i_pv] = true;
        }

        // tete = FORWARD (PaveSet.cxx L92-128).
        let mut found = false;
        let mut l1: Vec<Pave> = Vec::new();
        let mut l2: Vec<Pave> = Vec::new();
        for pv in &sorted {
            if !found {
                let o = pv.vertex().orientation;
                if o == Orientation::Forward {
                    found = true;
                    l1.push(pv.clone());
                } else {
                    l2.push(pv.clone());
                }
            } else {
                l1.push(pv.clone());
            }
        }

        sorted_list.clear();
        sorted_list.append(&mut l1);
        sorted_list.append(&mut l2);
    }

    /// PaveSet.cxx L131-140: FUN_islook — |P(v1) - P(v2)| > 1.e-8.
    fn fun_islook(brep: &BRep, e: &Shape) -> bool {
        let (v1, v2) = top_exp_vertices(brep, e);
        let p1 = vertex_point_of(&v1);
        let p2 = vertex_point_of(&v2);
        let dp1p2 = p1.distance(p2);
        dp1p2.abs() > 1.0e-8
    }

    /// PaveSet.cxx L144-319: Prepare — add the edge vertices to the list of
    /// interferences; if an edge vertex VE is already in the list as
    /// interference VI: do not add VE in the list; if VI is INTERNAL, set
    /// VI orientation to VE orientation; remove VI from the list if VI is
    /// EXTERNAL or VE and VI have opposite orientations.
    pub(crate) fn prepare(&mut self, brep: &BRep) {
        if self.my_prepare_done {
            return;
        }

        // L183: bool isEd = BRep_Tool::Degenerated(myEdge).
        let is_ed = self.my_edge.as_edge().expect("PaveSet edge").degenerated;
        let mut edge_vertex_count = 0usize;

        if self.my_remove_pv {
            // L186-189: TopExp_Explorer EVexp(myEdge, TopAbs_VERTEX) — the
            // bound vertices of the FORWARD edge: first FORWARD, last
            // REVERSED.
            let bounds: Vec<(Shape, Orientation)> = {
                let (vfirst, vlast) = top_exp_vertices(brep, &self.my_edge);
                vec![(vfirst, Orientation::Forward), (vlast, Orientation::Reversed)]
            };
            for (ve, ve_ori) in &bounds {
                // L194-195: VEori / VEbound.
                let ve_bound = *ve_ori == Orientation::Forward || *ve_ori == Orientation::Reversed;

                let mut edge_vertex_index = 0usize;
                let mut add_ve = true;

                // L200: bool add = false; (ofv)
                let mut add = false;

                let mut it = 0usize;
                while it < self.my_vertices.len() {
                    // L205-209: skip the edge vertices inserted at the head
                    // of the list.
                    edge_vertex_index += 1;
                    if edge_vertex_index <= edge_vertex_count {
                        it += 1;
                        continue;
                    }

                    // L212-228: PV = parametrized vertex, VI = interference
                    // vertex; hasVSD stays false on this carrier (ChFi3d
                    // fills no SameDomain table).
                    let vi_same_ve = shapes_same(&self.my_vertices[it].my_vertex, ve);
                    let vsd_same_ve = false;
                    let same_vertex_processing = (vi_same_ve || vsd_same_ve) && !is_ed;

                    if same_vertex_processing {
                        let vi_ori = self.my_vertices[it].my_vertex.orientation;
                        // L233: if (VEbound || vsd_sameve).
                        if ve_bound || vsd_same_ve {
                            match vi_ori {
                                Orientation::External => {
                                    // L238-240: remove VI.
                                    self.my_vertices.remove(it);
                                }
                                Orientation::Internal => {
                                    // L242-244: VI.Orientation(VEori).
                                    self.my_vertices[it].my_vertex.orientation = *ve_ori;
                                }
                                Orientation::Forward | Orientation::Reversed => {
                                    // L246-259 (ofv): remove on opposite
                                    // orientation.
                                    if vi_ori != *ve_ori {
                                        self.my_vertices.remove(it);
                                        let islook = Self::fun_islook(brep, &self.my_edge);
                                        if (ve_bound && (vsd_same_ve || vi_same_ve)) && islook {
                                            add = true; // ofv
                                        }
                                    }
                                }
                            }
                        }
                        // L262-266 (ofv): addVE = add.
                        add_ve = add;
                        break;
                    }
                    it += 1;
                }

                // L269-276: if VE not found in the list, add it.
                if add_ve {
                    // L271: parVE = BRep_Tool::Parameter(VE, myEdge).
                    let par_ve = brep_tool_parameter(brep, ve, &self.my_edge);
                    // L272-273: myVertices.Prepend(new TopOpeBRepBuild_Pave(VE,
                    // parVE, true)).
                    let new_pv = Pave::new(ve, par_ve, true);
                    self.my_vertices.insert(0, new_pv);
                    edge_vertex_count += 1;
                }
            }
        } // myRemovePV

        // L280-294.
        let ll = self.my_vertices.len();
        if ll == edge_vertex_count {
            // if no more interferences vertices, clear the list
            self.my_vertices.clear();
        } else if ll >= 2 {
            // sort the parametrized vertices on Parameter() value.
            let list = std::mem::take(&mut self.my_vertices);
            Self::sort_pave(list, &mut self.my_vertices);
        }

        self.my_prepare_done = true;
    }

    /// PaveSet.cxx L323-330: InitLoop — Prepare runs here (the OCCT
    /// `if (!myPrepareDone) Prepare()` gate).
    pub(crate) fn init_loop(&mut self, brep: &BRep) {
        if !self.my_prepare_done {
            self.prepare(brep);
        }
        self.my_vertices_it = 0;
    }

    /// PaveSet.cxx L334-338: MoreLoop.
    pub(crate) fn more_loop(&self) -> bool {
        self.my_vertices_it < self.my_vertices.len()
    }

    /// PaveSet.cxx L342-345: NextLoop.
    pub(crate) fn next_loop(&mut self) {
        self.my_vertices_it += 1;
    }

    /// PaveSet.cxx L349-352: Loop.
    pub(crate) fn loop_(&self) -> Pave {
        self.my_vertices[self.my_vertices_it].clone()
    }

    /// PaveSet.cxx L363-443: HasEqualParameters.
    pub(crate) fn has_equal_parameters(&mut self) -> bool {
        self.my_has_equal_parameters = false;

        for it1 in 0..self.my_vertices.len() {
            if self.my_has_equal_parameters {
                break;
            }
            let v1 = self.my_vertices[it1].vertex().clone();
            let p1 = self.my_vertices[it1].parameter();
            for it2 in 0..self.my_vertices.len() {
                if self.my_has_equal_parameters {
                    break;
                }
                let v2 = self.my_vertices[it2].vertex().clone();
                if shapes_equal(&v2, &v1) {
                    continue;
                }
                let p2 = self.my_vertices[it2].parameter();
                let d = (p1 - p2).abs();
                if d < PCONFUSION {
                    self.my_has_equal_parameters = true;
                    self.my_equal_parameters = p1;
                }
            }
        }

        if !self.my_has_equal_parameters {
            // L400-415: rd = the edge carries a 3d curve; f = the edge
            // range first (the BRep_Tool::Curve first output = the edge
            // range).
            let ed = self.my_edge.as_edge().expect("PaveSet edge");
            let rd = ed.curve.is_some();
            let f = if rd { ed.range[0] } else { 0.0 };
            if rd {
                for it1 in 0..self.my_vertices.len() {
                    if self.my_has_equal_parameters {
                        break;
                    }
                    let p1 = self.my_vertices[it1].parameter();
                    let d = (p1 - f).abs();
                    if d < PCONFUSION {
                        self.my_has_equal_parameters = true;
                        self.my_equal_parameters = f;
                    }
                }
            }
        }

        self.my_has_equal_parameters
    }

    /// PaveSet.cxx L447-458: EqualParameters.
    pub(crate) fn equal_parameters(&self) -> f64 {
        if self.my_has_equal_parameters {
            return self.my_equal_parameters;
        }
        0.0 // windowsNT
    }

    /// PaveSet.cxx L462-490: ClosedVertices.
    pub(crate) fn closed_vertices(&mut self) -> bool {
        if self.my_vertices.is_empty() {
            return false;
        }

        let mut vmin: Option<Shape> = None;
        let mut vmax: Option<Shape> = None;
        let mut parmin = f64::INFINITY;
        let mut parmax = f64::NEG_INFINITY;
        for pv in &self.my_vertices {
            let v = pv.vertex().clone();
            let par = pv.parameter();
            if par > parmax {
                vmax = Some(v.clone());
                parmax = par;
            }
            if par < parmin {
                vmin = Some(v.clone());
                parmin = par;
            }
        }

        // L488: myClosed = Vmin.IsSame(Vmax).
        self.my_closed = match (&vmin, &vmax) {
            (Some(a), Some(b)) => shapes_same(a, b),
            _ => false,
        };
        self.my_closed
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_PaveClassifier (PaveClassifier.cxx L35-393).
// =========================================================================
pub(crate) struct PaveClassifier {
    /// OCCT: TopoDS_Edge myEdge.
    #[allow(dead_code)]
    my_edge: Shape,
    /// OCCT: bool myEdgePeriodic.
    my_edge_periodic: bool,
    /// OCCT: bool mySameParameters.
    my_same_parameters: bool,
    /// OCCT: bool myClosedVertices.
    my_closed_vertices: bool,
    /// OCCT: double myFirst.
    my_first: f64,
    /// OCCT: double myPeriod.
    my_period: f64,
    /// OCCT: double myP1 / myP2.
    my_p1: f64,
    my_p2: f64,
    /// OCCT: TopAbs_Orientation myO1 / myO2.
    my_o1: Orientation,
    my_o2: Orientation,
    /// OCCT: int myCas1 / myCas2.
    my_cas1: i32,
    my_cas2: i32,
}

impl PaveClassifier {
    /// PaveClassifier.cxx L35-98.
    pub(crate) fn new(brep: &BRep, e: &Shape) -> Self {
        let mut c = PaveClassifier {
            my_edge: e.clone(),
            my_edge_periodic: false,
            my_same_parameters: false,
            my_closed_vertices: false,
            my_first: 0.0,
            my_period: 0.0,
            my_p1: 0.0,
            my_p2: 0.0,
            my_o1: Orientation::Forward,
            my_o2: Orientation::Forward,
            my_cas1: 0,
            my_cas2: 0,
        };

        let ed = e.as_edge().expect("PaveClassifier edge");
        if !ed.degenerated {
            // L44-47: C = BRep_Tool::Curve(myEdge, loc, f, l) — the edge
            // curve + the edge range (f, l).
            let (f, l) = (ed.range[0], ed.range[1]);
            if let Some(curve) = &ed.curve {
                if curve.is_periodic() {
                    // L50-52: TopExp::Vertices(myEdge, v1, v2) — v1 FORWARD,
                    // v2 REVERSED.
                    let (v1, v2) = top_exp_vertices(brep, e);
                    if !v1.is_null() && !v2.is_null() {
                        // L55-64: the edge has vertices.
                        c.my_first = f;
                        let [fc, lc] = curve.default_domain();
                        c.my_period = lc - fc;
                        c.my_edge_periodic = shapes_same(&v1, &v2);
                        c.my_same_parameters = c.my_edge_periodic;
                        if c.my_same_parameters {
                            c.my_first = brep_tool_parameter(brep, &v1, e);
                        }
                    } else {
                        // L66-73: the edge has no vertices.
                        c.my_first = f;
                        c.my_period = l - f;
                        c.my_edge_periodic = true;
                        c.my_same_parameters = false;
                    }
                }
            }
        } // ! degenerated

        c
    }

    /// PaveClassifier.cxx L102-177: CompareOnNonPeriodic.
    fn compare_on_non_periodic(&mut self) -> TopAbsState {
        let mut state = TopAbsState::Unknown;
        let lower: bool;
        match self.my_o2 {
            Orientation::Forward => {
                lower = false;
            }
            Orientation::Reversed => {
                lower = true;
            }
            Orientation::Internal => {
                lower = false;
                state = TopAbsState::In;
            }
            Orientation::External => {
                lower = false;
                state = TopAbsState::Out;
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
                if lower {
                    state = TopAbsState::In;
                } else {
                    state = TopAbsState::Out;
                }
            } else if lower {
                state = TopAbsState::Out;
            } else {
                state = TopAbsState::In;
            }
        }

        state
    }

    /// PaveClassifier.cxx L181-217: AdjustCase.
    fn adjust_case(
        &self,
        p1: f64,
        o: Orientation,
        first: f64,
        period: f64,
        tol: f64,
        cas: &mut i32,
    ) -> f64 {
        let p2;
        if (p1 - first).abs() < tol {
            // p1 is first
            if o == Orientation::Reversed {
                p2 = p1 + period;
                *cas = 1;
            } else {
                p2 = p1;
                *cas = 2;
            }
        } else {
            // p1 is not on first
            let last = first + period;
            if (p1 - last).abs() < tol {
                // p1 is on last
                p2 = p1;
                *cas = 3;
            } else {
                // p1 is not on last
                p2 = in_period(p1, first, last);
                *cas = 4;
            }
        }
        p2
    }

    /// PaveClassifier.cxx L221-262: AdjustOnPeriodic.
    fn adjust_on_periodic(&mut self) {
        if !self.to_adjust_on_periodic() {
            return;
        }

        let tol = PCONFUSION;

        if self.my_same_parameters {
            let mut cas1 = 0;
            self.my_p1 = self.adjust_case(self.my_p1, self.my_o1, self.my_first, self.my_period, tol, &mut cas1);
            self.my_cas1 = cas1;
            let mut cas2 = 0;
            self.my_p2 = self.adjust_case(self.my_p2, self.my_o2, self.my_first, self.my_period, tol, &mut cas2);
            self.my_cas2 = cas2;
        } else if self.my_o1 != self.my_o2 {
            if self.my_o1 == Orientation::Forward {
                let mut cas2 = 0;
                self.my_p2 = self.adjust_case(self.my_p2, self.my_o2, self.my_p1, self.my_period, tol, &mut cas2);
                self.my_cas2 = cas2;
            }
            if self.my_o2 == Orientation::Forward {
                let mut cas1 = 0;
                self.my_p1 = self.adjust_case(self.my_p1, self.my_o1, self.my_p2, self.my_period, tol, &mut cas1);
                self.my_cas1 = cas1;
            }
        }
    }

    /// PaveClassifier.cxx L266-270: ToAdjustOnPeriodic.
    fn to_adjust_on_periodic(&self) -> bool {
        self.my_same_parameters || (self.my_o1 != self.my_o2)
    }

    /// PaveClassifier.cxx L274-309: CompareOnPeriodic.
    fn compare_on_periodic(&mut self) -> TopAbsState {
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

    /// PaveClassifier.cxx L313-362: Compare.
    pub(crate) fn compare(&mut self, l1: &Pave, l2: &Pave) -> TopAbsState {
        self.my_cas1 = 0; // debug
        self.my_cas2 = 0; // debug
        self.my_o1 = l1.vertex().orientation;
        self.my_o2 = l2.vertex().orientation;
        self.my_p1 = l1.parameter();
        self.my_p2 = l2.parameter();

        if self.my_edge_periodic && self.to_adjust_on_periodic() {
            self.adjust_on_periodic();
        }

        let state = if self.my_edge_periodic {
            self.compare_on_periodic()
        } else {
            self.compare_on_non_periodic()
        };

        state
    }

    /// PaveClassifier.cxx L366-375: SetFirstParameter.
    pub(crate) fn set_first_parameter(&mut self, p: f64) {
        self.my_first = p;
        self.my_same_parameters = true;
    }

    /// PaveClassifier.cxx L379-392: ClosedVertices(Closed).
    #[allow(dead_code)]
    pub(crate) fn closed_vertices(&mut self, closed: bool) {
        self.my_closed_vertices = closed;
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Area1dBuilder (Area1dBuilder.cxx L67-298) + the
// AreaBuilder base iteration/helpers (AreaBuilder.cxx L334-429, L59-121).
// The 1d area walk over the pave set: each interference pave (a block
// loop) groups with the bound pave (a shape loop) that bounds the piece
// valid for the state — the paves left in boundaryloops are never built.
// =========================================================================
pub(crate) struct PaveAreaBuilder {
    /// OCCT: myUNKNOWNRaise = false (ctor AreaBuilder.cxx L30-33: no raise
    /// if UNKNOWN state found).
    my_unknown_raise: bool,
    /// OCCT: myArea (list of list of pave loops).
    pub(crate) my_area: Vec<Vec<Pave>>,
    my_area_iterator: usize,
    my_loop_iterator: usize,
}

impl PaveAreaBuilder {
    /// AreaBuilder.cxx L30-33.
    pub(crate) fn new() -> Self {
        PaveAreaBuilder {
            my_unknown_raise: false,
            my_area: Vec::new(),
            my_area_iterator: 0,
            my_loop_iterator: 0,
        }
    }

    /// AreaBuilder.cxx L59-103: CompareLoopWithListOfLoop.
    fn compare_loop_with_list_of_loop(
        &self,
        lc: &mut PaveClassifier,
        l: &Pave,
        lol: &[Pave],
        what: LoopEnum,
    ) -> TopAbsState {
        let mut state = TopAbsState::Unknown;
        let mut totest: bool; // L must or not be tested

        if lol.is_empty() {
            return TopAbsState::Out;
        }

        for cur_l in lol {
            match what {
                LoopEnum::AnyLoop => totest = true,
                LoopEnum::Boundary => totest = cur_l.is_shape(),
                LoopEnum::Block => totest = !cur_l.is_shape(),
            }
            if totest {
                state = lc.compare(l, cur_l);
                if state == TopAbsState::Out {
                    // <L> is out of at least one Loop of <LOL> : stop to
                    // explore
                    break;
                }
            }
        }

        state
    }

    /// AreaBuilder.cxx L111-121: Atomize.
    fn atomize(&self, state: &mut TopAbsState, newstate: TopAbsState) {
        if self.my_unknown_raise {
            assert!(
                state != &TopAbsState::Unknown,
                "AreaBuilder : Position Unknown"
            );
        } else {
            *state = newstate;
        }
    }

    /// AreaBuilder.cxx L401-407: ADD_Loop_TO_LISTOFLoop.
    fn add_loop_to_list_of_loop(l: &Pave, lol: &mut Vec<Pave>) {
        lol.push(l.clone());
    }

    /// AreaBuilder.cxx L411-417: REM_Loop_FROM_LISTOFLoop.
    fn rem_loop_from_list_of_loop(ita: usize, a: &mut Vec<Pave>) {
        a.remove(ita);
    }

    /// Area1dBuilder.cxx L67-298: InitAreaBuilder.
    pub(crate) fn init_area_builder(
        &mut self,
        brep: &BRep,
        ls: &mut PaveSet,
        lc: &mut PaveClassifier,
        force_class: bool,
    ) {
        let mut state: TopAbsState;
        let mut loopinside: bool;
        let mut loopoutside: bool;

        // boundaryloops : list of boundary loops out of the areas.
        let mut boundaryloops: Vec<Pave> = Vec::new();

        self.my_area.clear(); // Clear the list of Area to be built

        ls.init_loop(brep);
        while ls.more_loop() {
            // process a new loop : L is the new current Loop
            let l = ls.loop_();
            let boundary_l = l.is_shape();

            // L = shape et ForceClass  : on traite L comme un block
            // L = shape et !ForceClass : on traite L comme un pur shape
            // L = !shape               : on traite L comme un block
            let traitercommeblock = !boundary_l || force_class;
            if !traitercommeblock {
                // the loop L is a boundary loop :
                // - try to insert it in an existing area, such as L is
                //   inside all the block loops. Only block loops of the
                //   area are compared.
                // - if L could not be inserted, store it in list of
                //   boundary loops.
                loopinside = false;
                let mut area_iter = 0usize;
                while area_iter < self.my_area.len() {
                    let a_area = self.my_area[area_iter].clone();
                    if a_area.is_empty() {
                        area_iter += 1;
                        continue;
                    }
                    state = self.compare_loop_with_list_of_loop(lc, &l, &a_area, LoopEnum::Block);
                    if state == TopAbsState::Unknown {
                        self.atomize(&mut state, TopAbsState::In);
                    }
                    loopinside = state == TopAbsState::In;
                    if loopinside {
                        break;
                    }
                    area_iter += 1;
                } // end of Area scan

                if loopinside {
                    // ADD_Loop_TO_LISTOFLoop(L, aArea, "IN, to current area")
                    self.my_area[area_iter].push(l.clone());
                } else if !loopinside {
                    // ADD_Loop_TO_LISTOFLoop(L, boundaryloops, "! IN, to
                    // boundaryloops")
                    boundaryloops.push(l.clone());
                }
            } // end of boundary loop
            else {
                // the loop L is a block loop
                // if L is IN theArea :
                //   - stop area scan, insert L in theArea.
                //   - remove from the area all the loops outside L
                //   - make a new area with them, unless they are all boundary
                //   - if they are all boundary put them back in boundaryLoops
                // else :
                //   - create a new area with L.
                //   - insert boundary loops that are IN the new area
                //     (and remove them from 'boundaryloops')
                loopinside = false;
                let mut area_iter = 0usize;
                while area_iter < self.my_area.len() {
                    let a_area = self.my_area[area_iter].clone();
                    if a_area.is_empty() {
                        area_iter += 1;
                        continue;
                    }
                    state = self.compare_loop_with_list_of_loop(lc, &l, &a_area, LoopEnum::AnyLoop);
                    if state == TopAbsState::Unknown {
                        self.atomize(&mut state, TopAbsState::In);
                    }
                    loopinside = state == TopAbsState::In;
                    if loopinside {
                        break;
                    }
                    area_iter += 1;
                } // end of Area scan

                if loopinside {
                    let mut all_shape = true;
                    let mut removed_loops: Vec<Pave> = Vec::new();
                    let mut loop_iter = 0usize;
                    while loop_iter < self.my_area[area_iter].len() {
                        let cur_l = self.my_area[area_iter][loop_iter].clone();
                        state = lc.compare(&cur_l, &l);
                        if state == TopAbsState::Unknown {
                            self.atomize(&mut state, TopAbsState::In); // not OUT
                        }
                        loopoutside = state == TopAbsState::Out;
                        if loopoutside {
                            // remove the loop from the area
                            // ADD_Loop_TO_LISTOFLoop(curL, removedLoops,
                            // "loopoutside = 1, area = removedLoops")
                            removed_loops.push(cur_l.clone());
                            all_shape = all_shape && cur_l.is_shape();
                            // REM_Loop_FROM_LISTOFLoop(LoopIter,
                            // AreaIter.ChangeValue(), "loop of cur. area,
                            // cur. area")
                            Self::rem_loop_from_list_of_loop(loop_iter, &mut self.my_area[area_iter]);
                        } else {
                            loop_iter += 1;
                        }
                    }
                    // insert the loop in the area
                    // ADD_Loop_TO_LISTOFLoop(L, aArea, "area = current")
                    self.my_area[area_iter].push(l.clone());
                    if !removed_loops.is_empty() {
                        if all_shape {
                            // ADD_LISTOFLoop_TO_LISTOFLoop(removedLoops,
                            // boundaryloops, "allShape = 1")
                            boundaryloops.extend(removed_loops);
                        } else {
                            // make a new area with the removed loops
                            // ADD_LISTOFLoop_TO_LISTOFLoop(removedLoops,
                            // myArea.Last(), "allShape = 0")
                            self.my_area.push(removed_loops);
                        }
                    }
                } // Loopinside == True
                else {
                    let mut ashapeinside: bool;
                    let mut ablockinside: bool;
                    self.my_area.push(Vec::new());
                    let new_area0_idx = self.my_area.len() - 1;
                    // ADD_Loop_TO_LISTOFLoop(L, newArea0, "new area")
                    self.my_area[new_area0_idx].push(l.clone());

                    let mut loop_iter = 0usize;
                    while loop_iter < boundaryloops.len() {
                        ashapeinside = false;
                        ablockinside = false;
                        let lb = boundaryloops[loop_iter].clone();
                        state = lc.compare(&lb, &l);
                        if state == TopAbsState::Unknown {
                            self.atomize(&mut state, TopAbsState::In);
                        }
                        ashapeinside = state == TopAbsState::In;
                        if ashapeinside {
                            state = lc.compare(&l, &lb);
                            if state == TopAbsState::Unknown {
                                self.atomize(&mut state, TopAbsState::In);
                            }
                            ablockinside = state == TopAbsState::In;
                        }
                        if ashapeinside && ablockinside {
                            // ADD_Loop_TO_LISTOFLoop(curL, newArea0,
                            // "ashapeinside && ablockinside, new area")
                            self.my_area[new_area0_idx].push(lb.clone());
                            // REM_Loop_FROM_LISTOFLoop(LoopIter,
                            // boundaryloops, "loop of boundaryloops,
                            // boundaryloops")
                            Self::rem_loop_from_list_of_loop(loop_iter, &mut boundaryloops);
                        } else {
                            loop_iter += 1;
                        }
                    } // end of boundaryloops scan
                } // Loopinside == False
            } // end of block loop
            // the OCCT for-loop third clause: LS.NextLoop().
            ls.next_loop();
        } // end of LoopSet LS scan

        // L297: InitArea().
        self.init_area();
    }

    /// AreaBuilder.cxx L334-340: InitArea.
    fn init_area(&mut self) {
        self.my_area_iterator = 0;
        self.init_loop();
    }

    /// AreaBuilder.cxx L344-348: MoreArea.
    pub(crate) fn more_area(&self) -> bool {
        self.my_area_iterator < self.my_area.len()
    }

    /// AreaBuilder.cxx L352-356: NextArea.
    pub(crate) fn next_area(&mut self) {
        self.my_area_iterator += 1;
        self.init_loop();
    }

    /// AreaBuilder.cxx L360-374: InitLoop.
    fn init_loop(&mut self) {
        if self.my_area_iterator < self.my_area.len() {
            self.my_loop_iterator = 0;
        } else {
            // Create an empty ListIteratorOfListOfLoop
            self.my_loop_iterator = 0;
        }
    }

    /// AreaBuilder.cxx L378-382: MoreLoop.
    pub(crate) fn more_loop(&self) -> bool {
        self.my_area_iterator < self.my_area.len()
            && self.my_loop_iterator < self.my_area[self.my_area_iterator].len()
    }

    /// AreaBuilder.cxx L386-389: NextLoop.
    pub(crate) fn next_loop(&mut self) {
        self.my_loop_iterator += 1;
    }

    /// AreaBuilder.cxx L393-397: Loop.
    pub(crate) fn loop_(&self) -> &Pave {
        &self.my_area[self.my_area_iterator][self.my_loop_iterator]
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_EdgeBuilder (EdgeBuilder.cxx L23-105) — the thin
// Area1dBuilder alias the MakeEdges walk consumes.
// =========================================================================
pub(crate) struct EdgeBuilder {
    pub(crate) area_builder: PaveAreaBuilder,
}

impl EdgeBuilder {
    /// EdgeBuilder.cxx L27-40: ctor -> InitEdgeBuilder -> InitAreaBuilder
    /// (ForceClass = false at the SplitEdge1 call site, Builder.cxx L1044).
    pub(crate) fn new(
        brep: &BRep,
        ls: &mut PaveSet,
        lc: &mut PaveClassifier,
        force_class: bool,
    ) -> Self {
        let mut ab = PaveAreaBuilder::new();
        ab.init_area_builder(brep, ls, lc, force_class);
        EdgeBuilder { area_builder: ab }
    }

    /// EdgeBuilder.cxx L45-48: InitEdge -> InitArea.
    pub(crate) fn init_edge(&mut self) {
        self.area_builder.init_area();
    }

    /// EdgeBuilder.cxx L52-56: MoreEdge -> MoreArea.
    pub(crate) fn more_edge(&self) -> bool {
        self.area_builder.more_area()
    }

    /// EdgeBuilder.cxx L60-63: NextEdge -> NextArea.
    pub(crate) fn next_edge(&mut self) {
        self.area_builder.next_area();
    }

    /// EdgeBuilder.cxx L67-70: InitVertex -> InitLoop.
    pub(crate) fn init_vertex(&mut self) {
        self.area_builder.init_loop();
    }

    /// EdgeBuilder.cxx L74-78: MoreVertex -> MoreLoop.
    pub(crate) fn more_vertex(&self) -> bool {
        self.area_builder.more_loop()
    }

    /// EdgeBuilder.cxx L82-85: NextVertex -> NextLoop.
    pub(crate) fn next_vertex(&mut self) {
        self.area_builder.next_loop();
    }

    /// EdgeBuilder.cxx L89-95: Vertex.
    pub(crate) fn vertex(&self) -> &Shape {
        self.area_builder.loop_().vertex()
    }

    /// EdgeBuilder.cxx L99-105: Parameter.
    pub(crate) fn parameter(&self) -> f64 {
        self.area_builder.loop_().parameter()
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Builder::FillVertexSet / FillVertexSetOnValue
// (Builder.cxx L1991-2055) — the pave walk over the edge point
// interferences: each pave carries the vertex with the transition
// orientation of the state being built.
// =========================================================================
impl TopOpeBRepBuildHBuilder {
    /// Builder.cxx L1991-1999 — the OCCT TopOpeBRepDS_PointIterator walk
    /// (the Match() filter keeps the POINT / VERTEX interferences,
    /// PointIterator.cxx L28-32).
    fn fill_vertex_set(
        &self,
        interfs: &[TopOpeBRepDSInterference],
        to_build: TopAbsState,
        pvs: &mut PaveSet,
    ) {
        for i in interfs {
            let (kind_g, index_g) = match i {
                TopOpeBRepDSInterference::CurvePoint(ci) => (ci.kind_g, ci.index_g),
                _ => continue,
            };
            self.fill_vertex_set_on_value(i, kind_g, index_g, to_build, pvs);
        }
    }

    /// Builder.cxx L2003-2055.
    fn fill_vertex_set_on_value(
        &self,
        interf: &TopOpeBRepDSInterference,
        kind_g: TopOpeBRepDSKind,
        ind: i32,
        to_build: TopAbsState,
        pvs: &mut PaveSet,
    ) {
        let ds = self.my_data_structure.as_ref().expect("Perform first");

        // L2010-2011: ind = index of new point or existing vertex;
        // bool ispoint = IT.IsPoint().
        let ispoint = kind_g == TopOpeBRepDSKind::Point;
        // L2014-2022: if (ispoint && ind <= NbPoints()) V = NewVertex(ind);
        // else V = myDataStructure->Shape(ind).
        let v: Shape = if ispoint && ind <= ds.side.points.len() as i32 {
            // L2016: V = NewVertex(ind) — the BuildVertices product.
            match self.my_new_vertices.get(&ind) {
                Some(v) => v.clone(),
                None => Shape::null(),
            }
        } else {
            // L2021: V = myDataStructure->Shape(ind).
            ds.shape(ind).clone()
        };
        // L2023: par = IT.Parameter().
        let par = interf.parameter();
        // L2024: ori = IT.Orientation(ToBuild) — the transition orientation
        // of the state (TopOpeBRepDS_Transition::Orientation(S, T)).
        let ori = interf.transition().orientation_for_state(to_build);

        // L2026-2029: bool keep = true (the commented-out EXTERNAL /
        // INTERNAL skip stays out).
        let keep = true;

        if keep {
            // L2031: myBuildTool.Orientation(V, ori).
            let v = shape_oriented(&v, ori);
            // L2032: PVS.Append(new TopOpeBRepBuild_Pave(V, par, false)).
            pvs.append(Pave::new(&v, par, false));
        }
    }
}

// =========================================================================
// OCCT TopOpeBRepBuild_Builder::SplitEdge / SplitEdge1 (Builder.cxx
// L925-937 / L941-1079).
// =========================================================================
impl TopOpeBRepBuildHBuilder {
    /// Builder.cxx L925-937 — the GetcontextSF2 debug switch is compiled
    /// out, so SplitEdge1 is the body.
    pub(crate) fn split_edge(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        e: &Shape,
        to_build1: TopAbsState,
        to_build2: TopAbsState,
    ) {
        self.split_edge1(brep, ds, e, to_build1, to_build2);
    }

    /// Builder.cxx L941-1079.
    fn split_edge1(
        &mut self,
        brep: &mut BRep,
        ds: &TopOpeBRepDSHDataStructure,
        eoriented: &Shape,
        to_build1: TopAbsState,
        to_build2: TopAbsState,
    ) {
        // L947-948: work on a FORWARD edge <Eforward>.
        let eforward = shape_oriented(eoriented, Orientation::Forward);

        // L950: bool tosplit = ToSplit(Eoriented, ToBuild1).
        let tosplit = self.to_split(ds, eoriented, to_build1);

        // L963-966.
        if !tosplit {
            return;
        }

        // L968-969: Reverse(ToBuild1, ToBuild2); Reverse(ToBuild2, ToBuild1)
        // — the results are unused in SplitEdge1.
        let _ = Self::builder_reverse(to_build1, to_build2);
        let _ = Self::builder_reverse(to_build2, to_build1);
        let connect_to1 = true;
        let connect_to2 = false;

        // L974-976: build the list of edges to split : LE1, LE2.
        let mut le1: Vec<Shape> = Vec::new();
        let mut le2: Vec<Shape> = Vec::new();
        le1.push(eforward.clone());
        self.find_same_domain(ds, &mut le1, &mut le2);

        // L1010: Make a PaveSet <PVS> on edge <Eforward>.
        let mut pvs = PaveSet::new(&eforward);

        // L1013-1014: Add the points/vertices found on edge <Eforward> in
        // <PVS> — TopOpeBRepDS_PointIterator EPIT(myDataStructure->
        // EdgePoints(Eforward)); the facade interference list is the
        // EdgePoints payload (the PointIterator Match() filter runs in
        // fill_vertex_set).  Index space: the facade shape_interferences map
        // is keyed by the 1-based add_shape index = bopds index + 1 (the
        // merge_solid precedent).
        let ishape = ds.bopds.index(&eforward);
        let interfs: Vec<TopOpeBRepDSInterference> = if ishape >= 0 {
            ds.shape_interferences(ishape as i32 + 1).to_vec()
        } else {
            Vec::new()
        };
        self.fill_vertex_set(&interfs, to_build1, &mut pvs);

        // L1016-1021: PaveClassifier VCL(Eforward); equalpar.
        let mut vcl = PaveClassifier::new(brep, &eforward);
        let equalpar = pvs.has_equal_parameters();
        if equalpar {
            vcl.set_first_parameter(pvs.equal_parameters());
        }

        // L1027: before return if PVS has no vertices, mark <Eforward> as
        // split <ToBuild1>.
        self.mark_split(&eforward, to_build1, true);

        // L1029-1041: PVS.InitLoop(); if (!PVS.MoreLoop()) return.
        pvs.init_loop(brep);
        if !pvs.more_loop() {
            return;
        }

        // L1044: build the new edges — TopOpeBRepBuild_EdgeBuilder
        // EBU(PVS, VCL) (the 1d area walk over the pave set).
        let mut ebu = EdgeBuilder::new(brep, &mut pvs, &mut vcl, false);

        // L1048-1049: EdgeList = ChangeMerged(Eforward, ToBuild1);
        // MakeEdges(Eforward, EBU, EdgeList).  The rcad form computes the
        // list then assigns (the change_merged borrow cannot span the
        // call); the OCCT ChangeMerged binding runs first.
        self.change_merged(&eforward, to_build1);
        let edge_list = self.make_edges(brep, &eforward, &mut ebu);
        *self.change_merged(&eforward, to_build1) = edge_list.clone();

        // L1053-1064: connect new edges as edges built <ToBuild1> on LE1
        // edge.
        for ecur in &le1 {
            // L1058: MarkSplit(Ecur, ToBuild1).
            self.mark_split(ecur, to_build1, true);
            // L1059-1063: ChangeSplit(Ecur, ToBuild1); if (ConnectTo1) EL =
            // EdgeList.
            self.change_split(ecur, to_build1);
            if connect_to1 {
                *self.change_split(ecur, to_build1) = edge_list.clone();
            }
        }

        // L1066-1077: connect new edges as edges built <ToBuild2> on LE2
        // edges.
        for ecur in &le2 {
            // L1071: MarkSplit(Ecur, ToBuild2).
            self.mark_split(ecur, to_build2, true);
            // L1072-1076: ChangeSplit(Ecur, ToBuild2); if (ConnectTo2) EL =
            // EdgeList.
            self.change_split(ecur, to_build2);
            if connect_to2 {
                *self.change_split(ecur, to_build2) = edge_list.clone();
            }
        }
    }

    /// OCCT TopOpeBRepBuild_Builder::MakeEdges (Merge.cxx L516-622) — one
    /// new edge per EdgeBuilder area (the pave pair), the vertices passing
    /// the betonnage (Merge.cxx L558-599).
    ///
    /// `an_edge` is the FORWARD source edge (OCCT anEdge); the copy
    /// (BuildTool.cxx L293-316 CopyEdge: EmptyCopied + Range(f,l)) plus the
    /// AddEdgeVertex / Parameter writes (BuildTool.cxx L915-918 /
    /// L950-993) map to the add_tedge primitive + the in-place TEdge edits
    /// (the hbuilder.rs precedent).
    fn make_edges(&mut self, brep: &mut BRep, an_edge: &Shape, edbu: &mut EdgeBuilder) -> Vec<Shape> {
        let mut l: Vec<Shape> = Vec::new();

        // L527: for (EDBU.InitEdge(); EDBU.MoreEdge(); EDBU.NextEdge()).
        edbu.init_edge();
        while edbu.more_edge() {
            // L530-539: 1 vertex sur edge courante => suppression edge.
            let mut nloop = 0;
            edbu.init_vertex();
            while edbu.more_vertex() {
                nloop += 1;
                edbu.next_vertex();
            }
            if nloop <= 1 {
                edbu.next_edge();
                continue;
            }

            // L541: myBuildTool.CopyEdge(anEdge, newEdge) — the empty copy
            // carrying the source curve, tolerance, pcurve representations
            // and range (BuildTool.cxx L293-316).
            let (curve, src_tol, src_pcurves, src_degenerated) = {
                let src = an_edge.as_edge().expect("MakeEdges source edge");
                (
                    src.curve.clone(),
                    src.tolerance,
                    src.pcurves.clone(),
                    src.degenerated,
                )
            };

            // L543: bool hasvertex = false.
            let mut hasvertex = false;
            // The betonnage-added (vertex, parameter) pairs — the OCCT
            // TopExp_Explorer scan over the vertices already added to
            // newEdge (Merge.cxx L562-591).
            let mut added: Vec<(Shape, f64)> = Vec::new();

            // L544: for (EDBU.InitVertex(); EDBU.MoreVertex();
            // EDBU.NextVertex()).
            edbu.init_vertex();
            while edbu.more_vertex() {
                // L546-547: TopoDS_Shape V = EDBU.Vertex(); Vori =
                // V.Orientation().
                let v = edbu.vertex().clone();
                let v_ori = v.orientation;

                // L549: bool hassd = myDataStructure->HasSameDomain(V) —
                // GAP carrier: the ChFi3d D6 facade fills no SameDomain
                // table, so the OCCT false path is what returns (the
                // L552-556 SameDomain-reference arm is unreachable).
                let hassd = false;
                if hassd {
                    // L552-556: on prend le vertex reference de V.
                }

                // L558-559: if (oriV != TopAbs_EXTERNAL).
                if v_ori != Orientation::External {
                    // L561: double parV = EDBU.Parameter().
                    let par_v = edbu.parameter();
                    // L561-591: betonnage — equafound scan over the
                    // vertices already added to newEdge.
                    let mut equafound = false;
                    for (ve, _par_ve) in &added {
                        let ori_ve = ve.orientation;
                        // L567-570: V.IsEqual(VE).
                        if shapes_equal(&v, ve) {
                            equafound = true;
                            break;
                        } else if ori_ve == Orientation::Forward
                            || ori_ve == Orientation::Reversed
                        {
                            // L572-580.
                            if v_ori == ori_ve {
                                equafound = true;
                                break;
                            }
                        } else if ori_ve == Orientation::Internal
                            || ori_ve == Orientation::External
                        {
                            // L581-590.
                            let par_ve = brep_tool_parameter(brep, ve, an_edge);
                            if par_v == par_ve {
                                equafound = true;
                                break;
                            }
                        }
                    }
                    if !equafound {
                        // L592-598: hasvertex = true; AddEdgeVertex(newEdge,
                        // V); Parameter(newEdge, V, parV).
                        hasvertex = true;
                        added.push((v, par_v));
                    }
                }
                edbu.next_vertex();
            } // loop on vertices of new edge newEdge

            // L617-620: if (hasvertex) L.Append(newEdge).  The rcad piece
            // carries the source curve / tolerance / pcurves (the CopyEdge
            // payload) with the surviving pave vertices and parameters (the
            // effective trim — the BuildTool::Parameter writes).
            if hasvertex {
                if added.len() == 2 {
                    let (va, vb) = (added[0].0.clone(), added[1].0.clone());
                    let (pa, pb) = (added[0].1, added[1].1);
                    let mut piece = brep.add_tedge(curve.clone(), va, vb, [pa, pb]);
                    {
                        let pd = brep.edge_mut_inplace(piece.clone());
                        pd.tolerance = src_tol;
                        pd.degenerated = src_degenerated;
                        for (fk, (pcv, t1, t2)) in &src_pcurves {
                            pd.pcurves.insert(*fk, (pcv.clone(), *t1, *t2));
                        }
                        pd.vertex_params.insert(added[0].0.ptr_id(), pa);
                        pd.vertex_params.insert(added[1].0.ptr_id(), pb);
                    }
                    // The OCCT copy inherits the source orientation
                    // (EmptyCopied) — anEdge is FORWARD at this call site.
                    piece.orientation = an_edge.orientation;
                    l.push(piece);
                } else {
                    // Architecture difference: the OCCT one-surviving-vertex
                    // append (a one-vertex TopoDS_Edge, Merge.cxx L617) has
                    // no rcad TEdge form (two bound vertices required); the
                    // areas of the pave walk pair one interference pave
                    // with one bound pave, both non-EXTERNAL after
                    // PaveSet::Prepare, so the arm is unreachable on this
                    // carrier.
                }
            }
            edbu.next_edge();
        } // loop on EDBU edges

        l
    }
}
