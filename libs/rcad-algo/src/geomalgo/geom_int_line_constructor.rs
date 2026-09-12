// OCCT GeomInt_LineConstructor (TKGeomAlgo/GeomInt) — 1:1 Rust translation.
//
// OCCT sources:
//   GeomInt_LineConstructor.hxx L20-69  (class, members, statics)
//   GeomInt_LineConstructor.lxx L21-67  (ctor, Load, IsDone, NbParts, Part)
//   GeomInt_LineConstructor.cxx L44-983 (GeomInt_Vertex L51-78, Perform
//     L114-670, TreatCircle L674-733, the file statics AdjustPeriodic
//     L737-816, Parameters L820-862, GLinePoint L866-891, RejectMicroCircle
//     L895-915, RejectDuplicates L926-983)
//
// Architecture differences (Rust vs C++):
//   - `occ::handle<Adaptor3d_TopolTool> myDom1/myDom2` are stored as the
//     `GeomAdaptor_Surface`-equivalent adaptors they are built over.  The rcad
//     `TopolTool<'a, S, ST>` borrows its adaptor, so storing the tool itself
//     would make the struct self-referential; the tool is rebuilt on demand by
//     [`GeomIntLineConstructor::dom1`]/[`dom2`].  This is observationally
//     identical because `Adaptor3d_TopolTool` is a pure function of the
//     adaptor domain (`Initialize(S)`, Adaptor3d_TopolTool.cxx L57-209) and the
//     rcad tool is the 1:1 translation of both `Initialize(S)` and `Classify`.
//   - `occ::handle<GeomAdaptor_Surface> myHS1/myHS2` map to
//     `hlr::contap::surface_adaptor::GeomSurfaceAdapter` (the landed 1:1
//     `GeomAdaptor_Surface`), whose `SurfaceAdapter` trait is
//     `Adaptor3d_HSurfaceTool`.
//   - `GeomInt_LineTool::{NbVertex, Vertex, FirstParameter, LastParameter}`
//     (GeomInt_LineTool.cxx L281-410) dispatch on `L->ArcType()`; rcad's
//     `IntPatchLine` carries all line kinds in one struct, so the dispatch
//     selects the same accessor (see [`nb_vertex`], [`vertex`],
//     [`first_parameter`], [`last_parameter`]).

use glam::DVec3;

use rcad_kernel::geom::{Curve3, CurveEval};
use rcad_kernel::precision::{INFINITE_VALUE, PCONFUSION, REAL_LAST};
use rcad_kernel::topods::State;

use crate::geomalgo::int_patch::cycy_common::in_period;
use crate::geomalgo::int_patch::elclib::{
    circle_value, ellipse_value, hyperbola_value, line_value, parabola_value,
};
use crate::geomalgo::int_patch::{IntPatchIType, IntPatchLine, IntPatchVertex};
use crate::geomalgo::int_surf::Quadric;
use crate::hlr::contap::geom_tool::GeomTool;
use crate::hlr::contap::surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter};
use crate::topalgo::adaptor3d::topol_tool::TopolTool;

/// OCCT GeomInt_LineConstructor.cxx L44: `static const double TwoPI = M_PI + M_PI`.
const TWO_PI: f64 = std::f64::consts::PI + std::f64::consts::PI;

// ============================================================================
// OCCT GeomInt_LineTool (GeomInt_LineTool.cxx L281-410) — the static line
// accessors GeomInt_LineConstructor dispatches through.
// ============================================================================

/// OCCT GeomInt_LineTool::NbVertex (cxx L281-295).
pub fn nb_vertex(l: &IntPatchLine) -> usize {
    match l.line_type {
        IntPatchIType::Analytic => l.nb_vertex(),
        IntPatchIType::Restriction => l.nb_vertex(),
        IntPatchIType::Walking => l.nb_vertex(),
        _ => l.nb_vertex(),
    }
}

/// OCCT GeomInt_LineTool::Vertex (cxx L299-313) — 1-based index.
pub fn vertex(l: &IntPatchLine, i: usize) -> &IntPatchVertex {
    match l.line_type {
        IntPatchIType::Analytic => l.vertex(i),
        IntPatchIType::Restriction => l.vertex(i),
        IntPatchIType::Walking => l.vertex(i),
        _ => l.vertex(i),
    }
}

/// OCCT GeomInt_LineTool::FirstParameter (cxx L317-366).
pub fn first_parameter(l: &IntPatchLine) -> f64 {
    let typl = l.line_type;
    match typl {
        IntPatchIType::Analytic => {
            if l.has_first_point() {
                return l.first_point().param_on_line;
            }
            let firstp = l.t_range[0];
            // OCCT: ALine->FirstParameter(included) reflects whether the bound
            // is included in the arc; rcad's IntPatchLine stores the closed
            // range, i.e. included == true (IntPatch_ALine::FirstParameter
            // returns myFirst when myIsArc is false).
            let _included = true;
            if !_included {
                return firstp + standard_epsilon(firstp);
            }
            firstp
        }
        IntPatchIType::Restriction => {
            if l.has_first_point() {
                l.first_point().param_on_line
            } else {
                -INFINITE_VALUE
            }
        }
        IntPatchIType::Walking => {
            if l.has_first_point() {
                l.first_point().param_on_line
            } else {
                1.
            }
        }
        _ => {
            if l.has_first_point() {
                return l.first_point().param_on_line;
            }
            match typl {
                IntPatchIType::Line | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
                    -INFINITE_VALUE
                }
                _ => 0.0,
            }
        }
    }
}

/// OCCT GeomInt_LineTool::LastParameter (cxx L371-410).
pub fn last_parameter(l: &IntPatchLine) -> f64 {
    let typl = l.line_type;
    match typl {
        IntPatchIType::Analytic => {
            if l.has_last_point() {
                return l.last_point().param_on_line;
            }
            let lastp = l.t_range[1];
            let _included = true;
            if !_included {
                return lastp - standard_epsilon(lastp);
            }
            lastp
        }
        IntPatchIType::Restriction => {
            if l.has_last_point() {
                l.last_point().param_on_line
            } else {
                INFINITE_VALUE
            }
        }
        IntPatchIType::Walking => {
            if l.has_last_point() {
                l.last_point().param_on_line
            } else {
                l.nb_points() as f64
            }
        }
        _ => {
            if l.has_last_point() {
                return l.last_point().param_on_line;
            }
            match typl {
                IntPatchIType::Line | IntPatchIType::Parabola | IntPatchIType::Hyperbola => {
                    INFINITE_VALUE
                }
                _ => 0.0,
            }
        }
    }
}

// ============================================================================
// OCCT GeomInt_LineConstructor.cxx L51-78: class GeomInt_Vertex
// ============================================================================

/// OCCT GeomInt_Vertex (cxx L51-78) — an IntPatch_Point wrapper providing the
/// sort by parameter on line.
#[derive(Debug, Clone, Default)]
pub struct GeomIntVertex {
    my_vertex: Option<IntPatchVertex>,
}

impl GeomIntVertex {
    /// OCCT GeomInt_Vertex() = default (cxx L55).
    pub fn new() -> Self {
        GeomIntVertex { my_vertex: None }
    }

    /// OCCT GeomInt_Vertex::SetVertex(theOther) (cxx L57-62).
    pub fn set_vertex(&mut self, the_other: &IntPatchVertex) {
        self.my_vertex = Some(the_other.clone());
        let a_new_param = in_period(the_other.param_on_line, 0.0, TWO_PI);
        self.set_parameter(a_new_param);
    }

    /// OCCT GeomInt_Vertex::SetParameter(theParam) (cxx L65).
    pub fn set_parameter(&mut self, the_param: f64) {
        self.my_vertex
            .as_mut()
            .expect("GeomInt_Vertex::SetParameter on empty vertex")
            .param_on_line = the_param;
    }

    /// OCCT GeomInt_Vertex::Getvertex() (cxx L68).
    pub fn get_vertex(&self) -> &IntPatchVertex {
        self.my_vertex
            .as_ref()
            .expect("GeomInt_Vertex::Getvertex on empty vertex")
    }

    /// OCCT GeomInt_Vertex::operator< (cxx L71-74).
    pub fn less_than(&self, the_other: &GeomIntVertex) -> bool {
        self.get_vertex().param_on_line < the_other.get_vertex().param_on_line
    }
}

/// `std::sort` on `GeomInt_Vertex` (cxx L71-74 operator<); the C++ sort is not
/// stable, hence `sort_unstable_by`.
fn sort_vertices(v: &mut [GeomIntVertex]) {
    v.sort_unstable_by(|a, b| {
        a.get_vertex()
            .param_on_line
            .partial_cmp(&b.get_vertex().param_on_line)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
}

// ============================================================================
// OCCT GeomInt_LineConstructor.hxx L27-65 / .lxx L21-67
// ============================================================================

/// OCCT GeomInt_LineConstructor (hxx L27-65).
pub struct GeomIntLineConstructor {
    // OCCT hxx L59: done
    done: bool,
    // OCCT hxx L60: NCollection_Sequence<double> seqp
    seqp: Vec<f64>,
    // OCCT hxx L61-62: myDom1 / myDom2 (see the file-header note).
    my_dom1: Option<GeomSurfaceAdapter>,
    my_dom2: Option<GeomSurfaceAdapter>,
    // OCCT hxx L63-64: myHS1 / myHS2
    my_hs1: Option<GeomSurfaceAdapter>,
    my_hs2: Option<GeomSurfaceAdapter>,
}

impl Default for GeomIntLineConstructor {
    fn default() -> Self {
        Self::new()
    }
}

impl GeomIntLineConstructor {
    /// OCCT GeomInt_LineConstructor() (lxx L21-24): done = false.
    pub fn new() -> Self {
        GeomIntLineConstructor {
            done: false,
            seqp: Vec::new(),
            my_dom1: None,
            my_dom2: None,
            my_hs1: None,
            my_hs2: None,
        }
    }

    /// OCCT GeomInt_LineConstructor::Load(D1, D2, S1, S2) (lxx L28-37).
    pub fn load(
        &mut self,
        d1: &GeomSurfaceAdapter,
        d2: &GeomSurfaceAdapter,
        s1: &GeomSurfaceAdapter,
        s2: &GeomSurfaceAdapter,
    ) {
        self.my_dom1 = Some(d1.clone());
        self.my_dom2 = Some(d2.clone());
        self.my_hs1 = Some(s1.clone());
        self.my_hs2 = Some(s2.clone());
    }

    /// OCCT GeomInt_LineConstructor::IsDone (lxx L41-44).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT GeomInt_LineConstructor::NbParts (lxx L48-55).
    pub fn nb_parts(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: GeomInt_LineConstructor::NbParts");
        }
        self.seqp.len() / 2
    }

    /// OCCT GeomInt_LineConstructor::Part(I, WFirst, WLast) (lxx L59-67) —
    /// 1-based index.
    pub fn part(&self, i: usize) -> (f64, f64) {
        if !self.done {
            panic!("StdFail_NotDone: GeomInt_LineConstructor::Part");
        }
        (self.seqp[2 * i - 2], self.seqp[2 * i - 1])
    }

    /// OCCT myDom1->Classify(...) — the stored `Adaptor3d_TopolTool` rebuilt
    /// over myDom1 (see the file-header architecture note).
    fn dom1(&self) -> TopolTool<'_, GeomSurfaceAdapter, GeomTool> {
        TopolTool::<GeomSurfaceAdapter, GeomTool>::new(
            self.my_dom1.as_ref().expect("GeomInt_LineConstructor::Load"),
        )
    }

    /// OCCT myDom2->Classify(...).
    fn dom2(&self) -> TopolTool<'_, GeomSurfaceAdapter, GeomTool> {
        TopolTool::<GeomSurfaceAdapter, GeomTool>::new(
            self.my_dom2.as_ref().expect("GeomInt_LineConstructor::Load"),
        )
    }

    fn hs1(&self) -> &GeomSurfaceAdapter {
        self.my_hs1.as_ref().expect("GeomInt_LineConstructor::Load")
    }

    fn hs2(&self) -> &GeomSurfaceAdapter {
        self.my_hs2.as_ref().expect("GeomInt_LineConstructor::Load")
    }

    // ========================================================================
    // OCCT GeomInt_LineConstructor::Perform (cxx L114-670)
    // ========================================================================

    /// OCCT GeomInt_LineConstructor::Perform(L) (cxx L114-670).
    pub fn perform(&mut self, l: &IntPatchLine) {
        // OCCT L118: constexpr double Tol = Precision::PConfusion() * 35.0;
        const A_TOL: f64 = PCONFUSION * 35.0;

        // OCCT L120: const IntPatch_IType typl = L->ArcType();
        let typl = l.line_type;
        if typl == IntPatchIType::Analytic {
            // OCCT L123-150: IntPatch_ALine branch.
            self.seqp.clear();
            let nbvtx = nb_vertex(l);
            // for (i = 1; i < nbvtx; i++)
            let mut i = 1usize;
            while i < nbvtx {
                let firstp = vertex(l, i).param_on_line;
                let lastp = vertex(l, i + 1).param_on_line;
                if firstp != lastp {
                    let pmid = (firstp + lastp) * 0.5;
                    // ALine->Value(pmid) (IntPatch_ALine::Value over the
                    // carried 3D curve).
                    let p_mid = l.curve.point_at(pmid);
                    let (u1, v1, u2, v2) =
                        parameters_two(self.hs1(), self.hs2(), p_mid);
                    let (u1, v1, u2, v2) =
                        adjust_periodic_pair(self.hs1(), self.hs2(), u1, v1, u2, v2);
                    let in1 = self.dom1().classify(glam::DVec2::new(u1, v1), A_TOL, true);
                    if in1 != State::Out {
                        let in2 = self.dom2().classify(glam::DVec2::new(u2, v2), A_TOL, true);
                        if in2 != State::Out {
                            self.seqp.push(firstp);
                            self.seqp.push(lastp);
                        }
                    }
                }
                i += 1;
            }
            self.done = true;
            return;
        } else if typl == IntPatchIType::Walking {
            // OCCT L152-329: IntPatch_WLine branch.
            self.seqp.clear();
            let nbvtx = nb_vertex(l);
            let mut i = 1usize;
            while i < nbvtx {
                let firstp = vertex(l, i).param_on_line;
                let lastp = vertex(l, i + 1).param_on_line;
                if firstp != lastp {
                    if lastp != firstp + 1. {
                        let pmid = ((firstp + lastp) / 2.) as i32;
                        // OCCT: WLine->Point(pmid) — 1-based (rcad's
                        // IntPatchLine::point is 0-based).
                        let p = &l.wline_pnts[(pmid - 1) as usize];
                        let (u1, v1, u2, v2) = (p.u1, p.v1, p.u2, p.v2);
                        let (u1, v1, u2, v2) =
                            adjust_periodic_pair(self.hs1(), self.hs2(), u1, v1, u2, v2);
                        let in1 = self.dom1().classify(glam::DVec2::new(u1, v1), A_TOL, true);
                        if in1 != State::Out {
                            let in2 =
                                self.dom2().classify(glam::DVec2::new(u2, v2), A_TOL, true);
                            if in2 != State::Out {
                                self.seqp.push(firstp);
                                self.seqp.push(lastp);
                            }
                        }
                    } else {
                        // OCCT L183: WLine->GetCreatingWay() == IntPatch_WLImpPrm
                        if l.wl_type == crate::geomalgo::int_patch::WLineType::ImpPrm {
                            // OCCT L202-224 (the fix #29972).
                            let a_p_first = &l.wline_pnts[(firstp as i32 - 1) as usize];
                            let a_p_last = &l.wline_pnts[(lastp as i32 - 1) as usize];
                            let (mut u1, mut v1, mut u2, mut v2) =
                                (a_p_first.u1, a_p_first.v1, a_p_first.u2, a_p_first.v2);
                            let (u1b, v1b, u2b, v2b) =
                                adjust_periodic_pair(self.hs1(), self.hs2(), u1, v1, u2, v2);
                            u1 = u1b;
                            v1 = v1b;
                            u2 = u2b;
                            v2 = v2b;
                            let (a_u21, a_v21, a_u22, a_v22) =
                                (a_p_last.u1, a_p_last.v1, a_p_last.u2, a_p_last.v2);
                            let (a_u21, a_v21, a_u22, a_v22) = adjust_periodic_pair(
                                self.hs1(),
                                self.hs2(),
                                a_u21,
                                a_v21,
                                a_u22,
                                a_v22,
                            );
                            u1 = 0.5 * (u1 + a_u21);
                            v1 = 0.5 * (v1 + a_v21);
                            u2 = 0.5 * (u2 + a_u22);
                            v2 = 0.5 * (v2 + a_v22);

                            let in1 =
                                self.dom1().classify(glam::DVec2::new(u1, v1), A_TOL, true);
                            if in1 != State::Out {
                                let in2 =
                                    self.dom2().classify(glam::DVec2::new(u2, v2), A_TOL, true);
                                if in2 != State::Out {
                                    self.seqp.push(firstp);
                                    self.seqp.push(lastp);
                                }
                            }
                        } else {
                            // OCCT L228-251.
                            let p_first = &l.wline_pnts[(firstp as i32 - 1) as usize];
                            let (u1a, v1a, u2a, v2a) =
                                (p_first.u1, p_first.v1, p_first.u2, p_first.v2);
                            let (u1, v1, u2, v2) =
                                adjust_periodic_pair(self.hs1(), self.hs2(), u1a, v1a, u2a, v2a);
                            let mut in1 =
                                self.dom1().classify(glam::DVec2::new(u1, v1), A_TOL, true);
                            if in1 != State::Out {
                                let mut in2 =
                                    self.dom2().classify(glam::DVec2::new(u2, v2), A_TOL, true);
                                if in2 != State::Out {
                                    let p_last = &l.wline_pnts[(lastp as i32 - 1) as usize];
                                    let (u1, v1, u2, v2) = adjust_periodic_pair(
                                        self.hs1(),
                                        self.hs2(),
                                        p_last.u1,
                                        p_last.v1,
                                        p_last.u2,
                                        p_last.v2,
                                    );
                                    in1 = self.dom1().classify(
                                        glam::DVec2::new(u1, v1),
                                        A_TOL,
                                        true,
                                    );
                                    if in1 != State::Out {
                                        in2 = self.dom2().classify(
                                            glam::DVec2::new(u2, v2),
                                            A_TOL,
                                            true,
                                        );
                                        if in2 != State::Out {
                                            self.seqp.push(firstp);
                                            self.seqp.push(lastp);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                i += 1;
            }
            // OCCT L265-326: the "7 connected segments" rejection.
            let a_nb_parts = self.seqp.len() / 2;
            if a_nb_parts > 1 {
                let a_st1 = self.hs1().get_type();
                let a_st2 = self.hs2().get_type();
                let mut b_cond = false;
                if a_st1 == crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane {
                    if a_st2 == crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfExtrusion
                        || a_st2
                            == crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfRevolution
                    {
                        b_cond = !b_cond;
                    }
                } else if a_st2 == crate::geomalgo::int_patch::GeomAbsSurfaceType::Plane {
                    if a_st1 == crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfExtrusion
                        || a_st1
                            == crate::geomalgo::int_patch::GeomAbsSurfaceType::SurfaceOfRevolution
                    {
                        b_cond = !b_cond;
                    }
                }
                if b_cond {
                    // OCCT L293-324: NCollection_IndexedMap<int> aMap.
                    let mut a_map: Vec<i32> = Vec::new();
                    let mut a_seq_tmp: Vec<f64> = Vec::new();
                    let a_nb = self.seqp.len();
                    let mut i = 1usize;
                    while i <= a_nb {
                        let lastp = self.seqp[i - 1];
                        let an_index = lastp as i32;
                        if !a_map.contains(&an_index) {
                            a_map.push(an_index);
                            a_seq_tmp.push(lastp);
                        } else {
                            let a_nb_tmp = a_seq_tmp.len();
                            a_seq_tmp.remove(a_nb_tmp - 1);
                        }
                        i += 1;
                    }
                    self.seqp.clear();
                    let a_nb = a_seq_tmp.len() / 2;
                    let mut i = 1usize;
                    while i <= a_nb {
                        let jx = 2 * i;
                        let firstp = a_seq_tmp[jx - 2];
                        let lastp = a_seq_tmp[jx - 1];
                        self.seqp.push(firstp);
                        self.seqp.push(lastp);
                        i += 1;
                    }
                }
            }
            self.done = true;
            return;
        } else if typl != IntPatchIType::Restriction {
            // OCCT L332-386: the GLine (Line / Circle / Ellipse / Parabola /
            // Hyperbola) branch.
            self.seqp.clear();
            if typl == IntPatchIType::Circle || typl == IntPatchIType::Ellipse {
                self.treat_circle(l, A_TOL);
                self.done = true;
                return;
            }
            // OCCT L345-382.
            let mut intrvtested = false;
            let nbvtx = nb_vertex(l);
            let mut i = 1usize;
            while i < nbvtx {
                let firstp = vertex(l, i).param_on_line;
                let lastp = vertex(l, i + 1).param_on_line;
                if (firstp - lastp).abs() > PCONFUSION {
                    intrvtested = true;
                    let pmid = (firstp + lastp) * 0.5;
                    let p_mid = g_line_point(typl, l, pmid);
                    let (u1, v1, u2, v2) = parameters_two(self.hs1(), self.hs2(), p_mid);
                    let (u1, v1, u2, v2) =
                        adjust_periodic_pair(self.hs1(), self.hs2(), u1, v1, u2, v2);
                    let in1 = self.dom1().classify(glam::DVec2::new(u1, v1), A_TOL, true);
                    if in1 != State::Out {
                        let in2 = self.dom2().classify(glam::DVec2::new(u2, v2), A_TOL, true);
                        if in2 != State::Out {
                            self.seqp.push(firstp);
                            self.seqp.push(lastp);
                        }
                    }
                }
                i += 1;
            }
            if !intrvtested {
                self.seqp.push(first_parameter(l));
                self.seqp.push(last_parameter(l));
            }
            self.done = true;
            return;
        }

        // OCCT L388-669: the IntPatch_Restriction branch.
        self.done = false;
        self.seqp.clear();
        let mut nbvtx = nb_vertex(l);
        if nbvtx == 0 {
            self.seqp.push(first_parameter(l));
            self.seqp.push(last_parameter(l));
            self.done = true;
            return;
        }

        // OCCT L400-401: NCollection_Sequence<GeomInt_ParameterAndOrientation>
        // seqpss; or1/or2 = TopAbs_FORWARD.
        let mut seqpss: Vec<ParameterAndOrientation> = Vec::new();
        let mut or1 = OrientationKind::Forward;
        let mut or2 = OrientationKind::Forward;

        let mut i = 1usize;
        while i <= nbvtx {
            let thevtx = vertex(l, i);
            let prm = thevtx.param_on_line;
            if thevtx.on_dom_s1 {
                match thevtx.transition_line_arc1 {
                    crate::geomalgo::int_patch::transitions::TypeTrans::In => {
                        or1 = OrientationKind::Forward
                    }
                    crate::geomalgo::int_patch::transitions::TypeTrans::Out => {
                        or1 = OrientationKind::Reversed
                    }
                    crate::geomalgo::int_patch::transitions::TypeTrans::Touch => {
                        or1 = OrientationKind::Internal
                    }
                    crate::geomalgo::int_patch::transitions::TypeTrans::Undecided => {
                        or1 = OrientationKind::Internal
                    }
                }
            } else {
                or1 = OrientationKind::Internal;
            }

            if thevtx.on_dom_s2 {
                match thevtx.transition_line_arc2 {
                    crate::geomalgo::int_patch::transitions::TypeTrans::In => {
                        or2 = OrientationKind::Forward
                    }
                    crate::geomalgo::int_patch::transitions::TypeTrans::Out => {
                        or2 = OrientationKind::Reversed
                    }
                    crate::geomalgo::int_patch::transitions::TypeTrans::Touch => {
                        or2 = OrientationKind::Internal
                    }
                    crate::geomalgo::int_patch::transitions::TypeTrans::Undecided => {
                        or2 = OrientationKind::Internal
                    }
                }
            } else {
                or2 = OrientationKind::Internal;
            }
            // OCCT L453-505: insert / accumulate into seqpss.
            let nbinserted = seqpss.len();
            let mut inserted = false;
            let mut j = 1usize;
            while j <= nbinserted {
                if (prm - seqpss[j - 1].parameter()).abs() <= A_TOL {
                    // accumulate
                    let valj = &mut seqpss[j - 1];
                    if or1 != OrientationKind::Internal {
                        if valj.orientation1() != OrientationKind::Internal {
                            if or1 != valj.orientation1() {
                                valj.set_orientation1(OrientationKind::Internal);
                            }
                        } else {
                            valj.set_orientation1(or1);
                        }
                    }
                    if or2 != OrientationKind::Internal {
                        if valj.orientation2() != OrientationKind::Internal {
                            if or2 != valj.orientation2() {
                                valj.set_orientation2(OrientationKind::Internal);
                            }
                        } else {
                            valj.set_orientation2(or2);
                        }
                    }
                    inserted = true;
                    break;
                }
                if prm < seqpss[j - 1].parameter() - A_TOL {
                    // insert before position j
                    seqpss.insert(j - 1, ParameterAndOrientation::new(prm, or1, or2));
                    inserted = true;
                    break;
                }
                j += 1;
            }
            if !inserted {
                seqpss.push(ParameterAndOrientation::new(prm, or1, or2));
            }
            i += 1;
        }

        // OCCT L508-523: determine the state at the beginning of line.
        let mut trim = false;
        let mut dans_s1 = false;
        let mut dans_s2 = false;

        nbvtx = seqpss.len();
        let mut i = 1usize;
        while i <= nbvtx {
            or1 = seqpss[i - 1].orientation1();
            if or1 != OrientationKind::Internal {
                trim = true;
                dans_s1 = or1 != OrientationKind::Forward;
                break;
            }
            i += 1;
        }

        if i > nbvtx {
            // OCCT L525-543.
            let nb_v = nb_vertex(l);
            let mut i = 1usize;
            while i <= nb_v {
                if !vertex(l, i).on_dom_s1 {
                    let v = vertex(l, i);
                    let ppcc = glam::DVec2::new(v.u1, v.v1);
                    if self.dom1().classify(ppcc, A_TOL, true) == State::Out {
                        self.done = true;
                        return;
                    }
                    break;
                }
                i += 1;
            }
            dans_s1 = true; // Keep in doubt
        }
        // OCCT L545-554.
        let mut i = 1usize;
        while i <= nbvtx {
            or2 = seqpss[i - 1].orientation2();
            if or2 != OrientationKind::Internal {
                trim = true;
                dans_s2 = or2 != OrientationKind::Forward;
                break;
            }
            i += 1;
        }
        if i > nbvtx {
            // OCCT L556-573.
            let nb_v = nb_vertex(l);
            let mut i = 1usize;
            while i <= nb_v {
                if !vertex(l, i).on_dom_s2 {
                    let v = vertex(l, i);
                    if self.dom2().classify(glam::DVec2::new(v.u2, v.v2), A_TOL, true) == State::Out {
                        self.done = true;
                        return;
                    }
                    break;
                }
                i += 1;
            }
            dans_s2 = true; // Keep in doubt
        }

        if !trim {
            // OCCT L575-581: necessarily dansS1 == dansS2 == true.
            self.seqp.push(first_parameter(l));
            self.seqp.push(last_parameter(l));
            self.done = true;
            return;
        }

        // OCCT L583-656: seqpss is peeled to create valid ends.
        let thefirst = first_parameter(l);
        let thelast = last_parameter(l);
        let mut firstp = thefirst;

        let mut i = 1usize;
        while i <= nbvtx {
            or1 = seqpss[i - 1].orientation1();
            or2 = seqpss[i - 1].orientation2();
            if dans_s1 && dans_s2 {
                if or1 == OrientationKind::Reversed {
                    dans_s1 = false;
                }
                if or2 == OrientationKind::Reversed {
                    dans_s2 = false;
                }
                if !dans_s1 || !dans_s2 {
                    let lastp = seqpss[i - 1].parameter();
                    let stofirst = firstp.max(thefirst);
                    let stolast = lastp.min(thelast);
                    if stolast > stofirst {
                        self.seqp.push(stofirst);
                        self.seqp.push(stolast);
                    }
                    if lastp > thelast {
                        break;
                    }
                }
            } else {
                if dans_s1 {
                    if or1 == OrientationKind::Reversed {
                        dans_s1 = false;
                    }
                } else if or1 == OrientationKind::Forward {
                    dans_s1 = true;
                }
                if dans_s2 {
                    if or2 == OrientationKind::Reversed {
                        dans_s2 = false;
                    }
                } else if or2 == OrientationKind::Forward {
                    dans_s2 = true;
                }
                if dans_s1 && dans_s2 {
                    firstp = seqpss[i - 1].parameter();
                }
            }
            i += 1;
        }
        // OCCT L658-668: finally to add.
        if dans_s1 && dans_s2 {
            let lastp = thelast;
            firstp = firstp.max(thefirst);
            if lastp > firstp {
                self.seqp.push(firstp);
                self.seqp.push(lastp);
            }
        }
        self.done = true;
    }

    // ========================================================================
    // OCCT GeomInt_LineConstructor::TreatCircle (cxx L674-733)
    // ========================================================================

    /// OCCT GeomInt_LineConstructor::TreatCircle(theLine, theTol) (cxx
    /// L674-733).
    fn treat_circle(&mut self, the_line: &IntPatchLine, the_tol: f64) {
        let a_type = the_line.line_type;
        if reject_micro_circle(the_line, a_type, the_tol) {
            return;
        }
        // OCCT L684-691: the vertex array [1, aNbVtx + 1], sorted.
        let a_nb_vtx = the_line.nb_vertex();
        let mut a_vtx_arr: Vec<GeomIntVertex> = Vec::with_capacity(a_nb_vtx + 1);
        for _ in 0..=a_nb_vtx {
            a_vtx_arr.push(GeomIntVertex::new());
        }
        let mut i = 1usize;
        while i <= a_nb_vtx {
            let v = the_line.vertex(i).clone();
            a_vtx_arr[i - 1].set_vertex(&v);
            i += 1;
        }
        sort_vertices(&mut a_vtx_arr[0..a_nb_vtx]);

        // OCCT L694-695: create last vertex.
        let a_min_prm = a_vtx_arr[0].get_vertex().param_on_line + TWO_PI;
        a_vtx_arr[a_nb_vtx].set_parameter(a_min_prm);

        reject_duplicates(&mut a_vtx_arr);

        // std::sort(aVtxArr.begin(), aVtxArr.end()) — the whole array.
        sort_vertices(&mut a_vtx_arr);

        // OCCT L704-732.
        let mut i = 1usize; // aVtxArr.Lower() == 1 (0-based index 0)
        while i <= a_vtx_arr.len() - 1 {
            let a_t1 = a_vtx_arr[i - 1].get_vertex().param_on_line;
            let a_t2 = a_vtx_arr[i].get_vertex().param_on_line;

            if a_t2 == REAL_LAST {
                break;
            }

            let a_tmid = (a_t1 + a_t2) * 0.5;
            let a_pmid = g_line_point(a_type, the_line, a_tmid);
            let (u1, v1, u2, v2) = parameters_two(self.hs1(), self.hs2(), a_pmid);
            let (u1, v1, u2, v2) = adjust_periodic_pair(self.hs1(), self.hs2(), u1, v1, u2, v2);

            let mut a_state = self.dom1().classify(glam::DVec2::new(u1, v1), the_tol, true);
            if a_state != State::Out {
                a_state = self.dom2().classify(glam::DVec2::new(u2, v2), the_tol, true);
                if a_state != State::Out {
                    self.seqp.push(a_t1);
                    self.seqp.push(a_t2);
                }
            }
            i += 1;
        }
    }
}

// ============================================================================
// OCCT GeomInt_LineConstructor.cxx L400 / L36-53: GeomInt_ParameterAndOrientation
// + the TopAbs_Orientation values it stores.  OCCT uses `TopAbs_Orientation`
// directly; rcad's `GeomIntParameterAndOrientation` carries the same four
// values in the local [`OrientationKind`] enum (rcad's `topods::Orientation`
// has no INTERNAL value — the FINAL/EXTERNAL pair is split differently — so the
// four `TopAbs_Orientation` states are encoded here as in OCCT).
// ============================================================================

/// OCCT `TopAbs_FORWARD` / `TopAbs_REVERSED` / `TopAbs_INTERNAL`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrientationKind {
    Forward,
    Reversed,
    Internal,
}

/// OCCT GeomInt_ParameterAndOrientation (GeomInt_ParameterAndOrientation.hxx).
#[derive(Debug, Clone)]
pub struct ParameterAndOrientation {
    my_parameter: f64,
    my_or1: OrientationKind,
    my_or2: OrientationKind,
}

impl ParameterAndOrientation {
    /// OCCT GeomInt_ParameterAndOrientation(Par, Or1, Or2).
    pub fn new(par: f64, or1: OrientationKind, or2: OrientationKind) -> Self {
        ParameterAndOrientation {
            my_parameter: par,
            my_or1: or1,
            my_or2: or2,
        }
    }
    /// OCCT Parameter().
    pub fn parameter(&self) -> f64 {
        self.my_parameter
    }
    /// OCCT Orientation1().
    pub fn orientation1(&self) -> OrientationKind {
        self.my_or1
    }
    /// OCCT Orientation2().
    pub fn orientation2(&self) -> OrientationKind {
        self.my_or2
    }
    /// OCCT SetOrientation1().
    pub fn set_orientation1(&mut self, o: OrientationKind) {
        self.my_or1 = o;
    }
    /// OCCT SetOrientation2().
    pub fn set_orientation2(&mut self, o: OrientationKind) {
        self.my_or2 = o;
    }
}

// ============================================================================
// OCCT GeomInt_LineConstructor.cxx L737-816: file-static AdjustPeriodic
// ============================================================================

/// OCCT GeomInt::AdjustPeriodic (GeomInt.cxx L21-48).
pub fn geom_int_adjust_periodic(
    the_par: f64,
    the_par_min: f64,
    the_par_max: f64,
    the_period: f64,
    the_new_par: &mut f64,
    the_offset: &mut f64,
    the_eps: f64,
) -> bool {
    *the_offset = 0.;
    *the_new_par = the_par;
    let b_min = the_par_min - the_par > the_eps;
    let b_max = the_par - the_par_max > the_eps;
    if b_min || b_max {
        let dp = if b_min {
            the_par_max - the_par
        } else {
            the_par_min - the_par
        };
        // modf(dp / thePeriod, &aNbPer) — the integral part.
        let a_nb_per = (dp / the_period).trunc();
        *the_offset = a_nb_per * the_period;
        *the_new_par += *the_offset;
    }
    *the_offset > 0.
}

/// OCCT GeomInt_LineConstructor.cxx L737-816: the file-static AdjustPeriodic
/// over the two adaptor surfaces.
fn adjust_periodic_pair(
    my_hs1: &GeomSurfaceAdapter,
    my_hs2: &GeomSurfaceAdapter,
    u1: f64,
    v1: f64,
    u2: f64,
    v2: f64,
) -> (f64, f64, f64, f64) {
    use crate::geomalgo::int_patch::GeomAbsSurfaceType as ST;
    // OCCT L744-764.
    let typs1 = my_hs1.get_type();
    let (my_hs1_is_u_periodic, my_hs1_is_v_periodic) = match typs1 {
        ST::Cylinder | ST::Cone | ST::Sphere => (true, false),
        ST::Torus => (true, true),
        // Case of periodic biparameters is processed upstream.
        _ => (false, false),
    };
    // OCCT L765-785.
    let typs2 = my_hs2.get_type();
    let (my_hs2_is_u_periodic, my_hs2_is_v_periodic) = match typs2 {
        ST::Cylinder | ST::Cone | ST::Sphere => (true, false),
        ST::Torus => (true, true),
        _ => (false, false),
    };
    let mut du = 0.0f64;
    let mut dv = 0.0f64;
    // OCCT L786-815.  Note: the OCCT source comments out `myHS->UPeriod()` and
    // uses `M_PI + M_PI` literally.
    let mut u1 = u1;
    let mut v1 = v1;
    let mut u2 = u2;
    let mut v2 = v2;
    if my_hs1_is_u_periodic {
        let lmf = TWO_PI;
        let f = my_hs1.first_u_parameter();
        let l = my_hs1.last_u_parameter();
        let mut new_u = 0.0;
        geom_int_adjust_periodic(u1, f, l, lmf, &mut new_u, &mut du, PCONFUSION);
        u1 = new_u;
    }
    if my_hs1_is_v_periodic {
        let lmf = TWO_PI;
        let f = my_hs1.first_v_parameter();
        let l = my_hs1.last_v_parameter();
        let mut new_v = 0.0;
        geom_int_adjust_periodic(v1, f, l, lmf, &mut new_v, &mut dv, PCONFUSION);
        v1 = new_v;
    }
    if my_hs2_is_u_periodic {
        let lmf = TWO_PI;
        let f = my_hs2.first_u_parameter();
        let l = my_hs2.last_u_parameter();
        let mut new_u = 0.0;
        geom_int_adjust_periodic(u2, f, l, lmf, &mut new_u, &mut du, PCONFUSION);
        u2 = new_u;
    }
    if my_hs2_is_v_periodic {
        let lmf = TWO_PI;
        let f = my_hs2.first_v_parameter();
        let l = my_hs2.last_v_parameter();
        let mut new_v = 0.0;
        geom_int_adjust_periodic(v2, f, l, lmf, &mut new_v, &mut dv, PCONFUSION);
        v2 = new_v;
    }
    let _ = (du, dv);
    (u1, v1, u2, v2)
}

// ============================================================================
// OCCT GeomInt_LineConstructor.cxx L820-891: file-static Parameters / GLinePoint
// ============================================================================

/// OCCT `static void Parameters(myHS1, myHS2, Ptref, U1, V1, U2, V2)` (cxx
/// L820-830) + the single-surface overload (cxx L834-862, the per-surface
/// `IntSurf_Quadric` branch + Standard_ConstructionError fallthrough).
fn parameters_one(my_hs: &GeomSurfaceAdapter, ptref: DVec3) -> (f64, f64) {
    use crate::geomalgo::int_patch::GeomAbsSurfaceType as ST;
    let quad1 = match my_hs.get_type() {
        ST::Plane => Quadric::from_plane(&my_hs.plane()),
        ST::Cylinder => Quadric::from_cylinder(&my_hs.cylinder()),
        ST::Cone => Quadric::from_cone(&my_hs.cone()),
        ST::Sphere => Quadric::from_sphere(&my_hs.sphere()),
        ST::Torus => Quadric::from_torus(&my_hs.torus()),
        _ => panic!("Standard_ConstructionError: GeomInt_LineConstructor::Parameters"),
    };
    quad1.parameters(ptref)
}

/// OCCT `static void Parameters(myHS1, myHS2, Ptref, U1, V1, U2, V2)`
/// (cxx L820-830).
fn parameters_two(
    my_hs1: &GeomSurfaceAdapter,
    my_hs2: &GeomSurfaceAdapter,
    ptref: DVec3,
) -> (f64, f64, f64, f64) {
    let (u1, v1) = parameters_one(my_hs1, ptref);
    let (u2, v2) = parameters_one(my_hs2, ptref);
    (u1, v1, u2, v2)
}

/// OCCT `static void GLinePoint(typl, GLine, aT, aP)` (cxx L866-891).
fn g_line_point(typl: IntPatchIType, g_line: &IntPatchLine, a_t: f64) -> DVec3 {
    match typl {
        IntPatchIType::Line => match &g_line.curve {
            Curve3::Line(l) => line_value(l, a_t),
            _ => panic!("Standard_ConstructionError: GeomInt_LineConstructor::GLinePoint"),
        },
        IntPatchIType::Circle => match &g_line.curve {
            Curve3::Circle(c) => circle_value(c, a_t),
            _ => panic!("Standard_ConstructionError: GeomInt_LineConstructor::GLinePoint"),
        },
        IntPatchIType::Ellipse => match &g_line.curve {
            Curve3::Ellipse(e) => ellipse_value(e, a_t),
            _ => panic!("Standard_ConstructionError: GeomInt_LineConstructor::GLinePoint"),
        },
        IntPatchIType::Hyperbola => match &g_line.curve {
            Curve3::Hyperbola(h) => hyperbola_value(h, a_t),
            _ => panic!("Standard_ConstructionError: GeomInt_LineConstructor::GLinePoint"),
        },
        IntPatchIType::Parabola => match &g_line.curve {
            Curve3::Parabola(p) => parabola_value(p, a_t),
            _ => panic!("Standard_ConstructionError: GeomInt_LineConstructor::GLinePoint"),
        },
        _ => panic!("Standard_ConstructionError: GeomInt_LineConstructor::GLinePoint"),
    }
}

// ============================================================================
// OCCT GeomInt_LineConstructor.cxx L895-915: file-static RejectMicroCircle
// ============================================================================

/// OCCT `static bool RejectMicroCircle(aGLine, aType, aTol3D)` (cxx
/// L895-915).
fn reject_micro_circle(g_line: &IntPatchLine, a_type: IntPatchIType, a_tol3d: f64) -> bool {
    let mut b_ret = false;
    if a_type == IntPatchIType::Circle {
        let a_r = match &g_line.curve {
            Curve3::Circle(c) => c.radius,
            _ => 0.0,
        };
        b_ret = a_r < a_tol3d;
    } else if a_type == IntPatchIType::Ellipse {
        let a_r = match &g_line.curve {
            Curve3::Ellipse(e) => e.major_radius,
            _ => 0.0,
        };
        b_ret = a_r < a_tol3d;
    }
    b_ret
}

// ============================================================================
// OCCT GeomInt_LineConstructor.cxx L926-983: file-static RejectDuplicates
// ============================================================================

/// OCCT `static void RejectDuplicates(theVtxArr)` (cxx L926-983).
///
/// index note: OCCT indexes the array 1..Upper; the loops below use 0-based
/// indices (`ii = i - 1`) with the same bounds.
fn reject_duplicates(the_vtx_arr: &mut [GeomIntVertex]) {
    // About the value aTolPC = 1000. * Precision::PConfusion(), see
    // IntPatch_GLine::ComputeVertexParameters(...) for more details.
    const A_TOL_PC: f64 = 1000. * PCONFUSION;

    let upper = the_vtx_arr.len();

    // OCCT: for (i = Lower; i <= Upper - 2; i++)
    let mut ii = 0usize;
    while ii + 3 <= upper {
        let a_prmi = the_vtx_arr[ii].get_vertex().param_on_line;

        if a_prmi == REAL_LAST {
            ii += 1;
            continue;
        }

        // OCCT: for (j = i + 1; j <= Upper - 1; j++)
        let mut jj = ii + 1;
        while jj + 2 <= upper {
            let a_prmj = the_vtx_arr[jj].get_vertex().param_on_line;

            if a_prmj - a_prmi < A_TOL_PC {
                the_vtx_arr[jj].set_parameter(REAL_LAST);
            } else {
                break;
            }
            jj += 1;
        }
        ii += 1;
    }

    // Find duplicates with the last element of the array.
    // OCCT: for (i = Upper - 1; i > Lower; i--) with 1-based i.
    let a_max_prm = the_vtx_arr[upper - 1].get_vertex().param_on_line;
    let mut i = upper - 1; // 1-based Upper - 1
    while i > 1 {
        let idx = i - 1; // 0-based
        let a_prmi = the_vtx_arr[idx].get_vertex().param_on_line;

        if a_prmi == REAL_LAST {
            i -= 1;
            continue;
        }

        if (a_max_prm - a_prmi) < A_TOL_PC {
            the_vtx_arr[idx].set_parameter(REAL_LAST);
        } else {
            break;
        }
        i -= 1;
    }
}

/// OCCT `Epsilon(theValue)` (Standard_Real.hxx L242-246) — the ULP of
/// `theValue` toward the infinity of the same sign.
fn standard_epsilon(x: f64) -> f64 {
    if x >= 0.0 {
        next_after(x, f64::INFINITY) - x
    } else {
        x - next_after(x, f64::NEG_INFINITY)
    }
}

/// OCCT `std::nextafter` — bit-level construction (see the kernel twin).
fn next_after(x: f64, to: f64) -> f64 {
    if x.is_nan() || to.is_nan() || x == to {
        return x;
    }
    if x == 0.0 {
        return if to > 0.0 { f64::from_bits(1) } else { -f64::from_bits(1) };
    }
    if (to > x) == (x > 0.0) {
        f64::from_bits(x.to_bits() + 1)
    } else {
        f64::from_bits(x.to_bits() - 1)
    }
}
