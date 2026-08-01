// OCCT IntPatch (TKGeomAlgo) — intersection of patches (surfaces + domains).
//
// IntPatch_Point: a point of the intersection of two patches, with the
// transition relative to the restriction arcs / vertices of both surfaces.
// IntPatch_IType / IntPatch_Line: the abstract intersection line types.
//
// rcad data-model notes:
// - OCCT arcS1/arcS2 (Handle(Adaptor2d_Curve2d)) and vS1/vS2
//   (Handle(Adaptor3d_HVertex)) are restriction entities of the surfaces'
//   domains. rcad stores optional indices that the caller resolves against
//   its domain representation (0 = none).
// - IntSurf_Transition maps to crate::topalgo::int_res2d::Transition.

pub mod i_type;
pub mod line;

pub use i_type::IntPatchIType;
pub use line::IntPatchLine;

use glam::DVec3;

use crate::topalgo::int_res2d::{Position, Situation, Transition, TypeTrans};
use crate::topalgo::int_surf::PntOn2S;

/// OCCT IntPatch_Point (IntPatch_Point.cxx / .hxx) — a point on an
/// intersection line, possibly on a restriction arc/vertex of either surface.
#[derive(Debug, Clone)]
pub struct IntPatchPoint {
    /// IntSurf_PntOn2S — the point with its UV on both surfaces.
    pt: PntOn2S,
    /// Parameter on the intersection line.
    para: f64,
    /// Tolerance of the point.
    tol: f64,
    /// True when the intersection is tangent at this point.
    tgt: bool,
    /// True when the point is multiple (more than two branches).
    mult: bool,
    /// Point lies on the domain boundary (arc) of surface 1.
    on_s1: bool,
    /// Point is a vertex of the domain of surface 1.
    vtx_on_s1: bool,
    /// Arc of the domain of surface 1 the point lies on.
    arc_s1: Option<usize>,
    /// Transition of the line relative to arc_s1.
    traline1: Transition,
    /// Transition of the point relative to arc_s1.
    tra1: Transition,
    /// Parameter of the point on arc_s1.
    prm1: f64,
    /// Point lies on the domain boundary (arc) of surface 2.
    on_s2: bool,
    /// Point is a vertex of the domain of surface 2.
    vtx_on_s2: bool,
    /// Arc of the domain of surface 2 the point lies on.
    arc_s2: Option<usize>,
    /// Transition of the line relative to arc_s2.
    traline2: Transition,
    /// Transition of the point relative to arc_s2.
    tra2: Transition,
    /// Parameter of the point on arc_s2.
    prm2: f64,
}

impl IntPatchPoint {
    pub fn new() -> Self {
        IntPatchPoint {
            pt: PntOn2S::new(),
            para: 0.0,
            tol: 0.0,
            tgt: false,
            mult: false,
            on_s1: false,
            vtx_on_s1: false,
            arc_s1: None,
            traline1: Transition::undecided(Position::Middle),
            tra1: Transition::undecided(Position::Middle),
            prm1: 0.0,
            on_s2: false,
            vtx_on_s2: false,
            arc_s2: None,
            traline2: Transition::undecided(Position::Middle),
            tra2: Transition::undecided(Position::Middle),
            prm2: 0.0,
        }
    }

    /// OCCT SetValue(Pt, Tol, Tangent) (IntPatch_Point.cxx L25-35).
    pub fn set_value(&mut self, pt: DVec3, tol: f64, tangent: bool) {
        self.on_s1 = false;
        self.on_s2 = false;
        self.vtx_on_s1 = false;
        self.vtx_on_s2 = false;
        self.mult = false;
        self.tgt = tangent;
        self.pt.set_value(pt, true, 0.0, 0.0);
        self.pt.set_value_uv(false, 0.0, 0.0);
        self.tol = tol;
    }

    /// OCCT SetValue(Pt).
    pub fn set_value_pnt(&mut self, pt: DVec3) {
        self.pt.set_value(pt, true, 0.0, 0.0);
        self.pt.set_value_uv(false, 0.0, 0.0);
    }

    /// OCCT SetValue(PntOn2S).
    pub fn set_value_pnt_on_2s(&mut self, p: &PntOn2S) {
        self.pt = p.clone();
    }

    /// OCCT SetTolerance.
    pub fn set_tolerance(&mut self, tol: f64) {
        self.tol = tol;
    }

    /// OCCT SetParameters(U1, V1, U2, V2).
    pub fn set_parameters(&mut self, u1: f64, v1: f64, u2: f64, v2: f64) {
        self.pt.set_value_uv(true, u1, v1);
        self.pt.set_value_uv(false, u2, v2);
    }

    /// OCCT SetParameter(Para) — parameter on the intersection line.
    pub fn set_parameter(&mut self, para: f64) {
        self.para = para;
    }

    /// OCCT SetVertex(OnFirst, V).
    pub fn set_vertex(&mut self, on_first: bool, v: Option<usize>) {
        if on_first {
            self.on_s1 = true;
            self.vtx_on_s1 = true;
            let _ = v; // vS1 — vertex id resolved by caller
        } else {
            self.on_s2 = true;
            self.vtx_on_s2 = true;
            let _ = v;
        }
    }

    /// OCCT SetArc(OnFirst, A, Param, TLine, TArc).
    pub fn set_arc(
        &mut self, on_first: bool, arc: Option<usize>, param: f64,
        t_line: Transition, t_arc: Transition,
    ) {
        if on_first {
            self.on_s1 = true;
            self.arc_s1 = arc;
            self.traline1 = t_line;
            self.tra1 = t_arc;
            self.prm1 = param;
        } else {
            self.on_s2 = true;
            self.arc_s2 = arc;
            self.traline2 = t_line;
            self.tra2 = t_arc;
            self.prm2 = param;
        }
    }

    /// OCCT SetMultiple.
    pub fn set_multiple(&mut self, is_mult: bool) {
        self.mult = is_mult;
    }

    /// OCCT ParameterOnLine.
    pub fn parameter_on_line(&self) -> f64 {
        self.para
    }

    /// OCCT Tolerance.
    pub fn tolerance(&self) -> f64 {
        self.tol
    }

    /// OCCT IsTangencyPoint.
    pub fn is_tangency_point(&self) -> bool {
        self.tgt
    }

    /// OCCT IsMultiple.
    pub fn is_multiple(&self) -> bool {
        self.mult
    }

    /// OCCT IsOnDomS1.
    pub fn is_on_dom_s1(&self) -> bool {
        self.on_s1
    }

    /// OCCT IsVertexOnS1.
    pub fn is_vertex_on_s1(&self) -> bool {
        self.vtx_on_s1
    }

    /// OCCT ParameterOnArc1.
    pub fn parameter_on_arc1(&self) -> f64 {
        self.prm1
    }

    /// OCCT IsOnDomS2.
    pub fn is_on_dom_s2(&self) -> bool {
        self.on_s2
    }

    /// OCCT IsVertexOnS2.
    pub fn is_vertex_on_s2(&self) -> bool {
        self.vtx_on_s2
    }

    /// OCCT ParameterOnArc2.
    pub fn parameter_on_arc2(&self) -> f64 {
        self.prm2
    }

    /// OCCT PntOn2S() — the underlying point.
    pub fn pnt_on_2s(&self) -> &PntOn2S {
        &self.pt
    }

    /// OCCT ReverseTransition (IntPatch_Point.cxx L77-137) — swap In/Out
    /// transitions when the line orientation is reversed.
    pub fn reverse_transition(&mut self) {
        if self.on_s1 {
            self.traline1 = reverse_type(&self.traline1);
            self.tra1 = reverse_type(&self.tra1);
        }
        if self.on_s2 {
            self.traline2 = reverse_type(&self.traline2);
            self.tra2 = reverse_type(&self.tra2);
        }
    }
}

impl Default for IntPatchPoint {
    fn default() -> Self {
        Self::new()
    }
}

/// OCCT IntPatch_Point::ReverseTransition — invert a In/Out transition.
fn reverse_type(t: &Transition) -> Transition {
    let mut out = t.clone();
    match t.transition_type() {
        TypeTrans::In => {
            out = Transition::in_out(t.is_tangent(), Position::Middle, TypeTrans::Out);
        }
        TypeTrans::Out => {
            out = Transition::in_out(t.is_tangent(), Position::Middle, TypeTrans::In);
        }
        TypeTrans::Touch => {
            // OCCT: TArc.SetValue(false, IntSurf_Out) — only in/out handled;
            // for Touch the transition is left unchanged by the OCCT switch.
            let _ = Situation::Unknown;
        }
        _ => {}
    }
    out
}
