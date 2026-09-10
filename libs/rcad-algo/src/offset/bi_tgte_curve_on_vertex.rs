// OCCT BiTgte_CurveOnVertex.cxx L1-253 (+ BiTgte_CurveOnVertex.hxx L42-117)
// — 1:1 translation.
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BiTgte/
//         BiTgte_CurveOnVertex.cxx / .hxx
//
// OCCT inheritance chain (hxx L42): BiTgte_CurveOnVertex -> Adaptor3d_Curve.
// Architecture difference: rcad has no Adaptor3d_Curve class hierarchy —
// the adaptor is carried as a plain struct with the OCCT method names; the
// methods that throw Standard_NotImplemented in OCCT keep the same panic.

use glam::DVec3;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{
    BezierCurve3, BSplineCurve3, Circle3, Curve3, Ellipse3, Hyperbola3, Line3, Parabola3,
};
use rcad_kernel::topo_shape::Shape;

use super::brep_offset_offset::GeomAbsShapeKind;
use super::bi_tgte_curve_on_edge::{ResD1, ResD2, ResD3};

use crate::feat::loc_ope_wires_on_shape::brep_tool_pnt;

// ---------------------------------------------------------------------------
// OCCT class (BiTgte_CurveOnVertex.hxx L42-117, cxx L38-253).
// ---------------------------------------------------------------------------

/// OCCT BiTgte_CurveOnVertex (BiTgte_CurveOnVertex.hxx L42-117).
#[derive(Clone)]
pub struct BiTgteCurveOnVertex {
    my_first: f64,  // OCCT: myFirst (hxx L104)
    my_last: f64,   // OCCT: myLast (hxx L105)
    my_pnt: DVec3,  // OCCT: myPnt (hxx L106)
}

impl Default for BiTgteCurveOnVertex {
    fn default() -> Self {
        Self::new()
    }
}

impl BiTgteCurveOnVertex {
    /// OCCT BiTgte_CurveOnVertex::BiTgte_CurveOnVertex() (cxx L40-44).
    pub fn new() -> Self {
        BiTgteCurveOnVertex {
            my_first: 0.0,
            my_last: 0.0,
            my_pnt: DVec3::ZERO,
        }
    }

    /// OCCT BiTgte_CurveOnVertex::BiTgte_CurveOnVertex(EonF, Vertex)
    /// (cxx L48-54).
    pub fn with_edge_vertex(the_eon_f: &Shape, the_vertex: &Shape) -> Self {
        let mut res = BiTgteCurveOnVertex::new();
        res.init(the_eon_f, the_vertex);
        res
    }

    /// OCCT BiTgte_CurveOnVertex::Init(EonF, V) (cxx L58-62).
    pub fn init(&mut self, eon_f: &Shape, v: &Shape) {
        // OCCT L60: BRep_Tool::Range(EonF, myFirst, myLast) — the rcad edge
        // range read (TEdgeData::range).
        match eon_f.as_edge() {
            Some(ed) => {
                self.my_first = ed.range[0];
                self.my_last = ed.range[1];
            }
            None => panic!("BiTgte_CurveOnVertex::Init: EonF is not an edge"),
        }
        // OCCT L61: myPnt = BRep_Tool::Pnt(V).
        self.my_pnt = brep_tool_pnt(v).expect("BiTgte_CurveOnVertex::Init: null vertex");
    }

    /// OCCT BiTgte_CurveOnVertex::FirstParameter() (cxx L66-69).
    pub fn first_parameter(&self) -> f64 {
        self.my_first
    }

    /// OCCT BiTgte_CurveOnVertex::LastParameter() (cxx L73-76).
    pub fn last_parameter(&self) -> f64 {
        self.my_last
    }

    /// OCCT BiTgte_CurveOnVertex::Continuity() (cxx L80-83).
    pub fn continuity(&self) -> GeomAbsShapeKind {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::NbIntervals(S) (cxx L87-90).
    pub fn nb_intervals(&self, _s: GeomAbsShapeKind) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Intervals(T, S) (cxx L94-97).
    pub fn intervals(&self, _t: &mut [f64], _s: GeomAbsShapeKind) {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Trim(F, L, Tol) (cxx L101-106).
    pub fn trim(&self, _first: f64, _last: f64, _tol: f64) -> Curve3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::IsClosed() (cxx L110-113).
    pub fn is_closed(&self) -> bool {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::IsPeriodic() (cxx L116-119).
    pub fn is_periodic(&self) -> bool {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Period() (cxx L124-127).
    pub fn period(&self) -> f64 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::EvalD0(theU) (cxx L131-134).
    pub fn eval_d0(&self, _the_u: f64) -> DVec3 {
        self.my_pnt
    }

    /// OCCT BiTgte_CurveOnVertex::EvalD1(theU) (cxx L138-141).
    pub fn eval_d1(&self, _the_u: f64) -> ResD1 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::EvalD2(theU) (cxx L145-148).
    pub fn eval_d2(&self, _the_u: f64) -> ResD2 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::EvalD3(theU) (cxx L152-155).
    pub fn eval_d3(&self, _the_u: f64) -> ResD3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::EvalDN(theU, theN) (cxx L159-162).
    pub fn eval_dn(&self, _the_u: f64, _the_n: i32) -> DVec3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Resolution(theR3d) (cxx L166-169).
    pub fn resolution(&self, _the_r3d: f64) -> f64 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::GetType() (cxx L173-176).
    pub fn get_type(&self) -> CurveType {
        CurveType::Other
    }

    /// OCCT BiTgte_CurveOnVertex::Line() (cxx L180-183).
    pub fn line(&self) -> Line3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Circle() (cxx L187-190).
    pub fn circle(&self) -> Circle3 {
        // OCCT L189: throw Standard_NoSuchObject.
        panic!("NoSuchObject: BiTgte_CurveOnVertex::Circle");
    }

    /// OCCT BiTgte_CurveOnVertex::Ellipse() (cxx L194-197).
    pub fn ellipse(&self) -> Ellipse3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Hyperbola() (cxx L201-204).
    pub fn hyperbola(&self) -> Hyperbola3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Parabola() (cxx L208-211).
    pub fn parabola(&self) -> Parabola3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Degree() (cxx L215-218).
    pub fn degree(&self) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::IsRational() (cxx L222-225).
    pub fn is_rational(&self) -> bool {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::NbPoles() (cxx L229-232).
    pub fn nb_poles(&self) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::NbKnots() (cxx L236-239).
    pub fn nb_knots(&self) -> i32 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::Bezier() (cxx L243-246).
    pub fn bezier(&self) -> BezierCurve3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }

    /// OCCT BiTgte_CurveOnVertex::BSpline() (cxx L250-253).
    pub fn bspline(&self) -> BSplineCurve3 {
        panic!("NotImplemented: BiTgte_CurveOnVertex");
    }
}
