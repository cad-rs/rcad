// OCCT HLRBRep_CLProps + HLRBRep_CLPropsATool (TKHLR) — the curve local
// properties engine instantiated over HLRBRep_Curve.
//
// HLRBRep_CLProps.hxx L25-29:
//   using HLRBRep_CLProps = GeomLProp_CLPropsBase<gp_Pnt2d, gp_Vec2d,
//       gp_Dir2d, const HLRBRep_Curve*, LProp_CurveUtils::ToolAccess<
//       HLRBRep_CLPropsATool>>;
//
// HLRBRep_CLPropsATool.hxx L29-69 + .lxx (the tool statics delegate to the
// HLRBRep_Curve methods; Continuity returns GeomAbs_C2).  The rcad tool is
// the `CLPropsCurve2d` impl for [`Curve`].

use rcad_kernel::base::geom_lprop::{CLPropsCurve2d, ClPropsBase};

use super::curve::Curve;

/// OCCT HLRBRep_CLProps — the 2D local-properties engine over the projected
/// curve.
pub type CLProps<'a> = ClPropsBase<'a, Curve<'a>>;

/// OCCT HLRBRep_CLPropsATool (the ToolAccess policy for HLRBRep_Curve).
impl CLPropsCurve2d for Curve<'_> {
    /// OCCT Tool::Value(C, U, P) = A->Value(U) (lxx).
    fn eval_d0(&self, u: f64) -> glam::DVec2 {
        self.value(u)
    }

    /// OCCT Tool::D1(C, U, P, V1) = A->D1(U, P, V1) (lxx).
    fn eval_d1(&self, u: f64) -> (glam::DVec2, glam::DVec2) {
        let mut p = glam::DVec2::ZERO;
        let mut v = glam::DVec2::ZERO;
        self.d1_2d(u, &mut p, &mut v);
        (p, v)
    }

    /// OCCT Tool::D2(C, U, P, V1, V2) = A->D2(U, P, V1, V2) (lxx).
    fn eval_d2(&self, u: f64) -> (glam::DVec2, glam::DVec2, glam::DVec2) {
        let mut p = glam::DVec2::ZERO;
        let mut v1 = glam::DVec2::ZERO;
        let mut v2 = glam::DVec2::ZERO;
        self.d2_2d(u, &mut p, &mut v1, &mut v2);
        (p, v1, v2)
    }

    /// OCCT Tool::D3(C, U, P, V1, V2, V3) = A->D3(U, P, V1, V2, V3) — the
    /// HLRBRep_Curve::D3 body is empty (cxx L406), so the out values keep
    /// their default-constructed (0,0) state, verbatim.
    fn eval_d3(&self, _u: f64) -> (glam::DVec2, glam::DVec2, glam::DVec2, glam::DVec2) {
        (glam::DVec2::ZERO, glam::DVec2::ZERO, glam::DVec2::ZERO, glam::DVec2::ZERO)
    }

    /// OCCT Tool::FirstParameter(C) = A->FirstParameter() (lxx).
    fn first_parameter(&self) -> f64 {
        Curve::first_parameter(self)
    }

    /// OCCT Tool::LastParameter(C) = A->LastParameter() (lxx).
    fn last_parameter(&self) -> f64 {
        Curve::last_parameter(self)
    }
}

// The OCCT ATool::Continuity returns GeomAbs_C2; the generic engine carries
// myCN = 4 from SetCurve and never queries the tool continuity (the same
// dead-member status as OCCT).
#[allow(unused)]
fn _continuity_parity(_c: &Curve<'_>) -> i32 {
    2 // GeomAbs_C2
}
