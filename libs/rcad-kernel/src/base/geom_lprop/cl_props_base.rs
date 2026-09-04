// OCCT GeomLProp_CLPropsBase + LProp_CurveUtils (TKGeomBase) — the generic
// curve local-properties engine: point, derivatives up to order 3, tangent,
// curvature, normal, centre of curvature.
//
// GeomLProp_CLProps.hxx L46-216 (the ClPropsBase template, parameterized by
// point/vector/direction types and a curve access policy) +
// LProp_CurveUtils.hxx L26-466 (the shared algorithms: EvalDerivatives,
// SetParameter, EnsureDeriv, IsTangentDefined, ComputeTangent, Tangent,
// ComputeCurvature/Curvature, ComputeNormal/Normal,
// ComputeCentreOfCurvature/CentreOfCurvature).
//
// The OCCT template parameters (Pnt, Vec, Dir, CurveType, Access) map to a
// single accessor trait [`CLPropsCurve2d`] for the 2D instantiation
// (gp_Pnt2d/gp_Vec2d/gp_Dir2d + LProp_CurveUtils::ToolAccess<Tool> — the
// tool's statics become the trait methods).  The HLRBRep_CLProps
// instantiation lives in rcad-algo `hlr/brep/cl_props.rs`.

use crate::base::geom_lprop::LPropStatus;

/// OCCT gp::Resolution() floor for the line-parameter test.
const GP_RESOLUTION: f64 = 1e-15;

/// The `LProp_CurveUtils::ToolAccess<Tool>` (or DirectAccess) policy for a
/// 2D curve: the D0..D3 evaluations and the parameter bounds.
pub trait CLPropsCurve2d {
    /// OCCT Tool::Value(C, U, P) / Access::D0.
    fn eval_d0(&self, u: f64) -> DVec2;
    /// OCCT Tool::D1(C, U, P, V1) / Access::D1.
    fn eval_d1(&self, u: f64) -> (DVec2, DVec2);
    /// OCCT Tool::D2(C, U, P, V1, V2) / Access::D2.
    fn eval_d2(&self, u: f64) -> (DVec2, DVec2, DVec2);
    /// OCCT Tool::D3(C, U, P, V1, V2, V3) / Access::D3.
    fn eval_d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2);
    /// OCCT Tool::FirstParameter(C) / Access::FirstParameter.
    fn first_parameter(&self) -> f64;
    /// OCCT Tool::LastParameter(C) / Access::LastParameter.
    fn last_parameter(&self) -> f64;
}

use glam::DVec2;

/// OCCT GeomLProp_CLPropsBase (2D instantiation: gp_Pnt2d/gp_Vec2d/gp_Dir2d).
pub struct ClPropsBase<'a, C: ?Sized> {
    /// OCCT myCurve.
    curve: Option<&'a C>,
    /// OCCT myU.
    u: f64,
    /// OCCT myDerOrder.
    der_order: i32,
    /// OCCT myCN.
    cn: i32,
    /// OCCT myLinTol.
    lin_tol: f64,
    /// OCCT myPnt.
    pnt: DVec2,
    /// OCCT myDerivArr[3].
    deriv: [DVec2; 3],
    /// OCCT myCurvature.
    curvature: f64,
    /// OCCT myTangentStatus.
    tangent_status: LPropStatus,
    /// OCCT mySignificantFirstDerivativeOrder.
    significant_first_derivative_order: i32,
}

impl<'a, C: CLPropsCurve2d + ?Sized> ClPropsBase<'a, C> {
    /// OCCT GeomLProp_CLPropsBase(N, Resolution) (hxx L96-105) — the curve
    /// is set later with SetCurve.
    pub fn new(n: i32, resolution: f64) -> Self {
        assert!((0..=3).contains(&n), "Standard_OutOfRange: CLProps(N)");
        ClPropsBase {
            curve: None,
            u: f64::MAX, // RealLast()
            der_order: n,
            cn: 0,
            lin_tol: resolution,
            pnt: DVec2::ZERO,
            deriv: [DVec2::ZERO; 3],
            curvature: 0.0,
            tangent_status: LPropStatus::Undecided,
            significant_first_derivative_order: 0,
        }
    }

    /// OCCT SetCurve(C) (hxx L122-126).
    pub fn set_curve(&mut self, c: &'a C) {
        self.curve = Some(c);
        self.cn = 4;
    }

    fn c(&self) -> &'a C {
        self.curve.expect("CLProps: no curve set")
    }

    /// OCCT SetParameter(U) via LProp_CurveUtils::SetParameter
    /// (LProp_CurveUtils.hxx L287-298 + EvalDerivatives L146-163).
    pub fn set_parameter(&mut self, u: f64) {
        self.u = u;
        // EvalDerivatives: one combined evaluation up to myDerOrder.
        match self.der_order {
            0 => self.pnt = self.c().eval_d0(u),
            1 => {
                let (p, v1) = self.c().eval_d1(u);
                self.pnt = p;
                self.deriv[0] = v1;
            }
            2 => {
                let (p, v1, v2) = self.c().eval_d2(u);
                self.pnt = p;
                self.deriv[0] = v1;
                self.deriv[1] = v2;
            }
            3 => {
                let (p, v1, v2, v3) = self.c().eval_d3(u);
                self.pnt = p;
                self.deriv[0] = v1;
                self.deriv[1] = v2;
                self.deriv[2] = v3;
            }
            _ => {}
        }
        self.tangent_status = LPropStatus::Undecided;
    }

    /// OCCT Value() (hxx L129).
    pub fn value(&self) -> DVec2 {
        self.pnt
    }

    /// OCCT D1() via EnsureDeriv (LProp_CurveUtils.hxx L309-322).
    pub fn d1(&mut self) -> DVec2 {
        if self.der_order < 1 {
            self.der_order = 1;
            let (p, v1) = self.c().eval_d1(self.u);
            self.pnt = p;
            self.deriv[0] = v1;
        }
        self.deriv[0]
    }

    /// OCCT D2() via EnsureDeriv.
    pub fn d2(&mut self) -> DVec2 {
        if self.der_order < 2 {
            self.der_order = 2;
            let (p, v1, v2) = self.c().eval_d2(self.u);
            self.pnt = p;
            self.deriv[0] = v1;
            self.deriv[1] = v2;
        }
        self.deriv[1]
    }

    /// OCCT D3() via EnsureDeriv.
    pub fn d3(&mut self) -> DVec2 {
        if self.der_order < 3 {
            self.der_order = 3;
            let (p, v1, v2, v3) = self.c().eval_d3(self.u);
            self.pnt = p;
            self.deriv[0] = v1;
            self.deriv[1] = v2;
            self.deriv[2] = v3;
        }
        self.deriv[2]
    }

    /// OCCT IsTangentDefined via LProp_CurveUtils::IsTangentDefined
    /// (LProp_CurveUtils.hxx L333-380).
    pub fn is_tangent_defined(&mut self) -> bool {
        if self.tangent_status == LPropStatus::Undefined {
            return false;
        }
        if self.tangent_status as i32 >= LPropStatus::Defined as i32 {
            return true;
        }

        let a_tol_sq = self.lin_tol * self.lin_tol;
        let mut an_order = 0;
        while an_order < 4 {
            an_order += 1;
            if self.cn >= an_order {
                let a_v = match an_order {
                    1 => self.d1(),
                    2 => self.d2(),
                    3 => self.d3(),
                    _ => {
                        self.tangent_status = LPropStatus::Undefined;
                        return false;
                    }
                };
                if a_v.length_squared() > a_tol_sq {
                    self.significant_first_derivative_order = an_order;
                    self.tangent_status = LPropStatus::Defined;
                    return true;
                }
            } else {
                self.tangent_status = LPropStatus::Undefined;
                return false;
            }
        }
        false
    }

    /// OCCT Tangent(D) via LProp_CurveUtils::Tangent + ComputeTangent
    /// (LProp_CurveUtils.hxx L390-402 + L174-219).  Returns None for the
    /// LProp_NotDefined throw.
    pub fn tangent(&mut self) -> Option<DVec2> {
        if !self.is_tangent_defined() {
            return None; // LProp_NotDefined
        }
        // ComputeTangent.
        if self.significant_first_derivative_order == 1 {
            return Some(self.deriv[0].normalize_or_zero()); // Dir(aV)
        }

        const THE_DIVISION_FACTOR: f64 = 1.0e-3;
        const THE_MIN_STEP: f64 = 1.0e-7;

        let an_usupremum = self.c().last_parameter();
        let an_uinfimum = self.c().first_parameter();

        let a_du = if an_usupremum >= f64::MAX || an_uinfimum <= f64::MIN {
            0.0
        } else {
            an_usupremum - an_uinfimum
        };

        let a_delta = (a_du * THE_DIVISION_FACTOR).max(THE_MIN_STEP);

        let mut a_v = self.deriv[(self.significant_first_derivative_order - 1) as usize];

        let an_other_u = if self.u - an_uinfimum < a_delta {
            self.u + a_delta
        } else {
            self.u - a_delta
        };

        let p1 = self.c().eval_d0(self.u.min(an_other_u));
        let p2 = self.c().eval_d0(self.u.max(an_other_u));

        let a_chord = p2 - p1;
        if a_v.dot(a_chord) < 0.0 {
            a_v = -a_v;
        }

        Some(a_v.normalize_or_zero()) // Dir(aV)
    }

    /// OCCT Curvature() via LProp_CurveUtils::Curvature + ComputeCurvature
    /// (LProp_CurveUtils.hxx L412-427 + L227-242).  Returns RealLast() when
    /// the significant derivative is of higher order.
    pub fn curvature(&mut self) -> f64 {
        let _an_is_defined = self.is_tangent_defined();
        if self.significant_first_derivative_order > 1 {
            return f64::MAX; // RealLast()
        }
        let the_d1 = self.deriv[0];
        let the_d2 = self.deriv[1];
        let a_tol_sq = self.lin_tol * self.lin_tol;
        let a_dd1 = the_d1.length_squared();
        let a_dd2 = the_d2.length_squared();

        if a_dd2 <= a_tol_sq {
            self.curvature = 0.0;
            return 0.0;
        }

        let a_n = {
            // CrossSquareMagnitude for 2D: (D1 x D2)^2 = (x1*y2 - y1*x2)^2.
            let cross = the_d1.x * the_d2.y - the_d1.y * the_d2.x;
            cross * cross
        };
        let a_t = a_n / a_dd1 / a_dd2;
        if a_t <= a_tol_sq {
            self.curvature = 0.0;
            return 0.0;
        }

        self.curvature = a_n.sqrt() / a_dd1 / a_dd1.sqrt();
        self.curvature
    }

    /// OCCT Normal(N) via LProp_CurveUtils::Normal + ComputeNormal
    /// (LProp_CurveUtils.hxx L435-442 + L249-254).  Returns None for the
    /// LProp_NotDefined throw.
    pub fn normal(&mut self) -> Option<DVec2> {
        let a_curvature = self.curvature();
        if a_curvature == f64::MAX || a_curvature.abs() <= self.lin_tol {
            return None; // LProp_NotDefined
        }
        let the_d1 = self.deriv[0];
        let the_d2 = self.deriv[1];
        // Vec aNorm = D2*(D1*D1) - D1*(D1*D2); Dir = Dir(aNorm).
        let a_norm = the_d2 * the_d1.length_squared() - the_d1 * the_d1.dot(the_d2);
        Some(a_norm.normalize_or_zero())
    }

    /// OCCT CentreOfCurvature(P) via LProp_CurveUtils::CentreOfCurvature +
    /// ComputeCentreOfCurvature (LProp_CurveUtils.hxx L452-464 + L263-274).
    /// Returns None for the LProp_NotDefined throw.
    pub fn centre_of_curvature(&mut self) -> Option<DVec2> {
        if self.curvature().abs() <= self.lin_tol {
            return None; // LProp_NotDefined
        }
        let the_pnt = self.pnt;
        let the_d1 = self.deriv[0];
        let the_d2 = self.deriv[1];
        let mut a_norm = the_d2 * the_d1.length_squared() - the_d1 * the_d1.dot(the_d2);
        a_norm = a_norm.normalize_or_zero();
        a_norm /= self.curvature;
        Some(the_pnt + a_norm) // thePnt.Translated(aNorm)
    }
}

// The gp::Resolution floor reference (used by the callers constructing the
// tolerance, Epsilon(1.) in HLRBRep_Curve::Tangent).
#[allow(unused)]
fn _parity() -> f64 {
    GP_RESOLUTION
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A unit-speed line y = t along X: D1 = (1, 0), higher derivatives 0.
    struct TestLine2d;

    impl CLPropsCurve2d for TestLine2d {
        fn eval_d0(&self, u: f64) -> DVec2 {
            DVec2::new(u, 0.0)
        }
        fn eval_d1(&self, u: f64) -> (DVec2, DVec2) {
            (DVec2::new(u, 0.0), DVec2::new(1.0, 0.0))
        }
        fn eval_d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
            (DVec2::new(u, 0.0), DVec2::new(1.0, 0.0), DVec2::ZERO)
        }
        fn eval_d3(&self, u: f64) -> (DVec2, DVec2, DVec2, DVec2) {
            (DVec2::new(u, 0.0), DVec2::new(1.0, 0.0), DVec2::ZERO, DVec2::ZERO)
        }
        fn first_parameter(&self) -> f64 {
            0.0
        }
        fn last_parameter(&self) -> f64 {
            2.0
        }
    }

    /// OCCT anchor: a straight line — the tangent is defined from D1, the
    /// curvature is 0.
    #[test]
    fn clprops_line_tangent_curvature() {
        let mut clp = ClPropsBase::<dyn CLPropsCurve2d>::new(2, 1e-12);
        clp.set_curve(&TestLine2d);
        clp.set_parameter(1.0);
        let p = clp.value();
        assert!((p.x - 1.0).abs() < 1e-12 && p.y.abs() < 1e-12);
        assert!(clp.is_tangent_defined());
        let d = clp.tangent().unwrap();
        assert!((d.x - 1.0).abs() < 1e-12 && d.y.abs() < 1e-12);
        assert!(clp.curvature().abs() < 1e-12);
    }
}
