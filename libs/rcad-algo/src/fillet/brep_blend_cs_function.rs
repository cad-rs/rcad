//! OCCT Blend_CSFunction (TKFillet/Blend) — 1:1 port of Blend_CSFunction.hxx
//! (L40-199) + Blend_CSFunction.cxx (L19-52).
//!
//! Architecture mapping: `class Blend_CSFunction : public Blend_AppFunction`
//! maps to a Rust subtrait of [`BlendAppFunction`].  The default-implemented
//! trait methods correspond to the OCCT methods that are implemented (not
//! deferred) in Blend_CSFunction.cxx.  The three OCCT `Section` overloads are
//! named `section`, `section_d1` and `section_d2` (same convention as
//! [`super::brep_blend_function::BlendAppFunction`]).

use glam::{DVec2, DVec3};

use super::brep_blend_function::BlendAppFunction;
use super::brep_blend_point::BlendPoint;

/// OCCT Blend_CSFunction — deferred class for a function used to compute a
/// blending surface between a surface and a curve, using a guide line
/// (Blend_CSFunction.hxx L40).  The vector X is w, U, V.
pub trait BlendCSFunction: BlendAppFunction {
    /// OCCT NbVariables() (Blend_CSFunction.cxx L21-24) — returns 3.
    fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations() — the number of equations of the function.
    fn nb_equations(&self) -> usize;

    /// OCCT Pnt1() (Blend_CSFunction.cxx L26-29) — delegates to PointOnC.
    fn pnt1(&self) -> DVec3 {
        self.point_on_c()
    }

    /// OCCT Pnt2() (Blend_CSFunction.cxx L31-34) — delegates to PointOnS.
    fn pnt2(&self) -> DVec3 {
        self.point_on_s()
    }

    /// OCCT PointOnS() — the point on the surface.
    fn point_on_s(&self) -> DVec3;

    /// OCCT PointOnC() — the point on the curve.
    fn point_on_c(&self) -> DVec3;

    /// OCCT Pnt2d() — the U,V coordinates of the point on the surface.
    fn pnt2d(&self) -> DVec2;

    /// OCCT ParameterOnC() — the parameter of the point on the curve.
    fn parameter_on_c(&self) -> f64;

    /// OCCT IsTangencyPoint() — true when it is not possible to compute the
    /// tangent vectors at PointOnS and/or PointOnC.
    fn is_tangency_point(&self) -> bool;

    /// OCCT TangentOnS() — the tangent vector at PointOnS, in 3d space.
    fn tangent_on_s(&self) -> DVec3;

    /// OCCT Tangent2d() — the tangent vector at PointOnS, in the parametric
    /// space of the surface.
    fn tangent_2d(&self) -> DVec2;

    /// OCCT TangentOnC() — the tangent vector at PointOnC, in 3d space.
    fn tangent_on_c(&self) -> DVec3;

    /// OCCT Tangent(U, V, TgS, NormS) — the tangent vector at the section,
    /// and the normal (of the surface) at these points.
    fn tangent(&self, u: f64, v: f64, tg_s: &mut DVec3, norm_s: &mut DVec3);

    /// OCCT TwistOnS() — always false (not redefined in the CS hierarchy).
    fn twist_on_s(&self) -> bool {
        false
    }

    /// OCCT TwistOnC() — always false (see TwistOnS).
    fn twist_on_c(&self) -> bool {
        false
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weigths, DWeigths, D2Weigths) (Blend_CSFunction.cxx L36-51) — the
    /// default returns false.
    #[allow(clippy::too_many_arguments)]
    fn section_d2(
        &mut self,
        _p: &BlendPoint,
        _poles: &mut [DVec3],
        _d_poles: &mut [DVec3],
        _d2_poles: &mut [DVec3],
        _poles_2d: &mut [DVec2],
        _d_poles_2d: &mut [DVec2],
        _d2_poles_2d: &mut [DVec2],
        _weigths: &mut [f64],
        _d_weigths: &mut [f64],
        _d2_weigths: &mut [f64],
    ) -> bool {
        false
    }

    /// OCCT GetMinimalDistance() (Blend_CSFunction.cxx L53-56) — throws
    /// Standard_NotImplemented in OCCT.
    fn get_minimal_distance(&self) -> f64 {
        panic!("Standard_NotImplemented: Blend_CSFunction::GetMinimalDistance");
    }
}
