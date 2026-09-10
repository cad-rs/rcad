//! OCCT Blend_SurfRstFunction (TKFillet/Blend) — 1:1 port of
//! Blend_SurfRstFunction.hxx (L40-223) + Blend_SurfRstFunction.cxx
//! (L19-34).
//!
//! Architecture mapping: `class Blend_SurfRstFunction : public
//! Blend_AppFunction` maps to a Rust subtrait of [`BlendAppFunction`].  The
//! default-implemented trait methods correspond to the OCCT methods that are
//! implemented (not deferred) in Blend_SurfRstFunction.cxx.

use glam::{DVec2, DVec3};

use super::brep_blend_function::BlendAppFunction;

/// OCCT Blend_SurfRstFunction — deferred class for a function used to compute
/// a blending surface between a surface and a pcurve on the other surface,
/// using a guide line (Blend_SurfRstFunction.hxx L34).  The vector X is
/// U, V, W: the parametric coordinates of the extremities of a section on
/// the surface and the curve.
pub trait BlendSurfRstFunction: BlendAppFunction {
    /// OCCT NbVariables() — returns 3 (deferred in OCCT).
    fn nb_variables(&self) -> usize;

    /// OCCT NbEquations() — the number of equations of the function.
    fn nb_equations(&self) -> usize;

    /// OCCT Pnt1() (Blend_SurfRstFunction.cxx L21-24) — delegates to
    /// PointOnS.
    fn pnt1(&self) -> DVec3 {
        self.point_on_s()
    }

    /// OCCT Pnt2() (Blend_SurfRstFunction.cxx L26-29) — delegates to
    /// PointOnRst.
    fn pnt2(&self) -> DVec3 {
        self.point_on_rst()
    }

    /// OCCT PointOnS() — the point on the surface.
    fn point_on_s(&self) -> DVec3;

    /// OCCT PointOnRst() — the point on the restriction.
    fn point_on_rst(&self) -> DVec3;

    /// OCCT Pnt2dOnS() — the U,V coordinates of the point on the surface.
    fn pnt2d_on_s(&self) -> DVec2;

    /// OCCT Pnt2dOnRst() (Blend_SurfRstFunction.hxx L123) — the U,V
    /// coordinates of the point on the curve on surface.
    fn pnt2d_on_rst(&self) -> DVec2;

    /// OCCT ParameterOnRst() — the parameter of the point on the
    /// restriction.
    fn parameter_on_rst(&self) -> f64;

    /// OCCT IsTangencyPoint() — true when it is not possible to compute the
    /// tangent vectors at PointOnS and/or PointOnRst.
    fn is_tangency_point(&self) -> bool;

    /// OCCT TangentOnS() — the tangent vector at PointOnS, in 3d space.
    fn tangent_on_s(&self) -> DVec3;

    /// OCCT Tangent2dOnS() — the tangent vector at PointOnS, in the
    /// parametric space of the surface.
    fn tangent_2d_on_s(&self) -> DVec2;

    /// OCCT TangentOnRst() — the tangent vector at PointOnRst, in 3d space.
    fn tangent_on_rst(&self) -> DVec3;

    /// OCCT Tangent2dOnRst() — the tangent vector at PointOnRst, in the
    /// parametric space of the surface.
    fn tangent_2d_on_rst(&self) -> DVec2;

    /// OCCT Decroch(Sol, NS, TgS) — the flag that indicates whether a
    /// decrochage (pull-away) occurs at the solution point.
    fn decroch(&self, sol: &[f64], ns: &mut DVec3, tg_s: &mut DVec3) -> bool;

    /// OCCT GetMinimalDistance() (Blend_SurfRstFunction.cxx L31-34) — throws
    /// Standard_NotImplemented in OCCT.
    fn get_minimal_distance(&self) -> f64 {
        panic!("Standard_NotImplemented: Blend_SurfRstFunction::GetMinimalDistance");
    }
}
