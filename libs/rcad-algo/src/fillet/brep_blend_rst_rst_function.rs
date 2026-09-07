//! OCCT Blend_RstRstFunction (TKFillet/Blend) — 1:1 port of
//! Blend_RstRstFunction.hxx (L42-233) + Blend_RstRstFunction.cxx (L19-31).
//!
//! Architecture mapping: `class Blend_RstRstFunction : public
//! Blend_AppFunction` maps to a Rust subtrait of [`BlendAppFunction`].  The
//! default-implemented trait methods correspond to the OCCT methods that are
//! implemented (not deferred) in Blend_RstRstFunction.cxx.

use glam::{DVec2, DVec3};

use super::brep_blend_function::BlendAppFunction;

/// OCCT Blend_RstRstFunction — deferred class for a function used to compute
/// a blending surface between two pcurves (restrictions), using a guide line
/// (Blend_RstRstFunction.hxx L37).  The vector X is W, U1, V1, U2, V2.
pub trait BlendRstRstFunction: BlendAppFunction {
    /// OCCT NbVariables() — returns 5 (deferred in OCCT).
    fn nb_variables(&self) -> usize;

    /// OCCT NbEquations() — the number of equations of the function.
    fn nb_equations(&self) -> usize;

    /// OCCT Pnt1() (Blend_RstRstFunction.cxx L21-24) — delegates to
    /// PointOnRst1.
    fn pnt1(&self) -> DVec3 {
        self.point_on_rst1()
    }

    /// OCCT Pnt2() (Blend_RstRstFunction.cxx L26-29) — delegates to
    /// PointOnRst2.
    fn pnt2(&self) -> DVec3 {
        self.point_on_rst2()
    }

    /// OCCT PointOnRst1() — the point on the first restriction.
    fn point_on_rst1(&self) -> DVec3;

    /// OCCT PointOnRst2() — the point on the second restriction.
    fn point_on_rst2(&self) -> DVec3;

    /// OCCT Pnt2dOnRst1() — the U,V coordinates of the point on the first
    /// restriction.
    fn pnt2d_on_rst1(&self) -> DVec2;

    /// OCCT Pnt2dOnRst2() — the U,V coordinates of the point on the second
    /// restriction.
    fn pnt2d_on_rst2(&self) -> DVec2;

    /// OCCT ParameterOnRst1() — the parameter of the point on the first
    /// restriction.
    fn parameter_on_rst1(&self) -> f64;

    /// OCCT ParameterOnRst2() — the parameter of the point on the second
    /// restriction.
    fn parameter_on_rst2(&self) -> f64;

    /// OCCT IsTangencyPoint() — true when it is not possible to compute the
    /// tangent vectors at PointOnRst1 and/or PointOnRst2.
    fn is_tangency_point(&self) -> bool;

    /// OCCT TangentOnRst1() — the tangent vector at PointOnRst1, in 3d space.
    fn tangent_on_rst1(&self) -> DVec3;

    /// OCCT Tangent2dOnRst1() — the tangent vector at PointOnRst1, in the
    /// parametric space of the first surface.
    fn tangent_2d_on_rst1(&self) -> DVec2;

    /// OCCT TangentOnRst2() — the tangent vector at PointOnRst2, in 3d space.
    fn tangent_on_rst2(&self) -> DVec3;

    /// OCCT Tangent2dOnRst2() — the tangent vector at PointOnRst2, in the
    /// parametric space of the second surface.
    fn tangent_2d_on_rst2(&self) -> DVec2;

    /// OCCT GetMinimalDistance() (Blend_RstRstFunction.cxx L31) — throws
    /// Standard_NotImplemented in OCCT.
    fn get_minimal_distance(&self) -> f64 {
        panic!("Standard_NotImplemented: Blend_RstRstFunction::GetMinimalDistance");
    }
}
