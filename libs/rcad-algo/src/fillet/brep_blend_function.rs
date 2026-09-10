//! OCCT Blend_AppFunction (TKFillet/Blend) — 1:1 port of Blend_AppFunction.hxx
//! (L44-193) + Blend_AppFunction.cxx (L20-23); and Blend_Function — 1:1 port
//! of Blend_Function.hxx (L42-119) + Blend_Function.cxx (L19-56).
//!
//! Architecture mappings: `class Blend_AppFunction : public
//! math_FunctionSetWithDerivatives` maps to a Rust subtrait of
//! [`FunctionSetWithDerivatives`]; `class Blend_Function : public
//! Blend_AppFunction` maps to a subtrait of [`BlendAppFunction`].  The
//! default-implemented trait methods correspond to the OCCT methods that are
//! implemented (not deferred) in Blend_AppFunction.cxx / Blend_Function.cxx.
//! The three OCCT `Section` overloads are named `section`, `section_d1` and
//! `section_d2`; the two `GetTolerance` overloads are `get_tolerance` and
//! `get_approx_tolerance`.

use glam::{DVec2, DVec3};

use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;

use super::brep_blend_point::BlendPoint;

/// OCCT Blend_AppFunction — deferred class for a function used to compute a
/// blending surface between two surfaces, using a guide line
/// (Blend_AppFunction.hxx L44).
pub trait BlendAppFunction: FunctionSetWithDerivatives {
    /// OCCT Set(Param) — sets the value of the parameter along the guide
    /// line; this determines the plane in which the solution has to be found.
    fn set_param(&mut self, param: f64);

    /// OCCT Set(First, Last) — sets the bounds of the parametric interval on
    /// the guide line.
    fn set_interval(&mut self, first: f64, last: f64);

    /// OCCT GetTolerance(Tolerance, Tol) — the parametric tolerance for each
    /// variable; Tol is the tolerance used in 3d space.
    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64);

    /// OCCT GetBounds(InfBound, SupBound) — the lowest / greatest values
    /// allowed for each variable.
    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]);

    /// OCCT IsSolution(Sol, Tol) — true if Sol is a zero of the function.
    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool;

    /// OCCT GetMinimalDistance() — the minimal Distance between two
    /// extremities of calculated sections.
    fn get_minimal_distance(&self) -> f64;

    /// OCCT Pnt1() — the point on the first support.
    fn pnt1(&self) -> DVec3;

    /// OCCT Pnt2() — the point on the second support.
    fn pnt2(&self) -> DVec3;

    /// OCCT IsRational() — true if the section is rational.
    fn is_rational(&self) -> bool;

    /// OCCT GetSectionSize() — the length of the maximum section.
    fn get_section_size(&self) -> f64;

    /// OCCT GetMinimalWeight(Weigths) — the minimal value of the weight for
    /// each pole of all sections.
    fn get_minimal_weight(&self, weigths: &mut [f64]);

    /// OCCT NbIntervals(S) — the number of intervals for continuity S.
    fn nb_intervals(&self, s: GeomAbsShape) -> usize;

    /// OCCT Intervals(T, S) — the parameters bounding the intervals of
    /// continuity S.
    fn intervals(&self, t: &mut [f64], s: GeomAbsShape);

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d).
    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    );

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d) — the
    /// tolerances to reach in approximation.
    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    );

    /// OCCT Knots(TKnots).
    fn knots(&mut self, tknots: &mut [f64]);

    /// OCCT Mults(TMults).
    fn mults(&mut self, tmults: &mut [i32]);

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths) —
    /// used for the first and last section; returns true if the derivatives
    /// are computed, false otherwise.
    #[allow(clippy::too_many_arguments)]
    fn section_d1(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool;

    /// OCCT Section(P, Poles, Poles2d, Weigths).
    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    );

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weigths, DWeigths, D2Weigths) — used for the first and last section;
    /// returns true if the derivatives are computed, false otherwise.
    #[allow(clippy::too_many_arguments)]
    fn section_d2(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        d2_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        d2_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
        d2_weigths: &mut [f64],
    ) -> bool;

    /// OCCT Resolution(IC2d, Tol, TolU, TolV).
    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64);

    /// OCCT Parameter(P) (Blend_AppFunction.cxx L20-23) — the parameter of
    /// the point P, used to impose the parameters in the approximation.
    fn parameter(&self, p: &BlendPoint) -> f64 {
        p.parameter()
    }
}

/// OCCT Blend_Function — deferred class for a function used to compute a
/// blending surface between two surfaces, using a guide line
/// (Blend_Function.hxx L42).
pub trait BlendFunction: BlendAppFunction {
    /// OCCT NbVariables() (Blend_Function.cxx L19-22) — returns 4.
    fn nb_variables(&self) -> usize {
        4
    }

    /// OCCT Pnt1() (Blend_Function.cxx L24-27) — delegates to PointOnS1.
    fn pnt1(&self) -> DVec3 {
        self.point_on_s1()
    }

    /// OCCT Pnt2() (Blend_Function.cxx L29-32) — delegates to PointOnS2.
    fn pnt2(&self) -> DVec3 {
        self.point_on_s2()
    }

    /// OCCT PointOnS1() — the point on the first surface, at parameter
    /// Sol(1), Sol(2).
    fn point_on_s1(&self) -> DVec3;

    /// OCCT PointOnS2() — the point on the second surface, at parameter
    /// Sol(3), Sol(4).
    fn point_on_s2(&self) -> DVec3;

    /// OCCT IsTangencyPoint() — true when it is not possible to compute the
    /// tangent vectors at PointOnS1 and/or PointOnS2.
    fn is_tangency_point(&self) -> bool;

    /// OCCT TangentOnS1() — the tangent vector at PointOnS1, in 3d space.
    fn tangent_on_s1(&self) -> DVec3;

    /// OCCT Tangent2dOnS1() — the tangent vector at PointOnS1, in the
    /// parametric space of the first surface.
    fn tangent_2d_on_s1(&self) -> DVec2;

    /// OCCT TangentOnS2() — the tangent vector at PointOnS2, in 3d space.
    fn tangent_on_s2(&self) -> DVec3;

    /// OCCT Tangent2dOnS2() — the tangent vector at PointOnS2, in the
    /// parametric space of the second surface.
    fn tangent_2d_on_s2(&self) -> DVec2;

    /// OCCT Tangent(U1, V1, U2, V2, TgFirst, TgLast, NormFirst, NormLast) —
    /// the tangent vector at the section, at the beginning and the end of
    /// the section, and the normal (of the surfaces) at these points.
    #[allow(clippy::too_many_arguments)]
    fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut DVec3,
        tg_last: &mut DVec3,
        norm_first: &mut DVec3,
        norm_last: &mut DVec3,
    );

    /// OCCT TwistOnS1() (Blend_Function.cxx L34-37) — returns false.
    fn twist_on_s1(&self) -> bool {
        false
    }

    /// OCCT TwistOnS2() (Blend_Function.cxx L39-42) — returns false.
    fn twist_on_s2(&self) -> bool {
        false
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weigths, DWeigths, D2Weigths) (Blend_Function.cxx L44-56) — the
    /// default returns false.
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
}
