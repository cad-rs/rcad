//! IntCurveSurface_HInter.hxx/.cxx (L1-582) — the curve × surface
//! intersection entry class over adaptors.  The OCCT body is a thin shell
//! delegating to IntCurveSurface_InterImpl (Inter.pxx) and
//! IntCurveSurface_InterUtils (InterUtils.pxx) through lambdas; here the
//! lambdas are the [`HInterHost`] implementation at the bottom of this
//! file.

use rcad_kernel::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};
use rcad_kernel::math::bnd::BoundSortBox;

use crate::geomalgo::int_imp::int_cs::IntCS;
use crate::geomalgo::int_imp::zer_cs_par_func::ZerCSParFunc;
use crate::geomalgo::int_patch::int_conic_quad::IntConicQuad;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};

use super::inter_impl;
use super::inter_utils::{compute_append_point, do_new_bounds, do_surface};
use super::quad_curv_exact::TheQuadCurvExactHInter;
use super::{
    HCurveTool, HInterHost, HSurfaceTool, Intersection, IntersectionPoint, SurfaceType,
};

/// OCCT IntCurveSurface_HInter (HInter.hxx L47-208) — extends
/// IntCurveSurface_Intersection (the `base` subobject).
#[derive(Debug, Clone, Default)]
pub struct HInter {
    pub base: Intersection,
}

impl HInter {
    /// OCCT IntCurveSurface_HInter() (cxx L56).
    pub fn new() -> Self {
        HInter {
            base: Intersection::new(),
        }
    }

    /// OCCT Perform(Curve, Surface) (cxx L106-116) — compute the
    /// intersection decomposing the surface by C2 intervals.
    pub fn perform<C: ?Sized, CT, S, ST>(&mut self, curve: &C, surface: &S)
    where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform::<C, CT, S, ST, Self>(curve, surface, self);
    }

    /// OCCT Perform(Curve, Surface, U0, V0, U1, V1) (cxx L120-150) —
    /// perform with the given UV bounds.
    #[allow(clippy::too_many_arguments)]
    pub fn perform_bounds<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        surface: &S,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_bounds::<C, CT, S, ST, Self>(curve, surface, u1, v1, u2, v2, self);
    }

    /// OCCT Perform(Curve, Polygon, Surface) (cxx L154-168).
    pub fn perform_polygon<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        polygon: &ThePolygonOfHInter,
        surface: &S,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_polygon::<C, CT, S, ST, Self>(curve, polygon, surface, self);
    }

    /// OCCT Perform(Curve, Surface, Polyhedron) (cxx L172-190).
    pub fn perform_polyhedron<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        surface: &S,
        polyhedron: &ThePolyhedronOfHInter,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_polyhedron::<C, CT, S, ST, Self>(curve, surface, polyhedron, self);
    }

    /// OCCT Perform(Curve, Polygon, Surface, Polyhedron) (cxx L194-219).
    pub fn perform_polygon_polyhedron<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        polygon: &ThePolygonOfHInter,
        surface: &S,
        polyhedron: &ThePolyhedronOfHInter,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_polygon_polyhedron::<C, CT, S, ST, Self>(
            curve, polygon, surface, polyhedron, self,
        );
    }

    /// OCCT Perform(Curve, Polygon, Surface, Polyhedron, BndBSB)
    /// (cxx L223-251).
    pub fn perform_polygon_polyhedron_bsb<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        polygon: &ThePolygonOfHInter,
        surface: &S,
        polyhedron: &ThePolyhedronOfHInter,
        bnd_bsb: &mut BoundSortBox,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_polygon_polyhedron_bsb::<C, CT, S, ST, Self>(
            curve, polygon, surface, polyhedron, bnd_bsb, self,
        );
    }

    /// OCCT InternalPerform(Curve, Polygon, Surface, Polyhedron, U0, V0, U1,
    /// V1, BSB) (cxx L255-284).
    #[allow(clippy::too_many_arguments)]
    pub fn internal_perform_bsb<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        polygon: &ThePolygonOfHInter,
        surface: &S,
        polyhedron: &ThePolyhedronOfHInter,
        u0: f64,
        v0: f64,
        u1: f64,
        v1: f64,
        bsb: &mut BoundSortBox,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::internal_perform_bsb::<C, CT, S, ST, Self>(
            curve, polygon, surface, polyhedron, u0, v0, u1, v1, bsb, self,
        );
    }

    /// OCCT InternalPerform(Curve, Polygon, Surface, Polyhedron, U0, V0, U1,
    /// V1) (cxx L288-315).
    #[allow(clippy::too_many_arguments)]
    pub fn internal_perform<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        polygon: &ThePolygonOfHInter,
        surface: &S,
        polyhedron: &ThePolyhedronOfHInter,
        u0: f64,
        v0: f64,
        u1: f64,
        v1: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::internal_perform::<C, CT, S, ST, Self>(
            curve, polygon, surface, polyhedron, u0, v0, u1, v1, self,
        );
    }

    /// OCCT InternalPerformCurveQuadric(Curve, Surface) (cxx L319-330).
    pub fn internal_perform_curve_quadric<C: ?Sized, CT, S, ST>(&mut self, curve: &C, surface: &S)
    where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::internal_perform_curve_quadric::<C, CT, S, ST, TheQuadCurvExactHInter<S, ST, C, CT>, Self>(
            curve, surface, self,
        );
    }

    /// OCCT InternalPerform(Curve, Polygon, Surface, U1, V1, U2, V2)
    /// (cxx L334-365) — the 7-arg version; instantiates
    /// InternalPerformPolygonBounds.
    #[allow(clippy::too_many_arguments)]
    pub fn internal_perform_polygon_bounds<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        polygon: &ThePolygonOfHInter,
        surface: &S,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::internal_perform_polygon_bounds::<C, CT, S, ST, TheQuadCurvExactHInter<S, ST, C, CT>, Self>(
            curve, polygon, surface, u1, v1, u2, v2, self,
        );
    }

    /// OCCT PerformConicSurf(Line, Curve, Surface, U1, V1, U2, V2)
    /// (cxx L369-402).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_conic_line<C: ?Sized, CT, S, ST>(
        &mut self,
        line: &Line3,
        curve: &C,
        surface: &S,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_conic_surf_line::<C, CT, S, ST, Self>(line, curve, surface, u1, v1, u2, v2, self);
    }

    /// OCCT PerformConicSurf(Circle, ...) (cxx L406-433).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_conic_circle<C: ?Sized, CT, S, ST>(
        &mut self,
        circle: &Circle3,
        curve: &C,
        surface: &S,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_conic_surf_circle::<C, CT, S, ST, Self>(
            circle, curve, surface, u1, v1, u2, v2, self,
        );
    }

    /// OCCT PerformConicSurf(Ellipse, ...) (cxx L437-464).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_conic_ellipse<C: ?Sized, CT, S, ST>(
        &mut self,
        ellipse: &Ellipse3,
        curve: &C,
        surface: &S,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_conic_surf_ellipse::<C, CT, S, ST, Self>(
            ellipse, curve, surface, u1, v1, u2, v2, self,
        );
    }

    /// OCCT PerformConicSurf(Parab, ...) (cxx L468-500).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_conic_parabola<C: ?Sized, CT, S, ST>(
        &mut self,
        parab: &Parabola3,
        curve: &C,
        surface: &S,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_conic_surf_parabola::<C, CT, S, ST, Self>(
            parab, curve, surface, u1, v1, u2, v2, self,
        );
    }

    /// OCCT PerformConicSurf(Hypr, ...) (cxx L504-536).
    #[allow(clippy::too_many_arguments)]
    pub fn perform_conic_hyperbola<C: ?Sized, CT, S, ST>(
        &mut self,
        hyper: &Hyperbola3,
        curve: &C,
        surface: &S,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::perform_conic_surf_hyperbola::<C, CT, S, ST, Self>(
            hyper, curve, surface, u1, v1, u2, v2, self,
        );
    }

    /// OCCT AppendIntAna(Curve, Surface, IntAna) (cxx L540-550).
    pub fn append_int_ana<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        surface: &S,
        intana_conic_quad: &IntConicQuad,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        inter_impl::append_int_ana::<C, CT, S, ST, Self>(curve, surface, intana_conic_quad, self);
    }

    /// OCCT AppendPoint(Curve, w, Surface, u, v) (cxx L554-571).
    pub fn append_point<C: ?Sized, CT, S, ST>(
        &mut self,
        curve: &C,
        lw: f64,
        surface: &S,
        su: f64,
        sv: f64,
    ) where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        if let Some(a_point) = compute_append_point::<C, CT, S, ST>(curve, lw, surface, su, sv) {
            self.base.append(&a_point);
        }
    }

    /// OCCT AppendSegment(Curve, u0, u1, Surface) (cxx L575-581) — not
    /// implemented in OCCT either.
    pub fn append_segment<C: ?Sized, CT, S, ST>(&mut self, _u0: f64, _u1: f64)
    where
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    {
        // Not implemented
    }

    /// OCCT DoSurface(surface, u0, u1, v0, v1, pntsOnSurface, boxSurface,
    /// gap) (cxx L60-77) — the 50x50 sampling grid (row-major here).
    pub(crate) fn do_surface<S, ST: HSurfaceTool<Surface = S>>(
        surface: &S,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
        pnts_on_surface: &mut Vec<glam::DVec3>,
        box_surface: &mut rcad_kernel::math::bnd::BndBox,
        gap: &mut f64,
    ) {
        do_surface::<S, ST>(surface, u0, u1, v0, v1, pnts_on_surface, box_surface, gap);
    }

    /// OCCT DoNewBounds(surface, u0, u1, v0, v1, pntsOnSurface, X, Y, Z,
    /// Bounds) (cxx L81-102).
    pub(crate) fn do_new_bounds<S, ST: HSurfaceTool<Surface = S>>(
        surface: &S,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
        pnts_on_surface: &[glam::DVec3],
        x: &[f64; 3],
        y: &[f64; 3],
        z: &[f64; 3],
        bounds: &mut [f64; 4],
    ) {
        do_new_bounds::<S, ST>(surface, u0, u1, v0, v1, pnts_on_surface, x, y, z, bounds);
    }
}

impl<
        C: ?Sized,
        CT: HCurveTool<Curve = C> + crate::geomalgo::int_imp::CurveTool3d<Curve = C>,
        S,
        ST: HSurfaceTool<Surface = S> + crate::geomalgo::int_imp::PSurfaceTool<Surface = S>,
    > HInterHost<C, CT, S, ST> for HInter
{
    fn done_flag(&mut self) -> &mut bool {
        &mut self.base.done
    }

    fn is_parallel_flag(&mut self) -> &mut bool {
        &mut self.base.my_is_parallel
    }

    fn reset_fields(&mut self) {
        self.base.reset_fields();
    }

    fn append(&mut self, pt: &IntersectionPoint) {
        self.base.append(pt);
    }

    fn perform_bounds(&mut self, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64) {
        HInter::perform_bounds::<C, CT, S, ST>(self, c, s, u1, v1, u2, v2);
    }

    fn perform_conic_line(&mut self, line: &Line3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64) {
        HInter::perform_conic_line::<C, CT, S, ST>(self, line, c, s, u1, v1, u2, v2);
    }

    fn perform_conic_circle(&mut self, circle: &Circle3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64) {
        HInter::perform_conic_circle::<C, CT, S, ST>(self, circle, c, s, u1, v1, u2, v2);
    }

    fn perform_conic_ellipse(&mut self, ellipse: &Ellipse3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64) {
        HInter::perform_conic_ellipse::<C, CT, S, ST>(self, ellipse, c, s, u1, v1, u2, v2);
    }

    fn perform_conic_parabola(&mut self, parab: &Parabola3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64) {
        HInter::perform_conic_parabola::<C, CT, S, ST>(self, parab, c, s, u1, v1, u2, v2);
    }

    fn perform_conic_hyperbola(&mut self, hyper: &Hyperbola3, c: &C, s: &S, u1: f64, v1: f64, u2: f64, v2: f64) {
        HInter::perform_conic_hyperbola::<C, CT, S, ST>(self, hyper, c, s, u1, v1, u2, v2);
    }

    fn perform_polygon_polyhedron(
        &mut self,
        c: &C,
        p: &ThePolygonOfHInter,
        s: &S,
        ph: &ThePolyhedronOfHInter,
    ) {
        HInter::perform_polygon_polyhedron::<C, CT, S, ST>(self, c, p, s, ph);
    }

    fn internal_perform(
        &mut self,
        c: &C,
        p: &ThePolygonOfHInter,
        s: &S,
        ph: &ThePolyhedronOfHInter,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        HInter::internal_perform::<C, CT, S, ST>(self, c, p, s, ph, u1, v1, u2, v2);
    }

    fn internal_perform_bsb(
        &mut self,
        c: &C,
        p: &ThePolygonOfHInter,
        s: &S,
        ph: &ThePolyhedronOfHInter,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        bsb: &mut BoundSortBox,
    ) {
        HInter::internal_perform_bsb::<C, CT, S, ST>(self, c, p, s, ph, u1, v1, u2, v2, bsb);
    }

    fn internal_perform_bounds(&mut self, c: &C, p: &ThePolygonOfHInter, s: &S, u1: f64, v1: f64, u2: f64, v2: f64) {
        HInter::internal_perform_polygon_bounds::<C, CT, S, ST>(self, c, p, s, u1, v1, u2, v2);
    }

    fn internal_perform_curve_quadric(&mut self, c: &C, s: &S) {
        HInter::internal_perform_curve_quadric::<C, CT, S, ST>(self, c, s);
    }

    fn append_int_ana(&mut self, c: &C, s: &S, ana: &IntConicQuad) {
        HInter::append_int_ana::<C, CT, S, ST>(self, c, s, ana);
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::geomalgo::int_imp::{CurveTool3d, PSurfaceTool};
    use crate::geomalgo::int_curve_surface::{HCurveTool, HSurfaceTool, TransitionOnCurve};
    use glam::DVec3;
    use rcad_kernel::base::proj_lib::CurveType;
    use rcad_kernel::geom::{
        BezierCurve3, Circle3, ConicalSurface, CylindricalSurface, CurveEval, Ellipse3,
        Hyperbola3, Line3, Parabola3, Plane, SphericalSurface, ToroidalSurface,
    };
    use rcad_kernel::math::GeomAbsShape;

    // ======================================================================
    // Test adaptors — the OCCT adaptor equivalents used by the anchors:
    // a line, a quadratic Bezier, the plane z = 0 (S(u,v) = (u,v,0)) and a
    // cylinder (axis Z, u = atan2(y, x), v = z).
    // ======================================================================

    /// The line C(w) = origin + w * direction over [first, last].
    #[derive(Clone, Copy)]
    struct TestLineCurve {
        line: Line3,
        first: f64,
        last: f64,
    }

    impl TestLineCurve {
        fn point_at(&self, u: f64) -> DVec3 {
            self.line.origin + u * self.line.direction
        }
        fn derivative_at(&self, _u: f64) -> DVec3 {
            self.line.direction
        }
        fn derivative2_at(&self, _u: f64) -> DVec3 {
            DVec3::ZERO
        }
    }

    /// The quadratic Bezier over [0, 1].
    #[derive(Clone)]
    struct TestBezierCurve {
        bez: BezierCurve3,
    }

    impl TestBezierCurve {
        fn point_at(&self, u: f64) -> DVec3 {
            self.bez.point_at(u)
        }
        fn derivative_at(&self, u: f64) -> DVec3 {
            self.bez.derivative_at(u)
        }
        fn derivative2_at(&self, _u: f64) -> DVec3 {
            DVec3::ZERO
        }
        fn bez(&self) -> &BezierCurve3 {
            &self.bez
        }
    }

    /// The plane z = 0 with S(u, v) = (u, v, 0) over the stored window
    /// (the OCCT adaptor carries the parameter bounds so UTrim/VTrim can
    /// narrow them).
    #[derive(Clone)]
    struct TestPlaneSurf {
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    }

    impl TestPlaneSurf {
        fn point_at(&self, u: f64, v: f64) -> DVec3 {
            DVec3::new(u, v, 0.0)
        }
        fn d1_at(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
            (self.point_at(u, v), DVec3::X, DVec3::Y)
        }
    }

    /// The cylinder x^2 + y^2 = R^2 with v = z over the stored window.
    #[derive(Clone)]
    struct TestCylSurf {
        radius: f64,
        u0: f64,
        u1: f64,
        v0: f64,
        v1: f64,
    }

    impl TestCylSurf {
        fn point_at(&self, u: f64, v: f64) -> DVec3 {
            DVec3::new(self.radius * u.cos(), self.radius * u.sin(), v)
        }
        fn d1_at(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
            let r = self.radius;
            (
                self.point_at(u, v),
                DVec3::new(-r * u.sin(), r * u.cos(), 0.0),
                DVec3::Z,
            )
        }
    }

    // ------------------------------------------------------------------
    // HCurveTool + CurveTool3d over the test curves.  The unused
    // type-specific accessors panic like the OCCT
    // Standard_NoSuchObject raises.
    // ------------------------------------------------------------------

    /// Test tool over [`TestLineCurve`].
    struct LineCT;

    impl HCurveTool for LineCT {
        type Curve = TestLineCurve;

        fn first_parameter(c: &TestLineCurve) -> f64 {
            c.first
        }
        fn last_parameter(c: &TestLineCurve) -> f64 {
            c.last
        }
        fn continuity(_c: &TestLineCurve) -> GeomAbsShape {
            GeomAbsShape::C2
        }
        fn nb_intervals(_c: &TestLineCurve, _s: GeomAbsShape) -> usize {
            1
        }
        fn intervals(c: &TestLineCurve, t: &mut [f64], _s: GeomAbsShape) {
            t[0] = c.first;
            t[1] = c.last;
        }
        fn is_closed(_c: &TestLineCurve) -> bool {
            false
        }
        fn is_periodic(_c: &TestLineCurve) -> bool {
            false
        }
        fn period(_c: &TestLineCurve) -> f64 {
            panic!("Standard_NoSuchObject");
        }
        fn value(c: &TestLineCurve, u: f64) -> DVec3 {
            c.point_at(u)
        }
        fn d0(c: &TestLineCurve, u: f64) -> DVec3 {
            c.point_at(u)
        }
        fn d1(c: &TestLineCurve, u: f64) -> (DVec3, DVec3) {
            (c.point_at(u), c.derivative_at(u))
        }
        fn d2(c: &TestLineCurve, u: f64) -> (DVec3, DVec3, DVec3) {
            (c.point_at(u), c.derivative_at(u), c.derivative2_at(u))
        }
        fn d3(_c: &TestLineCurve, _u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
            panic!("Standard_NoSuchObject");
        }
        fn dn(_c: &TestLineCurve, _u: f64, _n: usize) -> DVec3 {
            panic!("Standard_NoSuchObject");
        }
        fn resolution(_c: &TestLineCurve, r3d: f64) -> f64 {
            r3d
        }
        fn get_type(_c: &TestLineCurve) -> CurveType {
            CurveType::Line
        }
        fn line(c: &TestLineCurve) -> Line3 {
            c.line
        }
        fn circle(_c: &TestLineCurve) -> Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(_c: &TestLineCurve) -> Ellipse3 {
            panic!("Standard_NoSuchObject");
        }
        fn hyperbola(_c: &TestLineCurve) -> Hyperbola3 {
            panic!("Standard_NoSuchObject");
        }
        fn parabola(_c: &TestLineCurve) -> Parabola3 {
            panic!("Standard_NoSuchObject");
        }
        fn bezier(_c: &TestLineCurve) -> &BezierCurve3 {
            panic!("Standard_NoSuchObject");
        }
        fn bspline(_c: &TestLineCurve) -> &rcad_kernel::geom::BSplineCurve3 {
            panic!("Standard_NoSuchObject");
        }
    }

    impl CurveTool3d for LineCT {
        type Curve = TestLineCurve;
        fn value(c: &TestLineCurve, u: f64) -> DVec3 {
            c.point_at(u)
        }
        fn d1(c: &TestLineCurve, u: f64) -> (DVec3, DVec3) {
            (c.point_at(u), c.derivative_at(u))
        }
        fn first_parameter(c: &TestLineCurve) -> f64 {
            <LineCT as HCurveTool>::first_parameter(c)
        }
        fn last_parameter(c: &TestLineCurve) -> f64 {
            <LineCT as HCurveTool>::last_parameter(c)
        }
        fn resolution(_c: &TestLineCurve, r3d: f64) -> f64 {
            r3d
        }
    }

    /// Test tool over [`TestBezierCurve`].
    struct BezierCT;

    impl HCurveTool for BezierCT {
        type Curve = TestBezierCurve;

        fn first_parameter(_c: &TestBezierCurve) -> f64 {
            0.0
        }
        fn last_parameter(_c: &TestBezierCurve) -> f64 {
            1.0
        }
        fn continuity(_c: &TestBezierCurve) -> GeomAbsShape {
            GeomAbsShape::C2
        }
        fn nb_intervals(_c: &TestBezierCurve, _s: GeomAbsShape) -> usize {
            1
        }
        fn intervals(_c: &TestBezierCurve, t: &mut [f64], _s: GeomAbsShape) {
            t[0] = 0.0;
            t[1] = 1.0;
        }
        fn is_closed(_c: &TestBezierCurve) -> bool {
            false
        }
        fn is_periodic(_c: &TestBezierCurve) -> bool {
            false
        }
        fn period(_c: &TestBezierCurve) -> f64 {
            panic!("Standard_NoSuchObject");
        }
        fn value(c: &TestBezierCurve, u: f64) -> DVec3 {
            c.point_at(u)
        }
        fn d0(c: &TestBezierCurve, u: f64) -> DVec3 {
            c.point_at(u)
        }
        fn d1(c: &TestBezierCurve, u: f64) -> (DVec3, DVec3) {
            (c.point_at(u), c.derivative_at(u))
        }
        fn d2(c: &TestBezierCurve, u: f64) -> (DVec3, DVec3, DVec3) {
            (c.point_at(u), c.derivative_at(u), c.derivative2_at(u))
        }
        fn d3(_c: &TestBezierCurve, _u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
            panic!("Standard_NoSuchObject");
        }
        fn dn(_c: &TestBezierCurve, _u: f64, _n: usize) -> DVec3 {
            panic!("Standard_NoSuchObject");
        }
        fn resolution(_c: &TestBezierCurve, r3d: f64) -> f64 {
            r3d
        }
        fn get_type(_c: &TestBezierCurve) -> CurveType {
            CurveType::Bezier
        }
        fn line(_c: &TestBezierCurve) -> Line3 {
            panic!("Standard_NoSuchObject");
        }
        fn circle(_c: &TestBezierCurve) -> Circle3 {
            panic!("Standard_NoSuchObject");
        }
        fn ellipse(_c: &TestBezierCurve) -> Ellipse3 {
            panic!("Standard_NoSuchObject");
        }
        fn hyperbola(_c: &TestBezierCurve) -> Hyperbola3 {
            panic!("Standard_NoSuchObject");
        }
        fn parabola(_c: &TestBezierCurve) -> Parabola3 {
            panic!("Standard_NoSuchObject");
        }
        fn bezier(c: &TestBezierCurve) -> &BezierCurve3 {
            c.bez()
        }
        fn bspline(_c: &TestBezierCurve) -> &rcad_kernel::geom::BSplineCurve3 {
            panic!("Standard_NoSuchObject");
        }
    }

    impl CurveTool3d for BezierCT {
        type Curve = TestBezierCurve;
        fn value(c: &TestBezierCurve, u: f64) -> DVec3 {
            c.point_at(u)
        }
        fn d1(c: &TestBezierCurve, u: f64) -> (DVec3, DVec3) {
            (c.point_at(u), c.derivative_at(u))
        }
        fn first_parameter(c: &TestBezierCurve) -> f64 {
            <BezierCT as HCurveTool>::first_parameter(c)
        }
        fn last_parameter(c: &TestBezierCurve) -> f64 {
            <BezierCT as HCurveTool>::last_parameter(c)
        }
        fn resolution(_c: &TestBezierCurve, r3d: f64) -> f64 {
            r3d
        }
    }

    // ------------------------------------------------------------------
    // HSurfaceTool + PSurfaceTool over the test surfaces.
    // ------------------------------------------------------------------

    macro_rules! impl_surface_tools {
        ($marker:ident, $surf:ty, $stype:expr,
         $is_u_closed:expr, $is_u_periodic:expr, $u_period:expr, $is_v_closed:expr,
         $is_v_periodic:expr, $v_period:expr, $plane:expr, $cylinder:expr) => {
            struct $marker;

            impl HSurfaceTool for $marker {
                type Surface = $surf;
                type BasisCurve = TestBasisCurve;
                type BasisSurface = TestBasisSurface;

                fn first_u_parameter(s: &$surf) -> f64 {
                    s.u0
                }
                fn first_v_parameter(s: &$surf) -> f64 {
                    s.v0
                }
                fn last_u_parameter(s: &$surf) -> f64 {
                    s.u1
                }
                fn last_v_parameter(s: &$surf) -> f64 {
                    s.v1
                }
                fn nb_u_intervals(_s: &$surf, _sh: GeomAbsShape) -> usize {
                    1
                }
                fn nb_v_intervals(_s: &$surf, _sh: GeomAbsShape) -> usize {
                    1
                }
                fn u_intervals(s: &$surf, tab: &mut [f64], _sh: GeomAbsShape) {
                    tab[0] = s.u0;
                    tab[1] = s.u1;
                }
                fn v_intervals(s: &$surf, tab: &mut [f64], _sh: GeomAbsShape) {
                    tab[0] = s.v0;
                    tab[1] = s.v1;
                }
                fn is_u_closed(_s: &$surf) -> bool {
                    $is_u_closed
                }
                fn is_v_closed(_s: &$surf) -> bool {
                    $is_v_closed
                }
                fn is_u_periodic(_s: &$surf) -> bool {
                    $is_u_periodic
                }
                fn u_period(_s: &$surf) -> f64 {
                    $u_period
                }
                fn is_v_periodic(_s: &$surf) -> bool {
                    $is_v_periodic
                }
                fn v_period(_s: &$surf) -> f64 {
                    $v_period
                }
                fn value(s: &$surf, u: f64, v: f64) -> DVec3 {
                    s.point_at(u, v)
                }
                fn d0(s: &$surf, u: f64, v: f64) -> DVec3 {
                    s.point_at(u, v)
                }
                fn d1(s: &$surf, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
                    s.d1_at(u, v)
                }
                fn d2(
                    s: &$surf,
                    u: f64,
                    v: f64,
                ) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3) {
                    let (p, du, dv) = s.d1_at(u, v);
                    (p, du, dv, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO)
                }
                fn dn(_s: &$surf, _u: f64, _v: f64, _nu: usize, _nv: usize) -> DVec3 {
                    panic!("Standard_NoSuchObject");
                }
                fn u_resolution(_s: &$surf, r3d: f64) -> f64 {
                    r3d
                }
                fn v_resolution(_s: &$surf, r3d: f64) -> f64 {
                    r3d
                }
                fn get_type(_s: &$surf) -> SurfaceType {
                    $stype
                }
                fn plane(s: &$surf) -> Plane {
                    let _ = s;
                    $plane
                }
                fn cylinder(s: &$surf) -> CylindricalSurface {
                    let _ = s;
                    $cylinder
                }
                fn cone(_s: &$surf) -> ConicalSurface {
                    panic!("Standard_NoSuchObject");
                }
                fn torus(_s: &$surf) -> ToroidalSurface {
                    panic!("Standard_NoSuchObject");
                }
                fn sphere(_s: &$surf) -> SphericalSurface {
                    panic!("Standard_NoSuchObject");
                }
                fn axe_of_revolution(_s: &$surf) -> (DVec3, DVec3) {
                    panic!("Standard_NoSuchObject");
                }
                fn direction(_s: &$surf) -> DVec3 {
                    panic!("Standard_NoSuchObject");
                }
                fn basis_curve(_s: &$surf) -> TestBasisCurve {
                    panic!("Standard_NoSuchObject");
                }
                fn basis_surface(_s: &$surf) -> TestBasisSurface {
                    panic!("Standard_NoSuchObject");
                }
                fn offset_value(_s: &$surf) -> f64 {
                    panic!("Standard_NoSuchObject");
                }
                fn u_trim(s: &$surf, first: f64, last: f64, _tol: f64) -> $surf {
                    // OCCT GeomAdaptor_Surface::UTrim — the same surface,
                    // narrowed U window.
                    let mut t = s.clone();
                    t.u0 = first;
                    t.u1 = last;
                    t
                }
                fn v_trim(s: &$surf, first: f64, last: f64, _tol: f64) -> $surf {
                    let mut t = s.clone();
                    t.v0 = first;
                    t.v1 = last;
                    t
                }
                fn bezier(_s: &$surf) -> &rcad_kernel::geom::BezierSurface {
                    panic!("Standard_NoSuchObject");
                }
                fn bspline(_s: &$surf) -> &rcad_kernel::geom::BSplineSurface {
                    panic!("Standard_NoSuchObject");
                }
                fn nb_u_poles(_s: &$surf) -> usize {
                    panic!("Standard_NoSuchObject");
                }
                fn nb_v_poles(_s: &$surf) -> usize {
                    panic!("Standard_NoSuchObject");
                }
                fn nb_u_knots(_s: &$surf) -> usize {
                    panic!("Standard_NoSuchObject");
                }
                fn nb_v_knots(_s: &$surf) -> usize {
                    panic!("Standard_NoSuchObject");
                }
                fn u_degree(_s: &$surf) -> usize {
                    panic!("Standard_NoSuchObject");
                }
                fn v_degree(_s: &$surf) -> usize {
                    panic!("Standard_NoSuchObject");
                }
            }

            impl PSurfaceTool for $marker {
                type Surface = $surf;
                fn value(s: &$surf, u: f64, v: f64) -> DVec3 {
                    s.point_at(u, v)
                }
                fn d1(s: &$surf, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
                    s.d1_at(u, v)
                }
                fn first_u_parameter(s: &$surf) -> f64 {
                    s.u0
                }
                fn last_u_parameter(s: &$surf) -> f64 {
                    s.u1
                }
                fn first_v_parameter(s: &$surf) -> f64 {
                    s.v0
                }
                fn last_v_parameter(s: &$surf) -> f64 {
                    s.v1
                }
                fn u_resolution(_s: &$surf, r3d: f64) -> f64 {
                    r3d
                }
                fn v_resolution(_s: &$surf, r3d: f64) -> f64 {
                    r3d
                }
            }
        };
    }

    impl_surface_tools!(
        PlaneST,
        TestPlaneSurf,
        SurfaceType::Plane,
        false,
        false,
        unreachable!(),
        false,
        false,
        unreachable!(),
        Plane {
            origin: DVec3::ZERO,
            normal: DVec3::Z,
            u_dir: DVec3::X,
            v_dir: DVec3::Y,
        },
        panic!("Standard_NoSuchObject")
    );

    impl_surface_tools!(
        CylST,
        TestCylSurf,
        SurfaceType::Cylinder,
        true,
        true,
        std::f64::consts::TAU,
        false,
        false,
        unreachable!(),
        panic!("Standard_NoSuchObject"),
        CylindricalSurface::new_with_ref_dir(DVec3::ZERO, DVec3::Z, 2.0, DVec3::X)
    );

    /// The OCCT Adaptor3d_Curve / Adaptor3d_Surface stubs — never reached by
    /// the anchors (only the EstLim* basis-surface paths use them).
    struct TestBasisCurve;
    impl crate::geomalgo::int_curve_surface::Adaptor3dCurveBasis for TestBasisCurve {
        fn get_type(&self) -> CurveType {
            panic!("Standard_NoSuchObject");
        }
        fn value(&self, _u: f64) -> DVec3 {
            panic!("Standard_NoSuchObject");
        }
        fn line(&self) -> Line3 {
            panic!("Standard_NoSuchObject");
        }
        fn parabola(&self) -> Parabola3 {
            panic!("Standard_NoSuchObject");
        }
        fn hyperbola(&self) -> Hyperbola3 {
            panic!("Standard_NoSuchObject");
        }
    }

    struct TestBasisSurface;
    impl crate::geomalgo::int_curve_surface::Adaptor3dSurfaceBasis<TestBasisCurve> for TestBasisSurface {
        fn get_type(&self) -> SurfaceType {
            panic!("Standard_NoSuchObject");
        }
        fn plane(&self) -> Plane {
            panic!("Standard_NoSuchObject");
        }
        fn cylinder(&self) -> CylindricalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn cone(&self) -> ConicalSurface {
            panic!("Standard_NoSuchObject");
        }
        fn direction(&self) -> DVec3 {
            panic!("Standard_NoSuchObject");
        }
        fn basis_curve(&self) -> TestBasisCurve {
            panic!("Standard_NoSuchObject");
        }
    }

    /// OCCT anchor: a line piercing the plane through the (Lin, Pln, Tolang)
    /// IntAna path.  C(w) = (0, 0.25, 1) + w * (1, 0, -1)/sqrt(2) over
    /// [-5, 5] hits z = 0 at w = sqrt(2), point (1, 0.25, 0), u = 1,
    /// v = 0.25, transition In (the gp_Lin parameter is the arc length
    /// along the unit direction).
    #[test]
    fn hinter_line_plane_intana_anchor() {
        let curve = TestLineCurve {
            line: Line3 {
                origin: DVec3::new(0.0, 0.25, 1.0),
                direction: DVec3::new(1.0, 0.0, -1.0).normalize(),
            },
            first: -5.0,
            last: 5.0,
        };
        let surface = TestPlaneSurf {
            u0: -2.0,
            u1: 2.0,
            v0: -2.0,
            v1: 2.0,
        };

        let mut h = HInter::new();
        h.perform::<TestLineCurve, LineCT, TestPlaneSurf, PlaneST>(&curve, &surface);

        assert!(h.base.is_done());
        assert!(!h.base.is_parallel());
        assert_eq!(h.base.nb_points(), 1, "points: {:?}", h.base.lpnt);
        let pt = h.base.point(1);
        assert!((pt.w() - std::f64::consts::SQRT_2).abs() < 1e-9, "w={}", pt.w());
        assert!((pt.u() - 1.0).abs() < 1e-9, "u={}", pt.u());
        assert!((pt.v() - 0.25).abs() < 1e-9, "v={}", pt.v());
        assert!(pt.pnt().distance(DVec3::new(1.0, 0.25, 0.0)) < 1e-9);
        assert_eq!(pt.transition(), TransitionOnCurve::In);
    }

    /// OCCT anchor: a line through the cylinder via the implicit
    /// (Lin, IntAna_Quadric) path — two crossings, In then Out.
    #[test]
    fn hinter_line_cylinder_intana_anchor() {
        let curve = TestLineCurve {
            line: Line3 {
                origin: DVec3::new(-10.0, 0.0, 0.5),
                direction: DVec3::X,
            },
            first: -10.0,
            last: 12.0,
        };
        let surface = TestCylSurf {
            radius: 2.0,
            u0: 0.0,
            u1: std::f64::consts::TAU,
            v0: -1.0,
            v1: 1.0,
        };

        let mut h = HInter::new();
        h.perform::<TestLineCurve, LineCT, TestCylSurf, CylST>(&curve, &surface);

        assert!(h.base.is_done());
        assert_eq!(h.base.nb_points(), 2, "points: {:?}", h.base.lpnt);
        // The IntAna path appends in the polynomial root order (not sorted);
        // match the two crossings order-agnostically: entry w = 8 (u = PI,
        // In) and exit w = 12 (u = 0, Out), from the location (-10, 0, 0.5).
        let ws = [h.base.point(1).w(), h.base.point(2).w()];
        let entry = if ws[0] < ws[1] { h.base.point(1) } else { h.base.point(2) };
        let exit = if ws[0] < ws[1] { h.base.point(2) } else { h.base.point(1) };
        assert!((entry.w() - 8.0).abs() < 1e-9, "entry w={}", entry.w());
        assert!((exit.w() - 12.0).abs() < 1e-9, "exit w={}", exit.w());
        assert!((entry.u() - std::f64::consts::PI).abs() < 1e-9, "entry u={}", entry.u());
        assert!(exit.u().abs() < 1e-9, "exit u={}", exit.u());
        assert!((entry.v() - 0.5).abs() < 1e-9);
        assert_eq!(entry.transition(), TransitionOnCurve::In);
        assert_eq!(exit.transition(), TransitionOnCurve::Out);
    }

    /// The quadratic Bezier (-1,0,1) -> (0,0,-2) -> (1,0,1): z(w) =
    /// 6w^2 - 6w + 1 crosses z = 0 at w = (3 +- sqrt(3))/6, x = 2w - 1 =
    /// -+ sqrt(3)/3.  Exact quadric intersection (QuadCurvExactHInter
    /// route: non-conic curve x quadric surface).
    #[test]
    fn hinter_bezier_plane_quadcurv_anchor() {
        let bez = BezierCurve3 {
            control_points: vec![
                DVec3::new(-1.0, 0.0, 1.0),
                DVec3::new(0.0, 0.0, -2.0),
                DVec3::new(1.0, 0.0, 1.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let curve = TestBezierCurve { bez };
        let surface = TestPlaneSurf {
            u0: -2.0,
            u1: 2.0,
            v0: -2.0,
            v1: 2.0,
        };

        let mut h = HInter::new();
        h.perform::<TestBezierCurve, BezierCT, TestPlaneSurf, PlaneST>(&curve, &surface);

        assert!(h.base.is_done());
        assert_eq!(h.base.nb_points(), 2, "points: {:?}", h.base.lpnt);
        let s3 = 3.0f64.sqrt();
        let w1 = (3.0 - s3) / 6.0;
        let w2 = (3.0 + s3) / 6.0;
        let p1 = h.base.point(1);
        let p2 = h.base.point(2);
        assert!((p1.w() - w1).abs() < 1e-6, "w1={} vs {}", p1.w(), w1);
        assert!((p2.w() - w2).abs() < 1e-6, "w2={} vs {}", p2.w(), w2);
        assert!((p1.u() + s3 / 3.0).abs() < 1e-6, "u1={}", p1.u());
        assert!((p2.u() - s3 / 3.0).abs() < 1e-6, "u2={}", p2.u());
        assert!(p1.v().abs() < 1e-9 && p2.v().abs() < 1e-9);
        assert_eq!(p1.transition(), TransitionOnCurve::In);
        assert_eq!(p2.transition(), TransitionOnCurve::Out);
    }

    /// The same Bezier x plane through the polygon/polyhedron interference
    /// and the IntCS exact refinement (PerformPolygonPolyhedron ->
    /// InternalPerform -> Interference + ZerCSParFunc + IntCS).
    #[test]
    fn hinter_bezier_plane_intcs_refine_anchor() {
        let bez = BezierCurve3 {
            control_points: vec![
                DVec3::new(-1.0, 0.0, 1.0),
                DVec3::new(0.0, 0.0, -2.0),
                DVec3::new(1.0, 0.0, 1.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let curve = TestBezierCurve { bez };
        let surface = TestPlaneSurf {
            u0: -2.0,
            u1: 2.0,
            v0: -2.0,
            v1: 2.0,
        };

        let polygon = crate::geomalgo::int_curv_surf::ThePolygonOfHInter::new_tool::<
            TestBezierCurve,
            BezierCT,
        >(&curve, 11);
        let polyhedron = crate::geomalgo::int_curv_surf::ThePolyhedronOfHInter::new_tool::<
            TestPlaneSurf,
            PlaneST,
        >(&surface, 6, 6, -2.0, -2.0, 2.0, 2.0);

        let mut h = HInter::new();
        h.perform_polygon_polyhedron::<TestBezierCurve, BezierCT, TestPlaneSurf, PlaneST>(
            &curve, &polygon, &surface, &polyhedron,
        );

        assert!(h.base.is_done());
        assert_eq!(h.base.nb_points(), 2, "points: {:?}", h.base.lpnt);
        let s3 = 3.0f64.sqrt();
        let w1 = (3.0 - s3) / 6.0;
        let w2 = (3.0 + s3) / 6.0;
        let p1 = h.base.point(1);
        let p2 = h.base.point(2);
        assert!((p1.w() - w1).abs() < 1e-5, "w1={} vs {}", p1.w(), w1);
        assert!((p2.w() - w2).abs() < 1e-5, "w2={} vs {}", p2.w(), w2);
        assert!((p1.u() + s3 / 3.0).abs() < 1e-5, "u1={}", p1.u());
        assert!((p2.u() - s3 / 3.0).abs() < 1e-5, "u2={}", p2.u());
    }
}
