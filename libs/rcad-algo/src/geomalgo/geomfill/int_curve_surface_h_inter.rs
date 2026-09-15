//! OCCT IntCurveSurface_HInter over the production adaptor pair (cxx L1-582)
//! as consumed by GeomFill_SectionPlacement / GeomFill_GuideTrihedronPlan:
//! `TheCurve = occ::handle<Adaptor3d_Curve>` and `TheSurface =
//! occ::handle<GeomAdaptor_Surface>`.
//!
//! The rcad IntCurveSurface engine
//! ([`crate::geomalgo::int_curve_surface::HInter`], the generic Inter.pxx /
//! InterUtils.pxx translation) is driven here through the production
//! host-tool markers — the OCCT template statics:
//! - [`TheCurveTool`] — IntCurveSurface_TheHCurveTool over the erased
//!   `Adaptor3d_Curve` handle (the [`CurveHost`] view bundle).
//! - [`TheSurfaceTool`] — Adaptor3d_HSurfaceTool over
//!   [`rcad_kernel::base::proj_lib::GeomSurfaceAdaptor`].
//!
//! Architecture difference (the rcad split of `Adaptor3d_Curve`): the OCCT
//! class carries the evaluation virtuals AND the geometry accessors
//! (GetType / Line / Circle / ... / Bezier / BSpline) on one class; the rcad
//! kernel splits them over `Adaptor3dCurve` and the companion trait
//! `Adaptor3dCurveGeom`, and the kernel `Curve3` payload answers the
//! reference-returning Bezier/BSpline downcasts.  [`CurveHost`] bundles the
//! views; the absent halves answer exactly the OCCT base-class values
//! (GetType -> GeomAbs_OtherCurve, the analytic downcasts ->
//! Standard_NoSuchObject).

use std::marker::PhantomData;

use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor3dCurve, Adaptor3dSurface, GeomAbsSurfaceType as KernelSurfaceType,
};
use rcad_kernel::base::proj_lib::geom_adaptor_curve::{Adaptor3dCurveGeom, GeomCurveAdaptor};
use rcad_kernel::base::proj_lib::geom_adaptor_surface::GeomSurfaceAdaptor;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{
    BezierCurve3, BezierSurface, BSplineCurve3, BSplineSurface, Circle3, ConicalSurface,
    CylindricalSurface, Curve3, Ellipse3, Hyperbola3, Line3, Parabola3, Plane, Point3,
    SphericalSurface, Surface3, ToroidalSurface, Vec3,
};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::int_curve_surface::{
    Adaptor3dCurveBasis, Adaptor3dSurfaceBasis, HCurveTool, HInter, HSurfaceTool,
    IntersectionPoint, SurfaceType,
};
use crate::geomalgo::int_imp::{CurveTool3d, PSurfaceTool};

/// OCCT `occ::handle<Adaptor3d_Curve>` — the erased curve adaptor handle
/// consumed by the HInter assembly, bundled with the rcad split adaptor
/// views (see the module doc).
#[derive(Clone, Copy)]
pub struct CurveHost<'a> {
    /// The evaluation / interval / resolution half (`Adaptor3d_Curve`).
    pub curve: &'a dyn Adaptor3dCurve,
    /// The GetType / Line / Circle / ... half (`Adaptor3d_Curve` geometry
    /// accessors), when the concrete adaptor provides it
    /// (the GeomAdaptor_Curve encodings).
    pub geom: Option<&'a dyn Adaptor3dCurveGeom>,
    /// The kernel `Geom_Curve` payload, when the caller holds it — answers
    /// the Bezier / BSpline downcasts (the reference-returning OCCT queries).
    pub curve3: Option<&'a Curve3>,
}

impl<'a> CurveHost<'a> {
    /// The erased `handle<Adaptor3d_Curve>` without the geometry views: the
    /// OCCT base-class answers (GetType -> OtherCurve, the analytic
    /// downcasts -> Standard_NoSuchObject).
    pub fn new(curve: &'a dyn Adaptor3dCurve) -> Self {
        CurveHost {
            curve,
            geom: None,
            curve3: None,
        }
    }
}

/// OCCT IntCurveSurface_TheHCurveTool (TheHCurveTool.hxx L41-191 + .cxx
/// L34-311) — the static tool over the erased `Adaptor3d_Curve` handle; the
/// OCCT template parameter `TheCurveTool` of the HInter instantiation.
pub struct TheCurveTool<'a>(PhantomData<&'a ()>);

impl<'a> HCurveTool for TheCurveTool<'a> {
    type Curve = CurveHost<'a>;

    /// OCCT FirstParameter(C) (hxx L46-49).
    fn first_parameter(c: &CurveHost<'a>) -> f64 {
        c.curve.first_parameter()
    }

    /// OCCT LastParameter(C) (hxx L51-54).
    fn last_parameter(c: &CurveHost<'a>) -> f64 {
        c.curve.last_parameter()
    }

    /// OCCT Continuity(C) (hxx L56-59).
    fn continuity(c: &CurveHost<'a>) -> GeomAbsShape {
        c.curve.continuity()
    }

    /// OCCT NbIntervals(C, S) (hxx L61-64).
    fn nb_intervals(c: &CurveHost<'a>, s: GeomAbsShape) -> usize {
        c.curve.nb_intervals(s)
    }

    /// OCCT Intervals(C, T, S) (hxx L67-72).
    fn intervals(c: &CurveHost<'a>, t: &mut [f64], s: GeomAbsShape) {
        let a_tab = c.curve.intervals(s);
        t[..a_tab.len()].copy_from_slice(&a_tab);
    }

    /// OCCT IsClosed(C) (hxx L74-77).
    fn is_closed(c: &CurveHost<'a>) -> bool {
        c.curve.is_closed()
    }

    /// OCCT IsPeriodic(C) (hxx L79-82).
    fn is_periodic(c: &CurveHost<'a>) -> bool {
        c.curve.is_periodic()
    }

    /// OCCT Period(C) (hxx L84-87).
    fn period(c: &CurveHost<'a>) -> f64 {
        c.curve.period()
    }

    /// OCCT Value(C, U) (hxx L89-92).
    fn value(c: &CurveHost<'a>, u: f64) -> Point3 {
        c.curve.value(u)
    }

    /// OCCT D0(C, U, P) (hxx L94-97).
    fn d0(c: &CurveHost<'a>, u: f64) -> Point3 {
        c.curve.value(u)
    }

    /// OCCT D1(C, U, P, V) (hxx L99-103).
    fn d1(c: &CurveHost<'a>, u: f64) -> (Point3, Vec3) {
        c.curve.d1(u)
    }

    /// OCCT D2(C, U, P, V1, V2) (hxx L108-112).
    fn d2(c: &CurveHost<'a>, u: f64) -> (Point3, Vec3, Vec3) {
        c.curve.d2(u)
    }

    /// OCCT D3(C, U, P, V1, V2, V3) (hxx L118-123).
    fn d3(c: &CurveHost<'a>, u: f64) -> (Point3, Vec3, Vec3, Vec3) {
        c.curve.d3(u)
    }

    /// OCCT DN(C, U, N) (hxx L127-130).
    fn dn(c: &CurveHost<'a>, u: f64, n: usize) -> Vec3 {
        c.curve.dn(u, n as i32)
    }

    /// OCCT Resolution(C, R3d) (hxx L134-137).
    fn resolution(c: &CurveHost<'a>, r3d: f64) -> f64 {
        c.curve.resolution(r3d)
    }

    /// OCCT GetType(C) (hxx L142-145).
    fn get_type(c: &CurveHost<'a>) -> CurveType {
        match c.geom {
            Some(g) => g.get_type(),
            None => c.curve.curve_type(),
        }
    }

    /// OCCT Line(C) (hxx L147-150) — Standard_NoSuchObject on the base class.
    fn line(c: &CurveHost<'a>) -> Line3 {
        c.geom
            .expect("Standard_NoSuchObject: IntCurveSurface_TheHCurveTool::Line")
            .line()
    }

    /// OCCT Circle(C) (hxx L152-155).
    fn circle(c: &CurveHost<'a>) -> Circle3 {
        c.geom
            .expect("Standard_NoSuchObject: IntCurveSurface_TheHCurveTool::Circle")
            .circle()
    }

    /// OCCT Ellipse(C) (hxx L157-160).
    fn ellipse(c: &CurveHost<'a>) -> Ellipse3 {
        c.geom
            .expect("Standard_NoSuchObject: IntCurveSurface_TheHCurveTool::Ellipse")
            .ellipse()
    }

    /// OCCT Hyperbola(C) (hxx L162-165).
    fn hyperbola(c: &CurveHost<'a>) -> Hyperbola3 {
        c.geom
            .expect("Standard_NoSuchObject: IntCurveSurface_TheHCurveTool::Hyperbola")
            .hyperbola()
    }

    /// OCCT Parabola(C) (hxx L167-170).
    fn parabola(c: &CurveHost<'a>) -> Parabola3 {
        c.geom
            .expect("Standard_NoSuchObject: IntCurveSurface_TheHCurveTool::Parabola")
            .parabola()
    }

    /// OCCT Bezier(C) (hxx L172-175) — Standard_NoSuchObject on the base
    /// class; answered from the kernel payload when the caller holds it.
    fn bezier(c: &CurveHost<'a>) -> &'a BezierCurve3 {
        match c.curve3 {
            Some(Curve3::Bezier(b)) => b,
            _ => panic!("Standard_NoSuchObject: IntCurveSurface_TheHCurveTool::Bezier"),
        }
    }

    /// OCCT BSpline(C) (hxx L177-180).
    fn bspline(c: &CurveHost<'a>) -> &'a BSplineCurve3 {
        match c.curve3 {
            Some(Curve3::BSpline(b)) => b,
            _ => panic!("Standard_NoSuchObject: IntCurveSurface_TheHCurveTool::BSpline"),
        }
    }
}

impl<'a> CurveTool3d for TheCurveTool<'a> {
    type Curve = CurveHost<'a>;

    fn value(c: &CurveHost<'a>, u: f64) -> Point3 {
        <Self as HCurveTool>::value(c, u)
    }

    fn d1(c: &CurveHost<'a>, u: f64) -> (Point3, Vec3) {
        <Self as HCurveTool>::d1(c, u)
    }

    fn first_parameter(c: &CurveHost<'a>) -> f64 {
        <Self as HCurveTool>::first_parameter(c)
    }

    fn last_parameter(c: &CurveHost<'a>) -> f64 {
        <Self as HCurveTool>::last_parameter(c)
    }

    fn resolution(c: &CurveHost<'a>, r3d: f64) -> f64 {
        <Self as HCurveTool>::resolution(c, r3d)
    }
}

/// OCCT Adaptor3d_HSurfaceTool (hxx L40-293 + .cxx L27-115) over
/// [`GeomSurfaceAdaptor`] — the static tool over the `GeomAdaptor_Surface`
/// the caller holds; the OCCT template parameter `TheSurfaceTool` of the
/// HInter instantiation.  Every accessor forwards to the adaptor instance
/// exactly as the OCCT static class forwards to `Adaptor3d_Surface`.
pub struct TheSurfaceTool;

impl HSurfaceTool for TheSurfaceTool {
    type Surface = GeomSurfaceAdaptor;
    /// OCCT BasisCurve(S) — the GeomAdaptor_Curve view of the basis curve.
    type BasisCurve = GeomCurveAdaptor;
    /// OCCT BasisSurface(S) — the GeomAdaptor_Surface view of the basis
    /// surface.
    type BasisSurface = GeomSurfaceAdaptor;

    /// OCCT FirstUParameter(S) (hxx L45-48).
    fn first_u_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        s.first_u_parameter()
    }

    /// OCCT FirstVParameter(S) (hxx L50-53).
    fn first_v_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        s.first_v_parameter()
    }

    /// OCCT LastUParameter(S) (hxx L55-58).
    fn last_u_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        s.last_u_parameter()
    }

    /// OCCT LastVParameter(S) (hxx L60-63).
    fn last_v_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        s.last_v_parameter()
    }

    /// OCCT NbUIntervals(S, Sh) (hxx L65-68).
    fn nb_u_intervals(s: &GeomSurfaceAdaptor, sh: GeomAbsShape) -> usize {
        s.nb_u_intervals(sh)
    }

    /// OCCT NbVIntervals(S, Sh) (hxx L70-73).
    fn nb_v_intervals(s: &GeomSurfaceAdaptor, sh: GeomAbsShape) -> usize {
        s.nb_v_intervals(sh)
    }

    /// OCCT UIntervals(S, Tab, Sh) (hxx L75-80).
    fn u_intervals(s: &GeomSurfaceAdaptor, tab: &mut [f64], sh: GeomAbsShape) {
        let a_tab = s.u_intervals(sh);
        tab[..a_tab.len()].copy_from_slice(&a_tab);
    }

    /// OCCT VIntervals(S, Tab, Sh) (hxx L82-87).
    fn v_intervals(s: &GeomSurfaceAdaptor, tab: &mut [f64], sh: GeomAbsShape) {
        let a_tab = s.v_intervals(sh);
        tab[..a_tab.len()].copy_from_slice(&a_tab);
    }

    /// OCCT UTrim(S, First, Last, Tol) (hxx L90-96) — the narrowed adaptor.
    fn u_trim(s: &GeomSurfaceAdaptor, first: f64, last: f64, tol: f64) -> GeomSurfaceAdaptor {
        s.u_trim_of(first, last, tol)
    }

    /// OCCT VTrim(S, First, Last, Tol) (hxx L99-105).
    fn v_trim(s: &GeomSurfaceAdaptor, first: f64, last: f64, tol: f64) -> GeomSurfaceAdaptor {
        s.v_trim_of(first, last, tol)
    }

    /// OCCT IsUClosed(S) (hxx L107-110).
    fn is_u_closed(s: &GeomSurfaceAdaptor) -> bool {
        s.is_u_closed()
    }

    /// OCCT IsVClosed(S) (hxx L112-115).
    fn is_v_closed(s: &GeomSurfaceAdaptor) -> bool {
        s.is_v_closed()
    }

    /// OCCT IsUPeriodic(S) (hxx L117-120).
    fn is_u_periodic(s: &GeomSurfaceAdaptor) -> bool {
        s.is_u_periodic()
    }

    /// OCCT UPeriod(S) (hxx L122-125).
    fn u_period(s: &GeomSurfaceAdaptor) -> f64 {
        s.u_period()
    }

    /// OCCT IsVPeriodic(S) (hxx L127-130).
    fn is_v_periodic(s: &GeomSurfaceAdaptor) -> bool {
        s.is_v_periodic()
    }

    /// OCCT VPeriod(S) (hxx L132-135).
    fn v_period(s: &GeomSurfaceAdaptor) -> f64 {
        s.v_period()
    }

    /// OCCT Value(S, U, V) (hxx L137-142).
    fn value(s: &GeomSurfaceAdaptor, u: f64, v: f64) -> Point3 {
        s.value(u, v)
    }

    /// OCCT D0(S, U, V, P) (hxx L144-150).
    fn d0(s: &GeomSurfaceAdaptor, u: f64, v: f64) -> Point3 {
        s.value(u, v)
    }

    /// OCCT D1(S, U, V, P, D1U, D1V) (hxx L152-160).
    fn d1(s: &GeomSurfaceAdaptor, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        s.d1(u, v)
    }

    /// OCCT D2(S, U, V, P, D1U, D1V, D2U, D2V, D2UV) (hxx L162-173).
    fn d2(
        s: &GeomSurfaceAdaptor,
        u: f64,
        v: f64,
    ) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        s.d2(u, v)
    }

    /// OCCT DN(S, U, V, NU, NV) (hxx L203-210).
    fn dn(s: &GeomSurfaceAdaptor, u: f64, v: f64, nu: usize, nv: usize) -> Vec3 {
        s.dn(u, v, nu as i32, nv as i32)
    }

    /// OCCT UResolution(S, R3d) (hxx L212-215).
    fn u_resolution(s: &GeomSurfaceAdaptor, r3d: f64) -> f64 {
        s.u_resolution(r3d)
    }

    /// OCCT VResolution(S, R3d) (hxx L217-220).
    fn v_resolution(s: &GeomSurfaceAdaptor, r3d: f64) -> f64 {
        s.v_resolution(r3d)
    }

    /// OCCT GetType(S) (hxx L222-225) — the GeomAdaptor_Surface type tag.
    fn get_type(s: &GeomSurfaceAdaptor) -> SurfaceType {
        match Adaptor3dSurface::get_type(s) {
            KernelSurfaceType::Plane => SurfaceType::Plane,
            KernelSurfaceType::Cylinder => SurfaceType::Cylinder,
            KernelSurfaceType::Cone => SurfaceType::Cone,
            KernelSurfaceType::Sphere => SurfaceType::Sphere,
            KernelSurfaceType::Torus => SurfaceType::Torus,
            KernelSurfaceType::BezierSurface => SurfaceType::BezierSurface,
            KernelSurfaceType::BSplineSurface => SurfaceType::BSplineSurface,
            KernelSurfaceType::SurfaceOfRevolution => SurfaceType::SurfaceOfRevolution,
            KernelSurfaceType::SurfaceOfExtrusion => SurfaceType::SurfaceOfExtrusion,
            KernelSurfaceType::OffsetSurface => SurfaceType::OffsetSurface,
            KernelSurfaceType::OtherSurface => SurfaceType::OtherSurface,
        }
    }

    /// OCCT Plane(S) (hxx L227-230).
    fn plane(s: &GeomSurfaceAdaptor) -> Plane {
        s.plane()
    }

    /// OCCT Cylinder(S) (hxx L229-232).
    fn cylinder(s: &GeomSurfaceAdaptor) -> CylindricalSurface {
        s.cylinder()
    }

    /// OCCT Cone(S) (hxx L234-237).
    fn cone(s: &GeomSurfaceAdaptor) -> ConicalSurface {
        s.cone()
    }

    /// OCCT Torus(S) (hxx L236-239).
    fn torus(s: &GeomSurfaceAdaptor) -> ToroidalSurface {
        s.torus()
    }

    /// OCCT Sphere(S) (hxx L238-241).
    fn sphere(s: &GeomSurfaceAdaptor) -> SphericalSurface {
        s.sphere()
    }

    /// OCCT AxeOfRevolution(S) (hxx L253-256).
    fn axe_of_revolution(s: &GeomSurfaceAdaptor) -> (Point3, Vec3) {
        s.axe_of_revolution()
    }

    /// OCCT Direction(S) (hxx L258-261).
    fn direction(s: &GeomSurfaceAdaptor) -> Vec3 {
        s.direction()
    }

    /// OCCT BasisCurve(S) (hxx L263-266).
    fn basis_curve(s: &GeomSurfaceAdaptor) -> GeomCurveAdaptor {
        s.basis_curve_of()
    }

    /// OCCT BasisSurface(S) (hxx L268-271).
    fn basis_surface(s: &GeomSurfaceAdaptor) -> GeomSurfaceAdaptor {
        s.basis_surface_of()
    }

    /// OCCT OffsetValue(S) (hxx L273-276).
    fn offset_value(s: &GeomSurfaceAdaptor) -> f64 {
        s.offset_value()
    }

    /// OCCT Adaptor3d_Surface::NbUPoles() — the NbSamplesU input.
    fn nb_u_poles(s: &GeomSurfaceAdaptor) -> usize {
        s.nb_u_poles_of()
    }

    /// OCCT Adaptor3d_Surface::NbVPoles().
    fn nb_v_poles(s: &GeomSurfaceAdaptor) -> usize {
        s.nb_v_poles_of()
    }

    /// OCCT Adaptor3d_Surface::NbUKnots().
    fn nb_u_knots(s: &GeomSurfaceAdaptor) -> usize {
        s.nb_u_knots_of()
    }

    /// OCCT Adaptor3d_Surface::NbVKnots().
    fn nb_v_knots(s: &GeomSurfaceAdaptor) -> usize {
        s.nb_v_knots_of()
    }

    /// OCCT Adaptor3d_Surface::UDegree().
    fn u_degree(s: &GeomSurfaceAdaptor) -> usize {
        s.u_degree_of()
    }

    /// OCCT Adaptor3d_Surface::VDegree().
    fn v_degree(s: &GeomSurfaceAdaptor) -> usize {
        s.v_degree_of()
    }

    /// OCCT Adaptor3d_Surface::Bezier().
    fn bezier(s: &GeomSurfaceAdaptor) -> &BezierSurface {
        match s.surface() {
            rcad_kernel::geom::Surface3::Bezier(b) => b,
            _ => panic!("Standard_NoSuchObject: Adaptor3d_HSurfaceTool::Bezier"),
        }
    }

    /// OCCT Adaptor3d_Surface::BSpline().
    fn bspline(s: &GeomSurfaceAdaptor) -> &BSplineSurface {
        match s.surface() {
            rcad_kernel::geom::Surface3::BSpline(b) => b,
            _ => panic!("Standard_NoSuchObject: Adaptor3d_HSurfaceTool::BSpline"),
        }
    }
}

impl PSurfaceTool for TheSurfaceTool {
    type Surface = GeomSurfaceAdaptor;

    fn value(s: &GeomSurfaceAdaptor, u: f64, v: f64) -> Point3 {
        <TheSurfaceTool as HSurfaceTool>::value(s, u, v)
    }

    fn d1(s: &GeomSurfaceAdaptor, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        <TheSurfaceTool as HSurfaceTool>::d1(s, u, v)
    }

    fn first_u_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        <TheSurfaceTool as HSurfaceTool>::first_u_parameter(s)
    }

    fn last_u_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        <TheSurfaceTool as HSurfaceTool>::last_u_parameter(s)
    }

    fn first_v_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        <TheSurfaceTool as HSurfaceTool>::first_v_parameter(s)
    }

    fn last_v_parameter(s: &GeomSurfaceAdaptor) -> f64 {
        <TheSurfaceTool as HSurfaceTool>::last_v_parameter(s)
    }

    fn u_resolution(s: &GeomSurfaceAdaptor, r3d: f64) -> f64 {
        <TheSurfaceTool as HSurfaceTool>::u_resolution(s, r3d)
    }

    fn v_resolution(s: &GeomSurfaceAdaptor, r3d: f64) -> f64 {
        <TheSurfaceTool as HSurfaceTool>::v_resolution(s, r3d)
    }
}

/// OCCT Adaptor3d_Curve::GetType() over the GeomAdaptor_Curve basis-curve
/// view — the Adaptor3dSurfaceBasis input of the EstLim* offset paths.
impl Adaptor3dCurveBasis for GeomCurveAdaptor {
    fn get_type(&self) -> CurveType {
        Adaptor3dCurveGeom::get_type(self)
    }

    fn value(&self, u: f64) -> Point3 {
        Adaptor3dCurve::value(self, u)
    }

    fn line(&self) -> Line3 {
        GeomCurveAdaptor::line(self)
    }

    fn parabola(&self) -> Parabola3 {
        GeomCurveAdaptor::parabola(self)
    }

    fn hyperbola(&self) -> Hyperbola3 {
        GeomCurveAdaptor::hyperbola(self)
    }
}

/// OCCT Adaptor3d_Surface accessors over the GeomAdaptor_Surface
/// basis-surface view — the Adaptor3dSurfaceBasis input of the EstLim*
/// offset paths.
impl Adaptor3dSurfaceBasis<GeomCurveAdaptor> for GeomSurfaceAdaptor {
    fn get_type(&self) -> SurfaceType {
        match Adaptor3dSurface::get_type(self) {
            KernelSurfaceType::Plane => SurfaceType::Plane,
            KernelSurfaceType::Cylinder => SurfaceType::Cylinder,
            KernelSurfaceType::Cone => SurfaceType::Cone,
            KernelSurfaceType::Sphere => SurfaceType::Sphere,
            KernelSurfaceType::Torus => SurfaceType::Torus,
            KernelSurfaceType::BezierSurface => SurfaceType::BezierSurface,
            KernelSurfaceType::BSplineSurface => SurfaceType::BSplineSurface,
            KernelSurfaceType::SurfaceOfRevolution => SurfaceType::SurfaceOfRevolution,
            KernelSurfaceType::SurfaceOfExtrusion => SurfaceType::SurfaceOfExtrusion,
            KernelSurfaceType::OffsetSurface => SurfaceType::OffsetSurface,
            KernelSurfaceType::OtherSurface => SurfaceType::OtherSurface,
        }
    }

    fn plane(&self) -> Plane {
        GeomSurfaceAdaptor::plane(self)
    }

    fn cylinder(&self) -> CylindricalSurface {
        GeomSurfaceAdaptor::cylinder(self)
    }

    fn cone(&self) -> ConicalSurface {
        GeomSurfaceAdaptor::cone(self)
    }

    fn direction(&self) -> Vec3 {
        GeomSurfaceAdaptor::direction(self)
    }

    fn basis_curve(&self) -> GeomCurveAdaptor {
        GeomSurfaceAdaptor::basis_curve_of(self)
    }
}

/// OCCT IntCurveSurface_HInter (HInter.hxx L47-208) over the production
/// adaptor pair — extends IntCurveSurface_Intersection (the engine
/// [`HInter`], the `base` subobject).
#[derive(Debug, Clone, Default)]
pub struct IntCurveSurfaceHInter {
    /// OCCT IntCurveSurface_Intersection base subobject.
    pub base: HInter,
}

impl IntCurveSurfaceHInter {
    /// OCCT IntCurveSurface_HInter() (cxx L56).
    pub fn new() -> Self {
        IntCurveSurfaceHInter {
            base: HInter::new(),
        }
    }

    /// OCCT Perform(Curve, Surface) (cxx L106-116) over the adaptor pair —
    /// instantiates the Inter.pxx templates with TheCurveTool =
    /// IntCurveSurface_TheHCurveTool over the erased Adaptor3d_Curve handle
    /// and TheSurfaceTool = Adaptor3d_HSurfaceTool over the
    /// GeomAdaptor_Surface (the GeomFill_SectionPlacement.cxx L471-472
    /// statement).
    pub fn perform_adaptors(&mut self, curve: &CurveHost, surface: &GeomSurfaceAdaptor) {
        fn drive<'a>(
            this: &mut IntCurveSurfaceHInter,
            curve: &CurveHost<'a>,
            surface: &GeomSurfaceAdaptor,
        ) {
            this.base.perform::<CurveHost<'a>, TheCurveTool<'a>, GeomSurfaceAdaptor, TheSurfaceTool>(
                curve, surface,
            );
        }
        drive(self, curve, surface);
    }

    /// OCCT Perform(Curve, Surface) over the kernel geometry encodings —
    /// the OCCT callers lift their Geom_Curve / Geom_Surface through
    /// GeomAdaptor_Curve / GeomAdaptor_Surface before the Perform
    /// (GeomFill_SectionPlacement.cxx L469-470
    /// GeomFill_GuideTrihedronPlan.cxx Init); this entry performs the same
    /// lift.
    pub fn perform(&mut self, curve: &Curve3, surface: &Surface3) {
        let a_curve_adaptor = GeomCurveAdaptor::new(curve.clone());
        let a_surface_adaptor = GeomSurfaceAdaptor::new(surface.clone());
        let a_host = CurveHost {
            curve: &a_curve_adaptor,
            geom: Some(&a_curve_adaptor),
            curve3: Some(curve),
        };
        self.perform_adaptors(&a_host, &a_surface_adaptor);
    }

    /// OCCT IsDone() (Intersection.cxx L33-36).
    pub fn is_done(&self) -> bool {
        self.base.base.is_done()
    }

    /// OCCT NbPoints() (Intersection.cxx L45-52).
    pub fn nb_points(&self) -> usize {
        self.base.base.nb_points()
    }

    /// OCCT Point(N) (Intersection.cxx L65-72) — 1-based.
    pub fn point(&self, j: usize) -> &IntersectionPoint {
        self.base.base.point(j)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::geom::Line3;

    /// OCCT anchor (the GeomFill_SectionPlacement.cxx L468-472 statement
    /// shape): a line adaptor piercing the plane z = 0.  C(w) = (0, 0.25, 1)
    /// + w * (1, 0, -1)/sqrt(2) over [-5, 5] hits z = 0 at w = sqrt(2),
    /// point (1, 0.25, 0), u = 1, v = 0.25.
    #[test]
    fn hinter_adaptor_pair_line_plane_anchor() {
        let line = Line3::new(DVec3::new(0.0, 0.25, 1.0), DVec3::new(1.0, 0.0, -1.0).normalize());
        let curve = Curve3::Line(line);
        let a_curve_adaptor = GeomCurveAdaptor::with_range(curve.clone(), -5.0, 5.0);
        let a_host = CurveHost {
            curve: &a_curve_adaptor,
            geom: Some(&a_curve_adaptor),
            curve3: Some(&curve),
        };
        let a_surface_adaptor = GeomSurfaceAdaptor::new(Surface3::Plane(Plane::new(
            DVec3::ZERO,
            DVec3::Z,
        )));

        let mut intersector = IntCurveSurfaceHInter::new();
        intersector.perform_adaptors(&a_host, &a_surface_adaptor);

        assert!(intersector.is_done());
        assert_eq!(intersector.nb_points(), 1);
        let pt = intersector.point(1);
        assert!((pt.w() - std::f64::consts::SQRT_2).abs() < 1e-9, "w={}", pt.w());
        assert!((pt.u() - 1.0).abs() < 1e-9, "u={}", pt.u());
        assert!((pt.v() - 0.25).abs() < 1e-9, "v={}", pt.v());
        assert!(pt.pnt().distance(DVec3::new(1.0, 0.25, 0.0)) < 1e-9);
    }

    /// The same anchor perturbed: the plane dropped to z = -0.5 moves the
    /// crossing to w = 1.5*sqrt(2) — the intersector tracks the geometry.
    #[test]
    fn hinter_adaptor_pair_line_plane_perturbed() {
        let line = Line3::new(DVec3::new(0.0, 0.25, 1.0), DVec3::new(1.0, 0.0, -1.0).normalize());
        let curve = Curve3::Line(line);
        let a_curve_adaptor = GeomCurveAdaptor::with_range(curve.clone(), -5.0, 5.0);
        let a_host = CurveHost {
            curve: &a_curve_adaptor,
            geom: Some(&a_curve_adaptor),
            curve3: Some(&curve),
        };
        let a_surface_adaptor = GeomSurfaceAdaptor::new(Surface3::Plane(Plane::new(
            DVec3::new(0.0, 0.0, -0.5),
            DVec3::Z,
        )));

        let mut intersector = IntCurveSurfaceHInter::new();
        intersector.perform_adaptors(&a_host, &a_surface_adaptor);

        assert!(intersector.is_done());
        assert_eq!(intersector.nb_points(), 1);
        let pt = intersector.point(1);
        assert!(
            (pt.w() - 1.5 * std::f64::consts::SQRT_2).abs() < 1e-9,
            "w={}",
            pt.w()
        );
    }

    /// The kernel-geometry entry (the GeomFill callers' lift): a BSpline
    /// guide against the plane z = 0 — the polygon route over the erased
    /// encoding with the payload view present.  The quadratic Bezier
    /// (-1,0,1)-(0,0,-2)-(1,0,1) crosses z = 0 at w = (3 +- sqrt(3))/6.
    #[test]
    fn hinter_kernel_entry_bezier_plane() {
        let curve = Curve3::Bezier(rcad_kernel::geom::BezierCurve3 {
            control_points: vec![
                DVec3::new(-1.0, 0.0, 1.0),
                DVec3::new(0.0, 0.0, -2.0),
                DVec3::new(1.0, 0.0, 1.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        let surface = Surface3::Plane(Plane::new(DVec3::ZERO, DVec3::Z));

        let mut intersector = IntCurveSurfaceHInter::new();
        intersector.perform(&curve, &surface);

        assert!(intersector.is_done());
        assert_eq!(intersector.nb_points(), 2, "points: {:?}", intersector.base.base.lpnt);
        let s3 = 3.0f64.sqrt();
        let w1 = (3.0 - s3) / 6.0;
        let w2 = (3.0 + s3) / 6.0;
        let p1 = intersector.point(1);
        let p2 = intersector.point(2);
        assert!((p1.w() - w1).abs() < 1e-6, "w1={} vs {}", p1.w(), w1);
        assert!((p2.w() - w2).abs() < 1e-6, "w2={} vs {}", p2.w(), w2);
        assert!(p1.pnt().z.abs() < 1e-9 && p2.pnt().z.abs() < 1e-9);
    }
}
