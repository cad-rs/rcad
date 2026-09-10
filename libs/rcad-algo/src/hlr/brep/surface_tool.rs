// OCCT HLRBRep_LineTool + HLRBRep_SurfaceTool (TKHLR) — the static tools
// over the HLR line (gp_Lin) and the HLR face surface (HLRBRep_Surface),
// bound into the IntCurveSurface / IntImp tool traits.
//
// HLRBRep_LineTool.hxx L47-177 + .lxx (infinite bounds, CN continuity,
// ElCLib evaluations, empty Intervals/Poles/Knots bodies, NbSamples = 3).
// HLRBRep_SurfaceTool.hxx L41-175 + .lxx/.cxx L14-157 (every member a
// one-line forward to the HLRBRep_Surface adaptor; NbSamplesU/V dispatch).

use glam::DVec3;
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3, Point3, Vec3};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::int_curve_surface::{Adaptor3dCurveBasis, Adaptor3dSurfaceBasis, HCurveTool, HSurfaceTool};
use crate::geomalgo::int_imp::{CurveTool3d, PSurfaceTool};
use crate::geomalgo::int_patch::GeomAbsSurfaceType;

use super::surface::Surface;

// ---------------------------------------------------------------------------
// HLRBRep_LineTool
// ---------------------------------------------------------------------------

/// OCCT HLRBRep_LineTool — the static tool over gp_Lin.
#[derive(Debug, Clone, Copy, Default)]
pub struct LineTool;

impl HCurveTool for LineTool {
    type Curve = Line3;

    /// OCCT FirstParameter(C) — RealFirst().
    fn first_parameter(_c: &Line3) -> f64 {
        f64::MIN
    }
    /// OCCT LastParameter(C) — RealLast().
    fn last_parameter(_c: &Line3) -> f64 {
        f64::MAX
    }
    /// OCCT Continuity(C) — GeomAbs_CN.
    fn continuity(_c: &Line3) -> GeomAbsShape {
        GeomAbsShape::CN
    }
    /// OCCT NbIntervals(C, S) — 1.
    fn nb_intervals(_c: &Line3, _s: GeomAbsShape) -> usize {
        1
    }
    /// OCCT Intervals(C, T, S) — empty body verbatim.
    fn intervals(_c: &Line3, _t: &mut [f64], _s: GeomAbsShape) {}
    /// OCCT IsClosed(C).
    fn is_closed(_c: &Line3) -> bool {
        false
    }
    /// OCCT IsPeriodic(C).
    fn is_periodic(_c: &Line3) -> bool {
        false
    }
    /// OCCT Period(C).
    fn period(_c: &Line3) -> f64 {
        0.0
    }
    /// OCCT Value(C, U) / D0 — ElCLib::Value.
    fn value(c: &Line3, u: f64) -> Point3 {
        rcad_kernel::math::el::elclib_line_value(u, c.origin, c.direction)
    }
    fn d0(c: &Line3, u: f64) -> Point3 {
        <Self as HCurveTool>::value(c, u)
    }
    /// OCCT D1(C, U, P, V) — ElCLib::D1.
    fn d1(c: &Line3, u: f64) -> (Point3, Vec3) {
        rcad_kernel::math::el::elclib_line_d1(u, c.origin, c.direction)
    }
    /// OCCT D2(C, U, P, V1, V2) — ElCLib::D1 + V2 = 0.
    fn d2(c: &Line3, u: f64) -> (Point3, Vec3, Vec3) {
        let (p, v1) = <Self as HCurveTool>::d1(c, u);
        (p, v1, Vec3::ZERO)
    }
    /// OCCT D3 — ElCLib::D1 + V2 = V3 = 0.
    fn d3(c: &Line3, u: f64) -> (Point3, Vec3, Vec3, Vec3) {
        let (p, v1) = <Self as HCurveTool>::d1(c, u);
        (p, v1, Vec3::ZERO, Vec3::ZERO)
    }
    /// OCCT DN(C, U, N) — the direction for N == 1, zero beyond.
    fn dn(c: &Line3, _u: f64, n: usize) -> Vec3 {
        if n == 1 {
            c.direction
        } else {
            Vec3::ZERO
        }
    }
    /// OCCT Resolution(C, R3D) — R3D.
    fn resolution(_c: &Line3, r3d: f64) -> f64 {
        r3d
    }
    /// OCCT GetType(C) — GeomAbs_Line.
    fn get_type(_c: &Line3) -> CurveType {
        CurveType::Line
    }
    /// OCCT Line(C) — identity.
    fn line(c: &Line3) -> Line3 {
        c.clone()
    }
    /// OCCT Circle(C) — default-constructed gp_Circ (never invoked for a
    /// line); the rcad panic is the unreachable-path marker.
    fn circle(_c: &Line3) -> Circle3 {
        panic!("Standard_NoSuchObject: LineTool::Circle");
    }
    fn ellipse(_c: &Line3) -> Ellipse3 {
        panic!("Standard_NoSuchObject: LineTool::Ellipse");
    }
    fn hyperbola(_c: &Line3) -> Hyperbola3 {
        panic!("Standard_NoSuchObject: LineTool::Hyperbola");
    }
    fn parabola(_c: &Line3) -> Parabola3 {
        panic!("Standard_NoSuchObject: LineTool::Parabola");
    }
    /// OCCT Bezier/BSpline(C) — null handles (unreachable for a line).
    fn bezier(_c: &Line3) -> &rcad_kernel::geom::BezierCurve3 {
        panic!("Standard_NoSuchObject: LineTool::Bezier");
    }
    fn bspline(_c: &Line3) -> &rcad_kernel::geom::BSplineCurve3 {
        panic!("Standard_NoSuchObject: LineTool::BSpline");
    }
    /// OCCT NbSamples(C, U0, U1) — 3 (lxx).
    fn nb_samples(_c: &Line3, _u0: f64, _u1: f64) -> usize {
        3
    }
}

impl CurveTool3d for LineTool {
    type Curve = Line3;
    /// OCCT TheCurveTool::Value(C, U).
    fn value(c: &Line3, u: f64) -> DVec3 {
        <Self as HCurveTool>::value(c, u)
    }
    /// OCCT TheCurveTool::D1(C, U, P, T).
    fn d1(c: &Line3, u: f64) -> (DVec3, DVec3) {
        <Self as HCurveTool>::d1(c, u)
    }
    fn first_parameter(c: &Line3) -> f64 {
        <Self as HCurveTool>::first_parameter(c)
    }
    fn last_parameter(c: &Line3) -> f64 {
        <Self as HCurveTool>::last_parameter(c)
    }
    fn resolution(c: &Line3, r3d: f64) -> f64 {
        <Self as HCurveTool>::resolution(c, r3d)
    }
}

// ---------------------------------------------------------------------------
// HLRBRep_SurfaceTool
// ---------------------------------------------------------------------------

/// The null basis-curve adaptor (OCCT BRepAdaptor_Surface::BasisCurve
/// raises Standard_NoSuchObject for the surfaces the HLR chain feeds in).
#[derive(Debug, Clone, Copy, Default)]
pub struct NullBasisCurve;

impl Adaptor3dCurveBasis for NullBasisCurve {
    fn get_type(&self) -> rcad_kernel::base::proj_lib::CurveType {
        panic!("Standard_NoSuchObject: BasisCurve");
    }
    fn value(&self, _u: f64) -> Point3 {
        panic!("Standard_NoSuchObject: BasisCurve");
    }
    fn line(&self) -> Line3 {
        panic!("Standard_NoSuchObject: BasisCurve");
    }
    fn parabola(&self) -> rcad_kernel::geom::Parabola3 {
        panic!("Standard_NoSuchObject: BasisCurve");
    }
    fn hyperbola(&self) -> Hyperbola3 {
        panic!("Standard_NoSuchObject: BasisCurve");
    }
}

/// The null basis-surface adaptor (BRepAdaptor_Surface::BasisSurface).
#[derive(Debug, Clone, Copy, Default)]
pub struct NullBasisSurface;

impl Adaptor3dSurfaceBasis<NullBasisCurve> for NullBasisSurface {
    fn get_type(&self) -> GeomAbsSurfaceType {
        panic!("Standard_NoSuchObject: BasisSurface");
    }
    fn plane(&self) -> rcad_kernel::geom::Plane {
        panic!("Standard_NoSuchObject: BasisSurface");
    }
    fn cylinder(&self) -> rcad_kernel::geom::CylindricalSurface {
        panic!("Standard_NoSuchObject: BasisSurface");
    }
    fn cone(&self) -> rcad_kernel::geom::ConicalSurface {
        panic!("Standard_NoSuchObject: BasisSurface");
    }
    fn direction(&self) -> Vec3 {
        panic!("Standard_NoSuchObject: BasisSurface");
    }
    fn basis_curve(&self) -> NullBasisCurve {
        NullBasisCurve
    }
}

/// OCCT HLRBRep_SurfaceTool — the static tool over HLRBRep_Surface; every
/// member forwards to the surface adaptor.
#[derive(Debug, Clone, Copy, Default)]
pub struct SurfaceTool<'a> {
    _life: std::marker::PhantomData<&'a ()>,
}

impl<'a> SurfaceTool<'a> {
    /// OCCT UTrim forwards to BRepAdaptor_Surface::UTrim — the narrowed
    /// window on the same surface (GeomAdaptor semantics).
    pub fn u_trim(s: &Surface<'a>, first: f64, last: f64, tol: f64) -> Surface<'a> {
        s.u_trim(first, last, tol)
    }
    /// OCCT VTrim.
    pub fn v_trim(s: &Surface<'a>, first: f64, last: f64, tol: f64) -> Surface<'a> {
        s.v_trim(first, last, tol)
    }
}

impl<'a> HSurfaceTool for SurfaceTool<'a> {
    type Surface = Surface<'a>;
    type BasisCurve = NullBasisCurve;
    type BasisSurface = NullBasisSurface;

    fn first_u_parameter(s: &Surface<'a>) -> f64 {
        s.first_u_parameter()
    }
    fn first_v_parameter(s: &Surface<'a>) -> f64 {
        s.first_v_parameter()
    }
    fn last_u_parameter(s: &Surface<'a>) -> f64 {
        s.last_u_parameter()
    }
    fn last_v_parameter(s: &Surface<'a>) -> f64 {
        s.last_v_parameter()
    }
    fn nb_u_intervals(_s: &Surface<'a>, _sh: GeomAbsShape) -> usize {
        1
    }
    fn nb_v_intervals(_s: &Surface<'a>, _sh: GeomAbsShape) -> usize {
        1
    }
    fn u_intervals(s: &Surface<'a>, tab: &mut [f64], _sh: GeomAbsShape) {
        tab[0] = s.first_u_parameter();
        tab[1] = s.last_u_parameter();
    }
    fn v_intervals(s: &Surface<'a>, tab: &mut [f64], _sh: GeomAbsShape) {
        tab[0] = s.first_v_parameter();
        tab[1] = s.last_v_parameter();
    }
    fn u_trim(s: &Surface<'a>, first: f64, last: f64, tol: f64) -> Surface<'a> {
        Self::u_trim(s, first, last, tol)
    }
    fn v_trim(s: &Surface<'a>, first: f64, last: f64, tol: f64) -> Surface<'a> {
        Self::v_trim(s, first, last, tol)
    }
    fn is_u_closed(s: &Surface<'a>) -> bool {
        Surface::is_u_closed(s)
    }
    fn is_v_closed(s: &Surface<'a>) -> bool {
        Surface::is_v_closed(s)
    }
    fn is_u_periodic(s: &Surface<'a>) -> bool {
        s.is_u_periodic()
    }
    fn u_period(s: &Surface<'a>) -> f64 {
        s.u_period()
    }
    fn is_v_periodic(s: &Surface<'a>) -> bool {
        s.is_v_periodic()
    }
    fn v_period(s: &Surface<'a>) -> f64 {
        s.v_period()
    }
    fn value(s: &Surface<'a>, u: f64, v: f64) -> Point3 {
        s.value(u, v)
    }
    fn d0(s: &Surface<'a>, u: f64, v: f64) -> Point3 {
        s.d0(u, v)
    }
    fn d1(s: &Surface<'a>, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        s.d1(u, v)
    }
    fn d2(s: &Surface<'a>, u: f64, v: f64) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        s.d2(u, v)
    }
    fn dn(_s: &Surface<'a>, _u: f64, _v: f64, _nu: usize, _nv: usize) -> Vec3 {
        // BRepAdaptor_Surface::DN on the analytic supports is not exercised
        // by the HLR chain (the walking uses D1/D2 only).
        panic!("Standard_NotImplemented: SurfaceTool::DN");
    }
    fn u_resolution(s: &Surface<'a>, r3d: f64) -> f64 {
        s.my_surface().u_resolution(r3d)
    }
    fn v_resolution(s: &Surface<'a>, r3d: f64) -> f64 {
        s.my_surface().v_resolution(r3d)
    }
    fn get_type(s: &Surface<'a>) -> GeomAbsSurfaceType {
        s.get_type()
    }
    fn plane(s: &Surface<'a>) -> rcad_kernel::geom::Plane {
        s.plane()
    }
    fn cylinder(s: &Surface<'a>) -> rcad_kernel::geom::CylindricalSurface {
        s.cylinder()
    }
    fn cone(s: &Surface<'a>) -> rcad_kernel::geom::ConicalSurface {
        s.cone()
    }
    fn torus(s: &Surface<'a>) -> rcad_kernel::geom::ToroidalSurface {
        s.torus()
    }
    fn sphere(s: &Surface<'a>) -> rcad_kernel::geom::SphericalSurface {
        s.sphere()
    }
    fn axe_of_revolution(s: &Surface<'a>) -> (Point3, Vec3) {
        s.axis()
    }
    fn direction(_s: &Surface<'a>) -> Vec3 {
        panic!("Standard_NoSuchObject: SurfaceTool::Direction");
    }
    fn basis_curve(_s: &Surface<'a>) -> NullBasisCurve {
        panic!("Standard_NoSuchObject: SurfaceTool::BasisCurve");
    }
    fn basis_surface(_s: &Surface<'a>) -> NullBasisSurface {
        panic!("Standard_NoSuchObject: SurfaceTool::BasisSurface");
    }
    fn offset_value(_s: &Surface<'a>) -> f64 {
        panic!("Standard_NoSuchObject: SurfaceTool::OffsetValue");
    }
    fn nb_u_poles(s: &Surface<'a>) -> usize {
        s.my_surface().nb_u_poles()
    }
    fn nb_v_poles(s: &Surface<'a>) -> usize {
        s.my_surface().nb_v_poles()
    }
    fn nb_u_knots(s: &Surface<'a>) -> usize {
        s.my_surface().nb_u_knots()
    }
    fn nb_v_knots(s: &Surface<'a>) -> usize {
        s.my_surface().nb_v_knots()
    }
    fn u_degree(s: &Surface<'a>) -> usize {
        s.my_surface().u_degree()
    }
    fn v_degree(s: &Surface<'a>) -> usize {
        s.my_surface().v_degree()
    }
    fn bezier<'b>(_s: &'b Surface<'a>) -> &'b rcad_kernel::geom::BezierSurface {
        panic!("Standard_NoSuchObject: SurfaceTool::Bezier");
    }
    fn bspline<'b>(_s: &'b Surface<'a>) -> &'b rcad_kernel::geom::BSplineSurface {
        panic!("Standard_NoSuchObject: SurfaceTool::BSpline");
    }
}

impl<'a> PSurfaceTool for SurfaceTool<'a> {
    type Surface = Surface<'a>;
    fn value(s: &Surface<'a>, u: f64, v: f64) -> DVec3 {
        s.value(u, v)
    }
    fn d1(s: &Surface<'a>, u: f64, v: f64) -> (DVec3, DVec3, DVec3) {
        s.d1(u, v)
    }
    fn first_u_parameter(s: &Surface<'a>) -> f64 {
        s.first_u_parameter()
    }
    fn last_u_parameter(s: &Surface<'a>) -> f64 {
        s.last_u_parameter()
    }
    fn first_v_parameter(s: &Surface<'a>) -> f64 {
        s.first_v_parameter()
    }
    fn last_v_parameter(s: &Surface<'a>) -> f64 {
        s.last_v_parameter()
    }
    fn u_resolution(s: &Surface<'a>, r3d: f64) -> f64 {
        s.my_surface().u_resolution(r3d)
    }
    fn v_resolution(s: &Surface<'a>, r3d: f64) -> f64 {
        s.my_surface().v_resolution(r3d)
    }
}
