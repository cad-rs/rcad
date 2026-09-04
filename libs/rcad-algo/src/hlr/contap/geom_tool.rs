// OCCT GeomAdaptor pairing — the `GeomTool` is the HSurfaceTool (the OCCT
// `Adaptor3d_HSurfaceTool` template argument) over the
// [`GeomSurfaceAdapter`] (the `GeomAdaptor_Surface`).  It feeds the generic
// Adaptor3d TopolTool (`topalgo::adaptor3d::topol_tool::TopolTool`) so the
// Contap engines run on an adaptor surface, as in OCCT.

use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{
    BezierSurface, Circle3, ConicalSurface, CylindricalSurface, Hyperbola3, Line3, Parabola3,
    Plane, Point3, SphericalSurface, ToroidalSurface, Vec3,
};

use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::int_curve_surface::{
    Adaptor3dCurveBasis, Adaptor3dSurfaceBasis, HSurfaceTool,
};
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::hlr::contap::surface_adaptor::{GeomSurfaceAdapter, SurfaceAdapter};

/// The basis-curve adaptor placeholder (the restriction curves of the
/// TopolTool carry their own adaptors; the Geom basis curve of an
/// extrusion/revolution surface is not modelled here yet — Stage 3b).
#[derive(Debug, Clone, Copy, Default)]
pub struct GeomBasisCurve;

impl Adaptor3dCurveBasis for GeomBasisCurve {
    fn get_type(&self) -> CurveType {
        panic!("Standard_NoSuchObject: Geom basis curve");
    }
    fn value(&self, _u: f64) -> Point3 {
        panic!("Standard_NoSuchObject: Geom basis curve");
    }
    fn line(&self) -> Line3 {
        panic!("Standard_NoSuchObject: Geom basis curve");
    }
    fn parabola(&self) -> Parabola3 {
        panic!("Standard_NoSuchObject: Geom basis curve");
    }
    fn hyperbola(&self) -> Hyperbola3 {
        panic!("Standard_NoSuchObject: Geom basis curve");
    }
}

/// The basis-surface adaptor placeholder (offset surfaces not modelled).
#[derive(Debug, Clone, Copy, Default)]
pub struct GeomBasisSurface;

impl Adaptor3dSurfaceBasis<GeomBasisCurve> for GeomBasisSurface {
    fn get_type(&self) -> GeomAbsSurfaceType {
        panic!("Standard_NoSuchObject: Geom basis surface");
    }
    fn plane(&self) -> Plane {
        panic!("Standard_NoSuchObject: Geom basis surface");
    }
    fn cylinder(&self) -> CylindricalSurface {
        panic!("Standard_NoSuchObject: Geom basis surface");
    }
    fn cone(&self) -> ConicalSurface {
        panic!("Standard_NoSuchObject: Geom basis surface");
    }
    fn direction(&self) -> Vec3 {
        panic!("Standard_NoSuchObject: Geom basis surface");
    }
    fn basis_curve(&self) -> GeomBasisCurve {
        GeomBasisCurve
    }
}

/// OCCT `GeomAdaptor_Surface` + `Adaptor3d_HSurfaceTool` pairing over the
/// kernel surface.
#[derive(Debug, Clone, Copy, Default)]
pub struct GeomTool;

impl HSurfaceTool for GeomTool {
    type Surface = GeomSurfaceAdapter;
    type BasisCurve = GeomBasisCurve;
    type BasisSurface = GeomBasisSurface;

    fn first_u_parameter(s: &GeomSurfaceAdapter) -> f64 {
        s.first_u_parameter()
    }
    fn first_v_parameter(s: &GeomSurfaceAdapter) -> f64 {
        s.first_v_parameter()
    }
    fn last_u_parameter(s: &GeomSurfaceAdapter) -> f64 {
        s.last_u_parameter()
    }
    fn last_v_parameter(s: &GeomSurfaceAdapter) -> f64 {
        s.last_v_parameter()
    }
    fn nb_u_intervals(_s: &GeomSurfaceAdapter, _sh: GeomAbsShape) -> usize {
        1
    }
    fn nb_v_intervals(_s: &GeomSurfaceAdapter, _sh: GeomAbsShape) -> usize {
        1
    }
    fn u_intervals(s: &GeomSurfaceAdapter, tab: &mut [f64], _sh: GeomAbsShape) {
        tab[0] = s.first_u_parameter();
        tab[1] = s.last_u_parameter();
    }
    fn v_intervals(s: &GeomSurfaceAdapter, tab: &mut [f64], _sh: GeomAbsShape) {
        tab[0] = s.first_v_parameter();
        tab[1] = s.last_v_parameter();
    }
    fn u_trim(s: &GeomSurfaceAdapter, first: f64, last: f64, _tol: f64) -> GeomSurfaceAdapter {
        // OCCT GeomAdaptor_Surface::UTrim narrows the same-surface window.
        let d = s.uv_window();
        s.clone_with_window([first.min(last), last.max(first), d[2], d[3]])
    }
    fn v_trim(s: &GeomSurfaceAdapter, first: f64, last: f64, _tol: f64) -> GeomSurfaceAdapter {
        let d = s.uv_window();
        s.clone_with_window([d[0], d[1], first.min(last), last.max(first)])
    }
    fn is_u_closed(s: &GeomSurfaceAdapter) -> bool {
        matches!(
            s.get_type(),
            GeomAbsSurfaceType::Cylinder | GeomAbsSurfaceType::Torus
        )
    }
    fn is_v_closed(s: &GeomSurfaceAdapter) -> bool {
        s.get_type() == GeomAbsSurfaceType::Torus
    }
    fn is_u_periodic(s: &GeomSurfaceAdapter) -> bool {
        s.is_u_periodic()
    }
    fn u_period(s: &GeomSurfaceAdapter) -> f64 {
        s.u_period()
    }
    fn is_v_periodic(s: &GeomSurfaceAdapter) -> bool {
        s.is_v_periodic()
    }
    fn v_period(s: &GeomSurfaceAdapter) -> f64 {
        s.v_period()
    }
    fn value(s: &GeomSurfaceAdapter, u: f64, v: f64) -> Point3 {
        s.value(u, v)
    }
    fn d0(s: &GeomSurfaceAdapter, u: f64, v: f64) -> Point3 {
        s.value(u, v)
    }
    fn d1(s: &GeomSurfaceAdapter, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        s.d1(u, v)
    }
    fn d2(
        s: &GeomSurfaceAdapter,
        u: f64,
        v: f64,
    ) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        s.d2(u, v)
    }
    fn dn(s: &GeomSurfaceAdapter, u: f64, v: f64, _nu: usize, _nv: usize) -> Vec3 {
        // The analytic DN is not used by the Contap chain; the SurfaceEval
        // interface carries only D1/D2 (documented).
        let _ = (s, u, v);
        panic!("Standard_NotImplemented: GeomTool::DN");
    }
    fn u_resolution(s: &GeomSurfaceAdapter, r3d: f64) -> f64 {
        s.u_resolution(r3d)
    }
    fn v_resolution(s: &GeomSurfaceAdapter, r3d: f64) -> f64 {
        s.v_resolution(r3d)
    }
    fn get_type(s: &GeomSurfaceAdapter) -> GeomAbsSurfaceType {
        s.get_type()
    }
    fn plane(s: &GeomSurfaceAdapter) -> Plane {
        s.plane()
    }
    fn cylinder(s: &GeomSurfaceAdapter) -> CylindricalSurface {
        s.cylinder()
    }
    fn cone(s: &GeomSurfaceAdapter) -> ConicalSurface {
        s.cone()
    }
    fn torus(s: &GeomSurfaceAdapter) -> ToroidalSurface {
        s.torus()
    }
    fn sphere(s: &GeomSurfaceAdapter) -> SphericalSurface {
        s.sphere()
    }
    fn axe_of_revolution(s: &GeomSurfaceAdapter) -> (Point3, Vec3) {
        match s.get_type() {
            GeomAbsSurfaceType::Cylinder => {
                let c = s.cylinder();
                (c.origin, c.axis)
            }
            GeomAbsSurfaceType::Cone => {
                let c = s.cone();
                (c.apex, c.axis)
            }
            GeomAbsSurfaceType::Sphere => {
                let c = s.sphere();
                (c.center, c.axis)
            }
            GeomAbsSurfaceType::Torus => {
                let c = s.torus();
                (c.center, c.axis)
            }
            _ => panic!("Standard_NoSuchObject: GeomTool::AxeOfRevolution"),
        }
    }
    fn direction(s: &GeomSurfaceAdapter) -> Vec3 {
        match s.get_type() {
            GeomAbsSurfaceType::SurfaceOfRevolution | GeomAbsSurfaceType::SurfaceOfExtrusion => {
                // The extrusion/revolution direction is not modelled on the
                // kernel wrapper yet (Stage 3b).
                panic!("Standard_NoSuchObject: GeomTool::Direction");
            }
            _ => panic!("Standard_NoSuchObject: GeomTool::Direction"),
        }
    }
    fn basis_curve(_s: &GeomSurfaceAdapter) -> GeomBasisCurve {
        panic!("Standard_NoSuchObject: GeomTool::BasisCurve");
    }
    fn basis_surface(_s: &GeomSurfaceAdapter) -> GeomBasisSurface {
        panic!("Standard_NoSuchObject: GeomTool::BasisSurface");
    }
    fn offset_value(_s: &GeomSurfaceAdapter) -> f64 {
        panic!("Standard_NoSuchObject: GeomTool::OffsetValue");
    }
    fn nb_u_poles(s: &GeomSurfaceAdapter) -> usize {
        s.nb_u_poles()
    }
    fn nb_v_poles(s: &GeomSurfaceAdapter) -> usize {
        s.nb_v_poles()
    }
    fn nb_u_knots(s: &GeomSurfaceAdapter) -> usize {
        s.nb_u_knots()
    }
    fn nb_v_knots(s: &GeomSurfaceAdapter) -> usize {
        s.nb_v_knots()
    }
    fn u_degree(s: &GeomSurfaceAdapter) -> usize {
        s.u_degree()
    }
    fn v_degree(s: &GeomSurfaceAdapter) -> usize {
        s.v_degree()
    }
    fn bezier(_s: &GeomSurfaceAdapter) -> &BezierSurface {
        panic!("Standard_NoSuchObject: GeomTool::Bezier");
    }
    fn bspline(_s: &GeomSurfaceAdapter) -> &rcad_kernel::geom::BSplineSurface {
        panic!("Standard_NoSuchObject: GeomTool::BSpline");
    }
}

// The Circle3 reference of the OCCT gp_Circ reconstruction (Contap_ContAna)
// stays shape-checked through the cont_ana tests.
#[allow(unused)]
fn _parity(c: Circle3) -> Circle3 {
    c
}
