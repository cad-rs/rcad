//! The tool specialisations binding the bop `BRepAdaptorCurve` /
//! `BRepAdaptorSurface` into the IntCurveSurface_HInter assembly
//! (OCCT IntCurveSurface_HInter.cxx instantiates the Inter.pxx templates
//! with TheCurveTool = IntCurveSurface_TheHCurveTool over
//! `Adaptor3d_Curve` and TheSurfaceTool = Adaptor3d_HSurfaceTool over
//! `Adaptor3d_Surface`; the concrete adaptor here is the bop
//! BRepAdaptor pair).  Also hosts the [`HInterHost`] implementation for
//! [`IntCurveSurfaceHInter`] — the OCCT lambda bundle of HInter.cxx.

use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{
    BezierCurve3, BSplineCurve3, Circle3, ConicalSurface, Curve3, CurveEval, CylindricalSurface,
    Ellipse3, Hyperbola3, Line3, Parabola3, Plane, Point3, SphericalSurface, Surface3, SurfaceEval,
    ToroidalSurface, Vec3,
};
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::int_curve_surface::{
    Adaptor3dCurveBasis, Adaptor3dSurfaceBasis, HCurveTool, HSurfaceTool, HInterHost,
};
use crate::geomalgo::int_curve_surface::quad_curv_exact::TheQuadCurvExactHInter;
use crate::geomalgo::int_imp::{CurveTool3d, PSurfaceTool};
use crate::geomalgo::int_patch::GeomAbsSurfaceType;
use crate::geomalgo::int_patch::int_conic_quad::IntConicQuad;
use crate::geomalgo::int_curv_surf::{ThePolygonOfHInter, ThePolyhedronOfHInter};

use super::bean_face_intersector::{BRepAdaptorCurve, BRepAdaptorSurface, IntCurveSurfaceHInter};

/// OCCT Adaptor3d_Curve over the rcad `Curve3` — the basis-curve adaptor
/// type handed out by [`BRepAdaptorSurfaceTool::basis_curve`] (the
/// Revolution / Extrusion profile).
#[derive(Debug, Clone)]
pub struct BRepBasisCurve(pub Curve3);

impl Adaptor3dCurveBasis for BRepBasisCurve {
    fn get_type(&self) -> CurveType {
        curve_type_of(&self.0)
    }
    fn value(&self, u: f64) -> Point3 {
        CurveEval::point_at(&self.0, u)
    }
    fn line(&self) -> Line3 {
        match &self.0 {
            Curve3::Line(l) => l.clone(),
            _ => panic!("Standard_NoSuchObject: BasisCurve::Line"),
        }
    }
    fn parabola(&self) -> Parabola3 {
        match &self.0 {
            Curve3::Parabola(p) => *p,
            _ => panic!("Standard_NoSuchObject: BasisCurve::Parabola"),
        }
    }
    fn hyperbola(&self) -> Hyperbola3 {
        match &self.0 {
            Curve3::Hyperbola(h) => *h,
            _ => panic!("Standard_NoSuchObject: BasisCurve::Hyperbola"),
        }
    }
}

/// OCCT Adaptor3d_Surface over the rcad `Surface3` — the basis-surface
/// adaptor type handed out by [`BRepAdaptorSurfaceTool::basis_surface`]
/// (the Offset basis).
#[derive(Debug, Clone)]
pub struct BRepBasisSurface(pub Surface3);

impl Adaptor3dSurfaceBasis<BRepBasisCurve> for BRepBasisSurface {
    fn get_type(&self) -> GeomAbsSurfaceType {
        surface_type_of(&self.0)
    }
    fn plane(&self) -> Plane {
        match &self.0 {
            Surface3::Plane(p) => *p,
            _ => panic!("Standard_NoSuchObject: BasisSurface::Plane"),
        }
    }
    fn cylinder(&self) -> CylindricalSurface {
        match &self.0 {
            Surface3::Cylinder(c) => *c,
            _ => panic!("Standard_NoSuchObject: BasisSurface::Cylinder"),
        }
    }
    fn cone(&self) -> ConicalSurface {
        match &self.0 {
            Surface3::Cone(c) => *c,
            _ => panic!("Standard_NoSuchObject: BasisSurface::Cone"),
        }
    }
    fn direction(&self) -> Vec3 {
        match &self.0 {
            Surface3::LinearExtrusion(e) => e.direction,
            _ => panic!("Standard_NoSuchObject: BasisSurface::Direction"),
        }
    }
    fn basis_curve(&self) -> BRepBasisCurve {
        match &self.0 {
            Surface3::Revolution(r) => BRepBasisCurve((*r.profile).clone()),
            Surface3::LinearExtrusion(e) => BRepBasisCurve((*e.profile).clone()),
            _ => panic!("Standard_NoSuchObject: BasisSurface::BasisCurve"),
        }
    }
}

/// The rcad `Curve3` classification (OCCT GeomAbs_CurveType).
fn curve_type_of(c: &Curve3) -> CurveType {
    match c {
        Curve3::Line(_) => CurveType::Line,
        Curve3::Circle(_) => CurveType::Circle,
        Curve3::Ellipse(_) => CurveType::Ellipse,
        Curve3::Hyperbola(_) => CurveType::Hyperbola,
        Curve3::Parabola(_) => CurveType::Parabola,
        Curve3::Bezier(_) => CurveType::Bezier,
        Curve3::BSpline(_) => CurveType::BSpline,
        _ => CurveType::Other,
    }
}

/// The rcad `Surface3` classification (OCCT GeomAbs_SurfaceType).
fn surface_type_of(s: &Surface3) -> GeomAbsSurfaceType {
    match s {
        Surface3::Plane(_) => GeomAbsSurfaceType::Plane,
        Surface3::Cylinder(_) => GeomAbsSurfaceType::Cylinder,
        Surface3::Cone(_) => GeomAbsSurfaceType::Cone,
        Surface3::Sphere(_) => GeomAbsSurfaceType::Sphere,
        Surface3::Torus(_) => GeomAbsSurfaceType::Torus,
        Surface3::Bezier(_) => GeomAbsSurfaceType::BezierSurface,
        Surface3::BSpline(_) => GeomAbsSurfaceType::BSplineSurface,
        Surface3::Revolution(_) => GeomAbsSurfaceType::SurfaceOfRevolution,
        Surface3::LinearExtrusion(_) => GeomAbsSurfaceType::SurfaceOfExtrusion,
        Surface3::Offset(_) => GeomAbsSurfaceType::OffsetSurface,
        _ => GeomAbsSurfaceType::OtherSurface,
    }
}

/// OCCT IntCurveSurface_TheHCurveTool over `BRepAdaptor_Curve`.
pub struct BRepAdaptorCurveTool;

impl HCurveTool for BRepAdaptorCurveTool {
    type Curve = BRepAdaptorCurve;

    fn first_parameter(c: &BRepAdaptorCurve) -> f64 {
        c.first_parameter()
    }
    fn last_parameter(c: &BRepAdaptorCurve) -> f64 {
        c.last_parameter()
    }
    /// OCCT BRepAdaptor_Curve::Continuity — the underlying geometry is C
    /// infinity for the analytic types (GeomAbs_CN); a C2-rated BSpline is
    /// likewise never split by the C2 interval queries.
    fn continuity(_c: &BRepAdaptorCurve) -> GeomAbsShape {
        GeomAbsShape::CN
    }
    fn nb_intervals(_c: &BRepAdaptorCurve, _s: GeomAbsShape) -> usize {
        1
    }
    fn intervals(c: &BRepAdaptorCurve, t: &mut [f64], _s: GeomAbsShape) {
        t[0] = c.first_parameter();
        t[1] = c.last_parameter();
    }
    fn is_closed(c: &BRepAdaptorCurve) -> bool {
        CurveEval::point_at(c.curve(), c.first_parameter())
            .distance(CurveEval::point_at(c.curve(), c.last_parameter()))
            < rcad_kernel::precision::CONFUSION
    }
    fn is_periodic(c: &BRepAdaptorCurve) -> bool {
        c.is_periodic()
    }
    fn period(c: &BRepAdaptorCurve) -> f64 {
        c.period()
    }
    fn value(c: &BRepAdaptorCurve, u: f64) -> Point3 {
        CurveEval::point_at(c.curve(), u)
    }
    fn d0(c: &BRepAdaptorCurve, u: f64) -> Point3 {
        CurveEval::point_at(c.curve(), u)
    }
    fn d1(c: &BRepAdaptorCurve, u: f64) -> (Point3, Vec3) {
        (
            CurveEval::point_at(c.curve(), u),
            CurveEval::derivative_at(c.curve(), u),
        )
    }
    fn d2(c: &BRepAdaptorCurve, u: f64) -> (Point3, Vec3, Vec3) {
        (
            CurveEval::point_at(c.curve(), u),
            CurveEval::derivative_at(c.curve(), u),
            CurveEval::derivative2_at(c.curve(), u),
        )
    }
    fn d3(_c: &BRepAdaptorCurve, _u: f64) -> (Point3, Vec3, Vec3, Vec3) {
        panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::D3");
    }
    fn dn(c: &BRepAdaptorCurve, u: f64, n: usize) -> Vec3 {
        let _ = (c, u);
        match n {
            1 => CurveEval::derivative_at(c.curve(), u),
            2 => CurveEval::derivative2_at(c.curve(), u),
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::DN"),
        }
    }
    fn resolution(c: &BRepAdaptorCurve, r3d: f64) -> f64 {
        c.resolution(r3d)
    }
    fn get_type(c: &BRepAdaptorCurve) -> CurveType {
        curve_type_of(c.curve())
    }
    fn line(c: &BRepAdaptorCurve) -> Line3 {
        match c.curve() {
            Curve3::Line(l) => *l,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::Line"),
        }
    }
    fn circle(c: &BRepAdaptorCurve) -> Circle3 {
        match c.curve() {
            Curve3::Circle(c2) => *c2,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::Circle"),
        }
    }
    fn ellipse(c: &BRepAdaptorCurve) -> Ellipse3 {
        match c.curve() {
            Curve3::Ellipse(e) => *e,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::Ellipse"),
        }
    }
    fn hyperbola(c: &BRepAdaptorCurve) -> Hyperbola3 {
        match c.curve() {
            Curve3::Hyperbola(h) => *h,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::Hyperbola"),
        }
    }
    fn parabola(c: &BRepAdaptorCurve) -> Parabola3 {
        match c.curve() {
            Curve3::Parabola(p) => *p,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::Parabola"),
        }
    }
    fn bezier(c: &BRepAdaptorCurve) -> &BezierCurve3 {
        match c.curve() {
            Curve3::Bezier(b) => b,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::Bezier"),
        }
    }
    fn bspline(c: &BRepAdaptorCurve) -> &BSplineCurve3 {
        match c.curve() {
            Curve3::BSpline(b) => b,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorCurveTool::BSpline"),
        }
    }
}

impl CurveTool3d for BRepAdaptorCurveTool {
    type Curve = BRepAdaptorCurve;
    fn value(c: &BRepAdaptorCurve, u: f64) -> glam::DVec3 {
        <Self as HCurveTool>::value(c, u)
    }
    fn d1(c: &BRepAdaptorCurve, u: f64) -> (glam::DVec3, glam::DVec3) {
        <Self as HCurveTool>::d1(c, u)
    }
    fn first_parameter(c: &BRepAdaptorCurve) -> f64 {
        <Self as HCurveTool>::first_parameter(c)
    }
    fn last_parameter(c: &BRepAdaptorCurve) -> f64 {
        <Self as HCurveTool>::last_parameter(c)
    }
    fn resolution(c: &BRepAdaptorCurve, r3d: f64) -> f64 {
        <Self as HCurveTool>::resolution(c, r3d)
    }
}

/// OCCT Adaptor3d_HSurfaceTool over `BRepAdaptor_Surface`.
pub struct BRepAdaptorSurfaceTool;

impl HSurfaceTool for BRepAdaptorSurfaceTool {
    type Surface = BRepAdaptorSurface;
    type BasisCurve = BRepBasisCurve;
    type BasisSurface = BRepBasisSurface;

    fn first_u_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.first_u_parameter()
    }
    fn first_v_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.first_v_parameter()
    }
    fn last_u_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.last_u_parameter()
    }
    fn last_v_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.last_v_parameter()
    }
    /// OCCT BRepAdaptor_Surface::NbUIntervals(GeomAbs_C2) — 1 for the
    /// analytic surfaces; the BSpline reports one interval per distinct
    /// knot span.
    fn nb_u_intervals(s: &BRepAdaptorSurface, _sh: GeomAbsShape) -> usize {
        match s.surface() {
            Surface3::BSpline(bs) => {
                crate::geomalgo::int_curve_surface::distinct_knots(&bs.knots_u)
                    .len()
                    .saturating_sub(1)
                    .max(1)
            }
            _ => 1,
        }
    }
    fn nb_v_intervals(s: &BRepAdaptorSurface, _sh: GeomAbsShape) -> usize {
        match s.surface() {
            Surface3::BSpline(bs) => {
                crate::geomalgo::int_curve_surface::distinct_knots(&bs.knots_v)
                    .len()
                    .saturating_sub(1)
                    .max(1)
            }
            _ => 1,
        }
    }
    fn u_intervals(s: &BRepAdaptorSurface, tab: &mut [f64], _sh: GeomAbsShape) {
        match s.surface() {
            Surface3::BSpline(bs) => {
                let knots = crate::geomalgo::int_curve_surface::distinct_knots(&bs.knots_u);
                for (i, k) in knots.iter().enumerate() {
                    tab[i] = *k;
                }
            }
            _ => {
                tab[0] = s.first_u_parameter();
                tab[1] = s.last_u_parameter();
            }
        }
    }
    fn v_intervals(s: &BRepAdaptorSurface, tab: &mut [f64], _sh: GeomAbsShape) {
        match s.surface() {
            Surface3::BSpline(bs) => {
                let knots = crate::geomalgo::int_curve_surface::distinct_knots(&bs.knots_v);
                for (i, k) in knots.iter().enumerate() {
                    tab[i] = *k;
                }
            }
            _ => {
                tab[0] = s.first_v_parameter();
                tab[1] = s.last_v_parameter();
            }
        }
    }
    /// OCCT BRepAdaptor_Surface::UTrim — the same surface, narrowed U window.
    fn u_trim(s: &BRepAdaptorSurface, first: f64, last: f64, _tol: f64) -> BRepAdaptorSurface {
        BRepAdaptorSurface::with_uv_bounds(
            s.surface().clone(),
            [first, last, s.first_v_parameter(), s.last_v_parameter()],
        )
    }
    fn v_trim(s: &BRepAdaptorSurface, first: f64, last: f64, _tol: f64) -> BRepAdaptorSurface {
        BRepAdaptorSurface::with_uv_bounds(
            s.surface().clone(),
            [s.first_u_parameter(), s.last_u_parameter(), first, last],
        )
    }
    fn is_u_closed(s: &BRepAdaptorSurface) -> bool {
        SurfaceEval::is_u_closed(s.surface())
    }
    fn is_v_closed(s: &BRepAdaptorSurface) -> bool {
        SurfaceEval::is_v_closed(s.surface())
    }
    fn is_u_periodic(s: &BRepAdaptorSurface) -> bool {
        s.is_u_periodic()
    }
    fn u_period(s: &BRepAdaptorSurface) -> f64 {
        s.u_period()
    }
    fn is_v_periodic(s: &BRepAdaptorSurface) -> bool {
        s.is_v_periodic()
    }
    fn v_period(s: &BRepAdaptorSurface) -> f64 {
        s.v_period()
    }
    fn value(s: &BRepAdaptorSurface, u: f64, v: f64) -> Point3 {
        SurfaceEval::point_at(s.surface(), u, v)
    }
    fn d0(s: &BRepAdaptorSurface, u: f64, v: f64) -> Point3 {
        SurfaceEval::point_at(s.surface(), u, v)
    }
    fn d1(s: &BRepAdaptorSurface, u: f64, v: f64) -> (Point3, Vec3, Vec3) {
        SurfaceEval::derivatives(s.surface(), u, v)
    }
    fn d2(
        s: &BRepAdaptorSurface,
        u: f64,
        v: f64,
    ) -> (Point3, Vec3, Vec3, Vec3, Vec3, Vec3) {
        // kernel SurfaceEval::derivatives2 = (P, D1U, D1V, D2UU, D2UV, D2VV);
        // the OCCT D2 order is (P, D1U, D1V, D2U, D2V, D2UV).
        let (p, du, dv, duu, duv, dvv) = SurfaceEval::derivatives2(s.surface(), u, v);
        (p, du, dv, duu, dvv, duv)
    }
    fn dn(_s: &BRepAdaptorSurface, _u: f64, _v: f64, _nu: usize, _nv: usize) -> Vec3 {
        panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::DN");
    }
    /// OCCT GeomAdaptor_Surface::UResolution (GeomAdaptor_Surface.cxx
    /// L1818-1898) — the per-type parameter tolerance.
    fn u_resolution(s: &BRepAdaptorSurface, r3d: f64) -> f64 {
        match s.surface() {
            Surface3::Plane(_) => r3d,
            Surface3::Cylinder(c) => {
                let r = c.radius;
                if r > rcad_kernel::precision::CONFUSION {
                    let res = r3d / (2.0 * r);
                    if res <= 1.0 {
                        2.0 * res.asin()
                    } else {
                        2.0 * std::f64::consts::PI
                    }
                } else {
                    0.0
                }
            }
            Surface3::Sphere(sp) => {
                let r = sp.radius;
                if r > rcad_kernel::precision::CONFUSION {
                    let res = r3d / (2.0 * r);
                    if res <= 1.0 {
                        2.0 * res.asin()
                    } else {
                        2.0 * std::f64::consts::PI
                    }
                } else {
                    0.0
                }
            }
            Surface3::Cone(c) => {
                // OCCT L1848-1864: unbounded V domain -> Precision::Parametric;
                // otherwise R3d / max(VIso radius at VFirst, VIso at VLast).
                let r1 = cone_iso_radius(c, s.first_v_parameter());
                let r2 = cone_iso_radius(c, s.last_v_parameter());
                let r = r1.max(r2);
                if r > rcad_kernel::precision::CONFUSION {
                    r3d / r
                } else {
                    0.0
                }
            }
            Surface3::Torus(t) => {
                let r = t.major_radius + t.minor_radius;
                if r > rcad_kernel::precision::CONFUSION {
                    let res = r3d / (2.0 * r);
                    if res <= 1.0 {
                        2.0 * res.asin()
                    } else {
                        2.0 * std::f64::consts::PI
                    }
                } else {
                    0.0
                }
            }
            _ => r3d * 0.01, // OCCT default: Precision::Parametric(R3d)
        }
    }
    /// OCCT GeomAdaptor_Surface::VResolution (GeomAdaptor_Surface.cxx
    /// L1900-1958).
    fn v_resolution(s: &BRepAdaptorSurface, r3d: f64) -> f64 {
        match s.surface() {
            Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Cone(_) => r3d,
            Surface3::Sphere(sp) => {
                let r = sp.radius;
                if r > rcad_kernel::precision::CONFUSION {
                    let res = r3d / (2.0 * r);
                    if res <= 1.0 {
                        2.0 * res.asin()
                    } else {
                        2.0 * std::f64::consts::PI
                    }
                } else {
                    0.0
                }
            }
            Surface3::Torus(t) => {
                // OCCT L1909-1920: R = MinorRadius.
                let r = t.minor_radius;
                if r > rcad_kernel::precision::CONFUSION {
                    let res = r3d / (2.0 * r);
                    if res <= 1.0 {
                        2.0 * res.asin()
                    } else {
                        2.0 * std::f64::consts::PI
                    }
                } else {
                    0.0
                }
            }
            _ => r3d * 0.01,
        }
    }
    fn get_type(s: &BRepAdaptorSurface) -> GeomAbsSurfaceType {
        surface_type_of(s.surface())
    }
    fn plane(s: &BRepAdaptorSurface) -> Plane {
        s.plane()
    }
    fn cylinder(s: &BRepAdaptorSurface) -> CylindricalSurface {
        s.cylinder()
    }
    fn cone(s: &BRepAdaptorSurface) -> ConicalSurface {
        s.cone()
    }
    fn torus(s: &BRepAdaptorSurface) -> ToroidalSurface {
        s.torus()
    }
    fn sphere(s: &BRepAdaptorSurface) -> SphericalSurface {
        s.sphere()
    }
    fn axe_of_revolution(s: &BRepAdaptorSurface) -> (Point3, Vec3) {
        match s.surface() {
            Surface3::Revolution(r) => (r.axis_origin, r.axis_dir),
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::AxeOfRevolution"),
        }
    }
    fn direction(s: &BRepAdaptorSurface) -> Vec3 {
        match s.surface() {
            Surface3::LinearExtrusion(e) => e.direction,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::Direction"),
        }
    }
    fn basis_curve(s: &BRepAdaptorSurface) -> BRepBasisCurve {
        match s.surface() {
            Surface3::Revolution(r) => BRepBasisCurve((*r.profile).clone()),
            Surface3::LinearExtrusion(e) => BRepBasisCurve((*e.profile).clone()),
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::BasisCurve"),
        }
    }
    fn basis_surface(s: &BRepAdaptorSurface) -> BRepBasisSurface {
        match s.surface() {
            Surface3::Offset(o) => BRepBasisSurface((*o.basis).clone()),
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::BasisSurface"),
        }
    }
    fn offset_value(s: &BRepAdaptorSurface) -> f64 {
        match s.surface() {
            Surface3::Offset(o) => o.offset_distance,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::OffsetValue"),
        }
    }
    fn nb_u_poles(s: &BRepAdaptorSurface) -> usize {
        match s.surface() {
            Surface3::BSpline(bs) => bs.control_points.len(),
            Surface3::Bezier(bz) => bz.control_points.len(),
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::NbUPoles"),
        }
    }
    fn nb_v_poles(s: &BRepAdaptorSurface) -> usize {
        match s.surface() {
            Surface3::BSpline(bs) => bs.control_points.first().map_or(0, |r| r.len()),
            Surface3::Bezier(bz) => bz.control_points.first().map_or(0, |r| r.len()),
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::NbVPoles"),
        }
    }
    fn nb_u_knots(s: &BRepAdaptorSurface) -> usize {
        match s.surface() {
            Surface3::BSpline(bs) => {
                crate::geomalgo::int_curve_surface::distinct_knots(&bs.knots_u).len()
            }
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::NbUKnots"),
        }
    }
    fn nb_v_knots(s: &BRepAdaptorSurface) -> usize {
        match s.surface() {
            Surface3::BSpline(bs) => {
                crate::geomalgo::int_curve_surface::distinct_knots(&bs.knots_v).len()
            }
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::NbVKnots"),
        }
    }
    fn u_degree(s: &BRepAdaptorSurface) -> usize {
        s.u_degree()
    }
    fn v_degree(s: &BRepAdaptorSurface) -> usize {
        s.v_degree()
    }
    fn bezier(s: &BRepAdaptorSurface) -> &rcad_kernel::geom::BezierSurface {
        match s.surface() {
            Surface3::Bezier(b) => b,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::Bezier"),
        }
    }
    fn bspline(s: &BRepAdaptorSurface) -> &rcad_kernel::geom::BSplineSurface {
        match s.surface() {
            Surface3::BSpline(b) => b,
            _ => panic!("Standard_NoSuchObject: BRepAdaptorSurfaceTool::BSpline"),
        }
    }
}

impl PSurfaceTool for BRepAdaptorSurfaceTool {
    type Surface = BRepAdaptorSurface;
    fn value(s: &BRepAdaptorSurface, u: f64, v: f64) -> glam::DVec3 {
        SurfaceEval::point_at(s.surface(), u, v)
    }
    fn d1(s: &BRepAdaptorSurface, u: f64, v: f64) -> (glam::DVec3, glam::DVec3, glam::DVec3) {
        SurfaceEval::derivatives(s.surface(), u, v)
    }
    fn first_u_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.first_u_parameter()
    }
    fn last_u_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.last_u_parameter()
    }
    fn first_v_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.first_v_parameter()
    }
    fn last_v_parameter(s: &BRepAdaptorSurface) -> f64 {
        s.last_v_parameter()
    }
    fn u_resolution(s: &BRepAdaptorSurface, r3d: f64) -> f64 {
        <Self as HSurfaceTool>::u_resolution(s, r3d)
    }
    fn v_resolution(s: &BRepAdaptorSurface, r3d: f64) -> f64 {
        <Self as HSurfaceTool>::v_resolution(s, r3d)
    }
}

/// OCCT GeomAdaptor_Surface::UResolution cone helper — the radius of the
/// v-isoparametric circle (the distance of the surface point at u = 0 from
/// the cone axis).
fn cone_iso_radius(c: &ConicalSurface, v: f64) -> f64 {
    let p = SurfaceEval::point_at(&Surface3::Cone(*c), 0.0, v);
    (p - c.apex).cross(c.axis).length()
}

/// OCCT IntCurveSurface_HInter::Perform(Curve, Surface) (HInter.cxx
/// L106-116) — the Inter.pxx entry: decompose the surface by C2 intervals
/// and run PerformBounds per interval.
pub(crate) fn run_perform(
    curve: &BRepAdaptorCurve,
    surface: &BRepAdaptorSurface,
    host: &mut IntCurveSurfaceHInter,
) {
    crate::geomalgo::int_curve_surface::inter_impl::perform::<
        BRepAdaptorCurve,
        BRepAdaptorCurveTool,
        BRepAdaptorSurface,
        BRepAdaptorSurfaceTool,
        IntCurveSurfaceHInter,
    >(curve, surface, host);
}

/// The OCCT lambda bundle of IntCurveSurface_HInter.cxx (the `[this](...)`
/// callbacks the HInter member functions pass into the Inter.pxx
/// templates), forwarding into the shared IntCurveSurface_InterImpl engine
/// exactly as the OCCT member functions do.
impl HInterHost<BRepAdaptorCurve, BRepAdaptorCurveTool, BRepAdaptorSurface, BRepAdaptorSurfaceTool>
    for IntCurveSurfaceHInter
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

    fn append(
        &mut self,
        pt: &crate::geomalgo::int_curve_surface::IntersectionPoint,
    ) {
        self.base.append(pt);
    }

    fn perform_bounds(
        &mut self,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::perform_bounds::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(c, s, u1, v1, u2, v2, self);
    }

    fn perform_conic_line(
        &mut self,
        line: &Line3,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::perform_conic_surf_line::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(line, c, s, u1, v1, u2, v2, self);
    }

    fn perform_conic_circle(
        &mut self,
        circle: &Circle3,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::perform_conic_surf_circle::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(circle, c, s, u1, v1, u2, v2, self);
    }

    fn perform_conic_ellipse(
        &mut self,
        ellipse: &Ellipse3,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::perform_conic_surf_ellipse::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(ellipse, c, s, u1, v1, u2, v2, self);
    }

    fn perform_conic_parabola(
        &mut self,
        parab: &Parabola3,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::perform_conic_surf_parabola::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(parab, c, s, u1, v1, u2, v2, self);
    }

    fn perform_conic_hyperbola(
        &mut self,
        hyper: &Hyperbola3,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::perform_conic_surf_hyperbola::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(hyper, c, s, u1, v1, u2, v2, self);
    }

    fn perform_polygon_polyhedron(
        &mut self,
        c: &BRepAdaptorCurve,
        p: &ThePolygonOfHInter,
        s: &BRepAdaptorSurface,
        ph: &ThePolyhedronOfHInter,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::perform_polygon_polyhedron::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(c, p, s, ph, self);
    }

    fn internal_perform(
        &mut self,
        c: &BRepAdaptorCurve,
        p: &ThePolygonOfHInter,
        s: &BRepAdaptorSurface,
        ph: &ThePolyhedronOfHInter,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::internal_perform::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(c, p, s, ph, u1, v1, u2, v2, self);
    }

    fn internal_perform_bsb(
        &mut self,
        c: &BRepAdaptorCurve,
        p: &ThePolygonOfHInter,
        s: &BRepAdaptorSurface,
        ph: &ThePolyhedronOfHInter,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        bsb: &mut rcad_kernel::math::bnd::BoundSortBox,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::internal_perform_bsb::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(c, p, s, ph, u1, v1, u2, v2, bsb, self);
    }

    fn internal_perform_bounds(
        &mut self,
        c: &BRepAdaptorCurve,
        p: &ThePolygonOfHInter,
        s: &BRepAdaptorSurface,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::internal_perform_polygon_bounds::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            TheQuadCurvExactHInter<BRepAdaptorSurface, BRepAdaptorSurfaceTool, BRepAdaptorCurve, BRepAdaptorCurveTool>,
            Self,
        >(c, p, s, u1, v1, u2, v2, self);
    }

    fn internal_perform_curve_quadric(
        &mut self,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::internal_perform_curve_quadric::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            TheQuadCurvExactHInter<BRepAdaptorSurface, BRepAdaptorSurfaceTool, BRepAdaptorCurve, BRepAdaptorCurveTool>,
            Self,
        >(c, s, self);
    }

    fn append_int_ana(
        &mut self,
        c: &BRepAdaptorCurve,
        s: &BRepAdaptorSurface,
        ana: &IntConicQuad,
    ) {
        crate::geomalgo::int_curve_surface::inter_impl::append_int_ana::<
            BRepAdaptorCurve,
            BRepAdaptorCurveTool,
            BRepAdaptorSurface,
            BRepAdaptorSurfaceTool,
            Self,
        >(c, s, ana, self);
    }
}
