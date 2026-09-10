//! Adaptor interface encodings shared by the ProjLib / Approx translations.
//!
//! OCCT passes `occ::handle<Adaptor3d_Surface>`, `occ::handle<Adaptor3d_Curve>`
//! and `occ::handle<Adaptor2d_Curve2d>` (TKG3d / TKG2d) through the whole
//! ProjLib and Approx packages and reaches them either directly or via the
//! static `Adaptor3d_HSurfaceTool` facade.  The rcad encoding keeps that
//! shape:
//! - each abstract adaptor class becomes an object-safe trait whose methods
//!   mirror the OCCT virtuals used by the consuming translations,
//! - `occ::handle<T>` maps to `Arc<dyn Trait>` (shared ownership; a handle
//!   copy in OCCT is shallow, `Arc::clone` here),
//! - `Adaptor3d_HSurfaceTool` one-liners collapse into direct trait calls
//!   (same encoding as `hlr::contap::surface_adaptor`).
//!
//! [`CurveOnSurface`] is the `Adaptor3d_CurveOnSurface` encoding.  Its
//! evaluation (Value/D1/D2 compositions) follows the OCCT definition; the
//! interval machinery is a GAP carrier (staged) — see the struct docs.

use std::sync::Arc;

use glam::{DVec2, DVec3};

use crate::geom::Curve2dEval;

use crate::math::GeomAbsShape;

use super::CurveType;

// ---------------------------------------------------------------------------
// GeomAbs_SurfaceType — mirrors GeomAbs_SurfaceType
// ---------------------------------------------------------------------------

/// Type of the adapted surface.
///
/// OCCT: `GeomAbs_SurfaceType`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomAbsSurfaceType {
    Plane,
    Cylinder,
    Cone,
    Sphere,
    Torus,
    BezierSurface,
    BSplineSurface,
    SurfaceOfRevolution,
    SurfaceOfExtrusion,
    OffsetSurface,
    OtherSurface,
}

// ---------------------------------------------------------------------------
// Adaptor2d_Curve2d — the 2D adaptor instance interface
// ---------------------------------------------------------------------------

/// OCCT Adaptor2d_Curve2d (TKG2d) — the instance interface with the methods
/// the ProjLib / Approx translations use.  `ProjLib_CompProjectedCurve`
/// inherits from this class in OCCT; the geomalgo translation implements
/// this trait for it.
pub trait Adaptor2dCurve2d {
    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64;
    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64;
    /// OCCT Value(U).
    fn value(&self, u: f64) -> DVec2;
    /// OCCT D0(U, P).
    fn d0(&self, u: f64) -> DVec2;
    /// OCCT D1(U, P, V).
    fn d1(&self, u: f64) -> (DVec2, DVec2);
    /// OCCT D2(U, P, V1, V2).
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2);
    /// OCCT Continuity().
    fn continuity(&self) -> GeomAbsShape;
    /// OCCT GetType().
    fn get_type(&self) -> CurveType;
    /// OCCT Line() — valid when GetType() == Line.
    fn line(&self) -> crate::geom::Line2d;
    /// OCCT BSpline() — the underlying bspline when GetType() == BSpline
    /// (the OCCT downcast to Geom2d_BSplineCurve; None models the null
    /// handle of the non-bspline case).
    fn bspline(&self) -> Option<crate::geom::BSplineCurve2>;
    /// OCCT Bezier() — the underlying bezier when GetType() == Bezier
    /// (the OCCT downcast to Geom2d_BezierCurve; None models the null
    /// handle of the non-bezier case).
    fn bezier(&self) -> Option<crate::geom::BezierCurve2>;
    /// OCCT Trim(FirstParam, LastParam, Tol) — returns a curve equivalent of
    /// `self` restricted to [First, Last].
    fn trim(&self, first_param: f64, last_param: f64, tol: f64) -> Arc<dyn Adaptor2dCurve2d>;
}

/// OCCT `occ::handle<Adaptor2d_Curve2d>` — shared ownership.
pub type Curve2dHandle = Arc<dyn Adaptor2dCurve2d>;

// ---------------------------------------------------------------------------
// Adaptor3d_Curve — the 3D adaptor instance interface
// ---------------------------------------------------------------------------

/// OCCT Adaptor3d_Curve (TKG3d) — the instance interface with the methods the
/// ProjLib translations use (ProjLib_CompProjectedCurve.cxx, ProjLib_PrjFunc).
pub trait Adaptor3dCurve {
    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64;
    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64;
    /// OCCT Value(U) / D0(U, P).
    fn value(&self, u: f64) -> DVec3;
    /// OCCT D1(U, P, V1).
    fn d1(&self, u: f64) -> (DVec3, DVec3);
    /// OCCT D2(U, P, V1, V2).
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3);
    /// OCCT Continuity().
    fn continuity(&self) -> GeomAbsShape;
    /// OCCT NbIntervals(S).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize;
    /// OCCT Intervals(T, S) — the S-discontinuity parameters.
    fn intervals(&self, s: GeomAbsShape) -> Vec<f64>;
    /// OCCT Trim(First, Last, Tol).
    fn trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurve>;
}

/// OCCT `occ::handle<Adaptor3d_Curve>` — shared ownership.
pub type CurveHandle = Arc<dyn Adaptor3dCurve>;

// ---------------------------------------------------------------------------
// Adaptor3d_Surface — the surface adaptor instance interface
// ---------------------------------------------------------------------------

/// OCCT Adaptor3d_Surface (TKG3d) — the instance interface with the methods
/// the ProjLib / Approx translations use (a subset of
/// Adaptor3d_HSurfaceTool.hxx).  Derivative out-parameters map to tuple slots
/// in OCCT declaration order.
pub trait Adaptor3dSurface {
    /// OCCT FirstUParameter().
    fn first_u_parameter(&self) -> f64;
    /// OCCT LastUParameter().
    fn last_u_parameter(&self) -> f64;
    /// OCCT FirstVParameter().
    fn first_v_parameter(&self) -> f64;
    /// OCCT LastVParameter().
    fn last_v_parameter(&self) -> f64;
    /// OCCT Value(U, V) / D0(U, V, P).
    fn value(&self, u: f64, v: f64) -> DVec3;
    /// OCCT D1(U, V, P, D1U, D1V).
    fn d1(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3);
    /// OCCT D2(U, V, P, D1U, D1V, D2U, D2V, D2UV).
    fn d2(&self, u: f64, v: f64) -> (DVec3, DVec3, DVec3, DVec3, DVec3, DVec3);
    /// OCCT D3(U, V, P, D1U, D1V, D2U, D2V, D2UV, D3U, D3V, D3UUV, D3UVV).
    #[allow(clippy::type_complexity)]
    fn d3(
        &self,
        u: f64,
        v: f64,
    ) -> (
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
        DVec3,
    );
    /// OCCT UResolution(R3d).
    fn u_resolution(&self, r3d: f64) -> f64;
    /// OCCT VResolution(R3d).
    fn v_resolution(&self, r3d: f64) -> f64;
    /// OCCT GetType().
    fn get_type(&self) -> GeomAbsSurfaceType;
    /// OCCT IsUPeriodic().
    fn is_u_periodic(&self) -> bool;
    /// OCCT UPeriod().
    fn u_period(&self) -> f64;
    /// OCCT IsVPeriodic().
    fn is_v_periodic(&self) -> bool;
    /// OCCT VPeriod().
    fn v_period(&self) -> f64;
    /// OCCT UContinuity().
    fn u_continuity(&self) -> GeomAbsShape;
    /// OCCT VContinuity().
    fn v_continuity(&self) -> GeomAbsShape;
    /// OCCT NbUIntervals(S).
    fn nb_u_intervals(&self, s: GeomAbsShape) -> usize;
    /// OCCT NbVIntervals(S).
    fn nb_v_intervals(&self, s: GeomAbsShape) -> usize;
    /// OCCT UIntervals(T, S).
    fn u_intervals(&self, s: GeomAbsShape) -> Vec<f64>;
    /// OCCT VIntervals(T, S).
    fn v_intervals(&self, s: GeomAbsShape) -> Vec<f64>;
    /// OCCT UTrim(U1, U2, Eps) — the surface restricted in U.
    /// GAP (staged): the Geom_RectangularTrimmedSurface adaptor family is not
    /// translated yet.
    fn u_trim(&self, _u1: f64, _u2: f64, _eps: f64) -> Arc<dyn Adaptor3dSurface> {
        panic!("GAP: Adaptor3d_Surface::UTrim not translated");
    }
    /// OCCT VTrim(V1, V2, Eps) — the surface restricted in V.
    /// GAP (staged): see [`Adaptor3dSurface::u_trim`].
    fn v_trim(&self, _v1: f64, _v2: f64, _eps: f64) -> Arc<dyn Adaptor3dSurface> {
        panic!("GAP: Adaptor3d_Surface::VTrim not translated");
    }
    /// OCCT ShallowCopy().
    fn shallow_copy(&self) -> Arc<dyn Adaptor3dSurface>;
    /// The kernel surface when the adaptor wraps one — the routing bridge for
    /// the landed Surface3-based Extrema engine (same bridge as
    /// `hlr::contap::surface_adaptor::SurfaceAdapter::kernel_surface`).
    fn kernel_surface(&self) -> Option<&crate::geom::Surface3>;
}

/// OCCT `occ::handle<Adaptor3d_Surface>` — shared ownership.
pub type SurfaceHandle = Arc<dyn Adaptor3dSurface>;

// ---------------------------------------------------------------------------
// Geom2dAdaptor_Curve — the 2D adaptor over a kernel curve
// ---------------------------------------------------------------------------

/// OCCT Geom2dAdaptor_Curve (TKGeomBase/Geom2dAdaptor) — the Adaptor2d_Curve2d
/// instance over a Geom2d_Curve (rcad's `Curve2d`), in the unrestricted
/// constructor form `Geom2dAdaptor_Curve(C)`.
///
/// The adaptor looks through `Curve2d::Trimmed` layers for the type accessors
/// (the OCCT adaptor holds the basis curve) and evaluates through the trimmed
/// parameterization.
pub struct Geom2dCurveAdaptor {
    /// OCCT: handle(Geom2d_Curve) myCurve.
    pub curve: crate::geom::Curve2d,
}

impl Geom2dCurveAdaptor {
    /// OCCT Geom2dAdaptor_Curve(C).
    pub fn new(curve: crate::geom::Curve2d) -> Self {
        Geom2dCurveAdaptor { curve }
    }
}

impl Adaptor2dCurve2d for Geom2dCurveAdaptor {
    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64 {
        self.curve.default_domain()[0]
    }

    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64 {
        self.curve.default_domain()[1]
    }

    /// OCCT Value(U).
    fn value(&self, u: f64) -> DVec2 {
        self.curve.point_at(u)
    }

    /// OCCT D0(U, P).
    fn d0(&self, u: f64) -> DVec2 {
        self.curve.point_at(u)
    }

    /// OCCT D1(U, P, V).
    fn d1(&self, u: f64) -> (DVec2, DVec2) {
        (self.curve.point_at(u), self.curve.derivative_at(u))
    }

    /// OCCT D2(U, P, V1, V2).
    fn d2(&self, u: f64) -> (DVec2, DVec2, DVec2) {
        (
            self.curve.point_at(u),
            self.curve.derivative_at(u),
            self.curve.derivative2_at(u),
        )
    }

    /// OCCT Geom2dAdaptor_Curve::Continuity() — CN for the elementary kinds
    /// (matching the codebase-wide encoding of non-composite rcad curves);
    /// the OCCT BSpline LocalContinuity walk and the Offset basis-continuity
    /// shift are staged.
    fn continuity(&self) -> GeomAbsShape {
        GeomAbsShape::CN
    }

    /// OCCT GetType() — of the basis curve (looking through Trimmed).
    fn get_type(&self) -> CurveType {
        curve2d_type_of(self.curve.inner())
    }

    /// OCCT Line() — valid when GetType() == Line.
    fn line(&self) -> crate::geom::Line2d {
        match self.curve.inner() {
            crate::geom::Curve2d::Line(l) => *l,
            _ => unreachable!("Geom2dCurveAdaptor::Line on a non-line curve"),
        }
    }

    /// OCCT BSpline() — the downcast to Geom2d_BSplineCurve.
    fn bspline(&self) -> Option<crate::geom::BSplineCurve2> {
        match self.curve.inner() {
            crate::geom::Curve2d::BSpline(b) => Some(b.clone()),
            _ => None,
        }
    }

    /// OCCT Bezier() — the downcast to Geom2d_BezierCurve.
    fn bezier(&self) -> Option<crate::geom::BezierCurve2> {
        match self.curve.inner() {
            crate::geom::Curve2d::Bezier(b) => Some(b.clone()),
            _ => None,
        }
    }

    /// OCCT Geom2dAdaptor_Curve::Trim(First, Last, Tol) — restricts the
    /// parameter range (the Geom2d_TrimmedCurve-like restriction).
    fn trim(&self, first_param: f64, last_param: f64, _tol: f64) -> Arc<dyn Adaptor2dCurve2d> {
        Arc::new(Geom2dCurveAdaptor {
            curve: crate::geom::Curve2d::Trimmed(crate::geom::TrimmedCurve2 {
                curve: Box::new(self.curve.clone()),
                t_min: first_param,
                t_max: last_param,
            }),
        })
    }
}

/// GeomAbs_CurveType of a kernel `Curve2d` kind (the Geom2dAdaptor_Curve type
/// dispatch).
fn curve2d_type_of(c: &crate::geom::Curve2d) -> CurveType {
    use crate::geom::Curve2d;
    match c {
        Curve2d::Line(_) => CurveType::Line,
        Curve2d::Circle(_) => CurveType::Circle,
        Curve2d::Ellipse(_) => CurveType::Ellipse,
        Curve2d::Parabola(_) => CurveType::Parabola,
        Curve2d::Hyperbola(_) => CurveType::Hyperbola,
        Curve2d::Bezier(_) => CurveType::Bezier,
        Curve2d::BSpline(_) => CurveType::BSpline,
        _ => CurveType::Other,
    }
}

// ---------------------------------------------------------------------------
// Adaptor3d_HSurfaceTool::IsSurfG1 — GAP carrier
// ---------------------------------------------------------------------------

/// OCCT Adaptor3d_HSurfaceTool::IsSurfG1(S, AlongU, AngleTol)
/// (Adaptor3d_HSurfaceTool.cxx) — checks the G1 continuity across the
/// surface's u/v joints.  GAP (staged): the underlying joint analysis is not
/// translated; every call reports "not G1" which is the conservative answer
/// that drives the Approx_CurveOnSurface tolerance clamping.
pub fn is_surf_g1(_s: &dyn Adaptor3dSurface, _along_u: bool, _angle_tol: f64) -> bool {
    false
}

// ---------------------------------------------------------------------------
// Adaptor3d_CurveOnSurface — the 3D curve-on-surface adaptor
// ---------------------------------------------------------------------------

/// OCCT Adaptor3d_CurveOnSurface (TKG3d) — a 3D curve defined as the image of
/// a 2D curve on a surface.
///
/// The evaluation members are the OCCT definition: the curve point is the
/// surface evaluated at the 2D curve point, and the derivatives chain by the
/// chain rule (Adaptor3d_CurveOnSurface.cxx Value/D1/D2).
///
/// GAP (staged): the interval machinery (NbIntervals/Intervals, which walks
/// the basis-curve and surface discontinuities) and the trimming support are
/// not translated; [`Adaptor3dCurve::nb_intervals`] reports a single interval
/// and [`Adaptor3dCurve::trim`] panics until the body lands.
pub struct CurveOnSurface {
    /// OCCT: const handle(Adaptor2d_Curve2d) my2dCurve.
    pub my2d_curve: Curve2dHandle,
    /// OCCT: handle(Adaptor3d_Surface) mySurface.
    pub my_surface: SurfaceHandle,
}

impl CurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface(C2D, S).
    pub fn new(the2d_curve: Curve2dHandle, the_surface: SurfaceHandle) -> Self {
        CurveOnSurface {
            my2d_curve: the2d_curve,
            my_surface: the_surface,
        }
    }

    /// OCCT GetCurve() — the 2D basis curve.
    pub fn get_curve(&self) -> &Curve2dHandle {
        &self.my2d_curve
    }

    /// OCCT GetSurface() — the surface.
    pub fn get_surface(&self) -> &SurfaceHandle {
        &self.my_surface
    }
}

impl Adaptor3dCurve for CurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface::FirstParameter.
    fn first_parameter(&self) -> f64 {
        self.my2d_curve.first_parameter()
    }

    /// OCCT Adaptor3d_CurveOnSurface::LastParameter.
    fn last_parameter(&self) -> f64 {
        self.my2d_curve.last_parameter()
    }

    /// OCCT Adaptor3d_CurveOnSurface::Value(U): S->Value(C2D->Value(U)).
    fn value(&self, u: f64) -> DVec3 {
        let p2d = self.my2d_curve.value(u);
        self.my_surface.value(p2d.x, p2d.y)
    }

    /// OCCT Adaptor3d_CurveOnSurface::D1(U, P, V):
    ///   S->D1(..., P, DS1u, DS1v); V = DS1u*x' + DS1v*y'.
    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        let p2d = self.my2d_curve.value(u);
        let (p, ds1_u, ds1_v) = self.my_surface.d1(p2d.x, p2d.y);
        let (_p2d, d2d) = self.my2d_curve.d1(u);
        let v = ds1_u * d2d.x + ds1_v * d2d.y;
        (p, v)
    }

    /// OCCT Adaptor3d_CurveOnSurface::D2(U, P, V1, V2):
    ///   V1 = DS1u*x' + DS1v*y',
    ///   V2 = DS2u*x'*x' + DS1u*x'' + 2*DS2uv*x'*y' + DS2v*y'*y' + DS1v*y''.
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        let p2d = self.my2d_curve.value(u);
        let (p, ds1_u, ds1_v, ds2_u, ds2_v, ds2_uv) = self.my_surface.d2(p2d.x, p2d.y);
        let (_p2d, d2d, dd2d) = self.my2d_curve.d2(u);
        let v1 = ds1_u * d2d.x + ds1_v * d2d.y;
        let v2 = ds2_u * d2d.x * d2d.x
            + ds1_u * dd2d.x
            + 2.0 * ds2_uv * d2d.x * d2d.y
            + ds2_v * d2d.y * d2d.y
            + ds1_v * dd2d.y;
        (p, v1, v2)
    }

    /// OCCT Adaptor3d_CurveOnSurface::Continuity — minimum of the basis curve
    /// and the two surface continuities.
    fn continuity(&self) -> GeomAbsShape {
        let cont_c = self.my2d_curve.continuity();
        let cont_su = self.my_surface.u_continuity();
        let cont_c = if (cont_su as u8) < (cont_c as u8) {
            cont_su
        } else {
            cont_c
        };
        let cont_sv = self.my_surface.v_continuity();
        if (cont_sv as u8) < (cont_c as u8) {
            cont_sv
        } else {
            cont_c
        }
    }

    /// GAP (staged): Adaptor3d_CurveOnSurface::NbIntervals — reports a single
    /// interval.
    fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        1
    }

    /// GAP (staged): Adaptor3d_CurveOnSurface::Intervals — reports the whole
    /// domain as a single interval.
    fn intervals(&self, _s: GeomAbsShape) -> Vec<f64> {
        vec![self.first_parameter(), self.last_parameter()]
    }

    /// GAP (staged): Adaptor3d_CurveOnSurface::Trim.
    fn trim(&self, _first: f64, _last: f64, _tol: f64) -> Arc<dyn Adaptor3dCurve> {
        panic!("GAP: Adaptor3d_CurveOnSurface::Trim not translated");
    }
}
