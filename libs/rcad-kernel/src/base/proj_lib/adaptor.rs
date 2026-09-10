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

use crate::core::precision;
use crate::geom::Curve2dEval;
use crate::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};

use crate::math::GeomAbsShape;

use super::CurveType;
use super::proj_lib_projected_curve::Adaptor3dCurveGeom;

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
    /// OCCT Adaptor2d_Curve2d::IsClosed() — the base class raises
    /// Standard_NoSuchObject (Adaptor2d_Curve2d.cxx).
    fn is_closed(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::IsClosed")
    }
    /// OCCT Adaptor2d_Curve2d::IsPeriodic() — consumed by
    /// Adaptor3d_CurveOnSurface::IsPeriodic (Adaptor3d_CurveOnSurface.cxx
    /// L1151-1160).
    fn is_periodic(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::IsPeriodic")
    }
    /// OCCT Adaptor2d_Curve2d::Period() — the base class raises
    /// Standard_NoSuchObject (Adaptor2d_Curve2d.cxx).
    fn period(&self) -> f64 {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Period")
    }
    /// OCCT Adaptor2d_Curve2d::Resolution(R3d) — the base class raises
    /// Standard_NoSuchObject (Adaptor2d_Curve2d.hxx L124).
    fn resolution(&self, _r3d: f64) -> f64 {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Resolution")
    }
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
    /// OCCT Adaptor3d_Curve::ShallowCopy() (Adaptor3d_Curve.cxx L38-42) — the
    /// base class throws Standard_NotImplemented.
    fn shallow_copy(&self) -> Arc<dyn Adaptor3dCurve> {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::ShallowCopy")
    }
    /// OCCT Adaptor3d_Curve::IsClosed() (Adaptor3d_Curve.cxx L89-93).
    fn is_closed(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::IsClosed")
    }
    /// OCCT Adaptor3d_Curve::IsPeriodic() (Adaptor3d_Curve.cxx L96-100).
    fn is_periodic(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::IsPeriodic")
    }
    /// OCCT Adaptor3d_Curve::Period() (Adaptor3d_Curve.cxx L103-108).
    fn period(&self) -> f64 {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::Period")
    }
    /// OCCT Adaptor3d_Curve::Resolution(R3d) (Adaptor3d_Curve.cxx L111-115).
    fn resolution(&self, _r3d: f64) -> f64 {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::Resolution")
    }
    /// OCCT Adaptor3d_Curve::D3(U, P, V1, V2, V3) — the EvalD3 composition
    /// (Adaptor3d_Curve.hxx L116-124); the base EvalD3 throws
    /// Standard_NotImplemented (Adaptor3d_Curve.cxx L225-232).
    fn d3(&self, _u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::EvalD3")
    }
    /// OCCT Adaptor3d_Curve::DN(U, N) — the EvalDN composition
    /// (Adaptor3d_Curve.hxx L126-130).  The rcad default keeps the
    /// established order-dispatch encoding: N=1 delegates to D1, N=2 to D2
    /// and the remaining orders raise the OCCT Standard_NotImplemented of
    /// the base EvalDN (Adaptor3d_Curve.cxx L235-244).
    fn dn(&self, u: f64, n: i32) -> DVec3 {
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            _ => panic!("Standard_NotImplemented: Adaptor3d_Curve::EvalDN"),
        }
    }
    /// OCCT Adaptor3d_Curve::Degree() (Adaptor3d_Curve.cxx L160-164).
    fn degree(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::Degree")
    }
    /// OCCT Adaptor3d_Curve::IsRational() (Adaptor3d_Curve.cxx L167-171).
    fn is_rational(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::IsRational")
    }
    /// OCCT Adaptor3d_Curve::NbPoles() (Adaptor3d_Curve.cxx L174-178).
    fn nb_poles(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::NbPoles")
    }
    /// OCCT Adaptor3d_Curve::NbKnots() (Adaptor3d_Curve.cxx L181-185).
    fn nb_knots(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::NbKnots")
    }
    /// OCCT Adaptor3d_Curve::Bezier() (Adaptor3d_Curve.cxx L188-192) — the
    /// Geom_BezierCurve payload.
    fn bezier(&self) -> crate::geom::BezierCurve3 {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::Bezier")
    }
    /// OCCT Adaptor3d_Curve::BSpline() (Adaptor3d_Curve.cxx L195-199) — the
    /// Geom_BSplineCurve payload.
    fn bspline(&self) -> crate::geom::BSplineCurve3 {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::BSpline")
    }
    /// OCCT Adaptor3d_Curve::OffsetCurve() (Adaptor3d_Curve.cxx L202-206) —
    /// the Geom_OffsetCurve payload.
    fn offset_curve(&self) -> crate::geom::OffsetCurve3 {
        panic!("Standard_NotImplemented: Adaptor3d_Curve::OffsetCurve")
    }
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
    /// The base class throws Standard_NotImplemented
    /// (Adaptor3d_Surface.cxx L123-131).
    fn u_trim(&self, _u1: f64, _u2: f64, _eps: f64) -> Arc<dyn Adaptor3dSurface> {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::UTrim")
    }
    /// OCCT VTrim(V1, V2, Eps) — the surface restricted in V.
    /// The base class throws Standard_NotImplemented
    /// (Adaptor3d_Surface.cxx L134-142).
    fn v_trim(&self, _v1: f64, _v2: f64, _eps: f64) -> Arc<dyn Adaptor3dSurface> {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::VTrim")
    }
    /// OCCT Adaptor3d_Surface::IsUClosed() (Adaptor3d_Surface.cxx L145-149).
    fn is_u_closed(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::IsUClosed")
    }
    /// OCCT Adaptor3d_Surface::IsVClosed() (Adaptor3d_Surface.cxx L152-156).
    fn is_v_closed(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::IsVClosed")
    }
    /// OCCT Adaptor3d_Surface::UDegree() (Adaptor3d_Surface.cxx).
    fn u_degree(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::UDegree")
    }
    /// OCCT Adaptor3d_Surface::NbUPoles() (Adaptor3d_Surface.cxx).
    fn nb_u_poles(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::NbUPoles")
    }
    /// OCCT Adaptor3d_Surface::VDegree() (Adaptor3d_Surface.cxx).
    fn v_degree(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::VDegree")
    }
    /// OCCT Adaptor3d_Surface::NbVPoles() (Adaptor3d_Surface.cxx).
    fn nb_v_poles(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::NbVPoles")
    }
    /// OCCT Adaptor3d_Surface::NbUKnots() (Adaptor3d_Surface.cxx).
    fn nb_u_knots(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::NbUKnots")
    }
    /// OCCT Adaptor3d_Surface::NbVKnots() (Adaptor3d_Surface.cxx).
    fn nb_v_knots(&self) -> usize {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::NbVKnots")
    }
    /// OCCT Adaptor3d_Surface::IsURational() (Adaptor3d_Surface.cxx).
    fn is_u_rational(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::IsURational")
    }
    /// OCCT Adaptor3d_Surface::IsVRational() (Adaptor3d_Surface.cxx).
    fn is_v_rational(&self) -> bool {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::IsVRational")
    }
    /// OCCT Adaptor3d_Surface::Bezier() — the Geom_BezierSurface payload.
    fn bezier(&self) -> crate::geom::BezierSurface {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::Bezier")
    }
    /// OCCT Adaptor3d_Surface::BSpline() — the Geom_BSplineSurface payload.
    fn bspline(&self) -> crate::geom::BSplineSurface {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::BSpline")
    }
    /// OCCT Adaptor3d_Surface::Direction() — the extrusion direction
    /// (gp_Dir payload).
    fn direction(&self) -> DVec3 {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::Direction")
    }
    /// OCCT Adaptor3d_Surface::BasisCurve() — the basis curve adaptor of an
    /// extrusion / revolution surface.
    fn basis_curve(&self) -> Arc<dyn Adaptor3dCurve> {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::BasisCurve")
    }
    /// OCCT Adaptor3d_Surface::BasisSurface() — the basis surface adaptor of
    /// an offset surface.
    fn basis_surface(&self) -> Arc<dyn Adaptor3dSurface> {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::BasisSurface")
    }
    /// OCCT Adaptor3d_Surface::OffsetValue().
    fn offset_value(&self) -> f64 {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::OffsetValue")
    }
    /// OCCT Adaptor3d_Surface::DN(U, V, Nu, Nv) — the EvalDN composition;
    /// the base EvalDN throws Standard_NotImplemented
    /// (Adaptor3d_Surface.cxx).
    fn dn(&self, _u: f64, _v: f64, _nu: i32, _nv: i32) -> DVec3 {
        panic!("Standard_NotImplemented: Adaptor3d_Surface::EvalDN")
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
    /// OCCT: Standard_Real myFirst — the restricted first parameter
    /// (Geom2dAdaptor_Curve::Load(UFirst, ULast) form).
    pub first: f64,
    /// OCCT: Standard_Real myLast.
    pub last: f64,
}

impl Geom2dCurveAdaptor {
    /// OCCT Geom2dAdaptor_Curve(C) -> Load(C): the curve natural domain
    /// (Geom2dAdaptor_Curve.hxx Load(theCurve)).
    pub fn new(curve: crate::geom::Curve2d) -> Self {
        let [first, last] = curve.default_domain();
        Geom2dCurveAdaptor {
            curve,
            first,
            last,
        }
    }

    /// OCCT Geom2dAdaptor_Curve(C, First, Last) -> Load(C, First, Last):
    /// the restricted parameter window (the form used by
    /// BRepAdaptor_Curve.cxx L202-203: HC->Load(PC, pf, pl)).
    pub fn with_range(curve: crate::geom::Curve2d, first: f64, last: f64) -> Self {
        Geom2dCurveAdaptor {
            curve,
            first,
            last,
        }
    }
}

impl Adaptor2dCurve2d for Geom2dCurveAdaptor {
    /// OCCT FirstParameter() — the restricted first parameter.
    fn first_parameter(&self) -> f64 {
        self.first
    }

    /// OCCT LastParameter() — the restricted last parameter.
    fn last_parameter(&self) -> f64 {
        self.last
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
            first: first_param,
            last: last_param,
        })
    }

    /// OCCT Geom2dAdaptor_Curve::IsClosed (Geom2dAdaptor_Curve.cxx
    /// L588-602): the endpoint-distance test on the restricted domain.  The
    /// OCCT Precision::IsPositiveInfinite / IsNegativeInfinite guards ride
    /// the IEEE infinity encoding (see the Precision::Infinite coordination
    /// batch note in docs).
    fn is_closed(&self) -> bool {
        let last = self.last;
        let first = self.first;
        if last != f64::INFINITY && first != f64::NEG_INFINITY {
            let pd = self.value(first);
            let pf = self.value(last);
            return (pf - pd).length() <= precision::CONFUSION;
        }
        false
    }

    /// OCCT Geom2dAdaptor_Curve::IsPeriodic (Geom2dAdaptor_Curve.cxx
    /// L604-607): the basis curve periodicity.
    fn is_periodic(&self) -> bool {
        self.curve.is_periodic()
    }

    /// OCCT Geom2dAdaptor_Curve::Period (Geom2dAdaptor_Curve.cxx L611-614):
    /// myCurve->LastParameter() - myCurve->FirstParameter().
    fn period(&self) -> f64 {
        self.curve.default_domain()[1] - self.curve.default_domain()[0]
    }

    /// OCCT Geom2dAdaptor_Curve::Resolution(Ruv)
    /// (Geom2dAdaptor_Curve.cxx L1186-1222) — the per-type parametric
    /// resolution.
    fn resolution(&self, ruv: f64) -> f64 {
        use crate::geom::Curve2d;
        match self.curve.inner() {
            Curve2d::Line(_) => ruv,
            Curve2d::Circle(c) => {
                let r = c.radius;
                if r > ruv / 2.0 {
                    2.0 * (ruv / (2.0 * r)).asin()
                } else {
                    2.0 * std::f64::consts::PI
                }
            }
            Curve2d::Ellipse(e) => ruv / e.major_radius,
            // OCCT L1206-1214: Geom2d_BezierCurve / Geom2d_BSplineCurve
            // Resolution — the BSplCLib pole-magnitude engine; the rcad 2D
            // poles ride the 3D engine through the z=0 embedding (pole
            // magnitudes are invariant).
            Curve2d::BSpline(b) => {
                let poles3: Vec<glam::DVec3> = b
                    .control_points
                    .iter()
                    .map(|p| glam::DVec3::new(p.x, p.y, 0.0))
                    .collect();
                crate::math::bspl::bspl_curve_resolution(
                    &poles3,
                    Some(&b.weights[..]),
                    &b.knots,
                    b.degree,
                    ruv,
                )
            }
            Curve2d::Bezier(b) => {
                let poles3: Vec<glam::DVec3> = b
                    .control_points
                    .iter()
                    .map(|p| glam::DVec3::new(p.x, p.y, 0.0))
                    .collect();
                let degree = b.control_points.len() - 1;
                let mut flat_knots = vec![0.0; degree + 1];
                flat_knots.extend(std::iter::repeat_n(1.0, degree + 1));
                crate::math::bspl::bspl_curve_resolution(
                    &poles3,
                    Some(&b.weights[..]),
                    &flat_knots,
                    degree,
                    ruv,
                )
            }
            // OCCT default arm L1216-1218: Precision::Parametric(Ruv).
            _ => precision::parametric_default(ruv),
        }
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
    /// OCCT: GeomAbs_CurveType myType — set by the Load/EvalKPart pair; the
    /// constructor default is GeomAbs_OtherCurve
    /// (Adaptor3d_CurveOnSurface.cxx L873-876).  GAP (staged): the EvalKPart
    /// refinement (Adaptor3d_CurveOnSurface.cxx L1553-1830) is not
    /// translated, so the type keeps the constructor value.
    pub my_type: CurveType,
}

impl CurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface(C2D, S) (Adaptor3d_CurveOnSurface.cxx
    /// L881-889): myType = GeomAbs_OtherCurve, myIntCont = GeomAbs_CN,
    /// Load(S), Load(C).
    pub fn new(the2d_curve: Curve2dHandle, the_surface: SurfaceHandle) -> Self {
        CurveOnSurface {
            my2d_curve: the2d_curve,
            my_surface: the_surface,
            my_type: CurveType::Other,
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

    /// GAP (staged): Adaptor3d_CurveOnSurface::Trim — the OCCT body
    /// (Adaptor3d_CurveOnSurface.cxx L1133-1141) rebuilds the adaptor with
    /// Load(mySurface) + Load(myCurve->Trim(First, Last, Tol)).
    fn trim(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurve> {
        let hcs = CurveOnSurface::new(self.my2d_curve.trim(first, last, tol), self.my_surface.clone());
        Arc::new(hcs)
    }

    /// OCCT Adaptor3d_CurveOnSurface::IsClosed (Adaptor3d_CurveOnSurface.cxx
    /// L1145-1148): myCurve->IsClosed().
    fn is_closed(&self) -> bool {
        self.my2d_curve.is_closed()
    }

    /// OCCT Adaptor3d_CurveOnSurface::IsPeriodic
    /// (Adaptor3d_CurveOnSurface.cxx L1151-1160): true for the Circle /
    /// Ellipse myType, else myCurve->IsPeriodic().
    fn is_periodic(&self) -> bool {
        if self.my_type == CurveType::Circle || self.my_type == CurveType::Ellipse {
            return true;
        }
        self.my2d_curve.is_periodic()
    }

    /// OCCT Adaptor3d_CurveOnSurface::Period
    /// (Adaptor3d_CurveOnSurface.cxx L1163-1170): 2*pi for the Circle /
    /// Ellipse myType, else myCurve->Period().
    fn period(&self) -> f64 {
        if self.my_type == CurveType::Circle || self.my_type == CurveType::Ellipse {
            return std::f64::consts::TAU;
        }
        self.my2d_curve.period()
    }

    /// OCCT Adaptor3d_CurveOnSurface::Resolution
    /// (Adaptor3d_CurveOnSurface.cxx L1364-1371):
    /// myCurve->Resolution(min(UResolution(R3d), VResolution(R3d))).
    fn resolution(&self, r3d: f64) -> f64 {
        let ru = self.my_surface.u_resolution(r3d);
        let rv = self.my_surface.v_resolution(r3d);
        self.my2d_curve.resolution(ru.min(rv))
    }
}

impl Adaptor3dCurveGeom for CurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface::GetType
    /// (Adaptor3d_CurveOnSurface.cxx L1374-1377): returns myType.
    fn get_type(&self) -> CurveType {
        self.my_type
    }

    /// OCCT Adaptor3d_CurveOnSurface::Line (L1381-1387) — the myLin payload
    /// cached by EvalKPart; the Standard_NoSuchObject raise fires when
    /// myType != GeomAbs_Line.
    fn line(&self) -> Line3 {
        panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Line(): curve is not a line")
    }

    /// OCCT Adaptor3d_CurveOnSurface::Circle (L1390-1396) — the myCirc
    /// payload cached by EvalKPart; raises when myType != GeomAbs_Circle.
    fn circle(&self) -> Circle3 {
        panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Line(): curve is not a circle")
    }

    /// OCCT Adaptor3d_CurveOnSurface::Ellipse (L1399-1403) —
    /// to3d(mySurface->Plane(), myCurve->Ellipse()); unreachable with the
    /// staged EvalKPart (see the struct docs).
    fn ellipse(&self) -> Ellipse3 {
        panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Ellipse")
    }

    /// OCCT Adaptor3d_CurveOnSurface::Hyperbola (L1406-1410) — see Ellipse.
    fn hyperbola(&self) -> Hyperbola3 {
        panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Hyperbola")
    }

    /// OCCT Adaptor3d_CurveOnSurface::Parabola (L1413-1417) — see Ellipse.
    fn parabola(&self) -> Parabola3 {
        panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Parabola")
    }

    /// OCCT Adaptor3d_CurveOnSurface::Trim — see the Adaptor3dCurve trim.
    fn trim_geom(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurveGeom> {
        let hcs = CurveOnSurface::new(self.my2d_curve.trim(first, last, tol), self.my_surface.clone());
        Arc::new(hcs)
    }
}
