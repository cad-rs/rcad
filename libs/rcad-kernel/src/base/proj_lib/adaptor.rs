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
    /// Standard_NotImplemented (Adaptor2d_Curve2d.cxx L136-140); the
    /// Geom2dAdaptor_Curve override is the per-type dispatch
    /// (Geom2dAdaptor_Curve.cxx L1186-1226).
    fn resolution(&self, _r3d: f64) -> f64 {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Resolution")
    }
    /// OCCT Adaptor2d_Curve2d::NbIntervals(S) — the base class raises
    /// Standard_NotImplemented (Adaptor2d_Curve2d.cxx L64-68).
    fn nb_intervals(&self, _s: GeomAbsShape) -> usize {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::NbIntervals")
    }
    /// OCCT Adaptor2d_Curve2d::Intervals(T, S) — the base class raises
    /// Standard_NotImplemented (Adaptor2d_Curve2d.cxx L72-77).
    fn intervals(&self, _s: GeomAbsShape) -> Vec<f64> {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Intervals")
    }
    /// OCCT Adaptor2d_Curve2d::Circle() — the gp_Circ2d payload; the base
    /// class raises Standard_NotImplemented (Adaptor2d_Curve2d.cxx L188-192).
    fn circle(&self) -> crate::geom::Circle2d {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Circle")
    }
    /// OCCT Adaptor2d_Curve2d::Ellipse() — the gp_Elips2d payload; the base
    /// class raises Standard_NotImplemented (Adaptor2d_Curve2d.cxx L195-199).
    fn ellipse(&self) -> crate::geom::Ellipse2d {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Ellipse")
    }
    /// OCCT Adaptor2d_Curve2d::Hyperbola() — the gp_Hypr2d payload; the base
    /// class raises Standard_NotImplemented (Adaptor2d_Curve2d.cxx L202-206).
    fn hyperbola(&self) -> crate::geom::Hyperbola2d {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Hyperbola")
    }
    /// OCCT Adaptor2d_Curve2d::Parabola() — the gp_Parab2d payload; the base
    /// class raises Standard_NotImplemented (Adaptor2d_Curve2d.cxx L209-213).
    fn parabola(&self) -> crate::geom::Parabola2d {
        panic!("Standard_NotImplemented: Adaptor2d_Curve2d::Parabola")
    }
    /// The kernel curve when the adaptor wraps one — the routing bridge for
    /// the `occ::down_cast<Geom2dAdaptor_Curve>` of the GeomLib / Approx
    /// consumers (the same shape as
    /// [`Adaptor3dSurface::kernel_surface`]): Some when the adaptor is a
    /// `Geom2dAdaptor_Curve` over a kernel `Curve2d`, None when the OCCT
    /// down-cast would leave a null handle.
    fn kernel_curve2d(&self) -> Option<&crate::geom::Curve2d> {
        None
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

    /// OCCT Geom2dAdaptor_Curve::Continuity()
    /// (Geom2dAdaptor_Curve.cxx L368-399): LocalContinuity for the BSpline
    /// kind, the shifted basis continuity for the offset kind, CN for the
    /// elementary kinds.
    fn continuity(&self) -> GeomAbsShape {
        use crate::geom::Curve2d;
        match self.curve.inner() {
            // OCCT L370-373: LocalContinuity(myFirst, myLast).
            Curve2d::BSpline(_) => self.local_continuity(self.first, self.last),
            Curve2d::Offset(_) => {
                // OCCT L374-392: the basis continuity shifted down by one
                // degree (Geom2d_OffsetCurve::GetBasisCurveContinuity); the
                // rcad basis carries the Geom2d elementary default (CN), so
                // the shifted value stays CN — the OCCT shift arms kept in
                // shape below.
                let base = GeomAbsShape::CN;
                match base {
                    GeomAbsShape::CN => GeomAbsShape::CN,
                    GeomAbsShape::C3 => GeomAbsShape::C2,
                    GeomAbsShape::C2 => GeomAbsShape::C1,
                    GeomAbsShape::C1 => GeomAbsShape::C0,
                    _ => panic!("Standard_NoSuchObject: Geom2dAdaptor_Curve::Continuity"),
                }
            }
            // OCCT L394-397: the OtherCurve kind raises.
            Curve2d::CircleInvolute(_)
            | Curve2d::ArchimedeanSpiral(_)
            | Curve2d::LogarithmicSpiral(_)
            | Curve2d::SineWave(_) => {
                panic!("Standard_NoSuchObject: Geom2dAdaptor_Curve::Continuity")
            }
            // OCCT L399: the elementary kinds answer CN.
            _ => GeomAbsShape::CN,
        }
    }

    /// OCCT GetType() — of the basis curve (looking through Trimmed).
    fn get_type(&self) -> CurveType {
        curve2d_type_of(self.curve.inner())
    }

    /// The kernel `Geom2d_Curve` wrapped by this `Geom2dAdaptor_Curve` — the
    /// `Curve()` accessor of the GeomLib / Approx down-cast sites.
    fn kernel_curve2d(&self) -> Option<&crate::geom::Curve2d> {
        Some(&self.curve)
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

    /// OCCT Geom2dAdaptor_Curve::Trim(First, Last, Tol) (Geom2dAdaptor_Curve.cxx
    /// L577-584): `return new Geom2dAdaptor_Curve(myCurve, First, Last)` — a new
    /// adaptor over the SAME basis curve with the restricted range.  The range
    /// lives in myFirst/myLast; it is NOT wrapped into a Geom2d_TrimmedCurve
    /// (which clamps on evaluation, so successive trims would nest clamps and
    /// saturate the parameter — the OrderedApprox path trims once per interval).
    fn trim(&self, first_param: f64, last_param: f64, _tol: f64) -> Arc<dyn Adaptor2dCurve2d> {
        Arc::new(Geom2dCurveAdaptor {
            curve: self.curve.clone(),
            first: first_param,
            last: last_param,
        })
    }

    /// OCCT Geom2dAdaptor_Curve::IsClosed (Geom2dAdaptor_Curve.cxx
    /// L588-600): the endpoint-distance test on the restricted domain.
    /// OCCT L590: !Precision::IsPositiveInfinite(myLast) &&
    /// !Precision::IsNegativeInfinite(myFirst) (Precision.hxx L357-367,
    /// threshold 0.5 * Precision::Infinite()).
    fn is_closed(&self) -> bool {
        let last = self.last;
        let first = self.first;
        if !precision::is_positive_infinite_value(last)
            && !precision::is_negative_infinite_value(first)
        {
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

    /// OCCT Geom2dAdaptor_Curve::Circle() (Geom2dAdaptor_Curve.cxx
    /// L1237-1243) — the gp_Circ2d payload; raises on a type mismatch.
    fn circle(&self) -> crate::geom::Circle2d {
        match self.curve.inner() {
            crate::geom::Curve2d::Circle(c) => *c,
            _ => panic!("Standard_NoSuchObject: Geom2dAdaptor_Curve::Circle() - curve is not a Circle"),
        }
    }

    /// OCCT Geom2dAdaptor_Curve::Ellipse() (Geom2dAdaptor_Curve.cxx
    /// L1246-1252) — the gp_Elips2d payload.
    fn ellipse(&self) -> crate::geom::Ellipse2d {
        match self.curve.inner() {
            crate::geom::Curve2d::Ellipse(e) => *e,
            _ => panic!("Standard_NoSuchObject: Geom2dAdaptor_Curve::Ellipse() - curve is not an Ellipse"),
        }
    }

    /// OCCT Geom2dAdaptor_Curve::Hyperbola() (Geom2dAdaptor_Curve.cxx
    /// L1255-1261) — the gp_Hypr2d payload.
    fn hyperbola(&self) -> crate::geom::Hyperbola2d {
        match self.curve.inner() {
            crate::geom::Curve2d::Hyperbola(h) => *h,
            _ => panic!("Standard_NoSuchObject: Geom2dAdaptor_Curve::Hyperbola() - curve is not an Hyperbola"),
        }
    }

    /// OCCT Geom2dAdaptor_Curve::Parabola() (Geom2dAdaptor_Curve.cxx
    /// L1264-1270) — the gp_Parab2d payload.
    fn parabola(&self) -> crate::geom::Parabola2d {
        match self.curve.inner() {
            crate::geom::Curve2d::Parabola(p) => *p,
            _ => panic!("Standard_NoSuchObject: Geom2dAdaptor_Curve::Parabola() - curve is not a Parabola"),
        }
    }

    /// OCCT Geom2dAdaptor_Curve::NbIntervals(S)
    /// (Geom2dAdaptor_Curve.cxx L409-486).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        use crate::geom::Curve2d;
        let inner = self.curve.inner();
        if let Curve2d::BSpline(a_bspline) = inner {
            // OCCT L414-417.
            if (!bspline_is_periodic(a_bspline) && s <= self.continuity()) || s == GeomAbsShape::C0 {
                return 1;
            }

            // OCCT L419-438.
            let a_degree = a_bspline.degree;
            let a_cont = match s {
                GeomAbsShape::C1 => 1,
                GeomAbsShape::C2 => 2,
                GeomAbsShape::C3 => 3,
                GeomAbsShape::CN => a_degree as i32,
                // OCCT L437: the G1/G2 arms raise Standard_DomainError
                // (the rcad GeomAbsShape carries no G kinds).
                _ => panic!("Standard_DomainError: Geom2dAdaptor_Curve::NbIntervals()"),
            };

            // OCCT L440.
            let an_eps = self
                .resolution(precision::CONFUSION)
                .min(precision::p_confusion());

            // OCCT L442-450: BSplCLib::Intervals with the null output array
            // returns the count; the rcad helper returns count + 1 entries.
            let (tk, tm) = knots_mults_of(&a_bspline.knots);
            let out = crate::math::bspl_lib::intervals(
                &tk,
                &tm,
                a_degree,
                bspline_is_periodic(a_bspline),
                a_cont,
                self.first,
                self.last,
                an_eps,
            );
            out.len() - 1
        } else if let Curve2d::Offset(an_offset) = inner {
            // OCCT L453-480.
            let base_s = match s {
                // OCCT L459-462: the G1/G2 arms raise Standard_DomainError.
                GeomAbsShape::C0 => GeomAbsShape::C1,
                GeomAbsShape::C1 => GeomAbsShape::C2,
                GeomAbsShape::C2 => GeomAbsShape::C3,
                _ => GeomAbsShape::CN,
            };
            let an_adaptor = Geom2dCurveAdaptor::new((*an_offset.basis).clone());
            an_adaptor.nb_intervals(base_s)
        } else {
            // OCCT L482-485.
            1
        }
    }

    /// OCCT Geom2dAdaptor_Curve::Intervals(T, S)
    /// (Geom2dAdaptor_Curve.cxx L490-570).
    fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        use crate::geom::Curve2d;
        let inner = self.curve.inner();
        if let Curve2d::BSpline(a_bspline) = inner {
            // OCCT L495-500.
            if (!bspline_is_periodic(a_bspline) && s <= self.continuity()) || s == GeomAbsShape::C0 {
                return vec![self.first, self.last];
            }

            // OCCT L502-521.
            let a_degree = a_bspline.degree;
            let a_cont = match s {
                GeomAbsShape::C1 => 1,
                GeomAbsShape::C2 => 2,
                GeomAbsShape::C3 => 3,
                GeomAbsShape::CN => a_degree as i32,
                _ => panic!("Standard_DomainError: Geom2dAdaptor_Curve::Intervals()"),
            };

            // OCCT L523.
            let an_eps = self
                .resolution(precision::CONFUSION)
                .min(precision::p_confusion());

            // OCCT L525-533.
            let (tk, tm) = knots_mults_of(&a_bspline.knots);
            crate::math::bspl_lib::intervals(
                &tk,
                &tm,
                a_degree,
                bspline_is_periodic(a_bspline),
                a_cont,
                self.first,
                self.last,
                an_eps,
            )
        } else if let Curve2d::Offset(an_offset) = inner {
            // OCCT L534-563.
            let base_s = match s {
                GeomAbsShape::C0 => GeomAbsShape::C1,
                GeomAbsShape::C1 => GeomAbsShape::C2,
                GeomAbsShape::C2 => GeomAbsShape::C3,
                _ => GeomAbsShape::CN,
            };
            let an_adaptor = Geom2dCurveAdaptor::new((*an_offset.basis).clone());
            // OCCT L578-580: anAdaptor.Intervals(T, BaseS) with the range
            // boundaries rewritten from the restricted window.
            let my_nb_intervals = an_adaptor.nb_intervals(base_s);
            let mut t = an_adaptor.intervals(base_s);
            t[0] = self.first;
            t[my_nb_intervals] = self.last;
            t
        } else {
            // OCCT L565-569.
            vec![self.first, self.last]
        }
    }
}

/// The (Knots, Multiplicities) pair of an rcad flat knot vector
/// (run-length compression; OCCT stores the pair directly).
fn knots_mults_of(flat: &[f64]) -> (Vec<f64>, Vec<i32>) {
    let mut knots: Vec<f64> = Vec::new();
    let mut mults: Vec<i32> = Vec::new();
    for &k in flat {
        match knots.last() {
            Some(&last) if last == k => {
                let n = mults.len();
                mults[n - 1] += 1;
            }
            _ => {
                knots.push(k);
                mults.push(1);
            }
        }
    }
    (knots, mults)
}

/// OCCT Geom2d_BSplineCurve::IsPeriodic() in the rcad encoding — the rcad 2D
/// BSpline payload carries a clamped flat knot vector without a periodicity
/// flag, so the adaptor answers false (architecture mapping).
fn bspline_is_periodic(_b: &crate::geom::BSplineCurve2) -> bool {
    false
}

impl Geom2dCurveAdaptor {
    /// OCCT Geom2dAdaptor_Curve::LocalContinuity(U1, U2)
    /// (Geom2dAdaptor_Curve.cxx L145-227) — the BSpline continuity between
    /// two parameters: C(degree - max knot multiplicity in (U1, U2)).
    fn local_continuity(&self, u1: f64, u2: f64) -> GeomAbsShape {
        let crate::geom::Curve2d::BSpline(a_bspline) = self.curve.inner() else {
            // OCCT L147: Standard_NoSuchObject_Raise_if.
            panic!("Standard_NoSuchObject: Geom2dAdaptor_Curve::LocalContinuity");
        };
        let (tk, tm) = knots_mults_of(&a_bspline.knots);
        let nb = tk.len() as i32; // OCCT L149: aBSpline->NbKnots().
        let mut index1 = 0i32;
        let mut index2 = 0i32;
        let mut new_first = 0.0f64;
        let mut new_last = 0.0f64;
        // OCCT L155-172: BSplCLib::LocateParameter for U1 and U2.
        crate::math::bspl_lib::locate_parameter_knots_mults(
            a_bspline.degree,
            &tk,
            &tm,
            u1,
            bspline_is_periodic(a_bspline),
            1,
            nb,
            &mut index1,
            &mut new_first,
        );
        crate::math::bspl_lib::locate_parameter_knots_mults(
            a_bspline.degree,
            &tk,
            &tm,
            u2,
            bspline_is_periodic(a_bspline),
            1,
            nb,
            &mut index2,
            &mut new_last,
        );
        let a_periodic = bspline_is_periodic(a_bspline);
        // OCCT L173-179.
        if (new_first - crate::math::bspl_lib::at(&tk, index1 + 1)).abs() < precision::p_confusion()
        {
            if index1 < nb {
                index1 += 1;
            }
        }
        // OCCT L180-183.
        if (new_last - crate::math::bspl_lib::at(&tk, index2)).abs() < precision::p_confusion() {
            index2 -= 1;
        }
        let mut mult_max;
        // OCCT L185-189: beware of periodic curves.
        if a_periodic && index1 == nb {
            index1 = 1;
        }

        // OCCT L191-206.
        if (index2 - index1 <= 0) && !a_periodic {
            mult_max = 100; // CN between 2 consecutive nodes
        } else {
            mult_max = crate::math::bspl_lib::ati(&tm, index1 + 1);
            for i in index1 + 1..=index2 {
                let m = crate::math::bspl_lib::ati(&tm, i);
                if m > mult_max {
                    mult_max = m;
                }
            }
            mult_max = a_bspline.degree as i32 - mult_max;
        }
        // OCCT L207-226.
        if mult_max <= 0 {
            GeomAbsShape::C0
        } else if mult_max == 1 {
            GeomAbsShape::C1
        } else if mult_max == 2 {
            GeomAbsShape::C2
        } else if mult_max == 3 {
            GeomAbsShape::C3
        } else {
            GeomAbsShape::CN
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
// Adaptor3d_Surface geometry payloads — the kernel_surface() routing
// ---------------------------------------------------------------------------

// The OCCT virtual payload accessors (Adaptor3d_Surface.cxx L210-243: Plane /
// Cylinder / Cone / Sphere / Torus, each raising Standard_NoSuchObject on the
// base class) cannot ride the [`Adaptor3dSurface`] trait as methods: the
// companion `Adaptor3dSurfaceGeom` subtrait already declares the same names
// and the trait-object receivers would turn ambiguous.  The free functions
// below keep the same encoding as the other HSurfaceTool one-liners (direct
// routing through the `kernel_surface()` bridge); a missing or mismatched
// payload preserves the OCCT raise.

/// OCCT Adaptor3d_Surface::Plane() — the gp_Pln payload
/// (Adaptor3d_Surface.cxx L210-215).
pub fn surface_plane(s: &dyn Adaptor3dSurface) -> crate::geom::Plane {
    match s.kernel_surface() {
        Some(crate::geom::Surface3::Plane(g)) => *g,
        _ => panic!("Standard_NoSuchObject: Adaptor3d_Surface::Plane"),
    }
}

/// OCCT Adaptor3d_Surface::Cylinder() — the gp_Cylinder payload
/// (Adaptor3d_Surface.cxx L217-222).
pub fn surface_cylinder(s: &dyn Adaptor3dSurface) -> crate::geom::CylindricalSurface {
    match s.kernel_surface() {
        Some(crate::geom::Surface3::Cylinder(g)) => *g,
        _ => panic!("Standard_NoSuchObject: Adaptor3d_Surface::Cylinder"),
    }
}

/// OCCT Adaptor3d_Surface::Cone() — the gp_Cone payload
/// (Adaptor3d_Surface.cxx L224-229).
pub fn surface_cone(s: &dyn Adaptor3dSurface) -> crate::geom::ConicalSurface {
    match s.kernel_surface() {
        Some(crate::geom::Surface3::Cone(g)) => *g,
        _ => panic!("Standard_NoSuchObject: Adaptor3d_Surface::Cone"),
    }
}

/// OCCT Adaptor3d_Surface::Sphere() — the gp_Sphere payload
/// (Adaptor3d_Surface.cxx L231-236).
pub fn surface_sphere(s: &dyn Adaptor3dSurface) -> crate::geom::SphericalSurface {
    match s.kernel_surface() {
        Some(crate::geom::Surface3::Sphere(g)) => *g,
        _ => panic!("Standard_NoSuchObject: Adaptor3d_Surface::Sphere"),
    }
}

/// OCCT Adaptor3d_Surface::Torus() — the gp_Torus payload
/// (Adaptor3d_Surface.cxx L238-243).
pub fn surface_torus(s: &dyn Adaptor3dSurface) -> crate::geom::ToroidalSurface {
    match s.kernel_surface() {
        Some(crate::geom::Surface3::Torus(g)) => *g,
        _ => panic!("Standard_NoSuchObject: Adaptor3d_Surface::Torus"),
    }
}

// ---------------------------------------------------------------------------
// Adaptor3d_CurveOnSurface — the 3D curve-on-surface adaptor
// ---------------------------------------------------------------------------

// The file-static `to3d` set of Adaptor3d_CurveOnSurface.cxx L57-98: the
// plane lift of the 2D gp payloads (Pnt2d / Vec2d / Ax22d / Circ2d / Elips2d
// / Hypr2d / Parab2d).

/// OCCT Adaptor3d_CurveOnSurface.cxx L57-60: `to3d(Pl, P)` —
/// ElSLib::Value(P.X(), P.Y(), Pl).
fn to3d_pnt(pl: &crate::geom::Plane, p: DVec2) -> DVec3 {
    crate::math::el::elslib_plane_value(p.x, p.y, pl.origin, pl.u_dir, pl.v_dir)
}

/// OCCT Adaptor3d_CurveOnSurface.cxx L62-70: `to3d(Pl, V)` —
/// V = XDirection*V.X + YDirection*V.Y.
fn to3d_vec(pl: &crate::geom::Plane, v: DVec2) -> DVec3 {
    let vx = pl.u_dir * v.x;
    let vy = pl.v_dir * v.y;
    vx + vy
}

/// OCCT Adaptor3d_CurveOnSurface.cxx L72-78: `to3d(Pl, A)` —
/// gp_Ax2(P, VX.Crossed(VY), VX); the gp_Ax2 3-arg constructor
/// (gp_Ax2.hxx L73-80) orthogonalizes the X direction and derives Y.
fn to3d_ax22d(
    pl: &crate::geom::Plane,
    location: DVec2,
    x_dir: DVec2,
    y_dir: DVec2,
) -> super::elslib_iso::Ax2View {
    let p = to3d_pnt(pl, location);
    let vx = to3d_vec(pl, x_dir);
    let vy = to3d_vec(pl, y_dir);
    let n = vx.cross(vy).normalize_or_zero();
    super::elslib_iso::Ax2View {
        location: p,
        direction: n,
        // OCCT: vxdir.CrossCross(Vx, N) = N ^ (Vx ^ N).
        x_direction: n.cross(vx.cross(n)).normalize_or_zero(),
    }
}

/// OCCT Adaptor3d_CurveOnSurface.cxx L80-83: `to3d(Pl, C)` —
/// gp_Circ(to3d(Pl, C.Axis()), C.Radius()).
fn to3d_circ(pl: &crate::geom::Plane, c: crate::geom::Circle2d) -> Circle3 {
    super::elslib_iso::circ_from_ax2(&to3d_ax22d(pl, c.center, c.x_dir, c.y_dir), c.radius)
}

/// OCCT Adaptor3d_CurveOnSurface.cxx L85-88: `to3d(Pl, E)` —
/// gp_Elips(to3d(Pl, E.Axis()), E.MajorRadius(), E.MinorRadius()).
fn to3d_elips(pl: &crate::geom::Plane, e: crate::geom::Ellipse2d) -> Ellipse3 {
    let axes = to3d_ax22d(pl, e.center, e.major_dir, e.minor_dir);
    Ellipse3 {
        center: axes.location,
        normal: axes.direction,
        major_dir: axes.x_direction,
        major_radius: e.major_radius,
        minor_radius: e.minor_radius,
    }
}

/// OCCT Adaptor3d_CurveOnSurface.cxx L90-93: `to3d(Pl, H)` —
/// gp_Hypr(to3d(Pl, H.Axis()), H.MajorRadius(), H.MinorRadius()).  The rcad
/// Hyperbola2d carries no explicit minor direction, so the Ax22d Y direction
/// uses the codebase conic convention (-major_dir.y, major_dir.x).
fn to3d_hypr(pl: &crate::geom::Plane, h: crate::geom::Hyperbola2d) -> Hyperbola3 {
    let y_dir = DVec2::new(-h.major_dir.y, h.major_dir.x);
    let axes = to3d_ax22d(pl, h.center, h.major_dir, y_dir);
    Hyperbola3 {
        center: axes.location,
        normal: axes.direction,
        major_dir: axes.x_direction,
        semi_major: h.semi_major,
        semi_minor: h.semi_minor,
    }
}

/// OCCT Adaptor3d_CurveOnSurface.cxx L95-98: `to3d(Pl, P)` —
/// gp_Parab(to3d(Pl, P.Axis()), P.Focal()).
fn to3d_parab(pl: &crate::geom::Plane, p: crate::geom::Parabola2d) -> Parabola3 {
    let axes = to3d_ax22d(pl, p.origin, p.axis_dir, DVec2::new(-p.axis_dir.y, p.axis_dir.x));
    Parabola3 {
        vertex: axes.location,
        normal: axes.direction,
        axis_dir: axes.x_direction,
        focal_param: p.focal_param,
    }
}

// OCCT Adaptor3d_CurveOnSurface.cxx L1008-1041 — the static interval helper.

/// OCCT static AddIntervals (Adaptor3d_CurveOnSurface.cxx L1008-1041):
/// appends the roots of the equation to the sorted sequence of parameters
/// along the curve, keeping it sorted and avoiding repetitions (within
/// tolerance theTol).
fn add_intervals(the_parameters: &mut Vec<f64>, the_roots: &crate::math::root::FunctionRoots, the_tol: f64) {
    if !the_roots.is_done() || the_roots.is_all_null() {
        return;
    }

    let nsol = the_roots.nb_solutions();
    for i in 1..=nsol {
        let param = the_roots.value(i);
        if param - the_parameters[0] < the_tol {
            // skip param if equal to or less than theParameters(1)
            continue;
        }
        for j in 2..=the_parameters.len() {
            let a_delta = the_parameters[j - 1] - param;
            if a_delta > the_tol {
                the_parameters.insert(j - 1, param);
                break;
            } else if a_delta >= -the_tol {
                // param == theParameters(j) within Tol
                break;
            }
        }
    }
}

/// OCCT Adaptor3d_CurveOnSurface (TKG3d) — a 3D curve defined as the image of
/// a 2D curve on a surface.
///
/// The evaluation members are the OCCT definition: the curve point is the
/// surface evaluated at the 2D curve point, and the derivatives chain by the
/// chain rule (Adaptor3d_CurveOnSurface.cxx Value/D1/D2).  The Load pair
/// derives the analytic kind through EvalKPart
/// (Adaptor3d_CurveOnSurface.cxx L1552-1732) and caches the Line / Circle
/// payloads.
///
/// GAP (staged): the trimming support keeps the constructor-based encoding
/// (Load(mySurface) + Load(myCurve->Trim(...)), L1133-1141 — behaviour
/// equal), and Load(C) does not run EvalFirstLastSurf (L1736-1829): its
/// myFirstSurf / myLastSurf results feed only the boundary branches of the
/// OCCT EvalD1/D2/D3 (L1212-1341) which the rcad evaluation encoding does
/// not carry.
pub struct CurveOnSurface {
    /// OCCT: const handle(Adaptor2d_Curve2d) my2dCurve.
    pub my2d_curve: Curve2dHandle,
    /// OCCT: handle(Adaptor3d_Surface) mySurface.
    pub my_surface: SurfaceHandle,
    /// OCCT: GeomAbs_CurveType myType — set by the Load/EvalKPart pair; the
    /// constructor default is GeomAbs_OtherCurve
    /// (Adaptor3d_CurveOnSurface.cxx L890-897).
    pub my_type: CurveType,
    /// OCCT: gp_Circ myCirc — the analytic payload cached by EvalKPart
    /// (read through Circle(), L1390-1396); None models the non-circle state
    /// (the raise fires on the myType guard before the payload is read).
    pub my_circ: Option<Circle3>,
    /// OCCT: gp_Lin myLin — the analytic payload cached by EvalKPart (read
    /// through Line(), L1381-1387).
    pub my_lin: Option<Line3>,
    /// OCCT: GeomAbs_Shape myIntCont + handle(NCollection_HSequence<double>)
    /// myIntervals — the interval cache of NbIntervals (L1045-1115); the
    /// OCCT const_cast write maps to the mutex.
    my_intervals: std::sync::Mutex<Option<(GeomAbsShape, Vec<f64>)>>,
}

impl CurveOnSurface {
    /// OCCT Adaptor3d_CurveOnSurface(C, S) (Adaptor3d_CurveOnSurface.cxx
    /// L890-897): myType = GeomAbs_OtherCurve, myIntCont = GeomAbs_CN,
    /// Load(S), Load(C).
    pub fn new(the2d_curve: Curve2dHandle, the_surface: SurfaceHandle) -> Self {
        let mut cos = CurveOnSurface {
            my2d_curve: the2d_curve,
            my_surface: the_surface,
            my_type: CurveType::Other,
            my_circ: None,
            my_lin: None,
            my_intervals: std::sync::Mutex::new(None),
        };
        cos.load_surface();
        cos.load_curve();
        cos
    }

    /// OCCT GetCurve() — the 2D basis curve.
    pub fn get_curve(&self) -> &Curve2dHandle {
        &self.my2d_curve
    }

    /// OCCT GetSurface() — the surface.
    pub fn get_surface(&self) -> &SurfaceHandle {
        &self.my_surface
    }

    /// OCCT Adaptor3d_CurveOnSurface::Load(S)
    /// (Adaptor3d_CurveOnSurface.cxx L932-939).
    pub fn load_surface(&mut self) {
        // OCCT: mySurface = S; if (!myCurve.IsNull()) EvalKPart(); — the rcad
        // curve handle is non-optional.
        self.eval_k_part();
    }

    /// OCCT Adaptor3d_CurveOnSurface::Load(C)
    /// (Adaptor3d_CurveOnSurface.cxx L943-964).
    pub fn load_curve(&mut self) {
        // OCCT: myCurve = C; if (mySurface.IsNull()) return; — the rcad
        // surface handle is non-optional.

        self.eval_k_part();

        let mut s_type = self.my_surface.get_type();
        if s_type == GeomAbsSurfaceType::OffsetSurface {
            s_type = self.my_surface.basis_surface().get_type();
        }

        // OCCT L959-963: for the BSpline / extrusion / revolution surfaces
        // the Load runs EvalFirstLastSurf (L1736-1829) — GAP (staged), see
        // the struct docs; the branch is preserved as the recorded no-op.
        let _ = s_type;
    }

    /// OCCT Adaptor3d_CurveOnSurface::EvalKPart
    /// (Adaptor3d_CurveOnSurface.cxx L1552-1732) — derives the curve type and
    /// fills the analytic payloads for the plane-based and isoparametric
    /// line-on-quadric cases.
    pub fn eval_k_part(&mut self) {
        use super::elslib_iso::{
            circ_rotated, circ_with_direction_reversed, dir2d_is_opposite, dir2d_is_parallel,
            elslib_cone_u_iso, elslib_cone_v_iso, elslib_cylinder_u_iso, elslib_cylinder_v_iso,
            elslib_sphere_u_iso, elslib_sphere_v_iso, elslib_torus_u_iso, elslib_torus_v_iso,
            Ax3View,
        };

        // OCCT L1554.
        self.my_type = CurveType::Other;

        let s_ty = self.my_surface.get_type();
        let c_ty = self.my2d_curve.get_type();
        // OCCT L1558-1577: the plane branch.
        if s_ty == GeomAbsSurfaceType::Plane {
            self.my_type = c_ty;
            if self.my_type == CurveType::Circle {
                // OCCT L1563: myCirc = to3d(mySurface->Plane(), myCurve->Circle()).
                let pl = surface_plane(&*self.my_surface);
                let c2d = self.my2d_curve.circle();
                self.my_circ = Some(to3d_circ(&pl, c2d));
            } else if self.my_type == CurveType::Line {
                // OCCT L1565-1576.
                let (p_uv, d_uv) = self.my2d_curve.d1(0.0);
                let (p, d1u, d1v) = self.my_surface.d1(p_uv.x, p_uv.y);
                // OCCT: V.SetLinearForm(Duv.X(), D1U, Duv.Y(), D1V).
                let v = d1u * d_uv.x + d1v * d_uv.y;
                self.my_lin = Some(Line3::new(p, v));
            }
        } else if c_ty == CurveType::Line {
            // OCCT L1580-1582: gp_Dir2d D = myCurve->Line().Direction().
            let d = self.my2d_curve.line().direction;
            if dir2d_is_parallel(d, DVec2::X, precision::ANGULAR) {
                // OCCT L1584: Iso V.
                match s_ty {
                    GeomAbsSurfaceType::Sphere => {
                        // OCCT L1585-1604.
                        let p = self.my2d_curve.line().origin;
                        if ((p.y.abs() - std::f64::consts::FRAC_PI_2).abs())
                            >= precision::p_confusion()
                        {
                            self.my_type = CurveType::Circle;
                            let sph = surface_sphere(&*self.my_surface);
                            let axis = Ax3View::from_axes(sph.center, sph.axis, sph.ref_dir);
                            let mut a_circ = elslib_sphere_v_iso(&axis, sph.radius, p.y);
                            // OCCT: DRev = Axis.XDirection().Crossed(Axis.YDirection());
                            //       AxeRev(Axis.Location(), DRev); myCirc.Rotate(AxeRev, P.X()).
                            let d_rev = axis.x_direction.cross(axis.y_direction);
                            a_circ = circ_rotated(&a_circ, axis.location, d_rev, p.x);
                            if dir2d_is_opposite(d, DVec2::X, precision::ANGULAR) {
                                a_circ = circ_with_direction_reversed(&a_circ);
                            }
                            self.my_circ = Some(a_circ);
                        }
                    }
                    GeomAbsSurfaceType::Cylinder => {
                        // OCCT L1605-1621.
                        self.my_type = CurveType::Circle;
                        let cyl = surface_cylinder(&*self.my_surface);
                        let p = self.my2d_curve.line().origin;
                        let axis = ax3_view_of_cylinder(&cyl);
                        let mut a_circ = elslib_cylinder_v_iso(&axis, cyl.radius, p.y);
                        let d_rev = axis.x_direction.cross(axis.y_direction);
                        a_circ = circ_rotated(&a_circ, axis.location, d_rev, p.x);
                        if dir2d_is_opposite(d, DVec2::X, precision::ANGULAR) {
                            a_circ = circ_with_direction_reversed(&a_circ);
                        }
                        self.my_circ = Some(a_circ);
                    }
                    GeomAbsSurfaceType::Cone => {
                        // OCCT L1622-1638.
                        self.my_type = CurveType::Circle;
                        let cone = surface_cone(&*self.my_surface);
                        let p = self.my2d_curve.line().origin;
                        let axis = Ax3View::from_axes(cone.apex, cone.axis, cone.ref_dir);
                        let mut a_circ =
                            elslib_cone_v_iso(&axis, cone.radius, cone.half_angle_rad, p.y);
                        let d_rev = axis.x_direction.cross(axis.y_direction);
                        a_circ = circ_rotated(&a_circ, axis.location, d_rev, p.x);
                        if dir2d_is_opposite(d, DVec2::X, precision::ANGULAR) {
                            a_circ = circ_with_direction_reversed(&a_circ);
                        }
                        self.my_circ = Some(a_circ);
                    }
                    GeomAbsSurfaceType::Torus => {
                        // OCCT L1639-1655.
                        self.my_type = CurveType::Circle;
                        let tore = surface_torus(&*self.my_surface);
                        let p = self.my2d_curve.line().origin;
                        let axis = Ax3View::from_axes(tore.center, tore.axis, tore.ref_dir);
                        let mut a_circ =
                            elslib_torus_v_iso(&axis, tore.major_radius, tore.minor_radius, p.y);
                        let d_rev = axis.x_direction.cross(axis.y_direction);
                        a_circ = circ_rotated(&a_circ, axis.location, d_rev, p.x);
                        if dir2d_is_opposite(d, DVec2::X, precision::ANGULAR) {
                            a_circ = circ_with_direction_reversed(&a_circ);
                        }
                        self.my_circ = Some(a_circ);
                    }
                    _ => {}
                }
            } else if dir2d_is_parallel(d, DVec2::Y, precision::ANGULAR) {
                // OCCT L1657: Iso U.
                match s_ty {
                    GeomAbsSurfaceType::Sphere => {
                        // OCCT L1659-1684.
                        self.my_type = CurveType::Circle;
                        let sph = surface_sphere(&*self.my_surface);
                        let p = self.my2d_curve.line().origin;
                        let axis = Ax3View::from_axes(sph.center, sph.axis, sph.ref_dir);
                        // OCCT L1666: compute the iso 0.
                        let mut a_circ = elslib_sphere_u_iso(&axis, sph.radius, 0.0);
                        // OCCT L1669-1671: same-parametrization (circle
                        // rotation - Y offset); DRev = Axis.XDirection()
                        // .Crossed(Axis.Direction()).
                        let d_rev = axis.x_direction.cross(axis.direction);
                        a_circ = circ_rotated(&a_circ, axis.location, d_rev, p.y);
                        // OCCT L1674-1676: transform to iso U (= P.X());
                        // DRev = Axis.XDirection().Crossed(Axis.YDirection()).
                        let d_rev = axis.x_direction.cross(axis.y_direction);
                        a_circ = circ_rotated(&a_circ, axis.location, d_rev, p.x);
                        if dir2d_is_opposite(d, DVec2::Y, precision::ANGULAR) {
                            a_circ = circ_with_direction_reversed(&a_circ);
                        }
                        self.my_circ = Some(a_circ);
                    }
                    GeomAbsSurfaceType::Cylinder => {
                        // OCCT L1685-1698.
                        self.my_type = CurveType::Line;
                        let cyl = surface_cylinder(&*self.my_surface);
                        let p = self.my2d_curve.line().origin;
                        let axis = ax3_view_of_cylinder(&cyl);
                        let mut a_lin = elslib_cylinder_u_iso(&axis, cyl.radius, p.x);
                        // OCCT: Tr(myLin.Direction()); Tr.Multiply(P.Y());
                        //       myLin.Translate(Tr).
                        a_lin.origin += a_lin.direction * p.y;
                        if dir2d_is_opposite(d, DVec2::Y, precision::ANGULAR) {
                            // OCCT: myLin.Reverse().
                            a_lin.direction = -a_lin.direction;
                        }
                        self.my_lin = Some(a_lin);
                    }
                    GeomAbsSurfaceType::Cone => {
                        // OCCT L1699-1712.
                        self.my_type = CurveType::Line;
                        let cone = surface_cone(&*self.my_surface);
                        let p = self.my2d_curve.line().origin;
                        let axis = Ax3View::from_axes(cone.apex, cone.axis, cone.ref_dir);
                        let mut a_lin =
                            elslib_cone_u_iso(&axis, cone.radius, cone.half_angle_rad, p.x);
                        a_lin.origin += a_lin.direction * p.y;
                        if dir2d_is_opposite(d, DVec2::Y, precision::ANGULAR) {
                            a_lin.direction = -a_lin.direction;
                        }
                        self.my_lin = Some(a_lin);
                    }
                    GeomAbsSurfaceType::Torus => {
                        // OCCT L1713-1728.
                        self.my_type = CurveType::Circle;
                        let tore = surface_torus(&*self.my_surface);
                        let p = self.my2d_curve.line().origin;
                        let axis = Ax3View::from_axes(tore.center, tore.axis, tore.ref_dir);
                        let a_circ =
                            elslib_torus_u_iso(&axis, tore.major_radius, tore.minor_radius, p.x);
                        // OCCT L1720: myCirc.Rotate(myCirc.Axis(), P.Y()) —
                        // the circle's own axis (Location = center,
                        // Direction = normal).
                        let mut a_circ = circ_rotated(&a_circ, a_circ.center, a_circ.normal, p.y);
                        if dir2d_is_opposite(d, DVec2::Y, precision::ANGULAR) {
                            a_circ = circ_with_direction_reversed(&a_circ);
                        }
                        self.my_circ = Some(a_circ);
                    }
                    _ => {}
                }
            }
        }
    }
}

/// The gp_Ax3 frame view of a cylinder payload — the explicit-Y (possibly
/// left-handed) swept-lateral frame is honored (see
/// `CylindricalSurface::y_dir`).
fn ax3_view_of_cylinder(cyl: &crate::geom::CylindricalSurface) -> super::elslib_iso::Ax3View {
    match cyl.y_dir {
        Some(y) => super::elslib_iso::Ax3View::with_y_dir(cyl.origin, cyl.axis, cyl.ref_dir, y),
        None => super::elslib_iso::Ax3View::from_axes(cyl.origin, cyl.axis, cyl.ref_dir),
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

    /// OCCT Adaptor3d_CurveOnSurface::NbIntervals(S)
    /// (Adaptor3d_CurveOnSurface.cxx L1045-1115): the cached sequence for the
    /// recorded continuity, else the sorted union of the curve intervals and
    /// the surface-discontinuity crossings (the Adaptor3d_InterFunc +
    /// math_FunctionRoots solve).
    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        // OCCT L1047-1050.
        {
            let cached = self.my_intervals.lock().unwrap();
            if let Some((cached_cont, cached_intervals)) = cached.as_ref() {
                if *cached_cont == s {
                    return cached_intervals.len() - 1;
                }
            }
        }

        let nu = self.my_surface.nb_u_intervals(s);
        let nv = self.my_surface.nb_v_intervals(s);
        let nc = self.my2d_curve.nb_intervals(s);

        // OCCT L1058-1061: the TabU / TabV / TabC arrays — the rcad
        // Intervals() calls return the filled arrays directly.
        let nb_sample: i32 = 20;
        let tdeb = self.my2d_curve.first_parameter();
        let tfin = self.my2d_curve.last_parameter();

        // OCCT L1068: myCurve->Intervals(TabC, S).
        let tab_c = self.my2d_curve.intervals(s);

        let tol = precision::p_confusion() / 10.0; // OCCT L1070.

        // OCCT L1072-1079: the sorted sequence of parameters defining
        // continuity intervals; started with own intervals of curve and
        // completed by additional points coming from surface
        // discontinuities.
        let mut a_intervals: Vec<f64> = Vec::with_capacity(nc + 1);
        for i in 1..=nc + 1 {
            a_intervals.push(tab_c[i - 1]);
        }

        // OCCT L1081-1091.
        if nu > 1 {
            let tab_u = self.my_surface.u_intervals(s);
            for iu in 2..=nu {
                let u = tab_u[iu - 1];
                let mut func =
                    super::adaptor3d_interfunc::Adaptor3dInterFunc::new(self.my2d_curve.clone(), u, 1);
                let resol = crate::math::root::FunctionRoots::new(
                    &mut func, tdeb, tfin, nb_sample, tol, tol, tol, 0.0,
                );
                add_intervals(&mut a_intervals, &resol, tol);
            }
        }
        // OCCT L1092-1102.
        if nv > 1 {
            let tab_v = self.my_surface.v_intervals(s);
            for iv in 2..=nv {
                let v = tab_v[iv - 1];
                let mut func =
                    super::adaptor3d_interfunc::Adaptor3dInterFunc::new(self.my2d_curve.clone(), v, 2);
                let resol = crate::math::root::FunctionRoots::new(
                    &mut func, tdeb, tfin, nb_sample, tol, tol, tol, 0.0,
                );
                add_intervals(&mut a_intervals, &resol, tol);
            }
        }

        // OCCT L1104-1110: for case intervals==1 and first point == last
        // point SequenceOfReal contains only one value, therefore it is
        // necessary to add second value into aIntervals which will be equal
        // first value.
        if a_intervals.len() == 1 {
            let v = a_intervals[0];
            a_intervals.push(v);
        }

        // OCCT L1112-1114: the const_cast cache write.
        let n = a_intervals.len() - 1;
        *self.my_intervals.lock().unwrap() = Some((s, a_intervals));
        n
    }

    /// OCCT Adaptor3d_CurveOnSurface::Intervals(T, S)
    /// (Adaptor3d_CurveOnSurface.cxx L1119-1129): the cached sequence of
    /// NbIntervals(S) + 1 bounds.  The OCCT Standard_ASSERT_RAISE on the
    /// caller buffer length has no counterpart in the Vec-returning rcad
    /// encoding (the returned sequence IS the filled array).
    fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        self.nb_intervals(s); // OCCT L1121.
        let cached = self.my_intervals.lock().unwrap();
        match cached.as_ref() {
            Some((_, v)) => v.clone(),
            None => unreachable!("Adaptor3d_CurveOnSurface::Intervals: NbIntervals must cache"),
        }
    }

    /// OCCT Adaptor3d_CurveOnSurface::Trim
    /// (Adaptor3d_CurveOnSurface.cxx L1133-1141): the rebuilt adaptor runs
    /// Load(mySurface) + Load(myCurve->Trim(First, Last, Tol)) — the rcad
    /// constructor performs the same Load pair.
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
        if self.my_type != CurveType::Line {
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Line(): curve is not a line");
        }
        self.my_lin.expect("EvalKPart must cache myLin for GeomAbs_Line")
    }

    /// OCCT Adaptor3d_CurveOnSurface::Circle (L1390-1396) — the myCirc
    /// payload cached by EvalKPart; raises when myType != GeomAbs_Circle.
    fn circle(&self) -> Circle3 {
        if self.my_type != CurveType::Circle {
            panic!("Standard_NoSuchObject: Adaptor3d_CurveOnSurface::Line(): curve is not a circle");
        }
        self.my_circ.expect("EvalKPart must cache myCirc for GeomAbs_Circle")
    }

    /// OCCT Adaptor3d_CurveOnSurface::Ellipse (L1399-1403) —
    /// to3d(mySurface->Plane(), myCurve->Ellipse()).
    fn ellipse(&self) -> Ellipse3 {
        let pl = surface_plane(&*self.my_surface);
        to3d_elips(&pl, self.my2d_curve.ellipse())
    }

    /// OCCT Adaptor3d_CurveOnSurface::Hyperbola (L1406-1410) —
    /// to3d(mySurface->Plane(), myCurve->Hyperbola()).
    fn hyperbola(&self) -> Hyperbola3 {
        let pl = surface_plane(&*self.my_surface);
        to3d_hypr(&pl, self.my2d_curve.hyperbola())
    }

    /// OCCT Adaptor3d_CurveOnSurface::Parabola (L1413-1417) —
    /// to3d(mySurface->Plane(), myCurve->Parabola()).
    fn parabola(&self) -> Parabola3 {
        let pl = surface_plane(&*self.my_surface);
        to3d_parab(&pl, self.my2d_curve.parabola())
    }

    /// OCCT Adaptor3d_CurveOnSurface::Trim — see the Adaptor3dCurve trim.
    fn trim_geom(&self, first: f64, last: f64, tol: f64) -> Arc<dyn Adaptor3dCurveGeom> {
        let hcs = CurveOnSurface::new(self.my2d_curve.trim(first, last, tol), self.my_surface.clone());
        Arc::new(hcs)
    }
}
