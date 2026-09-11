//! OCCT Extrema_CurveTool (TKGeomBase/Extrema/Extrema_CurveTool.hxx L38-142) —
//! the static facade over `Adaptor3d_Curve` used by the whole Extrema package
//! (`Extrema_ExtCC`, `Extrema_ExtElC`, `Extrema_GGenExtCC`, `Extrema_ExtCS`).
//!
//! rcad encoding (architecture difference, annotated):
//! OCCT `Adaptor3d_Curve` is one class carrying both the evaluation virtuals
//! and the type accessors (GetType/Line/Circle/.../IsPeriodic/Period/
//! Resolution/IsClosed).  The rcad kernel splits that class over two traits:
//! `base::proj_lib::adaptor::Adaptor3dCurve` (evaluation, interval machinery)
//! and `base::proj_lib::proj_lib_projected_curve::Adaptor3dCurveGeom`
//! (GetType + the analytic downcasts), and neither carries
//! IsPeriodic/Period/Resolution/IsClosed.
//!
//! `ExtremaCurveTool` restores the OCCT shape: it is the rcad translation of
//! `Extrema_CurveTool`, i.e. the facade the Extrema algorithms call instead of
//! touching the adaptor directly.  [`CurveToolHandle`] is the adapter that
//! routes the facade onto an rcad adaptor handle; the queries the rcad adaptor
//! traits do not expose (periodicity, period, resolution, closedness) are
//! supplied by the caller at construction — the OCCT value is filled in by the
//! caller that owns the underlying `Geom_Curve` (see
//! `proj_lib::proj_lib_projected_curve`'s `GeomCurveAdaptor` companion).
//!
//! Interface request (see the delivery report): enrich
//! `base::proj_lib::adaptor::Adaptor3dCurve` with `get_type`, `is_periodic`,
//! `period`, `is_closed`, `resolution` and the `Adaptor3dCurveGeom` downcasts so
//! this facade collapses back into the occt::handle call sites.

use glam::DVec3;

use crate::base::extrema_ext_pc::BSplineView;
use crate::base::proj_lib::adaptor::Adaptor3dCurve;
use crate::base::proj_lib::proj_lib_projected_curve::Adaptor3dCurveGeom;
use crate::base::proj_lib::CurveType;
use crate::geom::{Circle3, Ellipse3, Hyperbola3, Line3, Parabola3};
use crate::math::GeomAbsShape;

/// OCCT Extrema_CurveTool (Extrema_CurveTool.hxx L38-142) — the static facade.
/// In OCCT every member is a one-liner forwarding to `Adaptor3d_Curve`; here
/// the facade is the trait itself so the Extrema algorithms can be written
/// against a single, tool-shaped interface exactly as the OCCT templates are.
pub trait ExtremaCurveTool {
    /// OCCT Extrema_CurveTool::FirstParameter (hxx L43).
    fn first_parameter(&self) -> f64;
    /// OCCT Extrema_CurveTool::LastParameter (hxx L45).
    fn last_parameter(&self) -> f64;
    /// OCCT Extrema_CurveTool::Continuity (hxx L47).
    fn continuity(&self) -> GeomAbsShape;
    /// OCCT Extrema_CurveTool::NbIntervals (hxx L51-54).
    fn nb_intervals(&self, s: GeomAbsShape) -> i32;
    /// OCCT Extrema_CurveTool::Intervals (hxx L59-64).
    fn intervals(&self, s: GeomAbsShape) -> Vec<f64>;
    /// OCCT Extrema_CurveTool::IsPeriodic (hxx L71).
    fn is_periodic(&self) -> bool;
    /// OCCT Extrema_CurveTool::Period (hxx L73).
    fn period(&self) -> f64;
    /// OCCT Extrema_CurveTool::Resolution (hxx L75-78).
    fn resolution(&self, r3d: f64) -> f64;
    /// OCCT Extrema_CurveTool::GetType (hxx L80).
    fn get_type(&self) -> CurveType;
    /// OCCT `Adaptor3d_Curve::IsClosed()`.
    fn is_closed(&self) -> bool;
    /// OCCT Extrema_CurveTool::Value (hxx L82).
    fn value(&self, u: f64) -> DVec3;
    /// OCCT Extrema_CurveTool::D0 (hxx L84-87).
    fn d0(&self, u: f64) -> DVec3 {
        self.value(u)
    }
    /// OCCT Extrema_CurveTool::D1 (hxx L89-92).
    fn d1(&self, u: f64) -> (DVec3, DVec3);
    /// OCCT Extrema_CurveTool::D2 (hxx L94-101).
    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3);
    /// OCCT Extrema_CurveTool::DN (hxx L113-116) — the `Adaptor3d_Curve::DN`
    /// default implementation dispatches on the order.
    fn dn(&self, u: f64, n: i32) -> DVec3 {
        match n {
            1 => self.d1(u).1,
            2 => self.d2(u).2,
            _ => panic!("Standard_NotImplemented: Adaptor3d_Curve::DN"),
        }
    }
    /// OCCT Extrema_CurveTool::Line (hxx L118).
    fn line(&self) -> Line3;
    /// OCCT Extrema_CurveTool::Circle (hxx L120).
    fn circle(&self) -> Circle3;
    /// OCCT Extrema_CurveTool::Ellipse (hxx L122).
    fn ellipse(&self) -> Ellipse3;
    /// OCCT Extrema_CurveTool::Hyperbola (hxx L124).
    fn hyperbola(&self) -> Hyperbola3;
    /// OCCT Extrema_CurveTool::Parabola (hxx L126).
    fn parabola(&self) -> Parabola3;
    /// OCCT `Extrema_CurveTool::Bezier(theC)` (hxx L136) — the
    /// `Geom_BezierCurve` pole count `Extrema_GGExtPC::Perform` samples with
    /// (`mysample = TheCurveTool::Bezier(aCurve)->NbPoles() * 2`,
    /// Extrema_GGExtPC.hxx L184-186).
    fn bezier_nb_poles(&self) -> usize;
    /// OCCT `Extrema_CurveTool::BSpline(theC)` (hxx L138) — the
    /// `Geom_BSplineCurve` knot/degree view `Extrema_GGExtPC::Perform` walks
    /// (Extrema_GGExtPC.hxx L190-192).
    fn bspline(&self) -> BSplineView;
    /// OCCT `GCPnts_AbscissaPoint::Length(C)` (CPnts_AbscissaPoint.cxx
    /// L103-106) — the arc length of the adaptor over its whole domain.
    ///
    /// rcad encoding: answers `-1` when the adaptor carries no kernel curve
    /// (the arc-length engine is Curve3-based); `Extrema_GGenExtCC::Perform`
    /// then keeps `indmax == -1`, the state OCCT reaches in the same place
    /// whenever a curve parameter is infinite
    /// (Extrema_GGenExtCC.hxx L476-491).
    fn abscissa_length(&self) -> f64 {
        -1.0
    }
}

/// The rcad adapter routing [`ExtremaCurveTool`] onto an rcad adaptor handle
/// (the OCCT `Adaptor3d_Curve&` a tool call receives).
pub struct CurveToolHandle<'a> {
    /// The evaluation/interval half of the adaptor
    /// (`Adaptor3d_Curve`).
    pub curve: &'a dyn Adaptor3dCurve,
    /// The type/analytic-downcast half of the adaptor
    /// (`Adaptor3dCurveGeom`), when the caller holds it.
    pub geom: Option<&'a dyn Adaptor3dCurveGeom>,
    /// OCCT `Adaptor3d_Curve::GetType()` when `geom` is absent.
    pub curve_type: CurveType,
    /// OCCT `Adaptor3d_Curve::IsPeriodic()`.
    pub is_periodic: bool,
    /// OCCT `Adaptor3d_Curve::Period()`.
    pub period: f64,
    /// OCCT `Adaptor3d_Curve::Resolution(R3d)`.
    pub resolution: f64,
    /// OCCT `Adaptor3d_Curve::IsClosed()`.
    pub is_closed: bool,
    /// The kernel curve the adaptor wraps, when the caller holds it — used by
    /// `GCPnts_AbscissaPoint::Length(Adaptor3d_Curve)` through the kernel
    /// `base::gcpnts::abscissa_point::arc_length` engine.
    pub curve3: Option<&'a crate::geom::Curve3>,
}

impl<'a> CurveToolHandle<'a> {
    /// The facade over an adaptor handle that also carries the type half.
    pub fn with_geom(curve: &'a dyn Adaptor3dCurve, geom: &'a dyn Adaptor3dCurveGeom) -> Self {
        CurveToolHandle {
            curve,
            geom: Some(geom),
            curve_type: geom.get_type(),
            is_periodic: false,
            period: 0.0,
            resolution: 1.0,
            is_closed: false,
            curve3: None,
        }
    }

    /// The facade over the rcad encoding of OCCT `GeomAdaptor_Curve(C, First,
    /// Last)`: the adaptor half answers the evaluation, the geometry half the
    /// type/downcast, and the OCCT queries the rcad adaptor traits do not
    /// carry (`Geom_Curve::IsPeriodic` / `Period` / `IsClosed`) are taken from
    /// the kernel curve itself.
    pub fn for_curve3(
        curve3: &'a crate::geom::Curve3,
        curve: &'a dyn Adaptor3dCurve,
        geom: &'a dyn Adaptor3dCurveGeom,
    ) -> Self {
        // OCCT Geom_Circle::IsPeriodic() / Period() / IsClosed().
        let a_period = match curve3 {
            crate::geom::Curve3::Circle(_) | crate::geom::Curve3::Ellipse(_) => {
                std::f64::consts::PI + std::f64::consts::PI
            }
            _ => 0.0,
        };
        CurveToolHandle {
            curve,
            geom: Some(geom),
            curve_type: geom.get_type(),
            is_periodic: a_period > 0.0,
            period: a_period,
            resolution: 1.0,
            is_closed: a_period > 0.0,
            curve3: Some(curve3),
        }
    }

    /// The facade over an evaluation-only adaptor handle: the OCCT queries the
    /// rcad `Adaptor3dCurve` trait does not carry are supplied by the caller.
    pub fn new(
        curve: &'a dyn Adaptor3dCurve,
        curve_type: CurveType,
        is_periodic: bool,
        period: f64,
        resolution: f64,
        is_closed: bool,
    ) -> Self {
        CurveToolHandle {
            curve,
            geom: None,
            curve_type,
            is_periodic,
            period,
            resolution,
            is_closed,
            curve3: None,
        }
    }

    /// The default facade for an evaluation-only handle whose extra queries are
    /// unknown — every query answers the OCCT value of an adaptor that is not
    /// periodic, not closed and has unit resolution, and the type is the
    /// `GeomAbs_OtherCurve` OCCT reports for a generic adaptor
    /// (`Adaptor3d_CurveOnSurface::GetType()`).
    pub fn other(curve: &'a dyn Adaptor3dCurve) -> Self {
        CurveToolHandle::new(curve, CurveType::Other, false, 0.0, 1.0, false)
    }

    /// OCCT `GCPnts_AbscissaPoint::Length(C)` (CPnts_AbscissaPoint.cxx L103-106)
    /// — the arc length of the adaptor over its whole domain.
    ///
    /// rcad encoding: when the adaptor wraps a kernel curve the kernel
    /// `base::gcpnts::abscissa_point::arc_length` engine computes the length
    /// (the existing rcad GCPnts machinery, reused not retranslated); when it
    /// does not, the query answers the sentinel `-1` so the caller's
    /// "lengths are unavailable" branch runs — the branch OCCT itself takes
    /// whenever a parameter is infinite (`Extrema_GGenExtCC::Perform`
    /// cxx L476-491).
    pub fn abscissa_length(&self) -> f64 {
        match self.curve3 {
            Some(c) => crate::base::gcpnts::abscissa_point::arc_length(
                c,
                self.first_parameter(),
                self.last_parameter(),
            ),
            None => -1.0,
        }
    }
}

impl ExtremaCurveTool for CurveToolHandle<'_> {
    fn first_parameter(&self) -> f64 {
        self.curve.first_parameter()
    }

    fn last_parameter(&self) -> f64 {
        self.curve.last_parameter()
    }

    fn continuity(&self) -> GeomAbsShape {
        self.curve.continuity()
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> i32 {
        self.curve.nb_intervals(s) as i32
    }

    fn intervals(&self, s: GeomAbsShape) -> Vec<f64> {
        self.curve.intervals(s)
    }

    fn is_periodic(&self) -> bool {
        self.is_periodic
    }

    fn period(&self) -> f64 {
        self.period
    }

    fn resolution(&self, _r3d: f64) -> f64 {
        self.resolution
    }

    fn get_type(&self) -> CurveType {
        match self.geom {
            Some(g) => g.get_type(),
            None => self.curve_type,
        }
    }

    fn is_closed(&self) -> bool {
        self.is_closed
    }

    fn value(&self, u: f64) -> DVec3 {
        self.curve.value(u)
    }

    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        self.curve.d1(u)
    }

    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        self.curve.d2(u)
    }

    fn dn(&self, u: f64, n: i32) -> DVec3 {
        self.curve_dn(u, n)
    }

    fn line(&self) -> Line3 {
        self.geom.expect("Extrema_CurveTool: Line() on a non-line adaptor").line()
    }

    fn circle(&self) -> Circle3 {
        self.geom.expect("Extrema_CurveTool: Circle() on a non-circle adaptor").circle()
    }

    fn ellipse(&self) -> Ellipse3 {
        self.geom.expect("Extrema_CurveTool: Ellipse() on a non-ellipse adaptor").ellipse()
    }

    fn hyperbola(&self) -> Hyperbola3 {
        self.geom.expect("Extrema_CurveTool: Hyperbola() on a non-hyperbola adaptor").hyperbola()
    }

    fn parabola(&self) -> Parabola3 {
        self.geom.expect("Extrema_CurveTool: Parabola() on a non-parabola adaptor").parabola()
    }

    /// OCCT `Extrema_CurveTool::Bezier(theC)` (hxx L136) — `theC.Bezier()`,
    /// i.e. `Adaptor3d_Curve::Bezier()`, which raises
    /// `Standard_NoSuchObject` for a non-Bezier adaptor.
    fn bezier_nb_poles(&self) -> usize {
        match self.curve3 {
            Some(crate::geom::Curve3::Bezier(b)) => b.control_points.len(),
            _ => panic!("Standard_NoSuchObject: Extrema_CurveTool::Bezier on a non-bezier adaptor"),
        }
    }

    /// OCCT `Extrema_CurveTool::BSpline(theC)` (hxx L138) — `theC.BSpline()`,
    /// i.e. `Adaptor3d_Curve::BSpline()`, which raises
    /// `Standard_NoSuchObject` for a non-BSpline adaptor.
    fn bspline(&self) -> BSplineView {
        match self.curve3 {
            Some(crate::geom::Curve3::BSpline(bs)) => {
                // Reduced knot view out of the flat (multiplicity-expanded)
                // vector; see the BSplineView architecture-glue note.
                let mut knots = vec![bs.knots[0]];
                for &k in &bs.knots[1..] {
                    if (k - knots[knots.len() - 1]).abs() > 1e-12 {
                        knots.push(k);
                    }
                }
                BSplineView {
                    first_u_knot_index: 1,
                    last_u_knot_index: knots.len(),
                    knots,
                    degree: bs.degree,
                }
            }
            _ => panic!("Standard_NoSuchObject: Extrema_CurveTool::BSpline on a non-bspline adaptor"),
        }
    }
}

impl CurveToolHandle<'_> {
    /// OCCT `Adaptor3d_Curve::DN` inherited implementation (the adaptor
    /// dispatches on the derivative order; higher orders are not implemented).
    fn curve_dn(&self, u: f64, n: i32) -> DVec3 {
        match n {
            1 => self.curve.d1(u).1,
            2 => self.curve.d2(u).2,
            _ => panic!("Standard_NotImplemented: Adaptor3d_Curve::DN"),
        }
    }
}
