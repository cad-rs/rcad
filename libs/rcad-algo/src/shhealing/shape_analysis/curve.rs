//! OCCT ShapeAnalysis package class (TKShHealing): `ShapeAnalysis_Curve`
//! (`ShapeAnalysis_Curve.hxx` + `ShapeAnalysis_Curve.cxx` L1-1490).
//!
//! Tool for analysing curves: point projection (grid + Newton), parameter
//! range validation, pcurve bounding boxes, seam selection, planarity and
//! periodicity queries.
//!
//! Architecture bridges (numbering follows the W1 edge.rs style):
//! 1. OCCT `occ::handle<Geom_Curve>` / `handle<Geom2d_Curve>` -> `Curve3` /
//!    `Curve2d` values (clone-on-handle). `IsKind(Geom_BoundedCurve)` -> the
//!    BSpline/Bezier/Trimmed match; `IsClosed()`/`IsPeriodic()` -> the
//!    kernel `CurveEval` re-hosts.
//! 2. OCCT `Adaptor3d_Curve` (the `Project(Adaptor3d_Curve)` overload) ->
//!    the local [`Adaptor3dCurve`] trait; the concrete instances are the
//!    `brep_lib_validate_edge` re-hosts (`GeomAdaptorCurve`,
//!    `Adaptor3dCurveOnSurface`).
//! 3. `GeomAdaptor_Curve(C, u1, u2)` + `Load(C, u1, u2)` -> the immutable
//!    re-host re-construction (the rcad adaptor is a value).
//! 4. `Extrema_ExtPC(P, C)` / `Extrema_LocateExtPC(P, C, U0, Umin, Usup,
//!    TolF)` -> the kernel real `base::extrema_ext_pc::ExtremaExtPC` /
//!    `base::extrema_locate_ext_pc::LocateExtPC` over the kernel
//!    adaptor-stack tool (`GeomCurveAdaptor` / `CurveOnSurface` seen through
//!    the `CurveToolHandle` view — the concrete adaptor instance the OCCT
//!    ctors receive); for the re-host carrying no routed parts (the
//!    `CurveOnPlane` stand-in) the GAP below applies and the OCCT
//!    `!IsDone()` fallback (ProjectOnSegments, the `Project` call of
//!    NextProject) stays reachable. `Extrema_ExtPC::IsMin(i)` -> the kernel
//!    ExtremaExtPC carries the minima (the IsMin filter of ProjectAct).
//! 5. `ElCLib` -> the kernel `math::el` + `geomalgo::int_patch::elclib`
//!    re-hosts; `ElCLib::AdjustPeriodic` -> the local re-host below (the
//!    kernel topo/brep_lib/make_edge2d.rs precedent).
//! 6. `Geom2dAdaptor_Curve` / `Geom2dInt_Geom2dCurveTool` -> the
//!    `geomalgo::geom2d_int::Curve2dAdaptor` re-host; the C2 interval count
//!    of FillBndBox -> the local re-host below (the flat-knot walk of the
//!    brep_fill_sweep_c.rs precedent).
//! 7. `ShapeExtend_ComplexCurve` (IsPlanar arm) -> GAP: not translated yet
//!    (W1-3 delivered the status/MsgRegistrator/WireData side only); the
//!    OCCT failure path (`return false`) is preserved.
//! 8. `Standard_NoSuchObject` raises of the adaptor accessors -> GAP on the
//!    `Adaptor3d_CurveOnSurface` conic lifting (`EvalKPart`,
//!    Adaptor3d_CurveOnSurface.cxx): the rcad `curve3d()` returns None and
//!    the ProjectAct conic arms panic with the GAP anchor (the OCCT raise
//!    propagates out of the switch the same way).

use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::base::extrema_curve_tool::CurveToolHandle;
use rcad_kernel::base::extrema_ext_pc::ExtremaExtPC;
use rcad_kernel::base::extrema_locate_ext_pc::LocateExtPC;
use rcad_kernel::base::proj_lib::geom_adaptor_curve::GeomCurveAdaptor;
use rcad_kernel::base::proj_lib::geom_adaptor_surface::GeomSurfaceAdaptor;
use rcad_kernel::base::proj_lib::{CurveOnSurface, Geom2dCurveAdaptor};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Surface3};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::precision::{is_infinite_value, CONFUSION, INFINITE_VALUE, PCONFUSION};

use crate::geomalgo::geom2d_int::{curve2d_type_of, Curve2dAdaptor, Curve2dType};
use crate::topalgo::brep_lib_validate_edge::{Adaptor3dCurveOnSurface, GeomAdaptorCurve};

// OCCT Standard_Real.hxx L176-179 / L182-185.
const REAL_LAST: f64 = f64::MAX;
// OCCT gp.hxx L59-60: gp::Resolution() = RealSmall() = DBL_MIN.
const REAL_SMALL: f64 = f64::MIN_POSITIVE;

// ---------------------------------------------------------------------------
// Adaptor3d_Curve trait (architecture bridge #2)
// ---------------------------------------------------------------------------

/// OCCT GeomAbs_CurveType — the `Adaptor3d_Curve::GetType()` dispatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdaptorCurveKind {
    Line,
    Circle,
    Ellipse,
    Hyperbola,
    Parabola,
    Bezier,
    BSpline,
    Offset,
    Other,
}

/// OCCT Adaptor3d_Curve — the members ShapeAnalysis_Curve consumes.
/// `curve3d()` carries the per-type accessors (`Line()`/`Circle()`/...);
/// None models the Standard_NoSuchObject raise (bridge #8).
pub trait Adaptor3dCurve: CurveEval {
    /// OCCT FirstParameter().
    fn first_parameter(&self) -> f64;
    /// OCCT LastParameter().
    fn last_parameter(&self) -> f64;
    /// OCCT Value(U).
    fn value(&self, the_u: f64) -> DVec3;
    /// OCCT Resolution(R3d).
    fn resolution(&self, the_r3d: f64) -> f64;
    /// OCCT IsClosed().
    fn is_closed(&self) -> bool;
    /// OCCT IsKind(STANDARD_TYPE(Geom_BoundedCurve)).
    fn is_kind_bounded(&self) -> bool;
    /// OCCT GetType().
    fn get_type(&self) -> AdaptorCurveKind;
    /// OCCT Line()/Circle()/Ellipse()/Hyperbola()/Parabola()/BSpline()/
    /// Bezier() — the underlying curve value when the adaptor wraps one
    /// (None = the Standard_NoSuchObject raise, bridge #8).
    fn curve3d(&self) -> Option<Curve3>;
    /// The (pcurve, surface) pair of the `Adaptor3d_CurveOnSurface`
    /// instance — the concrete adaptor shape the OCCT Extrema ctors receive
    /// (bridge #4 tool routing below). The default None keeps the
    /// `GeomAdaptor_Curve`-like re-hosts on the `curve3d()` route.
    fn curve_on_surface_parts(&self) -> Option<(Curve2d, Surface3)> {
        None
    }
}

impl Adaptor3dCurve for GeomAdaptorCurve {
    fn first_parameter(&self) -> f64 {
        GeomAdaptorCurve::first_parameter(self)
    }
    fn last_parameter(&self) -> f64 {
        GeomAdaptorCurve::last_parameter(self)
    }
    fn value(&self, the_u: f64) -> DVec3 {
        GeomAdaptorCurve::value(self, the_u)
    }
    fn resolution(&self, the_r3d: f64) -> f64 {
        GeomAdaptorCurve::resolution(self, the_r3d)
    }
    fn is_closed(&self) -> bool {
        CurveEval::is_closed(self.curve())
    }
    fn is_kind_bounded(&self) -> bool {
        matches!(
            self.curve(),
            Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)
        )
    }
    fn get_type(&self) -> AdaptorCurveKind {
        match self.curve() {
            Curve3::Line(_) => AdaptorCurveKind::Line,
            Curve3::Circle(_) => AdaptorCurveKind::Circle,
            Curve3::Ellipse(_) => AdaptorCurveKind::Ellipse,
            Curve3::Hyperbola(_) => AdaptorCurveKind::Hyperbola,
            Curve3::Parabola(_) => AdaptorCurveKind::Parabola,
            Curve3::Bezier(_) => AdaptorCurveKind::Bezier,
            Curve3::BSpline(_) => AdaptorCurveKind::BSpline,
            Curve3::Offset(_) => AdaptorCurveKind::Offset,
            _ => AdaptorCurveKind::Other,
        }
    }
    fn curve3d(&self) -> Option<Curve3> {
        Some(self.curve().clone())
    }
}

impl Adaptor3dCurve for Adaptor3dCurveOnSurface {
    fn first_parameter(&self) -> f64 {
        Adaptor3dCurveOnSurface::first_parameter(self)
    }
    fn last_parameter(&self) -> f64 {
        Adaptor3dCurveOnSurface::last_parameter(self)
    }
    fn value(&self, the_u: f64) -> DVec3 {
        Adaptor3dCurveOnSurface::value(self, the_u)
    }
    fn resolution(&self, the_r3d: f64) -> f64 {
        Adaptor3dCurveOnSurface::resolution(self, the_r3d)
    }
    fn is_closed(&self) -> bool {
        // Adaptor3d_CurveOnSurface::IsClosed() — the ends coincide.
        let f = self.first_parameter();
        let l = self.last_parameter();
        self.value(f).distance_squared(self.value(l)) <= CONFUSION * CONFUSION
    }
    fn is_kind_bounded(&self) -> bool {
        // Adaptor3d_CurveOnSurface is a bounded curve (the pcurve range).
        true
    }
    fn get_type(&self) -> AdaptorCurveKind {
        // Adaptor3d_CurveOnSurface::GetType() — the pcurve type (myType).
        let pcurve = self.pcurve_ref();
        match curve2d_type_of(&pcurve) {
            Curve2dType::Line => AdaptorCurveKind::Line,
            Curve2dType::Circle => AdaptorCurveKind::Circle,
            Curve2dType::Ellipse => AdaptorCurveKind::Ellipse,
            Curve2dType::Hyperbola => AdaptorCurveKind::Hyperbola,
            Curve2dType::Parabola => AdaptorCurveKind::Parabola,
            Curve2dType::BSplineCurve => AdaptorCurveKind::BSpline,
            Curve2dType::BezierCurve => AdaptorCurveKind::Bezier,
            Curve2dType::OffsetCurve => AdaptorCurveKind::Offset,
            _ => AdaptorCurveKind::Other,
        }
    }
    fn curve3d(&self) -> Option<Curve3> {
        // GAP (bridge #8): the conic lifting of
        // Adaptor3d_CurveOnSurface::EvalKPart (Adaptor3d_CurveOnSurface.cxx)
        // is untranslated; the OCCT Standard_NoSuchObject path is kept.
        None
    }
    fn curve_on_surface_parts(&self) -> Option<(Curve2d, Surface3)> {
        Some(self.curve_on_surface_values())
    }
}

impl Adaptor3dCurveOnSurface {
    /// The underlying pcurve (bridge: the GetType dispatch reads the pcurve
    /// type; the value accessor models the myCurve handle).
    fn pcurve_ref(&self) -> Curve2d {
        self.curve_on_surface_values().0
    }
}

// ---------------------------------------------------------------------------
// The Extrema tool routing over the local Adaptor3dCurve re-hosts (bridge #4)
// ---------------------------------------------------------------------------

/// The kernel adaptor-stack instance the OCCT `Extrema_ExtPC` /
/// `Extrema_LocateExtPC` ctors receive for the `Adaptor3d_Curve&` argument.
/// Architecture glue: the kernel adaptor owns its curve values, and the
/// declaration order (backing, then the `CurveToolHandle` view, then the
/// Extrema object borrowing it) models the OCCT stack scoping of
/// `GeomAdaptor_Curve` / `Adaptor3d_CurveOnSurface`.
enum AdaptorCurveBacking {
    /// The GeomAdaptor_Curve instance (`curve3d()` route).
    Curve(GeomCurveAdaptor),
    /// The Adaptor3d_CurveOnSurface instance (`curve_on_surface_parts()`
    /// route).
    OnSurface(CurveOnSurface),
}

/// OCCT: the Extrema ctors take the `Adaptor3d_Curve&` itself — the
/// `Extrema_CurveTool` statics dispatch on it. The rcad routing builds the
/// matching kernel adaptor instance: `GeomAdaptor_Curve(C, First, Last)`
/// for the `curve3d()` re-hosts, `Adaptor3d_CurveOnSurface` over
/// `Geom2dAdaptor_Curve(C, First, Last)` + `GeomAdaptor_Surface(S)` for the
/// curve-on-surface re-host. None keeps the bridge #4 GAP exit (the
/// re-host carries no routed parts).
fn adaptor_curve_backing<C: Adaptor3dCurve>(the_curve: &C) -> Option<AdaptorCurveBacking> {
    let u_min = Adaptor3dCurve::first_parameter(the_curve);
    let u_max = Adaptor3dCurve::last_parameter(the_curve);
    if let Some(c3) = the_curve.curve3d() {
        // The GeomAdaptor_Curve instance.
        Some(AdaptorCurveBacking::Curve(GeomCurveAdaptor::with_range(
            c3, u_min, u_max,
        )))
    } else if let Some((c2d, surface)) = the_curve.curve_on_surface_parts() {
        // The Adaptor3d_CurveOnSurface instance.
        Some(AdaptorCurveBacking::OnSurface(CurveOnSurface::new(
            Arc::new(Geom2dCurveAdaptor::with_range(c2d, u_min, u_max)),
            Arc::new(GeomSurfaceAdaptor::new(surface)),
        )))
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Statics of ShapeAnalysis_Curve.cxx
// ---------------------------------------------------------------------------

/// OCCT ProjectOnSegments(theCurve, thePoint, theSegmentCount, aStartParam,
/// anEndParam, aProjDistance, aProjPoint, aProjParam) (cxx L81-124) —
/// projects a point onto a curve by evaluating the curve at several points
/// within the parameter range.
fn project_on_segments<C: Adaptor3dCurve>(
    the_curve: &C,
    the_point: DVec3,
    the_segment_count: i32,
    a_start_param: &mut f64,
    an_end_param: &mut f64,
    a_proj_distance: &mut f64,
    a_proj_point: &mut DVec3,
    a_proj_param: &mut f64,
) {
    // We consider <nbseg> points on [uMin,uMax]
    // Which is the closest? And what is the new interval?
    // (it cannot overflow the old one)

    if the_segment_count <= 0 {
        return; // No segments to project on
    }

    let a_param_step = (*an_end_param - *a_start_param) / the_segment_count as f64;
    let mut a_min_sq_distance = *a_proj_distance * *a_proj_distance;
    let mut a_has_changed = false;
    for i in 0..=the_segment_count {
        let a_current_param = *a_start_param + (a_param_step * i as f64);
        let a_current_point = the_curve.value(a_current_param);
        let a_current_sq_distance = a_current_point.distance_squared(the_point);
        if a_current_sq_distance < a_min_sq_distance {
            a_min_sq_distance = a_current_sq_distance;
            *a_proj_point = a_current_point;
            *a_proj_param = a_current_param;
            a_has_changed = true;
        }
    }
    if a_has_changed {
        *a_proj_distance = a_min_sq_distance.sqrt();
    }

    *an_end_param = (*an_end_param).min(*a_proj_param + a_param_step);
    *a_start_param = (*a_start_param).max(*a_proj_param - a_param_step);
}

/// OCCT ElCLib::AdjustPeriodic(UFirst, ULast, Preci, U1, U2) (the kernel
/// topo/brep_lib/make_edge2d.rs re-host, ElCLib.cxx L115-146) — the local
/// copy for the ValidateRange arm.
fn elclib_adjust_periodic(u_first: f64, u_last: f64, preci: f64, u1: &mut f64, u2: &mut f64) {
    let a_period = u_last - u_first;
    // OCCT L129-136: bring U1 into [UFirst - Preci, ULast + Preci).
    while *u1 < u_first - preci {
        *u1 += a_period;
    }
    while *u1 >= u_last + preci {
        *u1 -= a_period;
    }
    // OCCT L138-144: U2 = U1 + the smallest positive dy giving
    // U1 <= U2 < U1 + Period.
    let mut dy = *u2 - *u1;
    while dy < 0.0 {
        dy += a_period;
    }
    while dy >= a_period {
        dy -= a_period;
    }
    *u2 = *u1 + dy;
}

// ---------------------------------------------------------------------------
// ShapeAnalysis_Curve
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_Curve (hxx) — the curve analysis tool (stateless).
#[derive(Debug, Clone, Default)]
pub struct ShapeAnalysisCurve;

impl ShapeAnalysisCurve {
    /// OCCT Project(C3D, P3D, preci, proj, param, AdjustToEnds)
    /// (cxx L126-140) — the full-range form.
    pub fn project(
        &self,
        c3d: &Curve3,
        p3d: DVec3,
        preci: f64,
        proj: &mut DVec3,
        param: &mut f64,
        adjust_to_ends: bool,
    ) -> f64 {
        let u_min = CurveEval::default_domain(c3d)[0];
        let u_max = CurveEval::default_domain(c3d)[1];
        if u_min < u_max {
            self.project_cf_cl(c3d, p3d, preci, proj, param, u_min, u_max, adjust_to_ends)
        } else {
            self.project_cf_cl(c3d, p3d, preci, proj, param, u_max, u_min, adjust_to_ends)
        }
    }

    /// OCCT Project(C3D, P3D, preci, proj, param, cf, cl, AdjustToEnds)
    /// (cxx L147-203) — the restricted-range form.
    #[allow(clippy::too_many_arguments)]
    pub fn project_cf_cl(
        &self,
        c3d: &Curve3,
        p3d: DVec3,
        preci: f64,
        proj: &mut DVec3,
        param: &mut f64,
        cf: f64,
        cl: f64,
        adjust_to_ends: bool,
    ) -> f64 {
        let mut u_min = if cf < cl { cf } else { cl };
        let mut u_max = if cf < cl { cl } else { cf };

        let mut gac = GeomAdaptorCurve::new(c3d.clone(), u_min, u_max);
        if matches!(
            c3d,
            Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)
        ) {
            // :j8 abv 10 Dec 98: tr10_r0501_db.stp #9423: protection against
            // densing of points near one end
            let prec = if adjust_to_ends { preci } else { CONFUSION };
            let low_bound = gac.value(u_min);
            let hig_bound = gac.value(u_max);
            let mut distmin = low_bound.distance(p3d);
            if distmin <= prec {
                *param = u_min;
                *proj = low_bound;
                return distmin;
            }
            distmin = hig_bound.distance(p3d);
            if distmin <= prec {
                *param = u_max;
                *proj = hig_bound;
                return distmin;
            }
            let _ = distmin;
        }

        if !CurveEval::is_closed(c3d) {
            // modified by rln on 16/12/97 after CSR# PRO11641 entity 20767
            // modified by pdn on 01.07.98 after BUC60195 entity 1952
            // (std::min() added)
            let delta = gac.resolution(preci).min((u_max - u_min) * 0.1);
            u_min -= delta;
            u_max += delta;
            // GAC.Load(C3D, uMin, uMax) — the re-host re-construction
            // (bridge #3).
            gac = GeomAdaptorCurve::new(c3d.clone(), u_min, u_max);
        }

        self.project_act(&gac, p3d, preci, proj, param)
    }

    /// OCCT Project(C3D: Adaptor3d_Curve, P3D, preci, proj, param,
    /// AdjustToEnds) (cxx L205-263).
    pub fn project_adaptor<C: Adaptor3dCurve>(
        &self,
        c3d: &C,
        p3d: DVec3,
        preci: f64,
        proj: &mut DVec3,
        param: &mut f64,
        adjust_to_ends: bool,
    ) -> f64 {
        let u_min = c3d.first_parameter();
        let u_max = c3d.last_parameter();

        if is_infinite_value(u_min) && is_infinite_value(u_max) {
            return self.project_act(c3d, p3d, preci, proj, param);
        }

        let mut distmin_l = INFINITE_VALUE;
        let mut distmin_h = INFINITE_VALUE;
        // :j8 abv 10 Dec 98: tr10_r0501_db.stp #9423
        let prec = if adjust_to_ends { preci } else { CONFUSION };
        let low_bound = c3d.value(u_min);
        let hig_bound = c3d.value(u_max);
        distmin_l = low_bound.distance(p3d);
        distmin_h = hig_bound.distance(p3d);

        if distmin_l <= prec {
            *param = u_min;
            *proj = low_bound;
            return distmin_l;
        }

        if distmin_h <= prec {
            *param = u_max;
            *proj = hig_bound;
            return distmin_h;
        }

        let dist_proj = self.project_act(c3d, p3d, preci, proj, param);
        if dist_proj < distmin_l + CONFUSION && dist_proj < distmin_h + CONFUSION {
            return dist_proj;
        }

        if distmin_l < distmin_h {
            *param = u_min;
            *proj = low_bound;
            return distmin_l;
        }
        *param = u_max;
        *proj = hig_bound;
        distmin_h
    }

    /// OCCT ProjectAct(theCurve, thePoint, theTolerance, theProjPoint,
    /// theProjParam) (cxx L265-502) — the Newton algo core.
    pub fn project_act<C: Adaptor3dCurve>(
        &self,
        the_curve: &C,
        the_point: DVec3,
        the_tolerance: f64,
        the_proj_point: &mut DVec3,
        the_proj_param: &mut f64,
    ) -> f64 {
        let mut ok = false;
        *the_proj_param = 0.0;
        // OCC_CATCH_SIGNALS (bridge #8): the OCCT try/catch -> the OK=false
        // exit is preserved at the GAP site below.
        // Extrema_ExtPC aCurveExtrema(thePoint, theCurve) (cxx L277) — the
        // two-arg ctor over the adaptor, the default theTolF is 1.0e-10.
        let a_backing = adaptor_curve_backing(the_curve);
        let a_ext_pc_tool = match a_backing.as_ref() {
            Some(AdaptorCurveBacking::Curve(a_gac)) => {
                Some(CurveToolHandle::for_curve3(&a_gac.curve, a_gac, a_gac))
            }
            Some(AdaptorCurveBacking::OnSurface(a_cos)) => {
                Some(CurveToolHandle::with_geom(a_cos, a_cos))
            }
            None => {
                // GAP (bridge #4): the re-host carries no routed parts (the
                // CurveOnPlane stand-in); the OCCT `!IsDone() -> OK=false`
                // exit is preserved and the ProjectOnSegments /
                // Extrema_LocateExtPC fallback below stays reachable.
                None
            }
        };
        let a_curve_extrema = a_ext_pc_tool
            .as_ref()
            .map(|a_tool| ExtremaExtPC::new_point_curve(the_point, a_tool, 1.0e-10));
        let mut a_min_extrema_distance = REAL_LAST;
        let mut a_min_extrema_index = 0usize;
        if let Some(a_curve_extrema) = a_curve_extrema.as_ref() {
            if a_curve_extrema.is_done() && a_curve_extrema.nb_ext() > 0 {
                for i in 1..=a_curve_extrema.nb_ext() {
                    // OCCT IsMin(i): the kernel ExtPC stores minima only
                    // (bridge #4).
                    let a_current_distance = a_curve_extrema.square_distance(i);
                    if a_current_distance < a_min_extrema_distance {
                        a_min_extrema_distance = a_current_distance;
                        a_min_extrema_index = i;
                    }
                }

                if a_min_extrema_index != 0 {
                    *the_proj_param = a_curve_extrema.point(a_min_extrema_index).param;
                    *the_proj_point = a_curve_extrema.point(a_min_extrema_index).point;
                    ok = true;
                }
            }
        }

        // szv#4:S4163:12Mar99 moved
        let mut u_min = the_curve.first_parameter();
        let mut u_max = the_curve.last_parameter();
        let mut an_is_closed_curve = false;
        let mut a_curve_period = 0.0f64;
        // Distance between the point and the projection point.
        let mut a_proj_distance = INFINITE_VALUE;
        let mut a_mod_min = INFINITE_VALUE;

        // Remember the computed values.
        // These values will be used in case the projection is not successful.
        let a_computed_param = *the_proj_param;
        let a_computed_proj = *the_proj_point;

        // PTV 29.05.2002 remember the old solution, cause it could be better
        let mut an_is_have_old_solution = false;
        let mut an_old_param = 0.0f64;
        let mut an_old_proj = DVec3::ZERO;
        if ok {
            an_is_have_old_solution = true;
            an_old_proj = *the_proj_point;
            an_old_param = *the_proj_param;
            a_proj_distance = the_proj_point.distance(the_point);
            a_mod_min = a_proj_distance;
            if a_proj_distance > the_tolerance {
                ok = false;
            }
            if Adaptor3dCurve::is_closed(the_curve) {
                an_is_closed_curve = true;
                a_curve_period = u_max - u_min; // szv#4:S4163:12Mar99 optimized
            }
        }

        if !ok {
            // BUG NICOLAS - If the point is on the curve 0 Solutions.
            // This works with ElCLib

            // Generally speaking, we try to ALWAYS return a result that's
            // NOT EVEN GOOD. The caller can then decide what to do.
            *the_proj_param = 0.0;

            match the_curve.get_type() {
                AdaptorCurveKind::Circle => {
                    let Some(curve3) = the_curve.curve3d() else {
                        panic!(
                            "GAP: Adaptor3d_CurveOnSurface::Circle() (the EvalKPart lifting, \
                             Adaptor3d_CurveOnSurface.cxx) is not translated — \
                             see curve.rs header bridge #8"
                        );
                    };
                    let Curve3::Circle(a_circ) = curve3 else {
                        panic!("GAP: the adaptor kind/curve mismatch (bridge #8)");
                    };
                    *the_proj_point = a_circ.center;
                    if a_circ.radius <= REAL_SMALL
                        || the_point.distance_squared(*the_proj_point) <= REAL_SMALL
                    {
                        *the_proj_param = the_curve.first_parameter();
                        *the_proj_point = *the_proj_point
                            + a_circ.x_dir * a_circ.radius;
                    } else {
                        *the_proj_param =
                            crate::geomalgo::int_patch::elclib::circle_parameter(&a_circ, the_point);
                        *the_proj_point =
                            crate::geomalgo::int_patch::elclib::circle_value(&a_circ, *the_proj_param);
                    }
                    an_is_closed_curve = true;
                    a_curve_period = 2.0 * std::f64::consts::PI;
                }
                AdaptorCurveKind::Hyperbola => {
                    let Some(curve3) = the_curve.curve3d() else {
                        panic!(
                            "GAP: Adaptor3d_CurveOnSurface::Hyperbola() (the EvalKPart lifting, \
                             Adaptor3d_CurveOnSurface.cxx) is not translated — \
                             see curve.rs header bridge #8"
                        );
                    };
                    let Curve3::Hyperbola(a_hypo) = curve3 else {
                        panic!("GAP: the adaptor kind/curve mismatch (bridge #8)");
                    };
                    *the_proj_param =
                        crate::geomalgo::int_patch::elclib::hyperbola_parameter(&a_hypo, the_point);
                    *the_proj_point =
                        crate::geomalgo::int_patch::elclib::hyperbola_value(&a_hypo, *the_proj_param);
                }
                AdaptorCurveKind::Parabola => {
                    let Some(curve3) = the_curve.curve3d() else {
                        panic!(
                            "GAP: Adaptor3d_CurveOnSurface::Parabola() (the EvalKPart lifting, \
                             Adaptor3d_CurveOnSurface.cxx) is not translated — \
                             see curve.rs header bridge #8"
                        );
                    };
                    let Curve3::Parabola(a_para) = curve3 else {
                        panic!("GAP: the adaptor kind/curve mismatch (bridge #8)");
                    };
                    *the_proj_param =
                        crate::geomalgo::int_patch::elclib::parabola_parameter(&a_para, the_point);
                    *the_proj_point =
                        crate::geomalgo::int_patch::elclib::parabola_value(&a_para, *the_proj_param);
                }
                AdaptorCurveKind::Line => {
                    let Some(curve3) = the_curve.curve3d() else {
                        panic!(
                            "GAP: Adaptor3d_CurveOnSurface::Line() (the EvalKPart lifting, \
                             Adaptor3d_CurveOnSurface.cxx) is not translated — \
                             see curve.rs header bridge #8"
                        );
                    };
                    let Curve3::Line(a_lin) = curve3 else {
                        panic!("GAP: the adaptor kind/curve mismatch (bridge #8)");
                    };
                    *the_proj_param =
                        crate::geomalgo::int_patch::elclib::line_parameter(&a_lin, the_point);
                    *the_proj_point =
                        crate::geomalgo::int_patch::elclib::line_value(&a_lin, *the_proj_param);
                }
                AdaptorCurveKind::Ellipse => {
                    let Some(curve3) = the_curve.curve3d() else {
                        panic!(
                            "GAP: Adaptor3d_CurveOnSurface::Ellipse() (the EvalKPart lifting, \
                             Adaptor3d_CurveOnSurface.cxx) is not translated — \
                             see curve.rs header bridge #8"
                        );
                    };
                    let Curve3::Ellipse(a_ell) = curve3 else {
                        panic!("GAP: the adaptor kind/curve mismatch (bridge #8)");
                    };
                    *the_proj_param =
                        crate::geomalgo::int_patch::elclib::ellipse_parameter(&a_ell, the_point);
                    *the_proj_point =
                        crate::geomalgo::int_patch::elclib::ellipse_value(&a_ell, *the_proj_param);
                    an_is_closed_curve = true;
                    a_curve_period = 2.0 * std::f64::consts::PI;
                }
                _ => {
                    // Bspline or something, no easy solution.
                    // Here this algorithm tries to find the closest point on
                    // the curve by evaluating the distance between the point
                    // and the curve at several points within the parameter
                    // range.
                    a_proj_distance = INFINITE_VALUE;
                    project_on_segments(
                        the_curve,
                        the_point,
                        25,
                        &mut u_min,
                        &mut u_max,
                        &mut a_proj_distance,
                        the_proj_point,
                        the_proj_param,
                    );
                    if a_proj_distance <= the_tolerance {
                        return a_proj_distance;
                    }

                    // OCCT L421-426: Extrema_LocateExtPC aProjector(thePoint,
                    // theCurve, theProjParam /*U0*/, uMin, uMax,
                    // theTolerance /*TolU*/).
                    if let Some(a_backing) = adaptor_curve_backing(the_curve) {
                        let a_projector_tool = match &a_backing {
                            AdaptorCurveBacking::Curve(a_gac) => {
                                CurveToolHandle::for_curve3(&a_gac.curve, a_gac, a_gac)
                            }
                            AdaptorCurveBacking::OnSurface(a_cos) => {
                                CurveToolHandle::with_geom(a_cos, a_cos)
                            }
                        };
                        let a_projector = LocateExtPC::new_point_curve_seed_ranged(
                            the_point,
                            &a_projector_tool,
                            *the_proj_param, // U0
                            u_min,
                            u_max,
                            the_tolerance, // TolU
                        );
                        // OCCT L427-434.
                        if a_projector.is_done() {
                            *the_proj_param = a_projector.point().param;
                            *the_proj_point = a_projector.point().point;
                            let a_dist_newton = the_point.distance(*the_proj_point);
                            if a_dist_newton < a_mod_min {
                                return a_dist_newton;
                            }
                        }
                    }

                    // Here we are trying to find the closest point on the
                    // curve by repeatedly evaluating the distance between
                    // the point and the curve at several points within the
                    // parameter range.  The particular number of segments to
                    // probe was chosen to be 40, 20, 25, and 40 again.
                    for a_segment_count in [40, 20, 25, 40] {
                        project_on_segments(
                            the_curve,
                            the_point,
                            a_segment_count,
                            &mut u_min,
                            &mut u_max,
                            &mut a_proj_distance,
                            the_proj_point,
                            the_proj_param,
                        );
                        if a_proj_distance <= the_tolerance {
                            return a_proj_distance;
                        }
                    }

                    // Did not find a point on the curve that is closer than
                    // the tolerance. So, we return the closest point found
                    // so far.
                    if a_proj_distance > a_mod_min {
                        a_proj_distance = a_mod_min;
                        *the_proj_param = a_computed_param;
                        *the_proj_point = a_computed_proj;
                    }

                    return a_proj_distance;
                }
            }
        }

        if an_is_closed_curve && (*the_proj_param < u_min || *the_proj_param > u_max) {
            *the_proj_param += shape_analysis_adjust_by_period(
                *the_proj_param,
                0.5 * (u_min + u_max),
                a_curve_period,
            );
        }

        if an_is_have_old_solution {
            // PTV 29.05.2002 Compare old solution and new;
            let an_old_dist = an_old_proj.distance_squared(the_point);
            let a_new_dist = the_proj_point.distance_squared(the_point);
            if an_old_dist < a_new_dist {
                *the_proj_point = an_old_proj;
                *the_proj_param = an_old_param;
            }
        }
        the_proj_point.distance(the_point)
    }

    /// OCCT NextProject(paramPrev, C3D, P3D, preci, proj, param, cf, cl,
    /// AdjustToEnds) (cxx L504-560) — Newton algo for projecting point on
    /// curve (S4030).
    #[allow(clippy::too_many_arguments)]
    pub fn next_project(
        &self,
        param_prev: f64,
        c3d: &Curve3,
        p3d: DVec3,
        preci: f64,
        proj: &mut DVec3,
        param: &mut f64,
        cf: f64,
        cl: f64,
        adjust_to_ends: bool,
    ) -> f64 {
        let mut u_min = if cf < cl { cf } else { cl };
        let mut u_max = if cf < cl { cl } else { cf };
        let mut distmin = INFINITE_VALUE;
        let mut gac = GeomAdaptorCurve::new(c3d.clone(), u_min, u_max);
        if matches!(
            c3d,
            Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)
        ) {
            // :j8 abv 10 Dec 98: tr10_r0501_db.stp #9423
            let prec = if adjust_to_ends { preci } else { CONFUSION };
            let low_bound = gac.value(u_min);
            let hig_bound = gac.value(u_max);
            distmin = low_bound.distance(p3d);
            if distmin <= prec {
                *param = u_min;
                *proj = low_bound;
                return distmin;
            }
            distmin = hig_bound.distance(p3d);
            if distmin <= prec {
                *param = u_max;
                *proj = hig_bound;
                return distmin;
            }
        }

        if !CurveEval::is_closed(c3d) {
            // modified by rln on 16/12/97; modified by pdn on 01.07.98
            let delta = gac.resolution(preci).min((u_max - u_min) * 0.1);
            u_min -= delta;
            u_max += delta;
            // GAC.Load(C3D, uMin, uMax) (bridge #3).
            gac = GeomAdaptorCurve::new(c3d.clone(), u_min, u_max);
        }
        self.next_project_adaptor(param_prev, &gac, p3d, preci, proj, param)
    }

    /// OCCT NextProject(paramPrev, C3D: Adaptor3d_Curve, P3D, preci, proj,
    /// param) (cxx L561-580).
    pub fn next_project_adaptor<C: Adaptor3dCurve>(
        &self,
        param_prev: f64,
        c3d: &C,
        p3d: DVec3,
        preci: f64,
        proj: &mut DVec3,
        param: &mut f64,
    ) -> f64 {
        let u_min = c3d.first_parameter();
        let u_max = c3d.last_parameter();

        // OCCT L571: Extrema_LocateExtPC aProjector(P3D, C3D, paramPrev
        // /*U0*/, uMin, uMax, preci /*TolU*/).
        if let Some(a_backing) = adaptor_curve_backing(c3d) {
            let a_projector_tool = match &a_backing {
                AdaptorCurveBacking::Curve(a_gac) => {
                    CurveToolHandle::for_curve3(&a_gac.curve, a_gac, a_gac)
                }
                AdaptorCurveBacking::OnSurface(a_cos) => CurveToolHandle::with_geom(a_cos, a_cos),
            };
            let a_projector = LocateExtPC::new_point_curve_seed_ranged(
                p3d,
                &a_projector_tool,
                param_prev, // U0
                u_min,
                u_max,
                preci, // TolU
            );
            // OCCT L572-576.
            if a_projector.is_done() {
                *param = a_projector.point().param;
                *proj = a_projector.point().point;
                return p3d.distance(*proj);
            }
        }
        self.project_adaptor(c3d, p3d, preci, proj, param, false)
    }

    /// OCCT ValidateRange(theCurve, First, Last, preci) (cxx L586-739) —
    /// copied from StepToTopoDS_GeometricTuul::UpdateParam3d (Aug 2001).
    pub fn validate_range(&self, the_curve: &Curve3, first: &mut f64, last: &mut f64, preci: f64) -> bool {
        // First et/ou Last peuvent etre en dehors des bornes naturelles de
        // la courbe. On donnera alors la valeur en bout a First et/ou Last

        let cf = CurveEval::default_domain(the_curve)[0];
        let cl = CurveEval::default_domain(the_curve)[1];
        let bounded = matches!(
            the_curve,
            Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)
        );
        let is_closed = CurveEval::is_closed(the_curve);

        if bounded && !is_closed {
            if *first < cf {
                *first = cf;
            } else if *first > cl {
                *first = cl;
            }
            if *last < cf {
                *last = cf;
            } else if *last > cl {
                *last = cl;
            }
        }

        // 15.11.2002 PTV OCC966
        if self.is_periodic(the_curve) {
            // :a7 abv 11 Feb 98: preci -> PConfusion()
            elclib_adjust_periodic(cf, cl, PCONFUSION, first, last);
        } else if *first < *last {
            // nothing to fix
        } else if is_closed {
            // l'un des points projecte se trouve sur l'origine du
            // parametrage de la courbe 3D.

            // Last = cf au lieu de Last = cl
            if (*last - cf).abs() < PCONFUSION {
                *last = cl;
                // First = cl au lieu de First = cf
            } else if (*first - cl).abs() < PCONFUSION {
                *first = cf;
                // on se trouve dans un cas ou l origine est traversee
                // illegal sur une courbe fermee non periodique
                // on inverse quand meme les parametres !!!!!!
            } else {
                // : S4136 abv 20 Apr 99: r0701_ug.stp #6230: add check in 3d
                if CurveEval::point_at(the_curve, *first)
                    .distance(CurveEval::point_at(the_curve, cf))
                    < preci
                {
                    *first = cf;
                }
                if CurveEval::point_at(the_curve, *last)
                    .distance(CurveEval::point_at(the_curve, cl))
                    < preci
                {
                    *last = cl;
                }
                if *first > *last {
                    let tmp = *first;
                    *first = *last;
                    *last = tmp;
                }
            }
        }
        // The curve is closed within the 3D tolerance
        else if let Curve3::BSpline(a_bspline) = the_curve {
            let start = CurveEval::point_at(the_curve, cf);
            let end = CurveEval::point_at(the_curve, cl);
            let _ = &mut (start, end);
            let start_point = a_bspline
                .control_points
                .first()
                .copied()
                .unwrap_or(DVec3::ZERO);
            let end_point = a_bspline
                .control_points
                .last()
                .copied()
                .unwrap_or(DVec3::ZERO);
            if start_point.distance(end_point) <= preci {
                // : S4136 <= BRepAPI::Precision()
                // Last = cf au lieu de Last = cl
                if (*last - cf).abs() < PCONFUSION {
                    *last = cl;
                    // First = cl au lieu de First = cf
                } else if (*first - cl).abs() < PCONFUSION {
                    *first = cf;
                } else {
                    let tmp = *first;
                    *first = *last;
                    *last = tmp;
                }
            }
            // abv 15.03.00 #72 bm1_pe_t4 protection of exceptions in draw
            else if *first > *last {
                *first = CurveEval::reversed_parameter(the_curve, *first);
                *last = CurveEval::reversed_parameter(the_curve, *last);
                // theCurve->Reverse() — the rcad curve value is immutable;
                // the caller re-anchors the reversed curve (architecture
                // bridge: the Geom value Reverse is a no-op here, the
                // parameter swap carries the OCCT effect on [First, Last]).
            }
            // : j9 abv 11 Dec 98: PRO7747 #4875, after :j8: else
            if *first == *last {
                // gka 10.07.1998 file PRO7656 entity 33334
                *first = cf;
                *last = cl;
                return false;
            }
        } else {
            // abv 15.03.00 #72 bm1_pe_t4 protection of exceptions in draw
            if *first > *last {
                *first = CurveEval::reversed_parameter(the_curve, *first);
                *last = CurveEval::reversed_parameter(the_curve, *last);
                // theCurve->Reverse() — see the BSpline arm note.
            }
            // pdn 11.01.99 #144 bm1_pe_t4 protection of exceptions in draw
            if *first == *last {
                *first -= PCONFUSION;
                *last += PCONFUSION;
            }
            return false;
        }
        true
    }

    /// OCCT FillBndBox(C2d, First, Last, NPoints, Exact, Box) (cxx L788-850)
    /// — WORK-AROUND for methods from BndLib which give not exact bounds.
    pub fn fill_bnd_box(
        &self,
        c2d: &Curve2d,
        first: f64,
        last: f64,
        n_points: i32,
        exact: bool,
        box2d: &mut BndBox2d,
    ) {
        if !exact {
            let nseg = if n_points < 2 { 1 } else { n_points - 1 };
            let step = (last - first) / nseg as f64;
            for i in 0..=nseg {
                let par = first + i as f64 * step;
                let pnt = Curve2dEval::point_at(c2d, par);
                box2d.add_point(pnt);
            }
            return;
        }

        // We should solve the task on intervals of C2 continuity.
        // Geom2dAdaptor_Curve anAC(C2d, First, Last) — the restricted-range
        // C2 interval count (bridge #6).
        let nb_int = geom2d_adaptor_nb_intervals_c2(c2d, first, last);
        // If we have only 1 interval then use input NPoints parameter to get
        // samples.
        let nb_samples = if nb_int < 2 { n_points - 1 } else { nb_int };
        let mut a_params = vec![0.0f64; (nb_samples + 1) as usize];
        if nb_samples == nb_int {
            // anAC.Intervals(aParams, GeomAbs_C2): the C2 boundaries in
            // [First, Last] (1-based: aParams(1) = First, aParams(N+1) =
            // Last).
            fill_c2_intervals(c2d, first, last, &mut a_params);
        } else {
            let step = (last - first) / nb_samples as f64;
            for i in 0..=nb_samples {
                a_params[i as usize] = first + i as f64 * step;
            }
        }
        for i in 1..=nb_samples + 1 {
            let a_par1 = a_params[(i - 1) as usize];
            let a_pnt = Curve2dEval::point_at(c2d, a_par1);
            box2d.add_point(a_pnt);
            if i <= nb_samples {
                let a_par2 = a_params[i as usize];
                let par = (a_par1 + a_par2) * 0.5;
                let mut pextr = DVec2::ZERO;
                let mut parextr = par;
                if search_for_extremum(c2d, a_par1, a_par2, DVec2::new(1.0, 0.0), &mut parextr, &mut pextr)
                {
                    box2d.add_point(pextr);
                }
                parextr = par;
                if search_for_extremum(c2d, a_par1, a_par2, DVec2::new(0.0, 1.0), &mut parextr, &mut pextr)
                {
                    box2d.add_point(pextr);
                }
            }
        }
    }

    /// OCCT SelectForwardSeam(C1, C2) (cxx L852-975).
    pub fn select_forward_seam(&self, c1: &Curve2d, c2: &Curve2d) -> i32 {
        //  SelectForward est destine a devenir un outil distinct
        //  Il est sans doute optimisable !

        let mut the_curve_indice = 0i32;

        // occ::handle<Geom2d_Line> L1 = down_cast<Geom2d_Line>(C1) — the
        // Line2d value; the BoundedCurve fallback builds a line from the
        // end points.
        let l1: Option<(DVec2, DVec2)> = match c1 {
            Curve2d::Line(l) => Some((l.origin, l.direction)),
            Curve2d::BSpline(_) | Curve2d::Bezier(_) | Curve2d::Trimmed(_) => {
                // gp_Vec2d VecBC1(StartBC1, EndBC1).
                let start_bc1 = Curve2dEval::point_at(c1, curve2d_first_param(c1));
                let end_bc1 = Curve2dEval::point_at(c1, curve2d_last_param(c1));
                let vec_bc1 = end_bc1 - start_bc1;
                if vec_bc1.length_squared() < REAL_SMALL {
                    return the_curve_indice;
                }
                Some((start_bc1, vec_bc1))
            }
            _ => None,
        };
        let Some((l1_loc, l1_dir)) = l1 else {
            return the_curve_indice;
        };

        let l2: Option<(DVec2, DVec2)> = match c2 {
            Curve2d::Line(l) => Some((l.origin, l.direction)),
            Curve2d::BSpline(_) | Curve2d::Bezier(_) | Curve2d::Trimmed(_) => {
                let start_bc2 = Curve2dEval::point_at(c2, curve2d_first_param(c2));
                let end_bc2 = Curve2dEval::point_at(c2, curve2d_last_param(c2));
                let vec_bc2 = end_bc2 - start_bc2;
                if vec_bc2.length_squared() < REAL_SMALL {
                    return the_curve_indice;
                }
                Some((start_bc2, vec_bc2))
            }
            _ => None,
        };
        let Some((l2_loc, l2_dir)) = l2 else {
            return the_curve_indice;
        };

        let mut udir_pos = false;
        let mut udir_neg = false;
        let mut vdir_pos = false;
        let mut vdir_neg = false;

        // gp_Dir2d theDir = L1->Direction().
        let the_dir = l1_dir.normalize_or_zero();
        // gp_Pnt2d theLoc1 = L1->Location(); theLoc2 = L2->Location().
        let the_loc1 = l1_loc;
        let the_loc2 = l2_loc;

        if the_dir.x > 0.0 {
            udir_pos = true; // szv#4:S4163:12Mar99 Udir unused
        } else if the_dir.x < 0.0 {
            udir_neg = true;
        } else if the_dir.y > 0.0 {
            vdir_pos = true;
        } else if the_dir.y < 0.0 {
            vdir_neg = true;
        }

        if vdir_pos {
            // max of Loc1.X() Loc2.X()
            if the_loc1.x > the_loc2.x {
                the_curve_indice = 1;
            } else {
                the_curve_indice = 2;
            }
        } else if vdir_neg {
            if the_loc1.x > the_loc2.x {
                the_curve_indice = 2;
            } else {
                the_curve_indice = 1;
            }
        } else if udir_pos {
            // min of Loc1.X() Loc2.X()
            if the_loc1.y < the_loc2.y {
                the_curve_indice = 1;
            } else {
                the_curve_indice = 2;
            }
        } else if udir_neg {
            if the_loc1.y < the_loc2.y {
                the_curve_indice = 2;
            } else {
                the_curve_indice = 1;
            }
        }

        the_curve_indice
    }

    /// OCCT IsPlanar(pnts, Normal, preci) (cxx L1097-1172) — detects if
    /// points lie in some plane and returns normal.
    pub fn is_planar_pnts(pnts: &[DVec3], normal: &mut DVec3, preci: f64) -> bool {
        let precision = if preci > 0.0 { preci } else { CONFUSION };
        let no_norm = normal.length_squared() == 0.0;

        if pnts.len() < 3 {
            let n1 = pnts[0] - pnts[1];
            if no_norm {
                *normal = get_any_normal(n1);
                return true;
            }
            return n1.dot(*normal).abs() < CONFUSION;
        }

        let mut a_max_dir = DVec3::ZERO;
        if no_norm {
            // define a center point
            let mut a_center = DVec3::ZERO;
            for p in pnts {
                a_center += *p;
            }
            a_center /= pnts.len() as f64;

            a_max_dir = pnts[0] - a_center;
            *normal = (pnts[pnts.len() - 1] - a_center).cross(a_max_dir);

            for i in 0..pnts.len() - 1 {
                let a_tmp_dir = pnts[i + 1] - a_center;
                if a_tmp_dir.length_squared() > a_max_dir.length_squared() {
                    a_max_dir = a_tmp_dir;
                }

                let mut a_delta = (pnts[i] - a_center).cross(pnts[i + 1] - a_center);
                if normal.dot(a_delta) < 0.0 {
                    a_delta *= -1.0;
                }
                *normal += a_delta;
            }
        }

        // check if points are linear
        let nrm = normal.length();
        if nrm < CONFUSION {
            *normal = get_any_normal(a_max_dir);
            return true;
        }
        *normal /= nrm;

        let mut mind = REAL_LAST;
        let mut maxd = -REAL_LAST;
        for p in pnts {
            let dev = p.dot(*normal);
            if dev < mind {
                mind = dev;
            }
            if dev > maxd {
                maxd = dev;
            }
        }

        (maxd - mind) <= precision
    }

    /// OCCT IsPlanar(curve, Normal, preci) (cxx L1175-1256).
    pub fn is_planar(curve: &Curve3, normal: &mut DVec3, preci: f64) -> bool {
        let precision = if preci > 0.0 { preci } else { CONFUSION };
        let no_norm = normal.length_squared() == 0.0;

        match curve {
            Curve3::Line(line) => {
                // gp_XYZ N1 = Line->Position().Direction().XYZ().
                let n1 = line.direction;
                if no_norm {
                    *normal = get_any_normal(n1);
                    return true;
                }
                return n1.dot(*normal).abs() < CONFUSION;
            }
            Curve3::Circle(_)
            | Curve3::Ellipse(_)
            | Curve3::Hyperbola(_)
            | Curve3::Parabola(_) => {
                // Geom_Conic: gp_XYZ N1 = Conic->Axis().Direction().XYZ().
                let n1 = match curve {
                    Curve3::Circle(c) => c.normal,
                    Curve3::Ellipse(e) => e.normal,
                    Curve3::Hyperbola(h) => h.normal,
                    Curve3::Parabola(p) => p.normal,
                    _ => unreachable!(),
                };
                if no_norm {
                    *normal = n1;
                    return true;
                }
                let a_vec_mul = n1.cross(*normal);
                return a_vec_mul.length_squared() < CONFUSION * CONFUSION;
            }
            Curve3::Trimmed(trimmed) => {
                return Self::is_planar(&trimmed.curve, normal, precision);
            }
            Curve3::Offset(offset) => {
                // Geom_OffsetCurve::BasisCurve().
                return Self::is_planar(&offset.basis, normal, precision);
            }
            Curve3::BSpline(bspline) => {
                return Self::is_planar_pnts(&bspline.control_points, normal, precision);
            }
            Curve3::Bezier(bezier) => {
                return Self::is_planar_pnts(&bezier.control_points, normal, precision);
            }
            _ => {
                // GAP (bridge #7): ShapeExtend_ComplexCurve is not
                // translated yet; the OCCT failure path (`return false`)
                // is preserved.
                return false;
            }
        }
    }

    /// OCCT GetSamplePoints(curve, first, last, seq) (cxx L1258-1315) — the
    /// 3D form.
    pub fn get_sample_points(
        curve: &Curve3,
        first: f64,
        last: f64,
        seq: &mut Vec<DVec3>,
    ) -> bool {
        let c_first = CurveEval::default_domain(curve)[0];
        let c_last = CurveEval::default_domain(curve)[1];
        let adelta = c_last - c_first;
        if adelta == 0.0 {
            return false;
        }

        let a_k = ((last - first) / adelta).ceil() as i32;
        let mut nbp = 100 * a_k;
        match curve {
            Curve3::Line(_) => {
                nbp = 2;
            }
            Curve3::Circle(_) => {
                nbp = 360 * a_k;
            }
            Curve3::BSpline(a_bspl) => {
                // aBspl->NbKnots(): the distinct-knot count — the flat
                // vector carries the multiplicities; deduplicate.
                let mut nb_knots = 0usize;
                let mut prev = f64::NAN;
                for k in &a_bspl.knots {
                    if *k != prev {
                        nb_knots += 1;
                        prev = *k;
                    }
                }
                nbp = (nb_knots as i32) * (a_bspl.degree as i32) * a_k;
                if (nbp as f64) < 2.0 {
                    nbp = 2;
                }
            }
            Curve3::Bezier(a_b) => {
                nbp = 3 + a_b.control_points.len() as i32;
            }
            Curve3::Offset(a_c) => {
                return Self::get_sample_points(&a_c.basis, first, last, seq);
            }
            Curve3::Trimmed(a_c) => {
                return Self::get_sample_points(&a_c.curve, first, last, seq);
            }
            _ => {}
        }

        let step = (last - first) / (nbp - 1) as f64;
        for i in 0..nbp - 1 {
            seq.push(CurveEval::point_at(curve, first + step * i as f64));
        }
        seq.push(CurveEval::point_at(curve, last));
        true
    }

    /// OCCT GetSamplePoints(curve, first, last, seq) (cxx L1317-1420) — the
    /// 2D form.
    pub fn get_sample_points_2d(
        curve: &Curve2d,
        first: f64,
        last: f64,
        seq: &mut Vec<DVec2>,
    ) -> bool {
        //: abv 05.06.02: TUBE.stp
        // Use the same distribution of points as BRepTopAdaptor_FClass2d for
        // consistency
        // Geom2dAdaptor_Curve C(curve, first, last).
        let mut nbs = Curve2dAdaptor::nb_samples(curve);
        //-- Attention aux bsplines rationnelles de degree 3.
        if nbs > 2 {
            nbs *= 4;
        }
        let step = (last - first) / (nbs - 1) as f64;
        for i in 0..nbs - 1 {
            seq.push(Curve2dEval::point_at(curve, first + step * i as f64));
        }
        seq.push(Curve2dEval::point_at(curve, last));
        true
    }

    /// OCCT IsClosed(theCurve, preci) (cxx L1424-1447).
    pub fn is_closed(the_curve: &Curve3, preci: f64) -> bool {
        if CurveEval::is_closed(the_curve) {
            return true;
        }

        let prec = preci.max(CONFUSION);

        let f = CurveEval::default_domain(the_curve)[0];
        let l = CurveEval::default_domain(the_curve)[1];

        if is_infinite_value(f) || is_infinite_value(l) {
            return false;
        }

        let a_closed_val =
            CurveEval::point_at(the_curve, f).distance_squared(CurveEval::point_at(the_curve, l));
        let preci2 = prec * prec;

        a_closed_val <= preci2
    }

    /// OCCT IsPeriodic(theCurve) (cxx L1450-1469) — ask IsPeriodic on
    /// BasisCurve.
    pub fn is_periodic(&self, the_curve: &Curve3) -> bool {
        let mut a_tmp_curve = the_curve.clone();
        loop {
            match &a_tmp_curve {
                Curve3::Offset(off) => a_tmp_curve = (*off.basis).clone(),
                Curve3::Trimmed(trm) => a_tmp_curve = (*trm.curve).clone(),
                _ => break,
            }
        }
        CurveEval::is_periodic(&a_tmp_curve)
    }

    /// OCCT IsPeriodic(theCurve: Geom2d_Curve) (cxx L1471-1490).
    pub fn is_periodic_2d(&self, the_curve: &Curve2d) -> bool {
        let mut a_tmp_curve = the_curve.clone();
        loop {
            match &a_tmp_curve {
                Curve2d::Offset(off) => a_tmp_curve = (*off.basis).clone(),
                Curve2d::Trimmed(trm) => a_tmp_curve = (*trm.curve).clone(),
                _ => break,
            }
        }
        Curve2dEval::is_periodic(&a_tmp_curve)
    }
}

/// OCCT ShapeAnalysis::AdjustByPeriod(Val, ToVal, Period) (ShapeAnalysis.cxx
/// L48-62) — the local copy used by ProjectAct (the statics module carries
/// the 1:1 translation).
pub(crate) fn shape_analysis_adjust_by_period(the_val: f64, to_val: f64, period: f64) -> f64 {
    let diff = the_val - to_val;
    let d = diff.abs();
    let p = period.abs();
    if d <= 0.5 * p {
        return 0.0;
    }
    if p < 1e-100 {
        return diff;
    }
    (if diff > 0.0 { -p } else { p }) * (d / p + 0.5).floor()
}

/// OCCT SearchForExtremum(C2d, First, Last, dir, par, res) (cxx L741-786) —
/// search for extremum using Newton.
fn search_for_extremum(
    c2d: &Curve2d,
    first: f64,
    last: f64,
    dir: DVec2,
    par: &mut f64,
    res: &mut DVec2,
) -> bool {
    let mut nb_out = 0i32;
    for _ in 0..10 {
        let prevpar = *par;

        // C2d->D2(par, res, D1, D2) — the Curve2dAdaptor D2 re-host.
        let (r, d1, d2) = Curve2dAdaptor::d2(c2d, *par);
        *res = r;
        let det = d2.dot(dir);
        if det.abs() < 1e-10 {
            return true;
        }

        *par -= d1.dot(dir) / det;
        if (*par - prevpar).abs() < PCONFUSION {
            return true;
        }

        if *par < first {
            nb_out += 1;
            if nb_out > 3 || prevpar == first {
                return false;
            }
            *par = first;
        }
        if *par > last {
            nb_out += 1;
            if nb_out > 3 || prevpar == last {
                return false;
            }
            *par = last;
        }
    }
    true
}

/// OCCT GetAnyNormal(orig) (cxx L978-1000).
fn get_any_normal(orig: DVec3) -> DVec3 {
    if orig.z.abs() < CONFUSION {
        DVec3::new(0.0, 0.0, 1.0)
    } else {
        let mut norm = DVec3::new(orig.z, 0.0, -orig.x);
        let nrm = norm.length();
        if nrm < CONFUSION {
            norm = DVec3::new(0.0, 0.0, 1.0);
        } else {
            norm /= nrm;
        }
        norm
    }
}

/// OCCT Geom2dAdaptor_Curve::NbIntervals(GeomAbs_C2) (Geom2dAdaptor_Curve.cxx
/// L409-450) — the restricted-range C2 count over the pcurve (bridge #6).
/// BSpline: aCont = 2, the flat-knot multiplicity walk (the
/// brep_fill_sweep_c.rs re-host precedent); other types: 1.
fn geom2d_adaptor_nb_intervals_c2(c2d: &Curve2d, first: f64, last: f64) -> i32 {
    if let Curve2d::BSpline(bs) = c2d {
        let a_cont = 2i32;
        let mut intervals = 0i32;
        let mut i = 0usize;
        let n = bs.knots.len();
        while i < n {
            let mut j = i;
            while j < n && (bs.knots[j] - bs.knots[i]).abs() < 1.0e-12 {
                j += 1;
            }
            let mult = (j - i) as i32;
            let k = bs.knots[i];
            if k > first + 1.0e-12 && k < last - 1.0e-12 && mult > a_cont {
                intervals += 1;
            }
            i = j;
        }
        return intervals + 1;
    }
    1
}

/// OCCT Geom2dAdaptor_Curve::Intervals(T, GeomAbs_C2) — fills the C2
/// boundaries in [First, Last]: T(1) = First, the qualifying interior knots,
/// T(N+1) = Last (the same flat-knot walk).
fn fill_c2_intervals(c2d: &Curve2d, first: f64, last: f64, t: &mut [f64]) {
    let n = t.len() - 1;
    t[0] = first;
    t[n] = last;
    if let Curve2d::BSpline(bs) = c2d {
        let a_cont = 2i32;
        let mut idx = 1usize;
        let mut i = 0usize;
        let len = bs.knots.len();
        while i < len && idx < n {
            let mut j = i;
            while j < len && (bs.knots[j] - bs.knots[i]).abs() < 1.0e-12 {
                j += 1;
            }
            let mult = (j - i) as i32;
            let k = bs.knots[i];
            if k > first + 1.0e-12 && k < last - 1.0e-12 && mult > a_cont {
                t[idx] = k;
                idx += 1;
            }
            i = j;
        }
    }
}

/// The Geom2d first/last parameter of a bounded pcurve (the
/// Geom2d_BoundedCurve StartPoint/EndPoint parameters of SelectForwardSeam).
fn curve2d_first_param(c: &Curve2d) -> f64 {
    Curve2dEval::default_domain(c)[0]
}

fn curve2d_last_param(c: &Curve2d) -> f64 {
    Curve2dEval::default_domain(c)[1]
}
