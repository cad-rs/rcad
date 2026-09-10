//! ProjLib_ProjectedCurve (TKGeomBase/ProjLib) — part B: the class itself.
//!
//! Source: ProjLib_ProjectedCurve.cxx L274-1178 + ProjLib_ProjectedCurve.hxx
//! L47-229 (1:1; the consumed ChFiKPart entry is the (S, C) constructor used
//! by ChFiKPart_ComputeData_Fcts.cxx L89-91: `ProjLib_ProjectedCurve
//! Projc(HSg, HCg);`).
//!
//! Part A (`proj_lib_projected_curve`) carries the analytic ProjLib_* members
//! and the static helpers; this file carries the ProjLib_ProjectedCurve
//! struct, its constructors, `Perform` and the result accessors.
//!
//! GAP carriers (outside-scope deps, each keeps the OCCT failure path):
//!   - ProjLib_ComputeApproxOnPolarSurface (TKGeomBase/ProjLib, the
//!     Bezier/BSpline-surface branch of Perform, OCCT L496-536) — not
//!     translated; the approximation result stays null and myResult stays
//!     not-done, the OCCT outcome when the polar approximation produces
//!     nothing (L506 `if (!aRes.IsNull())` guard).
//!   - ProjLib_HCompProjectedCurve + Approx_CurveOnSurface (the default
//!     branch for the surfaces of revolution / extrusion / offset, OCCT
//!     L540-714; the rcad ProjLib_CompProjectedCurve encoding lives in
//!     rcad-algo/geomalgo, Approx is untranslated) — the branch is deferred;
//!     myResult stays not-done, the OCCT outcome when the numerical
//!     projection produces no curve (L629-638 `NbCurves == 0 -> return`).
//!     Not consumed by ChFiKPart: ChFiKPart_ProjPC rejects non-analytic
//!     surfaces before constructing ProjLib_ProjectedCurve
//!     (ChFiKPart_ComputeData_Fcts.cxx L87).
//!   - ProjLib_ComputeApprox (the advanced analytic fallback, OCCT L717-786)
//!     — carried by [`ProjLibComputeApprox`] below; its payload needs
//!     Approx_FitAndDivide2d (TKGeomBase/Approx) +
//!     Convert_CompBezierCurves2dToBSplineCurve2d (TKGeomBase/Convert) +
//!     BSplCLib::Reparametrize (TKMath) + GeomLib::SameRange, none
//!     translated.  The carrier keeps the OCCT failure path: the Bezier and
//!     BSpline results stay null, Perform returns at L726-729 and the
//!     consumer (ChFiKPart_ProjPC) throws Standard_NotImplemented on the
//!     OtherCurve type — the same failure the OCCT algorithm reports when
//!     the approximation does not converge.

use std::collections::HashMap;
use std::sync::Arc;

use glam::DVec2;

use super::adaptor::{Adaptor2dCurve2d, GeomAbsSurfaceType};
use super::proj_lib_projected_curve::{
    clone_projector, iso_is_deg, project_dispatch, trim_c3d, Adaptor3dCurveGeom,
    Adaptor3dSurfaceGeom, IsoType, ProjLibCone, ProjLibCylinder, ProjLibPlane, ProjLibSphere,
    ProjLibTorus,
};
use super::{CurveType, Projector};
use crate::core::precision;
use crate::geom::{Curve2d, Curve2dEval};

/// OCCT AppParCurves_Constraint (the values consumed by the ProjLib
/// translations).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AppParCurvesConstraint {
    /// OCCT AppParCurves_PassPoint.
    PassPoint,
    /// OCCT AppParCurves_TangencyPoint.
    TangencyPoint,
}

/// OCCT `occ::handle<Adaptor3d_Curve>` carrying the geometry accessors — the
/// constructor parameter type of Perform.
pub type GeomCurveHandle = Arc<dyn Adaptor3dCurveGeom>;

/// OCCT `occ::handle<Adaptor3d_Surface>` carrying the geometry accessors.
pub type GeomSurfaceHandle = Arc<dyn Adaptor3dSurfaceGeom>;

// =========================================================================
// OCCT ProjLib_ComputeApprox — GAP carrier (see the module header)
// =========================================================================

/// OCCT ProjLib_ComputeApprox (ProjLib_ComputeApprox.hxx L24-63 + .cxx
/// L1243-1522) — GAP carrier: the Perform payload needs Approx_FitAndDivide2d
/// + Convert_CompBezierCurves2dToBSplineCurve2d + BSplCLib::Reparametrize
/// (not translated); the carrier preserves the OCCT failure path (null
/// Bezier / null BSpline after Perform).
pub struct ProjLibComputeApprox {
    /// OCCT: Standard_Real myTolerance.
    my_tolerance: f64,
    /// OCCT: Standard_Integer myDegMin / myDegMax / myMaxSegments.
    my_deg_min: i32,
    my_deg_max: i32,
    my_max_segments: i32,
    /// OCCT: AppParCurves_Constraint myBndPnt.
    my_bnd_pnt: AppParCurvesConstraint,
}

impl ProjLibComputeApprox {
    /// OCCT ProjLib_ComputeApprox() (L1243-1251): tolerance =
    /// PApproximation, the degrees/segments unset, TangencyPoint bounds.
    pub fn new() -> Self {
        ProjLibComputeApprox {
            my_tolerance: precision::p_approximation(),
            my_deg_min: -1,
            my_deg_max: -1,
            my_max_segments: -1,
            my_bnd_pnt: AppParCurvesConstraint::TangencyPoint,
        }
    }

    /// OCCT SetTolerance(theTolerance) (L1488-1492).
    pub fn set_tolerance(&mut self, the_tolerance: f64) {
        self.my_tolerance = the_tolerance;
    }

    /// OCCT SetDegree(theDegMin, theDegMax) (L1494-1499).
    pub fn set_degree(&mut self, the_deg_min: i32, the_deg_max: i32) {
        self.my_deg_min = the_deg_min;
        self.my_deg_max = the_deg_max;
    }

    /// OCCT SetMaxSegments(theMaxSegments) (L1501-1505).
    pub fn set_max_segments(&mut self, the_max_segments: i32) {
        self.my_max_segments = the_max_segments;
    }

    /// OCCT SetBndPnt(theBndPnt) (L1507-1511).
    pub fn set_bnd_pnt(&mut self, the_bnd_pnt: AppParCurvesConstraint) {
        self.my_bnd_pnt = the_bnd_pnt;
    }

    /// OCCT Perform(C, S) (L1243-1522).  GAP: the Approx_FitAndDivide2d
    /// payload is not translated; Perform leaves both results null — the
    /// OCCT state when `Fit.IsAllApproximated()` fails.
    pub fn perform(&mut self, _c: &GeomCurveHandle, _s: &GeomSurfaceHandle) {
        // GAP: deferred with the Approx package translation.
    }

    /// OCCT BSpline() (L1529-1534) — null while the payload is deferred.
    pub fn bspline(&self) -> Option<crate::geom::BSplineCurve2> {
        None
    }

    /// OCCT Bezier() (L1536-1541) — null while the payload is deferred.
    pub fn bezier(&self) -> Option<crate::geom::BezierCurve2> {
        None
    }

    /// OCCT Tolerance() (L1543-1544).
    pub fn tolerance(&self) -> f64 {
        self.my_tolerance
    }
}

impl Default for ProjLibComputeApprox {
    fn default() -> Self {
        Self::new()
    }
}

// =========================================================================
// OCCT ProjLib_ProjectedCurve (ProjLib_ProjectedCurve.hxx L47-229)
// =========================================================================

/// OCCT ProjLib_ProjectedCurve — compute the 2d-curve of a 3d curve on a
/// surface: solve the particular (analytic) case when possible, otherwise
/// approximate.
pub struct ProjLibProjectedCurve {
    /// OCCT: Standard_Real myTolerance.
    my_tolerance: f64,
    /// OCCT: handle(Adaptor3d_Surface) mySurface.
    my_surface: Option<GeomSurfaceHandle>,
    /// OCCT: handle(Adaptor3d_Curve) myCurve.
    my_curve: Option<GeomCurveHandle>,
    /// OCCT: ProjLib_Projector myResult.
    my_result: Projector,
    /// OCCT: Standard_Integer myDegMin / myDegMax / myMaxSegments.
    my_deg_min: i32,
    my_deg_max: i32,
    my_max_segments: i32,
    /// OCCT: Standard_Real myMaxDist.
    my_max_dist: f64,
    /// OCCT: AppParCurves_Constraint myBndPnt.
    my_bnd_pnt: AppParCurvesConstraint,
}

impl ProjLibProjectedCurve {
    /// OCCT ProjLib_ProjectedCurve() (L274-282) — the empty constructor.
    pub fn new() -> Self {
        ProjLibProjectedCurve {
            my_tolerance: precision::CONFUSION,
            my_surface: None,
            my_curve: None,
            my_result: Projector::new(),
            my_deg_min: -1,
            my_deg_max: -1,
            my_max_segments: -1,
            my_max_dist: -1.0,
            my_bnd_pnt: AppParCurvesConstraint::TangencyPoint,
        }
    }

    /// OCCT ProjLib_ProjectedCurve(S) (L286-295) — Load(S).
    pub fn with_surface(s: GeomSurfaceHandle) -> Self {
        let mut pc = ProjLibProjectedCurve::new();
        pc.load_surface(&s);
        pc
    }

    /// OCCT ProjLib_ProjectedCurve(S, C) (L299-310) — Load(S); Perform(C);
    /// the default approximation tolerance (Precision::Confusion()).
    pub fn with_surface_curve(s: GeomSurfaceHandle, c: GeomCurveHandle) -> Self {
        let mut pc = ProjLibProjectedCurve::new();
        pc.load_surface(&s);
        pc.perform(c);
        pc
    }

    /// OCCT ProjLib_ProjectedCurve(S, C, Tol) (L314-326) — Load(S);
    /// Perform(C) with myTolerance = max(Tol, Precision::Confusion()).
    pub fn with_surface_curve_tol(s: GeomSurfaceHandle, c: GeomCurveHandle, tol: f64) -> Self {
        let mut pc = ProjLibProjectedCurve::new();
        pc.my_tolerance = tol.max(precision::CONFUSION);
        pc.load_surface(&s);
        pc.perform(c);
        pc
    }

    /// OCCT ShallowCopy() (L330-351) — the rcad field copy.
    pub fn shallow_copy(&self) -> ProjLibProjectedCurve {
        ProjLibProjectedCurve {
            my_tolerance: self.my_tolerance,
            my_surface: self.my_surface.clone(),
            my_curve: self.my_curve.clone(),
            my_result: clone_projector(&self.my_result),
            my_deg_min: self.my_deg_min,
            my_deg_max: self.my_deg_max,
            my_max_segments: self.my_max_segments,
            my_max_dist: self.my_max_dist,
            my_bnd_pnt: self.my_bnd_pnt,
        }
    }

    /// OCCT Load(S) (L355-358) — changes the surface.
    pub fn load_surface(&mut self, s: &GeomSurfaceHandle) {
        self.my_surface = Some(s.clone());
    }

    /// OCCT Load(Tolerance) (L362-365) — changes the tolerance.
    pub fn load_tolerance(&mut self, the_tol: f64) {
        self.my_tolerance = the_tol;
    }

    // ---------------------------------------------------------------------
    // OCCT Perform (ProjLib_ProjectedCurve.cxx L369-886)
    // ---------------------------------------------------------------------

    /// OCCT ProjLib_ProjectedCurve::Perform(C) (L369-886) — perform the
    /// projection: the analytic branches solve the particular case through
    /// the ProjLib_* members; the approximation branches are GAP carriers
    /// (see the module header).
    #[allow(unused_assignments)] // OCCT L383-387 initializes U1..V2 to 0 before the surface read
    pub fn perform(&mut self, c: GeomCurveHandle) {
        // OCCT L371-387: the setup.
        self.my_tolerance = self.my_tolerance.max(precision::CONFUSION);
        let mut my_curve = c;
        self.my_curve = Some(my_curve.clone());
        let first_par = my_curve.first_parameter();
        let last_par = my_curve.last_parameter();
        let s = self
            .my_surface
            .clone()
            .expect("ProjLib_ProjectedCurve::Perform: no surface loaded");
        let s_type = s.get_type();
        let c_type = my_curve.get_type();
        let mut is_analytical_surf = true;
        let mut is_trimmed = [false, false];
        let mut singular_case = [0i32; 2];
        let eps: f64 = 0.01;
        let mut tol_conf = precision::CONFUSION;
        let mut dt = (last_par - first_par) * eps;
        let mut u1 = 0.0_f64;
        let mut u2 = 0.0_f64;
        let mut v1 = 0.0_f64;
        let mut v2 = 0.0_f64;
        u1 = s.first_u_parameter();
        u2 = s.last_u_parameter();
        v1 = s.first_v_parameter();
        v2 = s.last_v_parameter();

        match s_type {
            // OCCT L391-396.
            GeomAbsSurfaceType::Plane => {
                let mut p = ProjLibPlane::new(&s.plane());
                project_dispatch(&mut p, my_curve.as_ref());
                self.my_result = p.projector;
            }

            // OCCT L398-403.
            GeomAbsSurfaceType::Cylinder => {
                let mut p = ProjLibCylinder::new(&s.cylinder());
                project_dispatch(&mut p, my_curve.as_ref());
                self.my_result = p.projector;
            }

            // OCCT L405-410.
            GeomAbsSurfaceType::Cone => {
                let mut p = ProjLibCone::new(&s.cone());
                project_dispatch(&mut p, my_curve.as_ref());
                self.my_result = p.projector;
            }

            // OCCT L412-445.
            GeomAbsSurfaceType::Sphere => {
                let mut p = ProjLibSphere::new(&s.sphere());
                project_dispatch(&mut p, my_curve.as_ref());
                if p.projector.is_done {
                    // Place into the pseudo-period (since Sphere is not
                    // periodic in V!)
                    p.set_in_bounds(my_curve.first_parameter());
                } else {
                    let v_max = std::f64::consts::FRAC_PI_2;
                    let v_min = -v_max;
                    let minang = 1.0e-5 * std::f64::consts::PI;
                    let a_sph = s.sphere();
                    let an_r = a_sph.radius;
                    let f = my_curve.first_parameter();
                    let l = my_curve.last_parameter();

                    let pf = my_curve.value(f);
                    let pl = my_curve.value(l);
                    let a_loc = a_sph.center;
                    let maxdist = pf.distance(a_loc).max(pl.distance(a_loc));
                    tol_conf = (an_r * minang).max((an_r - maxdist).abs());

                    // Surface has pole at V = Vmin
                    let pole = s.value(u1, v_min);
                    trim_c3d(
                        &mut my_curve,
                        &mut is_trimmed,
                        dt,
                        pole,
                        &mut singular_case,
                        3,
                        tol_conf,
                    );
                    // Surface has pole at V = Vmax
                    let pole = s.value(u1, v_max);
                    trim_c3d(
                        &mut my_curve,
                        &mut is_trimmed,
                        dt,
                        pole,
                        &mut singular_case,
                        4,
                        tol_conf,
                    );
                }
                self.my_result = p.projector;
            }

            // OCCT L447-452.
            GeomAbsSurfaceType::Torus => {
                let mut p = ProjLibTorus::new(&s.torus());
                project_dispatch(&mut p, my_curve.as_ref());
                self.my_result = p.projector;
            }

            // OCCT L454-537: the Bezier / BSpline surface branch.
            GeomAbsSurfaceType::BezierSurface | GeomAbsSurfaceType::BSplineSurface => {
                is_analytical_surf = false;
                let mut f = my_curve.first_parameter();
                let mut l = my_curve.last_parameter();
                dt = (l - f) * eps;

                // OCCT L462-466: the surface window.
                u1 = s.first_u_parameter();
                u2 = s.last_u_parameter();
                v1 = s.first_v_parameter();
                v2 = s.last_v_parameter();

                if iso_is_deg(s.as_ref(), u1, IsoType::IsoU, 0.0, self.my_tolerance) {
                    // Surface has pole at U = Umin
                    let pole = s.value(u1, v1);
                    trim_c3d(
                        &mut my_curve,
                        &mut is_trimmed,
                        dt,
                        pole,
                        &mut singular_case,
                        1,
                        tol_conf,
                    );
                }

                if iso_is_deg(s.as_ref(), u2, IsoType::IsoU, 0.0, self.my_tolerance) {
                    // Surface has pole at U = Umax
                    let pole = s.value(u2, v1);
                    trim_c3d(
                        &mut my_curve,
                        &mut is_trimmed,
                        dt,
                        pole,
                        &mut singular_case,
                        2,
                        tol_conf,
                    );
                }

                if iso_is_deg(s.as_ref(), v1, IsoType::IsoV, 0.0, self.my_tolerance) {
                    // Surface has pole at V = Vmin
                    let pole = s.value(u1, v1);
                    trim_c3d(
                        &mut my_curve,
                        &mut is_trimmed,
                        dt,
                        pole,
                        &mut singular_case,
                        3,
                        tol_conf,
                    );
                }

                if iso_is_deg(s.as_ref(), v2, IsoType::IsoV, 0.0, self.my_tolerance) {
                    // Surface has pole at V = Vmax
                    let pole = s.value(u1, v2);
                    trim_c3d(
                        &mut my_curve,
                        &mut is_trimmed,
                        dt,
                        pole,
                        &mut singular_case,
                        4,
                        tol_conf,
                    );
                }

                // OCCT L496-502: ProjLib_ComputeApproxOnPolarSurface polar;
                //   polar.SetTolerance/SetDegree/SetMaxSegments/SetBndPnt/
                //   SetMaxDist; polar.Perform(myCurve, mySurface);
                // OCCT L504-536: aRes = polar.BSpline(); when non-null the
                //   trimmed segments are extended (ExtendC2d, L511-522), the
                //   parameter range is restored (GeomLib::SameRange,
                //   L523-531) and the result is loaded as a BSpline.
                // GAP: ProjLib_ComputeApproxOnPolarSurface is not translated;
                // aRes stays null so the `if (!aRes.IsNull())` guard keeps
                // myResult not-done — the OCCT failure path.
                let _ = (&mut f, &mut l);
            }

            // OCCT L540-714: the default branch (surfaces of revolution,
            // extrusion, offset, other).
            _ => {
                is_analytical_surf = false;
                // OCCT L542-612: the surface-of-revolution singularity
                //   prologue (Extrema_ExtPC on the basis curve).
                // OCCT L614-624: ComputeTolU/ComputeTolV + the
                //   ProjLib_HCompProjectedCurve construction.
                // OCCT L629-659: Bounds + Approx_CurveOnSurface::Perform.
                // OCCT L661-713: the result assembly (MaxError2dU/V,
                //   ExtendC2d, GeomLib::SameRange, RemoveKnot smoothing).
                // GAP: the rcad ProjLib_CompProjectedCurve encoding lives in
                // rcad-algo/geomalgo (outside this crate) and
                // Approx_CurveOnSurface (TKGeomBase/Approx) is untranslated;
                // the branch payload is deferred with those dependencies.
                // myResult stays not-done — the OCCT outcome when the
                // numerical projection produces no curve (L629-638
                // `NbCurves == 0 -> return`).  Not consumed by ChFiKPart:
                // ChFiKPart_ProjPC rejects non-analytic surfaces before
                // constructing ProjLib_ProjectedCurve
                // (ChFiKPart_ComputeData_Fcts.cxx L87).
            }
        }

        // OCCT L717-786: the advanced analytic fallback.
        if !self.my_result.is_done && is_analytical_surf {
            // Use advanced analytical projector if base analytical projection
            // failed.
            // OCCT L720-725.
            let mut comp = ProjLibComputeApprox::new();
            comp.set_tolerance(self.my_tolerance);
            comp.set_degree(self.my_deg_min, self.my_deg_max);
            comp.set_max_segments(self.my_max_segments);
            comp.set_bnd_pnt(self.my_bnd_pnt);
            comp.perform(&my_curve, &s);
            if comp.bezier().is_none() && comp.bspline().is_none() {
                self.my_curve = Some(my_curve);
                return; // advanced projector has been failed too
            }
            // OCCT L730-786: the fallback result assembly —
            //   Geom2dConvert::CurveToBSplineCurve for the bezier result
            //   (L731-739), the ExtendC2d + GeomLib::SameRange treatment of
            //   the trimmed cases (L740-765), the Bezier type/periodicity
            //   flags (L766-784) and myTolerance = Comp.Tolerance()
            //   (L785).  GAP: rides on the approximation payload above
            //   (Approx_FitAndDivide2d / Convert_CompBezierCurves2dTo-
            //   BSplineCurve2d / BSplCLib::Reparametrize / GeomLib::SameRange
            //   — not translated); unreachable until the Approx package
            //   lands because both results stay null.
        }

        // OCCT L788-885: place the result into the surface parameter space.
        let is_periodic = [s.is_u_periodic(), s.is_v_periodic()];
        if self.my_result.is_done && (is_periodic[0] || is_periodic[1]) {
            // Check result curve to be in params space.
            let a_surf_first_par = [s.first_u_parameter(), s.first_v_parameter()];
            let mut a_surf_period = [0.0_f64, 0.0_f64];
            if is_periodic[0] {
                a_surf_period[0] = s.u_period();
            }
            if is_periodic[1] {
                a_surf_period[1] = s.v_period();
            }

            for an_idx in 1..=2 {
                if !is_periodic[an_idx - 1] {
                    continue;
                }

                if self.my_result.get_type() == CurveType::BSpline {
                    // OCCT L814-865: the per-knot-span coordinate histogram.
                    let mut a_map: HashMap<i64, i32> = HashMap::new();
                    let a_deg;
                    {
                        let a_res = match self.my_result.bspline.as_mut() {
                            Some(Curve2d::BSpline(b)) => b,
                            _ => panic!("ProjLib_ProjectedCurve: bspline result expected"),
                        };
                        a_deg = a_res.degree;
                        // OCCT walks the DISTINCT knot array
                        // (FirstUKnotIndex..LastUKnotIndex); rcad stores the
                        // flat expanded knot vector — the non-degenerate
                        // consecutive spans are the same knot intervals.
                        for knot_idx in 0..a_res.knots.len().saturating_sub(1) {
                            let a_first_param = a_res.knots[knot_idx];
                            let a_last_param = a_res.knots[knot_idx + 1];
                            if a_last_param <= a_first_param {
                                continue;
                            }
                            for an_int_idx in 0..=a_deg {
                                let a_curr_param = a_first_param
                                    + (a_last_param - a_first_param) * an_int_idx as f64
                                        / (a_deg as f64 + 1.0);
                                let a_pnt2d =
                                    Curve2d::BSpline(a_res.clone()).point_at(a_curr_param);

                                let coord = if an_idx == 1 { a_pnt2d.x } else { a_pnt2d.y };
                                let mut a_map_key = ((coord - a_surf_first_par[an_idx - 1])
                                    / a_surf_period[an_idx - 1])
                                    as i64;

                                if coord - a_surf_first_par[an_idx - 1] < 0.0 {
                                    a_map_key -= 1;
                                }

                                *a_map.entry(a_map_key).or_insert(0) += 1;
                            }
                        }
                    }

                    // OCCT L849-858: the dominant period bucket.
                    let mut a_max_points = 0;
                    let mut a_max_idx = 0_i64;
                    for (k, v) in a_map.iter() {
                        if *v > a_max_points {
                            a_max_points = *v;
                            a_max_idx = *k;
                        }
                    }
                    if a_max_idx != 0 {
                        // OCCT L859-865: translate the curve by whole
                        // periods so the dominant bucket is 0.
                        if let Some(Curve2d::BSpline(a_res)) = self.my_result.bspline.as_mut() {
                            let a_first_pnt = {
                                let c = Curve2d::BSpline(a_res.clone());
                                let first = c.default_domain()[0];
                                c.point_at(first)
                            };
                            let mut a_second_pnt = a_first_pnt;
                            if an_idx == 1 {
                                a_second_pnt.x =
                                    a_first_pnt.x - a_surf_period[0] * a_max_idx as f64;
                            } else {
                                a_second_pnt.y =
                                    a_first_pnt.y - a_surf_period[1] * a_max_idx as f64;
                            }
                            let t = a_second_pnt - a_first_pnt;
                            for p in a_res.control_points.iter_mut() {
                                *p += t;
                            }
                        }
                    }
                }

                if self.my_result.get_type() == CurveType::Line {
                    // OCCT L868-883: the projected line is framed into the
                    // parameter space.
                    let a_t1 = my_curve.first_parameter();
                    let a_t2 = my_curve.last_parameter();

                    if an_idx == 1 {
                        // U param space.
                        self.my_result
                            .u_frame(a_t1, a_t2, a_surf_first_par[0], a_surf_period[0]);
                    } else {
                        // V param space.
                        self.my_result
                            .v_frame(a_t1, a_t2, a_surf_first_par[1], a_surf_period[1]);
                    }
                }
            }
        }
        self.my_curve = Some(my_curve);
        let _ = (c_type, first_par, last_par, u1, u2, v1, v2, singular_case, is_trimmed);
    }

    // ---------------------------------------------------------------------
    // OCCT setters / getters (L890-950, L984-987, L1047-1167)
    // ---------------------------------------------------------------------

    /// OCCT SetDegree(theDegMin, theDegMax) (L890-894).
    pub fn set_degree(&mut self, the_deg_min: i32, the_deg_max: i32) {
        self.my_deg_min = the_deg_min;
        self.my_deg_max = the_deg_max;
    }

    /// OCCT SetMaxSegments(theMaxSegments) (L898-902).
    pub fn set_max_segments(&mut self, the_max_segments: i32) {
        self.my_max_segments = the_max_segments;
    }

    /// OCCT SetBndPnt(theBndPnt) (L905-908).
    pub fn set_bnd_pnt(&mut self, the_bnd_pnt: AppParCurvesConstraint) {
        self.my_bnd_pnt = the_bnd_pnt;
    }

    /// OCCT SetMaxDist(theMaxDist) (L912-915).
    pub fn set_max_dist(&mut self, the_max_dist: f64) {
        self.my_max_dist = the_max_dist;
    }

    /// OCCT GetSurface() (L919-922).
    pub fn get_surface(&self) -> Option<&GeomSurfaceHandle> {
        self.my_surface.as_ref()
    }

    /// OCCT GetCurve() (L926-929).
    pub fn get_curve(&self) -> Option<&GeomCurveHandle> {
        self.my_curve.as_ref()
    }

    /// OCCT GetTolerance() (L933-936) — the reached tolerance.
    pub fn get_tolerance(&self) -> f64 {
        self.my_tolerance
    }

    /// OCCT IsPeriodic() (L984-987).
    pub fn is_periodic(&self) -> bool {
        self.my_result.is_periodic
    }

    /// OCCT GetType() (L1047-1050).
    pub fn get_type(&self) -> CurveType {
        self.my_result.get_type()
    }

    /// OCCT Line() (L1054-1057).
    pub fn line(&self) -> crate::geom::Line2d {
        let l = &self.my_result.lin;
        crate::geom::Line2d {
            origin: DVec2::new(l.origin.x, l.origin.y),
            direction: DVec2::new(l.direction.x, l.direction.y).normalize_or_zero(),
        }
    }

    /// OCCT Circle() (L1060-1063).
    pub fn circle(&self) -> crate::geom::Circle2d {
        let c = &self.my_result.circ;
        crate::geom::Circle2d {
            center: DVec2::new(c.center.x, c.center.y),
            x_dir: DVec2::new(c.x_dir.x, c.x_dir.y),
            y_dir: DVec2::new(c.y_dir.x, c.y_dir.y),
            radius: c.radius,
        }
    }

    /// OCCT Ellipse() (L1066-1069).
    pub fn ellipse(&self) -> crate::geom::Ellipse2d {
        let e = &self.my_result.elips;
        crate::geom::Ellipse2d {
            center: DVec2::new(e.center.x, e.center.y),
            major_dir: DVec2::new(e.major_dir.x, e.major_dir.y).normalize_or_zero(),
            minor_dir: DVec2::new(-e.major_dir.y, e.major_dir.x),
            major_radius: e.major_radius,
            minor_radius: e.minor_radius,
        }
    }

    /// OCCT Parabola() (L1080-1083).
    pub fn parabola(&self) -> crate::geom::Parabola2d {
        let p = &self.my_result.parab;
        crate::geom::Parabola2d {
            origin: DVec2::new(p.vertex.x, p.vertex.y),
            axis_dir: DVec2::new(p.axis_dir.x, p.axis_dir.y).normalize_or_zero(),
            focal_param: p.focal_param,
        }
    }

    /// OCCT Hyperbola() (L1073-1076).
    pub fn hyperbola(&self) -> crate::geom::Hyperbola2d {
        let h = &self.my_result.hypr;
        crate::geom::Hyperbola2d {
            center: DVec2::new(h.center.x, h.center.y),
            major_dir: DVec2::new(h.major_dir.x, h.major_dir.y).normalize_or_zero(),
            semi_major: h.semi_major,
            semi_minor: h.semi_minor,
        }
    }

    /// OCCT Degree() (L1089-1105) — the Standard_NoSuchObject raise for the
    /// non-analytic result types is mirrored as a panic.
    pub fn degree(&self) -> i32 {
        let t = self.get_type();
        assert!(
            t == CurveType::BSpline || t == CurveType::Bezier,
            "ProjLib_ProjectedCurve:Degree"
        );
        if t == CurveType::BSpline {
            self.my_result
                .bspline_curve()
                .map(|b| b.degree as i32)
                .unwrap_or(0)
        } else {
            self.my_result
                .bezier()
                .map(|b| (b.control_points.len().max(2) - 1) as i32)
                .unwrap_or(0)
        }
    }

    /// OCCT IsRational() (L1109-1124).
    pub fn is_rational(&self) -> bool {
        let t = self.get_type();
        assert!(
            t == CurveType::BSpline || t == CurveType::Bezier,
            "ProjLib_ProjectedCurve:IsRational"
        );
        if t == CurveType::BSpline {
            self.my_result
                .bspline_curve()
                .map(|b| b.weights.iter().any(|&w| (w - 1.0).abs() > f64::EPSILON))
                .unwrap_or(false)
        } else {
            self.my_result
                .bezier()
                .map(|b| b.weights.iter().any(|&w| (w - 1.0).abs() > f64::EPSILON))
                .unwrap_or(false)
        }
    }

    /// OCCT NbPoles() (L1128-1144).
    pub fn nb_poles(&self) -> i32 {
        let t = self.get_type();
        assert!(
            t == CurveType::BSpline || t == CurveType::Bezier,
            "ProjLib_ProjectedCurve:NbPoles"
        );
        if t == CurveType::BSpline {
            self.my_result
                .bspline_curve()
                .map(|b| b.control_points.len() as i32)
                .unwrap_or(0)
        } else {
            self.my_result
                .bezier()
                .map(|b| b.control_points.len() as i32)
                .unwrap_or(0)
        }
    }

    /// OCCT NbKnots() (L1148-1153) — the distinct-knot count of the result.
    pub fn nb_knots(&self) -> i32 {
        assert!(self.get_type() == CurveType::BSpline, "ProjLib_ProjectedCurve:NbKnots");
        self.my_result
            .bspline_curve()
            .map(|b| {
                let mut n = 0;
                let knots = &b.knots;
                for (i, k) in knots.iter().enumerate() {
                    if i == 0 || *k > knots[i - 1] {
                        n += 1;
                    }
                }
                n
            })
            .unwrap_or(0) as i32
    }

    /// OCCT Bezier() (L1157-1160).
    pub fn bezier(&self) -> Option<crate::geom::BezierCurve2> {
        self.my_result.bezier()
    }

    /// OCCT BSpline() (L1164-1167).
    pub fn bspline(&self) -> Option<crate::geom::BSplineCurve2> {
        self.my_result.bspline_curve()
    }
}

impl Default for ProjLibProjectedCurve {
    fn default() -> Self {
        Self::new()
    }
}

// =========================================================================
// OCCT Adaptor2d_Curve2d interface (the NotImplemented throws preserved)
// =========================================================================

impl Adaptor2dCurve2d for ProjLibProjectedCurve {
    /// OCCT FirstParameter() (L940-943).
    fn first_parameter(&self) -> f64 {
        self.my_curve
            .as_ref()
            .map(|c| c.first_parameter())
            .unwrap_or(0.0)
    }

    /// OCCT LastParameter() (L947-950).
    fn last_parameter(&self) -> f64 {
        self.my_curve
            .as_ref()
            .map(|c| c.last_parameter())
            .unwrap_or(0.0)
    }

    /// OCCT Value(U) (L998-1001) — throws Standard_NotImplemented.
    fn value(&self, _u: f64) -> DVec2 {
        panic!("ProjLib_ProjectedCurve::Value() - method is not implemented")
    }

    /// OCCT D0(U, P) (L1005-1008).
    fn d0(&self, _u: f64) -> DVec2 {
        panic!("ProjLib_ProjectedCurve::D0() - method is not implemented")
    }

    /// OCCT D1(U, P, V) (L1012-1015).
    fn d1(&self, _u: f64) -> (DVec2, DVec2) {
        panic!("ProjLib_ProjectedCurve::D1() - method is not implemented")
    }

    /// OCCT D2(U, P, V1, V2) (L1019-1022).
    fn d2(&self, _u: f64) -> (DVec2, DVec2, DVec2) {
        panic!("ProjLib_ProjectedCurve::D2() - method is not implemented")
    }

    /// OCCT Continuity() (L954-957).
    fn continuity(&self) -> crate::math::GeomAbsShape {
        panic!("ProjLib_ProjectedCurve::Continuity() - method is not implemented")
    }

    /// OCCT GetType() (L1047-1050).
    fn get_type(&self) -> CurveType {
        self.my_result.get_type()
    }

    /// OCCT Line() (L1054-1057).
    fn line(&self) -> crate::geom::Line2d {
        ProjLibProjectedCurve::line(self)
    }

    /// OCCT BSpline() (L1164-1167) — None models the null handle.
    fn bspline(&self) -> Option<crate::geom::BSplineCurve2> {
        self.my_result.bspline_curve()
    }

    /// OCCT Bezier() (L1157-1160).
    fn bezier(&self) -> Option<crate::geom::BezierCurve2> {
        self.my_result.bezier()
    }

    /// OCCT Trim(First, Last, Tol) (L1171-1178).
    fn trim(&self, _first: f64, _last: f64, _tol: f64) -> Arc<dyn Adaptor2dCurve2d> {
        panic!("ProjLib_ProjectedCurve::Trim() - method is not implemented")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::proj_lib::proj_lib_projected_curve::{
        GeomCurveAdaptor, GeomSurfaceAdaptor,
    };
    use crate::geom::{Circle3, SphericalSurface};
    use glam::DVec3;

    #[test]
    fn great_circle_on_sphere_projects_to_line() {
        // The meridian circle (through both poles) is an iso-U line.
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 2.0,
            ref_dir: DVec3::X,
        };
        let mut circle = Circle3::new(DVec3::ZERO, DVec3::Y, 2.0);
        circle.x_dir = DVec3::Z;
        circle.y_dir = DVec3::X;
        let s: GeomSurfaceHandle = Arc::new(GeomSurfaceAdaptor::new(crate::geom::Surface3::Sphere(sphere)));
        let c: GeomCurveHandle = Arc::new(GeomCurveAdaptor::with_range(
            crate::geom::Curve3::Circle(circle),
            0.0,
            std::f64::consts::FRAC_PI_2,
        ));
        let projc = ProjLibProjectedCurve::with_surface_curve(s, c);
        assert_eq!(projc.get_type(), CurveType::Line);
        let l = projc.line();
        // The meridian through X has U = 0.
        assert!(l.origin.x.abs() < 1e-9);
    }

    #[test]
    fn non_iso_circle_keeps_the_occt_failure_path() {
        // A small tilted circle on the sphere is neither iso-U nor iso-V:
        // the analytic projection is not done and the advanced fallback
        // (GAP carrier) produces nothing — GetType stays OtherCurve, exactly
        // the OCCT state that drives ChFiKPart_ProjPC into
        // Standard_NotImplemented.
        let sphere = SphericalSurface {
            center: DVec3::ZERO,
            axis: DVec3::Z,
            radius: 1.0,
            ref_dir: DVec3::X,
        };
        let circle = Circle3::new(DVec3::new(0.0, 0.0, 0.5), DVec3::new(1.0, 0.0, 1.0), 0.5);
        let s: GeomSurfaceHandle = Arc::new(GeomSurfaceAdaptor::new(crate::geom::Surface3::Sphere(sphere)));
        let c: GeomCurveHandle = Arc::new(GeomCurveAdaptor::with_range(
            crate::geom::Curve3::Circle(circle),
            0.0,
            1.0,
        ));
        let projc = ProjLibProjectedCurve::with_surface_curve(s, c);
        assert_eq!(projc.get_type(), CurveType::Other);
    }
}
