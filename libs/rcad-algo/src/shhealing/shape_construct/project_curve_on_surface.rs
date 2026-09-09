//! OCCT ShapeConstruct_ProjectCurveOnSurface (TKShHealing
//! ShapeConstruct_ProjectCurveOnSurface.hxx L52-298 + .cxx L639-1102) —
//! 1:1 translation.
//!
//! This tool provides a method for computing pcurve by projecting a 3d curve
//! onto a surface.  Projection is done by 23 or more points (this number is
//! changed for B-Splines according to the rule: the total number of the
//! points is not less than number of spans * (degree + 1); it is increased
//! recursively starting with 23 and is added with 22 until the condition is
//! fulfilled).  Isoparametric cases (if curve corresponds to U=const or
//! V=const on the surface) are recognized with the given precision.
//!
//! The protected approximation members (approxPCurve and below, .cxx
//! L1104-2785) live in `project_curve_on_surface_approx.rs`; the file-local
//! helpers of the .cxx anonymous namespace live in
//! `project_curve_on_surface_ns.rs`.

use glam::{DVec2, DVec3};
use rcad_kernel::base::geom_proj_lib::project_on_plane::project_on_plane;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve::GeomCurveAdaptor;
use rcad_kernel::base::proj_lib::proj_lib_projected_curve_b::{
    GeomCurveHandle, ProjLibProjectedCurve,
};
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::geom::{
    BSplineCurve2, BSplineCurve3, Curve2d, Curve3, CurveEval, Line2d, Line3, Surface3,
    SurfaceEval, TrimmedCurve3,
};
use rcad_kernel::precision::{CONFUSION, PCONFUSION, SQUARE_CONFUSION};
use std::sync::Arc;

use crate::shhealing::shape_extend::status::{decode_status, encode_status, ShapeExtendStatus};

use super::array1::Array1;
use super::gap_deps::{surface_is_cn_u, surface_is_cn_v, surface_u_period, surface_v_period,
    ShapeAnalysisSurface};
use super::project_curve_on_surface_ns::{
    fix_periodicity_troubles, generate_curve_points, is_bspline_curve_invalid,
    SurfaceProjectorWithCache, THE_NCONTROL,
};

/// OCCT .hxx L55: `using ArrayOfPnt   = NCollection_Array1<gp_Pnt>`.
pub type ArrayOfPnt = Array1<DVec3>;
/// OCCT .hxx L57: `using ArrayOfPnt2d = NCollection_Array1<gp_Pnt2d>`.
pub type ArrayOfPnt2d = Array1<DVec2>;
/// OCCT .hxx L58: `using ArrayOfReal  = NCollection_Array1<double>`.
pub type ArrayOfReal = Array1<f64>;

/// OCCT `theC3D->IsKind(STANDARD_TYPE(Geom_BoundedCurve))` over the rcad
/// `Curve3` enum (arch. diff.): bounded curves are Bezier, BSpline and
/// TrimmedCurve; the analytic conics/lines are unbounded.
pub(crate) fn is_bounded_curve(the_c3d: &Curve3) -> bool {
    matches!(
        the_c3d,
        Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)
    )
}

/// OCCT `theSurf->IsKind(STANDARD_TYPE(Geom_SphericalSurface))` (arch. diff.).
pub(crate) fn is_spherical_surface(the_surf: &Surface3) -> bool {
    matches!(the_surf, Surface3::Sphere(_))
}

/// OCCT .hxx L52-298: `class ShapeConstruct_ProjectCurveOnSurface`.
///
/// This tool provides a method for computing pcurve by projecting
/// 3d curve onto a surface.
pub struct ProjectCurveOnSurface {
    /// OCCT mySurf — surface to project on (handle; None = null handle).
    pub(super) my_surf: Option<ShapeAnalysisSurface>,
    /// OCCT myPreci — current precision.
    pub(super) my_preci: f64,
    /// OCCT myStatus — operation status.
    pub(super) my_status: i32,
    /// OCCT myAdjustOverDegen — seam adjustment flag.
    pub(super) my_adjust_over_degen: i32,
    /// OCCT myCache — cached 3D/2D point pairs for projection optimization
    /// (0-based 2-element cache array).
    pub(super) my_cache: Array1<super::project_curve_on_surface_ns::CachePoint>,
}

impl ProjectCurveOnSurface {
    /// OCCT .cxx L641-646: empty constructor.
    pub fn new() -> Self {
        ProjectCurveOnSurface {
            my_surf: None,
            my_preci: CONFUSION,
            my_status: encode_status(ShapeExtendStatus::Ok),
            my_adjust_over_degen: 1,
            my_cache: Array1::empty(),
        }
    }

    /// OCCT .cxx L650-654: Init(theSurf: Geom_Surface, thePreci).
    pub fn init(&mut self, the_surf: Surface3, the_preci: f64) {
        self.init_sa(ShapeAnalysisSurface::new(the_surf), the_preci);
    }

    /// OCCT .cxx L658-663: Init(theSurf: ShapeAnalysis_Surface, thePreci)
    /// (the second OCCT Init overload; renamed `init_sa` — Rust has no
    /// overloading).
    pub fn init_sa(&mut self, the_surf: ShapeAnalysisSurface, the_preci: f64) {
        self.set_surface_sa(the_surf);
        self.set_precision(the_preci);
    }

    /// OCCT .cxx L667-670: SetSurface(theSurf: Geom_Surface).
    pub fn set_surface(&mut self, the_surf: Surface3) {
        self.set_surface_sa(ShapeAnalysisSurface::new(the_surf));
    }

    /// OCCT .cxx L674-683: SetSurface(theSurf: ShapeAnalysis_Surface) (the
    /// second OCCT overload; renamed `set_surface_sa` — Rust has no
    /// overloading).
    ///
    /// OCCT first compares `mySurf == theSurf` — a handle identity test
    /// (arch. diff.: the rcad carrier has no object identity), so the
    /// fast-path return is unrepresentable and the assignment always runs.
    pub fn set_surface_sa(&mut self, the_surf: ShapeAnalysisSurface) {
        self.my_surf = Some(the_surf);
        self.my_cache = Array1::empty();
    }

    /// OCCT .cxx L687-690: SetPrecision(thePreci).
    pub fn set_precision(&mut self, the_preci: f64) {
        self.my_preci = the_preci;
    }

    /// OCCT .cxx L694-697: AdjustOverDegenMode() — returns (modifiable) the
    /// flag specifying to which side of the parametrical space the part of
    /// the pcurve which lies on the seam is adjusted.
    pub fn adjust_over_degen_mode(&mut self) -> &mut i32 {
        &mut self.my_adjust_over_degen
    }

    /// OCCT .cxx L701-704: Status(theStatus) — the status of the last Perform.
    pub fn status(&self, the_status: ShapeExtendStatus) -> bool {
        decode_status(self.my_status, the_status)
    }

    /// OCCT .cxx L708-775: Perform(theC3D, theFirst, theLast, theC2D,
    /// theTolFirst, theTolLast).
    ///
    /// Computes the projection of a 3d curve onto a surface using the
    /// specialized algorithm.  Returns False if the projector fails,
    /// otherwise, if the pcurve is computed successfully, returns True.  The
    /// output curve 2D is guaranteed to be same-parameter with the input
    /// curve 3D on the interval [theFirst, theLast].
    ///
    /// OCCT defaults theTolFirst/theTolLast to Precision::Confusion(); Rust
    /// has no default arguments — the callers pass CONFUSION.
    ///
    /// Arch. diff.: the OCCT `theC3D.IsNull()` arm (FAIL1 on a null curve
    /// handle) is unrepresentable in a `&Curve3` signature; the mySurf null
    /// arm is kept.
    pub fn perform(
        &mut self,
        the_c3d: &Curve3,
        the_first: f64,
        the_last: f64,
        the_c2d: &mut Option<Curve2d>,
        the_tol_first: f64,
        the_tol_last: f64,
    ) -> bool {
        self.my_status = encode_status(ShapeExtendStatus::Ok);

        if self.my_surf.is_none() {
            *the_c2d = None;
            self.my_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }

        // Try analytical projection first
        let a_curve3d_trim = if !is_bounded_curve(the_c3d) {
            Curve3::Trimmed(TrimmedCurve3::new(the_c3d.clone(), the_first, the_last))
        } else {
            the_c3d.clone()
        };

        *the_c2d = self.project_analytic(&a_curve3d_trim);
        if the_c2d.is_some() {
            self.my_status |= encode_status(ShapeExtendStatus::Done1);
            return true;
        }

        // Handle problematic B-spline curves with uneven parameterization
        let mut a_bspline: Option<BSplineCurve3> = None;
        if is_bspline_curve_invalid(the_c3d, the_first, the_last, &mut a_bspline) {
            // Use ProjLib for curves with highly uneven parameterization speed
            if self.perform_by_proj_lib(the_c3d, the_first, the_last, the_c2d) {
                return self.status(ShapeExtendStatus::Done);
            }
        }

        // Generate curve points
        let mut a_points = ArrayOfPnt::empty();
        let mut a_params = ArrayOfReal::empty();
        let a_nb_pnt = generate_curve_points(
            the_c3d,
            the_first,
            the_last,
            THE_NCONTROL,
            &mut a_points,
            &mut a_params,
        );

        // Approximate pcurve
        let mut a_points2d = ArrayOfPnt2d::new(1, a_nb_pnt as usize, DVec2::ZERO);
        for i in 1..=a_nb_pnt {
            a_points2d.set_value(i as usize, DVec2::ZERO);
        }

        self.approx_p_curve(
            a_nb_pnt,
            the_c3d,
            the_tol_first,
            the_tol_last,
            &mut a_points,
            &mut a_params,
            &mut a_points2d,
            the_c2d,
        );

        let a_nb_pini = a_points.length();
        if the_c2d.is_some() {
            self.my_status |= encode_status(ShapeExtendStatus::Done2);
            return true;
        }

        // Interpolate the result
        *the_c2d = self.interpolate_p_curve(a_nb_pini as i32, &a_points2d, &a_params);

        self.my_status |= encode_status(if the_c2d.is_none() {
            ShapeExtendStatus::Fail1
        } else {
            ShapeExtendStatus::Done2
        });
        self.status(ShapeExtendStatus::Done)
    }

    /// OCCT .cxx L779-845: PerformByProjLib(theC3D, theFirst, theLast, theC2D).
    ///
    /// Computes the projection of a 3d curve onto a surface using the
    /// standard algorithm from ProjLib.  Returns False if the standard
    /// projector fails or raises an exception or cuts the curve by the
    /// parametrical bounds of the surface.
    pub fn perform_by_proj_lib(
        &mut self,
        the_c3d: &Curve3,
        the_first: f64,
        the_last: f64,
        the_c2d: &mut Option<Curve2d>,
    ) -> bool {
        *the_c2d = None;
        let a_gas = match self.my_surf.as_ref() {
            Some(s) => s.adaptor3d(),
            None => {
                self.my_status = encode_status(ShapeExtendStatus::Fail1);
                return false;
            }
        };

        // try { OCC_CATCH_SIGNALS ... } — the OCCT exception arm maps to the
        // catch_unwind Err arm (FAIL3 + nullify).
        let attempted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let a_gac: GeomCurveHandle = Arc::new(GeomCurveAdaptor::with_range(
                the_c3d.clone(),
                the_first,
                the_last,
            ));
            let a_projector = ProjLibProjectedCurve::with_surface_curve(a_gas, a_gac);

            let a_c2d: Option<Curve2d> = match a_projector.get_type() {
                CurveType::Line => Some(Curve2d::Line(a_projector.line())),
                CurveType::Circle => Some(Curve2d::Circle(a_projector.circle())),
                CurveType::Ellipse => Some(Curve2d::Ellipse(a_projector.ellipse())),
                CurveType::Parabola => Some(Curve2d::Parabola(a_projector.parabola())),
                CurveType::Hyperbola => Some(Curve2d::Hyperbola(a_projector.hyperbola())),
                CurveType::BSpline => a_projector.bspline().map(Curve2d::BSpline),
                _ => None,
            };
            a_c2d
        }));

        match attempted {
            Ok(Some(a_c2d)) => {
                *the_c2d = Some(a_c2d);
                self.my_status = encode_status(ShapeExtendStatus::Done1);
                true
            }
            Ok(None) => {
                self.my_status = encode_status(ShapeExtendStatus::Fail2);
                false
            }
            Err(_an_exception) => {
                self.my_status = encode_status(ShapeExtendStatus::Fail3);
                *the_c2d = None;
                false
            }
        }
    }

    /// OCCT .cxx L849-898: projectAnalytic(theC3D) — performs analytical
    /// projection for the plane-surface special cases.
    pub fn project_analytic(&self, the_c3d: &Curve3) -> Option<Curve2d> {
        let mut a_result: Option<Curve2d> = None;

        let a_surf = self.my_surf.as_ref()?.surface().clone();
        let mut a_plane = match &a_surf {
            Surface3::Plane(p) => Some(*p),
            _ => None,
        };

        if a_plane.is_none() {
            match &a_surf {
                Surface3::Trimmed(a_rts) => {
                    if let Surface3::Plane(p) = a_rts.basis.as_ref() {
                        a_plane = Some(*p);
                    }
                }
                Surface3::Offset(an_os) => {
                    if let Surface3::Plane(p) = an_os.basis.as_ref() {
                        a_plane = Some(*p);
                    }
                }
                _ => {}
            }
        }

        if let Some(a_plane) = a_plane {
            // OCCT GeomProjLib::ProjectOnPlane(theC3D, aPlane,
            // aPlane->Position().Direction(), true); the plane placement
            // direction is the plane normal.
            let a_proj_on_plane =
                project_on_plane(the_c3d, the_c3d.default_domain(), &a_plane, a_plane.normal)?;
            let a_gas = self.my_surf.as_ref().unwrap().adaptor3d();
            let a_hc: GeomCurveHandle = Arc::new(GeomCurveAdaptor::new(a_proj_on_plane));
            let a_proj = ProjLibProjectedCurve::with_surface_curve(a_gas, a_hc);

            // OCCT aResult = Geom2dAdaptor::MakeCurve(aProj)
            // (Geom2dAdaptor.cxx L33-117): the conic arms return the bare
            // type curve, the BSpline arm wraps a Geom2d_TrimmedCurve over
            // [First, Last]; the BasisCurve unwrap below strips that trim, so
            // the net BSpline result equals the projector's BSpline.
            a_result = match a_proj.get_type() {
                CurveType::Line => Some(Curve2d::Line(a_proj.line())),
                CurveType::Circle => Some(Curve2d::Circle(a_proj.circle())),
                CurveType::Ellipse => Some(Curve2d::Ellipse(a_proj.ellipse())),
                CurveType::Parabola => Some(Curve2d::Parabola(a_proj.parabola())),
                CurveType::Hyperbola => Some(Curve2d::Hyperbola(a_proj.hyperbola())),
                CurveType::BSpline => a_proj.bspline().map(Curve2d::BSpline),
                _ => None,
            };
            if a_result.is_none() {
                return a_result;
            }

            // OCCT L888-892: if (aResult->IsKind(Geom2d_TrimmedCurve))
            // aResult = aTC->BasisCurve() — covered by the direct BSpline arm
            // above (arch. diff. note).

            return a_result;
        }

        a_result
    }

    /// OCCT .cxx L902-1102: getLine(thePoints, theParams, thePoints2d, theTol,
    /// theIsRecompute, theIsFromCache) — tries to approximate the 3D curve by
    /// a Geom2d_Line or a degree-1 Geom2d_BSplineCurve with the specified
    /// tolerance.  Returns the resulting 2D curve or None if line fitting
    /// failed.
    pub(crate) fn get_line(
        &self,
        the_points: &ArrayOfPnt,
        the_params: &ArrayOfReal,
        the_points2d: &mut ArrayOfPnt2d,
        the_tol: f64,
        the_is_recompute: &mut bool,
        the_is_from_cache: &mut bool,
    ) -> Option<Curve2d> {
        let a_nb = the_points.length();
        let mut a_p: [DVec3; 4] = [DVec3::ZERO; 4];
        a_p[0] = *the_points.value(1);
        a_p[1] = *the_points.value(2);
        a_p[2] = *the_points.value(a_nb - 1);
        a_p[3] = *the_points.value(a_nb);

        let mut a_p2d: [DVec2; 4] = [DVec2::ZERO; 4];

        let mut a_tol2 = the_tol * the_tol;
        let mut a_tol_working = the_tol;
        let is_periodic_u = self
            .my_surf
            .as_ref()
            .map(|s| s.surface().is_u_periodic())
            .unwrap_or(false);
        let is_periodic_v = self
            .my_surf
            .as_ref()
            .map(|s| s.surface().is_v_periodic())
            .unwrap_or(false);

        // Protection against bad tolerance shapes
        if a_tol2 > 1.0 {
            a_tol_working = CONFUSION;
            a_tol2 = a_tol_working * a_tol_working;
        }
        if a_tol2 < SQUARE_CONFUSION {
            a_tol2 = SQUARE_CONFUSION;
        }
        let an_old_tol2 = a_tol2;

        let mut a_saved_point_num: i32 = -1;
        let mut a_saved_point = DVec2::ZERO;

        // Create projector with B-spline surface pole cache optimization
        let a_projector = SurfaceProjectorWithCache::new(self.my_surf.as_ref());

        // OCCT L944-966: helper lambda to project a point with cache lookup.
        // `the_tol2` is passed per call (OCCT captures it by reference; the
        // loop body mutates it between calls).
        let mut project_point = |the_point: DVec3,
                                 the_result: &mut DVec2,
                                 the_index: i32,
                                 the_tol2: f64| {
            // Try existing endpoint cache first
            let a_nb_cache = self.my_cache.length();
            for j in 0..a_nb_cache {
                let a_cache_pnt = self.my_cache.value(j);
                if a_cache_pnt.first.distance_squared(the_point) < the_tol2 {
                    *the_result = a_projector.next_value_of_uv(
                        a_cache_pnt.second,
                        the_point,
                        a_tol_working,
                        the_tol2,
                        a_tol_working,
                    );
                    a_saved_point_num = the_index;
                    a_saved_point = a_cache_pnt.second;
                    if the_index == 0 {
                        *the_is_from_cache = true;
                    }
                    return;
                }
            }

            // Fall back to full projection (with B-spline pole cache
            // optimization inside)
            *the_result = a_projector.value_of_uv(the_point, a_tol_working, the_tol2);
        };

        // Project first and last points
        let mut i = 0usize;
        while i < 4 {
            let mut res = a_p2d[i];
            project_point(a_p[i], &mut res, i as i32, a_tol2);
            a_p2d[i] = res;

            let a_dist = a_projector.gap();
            let a_cur_dist = a_dist * a_dist;
            if a_tol2 < a_cur_dist {
                a_tol2 = a_cur_dist;
            }
            i += 3;
        }

        if is_periodic_u || is_periodic_v {
            // Compute second and last but one c2d points
            let mut i = 1usize;
            while i < 3 {
                let mut res = a_p2d[i];
                project_point(a_p[i], &mut res, i as i32, a_tol2);
                a_p2d[i] = res;

                let a_dist = a_projector.gap();
                let a_cur_dist = a_dist * a_dist;
                if a_tol2 < a_cur_dist {
                    a_tol2 = a_cur_dist;
                }
                i += 1;
            }

            if is_periodic_u {
                *the_is_recompute = fix_periodicity_troubles(
                    &mut a_p2d,
                    1,
                    surface_u_period(self.my_surf.as_ref().unwrap().surface()),
                    a_saved_point_num,
                    a_saved_point.x,
                );
            }

            if is_periodic_v {
                *the_is_recompute = fix_periodicity_troubles(
                    &mut a_p2d,
                    2,
                    surface_v_period(self.my_surf.as_ref().unwrap().surface()),
                    a_saved_point_num,
                    a_saved_point.y,
                );
            }
        }

        if *the_is_recompute && is_spherical_surface(self.my_surf.as_ref().unwrap().surface()) {
            return None;
        }

        the_points2d.set_value(1, a_p2d[0]);
        the_points2d.set_value(a_nb, a_p2d[3]);

        // Restore old tolerance
        a_tol2 = an_old_tol2;

        let d_par = *the_params.value(a_nb) - *the_params.value(1);
        if d_par.abs() < PCONFUSION {
            return None;
        }

        let a_vec0 = a_p2d[3] - a_p2d[0]; // OCCT gp_Vec2d(aP2d[0], aP2d[3])
        let a_vec = a_vec0 / d_par;
        let a_surf = self.my_surf.as_ref().unwrap().surface().clone();
        let mut is_normal_check = surface_is_cn_u(&a_surf, 1) && surface_is_cn_v(&a_surf, 1);

        if is_normal_check {
            for i in 1..=a_nb {
                let a_cur_point =
                    a_p2d[0] + a_vec * (*the_params.value(i) - *the_params.value(1));
                let (a_cur_p, a_du, a_dv) = a_surf.derivatives(a_cur_point.x, a_cur_point.y);
                let a_normal_vec = a_du.cross(a_dv);
                if a_normal_vec.length_squared() < SQUARE_CONFUSION {
                    is_normal_check = false;
                    break;
                }
                let a_normal_line = Line3::new(a_cur_p, a_normal_vec);
                let a_dist = a_normal_line.distance(*the_points.value(i));
                if a_dist > a_tol_working {
                    return None;
                }
            }
        }

        if !is_normal_check {
            let a_first_point_dist = a_surf
                .point_at(a_p2d[0].x, a_p2d[0].y)
                .distance_squared(*the_points.value(1));
            a_tol2 = a_tol2.max(a_tol2 * 2.0 * a_first_point_dist);
            for i in 2..a_nb {
                let a_cur_point =
                    a_p2d[0] + a_vec * (*the_params.value(i) - *the_params.value(1));
                let a_cur_p = a_surf.point_at(a_cur_point.x, a_cur_point.y);
                let a_dist1 = a_cur_p.distance_squared(*the_points.value(i));

                if (a_first_point_dist - a_dist1).abs() > a_tol2 {
                    return None;
                }
            }
        }

        // Check if pcurve can be represented by Geom2d_Line
        let a_l_length = a_vec0.length();
        if (a_l_length - d_par).abs() <= PCONFUSION {
            let a_dir_l = a_vec0 / a_l_length;
            let a_pl = a_p2d[0] - *the_params.value(1) * a_dir_l;
            return Some(Curve2d::Line(Line2d::new(a_pl, a_dir_l)));
        }

        // Create straight bspline
        // OCCT L1089-1101: poles (1,2), knots (1,2), mults (2,2), degree 1.
        // Architecture bridge: OCCT passes distinct knots + multiplicities;
        // the rcad BSplineCurve2 stores the flat expanded knot vector
        // [k1, k1, k2, k2].
        let a_k1 = *the_params.value(1);
        let a_k2 = *the_params.value(the_params.length());
        let a_c2d = Curve2d::BSpline(BSplineCurve2 {
            degree: 1,
            knots: vec![a_k1, a_k1, a_k2, a_k2],
            control_points: vec![a_p2d[0], a_p2d[3]],
            weights: vec![1.0, 1.0],
        });
        Some(a_c2d)
    }

    // ------------------------------------------------------------------
    // OCCT .cxx L1104-2785 (approxPCurve and the protected members) —
    // translated in project_curve_on_surface_approx.rs.
    // ------------------------------------------------------------------
}
