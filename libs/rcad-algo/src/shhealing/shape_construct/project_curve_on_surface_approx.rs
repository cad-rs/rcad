//! OCCT ShapeConstruct_ProjectCurveOnSurface.cxx L1104-2785 — the protected
//! approximation members, continued from `project_curve_on_surface.rs`:
//! - `approxPCurve` (L1106-1906)
//! - `correctExtremity` (L1910-2035)
//! - `insertAdditionalPointOrAdjust` (L2039-2148)
//! - `interpolatePCurve` (L2152-2216)
//! - `approximatePCurve` (L2220-2287)
//! - `checkPoints` (L2291-2369)
//! - `checkPoints2d` (L2373-2457)
//! - `isAnIsoparametric` (L2461-2785)

use glam::{DVec2, DVec3};
use rcad_kernel::base::geom_api::geom2d_interpolate::Geom2dInterpolate;
use rcad_kernel::geom::{BSplineCurve2, Curve2d, Curve3, CurveEval, Surface3, SurfaceEval};
use rcad_kernel::math::el::in_period;
use rcad_kernel::math::gp::GP_RESOLUTION;
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::precision::{is_infinite_value, CONFUSION, PCONFUSION};

use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};

use super::array1::Array1;
use super::gap_deps::{
    int_curve_int_conic_conic_lin_lin, shape_analysis_adjust_to_period,
    shape_analysis_curve_next_project, shape_analysis_curve_project, ShapeAnalysisSurface,
};
use super::points_to_bspline::PointsToBSpline;
use super::project_curve_on_surface::is_spherical_surface;
use super::project_curve_on_surface::{ArrayOfPnt, ArrayOfPnt2d, ArrayOfReal, ProjectCurveOnSurface};
use super::project_curve_on_surface_ns::{
    adjust_second_to_first_point, p2d_coord, p2d_set_coord, project_degenerated_points, CachePoint,
};

/// OCCT `RealLast()` (Standard_Real.hxx).
const REAL_LAST: f64 = f64::MAX;

/// 1-based coordinate setter for the `IndCoord` (1=U, 2=V) OCCT idiom.
fn set_coord(the_p: &mut DVec2, the_ind_coord: usize, the_val: f64) {
    p2d_set_coord(the_p, the_ind_coord, the_val)
}

impl ProjectCurveOnSurface {
    /// OCCT .cxx L1106-1906: approxPCurve(theNbPnt, theC3D, theTolFirst,
    /// theTolLast, thePoints, theParams, thePoints2d, theC2D) — the main
    /// approximation routine for pcurve computation.
    pub(crate) fn approx_p_curve(
        &mut self,
        the_nb_pnt: i32,
        the_c3d: &Curve3,
        the_tol_first: f64,
        the_tol_last: f64,
        the_points: &mut ArrayOfPnt,
        the_params: &mut ArrayOfReal,
        the_points2d: &mut ArrayOfPnt2d,
        the_c2d: &mut Option<Curve2d>,
    ) -> bool {
        let a_surf_sa: &ShapeAnalysisSurface = self.my_surf.as_ref().unwrap();

        // For performance, first try to handle typical case when pcurve is straight
        let mut is_recompute = false;
        let mut is_from_cache_line = false;

        *the_c2d = self.get_line(
            the_points,
            the_params,
            the_points2d,
            self.my_preci,
            &mut is_recompute,
            &mut is_from_cache_line,
        );
        if the_c2d.is_some() {
            // Fill cache
            let mut a_change_cycle = false;
            if !self.my_cache.is_empty()
                && self.my_cache.value(0).first.distance(*the_points.value(1))
                    > self
                        .my_cache
                        .value(0)
                        .first
                        .distance(*the_points.value(the_nb_pnt as usize))
                && self
                    .my_cache
                    .value(0)
                    .first
                    .distance(*the_points.value(the_nb_pnt as usize))
                    < CONFUSION
            {
                a_change_cycle = true;
            }

            self.my_cache.resize(
                0,
                1,
                CachePoint {
                    first: DVec3::ZERO,
                    second: DVec2::ZERO,
                },
            );
            if a_change_cycle {
                *self.my_cache.change_value(0) = CachePoint {
                    first: *the_points.value(1),
                    second: *the_points2d.value(1),
                };
                *self.my_cache.change_value(1) = CachePoint {
                    first: *the_points.value(the_nb_pnt as usize),
                    second: *the_points2d.value(the_nb_pnt as usize),
                };
            } else {
                *self.my_cache.change_value(0) = CachePoint {
                    first: *the_points.value(the_nb_pnt as usize),
                    second: *the_points2d.value(the_nb_pnt as usize),
                };
                *self.my_cache.change_value(1) = CachePoint {
                    first: *the_points.value(1),
                    second: *the_points2d.value(1),
                };
            }
            return true;
        }

        let mut is_done = true;

        // Test if the curve 3d is a boundary of the surface
        let mut iso_type_u = false;
        let mut iso_par2d3d = false;
        let mut p1_on_iso = false;
        let mut p2_on_iso = false;
        let mut value_p1 = DVec2::ZERO;
        let mut value_p2 = DVec2::ZERO;
        let mut c_iso: Option<Curve3> = None;
        let mut t1 = 0.0f64;
        let mut t2 = 0.0f64;
        let mut iso_closed: bool;

        // OCCT L1153-1158: sType = mySurf->Surface()->DynamicType();
        // isAnalytic = true, set false for Bezier/BSpline.
        let is_analytic = !matches!(
            a_surf_sa.surface(),
            Surface3::Bezier(_) | Surface3::BSpline(_)
        );

        // OCCT L1160-1161: mySurf->Surface()->Bounds(uf, ul, vf, vl).
        let [uf, ul, vf, vl] = a_surf_sa.bounds();
        iso_closed = false;

        let mut pout = ArrayOfReal::new(1, the_nb_pnt as usize, 0.0);
        for i in 1..=the_nb_pnt {
            pout.set_value(i as usize, 0.0);
        }

        let iso_param = self.is_an_isoparametric(
            the_nb_pnt,
            the_points,
            the_params,
            &mut iso_type_u,
            &mut p1_on_iso,
            &mut value_p1,
            &mut p2_on_iso,
            &mut value_p2,
            &mut iso_par2d3d,
            &mut c_iso,
            &mut t1,
            &mut t2,
            &mut pout,
        );

        // Projection of the points on surfaces
        let mut p3d;
        let mut p2d = DVec2::ZERO;
        let mut iso_value = 0.0f64;
        let mut iso_par1 = 0.0f64;
        let mut iso_par2 = 0.0f64;
        let mut t_par = 0.0f64;

        if iso_param {
            let parf;
            let parl;
            if iso_type_u {
                iso_value = value_p1.x;
                iso_par1 = value_p1.y;
                iso_par2 = value_p2.y;
                iso_closed = a_surf_sa.is_v_closed(self.my_preci);
                parf = vf;
                parl = vl;
            } else {
                iso_value = value_p1.y;
                iso_par1 = value_p1.x;
                iso_par2 = value_p2.x;
                iso_closed = a_surf_sa.is_u_closed(self.my_preci);
                parf = uf;
                parl = ul;
            }
            if !iso_par2d3d && !is_analytic {
                // OCCT L1212-1213: Cf = cIso->FirstParameter(); Cl = cIso->LastParameter().
                let ci = c_iso.as_ref().unwrap();
                let mut cf = ci.default_domain()[0];
                let mut cl = ci.default_domain()[1];
                if is_infinite_value(cf) {
                    cf = -1000.0;
                }
                if is_infinite_value(cl) {
                    cl = 1000.0;
                }
                let _ = (&mut cf, &mut cl);

                let tdeb = *pout.value(2);

                if iso_closed && (iso_par1 == parf || iso_par1 == parl) {
                    if (tdeb - parf).abs() < (tdeb - parl).abs() {
                        iso_par1 = parf;
                    } else {
                        iso_par1 = parl;
                    }
                    if iso_type_u {
                        value_p1.y = iso_par1;
                    } else {
                        value_p1.x = iso_par1;
                    }
                }
                if iso_closed && (iso_par2 == parf || iso_par2 == parl) {
                    let tfin = *pout.value((the_nb_pnt - 1) as usize);
                    if (tfin - parf).abs() < (tfin - parl).abs() {
                        iso_par2 = parf;
                    } else {
                        iso_par2 = parl;
                    }
                    if iso_type_u {
                        value_p2.y = iso_par2;
                    } else {
                        value_p2.x = iso_par2;
                    }
                }

                if !iso_closed {
                    if ((tdeb - iso_par1).abs() > (tdeb - iso_par2).abs())
                        && ((*pout.value((the_nb_pnt - 1) as usize) - iso_par2).abs()
                            > (*pout.value((the_nb_pnt - 1) as usize) - iso_par1).abs())
                    {
                        let value_tmp = value_p1;
                        value_p1 = value_p2;
                        value_p2 = value_tmp;
                        if iso_type_u {
                            iso_value = value_p1.x;
                            iso_par1 = value_p1.y;
                            iso_par2 = value_p2.y;
                        } else {
                            iso_value = value_p1.y;
                            iso_par1 = value_p1.x;
                            iso_par2 = value_p2.x;
                        }
                    }
                }
            }
        }

        let up = ul - uf;
        let vp = vl - vf;
        let mut gap = self.my_preci;
        let mut a_change_cycle = false;

        let mut is_from_cache = false;
        let mut a_saved_point = DVec2::ZERO;

        if !self.my_cache.is_empty()
            && self.my_cache.value(0).first.distance(*the_points.value(1))
                > self
                    .my_cache
                    .value(0)
                    .first
                    .distance(*the_points.value(the_nb_pnt as usize))
        {
            if self
                .my_cache
                .value(0)
                .first
                .distance(*the_points.value(the_nb_pnt as usize))
                < CONFUSION
            {
                a_change_cycle = true;
            }
        }

        let mut need_resolve_u_jump = false;
        let mut need_resolve_v_jump = false;
        let mut prev_p3d = DVec3::ZERO;
        let mut prev_p2d = DVec2::ZERO;

        for ii in 1..=the_nb_pnt {
            let a_pnt_index = if a_change_cycle {
                the_nb_pnt - ii + 1
            } else {
                ii
            };
            p3d = *the_points.value(a_pnt_index as usize);

            if iso_param {
                if iso_par2d3d {
                    if iso_par2 > iso_par1 {
                        t_par = *the_params.value(a_pnt_index as usize);
                    } else {
                        t_par = t1 + t2 - *the_params.value(a_pnt_index as usize);
                    }
                } else if !is_analytic {
                    if a_pnt_index == 1 {
                        t_par = iso_par1;
                    } else if a_pnt_index == the_nb_pnt {
                        t_par = iso_par2;
                    } else {
                        t_par = *pout.value(a_pnt_index as usize);
                    }
                }

                if !iso_par2d3d && is_analytic {
                    if a_pnt_index == 1 {
                        p2d = value_p1;
                    } else if a_pnt_index == the_nb_pnt {
                        p2d = value_p2;
                    } else {
                        p2d = a_surf_sa.next_value_of_uv(
                            p2d,
                            p3d,
                            self.my_preci,
                            CONFUSION + 1000.0 * gap,
                        );
                        gap = a_surf_sa.gap();
                    }
                } else {
                    if iso_type_u {
                        p2d = DVec2::new(iso_value, t_par);
                    } else {
                        p2d = DVec2::new(t_par, iso_value);
                    }
                }
            } else {
                if a_pnt_index == 1 && p1_on_iso {
                    p2d = value_p1;
                } else if a_pnt_index == the_nb_pnt && p2_on_iso {
                    p2d = value_p2;
                } else {
                    if a_pnt_index == 1 || a_pnt_index == the_nb_pnt {
                        if !is_recompute {
                            p2d = *the_points2d.value(a_pnt_index as usize);
                            gap = a_surf_sa.gap();
                            if a_pnt_index == 1 {
                                is_from_cache = is_from_cache_line;
                                a_saved_point = p2d;
                            }
                            continue;
                        } else {
                            let mut j = 0usize;
                            let a_nb_cache = self.my_cache.length();
                            while j < a_nb_cache {
                                let a_cache_pnt = self.my_cache.value(j);
                                if a_cache_pnt.first.distance_squared(p3d)
                                    < self.my_preci * self.my_preci
                                {
                                    p2d = a_surf_sa.next_value_of_uv(
                                        a_cache_pnt.second,
                                        p3d,
                                        self.my_preci,
                                        CONFUSION + gap,
                                    );
                                    if a_pnt_index == 1 {
                                        is_from_cache = true;
                                        a_saved_point = a_cache_pnt.second;
                                    }
                                    break;
                                }
                                j += 1;
                            }
                            if j >= a_nb_cache {
                                p2d = a_surf_sa.value_of_uv(p3d, self.my_preci);
                            }
                        }
                    } else {
                        p2d = a_surf_sa.next_value_of_uv(
                            p2d,
                            p3d,
                            self.my_preci,
                            CONFUSION + 1000.0 * gap,
                        );
                    }
                    gap = a_surf_sa.gap();
                }
            }
            the_points2d.set_value(a_pnt_index as usize, p2d);

            if the_nb_pnt > 23 && ii > 2 && ii < the_nb_pnt {
                if (p2d.x - prev_p2d.x).abs() > 0.95 * up
                    && prev_p3d.distance(p3d) < self.my_preci
                    && !a_surf_sa.is_u_closed(self.my_preci)
                    && a_surf_sa.nb_singularities(self.my_preci) > 0
                    && matches!(a_surf_sa.surface(), Surface3::BSpline(_))
                {
                    need_resolve_u_jump = true;
                }
                if (p2d.y - prev_p2d.y).abs() > 0.95 * vp
                    && prev_p3d.distance(p3d) < self.my_preci
                    && !a_surf_sa.is_v_closed(self.my_preci)
                    && a_surf_sa.nb_singularities(self.my_preci) > 0
                    && matches!(a_surf_sa.surface(), Surface3::BSpline(_))
                {
                    need_resolve_v_jump = true;
                }
            }
            prev_p3d = p3d;
            prev_p2d = p2d;
            if ii > 1 {
                if a_change_cycle {
                    let next = *the_points2d.value((a_pnt_index + 1) as usize);
                    p2d = DVec2::new(2.0 * p2d.x - next.x, 2.0 * p2d.y - next.y);
                } else {
                    let prev_stored = *the_points2d.value((a_pnt_index - 1) as usize);
                    p2d = DVec2::new(2.0 * p2d.x - prev_stored.x, 2.0 * p2d.y - prev_stored.y);
                }
            }
        }

        if !iso_par2d3d {
            project_degenerated_points(
                a_surf_sa,
                the_nb_pnt,
                the_points,
                the_points2d,
                self.my_preci,
                true,
            );
            project_degenerated_points(
                a_surf_sa,
                the_nb_pnt,
                the_points,
                the_points2d,
                self.my_preci,
                false,
            );
        }

        // Check extremities for singularities
        let a_point_first = *the_points.first();
        let a_point_last = *the_points.last();
        let a_tol_first = if the_tol_first < 0.0 {
            CONFUSION
        } else {
            the_tol_first
        };
        let a_tol_last = if the_tol_last < 0.0 {
            CONFUSION
        } else {
            the_tol_last
        };

        let mut i = 1i32;
        loop {
            let mut a_preci = 0.0f64;
            let mut a_first_par = 0.0f64;
            let mut a_last_par = 0.0f64;
            let mut a_p3d = DVec3::ZERO;
            let mut a_first_p2d = DVec2::ZERO;
            let mut a_last_p2d = DVec2::ZERO;
            let mut is_u_iso = false;
            if !a_surf_sa.singularity(
                i,
                &mut a_preci,
                &mut a_p3d,
                &mut a_first_p2d,
                &mut a_last_p2d,
                &mut a_first_par,
                &mut a_last_par,
                &mut is_u_iso,
            ) {
                break;
            }
            if a_preci <= CONFUSION && a_point_first.distance(a_p3d) <= a_tol_first {
                self.correct_extremity(
                    the_c3d,
                    the_params,
                    the_points2d,
                    true,
                    a_first_p2d,
                    is_u_iso,
                );
            }
            if a_preci <= CONFUSION && a_point_last.distance(a_p3d) <= a_tol_last {
                self.correct_extremity(
                    the_c3d,
                    the_params,
                    the_points2d,
                    false,
                    a_first_p2d,
                    is_u_iso,
                );
            }
            i += 1;
        }

        // Handle U-closed surfaces
        let tol_on_u_period = CONFUSION * up;
        let tol_on_v_period = CONFUSION * vp;

        if a_surf_sa.is_u_closed(self.my_preci) || need_resolve_u_jump {
            let mut first_x = the_points2d.value(1).x;
            if !is_from_cache {
                while first_x < uf {
                    first_x += up;
                    the_points2d.change_value(1).x = first_x;
                }
                while first_x > ul {
                    first_x -= up;
                    the_points2d.change_value(1).x = first_x;
                }
            }

            if a_surf_sa.surface().is_u_periodic() && is_from_cache {
                let mut a_min_param = uf;
                let mut a_max_param = ul;
                while a_min_param > a_saved_point.x {
                    a_min_param -= up;
                    a_max_param -= up;
                }
                while a_max_param < a_saved_point.x {
                    a_min_param += up;
                    a_max_param += up;
                }
                let a_shift = shape_analysis_adjust_to_period(first_x, a_min_param, a_max_param);
                first_x += a_shift;
                the_points2d.change_value(1).x = first_x;
            }
            let mut prev_x = first_x;

            let mut min_x = first_x;
            let mut max_x = first_x;
            let mut to_adjust = false;

            let mut a_pnt_iter = 2usize;
            while a_pnt_iter <= the_points2d.length() {
                let mut cur_x = the_points2d.value(a_pnt_iter).x;
                let dist2d = (cur_x - prev_x).abs();
                if dist2d > up / 2.0 {
                    self.insert_additional_point_or_adjust(
                        &mut to_adjust,
                        1,
                        up,
                        tol_on_u_period,
                        &mut cur_x,
                        prev_x,
                        the_c3d,
                        &mut a_pnt_iter,
                        the_points,
                        the_params,
                        the_points2d,
                    );
                }
                prev_x = cur_x;
                if min_x > cur_x {
                    min_x = cur_x;
                } else if max_x < cur_x {
                    max_x = cur_x;
                }
                a_pnt_iter += 1;
            }

            if !is_from_cache {
                let mid_x = 0.5 * (min_x + max_x);
                let mut shift_x = 0.0f64;
                if mid_x > ul {
                    shift_x = -up;
                } else if mid_x < uf {
                    shift_x = up;
                }
                if shift_x != 0.0 {
                    for a_pnt_iter in 1..=the_points2d.length() {
                        let v = the_points2d.value(a_pnt_iter).x + shift_x;
                        the_points2d.change_value(a_pnt_iter).x = v;
                    }
                }
            }
        }

        // Handle V-closed surfaces
        if a_surf_sa.is_v_closed(self.my_preci)
            || need_resolve_v_jump
            || is_spherical_surface(a_surf_sa.surface())
        {
            let mut first_y = the_points2d.value(1).y;
            if !is_from_cache {
                while first_y < vf {
                    first_y += vp;
                    the_points2d.change_value(1).y = first_y;
                }
                while first_y > vl {
                    first_y -= vp;
                    the_points2d.change_value(1).y = first_y;
                }
            }

            if a_surf_sa.surface().is_v_periodic() && is_from_cache {
                let mut a_min_param = vf;
                let mut a_max_param = vl;
                while a_min_param > a_saved_point.y {
                    a_min_param -= vp;
                    a_max_param -= vp;
                }
                while a_max_param < a_saved_point.y {
                    a_min_param += vp;
                    a_max_param += vp;
                }
                let a_shift = shape_analysis_adjust_to_period(first_y, a_min_param, a_max_param);
                first_y += a_shift;
                the_points2d.change_value(1).y = first_y;
            }
            let mut prev_y = first_y;

            let mut min_y = first_y;
            let mut max_y = first_y;
            let mut to_adjust = false;

            let mut a_pnt_iter = 2usize;
            while a_pnt_iter <= the_points2d.length() {
                let mut cur_y = the_points2d.value(a_pnt_iter).y;
                let dist2d = (cur_y - prev_y).abs();
                if dist2d > vp / 2.0 {
                    self.insert_additional_point_or_adjust(
                        &mut to_adjust,
                        2,
                        vp,
                        tol_on_v_period,
                        &mut cur_y,
                        prev_y,
                        the_c3d,
                        &mut a_pnt_iter,
                        the_points,
                        the_params,
                        the_points2d,
                    );
                }
                prev_y = cur_y;
                if min_y > cur_y {
                    min_y = cur_y;
                } else if max_y < cur_y {
                    max_y = cur_y;
                }
                a_pnt_iter += 1;
            }

            if !is_from_cache {
                let mid_y = 0.5 * (min_y + max_y);
                let mut shift_y = 0.0f64;
                if mid_y > vl {
                    shift_y = -vp;
                } else if mid_y < vf {
                    shift_y = vp;
                }
                if shift_y != 0.0 {
                    for a_pnt_iter in 1..=the_points2d.length() {
                        let v = the_points2d.value(a_pnt_iter).y + shift_y;
                        the_points2d.change_value(a_pnt_iter).y = v;
                    }
                }
            }
        }

        // Handle V-closed seam adjustment
        if a_surf_sa.is_v_closed(self.my_preci) || is_spherical_surface(a_surf_sa.surface()) {
            for a_pnt_iter in 2..=the_points2d.length() {
                let dist2d =
                    (the_points2d.value(a_pnt_iter).y - the_points2d.value(a_pnt_iter - 1).y)
                        .abs();
                if dist2d > vp / 2.0 {
                    let mut prev_on_first;
                    let mut prev_on_last;
                    let mut curr_on_first;
                    let mut curr_on_last;

                    let dist_prev_vf = (the_points2d.value(a_pnt_iter - 1).y - vf).abs();
                    let dist_prev_vl = (the_points2d.value(a_pnt_iter - 1).y - vl).abs();
                    let dist_curr_vf = (the_points2d.value(a_pnt_iter).y - vf).abs();
                    let dist_curr_vl = (the_points2d.value(a_pnt_iter).y - vl).abs();

                    let mut the_min = dist_prev_vf;
                    prev_on_first = true;
                    prev_on_last = false;
                    curr_on_first = false;
                    curr_on_last = false;
                    if dist_prev_vl < the_min {
                        the_min = dist_prev_vl;
                        prev_on_first = false;
                        prev_on_last = true;
                    }
                    if dist_curr_vf < the_min {
                        the_min = dist_curr_vf;
                        prev_on_first = false;
                        prev_on_last = false;
                        curr_on_first = true;
                    }
                    if dist_curr_vl < the_min {
                        // OCCT L1719-1725: theMin is not updated here (kept
                        // verbatim).
                        prev_on_first = false;
                        prev_on_last = false;
                        curr_on_first = false;
                        curr_on_last = true;
                    }

                    if prev_on_first {
                        the_points2d.change_value(a_pnt_iter - 1).y = vf;
                    } else if prev_on_last {
                        the_points2d.change_value(a_pnt_iter - 1).y = vl;
                    } else if curr_on_first {
                        the_points2d.change_value(a_pnt_iter).y = vf;
                    } else if curr_on_last {
                        the_points2d.change_value(a_pnt_iter).y = vl;
                    }
                }
            }
        }

        // Handle AdjustOverDegen
        if self.my_adjust_over_degen != -1 {
            if a_surf_sa.is_u_closed(self.my_preci) {
                // OCCT calls IsDegenerated for effect (the result is unused).
                let _ = a_surf_sa.is_degenerated(DVec3::new(0.0, 0.0, 0.0), self.my_preci);
                if a_surf_sa.nb_singularities(self.my_preci) > 0 {
                    let mut prev_x = 0.0f64;
                    let mut on_bound = 0i32;
                    let mut prev_on_bound = 0i32;
                    let mut ind = 1usize;
                    let mut start = true;
                    while ind <= the_points2d.length() {
                        let cur_x = the_points2d.value(ind).x;
                        if a_surf_sa.is_degenerated(*the_points.value(ind), CONFUSION) {
                            ind += 1;
                            continue;
                        }
                        on_bound = (((cur_x - 0.5 * (ul + uf)).abs() - up / 2.0).abs()
                            <= PCONFUSION) as i32;
                        if !start && ((cur_x - prev_x).abs() - up / 2.0).abs() <= 0.01 * up {
                            break;
                        }
                        start = false;
                        prev_x = cur_x;
                        prev_on_bound = on_bound;
                        ind += 1;
                    }
                    if ind <= the_points2d.length() {
                        let prev_x2 = if self.my_adjust_over_degen != 0 { uf } else { ul };
                        let d_u = up / 2.0 + PCONFUSION;
                        if prev_on_bound != 0 {
                            the_points2d.change_value(ind - 1).x = prev_x2;
                            let mut j = ind as i64 - 2;
                            while j > 0 {
                                let mut cur_x = the_points2d.value(j as usize).x;
                                while cur_x < prev_x2 - d_u {
                                    cur_x += up;
                                    the_points2d.change_value(j as usize).x = cur_x;
                                }
                                while cur_x > prev_x2 + d_u {
                                    cur_x -= up;
                                    the_points2d.change_value(j as usize).x = cur_x;
                                }
                                j -= 1;
                            }
                        } else if on_bound != 0 {
                            the_points2d.change_value(ind).x = prev_x2;
                            let mut j = ind + 1;
                            while j <= the_points2d.length() {
                                let mut cur_x = the_points2d.value(j).x;
                                while cur_x < prev_x2 - d_u {
                                    cur_x += up;
                                    the_points2d.change_value(j).x = cur_x;
                                }
                                while cur_x > prev_x2 + d_u {
                                    cur_x -= up;
                                    the_points2d.change_value(j).x = cur_x;
                                }
                                j += 1;
                            }
                        }
                        self.my_status |= encode_status(ShapeExtendStatus::Done4);
                    }
                }
            } else if a_surf_sa.is_v_closed(self.my_preci) {
                // OCCT calls IsDegenerated for effect (the result is unused).
                let _ = a_surf_sa.is_degenerated(DVec3::new(0.0, 0.0, 0.0), self.my_preci);
                if a_surf_sa.nb_singularities(self.my_preci) > 0 {
                    let mut prev_y = 0.0f64;
                    let mut on_bound = 0i32;
                    let mut prev_on_bound = 0i32;
                    let mut ind = 1usize;
                    let mut start = true;
                    while ind <= the_points2d.length() {
                        let cur_y = the_points2d.value(ind).y;
                        if a_surf_sa.is_degenerated(*the_points.value(ind), CONFUSION) {
                            ind += 1;
                            continue;
                        }
                        on_bound = (((cur_y - 0.5 * (vl + vf)).abs() - vp / 2.0).abs()
                            <= PCONFUSION) as i32;
                        if !start && ((cur_y - prev_y).abs() - vp / 2.0).abs() <= 0.01 * vp {
                            break;
                        }
                        start = false;
                        prev_y = cur_y;
                        prev_on_bound = on_bound;
                        ind += 1;
                    }
                    if ind <= the_points2d.length() {
                        let prev_y2 = if self.my_adjust_over_degen != 0 { vf } else { vl };
                        let d_v = vp / 2.0 + PCONFUSION;
                        if prev_on_bound != 0 {
                            the_points2d.change_value(ind - 1).y = prev_y2;
                            let mut j = ind as i64 - 2;
                            while j > 0 {
                                let mut cur_y = the_points2d.value(j as usize).y;
                                while cur_y < prev_y2 - d_v {
                                    cur_y += vp;
                                    the_points2d.change_value(j as usize).y = cur_y;
                                }
                                while cur_y > prev_y2 + d_v {
                                    cur_y -= vp;
                                    the_points2d.change_value(j as usize).y = cur_y;
                                }
                                j -= 1;
                            }
                        } else if on_bound != 0 {
                            the_points2d.change_value(ind).y = prev_y2;
                            let mut j = ind + 1;
                            while j <= the_points2d.length() {
                                let mut cur_y = the_points2d.value(j).y;
                                while cur_y < prev_y2 - d_v {
                                    cur_y += vp;
                                    the_points2d.change_value(j).y = cur_y;
                                }
                                while cur_y > prev_y2 + d_v {
                                    cur_y -= vp;
                                    the_points2d.change_value(j).y = cur_y;
                                }
                                j += 1;
                            }
                        }
                        self.my_status |= encode_status(ShapeExtendStatus::Done4);
                    }
                }
            }
        }

        // Fill cache
        self.my_cache.resize(
            0,
            1,
            CachePoint {
                first: DVec3::ZERO,
                second: DVec2::ZERO,
            },
        );
        if a_change_cycle {
            *self.my_cache.change_value(0) = CachePoint {
                first: *the_points.value(1),
                second: *the_points2d.value(1),
            };
            *self.my_cache.change_value(1) = CachePoint {
                first: *the_points.last(),
                second: *the_points2d.last(),
            };
        } else {
            *self.my_cache.change_value(0) = CachePoint {
                first: *the_points.last(),
                second: *the_points2d.last(),
            };
            *self.my_cache.change_value(1) = CachePoint {
                first: *the_points.value(1),
                second: *the_points2d.value(1),
            };
        }

        is_done
    }

    /// OCCT .cxx L1910-2035: correctExtremity(theC3D, theParams, thePoints2d,
    /// theIsFirstPoint, thePointOnIsoLine, theIsUIso) — corrects the extremity
    /// point near a singularity.
    pub(crate) fn correct_extremity(
        &self,
        the_c3d: &Curve3,
        the_params: &ArrayOfReal,
        the_points2d: &mut ArrayOfPnt2d,
        the_is_first_point: bool,
        the_point_on_iso_line: DVec2,
        the_is_u_iso: bool,
    ) {
        let a_surf_sa = self.my_surf.as_ref().unwrap();
        let nb_pnt = the_points2d.length();
        let ind_coord = if the_is_u_iso { 2usize } else { 1usize };
        let singularity_coord = p2d_coord(the_point_on_iso_line, 3 - ind_coord);
        let end_point = if the_is_first_point {
            *the_points2d.value(1)
        } else {
            *the_points2d.value(nb_pnt)
        };
        let finish_coord = p2d_coord(end_point, 3 - ind_coord);

        let mut a_dir = if the_is_u_iso { DVec2::Y } else { DVec2::X };
        let an_iso_line = rcad_kernel::geom::Line2d::new(end_point, a_dir);

        // OCCT IntRes2d_Domain Dom1, Dom2 — default-constructed (unbounded);
        // the rcad line-line intersection bridge carries no domain argument
        // (arch. diff. note on `int_curve_int_conic_conic_lin_lin`).

        let is_periodic = if the_is_u_iso {
            a_surf_sa.surface().is_v_periodic()
        } else {
            a_surf_sa.surface().is_u_periodic()
        };

        let mut first_point_of_line;
        let mut second_point_of_line;
        let mut finish_param;
        let mut first_param;
        let mut second_param;

        if the_is_first_point {
            first_point_of_line = *the_points2d.value(3);
            second_point_of_line = *the_points2d.value(2);
            finish_param = *the_params.value(1);
            first_param = *the_params.value(3);
            second_param = *the_params.value(2);
        } else {
            first_point_of_line = *the_points2d.value(nb_pnt - 2);
            second_point_of_line = *the_points2d.value(nb_pnt - 1);
            finish_param = *the_params.value(nb_pnt);
            first_param = *the_params.value(nb_pnt - 2);
            second_param = *the_params.value(nb_pnt - 1);
        }

        if singularity_coord > finish_coord
            && p2d_coord(second_point_of_line, 3 - ind_coord) > finish_coord
        {
            return;
        }
        if singularity_coord < finish_coord
            && p2d_coord(second_point_of_line, 3 - ind_coord) < finish_coord
        {
            return;
        }

        {
            let a_prev_dist =
                (p2d_coord(second_point_of_line, ind_coord) - p2d_coord(first_point_of_line, ind_coord))
                    .abs();
            let a_cur_dist =
                (p2d_coord(end_point, ind_coord) - p2d_coord(second_point_of_line, ind_coord)).abs();
            if a_cur_dist <= 2.0 * a_prev_dist {
                return;
            }
        }

        let mut finish_point = if the_is_u_iso {
            DVec2::new(finish_coord, second_point_of_line.y)
        } else {
            DVec2::new(second_point_of_line.x, finish_coord)
        };

        loop {
            if (p2d_coord(second_point_of_line, 3 - ind_coord) - finish_coord).abs()
                <= 2.0 * PCONFUSION
            {
                break;
            }

            let a_vec = second_point_of_line - first_point_of_line;
            let a_sq_magnitude = a_vec.length_squared();
            if a_sq_magnitude <= 1.0e-32 {
                break;
            }
            // OCCT aDir.SetCoord(aVec.X(), aVec.Y()) — gp_Dir2d keeps the unit
            // invariant; Line2d::new normalizes the same way.
            a_dir = a_vec;

            let a_line = rcad_kernel::geom::Line2d::new(first_point_of_line, a_dir);
            if let Some(int_value) = int_curve_int_conic_conic_lin_lin(&an_iso_line, &a_line) {
                // OCCT IntRes2d_IntersectionPoint IntPoint = Intersector.Point(1);
                // FinishPoint = IntPoint.Value();
                finish_point = int_value;
            } else {
                finish_point = if the_is_u_iso {
                    DVec2::new(finish_coord, second_point_of_line.y)
                } else {
                    DVec2::new(second_point_of_line.x, finish_coord)
                };
            }

            let prev_point = first_point_of_line;
            first_point_of_line = second_point_of_line;
            first_param = second_param;
            second_param = (first_param + finish_param) / 2.0;
            if (second_param - first_param).abs() <= 2.0 * PCONFUSION {
                break;
            }
            let a_p3d = the_c3d.point_at(second_param);
            second_point_of_line =
                a_surf_sa.next_value_of_uv(first_point_of_line, a_p3d, self.my_preci, CONFUSION);
            if is_periodic {
                adjust_second_to_first_point(
                    first_point_of_line,
                    &mut second_point_of_line,
                    a_surf_sa.surface(),
                );
            }

            let a_prev_dist =
                (p2d_coord(first_point_of_line, ind_coord) - p2d_coord(prev_point, ind_coord))
                    .abs();
            let a_cur_dist = (p2d_coord(second_point_of_line, ind_coord)
                - p2d_coord(first_point_of_line, ind_coord))
            .abs();
            if a_cur_dist > 2.0 * a_prev_dist {
                break;
            }
        }

        if the_is_first_point {
            *the_points2d.change_value(1) = finish_point;
        } else {
            *the_points2d.change_value(nb_pnt) = finish_point;
        }
    }

    /// OCCT .cxx L2039-2148: insertAdditionalPointOrAdjust(theToAdjust,
    /// theIndCoord, thePeriod, theTolOnPeriod, theCurCoord, thePrevCoord,
    /// theC3D, theIndex, thePoints, theParams, thePoints2d) — inserts an
    /// additional point or adjusts the coordinate to handle period jumps.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn insert_additional_point_or_adjust(
        &self,
        the_to_adjust: &mut bool,
        the_ind_coord: usize,
        the_period: f64,
        the_tol_on_period: f64,
        the_cur_coord: &mut f64,
        the_prev_coord: f64,
        the_c3d: &Curve3,
        the_index: &mut usize,
        the_points: &mut ArrayOfPnt,
        the_params: &mut ArrayOfReal,
        the_points2d: &mut ArrayOfPnt2d,
    ) {
        let a_surf_sa = self.my_surf.as_ref().unwrap();
        let corrected_cur_coord = in_period(
            *the_cur_coord,
            the_prev_coord - the_period / 2.0,
            the_prev_coord + the_period / 2.0,
        );

        if !*the_to_adjust {
            let cur_par = *the_params.value(*the_index);
            let prev_par = *the_params.value(*the_index - 1);
            let mut mid_par = (prev_par + cur_par) / 2.0;
            let mut mid_p3d = the_c3d.point_at(mid_par);
            let mut mid_p2d = a_surf_sa.value_of_uv(mid_p3d, self.my_preci);
            let mut mid_coord = p2d_coord(mid_p2d, the_ind_coord);
            mid_coord = in_period(
                mid_coord,
                the_prev_coord - the_period / 2.0,
                the_prev_coord + the_period / 2.0,
            );
            let mut first_coord = the_prev_coord;
            let mut last_coord = corrected_cur_coord;
            if last_coord < first_coord {
                let tmp = first_coord;
                first_coord = last_coord;
                last_coord = tmp;
            }
            if last_coord - first_coord <= the_tol_on_period {
                *the_to_adjust = true;
            } else if first_coord <= mid_coord && mid_coord <= last_coord {
                *the_to_adjust = true;
            } else {
                let mut success = true;
                let mut first_t = prev_par;
                let mut last_t = cur_par;
                mid_coord = p2d_coord(mid_p2d, the_ind_coord);
                while (mid_coord - the_prev_coord).abs() >= the_period / 2.0 - the_tol_on_period
                    || (*the_cur_coord - mid_coord).abs() >= the_period / 2.0 - the_tol_on_period
                {
                    if mid_par - first_t <= PCONFUSION || last_t - mid_par <= PCONFUSION {
                        success = false;
                        break;
                    }
                    if (mid_coord - the_prev_coord).abs() >= the_period / 2.0 - the_tol_on_period
                    {
                        last_t = (first_t + last_t) / 2.0;
                    } else {
                        first_t = (first_t + last_t) / 2.0;
                    }
                    mid_par = (first_t + last_t) / 2.0;
                    mid_p3d = the_c3d.point_at(mid_par);
                    mid_p2d = a_surf_sa.value_of_uv(mid_p3d, self.my_preci);
                    mid_coord = p2d_coord(mid_p2d, the_ind_coord);
                }
                if success {
                    // Insert additional point - need to resize arrays
                    let a_new_length = the_points.length() + 1;
                    let mut a_new_points = ArrayOfPnt::new(1, a_new_length, DVec3::ZERO);
                    let mut a_new_params = ArrayOfReal::new(1, a_new_length, 0.0);
                    let mut a_new_points2d = ArrayOfPnt2d::new(1, a_new_length, DVec2::ZERO);

                    for i in 1..*the_index {
                        a_new_points.set_value(i, *the_points.value(i));
                        a_new_params.set_value(i, *the_params.value(i));
                        a_new_points2d.set_value(i, *the_points2d.value(i));
                    }
                    a_new_points.set_value(*the_index, mid_p3d);
                    a_new_params.set_value(*the_index, mid_par);
                    a_new_points2d.set_value(*the_index, mid_p2d);
                    for i in *the_index..=the_points.length() {
                        a_new_points.set_value(i + 1, *the_points.value(i));
                        a_new_params.set_value(i + 1, *the_params.value(i));
                        a_new_points2d.set_value(i + 1, *the_points2d.value(i));
                    }

                    *the_points = a_new_points;
                    *the_params = a_new_params;
                    *the_points2d = a_new_points2d;
                    *the_index += 1;
                } else {
                    *the_to_adjust = true;
                }
            }
        }
        if *the_to_adjust {
            *the_cur_coord = corrected_cur_coord;
            set_coord(
                the_points2d.change_value(*the_index),
                the_ind_coord,
                *the_cur_coord,
            );
        }
    }

    /// OCCT .cxx L2152-2216: interpolatePCurve(theNbPnt, thePoints2d,
    /// theParams) — interpolates the 2D curve from points.
    pub(crate) fn interpolate_p_curve(
        &self,
        the_nb_pnt: i32,
        the_points2d: &ArrayOfPnt2d,
        the_params: &ArrayOfReal,
    ) -> Option<Curve2d> {
        let the_tolerance2d = self.my_preci / (100.0 * the_nb_pnt as f64);

        // try { OCC_CATCH_SIGNALS ... } — the OCCT exception arm maps to the
        // catch_unwind Err arm (null curve).
        let attempted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            // Convert to HArrays for Geom2dAPI_Interpolate
            let mut a_pnts2d: Vec<DVec2> = Vec::new();
            let mut a_params: Vec<f64> = Vec::new();

            for i in 1..=the_nb_pnt {
                a_pnts2d.push(*the_points2d.value(i as usize));
                a_params.push(*the_params.value(i as usize));
            }

            // Check for coincident points
            let mut a_tmp_pnts2d = the_points2d.clone();
            let mut a_tmp_params = the_params.clone();
            let mut a_preci = the_tolerance2d;
            self.check_points2d(&mut a_tmp_pnts2d, &mut a_tmp_params, &mut a_preci);

            if a_tmp_pnts2d.length() < 2 {
                return None;
            }

            // Rebuild HArrays if points were removed
            if a_tmp_pnts2d.length() != the_nb_pnt as usize {
                a_pnts2d = (1..=a_tmp_pnts2d.length())
                    .map(|i| *a_tmp_pnts2d.value(i))
                    .collect();
                a_params = (1..=a_tmp_params.length())
                    .map(|i| *a_tmp_params.value(i))
                    .collect();
            }

            let mut my_inter_pol2d =
                Geom2dInterpolate::new_with_params(a_pnts2d, a_params, false, a_preci);
            my_inter_pol2d.perform();
            if my_inter_pol2d.is_done() {
                Some(Curve2d::BSpline(my_inter_pol2d.curve()))
            } else {
                None
            }
        }));
        match attempted {
            Ok(a_c2d) => a_c2d,
            Err(_an_exception) => None, // aC2D.Nullify()
        }
    }

    /// OCCT .cxx L2220-2287: approximatePCurve(thePoints2d, theParams) —
    /// approximates the 2D curve from points.
    // (Standard_EXPORT member; no internal caller in the .cxx — the warning
    // is silenced, the method stays in the OCCT API surface.)
    #[allow(dead_code)]
    pub(crate) fn approximate_p_curve(
        &self,
        the_points2d: &ArrayOfPnt2d,
        the_params: &ArrayOfReal,
    ) -> Option<Curve2d> {
        let the_tolerance2d = self.my_preci;

        // try { OCC_CATCH_SIGNALS ... }
        let attempted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut a_tmp_pnts2d = the_points2d.clone();
            let mut a_tmp_params = the_params.clone();
            let mut a_preci = the_tolerance2d;
            self.check_points2d(&mut a_tmp_pnts2d, &mut a_tmp_params, &mut a_preci);

            let number_pnt = a_tmp_pnts2d.length();
            if number_pnt < 2 {
                return None;
            }

            let mut points3d = ArrayOfPnt::new(1, number_pnt, DVec3::ZERO);
            let mut params = ArrayOfReal::new(1, number_pnt, 0.0);

            for i in 1..=number_pnt {
                let pnt2d = *a_tmp_pnts2d.value(i);
                points3d.set_value(i, DVec3::new(pnt2d.x, pnt2d.y, 0.0));
                params.set_value(i, *a_tmp_params.value(i));
            }

            // OCCT L2252: GeomAPI_PointsToBSpline appr(points3d, params, 1, 10,
            // GeomAbs_C1, theTolerance2d).
            let appr = PointsToBSpline::with_points_params(
                &points3d,
                &params,
                1,
                10,
                GeomAbsShape::C1,
                the_tolerance2d,
            );
            let crv3d = appr.curve();

            let nb_poles = crv3d.control_points.len();
            let mut poles2d = Vec::with_capacity(nb_poles);

            for i in 0..nb_poles {
                poles2d.push(DVec2::new(
                    crv3d.control_points[i].x,
                    crv3d.control_points[i].y,
                ));
            }

            // OCCT L2264-2273: weights = crv3d->WeightsArray();
            // knots = crv3d->Knots(); multiplicities = crv3d->Multiplicities();
            // aC2D = new Geom2d_BSplineCurve(poles2d, weights, knots,
            // multiplicities, degree, periodic).  The rcad BSplineCurve2
            // stores the flat knot vector (the OCCT Knots() array is flat as
            // well) and has no periodic flag (arch. diff.); the
            // PointsToBSpline result is non-periodic.
            Some(Curve2d::BSpline(BSplineCurve2 {
                degree: crv3d.degree,
                knots: crv3d.knots.clone(),
                control_points: poles2d,
                weights: crv3d.weights.clone(),
            }))
        }));
        match attempted {
            Ok(a_c2d) => a_c2d,
            Err(_an_exception) => None, // aC2D.Nullify()
        }
    }

    /// OCCT .cxx L2291-2369: checkPoints(thePoints, theParams, thePreci) —
    /// checks and removes coincident 3D points.
    // (Standard_EXPORT member; no internal caller in the .cxx.)
    #[allow(dead_code)]
    pub(crate) fn check_points(
        &self,
        the_points: &mut ArrayOfPnt,
        the_params: &mut ArrayOfReal,
        the_preci: &mut f64,
    ) {
        let first_elem = the_points.lower();
        let last_elem = the_points.upper();

        let mut nb_pnt_dropped = 0;
        let mut last_valid = first_elem;

        let mut tmp_param = Array1::<i32>::new(first_elem, last_elem, 1);
        for i in first_elem..=last_elem {
            tmp_param.set_value(i, 1);
        }

        let mut dist_min2 = REAL_LAST;
        let mut prev = *the_points.value(last_valid);

        for i in first_elem + 1..=last_elem {
            let curr = *the_points.value(i);
            let cur_dist2 = prev.distance_squared(curr);
            if cur_dist2 < GP_RESOLUTION {
                nb_pnt_dropped += 1;
                if i == last_elem {
                    tmp_param.set_value(last_valid, 0);
                } else {
                    tmp_param.set_value(i, 0);
                }
            } else {
                if cur_dist2 < dist_min2 {
                    dist_min2 = cur_dist2;
                }
                last_valid = i;
                prev = curr;
            }
        }

        if dist_min2 < REAL_LAST {
            *the_preci = 0.9 * dist_min2.sqrt();
        }

        if nb_pnt_dropped == 0 {
            return;
        }

        let new_last = last_elem - nb_pnt_dropped;
        if new_last + 1 < first_elem + 2 {
            return;
        }

        let mut new_pnts = ArrayOfPnt::new(first_elem, new_last, DVec3::ZERO);
        let mut new_params = ArrayOfReal::new(first_elem, new_last, 0.0);
        let mut new_curr = first_elem;

        for i in first_elem..=last_elem {
            if *tmp_param.value(i) == 1 {
                new_pnts.set_value(new_curr, *the_points.value(i));
                new_params.set_value(new_curr, *the_params.value(i));
                new_curr += 1;
            }
        }

        *the_points = new_pnts;
        *the_params = new_params;
    }

    /// OCCT .cxx L2373-2457: checkPoints2d(thePoints2d, theParams, thePreci) —
    /// checks and removes coincident 2D points.
    pub(crate) fn check_points2d(
        &self,
        the_points2d: &mut ArrayOfPnt2d,
        the_params: &mut ArrayOfReal,
        the_preci: &mut f64,
    ) {
        let first_elem = the_points2d.lower();
        let last_elem = the_points2d.upper();

        let mut nb_pnt_dropped = 0;
        let mut last_valid = first_elem;

        let mut tmp_param = Array1::<i32>::new(first_elem, last_elem, 1);
        for i in first_elem..=last_elem {
            tmp_param.set_value(i, 1);
        }

        let mut dist_min2 = REAL_LAST;
        let mut prev = *the_points2d.value(last_valid);

        for i in first_elem + 1..=last_elem {
            let curr = *the_points2d.value(i);
            let cur_dist2 = prev.distance_squared(curr);
            if cur_dist2 < GP_RESOLUTION {
                nb_pnt_dropped += 1;
                if i == last_elem {
                    tmp_param.set_value(last_valid, 0);
                } else {
                    tmp_param.set_value(i, 0);
                }
            } else {
                if cur_dist2 < dist_min2 {
                    dist_min2 = cur_dist2;
                }
                last_valid = i;
                prev = curr;
            }
        }

        if dist_min2 < REAL_LAST {
            *the_preci = 0.9 * dist_min2.sqrt();
        }

        if nb_pnt_dropped == 0 {
            return;
        }

        let mut new_last = last_elem - nb_pnt_dropped;
        if new_last + 1 < first_elem + 2 {
            // Create minimal length pcurve
            tmp_param.set_value(first_elem, 1);
            tmp_param.set_value(last_elem, 1);
            let mut last_pnt = *the_points2d.value(last_elem);
            last_pnt += DVec2::new(*the_preci, *the_preci);
            *the_points2d.change_value(last_elem) = last_pnt;
            new_last = first_elem + 1;
        }

        let mut new_pnts = ArrayOfPnt2d::new(first_elem, new_last, DVec2::ZERO);
        let mut new_params = ArrayOfReal::new(first_elem, new_last, 0.0);
        let mut new_curr = first_elem;

        for i in first_elem..=last_elem {
            if *tmp_param.value(i) == 1 {
                new_pnts.set_value(new_curr, *the_points2d.value(i));
                new_params.set_value(new_curr, *the_params.value(i));
                new_curr += 1;
            }
        }

        *the_points2d = new_pnts;
        *the_params = new_params;
    }

    /// OCCT .cxx L2461-2785: isAnIsoparametric(theNbPnt, thePoints, theParams,
    /// theIsTypeU, theP1OnIso, theValueP1, theP2OnIso, theValueP2,
    /// theIsoPar2d3d, theCIso, theT1, theT2, theParamsOut) — detects if the
    /// curve is isoparametric (U=const or V=const).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn is_an_isoparametric(
        &self,
        the_nb_pnt: i32,
        the_points: &ArrayOfPnt,
        the_params: &ArrayOfReal,
        the_is_type_u: &mut bool,
        the_p1_on_iso: &mut bool,
        the_value_p1: &mut DVec2,
        the_p2_on_iso: &mut bool,
        the_value_p2: &mut DVec2,
        the_iso_par2d3d: &mut bool,
        the_c_iso: &mut Option<Curve3>,
        the_t1: &mut f64,
        the_t2: &mut f64,
        the_params_out: &mut ArrayOfReal,
    ) -> bool {
        let a_surf_sa = self.my_surf.as_ref().unwrap();

        // try { OCC_CATCH_SIGNALS ... } — OCCT returns false on exception;
        // the rcad arm mirrors it through catch_unwind.
        let attempted = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| -> bool {
            let prec = CONFUSION;

            let mut iso_param = false;
            *the_iso_par2d3d = false;

            // OCCT L2485: mySurf->Bounds(U1, U2, V1, V2).
            let [mut u1, mut u2, mut v1, mut v2] = a_surf_sa.bounds();

            // OCCT L2487-2492: rectangular trimmed surfaces re-read their own
            // bounds.
            if let Surface3::Trimmed(s_trim) = a_surf_sa.surface() {
                let trim = s_trim.trim;
                u1 = trim[0];
                u2 = trim[1];
                v1 = trim[2];
                v2 = trim[3];
            }

            let mut pt = DVec3::ZERO;
            let mut mpt = [0i32; 2];
            let mut t = 0.0f64;
            let mut tpar = [0.0f64; 2];
            let mut iso_value = 0.0f64;
            let mut mindist2;
            let mut mind2 = [0.0f64; 2];
            mindist2 = 4.0 * prec * prec;
            mind2[0] = mindist2;
            mind2[1] = mindist2;

            *the_p1_on_iso = false;
            *the_p2_on_iso = false;

            for j in 1..=4i32 {
                let mut iso_val = 0.0f64;
                let iso_u;
                let c_i: Option<Curve3>;
                let tt1: f64;
                let tt2: f64;

                // OCCT Bnd_Box* selection.
                let a_box;
                match j {
                    1 => {
                        if is_infinite_value(u1) {
                            continue;
                        }
                        c_i = a_surf_sa.u_iso(u1);
                        iso_u = true;
                        iso_val = u1;
                        a_box = a_surf_sa.get_box_uf();
                    }
                    2 => {
                        if is_infinite_value(u2) {
                            continue;
                        }
                        c_i = a_surf_sa.u_iso(u2);
                        iso_u = true;
                        iso_val = u2;
                        a_box = a_surf_sa.get_box_ul();
                    }
                    3 => {
                        if is_infinite_value(v1) {
                            continue;
                        }
                        c_i = a_surf_sa.v_iso(v1);
                        iso_u = false;
                        iso_val = v1;
                        a_box = a_surf_sa.get_box_vf();
                    }
                    _ => {
                        if is_infinite_value(v2) {
                            continue;
                        }
                        c_i = a_surf_sa.v_iso(v2);
                        iso_u = false;
                        iso_val = v2;
                        a_box = a_surf_sa.get_box_vl();
                    }
                }
                let c_i = match c_i {
                    Some(c) => c,
                    None => continue,
                };

                if iso_u {
                    tt1 = v1;
                    tt2 = v2;
                } else {
                    tt1 = u1;
                    tt2 = u2;
                }

                let ext1 = c_i.point_at(tt1);
                let ext2 = c_i.point_at(tt2);

                let extmi = c_i.point_at((tt1 + tt2) / 2.0);
                // gp_Pnt::IsEqual(P, prec) — distance <= prec.
                if ext1.distance(ext2) <= prec && ext1.distance(extmi) <= prec {
                    continue;
                }

                let mut pt_eq_ext1 = false;
                let mut pt_eq_ext2 = false;

                let mut currd2 = [0.0f64; 2];
                let mut tp = [0.0f64; 2];
                let mut mp = [0i32; 2];

                for i in 0..2 {
                    mp[i] = 0;
                    let k = if i == 0 { 1usize } else { the_nb_pnt as usize };

                    currd2[i] = the_points.value(k).distance_squared(ext1);
                    if currd2[i] <= prec * prec && !pt_eq_ext1 {
                        mp[i] = 1;
                        tp[i] = tt1;
                        pt_eq_ext1 = true;
                        continue;
                    }

                    currd2[i] = the_points.value(k).distance_squared(ext2);
                    if currd2[i] <= prec * prec && !pt_eq_ext2 {
                        mp[i] = 2;
                        tp[i] = tt2;
                        pt_eq_ext2 = true;
                        continue;
                    }

                    if is_spherical_surface(a_surf_sa.surface()) && !iso_u {
                        continue;
                    }

                    if a_box.is_out_point(*the_points.value(k)) {
                        continue;
                    }

                    let mut cf = c_i.default_domain()[0];
                    let mut cl = c_i.default_domain()[1];
                    if is_infinite_value(cf) {
                        cf = -1000.0;
                    }
                    if is_infinite_value(cl) {
                        cl = 1000.0;
                    }

                    // OCCT L2634-2635: sac.Project(cI, thePoints(k), prec, pt, t, Cf, Cl).
                    let (proj, param, dist) =
                        shape_analysis_curve_project(&c_i, *the_points.value(k), prec, cf, cl);
                    pt = proj;
                    t = param;
                    currd2[i] = dist * dist;
                    if dist <= prec && t >= cf && t <= cl {
                        mp[i] = 3;
                        tp[i] = t;
                    }
                }

                if mp[0] > 0 && mp[1] > 0 && (tp[0] - tp[1]).abs() < PCONFUSION {
                    continue;
                }

                if mp[0] > 0 && (!*the_p1_on_iso || currd2[0] < mind2[0]) {
                    *the_p1_on_iso = true;
                    mind2[0] = currd2[0];
                    if iso_u {
                        *the_value_p1 = DVec2::new(iso_val, tp[0]);
                    } else {
                        *the_value_p1 = DVec2::new(tp[0], iso_val);
                    }
                }

                if mp[1] > 0 && (!*the_p2_on_iso || currd2[1] < mind2[1]) {
                    *the_p2_on_iso = true;
                    mind2[1] = currd2[1];
                    if iso_u {
                        *the_value_p2 = DVec2::new(iso_val, tp[1]);
                    } else {
                        *the_value_p2 = DVec2::new(tp[1], iso_val);
                    }
                }

                if mp[0] <= 0 || mp[1] <= 0 {
                    continue;
                }

                let md2 = currd2[0] + currd2[1];
                if mindist2 <= md2 {
                    continue;
                }

                mindist2 = md2;
                mpt[0] = mp[0];
                mpt[1] = mp[1];
                tpar[0] = tp[0];
                tpar[1] = tp[1];
                *the_is_type_u = iso_u;
                iso_value = iso_val;
                *the_c_iso = Some(c_i);
                *the_t1 = tt1;
                *the_t2 = tt2;
            }

            if mpt[0] > 0 && mpt[1] > 0 {
                *the_p1_on_iso = true;
                *the_p2_on_iso = true;
                if *the_is_type_u {
                    *the_value_p1 = DVec2::new(iso_value, tpar[0]);
                    *the_value_p2 = DVec2::new(iso_value, tpar[1]);
                } else {
                    *the_value_p1 = DVec2::new(tpar[0], iso_value);
                    *the_value_p2 = DVec2::new(tpar[1], iso_value);
                }

                if mpt[0] != 3 && mpt[1] != 3 {
                    *the_iso_par2d3d = true;
                    let c_iso_ref = the_c_iso.as_ref().unwrap();
                    let mut i = 2i32;
                    while i < the_nb_pnt && *the_iso_par2d3d {
                        if tpar[1] > tpar[0] {
                            t = *the_params.value(i as usize);
                        } else {
                            t = *the_t1 + *the_t2 - *the_params.value(i as usize);
                        }
                        pt = c_iso_ref.point_at(t);
                        // OCCT: if (!thePoints(i).IsEqual(pt, prec)) theIsoPar2d3d = false.
                        if the_points.value(i as usize).distance(pt) > prec {
                            *the_iso_par2d3d = false;
                        }
                        i += 1;
                    }
                }

                if *the_iso_par2d3d {
                    iso_param = true;
                } else {
                    let mut prev_param = tpar[0];
                    let c_iso_ref = the_c_iso.as_ref().unwrap();
                    let mut cf = c_iso_ref.default_domain()[0];
                    let mut cl = c_iso_ref.default_domain()[1];
                    let mut iso_by_distance = true;
                    if is_infinite_value(cf) {
                        cf = -1000.0;
                    }
                    if is_infinite_value(cl) {
                        cl = 1000.0;
                    }

                    let mut i = 2i32;
                    while i < the_nb_pnt && iso_by_distance {
                        // OCCT L2758-2759: sac.NextProject(prevParam, theCIso,
                        // thePoints(i), prec, pt, t, Cf, Cl, false).
                        let (proj, param, dist) = shape_analysis_curve_next_project(
                            prev_param,
                            c_iso_ref,
                            *the_points.value(i as usize),
                            prec,
                            cf,
                            cl,
                        );
                        pt = proj;
                        t = param;
                        prev_param = t;
                        the_params_out.set_value(i as usize, t);
                        if dist > prec || t < cf || t > cl {
                            iso_by_distance = false;
                        }
                        i += 1;
                    }
                    if iso_by_distance {
                        iso_param = true;
                    }
                }
            }
            iso_param
        }));
        match attempted {
            Ok(res) => res,
            Err(_an_exception) => false,
        }
    }
}
