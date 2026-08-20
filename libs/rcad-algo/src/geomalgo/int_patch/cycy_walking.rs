//! Walking-line helpers for the IntCyCy (cylinder-cylinder) numeric engine —
//! 1:1 translation of OCCT `IntPatch_ImpImpIntersection.cxx`:
//!   - AddPointIntoWL (L5615-5817)
//!   - CyCyNoGeometric (L6573-7880)
//!
//! rcad data-model notes:
//!   - `IntSurf_Quadric` -> `Quadric`; `stCoeffsValue` -> `StCoeffsValue`.
//!   - `gp_Pnt2d` -> `[f64; 2]`.
//!   - `IntSurf_LineOn2S` -> `cycy_common::WLine` (1-based semantics).
//!   - `IntSurf_PntOn2S` -> `WLinePnt { p3d, u1, v1, u2, v2 }`.
//!   - `Epsilon(1.0)` = DBL_EPSILON = `f64::EPSILON`.

use rcad_kernel::geom::CylindricalSurface;
use rcad_kernel::precision::{PCONFUSION, SQUARE_CONFUSION};

use super::cycy_coeffs::{
    StCoeffsValue, cyl_cyl_compute_parameters, cyl_cyl_compute_parameters_u2,
    cyl_cyl_compute_parameters_v,
};
use super::cycy_common::{
    Mat35, WLine, inscribe_interval, inscribe_point, is_equal, precision_infinite, real_first,
};
use super::cycy_boundaries::{
    WorkWithBoundaries, critical_points_computing, seek_additional_points,
};
use super::imp_imp_intersection::IntStatus;
use super::{IntPatchLine, IntPatchPoint, WLineType};
use super::WLinePnt;
use crate::geomalgo::int_surf::quadric::Quadric;

/// OCCT IntSurf_PntOn2S::IsSame(Other, Tol3D, Tol2D) (IntSurf_PntOn2S.cxx).
fn is_same(p: &WLinePnt, other: &WLinePnt, tol_3d: f64, tol_2d: f64) -> bool {
    if p.p3d.distance_squared(other.p3d) > tol_3d * tol_3d {
        return false;
    }
    if tol_2d < 0.0 {
        // We need not compare 2D-coordinates of the points.
        return true;
    }
    if (p.u1 - other.u1).abs() > tol_2d || (p.v1 - other.v1).abs() > tol_2d {
        return false;
    }
    if (p.u2 - other.u2).abs() > tol_2d || (p.v2 - other.v2).abs() > tol_2d {
        return false;
    }
    true
}

/// OCCT AddPointIntoWL (L5615-5817).
/// Surf1 is the surface whose U-parameter is variable.  If theFlBefore == TRUE
/// the U1-parameter of the added point may be less than the U1-parameter of the
/// previously added point (in general U1-parameter is always increased).
/// If theOnlyCheck == TRUE no point is added to the line.
#[allow(clippy::too_many_arguments)]
pub fn add_point_into_wl(
    quad1: &Quadric,
    quad2: &Quadric,
    coeffs: &StCoeffsValue,
    is_reverse: bool,
    is_precise: bool,
    pnt_on_surf1: [f64; 2],
    pnt_on_surf2: [f64; 2],
    uf_surf1: f64,
    ul_surf1: f64,
    uf_surf2: f64,
    ul_surf2: f64,
    vf_surf1: f64,
    vl_surf1: f64,
    vf_surf2: f64,
    vl_surf2: f64,
    period_of_surf1: f64,
    line: &mut WLine,
    wl_index: usize,
    tol_3d: f64,
    tol_2d: f64,
    fl_before: bool,
    only_check: bool,
) -> bool {
    // Check if the point is in the domain or can be inscribed in the domain
    // after adjusting.
    let a_pt1 = quad1.value(pnt_on_surf1[0], pnt_on_surf1[1]);
    let a_pt2 = quad2.value(pnt_on_surf2[0], pnt_on_surf2[1]);

    let mut a_u1par = pnt_on_surf1[0];

    // aU1par always increases.  Therefore, we must reduce its value in order to
    // continue creation of WLine.
    let fl_force = a_u1par > 0.5 * (uf_surf1 + ul_surf1);
    if !inscribe_point(uf_surf1, ul_surf1, &mut a_u1par, tol_2d, period_of_surf1, fl_force) {
        return false;
    }

    if (line.nb_points() > 0)
        && ((ul_surf1 - uf_surf1) >= (period_of_surf1 - tol_2d))
        && (((a_u1par + period_of_surf1 - ul_surf1) <= tol_2d)
            || ((a_u1par - period_of_surf1 - uf_surf1) >= tol_2d))
    {
        // aU1par can be adjusted to both theUlSurf1 and theUfSurf1 with equal
        // possibilities.  This fragment allows choosing the correct parameter.
        #[allow(unused_assignments)]
        let mut a_u1 = 0.0;
        let mut _a_v1 = 0.0;
        let plast = line.value(line.nb_points());
        if is_reverse {
            a_u1 = plast.u2;
            _a_v1 = plast.v2;
        } else {
            a_u1 = plast.u1;
            _a_v1 = plast.v1;
        }
        let a_delta = a_u1 - a_u1par;
        if 2.0 * a_delta.abs() > period_of_surf1 {
            a_u1par += a_delta.signum() * period_of_surf1;
        }
    }

    let mut a_u2par = pnt_on_surf2[0];
    if !inscribe_point(uf_surf2, ul_surf2, &mut a_u2par, tol_2d, period_of_surf1, false) {
        return false;
    }

    let a_v1par = pnt_on_surf1[1];
    if (a_v1par - vl_surf1 > tol_2d) || (vf_surf1 - a_v1par > tol_2d) {
        return false;
    }
    let a_v2par = pnt_on_surf2[1];
    if (a_v2par - vl_surf2 > tol_2d) || (vf_surf2 - a_v2par > tol_2d) {
        return false;
    }

    // Get intersection point and add it in the WL.
    let mid = 0.5 * (a_pt1 + a_pt2);
    let a_pnt = if is_reverse {
        WLinePnt { p3d: mid, u1: a_u2par, v1: a_v2par, u2: a_u1par, v2: a_v1par }
    } else {
        WLinePnt { p3d: mid, u1: a_u1par, v1: a_v1par, u2: a_u2par, v2: a_v2par }
    };

    let mut a_nb_pnts = line.nb_points();
    if a_nb_pnts > 0 {
        let a_plast = line.value(a_nb_pnts).clone();
        let (a_ul, _a_vl) = if is_reverse { (a_plast.u2, a_plast.v2) } else { (a_plast.u1, a_plast.v1) };

        if !fl_before && a_u1par <= a_ul {
            // Parameter value must be increased if theFlBefore == FALSE.
            a_u1par += period_of_surf1;
            // The condition is the same as in InscribePoint(...).
            if (uf_surf1 - a_u1par > tol_2d) || (a_u1par - ul_surf1 > tol_2d) {
                // New aU1par is out of target interval; go back to old value.
                return false;
            }
        }

        if only_check {
            return true;
        }

        // theTol2D is the minimal step along the parameter changed; therefore,
        // if we apply this minimal step two neighbor points will be always
        // "same".  Consequently, we should reduce tolerance for IsSame checking.
        let a_d_tol = 1.0 - f64::EPSILON;
        if is_same(&a_pnt, &a_plast, tol_3d * a_d_tol, tol_2d * a_d_tol) {
            line.remove_point(a_nb_pnts);
        }
    }

    if only_check {
        return true;
    }

    line.append(a_pnt);

    if !is_precise {
        return true;
    }

    // Try to precise existing WLine.
    a_nb_pnts = line.nb_points();
    if a_nb_pnts >= 3 {
        let (a_u3, a_u2, a_u1) = if is_reverse {
            let p3 = line.value(a_nb_pnts);
            let p2 = line.value(a_nb_pnts - 1);
            let p1 = line.value(a_nb_pnts - 2);
            (p3.u2, p2.u2, p1.u2)
        } else {
            let p3 = line.value(a_nb_pnts);
            let p2 = line.value(a_nb_pnts - 1);
            let p1 = line.value(a_nb_pnts - 2);
            (p3.u1, p2.u1, p1.u1)
        };

        let a_step_prev = a_u2 - a_u1;
        let a_step = a_u3 - a_u2;
        let a_delta_step = (a_step_prev / a_step) as i32;

        if (1 < a_delta_step) && (a_delta_step < 2000) {
            // Add new points in case of non-uniform distribution of existing points.
            super::cycy_boundaries::seek_additional_points(
                quad1,
                quad2,
                line,
                coeffs,
                wl_index,
                a_delta_step as usize,
                a_nb_pnts - 2,
                a_nb_pnts - 1,
                tol_2d,
                period_of_surf1,
                is_reverse,
            );
        }
    }

    true
}

/// OCCT `gp_Lin::Distance(Other)` — distance between two infinite 3D lines.
fn line_line_distance(loc1: glam::DVec3, dir1: glam::DVec3, loc2: glam::DVec3, dir2: glam::DVec3) -> f64 {
    let n = dir1.cross(dir2);
    let sq = n.length_squared();
    let d = loc2 - loc1;
    if sq < f64::MIN_POSITIVE {
        // Parallel lines: distance from loc2 to line1.
        return d.cross(dir1).length() / dir1.length().max(f64::MIN_POSITIVE);
    }
    (d.dot(n)).abs() / n.length()
}

/// OCCT gp_Ax1::IsNormal(Other, AngularTolerance).
fn axes_are_normal(dir1: glam::DVec3, dir2: glam::DVec3, tol_ang: f64) -> bool {
    (dir1.angle_between(dir2) - std::f64::consts::FRAC_PI_2).abs() <= tol_ang
}

/// OCCT WLFStatus enum (L6755-6762) — status of the WLine under construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WlfStatus {
    /// No points have been added in WL.
    Absent = 0,
    /// WL contains at least one point.
    Exist = 1,
    /// WL has been finished in some critical point; we should start a new line.
    Broken = 2,
}

/// OCCT IntPatch_Point PntOn2S compare helper (IntSurf_PntOn2S::IsSame).
fn pnt2s_is_same(a: &IntPatchPoint, b: &IntPatchPoint, tol_3d: f64) -> bool {
    let pa = WLinePnt { p3d: a.p1, u1: a.u1, v1: a.v1, u2: a.u2, v2: a.v2 };
    let pb = WLinePnt { p3d: b.p1, u1: b.u1, v1: b.v1, u2: b.u2, v2: b.v2 };
    is_same(&pa, &pb, tol_3d, -1.0)
}

/// OCCT CyCyNoGeometric (L6573-7880) — the walking-line generator for the
/// general (non-analytic) cylinder-cylinder intersection.
///
/// Returns IntStatus; on success the resulting WLines are appended to `slin`
/// and isolated tangent points to `spnt`.
#[allow(clippy::too_many_lines, unused_assignments)]
pub fn cy_cy_no_geometric(
    cyl1: &CylindricalSurface,
    cyl2: &CylindricalSurface,
    bw: &WorkWithBoundaries,
    ranges: &mut [super::cycy_common::BndRange; 2],
    nb_of_ranges: usize,
    is_empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> IntStatus {
    let uv1 = bw.uv_s1();
    let uv2 = bw.uv_s2();
    let (a_usurf1f, a_vsurf1f, a_usurf1l, a_vsurf1l) = (uv1[0], uv1[1], uv1[2], uv1[3]);
    let (a_usurf2f, a_vsurf2f, a_usurf2l, a_vsurf2l) = (uv2[0], uv2[1], uv2[2], uv2[3]);

    let mut a_range_s1 = super::cycy_common::BndRange::new();
    let mut a_range_s2 = super::cycy_common::BndRange::new();
    bw.boundary_estimation(cyl1, cyl2, &mut a_range_s1, &mut a_range_s2);
    if a_range_s1.is_void() || a_range_s2.is_void() {
        return IntStatus::OK;
    }

    {
        // We should return fail status from the intersector if the result
        // should be an infinite curve of non-analytical type.  The limit for
        // the extent is the radius divided by 1e+2 and multiplied by 1e+7.
        // Thus, taking into account the number of valuable digits (15), we
        // provide reliable computations with an error not exceeding R/100.
        let a_f = 1.0e+5;
        let a_max_v1_range = a_f * cyl1.radius;
        let a_max_v2_range = a_f * cyl2.radius;
        if (a_range_s1.delta() > a_max_v1_range) || (a_range_s2.delta() > a_max_v2_range) {
            return IntStatus::InfiniteSectionCurve;
        }
    }

    // Checking parameters of cylinders in order to define "good intersection":
    // axes are almost perpendicular and one radius is much smaller than the
    // other and the small cylinder is "inside" the big one.
    let mut is_good_intersection = false;
    let mut an_optdu = 0.0;
    loop {
        let a_to_much_coeff = 3.0;
        let a_crit_angle = std::f64::consts::PI / 18.0; // 10 degree
        let an_r1 = cyl1.radius;
        let an_r2 = cyl2.radius;
        let (mut an_rmin, mut an_rmax) = (0.0, 0.0);
        // Radius criterion
        if an_r1 > a_to_much_coeff * an_r2 {
            an_rmax = an_r1;
            an_rmin = an_r2;
        } else if an_r2 > a_to_much_coeff * an_r1 {
            an_rmax = an_r2;
            an_rmin = an_r1;
        } else {
            break;
        }
        // Angle criterion
        let an_ax1 = cyl1.axis.normalize_or_zero();
        let an_ax2 = cyl2.axis.normalize_or_zero();
        if !axes_are_normal(an_ax1, an_ax2, a_crit_angle) {
            break;
        }
        // Placement criterion
        let a_dist = line_line_distance(cyl1.origin, an_ax1, cyl2.origin, an_ax2);
        if a_dist > an_rmax / 2.0 {
            break;
        }
        is_good_intersection = true;
        // Estimation of "optimal" du
        // Relative deflection, absolute deflection is Rmin*aDeflection
        let a_deflection = 0.001;
        let mut a_nb_p = 3;
        if an_rmin * a_deflection > 1.0e-3 {
            let mut an_angle = 1.0e0 - a_deflection;
            an_angle = 2.0e0 * an_angle.acos();
            a_nb_p = (2.0 * std::f64::consts::PI / an_angle) as i32 + 1;
        }
        an_optdu = 2.0 * std::f64::consts::FRAC_PI_2 / (a_nb_p - 1) as f64;
        break;
    }

    let an_equation_coeffs = bw.si_coeffs();
    let a_quad1 = bw.get_q_surface(1);
    let a_quad2 = bw.get_q_surface(2);
    let is_reversed = bw.is_reversed();
    let a_tol_2d = bw.get_2d_tolerance();
    let a_tol_3d = bw.get_3d_tolerance();
    let a_period = 2.0 * std::f64::consts::PI;
    let mut a_nb_max_points = 1000usize;
    let mut a_nb_min_points = 200usize;
    let du;
    if is_good_intersection {
        du = an_optdu;
        a_nb_max_points = 200;
        a_nb_min_points = 50;
    } else {
        du = 2.0 * std::f64::consts::PI / a_nb_max_points as f64;
    }
    let a_nb_pts = (((a_usurf1l - a_usurf1f) / du) as usize + 1)
        .min((20.0 * cyl1.radius) as usize);
    let a_nb_points = a_nb_min_points.max(a_nb_pts).min(a_nb_max_points);
    let a_step_min = a_tol_2d.max(PCONFUSION);
    let a_step_max = if (a_usurf1l - a_usurf1f) > std::f64::consts::PI / 100.0 {
        (a_usurf1l - a_usurf1f) / a_nb_points as f64
    } else {
        a_usurf1l - a_usurf1f
    };

    // The main idea of the algorithm is to change U1-parameter (U-parameter of
    // theCyl1) from aU1f to aU1l with some (adaptive) step and to obtain the
    // set of intersection points.
    for i in 0..nb_of_ranges {
        if ranges[i].is_void() {
            continue;
        }
        inscribe_interval(a_usurf1f, a_usurf1l, &mut ranges[i], a_tol_2d, a_period);
    }

    let range1 = ranges[1];
    if ranges[0].union(&range1) {
        // Works only if (theNbOfRanges == 2).
        ranges[1].set_void();
    }

    // Critical points are the values of U1-parameter in the points where WL
    // must be decomposed.  When U1 goes through a critical point its value is
    // set up to this parameter forcefully and the intersection point is added
    // in the line.  After that, the WL is broken (next U1 value will
    // correspond to the new WL).
    const A_NB_CRIT_POINTS_MAX: usize = 12;
    let mut an_u1crit = [precision_infinite(); A_NB_CRIT_POINTS_MAX];
    let mut a_nb_crit_points = A_NB_CRIT_POINTS_MAX;
    critical_points_computing(
        an_equation_coeffs,
        a_usurf1f,
        a_usurf1l,
        a_usurf2f,
        a_usurf2l,
        a_period,
        a_tol_2d,
        &mut a_nb_crit_points,
        &mut an_u1crit,
    );

    let a_nb_wlines = 2usize;
    for a_cur_interval in 0..nb_of_ranges {
        // Process every continuous region.
        let mut is_added_into_wl = [false; 2];
        let (mut an_uf, an_ul) = match ranges[a_cur_interval].get_bounds() {
            Some(b) => b,
            None => continue,
        };

        let is_delta_period = is_equal(an_ul - an_uf, a_period);

        // Inscribe and sort critical points.
        for i in 0..a_nb_crit_points {
            inscribe_point(an_uf, an_ul, &mut an_u1crit[i], 0.0, a_period, false);
        }
        an_u1crit[..a_nb_crit_points].sort_by(|a, b| a.partial_cmp(b).unwrap());

        while an_uf < an_ul {
            // Change the value of U-parameter on the 1st surface from anUf to
            // anUl (anUf will be modified in the cycle body).  Step is computed
            // adaptively (see comments below).
            let mut a_u2 = [0.0; 2];
            let mut a_v1 = [0.0; 2];
            let mut a_v2 = [0.0; 2];
            let mut a_wl_find_status = [WlfStatus::Absent; 2];
            let mut a_v1_prev = [0.0; 2];
            let mut a_v2_prev = [0.0; 2];
            let mut an_u_expect = [0.0; 2];
            let mut is_adding_wl_enabled = [true; 2];
            let mut a_w_line = [WLine::new(), WLine::new()];

            for i in 0..a_nb_wlines {
                a_wl_find_status[i] = WlfStatus::Absent;
                is_adding_wl_enabled[i] = true;
                a_u2[i] = 0.0;
                a_v1[i] = 0.0;
                a_v2[i] = 0.0;
                a_v1_prev[i] = 0.0;
                a_v2_prev[i] = 0.0;
                an_u_expect[i] = an_uf;
            }

            let mut a_critical_delta = [0.0; A_NB_CRIT_POINTS_MAX];
            for a_crit_pid in 0..a_nb_crit_points {
                a_critical_delta[a_crit_pid] = an_uf - an_u1crit[a_crit_pid];
            }

            let mut an_u1 = an_uf;
            let a_min_critical_param = an_uf;
            let mut an_u1_prev = an_uf;
            let mut is_first = true;

            while an_u1 <= an_ul {
                for i in 0..a_nb_crit_points {
                    if (an_u1 - an_u1crit[i]) * a_critical_delta[i] < 0.0 {
                        // WL has gone through i-th critical point.
                        an_u1 = an_u1crit[i];
                        for j in 0..a_nb_wlines {
                            a_wl_find_status[j] = WlfStatus::Broken;
                            an_u_expect[j] = an_u1;
                        }
                        break;
                    }
                }

                if is_equal(an_u1, an_ul) {
                    for i in 0..a_nb_wlines {
                        a_wl_find_status[i] = WlfStatus::Broken;
                        an_u_expect[i] = an_u1;
                        if is_delta_period {
                            // If isAddedIntoWL[i] == TRUE the WLine contains
                            // only one point (which was the end point of the
                            // previous WLine).  If we add the point found on
                            // the current step, the WLine will contain only two
                            // points, both equal to the points found earlier;
                            // the new WLine would repeat the existing one.
                            // Therefore, forbid building a new line in this case.
                            is_adding_wl_enabled[i] = !is_added_into_wl[i];
                        } else {
                            is_adding_wl_enabled[i] =
                                (a_tol_2d >= (an_u_expect[i] - an_u1))
                                    || (a_wl_find_status[i] == WlfStatus::Absent);
                        }
                    }
                } else {
                    for i in 0..a_nb_wlines {
                        is_adding_wl_enabled[i] =
                            (a_tol_2d >= (an_u_expect[i] - an_u1))
                                || (a_wl_find_status[i] == WlfStatus::Absent);
                    }
                }

                for i in 0..a_nb_wlines {
                    let a_nb_pnts_wl = a_w_line[i].nb_points();

                    if (a_wl_find_status[i] == WlfStatus::Broken)
                        || (a_wl_find_status[i] == WlfStatus::Absent)
                    {
                        // Begin and end of WLine must be on boundary point or
                        // on seam-edge strictly (if it is possible).
                        let mut a_tol = a_tol_2d;
                        cyl_cyl_compute_parameters_u2(
                            an_u1,
                            i as i32,
                            an_equation_coeffs,
                            &mut a_u2[i],
                            Some(&mut a_tol),
                        );
                        inscribe_point(a_usurf2f, a_usurf2l, &mut a_u2[i], a_tol_2d, a_period, false);
                        a_tol = a_tol.max(a_tol_2d);
                        if a_u2[i].abs() <= a_tol {
                            a_u2[i] = 0.0;
                        } else if (a_u2[i] - a_period).abs() <= a_tol {
                            a_u2[i] = a_period;
                        } else if (a_u2[i] - a_usurf2f).abs() <= a_tol {
                            a_u2[i] = a_usurf2f;
                        } else if (a_u2[i] - a_usurf2l).abs() <= a_tol {
                            a_u2[i] = a_usurf2l;
                        }
                    } else {
                        cyl_cyl_compute_parameters_u2(
                            an_u1,
                            i as i32,
                            an_equation_coeffs,
                            &mut a_u2[i],
                            None,
                        );
                        inscribe_point(a_usurf2f, a_usurf2l, &mut a_u2[i], a_tol_2d, a_period, false);
                    }

                    if a_nb_pnts_wl == 0 {
                        // The line has not contained any points yet.
                        if ((a_usurf2f + a_period - a_usurf2l) <= 2.0 * a_tol_2d)
                            && ((a_u2[i] - a_usurf2f).abs() < a_tol_2d
                                || (a_u2[i] - a_usurf2l).abs() < a_tol_2d)
                        {
                            // In this case aU2[i] can have two values: current
                            // aU2[i] or aU2[i]+aPeriod (aU2[i]-aPeriod).  It is
                            // necessary to choose the correct value.
                            let mut is_increasing = true;
                            super::cycy_coeffs::cyl_cyl_monotonicity(
                                an_u1 + a_step_min,
                                i,
                                an_equation_coeffs,
                                a_period,
                                &mut is_increasing,
                            );
                            // If U2(U1) is increasing and U2 is considered to be
                            // equal aUSurf2l, then after the next step (when U1
                            // is increased) U2 will increase too and we will go
                            // out of the surface boundary.  Therefore, if
                            // U2(U1) is increasing then U2 must be equal aUSurf2f.
                            // Analogically for decreasing.
                            if is_increasing {
                                a_u2[i] = a_usurf2f;
                            } else {
                                a_u2[i] = a_usurf2l;
                            }
                        }
                    } else if ((a_usurf2l - a_usurf2f) >= a_period)
                        && ((a_u2[i] - a_usurf2f).abs() < a_tol_2d
                            || (a_u2[i] - a_usurf2l).abs() < a_tol_2d)
                    {
                        // End of the line.
                        let plast = a_w_line[i].value(a_nb_pnts_wl).clone();
                        let (a_u2_prev, _a_v2_prev) =
                            if is_reversed { (plast.u1, plast.v1) } else { (plast.u2, plast.v2) };
                        if 2.0 * (a_u2_prev - a_u2[i]).abs() > a_period {
                            if a_u2_prev > a_u2[i] {
                                a_u2[i] += a_period;
                            } else {
                                a_u2[i] -= a_period;
                            }
                        }
                    }

                    cyl_cyl_compute_parameters_v(
                        an_u1,
                        a_u2[i],
                        an_equation_coeffs,
                        &mut a_v1[i],
                        &mut a_v2[i],
                    );

                    if is_first {
                        a_v1_prev[i] = a_v1[i];
                        a_v2_prev[i] = a_v2[i];
                    }
                }

                is_first = false;

                // Looking for points in WLine.
                let mut is_broken = false;
                for i in 0..a_nb_wlines {
                    if !is_adding_wl_enabled[i] {
                        let mut is_bound_intersect = false;
                        if (a_v1[i] - a_vsurf1f).abs() <= a_tol_2d
                            || (a_v1[i] - a_vsurf1f) * (a_v1_prev[i] - a_vsurf1f) < 0.0
                        {
                            is_bound_intersect = true;
                        } else if (a_v1[i] - a_vsurf1l).abs() <= a_tol_2d
                            || (a_v1[i] - a_vsurf1l) * (a_v1_prev[i] - a_vsurf1l) < 0.0
                        {
                            is_bound_intersect = true;
                        } else if (a_v2[i] - a_vsurf2f).abs() <= a_tol_2d
                            || (a_v2[i] - a_vsurf2f) * (a_v2_prev[i] - a_vsurf2f) < 0.0
                        {
                            is_bound_intersect = true;
                        } else if (a_v2[i] - a_vsurf2l).abs() <= a_tol_2d
                            || (a_v2[i] - a_vsurf2l) * (a_v2_prev[i] - a_vsurf2l) < 0.0
                        {
                            is_bound_intersect = true;
                        }

                        if a_wl_find_status[i] == WlfStatus::Broken {
                            is_broken = true;
                        }

                        if !is_bound_intersect {
                            continue;
                        }
                        an_u_expect[i] = an_u1;
                    }

                    // True if the current point is already in the domain.
                    let is_inscribe = (a_usurf2f - a_u2[i]) <= a_tol_2d
                        && (a_u2[i] - a_usurf2l) <= a_tol_2d
                        && (a_vsurf1f - a_v1[i]) <= a_tol_2d
                        && (a_v1[i] - a_vsurf1l) <= a_tol_2d
                        && (a_vsurf2f - a_v2[i]) <= a_tol_2d
                        && (a_v2[i] - a_vsurf2l) <= a_tol_2d;

                    // isVIntersect == TRUE if the intersection line intersects
                    // two (!) V-bounds of a cylinder (1st or 2nd — no matter).
                    let is_v_intersect = ((a_vsurf1f - a_v1[i]) * (a_vsurf1f - a_v1_prev[i]) < f64::MIN_POSITIVE)
                        && ((a_vsurf1l - a_v1[i]) * (a_vsurf1l - a_v1_prev[i]) < f64::MIN_POSITIVE)
                        || ((a_vsurf2f - a_v2[i]) * (a_vsurf2f - a_v2_prev[i]) < f64::MIN_POSITIVE)
                        && ((a_vsurf2l - a_v2[i]) * (a_vsurf2l - a_v2_prev[i]) < f64::MIN_POSITIVE);

                    // isFound1 == TRUE if the intersection line intersects
                    // V-bounds (First or Last — no matter) of the 1st cylinder;
                    // isFound2 likewise for the 2nd cylinder.
                    let mut is_found1 = false;
                    let mut is_found2 = false;
                    let mut is_force = false;

                    if a_wl_find_status[i] == WlfStatus::Absent {
                        if ((a_usurf2l - a_usurf2f) >= a_period) && ((an_u1 - a_usurf1l).abs() < a_tol_2d) {
                            is_force = true;
                        }
                    }

                    bw.add_boundary_point(
                        &mut a_w_line[i],
                        an_u1,
                        an_u1_prev,
                        a_min_critical_param,
                        a_u2[i],
                        a_v1[i],
                        a_v1_prev[i],
                        a_v2[i],
                        a_v2_prev[i],
                        i,
                        is_force,
                        &mut is_found1,
                        &mut is_found2,
                    );

                    let is_prev_v_bound = !is_v_intersect
                        && ((a_v1_prev[i] - a_vsurf1f).abs() <= a_tol_2d
                            || (a_v1_prev[i] - a_vsurf1l).abs() <= a_tol_2d
                            || (a_v2_prev[i] - a_vsurf2f).abs() <= a_tol_2d
                            || (a_v2_prev[i] - a_vsurf2l).abs() <= a_tol_2d);

                    a_v1_prev[i] = a_v1[i];
                    a_v2_prev[i] = a_v2[i];

                    if (a_wl_find_status[i] == WlfStatus::Exist)
                        && (is_found1 || is_found2)
                        && !is_prev_v_bound
                    {
                        a_wl_find_status[i] = WlfStatus::Broken; // start a new line
                    } else if is_inscribe {
                        if (a_wl_find_status[i] == WlfStatus::Absent) && (is_found1 || is_found2) {
                            a_wl_find_status[i] = WlfStatus::Exist;
                        }
                        if (a_wl_find_status[i] != WlfStatus::Broken)
                            || (a_w_line[i].nb_points() >= 1)
                            || is_equal(an_u1, an_ul)
                        {
                            if a_w_line[i].nb_points() > 0 {
                                let plast = a_w_line[i].value(a_w_line[i].nb_points()).clone();
                                let (a_u2_p, _a_v2_p) =
                                    if is_reversed { (plast.u1, plast.v1) } else { (plast.u2, plast.v2) };
                                let a_delta = a_u2[i] - a_u2_p;
                                if 2.0 * a_delta.abs() > a_period {
                                    if a_delta > 0.0 {
                                        a_u2[i] -= a_period;
                                    } else {
                                        a_u2[i] += a_period;
                                    }
                                }
                            }

                            if add_point_into_wl(
                                a_quad1,
                                a_quad2,
                                an_equation_coeffs,
                                is_reversed,
                                true,
                                [an_u1, a_v1[i]],
                                [a_u2[i], a_v2[i]],
                                a_usurf1f,
                                a_usurf1l,
                                a_usurf2f,
                                a_usurf2l,
                                a_vsurf1f,
                                a_vsurf1l,
                                a_vsurf2f,
                                a_vsurf2l,
                                a_period,
                                &mut a_w_line[i],
                                i,
                                a_tol_3d,
                                a_tol_2d,
                                is_force,
                                false,
                            ) {
                                if a_wl_find_status[i] == WlfStatus::Absent {
                                    a_wl_find_status[i] = WlfStatus::Exist;
                                }
                            } else if !is_found1 && !is_found2 {
                                // We do not add any point while doing this iteration.
                                if a_wl_find_status[i] == WlfStatus::Exist {
                                    a_wl_find_status[i] = WlfStatus::Broken;
                                }
                            }
                        }
                    } else {
                        // We do not add any point while doing this iteration.
                        if a_wl_find_status[i] == WlfStatus::Exist {
                            a_wl_find_status[i] = WlfStatus::Broken;
                        }
                    }

                    if a_wl_find_status[i] == WlfStatus::Broken {
                        is_broken = true;
                    }
                }

                if is_broken {
                    // Current lines are filled; go to the next lines.
                    an_uf = an_u1;

                    let mut is_added = true;
                    for i in 0..a_nb_wlines {
                        if is_adding_wl_enabled[i] {
                            continue;
                        }
                        is_added = false;

                        let mut is_found1 = false;
                        let mut is_found2 = false;
                        bw.add_boundary_point(
                            &mut a_w_line[i],
                            an_u1,
                            an_u1_prev,
                            a_min_critical_param,
                            a_u2[i],
                            a_v1[i],
                            a_v1_prev[i],
                            a_v2[i],
                            a_v2_prev[i],
                            i,
                            false,
                            &mut is_found1,
                            &mut is_found2,
                        );
                        if is_found1 || is_found2 {
                            is_added = true;
                        }

                        if a_w_line[i].nb_points() > 0 {
                            let plast = a_w_line[i].value(a_w_line[i].nb_points()).clone();
                            let (a_u2_p, _a_v2_p) =
                                if is_reversed { (plast.u1, plast.v1) } else { (plast.u2, plast.v2) };
                            let a_delta = a_u2[i] - a_u2_p;
                            if 2.0 * a_delta.abs() > a_period {
                                if a_delta > 0.0 {
                                    a_u2[i] -= a_period;
                                } else {
                                    a_u2[i] += a_period;
                                }
                            }
                        }

                        if add_point_into_wl(
                            a_quad1,
                            a_quad2,
                            an_equation_coeffs,
                            is_reversed,
                            true,
                            [an_u1, a_v1[i]],
                            [a_u2[i], a_v2[i]],
                            a_usurf1f,
                            a_usurf1l,
                            a_usurf2f,
                            a_usurf2l,
                            a_vsurf1f,
                            a_vsurf1l,
                            a_vsurf2f,
                            a_vsurf2l,
                            a_period,
                            &mut a_w_line[i],
                            i,
                            a_tol_3d,
                            a_tol_2d,
                            false,
                            false,
                        ) {
                            is_added = true;
                        }
                    }

                    if !is_added {
                        // Before breaking the WL, we must complete it correctly
                        // (e.g. to prolong to the surface boundary).  Therefore,
                        // we take the point last added in some WL (having maximal
                        // U1-parameter) and try to add it in the current WL.
                        let mut an_umax_added = real_first();
                        {
                            let mut is_changed = false;
                            for i in 0..a_nb_wlines {
                                if (a_wl_find_status[i] == WlfStatus::Absent)
                                    || (a_w_line[i].nb_points() == 0)
                                {
                                    continue;
                                }
                                let plast = a_w_line[i].value(a_w_line[i].nb_points()).clone();
                                let (a_u1_c, _a_v1_c) =
                                    if is_reversed { (plast.u2, plast.v2) } else { (plast.u1, plast.v1) };
                                an_umax_added = an_umax_added.max(a_u1_c);
                                is_changed = true;
                            }
                            if !is_changed {
                                // If anUmaxAdded were not changed in the previous
                                // cycle then we would break existing WLines.
                                break;
                            }
                        }

                        for i in 0..a_nb_wlines {
                            if is_adding_wl_enabled[i] {
                                continue;
                            }
                            cyl_cyl_compute_parameters(
                                an_umax_added,
                                i as i32,
                                an_equation_coeffs,
                                &mut a_u2[i],
                                &mut a_v1[i],
                                &mut a_v2[i],
                            );
                            add_point_into_wl(
                                a_quad1,
                                a_quad2,
                                an_equation_coeffs,
                                is_reversed,
                                true,
                                [an_umax_added, a_v1[i]],
                                [a_u2[i], a_v2[i]],
                                a_usurf1f,
                                a_usurf1l,
                                a_usurf2f,
                                a_usurf2l,
                                a_vsurf1f,
                                a_vsurf1l,
                                a_vsurf2f,
                                a_vsurf2l,
                                a_period,
                                &mut a_w_line[i],
                                i,
                                a_tol_3d,
                                a_tol_2d,
                                false,
                                false,
                            );
                        }
                    }
                    break;
                }

                // Step computing.
                {
                    // Step of the aU1-parameter is computed adaptively.  The
                    // algorithm aims to provide given aDeltaV1 and aDeltaV2
                    // values (if possible because the intersection line can go
                    // along a V-isoline) in every iteration.
                    let a_delta_v1 = a_range_s1.delta() / a_nb_points as f64;
                    let a_delta_v2 = a_range_s2.delta() / a_nb_points as f64;

                    let mut a_matr = Mat35::new();
                    let mut a_min_uexp = super::cycy_common::real_last();
                    for i in 0..a_nb_wlines {
                        if a_tol_2d < (an_u_expect[i] - an_u1) {
                            continue;
                        }
                        if a_wl_find_status[i] == WlfStatus::Absent {
                            an_u_expect[i] += a_step_max;
                            a_min_uexp = a_min_uexp.min(an_u_expect[i]);
                            continue;
                        }
                        if is_good_intersection {
                            // Use constant step.
                            an_u_expect[i] += a_step_max;
                            a_min_uexp = a_min_uexp.min(an_u_expect[i]);
                            continue;
                        }

                        let mut a_step_tmp = a_step_max;

                        let (a_sin_u1, a_cos_u1) = an_u1.sin_cos();
                        let (a_sin_u2, a_cos_u2) = a_u2[i].sin_cos();

                        a_matr.set_col(1, an_equation_coeffs.m_vec_c1.to_array());
                        a_matr.set_col(2, an_equation_coeffs.m_vec_c2.to_array());
                        a_matr.set_col(
                            3,
                            (an_equation_coeffs.m_vec_a1 * a_sin_u1 - an_equation_coeffs.m_vec_b1 * a_cos_u1)
                                .to_array(),
                        );
                        a_matr.set_col(
                            4,
                            (an_equation_coeffs.m_vec_a2 * a_sin_u2 - an_equation_coeffs.m_vec_b2 * a_cos_u2)
                                .to_array(),
                        );
                        a_matr.set_col(
                            5,
                            (an_equation_coeffs.m_vec_a1 * a_cos_u1
                                + an_equation_coeffs.m_vec_b1 * a_sin_u1
                                + an_equation_coeffs.m_vec_a2 * a_cos_u2
                                + an_equation_coeffs.m_vec_b2 * a_sin_u2
                                + an_equation_coeffs.m_vec_d)
                                .to_array(),
                        );

                        // The main idea is to solve the linearized system (2)
                        // (see description to ComputationMethods class) in order
                        // to find the new U1-value to provide the new value V1
                        // or V2, which differs from the current one by aDeltaV1
                        // or aDeltaV2 respectively.  While linearizing, the
                        // following Taylor formulas are used:
                        //     cos(x0+dx) = cos(x0) - sin(x0)*dx
                        //     sin(x0+dx) = sin(x0) + cos(x0)*dx

                        if !super::cycy_boundaries::step_computing(
                            &a_matr,
                            a_v1[i],
                            a_v2[i],
                            a_delta_v1,
                            a_delta_v2,
                            &mut a_step_tmp,
                        ) {
                            // To avoid cycling-up.
                            an_u_expect[i] += a_step_max;
                            a_min_uexp = a_min_uexp.min(an_u_expect[i]);
                            continue;
                        }

                        if a_step_tmp < a_step_min {
                            a_step_tmp = a_step_min;
                        }
                        if a_step_tmp > a_step_max {
                            a_step_tmp = a_step_max;
                        }

                        an_u_expect[i] = an_u1 + a_step_tmp;
                        a_min_uexp = a_min_uexp.min(an_u_expect[i]);
                    }

                    an_u1_prev = an_u1;
                    an_u1 = a_min_uexp;
                }

                if PCONFUSION >= (an_ul - an_u1) {
                    an_u1 = an_ul;
                }

                an_uf = an_u1;

                for i in 0..a_nb_wlines {
                    if a_w_line[i].nb_points() != 1 {
                        is_added_into_wl[i] = false;
                    }
                    if an_u1 == an_ul {
                        // Strictly equal.  Tolerance is considered above.
                        an_u_expect[i] = an_ul;
                    }
                }
            }

            for i in 0..a_nb_wlines {
                if (a_w_line[i].nb_points() == 1) && (!is_added_into_wl[i]) {
                    *is_empty = false;
                    let p1 = a_w_line[i].value(1);
                    let (u1, v1, u2, v2) = (p1.u1, p1.v1, p1.u2, p1.v2);
                    // OCCT IntPatch_Point::SetParameter(u1) — the parameter on
                    // the (nonexistent) line; rcad IntPatchPoint drops it.
                    let a_p = IntPatchPoint {
                        p1: p1.p3d,
                        p2: p1.p3d,
                        u1,
                        v1,
                        u2,
                        v2,
                        tolerance: a_tol_3d,
                    };

                    // Check whether the added point exists.  It is enough to
                    // check the last point.
                    let same_as_last = match spnt.last() {
                        Some(last) => pnt2s_is_same(last, &a_p, rcad_kernel::precision::CONFUSION),
                        None => false,
                    };
                    if spnt.is_empty() || !same_as_last {
                        spnt.push(a_p);
                    }
                } else if a_w_line[i].nb_points() > 1 {
                    let mut is_good = true;
                    if a_w_line[i].nb_points() == 2 {
                        let a_pf = a_w_line[i].value(1).clone();
                        let a_pl = a_w_line[i].value(2).clone();
                        if is_same(&a_pf, &a_pl, rcad_kernel::precision::CONFUSION, -1.0) {
                            is_good = false;
                        }
                    } else if a_w_line[i].nb_points() > 2 {
                        // Sometimes points of the WLine are distributed linearly
                        // and uniformly.  However, such a position of the points
                        // does not always describe the real intersection curve:
                        // the real tangents at the ends of the intersection
                        // curve can significantly deviate from this "line"
                        // direction.  Here we process this case by inserting
                        // additional points at the beginning/end of the WLine to
                        // make it more precise (see issue #30082).
                        let a_sq_tol_3d = a_tol_3d * a_tol_3d;
                        for j in 0..2 {
                            // If j == 0 ==> add point at the begin of the WLine.
                            // If j == 1 ==> add point at the end of the WLine.
                            loop {
                                if a_w_line[i].nb_points() >= a_nb_max_points {
                                    break;
                                }
                                // Take 1st and 2nd point to compute the "line"
                                // direction.  For our convenience, make the 2nd
                                // point be the end of the WLine because it will
                                // be used for computation of the normals to the
                                // surfaces.
                                let an_idx1 = if j != 0 { a_w_line[i].nb_points() - 1 } else { 2 };
                                let an_idx2 = if j != 0 { a_w_line[i].nb_points() } else { 1 };
                                let a_p1 = a_w_line[i].value(an_idx1).p3d;
                                let a_p2 = a_w_line[i].value(an_idx2).p3d;
                                let a_dir = a_p2 - a_p1;
                                if a_dir.length_squared() < a_sq_tol_3d {
                                    break;
                                }
                                // Compute tangent in the first/last point of the
                                // WLine.  The flag "isReversed" is not taken into
                                // account because the strict direction of the
                                // tangent is not important here (we are
                                // interested in the tangent line itself).
                                let a_n1 = a_quad1.normale(a_p2);
                                let a_n2 = a_quad2.normale(a_p2);
                                let a_tg = a_n1.cross(a_n2);
                                if a_tg.length_squared() < SQUARE_CONFUSION {
                                    // Tangent zone.
                                    break;
                                }
                                // Check of the bending.
                                let mut an_angle = a_dir.angle_between(a_tg);
                                if an_angle > std::f64::consts::FRAC_PI_2 {
                                    an_angle -= std::f64::consts::PI;
                                }
                                if an_angle.abs() > 0.25 {
                                    // ~ 14 deg.
                                    let a_nb_pnts_prev = a_w_line[i].nb_points();
                                    seek_additional_points(
                                        a_quad1,
                                        a_quad2,
                                        &mut a_w_line[i],
                                        an_equation_coeffs,
                                        i,
                                        3,
                                        an_idx1,
                                        an_idx2,
                                        a_tol_2d,
                                        a_period,
                                        is_reversed,
                                    );
                                    if a_w_line[i].nb_points() == a_nb_pnts_prev {
                                        // No points have been added.  Exit from the loop.
                                        break;
                                    }
                                } else {
                                    // Good result has been achieved.  Exit from the loop.
                                    break;
                                }
                            }
                        }
                    }

                    if is_good {
                        *is_empty = false;
                        is_added_into_wl[i] = true;
                        let nb_line = a_w_line[i].nb_points();
                        seek_additional_points(
                            a_quad1,
                            a_quad2,
                            &mut a_w_line[i],
                            an_equation_coeffs,
                            i,
                            a_nb_points,
                            1,
                            nb_line,
                            a_tol_2d,
                            a_period,
                            is_reversed,
                        );
                        // OCCT L7622: aWLine[i]->ComputeVertexParameters(aTol3D).
                        // The IntCyCy WLines are created without SetPeriod, so
                        // the period array is all zeros (IntPatch_WLine.cxx
                        // L44-61); the vertex list (start/end + seam copies) is
                        // what GeomInt_LineConstructor's WLine path iterates.
                        let mut a_wl = IntPatchLine::walking(
                            std::mem::take(&mut a_w_line[i]).into_points(),
                            WLineType::ImpImp,
                        );
                        a_wl.compute_vertex_parameters_wline(a_tol_3d, [0.0, 0.0, 0.0, 0.0]);
                        slin.push(a_wl);
                    }
                } else {
                    is_added_into_wl[i] = false;
                }
            }
        }
    }

    // Delete the points in theSPnt which lie in at least one of the lines in
    // theSlin.
    let mut a_nb_pnt = 0usize;
    while a_nb_pnt < spnt.len() {
        for a_nb_lin in 0..slin.len() {
            let a_w_line1 = &slin[a_nb_lin];
            let Some(points) = (if a_w_line1.wline_pnts.is_empty() { None } else { Some(&a_w_line1.wline_pnts) })
            else {
                continue;
            };
            let a_pnt_f_wl1 = &points[0];
            let a_pnt_l_wl1 = &points[points.len() - 1];
            let a_pnt_cur_p = &spnt[a_nb_pnt];
            let a_pnt_cur = WLinePnt {
                p3d: a_pnt_cur_p.p1,
                u1: a_pnt_cur_p.u1,
                v1: a_pnt_cur_p.v1,
                u2: a_pnt_cur_p.u2,
                v2: a_pnt_cur_p.v2,
            };
            if is_same(&a_pnt_cur, a_pnt_f_wl1, a_tol_3d, -1.0)
                || is_same(&a_pnt_cur, a_pnt_l_wl1, a_tol_3d, -1.0)
            {
                spnt.remove(a_nb_pnt);
                a_nb_pnt = a_nb_pnt.wrapping_sub(1);
                break;
            }
        }
        a_nb_pnt += 1;
    }

    // Try to add new points in the neighborhood of an existing point.
    let mut a_nb_pnt = 0usize;
    while a_nb_pnt < spnt.len() {
        // The standard algorithm (implemented above) could not find any
        // continuous curve in the neighborhood of aPnt2S (e.g. because this
        // curve is too small; see tests/bugs/modalg_5/bug25292_35 and _36).
        // Here we try to find several new points nearer to aPnt2S.  The
        // algorithm below tries to find two points in every interval
        // [u1 - aStepMax, u1] and [u1, u1 + aStepMax]; every new point will be
        // at the maximal distance from u1.  If these two points exist they will
        // be joined by the intersection curve.
        let a_pnt2_s = spnt[a_nb_pnt].clone();
        let (u1, _v1, u2, _v2) = (a_pnt2_s.u1, a_pnt2_s.v1, a_pnt2_s.u2, a_pnt2_s.v2);

        let mut a_w_line = WLine::new();

        // Define the index of the WLine which the point aPnt2S lies in.
        let mut an_index = 0usize;

        let (mut an_uf, mut an_ul, a_cur_u2) = if is_reversed {
            ((u2 - a_step_max).max(a_usurf1f), (u2 + a_step_max).min(a_usurf1l), u1)
        } else {
            ((u1 - a_step_max).max(a_usurf1f), (u1 + a_step_max).min(a_usurf1l), u2)
        };

        let an_uinf = an_uf;
        let an_usup = an_ul;
        let an_umid = 0.5 * (an_uf + an_ul);

        {
            // Find the value of the anIndex variable.
            let mut a_delta = real_first();
            for i in 0..a_nb_wlines {
                let mut an_u2t = 0.0;
                if !cyl_cyl_compute_parameters_u2(
                    an_umid,
                    i as i32,
                    an_equation_coeffs,
                    &mut an_u2t,
                    None,
                ) {
                    continue;
                }
                let mut a_du2 = (an_u2t - a_cur_u2).abs() % a_period;
                a_du2 = a_du2.min((a_du2 - a_period).abs());
                if a_du2 < a_delta {
                    a_delta = a_du2;
                    an_index = i;
                }
            }
        }

        // Bisection method is used in order to find every new point.
        let mut an_added_par = [if is_reversed { u2 } else { u1 }, if is_reversed { u2 } else { u1 }];

        for a_par_id in 0..2 {
            if a_par_id == 0 {
                an_uf = an_uinf;
                an_ul = an_umid;
            } else {
                an_uf = an_umid;
                an_ul = an_usup;
            }

            while (an_ul - an_uf).abs() > a_step_min {
                let an_uc = 0.5 * (an_uf + an_ul);
                let mut a_u2 = 0.0;
                let mut a_v1 = 0.0;
                let mut a_v2 = 0.0;
                let mut is_done = cyl_cyl_compute_parameters(
                    an_uc,
                    an_index as i32,
                    an_equation_coeffs,
                    &mut a_u2,
                    &mut a_v1,
                    &mut a_v2,
                );

                if is_done {
                    if (a_v1 - a_vsurf1f).abs() <= a_tol_2d {
                        a_v1 = a_vsurf1f;
                    }
                    if (a_v1 - a_vsurf1l).abs() <= a_tol_2d {
                        a_v1 = a_vsurf1l;
                    }
                    if (a_v2 - a_vsurf2f).abs() <= a_tol_2d {
                        a_v2 = a_vsurf2f;
                    }
                    if (a_v2 - a_vsurf2l).abs() <= a_tol_2d {
                        a_v2 = a_vsurf2l;
                    }
                    is_done = add_point_into_wl(
                        a_quad1,
                        a_quad2,
                        an_equation_coeffs,
                        is_reversed,
                        true,
                        [an_uc, a_v1],
                        [a_u2, a_v2],
                        a_usurf1f,
                        a_usurf1l,
                        a_usurf2f,
                        a_usurf2l,
                        a_vsurf1f,
                        a_vsurf1l,
                        a_vsurf2f,
                        a_vsurf2l,
                        a_period,
                        &mut a_w_line,
                        an_index,
                        a_tol_3d,
                        a_tol_2d,
                        false,
                        true,
                    );
                }

                if is_done {
                    an_added_par[0] = an_added_par[0].min(an_uc);
                    an_added_par[1] = an_added_par[1].max(an_uc);
                    if a_par_id == 0 {
                        an_ul = an_uc;
                    } else {
                        an_uf = an_uc;
                    }
                } else {
                    if a_par_id == 0 {
                        an_uf = an_uc;
                    } else {
                        an_ul = an_uc;
                    }
                }
            }
        }

        // Fill aWLine with additional points.
        if an_added_par[1] - an_added_par[0] > a_step_min {
            for a_par_id in 0..2 {
                let mut a_u2 = 0.0;
                let mut a_v1 = 0.0;
                let mut a_v2 = 0.0;
                cyl_cyl_compute_parameters(
                    an_added_par[a_par_id],
                    an_index as i32,
                    an_equation_coeffs,
                    &mut a_u2,
                    &mut a_v1,
                    &mut a_v2,
                );
                add_point_into_wl(
                    a_quad1,
                    a_quad2,
                    an_equation_coeffs,
                    is_reversed,
                    true,
                    [an_added_par[a_par_id], a_v1],
                    [a_u2, a_v2],
                    a_usurf1f,
                    a_usurf1l,
                    a_usurf2f,
                    a_usurf2l,
                    a_vsurf1f,
                    a_vsurf1l,
                    a_vsurf2f,
                    a_vsurf2l,
                    a_period,
                    &mut a_w_line,
                    an_index,
                    a_tol_3d,
                    a_tol_2d,
                    false,
                    false,
                );
            }

            let nb_line = a_w_line.nb_points();
            seek_additional_points(
                a_quad1,
                a_quad2,
                &mut a_w_line,
                an_equation_coeffs,
                an_index,
                a_nb_min_points,
                1,
                nb_line,
                a_tol_2d,
                a_period,
                is_reversed,
            );
            // OCCT L7868: aWLine->ComputeVertexParameters(aTol3D) — see the
            // comment at the other flush point for the zero period array.
            let mut a_wl = IntPatchLine::walking(a_w_line.into_points(), WLineType::ImpImp);
            a_wl.compute_vertex_parameters_wline(a_tol_3d, [0.0, 0.0, 0.0, 0.0]);
            slin.push(a_wl);
            spnt.remove(a_nb_pnt);
            a_nb_pnt = a_nb_pnt.wrapping_sub(1);
        }
        a_nb_pnt += 1;
    }

    IntStatus::OK
}
