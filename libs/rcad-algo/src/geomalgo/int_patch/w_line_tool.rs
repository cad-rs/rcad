//! OCCT IntPatch_WLineTool — joining and extending Walking-lines
//! (IntPatch_WLineTool.cxx).
//!
//! 1:1 Rust translation of JoinWLines (L1605-1838) and ExtendTwoWLines
//! (L1874-2153) plus their static helpers.
//!
//! rcad data-model notes:
//!   - IntPatch_WLine -> IntPatchLine { line_type: Walking, wline_pnts, ... }.
//!   - IntSurf_PntOn2S -> WLinePnt { p3d, u1, v1, u2, v2 }.
//!   - Bnd_Box2d -> [f64; 4] = [u_min, v_min, u_max, v_max].

use super::{IntPatchLine, WLinePnt};
use glam::DVec3;
use rcad_kernel::geom::{Surface3, SurfaceEval};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};

/// OCCT IntPatch_WLineTool::myMaxConcatAngle = M_PI / 6.
const MY_MAX_CONCAT_ANGLE: f64 = std::f64::consts::FRAC_PI_6;

// OCCT L32-39: check-result bit-mask.
const WT_EN_ALL: u32 = 0x00;
const WT_DIS_LAST_LAST: u32 = 0x01;
const WT_DIS_LAST_FIRST: u32 = 0x02;
const WT_DIS_FIRST_LAST: u32 = 0x04;
const WT_DIS_FIRST_FIRST: u32 = 0x08;

/// OCCT IntPatchWT_WLsConnectionType (L41-47).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WlsConnectionType {
    NotConnected,
    Singular,
    Common,
    ReqExtend,
}

/// OCCT-aligned: IntPatch_WLineTool::MinMax (L56-64).
fn min_max(a: &mut f64, b: &mut f64) {
    if *a > *b {
        std::mem::swap(a, b);
    }
}

/// OCCT-aligned: Bnd_Range::IsIntersected (Bnd_Range.cxx L65-128).
/// Whether [r1, r2] (myFirst=r1, myLast=r2 — never void in
/// CheckArgumentsToExtend) contains p or a periodic copy of it.  Returns
/// truthy for In or Boundary, falsy for Out.
fn range_is_intersected(r1: f64, r2: f64, p: f64, period: f64) -> bool {
    let a_df = r1 - p;
    let a_dl = r2 - p;
    let a_period = period.abs();
    if a_period <= f64::MIN_POSITIVE {
        let a_delta = a_df * a_dl;
        if a_delta.abs() <= f64::EPSILON {
            return true; // myFirst or myLast lies ON the value: Boundary
        }
        return a_delta < 0.0; // In
    }
    let a_val1 = a_df / a_period;
    let a_val2 = a_dl / a_period;
    let a_par1 = a_val1.floor() as i64;
    let a_par2 = a_val2.floor() as i64;
    if a_par1 != a_par2 {
        return true; // In, or Boundary when myLast lies on the seam-edge
    }
    // a_par1 == a_par2: truthy only when myFirst lies on the seam-edge.
    (a_val1 - a_par1 as f64).abs() <= f64::EPSILON
}

/// OCCT IntPatch_WLineTool::JoinWLines (L1605-1838).
///
/// Joins consecutive Walking-lines that share an endpoint (for two cylindrical
/// surfaces) into a single line.
pub fn join_w_lines(slin: &mut Vec<IntPatchLine>, s1: &Surface3, s2: &Surface3, _tol_3d: f64) {
    if slin.is_empty() {
        return;
    }
    let (r1, r2) = match (s1, s2) {
        (Surface3::Cylinder(c1), Surface3::Cylinder(c2)) => (c1.radius, c2.radius),
        _ => (1.0, 1.0),
    };
    let a_min_rad = 1.0e-3 * r1.min(r2);

    let an_arr_periods = [period(s1, true), period(s1, false), period(s2, true), period(s2, false)];
    let an_arr_f_bonds = [first_u(s1), first_v(s1), first_u(s2), first_v(s2)];
    let an_arr_l_bonds = [last_u(s1), last_v(s1), last_u(s2), last_v(s2)];

    let mut a_n1 = 0usize;
    while a_n1 < slin.len() {
        let wl1 = slin[a_n1].clone();
        if wl1.line_type != super::IntPatchIType::Walking || wl1.wline_pnts.is_empty() {
            a_n1 += 1;
            continue;
        }
        let a_nb_pnts_wl1 = wl1.wline_pnts.len();
        let a_pnt_f_wl1 = wl1.wline_pnts[0];
        let a_pnt_l_wl1 = wl1.wline_pnts[a_nb_pnts_wl1 - 1];

        let mut list_fc: Vec<usize> = Vec::new();
        let mut list_lc: Vec<usize> = Vec::new();
        let mut is_first_connected = false;
        let mut is_last_connected = false;

        for a_n2 in 0..slin.len() {
            if a_n2 == a_n1 {
                continue;
            }
            let wl2 = slin[a_n2].clone();
            if wl2.line_type != super::IntPatchIType::Walking || wl2.wline_pnts.is_empty() {
                continue;
            }
            let a_nb_pnts_wl2 = wl2.wline_pnts.len();
            let a_pnt_f_wl2 = wl2.wline_pnts[0];
            let a_pnt_l_wl2 = wl2.wline_pnts[a_nb_pnts_wl2 - 1];

            is_first_connected = false;
            is_last_connected = false;

            let mut a_sq_dist_f = a_pnt_f_wl1.p3d.distance_squared(a_pnt_f_wl2.p3d);
            let mut a_sq_dist_l = a_pnt_f_wl1.p3d.distance_squared(a_pnt_l_wl2.p3d);
            let a_sq_min_f_dist = a_sq_dist_f.min(a_sq_dist_l);
            if a_sq_min_f_dist < CONFUSION * CONFUSION {
                let is_fm = a_sq_dist_f < a_sq_dist_l;
                let a_pt1 = wl1.wline_pnts[1];
                let a_pt2 = if is_fm { wl2.wline_pnts[1] } else { wl2.wline_pnts[a_nb_pnts_wl2 - 2] };
                if !is_seam_or_bound(&a_pt1, &a_pt2, &a_pnt_f_wl1, &an_arr_periods, &an_arr_f_bonds, &an_arr_l_bonds) {
                    is_first_connected = true;
                }
            }

            a_sq_dist_f = a_pnt_l_wl1.p3d.distance_squared(a_pnt_f_wl2.p3d);
            a_sq_dist_l = a_pnt_l_wl1.p3d.distance_squared(a_pnt_l_wl2.p3d);
            let a_sq_min_l_dist = a_sq_dist_f.min(a_sq_dist_l);
            if a_sq_min_l_dist < CONFUSION * CONFUSION {
                let is_fm = a_sq_dist_f < a_sq_dist_l;
                let a_pt1 = wl1.wline_pnts[a_nb_pnts_wl1 - 2];
                let a_pt2 = if is_fm { wl2.wline_pnts[1] } else { wl2.wline_pnts[a_nb_pnts_wl2 - 2] };
                if !is_seam_or_bound(&a_pt1, &a_pt2, &a_pnt_l_wl1, &an_arr_periods, &an_arr_f_bonds, &an_arr_l_bonds) {
                    is_last_connected = true;
                }
            }

            if is_first_connected && is_last_connected {
                if a_sq_min_f_dist < a_sq_min_l_dist {
                    list_fc.push(a_n2);
                } else {
                    list_lc.push(a_n2);
                }
            } else if is_first_connected {
                list_fc.push(a_n2);
            } else if is_last_connected {
                list_lc.push(a_n2);
            }
        }

        is_first_connected = list_fc.len() == 1;
        is_last_connected = list_lc.len() == 1;
        if !(is_first_connected || is_last_connected) {
            a_n1 += 1;
            continue;
        }

        let an_index_wl2 = if is_first_connected { list_fc[0] } else { list_lc[0] };
        let wl2 = slin[an_index_wl2].clone();
        let a_nb_pnts_wl2 = wl2.wline_pnts.len();
        let a_pnt_f_wl2 = wl2.wline_pnts[0];
        let a_pnt_l_wl2 = wl2.wline_pnts[a_nb_pnts_wl2 - 1];

        if is_first_connected {
            let a_sq_dist_f = a_pnt_f_wl1.p3d.distance_squared(a_pnt_f_wl2.p3d);
            let a_sq_dist_l = a_pnt_f_wl1.p3d.distance_squared(a_pnt_l_wl2.p3d);
            let is_fm = a_sq_dist_f < a_sq_dist_l;
            let a_pt1 = wl1.wline_pnts[1];
            let a_pt2 = if is_fm { wl2.wline_pnts[1] } else { wl2.wline_pnts[a_nb_pnts_wl2 - 2] };
            if !check_arguments_to_join(
                a_pnt_f_wl1.p3d,
                a_pt1.p3d,
                a_pnt_f_wl1.p3d,
                a_pt2.p3d,
                a_min_rad,
            ) {
                a_n1 += 1;
                continue;
            }
            // First-First or First-Last connection: prepend wl2's points.
            let mut joined: Vec<WLinePnt> = Vec::new();
            if is_fm {
                for pt in wl2.wline_pnts.iter() {
                    joined.push(*pt);
                }
            } else {
                for pt in wl2.wline_pnts.iter().rev() {
                    joined.push(*pt);
                }
            }
            joined.extend_from_slice(&wl1.wline_pnts);
            slin[a_n1].wline_pnts = joined;
        } else {
            let a_sq_dist_f = a_pnt_l_wl1.p3d.distance_squared(a_pnt_f_wl2.p3d);
            let a_sq_dist_l = a_pnt_l_wl1.p3d.distance_squared(a_pnt_l_wl2.p3d);
            let is_fm = a_sq_dist_f < a_sq_dist_l;
            let a_pt1 = wl1.wline_pnts[a_nb_pnts_wl1 - 2];
            let a_pt2 = if is_fm { wl2.wline_pnts[1] } else { wl2.wline_pnts[a_nb_pnts_wl2 - 2] };
            if !check_arguments_to_join(
                a_pnt_l_wl1.p3d,
                a_pt1.p3d,
                a_pnt_l_wl1.p3d,
                a_pt2.p3d,
                a_min_rad,
            ) {
                a_n1 += 1;
                continue;
            }
            // Last-First or Last-Last connection: append wl2's points.
            let mut joined = wl1.wline_pnts.clone();
            if is_fm {
                for pt in wl2.wline_pnts.iter() {
                    joined.push(*pt);
                }
            } else {
                for pt in wl2.wline_pnts.iter().rev() {
                    joined.push(*pt);
                }
            }
            slin[a_n1].wline_pnts = joined;
        }

        slin[a_n1].vertices.clear();
        slin.remove(an_index_wl2);
    }
}

/// OCCT static IsSeamOrBound: whether the line segment from aPt1 to aPt2 (both
/// near aPnt) lies on a seam or domain boundary.
#[allow(clippy::too_many_arguments)]
fn is_seam_or_bound(
    a_pt1: &WLinePnt,
    a_pt2: &WLinePnt,
    a_pnt: &WLinePnt,
    arr_periods: &[f64; 4],
    arr_f_bonds: &[f64; 4],
    arr_l_bonds: &[f64; 4],
) -> bool {
    let pars1 = [a_pt1.u1, a_pt1.v1, a_pt1.u2, a_pt1.v2];
    let pars2 = [a_pt2.u1, a_pt2.v1, a_pt2.u2, a_pt2.v2];
    let parsv = [a_pnt.u1, a_pnt.v1, a_pnt.u2, a_pnt.v2];
    for i in 0..4 {
        if arr_periods[i] == 0.0 {
            continue;
        }
        if (pars1[i] - parsv[i]).abs() > 0.5 * arr_periods[i]
            || (pars2[i] - parsv[i]).abs() > 0.5 * arr_periods[i]
        {
            return true;
        }
    }
    for i in 0..4 {
        if arr_periods[i] != 0.0 {
            continue;
        }
        if (pars1[i] - arr_f_bonds[i]).abs() < 1e-9 && (pars2[i] - arr_f_bonds[i]).abs() < 1e-9 {
            return true;
        }
        if (pars1[i] - arr_l_bonds[i]).abs() < 1e-9 && (pars2[i] - arr_l_bonds[i]).abs() < 1e-9 {
            return true;
        }
    }
    false
}

/// OCCT static CheckArgumentsToJoin (L1074-1119).
#[allow(clippy::too_many_arguments)]
fn check_arguments_to_join(
    the_pnt: DVec3,
    the_p1: DVec3,
    the_p2: DVec3,
    the_p3: DVec3,
    _the_min_rad: f64,
) -> bool {
    // OCCT computes the curvature radius of the intersection line; when it
    // exceeds theMinRad, joining is allowed.
    let _ = the_pnt;
    // rcad: the line's curvature radius is not computed here; fall back to the
    // polygon-smoothness check (OCCT L1108-1118).
    let a_v12f = the_p2 - the_p1;
    let a_v12l = the_p3 - the_p2;
    if a_v12f.angle_between(a_v12l) > MY_MAX_CONCAT_ANGLE {
        return false;
    }
    let a_v13 = the_p3 - the_p1;
    let a_sq13 = a_v13.length_squared();
    let a_cross = a_v12f.cross(a_v13);
    a_cross.length_squared() < 1.0e-4 * a_sq13 * a_sq13
}

/// OCCT-aligned: IntPatch_WLineTool::ExtendTwoWLines (L1874-2153).
///
/// For pairs of WLines that meet at an endpoint (within theToler3D), extends one
/// of them through the shared point so that the pair forms a single smooth
/// curve.  theBoxS1/theBoxS2 are the surface UV-domain rectangles; theListOf
/// critical points are the cone apexes / sphere poles.
#[allow(clippy::too_many_lines)]
pub fn extend_two_w_lines(
    slin: &mut Vec<IntPatchLine>,
    s1: &Surface3,
    s2: &Surface3,
    tol_3d: f64,
    arr_periods: &[f64; 4],
    box_s1: [f64; 4],
    box_s2: [f64; 4],
    list_of_critical_points: &[DVec3],
) {
    if slin.len() < 2 {
        return;
    }

    let mut has_been_joined_counter = 0usize;
    let mut a_num_of_line1 = 0usize;
    while a_num_of_line1 < slin.len() {
        if has_been_joined_counter > 0 {
            a_num_of_line1 -= 1;
        }
        has_been_joined_counter = 0;

        let wl1 = slin[a_num_of_line1].clone();
        if wl1.line_type != super::IntPatchIType::Walking || wl1.wline_pnts.is_empty() {
            a_num_of_line1 += 1;
            continue;
        }
        let a_nb_pnts_wl1 = wl1.wline_pnts.len();

        // OCCT L1910-1918: the first vertex must be at parameter 1 and the last
        // vertex at NbPnts.
        if !wl1_first_last_vertex_ok(&wl1) {
            a_num_of_line1 += 1;
            continue;
        }

        let a_pnt_f_wl1 = wl1.wline_pnts[0];
        let a_pnt_fp1_wl1 = wl1.wline_pnts[1];
        let a_pnt_l_wl1 = wl1.wline_pnts[a_nb_pnts_wl1 - 1];
        let a_pnt_lm1_wl1 = wl1.wline_pnts[a_nb_pnts_wl1 - 2];

        if is_need_skip_wl(&wl1, &box_s1, &box_s2, arr_periods) {
            a_num_of_line1 += 1;
            continue;
        }

        let mut a_check_result = WT_EN_ALL;

        // OCCT L1939-2009: build the check-result mask from the endpoint
        // coincidences of wl1 with the other lines and the critical points.
        for a_num_of_line2 in (a_num_of_line1 + 1)..slin.len() {
            let wl2 = slin[a_num_of_line2].clone();
            if wl2.line_type != super::IntPatchIType::Walking || wl2.wline_pnts.is_empty() {
                continue;
            }
            let a_pnt_f_wl2 = wl2.wline_pnts[0];
            let a_pnt_l_wl2 = wl2.wline_pnts[wl2.wline_pnts.len() - 1];

            if !(pnt_is_same(&a_pnt_f_wl1, &a_pnt_f_wl2, tol_3d, PCONFUSION)
                || pnt_is_same(&a_pnt_f_wl1, &a_pnt_l_wl2, tol_3d, PCONFUSION))
            {
                if pnt_is_same(&a_pnt_f_wl1, &a_pnt_f_wl2, tol_3d, -1.0)
                    || pnt_is_same(&a_pnt_f_wl1, &a_pnt_l_wl2, tol_3d, -1.0)
                {
                    a_check_result |= WT_DIS_FIRST_FIRST | WT_DIS_FIRST_LAST;
                }
            }

            if !(pnt_is_same(&a_pnt_l_wl1, &a_pnt_f_wl2, tol_3d, PCONFUSION)
                || pnt_is_same(&a_pnt_l_wl1, &a_pnt_l_wl2, tol_3d, PCONFUSION))
            {
                if pnt_is_same(&a_pnt_l_wl1, &a_pnt_f_wl2, tol_3d, -1.0)
                    || pnt_is_same(&a_pnt_l_wl1, &a_pnt_l_wl2, tol_3d, -1.0)
                {
                    a_check_result |= WT_DIS_LAST_FIRST | WT_DIS_LAST_LAST;
                }
            }

            for pt in list_of_critical_points.iter() {
                if (a_check_result & (WT_DIS_FIRST_FIRST | WT_DIS_FIRST_LAST)) == 0 {
                    if pt.distance_squared(a_pnt_f_wl1.p3d) < CONFUSION {
                        a_check_result |= WT_DIS_FIRST_FIRST | WT_DIS_FIRST_LAST;
                    }
                }
                if (a_check_result & (WT_DIS_LAST_FIRST | WT_DIS_LAST_LAST)) == 0 {
                    if pt.distance_squared(a_pnt_l_wl1.p3d) < CONFUSION {
                        a_check_result |= WT_DIS_LAST_FIRST | WT_DIS_LAST_LAST;
                    }
                }
                if (a_check_result & (WT_DIS_FIRST_FIRST | WT_DIS_LAST_FIRST)) == 0 {
                    if pt.distance_squared(a_pnt_f_wl2.p3d) < CONFUSION {
                        a_check_result |= WT_DIS_FIRST_FIRST | WT_DIS_LAST_FIRST;
                    }
                }
                if (a_check_result & (WT_DIS_FIRST_LAST | WT_DIS_LAST_LAST)) == 0 {
                    if pt.distance_squared(a_pnt_l_wl2.p3d) < CONFUSION {
                        a_check_result |= WT_DIS_FIRST_LAST | WT_DIS_LAST_LAST;
                    }
                }
            }
        }

        if a_check_result == (WT_DIS_FIRST_FIRST | WT_DIS_FIRST_LAST | WT_DIS_LAST_FIRST | WT_DIS_LAST_LAST) {
            a_num_of_line1 += 1;
            continue;
        }

        // OCCT L2018-2151: try the four extension directions.  OCCT re-evaluates
        // theSlin.Length() in the loop condition and decrements aNumOfLine2 after
        // a join (L2148-2149), so the line shifted into the removed slot is
        // re-checked.
        let mut a_num_of_line2 = a_num_of_line1 + 1;
        while a_num_of_line2 < slin.len() {
            let wl2 = slin[a_num_of_line2].clone();
            if wl2.line_type != super::IntPatchIType::Walking || wl2.wline_pnts.is_empty() {
                a_num_of_line2 += 1;
                continue;
            }
            if !wl1_first_last_vertex_ok(&wl2) {
                a_num_of_line2 += 1;
                continue;
            }
            if is_need_skip_wl(&wl2, &box_s1, &box_s2, arr_periods) {
                a_num_of_line2 += 1;
                continue;
            }

            let mut has_been_joined = false;
            let a_nb_pnts_wl2 = wl2.wline_pnts.len();
            let a_pnt_f_wl2 = wl2.wline_pnts[0];
            let a_pnt_fp1_wl2 = wl2.wline_pnts[1];
            let a_pnt_l_wl2 = wl2.wline_pnts[a_nb_pnts_wl2 - 1];
            let a_pnt_lm1_wl2 = wl2.wline_pnts[a_nb_pnts_wl2 - 2];

            let wl1_mut_ref = a_num_of_line1;
            let wl2_idx = a_num_of_line2;

            if (a_check_result & WT_DIS_FIRST_FIRST) == 0 {
                // First/First
                let a_vec1 = a_pnt_fp1_wl1.p3d - a_pnt_f_wl1.p3d;
                let a_vec2 = a_pnt_f_wl2.p3d - a_pnt_fp1_wl2.p3d;
                let a_vec3 = a_pnt_f_wl1.p3d - a_pnt_f_wl2.p3d;
                extend_two_wl(
                    s1, s2, wl1_mut_ref, wl2_idx, slin,
                    a_pnt_f_wl1, a_pnt_f_wl2, a_vec1, a_vec2, a_vec3,
                    &box_s1, &box_s2, tol_3d, arr_periods,
                    WT_DIS_FIRST_LAST | WT_DIS_LAST_FIRST,
                    &mut a_check_result, &mut has_been_joined,
                    ExtendMode::FirstFirst,
                );
            }

            if (a_check_result & WT_DIS_FIRST_LAST) == 0 {
                // First/Last
                let a_vec1 = a_pnt_fp1_wl1.p3d - a_pnt_f_wl1.p3d;
                let a_vec2 = a_pnt_l_wl2.p3d - a_pnt_lm1_wl2.p3d;
                let a_vec3 = a_pnt_f_wl1.p3d - a_pnt_l_wl2.p3d;
                extend_two_wl(
                    s1, s2, wl1_mut_ref, wl2_idx, slin,
                    a_pnt_f_wl1, a_pnt_l_wl2, a_vec1, a_vec2, a_vec3,
                    &box_s1, &box_s2, tol_3d, arr_periods,
                    WT_DIS_LAST_LAST,
                    &mut a_check_result, &mut has_been_joined,
                    ExtendMode::FirstLast,
                );
            }

            if (a_check_result & WT_DIS_LAST_FIRST) == 0 {
                // Last/First
                let a_vec1 = a_pnt_l_wl1.p3d - a_pnt_lm1_wl1.p3d;
                let a_vec2 = a_pnt_fp1_wl2.p3d - a_pnt_f_wl2.p3d;
                let a_vec3 = a_pnt_f_wl2.p3d - a_pnt_l_wl1.p3d;
                extend_two_wl(
                    s1, s2, wl1_mut_ref, wl2_idx, slin,
                    a_pnt_l_wl1, a_pnt_f_wl2, a_vec1, a_vec2, a_vec3,
                    &box_s1, &box_s2, tol_3d, arr_periods,
                    WT_DIS_LAST_LAST,
                    &mut a_check_result, &mut has_been_joined,
                    ExtendMode::LastFirst,
                );
            }

            if (a_check_result & WT_DIS_LAST_LAST) == 0 {
                // Last/Last
                let a_vec1 = a_pnt_l_wl1.p3d - a_pnt_lm1_wl1.p3d;
                let a_vec2 = a_pnt_lm1_wl2.p3d - a_pnt_l_wl2.p3d;
                let a_vec3 = a_pnt_l_wl2.p3d - a_pnt_l_wl1.p3d;
                extend_two_wl(
                    s1, s2, wl1_mut_ref, wl2_idx, slin,
                    a_pnt_l_wl1, a_pnt_l_wl2, a_vec1, a_vec2, a_vec3,
                    &box_s1, &box_s2, tol_3d, arr_periods,
                    WT_DIS_LAST_LAST,
                    &mut a_check_result, &mut has_been_joined,
                    ExtendMode::LastLast,
                );
            }

            if has_been_joined {
                has_been_joined_counter += 1;
                slin.remove(wl2_idx);
                // OCCT L2148-2149: aNumOfLine2-- — stay on the same index so the
                // line shifted into the removed slot is re-checked.
                continue;
            }
            a_num_of_line2 += 1;
        }
        a_num_of_line1 += 1;
    }
}

/// Which extension mode the four ExtendTwoWL* functions implement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExtendMode {
    FirstFirst,
    FirstLast,
    LastFirst,
    LastLast,
}

/// OCCT-aligned: IntPatch_WLineTool::ExtendTwoWLFirstFirst/FirstLast/
/// LastFirst/LastLast (L1126-1464): shared by the four modes; extends wl1
/// through its start/end point and joins wl2.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn extend_two_wl(
    s1: &Surface3,
    s2: &Surface3,
    wl1_idx: usize,
    wl2_idx: usize,
    slin: &mut Vec<IntPatchLine>,
    pt_wl1: WLinePnt,
    pt_wl2: WLinePnt,
    vec1: DVec3,
    vec2: DVec3,
    vec3: DVec3,
    box_s1: &[f64; 4],
    box_s2: &[f64; 4],
    tol_3d: f64,
    arr_periods: &[f64; 4],
    check_mask: u32,
    check_result: &mut u32,
    has_been_joined: &mut bool,
    mode: ExtendMode,
) {
    let mut a_p_on2s = WLinePnt {
        p3d: DVec3::ZERO,
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
    };
    let a_check_res = check_arguments_to_extend(
        s1, s2, &pt_wl1, &pt_wl2, &mut a_p_on2s, vec1, vec2, vec3,
        box_s1, box_s2, tol_3d, arr_periods,
    );

    if a_check_res != WlsConnectionType::NotConnected {
        *check_result |= check_mask;
    } else {
        return;
    }

    // OCCT: AdjustPointAndVertex(ptWL1, periods, aPOn2S); ExtendFirst/Last(wl1);
    //       AdjustPointAndVertex(ptWL2, periods, aPOn2S); ExtendFirst/Last(wl2).
    match mode {
        ExtendMode::FirstFirst => {
            adjust_uv(&pt_wl1, arr_periods, &mut a_p_on2s);
            extend_first(slin, wl1_idx, a_p_on2s);
            adjust_uv(&pt_wl2, arr_periods, &mut a_p_on2s);
            extend_first(slin, wl2_idx, a_p_on2s);
        }
        ExtendMode::FirstLast => {
            adjust_uv(&pt_wl1, arr_periods, &mut a_p_on2s);
            extend_first(slin, wl1_idx, a_p_on2s);
            adjust_uv(&pt_wl2, arr_periods, &mut a_p_on2s);
            extend_last(slin, wl2_idx, a_p_on2s);
        }
        ExtendMode::LastFirst => {
            adjust_uv(&pt_wl1, arr_periods, &mut a_p_on2s);
            extend_last(slin, wl1_idx, a_p_on2s);
            adjust_uv(&pt_wl2, arr_periods, &mut a_p_on2s);
            extend_first(slin, wl2_idx, a_p_on2s);
        }
        ExtendMode::LastLast => {
            adjust_uv(&pt_wl1, arr_periods, &mut a_p_on2s);
            extend_last(slin, wl1_idx, a_p_on2s);
            adjust_uv(&pt_wl2, arr_periods, &mut a_p_on2s);
            extend_last(slin, wl2_idx, a_p_on2s);
        }
    }

    if *has_been_joined || a_check_res == WlsConnectionType::Singular {
        return;
    }

    // OCCT L1175-1209 (FirstFirst): remove the duplicated first vertices, join
    // wl2's points to the front of wl1, re-index the vertices.
    // The other modes follow the same pattern with the respective ends.
    let (wl1_n, wl2_n) = {
        let n1 = slin[wl1_idx].wline_pnts.len();
        let n2 = slin[wl2_idx].wline_pnts.len();
        (n1, n2)
    };
    // Remove the first/last vertex (index 1 or N) of wl1 and wl2 as OCCT does.
    let (rm_wl1_front, rm_wl1_back, rm_wl2_front, rm_wl2_back) = match mode {
        ExtendMode::FirstFirst => (true, false, true, false),
        ExtendMode::FirstLast => (true, false, false, true),
        ExtendMode::LastFirst => (false, true, true, false),
        ExtendMode::LastLast => (false, true, false, true),
    };
    remove_endpoint_vertex(slin, wl1_idx, rm_wl1_front, rm_wl1_back);
    remove_endpoint_vertex(slin, wl2_idx, rm_wl2_front, rm_wl2_back);

    // Join the points.
    let (p1, p2) = (wl1_n, wl2_n);
    match mode {
        ExtendMode::FirstFirst => {
            // OCCT L1187-1192: for (aNPt = 2; aNPt <= aNbPts; aNPt++)
            //   theWLine1->Curve()->InsertBefore(1, theWLine2->Point(aNPt)).
            // Each InsertBefore(1) shifts right, so wl2's points 2..N land
            // REVERSED at the front: [wl2_N, ..., wl2_2, wl1_1, ..., wl1_N1].
            let mut joined: Vec<WLinePnt> = Vec::new();
            for pt in slin[wl2_idx].wline_pnts.iter().skip(1).rev() {
                joined.push(*pt);
            }
            for pt in slin[wl1_idx].wline_pnts.iter() {
                joined.push(*pt);
            }
            // OCCT L1194-1199: wl1 vertices stay in place, shifted by +(N2-1).
            for v in slin[wl1_idx].vertices.iter_mut() {
                v.param_on_line += (p2 - 1) as f64;
            }
            // OCCT L1201-1207: for (aNVtx = NbVertex; aNVtx >= 1; aNVtx--)
            //   SetParameter(N2 - cur + 1); AddVertex(v, true) = svtx.Prepend.
            // Iterating NbVertex..1 and prepending each leaves the wl2 vertices
            // at the front in their ORIGINAL order (last iteration = Vertex(1)
            // lands first).
            let wl2_verts = std::mem::take(&mut slin[wl2_idx].vertices);
            for mut v in wl2_verts.into_iter().rev() {
                v.param_on_line = p2 as f64 - v.param_on_line + 1.0;
                slin[wl1_idx].vertices.insert(0, v);
            }
            slin[wl1_idx].wline_pnts = joined;
        }
        ExtendMode::FirstLast => {
            // OCCT L1278-1283: for (aNPt = aNbPts - 1; aNPt >= 1; aNPt--)
            //   theWLine1->Curve()->InsertBefore(1, theWLine2->Point(aNPt)).
            // Each InsertBefore(1) shifts right, so the last insertion (Point(1))
            // lands first — net ASCENDING [wl2_1, ..., wl2_{N2-1}] at the front
            // (the connection point aNbPts is dropped).
            let mut joined: Vec<WLinePnt> = Vec::new();
            for pt in slin[wl2_idx].wline_pnts.iter().take(wl2_n - 1) {
                joined.push(*pt);
            }
            for pt in slin[wl1_idx].wline_pnts.iter() {
                joined.push(*pt);
            }
            // OCCT L1285-1290: wl1 vertices stay in place, shifted by +(N2-1).
            for v in slin[wl1_idx].vertices.iter_mut() {
                v.param_on_line += (p2 - 1) as f64;
            }
            // OCCT L1292-1296: wl2 vertices iterated aNVtx = NbVertex..1 and
            // prepended WITHOUT re-indexing their parameters.
            let wl2_verts = std::mem::take(&mut slin[wl2_idx].vertices);
            for v in wl2_verts.into_iter().rev() {
                slin[wl1_idx].vertices.insert(0, v);
            }
            slin[wl1_idx].wline_pnts = joined;
        }
        ExtendMode::LastFirst => {
            // wl2's points (2..N) appended to wl1.
            let mut joined = slin[wl1_idx].wline_pnts.clone();
            for pt in slin[wl2_idx].wline_pnts.iter().skip(1) {
                joined.push(*pt);
            }
            let wl2_verts = std::mem::take(&mut slin[wl2_idx].vertices);
            for mut v in wl2_verts {
                v.param_on_line = p1 as f64 + v.param_on_line - 1.0;
                slin[wl1_idx].vertices.push(v);
            }
            slin[wl1_idx].wline_pnts = joined;
        }
        ExtendMode::LastLast => {
            // OCCT L1448-1453: for (aNPt = NbPnts(wl2) - 1; aNPt >= 1; aNPt--)
            //   theWLine1->Curve()->Add(theWLine2->Point(aNPt)).
            let mut joined = slin[wl1_idx].wline_pnts.clone();
            for pt in slin[wl2_idx].wline_pnts.iter().rev().skip(1) {
                joined.push(*pt);
            }
            // OCCT L1455-1461: wl2 vertices iterated NbVertex..1, SetParameter
            // (NbPnts1+NbPnts2 - cur), appended — reverse iteration yields the
            // wl2 vertices reversed at the end.
            let wl2_verts = std::mem::take(&mut slin[wl2_idx].vertices);
            for mut v in wl2_verts.into_iter().rev() {
                v.param_on_line = p1 as f64 + p2 as f64 - v.param_on_line;
                slin[wl1_idx].vertices.push(v);
            }
            slin[wl1_idx].wline_pnts = joined;
        }
    }
    // wl2 is removed by the caller.
    let _ = (wl1_n, wl2_n, p1, p2);
    *has_been_joined = true;
}

/// OCCT-aligned: IntPatch_WLineTool::ExtendTwoWLines L1910-1918 (exact
/// integer param checks): Vertex(1).ParameterOnLine() == 1 and
/// Vertex(NbVertex).ParameterOnLine() == NbPnts.
fn wl1_first_last_vertex_ok(wl: &IntPatchLine) -> bool {
    if wl.vertices.is_empty() || wl.wline_pnts.is_empty() {
        return false;
    }
    // OCCT L1910-1918: Vertex(1).ParameterOnLine() == 1 and
    // Vertex(NbVertex).ParameterOnLine() == NbPnts (exact integer params).
    wl.vertices[0].param_on_line == 1.0
        && wl.vertices[wl.vertices.len() - 1].param_on_line == wl.wline_pnts.len() as f64
}

/// OCCT-aligned: IntPatch_WLineTool::ExtendTwoWL* L1175-1185 etc (the
/// while (Vertex(1)/Vertex(NbVertex).ParameterOnLine() == aPrm) RemoveVertex
/// loops).  Remove the first and/or last endpoint vertex run.
fn remove_endpoint_vertex(slin: &mut Vec<IntPatchLine>, idx: usize, front: bool, back: bool) {
    let verts = &mut slin[idx].vertices;
    // OCCT L1175-1185 / L1266-1276 / L1355-1365 / L1436-1446: removes the
    // contiguous run at the front (all sharing the first vertex's parameter) or
    // at the back (all sharing the last vertex's parameter), matching the
    // while (Vertex(1)/Vertex(NbVertex).ParameterOnLine() == aPrm) RemoveVertex
    // loops.  The run is bounded by the FIRST endpoint value only — vertices
    // further along with the same parameter are kept.
    if front {
        if let Some(prm) = verts.first().map(|v| v.param_on_line) {
            let mut k = 0;
            // OCCT: while (Vertex(1).ParameterOnLine() == aPrm) RemoveVertex(1);
            while k < verts.len() && verts[k].param_on_line == prm {
                k += 1;
            }
            verts.drain(..k);
        }
    }
    if back {
        if let Some(prm) = verts.last().map(|v| v.param_on_line) {
            let mut k = verts.len();
            // OCCT: while (Vertex(NbVertex).ParameterOnLine() == aPrm)
            //   RemoveVertex(NbVertex);
            while k > 0 && verts[k - 1].param_on_line == prm {
                k -= 1;
            }
            verts.truncate(k);
        }
    }
}

/// OCCT-aligned: IntPatch_WLineTool::IsNeedSkipWL (L1842-1867).
fn is_need_skip_wl(wl: &IntPatchLine, box_s1: &[f64; 4], box_s2: &[f64; 4], arr_periods: &[f64; 4]) -> bool {
    let a_nb_vtx = wl.vertices.len();
    // OCCT L1850-1860: Vertex params are 1-based point indices; Point(pmid)
    // indexes wline_pnts at pmid-1 (rcad Vec is 0-based).  No bounds guards in
    // OCCT — the callers guarantee Vertex(1)@1 .. Vertex(NbVertex)@NbPnts.
    for i in 0..a_nb_vtx.saturating_sub(1) {
        let a_firstp = wl.vertices[i].param_on_line;
        let a_lastp = wl.vertices[i + 1].param_on_line;
        let pmid = ((a_firstp + a_lastp) / 2.0) as usize;
        let a_pmid = &wl.wline_pnts[pmid - 1];
        if is_out_of_domain(box_s1, box_s2, a_pmid, arr_periods) {
            return true;
        }
    }
    false
}

/// OCCT-aligned: IntPatch_WLineTool::IsOutOfDomain (L891-911).
fn is_out_of_domain(box_s1: &[f64; 4], box_s2: &[f64; 4], p_on2s: &WLinePnt, arr_periods: &[f64; 4]) -> bool {
    let a_u1 = in_period(p_on2s.u1, box_s1[0], box_s1[0] + arr_periods[0]);
    let a_v1 = in_period(p_on2s.v1, box_s1[1], box_s1[1] + arr_periods[1]);
    let a_u2 = in_period(p_on2s.u2, box_s2[0], box_s2[0] + arr_periods[2]);
    let a_v2 = in_period(p_on2s.v2, box_s2[1], box_s2[1] + arr_periods[3]);
    is_out_rect(a_u1, a_v1, box_s1) || is_out_rect(a_u2, a_v2, box_s2)
}

/// OCCT-aligned: ElCLib::InPeriod (TKMath ElCLib.cxx L424-448) — the point
/// in the period interval.
fn in_period(par: f64, min: f64, max: f64) -> f64 {
    let period = max - min;
    if period <= 0.0 || period.is_infinite() {
        return par;
    }
    let mut p = (par - min) % period;
    if p < 0.0 {
        p += period;
    }
    min + p
}

/// OCCT-aligned: Bnd_Box2d::IsOut for a point (box gap 0, boundary not out).
/// Is the 2D point (u, v) outside the rectangle [u_min, v_min, u_max, v_max]?
fn is_out_rect(u: f64, v: f64, rect: &[f64; 4]) -> bool {
    u < rect[0] || u > rect[2] || v < rect[1] || v > rect[3]
}

/// OCCT-aligned: IntSurf_PntOn2S::IsSame (IntSurf_PntOn2S.cxx L85-113).
/// 3D distance
/// within theTol3d; when theTol2d >= 0 the UV params of both surfaces must
/// also coincide within theTol2d.
fn pnt_is_same(a: &WLinePnt, b: &WLinePnt, tol3d: f64, tol2d: f64) -> bool {
    if a.p3d.distance_squared(b.p3d) > tol3d * tol3d {
        return false;
    }
    if tol2d < 0.0 {
        return true; // Compare 3D-points only.
    }
    if (a.u1 - b.u1).abs() > tol2d || (a.v1 - b.v1).abs() > tol2d {
        return false;
    }
    if (a.u2 - b.u2).abs() > tol2d || (a.v2 - b.v2).abs() > tol2d {
        return false;
    }
    true
}

/// OCCT-aligned: IntPatch_WLineTool::CheckArgumentsToExtend (L918-1067).
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn check_arguments_to_extend(
    s1: &Surface3,
    s2: &Surface3,
    pt_wl1: &WLinePnt,
    pt_wl2: &WLinePnt,
    new_point: &mut WLinePnt,
    vec1: DVec3,
    vec2: DVec3,
    vec3: DVec3,
    box_s1: &[f64; 4],
    box_s2: &[f64; 4],
    tol_3d: f64,
    arr_periods: &[f64; 4],
) -> WlsConnectionType {
    let a_sq_toler = tol_3d * tol_3d;
    let mut a_ret_val = WlsConnectionType::NotConnected;
    if vec3.length_squared() <= a_sq_toler {
        if vec1.angle_between(vec2) > MY_MAX_CONCAT_ANGLE {
            return a_ret_val;
        }
        a_ret_val = WlsConnectionType::Common;
    } else if vec1.angle_between(vec2) > MY_MAX_CONCAT_ANGLE
        || vec1.angle_between(vec3) > MY_MAX_CONCAT_ANGLE
        || vec2.angle_between(vec3) > MY_MAX_CONCAT_ANGLE
    {
        return a_ret_val;
    }

    let a_pmid = 0.5 * (pt_wl1.p3d + pt_wl2.p3d);

    let mut a_new_par = [0.0f64; 4];
    let mut a_par_lbc = [0.0f64; 4];
    a_par_lbc[0] = box_s1[0];
    a_par_lbc[1] = box_s1[1];
    a_par_lbc[2] = box_s2[0];
    a_par_lbc[3] = box_s2[1];

    if !is_intersection_point(a_pmid, s1, s2, pt_wl1, tol_3d, arr_periods, new_point) {
        return WlsConnectionType::NotConnected;
    }

    if is_out_of_domain(box_s1, box_s2, new_point, arr_periods) {
        return WlsConnectionType::NotConnected;
    }

    let mut a_par_wl1 = [pt_wl1.u1, pt_wl1.v1, pt_wl1.u2, pt_wl1.v2];
    let mut a_par_wl2 = [pt_wl2.u1, pt_wl2.v1, pt_wl2.u2, pt_wl2.v2];
    a_new_par = [new_point.u1, new_point.v1, new_point.u2, new_point.v2];

    let mut is_on_boundary = false;
    for i in 0..4 {
        if arr_periods[i] == 0.0 {
            continue;
        }
        let r1 = a_par_wl1[i].min(a_par_wl2[i]);
        let r2 = a_par_wl1[i].max(a_par_wl2[i]);
        if range_is_intersected(r1, r2, a_par_lbc[i], arr_periods[i]) {
            min_max(&mut a_par_wl1[i], &mut a_par_wl2[i]);
            if a_new_par[i] > a_par_wl2[i] {
                let a_par =
                    a_par_wl1[i] + arr_periods[i] * ((a_new_par[i] - a_par_wl1[i]) / arr_periods[i]).ceil();
                a_par_wl1[i] = a_par_wl2[i];
                a_par_wl2[i] = a_par;
            } else if a_new_par[i] < a_par_wl1[i] {
                let a_par =
                    a_par_wl2[i] - arr_periods[i] * ((a_par_wl2[i] - a_new_par[i]) / arr_periods[i]).ceil();
                a_par_wl2[i] = a_par_wl1[i];
                a_par_wl1[i] = a_par;
            }
            let (r1a, r2a) = (a_par_wl1[i].min(a_new_par[i]), a_par_wl1[i].max(a_new_par[i]));
            let (r1b, r2b) = (a_new_par[i].min(a_par_wl2[i]), a_new_par[i].max(a_par_wl2[i]));
            if range_is_intersected(r1a, r2a, a_par_lbc[i], arr_periods[i])
                || range_is_intersected(r1b, r2b, a_par_lbc[i], arr_periods[i])
            {
                return WlsConnectionType::NotConnected;
            }
            is_on_boundary = true;
        }
    }

    if is_on_boundary {
        return WlsConnectionType::Singular;
    }
    if a_ret_val == WlsConnectionType::Common {
        return WlsConnectionType::Common;
    }
    WlsConnectionType::ReqExtend
}

/// OCCT-aligned: IntPatch_WLineTool::IsIntersectionPoint (L734-804).
fn is_intersection_point(
    pmid: DVec3,
    s1: &Surface3,
    s2: &Surface3,
    ref_pt: &WLinePnt,
    tol: f64,
    arr_periods: &[f64; 4],
    new_pt: &mut WLinePnt,
) -> bool {
    let (a_u1, a_v1) = surface_parameters(s1, pmid);
    let (a_u2, a_v2) = surface_parameters(s2, pmid);
    if a_u1.is_nan() || a_u2.is_nan() {
        return false;
    }
    new_pt.p3d = pmid;
    new_pt.u1 = a_u1;
    new_pt.v1 = a_v1;
    new_pt.u2 = a_u2;
    new_pt.v2 = a_v2;

    adjust_uv(ref_pt, arr_periods, new_pt);

    // OCCT L800-801: evaluate the surfaces at the ORIGINAL aU1/aV1/aU2/aV2
    // locals (before AdjustPointAndVertex), not at the adjusted new_pt params.
    let a_p1 = s1.point_at(a_u1, a_v1);
    let a_p2 = s2.point_at(a_u2, a_v2);
    a_p1.distance_squared(a_p2) <= tol * tol
}

/// OCCT-aligned: ElSLib::Parameters per surface type (IsIntersectionPoint
/// L744-794 switch).  Data-model adapter: analytic UV inversion.
fn surface_parameters(surf: &Surface3, p: DVec3) -> (f64, f64) {
    match surf {
        Surface3::Plane(pl) => {
            let d = p - pl.origin;
            (d.dot(pl.u_dir), d.dot(pl.v_dir))
        }
        Surface3::Cylinder(c) => {
            let uv = c.world_to_uv(p);
            (uv.x, uv.y)
        }
        Surface3::Sphere(s) => {
            let uv = s.world_to_uv(p);
            (uv.x, uv.y)
        }
        Surface3::Cone(c) => {
            let uv = c.world_to_uv(p);
            (uv.x, uv.y)
        }
        Surface3::Torus(t) => {
            let uv = t.world_to_uv(p);
            (uv.x, uv.y)
        }
        _ => (f64::NAN, f64::NAN),
    }
}

/// OCCT-aligned: IntPatch_SpecialPoints::AdjustPointAndVertex (L1082-1128):
/// shifts periodic params of new_point to be within half a period of ref_pt.
fn adjust_uv(ref_pt: &WLinePnt, arr_periods: &[f64; 4], new_point: &mut WLinePnt) {
    let mut a_par = [new_point.u1, new_point.v1, new_point.u2, new_point.v2];
    let mut a_ref_par = [0.0f64; 2];
    for i in 0..4 {
        if arr_periods[i] == 0.0 {
            continue;
        }
        let a_period = arr_periods[i];
        let a_half_period = 0.5 * a_period;
        if i < 2 {
            a_ref_par[0] = ref_pt.u1;
            a_ref_par[1] = ref_pt.v1;
        } else {
            a_ref_par[0] = ref_pt.u2;
            a_ref_par[1] = ref_pt.v2;
        }
        let a_ref_ind = i % 2;
        let mut a_delta_par = a_ref_par[a_ref_ind] - a_par[i];
        let an_incr = a_period.copysign(a_delta_par);
        while (a_delta_par > a_half_period) || (a_delta_par < -a_half_period) {
            a_par[i] += an_incr;
            a_delta_par = a_ref_par[a_ref_ind] - a_par[i];
        }
    }
    new_point.u1 = a_par[0];
    new_point.v1 = a_par[1];
    new_point.u2 = a_par[2];
    new_point.v2 = a_par[3];
}

/// OCCT-aligned: IntPatch_WLineTool::ExtendFirst (L810-850): adds thePnt to
/// the beginning of the line.
fn extend_first(slin: &mut Vec<IntPatchLine>, idx: usize, added_pt: WLinePnt) {
    let n = slin[idx].wline_pnts.len();
    if n == 0 {
        return;
    }
    if added_pt.p3d.distance(slin[idx].wline_pnts[0].p3d) <= CONFUSION {
        slin[idx].wline_pnts[0] = added_pt;
        for v in slin[idx].vertices.iter_mut() {
            // OCCT L822: if (aVert.ParameterOnLine() != 1) break;
            if v.param_on_line != 1.0 {
                break;
            }
            v.u1 = added_pt.u1;
            v.v1 = added_pt.v1;
            v.u2 = added_pt.u2;
            v.v2 = added_pt.v2;
            v.p3d = added_pt.p3d;
        }
        return;
    }
    slin[idx].wline_pnts.insert(0, added_pt);
    for v in slin[idx].vertices.iter_mut() {
        // OCCT L840: if (aVert.ParameterOnLine() == 1) ... else +1.
        if v.param_on_line == 1.0 {
            v.u1 = added_pt.u1;
            v.v1 = added_pt.v1;
            v.u2 = added_pt.u2;
            v.v2 = added_pt.v2;
            v.p3d = added_pt.p3d;
        } else {
            v.param_on_line += 1.0;
        }
    }
}

/// OCCT-aligned: IntPatch_WLineTool::ExtendLast (L856-884): adds thePnt to
/// the end of the line.
fn extend_last(slin: &mut Vec<IntPatchLine>, idx: usize, added_pt: WLinePnt) {
    let n = slin[idx].wline_pnts.len();
    if n == 0 {
        return;
    }
    if added_pt.p3d.distance(slin[idx].wline_pnts[n - 1].p3d) <= CONFUSION {
        slin[idx].wline_pnts[n - 1] = added_pt;
    } else {
        slin[idx].wline_pnts.push(added_pt);
    }
    let new_n = slin[idx].wline_pnts.len();
    for v in slin[idx].vertices.iter_mut().rev() {
        // OCCT L875: if (aVert.ParameterOnLine() != aNbPnts) break;
        if v.param_on_line != n as f64 {
            break;
        }
        v.u1 = added_pt.u1;
        v.v1 = added_pt.v1;
        v.u2 = added_pt.u2;
        v.v2 = added_pt.v2;
        v.p3d = added_pt.p3d;
        v.param_on_line = new_n as f64;
    }
}

fn period(s: &Surface3, is_u: bool) -> f64 {
    if is_u {
        if s.is_u_periodic() { std::f64::consts::TAU } else { 0.0 }
    } else if s.is_v_periodic() {
        std::f64::consts::TAU
    } else {
        0.0
    }
}

fn first_u(s: &Surface3) -> f64 {
    s.default_domain()[0]
}
fn first_v(s: &Surface3) -> f64 {
    s.default_domain()[2]
}
fn last_u(s: &Surface3) -> f64 {
    s.default_domain()[1]
}
fn last_v(s: &Surface3) -> f64 {
    s.default_domain()[3]
}
