//! WorkWithBoundaries and the cylinder-cylinder boundary helpers — 1:1
//! translation of OCCT `IntPatch_ImpImpIntersection.cxx`:
//!   - MinMax / ExtremaLineLine (L4400-4424)
//!   - VBoundaryPrecise (L4443-4530)
//!   - DeltaU1Computing (L4498-4516)
//!   - StepComputing (L4531-4682)
//!   - WorkWithBoundaries class (L4256-4387)
//!   - AddBoundaryPoint (L5818-6033)
//!   - SeekAdditionalPoints (L6034-6165)
//!   - BoundariesComputing (L6166-6382)
//!   - CriticalPointsComputing (L6383-6512)
//!   - BoundaryEstimation (L6513-6572)
//!
//! rcad data-model notes:
//!   - `Bnd_Box2d` -> `[f64; 4]` = `[u_min, v_min, u_max, v_max]`.
//!   - `IntSurf_LineOn2S` -> `cycy_common::WLine` (1-based semantics).
//!   - `math_Matrix`/`math_Vector` -> `cycy_common::Mat3`/`Mat35` and `DVec3`.
//!   - `MathRoot::Brent` -> `cycy_common::brent_root`.

use glam::DVec3;
use rcad_kernel::geom::CylindricalSurface;
use rcad_kernel::precision::{ANGULAR, CONFUSION, PCONFUSION};

use super::cycy_coeffs::{StCoeffsValue, cyl_cyl_compute_parameters};
use super::cycy_common::{
    A_NUL_VALUE, BndRange, BrentConfig, Mat3, Mat35, SolverStatus, WLine, brent_root,
    exclude_near_elements, inscribe_point, precision_infinite, precision_is_infinite, real_last,
};
use super::WLinePnt;
use crate::topalgo::int_surf::quadric::Quadric;

/// OCCT MinMax (L4400-4410): replaces theParMIN = MIN, theParMAX = MAX.
pub fn min_max(a: &mut f64, b: &mut f64) {
    if *a > *b {
        let aux = *b;
        *b = *a;
        *a = aux;
    }
}

/// OCCT ExtremaLineLine (L4412-4424) — computes extrema between two lines.
/// d1/d2 are the line directions, l1l2 = Loc2 - Loc1, returns parameters on
/// each line.
pub fn extrema_line_line(
    d1: DVec3,
    d2: DVec3,
    l1l2: DVec3,
    cos_a: f64,
    sq_sin_a: f64,
    par1: &mut f64,
    par2: &mut f64,
) {
    let a_d1l = d1.dot(l1l2);
    let a_d2l = d2.dot(l1l2);
    *par1 = (a_d1l - cos_a * a_d2l) / sq_sin_a;
    *par2 = (cos_a * a_d1l - a_d2l) / sq_sin_a;
}

/// OCCT VBoundaryPrecise (L4443-4530).
/// By default V1 and V2 are considered to increase when U1 increases; if that
/// is not the case, new V1set and/or V2set must be computed as
/// [V_current - DeltaV].  This function processes this case.
pub fn v_boundary_precise(
    matr: &Mat35,
    v1_after_decr_by_delta: f64,
    v2_after_decr_by_delta: f64,
    v1_set: &mut f64,
    v2_set: &mut f64,
) {
    // Now we are going to define if V1 (and V2) increases
    // (or decreases) when U1 will increase.
    let mut a_syst = Mat3::new();
    a_syst.set_col(1, matr.col(1));
    a_syst.set_col(2, matr.col(2));
    a_syst.set_col(3, matr.col(4));
    // We have the system (see comment to StepComputing):
    //     {a11*dV1 + a12*dV2 + a14*dU2 = -a13*dU1
    //     {a21*dV1 + a22*dV2 + a24*dU2 = -a23*dU1
    //     {a31*dV1 + a32*dV2 + a34*dU2 = -a33*dU1
    let a_det = a_syst.determinant();
    a_syst.set_col(1, matr.col(3));
    let a_det1 = a_syst.determinant();
    a_syst.set_col(1, matr.col(1));
    a_syst.set_col(2, matr.col(3));
    let a_det2 = a_syst.determinant();
    // Now, dV1 = -dU1*aDet1/aDet, dV2 = -dU1*aDet2/aDet.
    if a_det * a_det1 > 0.0 {
        *v1_set = v1_after_decr_by_delta;
    }
    if a_det * a_det2 > 0.0 {
        *v2_set = v2_after_decr_by_delta;
    }
}

/// OCCT DeltaU1Computing (L4498-4516) — computes new step for U1 parameter.
pub fn delta_u1_computing(syst: &Mat3, free: &[f64; 3], delta_u1_found: &mut f64) -> bool {
    let a_det = syst.determinant();
    if a_det.abs() > A_NUL_VALUE {
        let mut a_syst1 = syst.clone();
        a_syst1.set_col(2, *free);
        *delta_u1_found = (a_syst1.determinant() / a_det).abs();
        return true;
    }
    false
}

/// OCCT StepComputing (L4531-4682).
/// theMatr must have 3*5-dimension strictly.  For the system
///     {a11*V1+a12*V2+a13*dU1+a14*dU2=b1;
///     {a21*V1+a22*V2+a23*dU1+a24*dU2=b2;
///     {a31*V1+a32*V2+a33*dU1+a34*dU2=b3;
/// theMatr is (a11 a12 a13 a14 b1) / (a21 ... b2) / (a31 ... b3).
pub fn step_computing(
    matr: &Mat35,
    v1_cur: f64,
    v2_cur: f64,
    delta_v1: f64,
    delta_v2: f64,
    delta_u1_found: &mut f64,
) -> bool {
    let mut is_success = false;
    *delta_u1_found = real_last();
    let mut a_syst = Mat3::new();
    let mut a_free = [0.0; 3];

    // By default, increasing V1(U1) and V2(U1) functions is considered.
    let mut a_v1_set = v1_cur + delta_v1;
    let mut a_v2_set = v2_cur + delta_v2;
    v_boundary_precise(
        matr,
        v1_cur - delta_v1,
        v2_cur - delta_v2,
        &mut a_v1_set,
        &mut a_v2_set,
    );

    a_syst.set_col(2, matr.col(3));
    a_syst.set_col(3, matr.col(4));

    for i in 0..2 {
        if i == 0 {
            // V1 is known
            a_syst.set_col(1, matr.col(2));
            let col5 = matr.col(5);
            let col1 = matr.col(1);
            for r in 0..3 {
                a_free[r] = col5[r] - a_v1_set * col1[r];
            }
        } else {
            // i==1 => V2 is known
            a_syst.set_col(1, matr.col(1));
            let col5 = matr.col(5);
            let col2 = matr.col(2);
            for r in 0..3 {
                a_free[r] = col5[r] - a_v2_set * col2[r];
            }
        }

        let mut a_new_du = *delta_u1_found;
        if delta_u1_computing(&a_syst, &a_free, &mut a_new_du) {
            is_success = true;
            if a_new_du < *delta_u1_found {
                *delta_u1_found = a_new_du;
            }
        }
    }

    if !is_success {
        let col5 = matr.col(5);
        let col1 = matr.col(1);
        let col2 = matr.col(2);
        for r in 0..3 {
            a_free[r] = col5[r] - a_v1_set * col1[r] - a_v2_set * col2[r];
        }
        // (OCCT builds a dead 3x2 aSyst1 here that is never used.)

        // Now we have overdetermined system.
        let a_det1 = matr.get(1, 3) * matr.get(2, 4) - matr.get(2, 3) * matr.get(1, 4);
        let a_det2 = matr.get(1, 3) * matr.get(3, 4) - matr.get(3, 3) * matr.get(1, 4);
        let a_det3 = matr.get(2, 3) * matr.get(3, 4) - matr.get(3, 3) * matr.get(2, 4);
        let an_abs_d1 = a_det1.abs();
        let an_abs_d2 = a_det2.abs();
        let an_abs_d3 = a_det3.abs();

        if an_abs_d1 >= an_abs_d2 {
            if an_abs_d1 >= an_abs_d3 {
                // Det1
                if an_abs_d1 <= A_NUL_VALUE {
                    return is_success;
                }
                *delta_u1_found =
                    (a_free[0] * matr.get(2, 4) - a_free[1] * matr.get(1, 4)).abs() / an_abs_d1;
                is_success = true;
            } else {
                // Det3
                if an_abs_d3 <= A_NUL_VALUE {
                    return is_success;
                }
                *delta_u1_found =
                    (a_free[1] * matr.get(3, 4) - a_free[2] * matr.get(2, 4)).abs() / an_abs_d3;
                is_success = true;
            }
        } else {
            if an_abs_d2 >= an_abs_d3 {
                // Det2
                if an_abs_d2 <= A_NUL_VALUE {
                    return is_success;
                }
                *delta_u1_found =
                    (a_free[0] * matr.get(3, 4) - a_free[2] * matr.get(1, 4)).abs() / an_abs_d2;
                is_success = true;
            } else {
                // Det3
                if an_abs_d3 <= A_NUL_VALUE {
                    return is_success;
                }
                *delta_u1_found =
                    (a_free[1] * matr.get(3, 4) - a_free[2] * matr.get(2, 4)).abs() / an_abs_d3;
                is_success = true;
            }
        }
    }

    is_success
}

/// OCCT WorkWithBoundaries (L4256-4387).
pub struct WorkWithBoundaries<'a> {
    my_quad1: &'a Quadric,
    my_quad2: &'a Quadric,
    my_coeffs: &'a StCoeffsValue,
    my_uv_surf1: [f64; 4],
    my_uv_surf2: [f64; 4],
    #[allow(dead_code)]
    my_nb_wlines: usize,
    my_period: f64,
    my_tol_3d: f64,
    my_tol_2d: f64,
    my_is_reverse: bool,
}

/// OCCT WorkWithBoundaries::SearchBoundType (L4260-4264).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchBoundType {
    SearchNone = 0,
    SearchV1 = 1,
    SearchV2 = 2,
}

/// OCCT WorkWithBoundaries::StPInfo (L4267-4289) — a boundary point candidate.
#[derive(Debug, Clone)]
pub struct StPInfo {
    /// Equal to 0 for 1st surface, non-zero for 2nd one.
    pub my_surf_id: usize,
    pub my_u1: f64,
    pub my_v1: f64,
    pub my_u2: f64,
    pub my_v2: f64,
}

impl StPInfo {
    pub fn new() -> Self {
        StPInfo {
            my_surf_id: 0,
            my_u1: real_last(),
            my_v1: real_last(),
            my_u2: real_last(),
            my_v2: real_last(),
        }
    }
}

impl Default for StPInfo {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> WorkWithBoundaries<'a> {
    /// OCCT WorkWithBoundaries ctor (L4294-4322).
    pub fn new(
        quad1: &'a Quadric,
        quad2: &'a Quadric,
        coeffs: &'a StCoeffsValue,
        uv_surf1: [f64; 4],
        uv_surf2: [f64; 4],
        nb_wlines: usize,
        period: f64,
        tol_3d: f64,
        tol_2d: f64,
        is_reverse: bool,
    ) -> WorkWithBoundaries<'a> {
        WorkWithBoundaries {
            my_quad1: quad1,
            my_quad2: quad2,
            my_coeffs: coeffs,
            my_uv_surf1: uv_surf1,
            my_uv_surf2: uv_surf2,
            my_nb_wlines: nb_wlines,
            my_period: period,
            my_tol_3d: tol_3d,
            my_tol_2d: tol_2d,
            my_is_reverse: is_reverse,
        }
    }

    /// OCCT SICoeffs() (L4324).
    pub fn si_coeffs(&self) -> &StCoeffsValue {
        self.my_coeffs
    }

    /// OCCT GetQSurface(theIdx) (L4332-4340).
    pub fn get_q_surface(&self, idx: usize) -> &Quadric {
        if idx <= 1 {
            self.my_quad1
        } else {
            self.my_quad2
        }
    }

    /// OCCT IsReversed() (L4342).
    pub fn is_reversed(&self) -> bool {
        self.my_is_reverse
    }

    /// OCCT Get2dTolerance() (L4346).
    pub fn get_2d_tolerance(&self) -> f64 {
        self.my_tol_2d
    }

    /// OCCT Get3dTolerance() (L4350).
    pub fn get_3d_tolerance(&self) -> f64 {
        self.my_tol_3d
    }

    /// OCCT UVS1() (L4354).
    pub fn uv_s1(&self) -> [f64; 4] {
        self.my_uv_surf1
    }

    /// OCCT UVS2() (L4358).
    pub fn uv_s2(&self) -> [f64; 4] {
        self.my_uv_surf2
    }

    /// OCCT WorkWithBoundaries::AddBoundaryPoint (L5818-6033).
    ///
    /// theU1Prev together with theU1 form the bracket across which the straddle
    /// on V_bound was detected, and in which the V-boundary crossing is resolved
    /// by a 1D branch-aware Brent root-finder on the analytical curve.
    #[allow(clippy::too_many_arguments)]
    pub fn add_boundary_point(
        &self,
        the_wl: &mut WLine,
        the_u1: f64,
        the_u1_prev: f64,
        the_u1_min: f64,
        the_u2: f64,
        the_v1: f64,
        the_v1_prev: f64,
        the_v2: f64,
        the_v2_prev: f64,
        the_wl_index: usize,
        the_fl_force: bool,
        is_the_found1: &mut bool,
        is_the_found2: &mut bool,
    ) {
        let (a_usurf1f, a_vsurf1f, a_usurf1l, a_vsurf1l) =
            (self.my_uv_surf1[0], self.my_uv_surf1[1], self.my_uv_surf1[2], self.my_uv_surf1[3]);
        let (a_usurf2f, a_vsurf2f, a_usurf2l, a_vsurf2l) =
            (self.my_uv_surf2[0], self.my_uv_surf2[1], self.my_uv_surf2[2], self.my_uv_surf2[3]);

        let a_size = 4;
        let an_arr_vzad = [a_vsurf1f, a_vsurf1l, a_vsurf2f, a_vsurf2l];

        let mut a_uv_point: [StPInfo; 4] = [
            StPInfo::new(),
            StPInfo::new(),
            StPInfo::new(),
            StPInfo::new(),
        ];

        // Branch-aware 1D root-finder for V(U1) = V_bound on the analytical
        // intersection curve, in place of the original rank-deficient 3x3 Newton
        // (SearchOnVBounds).  The intersection curve is 1D in U1 on each arccos
        // branch, so the correct numerical tool is a bracketed 1D root finder.
        // Brent converges robustly even when the curve grazes V_bound at a
        // near-tangent (two close crossings on opposite sides of a V-extremum) —
        // the configuration that makes the 3x3 Newton Jacobian rank-deficient.
        // The only tolerance is Precision::PConfusion() (OCCT standard
        // parametric equality) for termination on U.
        let coeffs = self.my_coeffs;
        let find_v_bound_crossing = |is_v1: bool,
                                     v_bound: f64,
                                     u_lo: f64,
                                     u_hi: f64,
                                     u1_star: &mut f64|
         -> bool {
            if !(u_lo < u_hi) {
                return false;
            }
            let mut func = |x: f64| -> Option<f64> {
                let mut a_u2 = 0.0;
                let mut a_v1 = 0.0;
                let mut a_v2 = 0.0;
                if !cyl_cyl_compute_parameters(x, the_wl_index as i32, coeffs, &mut a_u2, &mut a_v1, &mut a_v2)
                {
                    return None;
                }
                Some((if is_v1 { a_v1 } else { a_v2 }) - v_bound)
            };
            let a_flo = match func(u_lo) {
                Some(v) => v,
                None => return false,
            };
            let a_fhi = match func(u_hi) {
                Some(v) => v,
                None => return false,
            };
            if a_flo == 0.0 {
                *u1_star = u_lo;
                return true;
            }
            if a_fhi == 0.0 {
                *u1_star = u_hi;
                return true;
            }
            if a_flo * a_fhi > 0.0 {
                return false; // no bracketed crossing
            }
            let mut a_cfg = BrentConfig::new();
            a_cfg.x_tolerance = PCONFUSION;
            let a_res = brent_root(&mut func, u_lo, u_hi, &a_cfg);
            if a_res.status != SolverStatus::Ok || a_res.root.is_none() {
                return false;
            }
            *u1_star = a_res.root.unwrap();
            true
        };

        let mut an_id_surf = 0usize;
        while an_id_surf < 4 {
            let a_vf = if an_id_surf == 0 { the_v1 } else { the_v2 };
            let a_vl = if an_id_surf == 0 { the_v1_prev } else { the_v2_prev };
            let a_ts = if an_id_surf == 0 { SearchBoundType::SearchV1 } else { SearchBoundType::SearchV2 };

            for an_id_bound in 0..2 {
                let an_index = an_id_surf + an_id_bound;
                a_uv_point[an_index].my_surf_id = an_id_surf;

                if (a_vf - an_arr_vzad[an_index]).abs() > self.my_tol_2d
                    && (a_vf - an_arr_vzad[an_index]) * (a_vl - an_arr_vzad[an_index]) > 0.0
                {
                    continue;
                }

                // Segment [aVf, aVl] intersects at least one V-boundary.
                let a_u_lo = the_u1_prev.min(the_u1);
                let a_u_hi = the_u1_prev.max(the_u1);
                let a_res = find_v_bound_crossing(
                    a_ts == SearchBoundType::SearchV1,
                    an_arr_vzad[an_index],
                    a_u_lo,
                    a_u_hi,
                    &mut a_uv_point[an_index].my_u1,
                );

                // aUVPoint[anIndex].myU1 is considered to be nearer to theU1 than
                // to theU1+/-Period.
                if !a_res
                    || a_uv_point[an_index].my_u1 >= the_u1
                    || a_uv_point[an_index].my_u1 < the_u1_min
                {
                    // Intersection point is not found or out of the domain.
                    a_uv_point[an_index].my_u1 = real_last();
                    continue;
                }
                // Intersection point is found.
                let a_u1 = a_uv_point[an_index].my_u1;
                let mut a_u2 = the_u2;
                let mut a_v1 = the_v1;
                let mut a_v2 = the_v2;
                if !cyl_cyl_compute_parameters(
                    a_u1,
                    the_wl_index as i32,
                    coeffs,
                    &mut a_u2,
                    &mut a_v1,
                    &mut a_v2,
                ) {
                    // Found point is wrong.
                    a_uv_point[an_index].my_u1 = real_last();
                    continue;
                }
                // Point on true V-boundary.
                if a_ts == SearchBoundType::SearchV1 {
                    a_v1 = an_arr_vzad[an_index];
                } else {
                    a_v2 = an_arr_vzad[an_index];
                }
                a_uv_point[an_index].my_u1 = a_u1;
                a_uv_point[an_index].my_u2 = a_u2;
                a_uv_point[an_index].my_v1 = a_v1;
                a_uv_point[an_index].my_v2 = a_v2;
            }
            an_id_surf += 2;
        }

        // Sort with ascending U1-parameter.
        a_uv_point.sort_by(|a, b| a.my_u1.partial_cmp(&b.my_u1).unwrap());

        *is_the_found1 = false;
        *is_the_found2 = false;

        // Adding found points on boundary in the WLine.
        for i in 0..a_size {
            if a_uv_point[i].my_u1 == real_last() {
                break;
            }
            if !super::cycy_walking::add_point_into_wl(
                self.my_quad1,
                self.my_quad2,
                coeffs,
                self.my_is_reverse,
                false,
                [a_uv_point[i].my_u1, a_uv_point[i].my_v1],
                [a_uv_point[i].my_u2, a_uv_point[i].my_v2],
                a_usurf1f,
                a_usurf1l,
                a_usurf2f,
                a_usurf2l,
                a_vsurf1f,
                a_vsurf1l,
                a_vsurf2f,
                a_vsurf2l,
                self.my_period,
                the_wl,
                the_wl_index,
                self.my_tol_3d,
                self.my_tol_2d,
                the_fl_force,
                false,
            ) {
                continue;
            }
            if a_uv_point[i].my_surf_id == 0 {
                *is_the_found1 = true;
            } else {
                *is_the_found2 = true;
            }
        }
    }

    /// OCCT WorkWithBoundaries::BoundaryEstimation (L6513-6572) — rough
    /// estimation of the parameter range.
    pub fn boundary_estimation(
        &self,
        cy1: &CylindricalSurface,
        cy2: &CylindricalSurface,
        out_box_s1: &mut BndRange,
        out_box_s2: &mut BndRange,
    ) {
        let a_d1 = cy1.axis.normalize_or_zero();
        let a_d2 = cy2.axis.normalize_or_zero();
        let a_r1 = cy1.radius;
        let a_r2 = cy2.radius;

        // Consider a parallelogram whose edges are parallel to aD1 and aD2 and
        // whose altitudes equal 2*aR1 and 2*aR2 (diameters of the cylinders).
        let a_cos_a = a_d1.dot(a_d2);
        let a_sq_sin_a = a_d1.cross(a_d2).length_squared();

        // If sine is small then it can be compared with angle.
        if a_sq_sin_a < ANGULAR * ANGULAR {
            return;
        }

        // Half of delta V — the distance between projections of two opposite
        // parallelogram vertices (joined by the maximal diagonal) to the cylinder axis.
        let a_sin_a = a_sq_sin_a.sqrt();
        let an_abs_cos_a = a_cos_a.abs();
        let a_hdv1 = (a_r1 * an_abs_cos_a + a_r2) / a_sin_a;
        let a_hdv2 = (a_r2 * an_abs_cos_a + a_r1) / a_sin_a;

        // V-parameters of the intersection point of the axes.
        let mut a_v01 = 0.0;
        let mut a_v02 = 0.0;
        extrema_line_line(a_d1, a_d2, cy2.origin - cy1.origin, a_cos_a, a_sq_sin_a, &mut a_v01, &mut a_v02);

        out_box_s1.add(a_v01 - a_hdv1);
        out_box_s1.add(a_v01 + a_hdv1);
        out_box_s2.add(a_v02 - a_hdv2);
        out_box_s2.add(a_v02 + a_hdv2);

        out_box_s1.enlarge(CONFUSION);
        out_box_s2.enlarge(CONFUSION);

        let (_a_u1, a_v1, _a_u2, a_v2) =
            (self.my_uv_surf1[0], self.my_uv_surf1[1], self.my_uv_surf1[2], self.my_uv_surf1[3]);
        out_box_s1.common(&BndRange::with_bounds(a_v1, a_v2));

        let (_a_u1, a_v1, _a_u2, a_v2) =
            (self.my_uv_surf2[0], self.my_uv_surf2[1], self.my_uv_surf2[2], self.my_uv_surf2[3]);
        out_box_s2.common(&BndRange::with_bounds(a_v1, a_v2));
    }
}

/// OCCT WorkWithBoundaries::BoundariesComputing (L6166-6382) — computes the true
/// domain of the future intersection curve.  Returns false if there is no
/// solution for U1.
pub fn boundaries_computing(coeffs: &StCoeffsValue, period: f64, u_range: &mut [BndRange; 2]) -> bool {
    // We have the equation cos(U2-FI2) = B*cos(U1-FI1) + C, hence
    //     -1 <= B*cos(U1-FI1)+C <= 1.

    if coeffs.m_b > 0.0 {
        // -(1+C)/B <= cos(U1-FI1) <= (1-C)/B
        if coeffs.m_b + coeffs.m_c.abs() < -1.0 {
            // (1-C)/B < -1 or -(1+C)/B > 1 ==> No solution
            return false;
        } else if coeffs.m_b + coeffs.m_c.abs() <= 1.0 {
            // (1-C)/B >= 1 and -(1+C)/B <= -1 ==> U=[0;2*PI]+aFI1
            u_range[0].add(coeffs.m_fi1);
            u_range[0].add(period + coeffs.m_fi1);
        } else if (1.0 + coeffs.m_c <= coeffs.m_b) && (coeffs.m_b <= 1.0 - coeffs.m_c) {
            // U=[0;aDAngle]+aFI1 || U=[2*PI-aDAngle;2*PI]+aFI1,
            // where aDAngle = acos(-(myCoeffs.mC + 1) / myCoeffs.mB)
            let mut an_arg = -(coeffs.m_c + 1.0) / coeffs.m_b;
            if an_arg > 1.0 {
                an_arg = 1.0;
            }
            if an_arg < -1.0 {
                an_arg = -1.0;
            }
            let a_d_angle = an_arg.acos();
            u_range[0].add(coeffs.m_fi1);
            u_range[0].add(a_d_angle + coeffs.m_fi1);
            u_range[1].add(period - a_d_angle + coeffs.m_fi1);
            u_range[1].add(period + coeffs.m_fi1);
        } else if (1.0 - coeffs.m_c <= coeffs.m_b) && (coeffs.m_b <= 1.0 + coeffs.m_c) {
            // U=[aDAngle;2*PI-aDAngle]+aFI1, where aDAngle = acos((1 - mC)/mB)
            let mut an_arg = (1.0 - coeffs.m_c) / coeffs.m_b;
            if an_arg > 1.0 {
                an_arg = 1.0;
            }
            if an_arg < -1.0 {
                an_arg = -1.0;
            }
            let a_d_angle = an_arg.acos();
            u_range[0].add(a_d_angle + coeffs.m_fi1);
            u_range[0].add(period - a_d_angle + coeffs.m_fi1);
        } else if coeffs.m_b - coeffs.m_c.abs() >= 1.0 {
            // U=[aDAngle1;aDAngle2]+aFI1 || U=[2*PI-aDAngle2;2*PI-aDAngle1]+aFI1
            let mut an_arg1 = (1.0 - coeffs.m_c) / coeffs.m_b;
            let mut an_arg2 = -(coeffs.m_c + 1.0) / coeffs.m_b;
            if an_arg1 > 1.0 {
                an_arg1 = 1.0;
            }
            if an_arg1 < -1.0 {
                an_arg1 = -1.0;
            }
            if an_arg2 > 1.0 {
                an_arg2 = 1.0;
            }
            if an_arg2 < -1.0 {
                an_arg2 = -1.0;
            }
            let a_d_angle1 = an_arg1.acos();
            let a_d_angle2 = an_arg2.acos();
            u_range[0].add(a_d_angle1 + coeffs.m_fi1);
            u_range[0].add(a_d_angle2 + coeffs.m_fi1);
            u_range[1].add(period - a_d_angle2 + coeffs.m_fi1);
            u_range[1].add(period - a_d_angle1 + coeffs.m_fi1);
        } else {
            return false;
        }
    } else if coeffs.m_b < 0.0 {
        // (1-C)/B <= cos(U1-FI1) <= -(1+C)/B
        if coeffs.m_b + coeffs.m_c.abs() > 1.0 {
            // -(1+C)/B < -1 or (1-C)/B > 1 ==> No solutions
            return false;
        } else if -coeffs.m_b + coeffs.m_c.abs() <= 1.0 {
            // -(1+C)/B >= 1 and (1-C)/B <= -1 ==> U=[0;2*PI]+aFI1
            u_range[0].add(coeffs.m_fi1);
            u_range[0].add(period + coeffs.m_fi1);
        } else if (-coeffs.m_c - 1.0 <= coeffs.m_b) && (coeffs.m_b <= coeffs.m_c - 1.0) {
            // U=[0;aDAngle]+aFI1 || U=[2*PI-aDAngle;2*PI]+aFI1,
            // where aDAngle = acos((1 - myCoeffs.mC) / myCoeffs.mB)
            let mut an_arg = (1.0 - coeffs.m_c) / coeffs.m_b;
            if an_arg > 1.0 {
                an_arg = 1.0;
            }
            if an_arg < -1.0 {
                an_arg = -1.0;
            }
            let a_d_angle = an_arg.acos();
            u_range[0].add(coeffs.m_fi1);
            u_range[0].add(a_d_angle + coeffs.m_fi1);
            u_range[1].add(period - a_d_angle + coeffs.m_fi1);
            u_range[1].add(period + coeffs.m_fi1);
        } else if (coeffs.m_c - 1.0 <= coeffs.m_b) && (coeffs.m_b <= -coeffs.m_b - 1.0) {
            // U=[aDAngle;2*PI-aDAngle]+aFI1,
            // where aDAngle = acos(-(myCoeffs.mC + 1) / myCoeffs.mB)
            let mut an_arg = -(coeffs.m_c + 1.0) / coeffs.m_b;
            if an_arg > 1.0 {
                an_arg = 1.0;
            }
            if an_arg < -1.0 {
                an_arg = -1.0;
            }
            let a_d_angle = an_arg.acos();
            u_range[0].add(a_d_angle + coeffs.m_fi1);
            u_range[0].add(period - a_d_angle + coeffs.m_fi1);
        } else if -coeffs.m_b - coeffs.m_c.abs() >= 1.0 {
            // U=[aDAngle1;aDAngle2]+aFI1 || U=[2*PI-aDAngle2;2*PI-aDAngle1]+aFI1,
            // where aDAngle1 = acos(-(mC + 1)/mB), aDAngle2 = acos((1 - mC)/mB)
            let mut an_arg1 = -(coeffs.m_c + 1.0) / coeffs.m_b;
            let mut an_arg2 = (1.0 - coeffs.m_c) / coeffs.m_b;
            if an_arg1 > 1.0 {
                an_arg1 = 1.0;
            }
            if an_arg1 < -1.0 {
                an_arg1 = -1.0;
            }
            if an_arg2 > 1.0 {
                an_arg2 = 1.0;
            }
            if an_arg2 < -1.0 {
                an_arg2 = -1.0;
            }
            let a_d_angle1 = an_arg1.acos();
            let a_d_angle2 = an_arg2.acos();
            u_range[0].add(a_d_angle1 + coeffs.m_fi1);
            u_range[0].add(a_d_angle2 + coeffs.m_fi1);
            u_range[1].add(period - a_d_angle2 + coeffs.m_fi1);
            u_range[1].add(period - a_d_angle1 + coeffs.m_fi1);
        } else {
            return false;
        }
    } else {
        return false;
    }

    true
}

/// OCCT CriticalPointsComputing (L6383-6512).
/// theNbCritPointsMax contains the true number of critical points; it must be
/// initialized correctly (to the array length) before calling.
pub fn critical_points_computing(
    coeffs: &StCoeffsValue,
    u_surf1f: f64,
    u_surf1l: f64,
    u_surf2f: f64,
    u_surf2l: f64,
    period: f64,
    tol_2d: f64,
    nb_crit_points_max: &mut usize,
    u1crit: &mut [f64],
) {
    // [0...1] — U1 goes through the seam-edge of the first cylinder.
    // [2...3] — First and last U1 parameter.
    // [4...5] — U2 goes through the seam-edge of the second cylinder.
    // [6...9] — the intersection line goes through U-boundaries of the 2nd surface.
    // [10...11] — boundary of monotonicity interval of U2(U1).

    u1crit[0] = 0.0;
    u1crit[1] = period;
    u1crit[2] = u_surf1f;
    u1crit[3] = u_surf1l;

    let a_cos = coeffs.m_fi2.cos();
    let a_bsb = coeffs.m_b.abs();
    if (coeffs.m_c - a_bsb <= a_cos) && (a_cos <= coeffs.m_c + a_bsb) {
        let mut an_arg = (a_cos - coeffs.m_c) / coeffs.m_b;
        if an_arg > 1.0 {
            an_arg = 1.0;
        }
        if an_arg < -1.0 {
            an_arg = -1.0;
        }
        u1crit[4] = -an_arg.acos() + coeffs.m_fi1;
        u1crit[5] = an_arg.acos() + coeffs.m_fi1;
    }

    let mut a_sf = (u_surf2f - coeffs.m_fi2).cos();
    let mut a_sl = (u_surf2l - coeffs.m_fi2).cos();
    min_max(&mut a_sf, &mut a_sl);

    // In accordance with pure mathematics, theU1crit[6] and [8] must be
    // -Precision::Infinite() instead of used +Precision::Infinite().
    u1crit[6] = if ((a_sl - coeffs.m_c) / coeffs.m_b).abs() < 1.0 {
        -((a_sl - coeffs.m_c) / coeffs.m_b).acos() + coeffs.m_fi1
    } else {
        precision_infinite()
    };
    u1crit[7] = if ((a_sf - coeffs.m_c) / coeffs.m_b).abs() < 1.0 {
        -((a_sf - coeffs.m_c) / coeffs.m_b).acos() + coeffs.m_fi1
    } else {
        precision_infinite()
    };
    u1crit[8] = if ((a_sf - coeffs.m_c) / coeffs.m_b).abs() < 1.0 {
        ((a_sf - coeffs.m_c) / coeffs.m_b).acos() + coeffs.m_fi1
    } else {
        precision_infinite()
    };
    u1crit[9] = if ((a_sl - coeffs.m_c) / coeffs.m_b).abs() < 1.0 {
        ((a_sl - coeffs.m_c) / coeffs.m_b).acos() + coeffs.m_fi1
    } else {
        precision_infinite()
    };

    u1crit[10] = coeffs.m_fi1;
    u1crit[11] = std::f64::consts::PI + coeffs.m_fi1;

    // Preparative treatment of array.  This array must have failed to contain
    // negative infinity number.
    for i in 0..*nb_crit_points_max {
        if precision_is_infinite(u1crit[i]) {
            continue;
        }
        u1crit[i] = u1crit[i] % period;
        if u1crit[i] < 0.0 {
            u1crit[i] += period;
        }
    }

    // Here all not infinite elements of theU1crit are in [0, thePeriod) range.
    loop {
        u1crit[..*nb_crit_points_max].sort_by(|a, b| a.partial_cmp(b).unwrap());
        if !exclude_near_elements(
            &mut u1crit[..*nb_crit_points_max],
            *nb_crit_points_max,
            u_surf1f,
            u_surf1l,
            tol_2d,
        ) {
            break;
        }
    }

    // Here all not infinite elements in theU1crit are different and sorted.
    while *nb_crit_points_max > 0 {
        let an_b = u1crit[*nb_crit_points_max - 1];
        if precision_is_infinite(an_b) {
            *nb_crit_points_max -= 1;
            continue;
        }
        // 1st not infinite element is found.
        if *nb_crit_points_max == 1 {
            break;
        }
        // Here theNbCritPointsMax > 1.
        let an_a = u1crit[0];
        // Compare 1st and last significant elements of theU1crit;
        // they may still differ by period.
        if (an_b - an_a - period).abs() < tol_2d {
            // E.g. anA == 2.0e-17, anB == (thePeriod - 1.0e-18).
            u1crit[0] = (an_a + an_b - period) / 2.0;
            u1crit[*nb_crit_points_max - 1] = precision_infinite();
            *nb_crit_points_max -= 1;
        }
        break;
    }
}

/// OCCT SeekAdditionalPoints (L6034-6165) — inserts additional intersection
/// points between neighbor points, splitting every interval in the middle until
/// the line contains at least theMinNbPoints points.
#[allow(clippy::too_many_arguments, unused_assignments, unused_variables)]
pub fn seek_additional_points(
    quad1: &Quadric,
    quad2: &Quadric,
    line: &mut WLine,
    coeffs: &StCoeffsValue,
    wl_index: usize,
    min_nb_points: usize,
    start_point_on_line: usize,
    end_point_on_line: usize,
    tol_2d: f64,
    period_of_surf2: f64,
    is_reverse: bool,
) {
    let mut a_nb_points = end_point_on_line - start_point_on_line + 1;

    let mut a_min_delta_param = tol_2d;

    {
        let mut u1 = 0.0;
        let mut v1 = 0.0;
        let mut u2 = 0.0;
        let mut v2 = 0.0;
        if is_reverse {
            let p1 = line.value(start_point_on_line);
            u1 = p1.u2;
            v1 = p1.v2;
            let p2 = line.value(end_point_on_line);
            u2 = p2.u2;
            v2 = p2.v2;
        } else {
            let p1 = line.value(start_point_on_line);
            u1 = p1.u1;
            v1 = p1.v1;
            let p2 = line.value(end_point_on_line);
            u2 = p2.u1;
            v2 = p2.v1;
        }
        a_min_delta_param = ((u2 - u1).abs() / min_nb_points as f64).max(a_min_delta_param);
    }

    let mut a_last_point_index = end_point_on_line;
    let mut u1_prec = 0.0;
    let mut v1_prec = 0.0;
    let mut u2_prec = 0.0;
    let mut v2_prec = 0.0;

    let mut a_nb_points_prev = 0;
    loop {
        a_nb_points_prev = a_nb_points;
        let mut fp = start_point_on_line;
        while fp < a_last_point_index {
            let mut u1f = 0.0;
            let mut v1f = 0.0;
            let mut u1l = 0.0;
            let mut v1l = 0.0;
            let mut u2f = 0.0;
            let mut v2f = 0.0;
            let mut u2l = 0.0;
            let mut v2l = 0.0;

            let lp = fp + 1;

            if is_reverse {
                let pf = line.value(fp);
                u1f = pf.u2;
                v1f = pf.v2;
                let pl = line.value(lp);
                u1l = pl.u2;
                v1l = pl.v2;
                let pf = line.value(fp);
                u2f = pf.u1;
                v2f = pf.v1;
                let pl = line.value(lp);
                u2l = pl.u1;
                v2l = pl.v1;
            } else {
                let pf = line.value(fp);
                u1f = pf.u1;
                v1f = pf.v1;
                let pl = line.value(lp);
                u1l = pl.u1;
                v1l = pl.v1;
                let pf = line.value(fp);
                u2f = pf.u2;
                v2f = pf.v2;
                let pl = line.value(lp);
                u2l = pl.u2;
                v2l = pl.v2;
            }

            if (u1l - u1f).abs() <= a_min_delta_param {
                // Step is minimal; it is not necessary to divide it.
                fp = lp + 1;
                continue;
            }

            u1_prec = 0.5 * (u1f + u1l);

            if !cyl_cyl_compute_parameters(
                u1_prec,
                wl_index as i32,
                coeffs,
                &mut u2_prec,
                &mut v1_prec,
                &mut v2_prec,
            ) {
                fp = lp + 1;
                continue;
            }

            min_max(&mut u2f, &mut u2l);
            if !inscribe_point(u2f, u2l, &mut u2_prec, tol_2d, period_of_surf2, false) {
                fp = lp + 1;
                continue;
            }

            let a_p1 = quad1.value(u1_prec, v1_prec);
            let a_p2 = quad2.value(u2_prec, v2_prec);
            let a_p_int = 0.5 * (a_p1 + a_p2);

            let an_ip = WLinePnt {
                p3d: a_p_int,
                u1: if is_reverse { u2_prec } else { u1_prec },
                v1: if is_reverse { v2_prec } else { v1_prec },
                u2: if is_reverse { u1_prec } else { u2_prec },
                v2: if is_reverse { v1_prec } else { v2_prec },
            };
            line.insert_before(lp, an_ip);

            a_nb_points += 1;
            a_last_point_index += 1;
            fp = lp + 1;
        }

        if a_nb_points >= min_nb_points {
            return;
        }
        if a_nb_points == a_nb_points_prev {
            return;
        }
    }
}
