//! IntCyCy — cylinder-cylinder intersection — 1:1 translation of OCCT
//! `IntPatch_ImpImpIntersection.cxx`:
//!   - CyCyAnalyticalIntersect (L4845-5196)
//!   - IntCyCy (L7881-8106)
//!
//! rcad data-model notes:
//!   - `IntSurf_Quadric` -> `Quadric`; `IntAna_QuadQuadGeo` -> `QuadQuadGeo`.
//!   - `gp_Pnt2d` -> `[f64; 2]`; `Bnd_Box2d` -> `[f64; 4]` = [u_min, v_min, u_max, v_max].
//!   - `IntPatch_Point` -> `IntPatchPoint` (spnt) / `IntPatchVertex` (GLine vertices).
//!   - `IntSurf_TypeTrans`/`IntSurf_Situation` -> `transitions::{TypeTrans, Situation}`.

use rcad_kernel::geom::Curve3;

use super::cycy_coeffs::new_coeffs;
use super::cycy_common::{BndRange, inscribe_interval};
use super::cycy_walking::cy_cy_no_geometric;
use super::elclib as ecl;
use super::imp_imp_intersection::IntStatus;
use super::quad_quad_geo::{AnaResultType, QuadQuadGeo};
use super::transitions::{Situation, TypeTrans, Transition};
use super::cycy_boundaries::{WorkWithBoundaries, boundaries_computing};
use super::{IntPatchIType, IntPatchLine, IntPatchPoint, IntPatchVertex};
use crate::topalgo::int_surf::quadric::Quadric;

/// OCCT IntPatch_GLine(Elips, Tang, Trans1, Trans2) helper — build an
/// IntPatchLine carrying an ellipse (with transitions and vertices).
fn gline_ellipse(e: rcad_kernel::geom::Ellipse3, tang: bool, t1: TypeTrans, t2: TypeTrans) -> IntPatchLine {
    let mut line =
        IntPatchLine::analytic(IntPatchIType::Ellipse, Curve3::Ellipse(e), [0.0, std::f64::consts::TAU]);
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// OCCT CyCyAnalyticalIntersect (L4845-5196) — post-processes the analytic
/// IntAna_QuadQuadGeo result (Empty/Same/Point/Line/Ellipse) into GLine /
/// spnt entries.  Returns false when the analytic handler declines the result
/// (Circle/Parabola/Hyperbola/NoGeometricSolution), in which case the caller
/// falls through to the numeric engine.
#[allow(clippy::too_many_arguments, unused_assignments)]
fn cy_cy_analytical_intersect(
    quad1: &Quadric,
    quad2: &Quadric,
    inter: &QuadQuadGeo,
    tol: f64,
    empty: &mut bool,
    same: &mut bool,
    multpoint: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let cy1 = quad1.cylinder();
    let cy2 = quad2.cylinder();

    let typint = inter.type_inter();
    let nb_sol = inter.nb_solutions();
    *empty = false;
    *same = false;

    match typint {
        AnaResultType::Empty => {
            *empty = true;
        }
        AnaResultType::Same => {
            *same = true;
        }
        AnaResultType::Point => {
            let psol = inter.point(1);
            let (u1, v1) = quad1.parameters(psol);
            let (u2, v2) = quad2.parameters(psol);
            spnt.push(IntPatchPoint {
                p1: psol,
                p2: psol,
                u1,
                v1,
                u2,
                v2,
                tolerance: tol,
            });
        }
        AnaResultType::Line => {
            if nb_sol == 1 {
                // Cylinders are tangent to each other by line.
                let linsol = inter.line(1);
                let ptref = linsol.origin;

                // Radius-vectors.
                let crb1 = (cy1.origin - ptref).normalize_or_zero();
                let crb2 = (cy2.origin - ptref).normalize_or_zero();

                // Outer normal lines.
                let norm1 = quad1.normale(ptref);
                let norm2 = quad2.normale(ptref);
                let (mut situcyl1, mut situcyl2) = (Situation::Unknown, Situation::Unknown);

                if crb1.dot(crb2) < 0.0 {
                    // Centres of curvature are "opposed".  Normal and
                    // radius-vector of the 1st(!) cylinder are used for judging
                    // what the situation of the 2nd(!) cylinder is.
                    situcyl2 = if norm1.dot(crb1) > 0.0 { Situation::Inside } else { Situation::Outside };
                    situcyl1 = if norm2.dot(crb2) > 0.0 { Situation::Inside } else { Situation::Outside };
                } else {
                    if cy1.radius < cy2.radius {
                        situcyl2 = if norm1.dot(crb1) > 0.0 { Situation::Inside } else { Situation::Outside };
                        situcyl1 = if norm2.dot(crb2) > 0.0 { Situation::Outside } else { Situation::Inside };
                    } else {
                        situcyl2 = if norm1.dot(crb1) > 0.0 { Situation::Outside } else { Situation::Inside };
                        situcyl1 = if norm2.dot(crb2) > 0.0 { Situation::Inside } else { Situation::Outside };
                    }
                }

                let mut glig = IntPatchLine::analytic(
                    IntPatchIType::Line,
                    Curve3::Line(linsol),
                    [f64::NEG_INFINITY, f64::INFINITY],
                );
                glig.trans1 = Some(Transition::new_touch(true, situcyl1, false));
                glig.trans2 = Some(Transition::new_touch(true, situcyl2, false));
                slin.push(glig);
            } else {
                for i in 1..=nb_sol {
                    let linsol = inter.line(i);
                    let ptref = linsol.origin;
                    let lsd = linsol.direction;

                    // Theoretically, qwe = +/- 1.0.
                    let qwe = lsd.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
                    let (t1, t2) = if qwe > 0.00000001 {
                        (TypeTrans::Out, TypeTrans::In)
                    } else if qwe < -0.00000001 {
                        (TypeTrans::In, TypeTrans::Out)
                    } else {
                        (TypeTrans::Undecided, TypeTrans::Undecided)
                    };
                    let mut glig = IntPatchLine::analytic(
                        IntPatchIType::Line,
                        Curve3::Line(linsol),
                        [f64::NEG_INFINITY, f64::INFINITY],
                    );
                    glig.trans1 = Some(Transition::new_in_out(false, t1));
                    glig.trans2 = Some(Transition::new_in_out(false, t2));
                    slin.push(glig);
                }
            }
        }
        AnaResultType::Ellipse => {
            let elipsol = inter.ellipse();

            let pttang1 = ecl::conic_value(&Curve3::Ellipse(elipsol), 0.5 * std::f64::consts::PI);
            let pttang2 = ecl::conic_value(&Curve3::Ellipse(elipsol), 1.5 * std::f64::consts::PI);

            *multpoint = true;
            let mut pmult1 = IntPatchVertex::default();
            pmult1.set_value(pttang1, tol, true);
            pmult1.set_multiple(true);
            let mut pmult2 = IntPatchVertex::default();
            pmult2.set_value(pttang2, tol, true);
            pmult2.set_multiple(true);

            let (o_u1, o_v1) = quad1.parameters(pttang1);
            let (o_u2, o_v2) = quad2.parameters(pttang1);
            pmult1.set_parameters(o_u1, o_v1, o_u2, o_v2);
            let (o_u1, o_v1) = quad1.parameters(pttang2);
            let (o_u2, o_v2) = quad2.parameters(pttang2);
            pmult2.set_parameters(o_u1, o_v1, o_u2, o_v2);

            // Process the first ellipse.

            // Compute the transition of the line.
            let (ptref, tgt) = ecl::conic_d1(&Curve3::Ellipse(elipsol), 0.0);

            // Theoretically, qwe = +/- |Tgt|.
            let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
            let (t1, t2) = if qwe > 0.00000001 {
                (TypeTrans::Out, TypeTrans::In)
            } else if qwe < -0.00000001 {
                (TypeTrans::In, TypeTrans::Out)
            } else {
                (TypeTrans::Undecided, TypeTrans::Undecided)
            };

            // Transition computed at point 0 -> Trans2, Trans1 because here it
            // should be computed at PI.
            let mut glig = gline_ellipse(elipsol, false, t2, t1);

            {
                let a_p = ecl::conic_value(&Curve3::Ellipse(elipsol), 0.0);
                let mut a_ip = IntPatchVertex::default();
                a_ip.set_value(a_p, tol, false);
                a_ip.set_multiple(false);
                let (a_u1, a_v1) = quad1.parameters(a_p);
                let (a_u2, a_v2) = quad2.parameters(a_p);
                a_ip.set_parameters(a_u1, a_v1, a_u2, a_v2);
                a_ip.set_parameter(0.0);
                glig.add_vertex(a_ip.clone());
                glig.set_first_point(1);
                a_ip.set_parameter(2.0 * std::f64::consts::PI);
                glig.add_vertex(a_ip);
                glig.set_last_point(2);
            }

            pmult1.set_parameter(0.5 * std::f64::consts::PI);
            glig.add_vertex(pmult1.clone());
            pmult2.set_parameter(1.5 * std::f64::consts::PI);
            glig.add_vertex(pmult2.clone());

            slin.push(glig);

            // Process the second ellipse.
            let elipsol2 = inter.ellipse_n(2);

            let param1 = ecl::conic_parameter(&Curve3::Ellipse(elipsol2), pttang1);
            let param2 = ecl::conic_parameter(&Curve3::Ellipse(elipsol2), pttang2);
            let mut parampourtransition = 0.0;
            if param1 < param2 {
                pmult1.set_parameter(0.5 * std::f64::consts::PI);
                pmult2.set_parameter(1.5 * std::f64::consts::PI);
                parampourtransition = std::f64::consts::PI;
            } else {
                pmult1.set_parameter(1.5 * std::f64::consts::PI);
                pmult2.set_parameter(0.5 * std::f64::consts::PI);
                parampourtransition = 0.0;
            }

            // Compute the transitions of the line for the second line.
            let (ptref, tgt) = ecl::conic_d1(&Curve3::Ellipse(elipsol2), parampourtransition);

            // Theoretically, qwe = +/- |Tgt|.
            let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
            let (t1, t2) = if qwe > 0.00000001 {
                (TypeTrans::Out, TypeTrans::In)
            } else if qwe < -0.00000001 {
                (TypeTrans::In, TypeTrans::Out)
            } else {
                (TypeTrans::Undecided, TypeTrans::Undecided)
            };

            // The transition was computed at a point of this line.
            let mut glig = gline_ellipse(elipsol2, false, t1, t2);

            {
                let a_p = ecl::conic_value(&Curve3::Ellipse(elipsol2), 0.0);
                let mut a_ip = IntPatchVertex::default();
                a_ip.set_value(a_p, tol, false);
                a_ip.set_multiple(false);
                let (a_u1, a_v1) = quad1.parameters(a_p);
                let (a_u2, a_v2) = quad2.parameters(a_p);
                a_ip.set_parameters(a_u1, a_v1, a_u2, a_v2);
                a_ip.set_parameter(0.0);
                glig.add_vertex(a_ip.clone());
                glig.set_first_point(1);
                a_ip.set_parameter(2.0 * std::f64::consts::PI);
                glig.add_vertex(a_ip);
                glig.set_last_point(2);
            }

            glig.add_vertex(pmult1);
            glig.add_vertex(pmult2);

            slin.push(glig);
        }
        AnaResultType::Parabola | AnaResultType::Hyperbola => {
            panic!("IntCyCy(): Wrong intersection type!");
        }
        // Circle is useful when we will work with trimmed surfaces (two
        // cylinders can be tangent by their basises, e.g. circle).
        AnaResultType::Circle | AnaResultType::PointAndCircle | AnaResultType::NoGeometricSolution => {
            return false;
        }
    }

    true
}

/// OCCT IntCyCy (L7881-8106) — cylinder-cylinder intersection.
///
/// The analytic IntAna_QuadQuadGeo result is post-processed by
/// CyCyAnalyticalIntersect; when it declines (or returns NoGeometricSolution),
/// the general cylinder-cylinder intersection is computed numerically by
/// CyCyNoGeometric.
#[allow(clippy::too_many_arguments)]
pub fn int_cycy(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_3d: f64,
    tol_2d: f64,
    uv1: [f64; 4],
    uv2: [f64; 4],
    is_empty: &mut bool,
    is_same_surface: &mut bool,
    is_multiple_point: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> IntStatus {
    *is_empty = true;
    *is_same_surface = false;
    *is_multiple_point = false;
    slin.clear();
    spnt.clear();

    let a_cyl1 = quad1.cylinder();
    let a_cyl2 = quad2.cylinder();

    let mut an_inter = QuadQuadGeo::new();
    an_inter.perform_cylinder_cylinder(quad1, quad2, tol_3d);

    if !an_inter.is_done() {
        return IntStatus::Fail;
    }

    if an_inter.type_inter() != AnaResultType::NoGeometricSolution {
        if cy_cy_analytical_intersect(
            quad1,
            quad2,
            &an_inter,
            tol_3d,
            is_empty,
            is_same_surface,
            is_multiple_point,
            slin,
            spnt,
        ) {
            return IntStatus::OK;
        }

        // Analytical handler declined: discard anything it may have appended so
        // the numerical path starts from a clean output state.
        *is_empty = true;
        *is_same_surface = false;
        *is_multiple_point = false;
        slin.clear();
        spnt.clear();
    }

    // Here, the intersection line is not an analytical curve
    // (line, circle, ellipse, etc.).

    // aUSBou[0/1][0]=Uf, [1]=Ul; aVSBou likewise for V (filled but unused here,
    // matching OCCT — the V bounds are consumed inside CyCyNoGeometric via the
    // WorkWithBoundaries UV boxes).
    let a_us_bou = [[uv1[0], uv1[2]], [uv2[0], uv2[2]]];
    let _a_vs_bou = [[uv1[1], uv1[3]], [uv2[1], uv2[3]]];

    let a_period = 2.0 * std::f64::consts::PI;
    let a_nb_wlines = 2usize;

    let an_equation_coeffs1 = new_coeffs(&a_cyl1, &a_cyl2);
    let an_equation_coeffs2 = new_coeffs(&a_cyl2, &a_cyl1);

    // Boundaries.  The intersection result can include two non-connected
    // regions (see WorkWithBoundaries::BoundariesComputing).
    let a_nb_of_boundaries = 2usize;
    let mut an_u_range: [[BndRange; 2]; 2] =
        [[BndRange::new(), BndRange::new()], [BndRange::new(), BndRange::new()]];

    if !boundaries_computing(&an_equation_coeffs1, a_period, &mut an_u_range[0]) {
        return IntStatus::OK;
    }
    if !boundaries_computing(&an_equation_coeffs2, a_period, &mut an_u_range[1]) {
        return IntStatus::OK;
    }

    // anURange[*] can be in different periodic regions compared with the
    // First-Last surface (e.g. the surface is a full cylinder [0, 2*PI] but
    // anURange is [5, 7]).  Trivial common-range computation returns [5, 2*PI]
    // and its summary length is 2*PI-5 == 1.28... only — wrong.  This problem
    // can be solved by the following algorithm:
    //  1. split anURange[*] by the surface boundary;
    //  2. shift every new range to inscribe it in [Ufirst, Ulast] of the cylinder;
    //  3. consider only common ranges between [Ufirst, Ulast] and new ranges.
    let mut a_sum_range = [0.0, 0.0];
    for a_cid in 0..2usize {
        let mut a_list_of_rng = vec![an_u_range[a_cid][0], an_u_range[a_cid][1]];
        let a_split_arr = [a_us_bou[a_cid][0], a_us_bou[a_cid][1], 0.0];
        for a_s_ind in 0..3 {
            let a_lst_temp = a_list_of_rng.clone();
            a_list_of_rng.clear();
            for a_rng in a_lst_temp {
                a_rng.split(a_split_arr[a_s_ind], &mut a_list_of_rng, a_period);
            }
        }
        for mut a_curr_range in a_list_of_rng {
            let mut a_bound_r = BndRange::new();
            a_bound_r.add(a_us_bou[a_cid][0]);
            a_bound_r.add(a_us_bou[a_cid][1]);

            if !inscribe_interval(
                a_us_bou[a_cid][0],
                a_us_bou[a_cid][1],
                &mut a_curr_range,
                tol_2d,
                a_period,
            ) {
                // If aCurrRange does not have a common block with [Ufirst, Ulast]
                // of the cylinder then try to inscribe [Ufirst, Ulast] in the
                // boundaries of aCurrRange.
                let (a_f, a_l) = match a_curr_range.get_bounds() {
                    Some(b) => b,
                    None => continue,
                };
                if a_l < a_us_bou[a_cid][0] {
                    a_curr_range.shift(a_period);
                } else if a_f > a_us_bou[a_cid][1] {
                    a_curr_range.shift(-a_period);
                }
            }

            a_bound_r.common(&a_curr_range);

            let a_delta = a_bound_r.delta();
            if a_delta > 0.0 {
                a_sum_range[a_cid] += a_delta;
            }
        }
    }

    // The bigger range the bigger number of points in the Walking-line (WLine)
    // we will be able to add and consequently the more precise the intersection
    // line.  Every point of the WLine is determined as a function of the
    // U1-parameter, where U1 is the U-parameter on the 1st quadric.  Therefore,
    // we should use the quadric with the bigger range as the 1st parameter in
    // IntCyCy().  On the other hand, there is no point in reversing in case of
    // an analytical intersection (when the result is a line, ellipse, point...);
    // this result is independent of the arguments order.
    let is_to_reverse = a_sum_range[1] > a_sum_range[0];

    if is_to_reverse {
        let a_bound_work = WorkWithBoundaries::new(
            quad2,
            quad1,
            &an_equation_coeffs2,
            uv2,
            uv1,
            a_nb_wlines,
            a_period,
            tol_3d,
            tol_2d,
            true,
        );
        return cy_cy_no_geometric(
            &a_cyl2,
            &a_cyl1,
            &a_bound_work,
            &mut an_u_range[1],
            a_nb_of_boundaries,
            is_empty,
            slin,
            spnt,
        );
    } else {
        let a_bound_work = WorkWithBoundaries::new(
            quad1,
            quad2,
            &an_equation_coeffs1,
            uv1,
            uv2,
            a_nb_wlines,
            a_period,
            tol_3d,
            tol_2d,
            false,
        );
        return cy_cy_no_geometric(
            &a_cyl1,
            &a_cyl2,
            &a_bound_work,
            &mut an_u_range[0],
            a_nb_of_boundaries,
            is_empty,
            slin,
            spnt,
        );
    }
}
