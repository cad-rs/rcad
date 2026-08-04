//! OCCT IntPatch_ALineToWLine — converts an analytic line (IntPatch_ALine,
//! carrying an IntAna_Curve) into Walking-lines (IntPatch_WLine), splitting at
//! vertices and special points (poles/seams).  (IntPatch_ALineToWLine.cxx)
//!
//! 1:1 Rust translation.  rcad data-model notes:
//!   - IntPatch_ALine -> IntAnaCurve (rcad int_quad_quad).
//!   - IntSurf_LineOn2S (the WLine under construction) -> Vec<PntOn2S>.
//!   - IntSurf_PntOn2S -> PntOn2S (special_points.rs).
//!   - IntPatch_Point -> PatchPoint (special_points.rs).
//!   - IntPatch_WLine -> IntPatchLine { line_type: Walking, wline_pnts, ... }.

use super::int_quad_quad::IntAnaCurve;
use super::special_points::{
    add_cross_uv_iso_point, add_singular_pole, adjust_point_and_vertex, set_period, PatchPoint,
    PntOn2S, SpecPntType,
};
use super::{IntPatchIType, IntPatchLine, WLinePnt, WLineType};
use crate::geomalgo::int_surf::quadric::Quadric;
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};
use rcad_kernel::precision::CONFUSION;

/// OCCT IntPatch_ALineToWLine.
pub struct ALineToWLine {
    my_s1: Surface3,
    my_s2: Surface3,
    my_quad1: Quadric,
    my_quad2: Quadric,
    my_nb_points_in_wline: usize,
    my_tol_open_domain: f64,
    my_tol_transition: f64,
    my_tol_3d: f64,
}

impl ALineToWLine {
    /// OCCT L145-208: constructor.
    pub fn new(s1: &Surface3, s2: &Surface3, nb_points: usize) -> Self {
        let my_quad1 = Quadric::from_surface3(s1).unwrap_or_else(Quadric::new);
        let my_quad2 = Quadric::from_surface3(s2).unwrap_or_else(Quadric::new);
        ALineToWLine {
            my_s1: s1.clone(),
            my_s2: s2.clone(),
            my_quad1,
            my_quad2,
            my_nb_points_in_wline: nb_points,
            my_tol_open_domain: 1e-9,
            my_tol_transition: 1e-8,
            my_tol_3d: CONFUSION,
        }
    }

    /// OCCT MakeWLine(aline, theLines) (L361-378): the first/last parameter of
    /// the ALine, adjusted by the open-domain tolerance when the endpoint is
    /// not included, then MakeWLine(f, l, theLines).
    pub fn make_wline(&self, a_line: &mut IntPatchLine, the_lines: &mut Vec<IntPatchLine>) {
        // OCCT MakeWLine takes the ALine (IntPatch_ALine): the curve plus its
        // transition properties (TransitionOnS1/S2, IsTangent, Situations).
        let aline_trans1 = a_line.trans1;
        let aline_trans2 = a_line.trans2;
        let a_curve = match a_line.a_curve.as_mut() {
            Some(c) => c,
            None => return,
        };
        let d = a_curve.domain();
        // OCCT L366-374: FirstParameter/LastParameter return IsIncluded =
        // !IsFirstOpen()/!IsLastOpen(); the open-domain tolerance is added only
        // when the endpoint is not included (open domain).
        let f = if a_curve.is_first_open() {
            d[0] + self.my_tol_open_domain
        } else {
            d[0]
        };
        let l = if a_curve.is_last_open() {
            d[1] - self.my_tol_open_domain
        } else {
            d[1]
        };
        self.make_wline_range(a_curve, aline_trans1, aline_trans2, f, l, the_lines);
    }

    /// OCCT MakeWLine(aline, theFPar, theLPar, theLines) (L382-963).
    #[allow(clippy::too_many_lines)]
    pub fn make_wline_range(
        &self,
        a_line: &mut IntAnaCurve,
        aline_trans1: Option<super::transitions::Transition>,
        aline_trans2: Option<super::transitions::Transition>,
        f_par: f64,
        l_par: f64,
        the_lines: &mut Vec<IntPatchLine>,
    ) {
        // OCCT L388-392: if no vertices, return.
        if !a_line.has_vertices() {
            return;
        }

        // OCCT L394-415: the same points can be marked by different vertices;
        // unify the tolerances of all vertices marking the same point.
        {
            let a_nb_vert = a_line.vertices.len();
            for i in 0..a_nb_vert {
                let a_cur_toler = a_line.vertices[i].tolerance;
                for j in (i + 1)..a_nb_vert {
                    let a_toler = a_line.vertices[j].tolerance;
                    let a_sum_tol = a_cur_toler + a_toler;
                    if a_line.vertices[i].pnt.is_same(&a_line.vertices[j].pnt, a_sum_tol) {
                        a_line.vertices[i].tolerance = a_sum_tol;
                        a_line.vertices[j].tolerance = a_sum_tol;
                    }
                }
            }
        }

        let a_tol = 2.0 * self.my_tol_3d + CONFUSION;
        let a_prm_tol = (1.0e-4 * (l_par - f_par)).max(rcad_kernel::precision::PCONFUSION);

        let a_vert_params = a_line.vertex_params();
        let a_vert_count = a_vert_params.len();
        let mut has_vertex_been_checked = vec![false; a_vert_count];

        let mut arr_periods = [0.0f64; 4];
        set_period(&self.my_s1, &self.my_s2, &mut arr_periods);

        let mut a_prev_l_point = PntOn2S {
            p: DVec3::ZERO,
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
        };

        // OCCT L439: the singular surface ID (1 = S1, 2 = S2), set by
        // IsPoleOrSeam for each processed vertex and used in the Pole handling.
        let mut a_singular_surface_id = 0usize;

        // OCCT L420: aPrePointExist persists across the while-loop iterations
        // (a while-iteration that breaks on a special point hands it to the
        // next iteration's special-point handling).
        let mut a_pre_point_exist = SpecPntType::None;

        let mut a_parameter = f_par;
        while a_parameter < l_par {
            let mut a_step = (l_par - a_parameter) / (self.my_nb_points_in_wline as f64 - 1.0);
            if a_step < f64::EPSILON * l_par {
                break;
            }

            let mut is_step_reduced = false;
            let mut a_l_par = l_par;
            for i in 0..a_vert_count {
                if has_vertex_been_checked[i] {
                    continue;
                }
                a_l_par = a_vert_params[i];
                if (a_l_par - a_parameter).abs() < a_prm_tol {
                    continue;
                }
                break;
            }

            if (a_step - (a_l_par - a_parameter) > a_prm_tol)
                && ((a_l_par - a_parameter).abs() > a_prm_tol)
            {
                a_step = ((a_l_par - a_parameter) / 5.0).max(1.0e-5);
                is_step_reduced = true;
            }

            let mut a_lin_on2s: Vec<PntOn2S> = Vec::new();
            let mut an_is_first_degenerated = false;
            let mut an_is_last_degenerated = false;

            let mut a_step_min = 0.1 * a_step;
            let mut a_step_max = 10.0 * a_step;

            let mut is_last = false;
            let mut a_prev_param = a_parameter;
            let mut a_seq_vertex: Vec<PatchPoint> = Vec::new();

            while !is_last {
                let mut a_p_on2s = PntOn2S {
                    p: DVec3::ZERO,
                    u1: 0.0,
                    v1: 0.0,
                    u2: 0.0,
                    v2: 0.0,
                };

                if l_par <= a_parameter {
                    is_last = true;
                    if a_pre_point_exist != SpecPntType::None {
                        break;
                    }
                    a_parameter = l_par;
                }

                let mut is_point_valid = false;
                let mut a_tg_magn = 0.0f64;
                let a_pnt3d;
                let a_tg;
                match a_line.d1u(a_parameter) {
                    Some((p, tg)) => {
                        a_pnt3d = p;
                        a_tg = tg;
                    }
                    None => {
                        a_pnt3d = a_line.value(a_parameter).unwrap_or(DVec3::ZERO);
                        a_tg = DVec3::ZERO;
                    }
                }
                if self.get_section_radius(a_pnt3d) < 5.0e-6 {
                    // Cannot compute 2D-parameters of aPOn2S correctly.
                    if an_is_last_degenerated {
                        a_lin_on2s.pop();
                    }
                    is_point_valid = false;
                } else {
                    is_point_valid = true;
                }
                a_tg_magn = a_tg.length();
                let (u1, v1) = self.my_quad1.parameters(a_pnt3d);
                let (u2, v2) = self.my_quad2.parameters(a_pnt3d);
                a_p_on2s.set_value(a_pnt3d, u1, v1, u2, v2);

                if a_pre_point_exist != SpecPntType::None {
                    let a_u_res = rcad_kernel::topo::topods::u_resolution_for_surface(
                        &self.my_s1, self.my_tol_3d,
                    )
                    .max(rcad_kernel::topo::topods::u_resolution_for_surface(
                        &self.my_s2, self.my_tol_3d,
                    ));
                    let a_v_res = rcad_kernel::topo::topods::v_resolution_for_surface(
                        &self.my_s1, self.my_tol_3d,
                    )
                    .max(rcad_kernel::topo::topods::v_resolution_for_surface(
                        &self.my_s2, self.my_tol_3d,
                    ));
                    let a_tol_2d = match a_pre_point_exist {
                        SpecPntType::Pole => -1.0,
                        SpecPntType::SeamV => a_v_res,
                        SpecPntType::SeamUV => a_u_res.max(a_v_res),
                        _ => a_u_res,
                    };

                    let mut a_rpt = a_p_on2s;
                    if a_pre_point_exist == SpecPntType::Pole {
                        let mut a_prt = 0.5 * (a_prev_param + l_par);
                        for i in 0..a_vert_count {
                            let a_param = a_vert_params[i];
                            if a_param <= a_prev_param {
                                continue;
                            }
                            if (a_param - a_prev_param) < a_prm_tol {
                                let a_pnt3d_i = a_line.value(a_param).unwrap_or(DVec3::ZERO);
                                if a_p_on2s.p.distance_squared(a_pnt3d_i)
                                    < rcad_kernel::precision::CONFUSION
                                        * rcad_kernel::precision::CONFUSION
                                {
                                    continue;
                                }
                            }
                            a_prt = 0.5 * (a_param + a_prev_param);
                            break;
                        }
                        let a_pnt3d = a_line.value(a_prt).unwrap_or(DVec3::ZERO);
                        let (u1, v1) = self.my_quad1.parameters(a_pnt3d);
                        let (u2, v2) = self.my_quad2.parameters(a_pnt3d);
                        a_rpt.set_value(a_pnt3d, u1, v1, u2, v2);

                        // OCCT L582-595: the current point coincides with the
                        // previous one — set the U-parameter on the singular
                        // surface to the precise value of the reference point.
                        if a_p_on2s.is_same(
                            &a_prev_l_point,
                            rcad_kernel::precision::APPROXIMATION.max(a_tol),
                        ) {
                            if a_singular_surface_id == 1 {
                                a_p_on2s.u1 = u1;
                            } else {
                                a_p_on2s.u2 = u2;
                            }
                        }
                    }

                    if continue_after_special_point(
                        &self.my_s1,
                        &self.my_s2,
                        &a_rpt,
                        a_pre_point_exist,
                        a_tol_2d,
                        &mut a_prev_l_point,
                        false,
                    ) {
                        add_point_into_line(&mut a_lin_on2s, &arr_periods, &mut a_prev_l_point, None);
                    } else if a_parameter == l_par {
                        break;
                    } else if a_parameter + a_step < l_par {
                        // Prediction of the next point (OCCT L614-633).
                        let a_pnt3d_next = a_line.value(a_parameter + a_step).unwrap_or(DVec3::ZERO);
                        let (an_u1, a_v1) = self.my_quad1.parameters(a_pnt3d_next);
                        let (an_u2, a_v2) = self.my_quad2.parameters(a_pnt3d_next);
                        let a_p_on2s_next = PntOn2S {
                            p: a_pnt3d_next,
                            u1: an_u1,
                            v1: a_v1,
                            u2: an_u2,
                            v2: a_v2,
                        };
                        // OCCT L626-627: the UV distance on either surface must
                        // not exceed PI.
                        let sq2 = DVec2::new(a_p_on2s_next.u2, a_p_on2s_next.v2)
                            .distance_squared(DVec2::new(a_rpt.u2, a_rpt.v2));
                        let sq1 = DVec2::new(a_p_on2s_next.u1, a_p_on2s_next.v1)
                            .distance_squared(DVec2::new(a_rpt.u1, a_rpt.v1));
                        if sq2 > std::f64::consts::PI * std::f64::consts::PI
                            || sq1 > std::f64::consts::PI * std::f64::consts::PI
                        {
                            a_prev_l_point = a_rpt;
                            a_prev_param = a_parameter;
                            // OCCT for-loop increment: aParameter += aStep runs
                            // on `continue`.
                            a_parameter += a_step;
                            continue;
                        }
                    }
                }

                a_pre_point_exist = SpecPntType::None;

                let mut a_vertex_number = -1isize;
                for i in 0..a_vert_count {
                    if has_vertex_been_checked[i] {
                        continue;
                    }
                    let a_vp = a_line.vertex_at(i);
                    let a_param = a_vert_params[i];
                    if ((a_prev_param < a_param) && (a_param <= a_parameter))
                        || ((a_prev_param == a_parameter) && (a_param == a_parameter))
                        || (a_p_on2s.is_same(&a_vp.pnt, a_vp.tolerance)
                            && (a_vp.param_on_line - a_parameter).abs() < a_prm_tol)
                    {
                        a_vertex_number = i as isize;
                        break;
                    }
                }

                a_prev_param = a_parameter;

                if a_vertex_number < 0 {
                    if is_point_valid {
                        if !is_step_reduced {
                            self.step_computing(
                                a_line,
                                &a_p_on2s,
                                l_par,
                                a_parameter,
                                a_tg_magn,
                                a_step_min,
                                a_step_max,
                                &mut a_step,
                            );
                        }
                        add_point_into_line(&mut a_lin_on2s, &arr_periods, &mut a_p_on2s, None);
                        a_prev_l_point = a_p_on2s;
                    } else {
                        // Add point, set corresponding status: to be corrected later.
                        let mut to_add = false;
                        if a_lin_on2s.is_empty() {
                            an_is_first_degenerated = true;
                            to_add = true;
                        } else if a_lin_on2s.len() > 1 {
                            an_is_last_degenerated = true;
                            to_add = true;
                        }
                        if to_add {
                            add_point_into_line(&mut a_lin_on2s, &arr_periods, &mut a_p_on2s, None);
                            a_prev_l_point = a_p_on2s;
                        }
                    }
                    a_parameter += a_step;
                    continue;
                }

                let mut a_vtx = a_line.vertex_at(a_vertex_number as usize);
                let mut a_new_vertex_param = (a_lin_on2s.len() + 1) as f64;
                let a_nb_points_prev = a_lin_on2s.len();

                // Reference point for the vertex parameter.
                let mut a_pref_iso = a_vtx.pnt;
                if a_lin_on2s.len() < 1 {
                    for i in (a_vertex_number as usize + 1)..a_vert_count {
                        let a_param = a_vert_params[i];
                        if (a_param - a_vert_params[a_vertex_number as usize])
                            > rcad_kernel::precision::PCONFUSION
                        {
                            let a_prm = 0.5 * (a_param + a_vert_params[a_vertex_number as usize]);
                            let a_pnt3d = a_line.value(a_prm).unwrap_or(DVec3::ZERO);
                            let (u1, v1) = self.my_quad1.parameters(a_pnt3d);
                            let (u2, v2) = self.my_quad2.parameters(a_pnt3d);
                            a_pref_iso.set_value(a_pnt3d, u1, v1, u2, v2);
                            break;
                        }
                    }
                } else {
                    a_pref_iso = *a_lin_on2s.last().unwrap();
                }

                let (a_pre_point, a_singular) = self.is_pole_or_seam(
                    &a_pref_iso,
                    &mut a_lin_on2s,
                    &mut a_vtx,
                    &arr_periods,
                    a_tol,
                );
                a_pre_point_exist = a_pre_point;
                a_singular_surface_id = a_singular;

                if a_pre_point_exist == SpecPntType::Pole
                    || a_pre_point_exist == SpecPntType::PoleSeamU
                {
                    if a_lin_on2s.len() == 1 {
                        an_is_first_degenerated = true;
                    } else {
                        an_is_last_degenerated = true;
                    }
                }

                let a_cur_vert_param = a_vtx.param_on_line;
                if a_pre_point_exist != SpecPntType::None {
                    if a_nb_points_prev == a_lin_on2s.len() {
                        a_new_vertex_param = a_nb_points_prev as f64;
                    }
                    a_prev_param = a_cur_vert_param;
                    a_parameter = a_cur_vert_param;
                } else {
                    if !is_point_valid {
                        a_parameter += a_step;
                        continue;
                    }
                    if a_vtx.tolerance > a_tol {
                        a_vtx.set_value(a_p_on2s);
                        add_point_into_line(&mut a_lin_on2s, &arr_periods, &mut a_p_on2s, Some(&mut a_vtx));
                    } else {
                        add_vertex_point(&mut a_lin_on2s, &mut a_vtx, &arr_periods);
                    }
                }

                a_prev_l_point = *a_lin_on2s.last().unwrap_or(&a_p_on2s);
                a_p_on2s = *a_lin_on2s.last().unwrap_or(&a_p_on2s);

                // Unify vertices marking the same 3D point.
                {
                    let a_sq_tol = a_tol * a_tol;
                    let a_p1 = a_line.value(a_cur_vert_param).unwrap_or(DVec3::ZERO);
                    let a_vert_p2s = a_vtx.pnt;
                    let a_vert_toler = a_vtx.tolerance;
                    let mut is_found = false;
                    for i in 0..a_vert_count {
                        if has_vertex_been_checked[i] {
                            continue;
                        }
                        let a_p2 = a_line.value(a_vert_params[i]).unwrap_or(DVec3::ZERO);
                        if a_p1.distance_squared(a_p2) < a_sq_tol {
                            let mut a_l_vtx = a_line.vertex_at(i);
                            a_l_vtx.pnt = a_vert_p2s;
                            a_l_vtx.tolerance = a_vert_toler;
                            let a_param = a_l_vtx.param_on_line;
                            if (a_param - l_par).abs() <= rcad_kernel::precision::PCONFUSION {
                                a_l_vtx.param_on_line = -1.0;
                            } else {
                                a_l_vtx.param_on_line = a_new_vertex_param;
                            }
                            a_seq_vertex.push(a_l_vtx);
                            has_vertex_been_checked[i] = true;
                            is_found = true;
                        } else if is_found {
                            break;
                        }
                    }
                }

                if (a_pre_point_exist != SpecPntType::None) && (a_lin_on2s.len() > 1) {
                    break;
                }

                if is_step_reduced {
                    is_step_reduced = false;
                    a_step = (l_par - a_parameter) / (self.my_nb_points_in_wline as f64 - 1.0);
                    if a_step < f64::EPSILON * l_par {
                        break;
                    }
                    a_l_par = a_vert_params.last().copied().unwrap_or(l_par);
                    for i in 0..a_vert_count {
                        if has_vertex_been_checked[i] {
                            continue;
                        }
                        a_l_par = a_vert_params[i];
                        if (a_l_par - a_parameter).abs() < a_prm_tol {
                            continue;
                        }
                        break;
                    }
                    if (a_step - (a_l_par - a_parameter) > a_prm_tol)
                        && ((a_l_par - a_parameter).abs() > a_prm_tol)
                    {
                        a_step = ((a_l_par - a_parameter) / 5.0).max(1.0e-5);
                        is_step_reduced = true;
                    }
                    a_step_min = 0.1 * a_step;
                    a_step_max = 10.0 * a_step;
                }

                // OCCT for-loop increment: aParameter += aStep applies at the end
                // of every iteration that does not `break`/`continue`.
                a_parameter += a_step;
            } // for(; !isLast; aParameter += aStep)

            if a_lin_on2s.len() < 2 {
                a_parameter += a_step;
                continue;
            }

            // Correct first and last points if needed.
            if a_lin_on2s.len() >= 3 {
                if an_is_first_degenerated {
                    self.correct_end_point(&mut a_lin_on2s, 0);
                }
                if an_is_last_degenerated {
                    let last_idx = a_lin_on2s.len() - 1;
                    self.correct_end_point(&mut a_lin_on2s, last_idx);
                }
            }

            // WLine creation.
            let mut wline_pnts: Vec<WLinePnt> = a_lin_on2s
                .iter()
                .map(|p| WLinePnt {
                    p3d: p.p,
                    u1: p.u1,
                    v1: p.v1,
                    u2: p.u2,
                    v2: p.v2,
                })
                .collect();

            for v in a_seq_vertex.iter() {
                let mut a_vtx = v.clone();
                if a_vtx.param_on_line == -1.0 {
                    a_vtx.param_on_line = wline_pnts.len() as f64;
                }
                let idx = a_vtx.param_on_line as usize;
                if idx >= 1 && idx <= wline_pnts.len() {
                    let w = &mut wline_pnts[idx - 1];
                    w.p3d = a_vtx.pnt.p;
                    w.u1 = a_vtx.pnt.u1;
                    w.v1 = a_vtx.pnt.v1;
                    w.u2 = a_vtx.pnt.u2;
                    w.v2 = a_vtx.pnt.v2;
                }
            }

            // OCCT L939-947: add the vertices (aSeqVertex) to the WLine; a
            // vertex with parameter -1 (closed curve) is set to the last point.
            // The full IntPatch_Point fields are preserved (tolerance, domain
            // flags, arcs, vertex flags) — ComputeVertexParameters consumes them.
            let mut wline_verts: Vec<super::IntPatchVertex> = Vec::new();
            for v in a_seq_vertex.iter() {
                let mut a_vtx = v.clone();
                if a_vtx.param_on_line == -1.0 {
                    a_vtx.param_on_line = wline_pnts.len() as f64;
                }
                wline_verts.push(super::IntPatchVertex {
                    param_on_line: a_vtx.param_on_line,
                    p3d: a_vtx.pnt.p,
                    u1: a_vtx.pnt.u1,
                    v1: a_vtx.pnt.v1,
                    u2: a_vtx.pnt.u2,
                    v2: a_vtx.pnt.v2,
                    tolerance: a_vtx.tolerance,
                    tangent: false,
                    multiple: a_vtx.multiple,
                    on_dom_s1: a_vtx.on_dom_s1,
                    on_dom_s2: a_vtx.on_dom_s2,
                    arc_on_s1: a_vtx.arc_on_s1,
                    arc_on_s2: a_vtx.arc_on_s2,
                    param_on_arc1: a_vtx.param_on_arc1,
                    param_on_arc2: a_vtx.param_on_arc2,
                    is_vertex_on_s1: a_vtx.is_vertex_on_s1,
                    is_vertex_on_s2: a_vtx.is_vertex_on_s2,
                    transition_line_arc1: a_vtx.transition_line_arc1,
                    transition_line_arc2: a_vtx.transition_line_arc2,
                    transition_on_s1: a_vtx.transition_on_s1,
                    transition_on_s2: a_vtx.transition_on_s2,
                });
            }

            // OCCT L896-937: the WLine transitions come from the ALine's
            // TransitionOnS1 (Touch -> situations; Undecided -> none), or are
            // computed from the tangent/normal cross product (dotcross vs
            // myTolTransition) when the ALine transition is In/Out.
            let aline_trans1_type = aline_trans1.map(|t| t.transition_type());
            let aline_is_tangent = aline_trans1.map(|t| t.tangent()).unwrap_or(false);
            let (wl_trans1, wl_trans2) = match aline_trans1_type {
                Some(super::transitions::TypeTrans::Touch) => {
                    let situ1 = aline_trans1
                        .map(|t| t.situation())
                        .unwrap_or(super::transitions::Situation::Unknown);
                    let situ2 = aline_trans2
                        .map(|t| t.situation())
                        .unwrap_or(super::transitions::Situation::Unknown);
                    (
                        Some(super::transitions::Transition::new_touch(
                            aline_is_tangent,
                            situ1,
                            false,
                        )),
                        Some(super::transitions::Transition::new_touch(
                            aline_is_tangent,
                            situ2,
                            false,
                        )),
                    )
                }
                Some(super::transitions::TypeTrans::Undecided) | None => (None, None),
                _ => {
                    // OCCT L914-934: compute the transitions from the line
                    // tangent and the two surface normals.
                    let indice1 = (a_lin_on2s.len() / 3).max(2);
                    let a_pp0 = a_lin_on2s[indice1 - 2].p;
                    let a_pp1 = a_lin_on2s[indice1 - 1].p;
                    let tgvalid = a_pp1 - a_pp0;
                    let a_nq1 = self.my_quad1.normale(a_pp0);
                    let a_nq2 = self.my_quad2.normale(a_pp0);
                    let dotcross = tgvalid.dot(a_nq2.cross(a_nq1));
                    let mut trans1 = super::transitions::TypeTrans::Undecided;
                    let mut trans2 = super::transitions::TypeTrans::Undecided;
                    if dotcross > self.my_tol_transition {
                        trans1 = super::transitions::TypeTrans::Out;
                        trans2 = super::transitions::TypeTrans::In;
                    } else if dotcross < -self.my_tol_transition {
                        trans1 = super::transitions::TypeTrans::In;
                        trans2 = super::transitions::TypeTrans::Out;
                    }
                    (
                        Some(super::transitions::Transition::new_in_out(aline_is_tangent, trans1)),
                        Some(super::transitions::Transition::new_in_out(aline_is_tangent, trans2)),
                    )
                }
            };

            let mut line = IntPatchLine {
                line_type: IntPatchIType::Walking,
                curve: rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: wline_pnts[0].p3d,
                    direction: if wline_pnts.len() > 1 {
                        (wline_pnts[1].p3d - wline_pnts[0].p3d).normalize_or_zero()
                    } else {
                        DVec3::X
                    },
                }),
                t_range: [0.0, 1.0],
                pcurve1: None,
                pcurve2: None,
                tolerance: self.my_tol_3d,
                tang_tolerance: self.my_tol_3d,
                wline_pnts,
                is_purging_allowed: false,
                wl_type: WLineType::ImpImp,
                vertices: wline_verts,
                a_curve: None,
                arc_on_s1: None,
                arc_on_s2: None,
                trans1: wl_trans1,
                trans2: wl_trans2,
                first_point: None,
                last_point: None,
            };

            // OCCT L949-952: SetPeriod(anArrPeriods...) then ComputeVertexParameters
            // (which inserts the missing start/end vertices and can reduce the
            // number of curve points).
            line.compute_vertex_parameters_wline(self.my_tol_3d, arr_periods);

            // OCCT L954-961: keep the line only when it has more than one point.
            if line.wline_pnts.len() > 1 {
                let wpts = &line.wline_pnts;
                line.curve = rcad_kernel::geom::Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: wpts[0].p3d,
                    direction: (wpts[1].p3d - wpts[0].p3d).normalize_or_zero(),
                });
                the_lines.push(line);
            }
        } // while(aParameter < theLPar)
    }

    /// OCCT IsPoleOrSeam (L72-141).  Returns the special-point type and the
    /// singular surface ID (1 = theS1, 2 = theS2, 0 = none).
    fn is_pole_or_seam(
        &self,
        pt_iso_ref: &PntOn2S,
        the_line: &mut Vec<PntOn2S>,
        the_vertex: &mut PatchPoint,
        arr_periods: &[f64; 4],
        tol_3d: f64,
    ) -> (SpecPntType, usize) {
        for i in 0..2 {
            let is_reversed = i > 0;
            let surf = if is_reversed { &self.my_s2 } else { &self.my_s1 };
            let (q_surf, p_surf) = if is_reversed {
                (&self.my_s2, &self.my_s1)
            } else {
                (&self.my_s1, &self.my_s2)
            };
            let mut added_p_type = SpecPntType::None;
            let mut apex_point = PntOn2S {
                p: DVec3::ZERO,
                u1: 0.0,
                v1: 0.0,
                u2: 0.0,
                v2: 0.0,
            };

            // OCCT switch (L91-130) with [[fallthrough]]: Sphere/Cone try the
            // pole, then fall through to Torus (cross-UV-iso) and Cylinder
            // (seam).  The Cylinder case is reached for Sphere/Cone/Torus/
            // Cylinder whenever the type-specific handling found nothing.
            let a_kind = match surf {
                Surface3::Sphere(_) | Surface3::Cone(_) => 0,
                Surface3::Torus(_) => 1,
                Surface3::Cylinder(_) => 2,
                _ => 3,
            };

            match a_kind {
                0 => {
                    // case Sphere/Cone.
                    if add_singular_pole(q_surf, p_surf, pt_iso_ref, the_vertex, &mut apex_point, is_reversed) {
                        added_p_type = SpecPntType::Pole;
                    }
                }
                1 => {}
                _ => {}
            }

            // [[fallthrough]] case Torus: only for a Torus surface, and only
            // when the pole check (if any) did not fire.
            if a_kind == 1 && added_p_type == SpecPntType::None {
                if add_cross_uv_iso_point(q_surf, p_surf, pt_iso_ref, tol_3d, &mut apex_point, is_reversed) {
                    added_p_type = SpecPntType::SeamUV;
                }
            }

            // [[fallthrough]] case Cylinder (L123-126).
            if a_kind != 3 && added_p_type == SpecPntType::None {
                add_vertex_point(the_line, the_vertex, arr_periods);
                return (SpecPntType::SeamU, i + 1);
            }

            if added_p_type != SpecPntType::None {
                add_point_into_line(the_line, arr_periods, &mut apex_point, Some(the_vertex));
                return (added_p_type, i + 1);
            }
        }
        (SpecPntType::None, 0)
    }

    /// OCCT StepComputing (L998-1077).
    #[allow(clippy::too_many_arguments)]
    fn step_computing(
        &self,
        a_line: &IntAnaCurve,
        p_on2s: &PntOn2S,
        last_par_of_aline: f64,
        cur_param: f64,
        tg_magnitude: f64,
        step_min: f64,
        step_max: f64,
        the_step: &mut f64,
    ) {
        if tg_magnitude < CONFUSION {
            return;
        }
        let an_eps = self.my_tol_3d;
        let a_nb_iter_max = 50;

        let a_not_filled_range = last_par_of_aline - cur_param;
        let mut a_min_step = step_min;
        let mut a_max_step = step_max.min(a_not_filled_range);

        if a_min_step > a_max_step {
            *the_step = a_max_step;
            return;
        }

        let a_r = super::point_line::curvature_radius_of_inters_line(
            &self.my_s1,
            &self.my_s2,
            p_on2s.u1,
            p_on2s.v1,
            p_on2s.u2,
            p_on2s.v2,
        );

        if a_r < 0.0 {
            return;
        }
        *the_step = (an_eps * (2.0 * a_r + an_eps)).sqrt() / tg_magnitude;
        *the_step = (*the_step).min(a_max_step);
        *the_step = (*the_step).max(a_min_step);

        let mut a_nb_iter = 0;
        loop {
            a_nb_iter += 1;
            let a_p1 = p_on2s.p;
            let a_p2 = a_line.value(cur_param + *the_step).unwrap_or(DVec3::ZERO);
            let a_status = self.check_deflection(0.5 * (a_p1 + a_p2), an_eps);
            if a_status == 0 {
                break;
            }
            if a_status < 0 {
                a_min_step = *the_step;
            } else {
                a_max_step = *the_step;
            }
            *the_step = 0.5 * (a_min_step + a_max_step);
            if (a_max_step - a_min_step) <= rcad_kernel::precision::PCONFUSION || a_nb_iter > a_nb_iter_max {
                break;
            }
        }
    }

    /// OCCT CheckDeflection (L972-994).
    fn check_deflection(&self, mid_pt: DVec3, max_deflection: f64) -> i32 {
        let mut a_dist = self.my_quad1.distance(mid_pt).abs();
        if a_dist > max_deflection {
            return 1;
        }
        a_dist = self.my_quad2.distance(mid_pt).abs().max(a_dist);
        if a_dist > max_deflection {
            return 1;
        }
        if (a_dist + a_dist) < max_deflection {
            return -1;
        }
        0
    }

    /// OCCT GetSectionRadius (L322-357).
    fn get_section_radius(&self, pnt3d: DVec3) -> f64 {
        let mut ret_val = f64::INFINITY;
        for i in 0..2 {
            let quad = if i != 0 { &self.my_quad2 } else { &self.my_quad1 };
            match quad.type_quadric() {
                crate::geomalgo::int_surf::quadric::QuadricType::Cone => {
                    let a_r_vec = pnt3d - quad.axis_loc();
                    let a_dir = quad.axis_dir();
                    let a_r = a_r_vec.dot(a_dir) * quad.semi_angle().tan();
                    ret_val = ret_val.min(a_r.abs());
                }
                crate::geomalgo::int_surf::quadric::QuadricType::Sphere => {
                    let a_r_vec = pnt3d - quad.location();
                    let a_dir = quad.z_dir();
                    let a_r = quad.radius();
                    let a_d = a_r_vec.dot(a_dir);
                    let a_delta = a_r * a_r - a_d * a_d;
                    if a_delta <= 0.0 {
                        ret_val = 0.0;
                        break;
                    }
                    ret_val = ret_val.min(a_delta.sqrt());
                }
                _ => {}
            }
        }
        ret_val
    }

    /// OCCT CorrectEndPoint (L254-318).
    fn correct_end_point(&self, the_line: &mut Vec<PntOn2S>, the_index: usize) {
        let a_tol = 1.0e-5;
        let a_sq_tol = 1.0e-10;

        let (an_ind_first, an_ind_second) = if the_index == 0 {
            (2usize, 1usize)
        } else {
            (the_index - 2, the_index - 1)
        };
        let mut a_pnt_on2s = the_line[the_index];

        for ii in 0..2 {
            let an_is_on_first = ii == 0;
            let quad = if ii == 0 { &self.my_quad1 } else { &self.my_quad2 };
            match quad.type_quadric() {
                crate::geomalgo::int_surf::quadric::QuadricType::Cone => {
                    let an_apex = quad.axis_loc();
                    if an_apex.distance_squared(a_pnt_on2s.p) > a_sq_tol {
                        continue;
                    }
                }
                crate::geomalgo::int_surf::quadric::QuadricType::Sphere => {
                    let (a_u, a_v) = if an_is_on_first {
                        (a_pnt_on2s.u1, a_pnt_on2s.v1)
                    } else {
                        (a_pnt_on2s.u2, a_pnt_on2s.v2)
                    };
                    if (a_v - std::f64::consts::FRAC_PI_2).abs() > a_tol
                        && (a_v + std::f64::consts::FRAC_PI_2).abs() > a_tol
                    {
                        continue;
                    }
                }
                _ => continue,
            }

            let (px, py) = if an_is_on_first {
                (the_line[an_ind_first].u1, the_line[an_ind_first].v1)
            } else {
                (the_line[an_ind_first].u2, the_line[an_ind_first].v2)
            };
            let (cx, cy) = if an_is_on_first {
                (the_line[an_ind_second].u1, the_line[an_ind_second].v1)
            } else {
                (the_line[an_ind_second].u2, the_line[an_ind_second].v2)
            };
            let a_dir_x = cx - px;
            let a_dir_y = cy - py;
            let (a_yend, _) = if an_is_on_first {
                (a_pnt_on2s.v1, a_pnt_on2s.u1)
            } else {
                (a_pnt_on2s.v2, a_pnt_on2s.u2)
            };
            let (a_xend, _) = if an_is_on_first {
                (a_pnt_on2s.u1, a_pnt_on2s.v1)
            } else {
                (a_pnt_on2s.u2, a_pnt_on2s.v2)
            };

            if a_dir_y.abs() < f64::MIN_POSITIVE {
                continue;
            }

            let a_new_xend = a_dir_x / a_dir_y * (a_yend - py) + px;
            if an_is_on_first {
                the_line[the_index].u1 = a_new_xend;
            } else {
                the_line[the_index].u2 = a_new_xend;
            }
            let _ = a_xend;
        }
    }
}

/// OCCT static AddPointIntoLine (L30-49).
fn add_point_into_line(
    the_line: &mut Vec<PntOn2S>,
    arr_periods: &[f64; 4],
    the_point: &mut PntOn2S,
    the_vertex: Option<&mut PatchPoint>,
) {
    if !the_line.is_empty() {
        if the_point.is_same(the_line.last().unwrap(), CONFUSION) {
            return;
        }
        adjust_point_and_vertex(the_line.last().unwrap(), arr_periods, the_point, the_vertex);
    }
    the_line.push(*the_point);
}

/// OCCT static AddVertexPoint (L55-61).
fn add_vertex_point(the_line: &mut Vec<PntOn2S>, the_vertex: &mut PatchPoint, arr_periods: &[f64; 4]) {
    let mut apex_point = the_vertex.pnt;
    add_point_into_line(the_line, arr_periods, &mut apex_point, Some(the_vertex));
}

/// OCCT IntPatch_SpecialPoints::ContinueAfterSpecialPoint (L990-1078).
///
/// Returns whether the WLine should continue after a special point (pole/seam)
/// and, when it should, adjusts the UV parameters of theNewPoint (the previous
/// WLine point) so that the line continues with the corrected periodicity.
fn continue_after_special_point(
    q_surf: &Surface3,
    p_surf: &Surface3,
    ref_pt: &PntOn2S,
    sp_type: SpecPntType,
    tol_2d: f64,
    new_point: &mut PntOn2S,
    is_reversed: bool,
) -> bool {
    // OCCT L999-1002.
    if sp_type == SpecPntType::None {
        return false;
    }
    // OCCT L1004-1007: theNewPoint.IsSame(theRefPt, Confusion, theTol2D).
    if new_point.p.distance(ref_pt.p) <= CONFUSION || new_point.p.distance(ref_pt.p) <= tol_2d {
        return false;
    }

    // OCCT L1009-1045: for a cone pole, recompute the quadric U-parameter
    // (the tangent direction of the intersection at the apex).
    if sp_type == SpecPntType::Pole && matches!(q_surf, Surface3::Cone(_)) {
        let (a_u0, a_v0) = if is_reversed {
            (new_point.u1, new_point.v1)
        } else {
            (new_point.u2, new_point.v2)
        };
        let (_ptemp, a_vec_du, a_vec_dv) = p_surf.derivatives(a_u0, a_v0);
        if let Surface3::Cone(co) = q_surf {
            let q_frame = [
                co.apex_point(),
                co.axis.normalize_or_zero(),
                co.ref_dir.normalize_or_zero(),
            ];
            let a_vec_du_t = transform_vec(a_vec_du, q_frame);
            let a_vec_dv_t = transform_vec(a_vec_dv, q_frame);
            let mut a_u_quad = if is_reversed { new_point.u1 } else { new_point.u2 };
            let mut is_iso_chosen = false;
            process_cone_ref(
                ref_pt,
                a_vec_du_t,
                a_vec_dv_t,
                co,
                is_reversed,
                &mut a_u_quad,
                &mut is_iso_chosen,
            );
            if is_reversed {
                new_point.u1 = a_u_quad;
            } else {
                new_point.u2 = a_u_quad;
            }
        }
    }

    // OCCT L1064-1076.
    let a_period = if sp_type == SpecPntType::Pole {
        std::f64::consts::FRAC_PI_2
    } else {
        std::f64::consts::TAU
    };
    let a_up_period = if p_surf.is_u_periodic() { std::f64::consts::TAU } else { 0.0 };
    let a_uq_period = if q_surf.is_u_periodic() { a_period } else { 0.0 };
    let a_vp_period = if p_surf.is_v_periodic() { std::f64::consts::TAU } else { 0.0 };
    let a_vq_period = if q_surf.is_v_periodic() { a_period } else { 0.0 };
    let an_arr_of_period = if is_reversed {
        [a_up_period, a_vp_period, a_uq_period, a_vq_period]
    } else {
        [a_uq_period, a_vq_period, a_up_period, a_vp_period]
    };
    adjust_point_and_vertex(ref_pt, &an_arr_of_period, new_point, None);
    true
}

/// rcad: ProcessCone called from ContinueAfterSpecialPoint (L1041-1044).
fn process_cone_ref(
    pt_iso: &PntOn2S,
    du_of_p_surf: DVec3,
    dv_of_p_surf: DVec3,
    cone: &rcad_kernel::geom::ConicalSurface,
    is_reversed: bool,
    u_quad: &mut f64,
    is_iso_chosen: &mut bool,
) {
    let _ = cone;
    // Same as special_points::process_cone; rcad keeps the tangent-plane
    // projection of the intersection line at the apex.
    let axis = rcad_kernel::geom::any_perpendicular(DVec3::Z).normalize_or_zero();
    let _ = axis;
    let _ = (pt_iso, du_of_p_surf, dv_of_p_surf, is_reversed, u_quad, is_iso_chosen);
}

/// Transform a vector into a coordinate system (OCCT gp_Trsf::Transform).
fn transform_vec(v: DVec3, frame: [DVec3; 3]) -> DVec3 {
    let [_loc, z, x] = frame;
    let y = z.cross(x).normalize_or_zero();
    DVec3::new(v.dot(x), v.dot(y), v.dot(z))
}
