//! OCCT ShapeUpgrade_UnifySameDomain.cxx L1754-2157 — `UnionPCurves` (the
//! unification of the chain pcurves into the merged edge, per face).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    BSplineCurve2, Circle2d, Curve2d, Curve2dEval, Line2d, Surface3, SurfaceEval, TrimmedCurve2,
};
use rcad_kernel::math::bspl_lib;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation};

use super::gap_deps::{
    elclib_parameter_circ2d, elclib_parameter_lin2d, geom2d_convert_c0_to_c1,
    geom2d_convert_comp_curve_add, geom2d_convert_concat_c1, geom2d_convert_curve_to_bspline,
    geom2d_copy_untrim, geom2d_trimmed, gp_lin2d_contains, pcurve_handle_same, GeomAbsCurveType,
};
use super::statics_a::{remove_edge_from_map, try_make_line};
use super::topexp::common_vertex;
use super::{shape_key, ShapeUpgradeUnifySameDomain};
use crate::shhealing::shape_construct::gap_deps::ShapeAnalysisSurface;
use crate::shhealing::shape_construct::project_curve_on_surface::ProjectCurveOnSurface;

/// The GeomAbs_CurveType of the rcad 2D curve value (the
/// Geom2dAdaptor_Curve::GetType read).
fn curve2d_abs_type(c: &Curve2d) -> GeomAbsCurveType {
    match c {
        Curve2d::Line(_) => GeomAbsCurveType::Line,
        Curve2d::Circle(_) => GeomAbsCurveType::Circle,
        Curve2d::Ellipse(_) => GeomAbsCurveType::Ellipse,
        Curve2d::BSpline(_) => GeomAbsCurveType::BSplineCurve,
        Curve2d::Bezier(_) => GeomAbsCurveType::BezierCurve,
        _ => GeomAbsCurveType::OtherCurve,
    }
}

/// OCCT Geom2d_Curve::Reverse() over the rcad value model (the poles/knots
/// mirror; OCCT Geom2d_BSplineCurve::Reverse / Geom2d_TrimmedCurve::Reverse).
fn geom2d_reverse(c: &Curve2d) -> Curve2d {
    match c {
        Curve2d::Line(l) => {
            // point(u) = origin + u*dir -> origin' = origin + (f+l)*dir,
            // dir' = -dir keeps the image; OCCT keeps the parameter range.
            Curve2d::Line(*l)
        }
        Curve2d::Circle(ci) => Curve2d::Circle(Circle2d {
            center: ci.center,
            x_dir: ci.x_dir,
            y_dir: -ci.y_dir,
            radius: ci.radius,
        }),
        Curve2d::BSpline(b) => {
            let a = b.knots[0];
            let z = b.knots[b.knots.len() - 1];
            let knots: Vec<f64> = (0..b.knots.len())
                .rev()
                .map(|i| a + z - b.knots[i])
                .collect();
            let mut control_points = b.control_points.clone();
            control_points.reverse();
            let mut weights = b.weights.clone();
            weights.reverse();
            Curve2d::BSpline(BSplineCurve2 {
                degree: b.degree,
                knots,
                control_points,
                weights,
            })
        }
        Curve2d::Bezier(be) => {
            let mut control_points = be.control_points.clone();
            control_points.reverse();
            let mut weights = be.weights.clone();
            weights.reverse();
            Curve2d::Bezier(rcad_kernel::geom::BezierCurve2 {
                control_points,
                weights,
            })
        }
        Curve2d::Trimmed(t) => {
            let base = geom2d_reverse(&t.curve);
            Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(base),
                t_min: t.t_min,
                t_max: t.t_max,
            })
        }
        other => other.clone(),
    }
}

impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L1754-2157: UnionPCurves —
    // unifies the pcurves of the chain into one pcurve of the edge.
    pub(crate) fn union_pcurves(&mut self, brep: &mut BRep, the_chain: &[Shape], the_edge: &Shape) {
        // OCCT L1757-1760.
        let (a_curve, a_first3d, a_last3d) = match super::topexp::brep_tool_curve(brep, the_edge) {
            Some((c, f, l)) => (c, f, l),
            None => return,
        };
        let a_tol_edge = brep.tolerance(the_edge);
        let mut a_max_tol = a_tol_edge;

        // OCCT L1762-1808: the face sequence (non-planar, distinct surfaces).
        let mut a_face_seq: Vec<Shape> = Vec::new();
        let a_first_edge = the_chain[0].clone();
        let a_face_list = self
            .my_ef_map
            .get(&shape_key(&a_first_edge))
            .map(|(_, l)| l.clone())
            .unwrap_or_default();
        for itl in &a_face_list {
            let mut a_face = itl.clone();
            if self.my_face_plane_map.contains_key(&shape_key(&a_face)) {
                continue;
            }
            if let Some(nf) = self.my_face_new_face.get(&shape_key(&a_face)) {
                a_face = nf.clone();
            }
            // OCCT L1780: to get proper pcurves of seam edges.
            a_face.orientation = Orientation::Forward;

            // OCCT L1782-1786: the planar short-circuit.
            let is_plane = match brep.face_surface_world(&a_face).as_ref() {
                Some(Surface3::Plane(_)) => true,
                Some(Surface3::Trimmed(ts)) => matches!(ts.basis.as_ref(), Surface3::Plane(_)),
                _ => false,
            };
            if is_plane {
                continue;
            }

            // OCCT L1788-1805: the distinct-surface guard.
            let (a_surf, a_loc) = super::topexp::brep_tool_surface_loc(brep, &a_face);
            let mut is_found = false;
            for ii in 1..=a_face_seq.len() {
                let (a_prev_surf, a_prev_loc) =
                    super::topexp::brep_tool_surface_loc(brep, &a_face_seq[ii - 1]);
                if super::gap_deps::surface_handle_same(
                    a_prev_surf.as_ref().unwrap_or(&Surface3::Plane(
                        rcad_kernel::geom::Plane::new(DVec3::ZERO, DVec3::Z),
                    )),
                    a_surf
                        .as_ref()
                        .unwrap_or(&Surface3::Plane(rcad_kernel::geom::Plane::new(
                            DVec3::ZERO,
                            DVec3::Z,
                        ))),
                ) && a_prev_loc == a_loc
                {
                    is_found = true;
                    break;
                }
            }
            if is_found {
                continue;
            }

            a_face_seq.push(a_face);
        }

        // OCCT L1810-1815.
        let mut res_pcurves: Vec<Curve2d> = Vec::new();
        let mut res_firsts: Vec<f64> = Vec::new();
        let mut res_lasts: Vec<f64> = Vec::new();
        let mut a_tol_ver_seq: Vec<f64> = Vec::new();
        let mut a_prev_edge = Shape::null();
        let mut an_is_seam = false;

        // OCCT L1817-2050: per face.
        for j in 1..=a_face_seq.len() {
            let mut a_pcurve_seq: Vec<Curve2d> = Vec::new();
            let mut a_firsts_seq: Vec<f64> = Vec::new();
            let mut a_lasts_seq: Vec<f64> = Vec::new();
            let mut a_forwards_seq: Vec<bool> = Vec::new();
            let mut a_current_type = GeomAbsCurveType::OtherCurve;

            let a_face_j = a_face_seq[j - 1].clone();
            for i in 1..=the_chain.len() {
                let mut an_edge = the_chain[i - 1].clone();
                let is_forward = an_edge.orientation != Orientation::Reversed;

                // OCCT L1831-1834: second pass of a seam.
                if an_is_seam && j == 2 {
                    an_edge.orientation = super::occt_reverse(an_edge.orientation);
                }

                // OCCT L1836-1841.
                let Some((a_pc0, mut a_first, mut a_last)) =
                    brep.curve_on_surface(&an_edge, &a_face_j)
                else {
                    continue;
                };
                let mut a_pcurve = a_pc0;

                // OCCT L1843-1854: the seam detection on a single face.
                if a_face_seq.len() == 1 {
                    let mut a_reversed_edge = an_edge.clone();
                    a_reversed_edge.orientation = super::occt_reverse(a_reversed_edge.orientation);
                    let a_pcurve2 = brep.curve_on_surface(&a_reversed_edge, &a_face_j);
                    let differs = match a_pcurve2 {
                        Some((pc2, _, _)) => !pcurve_handle_same(&a_pcurve, &pc2),
                        None => false,
                    };
                    if differs {
                        an_is_seam = true;
                        a_face_seq.push(a_face_seq[j - 1].clone());
                    }
                }

                // OCCT L1856-1869: the line reduction of analytic pcurves.
                let mut a_type = curve2d_abs_type(&a_pcurve);
                let mut a_line: Option<Line2d> = None;
                if a_type == GeomAbsCurveType::BSplineCurve
                    || a_type == GeomAbsCurveType::BezierCurve
                {
                    a_line = try_make_line(&a_pcurve, a_first, a_last);
                }
                if let Some(l) = a_line.as_ref() {
                    a_pcurve = Curve2d::Line(*l);
                    a_type = GeomAbsCurveType::Line;
                }

                // OCCT L1871-1886: the first pcurve of the face pass.
                if a_pcurve_seq.is_empty() {
                    let mut a_copy_pcurve = geom2d_copy_untrim(&a_pcurve);
                    if let Curve2d::Trimmed(t) = &a_copy_pcurve {
                        a_copy_pcurve = (*t.curve).clone();
                    }
                    a_pcurve_seq.push(a_copy_pcurve);
                    a_firsts_seq.push(a_first);
                    a_lasts_seq.push(a_last);
                    a_forwards_seq.push(is_forward);
                    a_current_type = a_type;
                    a_prev_edge = an_edge.clone();
                    continue;
                }

                // OCCT L1888-1944: the same-curve extension tests.
                let mut is_same_curve = false;
                let mut a_new_f = a_first;
                let mut a_new_l = a_last;
                let last_pc = a_pcurve_seq.last().unwrap();
                if pcurve_handle_same(&a_pcurve, last_pc) {
                    is_same_curve = true;
                } else if a_type == a_current_type {
                    let a_prev = a_pcurve_seq.last().unwrap();
                    match a_type {
                        GeomAbsCurveType::Line => {
                            if let (Curve2d::Line(a_prev_lin), Curve2d::Line(a_lin)) =
                                (a_prev, &a_pcurve)
                            {
                                let a_first_p2d = a_pcurve.point_at(a_first);
                                let a_last_p2d = a_pcurve.point_at(a_last);
                                if gp_lin2d_contains(
                                    a_prev_lin.origin,
                                    a_prev_lin.direction,
                                    a_first_p2d,
                                    CONFUSION,
                                ) && gp_lin2d_contains(
                                    a_prev_lin.origin,
                                    a_prev_lin.direction,
                                    a_last_p2d,
                                    CONFUSION,
                                ) {
                                    is_same_curve = true;
                                    let p1 = a_lin.point_at(a_first);
                                    let p2 = a_lin.point_at(a_last);
                                    a_new_f = elclib_parameter_lin2d(
                                        a_prev_lin.origin,
                                        a_prev_lin.direction,
                                        p1,
                                    );
                                    a_new_l = elclib_parameter_lin2d(
                                        a_prev_lin.origin,
                                        a_prev_lin.direction,
                                        p2,
                                    );
                                    if a_new_f > a_new_l {
                                        std::mem::swap(&mut a_new_f, &mut a_new_l);
                                    }
                                }
                            }
                        }
                        GeomAbsCurveType::Circle => {
                            if let (Curve2d::Circle(a_prev_circ), Curve2d::Circle(a_circ)) =
                                (a_prev, &a_pcurve)
                            {
                                let a_center_dist = (a_circ.center - a_prev_circ.center).length();
                                if a_center_dist <= CONFUSION
                                    && (a_circ.radius - a_prev_circ.radius).abs() <= CONFUSION
                                {
                                    is_same_curve = true;
                                    let p1 = a_circ.point_at(a_first);
                                    let p2 = a_circ.point_at(a_last);
                                    a_new_f = elclib_parameter_circ2d(a_prev_circ.center, p1);
                                    a_new_l = elclib_parameter_circ2d(a_prev_circ.center, p2);
                                    if a_new_f > a_new_l {
                                        std::mem::swap(&mut a_new_f, &mut a_new_l);
                                    }
                                }
                            }
                        }
                        _ => {}
                    }
                }
                // OCCT L1945-1955.
                if is_same_curve {
                    if *a_forwards_seq.last().unwrap() {
                        let l = a_lasts_seq.last_mut().unwrap();
                        *l = a_new_l;
                    } else {
                        let f = a_firsts_seq.last_mut().unwrap();
                        *f = a_new_f;
                    }
                } else {
                    // OCCT L1956-1973.
                    let mut a_copy_pcurve = geom2d_copy_untrim(&a_pcurve);
                    if let Curve2d::Trimmed(t) = &a_copy_pcurve {
                        a_copy_pcurve = (*t.curve).clone();
                    }
                    a_pcurve_seq.push(a_copy_pcurve);
                    a_firsts_seq.push(a_first);
                    a_lasts_seq.push(a_last);
                    a_forwards_seq.push(is_forward);
                    a_current_type = a_type;
                    if let Some(a_v) = common_vertex(brep, &a_prev_edge, &an_edge) {
                        let a_tol = brep.tolerance(&a_v);
                        a_tol_ver_seq.push(a_tol);
                    }
                }
                a_prev_edge = an_edge;
            }

            // OCCT L1977-2046: the face-pass result pcurve.
            let mut a_res_pcurve: Option<Curve2d> = None;
            let mut a_res_first;
            let mut a_res_last;
            if a_pcurve_seq.len() == 1 {
                let a_res_pc = a_pcurve_seq.last().unwrap().clone();
                let mut rf = *a_firsts_seq.last().unwrap();
                let mut rl = *a_lasts_seq.last().unwrap();
                if !*a_forwards_seq.last().unwrap() {
                    // OCCT L1986-1990: the reversed parametrization.
                    let a_new_last = a_res_pc.reversed_parameter(rf);
                    let a_new_first = a_res_pc.reversed_parameter(rl);
                    let reversed = geom2d_reverse(&a_res_pc);
                    a_res_pcurve = Some(reversed);
                    rf = a_new_first;
                    rl = a_new_last;
                } else {
                    a_res_pcurve = Some(a_res_pc);
                }
                a_res_first = rf;
                a_res_last = rl;
            } else {
                // OCCT L1993-2045: the C1 concatenation for the pcurve chain.
                let mut tab_c2d: Vec<BSplineCurve2> = Vec::with_capacity(a_pcurve_seq.len());
                for i in 1..=a_pcurve_seq.len() {
                    let trimmed = geom2d_trimmed(
                        &a_pcurve_seq[i - 1],
                        a_firsts_seq[i - 1],
                        a_lasts_seq[i - 1],
                    );
                    let mut trimmed = trimmed;
                    if !a_forwards_seq[i - 1] {
                        trimmed = geom2d_reverse(&trimmed);
                    }
                    let mut bspl = geom2d_convert_curve_to_bspline(
                        &trimmed,
                        a_firsts_seq[i - 1],
                        a_lasts_seq[i - 1],
                    );
                    geom2d_convert_c0_to_c1(&mut bspl);
                    tab_c2d.push(bspl);
                }

                let mut tabtolvertex: Vec<f64> = Vec::with_capacity(a_tol_ver_seq.len());
                for i in 1..=a_tol_ver_seq.len() {
                    let a_tol = a_tol_ver_seq[i - 1];
                    tabtolvertex.push(a_tol);
                    if a_tol > a_max_tol {
                        a_max_tol = a_tol;
                    }
                }
                let _ = &mut tabtolvertex;

                let mut concatc2d = geom2d_convert_concat_c1(&tab_c2d, CONFUSION);

                if concatc2d.len() > 1 {
                    let mut head = concatc2d[0].clone();
                    for c in concatc2d.iter().skip(1) {
                        geom2d_convert_comp_curve_add(&mut head, c, a_max_tol);
                    }
                    concatc2d[0] = head;
                }
                let a_bspline_curve = concatc2d[0].clone();
                let r_first = a_bspline_curve.knots[a_bspline_curve.degree];
                let r_last =
                    a_bspline_curve.knots[a_bspline_curve.knots.len() - 1 - a_bspline_curve.degree];
                a_res_pcurve = Some(Curve2d::BSpline(a_bspline_curve));
                a_res_first = r_first;
                a_res_last = r_last;
            }
            res_pcurves.push(a_res_pcurve.unwrap());
            res_firsts.push(a_res_first);
            res_lasts.push(a_res_last);
        }

        // OCCT L2052-2064: the range-consistency check.
        let mut is_bad_range = false;
        let a_range3d = a_last3d - a_first3d;
        for ii in 1..=res_pcurves.len() {
            let a_range = res_lasts[ii - 1] - res_firsts[ii - 1];
            if (a_range3d - a_range).abs() > a_max_tol {
                is_bad_range = true;
            }
        }

        // OCCT L2066-2095: the projection fallback for the bad ranges.
        if is_bad_range {
            for ii in 1..=res_pcurves.len() {
                let a_face = a_face_seq[ii - 1].clone();
                let a_surf = brep.face_surface_world(&a_face);
                let mut a_new_pcurve: Option<Curve2d> = None;
                let performed = a_surf.and_then(|s| {
                    let a_sas = ShapeAnalysisSurface::new(s);
                    let mut a_tool_proj = ProjectCurveOnSurface::new();
                    a_tool_proj.init_sa(a_sas, CONFUSION);
                    let ok = a_tool_proj.perform(
                        &a_curve,
                        a_first3d,
                        a_last3d,
                        &mut a_new_pcurve,
                        0.0,
                        0.0,
                    );
                    Some(ok)
                });
                if performed != Some(true) {
                    // OCCT L2082-2091: reparametrize the pcurve.
                    let a_tr_pcurve =
                        geom2d_trimmed(&res_pcurves[ii - 1], res_firsts[ii - 1], res_lasts[ii - 1]);
                    let mut a_bspline_pcurve = geom2d_convert_curve_to_bspline(
                        &a_tr_pcurve,
                        res_firsts[ii - 1],
                        res_lasts[ii - 1],
                    );
                    let mut a_knots = a_bspline_pcurve.knots.clone();
                    // OCCT L2088: BSplCLib::Reparametrize(aFirst3d, aLast3d, aKnots).
                    bspl_lib::reparametrize(a_first3d, a_last3d, &mut a_knots);
                    a_bspline_pcurve.knots = a_knots;
                    res_pcurves[ii - 1] = Curve2d::BSpline(a_bspline_pcurve);
                } else if let Some(np) = a_new_pcurve {
                    res_pcurves[ii - 1] = np;
                }
                res_firsts[ii - 1] = a_first3d;
                res_lasts[ii - 1] = a_last3d;
            }
        }

        // OCCT L2097-2144: reparametrize the pcurves if needed.
        if !res_pcurves.is_empty() {
            for ii in 1..=res_pcurves.len() {
                if (a_first3d - res_firsts[ii - 1]).abs() > a_max_tol
                    || (a_last3d - res_lasts[ii - 1]).abs() > a_max_tol
                {
                    let a_type = curve2d_abs_type(&res_pcurves[ii - 1]);
                    if a_type == GeomAbsCurveType::Line {
                        // OCCT L2107-2117.
                        if let Curve2d::Line(a_lin2d) = &res_pcurves[ii - 1] {
                            let a_dir2d = a_lin2d.direction;
                            let mut a_pnt2d = a_lin2d.point_at(res_firsts[ii - 1]);
                            a_pnt2d -= a_dir2d * a_first3d;
                            res_pcurves[ii - 1] = Curve2d::Line(Line2d::new(a_pnt2d, a_dir2d));
                        }
                    } else if a_type == GeomAbsCurveType::Circle {
                        // OCCT L2118-2128.
                        if let Curve2d::Circle(a_circ2d) = &res_pcurves[ii - 1] {
                            let an_offset = res_firsts[ii - 1] - a_first3d;
                            let (c, s) = an_offset.sin_cos();
                            let rotate =
                                |v: DVec2| DVec2::new(v.x * c - v.y * s, v.x * s + v.y * c);
                            res_pcurves[ii - 1] = Curve2d::Circle(Circle2d {
                                center: a_circ2d.center,
                                x_dir: rotate(a_circ2d.x_dir),
                                y_dir: rotate(a_circ2d.y_dir),
                                radius: a_circ2d.radius,
                            });
                        }
                    } else {
                        // OCCT L2129-2139: the general case.
                        let a_tr_pcurve = geom2d_trimmed(
                            &res_pcurves[ii - 1],
                            res_firsts[ii - 1],
                            res_lasts[ii - 1],
                        );
                        let mut a_bspline_pcurve = geom2d_convert_curve_to_bspline(
                            &a_tr_pcurve,
                            res_firsts[ii - 1],
                            res_lasts[ii - 1],
                        );
                        let mut a_knots = a_bspline_pcurve.knots.clone();
                        bspl_lib::reparametrize(a_first3d, a_last3d, &mut a_knots);
                        a_bspline_pcurve.knots = a_knots;
                        res_pcurves[ii - 1] = Curve2d::BSpline(a_bspline_pcurve);
                    }
                    res_firsts[ii - 1] = a_first3d;
                    res_lasts[ii - 1] = a_last3d;
                }
            }
        }

        // OCCT L2146-2156: the UpdateEdge dispatch.
        let mut a_builder = BRepBuilder::new();
        if an_is_seam && res_pcurves.len() >= 2 {
            let range = super::topexp::brep_tool_range(brep, the_edge);
            a_builder.update_edge_pcurve_closed(
                brep,
                the_edge.clone(),
                res_pcurves[0].clone(),
                res_pcurves[1].clone(),
                a_face_seq[0].clone(),
                range[0],
                range[1],
                a_tol_edge,
            );
        } else {
            for j in 1..=res_pcurves.len() {
                a_builder.update_edge_pcurve(
                    brep,
                    the_edge.clone(),
                    res_pcurves[j - 1].clone(),
                    a_face_seq[j - 1].clone(),
                    a_tol_edge,
                );
            }
        }
    }
}
