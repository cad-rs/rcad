// OCCT IntPatch_ImpImpIntersection restriction post-processing:
//   PutPointsOnLine (L439-724), MultiplePoint (L726-916),
//   PointOnSecondDom (L918-1079), FindLine (L1081-1392),
//   SingleLine (L1402-1484), ProcessSegments (L1486-1885),
//   IsRLineGood (L1933-2027), ProcessRLine (L2029-2343).
//
// 1:1 Rust translation.
//
// rcad data-model notes:
//   - NCollection_Sequence<handle<IntPatch_Line>> slin -> Vec<IntPatchLine>.
//   - IntPatch_GLine / IntPatch_ALine / IntPatch_RLine vertices all map to
//     IntPatchLine.vertices (GLine/RLine) / IntPatchLine.a_curve.vertices
//     (ALine).  The AddVertex/Replace/NbVertex/Vertex ops are provided by the
//     line methods in int_patch::mod.rs.
//   - The 2D boundary arc (Adaptor2d_Curve2d) is a Curve2d; the path-point
//     sequence is Vec<PathPoint> from so_on_bounds.rs.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3};

use super::so_on_bounds::{Domain, PathPoint};
use super::transitions::{make_transition, recadre, Transition, TypeTrans};
use super::{IntPatchIType, IntPatchLine, IntPatchVertex};
use crate::geomalgo::int_surf::quadric::Quadric;

/// OCCT gp_Vec::Angle(gp_Vec) — the angular value between two vectors in
/// [0, PI], computed exactly as gp_Dir::Angle (gp_Dir.cxx L27-54).
fn gp_vec_angle(a: DVec3, b: DVec3) -> f64 {
    let a = a.normalize_or_zero();
    let b = b.normalize_or_zero();
    let cosinus = a.dot(b);
    if cosinus > -0.70710678118655 && cosinus < 0.70710678118655 {
        cosinus.acos()
    } else {
        let sinus = a.cross(b).length();
        if cosinus < 0.0 {
            std::f64::consts::PI - sinus.asin()
        } else {
            sinus.asin()
        }
    }
}

/// OCCT gp_Vec::IsParallel (gp_Vec.hxx L142-146) — angle-based definition
/// (opposite directions count as parallel).
fn gp_vec_is_parallel(a: DVec3, b: DVec3, ang_tol: f64) -> bool {
    let an_ang = gp_vec_angle(a, b);
    an_ang <= ang_tol || std::f64::consts::PI - an_ang <= ang_tol
}

/// OCCT IntPatch_IType of a line (ArcType).
fn arc_type_of(line: &IntPatchLine) -> IntPatchIType {
    line.line_type
}

/// OCCT GLine/ALine::AddVertex equivalent.
fn add_vertex(line: &mut IntPatchLine, v: IntPatchVertex) {
    if line.line_type == IntPatchIType::Analytic {
        if let Some(ac) = line.a_curve.as_mut() {
            let pp = vertex_to_patch_point(&v);
            ac.vertices.push(pp);
        }
    } else {
        line.vertices.push(v);
    }
}

/// OCCT GLine/ALine::Replace(Index, Pnt) equivalent — 1-based index.
fn replace_vertex(line: &mut IntPatchLine, index: usize, v: IntPatchVertex) {
    if line.line_type == IntPatchIType::Analytic {
        if let Some(ac) = line.a_curve.as_mut() {
            ac.vertices[index - 1] = vertex_to_patch_point(&v);
        }
    } else {
        line.vertices[index - 1] = v;
    }
}

/// OCCT GLine/ALine::NbVertex() equivalent.
fn nb_vertex(line: &IntPatchLine) -> usize {
    if line.line_type == IntPatchIType::Analytic {
        line.a_curve.as_ref().map(|c| c.vertices.len()).unwrap_or(0)
    } else {
        line.vertices.len()
    }
}

/// OCCT GLine/ALine::Vertex(Index) equivalent — 1-based.
fn vertex_of(line: &IntPatchLine, index: usize) -> IntPatchVertex {
    if line.line_type == IntPatchIType::Analytic {
        let pp = &line.a_curve.as_ref().unwrap().vertices[index - 1];
        patch_point_to_vertex(pp)
    } else {
        line.vertices[index - 1].clone()
    }
}

fn vertex_to_patch_point(v: &IntPatchVertex) -> super::special_points::PatchPoint {
    super::special_points::PatchPoint {
        pnt: super::special_points::PntOn2S {
            p: v.p3d,
            u1: v.u1,
            v1: v.v1,
            u2: v.u2,
            v2: v.v2,
        },
        param_on_line: v.param_on_line,
        tolerance: v.tolerance,
        multiple: v.multiple,
        on_dom_s1: v.on_dom_s1,
        on_dom_s2: v.on_dom_s2,
        arc_on_s1: v.arc_on_s1.clone(),
        arc_on_s2: v.arc_on_s2.clone(),
        param_on_arc1: v.param_on_arc1,
        param_on_arc2: v.param_on_arc2,
        is_vertex_on_s1: v.is_vertex_on_s1,
        is_vertex_on_s2: v.is_vertex_on_s2,
        transition_line_arc1: v.transition_line_arc1,
        transition_line_arc2: v.transition_line_arc2,
        transition_on_s1: v.transition_on_s1,
        transition_on_s2: v.transition_on_s2,
    }
}

fn patch_point_to_vertex(pp: &super::special_points::PatchPoint) -> IntPatchVertex {
    IntPatchVertex {
        param_on_line: pp.param_on_line,
        p3d: pp.pnt.p,
        u1: pp.pnt.u1,
        v1: pp.pnt.v1,
        u2: pp.pnt.u2,
        v2: pp.pnt.v2,
        tolerance: pp.tolerance,
        tangent: false,
        multiple: pp.multiple,
        on_dom_s1: pp.on_dom_s1,
        on_dom_s2: pp.on_dom_s2,
        arc_on_s1: pp.arc_on_s1.clone(),
        arc_on_s2: pp.arc_on_s2.clone(),
        param_on_arc1: pp.param_on_arc1,
        param_on_arc2: pp.param_on_arc2,
        is_vertex_on_s1: pp.is_vertex_on_s1,
        is_vertex_on_s2: pp.is_vertex_on_s2,
        transition_line_arc1: pp.transition_line_arc1,
        transition_line_arc2: pp.transition_line_arc2,
        transition_on_s1: pp.transition_on_s1,
        transition_on_s2: pp.transition_on_s2,
    }
}

/// OCCT ElCLib::Parameter(gp_Lin, gp_Pnt) — parameter of a point on a 3D line.
fn elclib_line_parameter(line: &rcad_kernel::geom::Line3, p: DVec3) -> f64 {
    (p - line.origin).dot(line.direction.normalize_or_zero())
}

/// OCCT ElCLib::Value(para, gp_Lin).
fn elclib_line_value(line: &rcad_kernel::geom::Line3, t: f64) -> DVec3 {
    line.point_at(t)
}

/// OCCT ElCLib::Parameter(gp_Circ, gp_Pnt) — angle parameter of a point on a
/// circle.  OCCT CircleParameter (ElCLib.cxx L1199-1222) computes
/// Teta = XDirection.AngleWithRef(aVProj, dir), which returns the angle in
/// [0, 2*PI) (normalizeAngle), so a negative atan2 angle is shifted by 2*PI.
fn elclib_circle_parameter(circ: &rcad_kernel::geom::Circle3, p: DVec3) -> f64 {
    let d = p - circ.center;
    let x = d.dot(circ.x_dir.normalize_or_zero());
    let y = d.dot(circ.y_dir.normalize_or_zero());
    let mut t = y.atan2(x);
    if t < 0.0 {
        t += std::f64::consts::TAU;
    }
    t
}

/// OCCT ElCLib::Value(para, gp_Circ).
fn elclib_circle_value(circ: &rcad_kernel::geom::Circle3, t: f64) -> DVec3 {
    circ.point_at(t)
}

/// OCCT ElCLib::DN(para, gp_Circ, 1) — first derivative of the circle.
fn elclib_circle_d1(circ: &rcad_kernel::geom::Circle3, t: f64) -> DVec3 {
    let x = circ.x_dir.normalize_or_zero();
    let y = circ.y_dir.normalize_or_zero();
    circ.radius * (-t.sin() * x + t.cos() * y)
}

/// OCCT FindLine (L1081-1392).  Projects Psurf onto the intersection lines.
#[allow(clippy::too_many_arguments)]
fn find_line(
    psurf_in: &mut DVec3,
    slin: &[IntPatchLine],
    tol: f64,
    the_l_params: &mut Vec<f64>,
    v_tgt_int: &mut DVec3,
    the_line_idx: &mut usize,
    only_this_line: usize,
    the_arc: &Curve2d,
    the_parameter_on_arc: &mut f64,
    the_point_on_arc: &mut DVec3,
    quad_surf1: &Quadric,
    quad_surf2: &Quadric,
    the_output_toler: &mut f64,
) -> bool {
    let mut psurf = *psurf_in;
    if quad_surf1.distance(psurf) > tol || quad_surf2.distance(psurf) > tol {
        return false;
    }

    let a_sq_tol = tol * tol;
    let mut a_sq_dist_min = f64::MAX;
    let mut a_sq_dist;
    let mut para;
    let mut lower;
    let mut upper;
    let mut pt;
    let mut typarc;

    let mut a_para_int = f64::MAX;
    let mut nblin = slin.len();
    let mut i = 0;
    while i < nblin {
        if only_this_line != 0 {
            i = only_this_line - 1;
            nblin = 0;
        }
        let lin = &slin[i];
        typarc = arc_type_of(lin);
        if typarc == IntPatchIType::Analytic {
            let ac = lin.a_curve.as_ref().unwrap();
            lower = ac.domain()[0];
            upper = ac.domain()[1];
        } else {
            if lin.has_first_point() {
                lower = lin.first_point().parameter_on_line();
            } else {
                lower = f64::NEG_INFINITY;
            }
            if lin.has_last_point() {
                upper = lin.last_point().parameter_on_line();
            } else {
                upper = f64::INFINITY;
            }
        }

        match typarc {
            IntPatchIType::Line => {
                if let Curve3::Line(l) = &lin.curve {
                    para = elclib_line_parameter(l, psurf);
                    if para <= upper && para >= lower {
                        // OCCT: gp_Lin::Direction() is a unit vector, so
                        // ElCLib::Value(para, Lin) = Loc + para * Dir with the
                        // SAME normalized direction ElCLib::Parameter used.  A
                        // raw Line3 may carry a non-unit direction, so the point
                        // must be rebuilt with the normalized direction, not
                        // line.point_at(t) (which uses the raw direction).
                        let dir_n = l.direction.normalize_or_zero();
                        pt = l.origin + dir_n * para;
                        a_sq_dist = psurf.distance_squared(pt);
                        if (a_sq_dist < a_sq_tol) && (a_sq_dist < a_sq_dist_min) {
                            a_sq_dist_min = a_sq_dist;
                            a_para_int = para;
                            *the_line_idx = i;
                        }
                    }
                }
            }
            IntPatchIType::Circle => {
                if let Curve3::Circle(c) = &lin.curve {
                    para = elclib_circle_parameter(c, psurf);
                    if (para <= upper && para >= lower)
                        || (para + std::f64::consts::TAU <= upper
                            && para + std::f64::consts::TAU >= lower)
                        || (para - std::f64::consts::TAU <= upper
                            && para - std::f64::consts::TAU >= lower)
                    {
                        pt = elclib_circle_value(c, para);
                        a_sq_dist = psurf.distance_squared(pt);
                        if (a_sq_dist < a_sq_tol) && (a_sq_dist < a_sq_dist_min) {
                            a_sq_dist_min = a_sq_dist;
                            a_para_int = para;
                            *the_line_idx = i;
                        }
                    }
                }
            }
            IntPatchIType::Ellipse => {
                if let Curve3::Ellipse(e) = &lin.curve {
                    para = elclib_ellipse_parameter(e, psurf);
                    if (para <= upper && para >= lower)
                        || (para + std::f64::consts::TAU <= upper
                            && para + std::f64::consts::TAU >= lower)
                        || (para - std::f64::consts::TAU <= upper
                            && para - std::f64::consts::TAU >= lower)
                    {
                        pt = e.point_at(para);
                        a_sq_dist = psurf.distance_squared(pt);
                        if (a_sq_dist < a_sq_tol) && (a_sq_dist < a_sq_dist_min) {
                            a_sq_dist_min = a_sq_dist;
                            a_para_int = para;
                            *the_line_idx = i;
                        }
                    }
                }
            }
            IntPatchIType::Parabola => {
                if let Curve3::Parabola(par) = &lin.curve {
                    para = elclib_parabola_parameter(par, psurf);
                    if para <= upper && para >= lower {
                        let mut amelioration = 0;
                        loop {
                            let parabis = para + 0.0000001;
                            pt = par.point_at(para);
                            a_sq_dist = psurf.distance_squared(pt);
                            let ptbis = par.point_at(parabis);
                            let distbis = psurf.distance(ptbis);
                            let a_dist = a_sq_dist.sqrt();
                            let ddist = distbis - a_dist;
                            if (a_sq_dist < a_sq_tol) && (a_sq_dist < a_sq_dist_min) {
                                a_sq_dist_min = a_sq_dist;
                                a_para_int = para;
                                *the_line_idx = i;
                            }
                            if a_sq_dist < rcad_kernel::precision::square_p_confusion() {
                                amelioration = 100;
                            }
                            if ddist > 1.0e-9 || ddist < -1.0e-9 {
                                para = para - a_dist * (parabis - para) / ddist;
                            } else {
                                amelioration = 100;
                            }
                            amelioration += 1;
                            if amelioration >= 5 {
                                break;
                            }
                        }
                    }
                }
            }
            IntPatchIType::Hyperbola => {
                if let Curve3::Hyperbola(h) = &lin.curve {
                    para = elclib_hyperbola_parameter(h, psurf);
                    if para <= upper && para >= lower {
                        pt = h.point_at(para);
                        a_sq_dist = psurf.distance_squared(pt);
                        if (a_sq_dist < a_sq_tol) && (a_sq_dist < a_sq_dist_min) {
                            a_sq_dist_min = a_sq_dist;
                            a_para_int = para;
                            *the_line_idx = i;
                        }
                    }
                }
            }
            IntPatchIType::Analytic => {
                let ac = lin.a_curve.as_ref().unwrap();
                let a_l_params = ac.find_parameter(psurf);
                if !a_l_params.is_empty() {
                    a_sq_dist = f64::MAX;
                    for &p in &a_l_params {
                        let pt = ac.value(p).unwrap_or(DVec3::ZERO);
                        let a_sq_d = psurf.distance_squared(pt);
                        if a_sq_d < a_sq_dist {
                            a_sq_dist = a_sq_d;
                        }
                    }
                    if a_sq_dist < a_sq_dist_min {
                        a_sq_dist_min = a_sq_dist;
                        *the_l_params = a_l_params;
                        *the_line_idx = i;
                    }
                } else {
                    // The point was not found by direct projection; try the
                    // 2D restriction intersection.
                    let mut copie_psurf = psurf;
                    let mut theparamonarc = *the_parameter_on_arc;
                    let mut theparam = 0.0;
                    let intersect_ok = intersection_with_an_arc(
                        &mut copie_psurf,
                        ac,
                        &mut theparam,
                        the_arc,
                        &mut theparamonarc,
                        the_point_on_arc,
                        quad_surf1,
                        lower,
                        upper,
                    );
                    a_sq_dist = copie_psurf.distance_squared(psurf);
                    if intersect_ok {
                        if a_sq_dist < a_sq_tol {
                            *the_parameter_on_arc = theparamonarc;
                            psurf = *the_point_on_arc;
                            a_sq_dist_min = a_sq_dist;
                            the_l_params.push(theparam);
                            *the_line_idx = i;
                        }
                    }
                }
            }
            _ => {}
        }
        i += 1;
    }

    if a_sq_dist_min == f64::MAX {
        return false;
    }

    *the_output_toler = (*the_output_toler).max(a_sq_dist_min.sqrt());

    typarc = arc_type_of(&slin[*the_line_idx]);

    // Computation of the tangent vector.
    match typarc {
        IntPatchIType::Line => {
            the_l_params.push(a_para_int);
            if let Curve3::Line(l) = &slin[*the_line_idx].curve {
                *v_tgt_int = l.direction;
            }
        }
        IntPatchIType::Circle => {
            the_l_params.push(a_para_int);
            if let Curve3::Circle(c) = &slin[*the_line_idx].curve {
                *v_tgt_int = elclib_circle_d1(c, a_para_int);
            }
        }
        IntPatchIType::Ellipse => {
            the_l_params.push(a_para_int);
            if let Curve3::Ellipse(e) = &slin[*the_line_idx].curve {
                *v_tgt_int = e.derivative_at(a_para_int);
            }
        }
        IntPatchIType::Parabola => {
            the_l_params.push(a_para_int);
            if let Curve3::Parabola(p) = &slin[*the_line_idx].curve {
                *v_tgt_int = p.derivative_at(a_para_int);
            }
        }
        IntPatchIType::Hyperbola => {
            the_l_params.push(a_para_int);
            if let Curve3::Hyperbola(h) = &slin[*the_line_idx].curve {
                *v_tgt_int = h.derivative_at(a_para_int);
            }
        }
        IntPatchIType::Analytic => {
            let ac = slin[*the_line_idx].a_curve.as_ref().unwrap();
            match ac.d1u(the_l_params[the_l_params.len() - 1]) {
                Some((_pt, tg)) => {
                    *v_tgt_int = tg;
                }
                None => {
                    *v_tgt_int = DVec3::ZERO;
                }
            }
        }
        _ => {}
    }
    true
}

/// OCCT ElCLib::Parameter(gp_Elips, gp_Pnt) — approximate angle parameter.
/// OCCT ElCLib::EllipseParameter (ElCLib.cxx L1226-1249): the eccentric
/// anomaly.  The vector Om = NX*X + NY*(MajorRadius/MinorRadius)*Y is built
/// by scaling the Y component by the radius ratio, then Teta =
/// AngleWithRef(X, Om, Direction) = atan2(NY*(Major/Minor), NX).  The radius
/// ratio scaling is essential — without it the parameter of a point on the
/// ellipse is wrong (only valid for a circle).  The angle is normalized to
/// [0, 2*PI) (OCCT normalizeAngle, ElCLib.cxx L56-68) exactly like the
/// circle parameter, so the SOnBounds vertex parameters stay in the curve
/// domain [0, 2*PI].
fn elclib_ellipse_parameter(e: &rcad_kernel::geom::Ellipse3, p: DVec3) -> f64 {
    let d = p - e.center;
    let x = d.dot(e.major_dir.normalize_or_zero());
    let y = d.dot(e.normal.cross(e.major_dir).normalize_or_zero());
    let mut t = (y * e.major_radius / e.minor_radius).atan2(x);
    if t < 0.0 {
        t += std::f64::consts::TAU;
    }
    t
}

/// OCCT ElCLib::ParabolaParameter (ElCLib.cxx L1269-1272):
///   `t = (P - vertex) . YDirection`, where Y = N x X is the cross-axis direction.
fn elclib_parabola_parameter(p: &rcad_kernel::geom::Parabola3, pt: DVec3) -> f64 {
    let dir_perp = p.normal.cross(p.axis_dir).normalize_or_zero();
    (pt - p.vertex).dot(dir_perp)
}

/// OCCT ElCLib::HyperbolaParameter (ElCLib.cxx L1253-1265):
///   `sht = (P - center) . YDirection / MinorRadius;  t = asinh(sht)`.
fn elclib_hyperbola_parameter(h: &rcad_kernel::geom::Hyperbola3, pt: DVec3) -> f64 {
    let minor_dir = h.normal.cross(h.major_dir).normalize_or_zero();
    let sht = (pt - h.center).dot(minor_dir) / h.semi_minor;
    sht.asinh()
}

/// OCCT IntersectionWithAnArc (L127-314): coarse parameter search over the
/// ALine domain followed by a Newton iteration that matches the ALine UV
/// against the boundary arc parameter.
#[allow(clippy::too_many_arguments)]
fn intersection_with_an_arc(
    psurf: &mut DVec3,
    alin: &super::int_quad_quad::IntAnaCurve,
    para: &mut f64,
    the_arc: &Curve2d,
    the_parameter_on_arc: &mut f64,
    the_point_on_arc: &mut DVec3,
    quad_surf: &Quadric,
    u0alin: f64,
    u1alin: f64,
) -> bool {
    // OCCT L142: dtheta = (u1alin - u0alin) * 0.01.
    let dtheta = (u1alin - u0alin) * 0.01;
    // OCCT L143-147: du = 1e-9; if (du >= dtheta) du = dtheta / 2.
    let mut du = 0.000000001;
    if du >= dtheta {
        du = dtheta / 2.0;
    }
    let mut distmin = f64::MAX;
    let mut thetamin = 0.0;
    let mut theparameteronarc = *the_parameter_on_arc;

    // OCCT L153-162: coarse search of the point of the ALine closest to PSurf.
    let mut _theta = u0alin + dtheta;
    while _theta <= u1alin - dtheta {
        let p = alin.value(_theta).unwrap_or(DVec3::ZERO);
        let d = p.distance(*psurf);
        if d < distmin {
            thetamin = _theta;
            distmin = d;
        }
        _theta += dtheta;
    }

    // OCCT L164-176: initial distance.
    let mut bestpara = 0.0;
    let mut besttheta = 0.0;
    let mut bestdist = 0.0;
    let mut distinit = 0.0;
    {
        let pp0 = alin.value(thetamin).unwrap_or(DVec3::ZERO);
        let (ua0, va0) = quad_surf.parameters(pp0);
        let p2d = the_arc.point_at(theparameteronarc);
        let pa_pr = DVec2::new(ua0 - p2d.x, va0 - p2d.y);
        distinit = pa_pr.length();
    }
    let mut theta = thetamin;
    // OCCT L179-182.
    let mut nbiter = 0;
    let drmax = (the_arc.default_domain()[1] - the_arc.default_domain()[0]) * 0.05;
    let damax = (u1alin - u0alin) * 0.05;
    bestdist = f64::MAX;

    loop {
        let pp0 = alin.value(theta).unwrap_or(DVec3::ZERO);
        let pp1 = alin.value(theta + du).unwrap_or(DVec3::ZERO);
        let (ua0, va0) = quad_surf.parameters(pp0);
        let (ua1, va1) = quad_surf.parameters(pp1);
        let d1a = DVec2::new((ua1 - ua0) / du, (va1 - va0) / du);
        let p2d = the_arc.point_at(theparameteronarc);
        let d2d = the_arc.derivative_at(theparameteronarc);
        let pa_pr = DVec2::new(ua0 - p2d.x, va0 - p2d.y);

        let pbd = pa_pr.length();
        if bestdist > pbd {
            bestdist = pbd;
            bestpara = theparameteronarc;
            besttheta = theta;
        }

        let d1a = DVec2::new(-d1a.x, -d1a.y);
        let d = d1a.x * d2d.y - d1a.y * d2d.x;
        let mut da = (-pa_pr.x) * d2d.y - (-pa_pr.y) * d2d.x;
        let mut dr = d1a.x * (-pa_pr.y) - d1a.y * (-pa_pr.x);
        if d.abs() > 1e-15 {
            da /= d;
            dr /= d;
        } else {
            // OCCT L223-248: fallback when the Jacobian is null.
            if pa_pr.x.abs() > pa_pr.y.abs() {
                let mut xx = pa_pr.x;
                xx *= 0.5;
                if d1a.x != 0.0 {
                    da = -xx / d1a.x;
                }
                if d2d.x != 0.0 {
                    dr = -xx / d2d.x;
                }
            } else {
                let mut yy = pa_pr.y;
                yy *= 0.5;
                if d1a.y != 0.0 {
                    da = -yy / d1a.y;
                }
                if d2d.y != 0.0 {
                    dr = -yy / d2d.y;
                }
            }
        }

        // OCCT L253-268: clamp the increments.
        if da < -damax {
            da = -damax;
        } else if da > damax {
            da = damax;
        }
        if dr < -drmax {
            dr = -drmax;
        } else if dr > drmax {
            dr = drmax;
        }

        // OCCT L270-279: converged.
        if da.abs() < 1e-10 && dr.abs() < 1e-10 {
            *para = theta;
            *psurf = alin.value(*para).unwrap_or(DVec3::ZERO);
            *the_parameter_on_arc = theparameteronarc;
            *the_point_on_arc = alin.value(*para).unwrap_or(DVec3::ZERO);
            return true;
        }
        // OCCT L282-299: step.
        theta += da;
        theparameteronarc += dr;
        let arc_dom = the_arc.default_domain();
        if theparameteronarc > arc_dom[1] {
            theparameteronarc = arc_dom[1];
        }
        if theparameteronarc < arc_dom[0] {
            theparameteronarc = arc_dom[0];
        }
        if theta < u0alin {
            theta = u0alin;
        }
        if theta > u1alin - du {
            theta = u1alin - du - du;
        }
        nbiter += 1;
        if nbiter >= 20 {
            break;
        }
    }

    // OCCT L303-311: fall back to the best sample found.
    if bestdist < distinit {
        *para = besttheta;
        *psurf = alin.value(*para).unwrap_or(DVec3::ZERO);
        *the_parameter_on_arc = bestpara;
        *the_point_on_arc = alin.value(*para).unwrap_or(DVec3::ZERO);
        return true;
    }
    false
}

/// OCCT SingleLine (L1402-1484).  Projects Psurf onto the single line.
fn single_line(
    psurf: DVec3,
    lin: &IntPatchLine,
    tol: f64,
    paraint: &mut f64,
    v_tgt_int: &mut DVec3,
) -> bool {
    let typarc = arc_type_of(lin);
    let mut parproj = 0.0;
    let mut tgint = DVec3::ZERO;
    let mut ptproj = DVec3::ZERO;
    let mut retvalue = false;

    match typarc {
        IntPatchIType::Line => {
            if let Curve3::Line(l) = &lin.curve {
                parproj = elclib_line_parameter(l, psurf);
                // OCCT gp_Lin::Direction() is a unit vector — the projected
                // point uses the same normalized direction as the parameter
                // (see find_line Line branch).
                ptproj = l.origin + l.direction.normalize_or_zero() * parproj;
                tgint = l.direction;
            }
        }
        IntPatchIType::Circle => {
            if let Curve3::Circle(c) = &lin.curve {
                parproj = elclib_circle_parameter(c, psurf);
                ptproj = elclib_circle_value(c, parproj);
                tgint = elclib_circle_d1(c, parproj);
            }
        }
        IntPatchIType::Ellipse => {
            if let Curve3::Ellipse(e) = &lin.curve {
                parproj = elclib_ellipse_parameter(e, psurf);
                ptproj = e.point_at(parproj);
                tgint = e.derivative_at(parproj);
            }
        }
        IntPatchIType::Parabola => {
            if let Curve3::Parabola(p) = &lin.curve {
                parproj = elclib_parabola_parameter(p, psurf);
                ptproj = p.point_at(parproj);
                tgint = p.derivative_at(parproj);
            }
        }
        IntPatchIType::Hyperbola => {
            if let Curve3::Hyperbola(h) = &lin.curve {
                parproj = elclib_hyperbola_parameter(h, psurf);
                ptproj = h.point_at(parproj);
                tgint = h.derivative_at(parproj);
            }
        }
        IntPatchIType::Analytic => {
            let ac = lin.a_curve.as_ref().unwrap();
            let a_l_params = ac.find_parameter(psurf);
            if !a_l_params.is_empty() {
                ptproj = psurf;
                parproj = a_l_params[a_l_params.len() - 1];
                match ac.d1u(parproj) {
                    Some((_pt, tg)) => {
                        tgint = tg;
                    }
                    None => {
                        tgint = DVec3::ZERO;
                    }
                }
            } else {
                return false;
            }
        }
        _ => {}
    }

    if psurf.distance(ptproj) <= tol {
        *paraint = parproj;
        *v_tgt_int = tgint;
        retvalue = true;
    } else {
        retvalue = false;
    }
    retvalue
}

/// OCCT MultiplePoint (L726-916).
#[allow(clippy::too_many_arguments)]
fn multiple_point(
    listpnt: &[PathPoint],
    domain: &Domain,
    quad_surf: &Quadric,
    normale: DVec3,
    slin: &mut [IntPatchLine],
    done: &mut [i32],
    used_line: &mut [i32],
    index: usize,
    on_first: bool,
    the_toler: f64,
) -> bool {
    let mut localdone = done.to_vec();
    let nblin = slin.len();
    let nbpnt = listpnt.len();
    let currentpointonrst = &listpnt[index - 1];
    let point = currentpointonrst.value();
    let mut retvalue = true;

    let mut ii = 0;
    while ii < nblin {
        let the_type = slin[ii].line_type;
        let nbvtx = nb_vertex(&slin[ii]);
        let mut jj = 1usize;
        while jj <= nbvtx {
            let mut intpt = vertex_of(&slin[ii], jj);
            if intpt.multiple
                && ((on_first && !intpt.on_dom_s1) || (!on_first && !intpt.on_dom_s2))
            {
                if point.distance(intpt.p3d) <= intpt.tolerance {
                    retvalue = false;
                    let mut paraint = 0.0;
                    let mut v_tgt_int = DVec3::ZERO;
                    if !single_line(point, &slin[ii].clone(), intpt.tolerance, &mut paraint, &mut v_tgt_int)
                    {
                        return false;
                    }
                    let goon;
                    // OCCT L808-814: SetVertex when the point is a domain vertex.
                    if !currentpointonrst.is_new() {
                        goon = true;
                        intpt.set_vertex(on_first);
                    } else {
                        goon = false;
                    }
                    let currentarc = currentpointonrst.arc().clone();
                    let currentparameter = currentpointonrst.parameter();
                    let p2d = currentarc.point_at(currentparameter);
                    let d2d = currentarc.derivative_at(currentparameter);
                    let (_ptbid, d1u, d1v) = quad_surf.d1(p2d.x, p2d.y);
                    let v_tgrst = d2d.x * d1u + d2d.y * d1v;

                    let mut transline = Transition::new();
                    let mut transarc = Transition::new();
                    if normale.length_squared() < 1e-16 {
                        transline.set_value_in_out(true, TypeTrans::Undecided);
                        transarc.set_value_in_out(true, TypeTrans::Undecided);
                    } else {
                        make_transition(v_tgt_int, v_tgrst, normale, &mut transline, &mut transarc);
                    }

                    intpt.p3d = point;
                    intpt.set_arc(on_first, currentarc, currentparameter, transline, transarc);
                    intpt.tolerance = the_toler;

                    replace_vertex(&mut slin[ii], jj, intpt.clone());
                    localdone[index - 1] = 1;
                    if goon {
                        for k in index..nbpnt {
                            if done[k] != 1 {
                                let otherpt = &listpnt[k];
                                if !otherpt.is_new() {
                                    let vtxbis = otherpt.vertex();
                                    if domain.identical(currentpointonrst.vertex(), vtxbis) {
                                        let oarc = otherpt.arc().clone();
                                        let oparam = otherpt.parameter();
                                        // OCCT L868-869: Vtgrst uses the outer
                                        // d1u/d1v (surface derivatives at the
                                        // current point) with the other arc's d2d.
                                        let od2d = oarc.derivative_at(oparam);
                                        let ov_tgrst = od2d.x * d1u + od2d.y * d1v;
                                        let mut otransline = Transition::new();
                                        let mut otransarc = Transition::new();
                                        if normale.length_squared() < 1e-16 {
                                            otransline.set_value_in_out(true, TypeTrans::Undecided);
                                            otransarc.set_value_in_out(true, TypeTrans::Undecided);
                                        } else {
                                            make_transition(
                                                v_tgt_int,
                                                ov_tgrst,
                                                normale,
                                                &mut otransline,
                                                &mut otransarc,
                                            );
                                        }
                                        let mut ointpt = intpt.clone();
                                        // OCCT L859-864: SetVertex on the shared vertex.
                                        ointpt.set_vertex(on_first);
                                        ointpt.set_arc(on_first, oarc, oparam, otransline, otransarc);
                                        ointpt.tolerance = the_toler;
                                        add_vertex(&mut slin[ii], ointpt);
                                        used_line[ii] = 1;
                                        retvalue = true;
                                        localdone[k] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            jj += 1;
        }
        ii += 1;
    }

    for ii in 0..nbpnt {
        done[ii] = localdone[ii];
    }

    retvalue
}

/// OCCT PointOnSecondDom (L918-1079).
#[allow(clippy::too_many_arguments)]
fn point_on_second_dom(
    listpnt: &[PathPoint],
    domain: &Domain,
    quad_surf: &Quadric,
    normale: DVec3,
    v_tgt_int: DVec3,
    lin: &mut IntPatchLine,
    done: &mut [i32],
    index: usize,
    the_toler: f64,
) -> bool {
    let currentpointonrst = &listpnt[index - 1];
    let mut retvalue = true;
    let nbpnt = listpnt.len();

    // OCCT L962-963, L1069-1076: nbvtx is re-read after each iteration because
    // the k-loop may AddVertex to the same line.
    let mut jj = 1usize;
    let mut nbvtx = nb_vertex(lin);
    while jj <= nbvtx {
        let mut intpt = vertex_of(lin, jj);
        if !intpt.on_dom_s2 {
            if currentpointonrst.value().distance(intpt.p3d) <= intpt.tolerance {
                // OCCT L977: Retvalue = false (a vertex of the second domain
                // already matches; the caller must not place a new one).
                retvalue = false;
                let goon = !currentpointonrst.is_new();
                let currentarc = currentpointonrst.arc().clone();
                let currentparameter = currentpointonrst.parameter();
                let p2d = currentarc.point_at(currentparameter);
                let d2d = currentarc.derivative_at(currentparameter);
                let (_ptbid, d1u, d1v) = quad_surf.d1(p2d.x, p2d.y);
                let v_tgrst = d2d.x * d1u + d2d.y * d1v;
                let mut transline = Transition::new();
                let mut transarc = Transition::new();
                if normale.length_squared() < 1e-16 {
                    transline.set_value_in_out(true, TypeTrans::Undecided);
                    transarc.set_value_in_out(true, TypeTrans::Undecided);
                } else {
                    make_transition(v_tgt_int, v_tgrst, normale, &mut transline, &mut transarc);
                }
                // OCCT L978-982: SetVertex(false) when the point is a domain vertex.
                if !currentpointonrst.is_new() {
                    intpt.set_vertex(false);
                }
                intpt.set_arc(false, currentarc, currentparameter, transline, transarc);
                intpt.tolerance = the_toler;
                replace_vertex(lin, jj, intpt.clone());
                done[index - 1] = 1;

                if goon {
                    for k in index..nbpnt {
                        if done[k] != 1 {
                            let otherpt = &listpnt[k];
                            if !otherpt.is_new() {
                                let vtxbis = otherpt.vertex();
                                if domain.identical(currentpointonrst.vertex(), vtxbis) {
                                    let oarc = otherpt.arc().clone();
                                    let oparam = otherpt.parameter();
                                    // OCCT L1030-1031: Vtgrst uses the outer
                                    // d1u/d1v with the other arc's d2d.
                                    let od2d = oarc.derivative_at(oparam);
                                    let ov_tgrst = od2d.x * d1u + od2d.y * d1v;
                                    let mut otransline = Transition::new();
                                    let mut otransarc = Transition::new();
                                    if normale.length_squared() < 1e-16 {
                                        otransline.set_value_in_out(true, TypeTrans::Undecided);
                                        otransarc.set_value_in_out(true, TypeTrans::Undecided);
                                    } else {
                                        make_transition(
                                            v_tgt_int,
                                            ov_tgrst,
                                            normale,
                                            &mut otransline,
                                            &mut otransarc,
                                        );
                                    }
                                    let mut ointpt = intpt.clone();
                                    // OCCT L1022-1027: SetVertex(false) on the
                                    // shared domain vertex.
                                    ointpt.set_vertex(false);
                                    ointpt.set_arc(false, oarc, oparam, otransline, otransarc);
                                    ointpt.tolerance = the_toler;
                                    add_vertex(lin, ointpt);
                                    done[k] = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
        jj += 1;
        nbvtx = nb_vertex(lin);
    }

    retvalue
}

/// OCCT PutPointsOnLine (L439-724).  Traitement des points de depart: replaces
/// the boundary points onto the intersection lines with the correct
/// transition.
#[allow(clippy::too_many_arguments)]
pub fn put_points_on_line(
    s1: &Surface3,
    s2: &Surface3,
    listpnt: &[PathPoint],
    slin: &mut Vec<IntPatchLine>,
    on_first: bool,
    domain: &Domain,
    quad_surf: &Quadric,
    other_quad: &Quadric,
    multpoint: bool,
    tolarc: f64,
) {
    let nbpnt = listpnt.len();
    let nblin = slin.len();

    if slin.is_empty() || nbpnt == 0 {
        return;
    }

    let mut done = vec![0i32; nbpnt];
    let mut used_line = vec![0i32; nblin];

    for i in 0..nbpnt {
        if done[i] != 1 {
            let mut currentpointonrst = listpnt[i].clone();
            let mut psurf = currentpointonrst.value();
            let mut tolerance = currentpointonrst.tolerance();

            // Search first for a match with a "multiple point".
            for lu in used_line.iter_mut() {
                *lu = 0;
            }
            let mut goon = true;
            if multpoint {
                let normale = quad_surf.normale(psurf);
                let currentarc = currentpointonrst.arc().clone();
                let currentparameter = currentpointonrst.parameter();
                let p2d = currentarc.point_at(currentparameter);
                let d2d = currentarc.derivative_at(currentparameter);
                let (_ptbid, d1u, d1v) = quad_surf.d1(p2d.x, p2d.y);
                let v_tgrst = d2d.x * d1u + d2d.y * d1v;
                let _ = v_tgrst;
                goon = multiple_point(
                    listpnt,
                    domain,
                    quad_surf,
                    normale,
                    slin,
                    &mut done,
                    &mut used_line,
                    i + 1,
                    on_first,
                    tolarc,
                );
            }
            if goon {
                let mut linefound = false;

                for indiceline in 0..slin.len() {
                    if used_line[indiceline] != 0 {
                        continue;
                    }
                    let mut linenumber = indiceline;

                    // Points may have been moved; retake the original point.
                    currentpointonrst = listpnt[i].clone();
                    let currentarc = currentpointonrst.arc().clone();
                    let currentparameter = currentpointonrst.parameter();
                    psurf = currentpointonrst.value();
                    tolerance = currentpointonrst.tolerance();

                    // OCC4455: enlarge tolerance from the vertex resolution.
                    if !currentpointonrst.is_new() {
                        let a_vtx_tol =
                            domain.vertex_tolerance(currentpointonrst.vertex(), &currentarc);
                        let a_tol_ang = 0.01 * tolerance;
                        tolerance = tolerance.max(a_vtx_tol);
                        let a_norm1 = quad_surf.normale(psurf);
                        let a_norm2 = other_quad.normale(psurf);
                        if a_norm1.length() > f64::MIN_POSITIVE && a_norm2.length() > f64::MIN_POSITIVE
                        {
                            if gp_vec_is_parallel(a_norm1, a_norm2, a_tol_ang) {
                                tolerance = tolerance.sqrt();
                            }
                        }
                    }

                    let mut pointonarc = DVec3::ZERO;
                    let mut v_tgt_int = DVec3::ZERO;
                    let mut a_l_params: Vec<f64> = Vec::new();
                    let mut a_vert_tol = tolarc;
                    let mut param_on_arc = currentparameter;
                    linefound = find_line(
                        &mut psurf,
                        slin,
                        tolerance,
                        &mut a_l_params,
                        &mut v_tgt_int,
                        &mut linenumber,
                        indiceline + 1,
                        &currentarc,
                        &mut param_on_arc,
                        &mut pointonarc,
                        quad_surf,
                        other_quad,
                        &mut a_vert_tol,
                    );
                    let linenumber = linenumber; // &mut -> value for indexing below
                    if linefound {
                        let normale = quad_surf.normale(psurf);
                        let currentarc = currentpointonrst.arc().clone();
                        let p2d = currentarc.point_at(currentparameter);
                        let d2d = currentarc.derivative_at(currentparameter);
                        let (_ptbid, d1u, d1v) = quad_surf.d1(p2d.x, p2d.y);
                        let v_tgrst = d2d.x * d1u + d2d.y * d1v;

                        let the_type = slin[linenumber].line_type;

                        if !on_first {
                            // Match between the point on the first domain and
                            // the point on the second domain.
                            goon = point_on_second_dom(
                                listpnt,
                                domain,
                                quad_surf,
                                normale,
                                v_tgt_int,
                                &mut slin[linenumber],
                                &mut done,
                                i + 1,
                                a_vert_tol,
                            );
                        }

                        if goon {
                            let mut solpnt = IntPatchVertex::default();
                            solpnt.set_value(psurf, a_vert_tol, false);

                            let u1 = p2d.x;
                            let v1 = p2d.y;
                            let (u2, v2) = other_quad.parameters(psurf);

                            if on_first {
                                let (mut u1, mut v1, mut u2, mut v2) = (u1, v1, u2, v2);
                                recadre(s1, s2, &mut u1, &mut v1, &mut u2, &mut v2);
                                solpnt.set_parameters(u1, v1, u2, v2);
                            } else {
                                let (mut u1, mut v1, mut u2, mut v2) = (u1, v1, u2, v2);
                                recadre(s1, s2, &mut u2, &mut v2, &mut u1, &mut v1);
                                solpnt.set_parameters(u2, v2, u1, v1);
                            }

                            // OCCT L630-638: SetVertex when the point is a domain
                            // vertex of the current surface.
                            if !currentpointonrst.is_new() {
                                solpnt.set_vertex(on_first);
                            }

                            let mut transline = Transition::new();
                            let mut transarc = Transition::new();
                            if normale.length_squared() < 1e-16 {
                                transline.set_value_in_out(true, TypeTrans::Undecided);
                                transarc.set_value_in_out(true, TypeTrans::Undecided);
                            } else {
                                make_transition(v_tgt_int, v_tgrst, normale, &mut transline, &mut transarc);
                            }
                            solpnt.set_arc(on_first, currentarc, currentparameter, transline, transarc);

                            for &par in &a_l_params {
                                solpnt.set_parameter(par);
                                add_vertex(&mut slin[linenumber], solpnt.clone());
                            }

                            done[i] = 1;

                            // OCCT L666-712: the same solpnt is placed on the
                            // other path points that share the current point's
                            // domain vertex.
                            if goon && !currentpointonrst.is_new() {
                                let vtx = currentpointonrst.vertex();
                                for k in (i + 1)..nbpnt {
                                    if done[k] != 1 {
                                        let otherpt = &listpnt[k];
                                        if !otherpt.is_new() {
                                            let vtxbis = otherpt.vertex();
                                            if domain.identical(vtx, vtxbis) {
                                                // OCCT L681-698: reuse solpnt, only
                                                // the arc changes (d1u/d1v from the
                                                // current point).  The shared
                                                // domain vertex is kept.
                                                solpnt.set_tolerance(tolarc);
                                                solpnt.set_vertex(on_first);
                                                let karc = otherpt.arc().clone();
                                                let kparam = otherpt.parameter();
                                                let kd2d = karc.derivative_at(kparam);
                                                let k_v_tgrst = kd2d.x * d1u + kd2d.y * d1v;
                                                let mut ktransline = Transition::new();
                                                let mut ktransarc = Transition::new();
                                                if normale.length_squared() < 1e-16 {
                                                    ktransline
                                                        .set_value_in_out(true, TypeTrans::Undecided);
                                                    ktransarc
                                                        .set_value_in_out(true, TypeTrans::Undecided);
                                                } else {
                                                    make_transition(
                                                        v_tgt_int,
                                                        k_v_tgrst,
                                                        normale,
                                                        &mut ktransline,
                                                        &mut ktransarc,
                                                    );
                                                }
                                                solpnt.set_arc(
                                                    on_first,
                                                    karc,
                                                    kparam,
                                                    ktransline,
                                                    ktransarc,
                                                );
                                                add_vertex(&mut slin[linenumber], solpnt.clone());
                                                done[k] = 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        done[i] = 1;
                    }
                }
            }
        }
    }
}

/// OCCT ProcessSegments (L1486-1885).  Converts the solution segments into
/// IntPatch_RLine restriction lines.
#[allow(clippy::too_many_arguments)]
pub fn process_segments(
    listedg: &[super::so_on_bounds::Segment],
    slin: &mut Vec<IntPatchLine>,
    quad1: &Quadric,
    quad2: &Quadric,
    on_first: bool,
    tol_arc: f64,
) {
    for seg in listedg {
        let mut edge_degenere = false;
        let arc_ref = seg.curve().clone();

        let mut dofirst = false;
        let mut dolast = false;
        let mut procf = false;
        let mut procl = false;

        let mut paramf = 0.0;
        let mut paraml = 0.0;
        let mut pstartf: Option<&PathPoint> = None;
        let mut pstartl: Option<&PathPoint> = None;

        if seg.has_first_point() {
            dofirst = true;
            pstartf = Some(seg.first_point());
            paramf = pstartf.unwrap().parameter();
        }
        if seg.has_last_point() {
            dolast = true;
            pstartl = Some(seg.last_point());
            paraml = pstartl.unwrap().parameter();
        }

        let mut rline = IntPatchLine::analytic(
            IntPatchIType::Restriction,
            Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::X,
            }),
            [paramf, paraml],
        );
        if on_first {
            rline.set_arc_on_s1(arc_ref.clone());
        } else {
            rline.set_arc_on_s2(arc_ref.clone());
        }

        if dofirst && dolast {
            // Determination of the transition of the line.
            let p2d = arc_ref.point_at(0.5 * (paramf + paraml));
            let d2d = arc_ref.derivative_at(0.5 * (paramf + paraml));
            let (_valpt, d1u, d1v) = if on_first {
                quad1.d1(p2d.x, p2d.y)
            } else {
                quad2.d1(p2d.x, p2d.y)
            };
            let tgline = d2d.x * d1u + d2d.y * d1v;

            if d1u.length() < 1e-7 {
                // Degenerate edge?
                edge_degenere = true;
                for edg in 0..=10 {
                    let ep2d = arc_ref.point_at(paramf + (paraml - paramf) * edg as f64 * 0.1);
                    let ed2d = arc_ref.derivative_at(paramf + (paraml - paramf) * edg as f64 * 0.1);
                    let (_evalpt, ed1u, _ed1v) = if on_first {
                        quad1.d1(ep2d.x, ep2d.y)
                    } else {
                        quad2.d1(ep2d.x, ep2d.y)
                    };
                    let _ = ed2d;
                    if ed1u.length() > 1e-7 {
                        edge_degenere = false;
                    }
                }
                rline = IntPatchLine::analytic(
                    IntPatchIType::Restriction,
                    Curve3::Line(rcad_kernel::geom::Line3 {
                        origin: DVec3::ZERO,
                        direction: DVec3::X,
                    }),
                    [paramf, paraml],
                );
                if on_first {
                    rline.set_arc_on_s1(arc_ref.clone());
                } else {
                    rline.set_arc_on_s2(arc_ref.clone());
                }
            } else {
                let norm2 = quad2.normale(_valpt);
                let norm1 = quad1.normale(_valpt);
                let trans1;
                let trans2;
                if tgline.cross(norm2).dot(norm1) > 0.000000001 {
                    trans1 = TypeTrans::Out;
                    trans2 = TypeTrans::In;
                } else if tgline.cross(norm2).dot(norm1) < -0.000000001 {
                    trans1 = TypeTrans::In;
                    trans2 = TypeTrans::Out;
                } else {
                    trans1 = TypeTrans::Undecided;
                    trans2 = TypeTrans::Undecided;
                }
                rline = IntPatchLine::analytic(
                    IntPatchIType::Restriction,
                    Curve3::Line(rcad_kernel::geom::Line3 {
                        origin: DVec3::ZERO,
                        direction: DVec3::X,
                    }),
                    [paramf, paraml],
                );
                rline.trans1 = Some(Transition::new_in_out(false, trans1));
                rline.trans2 = Some(Transition::new_in_out(false, trans2));
                if on_first {
                    rline.set_arc_on_s1(arc_ref.clone());
                } else {
                    rline.set_arc_on_s2(arc_ref.clone());
                }
            }
        } else {
            rline = IntPatchLine::analytic(
                IntPatchIType::Restriction,
                Curve3::Line(rcad_kernel::geom::Line3 {
                    origin: DVec3::ZERO,
                    direction: DVec3::X,
                }),
                [paramf, paraml],
            );
            if on_first {
                rline.set_arc_on_s1(arc_ref.clone());
            } else {
                rline.set_arc_on_s2(arc_ref.clone());
            }
        }

        if dofirst || dolast {
            let nblines = slin.len();
            for j in 0..nblines {
                let typ = slin[j].line_type;
                let nbpts = nb_vertex(&slin[j]);
                for k in 1..=nbpts {
                    let mut ptvtx = vertex_of(&slin[j], k);

                    if !edge_degenere && dofirst {
                        if ptvtx.p3d.distance(pstartf.unwrap().value()) <= tol_arc {
                            ptvtx.set_multiple(true);
                            ptvtx.set_tolerance(tol_arc);
                            replace_vertex(&mut slin[j], k, ptvtx.clone());
                            let mut newptvtx = ptvtx.clone();
                            newptvtx.set_parameter(paramf);

                            let p2d = arc_ref.point_at(paramf);
                            let d2d = arc_ref.derivative_at(paramf);
                            let (_valpt, d1u, d1v) = if on_first {
                                quad1.d1(p2d.x, p2d.y)
                            } else {
                                quad2.d1(p2d.x, p2d.y)
                            };
                            let tgline = d2d.x * d1u + d2d.y * d1v;

                            if ptvtx.on_dom_s1 {
                                if let Some(thearc) = ptvtx.arc_on_s1.clone() {
                                    let ap2d = thearc.point_at(ptvtx.parameter_on_arc1());
                                    let ad2d = thearc.derivative_at(ptvtx.parameter_on_arc1());
                                    let (_avalpt, ad1u, ad1v) = quad1.d1(ap2d.x, ap2d.y);
                                    let tgarc = ad2d.x * ad1u + ad2d.y * ad1v;
                                    let norm1 = ad1u.cross(ad1v);
                                    let mut trest = Transition::new();
                                    let mut tarc = Transition::new();
                                    if norm1.length_squared() < 1e-16 {
                                        trest.set_value_in_out(true, TypeTrans::Undecided);
                                        tarc.set_value_in_out(true, TypeTrans::Undecided);
                                    } else {
                                        make_transition(tgline, tgarc, norm1, &mut trest, &mut tarc);
                                    }
                                    newptvtx.set_arc(true, thearc, ptvtx.parameter_on_arc1(), trest, tarc);
                                }
                            }
                            if ptvtx.on_dom_s2 {
                                if let Some(thearc) = ptvtx.arc_on_s2.clone() {
                                    let ap2d = thearc.point_at(ptvtx.parameter_on_arc2());
                                    let ad2d = thearc.derivative_at(ptvtx.parameter_on_arc2());
                                    let (_avalpt, ad1u, ad1v) = quad2.d1(ap2d.x, ap2d.y);
                                    let tgarc = ad2d.x * ad1u + ad2d.y * ad1v;
                                    let norm2 = ad1u.cross(ad1v);
                                    let mut trest = Transition::new();
                                    let mut tarc = Transition::new();
                                    if norm2.length_squared() < 1e-16 {
                                        trest.set_value_in_out(true, TypeTrans::Undecided);
                                        tarc.set_value_in_out(true, TypeTrans::Undecided);
                                    } else {
                                        make_transition(tgline, tgarc, norm2, &mut trest, &mut tarc);
                                    }
                                    newptvtx.set_arc(false, thearc, ptvtx.parameter_on_arc2(), trest, tarc);
                                }
                            }

                            add_vertex(&mut rline, newptvtx.clone());
                            if !procf {
                                procf = true;
                                rline.set_first_point(rline.vertices.len());
                            }
                        }
                    }
                    if !edge_degenere && dolast {
                        if ptvtx.p3d.distance(pstartl.unwrap().value()) <= tol_arc {
                            ptvtx.set_multiple(true);
                            ptvtx.set_tolerance(tol_arc);
                            replace_vertex(&mut slin[j], k, ptvtx.clone());
                            let mut newptvtx = ptvtx.clone();
                            newptvtx.set_parameter(paraml);

                            let p2d = arc_ref.point_at(paraml);
                            let d2d = arc_ref.derivative_at(paraml);
                            let (_valpt, d1u, d1v) = if on_first {
                                quad1.d1(p2d.x, p2d.y)
                            } else {
                                quad2.d1(p2d.x, p2d.y)
                            };
                            let tgline = d2d.x * d1u + d2d.y * d1v;

                            if ptvtx.on_dom_s1 {
                                if let Some(thearc) = ptvtx.arc_on_s1.clone() {
                                    let ap2d = thearc.point_at(ptvtx.parameter_on_arc1());
                                    let ad2d = thearc.derivative_at(ptvtx.parameter_on_arc1());
                                    let (_avalpt, ad1u, ad1v) = quad1.d1(ap2d.x, ap2d.y);
                                    let tgarc = ad2d.x * ad1u + ad2d.y * ad1v;
                                    let norm1 = ad1u.cross(ad1v);
                                    let mut trest = Transition::new();
                                    let mut tarc = Transition::new();
                                    if norm1.length_squared() < 1e-16 {
                                        trest.set_value_in_out(true, TypeTrans::Undecided);
                                        tarc.set_value_in_out(true, TypeTrans::Undecided);
                                    } else {
                                        make_transition(tgline, tgarc, norm1, &mut trest, &mut tarc);
                                    }
                                    newptvtx.set_arc(true, thearc, ptvtx.parameter_on_arc1(), trest, tarc);
                                }
                            }
                            if ptvtx.on_dom_s2 {
                                if let Some(thearc) = ptvtx.arc_on_s2.clone() {
                                    let ap2d = thearc.point_at(ptvtx.parameter_on_arc2());
                                    let ad2d = thearc.derivative_at(ptvtx.parameter_on_arc2());
                                    let (_avalpt, ad1u, ad1v) = quad2.d1(ap2d.x, ap2d.y);
                                    let tgarc = ad2d.x * ad1u + ad2d.y * ad1v;
                                    let norm2 = ad1u.cross(ad1v);
                                    let mut trest = Transition::new();
                                    let mut tarc = Transition::new();
                                    if norm2.length_squared() < 1e-16 {
                                        trest.set_value_in_out(true, TypeTrans::Undecided);
                                        tarc.set_value_in_out(true, TypeTrans::Undecided);
                                    } else {
                                        make_transition(tgline, tgarc, norm2, &mut trest, &mut tarc);
                                    }
                                    newptvtx.set_arc(false, thearc, ptvtx.parameter_on_arc2(), trest, tarc);
                                }
                            }

                            add_vertex(&mut rline, newptvtx.clone());
                            if !procl {
                                procl = true;
                                rline.set_last_point(rline.vertices.len());
                            }
                        }
                    }
                }
                if procf {
                    dofirst = false;
                }
                if procl {
                    dolast = false;
                }
            }
        }

        // If the first/last point was not found on a line, still place it on
        // the restriction solution.
        if dofirst {
            let mut ptvtx = IntPatchVertex::default();
            let pstartf = pstartf.unwrap();
            ptvtx.set_value(pstartf.value(), pstartf.tolerance(), false);
            let (u1, v1) = quad1.parameters(pstartf.value());
            let (u2, v2) = quad2.parameters(pstartf.value());
            ptvtx.set_parameters(u1, v1, u2, v2);
            ptvtx.set_parameter(paramf);
            // OCCT L1853-1859: SetVertex + set the arc when the point is on a
            // domain vertex.
            if !pstartf.is_new() {
                ptvtx.set_vertex(on_first);
                ptvtx.set_arc(
                    on_first,
                    pstartf.arc().clone(),
                    pstartf.parameter(),
                    Transition::new(),
                    Transition::new(),
                );
            }
            add_vertex(&mut rline, ptvtx.clone());
            rline.set_first_point(rline.vertices.len());
        }
        if dolast {
            let mut ptvtx = IntPatchVertex::default();
            let pstartl = pstartl.unwrap();
            ptvtx.set_value(pstartl.value(), pstartl.tolerance(), false);
            let (u1, v1) = quad1.parameters(pstartl.value());
            let (u2, v2) = quad2.parameters(pstartl.value());
            ptvtx.set_parameters(u1, v1, u2, v2);
            ptvtx.set_parameter(paraml);
            // OCCT L1871-1877: SetVertex + set the arc when the point is on a
            // domain vertex.
            if !pstartl.is_new() {
                ptvtx.set_vertex(on_first);
                ptvtx.set_arc(
                    on_first,
                    pstartl.arc().clone(),
                    pstartl.parameter(),
                    Transition::new(),
                    Transition::new(),
                );
            }
            add_vertex(&mut rline, ptvtx.clone());
            rline.set_last_point(rline.vertices.len());
        }
        slin.push(rline);
    }
}

/// OCCT SquareDistance (L1897-1931) — distance from a GLine to a point.  For
/// the Line/Circle types the analytic quadric distance is used; for the other
/// curve types (Ellipse/Parabola/Hyperbola) OCCT runs Extrema_ExtPC (rcad:
/// base::extrema::ExtPC) over the curve's parameter range.
fn square_distance_gline(line: &IntPatchLine, p: DVec3) -> f64 {
    match &line.curve {
        Curve3::Line(l) => {
            let pp = l.origin + l.direction.normalize_or_zero()
                * (p - l.origin).dot(l.direction.normalize_or_zero());
            pp.distance_squared(p)
        }
        Curve3::Circle(c) => {
            let d = p - c.center;
            let radial = d - c.normal.normalize_or_zero()
                * d.dot(c.normal.normalize_or_zero());
            let pp = c.center + radial.normalize_or_zero() * c.radius;
            pp.distance_squared(p)
        }
        _ => {
            // OCCT IsRLineGood L1978-1983: anExtr is initialized with the
            // GLine curve's FirstParameter/LastParameter.
            let d = line.curve.default_domain();
            let ext = rcad_kernel::base::extrema::ExtPC::new(p, &line.curve, 1e-7, d[0], d[1]);
            if !ext.is_done() || ext.nb_ext() == 0 {
                return f64::MAX;
            }
            let mut sq = ext.square_distance(1);
            for i in 2..=ext.nb_ext() {
                sq = sq.min(ext.square_distance(i));
            }
            sq
        }
    }
}

/// OCCT IsRLineGood (L1933-2027) — the RLine is discarded when all its points
/// lie on the GLine within the tolerance.
fn is_r_line_good(
    quad1: &Quadric,
    quad2: &Quadric,
    gline: &IntPatchLine,
    rline: &IntPatchLine,
    the_tol: f64,
) -> bool {
    let a_sq_tol = the_tol * the_tol;
    let a_nb_pnts_m1 = if rline.vertices.len() > 0 {
        rline.vertices.len() - 1
    } else {
        0
    };

    if a_nb_pnts_m1 < 1 {
        return false;
    }

    if a_nb_pnts_m1 == 1 {
        let a_p1 = rline.vertices[0].p3d;
        let a_p2 = rline.vertices[1].p3d;
        if a_p1.distance_squared(a_p2) < a_sq_tol {
            // RLine is degenerated.
            return false;
        }
        let a_pmid;
        if rline.is_arc_on_s1() {
            let an_ac2d = rline.arc_on_s1().unwrap();
            // OCCT uses the bounded arc parameters (anAC2d->FirstParameter()/
            // LastParameter()); the rcad FF arcs are unbounded 2D lines, so the
            // range is taken from the RLine first/last points instead.
            let a_par_f = if rline.has_first_point() {
                rline.first_point().parameter_on_line()
            } else {
                an_ac2d.default_domain()[0]
            };
            let a_par_l = if rline.has_last_point() {
                rline.last_point().parameter_on_line()
            } else {
                an_ac2d.default_domain()[1]
            };
            let a_p2d = an_ac2d.point_at(0.5 * (a_par_f + a_par_l));
            a_pmid = quad1.value(a_p2d.x, a_p2d.y);
        } else {
            let an_ac2d = rline.arc_on_s2().unwrap();
            let a_par_f = if rline.has_first_point() {
                rline.first_point().parameter_on_line()
            } else {
                an_ac2d.default_domain()[0]
            };
            let a_par_l = if rline.has_last_point() {
                rline.last_point().parameter_on_line()
            } else {
                an_ac2d.default_domain()[1]
            };
            let a_p2d = an_ac2d.point_at(0.5 * (a_par_f + a_par_l));
            a_pmid = quad2.value(a_p2d.x, a_p2d.y);
        }
        let a_sq_dist = square_distance_gline(gline, a_pmid);
        return a_sq_dist > a_sq_tol;
    }

    for i in 1..a_nb_pnts_m1 {
        let a_p = rline.vertices[i].p3d;
        let a_sq_dist = square_distance_gline(gline, a_p);
        if a_sq_dist > a_sq_tol {
            return true;
        }
    }
    false
}

/// OCCT ProcessRLine (L2029-2343).  Places the "multiple" points of the other
/// intersection lines onto the restriction lines.
pub fn process_r_line(
    slin: &mut Vec<IntPatchLine>,
    quad1: &Quadric,
    quad2: &Quadric,
    _tol_arc: f64,
    the_is_req_to_keep_r_line: bool,
) {
    let tol_arc = (100.0 * _tol_arc).min(0.1);

    let mut i = 0usize;
    while i < slin.len() {
        let typ1 = slin[i].line_type;
        let mut has_to_delete_r_line = false;
        if typ1 == IntPatchIType::Restriction {
            let mut seq_pnt3d: Vec<DVec3> = Vec::new();
            let mut seq_real: Vec<f64> = Vec::new();

            let rline_i = slin[i].clone();
            for j in 0..slin.len() {
                let nbpt = seq_pnt3d.len();
                let typ2 = slin[j].line_type;
                if typ2 != IntPatchIType::Restriction {
                    let (on_first, arcref) = if rline_i.is_arc_on_s1() {
                        (true, rline_i.arc_on_s1().unwrap().clone())
                    } else if rline_i.is_arc_on_s2() {
                        (false, rline_i.arc_on_s2().unwrap().clone())
                    } else {
                        continue;
                    };
                    let paramf = if rline_i.has_first_point() {
                        rline_i.first_point().parameter_on_line()
                    } else {
                        f64::NEG_INFINITY
                    };
                    let paraml = if rline_i.has_last_point() {
                        rline_i.last_point().parameter_on_line()
                    } else {
                        f64::INFINITY
                    };

                    let nbvtx = nb_vertex(&slin[j]);

                    // Edge degenerate check.
                    let mut edge_degenere = true;
                    let mut edg = 0;
                    while edge_degenere && edg <= 10 {
                        let ep2d = arcref.point_at(paramf + (paraml - paramf) * edg as f64 * 0.1);
                        let (_evalpt, ed1u, _ed1v) = if on_first {
                            quad1.d1(ep2d.x, ep2d.y)
                        } else {
                            quad2.d1(ep2d.x, ep2d.y)
                        };
                        if ed1u.length() > 1e-7 {
                            edge_degenere = false;
                        }
                        edg += 1;
                    }

                    let mut k = 1usize;
                    while !edge_degenere && k <= nbvtx {
                        let mut ptvtx = vertex_of(&slin[j], k);
                        if (on_first && !ptvtx.on_dom_s1) || (!on_first && !ptvtx.on_dom_s2) {
                            let mut project = true;
                            let mut keeppoint = false;
                            let toproj = ptvtx.p3d;

                            let mut jj = 0usize;
                            while jj < nbpt {
                                if toproj.distance(seq_pnt3d[jj]) < _tol_arc {
                                    project = false;
                                    break;
                                }
                                jj += 1;
                            }
                            if project {
                                let (u, v) = if on_first {
                                    ptvtx.parameters_on_s1()
                                } else {
                                    ptvtx.parameters_on_s2()
                                };
                                let p2d_input = DVec2::new(u, v);
                                let (proj_ok, paramproj, p2d_proj) =
                                    project_2d(&arcref, p2d_input, paramf, paraml);
                                if proj_ok {
                                    let ptproj = if on_first {
                                        quad1.value(p2d_proj.x, p2d_proj.y)
                                    } else {
                                        quad2.value(p2d_proj.x, p2d_proj.y)
                                    };
                                    if toproj.distance(ptproj) <= 100.0 * tol_arc
                                        && paramproj >= paramf
                                        && paramproj <= paraml
                                    {
                                        let mut newptvtx = ptvtx.clone();
                                        newptvtx.set_parameter(paramproj);
                                        keeppoint = true;
                                        seq_pnt3d.push(toproj);
                                        seq_real.push(paramproj);

                                        // Verify that the restriction carries this vertex.
                                        for ri in 0..slin.len() {
                                            if slin[ri].line_type == IntPatchIType::Restriction {
                                                if on_first && slin[ri].is_arc_on_s1() {
                                                    if super::so_on_bounds::curves_same(&arcref, slin[ri].arc_on_s1().unwrap()) {
                                                        add_vertex(&mut slin[ri], newptvtx.clone());
                                                    }
                                                } else if !on_first && slin[ri].is_arc_on_s2() {
                                                    if super::so_on_bounds::curves_same(&arcref, slin[ri].arc_on_s2().unwrap()) {
                                                        add_vertex(&mut slin[ri], newptvtx.clone());
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            } else {
                                keeppoint = true;
                                let mut newptvtx = ptvtx.clone();
                                if jj >= 1 && jj <= seq_real.len() {
                                    newptvtx.set_parameter(seq_real[jj - 1]);
                                }
                            }
                            if keeppoint {
                                ptvtx.set_multiple(true);
                                ptvtx.set_tolerance(_tol_arc);
                                let mut newptvtx = ptvtx.clone();
                                newptvtx.set_multiple(true);
                                replace_vertex(&mut slin[j], k, ptvtx.clone());

                                let mut tgrest = DVec3::ZERO;
                                if ptvtx.on_dom_s1 || ptvtx.on_dom_s2 {
                                    let p2d = arcref.point_at(newptvtx.parameter_on_line());
                                    let d2d = arcref.derivative_at(newptvtx.parameter_on_line());
                                    if on_first {
                                        // donc OnDomS2
                                        let (_valpt, d1u, d1v) = quad1.d1(p2d.x, p2d.y);
                                        tgrest = d2d.x * d1u + d2d.y * d1v;
                                        if let Some(thearc) = ptvtx.arc_on_s2.clone() {
                                            let ap2d = thearc.point_at(ptvtx.parameter_on_arc2());
                                            let ad2d = thearc.derivative_at(ptvtx.parameter_on_arc2());
                                            let (_avalpt, ad1u, ad1v) = quad2.d1(ap2d.x, ap2d.y);
                                            let tgarc = ad2d.x * ad1u + ad2d.y * ad1v;
                                            let norm = ad1u.cross(ad1v);
                                            let mut trest = Transition::new();
                                            let mut tarc = Transition::new();
                                            if norm.length_squared() < 1e-16 {
                                                trest.set_value_in_out(true, TypeTrans::Undecided);
                                                tarc.set_value_in_out(true, TypeTrans::Undecided);
                                            } else {
                                                make_transition(tgrest, tgarc, norm, &mut trest, &mut tarc);
                                            }
                                            newptvtx.set_arc(false, thearc, ptvtx.parameter_on_arc2(), trest, tarc);
                                        }
                                    } else {
                                        // donc OnDomS1
                                        let (_valpt, d1u, d1v) = quad2.d1(p2d.x, p2d.y);
                                        tgrest = d2d.x * d1u + d2d.y * d1v;
                                        if let Some(thearc) = ptvtx.arc_on_s1.clone() {
                                            let ap2d = thearc.point_at(ptvtx.parameter_on_arc1());
                                            let ad2d = thearc.derivative_at(ptvtx.parameter_on_arc1());
                                            let (_avalpt, ad1u, ad1v) = quad1.d1(ap2d.x, ap2d.y);
                                            let tgarc = ad2d.x * ad1u + ad2d.y * ad1v;
                                            let norm = ad1u.cross(ad1v);
                                            let mut trest = Transition::new();
                                            let mut tarc = Transition::new();
                                            if norm.length_squared() < 1e-16 {
                                                trest.set_value_in_out(true, TypeTrans::Undecided);
                                                tarc.set_value_in_out(true, TypeTrans::Undecided);
                                            } else {
                                                make_transition(tgrest, tgarc, norm, &mut trest, &mut tarc);
                                            }
                                            newptvtx.set_arc(true, thearc, ptvtx.parameter_on_arc1(), trest, tarc);
                                        }
                                    }
                                }
                                add_vertex(&mut slin[i], newptvtx);
                            }
                        }
                        k += 1;
                    }

                    if !the_is_req_to_keep_r_line {
                        // Discard the RLine when it is redundant with the GLine.
                        let a_g_line_type = slin[j].line_type;
                        let is_gline = a_g_line_type != IntPatchIType::Analytic
                            && a_g_line_type != IntPatchIType::Restriction;
                        if is_gline {
                            let gline = slin[j].clone();
                            let rline = slin[i].clone();
                            has_to_delete_r_line =
                                !is_r_line_good(quad1, quad2, &gline, &rline, tol_arc);
                        }
                        if has_to_delete_r_line {
                            break;
                        }
                    }
                }
            }
        }

        if has_to_delete_r_line {
            slin.remove(i);
            // i stays (OCCT: i-- then i++ in the loop).
        } else {
            i += 1;
        }
    }
}

/// OCCT IntPatch_HInterTool::Project (IntPatch_HInterTool.cxx L271-304) —
/// projects a 2D point onto the arc.  OCCT runs Extrema_EPCOfExtPC2d (a generic
/// point-to-2D-curve extrema solver with Nbu=20, epsX=1e-8, Tol=1e-5) on the
/// arc's bounded parameter range; rcad uses the ExtPC2d translation
/// (Extrema_ExtPC2d) which returns the closest point on the arc over
/// [u_inf, u_sup].
fn project_2d(c: &Curve2d, p: DVec2, u_inf: f64, u_sup: f64) -> (bool, f64, DVec2) {
    const TOL: f64 = 1.0e-5; // IntPatch_HInterTool::Project Tol
    let ext = rcad_kernel::base::extrema::ExtPC2d::new(p, c, TOL, u_inf, u_sup);
    if !ext.is_done() || ext.nb_ext() == 0 {
        return (false, 0.0, p);
    }
    let mut indexmin = 1usize;
    let mut dist2 = ext.square_distance(1);
    for i in 2..=ext.nb_ext() {
        if ext.square_distance(i) < dist2 {
            indexmin = i;
            dist2 = ext.square_distance(i);
        }
    }
    let pnt = ext.point(indexmin);
    (true, pnt.param, pnt.point)
}
