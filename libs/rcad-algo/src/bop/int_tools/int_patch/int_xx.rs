// OCCT IntPatch_ImpImpIntersection IntXX functions (L2345-9650) — the 15
// analytic pair intersection routines, translated 1:1 as free functions.
//
// Each IntXX takes the two quadrics (plus tolerances / reversal), the
// reference-solution flags (Empty / Same / Multpoint) as out-params, and
// appends to slin (Vec<IntPatchLine>) / spnt (Vec<IntPatchPoint>).
//
// rcad data-model notes:
//   - IntAna_QuadQuadGeo -> QuadQuadGeo (int_patch::quad_quad_geo).
//   - IntPatch_GLine (conic with transitions) -> IntPatchLine built by
//     gline_line/gline_circle/... which set line_type + curve + trans1/trans2.
//   - IntPatch_ALine (IntAna_Curve) -> IntPatchLine with a_curve set.
//   - IntSurf_TypeTrans / IntSurf_Situation -> transitions::{TypeTrans, Situation}.

use glam::DVec3;

use super::elclib as ecl;
use super::transitions::{Situation, Transition, TypeTrans};
use super::{IntPatchIType, IntPatchLine, IntPatchPoint};
use super::quad_quad_geo::{AnaResultType, QuadQuadGeo};
use super::special_points::{PatchPoint, PntOn2S};
use crate::topalgo::int_surf::quadric::Quadric;
use rcad_kernel::geom::{Curve3, CurveEval, Ellipse3, Hyperbola3, Line3, Parabola3};

/// OCCT IntPatch_GLine(Lin, Tang, Trans1, Trans2) — line with In/Out transitions.
fn gline_line(l: Line3, tang: bool, t1: TypeTrans, t2: TypeTrans) -> IntPatchLine {
    let mut line = IntPatchLine::analytic(
        IntPatchIType::Line,
        Curve3::Line(l),
        [f64::NEG_INFINITY, f64::INFINITY],
    );
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// OCCT IntPatch_GLine(Lin, Tang, Situ1, Situ2) — line with Touch transitions.
fn gline_line_touch(l: Line3, tang: bool, s1: Situation, s2: Situation) -> IntPatchLine {
    let mut line = IntPatchLine::analytic(
        IntPatchIType::Line,
        Curve3::Line(l),
        [f64::NEG_INFINITY, f64::INFINITY],
    );
    line.trans1 = Some(Transition::new_touch(tang, s1, false));
    line.trans2 = Some(Transition::new_touch(tang, s2, false));
    line
}

/// OCCT IntPatch_GLine(Circ, Tang, Trans1, Trans2).
fn gline_circle(c: rcad_kernel::geom::Circle3, tang: bool, t1: TypeTrans, t2: TypeTrans) -> IntPatchLine {
    let mut line = IntPatchLine::analytic(
        IntPatchIType::Circle,
        Curve3::Circle(c),
        [0.0, std::f64::consts::TAU],
    );
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// OCCT IntPatch_GLine(Circ, Tang, Situ1, Situ2).
fn gline_circle_touch(c: rcad_kernel::geom::Circle3, tang: bool, s1: Situation, s2: Situation) -> IntPatchLine {
    let mut line = IntPatchLine::analytic(
        IntPatchIType::Circle,
        Curve3::Circle(c),
        [0.0, std::f64::consts::TAU],
    );
    line.trans1 = Some(Transition::new_touch(tang, s1, false));
    line.trans2 = Some(Transition::new_touch(tang, s2, false));
    line
}

/// OCCT IntPatch_GLine(Elips, Tang, Trans1, Trans2).
fn gline_ellipse(e: Ellipse3, tang: bool, t1: TypeTrans, t2: TypeTrans) -> IntPatchLine {
    let mut line = IntPatchLine::analytic(
        IntPatchIType::Ellipse,
        Curve3::Ellipse(e),
        [0.0, std::f64::consts::TAU],
    );
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// OCCT IntPatch_GLine(Parab, Tang, Trans1, Trans2).
fn gline_parabola(p: Parabola3, tang: bool, t1: TypeTrans, t2: TypeTrans) -> IntPatchLine {
    let dom = p.default_domain();
    let mut line = IntPatchLine::analytic(IntPatchIType::Parabola, Curve3::Parabola(p), dom);
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// OCCT IntPatch_GLine(Hypr, Tang, Trans1, Trans2).
fn gline_hyperbola(h: Hyperbola3, tang: bool, t1: TypeTrans, t2: TypeTrans) -> IntPatchLine {
    let dom = h.default_domain();
    let mut line = IntPatchLine::analytic(IntPatchIType::Hyperbola, Curve3::Hyperbola(h), dom);
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// Build an IntPatch_Point (PatchPoint) for a 3D point on both quadrics.
fn make_point(p: DVec3, tol: f64, tangent: bool, q1: &Quadric, q2: &Quadric) -> PatchPoint {
    let (u1, v1) = q1.parameters(p);
    let (u2, v2) = q2.parameters(p);
    PatchPoint {
        pnt: PntOn2S { p, u1, v1, u2, v2 },
        param_on_line: 0.0,
        tolerance: tol,
        multiple: false,
        on_dom_s1: true,
        on_dom_s2: true,
        arc_on_s1: None,
        arc_on_s2: None,
        param_on_arc1: 0.0,
        param_on_arc2: 0.0,
    }
}

/// OCCT IntPP (L3106-3150) — Plane/Plane.
pub fn int_pp(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_ang: f64,
    tol_tang: f64,
    same: &mut bool,
    slin: &mut Vec<IntPatchLine>,
) -> bool {
    let pl1 = quad1.plane();
    let pl2 = quad2.plane();

    let mut inter = QuadQuadGeo::new();
    inter.perform_plane_plane(quad1, quad2, tol_ang, tol_tang);
    let _ = (pl1, pl2);
    if !inter.is_done() {
        return false;
    }
    *same = false;
    let typint = inter.type_inter();
    if typint == AnaResultType::Same {
        // coincident faces
        *same = true;
    } else if typint != AnaResultType::Empty {
        // an intersection line
        let linsol = inter.line(1);
        let discri = linsol.direction.dot(quad2.normale(linsol.origin).cross(quad1.normale(linsol.origin)));
        let (t1, t2) = if discri > 0.0 {
            (TypeTrans::Out, TypeTrans::In)
        } else {
            (TypeTrans::In, TypeTrans::Out)
        };
        slin.push(gline_line(linsol, false, t1, t2));
    }
    true
}

/// OCCT IntPCy (L3157-3345) — Plane/Cylinder (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_pcy(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_ang: f64,
    tol_tang: f64,
    reversed: bool,
    empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    h: f64,
) -> bool {
    let (pl, cy) = if !reversed {
        (quad1, quad2)
    } else {
        (quad2, quad1)
    };
    let _ = (pl, cy);

    let mut inter = QuadQuadGeo::new();
    inter.perform_plane_cylinder(pl, cy, tol_ang, tol_tang, h);
    if !inter.is_done() {
        return false;
    }
    let typint = inter.type_inter();
    let nb_sol = inter.nb_solutions();
    *empty = false;

    match typint {
        AnaResultType::Empty => {
            *empty = true;
        }
        AnaResultType::Line => {
            let linsol = inter.line(1);
            let orig = linsol.origin;
            if nb_sol == 1 {
                // tangency line
                let test_curvature = orig - cy.axis_loc();
                let (normp, normcyl) = if !reversed {
                    (quad1.normale(orig), quad2.normale(orig))
                } else {
                    (quad2.normale(orig), quad1.normale(orig))
                };
                let (situcyl, situp) = if normp.dot(test_curvature) > 0.0 {
                    let s = if normp.dot(normcyl) > 0.0 { Situation::Inside } else { Situation::Outside };
                    (Situation::Outside, s)
                } else {
                    let s = if normp.dot(normcyl) > 0.0 { Situation::Outside } else { Situation::Inside };
                    (Situation::Inside, s)
                };
                let glig = if !reversed {
                    gline_line_touch(linsol, true, situp, situcyl)
                } else {
                    gline_line_touch(linsol, true, situcyl, situp)
                };
                slin.push(glig);
            } else {
                // two lines: determine each transition
                if linsol.direction.dot(quad2.normale(orig).cross(quad1.normale(orig))) > 0.0 {
                    slin.push(gline_line(linsol, false, TypeTrans::Out, TypeTrans::In));
                } else {
                    slin.push(gline_line(linsol, false, TypeTrans::In, TypeTrans::Out));
                }
                let linsol2 = inter.line(2);
                let orig2 = linsol2.origin;
                if linsol2.direction.dot(quad2.normale(orig2).cross(quad1.normale(orig2))) > 0.0 {
                    slin.push(gline_line(linsol2, false, TypeTrans::Out, TypeTrans::In));
                } else {
                    slin.push(gline_line(linsol2, false, TypeTrans::In, TypeTrans::Out));
                }
            }
        }
        AnaResultType::Circle => {
            let cirsol = inter.circle();
            let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
            if tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref))) > 0.0 {
                slin.push(gline_circle(cirsol, false, TypeTrans::Out, TypeTrans::In));
            } else {
                slin.push(gline_circle(cirsol, false, TypeTrans::In, TypeTrans::Out));
            }
        }
        AnaResultType::Ellipse => {
            let elipsol = inter.ellipse();
            let (ptref, tgt) = ecl::ellipse_d1(&elipsol, 0.0);
            if tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref))) > 0.0 {
                slin.push(gline_ellipse(elipsol, false, TypeTrans::Out, TypeTrans::In));
            } else {
                slin.push(gline_ellipse(elipsol, false, TypeTrans::In, TypeTrans::Out));
            }
        }
        _ => {
            return false; // should not happen
        }
    }
    true
}

/// OCCT IntPSp (L3352-3439) — Plane/Sphere (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_psp(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_ang: f64,
    tol_tang: f64,
    reversed: bool,
    empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let (pl, sp) = if !reversed { (quad1, quad2) } else { (quad2, quad1) };
    let _ = (pl, sp);

    let mut inter = QuadQuadGeo::new();
    inter.perform_plane_sphere(pl, sp);
    if !inter.is_done() {
        return false;
    }
    let typint = inter.type_inter();
    *empty = false;

    match typint {
        AnaResultType::Empty => {
            *empty = true;
        }
        AnaResultType::Point => {
            let psol = inter.point(1);
            let (u1, v1) = quad1.parameters(psol);
            let (u2, v2) = quad2.parameters(psol);
            let mut ptsol = PatchPoint {
                pnt: PntOn2S { p: psol, u1, v1, u2, v2 },
                param_on_line: 0.0,
                tolerance: tol_tang,
                multiple: false,
                on_dom_s1: true,
                on_dom_s2: true,
                arc_on_s1: None,
                arc_on_s2: None,
                param_on_arc1: 0.0,
                param_on_arc2: 0.0,
            };
            let _ = &mut ptsol;
            spnt.push(IntPatchPoint {
                p1: psol,
                p2: psol,
                u1,
                v1,
                u2,
                v2,
                tolerance: tol_tang,
            });
        }
        AnaResultType::Circle => {
            let cirsol = inter.circle();
            let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
            if tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref))) > 0.0 {
                slin.push(gline_circle(cirsol, false, TypeTrans::Out, TypeTrans::In));
            } else {
                slin.push(gline_circle(cirsol, false, TypeTrans::In, TypeTrans::Out));
            }
        }
        _ => {
            return false;
        }
    }
    let _ = tol_ang;
    true
}

/// OCCT IntPTo (L3790-3857) — Plane/Torus (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_pto(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_tang: f64,
    b_reversed: bool,
    b_empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
) -> bool {
    let a_pln = if b_reversed { quad2.plane() } else { quad1.plane() };
    let a_torus = if b_reversed { quad1.torus() } else { quad2.torus() };

    let mut inter = QuadQuadGeo::new();
    inter.perform_plane_torus(quad1, quad2, tol_tang);
    let mut b_ret = inter.is_done();
    if !b_ret {
        return b_ret;
    }
    let typint = inter.type_inter();
    let nb_sol = inter.nb_solutions();
    *b_empty = false;

    match typint {
        AnaResultType::Empty => {
            *b_empty = true;
        }
        AnaResultType::Circle => {
            for i in 1..=nb_sol {
                let mut a_c = inter.circle_n(i as i32);
                // OCCT: if the plane is not normal to the torus axis, adjust the
                // circle to the torus seam.
                let pln_norm = a_pln.normal.normalize_or_zero();
                let tor_axis = a_torus.axis.normalize_or_zero();
                if pln_norm.cross(tor_axis).length() > rcad_kernel::precision::ANGULAR {
                    let quad_dz = if b_reversed { quad1.z_dir() } else { quad2.z_dir() };
                    let quad_dx = if b_reversed { quad1.x_dir() } else { quad2.x_dir() };
                    ecl::adjust_to_seam_quadric(&mut a_c, a_torus.center, quad_dz, quad_dx);
                }
                let (ptref, tgt) = ecl::circle_d1(&a_c, 0.0);
                if tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref))) > 0.0 {
                    slin.push(gline_circle(a_c, false, TypeTrans::Out, TypeTrans::In));
                } else {
                    slin.push(gline_circle(a_c, false, TypeTrans::In, TypeTrans::Out));
                }
            }
        }
        _ => {
            b_ret = false;
        }
    }
    b_ret
}

/// OCCT IntPCo (L3446-3786) — Plane/Cone (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_pco(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_ang: f64,
    tol_tang: f64,
    reversed: bool,
    empty: &mut bool,
    multpoint: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let (pl, co) = if !reversed { (quad1, quad2) } else { (quad2, quad1) };
    let apex = co.cone().apex_point();
    let _ = (pl, co);

    let mut inter = QuadQuadGeo::new();
    inter.perform_plane_cone(pl, co, tol_ang, tol_tang);
    if !inter.is_done() {
        return false;
    }
    let typint = inter.type_inter();
    let nb_sol = inter.nb_solutions();
    *empty = false;

    match typint {
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
                tolerance: tol_tang,
            });
        }
        AnaResultType::Line => {
            let mut linsol = inter.line(1);
            let co_dir = co.cone().axis_dir();
            if linsol.direction.dot(co_dir) < 0.0 {
                linsol.direction = -linsol.direction;
            }
            let para = ecl::line_parameter(&linsol, apex);
            let ptbid = ecl::line_value(&linsol, para + 5.0);
            let (u1, v1) = quad1.parameters(apex);
            let (u2, v2) = quad2.parameters(apex);

            if nb_sol == 1 {
                // tangency line
                let mut ptsol = make_point(apex, tol_tang, false, quad1, quad2);
                ptsol.param_on_line = para;
                let ptbid2 = apex + 5.0 * co_dir;
                let test_curvature = ptbid - ptbid2;
                let (normp, normco) = if !reversed {
                    (quad1.normale(ptbid), quad2.normale(ptbid))
                } else {
                    (quad2.normale(ptbid), quad1.normale(ptbid))
                };
                let (situco, situco_otherside, situp, situp_otherside);
                if normp.dot(test_curvature) > 0.0 {
                    situco = Situation::Outside;
                    situco_otherside = Situation::Inside;
                    if normp.dot(normco) > 0.0 {
                        situp = Situation::Inside;
                        situp_otherside = Situation::Outside;
                    } else {
                        situp = Situation::Outside;
                        situp_otherside = Situation::Inside;
                    }
                } else {
                    situco = Situation::Inside;
                    situco_otherside = Situation::Outside;
                    if normp.dot(normco) > 0.0 {
                        situp = Situation::Outside;
                        situp_otherside = Situation::Inside;
                    } else {
                        situp = Situation::Inside;
                        situp_otherside = Situation::Outside;
                    }
                }
                // Apex -> Cone.Direction
                let mut glig = if !reversed {
                    gline_line_touch(linsol, true, situp, situco)
                } else {
                    gline_line_touch(linsol, true, situco, situp)
                };
                add_vertex_int(&mut glig, ptsol.clone());
                glig.first_point = Some(glig.vertices.len());
                slin.push(glig);
                // -Cone.Direction <--- Apex
                linsol.direction = -linsol.direction;
                let mut glig2 = if !reversed {
                    gline_line_touch(linsol, true, situp_otherside, situco_otherside)
                } else {
                    gline_line_touch(linsol, true, situco_otherside, situp_otherside)
                };
                add_vertex_int(&mut glig2, ptsol.clone());
                glig2.first_point = Some(glig2.vertices.len());
                slin.push(glig2);
            } else {
                // Two lines: determine the transitions of each.  Each line is
                // oriented along the cone axis; the transitions of the two
                // lines are inverse of each other, so compute only the first.
                // OCCT L3605: if (dir.DotCross(N2, N1) > 0.) Out/In else In/Out.
                let (mut trans1, mut trans2) = if linsol
                    .direction
                    .dot(quad2.normale(ptbid).cross(quad1.normale(ptbid)))
                    > 0.0
                {
                    (TypeTrans::Out, TypeTrans::In)
                } else {
                    (TypeTrans::In, TypeTrans::Out)
                };
                *multpoint = true;
                // Ligne 1
                let mut ptsol = make_point(apex, tol_tang, false, quad1, quad2);
                ptsol.param_on_line = para;
                ptsol.multiple = true;
                let mut l1 = gline_line(linsol, false, trans1, trans2);
                add_vertex_int(&mut l1, ptsol.clone());
                l1.first_point = Some(l1.vertices.len());
                slin.push(l1);
                // Other side: transitions stay the same.
                linsol.direction = -linsol.direction;
                let p1r = ecl::line_parameter(&linsol, apex);
                ptsol.param_on_line = p1r;
                let mut l1r = gline_line(linsol, false, trans1, trans2);
                add_vertex_int(&mut l1r, ptsol.clone());
                l1r.first_point = Some(l1r.vertices.len());
                slin.push(l1r);
                // Ligne 2
                let mut linsol2 = inter.line(2);
                if linsol2.direction.dot(co_dir) < 0.0 {
                    linsol2.direction = -linsol2.direction;
                }
                let p2 = ecl::line_parameter(&linsol2, apex);
                let ptbid2 = ecl::line_value(&linsol2, p2 + 5.0);
                let (trans1b, trans2b) = if linsol2
                    .direction
                    .dot(quad2.normale(ptbid2).cross(quad1.normale(ptbid2)))
                    > 0.0
                {
                    (TypeTrans::Out, TypeTrans::In)
                } else {
                    (TypeTrans::In, TypeTrans::Out)
                };
                trans1 = trans1b;
                trans2 = trans2b;
                ptsol.param_on_line = p2;
                let mut l2 = gline_line(linsol2, false, trans1, trans2);
                add_vertex_int(&mut l2, ptsol.clone());
                l2.first_point = Some(l2.vertices.len());
                slin.push(l2);
                // Other side.
                linsol2.direction = -linsol2.direction;
                let p2r = ecl::line_parameter(&linsol2, apex);
                ptsol.param_on_line = p2r;
                let mut l2r = gline_line(linsol2, false, trans1, trans2);
                add_vertex_int(&mut l2r, ptsol.clone());
                l2r.first_point = Some(l2r.vertices.len());
                slin.push(l2r);
                let _ = (u1, v1, u2, v2);
            }
        }
        AnaResultType::Circle => {
            let mut cirsol = inter.circle();
            let co_frame_z = co.cone().axis_dir();
            let co_frame_x = rcad_kernel::geom::any_perpendicular(co_frame_z).normalize_or_zero();
            ecl::adjust_to_seam_quadric(&mut cirsol, co.cone().apex, co_frame_z, co_frame_x);
            let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
            if tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref))) > 0.0 {
                slin.push(gline_circle(cirsol, false, TypeTrans::Out, TypeTrans::In));
            } else {
                slin.push(gline_circle(cirsol, false, TypeTrans::In, TypeTrans::Out));
            }
        }
        AnaResultType::Ellipse => {
            let elipsol = inter.ellipse();
            let (ptref, tgt) = ecl::ellipse_d1(&elipsol, 0.0);
            if tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref))) > 0.0 {
                slin.push(gline_ellipse(elipsol, false, TypeTrans::Out, TypeTrans::In));
            } else {
                slin.push(gline_ellipse(elipsol, false, TypeTrans::In, TypeTrans::Out));
            }
        }
        AnaResultType::Parabola => {
            let parabsol = inter.parabola();
            // OCCT: Tgtorig = Parab.YAxis().Direction().  rcad Parabola3 tangent
            // at the vertex (t=0) is normal x axis_dir (the local Y axis).
            let tgt_orig = parabsol.normal.normalize_or_zero()
                .cross(parabsol.axis_dir.normalize_or_zero())
                .normalize_or_zero();
            let ptran = tgt_orig
                .dot(quad2.normale(parabsol.vertex).cross(quad1.normale(parabsol.vertex)));
            let (t1, t2) = transition_from_scalar(ptran);
            slin.push(gline_parabola(parabsol, false, t1, t2));
        }
        AnaResultType::Hyperbola => {
            for i in 1..=2 {
                let hyprsol = inter.hyperbola_n(i as i32);
                // OCCT: tophypr = ElCLib::Value(MajorRadius, XAxis) = the vertex
                // on the transverse axis; Tgttop = YAxis().Direction().
                let tophypr = hyprsol.center + hyprsol.semi_major * hyprsol.major_dir.normalize_or_zero();
                let major = hyprsol.major_dir.normalize_or_zero();
                let normal = hyprsol.normal.normalize_or_zero();
                let tgttop = normal.cross(major).normalize_or_zero();
                let qwe = tgttop.dot(quad2.normale(tophypr).cross(quad1.normale(tophypr)));
                let (t1, t2) = transition_from_scalar(qwe);
                slin.push(gline_hyperbola(hyprsol, false, t1, t2));
            }
        }
        _ => {
            return false;
        }
    }
    true
}

/// OCCT: qwe > 1e-8 -> Out/In, < -1e-8 -> In/Out, else Undecided/Undecided.
fn transition_from_scalar(qwe: f64) -> (TypeTrans, TypeTrans) {
    if qwe > 1.0e-8 {
        (TypeTrans::Out, TypeTrans::In)
    } else if qwe < -1.0e-8 {
        (TypeTrans::In, TypeTrans::Out)
    } else {
        (TypeTrans::Undecided, TypeTrans::Undecided)
    }
}

/// Add an IntPatch_Point (as PatchPoint vertex) to a GLine/ALine.
fn add_vertex_int(line: &mut IntPatchLine, v: super::special_points::PatchPoint) {
    if line.line_type == IntPatchIType::Analytic {
        if let Some(ac) = line.a_curve.as_mut() {
            ac.vertices.push(v);
        }
    } else {
        line.vertices.push(super::IntPatchVertex {
            param_on_line: v.param_on_line,
            p3d: v.pnt.p,
            u1: v.pnt.u1,
            v1: v.pnt.v1,
            u2: v.pnt.u2,
            v2: v.pnt.v2,
            tolerance: v.tolerance,
            multiple: v.multiple,
            on_dom_s1: v.on_dom_s1,
            on_dom_s2: v.on_dom_s2,
            arc_on_s1: v.arc_on_s1,
            arc_on_s2: v.arc_on_s2,
            param_on_arc1: v.param_on_arc1,
            param_on_arc2: v.param_on_arc2,
        });
    }
}
