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
use crate::geomalgo::int_surf::quadric::Quadric;
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
    // OCCT Geom_Parabola FirstParameter/LastParameter = -/+Precision::Infinite().
    let dom = [
        -rcad_kernel::precision::INFINITE_VALUE,
        rcad_kernel::precision::INFINITE_VALUE,
    ];
    let mut line = IntPatchLine::analytic(IntPatchIType::Parabola, Curve3::Parabola(p), dom);
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// OCCT IntPatch_GLine(Hypr, Tang, Trans1, Trans2).
fn gline_hyperbola(h: Hyperbola3, tang: bool, t1: TypeTrans, t2: TypeTrans) -> IntPatchLine {
    // OCCT Geom_Hyperbola FirstParameter/LastParameter = -/+Precision::Infinite().
    let dom = [
        -rcad_kernel::precision::INFINITE_VALUE,
        rcad_kernel::precision::INFINITE_VALUE,
    ];
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
        is_vertex_on_s1: false,
        is_vertex_on_s2: false,
        transition_line_arc1: super::transitions::TypeTrans::Undecided,
        transition_line_arc2: super::transitions::TypeTrans::Undecided,
        transition_on_s1: super::transitions::TypeTrans::Undecided,
        transition_on_s2: super::transitions::TypeTrans::Undecided,
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
            // OCCT L3298-3300: cirsol = inter.Circle(1); AdjustToSeam(Cy, cirsol).
            let mut cirsol = inter.circle();
            ecl::adjust_to_seam_quadric(&mut cirsol, cy.axis_loc(), cy.z_dir(), cy.x_dir());
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
                is_vertex_on_s1: false,
                is_vertex_on_s2: false,
                transition_line_arc1: super::transitions::TypeTrans::Undecided,
                transition_line_arc2: super::transitions::TypeTrans::Undecided,
                transition_on_s1: super::transitions::TypeTrans::Undecided,
                transition_on_s2: super::transitions::TypeTrans::Undecided,
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
            let co_frame_x = co.cone().ref_dir.normalize_or_zero();
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

/// Reconstruct a Surface3 from a Quadric.
fn quadric_to_surface3(quad: &Quadric) -> rcad_kernel::geom::Surface3 {
    use crate::geomalgo::int_surf::quadric::QuadricType;
    match quad.type_quadric() {
        QuadricType::Plane => rcad_kernel::geom::Surface3::Plane(quad.plane()),
        QuadricType::Cylinder => rcad_kernel::geom::Surface3::Cylinder(quad.cylinder()),
        QuadricType::Sphere => rcad_kernel::geom::Surface3::Sphere(quad.sphere()),
        QuadricType::Cone => rcad_kernel::geom::Surface3::Cone(quad.cone()),
        QuadricType::Torus => rcad_kernel::geom::Surface3::Torus(quad.torus()),
        QuadricType::Other => rcad_kernel::geom::Surface3::Plane(quad.plane()),
    }
}

/// Build an IntPatch_ALine (IntPatchLine with a_curve + transitions) from an
/// IntAna_Curve.
fn aline(
    curvsol: super::int_quad_quad::IntAnaCurve,
    tang: bool,
    t1: TypeTrans,
    t2: TypeTrans,
) -> IntPatchLine {
    let d = curvsol.domain();
    let mut line = IntPatchLine::analytic(
        IntPatchIType::Analytic,
        Curve3::Line(Line3 {
            origin: DVec3::ZERO,
            direction: DVec3::X,
        }),
        [d[0], d[1]],
    );
    line.a_curve = Some(curvsol);
    line.trans1 = Some(Transition::new_in_out(tang, t1));
    line.trans2 = Some(Transition::new_in_out(tang, t2));
    line
}

/// Shared NoGeometricSolution fallback (IntCySp L8263-8362, IntCyCo L8465-8593,
/// IntCoCo L9022-9139, IntCoSp L9349-9461): run IntAna_IntQuadQuad; append its
/// points to spnt and its curves as ALine (with transitions + ProcessBounds).
#[allow(clippy::too_many_arguments)]
fn int_quad_quad_fallback(
    quad1: &Quadric,
    quad2: &Quadric,
    explicit: &rcad_kernel::geom::Surface3,
    other: &Quadric,
    tol: f64,
    empty: &mut bool,
    multpoint: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
    use_explore_curve: bool,
    qwe_tol: f64,
    check_tangent_mag: bool,
) -> bool {
    let other_surf = quadric_to_surface3(other);
    let Some(other_quad) = super::int_quad_quad::IntAnaQuadric::from_surface3(&other_surf) else {
        return false;
    };
    let mut anaint = super::int_quad_quad::IntQuadQuad::new();
    match explicit {
        rcad_kernel::geom::Surface3::Cylinder(c) => {
            anaint.perform_cylinder(c, &other_quad);
        }
        rcad_kernel::geom::Surface3::Cone(c) => {
            anaint.perform_cone(c, &other_quad);
        }
        _ => return false,
    }
    if !anaint.is_done() {
        return false;
    }
    if anaint.nb_points() == 0 && anaint.nb_curves() == 0 {
        *empty = true;
        return true;
    }
    // Points.
    let nb_pnt = anaint.nb_points();
    for i in 1..=nb_pnt {
        let psol = anaint.point(i);
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
    // Curves.  IntCyCo uses ExploreCurve to split each curve at the cone apex.
    let nb_curv = anaint.nb_curves();
    for i in 0..nb_curv {
        let base_curve = anaint.curve(i).unwrap().clone();
        let mut curves: Vec<super::int_quad_quad::IntAnaCurve> = Vec::new();
        if use_explore_curve {
            // OCCT IntCyCo NoGeometricSolution (L8512): ExploreCurve(Co, aC, ...)
            // splits at the CONE apex.  The cone is the other quadric here (the
            // explicit surface may be the cylinder).
            let cone = if let rcad_kernel::geom::Surface3::Cone(c) = other_surf {
                c
            } else if let rcad_kernel::geom::Surface3::Cone(c) = explicit {
                *c
            } else {
                return false;
            };
            explore_curve(&cone, &base_curve, 10.0 * tol, &mut curves);
        } else {
            curves.push(base_curve);
        }
        for curvsol in curves {
            let d = curvsol.domain();
            let (first, last) = (d[0], d[1]);
            let firstp = !curvsol.is_first_open();
            let lastp = !curvsol.is_last_open();
            let ptf = if firstp { curvsol.value(first).unwrap_or(DVec3::ZERO) } else { DVec3::ZERO };
            let ptl = if lastp { curvsol.value(last).unwrap_or(DVec3::ZERO) } else { DVec3::ZERO };
            // Find a parameter where the tangent is valid.
            let mut para = last;
            let mut kount = 1;
            let mut tgfound = false;
            let mut ptvalid = DVec3::ZERO;
            let mut tgvalid = DVec3::ZERO;
            while !tgfound {
                para = (1.123 * first + para) / 2.123;
                match curvsol.d1u(para) {
                    Some((pv, tv)) => {
                        ptvalid = pv;
                        tgvalid = tv;
                        tgfound = true;
                        // OCCT IntCoCo NoGeometricSolution (L9080-9084): a
                        // near-zero tangent means the normals are meaningless
                        // there — retry with the next parameter.
                        if check_tangent_mag && tgvalid.length_squared() < 1.0e-14 {
                            tgfound = false;
                        }
                    }
                    None => {
                        tgfound = false;
                    }
                }
                if !tgfound {
                    kount += 1;
                    tgfound = kount > 5;
                }
            }
            let (trans1, trans2);
            let mut kept = false;
            if kount <= 5 {
                let qwe = tgvalid.dot(quad2.normale(ptvalid).cross(quad1.normale(ptvalid)));
                (trans1, trans2) = if qwe > qwe_tol {
                    (TypeTrans::Out, TypeTrans::In)
                } else if qwe < -qwe_tol {
                    (TypeTrans::In, TypeTrans::Out)
                } else {
                    (TypeTrans::Undecided, TypeTrans::Undecided)
                };
                kept = true;
            } else {
                ptvalid = curvsol.value(para).unwrap_or(DVec3::ZERO);
                (trans1, trans2) = (TypeTrans::Undecided, TypeTrans::Undecided);
                kept = true;
            }
            if kept {
                let mut alig = aline(curvsol.clone(), false, trans1, trans2);
                let mut n_firstp = !firstp;
                let mut n_lastp = !lastp;
                if let Some(ac) = alig.a_curve.as_mut() {
                    super::imp_imp_intersection::process_bounds(
                        ac,
                        slin,
                        quad1,
                        quad2,
                        &mut n_firstp,
                        ptf,
                        first,
                        &mut n_lastp,
                        ptl,
                        last,
                        multpoint,
                        tol,
                    );
                }
                slin.push(alig);
            }
        }
    }
    true
}

/// OCCT IntSpSp (L9473-9554) — Sphere/Sphere.
#[allow(clippy::too_many_arguments)]
pub fn int_spsp(
    quad1: &Quadric,
    quad2: &Quadric,
    tol: f64,
    empty: &mut bool,
    same: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let mut inter = QuadQuadGeo::new();
    inter.perform_sphere_sphere(quad1, quad2, tol);
    if !inter.is_done() {
        return false;
    }
    let typint = inter.type_inter();
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
        AnaResultType::Circle => {
            let cirsol = inter.circle();
            let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
            let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
            let (t1, t2) = transition_from_scalar(qwe);
            slin.push(gline_circle(cirsol, false, t1, t2));
        }
        _ => {
            return false;
        }
    }
    true
}

/// OCCT TreatResultTorus (L9648-9714) — shared torus result handling.
fn treat_result_torus(
    quad1: &Quadric,
    quad2: &Quadric,
    inter: &QuadQuadGeo,
    b_empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
) -> bool {
    let b_ret = inter.is_done();
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
                // OCCT: AdjustToSeam(Torus1, aC) when both quadrics are tori.
                if quad1.type_quadric() == quad2.type_quadric() {
                    ecl::adjust_to_seam_quadric(&mut a_c, quad1.axis_loc(), quad1.z_dir(), quad1.x_dir());
                }
                let (ptref, tgt) = ecl::circle_d1(&a_c, 0.0);
                let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
                let (t1, t2) = transition_from_scalar(qwe);
                slin.push(gline_circle(a_c, false, t1, t2));
            }
        }
        _ => {
            return false;
        }
    }
    b_ret
}

/// OCCT IntCyTo (L9564-9578) — Cylinder/Torus (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_cyto(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_tang: f64,
    b_reversed: bool,
    b_empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
) -> bool {
    let (a_cyl, a_torus) = if b_reversed { (quad2, quad1) } else { (quad1, quad2) };
    let _ = (a_cyl, a_torus);
    let mut inter = QuadQuadGeo::new();
    inter.perform_cylinder_torus(quad1, quad2, tol_tang);
    treat_result_torus(quad1, quad2, &inter, b_empty, slin)
}

/// OCCT IntCoTo (L9582-9596) — Cone/Torus (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_coto(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_tang: f64,
    b_reversed: bool,
    b_empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
) -> bool {
    let _ = (b_reversed, quad1, quad2);
    let mut inter = QuadQuadGeo::new();
    inter.perform_cone_torus(quad1, quad2, tol_tang);
    treat_result_torus(quad1, quad2, &inter, b_empty, slin)
}

/// OCCT IntSpTo (L9600-9614) — Sphere/Torus (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_spto(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_tang: f64,
    b_reversed: bool,
    b_empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
) -> bool {
    let _ = (b_reversed, quad1, quad2);
    let mut inter = QuadQuadGeo::new();
    inter.perform_sphere_torus(quad1, quad2, tol_tang);
    treat_result_torus(quad1, quad2, &inter, b_empty, slin)
}

/// OCCT IntToTo (L9618-9644) — Torus/Torus.
#[allow(clippy::too_many_arguments)]
pub fn int_toto(
    quad1: &Quadric,
    quad2: &Quadric,
    tol_tang: f64,
    b_same_surf: &mut bool,
    b_empty: &mut bool,
    slin: &mut Vec<IntPatchLine>,
) -> bool {
    let mut inter = QuadQuadGeo::new();
    inter.perform_torus_torus(quad1, quad2, tol_tang);
    let mut b_ret = inter.is_done();
    if b_ret {
        if inter.type_inter() == AnaResultType::Same {
            *b_empty = false;
            *b_same_surf = true;
        } else {
            b_ret = treat_result_torus(quad1, quad2, &inter, b_empty, slin);
        }
    }
    b_ret
}

/// OCCT IntCySp (L8107-8370) — Cylinder/Sphere (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_cysp(
    quad1: &Quadric,
    quad2: &Quadric,
    tol: f64,
    reversed: bool,
    empty: &mut bool,
    multpoint: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let (cy, sp) = if !reversed { (quad1, quad2) } else { (quad2, quad1) };

    let mut inter = QuadQuadGeo::new();
    inter.perform_cylinder_sphere(cy, sp, tol);
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
        AnaResultType::Circle => {
            let cirsol = inter.circle();
            let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
            if nb_sol == 1 {
                // Tangent circle: use Situation transitions.
                let test_curvature = ptref - sp.sphere().center;
                let (normcyl, normsp) = if !reversed {
                    (quad1.normale(ptref), quad2.normale(ptref))
                } else {
                    (quad2.normale(ptref), quad1.normale(ptref))
                };
                let (situcyl, situsp) = if normcyl.dot(test_curvature) > 0.0 {
                    let s = if normsp.dot(normcyl) > 0.0 { Situation::Inside } else { Situation::Outside };
                    (s, Situation::Outside)
                } else {
                    let s = if normsp.dot(normcyl) > 0.0 { Situation::Outside } else { Situation::Inside };
                    (s, Situation::Inside)
                };
                let glig = if !reversed {
                    gline_circle_touch(cirsol, true, situcyl, situsp)
                } else {
                    gline_circle_touch(cirsol, true, situsp, situcyl)
                };
                slin.push(glig);
            } else {
                // Two circles: TypeTrans transitions.
                let qwe1 = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
                let (t1, t2) = if qwe1 > 0.0 {
                    (TypeTrans::Out, TypeTrans::In)
                } else {
                    (TypeTrans::In, TypeTrans::Out)
                };
                slin.push(gline_circle(cirsol, false, t1, t2));
                let cirsol2 = inter.circle_n(2);
                let (ptref2, tgt2) = ecl::circle_d1(&cirsol2, 0.0);
                let qwe2 = tgt2.dot(quad2.normale(ptref2).cross(quad1.normale(ptref2)));
                let (t1, t2) = transition_from_scalar_1e7(qwe2);
                slin.push(gline_circle(cirsol2, false, t1, t2));
            }
        }
        AnaResultType::NoGeometricSolution => {
            // OCCT: IntAna_IntQuadQuad(Cy, Sp).  rcad: IntQuadQuad with the
            // cylinder as the explicit surface, the other as the quadric.
            return int_quad_quad_fallback(
                quad1,
                quad2,
                &rcad_kernel::geom::Surface3::Cylinder(cy.cylinder()),
                sp,
                tol,
                empty,
                multpoint,
                slin,
                spnt,
                false,
                1.0e-8,
                false,
            );
        }
        _ => {
            return false;
        }
    }
    true
}

/// OCCT IntCoSp (L9179-9470) — Cone/Sphere (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_cosp(
    quad1: &Quadric,
    quad2: &Quadric,
    tol: f64,
    reversed: bool,
    empty: &mut bool,
    multpoint: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let (co, sp) = if !reversed { (quad1, quad2) } else { (quad2, quad1) };
    let _ = (co, sp);

    let mut inter = QuadQuadGeo::new();
    inter.perform_sphere_cone(sp, co, tol);
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
        AnaResultType::Point => {
            let apex = co.cone().apex_point();
            let paramapex = ecl::line_parameter_of_axis(co.cone().apex, co.cone().axis_dir(), apex);
            for i in 1..=nb_sol {
                let ptcontact = inter.point(i as i32);
                let param = ecl::line_parameter_of_axis(co.cone().apex, co.cone().axis_dir(), ptcontact);
                let (u1, v1) = quad1.parameters(ptcontact);
                let (u2, v2) = quad2.parameters(ptcontact);
                if apex.distance(ptcontact) <= tol {
                    spnt.push(IntPatchPoint {
                        p1: ptcontact,
                        p2: ptcontact,
                        u1,
                        v1,
                        u2,
                        v2,
                        tolerance: tol,
                    });
                } else if param >= paramapex {
                    spnt.push(IntPatchPoint {
                        p1: ptcontact,
                        p2: ptcontact,
                        u1,
                        v1,
                        u2,
                        v2,
                        tolerance: tol,
                    });
                }
            }
        }
        AnaResultType::Circle => {
            for i in 1..=nb_sol {
                let cirsol = inter.circle_n(i as i32);
                let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
                let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
                let (t1, t2) = transition_from_scalar(qwe);
                slin.push(gline_circle(cirsol, false, t1, t2));
            }
        }
        AnaResultType::PointAndCircle => {
            let apex = co.cone().apex_point();
            let paramapex = ecl::line_parameter_of_axis(co.cone().apex, co.cone().axis_dir(), apex);
            // The point is necessarily the apex.
            let (u1, v1) = quad1.parameters(apex);
            let (u2, v2) = quad2.parameters(apex);
            spnt.push(IntPatchPoint {
                p1: apex,
                p2: apex,
                u1,
                v1,
                u2,
                v2,
                tolerance: tol,
            });
            let mut cirsol = inter.circle();
            let param = ecl::line_parameter_of_axis(co.cone().apex, co.cone().axis_dir(), cirsol.center);
            let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
            let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
            let (t1, t2) = if param >= paramapex {
                if qwe > rcad_kernel::precision::PCONFUSION {
                    (TypeTrans::Out, TypeTrans::In)
                } else if qwe < -rcad_kernel::precision::PCONFUSION {
                    (TypeTrans::In, TypeTrans::Out)
                } else {
                    (TypeTrans::Undecided, TypeTrans::Undecided)
                }
            } else if qwe < -rcad_kernel::precision::PCONFUSION {
                (TypeTrans::Out, TypeTrans::In)
            } else if qwe > rcad_kernel::precision::PCONFUSION {
                (TypeTrans::In, TypeTrans::Out)
            } else {
                (TypeTrans::Undecided, TypeTrans::Undecided)
            };
            slin.push(gline_circle(cirsol, false, t1, t2));
        }
        AnaResultType::NoGeometricSolution => {
            return int_quad_quad_fallback(
                quad1,
                quad2,
                &rcad_kernel::geom::Surface3::Cone(co.cone()),
                sp,
                tol,
                empty,
                multpoint,
                slin,
                spnt,
                false,
                1.0e-9,
                false,
            );
        }
        _ => {
            return false;
        }
    }
    true
}

/// OCCT IntCoCo (L8679-9175) — Cone/Cone.
#[allow(clippy::too_many_arguments)]
pub fn int_coco(
    quad1: &Quadric,
    quad2: &Quadric,
    tol: f64,
    empty: &mut bool,
    same: &mut bool,
    multpoint: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let co1 = quad1.cone();
    let co2 = quad2.cone();
    let apex1 = co1.apex_point();
    let apex2 = co2.apex_point();

    let mut inter = QuadQuadGeo::new();
    inter.perform_cone_cone(quad1, quad2, tol);
    if !inter.is_done() {
        return false;
    }
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
        AnaResultType::Line => {
            let (mut trans1, mut trans2) = (TypeTrans::In, TypeTrans::Out);
            if nb_sol == 1 {
                // Tangency line.
                let linsol = inter.line(1);
                let para = ecl::line_parameter(&linsol, apex1);
                let ptbid = ecl::line_value(&linsol, para + 5.0);
                let (u1, v1) = quad1.parameters(apex1);
                let (u2, v2) = quad2.parameters(apex1);
                let mut a_ptsol = make_point(apex1, tol, false, quad1, quad2);
                a_ptsol.param_on_line = para;
                let (mut u1, mut v1, mut u2, mut v2) = (u1, v1, u2, v2);
                let _ = (&mut u1, &mut v1, &mut u2, &mut v2);
                let norm_c1 = quad1.normale(ptbid);
                let norm_c2 = quad2.normale(ptbid);
                let a_dot = norm_c1.dot(norm_c2);
                let (mut situ_c1, mut situ_c2) = if a_dot < 0.0 {
                    (Situation::Outside, Situation::Outside)
                } else {
                    // Use distance from ptbid to each cone axis.
                    let a_l_ax1 = Line3 { origin: apex1, direction: co1.axis_dir() };
                    let a_l_ax2 = Line3 { origin: apex2, direction: co2.axis_dir() };
                    let a_r1 = point_line_distance(ptbid, &a_l_ax1);
                    let a_r2 = point_line_distance(ptbid, &a_l_ax2);
                    if a_r1 > a_r2 {
                        (Situation::Outside, Situation::Inside)
                    } else {
                        (Situation::Inside, Situation::Outside)
                    }
                };
                // 1
                let mut glig = gline_line_touch(linsol, true, situ_c1, situ_c2);
                add_vertex_int(&mut glig, a_ptsol.clone());
                glig.first_point = Some(glig.vertices.len());
                slin.push(glig);
                // 2
                let mut linsol_r = linsol;
                linsol_r.direction = -linsol_r.direction;
                let para_r = ecl::line_parameter(&linsol_r, apex1);
                a_ptsol.param_on_line = para_r;
                // OCCT L8778-8785: the reversed tangency line carries the
                // situations swapped (situC2, situC1).
                let mut glig2 = gline_line_touch(linsol_r, true, situ_c2, situ_c1);
                add_vertex_int(&mut glig2, a_ptsol.clone());
                glig2.first_point = Some(glig2.vertices.len());
                slin.push(glig2);
            } else if nb_sol == 2 {
                for i in 1..=2 {
                    let linsol = inter.line(i as i32);
                    let para = ecl::line_parameter(&linsol, apex1);
                    let ptbid = ecl::line_value(&linsol, para + 5.0);
                    let (u1, v1) = quad1.parameters(apex1);
                    let (u2, v2) = quad2.parameters(apex1);
                    let _ = (u1, v1, u2, v2);
                    trans1 = TypeTrans::In;
                    trans2 = TypeTrans::Out;
                    if linsol.direction.dot(quad2.normale(ptbid).cross(quad1.normale(ptbid))) > 0.0 {
                        trans1 = TypeTrans::Out;
                        trans2 = TypeTrans::In;
                    }
                    *multpoint = true;
                    let mut a_ptsol = make_point(apex1, tol, false, quad1, quad2);
                    a_ptsol.param_on_line = para;
                    a_ptsol.multiple = true;
                    // 1,3
                    let mut glig = gline_line(linsol, false, trans1, trans2);
                    add_vertex_int(&mut glig, a_ptsol.clone());
                    glig.first_point = Some(glig.vertices.len());
                    slin.push(glig);
                    // 2,4
                    let mut linsol_r = linsol;
                    linsol_r.direction = -linsol_r.direction;
                    let para_r = ecl::line_parameter(&linsol_r, apex1);
                    a_ptsol.param_on_line = para_r;
                    let mut glig2 = gline_line(linsol_r, false, trans1, trans2);
                    add_vertex_int(&mut glig2, a_ptsol.clone());
                    glig2.first_point = Some(glig2.vertices.len());
                    slin.push(glig2);
                }
            }
        }
        AnaResultType::Point => {
            let paramapex1 = ecl::line_parameter_of_axis(co1.apex, co1.axis_dir(), apex1);
            let paramapex2 = ecl::line_parameter_of_axis(co2.apex, co2.axis_dir(), apex2);
            for i in 1..=nb_sol {
                let ptcontact = inter.point(i as i32);
                let param1 = ecl::line_parameter_of_axis(co1.apex, co1.axis_dir(), ptcontact);
                let param2 = ecl::line_parameter_of_axis(co2.apex, co2.axis_dir(), ptcontact);
                let (u1, v1) = quad1.parameters(ptcontact);
                let (u2, v2) = quad2.parameters(ptcontact);
                if apex1.distance(ptcontact) <= tol && apex2.distance(ptcontact) <= tol {
                    spnt.push(IntPatchPoint {
                        p1: ptcontact,
                        p2: ptcontact,
                        u1,
                        v1,
                        u2,
                        v2,
                        tolerance: tol,
                    });
                } else if param1 >= paramapex1 && param2 >= paramapex2 {
                    spnt.push(IntPatchPoint {
                        p1: ptcontact,
                        p2: ptcontact,
                        u1,
                        v1,
                        u2,
                        v2,
                        tolerance: tol,
                    });
                }
            }
        }
        AnaResultType::Circle => {
            for i in 1..=nb_sol {
                let cirsol = inter.circle_n(i as i32);
                let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
                let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
                let (t1, t2) = transition_from_scalar(qwe);
                let mut glig = gline_circle(cirsol, false, t1, t2);
                if inter.has_common_gen() {
                    let a_p_char = inter.p_char();
                    let (u1, v1) = quad1.parameters(a_p_char);
                    let (u2, v2) = quad2.parameters(a_p_char);
                    let mut a_ptsol = make_point(a_p_char, tol, false, quad1, quad2);
                    let _ = (u1, v1, u2, v2);
                    a_ptsol.param_on_line = 0.0;
                    add_vertex_int(&mut glig, a_ptsol);
                }
                slin.push(glig);
            }
        }
        AnaResultType::Ellipse => {
            let elipsol = inter.ellipse();
            let (ptref, tgt) = ecl::ellipse_d1(&elipsol, 0.0);
            let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
            let (t1, t2) = transition_from_scalar(qwe);
            let mut glig = gline_ellipse(elipsol, false, t1, t2);
            if inter.has_common_gen() {
                let a_p_char = inter.p_char();
                let mut a_ptsol = make_point(a_p_char, tol, false, quad1, quad2);
                a_ptsol.param_on_line = 0.0;
                add_vertex_int(&mut glig, a_ptsol);
            }
            slin.push(glig);
        }
        AnaResultType::Hyperbola => {
            for i in 1..=2 {
                let hyprsol = inter.hyperbola_n(i as i32);
                let tophypr = hyprsol.center + hyprsol.semi_major * hyprsol.major_dir.normalize_or_zero();
                let major = hyprsol.major_dir.normalize_or_zero();
                let normal = hyprsol.normal.normalize_or_zero();
                let tgttop = normal.cross(major).normalize_or_zero();
                let qwe = tgttop.dot(quad2.normale(tophypr).cross(quad1.normale(tophypr)));
                let (t1, t2) = transition_from_scalar(qwe);
                let mut glig = gline_hyperbola(hyprsol, false, t1, t2);
                if inter.has_common_gen() {
                    let a_p_char = inter.p_char();
                    let mut a_ptsol = make_point(a_p_char, tol, false, quad1, quad2);
                    a_ptsol.param_on_line = 0.0;
                    add_vertex_int(&mut glig, a_ptsol);
                }
                slin.push(glig);
            }
        }
        AnaResultType::Parabola => {
            let parabsol = inter.parabola();
            let tgt_orig = parabsol.normal.normalize_or_zero()
                .cross(parabsol.axis_dir.normalize_or_zero())
                .normalize_or_zero();
            let ptran = tgt_orig
                .dot(quad2.normale(parabsol.vertex).cross(quad1.normale(parabsol.vertex)));
            let (t1, t2) = transition_from_scalar(ptran);
            let mut glig = gline_parabola(parabsol, false, t1, t2);
            if inter.has_common_gen() {
                let a_p_char = inter.p_char();
                let mut a_ptsol = make_point(a_p_char, tol, false, quad1, quad2);
                a_ptsol.param_on_line = 0.0;
                add_vertex_int(&mut glig, a_ptsol);
            }
            slin.push(glig);
        }
        AnaResultType::NoGeometricSolution => {
            let ok = int_quad_quad_fallback(
                quad1,
                quad2,
                &rcad_kernel::geom::Surface3::Cone(co1),
                quad2,
                tol,
                empty,
                multpoint,
                slin,
                spnt,
                false,
                1.0e-9,
                true,
            );
            if !ok {
                return false;
            }
        }
        _ => {
            return false;
        }
    }

    // OCCT L9147-9172: common generatrix through the apexes.
    if inter.has_common_gen() {
        let a_p_char = inter.p_char();
        let linsol = Line3 {
            origin: apex1,
            direction: (apex2 - apex1).normalize_or_zero(),
        };
        let mut glig = gline_line(linsol, true, TypeTrans::Undecided, TypeTrans::Undecided);
        let (u1, v1) = quad1.parameters(a_p_char);
        let (u2, v2) = quad2.parameters(a_p_char);
        let mut a_ptsol = make_point(a_p_char, tol, false, quad1, quad2);
        let _ = (u1, v1, u2, v2);
        let para = ecl::line_parameter(&linsol, a_p_char);
        a_ptsol.param_on_line = para;
        add_vertex_int(&mut glig, a_ptsol);
        slin.push(glig);
    }
    true
}

/// OCCT: distance from a point to an infinite line (gce_MakeLin proximity).
fn point_line_distance(p: DVec3, l: &Line3) -> f64 {
    let d = p - l.origin;
    let proj = d.dot(l.direction.normalize_or_zero());
    (d - l.direction.normalize_or_zero() * proj).length()
}

/// OCCT IntCyCo (L8374-8602) — Cylinder/Cone (and reverse).
#[allow(clippy::too_many_arguments)]
pub fn int_cyco(
    quad1: &Quadric,
    quad2: &Quadric,
    tol: f64,
    reversed: bool,
    empty: &mut bool,
    multpoint: &mut bool,
    slin: &mut Vec<IntPatchLine>,
    spnt: &mut Vec<IntPatchPoint>,
) -> bool {
    let (cy, co) = if !reversed { (quad1, quad2) } else { (quad2, quad1) };

    let mut inter = QuadQuadGeo::new();
    // OCCT L8395-8405: IntAna_QuadQuadGeo inter(Cy, Co, Tol) with Cy the
    // cylinder and Co the cone, in the reordered order.
    inter.perform_cylinder_cone(cy, co, tol);
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
        AnaResultType::Circle => {
            for j in 1..=2 {
                let cirsol = inter.circle_n(j as i32);
                let (ptref, tgt) = ecl::circle_d1(&cirsol, 0.0);
                let qwe = tgt.dot(quad2.normale(ptref).cross(quad1.normale(ptref)));
                let (t1, t2) = transition_from_scalar(qwe);
                slin.push(gline_circle(cirsol, false, t1, t2));
            }
        }
        AnaResultType::NoGeometricSolution => {
            // OCCT L8465-8593: IntAna_IntQuadQuad(Cy, Co) with ExploreCurve
            // splitting each curve at the cone apex.
            return int_quad_quad_fallback(
                quad1,
                quad2,
                &rcad_kernel::geom::Surface3::Cylinder(cy.cylinder()),
                co,
                tol,
                empty,
                multpoint,
                slin,
                spnt,
                true,
                1.0e-8,
                false,
            );
        }
        _ => {
            return false;
        }
    }
    true
}

/// OCCT ExploreCurve (L8608-8675) — split theCrv at the cone apex points.
fn explore_curve(
    the_co: &rcad_kernel::geom::ConicalSurface,
    the_crv: &super::int_quad_quad::IntAnaCurve,
    the_tol: f64,
    the_lc: &mut Vec<super::int_quad_quad::IntAnaCurve>,
) -> bool {
    let a_sq_tol = the_tol * the_tol;
    let a_papx = the_co.apex_point();
    let d0 = the_crv.domain();
    let (mut a_t1, a_t2) = (d0[0], d0[1]);

    the_lc.clear();
    let a_l_params = the_crv.find_parameter(a_papx);
    if a_l_params.is_empty() {
        the_lc.push(the_crv.clone());
        return false;
    }

    for &a_prm in &a_l_params {
        let mut a_prm = a_prm;
        if a_prm - a_t1 < rcad_kernel::precision::PCONFUSION {
            continue;
        }
        let mut is_last = false;
        if a_t2 - a_prm < rcad_kernel::precision::PCONFUSION {
            a_prm = a_t2;
            is_last = true;
        }
        let a_p = the_crv.value(a_prm).unwrap_or(DVec3::ZERO);
        let a_sq_d = a_p.distance_squared(a_papx);
        if a_sq_d < a_sq_tol {
            let mut a_c1 = the_crv.clone();
            a_c1.set_domain(a_t1, a_prm);
            a_t1 = a_prm;
            the_lc.push(a_c1);
        }
        if is_last {
            break;
        }
    }

    if the_lc.is_empty() {
        the_lc.push(the_crv.clone());
        return false;
    }
    if a_t2 - a_t1 > rcad_kernel::precision::PCONFUSION {
        let mut a_c1 = the_crv.clone();
        a_c1.set_domain(a_t1, a_t2);
        the_lc.push(a_c1);
    }
    true
}

/// OCCT: qwe > 1e-7 -> Out/In, < -1e-7 -> In/Out, else Undecided (IntCySp L8243).
fn transition_from_scalar_1e7(qwe: f64) -> (TypeTrans, TypeTrans) {
    if qwe > 1.0e-7 {
        (TypeTrans::Out, TypeTrans::In)
    } else if qwe < -1.0e-7 {
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
            tangent: false,
            multiple: v.multiple,
            on_dom_s1: v.on_dom_s1,
            on_dom_s2: v.on_dom_s2,
            arc_on_s1: v.arc_on_s1,
            arc_on_s2: v.arc_on_s2,
            param_on_arc1: v.param_on_arc1,
            param_on_arc2: v.param_on_arc2,
            is_vertex_on_s1: v.is_vertex_on_s1,
            is_vertex_on_s2: v.is_vertex_on_s2,
            transition_line_arc1: v.transition_line_arc1,
            transition_line_arc2: v.transition_line_arc2,
            transition_on_s1: v.transition_on_s1,
            transition_on_s2: v.transition_on_s2,
        });
    }
}
