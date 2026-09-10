//! OCCT BRepFill_Sweep.cxx, part B2 (the Filling static, cxx L867-1337) —
//! the second half of the part-B statics, split from [`super::
//! brep_fill_sweep_b`] to keep every file under 2,000 lines (AGENTS.md
//! File Editing Rule 5).  The GeomPlate deformation block L1003-1179 is
//! commented out in the OCCT source itself and carries no code.

use std::collections::HashMap;

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::{p_confusion, CONFUSION, SQUARE_CONFUSION};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d, Surface3, SurfaceEval, TrimmedCurve2};
use rcad_kernel::topo::topods::{BRep, BRepBuilder, GeomAbsShape, Orientation, Shape};

use crate::brep_fill::generator::ShapeKey;

use super::brep_fill_sweep_b::{
    build_edge_c3d, build_face, curve2d_param_bounds, gp_vec2d_angle, gp_vec_angle,
    point_in_ax2_frame, set_vertex_tolerance, update_edge_pcurve2_on_surf,
    update_edge_pcurve_on_surf, GpAx2,
};

// ---------------------------------------------------------------------------
// Filling (cxx L867-1337)
// ---------------------------------------------------------------------------

/// OCCT static Filling (cxx L867-1337) — construct the faces of filling.
#[allow(clippy::too_many_arguments)]
pub(super) fn filling(
    brep: &mut BRep,
    ef: &Shape,
    f1: &Shape,
    el: &Shape,
    f2: &Shape,
    eemap: &mut HashMap<ShapeKey, Shape>,
    tol: f64,
    axe: &GpAx2,
    tangent_on_part1: DVec3,
    aux1: &mut Shape,
    aux2: &mut Shape,
    result: &mut Shape,
) -> bool {
    let mut b = BRepBuilder::new();
    let _ = tol;
    let mut with_e3;
    let mut with_e4;

    // Return constraints
    let mut e1: Shape;
    let mut e2: Shape;
    let mut e3: Shape;
    let mut e4: Shape;
    e1 = ef.clone();
    e2 = el.clone();

    // TopExp::Vertices(E1, Vf, Vl);
    let (mut vf, mut vl) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(&e1);
    vf.orientation = Orientation::Forward;
    vl.orientation = Orientation::Forward;

    // TopExp::Vertices(E2, V1, V2);
    let (mut v1, mut v2) = crate::brep_fill::brep_fill_pipe::top_exp_vertices(&e2);
    v1.orientation = Orientation::Reversed;
    v2.orientation = Orientation::Reversed;

    // B.MakeEdge(E3); B.MakeEdge(E4);
    e3 = b.add_edge(brep, None, Shape::null(), Shape::null(), [0.0, 0.0]);
    e4 = b.add_edge(brep, None, Shape::null(), Shape::null(), [0.0, 0.0]);

    with_e3 = false;
    with_e4 = false;

    if !aux1.is_null() && !vf.is_same(&v1) {
        e3 = aux1.clone();
        with_e3 = true;
    }

    if vf.is_same(&vl) {
        e4 = e3.clone();
        // E4.Reverse();
        e4.orientation = match e4.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
        with_e4 = with_e3;
    } else if !aux2.is_null() && !vl.is_same(&v2) {
        e4 = aux2.clone();
        with_e4 = true;
    }

    // Construction of a surface of revolution
    // occ::handle<Geom_Curve> Prof1, Prof2;
    let ed1 = brep.edge(e1.clone());
    let prof1 = ed1
        .curve
        .clone()
        .expect("Filling: E1 without 3d curve");
    let (f1p, l1p) = (ed1.range[0], ed1.range[1]);
    let ed2 = brep.edge(e2.clone());
    let prof2 = ed2
        .curve
        .clone()
        .expect("Filling: E2 without 3d curve");
    let (f2p, l2p) = (ed2.range[0], ed2.range[1]);

    // Indeed, both Prof1 and Prof2 are the same curves but in different
    // positions ... bool bSameCurveDomain =
    // EF.Orientation() != EL.Orientation();
    let b_same_curve_domain = ef.orientation != el.orientation;

    // Choose the angle of the opening — gp_Trsf aTf; aTf.SetTransformation(Axe);
    // Choose the furthest point from the "center of revolution" to provide
    // correct angle measurement.
    let a_prm = [f1p, 0.5 * (f1p + l1p), l1p];
    let a_p1 = [
        point_in_ax2_frame(axe, prof1.point_at(a_prm[0])),
        point_in_ax2_frame(axe, prof1.point_at(a_prm[1])),
        point_in_ax2_frame(axe, prof1.point_at(a_prm[2])),
    ];

    let mut a_max_idx = -1i32;
    let mut a_max_dist = f64::MIN;
    for i in 0..3 {
        let a_dist = a_p1[i].x * a_p1[i].x + a_p1[i].z * a_p1[i].z;
        if a_dist > a_max_dist {
            a_max_dist = a_dist;
            a_max_idx = i as i32;
        }
    }

    let a_prm2 = [f2p, 0.5 * (f2p + l2p), l2p];
    let a_p2 = point_in_ax2_frame(
        axe,
        prof2.point_at(a_prm2[if b_same_curve_domain {
            a_max_idx as usize
        } else {
            2 - a_max_idx as usize
        }]),
    );
    let a_v1 = DVec2::new(a_p1[a_max_idx as usize].z, a_p1[a_max_idx as usize].x);
    let a_v2 = DVec2::new(a_p2.z, a_p2.x);
    let gp_resolution = rcad_kernel::math::gp::GP_RESOLUTION;
    if a_v1.length_squared() <= gp_resolution || a_v2.length_squared() <= gp_resolution {
        return false;
    }

    // Angle = aV1.Angle(aV2);
    let mut angle = gp_vec2d_angle(a_v1, a_v2);

    // gp_Ax1 axe(Axe.Location(), Axe.YDirection());
    let axe1 = (axe.location, axe.y_direction);

    if angle < 0.0 {
        angle = -angle;
        // axe.Reverse();
        let axe1_rev = (axe1.0, -axe1.1);
        return filling_continue(
            brep, ef, f1, el, f2, eemap, tol, axe, axe1_rev, tangent_on_part1,
            aux1, aux2, result,
            e1, e2, e3, e4, vf, vl, v1, v2, with_e3, with_e4,
            prof1, prof2, f1p, l1p, f2p, l2p, a_prm, a_max_idx, angle,
        );
    }

    filling_continue(
        brep, ef, f1, el, f2, eemap, tol, axe, axe1, tangent_on_part1,
        aux1, aux2, result,
        e1, e2, e3, e4, vf, vl, v1, v2, with_e3, with_e4,
        prof1, prof2, f1p, l1p, f2p, l2p, a_prm, a_max_idx, angle,
    )
}

/// The continuation of OCCT Filling after the angle opening choice (the
/// OCCT body continues in the same frame; the rcad split carries the
/// locals — no control-flow change, the split is at the `axe.Reverse()`
/// conditional tail only).
#[allow(clippy::too_many_arguments)]
fn filling_continue(
    brep: &mut BRep,
    _ef: &Shape,
    f1: &Shape,
    _el: &Shape,
    f2: &Shape,
    eemap: &mut HashMap<ShapeKey, Shape>,
    _tol: f64,
    _axe: &GpAx2,
    axe1: (DVec3, DVec3),
    tangent_on_part1: DVec3,
    aux1: &mut Shape,
    aux2: &mut Shape,
    result: &mut Shape,
    mut e1: Shape,
    mut e2: Shape,
    mut e3: Shape,
    mut e4: Shape,
    mut vf: Shape,
    mut vl: Shape,
    mut v1: Shape,
    mut v2: Shape,
    mut with_e3: bool,
    mut with_e4: bool,
    prof1: Curve3,
    _prof2: Curve3,
    f1p: f64,
    l1p: f64,
    f2p: f64,
    l2p: f64,
    a_prm: [f64; 3],
    a_max_idx: i32,
    mut angle: f64,
) -> bool {
    // occ::handle<Geom_SurfaceOfRevolution> Rev = new
    //   (Geom_SurfaceOfRevolution)(Prof1, axe);
    let rev = Surface3::Revolution(rcad_kernel::geom::RevolutionSurface {
        profile: Box::new(prof1.clone()),
        axis_origin: axe1.0,
        axis_dir: axe1.1,
    });

    // occ::handle<Geom_Surface> Surf = new (Geom_RectangularTrimmedSurface)
    //   (Rev, 0, Angle, f1, l1);
    let mut surf = Surface3::Trimmed(rcad_kernel::geom::TrimmedSurface {
        basis: Box::new(rev),
        trim: [0.0, angle, f1p, l1p],
    });

    // Control the direction of the rotation
    let mut to_reverse_result = false;
    // gp_Vec d1u; d1u = Surf->DN(0, aPrm[aMaxIdx], 1, 0);
    let (_p, d1u, _dv) = surf.derivatives(0.0, a_prm[a_max_idx as usize]);
    if gp_vec_angle(d1u, tangent_on_part1) > std::f64::consts::PI / 2.0 {
        // Invert everything
        to_reverse_result = true;
    }

    // occ::handle<Geom2d_Line> L; gp_Pnt2d P2d(0., 0.);
    // L = new (Geom2d_Line)(P2d, gp::DY2d()); C1 = new (Geom2d_TrimmedCurve)(L, f1, l1);
    let mut c1 = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(Line2d::new(
            DVec2::new(0.0, 0.0),
            DVec2::new(0.0, 1.0),
        ))),
        t_min: f1p,
        t_max: l1p,
    });

    // P2d.SetCoord(Angle, 0.); L = new (Geom2d_Line)(P2d, gp::DY2d());
    // C2 = new (Geom2d_TrimmedCurve)(L, f1, l1);
    let mut c2 = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(Line2d::new(
            DVec2::new(angle, 0.0),
            DVec2::new(0.0, 1.0),
        ))),
        t_min: f1p,
        t_max: l1p,
    });

    // It is required to control the direction and the range.
    // C2->D0(f1, P2d); Surf->D0(P2d.X(), P2d.Y(), P1);
    let mut p2d = c2.point_at(f1p);
    let mut p1 = surf.point_at(p2d.x, p2d.y);
    // C2->D0(l1, P2d); Surf->D0(P2d.X(), P2d.Y(), P2);
    p2d = c2.point_at(l1p);
    let p2 = surf.point_at(p2d.x, p2d.y);
    // P = BRep_Tool::Pnt(V1);
    let p = brep.vertex(v1.clone()).point;
    if p.distance(p2) + CONFUSION < p.distance(p1) {
        // E2 is parsed in the direction opposite to E1
        // C2->Reverse();
        c2 = reversed_trimmed_line_basis(&c2);
        let aux = v2.clone();
        v2 = v1.clone();
        v1 = aux;
    }
    // GeomLib::SameRange(PConfusion(), C2, C2->FirstParameter(),
    //                    C2->LastParameter(), f2, l2, C3); C2 = C3;
    let (c2_first, c2_last) = curve2d_param_bounds(&c2);
    let c3 = crate::geomalgo::geom_lib_same_range::same_range(
        p_confusion(),
        &c2,
        c2_first,
        c2_last,
        f2p,
        l2p,
    );
    c2 = c3;

    // P1 = BRep_Tool::Pnt(Vf); P2 = BRep_Tool::Pnt(V1); (pointu_f comment)
    let _p1 = brep.vertex(vf.clone()).point;
    let _p2 = brep.vertex(v1.clone()).point;
    // P1 = BRep_Tool::Pnt(Vl); P2 = BRep_Tool::Pnt(V2); (pointu_l comment)
    let _p1 = brep.vertex(vl.clone()).point;
    let _p2 = brep.vertex(v2.clone()).point;

    // P2d.SetCoord(0., f1); L = new (Geom2d_Line)(P2d, gp::DX2d());
    // C3 = new (Geom2d_TrimmedCurve)(L, 0, Angle);
    let c3 = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(Line2d::new(
            DVec2::new(0.0, f1p),
            DVec2::new(1.0, 0.0),
        ))),
        t_min: 0.0,
        t_max: angle,
    });
    // (C3/C4 are consumed through the UpdateEdge calls below)
    let c3 = c3;

    // P2d.SetCoord(0., l1); L = new (Geom2d_Line)(P2d, gp::DX2d());
    // C4 = new (Geom2d_TrimmedCurve)(L, 0, Angle);
    let c4 = Curve2d::Trimmed(TrimmedCurve2 {
        curve: Box::new(Curve2d::Line(Line2d::new(
            DVec2::new(0.0, l1p),
            DVec2::new(1.0, 0.0),
        ))),
        t_min: 0.0,
        t_max: angle,
    });

    // (the GeomPlate deformation block L1003-1179 is commented out in OCCT)

    // Update des Edges
    // TopLoc_Location Loc; occ::handle<Geom_Curve> C3d;
    // B.UpdateEdge(E1, C1, Surf, Loc, Precision::Confusion());
    update_edge_pcurve_on_surf(brep, &e1, &c1, &surf, 0, CONFUSION);
    // B.UpdateEdge(E2, C2, Surf, Loc, Precision::Confusion());
    update_edge_pcurve_on_surf(brep, &e2, &c2, &surf, 0, CONFUSION);

    let c3 = c3;
    if e3.is_same(&e4) {
        if !with_e3 {
            // C3d = Surf->VIso(f1);
            let c3d = crate::brep_fill::brep_fill_sweep::surface_viso(&surf, f1p);
            // E3 = BuildEdge(C3d, C3, Surf, Vf, V1, 0, Angle, PConfusion);
            e3 = build_edge_c3d(brep, &Some(c3d), &c3, &surf, &vf, &v1, 0.0, angle, CONFUSION);
        } else {
            // BRepAdaptor_Curve aCurve(E3); AngleOld = aCurve.LastParameter();
            let angle_old = brep.edge(e3.clone()).range[1];
            if angle > angle_old {
                // B.Range(E3, 0, Angle);
                BRepBuilder::new().set_edge_range(brep, e3.clone(), 0.0, angle);
                // TopoDS_Vertex V(TopExp::LastVertex(E3));
                let ed_last = brep.edge(e3.clone()).last.clone();
                // TVlast->Tolerance(Precision::Confusion());
                set_vertex_tolerance(brep, &ed_last, CONFUSION);
            }
        }

        // B.UpdateEdge(E3, C3, C4, Surf, Loc, Precision::Confusion());
        update_edge_pcurve2_on_surf(brep, &e3, &c3, &c4, &surf, 0, CONFUSION);
        // E4 = E3; E4.Reverse();
        e4 = e3.clone();
        e4.orientation = match e4.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
    } else {
        if !with_e3 {
            let c3d = crate::brep_fill::brep_fill_sweep::surface_viso(&surf, f1p);
            e3 = build_edge_c3d(brep, &Some(c3d), &c3, &surf, &vf, &v1, 0.0, angle, CONFUSION);
        } else {
            let angle_old = brep.edge(e3.clone()).range[1];
            if angle > angle_old {
                BRepBuilder::new().set_edge_range(brep, e3.clone(), 0.0, angle);
                let ed_last = brep.edge(e3.clone()).last.clone();
                set_vertex_tolerance(brep, &ed_last, CONFUSION);
            }
            // B.UpdateEdge(E3, C3, Surf, Loc, Precision::Confusion());
            update_edge_pcurve_on_surf(brep, &e3, &c3, &surf, 0, CONFUSION);
        }

        if !with_e4 {
            // C3d = Surf->VIso(l1);
            let c3d = crate::brep_fill::brep_fill_sweep::surface_viso(&surf, l1p);
            // E4 = BuildEdge(C3d, C4, Surf, Vl, V2, 0, Angle, PConfusion);
            e4 = build_edge_c3d(brep, &Some(c3d), &c4, &surf, &vl, &v2, 0.0, angle, CONFUSION);
        } else {
            let angle_old = brep.edge(e4.clone()).range[1];
            if angle > angle_old {
                BRepBuilder::new().set_edge_range(brep, e4.clone(), 0.0, angle);
                let ed_last = brep.edge(e4.clone()).last.clone();
                set_vertex_tolerance(brep, &ed_last, CONFUSION);
            }
            // B.UpdateEdge(E4, C4, Surf, Loc, Precision::Confusion());
            update_edge_pcurve_on_surf(brep, &e4, &c4, &surf, 0, CONFUSION);
        }
    }

    // Construct face
    build_face(brep, &surf, &e1, &e3, &e2, &e4, eemap, false, false, result);

    // Set the continuities.
    // B.Continuity(E1, TopoDS::Face(F1), Result, GeomAbs_G1);
    BRepBuilder::new().continuity(brep, &e1, f1, result, GeomAbsShape::G1);
    // B.Continuity(E2, TopoDS::Face(F2), Result, GeomAbs_G1);
    BRepBuilder::new().continuity(brep, &e2, f2, result, GeomAbsShape::G1);

    // Render the calculated borders.
    if !brep.edge(e3.clone()).degenerated {
        // Aux1 = E3;
        *aux1 = e3.clone();
    } else {
        // B.MakeEdge(Aux1); // Nullify
        *aux1 = BRepBuilder::new().add_edge(brep, None, Shape::null(), Shape::null(), [0.0, 0.0]);
    }

    if !brep.edge(e4.clone()).degenerated {
        // Aux2 = E4;
        *aux2 = e4.clone();
    } else {
        // B.MakeEdge(Aux2);
        *aux2 = BRepBuilder::new().add_edge(brep, None, Shape::null(), Shape::null(), [0.0, 0.0]);
    }

    // Set the orientation
    // Provide correct normals computation (the normal will be computed not
    // in singularity point definitely).
    let mut angle = f64::MIN;
    for i in 0..3 {
        // gp_Vec D1U, D1V, N1, N2;
        // C1->D0(aPrm[i], P2d); Surf->D1(P2d.X(), P2d.Y(), P, D1U, D1V);
        let p2d = c1.point_at(a_prm[i]);
        let (_p, d1u, d1v) = surf.derivatives(p2d.x, p2d.y);
        // N1 = D1U ^ D1V;
        let n1 = d1u.cross(d1v);

        if n1.length_squared() < SQUARE_CONFUSION {
            continue;
        }

        // C1 = BRep_Tool::CurveOnSurface(E1, TopoDS::Face(F1), f2, l2);
        let (c1_on_f1, _tf, _tl) =
            crate::brep_algo::tool::brep_tool_curve_on_surface(&e1, f1).expect("C1 on F1");
        let _ = &mut c1;
        // C1->D0(aPrm[i], P2d);
        let p2d = c1_on_f1.point_at(a_prm[i]);
        // occ::handle<BRepAdaptor_Surface> AS = new BRepAdaptor_Surface(Face(F1));
        // AS->D1(P2d.X(), P2d.Y(), P, D1U, D1V);
        let as_surf = rcad_kernel::topo::topods::face_surface_value(brep, f1)
            .expect("Filling: F1 without surface");
        let (_p, d1u, d1v) = as_surf.derivatives(p2d.x, p2d.y);
        // N2 = D1U ^ D1V;
        let n2 = d1u.cross(d1v);

        if n2.length_squared() < SQUARE_CONFUSION {
            continue;
        }

        // Angle = N1.Angle(N2);
        angle = gp_vec_angle(n1, n2);

        break;
    }

    if angle == f64::MIN {
        return false;
    }

    // if ((F1.Orientation() == TopAbs_REVERSED) ^ (Angle > M_PI / 2))
    if (f1.orientation == Orientation::Reversed) ^ (angle > std::f64::consts::PI / 2.0) {
        result.orientation = Orientation::Reversed;
    } else {
        result.orientation = Orientation::Forward;
    }

    if to_reverse_result {
        // Result.Reverse();
        result.orientation = match result.orientation {
            Orientation::Forward => Orientation::Reversed,
            Orientation::Reversed => Orientation::Forward,
            o => o,
        };
    }

    true
}

/// OCCT Geom2d_TrimmedCurve::Reverse over the trimmed-line C2 (the basis
/// line direction is reversed, the bounds swap).
fn reversed_trimmed_line_basis(c2: &Curve2d) -> Curve2d {
    match c2 {
        Curve2d::Trimmed(t) => {
            let basis = match t.curve.as_ref() {
                Curve2d::Line(l) => Curve2d::Line(Line2d {
                    origin: l.origin,
                    direction: -l.direction,
                }),
                _ => panic!("GAP: Geom2d_TrimmedCurve::Reverse basis — Filling"),
            };
            Curve2d::Trimmed(TrimmedCurve2 {
                curve: Box::new(basis),
                t_min: t.t_max,
                t_max: t.t_min,
            })
        }
        _ => panic!("GAP: Geom2d_TrimmedCurve::Reverse — Filling"),
    }
}
