// OCCT BRepOffset_Inter2d.cxx L1193-2394 — 1:1 translation (part B: the
// class entry points ExtentEdge / Compute / ConnexIntByInt /
// ConnexIntByIntInVert / FuseVertices and the file statics UpdateVertex /
// MakeChain; the file statics through ExtendPCurve live in
// brep_offset_inter2d.rs).
//
// Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffset/
//         BRepOffset_Inter2d.cxx / .hxx
//
// The architecture-difference list is the header of
// brep_offset_inter2d.rs (items #21-#34).

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{
    Curve2d, Curve2dEval, Curve3, CurveEval, Line3, Surface3, SurfaceEval, TrimmedCurve3,
};
use rcad_kernel::topods::{CurveRepresentation, GeomAbsShape, Orientation, ShapeType, TShape};
use rcad_kernel::topo_shape::Shape;

use crate::brep_algo::as_des::BRepAlgoAsDes;
use crate::brep_algo::image::BRepAlgoImage;
use crate::brep_algo::tool::{
    brep_tool_pnt, brep_tool_range, brep_tool_tolerance, builder_range_edge, empty_copied,
    explorer, oriented, shape_key, top_abs_reverse, ShapeKey,
};
use crate::feat::loc_ope_wires_on_shape_b::brep_tool_degenerated;

use super::brep_offset_inter2d::*;
use super::brep_offset_offset_b::BRepOffsetOffset;

// ---------------------------------------------------------------------------
// GAP leaves local to this part.
// ---------------------------------------------------------------------------

/// OCCT theSurf->Transformed(theLoc.Transformation()) — GAP (architecture
/// difference #32: the identity-location convention carries the surface
/// through; the location bake is not translated).
fn geom_surface_transformed_loc(_the_s: &Surface3) -> Surface3 {
    panic!("GAP: Geom_Surface::Transformed (location bake not translated)");
}

/// OCCT BRep_CurveRepresentation::Location() — the identity-location
/// convention (architecture difference #32).
fn curve_rep_location(_the_rep: &CurveRepresentation) -> u32 {
    0
}

// =========================================================================
// OCCT BRepOffset_Inter2d — the class (BRepOffset_Inter2d.hxx L37-119).
// =========================================================================

/// OCCT BRepOffset_Inter2d (BRepOffset_Inter2d.hxx L37-119) — the static-only
/// class computing the intersections between the edges stored in AsDes as
/// descendants of a face; the results are stored as SD in AsDes.
pub struct BRepOffsetInter2d;

impl BRepOffsetInter2d {
    /// OCCT BRepOffset_Inter2d::ExtentEdge(E, NE, theOffset)
    /// (BRepOffset_Inter2d.cxx L1193-1674) — extents the edge; NE receives
    /// the extended copy.  Returns false on the projection failure paths.
    pub fn extent_edge(e: &Shape, ne: &mut Shape, the_offset: f64) -> bool {
        // OCCT L1195: BRepLib::BuildCurve3d(E); (commented out in OCCT).

        // OCCT L1197-1204.
        let a2_offset = 2. * the_offset.abs();
        let (an_ef, an_el) = brep_tool_range(e);
        // OCCT L1197/1205: NE = TopoDS::Edge(E.EmptyCopied()).
        *ne = empty_copied(e);

        // OCCT L1212-1216.
        let mut nb_pcurves = 0usize;
        let mut first_par_on_pc = -f64::MAX; // OCCT RealFirst().
        let mut last_par_on_pc = f64::MAX; // OCCT RealLast().
        let mut min_pc: Option<Curve2d> = None;
        let mut min_surf: Option<Surface3> = None;
        let mut min_loc: u32 = 0;
        // Architecture difference #32: the OCCT `theCurve == MinPC` handle
        // identity maps to the representation index.
        let mut min_rep_index: Option<usize> = None;

        // OCCT L1218-1402: the representation walk.
        if let TShape::Edge(ed) = Arc::make_mut(&mut ne.data) {
            for i_rep in 0..ed.representations.len() {
                // OCCT L1224: CurveRep->IsCurveOnSurface().
                let is_curve_on_surface = matches!(
                    ed.representations[i_rep],
                    CurveRepresentation::CurveOnSurface { .. }
                        | CurveRepresentation::CurveOnClosedSurface { .. }
                );
                if !is_curve_on_surface {
                    continue;
                }
                // OCCT L1226: NbPCurves++.
                nb_pcurves += 1;
                // OCCT L1227: theCurve = CurveRep->PCurve().
                let mut the_curve: Curve2d = match &ed.representations[i_rep] {
                    CurveRepresentation::CurveOnSurface { pcurve, .. } => pcurve.clone(),
                    CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => pcurve1.clone(),
                    _ => unreachable!(),
                };
                // OCCT L1228-1229.
                let mut first_par = the_curve.default_domain()[0];
                let mut last_par = the_curve.default_domain()[1];

                // OCCT L1231-1249: the bounded-pcurve extension.
                let is_bounded = matches!(
                    the_curve,
                    Curve2d::BSpline(_) | Curve2d::Bezier(_) | Curve2d::Trimmed(_)
                );
                if is_bounded && (first_par > an_ef - a2_offset || last_par < an_el + a2_offset) {
                    let mut new_pcurve: Curve2d = the_curve.clone();
                    if extend_pcurve(&the_curve, an_ef, an_el, a2_offset, &mut new_pcurve) {
                        // OCCT L1237: CurveRep->PCurve(NewPCurve).
                        match &mut ed.representations[i_rep] {
                            CurveRepresentation::CurveOnSurface { pcurve, .. } => {
                                *pcurve = new_pcurve.clone();
                            }
                            CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => {
                                *pcurve1 = new_pcurve.clone();
                            }
                            _ => {}
                        }
                        // OCCT L1238-1239.
                        first_par = new_pcurve.default_domain()[0];
                        last_par = new_pcurve.default_domain()[1];
                        // OCCT L1240-1247: the closed-surface second pcurve.
                        if let CurveRepresentation::CurveOnClosedSurface { pcurve2, .. } =
                            &mut ed.representations[i_rep]
                        {
                            let mut new2: Curve2d = pcurve2.clone();
                            if extend_pcurve(pcurve2, an_ef, an_el, a2_offset, &mut new2) {
                                // OCCT L1245: CurveRep->PCurve2(NewPCurve).
                                *pcurve2 = new2;
                            }
                        }
                    }
                } else if the_curve.is_periodic() {
                    // OCCT L1250-1256.
                    let mut delta = (curve2d_period(&the_curve) - (an_el - an_ef)) * 0.5;
                    delta *= 0.95;
                    first_par = an_ef - delta;
                    last_par = an_el + delta;
                } else if the_curve.is_closed() {
                    // OCCT L1257-1260.
                    last_par -= 0.05 * (last_par - first_par);
                }

                // OCCT L1263: theCurve = CurveRep->PCurve() — re-read (the
                // possibly extended one).
                the_curve = match &ed.representations[i_rep] {
                    CurveRepresentation::CurveOnSurface { pcurve, .. } => pcurve.clone(),
                    CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => pcurve1.clone(),
                    _ => unreachable!(),
                };
                // OCCT L1264: theSurf = CurveRep->Surface() — GAP
                // (architecture difference #32).
                let the_surf = curve_rep_surface(&ed.representations[i_rep]);
                // OCCT L1265-1266: the surface bounds.
                let bounds = the_surf.default_domain();
                let (umin, umax, vmin, vmax) = (bounds[0], bounds[1], bounds[2], bounds[3]);

                // OCCT L1267-1291: the boundary lines on the finite bounds.
                let mut bound_lines: Vec<Curve2d> = Vec::new();
                if !vmin.is_infinite() {
                    bound_lines.push(Curve2d::Line(rcad_kernel::geom::Line2d::new(
                        DVec2::new(0., vmin),
                        DVec2::X,
                    )));
                }
                if !umin.is_infinite() {
                    bound_lines.push(Curve2d::Line(rcad_kernel::geom::Line2d::new(
                        DVec2::new(umin, 0.),
                        DVec2::Y,
                    )));
                }
                if !vmax.is_infinite() {
                    bound_lines.push(Curve2d::Line(rcad_kernel::geom::Line2d::new(
                        DVec2::new(0., vmax),
                        DVec2::X,
                    )));
                }
                if !umax.is_infinite() {
                    bound_lines.push(Curve2d::Line(rcad_kernel::geom::Line2d::new(
                        DVec2::new(umax, 0.),
                        DVec2::Y,
                    )));
                }

                // OCCT L1293-1337: the intersection parameters with the
                // boundaries.
                let mut params: Vec<f64> = Vec::new();
                let p_confusion = rcad_kernel::precision::PCONFUSION;
                for i in 1..=bound_lines.len() {
                    let mut int_cc = Geom2dIntGInter {
                        base: crate::geomalgo::int_res2d::IntersectionBase::new(),
                    };
                    int_cc.perform_natural(&the_curve, &bound_lines[i - 1], p_confusion, p_confusion);
                    if int_cc.is_done() {
                        for j in 1..=int_cc.nb_points() {
                            let ip = int_cc.point(j);
                            let a_point = ip.value();
                            if a_point.x >= umin && a_point.x <= umax && a_point.y >= vmin && a_point.y <= vmax {
                                params.push(ip.param_on_first());
                            }
                        }
                        for j in 1..=int_cc.nb_segments() {
                            let is = int_cc.segment(j);
                            if is.has_first_point() {
                                let ip = is.first_point();
                                let a_point = ip.value();
                                if a_point.x >= umin && a_point.x <= umax && a_point.y >= vmin && a_point.y <= vmax {
                                    params.push(ip.param_on_first());
                                }
                            }
                            if is.has_last_point() {
                                let ip = is.last_point();
                                let a_point = ip.value();
                                if a_point.x >= umin && a_point.x <= umax && a_point.y >= vmin && a_point.y <= vmax {
                                    params.push(ip.param_on_first());
                                }
                            }
                        }
                    }
                }

                // OCCT L1338-1379: clamp the range to the boundary hits.
                if !params.is_empty() {
                    if params.len() == 1 {
                        let pnt_first = the_curve.point_at(first_par);
                        if pnt_first.x >= umin && pnt_first.x <= umax && pnt_first.y >= vmin && pnt_first.y <= vmax {
                            if last_par > params[0] {
                                last_par = params[0];
                            }
                        } else if first_par < params[0] {
                            first_par = params[0];
                        }
                    } else {
                        let mut fpar = f64::MAX;
                        let mut lpar = -f64::MAX;
                        for i in 1..=params.len() {
                            if params[i - 1] < fpar {
                                fpar = params[i - 1];
                            }
                            if params[i - 1] > lpar {
                                lpar = params[i - 1];
                            }
                        }
                        if first_par < fpar {
                            first_par = fpar;
                        }
                        if last_par > lpar {
                            last_par = lpar;
                        }
                    }
                }
                // OCCT L1381: end of check — BRep_GCurve::SetRange.
                match &mut ed.representations[i_rep] {
                    CurveRepresentation::CurveOnSurface { range, .. } => {
                        *range = [first_par, last_par];
                    }
                    CurveRepresentation::CurveOnClosedSurface { range, .. } => {
                        *range = [first_par, last_par];
                    }
                    _ => {}
                }

                // OCCT L1386-1400: update FirstParOnPC and LastParOnPC.
                if first_par > first_par_on_pc {
                    first_par_on_pc = first_par;
                    min_pc = Some(the_curve.clone());
                    min_surf = Some(the_surf.clone());
                    min_loc = curve_rep_location(&ed.representations[i_rep]);
                    min_rep_index = Some(i_rep);
                }
                if last_par < last_par_on_pc {
                    last_par_on_pc = last_par;
                    min_pc = Some(the_curve);
                    min_surf = Some(the_surf);
                    min_loc = curve_rep_location(&ed.representations[i_rep]);
                    min_rep_index = Some(i_rep);
                }
            }
        }

        // OCCT L1404-1405: the 3D curve of NE.
        let c3d_read = match ne.as_edge() {
            Some(ed) => ed.curve.clone().map(|c| (c, ed.range[0], ed.range[1])),
            None => None,
        };
        let mut f = c3d_read.as_ref().map(|(_, fr, _)| *fr).unwrap_or(0.0);
        let mut l = c3d_read.as_ref().map(|(_, _, lr)| *lr).unwrap_or(0.0);
        let mut c3d: Option<Curve3> = c3d_read.map(|(c, _, _)| c);

        if nb_pcurves != 0 {
            // OCCT L1408: MinLoc = E.Location() * MinLoc — the
            // identity-location convention (architecture difference #32).
            if let Some(c3d_ref) = &c3d {
                let min_pc = min_pc.as_ref().expect("MinPC");
                let min_surf = min_surf.as_ref().expect("MinSurf");
                // OCCT L1411-1415.
                if min_pc.is_closed() {
                    f = first_par_on_pc;
                    l = last_par_on_pc;
                } else if c3d_ref.is_periodic() {
                    // OCCT L1416-1422.
                    let mut delta = (curve3_period(c3d_ref) - (l - f)) * 0.5;
                    delta *= 0.95;
                    f -= delta;
                    l += delta;
                } else if c3d_ref.is_closed() {
                    // OCCT L1423-1426.
                    l -= 0.05 * (l - f);
                } else {
                    // OCCT L1427-1462: the projector block.
                    f = first_par_on_pc;
                    l = last_par_on_pc;
                    // OCCT L1431: the MinLoc.Transformation() — the
                    // identity-location convention.
                    if !first_par_on_pc.is_infinite() {
                        let p2d1 = min_pc.point_at(first_par_on_pc);
                        let p1 = min_surf.point_at(p2d1.x, p2d1.y);
                        let projector = GeomAPIProjectPointOnCurve::init_point_curve(p1, c3d_ref);
                        if projector.nb_points() > 0 {
                            f = projector.lower_distance_parameter();
                        }
                    }
                    if !last_par_on_pc.is_infinite() {
                        let p2d2 = min_pc.point_at(last_par_on_pc);
                        let p2 = min_surf.point_at(p2d2.x, p2d2.y);
                        let projector = GeomAPIProjectPointOnCurve::init_point_curve(p2, c3d_ref);
                        if projector.nb_points() > 0 {
                            l = projector.lower_distance_parameter();
                        }
                    }
                }
                // OCCT L1464-1468.
                builder_range_edge(ne, f, l);
                if !f.is_infinite() && !l.is_infinite() {
                    brep_lib_same_parameter(ne, rcad_kernel::precision::CONFUSION);
                }
            } else if !brep_tool_degenerated(e) {
                // OCCT L1470-1603: no 3d curve.
                // OCCT L1472: MinSurf = MinSurf->Transformed(
                //             MinLoc.Transformation()) — the
                // identity-location convention carries the surface through.
                let min_surf = match min_surf {
                    Some(s) => s,
                    None => panic!("ExtentEdge: MinSurf is null"),
                };
                let min_pc = min_pc.as_ref().expect("MinPC");
                let mut max_deviation = 0.;
                if first_par_on_pc.is_infinite() || last_par_on_pc.is_infinite() {
                    // OCCT L1476-1501: the line construction.
                    if let Curve2d::Line(the_line) = min_pc {
                        let mut is_line = false;
                        if matches!(min_surf, Surface3::Plane(_)) {
                            is_line = true;
                        } else if matches!(min_surf, Surface3::Cylinder(_) | Surface3::Cone(_)) {
                            // OCCT L1486-1491: the line direction parallel
                            // to gp::DY2d().
                            if gp_dir2d_is_parallel(
                                the_line.direction,
                                DVec2::Y,
                                rcad_kernel::precision::ANGULAR,
                            ) {
                                is_line = true;
                            }
                        }
                        if is_line {
                            // OCCT L1493-1499: the two line points at 0/1.
                            let p2d1 = the_line.point_at(0.);
                            let p2d2 = the_line.point_at(1.);
                            let p1 = min_surf.point_at(p2d1.x, p2d1.y);
                            let p2 = min_surf.point_at(p2d2.x, p2d2.y);
                            let a_vec = p2 - p1;
                            c3d = Some(Curve3::Line(Line3::new(p1, a_vec)));
                        }
                    }
                } else {
                    // OCCT L1503-1524: the Adaptor3d_CurveOnSurface +
                    // GeomLib::BuildCurve3d — GAP (architecture difference
                    // #29).
                    let con_s = Adaptor3dCurveOnSurface::new(min_pc, &min_surf);
                    let _continuity = GeomAbsShape::C1; // OCCT GeomAbs_C1.
                    let _max_degree = 14i32;
                    let _max_segment = evaluate_max_segment(&con_s);
                    // OCCT L1510: double /*max_deviation,*/ average_deviation;
                    let mut average_deviation = 0.;
                    geom_lib_build_curve3d(&mut c3d, &mut max_deviation, &mut average_deviation);
                }
                // OCCT L1525: BB.UpdateEdge(NE, C3d, max_deviation).
                builder_update_edge_curve(ne, c3d.clone(), max_deviation);
                // OCCT L1527-1593: the second representation pass.
                let mut projection_success = true;
                if nb_pcurves > 1 {
                    if let TShape::Edge(ed) = Arc::make_mut(&mut ne.data) {
                        for i_rep in 0..ed.representations.len() {
                            let is_curve_on_surface = matches!(
                                ed.representations[i_rep],
                                CurveRepresentation::CurveOnSurface { .. }
                                    | CurveRepresentation::CurveOnClosedSurface { .. }
                            );
                            if !is_curve_on_surface {
                                continue;
                            }
                            let the_curve: Curve2d = match &ed.representations[i_rep] {
                                CurveRepresentation::CurveOnSurface { pcurve, .. } => pcurve.clone(),
                                CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => pcurve1.clone(),
                                _ => unreachable!(),
                            };
                            // OCCT L1539-1540: the surface read — GAP
                            // (architecture difference #32).
                            let the_surf = curve_rep_surface(&ed.representations[i_rep]);
                            // OCCT L1540: theLoc = CurveRep->Location() —
                            // the identity-location convention.
                            let the_loc = curve_rep_location(&ed.representations[i_rep]);
                            // OCCT L1541: the handle identity.
                            if Some(i_rep) == min_rep_index && the_loc == min_loc {
                                continue;
                            }
                            // OCCT L1545-1546: BRep_GCurve First()/Last().
                            let (first_par, last_par) = match &ed.representations[i_rep] {
                                CurveRepresentation::CurveOnSurface { range, .. } => (range[0], range[1]),
                                CurveRepresentation::CurveOnClosedSurface { range, .. } => (range[0], range[1]),
                                _ => unreachable!(),
                            };
                            // OCCT L1547-1548.
                            if (first_par - first_par_on_pc).abs() > rcad_kernel::precision::PCONFUSION
                                || (last_par - last_par_on_pc).abs() > rcad_kernel::precision::PCONFUSION
                            {
                                // OCCT L1550-1551: theLoc/theSurf — the
                                // transformed surface — GAP.
                                let the_surf = geom_surface_transformed_loc(&the_surf);
                                // OCCT L1553-1572.
                                if let Curve2d::Line(the_line) = &the_curve {
                                    if matches!(
                                        the_surf,
                                        Surface3::BSpline(_) | Surface3::Bezier(_) | Surface3::Trimmed(_)
                                    ) {
                                        let the_dir = the_line.direction;
                                        if gp_dir2d_is_parallel(the_dir, DVec2::X, rcad_kernel::precision::ANGULAR)
                                            || gp_dir2d_is_parallel(the_dir, DVec2::Y, rcad_kernel::precision::ANGULAR)
                                        {
                                            // OCCT L1560-1561: the bounds.
                                            let bounds = the_surf.default_domain();
                                            let (u1, u2, v1, v2) = (bounds[0], bounds[1], bounds[2], bounds[3]);
                                            let origin = the_line.origin;
                                            if (origin.x - u1).abs() <= rcad_kernel::precision::CONFUSION
                                                || (origin.x - u2).abs() <= rcad_kernel::precision::CONFUSION
                                                || (origin.y - v1).abs() <= rcad_kernel::precision::CONFUSION
                                                || (origin.y - v2).abs() <= rcad_kernel::precision::CONFUSION
                                            {
                                                // OCCT L1568-1569.
                                                brep_lib_same_parameter(ne, rcad_kernel::precision::CONFUSION);
                                                break;
                                            }
                                        }
                                    }
                                }
                                // OCCT L1573-1589.
                                if c3d.is_some() && first_par_on_pc < last_par_on_pc {
                                    let proj_pcurve = geom_proj_lib_curve2d(
                                        c3d.as_ref().expect("C3d"),
                                        first_par_on_pc,
                                        last_par_on_pc,
                                        &the_surf,
                                    );
                                    match proj_pcurve {
                                        None => projection_success = false,
                                        Some(proj) => {
                                            match &mut ed.representations[i_rep] {
                                                CurveRepresentation::CurveOnSurface { pcurve, .. } => {
                                                    *pcurve = proj;
                                                }
                                                CurveRepresentation::CurveOnClosedSurface { pcurve1, .. } => {
                                                    *pcurve1 = proj;
                                                }
                                                _ => {}
                                            }
                                        }
                                    }
                                } else {
                                    // OCCT L1588: return false.
                                    return false;
                                }
                            }
                        }
                    }
                }
                // OCCT L1594-1602.
                if projection_success {
                    builder_range_edge(ne, first_par_on_pc, last_par_on_pc);
                } else {
                    // OCCT L1600: BB.Range(NE, FirstParOnPC, LastParOnPC,
                    // true) — the only2d form (the rcad range write is the
                    // same leaf; the only2d flag carries no extra state).
                    builder_range_edge(ne, first_par_on_pc, last_par_on_pc);
                    brep_lib_same_parameter(ne, rcad_kernel::precision::CONFUSION);
                }
            }
        } else {
            // OCCT L1605-1672: no pcurves.
            let c3d = match &c3d {
                Some(c) => c.clone(),
                None => panic!("ExtentEdge: C3d is null (the OCCT would dereference the null handle)"),
            };
            let mut first_par = c3d.default_domain()[0];
            let mut last_par = c3d.default_domain()[1];

            // OCCT L1610-1658.
            let is_bounded = matches!(
                c3d,
                Curve3::BSpline(_) | Curve3::Bezier(_) | Curve3::Trimmed(_)
            );
            if is_bounded && (first_par > an_ef - a2_offset || last_par < an_el + a2_offset) {
                // OCCT L1613: the trimmed copy.
                let a_tr_curve = Curve3::Trimmed(TrimmedCurve3::new(c3d.clone(), first_par, last_par));
                // OCCT L1623: the comp-curve engine — GAP (architecture
                // difference #27).
                let mut a_comp_curve = GeomConvertCompCurveToBSplineCurve::new(&a_tr_curve);
                let a_tol = rcad_kernel::precision::CONFUSION;
                let a_delta = a2_offset.max(1.);

                // OCCT L1627-1639: the begin prolongation.
                if first_par > an_ef - a2_offset {
                    let a_p_bnd = c3d.point_at(first_par);
                    let a_v_bnd = c3d.derivative_at(first_par);
                    let a_d_bnd = a_v_bnd.normalize_or_zero();
                    let a_p_beg = a_p_bnd - a_delta * a_d_bnd;
                    let a_lin = Curve3::Line(Line3::new(a_p_beg, a_d_bnd));
                    let a_segment = Curve3::Trimmed(TrimmedCurve3::new(a_lin, 0., a_delta));
                    if !a_comp_curve.add(&a_segment, a_tol) {
                        // OCCT L1637: return true.
                        return true;
                    }
                }

                // OCCT L1641-1652: the end prolongation.
                if last_par < an_el + a2_offset {
                    let a_p_beg = c3d.point_at(last_par);
                    let a_v_bnd = c3d.derivative_at(last_par);
                    let a_d_bnd = a_v_bnd.normalize_or_zero();
                    let a_lin = Curve3::Line(Line3::new(a_p_beg, a_d_bnd));
                    let a_segment = Curve3::Trimmed(TrimmedCurve3::new(a_lin, 0., a_delta));
                    if !a_comp_curve.add(&a_segment, a_tol) {
                        // OCCT L1650: return true.
                        return true;
                    }
                }

                // OCCT L1654-1657.
                let c3d_new = a_comp_curve.bspline_curve();
                first_par = c3d_new.default_domain()[0];
                last_par = c3d_new.default_domain()[1];
                builder_update_edge_curve(ne, Some(c3d_new), rcad_kernel::precision::CONFUSION);
            } else if c3d.is_periodic() {
                // OCCT L1659-1665.
                let mut delta = (curve3_period(&c3d) - (an_el - an_ef)) * 0.5;
                delta *= 0.95;
                first_par = an_ef - delta;
                last_par = an_el + delta;
            } else if c3d.is_closed() {
                // OCCT L1666-1669.
                last_par -= 0.05 * (last_par - first_par);
            }

            // OCCT L1671.
            builder_range_edge(ne, first_par, last_par);
        }
        // OCCT L1673.
        true
    }

    /// OCCT BRepOffset_Inter2d::Compute(AsDes, F, NewEdges, Tol,
    /// theEdgeIntEdges, theDMVV, theRange) (BRepOffset_Inter2d.cxx
    /// L1725-1826).
    pub fn compute(
        as_des: &mut BRepAlgoAsDes,
        f: &Shape,
        new_edges: &IndexedShapeMap,
        tol: f64,
        the_edge_int_edges: &HashMap<ShapeKey, Vec<Shape>>,
        the_dmvv: &mut DmvvMap,
        the_range: (),
    ) {
        // OCCT L1738-1743: the edges of the face (NCollection_Map).
        let mut edges_of_face: HashSet<ShapeKey> = HashSet::new();
        for e in explorer(f, ShapeType::Edge, ShapeType::Shape) {
            edges_of_face.insert(shape_key(&e));
        }

        // OCCT L1755: LE = AsDes->Descendant(F) — the rcad borrow-split
        // clone (the AsDes is mutated inside the loop).
        let le = as_des.descendant(f).to_vec();

        // OCCT L1758: BRepAdaptor_Surface BAsurf(F).
        let b_asurf = BRepAdaptorSurface::new(f);
        // OCCT L1760-1766: the Message_ProgressScope guards — the rcad () range
        // (architecture difference #33).
        let _ = the_range;

        // OCCT L1757: int j, i = 1.
        let mut i = 1usize;
        for i1 in 0..le.len() {
            // OCCT L1767: const TopoDS_Edge& E1 = TopoDS::Edge(it1LE.Value()).
            let e1 = &le[i1];
            let mut j = 1usize;
            let mut j1 = 0usize;
            // OCCT L1771: while (j < i && it2LE.More()).
            while j < i && j1 < le.len() {
                let e2 = &le[j1];

                // OCCT L1775-1807: the theEdgeIntEdges filter.
                let mut to_intersect = true;
                let k1 = shape_key(e1);
                if let Some(a_elist) = the_edge_int_edges.get(&k1) {
                    for itedges in a_elist {
                        if e2.is_same(itedges) {
                            to_intersect = false;
                        }
                    }

                    if to_intersect {
                        for an_edge in a_elist {
                            if let Some(a_elist2) = the_edge_int_edges.get(&shape_key(an_edge)) {
                                for itedges2 in a_elist2 {
                                    if e2.is_same(itedges2) {
                                        to_intersect = false;
                                    }
                                }
                            }
                        }
                    }
                }

                // OCCT L1813-1820: the intersection condition.
                let k2 = shape_key(e2);
                if to_intersect
                    && (!edges_of_face.contains(&k1) || !edges_of_face.contains(&k2))
                    && (new_edges.contains_key(&k1) || new_edges.contains_key(&k2))
                {
                    // OCCT L1817-1818: EdgeInter(F.Oriented(FORWARD), ...).
                    let a_local_shape = oriented(f, Orientation::Forward);
                    edge_inter(
                        &a_local_shape,
                        &b_asurf,
                        e1,
                        e2,
                        as_des,
                        tol,
                        true,
                        the_dmvv,
                    );
                }
                // OCCT L1821-1822: it2LE.Next(); j++.
                j1 += 1;
                j += 1;
            }
            // OCCT L1824: i++.
            i += 1;
        }
    }

    /// OCCT BRepOffset_Inter2d::ConnexIntByInt(FI, OFI, MES, Build, theAsDes,
    /// AsDes2d, Offset, Tol, Analyse, FacesWithVerts, theImageVV,
    /// theEdgeIntEdges, theDMVV, theRange) (BRepOffset_Inter2d.cxx
    /// L1830-2094).
    pub fn connex_int_by_int(
        fi: &Shape,
        ofi: &mut BRepOffsetOffset,
        mes: &mut HashMap<ShapeKey, Shape>,
        build: &HashMap<ShapeKey, Shape>,
        the_as_des: &BRepAlgoAsDes,
        as_des2d: &mut BRepAlgoAsDes,
        offset: f64,
        tol: f64,
        analyse: &BRepOffsetAnalyse,
        faces_with_verts: &mut IndexedShapeMap,
        the_image_vv: &mut BRepAlgoImage,
        the_edge_int_edges: &mut HashMap<ShapeKey, Vec<Shape>>,
        the_dmvv: &mut DmvvMap,
        the_range: (),
    ) -> bool {
        // OCCT L1851-1854: the Message_ProgressScope — the rcad () range
        // (architecture difference #33).
        let _ = the_range;

        // OCCT L1849-1850: MVE = MapVertexEdges(FI) — GAP leaf
        // (architecture difference #23).
        let mut mve: HashMap<ShapeKey, Vec<Shape>> = HashMap::new();
        brep_offset_tool_map_vertex_edges(fi, &mut mve);

        //---------------------
        // OCCT L1855-1899: Extension of edges.
        //---------------------
        let mut ne = Shape::null();
        for (_, l) in mve.iter() {
            // OCCT L1868-1877: YaBuild.
            let mut ya_build = false;
            for itl in l {
                ya_build = build.contains_key(&shape_key(itl));
                if ya_build {
                    break;
                }
            }
            if ya_build {
                for itl in l {
                    let ei = itl;
                    // OCCT L1883-1886: the F/R orientation guard.
                    if ei.orientation != Orientation::Forward && ei.orientation != Orientation::Reversed {
                        continue;
                    }
                    // OCCT L1887-1888: OE = OFI.Generated(EI).
                    let a_local_shape = ofi.generated(ei);
                    let oe = &a_local_shape;
                    // OCCT L1889-1896.
                    if !mes.contains_key(&shape_key(oe)) && !build.contains_key(&shape_key(ei)) {
                        if !Self::extent_edge(oe, &mut ne, offset) {
                            return false;
                        }
                        mes.insert(shape_key(oe), ne.clone());
                    }
                }
            }
        }

        // OCCT L1901-1905: FIO.
        let mut fio = ofi.face();
        if let Some(v) = mes.get(&shape_key(&fio)) {
            fio = v.clone();
        }

        // OCCT L1907: BRepAdaptor_Surface BAsurf(FIO).
        let b_asurf = BRepAdaptorSurface::new(&fio);

        // OCCT L1909-2092: the wire walk.
        let fi_forward = oriented(fi, Orientation::Forward);
        for exp_current in explorer(&fi_forward, ShapeType::Wire, ShapeType::Shape) {
            let w = &exp_current;
            let mut wexp = BRepToolsWireExplorer::new();
            let mut end = false;
            // OCCT L1921-1923: the FORWARD-oriented wire/face.
            let a_local_wire = oriented(w, Orientation::Forward);
            let a_local_face = oriented(fi, Orientation::Forward);
            wexp.init(&a_local_wire, &a_local_face);
            if !wexp.more() {
                // OCCT L1926: the empty-wire protection.
                continue;
            }
            let first_e = wexp.current();
            let mut cur_e = first_e.clone();
            // OCCT L1929: the (unused in the OCCT source) Edges indexed map.
            let _edges: IndexedShapeMap = IndexedShapeMap::new();

            while !end {
                wexp.next();
                let mut next_e = if wexp.more() {
                    wexp.current()
                } else {
                    // OCCT L1940-1941.
                    end = true;
                    first_e.clone()
                };
                if cur_e.is_same(&next_e) {
                    continue;
                }

                // OCCT L1948: Vref = CommonVertex(CurE, NextE).
                let vref = common_vertex(&cur_e, &next_e);

                // OCCT L1950-1951: the Analyse replacement (GAP carrier).
                cur_e = analyse.edge_replacement(fi, &cur_e);
                next_e = analyse.edge_replacement(fi, &next_e);

                // OCCT L1953-1956: CEO / NEO.
                let ceo = ofi.generated(&cur_e);
                let neo = ofi.generated(&next_e);

                //------------------------------------------
                // OCCT L1960-2002: the image selection.
                //------------------------------------------
                // OCCT L1960: LV1, LV2 — declared but unused in the OCCT
                // source.
                let _lv1: Vec<Shape> = Vec::new();
                let _lv2: Vec<Shape> = Vec::new();
                let mut do_inter = true;
                let mut an_or1 = Orientation::External;
                let mut an_or2 = Orientation::External;

                let mut a_choice = 0i32;
                let mut ne1 = Shape::null();
                let mut ne2 = Shape::null();
                let _ = (&ne1, &ne2); // the OCCT default nulls (L1962).
                let mut ne1_seq: Vec<Shape> = Vec::new();
                let mut ne2_seq: Vec<Shape> = Vec::new();
                if build.contains_key(&shape_key(&cur_e)) && build.contains_key(&shape_key(&next_e)) {
                    // OCCT L1967-1976: aChoice = 1.
                    a_choice = 1;
                    ne1 = build[&shape_key(&cur_e)].clone();
                    ne2 = build[&shape_key(&next_e)].clone();
                    get_edges_oriented_in_face(&ne1, &fio, the_as_des, &mut ne1_seq);
                    get_edges_oriented_in_face(&ne2, &fio, the_as_des, &mut ne2_seq);
                    an_or1 = Orientation::Reversed;
                    an_or2 = Orientation::Forward;
                } else if build.contains_key(&shape_key(&cur_e)) && mes.contains_key(&shape_key(&neo)) {
                    // OCCT L1977-1987: aChoice = 2.
                    a_choice = 2;
                    ne1 = build[&shape_key(&cur_e)].clone();
                    ne2 = mes[&shape_key(&neo)].clone();
                    ne2.orientation = next_e.orientation;
                    get_edges_oriented_in_face(&ne1, &fio, the_as_des, &mut ne1_seq);
                    ne2_seq.push(ne2.clone());
                    an_or1 = Orientation::Reversed;
                    an_or2 = Orientation::Forward;
                } else if build.contains_key(&shape_key(&next_e)) && mes.contains_key(&shape_key(&ceo)) {
                    // OCCT L1988-1998: aChoice = 3.
                    a_choice = 3;
                    ne1 = build[&shape_key(&next_e)].clone();
                    ne2 = mes[&shape_key(&ceo)].clone();
                    ne2.orientation = cur_e.orientation;
                    get_edges_oriented_in_face(&ne1, &fio, the_as_des, &mut ne1_seq);
                    ne2_seq.push(ne2.clone());
                    an_or1 = Orientation::Forward;
                    an_or2 = Orientation::Reversed;
                } else {
                    do_inter = false;
                }
                if do_inter {
                    //------------------------------------
                    // OCCT L2005-2028: NE1/NE2 can be a compound of edges.
                    //------------------------------------
                    let mut b_coincide = false;
                    let (a_e1, a_e2) = if a_choice == 1 || a_choice == 2 {
                        (ne1_seq[ne1_seq.len() - 1].clone(), ne2_seq[0].clone())
                    } else {
                        (ne1_seq[0].clone(), ne2_seq[ne2_seq.len() - 1].clone())
                    };

                    if a_e1.orientation == Orientation::Reversed {
                        an_or1 = top_abs_reverse(an_or1);
                    }
                    if a_e2.orientation == Orientation::Reversed {
                        an_or2 = top_abs_reverse(an_or2);
                    }

                    // OCCT L2030-2042.
                    ref_edge_inter(
                        &fio,
                        &b_asurf,
                        &a_e1,
                        &a_e2,
                        an_or1,
                        an_or2,
                        as_des2d,
                        tol,
                        true,
                        &vref,
                        the_image_vv,
                        the_dmvv,
                        &mut b_coincide,
                    );

                    // OCCT L2044-2063: the symmetric theEdgeIntEdges record.
                    let k_e1 = shape_key(&a_e1);
                    let k_e2 = shape_key(&a_e2);
                    match the_edge_int_edges.get_mut(&k_e1) {
                        Some(list) => list.push(a_e2.clone()),
                        None => {
                            the_edge_int_edges.insert(k_e1, vec![a_e2.clone()]);
                        }
                    }
                    match the_edge_int_edges.get_mut(&k_e2) {
                        Some(list) => list.push(a_e1.clone()),
                        None => {
                            the_edge_int_edges.insert(k_e2, vec![a_e1.clone()]);
                        }
                    }

                    // OCCT L2066-2071.
                    if build.contains_key(&shape_key(&vref)) {
                        faces_with_verts.insert(shape_key(fi), fi.clone());
                    }
                } else {
                    // OCCT L2073-2089.
                    let v = common_vertex(&ceo, &neo);
                    if !v.is_null() {
                        // OCCT L2078-2082 (the MES edge is the shared TShape;
                        // the rcad mutation targets the local map copy —
                        // architecture difference #22).
                        if mes.contains_key(&shape_key(&ceo)) {
                            let mut oe = mes[&shape_key(&ceo)].clone();
                            update_vertex(&v, &ceo, &mut oe, tol);
                            as_des2d.add(&oe, &v);
                        }
                        // OCCT L2083-2087.
                        if mes.contains_key(&shape_key(&neo)) {
                            let mut oe = mes[&shape_key(&neo)].clone();
                            update_vertex(&v, &neo, &mut oe, tol);
                            as_des2d.add(&oe, &v);
                        }
                    }
                }
                // OCCT L2090: CurE = wexp.Current().
                cur_e = wexp.current();
            }
        }
        // OCCT L2093.
        true
    }

    /// OCCT BRepOffset_Inter2d::ConnexIntByIntInVert(FI, OFI, MES, Build,
    /// AsDes, AsDes2d, Tol, Analyse, theDMVV, theRange)
    /// (BRepOffset_Inter2d.cxx L2100-2307).
    pub fn connex_int_by_int_in_vert(
        fi: &Shape,
        ofi: &mut BRepOffsetOffset,
        mes: &mut HashMap<ShapeKey, Shape>,
        build: &HashMap<ShapeKey, Shape>,
        as_des: &BRepAlgoAsDes,
        as_des2d: &mut BRepAlgoAsDes,
        tol: f64,
        analyse: &BRepOffsetAnalyse,
        the_dmvv: &mut DmvvMap,
        the_range: (),
    ) {
        // OCCT L2110: theRange (architecture difference #33).
        let _ = the_range;

        // OCCT L2113-2117: FIO.
        let mut fio = ofi.face();
        if let Some(v) = mes.get(&shape_key(&fio)) {
            fio = v.clone();
        }

        // OCCT L2119-2126: aME — the descendants of FIO.
        let mut a_me: HashSet<ShapeKey> = HashSet::new();
        let a_le = as_des.descendant(&fio).to_vec();
        for a_e in &a_le {
            a_me.insert(shape_key(a_e));
        }

        // OCCT L2128: BRepAdaptor_Surface BAsurf(FIO).
        let b_asurf = BRepAdaptorSurface::new(&fio);

        // OCCT L2131-2306: the wire walk.
        let fi_forward = oriented(fi, Orientation::Forward);
        for exp_current in explorer(&fi_forward, ShapeType::Wire, ShapeType::Shape) {
            let w = &exp_current;
            let mut wexp = BRepToolsWireExplorer::new();
            let mut end = false;
            // OCCT L2144-2146.
            let a_local_wire = oriented(w, Orientation::Forward);
            let a_local_face = oriented(fi, Orientation::Forward);
            wexp.init(&a_local_wire, &a_local_face);
            if !wexp.more() {
                // OCCT L2149: the empty-wire protection.
                continue;
            }

            let first_e = wexp.current();
            let mut cur_e = first_e.clone();
            while !end {
                wexp.next();
                let mut next_e = if wexp.more() {
                    wexp.current()
                } else {
                    end = true;
                    first_e.clone()
                };
                if cur_e.is_same(&next_e) {
                    continue;
                }

                // OCCT L2170-2175.
                let vref = common_vertex(&cur_e, &next_e);
                if !build.contains_key(&shape_key(&vref)) {
                    cur_e = next_e;
                    continue;
                }

                // OCCT L2177-2178.
                cur_e = analyse.edge_replacement(fi, &cur_e);
                next_e = analyse.edge_replacement(fi, &next_e);

                // OCCT L2180-2183: CEO / NEO.
                let ceo = ofi.generated(&cur_e);
                let neo = ofi.generated(&next_e);

                // OCCT L2185-2207: the image selection.
                let an_or1 = Orientation::External;
                let an_or2 = Orientation::External;
                let ne1: Shape;
                let ne2: Shape;
                if build.contains_key(&shape_key(&cur_e)) && build.contains_key(&shape_key(&next_e)) {
                    ne1 = build[&shape_key(&cur_e)].clone();
                    ne2 = build[&shape_key(&next_e)].clone();
                } else if build.contains_key(&shape_key(&cur_e)) && mes.contains_key(&shape_key(&neo)) {
                    ne1 = build[&shape_key(&cur_e)].clone();
                    ne2 = mes[&shape_key(&neo)].clone();
                } else if build.contains_key(&shape_key(&next_e)) && mes.contains_key(&shape_key(&ceo)) {
                    ne1 = build[&shape_key(&next_e)].clone();
                    ne2 = mes[&shape_key(&ceo)].clone();
                } else {
                    // OCCT L2203-2207.
                    cur_e = wexp.current();
                    continue;
                }

                // OCCT L2212: NE3 = Build(Vref).
                let ne3 = build[&shape_key(&vref)].clone();

                // OCCT L2214-2303.
                for exp2_current in explorer(&ne3, ShapeType::Edge, ShapeType::Shape) {
                    let a_e3 = &exp2_current;
                    if !a_me.contains(&shape_key(a_e3)) {
                        continue;
                    }

                    // OCCT L2222-2245: the intersection with the first edge.
                    for exp1_current in explorer(&ne1, ShapeType::Edge, ShapeType::Shape) {
                        let a_e1 = &exp1_current;
                        let mut an_empty_image = BRepAlgoImage::new();
                        let mut b_coincide = false;
                        ref_edge_inter(
                            &fio,
                            &b_asurf,
                            a_e1,
                            a_e3,
                            an_or1,
                            an_or2,
                            as_des2d,
                            tol,
                            true,
                            &vref,
                            &mut an_empty_image,
                            the_dmvv,
                            &mut b_coincide,
                        );
                        if b_coincide {
                            // OCCT L2242-2243: trim the edge E3 the same way
                            // as E1 (the borrow-split clone of the
                            // descendants).
                            let mut a_e3_mut = a_e3.clone();
                            let mut desc = as_des2d.descendant(a_e1).to_vec();
                            store(&mut a_e3_mut, &mut desc, tol, true, as_des2d, the_dmvv);
                        }
                    }

                    // OCCT L2247-2270: the intersection with the second edge.
                    for exp1_current in explorer(&ne2, ShapeType::Edge, ShapeType::Shape) {
                        let a_e2 = &exp1_current;
                        let mut an_empty_image = BRepAlgoImage::new();
                        let mut b_coincide = false;
                        ref_edge_inter(
                            &fio,
                            &b_asurf,
                            a_e2,
                            a_e3,
                            an_or1,
                            an_or2,
                            as_des2d,
                            tol,
                            true,
                            &vref,
                            &mut an_empty_image,
                            the_dmvv,
                            &mut b_coincide,
                        );
                        if b_coincide {
                            // OCCT L2267-2268: trim E3 the same way as E2.
                            let mut a_e3_mut = a_e3.clone();
                            let mut desc = as_des2d.descendant(a_e2).to_vec();
                            store(&mut a_e3_mut, &mut desc, tol, true, as_des2d, the_dmvv);
                        }
                    }

                    // OCCT L2272-2302: the edges generated from the vertex
                    // among themselves.
                    let ne3_edges = explorer(&ne3, ShapeType::Edge, ShapeType::Shape);
                    // OCCT L2274-2280: skip the edges up to aE3.
                    let mut start = ne3_edges.len();
                    for (idx, e) in ne3_edges.iter().enumerate() {
                        if a_e3.is_same(e) {
                            start = idx + 1;
                            break;
                        }
                    }
                    for a_e3_next in ne3_edges.iter().skip(start) {
                        if a_me.contains(&shape_key(a_e3_next)) {
                            let mut an_empty_image = BRepAlgoImage::new();
                            let mut b_coincide = false;
                            ref_edge_inter(
                                &fio,
                                &b_asurf,
                                a_e3_next,
                                a_e3,
                                an_or1,
                                an_or2,
                                as_des2d,
                                tol,
                                true,
                                &vref,
                                &mut an_empty_image,
                                the_dmvv,
                                &mut b_coincide,
                            );
                        }
                    }
                }
                // OCCT L2304: CurE = wexp.Current().
                cur_e = wexp.current();
            }
        }
    }

    /// OCCT BRepOffset_Inter2d::FuseVertices(theDMVV, theAsDes, theImageVV)
    /// (BRepOffset_Inter2d.cxx L2335-2394).
    pub fn fuse_vertices(
        the_dmvv: &DmvvMap,
        the_as_des: &mut BRepAlgoAsDes,
        the_image_vv: &mut BRepAlgoImage,
    ) -> bool {
        // OCCT L2343.
        let mut a_mv_done: HashSet<ShapeKey> = HashSet::new();
        // OCCT L2344-2345: for i in 1..=theDMVV.Extent().
        for i in 1..=the_dmvv.len() {
            // OCCT L2347: aV = FindKey(i).
            let a_v = the_dmvv[i - 1].0.clone();

            // OCCT L2349-2351: the chain of vertices.
            let mut a_lv_chain: Vec<Shape> = Vec::new();
            make_chain(&a_v, the_dmvv, &mut a_mv_done, &mut a_lv_chain);

            // OCCT L2353-2356.
            if a_lv_chain.len() < 2 {
                continue;
            }

            // OCCT L2358-2362: the new vertex.
            let a_v_new = bop_tools_algo_tools_make_vertex(&a_lv_chain);
            let mut a_v_new_int = oriented(&a_v_new, Orientation::Internal);

            // OCCT L2364-2391.
            for a_v_old in &a_lv_chain {
                // OCCT L2369: aVOldInt = aVOld.Oriented(INTERNAL).
                let a_v_old_int = oriented(a_v_old, Orientation::Internal);
                // OCCT L2370: the ascendants (the rcad borrow-split clone;
                // the AsDes is replaced inside the loop).
                let mut a_le = the_as_des.ascendant(a_v_old).to_vec();

                for a_e in a_le.iter_mut() {
                    // OCCT L2375-2376.
                    let a_tol_e = brep_tool_tolerance(a_e);
                    // OCCT L2378-2381: BRep_Tool::Parameter(aVOldInt, aE, aT).
                    let a_t = match brep_tool_parameter_bool(&a_v_old_int, a_e) {
                        Some(t) => t,
                        None => return false,
                    };
                    // OCCT L2382: the ascendant edge is the shared TShape;
                    // the rcad mutation targets the local clone copy
                    // (architecture difference #22).
                    builder_update_vertex_on_edge(&mut a_v_new_int, a_t, a_e, a_tol_e);
                }
                // OCCT L2384-2385.
                the_as_des.replace(a_v_old, &a_v_new);
                // OCCT L2386-2390.
                if the_image_vv.is_image(a_v_old) {
                    let a_pro_vertex = the_image_vv.image_from(a_v_old).clone();
                    the_image_vv.add(&a_pro_vertex, &oriented(&a_v_new, Orientation::Forward));
                }
            }
        }
        // OCCT L2393.
        true
    }
}

/// OCCT UpdateVertex(V, OE, NE, TolConf) (BRepOffset_Inter2d.cxx L1680-1721)
/// — the file static.
fn update_vertex(v: &Shape, oe: &Shape, ne: &mut Shape, tol_conf: f64) -> bool {
    // OCCT L1682-1683: the adaptors.
    let oc = BRepAdaptorCurve::new(oe, oe);
    let nc = BRepAdaptorCurve::new(ne, ne);
    // OCCT L1684-1691.
    let of = oc.first_parameter();
    let ol = oc.last_parameter();
    let nf = nc.first_parameter();
    let nl = nc.last_parameter();
    let mut u = 0.0f64;
    let par_tol = rcad_kernel::precision::PCONFUSION;
    let p = brep_tool_pnt(v).unwrap_or(DVec3::ZERO);
    let mut ok = false;

    // OCCT L1693-1700.
    if p.distance(oc.value(of)) < tol_conf {
        if of >= nf + par_tol && of <= nl + par_tol && p.distance(nc.value(of)) < tol_conf {
            ok = true;
            u = of;
        }
    }
    // OCCT L1701-1708.
    if p.distance(oc.value(ol)) < tol_conf {
        if ol >= nf + par_tol && ol <= nl + par_tol && p.distance(nc.value(ol)) < tol_conf {
            ok = true;
            u = ol;
        }
    }
    // OCCT L1709-1719.
    if ok {
        // OCCT L1712-1716: EE = NE.Oriented(FORWARD) — the orientation copy
        // shares the TShape; the rcad mutation targets ne directly.
        // UpdateVertex(V.Oriented(INTERNAL), U, NE, Tolerance(NE)).
        let mut v = v.clone();
        let tol_ne = brep_tool_tolerance(ne);
        builder_update_vertex_on_edge(&mut v, u, ne, tol_ne);
    }
    // OCCT L1720.
    ok
}

/// OCCT MakeChain(theV, theDMVV, theMDone, theChain)
/// (BRepOffset_Inter2d.cxx L2311-2331) — the file static.
fn make_chain(
    the_v: &Shape,
    the_dmvv: &DmvvMap,
    the_m_done: &mut HashSet<ShapeKey>,
    the_chain: &mut Vec<Shape>,
) {
    // OCCT L2318: theMDone.Add — false when already done.
    if the_m_done.insert(shape_key(the_v)) {
        the_chain.push(the_v.clone());
        // OCCT L2321-2329.
        if let Some((_, p_lv)) = the_dmvv.get(&shape_key(the_v)) {
            for a_it in p_lv {
                make_chain(a_it, the_dmvv, the_m_done, the_chain);
            }
        }
    }
}
