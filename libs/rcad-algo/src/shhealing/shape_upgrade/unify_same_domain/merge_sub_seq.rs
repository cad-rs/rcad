//! OCCT ShapeUpgrade_UnifySameDomain.cxx L2164-2675 — `MergeSubSeq`
//! (L2164-2501), `IsMergingPossible` (L2508-2638), `GetLineEdgePoints`
//! (L2642-2675).
//!
//! Signature note: the OCCT const `BRep&` parameters become `&mut BRep` in
//! the rcad bridge — the shared explorer re-host registers composed
//! locations in the BRep table (the brep_tool.rs contract).

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Circle3, CurveEval, Line3};
use rcad_kernel::precision::{p_confusion, square_p_confusion, CONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation, ShapeType, TShape};
use std::sync::Arc;

use super::gap_deps::{brep_lib_make_edge, BRepAdaptorCurve, GeomAbsCurveType};
use super::statics_a::is_linear;
use super::statics_b::find_closest_points;
use super::topexp::{first_vertex, is_edge_degenerated, last_vertex, occt_is_same_shape};
use super::{occt_reverse, IndexedDataMapOfShapeListOfShape, ShapeUpgradeUnifySameDomain};
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;

impl ShapeUpgradeUnifySameDomain {
    // OCCT ShapeUpgrade_UnifySameDomain.cxx L2164-2501: MergeSubSeq —
    // merges a sequence of edges into one edge if possible.
    pub(crate) fn merge_sub_seq(
        &mut self,
        brep: &mut BRep,
        the_chain: &[Shape],
        the_vfmap: &IndexedDataMapOfShapeListOfShape,
        out_edge: &mut Shape,
    ) -> bool {
        // OCCT L2171-2172.
        let sae = ShapeAnalysisEdge::new();
        let mut b = BRepBuilder::new();
        // OCCT L2174-2178.
        let mut is_union_of_lines_possible = true;
        let mut is_union_of_circles_possible = true;
        for j in 1..the_chain.len() {
            let edge1 = the_chain[j - 1].clone();
            let edge2 = the_chain[j].clone();

            // OCCT L2184-2240: the double-degenerated 2D construction.
            if is_edge_degenerated(brep, &edge1) && is_edge_degenerated(brep, &edge2) {
                // OCCT L2186-2207.
                let edge_first = the_chain[0].clone();
                let edge_last = the_chain[the_chain.len() - 1].clone();
                let mut common_face = Shape::null();
                let mut min_sq_dist = 0.0f64;
                let mut or_of_e1 = Orientation::Forward;
                let mut or_of_e2 = Orientation::Forward;
                let mut ind_on_e1 = 0i32;
                let mut ind_on_e2 = 0i32;
                let mut points_on_edge1 = [DVec2::ZERO; 2];
                let mut points_on_edge2 = [DVec2::ZERO; 2];
                if !find_closest_points(
                    brep,
                    &edge_first,
                    &edge_last,
                    the_vfmap,
                    &mut common_face,
                    &mut min_sq_dist,
                    &mut or_of_e1,
                    &mut or_of_e2,
                    &mut ind_on_e1,
                    &mut ind_on_e2,
                    &mut points_on_edge1,
                    &mut points_on_edge2,
                ) {
                    return false;
                }

                // OCCT L2209-2211: the extremities of the future edge.
                let ind_on_e1 = 1 - ind_on_e1;
                let ind_on_e2 = 1 - ind_on_e2;

                // OCCT L2213-2222.
                let mut start_point = points_on_edge1[ind_on_e1 as usize];
                let mut end_point = points_on_edge2[ind_on_e2 as usize];
                if (or_of_e1 == Orientation::Forward && ind_on_e1 == 1)
                    || (or_of_e1 == Orientation::Reversed && ind_on_e1 == 0)
                {
                    std::mem::swap(&mut start_point, &mut end_point);
                }

                // OCCT L2224: GC_MakeLine2d(StartPoint, EndPoint).
                let a_line = super::gap_deps::gc_make_line2d(start_point, end_point);

                // OCCT L2226-2229.
                let a_vertex = first_vertex(brep, &edge_first, false);
                let mut start_vertex = a_vertex.clone();
                start_vertex.orientation = Orientation::Forward;
                let mut end_vertex = a_vertex.clone();
                end_vertex.orientation = Orientation::Reversed;

                // OCCT L2231-2237.
                let new_edge = make_empty_edge(brep);
                b.update_edge_pcurve(
                    brep,
                    new_edge.clone(),
                    rcad_kernel::geom::Curve2d::Line(a_line),
                    common_face.clone(),
                    CONFUSION,
                );
                b.set_edge_range(brep, new_edge.clone(), 0.0, start_point.distance(end_point));
                b.add_to_edge(brep, new_edge.clone(), start_vertex);
                b.add_to_edge(brep, new_edge.clone(), end_vertex);
                b.set_edge_degenerated(brep, new_edge.clone(), true);
                *out_edge = new_edge;
                return true;
            }

            // OCCT L2242-2252: the curves (the trimmed bases are already
            // unwrapped by the BRep_Tool::Curve re-host; the range applied).
            let c3d1 = super::topexp::brep_tool_curve(brep, &edge1);
            let c3d2 = super::topexp::brep_tool_curve(brep, &edge2);
            let a_ba_curve1 = BRepAdaptorCurve::new(brep, &edge1);
            let a_ba_curve2 = BRepAdaptorCurve::new(brep, &edge2);

            // OCCT L2249-2252.
            let (Some(_c3d1), Some(_c3d2)) = (c3d1, c3d2) else {
                return false;
            };

            // OCCT L2264-2276.
            let mut a_dir1 = DVec3::ONE;
            let mut a_dir2 = DVec3::ONE;
            let lin1 = is_linear(&a_ba_curve1, &mut a_dir1);
            let lin2 = is_linear(&a_ba_curve2, &mut a_dir2);
            if lin1 && lin2 {
                if !super::gap_deps::dir_is_parallel_3d(a_dir1, a_dir2, self.my_ang_tol) {
                    is_union_of_lines_possible = false;
                }
            } else {
                is_union_of_lines_possible = false;
            }

            // OCCT L2277-2291.
            let circle_pair = matches!(
                (a_ba_curve1.get_type(), a_ba_curve2.get_type()),
                (GeomAbsCurveType::Circle, GeomAbsCurveType::Circle)
            );
            if circle_pair {
                if let (Some(c1), Some(c2)) = (a_ba_curve1.circle(), a_ba_curve2.circle()) {
                    let p01 = c1.center;
                    let p02 = c2.center;
                    if p01.distance(p02) > CONFUSION {
                        is_union_of_circles_possible = false;
                    }
                }
            } else {
                is_union_of_circles_possible = false;
            }
        }
        // OCCT L2293-2296.
        if is_union_of_lines_possible && is_union_of_circles_possible {
            return false;
        }

        // OCCT L2298-2335: the union of lines.
        if is_union_of_lines_possible {
            let mut v = [
                sae.first_vertex(brep, &the_chain[0]),
                sae.last_vertex(brep, &the_chain[the_chain.len() - 1]),
            ];
            let pv1 = brep.vertex_position(&v[0]);
            let pv2 = brep.vertex_position(&v[1]);
            let vec = pv2 - pv1;
            // OCCT L2307-2322: the safe-input vertex copies.
            if self.my_safe_input_mode {
                for k in 0..2 {
                    if !self.my_context.is_recorded(&v[k]) {
                        let vcopy = brep.empty_copy(v[k].clone());
                        self.my_context.replace(brep, &v[k], &vcopy);
                        v[k] = vcopy;
                    } else {
                        v[k] = self.my_context.apply(brep, &v[k], ShapeType::Shape);
                    }
                }
            }
            // OCCT L2323-2325: the trimmed line over [0, dist].
            let dist = pv1.distance(pv2);
            let line = Line3 {
                origin: pv1,
                direction: vec.normalize_or_zero(),
            };
            let tc = rcad_kernel::geom::Curve3::Line(line);
            // OCCT L2326-2331.
            let e = brep_lib_make_edge(brep, &tc, &v[0], &v[1], 0.0, dist);
            b.update_vertex_on_edge(brep, v[0].clone(), 0.0, e.clone(), 0.0);
            b.update_vertex_on_edge(brep, v[1].clone(), dist, e.clone(), 0.0);
            self.union_pcurves(brep, the_chain, &e);
            *out_edge = e;
            return true;
        }

        // OCCT L2337-2462: the union of circles.
        if is_union_of_circles_possible {
            let fe = the_chain[0].clone();
            let Some((c3d, f, l)) = super::topexp::brep_tool_curve(brep, &fe) else {
                return false;
            };
            // OCCT L2343-2348: the circle basis.
            let cir = match &c3d {
                rcad_kernel::geom::Curve3::Circle(c) => Some(*c),
                _ => None,
            };
            let Some(cir) = cir else {
                return false;
            };
            let _ = (f, l);

            let mut v = [
                sae.first_vertex(brep, &fe),
                sae.last_vertex(brep, &the_chain[the_chain.len() - 1]),
            ];
            let mut is_closed = occt_is_same_shape(&v[0], &v[1]);
            if !is_closed {
                // OCCT L2356-2367: the point-equality closedness check.
                let a_p0 = brep.vertex_position(&v[0]);
                let a_p1 = brep.vertex_position(&v[1]);
                let a_tol = brep.tolerance(&v[0]).max(brep.tolerance(&v[1]));
                if a_p0.distance_squared(a_p1) < a_tol * a_tol {
                    is_closed = true;
                    v[1] = v[0].clone();
                    v[1].orientation = occt_reverse(v[1].orientation);
                }
            }
            let e;
            if is_closed {
                // OCCT L2369-2407: the closed chain.
                let a_def = BRepAdaptorCurve::new(brep, &fe);
                let fp;
                let lp;
                if fe.orientation == Orientation::Forward {
                    fp = a_def.first_parameter();
                    lp = a_def.last_parameter();
                } else {
                    fp = a_def.last_parameter();
                    lp = a_def.first_parameter();
                }
                if fp.abs() < p_confusion() {
                    // OCCT L2387-2391.
                    let ne = make_full_circle_edge(brep, &cir);
                    b.add_to_edge(brep, ne.clone(), v[0].clone());
                    b.add_to_edge(brep, ne.clone(), v[1].clone());
                    let mut ne = ne;
                    ne.orientation = fe.orientation;
                    e = ne;
                } else {
                    // OCCT L2394-2406: GC_MakeCircle over the three points.
                    let Some(cir1) = super::gap_deps::gc_make_circle_3_points(
                        a_def.value(fp),
                        a_def.value((fp + lp) * 0.5),
                        a_def.value(lp),
                    ) else {
                        return false;
                    };
                    let ne = make_full_circle_edge(brep, &cir1);
                    b.add_to_edge(brep, ne.clone(), v[0].clone());
                    b.add_to_edge(brep, ne.clone(), v[1].clone());
                    e = ne;
                }
            } else {
                // OCCT L2408-2458: the open chain.
                let param_first = brep
                    .parameter_on_edge(&v[0], &fe, &Shape::null())
                    .unwrap_or(0.0);
                let vertex_last_on_fe = sae.last_vertex(brep, &fe);
                let mut param_last = brep
                    .parameter_on_edge(&vertex_last_on_fe, &fe, &Shape::null())
                    .unwrap_or(0.0);

                // OCCT L2414-2429: the safe-input vertex copies.
                if self.my_safe_input_mode {
                    for k in 0..2 {
                        if !self.my_context.is_recorded(&v[k]) {
                            let vcopy = brep.empty_copy(v[k].clone());
                            self.my_context.replace(brep, &v[k], &vcopy);
                            v[k] = vcopy;
                        } else {
                            v[k] = self.my_context.apply(brep, &v[k], ShapeType::Shape);
                        }
                    }
                }

                let point_first = brep.vertex_position(&v[0]);
                // OCCT L2432-2435: the parametrization fold.
                while (param_last - param_first).abs() > 7.0 * std::f64::consts::PI / 8.0 {
                    param_last = (param_first + param_last) / 2.0;
                }
                let a_def = BRepAdaptorCurve::new(brep, &fe);
                let point_last = a_def.value(param_last);
                let origin = cir.center;
                let dir1 = (point_first - origin).normalize_or_zero();
                let dir2 = (point_last - origin).normalize_or_zero();
                let vdir = dir1.cross(dir2).normalize_or_zero();
                // OCCT L2443: the circle over the (Origin, Vdir, Dir1) frame.
                let a_new_circle = Circle3::new(origin, vdir, cir.radius);
                let point_last_in_chain = brep.vertex_position(&v[1]);
                let dir_last_in_chain = (point_last_in_chain - origin).normalize_or_zero();
                let mut lpar =
                    super::gap_deps::dir_angle_with_ref_3d(dir1, dir_last_in_chain, vdir);
                if lpar < 0.0 {
                    lpar += 2.0 * std::f64::consts::PI;
                }

                // OCCT L2452-2457.
                let tc = rcad_kernel::geom::Curve3::Circle(a_new_circle);
                let ne = brep_lib_make_edge(brep, &tc, &v[0], &v[1], 0.0, lpar);
                b.update_vertex_on_edge(brep, v[0].clone(), 0.0, ne.clone(), 0.0);
                b.update_vertex_on_edge(brep, v[1].clone(), lpar, ne.clone(), 0.0);
                e = ne;
            }
            self.union_pcurves(brep, the_chain, &e);
            *out_edge = e;
            return true;
        }

        // OCCT L2463-2499: the BSpline/Bezier gluing (myConcatBSplines).
        if the_chain.len() > 1 && self.my_concat_bsplines {
            let vf = sae.first_vertex(brep, &the_chain[0]);
            let vl = sae.last_vertex(brep, &the_chain[the_chain.len() - 1]);
            let mut need_union = true;
            for j in 1..=the_chain.len() {
                let edge = the_chain[j - 1].clone();
                let Some((c3d, _, _)) = super::topexp::brep_tool_curve(brep, &edge) else {
                    continue;
                };
                let is_bspl_or_bez = matches!(
                    &c3d,
                    rcad_kernel::geom::Curve3::BSpline(_) | rcad_kernel::geom::Curve3::Bezier(_)
                );
                if is_bspl_or_bez {
                    continue;
                }
                need_union = false;
                break;
            }
            if need_union {
                *out_edge = super::statics_b::glue_edges_with_3d_curves(brep, the_chain, &vf, &vl);
                return true;
            }
        }
        // OCCT L2500.
        false
    }
}

/// OCCT BRep_Builder::MakeEdge(E) (BRep_Builder.cxx L875-887) — the empty
/// TEdge (default tolerance, SameParameter/SameRange flags set).
pub(crate) fn make_empty_edge(brep: &mut BRep) -> Shape {
    use rcad_kernel::topods::tshape_flags;
    let tshape = Arc::new(TShape::Edge(rcad_kernel::topods::TEdgeData {
        my_shapes: Vec::new(),
        flags: tshape_flags::DEFAULT,
        curve: None,
        first: Shape::null(),
        last: Shape::null(),
        range: [0.0, 0.0],
        degenerated: false,
        pcurves: indexmap::IndexMap::new(),
        representations: Vec::new(),
        vertex_params: std::collections::HashMap::new(),
        tolerance: CONFUSION,
        same_parameter: true,
        same_range: true,
    }));
    let index = brep.tshapes.len();
    brep.tshapes.push(tshape);
    Shape {
        data: brep.tshapes[index].clone(),
        index,
        orientation: Orientation::Forward,
        location: 0,
    }
}

/// The full-circle edge constructor (the B.MakeEdge(E, Cir, Confusion())
/// form of cxx L2387/L2403).
pub(crate) fn make_full_circle_edge(brep: &mut BRep, cir: &Circle3) -> Shape {
    let e = make_empty_edge(brep);
    // SAFETY: fresh edge, no aliases (the in-place curve set mirroring the
    // BRep_Builder::MakeEdge(E, C, Tol) payload write).
    let ptr = Arc::as_ptr(&brep.tshapes[e.index]) as *mut TShape;
    let ts = unsafe { &mut *ptr };
    if let TShape::Edge(ed) = ts {
        ed.curve = Some(rcad_kernel::geom::Curve3::Circle(*cir));
        ed.tolerance = ed.tolerance.max(CONFUSION);
    }
    e
}

/// OCCT static IsMergingPossible (cxx L2508-2638).
#[allow(clippy::too_many_arguments)]
pub fn is_merging_possible(
    brep: &mut BRep,
    edge1: &Shape,
    edge2: &Shape,
    the_ang_tol: f64,
    the_lin_tol: f64,
    avoid_edge_vrt: &super::MapOfShape,
    the_line_direction_ok: bool,
    the_first_point: DVec3,
    the_direction_vec: DVec3,
    the_vfmap: &IndexedDataMapOfShapeListOfShape,
) -> bool {
    // OCCT L2521-2522.
    let is_deg_e1 = is_edge_degenerated(brep, edge1);
    let is_deg_e2 = is_edge_degenerated(brep, edge2);

    if is_deg_e1 && is_deg_e2 {
        // OCCT L2526-2552: the 2D connection-point check.
        let mut common_face = Shape::null();
        let mut min_sq_dist = 0.0f64;
        let mut or_of_e1 = Orientation::Forward;
        let mut or_of_e2 = Orientation::Forward;
        let mut ind_on_e1 = 0i32;
        let mut ind_on_e2 = 0i32;
        let mut points_on_edge1 = [DVec2::ZERO; 2];
        let mut points_on_edge2 = [DVec2::ZERO; 2];
        if !find_closest_points(
            brep,
            edge1,
            edge2,
            the_vfmap,
            &mut common_face,
            &mut min_sq_dist,
            &mut or_of_e1,
            &mut or_of_e2,
            &mut ind_on_e1,
            &mut ind_on_e2,
            &mut points_on_edge1,
            &mut points_on_edge2,
        ) {
            return false;
        }
        // OCCT L2547-2552.
        return min_sq_dist <= square_p_confusion();
    } else if is_deg_e1 || is_deg_e2 {
        // OCCT L2554-2557.
        return false;
    }

    // OCCT L2559-2563.
    let cv = last_vertex(brep, edge1, true);
    if cv.is_null() || avoid_edge_vrt.contains_key(&super::shape_key(&cv)) {
        return false;
    }

    // OCCT L2565-2577.
    let ade1 = BRepAdaptorCurve::new(brep, edge1);
    let ade2 = BRepAdaptorCurve::new(brep, edge2);
    let t1 = ade1.get_type();
    let t2 = ade2.get_type();
    if t1 == GeomAbsCurveType::Circle && t2 == GeomAbsCurveType::Circle {
        if let (Some(c1), Some(c2)) = (ade1.circle(), ade2.circle()) {
            if c1.center.distance(c2.center) > CONFUSION {
                return false;
            }
        }
    }

    // OCCT L2579-2586.
    let mut a_dir1 = DVec3::ONE;
    let mut a_dir2 = DVec3::ONE;
    let lin_pair = is_linear(&ade1, &mut a_dir1) && is_linear(&ade2, &mut a_dir2);
    if !lin_pair
        && ((t1 != GeomAbsCurveType::BezierCurve && t1 != GeomAbsCurveType::BSplineCurve)
            || (t2 != GeomAbsCurveType::BezierCurve && t2 != GeomAbsCurveType::BSplineCurve))
        && t1 != t2
    {
        return false;
    }

    // OCCT L2588-2608: the tangent-direction agreement at the junction.
    let (_p1, mut diff1) = if edge1.orientation == Orientation::Forward {
        ade1.d1(ade1.last_parameter())
    } else {
        let (p, d) = ade1.d1(ade1.first_parameter());
        (p, -d)
    };
    let (_p2, mut diff2) = if edge2.orientation == Orientation::Forward {
        ade2.d1(ade2.first_parameter())
    } else {
        let (p, d) = ade2.d1(ade2.last_parameter());
        (p, -d)
    };
    if super::gap_deps::dir_angle_3d(diff1.normalize_or_zero(), diff2.normalize_or_zero())
        > the_ang_tol
    {
        return false;
    }

    // OCCT L2615-2635: the accumulated deflection/angle checks.
    if the_line_direction_ok && t2 == GeomAbsCurveType::Line {
        let a_last = if edge2.orientation == Orientation::Forward {
            ade2.last_parameter()
        } else {
            ade2.first_parameter()
        };
        let a_cur_v = ade2.value(a_last) - the_first_point;
        let a_dd = the_direction_vec.cross(a_cur_v).length_squared();
        if a_dd > the_lin_tol * the_lin_tol {
            return false;
        }
        if super::gap_deps::dir_angle_3d(
            the_direction_vec.normalize_or_zero(),
            a_cur_v.normalize_or_zero(),
        ) > the_ang_tol
            || super::gap_deps::dir_angle_3d(diff2.normalize_or_zero(), a_cur_v.normalize_or_zero())
                > the_ang_tol
        {
            return false;
        }
    }
    let _ = &mut diff1;

    // OCCT L2637.
    true
}

/// OCCT static GetLineEdgePoints (cxx L2642-2675).
pub fn get_line_edge_points(
    brep: &BRep,
    the_inp_edge: &Shape,
    the_first_point: &mut DVec3,
    the_direction_vec: &mut DVec3,
) -> bool {
    // OCCT L2646-2651.
    let Some((a_cur, f, l)) = super::topexp::brep_tool_curve(brep, the_inp_edge) else {
        return false;
    };
    // OCCT L2653-2657: the trimmed basis.
    let a_cur = match &a_cur {
        rcad_kernel::geom::Curve3::Trimmed(t) => (*t.curve).clone(),
        other => other.clone(),
    };
    // OCCT L2659-2662: Geom_Line only.
    let rcad_kernel::geom::Curve3::Line(_) = &a_cur else {
        return false;
    };
    // OCCT L2664-2669.
    let (f, l) = if the_inp_edge.orientation == Orientation::Reversed {
        (l, f)
    } else {
        (f, l)
    };
    *the_first_point = a_cur.point_at(f);
    let a_lp = a_cur.point_at(l);
    *the_direction_vec = (a_lp - *the_first_point).normalize_or_zero();
    true
}
