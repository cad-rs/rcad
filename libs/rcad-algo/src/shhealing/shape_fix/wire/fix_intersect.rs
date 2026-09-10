//! OCCT `ShapeFix_Wire` — the `FixIntersectingEdges(num)` and
//! `FixIntersectingEdges(num1, num2)` methods (`ShapeFix_Wire.cxx`
//! L2966-3489).  The impl submodule of [`super::ShapeFixWire`] (the OCCT
//! continued-file convention); split out of `fix_adv.rs` for the
//! 2000-line file budget.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::Curve3;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, TShape};

use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};
use crate::shhealing::shape_fix::split_tool::ShapeFixSplitTool;
use crate::shhealing::shape_fix::wire::wire_statics::{
    compute_local_deviation, param_on_first, param_on_second, pcurve_range_on_face,
};
use crate::shhealing::shape_fix::wire::{
    brep_tool_pnt, brep_tool_same_parameter, brep_tool_tolerance, REAL_LAST,
    ShapeFixWire,
};

/// OCCT `gp::Resolution()` (gp.hxx L59-60) — the minimal positive double.
pub(crate) const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

impl ShapeFixWire {
    // OCCT ShapeFix_Wire.cxx L2966-3265 — FixIntersectingEdges(num).
    /// OCCT ShapeFix_Wire::FixIntersectingEdges(num) (cxx L2966-3265): tests
    /// if two consequent edges are intersecting and fixes the intersection
    /// by increasing the tolerance of the vertex between the edges, shifting
    /// this vertex to the intersection point, or cutting the edges to the
    /// intersection point.
    pub fn fix_intersecting_edges(&mut self, brep: &mut BRep, num: i32) -> bool {
        self.fix_intersecting_edges_impl(brep, num)
    }

    /// The FixIntersectingEdges(num) body (the shared implementation; kept
    /// separate so the num1/num2 form is its own documented function).
    fn fix_intersecting_edges_impl(&mut self, brep: &mut BRep, num: i32) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() || self.nb_edges() < 2 {
            return false;
        }

        // analysis
        let mut points2d: Vec<DVec2> = Vec::new();
        let mut points3d: Vec<DVec3> = Vec::new();
        let mut errors: Vec<f64> = Vec::new();
        self.my_analyzer.check_intersecting_edges_full(
            brep,
            num,
            &mut points2d,
            &mut points3d,
            &mut errors,
        );
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail) {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
        }
        if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
            return false;
        }

        // rln 03/02/98: CSR#BUC50004 entity 56 (to avoid later inserting a
        // lacking edge)

        // action: increase the tolerance of the vertex

        let nb_now = self.my_analyzer.wire_data().unwrap().nb_edges();
        let n2 = if num > 0 { num } else { nb_now };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_now };
        let mut e1 = self.my_analyzer.wire_data().unwrap().edge(n1);
        let mut e2 = self.my_analyzer.wire_data().unwrap().edge(n2);
        if self.base.my_context.is_some() {
            e1 = self.context_apply(brep, &e1);
            e2 = self.context_apply(brep, &e2);
        }

        let is_forward1 = e1.orientation == Orientation::Forward;
        let is_forward2 = e2.orientation == Orientation::Forward;
        let face = self.face();
        let (_, mut a1, mut b1) = pcurve_range_on_face(brep, &e1, &face);
        let (_, mut a2, mut b2) = pcurve_range_on_face(brep, &e2, &face);

        let sae = ShapeAnalysisEdge::new();
        let mut vp = sae.first_vertex(brep, &e1);
        let mut v1 = sae.last_vertex(brep, &e1);
        let mut v2 = sae.first_vertex(brep, &e2);
        let mut vn = sae.last_vertex(brep, &e2);

        let mut tol = brep_tool_tolerance(&v1);
        let mut pnt = brep_tool_pnt(&v1);

        let mut prev_range1 = REAL_LAST;
        let mut prev_range2 = REAL_LAST;
        let mut cut_edge1 = false;
        let mut cut_edge2 = false;
        let mut is_cut_line = false;
        let mut is_changed_edge = false;

        let mut b = BRepBuilder::new();

        let nb = points3d.len();
        for i in 1..=nb {
            let ip = &points2d[i - 1];
            let param1 = if num == 1 {
                param_on_second(ip)
            } else {
                param_on_first(ip)
            };
            let param2 = if num == 1 {
                param_on_first(ip)
            } else {
                param_on_second(ip)
            };

            let new_range1 = ((if is_forward1 { a1 } else { b1 }) - param1).abs();
            let new_range2 = ((if is_forward2 { b2 } else { a2 }) - param2).abs();
            if new_range1 > prev_range1 && new_range2 > prev_range2 {
                continue;
            }

            let pint = points3d[i - 1];
            let rad = errors[i - 1];
            let mut newtol = 1.0001 * (pnt.distance(pint) + rad);

            //: r8 abv 12 Apr 99: try increasing the tolerance of the edge

            let mut loc_may_edit = self.my_topo_mode;
            // Always try to modify the tolerance firstly as a better solution
            if /* ! myTopoMode &&*/ newtol > tol {
                let te1 = rad
                    + compute_local_deviation(
                        brep,
                        &e1,
                        pint,
                        pnt,
                        param1,
                        if is_forward1 { b1 } else { a1 },
                        &face,
                    );
                let te2 = rad
                    + compute_local_deviation(
                        brep,
                        &e2,
                        pint,
                        pnt,
                        if is_forward2 { a2 } else { b2 },
                        param2,
                        &face,
                    );
                let maxte = te1.max(te2);
                if maxte < self.base.my_max_tol && maxte < newtol {
                    if brep_tool_tolerance(&e1) < te1 || brep_tool_tolerance(&e2) < te2 {
                        // Make a copy of the edges.
                        if self.base.my_context.is_some() {
                            is_changed_edge = true; // To avoid double copying of vertexes.

                            // Intersection point of two base edges.
                            let mut vv1 = {
                                let ctx = self.base.my_context.as_mut().unwrap();
                                ctx.copy_vertex(brep, &v1, -1.0)
                            };

                            let mut vvp = vp.clone();
                            let mut vvn = vn.clone();
                            if vp.is_same(&vn) {
                                // Should modify only one vertex.
                                vvp = {
                                    let ctx = self.base.my_context.as_mut().unwrap();
                                    ctx.copy_vertex(brep, &vp, -1.0)
                                };
                                vvn = vvp.clone();
                            } else {
                                vvp = {
                                    let ctx = self.base.my_context.as_mut().unwrap();
                                    ctx.copy_vertex(brep, &vp, -1.0)
                                };
                                vvn = {
                                    let ctx = self.base.my_context.as_mut().unwrap();
                                    ctx.copy_vertex(brep, &vn, -1.0)
                                };
                            }

                            let ee1 = sae_copy_replace(brep, &e1, &vvp, &vv1);
                            let ee2 = sae_copy_replace(brep, &e2, &vv1, &vvn);

                            if let Some(ctx) = self.base.my_context.as_mut() {
                                ctx.replace(brep, &e1, &ee1);
                                ctx.replace(brep, &e2, &ee2);
                            }

                            self.update_wire(brep);
                            e1 = self.my_analyzer.wire_data().unwrap().edge(n1);
                            e2 = self.my_analyzer.wire_data().unwrap().edge(n2);
                            vp = sae.first_vertex(brep, &e1);
                            v1 = sae.last_vertex(brep, &e1);
                            v2 = sae.first_vertex(brep, &e2);
                            vn = sae.last_vertex(brep, &e2);
                            vv1 = v1.clone();
                        }

                        b.update_edge_tolerance(brep, e1.clone(), 1.000001 * te1);
                        b.update_vertex_tolerance(brep, sae.first_vertex(brep, &e1), 1.000001 * te1);
                        b.update_vertex_tolerance(brep, sae.last_vertex(brep, &e1), 1.000001 * te1);
                        b.update_edge_tolerance(brep, e2.clone(), 1.000001 * te2);
                        b.update_vertex_tolerance(brep, sae.first_vertex(brep, &e2), 1.000001 * te2);
                        b.update_vertex_tolerance(brep, sae.last_vertex(brep, &e2), 1.000001 * te2);

                        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done6);
                        loc_may_edit = false;
                    }
                    newtol = 1.000001 * maxte;
                }
            }

            if loc_may_edit || newtol <= self.base.my_max_tol {
                prev_range1 = new_range1;
                prev_range2 = new_range2;
                if loc_may_edit {
                    newtol = 1.0001 * (pnt.distance(pint) + rad);
                    //: j6 abv 7 Dec 98: ProSTEP TR10 r0601_id.stp #57676 &
                    // #58586: do not cut edges because of the influence on
                    // adjacent faces
                    let mut a_tool = ShapeFixSplitTool;

                    if !a_tool.cut_edge(
                        brep,
                        &e1,
                        if is_forward1 { a1 } else { b1 },
                        param1,
                        &face,
                        &mut is_cut_line,
                    ) {
                        if v1.is_same(&vp) {
                            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done3);
                        } else {
                            loc_may_edit = false;
                        }
                    } else {
                        cut_edge1 = true; //: h4
                    }

                    if !a_tool.cut_edge(
                        brep,
                        &e2,
                        if is_forward2 { b2 } else { a2 },
                        param2,
                        &face,
                        &mut is_cut_line,
                    ) {
                        if v2.is_same(&vn) {
                            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done4);
                        } else {
                            loc_may_edit = false;
                        }
                    } else {
                        cut_edge2 = true; //: h4
                    }
                }

                let same_parameter_ok =
                    brep_tool_same_parameter(&e1) && brep_tool_same_parameter(&e2);
                if loc_may_edit
                    && new_range1 <= prev_range1
                    && new_range2 <= prev_range2
                    && same_parameter_ok
                {
                    self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);
                    pnt = pint;
                    if tol <= rad {
                        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
                        tol = 1.001 * rad;
                    }
                } else if is_cut_line {
                    self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);
                    pnt = pint;
                    if tol <= rad {
                        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
                        tol = 1.001 * rad;
                    }
                } else {
                    // else increase tolerance
                    if tol < newtol {
                        // rln 07.04.99 CCI60005-brep.igs
                        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
                        tol = newtol;
                    }
                }
            } else {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
            }
        }

        if !self.last_fix_status(ShapeExtendStatus::Done) {
            return false;
        }

        if is_changed_edge {
            b.update_vertex_point(brep, v1.clone(), pnt, tol);
            b.update_vertex_point(brep, v2.clone(), pnt, tol);
        } else if self.base.my_context.is_some() {
            if v1.is_same(&v2) {
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.copy_vertex_at(brep, &v1, pnt, tol);
                }
            } else {
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.copy_vertex_at(brep, &v1, pnt, tol);
                    ctx.copy_vertex_at(brep, &v2, pnt, tol);
                }
            }
        } else {
            b.update_vertex_point(brep, v1.clone(), pnt, tol);
            b.update_vertex_point(brep, v2.clone(), pnt, tol);
        }

        //: h4: make edges SP (after all cuts: t4mug.stp #3730+#6460)
        if cut_edge1 {
            if self.base.my_context.is_some() {
                e1 = self.context_apply(brep, &e1);
            }
            self.my_fix_edge.fix_same_parameter(brep, &e1, 0.0);
        }
        if cut_edge2 && !is_cut_line {
            if self.base.my_context.is_some() {
                e2 = self.context_apply(brep, &e2);
            }
            self.my_fix_edge.fix_same_parameter(brep, &e2, 0.0);
        }
        if cut_edge1 || cut_edge2 {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done7);
        }
        if !self.base.my_shape.is_null() {
            // Edges were intersecting, corrected
            self.base.send_warning_own(&MessageMsg::from_key("FixAdvWire.FixIntersection.MSG10"));
        }
        let _ = (&mut a1, &mut b1, &mut a2, &mut b2);
        true
    }

    // OCCT ShapeFix_Wire.cxx L3271-3489 — FixIntersectingEdges(num1, num2).
    /// OCCT ShapeFix_Wire::FixIntersectingEdges(num1, num2) (cxx
    /// L3271-3489): tests if the two edges num1 and num2 are intersecting
    /// and fixes the intersection by increasing the tolerance of the vertex
    /// nearest to the intersection point.
    pub fn fix_intersecting_edges_pair(&mut self, brep: &mut BRep, num1: i32, num2: i32) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }
        let mut points2d: Vec<DVec2> = Vec::new();
        let mut points3d: Vec<DVec3> = Vec::new();
        let mut errors: Vec<f64> = Vec::new();
        self.my_analyzer.check_intersecting_edges_full_pair(
            brep,
            num1,
            num2,
            &mut points2d,
            &mut points3d,
            &mut errors,
        );
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail) {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
        }
        if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
            return false;
        }
        let mut vertex_points = [DVec3::ZERO; 4];
        let mut vertex_tolers = [0.0f64; 4];
        let mut new_tolers = [0.0f64; 4];
        let mut vertices: [Option<Shape>; 4] = [None, None, None, None];

        let sbwd_nb = self.my_analyzer.wire_data().unwrap().nb_edges();
        let n2 = if num1 > 0 { num1 } else { sbwd_nb };
        let n1 = if num2 > 1 { num2 } else { sbwd_nb };
        if n1 == n2 {
            return false;
        }

        let edge1 = self.my_analyzer.wire_data().unwrap().edge(n1);
        let edge2 = self.my_analyzer.wire_data().unwrap().edge(n2);

        let sae = ShapeAnalysisEdge::new();
        vertices[0] = Some(sae.first_vertex(brep, &edge1));
        vertices[1] = Some(sae.last_vertex(brep, &edge1));
        vertices[2] = Some(sae.first_vertex(brep, &edge2));
        vertices[3] = Some(sae.last_vertex(brep, &edge2));

        for i in 0..4 {
            let vtx = vertices[i].as_ref().unwrap();
            vertex_points[i] = brep_tool_pnt(vtx);
            vertex_tolers[i] = brep_tool_tolerance(vtx);
        }

        let mut a_new_tol_edge1 = 0.0f64;
        let mut a_new_tol_edge2 = 0.0f64;
        let nb = points3d.len();
        for i in 1..=nb {
            let pint = points3d[i - 1];

            // searching for the nearest vertices to the intersection point
            let mut a_vtx1_param = 0.0f64;
            let mut a_vtx2_param = 0.0f64;
            let mut a_min_dist = REAL_LAST;
            let mut a_nearest_vertex = DVec3::ZERO;
            let mut a_necessary_vtx_tole = 0.0f64;
            for a_vc1 in 1..=2 {
                for a_vc2 in 3..=4 {
                    let a_vtx_ip_dist = pint.distance(vertex_points[a_vc1 - 1]);
                    let a_vtx_vtx_dist = vertex_points[a_vc1 - 1].distance(vertex_points[a_vc2 - 1]);
                    if a_min_dist > a_vtx_ip_dist && a_vtx_ip_dist > a_vtx_vtx_dist {
                        a_necessary_vtx_tole = a_vtx_vtx_dist;
                        a_nearest_vertex = vertex_points[a_vc1 - 1];
                        a_min_dist = a_vtx_ip_dist;
                        a_vtx1_param = brep_tool_parameter(
                            brep,
                            vertices[(a_vc1 - 1) as usize].as_ref().unwrap(),
                            &edge1,
                        );
                        a_vtx2_param = brep_tool_parameter(
                            brep,
                            vertices[(a_vc2 - 1) as usize].as_ref().unwrap(),
                            &edge2,
                        );
                    }
                }
            }

            // calculation of the necessary tolerances of the edges
            let ip = &points2d[i - 1];
            let param1 = param_on_first(ip);
            let param2 = param_on_second(ip);
            let (mut a_curve1, loc1, f, l) = curve_with_loc(brep, &edge1);
            let (mut a_curve2, loc2, _f2, _l2) = curve_with_loc(brep, &edge2);
            let _ = (f, l);

            // if aMinDist is lower than the resolution then the intersection
            // point lies inside the vertex
            if a_min_dist < GP_RESOLUTION {
                continue;
            }

            let mut a_max_edge_tol1 = 0.0f64;
            let mut a_max_edge_tol2 = 0.0f64;
            if a_min_dist < REAL_LAST && a_curve1.is_some() && a_curve2.is_some() {
                // OCCT L3374: gp_Lin aLig(aNearestVertex, gp_Vec(aNearestVertex, pint)).
                let a_lig_dir = (pint - a_nearest_vertex).normalize_or_zero();
                let mut du1 = 0.05 * (param1 - a_vtx1_param);
                let mut du2 = 0.05 * (param2 - a_vtx2_param);
                let tole1 = brep_tool_tolerance(&edge1);
                let tole2 = brep_tool_tolerance(&edge2);
                for a_points_c in 2..19 {
                    let u = a_vtx1_param + a_points_c as f64 * du1;
                    let p1 = rcad_kernel::geom::CurveEval::point_at(a_curve1.as_ref().unwrap(), u);
                    let p1 = apply_loc(brep, loc1, p1);
                    let d1 = point_line_distance(p1, a_nearest_vertex, a_lig_dir) * 2.0000001;
                    if d1 > tole1 && d1 > a_max_edge_tol1 {
                        a_max_edge_tol1 = d1;
                    }

                    let u = a_vtx2_param + a_points_c as f64 * du2;
                    let p2 = rcad_kernel::geom::CurveEval::point_at(a_curve2.as_ref().unwrap(), u);
                    let p2 = apply_loc(brep, loc2, p2);
                    let d2 = point_line_distance(p2, a_nearest_vertex, a_lig_dir) * 2.0000001;
                    if d2 > tole2 && d2 > a_max_edge_tol2 {
                        a_max_edge_tol2 = d2;
                    }
                }
                if a_max_edge_tol1 == 0.0 && a_max_edge_tol2 == 0.0 {
                    continue;
                }
                // if the vertices are farther than the tolerances then we do
                // not need to increase the edge tolerance
                if a_necessary_vtx_tole > a_max_edge_tol1.max(tole1)
                    || a_necessary_vtx_tole > a_max_edge_tol2.max(tole2)
                {
                    a_max_edge_tol1 = 0.0;
                    a_max_edge_tol2 = 0.0;
                }
                let _ = &mut a_curve1;
                let _ = &mut a_curve2;
                du1 += 0.0;
                du2 += 0.0;
            }

            let rad = errors[i - 1];
            let mut fin_tol = REAL_LAST;
            let mut rank = 1usize;
            for j in 1..=4 {
                let newtol = 1.0001 * (pint.distance(vertex_points[j - 1]) + rad);
                if newtol < fin_tol {
                    rank = j;
                    fin_tol = newtol;
                }
            }
            if fin_tol <= self.base.my_max_tol {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
                if new_tolers[rank - 1] < fin_tol {
                    if a_max_edge_tol1.max(a_max_edge_tol2) < fin_tol
                        && (a_max_edge_tol1 > 0.0 || a_max_edge_tol2 > 0.0)
                    {
                        a_new_tol_edge1 = a_new_tol_edge1.max(a_max_edge_tol1);
                        a_new_tol_edge2 = a_new_tol_edge2.max(a_max_edge_tol2);
                    } else {
                        new_tolers[rank - 1] = fin_tol;
                    }
                }
            } else {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
            }
        }

        let mut b = BRepBuilder::new();
        // update of the edge tolerances
        if a_new_tol_edge1 > 0.0 {
            for i in 1..=2 {
                if a_new_tol_edge1 > vertex_tolers[i - 1].max(new_tolers[i - 1]) {
                    new_tolers[i - 1] = a_new_tol_edge1;
                }
            }
            b.update_edge_tolerance(brep, edge1.clone(), a_new_tol_edge1);
        }
        if a_new_tol_edge2 > 0.0 {
            for i in 3..=4 {
                if a_new_tol_edge2 > vertex_tolers[i - 1].max(new_tolers[i - 1]) {
                    new_tolers[i - 1] = a_new_tol_edge2;
                }
            }
            b.update_edge_tolerance(brep, edge2.clone(), a_new_tol_edge2);
        }

        // update of the vertex tolerances
        for i in 1..=4 {
            if new_tolers[i - 1] > 0.0 {
                let v = vertices[i - 1].clone().unwrap();
                b.update_vertex_tolerance(brep, v, new_tolers[i - 1]);
            }
        }

        if !self.base.my_shape.is_null() {
            // Edges were intersecting, corrected
            self.base.send_warning_own(&MessageMsg::from_key("FixAdvWire.FixIntersection.MSG10"));
        }
        true
    }
}

// ---------------------------------------------------------------------------
// File-local helpers.
// ---------------------------------------------------------------------------

/// OCCT ShapeAnalysis_Edge::CopyReplaceVertices usage site — the
/// ShapeBuild_Edge tool call (kept for the OCCT statement shape).
fn sae_copy_replace(brep: &mut BRep, edge: &Shape, v1: &Shape, v2: &Shape) -> Shape {
    ShapeBuildEdge.copy_replace_vertices(brep, edge, v1, v2)
}

/// OCCT BRep_Tool::Parameter(vertex, edge).
fn brep_tool_parameter(brep: &BRep, v: &Shape, e: &Shape) -> f64 {
    match brep.tshapes[e.index].as_ref() {
        TShape::Edge(ed) => *ed.vertex_params.get(&v.ptr_id()).unwrap_or(&0.0),
        _ => 0.0,
    }
}

/// OCCT `BRep_Tool::Curve(edge, L, f, l)` — the curve and its location.
fn curve_with_loc(brep: &BRep, e: &Shape) -> (Option<Curve3>, u32, f64, f64) {
    match brep.tshapes[e.index].as_ref() {
        TShape::Edge(ed) => (ed.curve.clone(), e.location, ed.range[0], ed.range[1]),
        _ => (None, 0, 0.0, 0.0),
    }
}

/// OCCT `P.Transform(L.Transformation())`.
fn apply_loc(brep: &BRep, loc: u32, p: DVec3) -> DVec3 {
    if loc != 0 {
        brep.get_location(loc).transform_point3(p)
    } else {
        p
    }
}

/// OCCT `gp_Lin::Distance(P)` — the point-to-line distance.
fn point_line_distance(p: DVec3, origin: DVec3, dir: DVec3) -> f64 {
    (p - origin).cross(dir).length()
}

