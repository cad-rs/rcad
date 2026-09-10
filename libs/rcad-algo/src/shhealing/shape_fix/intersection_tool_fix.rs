//! OCCT `ShapeFix_IntersectionTool` — the fixing drivers (the impl
//! submodule of [`super::intersection_tool::ShapeFixIntersectionTool`], the
//! OCCT continued-file convention): `FixSelfIntersectWire` (cxx L1029-1831)
//! and `FixIntersectingWires` (cxx L1835-2530) — split out of
//! `intersection_tool.rs` for the 2000-line file budget.

use std::collections::HashMap;

use rcad_kernel::geom::{Curve2d, Curve2dEval};
use rcad_kernel::math::bnd::BndBox2d;
use rcad_kernel::precision::CONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType};

use crate::geomalgo::int_res2d::{Domain, IntersectionPoint, Transition, Position};
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes, topexp_explorer};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::intersection_tool::{
    create_boxes2d, get_point_on_edge, select_int_pnt, shape_is_equal, Geom2dIntGInterGap,
    ShapeBoxes2d,
};
use crate::shhealing::shape_fix::intersection_tool::ShapeFixIntersectionTool;
use crate::shhealing::shape_fix::split_tool::{brep_loc_transform, brep_tool_surface};
use crate::shhealing::shape_fix::wire::{
    brep_tool_degenerated, brep_tool_pnt, brep_tool_same_parameter, brep_tool_tolerance,
    shape_oriented,
};

impl ShapeFixIntersectionTool {
    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L1029-1831 — FixSelfIntersectWire.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::FixSelfIntersectWire (cxx L1029-1831).
    pub fn fix_self_intersect_wire(
        &mut self,
        brep: &mut BRep,
        sewd: &mut WireData,
        face: &Shape,
        nb_split: &mut i32,
        nb_cut: &mut i32,
        nb_removed: &mut i32,
    ) -> bool {
        // OCCT L1035-1038.
        if self.my_context.is_none() || face.is_null() {
            return false;
        }

        // double area2d = ShapeAnalysis::TotCross2D(sewd,face);
        // if(area2d<Precision::PConfusion()*Precision::PConfusion()) return false; //gka
        // 06.09.04 BUG 6555

        // OCCT L1044-1052.
        let ctx = self.my_context.as_mut().unwrap();
        let sf = ctx.apply(brep, face, ShapeType::Shape);
        let mut max_tol_vert = 0.0f64;
        for exp in topexp_explorer(brep, &sf, ShapeType::Vertex) {
            let tol_v = brep_tool_tolerance(&exp);
            max_tol_vert = max_tol_vert.max(tol_v);
        }
        max_tol_vert = max_tol_vert.min(self.my_max_tol);
        let sae = ShapeAnalysisEdge::new();

        // OCCT L1054: step 1 : intersection of adjacent edges

        // OCCT L1056-1059: step 2 : intersection of non-adjacent edges
        let mut boxes: ShapeBoxes2d = HashMap::new();
        let _ = create_boxes2d(brep, sewd, face, &mut boxes);
        let sas = ShapeAnalysisSurface::new(brep_tool_surface(brep, face).unwrap());

        // OCCT L1061-1064.
        *nb_split = 0;
        *nb_cut = 0;
        let mut nb_replaced = 0i32;
        let mut is_done;

        // OCCT L1065.
        let mut num1 = 1i32;
        while num1 < sewd.nb_edges() && *nb_split < 30 {
            // for each edge from first wire
            let mut num2 = num1 + 2;
            while num2 <= sewd.nb_edges() && *nb_split < 30 {
                // for each edge from second wire
                // OCCT L1071-1074.
                if num1 == 1 && num2 == sewd.nb_edges() {
                    num2 += 1;
                    continue;
                }
                // OCCT L1075-1080.
                let mut edge1 = sewd.edge(num1);
                let mut edge2 = sewd.edge(num2);
                if edge1.is_same(&edge2) {
                    num2 += 1;
                    continue;
                }
                // OCCT L1081-1084.
                if brep_tool_degenerated(&edge1) || brep_tool_degenerated(&edge2) {
                    num2 += 1;
                    continue;
                }
                // OCCT L1085-1088.
                if !boxes.contains_key(&(edge1.ptr_id(), edge1.location))
                    || !boxes.contains_key(&(edge2.ptr_id(), edge2.location))
                {
                    num2 += 1;
                    continue;
                }
                // OCCT L1089-1091.
                let b1 = boxes.get(&(edge1.ptr_id(), edge1.location)).unwrap().clone();
                let b2 = boxes.get(&(edge2.ptr_id(), edge2.location)).unwrap().clone();
                if !b1.is_out_box(&b2) {
                    // intersection is possible...
                    // OCCT L1094-1103.
                    let mut a1 = 0.0f64;
                    let mut b1p = 0.0f64;
                    let mut a2 = 0.0f64;
                    let mut b2p = 0.0f64;
                    let mut crv1: Option<Curve2d> = None;
                    let mut crv2: Option<Curve2d> = None;
                    if !sae.pcurve_face(brep, &edge1, face, &mut crv1, &mut a1, &mut b1p, false) {
                        return false;
                    }
                    if !sae.pcurve_face(brep, &edge2, face, &mut crv2, &mut a2, &mut b2p, false) {
                        return false;
                    }
                    let (b1v, b2v) = (b1p, b2p);
                    // OCCT L1104-1109.
                    let tolint = 1.0e-10f64;
                    let crv1 = crv1.unwrap();
                    let crv2 = crv2.unwrap();
                    let d1 = Domain::bounded(
                        Curve2dEval::point_at(&crv1, a1),
                        a1,
                        tolint,
                        Curve2dEval::point_at(&crv1, b1v),
                        b1v,
                        tolint,
                    );
                    let d2 = Domain::bounded(
                        Curve2dEval::point_at(&crv2, a2),
                        a2,
                        tolint,
                        Curve2dEval::point_at(&crv2, b2v),
                        b2v,
                        tolint,
                    );
                    let mut inter = Geom2dIntGInterGap::new();
                    inter.perform(&crv1, &d1, &crv2, &d2, tolint, tolint);
                    // OCCT L1110-1113.
                    if !inter.is_done() {
                        num2 += 1;
                        continue;
                    }
                    // OCCT L1115: intersection is point
                    if inter.nb_points() > 0 && inter.nb_points() < 3 {
                        // OCCT L1117-1119.
                        let mut ip = IntersectionPoint::empty();
                        let mut tr1 = Transition::empty();
                        let mut tr2 = Transition::empty();
                        select_int_pnt(&inter, &mut ip, &mut tr1, &mut tr2);
                        // OCCT L1120-1240.
                        if tr1.position_on_curve() == Position::Middle
                            && tr2.position_on_curve() == Position::Middle
                        {
                            let param1 = ip.param_on_first();
                            let param2 = ip.param_on_second();
                            let pi1 = get_point_on_edge(brep, &edge1, &sas, &crv1, param1);
                            let pi2 = get_point_on_edge(brep, &edge2, &sas, &crv2, param2);
                            let mut b_builder = BRepBuilder::new();
                            let mut v = Shape::null();
                            let mut tol_v = 0.0f64;
                            // OCCT L1129: analysis for edge1
                            let mut modif_e1 = false;
                            let vf1 = sae.first_vertex(brep, &edge1);
                            let pvf1 = brep_tool_pnt(&vf1);
                            let vl1 = sae.last_vertex(brep, &edge1);
                            let pvl1 = brep_tool_pnt(&vl1);
                            let mut dist1 = pi1.distance(pvf1);
                            let mut dist2 = pi1.distance(pvl1);
                            let mut distmin = dist1.min(dist2);
                            if dist1 != dist2 && distmin < max_tol_vert {
                                if dist1 < dist2 {
                                    tol_v = (dist1 * 1.00001f64).max(brep_tool_tolerance(&vf1));
                                    b_builder.update_vertex_tolerance(brep, vf1.clone(), tol_v);
                                    v = vf1.clone();
                                } else {
                                    tol_v = (dist2 * 1.00001f64).max(brep_tool_tolerance(&vl1));
                                    b_builder.update_vertex_tolerance(brep, vl1.clone(), tol_v);
                                    v = vl1.clone();
                                }

                                let dista = (a1 - param1).abs();
                                let distb = (b1v - param1).abs();
                                let mut is_cut_line = false;
                                modif_e1 = self.cut_edge(
                                    brep,
                                    &mut edge1,
                                    if dista > distb { a1 } else { b1v },
                                    param1,
                                    face,
                                    &mut is_cut_line,
                                );
                                if modif_e1 {
                                    *nb_cut += 1;
                                }
                                // not needed split edge, if one of parts is too small
                                modif_e1 = modif_e1 || distmin < CONFUSION;
                            }
                            // OCCT L1164: analysis for edge2
                            let mut modif_e2 = false;
                            let vf2 = sae.first_vertex(brep, &edge2);
                            let pvf2 = brep_tool_pnt(&vf2);
                            let vl2 = sae.last_vertex(brep, &edge2);
                            let pvl2 = brep_tool_pnt(&vl2);
                            dist1 = pi2.distance(pvf2);
                            dist2 = pi2.distance(pvl2);
                            distmin = dist1.min(dist2);
                            if dist1 != dist2 && distmin < max_tol_vert {
                                if dist1 < dist2 {
                                    tol_v = (dist1 * 1.00001f64).max(brep_tool_tolerance(&vf2));
                                    b_builder.update_vertex_tolerance(brep, vf2.clone(), tol_v);
                                    v = vf2.clone();
                                } else {
                                    tol_v = (dist2 * 1.00001f64).max(brep_tool_tolerance(&vl2));
                                    b_builder.update_vertex_tolerance(brep, vl2.clone(), tol_v);
                                    v = vl2.clone();
                                }

                                let dista = (a2 - param2).abs();
                                let distb = (b2v - param2).abs();
                                let mut is_cut_line = false;
                                modif_e2 = self.cut_edge(
                                    brep,
                                    &mut edge2,
                                    if dista > distb { a2 } else { b2v },
                                    param2,
                                    face,
                                    &mut is_cut_line,
                                );
                                if modif_e2 {
                                    *nb_cut += 1;
                                }
                                // not needed split edge, if one of parts is too small
                                modif_e2 = modif_e2 || distmin < CONFUSION;
                            }
                            // OCCT L1199-1207.
                            if modif_e1 && !modif_e2 {
                                if self.split_edge1(
                                    brep, sewd, face, num2, param2, &v, tol_v, &mut boxes,
                                ) {
                                    *nb_split += 1;
                                    num2 -= 1;
                                    num2 += 1;
                                    continue;
                                }
                            }
                            // OCCT L1208-1216.
                            if !modif_e1 && modif_e2 {
                                if self.split_edge1(
                                    brep, sewd, face, num1, param1, &v, tol_v, &mut boxes,
                                ) {
                                    *nb_split += 1;
                                    num1 -= 1;
                                    break;
                                }
                            }
                            // OCCT L1217-1239.
                            if !modif_e1 && !modif_e2 {
                                let p0 = glam::DVec3::new(
                                    (pi1.x + pi2.x) / 2.0,
                                    (pi1.y + pi2.y) / 2.0,
                                    (pi1.z + pi2.z) / 2.0,
                                );
                                tol_v = ((pi1.distance(pi2) / 2.0) * 1.00001f64).max(CONFUSION);
                                let made = {
                                    let mut b = BRepBuilder::new();
                                    b.add_vertex(brep, p0, tol_v)
                                };
                                v = made;
                                max_tol_vert = max_tol_vert.max(tol_v);
                                let is_edge_split2 = self.split_edge1(
                                    brep, sewd, face, num2, param2, &v, tol_v, &mut boxes,
                                );
                                if is_edge_split2 {
                                    *nb_split += 1;
                                    num2 -= 1;
                                }
                                if self.split_edge1(
                                    brep, sewd, face, num1, param1, &v, tol_v, &mut boxes,
                                ) {
                                    *nb_split += 1;
                                    num1 -= 1;
                                    break;
                                }
                                if is_edge_split2 {
                                    num2 += 1;
                                    continue;
                                }
                            }
                        }
                        // OCCT L1241-1259.
                        if tr1.position_on_curve() == Position::Middle
                            && tr2.position_on_curve() != Position::Middle
                        {
                            // find needed vertex from edge2 and split edge1 using it
                            let param1 = ip.param_on_first();
                            if self.find_vert_and_split_edge(
                                brep,
                                param1,
                                &edge1,
                                &edge2,
                                &crv1,
                                &mut max_tol_vert,
                                &mut num1,
                                sewd,
                                face,
                                &mut boxes,
                                false,
                            ) {
                                *nb_split += 1;
                                break;
                            }
                        }
                        // OCCT L1260-1278.
                        if tr1.position_on_curve() != Position::Middle
                            && tr2.position_on_curve() == Position::Middle
                        {
                            // find needed vertex from edge1 and split edge2 using it
                            let param2 = ip.param_on_second();
                            if self.find_vert_and_split_edge(
                                brep,
                                param2,
                                &edge2,
                                &edge1,
                                &crv2,
                                &mut max_tol_vert,
                                &mut num2,
                                sewd,
                                face,
                                &mut boxes,
                                false,
                            ) {
                                *nb_split += 1;
                                num2 += 1;
                                continue;
                            }
                        }
                        // OCCT L1279-1286.
                        if tr1.position_on_curve() != Position::Middle
                            && tr2.position_on_curve() != Position::Middle
                        {
                            // union vertexes
                            if self.union_vertexes(brep, sewd, &mut edge1, &mut edge2, num2, &mut boxes, &b2) {
                                nb_replaced += 1; // gka 06.09.04
                            }
                        }
                    }
                    // OCCT L1288-1289: intersection is segment
                    if inter.nb_segments() == 1 {
                        // OCCT L1291-1314.
                        let is = inter.segment(1);
                        if is.has_first_point() && is.has_last_point() {
                            let mut is_modified1 = false;
                            let mut is_modified2 = false;
                            let mut new_v = Shape::null();
                            let mut newtol = 0.0f64;
                            let ipf = is.first_point();
                            let p11 = ipf.param_on_first();
                            let p21 = ipf.param_on_second();
                            let ipl = is.last_point();
                            let p12 = ipl.param_on_first();
                            let p22 = ipl.param_on_second();
                            let pnt11 = get_point_on_edge(brep, &edge1, &sas, &crv1, p11);
                            let pnt12 = get_point_on_edge(brep, &edge1, &sas, &crv1, p12);
                            let pnt21 = get_point_on_edge(brep, &edge2, &sas, &crv2, p21);
                            let pnt22 = get_point_on_edge(brep, &edge2, &sas, &crv2, p22);
                            // next string commented by skl 29.12.2004 for OCC7624
                            // if( Pnt11.Distance(Pnt21)>myPreci || Pnt12.Distance(Pnt22)>myPreci ) continue;
                            // OCCT L1311-1314.
                            if pnt11.distance(pnt21) > max_tol_vert || pnt12.distance(pnt22) > max_tol_vert {
                                num2 += 1;
                                continue;
                            }
                            // OCCT L1315: analysis for edge1
                            let v1 = sae.first_vertex(brep, &edge1);
                            let pv1 = brep_tool_pnt(&v1);
                            let v2 = sae.last_vertex(brep, &edge1);
                            let pv2 = brep_tool_pnt(&v2);
                            let mut dist1 = pnt11.distance(pv1);
                            let mut dist2 = pnt12.distance(pv1);
                            let mut maxdist = dist1.max(dist2);
                            let mut pdist;
                            if edge1.orientation == Orientation::Reversed {
                                pdist = (b1v - p11).abs().max((b1v - p12).abs());
                            } else {
                                pdist = (a1 - p11).abs().max((a1 - p12).abs());
                            }
                            // OCCT L1335-1341.
                            if maxdist < max_tol_vert || pdist < (b1v - a1).abs() * 0.01 {
                                // if(maxdist<maxtol || pdist<std::abs(b1-a1)*0.01) {
                                newtol = maxdist;
                                new_v = v1.clone();
                                is_modified1 = true;
                            }
                            // OCCT L1342-1362.
                            dist1 = pnt11.distance(pv2);
                            dist2 = pnt12.distance(pv2);
                            maxdist = dist1.max(dist2);
                            if edge1.orientation == Orientation::Reversed {
                                pdist = (a1 - p11).abs().max((a1 - p12).abs());
                            } else {
                                pdist = (b1v - p11).abs().max((b1v - p12).abs());
                            }
                            // if(maxdist<maxtol || pdist<std::abs(b1-a1)*0.01) {
                            if maxdist < max_tol_vert || pdist < (b1v - a1).abs() * 0.01 {
                                if (is_modified1 && maxdist < newtol) || !is_modified1 {
                                    newtol = maxdist;
                                    new_v = v2.clone();
                                    is_modified1 = true;
                                }
                            }
                            // OCCT L1363-1398.
                            if is_modified1 {
                                // cut edge1 and update tolerance NewV
                                let dista = (a1 - p11).abs() + (a1 - p12).abs();
                                let distb = (b1v - p11).abs() + (b1v - p12).abs();
                                let pend;
                                let cut;
                                if dista > distb {
                                    pend = a1;
                                } else {
                                    pend = b1v;
                                }
                                if (pend - p11).abs() > (pend - p12).abs() {
                                    cut = p12;
                                } else {
                                    cut = p11;
                                }
                                let mut is_cut_line = false;
                                if self.cut_edge(brep, &mut edge1, pend, cut, face, &mut is_cut_line) {
                                    *nb_cut += 1;
                                }
                                if newtol > brep_tool_tolerance(&new_v) {
                                    let mut b = BRepBuilder::new();
                                    b.update_vertex_tolerance(brep, new_v.clone(), newtol);
                                } else {
                                    newtol = brep_tool_tolerance(&new_v);
                                }
                            }
                            // OCCT L1399: analysis for edge2
                            let v12 = sae.first_vertex(brep, &edge2);
                            let pv12 = brep_tool_pnt(&v12);
                            let v22 = sae.last_vertex(brep, &edge2);
                            let pv22 = brep_tool_pnt(&v22);
                            dist1 = pnt21.distance(pv12);
                            dist2 = pnt22.distance(pv12);
                            maxdist = dist1.max(dist2);
                            if edge2.orientation == Orientation::Reversed {
                                pdist = (b2v - p21).abs().max((b2v - p22).abs());
                            } else {
                                pdist = (a2 - p21).abs().max((a2 - p22).abs());
                            }
                            // if(maxdist<maxtol || pdist<std::abs(b2-a2)*0.01) {
                            // OCCT L1418-1424.
                            if maxdist < max_tol_vert || pdist < (b2v - a2).abs() * 0.01 {
                                newtol = maxdist;
                                new_v = v12.clone();
                                is_modified2 = true;
                            }
                            // OCCT L1425-1445.
                            dist1 = pnt21.distance(pv22);
                            dist2 = pnt22.distance(pv22);
                            maxdist = dist1.max(dist2);
                            if edge2.orientation == Orientation::Reversed {
                                pdist = (a2 - p21).abs().max((a2 - p22).abs());
                            } else {
                                pdist = (b2v - p21).abs().max((b2v - p22).abs());
                            }
                            // if(maxdist<maxtol || pdist<std::abs(b2-a2)*0.01) {
                            if maxdist < max_tol_vert || pdist < (b2v - a2).abs() * 0.01 {
                                if (is_modified2 && maxdist < newtol) || !is_modified2 {
                                    newtol = maxdist;
                                    new_v = v22.clone();
                                    is_modified2 = true;
                                }
                            }
                            // OCCT L1446-1481.
                            if is_modified2 {
                                // cut edge1 and update tolerance NewV
                                let dista = (a2 - p21).abs() + (a2 - p22).abs();
                                let distb = (b2v - p21).abs() + (b2v - p22).abs();
                                let pend;
                                let cut;
                                if dista > distb {
                                    pend = a2;
                                } else {
                                    pend = b2v;
                                }
                                if (pend - p21).abs() > (pend - p22).abs() {
                                    cut = p22;
                                } else {
                                    cut = p21;
                                }
                                let mut is_cut_line = false;
                                if self.cut_edge(brep, &mut edge2, pend, cut, face, &mut is_cut_line) {
                                    *nb_cut += 1;
                                }
                                if newtol > brep_tool_tolerance(&new_v) {
                                    let mut b = BRepBuilder::new();
                                    b.update_vertex_tolerance(brep, new_v.clone(), newtol);
                                } else {
                                    newtol = brep_tool_tolerance(&new_v);
                                }
                            }

                            // OCCT L1483-1491.
                            if is_modified1 && !is_modified2 {
                                if self.split_edge2(
                                    brep, sewd, face, num2, p21, p22, &new_v, newtol, &mut boxes,
                                ) {
                                    *nb_split += 1;
                                    num2 -= 1;
                                    num2 += 1;
                                    continue;
                                }
                            }
                            // OCCT L1492-1500.
                            if !is_modified1 && is_modified2 {
                                if self.split_edge2(
                                    brep, sewd, face, num1, p11, p12, &new_v, newtol, &mut boxes,
                                ) {
                                    *nb_split += 1;
                                    num1 -= 1;
                                    break;
                                }
                            }
                            // OCCT L1501-1823.
                            if !is_modified1 && !is_modified2 {
                                let param1 = (p11 + p12) / 2.0;
                                let param2 = (p21 + p22) / 2.0;
                                let pnt10 = get_point_on_edge(brep, &edge1, &sas, &crv1, param1);
                                let pnt20 = get_point_on_edge(brep, &edge2, &sas, &crv2, param2);
                                let p0 = glam::DVec3::new(
                                    (pnt10.x + pnt20.x) / 2.0,
                                    (pnt10.y + pnt20.y) / 2.0,
                                    (pnt10.z + pnt20.z) / 2.0,
                                );
                                dist1 = pnt11.distance(p0).max(pnt12.distance(p0));
                                dist2 = pnt21.distance(p0).max(pnt22.distance(p0));
                                let mut tol_v = dist1.max(dist2);
                                tol_v = tol_v.max(pnt10.distance(pnt20)) * 1.00001f64;
                                let fix_segment = true;
                                if tol_v < max_tol_vert {
                                    // create new vertex and split each intersecting edge on two edges
                                    let mut b = BRepBuilder::new();
                                    let made = b.add_vertex(brep, pnt10, tol_v);
                                    new_v = made;
                                    if self.split_edge2(
                                        brep, sewd, face, num2, p21, p22, &new_v, tol_v, &mut boxes,
                                    ) {
                                        *nb_split += 1;
                                        num2 -= 1;
                                    }
                                    if self.split_edge2(
                                        brep, sewd, face, num1, p11, p12, &new_v, tol_v, &mut boxes,
                                    ) {
                                        *nb_split += 1;
                                        num1 -= 1;
                                        break;
                                    }
                                } else if fix_segment {
                                    // if( std::abs(p12-p11)>std::abs(b1-a1)/2 || std::abs(p22-p21)>std::abs(b2-a2)/2 )
                                    // {
                                    //  segment is big and we have to split each intersecting edge
                                    //  on 3 edges --> middle edge - edge based on segment
                                    //  after we can remove edges made from segment
                                    // OCCT L1538-1551.
                                    let p01 = glam::DVec3::new(
                                        (pnt11.x + pnt21.x) / 2.0,
                                        (pnt11.y + pnt21.y) / 2.0,
                                        (pnt11.z + pnt21.z) / 2.0,
                                    );
                                    let p02 = glam::DVec3::new(
                                        (pnt12.x + pnt22.x) / 2.0,
                                        (pnt12.y + pnt22.y) / 2.0,
                                        (pnt12.z + pnt22.z) / 2.0,
                                    );
                                    let mut tol_v1 = pnt11.distance(p01).max(pnt21.distance(p01));
                                    tol_v1 = tol_v1.max(CONFUSION) * 1.00001f64;
                                    let mut tol_v2 = pnt12.distance(p02).max(pnt22.distance(p02));
                                    tol_v2 = tol_v2.max(CONFUSION) * 1.00001f64;
                                    if tol_v1 > max_tol_vert || tol_v2 > max_tol_vert {
                                        num2 += 1;
                                        continue;
                                    }
                                    let mut new_v1 = Shape::null();
                                    let mut new_v2 = Shape::null();
                                    /* (the commented OCCT block L1553-1572) */
                                    let mut tmp_e = Shape::null();
                                    // OCCT L1573: split edge1
                                    let mut akey1 = 0i32;
                                    let mut akey2 = 0i32;
                                    let mut new_tolerance;
                                    // analysis fo P01
                                    // OCCT L1578-1591.
                                    new_tolerance = tol_v1.max(brep_tool_tolerance(&v1));
                                    if p01.distance(pv1) < new_tolerance {
                                        let mut b = BRepBuilder::new();
                                        let made =
                                            b.add_vertex(brep, brep_tool_pnt(&v1), new_tolerance);
                                        new_v1 = made;
                                        new_v1.orientation = v1.orientation;
                                        akey1 += 1;
                                    }
                                    new_tolerance = tol_v1.max(brep_tool_tolerance(&v2));
                                    if p01.distance(pv2) < new_tolerance {
                                        let mut b = BRepBuilder::new();
                                        let made =
                                            b.add_vertex(brep, brep_tool_pnt(&v2), new_tolerance);
                                        new_v1 = made;
                                        new_v1.orientation = v2.orientation;
                                        akey1 += 1;
                                    }
                                    // analysis fo P02
                                    // OCCT L1592-1606.
                                    new_tolerance = tol_v2.max(brep_tool_tolerance(&v1));
                                    if p02.distance(pv1) < new_tolerance {
                                        let mut b = BRepBuilder::new();
                                        let made =
                                            b.add_vertex(brep, brep_tool_pnt(&v1), new_tolerance);
                                        new_v2 = made;
                                        new_v2.orientation = v1.orientation;
                                        akey2 += 1;
                                    }
                                    new_tolerance = tol_v2.max(brep_tool_tolerance(&v2));
                                    if p02.distance(pv2) < new_tolerance {
                                        let mut b = BRepBuilder::new();
                                        let made =
                                            b.add_vertex(brep, brep_tool_pnt(&v2), new_tolerance);
                                        new_v2 = made;
                                        new_v2.orientation = v2.orientation;
                                        akey2 += 1;
                                    }
                                    // OCCT L1607-1610.
                                    if akey1 > 1 || akey2 > 1 {
                                        num2 += 1;
                                        continue;
                                    }
                                    let mut dnum1 = 0i32;
                                    let mut numseg1 = num1;
                                    // prepare vertices
                                    // OCCT L1612-1620.
                                    if akey1 == 0 {
                                        let mut b = BRepBuilder::new();
                                        new_v1 = b.add_vertex(brep, p01, tol_v1);
                                    }
                                    if akey2 == 0 {
                                        let mut b = BRepBuilder::new();
                                        new_v2 = b.add_vertex(brep, p02, tol_v2);
                                    }
                                    // split
                                    // OCCT L1621-1669.
                                    if akey1 == 0 && akey2 > 0 {
                                        if self.split_edge1(
                                            brep, sewd, face, num1, p11, &new_v1, tol_v1, &mut boxes,
                                        ) {
                                            *nb_split += 1;
                                            dnum1 = 1;
                                            numseg1 = num1 + 1;
                                        }
                                    }
                                    if akey1 > 0 && akey2 == 0 {
                                        if self.split_edge1(
                                            brep, sewd, face, num1, p12, &new_v2, tol_v2, &mut boxes,
                                        ) {
                                            *nb_split += 1;
                                            dnum1 = 1;
                                            numseg1 = num1;
                                        }
                                    }
                                    if akey1 == 0 && akey2 == 0 {
                                        if self.split_edge1(
                                            brep, sewd, face, num1, p11, &new_v1, tol_v1, &mut boxes,
                                        ) {
                                            *nb_split += 1;
                                            dnum1 = 1;
                                        }
                                        tmp_e = sewd.edge(num1);
                                        let mut a = 0.0f64;
                                        let mut b = 0.0f64;
                                        let mut c2d: Option<Curve2d> = None;
                                        sae.pcurve_face(brep, &tmp_e, face, &mut c2d, &mut a, &mut b, false);
                                        if (a - p12) * (b - p12) > 0.0 {
                                            // p12 - external for [a,b] => split next edge
                                            if self.split_edge1(
                                                brep, sewd, face, num1 + 1, p12, &new_v2, tol_v2,
                                                &mut boxes,
                                            ) {
                                                *nb_split += 1;
                                                dnum1 += 1;
                                                numseg1 = num1 + 1;
                                            }
                                        } else {
                                            if self.split_edge1(
                                                brep, sewd, face, num1, p12, &new_v2, tol_v2,
                                                &mut boxes,
                                            ) {
                                                *nb_split += 1;
                                                dnum1 += 1;
                                                numseg1 = num1 + 1;
                                            }
                                        }
                                    }
                                    // SegE = sewd->Edge(numseg1); // get edge from segment
                                    //  split edge2
                                    //  replace vertices if it is necessary
                                    // OCCT L1670-1674.
                                    let sbe = ShapeBuildEdge;
                                    akey1 = 0;
                                    akey2 = 0;

                                    // OCCT L1676-1697.
                                    if p01.distance(pv12) < tol_v1 {
                                        tol_v1 += p01.distance(pv12);
                                        let mut b = BRepBuilder::new();
                                        b.update_vertex_tolerance(brep, new_v1.clone(), tol_v1);
                                        if v12.orientation == new_v1.orientation {
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v12, &new_v1);
                                            // v12 = new_v1 (the later reads go
                                            // through the fresh vertex)
                                        } else {
                                            let new_v1r = shape_oriented(&new_v1, Orientation::Reversed);
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v12, &new_v1r);
                                        }
                                        nb_replaced += 1; // gka 06.09.04
                                        let new_e = sbe.copy_replace_vertices(brep, &edge2, &new_v1, &v22);
                                        let ctx = self.my_context.as_mut().unwrap();
                                        ctx.replace(brep, &edge2, &new_e);
                                        sewd.set_edge(&new_e, num2 + dnum1);
                                        boxes.insert(
                                            (new_e.ptr_id(), new_e.location),
                                            b2.clone(),
                                        ); // update boxes
                                        edge2 = new_e;
                                        akey1 = 1;
                                    }
                                    // OCCT L1698-1719.
                                    if p01.distance(pv22) < tol_v1 {
                                        tol_v1 += p01.distance(pv22);
                                        let mut b = BRepBuilder::new();
                                        b.update_vertex_tolerance(brep, new_v1.clone(), tol_v1);
                                        if v22.orientation == new_v1.orientation {
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v22, &new_v1);
                                        } else {
                                            let new_v1r = shape_oriented(&new_v1, Orientation::Reversed);
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v22, &new_v1r);
                                        }
                                        nb_replaced += 1; // gka 06.09.04
                                        let new_e = sbe.copy_replace_vertices(brep, &edge2, &v12, &new_v1);
                                        let ctx = self.my_context.as_mut().unwrap();
                                        ctx.replace(brep, &edge2, &new_e);
                                        sewd.set_edge(&new_e, num2 + dnum1);
                                        boxes.insert(
                                            (new_e.ptr_id(), new_e.location),
                                            b2.clone(),
                                        ); // update boxes
                                        edge2 = new_e;
                                        akey1 = 2;
                                    }
                                    // OCCT L1720-1741.
                                    if p02.distance(pv12) < tol_v2 {
                                        tol_v2 += p02.distance(pv12);
                                        let mut b = BRepBuilder::new();
                                        b.update_vertex_tolerance(brep, new_v2.clone(), tol_v2);
                                        if v12.orientation == new_v2.orientation {
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v12, &new_v2);
                                        } else {
                                            let new_v2r = shape_oriented(&new_v2, Orientation::Reversed);
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v12, &new_v2r);
                                        }
                                        nb_replaced += 1; // gka 06.09.04
                                        let new_e = sbe.copy_replace_vertices(brep, &edge2, &new_v2, &v22);
                                        let ctx = self.my_context.as_mut().unwrap();
                                        ctx.replace(brep, &edge2, &new_e);
                                        sewd.set_edge(&new_e, num2 + dnum1);
                                        boxes.insert(
                                            (new_e.ptr_id(), new_e.location),
                                            b2.clone(),
                                        ); // update boxes
                                        edge2 = new_e;
                                        akey2 = 1;
                                    }
                                    // OCCT L1742-1763.
                                    if p02.distance(pv22) < tol_v2 {
                                        tol_v2 += p02.distance(pv22);
                                        let mut b = BRepBuilder::new();
                                        b.update_vertex_tolerance(brep, new_v2.clone(), tol_v2);
                                        if v22.orientation == new_v2.orientation {
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v22, &new_v2);
                                        } else {
                                            let new_v2r = shape_oriented(&new_v2, Orientation::Reversed);
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &v22, &new_v2r);
                                        }
                                        nb_replaced += 1; // gka 06.09.04
                                        let new_e = sbe.copy_replace_vertices(brep, &edge2, &v12, &new_v2);
                                        let ctx = self.my_context.as_mut().unwrap();
                                        ctx.replace(brep, &edge2, &new_e);
                                        sewd.set_edge(&new_e, num2 + dnum1);
                                        boxes.insert(
                                            (new_e.ptr_id(), new_e.location),
                                            b2.clone(),
                                        ); // update boxes
                                        edge2 = new_e;
                                        akey2 = 2;
                                    }
                                    // OCCT L1764: split
                                    let mut dnum2 = 0i32;
                                    let mut numseg2 = num2 + dnum1;
                                    // OCCT L1765-1815.
                                    if akey1 == 0 && akey2 > 0 {
                                        if self.split_edge1(
                                            brep, sewd, face, num2 + dnum1, p21, &new_v1, tol_v1,
                                            &mut boxes,
                                        ) {
                                            *nb_split += 1;
                                            dnum2 = 1;
                                            // numseg2=num2+dnum1+1;
                                            numseg2 = num2 + dnum1;
                                        }
                                    }
                                    if akey1 > 0 && akey2 == 0 {
                                        if self.split_edge1(
                                            brep, sewd, face, num2 + dnum1, p22, &new_v2, tol_v2,
                                            &mut boxes,
                                        ) {
                                            *nb_split += 1;
                                            dnum2 = 1;
                                            // numseg2=num2+dnum1;
                                            numseg2 = num2 + dnum1 + 1;
                                        }
                                    }
                                    if akey1 == 0 && akey2 == 0 {
                                        if self.split_edge1(
                                            brep, sewd, face, num2 + dnum1, p21, &new_v1, tol_v1,
                                            &mut boxes,
                                        ) {
                                            *nb_split += 1;
                                            dnum2 = 1;
                                        }
                                        tmp_e = sewd.edge(num2 + dnum1);
                                        let mut a = 0.0f64;
                                        let mut b = 0.0f64;
                                        let mut c2d: Option<Curve2d> = None;
                                        sae.pcurve_face(brep, &tmp_e, face, &mut c2d, &mut a, &mut b, false);
                                        if (a - p22) * (b - p22) > 0.0 {
                                            // p22 - external for [a,b] => split next edge
                                            if self.split_edge1(
                                                brep, sewd, face, num2 + dnum1 + dnum2, p22, &new_v2,
                                                tol_v2, &mut boxes,
                                            ) {
                                                *nb_split += 1;
                                                numseg2 = num2 + dnum1 + dnum2;
                                                dnum2 += 1;
                                            }
                                        } else {
                                            if self.split_edge1(
                                                brep, sewd, face, num2 + dnum1, p22, &new_v2, tol_v2,
                                                &mut boxes,
                                            ) {
                                                *nb_split += 1;
                                                dnum2 += 1;
                                                numseg2 = num2 + dnum1 + 1;
                                            }
                                        }
                                    }
                                    // remove segment
                                    // OCCT L1816-1819.
                                    sewd.remove(numseg2);
                                    sewd.remove(numseg1);
                                    *nb_removed += 2;
                                    // num1--;
                                    // break;
                                }
                            }
                        }
                    }
                }
                num2 += 1;
            }
            num1 += 1;
        }
        // OCCT L1829-1830.
        is_done = (*nb_split != 0 || *nb_cut != 0 || nb_replaced != 0 || *nb_removed != 0);
        is_done
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_IntersectionTool.cxx L1835-2530 — FixIntersectingWires.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_IntersectionTool::FixIntersectingWires (cxx L1835-2530).
    pub fn fix_intersecting_wires(&mut self, brep: &mut BRep, face: &mut Shape) -> bool {
        // OCCT L1837-1840.
        if self.my_context.is_none() || face.is_null() {
            return false;
        }
        // TopoDS_Shape S = context->Apply(face);
        // TopoDS_Shape SF = TopoDS::Face(S);
        // OCCT L1843-1858.
        let mut sf = face.clone();
        let ori = face.orientation;
        let mut seq_wir: Vec<Shape> = Vec::new();
        let mut seq_nm_shapes: Vec<Shape> = Vec::new();
        for v in iter_subshapes(brep, &sf, false, true) {
            if v.shape_type() != ShapeType::Wire
                || (v.orientation != Orientation::Forward && v.orientation != Orientation::Reversed)
            {
                seq_nm_shapes.push(v);
                continue;
            }
            seq_wir.push(v);
        }
        // OCCT L1859-1862.
        if (seq_wir.len() as i32) < 2 {
            return false; // gka 06.09.04
        }

        // OCCT L1864-1872.
        let mut max_tol_vert = 0.0f64;
        for exp in topexp_explorer(brep, &sf, ShapeType::Vertex) {
            let tol_v = brep_tool_tolerance(&exp);
            max_tol_vert = max_tol_vert.max(tol_v);
        }
        let mut is_done = false; // gka 06.09.04
        let sae = ShapeAnalysisEdge::new();
        let sas = ShapeAnalysisSurface::new(brep_tool_surface(brep, face).unwrap());

        // OCCT L1874-1886: precompute edge boxes for all wires.
        let mut a_seq_wir_edge_boxes: Vec<ShapeBoxes2d> = Vec::new();
        let mut a_seq_wir_boxes: Vec<BndBox2d> = Vec::new();
        for n in 1..=(seq_wir.len() as i32) {
            let a_wire = seq_wir[(n - 1) as usize].clone();
            let mut a_sewd = WireData::new_from_wire(brep, &a_wire, true, false);
            let mut a_boxes: ShapeBoxes2d = HashMap::new();
            let a_total_box = create_boxes2d(brep, &mut a_sewd, face, &mut a_boxes);
            a_seq_wir_edge_boxes.push(a_boxes);
            a_seq_wir_boxes.push(a_total_box);
        }

        // OCCT L1888-2505.
        for n1 in 1..=((seq_wir.len() as i32) - 1) {
            let mut wire1 = seq_wir[(n1 - 1) as usize].clone();
            let mut sewd1 = WireData::new_from_wire(brep, &wire1, true, false);
            for n2 in (n1 + 1)..=(seq_wir.len() as i32) {
                let mut wire2 = seq_wir[(n2 - 1) as usize].clone();
                let mut sewd2 = WireData::new_from_wire(brep, &wire2, true, false);
                let a_box1 = a_seq_wir_boxes[(n1 - 1) as usize].clone();
                let a_box2 = a_seq_wir_boxes[(n2 - 1) as usize].clone();
                // OCCT L1902-1905.
                if !a_box1.is_void() && !a_box2.is_void() && a_box1.is_out_box(&a_box2) {
                    continue;
                }
                // OCCT L1906-1909: detect possible intersections:
                let mut nb_modif = 0i32;
                let mut nb_replaced = 0i32; // gka 06.09.04
                let mut has_modif_wire = false; // gka 06.09.04
                let mut num1 = 1i32;
                while num1 <= sewd1.nb_edges() && nb_modif < 30 {
                    // for each edge from first wire
                    // OCCT L1913.
                    let mut edge1 = sewd1.edge(num1); // gka 06.09.04

                    let mut num2 = 1i32;
                    while num2 <= sewd2.nb_edges() && nb_modif < 30 {
                        // for each edge from second wire
                        // OCCT L1918.
                        let mut edge2 = sewd2.edge(num2);
                        // OCCT L1919-1922.
                        if edge1.is_same(&edge2) {
                            num2 += 1;
                            continue;
                        }
                        // OCCT L1923-1926.
                        if brep_tool_degenerated(&edge1) || brep_tool_degenerated(&edge2) {
                            num2 += 1;
                            continue;
                        }
                        // OCCT L1927-1930.
                        let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                        if !boxes1.contains_key(&(edge1.ptr_id(), edge1.location)) {
                            num2 += 1;
                            continue;
                        }
                        let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                        if !boxes2.contains_key(&(edge2.ptr_id(), edge2.location)) {
                            num2 += 1;
                            continue;
                        }
                        // OCCT L1931-1933.
                        let b1 = a_seq_wir_edge_boxes[(n1 - 1) as usize]
                            .get(&(edge1.ptr_id(), edge1.location))
                            .unwrap()
                            .clone();
                        let b2 = a_seq_wir_edge_boxes[(n2 - 1) as usize]
                            .get(&(edge2.ptr_id(), edge2.location))
                            .unwrap()
                            .clone();
                        if !b1.is_out_box(&b2) {
                            // intersection is possible...
                            // OCCT L1936-1945.
                            let mut a1 = 0.0f64;
                            let mut b1p = 0.0f64;
                            let mut a2 = 0.0f64;
                            let mut b2p = 0.0f64;
                            let mut crv1: Option<Curve2d> = None;
                            let mut crv2: Option<Curve2d> = None;
                            if !sae.pcurve_face(brep, &edge1, face, &mut crv1, &mut a1, &mut b1p, false) {
                                num2 += 1;
                                continue; // return false; gka 06.09.04
                            }
                            if !sae.pcurve_face(brep, &edge2, face, &mut crv2, &mut a2, &mut b2p, false) {
                                num2 += 1;
                                continue; // return false;gka 06.09.04
                            }
                            let (b1v, b2v) = (b1p, b2p);
                            // OCCT L1946-1951.
                            let tolint = 1.0e-10f64;
                            let crv1 = crv1.unwrap();
                            let crv2 = crv2.unwrap();
                            let d1 = Domain::bounded(
                                Curve2dEval::point_at(&crv1, a1),
                                a1,
                                tolint,
                                Curve2dEval::point_at(&crv1, b1v),
                                b1v,
                                tolint,
                            );
                            let d2 = Domain::bounded(
                                Curve2dEval::point_at(&crv2, a2),
                                a2,
                                tolint,
                                Curve2dEval::point_at(&crv2, b2v),
                                b2v,
                                tolint,
                            );
                            let mut inter = Geom2dIntGInterGap::new();
                            inter.perform(&crv1, &d1, &crv2, &d2, tolint, tolint);
                            // OCCT L1952-1955.
                            if !inter.is_done() {
                                num2 += 1;
                                continue;
                            }
                            // OCCT L1957: intersection is point
                            if inter.nb_points() > 0 && inter.nb_points() < 3 {
                                // OCCT L1959-1961.
                                let mut ip = IntersectionPoint::empty();
                                let mut tr1 = Transition::empty();
                                let mut tr2 = Transition::empty();
                                select_int_pnt(&inter, &mut ip, &mut tr1, &mut tr2);
                                // OCCT L1962-1994.
                                if tr1.position_on_curve() == Position::Middle
                                    && tr2.position_on_curve() == Position::Middle
                                {
                                    // create new vertex and split both edges
                                    let param1 = ip.param_on_first();
                                    let param2 = ip.param_on_second();
                                    let pi1 = get_point_on_edge(brep, &edge1, &sas, &crv1, param1);
                                    let pi2 = get_point_on_edge(brep, &edge2, &sas, &crv2, param2);
                                    let p0 = glam::DVec3::new(
                                        (pi1.x + pi2.x) / 2.0,
                                        (pi1.y + pi2.y) / 2.0,
                                        (pi1.z + pi2.z) / 2.0,
                                    );
                                    let tol_v =
                                        ((pi1.distance(pi2) / 2.0) * 1.00001f64).max(CONFUSION);
                                    let mut b = BRepBuilder::new();
                                    let made = b.add_vertex(brep, p0, tol_v);
                                    let v = made;
                                    max_tol_vert = max_tol_vert.max(tol_v);
                                    let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                    let is_split_edge2 = self.split_edge1(
                                        brep, &mut sewd2, face, num2, param2, &v, tol_v, boxes2,
                                    );
                                    if is_split_edge2 {
                                        nb_modif += 1;
                                        num2 -= 1;
                                    }
                                    let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                                    if self.split_edge1(
                                        brep, &mut sewd1, face, num1, param1, &v, tol_v, boxes1,
                                    ) {
                                        nb_modif += 1;
                                        num1 -= 1;
                                        break;
                                    }
                                    if is_split_edge2 {
                                        num2 += 1;
                                        continue;
                                    }
                                }
                                // OCCT L1995-2014.
                                if tr1.position_on_curve() == Position::Middle
                                    && tr2.position_on_curve() != Position::Middle
                                {
                                    // find needed vertex from edge2 and split edge1 using it
                                    let param1 = ip.param_on_first();
                                    let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                                    if self.find_vert_and_split_edge(
                                        brep,
                                        param1,
                                        &edge1,
                                        &edge2,
                                        &crv1,
                                        &mut max_tol_vert,
                                        &mut num1,
                                        &mut sewd1,
                                        face,
                                        boxes1,
                                        true,
                                    ) {
                                        nb_modif += 1;
                                        break;
                                    }
                                }
                                // OCCT L2015-2034.
                                if tr1.position_on_curve() != Position::Middle
                                    && tr2.position_on_curve() == Position::Middle
                                {
                                    // find needed vertex from edge1 and split edge2 using it
                                    let param2 = ip.param_on_second();
                                    let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                    if self.find_vert_and_split_edge(
                                        brep,
                                        param2,
                                        &edge2,
                                        &edge1,
                                        &crv2,
                                        &mut max_tol_vert,
                                        &mut num2,
                                        &mut sewd2,
                                        face,
                                        boxes2,
                                        true,
                                    ) {
                                        nb_modif += 1;
                                        num2 += 1;
                                        continue;
                                    }
                                }
                                // OCCT L2035-2043.
                                if tr1.position_on_curve() != Position::Middle
                                    && tr2.position_on_curve() != Position::Middle
                                {
                                    // union vertexes
                                    let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                    if self.union_vertexes(
                                        brep, &mut sewd2, &mut edge1, &mut edge2, num2, boxes2, &b2,
                                    ) {
                                        nb_replaced += 1; // gka 06.09.04
                                    }
                                }
                            }
                            // OCCT L2045.
                            has_modif_wire = has_modif_wire || nb_modif != 0 || nb_replaced != 0;
                            // OCCT L2047: intersection is segment
                            if inter.nb_segments() == 1 {
                                // OCCT L2049-2066.
                                let is = inter.segment(1);
                                if is.has_first_point() && is.has_last_point() {
                                    let mut is_modified1 = false;
                                    let mut is_modified2 = false;
                                    let mut new_v = Shape::null();
                                    let mut newtol = 0.0f64;
                                    let ipf = is.first_point();
                                    let p11 = ipf.param_on_first();
                                    let p21 = ipf.param_on_second();
                                    let ipl = is.last_point();
                                    let p12 = ipl.param_on_first();
                                    let p22 = ipl.param_on_second();
                                    let pnt11 = get_point_on_edge(brep, &edge1, &sas, &crv1, p11);
                                    let pnt12 = get_point_on_edge(brep, &edge1, &sas, &crv1, p12);
                                    let pnt21 = get_point_on_edge(brep, &edge2, &sas, &crv2, p21);
                                    let pnt22 = get_point_on_edge(brep, &edge2, &sas, &crv2, p22);

                                    // OCCT L2068: analysis for edge1
                                    let v1 = sae.first_vertex(brep, &edge1);
                                    let pv1 = brep_tool_pnt(&v1);
                                    let v2 = sae.last_vertex(brep, &edge1);
                                    let pv2 = brep_tool_pnt(&v2);
                                    let mut dist1 = pnt11.distance(pv1);
                                    let mut dist2 = pnt12.distance(pv1);
                                    let mut maxdist = dist1.max(dist2);
                                    let mut pdist;
                                    if edge1.orientation == Orientation::Reversed {
                                        pdist = (b1v - p11).abs().max((b1v - p12).abs());
                                    } else {
                                        pdist = (a1 - p11).abs().max((a1 - p12).abs());
                                    }
                                    // OCCT L2085-2090.
                                    if maxdist < max_tol_vert || pdist < (b1v - a1).abs() * 0.01 {
                                        newtol = maxdist;
                                        new_v = v1.clone();
                                        is_modified1 = true;
                                    }
                                    // OCCT L2091-2110.
                                    dist1 = pnt11.distance(pv2);
                                    dist2 = pnt12.distance(pv2);
                                    maxdist = dist1.max(dist2);
                                    if edge1.orientation == Orientation::Reversed {
                                        pdist = (a1 - p11).abs().max((a1 - p12).abs());
                                    } else {
                                        pdist = (b1v - p11).abs().max((b1v - p12).abs());
                                    }
                                    if maxdist < max_tol_vert || pdist < (b1v - a1).abs() * 0.01 {
                                        if (is_modified1 && maxdist < newtol) || !is_modified1 {
                                            newtol = maxdist;
                                            new_v = v2.clone();
                                            is_modified1 = true;
                                        }
                                    }
                                    // OCCT L2111-2143.
                                    if is_modified1 {
                                        // cut edge1 and update tolerance NewV
                                        let dista = (a1 - p11).abs() + (a1 - p12).abs();
                                        let distb = (b1v - p11).abs() + (b1v - p12).abs();
                                        let pend;
                                        let cut;
                                        if dista > distb {
                                            pend = a1;
                                        } else {
                                            pend = b1v;
                                        }
                                        if (pend - p11).abs() > (pend - p12).abs() {
                                            cut = p12;
                                        } else {
                                            cut = p11;
                                        }
                                        let mut is_cut_line = false;
                                        if !self.cut_edge(brep, &mut edge1, pend, cut, face, &mut is_cut_line) {
                                            is_modified1 = false;
                                            num2 += 1;
                                            continue;
                                        }
                                        if newtol > brep_tool_tolerance(&new_v) {
                                            let mut b = BRepBuilder::new();
                                            b.update_vertex_tolerance(brep, new_v.clone(), newtol * 1.00001f64);
                                        }
                                    }

                                    // OCCT L2145: analysis for edge2
                                    let v12 = sae.first_vertex(brep, &edge2);
                                    let pv12 = brep_tool_pnt(&v12);
                                    let v22 = sae.last_vertex(brep, &edge2);
                                    let pv22 = brep_tool_pnt(&v22);
                                    dist1 = pnt21.distance(pv12);
                                    dist2 = pnt22.distance(pv12);
                                    maxdist = dist1.max(dist2);
                                    if edge2.orientation == Orientation::Reversed {
                                        pdist = (b2v - p21).abs().max((b2v - p22).abs());
                                    } else {
                                        pdist = (a2 - p21).abs().max((a2 - p22).abs());
                                    }
                                    // OCCT L2161-2166.
                                    if maxdist < max_tol_vert || pdist < (b2v - a2).abs() * 0.01 {
                                        newtol = maxdist;
                                        new_v = v12.clone();
                                        is_modified2 = true;
                                    }
                                    // OCCT L2167-2186.
                                    dist1 = pnt21.distance(pv22);
                                    dist2 = pnt22.distance(pv22);
                                    maxdist = dist1.max(dist2);
                                    if edge2.orientation == Orientation::Reversed {
                                        pdist = (a2 - p21).abs().max((a2 - p22).abs());
                                    } else {
                                        pdist = (b2v - p21).abs().max((b2v - p22).abs());
                                    }
                                    if maxdist < max_tol_vert || pdist < (b2v - a2).abs() * 0.01 {
                                        if (is_modified2 && maxdist < newtol) || !is_modified2 {
                                            newtol = maxdist;
                                            new_v = v22.clone();
                                            is_modified2 = true;
                                        }
                                    }
                                    // OCCT L2187-2219.
                                    if is_modified2 {
                                        // cut edge1 and update tolerance NewV
                                        let dista = (a2 - p21).abs() + (a2 - p22).abs();
                                        let distb = (b2v - p21).abs() + (b2v - p22).abs();
                                        let pend;
                                        let cut;
                                        if dista > distb {
                                            pend = a2;
                                        } else {
                                            pend = b2v;
                                        }
                                        if (pend - p21).abs() > (pend - p22).abs() {
                                            cut = p22;
                                        } else {
                                            cut = p21;
                                        }
                                        let mut is_cut_line = false;
                                        if !self.cut_edge(brep, &mut edge2, pend, cut, face, &mut is_cut_line) {
                                            is_modified2 = false;
                                            num2 += 1;
                                            continue;
                                        }
                                        if newtol > brep_tool_tolerance(&new_v) {
                                            let mut b = BRepBuilder::new();
                                            b.update_vertex_tolerance(brep, new_v.clone(), newtol * 1.00001f64);
                                        }
                                    }

                                    // OCCT L2221-2228.
                                    if is_modified1 || is_modified2 {
                                        // necessary to make intersect with the same pair of the edges
                                        // once again with modified ranges
                                        num2 -= 1;
                                        has_modif_wire = true; // gka 06.09.04
                                        num2 += 1;
                                        continue;
                                    } else {
                                        // OCCT L2229-2450: create new vertex and split edge1 and edge2
                                        // using it
                                        if (p12 - p11).abs() > (b1v - a1).abs() / 2.0
                                            || (p22 - p21).abs() > (b2v - a2).abs() / 2.0
                                        {
                                            // segment is big and we have to split each intersecting edge
                                            // on 3 edges --> middle edge - edge based on segment
                                            // OCCT L2237-2250.
                                            let p01 = glam::DVec3::new(
                                                (pnt11.x + pnt21.x) / 2.0,
                                                (pnt11.y + pnt21.y) / 2.0,
                                                (pnt11.z + pnt21.z) / 2.0,
                                            );
                                            let p02 = glam::DVec3::new(
                                                (pnt12.x + pnt22.x) / 2.0,
                                                (pnt12.y + pnt22.y) / 2.0,
                                                (pnt12.z + pnt22.z) / 2.0,
                                            );
                                            let mut tol_v1 =
                                                pnt11.distance(p01).max(pnt21.distance(p01));
                                            tol_v1 = tol_v1.max(CONFUSION) * 1.00001f64;
                                            let mut tol_v2 =
                                                pnt12.distance(p02).max(pnt22.distance(p02));
                                            tol_v2 = tol_v2.max(CONFUSION) * 1.00001f64;
                                            if tol_v1 > max_tol_vert || tol_v2 > max_tol_vert {
                                                num2 += 1;
                                                continue;
                                            }

                                            // OCCT L2252-2253.
                                            has_modif_wire = true; // gka 06.09.04
                                            let mut new_v1 = Shape::null();
                                            let mut new_v2 = Shape::null();
                                            let mut tmp_e = Shape::null();
                                            // split edge1
                                            let mut akey1 = 0i32;
                                            let mut akey2 = 0i32;
                                            // analysis fo P01
                                            // OCCT L2258-2275.
                                            if p01.distance(pv1) < tol_v1.max(brep_tool_tolerance(&v1)) {
                                                new_v1 = v1.clone();
                                                if tol_v1 > brep_tool_tolerance(&v1) {
                                                    let mut b = BRepBuilder::new();
                                                    b.update_vertex_tolerance(brep, new_v1.clone(), tol_v1);
                                                }
                                                akey1 += 1;
                                            }
                                            if p01.distance(pv2) < tol_v1.max(brep_tool_tolerance(&v2)) {
                                                new_v1 = v2.clone();
                                                if tol_v1 > brep_tool_tolerance(&v2) {
                                                    let mut b = BRepBuilder::new();
                                                    b.update_vertex_tolerance(brep, new_v1.clone(), tol_v1);
                                                }
                                                akey1 += 1;
                                            }
                                            // analysis fo P02
                                            // OCCT L2276-2294.
                                            if p02.distance(pv1) < tol_v2.max(brep_tool_tolerance(&v1)) {
                                                new_v2 = v1.clone();
                                                if tol_v2 > brep_tool_tolerance(&v1) {
                                                    let mut b = BRepBuilder::new();
                                                    b.update_vertex_tolerance(brep, new_v2.clone(), tol_v2);
                                                }
                                                akey2 += 1;
                                            }
                                            if p02.distance(pv2) < tol_v2.max(brep_tool_tolerance(&v2)) {
                                                new_v2 = v2.clone();
                                                if tol_v2 > brep_tool_tolerance(&v2) {
                                                    let mut b = BRepBuilder::new();
                                                    b.update_vertex_tolerance(brep, new_v2.clone(), tol_v2);
                                                }
                                                akey2 += 1;
                                            }
                                            // OCCT L2295-2298.
                                            if akey1 > 1 || akey2 > 1 {
                                                num2 += 1;
                                                continue;
                                            }
                                            // prepare vertices
                                            // OCCT L2299-2307.
                                            if akey1 == 0 {
                                                let mut b = BRepBuilder::new();
                                                new_v1 = b.add_vertex(brep, p01, tol_v1);
                                            }
                                            if akey2 == 0 {
                                                let mut b = BRepBuilder::new();
                                                new_v2 = b.add_vertex(brep, p02, tol_v2);
                                            }
                                            // split
                                            // OCCT L2308-2346.
                                            let mut numseg1 = num1;
                                            if akey1 == 0 && akey2 > 0 {
                                                let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd1, face, num1, p11, &new_v1, tol_v1,
                                                    boxes1,
                                                ) {
                                                    nb_modif += 1;
                                                    numseg1 = num1 + 1;
                                                }
                                            }
                                            if akey1 > 0 && akey2 == 0 {
                                                let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd1, face, num1, p12, &new_v2, tol_v2,
                                                    boxes1,
                                                ) {
                                                    nb_modif += 1;
                                                    numseg1 = num1;
                                                }
                                            }
                                            if akey1 == 0 && akey2 == 0 {
                                                let mut num1split2 = num1; // what edge to split by point 2
                                                let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd1, face, num1, p11, &new_v1, tol_v1,
                                                    boxes1,
                                                ) {
                                                    nb_modif += 1;
                                                    tmp_e = sewd1.edge(num1);
                                                    let mut a = 0.0f64;
                                                    let mut b = 0.0f64;
                                                    let mut c2d: Option<Curve2d> = None;
                                                    sae.pcurve_face(
                                                        brep, &tmp_e, face, &mut c2d, &mut a, &mut b,
                                                        false,
                                                    );
                                                    if (a - p12) * (b - p12) > 0.0 {
                                                        // p12 - external for [a,b] => split next edge
                                                        num1split2 += 1;
                                                    }
                                                }
                                                let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd1, face, num1split2, p12, &new_v2,
                                                    tol_v2, boxes1,
                                                ) {
                                                    nb_modif += 1;
                                                    numseg1 = num1 + 1;
                                                }
                                            }
                                            // OCCT L2347: SegE = sewd1->Edge(numseg1); get edge from segment
                                            let mut seg_e = sewd1.edge(numseg1);
                                            // split edge2
                                            // replace vertices if it is necessary
                                            // OCCT L2349-2399.
                                            let sbe = ShapeBuildEdge;
                                            akey1 = 0;
                                            akey2 = 0;
                                            if p01.distance(pv12) < tol_v1 {
                                                tol_v1 += p01.distance(pv12);
                                                let mut b = BRepBuilder::new();
                                                b.update_vertex_tolerance(brep, new_v1.clone(), tol_v1);
                                                // V12 = NewV1 — the OCCT updates its local
                                                // copy; the later reads use the fresh vertex.
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                let new_e =
                                                    sbe.copy_replace_vertices(brep, &edge2, &new_v1, &v22);
                                                let ctx = self.my_context.as_mut().unwrap();
                                                ctx.replace(brep, &edge2, &new_e);
                                                sewd2.set_edge(&new_e, num2);
                                                boxes2.insert((new_e.ptr_id(), new_e.location), b2.clone());
                                                edge2 = new_e;
                                                akey1 = 1;
                                            }
                                            if p01.distance(pv22) < tol_v1 {
                                                tol_v1 += p01.distance(pv22);
                                                let mut b = BRepBuilder::new();
                                                b.update_vertex_tolerance(brep, new_v1.clone(), tol_v1);
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                let new_e =
                                                    sbe.copy_replace_vertices(brep, &edge2, &v12, &new_v1);
                                                let ctx = self.my_context.as_mut().unwrap();
                                                ctx.replace(brep, &edge2, &new_e);
                                                sewd2.set_edge(&new_e, num2);
                                                boxes2.insert((new_e.ptr_id(), new_e.location), b2.clone());
                                                edge2 = new_e;
                                                akey1 = 2;
                                            }
                                            if p02.distance(pv12) < tol_v2 {
                                                tol_v2 += p02.distance(pv12);
                                                let mut b = BRepBuilder::new();
                                                b.update_vertex_tolerance(brep, new_v2.clone(), tol_v2);
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                let new_e =
                                                    sbe.copy_replace_vertices(brep, &edge2, &new_v2, &v22);
                                                let ctx = self.my_context.as_mut().unwrap();
                                                ctx.replace(brep, &edge2, &new_e);
                                                sewd2.set_edge(&new_e, num2);
                                                boxes2.insert((new_e.ptr_id(), new_e.location), b2.clone());
                                                edge2 = new_e;
                                                akey2 = 1;
                                            }
                                            if p02.distance(pv22) < tol_v2 {
                                                tol_v2 += p02.distance(pv22);
                                                let mut b = BRepBuilder::new();
                                                b.update_vertex_tolerance(brep, new_v2.clone(), tol_v2);
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                let new_e =
                                                    sbe.copy_replace_vertices(brep, &edge2, &v12, &new_v2);
                                                let ctx = self.my_context.as_mut().unwrap();
                                                ctx.replace(brep, &edge2, &new_e);
                                                sewd2.set_edge(&new_e, num2);
                                                boxes2.insert((new_e.ptr_id(), new_e.location), b2.clone());
                                                edge2 = new_e;
                                                akey2 = 2;
                                            }
                                            // split
                                            // OCCT L2400-2439.
                                            let mut numseg2 = num2;
                                            if akey1 == 0 && akey2 > 0 {
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd2, face, num2, p21, &new_v1, tol_v1,
                                                    boxes2,
                                                ) {
                                                    nb_modif += 1;
                                                    numseg2 = num2 + 1;
                                                }
                                            }
                                            if akey1 > 0 && akey2 == 0 {
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd2, face, num2, p22, &new_v2, tol_v2,
                                                    boxes2,
                                                ) {
                                                    nb_modif += 1;
                                                    numseg2 = num2;
                                                }
                                            }
                                            if akey1 == 0 && akey2 == 0 {
                                                let mut num2split2 = num2;
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd2, face, num2, p21, &new_v1, tol_v1,
                                                    boxes2,
                                                ) {
                                                    nb_modif += 1;
                                                    numseg2 = num2 + 1;
                                                    tmp_e = sewd2.edge(num2);
                                                    let mut a = 0.0f64;
                                                    let mut b = 0.0f64;
                                                    let mut c2d: Option<Curve2d> = None;
                                                    sae.pcurve_face(
                                                        brep, &tmp_e, face, &mut c2d, &mut a, &mut b,
                                                        false,
                                                    );
                                                    if (a - p22) * (b - p22) > 0.0 {
                                                        // p22 - external for [a,b] => split next edge
                                                        num2split2 += 1;
                                                    }
                                                }
                                                let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                                if self.split_edge1(
                                                    brep, &mut sewd2, face, num2split2, p22, &new_v2,
                                                    tol_v2, boxes2,
                                                ) {
                                                    nb_modif += 1;
                                                    numseg2 = num2 + 1;
                                                }
                                            }
                                            // OCCT L2440-2449.
                                            tmp_e = sewd2.edge(numseg2);
                                            let seg_box = a_seq_wir_edge_boxes[(n1 - 1) as usize]
                                                .get(&(seg_e.ptr_id(), seg_e.location))
                                                .cloned();
                                            if let Some(sb) = seg_box {
                                                a_seq_wir_edge_boxes[(n2 - 1) as usize]
                                                    .insert((tmp_e.ptr_id(), tmp_e.location), sb);
                                            }
                                            if !sae
                                                .first_vertex(brep, &seg_e)
                                                .is_same(&sae.first_vertex(brep, &tmp_e))
                                            {
                                                // OCCT SegE.Reverse() — TopAbs::Reverse.
                                                seg_e.orientation = match seg_e.orientation {
                                                    Orientation::Forward => Orientation::Reversed,
                                                    Orientation::Reversed => Orientation::Forward,
                                                    o => o,
                                                };
                                            }
                                            let ctx = self.my_context.as_mut().unwrap();
                                            ctx.replace(brep, &tmp_e, &seg_e);
                                            sewd2.set_edge(&seg_e, numseg2);
                                            num1 -= 1;
                                            break;
                                        } else {
                                            // OCCT L2451-2479: split each intersecting edge on two edges
                                            let p0 = glam::DVec3::new(
                                                (pnt11.x + pnt12.x) / 2.0,
                                                (pnt11.y + pnt12.y) / 2.0,
                                                (pnt11.z + pnt12.z) / 2.0,
                                            );
                                            let param1 = (p11 + p12) / 2.0;
                                            let param2 = (p21 + p22) / 2.0;
                                            let pnt10 = get_point_on_edge(brep, &edge1, &sas, &crv1, param1);
                                            let pnt20 = get_point_on_edge(brep, &edge2, &sas, &crv2, param2);
                                            dist1 = pnt11.distance(p0).max(pnt12.distance(pnt10));
                                            dist2 = pnt21.distance(p0).max(pnt22.distance(pnt10));
                                            let mut tol_v = dist1.max(dist2);
                                            tol_v = tol_v.max(pnt10.distance(pnt20)) * 1.00001f64;
                                            let mut b = BRepBuilder::new();
                                            new_v = b.add_vertex(brep, pnt10, tol_v);
                                            max_tol_vert = max_tol_vert.max(tol_v);
                                            has_modif_wire = true;
                                            let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                                            if self.split_edge2(
                                                brep, &mut sewd2, face, num2, p21, p22, &new_v, tol_v,
                                                boxes2,
                                            ) {
                                                nb_modif += 1;
                                                num2 -= 1;
                                            }
                                            let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                                            if self.split_edge2(
                                                brep, &mut sewd1, face, num1, p11, p12, &new_v, tol_v,
                                                boxes1,
                                            ) {
                                                nb_modif += 1;
                                                num1 -= 1;
                                                break;
                                            }
                                        }
                                    }
                                }
                            } // end if(Inter.NbSegments()==1)
                        }
                        num2 += 1;
                    }
                    num1 += 1;
                }
                // OCCT L2486-2503.
                if has_modif_wire {
                    is_done = true;
                    seq_wir[(n1 - 1) as usize] = sewd1.wire(brep);
                    let ctx = self.my_context.as_mut().unwrap();
                    ctx.replace(brep, &wire1, &seq_wir[(n1 - 1) as usize]);
                    wire1 = seq_wir[(n1 - 1) as usize].clone();
                    // recompute boxes for wire1
                    let boxes1 = &mut a_seq_wir_edge_boxes[(n1 - 1) as usize];
                    boxes1.clear();
                    let a_new_box1 = create_boxes2d(brep, &mut sewd1, face, boxes1);
                    a_seq_wir_boxes[(n1 - 1) as usize] = a_new_box1;
                    seq_wir[(n2 - 1) as usize] = sewd2.wire(brep);
                    let ctx = self.my_context.as_mut().unwrap();
                    ctx.replace(brep, &wire2, &seq_wir[(n2 - 1) as usize]);
                    wire2 = seq_wir[(n2 - 1) as usize].clone();
                    // recompute boxes for wire2
                    let boxes2 = &mut a_seq_wir_edge_boxes[(n2 - 1) as usize];
                    boxes2.clear();
                    let a_new_box2 = create_boxes2d(brep, &mut sewd2, face, boxes2);
                    a_seq_wir_boxes[(n2 - 1) as usize] = a_new_box2;
                }
            }
        }

        // OCCT L2507-2528.
        if is_done {
            // update face
            let mut b_builder = BRepBuilder::new();
            let empty_copied = brep.empty_copied(face);
            let mut newface = empty_copied;
            newface.orientation = Orientation::Forward;
            for i in 1..=(seq_wir.len() as i32) {
                let wire = seq_wir[(i - 1) as usize].clone();
                builder_add(brep, &newface, &wire);
            }
            for i in 1..=(seq_nm_shapes.len() as i32) {
                let a_nms = seq_nm_shapes[(i - 1) as usize].clone();
                builder_add(brep, &newface, &a_nms);
            }
            newface.orientation = ori;
            let ctx = self.my_context.as_mut().unwrap();
            ctx.replace(brep, face, &newface);
            *face = newface;
            let _ = &mut b_builder;
        }
        // OCCT L2529.
        is_done
    }
}
