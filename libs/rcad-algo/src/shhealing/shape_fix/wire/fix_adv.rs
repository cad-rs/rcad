//! OCCT `ShapeFix_Wire` — the fixing methods of advanced level
//! (`ShapeFix_Wire.cxx` L1351-4478): `FixReorder(wi)` / `FixSmall(num)` /
//! `FixConnected(num)` / `FixSeam` / `FixShifted` / `FixDegenerated(num)` /
//! `FixSelfIntersectingEdge` / `FixIntersectingEdges(num)` /
//! `FixIntersectingEdges(num1,num2)` / `FixLacking(num,force)` /
//! `FixNotchedEdges` / `FixTails` / `UpdateWire` / `FixDummySeam`.  The impl
//! submodule of [`super::ShapeFixWire`] (the OCCT continued-file
//! convention); the file statics live in `wire_statics.rs`.

use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Curve2d, Curve2dEval, Curve3, CurveEval, Line2d};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation, ShapeType, TShape};

use crate::shhealing::shape_analysis::analysis;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::transfer_parameters_proj::ShapeAnalysisTransferParametersProj;
use crate::shhealing::shape_analysis::wire_order::ShapeAnalysisWireOrder;
use crate::shhealing::shape_build::brep_tool::topexp_explorer;
use crate::shhealing::shape_build::edge::{builder_range_on_face, ShapeBuildEdge};
use crate::shhealing::shape_build::vertex::ShapeBuildVertex;
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::wire::wire_statics::{
    compute_local_deviation, copy_reverse_pcurves, param_on_first, param_on_second,
    pcurve_range_on_face, remove_loop, remove_loop_split, try_bending_pcurve,
    update_edge_uv_points,
};
use crate::shhealing::shape_fix::wire::{
    brep_tool_degenerated, brep_tool_pnt, brep_tool_same_parameter, brep_tool_tolerance,
    shape_free, shape_oriented, REAL_LAST, ShapeFixWire,
};

/// OCCT gp_XY cross product (the `^` operator).
fn cross2d(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

impl ShapeFixWire {
    // OCCT ShapeFix_Wire.cxx L1351-1400 — FixReorder(wi).
    /// OCCT ShapeFix_Wire::FixReorder(wi) (cxx L1351-1400): reorders edges
    /// in the wire as determined by WireOrder that should be filled and
    /// computed before.
    pub fn fix_reorder_ordered(&mut self, brep: &mut BRep, wi: &ShapeAnalysisWireOrder) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() {
            return false;
        }

        let status = wi.status();
        if status == 0 {
            return false;
        }
        if status <= -10 {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
            return false;
        }

        let nb = self.my_analyzer.wire_data().unwrap().nb_edges();
        if nb != wi.nb_edges() {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
            return false;
        }
        // D abord on protege
        for i in 1..=nb {
            if wi.ordered(i) == 0 {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail3);
                return false;
            }
        }

        let mut newedges: Vec<Shape> = Vec::new();
        for i in 1..=nb {
            let e = self.my_analyzer.wire_data().unwrap().edge(wi.ordered(i));
            newedges.push(e);
        }
        for i in 1..=nb {
            let e = newedges[(i - 1) as usize].clone();
            self.my_analyzer.wire_data_mut().unwrap().set_edge(&e, i);
        }

        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
        true
    }

    // OCCT ShapeFix_Wire.cxx L1404-1472 — FixSmall(num, lockvtx, precsmall).
    /// OCCT ShapeFix_Wire::FixSmall(num, lockvtx, precsmall) (cxx
    /// L1404-1472): fixes the Null Length Edge to be removed.
    pub fn fix_small_edge(&mut self, brep: &mut BRep, num: i32, lockvtx: bool, precsmall: f64) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() <= 1 {
            return false;
        }

        // analysis:
        self.my_analyzer.check_small(brep, num, precsmall);
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail) {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
            //: n2    return false;
        }

        if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
            return false;
        }

        // OUI cette edge est NULLE

        if self.my_analyzer.last_check_status(ShapeExtendStatus::Done2) {
            // edge is small, but vertices are not the same..
            if lockvtx || !self.my_topo_mode {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
                return false;
            }
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);
        } else {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
        }

        // action: remove edge
        let n = if num > 0 { num } else { self.nb_edges() };
        if self.base.my_context.is_some() {
            let e = self.my_analyzer.wire_data().unwrap().edge(n);
            if let Some(ctx) = self.base.my_context.as_mut() {
                ctx.remove(brep, &e);
            }
        }
        // Small edge(s) removed
        self.base.send_warning(
            &self
                .my_analyzer
                .wire_data()
                .map(|w| w.edge(n))
                .unwrap_or_else(Shape::null),
            &MessageMsg::from_key("FixAdvWire.FixSmall.MSG0"),
        );
        self.my_analyzer.wire_data_mut().unwrap().remove(n);

        // call FixConnected in the case if vertices of the small edge were
        // not the same
        if self.last_fix_status(ShapeExtendStatus::Done2) {
            let mut sav_last_fix_status = self.my_last_fix_status;
            // #43 rln 20.11.98 S4054 CTS18544 entity 21734 removing last edge
            let n_now = if n <= self.nb_edges() { n } else { 1 };
            self.fix_connected_edge(brep, n_now, precsmall, true);
            if self.last_fix_status(ShapeExtendStatus::Fail) {
                sav_last_fix_status |= encode_status(ShapeExtendStatus::Fail3);
            }
            self.my_last_fix_status = sav_last_fix_status;
        }

        true
    }

    // OCCT ShapeFix_Wire.cxx L1476-1611 — FixConnected(num, prec, update).
    /// OCCT ShapeFix_Wire::FixConnected(num, prec, theUpdateWire) (cxx
    /// L1476-1611): fixes connected edges (preceding and current); forces
    /// the vertices (end of preceding — begin of current) to be the same
    /// one.
    pub fn fix_connected_edge(&mut self, brep: &mut BRep, num: i32, prec: f64, the_update_wire: bool) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() <= 0 {
            return false;
        }

        // analysis
        let prec_used = if prec < 0.0 { self.base.my_max_tol } else { prec };
        self.my_analyzer.check_connected(brep, num, prec_used);
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail) {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
        }
        if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
            return false;
        }

        // action: replacing vertex
        let nb_edges = self.my_analyzer.wire_data().unwrap().nb_edges();
        let n2 = if num > 0 { num } else { nb_edges };
        let n1 = if n2 > 1 { n2 - 1 } else { nb_edges };
        let e1 = self.my_analyzer.wire_data().unwrap().edge(n1);
        let e2 = self.my_analyzer.wire_data().unwrap().edge(n2);

        let sae = ShapeAnalysisEdge::new();
        let mut v1 = sae.last_vertex(brep, &e1);
        let mut v2 = sae.first_vertex(brep, &e2);
        let v_null = Shape::null();
        let mut v = v_null.clone();

        if self.my_analyzer.last_check_status(ShapeExtendStatus::Done1) {
            // absolutely confused
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
            // #40 rln 18.11.98 S4054 BUC60035 entity 2393 (2-nd sub-curve is
            // edge with the same vertex)
            let last_v2 = sae.last_vertex(brep, &e2);
            if v2.is_same(&last_v2) {
                v = v2.clone();
                if self.base.my_context.is_some() {
                    let vo = shape_oriented(&v, v1.orientation);
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &v1, &vo);
                    }
                }
            } else {
                v = v1.clone();
                if self.base.my_context.is_some() {
                    let vo = shape_oriented(&v, v2.orientation);
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &v2, &vo);
                    }
                }
            }
        } else {
            // on moyenne ...
            if self.my_analyzer.last_check_status(ShapeExtendStatus::Done2) {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);
            } else {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done3);
            }
            let sbv = ShapeBuildVertex;
            v = sbv.combine_vertex(brep, &v1, &v2, 1.0001);
            if self.base.my_context.is_some() {
                let ctx = self.base.my_context.as_mut().unwrap();
                let vo1 = shape_oriented(&v, v1.orientation);
                ctx.replace(brep, &v1, &vo1);
                let vo2 = shape_oriented(&v, v2.orientation);
                ctx.replace(brep, &v2, &vo2);
            }
        }

        // replace vertices to a new one
        let sbe = ShapeBuildEdge;
        if self.my_analyzer.wire_data().unwrap().nb_edges() < 2 {
            if shape_free(brep, &e2) && self.my_topo_mode {
                let mut b = BRepBuilder::new();
                let fv = sae.first_vertex(brep, &e2);
                let lv = sae.last_vertex(brep, &e2);
                b.remove_from_edge(brep, e2.clone(), fv);
                b.remove_from_edge(brep, e2.clone(), lv);
                b.add_to_edge(brep, e2.clone(), shape_oriented(&v, Orientation::Forward));
                b.add_to_edge(brep, e2.clone(), shape_oriented(&v, Orientation::Reversed));
            } else {
                let tmp_e = sbe.copy_replace_vertices(brep, &e2, &v, &v);
                self.my_analyzer.wire_data_mut().unwrap().set_edge(&tmp_e, n2);
                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &e2, &tmp_e);
                    }
                }
            }
        } else {
            if shape_free(brep, &e2) && shape_free(brep, &e1) && self.my_topo_mode {
                let mut b = BRepBuilder::new();
                let fv2 = sae.first_vertex(brep, &e2);
                b.remove_from_edge(brep, e2.clone(), fv2);
                b.add_to_edge(brep, e2.clone(), shape_oriented(&v, Orientation::Forward));
                let lv1 = sae.last_vertex(brep, &e1);
                let lv2 = sae.last_vertex(brep, &e2);
                let fv2b = sae.first_vertex(brep, &e2);
                if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done1)
                    || fv2b.is_same(&lv2)
                {
                    b.remove_from_edge(brep, e1.clone(), lv1);
                    b.add_to_edge(brep, e1.clone(), shape_oriented(&v, Orientation::Reversed));
                }
            } else {
                let tmp_e2 = sbe.copy_replace_vertices(brep, &e2, &v, &v_null);
                self.my_analyzer.wire_data_mut().unwrap().set_edge(&tmp_e2, n2);
                if self.base.my_context.is_some() {
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &e2, &tmp_e2);
                    }
                }
                let lv2 = sae.last_vertex(brep, &e2);
                let fv2 = sae.first_vertex(brep, &e2);
                if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done1) || fv2.is_same(&lv2)
                {
                    let tmp_e1 = sbe.copy_replace_vertices(brep, &e1, &v_null, &v);
                    self.my_analyzer.wire_data_mut().unwrap().set_edge(&tmp_e1, n1);
                    if self.base.my_context.is_some() {
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &e1, &tmp_e1);
                        }
                    }
                }
            }
        }

        // Optionally update wire data with context replacements
        if the_update_wire && self.base.my_context.is_some() {
            self.update_wire(brep);
        }

        let _ = &mut v1;
        let _ = &mut v2;
        true
    }

    // OCCT ShapeFix_Wire.cxx L1615-1637 — FixSeam(num).
    /// OCCT ShapeFix_Wire::FixSeam(num) (cxx L1615-1637): fixes a seam edge;
    /// a seam edge has two pcurves, one for forward, one for reversed; the
    /// forward pcurve must be set as first.
    pub fn fix_seam(&mut self, brep: &mut BRep, num: i32) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let mut c1: Option<Curve2d> = None;
        let mut c2: Option<Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        if !self
            .my_analyzer
            .check_seam_full(brep, num, &mut c1, &mut c2, &mut cf, &mut cl)
        {
            return false;
        }

        let mut b = BRepBuilder::new();
        let nb = self.nb_edges();
        let e = self.my_analyzer.wire_data().unwrap().edge(if num > 0 { num } else { nb });
        let face = self.face();
        //: S4136: BRep_Tool::Tolerance(E) — OCCT passes 0.
        if let (Some(c2v), Some(c1v)) = (c2.clone(), c1.clone()) {
            b.update_edge_pcurve_closed(brep, e.clone(), c2v, c1v, face.clone(), cf, cl, 0.0);
        }
        builder_range_on_face(brep, &e, &face, cf, cl);
        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);

        true
    }

    // OCCT ShapeFix_Wire.cxx L1661-2126 — FixShifted.
    /// OCCT ShapeFix_Wire::FixShifted() (cxx L1661-2126): fixes parametric
    /// curves which may be shifted to the whole parametric range of a
    /// closed surface as a result of recomputing from the 3d
    /// representation.
    pub fn fix_shifted(&mut self, brep: &mut BRep) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        // #78 rln 12.03.99 S4135: checking spatial closure with Precision
        let (uclosed, vclosed, is_v_crv_closed, mut v_range) = {
            let precision = self.base.my_precision;
            let uc = self
                .my_analyzer
                .my_surf
                .as_mut()
                .map(|surf| surf.is_u_closed(precision))
                .unwrap_or(false);
            // #67 rln 01.03.99 S4135: ims010.igs entity D11900 (2D contour
            // is 2*PI higher than the V range [-pi/2,p/2])
            let vc = self
                .my_analyzer
                .my_surf
                .as_mut()
                .map(|surf| {
                    surf.is_v_closed(precision)
                        || matches!(surf.surface(), rcad_kernel::geom::Surface3::Sphere(_))
                })
                .unwrap_or(false);

            // PTV 26.06.2002 begin — CTS18546-2.igs entity 2222: the base
            // curve is periodic and the 2dcurve is shifted
            let mut is_v_crv_closed = false;
            let mut v_range = 1.0f64;
            let rev_profile: Option<Curve3> = match self.my_analyzer.surf().map(|s| s.surface()) {
                Some(rcad_kernel::geom::Surface3::Revolution(rev)) => {
                    Some((*rev.profile).clone())
                }
                _ => None,
            };
            if let Some(mut a_base_crv) = rev_profile {
                loop {
                    match &a_base_crv {
                        Curve3::Offset(off) => a_base_crv = (*off.basis).clone(),
                        Curve3::Trimmed(trm) => a_base_crv = (*trm.curve).clone(),
                        _ => break,
                    }
                }
                if CurveEval::is_periodic(&a_base_crv) {
                    v_range = curve_period_of(&a_base_crv);
                    is_v_crv_closed = true;
                }
            }
            // PTV 26.06.2002 end
            (uc, vc, is_v_crv_closed, v_range)
        };
        if !uclosed && !vclosed {
            return false;
        }

        let (mut u_range, suf, sul, svf, svl, sumid, svmid) = {
            let surf = self.my_analyzer.surf().unwrap();
            let (mut su_f, mut su_l, mut sv_f, mut sv_l) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            surf.bounds(&mut su_f, &mut su_l, &mut sv_f, &mut sv_l);
            let ur = if uclosed { (su_l - su_f).abs() } else { REAL_LAST };
            if !is_v_crv_closed {
                v_range = if vclosed { (sv_l - sv_f).abs() } else { REAL_LAST };
            }
            (ur, su_f, su_l, sv_f, sv_l, 0.5 * (su_f + su_l), 0.5 * (sv_f + sv_l))
        };
        let utol = 0.2 * u_range;
        let vtol = 0.2 * v_range;

        let sbwd_oring: WireData = {
            // OCCT: sbwdOring = WireData() (the shared handle); the rcad
            // filter loop reads the original and builds the filtered copy.
            let mut filtered = WireData::new();
            let nb = self.my_analyzer.wire_data().unwrap().nb_edges();
            for i in 1..=nb {
                let e1 = self.my_analyzer.wire_data().unwrap().edge(i);
                let face = self.face();
                let sae0 = ShapeAnalysisEdge::new();
                if brep_tool_degenerated(&e1) && !sae0.has_pcurve_face(brep, &e1, &face) {
                    continue;
                }
                filtered.add_edge(&e1, 0);
            }
            filtered
        };
        let sae = ShapeAnalysisEdge::new();
        let mut sbwd = WireData::new();
        for i in 1..=sbwd_oring.nb_edges() {
            let e1 = sbwd_oring.edge(i);
            sbwd.add_edge(&e1, 0);
        }

        let sbe = ShapeBuildEdge;
        let nb = sbwd.nb_edges();
        let mut end = (nb == 0);
        let mut degstop = false;
        let mut stop = nb;
        let mut degn2 = 0;
        let mut pdeg = DVec3::ZERO;
        // pdn 17.12.98 r0901_ec 38237 to shift wire at 0
        let mut box2d = rcad_kernel::math::bnd::BndBox2d::new();
        let mut n2 = 1;
        let mut n1 = nb;
        while !end {
            if n2 > nb {
                n2 = 1;
            }
            if n2 == stop {
                end = true;
            }

            let e1 = sbwd.edge(n1);
            let e2 = sbwd.edge(n2);

            if brep_tool_degenerated(&e1) || brep_tool_degenerated(&e2) {
                if !degstop {
                    stop = n2;
                    degstop = true;
                }
                n1 = n2;
                n2 += 1;
                continue;
            }

            let v = sae.first_vertex(brep, &e2);
            if v.is_null() {
                n1 = n2;
                n2 += 1;
                continue;
            }

            let p = brep_tool_pnt(&v);

            let mut a1 = 0.0f64;
            let mut b1 = 0.0f64;
            let mut a2 = 0.0f64;
            let mut b2 = 0.0f64;
            let mut c2d1: Option<Curve2d> = None;
            let mut c2d2: Option<Curve2d> = None;

            //: abv 29.08.01: torCuts.sat: distinguish degeneration by U and
            // by V; only the corresponding move is prohibited
            let mut is_deg = 0;
            let mut deg_p1 = DVec2::ZERO;
            let mut deg_p2 = DVec2::ZERO;
            let mut deg_t1 = 0.0f64;
            let mut deg_t2 = 0.0f64;
            {
                // OCCT: std::max(Precision(), BRep_Tool::Tolerance(V)).
                let preci = self.base.my_precision.max(brep_tool_tolerance(&v));
                let surf_mut = self.my_analyzer.my_surf.as_mut().unwrap();
                if surf_mut.degenerated_values(
                    p,
                    preci,
                    &mut deg_p1,
                    &mut deg_p2,
                    &mut deg_t1,
                    &mut deg_t2,
                    true,
                ) {
                    is_deg = if (deg_p1.x - deg_p2.x).abs() > (deg_p1.y - deg_p2.y).abs() {
                        1
                    } else {
                        2
                    };
                }
            }

            // abv 23 Feb 00: UKI60107-6 210: additional check for the
            // near-degenerated case.  smh#15 PRO19800. Check if the surface
            // is a surface of revolution.
            let is_revolution = matches!(
                self.my_analyzer.surf().unwrap().surface(),
                rcad_kernel::geom::Surface3::Revolution(_)
            );
            if is_revolution {
                let face = self.face();
                if is_deg == 0 && !vclosed {
                    if c2d1.is_none()
                        && !sae.pcurve_face(brep, &e1, &face, &mut c2d1, &mut a1, &mut b1, true)
                    {
                        n1 = n2;
                        n2 += 1;
                        continue;
                    }
                    let c = c2d1.as_ref().unwrap();
                    let p1 = DVec2::new(suf, Curve2dEval::point_at(c, b1).y);
                    let p2 = DVec2::new(sul, Curve2dEval::point_at(c, b1).y);
                    let max_tol = self.base.my_max_tol;
                    let surf_mut = self.my_analyzer.my_surf.as_mut().unwrap();
                    if surf_mut.is_degenerated_uv(p1, p2, max_tol, 10.0)
                        && !surf_mut.is_degenerated_uv(
                            Curve2dEval::point_at(c, a1),
                            Curve2dEval::point_at(c, b1),
                            max_tol,
                            10.0,
                        )
                    {
                        // abv 31.07.00: trj4_pm1-ec-214.stp #31274: still
                        // allow work if the edge already exists
                        is_deg = 1;
                    }
                }
                if is_deg == 0 && !uclosed {
                    if c2d1.is_none()
                        && !sae.pcurve_face(brep, &e1, &face, &mut c2d1, &mut a1, &mut b1, true)
                    {
                        n1 = n2;
                        n2 += 1;
                        continue;
                    }
                    let c = c2d1.as_ref().unwrap();
                    let p1 = DVec2::new(Curve2dEval::point_at(c, b1).x, svf);
                    let p2 = DVec2::new(Curve2dEval::point_at(c, b1).x, svl);
                    let max_tol = self.base.my_max_tol;
                    let surf_mut = self.my_analyzer.my_surf.as_mut().unwrap();
                    if surf_mut.is_degenerated_uv(p1, p2, max_tol, 10.0)
                        && !surf_mut.is_degenerated_uv(
                            Curve2dEval::point_at(c, a1),
                            Curve2dEval::point_at(c, b1),
                            max_tol,
                            10.0,
                        )
                    {
                        is_deg = 2;
                    }
                }
            }

            if is_deg != 0 {
                if !degstop {
                    stop = n2;
                    degstop = true;
                }

                if degn2 == 0 {
                    degn2 = n2;
                    pdeg = p;
                } else {
                    if pdeg.distance_squared(p) < self.base.my_precision * self.base.my_precision {
                        degn2 = n2;
                        // if ( stop < n2 ) { stop = n2; degstop = true; }
                    } else {
                        let mut ax1 = 0.0f64;
                        let mut bx1 = 0.0f64;
                        let mut ax2 = 0.0f64;
                        let mut bx2 = 0.0f64;
                        let mut cx1: Option<Curve2d> = None;
                        let mut cx2: Option<Curve2d> = None;
                        let face = self.face();
                        let e_prev = sbwd.edge(if degn2 > 1 { degn2 - 1 } else { nb });
                        let e_degn = sbwd.edge(degn2);
                        let ok = !{
                            (c2d1.is_none()
                                && !sae.pcurve_face(brep, &e1, &face, &mut c2d1, &mut a1, &mut b1, true))
                                || (c2d2.is_none()
                                    && !sae.pcurve_face(brep, &e2, &face, &mut c2d2, &mut a2, &mut b2, true))
                                || !sae.pcurve_face(brep, &e_prev, &face, &mut cx1, &mut ax1, &mut bx1, true)
                                || !sae.pcurve_face(brep, &e_degn, &face, &mut cx2, &mut ax2, &mut bx2, true)
                        };
                        if !ok {
                            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
                            n1 = n2;
                            n2 += 1;
                            continue;
                        }
                        let pd1 = Curve2dEval::point_at(cx1.as_ref().unwrap(), bx1);
                        let pd2 = Curve2dEval::point_at(cx2.as_ref().unwrap(), ax2);
                        let pn1 = Curve2dEval::point_at(c2d1.as_ref().unwrap(), b1);
                        let pn2 = Curve2dEval::point_at(c2d2.as_ref().unwrap(), a2);
                        let mut x = DVec2::ZERO; // shift vector
                        let period;
                        if uclosed {
                            x.x = 1.0;
                            period = u_range;
                        } else {
                            x.y = 1.0;
                            period = v_range;
                        }
                        let rot1 = cross2d(pn1 - pd2, x);
                        let rot2 = cross2d(pd1 - pn2, x);
                        let scld = (pd2 - pd1).dot(x);
                        let scln = (pn2 - pn1).dot(x);
                        if rot1 * rot2 < -PCONFUSION
                            && scld * scln < -PCONFUSION
                            && scln.abs() > 0.1 * period
                            && scld.abs() > 0.1 * period
                            && rot1 * scld > PCONFUSION
                            && rot2 * scln > PCONFUSION
                        {
                            // abv 02 Mar 00: trying more sophisticated
                            // analysis (ie_exhaust-A.stp #37520)
                            let c2d2r = c2d2.as_ref().unwrap();
                            let cx1r = cx1.as_ref().unwrap();
                            let cx2r = cx2.as_ref().unwrap();
                            let c2d1r = c2d1.as_ref().unwrap();
                            let sign = if rot2 > 0.0 { 1.0 } else { -1.0 };
                            let deep1 = (pn2.dot(x))
                                .min((pd1.dot(x)))
                                .min((Curve2dEval::point_at(c2d2r, b2).dot(x)))
                                .min((Curve2dEval::point_at(cx1r, ax1).dot(x)))
                                .min((Curve2dEval::point_at(c2d2r, 0.5 * (a2 + b2)).dot(x)))
                                .min((Curve2dEval::point_at(cx1r, 0.5 * (ax1 + bx1)).dot(x)))
                                * sign;
                            let deep2 = (pn1.dot(x))
                                .max((pd2.dot(x)))
                                .max((Curve2dEval::point_at(c2d1r, a1).dot(x)))
                                .max((Curve2dEval::point_at(cx2r, bx2).dot(x)))
                                .max((Curve2dEval::point_at(c2d1r, 0.5 * (a1 + b1)).dot(x)))
                                .max((Curve2dEval::point_at(cx2r, 0.5 * (ax2 + bx2)).dot(x)))
                                * sign;
                            let deep = deep2 - deep1; // estimated current size of wire by x
                            // pdn 30 Oct 00: trying the correct period
                            // [0,period] (trj5_k1-tc-203.stp #4698)
                            let dx = analysis::adjust_to_period(
                                deep,
                                PCONFUSION,
                                period + PCONFUSION,
                            );
                            let dx = if scld > 0.0 { -dx } else { dx };
                            x *= dx;
                            // OCCT: Shift.SetTranslation(x); the transformed
                            // pcurve replaces the stored one.
                            let mut k = degn2;
                            loop {
                                if k > nb {
                                    k = 1;
                                }
                                if k == n2 {
                                    break;
                                }
                                let edge = sbwd.edge(k);
                                let face = self.face();
                                let mut cxx: Option<Curve2d> = None;
                                let mut pax = 0.0f64;
                                let mut pbx = 0.0f64;
                                if !sae.pcurve_face(brep, &edge, &face, &mut cxx, &mut pax, &mut pbx, true)
                                {
                                    k += 1;
                                    continue;
                                }
                                // skl 15.05.2002 for OCC208 (if few edges
                                // have reference to one pcurve)
                                let cxx_new =
                                    rcad_kernel::geom::translate_curve2d(cxx.as_ref().unwrap(), x);
                                sbe.replace_pcurve(brep, &edge, &cxx_new, &face);
                                update_edge_uv_points(brep, &edge, &face);
                                k += 1;
                            }
                            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
                            n1 = n2;
                            n2 += 1;
                            continue;
                        }
                        // degn2 = n2; pdeg = p; // ie_exhaust-A.stp #37520
                    }
                }
                //: abv 29.08.01: torCuts.sat:      continue;
            }

            let face = self.face();
            if (c2d1.is_none() && !sae.pcurve_face(brep, &e1, &face, &mut c2d1, &mut a1, &mut b1, true))
                || (c2d2.is_none()
                    && !sae.pcurve_face(brep, &e2, &face, &mut c2d2, &mut a2, &mut b2, true))
            {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
                n1 = n2;
                n2 += 1;
                continue;
            }
            let p2d1 = Curve2dEval::point_at(c2d1.as_ref().unwrap(), b1);
            let p2d2 = Curve2dEval::point_at(c2d2.as_ref().unwrap(), a2);
            box2d.add_point(p2d1);

            let mut du = 0.0f64;
            let mut dv = 0.0f64;
            if uclosed && is_deg != 1 {
                let dx = (p2d2.x - p2d1.x).abs();
                if dx > u_range - utol {
                    du = analysis::adjust_by_period(p2d2.x, p2d1.x, u_range);
                } else if dx > utol && stop == nb {
                    stop = n2; //: abv 29.08.01: torCuts2.stp
                }
            }
            if vclosed && is_deg != 2 {
                let dy = (p2d2.y - p2d1.y).abs();
                if dy > v_range - vtol {
                    dv = analysis::adjust_by_period(p2d2.y, p2d1.y, v_range);
                } else if dy > vtol && stop == nb {
                    stop = n2;
                }
            }
            if du == 0.0 && dv == 0.0 {
                n1 = n2;
                n2 += 1;
                continue;
            }

            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
            // OCCT: Shift.SetTranslation(gp_Vec2d(du, dv)); skl 15.05.2002
            // for OCC208 (if few edges have reference to one pcurve)
            let c2d2_new =
                rcad_kernel::geom::translate_curve2d(c2d2.as_ref().unwrap(), DVec2::new(du, dv));
            sbe.replace_pcurve(brep, &e2, &c2d2_new, &face);
            update_edge_uv_points(brep, &e2, &face);

            n1 = n2;
            n2 += 1;
        }
        if box2d.is_void() {
            return false; // #3 smh 01.04.99. S4163: Overflow, when box is void.
        }

        let (umin, vmin, umax, vmax) = box2d.get().unwrap_or((0.0, 0.0, 0.0, 0.0));
        if (umin + umax - suf - sul).abs() < u_range
            && (vmin + vmax - svf - svl).abs() < v_range
            && !self.last_fix_status(ShapeExtendStatus::Done)
        {
            return false;
        }

        box2d.set_void();
        for n in 1..=nb {
            let edge = sbwd.edge(n);
            let face = self.face();
            let mut c2d: Option<Curve2d> = None;
            let mut a = 0.0f64;
            let mut b = 0.0f64;
            if !sae.pcurve_face(brep, &edge, &face, &mut c2d, &mut a, &mut b, true) {
                continue;
            }
            box2d.add_point(Curve2dEval::point_at(c2d.as_ref().unwrap(), a));
            box2d.add_point(Curve2dEval::point_at(
                c2d.as_ref().unwrap(),
                0.5 * (a + b),
            ));
        }
        let (umin, vmin, umax, vmax) = box2d.get().unwrap_or((0.0, 0.0, 0.0, 0.0));

        let mut du = 0.0f64;
        let mut dv = 0.0f64;

        if uclosed {
            let umid = 0.5 * (umin + umax);
            // PTV 26.06.2002 xloop torus-apple iges face mode
            du = analysis::adjust_by_period(umid, sumid, u_range);
        }
        if vclosed {
            let vmid = 0.5 * (vmin + vmax);
            // PTV 26.06.2002 xloop torus-apple iges face mode
            dv = analysis::adjust_by_period(vmid, svmid, v_range);
        }

        if du == 0.0 && dv == 0.0 {
            return true;
        }

        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);

        // OCCT: Shift.SetTranslation(gp_Vec2d(du, dv)); skl 15.05.2002 for
        // OCC208 (if few edges have reference to one pcurve)
        let shift = DVec2::new(du, dv);
        for n in 1..=sbwd_oring.nb_edges() {
            let ed = sbwd_oring.edge(n);
            let face = self.face();
            let mut c2d: Option<Curve2d> = None;
            let mut a = 0.0f64;
            let mut b = 0.0f64;
            if !sae.pcurve_face(brep, &ed, &face, &mut c2d, &mut a, &mut b, true) {
                continue;
            }
            let c2d_new = rcad_kernel::geom::translate_curve2d(c2d.as_ref().unwrap(), shift);
            sbe.replace_pcurve(brep, &ed, &c2d_new, &face);
            update_edge_uv_points(brep, &ed, &face);
        }
        true
    }

    // OCCT ShapeFix_Wire.cxx L2130-2204 — FixDegenerated(num).
    /// OCCT ShapeFix_Wire::FixDegenerated(num) (cxx L2130-2204): fixes the
    /// degenerated edge; checks the num-th edge or the point between the
    /// (num-1)-th and num-th edges for a singularity on the supporting
    /// surface.
    pub fn fix_degenerated_edge(&mut self, brep: &mut BRep, num: i32) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        // analysis
        let mut p2d1 = DVec2::ZERO;
        let mut p2d2 = DVec2::ZERO;
        self.my_analyzer
            .check_degenerated_full(brep, num, &mut p2d1, &mut p2d2);
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail1) {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
        }
        //: abv 29.08.01: torHalf2.sat: if the edge was encoded as degenerated
        // but has no pcurve and no singularity is found at that point,
        // remove it
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail2) {
            self.my_analyzer.wire_data_mut().unwrap().remove(num);
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done3);
            return true;
        }
        if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
            return false;
        }

        // action: create the degenerated edge and insert it (or replace)

        let vect2d = p2d2 - p2d1;
        let dir2d = vect2d.normalize_or_zero();
        let line2d = Curve2d::Line(Line2d::new(p2d1, dir2d));

        let mut b = BRepBuilder::new();
        let deg_edge = b.add_edge(brep, None, Shape::null(), Shape::null(), [0.0, 0.0]);
        b.set_edge_degenerated(brep, deg_edge.clone(), true);
        let face = self.face();
        b.update_edge_pcurve(brep, deg_edge.clone(), line2d, face.clone(), CONFUSION);
        builder_range_on_face(brep, &deg_edge, &face, 0.0, vect2d.length());

        let nb = self.my_analyzer.wire_data().unwrap().nb_edges();
        let n2 = if num > 0 { num } else { nb };
        let n1 = if n2 > 1 { n2 - 1 } else { nb };

        let lack = self.my_analyzer.last_check_status(ShapeExtendStatus::Done1);
        let n3 = if lack {
            n2
        } else if n2 < nb {
            n2 + 1
        } else {
            1
        };

        let sae = ShapeAnalysisEdge::new();
        let mut v1 = sae.last_vertex(brep, &self.my_analyzer.wire_data().unwrap().edge(n1));
        let mut v2 = sae.first_vertex(brep, &self.my_analyzer.wire_data().unwrap().edge(n3));

        v1.orientation = Orientation::Forward;
        v2.orientation = Orientation::Reversed;
        b.add_to_edge(brep, deg_edge.clone(), v1);
        b.add_to_edge(brep, deg_edge.clone(), v2);
        let mut deg_edge_o = deg_edge.clone();
        deg_edge_o.orientation = Orientation::Forward;

        if lack {
            self.my_analyzer.wire_data_mut().unwrap().add_edge(&deg_edge_o, n2);
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
        } else {
            self.my_analyzer.wire_data_mut().unwrap().set_edge(&deg_edge_o, n2);
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);
        }

        //  commented to avoid extra messages
        //  SendWarning ( degEdge, Message_Msg ( "FixWire.FixDegenerated.MSG0" ) );

        true
    }

    // OCCT ShapeFix_Wire.cxx L2708-2913 — FixSelfIntersectingEdge(num).
    /// OCCT ShapeFix_Wire::FixSelfIntersectingEdge(num) (cxx L2708-2913):
    /// detects and fixes the self-intersecting pcurve of the edge num.
    pub fn fix_self_intersecting_edge(&mut self, brep: &mut BRep, num: i32) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        // analysis
        let mut points2d: Vec<DVec2> = Vec::new();
        let mut points3d: Vec<DVec3> = Vec::new();
        self.my_analyzer
            .check_self_intersecting_edge_full(brep, num, &mut points2d, &mut points3d);
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail) {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
        }
        if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
            return false;
        }

        // action: increase the tolerance of the vertex

        let nb = self.nb_edges();
        let mut e = self.my_analyzer.wire_data().unwrap().edge(if num > 0 { num } else { nb });

        let sae = ShapeAnalysisEdge::new();
        let v1 = sae.first_vertex(brep, &e);
        let v2 = sae.last_vertex(brep, &e);
        let mut tol1 = brep_tool_tolerance(&v1);
        let mut tol2 = brep_tool_tolerance(&v2);
        let pnt1 = brep_tool_pnt(&v1);
        let pnt2 = brep_tool_pnt(&v2);

        // cycle is to verify the fix in case of RemoveLoop
        let tolfact = 0.1; // factor for shifting by parameter in RemoveLoop
        let mut f2d = 0.0f64;
        let mut l2d = 0.0f64;
        let mut c2d: Option<Curve2d> = None;
        let mut newtol = 0.0f64; // = Precision();

        if self.my_remove_loop_mode < 1 {
            for _iter in 0..30 {
                let mut loop_removed = false;
                let mut prev_first = 0.0f64;
                let mut prev_last = 0.0f64;
                let mut i = 1usize;
                while i <= points2d.len() {
                    let pint = points3d[i - 1];
                    let dist21 = pnt1.distance_squared(pint);
                    let dist22 = pnt2.distance_squared(pint);
                    if dist21 < tol1 * tol1 || dist22 < tol2 * tol2 {
                        i += 1;
                        continue;
                    }
                    newtol = 1.001 * dist21.min(dist22).sqrt(); //: f8

                    //: k3 abv 24 Dec 98: BUC50070 #26682 and #30087: try to
                    // remove the loop
                    if self.my_geom_mode {
                        let face = self.face();
                        if c2d.is_none() {
                            sae.pcurve_face(brep, &e, &face, &mut c2d, &mut f2d, &mut l2d, false);
                        }
                        let firstpar = param_on_first(&points2d[i - 1]);
                        let lastpar = param_on_second(&points2d[i - 1]);
                        if firstpar > prev_first && lastpar < prev_last {
                            i += 1;
                            continue;
                        }
                        let max_tol = self.base.my_max_tol;
                        let preci = self.base.my_precision;
                        let remove3d = self.my_remove_loop_mode == 0;
                        if remove_loop(
                            brep,
                            &mut e,
                            &face,
                            &points2d[i - 1],
                            tolfact,
                            max_tol.min(newtol.max(preci)),
                            remove3d,
                        ) {
                            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done4);
                            loop_removed = true;
                            prev_first = firstpar;
                            prev_last = lastpar;
                            // repeat of fix on that edge required (to be done
                            // by the caller)
                            i += 1;
                            continue;
                        }
                    }
                    if newtol < self.base.my_max_tol {
                        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
                        let mut b = BRepBuilder::new();
                        if dist21 < dist22 {
                            tol1 = newtol;
                            b.update_vertex_tolerance(brep, v1.clone(), newtol);
                        } else {
                            tol2 = newtol;
                            b.update_vertex_tolerance(brep, v2.clone(), newtol);
                        }
                    } else {
                        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
                    }
                    i += 1;
                }

                // after RemoveLoop, check that the self-intersection
                // disappeared
                if loop_removed {
                    let mut pnts2d: Vec<DVec2> = Vec::new();
                    let mut pnts3d: Vec<DVec3> = Vec::new();
                    self.my_analyzer
                        .check_self_intersecting_edge_full(brep, num, &mut pnts2d, &mut pnts3d);
                    if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
                        break;
                    }
                    points3d = pnts3d;
                    points2d = pnts2d;
                    let mut b = BRepBuilder::new();
                    let face = self.face();
                    if let Some(c2d_v) = c2d.clone() {
                        b.update_edge_pcurve(brep, e.clone(), c2d_v, face.clone(), 0.0);
                    }
                    builder_range_on_face(brep, &e, &face, f2d, l2d);
                    // newtol+=Precision();
                } else {
                    break;
                }
            }
        }

        //===============================================
        // RemoveLoopMode = 1 , insert vertex
        //===============================================
        if self.my_remove_loop_mode == 1 {
            // after fixing there will be nb+1 edges
            let mut loop_removed;
            // create a sequence of the resulting edges
            let mut ttss: Vec<Shape> = Vec::new();

            loop_removed = false;
            //: k3 abv 24 Dec 98: BUC50070 #26682 and #30087: try to remove
            // the loop
            if self.my_geom_mode {
                let face = self.face();
                if c2d.is_none() {
                    sae.pcurve_face(brep, &e, &face, &mut c2d, &mut f2d, &mut l2d, false);
                }
                let mut e1: Option<Shape> = None;
                let mut e2v: Option<Shape> = None;
                if remove_loop_split(brep, &e, &face, &points2d[0], &mut e1, &mut e2v) {
                    self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done4);
                    loop_removed = true;
                    if let Some(e1s) = &e1 {
                        ttss.push(e1s.clone());
                        newtol = brep_tool_tolerance(e1s)
                            .max(brep_tool_tolerance(e2v.as_ref().unwrap()));
                    } else {
                        newtol = brep_tool_tolerance(e2v.as_ref().unwrap());
                    }
                }
                ttss.push(e2v.unwrap_or_else(Shape::null));
            }

            if newtol > self.base.my_max_tol {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
            }

            let mut sewd = WireData::new();
            for s in &ttss {
                sewd.add_edge(s, 0);
            }
            let sewd_wire = sewd.wire(brep);
            if self.base.my_context.is_some() {
                if let Some(ctx) = self.base.my_context.as_mut() {
                    ctx.replace(brep, &e, &sewd_wire);
                }
                self.update_wire(brep);
            } else {
                let nb_now = self.nb_edges();
                let n = if num > 0 { num } else { nb_now };
                self.my_analyzer.wire_data_mut().unwrap().remove(n);
                self.my_analyzer.wire_data_mut().unwrap().add_edge(&sewd_wire, n);
            }
            if loop_removed {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done8);
            }
        }

        if self.last_fix_status(ShapeExtendStatus::Done) && !self.base.my_shape.is_null() {
            // Edge was self-intersecting, corrected
            self.base.send_warning(
                &e,
                &MessageMsg::from_key("FixAdvWire.FixIntersection.MSG5"),
            );
        }

        self.last_fix_status(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire.cxx L3617-3973 — FixLacking(num, force).
    /// OCCT ShapeFix_Wire::FixLacking(num, force) (cxx L3617-3973): fixes
    /// the lacking edge; tests if two adjacent edges are disconnected in 2d
    /// (while connected in 3d), and in that case either increases the
    /// tolerance of the vertex or adds a new edge (straight in 2d space) to
    /// close the wire in 2d.
    pub fn fix_lacking_edge(&mut self, brep: &mut BRep, num: i32, force: bool) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        //=============
        // First phase: analysis whether the problem (gap) exists
        let mut p2d1 = DVec2::ZERO;
        let mut p2d2 = DVec2::ZERO;
        self.my_analyzer
            .check_lacking_full(brep, num, if force { self.base.my_precision } else { 0.0 }, &mut p2d1, &mut p2d2);
        if self.my_analyzer.last_check_status(ShapeExtendStatus::Fail) {
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
        }
        if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done) {
            return false;
        }

        //=============
        // Second phase: collection of data necessary for further analysis

        let sbwd_nb = self.my_analyzer.wire_data().unwrap().nb_edges();
        let n2 = if num > 0 { num } else { sbwd_nb };
        let n1 = if n2 > 1 { n2 - 1 } else { sbwd_nb };
        let e1 = self.my_analyzer.wire_data().unwrap().edge(n1);
        let e2 = self.my_analyzer.wire_data().unwrap().edge(n2);

        let sae = ShapeAnalysisEdge::new();
        let v1 = sae.last_vertex(brep, &e1);
        let v2 = sae.first_vertex(brep, &e2);
        let tol = brep_tool_tolerance(&v1).max(brep_tool_tolerance(&v2));

        let prec = self.base.my_precision;
        let dist2d = self.my_analyzer.max_distance_2d();
        let mut inctol = self.my_analyzer.max_distance_3d();

        let face = self.my_analyzer.face().clone();

        let mut p3d1 = DVec3::ZERO;
        let mut p3d2 = DVec3::ZERO;
        let mut tol1 = CONFUSION;
        let mut tol2 = CONFUSION; // SK

        //=============
        //: s2 abv 21 Apr 99: Speculation: try bending pcurves
        let mut bendtol1 = 0.0f64;
        let mut bendtol2 = 0.0f64;
        let mut bendc1: Option<Curve2d> = None;
        let mut bendc2: Option<Curve2d> = None;
        let mut bendf1 = 0.0f64;
        let mut bendl1 = 0.0f64;
        let mut bendf2 = 0.0f64;
        let mut bendl2 = 0.0f64;
        if self.my_geom_mode
            && !brep.is_edge_closed_on_face(&e1, &face)
            && !brep.is_edge_closed_on_face(&e2, &face)
        {
            let p2d = 0.5 * (p2d1 + p2d2);
            let mut ok1 = try_bending_pcurve(
                brep,
                &e1,
                &face,
                p2d,
                e1.orientation == Orientation::Forward,
                &mut bendc1,
                &mut bendf1,
                &mut bendl1,
                &mut bendtol1,
            );
            let mut ok2 = try_bending_pcurve(
                brep,
                &e2,
                &face,
                p2d,
                e2.orientation == Orientation::Reversed,
                &mut bendc2,
                &mut bendf2,
                &mut bendl2,
                &mut bendtol2,
            );
            if ok1 && !ok2 {
                bendtol2 = brep_tool_tolerance(&e2);
                ok1 = try_bending_pcurve(
                    brep,
                    &e1,
                    &face,
                    p2d2,
                    e1.orientation == Orientation::Forward,
                    &mut bendc1,
                    &mut bendf1,
                    &mut bendl1,
                    &mut bendtol1,
                );
            } else if !ok1 && ok2 {
                bendtol1 = brep_tool_tolerance(&e1);
                ok2 = try_bending_pcurve(
                    brep,
                    &e2,
                    &face,
                    p2d1,
                    e2.orientation == Orientation::Forward,
                    &mut bendc2,
                    &mut bendf2,
                    &mut bendl2,
                    &mut bendtol2,
                );
            }
            if !ok1 && !ok2 {
                bendc1 = None;
            }
        }

        //=============
        // Third phase: analyse how to fix the problem

        // selector of solutions
        let mut do_increase = false; // increase tolerance
        let mut do_add_long = false; // add long 3d edge in replacement of a vertex
        let mut do_add_closed = false; // add closed 3d edge
        let mut do_add_degen = false; // add degenerated edge
        let mut do_bend = false; //: s2 bend pcurves

        // if bending is OK with the existing tolerances of the edges, take it
        if bendc1.is_some()
            && bendc2.is_some()
            && ((bendtol1 < brep_tool_tolerance(&e1) && bendtol2 < brep_tool_tolerance(&e2))
                || (inctol < prec && bendtol1 < inctol && bendtol2 < inctol))
        {
            do_bend = true;

            // is it OK just to increase the tolerance (to a value less than
            // preci)?
        } else if inctol < prec {
            do_increase = true;

            // If the increase is not OK or forced, try to find other
            // solutions (adding an edge)
        } else if !brep_tool_degenerated(&e2) && !brep_tool_degenerated(&e1) {
            // analyze the 3d space between the edges: is it enough to add a
            // long 3d edge?
            if self.my_topo_mode {
                let mut c3d: Option<Curve3> = None;
                let mut a = 0.0f64;
                let mut b = 0.0f64;
                if !sae.curve3d(brep, &e1, &mut c3d, &mut a, &mut b, true) {
                    // cannot work
                    self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
                    return false;
                }
                p3d1 = rcad_kernel::geom::CurveEval::point_at(c3d.as_ref().unwrap(), b);
                let dist2d3d1 = p3d1.distance(surface_value_of(self, p2d1));
                if !sae.curve3d(brep, &e2, &mut c3d, &mut a, &mut b, true) {
                    // cannot work
                    self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail1);
                    return false;
                }
                p3d2 = rcad_kernel::geom::CurveEval::point_at(c3d.as_ref().unwrap(), a);
                let dist2d3d2 = p3d2.distance(surface_value_of(self, p2d2));

                tol1 = brep_tool_tolerance(&e1).max(dist2d3d1);
                tol2 = brep_tool_tolerance(&e2).max(dist2d3d2);
                //: c5 abv 26 Feb 98: CTS17806 #44418
                let tol0 = tol1 + tol2;
                let dist3d2 = p3d1.distance_squared(p3d2);

                // is it OK to add a long 3d edge?
                if !self.my_analyzer.last_check_status(ShapeExtendStatus::Done2)
                    //: 81 abv 20 Jan 98: don`t add back-going edges (zigzags)
                    && dist3d2 > 1.25 * tol0 * tol0
                    && (force || dist3d2 > prec * prec || inctol > self.base.my_max_tol)
                {
                    do_add_long = true;
                }
            }

            //: h6 abv 25 Jun 98: BUC40132 6361: try to increase tol up to
            // MaxTol if not add
            if !do_add_long
                && inctol < self.base.my_max_tol
                && !self
                    .my_analyzer
                    .my_surf
                    .as_mut()
                    .unwrap()
                    .is_degenerated_uv(p2d1, p2d2, 2.0 * tol, 10.0)
            {
                //: p7
                if bendc1.is_some()
                    && bendc2.is_some()
                    && bendtol1 < inctol
                    && bendtol2 < inctol
                {
                    do_bend = true;
                } else {
                    do_increase = true;
                }
            } else if !do_add_long {
                // else try to add either a degenerated or a closed edge
                let pv = (brep_tool_pnt(&v1) + brep_tool_pnt(&v2)) * 0.5;
                let pm = surface_value_of(self, 0.5 * (p2d1 + p2d2));

                let dist = pv.distance(pm);
                if dist <= tol {
                    do_add_degen = true;
                } else if self.my_topo_mode {
                    do_add_closed = true;
                } else if dist <= self.base.my_max_tol {
                    //: r7 abv 12 Apr 99: t3d_opt.stp #14245 after S4136
                    do_add_degen = true;
                    do_increase = true;
                    inctol = dist;
                }
            }
        } else if !brep_tool_degenerated(&e2) && brep_tool_degenerated(&e1) {
            // create the new degenerated edge and replace E1 by the new edge
        } else if brep_tool_degenerated(&e2) && !brep_tool_degenerated(&e1) {
            // create the new degenerated edge and replace E2 by the new edge
        }

        //=============
        // Third phase - do the fixes
        let mut b = BRepBuilder::new();

        // add edge
        if do_add_long || do_add_degen || do_add_closed {
            // construct the new vertices
            let new_v1: Shape;
            let new_v2: Shape;
            if do_add_long {
                let nv1 = b.add_vertex(brep, p3d1, CONFUSION);
                let mut nv1r = nv1.clone();
                nv1r.orientation = Orientation::Reversed;
                let nv2 = b.add_vertex(brep, p3d2, CONFUSION);
                b.update_vertex_tolerance(brep, nv1r.clone(), 1.001 * tol1);
                b.update_vertex_tolerance(brep, nv2.clone(), 1.001 * tol2);
                new_v1 = nv1r;
                new_v2 = nv2;
            } else {
                new_v1 = v1.clone();
                new_v2 = v2.clone();
            }

            // prepare the new edge
            let edge = b.add_edge(brep, None, Shape::null(), Shape::null(), [0.0, 0.0]);
            if do_add_degen {
                b.set_edge_degenerated(brep, edge.clone(), true); // sln: do it before adding curve
            }
            let v12 = p2d2 - p2d1;
            let the_line2d = Curve2d::Line(Line2d::new(p2d1, v12.normalize_or_zero()));
            b.update_edge_pcurve(brep, edge.clone(), the_line2d, face.clone(), CONFUSION);
            builder_range_on_face(brep, &edge, &face, 0.0, dist2d);
            b.add_to_edge(brep, edge.clone(), shape_oriented(&new_v1, Orientation::Forward));
            b.add_to_edge(brep, edge.clone(), shape_oriented(&new_v2, Orientation::Reversed));
            let sbe = ShapeBuildEdge;
            if !do_add_degen && !sbe.build_curve3d(brep, &edge) {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail3);
                return false;
            }

            // if a long edge is added, replace the vertices of the adjacent
            // edges
            if do_add_long {
                // replace the 1st edge (n1==n2 - special case: the wire
                // consists of one edge)
                let v_null = Shape::null();
                let first_arg = if n1 == n2 { &new_v2 } else { &v_null };
                let edge1 = sbe.copy_replace_vertices(brep, &e1, first_arg, &new_v1);
                self.my_analyzer.wire_data_mut().unwrap().set_edge(&edge1, n1);
                if self.base.my_context.is_some() {
                    let ctx = self.base.my_context.as_mut().unwrap();
                    ctx.replace(brep, &e1, &edge1);
                    // actually, this will occur only in the context of a
                    // single face; hence, recording to ReShape is rather for
                    // tracking modifications than for keeping sharing
                    ctx.replace(brep, &v1, &shape_oriented(&new_v1, v1.orientation));
                    if !v1.is_same(&v2) {
                        ctx.replace(brep, &v2, &shape_oriented(&new_v2, v2.orientation));
                    }
                }
                // replace the 2nd edge
                if n1 != n2 {
                    let edge2 = sbe.copy_replace_vertices(brep, &e2, &new_v2, &v_null);
                    self.my_analyzer.wire_data_mut().unwrap().set_edge(&edge2, n2);
                    if self.base.my_context.is_some() {
                        let ctx = self.base.my_context.as_mut().unwrap();
                        ctx.replace(brep, &e2, &edge2);
                    }
                }
                if self.base.my_context.is_some() {
                    self.update_wire(brep);
                }
            }

            // insert the new edge
            if do_add_degen {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done3);
            } else if !do_add_long {
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done4);
            }
            self.my_analyzer.wire_data_mut().unwrap().add_edge(&edge, n2);
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);
        }
        // else try to increase tol up to MaxTol
        else if inctol > tol && inctol < self.base.my_max_tol {
            if bendc1.is_some() && bendc2.is_some() && bendtol1 < inctol && bendtol2 < inctol {
                do_bend = true;
            } else {
                do_increase = true;
            }
        }

        // bend pcurves
        if do_bend {
            //: s2 abv 21 Apr 99
            if let Some(bc1) = bendc1.clone() {
                b.update_edge_pcurve(brep, e1.clone(), bc1, face.clone(), bendtol1);
            }
            builder_range_on_face(brep, &e1, &face, bendf1, bendl1);
            if let Some(bc2) = bendc2.clone() {
                b.update_edge_pcurve(brep, e2.clone(), bc2, face.clone(), bendtol2);
            }
            builder_range_on_face(brep, &e2, &face, bendf2, bendl2);
            b.update_vertex_tolerance(brep, sae.first_vertex(brep, &e1), bendtol1);
            b.update_vertex_tolerance(brep, sae.last_vertex(brep, &e1), bendtol1);
            b.update_vertex_tolerance(brep, sae.first_vertex(brep, &e2), bendtol2);
            b.update_vertex_tolerance(brep, sae.last_vertex(brep, &e2), bendtol2);
            //: s3 abv 22 Apr 99: PRO7187 #11534: self-intersection not
            // detected until the curve is bent (!)
            self.fix_self_intersecting_edge(brep, n1);
            self.fix_self_intersecting_edge(brep, n2);
            self.fix_intersecting_edges(brep, n2); // skl 24.04.2003 for OCC58
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done5);
        }

        // increase the vertex tolerance
        if do_increase {
            b.update_vertex_tolerance(brep, v1.clone(), 1.001 * inctol);
            b.update_vertex_tolerance(brep, v2.clone(), 1.001 * inctol);
            self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
        }

        if self.last_fix_status(ShapeExtendStatus::Done) {
            return true;
        }

        self.my_last_fix_status |= encode_status(ShapeExtendStatus::Fail2);
        false
    }

    // OCCT ShapeFix_Wire.cxx L3977-4104 — FixNotchedEdges.
    /// OCCT ShapeFix_Wire::FixNotchedEdges() (cxx L3977-4104).
    pub fn fix_notched_edges(&mut self, brep: &mut BRep) -> bool {
        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let face = self.face();
        if self.base.my_context.is_some() {
            self.update_wire(brep);
        }

        let mut i = 1;
        while i <= self.nb_edges() && self.nb_edges() > 2 {
            let mut param = 0.0f64;
            let mut to_remove = 0i32;
            if self
                .my_analyzer
                .check_notched_edges(brep, i, &mut to_remove, &mut param, self.base.my_min_tol)
            {
                let nb_now = self.nb_edges();
                let n2 = if i > 0 { i } else { nb_now };
                let n1 = if n2 > 1 { n2 - 1 } else { nb_now };
                let is_remove_first = n1 == to_remove;
                let to_split = if n2 == to_remove { n1 } else { n2 };
                let split_e = self.my_analyzer.wire_data().unwrap().edge(to_split);
                let sae = ShapeAnalysisEdge::new();
                let mut c2d: Option<Curve2d> = None;
                let mut a = 0.0f64;
                let mut b = 0.0f64;
                sae.pcurve_face(brep, &split_e, &face, &mut c2d, &mut a, &mut b, true);
                let sbe = ShapeBuildEdge;
                let orient = split_e.orientation;

                // check whether the whole edges should be removed - this is
                // the case when the split point coincides with the end of
                // the edge; for closed edges the split point may fall at the
                // other end (see issue #0029780)
                if (param - if is_remove_first { b } else { a }).abs() <= PCONFUSION
                    || (sae.is_closed3d(brep, &split_e)
                        && (param - if is_remove_first { a } else { b }).abs() <= PCONFUSION)
                {
                    self.fix_dummy_seam(brep, n1);
                    // The seam edge is removed from the list. So, need to
                    // step back to avoid missing the edge processing.
                    i -= 1;
                } else {
                    // perform splitting of the edge and adding to wire

                    // pdn check if it is necessary
                    if ((if is_remove_first { a } else { b }) - param).abs() < PCONFUSION {
                        i += 1;
                        continue;
                    }

                    let mut transfer_parameters = ShapeAnalysisTransferParametersProj::new();
                    transfer_parameters.base.my_max_tolerance = self.base.my_max_tol;
                    transfer_parameters.init(brep, &split_e, &face);
                    let (first, last) = if a < b { (a, b) } else { (b, a) };
                    let surf_pt = self
                        .my_analyzer
                        .surf()
                        .unwrap()
                        .value(Curve2dEval::point_at(c2d.as_ref().unwrap(), param));
                    let mut b = BRepBuilder::new();
                    let vnew = b.add_vertex(brep, surf_pt, CONFUSION);
                    let mut we = split_e.clone();
                    we.orientation = Orientation::Forward;
                    let vnew_rev = shape_oriented(&vnew, Orientation::Reversed);
                    let first_v = sae.first_vertex(brep, &we);
                    let new_e1 = sbe.copy_replace_vertices(brep, &we, &first_v, &vnew_rev);
                    sbe.copy_pcurves(brep, &new_e1, &we);
                    transfer_parameters.transfer_range(brep, &new_e1, first, param, true);
                    b.set_edge_same_range(brep, new_e1.clone(), false);
                    b.set_edge_same_parameter(brep, new_e1.clone(), false);
                    let vnew_fwd = shape_oriented(&vnew, Orientation::Forward);
                    let last_v = sae.last_vertex(brep, &we);
                    let new_e2 = sbe.copy_replace_vertices(brep, &we, &vnew_fwd, &last_v);
                    sbe.copy_pcurves(brep, &new_e2, &we);
                    transfer_parameters.transfer_range(brep, &new_e2, param, last, true);
                    b.set_edge_same_range(brep, new_e2.clone(), false);
                    b.set_edge_same_parameter(brep, new_e2.clone(), false);

                    if self.base.my_context.is_some() {
                        let mut wire = b.make_wire(brep);
                        b.add_to_wire(brep, wire.clone(), new_e1.clone());
                        b.add_to_wire(brep, wire.clone(), new_e2.clone());
                        wire = wire;
                        if let Some(ctx) = self.base.my_context.as_mut() {
                            ctx.replace(brep, &we, &wire);
                        }
                    }

                    let mut new_e1 = new_e1;
                    let mut new_e2 = new_e2;
                    new_e1.orientation = orient;
                    new_e2.orientation = orient;
                    if orient == Orientation::Reversed {
                        std::mem::swap(&mut new_e1, &mut new_e2);
                    }

                    let is_remove_last = (n1 == self.nb_edges()) && (n2 == 1);
                    self.my_analyzer
                        .wire_data_mut()
                        .unwrap()
                        .set_edge(&new_e1, to_split);
                    let add_at = if to_split == self.nb_edges() { 0 } else { to_split + 1 };
                    self.my_analyzer.wire_data_mut().unwrap().add_edge(&new_e2, add_at);

                    self.fix_dummy_seam(brep, if is_remove_last { self.nb_edges() } else { to_remove });
                    self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done2);
                }

                i -= 1;
                if self.base.my_context.is_some() {
                    // skl 07.03.2002 for OCC180
                    self.update_wire(brep);
                }
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done1);
            }
            i += 1;
        }
        self.my_status_notches = self.my_last_fix_status;
        self.last_fix_status(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire.cxx L4213-4289 — FixDummySeam(num).
    /// OCCT ShapeFix_Wire::FixDummySeam(num) (cxx L4213-4289).
    pub fn fix_dummy_seam(&mut self, brep: &mut BRep, num: i32) {
        let sae = ShapeAnalysisEdge::new();
        let sbe = ShapeBuildEdge;
        let sbv = ShapeBuildVertex;
        let nb = self.nb_edges();
        let num1 = if num == nb { 1 } else { num + 1 };
        let e1 = self.my_analyzer.wire_data().unwrap().edge(num);
        let e2 = self.my_analyzer.wire_data().unwrap().edge(num1);
        let v1 = sae.first_vertex(brep, &e1);
        let v2 = sae.last_vertex(brep, &e2);
        let vm = sbv.combine_vertex(brep, &v1, &v2, 1.0001);

        // pnd defining if the new pcurves exist
        let to_remove = false;

        // creating the new edge with the pcurves and the new vertex
        let mut vs = sae.first_vertex(brep, &e2);
        if vs.is_same(&v1) || vs.is_same(&v2) {
            vs = vm.clone();
        }
        let new_edge = sbe.copy_replace_vertices(brep, &e2, &vs, &vm);
        copy_reverse_pcurves(brep, &new_edge, &e1, e1.orientation == e2.orientation);
        let mut b = BRepBuilder::new();
        b.set_edge_same_range(brep, new_edge.clone(), false);
        b.set_edge_same_parameter(brep, new_edge.clone(), false);

        if self.base.my_context.is_some() {
            let ctx = self.base.my_context.as_mut().unwrap();
            if to_remove {
                ctx.remove(brep, &e2);
                ctx.remove(brep, &e1);
            } else {
                ctx.replace(brep, &e2, &new_edge);
                let rev = shape_oriented(&new_edge, Orientation::Reversed);
                ctx.replace(brep, &e1, &rev);
            }
            ctx.replace(brep, &v1, &shape_oriented(&vm, v1.orientation));
            ctx.replace(brep, &v2, &shape_oriented(&vm, v2.orientation));
        }

        let nb_now = self.nb_edges();
        let next = if num1 == nb_now { 1 } else { num1 + 1 };
        let prev = if num > 1 { num - 1 } else { nb_now };
        let prev_e = self.my_analyzer.wire_data().unwrap().edge(prev);
        let next_e = self.my_analyzer.wire_data().unwrap().edge(next);

        let v_null = Shape::null();
        let tmp_e1 = sbe.copy_replace_vertices(brep, &prev_e, &v_null, &vm);
        self.my_analyzer.wire_data_mut().unwrap().set_edge(&tmp_e1, prev);
        if self.base.my_context.is_some() {
            if let Some(ctx) = self.base.my_context.as_mut() {
                ctx.replace(brep, &prev_e, &tmp_e1);
            }
        }

        let tmp_e1 = sbe.copy_replace_vertices(brep, &next_e, &vm, &v_null);
        self.my_analyzer.wire_data_mut().unwrap().set_edge(&tmp_e1, next);
        if self.base.my_context.is_some() {
            if let Some(ctx) = self.base.my_context.as_mut() {
                ctx.replace(brep, &next_e, &tmp_e1);
            }
        }

        // removing the edges from the wire
        let (n1r, n2r) = if num < num1 { (num, num1) } else { (num1, num) };
        self.my_analyzer.wire_data_mut().unwrap().remove(n2r);
        self.my_analyzer.wire_data_mut().unwrap().remove(n1r);
    }

    // OCCT ShapeFix_Wire.cxx L4293-4310 — UpdateWire.
    /// OCCT ShapeFix_Wire::UpdateWire (cxx L4293-4310): updates the WireData
    /// if some replacements are made; this is necessary for wires (unlike
    /// other shape types) since one edge can be present in the wire several
    /// times.
    pub fn update_wire(&mut self, brep: &mut BRep) {
        let mut i = 1;
        loop {
            let nb = self.my_analyzer.wire_data().unwrap().nb_edges();
            if i > nb {
                break;
            }
            let e = self.my_analyzer.wire_data().unwrap().edge(i);
            let s = self.context_apply(brep, &e);
            if s.is_same(&e) {
                i += 1;
                continue;
            }
            // OCCT: for (TopExp_Explorer exp(S, TopAbs_EDGE)...)
            // sbwd->Add(exp.Current(), i++); — inserts every sub-edge at the
            // incrementing position.
            let sub_edges = topexp_explorer(brep, &s, ShapeType::Edge);
            for se in &sub_edges {
                self.my_analyzer.wire_data_mut().unwrap().add_edge(se, i);
                i += 1;
            }
            self.my_analyzer.wire_data_mut().unwrap().remove(i);
            // OCCT: sbwd->Remove(i--) — the post-decrement cancels the loop
            // increment, so `i` stays.
        }
    }

    // OCCT ShapeFix_Wire.cxx L4314-4478 — FixTails.
    /// OCCT ShapeFix_Wire::FixTails (cxx L4314-4478).
    pub fn fix_tails(&mut self, brep: &mut BRep) -> bool {
        if self.my_max_tail_width < 0.0 || !self.is_ready() {
            return false;
        }

        self.my_last_fix_status = encode_status(ShapeExtendStatus::Ok);
        if self.base.my_context.is_some() {
            self.update_wire(brep);
        }
        let mut a_check_angle = true;
        let mut a_e_count = self.nb_edges();
        let mut a_en_ns = [a_e_count, 1];
        while a_e_count >= 2 && a_en_ns[1] <= a_e_count {
            let a_es = [
                self.my_analyzer.wire_data().unwrap().edge(a_en_ns[0]),
                self.my_analyzer.wire_data().unwrap().edge(a_en_ns[1]),
            ];
            let mut a_e_parts: [[Option<Shape>; 2]; 2] = [[None, None], [None, None]];
            let max_sine = if a_check_angle {
                self.my_max_tail_angle_sine
            } else {
                -1.0
            };
            // (the four output slots are destructured for the disjoint
            // mutable borrows)
            let [a_parts0, a_parts1] = &mut a_e_parts;
            let [a_p00, a_p01] = a_parts0;
            let [a_p10, a_p11] = a_parts1;
            if !self.my_analyzer.check_tail(
                brep,
                &a_es[0],
                &a_es[1],
                max_sine,
                self.my_max_tail_width,
                self.base.my_max_tol,
                a_p00,
                a_p01,
                a_p10,
                a_p11,
            ) {
                a_en_ns[0] = a_en_ns[1];
                a_en_ns[1] += 1;
                a_check_angle = true;
                continue;
            }

            // Provide not less than 1 edge in the result wire.
            let a_split_counts = [
                if a_e_parts[0][1].is_some() { 1 } else { 0 },
                if a_e_parts[1][1].is_some() { 1 } else { 0 },
            ];
            let a_remove_count = (if a_e_parts[0][0].is_some() { 1 } else { 0 })
                + (if a_e_parts[1][0].is_some() { 1 } else { 0 });
            if a_e_count + a_split_counts[0] + a_split_counts[1] < 1 + a_remove_count {
                a_en_ns[0] = a_en_ns[1];
                a_en_ns[1] += 1;
                a_check_angle = true;
                continue;
            }

            // Split the edges.
            for a_ei in 0..2 {
                if a_split_counts[a_ei] == 0 {
                    continue;
                }

                // Replace the edge by the wire of its parts in the shape.
                let a_e = a_es[a_ei].clone();
                if self.base.my_context.is_some() {
                    let mut b = BRepBuilder::new();
                    let mut a_ewire = b.make_wire(brep);
                    b.add_to_wire(brep, a_ewire.clone(), a_e_parts[a_ei][0].clone().unwrap());
                    b.add_to_wire(brep, a_ewire.clone(), a_e_parts[a_ei][1].clone().unwrap());
                    let a_fe = shape_oriented(&a_e, Orientation::Forward);
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.replace(brep, &a_fe, &a_ewire);
                    }
                    let _ = &mut a_ewire;
                }

                // Replace the edge by its parts in the edge wire.
                let a_orient = a_e.orientation;
                let p0 = a_e_parts[a_ei][0].as_mut().unwrap();
                p0.orientation = a_orient;
                let p1 = a_e_parts[a_ei][1].as_mut().unwrap();
                p1.orientation = a_orient;
                let a_first_pi = if a_orient != Orientation::Reversed { 0 } else { 1 };
                let a_add = if a_ei == 0 || a_en_ns[1] < a_en_ns[0] {
                    0
                } else {
                    a_split_counts[0]
                };
                self.my_analyzer
                    .wire_data_mut()
                    .unwrap()
                    .set_edge(a_e_parts[a_ei][a_first_pi].as_ref().unwrap(), a_en_ns[a_ei] + a_add);
                self.my_analyzer.wire_data_mut().unwrap().add_edge(
                    a_e_parts[a_ei][1 - a_first_pi].as_ref().unwrap(),
                    a_en_ns[a_ei] + 1 + a_add,
                );
            }

            // Remove the tail.
            if a_remove_count == 2 {
                a_check_angle = true;
                self.fix_dummy_seam(
                    brep,
                    a_en_ns[0]
                        + a_split_counts[0]
                        + (if a_en_ns[0] < a_en_ns[1] { 0 } else { a_split_counts[1] }),
                );
                if self.base.my_context.is_some() {
                    self.update_wire(brep);
                }
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done);

                if a_split_counts[0] + a_split_counts[1] == 2 {
                    a_en_ns[0] = a_en_ns[1];
                    a_en_ns[1] += 1;
                    continue;
                }

                if a_split_counts[0] == a_split_counts[1] {
                    a_e_count -= 2;
                    if a_en_ns[1] >= 3 {
                        a_en_ns[0] -= 1;
                        a_en_ns[1] -= 1;
                    } else {
                        a_en_ns[0] = a_e_count;
                        a_en_ns[1] = 1;
                    }
                    a_check_angle = false;
                } else {
                    a_e_count -= 1;
                    if a_split_counts[0] != 0 {
                        a_en_ns[0] = if a_en_ns[0] <= a_e_count { a_en_ns[0] } else { a_e_count };
                    } else if a_en_ns[1] >= 3 {
                        a_en_ns[0] -= 1;
                        a_en_ns[1] -= 1;
                    } else {
                        a_en_ns[0] = a_e_count;
                        a_en_ns[1] = 1;
                    }
                }
            } else {
                a_check_angle = false;
                a_e_count -= 1;
                let a_ri = if a_e_parts[0][0].is_none() { 1 } else { 0 };
                if a_split_counts[a_ri] != 0 {
                    if a_ri == 0 {
                        if a_en_ns[1] >= 3 {
                            a_en_ns[0] -= 1;
                            a_en_ns[1] -= 1;
                        } else {
                            a_en_ns[0] = a_e_count;
                            a_en_ns[1] = 1;
                        }
                    } else {
                        a_en_ns[0] = if a_en_ns[1] > 1 { a_en_ns[0] } else { a_e_count };
                    }
                }
                let rm_at = a_en_ns[a_ri] + (if a_ri != 0 || a_split_counts[0] == 0 { 0 } else { 1 });
                self.my_analyzer.wire_data_mut().unwrap().remove(rm_at);
                if self.base.my_context.is_some() {
                    let fwd = shape_oriented(&a_es[a_ri], Orientation::Forward);
                    if let Some(ctx) = self.base.my_context.as_mut() {
                        ctx.remove(brep, &fwd);
                    }
                    self.update_wire(brep);
                }
                self.my_last_fix_status |= encode_status(ShapeExtendStatus::Done);
            }
        }
        self.my_status_notches = self.my_last_fix_status;
        crate::shhealing::shape_extend::status::decode_status(
            self.my_last_fix_status,
            ShapeExtendStatus::Done,
        )
    }
}

// ---------------------------------------------------------------------------
// File-local helpers.
// ---------------------------------------------------------------------------
fn surface_value_of(fw: &ShapeFixWire, p2d: DVec2) -> DVec3 {
    fw.my_analyzer
        .surf()
        .map(|s| s.value(p2d))
        .unwrap_or(DVec3::ZERO)
}

/// OCCT `Geom_Curve::Period()` for the revolution basis curve.
fn curve_period_of(c: &Curve3) -> f64 {
    let d = CurveEval::default_domain(c);
    d[1] - d[0]
}
