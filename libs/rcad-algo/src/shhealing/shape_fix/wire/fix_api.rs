//! OCCT `ShapeFix_Wire` — the fixing methods of API level
//! (`ShapeFix_Wire.cxx` L297-1343): `Perform` / `FixReorder(bool)` /
//! `FixSmall(bool,prec)` / `FixConnected(prec)` / `FixEdgeCurves` /
//! `FixDegenerated()` / `FixSelfIntersection` / `FixLacking(bool)` /
//! `FixClosed`.  The impl submodule of [`super::ShapeFixWire`] (the OCCT
//! continued-file convention).

use glam::DVec2;
use rcad_kernel::geom::CurveEval;
use rcad_kernel::precision::PCONFUSION;
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::geom::SurfaceEval;
use rcad_kernel::topods::{BRep, BRepBuilder, Orientation, ShapeType, TShape};

use crate::shhealing::shape_analysis::curve::ShapeAnalysisCurve;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::transfer_parameters_proj::ShapeAnalysisTransferParametersProj;
use crate::shhealing::shape_analysis::wire_order::ShapeAnalysisWireOrder;
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_extend::status::{encode_status, ShapeExtendStatus};
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::root::ShapeFixRoot;
use crate::shhealing::shape_fix::shape_fix::MessageProgressRange;
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::shhealing::shape_fix::wire::{brep_tool_pnt, brep_tool_same_parameter, ShapeFixWire};

impl ShapeFixWire {
    // OCCT ShapeFix_Wire.cxx L297-483 — Perform.
    /// OCCT ShapeFix_Wire::Perform(theProgress) (cxx L297-483): performs all
    /// the available fixes; if some fix is turned on or off explicitly by
    /// the Fix..Mode() flag, this fix is either called or not depending on
    /// that flag; else (i.e. if the flag is default) the fix is called
    /// depending on the situation.
    pub fn perform(&mut self, brep: &mut BRep, the_progress: MessageProgressRange) -> bool {
        self.clear_statuses();
        if !self.is_loaded() {
            return false;
        }

        // OCCT: myFixEdge->SetContext(Context()) — the shared handle; the
        // rcad value model passes a copy (the ShapeFix_Edge::SetContext
        // by-value precedent of the W3 tranche 1).
        if let Some(ctx) = self.base.my_context.as_mut() {
            self.my_fix_edge.set_context(ctx.clone());
        }

        let mut fixed = false;

        // FixReorder is first, because as a rule wire is required to be
        // ordered.  We shall analyze the order of edges in the wire and set
        // the appropriate status even if FixReorder should not be called
        // (if it is forbidden).

        let mut sawo = ShapeAnalysisWireOrder::new();
        // OCCT: ReorderOK = (CheckOrder(sawo, myClosedMode) == 0) — the rcad
        // check_order returns the Done flag, i.e. status != 0.
        let reorder_ok = !self.my_analyzer.check_order(
            brep,
            &mut sawo,
            true,
            false,
            self.my_closed_mode,
        );
        let mut reorder_ok = reorder_ok;
        if ShapeFixRoot::need_fix(self.my_fix_reorder_mode, !reorder_ok) {
            if self.fix_reorder(brep, false) {
                fixed = true;
            }
            reorder_ok = !self.status_reorder(ShapeExtendStatus::Fail);
        }

        if the_progress.user_break() {
            return false;
        }

        // FixSmall is allowed to change topology only if mode is set and
        // FixReorder did not failed.
        if ShapeFixRoot::need_fix(self.my_fix_small_mode, self.my_topo_mode) {
            if self.fix_small_all(
                brep,
                !self.my_topo_mode || !reorder_ok,
                self.base.my_min_tol,
            ) != 0
            {
                fixed = true;
                // retry reorder if necessary (after FixSmall it can work
                // better)
                if ShapeFixRoot::need_fix(self.my_fix_reorder_mode, !reorder_ok) {
                    self.fix_reorder(brep, false);
                    reorder_ok = !self.status_reorder(ShapeExtendStatus::Fail);
                }
            }
        }

        if the_progress.user_break() {
            return false;
        }

        if ShapeFixRoot::need_fix(self.my_fix_connected_mode, reorder_ok) {
            if self.fix_connected(brep, -1.0) {
                fixed = true;
            }
        }

        if the_progress.user_break() {
            return false;
        }

        if ShapeFixRoot::need_fix(self.my_fix_edge_curves_mode, true) {
            let sav_fix_shifted_mode = self.my_fix_shifted_mode;
            // turn out FixShifted if reorder not done
            if self.my_fix_shifted_mode == -1 && !reorder_ok {
                self.my_fix_shifted_mode = 0;
            }
            if self.fix_edge_curves(brep) {
                fixed = true;
            }
            self.my_fix_shifted_mode = sav_fix_shifted_mode;
        }

        if the_progress.user_break() {
            return false;
        }

        if ShapeFixRoot::need_fix(self.my_fix_degenerated_mode, true) {
            if self.fix_degenerated(brep) {
                fixed = true; // ?? if ! ReorderOK ??
            }
        }

        if the_progress.user_break() {
            return false;
        }

        // pdn - temporary to test
        if self.my_fix_tail_mode <= 0 && ShapeFixRoot::need_fix(self.my_fix_notched_edges_mode, reorder_ok) {
            fixed |= self.fix_notched_edges(brep);
            if fixed {
                self.fix_shifted(brep); // skl 07.03.2002 for OCC180
            }
        }

        if the_progress.user_break() {
            return false;
        }

        if self.my_fix_tail_mode != 0 {
            if self.fix_tails(brep) {
                fixed = true;
                self.fix_shifted(brep);
            }
        }

        if the_progress.user_break() {
            return false;
        }

        if ShapeFixRoot::need_fix(self.my_fix_self_intersection_mode, self.my_closed_mode) {
            let sav_fix_intersecting_edges_mode = self.my_fix_intersecting_edges_mode;
            // switch off FixIntEdges if reorder not done
            if self.my_fix_intersecting_edges_mode == -1 && !reorder_ok {
                self.my_fix_intersecting_edges_mode = 0;
            }
            if self.fix_self_intersection(brep) {
                fixed = true;
            }
            self.fix_reorder(brep, false);
            self.my_fix_intersecting_edges_mode = sav_fix_intersecting_edges_mode;
        }

        if the_progress.user_break() {
            return false;
        }

        if ShapeFixRoot::need_fix(self.my_fix_lacking_mode, reorder_ok) {
            if self.fix_lacking(brep, false) {
                fixed = true;
            }
        }

        if the_progress.user_break() {
            return false;
        }

        // TEMPORARILY without special mode !!!
        let nb_edges = self.nb_edges();
        for iedge in 1..=nb_edges {
            let edge = self.my_analyzer.wire_data().unwrap().edge(iedge);
            let face = self.face();
            if self
                .my_fix_edge
                .fix_vertex_tolerance_face(brep, &edge, &face)
            {
                fixed = true;
            }
        }

        if the_progress.user_break() {
            return false;
        }

        if self.base.my_context.is_some() {
            self.update_wire(brep);
        }

        fixed
    }

    // OCCT ShapeFix_Wire.cxx L487-534 — FixReorder(bool).
    /// OCCT ShapeFix_Wire::FixReorder(theModeBoth) (cxx L487-534): performs
    /// an analysis and reorders edges in the wire using class WireOrder;
    /// the flag determines the use of miscible mode if necessary.
    pub fn fix_reorder(&mut self, brep: &mut BRep, the_mode_both: bool) -> bool {
        self.my_status_reorder = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() {
            return false;
        }

        // fix in Both mode for bi-periodic surface
        let mut sawo = ShapeAnalysisWireOrder::new();
        let bi_periodic = self
            .my_analyzer
            .surf()
            .map(|surf| {
                let s = surf.surface();
                s.is_u_periodic() && s.is_v_periodic()
            })
            .unwrap_or(false);
        if bi_periodic && the_mode_both {
            self.my_analyzer
                .check_order(brep, &mut sawo, true, true, self.my_closed_mode);
        } else {
            self.my_analyzer
                .check_order(brep, &mut sawo, true, false, self.my_closed_mode);
        }

        self.fix_reorder_ordered(brep, &sawo);

        if self.last_fix_status(ShapeExtendStatus::Fail) {
            self.my_status_reorder |= encode_status(if self.last_fix_status(ShapeExtendStatus::Fail1) {
                ShapeExtendStatus::Fail1
            } else {
                ShapeExtendStatus::Fail2
            });
        }
        if !self.last_fix_status(ShapeExtendStatus::Done) {
            return false;
        }

        self.my_status_reorder |= encode_status(ShapeExtendStatus::Done1);
        if sawo.status() == 2 || sawo.status() == -2 {
            self.my_status_reorder |= encode_status(ShapeExtendStatus::Done2);
        }
        if sawo.status() < 0 {
            self.my_status_reorder |= encode_status(ShapeExtendStatus::Done3);
        }
        if sawo.status() == 3 {
            // only shifted
            self.my_status_reorder |= encode_status(ShapeExtendStatus::Done5);
        }
        true
    }

    // OCCT ShapeFix_Wire.cxx L538-553 — FixSmall(bool, prec).
    /// OCCT ShapeFix_Wire::FixSmall(lockvtx, precsmall) (cxx L538-553):
    /// applies FixSmall(num) to all edges in the wire.
    pub fn fix_small_all(&mut self, brep: &mut BRep, lockvtx: bool, precsmall: f64) -> i32 {
        self.my_status_small = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() {
            // OCCT returns false from the int function.
            return 0;
        }

        let mut i = self.nb_edges();
        while i > 0 {
            self.fix_small_edge(brep, i, lockvtx, precsmall);
            self.my_status_small |= self.my_last_fix_status;
            i -= 1;
        }

        if self.status_small(ShapeExtendStatus::Done) {
            1
        } else {
            0
        }
    }

    // OCCT ShapeFix_Wire.cxx L557-596 — FixConnected(prec).
    /// OCCT ShapeFix_Wire::FixConnected(prec) (cxx L557-596): applies
    /// FixConnected(num) to all edges in the wire; the connection between
    /// the first and last edges is treated only if the flag ClosedMode is
    /// True; if <prec> is -1 then MaxTolerance() is taken.
    pub fn fix_connected(&mut self, brep: &mut BRep, prec: f64) -> bool {
        self.my_status_connected = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() {
            return false;
        }

        let a_stop = if self.my_closed_mode { 0 } else { 1 };
        let mut a_i = self.nb_edges();
        while a_i > a_stop {
            self.fix_connected_edge(brep, a_i, prec, false);
            self.my_status_connected |= self.my_last_fix_status;
            // Refresh the edge that the next iteration will analyze as n1.
            if self.base.my_context.is_none() || a_i - 1 <= a_stop {
                a_i -= 1;
                continue;
            }
            let nb_now = self.my_analyzer.wire_data().unwrap().nb_edges();
            let a_n1_next = if a_i - 1 > 1 { a_i - 2 } else { nb_now };
            if a_n1_next == a_i || a_n1_next == a_i - 1 {
                a_i -= 1;
                continue;
            }
            let a_e_prev = self.my_analyzer.wire_data().unwrap().edge(a_n1_next);
            let a_refresh = self.context_apply(brep, &a_e_prev);
            if a_refresh.is_null()
                || a_refresh.is_same(&a_e_prev)
                || a_refresh.shape_type() != ShapeType::Edge
            {
                a_i -= 1;
                continue;
            }
            self.my_analyzer
                .wire_data_mut()
                .unwrap()
                .set_edge(&a_refresh, a_n1_next);
            a_i -= 1;
        }

        if self.base.my_context.is_some() {
            self.update_wire(brep);
        }

        self.status_connected(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire.cxx L600-1030 — FixEdgeCurves.
    /// OCCT ShapeFix_Wire::FixEdgeCurves() (cxx L600-1030): groups the fixes
    /// dealing with 3d and pcurves of the edges: FixReversed2d /
    /// FixRemovePCurve / FixAddPCurve / FixRemoveCurve3d / FixAddCurve3d /
    /// FixSeam / FixShifted / FixSameParameter / FixVertexTolerance.
    pub fn fix_edge_curves(&mut self, brep: &mut BRep) -> bool {
        self.my_status_edge_curves = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() {
            return false;
        }
        let is_ready = self.is_ready();

        let mut nb = self.my_analyzer.wire_data().unwrap().nb_edges();
        let face = self.face();
        // OCCT: theAdvFixEdge = myFixEdge; the handle is never null in rcad
        // (the value model), so the `myFixReversed2dMode = false` fallback
        // (implicit int conversion) is dead in OCCT too.
        let _ = &mut self.my_fix_edge;

        // fix revesred 2d / 3d curves
        if is_ready && ShapeFixRoot::need_fix(self.my_fix_reversed2d_mode, true) {
            for i in 1..=nb {
                let edge = self.my_analyzer.wire_data().unwrap().edge(i);
                self.my_fix_edge.fix_reversed2d_face(brep, &edge, &face);
                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done1);
                }
                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail1);
                }
            }
        }

        // add / remove pcurve
        if is_ready && ShapeFixRoot::need_fix(self.my_fix_remove_pcurve_mode, false) {
            for i in 1..=nb {
                let edge = self.my_analyzer.wire_data().unwrap().edge(i);
                self.my_fix_edge.fix_remove_pcurve_face(brep, &edge, &face);
                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done2);
                }
                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail2);
                }
            }
        }

        if is_ready && ShapeFixRoot::need_fix(self.my_fix_add_pcurve_mode, true) {
            let mut overdegen = 0; //: c0
            let mut i = 1;
            while i <= nb {
                let edge = self.my_analyzer.wire_data().unwrap().edge(i);
                let is_seam = self.my_analyzer.wire_data_mut().unwrap().is_seam(i);
                let precision = self.base.my_precision;
                self.my_fix_edge.fix_add_pcurve_face(brep, &edge, &face, is_seam, precision);
                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done3);
                }
                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail3);
                }

                // abv 24 Feb 00: trj3_s1-ac-214.stp #1631 etc.: try to split
                // the edge in singularity
                if !is_seam && self.my_fix_edge.status(ShapeExtendStatus::Done2) {
                    if self.base.my_context.is_some() {
                        let sbe = ShapeBuildEdge;
                        let mut e = self.my_analyzer.wire_data().unwrap().edge(i);
                        let sac = ShapeAnalysisCurve;
                        let mut c: Option<rcad_kernel::geom::Curve3> = None;
                        let mut a = 0.0f64;
                        let mut b = 0.0f64;
                        {
                            let (cv, fa, fb) = brep_tool_curve_range(brep, &e);
                            c = cv;
                            a = fa;
                            b = fb;
                        }
                        let min_tol = self.base.my_min_tol;
                        let nbs = self
                            .my_analyzer
                            .my_surf
                            .as_mut()
                            .map(|s| s.nb_singularities(min_tol))
                            .unwrap_or(0);
                        let mut seq: Vec<f64> = Vec::new();
                        for j in 1..=nbs {
                            let mut preci = 0.0f64;
                            let mut pd1 = DVec2::ZERO;
                            let mut pd2 = DVec2::ZERO;
                            let mut p3d = glam::DVec3::ZERO;
                            let mut par1 = 0.0f64;
                            let mut par2 = 0.0f64;
                            let mut tmp_uiso_deg = false;
                            let ok = self
                                .my_analyzer
                                .my_surf
                                .as_mut()
                                .map(|s| {
                                    s.singularity(
                                        j,
                                        &mut preci,
                                        &mut p3d,
                                        &mut pd1,
                                        &mut pd2,
                                        &mut par1,
                                        &mut par2,
                                        &mut tmp_uiso_deg,
                                    )
                                })
                                .unwrap_or(false);
                            if !ok {
                                continue;
                            }
                            let sac2 = ShapeAnalysisCurve;
                            let mut pr = glam::DVec3::ZERO;
                            let mut split = 0.0f64;
                            let dist = c.as_ref().map(|cv| {
                                sac2.project(cv, p3d, min_tol, &mut pr, &mut split, true)
                            });
                            if let Some(dist) = dist {
                                if dist < preci.max(min_tol) {
                                    if split - a > PCONFUSION && b - split > PCONFUSION {
                                        let mut k = 1usize;
                                        while k <= seq.len() {
                                            if split < seq[k - 1] - PCONFUSION {
                                                seq.insert(k - 1, split);
                                                break;
                                            } else if split < seq[k - 1] + PCONFUSION {
                                                break;
                                            }
                                            k += 1;
                                        }
                                        if k > seq.len() {
                                            seq.push(split);
                                        }
                                    }
                                }
                            }
                        }
                        if !seq.is_empty() {
                            // supposed that edge is SP
                            let is_fwd = e.orientation == Orientation::Forward;
                            e.orientation = Orientation::Forward;

                            // 10.04.2003 skl for using trimmed lines as
                            // pcurves
                            let sae = ShapeAnalysisEdge::new();
                            if brep_tool_same_parameter(
                                &self.my_analyzer.wire_data().unwrap().edge(i),
                            ) {
                                sbe.remove_pcurve_face(brep, &e, &face);
                            } else if sae.has_pcurve_face(brep, &e, &face) {
                                let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
                                let mut fp2d = 0.0f64;
                                let mut lp2d = 0.0f64;
                                if sae.pcurve_face(brep, &e, &face, &mut c2d, &mut fp2d, &mut lp2d, false) {
                                    if let Some(c2d) = &c2d {
                                        if !matches!(
                                            c2d,
                                            rcad_kernel::geom::Curve2d::Trimmed(_)
                                        ) {
                                            sbe.remove_pcurve_face(brep, &e, &face);
                                        }
                                    }
                                }
                            }

                            let mut b_builder = BRepBuilder::new();
                            let mut v1 = sae.first_vertex(brep, &e);
                            let mut v2 = sae.last_vertex(brep, &e);
                            let mut v;

                            let mut sw = WireData::new();
                            let mut a = a;
                            let mut b = b;
                            for k in 0..=seq.len() {
                                let split = if k < seq.len() { seq[k] } else { b };
                                if k < seq.len() {
                                    let pt = c
                                        .as_ref()
                                        .map(|cv| CurveEval::point_at(cv, split))
                                        .unwrap_or(glam::DVec3::ZERO);
                                    v = b_builder.add_vertex(
                                        brep,
                                        pt,
                                        brep_tool_tolerance(&e),
                                    );
                                    // try increase tolerance before splitting
                                    let mut a_dist =
                                        brep_tool_pnt(&v1).distance(brep_tool_pnt(&v));
                                    if a_dist < brep_tool_tolerance(&v1) * 1.01 {
                                        b_builder.update_vertex_tolerance(
                                            brep,
                                            v1.clone(),
                                            a_dist.max(brep_tool_tolerance(&v1)),
                                        );
                                        a = split;
                                        v1 = v.clone();
                                        continue;
                                    } else {
                                        a_dist = brep_tool_pnt(&v2).distance(brep_tool_pnt(&v));
                                        if a_dist < brep_tool_tolerance(&v2) * 1.01 {
                                            b_builder.update_vertex_tolerance(
                                                brep,
                                                v.clone(),
                                                a_dist.max(brep_tool_tolerance(&v2)),
                                            );
                                            b = split;
                                            v2 = v.clone();
                                            continue;
                                        }
                                    }
                                } else {
                                    v = v2.clone();
                                }

                                let new_edge = sbe.copy_replace_vertices(brep, &e, &v1, &v);
                                if brep_tool_same_parameter(
                                    &self.my_analyzer.wire_data().unwrap().edge(i),
                                ) {
                                    b_builder.set_edge_range(brep, new_edge.clone(), a, split);
                                    sw.add_edge(&new_edge, 0);
                                } else {
                                    let sftp = ShapeAnalysisTransferParametersProj::new_edge_face(
                                        brep, &e, &face,
                                    );
                                    sftp.transfer_range(brep, &new_edge, a, split, false);
                                    sw.add_edge(&new_edge, 0);
                                }
                                a = split;
                                v1 = v.clone();
                            }
                            if !is_fwd {
                                sw.reverse();
                                e.orientation = Orientation::Reversed;
                            }
                            let sw_wire = sw.wire(brep);
                            if let Some(ctx) = self.base.my_context.as_mut() {
                                ctx.replace(brep, &e, &sw_wire);
                            }
                            self.update_wire(brep);
                            nb = self.my_analyzer.wire_data().unwrap().nb_edges();
                            // OCCT: i-- then the for-loop i++ — net zero, so
                            // the while-increment is skipped.
                            continue;
                        }
                    }

                    overdegen = i;
                }
                i += 1;
            }

            //: c0 abv 20 Feb 98: treat case of curve going over degenerated
            // pole and seam
            if overdegen != 0
                && self
                    .my_analyzer
                    .my_surf
                    .as_mut()
                    .map(|s| s.is_u_closed(self.base.my_precision))
                    .unwrap_or(false)
            {
                let sbe = ShapeBuildEdge;
                let suf = self.my_analyzer.surf().unwrap();
                let (mut suf_f, mut sul, mut svf, mut svl) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                suf.bounds(&mut suf_f, &mut sul, &mut svf, &mut svl);
                let u_range = (sul - suf_f).abs();
                let mut vec = DVec2::ZERO;
                let sae = ShapeAnalysisEdge::new();
                let mut k = 1;
                while k <= nb {
                    let edge = self.my_analyzer.wire_data().unwrap().edge(k);
                    let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
                    let mut cf = 0.0f64;
                    let mut cl = 0.0f64;
                    if !sae.pcurve_face(brep, &edge, &face, &mut c2d, &mut cf, &mut cl, true) {
                        break;
                    }
                    let p_first = c2d
                        .as_ref()
                        .map(|c| rcad_kernel::geom::Curve2dEval::point_at(c, cf))
                        .unwrap_or(DVec2::ZERO);
                    let p_last = c2d
                        .as_ref()
                        .map(|c| rcad_kernel::geom::Curve2dEval::point_at(c, cl))
                        .unwrap_or(DVec2::ZERO);
                    vec += p_last - p_first;
                    k += 1;
                }
                if k > nb && ((vec.x.abs() - u_range).abs()) < 0.1 * u_range {
                    let over_edge = self.my_analyzer.wire_data().unwrap().edge(overdegen);
                    sbe.remove_pcurve_face(brep, &over_edge, &face);
                    *self.my_fix_edge.projector().adjust_over_degen_mode() = 0;
                    let is_seam = self
                        .my_analyzer
                        .wire_data_mut()
                        .unwrap()
                        .is_seam(overdegen);
                    let precision = self.base.my_precision;
                    self.my_fix_edge.fix_add_pcurve_face(
                        brep,
                        &over_edge,
                        &face,
                        is_seam,
                        precision,
                    );
                }
            }
        }

        // add / remove pcurve
        if is_ready && ShapeFixRoot::need_fix(self.my_fix_remove_curve3d_mode, false) {
            for i in 1..=nb {
                let edge = self.my_analyzer.wire_data().unwrap().edge(i);
                self.my_fix_edge.fix_remove_curve3d(brep, &edge);
                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done4);
                }
                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail4);
                }
            }
        }
        if ShapeFixRoot::need_fix(self.my_fix_add_curve3d_mode, true) {
            let mut i = 1;
            while i <= nb {
                let edge = self.my_analyzer.wire_data().unwrap().edge(i);
                self.my_fix_edge.fix_add_curve3d(brep, &edge);
                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done5);
                }
                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                    //: abv 29.08.01: Spatial_firex_lofting.sat: if 3d curve
                    // cannot be built because edge has no pcurves either,
                    // remove that edge
                    let (c2d_opt, first, last) =
                        brep_tool_curve_on_surface_first(brep, &edge);
                    let _ = &c2d_opt;
                    if c2d_opt.is_none() || (last - first).abs() < PCONFUSION {
                        // Incomplete edge (with no pcurves or 3d curve)
                        // removed
                        self.base.send_warning(
                            &edge,
                            &crate::shhealing::shape_extend::msg::MessageMsg::from_key(
                                "FixWire.FixCurve3d.Removed",
                            ),
                        );
                        self.my_analyzer.wire_data_mut().unwrap().remove(i);
                        i -= 1;
                        nb -= 1;
                        self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done5);
                        if i == nb {
                            self.fix_closed(brep, self.base.my_precision);
                        } else {
                            self.fix_connected_edge(brep, i + 1, self.base.my_precision, true);
                        }
                    }
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail5);
                }
                i += 1;
            }
        }

        // fix seam
        if is_ready && ShapeFixRoot::need_fix(self.my_fix_seam_mode, false) {
            for i in 1..=nb {
                self.fix_seam(brep, i);
                if self.last_fix_status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done6);
                }
                if self.last_fix_status(ShapeExtendStatus::Fail) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail6);
                }
            }
        }

        // fix shifted
        if is_ready && ShapeFixRoot::need_fix(self.my_fix_shifted_mode, true) {
            self.fix_shifted(brep);
            if self.last_fix_status(ShapeExtendStatus::Done) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done7);
            }
            if self.last_fix_status(ShapeExtendStatus::Fail) {
                self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail7);
            }
        }

        // fix same parameter
        if is_ready && ShapeFixRoot::need_fix(self.my_fix_same_parameter_mode, true) {
            for i in 1..=nb {
                // skl 28.10.2004 for OCC6366 - check SameRange
                let sae = ShapeAnalysisEdge::new();
                let edge_i = self.my_analyzer.wire_data().unwrap().edge(i);
                let (tmpc3d, first, last) = brep_tool_curve_range(brep, &edge_i);
                let mut first = first;
                let mut last = last;
                let _ = &tmpc3d;
                if sae.has_pcurve_face(brep, &edge_i, &face) {
                    let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
                    let mut fp2d = 0.0f64;
                    let mut lp2d = 0.0f64;
                    if sae.pcurve_face(brep, &edge_i, &face, &mut c2d, &mut fp2d, &mut lp2d, false) {
                        if (first - fp2d).abs() > PCONFUSION || (last - lp2d).abs() > PCONFUSION {
                            let mut b = BRepBuilder::new();
                            b.set_edge_same_range(brep, edge_i.clone(), false);
                        } else if let Some(c2d) = &c2d {
                            if !sae.check_pcurve_range(first, last, c2d) {
                                // Replace pcurve
                                let (s, l) = brep_face_surface_loc(brep, &face);
                                if let Some(s) = s {
                                    ShapeBuildEdge.remove_pcurve_surface_loc(
                                        brep, &edge_i, &s, l,
                                    );
                                }
                                let is_seam =
                                    self.my_analyzer.wire_data_mut().unwrap().is_seam(i);
                                let precision = self.base.my_precision;
                                self.my_fix_edge.fix_add_pcurve_face(
                                    brep, &edge_i, &face, is_seam, precision,
                                );
                                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                                    self.my_status_edge_curves |=
                                        encode_status(ShapeExtendStatus::Done3);
                                }
                                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                                    self.my_status_edge_curves |=
                                        encode_status(ShapeExtendStatus::Fail3);
                                }
                            }
                        }
                    }
                }
                let face2 = self.face();
                self.my_fix_edge
                    .fix_same_parameter_face(brep, &edge_i, &face2, 0.0);
                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done8);
                }
                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail8);
                }
            }
        }

        //: abv 10.06.02: porting C40 -> dev (CC670-12608.stp): moved from
        // Perform().  Update with face is needed for plane surfaces (w/o
        // stored pcurves)
        if ShapeFixRoot::need_fix(self.my_fix_vertex_tolerance_mode, true) {
            for i in 1..=nb {
                let edge = self.my_analyzer.wire_data().unwrap().edge(i);
                let face2 = self.face();
                self.my_fix_edge.fix_vertex_tolerance_face(brep, &edge, &face2);
                if self.my_fix_edge.status(ShapeExtendStatus::Done) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Done8);
                }
                if self.my_fix_edge.status(ShapeExtendStatus::Fail) {
                    self.my_status_edge_curves |= encode_status(ShapeExtendStatus::Fail8);
                }
            }
            if self.base.my_context.is_some() {
                self.update_wire(brep);
            }
        }

        self.status_edge_curves(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire.cxx L1034-1075 — FixDegenerated().
    /// OCCT ShapeFix_Wire::FixDegenerated() (cxx L1034-1075): applies
    /// FixDegenerated(num) to all edges in the wire.
    pub fn fix_degenerated(&mut self, brep: &mut BRep) -> bool {
        self.my_status_degenerated = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let mut lastcoded = -1;
        let mut prevcoded = 0;
        let stop = if self.my_closed_mode { 0 } else { 1 };
        let mut i = self.nb_edges();
        while i > stop {
            self.fix_degenerated_edge(brep, i);
            self.my_status_degenerated |= self.my_last_fix_status;
            //: r0 abv 19 Mar 99: PRO7226.stp #489490: remove duplicated
            // degenerated edges
            let coded = if self.last_fix_status(ShapeExtendStatus::Done2) { 1 } else { 0 };
            if lastcoded == -1 {
                lastcoded = coded;
            }
            if coded != 0 && (prevcoded != 0 || (i == 1 && lastcoded != 0)) && self.nb_edges() > 1 {
                self.my_analyzer.wire_data_mut().unwrap().remove(i);
                if prevcoded == 0 {
                    i = self.nb_edges();
                }
                let edge_at_i = self.my_analyzer.wire_data().unwrap().edge(i);
                let mut b = BRepBuilder::new();
                b.set_edge_degenerated(brep, edge_at_i, false);
                prevcoded = 0;
                // OCCT: B.Degenerated(sbwd->Edge(i++), false) — the
                // post-increment cancels the loop i--, so `i` stays.
            } else {
                prevcoded = coded;
                i -= 1;
            }
        }

        self.status_degenerated(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire.cxx L1083-1280 — FixSelfIntersection.
    /// OCCT ShapeFix_Wire::FixSelfIntersection() (cxx L1083-1280): applies
    /// FixSelfIntersectingEdge(num) and FixIntersectingEdges(num) to all
    /// edges in the wire and FixIntersectingEdges(num1, num2) for all pairs
    /// num1 and num2 such that num2 >= num1 + 2, and removes wrong edges if
    /// any.
    pub fn fix_self_intersection(&mut self, brep: &mut BRep) -> bool {
        self.my_status_self_intersection = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let mut nb = self.my_analyzer.wire_data().unwrap().nb_edges();

        if ShapeFixRoot::need_fix(self.my_fix_self_intersecting_edge_mode, true) {
            if self.my_remove_loop_mode < 1 {
                for num in 1..=nb {
                    self.fix_self_intersecting_edge(brep, num);
                    self.my_status_self_intersection |= self.my_last_fix_status;
                }
            } else if self.my_remove_loop_mode == 1 {
                let mut num = 1;
                while num <= nb {
                    self.fix_self_intersecting_edge(brep, num);
                    self.my_status_self_intersection |= self.my_last_fix_status;
                    if nb < self.my_analyzer.wire_data().unwrap().nb_edges() {
                        num -= 1;
                    }
                    nb = self.my_analyzer.wire_data().unwrap().nb_edges();
                    num += 1;
                }
                self.fix_closed(brep, self.base.my_precision);
            }
        }

        if ShapeFixRoot::need_fix(self.my_fix_intersecting_edges_mode, true) {
            let mut num = if self.my_closed_mode { 1 } else { 2 };
            while nb > 1 && num <= nb {
                self.fix_intersecting_edges(brep, num);
                if self.last_fix_status(ShapeExtendStatus::Fail1) {
                    self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Fail1);
                }
                if self.last_fix_status(ShapeExtendStatus::Fail2) {
                    self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Fail2);
                }
                if !self.last_fix_status(ShapeExtendStatus::Done) {
                    num += 1;
                    continue;
                }

                if self.last_fix_status(ShapeExtendStatus::Done1) {
                    self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done1);
                }
                if self.last_fix_status(ShapeExtendStatus::Done2) {
                    self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done2);
                }
                if self.last_fix_status(ShapeExtendStatus::Done6) {
                    self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done6);
                }

                if /* ! myTopoMode ||*/ nb < 3 {
                    // #86 rln 22.03.99 sim2.igs, entity 4292: After fixing of
                    // self-intersecting BRepCheck finds one more
                    // self-intersection not found by ShapeAnalysis
                    //%15 pdn 06.04.99 repeat until fixed CTS18546-2 entity 777

                    // if the tolerance was modified we should recheck the
                    // result, if it was enough
                    if self.last_fix_status(ShapeExtendStatus::Done7) {
                        // num--;
                        self.fix_intersecting_edges(brep, num);
                    }
                    num += 1;
                    continue;
                }

                if self.last_fix_status(ShapeExtendStatus::Done4) {
                    self.my_analyzer.wire_data_mut().unwrap().remove(num);
                }
                if self.last_fix_status(ShapeExtendStatus::Done3) {
                    let nb_now = self.my_analyzer.wire_data().unwrap().nb_edges();
                    self.my_analyzer
                        .wire_data_mut()
                        .unwrap()
                        .remove(if num > 1 { num - 1 } else { nb_now + num - 1 });
                }
                if self.last_fix_status(ShapeExtendStatus::Done4)
                    || self.last_fix_status(ShapeExtendStatus::Done3)
                {
                    self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done3);
                    num = if self.my_closed_mode { 1 } else { 2 };
                    nb = self.my_analyzer.wire_data().unwrap().nb_edges();
                } else {
                    // #86 rln 22.03.99
                    //%15 pdn 06.04.99 repeat until fixed CTS18546-2 entity 777
                    self.fix_intersecting_edges(brep, num);
                    // Always revisit the fixed edge
                    // num--;
                    num += 1;
                }
            }
            if self.base.my_context.is_some() {
                self.update_wire(brep);
            }
        }

        // pdn 17.03.99 S4135 to avoid regression fixing not adjacent
        // intersection
        if ShapeFixRoot::need_fix(self.my_fix_non_adjacent_intersecting_edges_mode, true) {
            let precision = self.base.my_precision;
            let mut i_tool =
                crate::shhealing::shape_fix::intersection_tool::ShapeFixIntersectionTool::new(
                    self.base.my_context.clone(),
                    precision,
                    1.0,
                );
            let face = self.my_analyzer.face().clone();
            let mut nb_split = 0;
            let mut nb_cut = 0;
            let mut nb_removed = 0;
            {
                let sbwd = self.my_analyzer.wire_data_mut().unwrap();
                if i_tool.fix_self_intersect_wire(
                    brep,
                    sbwd,
                    &face,
                    &mut nb_split,
                    &mut nb_cut,
                    &mut nb_removed,
                ) {
                    self.my_status_self_intersection |= encode_status(ShapeExtendStatus::Done5); // gka 06.09.04
                }
            }
            if nb_split > 0 || nb_removed > 0 {
                if nb_removed > 0 {
                    self.my_status_removed_segment = true;
                }
                // Load(sbwd); commented by skl 29.12.2004 for OCC7624,
                // instead this string inserted following three strings:
                // (the OCCT handle round-trip becomes a move-out/move-in of
                // the analyzer-owned wire data).
                let sbwd = self.my_analyzer.my_wire.take().unwrap();
                self.my_analyzer.load_wire_data(sbwd);
                if self.base.my_context.is_some() {
                    self.update_wire(brep);
                }
                self.base.my_shape = Shape::null();
            }
        }

        self.status_self_intersection(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire.cxx L1284-1300 — FixLacking(bool).
    /// OCCT ShapeFix_Wire::FixLacking(force) (cxx L1284-1300): applies
    /// FixLacking(num) to all edges in the wire.
    pub fn fix_lacking(&mut self, brep: &mut BRep, force: bool) -> bool {
        self.my_status_lacking = encode_status(ShapeExtendStatus::Ok);
        if !self.is_ready() {
            return false;
        }

        let start = if self.my_closed_mode { 1 } else { 2 };
        let nb = self.nb_edges();
        for i in start..=nb {
            self.fix_lacking_edge(brep, i, force);
            self.my_status_lacking |= self.my_last_fix_status;
        }

        self.status_lacking(ShapeExtendStatus::Done)
    }

    // OCCT ShapeFix_Wire.cxx L1304-1343 — FixClosed(prec).
    /// OCCT ShapeFix_Wire::FixClosed(prec) (cxx L1304-1343): fixes a wire to
    /// be well closed; performs FixConnected, FixDegenerated and FixLacking
    /// between the last and first edges (independently of the flag
    /// ClosedMode and the modes for these fixings).
    pub fn fix_closed(&mut self, brep: &mut BRep, prec: f64) -> bool {
        self.my_status_closed = encode_status(ShapeExtendStatus::Ok);
        if !self.is_loaded() || self.nb_edges() < 1 {
            return false;
        }

        self.fix_connected_edge(brep, 1, prec, true);
        if self.last_fix_status(ShapeExtendStatus::Done) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Done1);
        }
        if self.last_fix_status(ShapeExtendStatus::Fail) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Fail1);
        }

        self.fix_degenerated_edge(brep, 1);
        if self.last_fix_status(ShapeExtendStatus::Done) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Done2);
        }
        if self.last_fix_status(ShapeExtendStatus::Fail) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Fail2);
        }

        self.fix_lacking_edge(brep, 1, false);
        if self.last_fix_status(ShapeExtendStatus::Done) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Done3);
        }
        if self.last_fix_status(ShapeExtendStatus::Fail) {
            self.my_status_closed |= encode_status(ShapeExtendStatus::Fail3);
        }

        self.status_closed(ShapeExtendStatus::Done)
    }
}

// ---------------------------------------------------------------------------
// File-local BRep_Tool re-hosts (the fix_api bodies only).
// ---------------------------------------------------------------------------

/// OCCT BRep_Tool::Curve(edge, first, last): the 3d curve and its range.
pub(crate) fn brep_tool_curve_range(
    brep: &BRep,
    edge: &Shape,
) -> (Option<rcad_kernel::geom::Curve3>, f64, f64) {
    match brep.tshapes[edge.index].as_ref() {
        TShape::Edge(ed) => (ed.curve.clone(), ed.range[0], ed.range[1]),
        _ => (None, 0.0, 0.0),
    }
}

/// OCCT BRep_Tool::CurveOnSurface(edge, C, S, L, first, last): the first
/// pcurve representation of the edge.
pub(crate) fn brep_tool_curve_on_surface_first(
    brep: &BRep,
    edge: &Shape,
) -> (Option<rcad_kernel::geom::Curve2d>, f64, f64) {
    match brep.tshapes[edge.index].as_ref() {
        TShape::Edge(ed) => {
            for r in &ed.representations {
                match r {
                    rcad_kernel::topods::CurveRepresentation::CurveOnSurface {
                        pcurve,
                        range,
                        ..
                    } => return (Some(pcurve.clone()), range[0], range[1]),
                    rcad_kernel::topods::CurveRepresentation::CurveOnClosedSurface {
                        pcurve1,
                        range,
                        ..
                    } => return (Some(pcurve1.clone()), range[0], range[1]),
                    _ => continue,
                }
            }
            (None, 0.0, 0.0)
        }
        _ => (None, 0.0, 0.0),
    }
}

/// OCCT BRep_Tool::Surface(face, L): the face surface and its location.
pub(crate) fn brep_face_surface_loc(brep: &BRep, face: &Shape) -> (Option<rcad_kernel::geom::Surface3>, u32) {
    match brep.tshapes[face.index].as_ref() {
        TShape::Face(fd) => (fd.surface.clone(), fd.surface_location),
        _ => (None, 0),
    }
}
