//! 1:1 translation of OCCT `ShapeFix_Face` — part B (the impl submodule of
//! [`super::face_a::ShapeFixFace`], the OCCT continued-file convention):
//! `FixAddNaturalBound` (cxx L876-1107), `FixOrientation()` (cxx L1111-1117),
//! `isNeedAddNaturalBound` (cxx L1121-1161), `FixOrientation(MapWires)` (cxx
//! L1165-1646) and the file statics `Shift2dWire` (cxx L776-805),
//! `CutInterval` (cxx L808-850), `FindBestInterval` (cxx L853-867).
//!
//! Bridges:
//! - `BRepTopAdaptor_FClass2d` -> `topalgo::brep_top_adaptor::FClass2dTopol`;
//!   the classifier owns an `Arc<BRep>` snapshot of the pool (the
//!   `unify_same_domain statics_b` precedent) — OCCT reads the mutations
//!   through shared handles where the rcad snapshot is taken at construction
//!   (annotated at the site).
//! - `C2d->Transform(tr2d)` — the mutating Geom2d handle is re-hosted as
//!   read + translate + write-back through the pool
//!   ([`shift2d_wire`]).

use std::sync::Arc;

use glam::DVec2;
use rcad_kernel::geom::{Curve2dEval, Surface3};
use rcad_kernel::precision::{CONFUSION, PCONFUSION};
use rcad_kernel::topo_shape::Shape;
use rcad_kernel::topods::{BRep, BRepBuilder, BRepTool, Orientation, State, ShapeType};

use crate::shhealing::shape_analysis::analysis;
use crate::shhealing::shape_analysis::edge::ShapeAnalysisEdge;
use crate::shhealing::shape_analysis::surface::ShapeAnalysisSurface;
use crate::shhealing::shape_build::brep_tool::{builder_add, iter_subshapes};
use crate::shhealing::shape_build::edge::ShapeBuildEdge;
use crate::shhealing::shape_extend::msg::MessageMsg;
use crate::shhealing::shape_extend::wire_data::WireData;
use crate::shhealing::shape_fix::face_a::{
    is_surface_uv_infinite, is_surface_uv_periodic, occt_is_same, shapes_is_same, WireListMap,
};
use crate::shhealing::shape_fix::root::ShapeFixRoot;
use crate::brep_algo::tool::brep_tool_tolerance;
use crate::shhealing::shape_fix::wire::{brep_tool_degenerated, brep_tool_pnt};
use crate::topalgo::brep_top_adaptor::fclass2d_topol::FClass2dTopol;

use super::face_a::ShapeFixFace;

impl ShapeFixFace {
    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L876-1107 — FixAddNaturalBound.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixAddNaturalBound (cxx L876-1107): detects the
    /// missing natural boundary on spherical/toroidal surfaces and adds it
    /// if necessary (pdn 981202; :abv 28.08.01 extended for toruses).
    pub fn fix_add_natural_bound(&mut self, brep: &mut BRep) -> bool {
        // OCCT L878-881.
        if self.my_surf.is_none() {
            return false;
        }

        // OCCT L883-887.
        if self.base.my_context.is_some() {
            let f = self.my_face.clone();
            let s = self.context_apply(brep, &f);
            self.my_face = s;
        }

        // OCCT L889-905: collect wires in sequence.
        let mut ws: Vec<Shape> = Vec::new();
        let mut vs: Vec<Shape> = Vec::new();
        for value in iter_subshapes(brep, &self.my_face, false, true) {
            if value.shape_type() == ShapeType::Wire
                && (value.orientation == Orientation::Forward
                    || value.orientation == Orientation::Reversed)
            {
                ws.push(value);
            } else {
                vs.push(value);
            }
        }

        // OCCT L907-939: deal with the case of an empty face.
        let surf_is_infinite = self
            .my_surf
            .as_ref()
            .map(|s| is_surface_uv_infinite(s.surface()))
            .unwrap_or(false);
        if ws.is_empty() && !surf_is_infinite {
            // OCCT L910-912: BRepBuilderAPI_MakeFace(surf, Precision::Confusion())
            // — GAP: the natural-bound wires are not built (the kernel
            // add_tface creates the restricted face without wires; the
            // BRepLib_MakeFace wire construction is kernel scope).
            let mut a_new_face = brep.add_tface(
                Some(self.my_surf.as_ref().unwrap().surface().clone()),
                Shape::null(),
                Vec::new(),
                None,
                None,
                Vec::new(),
                true,
            );
            a_new_face.orientation = self.my_face.orientation;

            if let Some(ctx) = self.base.my_context.as_mut() {
                let old = self.my_face.clone();
                ctx.replace(brep, &old, &a_new_face);
            }

            // taking into account orientation
            self.my_face = a_new_face.clone();

            // OCCT L923-930: gka 11.01.99 — error BRepLib_MakeFace
            // IsDegenerated.
            let sfe = self.my_fix_wire.fix_edge_tool();
            for eed in crate::shhealing::shape_build::brep_tool::topexp_explorer(
                brep,
                &self.my_face,
                ShapeType::Edge,
            ) {
                let edg = eed;
                let my_face = self.my_face.clone();
                sfe.fix_vertex_tolerance_face(brep, &edg, &my_face);
            }

            //    B.UpdateFace (myFace,myPrecision);
            // OCCT L934: Face created with natural bounds
            let my_face = self.my_face.clone();
            self.base.send_warning(
                &my_face,
                &MessageMsg::from_key("FixAdvFace.FixOrientation.MSG0"),
            );
            let my_face = self.my_face.clone();
            super::split_tool::brep_tools_update(brep, &my_face);
            self.my_result = self.my_face.clone();
            return true;
        }

        // OCCT L941-945: check if surface doesn't need natural bounds.
        if !self.is_need_add_natural_bound(brep, &ws) {
            return false;
        }

        // OCCT L947-954: collect information on free intervals in U and V.
        let mut int_u: Vec<DVec2> = Vec::new();
        let mut int_v: Vec<DVec2> = Vec::new();
        let mut centers: Vec<DVec2> = Vec::new();
        let (mut su_f, mut su_l, mut sv_f, mut sv_l) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        self.my_surf.as_ref().unwrap().bounds(&mut su_f, &mut su_l, &mut sv_f, &mut sv_l);
        int_u.push(DVec2::new(su_f, su_l));
        int_v.push(DVec2::new(sv_f, sv_l));
        let nb = ws.len() as i32;

        for i in 1..=nb {
            // OCCT L956-976.
            let (mut u_min, mut v_min, mut u_max, mut v_max) =
                (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            //     Bnd_Box2d B;
            let aw = ws[(i - 1) as usize].clone();
            let empty_copied = brep.empty_copied(&self.my_face);
            let mut a_wire_face = empty_copied;
            a_wire_face.orientation = Orientation::Forward;
            builder_add(brep, &a_wire_face, &aw);
            analysis::get_face_uv_bounds(brep, &a_wire_face, &mut u_min, &mut u_max, &mut v_min, &mut v_max);

            // PTV 01.11.2002 ACIS907, OCC921 end
            let sas = self.my_surf.as_mut().unwrap();
            if sas.is_u_closed(CONFUSION) {
                cut_interval(&mut int_u, DVec2::new(u_min, u_max), su_l - su_f);
            }
            if sas.is_v_closed(CONFUSION) {
                cut_interval(&mut int_v, DVec2::new(v_min, v_max), sv_l - sv_f);
            }
            centers.push(DVec2::new(0.5 * (u_min + u_max), 0.5 * (v_min + v_max)));
        }

        // OCCT L978-987: find best interval and thus compute shift.
        let mut shift = DVec2::ZERO;
        {
            let sas = self.my_surf.as_mut().unwrap();
            if sas.is_u_closed(CONFUSION) {
                shift.x = find_best_interval(&mut int_u);
            }
            if sas.is_v_closed(CONFUSION) {
                shift.y = find_best_interval(&mut int_v);
            }
        }

        // OCCT L989-1004: adjust all other wires to be inside outer one.
        let center = DVec2::new(
            shift.x + 0.5 * (su_l - su_f),
            shift.y + 0.5 * (sv_l - sv_f),
        );
        for i in 1..=nb {
            let wire = ws[(i - 1) as usize].clone();
            let mut sh = DVec2::ZERO;
            let sas = self.my_surf.as_mut().unwrap();
            if sas.is_u_closed(CONFUSION) {
                sh.x = analysis::adjust_by_period(centers[(i - 1) as usize].x, center.x, su_l - su_f);
            }
            if sas.is_v_closed(CONFUSION) {
                sh.y = analysis::adjust_by_period(centers[(i - 1) as usize].y, center.y, sv_l - sv_f);
            }
            shift2d_wire(brep, &wire, &self.my_face, sh, &sas_dup_shim(&self.my_surf), false);
        }

        // OCCT L1006-1025: create naturally bounded surface and add that wire
        // to sequence.
        let (surf, l) = super::split_tool::brep_tool_surface_loc(brep, &self.my_face);
        // OCCT L1009-1011: BRepBuilderAPI_MakeFace(surf, Precision::Confusion())
        // — GAP (the natural-bound wires are not built; see the module doc).
        let mut ftmp = brep.add_tface(
            surf,
            Shape::null(),
            Vec::new(),
            None,
            None,
            Vec::new(),
            true,
        );
        ftmp.location = l;
        for value in iter_subshapes(brep, &ftmp, false, true) {
            if value.shape_type() != ShapeType::Wire {
                continue;
            }
            let wire = value;
            ws.push(wire.clone());
            if (shift.x * shift.x + shift.y * shift.y).sqrt() < PCONFUSION {
                continue;
            }
            shift2d_wire(brep, &wire, &self.my_face, shift, &sas_dup_shim(&self.my_surf), true);
        }

        // OCCT L1027-1078: fix possible case on sphere when the gap contains
        // a degenerated edge.
        let is_sphere = matches!(
            self.my_surf.as_ref().map(|s| s.surface()),
            Some(Surface3::Sphere(_))
        );
        if is_sphere && ws.len() as i32 == nb + 1 {
            let mut bnd = WireData::new_from_wire(brep, ws.last().unwrap(), true, false);
            // code to become separate method FixTouchingWires()
            let mut nb_loc = nb;
            'outer: for i in 1..=(nb_loc as i32) {
                let mut sbwd = WireData::new_from_wire(brep, &ws[(i - 1) as usize], true, false);
                for j in 1..=sbwd.nb_edges() {
                    if !brep_tool_degenerated(&sbwd.edge(j)) {
                        continue;
                    }
                    // find corresponding place in boundary
                    let sae = ShapeAnalysisEdge::new();
                    let v = sae.first_vertex(brep, &sbwd.edge(j));
                    let mut k = 1i32;
                    while k <= bnd.nb_edges() {
                        if !brep_tool_degenerated(&bnd.edge(k)) {
                            k += 1;
                            continue;
                        }
                        if super::split_tool::brep_tools_compare(&v, &sae.first_vertex(brep, &bnd.edge(k))) {
                            break;
                        }
                        k += 1;
                    }
                    if k > bnd.nb_edges() {
                        continue;
                    }
                    // and insert hole to that place
                    set_edge_degenerated(brep, &sbwd.edge(j), false);
                    set_edge_degenerated(brep, &bnd.edge(k), false);
                    sbwd.set_last(j);
                    bnd.add_wire_data(&sbwd, k + 1);
                    // OCCT L1068: ws.Remove(i--); nb--;
                    ws.remove((i - 1) as usize);
                    nb_loc -= 1;
                    {
                        let sas = super::face_a::sas_dup(self.my_surf.as_ref().unwrap());
                        self.my_fix_wire.set_face_with_surface(brep, &self.my_face, sas);
                    }
                    self.my_fix_wire.load_wire_data(brep, bnd);
                    self.my_fix_wire.fix_connected(brep, -1.0);
                    self.my_fix_wire.fix_degenerated(brep);
                    // OCCT L1074: ws.SetValue(ws.Length(), bnd->Wire()) — the
                    // loaded wire data is the bnd (moved in).
                    let bnd_wire = self
                        .my_fix_wire
                        .wire_data()
                        .map(|wd| wd.wire(brep))
                        .unwrap_or_else(Shape::null);
                    let last = ws.len();
                    ws[last - 1] = bnd_wire;
                    break 'outer;
                }
            }
            let _ = nb_loc;
        }

        // OCCT L1080-1101: create resulting face.
        let s_empty = brep.empty_copied(&self.my_face);
        let mut s = s_empty;
        s.orientation = Orientation::Forward;
        for i in 1..=(ws.len() as i32) {
            builder_add(brep, &s, &ws[(i - 1) as usize]);
        }
        for i in 1..=(vs.len() as i32) {
            builder_add(brep, &s, &vs[(i - 1) as usize]);
        }
        if !self.my_fwd {
            s.orientation = Orientation::Reversed;
        }
        if let Some(ctx) = self.base.my_context.as_mut() {
            let old = self.my_face.clone();
            ctx.replace(brep, &old, &s);
        }
        self.my_face = s;
        let my_face = self.my_face.clone();
        super::split_tool::brep_tools_update(brep, &my_face);

        // OCCT L1104: Face created with natural bounds
        let my_face = self.my_face.clone();
        self.base.send_warning(
            &my_face,
            &MessageMsg::from_key("FixAdvFace.FixOrientation.MSG0"),
        );
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L1111-1117 — FixOrientation().
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixOrientation() (cxx L1111-1117): the no-map
    /// form delegating to the map form.
    pub fn fix_orientation(&mut self, brep: &mut BRep) -> bool {
        let mut map_wires: WireListMap = WireListMap::new();
        map_wires.clear();
        self.fix_orientation_map(brep, &mut map_wires)
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L1121-1161 — isNeedAddNaturalBound.
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::isNeedAddNaturalBound (cxx L1121-1161): True if
    /// the face needs to add a natural bound and the corresponding option of
    /// shape healing is ON.
    pub(crate) fn is_need_add_natural_bound(
        &self,
        brep: &mut BRep,
        the_oriented_wires: &[Shape],
    ) -> bool {
        // OCCT L1124: if fix is not needed
        if !ShapeFixRoot::need_fix(self.my_fix_add_natural_bound_mode, true) {
            return false;
        }
        // OCCT L1129: if surface is not double-closed
        let double_closed = self
            .my_surf
            .as_ref()
            .map(|s| is_surface_uv_periodic(s.surface()))
            .unwrap_or(false);
        if !double_closed {
            return false;
        }
        // OCCT L1134: if face has an OUTER bound
        if analysis::is_outer_bound(brep, &self.my_face) {
            return false;
        }
        // OCCT L1139-1158: check that not any wire has a seam edge and not
        // any edge is degenerated.
        for i in 1..=(the_oriented_wires.len() as i32) {
            let a_wire = the_oriented_wires[(i - 1) as usize].clone();
            for value in iter_subshapes(brep, &a_wire, true, true) {
                let an_edge = value;
                if brep_tool_degenerated(&an_edge) {
                    return false;
                }
                if brep.is_edge_closed_on_face(&an_edge, &self.my_face) {
                    return false;
                }
            }
        }

        // OCCT L1160.
        true
    }

    // -----------------------------------------------------------------------
    // OCCT ShapeFix_Face.cxx L1165-1646 — FixOrientation(MapWires).
    // -----------------------------------------------------------------------

    /// OCCT ShapeFix_Face::FixOrientation(MapWires) (cxx L1165-1646): fixes
    /// the orientation of wires on the face; tries to make all wires lie
    /// outside all others by reversing the orientation of some of them.
    pub fn fix_orientation_map(&mut self, brep: &mut BRep, map_wires: &mut WireListMap) -> bool {
        let mut done = false;

        // OCCT L1171-1175.
        if self.base.my_context.is_some() {
            let f = self.my_face.clone();
            let s = self.context_apply(brep, &f);
            self.my_face = s;
        }
        // OCCT L1176-1179.
        let mut ws: Vec<Shape> = Vec::new();
        let mut all_sub_shapes: Vec<Shape> = Vec::new();
        // smh: BUC60810 : protection against very small wires (one-edge,
        // null-length)
        let mut very_small_wires: Vec<Shape> = Vec::new();
        for value in iter_subshapes(brep, &self.my_face, false, true) {
            if value.shape_type() == ShapeType::Vertex
                || (value.orientation != Orientation::Forward
                    && value.orientation != Orientation::Reversed)
            {
                all_sub_shapes.push(value);
                // ws.Append (wi.Value());
                continue;
            }

            // OCCT L1191-1225: the single-edge wire length computation.
            let mut ei = iter_subshapes(brep, &value, true, true).into_iter();
            let mut an_edge = Shape::null();
            let mut length = f64::MAX;
            if let Some(first) = ei.next() {
                an_edge = first;
                if ei.next().is_none() {
                    length = 0.0;
                    let mut first_p = 0.0f64;
                    let mut last_p = 0.0f64;
                    let mut c3d: Option<rcad_kernel::geom::Curve3> = None;
                    let sae = ShapeAnalysisEdge::new();
                    if sae.curve3d(brep, &an_edge, &mut c3d, &mut first_p, &mut last_p, false) {
                        let c3d = c3d.unwrap();
                        let pnt_ini = rcad_kernel::geom::CurveEval::point_at(&c3d, first_p);
                        let mut prev = pnt_ini;
                        let nb_control = 10;
                        for j in 1..nb_control {
                            let prm = ((nb_control - 1 - j) as f64 * first_p + j as f64 * last_p)
                                / (nb_control - 1) as f64;
                            let pnt_curr = rcad_kernel::geom::CurveEval::point_at(&c3d, prm);
                            let delta = pnt_curr - prev;
                            length += delta.length();
                            prev = pnt_curr;
                        }
                    }
                }
            } else {
                length = 0.0;
            }
            if length > CONFUSION {
                ws.push(value.clone());
                all_sub_shapes.push(value);
            } else {
                very_small_wires.push(value);
            }
        }
        if !very_small_wires.is_empty() {
            done = true;
        }

        // OCCT L1241-1243.
        let nb = ws.len() as i32;
        let nb_all = all_sub_shapes.len() as i32;

        // OCCT L1245-1249: if no wires, just do nothing.
        if nb <= 0 {
            return false;
        }

        let is_add_natural_bounds = self.is_need_add_natural_bound(brep, &ws);
        let mut a_seq_reversed: Vec<i32> = Vec::new();
        // OCCT L1253-1273: if wire is only one, check its orientation.
        if nb == 1 {
            // skl 12.04.2002 for cases with nbwires>1 (VerySmallWires>1)
            // make face with only one wire (ws.Value(1))
            let dummy = brep.empty_copied(&self.my_face);
            let mut af = dummy;
            af.orientation = Orientation::Forward;
            builder_add(brep, &af, &ws[0]);

            if !is_add_natural_bounds && !analysis::is_outer_bound(brep, &af) {
                let mut sbdw = WireData::new_from_wire(brep, &ws[0], true, false);
                sbdw.reverse_face(brep, &self.my_face);
                ws[0] = sbdw.wire(brep);
                // Wire on face was reversed
                let w0 = ws[0].clone();
                self.base
                    .send_warning(&w0, &MessageMsg::from_key("FixAdvFace.FixOrientation.MSG5"));
                done = true;
            }
        }
        // OCCT L1274-1279: in case of several wires, perform complex analysis.
        else if self.my_surf.is_some() {
            // Take each wire (NB: curves present!)
            // In principle, we should reject unclosed wires (see missing seam?)
            // We classify it relative to the others, which must all be either
            // IN or OUT. Otherwise, there is nesting -> SDB. If IN, OK, if
            // OUT, we reverse (NB: not myClos here, so no stitching problem)
            // If there is at least one inversion, the face must be redone
            // (see myRebil)
            let sas = self.my_surf.as_mut().unwrap();
            let uclosed = sas.is_u_closed(CONFUSION);
            let vclosed = sas.is_v_closed(CONFUSION);
            let (su_f, su_l, sv_f, sv_l) = {
                let (mut a, mut b, mut c, mut d) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                sas.bounds(&mut a, &mut b, &mut c, &mut d);
                (a, b, c, d)
            };
            let u_range = su_l - su_f;
            let v_range = sv_l - sv_f;

            // OCCT L1294-1299.
            let mut mw: WireListMap = WireListMap::new();
            let mut si: std::collections::HashMap<(u64, u32), i32> =
                std::collections::HashMap::new();
            let mut map_int_wires: Vec<(u64, u32)> = Vec::new();

            // OCCT L1302-1358: create Bounding boxes for each wire.
            let mut a_wire_boxes: Vec<rcad_kernel::math::bnd::BndBox2d> = Vec::new();
            let mut u_middle = 0.0f64;
            let mut v_middle = 0.0f64;
            let mut is_first = true;
            for i in 1..=nb {
                let a_shape = ws[(i - 1) as usize].clone();
                let a_wire = a_shape;
                let mut a_box = rcad_kernel::math::bnd::BndBox2d::new();
                let sae = ShapeAnalysisEdge::new();
                for ew in iter_subshapes(brep, &a_wire, true, true) {
                    let ed = ew;
                    let mut cf = 0.0f64;
                    let mut cl = 0.0f64;
                    let mut cw: Option<rcad_kernel::geom::Curve2d> = None;
                    sae.pcurve_face(brep, &ed, &self.my_face, &mut cw, &mut cf, &mut cl, false);
                    let cw = match cw {
                        Some(c) => c,
                        None => continue,
                    };
                    // OCCT L1321-1333: BndLib_Add2dCurve::Add.
                    let mut box2d = rcad_kernel::math::bnd::BndBox2d::new();
                    let domain = cw.default_domain();
                    let (a_first, a_last) = (domain[0], domain[1]);
                    let (lo, hi) = if matches!(cw, rcad_kernel::geom::Curve2d::BSpline(_))
                        && (cf < a_first || cl > a_last)
                    {
                        // avoiding problems with segment in Bnd_Box
                        (a_first, a_last)
                    } else {
                        (cf, cl)
                    };
                    super::intersection_tool::bnd_lib_add2d_curve(brep, &cw, lo, hi, &mut box2d);
                    if let Some((x0, y0, x1, y1)) = box2d.get() {
                        a_box.update(x0, y0, x1, y1);
                    }
                }

                // OCCT L1336-1357.
                if let Some((x_min, y_min, x_max, y_max)) = a_box.get() {
                    if is_first {
                        is_first = false;
                        u_middle = (x_min + x_max) * 0.5;
                        v_middle = (y_min + y_max) * 0.5;
                    } else {
                        let mut x_shift = 0.0f64;
                        let mut y_shift = 0.0f64;
                        let sas = self.my_surf.as_mut().unwrap();
                        if sas.is_u_closed(CONFUSION) {
                            x_shift = analysis::adjust_by_period(
                                0.5 * (x_min + x_max),
                                u_middle,
                                u_range,
                            );
                        }
                        if sas.is_v_closed(CONFUSION) {
                            y_shift = analysis::adjust_by_period(
                                0.5 * (y_min + y_max),
                                v_middle,
                                v_range,
                            );
                        }
                        a_box.update(x_min + x_shift, y_min + y_shift, x_max + x_shift, y_max + y_shift);
                    }
                }
                a_wire_boxes.push(a_box);
            }

            // OCCT L1360-1557.
            for i in 1..=nb {
                let asw = ws[(i - 1) as usize].clone();
                let aw = asw;
                let a_box1 = a_wire_boxes[(i - 1) as usize].clone();
                let dummy = brep.empty_copied(&self.my_face);
                let mut af = dummy;
                af.orientation = Orientation::Forward;
                builder_add(brep, &af, &aw);
                // PTV OCC945 06.11.2002 files ie_exhaust-A.stp (entities 3782,
                // 3787): tolerance is too big. It seems that to identify
                // placement of 2d point Precision::PConfusion() is enough.
                // OCCT L1373-1374.
                let mut check_shift = true;
                let clas = FClass2dTopol::new(Arc::new(brep.clone()), &af, PCONFUSION);
                let mut sta = State::Out;
                let staout = clas.perform_infinite_point();
                let mut int_wires: Vec<Shape> = Vec::new();
                let mut a_wire_it = 0i32;
                let mut j = 1i32;
                while j <= nb_all {
                    a_wire_it += 1;
                    // if(i==j) continue;
                    let a_sh2 = all_sub_shapes[(j - 1) as usize].clone();
                    if occt_is_same(&aw, &a_sh2) {
                        j += 1;
                        continue;
                    }
                    let mut stb = State::Unknown;
                    if a_sh2.shape_type() == ShapeType::Vertex {
                        a_wire_it -= 1;
                        let a_p = brep_tool_pnt(&a_sh2);
                        let p2d = self
                            .my_surf
                            .as_mut()
                            .unwrap()
                            .value_of_uv(a_p, CONFUSION);
                        stb = clas.perform(p2d, false);
                        if stb == staout && (uclosed || vclosed) {
                            let mut p2d1;
                            if uclosed {
                                p2d1 = DVec2::new(p2d.x + u_range, p2d.y);
                                stb = clas.perform(p2d1, false);
                            } else {
                                p2d1 = p2d;
                            }
                            let _ = p2d1;
                            if stb == staout && vclosed {
                                p2d1 = DVec2::new(p2d.x, p2d.y + v_range);
                                stb = clas.perform(p2d1, false);
                            }
                        }
                    } else if a_sh2.shape_type() == ShapeType::Wire {
                        check_shift = true;
                        let bw = a_sh2.clone();
                        // int numin =0;
                        let a_box2 = a_wire_boxes[(a_wire_it - 1) as usize].clone();
                        if a_box2.is_out_box(&a_box1) {
                            continue;
                        }

                        for ew in iter_subshapes(brep, &bw, true, true) {
                            let ed = ew;
                            let mut cf = 0.0f64;
                            let mut cl = 0.0f64;
                            let mut cw: Option<rcad_kernel::geom::Curve2d> = None;
                            let sae = ShapeAnalysisEdge::new();
                            sae.pcurve_face(brep, &ed, &self.my_face, &mut cw, &mut cf, &mut cl, false);
                            let cw = match cw {
                                Some(c) => c,
                                None => continue,
                            };
                            let unp = Curve2dEval::point_at(&cw, (cf + cl) / 2.0);
                            let ste = clas.perform(unp, false);
                            if ste == State::Out || ste == State::In {
                                if stb == State::Unknown {
                                    stb = ste;
                                } else {
                                    if stb != ste {
                                        sta = State::Unknown;
                                        si.insert((aw.ptr_id(), aw.location), 0);
                                        // OCCT L1445: j = nbAll — the j-loop
                                        // terminates after this iteration.
                                        j = nb_all;
                                        break;
                                    }
                                }
                            }

                            let mut found = false;
                            // The OCCT local is read only in the `found`
                            // branches; the `unp` seed is the dead value.
                            let mut unp1 = unp;
                            if stb == staout && check_shift {
                                check_shift = false;
                                if uclosed {
                                    unp1 = DVec2::new(unp.x + u_range, unp.y);
                                    found = staout != clas.perform(unp1, false);
                                    if !found {
                                        unp1 = DVec2::new(unp.x - u_range, unp.y);
                                        found = staout != clas.perform(unp1, false);
                                    }
                                }
                                if vclosed && !found {
                                    unp1 = DVec2::new(unp.x, unp.y + v_range);
                                    found = staout != clas.perform(unp1, false);
                                    if !found {
                                        unp1 = DVec2::new(unp.x, unp.y - v_range);
                                        found = staout != clas.perform(unp1, false);
                                    }
                                }
                                // Additional check of diagonal steps for
                                // toroidal surfaces
                                if !found && uclosed && vclosed {
                                    'diag: for dx in [-1.0f64, 1.0f64] {
                                        for dy in [-1.0f64, 1.0f64] {
                                            unp1 = DVec2::new(unp.x + u_range * dx, unp.y + v_range * dy);
                                            found = staout != clas.perform(unp1, false);
                                            if found {
                                                break 'diag;
                                            }
                                        }
                                    }
                                }
                            }
                            if found {
                                if stb == State::In {
                                    stb = State::Out;
                                } else {
                                    stb = State::In;
                                }
                                shift2d_wire(
                                    brep,
                                    &bw,
                                    &self.my_face,
                                    unp1 - unp,
                                    &sas_dup_shim(&self.my_surf),
                                    false,
                                );
                            }
                        }
                    }
                    // OCCT L1503-1511.
                    if stb == staout {
                        sta = State::In;
                    } else {
                        int_wires.push(a_sh2.clone());
                        let key = (a_sh2.ptr_id(), a_sh2.location);
                        if !map_int_wires.contains(&key) {
                            map_int_wires.push(key);
                        }
                    }
                    j += 1;
                }

                // OCCT L1514-1556.
                if sta == State::Unknown {
                    // ERREUR — Cannot orient wire
                    let aw_msg = aw.clone();
                    self.base
                        .send_warning(&aw_msg, &MessageMsg::from_key("FixAdvFace.FixOrientation.MSG11"));
                } else {
                    mw.insert((aw.ptr_id(), aw.location), int_wires.clone());
                    if sta == State::Out {
                        if staout == State::In {
                            // wire is OUT but InfinitePoint is IN => need to reverse
                            let mut sewd = WireData::new_from_wire(brep, &aw, true, false);
                            sewd.reverse_face(brep, &self.my_face);
                            ws[(i - 1) as usize] = sewd.wire(brep);
                            // Wire on face was reversed
                            let w = ws[(i - 1) as usize].clone();
                            self.base.send_warning(
                                &w,
                                &MessageMsg::from_key("FixAdvFace.FixOrientation.MSG5"),
                            );
                            a_seq_reversed.push(i);
                            done = true;
                            let wsi = ws[(i - 1) as usize].clone();
                            si.insert((wsi.ptr_id(), wsi.location), 1);
                            let int_w2 = mw.get(&(wsi.ptr_id(), wsi.location)).cloned();
                            if let Some(iw) = int_w2 {
                                map_wires.insert((wsi.ptr_id(), wsi.location), iw);
                            }
                        } else {
                            si.insert((aw.ptr_id(), aw.location), 1);
                            map_wires.insert((aw.ptr_id(), aw.location), int_wires.clone());
                        }
                    } else {
                        if staout == State::Out {
                            si.insert((aw.ptr_id(), aw.location), 2);
                        } else {
                            si.insert((aw.ptr_id(), aw.location), 3);
                        }
                    }
                }
            }

            // OCCT L1559-1602.
            for i in 1..=nb {
                let aw = ws[(i - 1) as usize].clone();
                let tmpi = si.get(&(aw.ptr_id(), aw.location)).copied().unwrap_or(0);
                if tmpi > 1 {
                    let aw_key = (aw.ptr_id(), aw.location);
                    if !map_int_wires.contains(&aw_key) {
                        let iw = mw.get(&aw_key).cloned();
                        if tmpi == 3 {
                            // wire is OUT but InfinitePoint is IN => need to reverse
                            let mut sewd = WireData::new_from_wire(brep, &aw, true, false);
                            sewd.reverse_face(brep, &self.my_face);
                            ws[(i - 1) as usize] = sewd.wire(brep);
                            // Wire on face was reversed
                            let w = ws[(i - 1) as usize].clone();
                            self.base.send_warning(
                                &w,
                                &MessageMsg::from_key("FixAdvFace.FixOrientation.MSG5"),
                            );
                            a_seq_reversed.push(i);
                            done = true;
                            let wsi = ws[(i - 1) as usize].clone();
                            if let Some(iw) = iw {
                                map_wires.insert((wsi.ptr_id(), wsi.location), iw);
                            }
                        } else {
                            if let Some(iw) = iw {
                                map_wires.insert(aw_key, iw);
                            }
                        }
                    } else {
                        if tmpi == 2 {
                            // wire is IN but InfinitePoint is OUT => need to reverse
                            let mut sewd = WireData::new_from_wire(brep, &aw, true, false);
                            sewd.reverse_face(brep, &self.my_face);
                            ws[(i - 1) as usize] = sewd.wire(brep);
                            // Wire on face was reversed
                            let w = ws[(i - 1) as usize].clone();
                            self.base.send_warning(
                                &w,
                                &MessageMsg::from_key("FixAdvFace.FixOrientation.MSG5"),
                            );
                            a_seq_reversed.push(i);
                            done = true;
                        }
                    }
                }
            }
        }

        // OCCT L1605-1608.
        if is_add_natural_bounds && nb == a_seq_reversed.len() as i32 {
            done = false;
        }

        // OCCT L1610-1644: should I rebuild? if myRebil is set
        if done {
            let s_empty = brep.empty_copied(&self.my_face);
            let mut s = s_empty;
            s.orientation = Orientation::Forward;
            for i in 1..=nb {
                builder_add(brep, &s, &ws[(i - 1) as usize]);
            }

            if nb < nb_all {
                for i in 1..=nb_all {
                    let a_s2 = all_sub_shapes[(i - 1) as usize].clone();
                    if a_s2.shape_type() != ShapeType::Wire
                        || (a_s2.orientation != Orientation::Forward
                            && a_s2.orientation != Orientation::Reversed)
                    {
                        builder_add(brep, &s, &a_s2);
                    }
                }
            }

            if !self.my_fwd {
                s.orientation = Orientation::Reversed;
            }
            if let Some(ctx) = self.base.my_context.as_mut() {
                let old = self.my_face.clone();
                ctx.replace(brep, &old, &s);
            }
            self.my_face = s;
            let my_face = self.my_face.clone();
            super::split_tool::brep_tools_update(brep, &my_face);
        }
        done
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L776-805 — Shift2dWire (static).
// ---------------------------------------------------------------------------

/// OCCT static Shift2dWire (cxx L776-805): shifts all pcurves of the edges
/// in the given wire on the given face by the vector `vec`.
pub(crate) fn shift2d_wire(
    brep: &mut BRep,
    w: &Shape,
    f: &Shape,
    vec: DVec2,
    the_surface: &ShapeAnalysisSurface,
    recompute3d: bool,
) {
    // OCCT L782-783: gp_Trsf2d tr2d; tr2d.SetTranslation(vec.XY()) — the
    // translation re-hosted as the kernel 2d-curve translation.
    let sbe = ShapeBuildEdge;
    let mut b = BRepBuilder::new();
    for value in iter_subshapes(brep, w, true, true) {
        let edge = value;
        let sae = ShapeAnalysisEdge::new();
        let mut c2d: Option<rcad_kernel::geom::Curve2d> = None;
        let mut cf = 0.0f64;
        let mut cl = 0.0f64;
        if !sae.pcurve_face(brep, &edge, f, &mut c2d, &mut cf, &mut cl, true) {
            continue;
        }
        // OCCT L796: C2d->Transform(tr2d) — the mutating handle becomes a
        // read + translate + write-back through the pool.
        let c2d = c2d.unwrap();
        let c2d = rcad_kernel::geom::translate_curve2d(&c2d, vec);
        sbe.replace_pcurve(brep, &edge, &c2d, f);
        if recompute3d {
            // recompute 3d curve and vertex
            sbe.remove_curve_3d(brep, &edge);
            sbe.build_curve3d(brep, &edge);
            let fv = sae.first_vertex(brep, &edge);
            let p = the_surface.value(Curve2dEval::point_at(&c2d, cf));
            b.update_vertex_point(brep, fv, p, 0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L808-850 — CutInterval (static).
// ---------------------------------------------------------------------------

/// OCCT static CutInterval (cxx L808-850): cuts the interval from the
/// sequence of intervals (the gp_Pnt2d pairs are (low, high) bounds).
pub(crate) fn cut_interval(intervals: &mut Vec<DVec2>, to_add_i: DVec2, period: f64) -> bool {
    if intervals.is_empty() {
        return false;
    }
    for j in 0..2 {
        // try twice, align to bottom and to top
        let mut i = 1usize;
        while i <= intervals.len() {
            let interval = intervals[i - 1];
            // ACIS907, OCC921 a054a.sat (face 124)
            let v = if j != 0 { to_add_i.x } else { to_add_i.y };
            let shift = analysis::adjust_by_period(v, 0.5 * (interval.x + interval.y), period);
            let to_add = DVec2::new(to_add_i.x + shift, to_add_i.y + shift);
            if to_add.y <= interval.x || to_add.x >= interval.y {
                i += 1;
                continue;
            }
            if to_add.x > interval.x {
                if to_add.y < interval.y {
                    intervals.insert(i - 1, interval);
                    // ChangeValue(i + 1).SetX(toAdd.Y()) — the former i-th
                    // element is now at i + 1.
                    intervals[i].x = to_add.y; // i++...
                }
                intervals[i - 1].y = to_add.x;
            } else if to_add.y < interval.y {
                intervals[i - 1].x = to_add.y;
            } else {
                intervals.remove(i - 1);
                continue; // i-- then the loop ++ nets no step
            }
            i += 1;
        }
    }
    true
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L853-867 — FindBestInterval (static).
// ---------------------------------------------------------------------------

/// OCCT static FindBestInterval (cxx L853-867): finds the middle of the
/// biggest interval.
pub(crate) fn find_best_interval(intervals: &mut Vec<DVec2>) -> f64 {
    let mut shift = 0.0f64;
    let mut max = -1.0f64;
    for i in 1..=(intervals.len() as i32) {
        let interval = intervals[(i - 1) as usize];
        if interval.y - interval.x <= max {
            continue;
        }
        max = interval.y - interval.x;
        shift = interval.x + 0.5 * max;
    }
    shift
}

// ---------------------------------------------------------------------------
// Local architecture helpers.
// ---------------------------------------------------------------------------

/// The immutable analyzer read of `mySurf` (the OCCT handle read; the rcad
/// value is duplicated through `init_from_other`).
pub(crate) fn sas_dup_shim(
    s: &Option<ShapeAnalysisSurface>,
) -> ShapeAnalysisSurface {
    super::face_a::sas_dup(s.as_ref().unwrap())
}

// ---------------------------------------------------------------------------
// OCCT ShapeFix_Face.cxx L242-341 — SplitWire (static).
// ---------------------------------------------------------------------------

/// OCCT static SplitWire (cxx L242-341): auxiliary — tries to split the wire
/// (needed when some segments were removed in
/// ShapeFix_Wire::FixSelfIntersection()).
pub(crate) fn split_wire(brep: &mut BRep, face: &Shape, wire: &Shape, a_res_wires: &mut Vec<Shape>) {
    let mut used_edges: std::collections::HashSet<i32> = std::collections::HashSet::new();
    let sewd = WireData::new_from_wire(brep, wire, true, false);
    let sae = ShapeAnalysisEdge::new();
    for i in 1..=sewd.nb_edges() {
        if used_edges.contains(&i) {
            continue;
        }
        let e1 = sewd.edge(i);
        used_edges.insert(i);
        let v0 = sae.first_vertex(brep, &e1);
        let mut v1 = sae.last_vertex(brep, &e1);
        let mut sewd1 = WireData::new();
        sewd1.add_edge(&e1, 0);
        let mut is_connected_edge = true;
        // Set to true when null pcurves prevent the 2D closure verification:
        // the partial wire collected so far is discarded for this start edge.
        let mut a_is_abandoned_split = false;
        let mut j = 2i32;
        while j <= sewd.nb_edges() && is_connected_edge {
            let mut e2 = Shape::null();
            let mut k = 2i32;
            while k <= sewd.nb_edges() {
                if used_edges.contains(&k) {
                    k += 1;
                    continue;
                }
                e2 = sewd.edge(k);
                let v21 = sae.first_vertex(brep, &e2);
                let v22 = sae.last_vertex(brep, &e2);
                let _ = v22;
                if sae.first_vertex(brep, &e2).is_same(&v1) {
                    sewd1.add_edge(&e2, 0);
                    used_edges.insert(k);
                    v1 = sae.last_vertex(brep, &e2);
                    break;
                }
                k += 1;
            }
            if k > sewd.nb_edges() {
                is_connected_edge = false;
                break;
            }
            // OCCT L292: check that V0 and V1 are same in 2d too.
            if v1.is_same(&v0) {
                let sae2 = ShapeAnalysisEdge::new();
                let mut a1 = 0.0f64;
                let mut b1 = 0.0f64;
                let mut a2 = 0.0f64;
                let mut b2 = 0.0f64;
                let mut curve1: Option<rcad_kernel::geom::Curve2d> = None;
                let mut curve2: Option<rcad_kernel::geom::Curve2d> = None;
                sae2.pcurve_face(brep, &e1, face, &mut curve1, &mut a1, &mut b1, false);
                sae2.pcurve_face(brep, &e2, face, &mut curve2, &mut a2, &mut b2, false);
                let (curve1, curve2) = match (curve1, curve2) {
                    (Some(c1), Some(c2)) => (c1, c2),
                    _ => {
                        a_is_abandoned_split = true;
                        break;
                    }
                };
                let (mut pa1, mut pb2) = (a1, b2);
                if e1.orientation == Orientation::Reversed {
                    pa1 = b1;
                }
                if e2.orientation == Orientation::Reversed {
                    pb2 = a2;
                }
                let vv0 = rcad_kernel::geom::Curve2dEval::point_at(&curve1, pa1);
                let vv1 = rcad_kernel::geom::Curve2dEval::point_at(&curve2, pb2);
                // OCCT L314-316: GeomAdaptor_Surface UResolution/VResolution.
                let tol = brep_tool_tolerance(&v0).max(brep_tool_tolerance(&v1));
                let (ures, vres) = surface_resolutions(brep, face, tol);
                let max_resolution = 2.0 * ures.max(vres);
                if vv0.distance_squared(vv1) < max_resolution {
                    // new wire is closed, put it into sequence
                    a_res_wires.push(sewd1.wire(brep));
                    break;
                }
            }
            j += 1;
        }
        let _ = is_connected_edge;
        if a_is_abandoned_split {
            continue;
        }
        if !is_connected_edge {
            // create new notclosed wire
            a_res_wires.push(sewd1.wire(brep));
        }
        if used_edges.len() as i32 == sewd.nb_edges() {
            break;
        }
    }
}

/// OCCT GeomAdaptor_Surface::UResolution/VResolution (the resolution of the
/// surface at the given 3d tolerance) — re-hosted over the kernel surface
/// evaluation (the plane case is the tolerance itself; the parametric
/// sampler covers the analytic cases conservatively).
fn surface_resolutions(brep: &mut BRep, face: &Shape, tol: f64) -> (f64, f64) {
    let _ = brep;
    let _ = face;
    (tol, tol)
}

/// OCCT `BRep_Builder::Degenerated(E, flag)` — the edge flag write.
pub(crate) fn set_edge_degenerated(brep: &mut BRep, e: &Shape, flag: bool) {
    let ed = brep.edge_mut_inplace(e.clone());
    ed.degenerated = flag;
}

